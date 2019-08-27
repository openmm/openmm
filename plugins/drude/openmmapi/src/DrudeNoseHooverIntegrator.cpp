/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2019 Stanford University and the Authors.           *
 * Authors: Andreas Krämer and Andrew C. Simmonett                            *
 * Contributors:                                                              *
 *                                                                            *
 * Permission is hereby granted, free of charge, to any person obtaining a    *
 * copy of this software and associated documentation files (the "Software"), *
 * to deal in the Software without restriction, including without limitation  *
 * the rights to use, copy, modify, merge, publish, distribute, sublicense,   *
 * and/or sell copies of the Software, and to permit persons to whom the      *
 * Software is furnished to do so, subject to the following conditions:       *
 *                                                                            *
 * The above copyright notice and this permission notice shall be included in *
 * all copies or substantial portions of the Software.                        *
 *                                                                            *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR *
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,   *
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL    *
 * THE AUTHORS, CONTRIBUTORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,    *
 * DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR      *
 * OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE  *
 * USE OR OTHER DEALINGS IN THE SOFTWARE.                                     *
 * -------------------------------------------------------------------------- */

#include "openmm/DrudeNoseHooverIntegrator.h"
#include "openmm/Context.h"
#include "openmm/OpenMMException.h"
#include "openmm/internal/ContextImpl.h"
#include "openmm/DrudeKernels.h"
#include "openmm/DrudeForce.h"
#include "openmm/CMMotionRemover.h"
#include "openmm/kernels.h"
#include <ctime>
#include <iostream>
#include <string>
#include <set>

using namespace OpenMM;
using std::string;
using std::vector;

DrudeNoseHooverIntegrator::DrudeNoseHooverIntegrator(double stepSize, System &system, double temperature, double drudeTemperature,
                                                     double collisionFrequency, double drudeCollisionFrequency,
                                                     int chainLength, int numMTS, int numYoshidaSuzuki) :
    NoseHooverIntegrator(stepSize) {
    const DrudeForce* drudeForce = NULL;
    for (int i = 0; i < system.getNumForces(); i++)
        if (dynamic_cast<const DrudeForce*>(&system.getForce(i)) != NULL) {
            if (drudeForce == NULL)
                drudeForce = dynamic_cast<const DrudeForce*>(&system.getForce(i));
            else
                throw OpenMMException("The System contains multiple DrudeForces");
        }
    if (drudeForce == NULL)
        throw OpenMMException("The System does not contain a DrudeForce");
    std::set<int> realParticlesSet;
    vector<int> realParticles, drudeParticles, drudeParents;
    vector<std::pair<int, int>> drudePairs;
    for (int i = 0; i < system.getNumParticles(); i++) {
        if (system.getParticleMass(i) > 0.0) realParticlesSet.insert(i);
    }
    for (int i = 0; i < drudeForce->getNumParticles(); i++) {
        int p, p1, p2, p3, p4;
        double charge, polarizability, aniso12, aniso34;
        drudeForce->getParticleParameters(i, p, p1, p2, p3, p4, charge, polarizability, aniso12, aniso34);
        realParticlesSet.erase(p);
        realParticlesSet.erase(p1);
        drudePairs.push_back({p,p1});
    }
    for(const auto &p : realParticlesSet) realParticles.push_back(p);

    addSubsystemThermostat(system, realParticles, drudePairs, temperature, drudeTemperature, collisionFrequency,
                           drudeCollisionFrequency, chainLength, numMTS, numYoshidaSuzuki);
}

DrudeNoseHooverIntegrator::~DrudeNoseHooverIntegrator() { }


void DrudeNoseHooverIntegrator::initialize(ContextImpl& contextRef) {
    NoseHooverIntegrator::initialize(contextRef);
    if (owner != NULL && &contextRef.getOwner() != owner)
        throw OpenMMException("This Integrator is already bound to a context");
    const DrudeForce* drudeForce = NULL;
    const System& system = contextRef.getSystem();
    bool hasCMMotionRemover = false;
    for (int i = 0; i < system.getNumForces(); i++){
        if (dynamic_cast<const DrudeForce*>(&system.getForce(i)) != NULL) {
            if (drudeForce == NULL)
                drudeForce = dynamic_cast<const DrudeForce*>(&system.getForce(i));
            else
                throw OpenMMException("The System contains multiple DrudeForces");
        }
        if (dynamic_cast<const CMMotionRemover*>(&system.getForce(i))) {
            hasCMMotionRemover = true;
        }
    }
    if (drudeForce == NULL)
        throw OpenMMException("The System does not contain a DrudeForce");
    if (!hasCMMotionRemover) {
        std::cout << "Warning: Did not find a center-of-mass motion remover in the system. "
                     "This is problematic when using Drude." << std::endl;
    }
    context = &contextRef;
    owner = &contextRef.getOwner();
}

double DrudeNoseHooverIntegrator::computeDrudeKineticEnergy() {
    double kE = 0.0;
    for (const auto &nhc: noseHooverChains){
        if (nhc.getThermostatedPairs().size() != 0) {
            kE += nhcKernel.getAs<NoseHooverChainKernel>().computeMaskedKineticEnergy(*context, nhc, true).second;
        }
    }
    return kE;
}

double DrudeNoseHooverIntegrator::computeTotalKineticEnergy() {
    return vvKernel.getAs<IntegrateVelocityVerletStepKernel>().computeKineticEnergy(*context, *this);
}
