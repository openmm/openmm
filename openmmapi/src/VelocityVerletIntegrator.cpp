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

#include "openmm/VelocityVerletIntegrator.h"
#include "openmm/Context.h"
#include "openmm/Force.h"
#include "openmm/System.h"
#include "openmm/NoseHooverChain.h"
#include "openmm/OpenMMException.h"
#include "openmm/CMMotionRemover.h"
#include "openmm/internal/ContextImpl.h"
#include "openmm/kernels.h"
#include <iostream>
#include <string>
#include <algorithm>

using namespace OpenMM;
using std::string;
using std::vector;

VelocityVerletIntegrator::VelocityVerletIntegrator(double stepSize):
    forcesAreValid(false)
{
    setStepSize(stepSize);
    setConstraintTolerance(1e-5);
}

VelocityVerletIntegrator::~VelocityVerletIntegrator() {}

double VelocityVerletIntegrator::propagateChain(double kineticEnergy, int chainID) {
    return nhcKernel.getAs<NoseHooverChainKernel>().propagateChain(*context, noseHooverChains.at(chainID), kineticEnergy, getStepSize());
}

int VelocityVerletIntegrator::addNoseHooverChainThermostat(System& system, double temperature, double collisionFrequency,
                                                           int chainLength, int numMTS, int numYoshidaSuzuki) {
    int numForces = system.getNumForces();
    std::vector<int> thermostatedParticles;
    std::vector<int> parentParticles;
    for(int particle = 0; particle < system.getNumParticles(); ++particle) {
        double mass = system.getParticleMass(particle);
        if ( (mass > 0) && (mass < 0.8) ){
            std::cout << "Warning: Found particles with mass between 0.0 and 0.8 dalton. Did you mean to make a DrudeVelocityVerletIntegrator instead? "
                         "The thermostat you are about to use will not treat these particles as Drude particles!" << std::endl;
        }
        if(system.getParticleMass(particle) > 0) {
            thermostatedParticles.push_back(particle);
        }
    }

    return addMaskedNoseHooverChainThermostat(system, thermostatedParticles, parentParticles, temperature, collisionFrequency, chainLength, numMTS, numYoshidaSuzuki);
}

int VelocityVerletIntegrator::addMaskedNoseHooverChainThermostat(System& system, const std::vector<int>& thermostatedParticles, const std::vector<int>& parentParticles,
                                         double temperature, double collisionFrequency,
                                         int chainLength, int numMTS, int numYoshidaSuzuki){
    const auto & mask = thermostatedParticles; // just an alias
    if (context) {
        throw OpenMMException("Nose-Hoover chains cannot be added after binding this integrator to a context.");
    }
    if ((parentParticles.size() != mask.size()) && (parentParticles.size() != 0)) {
        throw OpenMMException("The number of parent particles has to be either equal to the number of thermostated particles (to thermostat relative motion)"
                              " or zero (to thermostat absolute motion).");
    }
    // figure out the number of DOFs
    int nDOF = 3*mask.size();
    for (int constraintNum = 0; constraintNum < system.getNumConstraints(); constraintNum++) {
        int particle1, particle2;
        double distance;
        system.getConstraintParameters(constraintNum, particle1, particle2, distance);
        bool particle1_in_mask = (std::find(mask.begin(), mask.end(), particle1) != mask.end());
        bool particle2_in_mask = (std::find(mask.begin(), mask.end(), particle2) != mask.end());
        if ((system.getParticleMass(particle1) > 0) && (system.getParticleMass(particle2) > 0)){
            if ((particle1_in_mask && !particle2_in_mask) || (!particle1_in_mask && particle2_in_mask)){
                throw OpenMMException("Cannot add only one of particles " + std::to_string(particle1) + " and " + std::to_string(particle2)
                                        + " to NoseHooverChain, because they are connected by a constraint.");
            }
            if (particle1_in_mask && particle2_in_mask){
                nDOF -= 1;
            }
        }
    }
    int numForces = system.getNumForces();
    if (parentParticles.size() == 0){ // remove 3 degrees of freedom from thermostats that act on absolute motions
        for (int forceNum = 0; forceNum < numForces; ++forceNum) {
            if (dynamic_cast<CMMotionRemover*>(&system.getForce(forceNum))) nDOF -= 3;
        }
    }

    // make sure that thermostats do not overlap 
    for (auto &other_nhc: noseHooverChains) {
        for (auto &particle: mask){
            bool isParticleInOtherChain = (std::find(other_nhc.getThermostatedAtoms().begin(), 
                                                    other_nhc.getThermostatedAtoms().end(),
                                                    particle) == other_nhc.getThermostatedAtoms().begin());
            if (isParticleInOtherChain){
                throw OpenMMException("Found particle " + std::to_string(particle) + "in a different NoseHooverChain, "
                                      "but particles can only be thermostated by one thermostat.");
            }
        }
    }

    // create and add new chain 
    int chainID = noseHooverChains.size();
    auto nhc = NoseHooverChain(temperature, collisionFrequency, nDOF, chainLength,
                                    numMTS, numYoshidaSuzuki, chainID,
                                    thermostatedParticles, parentParticles);
    noseHooverChains.push_back(nhc);
    return chainID;
}

double VelocityVerletIntegrator::getTemperature(int chainID) const {
    if (chainID >= noseHooverChains.size()) {
        throw OpenMMException("Cannot get temperature for chainID " + std::to_string(chainID)
                + ". Only " + std::to_string(noseHooverChains.size()) + " have been registered with this integrator."
        );
    }
    return noseHooverChains.at(chainID).getDefaultTemperature();
}

void VelocityVerletIntegrator::setTemperature(double temperature, int chainID){
    if (chainID >= noseHooverChains.size()) {
        throw OpenMMException("Cannot set temperature for chainID " + std::to_string(chainID)
                + ". Only " + std::to_string(noseHooverChains.size()) + " have been registered with this integrator."
        );
    }
    noseHooverChains.at(chainID).setDefaultTemperature(temperature);

}

double VelocityVerletIntegrator::getCollisionFrequency(int chainID) const {
    if (chainID >= noseHooverChains.size()) {
        throw OpenMMException("Cannot get collision frequency for chainID " + std::to_string(chainID)
                + ". Only " + std::to_string(noseHooverChains.size()) + " have been registered with this integrator."
        );
    }
    return noseHooverChains.at(chainID).getDefaultCollisionFrequency();
} 
void VelocityVerletIntegrator::setCollisionFrequency(double frequency, int chainID){
    if (chainID >= noseHooverChains.size()) {
        throw OpenMMException("Cannot set collision frequency for chainID " + std::to_string(chainID)
                + ". Only " + std::to_string(noseHooverChains.size()) + " have been registered with this integrator."
        );
    }
    noseHooverChains.at(chainID).setDefaultCollisionFrequency(frequency);
}

double VelocityVerletIntegrator::computeKineticEnergy() {
    double kE = 0.0;
    if(noseHooverChains.size()) {
        for (const auto &nhc: noseHooverChains){
            if (nhc.getParentAtoms().size() == 0) {
                kE += nhcKernel.getAs<NoseHooverChainKernel>().computeMaskedKineticEnergy(*context, nhc, true);
            }
        }
    } else {
        kE = vvKernel.getAs<IntegrateVelocityVerletStepKernel>().computeKineticEnergy(*context, *this);
    }
    return kE;
}

double VelocityVerletIntegrator::computeHeatBathEnergy() {
    double energy = 0;
    for(auto &nhc : noseHooverChains) {
        energy += nhcKernel.getAs<NoseHooverChainKernel>().computeHeatBathEnergy(*context, nhc);
    }
    return energy;
}

void VelocityVerletIntegrator::initialize(ContextImpl& contextRef) {
    if (owner != NULL && &contextRef.getOwner() != owner)
        throw OpenMMException("This Integrator is already bound to a context");
    context = &contextRef;
    owner = &contextRef.getOwner();
    vvKernel = context->getPlatform().createKernel(IntegrateVelocityVerletStepKernel::Name(), contextRef);
    vvKernel.getAs<IntegrateVelocityVerletStepKernel>().initialize(contextRef.getSystem(), *this);
    nhcKernel = context->getPlatform().createKernel(NoseHooverChainKernel::Name(), contextRef);
    nhcKernel.getAs<NoseHooverChainKernel>().initialize();
    forcesAreValid = false;
}

void VelocityVerletIntegrator::cleanup() {
    vvKernel = Kernel();
    nhcKernel = Kernel();
}

vector<string> VelocityVerletIntegrator::getKernelNames() {
    std::vector<std::string> names;
    names.push_back(NoseHooverChainKernel::Name());
    names.push_back(IntegrateVelocityVerletStepKernel::Name());
    return names;
}

void VelocityVerletIntegrator::step(int steps) {
    if (context == NULL)
        throw OpenMMException("This Integrator is not bound to a context!");
    double scale, kineticEnergy;
    for (int i = 0; i < steps; ++i) {
        context->updateContextState();
        for(auto &nhc : noseHooverChains) {
            kineticEnergy = nhcKernel.getAs<NoseHooverChainKernel>().computeMaskedKineticEnergy(*context, nhc, false);
            scale = nhcKernel.getAs<NoseHooverChainKernel>().propagateChain(*context, nhc, kineticEnergy, getStepSize());
            nhcKernel.getAs<NoseHooverChainKernel>().scaleVelocities(*context, nhc, scale);
        }
        vvKernel.getAs<IntegrateVelocityVerletStepKernel>().execute(*context, *this, forcesAreValid);
        for(auto &nhc : noseHooverChains) {
            kineticEnergy = nhcKernel.getAs<NoseHooverChainKernel>().computeMaskedKineticEnergy(*context, nhc, false);
            scale = nhcKernel.getAs<NoseHooverChainKernel>().propagateChain(*context, nhc, kineticEnergy, getStepSize());
            nhcKernel.getAs<NoseHooverChainKernel>().scaleVelocities(*context, nhc, scale);
        }
    }
}
