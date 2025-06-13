/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2025 Stanford University and the Authors.           *
 * Authors: Peter Eastman                                                     *
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

#include "openmm/DPDIntegrator.h"
#include "openmm/Context.h"
#include "openmm/OpenMMException.h"
#include "openmm/internal/ContextImpl.h"
#include "openmm/kernels.h"
#include <string>

using namespace OpenMM;
using namespace std;

DPDIntegrator::DPDIntegrator(double temperature, double defaultFriction, double defaultCutoff, double stepSize) {
    setTemperature(temperature);
    setDefaultFriction(defaultFriction);
    setDefaultCutoff(defaultCutoff);
    setStepSize(stepSize);
    setConstraintTolerance(1e-5);
    setRandomNumberSeed(0);
}

void DPDIntegrator::initialize(ContextImpl& contextRef) {
    if (owner != NULL && &contextRef.getOwner() != owner)
        throw OpenMMException("This Integrator is already bound to a context");
    context = &contextRef;
    owner = &contextRef.getOwner();
    kernel = context->getPlatform().createKernel(IntegrateDPDStepKernel::Name(), contextRef);
    kernel.getAs<IntegrateDPDStepKernel>().initialize(contextRef.getSystem(), *this);
}

void DPDIntegrator::setTemperature(double temp) {
    if (temp < 0)
        throw OpenMMException("Temperature cannot be negative");
    temperature = temp;
}

void DPDIntegrator::setDefaultFriction(double friction) {
    if (friction < 0)
        throw OpenMMException("Friction cannot be negative");
    defaultFriction = friction;
}

void DPDIntegrator::setDefaultCutoff(double cutoff) {
    if (cutoff <= 0)
        throw OpenMMException("Cutoff must be positive");
    defaultCutoff = cutoff;
}

int DPDIntegrator::getParticleType(int index) const {
    if (particleType.find(index) == particleType.end())
        return 0;
    return particleType.at(index);
    
}

void DPDIntegrator::setParticleType(int index, int type) {
    particleType[index] = type;
}

const map<int, int>& DPDIntegrator::getParticleTypes() const {
    return particleType;
}

int DPDIntegrator::addTypePair(int type1, int type2, double friction, double cutoff) {
    pairs.push_back(TypePairInfo(type1, type2, friction, cutoff));
    return pairs.size()-1;
}

void DPDIntegrator::getTypePairParameters(int pairIndex, int& type1, int& type2, double& friction, double& cutoff) const {
    ASSERT_VALID_INDEX(pairIndex, pairs);
    type1 = pairs[pairIndex].type1;
    type2 = pairs[pairIndex].type2;
    friction = pairs[pairIndex].friction;
    cutoff = pairs[pairIndex].cutoff;
}

void DPDIntegrator::setTypePairParameters(int pairIndex, int type1, int type2, double friction, double cutoff) {
    ASSERT_VALID_INDEX(pairIndex, pairs);
    pairs[pairIndex].type1 = type1;
    pairs[pairIndex].type2 = type2;
    pairs[pairIndex].friction = friction;
    pairs[pairIndex].cutoff = cutoff;
}

void DPDIntegrator::cleanup() {
    kernel = Kernel();
}

vector<string> DPDIntegrator::getKernelNames() {
    std::vector<std::string> names;
    names.push_back(IntegrateDPDStepKernel::Name());
    return names;
}

double DPDIntegrator::computeKineticEnergy() {
    return kernel.getAs<IntegrateDPDStepKernel>().computeKineticEnergy(*context, *this);
}

bool DPDIntegrator::kineticEnergyRequiresForce() const {
    return false;
}

void DPDIntegrator::step(int steps) {
    if (context == NULL)
        throw OpenMMException("This Integrator is not bound to a context!");
    for (int i = 0; i < steps; ++i) {
        context->updateContextState();
        context->calcForcesAndEnergy(true, false, getIntegrationForceGroups());
        kernel.getAs<IntegrateDPDStepKernel>().execute(*context, *this);
    }
}
