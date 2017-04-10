/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2008-2009 Stanford University and the Authors.      *
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

#include "openmm/Force.h"
#include "openmm/OpenMMException.h"
#include "openmm/GBSAOBCForce.h"
#include "openmm/internal/AssertionUtilities.h"
#include "openmm/internal/GBSAOBCForceImpl.h"

using namespace OpenMM;

GBSAOBCForce::GBSAOBCForce() : nonbondedMethod(NoCutoff), cutoffDistance(1.0), solventDielectric(78.3), soluteDielectric(1.0), surfaceAreaEnergy(2.25936) {
}

int GBSAOBCForce::addParticle(double charge, double radius, double scalingFactor) {
    particles.push_back(ParticleInfo(charge, radius, scalingFactor));
    return particles.size()-1;
}

void GBSAOBCForce::getParticleParameters(int index, double& charge, double& radius, double& scalingFactor) const {
    ASSERT_VALID_INDEX(index, particles);
    charge = particles[index].charge;
    radius = particles[index].radius;
    scalingFactor = particles[index].scalingFactor;
}

void GBSAOBCForce::setParticleParameters(int index, double charge, double radius, double scalingFactor) {
    ASSERT_VALID_INDEX(index, particles);
    particles[index].charge = charge;
    particles[index].radius = radius;
    particles[index].scalingFactor = scalingFactor;
}

GBSAOBCForce::NonbondedMethod GBSAOBCForce::getNonbondedMethod() const {
    return nonbondedMethod;
}

void GBSAOBCForce::setNonbondedMethod(NonbondedMethod method) {
    if (method < 0 || method > 2)
        throw OpenMMException("GBSAOBCForce: Illegal value for nonbonded method");
    nonbondedMethod = method;
}

double GBSAOBCForce::getCutoffDistance() const {
    return cutoffDistance;
}

void GBSAOBCForce::setCutoffDistance(double distance) {
    cutoffDistance = distance;
}

ForceImpl* GBSAOBCForce::createImpl() const {
    return new GBSAOBCForceImpl(*this);
}

void GBSAOBCForce::updateParametersInContext(Context& context) {
    dynamic_cast<GBSAOBCForceImpl&>(getImplInContext(context)).updateParametersInContext(getContextImpl(context));
}
