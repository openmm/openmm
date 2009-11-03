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
#include "openmm/GBSAOBCSoftcoreForce.h"
#include "openmm/internal/GBSAOBCSoftcoreForceImpl.h"

using namespace OpenMM;

GBSAOBCSoftcoreForce::GBSAOBCSoftcoreForce() : solventDielectric(78.3), soluteDielectric(1.0), nonPolarPrefactor( 0.0054 ){
}

int GBSAOBCSoftcoreForce::addParticle(double charge, double radius, double scalingFactor, double nonPolarScalingFactor) {
    particles.push_back(ParticleInfo(charge, radius, scalingFactor, nonPolarScalingFactor));
    return particles.size()-1;
}

void GBSAOBCSoftcoreForce::getParticleParameters(int index, double& charge, double& radius, double& scalingFactor) const {
    charge                = particles[index].charge;
    radius                = particles[index].radius;
    scalingFactor         = particles[index].scalingFactor;
}

void GBSAOBCSoftcoreForce::getParticleParameters(int index, double& charge, double& radius, double& scalingFactor, double& nonPolarScalingFactor) const {
    charge                = particles[index].charge;
    radius                = particles[index].radius;
    scalingFactor         = particles[index].scalingFactor;
    nonPolarScalingFactor = particles[index].nonPolarScalingFactor;
}

void GBSAOBCSoftcoreForce::setParticleParameters(int index, double charge, double radius, double scalingFactor, double nonPolarScalingFactor) {
    particles[index].charge                 = charge;
    particles[index].radius                 = radius;
    particles[index].scalingFactor          = scalingFactor;
    particles[index].nonPolarScalingFactor  = nonPolarScalingFactor;
}

ForceImpl* GBSAOBCSoftcoreForce::createImpl() {
    return new GBSAOBCSoftcoreForceImpl(*this);
}
