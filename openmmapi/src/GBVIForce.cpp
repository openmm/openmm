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
#include "openmm/GBVIForce.h"
#include "openmm/internal/AssertionUtilities.h"
#include "openmm/internal/GBVIForceImpl.h"
#include <sstream>

using namespace OpenMM;

GBVIForce::GBVIForce() : nonbondedMethod(NoCutoff), cutoffDistance(1.0), solventDielectric(78.3), soluteDielectric(1.0),
                         scalingMethod(QuinticSpline), quinticLowerLimitFactor(0.8), quinticUpperBornRadiusLimit(5.0) {
}

int GBVIForce::addParticle(double charge, double radius, double gamma) {
    particles.push_back(ParticleInfo(charge, radius, gamma));
    return particles.size()-1;
}

void GBVIForce::getParticleParameters(int index, double& charge, double& radius, double& gamma) const {
    ASSERT_VALID_INDEX(index, particles);
    charge = particles[index].charge;
    radius = particles[index].radius;
    gamma  = particles[index].gamma;
}

void GBVIForce::setParticleParameters(int index, double charge, double radius, double gamma) {
    ASSERT_VALID_INDEX(index, particles);
    particles[index].charge = charge;
    particles[index].radius = radius;
    particles[index].gamma  = gamma;
}

GBVIForce::NonbondedMethod GBVIForce::getNonbondedMethod() const {
    return nonbondedMethod;
}

void GBVIForce::setNonbondedMethod(NonbondedMethod method) {
    nonbondedMethod = method;
}

double GBVIForce::getCutoffDistance() const {
    return cutoffDistance;
}

void GBVIForce::setCutoffDistance(double distance) {
    cutoffDistance = distance;
}

GBVIForce::BornRadiusScalingMethod GBVIForce::getBornRadiusScalingMethod( void ) const {
    return scalingMethod;
}

void GBVIForce::setBornRadiusScalingMethod(BornRadiusScalingMethod method) {
    scalingMethod = method;
}

double GBVIForce::getQuinticLowerLimitFactor( void ) const {
    return quinticLowerLimitFactor;
}

void GBVIForce::setQuinticLowerLimitFactor(double inputQuinticLowerLimitFactor ){
    quinticLowerLimitFactor = inputQuinticLowerLimitFactor;
}

double GBVIForce::getQuinticUpperBornRadiusLimit( void ) const {
    return quinticUpperBornRadiusLimit;
}

void GBVIForce::setQuinticUpperBornRadiusLimit(double inputQuinticUpperBornRadiusLimit){
    quinticUpperBornRadiusLimit = inputQuinticUpperBornRadiusLimit;
}

int GBVIForce::addBond(int particle1, int particle2, double bondLength) {
    bonds.push_back(BondInfo(particle1, particle2, bondLength));
    return bonds.size()-1;
}

void GBVIForce::setBondParameters( int index, int particle1, int particle2, double bondLength) {
    ASSERT_VALID_INDEX(index, bonds);
    bonds[index].particle1  = particle1;
    bonds[index].particle2  = particle2;
    bonds[index].bondLength = bondLength;
}

int GBVIForce::getNumBonds( void ) const {
   return (int) bonds.size();
}

void GBVIForce::getBondParameters(int index, int& bondIndex1, int& bondIndex2, double& bondLength) const {
    ASSERT_VALID_INDEX(index, bonds);
    bondIndex1 = bonds[index].particle1;
    bondIndex2 = bonds[index].particle2;
    bondLength = bonds[index].bondLength;
}

ForceImpl* GBVIForce::createImpl() const {
    return new GBVIForceImpl(*this);
}
