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
#include "openmm/GBVISoftcoreForce.h"
#include "openmm/internal/GBVISoftcoreForceImpl.h"
#include <sstream>

using namespace OpenMM;

GBVISoftcoreForce::GBVISoftcoreForce() :  nonbondedMethod(NoCutoff), cutoffDistance(1.0), solventDielectric(78.3), soluteDielectric(1.0),
               scalingMethod(QuinticSpline), alpha(1.0), beta(0.8), gamma(4.85), quinticLowerLimitFactor(0.8), quinticUpperBornRadiusLimit(5.0) {

}

void GBVISoftcoreForce::addParticle(double charge, double radius, double gamma, double bornRadiusScaleFactor) {
    particles.push_back(ParticleInfo(charge, radius, gamma, bornRadiusScaleFactor));
}

void GBVISoftcoreForce::getParticleParameters(int index, double& charge, double& radius, double& gamma ) const {
    charge = particles[index].charge;
    radius = particles[index].radius;
    gamma  = particles[index].gamma;
}

void GBVISoftcoreForce::getParticleParameters(int index, double& charge, double& radius, double& gamma, double& bornRadiusScaleFactor ) const {
    charge                = particles[index].charge;
    radius                = particles[index].radius;
    gamma                 = particles[index].gamma;
    bornRadiusScaleFactor = particles[index].bornRadiusScaleFactor;
}

void GBVISoftcoreForce::setParticleParameters(int index, double charge, double radius, double gamma, double bornRadiusScaleFactor) {
    particles[index].charge                = charge;
    particles[index].radius                = radius;
    particles[index].gamma                 = gamma;
    particles[index].bornRadiusScaleFactor = bornRadiusScaleFactor;
}

GBVISoftcoreForce::NonbondedSoftcoreMethod GBVISoftcoreForce::getNonbondedMethod() const {
    return nonbondedMethod;
}

void GBVISoftcoreForce::setNonbondedMethod(NonbondedSoftcoreMethod method) {
    nonbondedMethod = method;
}

double GBVISoftcoreForce::getCutoffDistance() const {
    return cutoffDistance;
}

void GBVISoftcoreForce::setCutoffDistance(double distance) {
    cutoffDistance = distance;
}

GBVISoftcoreForce::BornRadiusScalingSoftcoreMethod GBVISoftcoreForce::getBornRadiusScalingMethod( void ) const {
    return scalingMethod;
}

void GBVISoftcoreForce::setBornRadiusScalingMethod(BornRadiusScalingSoftcoreMethod method) {
    scalingMethod = method;
}

double GBVISoftcoreForce::getQuinticLowerLimitFactor( void ) const {
    return quinticLowerLimitFactor;
}

void GBVISoftcoreForce::setQuinticLowerLimitFactor(double inputQuinticLowerLimitFactor ){
    quinticLowerLimitFactor = inputQuinticLowerLimitFactor;
}

double GBVISoftcoreForce::getQuinticUpperBornRadiusLimit( void ) const {
    return quinticUpperBornRadiusLimit;
}

void GBVISoftcoreForce::setQuinticUpperBornRadiusLimit(double inputQuinticUpperBornRadiusLimit){
    quinticUpperBornRadiusLimit = inputQuinticUpperBornRadiusLimit;
}

int GBVISoftcoreForce::addBond(int particle1, int particle2, double bondLength) {
    bonds.push_back(BondInfo(particle1, particle2, bondLength));
    return bonds.size()-1;
}

void GBVISoftcoreForce::setBondParameters( int index, int particle1, int particle2, double bondLength) {
    bonds[index].particle1  = particle1;
    bonds[index].particle2  = particle2;
    bonds[index].bondLength = bondLength;
}

int GBVISoftcoreForce::getNumBonds( void ) const {
   return (int) bonds.size();
}

void GBVISoftcoreForce::getBondParameters(int index, int& bondIndex1, int& bondIndex2, double& bondLength) const {
    bondIndex1 = bonds[index].particle1;
    bondIndex2 = bonds[index].particle2;
    bondLength = bonds[index].bondLength;
}

ForceImpl* GBVISoftcoreForce::createImpl() {
    return new GBVISoftcoreForceImpl(*this);
}
