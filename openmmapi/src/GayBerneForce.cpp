/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2008-2016 Stanford University and the Authors.      *
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
#include "openmm/GayBerneForce.h"
#include "openmm/internal/AssertionUtilities.h"
#include "openmm/internal/GayBerneForceImpl.h"
#include <cmath>
#include <map>
#include <sstream>
#include <utility>

using namespace OpenMM;
using std::map;
using std::pair;
using std::set;
using std::string;
using std::stringstream;
using std::vector;

GayBerneForce::GayBerneForce() : nonbondedMethod(NoCutoff), cutoffDistance(1.0), switchingDistance(-1.0), useSwitchingFunction(false) {
}

GayBerneForce::NonbondedMethod GayBerneForce::getNonbondedMethod() const {
    return nonbondedMethod;
}

void GayBerneForce::setNonbondedMethod(NonbondedMethod method) {
    if (method < 0 || method > 2)
        throw OpenMMException("GayBerneForce: Illegal value for nonbonded method");
    nonbondedMethod = method;
}

double GayBerneForce::getCutoffDistance() const {
    return cutoffDistance;
}

void GayBerneForce::setCutoffDistance(double distance) {
    cutoffDistance = distance;
}

bool GayBerneForce::getUseSwitchingFunction() const {
    return useSwitchingFunction;
}

void GayBerneForce::setUseSwitchingFunction(bool use) {
    useSwitchingFunction = use;
}

double GayBerneForce::getSwitchingDistance() const {
    return switchingDistance;
}

void GayBerneForce::setSwitchingDistance(double distance) {
    switchingDistance = distance;
}

int GayBerneForce::addParticle(double sigma, double epsilon, int xparticle, int yparticle, double sx, double sy, double sz, double ex, double ey, double ez) {
    if (yparticle == -1 && (sy != sz || ey != ez))
        throw OpenMMException("GayBerneForce: yparticle is -1 for a particle that is not axially symmetric");
    if (xparticle == -1 && (sx != sz || ex != ez))
        throw OpenMMException("GayBerneForce: xparticle is -1 for a particle that is not spherical");
    if (xparticle == -1 && yparticle != -1)
        throw OpenMMException("GayBerneForce: xparticle cannot be -1 if yparticle is not also -1");
    particles.push_back(ParticleInfo(sigma, epsilon, xparticle, yparticle, sx, sy, sz, ex, ey, ez));
    return particles.size()-1;
}

void GayBerneForce::getParticleParameters(int index, double& sigma, double& epsilon, int& xparticle, int& yparticle, double& sx, double& sy, double& sz, double& ex, double& ey, double& ez) const {
    ASSERT_VALID_INDEX(index, particles);
    sigma = particles[index].sigma;
    epsilon = particles[index].epsilon;
    xparticle = particles[index].xparticle;
    yparticle = particles[index].yparticle;
    sx = particles[index].sx;
    sy = particles[index].sy;
    sz = particles[index].sz;
    ex = particles[index].ex;
    ey = particles[index].ey;
    ez = particles[index].ez;
}

void GayBerneForce::setParticleParameters(int index, double sigma, double epsilon, int xparticle, int yparticle, double sx, double sy, double sz, double ex, double ey, double ez) {
    ASSERT_VALID_INDEX(index, particles);
    if (yparticle == -1 && (sy != sz || ey != ez))
        throw OpenMMException("GayBerneForce: yparticle is -1 for a particle that is not axially symmetric");
    if (xparticle == -1 && (sx != sz || ex != ez))
        throw OpenMMException("GayBerneForce: xparticle is -1 for a particle that is not spherical");
    if (xparticle == -1 && yparticle != -1)
        throw OpenMMException("GayBerneForce: xparticle cannot be -1 if yparticle is not also -1");
    particles[index].sigma = sigma;
    particles[index].epsilon = epsilon;
    particles[index].xparticle = xparticle;
    particles[index].yparticle = yparticle;
    particles[index].sx = sx;
    particles[index].sy = sy;
    particles[index].sz = sz;
    particles[index].ex = ex;
    particles[index].ey = ey;
    particles[index].ez = ez;
}

int GayBerneForce::addException(int particle1, int particle2, double sigma, double epsilon, bool replace) {
    map<pair<int, int>, int>::iterator iter = exceptionMap.find(pair<int, int>(particle1, particle2));
    int newIndex;
    if (iter == exceptionMap.end())
        iter = exceptionMap.find(pair<int, int>(particle2, particle1));
    if (iter != exceptionMap.end()) {
        if (!replace) {
            stringstream msg;
            msg << "GayBerneForce: There is already an exception for particles ";
            msg << particle1;
            msg << " and ";
            msg << particle2;
            throw OpenMMException(msg.str());
        }
        exceptions[iter->second] = ExceptionInfo(particle1, particle2, sigma, epsilon);
        newIndex = iter->second;
        exceptionMap.erase(iter->first);
    }
    else {
        exceptions.push_back(ExceptionInfo(particle1, particle2, sigma, epsilon));
        newIndex = exceptions.size()-1;
    }
    exceptionMap[pair<int, int>(particle1, particle2)] = newIndex;
    return newIndex;
}

void GayBerneForce::getExceptionParameters(int index, int& particle1, int& particle2, double& sigma, double& epsilon) const {
    ASSERT_VALID_INDEX(index, exceptions);
    particle1 = exceptions[index].particle1;
    particle2 = exceptions[index].particle2;
    sigma = exceptions[index].sigma;
    epsilon = exceptions[index].epsilon;
}

void GayBerneForce::setExceptionParameters(int index, int particle1, int particle2, double sigma, double epsilon) {
    ASSERT_VALID_INDEX(index, exceptions);
    exceptions[index].particle1 = particle1;
    exceptions[index].particle2 = particle2;
    exceptions[index].sigma = sigma;
    exceptions[index].epsilon = epsilon;
}

ForceImpl* GayBerneForce::createImpl() const {
    return new GayBerneForceImpl(*this);
}

void GayBerneForce::updateParametersInContext(Context& context) {
    dynamic_cast<GayBerneForceImpl&>(getImplInContext(context)).updateParametersInContext(getContextImpl(context));
}
