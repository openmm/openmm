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
#include "openmm/CustomNonbondedForce.h"
#include "openmm/internal/CustomNonbondedForceImpl.h"
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

CustomNonbondedForce::CustomNonbondedForce(const string& energy) : energyExpression(energy), nonbondedMethod(NoCutoff), cutoffDistance(1.0) {
    periodicBoxVectors[0] = Vec3(2, 0, 0);
    periodicBoxVectors[1] = Vec3(0, 2, 0);
    periodicBoxVectors[2] = Vec3(0, 0, 2);
}

const string& CustomNonbondedForce::getEnergyFunction() const {
    return energyExpression;
}

void CustomNonbondedForce::setEnergyFunction(const std::string& energy) {
    energyExpression = energy;
}

CustomNonbondedForce::NonbondedMethod CustomNonbondedForce::getNonbondedMethod() const {
    return nonbondedMethod;
}

void CustomNonbondedForce::setNonbondedMethod(NonbondedMethod method) {
    nonbondedMethod = method;
}

double CustomNonbondedForce::getCutoffDistance() const {
    return cutoffDistance;
}

void CustomNonbondedForce::setCutoffDistance(double distance) {
    cutoffDistance = distance;
}

void CustomNonbondedForce::getPeriodicBoxVectors(Vec3& a, Vec3& b, Vec3& c) const {
    a = periodicBoxVectors[0];
    b = periodicBoxVectors[1];
    c = periodicBoxVectors[2];
}

void CustomNonbondedForce::setPeriodicBoxVectors(Vec3 a, Vec3 b, Vec3 c) {
    if (a[1] != 0.0 || a[2] != 0.0)
        throw OpenMMException("First periodic box vector must be parallel to x.");
    if (b[0] != 0.0 || b[2] != 0.0)
        throw OpenMMException("Second periodic box vector must be parallel to y.");
    if (c[0] != 0.0 || c[1] != 0.0)
        throw OpenMMException("Third periodic box vector must be parallel to z.");
    periodicBoxVectors[0] = a;
    periodicBoxVectors[1] = b;
    periodicBoxVectors[2] = c;
}

int CustomNonbondedForce::addParameter(const string& combiningRule) {
    combiningRules.push_back(combiningRule);
    return combiningRules.size()-1;
}

const string& CustomNonbondedForce::getParameterCombiningRule(int index) const {
    return combiningRules[index];
}

void CustomNonbondedForce::setParameterCombiningRule(int index, const string& combiningRule) {
    combiningRules[index] = combiningRule;
}

int CustomNonbondedForce::addParticle(const vector<double>& parameters) {
    particles.push_back(ParticleInfo(parameters));
    return particles.size()-1;
}

void CustomNonbondedForce::getParticleParameters(int index, std::vector<double>& parameters) const {
    parameters = particles[index].parameters;
}

void CustomNonbondedForce::setParticleParameters(int index, const vector<double>& parameters) {
    particles[index].parameters = parameters;
}

int CustomNonbondedForce::addException(int particle1, int particle2, const vector<double>& parameters, bool replace) {
    map<pair<int, int>, int>::iterator iter = exceptionMap.find(pair<int, int>(particle1, particle2));
    int newIndex;
    if (iter == exceptionMap.end())
        iter = exceptionMap.find(pair<int, int>(particle2, particle1));
    if (iter != exceptionMap.end()) {
        if (!replace) {
            stringstream msg;
            msg << "CustomNonbondedForce: There is already an exception for particles ";
            msg << particle1;
            msg << " and ";
            msg << particle2;
            throw OpenMMException(msg.str());
        }
        exceptions[iter->second] = ExceptionInfo(particle1, particle2, parameters);
        newIndex = iter->second;
        exceptionMap.erase(iter->first);
    }
    else {
        exceptions.push_back(ExceptionInfo(particle1, particle2, parameters));
        newIndex = exceptions.size()-1;
    }
    exceptionMap[pair<int, int>(particle1, particle2)] = newIndex;
    return newIndex;
}
void CustomNonbondedForce::getExceptionParameters(int index, int& particle1, int& particle2, vector<double>& parameters) const {
    particle1 = exceptions[index].particle1;
    particle2 = exceptions[index].particle2;
    parameters = exceptions[index].parameters;
}

void CustomNonbondedForce::setExceptionParameters(int index, int particle1, int particle2, const vector<double>& parameters) {
    exceptions[index].particle1 = particle1;
    exceptions[index].particle2 = particle2;
    exceptions[index].parameters = parameters;
}

ForceImpl* CustomNonbondedForce::createImpl() {
    return new CustomNonbondedForceImpl(*this);
}
