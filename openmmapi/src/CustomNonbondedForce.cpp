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

int CustomNonbondedForce::addParameter(const string& name, const string& combiningRule) {
    parameters.push_back(ParameterInfo(name, combiningRule));
    return parameters.size()-1;
}

const string& CustomNonbondedForce::getParameterName(int index) const {
    return parameters[index].name;
}

void CustomNonbondedForce::setParameterName(int index, const string& name) {
    parameters[index].name = name;
}

const string& CustomNonbondedForce::getParameterCombiningRule(int index) const {
    return parameters[index].combiningRule;
}

void CustomNonbondedForce::setParameterCombiningRule(int index, const string& combiningRule) {
    parameters[index].combiningRule = combiningRule;
}

int CustomNonbondedForce::addGlobalParameter(const string& name, double defaultValue) {
    globalParameters.push_back(GlobalParameterInfo(name, defaultValue));
    return globalParameters.size()-1;
}

const string& CustomNonbondedForce::getGlobalParameterName(int index) const {
    return globalParameters[index].name;
}

void CustomNonbondedForce::setGlobalParameterName(int index, const string& name) {
    globalParameters[index].name = name;
}

double CustomNonbondedForce::getGlobalParameterDefaultValue(int index) const {
    return globalParameters[index].defaultValue;
}

void CustomNonbondedForce::setGlobalParameterDefaultValue(int index, double defaultValue) {
    globalParameters[index].defaultValue = defaultValue;
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

int CustomNonbondedForce::addFunction(const std::string& name, const std::vector<double>& values, double min, double max, bool interpolating) {
    if (max <= min)
        throw OpenMMException("CustomNonbondedForce: max <= min for a tabulated function.");
    if (values.size() < 2)
        throw OpenMMException("CustomNonbondedForce: a tabulated function must have at least two points");
    functions.push_back(FunctionInfo(name, values, min, max, interpolating));
    return functions.size()-1;
}

void CustomNonbondedForce::getFunctionParameters(int index, std::string& name, std::vector<double>& values, double& min, double& max, bool& interpolating) const {
    name = functions[index].name;
    values = functions[index].values;
    min = functions[index].min;
    max = functions[index].max;
    interpolating = functions[index].interpolating;
}

void CustomNonbondedForce::setFunctionParameters(int index, const std::string& name, const std::vector<double>& values, double min, double max, bool interpolating) {
    if (max <= min)
        throw OpenMMException("CustomNonbondedForce: max <= min for a tabulated function.");
    if (values.size() < 2)
        throw OpenMMException("CustomNonbondedForce: a tabulated function must have at least two points");
    functions[index].name = name;
    functions[index].values = values;
    functions[index].min = min;
    functions[index].max = max;
    functions[index].interpolating = interpolating;
}

ForceImpl* CustomNonbondedForce::createImpl() {
    return new CustomNonbondedForceImpl(*this);
}
