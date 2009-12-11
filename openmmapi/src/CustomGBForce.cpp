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
#include "openmm/CustomGBForce.h"
#include "openmm/internal/CustomGBForceImpl.h"
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

CustomGBForce::CustomGBForce() : nonbondedMethod(NoCutoff), cutoffDistance(1.0) {
}

CustomGBForce::NonbondedMethod CustomGBForce::getNonbondedMethod() const {
    return nonbondedMethod;
}

void CustomGBForce::setNonbondedMethod(NonbondedMethod method) {
    nonbondedMethod = method;
}

double CustomGBForce::getCutoffDistance() const {
    return cutoffDistance;
}

void CustomGBForce::setCutoffDistance(double distance) {
    cutoffDistance = distance;
}

int CustomGBForce::addPerParticleParameter(const string& name) {
    parameters.push_back(PerParticleParameterInfo(name));
    return parameters.size()-1;
}

const string& CustomGBForce::getPerParticleParameterName(int index) const {
    return parameters[index].name;
}

void CustomGBForce::setPerParticleParameterName(int index, const string& name) {
    parameters[index].name = name;
}

int CustomGBForce::addGlobalParameter(const string& name, double defaultValue) {
    globalParameters.push_back(GlobalParameterInfo(name, defaultValue));
    return globalParameters.size()-1;
}

const string& CustomGBForce::getGlobalParameterName(int index) const {
    return globalParameters[index].name;
}

void CustomGBForce::setGlobalParameterName(int index, const string& name) {
    globalParameters[index].name = name;
}

double CustomGBForce::getGlobalParameterDefaultValue(int index) const {
    return globalParameters[index].defaultValue;
}

void CustomGBForce::setGlobalParameterDefaultValue(int index, double defaultValue) {
    globalParameters[index].defaultValue = defaultValue;
}

int CustomGBForce::addParticle(const vector<double>& parameters) {
    particles.push_back(ParticleInfo(parameters));
    return particles.size()-1;
}

void CustomGBForce::getParticleParameters(int index, std::vector<double>& parameters) const {
    parameters = particles[index].parameters;
}

void CustomGBForce::setParticleParameters(int index, const vector<double>& parameters) {
    particles[index].parameters = parameters;
}

int CustomGBForce::addComputedValue(const std::string& name, const std::string& expression, ComputationType type) {
    computedValues.push_back(CustomGBForce::ComputationInfo(name, expression, type));
    return computedValues.size()-1;
}

void CustomGBForce::getComputedValueParameters(int index, std::string& name, std::string& expression, ComputationType& type) const {
    name = computedValues[index].name;
    expression = computedValues[index].expression;
    type = computedValues[index].type;
}

void CustomGBForce::setComputedValueParameters(int index, const std::string& name, const std::string& expression, ComputationType type) {
    computedValues[index].name = name;
    computedValues[index].expression = expression;
    computedValues[index].type = type;
}

int CustomGBForce::addEnergyTerm(const std::string& expression, ComputationType type) {
    energyTerms.push_back(CustomGBForce::ComputationInfo("", expression, type));
    return energyTerms.size()-1;
}

void CustomGBForce::getEnergyTermParameters(int index, std::string& expression, ComputationType& type) const {
    expression = energyTerms[index].expression;
    type = energyTerms[index].type;
}

void CustomGBForce::setEnergyTermParameters(int index, const std::string& expression, ComputationType type) {
    energyTerms[index].expression = expression;
    energyTerms[index].type = type;
}

int CustomGBForce::addExclusion(int particle1, int particle2) {
    exclusions.push_back(ExclusionInfo(particle1, particle2));
    return exclusions.size()-1;
}
void CustomGBForce::getExclusionParticles(int index, int& particle1, int& particle2) const {
    particle1 = exclusions[index].particle1;
    particle2 = exclusions[index].particle2;
}

void CustomGBForce::setExclusionParticles(int index, int particle1, int particle2) {
    exclusions[index].particle1 = particle1;
    exclusions[index].particle2 = particle2;
}

int CustomGBForce::addFunction(const std::string& name, const std::vector<double>& values, double min, double max, bool interpolating) {
    if (max <= min)
        throw OpenMMException("CustomGBForce: max <= min for a tabulated function.");
    if (values.size() < 2)
        throw OpenMMException("CustomGBForce: a tabulated function must have at least two points");
    functions.push_back(FunctionInfo(name, values, min, max, interpolating));
    return functions.size()-1;
}

void CustomGBForce::getFunctionParameters(int index, std::string& name, std::vector<double>& values, double& min, double& max, bool& interpolating) const {
    name = functions[index].name;
    values = functions[index].values;
    min = functions[index].min;
    max = functions[index].max;
    interpolating = functions[index].interpolating;
}

void CustomGBForce::setFunctionParameters(int index, const std::string& name, const std::vector<double>& values, double min, double max, bool interpolating) {
    if (max <= min)
        throw OpenMMException("CustomGBForce: max <= min for a tabulated function.");
    if (values.size() < 2)
        throw OpenMMException("CustomGBForce: a tabulated function must have at least two points");
    functions[index].name = name;
    functions[index].values = values;
    functions[index].min = min;
    functions[index].max = max;
    functions[index].interpolating = interpolating;
}

ForceImpl* CustomGBForce::createImpl() {
    return new CustomGBForceImpl(*this);
}

