/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2008-2012 Stanford University and the Authors.      *
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
#include "openmm/internal/AssertionUtilities.h"
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

CustomNonbondedForce::CustomNonbondedForce(const string& energy) : energyExpression(energy), nonbondedMethod(NoCutoff), cutoffDistance(1.0),
    switchingDistance(-1.0), useSwitchingFunction(false), useLongRangeCorrection(false) {
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

bool CustomNonbondedForce::getUseSwitchingFunction() const {
    return useSwitchingFunction;
}

void CustomNonbondedForce::setUseSwitchingFunction(bool use) {
    useSwitchingFunction = use;
}

double CustomNonbondedForce::getSwitchingDistance() const {
    return switchingDistance;
}

void CustomNonbondedForce::setSwitchingDistance(double distance) {
    switchingDistance = distance;
}

bool CustomNonbondedForce::getUseLongRangeCorrection() const {
    return useLongRangeCorrection;
}

void CustomNonbondedForce::setUseLongRangeCorrection(bool use) {
    useLongRangeCorrection = use;
}

int CustomNonbondedForce::addPerParticleParameter(const string& name) {
    parameters.push_back(PerParticleParameterInfo(name));
    return parameters.size()-1;
}

const string& CustomNonbondedForce::getPerParticleParameterName(int index) const {
    ASSERT_VALID_INDEX(index, parameters);
    return parameters[index].name;
}

void CustomNonbondedForce::setPerParticleParameterName(int index, const string& name) {
    ASSERT_VALID_INDEX(index, parameters);
    parameters[index].name = name;
}

int CustomNonbondedForce::addGlobalParameter(const string& name, double defaultValue) {
    globalParameters.push_back(GlobalParameterInfo(name, defaultValue));
    return globalParameters.size()-1;
}

const string& CustomNonbondedForce::getGlobalParameterName(int index) const {
    ASSERT_VALID_INDEX(index, globalParameters);
    return globalParameters[index].name;
}

void CustomNonbondedForce::setGlobalParameterName(int index, const string& name) {
    ASSERT_VALID_INDEX(index, globalParameters);
    globalParameters[index].name = name;
}

double CustomNonbondedForce::getGlobalParameterDefaultValue(int index) const {
    ASSERT_VALID_INDEX(index, globalParameters);
    return globalParameters[index].defaultValue;
}

void CustomNonbondedForce::setGlobalParameterDefaultValue(int index, double defaultValue) {
    ASSERT_VALID_INDEX(index, globalParameters);
    globalParameters[index].defaultValue = defaultValue;
}

int CustomNonbondedForce::addParticle(const vector<double>& parameters) {
    particles.push_back(ParticleInfo(parameters));
    return particles.size()-1;
}

void CustomNonbondedForce::getParticleParameters(int index, std::vector<double>& parameters) const {
    ASSERT_VALID_INDEX(index, particles);
    parameters = particles[index].parameters;
}

void CustomNonbondedForce::setParticleParameters(int index, const vector<double>& parameters) {
    ASSERT_VALID_INDEX(index, particles);
    particles[index].parameters = parameters;
}

int CustomNonbondedForce::addExclusion(int particle1, int particle2) {
    exclusions.push_back(ExclusionInfo(particle1, particle2));
    return exclusions.size()-1;
}
void CustomNonbondedForce::getExclusionParticles(int index, int& particle1, int& particle2) const {
    ASSERT_VALID_INDEX(index, exclusions);
    particle1 = exclusions[index].particle1;
    particle2 = exclusions[index].particle2;
}

void CustomNonbondedForce::setExclusionParticles(int index, int particle1, int particle2) {
    ASSERT_VALID_INDEX(index, exclusions);
    exclusions[index].particle1 = particle1;
    exclusions[index].particle2 = particle2;
}

int CustomNonbondedForce::addFunction(const std::string& name, const std::vector<double>& values, double min, double max) {
    if (max <= min)
        throw OpenMMException("CustomNonbondedForce: max <= min for a tabulated function.");
    if (values.size() < 2)
        throw OpenMMException("CustomNonbondedForce: a tabulated function must have at least two points");
    functions.push_back(FunctionInfo(name, values, min, max));
    return functions.size()-1;
}

void CustomNonbondedForce::getFunctionParameters(int index, std::string& name, std::vector<double>& values, double& min, double& max) const {
    ASSERT_VALID_INDEX(index, functions);
    name = functions[index].name;
    values = functions[index].values;
    min = functions[index].min;
    max = functions[index].max;
}

void CustomNonbondedForce::setFunctionParameters(int index, const std::string& name, const std::vector<double>& values, double min, double max) {
    if (max <= min)
        throw OpenMMException("CustomNonbondedForce: max <= min for a tabulated function.");
    if (values.size() < 2)
        throw OpenMMException("CustomNonbondedForce: a tabulated function must have at least two points");
    ASSERT_VALID_INDEX(index, functions);
    functions[index].name = name;
    functions[index].values = values;
    functions[index].min = min;
    functions[index].max = max;
}

int CustomNonbondedForce::addInteractionGroup(const std::set<int>& set1, const std::set<int>& set2) {
    interactionGroups.push_back(InteractionGroupInfo(set1, set2));
    return interactionGroups.size()-1;
}

void CustomNonbondedForce::getInteractionGroupParameters(int index, std::set<int>& set1, std::set<int>& set2) const {
    ASSERT_VALID_INDEX(index, interactionGroups);
    set1 = interactionGroups[index].set1;
    set2 = interactionGroups[index].set2;
}

void CustomNonbondedForce::setInteractionGroupParameters(int index, const std::set<int>& set1, const std::set<int>& set2) {
    ASSERT_VALID_INDEX(index, interactionGroups);
    interactionGroups[index].set1 = set1;
    interactionGroups[index].set2 = set2;
}

ForceImpl* CustomNonbondedForce::createImpl() const {
    return new CustomNonbondedForceImpl(*this);
}

void CustomNonbondedForce::updateParametersInContext(Context& context) {
    dynamic_cast<CustomNonbondedForceImpl&>(getImplInContext(context)).updateParametersInContext(getContextImpl(context));
}
