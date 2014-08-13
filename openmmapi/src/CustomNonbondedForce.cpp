/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2008-2014 Stanford University and the Authors.      *
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

CustomNonbondedForce::CustomNonbondedForce(const CustomNonbondedForce& rhs) {
    // Copy everything and deep copy the tabulated functions
    energyExpression = rhs.energyExpression;
    nonbondedMethod = rhs.nonbondedMethod;
    cutoffDistance = rhs.cutoffDistance;
    switchingDistance = rhs.switchingDistance;
    useSwitchingFunction = rhs.useSwitchingFunction;
    useLongRangeCorrection = rhs.useLongRangeCorrection;
    parameters = rhs.parameters;
    globalParameters = rhs.globalParameters;
    particles = rhs.particles;
    exclusions = rhs.exclusions;
    interactionGroups = rhs.interactionGroups;
    for (vector<FunctionInfo>::const_iterator it = rhs.functions.begin(); it != rhs.functions.end(); it++)
        functions.push_back(FunctionInfo(it->name, it->function->Copy()));
}

CustomNonbondedForce::~CustomNonbondedForce() {
    for (int i = 0; i < (int) functions.size(); i++)
        delete functions[i].function;
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

void CustomNonbondedForce::createExclusionsFromBonds(const vector<pair<int, int> >& bonds, int bondCutoff) {
    if (bondCutoff < 1)
        return;
    vector<set<int> > exclusions(particles.size());
    vector<set<int> > bonded12(exclusions.size());
    for (int i = 0; i < (int) bonds.size(); ++i) {
        int p1 = bonds[i].first;
        int p2 = bonds[i].second;
        exclusions[p1].insert(p2);
        exclusions[p2].insert(p1);
        bonded12[p1].insert(p2);
        bonded12[p2].insert(p1);
    }
    for (int level = 0; level < bondCutoff-1; level++) {
        vector<set<int> > currentExclusions = exclusions;
        for (int i = 0; i < (int) particles.size(); i++) {
            for (set<int>::const_iterator iter = currentExclusions[i].begin(); iter != currentExclusions[i].end(); ++iter)
                exclusions[*iter].insert(bonded12[i].begin(), bonded12[i].end());
        }
    }
    for (int i = 0; i < (int) exclusions.size(); ++i)
        for (set<int>::const_iterator iter = exclusions[i].begin(); iter != exclusions[i].end(); ++iter)
            if (*iter < i)
                addExclusion(*iter, i);
}

int CustomNonbondedForce::addTabulatedFunction(const std::string& name, TabulatedFunction* function) {
    functions.push_back(FunctionInfo(name, function));
    return functions.size()-1;
}

const TabulatedFunction& CustomNonbondedForce::getTabulatedFunction(int index) const {
    ASSERT_VALID_INDEX(index, functions);
    return *functions[index].function;
}

TabulatedFunction& CustomNonbondedForce::getTabulatedFunction(int index) {
    ASSERT_VALID_INDEX(index, functions);
    return *functions[index].function;
}

const string& CustomNonbondedForce::getTabulatedFunctionName(int index) const {
    ASSERT_VALID_INDEX(index, functions);
    return functions[index].name;
}

int CustomNonbondedForce::addFunction(const std::string& name, const std::vector<double>& values, double min, double max) {
    functions.push_back(FunctionInfo(name, new Continuous1DFunction(values, min, max)));
    return functions.size()-1;
}

void CustomNonbondedForce::getFunctionParameters(int index, std::string& name, std::vector<double>& values, double& min, double& max) const {
    ASSERT_VALID_INDEX(index, functions);
    Continuous1DFunction* function = dynamic_cast<Continuous1DFunction*>(functions[index].function);
    if (function == NULL)
        throw OpenMMException("CustomNonbondedForce: function is not a Continuous1DFunction");
    name = functions[index].name;
    function->getFunctionParameters(values, min, max);
}

void CustomNonbondedForce::setFunctionParameters(int index, const std::string& name, const std::vector<double>& values, double min, double max) {
    ASSERT_VALID_INDEX(index, functions);
    Continuous1DFunction* function = dynamic_cast<Continuous1DFunction*>(functions[index].function);
    if (function == NULL)
        throw OpenMMException("CustomNonbondedForce: function is not a Continuous1DFunction");
    functions[index].name = name;
    function->setFunctionParameters(values, min, max);
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
