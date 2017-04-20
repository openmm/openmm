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
#include "openmm/CustomManyParticleForce.h"
#include "openmm/internal/AssertionUtilities.h"
#include "openmm/internal/CustomManyParticleForceImpl.h"
#include <cmath>
#include <map>
#include <set>
#include <sstream>
#include <utility>

using namespace OpenMM;
using namespace std;

CustomManyParticleForce::CustomManyParticleForce(int particlesPerSet, const string& energy) :
        particlesPerSet(particlesPerSet), energyExpression(energy), nonbondedMethod(NoCutoff), permutationMode(SinglePermutation), cutoffDistance(1.0), typeFilters(particlesPerSet) {
}

CustomManyParticleForce::~CustomManyParticleForce() {
    for (auto function : functions)
        delete function.function;
}

const string& CustomManyParticleForce::getEnergyFunction() const {
    return energyExpression;
}

void CustomManyParticleForce::setEnergyFunction(const string& energy) {
    energyExpression = energy;
}

CustomManyParticleForce::NonbondedMethod CustomManyParticleForce::getNonbondedMethod() const {
    return nonbondedMethod;
}

void CustomManyParticleForce::setNonbondedMethod(NonbondedMethod method) {
    if (method < 0 || method > 2)
        throw OpenMMException("CustomManyParticleForce: Illegal value for nonbonded method");
    nonbondedMethod = method;
}

CustomManyParticleForce::PermutationMode CustomManyParticleForce::getPermutationMode() const {
    return permutationMode;
}

void CustomManyParticleForce::setPermutationMode(PermutationMode mode) {
    permutationMode = mode;
}

double CustomManyParticleForce::getCutoffDistance() const {
    return cutoffDistance;
}

void CustomManyParticleForce::setCutoffDistance(double distance) {
    cutoffDistance = distance;
}

int CustomManyParticleForce::addPerParticleParameter(const string& name) {
    particleParameters.push_back(ParticleParameterInfo(name));
    return particleParameters.size()-1;
}

const string& CustomManyParticleForce::getPerParticleParameterName(int index) const {
    ASSERT_VALID_INDEX(index, particleParameters);
    return particleParameters[index].name;
}

void CustomManyParticleForce::setPerParticleParameterName(int index, const string& name) {
    ASSERT_VALID_INDEX(index, particleParameters);
    particleParameters[index].name = name;
}

int CustomManyParticleForce::addGlobalParameter(const string& name, double defaultValue) {
    globalParameters.push_back(GlobalParameterInfo(name, defaultValue));
    return globalParameters.size()-1;
}

const string& CustomManyParticleForce::getGlobalParameterName(int index) const {
    ASSERT_VALID_INDEX(index, globalParameters);
    return globalParameters[index].name;
}

void CustomManyParticleForce::setGlobalParameterName(int index, const string& name) {
    ASSERT_VALID_INDEX(index, globalParameters);
    globalParameters[index].name = name;
}

double CustomManyParticleForce::getGlobalParameterDefaultValue(int index) const {
    ASSERT_VALID_INDEX(index, globalParameters);
    return globalParameters[index].defaultValue;
}

void CustomManyParticleForce::setGlobalParameterDefaultValue(int index, double defaultValue) {
    ASSERT_VALID_INDEX(index, globalParameters);
    globalParameters[index].defaultValue = defaultValue;
}

int CustomManyParticleForce::addParticle(const vector<double>& parameters, int type) {
    particles.push_back(ParticleInfo(parameters, type));
    return particles.size()-1;
}

void CustomManyParticleForce::getParticleParameters(int index, vector<double>& parameters, int& type) const {
    ASSERT_VALID_INDEX(index, particles);
    parameters = particles[index].parameters;
    type = particles[index].type;
}

void CustomManyParticleForce::setParticleParameters(int index, const vector<double>& parameters, int type) {
    ASSERT_VALID_INDEX(index, particles);
    particles[index].parameters = parameters;
    particles[index].type = type;
}

int CustomManyParticleForce::addExclusion(int particle1, int particle2) {
    exclusions.push_back(ExclusionInfo(particle1, particle2));
    return exclusions.size()-1;
}
void CustomManyParticleForce::getExclusionParticles(int index, int& particle1, int& particle2) const {
    ASSERT_VALID_INDEX(index, exclusions);
    particle1 = exclusions[index].particle1;
    particle2 = exclusions[index].particle2;
}

void CustomManyParticleForce::setExclusionParticles(int index, int particle1, int particle2) {
    ASSERT_VALID_INDEX(index, exclusions);
    exclusions[index].particle1 = particle1;
    exclusions[index].particle2 = particle2;
}

void CustomManyParticleForce::createExclusionsFromBonds(const vector<pair<int, int> >& bonds, int bondCutoff) {
    if (bondCutoff < 1)
        return;
    vector<set<int> > exclusions(particles.size());
    vector<set<int> > bonded12(exclusions.size());
    for (auto& bond : bonds) {
        int p1 = bond.first;
        int p2 = bond.second;
        exclusions[p1].insert(p2);
        exclusions[p2].insert(p1);
        bonded12[p1].insert(p2);
        bonded12[p2].insert(p1);
    }
    for (int level = 0; level < bondCutoff-1; level++) {
        vector<set<int> > currentExclusions = exclusions;
        for (int i = 0; i < (int) particles.size(); i++)
            for (int j : currentExclusions[i])
                exclusions[j].insert(bonded12[i].begin(), bonded12[i].end());
    }
    for (int i = 0; i < (int) exclusions.size(); ++i)
        for (int j : exclusions[i])
            if (j < i)
                addExclusion(j, i);
}

void CustomManyParticleForce::getTypeFilter(int index, set<int>& types) const {
    if (index < 0 || index >= particlesPerSet)
        throw OpenMMException("CustomManyParticleForce: index to getTypeFilter out of range");
    types = typeFilters[index];
}

void CustomManyParticleForce::setTypeFilter(int index, const set<int>& types) {
    if (index < 0 || index >= particlesPerSet)
        throw OpenMMException("CustomManyParticleForce: index to setTypeFilter out of range");
    typeFilters[index] = types;
}

int CustomManyParticleForce::addTabulatedFunction(const string& name, TabulatedFunction* function) {
    functions.push_back(FunctionInfo(name, function));
    return functions.size()-1;
}

const TabulatedFunction& CustomManyParticleForce::getTabulatedFunction(int index) const {
    ASSERT_VALID_INDEX(index, functions);
    return *functions[index].function;
}

TabulatedFunction& CustomManyParticleForce::getTabulatedFunction(int index) {
    ASSERT_VALID_INDEX(index, functions);
    return *functions[index].function;
}

const string& CustomManyParticleForce::getTabulatedFunctionName(int index) const {
    ASSERT_VALID_INDEX(index, functions);
    return functions[index].name;
}

ForceImpl* CustomManyParticleForce::createImpl() const {
    return new CustomManyParticleForceImpl(*this);
}

void CustomManyParticleForce::updateParametersInContext(Context& context) {
    dynamic_cast<CustomManyParticleForceImpl&>(getImplInContext(context)).updateParametersInContext(getContextImpl(context));
}
