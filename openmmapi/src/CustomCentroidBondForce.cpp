/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2015-2016 Stanford University and the Authors.      *
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
#include "openmm/CustomCentroidBondForce.h"
#include "openmm/internal/AssertionUtilities.h"
#include "openmm/internal/CustomCentroidBondForceImpl.h"
#include <cmath>
#include <map>
#include <set>
#include <sstream>
#include <utility>

using namespace OpenMM;
using std::map;
using std::pair;
using std::set;
using std::string;
using std::stringstream;
using std::vector;

CustomCentroidBondForce::CustomCentroidBondForce(int numGroups, const string& energy) : groupsPerBond(numGroups), energyExpression(energy), usePeriodic(false) {
}

CustomCentroidBondForce::~CustomCentroidBondForce() {
    for (auto function : functions)
        delete function.function;
}

const string& CustomCentroidBondForce::getEnergyFunction() const {
    return energyExpression;
}

void CustomCentroidBondForce::setEnergyFunction(const std::string& energy) {
    energyExpression = energy;
}

int CustomCentroidBondForce::addPerBondParameter(const string& name) {
    bondParameters.push_back(BondParameterInfo(name));
    return bondParameters.size()-1;
}

const string& CustomCentroidBondForce::getPerBondParameterName(int index) const {
    ASSERT_VALID_INDEX(index, bondParameters);
    return bondParameters[index].name;
}

void CustomCentroidBondForce::setPerBondParameterName(int index, const string& name) {
    ASSERT_VALID_INDEX(index, bondParameters);
    bondParameters[index].name = name;
}

int CustomCentroidBondForce::addGlobalParameter(const string& name, double defaultValue) {
    globalParameters.push_back(GlobalParameterInfo(name, defaultValue));
    return globalParameters.size()-1;
}

const string& CustomCentroidBondForce::getGlobalParameterName(int index) const {
    ASSERT_VALID_INDEX(index, globalParameters);
    return globalParameters[index].name;
}

void CustomCentroidBondForce::setGlobalParameterName(int index, const string& name) {
    ASSERT_VALID_INDEX(index, globalParameters);
    globalParameters[index].name = name;
}

double CustomCentroidBondForce::getGlobalParameterDefaultValue(int index) const {
    ASSERT_VALID_INDEX(index, globalParameters);
    return globalParameters[index].defaultValue;
}

void CustomCentroidBondForce::setGlobalParameterDefaultValue(int index, double defaultValue) {
    ASSERT_VALID_INDEX(index, globalParameters);
    globalParameters[index].defaultValue = defaultValue;
}

void CustomCentroidBondForce::addEnergyParameterDerivative(const string& name) {
    for (int i = 0; i < globalParameters.size(); i++)
        if (name == globalParameters[i].name) {
            energyParameterDerivatives.push_back(i);
            return;
        }
    throw OpenMMException(string("addEnergyParameterDerivative: Unknown global parameter '"+name+"'"));
}

const string& CustomCentroidBondForce::getEnergyParameterDerivativeName(int index) const {
    ASSERT_VALID_INDEX(index, energyParameterDerivatives);
    return globalParameters[energyParameterDerivatives[index]].name;
}

int CustomCentroidBondForce::addGroup(const vector<int>& particles, const vector<double>& weights) {
    if (particles.size() != weights.size() && weights.size() > 0)
        throw OpenMMException("CustomCentroidBondForce: wrong number of weights specified for a group.");
    groups.push_back(GroupInfo(particles, weights));
    return groups.size()-1;
}

void CustomCentroidBondForce::getGroupParameters(int index, vector<int>& particles, std::vector<double>& weights) const {
    ASSERT_VALID_INDEX(index, groups);
    particles = groups[index].particles;
    weights = groups[index].weights;
}

void CustomCentroidBondForce::setGroupParameters(int index, const vector<int>& particles, const vector<double>& weights) {
    ASSERT_VALID_INDEX(index, groups);
    if (particles.size() != weights.size() && weights.size() > 0)
        throw OpenMMException("CustomCentroidBondForce: wrong number of weights specified for a group.");
    groups[index].particles = particles;
    groups[index].weights = weights;
}

int CustomCentroidBondForce::addBond(const vector<int>& groups, const vector<double>& parameters) {
    if (groups.size() != groupsPerBond)
        throw OpenMMException("CustomCentroidBondForce: wrong number of groups specified for a bond.");
    bonds.push_back(BondInfo(groups, parameters));
    return bonds.size()-1;
}

void CustomCentroidBondForce::getBondParameters(int index, vector<int>& groups, std::vector<double>& parameters) const {
    ASSERT_VALID_INDEX(index, bonds);
    groups = bonds[index].groups;
    parameters = bonds[index].parameters;
}

void CustomCentroidBondForce::setBondParameters(int index, const vector<int>& groups, const vector<double>& parameters) {
    ASSERT_VALID_INDEX(index, bonds);
    if (groups.size() != groupsPerBond)
        throw OpenMMException("CustomCentroidBondForce: wrong number of groups specified for a bond.");
    bonds[index].groups = groups;
    bonds[index].parameters = parameters;
}

int CustomCentroidBondForce::addTabulatedFunction(const std::string& name, TabulatedFunction* function) {
    functions.push_back(FunctionInfo(name, function));
    return functions.size()-1;
}

const TabulatedFunction& CustomCentroidBondForce::getTabulatedFunction(int index) const {
    ASSERT_VALID_INDEX(index, functions);
    return *functions[index].function;
}

TabulatedFunction& CustomCentroidBondForce::getTabulatedFunction(int index) {
    ASSERT_VALID_INDEX(index, functions);
    return *functions[index].function;
}

const string& CustomCentroidBondForce::getTabulatedFunctionName(int index) const {
    ASSERT_VALID_INDEX(index, functions);
    return functions[index].name;
}

ForceImpl* CustomCentroidBondForce::createImpl() const {
    return new CustomCentroidBondForceImpl(*this);
}

void CustomCentroidBondForce::updateParametersInContext(Context& context) {
    dynamic_cast<CustomCentroidBondForceImpl&>(getImplInContext(context)).updateParametersInContext(getContextImpl(context));
}

void CustomCentroidBondForce::setUsesPeriodicBoundaryConditions(bool periodic) {
    usePeriodic = periodic;
}

bool CustomCentroidBondForce::usesPeriodicBoundaryConditions() const {
    return usePeriodic;
}
