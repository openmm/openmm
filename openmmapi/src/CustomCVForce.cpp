/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2008-2018 Stanford University and the Authors.      *
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

#include "openmm/OpenMMException.h"
#include "openmm/CustomCVForce.h"
#include "openmm/internal/AssertionUtilities.h"
#include "openmm/internal/CustomCVForceImpl.h"
#include <cmath>
#include <map>
#include <set>
#include <utility>

using namespace OpenMM;
using namespace std;

CustomCVForce::CustomCVForce(const string& energy) : energyExpression(energy) {
}

CustomCVForce::~CustomCVForce() {
    for (auto variable : variables)
        delete variable.variable;
    for (auto function : functions)
        delete function.function;
}

const string& CustomCVForce::getEnergyFunction() const {
    return energyExpression;
}

void CustomCVForce::setEnergyFunction(const std::string& energy) {
    energyExpression = energy;
}

int CustomCVForce::addCollectiveVariable(const std::string& name, Force* variable) {
    if (variables.size() >= 32)
        throw OpenMMException("CustomCVForce cannot have more than 32 collective variables");
    variables.push_back(VariableInfo(name, variable));
    return variables.size()-1;
}

const string& CustomCVForce::getCollectiveVariableName(int index) const {
    ASSERT_VALID_INDEX(index, variables);
    return variables[index].name;
}

Force& CustomCVForce::getCollectiveVariable(int index) {
    ASSERT_VALID_INDEX(index, variables);
    return *variables[index].variable;
}

const Force& CustomCVForce::getCollectiveVariable(int index) const {
    ASSERT_VALID_INDEX(index, variables);
    return *variables[index].variable;
}

int CustomCVForce::addGlobalParameter(const string& name, double defaultValue) {
    globalParameters.push_back(GlobalParameterInfo(name, defaultValue));
    return globalParameters.size()-1;
}

const string& CustomCVForce::getGlobalParameterName(int index) const {
    ASSERT_VALID_INDEX(index, globalParameters);
    return globalParameters[index].name;
}

void CustomCVForce::setGlobalParameterName(int index, const string& name) {
    ASSERT_VALID_INDEX(index, globalParameters);
    globalParameters[index].name = name;
}

double CustomCVForce::getGlobalParameterDefaultValue(int index) const {
    ASSERT_VALID_INDEX(index, globalParameters);
    return globalParameters[index].defaultValue;
}

void CustomCVForce::setGlobalParameterDefaultValue(int index, double defaultValue) {
    ASSERT_VALID_INDEX(index, globalParameters);
    globalParameters[index].defaultValue = defaultValue;
}

void CustomCVForce::addEnergyParameterDerivative(const string& name) {
    for (int i = 0; i < globalParameters.size(); i++)
        if (name == globalParameters[i].name) {
            energyParameterDerivatives.push_back(i);
            return;
        }
    throw OpenMMException(string("addEnergyParameterDerivative: Unknown global parameter '"+name+"'"));
}

const string& CustomCVForce::getEnergyParameterDerivativeName(int index) const {
    ASSERT_VALID_INDEX(index, energyParameterDerivatives);
    return globalParameters[energyParameterDerivatives[index]].name;
}

int CustomCVForce::addTabulatedFunction(const std::string& name, TabulatedFunction* function) {
    functions.push_back(FunctionInfo(name, function));
    return functions.size()-1;
}

const TabulatedFunction& CustomCVForce::getTabulatedFunction(int index) const {
    ASSERT_VALID_INDEX(index, functions);
    return *functions[index].function;
}

TabulatedFunction& CustomCVForce::getTabulatedFunction(int index) {
    ASSERT_VALID_INDEX(index, functions);
    return *functions[index].function;
}

const string& CustomCVForce::getTabulatedFunctionName(int index) const {
    ASSERT_VALID_INDEX(index, functions);
    return functions[index].name;
}

void CustomCVForce::getCollectiveVariableValues(Context& context, vector<double>& values) {
    dynamic_cast<CustomCVForceImpl&>(getImplInContext(context)).getCollectiveVariableValues(getContextImpl(context), values);
}

ForceImpl* CustomCVForce::createImpl() const {
    return new CustomCVForceImpl(*this);
}

Context& CustomCVForce::getInnerContext(Context& context) {
    return dynamic_cast<CustomCVForceImpl&>(getImplInContext(context)).getInnerContext();
}

void CustomCVForce::updateParametersInContext(Context& context) {
    dynamic_cast<CustomCVForceImpl&>(getImplInContext(context)).updateParametersInContext(getContextImpl(context));
}

bool CustomCVForce::usesPeriodicBoundaryConditions() const {
    for (auto& variable : variables)
        if (variable.variable->usesPeriodicBoundaryConditions())
            return true;
    return false;
}
