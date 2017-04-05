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
#include "openmm/CustomCompoundBondForce.h"
#include "openmm/internal/AssertionUtilities.h"
#include "openmm/internal/CustomCompoundBondForceImpl.h"
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

CustomCompoundBondForce::CustomCompoundBondForce(int numParticles, const string& energy) : particlesPerBond(numParticles), energyExpression(energy), usePeriodic(false) {
}


CustomCompoundBondForce::~CustomCompoundBondForce() {
    for (auto function : functions)
        delete function.function;
}

const string& CustomCompoundBondForce::getEnergyFunction() const {
    return energyExpression;
}

void CustomCompoundBondForce::setEnergyFunction(const std::string& energy) {
    energyExpression = energy;
}

int CustomCompoundBondForce::addPerBondParameter(const string& name) {
    bondParameters.push_back(BondParameterInfo(name));
    return bondParameters.size()-1;
}

const string& CustomCompoundBondForce::getPerBondParameterName(int index) const {
    ASSERT_VALID_INDEX(index, bondParameters);
    return bondParameters[index].name;
}

void CustomCompoundBondForce::setPerBondParameterName(int index, const string& name) {
    ASSERT_VALID_INDEX(index, bondParameters);
    bondParameters[index].name = name;
}

int CustomCompoundBondForce::addGlobalParameter(const string& name, double defaultValue) {
    globalParameters.push_back(GlobalParameterInfo(name, defaultValue));
    return globalParameters.size()-1;
}

const string& CustomCompoundBondForce::getGlobalParameterName(int index) const {
    ASSERT_VALID_INDEX(index, globalParameters);
    return globalParameters[index].name;
}

void CustomCompoundBondForce::setGlobalParameterName(int index, const string& name) {
    ASSERT_VALID_INDEX(index, globalParameters);
    globalParameters[index].name = name;
}

double CustomCompoundBondForce::getGlobalParameterDefaultValue(int index) const {
    ASSERT_VALID_INDEX(index, globalParameters);
    return globalParameters[index].defaultValue;
}

void CustomCompoundBondForce::setGlobalParameterDefaultValue(int index, double defaultValue) {
    ASSERT_VALID_INDEX(index, globalParameters);
    globalParameters[index].defaultValue = defaultValue;
}

void CustomCompoundBondForce::addEnergyParameterDerivative(const string& name) {
    for (int i = 0; i < globalParameters.size(); i++)
        if (name == globalParameters[i].name) {
            energyParameterDerivatives.push_back(i);
            return;
        }
    throw OpenMMException(string("addEnergyParameterDerivative: Unknown global parameter '"+name+"'"));
}

const string& CustomCompoundBondForce::getEnergyParameterDerivativeName(int index) const {
    ASSERT_VALID_INDEX(index, energyParameterDerivatives);
    return globalParameters[energyParameterDerivatives[index]].name;
}

int CustomCompoundBondForce::addBond(const vector<int>& particles, const vector<double>& parameters) {
    if (particles.size() != particlesPerBond)
        throw OpenMMException("CustomCompoundBondForce: wrong number of particles specified for a bond.");
    bonds.push_back(BondInfo(particles, parameters));
    return bonds.size()-1;
}

void CustomCompoundBondForce::getBondParameters(int index, vector<int>& particles, std::vector<double>& parameters) const {
    ASSERT_VALID_INDEX(index, bonds);
    particles = bonds[index].particles;
    parameters = bonds[index].parameters;
}

void CustomCompoundBondForce::setBondParameters(int index, const vector<int>& particles, const vector<double>& parameters) {
    ASSERT_VALID_INDEX(index, bonds);
    if (particles.size() != particlesPerBond)
        throw OpenMMException("CustomCompoundBondForce: wrong number of particles specified for a bond.");
    bonds[index].particles = particles;
    bonds[index].parameters = parameters;
}

int CustomCompoundBondForce::addTabulatedFunction(const std::string& name, TabulatedFunction* function) {
    functions.push_back(FunctionInfo(name, function));
    return functions.size()-1;
}

const TabulatedFunction& CustomCompoundBondForce::getTabulatedFunction(int index) const {
    ASSERT_VALID_INDEX(index, functions);
    return *functions[index].function;
}

TabulatedFunction& CustomCompoundBondForce::getTabulatedFunction(int index) {
    ASSERT_VALID_INDEX(index, functions);
    return *functions[index].function;
}

const string& CustomCompoundBondForce::getTabulatedFunctionName(int index) const {
    ASSERT_VALID_INDEX(index, functions);
    return functions[index].name;
}

int CustomCompoundBondForce::addFunction(const std::string& name, const std::vector<double>& values, double min, double max) {
    functions.push_back(FunctionInfo(name, new Continuous1DFunction(values, min, max)));
    return functions.size()-1;
}

void CustomCompoundBondForce::getFunctionParameters(int index, std::string& name, std::vector<double>& values, double& min, double& max) const {
    ASSERT_VALID_INDEX(index, functions);
    Continuous1DFunction* function = dynamic_cast<Continuous1DFunction*>(functions[index].function);
    if (function == NULL)
        throw OpenMMException("CustomCompoundBondForce: function is not a Continuous1DFunction");
    name = functions[index].name;
    function->getFunctionParameters(values, min, max);
}

void CustomCompoundBondForce::setFunctionParameters(int index, const std::string& name, const std::vector<double>& values, double min, double max) {
    ASSERT_VALID_INDEX(index, functions);
    Continuous1DFunction* function = dynamic_cast<Continuous1DFunction*>(functions[index].function);
    if (function == NULL)
        throw OpenMMException("CustomCompoundBondForce: function is not a Continuous1DFunction");
    functions[index].name = name;
    function->setFunctionParameters(values, min, max);
}

ForceImpl* CustomCompoundBondForce::createImpl() const {
    return new CustomCompoundBondForceImpl(*this);
}

void CustomCompoundBondForce::updateParametersInContext(Context& context) {
    dynamic_cast<CustomCompoundBondForceImpl&>(getImplInContext(context)).updateParametersInContext(getContextImpl(context));
}

void CustomCompoundBondForce::setUsesPeriodicBoundaryConditions(bool periodic) {
    usePeriodic = periodic;
}

bool CustomCompoundBondForce::usesPeriodicBoundaryConditions() const {
    return usePeriodic;
}
