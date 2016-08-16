/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2010-2016 Stanford University and the Authors.      *
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
#include "openmm/CustomTorsionForce.h"
#include "openmm/internal/AssertionUtilities.h"
#include "openmm/internal/CustomTorsionForceImpl.h"
#include <cmath>
#include <map>
#include <sstream>
#include <utility>

using namespace OpenMM;
using std::string;
using std::stringstream;
using std::vector;

CustomTorsionForce::CustomTorsionForce(const string& energy) : energyExpression(energy), usePeriodic(false) {
}

const string& CustomTorsionForce::getEnergyFunction() const {
    return energyExpression;
}

void CustomTorsionForce::setEnergyFunction(const std::string& energy) {
    energyExpression = energy;
}

int CustomTorsionForce::addPerTorsionParameter(const string& name) {
    parameters.push_back(TorsionParameterInfo(name));
    return parameters.size()-1;
}

const string& CustomTorsionForce::getPerTorsionParameterName(int index) const {
    ASSERT_VALID_INDEX(index, parameters);
    return parameters[index].name;
}

void CustomTorsionForce::setPerTorsionParameterName(int index, const string& name) {
    ASSERT_VALID_INDEX(index, parameters);
    parameters[index].name = name;
}

int CustomTorsionForce::addGlobalParameter(const string& name, double defaultValue) {
    globalParameters.push_back(GlobalParameterInfo(name, defaultValue));
    return globalParameters.size()-1;
}

const string& CustomTorsionForce::getGlobalParameterName(int index) const {
    ASSERT_VALID_INDEX(index, globalParameters);
    return globalParameters[index].name;
}

void CustomTorsionForce::setGlobalParameterName(int index, const string& name) {
    ASSERT_VALID_INDEX(index, globalParameters);
    globalParameters[index].name = name;
}

double CustomTorsionForce::getGlobalParameterDefaultValue(int index) const {
    ASSERT_VALID_INDEX(index, globalParameters);
    return globalParameters[index].defaultValue;
}

void CustomTorsionForce::setGlobalParameterDefaultValue(int index, double defaultValue) {
    ASSERT_VALID_INDEX(index, globalParameters);
    globalParameters[index].defaultValue = defaultValue;
}

void CustomTorsionForce::addEnergyParameterDerivative(const string& name) {
    for (int i = 0; i < globalParameters.size(); i++)
        if (name == globalParameters[i].name) {
            energyParameterDerivatives.push_back(i);
            return;
        }
    throw OpenMMException(string("addEnergyParameterDerivative: Unknown global parameter '"+name+"'"));
}

const string& CustomTorsionForce::getEnergyParameterDerivativeName(int index) const {
    ASSERT_VALID_INDEX(index, energyParameterDerivatives);
    return globalParameters[energyParameterDerivatives[index]].name;
}

int CustomTorsionForce::addTorsion(int particle1, int particle2, int particle3, int particle4, const vector<double>& parameters) {
    torsions.push_back(TorsionInfo(particle1, particle2, particle3, particle4, parameters));
    return torsions.size()-1;
}

void CustomTorsionForce::getTorsionParameters(int index, int& particle1, int& particle2, int& particle3, int& particle4, std::vector<double>& parameters) const {
    ASSERT_VALID_INDEX(index, torsions);
    particle1 = torsions[index].particle1;
    particle2 = torsions[index].particle2;
    particle3 = torsions[index].particle3;
    particle4 = torsions[index].particle4;
    parameters = torsions[index].parameters;
}

void CustomTorsionForce::setTorsionParameters(int index, int particle1, int particle2, int particle3, int particle4, const vector<double>& parameters) {
    ASSERT_VALID_INDEX(index, torsions);
    torsions[index].parameters = parameters;
    torsions[index].particle1 = particle1;
    torsions[index].particle2 = particle2;
    torsions[index].particle3 = particle3;
    torsions[index].particle4 = particle4;
}

ForceImpl* CustomTorsionForce::createImpl() const {
    return new CustomTorsionForceImpl(*this);
}

void CustomTorsionForce::updateParametersInContext(Context& context) {
    dynamic_cast<CustomTorsionForceImpl&>(getImplInContext(context)).updateParametersInContext(getContextImpl(context));
}

void CustomTorsionForce::setUsesPeriodicBoundaryConditions(bool periodic) {
    usePeriodic = periodic;
}

bool CustomTorsionForce::usesPeriodicBoundaryConditions() const {
    return usePeriodic;
}
