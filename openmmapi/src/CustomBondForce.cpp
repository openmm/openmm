/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2008-2024 Stanford University and the Authors.      *
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
#include "openmm/CustomBondForce.h"
#include "openmm/internal/AssertionUtilities.h"
#include "openmm/internal/CustomBondForceImpl.h"
#include <cmath>
#include <map>
#include <sstream>
#include <utility>

using namespace OpenMM;
using namespace std;

CustomBondForce::CustomBondForce(const string& energy) : energyExpression(energy), usePeriodic(false), numContexts(0) {
}

const string& CustomBondForce::getEnergyFunction() const {
    return energyExpression;
}

void CustomBondForce::setEnergyFunction(const std::string& energy) {
    energyExpression = energy;
}

int CustomBondForce::addPerBondParameter(const string& name) {
    parameters.push_back(BondParameterInfo(name));
    return parameters.size()-1;
}

const string& CustomBondForce::getPerBondParameterName(int index) const {
    ASSERT_VALID_INDEX(index, parameters);
    return parameters[index].name;
}

void CustomBondForce::setPerBondParameterName(int index, const string& name) {
    ASSERT_VALID_INDEX(index, parameters);
    parameters[index].name = name;
}

int CustomBondForce::addGlobalParameter(const string& name, double defaultValue) {
    globalParameters.push_back(GlobalParameterInfo(name, defaultValue));
    return globalParameters.size()-1;
}

const string& CustomBondForce::getGlobalParameterName(int index) const {
    ASSERT_VALID_INDEX(index, globalParameters);
    return globalParameters[index].name;
}

void CustomBondForce::setGlobalParameterName(int index, const string& name) {
    ASSERT_VALID_INDEX(index, globalParameters);
    globalParameters[index].name = name;
}

double CustomBondForce::getGlobalParameterDefaultValue(int index) const {
    ASSERT_VALID_INDEX(index, globalParameters);
    return globalParameters[index].defaultValue;
}

void CustomBondForce::setGlobalParameterDefaultValue(int index, double defaultValue) {
    ASSERT_VALID_INDEX(index, globalParameters);
    globalParameters[index].defaultValue = defaultValue;
}

void CustomBondForce::addEnergyParameterDerivative(const string& name) {
    for (int i = 0; i < globalParameters.size(); i++)
        if (name == globalParameters[i].name) {
            energyParameterDerivatives.push_back(i);
            return;
        }
    throw OpenMMException(string("addEnergyParameterDerivative: Unknown global parameter '"+name+"'"));
}

const string& CustomBondForce::getEnergyParameterDerivativeName(int index) const {
    ASSERT_VALID_INDEX(index, energyParameterDerivatives);
    return globalParameters[energyParameterDerivatives[index]].name;
}

int CustomBondForce::addBond(int particle1, int particle2, const vector<double>& parameters) {
    bonds.push_back(BondInfo(particle1, particle2, parameters));
    return bonds.size()-1;
}

void CustomBondForce::getBondParameters(int index, int& particle1, int& particle2, std::vector<double>& parameters) const {
    ASSERT_VALID_INDEX(index, bonds);
    particle1 = bonds[index].particle1;
    particle2 = bonds[index].particle2;
    parameters = bonds[index].parameters;
}

void CustomBondForce::setBondParameters(int index, int particle1, int particle2, const vector<double>& parameters) {
    ASSERT_VALID_INDEX(index, bonds);
    bonds[index].parameters = parameters;
    bonds[index].particle1 = particle1;
    bonds[index].particle2 = particle2;
    if (numContexts > 0) {
        firstChangedBond = min(index, firstChangedBond);
        lastChangedBond = max(index, lastChangedBond);
    }
}

ForceImpl* CustomBondForce::createImpl() const {
    if (numContexts == 0) {
        // Begin tracking changes to bonds.
        firstChangedBond = bonds.size();
        lastChangedBond = -1;
    }
    numContexts++;
    return new CustomBondForceImpl(*this);
}

void CustomBondForce::updateParametersInContext(Context& context) {
    dynamic_cast<CustomBondForceImpl&>(getImplInContext(context)).updateParametersInContext(getContextImpl(context), firstChangedBond, lastChangedBond);
    if (numContexts == 1) {
        // We just updated the only existing context for this force, so we can reset
        // the tracking of changed bonds.
        firstChangedBond = bonds.size();
        lastChangedBond = -1;
    }
}

void CustomBondForce::setUsesPeriodicBoundaryConditions(bool periodic) {
    usePeriodic = periodic;
}

bool CustomBondForce::usesPeriodicBoundaryConditions() const {
    return usePeriodic;
}
