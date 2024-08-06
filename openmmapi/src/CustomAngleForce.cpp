/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2010-2024 Stanford University and the Authors.      *
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
#include "openmm/CustomAngleForce.h"
#include "openmm/internal/AssertionUtilities.h"
#include "openmm/internal/CustomAngleForceImpl.h"
#include <cmath>
#include <map>
#include <sstream>
#include <utility>

using namespace OpenMM;
using namespace std;

CustomAngleForce::CustomAngleForce(const string& energy) : energyExpression(energy), usePeriodic(false), numContexts(0) {
}

const string& CustomAngleForce::getEnergyFunction() const {
    return energyExpression;
}

void CustomAngleForce::setEnergyFunction(const std::string& energy) {
    energyExpression = energy;
}

int CustomAngleForce::addPerAngleParameter(const string& name) {
    parameters.push_back(AngleParameterInfo(name));
    return parameters.size()-1;
}

const string& CustomAngleForce::getPerAngleParameterName(int index) const {
    ASSERT_VALID_INDEX(index, parameters);
    return parameters[index].name;
}

void CustomAngleForce::setPerAngleParameterName(int index, const string& name) {
    ASSERT_VALID_INDEX(index, parameters);
    parameters[index].name = name;
}

int CustomAngleForce::addGlobalParameter(const string& name, double defaultValue) {
    globalParameters.push_back(GlobalParameterInfo(name, defaultValue));
    return globalParameters.size()-1;
}

const string& CustomAngleForce::getGlobalParameterName(int index) const {
    ASSERT_VALID_INDEX(index, globalParameters);
    return globalParameters[index].name;
}

void CustomAngleForce::setGlobalParameterName(int index, const string& name) {
    ASSERT_VALID_INDEX(index, globalParameters);
    globalParameters[index].name = name;
}

double CustomAngleForce::getGlobalParameterDefaultValue(int index) const {
    ASSERT_VALID_INDEX(index, globalParameters);
    return globalParameters[index].defaultValue;
}

void CustomAngleForce::setGlobalParameterDefaultValue(int index, double defaultValue) {
    ASSERT_VALID_INDEX(index, globalParameters);
    globalParameters[index].defaultValue = defaultValue;
}

void CustomAngleForce::addEnergyParameterDerivative(const string& name) {
    for (int i = 0; i < globalParameters.size(); i++)
        if (name == globalParameters[i].name) {
            energyParameterDerivatives.push_back(i);
            return;
        }
    throw OpenMMException(string("addEnergyParameterDerivative: Unknown global parameter '"+name+"'"));
}

const string& CustomAngleForce::getEnergyParameterDerivativeName(int index) const {
    ASSERT_VALID_INDEX(index, energyParameterDerivatives);
    return globalParameters[energyParameterDerivatives[index]].name;
}

int CustomAngleForce::addAngle(int particle1, int particle2, int particle3, const vector<double>& parameters) {
    angles.push_back(AngleInfo(particle1, particle2, particle3, parameters));
    return angles.size()-1;
}

void CustomAngleForce::getAngleParameters(int index, int& particle1, int& particle2, int& particle3, std::vector<double>& parameters) const {
    ASSERT_VALID_INDEX(index, angles);
    particle1 = angles[index].particle1;
    particle2 = angles[index].particle2;
    particle3 = angles[index].particle3;
    parameters = angles[index].parameters;
}

void CustomAngleForce::setAngleParameters(int index, int particle1, int particle2, int particle3, const vector<double>& parameters) {
    ASSERT_VALID_INDEX(index, angles);
    angles[index].parameters = parameters;
    angles[index].particle1 = particle1;
    angles[index].particle2 = particle2;
    angles[index].particle3 = particle3;
    if (numContexts > 0) {
        firstChangedAngle = min(index, firstChangedAngle);
        lastChangedAngle = max(index, lastChangedAngle);
    }
}

ForceImpl* CustomAngleForce::createImpl() const {
    if (numContexts == 0) {
        // Begin tracking changes to angles.
        firstChangedAngle = angles.size();
        lastChangedAngle = -1;
    }
    numContexts++;
    return new CustomAngleForceImpl(*this);
}

void CustomAngleForce::updateParametersInContext(Context& context) {
    dynamic_cast<CustomAngleForceImpl&>(getImplInContext(context)).updateParametersInContext(getContextImpl(context), firstChangedAngle, lastChangedAngle);
    if (numContexts == 1) {
        // We just updated the only existing context for this force, so we can reset
        // the tracking of changed angles.
        firstChangedAngle = angles.size();
        lastChangedAngle = -1;
    }
}

void CustomAngleForce::setUsesPeriodicBoundaryConditions(bool periodic) {
    usePeriodic = periodic;
}

bool CustomAngleForce::usesPeriodicBoundaryConditions() const {
    return usePeriodic;
}
