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
#include "openmm/CustomExternalForce.h"
#include "openmm/internal/AssertionUtilities.h"
#include "openmm/internal/CustomExternalForceImpl.h"
#include <cmath>
#include <map>
#include <sstream>
#include <utility>

using namespace OpenMM;
using namespace std;

CustomExternalForce::CustomExternalForce(const string& energy) : energyExpression(energy), numContexts(0) {
}

const string& CustomExternalForce::getEnergyFunction() const {
    return energyExpression;
}

void CustomExternalForce::setEnergyFunction(const std::string& energy) {
    energyExpression = energy;
}

int CustomExternalForce::addPerParticleParameter(const string& name) {
    parameters.push_back(ParticleParameterInfo(name));
    return parameters.size()-1;
}

const string& CustomExternalForce::getPerParticleParameterName(int index) const {
    ASSERT_VALID_INDEX(index, parameters);
    return parameters[index].name;
}

void CustomExternalForce::setPerParticleParameterName(int index, const string& name) {
    ASSERT_VALID_INDEX(index, parameters);
    parameters[index].name = name;
}

int CustomExternalForce::addGlobalParameter(const string& name, double defaultValue) {
    globalParameters.push_back(GlobalParameterInfo(name, defaultValue));
    return globalParameters.size()-1;
}

const string& CustomExternalForce::getGlobalParameterName(int index) const {
    ASSERT_VALID_INDEX(index, globalParameters);
    return globalParameters[index].name;
}

void CustomExternalForce::setGlobalParameterName(int index, const string& name) {
    ASSERT_VALID_INDEX(index, globalParameters);
    globalParameters[index].name = name;
}

double CustomExternalForce::getGlobalParameterDefaultValue(int index) const {
    ASSERT_VALID_INDEX(index, globalParameters);
    return globalParameters[index].defaultValue;
}

void CustomExternalForce::setGlobalParameterDefaultValue(int index, double defaultValue) {
    ASSERT_VALID_INDEX(index, globalParameters);
    globalParameters[index].defaultValue = defaultValue;
}

int CustomExternalForce::addParticle(int particle, const vector<double>& parameters) {
    particles.push_back(ParticleInfo(particle, parameters));
    return particles.size()-1;
}

void CustomExternalForce::getParticleParameters(int index, int& particle, std::vector<double>& parameters) const {
    ASSERT_VALID_INDEX(index, particles);
    particle = particles[index].particle;
    parameters = particles[index].parameters;
}

void CustomExternalForce::setParticleParameters(int index, int particle, const vector<double>& parameters) {
    ASSERT_VALID_INDEX(index, particles);
    particles[index].parameters = parameters;
    particles[index].particle = particle;
    if (numContexts > 0) {
        firstChangedParticle = min(index, firstChangedParticle);
        lastChangedParticle = max(index, lastChangedParticle);
    }
}

ForceImpl* CustomExternalForce::createImpl() const {
    if (numContexts == 0) {
        // Begin tracking changes to particles.
        firstChangedParticle = particles.size();
        lastChangedParticle = -1;
    }
    numContexts++;    return new CustomExternalForceImpl(*this);
}

void CustomExternalForce::updateParametersInContext(Context& context) {
    dynamic_cast<CustomExternalForceImpl&>(getImplInContext(context)).updateParametersInContext(getContextImpl(context), firstChangedParticle, lastChangedParticle);
    if (numContexts == 1) {
        // We just updated the only existing context for this force, so we can reset
        // the tracking of changed particles.
        firstChangedParticle = particles.size();
        lastChangedParticle = -1;
    }
}

bool CustomExternalForce::usesPeriodicBoundaryConditions() const {
    return (energyExpression.find("periodicdistance") != string::npos);
}
