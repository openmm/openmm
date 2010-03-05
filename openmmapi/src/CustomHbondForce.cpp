/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2008-2010 Stanford University and the Authors.      *
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
#include "openmm/CustomHbondForce.h"
#include "openmm/internal/CustomHbondForceImpl.h"
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

CustomHbondForce::CustomHbondForce(const string& energy) : energyExpression(energy), nonbondedMethod(NoCutoff), cutoffDistance(1.0) {
}

const string& CustomHbondForce::getEnergyFunction() const {
    return energyExpression;
}

void CustomHbondForce::setEnergyFunction(const std::string& energy) {
    energyExpression = energy;
}

CustomHbondForce::NonbondedMethod CustomHbondForce::getNonbondedMethod() const {
    return nonbondedMethod;
}

void CustomHbondForce::setNonbondedMethod(NonbondedMethod method) {
    nonbondedMethod = method;
}

double CustomHbondForce::getCutoffDistance() const {
    return cutoffDistance;
}

void CustomHbondForce::setCutoffDistance(double distance) {
    cutoffDistance = distance;
}

int CustomHbondForce::addPerDonorParameter(const string& name) {
    donorParameters.push_back(PerPairParameterInfo(name));
    return donorParameters.size()-1;
}

const string& CustomHbondForce::getPerDonorParameterName(int index) const {
    return donorParameters[index].name;
}

void CustomHbondForce::setPerDonorParameterName(int index, const string& name) {
    donorParameters[index].name = name;
}

int CustomHbondForce::addPerAcceptorParameter(const string& name) {
    acceptorParameters.push_back(PerPairParameterInfo(name));
    return acceptorParameters.size()-1;
}

const string& CustomHbondForce::getPerAcceptorParameterName(int index) const {
    return acceptorParameters[index].name;
}

void CustomHbondForce::setPerAcceptorParameterName(int index, const string& name) {
    acceptorParameters[index].name = name;
}

int CustomHbondForce::addGlobalParameter(const string& name, double defaultValue) {
    globalParameters.push_back(GlobalParameterInfo(name, defaultValue));
    return globalParameters.size()-1;
}

const string& CustomHbondForce::getGlobalParameterName(int index) const {
    return globalParameters[index].name;
}

void CustomHbondForce::setGlobalParameterName(int index, const string& name) {
    globalParameters[index].name = name;
}

double CustomHbondForce::getGlobalParameterDefaultValue(int index) const {
    return globalParameters[index].defaultValue;
}

void CustomHbondForce::setGlobalParameterDefaultValue(int index, double defaultValue) {
    globalParameters[index].defaultValue = defaultValue;
}

int CustomHbondForce::addDonor(int d1, int d2, int d3, const vector<double>& parameters) {
    donors.push_back(GroupInfo(d1, d2, d3, parameters));
    return donors.size()-1;
}

void CustomHbondForce::getDonorParameters(int index, int& d1, int& d2, int&  d3, std::vector<double>& parameters) const {
    d1 = donors[index].p1;
    d2 = donors[index].p2;
    d3 = donors[index].p3;
    parameters = donors[index].parameters;
}

void CustomHbondForce::setDonorParameters(int index, int d1, int d2, int d3, const vector<double>& parameters) {
    donors[index].p1 = d1;
    donors[index].p2 = d2;
    donors[index].p3 = d3;
    donors[index].parameters = parameters;
}

int CustomHbondForce::addAcceptor(int a1, int a2, int a3, const vector<double>& parameters) {
    acceptors.push_back(GroupInfo(a1, a2, a3, parameters));
    return acceptors.size()-1;
}

void CustomHbondForce::getAcceptorParameters(int index, int& a1, int& a2, int& a3, std::vector<double>& parameters) const {
    a1 = acceptors[index].p1;
    a2 = acceptors[index].p2;
    a3 = acceptors[index].p3;
    parameters = acceptors[index].parameters;
}

void CustomHbondForce::setAcceptorParameters(int index, int a1, int a2, int a3, const vector<double>& parameters) {
    acceptors[index].p1 = a1;
    acceptors[index].p2 = a2;
    acceptors[index].p3 = a3;
    acceptors[index].parameters = parameters;
}

int CustomHbondForce::addExclusion(int donor, int acceptor) {
    exclusions.push_back(ExclusionInfo(donor, acceptor));
    return exclusions.size()-1;
}
void CustomHbondForce::getExclusionParticles(int index, int& donor, int& acceptor) const {
    donor = exclusions[index].donor;
    acceptor = exclusions[index].acceptor;
}

void CustomHbondForce::setExclusionParticles(int index, int donor, int acceptor) {
    exclusions[index].donor = donor;
    exclusions[index].acceptor = acceptor;
}

int CustomHbondForce::addFunction(const std::string& name, const std::vector<double>& values, double min, double max, bool interpolating) {
    if (max <= min)
        throw OpenMMException("CustomHbondForce: max <= min for a tabulated function.");
    if (values.size() < 2)
        throw OpenMMException("CustomHbondForce: a tabulated function must have at least two points");
    functions.push_back(FunctionInfo(name, values, min, max, interpolating));
    return functions.size()-1;
}

void CustomHbondForce::getFunctionParameters(int index, std::string& name, std::vector<double>& values, double& min, double& max, bool& interpolating) const {
    name = functions[index].name;
    values = functions[index].values;
    min = functions[index].min;
    max = functions[index].max;
    interpolating = functions[index].interpolating;
}

void CustomHbondForce::setFunctionParameters(int index, const std::string& name, const std::vector<double>& values, double min, double max, bool interpolating) {
    if (max <= min)
        throw OpenMMException("CustomHbondForce: max <= min for a tabulated function.");
    if (values.size() < 2)
        throw OpenMMException("CustomHbondForce: a tabulated function must have at least two points");
    functions[index].name = name;
    functions[index].values = values;
    functions[index].min = min;
    functions[index].max = max;
    functions[index].interpolating = interpolating;
}

ForceImpl* CustomHbondForce::createImpl() {
    return new CustomHbondForceImpl(*this);
}

