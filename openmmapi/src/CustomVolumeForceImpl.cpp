/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2025 Stanford University and the Authors.           *
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

#include "openmm/internal/CustomVolumeForceImpl.h"
#include "openmm/internal/ContextImpl.h"
#include "lepton/ParsedExpression.h"
#include "lepton/Parser.h"

using namespace OpenMM;
using namespace std;

CustomVolumeForceImpl::CustomVolumeForceImpl(const CustomVolumeForce& owner) : CustomCPPForceImpl(owner), owner(owner) {
    Lepton::ParsedExpression expression = Lepton::Parser::parse(owner.getEnergyFunction());
    energyExpression = expression.createCompiledExpression();
    map<string, double*> variableLocations;
    globalParameterNames.resize(owner.getNumGlobalParameters());
    globalValues.resize(owner.getNumGlobalParameters());
    for (int i = 0; i < owner.getNumGlobalParameters(); i++) {
        string name = owner.getGlobalParameterName(i);
        defaultParameters[name] = owner.getGlobalParameterDefaultValue(i);
        globalParameterNames[i] = name;
        variableLocations[name] = &globalValues[i];
    }
    variableLocations["v"] = &volume;
    variableLocations["ax"] = &a[0];
    variableLocations["bx"] = &b[0];
    variableLocations["by"] = &b[1];
    variableLocations["cx"] = &c[0];
    variableLocations["cy"] = &c[1];
    variableLocations["cz"] = &c[2];
    energyExpression.setVariableLocations(variableLocations);
}

void CustomVolumeForceImpl::initialize(ContextImpl& context) {
    CustomCPPForceImpl::initialize(context);
}

double CustomVolumeForceImpl::computeForce(ContextImpl& context, const vector<Vec3>& positions, vector<Vec3>& forces) {
    for (int i = 0; i < globalParameterNames.size(); i++)
        globalValues[i] = context.getParameter(globalParameterNames[i]);
    context.getPeriodicBoxVectors(a, b, c);
    volume = a[0]*b[1]*c[2];
    return energyExpression.evaluate();
}

map<string, double> CustomVolumeForceImpl::getDefaultParameters() {
    return defaultParameters;
}
