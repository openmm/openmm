/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2011-2017 Stanford University and the Authors.      *
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

#include "openmm/CustomIntegrator.h"
#include "openmm/Context.h"
#include "openmm/OpenMMException.h"
#include "openmm/internal/AssertionUtilities.h"
#include "openmm/internal/ContextImpl.h"
#include "openmm/kernels.h"
#include <set>
#include <string>

using namespace OpenMM;
using namespace std;

CustomIntegrator::CustomIntegrator(double stepSize) : globalsAreCurrent(true), forcesAreValid(false) {
    setStepSize(stepSize);
    setConstraintTolerance(1e-5);
    setRandomNumberSeed(0);
    kineticEnergy = "m*v*v/2";
}

CustomIntegrator::~CustomIntegrator() {
    for (auto function : functions)
        delete function.function;
}

void CustomIntegrator::initialize(ContextImpl& contextRef) {
    if (owner != NULL && &contextRef.getOwner() != owner)
        throw OpenMMException("This Integrator is already bound to a context");
    vector<std::string> variableList;
    set<std::string> variableSet;
    variableList.insert(variableList.end(), globalNames.begin(), globalNames.end());
    variableList.insert(variableList.end(), perDofNames.begin(), perDofNames.end());
    for (auto& name : variableList) {
        if (variableSet.find(name) != variableSet.end())
            throw OpenMMException("The Integrator defines two variables with the same name: "+name);
        variableSet.insert(name);
        if (contextRef.getParameters().find(name) != contextRef.getParameters().end())
            throw OpenMMException("The Integrator defines a variable with the same name as a Context parameter: "+name);
    }
    set<std::string> globalTargets;
    globalTargets.insert(globalNames.begin(), globalNames.end());
    globalTargets.insert("dt");
    for (auto& param : contextRef.getParameters())
        globalTargets.insert(param.first);
    for (int i = 0; i < computations.size(); i++) {
        if (computations[i].type == ComputeGlobal && globalTargets.find(computations[i].variable) == globalTargets.end())
            throw OpenMMException("Unknown global variable: "+computations[i].variable);
    }
    context = &contextRef;
    owner = &contextRef.getOwner();
    kernel = context->getPlatform().createKernel(IntegrateCustomStepKernel::Name(), contextRef);
    kernel.getAs<IntegrateCustomStepKernel>().initialize(contextRef.getSystem(), *this);
    kernel.getAs<IntegrateCustomStepKernel>().setGlobalVariables(contextRef, globalValues);
    for (int i = 0; i < (int) perDofValues.size(); i++) {
        if (perDofValues[i].size() == 1)
            perDofValues[i].resize(context->getSystem().getNumParticles(), perDofValues[i][0]);
        kernel.getAs<IntegrateCustomStepKernel>().setPerDofVariable(contextRef, i, perDofValues[i]);
    }
}

void CustomIntegrator::cleanup() {
    kernel = Kernel();
}

void CustomIntegrator::stateChanged(State::DataType changed) {
    forcesAreValid = false;
}

vector<string> CustomIntegrator::getKernelNames() {
    vector<string> names;
    names.push_back(IntegrateCustomStepKernel::Name());
    return names;
}

double CustomIntegrator::computeKineticEnergy() {
    return kernel.getAs<IntegrateCustomStepKernel>().computeKineticEnergy(*context, *this, forcesAreValid);
}

void CustomIntegrator::step(int steps) {
    if (context == NULL)
        throw OpenMMException("This Integrator is not bound to a context!");  
    globalsAreCurrent = false;
    for (int i = 0; i < steps; ++i) {
        kernel.getAs<IntegrateCustomStepKernel>().execute(*context, *this, forcesAreValid);
    }
}

int CustomIntegrator::addGlobalVariable(const string& name, double initialValue) {
    if (owner != NULL)
        throw OpenMMException("The integrator cannot be modified after it is bound to a context");
    globalNames.push_back(name);
    globalValues.push_back(initialValue);
    return globalNames.size()-1;
}

const string& CustomIntegrator::getGlobalVariableName(int index) const {
    ASSERT_VALID_INDEX(index, globalNames);
    return globalNames[index];
}

int CustomIntegrator::addPerDofVariable(const string& name, double initialValue) {
    if (owner != NULL)
        throw OpenMMException("The integrator cannot be modified after it is bound to a context");
    perDofNames.push_back(name);
    perDofValues.push_back(vector<Vec3>(1, Vec3(initialValue, initialValue, initialValue)));
    return perDofNames.size()-1;
}

const string& CustomIntegrator::getPerDofVariableName(int index) const {
    ASSERT_VALID_INDEX(index, perDofNames);
    return perDofNames[index];
}

double CustomIntegrator::getGlobalVariable(int index) const {
    ASSERT_VALID_INDEX(index, globalValues);
    if (owner != NULL && !globalsAreCurrent) {
        kernel.getAs<const IntegrateCustomStepKernel>().getGlobalVariables(*context, globalValues);
        globalsAreCurrent = true;
    }
    return globalValues[index];
}

double CustomIntegrator::getGlobalVariableByName(const string& name) const {
    for (int i = 0; i < (int) globalNames.size(); i++) {
        if (name == globalNames[i]) {
            return getGlobalVariable(i);
        }
    }
    throw OpenMMException("Illegal global variable name: "+name);
}

void CustomIntegrator::setGlobalVariable(int index, double value) {
    ASSERT_VALID_INDEX(index, globalValues);
    if (owner != NULL && !globalsAreCurrent) {
        kernel.getAs<IntegrateCustomStepKernel>().getGlobalVariables(*context, globalValues);
        globalsAreCurrent = true;
    }
    globalValues[index] = value;
    if (owner != NULL)
        kernel.getAs<IntegrateCustomStepKernel>().setGlobalVariables(*context, globalValues);
}

void CustomIntegrator::setGlobalVariableByName(const string& name, double value) {
    for (int i = 0; i < (int) globalNames.size(); i++)
        if (name == globalNames[i]) {
            setGlobalVariable(i, value);
            return;
        }
    throw OpenMMException("Illegal global variable name: "+name);
}

void CustomIntegrator::getPerDofVariable(int index, vector<Vec3>& values) const {
    ASSERT_VALID_INDEX(index, perDofValues);
    if (owner == NULL)
        values = perDofValues[index];
    else
        kernel.getAs<const IntegrateCustomStepKernel>().getPerDofVariable(*context, index, values);
}

void CustomIntegrator::getPerDofVariableByName(const string& name,  vector<Vec3>& values) const {
    for (int i = 0; i < (int) perDofNames.size(); i++) {
        if (name == perDofNames[i]) {
            getPerDofVariable(i, values);
            return;
        }
    }
    throw OpenMMException("Illegal per-DOF variable name: "+name);
}

void CustomIntegrator::setPerDofVariable(int index, const vector<Vec3>& values) {
    ASSERT_VALID_INDEX(index, perDofValues);
    if (owner != NULL && values.size() != context->getSystem().getNumParticles())
        throw OpenMMException("setPerDofVariable() called with wrong number of values");
    if (owner == NULL)
        perDofValues[index] = values;
    else
        kernel.getAs<IntegrateCustomStepKernel>().setPerDofVariable(*context, index, values);
}

void CustomIntegrator::setPerDofVariableByName(const string& name, const vector<Vec3>& value) {
    for (int i = 0; i < (int) perDofNames.size(); i++)
        if (name == perDofNames[i]) {
            setPerDofVariable(i, value);
            return;
        }
    throw OpenMMException("Illegal per-DOF variable name: "+name);
}

int CustomIntegrator::addComputeGlobal(const string& variable, const string& expression) {
    if (owner != NULL)
        throw OpenMMException("The integrator cannot be modified after it is bound to a context");
    computations.push_back(ComputationInfo(ComputeGlobal, variable, expression));
    return computations.size()-1;
}

int CustomIntegrator::addComputePerDof(const string& variable, const string& expression) {
    if (owner != NULL)
        throw OpenMMException("The integrator cannot be modified after it is bound to a context");
    computations.push_back(ComputationInfo(ComputePerDof, variable, expression));
    return computations.size()-1;
}

int CustomIntegrator::addComputeSum(const string& variable, const string& expression) {
    if (owner != NULL)
        throw OpenMMException("The integrator cannot be modified after it is bound to a context");
    computations.push_back(ComputationInfo(ComputeSum, variable, expression));
    return computations.size()-1;
}

int CustomIntegrator::addConstrainPositions() {
    if (owner != NULL)
        throw OpenMMException("The integrator cannot be modified after it is bound to a context");
    computations.push_back(ComputationInfo(ConstrainPositions, "", ""));
    return computations.size()-1;
}

int CustomIntegrator::addConstrainVelocities() {
    if (owner != NULL)
        throw OpenMMException("The integrator cannot be modified after it is bound to a context");
    computations.push_back(ComputationInfo(ConstrainVelocities, "", ""));
    return computations.size()-1;
}

int CustomIntegrator::addUpdateContextState() {
    if (owner != NULL)
        throw OpenMMException("The integrator cannot be modified after it is bound to a context");
    computations.push_back(ComputationInfo(UpdateContextState, "", ""));
    return computations.size()-1;
}

int CustomIntegrator::beginIfBlock(const string& expression) {
    if (owner != NULL)
        throw OpenMMException("The integrator cannot be modified after it is bound to a context");
    computations.push_back(ComputationInfo(IfBlockStart, "", expression));
    return computations.size()-1;
}

int CustomIntegrator::beginWhileBlock(const string& expression) {
    if (owner != NULL)
        throw OpenMMException("The integrator cannot be modified after it is bound to a context");
    computations.push_back(ComputationInfo(WhileBlockStart, "", expression));
    return computations.size()-1;
}

int CustomIntegrator::endBlock() {
    if (owner != NULL)
        throw OpenMMException("The integrator cannot be modified after it is bound to a context");
    computations.push_back(ComputationInfo(BlockEnd, "", ""));
    return computations.size()-1;
}

void CustomIntegrator::getComputationStep(int index, ComputationType& type, string& variable, string& expression) const {
    ASSERT_VALID_INDEX(index, computations);
    type = computations[index].type;
    variable = computations[index].variable;
    expression = computations[index].expression;
}

int CustomIntegrator::addTabulatedFunction(const std::string& name, TabulatedFunction* function) {
    functions.push_back(FunctionInfo(name, function));
    return functions.size()-1;
}

const TabulatedFunction& CustomIntegrator::getTabulatedFunction(int index) const {
    ASSERT_VALID_INDEX(index, functions);
    return *functions[index].function;
}

TabulatedFunction& CustomIntegrator::getTabulatedFunction(int index) {
    ASSERT_VALID_INDEX(index, functions);
    return *functions[index].function;
}

const string& CustomIntegrator::getTabulatedFunctionName(int index) const {
    ASSERT_VALID_INDEX(index, functions);
    return functions[index].name;
}

const string& CustomIntegrator::getKineticEnergyExpression() const {
    return kineticEnergy;
}

void CustomIntegrator::setKineticEnergyExpression(const string& expression) {
    kineticEnergy = expression;
}
