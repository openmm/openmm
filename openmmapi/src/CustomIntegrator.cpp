/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2011 Stanford University and the Authors.           *
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
#include "openmm/internal/ContextImpl.h"
#include "openmm/kernels.h"
#include <ctime>
#include <string>

using namespace OpenMM;
using std::string;
using std::vector;

CustomIntegrator::CustomIntegrator(double stepSize) : owner(NULL), globalsAreCurrent(true), forcesAreValid(false) {
    setStepSize(stepSize);
    setConstraintTolerance(1e-4);
    setRandomNumberSeed((int) time(NULL));
}

void CustomIntegrator::initialize(ContextImpl& contextRef) {
    if (owner != NULL && &contextRef.getOwner() != owner)
        throw OpenMMException("This Integrator is already bound to a context");
    context = &contextRef;
    owner = &contextRef.getOwner();
    kernel = context->getPlatform().createKernel(IntegrateCustomStepKernel::Name(), contextRef);
    dynamic_cast<IntegrateCustomStepKernel&>(kernel.getImpl()).initialize(contextRef.getSystem(), *this);
    dynamic_cast<IntegrateCustomStepKernel&>(kernel.getImpl()).setGlobalVariables(contextRef, globalValues);
    for (int i = 0; i < (int) perDofValues.size(); i++)
        dynamic_cast<IntegrateCustomStepKernel&>(kernel.getImpl()).setPerDofVariable(contextRef, i, perDofValues[i]);
}

void CustomIntegrator::stateChanged(State::DataType changed) {
    forcesAreValid = false;
}

vector<string> CustomIntegrator::getKernelNames() {
    std::vector<std::string> names;
    names.push_back(IntegrateCustomStepKernel::Name());
    return names;
}

void CustomIntegrator::step(int steps) {
    globalsAreCurrent = false;
    for (int i = 0; i < steps; ++i) {
        dynamic_cast<IntegrateCustomStepKernel&>(kernel.getImpl()).execute(*context, *this, forcesAreValid);
    }
}

int CustomIntegrator::addGlobalVariable(const std::string& name, double initialValue) {
    if (owner != NULL)
        throw OpenMMException("The integrator cannot be modified after it is bound to a context");
    globalNames.push_back(name);
    globalValues.push_back(initialValue);
    return globalNames.size()-1;
}

std::string CustomIntegrator::getGlobalVariableName(int index) const {
    if (index < 0 || index >= globalNames.size())
        throw OpenMMException("Index out of range");
    return globalNames[index];
}

int CustomIntegrator::addPerDofVariable(const std::string& name, double initialValue) {
    if (owner != NULL)
        throw OpenMMException("The integrator cannot be modified after it is bound to a context");
    perDofNames.push_back(name);
    perDofValues.push_back(vector<Vec3>(1, Vec3(initialValue, initialValue, initialValue)));
    return perDofNames.size()-1;
}

std::string CustomIntegrator::getPerDofVariableName(int index) const {
    if (index < 0 || index >= perDofNames.size())
        throw OpenMMException("Index out of range");
    return perDofNames[index];
}

double CustomIntegrator::getGlobalVariable(int index) const {
    if (index < 0 || index >= globalNames.size())
        throw OpenMMException("Index out of range");
    if (owner != NULL && !globalsAreCurrent) {
        dynamic_cast<const IntegrateCustomStepKernel&>(kernel.getImpl()).getGlobalVariables(*context, globalValues);
        globalsAreCurrent = true;
    }
    return globalValues[index];
}

void CustomIntegrator::setGlobalVariable(int index, double value) {
    if (index < 0 || index >= globalNames.size())
        throw OpenMMException("Index out of range");
    if (owner != NULL && !globalsAreCurrent) {
        dynamic_cast<IntegrateCustomStepKernel&>(kernel.getImpl()).getGlobalVariables(*context, globalValues);
        globalsAreCurrent = true;
    }
    globalValues[index] = value;
    dynamic_cast<IntegrateCustomStepKernel&>(kernel.getImpl()).setGlobalVariables(*context, globalValues);
}

void CustomIntegrator::getPerDofVariable(int index, std::vector<Vec3>& values) const {
    if (index < 0 || index >= perDofNames.size())
        throw OpenMMException("Index out of range");
    if (owner == NULL)
        values = perDofValues[index];
    else
        dynamic_cast<const IntegrateCustomStepKernel&>(kernel.getImpl()).getPerDofVariable(*context, index, values);
}

void CustomIntegrator::setPerDofVariable(int index, const std::vector<Vec3>& values) {
    if (index < 0 || index >= perDofNames.size())
        throw OpenMMException("Index out of range");
    if (owner == NULL)
        perDofValues[index] = values;
    else
        dynamic_cast<IntegrateCustomStepKernel&>(kernel.getImpl()).setPerDofVariable(*context, index, values);
}

int CustomIntegrator::addComputeGlobal(const std::string& variable, const std::string& expression) {
    if (owner != NULL)
        throw OpenMMException("The integrator cannot be modified after it is bound to a context");
    computations.push_back(ComputationInfo(ComputeGlobal, variable, expression));
    return computations.size()-1;
}

int CustomIntegrator::addComputePerDof(const std::string& variable, const std::string& expression) {
    if (owner != NULL)
        throw OpenMMException("The integrator cannot be modified after it is bound to a context");
    computations.push_back(ComputationInfo(ComputePerDof, variable, expression));
    return computations.size()-1;
}

int CustomIntegrator::addComputeSum(const std::string& variable, const std::string& expression) {
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

void CustomIntegrator::getComputationStep(int index, ComputationType& type, std::string& variable, std::string& expression) const {
    if (index < 0 || index >= computations.size())
        throw OpenMMException("Index out of range");
    type = computations[index].type;
    variable = computations[index].variable;
    expression = computations[index].expression;
}
