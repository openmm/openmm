/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2008-2009 Stanford University and the Authors.      *
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

#include "openmm/VariableVerletIntegrator.h"
#include "openmm/Context.h"
#include "openmm/internal/ContextImpl.h"
#include "openmm/kernels.h"
#include <limits>
#include <string>

using namespace OpenMM;
using std::string;
using std::vector;

VariableVerletIntegrator::VariableVerletIntegrator(double errorTol) : errorTol(errorTol) {
    setConstraintTolerance(1e-4);
}

void VariableVerletIntegrator::initialize(ContextImpl& contextRef) {
    context = &contextRef;
    kernel = context->getPlatform().createKernel(IntegrateVariableVerletStepKernel::Name(), contextRef);
    dynamic_cast<IntegrateVariableVerletStepKernel&>(kernel.getImpl()).initialize(contextRef.getSystem(), *this);
}

vector<string> VariableVerletIntegrator::getKernelNames() {
    std::vector<std::string> names;
    names.push_back(IntegrateVariableVerletStepKernel::Name());
    return names;
}

void VariableVerletIntegrator::step(int steps) {
    for (int i = 0; i < steps; ++i) {
        context->updateContextState();
        context->calcForcesAndEnergy(true, false);
        dynamic_cast<IntegrateVariableVerletStepKernel&>(kernel.getImpl()).execute(*context, *this, std::numeric_limits<double>::infinity());
    }
}

void VariableVerletIntegrator::stepTo(double time) {
    while (time > context->getTime()) {
        context->updateContextState();
        context->calcForcesAndEnergy(true, false);
        dynamic_cast<IntegrateVariableVerletStepKernel&>(kernel.getImpl()).execute(*context, *this, time);
    }
}
