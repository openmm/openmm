/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2008-2012 Stanford University and the Authors.      *
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
#include "openmm/OpenMMException.h"
#include "openmm/internal/ContextImpl.h"
#include "openmm/kernels.h"
#include <limits>
#include <string>

using namespace OpenMM;
using std::string;
using std::vector;

VariableVerletIntegrator::VariableVerletIntegrator(double errorTol) : errorTol(errorTol) {
    setConstraintTolerance(1e-5);
    setStepSize(0.0);
}

void VariableVerletIntegrator::initialize(ContextImpl& contextRef) {
    if (owner != NULL && &contextRef.getOwner() != owner)
        throw OpenMMException("This Integrator is already bound to a context");
    context = &contextRef;
    owner = &contextRef.getOwner();
    kernel = context->getPlatform().createKernel(IntegrateVariableVerletStepKernel::Name(), contextRef);
    kernel.getAs<IntegrateVariableVerletStepKernel>().initialize(contextRef.getSystem(), *this);
}

void VariableVerletIntegrator::cleanup() {
    kernel = Kernel();
}

vector<string> VariableVerletIntegrator::getKernelNames() {
    std::vector<std::string> names;
    names.push_back(IntegrateVariableVerletStepKernel::Name());
    return names;
}

double VariableVerletIntegrator::computeKineticEnergy() {
    return kernel.getAs<IntegrateVariableVerletStepKernel>().computeKineticEnergy(*context, *this);
}

void VariableVerletIntegrator::step(int steps) {
    if (context == NULL)
        throw OpenMMException("This Integrator is not bound to a context!");
    for (int i = 0; i < steps; ++i) {
        context->updateContextState();
        context->calcForcesAndEnergy(true, false);
        setStepSize(kernel.getAs<IntegrateVariableVerletStepKernel>().execute(*context, *this, std::numeric_limits<double>::infinity()));
    }
}

void VariableVerletIntegrator::stepTo(double time) {
    if (context == NULL)
        throw OpenMMException("This Integrator is not bound to a context!");  
    while (time > context->getTime()) {
        context->updateContextState();
        context->calcForcesAndEnergy(true, false);
        setStepSize(kernel.getAs<IntegrateVariableVerletStepKernel>().execute(*context, *this, time));
    }
}
