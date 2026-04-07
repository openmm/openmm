/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit.                   *
 * See https://openmm.org/development.                                        *
 *                                                                            *
 * Portions copyright (c) 2008-2020 Stanford University and the Authors.      *
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

#include "openmm/VerletIntegrator.h"
#include "openmm/Context.h"
#include "openmm/OpenMMException.h"
#include "openmm/System.h"
#include "openmm/internal/ContextImpl.h"
#include "openmm/internal/ForceImpl.h"
#include "openmm/kernels.h"
#include <string>

using namespace OpenMM;
using std::string;
using std::vector;

namespace {

bool systemUsesVerletPart1ContextUpdate(const System& system) {
    for (int i = 0; i < system.getNumForces(); i++) {
        if (system.getForce(i).usesVerletPart1ContextUpdate())
            return true;
    }
    return false;
}

bool contextHasApplyBussiThermostatKernel(ContextImpl& impl) {
    for (ForceImpl* fi : impl.getForceImpls()) {
        for (const string& kn : fi->getKernelNames()) {
            if (kn == ApplyBussiThermostatKernel::Name())
                return true;
        }
    }
    return false;
}

} // namespace

VerletIntegrator::VerletIntegrator(double stepSize) : useVerletPart1ContextUpdate_(false) {
    setStepSize(stepSize);
    setConstraintTolerance(1e-5);
}

void VerletIntegrator::initialize(ContextImpl& contextRef) {
    if (owner != NULL && &contextRef.getOwner() != owner)
        throw OpenMMException("This Integrator is already bound to a context");
    context = &contextRef;
    owner = &contextRef.getOwner();
    kernel = context->getPlatform().createKernel(IntegrateVerletStepKernel::Name(), contextRef);
    kernel.getAs<IntegrateVerletStepKernel>().initialize(contextRef.getSystem(), *this);
    useVerletPart1ContextUpdate_ = systemUsesVerletPart1ContextUpdate(contextRef.getSystem()) ||
            contextHasApplyBussiThermostatKernel(contextRef);
}

void VerletIntegrator::cleanup() {
    kernel = Kernel();
    useVerletPart1ContextUpdate_ = false;
}

vector<string> VerletIntegrator::getKernelNames() {
    return {IntegrateVerletStepKernel::Name()};
}

double VerletIntegrator::computeKineticEnergy() {
    return kernel.getAs<IntegrateVerletStepKernel>().computeKineticEnergy(*context, *this);
}

void VerletIntegrator::step(int steps) {
    if (context == NULL)
        throw OpenMMException("This Integrator is not bound to a context!");
    IntegrateVerletStepKernel& verletKernel = kernel.getAs<IntegrateVerletStepKernel>();
    for (int i = 0; i < steps; ++i) {
        if (useVerletPart1ContextUpdate_) {
            context->calcForcesAndEnergy(true, false, getIntegrationForceGroups());
            verletKernel.executePart1(*context, *this);
            context->setStepPhase(ContextImpl::STEP_PHASE_AFTER_VERLET_PART1);
            context->updateContextState();
            context->setStepPhase(ContextImpl::STEP_PHASE_NONE);
            verletKernel.executePart2(*context, *this);
        }
        else {
            context->updateContextState();
            context->calcForcesAndEnergy(true, false, getIntegrationForceGroups());
            verletKernel.execute(*context, *this);
        }
    }
}
