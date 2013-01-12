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

#include "openmm/RPMDIntegrator.h"
#include "openmm/Context.h"
#include "openmm/OpenMMException.h"
#include "openmm/internal/ContextImpl.h"
#include "openmm/RpmdKernels.h"
#include <ctime>
#include <string>

using namespace OpenMM;
using std::string;
using std::vector;

RPMDIntegrator::RPMDIntegrator(int numCopies, double temperature, double frictionCoeff, double stepSize) :
        numCopies(numCopies), forcesAreValid(false), hasSetPosition(false), hasSetVelocity(false), isFirstStep(true) {
    setTemperature(temperature);
    setFriction(frictionCoeff);
    setStepSize(stepSize);
    setConstraintTolerance(1e-4);
    setRandomNumberSeed((int) time(NULL));
}

void RPMDIntegrator::initialize(ContextImpl& contextRef) {
    if (owner != NULL && &contextRef.getOwner() != owner)
        throw OpenMMException("This Integrator is already bound to a context");
    if (contextRef.getSystem().getNumConstraints() > 0)
        throw OpenMMException("RPMDIntegrator cannot be used with Systems that include constraints");
    context = &contextRef;
    owner = &contextRef.getOwner();
    kernel = context->getPlatform().createKernel(IntegrateRPMDStepKernel::Name(), contextRef);
    kernel.getAs<IntegrateRPMDStepKernel>().initialize(contextRef.getSystem(), *this);
}

void RPMDIntegrator::cleanup() {
    kernel = Kernel();
}

void RPMDIntegrator::stateChanged(State::DataType changed) {
    forcesAreValid = false;
}

vector<string> RPMDIntegrator::getKernelNames() {
    std::vector<std::string> names;
    names.push_back(IntegrateRPMDStepKernel::Name());
    return names;
}

void RPMDIntegrator::setPositions(int copy, const vector<Vec3>& positions) {
    kernel.getAs<IntegrateRPMDStepKernel>().setPositions(copy, positions);
    hasSetPosition = true;
}

void RPMDIntegrator::setVelocities(int copy, const vector<Vec3>& velocities) {
    kernel.getAs<IntegrateRPMDStepKernel>().setVelocities(copy, velocities);
    hasSetVelocity = true;
}

State RPMDIntegrator::getState(int copy, int types, bool enforcePeriodicBox, int groups) {
    if (isFirstStep) {
        // Call setPositions() on the Context so it doesn't think the user is trying to
        // run a simulation without setting positions first.  These positions will
        // immediately get overwritten by the ones stored in this integrator.
        
        vector<Vec3> p(context->getSystem().getNumParticles(), Vec3());
        context->getOwner().setPositions(p);
        isFirstStep = false;
    }
    kernel.getAs<IntegrateRPMDStepKernel>().copyToContext(copy, *context);
    return context->getOwner().getState(types, enforcePeriodicBox, groups);
}

double RPMDIntegrator::computeKineticEnergy() {
    return kernel.getAs<IntegrateRPMDStepKernel>().computeKineticEnergy(*context, *this);
}

void RPMDIntegrator::step(int steps) {
    if (!hasSetPosition) {
        // Initialize the positions from the context.
        
        State s = context->getOwner().getState(State::Positions);
        for (int i = 0; i < numCopies; i++)
            setPositions(i, s.getPositions());
    }
    if (!hasSetVelocity) {
        // Initialize the velocities from the context.
        
        State s = context->getOwner().getState(State::Velocities);
        for (int i = 0; i < numCopies; i++)
            setVelocities(i, s.getVelocities());
    }
    if (isFirstStep) {
        // Call setPositions() on the Context so it doesn't think the user is trying to
        // run a simulation without setting positions first.  These positions will
        // immediately get overwritten by the ones stored in this integrator.
        
        vector<Vec3> p(context->getSystem().getNumParticles(), Vec3());
        context->getOwner().setPositions(p);
        isFirstStep = false;
    }
    for (int i = 0; i < steps; ++i) {
        kernel.getAs<IntegrateRPMDStepKernel>().execute(*context, *this, forcesAreValid);
        forcesAreValid = true;
    }
}
