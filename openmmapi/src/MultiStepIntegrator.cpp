/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2008-2018 Stanford University and the Authors.      *
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

#include "openmm/kernels.h"
#include "openmm/MultiStepIntegrator.h"
#include "openmm/OpenMMException.h"
#include "openmm/internal/ContextImpl.h"

using namespace OpenMM;

MultiStep::MultiStep(int fastGroups_, int slowGroups_, int slowPeriod_)
    : fastGroups(fastGroups_), slowGroups(slowGroups_), slowPeriod(slowPeriod_), slowStep(0) {};

int MultiStep::getGroups() {
    ++slowStep;
    slowStep %= slowPeriod;
    return slowStep == 0 ? slowGroups | fastGroups : fastGroups;
}

MultiStepVerletIntegrator::MultiStepVerletIntegrator(double stepSize, int fastGroups, int slowGroups, int slowPeriod)
    : VerletIntegrator(stepSize), MultiStep(fastGroups, slowGroups, slowPeriod) {}

void MultiStepVerletIntegrator::step(int steps) {

    if (context == NULL)
        throw OpenMMException("This Integrator is not bound to a context!");

    for (int i = 0; i < steps; ++i) {
        context->updateContextState();
        context->calcForcesAndEnergy(true, false, getGroups());
        kernel.getAs<IntegrateVerletStepKernel>().execute(*context, *this);
    }
}

MultiStepLangevinIntegrator::MultiStepLangevinIntegrator(double temperature, double frictionCoeff, double stepSize,
                                                         int fastGroups, int slowGroups, int slowPeriod)
    : LangevinIntegrator(temperature, frictionCoeff, stepSize), MultiStep(fastGroups, slowGroups, slowPeriod) {}

void MultiStepLangevinIntegrator::step(int steps) {

    if (context == NULL)
        throw OpenMMException("This Integrator is not bound to a context!");

    for (int i = 0; i < steps; ++i) {
        context->updateContextState();
        context->calcForcesAndEnergy(true, false, getGroups());
        kernel.getAs<IntegrateLangevinStepKernel>().execute(*context, *this);
    }
}

MultiStepVerletIntegrator3::MultiStepVerletIntegrator3(double stepSize, int fastGroups, int slowGroups)
    : VerletIntegrator(stepSize/2), fastGroups(fastGroups), slowGroups(slowGroups) {}

void MultiStepVerletIntegrator3::step(int steps) {

    if (context == NULL)
        throw OpenMMException("This Integrator is not bound to a context!");

    for (int i = 0; i < steps; ++i) {
        context->updateContextState();
        context->calcForcesAndEnergy(true, false, fastGroups);
        kernel.getAs<IntegrateVerletStepKernel>().execute(*context, *this);

        context->updateContextState();
        // The forces of "slowGroups" have to be scaled by 2!
        context->calcForcesAndEnergy(true, false, fastGroups | slowGroups);
        kernel.getAs<IntegrateVerletStepKernel>().execute(*context, *this);
    }
}