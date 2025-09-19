/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit.                   *
 * See https://openmm.org/development.                                        *
 *                                                                            *
 * Portions copyright (c) 2008-2025 Stanford University and the Authors.      *
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

#include "openmm/DrudeSCFIntegrator.h"
#include "openmm/DrudeKernels.h"
#include "openmm/Context.h"
#include "openmm/OpenMMException.h"
#include "openmm/internal/ContextImpl.h"
#include "openmm/internal/DrudeHelpers.h"
#include <cmath>
#include <ctime>
#include <string>

using namespace OpenMM;
using std::string;
using std::vector;

DrudeSCFIntegrator::DrudeSCFIntegrator(double stepSize) : DrudeIntegrator(stepSize)
{
    setDrudeTemperature(0.0);  // This is only used to initialize velocities for this integrator
    setStepSize(stepSize);
    setMinimizationErrorTolerance(1.0);
    setConstraintTolerance(1e-5);
    setMaxDrudeDistance(0.0);
}

void DrudeSCFIntegrator::initialize(ContextImpl& contextRef) {
    if (owner != NULL && &contextRef.getOwner() != owner)
        throw OpenMMException("This Integrator is already bound to a context");
    const DrudeForce* force = getDrudeForce(contextRef);
    if (getMaxDrudeDistance() != 0.0)
        throw OpenMMException("DrudeSCFIntegrator does not currently support setting max Drude distance");
    context = &contextRef;
    owner = &contextRef.getOwner();
    kernel = context->getPlatform().createKernel(IntegrateDrudeSCFStepKernel::Name(), contextRef);
    kernel.getAs<IntegrateDrudeSCFStepKernel>().initialize(contextRef.getSystem(), *this, *force);
}

void DrudeSCFIntegrator::setMinimizationErrorTolerance(double tol) {
    if (tol <= 0)
        throw OpenMMException("Minimization error tolerance must be positive");
    tolerance = tol;
}

void DrudeSCFIntegrator::cleanup() {
    kernel = Kernel();
}

vector<string> DrudeSCFIntegrator::getKernelNames() {
    std::vector<std::string> names;
    names.push_back(IntegrateDrudeSCFStepKernel::Name());
    return names;
}

double DrudeSCFIntegrator::computeKineticEnergy() {
    return kernel.getAs<IntegrateDrudeSCFStepKernel>().computeKineticEnergy(*context, *this);
}

void DrudeSCFIntegrator::step(int steps) {
    if (context == NULL)
        throw OpenMMException("This Integrator is not bound to a context!");    
    for (int i = 0; i < steps; ++i) {
        context->updateContextState();
        context->calcForcesAndEnergy(true, false);
        kernel.getAs<IntegrateDrudeSCFStepKernel>().execute(*context, *this);
    }
}
