/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2008-2015 Stanford University and the Authors.      *
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

#include "openmm/DrudeLangevinIntegrator.h"
#include "openmm/Context.h"
#include "openmm/OpenMMException.h"
#include "openmm/internal/ContextImpl.h"
#include "openmm/DrudeKernels.h"
#include <ctime>
#include <string>

using namespace OpenMM;
using std::string;
using std::vector;

DrudeLangevinIntegrator::DrudeLangevinIntegrator(double temperature, double frictionCoeff, double drudeTemperature, double drudeFrictionCoeff, double stepSize) {
    setTemperature(temperature);
    setFriction(frictionCoeff);
    setDrudeTemperature(drudeTemperature);
    setDrudeFriction(drudeFrictionCoeff);
    setMaxDrudeDistance(0);
    setStepSize(stepSize);
    setConstraintTolerance(1e-5);
    setRandomNumberSeed(0);
}

double DrudeLangevinIntegrator::getMaxDrudeDistance() const {
    return maxDrudeDistance;
}

void DrudeLangevinIntegrator::setMaxDrudeDistance(double distance) {
    if (distance < 0)
        throw OpenMMException("setMaxDrudeDistance: Distance cannot be negative");
    maxDrudeDistance = distance;
}

void DrudeLangevinIntegrator::initialize(ContextImpl& contextRef) {
    if (owner != NULL && &contextRef.getOwner() != owner)
        throw OpenMMException("This Integrator is already bound to a context");
    const DrudeForce* force = NULL;
    const System& system = contextRef.getSystem();
    for (int i = 0; i < system.getNumForces(); i++)
        if (dynamic_cast<const DrudeForce*>(&system.getForce(i)) != NULL) {
            if (force == NULL)
                force = dynamic_cast<const DrudeForce*>(&system.getForce(i));
            else
                throw OpenMMException("The System contains multiple DrudeForces");
        }
    if (force == NULL)
        throw OpenMMException("The System does not contain a DrudeForce");
    context = &contextRef;
    owner = &contextRef.getOwner();
    kernel = context->getPlatform().createKernel(IntegrateDrudeLangevinStepKernel::Name(), contextRef);
    kernel.getAs<IntegrateDrudeLangevinStepKernel>().initialize(contextRef.getSystem(), *this, *force);
}

void DrudeLangevinIntegrator::cleanup() {
    kernel = Kernel();
}

void DrudeLangevinIntegrator::stateChanged(State::DataType changed) {
}

vector<string> DrudeLangevinIntegrator::getKernelNames() {
    std::vector<std::string> names;
    names.push_back(IntegrateDrudeLangevinStepKernel::Name());
    return names;
}

double DrudeLangevinIntegrator::computeKineticEnergy() {
    return kernel.getAs<IntegrateDrudeLangevinStepKernel>().computeKineticEnergy(*context, *this);
}

void DrudeLangevinIntegrator::step(int steps) {
    if (context == NULL)
        throw OpenMMException("This Integrator is not bound to a context!");    
    for (int i = 0; i < steps; ++i) {
        context->updateContextState();
        context->calcForcesAndEnergy(true, false);
        kernel.getAs<IntegrateDrudeLangevinStepKernel>().execute(*context, *this);
    }
}
