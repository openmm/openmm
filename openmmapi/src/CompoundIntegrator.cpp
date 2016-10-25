/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2015-2016 Stanford University and the Authors.      *
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

#include "openmm/CompoundIntegrator.h"
#include "openmm/Context.h"
#include "openmm/OpenMMException.h"
#include "openmm/internal/ContextImpl.h"
#include "openmm/internal/AssertionUtilities.h"
#include "openmm/kernels.h"
#include <string>

using namespace OpenMM;
using namespace std;

CompoundIntegrator::CompoundIntegrator() : currentIntegrator(0) {
}

CompoundIntegrator::~CompoundIntegrator() {
    for (int i = 0; i < integrators.size(); i++)
        delete integrators[i];
}

int CompoundIntegrator::getNumIntegrators() const {
    return integrators.size();
}

int CompoundIntegrator::addIntegrator(Integrator* integrator) {
    if (owner != NULL)
        throw OpenMMException("addIntegrator() cannot be called after a CustomIntegrator is already bound to a context");
    integrators.push_back(integrator);
    return integrators.size()-1;
}

Integrator& CompoundIntegrator::getIntegrator(int index) {
    ASSERT_VALID_INDEX(index, integrators);
    return *integrators[index];
}

const Integrator& CompoundIntegrator::getIntegrator(int index) const {
    ASSERT_VALID_INDEX(index, integrators);
    return *integrators[index];
}

int CompoundIntegrator::getCurrentIntegrator() const {
    return currentIntegrator;
}

void CompoundIntegrator::setCurrentIntegrator(int index) {
    if (index < 0 || index >= integrators.size())
        throw OpenMMException("Illegal index for setCurrentIntegrator()");
    currentIntegrator = index;
}

double CompoundIntegrator::getStepSize() const {
    return integrators[currentIntegrator]->getStepSize();
}

void CompoundIntegrator::setStepSize(double size) {
    integrators[currentIntegrator]->setStepSize(size);
}

double CompoundIntegrator::getConstraintTolerance() const {
    return integrators[currentIntegrator]->getConstraintTolerance();
}

void CompoundIntegrator::setConstraintTolerance(double tol) {
    integrators[currentIntegrator]->setConstraintTolerance(tol);
}

void CompoundIntegrator::step(int steps) {
    integrators[currentIntegrator]->step(steps);
}

void CompoundIntegrator::initialize(ContextImpl& context) {
    if (integrators.size() == 0)
        throw OpenMMException("CompoundIntegrator must contain at least one Integrator");
    for (int i = 0; i < integrators.size(); i++)
        integrators[i]->initialize(context);
}

void CompoundIntegrator::cleanup() {
    for (int i = 0; i < integrators.size(); i++)
        integrators[i]->cleanup();
}

vector<string> CompoundIntegrator::getKernelNames() {
    vector<string> kernels;
    for (int i = 0; i < integrators.size(); i++) {
        vector<string> integratorKernels = integrators[i]->getKernelNames();
        kernels.insert(kernels.end(), integratorKernels.begin(), integratorKernels.end());
    }
    return kernels;
}

void CompoundIntegrator::stateChanged(State::DataType changed) {
    for (int i = 0; i < integrators.size(); i++)
        integrators[i]->stateChanged(changed);
}

double CompoundIntegrator::computeKineticEnergy() {
    return integrators[currentIntegrator]->computeKineticEnergy();
}
