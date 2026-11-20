/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit.                   *
 * See https://openmm.org/development.                                        *
 *                                                                            *
 * Portions copyright (c) 2015 Stanford University and the Authors.           *
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

#include "ReferenceTests.h"
#include "TestVariableLangevinIntegrator.h"

double getStepSizeWithMasslessParticles(int numMassless) {
    // Two massive particles joined by a stretched harmonic bond, plus
    // numMassless massless particles that experience a force but cannot move.
    // The massless particles contribute nothing to the integration error, so
    // they must not affect the step size selected by the integrator.
    System system;
    system.addParticle(2.0);
    system.addParticle(2.0);
    HarmonicBondForce* bond = new HarmonicBondForce();
    bond->addBond(0, 1, 1.5, 1);
    system.addForce(bond);
    CustomExternalForce* external = new CustomExternalForce("x");
    for (int i = 0; i < numMassless; ++i) {
        system.addParticle(0.0);
        external->addParticle(2+i);
    }
    system.addForce(external);
    VariableLangevinIntegrator integrator(0.0, 0.0, 1e-5);
    Context context(system, integrator, platform);
    vector<Vec3> positions(2+numMassless);
    positions[0] = Vec3(-1, 0, 0);
    positions[1] = Vec3(1, 0, 0);
    for (int i = 0; i < numMassless; ++i)
        positions[2+i] = Vec3(0, 0, 0);
    context.setPositions(positions);
    integrator.step(1);
    return integrator.getStepSize();
}

void testMasslessParticlesDoNotAffectStepSize() {
    // Massless particles (e.g. virtual sites) used to dilute the error
    // estimate, making the integrator select too large a step size
    // (see issue #5566).  The step size must be identical with and
    // without them.
    double stepWithout = getStepSizeWithMasslessParticles(0);
    double stepWith = getStepSizeWithMasslessParticles(20);
    ASSERT_EQUAL_TOL(stepWithout, stepWith, 1e-5);
}

void runPlatformTests() {
    testMasslessParticlesDoNotAffectStepSize();
}
