/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit.                   *
 * See https://openmm.org/development.                                        *
 *                                                                            *
 * Portions copyright (c) 2015-2025 Stanford University and the Authors.      *
 * Authors: Evan Pretti                                                       *
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
#include "TestConstantPotentialForce.h"

void platformInitialize() {
}

void testGradientFiniteDifference(ConstantPotentialForce::ConstantPotentialMethod method, bool usePreconditioner) {
    // Ensures that computed forces match actual changes in energy with particle
    // perturbations, accounting for changes in electrode atom charges.

    // Finite differences with single precision PME have a lot of error, so we
    // only run this test on the reference platform, and the other platforms'
    // forces are compared to forces computed by the reference platform.

    System system;
    ConstantPotentialForce* force;
    vector<Vec3> positions;
    makeTestUpdateSystem(method, usePreconditioner, system, force, positions);
    force->setEwaldErrorTolerance(2e-6);
    VerletIntegrator integrator(0.001);
    Context context(system, integrator, platform);
    context.setPositions(positions);

    State state = context.getState(State::Energy | State::Forces);
    double energy = state.getPotentialEnergy();
    vector<Vec3> forces = state.getForces();

    double delta = 1e-3;
    for (int i = 0; i < forces.size(); i++) {
        for (int d = 0; d < 3; d++) {
            double refPos = positions[i][d];

            positions[i][d] = refPos - delta;
            context.setPositions(positions);
            vector<double> c;
            double energyL = context.getState(State::Energy).getPotentialEnergy();
            positions[i][d] = refPos + delta;
            context.setPositions(positions);
            double energyR = context.getState(State::Energy).getPotentialEnergy();
            positions[i][d] = refPos;

            ASSERT_EQUAL_TOL((energyR - energyL) / (2 * delta), -forces[i][d], 5e-4);
        }
    }
}

void runPlatformTests(ConstantPotentialForce::ConstantPotentialMethod method, bool usePreconditioner) {
    testEnergyConservation(method, usePreconditioner, 10);
    testGradientFiniteDifference(method, usePreconditioner);
}
