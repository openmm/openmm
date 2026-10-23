/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit.                   *
 * See https://openmm.org/development.                                        *
 *                                                                            *
 * Portions copyright (c) 2013 Stanford University and the Authors.           *
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

#include "ReferenceDrudeTests.h"
#include "TestDrudeLangevinIntegrator.h"
#include "openmm/CustomExternalForce.h"

void runPlatformTests() {
    // At zero temperature, the first step under constant forces has an analytic
    // solution for both the pair's center of mass and its relative motion.
    const double dt = 0.001;
    for (double friction : {0.0, 1e-16, 1.0}) {
        for (double drudeFriction : {0.0, 1e-16, 1.0}) {
            System system;
            for (double mass : {1.0, 2.0, 3.0})
                system.addParticle(mass);
            DrudeForce* drude = new DrudeForce();
            drude->addParticle(1, 0, -1, -1, -1, 0.1, 0.001, 1, 1);
            system.addForce(drude);
            CustomExternalForce* force = new CustomExternalForce("-f*x");
            force->addPerParticleParameter("f");
            force->addParticle(1, {2.0});
            force->addParticle(2, {3.0});
            system.addForce(force);
            DrudeLangevinIntegrator integrator(0, friction, 0, drudeFriction, dt);
            integrator.setMaxDrudeDistance(0);
            Context context(system, integrator, platform);
            context.setPositions(vector<Vec3>(3));
            context.setVelocities(vector<Vec3>(3));
            integrator.step(1);
            State state = context.getState(State::Positions | State::Velocities);
            // For 1e-16 friction, the correction to dt is below double precision.
            const double a = friction == 1 ? -expm1(-dt) : dt;
            const double b = drudeFriction == 1 ? -expm1(-dt) : dt;
            const double expected[] = {2*(a-b)/3, (2*a+b)/3, a};
            for (int i = 0; i < 3; i++) {
                for (int j = 0; j < 3; j++) {
                    ASSERT(std::isfinite(state.getPositions()[i][j]));
                    ASSERT(std::isfinite(state.getVelocities()[i][j]));
                }
                ASSERT_EQUAL_VEC(Vec3(expected[i], 0, 0), state.getVelocities()[i], 1e-12);
                ASSERT_EQUAL_VEC(Vec3(dt*expected[i], 0, 0), state.getPositions()[i], 1e-14);
            }
        }
    }
}
