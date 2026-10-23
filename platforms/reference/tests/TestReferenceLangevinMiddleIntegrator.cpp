/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit.                   *
 * See https://openmm.org/development.                                        *
 *                                                                            *
 * Portions copyright (c) 2019 Stanford University and the Authors.           *
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
#include "SimTKOpenMMUtilities.h"
#include "TestLangevinMiddleIntegrator.h"

void testRandomForceAxisOrder() {
    const vector<double> masses = {1.0, 0.0, 2.0, 3.0};
    const double temperature = 300.0;
    const double friction = 10.0;
    const double stepSize = 0.005;
    const double velocityScale = exp(-stepSize*friction);
    for (int seed : {1, 42, 12345}) {
        System system;
        vector<Vec3> expectedVelocities(masses.size());
        SimTKOpenMMUtilities::setRandomNumberSeed(seed);
        for (int i = 0; i < masses.size(); i++) {
            system.addParticle(masses[i]);
            if (masses[i] == 0.0)
                continue;
            const double noiseScale = sqrt(1-velocityScale*velocityScale)*sqrt(BOLTZ*temperature/masses[i]);
            // The Reference trajectory assigns successive random draws to z, y, then x.
            for (int axis = 2; axis >= 0; axis--)
                expectedVelocities[i][axis] = noiseScale*SimTKOpenMMUtilities::getNormallyDistributedRandomNumber();
        }

        LangevinMiddleIntegrator integrator(temperature, friction, stepSize);
        integrator.setRandomNumberSeed(seed);
        Context context(system, integrator, platform);
        context.setPositions(vector<Vec3>(masses.size()));
        context.setVelocities(vector<Vec3>(masses.size()));
        integrator.step(1);
        State state = context.getState(State::Positions | State::Velocities);
        for (int i = 0; i < masses.size(); i++) {
            ASSERT_EQUAL_VEC(expectedVelocities[i], state.getVelocities()[i], 1e-12);
            ASSERT_EQUAL_VEC(expectedVelocities[i]*(0.5*stepSize), state.getPositions()[i], 1e-12);
        }
    }
}

void runPlatformTests() {
    testRandomForceAxisOrder();
}
