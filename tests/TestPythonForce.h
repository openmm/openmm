/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit.                   *
 * See https://openmm.org/development.                                        *
 *                                                                            *
 * Portions copyright (c) 2025 Stanford University and the Authors.           *
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

#include "openmm/internal/AssertionUtilities.h"
#include "openmm/Context.h"
#include "openmm/PythonForce.h"
#include "openmm/Platform.h"
#include "openmm/VerletIntegrator.h"
#include "sfmt/SFMT.h"
#include <iostream>

using namespace OpenMM;
using namespace std;

void testForce() {
    class Computation : public PythonForceComputation {
        void compute(const State& state, double& energy, void* forces, bool forcesAreDouble) const {
            ASSERT_EQUAL(5.0, state.getParameters().at("a"));
            ASSERT_EQUAL(10.0, state.getParameters().at("b"));
            Vec3 a, b, c;
            state.getPeriodicBoxVectors(a, b, c);
            ASSERT_EQUAL(Vec3(2, 0, 0), a);
            ASSERT_EQUAL(Vec3(0.1, 2, 0), b);
            ASSERT_EQUAL(Vec3(0.1, 0.1, 2), c);
            energy = 25.0;
            int numParticles = state.getPositions().size();
            for (int i = 0; i < numParticles; i++) {
                Vec3 f = state.getPositions()[i]*2;
                if (forcesAreDouble)
                    ((Vec3*) forces)[i] = f;
                else {
                    ((float*) forces)[3*i] = (float) f[0];
                    ((float*) forces)[3*i+1] = (float) f[1];
                    ((float*) forces)[3*i+2] = (float) f[2];
                }
            }
        }
    };
    int numParticles = 5;
    System system;
    Vec3 a(2, 0, 0);
    Vec3 b(0.1, 2, 0);
    Vec3 c(0.1, 0.1, 2);
    system.setDefaultPeriodicBoxVectors(a, b, c);
    vector<Vec3> positions;
    OpenMM_SFMT::SFMT sfmt;
    init_gen_rand(0, sfmt);
    for (int i = 0; i < numParticles; i++) {
        system.addParticle(1.0);
        positions.push_back(Vec3(genrand_real2(sfmt), genrand_real2(sfmt), genrand_real2(sfmt)));
    }
    map<string, double> params;
    params["a"] = 5.0;
    params["b"] = 10.0;
    PythonForce* force = new PythonForce(new Computation(), params);
    ASSERT(!force->usesPeriodicBoundaryConditions());
    force->setUsesPeriodicBoundaryConditions(true);
    ASSERT(force->usesPeriodicBoundaryConditions());
    system.addForce(force);
    VerletIntegrator integrator(0.01);
    Context context(system, integrator, platform);
    context.setPositions(positions);
    State state = context.getState(State::Energy | State::Forces);
    ASSERT_EQUAL_TOL(25.0, state.getPotentialEnergy(), 1e-6);
    for (int i = 0; i < numParticles; i++)
        ASSERT_EQUAL_VEC(2*positions[i], state.getForces()[i], 1e-6)

    // Check that force groups are handled correctly.

    ASSERT_EQUAL_TOL(25.0, context.getState(State::Energy, false, 1).getPotentialEnergy(), 1e-6);
    ASSERT_EQUAL_TOL(0.0, context.getState(State::Energy, false, 2).getPotentialEnergy(), 1e-6);
}

void runPlatformTests();

int main(int argc, char* argv[]) {
    try {
        initializeTests(argc, argv);
        testForce();
        runPlatformTests();
    }
    catch(const exception& e) {
        cout << "exception: " << e.what() << endl;
        return 1;
    }
    cout << "Done" << endl;
    return 0;
}
