
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

#include "openmm/internal/AssertionUtilities.h"
#include "openmm/Context.h"
#include "openmm/NonbondedForce.h"
#include "openmm/System.h"
#include "openmm/LangevinIntegrator.h"
#include "sfmt/SFMT.h"
#include <iostream>
#include <vector>

using namespace OpenMM;
using namespace std;

void testConstraints() {
    const int numMolecules = 10;
    const int numParticles = numMolecules*3;
    const int numConstraints = numMolecules*3;
    const double temp = 100.0;
    System system;
    LangevinIntegrator integrator(temp, 2.0, 0.001);
    integrator.setConstraintTolerance(1e-5);
    NonbondedForce* forceField = new NonbondedForce();
    for (int i = 0; i < numMolecules; ++i) {
        system.addParticle(16.0);
        system.addParticle(1.0);
        system.addParticle(1.0);
        forceField->addParticle(-0.82, 0.317, 0.65);
        forceField->addParticle(0.41, 1.0, 0.0);
        forceField->addParticle(0.41, 1.0, 0.0);
        system.addConstraint(i*3, i*3+1, 0.1);
        system.addConstraint(i*3, i*3+2, 0.1);
        system.addConstraint(i*3+1, i*3+2, 0.163);
    }
    system.addForce(forceField);
    Context context(system, integrator, platform);
    vector<Vec3> positions(numParticles);
    vector<Vec3> velocities(numParticles);
    OpenMM_SFMT::SFMT sfmt;
    init_gen_rand(0, sfmt);

    for (int i = 0; i < numMolecules; ++i) {
        positions[i*3] = Vec3((i%4)*0.4, (i/4)*0.4, 0);
        positions[i*3+1] = positions[i*3]+Vec3(0.1, 0, 0);
        positions[i*3+2] = positions[i*3]+Vec3(-0.03333, 0.09428, 0);
        velocities[i*3] = Vec3(genrand_real2(sfmt)-0.5, genrand_real2(sfmt)-0.5, genrand_real2(sfmt)-0.5);
        velocities[i*3+1] = Vec3(genrand_real2(sfmt)-0.5, genrand_real2(sfmt)-0.5, genrand_real2(sfmt)-0.5);
        velocities[i*3+2] = Vec3(genrand_real2(sfmt)-0.5, genrand_real2(sfmt)-0.5, genrand_real2(sfmt)-0.5);
    }
    context.setPositions(positions);
    context.setVelocities(velocities);

    // Simulate it and see whether the constraints remain satisfied.

    for (int i = 0; i < 1000; ++i) {
        integrator.step(1);
        State state = context.getState(State::Positions | State::Forces);
        for (int j = 0; j < numConstraints; ++j) {
            int particle1, particle2;
            double distance;
            system.getConstraintParameters(j, particle1, particle2, distance);
            Vec3 p1 = state.getPositions()[particle1];
            Vec3 p2 = state.getPositions()[particle2];
            double dist = std::sqrt((p1[0]-p2[0])*(p1[0]-p2[0])+(p1[1]-p2[1])*(p1[1]-p2[1])+(p1[2]-p2[2])*(p1[2]-p2[2]));
            ASSERT_EQUAL_TOL(distance, dist, 1e-5);
        }
    }
}

void runPlatformTests();

int main(int argc, char* argv[]) {
    try {
        initializeTests(argc, argv);
        testConstraints();
        runPlatformTests();
    }
    catch(const exception& e) {
        cout << "exception: " << e.what() << endl;
        return 1;
    }
    cout << "Done" << endl;
    return 0;
}
