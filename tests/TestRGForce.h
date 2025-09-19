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
#include "openmm/RGForce.h"
#include "openmm/Context.h"
#include "openmm/NonbondedForce.h"
#include "openmm/System.h"
#include "openmm/VerletIntegrator.h"
#include "sfmt/SFMT.h"
#include <cmath>
#include <iostream>
#include <vector>

using namespace OpenMM;
using namespace std;

void testRG(bool allParticles) {
    const int numParticles = 30;
    System system;
    vector<Vec3> positions(numParticles);
    vector<int> particles;
    OpenMM_SFMT::SFMT sfmt;
    init_gen_rand(0, sfmt);
    for (int i = 0; i < numParticles; ++i) {
        system.addParticle(1.0);
        if (genrand_real2(sfmt) < 0.5 || allParticles)
            particles.push_back(i);
    }
    if (allParticles)
        system.addForce(new RGForce()); // Omitting the list of particles should mean all particles.
    else
        system.addForce(new RGForce(particles));
    VerletIntegrator integrator(0.001);
    Context context(system, integrator, platform);
    for (int i = 0; i < 10; i++) {
        // Set all particles to random positions.

        for (int j = 0; j < numParticles; j++)
            positions[j] = Vec3(genrand_real2(sfmt), genrand_real2(sfmt), genrand_real2(sfmt))*5;
        context.setPositions(positions);

        // Compute Rg.

        Vec3 center;
        for (int j : particles)
            center += positions[j];
        center /= particles.size();
        double sum = 0;
        for (int j : particles) {
            Vec3 v = positions[j]-center;
            sum += v.dot(v);
        }
        double rg = sqrt(sum/particles.size());

        // Compare to the value computed by the force.

        State state = context.getState(State::Energy | State::Forces);
        ASSERT_EQUAL_TOL(rg, state.getPotentialEnergy(), 1e-6);

        // Take a small step in the direction of the energy gradient and see whether the potential energy changes by the expected amount.

        const vector<Vec3>& forces = state.getForces();
        double norm = 0.0;
        for (int j = 0; j < (int) forces.size(); ++j)
            norm += forces[j].dot(forces[j]);
        norm = sqrt(norm);
        const double stepSize = 0.1;
        double step = 0.5*stepSize/norm;
        vector<Vec3> positions2(numParticles), positions3(numParticles);
        for (int j = 0; j < positions.size(); ++j) {
            Vec3 p = positions[j];
            Vec3 f = forces[j];
            positions2[j] = Vec3(p[0]-f[0]*step, p[1]-f[1]*step, p[2]-f[2]*step);
            positions3[j] = Vec3(p[0]+f[0]*step, p[1]+f[1]*step, p[2]+f[2]*step);
        }
        context.setPositions(positions2);
        State state2 = context.getState(State::Energy);
        context.setPositions(positions3);
        State state3 = context.getState(State::Energy);
        ASSERT_EQUAL_TOL(norm, (state2.getPotentialEnergy()-state3.getPotentialEnergy())/stepSize, 1e-4);
    }
}

void testEnergyConservation() {
    const int numParticles = 50;
    System system;
    vector<Vec3> positions(numParticles);
    vector<int> particles;
    OpenMM_SFMT::SFMT sfmt;
    init_gen_rand(0, sfmt);
    NonbondedForce* nb = new NonbondedForce(); // Add a nonbonded force to activate reordering on the GPU
    nb->setNonbondedMethod(NonbondedForce::CutoffNonPeriodic);
    system.addForce(nb);
    for (int i = 0; i < numParticles; ++i) {
        system.addParticle(2.0);
        nb->addParticle(0.0, 0.1, 0.01);
        positions[i] = Vec3(genrand_real2(sfmt), genrand_real2(sfmt), genrand_real2(sfmt))*5;
        if (genrand_real2(sfmt) < 0.5)
            particles.push_back(i);
    }
    system.addForce(new RGForce(particles));
    VerletIntegrator integrator(0.001);
    Context context(system, integrator, platform);
    context.setPositions(positions);
    context.setVelocitiesToTemperature(300.0, 0);
    integrator.step(5);
    State initialState = context.getState(State::Energy);
    double energy = initialState.getPotentialEnergy()+initialState.getKineticEnergy();
    for (int i = 0; i < 100; i++) {
        integrator.step(5);
        State state = context.getState(State::Energy);
        ASSERT_EQUAL_TOL(energy, state.getPotentialEnergy()+state.getKineticEnergy(), 1e-4);
    }
}

void runPlatformTests();

int main(int argc, char* argv[]) {
    try {
        initializeTests(argc, argv);
        testRG(true);
        testRG(false);
        testEnergyConservation();
        runPlatformTests();
    }
    catch(const exception& e) {
        cout << "exception: " << e.what() << endl;
        return 1;
    }
    cout << "Done" << endl;
    return 0;
}
