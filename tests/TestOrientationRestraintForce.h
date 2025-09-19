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
#include "openmm/NonbondedForce.h"
#include "openmm/OrientationRestraintForce.h"
#include "openmm/System.h"
#include "openmm/VerletIntegrator.h"
#include "sfmt/SFMT.h"
#include <cmath>
#include <iostream>
#include <vector>

using namespace OpenMM;
using namespace std;

void testOrientationRestraint() {
    const int numParticles = 20;
    const double k = 3.5;
    System system;
    vector<Vec3> referencePos(numParticles);
    vector<Vec3> positions(numParticles);
    vector<int> particles;
    OpenMM_SFMT::SFMT sfmt;
    init_gen_rand(0, sfmt);
    Vec3 center;
    for (int i = 0; i < numParticles; ++i) {
        system.addParticle(1.0);
        referencePos[i] = Vec3(genrand_real2(sfmt), genrand_real2(sfmt), genrand_real2(sfmt))*10;
        if (i%5 != 0) {
            particles.push_back(i);
            center += referencePos[i];
        }
    }
    center /= particles.size();
    OrientationRestraintForce* force = new OrientationRestraintForce(k, referencePos, particles);
    system.addForce(force);
    VerletIntegrator integrator(0.001);
    Context context(system, integrator, platform);

    // Randomly transform the reference positions and see if the energy is correct.

    for (int i = 0; i < 20; i++) {
        // Select a random axis, angle, and translation, and transform the particles from the reference positions.

        Vec3 translation(genrand_real2(sfmt), genrand_real2(sfmt), genrand_real2(sfmt));
        Vec3 axis(genrand_real2(sfmt), genrand_real2(sfmt), genrand_real2(sfmt));
        axis /= sqrt(axis.dot(axis));
        double angle = 3*genrand_real2(sfmt);
        for (int j : particles) {
            Vec3 p = referencePos[j]-center;
            Vec3 cross1 = axis.cross(p);
            p += sin(angle)*cross1 + (1-cos(angle))*axis.cross(cross1);
            p += translation;
            positions[j] = p;
        }
        context.setPositions(positions);
        State state = context.getState(State::Energy);
        double s = sin(angle/2);
        ASSERT_EQUAL_TOL(2*k*s*s, state.getPotentialEnergy(), 1e-6);
    }

    // Take a small step in the direction of the energy gradient and see whether the potential energy changes by the expected amount.

    for (int i = 0; i < 20; i++) {
        for (int j = 0; j < numParticles; ++j)
            positions[j] = Vec3(genrand_real2(sfmt), genrand_real2(sfmt), genrand_real2(sfmt))*10;
        context.setPositions(positions);
        vector<Vec3> forces = context.getState(State::Forces).getForces();
        double norm = 0.0;
        for (int i = 0; i < forces.size(); ++i)
            norm += forces[i].dot(forces[i]);
        norm = sqrt(norm);
        const double stepSize = 0.01;
        double step = 0.5*stepSize/norm;
        vector<Vec3> positions2(numParticles), positions3(numParticles);
        for (int i = 0; i < positions.size(); ++i) {
            Vec3 p = positions[i];
            Vec3 f = forces[i];
            positions2[i] = Vec3(p[0]-f[0]*step, p[1]-f[1]*step, p[2]-f[2]*step);
            positions3[i] = Vec3(p[0]+f[0]*step, p[1]+f[1]*step, p[2]+f[2]*step);
        }
        context.setPositions(positions2);
        State state2 = context.getState(State::Energy);
        context.setPositions(positions3);
        State state3 = context.getState(State::Energy);
        ASSERT_EQUAL_TOL(norm, (state2.getPotentialEnergy()-state3.getPotentialEnergy())/stepSize, 1e-3);
    }

    // When the current positions equal the reference positions, all forces should be zero.

    context.setPositions(referencePos);
    vector<Vec3> forces = context.getState(State::Forces).getForces();
    Vec3 zero;
    for (int i = 0; i < numParticles; i++)
        ASSERT_EQUAL_VEC(zero, forces[i], 1e-4);
    
    // Check that updateParametersInContext() works correctly.
    
    context.setPositions(positions);
    double e1 = context.getState(State::Energy).getPotentialEnergy();
    force->setK(2*k);
    force->updateParametersInContext(context);
    double e2 = context.getState(State::Energy).getPotentialEnergy();
    force->setReferencePositions(positions);
    force->updateParametersInContext(context);
    double e3 = context.getState(State::Energy).getPotentialEnergy();
    ASSERT_EQUAL_TOL(2*e1, e2, 1e-6);
    ASSERT_EQUAL_TOL(0.0, e3, 1e-6);
}

void testEnergyConservation() {
    const int numParticles = 50;
    System system;
    vector<Vec3> referencePos(numParticles);
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
        referencePos[i] = Vec3(genrand_real2(sfmt), genrand_real2(sfmt), genrand_real2(sfmt))*5;
        if (genrand_real2(sfmt) < 0.5)
            particles.push_back(i);
    }
    OrientationRestraintForce* force = new OrientationRestraintForce(10.0, referencePos, particles);
    system.addForce(force);
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

    // If we modify the reference positions, the energy should change.

    for (int i = 0; i < numParticles; ++i)
        referencePos[i] = Vec3(genrand_real2(sfmt), genrand_real2(sfmt), genrand_real2(sfmt))*5;
    force->setReferencePositions(referencePos);
    force->updateParametersInContext(context);
    State state2 = context.getState(State::Energy);
    double energy2 = state2.getPotentialEnergy()+state2.getKineticEnergy();
    ASSERT(fabs(energy-energy2) > 1e-3);

    // Make sure it's still conserved.

    for (int i = 0; i < 100; i++) {
        integrator.step(5);
        State state = context.getState(State::Energy);
        ASSERT_EQUAL_TOL(energy2, state.getPotentialEnergy()+state.getKineticEnergy(), 1e-4);
    }
}

void runPlatformTests();

int main(int argc, char* argv[]) {
    try {
        initializeTests(argc, argv);
        testOrientationRestraint();
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
