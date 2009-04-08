/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2008-2009 Stanford University and the Authors.      *
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

#include "System.h"


/**
 * This tests the Cuda implementation of BrownianIntegrator.
 */

#include "../../../tests/AssertionUtilities.h"
#include "OpenMMContext.h"
#include "CudaPlatform.h"
#include "HarmonicBondForce.h"
#include "NonbondedForce.h"
#include "System.h"
#include "BrownianIntegrator.h"
#include "../src/SimTKUtilities/SimTKOpenMMRealType.h"
#include "../src/sfmt/SFMT.h"
#include <iostream>
#include <vector>

using namespace OpenMM;
using namespace std;

const double TOL = 1e-5;

void testSingleBond() {
    CudaPlatform platform;
    System system(2, 0);
    system.setParticleMass(0, 2.0);
    system.setParticleMass(1, 2.0);
    double dt = 0.01;
    BrownianIntegrator integrator(0, 0.1, dt);
    HarmonicBondForce* forceField = new HarmonicBondForce(1);
    forceField->setBondParameters(0, 0, 1, 1.5, 1);
    system.addForce(forceField);
    OpenMMContext context(system, integrator, platform);
    vector<Vec3> positions(2);
    positions[0] = Vec3(-1, 0, 0);
    positions[1] = Vec3(1, 0, 0);
    context.setPositions(positions);
    
    // This is simply an overdamped harmonic oscillator, so compare it to the analytical solution.
    
    double rate = 2*1.0/0.1;
    for (int i = 0; i < 1000; ++i) {
        State state = context.getState(State::Positions | State::Velocities);
        double time = state.getTime();
        double expectedDist = 1.5+0.5*std::exp(-rate*time);
        ASSERT_EQUAL_VEC(Vec3(-0.5*expectedDist, 0, 0), state.getPositions()[0], 0.02);
        ASSERT_EQUAL_VEC(Vec3(0.5*expectedDist, 0, 0), state.getPositions()[1], 0.02);
        if (i > 0) {
            double expectedSpeed = -0.5*rate*std::exp(-rate*(time-0.5*dt));
            ASSERT_EQUAL_VEC(Vec3(-0.5*expectedSpeed, 0, 0), state.getVelocities()[0], 0.11);
            ASSERT_EQUAL_VEC(Vec3(0.5*expectedSpeed, 0, 0), state.getVelocities()[1], 0.11);
        }
        integrator.step(1);
    }
}

void testTemperature() {
    const int numParticles = 8;
    const int numBonds = numParticles-1;
    const double temp = 10.0;
    CudaPlatform platform;
    System system(numParticles, 0);
    BrownianIntegrator integrator(temp, 2.0, 0.01);
    HarmonicBondForce* forceField = new HarmonicBondForce(numBonds);
    for (int i = 0; i < numParticles; ++i)
        system.setParticleMass(i, 2.0);
    for (int i = 0; i < numBonds; ++i)
        forceField->setBondParameters(i, i, i+1, 1.0, 5.0);
    system.addForce(forceField);
    OpenMMContext context(system, integrator, platform);
    vector<Vec3> positions(numParticles);
    for (int i = 0; i < numParticles; ++i)
        positions[i] = Vec3(i, 0, 0);
    context.setPositions(positions);
    
    // Let it equilibrate.
    
    integrator.step(10000);
    
    // Now run it for a while and see if the temperature is correct.
    
    double pe = 0.0;
    const int steps = 50000;
    for (int i = 0; i < steps; ++i) {
        State state = context.getState(State::Energy);
        pe += state.getPotentialEnergy();
        integrator.step(1);
    }
    pe /= steps;
    double expected = 0.5*numBonds*BOLTZ*temp;
    ASSERT_EQUAL_TOL(expected, pe, 20*expected/std::sqrt((double) steps));
}

void testConstraints() {
    const int numParticles = 8;
    const int numConstraints = 5;
    const double temp = 20.0;
    CudaPlatform platform;
    System system(numParticles, numConstraints);
    BrownianIntegrator integrator(temp, 2.0, 0.001);
    integrator.setConstraintTolerance(1e-5);
    NonbondedForce* forceField = new NonbondedForce(numParticles, 0);
    for (int i = 0; i < numParticles; ++i) {
        system.setParticleMass(i, 10.0);
        forceField->setParticleParameters(i, (i%2 == 0 ? 0.2 : -0.2), 0.5, 5.0);
    }
    system.setConstraintParameters(0, 0, 1, 1.0);
    system.setConstraintParameters(1, 1, 2, 1.0);
    system.setConstraintParameters(2, 2, 3, 1.0);
    system.setConstraintParameters(3, 4, 5, 1.0);
    system.setConstraintParameters(4, 6, 7, 1.0);
    system.addForce(forceField);
    OpenMMContext context(system, integrator, platform);
    vector<Vec3> positions(numParticles);
    vector<Vec3> velocities(numParticles);
    init_gen_rand(0);
    for (int i = 0; i < numParticles; ++i) {
        positions[i] = Vec3(i, 0, 0);
        velocities[i] = Vec3(genrand_real2()-0.5, genrand_real2()-0.5, genrand_real2()-0.5);
    }
    context.setPositions(positions);
    context.setVelocities(velocities);
    
    // Simulate it and see whether the constraints remain satisfied.
    
    for (int i = 0; i < 1000; ++i) {
        State state = context.getState(State::Positions);
        for (int j = 0; j < numConstraints; ++j) {
            int particle1, particle2;
            double distance;
            system.getConstraintParameters(j, particle1, particle2, distance);
            Vec3 p1 = state.getPositions()[particle1];
            Vec3 p2 = state.getPositions()[particle2];
            double dist = std::sqrt((p1[0]-p2[0])*(p1[0]-p2[0])+(p1[1]-p2[1])*(p1[1]-p2[1])+(p1[2]-p2[2])*(p1[2]-p2[2]));
            ASSERT_EQUAL_TOL(distance, dist, 2e-5);
        }
        integrator.step(1);
    }
}

void testRandomSeed() {
    const int numParticles = 8;
    const double temp = 100.0;
    const double collisionFreq = 10.0;
    CudaPlatform platform;
    System system(numParticles, 0);
    BrownianIntegrator integrator(temp, 2.0, 0.001);
    NonbondedForce* forceField = new NonbondedForce(numParticles, 0);
    for (int i = 0; i < numParticles; ++i) {
        system.setParticleMass(i, 2.0);
        forceField->setParticleParameters(i, (i%2 == 0 ? 1.0 : -1.0), 1.0, 5.0);
    }
    system.addForce(forceField);
    vector<Vec3> positions(numParticles);
    vector<Vec3> velocities(numParticles);
    for (int i = 0; i < numParticles; ++i) {
        positions[i] = Vec3((i%2 == 0 ? 2 : -2), (i%4 < 2 ? 2 : -2), (i < 4 ? 2 : -2));
        velocities[i] = Vec3(0, 0, 0);
    }

    // Try twice with the same random seed.

    integrator.setRandomNumberSeed(5);
    OpenMMContext context(system, integrator, platform);
    context.setPositions(positions);
    context.setVelocities(velocities);
    integrator.step(10);
    State state1 = context.getState(State::Positions);
    context.reinitialize();
    context.setPositions(positions);
    context.setVelocities(velocities);
    integrator.step(10);
    State state2 = context.getState(State::Positions);

    // Try twice with a different random seed.

    integrator.setRandomNumberSeed(10);
    context.reinitialize();
    context.setPositions(positions);
    context.setVelocities(velocities);
    integrator.step(10);
    State state3 = context.getState(State::Positions);
    context.reinitialize();
    context.setPositions(positions);
    context.setVelocities(velocities);
    integrator.step(10);
    State state4 = context.getState(State::Positions);

    // Compare the results.

    for (int i = 0; i < numParticles; i++) {
        for (int j = 0; j < 3; j++) {
            ASSERT(state1.getPositions()[i][j] == state2.getPositions()[i][j]);
            ASSERT(state3.getPositions()[i][j] == state4.getPositions()[i][j]);
            ASSERT(state1.getPositions()[i][j] != state3.getPositions()[i][j]);
        }
    }
}

int main() {
    try {
        testSingleBond();
        testTemperature();
        testConstraints();
        testRandomSeed();
    }
    catch(const exception& e) {
        cout << "exception: " << e.what() << endl;
        return 1;
    }
    cout << "Done" << endl;
    return 0;
}
