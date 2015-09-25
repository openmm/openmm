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
#include "openmm/HarmonicBondForce.h"
#include "openmm/NonbondedForce.h"
#include "openmm/System.h"
#include "openmm/VariableVerletIntegrator.h"
#include "SimTKOpenMMRealType.h"
#include "sfmt/SFMT.h"
#include <iostream>
#include <vector>

using namespace OpenMM;
using namespace std;

const double TOL = 1e-5;

void testSingleBond() {
    System system;
    system.addParticle(2.0);
    system.addParticle(2.0);
    VariableVerletIntegrator integrator(1e-6);
    HarmonicBondForce* forceField = new HarmonicBondForce();
    forceField->addBond(0, 1, 1.5, 1);
    system.addForce(forceField);
    Context context(system, integrator, platform);
    vector<Vec3> positions(2);
    positions[0] = Vec3(-1, 0, 0);
    positions[1] = Vec3(1, 0, 0);
    context.setPositions(positions);

    // This is simply a harmonic oscillator, so compare it to the analytical solution.

    const double freq = 1.0;
    State state = context.getState(State::Energy);
    const double initialEnergy = state.getKineticEnergy()+state.getPotentialEnergy();
    for (int i = 0; i < 1000; ++i) {
        state = context.getState(State::Positions | State::Velocities | State::Energy);
        double time = state.getTime();
        double expectedDist = 1.5+0.5*std::cos(freq*time);
        ASSERT_EQUAL_VEC(Vec3(-0.5*expectedDist, 0, 0), state.getPositions()[0], 0.02);
        ASSERT_EQUAL_VEC(Vec3(0.5*expectedDist, 0, 0), state.getPositions()[1], 0.02);
        double expectedSpeed = -0.5*freq*std::sin(freq*time);
        ASSERT_EQUAL_VEC(Vec3(-0.5*expectedSpeed, 0, 0), state.getVelocities()[0], 0.02);
        ASSERT_EQUAL_VEC(Vec3(0.5*expectedSpeed, 0, 0), state.getVelocities()[1], 0.02);
        double energy = state.getKineticEnergy()+state.getPotentialEnergy();
        ASSERT_EQUAL_TOL(initialEnergy, energy, 0.05);
        integrator.step(1);
    }
    ASSERT(state.getTime() > 1.0);
}

void testConstraints() {
    const int numParticles = 8;
    const int numConstraints = 5;
    System system;
    VariableVerletIntegrator integrator(1e-5);
    integrator.setConstraintTolerance(1e-5);
    NonbondedForce* forceField = new NonbondedForce();
    for (int i = 0; i < numParticles; ++i) {
        system.addParticle(i%2 == 0 ? 5.0 : 10.0);
        forceField->addParticle((i%2 == 0 ? 0.2 : -0.2), 0.5, 5.0);
    }
    system.addConstraint(0, 1, 1.0);
    system.addConstraint(1, 2, 1.0);
    system.addConstraint(2, 3, 1.0);
    system.addConstraint(4, 5, 1.0);
    system.addConstraint(6, 7, 1.0);
    system.addForce(forceField);
    Context context(system, integrator, platform);
    vector<Vec3> positions(numParticles);
    vector<Vec3> velocities(numParticles);
    OpenMM_SFMT::SFMT sfmt;
    init_gen_rand(0, sfmt);

    for (int i = 0; i < numParticles; ++i) {
        positions[i] = Vec3(i/2, (i+1)/2, 0);
        velocities[i] = Vec3(genrand_real2(sfmt)-0.5, genrand_real2(sfmt)-0.5, genrand_real2(sfmt)-0.5);
    }
    context.setPositions(positions);
    context.setVelocities(velocities);

    // Simulate it and see whether the constraints remain satisfied.

    double initialEnergy = 0.0;
    for (int i = 0; i < 1000; ++i) {
        State state = context.getState(State::Positions | State::Energy | State::Velocities | State::Forces);
        for (int j = 0; j < numConstraints; ++j) {
            int particle1, particle2;
            double distance;
            system.getConstraintParameters(j, particle1, particle2, distance);
            Vec3 p1 = state.getPositions()[particle1];
            Vec3 p2 = state.getPositions()[particle2];
            double dist = std::sqrt((p1[0]-p2[0])*(p1[0]-p2[0])+(p1[1]-p2[1])*(p1[1]-p2[1])+(p1[2]-p2[2])*(p1[2]-p2[2]));
            ASSERT_EQUAL_TOL(distance, dist, 1e-4);
        }
        double energy = state.getKineticEnergy()+state.getPotentialEnergy();
        if (i == 1)
            initialEnergy = energy;
        else if (i > 1)
            ASSERT_EQUAL_TOL(initialEnergy, energy, 0.01);
        integrator.step(1);
    }
    double finalTime = context.getState(State::Positions).getTime();
    ASSERT(finalTime > 0.1);

    // Now try the stepTo() method.

    finalTime += 0.5;
    integrator.stepTo(finalTime);
    ASSERT_EQUAL(finalTime, context.getState(State::Positions).getTime());
}

void testConstrainedClusters() {
    const int numParticles = 7;
    System system;
    VariableVerletIntegrator integrator(1e-5);
    integrator.setConstraintTolerance(1e-5);
    NonbondedForce* forceField = new NonbondedForce();
    for (int i = 0; i < numParticles; ++i) {
        system.addParticle(i > 1 ? 1.0 : 10.0);
        forceField->addParticle((i%2 == 0 ? 0.2 : -0.2), 0.5, 5.0);
    }
    system.addConstraint(0, 1, 1.0);
    system.addConstraint(0, 2, 1.0);
    system.addConstraint(0, 3, 1.0);
    system.addConstraint(0, 4, 1.0);
    system.addConstraint(1, 5, 1.0);
    system.addConstraint(1, 6, 1.0);
    system.addConstraint(2, 3, sqrt(2.0));
    system.addConstraint(2, 4, sqrt(2.0));
    system.addConstraint(3, 4, sqrt(2.0));
    system.addConstraint(5, 6, sqrt(2.0));
    system.addForce(forceField);
    Context context(system, integrator, platform);
    vector<Vec3> positions(numParticles);
    positions[0] = Vec3(0, 0, 0);
    positions[1] = Vec3(1, 0, 0);
    positions[2] = Vec3(-1, 0, 0);
    positions[3] = Vec3(0, 1, 0);
    positions[4] = Vec3(0, 0, 1);
    positions[5] = Vec3(2, 0, 0);
    positions[6] = Vec3(1, 1, 0);
    vector<Vec3> velocities(numParticles);
    OpenMM_SFMT::SFMT sfmt;
    init_gen_rand(0, sfmt);

    for (int i = 0; i < numParticles; ++i)
        velocities[i] = Vec3(genrand_real2(sfmt)-0.5, genrand_real2(sfmt)-0.5, genrand_real2(sfmt)-0.5);
    context.setPositions(positions);
    context.setVelocities(velocities);

    // Simulate it and see whether the constraints remain satisfied.

    double initialEnergy = 0.0;
    for (int i = 0; i < 1000; ++i) {
        State state = context.getState(State::Positions | State::Energy | State::Velocities | State::Forces);
        for (int j = 0; j < system.getNumConstraints(); ++j) {
            int particle1, particle2;
            double distance;
            system.getConstraintParameters(j, particle1, particle2, distance);
            Vec3 p1 = state.getPositions()[particle1];
            Vec3 p2 = state.getPositions()[particle2];
            double dist = std::sqrt((p1[0]-p2[0])*(p1[0]-p2[0])+(p1[1]-p2[1])*(p1[1]-p2[1])+(p1[2]-p2[2])*(p1[2]-p2[2]));
            ASSERT_EQUAL_TOL(distance, dist, 2e-5);
        }
        double energy = state.getKineticEnergy()+state.getPotentialEnergy();
        if (i == 1)
            initialEnergy = energy;
        else if (i > 1)
            ASSERT_EQUAL_TOL(initialEnergy, energy, 0.01);
        integrator.step(1);
    }
    ASSERT(context.getState(State::Positions).getTime() > 0.1);
}

void testConstrainedMasslessParticles() {
    System system;
    system.addParticle(0.0);
    system.addParticle(1.0);
    system.addConstraint(0, 1, 1.5);
    vector<Vec3> positions(2);
    positions[0] = Vec3(-1, 0, 0);
    positions[1] = Vec3(1, 0, 0);
    VariableVerletIntegrator integrator(0.01);
    bool failed = false;
    try {
        // This should throw an exception.
        
        Context context(system, integrator, platform);
    }
    catch (exception& ex) {
        failed = true;
    }
    ASSERT(failed);
    
    // Now make both particles massless, which should work.
    
    system.setParticleMass(1, 0.0);
    Context context(system, integrator, platform);
    context.setPositions(positions);
    context.setVelocitiesToTemperature(300.0);
    integrator.step(1);
    State state = context.getState(State::Velocities);
    ASSERT_EQUAL(0.0, state.getVelocities()[0][0]);
}

void testArgonBox() {
    const int gridSize = 8;
    const double mass = 40.0;            // Ar atomic mass
    const double temp = 120.0;           // K
    const double epsilon = BOLTZ * temp; // L-J well depth for Ar
    const double sigma = 0.34;           // L-J size for Ar in nm
    const double density = 0.8;          // atoms / sigma^3
    double cellSize = sigma / pow(density, 0.333);
    double boxSize = gridSize * cellSize;
    double cutoff = 2.0 * sigma;

    // Create a box of argon atoms.
    
    System system;
    NonbondedForce* nonbonded = new NonbondedForce();
    vector<Vec3> positions;
    OpenMM_SFMT::SFMT sfmt;
    init_gen_rand(0, sfmt);
    const Vec3 half(0.5, 0.5, 0.5);
    for (int i = 0; i < gridSize; i++) {
        for (int j = 0; j < gridSize; j++) {
            for (int k = 0; k < gridSize; k++) {
                system.addParticle(mass);
                nonbonded->addParticle(0, sigma, epsilon);
                positions.push_back((Vec3(i, j, k) + half + Vec3(genrand_real2(sfmt), genrand_real2(sfmt), genrand_real2(sfmt))*0.1) * cellSize);
            }
        }
    }

    nonbonded->setNonbondedMethod(NonbondedForce::CutoffPeriodic);
    nonbonded->setCutoffDistance(cutoff);
    system.setDefaultPeriodicBoxVectors(Vec3(boxSize, 0, 0), Vec3(0, boxSize, 0), Vec3(0, 0, boxSize));
    system.addForce(nonbonded);

    VariableVerletIntegrator integrator(1e-5);
    Context context(system, integrator, platform);
    context.setPositions(positions);
    context.setVelocitiesToTemperature(temp);

    // Equilibrate.
    
    integrator.stepTo(1.0);

    // Simulate it and see whether energy remains constant.
    
    State state0 = context.getState(State::Energy);
    double initialEnergy = state0.getKineticEnergy() + state0.getPotentialEnergy();
    for (int i = 0; i < 20; i++) {
        double t = 1.0 + 0.05*(i+1);
        integrator.stepTo(t);
        State state = context.getState(State::Energy);
        double energy = state.getKineticEnergy() + state.getPotentialEnergy();
        ASSERT_EQUAL_TOL(initialEnergy, energy, 0.01);
    }
}

void runPlatformTests();

int main(int argc, char* argv[]) {
    try {
        initializeTests(argc, argv);
        testSingleBond();
        testConstraints();
        testConstrainedClusters();
        testConstrainedMasslessParticles();
        testArgonBox();
        runPlatformTests();
    }
    catch(const exception& e) {
        cout << "exception: " << e.what() << endl;
        return 1;
    }
    cout << "Done" << endl;
    return 0;
}

