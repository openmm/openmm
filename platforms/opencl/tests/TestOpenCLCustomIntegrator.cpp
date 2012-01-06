/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2008-2011 Stanford University and the Authors.      *
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

/**
 * This tests the OpenCL implementation of CustomIntegrator.
 */

#include "../../../tests/AssertionUtilities.h"
#include "openmm/Context.h"
#include "OpenCLPlatform.h"
#include "openmm/AndersenThermostat.h"
#include "openmm/HarmonicBondForce.h"
#include "openmm/NonbondedForce.h"
#include "openmm/System.h"
#include "openmm/CustomIntegrator.h"
#include "../src/SimTKUtilities/SimTKOpenMMRealType.h"
#include "sfmt/SFMT.h"
#include <iostream>
#include <vector>

using namespace OpenMM;
using namespace std;

const double TOL = 1e-5;

/**
 * Test a simple leapfrog integrator on a single bond.
 */
void testSingleBond() {
    OpenCLPlatform platform;
    System system;
    system.addParticle(2.0);
    system.addParticle(2.0);
    CustomIntegrator integrator(0.01);
    integrator.addComputePerDof("v", "v+dt*f/m");
    integrator.addComputePerDof("x", "x+dt*v");
    HarmonicBondForce* forceField = new HarmonicBondForce();
    forceField->addBond(0, 1, 1.5, 1);
    system.addForce(forceField);
    Context context(system, integrator, platform);
    vector<Vec3> positions(2);
    positions[0] = Vec3(-1, 0, 0);
    positions[1] = Vec3(1, 0, 0);
    context.setPositions(positions);
    
    // This is simply a harmonic oscillator, so compare it to the analytical solution.
    
    const double freq = 1.0;;
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
        ASSERT_EQUAL_TOL(initialEnergy, energy, 0.01);
        integrator.step(1);
    }
}

/**
 * Test an integrator that enforces constraints.
 */
void testConstraints() {
    const int numParticles = 8;
    const double temp = 500.0;
    OpenCLPlatform platform;
    System system;
    CustomIntegrator integrator(0.002);
    integrator.addPerDofVariable("oldx", 0);
    integrator.addComputePerDof("v", "v+dt*f/m");
    integrator.addComputePerDof("oldx", "x");
    integrator.addComputePerDof("x", "x+dt*v");
    integrator.addConstrainPositions();
    integrator.addComputePerDof("v", "(x-oldx)/dt");
    integrator.setConstraintTolerance(1e-5);
    NonbondedForce* forceField = new NonbondedForce();
    for (int i = 0; i < numParticles; ++i) {
        system.addParticle(i%2 == 0 ? 5.0 : 10.0);
        forceField->addParticle((i%2 == 0 ? 0.2 : -0.2), 0.5, 5.0);
    }
    for (int i = 0; i < numParticles-1; ++i)
        system.addConstraint(i, i+1, 1.0);
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
        State state = context.getState(State::Positions | State::Energy);
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
            ASSERT_EQUAL_TOL(initialEnergy, energy, 0.05);
        integrator.step(1);
    }
}

/**
 * Test an integrator that applies constraints directly to velocities.
 */
void testVelocityConstraints() {
    const int numParticles = 8;
    const double temp = 500.0;
    OpenCLPlatform platform;
    System system;
    CustomIntegrator integrator(0.002);
    integrator.addComputePerDof("v", "v+0.5*dt*f/m");
    integrator.addComputePerDof("x", "x+dt*v");
    integrator.addConstrainPositions();
    integrator.addComputePerDof("v", "v+0.5*dt*f/m");
    integrator.addConstrainVelocities();
    integrator.setConstraintTolerance(1e-5);
    NonbondedForce* forceField = new NonbondedForce();
    for (int i = 0; i < numParticles; ++i) {
        system.addParticle(i%2 == 0 ? 5.0 : 10.0);
        forceField->addParticle((i%2 == 0 ? 0.2 : -0.2), 0.5, 5.0);
    }
    for (int i = 0; i < numParticles-1; ++i)
        system.addConstraint(i, i+1, 1.0);
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
        State state = context.getState(State::Positions | State::Energy);
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
            ASSERT_EQUAL_TOL(initialEnergy, energy, 0.05);
        integrator.step(2);
    }
}

/**
 * Test an integrator with an AndersenThermostat to see if updateContextState()
 * is being handled correctly.
 */
void testWithThermostat() {
    const int numParticles = 8;
    const double temp = 100.0;
    const double collisionFreq = 10.0;
    const int numSteps = 10000;
    OpenCLPlatform platform;
    System system;
    CustomIntegrator integrator(0.005);
    integrator.addUpdateContextState();
    integrator.addComputePerDof("v", "v+dt*f/m");
    integrator.addComputePerDof("x", "x+dt*v");
    NonbondedForce* forceField = new NonbondedForce();
    for (int i = 0; i < numParticles; ++i) {
        system.addParticle(2.0);
        forceField->addParticle((i%2 == 0 ? 1.0 : -1.0), 1.0, 5.0);
    }
    system.addForce(forceField);
    AndersenThermostat* thermostat = new AndersenThermostat(temp, collisionFreq);
    system.addForce(thermostat);
    Context context(system, integrator, platform);
    vector<Vec3> positions(numParticles);
    for (int i = 0; i < numParticles; ++i)
        positions[i] = Vec3((i%2 == 0 ? 2 : -2), (i%4 < 2 ? 2 : -2), (i < 4 ? 2 : -2));
    context.setPositions(positions);
    
    // Let it equilibrate.
    
    integrator.step(10000);
    
    // Now run it for a while and see if the temperature is correct.
    
    double ke = 0.0;
    for (int i = 0; i < numSteps; ++i) {
        State state = context.getState(State::Energy);
        ke += state.getKineticEnergy();
        integrator.step(1);
    }
    ke /= numSteps;
    double expected = 0.5*numParticles*3*BOLTZ*temp;
    ASSERT_USUALLY_EQUAL_TOL(expected, ke, 6/std::sqrt((double) numSteps));
}

/**
 * Test a Monte Carlo integrator that uses global variables and depends on energy.
 */
void testMonteCarlo() {
    OpenCLPlatform platform;
    System system;
    system.addParticle(1.0);
    system.addParticle(1.0);
    CustomIntegrator integrator(0.1);
    const double kT = BOLTZ*300.0;
    integrator.addGlobalVariable("kT", kT);
    integrator.addGlobalVariable("oldE", 0);
    integrator.addGlobalVariable("accept", 0);
    integrator.addPerDofVariable("oldx", 0);
    integrator.addComputeGlobal("oldE", "energy");
    integrator.addComputePerDof("oldx", "x");
    integrator.addComputePerDof("x", "x+dt*gaussian");
    integrator.addComputeGlobal("accept", "step(exp((oldE-energy)/kT)-uniform)");
    integrator.addComputePerDof("x", "accept*x + (1-accept)*oldx");
    HarmonicBondForce* forceField = new HarmonicBondForce();
    forceField->addBond(0, 1, 2.0, 10.0);
    system.addForce(forceField);
    Context context(system, integrator, platform);
    vector<Vec3> positions(2);
    positions[0] = Vec3(-1, 0, 0);
    positions[1] = Vec3(1, 0, 0);
    context.setPositions(positions);
    
    // Compute the histogram of distances and see if it satisfies a Boltzmann distribution.
    
    const int numBins = 100;
    const double maxDist = 4.0;
    const int numIterations = 5000;
    vector<int> counts(numBins, 0);
    for (int i = 0; i < numIterations; ++i) {
        integrator.step(10);
        State state = context.getState(State::Positions);
        Vec3 delta = state.getPositions()[0]-state.getPositions()[1];
        double dist = sqrt(delta.dot(delta));
        if (dist < maxDist)
            counts[(int) (numBins*dist/maxDist)]++;
    }
    vector<double> expected(numBins, 0);
    double sum = 0;
    for (int i = 0; i < numBins; i++) {
        double dist = (i+0.5)*maxDist/numBins;
        expected[i] = dist*dist*exp(-5.0*(dist-2)*(dist-2)/kT);
        sum += expected[i];
    }
    for (int i = 0; i < numBins; i++)
        ASSERT_USUALLY_EQUAL_TOL((double) counts[i]/numIterations, expected[i]/sum, 0.01);
}

/**
 * Test the ComputeSum operation.
 */
void testSum() {
    const int numParticles = 200;
    const double boxSize = 10;
    OpenCLPlatform platform;
    System system;
    system.setDefaultPeriodicBoxVectors(Vec3(boxSize, 0, 0), Vec3(0, boxSize, 0), Vec3(0, 0, boxSize));
    NonbondedForce* nb = new NonbondedForce();
    system.addForce(nb);
    vector<Vec3> positions(numParticles);
    OpenMM_SFMT::SFMT sfmt;
    init_gen_rand(0, sfmt);
    for (int i = 0; i < numParticles; i++) {
        system.addParticle(1.5);
        nb->addParticle(i%2 == 0 ? 1 : -1, 0.1, 1);
        bool close = true;
        while (close) {
            positions[i] = Vec3(boxSize*genrand_real2(sfmt), boxSize*genrand_real2(sfmt), boxSize*genrand_real2(sfmt));
            close = false;
            for (int j = 0; j < i; ++j) {
                Vec3 delta = positions[i]-positions[j];
                if (delta.dot(delta) < 0.1)
                    close = true;
            }
        }
    }
    CustomIntegrator integrator(0.01);
    integrator.addGlobalVariable("ke", 0);
    integrator.addComputePerDof("v", "v+dt*f/m");
    integrator.addComputePerDof("x", "x+dt*v");
    integrator.addComputeSum("ke", "m*v*v/2");
    Context context(system, integrator, platform);
    context.setPositions(positions);
    
    // See if the sum is being computed correctly.
    
    State state = context.getState(State::Energy);
    const double initialEnergy = state.getKineticEnergy()+state.getPotentialEnergy();
    for (int i = 0; i < 100; ++i) {
        state = context.getState(State::Energy);
        ASSERT_EQUAL_TOL(state.getKineticEnergy(), integrator.getGlobalVariable(0), 1e-5);
        integrator.step(1);
    }
}

/**
 * Test an integrator that both uses and modifies a context parameter.
 */
void testParameter() {
    OpenCLPlatform platform;
    System system;
    system.addParticle(1.0);
    AndersenThermostat* thermostat = new AndersenThermostat(0.1, 0.1);
    system.addForce(thermostat);
    CustomIntegrator integrator(0.1);
    integrator.addGlobalVariable("temp", 0);
    integrator.addComputeGlobal("temp", "AndersenTemperature");
    integrator.addComputeGlobal("AndersenTemperature", "temp*2");
    Context context(system, integrator, platform);
    
    // See if the parameter is being used correctly.
    
    for (int i = 0; i < 10; i++) {
        integrator.step(1);
        ASSERT_EQUAL_TOL(context.getParameter("AndersenTemperature"), 0.1*(1<<(i+1)), 1e-5);
    }
}

int main() {
    try {
        testSingleBond();
//        testConstraints();
//        testVelocityConstraints();
        testWithThermostat();
        testMonteCarlo();
        testSum();
        testParameter();
    }
    catch(const exception& e) {
        cout << "exception: " << e.what() << endl;
        return 1;
    }
    cout << "Done" << endl;
    return 0;
}
