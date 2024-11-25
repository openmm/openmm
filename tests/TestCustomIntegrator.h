/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2008-2023 Stanford University and the Authors.      *
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

#ifdef WIN32
  #define _USE_MATH_DEFINES // Needed to get M_PI
#endif
#include "openmm/internal/AssertionUtilities.h"
#include "openmm/internal/CustomIntegratorUtilities.h"
#include "openmm/Context.h"
#include "openmm/AndersenThermostat.h"
#include "openmm/CustomAngleForce.h"
#include "openmm/CustomBondForce.h"
#include "openmm/CustomExternalForce.h"
#include "openmm/CustomIntegrator.h"
#include "openmm/HarmonicBondForce.h"
#include "openmm/MonteCarloBarostat.h"
#include "openmm/NonbondedForce.h"
#include "openmm/System.h"
#include "SimTKOpenMMRealType.h"
#include "sfmt/SFMT.h"
#include <algorithm>
#include <cmath>
#include <iostream>
#include <sstream>
#include <vector>

using namespace OpenMM;
using namespace std;

const double TOL = 1e-5;

/**
 * Test a simple leapfrog integrator on a single bond.
 */
void testSingleBond() {
    System system;
    system.addParticle(2.0);
    system.addParticle(2.0);
    const double dt = 0.01;
    CustomIntegrator integrator(dt);
    integrator.addComputePerDof("v", "v+dt*f/m");
    integrator.addComputePerDof("x", "x+dt*v");
    integrator.setKineticEnergyExpression("m*v1*v1/2; v1=v+0.5*dt*f/m");
    HarmonicBondForce* forceField = new HarmonicBondForce();
    forceField->addBond(0, 1, 1.5, 1);
    system.addForce(forceField);
    Context context(system, integrator, platform);
    vector<Vec3> positions(2);
    positions[0] = Vec3(-1, 0, 0);
    positions[1] = Vec3(1, 0, 0);
    context.setPositions(positions);
    vector<Vec3> velocities(2);
    velocities[0] = Vec3(-0.5*dt*0.5*0.5, 0, 0);
    velocities[1] = Vec3(0.5*dt*0.5*0.5, 0, 0);
    context.setVelocities(velocities);
    
    // This is simply a harmonic oscillator, so compare it to the analytical solution.
    
    const double freq = 1.0;;
    for (int i = 0; i < 1000; ++i) {
        State state = context.getState(State::Positions | State::Velocities | State::Energy);
        double time = state.getTime();
        double expectedDist = 1.5+0.5*std::cos(freq*time);
        ASSERT_EQUAL_VEC(Vec3(-0.5*expectedDist, 0, 0), state.getPositions()[0], 1e-4);
        ASSERT_EQUAL_VEC(Vec3(0.5*expectedDist, 0, 0), state.getPositions()[1], 1e-4);
        double expectedSpeed = -0.5*freq*std::sin(freq*(time-dt/2));
        ASSERT_EQUAL_VEC(Vec3(-0.5*expectedSpeed, 0, 0), state.getVelocities()[0], 1e-4);
        ASSERT_EQUAL_VEC(Vec3(0.5*expectedSpeed, 0, 0), state.getVelocities()[1], 1e-4);
        double energy = state.getKineticEnergy()+state.getPotentialEnergy();
        ASSERT_EQUAL_TOL(0.5*0.5*0.5, energy, 1e-4);
        integrator.step(1);
    }
}

/**
 * Test an integrator that enforces constraints.
 */
void testConstraints() {
    const int numParticles = 8;
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
            ASSERT_EQUAL_TOL(initialEnergy, energy, 0.01);
        integrator.step(1);
    }
}

/**
 * Test an integrator that applies constraints directly to velocities.
 */
void testVelocityConstraints() {
    const int numParticles = 10;
    System system;
    CustomIntegrator integrator(0.002);
    integrator.addPerDofVariable("x1", 0);
    integrator.addComputePerDof("v", "v+0.5*dt*f/m");
    integrator.addComputePerDof("x", "x+dt*v");
    integrator.addComputePerDof("x1", "x");
    integrator.addConstrainPositions();
    integrator.addComputePerDof("v", "v+0.5*dt*f/m+(x-x1)/dt");
    integrator.addConstrainVelocities();
    integrator.setConstraintTolerance(1e-5);
    NonbondedForce* forceField = new NonbondedForce();
    for (int i = 0; i < numParticles; ++i) {
        system.addParticle(i%2 == 0 ? 5.0 : 10.0);
        forceField->addParticle((i%2 == 0 ? 0.2 : -0.2), 0.5, 5.0);
    }
    
    // Constrain the first three particles with SHAKE.
    
    system.addConstraint(0, 1, 1.0);
    system.addConstraint(1, 2, 1.0);
    
    // Constrain the next three with SETTLE.
    
    system.addConstraint(3, 4, 1.0);
    system.addConstraint(5, 4, 1.0);
    system.addConstraint(3, 5, sqrt(2.0));
    
    // Constraint the rest with CCMA.
    
    for (int i = 6; i < numParticles-1; ++i)
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
        integrator.step(2);
        State state = context.getState(State::Positions | State::Velocities | State::Energy);
        for (int j = 0; j < system.getNumConstraints(); ++j) {
            int particle1, particle2;
            double distance;
            system.getConstraintParameters(j, particle1, particle2, distance);
            Vec3 p1 = state.getPositions()[particle1];
            Vec3 p2 = state.getPositions()[particle2];
            double dist = std::sqrt((p1[0]-p2[0])*(p1[0]-p2[0])+(p1[1]-p2[1])*(p1[1]-p2[1])+(p1[2]-p2[2])*(p1[2]-p2[2]));
            ASSERT_EQUAL_TOL(distance, dist, 2e-5);
            if (i > 0) {
                Vec3 v1 = state.getVelocities()[particle1];
                Vec3 v2 = state.getVelocities()[particle2];
                double vel = (v1-v2).dot(p1-p2);
                ASSERT_EQUAL_TOL(0.0, vel, 2e-5);
            }
        }
        double energy = state.getKineticEnergy()+state.getPotentialEnergy();
        if (i == 0)
            initialEnergy = energy;
        else if (i > 0)
            ASSERT_EQUAL_TOL(initialEnergy, energy, 0.01);
    }
}

void testConstrainedMasslessParticles() {
    System system;
    system.addParticle(0.0);
    system.addParticle(1.0);
    system.addConstraint(0, 1, 1.5);
    vector<Vec3> positions(2);
    positions[0] = Vec3(-1, 0, 0);
    positions[1] = Vec3(1, 0, 0);
    CustomIntegrator integrator(0.002);
    integrator.addPerDofVariable("oldx", 0);
    integrator.addComputePerDof("v", "v+dt*f/m");
    integrator.addComputePerDof("oldx", "x");
    integrator.addComputePerDof("x", "x+dt*v");
    integrator.addConstrainPositions();
    integrator.addComputePerDof("v", "(x-oldx)/dt");
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
    State state = context.getState(State::Velocities | State::Positions);
    ASSERT_EQUAL(0.0, state.getVelocities()[0][0]);
}

/**
 * Test an integrator with an AndersenThermostat to see if updateContextState()
 * is being handled correctly.
 */
void testWithThermostat() {
    const int numParticles = 8;
    const double temp = 100.0;
    const double collisionFreq = 10.0;
    const int numSteps = 5000;
    System system;
    CustomIntegrator integrator(0.003);
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
    context.setVelocitiesToTemperature(temp);
    
    // Let it equilibrate.
    
    integrator.step(10000);
    
    // Now run it for a while and see if the temperature is correct.
    
    double ke = 0.0;
    for (int i = 0; i < numSteps; ++i) {
        State state = context.getState(State::Energy);
        ke += state.getKineticEnergy();
        integrator.step(10);
    }
    ke /= numSteps;
    double expected = 0.5*numParticles*3*BOLTZ*temp;
    ASSERT_USUALLY_EQUAL_TOL(expected, ke, 0.1);
}

/**
 * Test a Monte Carlo integrator that uses global variables and depends on energy.
 */
void testMonteCarlo() {
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
    integrator.addComputePerDof("x", "select(accept, x, oldx)");
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
    System system;
    system.setDefaultPeriodicBoxVectors(Vec3(boxSize, 0, 0), Vec3(0, boxSize, 0), Vec3(0, 0, boxSize));
    NonbondedForce* nb = new NonbondedForce();
    system.addForce(nb);
    vector<Vec3> positions(numParticles);
    OpenMM_SFMT::SFMT sfmt;
    init_gen_rand(0, sfmt);
    for (int i = 0; i < numParticles; i++) {
        system.addParticle(i%10 == 0 ? 0.0 : 1.5);
        nb->addParticle(i%2 == 0 ? 0.1 : -0.1, 0.1, 1);
        bool close = true;
        while (close) {
            positions[i] = Vec3(boxSize*genrand_real2(sfmt), boxSize*genrand_real2(sfmt), boxSize*genrand_real2(sfmt));
            close = false;
            for (int j = 0; j < i; ++j) {
                Vec3 delta = positions[i]-positions[j];
                if (delta.dot(delta) < 1)
                    close = true;
            }
        }
    }
    CustomIntegrator integrator(0.005);
    integrator.addGlobalVariable("ke", 0);
    integrator.addGlobalVariable("temp", 0);
    integrator.addComputePerDof("v", "v+dt*f/m");
    integrator.addComputePerDof("x", "x+dt*v");
    integrator.addComputeSum("ke", "m*v*v/2");
    integrator.addComputeGlobal("temp", "ke+dt");
    Context context(system, integrator, platform);
    context.setPositions(positions);
    
    // See if the sum is being computed correctly.
    
    for (int i = 0; i < 100; ++i) {
        integrator.step(1);
        State state = context.getState(State::Energy);
        ASSERT_EQUAL_TOL(state.getKineticEnergy(), integrator.getGlobalVariable(0), 1e-5);
        ASSERT_EQUAL_TOL(integrator.getGlobalVariable(0)+integrator.getStepSize(), integrator.getGlobalVariable(1), 1e-5);
    }
}

/**
 * Test an integrator that both uses and modifies a context parameter.
 */
void testParameter() {
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
        ASSERT_EQUAL_TOL(context.getParameter("AndersenTemperature"), 0.1*(1<<(i+1)), 1e-10);
    }
    context.setParameter("AndersenTemperature", -5.0);
    integrator.step(1);
    ASSERT_EQUAL(-5.0, integrator.getGlobalVariable(0));
    ASSERT_EQUAL(-10.0, context.getParameter("AndersenTemperature"));
}

/**
 * Test random number distributions.
 */
void testRandomDistributions() {
    const int numParticles = 100;
    const int numBins = 20;
    const int numSteps = 100;
    System system;
    for (int i = 0; i < numParticles; i++)
        system.addParticle(1.0);
    CustomIntegrator integrator(0.1);
    integrator.addPerDofVariable("a", 0);
    integrator.addPerDofVariable("b", 0);
    integrator.addComputePerDof("a", "uniform");
    integrator.addComputePerDof("b", "gaussian");
    Context context(system, integrator, platform);
    
    // See if the random numbers are distributed correctly.
    
    vector<int> bins(numBins);
    double mean = 0.0;
    double var = 0.0;
    double skew = 0.0;
    double kurtosis = 0.0;
    vector<Vec3> values;
    for (int i = 0; i < numSteps; i++) {
        integrator.step(1);
        integrator.getPerDofVariable(0, values);
        for (int i = 0; i < numParticles; i++)
            for (int j = 0; j < 3; j++) {
                double v = values[i][j];
                ASSERT(v >= 0 && v < 1);
                bins[(int) (v*numBins)]++;
            }
        integrator.getPerDofVariable(1, values);
        for (int i = 0; i < numParticles; i++)
            for (int j = 0; j < 3; j++) {
                double v = values[i][j];
                mean += v;
                var += v*v;
                skew += v*v*v;
                kurtosis += v*v*v*v;
            }
    }
    
    // Check the distribution of uniform randoms.
    
    int numValues = numParticles*numSteps*3;
    double expected = numValues/(double) numBins;
    double tol = 4*sqrt(expected);
    for (int i = 0; i < numBins; i++)
        ASSERT(bins[i] >= expected-tol && bins[i] <= expected+tol);
    
    // Check the distribution of gaussian randoms.
    
    mean /= numValues;
    var /= numValues;
    skew /= numValues;
    kurtosis /= numValues;
    double c2 = var-mean*mean;
    double c3 = skew-3*var*mean+2*mean*mean*mean;
    double c4 = kurtosis-4*skew*mean-3*var*var+12*var*mean*mean-6*mean*mean*mean*mean;
    ASSERT_EQUAL_TOL(0.0, mean, 3.0/sqrt((double) numValues));
    ASSERT_EQUAL_TOL(1.0, c2, 3.0/pow(numValues, 1.0/3.0));
    ASSERT_EQUAL_TOL(0.0, c3, 3.0/pow(numValues, 1.0/4.0));
    ASSERT_EQUAL_TOL(0.0, c4, 3.0/pow(numValues, 1.0/4.0));
}

/**
 * Test getting and setting per-DOF variables.
 */
void testPerDofVariables() {
    const int numParticles = 200;
    const double boxSize = 10;
    System system;
    system.setDefaultPeriodicBoxVectors(Vec3(boxSize, 0, 0), Vec3(0, boxSize, 0), Vec3(0, 0, boxSize));
    NonbondedForce* nb = new NonbondedForce();
    system.addForce(nb);
    nb->setNonbondedMethod(NonbondedForce::CutoffNonPeriodic);
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
    integrator.addPerDofVariable("temp", 0);
    integrator.addPerDofVariable("pos", 0);
    integrator.addPerDofVariable("computed", 0);
    integrator.addComputePerDof("v", "v+dt*f/m");
    integrator.addComputePerDof("x", "x+dt*v");
    integrator.addComputePerDof("pos", "x");
    integrator.addComputePerDof("computed", "step(v)*log(x^4)");
    Context context(system, integrator, platform);
    context.setPositions(positions);
    vector<Vec3> initialValues(numParticles);
    for (int i = 0; i < numParticles; i++)
        initialValues[i] = Vec3(i+0.1, i+0.2, i+0.3);
    integrator.setPerDofVariable(0, initialValues);
    
    // Run a simulation, then query per-DOF values and see if they are correct.
    
    vector<Vec3> values;
    for (int i = 0; i < 100; ++i) {
        integrator.step(1);
        State state = context.getState(State::Positions | State::Velocities);
        integrator.getPerDofVariable(0, values);
        for (int j = 0; j < numParticles; j++)
            ASSERT_EQUAL_VEC(initialValues[j], values[j], 1e-5);
        integrator.getPerDofVariable(1, values);
        for (int j = 0; j < numParticles; j++)
            ASSERT_EQUAL_VEC(state.getPositions()[j], values[j], 1e-5);
        integrator.getPerDofVariable(2, values);
        for (int j = 0; j < numParticles; j++)
            for (int k = 0; k < 3; k++) {
                if (state.getVelocities()[j][k] < 0) {
                    ASSERT(values[j][k] == 0.0);
                }
                else {
                    double v = state.getPositions()[j][k];
                    ASSERT_EQUAL_TOL(log(v*v*v*v), values[j][k], 1e-5);
                }
            }
    }
}

/**
 * Test evaluating force groups separately.
 */
void testForceGroups() {
    System system;
    system.addParticle(2.0);
    system.addParticle(2.0);
    CustomIntegrator integrator(0.01);
    integrator.addPerDofVariable("outf", 0);
    integrator.addPerDofVariable("outf1", 0);
    integrator.addPerDofVariable("outf2", 0);
    integrator.addGlobalVariable("oute", 0);
    integrator.addGlobalVariable("oute1", 0);
    integrator.addGlobalVariable("oute2", 0);
    integrator.addComputePerDof("outf", "f");
    integrator.addComputePerDof("outf1", "f1");
    integrator.addComputePerDof("outf2", "f2");
    integrator.addComputeGlobal("oute", "energy");
    integrator.addComputeGlobal("oute1", "energy1");
    integrator.addComputeGlobal("oute2", "energy2");
    integrator.setIntegrationForceGroups((1<<1) + (1<<2));
    HarmonicBondForce* bonds = new HarmonicBondForce();
    bonds->addBond(0, 1, 1.5, 1.1);
    bonds->setForceGroup(1);
    system.addForce(bonds);
    NonbondedForce* nb = new NonbondedForce();
    nb->addParticle(0.2, 1, 0);
    nb->addParticle(0.2, 1, 0);
    nb->setForceGroup(2);
    system.addForce(nb);
    CustomExternalForce* external = new CustomExternalForce("x");
    external->addParticle(0);
    external->setForceGroup(3);
    system.addForce(external);
    Context context(system, integrator, platform);
    vector<Vec3> positions(2);
    positions[0] = Vec3(-1, 0, 0);
    positions[1] = Vec3(1, 0, 0);
    context.setPositions(positions);
    
    // See if the various forces are computed correctly.
    
    integrator.step(1);
    vector<Vec3> f, f1, f2;
    double e1 = 0.5*1.1*0.5*0.5;
    double e2 = 138.935456*0.2*0.2/2.0;
    integrator.getPerDofVariable(0, f);
    integrator.getPerDofVariable(1, f1);
    integrator.getPerDofVariable(2, f2);
    ASSERT_EQUAL_VEC(Vec3(1.1*0.5, 0, 0), f1[0], 1e-5);
    ASSERT_EQUAL_VEC(Vec3(-1.1*0.5, 0, 0), f1[1], 1e-5);
    ASSERT_EQUAL_VEC(Vec3(-138.935456*0.2*0.2/4.0, 0, 0), f2[0], 1e-5);
    ASSERT_EQUAL_VEC(Vec3(138.935456*0.2*0.2/4.0, 0, 0), f2[1], 1e-5);
    ASSERT_EQUAL_VEC(f1[0]+f2[0], f[0], 1e-5);
    ASSERT_EQUAL_VEC(f1[1]+f2[1], f[1], 1e-5);
    ASSERT_EQUAL_TOL(e1, integrator.getGlobalVariable(1), 1e-5);
    ASSERT_EQUAL_TOL(e2, integrator.getGlobalVariable(2), 1e-5);
    ASSERT_EQUAL_TOL(e1+e2, integrator.getGlobalVariable(0), 1e-5);
    
    // Make sure they also match the values returned by the Context.
    
    State s = context.getState(State::Forces | State::Energy, false, 6);
    State s1 = context.getState(State::Forces | State::Energy, false, 2);
    State s2 = context.getState(State::Forces | State::Energy, false, 4);
    vector<Vec3> c, c1, c2;
    c = s.getForces();
    c1 = s1.getForces();
    c2 = s2.getForces();
    ASSERT_EQUAL_VEC(f[0], c[0], 1e-5);
    ASSERT_EQUAL_VEC(f[1], c[1], 1e-5);
    ASSERT_EQUAL_VEC(f1[0], c1[0], 1e-5);
    ASSERT_EQUAL_VEC(f1[1], c1[1], 1e-5);
    ASSERT_EQUAL_VEC(f2[0], c2[0], 1e-5);
    ASSERT_EQUAL_VEC(f2[1], c2[1], 1e-5);
    ASSERT_EQUAL_TOL(s.getPotentialEnergy(), integrator.getGlobalVariable(0), 1e-5);
    ASSERT_EQUAL_TOL(s1.getPotentialEnergy(), integrator.getGlobalVariable(1), 1e-5);
    ASSERT_EQUAL_TOL(s2.getPotentialEnergy(), integrator.getGlobalVariable(2), 1e-5);
}

/**
 * Test a multiple time step r-RESPA integrator.
 */
void testRespa() {
    const int numParticles = 8;
    System system;
    system.setDefaultPeriodicBoxVectors(Vec3(4, 0, 0), Vec3(0, 4, 0), Vec3(0, 0, 4));
    CustomIntegrator integrator(0.002);
    integrator.addComputePerDof("v", "v+0.5*dt*f1/m");
    for (int i = 0; i < 2; i++) {
        integrator.addComputePerDof("v", "v+0.5*(dt/2)*f0/m");
        integrator.addComputePerDof("x", "x+(dt/2)*v");
        integrator.addComputePerDof("v", "v+0.5*(dt/2)*f0/m");
    }
    integrator.addComputePerDof("v", "v+0.5*dt*f1/m");
    HarmonicBondForce* bonds = new HarmonicBondForce();
    for (int i = 0; i < numParticles-2; i++)
        bonds->addBond(i, i+1, 1.0, 0.5);
    system.addForce(bonds);
    NonbondedForce* nb = new NonbondedForce();
    nb->setCutoffDistance(2.0);
    nb->setNonbondedMethod(NonbondedForce::Ewald);
    for (int i = 0; i < numParticles; ++i) {
        system.addParticle(i%2 == 0 ? 5.0 : 10.0);
        nb->addParticle((i%2 == 0 ? 0.2 : -0.2), 0.5, 5.0);
    }
    nb->setForceGroup(1);
    nb->setReciprocalSpaceForceGroup(0);
    system.addForce(nb);
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
    
    // Simulate it and monitor energy conservations.
    
    double initialEnergy = 0.0;
    for (int i = 0; i < 1000; ++i) {
        State state = context.getState(State::Energy);
        double energy = state.getKineticEnergy()+state.getPotentialEnergy();
        if (i == 1)
            initialEnergy = energy;
        else if (i > 1)
            ASSERT_EQUAL_TOL(initialEnergy, energy, 0.05);
        integrator.step(2);
    }
}

void testIfBlock() {
    System system;
    system.addParticle(2.0);
    system.addParticle(2.0);
    const double dt = 0.01;
    CustomIntegrator integrator(dt);
    integrator.addGlobalVariable("a", 0);
    integrator.addGlobalVariable("b", 0);
    integrator.addComputeGlobal("b", "1");
    integrator.beginIfBlock("a < 3.5");
    integrator.addComputeGlobal("b", "a+1");
    integrator.endBlock();
    Context context(system, integrator, platform);

    // Set "a" to 1.7 and verify that "b" gets set to a+1.

    integrator.setGlobalVariable(0, 1.7);
    integrator.step(1);
    ASSERT_EQUAL_TOL(2.7, integrator.getGlobalVariable(1), 1e-6);

    // Now set it to a value that should cause the block to be skipped.

    integrator.setGlobalVariable(0, 4.1);
    integrator.step(1);
    ASSERT_EQUAL_TOL(1.0, integrator.getGlobalVariable(1), 1e-6);
}

void testWhileBlock() {
    System system;
    system.addParticle(2.0);
    system.addParticle(2.0);
    const double dt = 0.01;
    CustomIntegrator integrator(dt);
    integrator.addGlobalVariable("a", 0);
    integrator.addGlobalVariable("b", 0);
    integrator.addComputeGlobal("b", "1");
    integrator.beginWhileBlock("b <= a");
    integrator.addComputeGlobal("b", "b+1");
    integrator.endBlock();
    Context context(system, integrator, platform);

    // Try a case where the loop should be skipped.

    integrator.setGlobalVariable(0, -3.3);
    integrator.step(1);
    ASSERT_EQUAL_TOL(1.0, integrator.getGlobalVariable(1), 1e-6);

    // In this case it should be executed exactly once.

    integrator.setGlobalVariable(0, 1.2);
    integrator.step(1);
    ASSERT_EQUAL_TOL(2.0, integrator.getGlobalVariable(1), 1e-6);

    // In this case, it should be executed several times.

    integrator.setGlobalVariable(0, 5.3);
    integrator.step(1);
    ASSERT_EQUAL_TOL(6.0, integrator.getGlobalVariable(1), 1e-6);
}

/**
 * Test modifying a global variable, then using it in a per-DOF computation.
 */
void testChangingGlobal() {
    System system;
    system.addParticle(1.0);
    CustomIntegrator integrator(0.1);
    integrator.addGlobalVariable("g", 0);
    integrator.addPerDofVariable("a", 0);
    integrator.addPerDofVariable("b", 0);
    integrator.addComputeGlobal("g", "g+1");
    integrator.addComputePerDof("a", "0.5");
    integrator.addComputePerDof("b", "a+g");
    Context context(system, integrator, platform);
    
    // See if everything is being calculated correctly..
    
    for (int i = 0; i < 10; i++) {
        integrator.step(1);
        ASSERT_EQUAL_TOL(i+1.0, integrator.getGlobalVariable(0), 1e-5);
        vector<Vec3> values;
        integrator.getPerDofVariable(1, values);
        ASSERT_EQUAL_VEC(Vec3(i+1.5, i+1.5, i+1.5), values[0], 1e-5);
    }
}

/**
 * Test steps that depend on derivatives of the energy with respect to parameters.
 */
void testEnergyParameterDerivatives() {
    System system;
    for (int i = 0; i < 3; i++)
        system.addParticle(1.0);
    
    // Create some custom forces that depend on parameters.
    
    CustomBondForce* bonds = new CustomBondForce("K*(A*r-r0)^2");
    system.addForce(bonds);
    bonds->addGlobalParameter("K", 2.0);
    bonds->addGlobalParameter("A", 1.0);
    bonds->addGlobalParameter("r0", 1.5);
    bonds->addEnergyParameterDerivative("K");
    bonds->addEnergyParameterDerivative("r0");
    bonds->addBond(0, 1);
    bonds->setForceGroup(0);
    CustomAngleForce* angles = new CustomAngleForce("K*(B*theta-theta0)^2");
    system.addForce(angles);
    angles->addGlobalParameter("K", 2.0);
    angles->addGlobalParameter("B", 1.0);
    angles->addGlobalParameter("theta0", M_PI/3);
    angles->addEnergyParameterDerivative("K");
    angles->addEnergyParameterDerivative("theta0");
    angles->addAngle(0, 1, 2);
    angles->setForceGroup(1);
    
    // Create an integrator that records parameter derivatives.
    
    CustomIntegrator integrator(0.1);
    integrator.addGlobalVariable("dEdK", 0.0);
    integrator.addGlobalVariable("dEdr0", 0.0);
    integrator.addPerDofVariable("dEdtheta0", 0.0);
    integrator.addGlobalVariable("dEdK_0", 0.0);
    integrator.addPerDofVariable("dEdr0_0", 0.0);
    integrator.addGlobalVariable("dEdtheta0_0", 0.0);
    integrator.addPerDofVariable("dEdK_1", 0.0);
    integrator.addGlobalVariable("dEdr0_1", 0.0);
    integrator.addGlobalVariable("dEdtheta0_1", 0.0);
    integrator.addComputeGlobal("dEdK", "deriv(energy, K)");
    integrator.addComputeGlobal("dEdr0", "deriv(energy, r0)");
    integrator.addComputePerDof("dEdtheta0", "deriv(energy, theta0)");
    integrator.addComputeGlobal("dEdK_0", "deriv(energy0, K)");
    integrator.addComputePerDof("dEdr0_0", "deriv(energy0, r0)");
    integrator.addComputeGlobal("dEdtheta0_0", "deriv(energy0, theta0)");
    integrator.addComputePerDof("dEdK_1", "deriv(energy1, K)");
    integrator.addComputeGlobal("dEdr0_1", "deriv(energy1, r0)");
    integrator.addComputeGlobal("dEdtheta0_1", "deriv(energy1, theta0)");
    
    // Create a Context.
    
    Context context(system, integrator, platform);
    vector<Vec3> positions(3);
    positions[0] = Vec3(0, 1, 0);
    positions[1] = Vec3(0, 0, 0);
    positions[2] = Vec3(1, 0, 0);
    context.setPositions(positions);
    
    // Check the results.
    
    integrator.step(1);
    vector<Vec3> values;
    double dEdK_0 = (1.0-1.5)*(1.0-1.5);
    double dEdK_1 = (M_PI/2-M_PI/3)*(M_PI/2-M_PI/3);
    ASSERT_EQUAL_TOL(dEdK_0, integrator.getGlobalVariableByName("dEdK_0"), 1e-5);
    integrator.getPerDofVariableByName("dEdK_1", values);
    ASSERT_EQUAL_TOL(dEdK_1, values[0][2], 1e-5);
    ASSERT_EQUAL_TOL(dEdK_0+dEdK_1, integrator.getGlobalVariableByName("dEdK"), 1e-5);
    double dEdr0 = -2.0*2.0*(1.0-1.5);
    integrator.getPerDofVariableByName("dEdr0_0", values);
    ASSERT_EQUAL_TOL(dEdr0, values[1][0], 1e-5);
    ASSERT_EQUAL_TOL(0.0, integrator.getGlobalVariableByName("dEdr0_1"), 1e-5);
    ASSERT_EQUAL_TOL(dEdr0, integrator.getGlobalVariableByName("dEdr0"), 1e-5);
    double dEdtheta0 = -2.0*2.0*(M_PI/2-M_PI/3);
    ASSERT_EQUAL_TOL(0.0, integrator.getGlobalVariableByName("dEdtheta0_0"), 1e-5);
    ASSERT_EQUAL_TOL(dEdtheta0, integrator.getGlobalVariableByName("dEdtheta0_1"), 1e-5);
    integrator.getPerDofVariableByName("dEdtheta0", values);
    ASSERT_EQUAL_TOL(dEdtheta0, values[2][1], 1e-5);
}

/**
 * Test an integrator that modifies the step size.
 */
void testChangeDT() {
    System system;
    system.addParticle(1.0);
    CustomIntegrator integrator(0.5);
    integrator.addGlobalVariable("dt_global", 0.0);
    integrator.addPerDofVariable("dt_dof", 0.0);
    integrator.addComputeGlobal("dt", "dt+1");
    integrator.addComputePerDof("dt_dof", "dt");
    integrator.addComputeGlobal("dt_global", "dt");
    Context context(system, integrator, platform);
    vector<Vec3> positions(1);
    positions[0] = Vec3(0, 0, 0);
    context.setPositions(positions);
    for (int i = 0; i < 5; i++) {
        integrator.step(1);
        double dt = 1.5+i;
        ASSERT_EQUAL_TOL(dt, integrator.getStepSize(), 1e-5);
        ASSERT_EQUAL_TOL(dt, integrator.getGlobalVariable(0), 1e-5);
        vector<Vec3> values;
        integrator.getPerDofVariable(0, values);
        ASSERT_EQUAL_VEC(Vec3(dt, dt, dt), values[0], 1e-5);
    }
    integrator.setStepSize(1.0);
    for (int i = 0; i < 5; i++) {
        integrator.step(1);
        double dt = 2.0+i;
        ASSERT_EQUAL_TOL(dt, integrator.getStepSize(), 1e-5);
        ASSERT_EQUAL_TOL(dt, integrator.getGlobalVariable(0), 1e-5);
        vector<Vec3> values;
        integrator.getPerDofVariable(0, values);
        ASSERT_EQUAL_VEC(Vec3(dt, dt, dt), values[0], 1e-5);
    }
}

/**
 * Test an integrator that uses a tabulated function.
 */
void testTabulatedFunction() {
    System system;
    system.addParticle(1.0);
    CustomIntegrator integrator(1.0);
    integrator.addGlobalVariable("global", 1.5);
    integrator.addPerDofVariable("dof", 0.0);
    integrator.addComputeGlobal("global", "fn(global)");
    integrator.addComputePerDof("dof", "fn(x)");
    vector<double> table;
    table.push_back(10.0);
    table.push_back(20.0);
    integrator.addTabulatedFunction("fn", new Continuous1DFunction(table, 1.0, 2.0));
    Context context(system, integrator, platform);
    vector<Vec3> positions(1);
    positions[0] = Vec3(1.2, 1.3, 1.4);
    context.setPositions(positions);
    integrator.step(1);
    ASSERT_EQUAL_TOL(15.0, integrator.getGlobalVariable(0), 1e-5);
    vector<Vec3> values;
    integrator.getPerDofVariable(0, values);
    ASSERT_EQUAL_VEC(Vec3(12.0, 13.0, 14.0), values[0], 1e-5);
}

/**
 * Test an integrator that alternates repeatedly between force groups.
 */
void testAlternatingGroups() {
    System system;
    system.addParticle(1.0);
    CustomExternalForce* force1 = new CustomExternalForce("-0.5*x");
    force1->addParticle(0);
    system.addForce(force1);
    CustomExternalForce* force2 = new CustomExternalForce("-0.8*y");
    force2->addParticle(0);
    force2->setForceGroup(1);
    system.addForce(force2);
    CustomIntegrator integrator(0.5);
    integrator.addGlobalVariable("savede1", 0.0);
    integrator.addGlobalVariable("savede2", 0.0);
    integrator.addGlobalVariable("savede3", 0.0);
    integrator.addGlobalVariable("savede4", 0.0);
    integrator.addPerDofVariable("savedf1", 0.0);
    integrator.addPerDofVariable("savedf2", 0.0);
    integrator.addPerDofVariable("savedf3", 0.0);
    integrator.addPerDofVariable("savedf4", 0.0);
    integrator.addComputeGlobal("savede1", "energy0");
    integrator.addComputeGlobal("savede2", "energy1");
    integrator.addComputePerDof("savedf1", "f0");
    integrator.addComputePerDof("savedf2", "f1");
    integrator.addComputeGlobal("savede3", "energy0");
    integrator.addComputeGlobal("savede4", "energy1");
    integrator.addComputePerDof("savedf3", "f0");
    integrator.addComputePerDof("savedf4", "f1");
    Context context(system, integrator, platform);
    vector<Vec3> positions(1);
    positions[0] = Vec3(1, 2, 3);
    context.setPositions(positions);
    integrator.step(1);
    vector<Vec3> f;
    for (int i = 0; i < 2; i++) {
        ASSERT_EQUAL_TOL(-0.5*1, integrator.getGlobalVariable(2*i), 1e-5);
        ASSERT_EQUAL_TOL(-0.8*2, integrator.getGlobalVariable(2*i+1), 1e-5);
        integrator.getPerDofVariable(2*i, f);
        ASSERT_EQUAL_VEC(Vec3(0.5, 0, 0), f[0], 1e-5);
        integrator.getPerDofVariable(2*i+1, f);
        ASSERT_EQUAL_VEC(Vec3(0, 0.8, 0), f[0], 1e-5);
    }
}

/**
 * Test that the forces are recomputed when updateContextState() modifies positions.
 */
void testUpdateContextState() {
    const int numParticles = 100;
    const double boxSize = 5.0;
    System system;
    system.setDefaultPeriodicBoxVectors(Vec3(boxSize, 0, 0), Vec3(0, boxSize, 0), Vec3(0, 0, boxSize));
    CustomIntegrator integrator(0.003);
    integrator.addPerDofVariable("force1", 0.0);
    integrator.addPerDofVariable("force2", 0.0);
    integrator.addComputePerDof("force1", "f");
    integrator.addUpdateContextState();
    integrator.addComputePerDof("force2", "f");
    NonbondedForce* nonbonded = new NonbondedForce();
    nonbonded->setNonbondedMethod(NonbondedForce::CutoffPeriodic);
    nonbonded->setCutoffDistance(2.0);
    for (int i = 0; i < numParticles; ++i) {
        system.addParticle(2.0);
        nonbonded->addParticle((i%2 == 0 ? 1.0 : -1.0), 1.0, 5.0);
    }
    system.addForce(nonbonded);
    system.addForce(new MonteCarloBarostat(1.0, 300.0, 1));
    Context context(system, integrator, platform);
    OpenMM_SFMT::SFMT sfmt;
    init_gen_rand(0, sfmt);
    vector<Vec3> positions(numParticles);
    for (int i = 0; i < numParticles; ++i)
        positions[i] = Vec3(genrand_real2(sfmt), genrand_real2(sfmt), genrand_real2(sfmt))*boxSize;
    context.setPositions(positions);
    
    // Make sure the forces change when the barostat accepts a step, and don't change
    // otherwise.
    
    for (int i = 0; i < 50; i++) {
        State state1 = context.getState(0);
        integrator.step(1);
        State state2 = context.getState(0);
        bool changed = (state1.getPeriodicBoxVolume() != state2.getPeriodicBoxVolume());
        vector<Vec3> f1, f2;
        integrator.getPerDofVariable(0, f1);
        integrator.getPerDofVariable(1, f2);
        bool different = false;
        for (int j = 0; j < numParticles; j++)
            if (f1[j] != f2[j])
                different = true;
        ASSERT_EQUAL(changed, different);
    }
}

/**
 * Test using expressions that involve vector functions.
 */
void testVectorFunctions() {
    const int numParticles = 8;
    System system;
    CustomIntegrator integrator(0.001);
    integrator.addGlobalVariable("sumy", 0.0);
    integrator.addPerDofVariable("angular", 0.0);
    integrator.addPerDofVariable("shuffle", 0.0);
    integrator.addPerDofVariable("multicross", 0.0);
    integrator.addPerDofVariable("maxplus", 0.0);
    integrator.addComputeSum("sumy", "x*vector(0, 1, 0)");
    integrator.addComputePerDof("angular", "cross(v, x)");
    integrator.addComputePerDof("shuffle", "dot(vector(_z(x), _x(x), _y(x)), v)");
    integrator.addComputePerDof("multicross", "cross(vector(1, 0, 0), cross(vector(0, 0, 1), vector(1, 0, 0)))");
    integrator.addComputePerDof("maxplus", "max(x, 0.1)+0.5");
    OpenMM_SFMT::SFMT sfmt;
    init_gen_rand(0, sfmt);
    vector<Vec3> positions(numParticles);
    vector<Vec3> velocities(numParticles);
    for (int i = 0; i < numParticles; ++i) {
        system.addParticle(1.0);
        positions[i] = Vec3(genrand_real2(sfmt)-0.5, genrand_real2(sfmt)-0.5, genrand_real2(sfmt)-0.5);
        velocities[i] = Vec3(genrand_real2(sfmt)-0.5, genrand_real2(sfmt)-0.5, genrand_real2(sfmt)-0.5);
    }
    Context context(system, integrator, platform);
    context.setPositions(positions);
    context.setVelocities(velocities);
    integrator.step(1);
    
    // See if the expressions were computed correctly.
    
    double sumy = 0;
    vector<Vec3> angular, shuffle, multicross, maxplus;
    integrator.getPerDofVariable(0, angular);
    integrator.getPerDofVariable(1, shuffle);
    integrator.getPerDofVariable(2, multicross);
    integrator.getPerDofVariable(3, maxplus);
    for (int i = 0; i < numParticles; i++) {
        ASSERT_EQUAL_VEC(velocities[i].cross(positions[i]), angular[i], 1e-5);
        ASSERT_EQUAL_VEC(Vec3(1, 1, 1)*velocities[i].dot(Vec3(positions[i][2], positions[i][0], positions[i][1])), shuffle[i], 1e-5);
        ASSERT_EQUAL_VEC(Vec3(0, 0, 1), multicross[i], 1e-5);
        ASSERT_EQUAL_VEC(Vec3(max(positions[i][0], 0.1)+0.5, max(positions[i][1], 0.1)+0.5, max(positions[i][2], 0.1)+0.5), maxplus[i], 1e-5);
        sumy += positions[i][1];
    }
    ASSERT_EQUAL_TOL(sumy, integrator.getGlobalVariable(0), 1e-5);
}

/**
 * This test records energies at multiple points during the step and checks that
 * they're correct.
 */
void testRecordEnergy() {
    const int numParticles = 8;
    System system;
    CustomIntegrator integrator(0.002);
    integrator.addGlobalVariable("startEnergy", 0);
    integrator.addGlobalVariable("endEnergy", 0);
    integrator.addUpdateContextState();
    integrator.addComputePerDof("v", "v+0.5*dt*f/m");
    integrator.addComputeGlobal("startEnergy", "energy");
    integrator.addComputePerDof("x", "x+dt*v");
    integrator.addComputeGlobal("endEnergy", "energy");
    integrator.addConstrainPositions();
    integrator.addComputePerDof("v", "v+0.5*dt*f/m");
    NonbondedForce* forceField = new NonbondedForce();
    for (int i = 0; i < numParticles; ++i) {
        system.addParticle(i%2 == 0 ? 5.0 : 10.0);
        forceField->addParticle((i%2 == 0 ? 0.2 : -0.2), 0.5, 5.0);
    }
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
    
    // Simulate it and see whether the energies are recorded correctly.
    
    for (int i = 0; i < 10; ++i) {
        double startEnergy = context.getState(State::Energy).getPotentialEnergy();
        integrator.step(1);
        double endEnergy = context.getState(State::Energy).getPotentialEnergy();
        ASSERT_EQUAL_TOL(startEnergy, integrator.getGlobalVariable(0), 1e-6);
        ASSERT_EQUAL_TOL(endEnergy, integrator.getGlobalVariable(1), 1e-6);
    }
}

void testInitialTemperature() {
    // Check temperature initialization for a collection of randomly placed particles
    const int numParticles = 50000;
    const int nDoF = 3 * numParticles;
    const double targetTemperature = 300;
    System system;
    OpenMM_SFMT::SFMT sfmt;
    init_gen_rand(0, sfmt);
    std::vector<Vec3> positions(numParticles);

    for (int i = 0; i < numParticles; i++) {
        system.addParticle(1.0);
        positions[i][0] = genrand_real2(sfmt);
        positions[i][1] = genrand_real2(sfmt);
        positions[i][2] = genrand_real2(sfmt);
    }

    CustomIntegrator integrator(0.001);
    Context context(system, integrator, platform);
    context.setPositions(positions);
    context.setVelocitiesToTemperature(targetTemperature);
    auto velocities = context.getState(State::Velocities).getVelocities();
    double kineticEnergy = 0;
    for(const auto &v : velocities) kineticEnergy += 0.5 * v.dot(v);
    double temperature = (2*kineticEnergy / (nDoF*BOLTZ));
    ASSERT_USUALLY_EQUAL_TOL(targetTemperature, temperature, 0.01);
}

void testCheckpoint() {
    // Test that integrator variables get loaded correctly from checkpoints.
    System system;
    system.addParticle(1.0);
    CustomIntegrator integrator(0.001);
    integrator.addGlobalVariable("a", 1.0);
    integrator.addPerDofVariable("b", 2.0);
    Context context(system, integrator, platform);
    vector<Vec3> positions(1, Vec3());
    context.setPositions(positions);
    integrator.setGlobalVariable(0, 5.0);
    vector<Vec3> b1(1, Vec3(1, 2, 3));
    integrator.setPerDofVariable(0, b1);
    stringstream checkpoint; 
    context.createCheckpoint(checkpoint);
    integrator.setGlobalVariable(0, 10.0);
    vector<Vec3> b2(1, Vec3(4, 5, 6));
    integrator.setPerDofVariable(0, b2);
    context.loadCheckpoint(checkpoint);
    ASSERT_EQUAL(5.0, integrator.getGlobalVariable(0));
    vector<Vec3> b3;
    integrator.getPerDofVariable(0, b3);
    ASSERT_EQUAL_VEC(b1[0], b3[0], 1e-6);
}

void testSaveParameters() {
    // Test that integrator variables get loaded correctly from States.
    System system;
    system.addParticle(1.0);
    CustomIntegrator integrator(0.001);
    integrator.addGlobalVariable("a", 1.0);
    integrator.addPerDofVariable("b", 2.0);
    Context context(system, integrator, platform);
    vector<Vec3> positions(1, Vec3());
    context.setPositions(positions);
    integrator.setGlobalVariable(0, 5.0);
    vector<Vec3> b1(1, Vec3(1, 2, 3));
    integrator.setPerDofVariable(0, b1);
    State savedState = context.getState(State::IntegratorParameters); 
    integrator.setGlobalVariable(0, 10.0);
    vector<Vec3> b2(1, Vec3(4, 5, 6));
    integrator.setPerDofVariable(0, b2);
    context.setState(savedState);
    ASSERT_EQUAL(5.0, integrator.getGlobalVariable(0));
    vector<Vec3> b3;
    integrator.getPerDofVariable(0, b3);
    ASSERT_EQUAL_VEC(b1[0], b3[0], 1e-6);
}

void testAnalyzeComputations() {
    System system;
    system.addParticle(1.0);
    CustomBondForce* bond = new CustomBondForce("scale*r");
    bond->addGlobalParameter("scale", 2.0);
    bond->setForceGroup(1);
    system.addForce(bond);

    // Create a complex integrator with lots of nested blocks and steps that use or invalidate
    // forces or energies.

    CustomIntegrator integrator(0.001);
    integrator.addGlobalVariable("color", 1.5);
    integrator.addPerDofVariable("z", 0);
    integrator.addComputeGlobal("color", "energy");           // 0
    integrator.beginIfBlock("color > 1.0");                   // 1
        integrator.addComputeGlobal("scale", "energy0");      // 2
    integrator.endBlock();                                    // 3
    integrator.beginIfBlock("scale < color");                 // 4
        integrator.addComputePerDof("v", "x");                // 5
    integrator.endBlock();                                    // 6
    integrator.addComputePerDof("z", "f1");                   // 7
    integrator.beginWhileBlock("energy2 > 0");                // 8
        integrator.beginIfBlock("color = 1");                 // 9
            integrator.addComputePerDof("v", "2*z");          // 10
        integrator.endBlock();                                // 11
        integrator.beginIfBlock("color = 2");                 // 12
            integrator.addComputeGlobal("color", "color+1");  // 13
            integrator.addUpdateContextState();               // 14
        integrator.endBlock();                                // 15
    integrator.endBlock();                                    // 16
    integrator.addComputePerDof("x", "x+f");                  // 17

    // Call analyzeComputations() and see if the results are what we expect.

    Context context(system, integrator, platform);
    ContextImpl* contextImpl = *reinterpret_cast<ContextImpl**>(&context);
    vector<vector<Lepton::ParsedExpression> > expressions;
    vector<CustomIntegratorUtilities::Comparison> comparisons;
    vector<int> blockEnd, forceGroup;
    vector<bool> invalidatesForces, needsForces, needsEnergy, computeBoth;
    map<string, Lepton::CustomFunction*> functions;
    CustomIntegratorUtilities::analyzeComputations(*contextImpl, integrator, expressions, comparisons, blockEnd, invalidatesForces,
            needsForces, needsEnergy, computeBoth, forceGroup, functions);
    ASSERT_EQUAL(3, blockEnd[1]);
    ASSERT_EQUAL(6, blockEnd[4]);
    ASSERT_EQUAL(16, blockEnd[8]);
    ASSERT_EQUAL(11, blockEnd[9]);
    ASSERT_EQUAL(15, blockEnd[12]);
    for (int i = 0; i < integrator.getNumComputations(); i++) {
        ASSERT_EQUAL(i == 2 || i == 14 || i == 17, invalidatesForces[i]);
        ASSERT_EQUAL(i == 7 || i == 17, needsForces[i]);
        ASSERT_EQUAL(i == 0 || i == 2 || i == 8, needsEnergy[i]);
        ASSERT_EQUAL(i == 17, computeBoth[i]);
        if (needsForces[i] || needsEnergy[i]) {
            int group = -1;
            if (i == 2)
                group = 0;
            else if (i == 7)
                group = 1;
            else if (i == 8)
                group = 2;
            ASSERT_EQUAL(group, forceGroup[i]);
        }
    }
}

/**
 * Make sure random numbers are computed correctly when steps get merged.
 */
void testMergedRandoms() {
    const int numParticles = 10;
    const int numSteps = 10;
    System system;
    for (int i = 0; i < numParticles; i++)
        system.addParticle(1.0);
    CustomIntegrator integrator(0.1);
    integrator.addPerDofVariable("dofUniform1", 0);
    integrator.addPerDofVariable("dofUniform2", 0);
    integrator.addPerDofVariable("dofGaussian1", 0);
    integrator.addPerDofVariable("dofGaussian2", 0);
    integrator.addGlobalVariable("globalUniform1", 0);
    integrator.addGlobalVariable("globalUniform2", 0);
    integrator.addGlobalVariable("globalGaussian1", 0);
    integrator.addGlobalVariable("globalGaussian2", 0);
    integrator.addComputePerDof("dofUniform1", "uniform");
    integrator.addComputePerDof("dofUniform2", "uniform");
    integrator.addComputePerDof("dofGaussian1", "gaussian");
    integrator.addComputePerDof("dofGaussian2", "gaussian");
    integrator.addComputeGlobal("globalUniform1", "uniform");
    integrator.addComputeGlobal("globalUniform2", "uniform");
    integrator.addComputeGlobal("globalGaussian1", "gaussian");
    integrator.addComputeGlobal("globalGaussian2", "gaussian");
    Context context(system, integrator, platform);
    
    // See if the random numbers are computed correctly.
    
    vector<Vec3> values1, values2;
    for (int i = 0; i < numSteps; i++) {
        integrator.step(1);
        integrator.getPerDofVariable(0, values1);
        integrator.getPerDofVariable(1, values2);
        for (int i = 0; i < numParticles; i++)
            for (int j = 0; j < 3; j++) {
                double v1 = values1[i][j];
                double v2 = values2[i][j];
                ASSERT(v1 >= 0 && v1 < 1);
                ASSERT(v2 >= 0 && v2 < 1);
                ASSERT(v1 != v2);
            }
        integrator.getPerDofVariable(2, values1);
        integrator.getPerDofVariable(3, values2);
        for (int i = 0; i < numParticles; i++)
            for (int j = 0; j < 3; j++) {
                double v1 = values1[i][j];
                double v2 = values2[i][j];
                ASSERT(v1 >= -10 && v1 < 10);
                ASSERT(v2 >= -10 && v2 < 10);
                ASSERT(v1 != v2);
            }
        double v1 = integrator.getGlobalVariable(0);
        double v2 = integrator.getGlobalVariable(1);
        ASSERT(v1 >= 0 && v1 < 1);
        ASSERT(v2 >= 0 && v2 < 1);
        ASSERT(v1 != v2);
        v1 = integrator.getGlobalVariable(2);
        v2 = integrator.getGlobalVariable(3);
        ASSERT(v1 >= -10 && v1 < 10);
        ASSERT(v2 >= -10 && v2 < 10);
        ASSERT(v1 != v2);
    }
}

void runPlatformTests();

int main(int argc, char* argv[]) {
    try {
        initializeTests(argc, argv);
        testSingleBond();
        testConstraints();
        testVelocityConstraints();
        testConstrainedMasslessParticles();
        testWithThermostat();
        testMonteCarlo();
        testSum();
        testParameter();
        testRandomDistributions();
        testPerDofVariables();
        testForceGroups();
        testRespa();
        testIfBlock();
        testWhileBlock();
        testChangingGlobal();
        testEnergyParameterDerivatives();
        testChangeDT();
        testTabulatedFunction();
        testAlternatingGroups();
        testUpdateContextState();
        testVectorFunctions();
        testRecordEnergy();
        testInitialTemperature();
        testCheckpoint();
        testSaveParameters();
        testAnalyzeComputations();
        testMergedRandoms();
        runPlatformTests();
    }
    catch(const exception& e) {
        cout << "exception: " << e.what() << endl;
        return 1;
    }
    cout << "Done" << endl;
    return 0;
}
