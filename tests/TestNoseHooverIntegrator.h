/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2019-2020 Stanford University and the Authors.      *
 * Authors: Andreas Krämer and Andrew C. Simmonett                            *
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
#include "openmm/CustomExternalForce.h"
#include "openmm/HarmonicBondForce.h"
#include "openmm/NonbondedForce.h"
#include "openmm/System.h"
#include "openmm/NoseHooverIntegrator.h"
#include "openmm/VirtualSite.h"
#include "SimTKOpenMMRealType.h"
#include "sfmt/SFMT.h"
#include <iostream>
#include <algorithm>
#include <numeric>
#include <vector>

using namespace OpenMM;
using namespace std;

const double TOL = 1e-5;

void testSingleBond() {
    System system;
    system.addParticle(2.0);
    system.addParticle(2.0);
    NoseHooverIntegrator integrator(0.01);
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
        ASSERT_EQUAL_TOL(initialEnergy, energy, 0.01);
        integrator.step(1);
    }
    ASSERT_EQUAL_TOL(10.0, context.getState(0).getTime(), 1e-5);
}

void testConstraints() {
    const int numParticles = 8;
    const int numConstraints = 5;
    System system;
    NoseHooverIntegrator integrator(0.0005);
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
        double energy = state.getPotentialEnergy()+state.getKineticEnergy();
        if (i == 1)
            initialEnergy = energy;
        else if (i > 1)
            ASSERT_EQUAL_TOL(initialEnergy, energy, 0.01);
        integrator.step(1);
    }
}

void testConstrainedClusters() {
    const int numParticles = 7;
    System system;
    NoseHooverIntegrator integrator(0.0005);
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
        double energy = state.getPotentialEnergy()+state.getKineticEnergy();
        if (i == 1)
            initialEnergy = energy;
        else if (i > 1)
            ASSERT_EQUAL_TOL(initialEnergy, energy, 0.01);
        integrator.step(1);
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
    NoseHooverIntegrator integrator(0.01);
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

/**
 * Make sure that virtual sites are updated correctly
 */
void testThreeParticleVirtualSite() {
    System system;
    system.addParticle(1.0);
    system.addParticle(1.0);
    system.addParticle(1.0);
    system.addParticle(0.0);
    system.setVirtualSite(3, new ThreeParticleAverageSite(0, 1, 2, 0.2, 0.3, 0.5));
    CustomExternalForce* forceField = new CustomExternalForce("-a*x");
    system.addForce(forceField);
    forceField->addPerParticleParameter("a");
    vector<double> params(1);
    params[0] = 0.1;
    forceField->addParticle(0, params);
    params[0] = 0.2;
    forceField->addParticle(1, params);
    params[0] = 0.3;
    forceField->addParticle(2, params);
    params[0] = 0.4;
    forceField->addParticle(3, params);
    NoseHooverIntegrator integrator(0.002);
    Context context(system, integrator, platform);
    vector<Vec3> positions(4);
    positions[0] = Vec3(0, 0, 0);
    positions[1] = Vec3(1, 0, 0);
    positions[2] = Vec3(0, 1, 0);
    context.setPositions(positions);
    context.applyConstraints(0.0001);
    for (int i = 0; i < 1000; i++) {
        State state = context.getState(State::Positions | State::Forces);
        const vector<Vec3>& pos = state.getPositions();
        ASSERT_EQUAL_VEC(pos[0]*0.2+pos[1]*0.3+pos[2]*0.5, pos[3], 1e-5);
        ASSERT_EQUAL_VEC(Vec3(0.1+0.4*0.2, 0, 0), state.getForces()[0], 1e-5);
        ASSERT_EQUAL_VEC(Vec3(0.2+0.4*0.3, 0, 0), state.getForces()[1], 1e-5);
        ASSERT_EQUAL_VEC(Vec3(0.3+0.4*0.5, 0, 0), state.getForces()[2], 1e-5);
        integrator.step(1);
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

    NoseHooverIntegrator integrator(300, 25, 0.001);
    Context context(system, integrator, platform);
    context.setPositions(positions);
    context.setVelocitiesToTemperature(targetTemperature);
    auto velocities = context.getState(State::Velocities).getVelocities();
    double kineticEnergy = 0;
    for(const auto &v : velocities) kineticEnergy += 0.5 * v.dot(v);
    double temperature = (2*kineticEnergy / (nDoF*BOLTZ));
    ASSERT_USUALLY_EQUAL_TOL(targetTemperature, temperature, 0.01);
}

void testHarmonicOscillator() {
    const double mass = 1.0;
    double temperature = 300;
    double frequency = 1;
    double mts = 1, ys = 1, chain_length = 3;

    System system;    
    system.addParticle(mass);
    vector<Vec3> positions(1);
    positions[0] = Vec3(0.5,0.5,0.5);
    vector<Vec3> velocities(1);
    velocities[0] = Vec3(0, 0, 0);
    auto harmonic_restraint = new CustomExternalForce("0.5*(x^2+y^2+z^2)");
    harmonic_restraint->addParticle(0);
    system.addForce(harmonic_restraint);
    NoseHooverIntegrator integrator(0.001);
    integrator.addThermostat(temperature, frequency, chain_length, mts, ys);
    Context context(system, integrator, platform);
    context.setPositions(positions);
    context.setVelocities(velocities);

    double mean_temperature=0;
    // equilibration
    integrator.step(2000);
    for (size_t i=0; i < 2500; i++){
        integrator.step(10);
        State state = context.getState(State::Energy | State::Positions | State::Velocities);
        double kinetic_energy = state.getKineticEnergy();
        double temp = kinetic_energy/(0.5*3*BOLTZ);
        mean_temperature = (i*mean_temperature + temp)/(i+1);
        double PE = state.getPotentialEnergy();
        double time = state.getTime();
        double energy = kinetic_energy + PE + integrator.computeHeatBathEnergy();
    }
    ASSERT_EQUAL_TOL(temperature, mean_temperature, 0.02);
}

int makeDimerBox(System& system, std::vector<Vec3>& positions, bool constrain=true, int numMolecules=20, double bondLength=0.1){
    double boxLength = 2; // nm
    Vec3 a(boxLength, 0.0, 0.0);
    Vec3 b(0.0, boxLength, 0.0);
    Vec3 c(0.0, 0.0, boxLength);
    double mass = 20;
    double bondForceConstant = 30000; //0.001;
    int numDOF = 0;
    NonbondedForce* forceField = new NonbondedForce();
    HarmonicBondForce* bondForce = new HarmonicBondForce();
    for(int molecule = 0; molecule < numMolecules; ++molecule) {
        int particle1 = system.addParticle(mass);
        int particle2 = system.addParticle(mass);
        forceField->addParticle(0.0, 0.1, 1.0);
        forceField->addParticle(0.0, 0.1, 1.0);
        forceField->addException(particle1, particle2, 0, 0, 0);
        bondForce->addBond(particle1, particle2, bondLength, bondForceConstant);
        numDOF += 6;
        if (constrain) {
            system.addConstraint(particle1, particle2, bondLength);
            numDOF -= 1;
        }
    }
    forceField->setCutoffDistance(.99*boxLength/2);
    forceField->setSwitchingDistance(.88*boxLength/2);
    forceField->setUseSwitchingFunction(true);
    forceField->setUseDispersionCorrection(false);
    forceField->setNonbondedMethod(NonbondedForce::CutoffPeriodic);
    system.addForce(forceField);
    system.addForce(bondForce);
    system.setDefaultPeriodicBoxVectors(a, b, c);
    OpenMM_SFMT::SFMT sfmt;
    init_gen_rand(0, sfmt);
    for (int i = 0; i < numMolecules; i++) {
        while (true) {
            Vec3 pos = Vec3(boxLength*genrand_real2(sfmt), boxLength*genrand_real2(sfmt), boxLength*genrand_real2(sfmt));
            Vec3 pos1 = pos + Vec3(0,0, bondLength/2);
            Vec3 pos2 = pos + Vec3(0,0,-bondLength/2);
            double minDist = 2*boxLength;
            for (int j = 0; j < i; j++) {
                Vec3 delta = pos1-positions[j];
                minDist = std::min(minDist, sqrt(delta.dot(delta)));
                delta = pos2-positions[j];
                minDist = std::min(minDist, sqrt(delta.dot(delta)));
            }
            if (minDist > 0.15) {
                positions[2*i+0] = pos1;
                positions[2*i+1] = pos2;
                break;
            }
        }
    }
    return numDOF;
}

void testDimerBox(bool constrain=true) {
    // Check conservation of system + bath energy for a harmonic oscillator
    int numMolecules = 20;
    double bondLength = 0.1;
    double bondLengthSquared = bondLength * bondLength;
    System system;
    std::vector<Vec3> positions(numMolecules*2);
    int numDOF = makeDimerBox(system, positions, constrain, numMolecules, bondLength);

    bool simpleConstruct = true;
    double temperature = 300; // kelvin
    double collisionFrequency = 200; // 1/ps
    int numMTS = 3;
    int numYS = 3;
    int chainLength = 5;
    auto integrator = simpleConstruct ? NoseHooverIntegrator(temperature, collisionFrequency, 0.001, chainLength, numMTS, numYS)
                                      : NoseHooverIntegrator(0.001);
    if (!simpleConstruct)
        integrator.addThermostat(temperature, collisionFrequency, chainLength, numMTS, numYS);
    Context context(system, integrator, platform);
    context.setPositions(positions);
    context.setVelocitiesToTemperature(temperature);

    int nSteps = 5000;
    double mean_temp = 0.0;
    std::vector<double> energies(nSteps);
    for (int i = 0; i < nSteps; ++i) {
        integrator.step(1);
        State state = context.getState(State::Energy | (constrain ? State::Positions : 0));
        if (constrain) {
            auto positions = state.getPositions();
            for(int i = 0; i < numMolecules; ++i) {
                Vec3 delta = positions[2*i+1] - positions[2*i];
                double dR2 = delta.dot(delta);
                ASSERT_EQUAL_TOL(bondLengthSquared, dR2, 1e-4);
            }
        }
        double KE = state.getKineticEnergy();
        double PE = state.getPotentialEnergy();
        double time = state.getTime();
        double instantaneous_temperature = 2 * KE / (BOLTZ * numDOF);
        mean_temp = (i*mean_temp + instantaneous_temperature)/(i+1);
        double energy = KE + PE + integrator.computeHeatBathEnergy();
        energies[i] = energy;
    }
    double sum = std::accumulate(energies.begin(), energies.end(), 0.0);
    double mean = sum / energies.size();
    double sq_sum = std::inner_product(energies.begin(), energies.end(), energies.begin(), 0.0);
    double std = std::sqrt(sq_sum / energies.size() - mean * mean);
    double relative_std = std / mean;
    // Check mean temperature
    ASSERT_USUALLY_EQUAL_TOL(temperature, mean_temp, 1e-2);
    // Check fluctuation of conserved (total bath + system) energy
    ASSERT_USUALLY_EQUAL_TOL(relative_std, 0, 5e-3);
}

void testCheckpoints() {
    // Create a system with Drude-like particles to be thermostated as a pair, as well as another
    // particle to be thermostated independently, to test all integrator features.
    double timeStep = 0.001;
    NoseHooverIntegrator integrator(timeStep), newIntegrator(timeStep);
    System system;
    double mass = 1;
    system.addParticle(8*mass);
    system.addParticle(mass);
    system.addParticle(5*mass);
    HarmonicBondForce* force = new HarmonicBondForce();
    force->addBond(0, 1, 0.1, 50.0);
    force->addBond(0, 2, 0.1, 50.0);
    system.addForce(force);
    double kineticEnergy = 1e6;
    double temperature=300, collisionFrequency=1, chainLength=3, numMTS=3, numYS=3;
    chainLength = 10;
    integrator.addSubsystemThermostat(std::vector<int>{2}, std::vector<std::pair<int,int>>{{0,1}},  temperature, collisionFrequency, temperature, collisionFrequency,
                                      chainLength, numMTS, numYS);
    newIntegrator.addSubsystemThermostat(std::vector<int>{2}, std::vector<std::pair<int,int>>{{0,1}},  temperature, collisionFrequency, temperature, collisionFrequency,
                                      chainLength, numMTS, numYS);
    Context context(system, integrator, platform);
    Context newContext(system, newIntegrator, platform);
    std::vector<Vec3> positions(3);
    std::vector<Vec3> velocities(3);
    positions[1] = {0.1, 0.0, 0.0};
    velocities[1] = {0.1,0.2,-0.2};
    positions[2] = {-0.1, 0.001, 0.001};
    velocities[2] = {-0.1,0.2,-0.2};
    context.setPositions(positions);
    context.setVelocities(velocities);

    // Run a short simulation and checkpoint..
    integrator.step(500);
    std::stringstream checkpoint; 
    context.createCheckpoint(checkpoint);

    // Now continue the simulation
    integrator.step(5);

    // And try the same, starting from the checkpoint
    newContext.loadCheckpoint(checkpoint);
    newIntegrator.step(5);

    State state1 = context.getState(State::Positions | State::Velocities);
    State state2 = newContext.getState(State::Positions | State::Velocities);
    ASSERT_EQUAL_VEC(state1.getPositions()[0], state2.getPositions()[0], 1e-6); 
    ASSERT_EQUAL_VEC(state1.getPositions()[1], state2.getPositions()[1], 1e-6); 
    ASSERT_EQUAL_VEC(state1.getVelocities()[0], state2.getVelocities()[0], 1e-6); 
    ASSERT_EQUAL_VEC(state1.getVelocities()[1], state2.getVelocities()[1], 1e-6); 
}

void testSaveParameters() {
    // Create a system with Drude-like particles to be thermostated as a pair, as well as another
    // particle to be thermostated independently, to test all integrator features.
    double timeStep = 0.001;
    NoseHooverIntegrator integrator(timeStep), newIntegrator(timeStep);
    System system;
    double mass = 1;
    system.addParticle(8*mass);
    system.addParticle(mass);
    system.addParticle(5*mass);
    HarmonicBondForce* force = new HarmonicBondForce();
    force->addBond(0, 1, 0.1, 50.0);
    force->addBond(0, 2, 0.1, 50.0);
    system.addForce(force);
    double kineticEnergy = 1e6;
    double temperature=300, collisionFrequency=1, chainLength=3, numMTS=3, numYS=3;
    chainLength = 10;
    integrator.addSubsystemThermostat(std::vector<int>{2}, std::vector<std::pair<int,int>>{{0,1}},  temperature, collisionFrequency, temperature, collisionFrequency,
                                      chainLength, numMTS, numYS);
    newIntegrator.addSubsystemThermostat(std::vector<int>{2}, std::vector<std::pair<int,int>>{{0,1}},  temperature, collisionFrequency, temperature, collisionFrequency,
                                      chainLength, numMTS, numYS);
    Context context(system, integrator, platform);
    Context newContext(system, newIntegrator, platform);
    std::vector<Vec3> positions(3);
    std::vector<Vec3> velocities(3);
    positions[1] = {0.1, 0.0, 0.0};
    velocities[1] = {0.1,0.2,-0.2};
    positions[2] = {-0.1, 0.001, 0.001};
    velocities[2] = {-0.1,0.2,-0.2};
    context.setPositions(positions);
    context.setVelocities(velocities);

    // Run a short simulation and save a state..
    integrator.step(500);
    State savedState = context.getState(State::Positions | State::Velocities | State::IntegratorParameters); 

    // Now continue the simulation
    integrator.step(5);

    // And try the same, starting from the state
    newContext.setState(savedState);
    newIntegrator.step(5);

    State state1 = context.getState(State::Positions | State::Velocities);
    State state2 = newContext.getState(State::Positions | State::Velocities);
    ASSERT_EQUAL_VEC(state1.getPositions()[0], state2.getPositions()[0], 1e-6); 
    ASSERT_EQUAL_VEC(state1.getPositions()[1], state2.getPositions()[1], 1e-6); 
    ASSERT_EQUAL_VEC(state1.getVelocities()[0], state2.getVelocities()[0], 1e-6); 
    ASSERT_EQUAL_VEC(state1.getVelocities()[1], state2.getVelocities()[1], 1e-6); 
}

void testAPIChangeNumParticles() {
    bool constrain = true;
    int numMolecules = 20;
    double bondLength = 0.1;
    double bondLengthSquared = bondLength * bondLength;
    System system;
    std::vector<Vec3> positions(numMolecules*2);
    int numDOF = makeDimerBox(system, positions, constrain, numMolecules, bondLength);

}

void testForceGroups() {
    System system;
    system.addParticle(1.0);
    NoseHooverIntegrator integrator(1, 1.0, 0.001);
    integrator.setIntegrationForceGroups(1<<1);
    CustomExternalForce* f1 = new CustomExternalForce("x");
    f1->addParticle(0);
    f1->setForceGroup(1);
    CustomExternalForce* f2 = new CustomExternalForce("y");
    f2->addParticle(0);
    f2->setForceGroup(2);
    system.addForce(f1);
    system.addForce(f2);
    Context context(system, integrator, platform);
    context.setPositions(vector<Vec3>(1));

    // Take one step and verify that the position was updated based only on f1.

    integrator.step(1);
    Vec3 pos = context.getState(State::Positions).getPositions()[0];
    ASSERT(pos[0] < 0);
    ASSERT(pos[1] == 0);
    ASSERT(pos[2] == 0);
}

void runPlatformTests();

int main(int argc, char* argv[]) {
    try {
        initializeTests(argc, argv);
        // Underlying integrator tests
        testSingleBond();
        testConstraints();
        testConstrainedClusters();
        testConstrainedMasslessParticles();
        testThreeParticleVirtualSite();
        testInitialTemperature();
        // Thermostat tests
        testHarmonicOscillator();
        bool constrain;
        constrain = false; testDimerBox(constrain);
        constrain = true; testDimerBox(constrain);
        testCheckpoints();
        testSaveParameters();
        testForceGroups();
        runPlatformTests();
    }
    catch(const exception& e) {
        cout << "exception: " << e.what() << endl;
        return 1;
    }
    cout << "Done" << endl;
    return 0;
}
