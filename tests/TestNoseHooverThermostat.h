/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2019 Stanford University and the Authors.           *
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
#include "openmm/NoseHooverChain.h"
#include "openmm/NoseHooverIntegrator.h"
#include "openmm/Context.h"
#include "openmm/State.h"
#include "openmm/HarmonicBondForce.h"
#include "openmm/NonbondedForce.h"
#include "openmm/CustomExternalForce.h"
#include "openmm/System.h"
#include "SimTKOpenMMRealType.h"
#include "sfmt/SFMT.h"
#include <iostream>
#include <sstream>
#include <iomanip>
#include <vector>
#include <algorithm>
#include <numeric>

using namespace OpenMM;
using namespace std;
#define DEBUG 0
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
#if DEBUG
        std::cout << " Time: " << std::setw(4) << time 
                  << " Temperature: " << std::setw(8) << std::setprecision(3) << temp
                  << " Mean Temp: " << std::setw(8) << std::setprecision(3) << mean_temperature
                  << " Kinet Energy: " << std::setw(12) << std::setprecision(8) << kinetic_energy
                  << " Poten Energy: " << std::setw(12) << std::setprecision(8) << PE
                  << " Total Energy: " << std::setw(12) << std::setprecision(8) << energy << std::endl;
#endif
    }
#if DEBUG
    std::cout << "Mean Temperature: " << mean_temperature << std::endl;;
#endif
    ASSERT_EQUAL_TOL(temperature, mean_temperature, 0.02);
}

void testDimerBox(bool constrain=true) {
    // Check conservation of system + bath energy for a harmonic oscillator
    System system;
    int numMolecules = 20;
    double boxLength = 2; // nm
    Vec3 a(boxLength, 0.0, 0.0);
    Vec3 b(0.0, boxLength, 0.0);
    Vec3 c(0.0, 0.0, boxLength);
    double mass = 20;
    double bondForceConstant = 30000; //0.001;
    double bondLength = 0.1;
    double bondLengthSquared = bondLength * bondLength;
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
    std::vector<Vec3> positions(numMolecules*2);
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
    bool simpleConstruct = true;
    double temperature = 300; // kelvin
    double collisionFrequency = 25; // 1/ps
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
    integrator.step(1500);

    int nSteps = 2000;
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
#if DEBUG
        if(i % 10 == 0) std::cout << " Time: " << std::setw(4) << time 
                                  << " Temperature: " << std::setw(8) << std::setprecision(3) << temperature
                                  << " Mean Temperature: " << std::setw(8) << std::setprecision(3) << mean_temp
                                  << " Total Energy: " << std::setw(12) << std::setprecision(8) << energy << std::endl;
#endif
        energies[i] = energy;
    }
    double sum = std::accumulate(energies.begin(), energies.end(), 0.0);
    double mean = sum / energies.size();
    double sq_sum = std::inner_product(energies.begin(), energies.end(), energies.begin(), 0.0);
    double std = std::sqrt(sq_sum / energies.size() - mean * mean);
    double relative_std = std / mean;
    //std::cout << "Relative STD of energy " << relative_std << std::endl;
    // Check mean temperature
    ASSERT_EQUAL_TOL(temperature, temperature, 1e-2);
    // Check fluctuation of conserved (total bath + system) energy
    ASSERT_EQUAL_TOL(relative_std, 0, 1e-4);
}

void testCheckpoints() {
    double timeStep = 0.001;
    NoseHooverIntegrator integrator(timeStep), newIntegrator(timeStep);
    System system;
    double mass = 1;
    system.addParticle(mass);
    system.addParticle(mass);
    NonbondedForce* force = new NonbondedForce();
    force->addParticle(0, 0.1, 1.0);
    force->addParticle(0, 0.1, 1.0);
    system.addForce(force);
    double kineticEnergy = 1e6;
    double temperature=300, collisionFrequency=1, chainLength=3, numMTS=3, numYS=3;
    chainLength = 10;
    integrator.addSubsystemThermostat(std::vector<int>(), std::vector<std::pair<int,int>>{{0,1}},  temperature, collisionFrequency, temperature, collisionFrequency,
                                      chainLength, numMTS, numYS);
    newIntegrator.addSubsystemThermostat(std::vector<int>(), std::vector<std::pair<int,int>>{{0,1}},  temperature, collisionFrequency, temperature, collisionFrequency,
                                      chainLength, numMTS, numYS);
    Context context(system, integrator, platform);
    Context newContext(system, newIntegrator, platform);
    std::vector<Vec3> positions(2);
    positions[1] = {0.1,0.1,0.1};
    context.setPositions(positions);
    integrator.step(300);
#if DEBUG
    std::cout << std::endl << std::endl << std::endl << "writing checkpoint";
#endif
    std::stringstream checkpoint; 
    context.createCheckpoint(checkpoint);
    
    State state = context.getState(State::Positions | State::Velocities);
    for (size_t i=0; i<100; i++){
        state = context.getState(State::Positions | State::Velocities);
#if DEBUG
        std::cout << "posvel" << state.getPositions()[0] << " " << state.getVelocities()[0] << std::endl;
#endif
        integrator.step(1);
    }
#if DEBUG
    std::cout << std::endl << std::endl << "loading checkpoint" << std::endl;
#endif
    newContext.loadCheckpoint(checkpoint);

    
    State state2 = context.getState(State::Positions | State::Velocities);
    for (size_t i=0; i<100; i++){
        state2 = newContext.getState(State::Positions | State::Velocities);
#if DEBUG
        std::cout << "posvel" << state2.getPositions()[0] << " " << state2.getVelocities()[0] << std::endl;
#endif
        newIntegrator.step(1);
    }

    for (int i=0; i<3;i++){
        ASSERT_EQUAL_VEC(state.getPositions()[0], state2.getPositions()[0], 1e-6); 
        ASSERT_EQUAL_VEC(state.getPositions()[1], state2.getPositions()[1], 1e-6); 
        ASSERT_EQUAL_VEC(state.getVelocities()[0], state2.getVelocities()[0], 1e-6); 
        ASSERT_EQUAL_VEC(state.getVelocities()[1], state2.getVelocities()[1], 1e-6); 
    }
}

void runPlatformTests();

int main(int argc, char* argv[]) {
    try {
        initializeTests(argc, argv);
        testHarmonicOscillator();
        bool constrain;
        constrain = false; testDimerBox(constrain);
        constrain = true; testDimerBox(constrain);
        testCheckpoints();
        runPlatformTests();
    }
    catch(const exception& e) {
        cout << "exception: " << e.what() << endl;
        return 1;
    }
    cout << "Done" << endl;
    return 0;
}
