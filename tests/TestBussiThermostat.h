/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit.                   *
 * See https://openmm.org/development.                                        *
 *                                                                            *
 * Portions copyright (c) 2025 Stanford University and the Authors.           *
 * Authors: Muhammad Hasyim                                                   *
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
#include "openmm/BussiThermostat.h"
#include "openmm/Context.h"
#include "openmm/NonbondedForce.h"
#include "openmm/OpenMMException.h"
#include "openmm/System.h"
#include "openmm/VerletIntegrator.h"
#include "SimTKOpenMMRealType.h"
#include "sfmt/SFMT.h"
#include <cmath>
#include <iostream>
#include <string>
#include <vector>
#include <set>

using namespace OpenMM;
using namespace std;

/**
 * Test that the BussiThermostat maintains the correct temperature.
 */
void testBussiTemperature() {
    const int numParticles = 8;
    const double temp = 300.0;
    const double tau = 1.0; // ps
    const int numSteps = 5000;
    System system;
    VerletIntegrator integrator(0.002);
    NonbondedForce* forceField = new NonbondedForce();
    for (int i = 0; i < numParticles; ++i) {
        system.addParticle(12.0);
        forceField->addParticle((i%2 == 0 ? 0.5 : -0.5), 0.3, 1.0);
    }
    system.addForce(forceField);
    BussiThermostat* thermostat = new BussiThermostat(temp, tau);
    system.addForce(thermostat);
    ASSERT(!thermostat->usesPeriodicBoundaryConditions());
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
    double expected = 0.5 * numParticles * 3 * BOLTZ * temp;
    ASSERT_USUALLY_EQUAL_TOL(expected, ke, 0.1);
}

/**
 * Test that the BussiThermostat can be applied to a subset of particles.
 */
void testBussiParticleSubset() {
    const int numParticles = 10;
    const double temp = 300.0;
    const double tau = 0.5;
    const int numSteps = 2000;
    System system;
    VerletIntegrator integrator(0.002);
    NonbondedForce* forceField = new NonbondedForce();
    for (int i = 0; i < numParticles; ++i) {
        system.addParticle(12.0);
        forceField->addParticle(0.0, 0.3, 1.0);
    }
    system.addForce(forceField);
    
    // Apply thermostat only to first half of particles
    set<int> thermostatParticles;
    for (int i = 0; i < numParticles/2; ++i)
        thermostatParticles.insert(i);
    
    BussiThermostat* thermostat = new BussiThermostat(temp, tau);
    thermostat->setParticles(thermostatParticles);
    thermostat->setApplyToAllParticles(false);
    system.addForce(thermostat);
    
    Context context(system, integrator, platform);
    vector<Vec3> positions(numParticles);
    vector<Vec3> velocities(numParticles);
    for (int i = 0; i < numParticles; ++i) {
        positions[i] = Vec3(i * 0.5, 0, 0);
        velocities[i] = Vec3(0, 0, 0); // Start with zero velocity
    }
    context.setPositions(positions);
    context.setVelocities(velocities);
    
    // Run and check that thermostat particles gain kinetic energy
    integrator.step(1000);
    
    State state = context.getState(State::Velocities);
    const vector<Vec3>& vels = state.getVelocities();
    
    // Thermostat particles should have non-zero velocities
    double keThermostat = 0.0;
    double keOther = 0.0;
    for (int i = 0; i < numParticles; ++i) {
        double v2 = vels[i][0]*vels[i][0] + vels[i][1]*vels[i][1] + vels[i][2]*vels[i][2];
        if (i < numParticles/2)
            keThermostat += 0.5 * 12.0 * v2;
        else
            keOther += 0.5 * 12.0 * v2;
    }
    
    // Thermostat particles should have significant kinetic energy
    ASSERT(keThermostat > 0.1 * 0.5 * (numParticles/2) * 3 * BOLTZ * temp);
}

/**
 * Test reservoir energy tracking.
 */
void testBussiReservoirEnergy() {
    const int numParticles = 8;
    const double temp = 100.0;
    const double tau = 1.0;
    System system;
    VerletIntegrator integrator(0.002);
    NonbondedForce* forceField = new NonbondedForce();
    for (int i = 0; i < numParticles; ++i) {
        system.addParticle(2.0);
        forceField->addParticle(0.0, 0.3, 1.0);
    }
    system.addForce(forceField);
    BussiThermostat* thermostat = new BussiThermostat(temp, tau);
    system.addForce(thermostat);
    
    Context context(system, integrator, platform);
    vector<Vec3> positions(numParticles);
    for (int i = 0; i < numParticles; ++i)
        positions[i] = Vec3(i * 0.5, 0, 0);
    context.setPositions(positions);
    context.setVelocitiesToTemperature(temp * 0.5); // Start below target temp
    
    // Run simulation
    integrator.step(1000);
    
    // Check that reservoir energy has been tracked
    // Since we started below target temp, thermostat should have added energy
    // (reservoir energy should be negative, meaning energy flowed from reservoir to system)
    double reservoirEnergy = context.getParameter(BussiThermostat::ReservoirEnergyTranslational());
    // Reservoir energy should be non-zero (thermostat exchanged energy with bath)
    ASSERT(std::abs(reservoirEnergy) > 1e-6);
}

/**
 * Test random seed reproducibility.
 */
void testBussiRandomSeed() {
    const int numParticles = 8;
    const double temp = 100.0;
    const double tau = 0.5;
    System system;
    VerletIntegrator integrator(0.002);
    NonbondedForce* forceField = new NonbondedForce();
    for (int i = 0; i < numParticles; ++i) {
        system.addParticle(2.0);
        forceField->addParticle(0.0, 0.3, 1.0);
    }
    system.addForce(forceField);
    BussiThermostat* thermostat = new BussiThermostat(temp, tau);
    thermostat->setRandomNumberSeed(42);
    system.addForce(thermostat);
    
    vector<Vec3> positions(numParticles);
    vector<Vec3> velocities(numParticles);
    for (int i = 0; i < numParticles; ++i) {
        positions[i] = Vec3(i * 0.5, 0, 0);
        velocities[i] = Vec3(0, 0, 0);
    }
    
    // Run twice with same seed
    Context context(system, integrator, platform);
    context.setPositions(positions);
    context.setVelocities(velocities);
    integrator.step(100);
    State state1 = context.getState(State::Velocities);
    
    context.reinitialize();
    context.setPositions(positions);
    context.setVelocities(velocities);
    integrator.step(100);
    State state2 = context.getState(State::Velocities);
    
    // Results should be identical with same seed
    for (int i = 0; i < numParticles; i++) {
        for (int j = 0; j < 3; j++) {
            ASSERT_EQUAL_TOL(state1.getVelocities()[i][j], state2.getVelocities()[i][j], 1e-10);
        }
    }
}

/**
 * Test that Bussi thermostat throws when kinetic energy is zero (HOOMD-compatible behavior).
 */
void testBussiZeroKineticEnergyThrows() {
    // One particle, no net force (isolated), so after Verlet part1 velocity stays zero
    const int numParticles = 1;
    System system;
    system.addParticle(1.0);
    NonbondedForce* nbf = new NonbondedForce();
    nbf->addParticle(0.0, 0.3, 1.0);
    system.addForce(nbf);
    BussiThermostat* thermostat = new BussiThermostat(300.0, 1.0);
    system.addForce(thermostat);
    VerletIntegrator integrator(0.002);
    Context context(system, integrator, platform);
    context.setPositions(std::vector<Vec3>(numParticles, Vec3(0, 0, 0)));
    context.setVelocities(std::vector<Vec3>(numParticles, Vec3(0, 0, 0)));
    bool threw = false;
    try {
        integrator.step(1);
    } catch (const OpenMMException& e) {
        std::string msg(e.what());
        if (msg.find("non-zero initial momenta") != std::string::npos)
            threw = true;
    }
    ASSERT(threw);
}

/**
 * Test that step order (Bussi after half-kick) produces reasonable temperature.
 */
void testBussiStepOrderTemperature() {
    const int numParticles = 8;
    const double temp = 200.0;
    const double tau = 0.5;
    const int numSteps = 3000;
    System system;
    VerletIntegrator integrator(0.002);
    NonbondedForce* forceField = new NonbondedForce();
    for (int i = 0; i < numParticles; ++i) {
        system.addParticle(10.0);
        forceField->addParticle((i % 2 == 0 ? 0.5 : -0.5), 0.3, 1.0);
    }
    system.addForce(forceField);
    BussiThermostat* thermostat = new BussiThermostat(temp, tau);
    system.addForce(thermostat);
    Context context(system, integrator, platform);
    std::vector<Vec3> positions(numParticles);
    for (int i = 0; i < numParticles; ++i)
        positions[i] = Vec3((i % 2) * 1.0, (i % 4 < 2 ? 1 : -1) * 1.0, (i < 4 ? 1 : -1) * 1.0);
    context.setPositions(positions);
    context.setVelocitiesToTemperature(temp);
    integrator.step(2000);
    double ke = 0.0;
    for (int i = 0; i < numSteps; ++i) {
        State state = context.getState(State::Energy);
        ke += state.getKineticEnergy();
        integrator.step(5);
    }
    ke /= numSteps;
    double expected = 0.5 * numParticles * 3 * BOLTZ * temp;
    ASSERT_USUALLY_EQUAL_TOL(expected, ke, 0.15);
}

void runPlatformTests();

int main(int argc, char* argv[]) {
    try {
        initializeTests(argc, argv);
        testBussiTemperature();
        testBussiParticleSubset();
        testBussiReservoirEnergy();
        testBussiRandomSeed();
        testBussiZeroKineticEnergyThrows();
        testBussiStepOrderTemperature();
        runPlatformTests();
    }
    catch(const exception& e) {
        cout << "exception: " << e.what() << endl;
        return 1;
    }
    cout << "Done" << endl;
    return 0;
}
