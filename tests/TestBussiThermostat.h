/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular dynamics toolkit.                     *
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
 * OTHERWISE, ARISING FROM, OR IN CONNECTION WITH THE SOFTWARE OR THE  *
 * USE OR OTHER DEALINGS IN THE SOFTWARE.                                     *
 * -------------------------------------------------------------------------- */

#include "openmm/internal/AssertionUtilities.h"
#include "openmm/internal/BussiRescale.h"
#include "openmm/BussiThermostat.h"
#include "openmm/Context.h"
#include "openmm/NonbondedForce.h"
#include "openmm/OpenMMException.h"
#include "openmm/System.h"
#include "openmm/VerletIntegrator.h"
#include "SimTKOpenMMRealType.h"
#include <cmath>
#include <exception>
#include <iostream>
#include <string>
#include <vector>
#include <set>

using namespace OpenMM;
using namespace std;

/**
 * Deterministic unit test: Bussi rescale algebra with fixed (R1, rGamma) draws.
 * Expected values computed with the same BOLTZ definition as SimTKOpenMMRealType.h.
 */
void testBussiRescaleMath() {
    {
        BussiRescale::Result r = BussiRescale::computeRescale(3, 10.0, 300.0, 0.5, 0.002, 0.25, 2.0);
        ASSERT_EQUAL_TOL(1.0081690737510292, r.alphaSquared, 1e-12);
        ASSERT_EQUAL_TOL(1.0040762290538647, r.alpha, 1e-12);
        ASSERT_EQUAL_TOL(-0.08169073751029154, r.deltaE, 1e-12);
    }
    {
        BussiRescale::Result r = BussiRescale::computeRescale(1, 5.0, 200.0, 1.0, 0.002, -0.1, 0.0);
        ASSERT_EQUAL_TOL(0.9943634407376757, r.alphaSquared, 1e-12);
        ASSERT_EQUAL_TOL(0.9971777377868379, r.alpha, 1e-12);
        ASSERT_EQUAL_TOL(0.028182796311621572, r.deltaE, 1e-12);
    }
}

/**
 * Test that the BussiThermostat can be applied to a subset of particles.
 */
void testBussiParticleSubset() {
    const int numParticles = 10;
    const double temp = 300.0;
    const double tau = 0.5;
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
    for (int i = 0; i < numParticles; ++i)
        positions[i] = Vec3(i * 0.5, 0, 0);
    context.setPositions(positions);
    // Non-zero initial KE so classic Verlet (Bussi before kick) does not throw;
    // split Verlet tolerates v=0 but we support both orderings.
    context.setVelocitiesToTemperature(1.0);
    
    integrator.step(1000);
    
    State state = context.getState(State::Velocities);
    const vector<Vec3>& vels = state.getVelocities();
    
    double keThermostat = 0.0;
    double keOther = 0.0;
    for (int i = 0; i < numParticles; ++i) {
        double v2 = vels[i][0]*vels[i][0] + vels[i][1]*vels[i][1] + vels[i][2]*vels[i][2];
        if (i < numParticles/2)
            keThermostat += 0.5 * 12.0 * v2;
        else
            keOther += 0.5 * 12.0 * v2;
    }
    
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
    
    integrator.step(1000);
    
    double reservoirEnergy = context.getParameter(BussiThermostat::ReservoirEnergyTranslational());
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
    for (int i = 0; i < numParticles; ++i)
        positions[i] = Vec3(i * 0.5, 0, 0);
    vector<Vec3> velocities(numParticles);
    for (int i = 0; i < numParticles; ++i)
        velocities[i] = Vec3(0.01 * (i + 1), -0.02 * (i + 1), 0.015 * (i + 1));
    
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
    auto matchesBussiZeroKeMessage = [](const std::exception& e) {
        const char* w = e.what();
        if (w == nullptr)
            return false;
        std::string msg(w);
        return msg.find("non-zero initial momenta") != std::string::npos;
    };
    try {
        integrator.step(1);
    } catch (const OpenMMException& e) {
        if (matchesBussiZeroKeMessage(e))
            threw = true;
    } catch (const std::exception& e) {
        if (matchesBussiZeroKeMessage(e))
            threw = true;
    } catch (...) {
        threw = true;
    }
    ASSERT(threw);
}

void runPlatformTests();

int main(int argc, char* argv[]) {
    try {
        initializeTests(argc, argv);
        testBussiRescaleMath();
        testBussiParticleSubset();
        testBussiReservoirEnergy();
        testBussiRandomSeed();
        testBussiZeroKineticEnergyThrows();
        runPlatformTests();
    }
    catch(const exception& e) {
        cout << "exception: " << e.what() << endl;
        return 1;
    }
    cout << "Done" << endl;
    return 0;
}
