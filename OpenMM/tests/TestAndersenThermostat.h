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
#include "openmm/AndersenThermostat.h"
#include "openmm/Context.h"
#include "openmm/NonbondedForce.h"
#include "openmm/System.h"
#include "openmm/VerletIntegrator.h"
#include "SimTKOpenMMRealType.h"
#include "sfmt/SFMT.h"
#include <iostream>
#include <vector>

using namespace OpenMM;
using namespace std;

void testTemperature() {
    const int numParticles = 8;
    const double temp = 100.0;
    const double collisionFreq = 10.0;
    const int numSteps = 5000;
    System system;
    VerletIntegrator integrator(0.003);
    NonbondedForce* forceField = new NonbondedForce();
    for (int i = 0; i < numParticles; ++i) {
        system.addParticle(2.0);
        forceField->addParticle((i%2 == 0 ? 1.0 : -1.0), 1.0, 5.0);
    }
    system.addForce(forceField);
    AndersenThermostat* thermostat = new AndersenThermostat(temp, collisionFreq);
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
    double expected = 0.5*numParticles*3*BOLTZ*temp;
    ASSERT_USUALLY_EQUAL_TOL(expected, ke, 0.1);
}

void testConstraints() {
    const int numParticles = 8;
    const double temp = 100.0;
    const double collisionFreq = 10.0;
    const int numSteps = 15000;
    System system;
    VerletIntegrator integrator(0.004);
    NonbondedForce* forceField = new NonbondedForce();
    for (int i = 0; i < numParticles; ++i) {
        system.addParticle(2.0);
        forceField->addParticle((i%2 == 0 ? 1.0 : -1.0), 1.0, 5.0);
    }
    system.addForce(forceField);
    system.addConstraint(0, 1, 1);
    system.addConstraint(1, 2, 1);
    system.addConstraint(2, 3, 1);
    system.addConstraint(3, 0, 1);
    system.addConstraint(4, 5, 1);
    system.addConstraint(5, 6, 1);
    system.addConstraint(6, 7, 1);
    system.addConstraint(7, 4, 1);
    AndersenThermostat* thermostat = new AndersenThermostat(temp, collisionFreq);
    system.addForce(thermostat);
    Context context(system, integrator, platform);
    vector<Vec3> positions(numParticles);
    positions[0] = Vec3(0, 0, 0);
    positions[1] = Vec3(1, 0, 0);
    positions[2] = Vec3(1, 1, 0);
    positions[3] = Vec3(0, 1, 0);
    positions[4] = Vec3(1, 0, 1);
    positions[5] = Vec3(1, 1, 1);
    positions[6] = Vec3(0, 1, 1);
    positions[7] = Vec3(0, 0, 1);
    context.setPositions(positions);
    context.setVelocitiesToTemperature(temp);

    // Let it equilibrate.

    integrator.step(5000);

    // Now run it for a while and see if the temperature is correct.

    double ke = 0.0;
    for (int i = 0; i < numSteps; ++i) {
        State state = context.getState(State::Energy);
        ke += state.getKineticEnergy();
        integrator.step(1);
    }
    ke /= numSteps;
    double expected = 0.5*(numParticles*3-system.getNumConstraints())*BOLTZ*temp;
    ASSERT_USUALLY_EQUAL_TOL(expected, ke, 0.1);
}

void testRandomSeed() {
    const int numParticles = 8;
    const double temp = 100.0;
    const double collisionFreq = 10.0;
    System system;
    VerletIntegrator integrator(0.01);
    NonbondedForce* forceField = new NonbondedForce();
    for (int i = 0; i < numParticles; ++i) {
        system.addParticle(2.0);
        forceField->addParticle((i%2 == 0 ? 1.0 : -1.0), 1.0, 5.0);
    }
    system.addForce(forceField);
    AndersenThermostat* thermostat = new AndersenThermostat(temp, collisionFreq);
    system.addForce(thermostat);
    vector<Vec3> positions(numParticles);
    vector<Vec3> velocities(numParticles);
    for (int i = 0; i < numParticles; ++i) {
        positions[i] = Vec3((i%2 == 0 ? 2 : -2), (i%4 < 2 ? 2 : -2), (i < 4 ? 2 : -2));
        velocities[i] = Vec3(0, 0, 0);
    }

    // Try twice with the same random seed.

    thermostat->setRandomNumberSeed(5);
    Context context(system, integrator, platform);
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

    thermostat->setRandomNumberSeed(10);
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

void runPlatformTests();

int main(int argc, char* argv[]) {
    try {
        initializeTests(argc, argv);
        testTemperature();
        testConstraints();
        testRandomSeed();
        runPlatformTests();
    }
    catch(const exception& e) {
        cout << "exception: " << e.what() << endl;
        return 1;
    }
    cout << "Done" << endl;
    return 0;
}
