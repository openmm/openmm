/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2015 Stanford University and the Authors.           *
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
#include "openmm/BrownianIntegrator.h"
#include "openmm/CompoundIntegrator.h"
#include "openmm/Context.h"
#include "openmm/HarmonicBondForce.h"
#include "openmm/LangevinIntegrator.h"
#include "openmm/System.h"
#include "openmm/VerletIntegrator.h"
#include "SimTKOpenMMRealType.h"
#include <iostream>
#include <vector>

using namespace OpenMM;
using namespace std;

const double TOL = 1e-5;

void testChangingIntegrator() {
    System system;
    system.addParticle(2.0);
    system.addParticle(2.0);
    HarmonicBondForce* bonds = new HarmonicBondForce();
    bonds->addBond(0, 1, 1.5, 1);
    system.addForce(bonds);
    CompoundIntegrator integrator;
    integrator.addIntegrator(new VerletIntegrator(0.01));
    integrator.addIntegrator(new LangevinIntegrator(300.0, 10.0, 0.011));
    integrator.addIntegrator(new BrownianIntegrator(300.0, 10.0, 0.012));
    Context context(system, integrator, platform);
    ASSERT_EQUAL(0, integrator.getCurrentIntegrator());
    vector<Vec3> positions(2);
    positions[0] = Vec3(-1, 0, 0);
    positions[1] = Vec3(1, 0, 0);
    for (int iteration = 0; iteration < 2; ++iteration) {
        context.setPositions(positions);

        // First integrate with the Verlet integrator and compare it to the analytical solution.

        const double freq = 1.0;
        State state = context.getState(State::Energy);
        const double initialEnergy = state.getKineticEnergy()+state.getPotentialEnergy();
        for (int i = 0; i < 100; ++i) {
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
        ASSERT_EQUAL_TOL(100*0.01, context.getState(0).getTime(), 1e-5);

        // Switch to the Langevin integrator and make sure that it heats up.

        integrator.setCurrentIntegrator(1);
        integrator.step(100);
        double ke = 0.0;
        for (int i = 0; i < 1000; ++i) {
            integrator.step(10);
            state = context.getState(State::Energy);
            ke += state.getKineticEnergy();
        }
        double expectedKE = 0.5*2*3*BOLTZ*300.0;
        ASSERT_USUALLY_EQUAL_TOL(expectedKE, ke/1000, 0.1);
        ASSERT_EQUAL_TOL(100*0.01+10100*0.011, context.getState(0).getTime(), 1e-5);
        
        // Now reinitialize the context and repeat all of these tests to make sure that works correctly.
        
        context.reinitialize();
        integrator.setCurrentIntegrator(0);
    }
}

void testChangingParameters() {
    System system;
    system.addParticle(1.0);
    CompoundIntegrator integrator;
    integrator.addIntegrator(new VerletIntegrator(0.01));
    integrator.addIntegrator(new LangevinIntegrator(300.0, 10.0, 0.02));
    integrator.addIntegrator(new BrownianIntegrator(300.0, 10.0, 0.03));
    
    // Try getting and setting the step size for different component integrators.
    
    for (int i = 0; i < 3; i++) {
        integrator.setCurrentIntegrator(i);
        ASSERT_EQUAL_TOL(0.01*(i+1), integrator.getStepSize(), 1e-7);
    }
    for (int i = 0; i < 3; i++) {
        integrator.setCurrentIntegrator(i);
        integrator.setStepSize(0.02*(i+1));
        ASSERT_EQUAL_TOL(0.02*(i+1), integrator.getStepSize(), 1e-7);
    }
    for (int i = 0; i < 3; i++) {
        integrator.setCurrentIntegrator(i);
        ASSERT_EQUAL_TOL(0.02*(i+1), integrator.getStepSize(), 1e-7);
    }
    
    // Try getting and setting the constraint tolerance for different component integrators.
    
    for (int i = 0; i < 3; i++) {
        integrator.setCurrentIntegrator(i);
        ASSERT_EQUAL_TOL(1e-5, integrator.getConstraintTolerance(), 1e-7);
    }
    for (int i = 0; i < 3; i++) {
        integrator.setCurrentIntegrator(i);
        integrator.setConstraintTolerance(1e-4*(i+1));
        ASSERT_EQUAL_TOL(1e-4*(i+1), integrator.getConstraintTolerance(), 1e-7);
    }
    for (int i = 0; i < 3; i++) {
        integrator.setCurrentIntegrator(i);
        ASSERT_EQUAL_TOL(1e-4*(i+1), integrator.getConstraintTolerance(), 1e-7);
    }
}

void testDifferentStepSizes() {
    System system;
    system.addParticle(2.0);
    system.addParticle(2.0);
    HarmonicBondForce* bonds = new HarmonicBondForce();
    bonds->addBond(0, 1, 1.5, 1);
    system.addForce(bonds);
    CompoundIntegrator integrator;
    integrator.addIntegrator(new VerletIntegrator(0.005));
    integrator.addIntegrator(new VerletIntegrator(0.01));
    Context context(system, integrator, platform);
    ASSERT_EQUAL(0, integrator.getCurrentIntegrator());
    vector<Vec3> positions(2);
    positions[0] = Vec3(-1, 0, 0);
    positions[1] = Vec3(1, 0, 0);
    context.setPositions(positions);

    // Integrate with the first Verlet integrator and compare it to the analytical solution.

    const double freq = 1.0;
    double expectedTime = 0;
    for (int i = 0; i < 100; ++i) {
        State state = context.getState(State::Positions);
        double time = state.getTime();
        ASSERT_EQUAL_TOL(expectedTime, time, 1e-5);
        double expectedDist = 1.5+0.5*std::cos(freq*time);
        ASSERT_EQUAL_VEC(Vec3(-0.5*expectedDist, 0, 0), state.getPositions()[0], 0.02);
        ASSERT_EQUAL_VEC(Vec3(0.5*expectedDist, 0, 0), state.getPositions()[1], 0.02);
        integrator.step(1);
        expectedTime += 0.005;
    }
    
    // Now switch to the second Verlet integrator which has a different step size.

    integrator.setCurrentIntegrator(1);
    for (int i = 0; i < 100; ++i) {
        State state = context.getState(State::Positions);
        double time = state.getTime();
        ASSERT_EQUAL_TOL(expectedTime, time, 1e-5);
        double expectedDist = 1.5+0.5*std::cos(freq*time);
        ASSERT_EQUAL_VEC(Vec3(-0.5*expectedDist, 0, 0), state.getPositions()[0], 0.02);
        ASSERT_EQUAL_VEC(Vec3(0.5*expectedDist, 0, 0), state.getPositions()[1], 0.02);
        integrator.step(1);
        expectedTime += 0.01;
    }
    
    // Finally, switch back to the first one again.

    integrator.setCurrentIntegrator(0);
    for (int i = 0; i < 100; ++i) {
        State state = context.getState(State::Positions);
        double time = state.getTime();
        ASSERT_EQUAL_TOL(expectedTime, time, 1e-5);
        double expectedDist = 1.5+0.5*std::cos(freq*time);
        ASSERT_EQUAL_VEC(Vec3(-0.5*expectedDist, 0, 0), state.getPositions()[0], 0.02);
        ASSERT_EQUAL_VEC(Vec3(0.5*expectedDist, 0, 0), state.getPositions()[1], 0.02);
        integrator.step(1);
        expectedTime += 0.005;
    }
}

void runPlatformTests();

int main(int argc, char* argv[]) {
    try {
        initializeTests(argc, argv);
        testChangingIntegrator();
        testChangingParameters();
        testDifferentStepSizes();
        runPlatformTests();
    }
    catch(const exception& e) {
        cout << "exception: " << e.what() << endl;
        return 1;
    }
    cout << "Done" << endl;
    return 0;
}
