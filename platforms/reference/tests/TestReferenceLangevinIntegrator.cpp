/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2008 Stanford University and the Authors.           *
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
 * This tests the reference implementation of LangevinIntegrator.
 */

#include "../../../tests/AssertionUtilities.h"
#include "OpenMMContext.h"
#include "ReferencePlatform.h"
#include "StandardMMForceField.h"
#include "System.h"
#include "LangevinIntegrator.h"
#include "../src/SimTKUtilities/SimTKOpenMMRealType.h"
#include <iostream>
#include <vector>

using namespace OpenMM;
using namespace std;

const double TOL = 1e-5;

void testSingleBond() {
    ReferencePlatform platform;
    System system(2, 0);
    system.setAtomMass(0, 2.0);
    system.setAtomMass(1, 2.0);
    LangevinIntegrator integrator(0, 0.1, 0.01);
    StandardMMForceField* forceField = new StandardMMForceField(2, 1, 0, 0, 0);
    forceField->setBondParameters(0, 0, 1, 1.5, 1);
    system.addForce(forceField);
    OpenMMContext context(system, integrator, platform);
    vector<Vec3> positions(2);
    positions[0] = Vec3(-1, 0, 0);
    positions[1] = Vec3(1, 0, 0);
    context.setPositions(positions);
    
    // This is simply a damped harmonic oscillator, so compare it to the analytical solution.
    
    double freq = std::sqrt(1-0.05*0.05);
    for (int i = 0; i < 1000; ++i) {
        State state = context.getState(State::Positions | State::Velocities);
        double time = state.getTime();
        double expectedDist = 1.5+0.5*std::exp(-0.05*time)*std::cos(freq*time);
        ASSERT_EQUAL_VEC(Vec3(-0.5*expectedDist, 0, 0), state.getPositions()[0], 0.02);
        ASSERT_EQUAL_VEC(Vec3(0.5*expectedDist, 0, 0), state.getPositions()[1], 0.02);
        double expectedSpeed = -0.5*std::exp(-0.05*time)*(0.05*std::cos(freq*time)+freq*std::sin(freq*time));
        ASSERT_EQUAL_VEC(Vec3(-0.5*expectedSpeed, 0, 0), state.getVelocities()[0], 0.02);
        ASSERT_EQUAL_VEC(Vec3(0.5*expectedSpeed, 0, 0), state.getVelocities()[1], 0.02);
        integrator.step(1);
    }
    
    // Not set the friction to a tiny value and see if it conserves energy.
    
    integrator.setFriction(5e-5);
    context.setPositions(positions);
    State state = context.getState(State::Energy);
    double initialEnergy = state.getKineticEnergy()+state.getPotentialEnergy();
    for (int i = 0; i < 1000; ++i) {
        state = context.getState(State::Energy);
        double energy = state.getKineticEnergy()+state.getPotentialEnergy();
        ASSERT_EQUAL_TOL(initialEnergy, energy, 0.01);
        integrator.step(1);
    }
}

void testTemperature() {
    const int numAtoms = 8;
    const double temp = 100.0;
    ReferencePlatform platform;
    System system(numAtoms, 0);
    LangevinIntegrator integrator(temp, 2.0, 0.01);
    StandardMMForceField* forceField = new StandardMMForceField(numAtoms, 0, 0, 0, 0);
    for (int i = 0; i < numAtoms; ++i) {
        system.setAtomMass(i, 2.0);
        forceField->setAtomParameters(i, (i%2 == 0 ? 1.0 : -1.0), 1.0, 5.0);
    }
    system.addForce(forceField);
    OpenMMContext context(system, integrator, platform);
    vector<Vec3> positions(numAtoms);
    for (int i = 0; i < numAtoms; ++i)
        positions[i] = Vec3((i%2 == 0 ? 2 : -2), (i%4 < 2 ? 2 : -2), (i < 4 ? 2 : -2));
    context.setPositions(positions);
    
    // Let it equilibrate.
    
    integrator.step(10000);
    
    // Now run it for a while and see if the temperature is correct.
    
    double ke = 0.0;
    for (int i = 0; i < 1000; ++i) {
        State state = context.getState(State::Energy);
        ke += state.getKineticEnergy();
        integrator.step(1);
    }
    ke /= 1000;
    double expected = 0.5*numAtoms*3*BOLTZ*temp;
    ASSERT_EQUAL_TOL(expected, ke, 3*expected/std::sqrt(1000.0));
}

int main() {
    try {
        testSingleBond();
        testTemperature();
    }
    catch(const exception& e) {
        cout << "exception: " << e.what() << endl;
        return 1;
    }
    cout << "Done" << endl;
    return 0;
}
