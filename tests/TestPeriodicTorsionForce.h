/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2008-2024 Stanford University and the Authors.      *
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
#include "openmm/PeriodicTorsionForce.h"
#include "openmm/System.h"
#include "openmm/VerletIntegrator.h"
#include "SimTKOpenMMRealType.h"
#include <iostream>
#include <vector>

using namespace OpenMM;
using namespace std;

const double TOL = 1e-5;

void testPeriodicTorsions() {
    System system;
    system.addParticle(1.0);
    system.addParticle(1.0);
    system.addParticle(1.0);
    system.addParticle(1.0);
    VerletIntegrator integrator(0.01);
    PeriodicTorsionForce* forceField = new PeriodicTorsionForce();
    forceField->addTorsion(0, 1, 2, 3, 2, PI_M/3, 1.1);
    system.addForce(forceField);
    ASSERT(!forceField->usesPeriodicBoundaryConditions());
    ASSERT(!system.usesPeriodicBoundaryConditions());
    Context context(system, integrator, platform);
    vector<Vec3> positions(4);
    positions[0] = Vec3(0, 1, 0);
    positions[1] = Vec3(0, 0, 0);
    positions[2] = Vec3(1, 0, 0);
    positions[3] = Vec3(1, 0, 2);
    context.setPositions(positions);
    State state = context.getState(State::Forces | State::Energy);
    {
        const vector<Vec3>& forces = state.getForces();
        double torque = -2*1.1*std::sin(2*PI_M/3);
        ASSERT_EQUAL_VEC(Vec3(0, 0, torque), forces[0], TOL);
        ASSERT_EQUAL_VEC(Vec3(0, 0.5*torque, 0), forces[3], TOL);
        ASSERT_EQUAL_VEC(Vec3(forces[0][0]+forces[1][0]+forces[2][0]+forces[3][0], forces[0][1]+forces[1][1]+forces[2][1]+forces[3][1], forces[0][2]+forces[1][2]+forces[2][2]+forces[3][2]), Vec3(0, 0, 0), TOL);
        ASSERT_EQUAL_TOL(1.1*(1+std::cos(2*PI_M/3)), state.getPotentialEnergy(), TOL);
    }
    
    // Try changing the torsion parameters and make sure it's still correct.
    
    forceField->setTorsionParameters(0, 0, 1, 2, 3, 3, PI_M/3.2, 1.3);
    forceField->updateParametersInContext(context);
    state = context.getState(State::Forces | State::Energy);
    {
        const vector<Vec3>& forces = state.getForces();
        double dtheta = (3*PI_M/2)-(PI_M/3.2);
        double torque = -3*1.3*std::sin(dtheta);
        ASSERT_EQUAL_VEC(Vec3(0, 0, torque), forces[0], TOL);
        ASSERT_EQUAL_VEC(Vec3(0, 0.5*torque, 0), forces[3], TOL);
        ASSERT_EQUAL_VEC(Vec3(forces[0][0]+forces[1][0]+forces[2][0]+forces[3][0], forces[0][1]+forces[1][1]+forces[2][1]+forces[3][1], forces[0][2]+forces[1][2]+forces[2][2]+forces[3][2]), Vec3(0, 0, 0), TOL);
        ASSERT_EQUAL_TOL(1.3*(1+std::cos(dtheta)), state.getPotentialEnergy(), TOL);
    }
}

void testPeriodic() {
    // Create a force that uses periodic boundary conditions.
    
    System system;
    system.addParticle(1.0);
    system.addParticle(1.0);
    system.addParticle(1.0);
    system.addParticle(1.0);
    system.setDefaultPeriodicBoxVectors(Vec3(3, 0, 0), Vec3(0, 3, 0), Vec3(0, 0, 3));
    VerletIntegrator integrator(0.01);
    PeriodicTorsionForce* torsions = new PeriodicTorsionForce();
    torsions->addTorsion(0, 1, 2, 3, 2, PI_M/3, 1.1);
    torsions->setUsesPeriodicBoundaryConditions(true);
    system.addForce(torsions);
    Context context(system, integrator, platform);
    vector<Vec3> positions(4);
    positions[0] = Vec3(0, 1, 0);
    positions[1] = Vec3(0, 0, 0);
    positions[2] = Vec3(1, 0, 0);
    positions[3] = Vec3(1, 0, 2);
    context.setPositions(positions);
    State state = context.getState(State::Forces | State::Energy);
    const vector<Vec3>& forces = state.getForces();
    double torque = -2*1.1*std::sin(2*PI_M/3);
    ASSERT_EQUAL_VEC(Vec3(0, 0, torque), forces[0], TOL);
    ASSERT_EQUAL_VEC(Vec3(0, -torque, 0), forces[3], TOL);
    ASSERT_EQUAL_VEC(Vec3(forces[0][0]+forces[1][0]+forces[2][0]+forces[3][0], forces[0][1]+forces[1][1]+forces[2][1]+forces[3][1], forces[0][2]+forces[1][2]+forces[2][2]+forces[3][2]), Vec3(0, 0, 0), TOL);
    ASSERT_EQUAL_TOL(1.1*(1+std::cos(2*PI_M/3)), state.getPotentialEnergy(), TOL);
}

void testParallelComputation() {
    System system;
    const int numParticles = 200;
    for (int i = 0; i < numParticles; i++)
        system.addParticle(1.0);
    PeriodicTorsionForce* force = new PeriodicTorsionForce();
    for (int i = 3; i < numParticles; i++)
        force->addTorsion(i-3, i-2, i-1, i, 2, 1.1, i);
    system.addForce(force);
    vector<Vec3> positions(numParticles);
    for (int i = 0; i < numParticles; i++)
        positions[i] = Vec3(i, i%2, i%3);
    VerletIntegrator integrator1(0.01);
    Context context1(system, integrator1, platform);
    context1.setPositions(positions);
    State state1 = context1.getState(State::Forces | State::Energy);
    VerletIntegrator integrator2(0.01);
    string deviceIndex = platform.getPropertyValue(context1, "DeviceIndex");
    map<string, string> props;
    props["DeviceIndex"] = deviceIndex+","+deviceIndex;
    Context context2(system, integrator2, platform, props);
    context2.setPositions(positions);
    State state2 = context2.getState(State::Forces | State::Energy);
    ASSERT_EQUAL_TOL(state1.getPotentialEnergy(), state2.getPotentialEnergy(), 1e-5);
    for (int i = 0; i < numParticles; i++)
        ASSERT_EQUAL_VEC(state1.getForces()[i], state2.getForces()[i], 1e-5);

    // Try updating some parameters and see if they still match.

    for (int i = 95; i < 102; i++) {
        int p1, p2, p3, p4, periodicity;
        double phase, k;
        force->getTorsionParameters(i, p1, p2, p3, p4, periodicity, phase, k);
        force->setTorsionParameters(i, p1, p2, p3, p4, periodicity, phase+0.1, 2*k);
    }
    force->updateParametersInContext(context1);
    force->updateParametersInContext(context2);
    State state3 = context1.getState(State::Energy);
    State state4 = context2.getState(State::Energy);
    ASSERT_EQUAL_TOL(state3.getPotentialEnergy(), state4.getPotentialEnergy(), 1e-5);
    ASSERT(fabs(state1.getPotentialEnergy()-state3.getPotentialEnergy()) > 0.1);
}

void runPlatformTests();

int main(int argc, char* argv[]) {
    try {
        initializeTests(argc, argv);
        testPeriodicTorsions();
        testPeriodic();
        runPlatformTests();
    }
    catch(const exception& e) {
        cout << "exception: " << e.what() << endl;
        return 1;
    }
    cout << "Done" << endl;
    return 0;
}
