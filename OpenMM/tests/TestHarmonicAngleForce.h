/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2008-2016 Stanford University and the Authors.s      *
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
#include "openmm/HarmonicAngleForce.h"
#include "openmm/System.h"
#include "openmm/VerletIntegrator.h"
#include "SimTKOpenMMRealType.h"
#include <iostream>
#include <vector>

using namespace OpenMM;
using namespace std;

const double TOL = 1e-5;

void testAngles() {
    System system;
    system.addParticle(1.0);
    system.addParticle(1.0);
    system.addParticle(1.0);
    system.addParticle(1.0);
    VerletIntegrator integrator(0.01);
    HarmonicAngleForce* forceField = new HarmonicAngleForce();
    forceField->addAngle(0, 1, 2, PI_M/3, 1.1);
    forceField->addAngle(1, 2, 3, PI_M/2, 1.2);
    system.addForce(forceField);
    ASSERT(!forceField->usesPeriodicBoundaryConditions());
    ASSERT(!system.usesPeriodicBoundaryConditions());
    Context context(system, integrator, platform);
    vector<Vec3> positions(4);
    positions[0] = Vec3(0, 1, 0);
    positions[1] = Vec3(0, 0, 0);
    positions[2] = Vec3(1, 0, 0);
    positions[3] = Vec3(2, 1, 0);
    context.setPositions(positions);
    State state = context.getState(State::Forces | State::Energy);
    {
        const vector<Vec3>& forces = state.getForces();
        double torque1 = 1.1*PI_M/6;
        double torque2 = 1.2*PI_M/4;
        ASSERT_EQUAL_VEC(Vec3(torque1, 0, 0), forces[0], TOL);
        ASSERT_EQUAL_VEC(Vec3(-0.5*torque2, 0.5*torque2, 0), forces[3], TOL); // reduced by sqrt(2) due to the bond length, another sqrt(2) due to the angle
        ASSERT_EQUAL_VEC(Vec3(forces[0][0]+forces[1][0]+forces[2][0]+forces[3][0], forces[0][1]+forces[1][1]+forces[2][1]+forces[3][1], forces[0][2]+forces[1][2]+forces[2][2]+forces[3][2]), Vec3(0, 0, 0), TOL);
        ASSERT_EQUAL_TOL(0.5*1.1*(PI_M/6)*(PI_M/6) + 0.5*1.2*(PI_M/4)*(PI_M/4), state.getPotentialEnergy(), TOL);
    }
    
    // Try changing the angle parameters and make sure it's still correct.
    
    forceField->setAngleParameters(0, 0, 1, 2, PI_M/3.1, 1.3);
    forceField->setAngleParameters(1, 1, 2, 3, PI_M/2.1, 1.4);
    forceField->updateParametersInContext(context);
    state = context.getState(State::Forces | State::Energy);
    {
        const vector<Vec3>& forces = state.getForces();
        double dtheta1 = (PI_M/2)-(PI_M/3.1);
        double dtheta2 = (3*PI_M/4)-(PI_M/2.1);
        double torque1 = 1.3*dtheta1;
        double torque2 = 1.4*dtheta2;
        ASSERT_EQUAL_VEC(Vec3(torque1, 0, 0), forces[0], TOL);
        ASSERT_EQUAL_VEC(Vec3(-0.5*torque2, 0.5*torque2, 0), forces[3], TOL); // reduced by sqrt(2) due to the bond length, another sqrt(2) due to the angle
        ASSERT_EQUAL_VEC(Vec3(forces[0][0]+forces[1][0]+forces[2][0]+forces[3][0], forces[0][1]+forces[1][1]+forces[2][1]+forces[3][1], forces[0][2]+forces[1][2]+forces[2][2]+forces[3][2]), Vec3(0, 0, 0), TOL);
        ASSERT_EQUAL_TOL(0.5*1.3*dtheta1*dtheta1 + 0.5*1.4*dtheta2*dtheta2, state.getPotentialEnergy(), TOL);
    }
}

void testPeriodic() {
    // Create a force that uses periodic boundary conditions.
    
    System system;
    system.addParticle(1.0);
    system.addParticle(1.0);
    system.addParticle(1.0);
    system.setDefaultPeriodicBoxVectors(Vec3(3, 0, 0), Vec3(0, 1.5, 0), Vec3(0, 0, 3));
    VerletIntegrator integrator(0.01);
    HarmonicAngleForce* angles = new HarmonicAngleForce();
    angles->addAngle(0, 1, 2, PI_M/3, 1.1);
    system.addForce(angles);
    angles->setUsesPeriodicBoundaryConditions(true);
    Context context(system, integrator, platform);
    vector<Vec3> positions(3);
    positions[0] = Vec3(0, 1, 0);
    positions[1] = Vec3(0, 0, 0);
    positions[2] = Vec3(1, 0, 0);
    context.setPositions(positions);
    State state = context.getState(State::Forces | State::Energy);
    const vector<Vec3>& forces = state.getForces();
    double torque = 1.1*PI_M/6;
    ASSERT_EQUAL_VEC(Vec3(2*torque, 0, 0), forces[0], TOL);
    ASSERT_EQUAL_VEC(Vec3(0, -torque, 0), forces[2], TOL);
    ASSERT_EQUAL_TOL(0.5*1.1*(PI_M/6)*(PI_M/6), state.getPotentialEnergy(), TOL);
}

void runPlatformTests();

int main(int argc, char* argv[]) {
    try {
        initializeTests(argc, argv);
        testAngles();
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
