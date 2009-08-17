
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
 * This tests all the different force terms in the reference implementation of CustomNonbondedForce.
 */

#include "../../../tests/AssertionUtilities.h"
#include "openmm/Context.h"
#include "ReferencePlatform.h"
#include "openmm/CustomNonbondedForce.h"
#include "openmm/System.h"
#include "openmm/VerletIntegrator.h"
#include <iostream>
#include <vector>

using namespace OpenMM;
using namespace std;

const double TOL = 1e-5;

void testSimpleExpression() {
    ReferencePlatform platform;
    System system;
    system.addParticle(1.0);
    system.addParticle(1.0);
    VerletIntegrator integrator(0.01);
    CustomNonbondedForce* forceField = new CustomNonbondedForce("-0.1*r^3");
    forceField->addParticle(vector<double>());
    forceField->addParticle(vector<double>());
    system.addForce(forceField);
    Context context(system, integrator, platform);
    vector<Vec3> positions(2);
    positions[0] = Vec3(0, 0, 0);
    positions[1] = Vec3(2, 0, 0);
    context.setPositions(positions);
    State state = context.getState(State::Forces | State::Energy);
    const vector<Vec3>& forces = state.getForces();
    double force = 0.1*3*(2*2);
    ASSERT_EQUAL_VEC(Vec3(-force, 0, 0), forces[0], TOL);
    ASSERT_EQUAL_VEC(Vec3(force, 0, 0), forces[1], TOL);
    ASSERT_EQUAL_TOL(-0.1*(2*2*2), state.getPotentialEnergy(), TOL);
}

void testParameters() {
    ReferencePlatform platform;
    System system;
    system.addParticle(1.0);
    system.addParticle(1.0);
    VerletIntegrator integrator(0.01);
    CustomNonbondedForce* forceField = new CustomNonbondedForce("scale*a*(r*b)^3");
    forceField->addParameter("a", "a1*a2");
    forceField->addParameter("b", "c+b1+b2");
    forceField->addGlobalParameter("scale", 3.0);
    forceField->addGlobalParameter("c", -1.0);
    vector<double> params(2);
    params[0] = 1.5;
    params[1] = 2.0;
    forceField->addParticle(params);
    params[0] = 2.0;
    params[1] = 3.0;
    forceField->addParticle(params);
    system.addForce(forceField);
    Context context(system, integrator, platform);
    vector<Vec3> positions(2);
    positions[0] = Vec3(0, 0, 0);
    positions[1] = Vec3(2, 0, 0);
    context.setPositions(positions);
    context.setParameter("scale", 1.0);
    context.setParameter("c", 0.0);
    State state = context.getState(State::Forces | State::Energy);
    vector<Vec3> forces = state.getForces();
    double force = -3.0*3*5.0*(10*10);
    ASSERT_EQUAL_VEC(Vec3(-force, 0, 0), forces[0], TOL);
    ASSERT_EQUAL_VEC(Vec3(force, 0, 0), forces[1], TOL);
    ASSERT_EQUAL_TOL(3.0*(10*10*10), state.getPotentialEnergy(), TOL);
    context.setParameter("scale", 1.5);
    context.setParameter("c", 1.0);
    state = context.getState(State::Forces | State::Energy);
    forces = state.getForces();
    force = -1.5*3.0*3*6.0*(12*12);
    ASSERT_EQUAL_VEC(Vec3(-force, 0, 0), forces[0], TOL);
    ASSERT_EQUAL_VEC(Vec3(force, 0, 0), forces[1], TOL);
    ASSERT_EQUAL_TOL(1.5*3.0*(12*12*12), state.getPotentialEnergy(), TOL);
}

void testExceptions() {
    ReferencePlatform platform;
    System system;
    VerletIntegrator integrator(0.01);
    CustomNonbondedForce* nonbonded = new CustomNonbondedForce("a*r");
    nonbonded->addParameter("a", "a1+a2");
    vector<double> params(1);
    vector<Vec3> positions(4);
    for (int i = 0; i < 4; i++) {
        system.addParticle(1.0);
        params[0] = i+1;
        nonbonded->addParticle(params);
        positions[i] = Vec3(i, 0, 0);
    }
    nonbonded->addException(0, 1, vector<double>());
    nonbonded->addException(1, 2, vector<double>());
    nonbonded->addException(2, 3, vector<double>());
    params[0] = 0.5;
    nonbonded->addException(0, 2, params);
    nonbonded->addException(1, 3, params);
    system.addForce(nonbonded);
    Context context(system, integrator, platform);
    context.setPositions(positions);
    State state = context.getState(State::Forces | State::Energy);
    const vector<Vec3>& forces = state.getForces();
    ASSERT_EQUAL_VEC(Vec3(0.5+1+4, 0, 0), forces[0], TOL);
    ASSERT_EQUAL_VEC(Vec3(0.5, 0, 0), forces[1], TOL);
    ASSERT_EQUAL_VEC(Vec3(-0.5, 0, 0), forces[2], TOL);
    ASSERT_EQUAL_VEC(Vec3(-(0.5+1+4), 0, 0), forces[3], TOL);
    ASSERT_EQUAL_TOL((1+4)*3+0.5*2+0.5*2, state.getPotentialEnergy(), TOL);
}

void testCutoff() {
    ReferencePlatform platform;
    System system;
    system.addParticle(1.0);
    system.addParticle(1.0);
    system.addParticle(1.0);
    VerletIntegrator integrator(0.01);
    CustomNonbondedForce* forceField = new CustomNonbondedForce("r");
    forceField->addParticle(vector<double>());
    forceField->addParticle(vector<double>());
    forceField->addParticle(vector<double>());
    forceField->setNonbondedMethod(CustomNonbondedForce::CutoffNonPeriodic);
    forceField->setCutoffDistance(2.5);
    system.addForce(forceField);
    Context context(system, integrator, platform);
    vector<Vec3> positions(3);
    positions[0] = Vec3(0, 0, 0);
    positions[1] = Vec3(0, 2, 0);
    positions[2] = Vec3(0, 3, 0);
    context.setPositions(positions);
    State state = context.getState(State::Forces | State::Energy);
    const vector<Vec3>& forces = state.getForces();
    ASSERT_EQUAL_VEC(Vec3(0, 1, 0), forces[0], TOL);
    ASSERT_EQUAL_VEC(Vec3(0, 0, 0), forces[1], TOL);
    ASSERT_EQUAL_VEC(Vec3(0, -1, 0), forces[2], TOL);
    ASSERT_EQUAL_TOL(2.0+1.0, state.getPotentialEnergy(), TOL);
}

void testPeriodic() {
    ReferencePlatform platform;
    System system;
    system.addParticle(1.0);
    system.addParticle(1.0);
    system.addParticle(1.0);
    VerletIntegrator integrator(0.01);
    CustomNonbondedForce* forceField = new CustomNonbondedForce("r");
    forceField->addParticle(vector<double>());
    forceField->addParticle(vector<double>());
    forceField->addParticle(vector<double>());
    forceField->setNonbondedMethod(CustomNonbondedForce::CutoffPeriodic);
    forceField->setCutoffDistance(2.0);
    forceField->setPeriodicBoxVectors(Vec3(4, 0, 0), Vec3(0, 4, 0), Vec3(0, 0, 4));
    system.addForce(forceField);
    Context context(system, integrator, platform);
    vector<Vec3> positions(3);
    positions[0] = Vec3(0, 0, 0);
    positions[1] = Vec3(0, 2.1, 0);
    positions[2] = Vec3(0, 3, 0);
    context.setPositions(positions);
    State state = context.getState(State::Forces | State::Energy);
    const vector<Vec3>& forces = state.getForces();
    ASSERT_EQUAL_VEC(Vec3(0, -2, 0), forces[0], TOL);
    ASSERT_EQUAL_VEC(Vec3(0, 2, 0), forces[1], TOL);
    ASSERT_EQUAL_VEC(Vec3(0, 0, 0), forces[2], TOL);
    ASSERT_EQUAL_TOL(1.9+1+0.9, state.getPotentialEnergy(), TOL);
}

int main() {
    try {
        testSimpleExpression();
        testParameters();
        testExceptions();
        testCutoff();
        testPeriodic();
    }
    catch(const exception& e) {
        cout << "exception: " << e.what() << endl;
        return 1;
    }
    cout << "Done" << endl;
    return 0;
}
