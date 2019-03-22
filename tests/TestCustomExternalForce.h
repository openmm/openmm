/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2008-2019 Stanford University and the Authors.      *
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
#include "openmm/CustomExternalForce.h"
#include "openmm/System.h"
#include "openmm/VerletIntegrator.h"
#include "SimTKOpenMMRealType.h"
#include <iostream>
#include <vector>

using namespace OpenMM;
using namespace std;

const double TOL = 1e-5;

void testForce() {
    System system;
    system.addParticle(1.0);
    system.addParticle(1.0);
    system.addParticle(1.0);
    VerletIntegrator integrator(0.01);
    CustomExternalForce* forceField = new CustomExternalForce("scale*(x+yscale*(y-y0)^2)");
    forceField->addPerParticleParameter("y0");
    forceField->addPerParticleParameter("yscale");
    forceField->addGlobalParameter("scale", 0.5);
    vector<double> parameters(2);
    parameters[0] = 0.5;
    parameters[1] = 2.0;
    forceField->addParticle(0, parameters);
    parameters[0] = 1.5;
    parameters[1] = 3.0;
    forceField->addParticle(2, parameters);
    system.addForce(forceField);
    ASSERT(!forceField->usesPeriodicBoundaryConditions());
    ASSERT(!system.usesPeriodicBoundaryConditions());
    Context context(system, integrator, platform);
    vector<Vec3> positions(3);
    positions[0] = Vec3(0, 2, 0);
    positions[1] = Vec3(0, 0, 1);
    positions[2] = Vec3(1, 0, 1);
    context.setPositions(positions);
    State state = context.getState(State::Forces | State::Energy);
    {
        const vector<Vec3>& forces = state.getForces();
        ASSERT_EQUAL_VEC(Vec3(-0.5, -0.5*2.0*2.0*1.5, 0), forces[0], TOL);
        ASSERT_EQUAL_VEC(Vec3(0, 0, 0), forces[1], TOL);
        ASSERT_EQUAL_VEC(Vec3(-0.5, 0.5*3.0*2.0*1.5, 0), forces[2], TOL);
        ASSERT_EQUAL_TOL(0.5*(1.0 + 2.0*1.5*1.5 + 3.0*1.5*1.5), state.getPotentialEnergy(), TOL);
    }
    
    // Try changing the parameters and make sure it's still correct.
    
    parameters[0] = 1.4;
    parameters[1] = 3.5;
    forceField->setParticleParameters(1, 2, parameters);
    forceField->updateParametersInContext(context);
    state = context.getState(State::Forces | State::Energy);
    {
        const vector<Vec3>& forces = state.getForces();
        ASSERT_EQUAL_VEC(Vec3(-0.5, -0.5*2.0*2.0*1.5, 0), forces[0], TOL);
        ASSERT_EQUAL_VEC(Vec3(0, 0, 0), forces[1], TOL);
        ASSERT_EQUAL_VEC(Vec3(-0.5, 0.5*3.5*2.0*1.4, 0), forces[2], TOL);
        ASSERT_EQUAL_TOL(0.5*(1.0 + 2.0*1.5*1.5 + 3.5*1.4*1.4), state.getPotentialEnergy(), TOL);
    }
}

void testManyParameters() {
    System system;
    system.addParticle(1.0);
    VerletIntegrator integrator(0.01);
    CustomExternalForce* forceField = new CustomExternalForce("xscale*(x-x0)^2+yscale*(y-y0)^2+zscale*(z-z0)^2");
    forceField->addPerParticleParameter("x0");
    forceField->addPerParticleParameter("y0");
    forceField->addPerParticleParameter("z0");
    forceField->addPerParticleParameter("xscale");
    forceField->addPerParticleParameter("yscale");
    forceField->addPerParticleParameter("zscale");
    vector<double> parameters(6);
    parameters[0] = 1.0;
    parameters[1] = 2.0;
    parameters[2] = 3.0;
    parameters[3] = 0.1;
    parameters[4] = 0.2;
    parameters[5] = 0.3;
    forceField->addParticle(0, parameters);
    system.addForce(forceField);
    Context context(system, integrator, platform);
    vector<Vec3> positions(1);
    positions[0] = Vec3(0, -1, 0);
    context.setPositions(positions);
    State state = context.getState(State::Forces | State::Energy);
    const vector<Vec3>& forces = state.getForces();
    ASSERT_EQUAL_VEC(Vec3(2*0.1*1.0, 2*0.2*3.0, 2*0.3*3.0), forces[0], TOL);
    ASSERT_EQUAL_TOL(0.1*1*1 + 0.2*3*3 + 0.3*3*3, state.getPotentialEnergy(), TOL);
}

void testPeriodic() {
    Vec3 vx(5, 0, 0);
    Vec3 vy(0, 6, 0);
    Vec3 vz(1, 2, 7);
    double x0 = 51, y0 = -17, z0 = 11.2;
    System system;
    system.setDefaultPeriodicBoxVectors(vx, vy, vz);
    system.addParticle(1.0);
    CustomExternalForce* force = new CustomExternalForce("periodicdistance(x, y, z, x0, y0, z0)^2");
    force->addPerParticleParameter("x0");
    force->addPerParticleParameter("y0");
    force->addPerParticleParameter("z0");
    vector<double> params(3);
    params[0] = x0;
    params[1] = y0;
    params[2] = z0;
    force->addParticle(0, params);
    system.addForce(force);
    ASSERT(force->usesPeriodicBoundaryConditions());
    ASSERT(system.usesPeriodicBoundaryConditions());
    VerletIntegrator integrator(0.01);
    Context context(system, integrator, platform);
    vector<Vec3> positions(1);
    positions[0] = Vec3(0, 2, 0);
    context.setPositions(positions);
    for (int i = 0; i < 100; i++) {
        State state = context.getState(State::Positions | State::Forces | State::Energy);

        // Apply periodic boundary conditions to the difference between the two positions.

        Vec3 delta = Vec3(x0, y0, z0)-state.getPositions()[0];
        delta -= vz*floor(delta[2]/vz[2]+0.5);
        delta -= vy*floor(delta[1]/vy[1]+0.5);
        delta -= vx*floor(delta[0]/vx[0]+0.5);

        // Verify that the force and energy are correct.

        ASSERT_EQUAL_VEC(delta*2, state.getForces()[0], TOL);
        ASSERT_EQUAL_TOL(delta.dot(delta), state.getPotentialEnergy(), TOL);
        integrator.step(1);
    }
}

void testZeroPeriodicDistance() {
    Vec3 vx(5, 0, 0);
    Vec3 vy(0, 6, 0);
    Vec3 vz(1, 2, 7);
    double x0 = 51, y0 = -17, z0 = 11.2;
    System system;
    system.setDefaultPeriodicBoxVectors(vx, vy, vz);
    system.addParticle(1.0);
    CustomExternalForce* force = new CustomExternalForce("periodicdistance(x, y, z, x0, y0, z0)^2");
    force->addPerParticleParameter("x0");
    force->addPerParticleParameter("y0");
    force->addPerParticleParameter("z0");
    vector<double> params(3);
    params[0] = x0;
    params[1] = y0;
    params[2] = z0;
    force->addParticle(0, params);
    system.addForce(force);
    ASSERT(force->usesPeriodicBoundaryConditions());
    ASSERT(system.usesPeriodicBoundaryConditions());
    VerletIntegrator integrator(0.01);
    Context context(system, integrator, platform);
    vector<Vec3> positions(1);
    positions[0] = Vec3(x0, y0, z0);
    context.setPositions(positions);

    State state = context.getState(State::Positions | State::Forces | State::Energy);
    vector<Vec3> forces = state.getForces();
    for (int i = 0; i < 3; i++)
        ASSERT_EQUAL(forces[0][i], forces[0][i]);
}

void testIllegalVariable() {
    System system;
    system.addParticle(1.0);
    CustomExternalForce* force = new CustomExternalForce("x+none");
    force->addParticle(0);
    system.addForce(force);
    VerletIntegrator integrator(0.001);
    bool threwException = false;
    try {
        Context(system, integrator, platform);
    }
    catch (const exception& e) {
        threwException = true;
    }
    ASSERT(threwException);
}

void testAtan2() {
    System system;
    system.addParticle(1.0);
    CustomExternalForce* force = new CustomExternalForce("atan2(x, y)");
    force->addParticle(0);
    system.addForce(force);
    VerletIntegrator integrator(0.01);
    Context context(system, integrator, platform);
    vector<Vec3> positions(1);
    positions[0] = Vec3(1.5, -2.1, 1.2);
    context.setPositions(positions);
    State state = context.getState(State::Energy);
    ASSERT_EQUAL_TOL(atan2(positions[0][0], positions[0][1]), state.getPotentialEnergy(), 1e-5);
}

void runPlatformTests();

int main(int argc, char* argv[]) {
    try {
        initializeTests(argc, argv);
        testForce();
        testManyParameters();
        testPeriodic();
        testZeroPeriodicDistance();
        testIllegalVariable();
        testAtan2();
        runPlatformTests();
    }
    catch(const exception& e) {
        cout << "exception: " << e.what() << endl;
        return 1;
    }
    cout << "Done" << endl;
    return 0;
}


