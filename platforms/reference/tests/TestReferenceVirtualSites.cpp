/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2012 Stanford University and the Authors.           *
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
 * This tests the reference implementation of virtual sites.
 */

#include "../../../tests/AssertionUtilities.h"
#include "openmm/Context.h"
#include "ReferencePlatform.h"
#include "openmm/CustomBondForce.h"
#include "openmm/CustomExternalForce.h"
#include "openmm/LangevinIntegrator.h"
#include "openmm/NonbondedForce.h"
#include "openmm/System.h"
#include "openmm/VerletIntegrator.h"
#include "openmm/VirtualSite.h"
#include "sfmt/SFMT.h"
#include <iostream>
#include <vector>

using namespace OpenMM;
using namespace std;

/**
 * Check that massless particles are handled correctly.
 */
void testMasslessParticle() {
    ReferencePlatform platform;
    System system;
    system.addParticle(0.0);
    system.addParticle(1.0);
    CustomBondForce* bonds = new CustomBondForce("-1/r");
    system.addForce(bonds);
    vector<double> params;
    bonds->addBond(0, 1, params);
    VerletIntegrator integrator(0.002);
    Context context(system, integrator, platform);
    vector<Vec3> positions(2);
    positions[0] = Vec3(0, 0, 0);
    positions[1] = Vec3(1, 0, 0);
    context.setPositions(positions);
    vector<Vec3> velocities(2);
    velocities[0] = Vec3(0, 0, 0);
    velocities[1] = Vec3(0, 1, 0);
    context.setVelocities(velocities);
    
    // The second particle should move in a circular orbit around the first one.
    // Compare it to the analytical solution.
    
    for (int i = 0; i < 1000; ++i) {
        State state = context.getState(State::Positions | State::Velocities | State::Forces);
        double time = state.getTime();
        ASSERT_EQUAL_VEC(Vec3(), state.getPositions()[0], 0.0);
        ASSERT_EQUAL_VEC(Vec3(), state.getVelocities()[0], 0.0);
        ASSERT_EQUAL_VEC(Vec3(cos(time), sin(time), 0), state.getPositions()[1], 0.01);
        ASSERT_EQUAL_VEC(Vec3(-sin(time), cos(time), 0), state.getVelocities()[1], 0.01);
        integrator.step(1);
    }
}

/**
 * Test a TwoParticleAverage virtual site.
 */
void testTwoParticleAverage() {
    ReferencePlatform platform;
    System system;
    system.addParticle(1.0);
    system.addParticle(1.0);
    system.addParticle(0.0);
    system.setVirtualSite(2, new VirtualSite::TwoParticleAverage(0, 1, 0.8, 0.2));
    CustomExternalForce* forceField = new CustomExternalForce("-a*x");
    system.addForce(forceField);
    forceField->addPerParticleParameter("a");
    vector<double> params(1);
    params[0] = 0.1;
    forceField->addParticle(0, params);
    params[0] = 0.2;
    forceField->addParticle(1, params);
    params[0] = 0.3;
    forceField->addParticle(2, params);
    LangevinIntegrator integrator(300.0, 0.1, 0.002);
    Context context(system, integrator, platform);
    vector<Vec3> positions(3);
    positions[0] = Vec3(0, 0, 0);
    positions[1] = Vec3(1, 0, 0);
    context.setPositions(positions);
    context.applyConstraints(0.0001);
    for (int i = 0; i < 1000; i++) {
        State state = context.getState(State::Positions | State::Forces);
        const vector<Vec3>& pos = state.getPositions();
        ASSERT_EQUAL_VEC(pos[0]*0.8+pos[1]*0.2, pos[2], 1e-10);
        ASSERT_EQUAL_VEC(Vec3(0.1+0.3*0.8, 0, 0), state.getForces()[0], 1e-10);
        ASSERT_EQUAL_VEC(Vec3(0.2+0.3*0.2, 0, 0), state.getForces()[1], 1e-10);
    }
}

/**
 * Test a ThreeParticleAverage virtual site.
 */
void testThreeParticleAverage() {
    ReferencePlatform platform;
    System system;
    system.addParticle(1.0);
    system.addParticle(1.0);
    system.addParticle(1.0);
    system.addParticle(0.0);
    system.setVirtualSite(3, new VirtualSite::ThreeParticleAverage(0, 1, 2, 0.2, 0.3, 0.5));
    CustomExternalForce* forceField = new CustomExternalForce("-a*x");
    system.addForce(forceField);
    forceField->addPerParticleParameter("a");
    vector<double> params(1);
    params[0] = 0.1;
    forceField->addParticle(0, params);
    params[0] = 0.2;
    forceField->addParticle(1, params);
    params[0] = 0.3;
    forceField->addParticle(2, params);
    params[0] = 0.4;
    forceField->addParticle(3, params);
    LangevinIntegrator integrator(300.0, 0.1, 0.002);
    Context context(system, integrator, platform);
    vector<Vec3> positions(4);
    positions[0] = Vec3(0, 0, 0);
    positions[1] = Vec3(1, 0, 0);
    positions[2] = Vec3(0, 1, 0);
    context.setPositions(positions);
    context.applyConstraints(0.0001);
    for (int i = 0; i < 1000; i++) {
        State state = context.getState(State::Positions | State::Forces);
        const vector<Vec3>& pos = state.getPositions();
        ASSERT_EQUAL_VEC(pos[0]*0.2+pos[1]*0.3+pos[2]*0.5, pos[3], 1e-10);
        ASSERT_EQUAL_VEC(Vec3(0.1+0.4*0.2, 0, 0), state.getForces()[0], 1e-10);
        ASSERT_EQUAL_VEC(Vec3(0.2+0.4*0.3, 0, 0), state.getForces()[1], 1e-10);
        ASSERT_EQUAL_VEC(Vec3(0.3+0.4*0.5, 0, 0), state.getForces()[2], 1e-10);
    }
}

/**
 * Test an OutOfPlane virtual site.
 */
void testOutOfPlane() {
    ReferencePlatform platform;
    System system;
    system.addParticle(1.0);
    system.addParticle(1.0);
    system.addParticle(1.0);
    system.addParticle(0.0);
    system.setVirtualSite(3, new VirtualSite::OutOfPlane(0, 1, 2, 0.3, 0.4, 0.5));
    CustomExternalForce* forceField = new CustomExternalForce("-a*x");
    system.addForce(forceField);
    forceField->addPerParticleParameter("a");
    vector<double> params(1);
    params[0] = 0.1;
    forceField->addParticle(0, params);
    params[0] = 0.2;
    forceField->addParticle(1, params);
    params[0] = 0.3;
    forceField->addParticle(2, params);
    params[0] = 0.4;
    forceField->addParticle(3, params);
    LangevinIntegrator integrator(300.0, 0.1, 0.002);
    Context context(system, integrator, platform);
    vector<Vec3> positions(4);
    positions[0] = Vec3(0, 0, 0);
    positions[1] = Vec3(1, 0, 0);
    positions[2] = Vec3(0, 1, 0);
    context.setPositions(positions);
    context.applyConstraints(0.0001);
    for (int i = 0; i < 1000; i++) {
        State state = context.getState(State::Positions | State::Forces);
        const vector<Vec3>& pos = state.getPositions();
        Vec3 v12 = pos[1]-pos[0];
        Vec3 v13 = pos[2]-pos[0];
        Vec3 cross = v12.cross(v13);
        ASSERT_EQUAL_VEC(pos[0]+v12*0.3+v13*0.4+cross*0.5, pos[3], 1e-10);
        const vector<Vec3>& f = state.getForces();
        ASSERT_EQUAL_VEC(Vec3(0.1+0.2+0.3+0.4, 0, 0), f[0]+f[1]+f[2], 1e-10);
        Vec3 f2(0.4*0.3, 0.4*0.5*v13[2], -0.4*0.5*v13[1]);
        Vec3 f3(0.4*0.4, -0.4*0.5*v12[2], 0.4*0.5*v12[1]);
        ASSERT_EQUAL_VEC(Vec3(0.1+0.4, 0, 0)-f2-f3, f[0], 1e-10);
        ASSERT_EQUAL_VEC(Vec3(0.2, 0, 0)+f2, f[1], 1e-10);
        ASSERT_EQUAL_VEC(Vec3(0.3, 0, 0)+f3, f[2], 1e-10);
    }
}

/**
 * Make sure that energy, linear momentum, and angular momentum are all conserved
 * when using virtual sites.
 */
void testConservationLaws() {
    ReferencePlatform platform;
    System system;
    NonbondedForce* forceField = new NonbondedForce();
    system.addForce(forceField);
    vector<Vec3> positions;
    
    // Create a linear molecule with a TwoParticleAverage virtual site.
    
    system.addParticle(1.0);
    system.addParticle(1.0);
    system.addParticle(0.0);
    system.setVirtualSite(2, new VirtualSite::TwoParticleAverage(0, 1, 0.4, 0.6));
    system.addConstraint(0, 1, 2.0);
    for (int i = 0; i < 3; i++) {
        forceField->addParticle(0, 1, 10);
        for (int j = 0; j < i; j++)
            forceField->addException(i, j, 0, 1, 0);
    }
    positions.push_back(Vec3(0, 0, 0));
    positions.push_back(Vec3(2, 0, 0));
    positions.push_back(Vec3());
    
    // Create a planar molecule with a ThreeParticleAverage virtual site.
    
    system.addParticle(1.0);
    system.addParticle(1.0);
    system.addParticle(1.0);
    system.addParticle(0.0);
    system.setVirtualSite(6, new VirtualSite::ThreeParticleAverage(3, 4, 5, 0.3, 0.5, 0.2));
    system.addConstraint(3, 4, 1.0);
    system.addConstraint(3, 5, 1.0);
    system.addConstraint(4, 5, sqrt(2.0));
    for (int i = 0; i < 4; i++) {
        forceField->addParticle(0, 1, 10);
        for (int j = 0; j < i; j++)
            forceField->addException(i+3, j+3, 0, 1, 0);
    }
    positions.push_back(Vec3(0, 0, 1));
    positions.push_back(Vec3(1, 0, 1));
    positions.push_back(Vec3(0, 1, 1));
    positions.push_back(Vec3());
    
    // Create a tetrahedral molecule with an OutOfPlane virtual site.
    
    system.addParticle(1.0);
    system.addParticle(1.0);
    system.addParticle(1.0);
    system.addParticle(0.0);
    system.setVirtualSite(10, new VirtualSite::OutOfPlane(7, 8, 9, 0.3, 0.5, 0.2));
    system.addConstraint(7, 8, 1.0);
    system.addConstraint(7, 9, 1.0);
    system.addConstraint(8, 9, sqrt(2.0));
    for (int i = 0; i < 4; i++) {
        forceField->addParticle(0, 1, 10);
        for (int j = 0; j < i; j++)
            forceField->addException(i+7, j+7, 0, 1, 0);
    }
    positions.push_back(Vec3(1, 0, -1));
    positions.push_back(Vec3(2, 0, -1));
    positions.push_back(Vec3(1, 1, -1));
    positions.push_back(Vec3());

    // Simulate it and check conservation laws.
    
    VerletIntegrator integrator(0.002);
    Context context(system, integrator, platform);
    context.setPositions(positions);
    context.applyConstraints(0.0001);
    int numParticles = system.getNumParticles();
    double initialEnergy;
    Vec3 initialMomentum, initialAngularMomentum;
    for (int i = 0; i < 1000; i++) {
        State state = context.getState(State::Positions | State::Velocities | State::Forces | State::Energy);
        const vector<Vec3>& pos = state.getPositions();
        const vector<Vec3>& vel = state.getVelocities();
        const vector<Vec3>& f = state.getForces();
        double energy = state.getPotentialEnergy();
        for (int i = 0; i < numParticles; i++) {
            Vec3 v = vel[i] + f[i]*0.5*integrator.getStepSize();
            energy += 0.5*system.getParticleMass(i)*v.dot(v);
        }
        if (i == 0)
            initialEnergy = energy;
        else
            ASSERT_EQUAL_TOL(initialEnergy, energy, 0.01);
        Vec3 momentum;
        for (int j = 0; j < numParticles; j++)
            momentum += vel[j]*system.getParticleMass(j);
        if (i == 0)
            initialMomentum = momentum;
        else
            ASSERT_EQUAL_VEC(initialMomentum, momentum, 1e-10);
        Vec3 angularMomentum;
        for (int j = 0; j < numParticles; j++)
            angularMomentum += pos[j].cross(vel[j])*system.getParticleMass(j);
        if (i == 0)
            initialAngularMomentum = angularMomentum;
        else
            ASSERT_EQUAL_VEC(initialAngularMomentum, angularMomentum, 1e-10);
        integrator.step(1);
    }
}

int main() {
    try {
        testMasslessParticle();
        testTwoParticleAverage();
        testThreeParticleAverage();
        testOutOfPlane();
        testConservationLaws();
    }
    catch(const exception& e) {
        cout << "exception: " << e.what() << endl;
        return 1;
    }
    cout << "Done" << endl;
    return 0;
}
