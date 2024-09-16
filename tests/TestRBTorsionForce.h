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
#include "openmm/RBTorsionForce.h"
#include "openmm/System.h"
#include "openmm/VerletIntegrator.h"
#include "SimTKOpenMMRealType.h"
#include <iostream>
#include <vector>

using namespace OpenMM;
using namespace std;

const double TOL = 1e-5;

void testRBTorsions() {
    System system;
    system.addParticle(1.0);
    system.addParticle(1.0);
    system.addParticle(1.0);
    system.addParticle(1.0);
    VerletIntegrator integrator(0.01);
    RBTorsionForce* forceField = new RBTorsionForce();
    forceField->addTorsion(0, 1, 2, 3, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6);
    system.addForce(forceField);
    ASSERT(!forceField->usesPeriodicBoundaryConditions());
    ASSERT(!system.usesPeriodicBoundaryConditions());
    Context context(system, integrator, platform);
    vector<Vec3> positions(4);
    positions[0] = Vec3(0, 1, 0);
    positions[1] = Vec3(0, 0, 0);
    positions[2] = Vec3(1, 0, 0);
    positions[3] = Vec3(1, 1, 1);
    context.setPositions(positions);
    State state = context.getState(State::Forces | State::Energy);
    {
        const vector<Vec3>& forces = state.getForces();
        double psi = 0.25*PI_M - PI_M;
        double torque = 0.0;
        for (int i = 1; i < 6; ++i) {
            double c = 0.1*(i+1);
            torque += -c*i*std::pow(std::cos(psi), i-1)*std::sin(psi);
        }
        ASSERT_EQUAL_VEC(Vec3(0, 0, torque), forces[0], TOL);
        ASSERT_EQUAL_VEC(Vec3(0, 0.5*torque, -0.5*torque), forces[3], TOL);
        ASSERT_EQUAL_VEC(Vec3(forces[0][0]+forces[1][0]+forces[2][0]+forces[3][0], forces[0][1]+forces[1][1]+forces[2][1]+forces[3][1], forces[0][2]+forces[1][2]+forces[2][2]+forces[3][2]), Vec3(0, 0, 0), TOL);
        double energy = 0.0;
        for (int i = 0; i < 6; ++i) {
            double c = 0.1*(i+1);
            energy += c*std::pow(std::cos(psi), i);
        }
        ASSERT_EQUAL_TOL(energy, state.getPotentialEnergy(), TOL);
    }
    
    // Try changing the torsion parameters and make sure it's still correct.
    
    forceField->setTorsionParameters(0, 0, 1, 2, 3, 0.11, 0.22, 0.33, 0.44, 0.55, 0.66);
    forceField->updateParametersInContext(context);
    state = context.getState(State::Forces | State::Energy);
    {
        const vector<Vec3>& forces = state.getForces();
        double psi = 0.25*PI_M - PI_M;
        double torque = 0.0;
        for (int i = 1; i < 6; ++i) {
            double c = 0.11*(i+1);
            torque += -c*i*std::pow(std::cos(psi), i-1)*std::sin(psi);
        }
        ASSERT_EQUAL_VEC(Vec3(0, 0, torque), forces[0], TOL);
        ASSERT_EQUAL_VEC(Vec3(0, 0.5*torque, -0.5*torque), forces[3], TOL);
        ASSERT_EQUAL_VEC(Vec3(forces[0][0]+forces[1][0]+forces[2][0]+forces[3][0], forces[0][1]+forces[1][1]+forces[2][1]+forces[3][1], forces[0][2]+forces[1][2]+forces[2][2]+forces[3][2]), Vec3(0, 0, 0), TOL);
        double energy = 0.0;
        for (int i = 0; i < 6; ++i) {
            double c = 0.11*(i+1);
            energy += c*std::pow(std::cos(psi), i);
        }
        ASSERT_EQUAL_TOL(energy, state.getPotentialEnergy(), TOL);
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
    RBTorsionForce* torsions = new RBTorsionForce();
    torsions->addTorsion(0, 1, 2, 3, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6);
    torsions->setUsesPeriodicBoundaryConditions(true);
    system.addForce(torsions);
    Context context(system, integrator, platform);
    vector<Vec3> positions(4);
    positions[0] = Vec3(0, 1, 0);
    positions[1] = Vec3(0, 0, 0);
    positions[2] = Vec3(1, 0, 0);
    positions[3] = Vec3(0, 0, 2);
    context.setPositions(positions);
    State state = context.getState(State::Forces | State::Energy);
    const vector<Vec3>& forces = state.getForces();
    double psi = 0.5*PI_M;
    double torque = 0.0;
    for (int i = 1; i < 6; ++i) {
        double c = 0.1*(i+1);
        torque += -c*i*std::pow(std::cos(psi), i-1)*std::sin(psi);
    }
    ASSERT_EQUAL_VEC(Vec3(0, 0, torque), forces[0], TOL);
    ASSERT_EQUAL_VEC(Vec3(0, -torque, 0), forces[3], TOL);
    ASSERT_EQUAL_VEC(Vec3(forces[0][0]+forces[1][0]+forces[2][0]+forces[3][0], forces[0][1]+forces[1][1]+forces[2][1]+forces[3][1], forces[0][2]+forces[1][2]+forces[2][2]+forces[3][2]), Vec3(0, 0, 0), TOL);
    double energy = 0.0;
    for (int i = 0; i < 6; ++i) {
        double c = 0.1*(i+1);
        energy += c*std::pow(std::cos(psi), i);
    }
    ASSERT_EQUAL_TOL(energy, state.getPotentialEnergy(), TOL);
}

void testParallelComputation() {
    System system;
    const int numParticles = 200;
    for (int i = 0; i < numParticles; i++)
        system.addParticle(1.0);
    RBTorsionForce* force = new RBTorsionForce();
    for (int i = 3; i < numParticles; i++)
        force->addTorsion(i-3, i-2, i-1, i, 2, 0.1*i, 0.5*i, i, 1, 1);
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
}

void runPlatformTests();

int main(int argc, char* argv[]) {
    try {
        initializeTests(argc, argv);
        testRBTorsions();
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
