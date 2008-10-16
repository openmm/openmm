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
 * This tests all the different force terms in the reference implementation of StandardMMForceField.
 */

#include "../../../tests/AssertionUtilities.h"
#include "OpenMMContext.h"
#include "CudaPlatform.h"
#include "StandardMMForceField.h"
#include "System.h"
#include "LangevinIntegrator.h"
#include "../src/SimTKUtilities/SimTKOpenMMRealType.h"
#include <iostream>
#include <vector>

using namespace OpenMM;
using namespace std;

const double TOL = 1e-5;

void testBonds() {
    CudaPlatform platform;
    System system(3, 0);
    LangevinIntegrator integrator(0.0, 0.1, 0.01);
    StandardMMForceField* forceField = new StandardMMForceField(3, 2, 0, 0, 0, 0);
    forceField->setBondParameters(0, 0, 1, 1.5, 0.8);
    forceField->setBondParameters(1, 1, 2, 1.2, 0.7);
    system.addForce(forceField);
    OpenMMContext context(system, integrator, platform);
    vector<Vec3> positions(3);
    positions[0] = Vec3(0, 2, 0);
    positions[1] = Vec3(0, 0, 0);
    positions[2] = Vec3(1, 0, 0);
    context.setPositions(positions);
    State state = context.getState(State::Forces | State::Energy);
    const vector<Vec3>& forces = state.getForces();
    ASSERT_EQUAL_VEC(Vec3(0, -0.8*0.5, 0), forces[0], TOL);
    ASSERT_EQUAL_VEC(Vec3(0.7*0.2, 0, 0), forces[2], TOL);
    ASSERT_EQUAL_VEC(Vec3(-forces[0][0]-forces[2][0], -forces[0][1]-forces[2][1], -forces[0][2]-forces[2][2]), forces[1], TOL);
    ASSERT_EQUAL_TOL(0.5*0.8*0.5*0.5 + 0.5*0.7*0.2*0.2, state.getPotentialEnergy(), TOL);
}

void testAngles() {
    CudaPlatform platform;
    System system(4, 0);
    LangevinIntegrator integrator(0.0, 0.1, 0.01);
    StandardMMForceField* forceField = new StandardMMForceField(4, 0, 2, 0, 0, 0);
    forceField->setAngleParameters(0, 0, 1, 2, PI_M/3, 1.1);
    forceField->setAngleParameters(1, 1, 2, 3, PI_M/2, 1.2);
    system.addForce(forceField);
    OpenMMContext context(system, integrator, platform);
    vector<Vec3> positions(4);
    positions[0] = Vec3(0, 1, 0);
    positions[1] = Vec3(0, 0, 0);
    positions[2] = Vec3(1, 0, 0);
    positions[3] = Vec3(2, 1, 0);
    context.setPositions(positions);
    State state = context.getState(State::Forces | State::Energy);
    const vector<Vec3>& forces = state.getForces();
    double torque1 = 1.1*PI_M/6;
    double torque2 = 1.2*PI_M/4;
    ASSERT_EQUAL_VEC(Vec3(torque1, 0, 0), forces[0], TOL);
    ASSERT_EQUAL_VEC(Vec3(-0.5*torque2, 0.5*torque2, 0), forces[3], TOL); // reduced by sqrt(2) due to the bond length, another sqrt(2) due to the angle
    ASSERT_EQUAL_VEC(Vec3(forces[0][0]+forces[1][0]+forces[2][0]+forces[3][0], forces[0][1]+forces[1][1]+forces[2][1]+forces[3][1], forces[0][2]+forces[1][2]+forces[2][2]+forces[3][2]), Vec3(0, 0, 0), TOL);
    ASSERT_EQUAL_TOL(0.5*1.1*(PI_M/6)*(PI_M/6) + 0.5*1.2*(PI_M/4)*(PI_M/4), state.getPotentialEnergy(), TOL);
}

void testPeriodicTorsions() {
    CudaPlatform platform;
    System system(4, 0);
    LangevinIntegrator integrator(0.0, 0.1, 0.01);
    StandardMMForceField* forceField = new StandardMMForceField(4, 0, 0, 1, 0, 0);
    forceField->setPeriodicTorsionParameters(0, 0, 1, 2, 3, 2, PI_M/3, 1.1);
    system.addForce(forceField);
    OpenMMContext context(system, integrator, platform);
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
    ASSERT_EQUAL_VEC(Vec3(0, 0.5*torque, 0), forces[3], TOL);
    ASSERT_EQUAL_VEC(Vec3(forces[0][0]+forces[1][0]+forces[2][0]+forces[3][0], forces[0][1]+forces[1][1]+forces[2][1]+forces[3][1], forces[0][2]+forces[1][2]+forces[2][2]+forces[3][2]), Vec3(0, 0, 0), TOL);
    ASSERT_EQUAL_TOL(1.1*(1+std::cos(2*PI_M/3)), state.getPotentialEnergy(), TOL);
}

void testRBTorsions() {
    CudaPlatform platform;
    System system(4, 0);
    LangevinIntegrator integrator(0.0, 0.1, 0.01);
    StandardMMForceField* forceField = new StandardMMForceField(4, 0, 0, 0, 1, 0);
    forceField->setRBTorsionParameters(0, 0, 1, 2, 3, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6);
    system.addForce(forceField);
    OpenMMContext context(system, integrator, platform);
    vector<Vec3> positions(4);
    positions[0] = Vec3(0, 1, 0);
    positions[1] = Vec3(0, 0, 0);
    positions[2] = Vec3(1, 0, 0);
    positions[3] = Vec3(1, 1, 1);
    context.setPositions(positions);
    State state = context.getState(State::Forces | State::Energy);
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

void testCoulomb() {
    CudaPlatform platform;
    System system(2, 0);
    LangevinIntegrator integrator(0.0, 0.1, 0.01);
    StandardMMForceField* forceField = new StandardMMForceField(2, 0, 0, 0, 0, 0);
    forceField->setAtomParameters(0, 0.5, 1, 0);
    forceField->setAtomParameters(1, -1.5, 1, 0);
    system.addForce(forceField);
    OpenMMContext context(system, integrator, platform);
    vector<Vec3> positions(2);
    positions[0] = Vec3(0, 0, 0);
    positions[1] = Vec3(2, 0, 0);
    context.setPositions(positions);
    State state = context.getState(State::Forces | State::Energy);
    const vector<Vec3>& forces = state.getForces();
    double force = 138.935485*(-0.75)/4.0;
    ASSERT_EQUAL_VEC(Vec3(-force, 0, 0), forces[0], TOL);
    ASSERT_EQUAL_VEC(Vec3(force, 0, 0), forces[1], TOL);
    ASSERT_EQUAL_TOL(138.935485*(-0.75)/2.0, state.getPotentialEnergy(), TOL);
}

void testLJ() {
    CudaPlatform platform;
    System system(2, 0);
    LangevinIntegrator integrator(0.0, 0.1, 0.01);
    StandardMMForceField* forceField = new StandardMMForceField(2, 0, 0, 0, 0, 0);
    forceField->setAtomParameters(0, 0, 1.2, 1);
    forceField->setAtomParameters(1, 0, 1.4, 2);
    system.addForce(forceField);
    OpenMMContext context(system, integrator, platform);
    vector<Vec3> positions(2);
    positions[0] = Vec3(0, 0, 0);
    positions[1] = Vec3(2, 0, 0);
    context.setPositions(positions);
    State state = context.getState(State::Forces | State::Energy);
    const vector<Vec3>& forces = state.getForces();
    double x = 1.3/2.0;
    double eps = SQRT_TWO;
    double force = 4.0*eps*(12*std::pow(x, 12.0)-6*std::pow(x, 6.0))/2.0;
    ASSERT_EQUAL_VEC(Vec3(-force, 0, 0), forces[0], TOL);
    ASSERT_EQUAL_VEC(Vec3(force, 0, 0), forces[1], TOL);
    ASSERT_EQUAL_TOL(4.0*eps*(std::pow(x, 12.0)-std::pow(x, 6.0)), state.getPotentialEnergy(), TOL);
}

void testExclusionsAnd14() {
    CudaPlatform platform;
    System system(5, 0);
    LangevinIntegrator integrator(0.0, 0.1, 0.01);
    StandardMMForceField* forceField = new StandardMMForceField(5, 4, 0, 0, 0, 2);
    forceField->setBondParameters(0, 0, 1, 1, 0);
    forceField->setBondParameters(1, 1, 2, 1, 0);
    forceField->setBondParameters(2, 2, 3, 1, 0);
    forceField->setBondParameters(3, 3, 4, 1, 0);
    system.addForce(forceField);
    for (int i = 1; i < 5; ++i) {
 
        // Test LJ forces
        
        vector<Vec3> positions(5);
        const double r = 1.0;
        for (int j = 0; j < 5; ++j) {
            forceField->setAtomParameters(j, 0, 1.5, 0);
            positions[j] = Vec3(0, j, 0);
        }
        forceField->setAtomParameters(0, 0, 1.5, 1);
        forceField->setAtomParameters(i, 0, 1.5, 1);
        forceField->setNonbonded14Parameters(0, 0, 3, 0, 1.5, i == 3 ? 0.5 : 0.0);
        forceField->setNonbonded14Parameters(1, 1, 4, 0, 1.5, 0.0);
        positions[i] = Vec3(r, 0, 0);
        OpenMMContext context(system, integrator, platform);
        context.setPositions(positions);
        State state = context.getState(State::Forces | State::Energy);
        const vector<Vec3>& forces = state.getForces();
        double x = 1.5/r;
        double eps = 1.0;
        double force = 4.0*eps*(12*std::pow(x, 12.0)-6*std::pow(x, 6.0))/r;
        double energy = 4.0*eps*(std::pow(x, 12.0)-std::pow(x, 6.0));
        if (i == 3) {
            force *= 0.5;
            energy *= 0.5;
        }
        if (i < 3) {
            force = 0;
            energy = 0;
        }
        ASSERT_EQUAL_VEC(Vec3(-force, 0, 0), forces[0], TOL);
        ASSERT_EQUAL_VEC(Vec3(force, 0, 0), forces[i], TOL);
        ASSERT_EQUAL_TOL(energy, state.getPotentialEnergy(), TOL);

        // Test Coulomb forces
        
        forceField->setAtomParameters(0, 2, 1.5, 0);
        forceField->setAtomParameters(i, 2, 1.5, 0);
        forceField->setNonbonded14Parameters(0, 0, 3, i == 3 ? 4/1.2 : 0, 1.5, 0);
        forceField->setNonbonded14Parameters(1, 1, 4, 0, 1.5, 0);
        OpenMMContext context2(system, integrator, platform);
        context2.setPositions(positions);
        state = context2.getState(State::Forces | State::Energy);
        const vector<Vec3>& forces2 = state.getForces();
        force = 138.935485*4/(r*r);
        energy = 138.935485*4/r;
        if (i == 3) {
            force /= 1.2;
            energy /= 1.2;
        }
        if (i < 3) {
            force = 0;
            energy = 0;
        }
        ASSERT_EQUAL_VEC(Vec3(-force, 0, 0), forces2[0], TOL);
        ASSERT_EQUAL_VEC(Vec3(force, 0, 0), forces2[i], TOL);
        ASSERT_EQUAL_TOL(energy, state.getPotentialEnergy(), TOL);
    }
}

void testCutoff() {
    CudaPlatform platform;
    System system(3, 0);
    LangevinIntegrator integrator(0.0, 0.1, 0.01);
    StandardMMForceField* forceField = new StandardMMForceField(3, 0, 0, 0, 0, 0);
    forceField->setAtomParameters(0, 1.0, 1, 0);
    forceField->setAtomParameters(1, 1.0, 1, 0);
    forceField->setAtomParameters(2, 1.0, 1, 0);
    forceField->setNonbondedMethod(StandardMMForceField::CutoffNonPeriodic);
    const double cutoff = 2.9;
    forceField->setCutoffDistance(cutoff);
    system.addForce(forceField);
    OpenMMContext context(system, integrator, platform);
    vector<Vec3> positions(3);
    positions[0] = Vec3(0, 0, 0);
    positions[1] = Vec3(0, 2, 0);
    positions[2] = Vec3(0, 3, 0);
    context.setPositions(positions);
    State state = context.getState(State::Forces | State::Energy);
    const vector<Vec3>& forces = state.getForces();
    const double eps = 78.3;
    const double krf = (1.0/(cutoff*cutoff*cutoff))*(eps-1.0)/(2.0*eps+1.0);
    const double crf = (1.0/cutoff)*(3.0*eps)/(2.0*eps+1.0);
    const double force1 = 138.935485*(1.0)*(0.25-2.0*krf*2.0);
    const double force2 = 138.935485*(1.0)*(1.0-2.0*krf*1.0);
    ASSERT_EQUAL_VEC(Vec3(0, -force1, 0), forces[0], TOL);
    ASSERT_EQUAL_VEC(Vec3(0, force1-force2, 0), forces[1], TOL);
    ASSERT_EQUAL_VEC(Vec3(0, force2, 0), forces[2], TOL);
    const double energy1 = 138.935485*(1.0)*(0.5+krf*4.0-crf);
    const double energy2 = 138.935485*(1.0)*(1.0+krf*1.0-crf);
    ASSERT_EQUAL_TOL(energy1+energy2, state.getPotentialEnergy(), TOL);
}

void testCutoff14() {
    CudaPlatform platform;
    System system(5, 0);
    LangevinIntegrator integrator(0.0, 0.1, 0.01);
    StandardMMForceField* forceField = new StandardMMForceField(5, 4, 0, 0, 0, 2);
    forceField->setBondParameters(0, 0, 1, 1, 0);
    forceField->setBondParameters(1, 1, 2, 1, 0);
    forceField->setBondParameters(2, 2, 3, 1, 0);
    forceField->setBondParameters(3, 3, 4, 1, 0);
    forceField->setNonbondedMethod(StandardMMForceField::CutoffNonPeriodic);
    const double cutoff = 3.5;
    forceField->setCutoffDistance(cutoff);
    system.addForce(forceField);
    OpenMMContext context(system, integrator, platform);
    vector<Vec3> positions(5);
    positions[0] = Vec3(0, 0, 0);
    positions[1] = Vec3(1, 0, 0);
    positions[2] = Vec3(2, 0, 0);
    positions[3] = Vec3(3, 0, 0);
    positions[4] = Vec3(4, 0, 0);
    for (int i = 1; i < 5; ++i) {
 
        // Test LJ forces
        
        forceField->setAtomParameters(0, 0, 1.5, 1);
        for (int j = 1; j < 5; ++j)
            forceField->setAtomParameters(j, 0, 1.5, 0);
        forceField->setAtomParameters(i, 0, 1.5, 1);
        forceField->setNonbonded14Parameters(0, 0, 3, 0, 1.5, i == 3 ? 0.5 : 0.0);
        forceField->setNonbonded14Parameters(1, 1, 4, 0, 1.5, 0.0);
        context.reinitialize();
        context.setPositions(positions);
        State state = context.getState(State::Forces | State::Energy);
        const vector<Vec3>& forces = state.getForces();
        double r = positions[i][0];
        double x = 1.5/r;
        double e = 1.0;
        double force = 4.0*e*(12*std::pow(x, 12.0)-6*std::pow(x, 6.0))/r;
        double energy = 4.0*e*(std::pow(x, 12.0)-std::pow(x, 6.0));
        if (i == 3) {
            force *= 0.5;
            energy *= 0.5;
        }
        if (i < 3 || r > cutoff) {
            force = 0;
            energy = 0;
        }
        ASSERT_EQUAL_VEC(Vec3(-force, 0, 0), forces[0], TOL);
        ASSERT_EQUAL_VEC(Vec3(force, 0, 0), forces[i], TOL);
        ASSERT_EQUAL_TOL(energy, state.getPotentialEnergy(), TOL);

        // Test Coulomb forces
        
        const double q = 0.7;
        forceField->setAtomParameters(0, q, 1.5, 0);
        forceField->setAtomParameters(i, q, 1.5, 0);
        forceField->setNonbonded14Parameters(0, 0, 3, i == 3 ? q*q/1.2 : 0, 1.5, 0);
        forceField->setNonbonded14Parameters(1, 1, 4, 0, 1.5, 0);
        context.reinitialize();
        context.setPositions(positions);
        state = context.getState(State::Forces | State::Energy);
        const vector<Vec3>& forces2 = state.getForces();
        const double eps = 78.3;
        const double krf = (1.0/(cutoff*cutoff*cutoff))*(eps-1.0)/(2.0*eps+1.0);
        const double crf = (1.0/cutoff)*(3.0*eps)/(2.0*eps+1.0);
        force = 138.935485*q*q*(1.0/(r*r)-2.0*krf*r);
        energy = 138.935485*q*q*(1.0/r+krf*r*r-crf);
        if (i == 3) {
            force /= 1.2;
            energy /= 1.2;
        }
        if (i < 3 || r > cutoff) {
            force = 0;
            energy = 0;
        }
        ASSERT_EQUAL_VEC(Vec3(-force, 0, 0), forces2[0], TOL);
        ASSERT_EQUAL_VEC(Vec3(force, 0, 0), forces2[i], TOL);
        ASSERT_EQUAL_TOL(energy, state.getPotentialEnergy(), TOL);
    }
}

void testPeriodic() {
    CudaPlatform platform;
    System system(3, 0);
    LangevinIntegrator integrator(0.0, 0.1, 0.01);
    StandardMMForceField* forceField = new StandardMMForceField(3, 1, 0, 0, 0, 0);
    forceField->setAtomParameters(0, 1.0, 1, 0);
    forceField->setAtomParameters(1, 1.0, 1, 0);
    forceField->setAtomParameters(2, 1.0, 1, 0);
    forceField->setBondParameters(0, 0, 1, 1.0, 0.0);
    forceField->setNonbondedMethod(StandardMMForceField::CutoffPeriodic);
    const double cutoff = 2.0;
    forceField->setCutoffDistance(cutoff);
    forceField->setPeriodicBoxSize(4.0, 4.0, 4.0);
    system.addForce(forceField);
    OpenMMContext context(system, integrator, platform);
    vector<Vec3> positions(3);
    positions[0] = Vec3(0, 0, 0);
    positions[1] = Vec3(2, 0, 0);
    positions[2] = Vec3(3, 0, 0);
    context.setPositions(positions);
    State state = context.getState(State::Forces | State::Energy);
    const vector<Vec3>& forces = state.getForces();
    const double eps = 78.3;
    const double krf = (1.0/(cutoff*cutoff*cutoff))*(eps-1.0)/(2.0*eps+1.0);
    const double crf = (1.0/cutoff)*(3.0*eps)/(2.0*eps+1.0);
    const double force = 138.935485*(1.0)*(1.0-2.0*krf*1.0);
    ASSERT_EQUAL_VEC(Vec3(force, 0, 0), forces[0], TOL);
    ASSERT_EQUAL_VEC(Vec3(-force, 0, 0), forces[1], TOL);
    ASSERT_EQUAL_VEC(Vec3(0, 0, 0), forces[2], TOL);
    ASSERT_EQUAL_TOL(2*138.935485*(1.0)*(1.0+krf*1.0-crf), state.getPotentialEnergy(), TOL);
}

int main() {
    try {
        testBonds();
        testAngles();
        testPeriodicTorsions();
        testRBTorsions();
        testCoulomb();
        testLJ();
        testExclusionsAnd14();
//        testCutoff();
//        testCutoff14();
//        testPeriodic();
    }
    catch(const exception& e) {
        cout << "exception: " << e.what() << endl;
        return 1;
    }
    cout << "Done" << endl;
    return 0;
}
