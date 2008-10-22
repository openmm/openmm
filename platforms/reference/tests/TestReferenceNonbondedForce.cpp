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
 * This tests all the different force terms in the reference implementation of NonbondedForce.
 */

#include "../../../tests/AssertionUtilities.h"
#include "OpenMMContext.h"
#include "ReferencePlatform.h"
#include "NonbondedForce.h"
#include "System.h"
#include "VerletIntegrator.h"
#include "../src/SimTKUtilities/SimTKOpenMMRealType.h"
#include "HarmonicBondForce.h"
#include <iostream>
#include <vector>

using namespace OpenMM;
using namespace std;

const double TOL = 1e-5;

void testCoulomb() {
    ReferencePlatform platform;
    System system(2, 0);
    VerletIntegrator integrator(0.01);
    NonbondedForce* forceField = new NonbondedForce(2, 0);
    forceField->setParticleParameters(0, 0.5, 1, 0);
    forceField->setParticleParameters(1, -1.5, 1, 0);
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
    ReferencePlatform platform;
    System system(2, 0);
    VerletIntegrator integrator(0.01);
    NonbondedForce* forceField = new NonbondedForce(2, 0);
    forceField->setParticleParameters(0, 0, 1.2, 1);
    forceField->setParticleParameters(1, 0, 1.4, 2);
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
    ReferencePlatform platform;
    System system(5, 0);
    VerletIntegrator integrator(0.01);
    HarmonicBondForce* bonds = new HarmonicBondForce(4);
    bonds->setBondParameters(0, 0, 1, 1, 0);
    bonds->setBondParameters(1, 1, 2, 1, 0);
    bonds->setBondParameters(2, 2, 3, 1, 0);
    bonds->setBondParameters(3, 3, 4, 1, 0);
    system.addForce(bonds);
    NonbondedForce* nonbonded = new NonbondedForce(5, 2);
    system.addForce(nonbonded);
    OpenMMContext context(system, integrator, platform);
    for (int i = 1; i < 5; ++i) {
 
        // Test LJ forces
        
        vector<Vec3> positions(5);
        const double r = 1.0;
        for (int j = 0; j < 5; ++j) {
            nonbonded->setParticleParameters(j, 0, 1.5, 0);
            positions[j] = Vec3(0, j, 0);
        }
        nonbonded->setParticleParameters(0, 0, 1.5, 1);
        nonbonded->setParticleParameters(i, 0, 1.5, 1);
        nonbonded->setNonbonded14Parameters(0, 0, 3, 0, 1.5, i == 3 ? 0.5 : 0.0);
        nonbonded->setNonbonded14Parameters(1, 1, 4, 0, 1.5, 0.0);
        positions[i] = Vec3(r, 0, 0);
        context.reinitialize();
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
        
        nonbonded->setParticleParameters(0, 2, 1.5, 0);
        nonbonded->setParticleParameters(i, 2, 1.5, 0);
        nonbonded->setNonbonded14Parameters(0, 0, 3, i == 3 ? 4/1.2 : 0, 1.5, 0);
        nonbonded->setNonbonded14Parameters(1, 1, 4, 0, 1.5, 0);
        context.reinitialize();
        context.setPositions(positions);
        state = context.getState(State::Forces | State::Energy);
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
    ReferencePlatform platform;
    System system(3, 0);
    VerletIntegrator integrator(0.01);
    NonbondedForce* forceField = new NonbondedForce(3, 0);
    forceField->setParticleParameters(0, 1.0, 1, 0);
    forceField->setParticleParameters(1, 1.0, 1, 0);
    forceField->setParticleParameters(2, 1.0, 1, 0);
    forceField->setNonbondedMethod(NonbondedForce::CutoffNonPeriodic);
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
    ReferencePlatform platform;
    System system(5, 0);
    VerletIntegrator integrator(0.01);
    HarmonicBondForce* bonds = new HarmonicBondForce(4);
    bonds->setBondParameters(0, 0, 1, 1, 0);
    bonds->setBondParameters(1, 1, 2, 1, 0);
    bonds->setBondParameters(2, 2, 3, 1, 0);
    bonds->setBondParameters(3, 3, 4, 1, 0);
    system.addForce(bonds);
    NonbondedForce* nonbonded = new NonbondedForce(5, 2);
    nonbonded->setNonbondedMethod(NonbondedForce::CutoffNonPeriodic);
    const double cutoff = 3.5;
    nonbonded->setCutoffDistance(cutoff);
    system.addForce(nonbonded);
    OpenMMContext context(system, integrator, platform);
    vector<Vec3> positions(5);
    positions[0] = Vec3(0, 0, 0);
    positions[1] = Vec3(1, 0, 0);
    positions[2] = Vec3(2, 0, 0);
    positions[3] = Vec3(3, 0, 0);
    positions[4] = Vec3(4, 0, 0);
    for (int i = 1; i < 5; ++i) {
 
        // Test LJ forces
        
        nonbonded->setParticleParameters(0, 0, 1.5, 1);
        for (int j = 1; j < 5; ++j)
            nonbonded->setParticleParameters(j, 0, 1.5, 0);
        nonbonded->setParticleParameters(i, 0, 1.5, 1);
        nonbonded->setNonbonded14Parameters(0, 0, 3, 0, 1.5, i == 3 ? 0.5 : 0.0);
        nonbonded->setNonbonded14Parameters(1, 1, 4, 0, 1.5, 0.0);
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
        nonbonded->setParticleParameters(0, q, 1.5, 0);
        nonbonded->setParticleParameters(i, q, 1.5, 0);
        nonbonded->setNonbonded14Parameters(0, 0, 3, i == 3 ? q*q/1.2 : 0, 1.5, 0);
        nonbonded->setNonbonded14Parameters(1, 1, 4, 0, 1.5, 0);
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
    ReferencePlatform platform;
    System system(3, 0);
    VerletIntegrator integrator(0.01);
    HarmonicBondForce* bonds = new HarmonicBondForce(1);
    bonds->setBondParameters(0, 0, 1, 1, 0);
    system.addForce(bonds);
    NonbondedForce* nonbonded = new NonbondedForce(3, 0);
    nonbonded->setParticleParameters(0, 1.0, 1, 0);
    nonbonded->setParticleParameters(1, 1.0, 1, 0);
    nonbonded->setParticleParameters(2, 1.0, 1, 0);
    nonbonded->setNonbondedMethod(NonbondedForce::CutoffPeriodic);
    const double cutoff = 2.0;
    nonbonded->setCutoffDistance(cutoff);
    nonbonded->setPeriodicBoxVectors(Vec3(4, 0, 0), Vec3(0, 4, 0), Vec3(0, 0, 4));
    system.addForce(nonbonded);
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
        testCoulomb();
        testLJ();
        testExclusionsAnd14();
        testCutoff();
        testCutoff14();
        testPeriodic();
    }
    catch(const exception& e) {
        cout << "exception: " << e.what() << endl;
        return 1;
    }
    cout << "Done" << endl;
    return 0;
}
