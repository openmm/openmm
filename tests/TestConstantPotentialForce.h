/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2008-2025 Stanford University and the Authors.      *
 * Authors: Peter Eastman, Evan Pretti                                        *
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
#include "ReferencePlatform.h"
#include "openmm/NonbondedForce.h"
#include "openmm/ConstantPotentialForce.h"
#include "openmm/System.h"
#include "openmm/VerletIntegrator.h"
#include "SimTKOpenMMRealType.h"
#include "sfmt/SFMT.h"
#include <iostream>
#include <iomanip>
#include <vector>

using namespace OpenMM;
using namespace std;

const double TOL = 1e-5;

void testCoulomb(bool periodicExceptions, bool triclinic) {
    // Ensures that the Coulomb energy and force computation for
    // ConstantPotentialForce without electrodes matches NonbondedForce.

    System refSystem;
    System testSystem;
    NonbondedForce* refForce = new NonbondedForce();
    ConstantPotentialForce* testForce = new ConstantPotentialForce();

    Vec3 a, b, c;
    if (triclinic) {
        a = Vec3(10, 0, 0);
        b = Vec3(1, 9, 0);
        c = Vec3(2, 3, 8);
    } else {
        a = Vec3(10, 0, 0);
        b = Vec3(0, 10, 0);
        c = Vec3(0, 0, 10);
    }
    refSystem.setDefaultPeriodicBoxVectors(a, b, c);
    testSystem.setDefaultPeriodicBoxVectors(a, b, c);

    vector<Vec3> positions;
    vector<pair<int, int>> bonds;
    for (int i = 0; i < 10; i++) {
        refSystem.addParticle(1);
        testSystem.addParticle(1);
        refForce->addParticle(i % 2 ? -1 : 1, 1, 0);
        testForce->addParticle(i % 2 ? -1 : 1);
        double f = 0.1 * i;
        positions.push_back(f * a + f * f * b + f * f * f * c);
        bonds.push_back(pair<int, int>(i, (i + 1) % 10));
    }
    refForce->createExceptionsFromBonds(bonds, 0.5, 0.0);
    testForce->createExceptionsFromBonds(bonds, 0.5);
    refForce->setExceptionsUsePeriodicBoundaryConditions(periodicExceptions);
    testForce->setExceptionsUsePeriodicBoundaryConditions(periodicExceptions);

    refForce->setNonbondedMethod(NonbondedForce::PME);

    refSystem.addForce(refForce);
    testSystem.addForce(testForce);

    ASSERT(testForce->usesPeriodicBoundaryConditions());
    ASSERT(testSystem.usesPeriodicBoundaryConditions());

    VerletIntegrator refIntegrator(0.001);
    VerletIntegrator testIntegrator(0.001);
    Context refContext(refSystem, refIntegrator, platform);
    Context testContext(testSystem, testIntegrator, platform);
    refContext.setPositions(positions);
    testContext.setPositions(positions);
    State refState = refContext.getState(State::Energy | State::Forces);
    State testState = testContext.getState(State::Energy | State::Forces);
    const vector<Vec3>& refForces = refState.getForces();
    const vector<Vec3>& testForces = testState.getForces();

    ASSERT_EQUAL_TOL(refState.getPotentialEnergy(), testState.getPotentialEnergy(), TOL);
    for (int i = 0; i < refForces.size(); i++) {
        ASSERT_EQUAL_VEC(refForces[i], testForces[i], TOL);
    }
}

void testCoulombOverlap() {
    // Ensures that two particles with an exception can overlap without
    // generating NaN energies and forces.

    System refSystem;
    System testSystem;
    NonbondedForce* refForce = new NonbondedForce();
    ConstantPotentialForce* testForce = new ConstantPotentialForce();

    refSystem.setDefaultPeriodicBoxVectors(Vec3(10, 0, 0), Vec3(0, 10, 0), Vec3(0, 0, 10));
    testSystem.setDefaultPeriodicBoxVectors(Vec3(10, 0, 0), Vec3(0, 10, 0), Vec3(0, 0, 10));

    for (int i = 0; i < 3; i++) {
        refSystem.addParticle(1);
        testSystem.addParticle(1);
    }
    refForce->addParticle(1, 1, 0);
    testForce->addParticle(1);
    refForce->addParticle(1, 1, 0);
    testForce->addParticle(1);
    refForce->addParticle(-2, 1, 0);
    testForce->addParticle(-2);
    refForce->addException(0, 1, 0, 1, 0);
    testForce->addException(0, 1, 0);

    refForce->setNonbondedMethod(NonbondedForce::PME);

    refSystem.addForce(refForce);
    testSystem.addForce(testForce);

    vector<Vec3> positions{
        Vec3(0, 0, 0),
        Vec3(0, 0, 0),
        Vec3(1, 2, 3)
    };

    VerletIntegrator refIntegrator(0.001);
    VerletIntegrator testIntegrator(0.001);
    Context refContext(refSystem, refIntegrator, platform);
    Context testContext(testSystem, testIntegrator, platform);
    refContext.setPositions(positions);
    testContext.setPositions(positions);
    State refState = refContext.getState(State::Energy | State::Forces);
    State testState = testContext.getState(State::Energy | State::Forces);
    const vector<Vec3>& refForces = refState.getForces();
    const vector<Vec3>& testForces = testState.getForces();

    ASSERT_EQUAL_TOL(refState.getPotentialEnergy(), testState.getPotentialEnergy(), TOL);
    for (int i = 0; i < refForces.size(); i++) {
        ASSERT_EQUAL_VEC(refForces[i], testForces[i], TOL);
    }
}

void testElectrodeNoOverlaps() {
    // Ensures that a particle cannot belong to more than one electrode.

    System system;
    ConstantPotentialForce* force = new ConstantPotentialForce();
    for (int i = 0; i < 3; i++) {
        system.addParticle(1);
        force->addParticle(0);
    }
    std::set<int> electrode1{0, 1};
    std::set<int> electrode2{1, 2};
    force->addElectrode(electrode1, 0, 0, 0);
    force->addElectrode(electrode2, 0, 0, 0);
    system.addForce(force);
    VerletIntegrator integrator(0.001);

    bool thrown = false;
    try {
        Context context(system, integrator, platform);
    }
    catch (const OpenMMException& exception) {
        thrown = true;
    }
    ASSERT(thrown);
}

void testElectrodeNoExceptions() {
    // Ensures that a particle cannot belong to both an electrode and an
    // exception.

    System system;
    ConstantPotentialForce* force = new ConstantPotentialForce();
    for (int i = 0; i < 3; i++) {
        system.addParticle(1);
        force->addParticle(0);
    }
    std::set<int> electrode{0, 1};
    force->addElectrode(electrode, 0, 0, 0);
    force->addException(1, 2, 0);
    system.addForce(force);
    VerletIntegrator integrator(0.001);

    bool thrown = false;
    try {
        Context context(system, integrator, platform);
    }
    catch (const OpenMMException& exception) {
        thrown = true;
    }
    ASSERT(thrown);
}

void testElectrodeMatrixNoMass() {
    // Ensures that we can use the matrix method when electrode particles have
    // no mass, but that we cannot when an electrode particle has mass.

    System system1;
    ConstantPotentialForce* force1 = new ConstantPotentialForce();
    system1.addParticle(0);
    system1.addParticle(0);
    system1.addParticle(1);
    force1->addParticle(0);
    force1->addParticle(0);
    force1->addParticle(0);
    std::set<int> electrode1{0, 1};
    force1->addElectrode(electrode1, 0, 0, 0);
    force1->setConstantPotentialMethod(ConstantPotentialForce::Matrix);
    system1.addForce(force1);
    VerletIntegrator integrator1(0.001);
    Context context1(system1, integrator1, platform);

    System system2;
    ConstantPotentialForce* force2 = new ConstantPotentialForce();
    system2.addParticle(0);
    system2.addParticle(1);
    system2.addParticle(1);
    force2->addParticle(0);
    force2->addParticle(0);
    force2->addParticle(0);
    std::set<int> electrode2{0, 1};
    force2->addElectrode(electrode2, 0, 0, 0);
    force2->setConstantPotentialMethod(ConstantPotentialForce::Matrix);
    system2.addForce(force2);
    VerletIntegrator integrator2(0.001);

    bool thrown = false;
    try {
        Context context2(system2, integrator2, platform);
    }
    catch (const OpenMMException& exception) {
        thrown = true;
    }
    ASSERT(thrown);
}

void runPlatformTests();

int main(int argc, char* argv[]) {
    try {
        initializeTests(argc, argv);
        testCoulomb(false, false);
        testCoulomb(true, false);
        testCoulomb(false, true);
        testCoulombOverlap();
        testElectrodeNoOverlaps();
        testElectrodeNoExceptions();
        testElectrodeMatrixNoMass();
        runPlatformTests();
    }
    catch(const exception& e) {
        cout << "exception: " << e.what() << endl;
        return 1;
    }
    cout << "Done" << endl;
    return 0;
}
