/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit.                   *
 * See https://openmm.org/development.                                        *
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
#include "openmm/HarmonicBondForce.h"
#include "openmm/NonbondedForce.h"
#include "openmm/CustomNonbondedForce.h"
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

const double TOL = 1e-4;

void testCoulomb(bool periodicExceptions) {
    // Ensures that the Coulomb energy and force computation for
    // ConstantPotentialForce without electrodes matches NonbondedForce.

    System refSystem;
    System testSystem;
    NonbondedForce* refForce = new NonbondedForce();
    ConstantPotentialForce* testForce = new ConstantPotentialForce();

    Vec3 a(10, 0, 0);
    Vec3 b(1, 9, 0);
    Vec3 c(2, 3, 8);
    refSystem.setDefaultPeriodicBoxVectors(a, b, c);
    testSystem.setDefaultPeriodicBoxVectors(a, b, c);

    vector<Vec3> positions;
    vector<pair<int, int>> bonds;
    for (int i = 0; i < 10; i++) {
        refSystem.addParticle(1);
        testSystem.addParticle(1);
        double testCharge = (i / 2 + 1) * (i % 2 ? -1 : 1);
        refForce->addParticle(testCharge, 1, 0);
        testForce->addParticle(testCharge);
        double f = 0.1 * i;
        positions.push_back(f * a + f * f * b + f * f * f * c);
        bonds.push_back(pair<int, int>(i, (i + 1) % 10));
    }
    refForce->createExceptionsFromBonds(bonds, 0.5, 0.0);
    testForce->createExceptionsFromBonds(bonds, 0.5);
    refForce->setExceptionsUsePeriodicBoundaryConditions(periodicExceptions);
    testForce->setExceptionsUsePeriodicBoundaryConditions(periodicExceptions);

    refForce->setNonbondedMethod(NonbondedForce::PME);
    refForce->setCutoffDistance(3.0);
    testForce->setCutoffDistance(3.0);

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

void testCoulombNonNeutral() {
    // Ensures that the energy of a non-neutral system at different cutoffs (and
    // thus different values of alpha) matches NonbondedForce.  Since the latter
    // includes a neutralizing plasma correction, this ensures that the same
    // correction is working in ConstantPotentialForce.

    System refSystem;
    System testSystem;
    NonbondedForce* refForce = new NonbondedForce();
    ConstantPotentialForce* testForce = new ConstantPotentialForce();

    Vec3 a, b, c;
    a = Vec3(10, 0, 0);
    b = Vec3(1, 9, 0);
    c = Vec3(2, 3, 8);
    refSystem.setDefaultPeriodicBoxVectors(a, b, c);
    testSystem.setDefaultPeriodicBoxVectors(a, b, c);

    vector<Vec3> positions;
    for (int i = 0; i < 10; i++) {
        refSystem.addParticle(1);
        testSystem.addParticle(1);
        double f = 0.1 * i;
        refForce->addParticle((i % 2 ? -1 : 1) + f * f, 1, 0);
        testForce->addParticle((i % 2 ? -1 : 1) + f * f);
        positions.push_back(f * a + f * f * b + f * f * f * c);
    }

    refForce->setNonbondedMethod(NonbondedForce::PME);

    refForce->setCutoffDistance(1.5);
    testForce->setCutoffDistance(1.5);

    refSystem.addForce(refForce);
    testSystem.addForce(testForce);

    VerletIntegrator refIntegrator(0.001);
    VerletIntegrator testIntegrator(0.001);
    Context refContext(refSystem, refIntegrator, platform);
    Context testContext(testSystem, testIntegrator, platform);
    refContext.setPositions(positions);
    testContext.setPositions(positions);

    {
        State refState = refContext.getState(State::Energy | State::Forces);
        State testState = testContext.getState(State::Energy | State::Forces);
        const vector<Vec3>& refForces = refState.getForces();
        const vector<Vec3>& testForces = testState.getForces();

        ASSERT_EQUAL_TOL(refState.getPotentialEnergy(), testState.getPotentialEnergy(), TOL);
        for (int i = 0; i < refForces.size(); i++) {
            ASSERT_EQUAL_VEC(refForces[i], testForces[i], TOL);
        }
    }

    refForce->setCutoffDistance(2.5);
    testForce->setCutoffDistance(2.5);

    refContext.reinitialize(true);
    testContext.reinitialize(true);

    {
        State refState = refContext.getState(State::Energy | State::Forces);
        State testState = testContext.getState(State::Energy | State::Forces);
        const vector<Vec3>& refForces = refState.getForces();
        const vector<Vec3>& testForces = testState.getForces();

        ASSERT_EQUAL_TOL(refState.getPotentialEnergy(), testState.getPotentialEnergy(), TOL);
        for (int i = 0; i < refForces.size(); i++) {
            ASSERT_EQUAL_VEC(refForces[i], testForces[i], TOL);
        }
    }
}

void testCoulombGaussian() {
    // Ensures that Gaussian charge interactions are computed correctly by
    // comparing a ConstantPotentialForce-containing system with Gaussian
    // charges to a NonbondedForce-containing system with point charges.

    System testSystem;
    ConstantPotentialForce* testForce = new ConstantPotentialForce();
    testSystem.setDefaultPeriodicBoxVectors(Vec3(10, 0, 0), Vec3(0, 10, 0), Vec3(0, 0, 10));
    testSystem.addParticle(0);
    testSystem.addParticle(0);
    testSystem.addParticle(0);
    testForce->addParticle(2);
    testForce->addParticle(0);
    testForce->addParticle(0);
    set<int> electrode1{1};
    testForce->addElectrode(electrode1, 0, 0.2, 0);
    set<int> electrode2{2};
    testForce->addElectrode(electrode2, 0, 0.3, 0);
    testSystem.addForce(testForce);

    vector<Vec3> positions{Vec3(0, 0, 0), Vec3(0.9, 0, 0), Vec3(0.6, 0.4, 0)};

    VerletIntegrator testIntegrator(0.001);
    Context testContext(testSystem, testIntegrator, platform);
    testContext.setPositions(positions);

    vector<double> charges;
    testForce->getCharges(testContext, charges);

    System refSystem;
    NonbondedForce* refForce = new NonbondedForce();
    refSystem.setDefaultPeriodicBoxVectors(Vec3(10, 0, 0), Vec3(0, 10, 0), Vec3(0, 0, 10));
    refSystem.addParticle(0);
    refSystem.addParticle(0);
    refSystem.addParticle(0);
    refForce->addParticle(charges[0], 1, 0);
    refForce->addParticle(charges[1], 1, 0);
    refForce->addParticle(charges[2], 1, 0);
    refForce->setNonbondedMethod(NonbondedForce::PME);
    refSystem.addForce(refForce);

    VerletIntegrator refIntegrator(0.001);
    Context refContext(refSystem, refIntegrator, platform);
    refContext.setPositions(positions);

    State testState = testContext.getState(State::Energy);
    State refState = refContext.getState(State::Energy);

    Vec3 x01 = positions[1] - positions[0];
    Vec3 x02 = positions[2] - positions[0];
    Vec3 x12 = positions[2] - positions[1];
    double r01 = sqrt(x01.dot(x01));
    double r02 = sqrt(x02.dot(x02));
    double r12 = sqrt(x12.dot(x12));

    double du_pair = -ONE_4PI_EPS0 * (
        charges[0] * charges[1] * erfc(r01 / 0.2) / r01 +
        charges[0] * charges[2] * erfc(r02 / 0.3) / r02 +
        charges[1] * charges[2] * erfc(r12 / sqrt(0.13)) / r12
    );
    double du_self = ONE_4PI_EPS0 * (charges[1] * charges[1] / 0.2 + charges[2] * charges[2] / 0.3) / sqrt(2 * M_PI);

    ASSERT_EQUAL_TOL(refState.getPotentialEnergy() + du_pair + du_self, testState.getPotentialEnergy(), TOL);
}

void testFiniteFieldNonPeriodic() {
    // Ensures that the electric field energy is based on an absolute offset
    // from the origin and not affected by periodic wrapping.

    System testSystem;
    Vec3 a, b, c;
    a = Vec3(10, 0, 0);
    b = Vec3(1, 9, 0);
    c = Vec3(2, 3, 8);
    testSystem.setDefaultPeriodicBoxVectors(a, b, c);
    testSystem.addParticle(0);

    Vec3 field(123, -456, 789);
    double charge = -2;

    ConstantPotentialForce* testForce = new ConstantPotentialForce();
    testForce->setEwaldErrorTolerance(1e-4);
    testForce->addParticle(charge);
    testSystem.addForce(testForce);

    VerletIntegrator baseIntegrator(0.001);
    Context baseContext(testSystem, baseIntegrator, platform);
    baseContext.setPositions(vector<Vec3>{Vec3(0, 0, 0)});
    double pmeBaseEnergy = baseContext.getState(State::Energy).getPotentialEnergy();

    testForce->setExternalField(field);

    Vec3 position(30, 20, -10);
    double referenceEnergy = pmeBaseEnergy - charge * field.dot(position);
    Vec3 referenceForce = charge * field;

    VerletIntegrator testIntegrator(0.001);
    Context testContext(testSystem, testIntegrator, platform);
    testContext.setPositions(vector<Vec3>{position});
    testIntegrator.step(1);

    State testState = testContext.getState(State::Energy | State::Forces);
    ASSERT_EQUAL_TOL(referenceEnergy, testState.getPotentialEnergy(), TOL);
    ASSERT_EQUAL_VEC(referenceForce, testState.getForces()[0], TOL);
}

void testElectrodesDisjoint() {
    // Ensures that a particle cannot belong to more than one electrode.

    System system;
    ConstantPotentialForce* force = new ConstantPotentialForce();
    for (int i = 0; i < 3; i++) {
        system.addParticle(1);
        force->addParticle(0);
    }
    set<int> electrode1{0, 1};
    set<int> electrode2{1, 2};
    force->addElectrode(electrode1, 0, 0.2, 0);
    force->addElectrode(electrode2, 0, 0.2, 0);
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

void testNoElectrodeExceptions() {
    // Ensures that a particle cannot belong to both an electrode and an
    // exception.

    System system;
    ConstantPotentialForce* force = new ConstantPotentialForce();
    for (int i = 0; i < 3; i++) {
        system.addParticle(1);
        force->addParticle(0);
    }
    set<int> electrode{0, 1};
    force->addElectrode(electrode, 0, 0.2, 0);
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
    set<int> electrode1{0, 1};
    force1->addElectrode(electrode1, 0, 0.2, 0);
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
    set<int> electrode2{0, 1};
    force2->addElectrode(electrode2, 0, 0.2, 0);
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

void testSmallSystems(ConstantPotentialForce::ConstantPotentialMethod method, bool usePreconditioner) {
    // Ensures that solver convergence is achieved on various small systems with
    // zero or one degrees of freedom.

    // Using one electrode particle with a total charge constraint, the solution
    // is uniquely determined.

    System system1;
    ConstantPotentialForce* force1 = new ConstantPotentialForce();
    system1.addParticle(0);
    system1.addParticle(1);
    force1->addParticle(0);
    force1->addParticle(2);
    set<int> electrode1{0};
    force1->addElectrode(electrode1, 0, 0.2, 0);
    force1->setConstantPotentialMethod(method);
    force1->setUsePreconditioner(usePreconditioner);
    force1->setUseChargeConstraint(true);
    force1->setChargeConstraintTarget(5);
    system1.addForce(force1);
    VerletIntegrator integrator1(0.001);
    Context context1(system1, integrator1, platform);
    vector<Vec3> positions1{Vec3(0, 0, 0), Vec3(0.5, 0, 0)};
    context1.setPositions(positions1);
    vector<double> charges1;
    force1->getCharges(context1, charges1);
    ASSERT_EQUAL(charges1.size(), 2);
    ASSERT_EQUAL_TOL(3, charges1[0], TOL);
    ASSERT_EQUAL(2, charges1[1]);

    // Using two electrode particles with a total charge constraint, there is
    // one degree of freedom.  Make a symmetric system so that we know the
    // solution.

    System system2;
    ConstantPotentialForce* force2 = new ConstantPotentialForce();
    system2.addParticle(0);
    system2.addParticle(0);
    system2.addParticle(1);
    force2->addParticle(0);
    force2->addParticle(0);
    force2->addParticle(2);
    set<int> electrode2{0, 1};
    force2->addElectrode(electrode2, 0, 0.2, 0);
    force2->setConstantPotentialMethod(method);
    force2->setUsePreconditioner(usePreconditioner);
    force2->setUseChargeConstraint(true);
    force2->setChargeConstraintTarget(5);
    system2.addForce(force2);
    VerletIntegrator integrator2(0.001);
    Context context2(system2, integrator2, platform);
    vector<Vec3> positions2{Vec3(-0.25, 0, 0), Vec3(0.25, 0, 0), Vec3(0, 0, 0)};
    context2.setPositions(positions2);
    vector<double> charges2;
    force2->getCharges(context2, charges2);
    ASSERT_EQUAL(charges2.size(), 3);
    ASSERT_EQUAL_TOL(1.5, charges2[0], TOL);
    ASSERT_EQUAL_TOL(1.5, charges2[1], TOL);
    ASSERT_EQUAL(2, charges2[2]);

    // Using no electrode particles without a total charge constraint, we should
    // simply recover the fixed charges.

    System system3;
    ConstantPotentialForce* force3 = new ConstantPotentialForce();
    system3.addParticle(1);
    system3.addParticle(1);
    force3->addParticle(2);
    force3->addParticle(3);
    force3->setConstantPotentialMethod(method);
    force3->setUsePreconditioner(usePreconditioner);
    force3->setUseChargeConstraint(false);
    system3.addForce(force3);
    VerletIntegrator integrator3(0.001);
    Context context3(system3, integrator3, platform);
    vector<Vec3> positions3{Vec3(0, 0, 0), Vec3(0.5, 0, 0)};
    context3.setPositions(positions3);
    vector<double> charges3;
    force3->getCharges(context3, charges3);
    ASSERT_EQUAL(charges3.size(), 2);
    ASSERT_EQUAL(2, charges3[0]);
    ASSERT_EQUAL(3, charges3[1]);

    // Using one electrode particle without a total charge constraint, there is
    // one degree of freedom.  Make a symmetric system so that we know the
    // solution.

    System system4;
    ConstantPotentialForce* force4 = new ConstantPotentialForce();
    system4.addParticle(0);
    system4.addParticle(1);
    system4.addParticle(1);
    force4->addParticle(0);
    force4->addParticle(1);
    force4->addParticle(-1);
    set<int> electrode4{0};
    force4->addElectrode(electrode4, 0, 0.2, 0);
    force4->setConstantPotentialMethod(method);
    force4->setUsePreconditioner(usePreconditioner);
    force4->setUseChargeConstraint(false);
    system4.addForce(force4);
    VerletIntegrator integrator4(0.001);
    Context context4(system4, integrator4, platform);
    vector<Vec3> positions4{Vec3(0, 0, 0), Vec3(-0.75, 0, 0), Vec3(0.75, 0, 0)};
    context4.setPositions(positions4);
    vector<double> charges4;
    force4->getCharges(context4, charges4);
    ASSERT_EQUAL(charges4.size(), 3);
    ASSERT_EQUAL_TOL(0, charges4[0], TOL);
    ASSERT_EQUAL(1, charges4[1]);
    ASSERT_EQUAL(-1, charges4[2]);
}

void testNoConstraintWithoutElectrode(ConstantPotentialForce::ConstantPotentialMethod method, bool usePreconditioner) {
    // If there are no electrodes (no charge degrees of freedom in the system),
    // the user should not be allowed to specify a total charge constraint.

    System system;
    ConstantPotentialForce* force = new ConstantPotentialForce();
    system.addParticle(1);
    system.addParticle(1);
    force->addParticle(1);
    force->addParticle(1);
    force->setConstantPotentialMethod(method);
    force->setUsePreconditioner(usePreconditioner);
    force->setUseChargeConstraint(true);
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

void testConstrainCharge(ConstantPotentialForce::ConstantPotentialMethod method, bool usePreconditioner) {
    // Ensures that the total charge constraint works correctly, including when
    // the system charge target is set to a non-neutral value and the
    // non-electrode particles carry a net charge.

    System system;
    ConstantPotentialForce* force = new ConstantPotentialForce();

    system.setDefaultPeriodicBoxVectors(Vec3(10, 0, 0), Vec3(0, 10, 0), Vec3(0, 0, 10));

    vector<Vec3> positions;
    set<int> electrode;
    double electrolyteCharge = 0;
    for (int i = 0; i < 10; i++) {
        system.addParticle(0);
        force->addParticle(0);
        positions.push_back(Vec3(i, 0, 0));
        electrode.insert(i);
    }
    for (int i = 0; i < 10; i++) {
        system.addParticle(1);
        double charge = (i % 2 ? -1 : 1) + 0.01 * i * i;
        force->addParticle(charge);
        positions.push_back(Vec3(i, 5, 5));
        electrolyteCharge += charge;
    }
    double chargeTarget = -2;

    force->addElectrode(electrode, 0, 0.2, 0);
    force->setConstantPotentialMethod(method);
    force->setUsePreconditioner(usePreconditioner);
    force->setUseChargeConstraint(true);
    force->setChargeConstraintTarget(chargeTarget);

    system.addForce(force);

    VerletIntegrator integrator(0.001);
    Context context(system, integrator, platform);
    context.setPositions(positions);
    vector<double> charges;
    force->getCharges(context, charges);

    ASSERT_EQUAL(20, charges.size());
    double electrodeCharge = 0;
    for (int i = 0; i < 10; i++) {
        electrodeCharge += charges[i];
    }
    ASSERT_EQUAL_TOL(chargeTarget - electrolyteCharge, electrodeCharge, TOL);
}

void makeTestUpdateSystem(ConstantPotentialForce::ConstantPotentialMethod method, bool usePreconditioner, System& system, ConstantPotentialForce*& force, vector<Vec3>& positions) {
    // Makes a reference system for tests below checking that the properties of
    // a system can be changed after context creation.

    system = System();
    system.setDefaultPeriodicBoxVectors(Vec3(3, 0, 0), Vec3(0, 3, 0), Vec3(0, 0, 3));
    force = new ConstantPotentialForce();
    system.addParticle(1);
    system.addParticle(1);
    system.addParticle(1);
    system.addParticle(0);
    system.addParticle(0);
    force->addParticle(1);
    force->addParticle(2);
    force->addParticle(-3);
    force->addParticle(0);
    force->addParticle(0);
    force->addException(0, 1, 1.5);
    force->addException(1, 2, 0);
    force->addElectrode({3}, 1, 0.02, 0.3);
    force->addElectrode({4}, 2, 0.04, 0.5);
    force->setConstantPotentialMethod(method);
    force->setUsePreconditioner(usePreconditioner);
    force->setUseChargeConstraint(true);
    force->setChargeConstraintTarget(5);
    force->setExternalField(Vec3(1, 2, 3));
    force->setCutoffDistance(1.1);
    system.addForce(force);
    
    positions.clear();
    positions.push_back(Vec3(0, 0, 0));
    positions.push_back(Vec3(0, 0, 1));
    positions.push_back(Vec3(0, 1, 1));
    positions.push_back(Vec3(1, 1, 1));
    positions.push_back(Vec3(2, 2, 2));
}

void testUpdate(ConstantPotentialForce::ConstantPotentialMethod method, bool usePreconditioner) {
    // Ensures that parameters can be updated and that the results match those
    // computed by a completely different system.

    System system1, system2;
    ConstantPotentialForce* force1, * force2;
    vector<Vec3> positions1, positions2;
    VerletIntegrator integrator1(0.001), integrator2(0.001);
    vector<double> charges1, charges2;
    makeTestUpdateSystem(method, usePreconditioner, system1, force1, positions1);
    makeTestUpdateSystem(method, usePreconditioner, system2, force2, positions2);

    // Make sure to get charges before updating so that the constant potential
    // method initializes: we are testing that it reinitializes with the update.
    Context context1(system1, integrator1, platform);
    context1.setPositions(positions1);
    force1->getCharges(context1, charges1);
    force1->setParticleParameters(1, 3);
    force1->setParticleParameters(2, -5);
    force1->setExceptionParameters(0, 0, 1, 1.75);
    force1->setElectrodeParameters(0, {3}, 3, 0.06, 0.7);
    force1->setChargeConstraintTarget(-7);
    force1->setExternalField(Vec3(-6, -4, -2));
    force1->updateParametersInContext(context1);
    force1->getCharges(context1, charges1);

    // Reinitialize the reference system completely to force an update.
    force2->setParticleParameters(1, 3);
    force2->setParticleParameters(2, -5);
    force2->setExceptionParameters(0, 0, 1, 1.75);
    force2->setElectrodeParameters(0, {3}, 3, 0.06, 0.7);
    force2->setChargeConstraintTarget(-7);
    force2->setExternalField(Vec3(-6, -4, -2));
    Context context2(system2, integrator2, platform);
    context2.setPositions(positions2);
    force2->getCharges(context2, charges2);

    {
        // Test charges.
        for (int i = 0; i < charges2.size(); i++) {
            ASSERT_EQUAL_TOL(charges2[i], charges1[i], TOL);
        }

        // Get states and test energies and forces.
        State state1 = context1.getState(State::Energy | State::Forces);
        State state2 = context2.getState(State::Energy | State::Forces);
        const vector<Vec3>& forces1 = state1.getForces();
        const vector<Vec3>& forces2 = state2.getForces();
        ASSERT_EQUAL_TOL(state2.getPotentialEnergy(), state1.getPotentialEnergy(), TOL);
        for (int i = 0; i < forces2.size(); i++) {
            ASSERT_EQUAL_VEC(forces2[i], forces1[i], TOL);
        }
    }

    // Now test updating the box.  Get PME parameters to apply to the second
    // context so that energies and forces can be compared.
    double alpha;
    int nx, ny, nz;
    force1->getPMEParametersInContext(context1, alpha, nx, ny, nz);

    context1.setPeriodicBoxVectors(Vec3(3.25, 0, 0), Vec3(0, 3.5, 0), Vec3(0, 0, 3.75));
    force1->getCharges(context1, charges1);

    force2->setPMEParameters(alpha, nx, ny, nz);
    system2.setDefaultPeriodicBoxVectors(Vec3(3.25, 0, 0), Vec3(0, 3.5, 0), Vec3(0, 0, 3.75));
    context2.reinitialize(false);
    context2.setPositions(positions2);
    force2->getCharges(context2, charges2);

    {
        // Test charges.
        for (int i = 0; i < charges2.size(); i++) {
            ASSERT_EQUAL_TOL(charges2[i], charges1[i], TOL);
        }

        // Get states and test energies and forces.
        State state1 = context1.getState(State::Energy | State::Forces);
        State state2 = context2.getState(State::Energy | State::Forces);
        const vector<Vec3>& forces1 = state1.getForces();
        const vector<Vec3>& forces2 = state2.getForces();
        ASSERT_EQUAL_TOL(state2.getPotentialEnergy(), state1.getPotentialEnergy(), TOL);
        for (int i = 0; i < forces2.size(); i++) {
            ASSERT_EQUAL_VEC(forces2[i], forces1[i], TOL);
        }
    }

    // Now test updating positions.  This checks that, e.g., a precomputed
    // matrix is invalidated and recomputed when electrode particle positions
    // are adjusted.
    positions1[0] += Vec3(0.1, 0.2, 0.3);
    positions1[4] += Vec3(0.4, 0.5, 0.6);
    context1.setPositions(positions1);
    force1->getCharges(context1, charges1);

    context2.reinitialize(false);
    context2.setPositions(positions1);
    force2->getCharges(context2, charges2);

    {
        // Test charges.
        for (int i = 0; i < charges2.size(); i++) {
            ASSERT_EQUAL_TOL(charges2[i], charges1[i], TOL);
        }

        // Get states and test energies and forces.
        State state1 = context1.getState(State::Energy | State::Forces);
        State state2 = context2.getState(State::Energy | State::Forces);
        const vector<Vec3>& forces1 = state1.getForces();
        const vector<Vec3>& forces2 = state2.getForces();
        ASSERT_EQUAL_TOL(state2.getPotentialEnergy(), state1.getPotentialEnergy(), TOL);
        for (int i = 0; i < forces2.size(); i++) {
            ASSERT_EQUAL_VEC(forces2[i], forces1[i], TOL);
        }
    }
}

void testParallelPlateCapacitorDoubleCell() {
    // Uses the zero field double-cell geometry to test an ideal parallel plate
    // capacitor with vacuum between the plates.

    double capacitance = EPSILON0 * 2.0 * 3.0 / 7.0;

    double potential = 100.0;

    System system;
    ConstantPotentialForce* force = new ConstantPotentialForce();

    system.setDefaultPeriodicBoxVectors(Vec3(2, 0, 0), Vec3(0, 3, 0), Vec3(0, 0, 16));

    vector<Vec3> positions;
    set<int> electrode1, electrode2;
    for (int iz = 0; iz < 4; iz++) {
        for (int iy = 0; iy < 12; iy++) {
            for (int ix = 0; ix < 8; ix++) {
                electrode1.insert(force->getNumParticles());
                system.addParticle(0);
                force->addParticle(0);
                positions.push_back(Vec3(0.25 * ix, 0.25 * iy, 0.25 * iz));

                electrode2.insert(force->getNumParticles());
                system.addParticle(0);
                force->addParticle(0);
                positions.push_back(Vec3(0.25 * ix, 0.25 * iy, 0.25 * iz + 8.0));
            }
        }
    }
    force->addElectrode(electrode1, 0, 0.2, 0);
    force->addElectrode(electrode2, potential, 0.2, 0);
    force->setUseChargeConstraint(true);

    system.addForce(force);

    VerletIntegrator integrator(0.001);
    Context context(system, integrator, platform);
    context.setPositions(positions);
    vector<double> charges;
    double q1, q2;

    force->getCharges(context, charges);

    q1 = 0;
    for (int ii : electrode1) {
        q1 += charges[ii];
    }

    q2 = 0;
    for (int ii : electrode2) {
        q2 += charges[ii];
    }

    // Charge on each electrode is doubled since we are using a double-cell.
    // Use a looser tolerance since the analytical formula is for uniform
    // plates with defined edges, not cubic arrays of Gaussian charges.
    ASSERT_EQUAL_TOL(-2.0 * potential * capacitance, q1, 1e-3);
    ASSERT_EQUAL_TOL(2.0 * potential * capacitance, q2, 1e-3);

    // Update potential to a new value and recompute.
    potential = -200.0;
    force->setElectrodeParameters(1, electrode2, potential, 0.2, 0);
    force->updateParametersInContext(context);
    force->getCharges(context, charges);

    q1 = 0;
    for (int ii : electrode1) {
        q1 += charges[ii];
    }

    q2 = 0;
    for (int ii : electrode2) {
        q2 += charges[ii];
    }

    ASSERT_EQUAL_TOL(-2.0 * potential * capacitance, q1, 1e-3);
    ASSERT_EQUAL_TOL(2.0 * potential * capacitance, q2, 1e-3);
}

void testParallelPlateCapacitorFiniteField() {
    // Uses the finite field geometry (still applying a potential to one
    // electrode) to test an ideal parallel plate capacitor.

    double capacitance = EPSILON0 * 2.0 * 3.0 / 11.0;

    double potential = 100.0;

    System system;
    ConstantPotentialForce* force = new ConstantPotentialForce();

    system.setDefaultPeriodicBoxVectors(Vec3(2, 0, 0), Vec3(0, 3, 0), Vec3(0, 0, 16));

    vector<Vec3> positions;
    set<int> electrode1, electrode2;
    for (int iz = 0; iz < 4; iz++) {
        for (int iy = 0; iy < 12; iy++) {
            for (int ix = 0; ix < 8; ix++) {
                electrode1.insert(force->getNumParticles());
                system.addParticle(0);
                force->addParticle(0);
                positions.push_back(Vec3(0.25 * ix, 0.25 * iy, 0.25 * iz));

                electrode2.insert(force->getNumParticles());
                system.addParticle(0);
                force->addParticle(0);
                positions.push_back(Vec3(0.25 * ix, 0.25 * iy, 0.25 * iz + 12.0));
            }
        }
    }
    force->addElectrode(electrode1, 0, 0.2, 0);
    force->addElectrode(electrode2, potential, 0.2, 0);
    force->setUseChargeConstraint(true);
    force->setExternalField(Vec3(0, 0, -potential / 16));

    system.addForce(force);

    VerletIntegrator integrator(0.001);
    Context context(system, integrator, platform);
    context.setPositions(positions);
    vector<double> charges;
    force->getCharges(context, charges);

    double q1 = 0;
    for (int ii : electrode1) {
        q1 += charges[ii];
    }

    double q2 = 0;
    for (int ii : electrode2) {
        q2 += charges[ii];
    }

    // Use a looser tolerance since the analytical formula is for uniform
    // plates with defined edges, not cubic arrays of Gaussian charges.
    ASSERT_EQUAL_TOL(-potential * capacitance, q1, 1e-3);
    ASSERT_EQUAL_TOL(potential * capacitance, q2, 1e-3);

    // Update potential to a new value and recompute.
    potential = -200.0;
    force->setElectrodeParameters(1, electrode2, potential, 0.2, 0);
    force->setExternalField(Vec3(0, 0, -potential / 16));
    force->updateParametersInContext(context);
    force->getCharges(context, charges);

    q1 = 0;
    for (int ii : electrode1) {
        q1 += charges[ii];
    }

    q2 = 0;
    for (int ii : electrode2) {
        q2 += charges[ii];
    }

    ASSERT_EQUAL_TOL(-potential * capacitance, q1, 1e-3);
    ASSERT_EQUAL_TOL(potential * capacitance, q2, 1e-3);
}

void makeTestReferenceSystem(bool freezeAll, System& system, ConstantPotentialForce*& constantPotentialForce, vector<Vec3>& positions) {
    // Generates a test system with SPC/E water in contact with a model
    // electrode.

    system = System();
    system.setDefaultPeriodicBoxVectors(Vec3(3, 0, 0), Vec3(0, 2, 0), Vec3(0, 0, 2));

    HarmonicBondForce* bondForce = new HarmonicBondForce();
    NonbondedForce* nonbondedForce = new NonbondedForce();
    constantPotentialForce = new ConstantPotentialForce();

    for (int i = 0; i < 300; i++) {
        system.addParticle((i >= 100 && i < 200) || freezeAll ? 0 : 195.08);
        nonbondedForce->addParticle(0, 0.25, 1);
        constantPotentialForce->addParticle(0);
    }

    // static const int bondCount = ...;
    // static const int bondIndices[...][2] = {...};
#include "conp_bond_indices.dat"

    for (int i = 0; i < bondCount; i++) {
        bondForce->addBond(bondIndices[i][0], bondIndices[i][1], 0.28284271247461906, 10000);
    }

    for (int i = 0; i < 222; i++) {
        int i1 = system.addParticle(15.99943);
        int i2 = system.addParticle(1.007947);
        int i3 = system.addParticle(1.007947);
        system.addConstraint(i1, i2, 0.1);
        system.addConstraint(i1, i3, 0.1);
        system.addConstraint(i2, i3, 0.1632980861841278);
        nonbondedForce->addParticle(0, 0.3165719505039882, 0.6497752);
        nonbondedForce->addParticle(0, 1, 0);
        nonbondedForce->addParticle(0, 1, 0);
        nonbondedForce->addException(i1, i2, 0, 1, 0);
        nonbondedForce->addException(i1, i3, 0, 1, 0);
        nonbondedForce->addException(i2, i3, 0, 1, 0);
        constantPotentialForce->addParticle(-0.8476);
        constantPotentialForce->addParticle(0.4238);
        constantPotentialForce->addParticle(0.4238);
        constantPotentialForce->addException(i1, i2, 0);
        constantPotentialForce->addException(i1, i3, 0);
        constantPotentialForce->addException(i2, i3, 0);
    }

    bondForce->setUsesPeriodicBoundaryConditions(true);
    system.addForce(bondForce);

    nonbondedForce->setNonbondedMethod(NonbondedForce::CutoffPeriodic);
    system.addForce(nonbondedForce);

    positions.clear();
    positions.resize(system.getNumParticles());

    // positions[...] = Vec3(..., ..., ...);
#include "conp_positions.dat"

    set<int> electrodeLeft;
    set<int> electrodeRight;
    for (int i = 0; i < 300; i++) {
        if (positions[i][0] < 1.5) {
            electrodeLeft.insert(i);
        } else {
            electrodeRight.insert(i);
        }
    }

    // Initialize the electrodes (parameters will be set later).
    constantPotentialForce->addElectrode(electrodeLeft, 0, 0, 0);
    constantPotentialForce->addElectrode(electrodeRight, 0, 0, 0);
    system.addForce(constantPotentialForce);
}

void testReferenceCharges(bool testThomasFermi, ConstantPotentialForce::ConstantPotentialMethod method, bool usePreconditioner) {
    // Make sure that charges solved for match values computed by an external
    // reference implementation.

    System testSystem;
    ConstantPotentialForce* testForce;
    vector<Vec3> positions;

    makeTestReferenceSystem(method == ConstantPotentialForce::Matrix, testSystem, testForce, positions);

    testForce->setEwaldErrorTolerance(5e-5);
    testForce->setConstantPotentialMethod(method);
    testForce->setUsePreconditioner(usePreconditioner);
    testForce->setUseChargeConstraint(true);
    testForce->setChargeConstraintTarget(0);
    testForce->setCutoffDistance(1);

    if (!testThomasFermi) {
        // The reference implementation (MetalWalls) doesn't support
        // simultaneous specification of these options.
        testForce->setExternalField(Vec3(100, 50, 25));
    }

    set<int> electrodeParticles;
    double potential, gaussianWidth, thomasFermiScale;
    testForce->getElectrodeParameters(0, electrodeParticles, potential, gaussianWidth, thomasFermiScale);
    testForce->setElectrodeParameters(0, electrodeParticles, 100, 0.05, testThomasFermi ? 0.50625 : 0.0);
    testForce->getElectrodeParameters(1, electrodeParticles, potential, gaussianWidth, thomasFermiScale);
    testForce->setElectrodeParameters(1, electrodeParticles, 200, 0.1, testThomasFermi ? 0.225 : 0.0);

    VerletIntegrator integrator(0.001);
    Context testContext(testSystem, integrator, platform);

    testContext.setPositions(positions);

    // Charges of electrode atoms should match reference values.
    vector<double> refCharges(300);
    if (testThomasFermi) {
        // refCharges[...] = ...;
#include "conp_ref_charges_tf.dat"
    } else {
        // refCharges[...] = ...;
#include "conp_ref_charges.dat"
    }
    vector<double> testCharges;
    testForce->getCharges(testContext, testCharges);
    for (int i = 0; i < 300; i++) {
        ASSERT_EQUAL_TOL(refCharges[i], testCharges[i], TOL);
    }
}

void testChargeUpdate(ConstantPotentialForce::ConstantPotentialMethod method, bool usePreconditioner) {
    // Make sure that charges get updated correctly before and after dynamics.

    System testSystem;
    ConstantPotentialForce* testForce;
    vector<Vec3> positions;

    makeTestReferenceSystem(method == ConstantPotentialForce::Matrix, testSystem, testForce, positions);

    testForce->setConstantPotentialMethod(method);
    testForce->setUsePreconditioner(usePreconditioner);
    testForce->setUseChargeConstraint(true);
    testForce->setChargeConstraintTarget(1);
    testForce->setExternalField(Vec3(10, 0, 0));

    set<int> electrodeParticles;
    double potential, gaussianWidth, thomasFermiScale;
    testForce->getElectrodeParameters(0, electrodeParticles, potential, gaussianWidth, thomasFermiScale);
    testForce->setElectrodeParameters(0, electrodeParticles, 1, 0.05, 1);
    testForce->getElectrodeParameters(1, electrodeParticles, potential, gaussianWidth, thomasFermiScale);
    testForce->setElectrodeParameters(1, electrodeParticles, 2, 0.06, 1.2);

    VerletIntegrator integrator(0.001);
    Context testContext(testSystem, integrator, platform);

    testContext.setPositions(positions);
    testContext.setVelocitiesToTemperature(300);

    // Charges of non-electrode atoms should exactly match charges set.
    vector<double> charges0;
    testForce->getCharges(testContext, charges0);
    for (int i = 300; i < charges0.size(); i++) {
        double refCharge;
        testForce->getParticleParameters(i, refCharge);
        ASSERT_EQUAL(refCharge, charges0[i]);
    }

    integrator.step(1);

    // Charges of electrode atoms should be changed, and charges of
    // non-electrode atoms should exactly match charges set.
    vector<double> charges1;
    testForce->getCharges(testContext, charges1);
    for (int i = 0; i < 300; i++) {
        ASSERT(charges1[i] != charges0[i]);
    }
    for (int i = 300; i < charges1.size(); i++) {
        double refCharge;
        testForce->getParticleParameters(i, refCharge);
        ASSERT_EQUAL(refCharge, charges1[i]);
    }
}

void testEnergyConservation(ConstantPotentialForce::ConstantPotentialMethod method, bool usePreconditioner, int stepCount) {
    // Do a short dynamics run and ensure that energy is conserved.

    System system;
    system.setDefaultPeriodicBoxVectors(Vec3(2, 0, 0), Vec3(0, 2, 0), Vec3(0, 0, 3));

    // WCA potential to keep particles from getting too close to each other.
    CustomNonbondedForce* nonbondedForce = new CustomNonbondedForce("1/r^12-2/r^6+1");
    nonbondedForce->setNonbondedMethod(CustomNonbondedForce::CutoffPeriodic);
    nonbondedForce->setCutoffDistance(1);

    ConstantPotentialForce* constantPotentialForce = new ConstantPotentialForce();
    constantPotentialForce->setCutoffDistance(1);
    constantPotentialForce->setUseChargeConstraint(true);
    constantPotentialForce->setChargeConstraintTarget(1);
    constantPotentialForce->setExternalField(Vec3(1, 2, 4));
    constantPotentialForce->setConstantPotentialMethod(method);
    constantPotentialForce->setUsePreconditioner(usePreconditioner);
    constantPotentialForce->setEwaldErrorTolerance(1e-4);

    for (int i = 0; i < 8; i++) {
        system.addParticle(method == ConstantPotentialForce::Matrix ? 0 : 1);
        constantPotentialForce->addParticle(0);
    }
    for (int i = 0; i < 4; i++) {
        system.addParticle(1);
    }
    for (int i = 0; i < 12; i++) {
        nonbondedForce->addParticle();
    }
    constantPotentialForce->addParticle(1);
    constantPotentialForce->addParticle(-1);
    constantPotentialForce->addParticle(-1);
    constantPotentialForce->addParticle(1);
    constantPotentialForce->addElectrode({0, 1, 2, 3}, -8, 0.05, 1);
    constantPotentialForce->addElectrode({4, 5, 6, 7}, 8, 0.1, 1.5);

    system.addForce(nonbondedForce);
    system.addForce(constantPotentialForce);

    vector<Vec3> positions{
        Vec3(0.0, 0.0, 0.0),
        Vec3(0.0, 1.0, 0.0),
        Vec3(1.0, 0.0, 0.0),
        Vec3(1.0, 1.0, 0.0),
        Vec3(0.0, 0.0, 2.0),
        Vec3(0.0, 1.0, 2.0),
        Vec3(1.0, 0.0, 2.0),
        Vec3(1.0, 1.0, 2.0),
        Vec3(0.5, 0.5, 1.0),
        Vec3(0.5, 1.5, 1.0),
        Vec3(1.5, 0.5, 1.0),
        Vec3(1.5, 1.5, 1.0)
    };

    VerletIntegrator integrator(0.002);
    Context context(system, integrator, platform);
    context.setPositions(positions);
    context.setVelocitiesToTemperature(1 / BOLTZ);

    State state1 = context.getState(State::Energy);
    double energy1 = state1.getPotentialEnergy() + state1.getKineticEnergy();

    integrator.step(stepCount);

    State state2 = context.getState(State::Energy);
    double energy2 = state2.getPotentialEnergy() + state2.getKineticEnergy();

    ASSERT_USUALLY_EQUAL_TOL(energy1, energy2, 5e-4);
}

void compareToReferencePlatform(System& system, ConstantPotentialForce* force, const vector<Vec3>& positions) {
    // Compares results from the current platform to results computed by the
    // reference platform.

    VerletIntegrator refIntegrator(0.001);
    VerletIntegrator testIntegrator(0.001);
    ReferencePlatform refPlatform;
    Context refContext(system, refIntegrator, refPlatform);
    Context testContext(system, testIntegrator, platform);
    refContext.setPositions(positions);
    testContext.setPositions(positions);

    State refState = refContext.getState(State::Energy | State::Forces);
    State testState = testContext.getState(State::Energy | State::Forces);
    const vector<Vec3>& refForces = refState.getForces();
    const vector<Vec3>& testForces = testState.getForces();
    vector<double> refCharges, testCharges;
    force->getCharges(refContext, refCharges);
    force->getCharges(testContext, testCharges);

    ASSERT_EQUAL_TOL(refState.getPotentialEnergy(), testState.getPotentialEnergy(), TOL);
    for (int i = 0; i < system.getNumParticles(); i++) {
        ASSERT_EQUAL_VEC(refForces[i], testForces[i], 3e-3);
    }
    for (int i = 0; i < system.getNumParticles(); i++) {
        ASSERT_EQUAL_TOL(refCharges[i], testCharges[i], TOL);
    }
}

void testCompareToReferencePlatform(ConstantPotentialForce::ConstantPotentialMethod method, bool usePreconditioner) {
    // Compares results between the current and reference platforms for a few
    // test systems.  This is only called in the runPlatformTests() for
    // platforms other than the reference platform.

    System system;
    ConstantPotentialForce* force;
    vector<Vec3> positions;

    makeTestUpdateSystem(method, usePreconditioner, system, force, positions);
    force->setEwaldErrorTolerance(5e-5);
    compareToReferencePlatform(system, force, positions);

    makeTestReferenceSystem(method == ConstantPotentialForce::Matrix, system, force, positions);
    force->setConstantPotentialMethod(method);
    force->setUsePreconditioner(usePreconditioner);
    force->setUseChargeConstraint(true);
    force->setChargeConstraintTarget(1);
    force->setExternalField(Vec3(10, 0, 0));
    force->setEwaldErrorTolerance(5e-5);
    set<int> electrodeParticles;
    double potential, gaussianWidth, thomasFermiScale;
    force->getElectrodeParameters(0, electrodeParticles, potential, gaussianWidth, thomasFermiScale);
    force->setElectrodeParameters(0, electrodeParticles, 1, 0.05, 1);
    force->getElectrodeParameters(1, electrodeParticles, potential, gaussianWidth, thomasFermiScale);
    force->setElectrodeParameters(1, electrodeParticles, 2, 0.06, 1.2);
    compareToReferencePlatform(system, force, positions);
}

void testLargeNeighborList(ConstantPotentialForce::ConstantPotentialMethod method, bool usePreconditioner) {
    // Runs a test where the initial neighbor list should overflow on GPU platforms.

    const int n = 9;
    const double l = 3.0;
    const double scale = l / n;

    System system;
    system.setDefaultPeriodicBoxVectors(Vec3(l, 0, 0), Vec3(0, l, 0), Vec3(0, 0, l));

    ConstantPotentialForce* force = new ConstantPotentialForce();
    force->setConstantPotentialMethod(method);
    force->setUsePreconditioner(usePreconditioner);
    force->setUseChargeConstraint(true);
    force->setCutoffDistance(1.4);
    force->setCGErrorTolerance(5e-4);
    system.addForce(force);

    vector<Vec3> positions;
    set<int> electrodeParticles;
    for (int ix = 0; ix < n; ix++) {
        for (int iy = 0; iy < n; iy++) {
            for (int iz = 0; iz < n; iz++) {
                positions.push_back(scale * Vec3(ix, iy, iz));
                positions.push_back(scale * Vec3(ix + 0.5, iy + 0.5, iz + 0.5));
                electrodeParticles.insert(system.addParticle(0.0));
                system.addParticle(0.0);
                force->addParticle(0.0);
                force->addParticle(1.0);
            }
        }
    }
    force->addElectrode(electrodeParticles, 0.0, 0.01, 0.0);

    VerletIntegrator integrator(0.001);
    Context context(system, integrator, platform);
    context.setPositions(positions);

    // Get charges: if the neighbor list is incomplete, they will not be uniformly equal to -1.
    vector<double> charges;
    force->getCharges(context, charges);
    for (int i : electrodeParticles) {
        ASSERT_EQUAL_TOL(-1.0, charges[i], 2e-3);
    }

    // Run again, this time doing an energy/force calculation before getting charges.
    context.reinitialize();
    context.setPositions(positions);
    context.getState(State::Energy | State::Forces);
    force->getCharges(context, charges);
    for (int i : electrodeParticles) {
        ASSERT_EQUAL_TOL(-1.0, charges[i], 2e-3);
    }
}

void platformInitialize();
void runPlatformTests(ConstantPotentialForce::ConstantPotentialMethod method, bool usePreconditioner);

void runMethodDependentTests(ConstantPotentialForce::ConstantPotentialMethod method, bool usePreconditioner) {
    testSmallSystems(method, usePreconditioner);
    testNoConstraintWithoutElectrode(method, usePreconditioner);
    testConstrainCharge(method, usePreconditioner);
    testUpdate(method, usePreconditioner);
    testReferenceCharges(false, method, usePreconditioner); // External field
    testReferenceCharges(true, method, usePreconditioner);  // Thomas-Fermi
    testChargeUpdate(method, usePreconditioner);
    runPlatformTests(method, usePreconditioner);
}

int main(int argc, char* argv[]) {
    try {
        initializeTests(argc, argv);
        platformInitialize();

        testCoulomb(false); // Non-periodic exceptions
        testCoulomb(true);  // Periodic exceptions
        testCoulombOverlap();
        testCoulombNonNeutral();
        testCoulombGaussian();
        testFiniteFieldNonPeriodic();

        testElectrodesDisjoint();
        testNoElectrodeExceptions();
        testElectrodeMatrixNoMass();

        testParallelPlateCapacitorDoubleCell();
        testParallelPlateCapacitorFiniteField();

        runMethodDependentTests(ConstantPotentialForce::Matrix, false); // Matrix inversion (usePreconditioner ignored)
        runMethodDependentTests(ConstantPotentialForce::CG, false);     // Conjugate gradient (not preconditioned)
        runMethodDependentTests(ConstantPotentialForce::CG, true);      // Conjugate gradient (preconditioned)
    }
    catch(const exception& e) {
        cout << "exception: " << e.what() << endl;
        return 1;
    }
    cout << "Done" << endl;
    return 0;
}
