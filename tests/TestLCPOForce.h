/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit.                   *
 * See https://openmm.org/development.                                        *
 *                                                                            *
 * Portions copyright (c) 2025 Stanford University and the Authors.           *
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
#include "openmm/LCPOForce.h"
#include "openmm/System.h"
#include "openmm/VerletIntegrator.h"
#include "ReferencePlatform.h"
#include "SimTKOpenMMRealType.h"
#include "sfmt/SFMT.h"
#include <iostream>

using namespace OpenMM;
using namespace std;

const double TOL = 5e-5;

void testSingleParticle() {
    System system;
    system.addParticle(1.0);

    LCPOForce* lcpo = new LCPOForce();
    lcpo->setSurfaceTension(5.0);
    lcpo->addParticle(0.2, 10.0, 11.0, 12.0, 13.0);
    system.addForce(lcpo);

    VerletIntegrator integrator(0.001);
    Context context(system, integrator, platform);
    context.setPositions(vector<Vec3>(1));
    State state = context.getState(State::Energy | State::Forces);
    const vector<Vec3>& forces = state.getForces();

    ASSERT_EQUAL_TOL(5.0 * 10.0 * 4.0 * PI_M * 0.2 * 0.2, state.getPotentialEnergy(), TOL);
    ASSERT_EQUAL_VEC(Vec3(0.0, 0.0, 0.0), forces[0], TOL);
}

void testTwoParticles() {
    System system;
    system.addParticle(1.0);
    system.addParticle(1.0);

    LCPOForce* lcpo = new LCPOForce();
    lcpo->setSurfaceTension(5.0);
    lcpo->addParticle(0.1, 10.0, 11.0, 12.0, 13.0);
    lcpo->addParticle(0.2, 14.0, 15.0, 16.0, 17.0);
    system.addForce(lcpo);

    Vec3 offset(0.48, 0.6, 0.64); // unit

    VerletIntegrator integrator(0.001);
    Context context(system, integrator, platform);

    for (int trial = 0; trial < 3; trial++) {
        double r = 0.15 + 0.1 * trial;
        context.setPositions(vector<Vec3>{Vec3(), r * offset});
        State state = context.getState(State::Energy | State::Forces);
        const vector<Vec3>& forces = state.getForces();
        ASSERT_EQUAL_TOL(5.0 * (
            // One-body terms
            10.0 * 4.0 * PI_M * 0.1 * 0.1 +
            14.0 * 4.0 * PI_M * 0.2 * 0.2 +
            // Two-body terms
            (r >= 0.3 ? 0.0 : 11.0 * 2.0 * PI_M * 0.1 * (0.1 - r / 2.0 - (0.1 * 0.1 - 0.2 * 0.2) / (2.0 * r))) +
            (r >= 0.3 ? 0.0 : 15.0 * 2.0 * PI_M * 0.2 * (0.2 - r / 2.0 - (0.2 * 0.2 - 0.1 * 0.1) / (2.0 * r)))),
            state.getPotentialEnergy(), TOL
        );

        double forceExpected = 5.0 * (
            // Two-body terms
            (r >= 0.3 ? 0.0 : 11.0 * 2.0 * PI_M * 0.1 * (-0.5 + (0.1 * 0.1 - 0.2 * 0.2) / (2.0 * r * r))) +
            (r >= 0.3 ? 0.0 : 15.0 * 2.0 * PI_M * 0.2 * (-0.5 + (0.2 * 0.2 - 0.1 * 0.1) / (2.0 * r * r))));
        ASSERT_EQUAL_VEC(forceExpected * offset, forces[0], TOL);
        ASSERT_EQUAL_VEC(forceExpected * -offset, forces[1], TOL);
    }
}

void testThreeParticles() {
    System system;
    system.addParticle(1.0);
    system.addParticle(1.0);
    system.addParticle(1.0);

    LCPOForce* lcpo = new LCPOForce();
    lcpo->setSurfaceTension(5.0);
    lcpo->addParticle(0.22, 10.0, 11.0, 12.0, 13.0);
    lcpo->addParticle(0.24, 14.0, 15.0, 16.0, 17.0);
    lcpo->addParticle(0.26, 18.0, 19.0, 20.0, 21.0);
    system.addForce(lcpo);

    VerletIntegrator integrator(0.001);
    Context context(system, integrator, platform);
    vector<Vec3> positions{
        Vec3(1.1, 1.2, 1.3),
        Vec3(1.4, 1.4, 1.4),
        Vec3(1.3, 1.0, 1.6)
    };
    context.setPositions(positions);
    State state = context.getState(State::Energy | State::Forces);
    const vector<Vec3>& forces = state.getForces();

    Vec3 x12 = positions[1] - positions[0];
    Vec3 x13 = positions[2] - positions[0];
    Vec3 x23 = positions[2] - positions[1];
    double r12 = sqrt(x12.dot(x12));
    double r13 = sqrt(x13.dot(x13));
    double r23 = sqrt(x23.dot(x23));
    double A12 = 2.0 * PI_M * 0.22 * (0.22 - r12 / 2.0 - (0.22 * 0.22 - 0.24 * 0.24) / (2.0 * r12));
    double A13 = 2.0 * PI_M * 0.22 * (0.22 - r13 / 2.0 - (0.22 * 0.22 - 0.26 * 0.26) / (2.0 * r13));
    double A21 = 2.0 * PI_M * 0.24 * (0.24 - r12 / 2.0 - (0.24 * 0.24 - 0.22 * 0.22) / (2.0 * r12));
    double A23 = 2.0 * PI_M * 0.24 * (0.24 - r23 / 2.0 - (0.24 * 0.24 - 0.26 * 0.26) / (2.0 * r23));
    double A31 = 2.0 * PI_M * 0.26 * (0.26 - r13 / 2.0 - (0.26 * 0.26 - 0.22 * 0.22) / (2.0 * r13));
    double A32 = 2.0 * PI_M * 0.26 * (0.26 - r23 / 2.0 - (0.26 * 0.26 - 0.24 * 0.24) / (2.0 * r23));
    double dA12 = 2.0 * PI_M * 0.22 * (-0.5 + (0.22 * 0.22 - 0.24 * 0.24) / (2.0 * r12 * r12)) / r12;
    double dA13 = 2.0 * PI_M * 0.22 * (-0.5 + (0.22 * 0.22 - 0.26 * 0.26) / (2.0 * r13 * r13)) / r13;
    double dA21 = 2.0 * PI_M * 0.24 * (-0.5 + (0.24 * 0.24 - 0.22 * 0.22) / (2.0 * r12 * r12)) / r12;
    double dA23 = 2.0 * PI_M * 0.24 * (-0.5 + (0.24 * 0.24 - 0.26 * 0.26) / (2.0 * r23 * r23)) / r23;
    double dA31 = 2.0 * PI_M * 0.26 * (-0.5 + (0.26 * 0.26 - 0.22 * 0.22) / (2.0 * r13 * r13)) / r13;
    double dA32 = 2.0 * PI_M * 0.26 * (-0.5 + (0.26 * 0.26 - 0.24 * 0.24) / (2.0 * r23 * r23)) / r23;

    ASSERT_EQUAL_TOL(5.0 * (
        // One-body terms
        10.0 * 4.0 * PI_M * 0.22 * 0.22 +
        14.0 * 4.0 * PI_M * 0.24 * 0.24 +
        18.0 * 4.0 * PI_M * 0.26 * 0.26 +
        // Two-body terms
        11.0 * (A12 + A13) + 15.0 * (A21 + A23) + 19.0 * (A31 + A32) +
        // Three-body terms
        12.0 * (A23 + A32) + 16.0 * (A13 + A31) + 20.0 * (A12 + A21) +
        13.0 * (A12 * A23 + A13 * A32) + 17.0 * (A21 * A13 + A23 * A31) + 21.0 * (A31 * A12 + A32 * A21)),
        state.getPotentialEnergy(), TOL
    );
    ASSERT_EQUAL_VEC(5.0 * (
        11.0 * (dA12 * -x12 + dA13 * -x13) + 15.0 * dA21 * -x12 + 19.0 * dA31 * -x13 +
        16.0 * (dA13 * -x13 + dA31 * -x13) + 20.0 * (dA12 * -x12 + dA21 * -x12) +
        13.0 * (dA12 * A23 * -x12 + dA13 * A32 * -x13) + 17.0 * (dA21 * A13 * -x12 + A21 * dA13 * -x13 + A23 * dA31 * -x13) + 21.0 * (dA31 * A12 * -x13 + A31 * dA12 * -x12 + A32 * dA21 * -x12)),
        -forces[0], TOL
    );
    ASSERT_EQUAL_VEC(5.0 * (
        11.0 * dA12 * x12 + 15.0 * (dA21 * x12 + dA23 * -x23) + 19.0 * dA32 * -x23 +
        12.0 * (dA23 * -x23 + dA32 * -x23) + 20.0 * (dA12 * x12 + dA21 * x12) +
        13.0 * (dA12 * A23 * x12 + A12 * dA23 * -x23 + A13 * dA32 * -x23) + 17.0 * (dA21 * A13 * x12 + dA23 * A31 * -x23) + 21.0 * (A31 * dA12 * x12 + dA32 * A21 * -x23 + A32 * dA21 * x12)),
        -forces[1], TOL
    );
    ASSERT_EQUAL_VEC(5.0 * (
        11.0 * dA13 * x13 + 15.0 * dA23 * x23 + 19.0 * (dA31 * x13 + dA32 * x23) +
        12.0 * (dA23 * x23 + dA32 * x23) + 16.0 * (dA13 * x13 + dA31 * x13) +
        13.0 * (A12 * dA23 * x23 + dA13 * A32 * x13 + A13 * dA32 * x23) + 17.0 * (A21 * dA13 * x13 + dA23 * A31 * x23 + A23 * dA31 * x13) + 21.0 * (dA31 * A12 * x13 + dA32 * A21 * x23)),
        -forces[2], TOL
    );
}

void testUsePeriodic(bool usePeriodic) {
    System system;
    system.addParticle(1.0);
    system.addParticle(1.0);

    LCPOForce* lcpo = new LCPOForce();
    lcpo->setSurfaceTension(5.0);
    lcpo->addParticle(0.1, 10.0, 11.0, 12.0, 13.0);
    lcpo->addParticle(0.2, 14.0, 15.0, 16.0, 17.0);
    lcpo->setUsesPeriodicBoundaryConditions(usePeriodic);
    system.addForce(lcpo);

    system.setDefaultPeriodicBoxVectors(Vec3(4.0, 0.0, 0.0), Vec3(0.1, 3.0, 0.0), Vec3(0.2, 0.3, 2.0));

    VerletIntegrator integrator(0.001);
    Context context(system, integrator, platform);

    double energyOneBodyExpected = 5.0 * (10.0 * 4.0 * PI_M * 0.1 * 0.1 + 14.0 * 4.0 * PI_M * 0.2 * 0.2);
    double energyTwoBodyExpected = 5.0 * (
        11.0 * 2.0 * PI_M * 0.1 * (0.1 - 0.25 / 2.0 - (0.1 * 0.1 - 0.2 * 0.2) / (2.0 * 0.25)) +
        15.0 * 2.0 * PI_M * 0.2 * (0.2 - 0.25 / 2.0 - (0.2 * 0.2 - 0.1 * 0.1) / (2.0 * 0.25)));
    double forceExpected = 5.0 * (
        11.0 * 2.0 * PI_M * 0.1 * (-0.5 + (0.1 * 0.1 - 0.2 * 0.2) / (2.0 * 0.25 * 0.25)) +
        15.0 * 2.0 * PI_M * 0.2 * (-0.5 + (0.2 * 0.2 - 0.1 * 0.1) / (2.0 * 0.25 * 0.25)));

    for (int image = -2; image <= 2; image++) {
        context.setPositions(vector<Vec3>{Vec3(0.0, 0.0, 0.0), Vec3(0.25, 0.0, 0.0) + image * Vec3(0.2, 0.3, 2.0)});
        State state = context.getState(State::Energy | State::Forces);
        const vector<Vec3>& forces = state.getForces();
        if (!image || usePeriodic) {
            ASSERT_EQUAL_TOL(energyOneBodyExpected + energyTwoBodyExpected, state.getPotentialEnergy(), TOL);
            ASSERT_EQUAL_VEC(forceExpected * Vec3(1.0, 0.0, 0.0), forces[0], TOL);
            ASSERT_EQUAL_VEC(forceExpected * Vec3(-1.0, 0.0, 0.0), forces[1], TOL);
        }
        else {
            ASSERT_EQUAL_TOL(energyOneBodyExpected, state.getPotentialEnergy(), TOL);
            ASSERT_EQUAL_VEC(Vec3(0.0, 0.0, 0.0), forces[0], TOL);
            ASSERT_EQUAL_VEC(Vec3(0.0, 0.0, 0.0), forces[1], TOL);
        }
    }
}

void testIncludeZeroScale() {
    // Particles with non-zero radii and zero LCPO coefficients should still
    // contribute to the LCPO energies of other particles.

    System system;
    system.addParticle(1.0);
    system.addParticle(1.0);
    system.addParticle(1.0);

    LCPOForce* lcpo = new LCPOForce();
    lcpo->setSurfaceTension(5.0);
    lcpo->addParticle(0.22, 10.0, 11.0, 12.0, 13.0);
    lcpo->addParticle(0.24, 0.0, 0.0, 0.0, 0.0);
    lcpo->addParticle(0.26, 0.0, 0.0, 0.0, 0.0);
    system.addForce(lcpo);

    VerletIntegrator integrator(0.001);
    Context context(system, integrator, platform);
    vector<Vec3> positions{
        Vec3(1.1, 1.2, 1.3),
        Vec3(1.4, 1.4, 1.4),
        Vec3(1.3, 1.0, 1.6)
    };
    context.setPositions(positions);
    State state = context.getState(State::Energy | State::Forces);
    const vector<Vec3>& forces = state.getForces();

    Vec3 x12 = positions[1] - positions[0];
    Vec3 x13 = positions[2] - positions[0];
    Vec3 x23 = positions[2] - positions[1];
    double r12 = sqrt(x12.dot(x12));
    double r13 = sqrt(x13.dot(x13));
    double r23 = sqrt(x23.dot(x23));
    double A12 = 2.0 * PI_M * 0.22 * (0.22 - r12 / 2.0 - (0.22 * 0.22 - 0.24 * 0.24) / (2.0 * r12));
    double A13 = 2.0 * PI_M * 0.22 * (0.22 - r13 / 2.0 - (0.22 * 0.22 - 0.26 * 0.26) / (2.0 * r13));
    double A23 = 2.0 * PI_M * 0.24 * (0.24 - r23 / 2.0 - (0.24 * 0.24 - 0.26 * 0.26) / (2.0 * r23));
    double A32 = 2.0 * PI_M * 0.26 * (0.26 - r23 / 2.0 - (0.26 * 0.26 - 0.24 * 0.24) / (2.0 * r23));
    double dA12 = 2.0 * PI_M * 0.22 * (-0.5 + (0.22 * 0.22 - 0.24 * 0.24) / (2.0 * r12 * r12)) / r12;
    double dA13 = 2.0 * PI_M * 0.22 * (-0.5 + (0.22 * 0.22 - 0.26 * 0.26) / (2.0 * r13 * r13)) / r13;
    double dA23 = 2.0 * PI_M * 0.24 * (-0.5 + (0.24 * 0.24 - 0.26 * 0.26) / (2.0 * r23 * r23)) / r23;
    double dA32 = 2.0 * PI_M * 0.26 * (-0.5 + (0.26 * 0.26 - 0.24 * 0.24) / (2.0 * r23 * r23)) / r23;

    ASSERT_EQUAL_TOL(5.0 * (10.0 * 4.0 * PI_M * 0.22 * 0.22 + 11.0 * (A12 + A13) + 12.0 * (A23 + A32) + 13.0 * (A12 * A23 + A13 * A32)), state.getPotentialEnergy(), TOL);
    ASSERT_EQUAL_VEC(5.0 * (11.0 * (dA12 * -x12 + dA13 * -x13) + 13.0 * (dA12 * A23 * -x12 + dA13 * A32 * -x13)), -forces[0], TOL);
    ASSERT_EQUAL_VEC(5.0 * (11.0 * dA12 * x12 + 12.0 * (dA23 * -x23 + dA32 * -x23) + 13.0 * (dA12 * A23 * x12 + A12 * dA23 * -x23 + A13 * dA32 * -x23)), -forces[1], TOL);
    ASSERT_EQUAL_VEC(5.0 * (11.0 * dA13 * x13 + 12.0 * (dA23 * x23 + dA32 * x23) + 13.0 * (A12 * dA23 * x23 + dA13 * A32 * x13 + A13 * dA32 * x23)), -forces[2], TOL);
}

void testExcludeZeroRadius() {
    // Particles with zero radii should be excluded.

    System system;
    system.addParticle(1.0);
    system.addParticle(1.0);
    system.addParticle(1.0);
    system.addParticle(1.0);

    LCPOForce* lcpo = new LCPOForce();
    lcpo->setSurfaceTension(5.0);
    lcpo->addParticle(0.1, 10.0, 11.0, 12.0, 13.0);
    lcpo->addParticle(0.2, 14.0, 15.0, 16.0, 17.0);
    lcpo->addParticle(0.0, 20.0, 20.0, 20.0, 20.0);
    lcpo->addParticle(0.0, 30.0, 30.0, 30.0, 30.0);
    system.addForce(lcpo);

    Vec3 offset(0.48, 0.6, 0.64); // unit

    VerletIntegrator integrator(0.001);
    Context context(system, integrator, platform);

    for (int trial = 0; trial < 3; trial++) {
        double r = 0.15 + 0.1 * trial;
        context.setPositions(vector<Vec3>{Vec3(), r * offset, r * offset * 0.3, r * offset * 0.7});
        State state = context.getState(State::Energy | State::Forces);
        const vector<Vec3>& forces = state.getForces();
        ASSERT_EQUAL_TOL(5.0 * (
            // One-body terms
            10.0 * 4.0 * PI_M * 0.1 * 0.1 +
            14.0 * 4.0 * PI_M * 0.2 * 0.2 +
            // Two-body terms
            (r >= 0.3 ? 0.0 : 11.0 * 2.0 * PI_M * 0.1 * (0.1 - r / 2.0 - (0.1 * 0.1 - 0.2 * 0.2) / (2.0 * r))) +
            (r >= 0.3 ? 0.0 : 15.0 * 2.0 * PI_M * 0.2 * (0.2 - r / 2.0 - (0.2 * 0.2 - 0.1 * 0.1) / (2.0 * r)))),
            state.getPotentialEnergy(), TOL
        );

        double forceExpected = 5.0 * (
            // Two-body terms
            (r >= 0.3 ? 0.0 : 11.0 * 2.0 * PI_M * 0.1 * (-0.5 + (0.1 * 0.1 - 0.2 * 0.2) / (2.0 * r * r))) +
            (r >= 0.3 ? 0.0 : 15.0 * 2.0 * PI_M * 0.2 * (-0.5 + (0.2 * 0.2 - 0.1 * 0.1) / (2.0 * r * r))));
        ASSERT_EQUAL_VEC(forceExpected * offset, forces[0], TOL);
        ASSERT_EQUAL_VEC(forceExpected * -offset, forces[1], TOL);
    }
}

void testAllParticlesExcluded() {
    // If all particles are excluded, LCPOForce should work correctly but give exactly zero for energy and forces.

    System system;
    system.addParticle(1.0);
    system.addParticle(1.0);
    system.addParticle(1.0);
    system.addParticle(1.0);

    LCPOForce* lcpo = new LCPOForce();
    lcpo->addParticle(0.0, 10.0, 11.0, 12.0, 13.0);
    lcpo->addParticle(0.0, 14.0, 15.0, 16.0, 17.0);
    lcpo->addParticle(0.0, 20.0, 20.0, 20.0, 20.0);
    lcpo->addParticle(0.0, 30.0, 30.0, 30.0, 30.0);
    system.addForce(lcpo);

    VerletIntegrator integrator(0.001);
    Context context(system, integrator, platform);
    context.setPositions(vector<Vec3>{Vec3(), Vec3(1.0, 0.0, 0.0), Vec3(0.0, 1.0, 0.0), Vec3(0.0, 0.0, 1.0)});
    State state = context.getState(State::Energy | State::Forces);
    const vector<Vec3>& forces = state.getForces();
    ASSERT_EQUAL(0.0, state.getPotentialEnergy());
    for (int i = 0; i < forces.size(); i++) {
        ASSERT_EQUAL(0.0, forces[i][0]);
        ASSERT_EQUAL(0.0, forces[i][1]);
        ASSERT_EQUAL(0.0, forces[i][2]);
    }
}

void testZeroSurfaceTension() {
    // If the surface tension is zero, LCPOForce should work correctly but give exactly zero for energy and forces.

    System system;
    system.addParticle(1.0);
    system.addParticle(1.0);
    system.addParticle(1.0);
    system.addParticle(1.0);

    LCPOForce* lcpo = new LCPOForce();
    lcpo->setSurfaceTension(0.0);
    lcpo->addParticle(1.0, 10.0, 11.0, 12.0, 13.0);
    lcpo->addParticle(1.0, 14.0, 15.0, 16.0, 17.0);
    lcpo->addParticle(1.0, 20.0, 20.0, 20.0, 20.0);
    lcpo->addParticle(1.0, 30.0, 30.0, 30.0, 30.0);
    system.addForce(lcpo);

    VerletIntegrator integrator(0.001);
    Context context(system, integrator, platform);
    context.setPositions(vector<Vec3>{Vec3(), Vec3(1.0, 0.0, 0.0), Vec3(0.0, 1.0, 0.0), Vec3(0.0, 0.0, 1.0)});
    State state = context.getState(State::Energy | State::Forces);
    const vector<Vec3>& forces = state.getForces();
    ASSERT_EQUAL(0.0, state.getPotentialEnergy());
    for (int i = 0; i < forces.size(); i++) {
        ASSERT_EQUAL(0.0, forces[i][0]);
        ASSERT_EQUAL(0.0, forces[i][1]);
        ASSERT_EQUAL(0.0, forces[i][2]);
    }
}

void makeSmallTestCase(System& system, vector<Vec3>& positions, double& energy) {
    system = System();
    system.addParticle(10.0);
    system.addParticle(10.0);
    system.addParticle(10.0);
    system.addParticle(10.0);
    system.addParticle(10.0);

    LCPOForce* lcpo = new LCPOForce();
    lcpo->setSurfaceTension(5.0);
    lcpo->addParticle(0.6, 10.0, 11.0, 12.0, 13.0);
    lcpo->addParticle(0.65, 14.0, 15.0, 16.0, 17.0);
    lcpo->addParticle(0.7, 18.0, 19.0, 20.0, 21.0);
    lcpo->addParticle(0.75, 22.0, 23.0, 24.0, 25.0);
    lcpo->addParticle(0.8, 26.0, 27.0, 28.0, 29.0);
    system.addForce(lcpo);

    positions = vector<Vec3>{
        Vec3(0.0, 0.3, 0.1),
        Vec3(0.9, 0.1, 0.2),
        Vec3(2.1, 0.2, -0.1),
        Vec3(1.5, 0.8, 0.0),
        Vec3(1.4, -0.8, -0.2)
    };

    double r[5][5], A[5][5];
    for (int i = 0; i < 5; i++) {
        for (int j = 0; j < 5; j++) {
            if (i != j) {
                Vec3 xij = positions[j] - positions[i];
                r[i][j] = sqrt(xij.dot(xij));
                double ri = 0.6 + 0.05 * i;
                double rj = 0.6 + 0.05 * j;
                A[i][j] = 2.0 * PI_M * ri * (ri - r[i][j] / 2.0 - (ri * ri - rj * rj) / (2.0 * r[i][j]));
            }
            else {
                r[i][j] = A[i][j] = 0.0;
            }
        }
    }

    energy = 4.0 * PI_M * (10.0 * 0.6 * 0.6 + 14.0 * 0.65 * 0.65 + 18.0 * 0.7 * 0.7 + 22.0 * 0.75 * 0.75 + 26.0 * 0.8 * 0.8);
    energy += 11.0 * A[0][1] + 15.0 * (A[1][0] + A[1][2] + A[1][3] + A[1][4]) + 19.0 * (A[2][1] + A[2][3] + A[2][4]) + 23.0 * (A[3][1] + A[3][2]) + 27.0 * (A[4][1] + A[4][2]);
    energy += 16.0 * (A[2][3] + A[2][4] + A[3][2] + A[4][2]) + 20.0 * (A[1][3] + A[1][4] + A[3][1] + A[4][1]) + 24.0 * (A[1][2] + A[2][1]) + 28.0 * (A[1][2] + A[2][1]);
    energy += 17.0 * (A[1][2] * (A[2][3] + A[2][4]) + A[1][3] * A[3][2] + A[1][4] * A[4][2]) + 21.0 * (A[2][1] * (A[1][3] + A[1][4]) + A[2][3] * A[3][1] + A[2][4] * A[4][1]) + 25.0 * (A[3][1] * A[1][2] + A[3][2] * A[2][1]) + 29.0 * (A[4][1] * A[1][2] + A[4][2] * A[2][1]);
    energy *= 5.0;
}

void makeAlanineDipeptideTestCase(int n, System& system, vector<Vec3>& positions, double& energy) {
    system = System();
    for (int i = 0; i < n * n * n; i++) {
        system.addParticle(1.008);
        system.addParticle(12.01);
        system.addParticle(1.008);
        system.addParticle(1.008);
        system.addParticle(12.01);
        system.addParticle(16.0);
        system.addParticle(14.01);
        system.addParticle(1.008);
        system.addParticle(12.01);
        system.addParticle(1.008);
        system.addParticle(12.01);
        system.addParticle(1.008);
        system.addParticle(1.008);
        system.addParticle(1.008);
        system.addParticle(12.01);
        system.addParticle(16.0);
        system.addParticle(14.01);
        system.addParticle(1.008);
        system.addParticle(12.01);
        system.addParticle(1.008);
        system.addParticle(1.008);
        system.addParticle(1.008);
    }

    LCPOForce * lcpo = new LCPOForce();
    lcpo->setSurfaceTension(1.0);
    for (int i = 0; i < n * n * n; i++) {
        lcpo->addParticle(0.0, 0.0, 0.0, 0.0, 0.0);
        lcpo->addParticle(0.31, 1.62939604, -0.58707796, -0.0027129056, 0.082274176);
        lcpo->addParticle(0.0, 0.0, 0.0, 0.0, 0.0);
        lcpo->addParticle(0.0, 0.0, 0.0, 0.0, 0.0);
        lcpo->addParticle(0.31, 0.147159648, -0.03977938, -4.6042828e-05, 0.00353025);
        lcpo->addParticle(0.3, 1.43433796, -0.3907856, -0.00283618716, 0.049670356);
        lcpo->addParticle(0.305, 0.85985384, -0.25635368, -0.000157837216, 0.024693968);
        lcpo->addParticle(0.0, 0.0, 0.0, 0.0, 0.0);
        lcpo->addParticle(0.31, 0.48844016, -0.151935684, -0.00042005268, 0.016666964);
        lcpo->addParticle(0.0, 0.0, 0.0, 0.0, 0.0);
        lcpo->addParticle(0.31, 1.62939604, -0.58707796, -0.0027129056, 0.082274176);
        lcpo->addParticle(0.0, 0.0, 0.0, 0.0, 0.0);
        lcpo->addParticle(0.0, 0.0, 0.0, 0.0, 0.0);
        lcpo->addParticle(0.0, 0.0, 0.0, 0.0, 0.0);
        lcpo->addParticle(0.31, 0.147159648, -0.03977938, -4.6042828e-05, 0.00353025);
        lcpo->addParticle(0.3, 1.43433796, -0.3907856, -0.00283618716, 0.049670356);
        lcpo->addParticle(0.305, 0.85985384, -0.25635368, -0.000157837216, 0.024693968);
        lcpo->addParticle(0.0, 0.0, 0.0, 0.0, 0.0);
        lcpo->addParticle(0.31, 1.62939604, -0.58707796, -0.0027129056, 0.082274176);
        lcpo->addParticle(0.0, 0.0, 0.0, 0.0, 0.0);
        lcpo->addParticle(0.0, 0.0, 0.0, 0.0, 0.0);
        lcpo->addParticle(0.0, 0.0, 0.0, 0.0, 0.0);
    }
    system.addForce(lcpo);

    positions.clear();
    for (int x = 0; x < n; x++) {
        for (int y = 0; y < n; y++) {
            for (int z = 0; z < n; z++) {
                Vec3 delta = 2.0 * Vec3(x, y, z);
                positions.push_back(delta + Vec3(0.5734, 0.1803, -0.1291));
                positions.push_back(delta + Vec3(0.6163, 0.1897, -0.0278));
                positions.push_back(delta + Vec3(0.7195, 0.2294, -0.02));
                positions.push_back(delta + Vec3(0.6097, 0.0889, 0.0122));
                positions.push_back(delta + Vec3(0.5296, 0.2862, 0.0558));
                positions.push_back(delta + Vec3(0.4736, 0.2411, 0.1572));
                positions.push_back(delta + Vec3(0.5179, 0.4116, 0.0106));
                positions.push_back(delta + Vec3(0.5847, 0.4262, -0.0666));
                positions.push_back(delta + Vec3(0.4547, 0.5358, 0.0679));
                positions.push_back(delta + Vec3(0.4001, 0.5041, 0.1534));
                positions.push_back(delta + Vec3(0.5605, 0.6416, 0.1053));
                positions.push_back(delta + Vec3(0.5941, 0.6982, 0.0151));
                positions.push_back(delta + Vec3(0.5201, 0.7186, 0.1671));
                positions.push_back(delta + Vec3(0.6541, 0.6078, 0.1553));
                positions.push_back(delta + Vec3(0.3495, 0.6083, -0.026));
                positions.push_back(delta + Vec3(0.2817, 0.7019, 0.0017));
                positions.push_back(delta + Vec3(0.3411, 0.5571, -0.1505));
                positions.push_back(delta + Vec3(0.3787, 0.4676, -0.1632));
                positions.push_back(delta + Vec3(0.2684, 0.6092, -0.264));
                positions.push_back(delta + Vec3(0.3257, 0.5949, -0.3466));
                positions.push_back(delta + Vec3(0.1753, 0.5611, -0.2821));
                positions.push_back(delta + Vec3(0.256, 0.7156, -0.2582));
            }
        }
    }

    // Reference energy value from Amber.
    energy = 1.5511 * 4.184 * n * n * n;
}

void makeGridTestCase(int n, System& system, vector<Vec3>& positions, double& energy) {
    OpenMM_SFMT::SFMT sfmt;
    init_gen_rand(0, sfmt);

    system = System();
    system.setDefaultPeriodicBoxVectors(Vec3(n, 0, 0), Vec3(0, n, 0), Vec3(0, 0, n));
    for (int i = 0; i < n * n * n; i++) {
        system.addParticle(1);
    }

    // Round random variates into {0, 1} to generate radius and position
    // parameters from a discrete set.  Since the forces are discontinuous at
    // the cutoff in LCPO, single precision Platform calculations might not
    // match the ReferencePlatform if any particles are positioned very close
    // to the sum of their radii from each other.  Also round the LCPO
    // parameters since otherwise the GPU context will take too long to try to
    // identify molecules.

    LCPOForce * lcpo = new LCPOForce();
    lcpo->setSurfaceTension(5.0);
    for (int i = 0; i < n * n * n; i++) {
        double ri = 0.6 + 0.2 * round(genrand_real2(sfmt));
        double p1i = round(2 * genrand_real2(sfmt)) / 2;
        double p2i = round(2 * genrand_real2(sfmt)) / 2;
        double p3i = round(2 * genrand_real2(sfmt)) / 2;
        double p4i = round(2 * genrand_real2(sfmt)) / 2;
        lcpo->addParticle(ri, p1i, p2i, p3i, p4i);
    }
    lcpo->setUsesPeriodicBoundaryConditions(true);
    system.addForce(lcpo);

    positions.clear();
    for (int x = 0; x < n; x++) {
        for (int y = 0; y < n; y++) {
            for (int z = 0; z < n; z++) {
                positions.push_back(Vec3(x + 0.1 * round(genrand_real2(sfmt)) - 0.05, y + 0.1 * round(genrand_real2(sfmt)) - 0.05, z + 0.1 * round(genrand_real2(sfmt)) - 0.05));
            }
        }
    }

    vector<double> r, p1, p2, p3, p4;
    for (int i = 0; i < n * n * n; i++) {
        double ri, p1i, p2i, p3i, p4i;
        lcpo->getParticleParameters(i, ri, p1i, p2i, p3i, p4i);
        r.push_back(ri);
        p1.push_back(p1i);
        p2.push_back(p2i);
        p3.push_back(p3i);
        p4.push_back(p4i);
    }

    energy = 0.0;

    for (int xi = 0; xi < n; xi++) {
        for (int yi = 0; yi < n; yi++) {
            for (int zi = 0; zi < n; zi++) {
                int i = (xi * n + yi) * n + zi;
                energy += p1[i] * 4.0 * M_PI * r[i] * r[i];

                // Check the local 3x3x3 grid around this particle for neighbors.
                vector<int> neighbors;
                for (int dx = -1; dx <= 1; dx++) {
                    for (int dy = -1; dy <= 1; dy++) {
                        for (int dz = -1; dz <= 1; dz++) {
                            int xj = xi + dx;
                            int yj = yi + dy;
                            int zj = zi + dz;

                            if (xj < 0) xj += n;
                            if (xj >= n) xj -= n;
                            if (yj < 0) yj += n;
                            if (yj >= n) yj -= n;
                            if (zj < 0) zj += n;
                            if (zj >= n) zj -= n;

                            int j = (xj * n + yj) * n + zj;
                            if (j != i) {
                                neighbors.push_back(j);
                            }
                        }
                    }
                }

                // Do the LCPO calculation manually.
                for (int j : neighbors) {
                    Vec3 xij = positions[j] - positions[i];
                    xij[0] -= n * floor(xij[0] / n + 0.5);
                    xij[1] -= n * floor(xij[1] / n + 0.5);
                    xij[2] -= n * floor(xij[2] / n + 0.5);
                    double rij = sqrt(xij.dot(xij));

                    if (rij >= r[i] + r[j]) {
                        continue;
                    }
                    energy += p2[i] * 2.0 * PI_M * r[i] * (r[i] - rij / 2.0 - (r[i] * r[i] - r[j] * r[j]) / (2.0 * rij));

                    for (int k : neighbors) {
                        if (k == j) {
                            continue;
                        }

                        Vec3 xik = positions[k] - positions[i];
                        Vec3 xjk = positions[k] - positions[j];
                        xik[0] -= n * floor(xik[0] / n + 0.5);
                        xik[1] -= n * floor(xik[1] / n + 0.5);
                        xik[2] -= n * floor(xik[2] / n + 0.5);
                        xjk[0] -= n * floor(xjk[0] / n + 0.5);
                        xjk[1] -= n * floor(xjk[1] / n + 0.5);
                        xjk[2] -= n * floor(xjk[2] / n + 0.5);
                        double rik = sqrt(xik.dot(xik));
                        double rjk = sqrt(xjk.dot(xjk));

                        if (rik >= r[i] + r[k] || rjk >= r[j] + r[k]) {
                            continue;
                        }
                        energy += p3[i] * 2.0 * PI_M * r[j] * (r[j] - rjk / 2.0 - (r[j] * r[j] - r[k] * r[k]) / (2.0 * rjk));
                        energy += p4[i] * 4.0 * PI_M * PI_M * r[i] * r[j] * (r[i] - rij / 2.0 - (r[i] * r[i] - r[j] * r[j]) / (2.0 * rij)) * (r[j] - rjk / 2.0 - (r[j] * r[j] - r[k] * r[k]) / (2.0 * rjk));
                    }
                }
            }
        }
    }

    energy *= 5.0;
}

void runEnergyForcesTestCase(const System& system, vector<Vec3>& positions, double refEnergy, bool doFiniteDifference) {
    VerletIntegrator integrator(0.001);
    Context context(system, integrator, platform);
    context.setPositions(positions);
    State state = context.getState(State::Energy | State::Forces);

    double testEnergy = state.getPotentialEnergy();
    const vector<Vec3>& forces = state.getForces();

    ASSERT_EQUAL_TOL(refEnergy, testEnergy, TOL);

    if (doFiniteDifference) {
        // Check forces with finite differences.
        double delta = 1e-4;
        for (int i = 0; i < forces.size(); i++) {
            for (int d = 0; d < 3; d++) {
                double refPos = positions[i][d];

                positions[i][d] = refPos - delta;
                context.setPositions(positions);
                double energyL = context.getState(State::Energy).getPotentialEnergy();

                positions[i][d] = refPos + delta;
                context.setPositions(positions);
                double energyR = context.getState(State::Energy).getPotentialEnergy();

                positions[i][d] = refPos;

                ASSERT_EQUAL_TOL((energyR - energyL) / (2 * delta), -forces[i][d], TOL);
            }
        }
    }
    else {
        // Calculate forces with the reference platform and compare.
        VerletIntegrator refIntegrator(0.001);
        ReferencePlatform refPlatform;
        Context refContext(system, refIntegrator, refPlatform);
        refContext.setPositions(positions);
        State refState = refContext.getState(State::Forces);
        const vector<Vec3>& refForces = refState.getForces();

        for (int i = 0; i < forces.size(); i++) {
            ASSERT_EQUAL_VEC(refForces[i], forces[i], 1e-3);
        }
    }
}

void testEnergyForces() {
    System system;
    vector<Vec3> positions;
    double energy;

    bool isReferencePlatform = (platform.getName() == "Reference");

    // 5 particles:
    makeSmallTestCase(system, positions, energy);
    runEnergyForcesTestCase(system, positions, energy, isReferencePlatform);

    // 22 particles:
    makeAlanineDipeptideTestCase(1, system, positions, energy);
    runEnergyForcesTestCase(system, positions, energy, isReferencePlatform);

    // 594 particles:
    makeAlanineDipeptideTestCase(3, system, positions, energy);
    runEnergyForcesTestCase(system, positions, energy, isReferencePlatform);

    if (!isReferencePlatform) {
        // 22000 particles:
        makeAlanineDipeptideTestCase(10, system, positions, energy);
        runEnergyForcesTestCase(system, positions, energy, false);

        // 125000 particles:
        makeGridTestCase(50, system, positions, energy);
        runEnergyForcesTestCase(system, positions, energy, false);
    }
}

void testUpdateInContext() {
    System system;
    vector<Vec3> positions;
    double energy;

    makeAlanineDipeptideTestCase(1, system, positions, energy);

    // Get an initial energy and forces.

    VerletIntegrator updateIntegrator(0.001);
    Context updateContext(system, updateIntegrator, platform);
    updateContext.setPositions(positions);
    State initialState = updateContext.getState(State::Energy | State::Forces);

    // Change some parameters on atoms whose radii are set to zero.

    LCPOForce* lcpo = dynamic_cast<LCPOForce*>(&system.getForce(0));
    lcpo->setParticleParameters(3, 0.0, 0.1, 0.2, 0.3, 0.4);
    lcpo->setParticleParameters(9, 0.0, 0.5, 0.6, 0.7, 0.8);
    lcpo->setParticleParameters(11, 0.0, 0.9, 1.0, 1.1, 1.2);
    lcpo->setParticleParameters(12, 0.0, 1.3, 1.4, 1.5, 1.6);
    lcpo->updateParametersInContext(updateContext);

    State updateState1 = updateContext.getState(State::Energy | State::Forces);

    // Change some parameters on atoms involved in the force.

    lcpo->setSurfaceTension(7.0);
    lcpo->setParticleParameters(4, 0.4, 1.7, 1.8, 1.9, 2.0);
    lcpo->setParticleParameters(6, 0.5, 2.1, 2.2, 2.3, 2.4);
    lcpo->setParticleParameters(10, 0.6, 2.5, 2.6, 2.7, 2.8);
    lcpo->updateParametersInContext(updateContext);

    State updateState2 = updateContext.getState(State::Energy | State::Forces);

    // Get an energy and forces from a freshly created context.

    VerletIntegrator referenceIntegrator(0.001);
    Context referenceContext(system, referenceIntegrator, platform);
    referenceContext.setPositions(positions);
    State referenceState = referenceContext.getState(State::Energy | State::Forces);

    // Check energies and forces.

    ASSERT_EQUAL_TOL(initialState.getPotentialEnergy(), updateState1.getPotentialEnergy(), TOL);
    ASSERT_EQUAL_TOL(referenceState.getPotentialEnergy(), updateState2.getPotentialEnergy(), TOL);

    for (int i = 0; i < system.getNumParticles(); i++) {
        ASSERT_EQUAL_VEC(initialState.getForces()[i], updateState1.getForces()[i], TOL);
        ASSERT_EQUAL_VEC(referenceState.getForces()[i], updateState2.getForces()[i], TOL);
    }
}

void testPeriodicShape(bool triclinic) {
    System system;
    vector<Vec3> positions;
    double energy;

    makeAlanineDipeptideTestCase(3, system, positions, energy);

    Vec3 a = Vec3(10.0, 0.0, 0.0);
    Vec3 b = triclinic ? Vec3(-1.0, 9.0, 0.0) : Vec3(0.0, 9.0, 0.0);
    Vec3 c = triclinic ? Vec3(2.0, -3.0, 8.0) : Vec3(0.0, 0.0, 8.0);
    system.setDefaultPeriodicBoxVectors(a, b, c);

    LCPOForce* lcpo = dynamic_cast<LCPOForce*>(&system.getForce(0));
    lcpo->setUsesPeriodicBoundaryConditions(true);

    VerletIntegrator integrator(0.001);
    Context context(system, integrator, platform);

    context.setPositions(positions);
    State referenceState = context.getState(State::Energy | State::Forces);

    // Translate particles at random by some small number of periodic images.

    OpenMM_SFMT::SFMT sfmt;
    init_gen_rand(0, sfmt);
    for (int i = 0; i < system.getNumParticles(); i++) {
        int ia = gen_rand32(sfmt) % 21 - 10;
        int ib = gen_rand32(sfmt) % 21 - 10;
        int ic = gen_rand32(sfmt) % 21 - 10;
        positions[i] += a * (double) ia + b * (double) ib + c * (double) ic;
    }

    context.setPositions(positions);
    State testState = context.getState(State::Energy | State::Forces);

    // Check energies and forces.

    ASSERT_EQUAL_TOL(referenceState.getPotentialEnergy(), testState.getPotentialEnergy(), TOL);

    for (int i = 0; i < system.getNumParticles(); i++) {
        ASSERT_EQUAL_VEC(referenceState.getForces()[i], testState.getForces()[i], 1e-4);
    }
}

void runPlatformTests();

int main(int argc, char* argv[]) {
    try {
        initializeTests(argc, argv);

        testSingleParticle();
        testTwoParticles();
        testThreeParticles();
        testUsePeriodic(false);
        testUsePeriodic(true);
        testIncludeZeroScale();
        testExcludeZeroRadius();
        testAllParticlesExcluded();
        testZeroSurfaceTension();
        testEnergyForces();
        testUpdateInContext();
        testPeriodicShape(false);
        testPeriodicShape(true);
        runPlatformTests();
    }
    catch(const exception& e) {
        cout << "exception: " << e.what() << endl;
        return 1;
    }
    cout << "Done" << endl;
    return 0;
}
