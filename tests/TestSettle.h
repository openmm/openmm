
/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit.                   *
 * See https://openmm.org/development.                                        *
 *                                                                            *
 * Portions copyright (c) 2008-2015 Stanford University and the Authors.      *
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
#include "openmm/NonbondedForce.h"
#include "openmm/OpenMMException.h"
#include "openmm/System.h"
#include "openmm/LangevinIntegrator.h"
#include "sfmt/SFMT.h"
#include <cmath>
#include <iostream>
#include <limits>
#include <string>
#include <vector>

using namespace OpenMM;
using namespace std;

struct TriangleExecutionResult {
    bool contextRejected;
    bool positionsFinite;
    bool velocitiesFinite;
};

bool areFinite(const vector<Vec3>& values) {
    for (const Vec3& value : values)
        for (int component = 0; component < 3; component++)
            if (!isfinite(value[component]))
                return false;
    return true;
}

TriangleExecutionResult executeTriangleAtConstraintSeam(double distance12) {
    TriangleExecutionResult result = {false, true, true};
    System system;
    for (int i = 0; i < 3; i++)
        system.addParticle(1.0);
    system.addConstraint(0, 1, 1.0);
    system.addConstraint(0, 2, 1.0);
    system.addConstraint(1, 2, distance12);
    LangevinIntegrator integrator(100.0, 2.0, 0.001);
    Context* context = NULL;
    try {
        context = new Context(system, integrator, platform);
    }
    catch (const OpenMMException& exception) {
        result.contextRejected = true;
        ASSERT(string(exception.what()).find("degenerate, inconsistent, or numerically unsafe triangle") != string::npos);
    }
    if (context != NULL) {
        vector<Vec3> positions = {Vec3(0, 0, 0), Vec3(1, 0, 0), Vec3(0, 1, 0)};
        vector<Vec3> velocities = {Vec3(0.1, 0.2, 0.3), Vec3(-0.2, 0.1, -0.1), Vec3(0.3, -0.1, 0.2)};
        context->setPositions(positions);
        context->setVelocities(velocities);
        context->applyConstraints(1e-5);
        result.positionsFinite = areFinite(context->getState(State::Positions).getPositions());
        context->applyVelocityConstraints(1e-5);
        result.velocitiesFinite = areFinite(context->getState(State::Velocities).getVelocities());
        delete context;
    }
    return result;
}

void testDegenerateTriangleExecution() {
    const double previousFloat = nextafterf(2.0f, 0.0f);
    const double roundingBoundary = 0.5*(2.0+previousFloat);
    const double lastRoundsDown = nextafter(roundingBoundary, 0.0);
    const double firstRoundsToCollinear = nextafter(roundingBoundary, 2.0);
    const double nextBelowTwo = nextafter(2.0, 0.0);
    const double nextAboveTwo = nextafter(2.0, numeric_limits<double>::infinity());
    ASSERT((float) lastRoundsDown == (float) previousFloat);
    ASSERT((float) firstRoundsToCollinear == 2.0f);
    ASSERT((float) nextBelowTwo == 2.0f);
    ASSERT((float) nextAboveTwo == 2.0f);

    TriangleExecutionResult exact = executeTriangleAtConstraintSeam(2.0);
    TriangleExecutionResult nextBelow = executeTriangleAtConstraintSeam(nextBelowTwo);
    TriangleExecutionResult nextAbove = executeTriangleAtConstraintSeam(nextAboveTwo);
    TriangleExecutionResult firstRounded = executeTriangleAtConstraintSeam(firstRoundsToCollinear);
    TriangleExecutionResult nondegenerate = executeTriangleAtConstraintSeam(lastRoundsDown);

    ASSERT(exact.contextRejected);
    ASSERT(nextBelow.contextRejected);
    ASSERT(nextAbove.contextRejected);
    ASSERT(firstRounded.contextRejected);
    ASSERT(!nondegenerate.contextRejected);
    ASSERT(nondegenerate.positionsFinite);
    ASSERT(nondegenerate.velocitiesFinite);
}

void testEmbeddedNearDegenerateTriangle() {
    const double previousFloat = nextafterf(2.0f, 0.0f);
    const double roundingBoundary = 0.5*(2.0+previousFloat);
    const double firstRoundsToCollinear = nextafter(roundingBoundary, 2.0);
    System system;
    for (int i = 0; i < 4; i++)
        system.addParticle(1.0);
    system.addConstraint(0, 1, 1.0);
    system.addConstraint(0, 2, 1.0);
    system.addConstraint(1, 2, firstRoundsToCollinear);
    system.addConstraint(0, 3, 1.0);
    LangevinIntegrator integrator(100.0, 2.0, 0.001);
    Context context(system, integrator, platform);
    const double halfBase = 0.5*firstRoundsToCollinear;
    const double height = sqrt(1.0-halfBase*halfBase);
    vector<Vec3> positions = {Vec3(0, height, 0), Vec3(-halfBase, 0, 0),
            Vec3(halfBase, 0, 0), Vec3(0, height, 1)};
    vector<Vec3> velocities = {Vec3(0.1, 0.2, 0.3), Vec3(-0.2, 0.1, -0.1),
            Vec3(0.3, -0.1, 0.2), Vec3(-0.1, 0.2, -0.2)};
    context.setPositions(positions);
    context.setVelocities(velocities);
    context.applyConstraints(1e-5);
    ASSERT(areFinite(context.getState(State::Positions).getPositions()));
    context.applyVelocityConstraints(1e-5);
    ASSERT(areFinite(context.getState(State::Velocities).getVelocities()));
}

void assertTriangleRejected(double distance01, double distance02, double distance12) {
    System system;
    for (int i = 0; i < 3; i++)
        system.addParticle(1.0);
    system.addConstraint(0, 1, distance01);
    system.addConstraint(0, 2, distance02);
    system.addConstraint(1, 2, distance12);
    LangevinIntegrator integrator(100.0, 2.0, 0.001);
    bool threwException = false;
    try {
        Context context(system, integrator, platform);
    }
    catch (const OpenMMException& exception) {
        threwException = true;
        ASSERT(string(exception.what()).find("degenerate, inconsistent, or numerically unsafe triangle") != string::npos);
    }
    ASSERT(threwException);
}

void testImpossibleTriangle() {
    assertTriangleRejected(1.0, 1.0, 3.0);
    assertTriangleRejected(1.0, 2.0, 4.0);
    assertTriangleRejected(1.0e-20, 1.0e-20, 3.0e-20);
    assertTriangleRejected(1.0e20, 1.0e20, 3.0e20);
    const double denorm = numeric_limits<double>::denorm_min();
    const double maximum = numeric_limits<double>::max();
    assertTriangleRejected(denorm, denorm, 3*denorm);
    assertTriangleRejected(maximum/3, maximum/3, maximum);
    const double scale = 1.0e20;
    const double epsilon = numeric_limits<double>::epsilon();
    assertTriangleRejected(scale, scale, 2*scale*(1+16*epsilon));
}

void testSettleArithmeticBoundaries() {
    const float largeArm = 1.0e20f;
    const float largeHalfBase = 0.5f*largeArm;
    const float largeRadicand = largeArm*largeArm-largeHalfBase*largeHalfBase;
    ASSERT(!isfinite(largeRadicand) || largeRadicand <= 0);
    assertTriangleRejected(1.0e20, 1.0e20, 1.0e20);

    const float subnormalArm = 1.0e-19f;
    const float subnormalHalfBase = 0.5f*subnormalArm;
    const float subnormalRadicand = subnormalArm*subnormalArm-subnormalHalfBase*subnormalHalfBase;
    ASSERT(subnormalRadicand > 0);
    ASSERT(!std::isnormal(subnormalRadicand));
    assertTriangleRejected(1.0e-19, 1.0e-19, 1.0e-19);

    const float smallArm = 1.0e-30f;
    const float smallHalfBase = 0.5f*smallArm;
    const float smallRadicand = smallArm*smallArm-smallHalfBase*smallHalfBase;
    ASSERT(!isfinite(smallRadicand) || smallRadicand <= 0);
    assertTriangleRejected(1.0e-30, 1.0e-30, 1.0e-30);
}

void testEmbeddedImpossibleTriangle() {
    System system;
    for (int i = 0; i < 6; i++)
        system.addParticle(1.0);
    system.addConstraint(4, 1, 1.0);
    system.addConstraint(5, 1, 1.0);
    system.addConstraint(5, 4, 3.0);
    system.addConstraint(4, 2, 1.0);
    LangevinIntegrator integrator(100.0, 2.0, 0.001);
    bool threwException = false;
    try {
        Context context(system, integrator, platform);
    }
    catch (const OpenMMException& exception) {
        threwException = true;
        ASSERT(string(exception.what()).find("degenerate, inconsistent, or numerically unsafe triangle") != string::npos);
    }
    ASSERT(threwException);
}

void assertTriangleAccepted(double distance01, double distance02, double distance12) {
    System system;
    for (int i = 0; i < 3; i++)
        system.addParticle(1.0);
    system.addConstraint(0, 1, distance01);
    system.addConstraint(0, 2, distance02);
    system.addConstraint(1, 2, distance12);
    LangevinIntegrator integrator(100.0, 2.0, 0.001);
    Context context(system, integrator, platform);
}

void testTriangleBoundaries() {
    assertTriangleAccepted(3.0, 4.0, 5.0);
}

void testPartialConstraintGraph() {
    System system;
    for (int i = 0; i < 6; i++)
        system.addParticle(1.0);
    for (int particle1 = 0; particle1 < 3; particle1++)
        for (int particle2 = 3; particle2 < 6; particle2++)
            system.addConstraint(particle1, particle2, 1.0);
    LangevinIntegrator integrator(100.0, 2.0, 0.001);
    Context context(system, integrator, platform);
}

void assertInvalidConstraintDistance(double distance) {
    System system;
    system.addParticle(1.0);
    system.addParticle(1.0);
    system.addConstraint(0, 1, distance);
    LangevinIntegrator integrator(100.0, 2.0, 0.001);
    bool threwException = false;
    try {
        Context context(system, integrator, platform);
    }
    catch (const OpenMMException& exception) {
        threwException = true;
        ASSERT(string(exception.what()).find("finite and positive") != string::npos);
    }
    ASSERT(threwException);
}

void testInvalidConstraintDistances() {
    assertInvalidConstraintDistance(0.0);
    assertInvalidConstraintDistance(-1.0);
    assertInvalidConstraintDistance(numeric_limits<double>::infinity());
    assertInvalidConstraintDistance(numeric_limits<double>::quiet_NaN());
}

void testConstraints() {
    const int numMolecules = 10;
    const int numParticles = numMolecules*3;
    const int numConstraints = numMolecules*3;
    const double temp = 100.0;
    System system;
    LangevinIntegrator integrator(temp, 2.0, 0.001);
    integrator.setConstraintTolerance(1e-5);
    NonbondedForce* forceField = new NonbondedForce();
    for (int i = 0; i < numMolecules; ++i) {
        system.addParticle(16.0);
        system.addParticle(1.0);
        system.addParticle(1.0);
        forceField->addParticle(-0.82, 0.317, 0.65);
        forceField->addParticle(0.41, 1.0, 0.0);
        forceField->addParticle(0.41, 1.0, 0.0);
        system.addConstraint(i*3, i*3+1, 0.1);
        system.addConstraint(i*3, i*3+2, 0.1);
        system.addConstraint(i*3+1, i*3+2, 0.163);
    }
    system.addForce(forceField);
    Context context(system, integrator, platform);
    vector<Vec3> positions(numParticles);
    vector<Vec3> velocities(numParticles);
    OpenMM_SFMT::SFMT sfmt;
    init_gen_rand(0, sfmt);

    for (int i = 0; i < numMolecules; ++i) {
        positions[i*3] = Vec3((i%4)*0.4, (i/4)*0.4, 0);
        positions[i*3+1] = positions[i*3]+Vec3(0.1, 0, 0);
        positions[i*3+2] = positions[i*3]+Vec3(-0.03333, 0.09428, 0);
        velocities[i*3] = Vec3(genrand_real2(sfmt)-0.5, genrand_real2(sfmt)-0.5, genrand_real2(sfmt)-0.5);
        velocities[i*3+1] = Vec3(genrand_real2(sfmt)-0.5, genrand_real2(sfmt)-0.5, genrand_real2(sfmt)-0.5);
        velocities[i*3+2] = Vec3(genrand_real2(sfmt)-0.5, genrand_real2(sfmt)-0.5, genrand_real2(sfmt)-0.5);
    }
    context.setPositions(positions);
    context.setVelocities(velocities);

    // Simulate it and see whether the constraints remain satisfied.

    for (int i = 0; i < 1000; ++i) {
        integrator.step(1);
        State state = context.getState(State::Positions | State::Forces);
        for (int j = 0; j < numConstraints; ++j) {
            int particle1, particle2;
            double distance;
            system.getConstraintParameters(j, particle1, particle2, distance);
            Vec3 p1 = state.getPositions()[particle1];
            Vec3 p2 = state.getPositions()[particle2];
            double dist = std::sqrt((p1[0]-p2[0])*(p1[0]-p2[0])+(p1[1]-p2[1])*(p1[1]-p2[1])+(p1[2]-p2[2])*(p1[2]-p2[2]));
            ASSERT_EQUAL_TOL(distance, dist, 1e-5);
        }
    }
}

void runPlatformTests();

int main(int argc, char* argv[]) {
    try {
        initializeTests(argc, argv);
        testImpossibleTriangle();
        testSettleArithmeticBoundaries();
        testEmbeddedImpossibleTriangle();
        testTriangleBoundaries();
        testDegenerateTriangleExecution();
        testEmbeddedNearDegenerateTriangle();
        testPartialConstraintGraph();
        testInvalidConstraintDistances();
        testConstraints();
        runPlatformTests();
    }
    catch(const exception& e) {
        cout << "exception: " << e.what() << endl;
        return 1;
    }
    cout << "Done" << endl;
    return 0;
}
