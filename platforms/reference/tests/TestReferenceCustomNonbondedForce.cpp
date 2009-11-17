
/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2008-2009 Stanford University and the Authors.      *
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
#include "../src/sfmt/SFMT.h"
#include "openmm/Context.h"
#include "ReferencePlatform.h"
#include "openmm/CustomNonbondedForce.h"
#include "openmm/NonbondedForce.h"
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
    CustomNonbondedForce* forceField = new CustomNonbondedForce("scale*a*(r*b)^3; a=a1*a2; b=c+b1+b2");
    forceField->addPerParticleParameter("a");
    forceField->addPerParticleParameter("b");
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

void testExclusions() {
    ReferencePlatform platform;
    System system;
    VerletIntegrator integrator(0.01);
    CustomNonbondedForce* nonbonded = new CustomNonbondedForce("a*r; a=a1+a2");
    nonbonded->addPerParticleParameter("a");
    vector<double> params(1);
    vector<Vec3> positions(4);
    for (int i = 0; i < 4; i++) {
        system.addParticle(1.0);
        params[0] = i+1;
        nonbonded->addParticle(params);
        positions[i] = Vec3(i, 0, 0);
    }
    nonbonded->addExclusion(0, 1);
    nonbonded->addExclusion(1, 2);
    nonbonded->addExclusion(2, 3);
    nonbonded->addExclusion(0, 2);
    nonbonded->addExclusion(1, 3);
    system.addForce(nonbonded);
    Context context(system, integrator, platform);
    context.setPositions(positions);
    State state = context.getState(State::Forces | State::Energy);
    const vector<Vec3>& forces = state.getForces();
    ASSERT_EQUAL_VEC(Vec3(1+4, 0, 0), forces[0], TOL);
    ASSERT_EQUAL_VEC(Vec3(0, 0, 0), forces[1], TOL);
    ASSERT_EQUAL_VEC(Vec3(0, 0, 0), forces[2], TOL);
    ASSERT_EQUAL_VEC(Vec3(-(1+4), 0, 0), forces[3], TOL);
    ASSERT_EQUAL_TOL((1+4)*3, state.getPotentialEnergy(), TOL);
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
    system.setPeriodicBoxVectors(Vec3(4, 0, 0), Vec3(0, 4, 0), Vec3(0, 0, 4));
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

void testTabulatedFunction(bool interpolating) {
    ReferencePlatform platform;
    System system;
    system.addParticle(1.0);
    system.addParticle(1.0);
    VerletIntegrator integrator(0.01);
    CustomNonbondedForce* forceField = new CustomNonbondedForce("fn(r)+1");
    forceField->addParticle(vector<double>());
    forceField->addParticle(vector<double>());
    vector<double> table;
    for (int i = 0; i < 21; i++)
        table.push_back(std::sin(0.25*i));
    forceField->addFunction("fn", table, 1.0, 6.0, interpolating);
    system.addForce(forceField);
    Context context(system, integrator, platform);
    vector<Vec3> positions(2);
    positions[0] = Vec3(0, 0, 0);
    double tol = 0.01;
    for (int i = 1; i < 30; i++) {
        double x = (7.0/30.0)*i;
        positions[1] = Vec3(x, 0, 0);
        context.setPositions(positions);
        State state = context.getState(State::Forces | State::Energy);
        const vector<Vec3>& forces = state.getForces();
        double force = (x < 1.0 || x > 6.0 ? 0.0 : -std::cos(x-1.0));
        double energy = (x < 1.0 || x > 6.0 ? 0.0 : std::sin(x-1.0))+1.0;
        ASSERT_EQUAL_VEC(Vec3(-force, 0, 0), forces[0], 0.1);
        ASSERT_EQUAL_VEC(Vec3(force, 0, 0), forces[1], 0.1);
        ASSERT_EQUAL_TOL(energy, state.getPotentialEnergy(), 0.02);
    }
}

void testCoulombLennardJones() {
    const int numMolecules = 300;
    const int numParticles = numMolecules*2;
    const double boxSize = 20.0;
    ReferencePlatform platform;

    // Create two systems: one with a NonbondedForce, and one using a CustomNonbondedForce to implement the same interaction.

    System standardSystem;
    System customSystem;
    for (int i = 0; i < numParticles; i++) {
        standardSystem.addParticle(1.0);
        customSystem.addParticle(1.0);
    }
    NonbondedForce* standardNonbonded = new NonbondedForce();
    CustomNonbondedForce* customNonbonded = new CustomNonbondedForce("4*eps*((sigma/r)^12-(sigma/r)^6)+138.935485*q/r; q=q1*q2; sigma=0.5*(sigma1+sigma2); eps=sqrt(eps1*eps2)");
    customNonbonded->addPerParticleParameter("q");
    customNonbonded->addPerParticleParameter("sigma");
    customNonbonded->addPerParticleParameter("eps");
    vector<Vec3> positions(numParticles);
    vector<Vec3> velocities(numParticles);
    init_gen_rand(0);
    vector<double> params(3);
    for (int i = 0; i < numMolecules; i++) {
        if (i < numMolecules/2) {
            standardNonbonded->addParticle(1.0, 0.2, 0.1);
            params[0] = 1.0;
            params[1] = 0.2;
            params[2] = 0.1;
            customNonbonded->addParticle(params);
            standardNonbonded->addParticle(-1.0, 0.1, 0.1);
            params[0] = -1.0;
            params[1] = 0.1;
            customNonbonded->addParticle(params);
        }
        else {
            standardNonbonded->addParticle(1.0, 0.2, 0.2);
            params[0] = 1.0;
            params[1] = 0.2;
            params[2] = 0.2;
            customNonbonded->addParticle(params);
            standardNonbonded->addParticle(-1.0, 0.1, 0.2);
            params[0] = -1.0;
            params[1] = 0.1;
            customNonbonded->addParticle(params);
        }
        positions[2*i] = Vec3(boxSize*genrand_real2(), boxSize*genrand_real2(), boxSize*genrand_real2());
        positions[2*i+1] = Vec3(positions[2*i][0]+1.0, positions[2*i][1], positions[2*i][2]);
        velocities[2*i] = Vec3(genrand_real2(), genrand_real2(), genrand_real2());
        velocities[2*i+1] = Vec3(genrand_real2(), genrand_real2(), genrand_real2());
        standardNonbonded->addException(2*i, 2*i+1, 0.0, 1.0, 0.0);
        customNonbonded->addExclusion(2*i, 2*i+1);
    }
    standardNonbonded->setNonbondedMethod(NonbondedForce::NoCutoff);
    customNonbonded->setNonbondedMethod(CustomNonbondedForce::NoCutoff);
    standardSystem.addForce(standardNonbonded);
    customSystem.addForce(customNonbonded);
    VerletIntegrator integrator1(0.01);
    VerletIntegrator integrator2(0.01);
    Context context1(standardSystem, integrator1, platform);
    Context context2(customSystem, integrator2, platform);
    context1.setPositions(positions);
    context2.setPositions(positions);
    context1.setVelocities(velocities);
    context2.setVelocities(velocities);
    State state1 = context1.getState(State::Forces | State::Energy);
    State state2 = context2.getState(State::Forces | State::Energy);
    ASSERT_EQUAL_TOL(state1.getPotentialEnergy(), state2.getPotentialEnergy(), 1e-4);
    for (int i = 0; i < numParticles; i++) {
        ASSERT_EQUAL_VEC(state1.getForces()[i], state2.getForces()[i], 1e-4);
    }
}

int main() {
    try {
        testSimpleExpression();
        testParameters();
        testExclusions();
        testCutoff();
        testPeriodic();
        testTabulatedFunction(true);
        testTabulatedFunction(false);
        testCoulombLennardJones();
    }
    catch(const exception& e) {
        cout << "exception: " << e.what() << endl;
        return 1;
    }
    cout << "Done" << endl;
    return 0;
}
