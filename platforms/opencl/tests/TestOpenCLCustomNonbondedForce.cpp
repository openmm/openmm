
/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2008-2013 Stanford University and the Authors.      *
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
 * This tests all the different force terms in the OpenCL implementation of CustomNonbondedForce.
 */

#include "openmm/internal/AssertionUtilities.h"
#include "sfmt/SFMT.h"
#include "openmm/Context.h"
#include "OpenCLPlatform.h"
#include "openmm/CustomNonbondedForce.h"
#include "openmm/NonbondedForce.h"
#include "openmm/System.h"
#include "openmm/VerletIntegrator.h"
#include <iostream>
#include <vector>

using namespace OpenMM;
using namespace std;

static OpenCLPlatform platform;

const double TOL = 1e-5;

void testSimpleExpression() {
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
    
    // Try changing the global parameters and make sure it's still correct.
    
    context.setParameter("scale", 1.5);
    context.setParameter("c", 1.0);
    state = context.getState(State::Forces | State::Energy);
    forces = state.getForces();
    force = -1.5*3.0*3*6.0*(12*12);
    ASSERT_EQUAL_VEC(Vec3(-force, 0, 0), forces[0], TOL);
    ASSERT_EQUAL_VEC(Vec3(force, 0, 0), forces[1], TOL);
    ASSERT_EQUAL_TOL(1.5*3.0*(12*12*12), state.getPotentialEnergy(), TOL);
    
    // Try changing the per-particle parameters and make sure it's still correct.
    
    params[0] = 1.6;
    params[1] = 2.1;
    forceField->setParticleParameters(0, params);
    params[0] = 1.9;
    params[1] = 2.8;
    forceField->setParticleParameters(1, params);
    forceField->updateParametersInContext(context);
    state = context.getState(State::Forces | State::Energy);
    forces = state.getForces();
    force = -1.5*1.6*1.9*3*5.9*(11.8*11.8);
    ASSERT_EQUAL_VEC(Vec3(-force, 0, 0), forces[0], TOL);
    ASSERT_EQUAL_VEC(Vec3(force, 0, 0), forces[1], TOL);
    ASSERT_EQUAL_TOL(1.5*1.6*1.9*(11.8*11.8*11.8), state.getPotentialEnergy(), TOL);
}

void testManyParameters() {
    System system;
    system.addParticle(1.0);
    system.addParticle(1.0);
    VerletIntegrator integrator(0.01);
    CustomNonbondedForce* forceField = new CustomNonbondedForce("(a1*a2+b1*b2+c1*c2+d1*d2+e1*e2)*r");
    forceField->addPerParticleParameter("a");
    forceField->addPerParticleParameter("b");
    forceField->addPerParticleParameter("c");
    forceField->addPerParticleParameter("d");
    forceField->addPerParticleParameter("e");
    vector<double> params(5);
    params[0] = 1.0;
    params[1] = 2.0;
    params[2] = 3.0;
    params[3] = 4.0;
    params[4] = 5.0;
    forceField->addParticle(params);
    params[0] = 1.1;
    params[1] = 1.2;
    params[2] = 1.3;
    params[3] = 1.4;
    params[4] = 1.5;
    forceField->addParticle(params);
    system.addForce(forceField);
    Context context(system, integrator, platform);
    vector<Vec3> positions(2);
    positions[0] = Vec3(0, 0, 0);
    positions[1] = Vec3(2, 0, 0);
    context.setPositions(positions);
    State state = context.getState(State::Forces | State::Energy);
    vector<Vec3> forces = state.getForces();
    double force = 1*1.1 + 2*1.2 + 3*1.3 + 4*1.4 + 5*1.5;
    ASSERT_EQUAL_VEC(Vec3(force, 0, 0), forces[0], TOL);
    ASSERT_EQUAL_VEC(Vec3(-force, 0, 0), forces[1], TOL);
    ASSERT_EQUAL_TOL(2*force, state.getPotentialEnergy(), TOL);
}

void testExclusions() {
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
    ASSERT_EQUAL_TOL((1+4)*3.0, state.getPotentialEnergy(), TOL);
}

void testCutoff() {
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
    system.setDefaultPeriodicBoxVectors(Vec3(4, 0, 0), Vec3(0, 4, 0), Vec3(0, 0, 4));
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

void testTabulatedFunction() {
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
    forceField->addFunction("fn", table, 1.0, 6.0);
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
    for (int i = 1; i < 20; i++) {
        double x = 0.25*i+1.0;
        positions[1] = Vec3(x, 0, 0);
        context.setPositions(positions);
        State state = context.getState(State::Energy);
        double energy = (x < 1.0 || x > 6.0 ? 0.0 : std::sin(x-1.0))+1.0;
        ASSERT_EQUAL_TOL(energy, state.getPotentialEnergy(), 1e-4);
    }
}

void testCoulombLennardJones() {
    const int numMolecules = 300;
    const int numParticles = numMolecules*2;
    const double boxSize = 20.0;

    // Create two systems: one with a NonbondedForce, and one using a CustomNonbondedForce to implement the same interaction.

    System standardSystem;
    System customSystem;
    for (int i = 0; i < numParticles; i++) {
        standardSystem.addParticle(1.0);
        customSystem.addParticle(1.0);
    }
    NonbondedForce* standardNonbonded = new NonbondedForce();
    CustomNonbondedForce* customNonbonded = new CustomNonbondedForce("4*eps*((sigma/r)^12-(sigma/r)^6)+138.935456*q/r; q=q1*q2; sigma=0.5*(sigma1+sigma2); eps=sqrt(eps1*eps2)");
    customNonbonded->addPerParticleParameter("q");
    customNonbonded->addPerParticleParameter("sigma");
    customNonbonded->addPerParticleParameter("eps");
    vector<Vec3> positions(numParticles);
    vector<Vec3> velocities(numParticles);
    OpenMM_SFMT::SFMT sfmt;
    init_gen_rand(0, sfmt);

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
        positions[2*i] = Vec3(boxSize*genrand_real2(sfmt), boxSize*genrand_real2(sfmt), boxSize*genrand_real2(sfmt));
        positions[2*i+1] = Vec3(positions[2*i][0]+1.0, positions[2*i][1], positions[2*i][2]);
        velocities[2*i] = Vec3(genrand_real2(sfmt), genrand_real2(sfmt), genrand_real2(sfmt));
        velocities[2*i+1] = Vec3(genrand_real2(sfmt), genrand_real2(sfmt), genrand_real2(sfmt));
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

void testParallelComputation() {
    System system;
    const int numParticles = 200;
    for (int i = 0; i < numParticles; i++)
        system.addParticle(1.0);
    CustomNonbondedForce* force = new CustomNonbondedForce("4*eps*((sigma/r)^12-(sigma/r)^6); sigma=0.5; eps=1");
    vector<double> params;
    for (int i = 0; i < numParticles; i++)
        force->addParticle(params);
    system.addForce(force);
    OpenMM_SFMT::SFMT sfmt;
    init_gen_rand(0, sfmt);
    vector<Vec3> positions(numParticles);
    for (int i = 0; i < numParticles; i++)
        positions[i] = Vec3(5*genrand_real2(sfmt), 5*genrand_real2(sfmt), 5*genrand_real2(sfmt));
    for (int i = 0; i < numParticles; ++i)
        for (int j = 0; j < i; ++j) {
            Vec3 delta = positions[i]-positions[j];
            if (delta.dot(delta) < 0.1)
                force->addExclusion(i, j);
        }
    VerletIntegrator integrator1(0.01);
    Context context1(system, integrator1, platform);
    context1.setPositions(positions);
    State state1 = context1.getState(State::Forces | State::Energy);
    VerletIntegrator integrator2(0.01);
    string deviceIndex = platform.getPropertyValue(context1, OpenCLPlatform::OpenCLDeviceIndex());
    map<string, string> props;
    props[OpenCLPlatform::OpenCLDeviceIndex()] = deviceIndex+","+deviceIndex;
    Context context2(system, integrator2, platform, props);
    context2.setPositions(positions);
    State state2 = context2.getState(State::Forces | State::Energy);
    ASSERT_EQUAL_TOL(state1.getPotentialEnergy(), state2.getPotentialEnergy(), 1e-5);
    for (int i = 0; i < numParticles; i++)
        ASSERT_EQUAL_VEC(state1.getForces()[i], state2.getForces()[i], 1e-5);
}

void testSwitchingFunction() {
    System system;
    system.addParticle(1.0);
    system.addParticle(1.0);
    VerletIntegrator integrator(0.01);
    CustomNonbondedForce* nonbonded = new CustomNonbondedForce("10/r^2");
    vector<double> params;
    nonbonded->addParticle(params);
    nonbonded->addParticle(params);
    nonbonded->setNonbondedMethod(CustomNonbondedForce::CutoffNonPeriodic);
    nonbonded->setCutoffDistance(2.0);
    nonbonded->setUseSwitchingFunction(true);
    nonbonded->setSwitchingDistance(1.5);
    system.addForce(nonbonded);
    Context context(system, integrator, platform);
    vector<Vec3> positions(2);
    positions[0] = Vec3(0, 0, 0);
    
    // Compute the interaction at various distances.
    
    for (double r = 1.0; r < 2.5; r += 0.1) {
        positions[1] = Vec3(r, 0, 0);
        context.setPositions(positions);
        State state = context.getState(State::Forces | State::Energy);
        
        // See if the energy is correct.
        
        double expectedEnergy = 10/(r*r);
        double switchValue;
        if (r <= 1.5)
            switchValue = 1;
        else if (r >= 2.0)
            switchValue = 0;
        else {
            double t = (r-1.5)/0.5;
            switchValue = 1+t*t*t*(-10+t*(15-t*6));
        }
        ASSERT_EQUAL_TOL(switchValue*expectedEnergy, state.getPotentialEnergy(), TOL);
        
        // See if the force is the gradient of the energy.
        
        double delta = 1e-3;
        positions[1] = Vec3(r-delta, 0, 0);
        context.setPositions(positions);
        double e1 = context.getState(State::Energy).getPotentialEnergy();
        positions[1] = Vec3(r+delta, 0, 0);
        context.setPositions(positions);
        double e2 = context.getState(State::Energy).getPotentialEnergy();
        ASSERT_EQUAL_TOL((e2-e1)/(2*delta), state.getForces()[0][0], 2e-3);
    }
}

void testLongRangeCorrection() {
    // Create a box of particles.

    int gridSize = 5;
    int numParticles = gridSize*gridSize*gridSize;
    double boxSize = gridSize*0.7;
    double cutoff = boxSize/3;
    System standardSystem;
    System customSystem;
    VerletIntegrator integrator1(0.01);
    VerletIntegrator integrator2(0.01);
    NonbondedForce* standardNonbonded = new NonbondedForce();
    CustomNonbondedForce* customNonbonded = new CustomNonbondedForce("4*eps*((sigma/r)^12-(sigma/r)^6); sigma=0.5*(sigma1+sigma2); eps=sqrt(eps1*eps2)");
    customNonbonded->addPerParticleParameter("sigma");
    customNonbonded->addPerParticleParameter("eps");
    vector<Vec3> positions(numParticles);
    int index = 0;
    vector<double> params1(2);
    params1[0] = 1.1;
    params1[1] = 0.5;
    vector<double> params2(2);
    params2[0] = 1;
    params2[1] = 1;
    for (int i = 0; i < gridSize; i++)
        for (int j = 0; j < gridSize; j++)
            for (int k = 0; k < gridSize; k++) {
                standardSystem.addParticle(1.0);
                customSystem.addParticle(1.0);
                if (index%2 == 0) {
                    standardNonbonded->addParticle(0, params1[0], params1[1]);
                    customNonbonded->addParticle(params1);
                }
                else {
                    standardNonbonded->addParticle(0, params2[0], params2[1]);
                    customNonbonded->addParticle(params2);
                }
                positions[index] = Vec3(i*boxSize/gridSize, j*boxSize/gridSize, k*boxSize/gridSize);
                index++;
            }
    standardNonbonded->setNonbondedMethod(NonbondedForce::CutoffPeriodic);
    customNonbonded->setNonbondedMethod(CustomNonbondedForce::CutoffPeriodic);
    standardNonbonded->setCutoffDistance(cutoff);
    customNonbonded->setCutoffDistance(cutoff);
    standardSystem.setDefaultPeriodicBoxVectors(Vec3(boxSize, 0, 0), Vec3(0, boxSize, 0), Vec3(0, 0, boxSize));
    customSystem.setDefaultPeriodicBoxVectors(Vec3(boxSize, 0, 0), Vec3(0, boxSize, 0), Vec3(0, 0, boxSize));
    standardNonbonded->setUseDispersionCorrection(true);
    customNonbonded->setUseLongRangeCorrection(true);
    standardNonbonded->setUseSwitchingFunction(true);
    customNonbonded->setUseSwitchingFunction(true);
    standardNonbonded->setSwitchingDistance(0.8*cutoff);
    customNonbonded->setSwitchingDistance(0.8*cutoff);
    standardSystem.addForce(standardNonbonded);
    customSystem.addForce(customNonbonded);

    // Compute the correction for the standard force.

    Context context1(standardSystem, integrator1, platform);
    context1.setPositions(positions);
    double standardEnergy1 = context1.getState(State::Energy).getPotentialEnergy();
    standardNonbonded->setUseDispersionCorrection(false);
    context1.reinitialize();
    context1.setPositions(positions);
    double standardEnergy2 = context1.getState(State::Energy).getPotentialEnergy();

    // Compute the correction for the custom force.

    Context context2(customSystem, integrator2, platform);
    context2.setPositions(positions);
    double customEnergy1 = context2.getState(State::Energy).getPotentialEnergy();
    customNonbonded->setUseLongRangeCorrection(false);
    context2.reinitialize();
    context2.setPositions(positions);
    double customEnergy2 = context2.getState(State::Energy).getPotentialEnergy();

    // See if they agree.

    ASSERT_EQUAL_TOL(standardEnergy1-standardEnergy2, customEnergy1-customEnergy2, 1e-4);
}

int main(int argc, char* argv[]) {
    try {
        if (argc > 1)
            platform.setPropertyDefaultValue("OpenCLPrecision", string(argv[1]));
        testSimpleExpression();
        testParameters();
        testManyParameters();
        testExclusions();
        testCutoff();
        testPeriodic();
        testTabulatedFunction();
        testCoulombLennardJones();
        testParallelComputation();
        testSwitchingFunction();
        testLongRangeCorrection();
    }
    catch(const exception& e) {
        cout << "exception: " << e.what() << endl;
        return 1;
    }
    cout << "Done" << endl;
    return 0;
}
