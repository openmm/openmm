/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2016 Stanford University and the Authors.           *
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

#ifdef WIN32
  #define _USE_MATH_DEFINES // Needed to get M_PI
#endif
#include "openmm/internal/AssertionUtilities.h"
#include "openmm/Context.h"
#include "openmm/GayBerneForce.h"
#include "openmm/NonbondedForce.h"
#include "openmm/System.h"
#include "openmm/VerletIntegrator.h"
#include "sfmt/SFMT.h"
#include <iostream>
#include <vector>

using namespace OpenMM;
using namespace std;

const double TOL = 1e-5;

void testPointParticles() {
    // For point particles, it should be identical to a standard Lennard-Jones force.

    const int numParticles = 10;
    const double sigma = 0.5;
    const double epsilon = 1.5;
    System system;
    GayBerneForce* gb = new GayBerneForce();
    NonbondedForce* nb = new NonbondedForce();
    system.addForce(gb);
    system.addForce(nb);
    gb->setForceGroup(1);
    vector<Vec3> positions;
    OpenMM_SFMT::SFMT sfmt;
    init_gen_rand(0, sfmt);
    for (int i = 0; i < numParticles; i++) {
        system.addParticle(1.0);
        gb->addParticle(sigma, epsilon, -1, -1, sigma, sigma, sigma, 1, 1, 1);
        nb->addParticle(0, sigma, epsilon);
        positions.push_back(Vec3(2.0*genrand_real2(sfmt), 2.0*genrand_real2(sfmt), 2.0*genrand_real2(sfmt)));
    }
    VerletIntegrator integ(0.001);

    // Compute forces and energy with each one and compare them.

    Context context(system, integ, platform);
    context.setPositions(positions);
    State state1 = context.getState(State::Forces | State::Energy, false, 1);
    State state2 = context.getState(State::Forces | State::Energy, false, 2);
    ASSERT_EQUAL_TOL(state1.getPotentialEnergy(), state2.getPotentialEnergy(), 1e-5);
    for (int i = 0; i < numParticles; i++)
        ASSERT_EQUAL_VEC(state1.getForces()[i], state2.getForces()[i], 1e-4);
}

void testEnergyScales() {
    // Create two Lennard-Jones particles for which the energy scale factors vary.

    const double sigma = 0.5;
    const double epsilon = 1.5;
    System system;
    for (int i = 0; i < 6; i++)
        system.addParticle(1.0);
    GayBerneForce* gb = new GayBerneForce();
    system.addForce(gb);
    gb->addParticle(sigma, epsilon, 1, 2, sigma, sigma, sigma, 1.1, 1.5, 1.8);
    gb->addParticle(1, 0, -1, -1, 1, 1, 1, 1, 1, 1);
    gb->addParticle(1, 0, -1, -1, 1, 1, 1, 1, 1, 1);
    gb->addParticle(sigma, epsilon, 4, 5, sigma, sigma, sigma, 1.2, 1.6, 1.7);
    gb->addParticle(1, 0, -1, -1, 1, 1, 1, 1, 1, 1);
    gb->addParticle(1, 0, -1, -1, 1, 1, 1, 1, 1, 1);
    vector<Vec3> positions(6);
    positions[0] = Vec3(0, 0, 0);
    positions[1] = Vec3(1, 0, 0);
    positions[2] = Vec3(0, 1, 0);
    positions[3] = Vec3(1, 0, 0);
    positions[4] = Vec3(2, 0, 0);
    positions[5] = Vec3(1, 1, 0);
    VerletIntegrator integ(0.001);
    Context context(system, integ, platform);
    context.setPositions(positions);

    // Depending on their relative orientations, the interaction should be equivalent
    // to LJ with different values of epsilon.

    double expectedEnergy = 4*epsilon*(pow(sigma, 12.0)-pow(sigma, 6.0));
    double expectedForce = 4*epsilon*(12*pow(sigma, 12.0)-6*pow(sigma, 6.0));
    double expectedScale = pow(2.0/(1/sqrt(1.1) + 1/sqrt(1.2)), 2.0);
    State state = context.getState(State::Forces | State::Energy);
    ASSERT_EQUAL_TOL(expectedEnergy*expectedScale, state.getPotentialEnergy(), 1e-5);
    ASSERT_EQUAL_VEC(Vec3(expectedForce*expectedScale, 0, 0), state.getForces()[3], 1e-5);
    positions[3] = Vec3(0, 1, 0);
    positions[4] = Vec3(1, 1, 0);
    positions[5] = Vec3(0, 2, 0);
    context.setPositions(positions);
    expectedScale = pow(2.0/(1/sqrt(1.5) + 1/sqrt(1.6)), 2.0);
    state = context.getState(State::Forces | State::Energy);
    ASSERT_EQUAL_TOL(expectedEnergy*expectedScale, state.getPotentialEnergy(), 1e-5);
    ASSERT_EQUAL_VEC(Vec3(0, expectedForce*expectedScale, 0), state.getForces()[3], 1e-5);
    positions[3] = Vec3(0, 1, 0);
    positions[4] = Vec3(1, 1, 0);
    positions[5] = Vec3(0, 1, 1);
    context.setPositions(positions);
    expectedScale = pow(2.0/(1/sqrt(1.5) + 1/sqrt(1.7)), 2.0);
    state = context.getState(State::Forces | State::Energy);
    ASSERT_EQUAL_TOL(expectedEnergy*expectedScale, state.getPotentialEnergy(), 1e-5);
    ASSERT_EQUAL_VEC(Vec3(0, expectedForce*expectedScale, 0), state.getForces()[3], 1e-5);
    
    // Modify their parameters and see if the result is still correct.
    
    double newSigma = 1.1*sigma;
    gb->setParticleParameters(0, newSigma, 1.5*epsilon, 1, 2, newSigma, newSigma, newSigma, 1.2, 1.6, 1.9);
    gb->setParticleParameters(3, newSigma, epsilon, 4, 5, newSigma, newSigma, newSigma, 1.3, 1.7, 1.8);
    gb->updateParametersInContext(context);
    double combinedEpsilon = sqrt(1.5)*epsilon;
    expectedEnergy = 4*combinedEpsilon*(pow(newSigma, 12.0)-pow(newSigma, 6.0));
    expectedForce = 4*combinedEpsilon*(12*pow(newSigma, 12.0)-6*pow(newSigma, 6.0));
    expectedScale = pow(2.0/(1/sqrt(1.6) + 1/sqrt(1.8)), 2.0);
    state = context.getState(State::Forces | State::Energy);
    ASSERT_EQUAL_TOL(expectedEnergy*expectedScale, state.getPotentialEnergy(), 1e-5);
    ASSERT_EQUAL_VEC(Vec3(0, expectedForce*expectedScale, 0), state.getForces()[3], 1e-5);
}

void testEnergyConservation() {
    // Create a box of ellipsoids and make sure a simulation conserves energy.
    // That verifies that forces and energies are consistent.

    const double boxSize = 3.0;
    System system;
    system.setDefaultPeriodicBoxVectors(Vec3(boxSize, 0, 0), Vec3(0, boxSize, 0), Vec3(0, 0, boxSize));
    GayBerneForce* gb = new GayBerneForce();
    system.addForce(gb);
    gb->setNonbondedMethod(GayBerneForce::CutoffPeriodic);
    gb->setCutoffDistance(1.0);
    gb->setUseSwitchingFunction(true);
    gb->setSwitchingDistance(0.9);
    vector<Vec3> positions;
    for (int x = 0; x < 3; x++) {
        for (int y = 0; y < 3; y++) {
            for (int z = 0; z < 3; z++) {
                int first = system.getNumParticles();
                system.addParticle(10.0);
                system.addParticle(1.0);
                system.addParticle(1.0);
                gb->addParticle(0.2, 10.0, first+1, first+2, 0.2, 0.25, 0.3, 0.9, 1.0, 1.1);
                gb->addParticle(1.0, 0.0, -1, -1, 1, 1, 1, 1, 1, 1);
                gb->addParticle(1.0, 0.0, -1, -1, 1, 1, 1, 1, 1, 1);
                positions.push_back(Vec3(x, y, z));
                positions.push_back(Vec3(x+0.1, y, z));
                positions.push_back(Vec3(x, y+0.1, z));
                system.addConstraint(first, first+1, 0.1);
                system.addConstraint(first, first+2, 0.1);
                system.addConstraint(first+1, first+2, 0.1*sqrt(2.0));
            }
        }
    }
    VerletIntegrator integ(0.0005);
    Context context(system, integ, platform);
    context.setPositions(positions);
    context.setVelocitiesToTemperature(300.0);
    double initialEnergy;
    for (int i = 0; i < 100; i++) {
        integ.step(5);
        State state = context.getState(State::Energy);
        double energy = state.getPotentialEnergy()+state.getKineticEnergy();
        if (i == 0)
            initialEnergy = energy;
        else
            ASSERT_EQUAL_TOL(initialEnergy, energy, 1e-4);
    }
}

void testExceptions() {
    // Create two Lennard-Jones particles for which the energy scale factors vary,
    // then override their interaction with an exception.

    const double sigma = 0.5;
    const double epsilon = 1.5;
    System system;
    for (int i = 0; i < 6; i++)
        system.addParticle(1.0);
    GayBerneForce* gb = new GayBerneForce();
    system.addForce(gb);
    gb->addParticle(sigma, epsilon, 1, 2, sigma, sigma, sigma, 1.1, 1.5, 1.8);
    gb->addParticle(1, 0, -1, -1, 1, 1, 1, 1, 1, 1);
    gb->addParticle(1, 0, -1, -1, 1, 1, 1, 1, 1, 1);
    gb->addParticle(sigma, epsilon, 4, 5, sigma, sigma, sigma, 1.2, 1.6, 1.7);
    gb->addParticle(1, 0, -1, -1, 1, 1, 1, 1, 1, 1);
    gb->addParticle(1, 0, -1, -1, 1, 1, 1, 1, 1, 1);
    gb->addException(0, 3, sigma, 3.5*epsilon);
    vector<Vec3> positions(6);
    positions[0] = Vec3(0, 0, 0);
    positions[1] = Vec3(1, 0, 0);
    positions[2] = Vec3(0, 1, 0);
    positions[3] = Vec3(1, 0, 0);
    positions[4] = Vec3(2, 0, 0);
    positions[5] = Vec3(1, 1, 0);
    VerletIntegrator integ(0.001);
    Context context(system, integ, platform);
    context.setPositions(positions);
    double expectedEnergy = 3.5*4*epsilon*(pow(sigma, 12.0)-pow(sigma, 6.0));
    double expectedForce = 3.5*4*epsilon*(12*pow(sigma, 12.0)-6*pow(sigma, 6.0));
    double expectedScale = pow(2.0/(1/sqrt(1.1) + 1/sqrt(1.2)), 2.0);
    State state = context.getState(State::Forces | State::Energy);
    ASSERT_EQUAL_TOL(expectedEnergy*expectedScale, state.getPotentialEnergy(), 1e-5);
    ASSERT_EQUAL_VEC(Vec3(expectedForce*expectedScale, 0, 0), state.getForces()[3], 1e-5);
    
    // Modify the exception and see if the results are still correct.

    gb->setExceptionParameters(0, 0, 3, sigma, 3.1*epsilon);
    gb->updateParametersInContext(context);
    expectedEnergy = 3.1*4*epsilon*(pow(sigma, 12.0)-pow(sigma, 6.0));
    expectedForce = 3.1*4*epsilon*(12*pow(sigma, 12.0)-6*pow(sigma, 6.0));
    state = context.getState(State::Forces | State::Energy);
    ASSERT_EQUAL_TOL(expectedEnergy*expectedScale, state.getPotentialEnergy(), 1e-5);
    ASSERT_EQUAL_VEC(Vec3(expectedForce*expectedScale, 0, 0), state.getForces()[3], 1e-5);
}

void runPlatformTests();

int main(int argc, char* argv[]) {
    try {
        initializeTests(argc, argv);
        testPointParticles();
        testEnergyScales();
        testEnergyConservation();
        testExceptions();
        runPlatformTests();
    }
    catch(const exception& e) {
        cout << "exception: " << e.what() << endl;
        return 1;
    }
    cout << "Done" << endl;
    return 0;
}
