/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
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
#include "ReferencePlatform.h"
#include "openmm/NonbondedForce.h"
#include "openmm/System.h"
#include "openmm/LangevinIntegrator.h"
#include "openmm/VerletIntegrator.h"
#include "openmm/internal/ContextImpl.h"
#include "openmm/internal/NonbondedForceImpl.h"
#include "SimTKOpenMMRealType.h"
#include "sfmt/SFMT.h"
#include <iostream>
#include <vector>

using namespace OpenMM;
using namespace std;

void testEwaldExact() {

//    Use a NaCl crystal to compare the calculated and Madelung energies

    const int numParticles = 1000;
    const double cutoff = 1.0;
    const double boxSize = 2.82;
    const double ewaldTol = 1e-5;

    System system;
    for (int i = 0; i < numParticles/2; i++)
        system.addParticle(22.99);
    for (int i = 0; i < numParticles/2; i++)
        system.addParticle(35.45);
    VerletIntegrator integrator(0.01);
    NonbondedForce* nonbonded = new NonbondedForce();
    for (int i = 0; i < numParticles/2; i++)
        nonbonded->addParticle(1.0, 1.0,0.0);
    for (int i = 0; i < numParticles/2; i++)
        nonbonded->addParticle(-1.0, 1.0,0.0);
    nonbonded->setNonbondedMethod(NonbondedForce::Ewald);
    nonbonded->setCutoffDistance(cutoff);
    system.setDefaultPeriodicBoxVectors(Vec3(boxSize, 0, 0), Vec3(0, boxSize, 0), Vec3(0, 0, boxSize));
    nonbonded->setEwaldErrorTolerance(ewaldTol);
    system.addForce(nonbonded);
    Context context(system, integrator, platform);
    vector<Vec3> positions(numParticles);
    #include "nacl_crystal.dat"
    context.setPositions(positions);

    State state = context.getState(State::Forces | State::Energy);
    const vector<Vec3>& forces = state.getForces();

//   The potential energy of an ion in a crystal is 
//   E = - (M*e^2/ 4*pi*epsilon0*a0),
//   where 
//   M            :    Madelung constant (dimensionless, for FCC cells such as NaCl it is 1.7476)
//   e            :    1.6022 × 10−19 C
//   4*pi*epsilon0:     1.112 × 10−10 C²/(J m)
//   a0           :    0.282 x 10-9 m (perfect cell)
// 
//   E is then the energy per pair of ions, so for our case
//   E has to be divided by 2 (per ion), multiplied by N(avogadro), multiplied by number of particles, and divided by 1000 for kJ
    double exactEnergy        = - (1.7476 * 1.6022e-19 * 1.6022e-19  * 6.02214e+23 * numParticles) / (1.112e-10 * 0.282e-9 * 2 * 1000);
    //cout << "exact\t\t: " << exactEnergy << endl;
    //cout << "calc\t\t: " << state.getPotentialEnergy() << endl;
    ASSERT_EQUAL_TOL(exactEnergy, state.getPotentialEnergy(), 100*ewaldTol);
}

void testEwaldPME(bool includeExceptions) {

//      Use amorphous NaCl system for the tests

    const int numParticles = 894;
    const double cutoff = 1.2;
    const double boxSize = 3.00646;
    double tol = 1e-5;

    ReferencePlatform reference;
    System system;
    NonbondedForce* nonbonded = new NonbondedForce();
    nonbonded->setNonbondedMethod(NonbondedForce::Ewald);
    nonbonded->setCutoffDistance(cutoff);
    nonbonded->setEwaldErrorTolerance(tol);

    for (int i = 0; i < numParticles/2; i++)
        system.addParticle(22.99);
    for (int i = 0; i < numParticles/2; i++)
        system.addParticle(35.45);
    for (int i = 0; i < numParticles/2; i++)
        nonbonded->addParticle(1.0, 1.0,0.0);
    for (int i = 0; i < numParticles/2; i++)
        nonbonded->addParticle(-1.0, 1.0,0.0);
    system.setDefaultPeriodicBoxVectors(Vec3(boxSize, 0, 0), Vec3(0, boxSize, 0), Vec3(0, 0, boxSize));
    system.addForce(nonbonded);

    vector<Vec3> positions(numParticles);
    #include "nacl_amorph.dat"
    if (includeExceptions) {
        // Add some exclusions.

        for (int i = 0; i < numParticles-1; i++) {
            Vec3 delta = positions[i]-positions[i+1];
            if (sqrt(delta.dot(delta)) < 0.5*cutoff)
                nonbonded->addException(i, i+1, i%2 == 0 ? 0.0 : 0.5, 1.0, 0.0);
        }
    }

//    (1)  Check whether this matches the Reference platform when using Ewald Method

    VerletIntegrator integrator1(0.01);
    VerletIntegrator integrator2(0.01);
    Context context(system, integrator1, platform);
    Context referenceContext(system, integrator2, reference);
    context.setPositions(positions);
    referenceContext.setPositions(positions);
    State state = context.getState(State::Forces | State::Energy);
    State referenceState = referenceContext.getState(State::Forces | State::Energy);
    tol = 1e-2;
    for (int i = 0; i < numParticles; i++) {
        ASSERT_EQUAL_VEC(referenceState.getForces()[i], state.getForces()[i], tol);
    }
    tol = 1e-5;
    ASSERT_EQUAL_TOL(referenceState.getPotentialEnergy(), state.getPotentialEnergy(), tol);
    if (!includeExceptions)
        ASSERT_EQUAL_TOL(-3.82047e+05, state.getPotentialEnergy(), 1e-5); // Value from Gromacs

//    (2) Check whether Ewald method is self-consistent

    double norm = 0.0;
    for (int i = 0; i < numParticles; ++i) {
        Vec3 f = state.getForces()[i];
        norm += f[0]*f[0] + f[1]*f[1] + f[2]*f[2];
    }

    norm = std::sqrt(norm);
    const double delta = 5e-3;
    double step = delta/norm;
    for (int i = 0; i < numParticles; ++i) {
        Vec3 p = positions[i];
        Vec3 f = state.getForces()[i];
        positions[i] = Vec3(p[0]-f[0]*step, p[1]-f[1]*step, p[2]-f[2]*step);
    }
    VerletIntegrator integrator3(0.01);
    Context context2(system, integrator3, platform);
    context2.setPositions(positions);

    tol = 1e-2;
    State state2 = context2.getState(State::Energy);
    ASSERT_EQUAL_TOL(norm, (state2.getPotentialEnergy()-state.getPotentialEnergy())/delta, tol)

//    (3)  Check whether this matches the Reference platform when using PME

    nonbonded->setNonbondedMethod(NonbondedForce::PME);
    context.reinitialize();
    referenceContext.reinitialize();
    context.setPositions(positions);
    referenceContext.setPositions(positions);
    state = context.getState(State::Forces | State::Energy);
    referenceState = referenceContext.getState(State::Forces | State::Energy);
    tol = 1e-2;
    for (int i = 0; i < numParticles; i++) {
        ASSERT_EQUAL_VEC(referenceState.getForces()[i], state.getForces()[i], tol);
    }
    tol = 1e-5;
    ASSERT_EQUAL_TOL(referenceState.getPotentialEnergy(), state.getPotentialEnergy(), tol);
    if (!includeExceptions)
        ASSERT_EQUAL_TOL(-3.82047e+05, state.getPotentialEnergy(), 1e-3); // Value from Gromacs

//    (4) Check whether PME method is self-consistent

    norm = 0.0;
    for (int i = 0; i < numParticles; ++i) {
        Vec3 f = state.getForces()[i];
        norm += f[0]*f[0] + f[1]*f[1] + f[2]*f[2];
    }

    norm = std::sqrt(norm);
    step = delta/norm;
    for (int i = 0; i < numParticles; ++i) {
        Vec3 p = positions[i];
        Vec3 f = state.getForces()[i];
        positions[i] = Vec3(p[0]-f[0]*step, p[1]-f[1]*step, p[2]-f[2]*step);
    }
    VerletIntegrator integrator4(0.01);
    Context context3(system, integrator4, platform);
    context3.setPositions(positions);

    tol = 1e-2;
    State state3 = context3.getState(State::Energy);
    ASSERT_EQUAL_TOL(norm, (state3.getPotentialEnergy()-state.getPotentialEnergy())/delta, tol)
}

void testTriclinic() {
    // Create a triclinic box containing eight particles.

    System system;
    system.setDefaultPeriodicBoxVectors(Vec3(2.5, 0, 0), Vec3(0.5, 3.0, 0), Vec3(0.7, 0.9, 3.5));
    for (int i = 0; i < 8; i++)
        system.addParticle(1.0);
    NonbondedForce* force = new NonbondedForce();
    system.addForce(force);
    force->setNonbondedMethod(NonbondedForce::PME);
    force->setCutoffDistance(1.0);
    force->setPMEParameters(3.45891, 32, 40, 48);
    for (int i = 0; i < 4; i++)
        force->addParticle(-1, 0.440104, 0.4184); // Cl parameters
    for (int i = 0; i < 4; i++)
        force->addParticle(1, 0.332840, 0.0115897); // Na parameters
    vector<Vec3> positions(8);
    positions[0] = Vec3(1.744, 2.788, 3.162);
    positions[1] = Vec3(1.048, 0.762, 2.340);
    positions[2] = Vec3(2.489, 1.570, 2.817);
    positions[3] = Vec3(1.027, 1.893, 3.271);
    positions[4] = Vec3(0.937, 0.825, 0.009);
    positions[5] = Vec3(2.290, 1.887, 3.352);
    positions[6] = Vec3(1.266, 1.111, 2.894);
    positions[7] = Vec3(0.933, 1.862, 3.490);

    // Compute the forces and energy.

    VerletIntegrator integ(0.001);
    Context context(system, integ, platform);
    context.setPositions(positions);
    State state = context.getState(State::Forces | State::Energy);

    // Compare them to values computed by Gromacs.

    double expectedEnergy = -963.370;
    vector<Vec3> expectedForce(8);
    expectedForce[0] = Vec3(4.25253e+01, -1.23503e+02, 1.22139e+02);
    expectedForce[1] = Vec3(9.74752e+01, 1.68213e+02, 1.93169e+02);
    expectedForce[2] = Vec3(-1.50348e+02, 1.29165e+02, 3.70435e+02);
    expectedForce[3] = Vec3(9.18644e+02, -3.52571e+00, -1.34772e+03);
    expectedForce[4] = Vec3(-1.61193e+02, 9.01528e+01, -7.12904e+01);
    expectedForce[5] = Vec3(2.82630e+02, 2.78029e+01, -3.72864e+02);
    expectedForce[6] = Vec3(-1.47454e+02, -2.14448e+02, -3.55789e+02);
    expectedForce[7] = Vec3(-8.82195e+02, -7.39132e+01, 1.46202e+03);
    for (int i = 0; i < 8; i++) {
        ASSERT_EQUAL_VEC(expectedForce[i], state.getForces()[i], 1e-4);
    }
    ASSERT_EQUAL_TOL(expectedEnergy, state.getPotentialEnergy(), 1e-4);
}

void testTriclinic2() {
    // Create a triclinic box containing a large molecule made up of randomly positioned particles and make sure the
    // results match the reference platform.

    if (platform.getName() == "Reference")
        return;
    const int numParticles = 1000;
    System system;
    system.setDefaultPeriodicBoxVectors(Vec3(3.2, 0, 0), Vec3(-1.1, 3.1, 0), Vec3(-1.1, -1.5, 2.7));
    vector<Vec3> positions(numParticles);
    OpenMM_SFMT::SFMT sfmt;
    init_gen_rand(0, sfmt);
    NonbondedForce* force = new NonbondedForce();
    for (int i = 0; i < numParticles; i++) {
        system.addParticle(1.0);
        force->addParticle(i%2 == 0 ? -1.0 : 1.0, 1.0, 0.0);
        positions[i] = Vec3(10*genrand_real2(sfmt)-2, 10*genrand_real2(sfmt)-2, 10*genrand_real2(sfmt)-2);
        if (i > 0) {
            Vec3 delta = positions[i-1]-positions[i];
            system.addConstraint(i-1, i, sqrt(delta.dot(delta)));
        }
    }
    system.addForce(force);
    force->setNonbondedMethod(NonbondedForce::PME);
    force->setCutoffDistance(1.0);
    force->setReciprocalSpaceForceGroup(1);
    force->setPMEParameters(2.62826, 27, 25, 24);

    // Compute the forces and energy.

    VerletIntegrator integ1(0.001);
    Context context1(system, integ1, platform);
    context1.setPositions(positions);
    VerletIntegrator integ2(0.001);
    Context context2(system, integ2, Platform::getPlatformByName("Reference"));
    context2.setPositions(positions);
    State state1 = context1.getState(State::Forces | State::Energy, false, 2);
    State state2 = context2.getState(State::Forces | State::Energy, false, 2);
    ASSERT_EQUAL_TOL(state2.getPotentialEnergy(), state1.getPotentialEnergy(), 1e-4);
    for (int i = 0; i < numParticles; i++)
        ASSERT_EQUAL_VEC(state2.getForces()[i], state1.getForces()[i], 1e-4);
}

void testErrorTolerance(NonbondedForce::NonbondedMethod method) {
    // Create a cloud of random point charges.

    const int numParticles = 51;
    const double boxWidth = 5.0;
    System system;
    system.setDefaultPeriodicBoxVectors(Vec3(boxWidth, 0, 0), Vec3(0, boxWidth, 0), Vec3(0, 0, boxWidth));
    NonbondedForce* force = new NonbondedForce();
    system.addForce(force);
    vector<Vec3> positions(numParticles);
    OpenMM_SFMT::SFMT sfmt;
    init_gen_rand(0, sfmt);

    for (int i = 0; i < numParticles; i++) {
        system.addParticle(1.0);
        force->addParticle(-1.0+i*2.0/(numParticles-1), 1.0, 0.0);
        positions[i] = Vec3(boxWidth*genrand_real2(sfmt), boxWidth*genrand_real2(sfmt), boxWidth*genrand_real2(sfmt));
    }
    force->setNonbondedMethod(method);

    // For various values of the cutoff and error tolerance, see if the actual error is reasonable.

    for (double cutoff = 1.0; cutoff < boxWidth/2; cutoff *= 1.2) {
        force->setCutoffDistance(cutoff);
        vector<Vec3> refForces;
        double norm = 0.0;
        for (double tol = 5e-5; tol < 1e-3; tol *= 2.0) {
            force->setEwaldErrorTolerance(tol);
            VerletIntegrator integrator(0.01);
            Context context(system, integrator, platform);
            context.setPositions(positions);
            State state = context.getState(State::Forces);
            if (refForces.size() == 0) {
                refForces = state.getForces();
                for (int i = 0; i < numParticles; i++)
                    norm += refForces[i].dot(refForces[i]);
                norm = sqrt(norm);
            }
            else {
                double diff = 0.0;
                for (int i = 0; i < numParticles; i++) {
                    Vec3 delta = refForces[i]-state.getForces()[i];
                    diff += delta.dot(delta);
                }
                diff = sqrt(diff)/norm;
                ASSERT(diff < 2*tol);
            }
            if (method == NonbondedForce::PME) {
                // See if the PME parameters were calculated correctly.

                double expectedAlpha, actualAlpha;
                int expectedSize[3], actualSize[3];
                NonbondedForceImpl::calcPMEParameters(system, *force, expectedAlpha, expectedSize[0], expectedSize[1], expectedSize[2], false);
                force->getPMEParametersInContext(context, actualAlpha, actualSize[0], actualSize[1], actualSize[2]);
                ASSERT_EQUAL_TOL(expectedAlpha, actualAlpha, 1e-5);
                for (int i = 0; i < 3; i++) {
                    ASSERT(actualSize[i] >= expectedSize[i]);
                    ASSERT(actualSize[i] < expectedSize[i]+10);
                }
            }
        }
    }
}

void testPMEParameters() {
    // Create a cloud of random point charges.

    const int numParticles = 51;
    const double boxWidth = 4.7;
    System system;
    system.setDefaultPeriodicBoxVectors(Vec3(boxWidth, 0, 0), Vec3(0, boxWidth, 0), Vec3(0, 0, boxWidth));
    NonbondedForce* force = new NonbondedForce();
    system.addForce(force);
    vector<Vec3> positions(numParticles);
    OpenMM_SFMT::SFMT sfmt;
    init_gen_rand(0, sfmt);

    for (int i = 0; i < numParticles; i++) {
        system.addParticle(1.0);
        force->addParticle(-1.0+i*2.0/(numParticles-1), 1.0, 0.0);
        positions[i] = Vec3(boxWidth*genrand_real2(sfmt), boxWidth*genrand_real2(sfmt), boxWidth*genrand_real2(sfmt));
    }
    force->setNonbondedMethod(NonbondedForce::PME);
    
    // Compute the energy with an error tolerance of 1e-3.

    force->setEwaldErrorTolerance(1e-3);
    VerletIntegrator integrator1(0.01);
    Context context1(system, integrator1, platform);
    context1.setPositions(positions);
    double energy1 = context1.getState(State::Energy).getPotentialEnergy();
    
    // Try again with an error tolerance of 1e-4.

    force->setEwaldErrorTolerance(1e-4);
    VerletIntegrator integrator2(0.01);
    Context context2(system, integrator2, platform);
    context2.setPositions(positions);
    double energy2 = context2.getState(State::Energy).getPotentialEnergy();
    
    // Now explicitly set the parameters.  These should match the values that were
    // used for tolerance 1e-3.

    force->setPMEParameters(2.49291157051793, 32, 32, 32);
    VerletIntegrator integrator3(0.01);
    Context context3(system, integrator3, platform);
    context3.setPositions(positions);
    double energy3 = context3.getState(State::Energy).getPotentialEnergy();
    ASSERT_EQUAL_TOL(energy1, energy3, 1e-6);
    ASSERT(fabs((energy1-energy2)/energy1) > 1e-5);
}

void runPlatformTests();

int main(int argc, char* argv[]) {
    try {
        initializeTests(argc, argv);
        testEwaldExact();
        testEwaldPME(false);
        testEwaldPME(true);
        testTriclinic();
        testTriclinic2();
        testErrorTolerance(NonbondedForce::Ewald);
        testErrorTolerance(NonbondedForce::PME);
        testPMEParameters();
        runPlatformTests();
    }
    catch(const exception& e) {
        cout << "exception: " << e.what() << endl;
        return 1;
    }
    cout << "Done" << endl;
    return 0;
}
