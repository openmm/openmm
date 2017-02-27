/* -------------------------------------------------------------------------- *
 *                                   OpenMMAmoeba                             *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2008-2012 Stanford University and the Authors.      *
 * Authors: Mark Friedrichs                                                   *
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
 * This tests the CUDA implementation of CudaAmoebaWcaDispersionForce.
 */

#include "openmm/internal/AssertionUtilities.h"
#include "openmm/Context.h"
#include "OpenMMAmoeba.h"
#include "openmm/System.h"
#include "openmm/AmoebaWcaDispersionForce.h"
#include "openmm/LangevinIntegrator.h"
#include <iostream>
#include <vector>
#include <stdlib.h>
#include <stdio.h>

#define ASSERT_EQUAL_TOL_MOD(expected, found, tol, testname) {double _scale_ = std::abs(expected) > 1.0 ? std::abs(expected) : 1.0; if (!(std::abs((expected)-(found))/_scale_ <= (tol))) {std::stringstream details; details << testname << " Expected "<<(expected)<<", found "<<(found); throwException(__FILE__, __LINE__, details.str());}};

#define ASSERT_EQUAL_VEC_MOD(expected, found, tol,testname) {ASSERT_EQUAL_TOL_MOD((expected)[0], (found)[0], (tol),(testname)); ASSERT_EQUAL_TOL_MOD((expected)[1], (found)[1], (tol),(testname)); ASSERT_EQUAL_TOL_MOD((expected)[2], (found)[2], (tol),(testname));};


using namespace OpenMM;
const double TOL = 1e-4;

extern "C" void registerAmoebaCudaKernelFactories();

void compareForcesEnergy(std::string& testName, double expectedEnergy, double energy,
                         const std::vector<Vec3>& expectedForces,
                         const std::vector<Vec3>& forces, double tolerance) {
    for (unsigned int ii = 0; ii < forces.size(); ii++) {
        ASSERT_EQUAL_VEC_MOD(expectedForces[ii], forces[ii], tolerance, testName);
    }
    ASSERT_EQUAL_TOL_MOD(expectedEnergy, energy, tolerance, testName);
}

// test Wca dispersion

void testWcaDispersionAmmonia() {

    std::string testName      = "testWcaDispersionAmmonia";

    int numberOfParticles     = 8;

    // Create the system.
    
    System system;
    AmoebaWcaDispersionForce* amoebaWcaDispersionForce = new AmoebaWcaDispersionForce();;

    amoebaWcaDispersionForce->setEpso(   4.6024000e-01);
    amoebaWcaDispersionForce->setEpsh(   5.6484000e-02);
    amoebaWcaDispersionForce->setRmino(  1.7025000e-01);
    amoebaWcaDispersionForce->setRminh(  1.3275000e-01);
    amoebaWcaDispersionForce->setDispoff(2.6000000e-02);
    amoebaWcaDispersionForce->setAwater( 3.3428000e+01);
    amoebaWcaDispersionForce->setSlevy(  1.0000000e+00);
    amoebaWcaDispersionForce->setShctd(  8.1000000e-01);

    // addParticle: radius, epsilon

    for (unsigned int ii = 0; ii < 2; ii++) {
        system.addParticle(  1.4007000e+01);
        amoebaWcaDispersionForce->addParticle(  1.8550000e-01,   4.3932000e-01);
    
        system.addParticle(  1.0080000e+00);
        amoebaWcaDispersionForce->addParticle(  1.3500000e-01,   8.3680000e-02);
    
        system.addParticle(  1.0080000e+00);
        amoebaWcaDispersionForce->addParticle(  1.3500000e-01,   8.3680000e-02);
    
        system.addParticle(  1.0080000e+00);
        amoebaWcaDispersionForce->addParticle(  1.3500000e-01,   8.3680000e-02);
    }

    std::vector<Vec3> positions(numberOfParticles);

    positions[0]              = Vec3(  1.5927280e-01,   1.7000000e-06,    1.6491000e-03);
    positions[1]              = Vec3(  2.0805540e-01,  -8.1258800e-02,    3.7282500e-02);
    positions[2]              = Vec3(  2.0843610e-01,   8.0953200e-02,    3.7462200e-02);
    positions[3]              = Vec3(  1.7280780e-01,   2.0730000e-04,   -9.8741700e-02);
    positions[4]              = Vec3( -1.6743680e-01,   1.5900000e-05,   -6.6149000e-03);
    positions[5]              = Vec3( -2.0428260e-01,   8.1071500e-02,    4.1343900e-02);
    positions[6]              = Vec3( -6.7308300e-02,   1.2800000e-05,    1.0623300e-02);
    positions[7]              = Vec3( -2.0426290e-01,  -8.1231400e-02,    4.1033500e-02);

    system.addForce(amoebaWcaDispersionForce);

    std::string platformName;
    platformName = "CUDA";
    LangevinIntegrator integrator(0.0, 0.1, 0.01);
    Context context(system, integrator, Platform::getPlatformByName(platformName));

    context.setPositions(positions);
    State state = context.getState(State::Forces | State::Energy);
    std::vector<Vec3> forces = state.getForces();
    double energy = state.getPotentialEnergy();

    // TINKER-computed values

    std::vector<Vec3> expectedForces(numberOfParticles);
    double expectedEnergy     =  -2.6981209e+01;

    expectedForces[0]         = Vec3(  4.7839388e+00,  -7.3510133e-04,  -5.0382764e-01);
    expectedForces[1]         = Vec3(  1.4657758e+00,   1.2431003e+00,  -6.7075886e-01);
    expectedForces[2]         = Vec3(  1.4563936e+00,  -1.2399917e+00,  -6.7443841e-01);
    expectedForces[3]         = Vec3(  2.1116744e+00,  -2.7407512e-03,   1.3271245e+00);
    expectedForces[4]         = Vec3( -4.7528440e+00,  -1.5148066e-03,   1.2653813e+00);
    expectedForces[5]         = Vec3( -1.1875619e+00,  -1.2866678e+00,  -3.9109060e-01);
    expectedForces[6]         = Vec3( -2.6885679e+00,  -4.3038639e-04,   3.3763583e-02);
    expectedForces[7]         = Vec3( -1.1888087e+00,   1.2889802e+00,  -3.8615387e-01);

    double tolerance          = 1.0e-04;
    compareForcesEnergy(testName, expectedEnergy, energy, expectedForces, forces, tolerance);
    
    // Try changing the particle parameters and make sure it's still correct.
    
    for (int i = 0; i < numberOfParticles; i++) {
        double radius, epsilon;
        amoebaWcaDispersionForce->getParticleParameters(i, radius, epsilon);
        amoebaWcaDispersionForce->setParticleParameters(i, 0.9*radius, 2.0*epsilon);
    }
    LangevinIntegrator integrator2(0.0, 0.1, 0.01);
    Context context2(system, integrator2, Platform::getPlatformByName(platformName));
    context2.setPositions(positions);
    State state1 = context.getState(State::Forces | State::Energy);
    State state2 = context2.getState(State::Forces | State::Energy);
    bool exceptionThrown = false;
    try {
        // This should throw an exception.
        compareForcesEnergy(testName, state1.getPotentialEnergy(), state2.getPotentialEnergy(), state1.getForces(), state2.getForces(), tolerance);
    }
    catch (std::exception ex) {
        exceptionThrown = true;
    }
    ASSERT(exceptionThrown);
    amoebaWcaDispersionForce->updateParametersInContext(context);
    state1 = context.getState(State::Forces | State::Energy);
    compareForcesEnergy(testName, state1.getPotentialEnergy(), state2.getPotentialEnergy(), state1.getForces(), state2.getForces(), tolerance);
}

int main(int argc, char* argv[]) {
    try {
        std::cout << "TestCudaAmoebaWcaDispersionForce running test..." << std::endl;
        registerAmoebaCudaKernelFactories();
        if (argc > 1)
            Platform::getPlatformByName("CUDA").setPropertyDefaultValue("Precision", std::string(argv[1]));

        // test Wca dispersion force using two ammonia molecules

        testWcaDispersionAmmonia();


    } catch(const std::exception& e) {
        std::cout << "exception: " << e.what() << std::endl;
        std::cout << "FAIL - ERROR.  Test failed." << std::endl;
        return 1;
    }
    std::cout << "Done" << std::endl;
    return 0;
}
