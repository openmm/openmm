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
 * This tests the CUDA implementation of CudaAmoebaVdwForce.
 */

#include "openmm/internal/AssertionUtilities.h"
#include "openmm/Context.h"
#include "OpenMMAmoeba.h"
#include "openmm/System.h"
#include "openmm/AmoebaVdwForce.h"
#include "openmm/LangevinIntegrator.h"
#include "sfmt/SFMT.h"
#include <iostream>
#include <vector>
#include <stdlib.h>
#include <stdio.h>

#define ASSERT_EQUAL_TOL_MOD(expected, found, tol, testname) {double _scale_ = std::abs(expected) > 1.0 ? std::abs(expected) : 1.0; if (!(std::abs((expected)-(found))/_scale_ <= (tol))) {std::stringstream details; details << testname << " Expected "<<(expected)<<", found "<<(found); throwException(__FILE__, __LINE__, details.str());}};

#define ASSERT_EQUAL_VEC_MOD(expected, found, tol,testname) {ASSERT_EQUAL_TOL_MOD((expected)[0], (found)[0], (tol),(testname)); ASSERT_EQUAL_TOL_MOD((expected)[1], (found)[1], (tol),(testname)); ASSERT_EQUAL_TOL_MOD((expected)[2], (found)[2], (tol),(testname));};


using namespace OpenMM;
using namespace std;

extern "C" void registerAmoebaCudaKernelFactories();

const double TOL = 1e-4;

void testVdw() {

    System system;
    int numberOfParticles          = 6;
    AmoebaVdwForce* amoebaVdwForce = new AmoebaVdwForce();
    std::string sigmaCombiningRule = std::string("CUBIC-MEAN");
    amoebaVdwForce->setSigmaCombiningRule(sigmaCombiningRule);

    std::string epsilonCombiningRule = std::string("HHG");
    amoebaVdwForce->setEpsilonCombiningRule(epsilonCombiningRule);
    for (int ii = 0; ii < numberOfParticles; ii++) {
        int indexIV;
        double mass, sigma, epsilon, reduction;
        std::vector< int > exclusions;
        if (ii == 0 || ii == 3) {
            mass        = 16.0;
            indexIV     = ii;
            sigma       = 1.70250E+00;
            epsilon     = 1.10000E-01;
            reduction   = 0.0;
        } else {
            mass        = 1.0;
            indexIV     = ii < 3 ? 0 : 3;
            sigma       = 1.32750E+00;
            epsilon     = 1.35000E-02;
            reduction   = 0.91;
        }

        // exclusions

        if (ii < 3) {
            exclusions.push_back (0);
            exclusions.push_back (1);
            exclusions.push_back (2);
        } else {
            exclusions.push_back (3);
            exclusions.push_back (4);
            exclusions.push_back (5);
        }
        system.addParticle(mass);
        amoebaVdwForce->addParticle(indexIV, sigma, epsilon, reduction);
        amoebaVdwForce->setParticleExclusions(ii, exclusions);
    }
    LangevinIntegrator integrator(0.0, 0.1, 0.01);

    std::vector<Vec3> positions(numberOfParticles);
    std::vector<Vec3> expectedForces(numberOfParticles);
    double expectedEnergy;

    positions[0]          = Vec3(-0.254893450E+02, -0.876646600E+01,  0.174761600E+01);
    positions[1]          = Vec3(-0.263489690E+02, -0.907798000E+01,  0.205385100E+01);
    positions[2]          = Vec3(-0.252491680E+02, -0.949411200E+01,  0.115017600E+01);
    positions[3]          = Vec3( 0.172827200E+01,  0.195873090E+02,  0.100059800E+01);
    positions[4]          = Vec3( 0.129370700E+01,  0.190112810E+02,  0.169576300E+01);
    positions[5]          = Vec3( 0.256122300E+01,  0.191601930E+02,  0.854382000E+00);

    double offset         = 27.0;
    for (int ii = 0; ii < 3; ii++) {
       positions[ii][0]      += offset;
       positions[ii][1]      += offset;
    }

    expectedForces[0]     = Vec3( -0.729561040E+03,  0.425828484E+04, -0.769114213E+03);
    expectedForces[1]     = Vec3(  0.181000041E+02,  0.328216639E+02, -0.126210511E+02);
    expectedForces[2]     = Vec3( -0.943743014E+00,  0.199728310E+02,  0.884567842E+00);
    expectedForces[3]     = Vec3(  0.615734500E+01, -0.747350431E+03,  0.264726489E+03);
    expectedForces[4]     = Vec3(  0.735772031E+03, -0.353310112E+04,  0.490066356E+03);
    expectedForces[5]     = Vec3( -0.295245970E+02, -0.306277797E+02,  0.260578506E+02);

    expectedEnergy        = 0.740688488E+03;

    system.addForce(amoebaVdwForce);
    std::string platformName;
    #define AngstromToNm 0.1    
    #define CalToJoule   4.184    
    for (int ii = 0; ii < numberOfParticles; ii++) {
        positions[ii][0] *= AngstromToNm;
        positions[ii][1] *= AngstromToNm;
        positions[ii][2] *= AngstromToNm;
    }
    for (int ii = 0; ii < amoebaVdwForce->getNumParticles();  ii++) {
        int indexIV;
        double sigma, epsilon, reduction;
        amoebaVdwForce->getParticleParameters(ii, indexIV, sigma, epsilon, reduction);
        sigma        *= AngstromToNm;
        epsilon      *= CalToJoule;
        amoebaVdwForce->setParticleParameters(ii, indexIV, sigma, epsilon, reduction);
    }
    platformName = "CUDA";
    Context context(system, integrator, Platform::getPlatformByName(platformName));

    context.setPositions(positions);
    State state                      = context.getState(State::Forces | State::Energy);
    std::vector<Vec3> forces         = state.getForces();
    const double conversion          = -AngstromToNm/CalToJoule;

    for (unsigned int ii = 0; ii < forces.size(); ii++) {
        forces[ii][0] *= conversion;
        forces[ii][1] *= conversion;
        forces[ii][2] *= conversion;
    }    
    expectedEnergy *= CalToJoule;
    double tolerance = 1.0e-03;
    for (unsigned int ii = 0; ii < forces.size(); ii++) {
        ASSERT_EQUAL_VEC(expectedForces[ii], forces[ii], tolerance);
    }
    ASSERT_EQUAL_TOL(expectedEnergy, state.getPotentialEnergy(), tolerance);
    
    // Try changing the particle parameters and make sure it's still correct.
    
    for (int i = 0; i < numberOfParticles; i++) {
        int indexIV;
        double mass, sigma, epsilon, reduction;
        amoebaVdwForce->getParticleParameters(i, indexIV, sigma, epsilon, reduction);
        amoebaVdwForce->setParticleParameters(i, indexIV, 0.9*sigma, 2.0*epsilon, 0.95*reduction);
    }
    LangevinIntegrator integrator2(0.0, 0.1, 0.01);
    Context context2(system, integrator2, Platform::getPlatformByName(platformName));
    context2.setPositions(positions);
    State state1 = context.getState(State::Forces | State::Energy);
    State state2 = context2.getState(State::Forces | State::Energy);
    bool exceptionThrown = false;
    try {
        // This should throw an exception.
        for (int i = 0; i < numberOfParticles; i++)
            ASSERT_EQUAL_VEC(state1.getForces()[i], state2.getForces()[i], tolerance);
    }
    catch (std::exception ex) {
        exceptionThrown = true;
    }
    ASSERT(exceptionThrown);
    amoebaVdwForce->updateParametersInContext(context);
    state1 = context.getState(State::Forces | State::Energy);
    for (int i = 0; i < numberOfParticles; i++)
        ASSERT_EQUAL_VEC(state1.getForces()[i], state2.getForces()[i], tolerance);
    ASSERT_EQUAL_TOL(state1.getPotentialEnergy(), state2.getPotentialEnergy(), tolerance);
}

void setupAndGetForcesEnergyVdwAmmonia(const std::string& sigmaCombiningRule, const std::string& epsilonCombiningRule, double cutoff,
                                       double boxDimension, std::vector<Vec3>& forces, double& energy) {

    // beginning of Vdw setup

    System system;
    AmoebaVdwForce* amoebaVdwForce        = new AmoebaVdwForce();;
    int numberOfParticles                 = 8;
    amoebaVdwForce->setSigmaCombiningRule(sigmaCombiningRule);
    amoebaVdwForce->setEpsilonCombiningRule(epsilonCombiningRule);
    amoebaVdwForce->setCutoff(cutoff);
    if (boxDimension > 0.0) {
        Vec3 a(boxDimension, 0.0, 0.0);
        Vec3 b(0.0, boxDimension, 0.0);
        Vec3 c(0.0, 0.0, boxDimension);
        system.setDefaultPeriodicBoxVectors(a, b, c);
        amoebaVdwForce->setNonbondedMethod(AmoebaVdwForce::CutoffPeriodic);
        amoebaVdwForce->setUseDispersionCorrection(1);
    } else {
        amoebaVdwForce->setNonbondedMethod(AmoebaVdwForce::NoCutoff);
        amoebaVdwForce->setUseDispersionCorrection(0);
    }

    // addParticle: ivIndex, radius, epsilon, reductionFactor

    system.addParticle(  1.4007000e+01);
    amoebaVdwForce->addParticle(0,   1.8550000e-01,   4.3932000e-01,   0.0000000e+00);

    system.addParticle(  1.0080000e+00);
    amoebaVdwForce->addParticle(0,   1.3500000e-01,   8.3680000e-02,   9.1000000e-01);

    system.addParticle(  1.0080000e+00);
    amoebaVdwForce->addParticle(0,   1.3500000e-01,   8.3680000e-02,   9.1000000e-01);

    system.addParticle(  1.0080000e+00);
    amoebaVdwForce->addParticle(0,   1.3500000e-01,   8.3680000e-02,   9.1000000e-01);

    system.addParticle(  1.4007000e+01);
    amoebaVdwForce->addParticle(4,   1.8550000e-01,   4.3932000e-01,   0.0000000e+00);

    system.addParticle(  1.0080000e+00);
    amoebaVdwForce->addParticle(4,   1.3500000e-01,   8.3680000e-02,   9.1000000e-01);

    system.addParticle(  1.0080000e+00);
    amoebaVdwForce->addParticle(4,   1.3500000e-01,   8.3680000e-02,   9.1000000e-01);

    system.addParticle(  1.0080000e+00);
    amoebaVdwForce->addParticle(4,   1.3500000e-01,   8.3680000e-02,   9.1000000e-01);

    // ParticleExclusions

    std::vector< int > exclusions;
    exclusions.resize(0);
    exclusions.push_back(0);
    exclusions.push_back(1);
    exclusions.push_back(2);
    exclusions.push_back(3);
    amoebaVdwForce->setParticleExclusions(0, exclusions);

    exclusions.resize(0);
    exclusions.push_back(1);
    exclusions.push_back(0);
    exclusions.push_back(2);
    exclusions.push_back(3);
    amoebaVdwForce->setParticleExclusions(1, exclusions);

    exclusions.resize(0);
    exclusions.push_back(2);
    exclusions.push_back(0);
    exclusions.push_back(1);
    exclusions.push_back(3);
    amoebaVdwForce->setParticleExclusions(2, exclusions);

    exclusions.resize(0);
    exclusions.push_back(3);
    exclusions.push_back(0);
    exclusions.push_back(1);
    exclusions.push_back(2);
    amoebaVdwForce->setParticleExclusions(3, exclusions);

    exclusions.resize(0);
    exclusions.push_back(4);
    exclusions.push_back(5);
    exclusions.push_back(6);
    exclusions.push_back(7);
    amoebaVdwForce->setParticleExclusions(4, exclusions);

    exclusions.resize(0);
    exclusions.push_back(5);
    exclusions.push_back(4);
    exclusions.push_back(6);
    exclusions.push_back(7);
    amoebaVdwForce->setParticleExclusions(5, exclusions);

    exclusions.resize(0);
    exclusions.push_back(6);
    exclusions.push_back(4);
    exclusions.push_back(5);
    exclusions.push_back(7);
    amoebaVdwForce->setParticleExclusions(6, exclusions);

    exclusions.resize(0);
    exclusions.push_back(7);
    exclusions.push_back(4);
    exclusions.push_back(5);
    exclusions.push_back(6);
    amoebaVdwForce->setParticleExclusions(7, exclusions);

    // end of Vdw setup

    std::vector<Vec3> positions(numberOfParticles);

    positions[0]              = Vec3(  1.5927280e-01,   1.7000000e-06,    1.6491000e-03);
    positions[1]              = Vec3(  2.0805540e-01,  -8.1258800e-02,    3.7282500e-02);
    positions[2]              = Vec3(  2.0843610e-01,   8.0953200e-02,    3.7462200e-02);
    positions[3]              = Vec3(  1.7280780e-01,   2.0730000e-04,   -9.8741700e-02);
    positions[4]              = Vec3( -1.6743680e-01,   1.5900000e-05,   -6.6149000e-03);
    positions[5]              = Vec3( -2.0428260e-01,   8.1071500e-02,    4.1343900e-02);
    positions[6]              = Vec3( -6.7308300e-02,   1.2800000e-05,    1.0623300e-02);
    positions[7]              = Vec3( -2.0426290e-01,  -8.1231400e-02,    4.1033500e-02);

    system.addForce(amoebaVdwForce);

    std::string platformName;
    platformName = "CUDA";
    LangevinIntegrator integrator(0.0, 0.1, 0.01);
    Context context(system, integrator, Platform::getPlatformByName(platformName));

    context.setPositions(positions);
    State state                      = context.getState(State::Forces | State::Energy);
    forces                           = state.getForces();
    energy                           = state.getPotentialEnergy();
}

void compareForcesEnergy(std::string& testName, double expectedEnergy, double energy,
                         std::vector<Vec3>& expectedForces,
                         std::vector<Vec3>& forces, double tolerance) {
    for (unsigned int ii = 0; ii < forces.size(); ii++) {
        ASSERT_EQUAL_VEC_MOD(expectedForces[ii], forces[ii], tolerance, testName);
    }
    ASSERT_EQUAL_TOL_MOD(expectedEnergy, energy, tolerance, testName);
}

// test VDW w/ sigmaRule=CubicMean and epsilonRule=HHG

void testVdwAmmoniaCubicMeanHhg() {

    std::string testName      = "testVdwAmmoniaCubicMeanHhg";

    int numberOfParticles     = 8;
    double boxDimension       = -1.0;
    double cutoff             = 9000000.0;
    std::vector<Vec3> forces;
    double energy;

    setupAndGetForcesEnergyVdwAmmonia("CUBIC-MEAN", "HHG", cutoff, boxDimension, forces, energy);
    std::vector<Vec3> expectedForces(numberOfParticles);

    double expectedEnergy     =  4.8012258e+00;

    expectedForces[0]         = Vec3(  2.9265247e+02,  -1.4507808e-02,  -6.9562123e+00);
    expectedForces[1]         = Vec3( -2.2451693e+00,   4.8143073e-01,  -2.0041494e-01);
    expectedForces[2]         = Vec3( -2.2440698e+00,  -4.7905450e-01,  -2.0125284e-01);
    expectedForces[3]         = Vec3( -1.0840394e+00,  -5.8531253e-04,   2.6934135e-01);
    expectedForces[4]         = Vec3( -5.6305662e+01,   1.4733908e-03,  -1.8083306e-01);
    expectedForces[5]         = Vec3(  1.6750145e+00,  -3.2448374e-01,  -1.8030914e-01);
    expectedForces[6]         = Vec3( -2.3412420e+02,   1.0754069e-02,   7.6287492e+00);
    expectedForces[7]         = Vec3(  1.6756544e+00,   3.2497316e-01,  -1.7906832e-01);

    double tolerance          = 1.0e-04;
    compareForcesEnergy(testName, expectedEnergy, energy, expectedForces, forces, tolerance);
}

// test VDW w/ sigmaRule=Arithmetic and epsilonRule=Arithmetic

void testVdwAmmoniaArithmeticArithmetic() {

    std::string testName      = "testVdwAmmoniaArithmeticArithmetic";

    int numberOfParticles     = 8;
    double boxDimension       = -1.0;
    double cutoff             = 9000000.0;

    std::vector<Vec3> forces;
    double energy;
    setupAndGetForcesEnergyVdwAmmonia("ARITHMETIC", "ARITHMETIC", cutoff, boxDimension, forces, energy);
    std::vector<Vec3> expectedForces(numberOfParticles);

    double expectedEnergy     =  4.2252403e+00;

    expectedForces[0]         = Vec3(  3.0603839e+02,  -1.5550310e-02,  -7.2661707e+00);
    expectedForces[1]         = Vec3( -2.7801357e+00,   5.8805051e-01,  -2.5907269e-01);
    expectedForces[2]         = Vec3( -2.7753968e+00,  -5.8440732e-01,  -2.5969111e-01);
    expectedForces[3]         = Vec3( -2.2496416e+00,  -1.1797440e-03,   5.5501757e-01);
    expectedForces[4]         = Vec3( -5.5077629e+01,   8.3417114e-04,  -3.3668921e-01);
    expectedForces[5]         = Vec3(  2.3752452e+00,  -4.6788669e-01,  -2.4907764e-01);
    expectedForces[6]         = Vec3( -2.4790697e+02,   1.1419770e-02,   8.0629999e+00);
    expectedForces[7]         = Vec3(  2.3761408e+00,   4.6871961e-01,  -2.4731607e-01);

    double tolerance          = 1.0e-04;
    compareForcesEnergy(testName, expectedEnergy, energy, expectedForces, forces, tolerance);
}

// test VDW w/ sigmaRule=Geometric and epsilonRule=Geometric

void testVdwAmmoniaGeometricGeometric() {

    std::string testName      = "testVdwAmmoniaGeometricGeometric";

    int numberOfParticles     = 8;
    double boxDimension       = -1.0;
    double cutoff             = 9000000.0;
    std::vector<Vec3> forces;
    double energy;
    setupAndGetForcesEnergyVdwAmmonia("GEOMETRIC", "GEOMETRIC", cutoff, boxDimension, forces, energy);
    std::vector<Vec3> expectedForces(numberOfParticles);

    double expectedEnergy     =  2.5249914e+00;

    expectedForces[0]         = Vec3(  2.1169631e+02,  -1.0710925e-02,  -4.3728025e+00);
    expectedForces[1]         = Vec3( -2.2585621e+00,   4.8409995e-01,  -2.0188344e-01);
    expectedForces[2]         = Vec3( -2.2551351e+00,  -4.8124855e-01,  -2.0246986e-01);
    expectedForces[3]         = Vec3( -1.7178028e+00,  -9.0851787e-04,   4.2466975e-01);
    expectedForces[4]         = Vec3( -4.8302147e+01,   9.6603376e-04,  -5.7972068e-01);
    expectedForces[5]         = Vec3(  1.8100634e+00,  -3.5214093e-01,  -1.9357207e-01);
    expectedForces[6]         = Vec3( -1.6078365e+02,   7.2117601e-03,   5.3180261e+00);
    expectedForces[7]         = Vec3(  1.8109211e+00,   3.5273117e-01,  -1.9224723e-01);

    double tolerance          = 1.0e-04;
    compareForcesEnergy(testName, expectedEnergy, energy, expectedForces, forces, tolerance);
}

void testVdwAmmoniaCubicMeanHarmonic() {

    std::string testName      = "testVdwAmmoniaCubicMeanHarmonic";

    int numberOfParticles     = 8;
    double boxDimension       = -1.0;
    double cutoff             = 9000000.0;
    std::vector<Vec3> forces;
    double energy;
    setupAndGetForcesEnergyVdwAmmonia("CUBIC-MEAN", "HARMONIC", cutoff, boxDimension, forces, energy);
    std::vector<Vec3> expectedForces(numberOfParticles);

    double expectedEnergy     =  4.1369069e+00;

    expectedForces[0]         = Vec3(  2.5854436e+02,  -1.2779529e-02,  -5.9041148e+00);
    expectedForces[1]         = Vec3( -2.0832419e+00,   4.4915831e-01,  -1.8266000e-01);
    expectedForces[2]         = Vec3( -2.0823991e+00,  -4.4699804e-01,  -1.8347141e-01);
    expectedForces[3]         = Vec3( -9.5914714e-01,  -5.2162026e-04,   2.3873165e-01);
    expectedForces[4]         = Vec3( -5.3724787e+01,   1.4838241e-03,  -2.8089191e-01);
    expectedForces[5]         = Vec3(  1.5074325e+00,  -2.9016397e-01,  -1.6385118e-01);
    expectedForces[6]         = Vec3( -2.0271029e+02,   9.2367947e-03,   6.6389988e+00);
    expectedForces[7]         = Vec3(  1.5080748e+00,   2.9058422e-01,  -1.6274118e-01);

    double tolerance          = 1.0e-04;
    compareForcesEnergy(testName, expectedEnergy, energy, expectedForces, forces, tolerance);
}

// test w/ cutoff=0.25 nm; single ixn between two particles (0 and 6); force nonzero on
// particle 4 due to reduction applied to NH
// the distance between 0 and 6 is ~ 0.235 so the ixn is in the tapered region

void testVdwTaper() {

    std::string testName      = "testVdwTaper";

    int numberOfParticles     = 8;
    double boxDimension       = 50.0;
    double cutoff             = 0.25;

    std::vector<Vec3> forces;
    double energy;
    setupAndGetForcesEnergyVdwAmmonia("CUBIC-MEAN", "HHG", cutoff, boxDimension, forces, energy);
    std::vector<Vec3> expectedForces(numberOfParticles);

    double expectedEnergy     =  3.5478444e+00;

    expectedForces[0]         = Vec3(  5.6710779e+02,  -2.7391004e-02,  -1.7867730e+01);
    expectedForces[1]         = Vec3( -0.0000000e+00,  -0.0000000e+00,  -0.0000000e+00);
    expectedForces[2]         = Vec3( -0.0000000e+00,  -0.0000000e+00,  -0.0000000e+00);
    expectedForces[3]         = Vec3( -0.0000000e+00,  -0.0000000e+00,  -0.0000000e+00);
    expectedForces[4]         = Vec3( -5.1039701e+01,   2.4651903e-03,   1.6080957e+00);
    expectedForces[5]         = Vec3( -0.0000000e+00,  -0.0000000e+00,  -0.0000000e+00);
    expectedForces[6]         = Vec3( -5.1606809e+02,   2.4925813e-02,   1.6259634e+01);
    expectedForces[7]         = Vec3( -0.0000000e+00,  -0.0000000e+00,  -0.0000000e+00);

    double tolerance          = 1.0e-04;
    compareForcesEnergy(testName, expectedEnergy, energy, expectedForces, forces, tolerance);
}

// test PBC

void testVdwPBC() {

    std::string testName      = "testVdwPBC";

    int numberOfParticles     = 8;
    double boxDimension       = 0.6;
    double cutoff             = 0.25;

    std::vector<Vec3> forces;
    double energy;
    setupAndGetForcesEnergyVdwAmmonia("CUBIC-MEAN", "HHG", cutoff, boxDimension, forces, energy);
    std::vector<Vec3> expectedForces(numberOfParticles);

    double expectedEnergy     =  1.4949141e+01;

    expectedForces[0]         = Vec3(  5.1453069e+02,   4.9751912e-01,  -1.2759570e+01);
    expectedForces[1]         = Vec3( -2.5622586e+02,  -4.6524265e+01,   2.4281465e+01);
    expectedForces[2]         = Vec3( -2.7538705e+02,   5.1831690e+01,   2.7367710e+01);
    expectedForces[3]         = Vec3( -0.0000000e+00,  -0.0000000e+00,  -0.0000000e+00);
    expectedForces[4]         = Vec3(  3.0883034e+02,  -5.8876974e+00,  -5.8286122e+01);
    expectedForces[5]         = Vec3(  1.1319359e+02,  -3.2047069e-01,   1.6181231e+00);
    expectedForces[6]         = Vec3( -5.1606809e+02,   2.4925813e-02,   1.6259634e+01);
    expectedForces[7]         = Vec3(  1.1112638e+02,   3.7829857e-01,   1.5187587e+00);

    // tolerance is higher here due to interpolation used in setting tapering coefficients;
    // if tapering turned off, then absolute difference < 2.0e-05

    double tolerance          = 5.0e-04;
    compareForcesEnergy(testName, expectedEnergy, energy, expectedForces, forces, tolerance);
}

// create box of 216 water molecules

void setupAndGetForcesEnergyVdwWater(const std::string& sigmaCombiningRule, const std::string& epsilonCombiningRule, double cutoff,
                                     double boxDimension, int includeVdwDispersionCorrection,
                                     std::vector<Vec3>& forces, double& energy) {

    // beginning of Vdw setup

    System system;
    AmoebaVdwForce* amoebaVdwForce        = new AmoebaVdwForce();;
    int numberOfParticles                 = 648;
    amoebaVdwForce->setSigmaCombiningRule(sigmaCombiningRule);
    amoebaVdwForce->setEpsilonCombiningRule(epsilonCombiningRule);
    amoebaVdwForce->setCutoff(cutoff);
    if (boxDimension > 0.0) {
        Vec3 a(boxDimension, 0.0, 0.0);
        Vec3 b(0.0, boxDimension, 0.0);
        Vec3 c(0.0, 0.0, boxDimension);
        system.setDefaultPeriodicBoxVectors(a, b, c);
        amoebaVdwForce->setNonbondedMethod(AmoebaVdwForce::CutoffPeriodic);
        amoebaVdwForce->setUseDispersionCorrection(includeVdwDispersionCorrection);
    } else {
        amoebaVdwForce->setNonbondedMethod(AmoebaVdwForce::NoCutoff);
        amoebaVdwForce->setUseDispersionCorrection(0);
    }

    // addParticle: ivIndex, radius, epsilon, reductionFactor

    int classIndex = 0;
    for (unsigned int ii = 0; ii < numberOfParticles; ii += 3) {

       system.addParticle(  1.5995000e+01);
       amoebaVdwForce->addParticle(ii, 1.7025000e-01,   4.6024000e-01,   0.0000000e+00);

       system.addParticle(  1.0080000e+00);
       amoebaVdwForce->addParticle(ii, 1.3275000e-01,   5.6484000e-02,   9.1000000e-01);

       system.addParticle(  1.0080000e+00);
       amoebaVdwForce->addParticle(ii, 1.3275000e-01,   5.6484000e-02,   9.1000000e-01);
   }

    // exclusions

    std::vector< int > exclusions(3);
    for (unsigned int ii = 0; ii < numberOfParticles; ii += 3) {
        exclusions[0] = ii;
        exclusions[1] = ii+1;
        exclusions[2] = ii+2;
        amoebaVdwForce->setParticleExclusions(ii, exclusions);
        amoebaVdwForce->setParticleExclusions(ii+1, exclusions);
        amoebaVdwForce->setParticleExclusions(ii+2, exclusions);
    }

    // end of Vdw setup

    // set positions

    std::vector<Vec3> positions(numberOfParticles);

    positions[0]             = Vec3(  8.0394300e-01,   5.8680350e-01,    4.9277700e-02);
    positions[1]             = Vec3(  7.5814940e-01,   5.0226660e-01,    4.0375900e-02);
    positions[2]             = Vec3(  8.2870560e-01,   6.0624400e-01,   -3.9707400e-02);
    positions[3]             = Vec3(  1.1484000e-02,  -8.8765990e-01,    6.4458520e-01);
    positions[4]             = Vec3(  9.5892500e-02,  -8.4464940e-01,    6.4052470e-01);
    positions[5]             = Vec3(  2.4723500e-02,  -9.7944710e-01,    6.1378930e-01);
    positions[6]             = Vec3( -6.5763670e-01,  -2.5260000e-02,    8.1046320e-01);
    positions[7]             = Vec3( -6.6454990e-01,   6.8992500e-02,    7.8963560e-01);
    positions[8]             = Vec3( -6.6845370e-01,  -4.0076000e-02,    9.1037470e-01);
    positions[9]             = Vec3(  6.5831270e-01,   8.5501500e-02,   -6.6685290e-01);
    positions[10]            = Vec3(  6.2600580e-01,   8.8732600e-02,   -5.7651320e-01);
    positions[11]            = Vec3(  6.1694860e-01,   5.3229000e-03,   -7.0543230e-01);
    positions[12]            = Vec3(  5.4954790e-01,   6.4357640e-01,    1.8420070e-01);
    positions[13]            = Vec3(  4.7740750e-01,   6.5609280e-01,    1.2079650e-01);
    positions[14]            = Vec3(  6.2544340e-01,   6.3485600e-01,    1.2346110e-01);
    positions[15]            = Vec3( -4.6646340e-01,  -8.5021310e-01,   -2.6526210e-01);
    positions[16]            = Vec3( -4.5053590e-01,  -8.3883300e-01,   -3.6069710e-01);
    positions[17]            = Vec3( -5.5653260e-01,  -8.7510810e-01,   -2.5955820e-01);
    positions[18]            = Vec3( -7.7550740e-01,  -4.6613180e-01,    4.9045930e-01);
    positions[19]            = Vec3( -7.3577510e-01,  -5.4400590e-01,    5.3107060e-01);
    positions[20]            = Vec3( -7.0755520e-01,  -4.1773140e-01,    4.4037930e-01);
    positions[21]            = Vec3( -2.8190600e-02,   7.4872450e-01,   -7.6855300e-01);
    positions[22]            = Vec3( -7.9443300e-02,   7.4463600e-01,   -6.8256160e-01);
    positions[23]            = Vec3(  1.7033100e-02,   8.3813000e-01,   -7.6365310e-01);
    positions[24]            = Vec3( -3.7112750e-01,  -2.2624390e-01,   -1.9170030e-01);
    positions[25]            = Vec3( -4.4236150e-01,  -2.4258640e-01,   -2.4723220e-01);
    positions[26]            = Vec3( -4.0233380e-01,  -2.2106530e-01,   -9.7227800e-02);
    positions[27]            = Vec3( -5.8120030e-01,  -5.6157220e-01,    8.3549400e-02);
    positions[28]            = Vec3( -6.6764500e-01,  -5.7119710e-01,    1.2970660e-01);
    positions[29]            = Vec3( -5.1434340e-01,  -5.5317060e-01,    1.5597670e-01);
    positions[30]            = Vec3(  8.5281410e-01,   4.9997870e-01,    3.4439320e-01);
    positions[31]            = Vec3(  8.8661040e-01,   4.7595500e-01,    4.3409810e-01);
    positions[32]            = Vec3(  7.6829200e-01,   4.5403270e-01,    3.3783460e-01);
    positions[33]            = Vec3(  6.2913000e-03,   3.9622090e-01,   -6.4448110e-01);
    positions[34]            = Vec3( -6.4546800e-02,   4.4539620e-01,   -6.0008300e-01);
    positions[35]            = Vec3(  7.0262000e-03,   3.1229330e-01,   -5.9892730e-01);
    positions[36]            = Vec3(  1.6883500e-02,   6.5824910e-01,    6.0982750e-01);
    positions[37]            = Vec3(  2.9114400e-02,   6.3714540e-01,    7.0403040e-01);
    positions[38]            = Vec3( -3.9569500e-02,   5.9419720e-01,    5.6714930e-01);
    positions[39]            = Vec3(  3.7393550e-01,   6.2909200e-01,    8.1318410e-01);
    positions[40]            = Vec3(  4.1500630e-01,   6.1010560e-01,    9.0110400e-01);
    positions[41]            = Vec3(  4.3953600e-01,   5.9208230e-01,    7.5268270e-01);
    positions[42]            = Vec3(  3.2500410e-01,   4.5615770e-01,   -2.5643980e-01);
    positions[43]            = Vec3(  3.7432790e-01,   4.5313140e-01,   -3.3754880e-01);
    positions[44]            = Vec3(  2.6987370e-01,   5.3785040e-01,   -2.4860760e-01);
    positions[45]            = Vec3(  5.6184630e-01,   5.2015900e-01,    6.3763990e-01);
    positions[46]            = Vec3(  5.6189080e-01,   5.6140190e-01,    5.5312940e-01);
    positions[47]            = Vec3(  5.4901540e-01,   4.2688810e-01,    6.2109450e-01);
    positions[48]            = Vec3( -8.7750980e-01,   6.9408570e-01,   -6.1784650e-01);
    positions[49]            = Vec3( -8.2179580e-01,   7.3187880e-01,   -5.4705510e-01);
    positions[50]            = Vec3( -9.0362240e-01,   7.7367480e-01,   -6.6488210e-01);
    positions[51]            = Vec3( -6.9406820e-01,   2.2491740e-01,    7.1940890e-01);
    positions[52]            = Vec3( -7.3674620e-01,   2.2091000e-01,    6.3486690e-01);
    positions[53]            = Vec3( -7.4149900e-01,   2.8970280e-01,    7.7200060e-01);
    positions[54]            = Vec3(  4.8285280e-01,  -1.8445100e-02,    3.1521130e-01);
    positions[55]            = Vec3(  5.5574910e-01,   2.4338500e-02,    2.7236750e-01);
    positions[56]            = Vec3(  4.1347360e-01,   5.0063500e-02,    3.2371450e-01);
    positions[57]            = Vec3( -2.2024800e-01,  -3.1071870e-01,    9.1706370e-01);
    positions[58]            = Vec3( -2.3195790e-01,  -4.0722320e-01,    9.2465160e-01);
    positions[59]            = Vec3( -2.8015290e-01,  -2.9349640e-01,    8.4209880e-01);
    positions[60]            = Vec3(  1.6893780e-01,   6.6734280e-01,   -2.4352040e-01);
    positions[61]            = Vec3(  1.9716270e-01,   7.5186390e-01,   -2.0536790e-01);
    positions[62]            = Vec3(  8.7430700e-02,   6.4225300e-01,   -1.9539020e-01);
    positions[63]            = Vec3( -9.0804840e-01,  -6.2437310e-01,   -8.8188300e-02);
    positions[64]            = Vec3( -8.6732940e-01,  -7.0428590e-01,   -4.8030200e-02);
    positions[65]            = Vec3( -8.3644480e-01,  -5.8139450e-01,   -1.3828190e-01);
    positions[66]            = Vec3( -8.6567760e-01,  -8.6537570e-01,    5.6295900e-02);
    positions[67]            = Vec3( -8.1778220e-01,  -9.4654890e-01,    8.4163600e-02);
    positions[68]            = Vec3( -9.4534460e-01,  -8.6858770e-01,    1.0560810e-01);
    positions[69]            = Vec3( -5.7716930e-01,  -2.6316670e-01,   -4.5880740e-01);
    positions[70]            = Vec3( -5.4569620e-01,  -3.1693230e-01,   -5.2720970e-01);
    positions[71]            = Vec3( -5.5496000e-01,  -1.7071220e-01,   -4.7392400e-01);
    positions[72]            = Vec3(  7.2367810e-01,  -8.4678300e-01,   -6.9502250e-01);
    positions[73]            = Vec3(  7.9899670e-01,  -8.9648580e-01,   -7.2759260e-01);
    positions[74]            = Vec3(  7.5075030e-01,  -8.1725850e-01,   -6.0600380e-01);
    positions[75]            = Vec3( -2.3769060e-01,  -6.2523350e-01,    1.2921080e-01);
    positions[76]            = Vec3( -1.8309420e-01,  -6.2163180e-01,    4.8693900e-02);
    positions[77]            = Vec3( -2.3929030e-01,  -5.3708810e-01,    1.6453540e-01);
    positions[78]            = Vec3(  8.3347800e-02,  -5.0189060e-01,    5.4317800e-01);
    positions[79]            = Vec3(  1.0917180e-01,  -5.7641330e-01,    4.8632230e-01);
    positions[80]            = Vec3(  1.4837200e-02,  -5.5084220e-01,    5.9546910e-01);
    positions[81]            = Vec3(  7.4250070e-01,  -2.7418580e-01,    8.3795900e-02);
    positions[82]            = Vec3(  6.8666720e-01,  -2.4554090e-01,    1.6206940e-01);
    positions[83]            = Vec3(  7.1516850e-01,  -3.6419530e-01,    7.2493400e-02);
    positions[84]            = Vec3( -2.5059100e-02,   8.6314620e-01,    2.2861410e-01);
    positions[85]            = Vec3(  9.6445000e-03,   9.0720400e-01,    1.4964290e-01);
    positions[86]            = Vec3(  4.5097900e-02,   8.7155360e-01,    2.9051950e-01);
    positions[87]            = Vec3(  4.7779490e-01,   9.0242640e-01,    8.2515620e-01);
    positions[88]            = Vec3(  4.3957480e-01,   8.0786830e-01,    8.2489220e-01);
    positions[89]            = Vec3(  4.6833310e-01,   9.2867710e-01,    9.1788160e-01);
    positions[90]            = Vec3(  8.2204140e-01,   9.0145630e-01,   -2.5081510e-01);
    positions[91]            = Vec3(  8.5191840e-01,   8.1397830e-01,   -2.2168590e-01);
    positions[92]            = Vec3(  7.6397810e-01,   9.2011290e-01,   -1.8137750e-01);
    positions[93]            = Vec3( -7.9443650e-01,   1.7601300e-01,    4.6436790e-01);
    positions[94]            = Vec3( -7.9212150e-01,   2.3533020e-01,    3.8657500e-01);
    positions[95]            = Vec3( -8.7057070e-01,   1.1288830e-01,    4.4595260e-01);
    positions[96]            = Vec3(  3.2425690e-01,   3.8214720e-01,   -8.2471120e-01);
    positions[97]            = Vec3(  2.8321830e-01,   4.2912450e-01,   -7.4875880e-01);
    positions[98]            = Vec3(  2.7681870e-01,   2.9837230e-01,   -8.2620080e-01);
    positions[99]            = Vec3(  7.5575820e-01,  -8.9620900e-01,    2.3680670e-01);
    positions[100]           = Vec3(  6.6600420e-01,  -8.7027760e-01,    2.7104280e-01);
    positions[101]           = Vec3(  8.1544110e-01,  -9.1190240e-01,    3.1149610e-01);
    positions[102]           = Vec3( -8.4248740e-01,   3.5007110e-01,   -4.4389740e-01);
    positions[103]           = Vec3( -7.5693800e-01,   3.9510690e-01,   -4.4710480e-01);
    positions[104]           = Vec3( -8.6984880e-01,   3.5457480e-01,   -5.3702920e-01);
    positions[105]           = Vec3(  3.8837250e-01,  -4.8496240e-01,    6.5322550e-01);
    positions[106]           = Vec3(  4.1237110e-01,  -4.0401080e-01,    7.0255980e-01);
    positions[107]           = Vec3(  3.0065040e-01,  -4.6399160e-01,    6.0513040e-01);
    positions[108]           = Vec3(  6.2063930e-01,  -5.0831230e-01,    4.9540430e-01);
    positions[109]           = Vec3(  6.8959700e-01,  -5.3506820e-01,    5.6328860e-01);
    positions[110]           = Vec3(  5.3663630e-01,  -5.1121830e-01,    5.4640900e-01);
    positions[111]           = Vec3(  7.0354670e-01,  -5.1748580e-01,   -7.3878700e-02);
    positions[112]           = Vec3(  7.8529450e-01,  -5.6535940e-01,   -9.5943500e-02);
    positions[113]           = Vec3(  6.7807440e-01,  -4.7921810e-01,   -1.6187590e-01);
    positions[114]           = Vec3( -4.4116790e-01,  -4.7749880e-01,    3.0876830e-01);
    positions[115]           = Vec3( -5.0645290e-01,  -4.1075220e-01,    3.1159470e-01);
    positions[116]           = Vec3( -4.6594720e-01,  -5.2568230e-01,    3.8755370e-01);
    positions[117]           = Vec3( -9.1937480e-01,  -5.8400000e-05,   -2.5359570e-01);
    positions[118]           = Vec3( -8.5894750e-01,  -7.0402500e-02,   -2.2230370e-01);
    positions[119]           = Vec3( -8.7441760e-01,   8.3170500e-02,   -2.3447490e-01);
    positions[120]           = Vec3(  5.0867290e-01,   2.3568780e-01,    5.5935510e-01);
    positions[121]           = Vec3(  4.1446460e-01,   2.6088930e-01,    5.8683440e-01);
    positions[122]           = Vec3(  5.1853820e-01,   1.4937830e-01,    5.8561390e-01);
    positions[123]           = Vec3( -4.6831090e-01,  -6.1465890e-01,   -1.6794620e-01);
    positions[124]           = Vec3( -4.8688540e-01,  -5.9611250e-01,   -7.4636500e-02);
    positions[125]           = Vec3( -4.9162010e-01,  -7.0497770e-01,   -1.8127910e-01);
    positions[126]           = Vec3( -3.1791800e-01,  -5.4450000e-03,   -3.6397680e-01);
    positions[127]           = Vec3( -2.2253910e-01,  -2.4457600e-02,   -3.5240990e-01);
    positions[128]           = Vec3( -3.6044390e-01,  -3.5065000e-02,   -2.8414310e-01);
    positions[129]           = Vec3(  1.0461140e-01,   2.6758700e-01,   -2.2684050e-01);
    positions[130]           = Vec3(  1.8426490e-01,   3.2453330e-01,   -2.3574350e-01);
    positions[131]           = Vec3(  1.0569370e-01,   2.3628020e-01,   -1.3834830e-01);
    positions[132]           = Vec3( -1.4119340e-01,   4.1653970e-01,   -2.7320250e-01);
    positions[133]           = Vec3( -5.2065100e-02,   3.6979030e-01,   -2.6662970e-01);
    positions[134]           = Vec3( -1.3834110e-01,   4.7690560e-01,   -1.9435870e-01);
    positions[135]           = Vec3( -7.6602450e-01,  -2.1216400e-01,   -1.9516640e-01);
    positions[136]           = Vec3( -8.0191290e-01,  -2.8391260e-01,   -1.3557910e-01);
    positions[137]           = Vec3( -7.4415500e-01,  -2.6044280e-01,   -2.8169590e-01);
    positions[138]           = Vec3( -1.3600310e-01,   1.9674000e-01,    2.0349610e-01);
    positions[139]           = Vec3( -1.6201050e-01,   2.8693750e-01,    2.3123820e-01);
    positions[140]           = Vec3( -2.1785650e-01,   1.4514420e-01,    1.9201990e-01);
    positions[141]           = Vec3(  6.2897820e-01,  -4.2302590e-01,   -7.6557210e-01);
    positions[142]           = Vec3(  6.2334100e-01,  -4.4471660e-01,   -6.7174140e-01);
    positions[143]           = Vec3(  6.3346670e-01,  -5.0696850e-01,   -8.0495300e-01);
    positions[144]           = Vec3(  9.1588260e-01,  -3.9845200e-02,    3.5189180e-01);
    positions[145]           = Vec3(  9.8891550e-01,  -4.4673900e-02,    2.9156120e-01);
    positions[146]           = Vec3(  8.4126090e-01,  -2.2841000e-03,    2.9707980e-01);
    positions[147]           = Vec3(  4.8470900e-01,  -8.2561400e-02,    6.0082980e-01);
    positions[148]           = Vec3(  3.9021850e-01,  -6.2932500e-02,    6.0195610e-01);
    positions[149]           = Vec3(  5.0563070e-01,  -7.9866200e-02,    5.0777230e-01);
    positions[150]           = Vec3( -7.2845180e-01,  -3.4650580e-01,    7.5973620e-01);
    positions[151]           = Vec3( -7.6073760e-01,  -3.6974690e-01,    6.7323450e-01);
    positions[152]           = Vec3( -7.1326740e-01,  -2.4916760e-01,    7.5651020e-01);
    positions[153]           = Vec3( -3.0896820e-01,  -3.8029640e-01,    6.5520670e-01);
    positions[154]           = Vec3( -3.5019560e-01,  -4.5571260e-01,    6.1040330e-01);
    positions[155]           = Vec3( -2.8479430e-01,  -3.2175460e-01,    5.7933340e-01);
    positions[156]           = Vec3( -6.2826700e-02,  -6.4315900e-02,   -6.8812300e-02);
    positions[157]           = Vec3( -7.4971500e-02,   1.9900000e-02,   -1.8191100e-02);
    positions[158]           = Vec3(  3.2478400e-02,  -8.8932300e-02,   -5.6413600e-02);
    positions[159]           = Vec3(  1.1667520e-01,  -6.6784990e-01,    1.1452860e-01);
    positions[160]           = Vec3(  6.5194200e-02,  -7.3080350e-01,    6.5294000e-02);
    positions[161]           = Vec3(  1.6133150e-01,  -6.1778770e-01,    4.7196600e-02);
    positions[162]           = Vec3(  8.8627400e-02,  -7.1850240e-01,    3.7581390e-01);
    positions[163]           = Vec3(  1.2356120e-01,  -8.0690930e-01,    3.7094210e-01);
    positions[164]           = Vec3(  9.2028600e-02,  -6.8313750e-01,    2.8412340e-01);
    positions[165]           = Vec3(  2.1347270e-01,   8.4107000e-03,    6.0413030e-01);
    positions[166]           = Vec3(  1.8845570e-01,   4.3251500e-02,    5.1600410e-01);
    positions[167]           = Vec3(  1.6789670e-01,  -7.6656800e-02,    6.0793520e-01);
    positions[168]           = Vec3(  1.8425700e-02,   3.0164400e-02,    8.4213210e-01);
    positions[169]           = Vec3( -7.1641800e-02,   4.1848500e-02,    8.7065260e-01);
    positions[170]           = Vec3(  4.4510400e-02,  -6.2982500e-02,    8.6373290e-01);
    positions[171]           = Vec3( -3.1486750e-01,  -1.9966860e-01,   -5.7954700e-01);
    positions[172]           = Vec3( -3.2321140e-01,  -1.4613590e-01,   -5.0133480e-01);
    positions[173]           = Vec3( -3.4769180e-01,  -1.4810900e-01,   -6.5567720e-01);
    positions[174]           = Vec3(  2.2013690e-01,  -4.8207100e-02,   -6.6169910e-01);
    positions[175]           = Vec3(  1.3676160e-01,  -9.4600100e-02,   -6.4525960e-01);
    positions[176]           = Vec3(  2.7051720e-01,  -1.2158460e-01,   -6.9535940e-01);
    positions[177]           = Vec3( -1.5721060e-01,  -2.0015580e-01,    4.8442010e-01);
    positions[178]           = Vec3( -7.4675400e-02,  -2.0952300e-01,    5.3560160e-01);
    positions[179]           = Vec3( -1.8522760e-01,  -1.0781560e-01,    5.0024110e-01);
    positions[180]           = Vec3(  5.4002730e-01,   6.3800500e-01,   -8.0040500e-01);
    positions[181]           = Vec3(  5.0366070e-01,   7.1545920e-01,   -7.5257350e-01);
    positions[182]           = Vec3(  5.1480770e-01,   5.5941670e-01,   -7.4903220e-01);
    positions[183]           = Vec3( -6.3383580e-01,   5.7282910e-01,   -1.7429980e-01);
    positions[184]           = Vec3( -6.0668100e-01,   4.7712900e-01,   -1.7677570e-01);
    positions[185]           = Vec3( -5.6638740e-01,   6.1288510e-01,   -2.2951390e-01);
    positions[186]           = Vec3( -2.0998170e-01,  -2.7747820e-01,    7.0579400e-02);
    positions[187]           = Vec3( -1.4055440e-01,  -3.0201380e-01,    1.3644740e-01);
    positions[188]           = Vec3( -1.6881700e-01,  -2.1818660e-01,    6.9733000e-03);
    positions[189]           = Vec3( -7.6400000e-04,   5.6326380e-01,    1.4175360e-01);
    positions[190]           = Vec3( -7.3688000e-02,   5.0031150e-01,    1.5514670e-01);
    positions[191]           = Vec3( -2.5553000e-02,   6.4733770e-01,    1.7711800e-01);
    positions[192]           = Vec3(  3.9595890e-01,  -1.9078420e-01,   -1.9708050e-01);
    positions[193]           = Vec3(  4.3887020e-01,  -1.5694200e-01,   -1.1582060e-01);
    positions[194]           = Vec3(  3.7635540e-01,  -1.1834040e-01,   -2.5323660e-01);
    positions[195]           = Vec3(  3.9638900e-02,  -2.4093090e-01,    8.9424300e-01);
    positions[196]           = Vec3( -4.9643600e-02,  -2.7156660e-01,    8.9962920e-01);
    positions[197]           = Vec3(  8.4318200e-02,  -2.7149840e-01,    9.7721820e-01);
    positions[198]           = Vec3( -5.9039370e-01,  -3.5975630e-01,   -7.1984370e-01);
    positions[199]           = Vec3( -5.3914870e-01,  -4.0214860e-01,   -7.8361060e-01);
    positions[200]           = Vec3( -6.8562580e-01,  -3.9051900e-01,   -7.4071320e-01);
    positions[201]           = Vec3(  4.7759800e-01,   3.2863960e-01,   -5.4274200e-02);
    positions[202]           = Vec3(  4.5034450e-01,   3.6680450e-01,   -1.4201230e-01);
    positions[203]           = Vec3(  4.3083410e-01,   3.8043410e-01,    1.5118500e-02);
    positions[204]           = Vec3(  1.8100450e-01,   1.6674000e-01,   -8.4907090e-01);
    positions[205]           = Vec3(  1.0479500e-01,   1.5721720e-01,   -9.0737790e-01);
    positions[206]           = Vec3(  1.7365410e-01,   9.7140100e-02,   -7.7842430e-01);
    positions[207]           = Vec3( -6.9841710e-01,   8.5211760e-01,    4.9956020e-01);
    positions[208]           = Vec3( -6.3194850e-01,   9.0336360e-01,    4.5467020e-01);
    positions[209]           = Vec3( -6.7863830e-01,   7.5666570e-01,    5.1268950e-01);
    positions[210]           = Vec3(  8.0356880e-01,  -7.6669620e-01,    5.6240980e-01);
    positions[211]           = Vec3(  8.9444390e-01,  -7.9421520e-01,    5.4379860e-01);
    positions[212]           = Vec3(  8.0061200e-01,  -7.1151420e-01,    6.3743510e-01);
    positions[213]           = Vec3( -2.3686380e-01,   4.4018650e-01,    2.7494630e-01);
    positions[214]           = Vec3( -2.1006750e-01,   4.1932880e-01,    3.6593160e-01);
    positions[215]           = Vec3( -3.2910900e-01,   4.6299420e-01,    2.7725190e-01);
    positions[216]           = Vec3(  7.3324180e-01,   9.1021100e-02,    8.6347740e-01);
    positions[217]           = Vec3(  6.4934460e-01,   5.3444800e-02,    8.7843600e-01);
    positions[218]           = Vec3(  7.1407590e-01,   1.8691830e-01,    8.6323690e-01);
    positions[219]           = Vec3(  3.6906600e-02,   1.4742360e-01,    4.0082880e-01);
    positions[220]           = Vec3( -1.0515300e-02,   1.4450010e-01,    4.8531790e-01);
    positions[221]           = Vec3( -3.6861400e-02,   1.5333190e-01,    3.3364650e-01);
    positions[222]           = Vec3(  5.7666790e-01,  -9.2075640e-01,    5.7305300e-01);
    positions[223]           = Vec3(  5.4452540e-01,  -9.3954290e-01,    6.5798160e-01);
    positions[224]           = Vec3(  6.7020160e-01,  -8.8052280e-01,    5.6852240e-01);
    positions[225]           = Vec3(  4.1616300e-01,  -2.3723450e-01,    7.8105700e-02);
    positions[226]           = Vec3(  4.4947640e-01,  -2.2465620e-01,    1.6469280e-01);
    positions[227]           = Vec3(  3.6093380e-01,  -3.1332780e-01,    7.1125100e-02);
    positions[228]           = Vec3( -1.9830990e-01,  -6.8678560e-01,   -7.6648560e-01);
    positions[229]           = Vec3( -1.1489950e-01,  -6.8356660e-01,   -8.2028210e-01);
    positions[230]           = Vec3( -2.0935090e-01,  -5.9618710e-01,   -7.3178710e-01);
    positions[231]           = Vec3( -4.3741650e-01,  -7.8889500e-01,    1.7785560e-01);
    positions[232]           = Vec3( -3.6424030e-01,  -7.2995610e-01,    1.5380490e-01);
    positions[233]           = Vec3( -5.0710310e-01,  -7.4066850e-01,    1.3917790e-01);
    positions[234]           = Vec3(  5.1605280e-01,   6.8521860e-01,    4.5545030e-01);
    positions[235]           = Vec3(  5.3920960e-01,   7.6750670e-01,    4.8965960e-01);
    positions[236]           = Vec3(  5.4441350e-01,   6.8153880e-01,    3.6305340e-01);
    positions[237]           = Vec3( -9.1377180e-01,   9.0412110e-01,   -8.0577110e-01);
    positions[238]           = Vec3( -8.6299150e-01,   9.8552780e-01,   -7.9463610e-01);
    positions[239]           = Vec3( -9.1270510e-01,   8.7715830e-01,   -9.0107170e-01);
    positions[240]           = Vec3( -5.6874630e-01,  -3.9330600e-02,    5.3540210e-01);
    positions[241]           = Vec3( -6.0667690e-01,   3.6619200e-02,    4.9922460e-01);
    positions[242]           = Vec3( -5.8307630e-01,  -4.4694300e-02,    6.3380260e-01);
    positions[243]           = Vec3(  1.0312020e-01,   2.2809180e-01,    5.7525600e-02);
    positions[244]           = Vec3(  1.8161800e-02,   2.2164820e-01,    1.0293620e-01);
    positions[245]           = Vec3(  1.4691520e-01,   3.0734480e-01,    9.4432600e-02);
    positions[246]           = Vec3( -5.3437690e-01,  -9.0689060e-01,   -7.7012560e-01);
    positions[247]           = Vec3( -6.0761130e-01,  -8.5593580e-01,   -8.0463440e-01);
    positions[248]           = Vec3( -5.5313680e-01,  -9.9745020e-01,   -8.0224750e-01);
    positions[249]           = Vec3(  1.7436730e-01,  -4.6935620e-01,   -7.7408150e-01);
    positions[250]           = Vec3(  1.3315640e-01,  -4.6856170e-01,   -6.8363440e-01);
    positions[251]           = Vec3(  2.3486700e-01,  -3.9970620e-01,   -7.7872930e-01);
    positions[252]           = Vec3(  5.0382310e-01,   8.6391330e-01,   -6.1751380e-01);
    positions[253]           = Vec3(  5.7851670e-01,   9.1774780e-01,   -6.3741940e-01);
    positions[254]           = Vec3(  5.2100060e-01,   8.2278060e-01,   -5.3449130e-01);
    positions[255]           = Vec3( -2.3461000e-03,   8.8439120e-01,   -3.5703750e-01);
    positions[256]           = Vec3( -4.5869800e-02,   9.2025060e-01,   -4.4264850e-01);
    positions[257]           = Vec3(  7.7568300e-02,   8.3812640e-01,   -3.7824790e-01);
    positions[258]           = Vec3( -1.6677150e-01,  -9.0353490e-01,   -5.6323410e-01);
    positions[259]           = Vec3( -1.5077930e-01,  -8.7448310e-01,   -6.5150250e-01);
    positions[260]           = Vec3( -2.5054260e-01,  -8.5746520e-01,   -5.4471400e-01);
    positions[261]           = Vec3( -1.0245710e-01,  -4.1390500e-01,    2.9240710e-01);
    positions[262]           = Vec3( -1.6375100e-01,  -3.5806090e-01,    3.3803800e-01);
    positions[263]           = Vec3( -3.4371600e-02,  -4.4188880e-01,    3.6032470e-01);
    positions[264]           = Vec3(  6.7721230e-01,  -9.2755000e-01,   -6.1695000e-03);
    positions[265]           = Vec3(  6.3209610e-01,  -8.4066740e-01,    1.5854000e-03);
    positions[266]           = Vec3(  7.2195780e-01,  -9.3506790e-01,    7.6821700e-02);
    positions[267]           = Vec3( -5.2597410e-01,   5.0741940e-01,    2.8142130e-01);
    positions[268]           = Vec3( -5.3172740e-01,   5.6506650e-01,    2.0013640e-01);
    positions[269]           = Vec3( -5.9533220e-01,   4.4193270e-01,    2.6673520e-01);
    positions[270]           = Vec3(  4.3852700e-02,  -7.1092730e-01,   -3.0056810e-01);
    positions[271]           = Vec3(  1.2232900e-02,  -6.7601300e-01,   -2.1679320e-01);
    positions[272]           = Vec3(  3.0039200e-02,  -8.0474130e-01,   -3.0050550e-01);
    positions[273]           = Vec3(  4.7537430e-01,   6.7956000e-03,   -8.8926760e-01);
    positions[274]           = Vec3(  4.4972180e-01,  -8.1937800e-02,   -8.5037740e-01);
    positions[275]           = Vec3(  3.9238110e-01,   5.4650000e-02,   -9.0978500e-01);
    positions[276]           = Vec3(  8.4526190e-01,  -3.2384610e-01,    4.4702430e-01);
    positions[277]           = Vec3(  8.5335920e-01,  -2.3860050e-01,    4.1507690e-01);
    positions[278]           = Vec3(  9.3799800e-01,  -3.6222940e-01,    4.5249690e-01);
    positions[279]           = Vec3( -8.5624140e-01,  -3.3540460e-01,   -5.3955060e-01);
    positions[280]           = Vec3( -8.9833150e-01,  -3.2177130e-01,   -6.2636700e-01);
    positions[281]           = Vec3( -7.6568080e-01,  -3.0076830e-01,   -5.3672910e-01);
    positions[282]           = Vec3( -4.0866080e-01,  -7.0070860e-01,    9.2586930e-01);
    positions[283]           = Vec3( -4.5043520e-01,  -7.7640050e-01,    9.7012510e-01);
    positions[284]           = Vec3( -3.2086210e-01,  -6.9414110e-01,    9.6526100e-01);
    positions[285]           = Vec3( -2.9612090e-01,   2.9021400e-01,   -4.6137730e-01);
    positions[286]           = Vec3( -3.0085180e-01,   1.9752840e-01,   -4.3159520e-01);
    positions[287]           = Vec3( -2.4502340e-01,   3.3756140e-01,   -3.9070450e-01);
    positions[288]           = Vec3( -8.4956240e-01,  -3.3051010e-01,    4.2215900e-02);
    positions[289]           = Vec3( -8.2077940e-01,  -3.9086690e-01,    1.1548590e-01);
    positions[290]           = Vec3( -9.3822180e-01,  -3.1618550e-01,    5.2894000e-02);
    positions[291]           = Vec3( -8.6464030e-01,   7.5345250e-01,    1.9545370e-01);
    positions[292]           = Vec3( -9.2073720e-01,   6.7584430e-01,    1.8998460e-01);
    positions[293]           = Vec3( -8.9310500e-01,   7.8515510e-01,    2.8077440e-01);
    positions[294]           = Vec3(  5.3248170e-01,   6.8435100e-02,   -1.1431070e-01);
    positions[295]           = Vec3(  6.1630600e-01,   5.7417300e-02,   -6.9794300e-02);
    positions[296]           = Vec3(  4.9275030e-01,   1.5234490e-01,   -7.9235100e-02);
    positions[297]           = Vec3( -3.0166400e-02,   3.6028840e-01,   -9.2023940e-01);
    positions[298]           = Vec3(  2.5390700e-02,   4.3355180e-01,   -9.4581010e-01);
    positions[299]           = Vec3( -1.2837900e-02,   3.5198820e-01,   -8.2331230e-01);
    positions[300]           = Vec3( -7.6094250e-01,  -7.4142570e-01,   -7.6415170e-01);
    positions[301]           = Vec3( -7.5826150e-01,  -6.8315050e-01,   -8.4024930e-01);
    positions[302]           = Vec3( -7.8169550e-01,  -6.8557300e-01,   -6.8728990e-01);
    positions[303]           = Vec3( -7.1618050e-01,  -8.6617600e-02,   -7.8297100e-01);
    positions[304]           = Vec3( -6.9164460e-01,  -1.6643810e-01,   -7.3660090e-01);
    positions[305]           = Vec3( -8.1169890e-01,  -8.6541300e-02,   -7.7380060e-01);
    positions[306]           = Vec3(  8.6280550e-01,  -2.8731190e-01,   -7.5013210e-01);
    positions[307]           = Vec3(  8.4297110e-01,  -2.0142080e-01,   -7.9688520e-01);
    positions[308]           = Vec3(  7.7553640e-01,  -3.3421630e-01,   -7.5754400e-01);
    positions[309]           = Vec3(  2.9607200e-02,  -6.7251560e-01,   -9.1368960e-01);
    positions[310]           = Vec3(  1.0909000e-03,  -6.2708430e-01,   -9.9528360e-01);
    positions[311]           = Vec3(  8.0161300e-02,  -5.9814710e-01,   -8.7106130e-01);
    positions[312]           = Vec3( -2.2829370e-01,   4.6661410e-01,    7.7985190e-01);
    positions[313]           = Vec3( -2.4730820e-01,   5.6404020e-01,    7.8763210e-01);
    positions[314]           = Vec3( -1.7899690e-01,   4.3324110e-01,    8.5622400e-01);
    positions[315]           = Vec3( -5.1323270e-01,  -2.6480150e-01,    7.2113100e-02);
    positions[316]           = Vec3( -4.4180310e-01,  -2.8480730e-01,    1.4166490e-01);
    positions[317]           = Vec3( -5.5826690e-01,  -3.4508980e-01,    5.6782300e-02);
    positions[318]           = Vec3(  3.5970320e-01,  -7.1101700e-01,   -8.5706800e-01);
    positions[319]           = Vec3(  3.5573750e-01,  -6.6123030e-01,   -7.7069560e-01);
    positions[320]           = Vec3(  2.9308100e-01,  -6.7738800e-01,   -9.1162920e-01);
    positions[321]           = Vec3( -8.6077820e-01,  -8.3187420e-01,    3.5264550e-01);
    positions[322]           = Vec3( -7.9919290e-01,  -8.9965630e-01,    3.8875110e-01);
    positions[323]           = Vec3( -8.4377450e-01,  -8.2428940e-01,    2.5657630e-01);
    positions[324]           = Vec3( -6.9407750e-01,   8.5240530e-01,   -4.8975260e-01);
    positions[325]           = Vec3( -6.0369970e-01,   8.2005830e-01,   -4.7948010e-01);
    positions[326]           = Vec3( -6.8257340e-01,   9.0158170e-01,   -5.7057020e-01);
    positions[327]           = Vec3( -8.6181560e-01,   2.1174420e-01,    3.2775000e-02);
    positions[328]           = Vec3( -9.5070390e-01,   2.5868190e-01,    2.6787700e-02);
    positions[329]           = Vec3( -8.8015990e-01,   1.3696510e-01,    9.1486900e-02);
    positions[330]           = Vec3( -6.7034530e-01,  -7.0959980e-01,    5.7197940e-01);
    positions[331]           = Vec3( -6.3447070e-01,  -7.7970770e-01,    5.1435410e-01);
    positions[332]           = Vec3( -7.1147280e-01,  -7.6230200e-01,    6.4084900e-01);
    positions[333]           = Vec3( -4.2433970e-01,   1.6353470e-01,   -7.5364040e-01);
    positions[334]           = Vec3( -3.3715920e-01,   1.3734360e-01,   -7.8660110e-01);
    positions[335]           = Vec3( -4.5203330e-01,   2.3873860e-01,   -8.1607320e-01);
    positions[336]           = Vec3( -4.2091960e-01,  -8.1633330e-01,   -5.3063920e-01);
    positions[337]           = Vec3( -4.2728590e-01,  -7.1806470e-01,   -5.4109270e-01);
    positions[338]           = Vec3( -4.5013260e-01,  -8.3810340e-01,   -6.1998700e-01);
    positions[339]           = Vec3(  6.0367930e-01,   3.3084920e-01,   -8.4465460e-01);
    positions[340]           = Vec3(  5.0455880e-01,   3.3698360e-01,   -8.4011240e-01);
    positions[341]           = Vec3(  6.2487550e-01,   2.4834360e-01,   -8.0607210e-01);
    positions[342]           = Vec3(  1.8546120e-01,  -6.3282200e-02,    5.1304500e-02);
    positions[343]           = Vec3(  2.8101390e-01,  -7.7771500e-02,    5.1163200e-02);
    positions[344]           = Vec3(  1.7127760e-01,   2.0996700e-02,    9.0574100e-02);
    positions[345]           = Vec3( -3.5029200e-02,  -7.9917400e-02,   -3.4468400e-01);
    positions[346]           = Vec3( -6.3903800e-02,  -6.0213300e-02,   -2.5206780e-01);
    positions[347]           = Vec3( -4.6785200e-02,  -1.7349570e-01,   -3.5772680e-01);
    positions[348]           = Vec3(  2.5567190e-01,   6.2355480e-01,    4.2852620e-01);
    positions[349]           = Vec3(  1.9093710e-01,   6.4505930e-01,    4.9102940e-01);
    positions[350]           = Vec3(  3.4540670e-01,   6.4937420e-01,    4.5902510e-01);
    positions[351]           = Vec3( -7.3742490e-01,  -8.7628820e-01,   -2.6411710e-01);
    positions[352]           = Vec3( -7.3220480e-01,  -9.1540050e-01,   -3.5104230e-01);
    positions[353]           = Vec3( -7.9968040e-01,  -9.2863850e-01,   -2.1682500e-01);
    positions[354]           = Vec3(  5.1017210e-01,  -2.7173980e-01,    7.9174500e-01);
    positions[355]           = Vec3(  5.1045830e-01,  -2.0746280e-01,    7.2138780e-01);
    positions[356]           = Vec3(  5.9967910e-01,  -3.0815350e-01,    7.9296320e-01);
    positions[357]           = Vec3(  6.1703300e-02,  -6.0490320e-01,   -5.4304490e-01);
    positions[358]           = Vec3(  6.5202000e-03,  -6.6388800e-01,   -5.9525970e-01);
    positions[359]           = Vec3(  6.2525700e-02,  -6.3466150e-01,   -4.5175130e-01);
    positions[360]           = Vec3( -5.0181950e-01,   6.8138390e-01,   -8.8794760e-01);
    positions[361]           = Vec3( -4.0469720e-01,   6.5541180e-01,   -8.8475300e-01);
    positions[362]           = Vec3( -5.4953810e-01,   6.3245150e-01,   -8.1669610e-01);
    positions[363]           = Vec3( -3.5708340e-01,   8.1787480e-01,    1.0372050e-01);
    positions[364]           = Vec3( -4.3575160e-01,   7.6657380e-01,    8.8357500e-02);
    positions[365]           = Vec3( -3.8126100e-01,   9.1312250e-01,    1.2894930e-01);
    positions[366]           = Vec3( -1.0889180e-01,   6.4289110e-01,   -1.1000150e-01);
    positions[367]           = Vec3( -9.5792300e-02,   6.5121590e-01,   -1.2915400e-02);
    positions[368]           = Vec3( -1.4253020e-01,   7.3532640e-01,   -1.2649680e-01);
    positions[369]           = Vec3( -8.0675190e-01,   3.8993580e-01,   -9.3061890e-01);
    positions[370]           = Vec3( -8.4285770e-01,   4.7693320e-01,   -9.5868770e-01);
    positions[371]           = Vec3( -7.4065520e-01,   4.1059110e-01,   -8.6270860e-01);
    positions[372]           = Vec3( -7.3221050e-01,  -8.3486000e-02,    1.8651540e-01);
    positions[373]           = Vec3( -6.6332990e-01,  -2.5838100e-02,    1.5155080e-01);
    positions[374]           = Vec3( -7.5939010e-01,  -1.4675440e-01,    1.1813700e-01);
    positions[375]           = Vec3(  6.1370510e-01,  -3.7510720e-01,   -2.9444790e-01);
    positions[376]           = Vec3(  5.3141590e-01,  -3.1971250e-01,   -2.8369080e-01);
    positions[377]           = Vec3(  6.7472620e-01,  -3.0544670e-01,   -3.2680390e-01);
    positions[378]           = Vec3(  2.8333090e-01,   7.0116700e-01,    6.3582400e-02);
    positions[379]           = Vec3(  2.3304950e-01,   7.8436370e-01,    8.8113000e-02);
    positions[380]           = Vec3(  2.1603670e-01,   6.3345680e-01,    4.3706900e-02);
    positions[381]           = Vec3(  3.4046290e-01,  -5.8425160e-01,   -5.8383960e-01);
    positions[382]           = Vec3(  4.2396660e-01,  -5.6867730e-01,   -5.4787780e-01);
    positions[383]           = Vec3(  2.7987870e-01,  -5.6273080e-01,   -5.1485370e-01);
    positions[384]           = Vec3(  4.8651200e-01,   3.9384650e-01,   -5.0852640e-01);
    positions[385]           = Vec3(  4.8954070e-01,   2.9830160e-01,   -5.1540010e-01);
    positions[386]           = Vec3(  5.7513360e-01,   4.2777280e-01,   -4.8094980e-01);
    positions[387]           = Vec3( -4.9931530e-01,  -8.6556710e-01,    4.1410020e-01);
    positions[388]           = Vec3( -4.0971070e-01,  -9.0364250e-01,    4.1539320e-01);
    positions[389]           = Vec3( -5.0187830e-01,  -8.1863570e-01,    3.2854240e-01);
    positions[390]           = Vec3( -9.2923250e-01,  -9.5140200e-02,    7.7175180e-01);
    positions[391]           = Vec3( -1.0068535e+00,  -4.9193300e-02,    8.1361050e-01);
    positions[392]           = Vec3( -8.5382270e-01,  -3.5167000e-02,    7.7988780e-01);
    positions[393]           = Vec3(  5.8200510e-01,  -2.7347380e-01,    3.2175080e-01);
    positions[394]           = Vec3(  5.9114530e-01,  -2.1232990e-01,    3.9188270e-01);
    positions[395]           = Vec3(  6.2697690e-01,  -3.5436570e-01,    3.5518080e-01);
    positions[396]           = Vec3( -4.3869270e-01,   7.1030180e-01,   -3.4435510e-01);
    positions[397]           = Vec3( -3.5798370e-01,   6.6801330e-01,   -3.8293170e-01);
    positions[398]           = Vec3( -3.9584820e-01,   7.8582280e-01,   -3.0015890e-01);
    positions[399]           = Vec3(  3.0315060e-01,   2.0553140e-01,    3.3518590e-01);
    positions[400]           = Vec3(  2.0466680e-01,   2.0029920e-01,    3.3800050e-01);
    positions[401]           = Vec3(  3.1784090e-01,   2.6138240e-01,    4.0966770e-01);
    positions[402]           = Vec3(  7.3144120e-01,   1.1861840e-01,    2.1872590e-01);
    positions[403]           = Vec3(  6.9245610e-01,   2.0755440e-01,    2.3848660e-01);
    positions[404]           = Vec3(  7.4250960e-01,   1.1063670e-01,    1.1673060e-01);
    positions[405]           = Vec3(  3.0774670e-01,  -6.7782260e-01,   -6.9330000e-02);
    positions[406]           = Vec3(  3.0161020e-01,  -7.5652530e-01,   -1.2627210e-01);
    positions[407]           = Vec3(  3.7612340e-01,  -6.9199170e-01,   -2.0688000e-03);
    positions[408]           = Vec3(  4.8241200e-02,   1.4991530e-01,   -4.8562930e-01);
    positions[409]           = Vec3(  7.0825700e-02,   1.7883510e-01,   -4.0076820e-01);
    positions[410]           = Vec3( -1.4581300e-02,   7.7868400e-02,   -4.8044320e-01);
    positions[411]           = Vec3(  2.6566210e-01,  -4.7972300e-02,   -3.9240060e-01);
    positions[412]           = Vec3(  2.5708940e-01,  -2.6958700e-02,   -4.8906580e-01);
    positions[413]           = Vec3(  1.8079360e-01,  -1.7099600e-02,   -3.5945650e-01);
    positions[414]           = Vec3(  7.3593670e-01,   3.2192010e-01,    6.3185000e-03);
    positions[415]           = Vec3(  7.5313070e-01,   3.1236830e-01,   -9.0780600e-02);
    positions[416]           = Vec3(  6.4125230e-01,   3.3242850e-01,    5.3072000e-03);
    positions[417]           = Vec3( -7.2074000e-03,  -2.1935180e-01,   -6.7044710e-01);
    positions[418]           = Vec3( -7.9916200e-02,  -2.2604130e-01,   -7.3330810e-01);
    positions[419]           = Vec3( -3.2871000e-03,  -2.9557560e-01,   -6.1702790e-01);
    positions[420]           = Vec3(  8.0182800e-01,   3.3340310e-01,   -2.5836160e-01);
    positions[421]           = Vec3(  8.9266890e-01,   3.1760310e-01,   -2.9990300e-01);
    positions[422]           = Vec3(  7.7135080e-01,   4.0881250e-01,   -3.1490320e-01);
    positions[423]           = Vec3( -3.1753700e-01,   3.7248900e-02,    5.0846140e-01);
    positions[424]           = Vec3( -3.3276340e-01,   1.2794660e-01,    5.4135580e-01);
    positions[425]           = Vec3( -4.0442920e-01,  -2.1535000e-03,    5.2164500e-01);
    positions[426]           = Vec3(  7.7089090e-01,  -1.7749490e-01,   -4.1090550e-01);
    positions[427]           = Vec3(  8.0919970e-01,  -9.9267700e-02,   -3.6080690e-01);
    positions[428]           = Vec3(  8.4794900e-01,  -2.2265030e-01,   -4.4286640e-01);
    positions[429]           = Vec3( -5.0985980e-01,   6.5271910e-01,    5.1660950e-01);
    positions[430]           = Vec3( -4.1891080e-01,   6.9500010e-01,    5.0933000e-01);
    positions[431]           = Vec3( -5.2072650e-01,   6.0609800e-01,    4.2889530e-01);
    positions[432]           = Vec3(  8.8931480e-01,  -1.5854900e-02,   -7.9057690e-01);
    positions[433]           = Vec3(  8.4049130e-01,   2.2454500e-02,   -7.1223150e-01);
    positions[434]           = Vec3(  8.6392620e-01,   4.6002000e-02,   -8.5696830e-01);
    positions[435]           = Vec3( -4.2632820e-01,  -5.4538160e-01,   -5.2698140e-01);
    positions[436]           = Vec3( -3.4047810e-01,  -5.2088280e-01,   -5.5637760e-01);
    positions[437]           = Vec3( -4.9107950e-01,  -5.2513960e-01,   -5.9520410e-01);
    positions[438]           = Vec3(  8.8830700e-01,   7.8506050e-01,    4.7420010e-01);
    positions[439]           = Vec3(  9.6737760e-01,   8.0796480e-01,    5.2210120e-01);
    positions[440]           = Vec3(  8.3449840e-01,   7.2694370e-01,    5.2968560e-01);
    positions[441]           = Vec3( -3.0889500e-02,  -5.4040860e-01,   -7.7446500e-02);
    positions[442]           = Vec3(  2.4910200e-02,  -4.7046460e-01,   -5.3187100e-02);
    positions[443]           = Vec3( -1.0937030e-01,  -5.1212170e-01,   -1.2642620e-01);
    positions[444]           = Vec3(  5.0722190e-01,  -8.0898340e-01,    3.3208510e-01);
    positions[445]           = Vec3(  5.1254280e-01,  -8.4333670e-01,    4.2962250e-01);
    positions[446]           = Vec3(  4.8459280e-01,  -7.1548850e-01,    3.3664280e-01);
    positions[447]           = Vec3(  7.0974400e-02,  -8.6268490e-01,   -7.2122900e-01);
    positions[448]           = Vec3(  8.8211100e-02,  -8.1266230e-01,   -7.9698760e-01);
    positions[449]           = Vec3(  1.4856180e-01,  -8.7440360e-01,   -6.6601020e-01);
    positions[450]           = Vec3( -2.7264270e-01,   8.2117820e-01,    4.0979220e-01);
    positions[451]           = Vec3( -1.8893860e-01,   7.8611730e-01,    4.4435560e-01);
    positions[452]           = Vec3( -2.7256440e-01,   8.1557060e-01,    3.0746650e-01);
    positions[453]           = Vec3( -2.3667600e-01,   7.0807760e-01,    9.0055470e-01);
    positions[454]           = Vec3( -1.7087350e-01,   7.0278860e-01,    9.7330650e-01);
    positions[455]           = Vec3( -2.2325560e-01,   8.0596230e-01,    8.7050690e-01);
    positions[456]           = Vec3(  6.0904540e-01,  -5.3471490e-01,   -5.1588800e-01);
    positions[457]           = Vec3(  6.6627390e-01,  -6.1177680e-01,   -4.9309950e-01);
    positions[458]           = Vec3(  6.1303950e-01,  -4.7414890e-01,   -4.3691960e-01);
    positions[459]           = Vec3( -6.9432470e-01,   5.5588670e-01,   -7.2750070e-01);
    positions[460]           = Vec3( -6.8524660e-01,   5.1427650e-01,   -6.4407660e-01);
    positions[461]           = Vec3( -7.7219850e-01,   6.0882800e-01,   -7.1352640e-01);
    positions[462]           = Vec3( -6.5544400e-01,   5.6801890e-01,    7.6654940e-01);
    positions[463]           = Vec3( -5.9853210e-01,   5.8150060e-01,    6.8630620e-01);
    positions[464]           = Vec3( -6.0728400e-01,   6.2604000e-01,    8.2970960e-01);
    positions[465]           = Vec3( -1.7725100e-01,  -7.5128040e-01,    4.8288320e-01);
    positions[466]           = Vec3( -1.1106490e-01,  -7.1604590e-01,    4.2681180e-01);
    positions[467]           = Vec3( -1.2808000e-01,  -8.2063050e-01,    5.2385060e-01);
    positions[468]           = Vec3(  5.0880810e-01,  -1.7782370e-01,   -5.5526690e-01);
    positions[469]           = Vec3(  4.7579150e-01,  -1.6757400e-01,   -4.6732050e-01);
    positions[470]           = Vec3(  6.0010540e-01,  -1.6566020e-01,   -5.4639700e-01);
    positions[471]           = Vec3(  7.9737120e-01,  -5.3326000e-03,   -2.2789800e-02);
    positions[472]           = Vec3(  7.5436910e-01,  -9.2537600e-02,   -2.7176000e-03);
    positions[473]           = Vec3(  8.4035540e-01,  -2.5845500e-02,   -1.0913300e-01);
    positions[474]           = Vec3(  2.4805290e-01,  -4.5182680e-01,   -2.5649240e-01);
    positions[475]           = Vec3(  2.6536400e-01,  -5.1313010e-01,   -1.8699050e-01);
    positions[476]           = Vec3(  2.8661880e-01,  -3.6531040e-01,   -2.2184290e-01);
    positions[477]           = Vec3(  8.9407190e-01,   6.4140150e-01,   -2.2838520e-01);
    positions[478]           = Vec3(  8.6394270e-01,   5.7649930e-01,   -2.9124340e-01);
    positions[479]           = Vec3(  9.8698980e-01,   6.3685520e-01,   -2.2087390e-01);
    positions[480]           = Vec3( -5.0297400e-01,   3.8595440e-01,   -9.1329410e-01);
    positions[481]           = Vec3( -4.2761710e-01,   4.4573350e-01,   -8.9961960e-01);
    positions[482]           = Vec3( -5.7109730e-01,   4.3620760e-01,   -9.6075240e-01);
    positions[483]           = Vec3( -5.7912630e-01,   3.1473530e-01,   -1.3174480e-01);
    positions[484]           = Vec3( -6.4359930e-01,   2.3775370e-01,   -1.3815260e-01);
    positions[485]           = Vec3( -5.1910350e-01,   2.9841740e-01,   -5.0386200e-02);
    positions[486]           = Vec3( -2.3287450e-01,  -4.5325250e-01,   -2.6295780e-01);
    positions[487]           = Vec3( -3.1705790e-01,  -5.0582880e-01,   -2.5755610e-01);
    positions[488]           = Vec3( -2.4988940e-01,  -3.6717760e-01,   -2.2456790e-01);
    positions[489]           = Vec3(  7.3902040e-01,   6.0596960e-01,    8.7531410e-01);
    positions[490]           = Vec3(  6.8510920e-01,   5.6584400e-01,    8.0394270e-01);
    positions[491]           = Vec3(  6.7885220e-01,   6.2492120e-01,    9.4854880e-01);
    positions[492]           = Vec3( -1.7342650e-01,  -4.4833620e-01,   -6.2689720e-01);
    positions[493]           = Vec3( -1.2120170e-01,  -4.7044370e-01,   -5.5178260e-01);
    positions[494]           = Vec3( -2.1756220e-01,  -3.7095040e-01,   -5.9046920e-01);
    positions[495]           = Vec3( -9.7909800e-02,   4.1047140e-01,    5.5154950e-01);
    positions[496]           = Vec3( -1.5085520e-01,   4.3078300e-01,    6.2923170e-01);
    positions[497]           = Vec3( -2.2121000e-02,   3.5854830e-01,    5.8628920e-01);
    positions[498]           = Vec3( -2.9249090e-01,   4.2502460e-01,   -7.0552100e-01);
    positions[499]           = Vec3( -2.3858950e-01,   3.8279980e-01,   -7.7129020e-01);
    positions[500]           = Vec3( -2.8537830e-01,   3.6194940e-01,   -6.3036240e-01);
    positions[501]           = Vec3( -4.1927200e-01,  -1.0765570e-01,   -8.1010100e-01);
    positions[502]           = Vec3( -4.5513170e-01,  -1.8389200e-01,   -8.5055730e-01);
    positions[503]           = Vec3( -4.6978720e-01,  -3.2915400e-02,   -8.4249770e-01);
    positions[504]           = Vec3( -8.3022800e-01,  -5.9366610e-01,   -5.2440890e-01);
    positions[505]           = Vec3( -8.3569020e-01,  -5.0053960e-01,   -5.4596070e-01);
    positions[506]           = Vec3( -7.7653500e-01,  -5.9680800e-01,   -4.4872510e-01);
    positions[507]           = Vec3(  4.7451900e-02,   2.4985900e-01,    7.1027380e-01);
    positions[508]           = Vec3(  5.2750000e-03,   2.6682820e-01,    8.0047760e-01);
    positions[509]           = Vec3(  9.2790500e-02,   1.6390540e-01,    7.2751450e-01);
    positions[510]           = Vec3(  9.8318300e-02,  -2.4834430e-01,    6.2217110e-01);
    positions[511]           = Vec3(  7.1376800e-02,  -2.3868900e-01,    7.1029050e-01);
    positions[512]           = Vec3(  1.0725160e-01,  -3.3946690e-01,    5.9525570e-01);
    positions[513]           = Vec3( -1.7389390e-01,   6.3857050e-01,   -4.3802350e-01);
    positions[514]           = Vec3( -1.0857550e-01,   7.0876020e-01,   -4.2023360e-01);
    positions[515]           = Vec3( -1.6180390e-01,   5.6775180e-01,   -3.7084920e-01);
    positions[516]           = Vec3( -8.3384410e-01,  -7.8320210e-01,    7.9714340e-01);
    positions[517]           = Vec3( -8.6597850e-01,  -8.7176550e-01,    7.7689410e-01);
    positions[518]           = Vec3( -9.1332720e-01,  -7.1912210e-01,    7.9807020e-01);
    positions[519]           = Vec3(  3.0122650e-01,   4.4099240e-01,    1.7747380e-01);
    positions[520]           = Vec3(  3.0879580e-01,   3.5962220e-01,    2.2668340e-01);
    positions[521]           = Vec3(  2.9198270e-01,   5.0655710e-01,    2.4655760e-01);
    positions[522]           = Vec3(  3.8346200e-01,  -2.8443150e-01,   -8.3961770e-01);
    positions[523]           = Vec3(  4.1227770e-01,  -2.9408340e-01,   -9.3409110e-01);
    positions[524]           = Vec3(  4.5498420e-01,  -3.3552520e-01,   -7.8643110e-01);
    positions[525]           = Vec3(  5.4535540e-01,   1.2249720e-01,   -4.0869350e-01);
    positions[526]           = Vec3(  6.0755050e-01,   1.6343320e-01,   -3.4805580e-01);
    positions[527]           = Vec3(  4.8362230e-01,   8.8573600e-02,   -3.4405000e-01);
    positions[528]           = Vec3(  1.3637990e-01,  -3.3186850e-01,    1.0338270e-01);
    positions[529]           = Vec3(  1.5761460e-01,  -2.5187340e-01,    1.5683210e-01);
    positions[530]           = Vec3(  7.8556700e-02,  -3.8461200e-01,    1.6118390e-01);
    positions[531]           = Vec3(  8.4245020e-01,   3.8084570e-01,   -6.9184990e-01);
    positions[532]           = Vec3(  9.0750590e-01,   3.9283710e-01,   -7.7288830e-01);
    positions[533]           = Vec3(  7.5053500e-01,   3.8878480e-01,   -7.2751780e-01);
    positions[534]           = Vec3(  2.7768360e-01,  -8.5899240e-01,   -5.3138620e-01);
    positions[535]           = Vec3(  2.8386750e-01,  -7.7018020e-01,   -5.6323660e-01);
    positions[536]           = Vec3(  3.4891330e-01,  -9.1242960e-01,   -5.6853820e-01);
    positions[537]           = Vec3(  2.6823810e-01,  -7.8504070e-01,    6.9926380e-01);
    positions[538]           = Vec3(  3.3824260e-01,  -8.3764610e-01,    7.3839250e-01);
    positions[539]           = Vec3(  3.0089590e-01,  -6.9098950e-01,    7.0290360e-01);
    positions[540]           = Vec3(  9.5946000e-02,   5.9757730e-01,    8.8417370e-01);
    positions[541]           = Vec3(  1.9084960e-01,   5.8892180e-01,    8.6811780e-01);
    positions[542]           = Vec3(  7.0090900e-02,   6.3001980e-01,    9.7622150e-01);
    positions[543]           = Vec3( -3.2687830e-01,  -9.5478000e-03,    2.1684540e-01);
    positions[544]           = Vec3( -3.2605730e-01,  -1.4225700e-02,    3.1463820e-01);
    positions[545]           = Vec3( -2.7582100e-01,  -8.4479800e-02,    1.8809180e-01);
    positions[546]           = Vec3( -8.3433230e-01,  -5.5202940e-01,    2.1864880e-01);
    positions[547]           = Vec3( -8.2396710e-01,  -5.4694370e-01,    3.1266070e-01);
    positions[548]           = Vec3( -9.3264700e-01,  -5.5452020e-01,    1.9359510e-01);
    positions[549]           = Vec3( -3.7479050e-01,   2.2505660e-01,    7.1205330e-01);
    positions[550]           = Vec3( -4.7509020e-01,   2.3675960e-01,    7.1906840e-01);
    positions[551]           = Vec3( -3.3344270e-01,   3.0911900e-01,    7.2096390e-01);
    positions[552]           = Vec3(  5.4909720e-01,  -6.8048160e-01,    7.2400200e-02);
    positions[553]           = Vec3(  6.0527360e-01,  -6.6696760e-01,    1.5170900e-01);
    positions[554]           = Vec3(  5.8614280e-01,  -6.1178520e-01,    1.7524700e-02);
    positions[555]           = Vec3( -2.3127640e-01,   9.0287820e-01,   -1.3411380e-01);
    positions[556]           = Vec3( -2.8615520e-01,   9.5668910e-01,   -1.9830460e-01);
    positions[557]           = Vec3( -2.9306830e-01,   8.7146310e-01,   -6.8234400e-02);
    positions[558]           = Vec3( -5.4794480e-01,   6.9927600e-02,    4.9211700e-02);
    positions[559]           = Vec3( -4.8467110e-01,   1.3673600e-02,    9.6662900e-02);
    positions[560]           = Vec3( -5.5944570e-01,   3.5041600e-02,   -4.0422400e-02);
    positions[561]           = Vec3( -4.0842490e-01,  -6.1610810e-01,    5.3013490e-01);
    positions[562]           = Vec3( -3.5055240e-01,  -6.7988460e-01,    4.9398580e-01);
    positions[563]           = Vec3( -4.6296070e-01,  -6.7880320e-01,    5.8633470e-01);
    positions[564]           = Vec3(  4.6585780e-01,   7.8746100e-01,   -1.2817710e-01);
    positions[565]           = Vec3(  5.3858490e-01,   8.3094890e-01,   -7.7410200e-02);
    positions[566]           = Vec3(  4.0552000e-01,   7.4979180e-01,   -6.1891900e-02);
    positions[567]           = Vec3( -1.6560700e-02,  -3.7062430e-01,   -3.6569060e-01);
    positions[568]           = Vec3( -9.0792700e-02,  -4.1378610e-01,   -3.2720710e-01);
    positions[569]           = Vec3(  5.1374900e-02,  -4.3774530e-01,   -3.4403280e-01);
    positions[570]           = Vec3( -5.9512760e-01,   1.7073000e-02,   -2.2772060e-01);
    positions[571]           = Vec3( -5.8225940e-01,   7.1421900e-02,   -3.0604790e-01);
    positions[572]           = Vec3( -6.2819960e-01,  -6.3276600e-02,   -2.6202260e-01);
    positions[573]           = Vec3( -4.7641750e-01,  -4.2323550e-01,    8.9604240e-01);
    positions[574]           = Vec3( -5.4796980e-01,  -4.1341290e-01,    8.3129580e-01);
    positions[575]           = Vec3( -4.6422920e-01,  -5.2061790e-01,    8.9834640e-01);
    positions[576]           = Vec3(  1.4489300e-02,  -8.9340740e-01,   -3.4831200e-02);
    positions[577]           = Vec3( -6.9252500e-02,  -9.1064710e-01,   -8.0576000e-02);
    positions[578]           = Vec3(  8.8146500e-02,  -9.1521790e-01,   -9.6184600e-02);
    positions[579]           = Vec3( -6.0237270e-01,   6.8170090e-01,    6.7672100e-02);
    positions[580]           = Vec3( -6.3353490e-01,   6.5944010e-01,   -1.4737000e-02);
    positions[581]           = Vec3( -6.7945440e-01,   7.0260310e-01,    1.2553510e-01);
    positions[582]           = Vec3( -7.9759390e-01,  -4.8566970e-01,   -8.8075620e-01);
    positions[583]           = Vec3( -7.6587590e-01,  -4.4277470e-01,   -9.5861500e-01);
    positions[584]           = Vec3( -8.8262650e-01,  -4.4354590e-01,   -8.6497650e-01);
    positions[585]           = Vec3(  6.0913180e-01,   7.5063640e-01,   -3.7944500e-01);
    positions[586]           = Vec3(  6.8958950e-01,   8.0236210e-01,   -3.5044320e-01);
    positions[587]           = Vec3(  5.5351750e-01,   7.5362410e-01,   -2.9669720e-01);
    positions[588]           = Vec3(  7.4485800e-01,   5.3041050e-01,   -4.4708420e-01);
    positions[589]           = Vec3(  7.0182180e-01,   6.1806940e-01,   -4.3652910e-01);
    positions[590]           = Vec3(  8.0156580e-01,   5.2857300e-01,   -5.2411300e-01);
    positions[591]           = Vec3( -6.9004280e-01,  -5.9012070e-01,   -2.9270410e-01);
    positions[592]           = Vec3( -7.1539690e-01,  -6.8384200e-01,   -2.8572180e-01);
    positions[593]           = Vec3( -5.9319910e-01,  -5.8219810e-01,   -2.7391860e-01);
    positions[594]           = Vec3( -2.0769030e-01,  -9.0263320e-01,    8.2559380e-01);
    positions[595]           = Vec3( -1.2326710e-01,  -9.0347650e-01,    7.7889800e-01);
    positions[596]           = Vec3( -2.4674410e-01,  -8.1114260e-01,    8.1400270e-01);
    positions[597]           = Vec3( -5.9770390e-01,  -2.5353030e-01,    3.7815410e-01);
    positions[598]           = Vec3( -5.7799760e-01,  -1.8503970e-01,    4.3781640e-01);
    positions[599]           = Vec3( -6.3056510e-01,  -1.9169960e-01,    3.0646360e-01);
    positions[600]           = Vec3(  2.5756560e-01,  -9.0983610e-01,   -2.2681580e-01);
    positions[601]           = Vec3(  3.3909840e-01,  -9.6122750e-01,   -2.1952540e-01);
    positions[602]           = Vec3(  2.5286730e-01,  -8.8095350e-01,   -3.1936560e-01);
    positions[603]           = Vec3( -5.7980030e-01,   4.5624440e-01,   -4.8053250e-01);
    positions[604]           = Vec3( -5.0283550e-01,   4.0235700e-01,   -4.7629430e-01);
    positions[605]           = Vec3( -5.5234760e-01,   5.5030960e-01,   -4.6445780e-01);
    positions[606]           = Vec3(  2.6417710e-01,   3.6149920e-01,    5.5726940e-01);
    positions[607]           = Vec3(  1.8510390e-01,   3.4649300e-01,    6.1450570e-01);
    positions[608]           = Vec3(  2.5328370e-01,   4.4988030e-01,    5.1753010e-01);
    positions[609]           = Vec3(  8.0108540e-01,  -7.3935090e-01,   -4.6186460e-01);
    positions[610]           = Vec3(  8.1693510e-01,  -7.9012220e-01,   -3.8198140e-01);
    positions[611]           = Vec3(  8.8304810e-01,  -6.9238110e-01,   -4.8097940e-01);
    positions[612]           = Vec3( -5.8628640e-01,   1.5133800e-02,   -5.2805090e-01);
    positions[613]           = Vec3( -5.0874980e-01,   5.8718200e-02,   -5.5884230e-01);
    positions[614]           = Vec3( -6.4503990e-01,   1.9133400e-02,   -6.0165090e-01);
    positions[615]           = Vec3(  7.6453220e-01,  -5.9994620e-01,    2.8797170e-01);
    positions[616]           = Vec3(  7.0859250e-01,  -5.4012040e-01,    3.3515640e-01);
    positions[617]           = Vec3(  7.9449730e-01,  -6.7260900e-01,    3.4844210e-01);
    positions[618]           = Vec3( -4.1271350e-01,   6.8162960e-01,   -6.2517570e-01);
    positions[619]           = Vec3( -3.4841290e-01,   6.1054470e-01,   -6.4194430e-01);
    positions[620]           = Vec3( -3.5808100e-01,   7.5958210e-01,   -6.0333290e-01);
    positions[621]           = Vec3( -2.3867290e-01,   5.9441400e-02,    9.1386800e-01);
    positions[622]           = Vec3( -2.9103650e-01,  -1.4337800e-02,    9.5259360e-01);
    positions[623]           = Vec3( -2.8602000e-01,   1.0405050e-01,    8.3648420e-01);
    positions[624]           = Vec3(  6.2908620e-01,  -6.6369160e-01,   -8.8313160e-01);
    positions[625]           = Vec3(  5.3309000e-01,  -6.6824080e-01,   -8.8386380e-01);
    positions[626]           = Vec3(  6.6687380e-01,  -7.2037270e-01,   -8.1674370e-01);
    positions[627]           = Vec3(  2.5101170e-01,  -8.8838680e-01,    2.2900940e-01);
    positions[628]           = Vec3(  2.4302200e-01,  -8.1686710e-01,    1.6969450e-01);
    positions[629]           = Vec3(  3.4457660e-01,  -8.9596990e-01,    2.4839760e-01);
    positions[630]           = Vec3( -9.1418940e-01,   8.0389630e-01,    7.8826000e-01);
    positions[631]           = Vec3( -8.3833600e-01,   7.4209380e-01,    7.8290720e-01);
    positions[632]           = Vec3( -9.9161100e-01,   7.5608500e-01,    8.2971860e-01);
    positions[633]           = Vec3(  7.9708930e-01,  -3.2882190e-01,    7.1789600e-01);
    positions[634]           = Vec3(  8.5609970e-01,  -2.5716920e-01,    7.4938090e-01);
    positions[635]           = Vec3(  7.9853320e-01,  -3.2248890e-01,    6.2155040e-01);
    positions[636]           = Vec3(  7.9743030e-01,  -6.0061740e-01,    7.6822330e-01);
    positions[637]           = Vec3(  8.2105340e-01,  -5.0895770e-01,    7.5902860e-01);
    positions[638]           = Vec3(  7.2970170e-01,  -6.0508550e-01,    8.3860140e-01);
    positions[639]           = Vec3( -1.1738970e-01,  -5.9305270e-01,    7.0381050e-01);
    positions[640]           = Vec3( -1.5290840e-01,  -6.5518590e-01,    6.3431800e-01);
    positions[641]           = Vec3( -1.6038250e-01,  -5.0776740e-01,    6.8496070e-01);
    positions[642]           = Vec3(  5.8567050e-01,   3.6131160e-01,    3.0656670e-01);
    positions[643]           = Vec3(  5.1450330e-01,   4.2381370e-01,    2.6162660e-01);
    positions[644]           = Vec3(  5.3597340e-01,   3.2574600e-01,    3.8272470e-01);
    positions[645]           = Vec3( -7.5114680e-01,   3.5944460e-01,    2.4369600e-01);
    positions[646]           = Vec3( -8.1938720e-01,   4.1907000e-01,    2.7265900e-01);
    positions[647]           = Vec3( -7.8315770e-01,   3.2541650e-01,    1.6165560e-01);

    system.addForce(amoebaVdwForce);

    std::string platformName;
    platformName = "CUDA";
    LangevinIntegrator integrator(0.0, 0.1, 0.01);
    Context context(system, integrator, Platform::getPlatformByName(platformName));

    context.setPositions(positions);
    State state                      = context.getState(State::Forces | State::Energy);
    forces                           = state.getForces();
    energy                           = state.getPotentialEnergy();
}

// test employing box of 216 water molecules w/ and w/o dispersion correction

void testVdwWater(int includeVdwDispersionCorrection) {


    std::string testName;
    if (includeVdwDispersionCorrection) {
        testName      = "testVdwWaterWithDispersionCorrection";
    } else {
        testName      = "testVdwWater";
    }

    int numberOfParticles     = 648;
    double boxDimension       = 1.8643;
    double cutoff             = 0.9;

    std::vector<Vec3> forces;
    double energy;
    setupAndGetForcesEnergyVdwWater("CUBIC-MEAN", "HHG", cutoff, boxDimension, includeVdwDispersionCorrection, forces, energy);
    std::vector<Vec3> expectedForces(numberOfParticles);

    // initialize expected energy and forces

    double expectedEnergy;
    if (includeVdwDispersionCorrection) {
        expectedEnergy     = 4.0108819792e+03;
    } else {
        expectedEnergy     = 4.0349101e+03;
    }

    expectedForces[0]         = Vec3(  2.3025909e+02,  -1.0422757e+01,  -2.1413965e+02);
    expectedForces[1]         = Vec3(  1.1261936e+02,   5.5882575e+02,   9.6539143e+01);
    expectedForces[2]         = Vec3( -8.8436857e+01,  -4.4737313e+01,   2.5242022e+02);
    expectedForces[3]         = Vec3(  6.3548886e+02,  -1.9582636e+02,  -1.2229882e+02);
    expectedForces[4]         = Vec3( -3.9513809e+02,  -1.3738635e+02,  -1.2488717e+02);
    expectedForces[5]         = Vec3(  6.8771170e+00,   9.1574345e+01,  -1.3865672e+00);
    expectedForces[6]         = Vec3(  3.6928792e+02,  -7.1025648e+01,   5.2550320e+02);
    expectedForces[7]         = Vec3(  1.5531325e+02,  -7.9256260e+02,   3.8119809e+02);
    expectedForces[8]         = Vec3(  1.5408049e+02,   1.4633285e+02,  -5.4573735e+02);
    expectedForces[9]         = Vec3( -3.2549641e+02,   2.8613802e+01,   1.8082150e+02);
    expectedForces[10]        = Vec3(  1.9968249e+02,  -7.7742415e+01,  -4.2188467e+02);
    expectedForces[11]        = Vec3(  4.6645952e+01,   2.1362909e+01,   5.7120135e+01);
    expectedForces[12]        = Vec3(  4.2971265e+01,  -4.2927865e+01,  -7.6319121e+02);
    expectedForces[13]        = Vec3(  2.1889551e+02,  -5.1806784e+01,   7.6637073e+01);
    expectedForces[14]        = Vec3( -3.0028991e+02,   9.0479382e+01,   1.3199449e+02);
    expectedForces[15]        = Vec3( -1.8041801e+00,  -1.2176813e+03,  -5.9679371e+02);
    expectedForces[16]        = Vec3( -1.4773339e+02,  -1.2582057e+02,   8.4242040e+02);
    expectedForces[17]        = Vec3(  6.8633312e+02,  -2.4864325e+01,   1.0675519e+01);
    expectedForces[18]        = Vec3(  4.8211733e+02,  -2.2044742e+02,   1.8105309e+02);
    expectedForces[19]        = Vec3( -2.1956736e+02,   5.4786741e+02,  -1.4827263e+02);
    expectedForces[20]        = Vec3( -1.0476100e+02,  -1.6922280e+02,   8.3220473e+01);
    expectedForces[21]        = Vec3(  5.5455350e+01,   1.7851224e+02,   6.0598643e+02);
    expectedForces[22]        = Vec3(  7.1915871e+00,  -1.6455057e+01,  -6.3929797e+00);
    expectedForces[23]        = Vec3( -2.4371864e+02,  -7.0765288e+02,  -1.7820148e+02);
    expectedForces[24]        = Vec3( -4.3715367e+02,   2.8707226e+02,   2.2558361e+02);
    expectedForces[25]        = Vec3(  2.3735557e+00,  -4.2064130e+00,   2.5933632e+01);
    expectedForces[26]        = Vec3(  1.0192415e+02,   5.3205842e+01,  -2.0158900e+02);
    expectedForces[27]        = Vec3( -5.1276812e+02,   2.4860179e+02,   2.5150125e+02);
    expectedForces[28]        = Vec3(  3.9608987e+02,  -4.0544816e+01,  -2.1579831e+02);
    expectedForces[29]        = Vec3( -2.5926005e+02,  -1.1691294e+02,  -4.7278417e+02);
    expectedForces[30]        = Vec3( -1.1036158e+02,   3.3229287e+01,   8.8351491e+01);
    expectedForces[31]        = Vec3( -4.1178603e+00,   1.2957632e+00,   1.2617868e+00);
    expectedForces[32]        = Vec3(  1.9401524e+02,   1.0030314e+02,   3.0386141e+01);
    expectedForces[33]        = Vec3(  1.0870981e+02,   2.1865882e+02,   6.7672725e+02);
    expectedForces[34]        = Vec3(  3.6812338e+01,  -2.9347382e-01,   1.2959790e+01);
    expectedForces[35]        = Vec3( -5.5738560e+01,   2.2268442e+02,  -1.4911862e+02);
    expectedForces[36]        = Vec3( -1.6659719e+02,   8.1941190e+01,   7.8132996e+01);
    expectedForces[37]        = Vec3( -1.1985659e+02,   7.1891763e+01,  -3.2559738e+02);
    expectedForces[38]        = Vec3(  1.4441653e+02,   4.1824194e+02,   3.4293255e+01);
    expectedForces[39]        = Vec3(  1.1307312e+02,  -3.4981802e+02,  -1.0676819e+02);
    expectedForces[40]        = Vec3( -1.2857487e+02,  -4.0798950e+01,  -1.7800009e+02);
    expectedForces[41]        = Vec3( -4.1244991e+02,   2.3439565e+02,   3.8546709e+02);
    expectedForces[42]        = Vec3(  2.1799653e+02,   2.9397348e+02,  -3.8955146e+02);
    expectedForces[43]        = Vec3( -9.4254462e+01,   5.0157547e+01,   1.3707146e+02);
    expectedForces[44]        = Vec3(  8.0520929e+02,  -1.0228045e+03,  -4.5813594e+01);
    expectedForces[45]        = Vec3(  5.3451125e+02,  -7.0725115e+02,  -2.3828883e+02);
    expectedForces[46]        = Vec3(  3.6991801e+02,  -1.0410466e+03,   8.4889413e+02);
    expectedForces[47]        = Vec3(  5.8096458e+01,   2.3859368e+02,   7.3104979e+01);
    expectedForces[48]        = Vec3( -1.3230378e+03,   7.1827816e+02,   9.4372946e+02);
    expectedForces[49]        = Vec3( -4.2398413e+02,  -3.6766481e+02,  -1.6666711e+02);
    expectedForces[50]        = Vec3(  8.6374766e+00,  -2.7870227e+02,   3.0621233e+02);
    expectedForces[51]        = Vec3( -2.4404981e+02,   1.0510491e+03,  -3.5202340e+02);
    expectedForces[52]        = Vec3(  1.7461155e+02,   1.3419948e+02,   5.0279733e+02);
    expectedForces[53]        = Vec3(  9.7089193e+01,  -1.4826120e+02,  -2.3426711e+02);
    expectedForces[54]        = Vec3( -1.0647432e+02,   2.3266543e+02,  -3.5025025e+02);
    expectedForces[55]        = Vec3( -1.9951717e+02,  -1.0804258e+02,   6.2114265e+01);
    expectedForces[56]        = Vec3(  2.5354284e+02,  -3.5074926e+02,  -2.9176946e+01);
    expectedForces[57]        = Vec3( -8.8328002e+02,  -1.8937053e+02,   2.0772189e+02);
    expectedForces[58]        = Vec3(  4.5809698e+01,  -1.8636513e+00,   1.6742910e+01);
    expectedForces[59]        = Vec3(  8.6670673e+01,   1.2653519e+02,   2.0046372e+02);
    expectedForces[60]        = Vec3( -9.6142275e+02,   1.1800289e+03,   1.0064471e+02);
    expectedForces[61]        = Vec3( -6.9751453e+01,  -1.9324691e+02,   2.2638891e+01);
    expectedForces[62]        = Vec3(  1.3671815e+02,   1.0129231e+01,  -6.1136892e+01);
    expectedForces[63]        = Vec3(  7.6032938e+02,  -1.9779542e+02,  -2.5706161e+01);
    expectedForces[64]        = Vec3( -9.8966607e+00,   3.6294058e+02,  -2.3247376e+02);
    expectedForces[65]        = Vec3( -1.3142673e+02,   9.9368000e+00,   1.4225069e+02);
    expectedForces[66]        = Vec3(  5.2691625e+01,  -4.0618331e+02,  -1.0942115e+02);
    expectedForces[67]        = Vec3(  5.8137742e+01,   2.1942441e+02,  -1.7191686e+02);
    expectedForces[68]        = Vec3(  1.2480778e+02,   2.4086799e+01,  -1.9588545e+02);
    expectedForces[69]        = Vec3(  2.4577512e+02,  -2.3516051e+01,   2.3764379e+02);
    expectedForces[70]        = Vec3(  4.2298797e+01,   8.6807583e+01,   2.9276129e+02);
    expectedForces[71]        = Vec3(  3.0309915e+01,  -3.4686002e+02,   1.0849989e+02);
    expectedForces[72]        = Vec3(  5.5170999e+02,  -1.5797913e+02,   9.8379834e+01);
    expectedForces[73]        = Vec3( -5.3706921e+02,   2.2236470e+02,   2.7548665e+02);
    expectedForces[74]        = Vec3( -2.7734936e+02,  -4.2949693e+02,  -8.2303469e+02);
    expectedForces[75]        = Vec3(  1.1749001e+03,   9.2734627e+02,  -2.6342291e+02);
    expectedForces[76]        = Vec3( -1.2345216e+02,  -6.6220297e+01,   1.1022778e+02);
    expectedForces[77]        = Vec3( -4.4861284e+01,  -7.4560959e+01,  -8.0791275e+01);
    expectedForces[78]        = Vec3(  1.7692476e+01,  -9.9067626e+02,  -3.3083363e+02);
    expectedForces[79]        = Vec3(  6.4276566e+01,   5.1865608e+02,   4.1241764e+02);
    expectedForces[80]        = Vec3(  6.5867999e+02,   2.0995465e+02,  -5.0989840e+02);
    expectedForces[81]        = Vec3( -5.2826435e+02,  -1.3719887e+02,   1.9730314e+02);
    expectedForces[82]        = Vec3(  2.3329314e+02,   5.1589944e+01,  -3.3942818e+02);
    expectedForces[83]        = Vec3( -2.5756193e+00,   1.3788257e+02,   1.3074778e+02);
    expectedForces[84]        = Vec3( -2.8925811e+01,   8.8185346e+01,   1.1547804e+02);
    expectedForces[85]        = Vec3( -3.3582542e+01,  -1.3079843e+02,   3.4774799e+02);
    expectedForces[86]        = Vec3( -5.9585326e+01,  -3.7973755e+01,   1.0206767e+01);
    expectedForces[87]        = Vec3( -7.1691223e+01,  -1.9547311e+02,   9.3856472e+02);
    expectedForces[88]        = Vec3(  1.4832977e+02,   4.2005795e+02,   2.7235044e+01);
    expectedForces[89]        = Vec3(  8.6221698e+00,  -1.3935906e+01,  -2.5473157e+00);
    expectedForces[90]        = Vec3(  3.0044348e+02,   3.4320406e+02,   2.5305998e+02);
    expectedForces[91]        = Vec3( -1.6694013e+02,   7.4645776e+02,   3.3266168e+01);
    expectedForces[92]        = Vec3(  1.6931562e+02,  -2.3972286e+01,  -3.0802865e+02);
    expectedForces[93]        = Vec3( -2.5931152e+02,  -9.7825375e+01,  -5.5598386e+02);
    expectedForces[94]        = Vec3( -7.8473519e+01,  -2.6244561e+02,   2.9622027e+02);
    expectedForces[95]        = Vec3(  1.5788473e+02,   2.9806605e+02,   1.8283700e+02);
    expectedForces[96]        = Vec3( -4.2033208e+02,   4.8745845e+02,   9.7957219e+01);
    expectedForces[97]        = Vec3(  1.3261733e+00,   3.7143158e+00,  -1.8329960e+00);
    expectedForces[98]        = Vec3(  7.4741051e+02,   1.0482696e+03,   1.7182445e+02);
    expectedForces[99]        = Vec3(  2.0366606e+02,   1.1586652e+02,   1.5482927e+03);
    expectedForces[100]       = Vec3(  5.7569884e+02,  -2.2068597e+02,  -1.9305186e+02);
    expectedForces[101]       = Vec3( -2.4259775e+02,  -6.1359858e+01,  -7.5753918e+01);
    expectedForces[102]       = Vec3(  2.4382837e+02,   1.3660693e+01,  -2.8717824e+02);
    expectedForces[103]       = Vec3( -4.1914381e+02,  -1.4506334e+02,   7.7692658e+01);
    expectedForces[104]       = Vec3(  9.5422158e+01,  -1.4261849e+01,   1.0122453e+02);
    expectedForces[105]       = Vec3( -6.1222151e+02,   6.1398147e+01,   3.2961693e+02);
    expectedForces[106]       = Vec3( -2.7982116e+02,  -3.7564793e+02,  -2.4318658e+02);
    expectedForces[107]       = Vec3(  8.2892383e+01,   1.1992522e+01,   2.4802899e+01);
    expectedForces[108]       = Vec3( -3.3626753e+02,  -7.8892295e+01,   8.1184535e+02);
    expectedForces[109]       = Vec3( -4.6396103e+01,   3.1165111e+01,  -3.5184949e+01);
    expectedForces[110]       = Vec3(  4.6295297e+02,  -8.0135823e+01,  -3.2963637e+02);
    expectedForces[111]       = Vec3(  5.5325945e+02,   3.8213691e+02,  -5.9449911e+02);
    expectedForces[112]       = Vec3( -6.1770559e+02,   2.2322880e+02,  -2.1022934e+01);
    expectedForces[113]       = Vec3(  2.3716825e+02,  -3.8085281e+02,   4.9618193e+02);
    expectedForces[114]       = Vec3(  2.9972570e+02,   3.3773759e+02,   3.4933501e+02);
    expectedForces[115]       = Vec3(  1.9913039e+02,  -3.2035626e+02,  -9.7987777e+01);
    expectedForces[116]       = Vec3( -2.3506456e+02,   4.0030561e+02,  -6.2604809e+02);
    expectedForces[117]       = Vec3(  6.0408777e+02,   5.6271523e+02,  -5.8742635e+02);
    expectedForces[118]       = Vec3( -5.0707285e+02,   7.8511351e+02,  -1.7561830e+02);
    expectedForces[119]       = Vec3(  5.6759945e+00,   6.9735233e+00,  -9.8853787e+00);
    expectedForces[120]       = Vec3(  1.1918531e+02,  -5.6200785e+02,   3.2333502e+02);
    expectedForces[121]       = Vec3(  5.1801315e+02,  -3.4460898e+02,   1.1488242e+02);
    expectedForces[122]       = Vec3(  7.5859703e+00,   6.9528951e+01,  -6.2437444e+00);
    expectedForces[123]       = Vec3(  1.1990899e+03,  -1.2507886e+01,   1.1906562e+03);
    expectedForces[124]       = Vec3(  2.6001536e+02,  -9.5041373e+01,  -4.2547252e+02);
    expectedForces[125]       = Vec3( -3.5601586e+01,   9.1460220e+02,   6.2306102e+02);
    expectedForces[126]       = Vec3( -3.8254176e+01,   2.4796243e+02,   3.0638741e+02);
    expectedForces[127]       = Vec3( -3.8798311e+02,   1.2532794e+02,  -3.1926601e+01);
    expectedForces[128]       = Vec3(  5.7672682e+01,   1.7487900e+02,  -1.0319581e+02);
    expectedForces[129]       = Vec3(  4.4276215e+02,  -1.5517993e+02,   3.8122873e+02);
    expectedForces[130]       = Vec3( -2.8119980e+02,  -2.6530936e+02,   4.3545814e+01);
    expectedForces[131]       = Vec3(  5.9375768e+00,   1.6335780e+01,  -3.6204644e+02);
    expectedForces[132]       = Vec3(  5.9733834e+02,  -2.9843602e+02,   1.0541092e+03);
    expectedForces[133]       = Vec3( -3.4568461e+02,   2.2794808e+02,  -8.4675886e+01);
    expectedForces[134]       = Vec3( -6.7882706e+01,  -4.5016799e+02,  -1.9387646e+02);
    expectedForces[135]       = Vec3(  5.6125478e+02,  -1.1741092e+03,   1.4214256e+02);
    expectedForces[136]       = Vec3(  1.2766648e+02,   1.2519441e+02,  -4.4359129e+02);
    expectedForces[137]       = Vec3( -2.7235507e+01,  -1.0871423e+01,   3.2406976e+01);
    expectedForces[138]       = Vec3( -1.2110689e+03,   1.3038667e+02,  -7.5331947e+02);
    expectedForces[139]       = Vec3(  3.2504988e+02,  -7.1577942e+02,  -2.4351454e+02);
    expectedForces[140]       = Vec3(  2.6466925e+02,   3.6878498e+02,  -6.2044711e+01);
    expectedForces[141]       = Vec3( -6.2647276e+02,  -4.7112562e+02,  -6.1885184e+01);
    expectedForces[142]       = Vec3(  3.0996699e+01,   3.1033732e+02,  -5.7493337e+02);
    expectedForces[143]       = Vec3(  2.4664962e+01,   7.8306665e+02,   3.7907751e+02);
    expectedForces[144]       = Vec3( -1.9847551e+02,  -3.0052499e+02,  -6.5192951e+01);
    expectedForces[145]       = Vec3( -5.0686773e+02,   1.1972933e+02,   3.6202910e+02);
    expectedForces[146]       = Vec3(  4.0004521e+02,  -4.3205645e+02,   2.8583074e+02);
    expectedForces[147]       = Vec3( -2.9447542e+01,   6.3434099e+02,  -6.5189179e+02);
    expectedForces[148]       = Vec3(  4.1588516e+02,  -1.5534873e+02,  -1.1031323e+01);
    expectedForces[149]       = Vec3(  1.2761543e+01,  -5.6871456e+01,   2.8032565e+02);
    expectedForces[150]       = Vec3( -4.7527199e+01,   6.2549647e+02,  -8.5007589e+02);
    expectedForces[151]       = Vec3(  1.6598949e+01,   1.1210199e+02,   2.0667266e+02);
    expectedForces[152]       = Vec3( -2.4012423e+00,  -6.1728962e+01,  -1.5057977e+01);
    expectedForces[153]       = Vec3( -3.6934056e+02,   2.3087063e+02,  -3.3268511e+02);
    expectedForces[154]       = Vec3(  1.2658460e+02,   4.1809355e+02,   1.9550635e+02);
    expectedForces[155]       = Vec3( -1.9528152e+02,  -1.9035919e+02,   1.4986704e+02);
    expectedForces[156]       = Vec3( -5.9608803e+00,   2.7258190e+02,   5.5318395e+02);
    expectedForces[157]       = Vec3( -1.0545084e+01,  -1.8397001e+00,  -3.3615498e+00);
    expectedForces[158]       = Vec3( -4.0827462e+02,  -6.4070993e+01,  -2.6614992e+02);
    expectedForces[159]       = Vec3( -1.5773827e+02,   4.7708702e+02,  -9.9622711e+02);
    expectedForces[160]       = Vec3(  8.3896122e+01,   2.5959601e+02,   1.7728081e+02);
    expectedForces[161]       = Vec3( -2.4859055e+02,   8.1617294e+01,   2.7905586e+02);
    expectedForces[162]       = Vec3(  1.3694593e+02,  -6.2640483e+02,  -2.5308878e+02);
    expectedForces[163]       = Vec3( -1.5156337e+02,   1.1567379e+02,   1.8384602e+02);
    expectedForces[164]       = Vec3( -1.4794785e+02,  -8.0381965e+01,   9.4028021e+02);
    expectedForces[165]       = Vec3( -3.5851252e+02,   2.6066489e+02,  -4.1225343e+01);
    expectedForces[166]       = Vec3(  1.0820389e+02,  -8.3485138e+01,   9.5047812e+01);
    expectedForces[167]       = Vec3(  2.0019411e+02,   5.0864719e+02,  -4.8662239e+01);
    expectedForces[168]       = Vec3( -3.0716707e+02,  -9.8820762e+02,   1.3660536e+01);
    expectedForces[169]       = Vec3(  8.5788590e+02,  -1.0335921e+02,  -2.2197815e+02);
    expectedForces[170]       = Vec3(  9.8072705e+00,   6.6224870e+02,  -1.1594055e+02);
    expectedForces[171]       = Vec3( -3.6572453e+01,   2.4062731e+02,   1.7831682e+02);
    expectedForces[172]       = Vec3( -2.1274187e+00,  -2.5471785e+02,  -2.5813525e+02);
    expectedForces[173]       = Vec3(  3.4763673e+02,  -2.0753062e+02,   7.5125390e+02);
    expectedForces[174]       = Vec3(  2.9010685e+01,  -4.1894498e+02,  -6.3847846e+02);
    expectedForces[175]       = Vec3(  3.0255245e+02,   2.5945698e+02,   2.3530662e+01);
    expectedForces[176]       = Vec3( -2.0827676e+01,   2.2568888e+01,   1.0890440e+01);
    expectedForces[177]       = Vec3(  2.2415739e+02,   3.5720682e+02,  -8.7208211e+01);
    expectedForces[178]       = Vec3( -2.9808276e+02,   7.0108444e+01,  -1.5361211e+02);
    expectedForces[179]       = Vec3(  2.3882234e+02,  -2.6764176e+02,  -1.7847229e+01);
    expectedForces[180]       = Vec3( -4.7380966e+02,   5.2547306e+01,   7.4446526e+02);
    expectedForces[181]       = Vec3(  4.2804072e+00,  -2.2456284e+02,  -2.0288206e+02);
    expectedForces[182]       = Vec3(  2.5474630e+00,   2.5211778e+01,   1.2675317e+01);
    expectedForces[183]       = Vec3( -1.4519676e+02,  -1.7268277e+02,  -8.4750014e+02);
    expectedForces[184]       = Vec3( -1.6850633e+02,   9.5223081e+02,  -2.5415253e+02);
    expectedForces[185]       = Vec3( -2.1600185e+02,  -1.6313340e+02,   1.8746411e+02);
    expectedForces[186]       = Vec3(  4.5719348e+01,  -5.8877840e+01,  -1.1199622e+02);
    expectedForces[187]       = Vec3( -9.0377199e+01,   2.2045926e+02,  -3.2532356e+02);
    expectedForces[188]       = Vec3( -1.4523032e+02,  -2.2127499e+02,   1.1526636e+02);
    expectedForces[189]       = Vec3(  1.9857975e+02,  -1.8170037e+02,   3.2016150e+02);
    expectedForces[190]       = Vec3(  1.5659789e+02,   5.2984420e+01,  -1.0121553e+02);
    expectedForces[191]       = Vec3(  1.3925504e+01,  -1.1627739e+02,   1.7106106e+00);
    expectedForces[192]       = Vec3(  1.7560737e+01,   3.1275548e+02,   1.4568181e+02);
    expectedForces[193]       = Vec3(  8.5129168e-01,   3.2688468e+01,  -2.0471785e+02);
    expectedForces[194]       = Vec3(  2.5640904e+02,  -1.8345358e+02,   3.3082812e+02);
    expectedForces[195]       = Vec3(  1.5159266e+02,  -7.9091822e+02,   7.6475475e+02);
    expectedForces[196]       = Vec3(  7.8546101e+02,   1.7887845e+02,  -7.4508377e+01);
    expectedForces[197]       = Vec3(  6.9010029e+00,   2.1336674e+01,  -7.8256484e+01);
    expectedForces[198]       = Vec3( -3.4033432e+01,  -1.6471994e+02,  -2.8721884e+02);
    expectedForces[199]       = Vec3( -1.2248359e+02,   5.3002993e+01,   3.5296174e+02);
    expectedForces[200]       = Vec3(  1.7497881e+02,   1.1744848e+02,   1.8127502e+02);
    expectedForces[201]       = Vec3( -1.1710882e+03,   8.7114752e+02,  -1.3353036e+02);
    expectedForces[202]       = Vec3(  2.7510368e+02,  -1.9393814e+02,   2.6671453e+02);
    expectedForces[203]       = Vec3(  7.3421287e+01,  -3.8409917e+01,  -1.0910368e+02);
    expectedForces[204]       = Vec3( -8.7394869e+02,  -1.2323430e+03,  -1.3693068e+02);
    expectedForces[205]       = Vec3(  2.5439149e+02,   2.4002077e+02,   3.0231859e+02);
    expectedForces[206]       = Vec3( -9.8167431e+01,   3.1515278e+02,  -2.5512038e+02);
    expectedForces[207]       = Vec3(  2.7501189e+02,  -4.9189774e+02,   4.5827877e+02);
    expectedForces[208]       = Vec3( -8.2333908e+02,  -7.0609881e+02,   2.9188070e+02);
    expectedForces[209]       = Vec3( -2.9200330e+02,   1.9955124e+02,  -7.9227826e+00);
    expectedForces[210]       = Vec3(  8.0151971e+02,   3.8996996e+02,  -1.9508756e+02);
    expectedForces[211]       = Vec3( -6.2168096e+01,   2.4323951e+01,   9.8845638e+01);
    expectedForces[212]       = Vec3(  2.2980428e+01,  -6.2829693e+02,  -7.5253154e+02);
    expectedForces[213]       = Vec3( -5.7502634e+02,   8.2310635e+02,   3.7934867e+02);
    expectedForces[214]       = Vec3( -8.9168889e+01,   5.0052167e+01,  -8.6850564e+01);
    expectedForces[215]       = Vec3(  2.7414821e+02,  -5.9725191e+01,  -3.7278778e+00);
    expectedForces[216]       = Vec3( -5.3230933e+02,   4.7698725e+02,  -2.5857751e+02);
    expectedForces[217]       = Vec3(  2.1615644e+02,   6.0001233e+01,  -1.1880182e+02);
    expectedForces[218]       = Vec3(  2.9837823e+01,  -3.9275454e+01,  -4.3913571e+01);
    expectedForces[219]       = Vec3( -5.6516561e+02,  -1.6539113e+02,   3.3185872e+02);
    expectedForces[220]       = Vec3( -1.5484833e+01,  -1.0085566e+01,  -2.8199778e+01);
    expectedForces[221]       = Vec3(  6.4295337e+02,  -2.7398997e+02,   8.2320867e+02);
    expectedForces[222]       = Vec3(  3.3816540e+02,  -1.8492046e+02,   8.6215770e+02);
    expectedForces[223]       = Vec3(  2.4367872e+02,   8.9527142e+01,  -5.9973680e+02);
    expectedForces[224]       = Vec3( -6.1601439e+02,  -5.3093862e+02,   4.1696221e+01);
    expectedForces[225]       = Vec3(  2.2521478e+02,  -2.9234267e+02,   2.9760244e+02);
    expectedForces[226]       = Vec3( -1.2672295e+02,   2.5753736e+01,  -1.6216252e+02);
    expectedForces[227]       = Vec3(  1.1136724e+02,   4.3446676e+00,  -1.4045447e+01);
    expectedForces[228]       = Vec3(  3.5554849e+02,   7.7446255e+00,   5.8722356e+02);
    expectedForces[229]       = Vec3( -7.8313952e+02,  -5.9206053e+01,   5.1181391e+02);
    expectedForces[230]       = Vec3( -9.5647729e+01,  -4.6045792e+02,  -3.1515210e+02);
    expectedForces[231]       = Vec3(  8.7273284e+01,   8.0310279e+02,  -1.2212040e+03);
    expectedForces[232]       = Vec3( -9.0853819e+02,  -7.6005570e+02,   1.8578569e+02);
    expectedForces[233]       = Vec3(  1.2061615e+02,  -3.2834560e+02,   5.1291192e+01);
    expectedForces[234]       = Vec3(  6.0149262e+02,   1.6081852e+03,  -1.0922015e+03);
    expectedForces[235]       = Vec3( -5.8087260e+01,  -2.6676779e+02,  -1.5006920e+02);
    expectedForces[236]       = Vec3( -1.0997771e+01,   1.3012815e+02,   5.8615633e+02);
    expectedForces[237]       = Vec3(  4.9821503e+02,  -9.5783644e+01,  -6.5708073e+02);
    expectedForces[238]       = Vec3( -4.9394992e+02,  -7.0445435e+02,  -1.5540228e+02);
    expectedForces[239]       = Vec3(  1.5710340e+01,   1.9170336e+02,   4.6707393e+02);
    expectedForces[240]       = Vec3( -1.2599557e+03,   5.6104975e+02,   6.5853210e+02);
    expectedForces[241]       = Vec3(  1.3032986e+01,  -2.7481865e+01,   8.0373029e+00);
    expectedForces[242]       = Vec3(  1.4371505e+02,  -5.2615293e+01,  -3.7929456e+02);
    expectedForces[243]       = Vec3(  6.0025739e+01,   1.4474905e+02,   3.6414330e+02);
    expectedForces[244]       = Vec3(  4.5364033e+02,   7.0785823e+01,  -2.9279112e+02);
    expectedForces[245]       = Vec3( -9.7012956e+01,  -8.2545229e+01,  -5.4567106e+01);
    expectedForces[246]       = Vec3( -3.6449281e+02,  -5.6669832e+02,  -3.5196038e+02);
    expectedForces[247]       = Vec3(  2.6571354e+02,  -2.5051716e+02,  -4.7199161e+01);
    expectedForces[248]       = Vec3( -4.7285701e+01,   1.7950093e+02,   8.0633971e+01);
    expectedForces[249]       = Vec3(  1.8035439e+02,   4.6977540e+02,   2.2589769e+02);
    expectedForces[250]       = Vec3(  5.9853399e+01,   1.5477006e+02,  -1.7011806e+02);
    expectedForces[251]       = Vec3( -2.6084741e+02,  -2.0538756e+02,   1.0269002e+02);
    expectedForces[252]       = Vec3(  4.0221722e+02,  -1.0368371e+02,  -1.3842460e+01);
    expectedForces[253]       = Vec3( -4.4264447e+02,  -3.0334810e+02,   1.6740779e+02);
    expectedForces[254]       = Vec3( -1.7373662e+02,   1.5020202e+02,  -3.3761418e+02);
    expectedForces[255]       = Vec3(  4.7508953e+01,  -5.4028784e+02,  -3.6291329e+01);
    expectedForces[256]       = Vec3(  5.4807168e+02,  -2.0293664e+02,   5.5114363e+02);
    expectedForces[257]       = Vec3( -4.9787370e+01,   4.5106329e+01,  -5.7524550e+01);
    expectedForces[258]       = Vec3( -5.1584935e+02,   1.5014186e+02,  -6.6369834e+02);
    expectedForces[259]       = Vec3( -6.0140298e+01,  -9.4127280e+01,   8.7495796e+01);
    expectedForces[260]       = Vec3(  8.0619588e+02,  -2.1066454e+02,  -5.9371508e+01);
    expectedForces[261]       = Vec3(  1.2817429e+02,  -2.7044116e+02,   5.3213275e+02);
    expectedForces[262]       = Vec3( -4.1467442e-01,  -1.3641138e+02,  -8.5576287e+01);
    expectedForces[263]       = Vec3( -5.6508528e+01,   2.4091504e+01,  -7.9236147e+01);
    expectedForces[264]       = Vec3(  1.4123799e+02,   2.1326444e+02,   1.9141445e+02);
    expectedForces[265]       = Vec3(  1.6579721e+02,  -3.2221503e+02,  -1.5197776e+02);
    expectedForces[266]       = Vec3( -2.5871222e+02,  -2.6345996e+02,  -1.1149413e+03);
    expectedForces[267]       = Vec3( -1.2859743e+02,  -3.9214753e+02,  -6.5874862e+02);
    expectedForces[268]       = Vec3(  1.7237157e+02,  -2.8729010e+02,   3.2000889e+02);
    expectedForces[269]       = Vec3(  6.6312203e+02,   3.5683926e+02,   9.7235616e+01);
    expectedForces[270]       = Vec3( -1.1860349e+02,  -4.3528895e+02,   1.1598672e+03);
    expectedForces[271]       = Vec3(  7.7135123e+01,  -2.1671305e+02,  -2.3088548e+02);
    expectedForces[272]       = Vec3(  7.2863708e+01,   5.2065158e+02,   1.5704752e+02);
    expectedForces[273]       = Vec3( -3.6727477e+02,  -6.0885463e+01,   4.1193960e+01);
    expectedForces[274]       = Vec3(  5.4858881e+01,   1.6886577e+02,  -1.1662392e+01);
    expectedForces[275]       = Vec3(  3.2070061e+01,  -1.6138709e+01,  -7.6648912e+00);
    expectedForces[276]       = Vec3(  2.5695892e+02,   7.7719531e+01,  -8.3154765e+02);
    expectedForces[277]       = Vec3( -3.3460202e+01,  -1.3077348e+02,   3.6110539e+01);
    expectedForces[278]       = Vec3( -4.1360402e+02,   2.7754335e+02,  -1.1139467e+02);
    expectedForces[279]       = Vec3(  1.0012775e+02,   1.4365237e+03,   2.6972544e+02);
    expectedForces[280]       = Vec3(  7.8750024e+02,  -2.3376123e+02,   9.5969232e+02);
    expectedForces[281]       = Vec3( -2.2493912e+02,  -1.1444344e+01,  -6.3713547e+01);
    expectedForces[282]       = Vec3( -9.4211966e+01,  -3.9286586e+02,  -7.0932607e+01);
    expectedForces[283]       = Vec3(  1.6192138e+02,   2.3010404e+02,  -2.1511498e+02);
    expectedForces[284]       = Vec3( -4.9644702e+02,  -1.1894844e+01,  -4.7458477e+02);
    expectedForces[285]       = Vec3( -1.0151600e+02,  -4.2789293e+02,   5.2698223e+02);
    expectedForces[286]       = Vec3(  1.1350409e+01,   1.5126445e+02,  -5.2740587e+01);
    expectedForces[287]       = Vec3( -4.7715521e+02,  -3.6476152e+02,  -5.4346958e+02);
    expectedForces[288]       = Vec3( -9.4019721e+01,  -2.4758228e+02,   4.2219712e+02);
    expectedForces[289]       = Vec3(  3.1525053e+01,   3.7927245e+02,  -2.4577041e+02);
    expectedForces[290]       = Vec3(  4.4004025e+02,  -9.4883983e+01,  -6.9075632e+01);
    expectedForces[291]       = Vec3( -2.6553328e+02,  -2.7059185e+02,   3.1145935e+02);
    expectedForces[292]       = Vec3(  1.1199463e+02,   9.0805870e+01,   8.0258804e+01);
    expectedForces[293]       = Vec3(  7.8056273e+01,  -3.5931761e+01,  -1.9249727e+02);
    expectedForces[294]       = Vec3( -1.2443366e+01,  -1.4631721e+02,   7.9733561e+01);
    expectedForces[295]       = Vec3( -3.2866900e+02,   1.1527926e+02,  -9.3135590e+01);
    expectedForces[296]       = Vec3(  7.5078643e+01,  -7.2653837e+02,  -1.0830908e+02);
    expectedForces[297]       = Vec3(  1.0361037e+02,   1.5478274e+02,   1.1293089e+03);
    expectedForces[298]       = Vec3( -2.4014663e+02,  -5.6771824e+02,   1.5441974e+02);
    expectedForces[299]       = Vec3( -6.8441034e+01,  -1.1226846e+02,  -5.4497343e+02);
    expectedForces[300]       = Vec3(  2.9803689e+02,   1.0649462e+03,   2.5408768e+02);
    expectedForces[301]       = Vec3(  5.6356599e+01,  -2.3035579e+02,   7.6928758e+01);
    expectedForces[302]       = Vec3(  1.0229494e+02,  -1.9615696e+02,  -3.2600199e+02);
    expectedForces[303]       = Vec3( -7.0312065e+01,  -3.1170055e+02,   5.6453888e+02);
    expectedForces[304]       = Vec3( -7.0752582e+01,   1.2273280e+02,  -1.6740503e+01);
    expectedForces[305]       = Vec3(  6.8972998e+02,  -2.8067032e+02,   6.6502848e+01);
    expectedForces[306]       = Vec3( -8.8407472e+02,   3.7818347e+02,  -1.2141945e+03);
    expectedForces[307]       = Vec3( -1.0445122e+02,  -4.4818836e+02,  -9.5227168e+00);
    expectedForces[308]       = Vec3(  8.3015802e+02,   5.0941309e+02,   3.2593793e+01);
    expectedForces[309]       = Vec3(  7.6887033e+02,   4.3490548e+02,  -1.0826101e+03);
    expectedForces[310]       = Vec3(  1.4352949e+02,  -3.9479802e+01,   1.7718790e+02);
    expectedForces[311]       = Vec3( -2.5204642e+02,  -3.3888239e+02,  -2.7365802e+02);
    expectedForces[312]       = Vec3( -4.4783121e+02,   4.4465580e+02,   9.2506336e+02);
    expectedForces[313]       = Vec3( -3.3390337e+01,  -5.0097602e+02,  -3.6812516e+02);
    expectedForces[314]       = Vec3( -4.1862341e+02,   1.9519764e+02,  -2.4571211e+02);
    expectedForces[315]       = Vec3( -1.4176845e+02,  -3.3767798e+01,   1.7800248e+02);
    expectedForces[316]       = Vec3( -4.7411046e+01,   2.3841278e+01,  -1.5183689e+01);
    expectedForces[317]       = Vec3(  1.6723204e+01,   1.4145797e+02,  -1.6930384e+01);
    expectedForces[318]       = Vec3( -8.4841167e+02,  -2.1282818e+02,   8.1614626e+01);
    expectedForces[319]       = Vec3(  6.5941519e+00,  -1.1472179e+02,  -2.3181002e+02);
    expectedForces[320]       = Vec3(  2.2464095e+01,  -5.7606215e+00,   5.2865780e+00);
    expectedForces[321]       = Vec3(  3.8063827e+02,   2.4597014e+02,  -4.1947008e+01);
    expectedForces[322]       = Vec3( -3.1531672e+02,   3.5038104e+02,  -3.2724191e+02);
    expectedForces[323]       = Vec3(  7.4741163e+01,   7.2296536e+01,   3.1768919e+02);
    expectedForces[324]       = Vec3(  7.9578966e+02,  -1.5345175e+02,  -7.1301288e+02);
    expectedForces[325]       = Vec3( -3.1387318e+01,   1.7528515e+01,  -3.4100670e+01);
    expectedForces[326]       = Vec3(  9.3853572e-01,   6.8496077e+00,   1.2234894e+01);
    expectedForces[327]       = Vec3( -1.3099690e+02,  -4.0135958e+02,  -4.3248135e+02);
    expectedForces[328]       = Vec3(  4.3988113e+02,  -1.5944416e+02,   5.4620163e+01);
    expectedForces[329]       = Vec3(  1.3509243e+01,   3.6990268e+00,  -2.6251838e-01);
    expectedForces[330]       = Vec3( -8.7459998e+01,  -6.4319025e+02,   2.2863360e+02);
    expectedForces[331]       = Vec3( -3.5181606e+02,   3.1131978e+02,   2.8037788e+02);
    expectedForces[332]       = Vec3(  1.9545884e+02,   4.3994512e+01,  -2.5158598e+02);
    expectedForces[333]       = Vec3(  1.8100917e+00,   1.9534738e+02,   1.7181161e+02);
    expectedForces[334]       = Vec3( -1.1846876e+02,   1.5368782e+02,   2.2261016e+02);
    expectedForces[335]       = Vec3(  1.4866155e+02,  -4.7983795e+02,   3.0361323e+02);
    expectedForces[336]       = Vec3( -7.7242468e+02,   2.0736079e+02,  -8.6661236e+02);
    expectedForces[337]       = Vec3( -7.9875679e+00,  -9.0946585e+02,  -7.1617627e+01);
    expectedForces[338]       = Vec3(  2.2543457e+02,   2.0215639e+02,   4.4285374e+02);
    expectedForces[339]       = Vec3( -2.6252013e+02,  -1.0336519e+02,  -2.6630451e+02);
    expectedForces[340]       = Vec3(  5.0907589e+02,  -1.2491931e+02,  -4.3015857e+01);
    expectedForces[341]       = Vec3( -6.0362692e+01,   1.0025664e+02,  -9.3389235e+01);
    expectedForces[342]       = Vec3(  5.8779874e+02,   3.7600030e+02,   2.5345454e+02);
    expectedForces[343]       = Vec3( -1.3661048e+02,   1.6498482e+02,  -2.6789925e+01);
    expectedForces[344]       = Vec3(  4.9840419e+01,  -1.2136315e+02,   2.5455320e+01);
    expectedForces[345]       = Vec3(  2.3324807e+02,  -3.6601512e+02,   5.6798114e+01);
    expectedForces[346]       = Vec3(  5.5853546e+01,  -4.9035168e+00,  -5.8726331e+02);
    expectedForces[347]       = Vec3( -3.8918158e+01,   3.0380495e+02,   1.0944883e+01);
    expectedForces[348]       = Vec3( -2.6092065e+02,   4.3443036e+02,  -1.3673862e+02);
    expectedForces[349]       = Vec3(  1.5195166e+02,   1.1879426e+01,  -1.1429493e+02);
    expectedForces[350]       = Vec3( -8.5687107e+02,  -1.7938881e+02,   4.5576455e+00);
    expectedForces[351]       = Vec3( -1.0117388e+03,  -3.2686610e+02,   2.9733049e+02);
    expectedForces[352]       = Vec3( -2.3946369e+02,   4.9753786e+02,   7.2681440e+02);
    expectedForces[353]       = Vec3(  3.9085837e+01,   7.4391931e+00,   3.4179196e+00);
    expectedForces[354]       = Vec3(  1.0462977e+03,   4.3639578e+02,  -5.7349371e+02);
    expectedForces[355]       = Vec3(  1.2715333e+02,  -5.7907885e+02,   5.7463165e+02);
    expectedForces[356]       = Vec3( -1.5968981e+02,   2.2295605e+01,   6.3554487e+01);
    expectedForces[357]       = Vec3( -1.9666111e+02,  -2.0286725e+02,  -8.3590992e+01);
    expectedForces[358]       = Vec3(  5.9985006e+00,   2.9764307e+01,   3.1304672e+01);
    expectedForces[359]       = Vec3(  1.0409792e+02,   4.3974987e+02,  -9.1530162e+02);
    expectedForces[360]       = Vec3(  2.9538724e+02,   1.4353710e+02,   2.1558768e+02);
    expectedForces[361]       = Vec3( -3.7974445e+02,  -9.0749144e+01,   1.4531502e+02);
    expectedForces[362]       = Vec3(  3.7089093e+02,   2.2523397e+02,  -2.9760698e+02);
    expectedForces[363]       = Vec3( -7.6884906e+01,  -2.7486958e+02,   3.5671043e+02);
    expectedForces[364]       = Vec3(  4.3400097e+02,   2.2275226e+02,   5.9201087e+01);
    expectedForces[365]       = Vec3(  2.1177620e+02,  -6.4990517e+02,  -1.8431623e+02);
    expectedForces[366]       = Vec3( -8.2234148e+01,   4.9900036e+02,   3.0906858e+02);
    expectedForces[367]       = Vec3( -1.4784097e+02,   1.3442463e+02,  -2.6811487e+02);
    expectedForces[368]       = Vec3(  2.2289541e+02,  -4.1933983e+02,   1.3646246e+01);
    expectedForces[369]       = Vec3( -2.1199481e+02,   6.2103281e+01,   1.6049768e+02);
    expectedForces[370]       = Vec3( -2.9429845e+01,  -1.3944740e+01,   1.6000003e+01);
    expectedForces[371]       = Vec3( -1.5093250e+02,  -1.9541345e+02,  -1.5879844e+02);
    expectedForces[372]       = Vec3(  2.2218133e+02,   7.5740668e+01,  -7.0754853e+02);
    expectedForces[373]       = Vec3( -4.7529615e+02,  -3.5259490e+02,   3.8025744e+02);
    expectedForces[374]       = Vec3(  4.9308495e+01,   1.2206152e+02,   4.2704595e+01);
    expectedForces[375]       = Vec3( -3.6408113e+02,   9.5905960e+02,   3.3861668e+02);
    expectedForces[376]       = Vec3(  1.5324093e+02,  -1.4375341e+02,  -9.1369900e+01);
    expectedForces[377]       = Vec3( -3.5343609e+02,  -4.4662438e+02,   3.1741831e+02);
    expectedForces[378]       = Vec3( -8.2721281e+02,  -6.4969313e+01,   4.1869320e+02);
    expectedForces[379]       = Vec3( -8.7419096e+00,  -4.6278343e+01,  -3.0259942e+01);
    expectedForces[380]       = Vec3(  1.6814871e+01,   4.5462601e+01,  -3.5204538e+01);
    expectedForces[381]       = Vec3(  1.9043143e+02,   5.1101135e+02,   2.9970918e+02);
    expectedForces[382]       = Vec3( -4.3547560e+02,  -7.7495770e+01,  -8.2311440e+01);
    expectedForces[383]       = Vec3(  1.3206505e+02,   4.0993476e+01,   1.5150623e+01);
    expectedForces[384]       = Vec3(  4.6163268e+01,  -2.7722119e+01,  -2.2996255e+02);
    expectedForces[385]       = Vec3( -4.7033552e+01,   1.5721545e+02,  -9.2562636e+01);
    expectedForces[386]       = Vec3( -2.4164903e+02,  -1.4771464e+02,  -5.1675638e+01);
    expectedForces[387]       = Vec3(  1.3445564e+03,   3.5387081e+02,  -3.7402713e+02);
    expectedForces[388]       = Vec3( -2.6290280e+02,   2.5025510e+02,   3.0570707e+01);
    expectedForces[389]       = Vec3( -4.4624214e+02,  -2.4919412e+02,   1.1336564e+03);
    expectedForces[390]       = Vec3(  3.1188483e+02,   6.4836125e+02,   4.5477683e+01);
    expectedForces[391]       = Vec3(  2.5570674e+02,  -2.7101115e+02,  -1.2579817e+02);
    expectedForces[392]       = Vec3( -3.4941759e+02,  -3.2508734e+01,  -5.6131515e+01);
    expectedForces[393]       = Vec3( -1.2734301e+02,  -2.0417657e+02,   5.7738004e+02);
    expectedForces[394]       = Vec3(  5.0435795e+01,  -9.5203527e+01,  -4.1237509e+00);
    expectedForces[395]       = Vec3( -6.8258463e+01,   1.7475914e+02,  -1.6422294e+02);
    expectedForces[396]       = Vec3(  2.6455575e+02,   2.6880230e+02,  -1.6607620e+01);
    expectedForces[397]       = Vec3( -3.5968167e+02,   6.6092937e+01,   1.5915445e+02);
    expectedForces[398]       = Vec3(  8.4092494e-01,  -5.9896731e+01,  -2.1856007e+01);
    expectedForces[399]       = Vec3( -1.9237984e+02,  -1.6506355e+02,   1.3370845e+02);
    expectedForces[400]       = Vec3(  4.7023043e+02,   1.2931257e+02,  -1.7013618e+02);
    expectedForces[401]       = Vec3(  9.4928843e+01,  -3.1613892e+02,  -4.6973091e+02);
    expectedForces[402]       = Vec3( -2.9180818e+02,   6.5496755e+02,  -2.9430870e+02);
    expectedForces[403]       = Vec3(  1.6670668e+02,  -2.4693913e+02,  -8.8775735e+01);
    expectedForces[404]       = Vec3( -1.4218604e+02,   2.3173349e+02,   3.9869740e+02);
    expectedForces[405]       = Vec3(  4.3469431e+02,  -2.4686568e+02,  -2.6838523e+02);
    expectedForces[406]       = Vec3(  1.0671130e+02,   3.9710106e+02,   2.6478100e+02);
    expectedForces[407]       = Vec3( -4.4281099e+02,  -3.2746126e+01,  -2.0719929e+02);
    expectedForces[408]       = Vec3(  1.1864952e+01,  -1.5901568e+02,   1.1865945e+01);
    expectedForces[409]       = Vec3( -6.0546863e+01,  -1.3488005e+02,  -2.9223986e+02);
    expectedForces[410]       = Vec3(  2.5845942e+01,   1.7545933e+02,  -1.4590654e+02);
    expectedForces[411]       = Vec3( -3.0776460e+02,   2.7351986e+02,  -2.0386483e+02);
    expectedForces[412]       = Vec3(  1.7301558e+02,   8.3012339e+01,   7.6018513e+02);
    expectedForces[413]       = Vec3(  1.2893147e+02,  -3.9937475e+00,   1.4197628e+01);
    expectedForces[414]       = Vec3( -3.9803615e+02,  -4.4825663e+02,   2.8616211e+01);
    expectedForces[415]       = Vec3( -2.3586809e+02,  -8.9306136e+01,   7.5265669e+02);
    expectedForces[416]       = Vec3(  8.3814512e+02,  -9.5563706e+00,   2.8738759e+02);
    expectedForces[417]       = Vec3( -3.8269399e+02,  -2.3650170e+02,   1.7330415e+01);
    expectedForces[418]       = Vec3(  7.2310675e+00,   2.4731653e+01,   1.9491126e+01);
    expectedForces[419]       = Vec3(  8.4485111e+01,   9.0560261e+01,  -1.4261777e+01);
    expectedForces[420]       = Vec3(  2.1829925e+02,  -8.3951651e+01,  -7.8910117e+02);
    expectedForces[421]       = Vec3( -2.4228018e+02,  -5.7097512e+01,   2.6447829e+02);
    expectedForces[422]       = Vec3(  6.9666665e+01,  -4.7313678e+02,   4.5559673e+02);
    expectedForces[423]       = Vec3(  1.0108369e+02,   4.6568610e+02,   2.5612541e+02);
    expectedForces[424]       = Vec3(  7.2334238e+01,  -1.5660040e+02,  -2.6278686e+02);
    expectedForces[425]       = Vec3(  1.0749954e+03,   2.4395658e+02,  -8.7308262e+01);
    expectedForces[426]       = Vec3(  5.4983082e+02,   5.6528577e+02,  -3.0489991e+02);
    expectedForces[427]       = Vec3( -2.1210132e+02,  -1.5871954e+02,  -1.6936452e+02);
    expectedForces[428]       = Vec3( -9.5473983e+01,   7.4879276e+01,   5.5720204e+01);
    expectedForces[429]       = Vec3(  5.6417182e+02,  -1.3478250e+02,  -4.8006673e+01);
    expectedForces[430]       = Vec3( -1.0055266e+02,  -8.3811289e+01,   6.5259748e+01);
    expectedForces[431]       = Vec3(  3.0553989e+01,   4.0927300e+02,   6.2242246e+02);
    expectedForces[432]       = Vec3( -6.0956353e+02,   9.0301312e+02,  -1.4063582e+01);
    expectedForces[433]       = Vec3(  3.2774559e+02,  -1.0849473e+02,  -8.8877712e+01);
    expectedForces[434]       = Vec3(  2.1371643e+02,  -7.2067980e+01,   2.4803863e+02);
    expectedForces[435]       = Vec3( -9.0636479e+01,   1.0673514e+03,   1.3005625e+02);
    expectedForces[436]       = Vec3( -3.5058624e+02,  -1.4968284e+02,   1.4340635e+02);
    expectedForces[437]       = Vec3(  3.4640813e+01,  -4.8035225e+01,   4.4066226e+01);
    expectedForces[438]       = Vec3( -2.3414523e+02,  -9.1359130e+01,   2.8535427e+02);
    expectedForces[439]       = Vec3( -2.7156547e+02,  -6.3426763e+01,   1.6968017e+01);
    expectedForces[440]       = Vec3( -2.8751166e+00,   4.3192169e+00,  -2.1401359e+00);
    expectedForces[441]       = Vec3( -2.0241419e+01,   4.2856138e+02,   7.3204868e+01);
    expectedForces[442]       = Vec3( -4.2166169e+01,  -2.2454013e+01,  -5.4998450e+01);
    expectedForces[443]       = Vec3(  2.6442026e+02,  -1.1670567e+02,   2.8879291e+02);
    expectedForces[444]       = Vec3( -4.9762729e+02,   4.7320314e+02,   3.2506251e+02);
    expectedForces[445]       = Vec3( -3.0994232e+02,   3.6728867e+02,  -6.8793662e+02);
    expectedForces[446]       = Vec3( -5.6328374e+00,  -2.3284824e+00,   1.9650486e+01);
    expectedForces[447]       = Vec3(  4.0145354e+02,   7.0143762e+02,   1.3039041e+02);
    expectedForces[448]       = Vec3(  1.3270713e+02,  -3.3050447e+02,   2.8630784e+02);
    expectedForces[449]       = Vec3( -3.8904836e+02,  -4.6125118e+01,  -3.8815388e+02);
    expectedForces[450]       = Vec3(  4.6293066e+02,  -2.4970158e+02,  -5.7231818e+01);
    expectedForces[451]       = Vec3( -4.7707142e+00,   1.7120822e-01,   4.5081299e+00);
    expectedForces[452]       = Vec3(  3.5968818e+01,  -1.1892248e+01,   1.2184792e+02);
    expectedForces[453]       = Vec3(  4.4939172e+02,   4.8377647e+02,   3.0529778e+02);
    expectedForces[454]       = Vec3( -2.9071570e+02,  -9.0611613e+01,  -2.5190431e+02);
    expectedForces[455]       = Vec3( -1.4344537e+02,  -1.3061997e+03,   3.8343947e+02);
    expectedForces[456]       = Vec3(  4.1580940e+02,  -3.6864043e+02,   6.6384066e+02);
    expectedForces[457]       = Vec3( -3.4869949e+02,   3.3262427e+02,  -8.3156358e+01);
    expectedForces[458]       = Vec3( -5.8565457e+00,  -5.1879987e+02,  -7.2717213e+02);
    expectedForces[459]       = Vec3( -1.0514593e+02,  -2.6844396e+02,   3.0849659e+02);
    expectedForces[460]       = Vec3( -1.2562153e+02,   6.5438031e+01,  -2.3521687e+02);
    expectedForces[461]       = Vec3(  8.3330016e+02,  -6.7961048e+02,  -7.2442064e+02);
    expectedForces[462]       = Vec3( -1.3328298e+02,   1.6143544e+02,  -3.2000493e+02);
    expectedForces[463]       = Vec3( -1.1035002e+02,  -8.6040775e+01,   2.0953986e+02);
    expectedForces[464]       = Vec3( -2.7702219e+02,  -1.2156319e+02,  -3.8897654e+02);
    expectedForces[465]       = Vec3(  3.6277432e+02,  -7.0422679e+02,  -8.2487472e+02);
    expectedForces[466]       = Vec3( -2.3947587e+02,  -4.3625964e+00,   4.8716989e+01);
    expectedForces[467]       = Vec3( -2.6210552e+02,   9.9730213e+01,  -2.4465307e+02);
    expectedForces[468]       = Vec3( -4.1344138e+01,  -3.8636310e+01,  -1.0720610e+01);
    expectedForces[469]       = Vec3(  2.3928584e+01,  -1.3030908e+01,  -1.7631168e+01);
    expectedForces[470]       = Vec3( -1.1873594e+02,  -3.3905287e+00,  -8.4575994e+01);
    expectedForces[471]       = Vec3(  4.3384846e+02,  -4.8660929e+02,  -2.4015737e+02);
    expectedForces[472]       = Vec3(  3.8974922e+01,   2.4401913e+02,  -1.2060132e+02);
    expectedForces[473]       = Vec3( -4.2239647e+02,  -8.6924630e+01,   5.7884013e+02);
    expectedForces[474]       = Vec3(  1.7222431e+02,  -1.9292167e+01,   4.4184770e+01);
    expectedForces[475]       = Vec3( -4.8871102e+01,   1.8185057e+02,  -1.3144294e+02);
    expectedForces[476]       = Vec3( -1.1467881e+02,  -1.8387284e+02,  -2.6935761e+01);
    expectedForces[477]       = Vec3(  4.4415601e+02,  -8.6415426e+02,  -2.2765803e+02);
    expectedForces[478]       = Vec3(  2.0893474e+02,   1.3529112e+02,   2.4108004e+02);
    expectedForces[479]       = Vec3( -1.2802570e+01,  -4.7089261e+00,  -7.0024810e+00);
    expectedForces[480]       = Vec3( -1.0878354e+02,   4.5030357e+02,  -3.9086153e+02);
    expectedForces[481]       = Vec3( -3.7659745e+01,  -4.9295745e+01,  -6.0028516e+01);
    expectedForces[482]       = Vec3(  1.8851809e+02,  -1.6863914e+02,   1.2031505e+02);
    expectedForces[483]       = Vec3(  2.2040686e+02,  -1.1953187e+03,   2.9106238e+02);
    expectedForces[484]       = Vec3( -1.1642847e+01,   5.6739977e+01,   3.7532413e+00);
    expectedForces[485]       = Vec3(  3.5013234e+00,   2.8236323e+01,  -1.1738866e+01);
    expectedForces[486]       = Vec3( -2.0051506e+03,  -5.3138994e+02,   4.0428126e+02);
    expectedForces[487]       = Vec3(  1.6887843e+02,   1.1687735e+02,  -9.6232177e+01);
    expectedForces[488]       = Vec3(  2.5030722e+02,  -3.6662537e+02,  -5.8595505e+01);
    expectedForces[489]       = Vec3( -1.5924273e+01,  -2.1761920e+02,  -3.4885425e+01);
    expectedForces[490]       = Vec3(  1.0423157e+02,   4.0767827e+01,   1.4036192e+02);
    expectedForces[491]       = Vec3(  5.0090776e+02,  -5.0455702e+01,  -4.2524815e+02);
    expectedForces[492]       = Vec3(  3.8037844e+02,   5.7460724e+02,   2.1243310e+02);
    expectedForces[493]       = Vec3( -1.1487653e+02,   3.3550952e+01,  -5.7050397e+01);
    expectedForces[494]       = Vec3(  1.8425697e+02,  -2.7223857e+02,  -2.9384020e+01);
    expectedForces[495]       = Vec3( -1.1115514e+02,  -4.4519957e+02,  -3.3628896e+02);
    expectedForces[496]       = Vec3(  3.8865823e+02,  -2.1372852e+02,  -7.8979988e+02);
    expectedForces[497]       = Vec3( -3.1740258e+02,   4.4824683e+02,  -5.0382680e+02);
    expectedForces[498]       = Vec3(  1.5521684e+02,  -1.5278429e+02,  -1.8418006e+02);
    expectedForces[499]       = Vec3( -2.2213191e+01,   7.1015095e+00,   1.0308820e+01);
    expectedForces[500]       = Vec3(  3.0011560e+01,   2.4584541e+02,  -5.3231694e+02);
    expectedForces[501]       = Vec3( -8.4865171e+02,  -2.0279022e+02,  -6.8435095e+02);
    expectedForces[502]       = Vec3(  7.0477549e+00,   2.7626774e+01,  -1.6246047e+01);
    expectedForces[503]       = Vec3( -7.1309648e+01,  -1.6054218e+02,  -3.1621746e+01);
    expectedForces[504]       = Vec3(  3.6317856e+02,   1.9055487e+02,   3.9196046e+01);
    expectedForces[505]       = Vec3(  1.4643265e+02,  -1.2295335e+03,  -2.2268806e+01);
    expectedForces[506]       = Vec3( -3.6638240e+02,  -6.5310612e+00,  -6.8077013e+02);
    expectedForces[507]       = Vec3(  1.9107813e+02,  -6.2176534e+02,   5.1826941e+02);
    expectedForces[508]       = Vec3(  1.7226933e+02,  -3.8939126e+02,  -7.4174355e+02);
    expectedForces[509]       = Vec3(  1.3541890e+02,   4.1350561e+02,  -2.6155396e+02);
    expectedForces[510]       = Vec3(  1.6452609e+02,  -4.1963597e+02,   1.6322800e+02);
    expectedForces[511]       = Vec3(  9.1248387e+01,  -1.3380408e+01,  -5.1326806e+02);
    expectedForces[512]       = Vec3(  1.2420564e+02,   8.9492432e+02,   2.8660661e+02);
    expectedForces[513]       = Vec3(  4.5609859e+02,   3.4252036e+01,  -2.3830457e+02);
    expectedForces[514]       = Vec3( -8.8282131e+01,  -1.4512876e+02,  -4.5703642e+01);
    expectedForces[515]       = Vec3( -6.9194714e+01,   5.5053415e+02,  -3.9527911e+02);
    expectedForces[516]       = Vec3( -1.5570798e+02,   4.1650792e-02,   3.0224667e+02);
    expectedForces[517]       = Vec3(  9.2007640e+01,   3.6378499e+02,  -2.2069415e+01);
    expectedForces[518]       = Vec3(  2.7996922e+02,  -2.0932322e+02,   6.2987796e+01);
    expectedForces[519]       = Vec3( -7.4985282e+01,   1.3508916e+02,   1.1453750e+02);
    expectedForces[520]       = Vec3(  1.2149786e+01,   4.2161047e+02,  -2.9521571e+02);
    expectedForces[521]       = Vec3(  1.9419524e+01,  -8.6619683e+01,  -9.4040609e+01);
    expectedForces[522]       = Vec3(  3.5495800e+01,   6.4028894e+01,   1.2320699e+02);
    expectedForces[523]       = Vec3( -5.6501302e+02,  -1.2497462e+02,   8.1633096e+02);
    expectedForces[524]       = Vec3( -3.1775230e+02,   1.6467221e+02,  -4.8242287e+01);
    expectedForces[525]       = Vec3( -1.8325245e+02,  -1.0751941e+02,   6.2172117e+02);
    expectedForces[526]       = Vec3(  1.8164299e+00,   8.0654784e-01,  -1.4484741e+01);
    expectedForces[527]       = Vec3(  2.6686986e+00,   1.5933475e+01,  -6.8843007e+01);
    expectedForces[528]       = Vec3( -1.3208642e+02,  -1.2281791e+02,   7.8417959e+01);
    expectedForces[529]       = Vec3( -4.8078193e+01,  -1.5140963e+02,   8.3605002e+01);
    expectedForces[530]       = Vec3(  8.5925332e+01,   1.9959621e+01,  -5.9034655e+01);
    expectedForces[531]       = Vec3( -3.1311626e+00,  -3.4009931e+01,  -1.5895685e+02);
    expectedForces[532]       = Vec3( -9.9694666e+01,   1.0296874e+00,   1.0482739e+02);
    expectedForces[533]       = Vec3(  2.8954747e+02,   1.2018666e+02,   2.2155006e+02);
    expectedForces[534]       = Vec3(  3.2588176e+02,   3.6622498e+01,   3.5801460e+02);
    expectedForces[535]       = Vec3( -6.5985617e+01,  -3.5501709e+02,   6.8149810e+01);
    expectedForces[536]       = Vec3( -4.7632753e+02,   2.7936656e+02,   1.5548558e+02);
    expectedForces[537]       = Vec3(  4.7724904e+02,   2.1647106e+02,   1.2003317e+02);
    expectedForces[538]       = Vec3( -1.6225405e+02,   1.4264998e+02,  -1.0113321e+02);
    expectedForces[539]       = Vec3( -3.2540487e+01,  -7.5643223e+01,   1.6054148e+01);
    expectedForces[540]       = Vec3(  4.3311991e+02,   5.8082595e+02,   1.9354276e+02);
    expectedForces[541]       = Vec3( -3.4924430e+02,  -7.0056069e+01,   1.0274560e+02);
    expectedForces[542]       = Vec3(  1.9441645e+02,  -2.0017354e+02,  -2.3280717e+02);
    expectedForces[543]       = Vec3( -4.1530380e+01,  -5.6394351e+02,   2.4472509e+02);
    expectedForces[544]       = Vec3( -1.6193396e+01,  -8.0431396e+01,  -2.9094018e+02);
    expectedForces[545]       = Vec3( -1.9414773e+01,   4.6180982e+01,   3.3123072e+01);
    expectedForces[546]       = Vec3( -3.9314039e+02,  -3.9874866e+02,   4.5308571e+02);
    expectedForces[547]       = Vec3( -7.1897482e+01,  -1.1940445e+02,  -2.7405931e+02);
    expectedForces[548]       = Vec3(  3.0646396e+02,   6.1235747e+01,  -1.5253270e+02);
    expectedForces[549]       = Vec3( -3.2480464e+02,   4.3056561e+02,   3.2532485e+01);
    expectedForces[550]       = Vec3(  1.2818655e+02,   7.4294994e-01,  -5.7650521e+00);
    expectedForces[551]       = Vec3( -1.7481437e+02,  -2.6203225e+02,  -9.6793481e+01);
    expectedForces[552]       = Vec3(  1.9259110e+02,   3.0171883e+02,   5.2403235e+02);
    expectedForces[553]       = Vec3( -7.8132123e+01,  -7.6569340e+00,  -1.1140240e+02);
    expectedForces[554]       = Vec3( -5.2504673e+02,  -4.3700574e+02,   4.3321261e+02);
    expectedForces[555]       = Vec3( -9.9785283e+02,   2.7139143e+02,  -4.1723547e+02);
    expectedForces[556]       = Vec3(  2.7115933e+02,  -8.7444541e+01,   1.0745103e+02);
    expectedForces[557]       = Vec3(  1.4348018e+02,   1.4013343e+02,  -4.0305733e+02);
    expectedForces[558]       = Vec3(  4.9491010e+02,   4.3040089e+02,  -3.5558583e+02);
    expectedForces[559]       = Vec3( -1.8510736e+02,   6.7129731e+01,  -2.1324364e+02);
    expectedForces[560]       = Vec3(  8.4478090e+01,   4.0986269e+01,   4.1401094e+02);
    expectedForces[561]       = Vec3(  4.4188419e+01,  -8.6520108e+02,   6.8555821e+02);
    expectedForces[562]       = Vec3( -4.5036974e+02,   2.2379481e+02,   6.2457971e+01);
    expectedForces[563]       = Vec3(  2.4325836e+02,   6.9725116e+01,   3.6266296e+01);
    expectedForces[564]       = Vec3(  1.5611925e+02,  -2.4471023e+02,   4.7857332e+02);
    expectedForces[565]       = Vec3( -3.5692822e+02,  -2.7248961e+02,  -1.8165146e+02);
    expectedForces[566]       = Vec3(  4.3093388e+02,   1.6944839e+02,  -4.4177741e+02);
    expectedForces[567]       = Vec3(  5.2118711e+02,  -1.5252783e+02,  -1.6964047e+02);
    expectedForces[568]       = Vec3(  1.3899083e+03,   3.8363654e+02,  -6.4623177e+02);
    expectedForces[569]       = Vec3( -1.4984858e+02,   4.5097170e+01,  -5.0257768e+01);
    expectedForces[570]       = Vec3( -8.1979709e+01,  -4.1632909e+01,  -5.1367825e+02);
    expectedForces[571]       = Vec3( -7.8645235e+00,   1.8105470e+01,   9.2902759e+01);
    expectedForces[572]       = Vec3(  1.2725548e+02,   1.3926229e+02,  -4.3752260e+01);
    expectedForces[573]       = Vec3( -2.7817289e+01,  -8.0900345e+01,  -4.3178498e+02);
    expectedForces[574]       = Vec3(  2.2115263e+02,  -7.8807799e+01,   8.8367391e+01);
    expectedForces[575]       = Vec3( -1.2756938e+02,   4.2170937e+02,  -6.8418245e+01);
    expectedForces[576]       = Vec3(  8.8370341e+01,  -1.4671130e+02,  -5.6211071e+02);
    expectedForces[577]       = Vec3(  6.7701453e+02,   2.1228010e+02,   2.3047236e+02);
    expectedForces[578]       = Vec3( -1.2694819e+02,  -7.2734890e+00,   1.1007893e+02);
    expectedForces[579]       = Vec3( -6.8945756e+02,   1.7359485e+02,  -2.5607776e+02);
    expectedForces[580]       = Vec3(  6.7323923e+00,   3.1574487e+02,   5.9741152e+02);
    expectedForces[581]       = Vec3(  2.1214031e+02,  -4.7565197e+01,  -9.6896159e+01);
    expectedForces[582]       = Vec3( -3.4876562e+02,   5.8335489e+01,  -1.7058451e+02);
    expectedForces[583]       = Vec3( -1.6516914e+02,  -3.8913162e+02,   5.8832025e+02);
    expectedForces[584]       = Vec3(  4.7753612e+01,  -6.6792213e+01,  -4.6492749e+01);
    expectedForces[585]       = Vec3( -4.7166578e+02,   7.1811831e+02,   7.0620222e+02);
    expectedForces[586]       = Vec3( -3.0570577e+02,  -1.7681907e+02,  -2.1227370e+02);
    expectedForces[587]       = Vec3(  1.8122997e+02,  -6.8318737e+01,  -3.4853368e+02);
    expectedForces[588]       = Vec3(  7.5752093e+01,   3.9702196e+02,  -9.0196404e+02);
    expectedForces[589]       = Vec3(  5.1185528e+02,  -7.7747186e+02,  -3.3969581e+02);
    expectedForces[590]       = Vec3( -3.1187606e+01,   5.9593198e+01,   7.7918649e+01);
    expectedForces[591]       = Vec3(  2.2741676e+02,   2.1705631e+02,   3.9920149e+02);
    expectedForces[592]       = Vec3(  7.0516894e+01,   3.3389665e+02,  -6.8443760e+00);
    expectedForces[593]       = Vec3( -9.9645406e+02,   2.8574127e+02,  -8.0946834e+02);
    expectedForces[594]       = Vec3(  1.0930460e+02,   1.5869363e+03,  -4.0704089e+02);
    expectedForces[595]       = Vec3( -3.0341438e+02,  -5.7494385e+00,   2.9655509e+02);
    expectedForces[596]       = Vec3(  8.7955203e+01,  -7.6811348e+01,  -6.9558652e+01);
    expectedForces[597]       = Vec3( -9.4378353e+01,   3.6253070e+02,   9.0703343e-01);
    expectedForces[598]       = Vec3( -4.7595525e+01,  -6.9148379e+02,  -4.6737971e+02);
    expectedForces[599]       = Vec3(  2.2339987e+02,  -2.5503335e+02,   2.8552040e+02);
    expectedForces[600]       = Vec3(  1.0309072e+02,  -2.1298204e+02,  -4.3393283e+02);
    expectedForces[601]       = Vec3( -2.4581477e+02,   2.5126386e+02,  -1.8679710e+02);
    expectedForces[602]       = Vec3(  7.7926618e+00,  -2.4947558e+01,   1.5962011e+02);
    expectedForces[603]       = Vec3(  7.0099476e+02,   5.9496129e+01,   2.0314640e+02);
    expectedForces[604]       = Vec3( -5.2583810e+01,   3.1237489e+01,  -4.0902937e+00);
    expectedForces[605]       = Vec3( -5.4435871e+01,  -7.7056693e+01,  -1.5137215e+01);
    expectedForces[606]       = Vec3( -7.8684074e+02,   8.1673860e+02,   5.0192222e+02);
    expectedForces[607]       = Vec3(  2.9859285e+02,   1.9661708e+02,  -1.9940993e+02);
    expectedForces[608]       = Vec3( -9.9358476e-01,  -3.4789219e+02,   1.7693106e+02);
    expectedForces[609]       = Vec3(  6.7880927e+02,   7.2768533e+01,   1.1973248e+03);
    expectedForces[610]       = Vec3( -1.2072838e+00,   1.1770157e+02,  -8.6765918e+01);
    expectedForces[611]       = Vec3( -4.3236695e+02,  -2.8583002e+02,   1.3459508e+02);
    expectedForces[612]       = Vec3( -6.1894609e+01,   4.4705141e+02,  -2.3264215e+02);
    expectedForces[613]       = Vec3( -2.6948720e+01,  -2.7688230e+01,   4.6300601e+01);
    expectedForces[614]       = Vec3(  4.3954815e+01,   6.1244278e+01,   1.0727211e+02);
    expectedForces[615]       = Vec3( -2.2049132e+02,  -9.7233594e+01,   8.3807949e+01);
    expectedForces[616]       = Vec3(  2.8631226e+02,  -1.3478640e+02,  -5.0675977e+02);
    expectedForces[617]       = Vec3( -1.1458029e+01,   6.3762007e+01,  -4.9643140e+01);
    expectedForces[618]       = Vec3(  7.7363305e+01,   1.1936816e+02,   3.9692323e+01);
    expectedForces[619]       = Vec3( -8.4199858e+01,   2.3567494e+02,   7.2885200e+01);
    expectedForces[620]       = Vec3( -7.4500104e+00,   3.5567340e+00,  -2.2676267e+01);
    expectedForces[621]       = Vec3( -7.9200057e+02,   3.2245192e+01,   1.1829644e+01);
    expectedForces[622]       = Vec3(  3.9959336e+02,   2.2447865e+02,  -3.1413233e+02);
    expectedForces[623]       = Vec3(  1.5810527e+02,  -2.2341997e+02,   2.3388212e+02);
    expectedForces[624]       = Vec3( -2.9248657e+02,  -1.0806250e+03,  -6.1268035e+01);
    expectedForces[625]       = Vec3(  6.5061442e+02,   1.2244756e+02,  -1.1976147e+02);
    expectedForces[626]       = Vec3( -1.7787710e+02,   3.7279354e+02,  -3.7911745e+02);
    expectedForces[627]       = Vec3(  3.7145927e+02,  -2.5203732e+02,  -2.2889342e+02);
    expectedForces[628]       = Vec3(  1.9440281e+02,  -2.3076142e+02,   7.2945279e+01);
    expectedForces[629]       = Vec3( -2.3603693e+02,  -1.1987203e+02,  -1.2000208e+02);
    expectedForces[630]       = Vec3( -8.4811726e+01,  -6.4500846e+02,  -5.1257210e+02);
    expectedForces[631]       = Vec3( -1.6183393e+01,   1.3980655e+01,  -2.3189702e+00);
    expectedForces[632]       = Vec3(  1.3526098e+02,   1.3800365e+02,  -9.4360476e+01);
    expectedForces[633]       = Vec3(  2.0983514e+01,   5.1156494e+02,  -4.3526557e+01);
    expectedForces[634]       = Vec3( -2.8666500e+02,  -5.8235691e+02,  -8.7425194e+01);
    expectedForces[635]       = Vec3( -1.7183604e+02,   1.0002944e+01,   6.5855812e+02);
    expectedForces[636]       = Vec3( -2.1735465e+02,   9.4664675e+02,   8.9132890e+02);
    expectedForces[637]       = Vec3(  6.3759184e+01,  -5.2436201e+02,   1.2455578e+02);
    expectedForces[638]       = Vec3(  3.2176621e+02,   1.8182348e+02,  -4.4863340e+02);
    expectedForces[639]       = Vec3( -8.4308252e+02,  -1.6268681e+02,   5.8113352e+02);
    expectedForces[640]       = Vec3(  6.1223389e+01,   3.6767423e+02,   5.9071999e+02);
    expectedForces[641]       = Vec3(  2.2546928e+02,  -2.2022805e+02,   8.3027103e+01);
    expectedForces[642]       = Vec3( -4.0385190e+02,   2.0141033e+02,   2.7296512e+01);
    expectedForces[643]       = Vec3(  8.4204485e+01,  -7.1102294e+01,   5.7642677e+01);
    expectedForces[644]       = Vec3(  6.1859901e+01,   1.4828523e+02,  -2.6764016e+02);
    expectedForces[645]       = Vec3( -6.7249932e+02,  -8.4613525e+01,  -3.9792609e+02);
    expectedForces[646]       = Vec3(  9.9116485e+01,  -3.7714583e+01,  -5.3966332e+01);
    expectedForces[647]       = Vec3(  2.0868628e+02,   2.9747206e+02,   3.3931416e+02);

    // tolerance is higher here due to interpolation used in setting tapering coefficients;
    // if tapering turned off, then absolute difference < 2.0e-05

    double tolerance          = 5.0e-04;
    compareForcesEnergy(testName, expectedEnergy, energy, expectedForces, forces, tolerance);

    // test sigma/epsilon rules for dispersion correction

    if (includeVdwDispersionCorrection) {

         std::vector<std::string> sigmaRules;
         std::vector<std::string> epsilonRules;
         std::vector<double> expectedEnergies;

         sigmaRules.push_back("ARITHMETIC");
         epsilonRules.push_back("ARITHMETIC");
         expectedEnergies.push_back(6.2137988e+03);

         sigmaRules.push_back("GEOMETRIC");
         epsilonRules.push_back("GEOMETRIC");
         expectedEnergies.push_back( 3.6358216e+03);

         sigmaRules.push_back("CUBIC-MEAN");
         epsilonRules.push_back("HARMONIC");
         expectedEnergies.push_back(3.2774624e+03);

         for (unsigned int ii = 0; ii < sigmaRules.size(); ii++) {
             setupAndGetForcesEnergyVdwWater(sigmaRules[ii], epsilonRules[ii], cutoff, boxDimension, includeVdwDispersionCorrection, forces, energy);
             testName    = "testVdwWaterWithDispersionCorrection_" + sigmaRules[ii] + '_' + epsilonRules[ii];
             ASSERT_EQUAL_TOL_MOD(expectedEnergies[ii], energy, tolerance, testName);
         }
 
    }
}

void testTriclinic() {
    System system;
    system.addParticle(1.0);
    system.addParticle(1.0);
    Vec3 a(3.1, 0, 0);
    Vec3 b(0.4, 3.5, 0);
    Vec3 c(-0.1, -0.5, 4.0);
    system.setDefaultPeriodicBoxVectors(a, b, c);
    LangevinIntegrator integrator(0.0, 0.1, 0.01);
    AmoebaVdwForce* vdw = new AmoebaVdwForce();
    vdw->setUseDispersionCorrection(false);
    vdw->addParticle(0, 0.5, 1.0, 0.0);
    vdw->addParticle(1, 0.5, 1.0, 0.0);
    vdw->setNonbondedMethod(AmoebaVdwForce::CutoffPeriodic);
    const double cutoff = 1.5;
    vdw->setCutoff(cutoff);
    system.addForce(vdw);
    Context context(system, integrator, Platform::getPlatformByName("CUDA"));
    vector<Vec3> positions(2);
    OpenMM_SFMT::SFMT sfmt;
    init_gen_rand(0, sfmt);
    for (int iteration = 0; iteration < 50; iteration++) {
        // Generate random positions for the two particles.

        positions[0] = a*genrand_real2(sfmt) + b*genrand_real2(sfmt) + c*genrand_real2(sfmt);
        positions[1] = a*genrand_real2(sfmt) + b*genrand_real2(sfmt) + c*genrand_real2(sfmt);
        context.setPositions(positions);

        // Loop over all possible periodic copies and find the nearest one.

        Vec3 delta;
        double distance2 = 100.0;
        for (int i = -1; i < 2; i++)
            for (int j = -1; j < 2; j++)
                for (int k = -1; k < 2; k++) {
                    Vec3 d = positions[1]-positions[0]+a*i+b*j+c*k;
                    if (d.dot(d) < distance2) {
                        delta = d;
                        distance2 = d.dot(d);
                    }
                }
        double distance = sqrt(distance2);

        // See if the energy is correct.

        State state = context.getState(State::Energy);
        if (distance >= cutoff) {
            ASSERT_EQUAL(0.0, state.getPotentialEnergy());
        }
        else if (distance < 0.9*cutoff) {
            const double energy = pow(1.07/(distance+0.07), 7.0)*(1.12/(pow(distance, 7.0)+0.12)-2);
            ASSERT_EQUAL_TOL(energy, state.getPotentialEnergy(), TOL);
        }
    }
}

int main(int argc, char* argv[]) {
    try {
        std::cout << "TestCudaAmoebaVdwForce running test..." << std::endl;
        registerAmoebaCudaKernelFactories();
        if (argc > 1)
            Platform::getPlatformByName("CUDA").setPropertyDefaultValue("Precision", std::string(argv[1]));

        testVdw();

        // tests using two ammonia molecules

        // test VDW w/ sigmaRule=CubicMean and epsilonRule=HHG

        testVdwAmmoniaCubicMeanHhg();

        // test VDW w/ sigmaRule=Arithmetic and epsilonRule=Arithmetic

        testVdwAmmoniaArithmeticArithmetic();

        // test VDW w/ sigmaRule=Geometric and epsilonRule=Geometric

        testVdwAmmoniaGeometricGeometric();

        // test VDW w/ sigmaRule=CubicMean and epsilonRule=Harmonic

        testVdwAmmoniaCubicMeanHarmonic();

        // test w/ cutoff=0.25 nm; single ixn between two particles (0 and 6); force nonzero on
        // particle 4 due to reduction applied to NH
        // the distance between 0 and 6 is ~ 0.235 so the ixn is in the tapered region

        testVdwTaper();

        // test PBC

        testVdwPBC();

        // tests based on box of water

        int includeVdwDispersionCorrection = 0;
        testVdwWater(includeVdwDispersionCorrection);

        // includes tests for various combinations of sigma/epsilon rules
        // when computing vdw dispersion correction
 
        includeVdwDispersionCorrection     = 1;
        testVdwWater(includeVdwDispersionCorrection);
        
        // test triclinic boxes

        testTriclinic();

    } catch(const std::exception& e) {
        std::cout << "exception: " << e.what() << std::endl;
        std::cout << "FAIL - ERROR.  Test failed." << std::endl;
        return 1;
    }
    std::cout << "Done" << std::endl;
    return 0;
}
