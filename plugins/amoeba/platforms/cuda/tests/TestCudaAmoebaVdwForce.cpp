/* -------------------------------------------------------------------------- *
 *                                   OpenMMAmoeba                             *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2008 Stanford University and the Authors.           *
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
 * This tests the Cuda implementation of CudaAmoebaVdwForce.
 */

#include "openmm/internal/AssertionUtilities.h"
#include "openmm/Context.h"
#include "OpenMMAmoeba.h"
#include "AmoebaTinkerParameterFile.h"
#include "openmm/System.h"
#include "openmm/AmoebaVdwForce.h"
#include "openmm/LangevinIntegrator.h"
#include <iostream>
#include <vector>
#include <stdlib.h>
#include <stdio.h>

#define ASSERT_EQUAL_TOL_MOD(expected, found, tol, testname) {double _scale_ = std::abs(expected) > 1.0 ? std::abs(expected) : 1.0; if (!(std::abs((expected)-(found))/_scale_ <= (tol))) {std::stringstream details; details << testname << " Expected "<<(expected)<<", found "<<(found); throwException(__FILE__, __LINE__, details.str());}};

#define ASSERT_EQUAL_VEC_MOD(expected, found, tol,testname) {ASSERT_EQUAL_TOL_MOD((expected)[0], (found)[0], (tol),(testname)); ASSERT_EQUAL_TOL_MOD((expected)[1], (found)[1], (tol),(testname)); ASSERT_EQUAL_TOL_MOD((expected)[2], (found)[2], (tol),(testname));};


using namespace OpenMM;
const double TOL = 1e-4;

void testVdw( FILE* log ) {

    System system;
    int numberOfParticles          = 6;
    AmoebaVdwForce* amoebaVdwForce = new AmoebaVdwForce();
    std::string sigmaCombiningRule = std::string("CUBIC-MEAN");
    amoebaVdwForce->setSigmaCombiningRule( sigmaCombiningRule );

    std::string epsilonCombiningRule = std::string("HHG");
    amoebaVdwForce->setEpsilonCombiningRule( epsilonCombiningRule );
    int classIndex = 0;
    for( int ii = 0; ii < numberOfParticles; ii++ ){
        int indexIV;
        double mass, sigma, epsilon, reduction;
        std::vector< int > exclusions;
        if( ii == 0 || ii == 3 ){
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

        if( ii < 3 ){
            exclusions.push_back ( 0 );
            exclusions.push_back ( 1 );
            exclusions.push_back ( 2 );
        } else {
            exclusions.push_back ( 3 );
            exclusions.push_back ( 4 );
            exclusions.push_back ( 5 );
        }
        system.addParticle(mass);
        amoebaVdwForce->addParticle( indexIV, classIndex, sigma, epsilon, reduction );
        amoebaVdwForce->setParticleExclusions( ii, exclusions );
    }
    LangevinIntegrator integrator(0.0, 0.1, 0.01);

    std::vector<Vec3> positions(numberOfParticles);
    std::vector<Vec3> expectedForces(numberOfParticles);
    double expectedEnergy;

    positions[0]          = Vec3( -0.254893450E+02, -0.876646600E+01,  0.174761600E+01 );
    positions[1]          = Vec3( -0.263489690E+02, -0.907798000E+01,  0.205385100E+01 );
    positions[2]          = Vec3( -0.252491680E+02, -0.949411200E+01,  0.115017600E+01 );
    positions[3]          = Vec3(  0.172827200E+01,  0.195873090E+02,  0.100059800E+01 );
    positions[4]          = Vec3(  0.129370700E+01,  0.190112810E+02,  0.169576300E+01 );
    positions[5]          = Vec3(  0.256122300E+01,  0.191601930E+02,  0.854382000E+00 );

    double offset         = 27.0;
    for( int ii = 0; ii < 3; ii++ ){
       positions[ii][0]      += offset;
       positions[ii][1]      += offset;
    }

    expectedForces[0]     = Vec3(  -0.729561040E+03,  0.425828484E+04, -0.769114213E+03 );
    expectedForces[1]     = Vec3(   0.181000041E+02,  0.328216639E+02, -0.126210511E+02 );
    expectedForces[2]     = Vec3(  -0.943743014E+00,  0.199728310E+02,  0.884567842E+00 );
    expectedForces[3]     = Vec3(   0.615734500E+01, -0.747350431E+03,  0.264726489E+03 );
    expectedForces[4]     = Vec3(   0.735772031E+03, -0.353310112E+04,  0.490066356E+03 );
    expectedForces[5]     = Vec3(  -0.295245970E+02, -0.306277797E+02,  0.260578506E+02 );

    expectedEnergy        = 0.740688488E+03;

    system.addForce(amoebaVdwForce);
    std::string platformName;
    #define AngstromToNm 0.1    
    #define CalToJoule   4.184    
    for( int ii = 0; ii < numberOfParticles; ii++ ){
        positions[ii][0] *= AngstromToNm;
        positions[ii][1] *= AngstromToNm;
        positions[ii][2] *= AngstromToNm;
    }
    for( int ii = 0; ii < amoebaVdwForce->getNumParticles();  ii++ ){
        int indexIV;
        int classIndex;
        double sigma, epsilon, reduction;
        amoebaVdwForce->getParticleParameters( ii, indexIV, classIndex, sigma, epsilon, reduction );
        sigma        *= AngstromToNm;
        epsilon      *= CalToJoule;
        amoebaVdwForce->setParticleParameters( ii, indexIV, classIndex, sigma, epsilon, reduction );
    }
    platformName = "Cuda";
    Context context(system, integrator, Platform::getPlatformByName( platformName ) );

    context.setPositions(positions);
    State state                      = context.getState(State::Forces | State::Energy);
    std::vector<Vec3> forces         = state.getForces();
    const double conversion          = -AngstromToNm/CalToJoule;

    for( unsigned int ii = 0; ii < forces.size(); ii++ ){
        forces[ii][0] *= conversion;
        forces[ii][1] *= conversion;
        forces[ii][2] *= conversion;
    }    
    expectedEnergy *= CalToJoule;

#ifdef AMOEBA_DEBUG
    if( log ){
        (void) fprintf( log, "computeAmoebaVdwForces: expected energy=%14.7e %14.7e\n", expectedEnergy, state.getPotentialEnergy() );
        for( unsigned int ii = 0; ii < forces.size(); ii++ ){
            (void) fprintf( log, "%6u [%14.7e %14.7e %14.7e]   [%14.7e %14.7e %14.7e]\n", ii,
                            expectedForces[ii][0], expectedForces[ii][1], expectedForces[ii][2], forces[ii][0], forces[ii][1], forces[ii][2] );
        }
        (void) fflush( log );
    }
#endif

    double tolerance = 1.0e-03;
    for( unsigned int ii = 0; ii < forces.size(); ii++ ){
        ASSERT_EQUAL_VEC( expectedForces[ii], forces[ii], tolerance );
    }
    ASSERT_EQUAL_TOL( expectedEnergy, state.getPotentialEnergy(), tolerance );
}

void setupAndGetForcesEnergyVdwAmmonia( const std::string& sigmaCombiningRule, const std::string& epsilonCombiningRule, double cutoff,
                                        double boxDimension, std::vector<Vec3>& forces, double& energy, FILE* log ){

    // beginning of Vdw setup

    System system;
    AmoebaVdwForce* amoebaVdwForce        = new AmoebaVdwForce();;
    int numberOfParticles                 = 8;
    amoebaVdwForce->setSigmaCombiningRule( sigmaCombiningRule );
    amoebaVdwForce->setEpsilonCombiningRule( epsilonCombiningRule );
    amoebaVdwForce->setUseNeighborList( 1 );
    amoebaVdwForce->setCutoff( cutoff );
    if( boxDimension > 0.0 ){
        Vec3 a( boxDimension, 0.0, 0.0 );
        Vec3 b( 0.0, boxDimension, 0.0 );
        Vec3 c( 0.0, 0.0, boxDimension );
        system.setDefaultPeriodicBoxVectors( a, b, c );
        amoebaVdwForce->setPBC( 1 );
    } else {
        amoebaVdwForce->setPBC( 0 );
    }

    // addParticle: ivIndex, radius, epsilon, reductionFactor

    int classIndex = 0;
    system.addParticle(   1.4007000e+01 );
    amoebaVdwForce->addParticle( 0, classIndex,   1.8550000e-01,   4.3932000e-01,   0.0000000e+00 );

    system.addParticle(   1.0080000e+00 );
    amoebaVdwForce->addParticle( 0, classIndex,   1.3500000e-01,   8.3680000e-02,   9.1000000e-01 );

    system.addParticle(   1.0080000e+00 );
    amoebaVdwForce->addParticle( 0, classIndex,   1.3500000e-01,   8.3680000e-02,   9.1000000e-01 );

    system.addParticle(   1.0080000e+00 );
    amoebaVdwForce->addParticle( 0, classIndex,   1.3500000e-01,   8.3680000e-02,   9.1000000e-01 );

    system.addParticle(   1.4007000e+01 );
    amoebaVdwForce->addParticle( 4, classIndex,   1.8550000e-01,   4.3932000e-01,   0.0000000e+00 );

    system.addParticle(   1.0080000e+00 );
    amoebaVdwForce->addParticle( 4, classIndex,   1.3500000e-01,   8.3680000e-02,   9.1000000e-01 );

    system.addParticle(   1.0080000e+00 );
    amoebaVdwForce->addParticle( 4, classIndex,   1.3500000e-01,   8.3680000e-02,   9.1000000e-01 );

    system.addParticle(   1.0080000e+00 );
    amoebaVdwForce->addParticle( 4, classIndex,   1.3500000e-01,   8.3680000e-02,   9.1000000e-01 );

    // ParticleExclusions

    std::vector< int > exclusions;
    exclusions.resize(0);
    exclusions.push_back( 0 );
    exclusions.push_back( 1 );
    exclusions.push_back( 2 );
    exclusions.push_back( 3 );
    amoebaVdwForce->setParticleExclusions( 0, exclusions );

    exclusions.resize(0);
    exclusions.push_back( 1 );
    exclusions.push_back( 0 );
    exclusions.push_back( 2 );
    exclusions.push_back( 3 );
    amoebaVdwForce->setParticleExclusions( 1, exclusions );

    exclusions.resize(0);
    exclusions.push_back( 2 );
    exclusions.push_back( 0 );
    exclusions.push_back( 1 );
    exclusions.push_back( 3 );
    amoebaVdwForce->setParticleExclusions( 2, exclusions );

    exclusions.resize(0);
    exclusions.push_back( 3 );
    exclusions.push_back( 0 );
    exclusions.push_back( 1 );
    exclusions.push_back( 2 );
    amoebaVdwForce->setParticleExclusions( 3, exclusions );

    exclusions.resize(0);
    exclusions.push_back( 4 );
    exclusions.push_back( 5 );
    exclusions.push_back( 6 );
    exclusions.push_back( 7 );
    amoebaVdwForce->setParticleExclusions( 4, exclusions );

    exclusions.resize(0);
    exclusions.push_back( 5 );
    exclusions.push_back( 4 );
    exclusions.push_back( 6 );
    exclusions.push_back( 7 );
    amoebaVdwForce->setParticleExclusions( 5, exclusions );

    exclusions.resize(0);
    exclusions.push_back( 6 );
    exclusions.push_back( 4 );
    exclusions.push_back( 5 );
    exclusions.push_back( 7 );
    amoebaVdwForce->setParticleExclusions( 6, exclusions );

    exclusions.resize(0);
    exclusions.push_back( 7 );
    exclusions.push_back( 4 );
    exclusions.push_back( 5 );
    exclusions.push_back( 6 );
    amoebaVdwForce->setParticleExclusions( 7, exclusions );

    // end of Vdw setup

    std::vector<Vec3> positions(numberOfParticles);

    positions[0]              = Vec3(   1.5927280e-01,   1.7000000e-06,    1.6491000e-03 );
    positions[1]              = Vec3(   2.0805540e-01,  -8.1258800e-02,    3.7282500e-02 );
    positions[2]              = Vec3(   2.0843610e-01,   8.0953200e-02,    3.7462200e-02 );
    positions[3]              = Vec3(   1.7280780e-01,   2.0730000e-04,   -9.8741700e-02 );
    positions[4]              = Vec3(  -1.6743680e-01,   1.5900000e-05,   -6.6149000e-03 );
    positions[5]              = Vec3(  -2.0428260e-01,   8.1071500e-02,    4.1343900e-02 );
    positions[6]              = Vec3(  -6.7308300e-02,   1.2800000e-05,    1.0623300e-02 );
    positions[7]              = Vec3(  -2.0426290e-01,  -8.1231400e-02,    4.1033500e-02 );

    system.addForce(amoebaVdwForce);

    std::string platformName;
    platformName = "Cuda";
    LangevinIntegrator integrator(0.0, 0.1, 0.01);
    Context context(system, integrator, Platform::getPlatformByName( platformName ) );

    context.setPositions(positions);
    State state                      = context.getState(State::Forces | State::Energy);
    forces                           = state.getForces();
    energy                           = state.getPotentialEnergy();
}

void compareForcesEnergy( std::string& testName, double expectedEnergy, double energy,
                          std::vector<Vec3>& expectedForces,
                          std::vector<Vec3>& forces, double tolerance, FILE* log ) {


#ifdef AMOEBA_DEBUG
    if( log ){
        (void) fprintf( log, "%s: expected energy=%14.7e %14.7e\n", testName.c_str(), expectedEnergy, state.getPotentialEnergy() );
        for( unsigned int ii = 0; ii < forces.size(); ii++ ){
            (void) fprintf( log, "%6u [%14.7e %14.7e %14.7e]   [%14.7e %14.7e %14.7e]\n", ii,
                            expectedForces[ii][0], expectedForces[ii][1], expectedForces[ii][2], forces[ii][0], forces[ii][1], forces[ii][2] );
        }
        (void) fflush( log );
    }
#endif

    for( unsigned int ii = 0; ii < forces.size(); ii++ ){
        ASSERT_EQUAL_VEC_MOD( expectedForces[ii], forces[ii], tolerance, testName );
    }
    ASSERT_EQUAL_TOL_MOD( expectedEnergy, energy, tolerance, testName );
}

// test VDW w/ sigmaRule=CubicMean and epsilonRule=HHG

void testVdwAmmoniaCubicMeanHhg( FILE* log ) {

    std::string testName      = "testVdwAmmoniaCubicMeanHhg";

    int numberOfParticles     = 8;
    double boxDimension       = -1.0;
    double cutoff             = 9000000.0;
    std::vector<Vec3> forces;
    double energy;

    setupAndGetForcesEnergyVdwAmmonia( "CUBIC-MEAN", "HHG", cutoff, boxDimension, forces, energy, log );
    std::vector<Vec3> expectedForces(numberOfParticles);

    double expectedEnergy     =  4.8012258e+00;

    expectedForces[0]         = Vec3(   2.9265247e+02,  -1.4507808e-02,  -6.9562123e+00 );
    expectedForces[1]         = Vec3(  -2.2451693e+00,   4.8143073e-01,  -2.0041494e-01 );
    expectedForces[2]         = Vec3(  -2.2440698e+00,  -4.7905450e-01,  -2.0125284e-01 );
    expectedForces[3]         = Vec3(  -1.0840394e+00,  -5.8531253e-04,   2.6934135e-01 );
    expectedForces[4]         = Vec3(  -5.6305662e+01,   1.4733908e-03,  -1.8083306e-01 );
    expectedForces[5]         = Vec3(   1.6750145e+00,  -3.2448374e-01,  -1.8030914e-01 );
    expectedForces[6]         = Vec3(  -2.3412420e+02,   1.0754069e-02,   7.6287492e+00 );
    expectedForces[7]         = Vec3(   1.6756544e+00,   3.2497316e-01,  -1.7906832e-01 );

    double tolerance          = 1.0e-04;
    compareForcesEnergy( testName, expectedEnergy, energy, expectedForces, forces, tolerance, log );
}

// test VDW w/ sigmaRule=Arithmetic and epsilonRule=Arithmetic

void testVdwAmmoniaArithmeticArithmetic( FILE* log ) {

    std::string testName      = "testVdwAmmoniaArithmeticArithmetic";

    int numberOfParticles     = 8;
    double boxDimension       = -1.0;
    double cutoff             = 9000000.0;

    std::vector<Vec3> forces;
    double energy;
    setupAndGetForcesEnergyVdwAmmonia( "ARITHMETIC", "ARITHMETIC", cutoff, boxDimension, forces, energy, log );
    std::vector<Vec3> expectedForces(numberOfParticles);

    double expectedEnergy     =  4.2252403e+00;

    expectedForces[0]         = Vec3(   3.0603839e+02,  -1.5550310e-02,  -7.2661707e+00 );
    expectedForces[1]         = Vec3(  -2.7801357e+00,   5.8805051e-01,  -2.5907269e-01 );
    expectedForces[2]         = Vec3(  -2.7753968e+00,  -5.8440732e-01,  -2.5969111e-01 );
    expectedForces[3]         = Vec3(  -2.2496416e+00,  -1.1797440e-03,   5.5501757e-01 );
    expectedForces[4]         = Vec3(  -5.5077629e+01,   8.3417114e-04,  -3.3668921e-01 );
    expectedForces[5]         = Vec3(   2.3752452e+00,  -4.6788669e-01,  -2.4907764e-01 );
    expectedForces[6]         = Vec3(  -2.4790697e+02,   1.1419770e-02,   8.0629999e+00 );
    expectedForces[7]         = Vec3(   2.3761408e+00,   4.6871961e-01,  -2.4731607e-01 );

    double tolerance          = 1.0e-04;
    compareForcesEnergy( testName, expectedEnergy, energy, expectedForces, forces, tolerance, log );
}

// test VDW w/ sigmaRule=Geometric and epsilonRule=Geometric

void testVdwAmmoniaGeometricGeometric( FILE* log ) {

    std::string testName      = "testVdwAmmoniaGeometricGeometric";

    int numberOfParticles     = 8;
    double boxDimension       = -1.0;
    double cutoff             = 9000000.0;
    std::vector<Vec3> forces;
    double energy;
    setupAndGetForcesEnergyVdwAmmonia( "GEOMETRIC", "GEOMETRIC", cutoff, boxDimension, forces, energy, log );
    std::vector<Vec3> expectedForces(numberOfParticles);

    double expectedEnergy     =  2.5249914e+00;

    expectedForces[0]         = Vec3(   2.1169631e+02,  -1.0710925e-02,  -4.3728025e+00 );
    expectedForces[1]         = Vec3(  -2.2585621e+00,   4.8409995e-01,  -2.0188344e-01 );
    expectedForces[2]         = Vec3(  -2.2551351e+00,  -4.8124855e-01,  -2.0246986e-01 );
    expectedForces[3]         = Vec3(  -1.7178028e+00,  -9.0851787e-04,   4.2466975e-01 );
    expectedForces[4]         = Vec3(  -4.8302147e+01,   9.6603376e-04,  -5.7972068e-01 );
    expectedForces[5]         = Vec3(   1.8100634e+00,  -3.5214093e-01,  -1.9357207e-01 );
    expectedForces[6]         = Vec3(  -1.6078365e+02,   7.2117601e-03,   5.3180261e+00 );
    expectedForces[7]         = Vec3(   1.8109211e+00,   3.5273117e-01,  -1.9224723e-01 );

    double tolerance          = 1.0e-04;
    compareForcesEnergy( testName, expectedEnergy, energy, expectedForces, forces, tolerance, log );
}

void testVdwAmmoniaCubicMeanHarmonic( FILE* log ) {

    std::string testName      = "testVdwAmmoniaCubicMeanHarmonic";

    int numberOfParticles     = 8;
    double boxDimension       = -1.0;
    double cutoff             = 9000000.0;
    std::vector<Vec3> forces;
    double energy;
    setupAndGetForcesEnergyVdwAmmonia( "CUBIC-MEAN", "HARMONIC", cutoff, boxDimension, forces, energy, log );
    std::vector<Vec3> expectedForces(numberOfParticles);

    double expectedEnergy     =  4.1369069e+00;

    expectedForces[0]         = Vec3(   2.5854436e+02,  -1.2779529e-02,  -5.9041148e+00 );
    expectedForces[1]         = Vec3(  -2.0832419e+00,   4.4915831e-01,  -1.8266000e-01 );
    expectedForces[2]         = Vec3(  -2.0823991e+00,  -4.4699804e-01,  -1.8347141e-01 );
    expectedForces[3]         = Vec3(  -9.5914714e-01,  -5.2162026e-04,   2.3873165e-01 );
    expectedForces[4]         = Vec3(  -5.3724787e+01,   1.4838241e-03,  -2.8089191e-01 );
    expectedForces[5]         = Vec3(   1.5074325e+00,  -2.9016397e-01,  -1.6385118e-01 );
    expectedForces[6]         = Vec3(  -2.0271029e+02,   9.2367947e-03,   6.6389988e+00 );
    expectedForces[7]         = Vec3(   1.5080748e+00,   2.9058422e-01,  -1.6274118e-01 );

    double tolerance          = 1.0e-04;
    compareForcesEnergy( testName, expectedEnergy, energy, expectedForces, forces, tolerance, log );
}

// test w/ cutoff=0.25 nm; single ixn between two particles (0 and 6); force nonzero on
// particle 4 due to reduction applied to NH
// the distance between 0 and 6 is ~ 0.235 so the ixn is in the tapered region

void testVdwTaper( FILE* log ) {

    std::string testName      = "testVdwTaper";

    int numberOfParticles     = 8;
    double boxDimension       = -1.0;
    double cutoff             = 0.25;

    std::vector<Vec3> forces;
    double energy;
    setupAndGetForcesEnergyVdwAmmonia( "CUBIC-MEAN", "HHG", cutoff, boxDimension, forces, energy, log );
    std::vector<Vec3> expectedForces(numberOfParticles);

    double expectedEnergy     =  3.5478444e+00;

    expectedForces[0]         = Vec3(   5.6710779e+02,  -2.7391004e-02,  -1.7867730e+01 );
    expectedForces[1]         = Vec3(  -0.0000000e+00,  -0.0000000e+00,  -0.0000000e+00 );
    expectedForces[2]         = Vec3(  -0.0000000e+00,  -0.0000000e+00,  -0.0000000e+00 );
    expectedForces[3]         = Vec3(  -0.0000000e+00,  -0.0000000e+00,  -0.0000000e+00 );
    expectedForces[4]         = Vec3(  -5.1039701e+01,   2.4651903e-03,   1.6080957e+00 );
    expectedForces[5]         = Vec3(  -0.0000000e+00,  -0.0000000e+00,  -0.0000000e+00 );
    expectedForces[6]         = Vec3(  -5.1606809e+02,   2.4925813e-02,   1.6259634e+01 );
    expectedForces[7]         = Vec3(  -0.0000000e+00,  -0.0000000e+00,  -0.0000000e+00 );

    double tolerance          = 1.0e-04;
    compareForcesEnergy( testName, expectedEnergy, energy, expectedForces, forces, tolerance, log );
}

// test PBC

void testVdwPBC( FILE* log ) {

    std::string testName      = "testVdwPBC";

    int numberOfParticles     = 8;
    double boxDimension       = 0.6;
    double cutoff             = 0.25;

    std::vector<Vec3> forces;
    double energy;
    setupAndGetForcesEnergyVdwAmmonia( "CUBIC-MEAN", "HHG", cutoff, boxDimension, forces, energy, log );
    std::vector<Vec3> expectedForces(numberOfParticles);

    double expectedEnergy     =  8.4385405e+00;

    expectedForces[0]         = Vec3(   5.1453069e+02,   4.9751912e-01,  -1.2759570e+01 );
    expectedForces[1]         = Vec3(  -2.5622586e+02,  -4.6524265e+01,   2.4281465e+01 );
    expectedForces[2]         = Vec3(  -2.7538705e+02,   5.1831690e+01,   2.7367710e+01 );
    expectedForces[3]         = Vec3(  -0.0000000e+00,  -0.0000000e+00,  -0.0000000e+00 );
    expectedForces[4]         = Vec3(   3.0883034e+02,  -5.8876974e+00,  -5.8286122e+01 );
    expectedForces[5]         = Vec3(   1.1319359e+02,  -3.2047069e-01,   1.6181231e+00 );
    expectedForces[6]         = Vec3(  -5.1606809e+02,   2.4925813e-02,   1.6259634e+01 );
    expectedForces[7]         = Vec3(   1.1112638e+02,   3.7829857e-01,   1.5187587e+00 );

    // tolerance is higher here due to interpolation used in setting tapering coefficients;
    // if tapering turned off, then absolute difference < 2.0e-05

    double tolerance          = 5.0e-04;
    compareForcesEnergy( testName, expectedEnergy, energy, expectedForces, forces, tolerance, log );
}

int main( int numberOfArguments, char* argv[] ) {

    try {
        std::cout << "TestCudaAmoebaVdwForce running test..." << std::endl;
        registerAmoebaCudaKernelFactories();

        FILE* log = NULL;

        testVdw( log );

        // tests using two ammonia molecules

        // test VDW w/ sigmaRule=CubicMean and epsilonRule=HHG

        testVdwAmmoniaCubicMeanHhg( log );

        // test VDW w/ sigmaRule=Arithmetic and epsilonRule=Arithmetic

        testVdwAmmoniaArithmeticArithmetic( log );

        // test VDW w/ sigmaRule=Geometric and epsilonRule=Geometric

        testVdwAmmoniaGeometricGeometric( log );

        // test VDW w/ sigmaRule=CubicMean and epsilonRule=Harmonic

        testVdwAmmoniaCubicMeanHarmonic( log );

        // test w/ cutoff=0.25 nm; single ixn between two particles (0 and 6); force nonzero on
        // particle 4 due to reduction applied to NH
        // the distance between 0 and 6 is ~ 0.235 so the ixn is in the tapered region

        testVdwTaper( log );

        // test PBC

        testVdwPBC( log );

    } catch(const std::exception& e) {
        std::cout << "exception: " << e.what() << std::endl;
        std::cout << "FAIL - ERROR.  Test failed." << std::endl;
        return 1;
    }
    std::cout << "Done" << std::endl;
    return 0;
}
