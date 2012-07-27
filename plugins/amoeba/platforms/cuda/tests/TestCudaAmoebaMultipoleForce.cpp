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
 * the rights to use, copy, modify, merge, publish, distribute, sublicense,  *
 * and/or sell copies of the Software, and to permit persons to whom the      *
 * Software is furnished to do so, subject to the following conditions:       *
 *                                                                            *
 * The above copyright notice and this permission notice shall be included in *
 * all copies or substantial portions of the Software.                        *
 *                                                                            *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR *
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,  *
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL    *
 * THE AUTHORS, CONTRIBUTORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,   *
 * DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR      *
 * OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE  *
 * USE OR OTHER DEALINGS IN THE SOFTWARE.                                     *
 * -------------------------------------------------------------------------- */

/**
 * This tests the Cuda implementation of CudaAmoebaMultipoleForce.
 */

#include "openmm/internal/AssertionUtilities.h"
#include "openmm/Context.h"
#include "OpenMMAmoeba.h"
#include "AmoebaTinkerParameterFile.h"
#include "openmm/System.h"
#include "openmm/AmoebaMultipoleForce.h"
#include "openmm/LangevinIntegrator.h"
#include <iostream>
#include <vector>
#include <stdlib.h>
#include <stdio.h>

#define ASSERT_EQUAL_TOL_MOD(expected, found, tol, testname) {double _scale_ = std::abs(expected) > 1.0 ? std::abs(expected) : 1.0; if (!(std::abs((expected)-(found))/_scale_ <= (tol))) {std::stringstream details; details << testname << " Expected "<<(expected)<<", found "<<(found); throwException(__FILE__, __LINE__, details.str());}};

#define ASSERT_EQUAL_VEC_MOD(expected, found, tol,testname) {ASSERT_EQUAL_TOL_MOD((expected)[0], (found)[0], (tol),(testname)); ASSERT_EQUAL_TOL_MOD((expected)[1], (found)[1], (tol),(testname)); ASSERT_EQUAL_TOL_MOD((expected)[2], (found)[2], (tol),(testname));};


using namespace OpenMM;
const double TOL = 1e-4;

void setupAndGetForcesEnergyMultipoleAmmonia( AmoebaMultipoleForce::AmoebaNonbondedMethod nonbondedMethod, AmoebaMultipoleForce::AmoebaPolarizationType polarizationType,
                                              double cutoff, int inputPmeGridDimension, std::vector<Vec3>& forces, double& energy, FILE* log ){

    // beginning of Multipole setup

    System system;

    double boxDimension                               = 0.6;
    Vec3 a( boxDimension, 0.0, 0.0 );
    Vec3 b( 0.0, boxDimension, 0.0 );
    Vec3 c( 0.0, 0.0, boxDimension );
    system.setDefaultPeriodicBoxVectors( a, b, c );

    AmoebaMultipoleForce* amoebaMultipoleForce        = new AmoebaMultipoleForce();;
    int numberOfParticles                             = 8;
    amoebaMultipoleForce->setNonbondedMethod( nonbondedMethod );
    amoebaMultipoleForce->setPolarizationType( polarizationType );
    amoebaMultipoleForce->setCutoffDistance( cutoff );
    amoebaMultipoleForce->setMutualInducedTargetEpsilon( 1.0e-06 );
    amoebaMultipoleForce->setMutualInducedMaxIterations( 500 );
    amoebaMultipoleForce->setAEwald( 1.4024714e+01 );
    amoebaMultipoleForce->setEwaldErrorTolerance( 1.0e-04 );

    std::vector<int> pmeGridDimension( 3 );
    pmeGridDimension[0] = pmeGridDimension[1] = pmeGridDimension[2] = inputPmeGridDimension;
    amoebaMultipoleForce->setPmeGridDimensions( pmeGridDimension );

    std::vector<double> nitrogenMolecularDipole(3);
    std::vector<double> nitrogenMolecularQuadrupole(9);

    nitrogenMolecularDipole[0]     =   8.3832254e-03;
    nitrogenMolecularDipole[1]     =   0.0000000e+00;
    nitrogenMolecularDipole[2]     =   3.4232474e-03;

    nitrogenMolecularQuadrupole[0] =  -4.0406249e-04;
    nitrogenMolecularQuadrupole[1] =   0.0000000e+00;
    nitrogenMolecularQuadrupole[2] =  -2.6883671e-04;
    nitrogenMolecularQuadrupole[3] =   0.0000000e+00;
    nitrogenMolecularQuadrupole[4] =   2.5463927e-04;
    nitrogenMolecularQuadrupole[5] =   0.0000000e+00;
    nitrogenMolecularQuadrupole[6] =  -2.6883671e-04;
    nitrogenMolecularQuadrupole[7] =   0.0000000e+00;
    nitrogenMolecularQuadrupole[8] =   1.4942322e-04;

    // first N

    system.addParticle( 1.4007000e+01 );
    amoebaMultipoleForce->addParticle(  -5.7960000e-01, nitrogenMolecularDipole, nitrogenMolecularQuadrupole, 2, 1, 2, 3,  3.9000000e-01,  3.1996314e-01,  1.0730000e-03 );

    // 3 H attached to first N

    std::vector<double> hydrogenMolecularDipole(3);
    std::vector<double> hydrogenMolecularQuadrupole(9);
    hydrogenMolecularDipole[0]     =  -1.7388763e-03;
    hydrogenMolecularDipole[1]     =   0.0000000e+00;
    hydrogenMolecularDipole[2]     =  -4.6837475e-03;

    hydrogenMolecularQuadrupole[0] =  -4.4253841e-05;
    hydrogenMolecularQuadrupole[1] =   0.0000000e+00;
    hydrogenMolecularQuadrupole[2] =   1.5429571e-05;
    hydrogenMolecularQuadrupole[3] =   0.0000000e+00;
    hydrogenMolecularQuadrupole[4] =   4.1798924e-05;
    hydrogenMolecularQuadrupole[5] =   0.0000000e+00;
    hydrogenMolecularQuadrupole[6] =   1.5429571e-05;
    hydrogenMolecularQuadrupole[7] =   0.0000000e+00;
    hydrogenMolecularQuadrupole[8] =   2.4549167e-06;

    system.addParticle( 1.0080000e+00 );
    system.addParticle( 1.0080000e+00 );
    system.addParticle( 1.0080000e+00 );
    amoebaMultipoleForce->addParticle(   1.9320000e-01, hydrogenMolecularDipole, hydrogenMolecularQuadrupole, 2, 0, 2, 3, 3.9000000e-01,  2.8135002e-01,  4.9600000e-04 );
    amoebaMultipoleForce->addParticle(   1.9320000e-01, hydrogenMolecularDipole, hydrogenMolecularQuadrupole, 2, 0, 1, 3, 3.9000000e-01,  2.8135002e-01,  4.9600000e-04 );
    amoebaMultipoleForce->addParticle(   1.9320000e-01, hydrogenMolecularDipole, hydrogenMolecularQuadrupole, 2, 0, 1, 2, 3.9000000e-01,  2.8135002e-01,  4.9600000e-04 );

    // second N

    system.addParticle(   1.4007000e+01 );
    amoebaMultipoleForce->addParticle(  -5.7960000e-01, nitrogenMolecularDipole, nitrogenMolecularQuadrupole, 2, 5, 6, 7,  3.9000000e-01,  3.1996314e-01,  1.0730000e-03 );

    // 3 H attached to second N

    system.addParticle(   1.0080000e+00 );
    system.addParticle(   1.0080000e+00 );
    system.addParticle(   1.0080000e+00 );
    amoebaMultipoleForce->addParticle(   1.9320000e-01, hydrogenMolecularDipole, hydrogenMolecularQuadrupole, 2, 4, 6, 7, 3.9000000e-01,  2.8135002e-01,  4.9600000e-04 );
    amoebaMultipoleForce->addParticle(   1.9320000e-01, hydrogenMolecularDipole, hydrogenMolecularQuadrupole, 2, 4, 5, 7, 3.9000000e-01,  2.8135002e-01,  4.9600000e-04 );
    amoebaMultipoleForce->addParticle(   1.9320000e-01, hydrogenMolecularDipole, hydrogenMolecularQuadrupole, 2, 4, 5, 6, 3.9000000e-01,  2.8135002e-01,  4.9600000e-04 );

    // CovalentMaps

    std::vector< int > covalentMap;
    covalentMap.resize(0);
    covalentMap.push_back( 1 );
    covalentMap.push_back( 2 );
    covalentMap.push_back( 3 );
    amoebaMultipoleForce->setCovalentMap( 0, static_cast<OpenMM::AmoebaMultipoleForce::CovalentType>(0), covalentMap );

    covalentMap.resize(0);
    covalentMap.push_back( 0 );
    covalentMap.push_back( 1 );
    covalentMap.push_back( 2 );
    covalentMap.push_back( 3 );
    amoebaMultipoleForce->setCovalentMap( 0, static_cast<OpenMM::AmoebaMultipoleForce::CovalentType>(4), covalentMap );

    covalentMap.resize(0);
    covalentMap.push_back( 0 );
    amoebaMultipoleForce->setCovalentMap( 1, static_cast<OpenMM::AmoebaMultipoleForce::CovalentType>(0), covalentMap );

    covalentMap.resize(0);
    covalentMap.push_back( 2 );
    covalentMap.push_back( 3 );
    amoebaMultipoleForce->setCovalentMap( 1, static_cast<OpenMM::AmoebaMultipoleForce::CovalentType>(1), covalentMap );

    covalentMap.resize(0);
    covalentMap.push_back( 0 );
    covalentMap.push_back( 1 );
    covalentMap.push_back( 2 );
    covalentMap.push_back( 3 );
    amoebaMultipoleForce->setCovalentMap( 1, static_cast<OpenMM::AmoebaMultipoleForce::CovalentType>(4), covalentMap );

    covalentMap.resize(0);
    covalentMap.push_back( 0 );
    amoebaMultipoleForce->setCovalentMap( 2, static_cast<OpenMM::AmoebaMultipoleForce::CovalentType>(0), covalentMap );

    covalentMap.resize(0);
    covalentMap.push_back( 1 );
    covalentMap.push_back( 3 );
    amoebaMultipoleForce->setCovalentMap( 2, static_cast<OpenMM::AmoebaMultipoleForce::CovalentType>(1), covalentMap );

    covalentMap.resize(0);
    covalentMap.push_back( 0 );
    covalentMap.push_back( 1 );
    covalentMap.push_back( 2 );
    covalentMap.push_back( 3 );
    amoebaMultipoleForce->setCovalentMap( 2, static_cast<OpenMM::AmoebaMultipoleForce::CovalentType>(4), covalentMap );

    covalentMap.resize(0);
    covalentMap.push_back( 0 );
    amoebaMultipoleForce->setCovalentMap( 3, static_cast<OpenMM::AmoebaMultipoleForce::CovalentType>(0), covalentMap );

    covalentMap.resize(0);
    covalentMap.push_back( 1 );
    covalentMap.push_back( 2 );
    amoebaMultipoleForce->setCovalentMap( 3, static_cast<OpenMM::AmoebaMultipoleForce::CovalentType>(1), covalentMap );

    covalentMap.resize(0);
    covalentMap.push_back( 0 );
    covalentMap.push_back( 1 );
    covalentMap.push_back( 2 );
    covalentMap.push_back( 3 );
    amoebaMultipoleForce->setCovalentMap( 3, static_cast<OpenMM::AmoebaMultipoleForce::CovalentType>(4), covalentMap );

    covalentMap.resize(0);
    covalentMap.push_back( 5 );
    covalentMap.push_back( 6 );
    covalentMap.push_back( 7 );
    amoebaMultipoleForce->setCovalentMap( 4, static_cast<OpenMM::AmoebaMultipoleForce::CovalentType>(0), covalentMap );

    covalentMap.resize(0);
    covalentMap.push_back( 4 );
    covalentMap.push_back( 5 );
    covalentMap.push_back( 6 );
    covalentMap.push_back( 7 );
    amoebaMultipoleForce->setCovalentMap( 4, static_cast<OpenMM::AmoebaMultipoleForce::CovalentType>(4), covalentMap );

    covalentMap.resize(0);
    covalentMap.push_back( 4 );
    amoebaMultipoleForce->setCovalentMap( 5, static_cast<OpenMM::AmoebaMultipoleForce::CovalentType>(0), covalentMap );

    covalentMap.resize(0);
    covalentMap.push_back( 6 );
    covalentMap.push_back( 7 );
    amoebaMultipoleForce->setCovalentMap( 5, static_cast<OpenMM::AmoebaMultipoleForce::CovalentType>(1), covalentMap );

    covalentMap.resize(0);
    covalentMap.push_back( 4 );
    covalentMap.push_back( 5 );
    covalentMap.push_back( 6 );
    covalentMap.push_back( 7 );
    amoebaMultipoleForce->setCovalentMap( 5, static_cast<OpenMM::AmoebaMultipoleForce::CovalentType>(4), covalentMap );

    covalentMap.resize(0);
    covalentMap.push_back( 4 );
    amoebaMultipoleForce->setCovalentMap( 6, static_cast<OpenMM::AmoebaMultipoleForce::CovalentType>(0), covalentMap );

    covalentMap.resize(0);
    covalentMap.push_back( 5 );
    covalentMap.push_back( 7 );
    amoebaMultipoleForce->setCovalentMap( 6, static_cast<OpenMM::AmoebaMultipoleForce::CovalentType>(1), covalentMap );

    covalentMap.resize(0);
    covalentMap.push_back( 4 );
    covalentMap.push_back( 5 );
    covalentMap.push_back( 6 );
    covalentMap.push_back( 7 );
    amoebaMultipoleForce->setCovalentMap( 6, static_cast<OpenMM::AmoebaMultipoleForce::CovalentType>(4), covalentMap );

    covalentMap.resize(0);
    covalentMap.push_back( 4 );
    amoebaMultipoleForce->setCovalentMap( 7, static_cast<OpenMM::AmoebaMultipoleForce::CovalentType>(0), covalentMap );

    covalentMap.resize(0);
    covalentMap.push_back( 5 );
    covalentMap.push_back( 6 );
    amoebaMultipoleForce->setCovalentMap( 7, static_cast<OpenMM::AmoebaMultipoleForce::CovalentType>(1), covalentMap );

    covalentMap.resize(0);
    covalentMap.push_back( 4 );
    covalentMap.push_back( 5 );
    covalentMap.push_back( 6 );
    covalentMap.push_back( 7 );
    amoebaMultipoleForce->setCovalentMap( 7, static_cast<OpenMM::AmoebaMultipoleForce::CovalentType>(4), covalentMap );

    AmoebaHarmonicBondForce* amoebaHarmonicBondForce  = new AmoebaHarmonicBondForce();
    int numberOfBonds                                 = 6; 

    // addParticle: particle1, particle2, length, quadraticK

    amoebaHarmonicBondForce->addBond( 0, 1,   0.0000000e+00,   0.0000000e+00 );
    amoebaHarmonicBondForce->addBond( 0, 2,   0.0000000e+00,   0.0000000e+00 );
    amoebaHarmonicBondForce->addBond( 0, 3,   0.0000000e+00,   0.0000000e+00 );

    amoebaHarmonicBondForce->addBond( 4, 5,   0.0000000e+00,   0.0000000e+00 );
    amoebaHarmonicBondForce->addBond( 4, 6,   0.0000000e+00,   0.0000000e+00 );
    amoebaHarmonicBondForce->addBond( 4, 7,   0.0000000e+00,   0.0000000e+00 );
    amoebaHarmonicBondForce->setAmoebaGlobalHarmonicBondCubic( -2.5500000e+01 ); 
    amoebaHarmonicBondForce->setAmoebaGlobalHarmonicBondQuartic( 3.7931250e+02 ); 
    system.addForce(amoebaHarmonicBondForce);

    std::vector<Vec3> positions(numberOfParticles);

    positions[0]              = Vec3(   1.5927280e-01,  1.7000000e-06,   1.6491000e-03 );
    positions[1]              = Vec3(   2.0805540e-01, -8.1258800e-02,   3.7282500e-02 );
    positions[2]              = Vec3(   2.0843610e-01,  8.0953200e-02,   3.7462200e-02 );
    positions[3]              = Vec3(   1.7280780e-01,  2.0730000e-04,  -9.8741700e-02 );
    positions[4]              = Vec3(  -1.6743680e-01,  1.5900000e-05,  -6.6149000e-03 );
    positions[5]              = Vec3(  -2.0428260e-01,  8.1071500e-02,   4.1343900e-02 );
    positions[6]              = Vec3(  -6.7308300e-02,  1.2800000e-05,   1.0623300e-02 );
    positions[7]              = Vec3(  -2.0426290e-01, -8.1231400e-02,   4.1033500e-02 );

    system.addForce(amoebaMultipoleForce);

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


//#define AMOEBA_DEBUG
#ifdef AMOEBA_DEBUG
    if( log ){
        double conversion = 1.0/4.184;
        (void) fprintf( log, "%s: expected energy=%14.7e %14.7e\n", testName.c_str(), conversion*expectedEnergy, conversion*energy );
        if( conversion != 1.0 )conversion *= -0.1;
        for( unsigned int ii = 0; ii < forces.size(); ii++ ){
            (void) fprintf( log, "%6u [%14.7e %14.7e %14.7e]   [%14.7e %14.7e %14.7e]\n", ii,
                            conversion*expectedForces[ii][0], conversion*expectedForces[ii][1], conversion*expectedForces[ii][2],
                            conversion*forces[ii][0], conversion*forces[ii][1], conversion*forces[ii][2] );
        }
        (void) fflush( log );
        conversion = 1.0;
        (void) fprintf( log, "\n%s: expected energy=%14.7e %14.7e no conversion\n", testName.c_str(), conversion*expectedEnergy, conversion*energy );
        if( conversion != 1.0 )conversion = -1.0;
        for( unsigned int ii = 0; ii < forces.size(); ii++ ){
            (void) fprintf( log, "%6u [%14.7e %14.7e %14.7e]   [%14.7e %14.7e %14.7e]\n", ii,
                            conversion*expectedForces[ii][0], conversion*expectedForces[ii][1], conversion*expectedForces[ii][2],
                            conversion*forces[ii][0], conversion*forces[ii][1], conversion*forces[ii][2] );
        }
        (void) fflush( log );
    }
#endif

    for( unsigned int ii = 0; ii < forces.size(); ii++ ){
        ASSERT_EQUAL_VEC_MOD( expectedForces[ii], forces[ii], tolerance, testName );
    }
    ASSERT_EQUAL_TOL_MOD( expectedEnergy, energy, tolerance, testName );
}

// test multipole direct polarization

void testMultipoleAmmoniaDirectPolarization( FILE* log ) {

    std::string testName      = "testMultipoleAmmoniaDirectPolarization";

    int numberOfParticles     = 8;
    int inputPmeGridDimension = 0;
    double cutoff             = 9000000.0;
    std::vector<Vec3> forces;
    double energy;

    setupAndGetForcesEnergyMultipoleAmmonia( AmoebaMultipoleForce::NoCutoff, AmoebaMultipoleForce::Direct, 
                                             cutoff, inputPmeGridDimension, forces, energy, log );
    std::vector<Vec3> expectedForces(numberOfParticles);

    double expectedEnergy     = -1.7428832e+01;

    expectedForces[0]         = Vec3(  -3.5574000e+02, -7.3919340e+00,  3.8989934e+01 );
    expectedForces[1]         = Vec3(   3.0368045e+01, -8.7325694e+00,  6.9731151e+00 );
    expectedForces[2]         = Vec3(   3.2358980e+01,  1.0234924e+01,  4.7203694e-01 );
    expectedForces[3]         = Vec3(   2.1439022e+01,  5.8998414e+00, -3.8355239e+01 );
    expectedForces[4]         = Vec3(  -1.8052760e+02, -1.0618455e+00, -7.0030146e+01 );
    expectedForces[5]         = Vec3(   4.2411304e+01, -1.6569222e+01,  1.9047581e+00 );
    expectedForces[6]         = Vec3(   3.6823677e+02,  7.7839986e-01,  5.8404590e+01 );
    expectedForces[7]         = Vec3(   4.1453480e+01,  1.6842405e+01,  1.6409513e+00 );

    double tolerance          = 1.0e-04;
    compareForcesEnergy( testName, expectedEnergy, energy, expectedForces, forces, tolerance, log );
}

// test multipole mutual polarization

void testMultipoleAmmoniaMutualPolarization( FILE* log ) {

    std::string testName      = "testMultipoleAmmoniaMutualPolarization";

    int numberOfParticles     = 8;
    int inputPmeGridDimension = 0;
    double cutoff             = 9000000.0;
    std::vector<Vec3> forces;
    double energy;

    setupAndGetForcesEnergyMultipoleAmmonia( AmoebaMultipoleForce::NoCutoff, AmoebaMultipoleForce::Mutual, 
                                             cutoff, inputPmeGridDimension, forces, energy, log );
    std::vector<Vec3> expectedForces(numberOfParticles);

    double expectedEnergy     = -1.7790449e+01;

    expectedForces[0]         = Vec3(  -3.7523158e+02,  -7.9806295e+00,   3.7464051e+01 );
    expectedForces[1]         = Vec3(   3.1352410e+01,  -9.4055551e+00,   8.5230415e+00 );
    expectedForces[2]         = Vec3(   3.3504923e+01,   1.1029935e+01,   1.5052263e+00 );
    expectedForces[3]         = Vec3(   2.3295507e+01,   6.3698827e+00,  -4.0403553e+01 );
    expectedForces[4]         = Vec3(  -1.9379275e+02,  -1.0903937e+00,  -7.3461740e+01 );
    expectedForces[5]         = Vec3(   4.3278067e+01,  -1.6906589e+01,   1.5721909e+00 );
    expectedForces[6]         = Vec3(   3.9529983e+02,   7.9661172e-01,   6.3499055e+01 );
    expectedForces[7]         = Vec3(   4.2293601e+01,   1.7186738e+01,   1.3017270e+00 );

    double tolerance          = 1.0e-04;
    compareForcesEnergy( testName, expectedEnergy, energy, expectedForces, forces, tolerance, log );
}

// test multipole mutual polarization using PME

void testMultipoleAmmoniaPMEMutualPolarization( FILE* log ) {

    std::string testName      = "testMultipoleAmmoniaMutualPolarization";

    int numberOfParticles     = 8;
    int inputPmeGridDimension = 80;
    double cutoff             = 0.28;
    std::vector<Vec3> forces;
    double energy;

    setupAndGetForcesEnergyMultipoleAmmonia( AmoebaMultipoleForce::PME, AmoebaMultipoleForce::Mutual, 
                                             cutoff, inputPmeGridDimension, forces, energy, log );
    std::vector<Vec3> expectedForces(numberOfParticles);

    double expectedEnergy     = -3.8843855e+00*4.184;

    expectedForces[0]         = Vec3(  -5.5121003e+02,  -1.7325162e+01,   8.1781763e+01 );
    expectedForces[1]         = Vec3(   2.5644285e+01,   3.4249929e+01,  -1.6040819e+01 );
    expectedForces[2]         = Vec3(   2.9751974e+01,  -3.0851062e+01,  -3.1741707e+01 );
    expectedForces[3]         = Vec3(   1.6147883e+02,   1.3744637e+01,  -2.6982931e+01 );
    expectedForces[4]         = Vec3(  -3.1289925e+02,  -3.3127430e+00,  -5.9526840e+01 );
    expectedForces[5]         = Vec3(   1.3046796e+02,  -3.8528176e+01,   4.9662804e+00 );
    expectedForces[6]         = Vec3(   4.6654349e+02,   2.5146526e+00,   7.4042527e+01 );
    expectedForces[7]         = Vec3(   1.2682825e+02,   3.9078836e+01,   3.8640742e+00 );

    double tolerance          = 1.0e-04;
    compareForcesEnergy( testName, expectedEnergy, energy, expectedForces, forces, tolerance, log );
}

void setupAndGetForcesEnergyMultipoleWater( AmoebaMultipoleForce::AmoebaNonbondedMethod nonbondedMethod, AmoebaMultipoleForce::AmoebaPolarizationType polarizationType,
                                            double cutoff, int inputPmeGridDimension, std::vector<Vec3>& forces, double& energy, FILE* log ){

    // beginning of Multipole setup

    System system;

    double boxDimension                               = 1.8643;
    Vec3 a( boxDimension, 0.0, 0.0 );
    Vec3 b( 0.0, boxDimension, 0.0 );
    Vec3 c( 0.0, 0.0, boxDimension );
    system.setDefaultPeriodicBoxVectors( a, b, c );

    AmoebaMultipoleForce* amoebaMultipoleForce        = new AmoebaMultipoleForce();;
    int numberOfParticles                             = 12;
    amoebaMultipoleForce->setNonbondedMethod( nonbondedMethod );
    amoebaMultipoleForce->setPolarizationType( polarizationType );
    amoebaMultipoleForce->setCutoffDistance( cutoff );
    amoebaMultipoleForce->setMutualInducedTargetEpsilon( 1.0e-06 );
    amoebaMultipoleForce->setMutualInducedMaxIterations( 500 );
    amoebaMultipoleForce->setAEwald( 5.4459052e+00 );
    amoebaMultipoleForce->setEwaldErrorTolerance( 1.0e-04 );

    std::vector<int> pmeGridDimension( 3 );
    pmeGridDimension[0] = pmeGridDimension[1] = pmeGridDimension[2] = inputPmeGridDimension;
    amoebaMultipoleForce->setPmeGridDimensions( pmeGridDimension );

    for( unsigned int jj = 0; jj < numberOfParticles; jj += 3 ){
        system.addParticle( 1.5995000e+01 );
        system.addParticle( 1.0080000e+00 );
        system.addParticle( 1.0080000e+00 );
    }

    std::vector<double> oxygenMolecularDipole(3);
    std::vector<double> oxygenMolecularQuadrupole(9);

    oxygenMolecularDipole[0]     =   0.0000000e+00;
    oxygenMolecularDipole[1]     =   0.0000000e+00;
    oxygenMolecularDipole[2]     =   7.5561214e-03;

    oxygenMolecularQuadrupole[0] =   3.5403072e-04;
    oxygenMolecularQuadrupole[1] =   0.0000000e+00;
    oxygenMolecularQuadrupole[2] =   0.0000000e+00;
    oxygenMolecularQuadrupole[3] =   0.0000000e+00;
    oxygenMolecularQuadrupole[4] =  -3.9025708e-04;
    oxygenMolecularQuadrupole[5] =   0.0000000e+00;
    oxygenMolecularQuadrupole[6] =   0.0000000e+00;
    oxygenMolecularQuadrupole[7] =   0.0000000e+00;
    oxygenMolecularQuadrupole[8] =   3.6226356e-05;

    std::vector<double> hydrogenMolecularDipole(3);
    std::vector<double> hydrogenMolecularQuadrupole(9);
    hydrogenMolecularDipole[0]     =  -2.0420949e-03;
    hydrogenMolecularDipole[1]     =   0.0000000e+00;
    hydrogenMolecularDipole[2]     =  -3.0787530e-03;

    hydrogenMolecularQuadrupole[0] =  -3.4284825e-05;
    hydrogenMolecularQuadrupole[1] =   0.0000000e+00;
    hydrogenMolecularQuadrupole[2] =  -1.8948597e-06;
    hydrogenMolecularQuadrupole[3] =   0.0000000e+00;
    hydrogenMolecularQuadrupole[4] =  -1.0024088e-04;
    hydrogenMolecularQuadrupole[5] =   0.0000000e+00;
    hydrogenMolecularQuadrupole[6] =  -1.8948597e-06;
    hydrogenMolecularQuadrupole[7] =   0.0000000e+00;
    hydrogenMolecularQuadrupole[8] =   1.3452570e-04;

    for( unsigned int jj = 0; jj < numberOfParticles; jj += 3 ){
        amoebaMultipoleForce->addParticle( -5.1966000e-01, oxygenMolecularDipole, oxygenMolecularQuadrupole, 1, jj+1, jj+2, -1,
                                            3.9000000e-01, 3.0698765e-01, 8.3700000e-04 );
        amoebaMultipoleForce->addParticle(  2.5983000e-01, hydrogenMolecularDipole, hydrogenMolecularQuadrupole, 0, jj, jj+2, -1, 
                                            3.9000000e-01, 2.8135002e-01, 4.9600000e-04 );
        amoebaMultipoleForce->addParticle(  2.5983000e-01, hydrogenMolecularDipole, hydrogenMolecularQuadrupole, 0, jj, jj+1, -1, 
                                            3.9000000e-01, 2.8135002e-01, 4.9600000e-04 );
    }

    // CovalentMaps

    std::vector< int > covalentMap;
    for( unsigned int jj = 0; jj < numberOfParticles; jj += 3 ){
        covalentMap.resize(0);
        covalentMap.push_back( jj+1 );
        covalentMap.push_back( jj+2 );
        amoebaMultipoleForce->setCovalentMap( jj, static_cast<OpenMM::AmoebaMultipoleForce::CovalentType>(0), covalentMap );

        covalentMap.resize(0);
        covalentMap.push_back( jj );
        covalentMap.push_back( jj+1 );
        covalentMap.push_back( jj+2 );
        amoebaMultipoleForce->setCovalentMap( jj,   static_cast<OpenMM::AmoebaMultipoleForce::CovalentType>(4), covalentMap );
        amoebaMultipoleForce->setCovalentMap( jj+1, static_cast<OpenMM::AmoebaMultipoleForce::CovalentType>(4), covalentMap );
        amoebaMultipoleForce->setCovalentMap( jj+2, static_cast<OpenMM::AmoebaMultipoleForce::CovalentType>(4), covalentMap );
    
        covalentMap.resize(0);
        covalentMap.push_back( jj );
        amoebaMultipoleForce->setCovalentMap( jj+1, static_cast<OpenMM::AmoebaMultipoleForce::CovalentType>(0), covalentMap );
        amoebaMultipoleForce->setCovalentMap( jj+2, static_cast<OpenMM::AmoebaMultipoleForce::CovalentType>(0), covalentMap );
    
        covalentMap.resize(0);
        covalentMap.push_back( jj+2 );
        amoebaMultipoleForce->setCovalentMap( jj+1, static_cast<OpenMM::AmoebaMultipoleForce::CovalentType>(1), covalentMap );
    
        covalentMap.resize(0);
        covalentMap.push_back( jj+1 );
        amoebaMultipoleForce->setCovalentMap( jj+2, static_cast<OpenMM::AmoebaMultipoleForce::CovalentType>(1), covalentMap );
    
    } 
 
    AmoebaHarmonicBondForce* amoebaHarmonicBondForce  = new AmoebaHarmonicBondForce();

    // addParticle: particle1, particle2, length, quadraticK

    for( unsigned int jj = 0; jj < numberOfParticles; jj += 3 ){
        amoebaHarmonicBondForce->addBond( jj, jj+1,   0.0000000e+00,   0.0000000e+00 );
        amoebaHarmonicBondForce->addBond( jj, jj+2,   0.0000000e+00,   0.0000000e+00 );
    }

    amoebaHarmonicBondForce->setAmoebaGlobalHarmonicBondCubic( -2.5500000e+01 ); 
    amoebaHarmonicBondForce->setAmoebaGlobalHarmonicBondQuartic( 3.7931250e+02 ); 
    system.addForce(amoebaHarmonicBondForce);

    std::vector<Vec3> positions(numberOfParticles);

    positions[0]              = Vec3(  -8.7387270e-01,   5.3220410e-01,    7.4214000e-03 );
    positions[1]              = Vec3(  -9.6050090e-01,   5.1173410e-01,   -2.2202700e-02 );
    positions[2]              = Vec3(  -8.5985900e-01,   4.9658230e-01,    1.0283390e-01 );
    positions[3]              = Vec3(   9.1767100e-02,  -7.8956650e-01,    4.3804200e-01 );
    positions[4]              = Vec3(   1.2333420e-01,  -7.0267430e-01,    4.2611550e-01 );
    positions[5]              = Vec3(   1.7267090e-01,  -8.2320810e-01,    4.8124750e-01 );
    positions[6]              = Vec3(   8.6290110e-01,   6.2153500e-02,    4.1280850e-01 );
    positions[7]              = Vec3(   8.6385200e-01,   1.2684730e-01,    3.3887060e-01 );
    positions[8]              = Vec3(   9.5063550e-01,   5.3173300e-02,    4.4799160e-01 );
    positions[9]              = Vec3(   5.0844930e-01,   2.8684740e-01,   -6.9293750e-01 );
    positions[10]             = Vec3(   6.0459330e-01,   3.0620510e-01,   -7.0100130e-01 );
    positions[11]             = Vec3(   5.0590640e-01,   1.8880920e-01,   -6.8813470e-01 );

    system.addForce(amoebaMultipoleForce);

    std::string platformName;
    platformName = "Cuda";
    LangevinIntegrator integrator(0.0, 0.1, 0.01);
    Context context(system, integrator, Platform::getPlatformByName( platformName ) );

    context.setPositions(positions);
    State state                      = context.getState(State::Forces | State::Energy);
    forces                           = state.getForces();
    energy                           = state.getPotentialEnergy();
}

// test multipole mutual polarization using PME

void testMultipoleWaterPMEDirectPolarization( FILE* log ) {

    std::string testName      = "testMultipoleWaterDirectPolarization";

    int numberOfParticles     = 12;
    int inputPmeGridDimension = 20;
    double cutoff             = 0.70;
    std::vector<Vec3> forces;
    double energy;

    setupAndGetForcesEnergyMultipoleWater( AmoebaMultipoleForce::PME, AmoebaMultipoleForce::Direct, 
                                            cutoff, inputPmeGridDimension, forces, energy, log );
    std::vector<Vec3> expectedForces(numberOfParticles);

    double expectedEnergy     = 6.4585115e-01;

    expectedForces[0]         = Vec3(  -1.2396731e+00,  -2.4231698e+01,   8.3348523e+00 );
    expectedForces[1]         = Vec3(  -3.3737276e+00,   9.9304523e+00,  -6.3917827e+00 );
    expectedForces[2]         = Vec3(   4.4062247e+00,   1.9518971e+01,  -4.6552873e+00 );
    expectedForces[3]         = Vec3(  -1.3128824e+00,  -1.2887339e+00,  -1.4473147e+00 );
    expectedForces[4]         = Vec3(   2.1137034e+00,   3.9457973e-01,   2.9269129e-01 );
    expectedForces[5]         = Vec3(   1.0271174e+00,   1.2039367e+00,   1.2112214e+00 );
    expectedForces[6]         = Vec3(  -3.2082903e+00,   1.4979371e+01,  -1.0274832e+00 );
    expectedForces[7]         = Vec3(  -1.1880320e+00,  -1.5177166e+01,   2.5525509e+00 );
    expectedForces[8]         = Vec3(   4.3607105e+00,  -7.0253274e+00,   2.9522580e-01 );
    expectedForces[9]         = Vec3(  -3.0175134e+00,   1.3607102e+00,   6.6883370e+00 );
    expectedForces[10]        = Vec3(   9.2036949e-01,  -1.4717629e+00,  -3.3362339e+00 );
    expectedForces[11]        = Vec3(   1.2523841e+00,  -1.9794292e+00,  -3.4670129e+00 );

    double tolerance          = 1.0e-03;
    compareForcesEnergy( testName, expectedEnergy, energy, expectedForces, forces, tolerance, log );
}

// test multipole mutual polarization using PME

void testMultipoleWaterPMEMutualPolarization( FILE* log ) {

    std::string testName      = "testMultipoleWaterMutualPolarization";

    int numberOfParticles     = 12;
    int inputPmeGridDimension = 20;
    double cutoff             = 0.70;
    std::vector<Vec3> forces;
    double energy;

    setupAndGetForcesEnergyMultipoleWater( AmoebaMultipoleForce::PME, AmoebaMultipoleForce::Mutual, 
                                            cutoff, inputPmeGridDimension, forces, energy, log );
    std::vector<Vec3> expectedForces(numberOfParticles);

    double expectedEnergy     =  6.5029855e-01;

    expectedForces[0]         = Vec3(  -1.2367386e+00,  -2.4197036e+01,   8.3256759e+00 );
    expectedForces[1]         = Vec3(  -3.3825187e+00,   9.9387618e+00,  -6.4200475e+00 );
    expectedForces[2]         = Vec3(   4.4108644e+00,   1.9486127e+01,  -4.6530661e+00 );
    expectedForces[3]         = Vec3(  -1.3129168e+00,  -1.2947383e+00,  -1.4438198e+00 );
    expectedForces[4]         = Vec3(   2.1144837e+00,   3.9590305e-01,   2.9040889e-01 );
    expectedForces[5]         = Vec3(   1.0287222e+00,   1.2100201e+00,   1.2103068e+00 );
    expectedForces[6]         = Vec3(  -3.2017550e+00,   1.4995985e+01,  -1.1036504e+00 );
    expectedForces[7]         = Vec3(  -1.2065398e+00,  -1.5192899e+01,   2.6233368e+00 );
    expectedForces[8]         = Vec3(   4.3698604e+00,  -7.0550315e+00,   3.4204565e-01 );
    expectedForces[9]         = Vec3(  -3.0082825e+00,   1.3575082e+00,   6.6901032e+00 );
    expectedForces[10]        = Vec3(   9.1775539e-01,  -1.4651882e+00,  -3.3322516e+00 );
    expectedForces[11]        = Vec3(   1.2467701e+00,  -1.9832979e+00,  -3.4684052e+00 );

    double tolerance          = 1.0e-03;
    compareForcesEnergy( testName, expectedEnergy, energy, expectedForces, forces, tolerance, log );
}

int main( int numberOfArguments, char* argv[] ) {

    try {
        std::cout << "TestCudaAmoebaMultipoleForce running test..." << std::endl;
        registerAmoebaCudaKernelFactories();

        FILE* log = NULL;
        //FILE* log = stderr;

        // tests using two ammonia molecules

        // test direct polarization, no cutoff

        testMultipoleAmmoniaDirectPolarization( log );

        // test mutual polarization, no cutoff

        testMultipoleAmmoniaMutualPolarization( log );

        // test multipole direct & mutual polarization using PME

        //testMultipoleAmmoniaPMEMutualPolarization( log );

        testMultipoleWaterPMEDirectPolarization( log );
        testMultipoleWaterPMEMutualPolarization( log );

    } catch(const std::exception& e) {
        std::cout << "exception: " << e.what() << std::endl;
        std::cout << "FAIL - ERROR.  Test failed." << std::endl;
        return 1;
    }
    std::cout << "Done" << std::endl;
    return 0;
}
