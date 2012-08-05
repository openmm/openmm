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

// setup for 2 ammonia molecules

static void setupAndGetForcesEnergyMultipoleAmmonia( AmoebaMultipoleForce::AmoebaPolarizationType polarizationType,
                                                     int includeCavityTerm, std::vector<Vec3>& forces, double& energy, FILE* log ){

    // beginning of Multipole setup

    System system;


    AmoebaMultipoleForce* amoebaMultipoleForce        = new AmoebaMultipoleForce();;
    int numberOfParticles                             = 8;

    amoebaMultipoleForce->setNonbondedMethod( AmoebaMultipoleForce::NoCutoff );
    amoebaMultipoleForce->setPolarizationType( polarizationType );
    amoebaMultipoleForce->setMutualInducedTargetEpsilon( 1.0e-06 );
    amoebaMultipoleForce->setMutualInducedMaxIterations( 500 );

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

    // covalent maps

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
    system.addForce(amoebaMultipoleForce);

    // GK force

    AmoebaGeneralizedKirkwoodForce* amoebaGeneralizedKirkwoodForce  = new AmoebaGeneralizedKirkwoodForce();
    amoebaGeneralizedKirkwoodForce->setSolventDielectric(   7.8300000e+01 );
    amoebaGeneralizedKirkwoodForce->setSoluteDielectric(    1.0000000e+00 );
    amoebaGeneralizedKirkwoodForce->setIncludeCavityTerm( includeCavityTerm );

    // addParticle: charge, radius, scalingFactor

    for( unsigned int ii = 0; ii < 2; ii++ ){
        amoebaGeneralizedKirkwoodForce->addParticle(  -5.7960000e-01,   1.5965000e-01,   6.9000000e-01 );
        amoebaGeneralizedKirkwoodForce->addParticle(   1.9320000e-01,   1.2360000e-01,   6.9000000e-01 );
        amoebaGeneralizedKirkwoodForce->addParticle(   1.9320000e-01,   1.2360000e-01,   6.9000000e-01 );
        amoebaGeneralizedKirkwoodForce->addParticle(   1.9320000e-01,   1.2360000e-01,   6.9000000e-01 );
    }
    system.addForce(amoebaGeneralizedKirkwoodForce);

    // 1-2 bonds needed 
/*
    AmoebaHarmonicBondForce* amoebaHarmonicBondForce  = new AmoebaHarmonicBondForce();

    // addBond: particle1, particle2, length, quadraticK

    amoebaHarmonicBondForce->addBond( 0, 1,   0.0000000e+00,   0.0000000e+00 );
    amoebaHarmonicBondForce->addBond( 0, 2,   0.0000000e+00,   0.0000000e+00 );
    amoebaHarmonicBondForce->addBond( 0, 3,   0.0000000e+00,   0.0000000e+00 );

    amoebaHarmonicBondForce->addBond( 4, 5,   0.0000000e+00,   0.0000000e+00 );
    amoebaHarmonicBondForce->addBond( 4, 6,   0.0000000e+00,   0.0000000e+00 );
    amoebaHarmonicBondForce->addBond( 4, 7,   0.0000000e+00,   0.0000000e+00 );
    amoebaHarmonicBondForce->setAmoebaGlobalHarmonicBondCubic( -2.5500000e+01 ); 
    amoebaHarmonicBondForce->setAmoebaGlobalHarmonicBondQuartic( 3.7931250e+02 ); 
    system.addForce(amoebaHarmonicBondForce);
*/
    std::vector<Vec3> positions(numberOfParticles);

    positions[0]              = Vec3(   1.5927280e-01,  1.7000000e-06,   1.6491000e-03 );
    positions[1]              = Vec3(   2.0805540e-01, -8.1258800e-02,   3.7282500e-02 );
    positions[2]              = Vec3(   2.0843610e-01,  8.0953200e-02,   3.7462200e-02 );
    positions[3]              = Vec3(   1.7280780e-01,  2.0730000e-04,  -9.8741700e-02 );
    positions[4]              = Vec3(  -1.6743680e-01,  1.5900000e-05,  -6.6149000e-03 );
    positions[5]              = Vec3(  -2.0428260e-01,  8.1071500e-02,   4.1343900e-02 );
    positions[6]              = Vec3(  -6.7308300e-02,  1.2800000e-05,   1.0623300e-02 );
    positions[7]              = Vec3(  -2.0426290e-01, -8.1231400e-02,   4.1033500e-02 );

    std::string platformName;
    platformName = "Cuda";
    LangevinIntegrator integrator(0.0, 0.1, 0.01);
    Context context(system, integrator, Platform::getPlatformByName( platformName ) );

    context.setPositions(positions);
    State state                      = context.getState(State::Forces | State::Energy);
    forces                           = state.getForces();
    energy                           = state.getPotentialEnergy();
}

// compare forces and energies 

static void compareForcesEnergy( std::string& testName, double expectedEnergy, double energy,
                                 std::vector<Vec3>& expectedForces,
                                 std::vector<Vec3>& forces, double tolerance, FILE* log ) {


//#define AMOEBA_DEBUG
#ifdef AMOEBA_DEBUG
    if( log ){
        double conversion = 1.0/4.184;
        double energyAbsDiff = fabs( expectedEnergy - energy );   
        double energyRelDiff =  2.0*energyAbsDiff/( fabs( expectedEnergy ) + fabs( energy ) + 1.0e-08 );   
        (void) fprintf( log, "%s: expected energy=%14.7e %14.7e  absDiff=%15.7e relDiff=%15.7e\n", testName.c_str(), conversion*expectedEnergy, conversion*energy,
                        conversion*energyAbsDiff, conversion*energyRelDiff );
        if( conversion != 1.0 )conversion *= -0.1;
        for( unsigned int ii = 0; ii < forces.size(); ii++ ){

            double expectedNorm = sqrt( expectedForces[ii][0]*expectedForces[ii][0] +
                                        expectedForces[ii][1]*expectedForces[ii][1] +
                                        expectedForces[ii][2]*expectedForces[ii][2] );

            double norm         = sqrt( forces[ii][0]*forces[ii][0] + forces[ii][1]*forces[ii][1] + forces[ii][2]*forces[ii][2] );
            double absDiff      = fabs( norm - expectedNorm );
            double relDiff      = 2.0*absDiff/(fabs( norm ) + fabs( expectedNorm ) + 1.0e-08);

            (void) fprintf( log, "%6u %15.7e %15.7e [%14.7e %14.7e %14.7e]   [%14.7e %14.7e %14.7e]\n", ii,
                            conversion*absDiff, conversion*relDiff,
                            conversion*expectedForces[ii][0], conversion*expectedForces[ii][1], conversion*expectedForces[ii][2],
                            conversion*forces[ii][0], conversion*forces[ii][1], conversion*forces[ii][2], conversion*expectedNorm, conversion*norm );
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

// test GK direct polarization for system comprised of two ammonia molecules

static void testGeneralizedKirkwoodAmmoniaDirectPolarization( FILE* log ) {

    std::string testName      = "testGeneralizedKirkwoodAmmoniaDirectPolarization";

    int numberOfParticles     = 8;
    std::vector<Vec3> forces;
    double energy;

    setupAndGetForcesEnergyMultipoleAmmonia( AmoebaMultipoleForce::Direct, 0, forces, energy, log );
    std::vector<Vec3> expectedForces(numberOfParticles);

    double expectedEnergy     = -7.6636680e+01;

    expectedForces[0]         = Vec3(  -6.9252994e+02,  -8.9085133e+00,   9.6489739e+01 );
    expectedForces[1]         = Vec3(   1.5593797e+02,  -6.0331931e+01,   1.5104507e+01 );
    expectedForces[2]         = Vec3(   1.5870088e+02,   6.1702809e+01,   6.7708985e+00 );
    expectedForces[3]         = Vec3(   1.4089885e+02,   7.5870617e+00,  -1.1362294e+02 );
    expectedForces[4]         = Vec3(  -1.8916205e+02,   2.1465549e-01,  -4.3433152e+02 );
    expectedForces[5]         = Vec3(   1.0208290e+01,   6.2676753e+01,   1.4987953e+02 );
    expectedForces[6]         = Vec3(   4.0621859e+02,   1.8962203e-01,   1.3021956e+02 );
    expectedForces[7]         = Vec3(   9.7274235e+00,  -6.3130458e+01,   1.4949024e+02 );

    double tolerance          = 1.0e-04;
    compareForcesEnergy( testName, expectedEnergy, energy, expectedForces, forces, tolerance, log );
}

// test GK mutual polarization for system comprised of two ammonia molecules

static void testGeneralizedKirkwoodAmmoniaMutualPolarization( FILE* log ) {

    std::string testName      = "testGeneralizedKirkwoodAmmoniaMutualPolarization";

    int numberOfParticles     = 8;
    std::vector<Vec3> forces;
    double energy;

    setupAndGetForcesEnergyMultipoleAmmonia( AmoebaMultipoleForce::Mutual, 0, forces, energy, log );
    std::vector<Vec3> expectedForces(numberOfParticles);

    double expectedEnergy     =  -7.8018875e+01;

    expectedForces[0]         = Vec3(  -7.6820301e+02,  -1.0102760e+01,   1.0094389e+02 );
    expectedForces[1]         = Vec3(   1.7037307e+02,  -7.5621857e+01,   2.3320365e+01 );
    expectedForces[2]         = Vec3(   1.7353828e+02,   7.7199741e+01,   1.3965379e+01 );
    expectedForces[3]         = Vec3(   1.5045244e+02,   8.5784569e+00,  -1.3377619e+02 );
    expectedForces[4]         = Vec3(  -2.1811615e+02,  -1.6818022e-01,  -4.6103163e+02 );
    expectedForces[5]         = Vec3(   6.2091942e+00,   7.6748687e+01,   1.5883463e+02 );
    expectedForces[6]         = Vec3(   4.8035662e+02,   4.9704902e-01,   1.3948083e+02 );
    expectedForces[7]         = Vec3(   5.3895456e+00,  -7.7131137e+01,   1.5826273e+02 );

    double tolerance          = 1.0e-04;
    compareForcesEnergy( testName, expectedEnergy, energy, expectedForces, forces, tolerance, log );
}

// test GK mutual polarization for system comprised of two ammonia molecules

static void testGeneralizedKirkwoodAmmoniaMutualPolarizationWithCavityTerm( FILE* log ) {

    std::string testName      = "testGeneralizedKirkwoodAmmoniaMutualPolarizationWithCavityTerm";

    int numberOfParticles     = 8;
    std::vector<Vec3> forces;
    double energy;

    setupAndGetForcesEnergyMultipoleAmmonia( AmoebaMultipoleForce::Mutual, 1, forces, energy, log );
    std::vector<Vec3> expectedForces(numberOfParticles);

    double expectedEnergy     = -6.0434582e+01;

    expectedForces[0]         = Vec3(  -7.8323218e+02,  -1.0097644e+01,   1.0256890e+02 );
    expectedForces[1]         = Vec3(   1.7078480e+02,  -7.1896701e+01,   2.0840172e+01 );
    expectedForces[2]         = Vec3(   1.7394089e+02,   7.3488594e+01,   1.1484648e+01 );
    expectedForces[3]         = Vec3(   1.5169364e+02,   8.5611299e+00,  -1.2968050e+02 );
    expectedForces[4]         = Vec3(  -2.1669693e+02,  -1.5926823e-01,  -4.6636274e+02 );
    expectedForces[5]         = Vec3(   8.7397444e+00,   7.3330990e+01,   1.6016898e+02 );
    expectedForces[6]         = Vec3(   4.8684950e+02,   4.8937161e-01,   1.4137061e+02 );
    expectedForces[7]         = Vec3(   7.9205382e+00,  -7.3716473e+01,   1.5960993e+02 );

    double tolerance          = 1.0e-04;
    compareForcesEnergy( testName, expectedEnergy, energy, expectedForces, forces, tolerance, log );
}

int main( int numberOfArguments, char* argv[] ) {

    try {
        std::cout << "TestCudaAmoebaMultipoleForce running test..." << std::endl;
        registerAmoebaCudaKernelFactories();

        FILE* log = NULL;

        // test direct and mutual polarization cases and
        // mutual polarization w/ the cavity term

        testGeneralizedKirkwoodAmmoniaDirectPolarization( log );
        testGeneralizedKirkwoodAmmoniaMutualPolarization( log );
        testGeneralizedKirkwoodAmmoniaMutualPolarizationWithCavityTerm( log );

    } catch(const std::exception& e) {
        std::cout << "exception: " << e.what() << std::endl;
        std::cout << "FAIL - ERROR.  Test failed." << std::endl;
        return 1;
    }
    std::cout << "Done" << std::endl;
    return 0;
}
