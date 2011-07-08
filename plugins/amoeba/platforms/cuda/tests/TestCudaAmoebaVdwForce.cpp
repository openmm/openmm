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

#include "../../../tests/AssertionUtilities.h"
#include "openmm/Context.h"
#include "OpenMMAmoeba.h"
#include "AmoebaTinkerParameterFile.h"
#include "openmm/System.h"
#include "openmm/AmoebaVdwForce.h"
#include "openmm/LangevinIntegrator.h"
#include <iostream>
#include <vector>

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

    for( int ii = 0; ii < numberOfParticles; ii++ ){
        int indexIV, indexClass;
        double mass, sigma, epsilon, reduction;
        std::vector< int > exclusions;
        if( ii == 0 || ii == 3 ){
            mass        = 16.0;
            indexIV     = ii;
            indexClass  = 70;
            sigma       = 1.70250E+00;
            epsilon     = 1.10000E-01;
            reduction   = 0.0;
        } else {
            mass        = 1.0;
            indexIV     = ii < 3 ? 0 : 3;
            indexClass  = 71;
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
        amoebaVdwForce->addParticle( indexIV, indexClass, sigma, epsilon, reduction );
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
        int indexIV, indexClass;
        double sigma, epsilon, reduction;
        amoebaVdwForce->getParticleParameters( ii, indexIV, indexClass, sigma, epsilon, reduction );
        sigma        *= AngstromToNm;
        epsilon      *= CalToJoule;
        amoebaVdwForce->setParticleParameters( ii, indexIV, indexClass, sigma, epsilon, reduction );
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


int main( int numberOfArguments, char* argv[] ) {

    try {
        std::cout << "TestCudaAmoebaVdwForce running test..." << std::endl;
        registerAmoebaCudaKernelFactories();

        FILE* log = NULL;
        testVdw( log );
    } catch(const std::exception& e) {
        std::cout << "exception: " << e.what() << std::endl;
        std::cout << "FAIL - ERROR.  Test failed." << std::endl;
        return 1;
    }
    std::cout << "Done" << std::endl;
    return 0;
}
