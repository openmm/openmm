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
 * This tests the Cuda implementation of CudaAmoebaTorsionForce.
 */

#include "../../../tests/AssertionUtilities.h"
#include "AmoebaTinkerParameterFile.h"
#include "openmm/Context.h"
#include "OpenMMAmoeba.h"
#include "openmm/System.h"
#include "openmm/LangevinIntegrator.h"
#include <iostream>
#include <vector>

using namespace OpenMM;

const double TOL = 1e-3;
#define PI_M               3.141592653589
#define RADIAN            57.29577951308

/* ---------------------------------------------------------------------------------------

   Compute cross product of two 3-vectors and place in 3rd vector

   vectorZ = vectorX x vectorY

   @param vectorX             x-vector
   @param vectorY             y-vector
   @param vectorZ             z-vector

   @return vector is vectorZ

   --------------------------------------------------------------------------------------- */
     
static void crossProductVector3( double* vectorX, double* vectorY, double* vectorZ ){

    vectorZ[0]  = vectorX[1]*vectorY[2] - vectorX[2]*vectorY[1];
    vectorZ[1]  = vectorX[2]*vectorY[0] - vectorX[0]*vectorY[2];
    vectorZ[2]  = vectorX[0]*vectorY[1] - vectorX[1]*vectorY[0];

    return;
}

static double dotVector3( double* vectorX, double* vectorY ){
    return vectorX[0]*vectorY[0] + vectorX[1]*vectorY[1] + vectorX[2]*vectorY[2];
}

static void computeAmoebaTorsionForce(int bondIndex,  std::vector<Vec3>& positions, AmoebaTorsionForce& amoebaTorsionForce,
                                      std::vector<Vec3>& forces, double* energy, FILE* log ) {

    int particle1, particle2, particle3, particle4;

    std::vector<double> torsion1;
    std::vector<double> torsion2;
    std::vector<double> torsion3;

    torsion1.resize(3);
    torsion2.resize(3);
    torsion3.resize(3);

    amoebaTorsionForce.getTorsionParameters(bondIndex, particle1, particle2, particle3, particle4, torsion1, torsion2, torsion3);

    std::vector< std::vector<double> > torsions;
    torsions.push_back( torsion1 );
    torsions.push_back( torsion2 );
    torsions.push_back( torsion3 );

    if( log ){
        (void) fprintf( log, "computeAmoebaTorsionForce: bond %d [%d %d %d %d]\n", 
                             bondIndex, particle1, particle2, particle3, particle4 );
        for( unsigned int ii = 0; ii < 3; ii++ ){
            (void) fprintf( log, "    [%10.3e %10.3e %10.3e]\n", torsions[ii][0], torsions[ii][1], torsions[ii][2] ); 
        }
        (void) fflush( log );
    }

    enum { BA, CB, DC, CA, DB, LastDeltaIndex };
    double deltaR[LastDeltaIndex][3];
    for( int ii = 0; ii < 3; ii++ ){
        deltaR[BA][ii] = positions[particle2][ii] - positions[particle1][ii];
        deltaR[CB][ii] = positions[particle3][ii] - positions[particle2][ii];
        deltaR[DC][ii] = positions[particle4][ii] - positions[particle3][ii];
        deltaR[CA][ii] = positions[particle3][ii] - positions[particle1][ii];
        deltaR[DB][ii] = positions[particle4][ii] - positions[particle2][ii];
    }   

    enum { Xt, Xu, Xtu, LastXtIndex };
    double crossProducts[LastXtIndex][3];
    crossProductVector3( deltaR[BA], deltaR[CB], crossProducts[Xt] );
    crossProductVector3( deltaR[CB], deltaR[DC], crossProducts[Xu] );
    crossProductVector3( crossProducts[Xt], crossProducts[Xu], crossProducts[Xtu] );
 
    double rT2   = dotVector3( crossProducts[Xt], crossProducts[Xt] );
    double rU2   = dotVector3( crossProducts[Xu], crossProducts[Xu] );
    double rTrU  = sqrt( rT2*rU2 );
    if( rTrU <= 0.0 ){
        return;
    }
 
    double rCB   = dotVector3( deltaR[CB], deltaR[CB] );
         rCB     = sqrt( rCB );
 
    // ---------------------------------------------------------------------------------------
 
    // cos(w), cos(2w), cos(3w), ... 
    // sin(w), sin(2w), sin(3w), ... 
  
    double cosine[6], sine[6];
   
    cosine[0]  = dotVector3( crossProducts[Xt], crossProducts[Xu] );
    cosine[0] /= rTrU;
 
    sine[0]    = dotVector3( deltaR[CB], crossProducts[Xtu] );
    sine[0]   /= (rCB*rTrU);
 
    for( int ii = 1; ii < 3; ii++ ){
        cosine[ii] = cosine[0]*cosine[ii-1] - sine[0]*  sine[ii-1];
          sine[ii] = cosine[0]*  sine[ii-1] + sine[0]*cosine[ii-1];
    }
 
    // ---------------------------------------------------------------------------------------
 
    // dEdPhi prefactor
  
    double dEdPhi = 0.0;
    for( int ii = 0; ii < 3; ii++ ){
        dEdPhi += torsions[ii][0]*((double) (ii+1))*( cosine[ii]*sin( torsions[ii][1] ) - sine[ii]*cos( torsions[ii][1] ) );
    }
 
    // ---------------------------------------------------------------------------------------
  
    // dEdtu[0]      = dEdT
    // dEdtu[1]      = dEdU
  
    // tempVector[0] == dEdA: dEdT x CB
    // tempVector[1] == dEdB: (CA x dEdT) + (dEdU x DC)
    // tempVector[2] == dEdC: (dEdT x BA) + (DB x dEdU)
    // tempVector[3] == dEdD: (dEdU x CB)
   
    double dEdtu[2][3];
    double tempVector[6][3];
 
    // dEdT & dEdU
  
    crossProductVector3( crossProducts[Xt], deltaR[CB], tempVector[0] );
    crossProductVector3( crossProducts[Xu], deltaR[CB], tempVector[1] );
    double norm[2] = { dEdPhi/(rT2*rCB ), -dEdPhi/(rU2*rCB ) };
    for( int jj = 0; jj < 2; jj++ ){
        for( int ii = 0; ii < 3; ii++ ){
            dEdtu[jj][ii] = norm[jj]*tempVector[jj][ii];
        }
    }
 
    // dEdA
   
    crossProductVector3( dEdtu[0], deltaR[CB], tempVector[0] );
 
    // dEdB
  
    crossProductVector3( deltaR[CA], dEdtu[0], tempVector[4] );
    crossProductVector3( dEdtu[1], deltaR[DC], tempVector[1] );
 
    // dEdC
 
    crossProductVector3( dEdtu[0], deltaR[BA], tempVector[5] );
    crossProductVector3( deltaR[DB], dEdtu[1], tempVector[2] );
    for( int jj = 0; jj < 3; jj++ ){
        tempVector[1][jj] += tempVector[4][jj];
        tempVector[2][jj] += tempVector[5][jj];
    }
 
    // dEdD
  
    crossProductVector3( dEdtu[1], deltaR[CB], tempVector[3] );
 
    // ---------------------------------------------------------------------------------------
     
    // accumulate forces and energy
  
    forces[particle1][0]       -= tempVector[0][0];
    forces[particle1][1]       -= tempVector[0][1];
    forces[particle1][2]       -= tempVector[0][2];

    forces[particle2][0]       -= tempVector[1][0];
    forces[particle2][1]       -= tempVector[1][1];
    forces[particle2][2]       -= tempVector[1][2];

    forces[particle3][0]       -= tempVector[2][0];
    forces[particle3][1]       -= tempVector[2][1];
    forces[particle3][2]       -= tempVector[2][2];

    forces[particle4][0]       -= tempVector[3][0];
    forces[particle4][1]       -= tempVector[3][1];
    forces[particle4][2]       -= tempVector[3][2];

    double energyTerm = 0.0;
    for( int ii = 0; ii < 3; ii++ ){
        energyTerm += torsions[ii][0]*( 1.0 + cosine[ii]*cos( torsions[ii][1] ) + sine[ii]*sin( torsions[ii][1] ) );
    }
    *energy    += energyTerm;

} 

static void computeAmoebaTorsionForces( Context& context, AmoebaTorsionForce& amoebaTorsionForce,
                                        std::vector<Vec3>& expectedForces, double* expectedEnergy, FILE* log ) {

    // get positions and zero forces

    State state                 = context.getState(State::Positions);
    std::vector<Vec3> positions = state.getPositions();
    expectedForces.resize( positions.size() );
    
    for( unsigned int ii = 0; ii < expectedForces.size(); ii++ ){
        expectedForces[ii][0] = expectedForces[ii][1] = expectedForces[ii][2] = 0.0;
    }

    // calculates forces/energy

    *expectedEnergy = 0.0;
    for( int ii = 0; ii < amoebaTorsionForce.getNumTorsions(); ii++ ){
        computeAmoebaTorsionForce(ii, positions, amoebaTorsionForce, expectedForces, expectedEnergy, log );
    }

    if( log ){
        (void) fprintf( log, "computeAmoebaTorsionForces: expected energy=%14.7e\n", *expectedEnergy );
        for( unsigned int ii = 0; ii < positions.size(); ii++ ){
            (void) fprintf( log, "%6u [%14.7e %14.7e %14.7e]\n", ii, expectedForces[ii][0], expectedForces[ii][1], expectedForces[ii][2] );
        }
        (void) fflush( log );
    }
    return;

}

void compareWithExpectedForceAndEnergy( Context& context, AmoebaTorsionForce& amoebaTorsionForce,
                                        double tolerance, const std::string& idString, FILE* log) {

    std::vector<Vec3> expectedForces;
    double expectedEnergy;
    computeAmoebaTorsionForces( context, amoebaTorsionForce, expectedForces, &expectedEnergy, log );
   
    State state                      = context.getState(State::Forces | State::Energy);
    const std::vector<Vec3> forces   = state.getForces();

    if( log ){
        (void) fprintf( log, "computeAmoebaTorsionForces: expected energy=%14.7e %14.7e\n", expectedEnergy, state.getPotentialEnergy() );
        for( unsigned int ii = 0; ii < forces.size(); ii++ ){
            (void) fprintf( log, "%6u [%14.7e %14.7e %14.7e]   [%14.7e %14.7e %14.7e]\n", ii,
                            expectedForces[ii][0], expectedForces[ii][1], expectedForces[ii][2], forces[ii][0], forces[ii][1], forces[ii][2] );
        }
        (void) fflush( log );
    }

    for( unsigned int ii = 0; ii < forces.size(); ii++ ){
        ASSERT_EQUAL_VEC( expectedForces[ii], forces[ii], tolerance );
    }
    ASSERT_EQUAL_TOL( expectedEnergy, state.getPotentialEnergy(), tolerance );
}

void testOneTorsion( FILE* log ) {

    System system;
    int numberOfParticles = 4;
    for( int ii = 0; ii < numberOfParticles; ii++ ){
        system.addParticle(1.0);
    }

    LangevinIntegrator integrator(0.0, 0.1, 0.01);

    AmoebaTorsionForce* amoebaTorsionForce = new AmoebaTorsionForce();

    std::vector<double> torsion1;
    torsion1.push_back(  0.619500000E+00 );
    torsion1.push_back(  0.000000000E+00 );

    std::vector<double> torsion2;
    torsion2.push_back( -0.202500000E+00 );
    torsion2.push_back(  0.180000000E+03 );

    std::vector<double> torsion3;
    torsion3.push_back(  0.175000000E-01 );
    torsion3.push_back(  0.000000000E+00 );
    amoebaTorsionForce->addTorsion(0, 1, 2, 3, torsion1, torsion2, torsion3 );

    system.addForce(amoebaTorsionForce);
    Context context(system, integrator, Platform::getPlatformByName( "Cuda"));

    std::vector<Vec3> positions(numberOfParticles);

    positions[0] = Vec3( 0.278860000E+01,  0.264630000E+01,  0.426300000E+00 );
    positions[1] = Vec3( 0.273400000E+01,  0.244300000E+01,  0.261400000E+00 );
    positions[2] = Vec3( 0.262660000E+01,  0.254130000E+01,  0.284200000E+00 );
    positions[3] = Vec3( 0.269130000E+01,  0.266390000E+01,  0.353100000E+00 );

    context.setPositions(positions);
    compareWithExpectedForceAndEnergy( context, *amoebaTorsionForce, TOL, "testOneTorsion", log );

}

int main( int numberOfArguments, char* argv[] ) {

    try {
        std::cout << "TestCudaAmoebaTorsionForce running test..." << std::endl;
        Platform::loadPluginsFromDirectory( Platform::getDefaultPluginsDirectory() );
        //FILE* log = stderr;
        FILE* log = NULL;
        //FILE* log = fopen( "AmoebaTorsionForce.log", "w" );;
        testOneTorsion( log );
        if( log && log != stderr )
            (void) fclose( log );

    }
    catch(const std::exception& e) {
        std::cout << "exception: " << e.what() << std::endl;
        std::cout << "FAIL - ERROR.  Test failed." << std::endl;
        return 1;
    }
    std::cout << "PASS - Test succeeded." << std::endl;
    return 0;
}
