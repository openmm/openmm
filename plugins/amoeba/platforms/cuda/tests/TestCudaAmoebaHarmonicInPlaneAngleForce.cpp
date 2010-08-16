/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
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
 * This tests the Cuda implementation of CudaAmoebaHarmonicInPlaneAngleForce.
 */

#include "../../../tests/AssertionUtilities.h"
#include "AmoebaTinkerParameterFile.h"
#include "openmm/Context.h"
#include "AmoebaOpenMM.h"
#include "openmm/System.h"
#include "openmm/LangevinIntegrator.h"
#include <iostream>
#include <vector>

using namespace OpenMM;

const double TOL = 1e-5;
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

static void getPrefactorsGivenInPlaneAngleCosine( double cosine, double idealInPlaneAngle, double quadraticK, double cubicK,
                                                  double quarticK, double penticK, double sexticK,
                                                  double* dEdR, double* energyTerm, FILE* log ) {

    double angle;
    if( cosine >= 1.0 ){
        angle = 0.0f;
    } else if( cosine <= -1.0 ){
        angle = RADIAN*PI_M;
    } else {
        angle = RADIAN*acos(cosine);
    }

    if( log ){
        (void) fprintf( log, "getPrefactorsGivenInPlaneAngleCosine: cosine=%10.3e angle=%10.3e ideal=%10.3e\n", cosine, angle, idealInPlaneAngle ); 
        (void) fflush( log );
    }

    double deltaIdeal         = angle - idealInPlaneAngle;
    double deltaIdeal2        = deltaIdeal*deltaIdeal;
    double deltaIdeal3        = deltaIdeal*deltaIdeal2;
    double deltaIdeal4        = deltaIdeal2*deltaIdeal2;
 
    // deltaIdeal = r - r_0
 
    *dEdR        = ( 2.0                        +
                     3.0*cubicK*  deltaIdeal    +
                     4.0*quarticK*deltaIdeal2   +
                     5.0*penticK* deltaIdeal3   +
                     6.0*sexticK* deltaIdeal4     );
 
    *dEdR       *= RADIAN*quadraticK*deltaIdeal;
 

    *energyTerm  = 1.0f + cubicK* deltaIdeal    +
                          quarticK*deltaIdeal2   +
                          penticK* deltaIdeal3   +
                          sexticK* deltaIdeal4;
    *energyTerm *= quadraticK*deltaIdeal2;

    return;
}

static void computeAmoebaHarmonicInPlaneAngleForce(int bondIndex,  std::vector<Vec3>& positions, AmoebaHarmonicInPlaneAngleForce& amoebaHarmonicInPlaneAngleForce,
                                                   std::vector<Vec3>& forces, double* energy, FILE* log ) {

    int particle1, particle2, particle3, particle4;
    double idealInPlaneAngle;
    double quadraticK;
    amoebaHarmonicInPlaneAngleForce.getAngleParameters(bondIndex, particle1, particle2, particle3, particle4, idealInPlaneAngle, quadraticK );

    double cubicK         = amoebaHarmonicInPlaneAngleForce.getAmoebaGlobalHarmonicInPlaneAngleCubic();
    double quarticK       = amoebaHarmonicInPlaneAngleForce.getAmoebaGlobalHarmonicInPlaneAngleQuartic();
    double penticK        = amoebaHarmonicInPlaneAngleForce.getAmoebaGlobalHarmonicInPlaneAnglePentic();
    double sexticK        = amoebaHarmonicInPlaneAngleForce.getAmoebaGlobalHarmonicInPlaneAngleSextic();

    if( log ){
        (void) fprintf( log, "computeAmoebaHarmonicInPlaneAngleForce: bond %d [%d %d %d %d] ang=%10.3f k=%10.3f [%10.3e %10.3e %10.3e %10.3e]\n", 
                             bondIndex, particle1, particle2, particle3, particle4, idealInPlaneAngle, quadraticK, cubicK, quarticK, penticK, sexticK );
        (void) fflush( log );
    }

    // T   = AD x CD
    // P   = B + T*delta
    // AP  = A - P
    // CP  = A - P
    // M   = CP x AP

    enum { AD, BD, CD, T, AP, P, CP, M, APxM, CPxM, ADxBD, BDxCD, TxCD, ADxT, dBxAD, CDxdB, LastDeltaAtomIndex };
 
    // AD   0
    // BD   1
    // CD   2 
    //  T   3
    // AP   4
    //  P   5
    // CP   6
    // M    7
    // APxM, CPxM, ADxBD, BDxCD, TxCD, ADxT, dBxAD, CDxdB, LastDeltaAtomIndex

    double deltaR[LastDeltaAtomIndex][3];
    for( int ii = 0; ii < 3; ii++ ){
        deltaR[AD][ii] = positions[particle1][ii] - positions[particle4][ii];
        deltaR[BD][ii] = positions[particle2][ii] - positions[particle4][ii];
        deltaR[CD][ii] = positions[particle3][ii] - positions[particle4][ii];
    }
    crossProductVector3( deltaR[AD], deltaR[CD], deltaR[T] );
 
    double rT2     = dotVector3( deltaR[T], deltaR[T] );
    double delta   = dotVector3( deltaR[T], deltaR[BD] );
         delta    *= -1.0/rT2;
 
    for( int ii = 0; ii < 3; ii++ ){
       deltaR[P][ii]  = positions[particle2][ii] + deltaR[T][ii]*delta;
       deltaR[AP][ii] = positions[particle1][ii] - deltaR[P][ii];
       deltaR[CP][ii] = positions[particle3][ii] - deltaR[P][ii];
    }   
 
    double rAp2 = dotVector3( deltaR[AP],  deltaR[AP] );
    double rCp2 = dotVector3( deltaR[CP],  deltaR[CP] );
    if( rAp2 <= 0.0 && rCp2 <= 0.0 ){
        if( log ){
            (void) fprintf( log, "computeAmoebaHarmonicInPlaneAngleForce:  rAp2 or rCp2 <= 0.0\n" );
            (void) fflush( log );
        }
        return;
    }

    crossProductVector3( deltaR[CP], deltaR[AP], deltaR[M] );
 
    double rm = dotVector3( deltaR[M], deltaR[M] );
         rm   = sqrt( rm );
    if( rm < 0.000001 ){
       rm = 0.000001;
    }
 
    double dot     = dotVector3( deltaR[AP], deltaR[CP] );
    double cosine  = dot/sqrt( rAp2*rCp2 );
 
    double dEdR;
    double energyTerm;
    getPrefactorsGivenInPlaneAngleCosine( cosine, idealInPlaneAngle, quadraticK, cubicK,
                                          quarticK, penticK, sexticK, &dEdR,  &energyTerm, log );
 
    double termA   = -dEdR/(rAp2*rm);
    double termC   =  dEdR/(rCp2*rm);
 
    crossProductVector3( deltaR[AP], deltaR[M], deltaR[APxM] );
    crossProductVector3( deltaR[CP], deltaR[M], deltaR[CPxM] );
 
    // forces will be gathered here
 
    enum { dA, dB, dC, dD, LastDIndex };
    double forceTerm[LastDIndex][3];
 
    for( int ii = 0; ii < 3; ii++ ){
       forceTerm[dA][ii] = deltaR[APxM][ii]*termA;
       forceTerm[dC][ii] = deltaR[CPxM][ii]*termC;
       forceTerm[dB][ii] = -1.0*( forceTerm[dA][ii] + forceTerm[dC][ii] );
    }
 
    double pTrT2  = dotVector3( forceTerm[dB], deltaR[T] );
         pTrT2   /= rT2;
 
    crossProductVector3( deltaR[CD], forceTerm[dB], deltaR[CDxdB] );
    crossProductVector3( forceTerm[dB], deltaR[AD], deltaR[dBxAD] );
 
    if( fabs( pTrT2 ) > 1.0e-08 ){
       double delta2 = delta*2.0;
 
       crossProductVector3( deltaR[BD], deltaR[CD], deltaR[BDxCD] );
       crossProductVector3( deltaR[T],  deltaR[CD], deltaR[TxCD]  );
       crossProductVector3( deltaR[AD], deltaR[BD], deltaR[ADxBD] );
       crossProductVector3( deltaR[AD], deltaR[T],  deltaR[ADxT]  );
       for( int ii = 0; ii < 3; ii++ ){
 
          double term           = deltaR[BDxCD][ii] + delta2*deltaR[TxCD][ii];
          forceTerm[dA][ii]  += delta*deltaR[CDxdB][ii] + term*pTrT2;
 
               term           = deltaR[ADxBD][ii] + delta2*deltaR[ADxT][ii];
          forceTerm[dC][ii]  += delta*deltaR[dBxAD][ii] + term*pTrT2;
 
          forceTerm[dD][ii]  = -( forceTerm[dA][ii] + forceTerm[dB][ii] + forceTerm[dC][ii] );
       }
    } else {
       for( int ii = 0; ii < 3; ii++ ){
 
          forceTerm[dA][ii] += delta*deltaR[CDxdB][ii];
          forceTerm[dC][ii] += delta*deltaR[dBxAD][ii];
 
          forceTerm[dD][ii]  = -( forceTerm[dA][ii] + forceTerm[dB][ii] + forceTerm[dC][ii] );
       }
    }
 
    // accumulate forces and energy
 
    *energy                    += energyTerm;
 
    forces[particle1][0]       -= forceTerm[0][0];
    forces[particle1][1]       -= forceTerm[0][1];
    forces[particle1][2]       -= forceTerm[0][2];

    forces[particle2][0]       -= forceTerm[1][0];
    forces[particle2][1]       -= forceTerm[1][1];
    forces[particle2][2]       -= forceTerm[1][2];

    forces[particle3][0]       -= forceTerm[2][0];
    forces[particle3][1]       -= forceTerm[2][1];
    forces[particle3][2]       -= forceTerm[2][2];

    forces[particle4][0]       -= forceTerm[3][0];
    forces[particle4][1]       -= forceTerm[3][1];
    forces[particle4][2]       -= forceTerm[3][2];

}

static void computeAmoebaHarmonicInPlaneAngleForces( Context& context, AmoebaHarmonicInPlaneAngleForce& amoebaHarmonicInPlaneAngleForce,
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
    for( int ii = 0; ii < amoebaHarmonicInPlaneAngleForce.getNumAngles(); ii++ ){
        computeAmoebaHarmonicInPlaneAngleForce(ii, positions, amoebaHarmonicInPlaneAngleForce, expectedForces, expectedEnergy, log );
    }

    if( log ){
        (void) fprintf( log, "computeAmoebaHarmonicInPlaneAngleForces: expected energy=%14.7e\n", *expectedEnergy );
        for( unsigned int ii = 0; ii < positions.size(); ii++ ){
            (void) fprintf( log, "%6u [%14.7e %14.7e %14.7e]\n", ii, expectedForces[ii][0], expectedForces[ii][1], expectedForces[ii][2] );
        }
        (void) fflush( log );
    }
    return;

}

void compareWithExpectedForceAndEnergy( Context& context, AmoebaHarmonicInPlaneAngleForce& amoebaHarmonicInPlaneAngleForce,
                                        double tolerance, const std::string& idString, FILE* log) {

    std::vector<Vec3> expectedForces;
    double expectedEnergy;
    computeAmoebaHarmonicInPlaneAngleForces( context, amoebaHarmonicInPlaneAngleForce, expectedForces, &expectedEnergy, log );
   
    State state                      = context.getState(State::Forces | State::Energy);
    const std::vector<Vec3> forces   = state.getForces();

    if( log ){
        (void) fprintf( log, "computeAmoebaHarmonicInPlaneAngleForces: expected energy=%14.7e %14.7e\n", expectedEnergy, state.getPotentialEnergy() );
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

void testOneAngle( FILE* log ) {

    System system;
    int numberOfParticles = 4;
    for( int ii = 0; ii < numberOfParticles; ii++ ){
        system.addParticle(1.0);
    }

    LangevinIntegrator integrator(0.0, 0.1, 0.01);

    AmoebaHarmonicInPlaneAngleForce* amoebaHarmonicInPlaneAngleForce = new AmoebaHarmonicInPlaneAngleForce();

    double angle      = 65.0;
    double quadraticK = 1.0;
    double cubicK     = 0.0e-01;
    double quarticK   = 0.0e-02;
    double penticK    = 0.0e-03;
    double sexticK    = 0.0e-04;
    amoebaHarmonicInPlaneAngleForce->addAngle(0, 1, 2, 3, angle, quadraticK);

    amoebaHarmonicInPlaneAngleForce->setAmoebaGlobalHarmonicInPlaneAngleCubic(cubicK);
    amoebaHarmonicInPlaneAngleForce->setAmoebaGlobalHarmonicInPlaneAngleQuartic(quarticK);
    amoebaHarmonicInPlaneAngleForce->setAmoebaGlobalHarmonicInPlaneAnglePentic(penticK);
    amoebaHarmonicInPlaneAngleForce->setAmoebaGlobalHarmonicInPlaneAngleSextic(sexticK);

    system.addForce(amoebaHarmonicInPlaneAngleForce);
    Context context(system, integrator, Platform::getPlatformByName( "Cuda"));

    std::vector<Vec3> positions(numberOfParticles);

    positions[0] = Vec3(0, 1, 0);
    positions[1] = Vec3(0, 0, 0);
    positions[2] = Vec3(0, 0, 1);
    positions[3] = Vec3(1, 1, 1);

    context.setPositions(positions);
    compareWithExpectedForceAndEnergy( context, *amoebaHarmonicInPlaneAngleForce, TOL, "testOneInPlaneAngle", log );

}

int main( int numberOfArguments, char* argv[] ) {

    try {
        std::cout << "TestCudaAmoebaHarmonicInPlaneAngleForce running test..." << std::endl;
        Platform::loadPluginsFromDirectory( Platform::getDefaultPluginsDirectory() );
        FILE* log = NULL;
        //FILE* log = fopen( "AmoebaHarmonicInPlaneAngleForce.log", "w" );;

        testOneAngle( NULL );
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
