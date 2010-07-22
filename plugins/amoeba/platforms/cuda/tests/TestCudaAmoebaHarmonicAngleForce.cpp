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
 * This tests the Cuda implementation of CudaAmoebaHarmonicAngleForce.
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
#define RADIAN_TO_DEGREE  57.29577951308
#define DEGREE_TO_RADIAN   0.01745329252
#define RADIAN_INVERSE     0.01745329252

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

static void getPrefactorsGivenAngleCosine( double cosine, double idealAngle, double quadraticK, double cubicK,
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
        (void) fprintf( log, "getPrefactorsGivenAngleCosine: cosine=%10.3e angle=%10.3e ideal=%10.3e\n", cosine, angle, idealAngle ); 
        (void) fflush( log );
    }

    double deltaIdeal         = angle - idealAngle;
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

static void computeAmoebaHarmonicAngleForce(int bondIndex,  std::vector<Vec3>& positions, AmoebaHarmonicAngleForce& amoebaHarmonicAngleForce,
                                             std::vector<Vec3>& forces, double* energy, FILE* log ) {

    int particle1, particle2, particle3;
    double idealAngle;
    double quadraticK;
    amoebaHarmonicAngleForce.getAngleParameters(bondIndex, particle1, particle2, particle3, idealAngle, quadraticK );

    double cubicK         = amoebaHarmonicAngleForce.getAmoebaGlobalHarmonicAngleCubic();
    double quarticK       = amoebaHarmonicAngleForce.getAmoebaGlobalHarmonicAngleQuartic();
    double penticK        = amoebaHarmonicAngleForce.getAmoebaGlobalHarmonicAnglePentic();
    double sexticK        = amoebaHarmonicAngleForce.getAmoebaGlobalHarmonicAngleSextic();

    if( log ){
        (void) fprintf( log, "computeAmoebaHarmonicAngleForce: bond %d [%d %d %d] ang=%10.3f k=%10.3f [%10.3e %10.3e %10.3e %10.3e]\n", 
                             bondIndex, particle1, particle2, particle3, idealAngle, quadraticK, cubicK, quarticK, penticK, sexticK );
        (void) fflush( log );
    }

    double deltaR[2][3];
    double r2_0 = 0.0;
    double r2_1 = 0.0;
    for( int ii = 0; ii < 3; ii++ ){

           deltaR[0][ii]    = positions[particle1][ii] - positions[particle2][ii];
           r2_0            += deltaR[0][ii]*deltaR[0][ii];

           deltaR[1][ii]    = positions[particle3][ii] - positions[particle2][ii];
           r2_1            += deltaR[1][ii]*deltaR[1][ii];

    }

    double pVector[3];
    crossProductVector3( deltaR[0], deltaR[1], pVector );
    double rp      = sqrt( pVector[0]*pVector[0] + pVector[1]*pVector[1] + pVector[2]*pVector[2] );
    if( rp < 1.0e-06 ){
       rp = 1.0e-06;
    }   
    double dot    = deltaR[0][0]*deltaR[1][0] + deltaR[0][1]*deltaR[1][1] + deltaR[0][2]*deltaR[1][2];
    double cosine = dot/sqrt(r2_0*r2_1);

    if( log ){
        (void) fprintf( log, "dot=%10.3e r2_0=%10.3e r2_1=%10.3e\n", dot, r2_0, r2_1 ); 
        (void) fflush( log );
    }

    double dEdR;
    double energyTerm;
    getPrefactorsGivenAngleCosine( cosine, idealAngle, quadraticK, cubicK,
                                   quarticK, penticK, sexticK, &dEdR, &energyTerm, log );

    double termA  = -dEdR/(r2_0*rp);
    double termC  =  dEdR/(r2_1*rp);

    double deltaCrossP[3][3];
    crossProductVector3( deltaR[0], pVector, deltaCrossP[0] );
    crossProductVector3( deltaR[1], pVector, deltaCrossP[2] );
    for( int ii = 0; ii < 3; ii++ ){
        deltaCrossP[0][ii] *= termA;
        deltaCrossP[2][ii] *= termC;
        deltaCrossP[1][ii]  = -1.0*(deltaCrossP[0][ii] + deltaCrossP[2][ii]);
    }

    forces[particle1][0]       += deltaCrossP[0][0];
    forces[particle1][1]       += deltaCrossP[0][1];
    forces[particle1][2]       += deltaCrossP[0][2];

    forces[particle2][0]       += deltaCrossP[1][0];
    forces[particle2][1]       += deltaCrossP[1][1];
    forces[particle2][2]       += deltaCrossP[1][2];

    forces[particle3][0]       += deltaCrossP[2][0];
    forces[particle3][1]       += deltaCrossP[2][1];
    forces[particle3][2]       += deltaCrossP[2][2];

    *energy                    += energyTerm;
}

static void computeAmoebaHarmonicAngleForces( Context& context, AmoebaHarmonicAngleForce& amoebaHarmonicAngleForce,
                                             std::vector<Vec3>& expectedForces, double* expectedEnergy, FILE* log ) {

    // get positions and zero forces

    State state = context.getState(State::Positions);
    std::vector<Vec3> positions = state.getPositions();
    expectedForces.resize( positions.size() );
    
    for( unsigned int ii = 0; ii < expectedForces.size(); ii++ ){
        expectedForces[ii][0] = expectedForces[ii][1] = expectedForces[ii][2] = 0.0;
    }

    // calculates forces/energy

    *expectedEnergy = 0.0;
    for( int ii = 0; ii < amoebaHarmonicAngleForce.getNumAngles(); ii++ ){
        computeAmoebaHarmonicAngleForce(ii, positions, amoebaHarmonicAngleForce, expectedForces, expectedEnergy, log );
    }

    if( log ){
        (void) fprintf( log, "computeAmoebaHarmonicAngleForces: expected energy=%14.7e\n", *expectedEnergy );
        for( unsigned int ii = 0; ii < positions.size(); ii++ ){
            (void) fprintf( log, "%6u [%14.7e %14.7e %14.7e]\n", ii, expectedForces[ii][0], expectedForces[ii][1], expectedForces[ii][2] );
        }
        (void) fflush( log );
    }
    return;

}

void compareWithExpectedForceAndEnergy( Context& context, AmoebaHarmonicAngleForce& amoebaHarmonicAngleForce,
                                        double tolerance, const std::string& idString, FILE* log) {

    std::vector<Vec3> expectedForces;
    double expectedEnergy;
    computeAmoebaHarmonicAngleForces( context, amoebaHarmonicAngleForce, expectedForces, &expectedEnergy, log );
   
    State state                      = context.getState(State::Forces | State::Energy);
    const std::vector<Vec3> forces   = state.getForces();

    if( log ){
        (void) fprintf( log, "computeAmoebaHarmonicAngleForces: expected energy=%14.7e %14.7e\n", expectedEnergy, state.getPotentialEnergy() );
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
    int numberOfParticles = 3;
    for( int ii = 0; ii < numberOfParticles; ii++ ){
        system.addParticle(1.0);
    }

    LangevinIntegrator integrator(0.0, 0.1, 0.01);

    AmoebaHarmonicAngleForce* amoebaHarmonicAngleForce = new AmoebaHarmonicAngleForce();

    double angle      = 100.0;
    double quadraticK = 1.0;
    double cubicK     = 1.0e-01;
    double quarticK   = 1.0e-02;
    double penticK    = 1.0e-03;
    double sexticK    = 1.0e-04;
    amoebaHarmonicAngleForce->addAngle(0, 1, 2, angle, quadraticK);

    amoebaHarmonicAngleForce->setAmoebaGlobalHarmonicAngleCubic(cubicK);
    amoebaHarmonicAngleForce->setAmoebaGlobalHarmonicAngleQuartic(quarticK);
    amoebaHarmonicAngleForce->setAmoebaGlobalHarmonicAnglePentic(penticK);
    amoebaHarmonicAngleForce->setAmoebaGlobalHarmonicAngleSextic(sexticK);

    system.addForce(amoebaHarmonicAngleForce);
    Context context(system, integrator, Platform::getPlatformByName( "Cuda"));

    std::vector<Vec3> positions(numberOfParticles);

    positions[0] = Vec3(0, 1, 0);
    positions[1] = Vec3(0, 0, 0);
    positions[2] = Vec3(0, 0, 1);

    context.setPositions(positions);
    compareWithExpectedForceAndEnergy( context, *amoebaHarmonicAngleForce, TOL, "testOneAngle", log );

}

int main( int numberOfArguments, char* argv[] ) {

    try {
        std::cout << "Running test..." << std::endl;
        std::string openmmPluginDirectory = "/home/friedrim/src/openmm/trunk/OpenMM/bin";
        Platform::loadPluginsFromDirectory( openmmPluginDirectory );
        FILE* log = fopen( "AmoebaHarmonicAngleForce.log", "w" );;

        testOneAngle( log );
        fclose( log );

    }
    catch(const std::exception& e) {
        std::cout << "exception: " << e.what() << std::endl;
        std::cout << "FAIL - ERROR.  Test failed." << std::endl;
        return 1;
    }
    std::cout << "PASS - Test succeeded." << std::endl;
    return 0;
}
