/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2008 Stanford University and the Authors.           *
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
 * This tests the Brook harmonic angle bond force/energy
 */

#include <vector>

#include "../../../tests/AssertionUtilities.h"
#include "BrookPlatform.h"
#include "OpenMMContext.h"
#include "HarmonicAngleForce.h"
#include "System.h"
#include "LangevinIntegrator.h"

#define PI_M               3.141592653589

using namespace OpenMM;
using namespace std;

const double TOL = 1e-5;

void testBrookAngles( FILE* log ){

// ---------------------------------------------------------------------------------------

   static const std::string methodName      = "testBrookAngles";
   static const int debug                   = 1;
   
   int PrintOn                              = 0; 
   int numberOfParticles                    = 4;
   double mass                              = 2.0;

// ---------------------------------------------------------------------------------------

   PrintOn = log ? PrintOn : 0;

   if( PrintOn ){
      (void) fprintf( log, "%s\n", methodName.c_str() ); (void) fflush( log );
   }

   BrookPlatform platform( 32, "cal", log );
   System system;
   for (int i = 0; i < numberOfParticles; i++)
       system.addParticle(1.0);
   LangevinIntegrator integrator( 0, 0.1, 0.01 );

   HarmonicAngleForce* forceField = new HarmonicAngleForce();

   // int atom1, int atom2, int atom3, double angle, double k
   forceField->addAngle(0, 1, 2, PI_M/3, 1.1);
   forceField->addAngle(1, 2, 3, PI_M/2, 1.2);
   system.addForce(forceField);

   OpenMMContext context(system, integrator, platform);
   vector<Vec3> positions(numberOfParticles);

   positions[0] = Vec3(0, 1, 0);
   positions[1] = Vec3(0, 0, 0);
   positions[2] = Vec3(1, 0, 0);
   positions[3] = Vec3(2, 1, 0);

   context.setPositions(positions);

   State state = context.getState( State::Forces | State::Energy );

   const vector<Vec3>& forces = state.getForces();
   if( PrintOn ){
      (void) fprintf( log, "Angle bond forces\n");
      for( int ii = 0; ii < numberOfParticles; ii++ ){
         (void) fprintf( log, "%d [%.5e %.5e %.5e]\n", ii, forces[ii][0], forces[ii][1], forces[ii][2] );
      }
      (void) fflush( log );
   }

   double tolerance  = 1.0e-03;
   double torque1    = 1.1*PI_M/6;
   double torque2    = 1.2*PI_M/4;
   ASSERT_EQUAL_VEC(Vec3(torque1, 0, 0), forces[0], tolerance);
   ASSERT_EQUAL_VEC(Vec3(-0.5*torque2, 0.5*torque2, 0), forces[3], tolerance); // reduced by sqrt(2) due to the bond length, another sqrt(2) due to the angle
   ASSERT_EQUAL_VEC(Vec3(forces[0][0]+forces[1][0]+forces[2][0]+forces[3][0],
                         forces[0][1]+forces[1][1]+forces[2][1]+forces[3][1],
                         forces[0][2]+forces[1][2]+forces[2][2]+forces[3][2]),
                         Vec3(0, 0, 0), tolerance);

   ASSERT_EQUAL_TOL(0.5*1.1*(PI_M/6)*(PI_M/6) + 0.5*1.2*(PI_M/4)*(PI_M/4), state.getPotentialEnergy(), tolerance);

   if( PrintOn ){
      (void) fprintf( log, "Angle bond forces ok tolerance=%.2e\n", tolerance ); (void) fflush( log );
   }
}

int main( ){

// ---------------------------------------------------------------------------------------

   static const std::string methodName      = "testBrookAngles";
   FILE* log                                = stdout;

// ---------------------------------------------------------------------------------------

   (void) fflush( stdout );
   (void) fflush( stderr );
   try {
      testBrookAngles( log );
    } catch( const exception& e ){
      (void) fprintf( log, "Exception %s %.s\n", methodName.c_str(),  e.what() ); (void) fflush( log );
      return 1;
   }   
   (void) fprintf( log, "\n%s done\n", methodName.c_str() ); (void) fflush( log );

   return 0;
}
