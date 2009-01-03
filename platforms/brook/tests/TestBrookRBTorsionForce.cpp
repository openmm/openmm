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
 * This tests the Brook RB torsion bond force/energy
 */

#include "../../../tests/AssertionUtilities.h"
#include "BrookPlatform.h"
#include "ReferencePlatform.h"
#include "OpenMMContext.h"
#include "RBTorsionForce.h"
#include "System.h"
#include "LangevinIntegrator.h"
#include <vector>

#define PI_M               3.141592653589

using namespace OpenMM;
using namespace std;

const double TOL = 1e-5;

void testBrookRBTorsions( FILE* log ){

// ---------------------------------------------------------------------------------------

   static const std::string methodName      = "RBTorsions";
   int PrintOn                              = 0;
   
   int numberOfParticles                    = 4;
   double mass                              = 2.0;

// ---------------------------------------------------------------------------------------

   PrintOn = log ? PrintOn : 0;

   if( PrintOn ){
      (void) fprintf( log, "%s\n", methodName.c_str() ); (void) fflush( log );
   }

   BrookPlatform platform( 32, "cal", log );
   //ReferencePlatform platform;
   System system( numberOfParticles, 0 ); 
   LangevinIntegrator integrator( 0, 0.1, 0.01 );

   RBTorsionForce* forceField = new RBTorsionForce( 1 );
   forceField->setTorsionParameters(0, 0, 1, 2, 3, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6);
   system.addForce(forceField);

   OpenMMContext context(system, integrator, platform);
   vector<Vec3> positions(numberOfParticles);

   positions[0] = Vec3(0, 1, 0);
   positions[1] = Vec3(0, 0, 0);
   positions[2] = Vec3(1, 0, 0);
   positions[3] = Vec3(1, 1, 1);

   context.setPositions(positions);

   State state = context.getState( State::Forces | State::Energy );

   const vector<Vec3>& forces = state.getForces();
   if( PrintOn ){
      (void) fprintf( log, "RB torsion bond forces\n");
      for( int ii = 0; ii < numberOfParticles; ii++ ){
         (void) fprintf( log, "%d [%.5e %.5e %.5e]\n", ii, forces[ii][0], forces[ii][1], forces[ii][2] );
      }
      (void) fflush( log );
   }

   double psi    = 0.25*PI_M - PI_M;
   double torque = 0.0;
   for (int i = 1; i < 6; ++i) {
      double c  = 0.1*(i+1);
      torque   += -c*i*std::pow(std::cos(psi), i-1)*std::sin(psi);
   }

   if( PrintOn ){
      (void) fprintf( log, "RB torsion bond expected forces\n");
      (void) fprintf( log, "0 [0.0 0.0 %.5e]\n", torque );
      (void) fprintf( log, "3 [0.0 %.5e %.5e]\n", 0.5*torque, -0.5*torque );
      (void) fflush( log );
   }

   double tolerance = 0.001;
   ASSERT_EQUAL_VEC(Vec3(0, 0, torque), forces[0], tolerance );
   ASSERT_EQUAL_VEC(Vec3(0, 0.5*torque, -0.5*torque), forces[3], tolerance );
   ASSERT_EQUAL_VEC(Vec3(forces[0][0]+forces[1][0]+forces[2][0]+forces[3][0],
                         forces[0][1]+forces[1][1]+forces[2][1]+forces[3][1],
                         forces[0][2]+forces[1][2]+forces[2][2]+forces[3][2]),
                         Vec3(0, 0, 0), tolerance);

   double energy = 0.0;
   for (int i = 0; i < 6; ++i) {
       double c = 0.1*(i+1);
       energy += c*std::pow(std::cos(psi), i);
   }
   ASSERT_EQUAL_TOL(energy, state.getPotentialEnergy(), tolerance );

   if( PrintOn ){
      (void) fprintf( log, "RB torsion bond forces ok tolerance=%.2e\n", tolerance); fflush( log );
   }
}

int main( ){

// ---------------------------------------------------------------------------------------

   static const std::string methodName      = "testBrookRbTorsion";
   FILE* log                                = stdout;

// ---------------------------------------------------------------------------------------

   (void) fflush( stdout );
   (void) fflush( stderr );
   try {
      testBrookRBTorsions( log );
    } catch( const exception& e ){
      (void) fprintf( log, "Exception %s %.s\n", methodName.c_str(),  e.what() ); (void) fflush( log );
      return 1;
   }   
   (void) fprintf( log, "\n%s done\n", methodName.c_str() ); (void) fflush( log );

   return 0;
}
