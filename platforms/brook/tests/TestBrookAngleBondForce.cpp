/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2009 Stanford University and the Authors.           *
 * Authors: Peter Eastman, Mark Friedrichs                                    *
 * Contributors:                                                              *
 *                                                                            *
 * This program is free software: you can redistribute it and/or modify       *
 * it under the terms of the GNU Lesser General Public License as published   *
 * by the Free Software Foundation, either version 3 of the License, or       *
 * (at your option) any later version.                                        *
 *                                                                            *
 * This program is distributed in the hope that it will be useful,            *
 * but WITHOUT ANY WARRANTY; without even the implied warranty of             *
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the              *
 * GNU Lesser General Public License for more details.                        *
 *                                                                            *
 * You should have received a copy of the GNU Lesser General Public License   *
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.      *
 * -------------------------------------------------------------------------- */

#include <vector>

#include "../../../tests/AssertionUtilities.h"
#include "BrookPlatform.h"
#include "openmm/Context.h"
#include "openmm/HarmonicAngleForce.h"
#include "openmm/System.h"
#include "openmm/LangevinIntegrator.h"

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
fprintf( stderr, "%s xxxx\n", methodName.c_str() ); fflush( stderr );
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

   Context context(system, integrator, platform);
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
