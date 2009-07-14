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

/**
 * This tests the Brook periodic torsion bond force/energy
 */

#include "../../../tests/AssertionUtilities.h"
#include "BrookPlatform.h"
#include "openmm/Context.h"
#include "openmm/PeriodicTorsionForce.h"
#include "openmm/System.h"
#include "openmm/LangevinIntegrator.h"
#include <vector>

#define PI_M               3.141592653589

using namespace OpenMM;
using namespace std;

const double TOL = 1e-5;



void testBrookPeriodicTorsions( FILE* log ){

// ---------------------------------------------------------------------------------------

   static const std::string methodName      = "PeriodicTorsions";
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
   system.addParticle(1.0);
   system.addParticle(1.0);
   system.addParticle(1.0);
   system.addParticle(1.0);
   LangevinIntegrator integrator( 0, 0.1, 0.01 );

   // int numParticles, int numBonds, int numAngles, int numPeriodicTorsions, int numRBTorsions

   PeriodicTorsionForce* forceField = new PeriodicTorsionForce(); 

   // int index, int atom1, int atom2, int atom3, double angle, double k
   forceField->addTorsion(0, 1, 2, 3, 2, PI_M/3, 1.1);
   system.addForce(forceField);

   Context context(system, integrator, platform);
   vector<Vec3> positions(numberOfParticles);

   positions[0] = Vec3(0, 1, 0);
   positions[1] = Vec3(0, 0, 0);
   positions[2] = Vec3(1, 0, 0);
   positions[3] = Vec3(1, 0, 2);

   context.setPositions(positions);

   State state = context.getState( State::Forces | State::Energy );

   const vector<Vec3>& forces = state.getForces();
   if( PrintOn ){
      (void) fprintf( log, "Periodic torsion bond forces\n");
      for( int ii = 0; ii < numberOfParticles; ii++ ){
         (void) fprintf( log, "%d [%.5e %.5e %.5e]\n", ii, forces[ii][0], forces[ii][1], forces[ii][2] );
      }
      (void) fflush( log );
   }

   double torque    = -2*1.1*std::sin(2*PI_M/3);
   double tolerance = 1.0e-03;
   ASSERT_EQUAL_VEC(Vec3(0, 0, torque), forces[0], tolerance );
   ASSERT_EQUAL_VEC(Vec3(0, 0.5*torque, 0), forces[3], tolerance );
   ASSERT_EQUAL_VEC(Vec3(forces[0][0]+forces[1][0]+forces[2][0]+forces[3][0],
                         forces[0][1]+forces[1][1]+forces[2][1]+forces[3][1],
                         forces[0][2]+forces[1][2]+forces[2][2]+forces[3][2]),
                         Vec3(0, 0, 0), tolerance );

   ASSERT_EQUAL_TOL(1.1*(1+std::cos(2*PI_M/3)), state.getPotentialEnergy(), tolerance );
   if( PrintOn ){
      (void) fprintf( log, "Periodic torsion bond forces ok -- tolerance=%.2e\n", tolerance); (void) fflush( log );
   }
}

int main( ){

// ---------------------------------------------------------------------------------------

   static const std::string methodName      = "testBrookPeriodicTorsions";
   FILE* log                                = stdout;

// ---------------------------------------------------------------------------------------

   (void) fflush( stdout );
   (void) fflush( stderr );
   try {
      testBrookPeriodicTorsions( log );
    } catch( const exception& e ){
      (void) fprintf( log, "Exception %s %.s\n", methodName.c_str(),  e.what() ); (void) fflush( log );
      return 1;
   }   
   (void) fprintf( log, "\n%s done\n", methodName.c_str() ); (void) fflush( log );

   return 0;
}
