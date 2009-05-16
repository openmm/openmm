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

#include "../../../tests/AssertionUtilities.h"
#include "BrookPlatform.h"
#include "openmm/OpenMMContext.h"
#include "openmm/HarmonicBondForce.h"
#include "openmm/System.h"
#include "openmm/LangevinIntegrator.h"
#include <vector>

using namespace OpenMM;
using namespace std;

const double TOL = 1e-5;

void testBrookBonds( FILE* log ){

// ---------------------------------------------------------------------------------------

   static const std::string methodName      = "testBrookBonds";
   int PrintOn                              = 0;
   int numberOfParticles                    = 3;
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
   LangevinIntegrator integrator(0, 0.1, 0.01);

   // int numParticles, int numBonds, int numAngles, int numPeriodicTorsions, int numRBTorsions

   HarmonicBondForce* forceField = new HarmonicBondForce(); 

   // ( index, atom1, atom2, length, k )
   forceField->addBond(0, 1, 1.5, 0.8);
   forceField->addBond(1, 2, 1.2, 0.7);
   system.addForce(forceField);

   OpenMMContext context(system, integrator, platform);
   vector<Vec3> positions(3);

   positions[0] = Vec3(0, 2, 0); 
   positions[1] = Vec3(0, 0, 0); 
   positions[2] = Vec3(1, 0, 0); 

   context.setPositions(positions);

   //(void) fprintf( log, "testBrookBonds:Calling getState\n");
   //(void) fflush( log );

   State state = context.getState( State::Forces | State::Energy );

   const vector<Vec3>& forces = state.getForces();
   if( PrintOn ){
      (void) fprintf( log, "Harmonic bond forces\n");
      for( int ii = 0; ii < numberOfParticles; ii++ ){
         (void) fprintf( log, "%d [%.5e %.5e %.5e]\n", ii, forces[ii][0], forces[ii][1], forces[ii][2] );
      }
      (void) fflush( log );
   }

   ASSERT_EQUAL_VEC(Vec3(0, -0.8*0.5, 0), forces[0], TOL);
   ASSERT_EQUAL_VEC(Vec3(0.7*0.2, 0, 0), forces[2], TOL);
   ASSERT_EQUAL_VEC(Vec3(-forces[0][0]-forces[2][0], -forces[0][1]-forces[2][1], -forces[0][2]-forces[2][2]), forces[1], TOL);

   if( PrintOn ){
      (void) fprintf( log, "Harmonic bond forces ok -- checking energy\n"); (void) fflush( log );
   }

   ASSERT_EQUAL_TOL( 0.5*0.8*0.5*0.5 + 0.5*0.7*0.2*0.2, state.getPotentialEnergy(), TOL);

   if( PrintOn ){
      (void) fprintf( log, "Harmonic bond energy ok\n"); (void) fflush( log );
   }

}

int main( ){

// ---------------------------------------------------------------------------------------

   static const std::string methodName      = "testBrookBonds";
   FILE* log                                = stdout;

// ---------------------------------------------------------------------------------------

   (void) fflush( stdout );
   (void) fflush( stderr );
   try {
      testBrookBonds( log );
    }  catch( const exception& e ){
      (void) fprintf( log, "Exception %s %.s\n", methodName.c_str(), e.what() ); (void) fflush( log );
      return 1;
   }   
   (void) fprintf( log, "\n%s done\n", methodName.c_str() ); (void) fflush( log );

   return 0;
}
