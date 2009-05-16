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
 * This tests the Brook RB torsion bond force/energy
 */

#include "../../../tests/AssertionUtilities.h"
#include "BrookPlatform.h"
#include "ReferencePlatform.h"
#include "openmm/OpenMMContext.h"
#include "openmm/RBTorsionForce.h"
#include "openmm/System.h"
#include "openmm/LangevinIntegrator.h"
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
   System system;
   system.addParticle(1.0);
   system.addParticle(1.0);
   system.addParticle(1.0);
   system.addParticle(1.0);
   LangevinIntegrator integrator( 0, 0.1, 0.01 );

   RBTorsionForce* forceField = new RBTorsionForce();
   forceField->addTorsion(0, 1, 2, 3, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6);
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
