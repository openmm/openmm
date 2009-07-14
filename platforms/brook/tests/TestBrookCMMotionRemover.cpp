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
#include "openmm/HarmonicBondForce.h"
#include "openmm/NonbondedForce.h"
#include "openmm/System.h"
#include "openmm/VerletIntegrator.h"
#include "openmm/CMMotionRemover.h"

#include "../src/sfmt/SFMT.h"
#include "../../reference/src/SimTKUtilities/SimTKOpenMMRealType.h"

using namespace OpenMM;
using namespace std;

const double TOL = 1e-5;

Vec3 calcCM(const vector<Vec3>& values, System& system) {
    Vec3 cm;
    for (int j = 0; j < system.getNumParticles(); ++j) {
        cm[0] += values[j][0]*system.getParticleMass(j);
        cm[1] += values[j][1]*system.getParticleMass(j);
        cm[2] += values[j][2]*system.getParticleMass(j);
    }
    return cm;
}

void testMotionRemoval( FILE* log ) {

// ---------------------------------------------------------------------------------------

   static const std::string methodName      = "testMotionRemoval";
   int PrintOn                              = 1;
   int numberOfParticles                    = 8;
   double mass                              = 2.0; 

// ---------------------------------------------------------------------------------------

   PrintOn = log ? PrintOn : 0;

   if( PrintOn ){
      (void) fprintf( log, "%s\n", methodName.c_str() ); (void) fflush( log );
   }

   //ReferencePlatform platform;
   BrookPlatform platform( 32, "cal", log );

   System system;
   VerletIntegrator integrator(0.001);
   HarmonicBondForce* bonds = new HarmonicBondForce();
   bonds->addBond(2, 3, 2.0, 0.5);
   system.addForce(bonds);

   NonbondedForce* nonbonded = new NonbondedForce();
   for (int i = 0; i < numberOfParticles; ++i) {
       system.addParticle((double) (i+1) );
       nonbonded->addParticle((i%2 == 0 ? 1.0 : -1.0), 1.0, 5.0);
   }
   system.addForce(nonbonded);

   CMMotionRemover* remover = new CMMotionRemover( 1 );
   system.addForce(remover);
   Context context(system, integrator, platform);
   vector<Vec3> positions(numberOfParticles);
   vector<Vec3> velocities(numberOfParticles);
   init_gen_rand(0);
   for (int i = 0; i < numberOfParticles; ++i) {
       positions[i] = Vec3((i%2 == 0 ? 2 : -2), (i%4 < 2 ? 2 : -2), (i < 4 ? 2 : -2));
       velocities[i] = Vec3(genrand_real2()-0.5, genrand_real2()-0.5, genrand_real2()-0.5);
   }
   context.setPositions(positions);
   context.setVelocities(velocities);

   // Now run it for a while and see if the center of mass remains fixed.

   Vec3 cmPos = calcCM(context.getState(State::Positions).getPositions(), system);
   for (int i = 0; i < 1000; ++i) {
       integrator.step(1);
       State state = context.getState(State::Positions | State::Velocities);
       Vec3 pos = calcCM(state.getPositions(), system);
       ASSERT_EQUAL_VEC(cmPos, pos, 1e-2);
       Vec3 vel = calcCM(state.getVelocities(), system);
       ASSERT_EQUAL_VEC(Vec3(0, 0, 0), vel, 1e-2);
   }
   if( PrintOn ){
      (void) fprintf( log, "%s ok\n", methodName.c_str() ); (void) fflush( log );
   }
}

int main( ){

// ---------------------------------------------------------------------------------------

   static const std::string methodName      = "RemoveCMMotion";
   FILE* log                                = stdout;

// ---------------------------------------------------------------------------------------

   (void) fflush( stdout );
   (void) fflush( stderr );
   try {
      testMotionRemoval( log );
    } catch( const exception& e ){
      (void) fprintf( log, "Exception %s %.s\n", methodName.c_str(),  e.what() ); (void) fflush( log );
      return 1;
   }   
   (void) fprintf( log, "\n%s done\n", methodName.c_str() ); (void) fflush( log );

   return 0;
}
