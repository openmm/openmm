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
#include "HarmonicBondForce.h"
#include "NonbondedForce.h"
#include "System.h"
#include "VerletIntegrator.h"
#include "CMMotionRemover.h"

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

   System system( numberOfParticles, 0 );
   VerletIntegrator integrator(0.001);
   HarmonicBondForce* bonds = new HarmonicBondForce(1);
   bonds->setBondParameters(0, 2, 3, 2.0, 0.5);
   system.addForce(bonds);

   NonbondedForce* nonbonded = new NonbondedForce();
   for (int i = 0; i < numberOfParticles; ++i) {
       system.setParticleMass(i, (double) (i+1) );
       nonbonded->addParticle((i%2 == 0 ? 1.0 : -1.0), 1.0, 5.0);
   }
   system.addForce(nonbonded);

   CMMotionRemover* remover = new CMMotionRemover( 1 );
   system.addForce(remover);
   OpenMMContext context(system, integrator, platform);
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
