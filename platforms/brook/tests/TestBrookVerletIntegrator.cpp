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
#include "CMMotionRemover.h"
#include "System.h"
#include "VerletIntegrator.h"
#include "../src/sfmt/SFMT.h"

#define PI_M               3.141592653589

using namespace OpenMM;
using namespace std;

const double TOL = 1e-5;

void testVerletSingleBond( FILE* log ){

// ---------------------------------------------------------------------------------------

   static const std::string methodName      = "testVerletSingleBond";
   int PrintOn                              = 1; 
   
   int numberOfParticles                    = 2;
   double mass                              = 2.0; 

// ---------------------------------------------------------------------------------------

   PrintOn = log ? PrintOn : 0;

   if( PrintOn ){
      (void) fprintf( log, "%s\n", methodName.c_str() ); (void) fflush( log );
   }   

   BrookPlatform platform( 32, "cal", log );

   System system( numberOfParticles, 0 ); 
   system.setParticleMass(0, 2.0);
   system.setParticleMass(1, 2.0);

   VerletIntegrator integrator(0.001);

   HarmonicBondForce* forceField = new HarmonicBondForce(1);
   forceField->setBondParameters(0, 0, 1, 1.5, 1); 
   system.addForce(forceField);

//   CMMotionRemover* remover = new CMMotionRemover();
//   system.addForce(remover);

   OpenMMContext context(system, integrator, platform);
   vector<Vec3> positions(2);
   positions[0] = Vec3(-1, 0, 0); 
   positions[1] = Vec3(1, 0, 0); 
   context.setPositions(positions);
      
   // This is simply a harmonic oscillator, so compare it to the analytical solution.
      
   const double freq = 1.0;;
   State state = context.getState(State::Energy);
   const double initialEnergy = state.getKineticEnergy()+state.getPotentialEnergy();
   if( PrintOn ){
      (void) fprintf( log, "%s Energy initialEnergy=%12.5e KE=%12.5e PE=%12.5e\n",
                      methodName.c_str(), initialEnergy,
                      state.getKineticEnergy(), state.getPotentialEnergy() );
       (void) fflush( log );
   }


   for (int i = 0; i < 1000; ++i) {

       state               = context.getState(State::Positions | State::Velocities | State::Energy);
       double time         = state.getTime();
       double expectedDist = 1.5+0.5*std::cos(freq*time);

       Vec3 position0 = state.getPositions()[0];
       Vec3 position1 = state.getPositions()[1];
       if( PrintOn > 1 ){
          (void) fprintf( log, "%s %d Pos expected=[%12.5e 0 0] actual=[%12.5e %12.5e %12.5e] [%12.5e %12.5e %12.5e]\n",
                          methodName.c_str(), i, -0.5*expectedDist, 
                          position0[0], position0[1], position0[2],
                          position1[0], position1[1], position1[2] );
          (void) fflush( log );
       }

       ASSERT_EQUAL_VEC(Vec3(-0.5*expectedDist, 0, 0), state.getPositions()[0], 0.02);
       ASSERT_EQUAL_VEC(Vec3(0.5*expectedDist, 0, 0), state.getPositions()[1], 0.02);

       double expectedSpeed = -0.5*freq*std::sin(freq*time);

       Vec3 velocity0 = state.getVelocities()[0];
       Vec3 velocity1 = state.getVelocities()[1];
       if( PrintOn > 1 ){
          (void) fprintf( log, "%s %d Vel expected=[%12.5e 0 0] actual=[%12.5e %12.5e %12.5e] [%12.5e %12.5e %12.5e]\n",
                          methodName.c_str(), i, -0.5*expectedSpeed, 
                          velocity0[0], velocity0[1], velocity0[2],
                          velocity1[0], velocity1[1], velocity1[2] );
          (void) fflush( log );
       }

       ASSERT_EQUAL_VEC(Vec3(-0.5*expectedSpeed, 0, 0), state.getVelocities()[0], 0.02);
       ASSERT_EQUAL_VEC(Vec3(0.5*expectedSpeed, 0, 0), state.getVelocities()[1], 0.02);

       double energy = state.getKineticEnergy()+state.getPotentialEnergy();

       if( PrintOn > 1 ){
          (void) fprintf( log, "%s %d Energy initialEnergy=%12.5e actual=%12.5e KE=%12.5e PE=%12.5e\n",
                          methodName.c_str(), i, initialEnergy, energy,
                          state.getKineticEnergy(), state.getPotentialEnergy() );
          (void) fflush( log );
       }

       ASSERT_EQUAL_TOL(initialEnergy, energy, 0.01);
       integrator.step(1);
   }   

   if( PrintOn ){
      (void) fprintf( log, "%s ok\n", methodName.c_str() );
      (void) fflush( log );
   }

}

void testVerletConstraints( FILE* log ){

// ---------------------------------------------------------------------------------------

  static const std::string methodName      = "testVerletConstraints";
  int PrintOn                              = 1; 
  
  const int numParticles                   = 8;
  const int numConstraints                 = numParticles/2;
  double mass                              = 2.0; 
  const double temp                        = 100.0;

// ---------------------------------------------------------------------------------------

   PrintOn = log ? PrintOn : 0;

   if( PrintOn ){
      (void) fprintf( log, "%s\n", methodName.c_str() ); (void) fflush( log );
   }   

   //ReferencePlatform platform;
   BrookPlatform platform( 32, "cal", log );
   System system(numParticles, numConstraints);
   VerletIntegrator integrator(0.001);
   integrator.setConstraintTolerance(1e-5);
   NonbondedForce* forceField = new NonbondedForce(numParticles, 0);
   for (int i = 0; i < numParticles; ++i) {
       system.setParticleMass(i, 10.0);
       forceField->setParticleParameters(i, (i%2 == 0 ? 0.2 : -0.2), 0.5, 5.0);
   }
   for (int i = 0; i < numConstraints; ++i)
       system.setConstraintParameters(i, 2*i, 2*i+1, 1.0);
   system.addForce(forceField);

   CMMotionRemover* remover = new CMMotionRemover();
   system.addForce(remover);

   OpenMMContext context(system, integrator, platform);
   vector<Vec3> positions(numParticles);
   vector<Vec3> velocities(numParticles);
   init_gen_rand(0);
   for (int i = 0; i < numParticles; ++i) {
       positions[i] = Vec3(i/2, (i+1)/2, 0);
       velocities[i] = Vec3(genrand_real2()-0.5, genrand_real2()-0.5, genrand_real2()-0.5);
   }
   context.setPositions(positions);
   context.setVelocities(velocities);
   
   // Simulate it and see whether the constraints remain satisfied.
   
   double initialEnergy = 0.0;
   for (int i = 0; i < 1000; ++i) {
      State state = context.getState(State::Positions | State::Energy);
      for (int j = 0; j < numConstraints; ++j) {
         int particle1, particle2;
         double distance;
         system.getConstraintParameters(j, particle1, particle2, distance);
         Vec3 p1 = state.getPositions()[particle1];
         Vec3 p2 = state.getPositions()[particle2];
         double dist = std::sqrt((p1[0]-p2[0])*(p1[0]-p2[0])+(p1[1]-p2[1])*(p1[1]-p2[1])+(p1[2]-p2[2])*(p1[2]-p2[2]));
          if( PrintOn ){
             (void) fprintf( log, "%s step=%d constraint=%d p[%d %d] d=%.5e exptd=%.5e [%.5e %.5e %.5e] [%.5e %.5e %.5e]e\n", 
                             methodName.c_str(), i, j, particle1, particle2, dist, distance,
                             p1[0], p1[1], p1[2], p2[0], p2[1], p2[2]); (void) fflush( log );
         }
         ASSERT_EQUAL_TOL(distance, dist, 2e-2);
       }

       double energy = state.getKineticEnergy()+state.getPotentialEnergy();
       if( PrintOn ){
          (void) fprintf( log, "%s %d e[%.5e %.5e] ke=%.5e pe=%.5e\n", 
                          methodName.c_str(), i, initialEnergy, energy, state.getKineticEnergy(), state.getPotentialEnergy() ); (void) fflush( log );
       }
       if (i == 1)
           initialEnergy = energy;
       else if (i > 1)
           ASSERT_EQUAL_TOL(initialEnergy, energy, 0.5);
       integrator.step(1);
   }
   if( PrintOn ){
      (void) fprintf( log, "%s ok\n", methodName.c_str() );
      (void) fflush( log );
   }
}

int main( ){

// ---------------------------------------------------------------------------------------

   static const std::string methodName      = "testBrookVerletIntegrator";
   FILE* log                                = stdout;

// ---------------------------------------------------------------------------------------

   (void) fflush( stdout );
   (void) fflush( stderr );
   try {
      testVerletSingleBond( log );
      testVerletConstraints( log );
    } catch( const exception& e ){
      (void) fprintf( log, "Exception %s %.s\n", methodName.c_str(),  e.what() ); (void) fflush( log );
      return 1;
   }   
   (void) fprintf( log, "\n%s done\n", methodName.c_str() ); (void) fflush( log );

   return 0;
}
