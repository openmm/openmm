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
 * This tests the Brook Langevin integrator
 */

#include <vector>

#include "../../../tests/AssertionUtilities.h"
#include "BrookPlatform.h"
#include "ReferencePlatform.h"
#include "OpenMMContext.h"
#include "HarmonicBondForce.h"
#include "NonBondedForce.h"
#include "CMMotionRemover.h"
#include "System.h"
#include "LangevinIntegrator.h"
#include "../src/sfmt/SFMT.h"
#include "../../reference/src/SimTKUtilities/SimTKOpenMMRealType.h"

using namespace OpenMM;
using namespace std;

const double TOL = 1e-5;

static OpenMMContext* testLangevinSingleBondSetup( int brookContext, LangevinIntegrator** outIntegrator, FILE* log ){

// ---------------------------------------------------------------------------------------

   static const std::string methodName      = "LangevinSingleBondSetup";
   int PrintOn                              = 1; 
   
   int numberOfParticles                    = 2;
   double mass                              = 2.0; 

// ---------------------------------------------------------------------------------------

   PrintOn = log ? PrintOn : 0;

   if( PrintOn ){
      (void) fprintf( log, "%s type=%s\n", methodName.c_str(), (brookContext?"Brook":"Reference") );
      (void) fflush( log );
   }   

   Platform* platform;
   if( brookContext ){
      platform = new BrookPlatform( 32, "cal", log );
      //platform = new BrookPlatform( 32, "cpu", log );
   } else {
      platform = new ReferencePlatform();
   }

   System* system = new System( numberOfParticles, 0);
   system->setParticleMass(0, mass );
   system->setParticleMass(1, mass );

   // double temperature, double frictionCoeff, double stepSize
   LangevinIntegrator* integrator = new LangevinIntegrator(0, 0.1, 0.001);
   integrator->setConstraintTolerance(1e-5);
   *outIntegrator                 = integrator;

   HarmonicBondForce* forceField  = new HarmonicBondForce(1);
   forceField->setBondParameters(0, 0, 1, 1.5, 1);
   system->addForce(forceField);

   OpenMMContext* context         = new OpenMMContext( *system, *integrator, *platform );

   vector<Vec3> positions(2);
   positions[0] = Vec3(-1, 0, 0);
   positions[1] = Vec3(1, 0, 0);
   context->setPositions(positions);
   
   return context;
}

void testLangevinSingleBond( FILE* log ){

// ---------------------------------------------------------------------------------------

   static const std::string methodName      = "LangevinSingleBond";
   int PrintOn                              = 1; 
   
// ---------------------------------------------------------------------------------------

   PrintOn = log ? PrintOn : 0;

   if( PrintOn ){
      (void) fprintf( log, "%s\n", methodName.c_str() ); (void) fflush( log );
   }   

   LangevinIntegrator* langevinIntegrator;
   OpenMMContext* context = testLangevinSingleBondSetup( 1, &langevinIntegrator, log ); 

   // This is simply a damped harmonic oscillator, so compare it to the analytical solution.
   
   double freq            = std::sqrt(1-0.05*0.05);
   int numberOfIterations = 1000;
   for (int i = 0; i < numberOfIterations; ++i) {
       State state = context->getState( State::Positions | State::Velocities );
       double time = state.getTime();
       double expectedDist = 1.5+0.5*std::exp(-0.05*time)*std::cos(freq*time);

       Vec3 pos1 = state.getPositions()[0];
       Vec3 pos2 = state.getPositions()[1];
       if( PrintOn > 1 ){
          (void) fprintf( log, "%s %d time=%.5e expD=%.5e pos=[%.5f %.5f %.5f] [%.5f %.5f %.5f] ", methodName.c_str(), i, time, -0.5*expectedDist, pos1[0], pos1[1], pos1[2], pos2[0], pos2[1], pos2[2] ); 
          (void) fflush( log );
       }

       ASSERT_EQUAL_VEC(Vec3(-0.5*expectedDist, 0, 0), state.getPositions()[0], 0.02);
       ASSERT_EQUAL_VEC(Vec3(0.5*expectedDist, 0, 0), state.getPositions()[1], 0.02);

       double expectedSpeed = -0.5*std::exp(-0.05*time)*(0.05*std::cos(freq*time)+freq*std::sin(freq*time));
       ASSERT_EQUAL_VEC(Vec3(-0.5*expectedSpeed, 0, 0), state.getVelocities()[0], 0.02);
       ASSERT_EQUAL_VEC(Vec3(0.5*expectedSpeed, 0, 0), state.getVelocities()[1], 0.02);

       Vec3 vel1 = state.getVelocities()[0];
       Vec3 vel2 = state.getVelocities()[1];
       if( PrintOn > 1 ){
          (void) fprintf( log, "expVel=%.5e vel=[%.5f %.5f %.5f] [%.5f %.5f %.5f]\n", -0.5*expectedSpeed, vel1[0], vel1[1], vel1[2], vel2[0], vel2[1], vel2[2] ); 
          (void) fflush( stdout );
       }

       langevinIntegrator->step(1);
   }
   
   if( PrintOn ){
      (void) fprintf( log, "%s 1 ok\n", methodName.c_str() ); fflush( log );
   }

   // Not set the friction to a tiny value and see if it conserves energy.
   
   langevinIntegrator->setFriction(5e-5);
   State state = context->getState(State::Energy);
   double potentialEnergy  = state.getPotentialEnergy();
   double kineticEnergy    = state.getKineticEnergy();
   double initialEnergy    = potentialEnergy + kineticEnergy;
   if( PrintOn ){
      (void) fprintf( log, "%s 2: initial energy: pot=%.5e ke=%.5e tot=%.5e\n", methodName.c_str(), potentialEnergy, kineticEnergy, initialEnergy );
      (void) fflush( log );
   }

   for (int i = 0; i < 1000; ++i) {

       state                   = context->getState(State::Energy);

       double potentialEnergy  = state.getPotentialEnergy();
       double kineticEnergy    = state.getKineticEnergy();
       double energy           = potentialEnergy + kineticEnergy;
       if( PrintOn > 1 ){
          (void) fprintf( log, "%s 2: energy: %d %.5e %.5e\n", methodName.c_str(), i, initialEnergy, energy, potentialEnergy, kineticEnergy );
          (void) fflush( log );
       }

       ASSERT_EQUAL_TOL( initialEnergy, energy, 0.01);

       langevinIntegrator->step(1);
   }

   if( PrintOn ){
      (void) fprintf( log, "%s 2 ok\n", methodName.c_str() ); fflush( log );
   }

   //delete langevinIntegrator;
   //delete context;
}

void testLangevinTemperature( FILE* log ){

// ---------------------------------------------------------------------------------------

   static const std::string methodName      = "LangevinTemperature";
   int PrintOn                              = 1; 
   const int numberOfParticles              = 8;
   double mass                              = 2.0; 
   const double temp                        = 100.0;

// ---------------------------------------------------------------------------------------

   PrintOn = log ? PrintOn : 0;

   if( PrintOn ){
      (void) fprintf( log, "%s\n", methodName.c_str() ); (void) fflush( log );
   }   

   BrookPlatform platform( 32, "cal", log );
   //ReferencePlatform platform;

   System system(numberOfParticles, 0);
   LangevinIntegrator integrator(temp, 0.2, 0.002);
   NonbondedForce* forceField = new NonbondedForce(numberOfParticles, 0);
   for (int i = 0; i < numberOfParticles; ++i){
      system.setParticleMass(i, mass );
      forceField->setParticleParameters(i, (i%2 == 0 ? 1.0 : -1.0), 1.0, 5.0);
   }
   system.addForce(forceField);

   CMMotionRemover* remover = new CMMotionRemover( 10 );
   system.addForce(remover);

   OpenMMContext context(system, integrator, platform);
   vector<Vec3> positions(numberOfParticles);
   for (int i = 0; i < numberOfParticles; ++i){
       positions[i] = Vec3((i%2 == 0 ? 2 : -2), (i%4 < 2 ? 2 : -2), (i < 4 ? 2 : -2));
   }
   context.setPositions(positions);
    
   // Let it equilibrate.
   
   integrator.step(10000);
   
   // Now run it for a while and see if the temperature is correct.
   
   double ke = 0.0;
   int steps = 1000;
   for( int i = 0; i < steps; ++i ){

       State state  = context.getState(State::Positions | State::Velocities | State::Energy);
       //State state  = context.getState(State::Energy);

       ke          += state.getKineticEnergy();
       if( PrintOn > 1 ){
          (void) fprintf( log, "%s %d KE=%12.5e ttl=%12.5e\n",
                          methodName.c_str(), i, state.getKineticEnergy(), ke );
          vector<Vec3> positions  = state.getPositions(); 
          vector<Vec3> velocities = state.getVelocities(); 
          double com[3] = { 0.0, 0.0, 0.0 };
          for( int ii = 0; ii < numberOfParticles; ii++ ){
             com[0] += velocities[ii][0];
             com[1] += velocities[ii][1];
             com[2] += velocities[ii][2];
             (void) fprintf( log, "   %d q[%12.5e %12.5e %12.5e] v[%12.5e %12.5e %12.5e]\n", ii,
                             positions[ii][0],  positions[ii][1],  positions[ii][2], 
                             velocities[ii][0], velocities[ii][1], velocities[ii][2] );
          }
          (void) fprintf( log, "VelCom[%12.5e %12.5e %12.5e]\n", com[0], com[1], com[2] );
          (void)fflush( log );
       }
       integrator.step(1);
   }
   ke               /= (double) steps;
   double expected   = 0.5*numberOfParticles*3.0*BOLTZ*temp;
   double tol        = 3*expected/std::sqrt(1000.0);

   double diff       = std::fabs( expected - ke );
   if( PrintOn ){
      (void) fprintf( log, "%s expected=%12.5e found=%12.5e diff=%12.5e tol=%12.5e\n", methodName.c_str(), expected, ke, diff, tol ); fflush( log );
   }

   ASSERT_EQUAL_TOL(expected, ke, 3*expected/std::sqrt(1000.0));
   if( PrintOn ){
      (void) fprintf( log, "%s ok\n", methodName.c_str(), expected, ke, diff, tol ); fflush( log );
   }

/*
/tests/AssertionUtilities.h
#define ASSERT_EQUAL_TOL(expected, found, tol){ 
   double _scale_ = std::fabs(expected) > 1.0 ? std::fabs(expected) : 1.0;
   if (std::fabs((expected)-(found))/_scale_ > (t    ol)) {std::stringstream details; details << "Expected "<<(expected)<<", found "<<(found); throwException(__FILE__, __LINE__, details.str());}};
    ASSERT_EQUAL_TOL(expected, ke, tol );
*/

}

void testLangevinConstraints( FILE* log ){

// ---------------------------------------------------------------------------------------

   static const std::string methodName      = "LangevinConstraints";
   int PrintOn                              = 1; 
   double mass                              = 1.0; 

// ---------------------------------------------------------------------------------------

   if( PrintOn ){
      (void) fprintf( log, "%s\n", methodName.c_str() );
      (void) fflush( log );
   }

   BrookPlatform platform( 32, "cal", log );

   const int numParticles   = 8;
   const int numConstraints = 4;
   const double temp        = 100.0;
   // ReferencePlatform platform;
   System system( numParticles, numConstraints );
   LangevinIntegrator integrator( temp, 2.0, 0.001 );
   integrator.setConstraintTolerance(1e-5);
   NonbondedForce* forceField = new NonbondedForce(numParticles, 0);
   for (int i = 0; i < numParticles; ++i) {
       system.setParticleMass(i, mass);
       forceField->setParticleParameters(i, (i%2 == 0 ? 0.2 : -0.2), 0.5, 5.0);
   }
   for (int i = 0; i < numConstraints; ++i){
       system.setConstraintParameters(i, 2*i, 2*i+1, 1.0);
   }
   system.addForce(forceField);

   CMMotionRemover* remover = new CMMotionRemover();
   system.addForce(remover);

   OpenMMContext context(system, integrator, platform);
   vector<Vec3> positions(numParticles);
   vector<Vec3> velocities(numParticles);
   init_gen_rand(0);

   for (int i = 0; i < numParticles; ++i) {
       positions[i]  = Vec3(i/2, (i+1)/2, 0);
       velocities[i] = Vec3(genrand_real2()-0.5, genrand_real2()-0.5, genrand_real2()-0.5);
   }
   context.setPositions(positions);
   context.setVelocities(velocities);
   
   // Simulate it and see whether the constraints remain satisfied.
   
   for (int i = 0; i < 1000; ++i) {
       State state = context.getState(State::Positions);
       for (int j = 0; j < numConstraints; ++j) {
          int particle1, particle2;
          double distance;
          system.getConstraintParameters(j, particle1, particle2, distance);
          Vec3 p1 = state.getPositions()[particle1];
          Vec3 p2 = state.getPositions()[particle2];
          double dist = std::sqrt((p1[0]-p2[0])*(p1[0]-p2[0])+(p1[1]-p2[1])*(p1[1]-p2[1])+(p1[2]-p2[2])*(p1[2]-p2[2]));
          if( PrintOn > 1 ){
             (void) fprintf( log, "%s %d %d dist=%12.5e %12.5e ok\n", methodName.c_str(), i, j, dist, fabs(dist-1.0) ); fflush( log );
          }
          ASSERT_EQUAL_TOL(1.0, dist, 2e-3);
       }
       integrator.step(1);
   }

   if( PrintOn ){
      (void) fprintf( log, "%s ok\n", methodName.c_str() ); fflush( log );
   }
}


int main( ){

// ---------------------------------------------------------------------------------------

   static const std::string methodName      = "testBrookLangevinIntegrator";
   FILE* log                                = stdout;

// ---------------------------------------------------------------------------------------

   (void) fflush( stdout );
   (void) fflush( stderr );
   try {
      testLangevinTemperature( log );
      testLangevinSingleBond( log );
      testLangevinConstraints( log );
    } catch( const exception& e ){
      (void) fprintf( log, "Exception %s %.s\n", methodName.c_str(),  e.what() ); (void) fflush( log );
      return 1;
   }   
   (void) fprintf( log, "\n%s done\n", methodName.c_str() ); (void) fflush( log );

   return 0;
}
