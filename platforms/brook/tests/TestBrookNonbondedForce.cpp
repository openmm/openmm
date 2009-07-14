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
 * This tests the Brook harmonic angle bond force/energy
 */

#include <vector>

#include "../../../tests/AssertionUtilities.h"
#include "BrookPlatform.h"
#include "openmm/Context.h"
#include "openmm/NonbondedForce.h"
#include "openmm/HarmonicBondForce.h"
#include "openmm/System.h"
#include "openmm/LangevinIntegrator.h"

#define PI_M               3.141592653589

using namespace OpenMM;
using namespace std;

const double TOL = 1e-5;

void testBrookCoulomb( FILE* log ){

// ---------------------------------------------------------------------------------------

   static const std::string methodName      = "Coulomb";
   int PrintOn                              = 0; 
   int numberOfParticles                    = 2;
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
   LangevinIntegrator integrator( 0, 0.1, 0.01 );

   // int index, double charge, double radius, double depth

   NonbondedForce* forceField = new NonbondedForce(); 
   forceField->addParticle(0.5, 1, 0);
   forceField->addParticle(-1.5, 1, 0);
   system.addForce(forceField);

   //(void) fprintf( log, "%s: Calling context\n",  methodName.c_str() );
   //(void) fflush( log );

   Context context(system, integrator, platform);
   vector<Vec3> positions(numberOfParticles);

   positions[0] = Vec3(0, 0, 0);
   positions[1] = Vec3(2, 0, 0);

   context.setPositions(positions);

   //(void) fprintf( log, "%s :Calling getState\n", methodName.c_str() );
   //(void) fflush( log );

   State state = context.getState( State::Forces | State::Energy );

   const vector<Vec3>& forces = state.getForces();
   if( PrintOn ){
      (void) fprintf( log, "\nCoulomb forces\n");
      for( int ii = 0; ii < numberOfParticles; ii++ ){
         (void) fprintf( log, "%d [%.5e %.5e %.5e]\n", ii, forces[ii][0], forces[ii][1], forces[ii][2] );
      }
      (void) fflush( log );
   }

   double force = 138.935485*(-0.75)/4.0;
   ASSERT_EQUAL_VEC(Vec3(-force, 0, 0), forces[0], TOL);
   ASSERT_EQUAL_VEC(Vec3(force, 0, 0), forces[1], TOL);
   ASSERT_EQUAL_TOL(138.935485*(-0.75)/2.0, state.getPotentialEnergy(), TOL);

   if( PrintOn ){
      (void) fprintf( log, "Coulomb forces ok\n"); fflush( log );
   }

   // delete forceField;

}

void testBrookLJ( FILE* log ){

// ---------------------------------------------------------------------------------------

   static const std::string methodName      = "LJ";
   
   int PrintOn                              = 0; 
   int numberOfParticles                    = 2;
   double mass                              = 2.0;

// ---------------------------------------------------------------------------------------

   PrintOn = log ? PrintOn : 0;

   if( PrintOn ){
      (void) fprintf( log, "%s\n", methodName.c_str() ); (void) fflush( log );
   }   

   BrookPlatform platform( 32, "cal", log );
   // ReferencePlatform platform;
   System system;
   system.addParticle(1.0);
   system.addParticle(1.0);
   LangevinIntegrator integrator( 0, 0.1, 0.01 );

   // int index, double charge, double radius, double depth

   NonbondedForce* forceField = new NonbondedForce(); 
   forceField->addParticle(0, 1.2, 1);
   forceField->addParticle(0, 1.4, 2);
   system.addForce(forceField);

   //(void) fprintf( log, "%s: Calling context\n",  methodName.c_str() );
   //(void) fflush( log );

   Context context(system, integrator, platform);
   vector<Vec3> positions(numberOfParticles);

   positions[0] = Vec3(0, 0, 0);
   positions[1] = Vec3(2, 0, 0);

   context.setPositions(positions);

   //(void) fprintf( log, "%s :Calling getState\n", methodName.c_str() );
   //(void) fflush( log );

   State state = context.getState( State::Forces | State::Energy );

   const vector<Vec3>& forces = state.getForces();
   if( PrintOn ){
      (void) fprintf( log, "LJ forces\n");
      for( int ii = 0; ii < numberOfParticles; ii++ ){
         (void) fprintf( log, "%d [%.5e %.5e %.5e]\n", ii, forces[ii][0], forces[ii][1], forces[ii][2] );
      }
      (void) fflush( log );
   }

   double x   = 1.3/2.0;
   double eps = sqrt( 2.0 );
   double force = 4.0*eps*(12*std::pow(x, 12.0)-6*std::pow(x, 6.0))/2.0;
   ASSERT_EQUAL_VEC(Vec3(-force, 0, 0), forces[0], TOL);
   ASSERT_EQUAL_VEC(Vec3(force, 0, 0), forces[1], TOL);
   ASSERT_EQUAL_TOL(4.0*eps*(std::pow(x, 12.0)-std::pow(x, 6.0)), state.getPotentialEnergy(), TOL);

   if( PrintOn ){
      (void) fprintf( log, "LJ forces ok\n"); fflush( log );
   }

}

void testBrookExclusionsAnd14( FILE* log ){

// ---------------------------------------------------------------------------------------

   static const std::string methodName      = "ExclusionsAnd14";
   int numberOfParticles                    = 5;
   int PrintOn                              = 0; 
   double mass                              = 2.0;

// ---------------------------------------------------------------------------------------

   PrintOn = log ? PrintOn : 0;

   if( PrintOn ){
      (void) fprintf( log, "%s\n", methodName.c_str() ); (void) fflush( log );
   }   

   BrookPlatform platform( 32, "cpu", log );
   //ReferencePlatform platform;
   System system;
   LangevinIntegrator integrator( 0, 0.1, 0.01 );

   // int index, double charge, double radius, double depth

   NonbondedForce* nonbonded = new NonbondedForce();
   for (int i = 0; i < numberOfParticles; i++) {
       system.addParticle(1.0);
       nonbonded->addParticle(0, 1.5, 0);
   }
   vector<pair<int, int> > bonds;
   bonds.push_back(pair<int, int>(0, 1));
   bonds.push_back(pair<int, int>(1, 2));
   bonds.push_back(pair<int, int>(2, 3));
   bonds.push_back(pair<int, int>(3, 4));
   nonbonded->createExceptionsFromBonds(bonds, 0.0, 0.0);
   int first14, second14;
   for (int i = 0; i < nonbonded->getNumExceptions(); i++) {
       int particle1, particle2;
       double chargeProd, sigma, epsilon;
       nonbonded->getExceptionParameters(i, particle1, particle2, chargeProd, sigma, epsilon);
       if ((particle1 == 0 && particle2 == 3) || (particle1 == 3 && particle2 == 0))
           first14 = i;
       if ((particle1 == 1 && particle2 == 4) || (particle1 == 4 && particle2 == 1))
           second14 = i;
   }
   system.addForce(nonbonded);

   //(void) fprintf( log, "%s: Calling context\n",  methodName.c_str() );
   //(void) fflush( log );

   Context context(system, integrator, platform);
   vector<Vec3> positions(numberOfParticles);

   const double r = 1.0;
   positions[0] = Vec3(0, 0, 0);
   for( int ii = 1; ii < numberOfParticles; ii++ ){
      positions[ii] = Vec3(r, 0, 0);
   }
   for( int ii = 1; ii < numberOfParticles; ii++ ){

      // Test LJ forces

       vector<Vec3> positions(5);
       const double r = 1.0;
       for (int j = 0; j < 5; ++j) {
           nonbonded->setParticleParameters(j, 0, 1.5, 0); 
           positions[j] = Vec3(0, j, 0); 
       }
       nonbonded->setParticleParameters(0, 0, 1.5, 1); 
       nonbonded->setParticleParameters(ii, 0, 1.5, 1); 
       nonbonded->setExceptionParameters(first14, 0, 3, 0, 1.5, ii == 3 ? 0.5 : 0.0);
       nonbonded->setExceptionParameters(second14, 1, 4, 0, 1.5, 0.0);
       positions[ii] = Vec3(r, 0, 0); 

      context.reinitialize();
      context.setPositions(positions);

      State state = context.getState( State::Forces | State::Energy );
      const vector<Vec3>& forces = state.getForces();
      double x = 1.5/r;
      double eps = 1.0;
      double force = 4.0*eps*(12*std::pow(x, 12.0)-6*std::pow(x, 6.0))/r;
      double energy = 4.0*eps*(std::pow(x, 12.0)-std::pow(x, 6.0));
      if( ii == 3 ){
         force *= 0.5;
         energy *= 0.5;
      }
      if( ii < 3 ){
         force = 0;
         energy = 0;
      }

      if( PrintOn ){
         (void) fprintf( log, "14 LJ forces ii=%d F=%.6e\n", ii, force );
         for( int jj = 0; jj < numberOfParticles; jj++ ){
            (void) fprintf( log, "%d [%.5e %.5e %.5e]\n", jj, forces[jj][0], forces[jj][1], forces[jj][2] );
         }
         (void) fflush( log );
      }

      ASSERT_EQUAL_VEC(Vec3(-force, 0, 0), forces[0], TOL);
      ASSERT_EQUAL_VEC(Vec3(force, 0, 0), forces[ii], TOL);

      if( PrintOn ){
         (void) fprintf( log, "14 LJ forces ok for index=%d\n\n", ii );
         (void) fflush( log );
      }

      ASSERT_EQUAL_TOL(energy, state.getPotentialEnergy(), TOL);

      // Test Coulomb forces

      nonbonded->setParticleParameters( 0, 2, 1.5, 0 );
      nonbonded->setParticleParameters( ii, 2, 1.5, 0 );
      nonbonded->setExceptionParameters( first14, 0, 3, ii == 3 ? 4/1.2 : 0, 1.5, 0 );
      nonbonded->setExceptionParameters( second14, 1, 4, 0, 1.5, 0 );

      context.reinitialize();

      context.setPositions(positions);

      state = context.getState( State::Forces | State::Energy );

      const vector<Vec3>& forces2 = state.getForces();

      force = 138.935485*4/(r*r);
      energy = 138.935485*4/r;
      if( ii == 3 ){
         force /= 1.2;
         energy /= 1.2;
      }
      if( ii < 3 ){
          force = 0;
          energy = 0;
      }

      if( PrintOn ){
         (void) fprintf( log, "14 Coulomb forces ii=%d F=%.6e\n", ii, force );
         for( int jj = 0; jj < numberOfParticles; jj++ ){
            (void) fprintf( log, "%d [%.5e %.5e %.5e]\n", jj, forces2[jj][0], forces2[jj][1], forces2[jj][2] );
         }
         (void) fflush( log );
      }

      ASSERT_EQUAL_VEC(Vec3(-force, 0, 0), forces2[0], TOL);
      ASSERT_EQUAL_VEC(Vec3(force, 0, 0), forces2[ii], TOL);

      if( PrintOn ){
         (void) fprintf( log, "14 Coulomb forces ok for index=%d\n\n", ii );
         (void) fflush( log );
      }
       ASSERT_EQUAL_TOL(energy, state.getPotentialEnergy(), TOL);
   }

   if( PrintOn ){
      (void) fprintf( log, "ExclusionsAnd14 ok\n"); fflush( log );
   }

}

int main( ){

// ---------------------------------------------------------------------------------------

   static const std::string methodName      = "testBrookNonBonded";
   FILE* log                                = stdout;

// ---------------------------------------------------------------------------------------

   (void) fflush( stdout );
   (void) fflush( stderr );
   try {
      testBrookCoulomb( log );
      testBrookLJ( log );
      testBrookExclusionsAnd14( log );
    } catch( const exception& e ){
      (void) fprintf( log, "Exception %s %.s\n", methodName.c_str(),  e.what() ); (void) fflush( log );
      return 1;
   }   
   (void) fprintf( log, "\n%s done\n", methodName.c_str() ); (void) fflush( log );

   return 0;
}
