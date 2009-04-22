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
 * This tests the Brook OBC force/energy
 */

#include <vector>
#include <fstream>
#include <iostream>

#include "../../../tests/AssertionUtilities.h"
#include "BrookPlatform.h"
#include "ReferencePlatform.h"
#include "openmm/OpenMMContext.h"
#include "openmm/GBSAOBCForce.h"
#include "openmm/System.h"
#include "openmm/LangevinIntegrator.h"
#include "../src/sfmt/SFMT.h"
#include "../../reference/src/SimTKUtilities/SimTKOpenMMRealType.h"

typedef std::vector<std::string> StringVector;
typedef StringVector::iterator StringVectorI;
typedef StringVector::const_iterator StringVectorCI;

using namespace OpenMM;
using namespace std;

const double TOL = 1e-5;

/**---------------------------------------------------------------------------------------

   Replacement of sorts for strtok() (static method) (Simbios)
   Used to parse parameter file lines

   Should be moved to Utilities file

   @param lineBuffer           string to tokenize
   @param delimiter            token delimter

   @return number of args; if return value equals maxTokens, then more tokens than allocated

   --------------------------------------------------------------------------------------- */

char* strsepLocal( char** lineBuffer, const char* delimiter ){

   // ---------------------------------------------------------------------------------------

   // static const std::string methodName = "\nTinkerParameterSet::strsep"

   char *s;
   const char *spanp;
   int c, sc;
   char *tok;

   // ---------------------------------------------------------------------------------------

   s = *lineBuffer;
   if( s == NULL ){
      return (NULL);
   }

   for( tok = s;; ){
      c = *s++;
      spanp = delimiter;
      do {
         if( (sc = *spanp++) == c ){
            if( c == 0 ){
               s = NULL;
            } else {
               s[-1] = 0;
            }
            *lineBuffer = s;
            return( tok );
         }
      } while( sc != 0 );
   }
}

/**---------------------------------------------------------------------------------------

   Tokenize a string (static method) (Simbios)

   @param lineBuffer           string to tokenize
   @param tokenArray           upon return vector of tokens
   @param delimiter            token delimter

   @return number of args

   --------------------------------------------------------------------------------------- */

int tokenizeString( char* lineBuffer, StringVector& tokenArray, const std::string delimiter ){

   // ---------------------------------------------------------------------------------------

   // static const std::string methodName = "\nSimTKOpenMMUtilities::tokenizeString";

   // ---------------------------------------------------------------------------------------

// (void) fprintf( stdout, "\nIn SimTKOpenMMUtilities::tokenizeString <%s>", lineBuffer );
// (void) fflush( stdout );

   char *ptr_c = NULL;

   for( ; (ptr_c = strsepLocal( &lineBuffer, delimiter.c_str() )) != NULL; ){
      if( *ptr_c ){
         tokenArray.push_back( std::string( ptr_c ) );
      }
   }

   return (int) tokenArray.size();
}

static OpenMMContext* testObcForceSetup( int numParticles, int brookContext, FILE* log ){

// ---------------------------------------------------------------------------------------

   static const std::string methodName      = "testObcForceSetup";
   
// ---------------------------------------------------------------------------------------

   Platform* platform;
   if( brookContext ){
      platform = new BrookPlatform( 32, "cal", log );
   } else {
      platform = new ReferencePlatform();
   }

   System* system                 = new System();
   for (int i = 0; i < numParticles; i++)
       system->addParticle(1.0);
   LangevinIntegrator* integrator = new LangevinIntegrator(0, 0.1, 0.01);
   GBSAOBCForce* forceField  = new GBSAOBCForce();

   for (int i = 0; i < numParticles; ++i){
       // charge radius scalingFactor
       forceField->addParticle(i%2 == 0 ? -1 : 1, 0.15, 1);
       //forceField->setParticleParameters(i, i%2 == 0 ? -1 : 1, 1.5, 1);
   }
   system->addForce(forceField);
   OpenMMContext* context = new OpenMMContext( *system, *integrator, *platform );
   
   // Set random positions for all the atoms.
   
   vector<Vec3> positions(numParticles);
   init_gen_rand(0);
   for (int i = 0; i < numParticles; ++i){
       positions[i] = Vec3(1.0*genrand_real2(), 1.0*genrand_real2(), 1.0*genrand_real2());
   }
   context->setPositions(positions);

   return context;
}

/**---------------------------------------------------------------------------------------

   Replacement of sorts for strtok() (static method) (Simbios)
   Used to parse parameter file lines

   Should be moved to Utilities file

   @param lineBuffer           string to tokenize
   @param delimiter            token delimter

   @return number of args; if return value equals maxTokens, then more tokens than allocated

   --------------------------------------------------------------------------------------- */

static char* localStrsep( char** lineBuffer, const char* delimiter ){

   // ---------------------------------------------------------------------------------------

   // static const std::string methodName = "\nTinkerParameterSet::strsep"

   char *s;
   const char *spanp;
   int c, sc;
   char *tok;

   // ---------------------------------------------------------------------------------------

   s = *lineBuffer;
   if( s == NULL ){
      return (NULL);
   }

   for( tok = s;; ){
      c = *s++;
      spanp = delimiter;
      do {
         if( (sc = *spanp++) == c ){
            if( c == 0 ){
               s = NULL;
            } else {
               s[-1] = 0;
            }
            *lineBuffer = s;
            return( tok );
         }
      } while( sc != 0 );
   }
}

static OpenMMContext* testObcForceFileSetup( std::string fileName, int brookContext, int* numParticles, FILE* log ){

// ---------------------------------------------------------------------------------------

   static const std::string methodName      = "testObcForceFileSetup";
   const int bufferSize                     = 1024;
   char buffer[bufferSize];

// ---------------------------------------------------------------------------------------

   Platform* platform;
   if( brookContext ){
      platform = new BrookPlatform( 32, "cal", log );
   } else {
      platform = new ReferencePlatform();
   }

   LangevinIntegrator* integrator = new LangevinIntegrator(0, 0.1, 0.01);
   FILE* readFile                 = fopen( fileName.c_str(), "r" );
   if( log ){
      if( !readFile ){
         (void) fprintf( log, "\n%s File=%s not opened.", methodName.c_str(), fileName.c_str() );
         (void) fflush( log );
         return NULL;
      } else {
         (void) fprintf( log, "\n%s File=%s opened.", methodName.c_str(), fileName.c_str() );
         (void) fflush( log );
      }
   }

   // atom count

   int lineCount         = 0;
   std::string delimiter = " ";
   (void) fgets( buffer, bufferSize, readFile );
   StringVector tokens;
   tokenizeString( buffer, tokens, delimiter );
   if( tokens.size() < 1 && log ){
      std::stringstream message;
      message << "\nProblem w/ line=" << lineCount << " <" << buffer << ">";
      (void) fprintf( log, "\n%s", message.str().c_str() );
      return NULL;
   }

   StringVectorI ii        = tokens.begin();
   int numberOfParticles   = atoi( (*ii).c_str() );

   if( log ){
      (void) fprintf( log, "\n%s %d\n", methodName.c_str(), numberOfParticles );
   }
   *numParticles           = numberOfParticles;
   lineCount++;

   System* system                 = new System();
   for (int i = 0; i < numberOfParticles; i++)
       system->addParticle(1.0);
   GBSAOBCForce* forceField       = new GBSAOBCForce();

   vector<Vec3> positions(numberOfParticles);
   int index                      = 0;
   for (int i = 0; i < numberOfParticles; ++i){

      (void) fgets( buffer, bufferSize, readFile );
      StringVector tokens;
      tokenizeString( buffer, tokens, delimiter );
      if( tokens.size() < 8 && log ){
         std::stringstream message;
         message << "\nProblem w/ line=" << lineCount << " <" << buffer << ">";
         (void) fprintf( log, "\n%s", message.str().c_str() );
         return NULL;
      }

      StringVectorI ii          = tokens.begin();
      ii++;
      float coordX              = (float) atof( (*ii).c_str() ); ii++;
      float coordY              = (float) atof( (*ii).c_str() ); ii++;
      float coordZ              = (float) atof( (*ii).c_str() ); ii++;
      float radius              = (float) atof( (*ii).c_str() ); ii++;
      float scalingFactor       = (float) atof( (*ii).c_str() ); ii++;
      float charge              = (float) atof( (*ii).c_str() ); ii++;
      float bornRadi            = (float) atof( (*ii).c_str() ); ii++;
   
      positions[index++]        = Vec3( coordX, coordY, coordZ );
      // charge radius scalingFactor
      forceField->addParticle( charge, radius, scalingFactor );

      if( log ){
         (void) fprintf( log, "%d [%.6f %.6f %.6f] q=%.6f rad=%.6f scl=%.6f bR=%.6f\n", i,
                         coordX, coordY, coordZ, charge, radius, scalingFactor, bornRadi );
      }
      lineCount++;
   }
   if( log ){
      (void) fflush( log );
   }
   system->addForce(forceField);
   OpenMMContext* context = new OpenMMContext( *system, *integrator, *platform );

   // Set positions for all the atoms.
   
   context->setPositions(positions);

   return context;
}

void testObcForce( FILE* log, char* testInputFileName ){

// ---------------------------------------------------------------------------------------

   static const std::string methodName      = "testObcForce";
   static const int maxPrint                = 10;
   static int PrintOn                       = 0;
   float epsilon                            = 1.0e-04f;
   int totalErrors                          = 0;
   int numParticles                         = 10;
   
// ---------------------------------------------------------------------------------------

   PrintOn = log ? PrintOn : 0;
   log     = PrintOn ? log : NULL;

   //OpenMMContext* context      = testObcForceSetup( numParticles, 0, log );
   //OpenMMContext* brookContext = testObcForceSetup( numParticles, 1, log );
   
   OpenMMContext* context      = testObcForceFileSetup( std::string( testInputFileName ), 0, &numParticles, log );
   //OpenMMContext* context      = NULL;
   OpenMMContext* brookContext = testObcForceFileSetup( std::string( testInputFileName ), 1, &numParticles, log );
   //OpenMMContext* brookContext = NULL;
   
   vector<Vec3> forces;
   if( context ){
      State state      = context->getState( State::Forces );
      forces           = state.getForces();
   }

   vector<Vec3> brookForces;
   if( brookContext ){
      State brookState   = brookContext->getState( State::Forces );
      brookForces        = brookState.getForces();
   }
   
   if( PrintOn > 1 ){
      (void) fprintf( log, "%s OBC forces [%s %s]\n", methodName.c_str(), (context?"Ref":""), (brookContext?"Brook":"") );
      for( int ii = 0; ii < numParticles; ii++ ){
         (void) fprintf( log, "%4d ", ii );
         if( context && brookContext ){
            double diff[3];
            double rdiff[3];
            int hit = 0;
            for( int jj = 0; jj < 3; jj++ ){
               diff[jj] = fabs( forces[ii][jj] -  brookForces[ii][jj] );
               if( forces[ii][jj] != 0.0 ){
                  rdiff[jj] = fabs( diff[jj]/forces[ii][jj] );
                  if( rdiff[jj] > 0.01 ){
                     hit++;
                  }
               } else {
                  rdiff[jj] = diff[jj];
               }
            }
            (void) fprintf( log, "diff[%8.3f %8.3f %8.3f] [%8.3f %8.3f %8.3f] ", rdiff[0], rdiff[1], rdiff[2], diff[0], diff[1], diff[2] );
            if( hit ){
               (void) fprintf( log, " XXX" );
            }
         }
         if( context ){
            (void) fprintf( log, "[%14.5e %14.5e %14.5e] ", forces[ii][0], forces[ii][1], forces[ii][2] );
         }
         if( brookContext ){
            (void) fprintf( log, "[%14.5e %14.5e %14.5e]", brookForces[ii][0], brookForces[ii][1], brookForces[ii][2] );
         }
         (void) fprintf( log, "\n" );
      }
      (void) fflush( log );
   }

   double tolerance = 1.0e-03;
   if( context && brookContext ){
      for (int i = 0; i < numParticles; ++i) {
          Vec3 f         = forces[i];
          Vec3 fBrook    = brookForces[i];
          ASSERT_EQUAL_VEC( f, fBrook, tolerance );
      }
   }

   if( PrintOn ){
      (void) fprintf( log, "testObcForce ok w/ tolerance=%.3e\n", tolerance );
      (void) fflush( log );
   }

}

void testObcSingleParticle( FILE* log ){

// ---------------------------------------------------------------------------------------

   int PrintOn                              = 0;
   static const std::string methodName      = "testObcSingleParticle";
   
// ---------------------------------------------------------------------------------------

   PrintOn = log ? PrintOn : 0;

   BrookPlatform platform( 32, "cal", log );
   //ReferencePlatform platform;

   System system;
   system.addParticle(2.0);
   LangevinIntegrator integrator(0, 0.1, 0.01);
   GBSAOBCForce* forceField = new GBSAOBCForce();
   forceField->addParticle(0.5, 0.15, 1);
   system.addForce(forceField);
   OpenMMContext context(system, integrator, platform);
   vector<Vec3> positions(1);
   positions[0] = Vec3(0, 0, 0);
   context.setPositions(positions);
   State state = context.getState(State::Energy);
   double bornRadius = 0.15-0.009; // dielectric offset
   double eps0 = EPSILON0;
   double bornEnergy = (-0.5*0.5/(8*PI_M*eps0))*(1.0/forceField->getSoluteDielectric()-1.0/forceField->getSolventDielectric())/bornRadius;
   double extendedRadius = bornRadius+0.14; // probe radius
   double nonpolarEnergy = CAL2JOULE*PI_M*0.0216*(10*extendedRadius)*(10*extendedRadius)*std::pow(0.15/bornRadius, 6.0); // Where did this formula come from?  Just copied it from CpuImplicitSolvent.cpp
   ASSERT_EQUAL_TOL((bornEnergy+nonpolarEnergy), state.getPotentialEnergy(), 0.01);

   if( PrintOn ){
      (void) fprintf( log, "%s ok\n", methodName.c_str() );
      (void) fflush( log );
   }
}

void testObcEConsistentForce( FILE* log ){

// ---------------------------------------------------------------------------------------

   static const std::string methodName      = "testObcEConsistentForce";
   static int PrintOn                       = 0;
   
// ---------------------------------------------------------------------------------------

   BrookPlatform platform( 32, "cal", log );

   //ReferencePlatform platform;
   const int numParticles = 10; 
   System system; 
   LangevinIntegrator integrator(0, 0.1, 0.01);
   GBSAOBCForce* forceField = new GBSAOBCForce();
   for (int i = 0; i < numParticles; ++i) {
       system.addParticle(1.0);
       forceField->addParticle(i%2 == 0 ? -1 : 1, 0.15, 1);
   }
   system.addForce(forceField);
   OpenMMContext context(system, integrator, platform);
   
   // Set random positions for all the atoms.
   
   vector<Vec3> positions(numParticles);
   init_gen_rand(0);
   for (int i = 0; i < numParticles; ++i)
       positions[i] = Vec3(5.0*genrand_real2(), 5.0*genrand_real2(), 5.0*genrand_real2());
   context.setPositions(positions);
   State state = context.getState(State::Forces | State::Energy);
   
   // Take a small step in the direction of the energy gradient.
   
   double norm = 0.0;
   for (int i = 0; i < numParticles; ++i) {
       Vec3 f = state.getForces()[i];
       norm += f[0]*f[0] + f[1]*f[1] + f[2]*f[2];
   }   
   norm = std::sqrt(norm);
   const double delta = 1e-3;
   double step = delta/norm;
   for (int i = 0; i < numParticles; ++i) {
       Vec3 p = positions[i];
       Vec3 f = state.getForces()[i];
       positions[i] = Vec3(p[0]-f[0]*step, p[1]-f[1]*step, p[2]-f[2]*step);
   }   
   context.setPositions(positions);
   
   // See whether the potential energy changed by the expected amount.
   
   State state2 = context.getState(State::Energy);
   ASSERT_EQUAL_TOL(norm, (state2.getPotentialEnergy()-state.getPotentialEnergy())/delta, 0.01)

   if( PrintOn ){
      (void) fprintf( log, "%s ok %.8e %.8e    %.8e %.8e\n", methodName.c_str(),
                      state2.getPotentialEnergy(), state.getPotentialEnergy(), norm,
                      (state2.getPotentialEnergy()-state.getPotentialEnergy())/delta ); (void) fflush( log );
    
   }

}

int main( int argc, char* argv[] ){

// ---------------------------------------------------------------------------------------

   static const std::string methodName      = "testBrookObcGbsa";
   FILE* log                                = stdout;

// ---------------------------------------------------------------------------------------

   (void) fflush( stdout ); (void) fflush( stderr );

   try {

      // skip test if input file name is not provided
      // but let user know

      if( argc > 1 ){
         testObcForce( log, argv[1] );
      } else {
         (void) fprintf( log, "%s missing input file (=ObcInfo.txt?): skipping testObcForce() argc=%d\n", methodName.c_str(), argc );
         for( int ii = 0; ii < argc; ii++ ){
            (void) fprintf( log, "%d <%s>\n", ii, argv[ii] );
         }
         (void) fflush( log );
      }

      testObcSingleParticle( log );
      testObcEConsistentForce( log );

   } catch( const exception& e ){
      (void) fprintf( log, "Exception %s %.s\n", methodName.c_str(),  e.what() ); (void) fflush( log );
      return 1;
   }   
   (void) fprintf( log, "\n%s done\n", methodName.c_str() ); (void) fflush( log );

   return 0;
}
