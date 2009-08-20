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
 * Tests:
 *    (1) the relative differences between the Cuda and Reference forces agree to within specified tolerance
 *    (2) energy and forces are consistent
 *    (3) energy conservation (Verlet)/thermal stability (Langevin)
 * 
 */


#include "../../../tests/AssertionUtilities.h"
#include "CudaPlatform.h"
#include "ReferencePlatform.h"

//#define FREE_ENERGY_BRANCH

#ifdef FREE_ENERGY_BRANCH
#include "openmm/OpenMMContext.h"
#else
#include "openmm/Context.h"
#endif

#include "openmm/HarmonicBondForce.h"
#include "openmm/HarmonicAngleForce.h"
#include "openmm/PeriodicTorsionForce.h"
#include "openmm/RBTorsionForce.h"
#include "openmm/GBSAOBCForce.h"
#include "openmm/NonbondedForce.h"
#include "openmm/CMMotionRemover.h"
#include "openmm/System.h"
#include "openmm/LangevinIntegrator.h"
#include "openmm/VariableLangevinIntegrator.h"
#include "openmm/VerletIntegrator.h"
#include "openmm/VariableVerletIntegrator.h"
#include "openmm/BrownianIntegrator.h"
#include "../src/sfmt/SFMT.h"

#include <ctime>
#include <vector>
#include <cfloat>

// force enums

#define MAX_PRINT 5
#define HARMONIC_BOND_FORCE                      1
#define HARMONIC_ANGLE_FORCE                     2
#define PERIODIC_TORSION_FORCE                   4
#define RB_TORSION_FORCE                         8
#define NB_FORCE                                 16
#define NB_EXCEPTION_FORCE                       32
#define GBSA_OBC_FORCE                           64

#define BOLTZMANN  (1.380658e-23)                /* (J/K) */
#define AVOGADRO   (6.0221367e23)               /* ()    */
#define RGAS       (BOLTZMANN*AVOGADRO)         /* (J/(mol K))  */
#define BOLTZ      (RGAS/1.0e+03)               /* (kJ/(mol K)) */

using namespace OpenMM;
using namespace std;

// the following are used in parsing parameter file

typedef std::vector<std::string> StringVector;
typedef StringVector::iterator StringVectorI;
typedef StringVector::const_iterator StringVectorCI;

// default return value from methods

static const int DefaultReturnValue               = 0;

/**---------------------------------------------------------------------------------------

   Find stats for vec3

   @param array                 array 
   @param statVector            vector of stats 

   @return 0

   --------------------------------------------------------------------------------------- */

static int findStatsForVec3( const std::vector<Vec3>& array, std::vector<double>& statVector ){

   // ---------------------------------------------------------------------------------------
   
   static const int STAT_AVG = 0;
   static const int STAT_STD = 1;
   static const int STAT_MIN = 2;
   static const int STAT_ID1 = 3;
   static const int STAT_MAX = 4;
   static const int STAT_ID2 = 5;
   static const int STAT_CNT = 6;

   //static const char* methodName  = "\nfindStatsForVec3: ";

   // ---------------------------------------------------------------------------------------

   statVector.resize( STAT_CNT + 1 );

   double avgValue   =  0.0;
   double stdValue   =  0.0;
   double minValue   =  1.0e+30;
   double maxValue   = -1.0e+30;
   int minValueIndex = 0;
   int maxValueIndex = 0;

   for( unsigned int ii = 0; ii < array.size(); ii++ ){

	   double norm2 = array[ii][0]*array[ii][0] + array[ii][1]*array[ii][1] + array[ii][2]*array[ii][2];
	   double norm  = std::sqrt( norm2 );

      avgValue    += norm;
      stdValue    += norm2;

      if( norm > maxValue ){
         maxValue       = norm;
         maxValueIndex  = ii;
      }
      if( norm < minValue ){
         minValue       = norm;
         minValueIndex  = ii;
      }
   }

   double count  = static_cast<double>(array.size());
   double iCount = count > 0.0 ? 1.0/count : 0.0;
  
   statVector[STAT_AVG] = avgValue*iCount;
   statVector[STAT_STD] = stdValue - avgValue*avgValue*count;
   if( count > 1.0 ){
      statVector[STAT_STD] = std::sqrt( stdValue/( count - 1.0 ) );
   }
   statVector[STAT_MIN] = minValue;
   statVector[STAT_ID1] = static_cast<double>(minValueIndex);
   statVector[STAT_MAX] = maxValue;
   statVector[STAT_ID2] = static_cast<double>(maxValueIndex);
   statVector[STAT_CNT] = count;

   return DefaultReturnValue;
}

/* ---------------------------------------------------------------------------------------

   Compute cross product of two 3-vectors and place in 3rd vector  -- helper method

   vectorZ = vectorX x vectorY

   @param vectorX             x-vector
   @param vectorY             y-vector
   @param vectorZ             z-vector

   @return vector is vectorZ

   --------------------------------------------------------------------------------------- */
     
void crossProductVector3D( double* vectorX, double* vectorY, double* vectorZ ){

   // ---------------------------------------------------------------------------------------

   // static const char* methodName = "crossProductVector3D";

   // ---------------------------------------------------------------------------------------

   vectorZ[0]  = vectorX[1]*vectorY[2] - vectorX[2]*vectorY[1];
   vectorZ[1]  = vectorX[2]*vectorY[0] - vectorX[0]*vectorY[2];
   vectorZ[2]  = vectorX[0]*vectorY[1] - vectorX[1]*vectorY[0];

   return;
}

/* ---------------------------------------------------------------------------------------

   Compute cross product of two 3-vectors and place in 3rd vector  -- helper method

   vectorZ = vectorX x vectorY

   @param vectorX             x-vector
   @param vectorY             y-vector
   @param vectorZ             z-vector

   @return vector is vectorZ

   --------------------------------------------------------------------------------------- */
     
void crossProductVector3F( float* vectorX, float* vectorY, float* vectorZ ){

   // ---------------------------------------------------------------------------------------

   // static const char* methodName = "crossProductVector3D";

   // ---------------------------------------------------------------------------------------

   vectorZ[0]  = vectorX[1]*vectorY[2] - vectorX[2]*vectorY[1];
   vectorZ[1]  = vectorX[2]*vectorY[0] - vectorX[0]*vectorY[2];
   vectorZ[2]  = vectorX[0]*vectorY[1] - vectorX[1]*vectorY[0];

   return;
}

/* ---------------------------------------------------------------------------------------

   Return nonzero if all entries in array targets match all entries in array bond (order unimportant)

   @param numberIndices       number of entries in array
   @param targets             array of numberIndices ints
   @param bond                array of numberIndices ints

   @return nonzero if all entries in targets match all entries in bond (order unimportant)

   --------------------------------------------------------------------------------------- */
     
static int checkBondIndices( int numberIndices, const int* targets, const int* bond ){

   // ---------------------------------------------------------------------------------------

   // static const char* methodName = "checkBondIndices";

   // ---------------------------------------------------------------------------------------

   for( int ii = 0; ii < numberIndices; ii++ ){
      int hit = 0;
      for( int jj = 0; jj < numberIndices && hit == 0; jj++ ){
         if( targets[ii] == bond[jj] )hit = 1;
      }
      if( hit == 0 )return 0;
   }

   return 1;
}

/**---------------------------------------------------------------------------------------

   Find stats for vec3

   @param array                 array 
   @param statVector            vector of stats 

   @return 0

   --------------------------------------------------------------------------------------- */

static int angleTestCalculate( double* vector1, double* vector2, FILE* log ){

   // ---------------------------------------------------------------------------------------
   
   //static const char* methodName  = "\nangleTest: ";

   // ---------------------------------------------------------------------------------------

   double crossProduct[3];

   float  vector1F[3];
   float  vector2F[3];
   float  crossProductF[3];

   for( int ii = 0; ii < 3; ii++ ){
      vector1F[ii] = static_cast<float>(vector1[ii]);
      vector2F[ii] = static_cast<float>(vector2[ii]);
   }

#define DOT3(u,v) ((u[0])*(v[0]) + (u[1])*(v[1]) + (u[2])*(v[2]))

   double dotProductD    = DOT3( vector1,  vector2 );
   double norm1D         = DOT3( vector1,  vector1 );
   double norm2D         = DOT3( vector2,  vector2 );
   dotProductD          /= sqrt( norm1D*norm2D );
   dotProductD           = dotProductD <  1.0 ? dotProductD :  1.0;
   dotProductD           = dotProductD > -1.0 ? dotProductD : -1.0;

   crossProductVector3D( vector1, vector2, crossProduct );
   double normCrossD     = DOT3( crossProduct, crossProduct );
   normCrossD           /= sqrt( norm1D*norm2D );

   //(void) fprintf( log, "D: dot=%14.7e norms=%14.7e %14.7e cross=%14.7e\n", dotProductD, norm1D, norm2D, normCrossD );

   double angleCosD      = acos( dotProductD );
   double angleSinD      = asin( normCrossD );

   // ----------------------------------------------------------------------------

   float  dotProductF    = DOT3( vector1F, vector2F );
   float  norm1F         = DOT3( vector1F, vector1F );
   float  norm2F         = DOT3( vector2F, vector2F );
   dotProductF          /= sqrt( norm1F*norm2F );
   dotProductF           = dotProductF <  1.0f ? dotProductF :  1.0f;
   dotProductF           = dotProductF > -1.0f ? dotProductF : -1.0f;
   crossProductVector3F( vector1F, vector2F, crossProductF );
   float  normCrossF     = DOT3( crossProductF, crossProductF );
   normCrossF           /= sqrt( norm1F*norm2F );

   float angleCosF       = acosf( dotProductF );
   float angleSinF       = asinf( normCrossF );

   // ----------------------------------------------------------------------------

   double deltaAngleCos  = fabs( angleCosD - static_cast<float>(angleCosF) ); 
   double deltaAngleSin  = fabs( angleSinD - static_cast<float>(angleSinF) ); 

   (void) fprintf( log, "%14.7e %14.7e %14.7e %14.7e %14.7e     %14.7e %14.7e %14.7e %14.7e %14.7e\n",
                   deltaAngleCos, dotProductD, dotProductF, angleCosD, angleCosF,
                   deltaAngleSin, normCrossD,  normCrossF,  angleSinD, angleSinF );

   return DefaultReturnValue;
}

/**---------------------------------------------------------------------------------------

   Find stats for vec3

   @param array                 array 
   @param statVector            vector of stats 

   @return 0

   --------------------------------------------------------------------------------------- */

static int angleTest( FILE* log ){

   // ---------------------------------------------------------------------------------------
   
   //static const char* methodName  = "\nangleTest: ";

   // ---------------------------------------------------------------------------------------

   double vector1[3];
   double vector2[3];

   double tempVector1[3];
   double tempVector2[3];

/*
Bpti atom 319 sdObc
ReferenceRbDihedralBond::calculateBondIxn

 Atm 327 [-0.404 0.604  -0.415]  Atm 318 [-0.358 0.487  -0.358]  Atm 319 [-0.299 0.391  -0.439]  Atm 320 [-0.262 0.3  -0.396] 
 Delta: [-0.046 0.117 -0.057 0.019054 0.138036 ] [0.059 -0.096 -0.081 0.019258 0.138773 ] [-0.037 0.091 -0.043 0.011499 0.107233 ]

 Cross: [-0.014949 -0.007089 -0.002487 ] [0.011499 0.005534 0.001817 ]
 k=30.334 a=0 m=-30.334 ang=-0.00962353 dotD=0.999954 sign=1
   dEdAngle=-0.583804 E=0.00280952 force factors: [289.436 -0.484422 -0.386125 -487.599 ] F=compute force; f=cumulative force
   F1[-4.32677 -2.05181 -0.719827 ] F2[-4.25779 -2.00384 -0.726432 ] F3[-5.67588 -2.74634 -0.879363 ] F4[-5.6069 -2.69837 -0.885968 ]
   f1[-4.32677 -2.05181 -0.719827 ] f2[26.0422 -32.83 -2.51618 ] f3[6.17743 2.98757 0.95879 ] f4[-5.78133 -2.78232 -0.91353 ]
*/

   vector1[0] = -0.014949;
   vector1[1] = -0.007089;
   vector1[2] = -0.002487;

   vector2[0] =  0.011499;
   vector2[1] =  0.005534;
   vector2[2] =  0.001817;

   vector1[0] = -1.0;
   vector1[1] =  0.0;
   vector1[2] =  0.0;

   vector2[0] =  0.0;
   vector2[1] =  0.0;
   vector2[2] =  1.0;

   double dotProductD    = DOT3( vector1,  vector2 );
   double norm1D         = DOT3( vector1,  vector1 );
   double norm2D         = DOT3( vector2,  vector2 );
   double target         = -1.0;
   double alpha          = (target - dotProductD)/(norm1D);
   double offset         = 1.0e-03;

   for( int ii = 1; ii < 100; ii++ ){
      double tempAlpha = alpha*(1.0 + static_cast<double>(ii)*offset );
      for( int jj = 0; jj < 3; jj++ ){
         tempVector1[jj] = vector1[jj];
         //tempVector2[jj] = vector2[jj] + vector1[jj]*tempAlpha;
         tempVector2[jj] = vector2[jj];
      }
      tempVector2[0] = offset/static_cast<double>(ii);
      
      angleTestCalculate( tempVector1, tempVector2, log );
   }

   return DefaultReturnValue;
}

/**---------------------------------------------------------------------------------------

   Find stats for double array

   @param array                   array 
   @param statVector              vector of stats 

   @return 0

   --------------------------------------------------------------------------------------- */

static int findStatsForDouble( const std::vector<double>& array, std::vector<double>& statVector ){

   // ---------------------------------------------------------------------------------------
   
   static const int STAT_AVG = 0;
   static const int STAT_STD = 1;
   static const int STAT_MIN = 2;
   static const int STAT_ID1 = 3;
   static const int STAT_MAX = 4;
   static const int STAT_ID2 = 5;
   static const int STAT_CNT = 6;

   //static const char* methodName  = "\nfindStatsForDouble: ";

   // ---------------------------------------------------------------------------------------

   statVector.resize( STAT_CNT + 1 );

   double avgValue   =  0.0;
   double stdValue   =  0.0;
   double minValue   =  1.0e+30;
   double maxValue   = -1.0e+30;
   int minValueIndex = 0;
   int maxValueIndex = 0;

   for( unsigned int ii = 0; ii < array.size(); ii++ ){

	   double norm  =  array[ii];

      avgValue    += norm;
      stdValue    += norm*norm;

      if( norm > maxValue ){
         maxValue       = norm;
         maxValueIndex  = ii;
      }
      if( norm < minValue ){
         minValue       = norm;
         minValueIndex  = ii;
      }
   }

   double count  = static_cast<double>(array.size());
   double iCount = count > 0.0 ? 1.0/count : 0.0;
  
   statVector[STAT_AVG] = avgValue*iCount;
   statVector[STAT_STD] = stdValue - avgValue*avgValue*count;
   if( count > 1.0 ){
      statVector[STAT_STD] = std::sqrt( stdValue/( count - 1.0 ) );
   }
   statVector[STAT_MIN] = minValue;
   statVector[STAT_ID1] = static_cast<double>(minValueIndex);
   statVector[STAT_MAX] = maxValue;
   statVector[STAT_ID2] = static_cast<double>(maxValueIndex);
   statVector[STAT_CNT] = count;

   return DefaultReturnValue;
}

/**---------------------------------------------------------------------------------------

 * Check constraints
 *
 * @param    context            OpenMM::Context used to get current positions
 * @param    system             OpenMM::System to be created
 * @param    tolerance          constraint tolerance
 * @param    log                log file
 *
 * @return   return number of violations

   --------------------------------------------------------------------------------------- */

#ifdef FREE_ENERGY_BRANCH
static int checkConstraints( const OpenMMContext& context, const System& system, double tolerance, FILE* log ) {
#else
static int checkConstraints( const Context& context, const System& system, double tolerance, FILE* log ) {
#endif
    
	 int violations                    = 0;
	 State state                       = context.getState(State::Positions);
	 const std::vector<Vec3>& pos      = state.getPositions();
	 for( int ii = 0; ii < system.getNumConstraints(); ii++ ){
	    int particle1;
	    int particle2;
	    double distance;
		 system.getConstraintParameters( ii, particle1, particle2, distance );
		 double actualDistance = sqrt(
 		                         (pos[particle2][0] - pos[particle1][0])*(pos[particle2][0] - pos[particle1][0]) +
 		                         (pos[particle2][1] - pos[particle1][1])*(pos[particle2][1] - pos[particle1][1]) +
 		                         (pos[particle2][2] - pos[particle1][2])*(pos[particle2][2] - pos[particle1][2]) );
       double delta          = fabs( actualDistance - distance );
		 if( delta > tolerance ){
		    violations++;
			 if( log ){
			    (void) fprintf( log, "CnstrViolation: %6d %6d [%6d %6d] %10.3e [%12.5e %12.5e] \n",
				                 ii, violations, particle1, particle2, delta, distance, actualDistance );
          }
       }
	 }
	 if( log ){
       (void) fprintf( log, "CnstrViolation: total violations=%d out of %d constraints; tolerance=%.3e.\n",
		                 violations, system.getNumConstraints(), tolerance );
    }

	 return violations;
}

/**---------------------------------------------------------------------------------------

   Sum forces

      @param context                  OpenMM::Context used to get current positions
      @param system                   OpenMM::System to be created
      @param forceSum                 on return, sum of forces
      @param step                     step index
      @param log                      log reference (stdlog in md.c)

      @return DefaultReturnValue

      --------------------------------------------------------------------------------------- */

#ifdef FREE_ENERGY_BRANCH
static int sumForces( const OpenMMContext& context, const System& system, double forceSum[3], int step, FILE* log ){  
#else
static int sumForces( const Context& context, const System& system, double forceSum[3], int step, FILE* log ){  
#endif

// ---------------------------------------------------------------------------------------
      
   double sum;

// ---------------------------------------------------------------------------------------

	 State state                       = context.getState(State::Forces);
	 const std::vector<Vec3>& forces   = state.getForces();

   // sum forces and track max value

   forceSum[0]       = forceSum[1] = forceSum[2] = 0.0;
	double forceMax   = -1.0;
	int forceMaxIndex = -1;
	for( int ii = 0; ii < system.getNumParticles(); ii++ ){
      forceSum[0]             += forces[ii][0];
      forceSum[1]             += forces[ii][1];
      forceSum[2]             += forces[ii][2];
		double forceMagnitude    = forces[ii][0]*forces[ii][0] + forces[ii][1]*forces[ii][1] + forces[ii][2]*forces[ii][2];
		if( forceMagnitude > forceMax ){
		   forceMax      = forceMagnitude;
			forceMaxIndex = ii;
      }
   }   

   if( 0 ){
      sum = fabs( forceSum[0] ) + fabs( forceSum[1] ) + fabs( forceSum[2] );
      (void) fprintf( log, "Force: Step=%d %.4e f[%.4e %.4e %.4e] Max=%.3e at index=%d\n", step, sum,
                      forceSum[0], forceSum[1], forceSum[2], sqrt( forceMax ), forceMaxIndex );
   }   

   return 0;
}

/**---------------------------------------------------------------------------------------

   Check kinetic energy

   @param numberOfAtoms            number of atoms
   @param nrdf                     number of degrees of freedom
   @param v                        velocities
   @param mass                     masses
   @param temperature              temperature
   @param step                     step index
   @param log                      log reference (stdlog in md.c)

   @return DefaultReturnValue if k.e. ~ temp;

   --------------------------------------------------------------------------------------- */

#ifdef FREE_ENERGY_BRANCH
static int checkKineticEnergy( const OpenMMContext& context, const System& system, double temperature, int step, FILE* log ){  
#else
static int checkKineticEnergy( const Context& context, const System& system, double temperature, int step, FILE* log ){  
#endif

   // ---------------------------------------------------------------------------------------

   double kineticEnergy;

   int status            = 0;
   int print             = 1;
   double cutoff         = 200.0;
   static double average = 0.0;
   static double stddev  = 0.0;
   static double count   = 0.0;

   // ---------------------------------------------------------------------------------------

   // calculate kineticEnergy

   State state           = context.getState(State::Energy);
	kineticEnergy         = 2.0*state.getKineticEnergy();

   int nrdf              = 3*system.getNumParticles() - system.getNumConstraints() - 3;
   kineticEnergy        /= (((double) BOLTZ)*((double) nrdf));
   if( print ){
      double averageL, stddevL;
      average  += kineticEnergy;
      stddev   += kineticEnergy*kineticEnergy;
      count    += 1.0;
      averageL  = average/count;
      stddevL   = stddev - averageL*averageL*count;
      if( stddevL > 0.0 && count > 1 ){
         stddevL = sqrt( stddevL/(count - 1.0 ) );
      }
      (void) (void) fprintf( log, "checkKineticEnergy: Step=%d T=%.3f avg=%g std=%.3f nrdf=%d\n", step, kineticEnergy, averageL, stddevL, nrdf);
   }

   // only check if calculated T is > specified T
/*
   if( (kineticEnergy - temperature) > cutoff ){
      (void) (void) fprintf( log, "checkKineticEnergy: ERROR Step=%d T=%.3f tpr-T=%.3f diff=%.3f cutoff=%.3f\n",
                      step, kineticEnergy, temperature, (kineticEnergy - temperature), cutoff );
       (void) fflush( NULL );
       status = 1;
    }
*/
    // ignore calculation prior to 2000 steps

    return 0;
}

/**---------------------------------------------------------------------------------------

   Replacement of sorts for strtok()
   Used to parse parameter file lines

   @param lineBuffer           string to tokenize
   @param delimiter            token delimter

   @return number of args

   --------------------------------------------------------------------------------------- */

char* strsepLocal( char** lineBuffer, const char* delimiter ){

   // ---------------------------------------------------------------------------------------

   // static const std::string methodName = "strsepLocal";

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
      c     = *s++;
      spanp = delimiter;
      do {
         if( (sc = *spanp++) == c ){
            if( c == 0 ){
               s = NULL;
            } else {
               s[-1] = 0;
            }
/*
            if( *s == '\n' ){ 
               *s = NULL;
            }
*/
            *lineBuffer = s;
            return( tok );
         }
      } while( sc != 0 );
   }
}

/**---------------------------------------------------------------------------------------

   Tokenize a string

   @param lineBuffer           string to tokenize
   @param tokenArray           upon return vector of tokens
   @param delimiter            token delimter

   @return number of tokens

   --------------------------------------------------------------------------------------- */

int tokenizeString( char* lineBuffer, StringVector& tokenArray, const std::string delimiter ){

   // ---------------------------------------------------------------------------------------

   // static const std::string methodName = "\nSimTKOpenMMUtilities::tokenizeString";

   // ---------------------------------------------------------------------------------------

   char *ptr_c = NULL;

   for( ; (ptr_c = strsepLocal( &lineBuffer, delimiter.c_str() )) != NULL; ){
      if( *ptr_c ){
/*
         char* endOfLine = ptr_c;
         while( endOfLine ){
printf( "%c", *endOfLine ); fflush( stdout );
            if( *endOfLine == '\n' )*endOfLine = '\0';
            endOfLine++;
         }  
*/
         tokenArray.push_back( std::string( ptr_c ) );
      }
   }

   return (int) tokenArray.size();
}

/**---------------------------------------------------------------------------------------

   Read a line from a file and tokenize into an array of strings

   @param filePtr              file to read from
   @param tokens               array of token strings
   @param lineCount            line count
   @param log                  optional file ptr for logging

   @return ptr to string containing line

   --------------------------------------------------------------------------------------- */

static char* readLine( FILE* filePtr, StringVector& tokens, int* lineCount, FILE* log ){

// ---------------------------------------------------------------------------------------

   //static const std::string methodName      = "readLine";
   
   std::string delimiter                    = " \n";
   const int bufferSize                     = 4096;
   char buffer[bufferSize];

// ---------------------------------------------------------------------------------------

   char* isNotEof = fgets( buffer, bufferSize, filePtr );
   if( isNotEof ){
      (*lineCount)++;
      tokenizeString( buffer, tokens, delimiter );
   }
   return isNotEof;

}

/**---------------------------------------------------------------------------------------

   Read particles count

   @param filePtr              file pointer to parameter file
   @param tokens               array of strings from first line of parameter file for this block of parameters
   @param system               System reference
   @param lineCount            used to track line entries read from parameter file
   @param log                  log file pointer -- may be NULL

   @return number of parameters read

   --------------------------------------------------------------------------------------- */

static int readParticles( FILE* filePtr, const StringVector& tokens, System& system, int* lineCount, FILE* log ){

// ---------------------------------------------------------------------------------------

   static const std::string methodName      = "readParticles";
   
// ---------------------------------------------------------------------------------------

   if( tokens.size() < 1 ){
      char buffer[1024];
      (void) sprintf( buffer, "%s no particles number entry???\n", methodName.c_str() );
      throwException(__FILE__, __LINE__, buffer );
      exit(-1);
   }

   int numberOfParticles = atoi( tokens[1].c_str() );
   if( log ){
      (void) fprintf( log, "%s particles=%d\n", methodName.c_str(), numberOfParticles );
   }

   return numberOfParticles;
}

/**---------------------------------------------------------------------------------------

   Read particle masses

   @param filePtr              file pointer to parameter file
   @param tokens               array of strings from first line of parameter file for this block of parameters
   @param system               System reference
   @param lineCount            used to track line entries read from parameter file
   @param log                  log file pointer -- may be NULL

   @return number of masses read

   --------------------------------------------------------------------------------------- */

static int readMasses( FILE* filePtr, const StringVector& tokens, System& system, int* lineCount, FILE* log ){

// ---------------------------------------------------------------------------------------

   static const std::string methodName      = "readMasses";
   
// ---------------------------------------------------------------------------------------

   if( tokens.size() < 1 ){
      char buffer[1024];
      (void) sprintf( buffer, "%s no particles number entry???\n", methodName.c_str() );
      throwException(__FILE__, __LINE__, buffer );
      exit(-1);
   }

   int numberOfParticles = atoi( tokens[1].c_str() );
   if( log ){
      (void) fprintf( log, "%s particle masses=%d\n", methodName.c_str(), numberOfParticles );
   }
   for( int ii = 0; ii < numberOfParticles; ii++ ){
      StringVector lineTokens;
      char* isNotEof = readLine( filePtr, lineTokens, lineCount, log );
      if( lineTokens.size() >= 1 ){
         int index   = atoi( lineTokens[0].c_str() );
         double mass = atof( lineTokens[1].c_str() );
         system.addParticle( mass );
      } else {
         char buffer[1024];
         (void) sprintf( buffer, "%s particle tokens incomplete at line=%d\n", methodName.c_str(), *lineCount );
         throwException(__FILE__, __LINE__, buffer );
         exit(-1);
      }
   }

   // diagnostics

   if( log ){
      static const unsigned int maxPrint   = MAX_PRINT;
      unsigned int arraySize               = static_cast<unsigned int>(system.getNumParticles());
      (void) fprintf( log, "%s: sample of masses\n", methodName.c_str() );
      for( unsigned int ii = 0; ii < arraySize && ii < maxPrint; ii++ ){
         (void) fprintf( log, "%6u %14.7e \n", ii, system.getParticleMass( ii ) );
      }
      if( arraySize > maxPrint ){
         for( unsigned int ii = arraySize - maxPrint; ii < arraySize; ii++ ){
            (void) fprintf( log, "%6u %14.7e \n", ii, system.getParticleMass( ii ) );
         }
      }
   }

   return system.getNumParticles();
}

/**---------------------------------------------------------------------------------------

   Read harmonic bond parameters

   @param filePtr              file pointer to parameter file
   @param forceFlag            flag signalling whether force is to be added to system
                               if force == 0 || forceFlag & HARMONIC_BOND_FORCE, then included
   @param tokens               array of strings from first line of parameter file for this block of parameters
   @param system               System reference
   @param lineCount            used to track line entries read from parameter file
   @param log                  log file pointer -- may be NULL

   @return number of bonds

   --------------------------------------------------------------------------------------- */

static int readHarmonicBondForce( FILE* filePtr, int forceFlag, const StringVector& tokens, System& system, int* lineCount, FILE* log ){

// ---------------------------------------------------------------------------------------

   static const std::string methodName      = "readHarmonicBondForce";
   
// ---------------------------------------------------------------------------------------

   if( tokens.size() < 1 ){
      char buffer[1024];
      (void) sprintf( buffer, "%s no bonds number entry???\n", methodName.c_str() );
      throwException(__FILE__, __LINE__, buffer );
      exit(-1);
   }

   HarmonicBondForce* bondForce = new HarmonicBondForce();
   if( forceFlag == 0 || forceFlag & HARMONIC_BOND_FORCE ){
      system.addForce( bondForce );
      if( log && forceFlag ){
         (void) fprintf( log, "harmonic bond force is being included.\n" );
      }
   } else if( forceFlag && log ){
      (void) fprintf( log, "harmonic bond force is not being included.\n" );
   }

   int numberOfBonds            = atoi( tokens[1].c_str() );
   if( log ){
      (void) fprintf( log, "%s number of HarmonicBondForce terms=%d\n", methodName.c_str(), numberOfBonds );
   }
   for( int ii = 0; ii < numberOfBonds; ii++ ){
      StringVector lineTokens;
      char* isNotEof = readLine( filePtr, lineTokens, lineCount, log );
      if( lineTokens.size() > 4 ){
         int index      = atoi( lineTokens[0].c_str() );
         int particle1  = atoi( lineTokens[1].c_str() );
         int particle2  = atoi( lineTokens[2].c_str() );
         double length  = atof( lineTokens[3].c_str() );
         double k       = atof( lineTokens[4].c_str() );
         bondForce->addBond( particle1, particle2, length, k );
      } else {
         (void) fprintf( log, "%s HarmonicBondForce tokens incomplete at line=%d\n", methodName.c_str(), *lineCount );
         exit(-1);
      }
   }

   // diagnostics

   if( log ){
      static const unsigned int maxPrint   = MAX_PRINT;
      unsigned int arraySize               = static_cast<unsigned int>(bondForce->getNumBonds());
      (void) fprintf( log, "%s: sample of HarmonicBondForce parameters\n", methodName.c_str() );
      for( unsigned int ii = 0; ii < arraySize && ii < maxPrint; ii++ ){
         int particle1, particle2;
         double length, k;
         bondForce->getBondParameters( ii, particle1, particle2, length, k ); 
         (void) fprintf( log, "%8d %8d %8d %14.7e %14.7e\n", ii, particle1, particle2, length, k);
      }
      if(  arraySize > maxPrint ){
         for( unsigned int ii = arraySize - maxPrint; ii < arraySize; ii++ ){
            int particle1, particle2;
            double length, k;
            bondForce->getBondParameters( ii, particle1, particle2, length, k ); 
            (void) fprintf( log, "%8d %8d %8d %14.7e %14.7e\n", ii, particle1, particle2, length, k);
         }
      }
   }

   return bondForce->getNumBonds();
}

/**---------------------------------------------------------------------------------------

   Read harmonic angle parameters

   @param filePtr              file pointer to parameter file
   @param forceFlag            flag signalling whether force is to be added to system
                               if force == 0 || forceFlag & HARMONIC_ANGLE_FORCE, then included
   @param tokens               array of strings from first line of parameter file for this block of parameters
   @param system               System reference
   @param lineCount            used to track line entries read from parameter file
   @param log                  log file pointer -- may be NULL

   @return number of bonds

   --------------------------------------------------------------------------------------- */

static int readHarmonicAngleForce( FILE* filePtr, int forceFlag, const StringVector& tokens,
                                   System& system, int* lineCount, FILE* log ){

// ---------------------------------------------------------------------------------------

   static const std::string methodName      = "readHarmonicAngleForce";
   
// ---------------------------------------------------------------------------------------

   if( tokens.size() < 1 ){
      char buffer[1024];
      (void) sprintf( buffer, "%s no angle bonds number entry???\n", methodName.c_str() );
      throwException(__FILE__, __LINE__, buffer );
   }

   HarmonicAngleForce* bondForce = new HarmonicAngleForce();
   if( forceFlag == 0 || forceFlag & HARMONIC_ANGLE_FORCE ){
      system.addForce( bondForce );
      if( log && forceFlag ){
         (void) fprintf( log, "harmonic angle force is being included.\n" );
      }
   } else if( forceFlag && log ){
      (void) fprintf( log, "harmonic angle force is not being included.\n" );
   }

   int numberOfAngles            = atoi( tokens[1].c_str() );
   if( log ){
      (void) fprintf( log, "%s number of HarmonicAngleForce terms=%d\n", methodName.c_str(), numberOfAngles );
   }
   for( int ii = 0; ii < numberOfAngles; ii++ ){
      StringVector lineTokens;
      char* isNotEof = readLine( filePtr, lineTokens, lineCount, log );
      if( lineTokens.size() > 5 ){
         int index      = atoi( lineTokens[0].c_str() );
         int particle1  = atoi( lineTokens[1].c_str() );
         int particle2  = atoi( lineTokens[2].c_str() );
         int particle3  = atoi( lineTokens[3].c_str() );
         double angle   = atof( lineTokens[4].c_str() );
         double k       = atof( lineTokens[5].c_str() );
         bondForce->addAngle( particle1, particle2, particle3, angle, k );
      } else {
         char buffer[1024];
         (void) sprintf( buffer, "%s HarmonicAngleForce tokens incomplete at line=%d\n", methodName.c_str(), *lineCount );
         throwException(__FILE__, __LINE__, buffer );
         exit(-1);
      }
   }

   // diagnostics

   if( log ){
      static const unsigned int maxPrint   = MAX_PRINT;
      unsigned int arraySize               = static_cast<unsigned int>(bondForce->getNumAngles());
      (void) fprintf( log, "%s: sample of HarmonicAngleForce parameters\n", methodName.c_str() );
      for( unsigned int ii = 0; ii < arraySize && ii < maxPrint; ii++ ){
         int particle1, particle2, particle3;
         double angle, k;
         bondForce->getAngleParameters( ii, particle1, particle2, particle3, angle, k ); 
         (void) fprintf( log, "%8d %8d %8d %8d %14.7e %14.7e\n", ii, particle1, particle2, particle3, angle, k);
      }
      if(  arraySize > maxPrint ){
         for( unsigned int ii = arraySize - maxPrint; ii < arraySize; ii++ ){
            int particle1, particle2, particle3;
            double angle, k;
            bondForce->getAngleParameters( ii, particle1, particle2, particle3, angle, k ); 
            (void) fprintf( log, "%8d %8d %8d %8d %14.7e %14.7e\n", ii, particle1, particle2, particle3, angle, k);
         }
      }
   }

   return bondForce->getNumAngles();
}

/**---------------------------------------------------------------------------------------

   Read PeriodicTorsionForce parameters

   @param filePtr              file pointer to parameter file
   @param forceFlag            flag signalling whether force is to be added to system
                               if force == 0 || forceFlag & PERIODIC_TORSION_FORCE, then included
   @param tokens               array of strings from first line of parameter file for this block of parameters
   @param system               System reference
   @param lineCount            used to track line entries read from parameter file
   @param log                  log file pointer -- may be NULL

   @return number of torsion bonds read

   --------------------------------------------------------------------------------------- */

static int readPeriodicTorsionForce( FILE* filePtr, int forceFlag, const StringVector& tokens, System& system, int* lineCount, FILE* log ){

// ---------------------------------------------------------------------------------------

   static const std::string methodName      = "readPeriodicTorsionForce";
   
// ---------------------------------------------------------------------------------------

   if( tokens.size() < 1 ){
      char buffer[1024];
      (void) sprintf( buffer, "%s no PeriodicTorsion bonds number entry???\n", methodName.c_str() );
      throwException(__FILE__, __LINE__, buffer );
      exit(-1);
   }

   PeriodicTorsionForce* bondForce = new PeriodicTorsionForce();
   if( forceFlag == 0 || forceFlag & PERIODIC_TORSION_FORCE ){
      system.addForce( bondForce );
      if( log && forceFlag ){
         (void) fprintf( log, "periodic torsion force is being included.\n" );
      }
   } else if( forceFlag && log ){
      (void) fprintf( log, "periodic torsion force is not being included.\n" );
   }

   int numberOfTorsions            = atoi( tokens[1].c_str() );
   if( log ){
      (void) fprintf( log, "%s number of PeriodicTorsionForce terms=%d\n", methodName.c_str(), numberOfTorsions );
   }
   for( int ii = 0; ii < numberOfTorsions; ii++ ){
      StringVector lineTokens;
      char* isNotEof = readLine( filePtr, lineTokens, lineCount, log );
      if( lineTokens.size() > 7 ){
         int index       = atoi( lineTokens[0].c_str() );
         int particle1   = atoi( lineTokens[1].c_str() );
         int particle2   = atoi( lineTokens[2].c_str() );
         int particle3   = atoi( lineTokens[3].c_str() );
         int particle4   = atoi( lineTokens[4].c_str() );
         int periodicity = atoi( lineTokens[5].c_str() );
         double phase    = atof( lineTokens[6].c_str() );
         double k        = atof( lineTokens[7].c_str() );
         bondForce->addTorsion( particle1, particle2, particle3, particle4, periodicity, phase, k );
      } else {
         char buffer[1024];
         (void) sprintf( buffer, "%s PeriodicTorsionForce tokens incomplete at line=%d\n", methodName.c_str(), *lineCount );
         throwException(__FILE__, __LINE__, buffer );
         exit(-1);
      }
   }

   // diagnostics

   if( log ){
      static const unsigned int maxPrint   = MAX_PRINT;
      unsigned int arraySize               = static_cast<unsigned int>(bondForce->getNumTorsions());
      (void) fprintf( log, "%s: sample of PeriodicTorsionForce parameters\n", methodName.c_str() );
      for( unsigned int ii = 0; ii < arraySize && ii < maxPrint; ii++ ){
         int particle1, particle2, particle3, particle4, periodicity;
         double phase, k;
         bondForce->getTorsionParameters( ii, particle1, particle2, particle3, particle4, periodicity, phase, k );
         (void) fprintf( log, "%8d %8d %8d %8d %8d %8d %14.7e %14.7e\n", ii, particle1, particle2, particle3, particle4, periodicity, phase, k );
      }
      if(  arraySize > maxPrint ){
         for( unsigned int ii = arraySize - maxPrint; ii < arraySize; ii++ ){
            int particle1, particle2, particle3, particle4, periodicity;
            double phase, k;
            bondForce->getTorsionParameters( ii, particle1, particle2, particle3, particle4, periodicity, phase, k );
            (void) fprintf( log, "%8d %8d %8d %8d %8d %8d %14.7e %14.7e\n", ii, particle1, particle2, particle3, particle4, periodicity, phase, k );
         }
      }
   }

   return bondForce->getNumTorsions();
}

/**---------------------------------------------------------------------------------------

   Read RBTorsionForce parameters

   @param filePtr              file pointer to parameter file
   @param forceFlag            flag signalling whether force is to be added to system
                               if force == 0 || forceFlag & RB_TORSION_FORCE, then included
   @param tokens               array of strings from first line of parameter file for this block of parameters
   @param system               System reference
   @param lineCount            used to track line entries read from parameter file
   @param log                  log file pointer -- may be NULL

   @return number of torsion bonds read

   --------------------------------------------------------------------------------------- */

static int readRBTorsionForce( FILE* filePtr, int forceFlag, const StringVector& tokens, System& system, int* lineCount, FILE* log ){

// ---------------------------------------------------------------------------------------

   static const std::string methodName      = "readRBTorsionForce";
   
// ---------------------------------------------------------------------------------------

   if( tokens.size() < 1 ){
      char buffer[1024];
      (void) sprintf( buffer, "%s no RBTorsion bonds number entry???\n", methodName.c_str() );
      throwException(__FILE__, __LINE__, buffer );
      exit(-1);
   }

   RBTorsionForce* bondForce       = new RBTorsionForce();
   if( forceFlag == 0 || forceFlag & RB_TORSION_FORCE ){
      system.addForce( bondForce );
      if( log && forceFlag ){
         (void) fprintf( log, "RB torsion force is being included.\n" );
      }
   } else if( forceFlag && log ){
      (void) fprintf( log, "RB torsion force is not being included.\n" );
   }

   int numberOfTorsions            = atoi( tokens[1].c_str() );
   if( log ){
      (void) fprintf( log, "%s number of RBTorsionForce terms=%d\n", methodName.c_str(), numberOfTorsions );
      (void) fflush( log );
   }

#if 1
static int nextIndex = 0;
static int parity    = 0;
int targets[12][4] = { 
                        { 315,318,319,320 },
                        { 315,318,319,321 },
                        { 327,318,319,320 },
                        { 327,318,319,321 },
                        { 319,318,327,325 },
                        { 319,318,327,328 },
                        { 318,319,321,322 },
                        { 318,319,321,323 },
                        { 320,319,321,322 },
                        { 320,319,321,323 },
                        { 319,321,323,324 },
                        { 319,321,323,325 } };



   for( int ii = 0; ii < numberOfTorsions; ii++ ){
      StringVector lineTokens;
      char* isNotEof = readLine( filePtr, lineTokens, lineCount, log );
      if( lineTokens.size() > 10 ){
         int index       = atoi( lineTokens[0].c_str() );
         int particle1   = atoi( lineTokens[1].c_str() );
         int particle2   = atoi( lineTokens[2].c_str() );
         int particle3   = atoi( lineTokens[3].c_str() );
         int particle4   = atoi( lineTokens[4].c_str() );
         double c0       = atof( lineTokens[5].c_str() );
         double c1       = atof( lineTokens[6].c_str() );
         double c2       = atof( lineTokens[7].c_str() );
         double c3       = atof( lineTokens[8].c_str() );
         double c4       = atof( lineTokens[9].c_str() );
         double c5       = atof( lineTokens[10].c_str() );


int bond[4] = { particle1, particle2, particle3, particle4 };
if( nextIndex >= 12 )nextIndex = 0;
int isBond = checkBondIndices( 4, targets[nextIndex], bond );
if( isBond ){
if( log )
(void) fprintf( log, "TGT %d %d [%d %d %d %d]\n", nextIndex, parity, targets[nextIndex][0], targets[nextIndex][1], targets[nextIndex][2], targets[nextIndex][3] );
         bondForce->addTorsion( particle1, particle2, particle3, particle4, c0, c1, c2, c3, c4, c5 );
} 

#else


         bondForce->addTorsion( particle1, particle2, particle3, particle4, c0, c1, c2, c3, c4, c5 );
#endif

      } else {
         char buffer[1024];
         (void) sprintf( buffer, "%s RBTorsionForce tokens incomplete at line=%d\n", methodName.c_str(), *lineCount );
         throwException(__FILE__, __LINE__, buffer );
         exit(-1);
      }
   }

#if 0
if( parity ){
   nextIndex++;
   parity = 0;
} else {
   parity = 1;
}
#endif

   // diagnostics

   if( log ){
      static const unsigned int maxPrint   = MAX_PRINT;
      unsigned int arraySize               = static_cast<unsigned int>(bondForce->getNumTorsions());
      (void) fprintf( log, "%s: sample of PeriodicTorsionForce parameters\n", methodName.c_str() );
      for( unsigned int ii = 0; ii < arraySize && ii < maxPrint; ii++ ){
         int particle1, particle2, particle3, particle4;
         double c0, c1, c2, c3, c4, c5;
         bondForce->getTorsionParameters( ii, particle1, particle2, particle3, particle4, c0, c1, c2, c3, c4, c5 );
         (void) fprintf( log, "%8d %8d %8d %8d %8d %14.7e %14.7e %14.7e %14.7e %14.7e %14.7e\n",
                         ii, particle1, particle2, particle3, particle4, c0, c1, c2, c3, c4, c5 );
      }
      if(  arraySize > maxPrint ){
         for( unsigned int ii = arraySize - maxPrint; ii < arraySize; ii++ ){
            int particle1, particle2, particle3, particle4;
            double c0, c1, c2, c3, c4, c5;
            bondForce->getTorsionParameters( ii, particle1, particle2, particle3, particle4, c0, c1, c2, c3, c4, c5 );
            (void) fprintf( log, "%8d %8d %8d %8d %8d %14.7e %14.7e %14.7e %14.7e %14.7e %14.7e\n",
                            ii, particle1, particle2, particle3, particle4, c0, c1, c2, c3, c4, c5 );
         }
      }
   }

   return bondForce->getNumTorsions();
}

/**---------------------------------------------------------------------------------------

   Read NonbondedExceptions parameters

   @param filePtr                     file pointer to parameter file
   @param includeNonbondedExceptions  if set, then include exceptions; otherwise set charge and epsilon to zero and
                                      sigma to 1
   @param tokens                      array of strings from first line of parameter file for this block of parameters
   @param nonbondedForce              NonBondedForce reference
   @param lineCount                   used to track line entries read from parameter file
   @param log                         log file pointer -- may be NULL

   @return number of parameters read

   --------------------------------------------------------------------------------------- */

static int readNonbondedExceptions( FILE* filePtr, int includeNonbondedExceptions, const StringVector& tokens, NonbondedForce& nonbondedForce, int* lineCount, FILE* log ){

// ---------------------------------------------------------------------------------------

   static const std::string methodName      = "readNonbondedExceptions";
   
// ---------------------------------------------------------------------------------------

   if( tokens.size() < 1 ){
      char buffer[1024];
      (void) sprintf( buffer, "%s no Nonbonded bonds number entry???\n", methodName.c_str() );
      throwException(__FILE__, __LINE__, buffer );
      exit(-1);
   }

   int numberOfExceptions           = atoi( tokens[1].c_str() );
   if( log ){
      (void) fprintf( log, "%s number of NonbondedExceptions terms=%d\n", methodName.c_str(), numberOfExceptions );
      (void) fflush( log );
   }
   for( int ii = 0; ii < numberOfExceptions; ii++ ){
      StringVector lineTokens;
      char* isNotEof = readLine( filePtr, lineTokens, lineCount, log );
      if( lineTokens.size() > 5 ){
         int index                   = atoi( lineTokens[0].c_str() );
         int particle1               = atoi( lineTokens[1].c_str() );
         int particle2               = atoi( lineTokens[2].c_str() );
         double charge               = includeNonbondedExceptions ? atof( lineTokens[3].c_str() ) : 0.0;
         double sigma                = includeNonbondedExceptions ? atof( lineTokens[4].c_str() ) : 1.0;
         double epsilon              = includeNonbondedExceptions ? atof( lineTokens[5].c_str() ) : 0.0;
         double softcoreLJLambda     = 1.0;
         if( lineTokens.size() > 4 ){ 
            softcoreLJLambda = includeNonbondedExceptions ? atof( lineTokens[4].c_str() ) : 1.0; 
         }    

#if 0
   if( log && ii < 2 )
      (void) fprintf( log, "************************ Setting q to zero ************************\n" );
   charge  = 0.0;
//   sigma   = 1.0;
//   epsilon = 0.0;
#endif
#ifdef FREE_ENERGY_BRANCH
         nonbondedForce.addException( particle1, particle2, charge, sigma, epsilon, softcoreLJLambda );
#else
         nonbondedForce.addException( particle1, particle2, charge, sigma, epsilon );
#endif
      } else if( log ){
         char buffer[1024];
         (void) sprintf( buffer, "%s readNonbondedExceptions tokens incomplete at line=%d\n", methodName.c_str(), *lineCount );
         throwException(__FILE__, __LINE__, buffer );
         exit(-1);
      }
   }

   // diagnostics

   if( log ){
      static const unsigned int maxPrint   = MAX_PRINT;
      unsigned int arraySize               = static_cast<unsigned int>(nonbondedForce.getNumExceptions());
      (void) fprintf( log, "%s: sample of NonbondedExceptions parameters\n", methodName.c_str() );
      for( unsigned int ii = 0; ii < arraySize && ii < maxPrint; ii++ ){
         int particle1, particle2;
#ifdef FREE_ENERGY_BRANCH
         double chargeProd, sigma, epsilon, softcoreLJLambda;
         nonbondedForce.getExceptionParameters( ii, particle1, particle2, chargeProd, sigma, epsilon, softcoreLJLambda );
         (void) fprintf( log, "%8d %8d %8d %14.7e %14.7e %14.7e %14.7e\n", ii, particle1, particle2, chargeProd, sigma, epsilon, softcoreLJLambda );
#else
         double chargeProd, sigma, epsilon;
         nonbondedForce.getExceptionParameters( ii, particle1, particle2, chargeProd, sigma, epsilon );
         (void) fprintf( log, "%8d %8d %8d %14.7e %14.7e %14.7e\n", ii, particle1, particle2, chargeProd, sigma, epsilon );
#endif
      }
      if( arraySize > maxPrint ){
         for( unsigned int ii = arraySize - maxPrint; ii < arraySize; ii++ ){
            int particle1, particle2;
#ifdef FREE_ENERGY_BRANCH
            double chargeProd, sigma, epsilon, softcoreLJLambda;
            nonbondedForce.getExceptionParameters( ii, particle1, particle2, chargeProd, sigma, epsilon, softcoreLJLambda );
            (void) fprintf( log, "%8d %8d %8d %14.7e %14.7e %14.7e %14.7e\n", ii, particle1, particle2, chargeProd, sigma, epsilon, softcoreLJLambda );
#else
            double chargeProd, sigma, epsilon;
            nonbondedForce.getExceptionParameters( ii, particle1, particle2, chargeProd, sigma, epsilon );
            (void) fprintf( log, "%8d %8d %8d %14.7e %14.7e %14.7e\n", ii, particle1, particle2, chargeProd, sigma, epsilon );
#endif
         }
      }
   }

   return nonbondedForce.getNumExceptions();
}

/**---------------------------------------------------------------------------------------

   Read NonbondedForce parameters

   @param filePtr              file pointer to parameter file
   @param tokens               array of strings from first line of parameter file for this block of parameters
   @param system               System reference
   @param lineCount            used to track line entries read from parameter file
   @param log                  log file pointer -- may be NULL

   @return number of parameters read

   --------------------------------------------------------------------------------------- */

static int readNonbondedForce( FILE* filePtr, int forceFlag, const StringVector& tokens,
                               System& system, int* lineCount, FILE* log ){

// ---------------------------------------------------------------------------------------

   static const std::string methodName      = "readNonbondedForce";
   
// ---------------------------------------------------------------------------------------

   if( tokens.size() < 1 ){
      char buffer[1024];
      (void) sprintf( buffer, "%s no Nonbonded bonds number entry???\n", methodName.c_str() );
      throwException(__FILE__, __LINE__, buffer );
      exit(-1);
   }

   NonbondedForce* nonbondedForce  = new NonbondedForce();

#ifdef FREE_ENERGY_BRANCH
//double lambda = 1.0;
double lambda = 0.05;
nonbondedForce->setSoftCoreLJLambda( lambda );
if( log ){
   (void) fprintf( log, "************************ Setting softCoreLJLambda=%.5f ************************\n", lambda );
}
#endif

   if( forceFlag == 0 || (forceFlag & NB_FORCE) || (forceFlag & NB_EXCEPTION_FORCE) ){
      system.addForce( nonbondedForce );
      if( log && forceFlag ){
         if( forceFlag & NB_FORCE ){
            (void) fprintf( log, "nonbonded force is being included.\n" );
         } else {
            (void) fprintf( log, "nonbonded force is not being included.\n" );
         }
         if( forceFlag & NB_EXCEPTION_FORCE ){
            (void) fprintf( log, "nonbonded exceptions are being included.\n" );
         } else {
            (void) fprintf( log, "nonbonded exceptions are not being included.\n" );
         }
      }
   } else if( forceFlag && log ){
      (void) fprintf( log, "nonbonded force is not being included.\n" );
   }

   int numberOfParticles           = atoi( tokens[1].c_str() );
   if( log ){
      (void) fprintf( log, "%s number of NonbondedForce terms=%d\n", methodName.c_str(), numberOfParticles );
      (void) fflush( log );
   }

   // get charge, sigma, epsilon for each particle

   int includeNonbonded = 1;
   if( forceFlag ){
      includeNonbonded = forceFlag & NB_FORCE;
   } 

   int includeNonbondedExceptions = 1;
   if( forceFlag ){
      includeNonbondedExceptions = forceFlag & NB_EXCEPTION_FORCE;
   } 

   for( int ii = 0; ii < numberOfParticles; ii++ ){
      StringVector lineTokens;
      char* isNotEof = readLine( filePtr, lineTokens, lineCount, log );
      if( lineTokens.size() > 3 ){
         int index                   = atoi( lineTokens[0].c_str() );
         double charge               = includeNonbonded ? atof( lineTokens[1].c_str() ) : 0.0;
         double sigma                = includeNonbonded ? atof( lineTokens[2].c_str() ) : 1.0;
         double epsilon              = includeNonbonded ? atof( lineTokens[3].c_str() ) : 0.0;
#ifdef FREE_ENERGY_BRANCH
         double softcoreLJLambda     = 1.0;
         if( lineTokens.size() > 4 ){ 
            softcoreLJLambda = includeNonbondedExceptions ? atof( lineTokens[4].c_str() ) : 1.0; 
         } else {    
            softcoreLJLambda     = ((ii % 3) == 0) ? lambda : 1.0;
         }
         nonbondedForce->addParticle( charge, sigma, epsilon, softcoreLJLambda );
#else
         nonbondedForce->addParticle( charge, sigma, epsilon );
#endif
      } else {
         char buffer[1024];
         (void) sprintf( buffer, "%s NonbondedForce tokens incomplete at line=%d\n", methodName.c_str(), *lineCount );
         throwException(__FILE__, __LINE__, buffer );
         exit(-1);
      }
   }

   // get cutoff distance, exceptions, periodic box, method,

   char* isNotEof                 = "1";
   int hits                       = 0;
   while( hits < 6 ){
      StringVector tokens;
      isNotEof = readLine( filePtr, tokens, lineCount, log );
      if( isNotEof && tokens.size() > 0 ){

         std::string field       = tokens[0];
         if( field.compare( "#" ) == 0 ){

            // skip
            if( log ){
                  (void) fprintf( log, "skip <%s>\n", field.c_str());
            }

         } else if( field.compare( "CutoffDistance" ) == 0 ){
            nonbondedForce->setCutoffDistance( atof( tokens[1].c_str() ) );
            hits++;
         } else if( field.compare( "RFDielectric" ) == 0 ){
            nonbondedForce->setReactionFieldDielectric( atof( tokens[1].c_str() ) );
            hits++;
         } else if( field.compare( "EwaldRTolerance" ) == 0 ){
            nonbondedForce->setEwaldErrorTolerance( atof( tokens[1].c_str() ) );
            hits++;
         } else if( field.compare( "NonbondedForceExceptions" ) == 0 ){
            readNonbondedExceptions( filePtr, includeNonbondedExceptions, tokens, *nonbondedForce, lineCount, log );
            hits++;
         } else if( field.compare( "Box" ) == 0 ){

            hits++;
            std::vector< Vec3 > box;
            box.resize( 3 );
            int xyzIndex = 0;
            int boxIndex = 0;
            for( int ii = 1; ii < 10; ii++ ){
               box[boxIndex][xyzIndex++] = atof( tokens[ii].c_str() );
               if( xyzIndex == 3 ){
                  xyzIndex = 0;
                  boxIndex++;
               }
            }
            nonbondedForce->setPeriodicBoxVectors( box[0], box[1], box[2] );

         } else if( field.compare( "NonbondedForceMethod" ) == 0 ){
            std::string nonbondedForceMethod = tokens[1];
            hits++;
            if( nonbondedForceMethod.compare( "NoCutoff" ) == 0 ){
               nonbondedForce->setNonbondedMethod( NonbondedForce::NoCutoff ); 
            } else if( nonbondedForceMethod.compare( "CutoffNonPeriodic" ) == 0 ){
               nonbondedForce->setNonbondedMethod( NonbondedForce::CutoffNonPeriodic ); 
            } else if( nonbondedForceMethod.compare( "CutoffPeriodic" ) == 0 ){
               nonbondedForce->setNonbondedMethod( NonbondedForce::CutoffPeriodic ); 
            } else if( nonbondedForceMethod.compare( "Ewald" ) == 0 ){
               nonbondedForce->setNonbondedMethod( NonbondedForce::Ewald ); 
            } else if( nonbondedForceMethod.compare( "PME" ) == 0 ){
               nonbondedForce->setNonbondedMethod( NonbondedForce::PME ); 
            } else {
               char buffer[1024];
               (void) sprintf( buffer, "nonbondedForce NonbondedForceMethod <%s> is not recognized at line=%d\n", nonbondedForceMethod.c_str(), lineCount );
               if( log ){
                  (void) fprintf( log, "%s", buffer ); (void) fflush( log );
               }
               throwException(__FILE__, __LINE__, buffer );
               exit(-1);
            }
         }
      }
   }

   // diagnostics

   if( log ){
      static const unsigned int maxPrint   = MAX_PRINT;
      unsigned int arraySize               = static_cast<unsigned int>(nonbondedForce->getNumParticles());

      (void) fprintf( log, "%s: nonbonded parameters\n", methodName.c_str() );

      // cutoff distance and box

      (void) fprintf( log, "CutoffDistance %14.7e\n", nonbondedForce->getCutoffDistance() );
      Vec3 a, b, c;
      nonbondedForce->getPeriodicBoxVectors( a, b, c);
      (void) fprintf( log, "Box [%14.7f %14.7f %14.7f]\n    [%14.7f %14.7f %14.7f]\n    [%14.7f %14.7f %14.7f]\n",
                      a[0], a[1], a[2], b[0], b[1], b[2], c[0], c[1], c[2] );
  
      // nonbond method

      std::string nonbondedForceMethod;
      switch( nonbondedForce->getNonbondedMethod() ){
          case NonbondedForce::NoCutoff:
              nonbondedForceMethod = "NoCutoff";
              break;
          case NonbondedForce::CutoffNonPeriodic:
              nonbondedForceMethod = "CutoffNonPeriodic";
              break;
          case NonbondedForce::CutoffPeriodic:
              nonbondedForceMethod = "CutoffPeriodic";
              break;
          case NonbondedForce::Ewald:
              nonbondedForceMethod = "Ewald";
              break;
          default:
              nonbondedForceMethod = "Unknown";
      }
      (void) fprintf( log, "NonbondedForceMethod=%s\n", nonbondedForceMethod.c_str() );
  
#ifdef FREE_ENERGY_BRANCH
      (void) fprintf( log, "charge, sigma, epsilon, softcoreLJLambda\n" );
#else
      (void) fprintf( log, "charge, sigma, epsilon\n" );
#endif

      for( unsigned int ii = 0; ii < arraySize && ii < maxPrint; ii++ ){
#ifdef FREE_ENERGY_BRANCH
         double charge, sigma, epsilon, softcoreLJLambda;
         nonbondedForce->getParticleParameters( ii, charge, sigma, epsilon, softcoreLJLambda );
         (void) fprintf( log, "%8d %14.7e %14.7e %14.7e %14.7e\n", ii, charge, sigma, epsilon, softcoreLJLambda );
#else
         double charge, sigma, epsilon;
         nonbondedForce->getParticleParameters( ii, charge, sigma, epsilon );
         (void) fprintf( log, "%8d %14.7e %14.7e %14.7e\n", ii, charge, sigma, epsilon );
#endif
      }
      if( arraySize > maxPrint ){
         for( unsigned int ii = arraySize - maxPrint; ii < arraySize; ii++ ){
#ifdef FREE_ENERGY_BRANCH
            double charge, sigma, epsilon, softcoreLJLambda;
            nonbondedForce->getParticleParameters( ii, charge, sigma, epsilon, softcoreLJLambda );
            (void) fprintf( log, "%8d %14.7e %14.7e %14.7e %14.7e\n", ii, charge, sigma, epsilon, softcoreLJLambda );
#else
            double charge, sigma, epsilon;
            nonbondedForce->getParticleParameters( ii, charge, sigma, epsilon );
            (void) fprintf( log, "%8d %14.7e %14.7e %14.7e\n", ii, charge, sigma, epsilon );
#endif
         }
      }
   }

   return nonbondedForce->getNumParticles();
}

/**---------------------------------------------------------------------------------------

   Read GBSAOBCForce parameters

   @param filePtr              file pointer to parameter file
   @param forceFlag            flag signalling whether force is to be added to system
                               if force == 0 || forceFlag & GBSA_OBC_FORCE, then included
   @param tokens               array of strings from first line of parameter file for this block of parameters
   @param system               System reference
   @param lineCount            used to track line entries read from parameter file
   @param log                  log file pointer -- may be NULL

   @return number of parameters read

   --------------------------------------------------------------------------------------- */

static int readGBSAOBCForce( FILE* filePtr, int forceFlag, const StringVector& tokens, System& system, int* lineCount, FILE* log ){

// ---------------------------------------------------------------------------------------

   static const std::string methodName      = "readGBSAOBCForce";
   
// ---------------------------------------------------------------------------------------

   if( tokens.size() < 1 ){
      char buffer[1024];
      (void) sprintf( buffer, "%s no GBSAOBC terms entry???\n", methodName.c_str() );
      throwException(__FILE__, __LINE__, buffer );
      exit(-1);
   }

   GBSAOBCForce* gbsaObcForce = new GBSAOBCForce();
   if( forceFlag == 0 || forceFlag & GBSA_OBC_FORCE ){
      system.addForce( gbsaObcForce );
      if( log && forceFlag ){
         (void) fprintf( log, "GBSA OBC force is being included.\n" );
      }
   } else if( forceFlag && log ){
      (void) fprintf( log, "GBSA OBC force is not being included.\n" );
   }

   int numberOfParticles           = atoi( tokens[1].c_str() );
   if( log ){
      (void) fprintf( log, "%s number of GBSAOBCForce terms=%d\n", methodName.c_str(), numberOfParticles );
      (void) fflush( log );
   }
   for( int ii = 0; ii < numberOfParticles; ii++ ){
      StringVector lineTokens;
      char* isNotEof = readLine( filePtr, lineTokens, lineCount, log );
      if( lineTokens.size() > 3 ){
         int index            = atoi( lineTokens[0].c_str() );
         double charge        = atof( lineTokens[1].c_str() );
         double radius        = atof( lineTokens[2].c_str() );
         double scalingFactor = atof( lineTokens[3].c_str() );
#ifdef FREE_ENERGY_BRANCH
double saScale;
if( (ii % 3 ) == 0 ){
   saScale = 0.0;
} else {
   saScale = 1.0;
}
   
         gbsaObcForce->addParticle( charge, radius, scalingFactor, saScale );
#else
         gbsaObcForce->addParticle( charge, radius, scalingFactor );
#endif
      } else {
         char buffer[1024];
         (void) sprintf( buffer, "%s GBSAOBCForce tokens incomplete at line=%d\n", methodName.c_str(), *lineCount );
         throwException(__FILE__, __LINE__, buffer );
         exit(-1);
      }
   }

   char* isNotEof                 = "1";
   int hits                       = 0;
   while( hits < 2 ){
      StringVector tokens;
      isNotEof = readLine( filePtr, tokens, lineCount, log );
      if( isNotEof && tokens.size() > 0 ){

         std::string field       = tokens[0];
         if( field.compare( "SoluteDielectric" ) == 0 ){
            gbsaObcForce->setSoluteDielectric( atof( tokens[1].c_str() ) );
            hits++;
         } else if( field.compare( "SolventDielectric" ) == 0 ){
            gbsaObcForce->setSolventDielectric( atof( tokens[1].c_str() ) );
            hits++;
         } else {
               char buffer[1024];
               (void) sprintf( buffer, "%s read past GBSA Obc block at line=%d\n", methodName.c_str(), lineCount );
               throwException(__FILE__, __LINE__, buffer );
               exit(-1);
         }
      } else {
         char buffer[1024];
         (void) sprintf( buffer, "%s invalid token count at line=%d?\n", methodName.c_str(), lineCount );
         throwException(__FILE__, __LINE__, buffer );
         exit(-1);
      }
   }

   if( log ){
      static const unsigned int maxPrint   = MAX_PRINT;
      unsigned int arraySize               = static_cast<unsigned int>(gbsaObcForce->getNumParticles());
      (void) fprintf( log, "%s: sample of GBSA OBC Force parameters; no. of particles=%d\n",
                      methodName.c_str(), gbsaObcForce->getNumParticles() );
      (void) fprintf( log, "solute/solvent dielectrics: [%10.4f %10.4f]\n",
                      gbsaObcForce->getSoluteDielectric(),  gbsaObcForce->getSolventDielectric() );

      for( unsigned int ii = 0; ii < arraySize && ii < maxPrint; ii++ ){
#ifdef FREE_ENERGY_BRANCH
         double charge, radius, scalingFactor, nonpolarScaleFactor;
         gbsaObcForce->getParticleParameters( ii, charge, radius, scalingFactor, nonpolarScaleFactor );
         (void) fprintf( log, "%8d  %14.7e %14.7e %14.7e %14.7e\n", ii, charge, radius, scalingFactor, nonpolarScaleFactor );
#else
         double charge, radius, scalingFactor;
         gbsaObcForce->getParticleParameters( ii, charge, radius, scalingFactor );
         (void) fprintf( log, "%8d  %14.7e %14.7e %14.7e\n", ii, charge, radius, scalingFactor );
#endif
      }
      if( arraySize > maxPrint ){
         for( unsigned int ii = arraySize - maxPrint; ii < arraySize; ii++ ){
#ifdef FREE_ENERGY_BRANCH
            double charge, radius, scalingFactor, nonpolarScaleFactor;
            gbsaObcForce->getParticleParameters( ii, charge, radius, scalingFactor, nonpolarScaleFactor );
            (void) fprintf( log, "%8d  %14.7e %14.7e %14.7e %14.7e\n", ii, charge, radius, scalingFactor, nonpolarScaleFactor );
#else
            double charge, radius, scalingFactor;
            gbsaObcForce->getParticleParameters( ii, charge, radius, scalingFactor );
            (void) fprintf( log, "%8d  %14.7e %14.7e %14.7e\n", ii, charge, radius, scalingFactor );
#endif
         }
      }
   }

   return gbsaObcForce->getNumParticles();
}

/**---------------------------------------------------------------------------------------

   Read Constraints

   @param filePtr              file pointer to parameter file
   @param tokens               array of strings from first line of parameter file for this block of parameters
   @param system               System reference
   @param lineCount            used to track line entries read from parameter file
   @param log                  log file pointer -- may be NULL

   @return number of parameters read

   --------------------------------------------------------------------------------------- */

static int readConstraints( FILE* filePtr, const StringVector& tokens, System& system, int* lineCount, FILE* log ){

// ---------------------------------------------------------------------------------------

   static const std::string methodName      = "readConstraints";
   
// ---------------------------------------------------------------------------------------

   if( tokens.size() < 1 ){
      char buffer[1024];
      (void) sprintf( buffer, "%s no Constraints terms entry???\n", methodName.c_str() );
      throwException(__FILE__, __LINE__, buffer );
      exit(-1);
   }

   int numberOfConstraints = atoi( tokens[1].c_str() );
   if( log ){
      (void) fprintf( log, "%s number of constraints=%d\n", methodName.c_str(), numberOfConstraints );
      (void) fflush( log );
   }
   for( int ii = 0; ii < numberOfConstraints; ii++ ){
      StringVector lineTokens;
      char* isNotEof = readLine( filePtr, lineTokens, lineCount, log );
      if( lineTokens.size() > 3 ){
         int index            = atoi( lineTokens[0].c_str() );
         int particle1        = atoi( lineTokens[1].c_str() );
         int particle2        = atoi( lineTokens[2].c_str() );
         double distance      = atof( lineTokens[3].c_str() );
         system.addConstraint( particle1, particle2, distance );
      } else {
         char buffer[1024];
         (void) sprintf( buffer, "%s constraint tokens incomplete at line=%d\n", methodName.c_str(), *lineCount );
         throwException(__FILE__, __LINE__, buffer );
         exit(-1);
      }
   }

   if( log ){
      static const unsigned int maxPrint   = MAX_PRINT;
      unsigned int arraySize               = static_cast<unsigned int>(system.getNumConstraints());
      (void) fprintf( log, "%s: sample constraints\n", methodName.c_str() );
      for( unsigned int ii = 0; ii < arraySize && ii < maxPrint; ii++ ){
         int particle1, particle2;
         double distance;
         system.getConstraintParameters( ii, particle1, particle2, distance ); 
         (void) fprintf( log, "%8d %8d %8d %14.7e\n", ii, particle1, particle2, distance );
      }
      if( arraySize > maxPrint ){
         for( unsigned int ii = arraySize - maxPrint; ii < arraySize; ii++ ){
            int particle1, particle2;
            double distance;
            system.getConstraintParameters( ii, particle1, particle2, distance ); 
            (void) fprintf( log, "%8d %8d %8d %14.7e\n", ii, particle1, particle2, distance );
         }
      }
   }

   return system.getNumConstraints();
}

/**---------------------------------------------------------------------------------------

   Read integrator

   @param filePtr              file pointer to parameter file
   @param tokens               array of strings from first line of parameter file for this block of parameters
   @param system               System reference
   @param lineCount            used to track line entries read from parameter file
   @param log                  log file pointer -- may be NULL

   @return integrator

   --------------------------------------------------------------------------------------- */

static Integrator* readIntegrator( FILE* filePtr, const StringVector& tokens, System& system, int* lineCount, FILE* log ){

// ---------------------------------------------------------------------------------------

   static const std::string methodName      = "readIntegrator";
   
// ---------------------------------------------------------------------------------------

   if( tokens.size() < 1 ){
      char buffer[1024];
      (void) sprintf( buffer, "%s integrator name missing?\n", methodName.c_str() );
      throwException(__FILE__, __LINE__, buffer );
      exit(-1);
   }

   std::string integratorName = tokens[1];
   if( log ){
      (void) fprintf( log, "%s integrator=%s\n", methodName.c_str(), integratorName.c_str() );
      (void) fflush( log );
   }

   // set number of parameters (lines to read)

   int readLines;
   if( integratorName.compare( "LangevinIntegrator" ) == 0 ){
      readLines = 5;
   } else if( integratorName.compare( "VariableLangevinIntegrator" ) == 0 ){
      readLines = 6;
   } else if( integratorName.compare( "VerletIntegrator" ) == 0 ){
      readLines = 2;
   } else if( integratorName.compare( "VariableVerletIntegrator" ) == 0 ){
      readLines = 3;
   } else if( integratorName.compare( "BrownianIntegrator" ) == 0 ){
      readLines = 5;
   } else {
      (void) fprintf( log, "%s integrator=%s not recognized.\n", methodName.c_str(), integratorName.c_str() );
      (void) fflush( log );
      exit(-1);
   }
   
   // read in parameters

   double stepSize               = 0.001;
   double constraintTolerance    = 1.0e-05;
   double temperature            = 300.0;
   double friction               = 0.01099;
   double errorTolerance         = 1.0e-05;
   int randomNumberSeed          = 1993;

   for( int ii = 0; ii < readLines; ii++ ){
      StringVector lineTokens;
      char* isNotEof = readLine( filePtr, lineTokens, lineCount, log );
      if( lineTokens.size() > 1 ){
         if( lineTokens[0].compare( "StepSize" ) == 0 ){
            stepSize            =  atof( lineTokens[1].c_str() );
         } else if( lineTokens[0].compare( "ConstraintTolerance" ) == 0 ){
            constraintTolerance =  atof( lineTokens[1].c_str() );
         } else if( lineTokens[0].compare( "Temperature" ) == 0 ){
            temperature         =  atof( lineTokens[1].c_str() );
         } else if( lineTokens[0].compare( "Friction" ) == 0 ){
            friction            =  atof( lineTokens[1].c_str() );
         } else if( lineTokens[0].compare( "ErrorTolerance" ) == 0 ){
            errorTolerance      =  atof( lineTokens[1].c_str() );
         } else if( lineTokens[0].compare( "RandomNumberSeed" ) == 0 ){
            randomNumberSeed    =  atoi( lineTokens[1].c_str() );
         } else {
            (void) fprintf( log, "%s integrator field=%s not recognized.\n", methodName.c_str(), lineTokens[0].c_str() );
            (void) fflush( log );
            exit(-1);
         }
      } else {
         char buffer[1024];
         (void) sprintf( buffer, "%s integrator parameters incomplete at line=%d\n", methodName.c_str(), *lineCount );
         throwException(__FILE__, __LINE__, buffer );
         exit(-1);
      }
   }

   // build integrator

   Integrator* returnIntegrator = NULL;

   if( integratorName.compare( "LangevinIntegrator" ) == 0 ){
      returnIntegrator = new LangevinIntegrator( temperature, friction, stepSize );
//      returnIntegrator->setRandomNumberSeed( randomNumberSeed );
   } else if( integratorName.compare( "VariableLangevinIntegrator" ) == 0 ){
      returnIntegrator = new VariableLangevinIntegrator( temperature, friction, errorTolerance );
      returnIntegrator->setStepSize( stepSize );
//      returnIntegrator->setRandomNumberSeed( randomNumberSeed );
   } else if( integratorName.compare( "VerletIntegrator" ) == 0 ){
      returnIntegrator = new VerletIntegrator( stepSize );
   } else if( integratorName.compare( "VariableVerletIntegrator" ) == 0 ){
      returnIntegrator = new VariableVerletIntegrator( errorTolerance );
      returnIntegrator->setStepSize( stepSize );
   } else if( integratorName.compare( "BrownianIntegrator" ) == 0 ){
      returnIntegrator = new BrownianIntegrator( temperature, friction, stepSize );
//      returnIntegrator->setRandomNumberSeed( randomNumberSeed );
   }
   returnIntegrator->setConstraintTolerance( constraintTolerance );
   
   if( log ){
      static const unsigned int maxPrint   = MAX_PRINT;
      (void) fprintf( log, "%s: parameters\n", methodName.c_str() );
      (void) fprintf( log, "StepSize=%14.7e constraint tolerance=%14.7e ", stepSize, constraintTolerance );
      if( integratorName.compare( "LangevinIntegrator" ) == 0 || 
          integratorName.compare( "BrownianIntegrator" ) == 0 ||  
          integratorName.compare( "VariableLangevinIntegrator" ) == 0 ){
          (void) fprintf( log, "Temperature=%14.7e friction=%14.7e seed=%d (seed may not be set!) ", temperature, friction, randomNumberSeed );
      }
      if( integratorName.compare( "VariableLangevinIntegrator" ) == 0 || 
          integratorName.compare( "VariableVerletIntegrator" ) == 0 ){
          (void) fprintf( log, "Error tolerance=%14.7e", errorTolerance);
      }
      (void) fprintf( log, "\n" );
   }

   return returnIntegrator;
}

/**---------------------------------------------------------------------------------------

   Read arrays of Vec3s (coordinates/velocities/forces/...)

   @param filePtr              file pointer to parameter file
   @param tokens               array of strings from first line of parameter file for this block of parameters
   @param coordinates          Vec3 array
   @param lineCount            used to track line entries read from parameter file
   @param log                  log file pointer -- may be NULL

   @return number of entries read

   --------------------------------------------------------------------------------------- */

static int readVec3( FILE* filePtr, const StringVector& tokens, std::vector<Vec3>& coordinates, int* lineCount, FILE* log ){

// ---------------------------------------------------------------------------------------

   static const std::string methodName      = "readVec3";
   
// ---------------------------------------------------------------------------------------

   if( tokens.size() < 1 ){
      char buffer[1024];
      (void) sprintf( buffer, "%s no Coordinates terms entry???\n", methodName.c_str() );
      throwException(__FILE__, __LINE__, buffer );
      exit(-1);
   }

   int numberOfCoordinates= atoi( tokens[1].c_str() );
   if( log ){
      (void) fprintf( log, "%s number of coordinates=%d\n", methodName.c_str(), numberOfCoordinates );
      (void) fflush( log );
   }
   for( int ii = 0; ii < numberOfCoordinates; ii++ ){
      StringVector lineTokens;
      char* isNotEof = readLine( filePtr, lineTokens, lineCount, log );
      if( lineTokens.size() > 3 ){
         int index            = atoi( lineTokens[0].c_str() );
         double xCoord        = atof( lineTokens[1].c_str() );
         double yCoord        = atof( lineTokens[2].c_str() );
         double zCoord        = atof( lineTokens[3].c_str() );
         coordinates.push_back( Vec3( xCoord, yCoord, zCoord ) );
      } else {
         char buffer[1024];
         (void) sprintf( buffer, "%s coordinates tokens incomplete at line=%d\n", methodName.c_str(), *lineCount );
         throwException(__FILE__, __LINE__, buffer );
         exit(-1);
      }
   }

   // diagnostics

   if( log ){
      static const unsigned int maxPrint = MAX_PRINT;
      (void) fprintf( log, "%s: sample of vec3: %u\n", methodName.c_str(), coordinates.size() );
      for( unsigned int ii = 0; ii < coordinates.size() && ii < maxPrint; ii++ ){
         (void) fprintf( log, "%6u [%14.7e %14.7e %14.7e]\n", ii,
                         coordinates[ii][0], coordinates[ii][1], coordinates[ii][2] );
      }
      if( coordinates.size() > maxPrint ){
         for( unsigned int ii = coordinates.size() - maxPrint; ii < coordinates.size(); ii++ ){
            (void) fprintf( log, "%6u [%14.7e %14.7e %14.7e]\n", ii,
                            coordinates[ii][0], coordinates[ii][1], coordinates[ii][2] );
         }
      }
   }

   return static_cast<int>(coordinates.size());
}

/**---------------------------------------------------------------------------------------

   Read parameter file

   @param inputParameterFile   input parameter file name
   @param system               system to which forces based on parameters are to be added
   @param coordinates          Vec3 array containing coordinates on output
   @param velocities           Vec3 array containing velocities on output
   @param inputLog             log file pointer -- may be NULL

   @return number of lines read

   --------------------------------------------------------------------------------------- */

Integrator* readParameterFile( const std::string& inputParameterFile, int forceFlag, System& system,
                               std::vector<Vec3>& coordinates, 
                               std::vector<Vec3>& velocities,
                               std::vector<Vec3>& forces, double* kineticEnergy, double* potentialEnergy,
                               FILE* inputLog ){

// ---------------------------------------------------------------------------------------

   static const std::string methodName      = "readParameterFile";
   int PrintOn                              = 1; 
   
// ---------------------------------------------------------------------------------------

   FILE* log;
	if( PrintOn == 0 && inputLog ){
      log = NULL;
   } else {
      log = inputLog;
   }

   if( log ){
      (void) fprintf( log, "%s\n", methodName.c_str() );
      (void) fflush( log );
   }   

   // open parameter file

   FILE* filePtr;
#ifdef _MSC_VER
   fopen_s( &filePtr, inputParameterFile.c_str(), "r" );
#else
   filePtr = fopen( inputParameterFile.c_str(), "r" );
#endif

   if( filePtr == NULL ){
      char buffer[1024];
      (void) sprintf( buffer, "Input parameter file=<%s> could not be opened -- aborting.\n", methodName.c_str(), inputParameterFile.c_str() );
      throwException(__FILE__, __LINE__, buffer );
      (void) fflush( stderr);
      exit(-1);
   } else if( log ){
      (void) fprintf( log, "Input parameter file=<%s> opened.\n", methodName.c_str(), inputParameterFile.c_str() );
   }

   int lineCount                  = 0;
   std::string version            = "0.1"; 
   char* isNotEof                 = "1";
   Integrator* returnIntegrator   = NULL;

   // loop over lines in file

   while( isNotEof ){

      // read line and continue if not EOF and tokens found on line

      StringVector tokens;
      isNotEof = readLine( filePtr, tokens, &lineCount, log );

      if( isNotEof && tokens.size() > 0 ){

         std::string field       = tokens[0];

         if( log ){
          (void) fprintf( log, "Field=<%s> at line=%d\n", field.c_str(), lineCount );
         }
 
         if( field.compare( "Version" ) == 0 ){
            if( tokens.size() > 1 ){
               version = tokens[1];
               if( log ){
                  (void) fprintf( log, "Version=<%s> at line=%d\n", version.c_str(), lineCount );
               }
            }
         } else if( field.compare( "Particles" ) == 0 ){
            readParticles( filePtr, tokens, system, &lineCount, log );
         } else if( field.compare( "Masses" ) == 0 ){
            readMasses( filePtr, tokens, system, &lineCount, log );
         } else if( field.compare( "NumberOfForces" ) == 0 ){
            // skip
         } else if( field.compare( "CMMotionRemover" ) == 0 ){
            int frequency = atoi( tokens[1].c_str() );
            system.addForce( new CMMotionRemover( frequency ) );
            if( log ){
               (void) fprintf( log, "CMMotionRemover added w/ frequency=%d at line=%d\n", frequency, lineCount );
            }
         } else if( field.compare( "HarmonicBondForce" ) == 0 ){
            readHarmonicBondForce( filePtr, forceFlag, tokens, system, &lineCount, log );
         } else if( field.compare( "HarmonicAngleForce" ) == 0 ){
            readHarmonicAngleForce( filePtr, forceFlag, tokens, system, &lineCount, log );
         } else if( field.compare( "PeriodicTorsionForce" ) == 0 ){
            readPeriodicTorsionForce( filePtr, forceFlag, tokens, system, &lineCount, log );
         } else if( field.compare( "RBTorsionForce" ) == 0 ){
            readRBTorsionForce( filePtr, forceFlag, tokens, system, &lineCount, log );
         } else if( field.compare( "NonbondedForce" ) == 0 ){
            readNonbondedForce( filePtr, forceFlag, tokens, system, &lineCount, log );
         } else if( field.compare( "GBSAOBCForce" ) == 0 ){
            readGBSAOBCForce( filePtr, forceFlag, tokens, system, &lineCount, log );
         } else if( field.compare( "Constraints" ) == 0 ){
            readConstraints( filePtr, tokens, system, &lineCount, log );
         } else if( field.compare( "Integrator" ) == 0 ){
            returnIntegrator = readIntegrator( filePtr, tokens, system, &lineCount, log );
         } else if( field.compare( "Positions" ) == 0 ){
            readVec3( filePtr, tokens, coordinates, &lineCount, log );
         } else if( field.compare( "Velocities" ) == 0 ){
            readVec3( filePtr, tokens, velocities, &lineCount, log );
         } else if( field.compare( "Forces" ) == 0 ){
            readVec3( filePtr, tokens, forces, &lineCount, log );
         } else if( field.compare( "KineticEnergy" ) == 0 ||
                    field.compare( "PotentialEnergy" ) == 0 ){
            double value = 0.0;
            if( tokens.size() > 1 ){
               value = atof( tokens[1].c_str() );
               if( log ){
                  (void) fprintf( log, "%s =%s\n", tokens[0].c_str(), tokens[1].c_str());
               }
            } else {
               char buffer[1024];
               (void) sprintf( buffer, "Missing energy for field=<%s> at line=%d\n", field.c_str(), lineCount );
               throwException(__FILE__, __LINE__, buffer );
               exit(-1);
            }
            if( field.compare( "KineticEnergy" ) == 0 ){
               *kineticEnergy    = value;
            } else {
               *potentialEnergy  = value;
            }
         } else {
            char buffer[1024];
            (void) sprintf( buffer, "Field=<%s> not recognized at line=%d\n", field.c_str(), lineCount );
            throwException(__FILE__, __LINE__, buffer );
            exit(-1);
         }
      }
   }
 
   if( log ){
      (void) fprintf( log, "Read %d lines from file=<%s>\n", lineCount, inputParameterFile.c_str() );
   }

   return returnIntegrator;
}

/** 
 * Check that energy and force are consistent
 * 
 * @return DefaultReturnValue or ErrorReturnValue
 *
 */

#ifdef FREE_ENERGY_BRANCH
static int checkEnergyForceConsistent( OpenMMContext& context, double delta, double tolerance, FILE* log ) {
#else
static int checkEnergyForceConsistent( Context& context, double delta, double tolerance, FILE* log ) {
#endif
    
// ---------------------------------------------------------------------------------------

   static const std::string methodName    = "checkEnergyForceConsistent";

// ---------------------------------------------------------------------------------------

   int returnStatus                       = 0;

   // get positions, forces and potential energy

   int types                              = State::Positions | State::Velocities | State::Forces | State::Energy;

   State state                            = context.getState( types );

   std::vector<Vec3> coordinates          = state.getPositions();
   std::vector<Vec3> velocities           = state.getVelocities();
   std::vector<Vec3> forces               = state.getForces();
   double kineticEnergy                   = state.getKineticEnergy();
   double potentialEnergy                 = state.getPotentialEnergy();

   // conpute norm of force

   double forceNorm         = 0.0;
   for( unsigned int ii = 0; ii < forces.size(); ii++ ){
      forceNorm += forces[ii][0]*forces[ii][0] + forces[ii][1]*forces[ii][1] + forces[ii][2]*forces[ii][2];
   }
   forceNorm = std::sqrt( forceNorm );

   if( forceNorm <= 0.0 ){
      if( log ){
         (void) fprintf( log, "%s norm of force is <= 0 norm=%.3e\n", methodName.c_str(), forceNorm );
         (void) fflush( log );
      }
      return returnStatus;
   }
 
   // take step in direction of energy gradient

   double step = delta/forceNorm;
   std::vector<Vec3> perturbedPositions; 
   perturbedPositions.resize( forces.size() );
   for( unsigned int ii = 0; ii < forces.size(); ii++ ){
      perturbedPositions[ii] = Vec3( coordinates[ii][0] - step*forces[ii][0], coordinates[ii][1] - step*forces[ii][1], coordinates[ii][2] - step*forces[ii][2] ); 
   }

   context.setPositions( perturbedPositions );

   // get new potential energy

   state    = context.getState( types );

   // report energies

   double perturbedPotentialEnergy = state.getPotentialEnergy();
   double deltaEnergy              = ( perturbedPotentialEnergy - potentialEnergy )/delta;
   double difference               = fabs( deltaEnergy - forceNorm );
   double denominator              = 0.5*( fabs( deltaEnergy ) + fabs( forceNorm ) ); 
   if( denominator > 0.0 ){
      difference /= denominator;
   }

   if( log ){
      (void) fprintf( log, "%s  difference=%14.7e dE=%14.7e Pe2/1 [%14.7e %14.7e] delta=%14.7e nrm=%14.7e\n",
                      methodName.c_str(), difference, deltaEnergy, perturbedPotentialEnergy,
                      potentialEnergy, delta, forceNorm );
      (void) fflush( log );
   }

   ASSERT( difference < tolerance );

   if( log ){
      (void) fprintf( log, "\n%s passed\n", methodName.c_str() );
      (void) fflush( log );
   }
   return returnStatus;

}

/**---------------------------------------------------------------------------------------
 * Get integrator
 * 
 * @param  integratorName       integratorName (VerletIntegrator, BrownianIntegrator, LangevinIntegrator, ...)
 * @param  timeStep             time step
 * @param  friction (ps)        friction
 * @param  temperature          temperature
 * @param  shakeTolerance       Shake tolerance
 * @param  errorTolerance       Error tolerance
 * @param  randomNumberSeed     seed
 *
 * @return DefaultReturnValue or ErrorReturnValue
 *
   --------------------------------------------------------------------------------------- */

Integrator* _getIntegrator( std::string& integratorName, double timeStep,
                            double friction, double temperature,
                            double shakeTolerance, double errorTolerance,
                            int randomNumberSeed, FILE* log ){

// ---------------------------------------------------------------------------------------

   static const std::string methodName      = "_getIntegrator";

// ---------------------------------------------------------------------------------------

    // Create an integrator 
    
    Integrator* integrator;

    if( integratorName.compare( "VerletIntegrator" ) == 0 ){
        integrator = new VerletIntegrator( timeStep );
    } else if( integratorName.compare( "VariableVerletIntegrator" ) == 0 ){
        integrator = new VariableVerletIntegrator( errorTolerance );
    } else if( integratorName.compare( "BrownianIntegrator" ) == 0 ){
        integrator = new BrownianIntegrator( temperature, friction, timeStep );
    } else if( integratorName.compare( "LangevinIntegrator" ) == 0 ){
        integrator                                = new LangevinIntegrator( temperature, friction, timeStep );
        LangevinIntegrator* langevinIntegrator    = dynamic_cast<LangevinIntegrator*>(integrator);
        if( randomNumberSeed <= 0 ){
           time_t zero = time(NULL);
           langevinIntegrator->setRandomNumberSeed(static_cast<int>(zero));
        } else { 
           langevinIntegrator->setRandomNumberSeed( randomNumberSeed );
        }
    } else if( integratorName.compare( "VariableLangevinIntegrator" ) == 0 ){
        integrator                                        = new VariableLangevinIntegrator( temperature, friction, errorTolerance );
        VariableLangevinIntegrator* langevinIntegrator    = dynamic_cast<VariableLangevinIntegrator*>(integrator);
        if( randomNumberSeed <= 0 ){
           time_t zero = time(NULL);
           langevinIntegrator->setRandomNumberSeed(static_cast<int>(zero));
        } else { 
           langevinIntegrator->setRandomNumberSeed( randomNumberSeed );
        }
    } else {
       char buffer[1024];
       (void) sprintf( buffer, "%s  integrator=<%s> not recognized.\n", methodName.c_str(), integratorName.c_str() );
       if( log ){
          (void) fprintf( log , "%s", buffer );
          (void) fflush( log );
       }
       throwException(__FILE__, __LINE__, buffer );
       return NULL;
    }    

    integrator->setConstraintTolerance( shakeTolerance );
    
   return integrator;
}

/**---------------------------------------------------------------------------------------
 * Get integrator type
 * 
 * @param  integrator 
 *
 * @return name or "NotFound"
 *
   --------------------------------------------------------------------------------------- */

static std::string _getIntegratorName( Integrator* integrator ){

// ---------------------------------------------------------------------------------------

//   static const std::string methodName      = "_getIntegratorName";

// ---------------------------------------------------------------------------------------

   // LangevinIntegrator

   try {
      LangevinIntegrator& langevinIntegrator = dynamic_cast<LangevinIntegrator&>(*integrator);
      return  "LangevinIntegrator";
   } catch( std::bad_cast ){
   }

   // VariableLangevinIntegrator

   try {
      VariableLangevinIntegrator& langevinIntegrator = dynamic_cast<VariableLangevinIntegrator&>(*integrator);
      return "VariableLangevinIntegrator";
   } catch( std::bad_cast ){
   }

   // VerletIntegrator

   try {
      VerletIntegrator& verletIntegrator = dynamic_cast<VerletIntegrator&>(*integrator);
      return "VerletIntegrator";
   } catch( std::bad_cast ){
   }
    
   // VariableVerletIntegrator

   try {
      VariableVerletIntegrator & variableVerletIntegrator = dynamic_cast<VariableVerletIntegrator&>(*integrator);
      return "VariableVerletIntegrator";
   } catch( std::bad_cast ){
   }
    
   // BrownianIntegrator

   try {
      BrownianIntegrator& brownianIntegrator = dynamic_cast<BrownianIntegrator&>(*integrator);
      return "BrownianIntegrator";
   } catch( std::bad_cast ){
   }
    
   return "NotFound";
}

/**---------------------------------------------------------------------------------------
 * Set velocities based on temperature
 * 
 * @param system       System reference -- retrieve particle masses
 * @param velocities   array of Vec3 for velocities (size must be set)
 * @param temperature  temperature
 * @param log          optional log reference
 *
 * @return DefaultReturnValue
 *
   --------------------------------------------------------------------------------------- */

static int _setVelocitiesBasedOnTemperature( const System& system, std::vector<Vec3>& velocities, double temperature, FILE* log ) {
    
// ---------------------------------------------------------------------------------------

   static const std::string methodName    = "setVelocitiesBasedOnTemperature";

   double randomValues[3];

// ---------------------------------------------------------------------------------------

   // set velocities based on temperature

   temperature   *= BOLTZ;
   double randMax = static_cast<double>(RAND_MAX);
   randMax        = 1.0/randMax;
   for( unsigned int ii = 0; ii < velocities.size(); ii++ ){
      double velocityScale      = std::sqrt( temperature/system.getParticleMass(ii) );
      randomValues[0]           = randMax*( (double) rand() );
      randomValues[1]           = randMax*( (double) rand() );
      randomValues[2]           = randMax*( (double) rand() );
      velocities[ii]            = Vec3( randomValues[0]*velocityScale, randomValues[1]*velocityScale, randomValues[2]*velocityScale );
   }

   return DefaultReturnValue;
}

/**---------------------------------------------------------------------------------------
 * Print Integrator info to log
 * 
 * @param integrator   integrator
 * @param log          optional log reference
 *
 * @return DefaultReturnValue
 *
   --------------------------------------------------------------------------------------- */

static int _printIntegratorInfo( Integrator* integrator, FILE* log ){
    
// ---------------------------------------------------------------------------------------

   //static const std::string methodName    = "_printIntegratorInfo";

// ---------------------------------------------------------------------------------------

   std::string integratorName           = _getIntegratorName( integrator );
   (void) fprintf( log, "Integrator=%s stepSize=%.3f ShakeTol=%.3e\n", 
                   integratorName.c_str(), integrator->getStepSize(), integrator->getConstraintTolerance() );

   // stochastic integrators (seed, friction, temperature)

   if( integratorName.compare( "LangevinIntegrator" ) == 0 || integratorName.compare( "VariableLangevinIntegrator" ) == 0 ||
       integratorName.compare( "BrownianIntegrator" ) == 0 ){

      double temperature = 300.0;
      double friction    = 100.0;
      int seed           = 0;

      if( integratorName.compare( "LangevinIntegrator" )                == 0 ){
         LangevinIntegrator* langevinIntegrator          = dynamic_cast<LangevinIntegrator*>(integrator);
         temperature                                     = langevinIntegrator->getTemperature();
         friction                                        = langevinIntegrator->getFriction();
         seed                                            = langevinIntegrator->getRandomNumberSeed();
      } else if( integratorName.compare( "VariableLangevinIntegrator" ) == 0 ){
         VariableLangevinIntegrator* langevinIntegrator  = dynamic_cast<VariableLangevinIntegrator*>(integrator);
         temperature                                     = langevinIntegrator->getTemperature();
         friction                                        = langevinIntegrator->getFriction();
         seed                                            = langevinIntegrator->getRandomNumberSeed();
      } else if( integratorName.compare( "BrownianIntegrator" )         == 0 ){
         BrownianIntegrator* brownianIntegrator          = dynamic_cast<BrownianIntegrator*>(integrator);
         temperature                                     = brownianIntegrator->getTemperature();
         friction                                        = brownianIntegrator->getFriction();
         seed                                            = brownianIntegrator->getRandomNumberSeed();
      }
   
      (void) fprintf( log, "T=%.3f friction=%.3f seed=%d\n", temperature, friction, seed );
   }

   // variable integrators -- error tolerance

   if( integratorName.compare( "VariableLangevinIntegrator" ) == 0 || integratorName.compare( "VariableVerletIntegrator" ) == 0 ){
      double errorTolerance = 0.0;
      if( integratorName.compare( "VariableLangevinIntegrator" ) == 0 ){
         VariableLangevinIntegrator* langevinIntegrator          = dynamic_cast<VariableLangevinIntegrator*>(integrator);
         errorTolerance                                          = langevinIntegrator->getErrorTolerance();
      } else {
         VariableVerletIntegrator* verletIntegrator              = dynamic_cast<VariableVerletIntegrator*>(integrator);
         errorTolerance                                          = verletIntegrator->getErrorTolerance();
      }
      (void) fprintf( log, "Error tolerance=%.3e\n", errorTolerance );
   }

   (void) fflush( log );

   return DefaultReturnValue;
}

/**---------------------------------------------------------------------------------------

   Set the velocities/positions of context2 to those of context1

   @param context1                 context1 
   @param context2                 context2 

   @return 0

   --------------------------------------------------------------------------------------- */

#ifdef FREE_ENERGY_BRANCH
static int _synchContexts( const OpenMMContext& context1, OpenMMContext& context2 ){
#else
static int _synchContexts( const Context& context1, Context& context2 ){
#endif

   // ---------------------------------------------------------------------------------------

   //static const char* methodName  = "\n_synchContexts: ";

   // ---------------------------------------------------------------------------------------

   const State state                       = context1.getState(State::Positions | State::Velocities);
   const std::vector<Vec3>& positions      = state.getPositions();
   const std::vector<Vec3>& velocities     = state.getVelocities();

   context2.setPositions( positions );
   context2.setVelocities( velocities );

   return DefaultReturnValue;
}

/**---------------------------------------------------------------------------------------
 * Get context
 * 
 * @param system          system
 * @param inputContext    input context -- if set, the newly created context is updated w/ positions & velocities
 * @param inputIntegrator input integrator for new context
 * @param platformName    name of platform( ReferencePlatform, CudaPlatform)
 * @param idString        diagnostic string (used in logging)
 * @param log             log file reference
 *
 * @return OpenMM context
 *
   --------------------------------------------------------------------------------------- */

Context* _getContext( System* system, Context* inputContext, Integrator* inputIntegrator, const std::string& platformName,
                      const std::string& idString, FILE* log ){

// ---------------------------------------------------------------------------------------

   static const std::string methodName      = "_getContext";

// ---------------------------------------------------------------------------------------

    // Create a context and initialize it.

    Context* context;
    ReferencePlatform referencePlatform;
    CudaPlatform gpuPlatform;
    if( platformName.compare( "ReferencePlatform" ) == 0 ){
       context = new Context( *system, *inputIntegrator, referencePlatform );
    } else {
       context = new Context( *system, *inputIntegrator, gpuPlatform );
    }

    if( log ){
       (void) fprintf( log, "%s Using Platform: %s\n", idString.c_str(), context->getPlatform().getName().c_str() );
       (void) fflush( log );
    }

    if( inputContext ){
       _synchContexts( *inputContext, *context );
    }

    return context;

}

/**---------------------------------------------------------------------------------------
      
   Get statistics of elements in array
   
   @param array               array to collect stats
   @param statistics          statistics of array
      index = 0   mean
      index = 1   stddev
      index = 2   min 
      index = 3   index of min value 
      index = 4   max
      index = 5   index of max value 
      index = 6   size of array

   @return DefaultReturnValue
      
   --------------------------------------------------------------------------------------- */

static int _getStatistics( const std::vector<double> & array,  std::vector<double> & statistics ){

   // ---------------------------------------------------------------------------------------

   static const char* methodName = "_getStatistics";

   static const int mean         = 0;
   static const int stddev       = 1;
   static const int min          = 2;
   static const int minIndex     = 3;
   static const int max          = 4;
   static const int maxIndex     = 5;
   static const int size         = 6;

   // ---------------------------------------------------------------------------------------

   // initialize stat array

   statistics.resize( 10 );
   for( unsigned int jj = 0; jj < statistics.size(); jj++ ){
      statistics[jj] = 0.0;
   }
   statistics[min] =  1.0e+30;
   statistics[max] = -1.0e+30;

   // collect stats

   int index       = 0;
   for( std::vector<double>::const_iterator ii = array.begin(); ii != array.end(); ii++ ){

      // first/second moments

      statistics[mean]     += *ii;
      statistics[stddev]   += (*ii)*(*ii);

      // min/max

      if( *ii < statistics[min] ){
         statistics[min]      = *ii;
         statistics[minIndex] = index;
      }
      if( *ii > statistics[max] ){
         statistics[max]      = *ii;
         statistics[maxIndex] = index;
      }
      index++;
   }

   // compute mean & std dev

   double arraySz      = (double) index;
   statistics[size]    = arraySz;
   if( index ){
      statistics[mean]   /= arraySz;
      statistics[stddev]  = statistics[stddev] - arraySz*statistics[mean]*statistics[mean];
      if( index > 1 ){
         statistics[stddev]  = std::sqrt( statistics[stddev] / ( arraySz - 1.0 ) );
      }
   }

   return DefaultReturnValue;
}

/**---------------------------------------------------------------------------------------
 * Check energy conservation
 * 
 * @param  context                context to run test on
 * @param  totalSimulationSteps   total number of simulation steps
 * @param  log                    log file reference
 *
 * @return DefaultReturnValue or ErrorReturnValue
 *
   --------------------------------------------------------------------------------------- */

#ifdef FREE_ENERGY_BRANCH
static int checkEnergyConservation( OpenMMContext& context, int totalSimulationSteps, FILE* log ) {
#else
static int checkEnergyConservation( Context& context, int totalSimulationSteps, FILE* log ) {
#endif
    
// ---------------------------------------------------------------------------------------

   static const std::string methodName             = "checkEnergyConservation";

   // initial T

   double initialTemperature                       = 300.0;

   // tolerance for thermostat

   double temperatureTolerance                     = 3.0;

   // tolerance for energy conservation test

   double energyTolerance                          = 0.05;

   std::string equilibrationIntegratorName         = "LangevinIntegrator";
   //std::string equilibrationIntegratorName         = "VerletIntegrator";
   int equilibrationTotalSteps                     = 10000;
   double equilibrationStepsBetweenReportsRatio    = 0.1;
   double equilibrationTimeStep                    = 0.002;
   double equilibrationFriction                    = 91.0;
   double equilibrationShakeTolerance              = 1.0e-05;
   double equilibrationErrorTolerance              = 1.0e-05;
   int equilibrationSeed                           = 1993;

   std::string simulationIntegratorName            = "VerletIntegrator";
   int simulationTotalSteps                        = totalSimulationSteps;
   double simulationStepsBetweenReportsRatio       = 0.01;
   double simulationTimeStep                       = 0.001;
   double simulationFriction                       = 91.0;
   double simulationShakeTolerance                 = 1.0e-06;
   double simulationErrorTolerance                 = 1.0e-05;
   int simulationSeed                              = 1993;

// ---------------------------------------------------------------------------------------

   int returnStatus                   = 0;
   clock_t     totalEquilibrationTime = 0;
   clock_t     totalSimulationTime    = 0;
   clock_t     cpuTime;

   int allTypes                       = State::Positions | State::Velocities | State::Forces | State::Energy;

   // set velocities based on temperature

   System& system                     = context.getSystem();
   int numberOfAtoms                  = system.getNumParticles();
   std::vector<Vec3> velocities; 
   //velocities.resize( numberOfAtoms );
   //_setVelocitiesBasedOnTemperature( system, velocities, initialTemperature, log );

   // get integrator for equilibration and context

   Integrator* integrator = _getIntegrator( equilibrationIntegratorName, equilibrationTimeStep,
                                            equilibrationFriction, initialTemperature,
                                            equilibrationShakeTolerance, equilibrationErrorTolerance, equilibrationSeed, log );

   if( log ){
      _printIntegratorInfo( integrator, log );
   }
#ifdef FREE_ENERGY_BRANCH
   OpenMMContext* equilibrationContext = _getContext( &system, &context, integrator, "CudaPlatform", "EquilibrationContext", log );
#else
   Context*       equilibrationContext = _getContext( &system, &context, integrator, "CudaPlatform", "EquilibrationContext", log );
#endif

   // equilibration loop

   int currentStep                        = 0;
   int equilibrationStepsBetweenReports   = static_cast<int>(static_cast<double>(equilibrationTotalSteps)*equilibrationStepsBetweenReportsRatio);
   if( equilibrationStepsBetweenReports < 1 )equilibrationStepsBetweenReports = 1;

   if( log ){
      (void) fprintf( log, "equilibrationTotalSteps=%d equilibrationStepsBetweenReports=%d ratio=%.4f\n", 
                      equilibrationTotalSteps, equilibrationStepsBetweenReports, equilibrationStepsBetweenReportsRatio);
      (void) fflush( log );
   }

   while( currentStep < equilibrationTotalSteps ){

      int nextStep = currentStep + equilibrationStepsBetweenReports;
      if( nextStep > equilibrationTotalSteps ){
         equilibrationStepsBetweenReports = equilibrationTotalSteps - currentStep;
      }

      // integrate

      cpuTime                 = clock();
      integrator->step(equilibrationStepsBetweenReports);
      totalEquilibrationTime += clock() - cpuTime;
      currentStep            += equilibrationStepsBetweenReports;

      // get energies and check for nans

      State state                            = equilibrationContext->getState( State::Energy );
   
      double kineticEnergy                   = state.getKineticEnergy();
      double potentialEnergy                 = state.getPotentialEnergy();
      double totalEnergy                     = kineticEnergy + potentialEnergy;
      if( log ){
         (void) fprintf( log, "%6d KE=%14.7e PE=%14.7e E=%14.7e\n", currentStep, kineticEnergy, potentialEnergy, totalEnergy );
         (void) fflush( log );
      }

#ifdef _MSC_VER
   #define isinf !_finite
   #define isnan _isnan
#endif
      if( isinf( totalEnergy ) || isnan( totalEnergy ) ){
         char buffer[1024];
         (void) sprintf( buffer, "%s Equilibration: nans detected at step %d -- aborting.\n", methodName.c_str(), currentStep );
         throwException(__FILE__, __LINE__, buffer );
         exit(-1);
      }

   }

   double kineticEnergy;
   double potentialEnergy;
   double totalEnergy;

   // report energies

   if( log ){

      State state                = equilibrationContext->getState( State::Energy );
      kineticEnergy              = state.getKineticEnergy();
      potentialEnergy            = state.getPotentialEnergy();
      totalEnergy                = kineticEnergy + potentialEnergy;
   
      double totalTime           = static_cast<double>(totalEquilibrationTime)/static_cast<double>(CLOCKS_PER_SEC);
      double timePerStep         = totalTime/static_cast<double>(equilibrationTotalSteps);
      double timePerStepPerAtom  = timePerStep/static_cast<double>(numberOfAtoms);
      (void) fprintf( log, "Final Equilibration energies: %6d  E=%14.7e [%14.7e %14.7e]  cpu time=%.3f time/step=%.3e time/step/atom=%.3e\n",
                      currentStep, (kineticEnergy + potentialEnergy), kineticEnergy, potentialEnergy,
                      totalTime, timePerStep, timePerStepPerAtom );
      (void) fflush( log );
   }

   // get simulation integrator & context

   Integrator* simulationIntegrator = _getIntegrator( simulationIntegratorName, simulationTimeStep,
                                                      simulationFriction, initialTemperature,
                                                      simulationShakeTolerance, simulationErrorTolerance,
                                                      simulationSeed, log );

   if( log ){
      _printIntegratorInfo( simulationIntegrator, log );
   }

#ifdef FREE_ENERGY_BRANCH
   OpenMMContext* simulationContext = _getContext( &system, equilibrationContext, simulationIntegrator, "CudaPlatform", "SimulationContext", log );
#else
   Context*       simulationContext = _getContext( &system, equilibrationContext, simulationIntegrator, "CudaPlatform", "SimulationContext", log );
#endif

   // create/initialize arrays used to track energies

   std::vector<double> stepIndexArray;
   std::vector<double> kineticEnergyArray;
   std::vector<double> potentialEnergyArray;
   std::vector<double> totalEnergyArray;

   State state                            = simulationContext->getState( State::Energy );
   kineticEnergy                          = state.getKineticEnergy();
   potentialEnergy                        = state.getPotentialEnergy();
   totalEnergy                            = kineticEnergy + potentialEnergy;

   stepIndexArray.push_back( 0.0 );
   kineticEnergyArray.push_back( kineticEnergy );
   potentialEnergyArray.push_back( potentialEnergy );
   totalEnergyArray.push_back( totalEnergy );

   // log

   if( log ){
      (void) fprintf( log, "Initial Simulation energies: E=%14.7e [%14.7e %14.7e]\n",
                      (kineticEnergy + potentialEnergy), kineticEnergy, potentialEnergy );
      (void) fflush( log );
   }

   /* -------------------------------------------------------------------------------------------------------------- */

   // prelude for simulation

   int simulationStepsBetweenReports   = static_cast<int>(static_cast<double>(simulationTotalSteps)*simulationStepsBetweenReportsRatio);
   if( simulationStepsBetweenReports < 1 )simulationStepsBetweenReports = 1;
   currentStep                       = 0;

   if( log ){
      (void) fprintf( log, "simulationTotalSteps=%d simulationStepsBetweenReports=%d ratio=%.4f\n", 
                      simulationTotalSteps, simulationStepsBetweenReports, simulationStepsBetweenReportsRatio );
      (void) fflush( log );
   }

   // main simulation loop

   while( currentStep < simulationTotalSteps ){

      // set step increment, perform integration, update step 

      int nextStep = currentStep + simulationStepsBetweenReports;
      if( nextStep > simulationTotalSteps ){
         simulationStepsBetweenReports = simulationTotalSteps - currentStep;
      }

      cpuTime              = clock();
      simulationIntegrator->step( simulationStepsBetweenReports );
      totalSimulationTime += clock() - cpuTime;

      currentStep         += simulationStepsBetweenReports;

      // record energies

      State state                            = simulationContext->getState( State::Energy );
      double kineticEnergy                   = state.getKineticEnergy();
      double potentialEnergy                 = state.getPotentialEnergy();
      double totalEnergy                     = kineticEnergy + potentialEnergy;

      stepIndexArray.push_back( (double) currentStep );
      kineticEnergyArray.push_back( kineticEnergy );
      potentialEnergyArray.push_back( potentialEnergy );
      totalEnergyArray.push_back( totalEnergy );

      // diagnostics & check for nans

      if( log ){
         (void) fprintf( log, "%6d KE=%14.7e PE=%14.7e E=%14.7e\n", currentStep, kineticEnergy, potentialEnergy, totalEnergy );
         (void) fflush( log );
      }

      if( isinf( totalEnergy ) || isnan( totalEnergy ) ){
         char buffer[1024];
         (void) sprintf( buffer, "%s Simulation: nans detected at step %d -- aborting.\n", methodName.c_str(), currentStep );
         throwException(__FILE__, __LINE__, buffer );
         exit(-1);
      }
   }

   state                                  = simulationContext->getState( State::Energy );
   kineticEnergy                          = state.getKineticEnergy();
   potentialEnergy                        = state.getPotentialEnergy();
   totalEnergy                            = kineticEnergy + potentialEnergy;

   // log times and energies

   if( log ){
      double totalTime           = static_cast<double>(totalSimulationTime)/static_cast<double>(CLOCKS_PER_SEC);
      double timePerStep         = totalTime/static_cast<double>(simulationTotalSteps);
      double timePerStepPerAtom  = timePerStep/static_cast<double>(numberOfAtoms);
      (void) fprintf( log, "Final Simulation: %6d  E=%14.7e [%14.7e %14.7e]  cpu time=%.3f time/step=%.3e time/step/atom=%.3e\n",
                      currentStep, (kineticEnergy + potentialEnergy), kineticEnergy, potentialEnergy,
                      totalTime, timePerStep, timePerStepPerAtom );
      (void) fflush( log );
   }

   // set dof

   double degreesOfFreedom  = static_cast<double>(3*numberOfAtoms - system.getNumConstraints() - 3 );
   double conversionFactor  = degreesOfFreedom*0.5*BOLTZ;
          conversionFactor  = 1.0/conversionFactor;

   // if Langevin or Brownian integrator, then check that temperature constant
   // else (Verlet integrator) check that energy drift is acceptable

   if( (simulationIntegratorName.compare( "LangevinIntegrator" ) == 0          ||
        simulationIntegratorName.compare( "VariableLangevinIntegrator" ) == 0  ||
        simulationIntegratorName.compare( "BrownianIntegrator" ) == 0) && numberOfAtoms > 0 ){

      // check that temperature constant

      // convert KE to temperature

      std::vector<double> temperature;
      for( std::vector<double>::const_iterator ii = kineticEnergyArray.begin(); ii != kineticEnergyArray.end(); ii++ ){
         temperature.push_back( (*ii)*conversionFactor );
      }

      // get temperature stats

      std::vector<double> temperatureStatistics;
      _getStatistics( temperature, temperatureStatistics );

      if( log ){
         (void) fprintf( log, "Simulation temperature results: mean=%14.7e stddev=%14.7e   min=%14.7e   %d max=%14.7e %d\n",
                         temperatureStatistics[0], temperatureStatistics[1], temperatureStatistics[2],
                         (int) (temperatureStatistics[3] + 0.001), temperatureStatistics[4],
                         (int) (temperatureStatistics[5] + 0.001) );

      }

      // check that <temperature> is within tolerance

      ASSERT_EQUAL_TOL( temperatureStatistics[0], initialTemperature, temperatureTolerance );

   } else {

      // total energy constant

      std::vector<double> statistics;
      _getStatistics( totalEnergyArray, statistics );

      std::vector<double> kineticEnergyStatistics;
      _getStatistics( kineticEnergyArray, kineticEnergyStatistics );
      double temperature  = kineticEnergyStatistics[0]*conversionFactor;
      double kT           = temperature*BOLTZ;

      // compute stddev in units of kT/dof/ns

      double stddevE      = statistics[1]/kT;
             stddevE     /= degreesOfFreedom;
             stddevE     /= simulationTotalSteps*simulationTimeStep*0.001;

      if( log ){
         (void) fprintf( log, "Simulation results: mean=%14.7e stddev=%14.7e  kT/dof/ns=%14.7e kT=%14.7e  min=%14.7e   %d max=%14.7e %d\n",
                         statistics[0], statistics[1], stddevE, kT, statistics[2], (int) (statistics[3] + 0.001), statistics[4], (int) (statistics[5] + 0.001) );
      }

      // check that energy fluctuation is within tolerance

      ASSERT_EQUAL_TOL( stddevE, 0.0, energyTolerance );

   }

   if( log ){
      (void) fprintf( log, "\n%s passed\n", methodName.c_str() );
      (void) fflush( log );
   }

   return returnStatus;

}

/**---------------------------------------------------------------------------------------

   Create context using content in parameter file (parameters/coordinates/velocities)

   @param parameterFileName    parameter file name
   @param forceFlag            flag controlling which forces are to be included
   @param platform             platform reference
   @param log                  FILE ptr; if NULL, diagnostic messages are not printed

   @return context

   --------------------------------------------------------------------------------------- */

#ifdef FREE_ENERGY_BRANCH
OpenMMContext* testSetup( std::string parameterFileName, int forceFlag, Platform& platform, std::vector<Vec3>& forces,
                          double* kineticEnergy, double* potentialEnergy, FILE* log ){
#else
Context* testSetup( std::string parameterFileName, int forceFlag, Platform& platform, std::vector<Vec3>& forces,
                    double* kineticEnergy, double* potentialEnergy, FILE* log ){
#endif

// ---------------------------------------------------------------------------------------

  static const std::string methodName      = "testSetup";
  double timeStep                          = 0.001; 
  double constraintTolerance               = 1.0e-05; 

// ---------------------------------------------------------------------------------------

   System* system  = new System();

   std::vector<Vec3> coordinates; 
   std::vector<Vec3> velocities; 

   // read parameters into system and coord/velocities into appropriate arrays

   Integrator* integrator =
      readParameterFile( parameterFileName, forceFlag, *system, coordinates, velocities,
                         forces, kineticEnergy, potentialEnergy, log );

#ifdef FREE_ENERGY_BRANCH
   OpenMMContext* context;
   context = new OpenMMContext( *system, *integrator, platform);
#else
   Context* context;
   context = new Context( *system, *integrator, platform);
#endif

   context->setPositions( coordinates );
   context->setVelocities( velocities );

   return context;

}

/**---------------------------------------------------------------------------------------

   Find stats for vec3

   @param array                 array 
   @param statVector              vector of stats 

   @return 0

   --------------------------------------------------------------------------------------- */

int compareForces( const std::vector<Vec3>& forceArray1, const std::string& f1Name, std::vector<double>& forceArray1Sum, std::vector<double>& forceArray1Stats,
                   const std::vector<Vec3>& forceArray2, const std::string& f2Name, std::vector<double>& forceArray2Sum, std::vector<double>& forceArray2Stats,
                   double* maxDelta, double* maxRelativeDelta, double forceTolerance, FILE* inputLog ){

// ---------------------------------------------------------------------------------------

  static const std::string methodName      = "compareForces";
  int PrintOn                              = 1; 

// ---------------------------------------------------------------------------------------

   FILE* log;
   if( PrintOn == 0 && inputLog ){
      log = NULL;
   } else {
      log = inputLog;
   } 

   if( log ){
      (void) fprintf( log, "%s\n", methodName.c_str() );
      (void) fflush( log );
   }   

   *maxDelta                           = -1.0e+30;
   *maxRelativeDelta                   = -1.0e+30;

   std::vector<double> forceArray1Norms;
   std::vector<double> forceArray2Norms;

   forceArray1Sum.resize( 3 );
   forceArray2Sum.resize( 3 );
   for( unsigned int ii = 0; ii < 3; ii++ ){
      forceArray1Sum[ii] = forceArray2Sum[ii] = 0.0;
   }

   for( unsigned int ii = 0; ii < forceArray2.size(); ii++ ){

      Vec3 f1                = forceArray1[ii];
      double normF1          = std::sqrt( (f1[0]*f1[0]) + (f1[1]*f1[1]) + (f1[2]*f1[2]) );
      forceArray1Norms.push_back( normF1 );
      forceArray1Sum[0]     += f1[0];
      forceArray1Sum[1]     += f1[1];
      forceArray1Sum[2]     += f1[2];

      Vec3 f2                = forceArray2[ii];
      double normF2          = std::sqrt( (f2[0]*f2[0]) + (f2[1]*f2[1]) + (f2[2]*f2[2]) );

      forceArray2Norms.push_back( normF2 );
      forceArray2Sum[0]     += f2[0];
      forceArray2Sum[1]     += f2[1];
      forceArray2Sum[2]     += f2[2];

      double delta           = std::sqrt( (f1[0]-f2[0])*(f1[0]-f2[0]) + (f1[1]-f2[1])*(f1[1]-f2[1]) + (f1[2]-f2[2])*(f1[2]-f2[2]) );
      double relativeDelta   = (delta*2.0)/(normF1+normF2);

      if( 0 && log && delta > forceTolerance ){
         (void) fprintf( log, "%5d delta=%13.6e relativeDelta=%13.6e  %s %13.6e [%14.7e %14.7e %14.7e]  %s %13.6e [%14.7e %14.7e %14.7e]\n", 
                         ii, delta, relativeDelta, f1Name.c_str(), normF1, f1[0], f1[1], f1[2],  f2Name.c_str(), normF2, f2[0], f2[1], f2[2] );
         (void) fflush( log );
      }
      if( log && *maxRelativeDelta < relativeDelta ){
         (void) fprintf( log, "%5d delta=%13.6e relativeDelta=%13.6e  %s %13.6e [%14.7e %14.7e %14.7e]  %s %13.6e [%14.7e %14.7e %14.7e]\n", 
                         ii, delta, relativeDelta, f1Name.c_str(), normF1, f1[0], f1[1], f1[2], f2Name.c_str(), normF2, f2[0], f2[1], f2[2] );
         (void) fflush( log );
         *maxRelativeDelta = relativeDelta;   
      }
      if( *maxDelta < delta ){
         *maxDelta = delta;
      }
   }

   findStatsForDouble( forceArray1Norms, forceArray1Stats );
   findStatsForDouble( forceArray2Norms, forceArray2Stats );

   return 0;
}

void testReferenceCudaForces( std::string parameterFileName, int forceFlag, FILE* inputLog ){

// ---------------------------------------------------------------------------------------

  static const std::string methodName      = "testReferenceCudaForces";
  int PrintOn                              = 1; 
  int compareParameterForces               = 0; 

  double forceTolerance                    = 0.01; 
  double energyTolerance                   = 0.01; 
  int numberOfSteps                        = 1; 
  int steps                                = 0; 
  int testOn                               = 0; 

// ---------------------------------------------------------------------------------------

   FILE* log;
   if( PrintOn == 0 && inputLog ){
      log = NULL;
   } else {
      log = inputLog;
   } 

   if( log ){
      (void) fprintf( log, "%s force tolerance=%.3e energy tolerance=%.3e step=%d\n",
                      methodName.c_str(), forceTolerance, energyTolerance, numberOfSteps );
      (void) fflush( log );
   }   

   ReferencePlatform referencePlatform;
   CudaPlatform cudaPlatform;

   double parameterKineticEnergy, parameterPotentialEnergy;

   std::vector<Vec3> parameterForces;
   std::vector<Vec3> parameterForces2;

#ifdef FREE_ENERGY_BRANCH
   OpenMMContext* referenceOpenMMContext = testSetup( parameterFileName, forceFlag, referencePlatform, 
                                                      parameterForces, &parameterKineticEnergy, &parameterPotentialEnergy, log );
   OpenMMContext* cudaOpenMMContext      = testSetup( parameterFileName, forceFlag,  cudaPlatform,
                                                      parameterForces2, &parameterKineticEnergy, &parameterPotentialEnergy, NULL );
#else
   Context* referenceOpenMMContext       = testSetup( parameterFileName, forceFlag, referencePlatform, 
                                                      parameterForces, &parameterKineticEnergy, &parameterPotentialEnergy, log );
   Context* cudaOpenMMContext            = testSetup( parameterFileName, forceFlag,  cudaPlatform,
                                                      parameterForces2, &parameterKineticEnergy, &parameterPotentialEnergy, log );
#endif

   Integrator& referenceIntegrator       = referenceOpenMMContext->getIntegrator();
   Integrator& cudaIntegrator            = cudaOpenMMContext->getIntegrator();

   // Run several steps and see if relative force difference is within tolerance
   
   for( int step = 0; step < numberOfSteps; step++ ){

      // pull info out of contexts

      int types                                       = State::Positions | State::Velocities | State::Forces | State::Energy;

      State cudaState                                 =      cudaOpenMMContext->getState( types );
      State referenceState                            = referenceOpenMMContext->getState( types );

      std::vector<Vec3> referenceCoordinates          = referenceState.getPositions();
      std::vector<Vec3> referenceVelocities           = referenceState.getVelocities();
      std::vector<Vec3> referenceForces               = referenceState.getForces();
      double referenceKineticEnergy                   = referenceState.getKineticEnergy();
      double referencePotentialEnergy                 = referenceState.getPotentialEnergy();

      std::vector<Vec3> cudaCoordinates               = cudaState.getPositions();
      std::vector<Vec3> cudaVelocities                = cudaState.getVelocities();
      std::vector<Vec3> cudaForces                    = cudaState.getForces();
      double cudaKineticEnergy                        = cudaState.getKineticEnergy();
      double cudaPotentialEnergy                      = cudaState.getPotentialEnergy();

      // diagnostics

      if( log ){
         //static const unsigned int maxPrint = MAX_PRINT;
         static const unsigned int maxPrint   = 1000000;

         // print x,y,z components separately, if formatType == 1
         // else print reference, cuda and parameter forces in blocks of 3

         static const unsigned int formatType = 1;

         (void) fprintf( log, "%s\n", methodName.c_str() );
         if( compareParameterForces ){
            (void) fprintf( log, "Kinetic   energies: r=%14.7e c=%14.7e, p=%14.7e\n", referenceKineticEnergy, cudaKineticEnergy, parameterKineticEnergy );
            (void) fprintf( log, "Potential energies: r=%14.7e c=%14.7e, p=%14.7e\n", referencePotentialEnergy, cudaPotentialEnergy, parameterPotentialEnergy );
            (void) fprintf( log, "Sample of forces: %u (r=reference, c=cuda, p=parameter) file forces\n", referenceForces.size() );
         } else {
            (void) fprintf( log, "Kinetic   energies: r=%14.7e c=%14.7e\n", referenceKineticEnergy, cudaKineticEnergy );
            (void) fprintf( log, "Potential energies: r=%14.7e c=%14.7e\n", referencePotentialEnergy, cudaPotentialEnergy );
            (void) fprintf( log, "Sample of forces: %u (r=reference, c=cuda) file forces\n", referenceForces.size() );
         }

         if( formatType == 1 ){
            (void) fprintf( log, "%s: atoms=%d [reference, cuda %s]\n", methodName.c_str(), referenceForces.size(), (compareParameterForces ? ", parameter" : "") );

            if( compareParameterForces ){
               for( unsigned int ii = 0; ii < referenceForces.size() && ii < maxPrint; ii++ ){
                  (void) fprintf( log, "%6u 0[%14.7e %14.7e %14.7e] 1[%14.7e %14.7e %14.7e] 2[%14.7e %14.7e %14.7e]\n", ii,
                                  referenceForces[ii][0], cudaForces[ii][0], parameterForces[ii][0],
                                  referenceForces[ii][1], cudaForces[ii][1], parameterForces[ii][1],
                                  referenceForces[ii][2], cudaForces[ii][2], parameterForces[ii][2] );
               }
               if( referenceForces.size() > maxPrint ){
                  for( unsigned int ii = referenceForces.size() - maxPrint; ii < referenceForces.size(); ii++ ){
                     (void) fprintf( log, "%6u 0[%14.7e %14.7e %14.7e] 1[%14.7e %14.7e %14.7e] 2[%14.7e %14.7e %14.7e]\n", ii,
                                     referenceForces[ii][0], cudaForces[ii][0], parameterForces[ii][0],
                                     referenceForces[ii][1], cudaForces[ii][1], parameterForces[ii][1],
                                     referenceForces[ii][2], cudaForces[ii][2], parameterForces[ii][2] );
                  }
               }
            } else {
               for( unsigned int ii = 0; ii < referenceForces.size() && ii < maxPrint; ii++ ){
                  (void) fprintf( log, "%6u 0[%14.7e %14.7e] 1[%14.7e %14.7e] 2[%14.7e %14.7e]\n", ii,
                                  referenceForces[ii][0], cudaForces[ii][0],
                                  referenceForces[ii][1], cudaForces[ii][1],
                                  referenceForces[ii][2], cudaForces[ii][2]  );
               }
               if( referenceForces.size() > maxPrint ){
                  for( unsigned int ii = referenceForces.size() - maxPrint; ii < referenceForces.size(); ii++ ){
                     (void) fprintf( log, "%6u 0[%14.7e %14.7e] 1[%14.7e %14.7e] 2[%14.7e %14.7e]\n", ii,
                                     referenceForces[ii][0], cudaForces[ii][0],
                                     referenceForces[ii][1], cudaForces[ii][1],
                                     referenceForces[ii][2], cudaForces[ii][2] );
                  }
               }
            }

         } else { 

            if( compareParameterForces ){
               for( unsigned int ii = 0; ii < referenceForces.size() && ii < maxPrint; ii++ ){
                  (void) fprintf( log, "%6u r[%14.7e %14.7e %14.7e] c[%14.7e %14.7e %14.7e] p[%14.7e %14.7e %14.7e]\n", ii,
                                  referenceForces[ii][0], referenceForces[ii][1], referenceForces[ii][2],
                                  cudaForces[ii][0], cudaForces[ii][1], cudaForces[ii][2],
                                  parameterForces[ii][0], parameterForces[ii][1], parameterForces[ii][2] );
               }
               if( referenceForces.size() > maxPrint ){
                  for( unsigned int ii = referenceForces.size() - maxPrint; ii < referenceForces.size(); ii++ ){
                     (void) fprintf( log, "%6u r[%14.7e %14.7e %14.7e] c[%14.7e %14.7e %14.7e] p[%14.7e %14.7e %14.7e]\n", ii,
                                     referenceForces[ii][0], referenceForces[ii][1], referenceForces[ii][2],
                                     cudaForces[ii][0], cudaForces[ii][1], cudaForces[ii][2],
                                     parameterForces[ii][0], parameterForces[ii][1], parameterForces[ii][2] );
                  }
               }
            } else {
               for( unsigned int ii = 0; ii < referenceForces.size() && ii < maxPrint; ii++ ){
                  (void) fprintf( log, "%6u r[%14.7e %14.7e %14.7e] c[%14.7e %14.7e %14.7e]\n", ii,
                                  referenceForces[ii][0], referenceForces[ii][1], referenceForces[ii][2],
                                  cudaForces[ii][0], cudaForces[ii][1], cudaForces[ii][2] );
               }
               if( referenceForces.size() > maxPrint ){
                  for( unsigned int ii = referenceForces.size() - maxPrint; ii < referenceForces.size(); ii++ ){
                     (void) fprintf( log, "%6u r[%14.7e %14.7e %14.7e] c[%14.7e %14.7e %14.7e]\n", ii,
                                     referenceForces[ii][0], referenceForces[ii][1], referenceForces[ii][2],
                                     cudaForces[ii][0], cudaForces[ii][1], cudaForces[ii][2] );
                  }
               }
            }
         }
      
      }

      // compare reference vs cuda forces

      double maxDeltaRefCud                          = -1.0e+30;
      double maxRelativeDeltaRefCud                  = -1.0e+30;
      double maxDeltaPrmCud                          = -1.0e+30;
      double maxRelativeDeltaPrmCud                  = -1.0e+30;

      std::vector<double> forceArray1Sum;
      std::vector<double> forceArray2Sum;
      std::vector<double> forceArray3Sum;

      std::vector<double> referenceForceStats;
      std::vector<double> cudaForceStats;
      std::vector<double> cudaForceStats1;
      std::vector<double> paramForceStats;

      compareForces( referenceForces, "fRef", forceArray1Sum, referenceForceStats,
                     cudaForces,      "fCud", forceArray2Sum, cudaForceStats, 
                     &maxDeltaRefCud, &maxRelativeDeltaRefCud, forceTolerance, log );
      
      (void) fflush( log );

      if( compareParameterForces ){

         // compare cuda & forces retreived from parameter file

         compareForces( parameterForces, "fPrm", forceArray3Sum, paramForceStats,
                        cudaForces,      "fCud", forceArray2Sum, cudaForceStats1,
                        &maxDeltaPrmCud, &maxRelativeDeltaPrmCud, forceTolerance, log );
      }

      if( log ){
         (void) fprintf( log, "Step=%d max delta=%14.7e maxRelativeDelta=%14.7e\n", step, maxDeltaRefCud, maxRelativeDeltaRefCud);
         (void) fprintf( log, "Reference force sum [%14.7e %14.7e %14.7e]\n", forceArray1Sum[0], forceArray1Sum[1], forceArray1Sum[2] );
         (void) fprintf( log, "Cuda      force sum [%14.7e %14.7e %14.7e]\n", forceArray2Sum[0], forceArray2Sum[1], forceArray2Sum[2] );
         if( compareParameterForces ){
            (void) fprintf( log, "Parameter force sum [%14.7e %14.7e %14.7e]\n", forceArray3Sum[0], forceArray3Sum[1], forceArray3Sum[2] );
         }

         (void) fprintf( log, "Reference force average=%14.7e stddev=%14.7e min=%14.7e at %6.0f max=%14.7e at %6.0f\n",
                         referenceForceStats[0], referenceForceStats[1], referenceForceStats[2], referenceForceStats[3],
                         referenceForceStats[4], referenceForceStats[5] );

         (void) fprintf( log, "     Cuda force average=%14.7e stddev=%14.7e min=%14.7e at %6.0f max=%14.7e at %6.0f\n",
                         cudaForceStats[0], cudaForceStats[1], cudaForceStats[2], cudaForceStats[3],
                         cudaForceStats[4], cudaForceStats[5] );

         if( compareParameterForces ){
            (void) fprintf( log, "    Param force average=%14.7e stddev=%14.7e min=%14.7e at %6.0f max=%14.7e at %6.0f\n",
                            paramForceStats[0], paramForceStats[1], paramForceStats[2], paramForceStats[3],
                            paramForceStats[4], paramForceStats[5] );
         }

         (void) fflush( log );
      }

      // check that relative force difference is small

      if( testOn ){
         ASSERT( maxRelativeDeltaRefCud < forceTolerance );

         // check energies

         ASSERT_EQUAL_TOL( referenceKineticEnergy,    cudaKineticEnergy,   energyTolerance );
         ASSERT_EQUAL_TOL( referencePotentialEnergy,  cudaPotentialEnergy, energyTolerance );
         if( compareParameterForces ){
            ASSERT_EQUAL_TOL( referencePotentialEnergy, parameterPotentialEnergy, energyTolerance );
         }
      }

/*
       double energy = state.getKineticEnergy()+state.getPotentialEnergy();
       if( PrintOn > 1 ){
          (void) fprintf( log, "%s %d e[%.5e %.5e] ke=%.5e pe=%.5e\n", 
                          methodName.c_str(), i, initialEnergy, energy, state.getKineticEnergy(), state.getPotentialEnergy() ); (void) fflush( log );
       }
       if( i == 1 ){
           initialEnergy = energy;
       } else if( i > 1 ){
           ASSERT_EQUAL_TOL(initialEnergy, energy, 0.5);
       }
*/
      if( steps ){
         cudaIntegrator.step( steps );
         _synchContexts( *cudaOpenMMContext, *referenceOpenMMContext );
      }

   }

   if( log ){
      if( testOn ){
         (void) fprintf( log, "\n%s tests passed\n", methodName.c_str() );
      } else {
         (void) fprintf( log, "\n%s tests off\n", methodName.c_str() );
      }
      (void) fflush( log );
   }
}

void testEnergyForcesConsistent( std::string parameterFileName, int forceFlag, double delta, FILE* inputLog ){

// ---------------------------------------------------------------------------------------

  static const std::string methodName      = "testEnergyForcesConsistent";
  int PrintOn                              = 1; 
  double tolerance                         = 0.01; 

// ---------------------------------------------------------------------------------------

   FILE* log;
   if( PrintOn == 0 && inputLog ){
      log = NULL;
   } else {
      log = inputLog;
   } 

   if( log ){
      (void) fprintf( log, "%s delta=%14.7e tolerance=%.3e\n",
                      methodName.c_str(), delta, tolerance );
      (void) fflush( log );
   }   

   ReferencePlatform referencePlatform;
   CudaPlatform cudaPlatform;

   double parameterKineticEnergy, parameterPotentialEnergy;

   std::vector<Vec3> parameterForces;
   std::vector<Vec3> parameterForces2;

#ifdef FREE_ENERGY_BRANCH
   OpenMMContext* referenceOpenMMContext = testSetup( parameterFileName, forceFlag, referencePlatform, 
                                                      parameterForces, &parameterKineticEnergy, &parameterPotentialEnergy, log );
   OpenMMContext* cudaOpenMMContext      = testSetup( parameterFileName, forceFlag,  cudaPlatform,
                                                      parameterForces2, &parameterKineticEnergy, &parameterPotentialEnergy, NULL );
#else
   Context* referenceOpenMMContext       = testSetup( parameterFileName, forceFlag, referencePlatform, 
                                                      parameterForces, &parameterKineticEnergy, &parameterPotentialEnergy, log );
   Context* cudaOpenMMContext            = testSetup( parameterFileName, forceFlag,  cudaPlatform,
                                                      parameterForces2, &parameterKineticEnergy, &parameterPotentialEnergy, NULL );
#endif

   if( log ){
      (void) fprintf( log, "%s Testing cuda platform\n", methodName.c_str() );
      (void) fflush( log );
   }   

   checkEnergyForceConsistent( *cudaOpenMMContext, delta, tolerance, log );

   if( log ){
      (void) fprintf( log, "%s Testing reference platform\n", methodName.c_str() );
      (void) fflush( log );
   }   

   checkEnergyForceConsistent( *referenceOpenMMContext, delta, tolerance, log );
}

void testEnergyConservation( std::string parameterFileName, int forceFlag, int totalSimulationSteps, FILE* inputLog ){

// ---------------------------------------------------------------------------------------

  static const std::string methodName      = "testEnergyConservation";
  int PrintOn                              = 1; 

// ---------------------------------------------------------------------------------------

   FILE* log;
   if( PrintOn == 0 && inputLog ){
      log = NULL;
   } else {
      log = inputLog;
   } 

   if( log ){
      (void) fprintf( log, "%s steps=%d\n", methodName.c_str(), totalSimulationSteps);
      (void) fflush( log );
   }   

   CudaPlatform cudaPlatform;

   double parameterKineticEnergy, parameterPotentialEnergy;

   std::vector<Vec3> parameterForces;
   std::vector<Vec3> parameterForces2;

#ifdef FREE_ENERGY_BRANCH
   OpenMMContext* cudaOpenMMContext      = testSetup( parameterFileName, forceFlag,  cudaPlatform,
                                                      parameterForces2, &parameterKineticEnergy, &parameterPotentialEnergy, log );
#else
   Context* cudaOpenMMContext            = testSetup( parameterFileName, forceFlag,  cudaPlatform,
                                                      parameterForces2, &parameterKineticEnergy, &parameterPotentialEnergy, log );
#endif

   if( log ){
      (void) fprintf( log, "%s Testing cuda platform\n", methodName.c_str() );
      (void) fflush( log );
   }   

   checkEnergyConservation( *cudaOpenMMContext, totalSimulationSteps, log );
}

/**---------------------------------------------------------------------------------------

   Print usage to screen and exit

   @param defaultParameterFileName   default parameter name

   @return 0

   --------------------------------------------------------------------------------------- */

int printUsage( std::string defaultParameterFileName ){

   (void) printf( "Usage:\nTestCudaUsingParameterFile\n" );

   (void) printf( "   -help this message\n" );

   (void) printf( "   -log log info to stdout for now\n" );
   (void) printf( "   -logFileName <log file name> (default=stdout)\n" );

   (void) printf( "\n" );
   (void) printf( "   -parameterFileName <parameter file name> (default=%s)\n", defaultParameterFileName.c_str() );

   (void) printf( "\n" );
   (void) printf( "   -checkEnergyForceConsistent do not check that force/energy are consistent\n" );
   (void) printf( "   +checkEnergyForceConsistent check that force/energy are consistent\n" );
   (void) printf( "   -delta <value> is size of perturbation used in numerically calculating force in checkEnergyForceConsistent test\n" );
   (void) printf( "                  default value is 1.0e-04\n" );

   (void) printf( "\n" );
   (void) printf( "   -checkEnergyConservation do not check that energy conservation\n" );
   (void) printf( "   +checkEnergyForceConsistent check energy conservation\n" );

   (void) printf( "\n" );
   (void) printf( "   -checkForces do not check that cuda/reference forces agree\n" );
   (void) printf( "   +checkForces check that cuda/reference forces agree\n" );
   (void) printf( "   +all include all forces (typically followed by -force entries)\n" );
   (void) printf( "   -force cr +force where force equals\n" );

   (void) printf( "   HarmonicBond \n" );
   (void) printf( "   HarmonicAngle\n" );
   (void) printf( "   PeriodicTorsion\n" );
   (void) printf( "   RBTorsion\n" );
   (void) printf( "   NB\n" );
   (void) printf( "   NbExceptions\n" );
   (void) printf( "   GbsaObc\n" );
   (void) printf( "   Note: multiple force entries are allowed.\n" );
   (void) printf( "   +force adds in the force; -force removes the force\n" );
   (void) printf( "   The arguments are case-insensitive\n" );
   (void) printf( "   The defaults is to include all forces represented in parameter file.\n" );
   (void) printf( "   Examples:\n\n" );
   (void) printf( "      To include all forces but the GBSA Obc force:\n" );
   (void) printf( "         TestCudaUsingParameterFile -parameterFileName %s +all -GbsaObc\n\n",  defaultParameterFileName.c_str() );
   (void) printf( "      To include only the harmonic bond force:\n" );
   (void) printf( "         TestCudaUsingParameterFile -parameterFileName %s +HarmonicBond\n\n",  defaultParameterFileName.c_str() );
   (void) printf( "      To include only the bond forces:\n" );
   (void) printf( "         TestCudaUsingParameterFile -parameterFileName %s +HarmonicBond +HarmonicAngle +PeriodicTorsion +RBTorsion\n\n", 
                            defaultParameterFileName.c_str() );

   exit(0);

   return 0;
}

/**---------------------------------------------------------------------------------------
 * Return forceEnum value if input argument matches one of the 
 * force names (HarmonicBond, HarmonicAngle, ...
 * The value returned is signed depending on whether the argument
 * contained a + or - (+HarmonicBond or -HarmonicBond)
 *
 * @param inputArgument   command-line argument
 * @param forceEnum       retrurn value
 *
 * @return 0 if argument is not a forceEnum argument
   --------------------------------------------------------------------------------------- */

int getForceOffset( const char* inputArgument, int* forceEnum ){

   int returnValue = -1;
   int sign;
   if( inputArgument[0] == '+' ){
      sign = 1;
   } else if( inputArgument[0] == '-' ){
      sign = -1;
   } else {
      return 0;
   }

   // skip over sign

   char* copyOfArgument = strdup( inputArgument );
   char* localArgument  = copyOfArgument + 1; 

#ifdef _MSC_VER
#define strcasecmp strcmp
#endif
   if( strcasecmp( localArgument, "HarmonicBond" ) == 0 ){
      returnValue =  HARMONIC_BOND_FORCE;
   } else  if( strcasecmp( localArgument, "HarmonicAngle" ) == 0 ){
      returnValue =  HARMONIC_ANGLE_FORCE;
   } else  if( strcasecmp( localArgument, "PeriodicTorsion" ) == 0 ){
      returnValue =  PERIODIC_TORSION_FORCE;
   } else  if( strcasecmp( localArgument, "RbTorsion" ) == 0 ){
      returnValue =  RB_TORSION_FORCE;
   } else  if( strcasecmp( localArgument, "Nb" ) == 0 ){
      returnValue =  NB_FORCE;
   } else  if( strcasecmp( localArgument, "NbExceptions" ) == 0 ){
      returnValue =  NB_EXCEPTION_FORCE;
   } else  if( strcasecmp( localArgument, "GbsaObc" ) == 0 ){
      returnValue =  GBSA_OBC_FORCE;
   }  

   free( copyOfArgument );
   if( returnValue > 0 ){
      *forceEnum = returnValue*sign;
      return 1;
   } else {
      return 0;
   }

}

// ---------------------------------------------------------------------------------------

int main( int numberOfArguments, char* argv[] ){

// ---------------------------------------------------------------------------------------

   static const std::string methodName               = "TestCudaFromFile";
   int checkForces                                   = 1;
   int checkEnergyForceConsistent                    = 0;
   int checkEnergyConservation                       = 0;

   // delta step used in EnergyForceConsistent calculations
   double delta                                      = 1.0e-04;

   // total number of steps to use in energy conservation test
   int simulationTotalSteps                          = 100000;

   FILE* log                                         = NULL;

// ---------------------------------------------------------------------------------------

   std::string defaultParameterFileName     = "OpenParameters.txt";

   if( numberOfArguments < 2 ){
      printUsage( defaultParameterFileName );
   }

   std::string parameterFileName            = defaultParameterFileName;
   int forceFlag                            = 0;
   int logFileNameIndex                     = -1;
   int forceEnum;

   // parse arguments

   for( int ii = 1; ii < numberOfArguments; ii++ ){
      if( strcasecmp( argv[ii], "-parameterFileName" ) == 0 ){
         parameterFileName          = argv[ii+1];
         ii++;
      } else if( strcasecmp( argv[ii], "-logFileName" ) == 0 ){
         logFileNameIndex           = ii + 1;
         ii++;
      } else if( strcasecmp( argv[ii], "-checkForces" ) == 0 ){
         checkForces                = 0;
      } else if( strcasecmp( argv[ii], "+checkForces" ) == 0 ){
         checkForces                = 1;
      } else if( strcasecmp( argv[ii], "-checkEnergyForceConsistent" ) == 0 ){
         checkEnergyForceConsistent = 0;
      } else if( strcasecmp( argv[ii], "+checkEnergyForceConsistent" ) == 0 ){
         checkEnergyForceConsistent = 1;
      } else if( strcasecmp( argv[ii], "-delta" ) == 0 ){
         delta                      = atof( argv[ii+1] );
         ii++;

      } else if( strcasecmp( argv[ii], "-checkEnergyConservation" ) == 0 ){
         checkEnergyConservation = 0;
      } else if( strcasecmp( argv[ii], "+checkEnergyConservation" ) == 0 ){
         checkEnergyConservation = 1;

      } else if( strcasecmp( argv[ii], "-simulationTotalSteps" ) == 0 ){
         simulationTotalSteps   = atoi( argv[ii+1] );
         ii++;
      } else if( strcasecmp( argv[ii], "+all" ) == 0 ){
         forceFlag += HARMONIC_BOND_FORCE + HARMONIC_ANGLE_FORCE + PERIODIC_TORSION_FORCE +
                      RB_TORSION_FORCE    + NB_FORCE             + NB_EXCEPTION_FORCE     +
                      GBSA_OBC_FORCE;
      } else if( strcasecmp( argv[ii], "-log" ) == 0 ){
         log = stdout;
      } else if( strcasecmp( argv[ii], "-help" ) == 0 ){
         printUsage( defaultParameterFileName );
      } else if( getForceOffset( argv[ii], &forceEnum ) ){
         forceFlag += forceEnum;
      } else {
         (void) printf( "Argument=<%s> not recognized -- aborting\n", argv[ii] );
         exit(-1);
      }
   }

   // open log file

   if( log && logFileNameIndex > -1 ){
#ifdef _MSC_VER
         fopen_s( &log, argv[logFileNameIndex], "w" );
#else
         log = fopen( argv[logFileNameIndex], "w" );
#endif
   }

   // log info

   if( log ){
      (void) fprintf( log, "Input arguments:\n" );
      for( int ii = 1; ii < numberOfArguments; ii++ ){
         (void) fprintf( log, "%3d arg=<%s>\n", ii, argv[ii] );
      }
      (void) fprintf( log, "parameter file=<%s>\n", parameterFileName.c_str() );
      (void) fprintf( log, "forceFlag=%d\n", forceFlag );
      (void) fprintf( log, "checkEnergyForceConsistent=%d delta=%14.4e\n", checkEnergyForceConsistent, delta );
      (void) fprintf( log, "checkEnergyConservation=%d steps=%d\n",    checkEnergyConservation, simulationTotalSteps);
      if( forceFlag ){
         (void) fprintf( log, "HarmonicBond     %1d\n", ((forceFlag & HARMONIC_BOND_FORCE)    ? 1 : 0) );
         (void) fprintf( log, "HarmonicAngle    %1d\n", ((forceFlag & HARMONIC_ANGLE_FORCE)   ? 1 : 0) );
         (void) fprintf( log, "PeriodicTorsion  %1d\n", ((forceFlag & PERIODIC_TORSION_FORCE) ? 1 : 0) );
         (void) fprintf( log, "RB Torsion       %1d\n", ((forceFlag & RB_TORSION_FORCE)       ? 1 : 0) );
         (void) fprintf( log, "NB               %1d\n", ((forceFlag & NB_FORCE)               ? 1 : 0) );
         (void) fprintf( log, "NB Exceptions    %1d\n", ((forceFlag & NB_EXCEPTION_FORCE)     ? 1 : 0) );
         (void) fprintf( log, "GBSA OBC         %1d\n", ((forceFlag & GBSA_OBC_FORCE)         ? 1 : 0) );
      } 
      (void) fflush( log );
   }

   // do tests

//angleTest( log );
//exit(0);

   // check forces

   if( checkForces ){
      try {
         testReferenceCudaForces( parameterFileName, forceFlag, log );
       } catch( const exception& e ){
         (void) fprintf( stderr, "Exception checkForces %s %.s\n", methodName.c_str(),  e.what() ); (void) fflush( stderr );
         return 1;
      }   
   }

   // check energy/force consistent

   if( checkEnergyForceConsistent ){
      try {
         testEnergyForcesConsistent( parameterFileName, forceFlag, delta, log );
       } catch( const exception& e ){
         (void) fprintf( stderr, "Exception checkEnergyForceConsistent %s %.s\n", methodName.c_str(),  e.what() ); (void) fflush( stderr );
         return 1;
      }   
   }

   // check energy conservation or thermal stability

   if( checkEnergyConservation ){
      try {
         testEnergyConservation( parameterFileName, forceFlag, simulationTotalSteps, log );
       } catch( const exception& e ){
         (void) fprintf( stderr, "Exception checkEnergyConservation %s %.s\n", methodName.c_str(),  e.what() ); (void) fflush( stderr );
         return 1;
      }   
   }

   if( log ){
      (void) fprintf( log, "\n%s done\n", methodName.c_str() ); (void) fflush( log );
   }

   return 0;
}
