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

#include "openmm/Context.h"

#include "openmm/HarmonicBondForce.h"
#include "openmm/HarmonicAngleForce.h"
#include "openmm/PeriodicTorsionForce.h"
#include "openmm/RBTorsionForce.h"
#include "openmm/GBSAOBCForce.h"
#include "openmm/GBVIForce.h"
#include "openmm/NonbondedForce.h"
#include "openmm/CMMotionRemover.h"
#include "openmm/System.h"
#include "openmm/LangevinIntegrator.h"
#include "openmm/VariableLangevinIntegrator.h"
#include "openmm/VerletIntegrator.h"
#include "openmm/VariableVerletIntegrator.h"
#include "openmm/BrownianIntegrator.h"

// free-energy plugin includes
#define	INCLUDE_FREE_ENERGY_PLUGIN
#ifdef INCLUDE_FREE_ENERGY_PLUGIN
#include "OpenMMFreeEnergy.h"
#include "openmm/freeEnergyKernels.h"
#include "ReferenceFreeEnergyKernelFactory.h"
#include "CudaFreeEnergyKernelFactory.h"
#endif

#include <ctime>
#include <vector>
#include <cfloat>
#include <cstring>
#include <cstdlib>
#include <typeinfo>
#include <sstream>

#ifdef _MSC_VER
   #define isinf !_finite
   #define isnan _isnan
#endif

// max entries to print for default output
#define MAX_PRINT 5

// force names

std::string HARMONIC_BOND_FORCE             = "HarmonicBond";
std::string HARMONIC_ANGLE_FORCE            = "HarmonicAngle"; 
std::string PERIODIC_TORSION_FORCE          = "PeriodicTorsion";
std::string RB_TORSION_FORCE                = "RbTorsion";

std::string NB_FORCE                        = "Nb";
std::string NB_SOFTCORE_FORCE               = "NbSoftcore";

std::string NB_EXCEPTION_FORCE              = "NbException";
std::string NB_EXCEPTION_SOFTCORE_FORCE     = "NbExceptionSoftcore";

std::string GBSA_OBC_FORCE                  = "Obc";
std::string GBSA_OBC_SOFTCORE_FORCE         = "ObcSoftcore";

std::string GBVI_FORCE                      = "GBVI";
std::string GBVI_SOFTCORE_FORCE             = "GBVISoftcore";

#define BOLTZMANN                     (1.380658e-23)               /* (J/K) */
#define AVOGADRO                      (6.0221367e23)               /* ()    */
#define RGAS                          (BOLTZMANN*AVOGADRO)         /* (J/(mol K))  */
#define BOLTZ                         (RGAS/1.0e+03)               /* (kJ/(mol K)) */

using namespace OpenMM;
using namespace std;

// the following are used in parsing parameter file

typedef std::vector<std::string> StringVector;
typedef StringVector::iterator StringVectorI;
typedef StringVector::const_iterator StringVectorCI;

typedef std::vector<std::vector<double> > VectorOfVectors;
typedef VectorOfVectors::iterator VectorOfVectorsI;
typedef VectorOfVectors::const_iterator VectorOfVectorsCI;

typedef std::map< std::string, VectorOfVectors > MapStringVectorOfVectors;
typedef MapStringVectorOfVectors::iterator MapStringVectorOfVectorsI;
typedef MapStringVectorOfVectors::const_iterator MapStringVectorOfVectorsCI;

typedef std::map< std::string, std::string > MapStringString;
typedef MapStringString::iterator MapStringStringI;
typedef MapStringString::const_iterator MapStringStringCI;

typedef std::map< std::string, int > MapStringInt;
typedef MapStringInt::iterator MapStringIntI;
typedef MapStringInt::const_iterator MapStringIntCI;

/* --------------------------------------------------------------------------------------- */
// internal routines

char* readLine( FILE* filePtr, StringVector& tokens, int* lineCount, FILE* log );
int readVec3( FILE* filePtr, const StringVector& tokens, std::vector<Vec3>& coordinates, int* lineCount, FILE* log );

/* --------------------------------------------------------------------------------------- */

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

/**
 * Write vec3 array to file
 *
 * @param    filePtr            file ptr to output data
 * @param    vect3Array         array to output
 *
 * @return   0
 */

static int writeFileVec3( FILE* filePtr, const std::vector<Vec3>& vect3Array ){

    for( unsigned int ii = 0; ii < vect3Array.size(); ii++ ){
       (void) fprintf( filePtr, "%8d  %14.7e %14.7e %14.7e\n", ii, 
                       vect3Array[ii][0], vect3Array[ii][1], vect3Array[ii][2] );
    }   

    return 0;
}

/**---------------------------------------------------------------------------------------

 * Write context to file
 *
 * @param    fileName           file name
 * @param    context            OpenMM::Context used to get current positions
 * @param    stateFlag          State::Positions | State::Velocities | State::Forces  | State::Energy
 * @param    log                log file
 *
 * @return   0

   --------------------------------------------------------------------------------------- */

static int writeContextToFile( std::string fileName, Context& context, int stateFlag, FILE* log ){

// ---------------------------------------------------------------------------------------

   static const std::string methodName      = "writeContextToFile";
   
// ---------------------------------------------------------------------------------------

   // open file

   FILE* filePtr;
#ifdef _MSC_VER
    fopen_s( &filePtr, fileName.c_str(), "w" );
#else
    filePtr = fopen( fileName.c_str(), "w" );
#endif

    if( filePtr == NULL ){
        char buffer[1024];
        (void) sprintf( buffer, "%s file=<%s> not opened.\n", methodName.c_str(), fileName.c_str());
        throwException(__FILE__, __LINE__, buffer );
        exit(-1);
    } else if( log ){
        (void) fprintf( log, "%s opened file %s.\n", methodName.c_str(), fileName.c_str());
    }
      
    State state = context.getState( stateFlag );

    if( stateFlag && State::Positions ){
       std::vector<Vec3> positions  = state.getPositions();
       (void) fprintf( filePtr, "Positions %u\n", positions.size() );
       writeFileVec3( filePtr, positions );
    }

    if( stateFlag && State::Velocities ){
       std::vector<Vec3> velocities = state.getVelocities();
       (void) fprintf( filePtr, "Velocities %u\n", velocities.size() );
       writeFileVec3( filePtr, velocities );
    }

    if( stateFlag && State::Forces ){
       std::vector<Vec3> forces     = state.getForces();
       (void) fprintf( filePtr, "Forces %u\n", forces.size() );
       writeFileVec3( filePtr, forces );
    }

    if( stateFlag && State::Energy ){
       (void) fprintf( filePtr, "KineticEnergy %14.7e\n", state.getKineticEnergy() );
       (void) fprintf( filePtr, "PotentialEnergy %14.7e\n", state.getPotentialEnergy() );
    }

    (void) fclose( filePtr );

    return 0;
}

/**---------------------------------------------------------------------------------------

 * Read context from file
 *
 * @param    fileName           file name
 * @param    context            OpenMM::Context to update
 * @param    stateFlag          State::Positions | State::Velocities | State::Forces  | State::Energy
 * @param    log                log file
 *
 * @return   0

   --------------------------------------------------------------------------------------- */

static int readContextFromFile( std::string fileName, Context& context, int stateFlag, FILE* log ){

// ---------------------------------------------------------------------------------------

   static const std::string methodName      = "readContextFromFile";
   
// ---------------------------------------------------------------------------------------

   // open file

   FILE* filePtr;
#ifdef _MSC_VER
    fopen_s( &filePtr, fileName.c_str(), "r" );
#else
    filePtr = fopen( fileName.c_str(), "r" );
#endif

    if( filePtr == NULL ){
        char buffer[1024];
        (void) sprintf( buffer, "%s file=<%s> not opened.\n", methodName.c_str(), fileName.c_str());
        throwException(__FILE__, __LINE__, buffer );
        exit(-1);
    } else if( log ){
        (void) fprintf( log, "%s opened file %s.\n", methodName.c_str(), fileName.c_str());
    }
      
    std::vector<Vec3> coordinates; 
    std::vector<Vec3> velocities; 
    std::vector<Vec3> forces; 
    double kineticEnergy, potentialEnergy;
    std::string version;

    int lineCount  = 0;
    char* isNotEof = "1";

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
                    kineticEnergy    = value;
                } else {
                    potentialEnergy  = value;
                }
            } else {
                char buffer[1024];
                (void) sprintf( buffer, "Field=<%s> not recognized at line=%d\n", field.c_str(), lineCount );
                throwException(__FILE__, __LINE__, buffer );
                exit(-1);
            }
        }
     }
  
    // close file

    (void) fclose( filePtr );

    System& system = context.getSystem();
    if( stateFlag & State::Positions ){
        if( system.getNumParticles() != coordinates.size() ){
            char buffer[1024];
            (void) sprintf( buffer, "%s: number of positions=%u does not agree w/ number in system=%d\n",
                            methodName.c_str(), coordinates.size(), system.getNumParticles() );
            throwException(__FILE__, __LINE__, buffer );
            exit(-1);
        } else if( log ){
            (void) fprintf( log, "%s setting positions from context file.\n", methodName.c_str() );
        }
        context.setPositions( coordinates );
    }

    if( stateFlag & State::Velocities ){
        if( system.getNumParticles() != velocities.size() ){
            char buffer[1024];
            (void) sprintf( buffer, "%s: number of velocities=%u does not agree w/ number in system=%d\n",
                            methodName.c_str(), velocities.size(), system.getNumParticles() );
            throwException(__FILE__, __LINE__, buffer );
            exit(-1);
        } else if( log ){
            (void) fprintf( log, "%s setting velocities from context file.\n", methodName.c_str() );
        }
        context.setVelocities( velocities );
    }

    return 0;
}

/**---------------------------------------------------------------------------------------

 * Check constraints
 *
 * @param    context            OpenMM::Context used to get current positions
 * @param    system             OpenMM::System to be created
 * @param    tolerance          constraint tolerance
 * @param    maxViolation       output max constraint violation
 * @param    log                log file
 *
 * @return   return number of violations

   --------------------------------------------------------------------------------------- */

static int checkConstraints( const Context& context, const System& system, double tolerance, 
                             double* maxViolation, FILE* log ) {
    
// ---------------------------------------------------------------------------------------
      
    int totalPrints                    = 0;
	 int violations                     = 0;

// ---------------------------------------------------------------------------------------

    *maxViolation                      = -1.0e-10;
	 State state                        = context.getState(State::Positions);
	 const std::vector<Vec3>& pos       = state.getPositions();
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
          if( delta > *maxViolation ){
             *maxViolation = delta;
          }
			 if( log && totalPrints++ < 10 ){
			    (void) fprintf( log, "CnstrViolation: %6d %6d particles[%6d %6d] delta=%10.3e d[%12.5e %12.5e] \n",
				                 ii, violations, particle1, particle2, delta, distance, actualDistance );
          }
       }
	 }
	 if( log && violations ){
       (void) fprintf( log, "CnstrViolation: total violations=%d out of %d constraints; maxViolation=%13.6e tolerance=%.3e.\n",
		                 violations, system.getNumConstraints(), *maxViolation, tolerance );
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

static int sumForces( const Context& context, const System& system, double forceSum[3], int step, FILE* log ){  

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

static int checkKineticEnergy( const Context& context, const System& system, double temperature, int step, FILE* log ){  

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

   // static const std::string methodName = "tokenizeString";

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

char* readLine( FILE* filePtr, StringVector& tokens, int* lineCount, FILE* log ){

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

   Read vector of double vectors

   @param filePtr              file pointer to parameter file
   @param tokens               array of strings from first line of parameter file for this block of parameters
   @param vectorOfVectors      output of vector of vectors
   @param lineCount            used to track line entries read from parameter file
   @param typeName             id of entries being read
   @param log                  log file pointer -- may be NULL

   @return number of entries read

   --------------------------------------------------------------------------------------- */

static int readVectorOfVectors( FILE* filePtr, const StringVector& tokens, std::vector< std::vector<double> >& vectorOfVectors, 
                                int* lineCount, std::string typeName, FILE* log ){

// ---------------------------------------------------------------------------------------

   static const std::string methodName      = "readVec3";
   
// ---------------------------------------------------------------------------------------

   if( tokens.size() < 1 ){
      char buffer[1024];
      (void) sprintf( buffer, "%s no Coordinates terms entry???\n", methodName.c_str() );
      throwException(__FILE__, __LINE__, buffer );
      exit(-1);
   }

   int numberToRead = atoi( tokens[1].c_str() );
   if( log ){
      (void) fprintf( log, "%s number of %s to read: %d\n", methodName.c_str(), typeName.c_str(), numberToRead );
      (void) fflush( log );
   }
   for( int ii = 0; ii < numberToRead; ii++ ){
      StringVector lineTokens;
      char* isNotEof = readLine( filePtr, lineTokens, lineCount, log );
      int tokenIndex = 0;
      if( lineTokens.size() > 1 ){
         int index            = atoi( lineTokens[tokenIndex++].c_str() );
         std::vector<double> nextEntry;
         for( unsigned int jj = 1; jj < lineTokens.size(); jj++ ){
             double value = atof( lineTokens[jj].c_str() );
             nextEntry.push_back( value );
         }
         vectorOfVectors.push_back( nextEntry );
      } else {
         char buffer[1024];
         (void) sprintf( buffer, "%s %s tokens incomplete at line=%d\n", methodName.c_str(), typeName.c_str(), *lineCount );
         throwException(__FILE__, __LINE__, buffer );
         exit(-1);
      }
   }

   // diagnostics

   if( log ){
      static const unsigned int maxPrint = MAX_PRINT;
      unsigned int   arraySize           = vectorOfVectors.size();
      (void) fprintf( log, "%s: sample of %s size=%u\n", methodName.c_str(), typeName.c_str(), arraySize );
      for( unsigned int ii = 0; ii < vectorOfVectors.size(); ii++ ){
         (void) fprintf( log, "%6u [", ii );
         for( unsigned int jj = 0; jj < vectorOfVectors[ii].size(); jj++ ){
            (void) fprintf( log, "%14.7e ", vectorOfVectors[ii][jj] );
         }
         (void) fprintf( log, "]\n" );

         // skip to end

         if( ii == maxPrint && (arraySize - maxPrint) > ii ){
            ii = arraySize - maxPrint - 1;
            if( ii < maxPrint )ii = maxPrint;
         } 
      }
   }

   return static_cast<int>(vectorOfVectors.size());
}

/**---------------------------------------------------------------------------------------
 * Set field if in map
 * 
 * @param  argumentMap            map to check
 * @param  fieldToCheck           key
 * @param  fieldToSet             field to set
 *
 * @return 1 if argument set, else 0
 *
   --------------------------------------------------------------------------------------- */

static int setStringFromMap( MapStringString& argumentMap, std::string fieldToCheck, std::string& fieldToSet ){

// ---------------------------------------------------------------------------------------

   static const std::string methodName             = "setStringFromMap";

// ---------------------------------------------------------------------------------------

   MapStringStringCI check = argumentMap.find( fieldToCheck );
   if( check != argumentMap.end() ){
      fieldToSet = (*check).second; 
      return 1;
   }
   return 0;
}

/**---------------------------------------------------------------------------------------
 * Set field if in map
 * 
 * @param  argumentMap            map to check
 * @param  fieldToCheck           key
 * @param  fieldToSet             field to set
 *
 * @return 1 if argument set, else 0
 *
   --------------------------------------------------------------------------------------- */

static int setIntFromMap( MapStringString& argumentMap, std::string fieldToCheck, int& fieldToSet ){

// ---------------------------------------------------------------------------------------

   static const std::string methodName             = "setIntFromMap";

// ---------------------------------------------------------------------------------------

   MapStringStringCI check = argumentMap.find( fieldToCheck );
   if( check != argumentMap.end() ){
      fieldToSet = atoi( (*check).second.c_str() ); 
      return 1;
   }
   return 0;
}

/**---------------------------------------------------------------------------------------

 * Set field if in map
 * 
 * @param  argumentMap            map to check
 * @param  fieldToCheck           key
 * @param  fieldToSet             field to set
 *
 * @return 1 if argument set, else 0
 *
   --------------------------------------------------------------------------------------- */

static int setDoubleFromMap( MapStringString& argumentMap, std::string fieldToCheck, double& fieldToSet ){

// ---------------------------------------------------------------------------------------

   static const std::string methodName             = "setDoubleFromMap";

// ---------------------------------------------------------------------------------------

   MapStringStringCI check = argumentMap.find( fieldToCheck );
   if( check != argumentMap.end() ){
      fieldToSet = atof( (*check).second.c_str() ); 
      return 1;
   }
   return 0;
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
      int tokenIndex = 0;
      if( lineTokens.size() >= 1 ){
         int index   = atoi( lineTokens[tokenIndex++].c_str() );
         double mass = atof( lineTokens[tokenIndex++].c_str() );
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
      for( unsigned int ii = 0; ii < arraySize; ii++ ){
         (void) fprintf( log, "%6u %14.7e \n", ii, system.getParticleMass( ii ) );
         if( ii == maxPrint ){
            ii = arraySize - maxPrint;
            if( ii < maxPrint )ii = maxPrint;
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

static int readHarmonicBondForce( FILE* filePtr, MapStringInt& forceMap, const StringVector& tokens, System& system, int* lineCount, FILE* log ){

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
   MapStringIntI forceActive    = forceMap.find( HARMONIC_BOND_FORCE );
   if( forceActive != forceMap.end() && (*forceActive).second ){
      system.addForce( bondForce );
      if( log ){
         (void) fprintf( log, "harmonic bond force is being included.\n" );
      }
   } else if( log ){
      (void) fprintf( log, "harmonic bond force is not being included.\n" );
   }

   int numberOfBonds            = atoi( tokens[1].c_str() );
   if( log ){
      (void) fprintf( log, "%s number of HarmonicBondForce terms=%d\n", methodName.c_str(), numberOfBonds );
   }
   for( int ii = 0; ii < numberOfBonds; ii++ ){
      StringVector lineTokens;
      char* isNotEof = readLine( filePtr, lineTokens, lineCount, log );
      int tokenIndex = 0;
      if( lineTokens.size() > 4 ){
         int index      = atoi( lineTokens[tokenIndex++].c_str() );
         int particle1  = atoi( lineTokens[tokenIndex++].c_str() );
         int particle2  = atoi( lineTokens[tokenIndex++].c_str() );
         double length  = atof( lineTokens[tokenIndex++].c_str() );
         double k       = atof( lineTokens[tokenIndex++].c_str() );
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
      for( unsigned int ii = 0; ii < arraySize; ii++ ){
         int particle1, particle2;
         double length, k;
         bondForce->getBondParameters( ii, particle1, particle2, length, k ); 
         (void) fprintf( log, "%8d %8d %8d %14.7e %14.7e\n", ii, particle1, particle2, length, k);
         if( ii == maxPrint ){
            ii = arraySize - maxPrint;
            if( ii < maxPrint )ii = maxPrint;
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

static int readHarmonicAngleForce( FILE* filePtr, MapStringInt& forceMap, const StringVector& tokens,
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
   MapStringIntI forceActive     = forceMap.find( HARMONIC_ANGLE_FORCE );
   if( forceActive != forceMap.end() && (*forceActive).second ){
      system.addForce( bondForce );
      if( log ){
         (void) fprintf( log, "harmonic angle force is being included.\n" );
      }
   } else if( log ){
      (void) fprintf( log, "harmonic angle force is not being included.\n" );
   }

   int numberOfAngles            = atoi( tokens[1].c_str() );
   if( log ){
      (void) fprintf( log, "%s number of HarmonicAngleForce terms=%d\n", methodName.c_str(), numberOfAngles );
   }
   for( int ii = 0; ii < numberOfAngles; ii++ ){
      StringVector lineTokens;
      char* isNotEof = readLine( filePtr, lineTokens, lineCount, log );
      int tokenIndex = 0;
      if( lineTokens.size() > 5 ){
         int index      = atoi( lineTokens[tokenIndex++].c_str() );
         int particle1  = atoi( lineTokens[tokenIndex++].c_str() );
         int particle2  = atoi( lineTokens[tokenIndex++].c_str() );
         int particle3  = atoi( lineTokens[tokenIndex++].c_str() );
         double angle   = atof( lineTokens[tokenIndex++].c_str() );
         double k       = atof( lineTokens[tokenIndex++].c_str() );
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
      for( unsigned int ii = 0; ii < arraySize; ii++ ){
         int particle1, particle2, particle3;
         double angle, k;
         bondForce->getAngleParameters( ii, particle1, particle2, particle3, angle, k ); 
         (void) fprintf( log, "%8d %8d %8d %8d %14.7e %14.7e\n", ii, particle1, particle2, particle3, angle, k);
         if( ii == maxPrint ){
            ii = arraySize - maxPrint;
            if( ii < maxPrint )ii = maxPrint;
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

static int readPeriodicTorsionForce( FILE* filePtr, MapStringInt& forceMap, const StringVector& tokens, System& system, int* lineCount, FILE* log ){

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
   MapStringIntI forceActive    = forceMap.find( PERIODIC_TORSION_FORCE );
   if( forceActive != forceMap.end() && (*forceActive).second ){
      system.addForce( bondForce );
      if( log ){
         (void) fprintf( log, "periodic torsion force is being included.\n" );
      }
   } else if( log ){
      (void) fprintf( log, "periodic torsion force is not being included.\n" );
   }

   int numberOfTorsions            = atoi( tokens[1].c_str() );
   if( log ){
      (void) fprintf( log, "%s number of PeriodicTorsionForce terms=%d\n", methodName.c_str(), numberOfTorsions );
   }
   for( int ii = 0; ii < numberOfTorsions; ii++ ){
      StringVector lineTokens;
      char* isNotEof = readLine( filePtr, lineTokens, lineCount, log );
      int tokenIndex = 0;
      if( lineTokens.size() > 7 ){
         int index       = atoi( lineTokens[tokenIndex++].c_str() );
         int particle1   = atoi( lineTokens[tokenIndex++].c_str() );
         int particle2   = atoi( lineTokens[tokenIndex++].c_str() );
         int particle3   = atoi( lineTokens[tokenIndex++].c_str() );
         int particle4   = atoi( lineTokens[tokenIndex++].c_str() );
         int periodicity = atoi( lineTokens[tokenIndex++].c_str() );
         double phase    = atof( lineTokens[tokenIndex++].c_str() );
         double k        = atof( lineTokens[tokenIndex++].c_str() );
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
      for( unsigned int ii = 0; ii < arraySize; ii++ ){
         int particle1, particle2, particle3, particle4, periodicity;
         double phase, k;
         bondForce->getTorsionParameters( ii, particle1, particle2, particle3, particle4, periodicity, phase, k );
         (void) fprintf( log, "%8d %8d %8d %8d %8d %8d %14.7e %14.7e\n", ii, particle1, particle2, particle3, particle4, periodicity, phase, k );
         if( ii == maxPrint ){
            ii = arraySize - maxPrint;
            if( ii < maxPrint )ii = maxPrint;
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

static int readRBTorsionForce( FILE* filePtr, MapStringInt& forceMap, const StringVector& tokens, System& system, int* lineCount, FILE* log ){

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
   MapStringIntI forceActive    = forceMap.find( RB_TORSION_FORCE );
   if( forceActive != forceMap.end() && (*forceActive).second ){
      system.addForce( bondForce );
      if( log ){
         (void) fprintf( log, "RB torsion force is being included.\n" );
      }
   } else if( log ){
      (void) fprintf( log, "RB torsion force is not being included.\n" );
   }

   int numberOfTorsions            = atoi( tokens[1].c_str() );
   if( log ){
      (void) fprintf( log, "%s number of RBTorsionForce terms=%d\n", methodName.c_str(), numberOfTorsions );
      (void) fflush( log );
   }

#if 0
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
#endif


   for( int ii = 0; ii < numberOfTorsions; ii++ ){
      StringVector lineTokens;
      char* isNotEof = readLine( filePtr, lineTokens, lineCount, log );
      int tokenIndex = 0;
      if( lineTokens.size() > 10 ){
         int index       = atoi( lineTokens[tokenIndex++].c_str() );
         int particle1   = atoi( lineTokens[tokenIndex++].c_str() );
         int particle2   = atoi( lineTokens[tokenIndex++].c_str() );
         int particle3   = atoi( lineTokens[tokenIndex++].c_str() );
         int particle4   = atoi( lineTokens[tokenIndex++].c_str() );
         double c0       = atof( lineTokens[tokenIndex++].c_str() );
         double c1       = atof( lineTokens[tokenIndex++].c_str() );
         double c2       = atof( lineTokens[tokenIndex++].c_str() );
         double c3       = atof( lineTokens[tokenIndex++].c_str() );
         double c4       = atof( lineTokens[tokenIndex++].c_str() );
         double c5       = atof( lineTokens[tokenIndex++].c_str() );


#if 0
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
      (void) fprintf( log, "%s: sample of RBTorsionForce parameters\n", methodName.c_str() );
      for( unsigned int ii = 0; ii < arraySize; ii++ ){
         int particle1, particle2, particle3, particle4;
         double c0, c1, c2, c3, c4, c5;
         bondForce->getTorsionParameters( ii, particle1, particle2, particle3, particle4, c0, c1, c2, c3, c4, c5 );
         (void) fprintf( log, "%8d %8d %8d %8d %8d %14.7e %14.7e %14.7e %14.7e %14.7e %14.7e\n",
                         ii, particle1, particle2, particle3, particle4, c0, c1, c2, c3, c4, c5 );
         if( ii == maxPrint ){
            ii = arraySize - maxPrint;
            if( ii < maxPrint )ii = maxPrint;
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
      int tokenIndex = 0;
      if( lineTokens.size() > 5 ){
         int index                   = atoi( lineTokens[tokenIndex++].c_str() );
         int particle1               = atoi( lineTokens[tokenIndex++].c_str() );
         int particle2               = atoi( lineTokens[tokenIndex++].c_str() );
         double charge               = includeNonbondedExceptions ? atof( lineTokens[tokenIndex++].c_str() ) : 0.0;
         double sigma                = includeNonbondedExceptions ? atof( lineTokens[tokenIndex++].c_str() ) : 1.0;
         double epsilon              = includeNonbondedExceptions ? atof( lineTokens[tokenIndex++].c_str() ) : 0.0;

#if 0
   if( log && ii < 2 )
      (void) fprintf( log, "************************ Setting q to zero ************************\n" );
   charge  = 0.0;
//   sigma   = 1.0;
//   epsilon = 0.0;
#endif
         nonbondedForce.addException( particle1, particle2, charge, sigma, epsilon );
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
      for( unsigned int ii = 0; ii < arraySize; ii++ ){
         int particle1, particle2;
         double chargeProd, sigma, epsilon;
         nonbondedForce.getExceptionParameters( ii, particle1, particle2, chargeProd, sigma, epsilon );
         (void) fprintf( log, "%8d %8d %8d %14.7e %14.7e %14.7e\n", ii, particle1, particle2, chargeProd, sigma, epsilon );
         if( ii == maxPrint ){
            ii = arraySize - maxPrint;
            if( ii < maxPrint )ii = maxPrint;
         }
      }
   }

   return nonbondedForce.getNumExceptions();
}

/**---------------------------------------------------------------------------------------

   Read NonbondedSoftcoreExceptions parameters

   @param filePtr                     file pointer to parameter file
   @param includeNonbondedSoftcoreExceptions 
                                      if set, then include exceptions; otherwise set charge and epsilon to zero and
                                      sigma to 1
   @param tokens                      array of strings from first line of parameter file for this block of parameters
   @param nonbondedForce              NonBondedForce reference
   @param lineCount                   used to track line entries read from parameter file
   @param log                         log file pointer -- may be NULL

   @return number of parameters read

   --------------------------------------------------------------------------------------- */

#ifdef INCLUDE_FREE_ENERGY_PLUGIN
static int readNonbondedSoftcoreExceptions( FILE* filePtr, int includeNonbondedSoftcoreExceptions,
                                            const StringVector& tokens, NonbondedSoftcoreForce& nonbondedForce,
                                            int* lineCount, FILE* log ){

// ---------------------------------------------------------------------------------------

   static const std::string methodName      = "readNonbondedSoftcoreExceptions";
   
// ---------------------------------------------------------------------------------------

   if( tokens.size() < 1 ){
      char buffer[1024];
      (void) sprintf( buffer, "%s no Nonbonded softcore exceptions ???\n", methodName.c_str() );
      throwException(__FILE__, __LINE__, buffer );
      exit(-1);
   }

   int numberOfExceptions           = atoi( tokens[1].c_str() );
   if( log ){
      (void) fprintf( log, "%s number of NonbondedSoftcoreExceptions terms=%d\n", methodName.c_str(), numberOfExceptions );
      (void) fflush( log );
   }
   for( int ii = 0; ii < numberOfExceptions; ii++ ){
      StringVector lineTokens;
      char* isNotEof = readLine( filePtr, lineTokens, lineCount, log );
      int tokenIndex = 0;
      if( lineTokens.size() > 5 ){
         int index                   = atoi( lineTokens[tokenIndex++].c_str() );
         int particle1               = atoi( lineTokens[tokenIndex++].c_str() );
         int particle2               = atoi( lineTokens[tokenIndex++].c_str() );
         double charge               = includeNonbondedSoftcoreExceptions ? atof( lineTokens[tokenIndex++].c_str() ) : 0.0;
         double sigma                = includeNonbondedSoftcoreExceptions ? atof( lineTokens[tokenIndex++].c_str() ) : 1.0;
         double epsilon              = includeNonbondedSoftcoreExceptions ? atof( lineTokens[tokenIndex++].c_str() ) : 0.0;
         double softcoreLJLambda     = 1.0;
         if( includeNonbondedSoftcoreExceptions && lineTokens.size() > tokenIndex ){
            softcoreLJLambda     = atof( lineTokens[tokenIndex++].c_str() );
         }
         nonbondedForce.addException( particle1, particle2, charge, sigma, epsilon, softcoreLJLambda );

      } else if( log ){
         char buffer[1024];
         (void) sprintf( buffer, "%s readNonbondedSoftcoreExceptions tokens incomplete at line=%d\n", methodName.c_str(), *lineCount );
         throwException(__FILE__, __LINE__, buffer );
         exit(-1);
      }
   }

   // diagnostics

   if( log ){
      static const unsigned int maxPrint   = MAX_PRINT;
      unsigned int arraySize               = static_cast<unsigned int>(nonbondedForce.getNumExceptions());
      (void) fprintf( log, "%s: sample of NonbondedSoftcoreExceptions parameters\n", methodName.c_str() );
      for( unsigned int ii = 0; ii < arraySize; ii++ ){
         int particle1, particle2;
         double chargeProd, sigma, epsilon, softcoreLJLambda;
         nonbondedForce.getExceptionParameters( ii, particle1, particle2, chargeProd, sigma, epsilon, softcoreLJLambda );
         (void) fprintf( log, "%8d %8d %8d %14.7e %14.7e %14.7e %14.7e\n", ii, particle1, particle2, chargeProd, sigma, epsilon, softcoreLJLambda );
         if( ii == maxPrint ){
            ii = arraySize - maxPrint;
            if( ii < maxPrint )ii = maxPrint;
         }
      }
   }

   return nonbondedForce.getNumExceptions();
}
#endif

/**---------------------------------------------------------------------------------------

   Set NonbondedForce method

   @param nonbondedForce       force method is to be set for (optional)
   @param nonbondedForceMethod nonbonded force method name
   @param log                  log file pointer -- may be NULL

   @return NonbondedForce enum

   --------------------------------------------------------------------------------------- */

static int setNonbondedForceMethod( NonbondedForce* nonbondedForce, std::string nonbondedForceMethod, FILE* log ){

// ---------------------------------------------------------------------------------------

   static const std::string methodName      = "setNonbondedForceMethod";
   
// ---------------------------------------------------------------------------------------

   if( log ){
      (void) fprintf( log, "%s Nonbonded force is being set %s.\n", methodName.c_str(), nonbondedForceMethod.c_str() );
   }

   NonbondedForce::NonbondedMethod method = NonbondedForce::NoCutoff;
   if( nonbondedForceMethod.compare( "NoCutoff" ) == 0 ){
      method = NonbondedForce::NoCutoff; 
   } else if( nonbondedForceMethod.compare( "CutoffNonPeriodic" ) == 0 ){
      method = NonbondedForce::CutoffNonPeriodic; 
   } else if( nonbondedForceMethod.compare( "CutoffPeriodic" ) == 0 ){
      method = NonbondedForce::CutoffPeriodic; 
   } else if( nonbondedForceMethod.compare( "Ewald" ) == 0 ){
      method = NonbondedForce::Ewald; 
   } else if( nonbondedForceMethod.compare( "PME" ) == 0 ){
      method = NonbondedForce::PME; 
   } else {
      char buffer[1024];
      (void) sprintf( buffer, "nonbondedForce NonbondedForceMethod <%s> is not recognized.\n", nonbondedForceMethod.c_str() );
      if( log ){
            (void) fprintf( log, "%s", buffer ); (void) fflush( log );
      }
      throwException(__FILE__, __LINE__, buffer );
      exit(-1);
   }

   if( nonbondedForce ){
      if( log ){
         (void) fprintf( log, "%s Nonbonded force is being set %s %d.\n", methodName.c_str(), nonbondedForceMethod.c_str(), method );
      }
      nonbondedForce->setNonbondedMethod( method ); 
   }

   return method;

}

/**---------------------------------------------------------------------------------------

   Set NonbondedSoftcoreForce method

   @param nonbondedSoftcoreForce       force method is to be set for (optional)
   @param nonbondedSoftcoreForceMethod nonbonded force method name
   @param log                  log file pointer -- may be NULL

   @return NonbondedSoftcoreForce enum

   --------------------------------------------------------------------------------------- */

#ifdef INCLUDE_FREE_ENERGY_PLUGIN
static int setNonbondedSoftcoreForceMethod( NonbondedSoftcoreForce* nonbondedSoftcoreForce, std::string nonbondedSoftcoreForceMethod, FILE* log ){

// ---------------------------------------------------------------------------------------

   static const std::string methodName      = "setNonbondedSoftcoreForceMethod";
   
// ---------------------------------------------------------------------------------------

   if( log ){
      (void) fprintf( log, "%s Nonbonded softcore force is being set %s.\n", methodName.c_str(), nonbondedSoftcoreForceMethod.c_str() );
   }

   NonbondedSoftcoreForce::NonbondedSoftcoreMethod method = NonbondedSoftcoreForce::NoCutoff;
   if( nonbondedSoftcoreForceMethod.compare( "NoCutoff" ) == 0 ){
      method = NonbondedSoftcoreForce::NoCutoff; 
   } else if( nonbondedSoftcoreForceMethod.compare( "CutoffNonPeriodic" ) == 0 ){
      method = NonbondedSoftcoreForce::CutoffNonPeriodic; 
   } else if( nonbondedSoftcoreForceMethod.compare( "CutoffPeriodic" ) == 0 ){
      method = NonbondedSoftcoreForce::CutoffPeriodic; 
   } else if( nonbondedSoftcoreForceMethod.compare( "Ewald" ) == 0 ){
      method = NonbondedSoftcoreForce::Ewald; 
   } else if( nonbondedSoftcoreForceMethod.compare( "PME" ) == 0 ){
      method = NonbondedSoftcoreForce::PME; 
   } else {
      char buffer[1024];
      (void) sprintf( buffer, "nonbondedSoftcoreForce NonbondedSoftcoreForceMethod <%s> is not recognized.\n", nonbondedSoftcoreForceMethod.c_str() );
      if( log ){
            (void) fprintf( log, "%s", buffer ); (void) fflush( log );
      }
      throwException(__FILE__, __LINE__, buffer );
      exit(-1);
   }

   if( nonbondedSoftcoreForce ){
      if( log ){
         (void) fprintf( log, "%s Nonbonded softcore force is being set %s %d.\n", methodName.c_str(), nonbondedSoftcoreForceMethod.c_str(), method );
      }
      nonbondedSoftcoreForce->setNonbondedMethod( method ); 
   }

   return method;

}
#endif

/**---------------------------------------------------------------------------------------

   Read NonbondedForce parameters

   @param filePtr              file pointer to parameter file
   @param tokens               array of strings from first line of parameter file for this block of parameters
   @param system               System reference
   @param lineCount            used to track line entries read from parameter file
   @param log                  log file pointer -- may be NULL

   @return number of parameters read

   --------------------------------------------------------------------------------------- */

static int readNonbondedForce( FILE* filePtr, MapStringInt& forceMap, const StringVector& tokens,
                               System& system, int* lineCount, MapStringString& inputArgumentMap, FILE* log ){

// ---------------------------------------------------------------------------------------

   static const std::string methodName      = "readNonbondedForce";
   
// ---------------------------------------------------------------------------------------

   if( tokens.size() < 1 ){
      char buffer[1024];
      (void) sprintf( buffer, "%s no Nonbonded bonds number entry???\n", methodName.c_str() );
      throwException(__FILE__, __LINE__, buffer );
      exit(-1);
   }

   NonbondedForce* nonbondedForce         = new NonbondedForce();

   MapStringIntI forceActive              = forceMap.find( NB_FORCE );
   MapStringIntI forceExceptionsActive    = forceMap.find( NB_EXCEPTION_FORCE );
   int includeNonbonded                   = ( forceActive != forceMap.end() && (*forceActive).second ) ? 1 : 0;
   int includeNonbondedExceptions         = ( forceExceptionsActive != forceMap.end() && (*forceExceptionsActive).second ) ? 1 : 0;
   if( includeNonbonded || includeNonbondedExceptions ){
      system.addForce( nonbondedForce );
      if( log ){
         if( includeNonbonded ){
            (void) fprintf( log, "nonbonded force is being included.\n" );
         } else {
            (void) fprintf( log, "nonbonded force is not being included.\n" );
         }
         if( includeNonbondedExceptions ){
            (void) fprintf( log, "nonbonded exceptions are being included.\n" );
         } else {
            (void) fprintf( log, "nonbonded exceptions are not being included.\n" );
         }
      }
   } else if( log ){
      (void) fprintf( log, "nonbonded force is not being included.\n" );
   }

   int numberOfParticles           = atoi( tokens[1].c_str() );
   if( log ){
      (void) fprintf( log, "%s number of NonbondedForce terms=%d\n", methodName.c_str(), numberOfParticles );
      (void) fflush( log );
   }

   // get charge, sigma, epsilon for each particle

   for( int ii = 0; ii < numberOfParticles; ii++ ){
      StringVector lineTokens;
      char* isNotEof = readLine( filePtr, lineTokens, lineCount, log );
      int tokenIndex = 0;
      if( lineTokens.size() > 3 ){
         int index                   = atoi( lineTokens[tokenIndex++].c_str() );
         double charge               = includeNonbonded ? atof( lineTokens[tokenIndex++].c_str() ) : 0.0;
         double sigma                = includeNonbonded ? atof( lineTokens[tokenIndex++].c_str() ) : 1.0;
         double epsilon              = includeNonbonded ? atof( lineTokens[tokenIndex++].c_str() ) : 0.0;
         nonbondedForce->addParticle( charge, sigma, epsilon );
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
   while( hits < 5 ){
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
         } else if( field.compare( "NonbondedForceMethod" ) == 0 ){
            setNonbondedForceMethod( nonbondedForce, tokens[1], log );
            hits++;
         }
      }
   }

   // overrides

   double cutoffDistance = nonbondedForce->getCutoffDistance( );
   if( setDoubleFromMap( inputArgumentMap, "nonbondedCutoffDistance", cutoffDistance ) ){
      nonbondedForce->setCutoffDistance( cutoffDistance );
   }
   
   double ewaldTolerance = nonbondedForce->getEwaldErrorTolerance( );
   if( setDoubleFromMap( inputArgumentMap, "nonbondedEwaldTolerance", ewaldTolerance ) ){
      nonbondedForce->setEwaldErrorTolerance( ewaldTolerance );
   }
   
   double rFDielectric = nonbondedForce->getReactionFieldDielectric( );
   if( setDoubleFromMap( inputArgumentMap, "nonbondedRFDielectric", rFDielectric ) ){
      nonbondedForce->setReactionFieldDielectric( rFDielectric );
   }
   
   std::string nonbondedMethod;
   if( setStringFromMap( inputArgumentMap, "nonbondedForceMethod", nonbondedMethod) ){
      setNonbondedForceMethod( nonbondedForce, nonbondedMethod, log );
   }

   // diagnostics

   if( log ){
      static const unsigned int maxPrint   = MAX_PRINT;
      unsigned int arraySize               = static_cast<unsigned int>(nonbondedForce->getNumParticles());

      (void) fprintf( log, "%s: nonbonded parameters\n", methodName.c_str() );

      // cutoff distance and box

      (void) fprintf( log, "CutoffDistance %14.7e\n", nonbondedForce->getCutoffDistance() );
  
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
          case NonbondedForce::PME:
              nonbondedForceMethod = "PME";
              break;
          default:
              nonbondedForceMethod = "Unknown";
      }
      (void) fprintf( log, "NonbondedForceMethod=%s\n", nonbondedForceMethod.c_str() );
  
      (void) fprintf( log, "charge, sigma, epsilon\n" );
      for( unsigned int ii = 0; ii < arraySize; ii++ ){
         double charge, sigma, epsilon;
         nonbondedForce->getParticleParameters( ii, charge, sigma, epsilon );
         (void) fprintf( log, "%8d %14.7e %14.7e %14.7e\n", ii, charge, sigma, epsilon );
         if( ii == maxPrint ){
            ii = arraySize - maxPrint;
            if( ii < maxPrint )ii = maxPrint;
         }
      }
   }

   return nonbondedForce->getNumParticles();
}

/**---------------------------------------------------------------------------------------

   Read NonbondedSoftcoreForce parameters

   @param filePtr              file pointer to parameter file
   @param tokens               array of strings from first line of parameter file for this block of parameters
   @param system               System reference
   @param lineCount            used to track line entries read from parameter file
   @param log                  log file pointer -- may be NULL

   @return number of parameters read

   --------------------------------------------------------------------------------------- */

#ifdef INCLUDE_FREE_ENERGY_PLUGIN
static int readNonbondedSoftcoreForce( FILE* filePtr, MapStringInt& forceMap, const StringVector& tokens,
                                       System& system, int* lineCount, MapStringString& inputArgumentMap, FILE* log ){

// ---------------------------------------------------------------------------------------

   static const std::string methodName      = "readNonbondedSoftcoreForce";
   
// ---------------------------------------------------------------------------------------

   if( tokens.size() < 1 ){
      char buffer[1024];
      (void) sprintf( buffer, "%s no Nonbonded softcore entries???\n", methodName.c_str() );
      throwException(__FILE__, __LINE__, buffer );
      exit(-1);
   }

   NonbondedSoftcoreForce* nonbondedSoftcore  = new NonbondedSoftcoreForce();

   MapStringIntI forceActive              = forceMap.find( NB_SOFTCORE_FORCE );
   MapStringIntI forceExceptionsActive    = forceMap.find( NB_EXCEPTION_SOFTCORE_FORCE );
   int includeNonbonded                   = ( forceActive != forceMap.end() && (*forceActive).second ) ? 1 : 0;
   int includeNonbondedExceptions         = ( forceExceptionsActive != forceMap.end() && (*forceExceptionsActive).second ) ? 1 : 0;
   if( includeNonbonded || includeNonbondedExceptions ){
      system.addForce( nonbondedSoftcore );
      if( log ){
         if( includeNonbonded ){
            (void) fprintf( log, "nonbonded softcore force is being included.\n" );
         } else {
            (void) fprintf( log, "nonbonded softcore force is not being included.\n" );
         }
         if( includeNonbondedExceptions ){
            (void) fprintf( log, "nonbonded softcore exceptions are being included.\n" );
         } else {
            (void) fprintf( log, "nonbonded softcore exceptions are not being included.\n" );
         }
      }
   } else if( log ){
      (void) fprintf( log, "nonbonded softcore force is not being included.\n" );
   }

   int numberOfParticles           = atoi( tokens[1].c_str() );
   if( log ){
      (void) fprintf( log, "%s number of NonbondedSoftcoreForce terms=%d\n", methodName.c_str(), numberOfParticles );
      (void) fflush( log );
   }

   // get charge, sigma, epsilon for each particle

   for( int ii = 0; ii < numberOfParticles; ii++ ){
      StringVector lineTokens;
      char* isNotEof = readLine( filePtr, lineTokens, lineCount, log );
      int tokenIndex = 0;
      if( lineTokens.size() > 3 ){
         int index                   = atoi( lineTokens[tokenIndex++].c_str() );
         double charge               = includeNonbonded ? atof( lineTokens[tokenIndex++].c_str() ) : 0.0;
         double sigma                = includeNonbonded ? atof( lineTokens[tokenIndex++].c_str() ) : 1.0;
         double epsilon              = includeNonbonded ? atof( lineTokens[tokenIndex++].c_str() ) : 0.0;
        
         double softcoreLJLambda     = 1.0;
         if( includeNonbonded && lineTokens.size() > tokenIndex ){
            softcoreLJLambda = atof( lineTokens[tokenIndex++].c_str() );
         }
         nonbondedSoftcore->addParticle( charge, sigma, epsilon, softcoreLJLambda );
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
   while( hits < 5 ){
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
            nonbondedSoftcore->setCutoffDistance( atof( tokens[1].c_str() ) );
            hits++;
         } else if( field.compare( "RFDielectric" ) == 0 ){
            nonbondedSoftcore->setReactionFieldDielectric( atof( tokens[1].c_str() ) );
            hits++;
         } else if( field.compare( "EwaldRTolerance" ) == 0 ){
            nonbondedSoftcore->setEwaldErrorTolerance( atof( tokens[1].c_str() ) );
            hits++;
         } else if( field.compare( "NonbondedSoftcoreForceExceptions" ) == 0 ){
            readNonbondedSoftcoreExceptions( filePtr, includeNonbondedExceptions, tokens, *nonbondedSoftcore, lineCount, log );
            hits++;
         } else if( field.compare( "NonbondedSoftcoreForceMethod" ) == 0 ){
            setNonbondedSoftcoreForceMethod( nonbondedSoftcore, tokens[1], log );
            hits++;
         }
      }
   }

   // overrides

   double cutoffDistance = nonbondedSoftcore->getCutoffDistance( );
   if( setDoubleFromMap( inputArgumentMap, "nonbondedCutoffDistance", cutoffDistance ) ){
      nonbondedSoftcore->setCutoffDistance( cutoffDistance );
   }
   
   double ewaldTolerance = nonbondedSoftcore->getEwaldErrorTolerance( );
   if( setDoubleFromMap( inputArgumentMap, "nonbondedEwaldTolerance", ewaldTolerance ) ){
      nonbondedSoftcore->setEwaldErrorTolerance( ewaldTolerance );
   }
   
   double rFDielectric = nonbondedSoftcore->getReactionFieldDielectric( );
   if( setDoubleFromMap( inputArgumentMap, "nonbondedRFDielectric", rFDielectric ) ){
      nonbondedSoftcore->setReactionFieldDielectric( rFDielectric );
   }
   
   std::string nonbondedMethod;
   if( setStringFromMap( inputArgumentMap, "nonbondedSoftcoreMethod", nonbondedMethod) ){
      setNonbondedSoftcoreForceMethod( nonbondedSoftcore, nonbondedMethod, log );
   }

   // diagnostics

   if( log ){
      static const unsigned int maxPrint   = MAX_PRINT;
      unsigned int arraySize               = static_cast<unsigned int>(nonbondedSoftcore->getNumParticles());

      (void) fprintf( log, "%s: nonbonded parameters\n", methodName.c_str() );

      // cutoff distance and box

      (void) fprintf( log, "CutoffDistance %14.7e\n", nonbondedSoftcore->getCutoffDistance() );
  
      // nonbond method

      std::string nonbondedSoftcoreMethod;
      switch( nonbondedSoftcore->getNonbondedMethod() ){
          case NonbondedForce::NoCutoff:
              nonbondedSoftcoreMethod = "NoCutoff";
              break;
          case NonbondedForce::CutoffNonPeriodic:
              nonbondedSoftcoreMethod = "CutoffNonPeriodic";
              break;
          case NonbondedForce::CutoffPeriodic:
              nonbondedSoftcoreMethod = "CutoffPeriodic";
              break;
          case NonbondedForce::Ewald:
              nonbondedSoftcoreMethod = "Ewald";
              break;
          case NonbondedForce::PME:
              nonbondedSoftcoreMethod = "PME";
              break;
          default:
              nonbondedSoftcoreMethod = "Unknown";
      }
      (void) fprintf( log, "NonbondedForceMethod=%s\n", nonbondedSoftcoreMethod.c_str() );
  
      (void) fprintf( log, "charge, sigma, epsilon\n" );
      for( unsigned int ii = 0; ii < arraySize; ii++ ){
         double charge, sigma, epsilon, softcoreLJLambda;
         nonbondedSoftcore->getParticleParameters( ii, charge, sigma, epsilon, softcoreLJLambda );
         (void) fprintf( log, "%8d %14.7e %14.7e %14.7e %14.7e\n", ii, charge, sigma, epsilon, softcoreLJLambda );
         if( ii == maxPrint ){
            ii = arraySize - maxPrint;
            if( ii < maxPrint )ii = maxPrint;
         }
      }
   }

   return nonbondedSoftcore->getNumParticles();
}
#endif

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

static int readGBSAOBCForce( FILE* filePtr, MapStringInt& forceMap, const StringVector& tokens, System& system, int* lineCount, FILE* log ){

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
   MapStringIntI forceActive    = forceMap.find( GBSA_OBC_FORCE );
   if( forceActive != forceMap.end() && (*forceActive).second ){
      system.addForce( gbsaObcForce );
      if( log ){
         (void) fprintf( log, "GBSA OBC force is being included.\n" );
      }
   } else if( log ){
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
      int tokenIndex = 0;
      if( lineTokens.size() > 3 ){
         int index            = atoi( lineTokens[tokenIndex++].c_str() );
         double charge        = atof( lineTokens[tokenIndex++].c_str() );
         double radius        = atof( lineTokens[tokenIndex++].c_str() );
         double scalingFactor = atof( lineTokens[tokenIndex++].c_str() );
         gbsaObcForce->addParticle( charge, radius, scalingFactor );
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

      for( unsigned int ii = 0; ii < arraySize; ii++ ){
         double charge, radius, scalingFactor;
         gbsaObcForce->getParticleParameters( ii, charge, radius, scalingFactor );
         (void) fprintf( log, "%8d  %14.7e %14.7e %14.7e\n", ii, charge, radius, scalingFactor );
         if( ii == maxPrint ){
            ii = arraySize - maxPrint;
            if( ii < maxPrint )ii = maxPrint;
         }
      }
   }

   return gbsaObcForce->getNumParticles();
}

/**---------------------------------------------------------------------------------------

   Read GBSAOBCSoftcoreForce parameters

   @param filePtr              file pointer to parameter file
   @param forceFlag            flag signalling whether force is to be added to system
                               if force == 0 || forceFlag & GBSA_OBCSoftcore_FORCE, then included
   @param tokens               array of strings from first line of parameter file for this block of parameters
   @param system               System reference
   @param lineCount            used to track line entries read from parameter file
   @param log                  log file pointer -- may be NULL

   @return number of parameters read

   --------------------------------------------------------------------------------------- */

#ifdef INCLUDE_FREE_ENERGY_PLUGIN
static int readGBSAOBCSoftcoreForce( FILE* filePtr, MapStringInt& forceMap, const StringVector& tokens, System& system, int* lineCount, FILE* log ){

// ---------------------------------------------------------------------------------------

   static const std::string methodName      = "readGBSAOBCSoftcoreForce";
   
// ---------------------------------------------------------------------------------------

   if( tokens.size() < 1 ){
      char buffer[1024];
      (void) sprintf( buffer, "%s no GBSAOBCSoftcore terms entry???\n", methodName.c_str() );
      throwException(__FILE__, __LINE__, buffer );
      exit(-1);
   }

   GBSAOBCSoftcoreForce* gbsaObcSoftcoreForce   = new GBSAOBCSoftcoreForce();
   MapStringIntI forceActive            = forceMap.find( GBSA_OBC_SOFTCORE_FORCE );
   if( forceActive != forceMap.end() && (*forceActive).second ){
      system.addForce( gbsaObcSoftcoreForce );
      if( log ){
         (void) fprintf( log, "GBSA OBCSoftcore force is being included.\n" );
      }
   } else if( log ){
      (void) fprintf( log, "GBSA OBCSoftcore force is not being included.\n" );
   }

   int numberOfParticles           = atoi( tokens[1].c_str() );
   if( log ){
      (void) fprintf( log, "%s number of GBSAOBCSoftcoreForce terms=%d\n", methodName.c_str(), numberOfParticles );
      (void) fflush( log );
   }
   for( int ii = 0; ii < numberOfParticles; ii++ ){
      StringVector lineTokens;
      char* isNotEof = readLine( filePtr, lineTokens, lineCount, log );
      int tokenIndex = 0;
      if( lineTokens.size() > 3 ){
         int index            = atoi( lineTokens[tokenIndex++].c_str() );
         double charge        = atof( lineTokens[tokenIndex++].c_str() );
         double radius        = atof( lineTokens[tokenIndex++].c_str() );
         double scalingFactor = atof( lineTokens[tokenIndex++].c_str() );
         double saScale       = atof( lineTokens[tokenIndex++].c_str() );
         gbsaObcSoftcoreForce->addParticle( charge, radius, scalingFactor, saScale );
      } else {
         char buffer[1024];
         (void) sprintf( buffer, "%s GBSAOBCSoftcoreForce tokens incomplete at line=%d\n", methodName.c_str(), *lineCount );
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
            gbsaObcSoftcoreForce->setSoluteDielectric( atof( tokens[1].c_str() ) );
            hits++;
         } else if( field.compare( "SolventDielectric" ) == 0 ){
            gbsaObcSoftcoreForce->setSolventDielectric( atof( tokens[1].c_str() ) );
            hits++;
         } else {
               char buffer[1024];
               (void) sprintf( buffer, "%s read past GBSA Obc softcore block at line=%d\n", methodName.c_str(), lineCount );
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
      unsigned int arraySize               = static_cast<unsigned int>(gbsaObcSoftcoreForce->getNumParticles());
      (void) fprintf( log, "%s: sample of GBSA OBCSoftcore Force parameters; no. of particles=%d\n",
                      methodName.c_str(), gbsaObcSoftcoreForce->getNumParticles() );
      (void) fprintf( log, "solute/solvent dielectrics: [%10.4f %10.4f]\n",
                      gbsaObcSoftcoreForce->getSoluteDielectric(),  gbsaObcSoftcoreForce->getSolventDielectric() );

      for( unsigned int ii = 0; ii < arraySize; ii++ ){
         double charge, radius, scalingFactor, nonpolarScaleFactor;
         gbsaObcSoftcoreForce->getParticleParameters( ii, charge, radius, scalingFactor, nonpolarScaleFactor );
         (void) fprintf( log, "%8d  %14.7e %14.7e %14.7e %14.7e\n", ii, charge, radius, scalingFactor, nonpolarScaleFactor );
         if( ii == maxPrint ){
            ii = arraySize - maxPrint;
            if( ii < maxPrint )ii = maxPrint;
         }
      }
   }

   return gbsaObcSoftcoreForce->getNumParticles();
}
#endif

/**---------------------------------------------------------------------------------------

   Read GBVIForce parameters

   @param filePtr              file pointer to parameter file
   @param forceFlag            flag signalling whether force is to be added to system
                               if force == 0 || forceFlag & GBVI_FORCE, then included
   @param tokens               array of strings from first line of parameter file for this block of parameters
   @param system               System reference
   @param lineCount            used to track line entries read from parameter file
   @param log                  log file pointer -- may be NULL

   @return number of parameters read

   --------------------------------------------------------------------------------------- */

static int readGBVIForceMod( FILE* filePtr, MapStringInt& forceMap, const StringVector& tokens, System& system, int* lineCount, FILE* log ){

// ---------------------------------------------------------------------------------------

   static const std::string methodName      = "readGBVIForce";
   
// ---------------------------------------------------------------------------------------

   if( tokens.size() < 1 ){
      char buffer[1024];
      (void) sprintf( buffer, "%s no GBVI terms entry???\n", methodName.c_str() );
      throwException(__FILE__, __LINE__, buffer );
      exit(-1);
   }

//   GBVIForce* gbviForce          = new GBVIForce();
   MapStringIntI forceActive     = forceMap.find( GBVI_FORCE );
   if( forceActive != forceMap.end() && (*forceActive).second ){
      forceMap[GBVI_SOFTCORE_FORCE] = 0;
//      system.addForce( gbviForce );
      if( log ){
         (void) fprintf( log, "GBVI force is being included & GBVI softcore force excluded.\n" );
      }
   } else if( log ){
      (void) fprintf( log, "GBVI force is not being included.\n" );
   }

   int numberOfParticles           = atoi( tokens[1].c_str() );
   if( log ){
      (void) fprintf( log, "%s number of GBVIForce terms=%d\n", methodName.c_str(), numberOfParticles );
      (void) fflush( log );
   }
   for( int ii = 0; ii < numberOfParticles; ii++ ){
      StringVector lineTokens;
      char* isNotEof = readLine( filePtr, lineTokens, lineCount, log );
      int tokenIndex = 0;
      if( lineTokens.size() > 3 ){
         int index            = atoi( lineTokens[tokenIndex++].c_str() );
         double charge        = atof( lineTokens[tokenIndex++].c_str() );
         double radius        = atof( lineTokens[tokenIndex++].c_str() );
         double gamma         = atof( lineTokens[tokenIndex++].c_str() );
//         gbviForce->addParticle( charge, radius, gamma );
      } else {
         char buffer[1024];
         (void) sprintf( buffer, "%s GBVIForce tokens incomplete at line=%d\n", methodName.c_str(), *lineCount );
         throwException(__FILE__, __LINE__, buffer );
         exit(-1);
      }
   }

   char* isNotEof                 = "1";
   int hits                       = 0;
   while( hits < 3 ){
      StringVector tokens;
      isNotEof = readLine( filePtr, tokens, lineCount, log );
      if( isNotEof && tokens.size() > 0 ){

         std::string field       = tokens[0];
         if( field.compare( "SoluteDielectric" ) == 0 ){
//            gbviForce->setSoluteDielectric( atof( tokens[1].c_str() ) );
            hits++;
         } else if( field.compare( "SolventDielectric" ) == 0 ){
//           gbviForce->setSolventDielectric( atof( tokens[1].c_str() ) );
            hits++;
         } else if( field.compare( "GBVIBonds" ) == 0 ){

            int numberOfBonds = atoi( tokens[1].c_str() );

            for( int ii = 0; ii < numberOfBonds; ii++ ){
               StringVector lineTokens;
               char* isNotEof = readLine( filePtr, lineTokens, lineCount, log );
               int tokenIndex = 0;
               if( lineTokens.size() > 3 ){
                  int index            = atoi( lineTokens[tokenIndex++].c_str() );
                  int atomI            = atoi( lineTokens[tokenIndex++].c_str() );
                  int atomJ            = atoi( lineTokens[tokenIndex++].c_str() );
                  double bondLength    = atof( lineTokens[tokenIndex++].c_str() );
//                  gbviForce->addBond( atomI, atomJ, bondLength );
               }
            }
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

#if 0
   if( log ){
      static const unsigned int maxPrint   = MAX_PRINT;
      unsigned int arraySize               = static_cast<unsigned int>(gbviForce->getNumParticles());
      (void) fprintf( log, "%s: sample of GBVI Force parameters; no. of particles=%d\n",
                      methodName.c_str(), gbviForce->getNumParticles() );
      (void) fprintf( log, "solute/solvent dielectrics: [%10.4f %10.4f]\n",
                      gbviForce->getSoluteDielectric(),  gbviForce->getSolventDielectric() );

      for( unsigned int ii = 0; ii < arraySize; ii++ ){
         double charge, radius, gamma;
         gbviForce->getParticleParameters( ii, charge, radius, gamma );
         (void) fprintf( log, "%8d  %14.7e %14.7e %14.7e\n", ii, charge, radius, gamma);
         if( ii == maxPrint ){
            ii = arraySize - maxPrint;
            if( ii < maxPrint )ii = maxPrint;
         }
      }
      arraySize               = static_cast<unsigned int>(gbviForce->getNumBonds());
      (void) fprintf( log, "%s: sample of GBVI: no. of bonds=%d\n",
                      methodName.c_str(), gbviForce->getNumBonds() );
      for( unsigned int ii = 0; ii < arraySize; ii++ ){
         int atomI, atomJ;
         double bondLength;
         gbviForce->getBondParameters( ii, atomI, atomJ, bondLength );
         (void) fprintf( log, "%8d %8d %8d %14.7e\n", ii, atomI, atomJ, bondLength );
         if( ii == maxPrint ){
            ii = arraySize - maxPrint;
            if( ii < maxPrint )ii = maxPrint;
         }
      }
   }

   return gbviForce->getNumParticles();
#endif
   return numberOfParticles;
}

/**---------------------------------------------------------------------------------------

   Read GBVIForce parameters

   @param filePtr              file pointer to parameter file
   @param forceFlag            flag signalling whether force is to be added to system
                               if force == 0 || forceFlag & GBVI_FORCE, then included
   @param tokens               array of strings from first line of parameter file for this block of parameters
   @param system               System reference
   @param lineCount            used to track line entries read from parameter file
   @param log                  log file pointer -- may be NULL

   @return number of parameters read

   --------------------------------------------------------------------------------------- */

static int readGBVIForce( FILE* filePtr, MapStringInt& forceMap, const StringVector& tokens, System& system, int* lineCount, FILE* log ){

// ---------------------------------------------------------------------------------------

   static const std::string methodName      = "readGBVIForce";
   
// ---------------------------------------------------------------------------------------

   if( tokens.size() < 1 ){
      char buffer[1024];
      (void) sprintf( buffer, "%s no GBVI terms entry???\n", methodName.c_str() );
      throwException(__FILE__, __LINE__, buffer );
      exit(-1);
   }

   GBVIForce* gbviForce          = new GBVIForce();
   MapStringIntI forceActive     = forceMap.find( GBVI_FORCE );
   if( forceActive != forceMap.end() && (*forceActive).second ){
      forceMap[GBVI_SOFTCORE_FORCE] = 0;
      system.addForce( gbviForce );
      if( log ){
         (void) fprintf( log, "GBVI force is being included & GBVI softcore force excluded.\n" );
      }
   } else if( log ){
      (void) fprintf( log, "GBVI force is not being included.\n" );
   }

   int numberOfParticles           = atoi( tokens[1].c_str() );
   if( log ){
      (void) fprintf( log, "%s number of GBVIForce terms=%d\n", methodName.c_str(), numberOfParticles );
      (void) fflush( log );
   }
   for( int ii = 0; ii < numberOfParticles; ii++ ){
      StringVector lineTokens;
      char* isNotEof = readLine( filePtr, lineTokens, lineCount, log );
      int tokenIndex = 0;
      if( lineTokens.size() > 3 ){
         int index            = atoi( lineTokens[tokenIndex++].c_str() );
         double charge        = atof( lineTokens[tokenIndex++].c_str() );
         double radius        = atof( lineTokens[tokenIndex++].c_str() );
         double gamma         = atof( lineTokens[tokenIndex++].c_str() );
         gbviForce->addParticle( charge, radius, gamma );
      } else {
         char buffer[1024];
         (void) sprintf( buffer, "%s GBVIForce tokens incomplete at line=%d\n", methodName.c_str(), *lineCount );
         throwException(__FILE__, __LINE__, buffer );
         exit(-1);
      }
   }

   char* isNotEof                 = "1";
   int hits                       = 0;
   while( hits < 3 ){
      StringVector tokens;
      isNotEof = readLine( filePtr, tokens, lineCount, log );
      if( isNotEof && tokens.size() > 0 ){

         std::string field       = tokens[0];
         if( field.compare( "SoluteDielectric" ) == 0 ){
            gbviForce->setSoluteDielectric( atof( tokens[1].c_str() ) );
            hits++;
         } else if( field.compare( "SolventDielectric" ) == 0 ){
            gbviForce->setSolventDielectric( atof( tokens[1].c_str() ) );
            hits++;
         } else if( field.compare( "GBVIBonds" ) == 0 ){

            int numberOfBonds = atoi( tokens[1].c_str() );

            for( int ii = 0; ii < numberOfBonds; ii++ ){
               StringVector lineTokens;
               char* isNotEof = readLine( filePtr, lineTokens, lineCount, log );
               int tokenIndex = 0;
               if( lineTokens.size() > 3 ){
                  int index            = atoi( lineTokens[tokenIndex++].c_str() );
                  int atomI            = atoi( lineTokens[tokenIndex++].c_str() );
                  int atomJ            = atoi( lineTokens[tokenIndex++].c_str() );
                  double bondLength    = atof( lineTokens[tokenIndex++].c_str() );
                  gbviForce->addBond( atomI, atomJ, bondLength );
               }
            }
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
      unsigned int arraySize               = static_cast<unsigned int>(gbviForce->getNumParticles());
      (void) fprintf( log, "%s: sample of GBVI Force parameters; no. of particles=%d\n",
                      methodName.c_str(), gbviForce->getNumParticles() );
      (void) fprintf( log, "solute/solvent dielectrics: [%10.4f %10.4f]\n",
                      gbviForce->getSoluteDielectric(),  gbviForce->getSolventDielectric() );

      for( unsigned int ii = 0; ii < arraySize; ii++ ){
         double charge, radius, gamma;
         gbviForce->getParticleParameters( ii, charge, radius, gamma );
         (void) fprintf( log, "%8d  %14.7e %14.7e %14.7e\n", ii, charge, radius, gamma);
         if( ii == maxPrint ){
            ii = arraySize - maxPrint;
            if( ii < maxPrint )ii = maxPrint;
         }
      }
      arraySize               = static_cast<unsigned int>(gbviForce->getNumBonds());
      (void) fprintf( log, "%s: sample of GBVI: no. of bonds=%d\n",
                      methodName.c_str(), gbviForce->getNumBonds() );
      for( unsigned int ii = 0; ii < arraySize; ii++ ){
         int atomI, atomJ;
         double bondLength;
         gbviForce->getBondParameters( ii, atomI, atomJ, bondLength );
         (void) fprintf( log, "%8d %8d %8d %14.7e\n", ii, atomI, atomJ, bondLength );
         if( ii == maxPrint ){
            ii = arraySize - maxPrint;
            if( ii < maxPrint )ii = maxPrint;
         }
      }
   }

   return gbviForce->getNumParticles();
}

/**---------------------------------------------------------------------------------------

   Read GBVISoftcoreForce parameters

   @param filePtr              file pointer to parameter file
   @param forceFlag            flag signalling whether force is to be added to system
                               if force == 0 || forceFlag & GBVISoftcore_FORCE, then included
   @param tokens               array of strings from first line of parameter file for this block of parameters
   @param system               System reference
   @param lineCount            used to track line entries read from parameter file
   @param log                  log file pointer -- may be NULL

   @return number of parameters read

   --------------------------------------------------------------------------------------- */

#ifdef INCLUDE_FREE_ENERGY_PLUGIN
static int readGBVISoftcoreForce( FILE* filePtr, MapStringInt& forceMap, const StringVector& tokens, System& system, int* lineCount, FILE* log ){

// ---------------------------------------------------------------------------------------

   static const std::string methodName      = "readGBVISoftcoreForce";
   
// ---------------------------------------------------------------------------------------

   if( tokens.size() < 1 ){
      char buffer[1024];
      (void) sprintf( buffer, "%s no GBVISoftcore terms entry???\n", methodName.c_str() );
      throwException(__FILE__, __LINE__, buffer );
      exit(-1);
   }

   GBVISoftcoreForce* gbviForce          = new GBVISoftcoreForce();
   MapStringIntI forceActive             = forceMap.find( GBVI_SOFTCORE_FORCE );
   if( forceActive != forceMap.end() && (*forceActive).second ){
      system.addForce( gbviForce );
      forceMap[GBVI_FORCE] = 0;
      if( log ){
         (void) fprintf( log, "GBVISoftcore force is being included and GBVI excluded.\n" );
      }
   } else if( log ){
      (void) fprintf( log, "GBVISoftcore force is not being included.\n" );
   }

   int numberOfParticles           = atoi( tokens[1].c_str() );
   if( log ){
      (void) fprintf( log, "%s number of GBVISoftcoreForce terms=%d\n", methodName.c_str(), numberOfParticles );
      (void) fflush( log );
   }
   for( int ii = 0; ii < numberOfParticles; ii++ ){
      StringVector lineTokens;
      char* isNotEof = readLine( filePtr, lineTokens, lineCount, log );
      int tokenIndex = 0;
      if( lineTokens.size() > 3 ){
         int index                    = atoi( lineTokens[tokenIndex++].c_str() );
         double charge                = atof( lineTokens[tokenIndex++].c_str() );
         double radius                = atof( lineTokens[tokenIndex++].c_str() );
         double gamma                 = atof( lineTokens[tokenIndex++].c_str() );
         double bornRadiusScaleFactor = 1.0;
         if( lineTokens.size() > tokenIndex ){
            bornRadiusScaleFactor = atof( lineTokens[tokenIndex++].c_str() );
         }
         gbviForce->addParticle( charge, radius, gamma, bornRadiusScaleFactor );
      } else {
         char buffer[1024];
         (void) sprintf( buffer, "%s GBVISoftcoreForce tokens incomplete at line=%d\n", methodName.c_str(), *lineCount );
         throwException(__FILE__, __LINE__, buffer );
         exit(-1);
      }
   }

   char* isNotEof                 = "1";
   int hits                       = 0;
   while( hits < 6 ){
      StringVector tokens;
      isNotEof = readLine( filePtr, tokens, lineCount, log );
      if( isNotEof && tokens.size() > 0 ){

         std::string field       = tokens[0];
         if( field.compare( "SoluteDielectric" ) == 0 ){
            gbviForce->setSoluteDielectric( atof( tokens[1].c_str() ) );
            hits++;
         } else if( field.compare( "SolventDielectric" ) == 0 ){
            gbviForce->setSolventDielectric( atof( tokens[1].c_str() ) );
            hits++;
         } else if( field.compare( "BornRadiusScalingMethod" ) == 0 ){
            int method = atoi( tokens[1].c_str() );
//method = 0;
//(void) fprintf( log, "%s: !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! BornRadiusScalingMethod forced to NoScale!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n", methodName.c_str() );


            if( method == 0 ){
                gbviForce->setBornRadiusScalingMethod( GBVISoftcoreForce::NoScaling );
            } else if( method == 1 ){
                gbviForce->setBornRadiusScalingMethod( GBVISoftcoreForce::Tanh );
            } else if( method == 2 ){
                gbviForce->setBornRadiusScalingMethod( GBVISoftcoreForce::QuinticSpline );
            } else {
               // not recognized force error
               (void) fprintf( log, "%s: BornRadiusScalingMethod id=%s not recognized.\n", methodName.c_str(), tokens[1].c_str() );
               (void) fprintf( stderr, "%s: BornRadiusScalingMethod id=%s not recognized.\n", methodName.c_str(), tokens[1].c_str() );
               (void) fflush( NULL );
               exit(0);
            }
            hits++;
         } else if( field.compare( "QuinticLowerLimitFactor" ) == 0 ){
            gbviForce->setQuinticLowerLimitFactor( atof( tokens[1].c_str() ) );
            hits++;
         } else if( field.compare( "QuinticUpperBornRadiusLimit" ) == 0 ){
            gbviForce->setQuinticUpperBornRadiusLimit( atof( tokens[1].c_str() ) );
            hits++;
         } else if( field.compare( "GBVISoftcoreBonds" ) == 0 ){

            int numberOfBonds = atoi( tokens[1].c_str() );

            for( int ii = 0; ii < numberOfBonds; ii++ ){
               StringVector lineTokens;
               char* isNotEof = readLine( filePtr, lineTokens, lineCount, log );
               int tokenIndex = 0;
               if( lineTokens.size() > 3 ){
                  int index            = atoi( lineTokens[tokenIndex++].c_str() );
                  int atomI            = atoi( lineTokens[tokenIndex++].c_str() );
                  int atomJ            = atoi( lineTokens[tokenIndex++].c_str() );
                  double bondLength    = atof( lineTokens[tokenIndex++].c_str() );
                  gbviForce->addBond( atomI, atomJ, bondLength );
               }
            }
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
      unsigned int arraySize               = static_cast<unsigned int>(gbviForce->getNumParticles());
      (void) fprintf( log, "%s: sample of GBVISoftcore Force parameters; no. of particles=%d\n",
                      methodName.c_str(), gbviForce->getNumParticles() );
      (void) fprintf( log, "solute/solvent dielectrics: [%10.4f %10.4f]\n",
                      gbviForce->getSoluteDielectric(),  gbviForce->getSolventDielectric() );
      (void) fprintf( log, "Born radius scaling method=%d: param[%10.4f %10.4f] [0=none, 1=tanh (not implemented), 2=quintic]\n",
                      gbviForce->getBornRadiusScalingMethod(),
                      gbviForce->getQuinticLowerLimitFactor(), gbviForce->getQuinticUpperBornRadiusLimit() );

      for( unsigned int ii = 0; ii < arraySize; ii++ ){
         double charge, radius, gamma, bornRadiusScaleFactor;
         gbviForce->getParticleParameters( ii, charge, radius, gamma, bornRadiusScaleFactor );
         (void) fprintf( log, "%8d  %14.7e %14.7e %14.7e %14.7e\n", ii, charge, radius, gamma, bornRadiusScaleFactor);
         if( ii == maxPrint ){
            ii = arraySize - maxPrint;
            if( ii < maxPrint )ii = maxPrint;
         }
      }
      arraySize               = static_cast<unsigned int>(gbviForce->getNumBonds());
      (void) fprintf( log, "%s: sample of GBVISoftcore: no. of bonds=%d\n",
                      methodName.c_str(), gbviForce->getNumBonds() );
      for( unsigned int ii = 0; ii < arraySize; ii++ ){
         int atomI, atomJ;
         double bondLength;
         gbviForce->getBondParameters( ii, atomI, atomJ, bondLength );
         (void) fprintf( log, "%8d %8d %8d %14.7e\n", ii, atomI, atomJ, bondLength );
         if( ii == maxPrint ){
            ii = arraySize - maxPrint;
            if( ii < maxPrint )ii = maxPrint;
         }
      }
   }

   return gbviForce->getNumParticles();
}
#endif

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
      int tokenIndex = 0;
      if( lineTokens.size() > 3 ){
         int index            = atoi( lineTokens[tokenIndex++].c_str() );
         int particle1        = atoi( lineTokens[tokenIndex++].c_str() );
         int particle2        = atoi( lineTokens[tokenIndex++].c_str() );
         double distance      = atof( lineTokens[tokenIndex++].c_str() );
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
         if( ii == maxPrint ){
            ii = arraySize - maxPrint;
            if( ii < maxPrint )ii = maxPrint;
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

int readVec3( FILE* filePtr, const StringVector& tokens, std::vector<Vec3>& coordinates, int* lineCount, FILE* log ){

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
      int tokenIndex = 0;
      if( lineTokens.size() > 3 ){
         int index            = atoi( lineTokens[tokenIndex++].c_str() );
         double xCoord        = atof( lineTokens[tokenIndex++].c_str() );
         double yCoord        = atof( lineTokens[tokenIndex++].c_str() );
         double zCoord        = atof( lineTokens[tokenIndex++].c_str() );
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
      for( unsigned int ii = 0; ii < coordinates.size(); ii++ ){
         (void) fprintf( log, "%6u [%14.7e %14.7e %14.7e]\n", ii,
                         coordinates[ii][0], coordinates[ii][1], coordinates[ii][2] );
         if( ii == maxPrint ){
            ii = coordinates.size() - maxPrint;
            if( ii < maxPrint )ii = maxPrint;
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

Integrator* readParameterFile( const std::string& inputParameterFile, MapStringInt& forceMap, System& system,
                               std::vector<Vec3>& coordinates, 
                               std::vector<Vec3>& velocities,
                               std::vector<Vec3>& forces, double* kineticEnergy, double* potentialEnergy,
                               MapStringVectorOfVectors& supplementary,
                               MapStringString& inputArgumentMap, FILE* inputLog ){

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
         } else if( field.compare( "Box" ) == 0 ){

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
            system.setDefaultPeriodicBoxVectors( box[0], box[1], box[2] );
            Vec3 a, b, c;
            system.getDefaultPeriodicBoxVectors( a, b, c);
            if( log ){
               (void) fprintf( log, "Box [%14.7f %14.7f %14.7f]\n    [%14.7f %14.7f %14.7f]\n    [%14.7f %14.7f %14.7f]\n",
                               a[0], a[1], a[2], b[0], b[1], b[2], c[0], c[1], c[2] );
            }

         } else if( field.compare( "CMMotionRemover" ) == 0 ){
            int frequency = atoi( tokens[1].c_str() );
            system.addForce( new CMMotionRemover( frequency ) );
            if( log ){
               (void) fprintf( log, "CMMotionRemover added w/ frequency=%d at line=%d\n", frequency, lineCount );
            }
         } else if( field.compare( "HarmonicBondForce" ) == 0 ){
            readHarmonicBondForce( filePtr, forceMap, tokens, system, &lineCount, log );
         } else if( field.compare( "HarmonicAngleForce" ) == 0 ){
            readHarmonicAngleForce( filePtr, forceMap, tokens, system, &lineCount, log );
         } else if( field.compare( "PeriodicTorsionForce" ) == 0 ){
            readPeriodicTorsionForce( filePtr, forceMap, tokens, system, &lineCount, log );
         } else if( field.compare( "RBTorsionForce" ) == 0 ){
            readRBTorsionForce( filePtr, forceMap, tokens, system, &lineCount, log );
         } else if( field.compare( "NonbondedForce" ) == 0 ){
            readNonbondedForce( filePtr, forceMap, tokens, system, &lineCount, inputArgumentMap, log );
#ifdef INCLUDE_FREE_ENERGY_PLUGIN
         } else if( field.compare( "NonbondedSoftcoreForce" ) == 0 ){
            readNonbondedSoftcoreForce( filePtr, forceMap, tokens, system, &lineCount, inputArgumentMap, log );
#endif
         } else if( field.compare( "GBSAOBCForce" ) == 0 ){
            readGBSAOBCForce( filePtr, forceMap, tokens, system, &lineCount, log );
#ifdef INCLUDE_FREE_ENERGY_PLUGIN
         } else if( field.compare( "GBSAOBCSoftcoreForce" ) == 0 ){
            readGBSAOBCSoftcoreForce( filePtr, forceMap, tokens, system, &lineCount, log );
#endif
         } else if( field.compare( "GBVIForce" ) == 0 ){
            readGBVIForce( filePtr, forceMap, tokens, system, &lineCount, log );
#ifdef INCLUDE_FREE_ENERGY_PLUGIN
         } else if( field.compare( "GBVISoftcoreForce" ) == 0 ){
            readGBVISoftcoreForce( filePtr, forceMap, tokens, system, &lineCount, log );
#endif
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
         } else if( field.compare( "GromacsHarmonicBondForce" )        == 0 ||
                    field.compare( "GromacsHarmonicAngleForce" )       == 0 ||
                    field.compare( "GromacsPeriodicTorsionForce" )     == 0 ||
                    field.compare( "GromacsRBTorsionForce" )           == 0 ||
                    field.compare( "GromacsNonbondedForceExceptions" ) == 0 ||
                    field.compare( "GromacsNonbondedForce" )           == 0 ){

            std::vector< std::vector<double> > vectorOfVectors;
            readVectorOfVectors( filePtr, tokens, vectorOfVectors, &lineCount, field, log );
            if( supplementary.find( field ) == supplementary.end() ){
                supplementary[field] = vectorOfVectors;
            }

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

   // close file

   (void) fclose( filePtr );
 
   if( log ){
      (void) fprintf( log, "Read %d lines from file=<%s>\n", lineCount, inputParameterFile.c_str() );
      (void) fflush( log );
   }

   return returnIntegrator;
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
//        seed                                            = brownianIntegrator->getRandomNumberSeed();
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

static int _synchContexts( const Context& context1, Context& context2 ){

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
 * @param deviceId        deviceId (Cuda only)
 * @param log             log file reference
 *
 * @return OpenMM context
 *
   --------------------------------------------------------------------------------------- */

Context* _getContext( System* system, Context* inputContext, Integrator* inputIntegrator, const std::string& platformName,
                      const std::string& idString, std::string deviceId, FILE* log ){

// ---------------------------------------------------------------------------------------

   static const std::string methodName      = "_getContext";

// ---------------------------------------------------------------------------------------

    // Create a context and initialize it.

    Context* context;
    ReferencePlatform referencePlatform;

    if( platformName.compare( "ReferencePlatform" ) == 0 ){
       context = new Context( *system, *inputIntegrator, Platform::getPlatformByName("Reference") );
    } else {
       context                    = new Context( *system, *inputIntegrator, Platform::getPlatformByName("Cuda") );
    }
    if( log ){
       (void) fprintf( log, "%s Using Platform: %s device=%s\n", idString.c_str(), context->getPlatform().getName().c_str(), deviceId.c_str() );
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

static int getForceStrings( System& system, StringVector& forceStringArray, FILE* log ){

    // print active forces and relevant parameters

    for( int ii = 0; ii < system.getNumForces(); ii++ ) {

        int hit                 = 0;
        Force& force            = system.getForce(ii);

        // bond

        if( !hit ){

            try {
               HarmonicBondForce& harmonicBondForce = dynamic_cast<HarmonicBondForce&>(force);
               forceStringArray.push_back( HARMONIC_BOND_FORCE );
               hit++;
            } catch( std::bad_cast ){
            }
        }
    
        // angle

        if( !hit ){
    
            try {
               HarmonicAngleForce& harmonicAngleForce = dynamic_cast<HarmonicAngleForce&>(force);
               forceStringArray.push_back( HARMONIC_ANGLE_FORCE );
               hit++;
            } catch( std::bad_cast ){
            }
        }

        // PeriodicTorsionForce
    
        if( !hit ){
    
            try {
               PeriodicTorsionForce & periodicTorsionForce = dynamic_cast<PeriodicTorsionForce&>(force);
               forceStringArray.push_back( PERIODIC_TORSION_FORCE );
               hit++;
            } catch( std::bad_cast ){
            }
        }
    
        // RBTorsionForce
    
        if( !hit ){
            try {
               RBTorsionForce& rBTorsionForce = dynamic_cast<RBTorsionForce&>(force);
               forceStringArray.push_back( RB_TORSION_FORCE );
               hit++;
            } catch( std::bad_cast ){
            }
        }
    
        // nonbonded
    
        if( !hit ){
            try {
               NonbondedForce& nbForce = dynamic_cast<NonbondedForce&>(force);
               std::stringstream nonbondedForceMethod;
               hit++;
               switch( nbForce.getNonbondedMethod() ){
                  case NonbondedForce::NoCutoff:
                      nonbondedForceMethod << "NoCutoff";
                      break;
                  case NonbondedForce::CutoffNonPeriodic:
                      nonbondedForceMethod << "CutoffNonPeriodic_Cut=";
                      nonbondedForceMethod << nbForce.getCutoffDistance();
                      break;
                  case NonbondedForce::CutoffPeriodic:
                      nonbondedForceMethod << "CutoffPeriodic_Cut=";
                      nonbondedForceMethod << nbForce.getCutoffDistance();
                      break;
                  case NonbondedForce::Ewald:
                      nonbondedForceMethod << "Ewald_Tol=";
                      nonbondedForceMethod << nbForce.getEwaldErrorTolerance();
                      break;
                  case NonbondedForce::PME:
                      nonbondedForceMethod << "PME";
                      break;
                  default:
                      nonbondedForceMethod << "Unknown";
               }
               forceStringArray.push_back( NB_FORCE + nonbondedForceMethod.str() );
               int nbExceptions = 0;
               for( int ii = 0; ii < nbForce.getNumExceptions() && nbExceptions == 0; ii++ ){
                    int particle1, particle2;
                    double chargeProd, sigma, epsilon;
                    nbForce.getExceptionParameters(ii, particle1, particle2, chargeProd, sigma, epsilon);
                    if( fabs( chargeProd ) > 0.0 || fabs( epsilon ) > 0.0 ){
                        nbExceptions = 1;
                    }
                }
                if( nbExceptions ){
                   forceStringArray.push_back( NB_EXCEPTION_FORCE );
                }
            } catch( std::bad_cast ){
            }
        } 

        // nonbonded softcore
    
#ifdef INCLUDE_FREE_ENERGY_PLUGIN
        if( !hit ){
            try {
               NonbondedSoftcoreForce& nbForce = dynamic_cast<NonbondedSoftcoreForce&>(force);
               std::stringstream nonbondedForceMethod;
               hit++;
               switch( nbForce.getNonbondedMethod() ){
                  case NonbondedSoftcoreForce::NoCutoff:
                      nonbondedForceMethod << "NoCutoff";
                      break;
                  case NonbondedForce::CutoffNonPeriodic:
                      nonbondedForceMethod << "CutoffNonPeriodic_Cut=";
                      nonbondedForceMethod << nbForce.getCutoffDistance();
                      break;
                  case NonbondedForce::CutoffPeriodic:
                      nonbondedForceMethod << "CutoffPeriodic_Cut=";
                      nonbondedForceMethod << nbForce.getCutoffDistance();
                      break;
                  case NonbondedForce::Ewald:
                      nonbondedForceMethod << "Ewald_Tol=";
                      nonbondedForceMethod << nbForce.getEwaldErrorTolerance();
                      break;
                  case NonbondedForce::PME:
                      nonbondedForceMethod << "PME";
                      break;
                  default:
                      nonbondedForceMethod << "Unknown";
               }
               forceStringArray.push_back( NB_SOFTCORE_FORCE + nonbondedForceMethod.str() );
               int nbExceptions = 0;
               for( int ii = 0; ii < nbForce.getNumExceptions() && nbExceptions == 0; ii++ ){
                    int particle1, particle2;
                    double chargeProd, sigma, epsilon;
                    nbForce.getExceptionParameters(ii, particle1, particle2, chargeProd, sigma, epsilon);
                    if( fabs( chargeProd ) > 0.0 || fabs( epsilon ) > 0.0 ){
                        nbExceptions = 1;
                    }
                }
                if( nbExceptions ){
                   forceStringArray.push_back( NB_EXCEPTION_SOFTCORE_FORCE );
                }
            } catch( std::bad_cast ){
            }
        } 
#endif

        // GBSA OBC
    
        if( !hit ){
            try {
               GBSAOBCForce& obcForce = dynamic_cast<GBSAOBCForce&>(force);
               forceStringArray.push_back( GBSA_OBC_FORCE );
               hit++;
            } catch( std::bad_cast ){
            }
        }
    
        // GBSA OBC softcore
    
#ifdef INCLUDE_FREE_ENERGY_PLUGIN
        if( !hit ){
            try {
               GBSAOBCSoftcoreForce& obcForce = dynamic_cast<GBSAOBCSoftcoreForce&>(force);
               forceStringArray.push_back( GBSA_OBC_SOFTCORE_FORCE );
               hit++;
            } catch( std::bad_cast ){
            }
        }
#endif
    
        // GBVI
    
        if( !hit ){
            try {
               GBVIForce& gbviForce = dynamic_cast<GBVIForce&>(force);
               forceStringArray.push_back( GBVI_FORCE );
               hit++;
            } catch( std::bad_cast ){
            }
        }
    
        // GBVI softcore
    
#ifdef INCLUDE_FREE_ENERGY_PLUGIN
        if( !hit ){
            try {
               GBVISoftcoreForce& gbviForce = dynamic_cast<GBVISoftcoreForce&>(force);
               forceStringArray.push_back( GBVI_SOFTCORE_FORCE );
               hit++;
            } catch( std::bad_cast ){
            }
        }
#endif
    
        // COM

        if( !hit ){
    
            try {
               CMMotionRemover& cMMotionRemover = dynamic_cast<CMMotionRemover&>(force);
               hit++;
            } catch( std::bad_cast ){
            }
        }

        if( !hit && log ){
           (void) fprintf( log, "   entry=%2d force not recognized XXXX\n", ii );
        }

    }

    return 0;
}

/** 
 * Check that energy and force are consistent
 * 
 * @return DefaultReturnValue or ErrorReturnValue
 *
 */

static int checkEnergyForceConsistent( Context& context, MapStringString& inputArgumentMap,
                                       FILE* log, FILE* summaryFile ) {
    
// ---------------------------------------------------------------------------------------

   int applyAssertion                     = 1;
   double delta                           = 1.0e-04;
   double tolerance                       = 0.01;
  
   static const std::string methodName    = "checkEnergyForceConsistent";

// ---------------------------------------------------------------------------------------

   setIntFromMap(    inputArgumentMap, "applyAssertion",          applyAssertion  );
   setDoubleFromMap( inputArgumentMap, "energyForceDelta",        delta           );
   setDoubleFromMap( inputArgumentMap, "energyForceTolerance",    tolerance       );

   StringVector forceStringArray;
   System system = context.getSystem();
   getForceStrings( system, forceStringArray, log );

   if( log ){
      (void) fprintf( log, "%s delta=%.3e tolerance=%.3e applyAssertion=%d\n", methodName.c_str(), delta, tolerance, applyAssertion );
      (void) fprintf( log, "\nForces:\n" );
      for( StringVectorCI ii = forceStringArray.begin(); ii != forceStringArray.end(); ii++ ){
         (void) fprintf( log, "   %s\n", (*ii).c_str() );
      }
      (void) fflush( log );
   }

   int returnStatus                       = 0;

   // get positions, forces and potential energy

   int types                              = State::Positions | State::Velocities | State::Forces | State::Energy;

   State state                            = context.getState( types );

   std::vector<Vec3> coordinates          = state.getPositions();
   std::vector<Vec3> velocities           = state.getVelocities();
   std::vector<Vec3> forces               = state.getForces();
   double kineticEnergy                   = state.getKineticEnergy();
   double potentialEnergy                 = state.getPotentialEnergy();

   // compute norm of force

   double forceNorm         = 0.0;
   for( unsigned int ii = 0; ii < forces.size(); ii++ ){

#if 0
(void) fprintf( log, "%6u x[%14.7e %14.7e %14.7e] f[%14.7e %14.7e %14.7e]\n", ii,
                coordinates[ii][0], coordinates[ii][1], coordinates[ii][2],
                forces[ii][0], forces[ii][1], forces[ii][2] );
#endif

      forceNorm += forces[ii][0]*forces[ii][0] + forces[ii][1]*forces[ii][1] + forces[ii][2]*forces[ii][2];
   }

   // check norm is not nan

   if( isinf( forceNorm ) || isnan( forceNorm ) ){ 
      if( log ){
         (void) fprintf( log, "%s norm of force is nan -- aborting.\n", methodName.c_str() );
         unsigned int hitNan = 0;
         for( unsigned int ii = 0; (ii < forces.size()) && (hitNan < 10); ii++ ){

            if( isinf( forces[ii][0] ) || isnan( forces[ii][0] ) ||
                isinf( forces[ii][1] ) || isnan( forces[ii][1] ) ||
                isinf( forces[ii][2] ) || isnan( forces[ii][2] ) )hitNan++;

            (void) fprintf( log, "%6u x[%14.7e %14.7e %14.7e] f[%14.7e %14.7e %14.7e]\n", ii,
                            coordinates[ii][0], coordinates[ii][1], coordinates[ii][2],
                            forces[ii][0], forces[ii][1], forces[ii][2] );
         }
         char buffer[1024];
         (void) sprintf( buffer, "%s : nans detected -- aborting.\n", methodName.c_str() );
         throwException(__FILE__, __LINE__, buffer );
      }    
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
   double denominator              = forceNorm; 
   if( denominator > 0.0 ){
      difference /= denominator;
   }

   if( log ){
      (void) fprintf( log, "%s  difference=%14.8e dE=%14.8e Pe2/1 [%16.10e %16.10e] delta=%10.4e nrm=%16.10e\n",
                      methodName.c_str(), difference, deltaEnergy, perturbedPotentialEnergy,
                      potentialEnergy, delta, forceNorm );
      (void) fflush( log );
   }
   if( summaryFile ){
      std::string forceString;
      if( forceStringArray.size() > 5 ){
         forceString = "All";
      } else {
         for( StringVectorCI ii = forceStringArray.begin(); ii != forceStringArray.end(); ii++ ){
            forceString += (*ii) + "_";
         }
      }
      if( forceString.size() < 1 ){
         forceString = "NA";
      }
      (void) fprintf( summaryFile, "EnergyForceConsistent %s\nForce %s\nCalculated %14.6e\nExpected %14.7e\nDiffNorm %14.7e\nE0 %14.7e\nE1 %14.7e\nForceNorm %14.7e\nDelta %14.7e\n",
                      context.getPlatform().getName().c_str(), forceString.c_str(), deltaEnergy, forceNorm, difference, potentialEnergy, perturbedPotentialEnergy, forceNorm, delta );
   }

   if( applyAssertion ){
      ASSERT( difference < tolerance );
      if( log ){
         (void) fprintf( log, "\n%s passed\n", methodName.c_str() );
         (void) fflush( log );
      }
   }
   return returnStatus;

}


/**---------------------------------------------------------------------------------------

   Find stats for vec3

   @param array                 array 
   @param statVector              vector of stats 

   @return 0

   --------------------------------------------------------------------------------------- */

int compareForces( const std::vector<Vec3>& forceArray1, const std::string& f1Name, std::vector<double>& forceArray1Sum, std::vector<double>& forceArray1Stats,
                   const std::vector<Vec3>& forceArray2, const std::string& f2Name, std::vector<double>& forceArray2Sum, std::vector<double>& forceArray2Stats,
                   double* maxDelta, int* maxDeltaIndex, double* maxRelativeDelta, int* maxRelativeDeltaIndex, double* maxDot, double forceTolerance, FILE* inputLog ){

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
   *maxDot                             = -1.0e+30;
   *maxDeltaIndex                      = -1;
   *maxRelativeDeltaIndex              = -1;

   std::vector<double> forceArray1Norms;
   std::vector<double> forceArray2Norms;

   forceArray1Sum.resize( 3 );
   forceArray2Sum.resize( 3 );
   for( unsigned int ii = 0; ii < 3; ii++ ){
      forceArray1Sum[ii] = forceArray2Sum[ii] = 0.0;
   }

   (void) fprintf( log, "   Id     delta  relDelta       dot %4s     norm                                          force    %4s     norm                                         force\n",
                   f1Name.c_str(), f2Name.c_str() ); 
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
      double dotProduct      = f1[0]*f2[0] + f1[1]*f2[1] + f1[2]*f2[2];
      if(  (normF1*normF2) > 0.0 ){
         dotProduct     /= (normF1*normF2);
         dotProduct      = 1.0 - dotProduct;
      } else {
         dotProduct      = 0.0;
      }

      double relativeDelta   = (delta*2.0)/(normF1+normF2);

      int print              = 0;
      if( delta > forceTolerance ){
         print++;
      }

      if( *maxRelativeDelta < relativeDelta ){
         print++;
         *maxRelativeDelta        = relativeDelta;   
         *maxRelativeDeltaIndex   = static_cast<int>(ii);
      }

      if( *maxDot < dotProduct ){
         *maxDot = dotProduct;   
         if( dotProduct > 1.0e-06 )print++;
      }

      if( *maxDelta < delta ){
         *maxDelta      = delta;
         *maxDeltaIndex = static_cast<int>(ii);
      }

      if( print && log ){
//         (void) fprintf( log, "%5d delta=%9.3e relDelta=%9.3e dot=%9.3e  %s %13.7e [%14.7e %14.7e %14.7e]  %s %13.7e [%14.7e %14.7e %14.7e]\n", 
//                         ii, delta, relativeDelta, dotProduct, f1Name.c_str(), normF1, f1[0], f1[1], f1[2], f2Name.c_str(), normF2, f2[0], f2[1], f2[2] );
         (void) fprintf( log, "%5d %9.3e %9.3e %9.3e %13.7e [%14.7e %14.7e %14.7e]    %13.7e [%14.7e %14.7e %14.7e] %s\n", 
                         ii, delta, relativeDelta, dotProduct, normF1, f1[0], f1[1], f1[2], normF2, f2[0], f2[1], f2[2], ((normF1 > 1.0e+06 || normF2 > 1.0e+06) ? "!!!" : "") );
         (void) fflush( log );
      }
   }

   findStatsForDouble( forceArray1Norms, forceArray1Stats );
   findStatsForDouble( forceArray2Norms, forceArray2Stats );

   return 0;
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

static int checkForcesDuringSimulation( int currentStep, Context& cudaContext, Context& referenceContext, FILE* log ) {
    
// ---------------------------------------------------------------------------------------

   static const std::string methodName             = "checkForcesDuringSimulation";

// ---------------------------------------------------------------------------------------

   _synchContexts( cudaContext, referenceContext );

   State referenceState                            = referenceContext.getState( State::Energy | State::Forces );
   double referenceKineticEnergy                   = referenceState.getKineticEnergy();
   double referencePotentialEnergy                 = referenceState.getPotentialEnergy();
   double referenceTotalEnergy                     = referenceKineticEnergy + referencePotentialEnergy;

   State cudaState                                 = cudaContext.getState( State::Energy | State::Forces );
   double cudaKineticEnergy                        = cudaState.getKineticEnergy();
   double cudaPotentialEnergy                      = cudaState.getPotentialEnergy();
   double cudaTotalEnergy                          = cudaKineticEnergy + cudaPotentialEnergy;

   (void) fprintf( log, "%6d PE=%14.7e %14.7e KE=%14.7e %14.7e E=%14.7e %14.7ed\n",
                   currentStep, referencePotentialEnergy, cudaPotentialEnergy,
                                referenceKineticEnergy,   cudaKineticEnergy, 
                                referenceTotalEnergy, cudaTotalEnergy );

   // compare reference vs cuda forces

   std::vector<Vec3> referenceForces               = referenceState.getForces();
   std::vector<Vec3> cudaForces                    = cudaState.getForces();

   double maxDeltaRefCud                           = -1.0e+30;
   double maxRelativeDeltaRefCud                   = -1.0e+30;
   double maxDotRefCud                             = -1.0e+30;
   double maxDeltaPrmCud                           = -1.0e+30;
   double maxRelativeDeltaPrmCud                   = -1.0e+30;
   double maxDotPrmCud                             = -1.0e+30;
   double forceTolerance                           = 1.0e-01;
   int maxDeltaIndex;
   int maxRelativeDeltaRefCudIndex;
   
   std::vector<double> forceArray1Sum;
   std::vector<double> forceArray2Sum;
   std::vector<double> forceArray3Sum;
   
   std::vector<double> referenceForceStats;
   std::vector<double> cudaForceStats;
   
   compareForces( referenceForces, "fRef", forceArray1Sum, referenceForceStats,
                  cudaForces,      "fCud", forceArray2Sum, cudaForceStats, 
                  &maxDeltaRefCud, &maxDeltaIndex, &maxRelativeDeltaRefCud, &maxRelativeDeltaRefCudIndex, &maxDotRefCud, forceTolerance, log );

   (void) fprintf( log, "MaxDelta=%13.7e at %d MaxRelativeDelta=%13.7e at %d maxDotRefCud=%14.6e\n",
                   maxDeltaRefCud, maxDeltaIndex, maxRelativeDeltaRefCud, maxRelativeDeltaRefCudIndex, maxDotRefCud );
   (void) fprintf( log, "Reference force average=%14.7e stddev=%14.7e min=%14.7e at %6.0f max=%14.7e at %6.0f\n",
                   referenceForceStats[0], referenceForceStats[1], referenceForceStats[2], referenceForceStats[3],
                   referenceForceStats[4], referenceForceStats[5] );

   (void) fprintf( log, "     Cuda force average=%14.7e stddev=%14.7e min=%14.7e at %6.0f max=%14.7e at %6.0f\n",
                   cudaForceStats[0], cudaForceStats[1], cudaForceStats[2], cudaForceStats[3],
                   cudaForceStats[4], cudaForceStats[5] );

   (void) fflush( log );

   return 0;

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

static int checkEnergyConservation( Context& context,  MapStringString& inputArgumentMap, FILE* log,
                                    FILE* summaryFile ) {
    
// ---------------------------------------------------------------------------------------

   static const std::string methodName             = "checkEnergyConservation";

   // tolerance for thermostat

   double temperatureTolerance                     = 3.0;

   // tolerance for energy conservation test

   double energyTolerance                          = 0.05;

   std::string equilibrationIntegratorName         = "LangevinIntegrator";
   //std::string equilibrationIntegratorName         = "VerletIntegrator";
   int equilibrationTotalSteps                     = 1000;
   double equilibrationStepsBetweenReportsRatio    = 0.1;
   double equilibrationTimeStep                    = 0.002;
   double equilibrationFriction                    = 91.0;
   double equilibrationShakeTolerance              = 1.0e-05;
   double equilibrationErrorTolerance              = 1.0e-05;
   double equilibrationTemperature                 = 300.0;
   int equilibrationSeed                           = 1993;
   int equilibrationWriteContext                   = 0;

   std::string simulationIntegratorName            = "VerletIntegrator";
   int simulationTotalSteps                        = 10000;
   double simulationStepsBetweenReportsRatio       = 0.01;
   double simulationTimeStep                       = 0.001;
   double simulationFriction                       = 91.0;
   double simulationShakeTolerance                 = 1.0e-06;
   double simulationErrorTolerance                 = 1.0e-05;
   double simulationTemperature                    = 300.0;
   int simulationSeed                              = 1993;
   int simulationWriteContext                      = 0;

   int applyAssertion                              = 1;
   std::string deviceId                            = "0";
   std::string runId                               = "RunId";

// ---------------------------------------------------------------------------------------

   setIntFromMap(    inputArgumentMap, "applyAssertion",                           applyAssertion                         );
   setStringFromMap( inputArgumentMap, "cudaDeviceId",                             deviceId                               );
   setStringFromMap( inputArgumentMap, "runId",                                    runId                                  );

   setStringFromMap( inputArgumentMap, "equilibrationIntegrator",                  equilibrationIntegratorName            );
   setIntFromMap(    inputArgumentMap, "equilibrationTotalSteps",                  equilibrationTotalSteps                );
   setDoubleFromMap( inputArgumentMap, "equilibrationStepsBetweenReportsRatio",    equilibrationStepsBetweenReportsRatio  );
   setDoubleFromMap( inputArgumentMap, "equilibrationTimeStep",                    equilibrationTimeStep                  );
   setDoubleFromMap( inputArgumentMap, "equilibrationFriction",                    equilibrationFriction                  );
   setDoubleFromMap( inputArgumentMap, "equilibrationShakeTolerance",              equilibrationShakeTolerance            );
   setDoubleFromMap( inputArgumentMap, "equilibrationErrorTolerance",              equilibrationErrorTolerance            );
   setDoubleFromMap( inputArgumentMap, "equilibrationTemperature",                 equilibrationTemperature               );
   setIntFromMap(    inputArgumentMap, "equilibrationSeed",                        equilibrationSeed                      );
   setIntFromMap(    inputArgumentMap, "equilibrationWriteContext",                equilibrationWriteContext              );

   setStringFromMap( inputArgumentMap, "simulationIntegrator",                     simulationIntegratorName               );
   setIntFromMap(    inputArgumentMap, "simulationTotalSteps",                     simulationTotalSteps                   );
   setDoubleFromMap( inputArgumentMap, "simulationStepsBetweenReportsRatio",       simulationStepsBetweenReportsRatio     );
   setDoubleFromMap( inputArgumentMap, "simulationTimeStep",                       simulationTimeStep                     );
   setDoubleFromMap( inputArgumentMap, "simulationFriction",                       simulationFriction                     );
   setDoubleFromMap( inputArgumentMap, "simulationShakeTolerance",                 simulationShakeTolerance               );
   setDoubleFromMap( inputArgumentMap, "simulationErrorTolerance",                 simulationErrorTolerance               );
   setDoubleFromMap( inputArgumentMap, "simulationTemperature",                    simulationTemperature                  );
   setIntFromMap(    inputArgumentMap, "simulationSeed",                           simulationSeed                         );
   setIntFromMap(    inputArgumentMap, "simulationWriteContext",                   simulationWriteContext                 );

   if( log ){
      (void) fprintf( log, "%s Equilbration: %s steps=%d ratioRport=%.2f timeStep=%.4f T=%8.3f friction=%8.3f\n"
                           "ShakeTol=%3e ErrorTol=%.3e seed=%d\n", methodName.c_str(),
                      equilibrationIntegratorName.c_str(), equilibrationTotalSteps, equilibrationStepsBetweenReportsRatio,
                      equilibrationTimeStep, equilibrationTemperature, equilibrationFriction,
                      equilibrationShakeTolerance, equilibrationErrorTolerance, equilibrationSeed );

      (void) fprintf( log, "%s Simulation: %s steps=%d ratioRport=%.2f timeStep=%.4f T=%8.3f friction=%8.3f\n"
                           "ShakeTol=%3e ErrorTol=%.3e seed=%d\n", methodName.c_str(),
                      simulationIntegratorName.c_str(), simulationTotalSteps, simulationStepsBetweenReportsRatio,
                      simulationTimeStep, simulationTemperature, simulationFriction,
                      simulationShakeTolerance, simulationErrorTolerance, simulationSeed );
      (void) fprintf( log, "deviceId=%s applyAssertion=%d\n", deviceId.c_str(), applyAssertion );
      (void) fflush( log );
   }

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
                                            equilibrationFriction, equilibrationTemperature,
                                            equilibrationShakeTolerance, equilibrationErrorTolerance, equilibrationSeed, log );

   if( log ){
      _printIntegratorInfo( integrator, log );
   }
   Context*       equilibrationContext = _getContext( &system, &context, integrator, "CudaPlatform", "EquilibrationContext", deviceId, log );

   // equilibration loop

   int constraintViolations               = 0;
   int constraintChecks                   = 0;

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

      // get energies, check for constraint violations and nans

      State state                            = equilibrationContext->getState( State::Energy | State::Forces );
   
      double kineticEnergy                   = state.getKineticEnergy();
      double potentialEnergy                 = state.getPotentialEnergy();
      double totalEnergy                     = kineticEnergy + potentialEnergy;
      double maxViolation;
      int violations                         = checkConstraints( *equilibrationContext, system, equilibrationShakeTolerance, &maxViolation, log );
      constraintViolations                  += violations;
      constraintChecks++;
      if( log ){
         (void) fprintf( log, "Equilibration: %6d KE=%14.7e PE=%14.7e E=%14.7e violations=%6d max=%13.6e totalViolation=%6d\n",
                         currentStep, kineticEnergy, potentialEnergy, totalEnergy, violations, maxViolation, constraintViolations );
         (void) fflush( log );
      }

      // compare reference and gpu forces, if violations found

      if( violations && log ){
         checkForcesDuringSimulation( currentStep, *equilibrationContext, context, log );
      }

      // output context?

      if( equilibrationWriteContext ){
         std::stringstream fileName;
         fileName << "EquilCnxt_" << runId << "_" << currentStep << ".txt";
         writeContextToFile( fileName.str(), *equilibrationContext, (State::Positions | State::Velocities | State::Forces  | State::Energy), log );
      }

      // nans

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
                                                      simulationFriction, simulationTemperature,
                                                      simulationShakeTolerance, simulationErrorTolerance,
                                                      simulationSeed, log );

   if( log ){
      _printIntegratorInfo( simulationIntegrator, log );
   }

   //delete equilibrationContext;
   Context*       simulationContext = _getContext( &system, equilibrationContext, simulationIntegrator, "CudaPlatform", "SimulationContext", deviceId, log );
   //Context*       simulationContext = _getContext( &system, equilibrationContext, simulationIntegrator, "ReferencePlatform", "SimulationContext", deviceId, log );
   //Context*       simulationContext = equilibrationContext;
   //Context*       simulationContext = _getContext( &system, &context, simulationIntegrator, "CudaPlatform", "SimulationContext", deviceId, log );
   //Context*       simulationContext = _getContext( &system, &context, simulationIntegrator, "ReferencePlatform", "SimulationContext", deviceId, log );
   //Context*       simulationContext = _getContext( &system, &context, simulationIntegrator, "ReferencePlatform", "SimulationContext", deviceId, log );

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
   currentStep                         = 0;

   if( log ){
      (void) fprintf( log, "simulationTotalSteps=%d simulationStepsBetweenReports=%d ratio=%.4f\n", 
                      simulationTotalSteps, simulationStepsBetweenReports, simulationStepsBetweenReportsRatio );
      (void) fflush( log );
   }

   // write initial context

   if( simulationWriteContext ){
      std::stringstream fileName;
      fileName << "SimulCnxt_" << runId << "_" << currentStep << ".txt";
      writeContextToFile( fileName.str(), *simulationContext, (State::Positions | State::Velocities | State::Forces  | State::Energy), log );
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
      double maxViolation;

      int violations                         = checkConstraints( *simulationContext, system, simulationShakeTolerance, &maxViolation, log );
      constraintViolations                  += violations;
      constraintChecks++;

      stepIndexArray.push_back( (double) currentStep );
      kineticEnergyArray.push_back( kineticEnergy );
      potentialEnergyArray.push_back( potentialEnergy );
      totalEnergyArray.push_back( totalEnergy );

      // diagnostics & check for nans

      if( log ){
         (void) fprintf( log, "Simulation: %6d KE=%14.7e PE=%14.7e E=%14.7e violations=%6d max=%13.6e totalViolation=%6d\n",
                         currentStep, kineticEnergy, potentialEnergy, totalEnergy, violations, maxViolation, constraintViolations );
         (void) fflush( log );
      }

      if( violations && log ){
         checkForcesDuringSimulation( currentStep, *simulationContext, context, log );
      }

      // output context?

      if( simulationWriteContext ){
         std::stringstream fileName;
         fileName << "SimulCnxt_" << runId << "_" << currentStep << ".txt";
         writeContextToFile( fileName.str(), *simulationContext, (State::Positions | State::Velocities | State::Forces  | State::Energy), log );
      }

      // check nans

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

   if( summaryFile ){
      (void) fprintf( summaryFile, "Platform %s\nIntegrator %s\nSteps %d\nTimeStepSize %14.7e\nAtoms %d\n",
// crashes???
//                      simulationContext->getPlatform().getName().c_str(), 
                      "Cuda", simulationIntegratorName.c_str(), simulationTotalSteps, simulationTimeStep, numberOfAtoms );
   }

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

      // summary info

      if( summaryFile ){
         double totalTime           = static_cast<double>(totalSimulationTime)/static_cast<double>(CLOCKS_PER_SEC);
         double timePerStep         = totalTime/static_cast<double>(simulationTotalSteps);
         (void) fprintf( summaryFile, "T %14.7e\nCalcT %14.7e\nStddevT %14.7e\nMinT %14.7e\nMaxT %14.7e\n",
                         simulationTemperature, temperatureStatistics[0], temperatureStatistics[1], temperatureStatistics[2],
                         temperatureStatistics[4] );

      }

      // check that <temperature> is within tolerance

      if( applyAssertion ){
         ASSERT_EQUAL_TOL( temperatureStatistics[0], simulationTemperature, temperatureTolerance );
      }

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

      // summary info

      if( summaryFile ){
         double totalTime           = static_cast<double>(totalSimulationTime)/static_cast<double>(CLOCKS_PER_SEC);
         double timePerStep         = totalTime/static_cast<double>(simulationTotalSteps);
         (void) fprintf( summaryFile, "DriftE %14.7e\nAvgE %14.7e\nStddevE %14.7e\n"
                         "Dof %d\nMinE %14.7e\nMinE_Idx %d\nMaxE %14.7e\nMax_E_Idx %d\n",
                          stddevE,                                      // drift
                          statistics[0], statistics[1],                 // mean & stddev
                          (3*numberOfAtoms - system.getNumConstraints() - 3), // dof
                          statistics[2], (int) (statistics[3] + 0.001), // min & index
                          statistics[4], (int) (statistics[5] + 0.001) ); // max index

      }

      // check that energy fluctuation is within tolerance
  
      if( applyAssertion ){
         ASSERT_EQUAL_TOL( stddevE, 0.0, energyTolerance );
      }

   }

   // summary info

   if( summaryFile ){
      double totalTime           = static_cast<double>(totalSimulationTime)/static_cast<double>(CLOCKS_PER_SEC);
      double timePerStep         = totalTime/static_cast<double>(simulationTotalSteps);
      (void) fprintf( summaryFile, "ConstraintViolations %d\nConstraintChecks %d\nWallTime %.3e\nWallTimePerStep %.3e\n",
                      constraintViolations, constraintChecks, totalTime, timePerStep );
       (void) fclose( summaryFile );
   }

   if( applyAssertion && log ){
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

Context* testSetup( std::string parameterFileName, MapStringInt& forceMap, Platform& platform, std::vector<Vec3>& forces,
                    double* kineticEnergy, double* potentialEnergy, MapStringVectorOfVectors& supplementary,
                     MapStringString& inputArgumentMap, FILE* log ){

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
      readParameterFile( parameterFileName, forceMap, *system, coordinates, velocities,
                         forces, kineticEnergy, potentialEnergy, supplementary, inputArgumentMap, log );

   Context* context;
   context = new Context( *system, *integrator, platform);

   StringVector forceStringArray;
   getForceStrings( *system, forceStringArray, log );

   if( log ){
      (void) fprintf( log, "\n%s Active Forces:\n", methodName.c_str() );
      for( StringVectorCI ii = forceStringArray.begin(); ii != forceStringArray.end(); ii++ ){
         (void) fprintf( log, "   %s\n", (*ii).c_str() );
      }
      (void) fflush( log );
   }

   // read context if present in inputArgumentMap

   MapStringStringI readContext = inputArgumentMap.find( "readContext" );
   if( readContext != inputArgumentMap.end() ){
      readContextFromFile( (*readContext).second, *context, (State::Positions | State::Velocities), log );
   } else {
      context->setPositions( coordinates );
      context->setVelocities( velocities );
   }

   return context;

}
void testReferenceCudaForces( std::string parameterFileName, MapStringInt& forceMap,
                              MapStringString& inputArgumentMap, FILE* inputLog, FILE* summaryFile ){

// ---------------------------------------------------------------------------------------

  static const std::string methodName      = "testReferenceCudaForces";
  int PrintOn                              = 1; 
  int compareParameterForces               = 0; 

  double forceTolerance                    = 0.01; 
  double energyTolerance                   = 0.01; 
  int numberOfSteps                        = 1; 
  int steps                                = 0; 
  int applyAssertion                       = 1;

// ---------------------------------------------------------------------------------------

   FILE* log;
   if( PrintOn == 0 && inputLog ){
      log = NULL;
   } else {
      log = inputLog;
   } 

   setIntFromMap(    inputArgumentMap, "applyAssertion",     applyAssertion    );

   if( log ){
      (void) fprintf( log, "%s force tolerance=%.3e energy tolerance=%.3e step=%d\n",
                      methodName.c_str(), forceTolerance, energyTolerance, numberOfSteps );
      (void) fflush( log );
   }   

   Platform& referencePlatform         =  Platform::getPlatformByName("Reference");
   Platform& cudaPlatform              = Platform::getPlatformByName("Cuda");;

   double parameterKineticEnergy, parameterPotentialEnergy;

   std::vector<Vec3> parameterForces;
   std::vector<Vec3> parameterForces2;
   MapStringVectorOfVectors supplementary;

   Context* referenceContext       = testSetup( parameterFileName, forceMap, referencePlatform, 
                                                parameterForces, &parameterKineticEnergy, &parameterPotentialEnergy,
                                                supplementary, inputArgumentMap, log );
   Context* cudaContext            = testSetup( parameterFileName, forceMap,  cudaPlatform,
                                                parameterForces2, &parameterKineticEnergy, &parameterPotentialEnergy,
                                                supplementary, inputArgumentMap, log );

   Integrator& referenceIntegrator = referenceContext->getIntegrator();
   Integrator& cudaIntegrator      = cudaContext->getIntegrator();

   // Run several steps and see if relative force difference is within tolerance
   
   int includeEnergy = 1;
   for( int step = 0; step < numberOfSteps; step++ ){

      // pull info out of contexts

      int types                                       = State::Positions | State::Velocities | State::Forces;
      if( includeEnergy ){
         types                                       |= State::Energy;
      }

      State cudaState                                 =      cudaContext->getState( types );
      State referenceState                            = referenceContext->getState( types );

      std::vector<Vec3> referenceCoordinates          = referenceState.getPositions();
      std::vector<Vec3> referenceVelocities           = referenceState.getVelocities();
      std::vector<Vec3> referenceForces               = referenceState.getForces();

      double referenceKineticEnergy;
      double referencePotentialEnergy;
      if( includeEnergy ){
         referenceKineticEnergy                       = referenceState.getKineticEnergy();
         referencePotentialEnergy                     = referenceState.getPotentialEnergy();
      } else {
         referenceKineticEnergy                       = 0.0;
         referencePotentialEnergy                     = 0.0;
      }

      std::vector<Vec3> cudaCoordinates               = cudaState.getPositions();
      std::vector<Vec3> cudaVelocities                = cudaState.getVelocities();
      std::vector<Vec3> cudaForces                    = cudaState.getForces();

      double cudaKineticEnergy                        = cudaState.getKineticEnergy();
      double cudaPotentialEnergy                      = cudaState.getPotentialEnergy();
      if( includeEnergy ){
         cudaKineticEnergy                            = cudaState.getKineticEnergy();
         cudaPotentialEnergy                          = cudaState.getPotentialEnergy();
      } else {
         cudaKineticEnergy                            = 0.0;
         cudaPotentialEnergy                          = 0.0;
      }

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
               for( unsigned int ii = 0; ii < referenceForces.size(); ii++ ){
                  (void) fprintf( log, "%6u 0[%14.7e %14.7e %14.7e] 1[%14.7e %14.7e %14.7e] 2[%14.7e %14.7e %14.7e]\n", ii,
                                  referenceForces[ii][0], cudaForces[ii][0], parameterForces[ii][0],
                                  referenceForces[ii][1], cudaForces[ii][1], parameterForces[ii][1],
                                  referenceForces[ii][2], cudaForces[ii][2], parameterForces[ii][2] );
                  if( ii == maxPrint ){
                      ii = referenceForces.size()- maxPrint;
                      if( ii < maxPrint )ii = maxPrint;
                  }
               }
            } else {
               for( unsigned int ii = 0; ii < referenceForces.size(); ii++ ){
                  (void) fprintf( log, "%6u 0[%14.7e %14.7e] 1[%14.7e %14.7e] 2[%14.7e %14.7e]\n", ii,
                                  referenceForces[ii][0], cudaForces[ii][0],
                                  referenceForces[ii][1], cudaForces[ii][1],
                                  referenceForces[ii][2], cudaForces[ii][2]  );
                  if( ii == maxPrint ){
                      ii = referenceForces.size() - maxPrint;
                      if( ii < maxPrint )ii = maxPrint;
                  }
               }
            }

         } else { 

            if( compareParameterForces ){
               for( unsigned int ii = 0; ii < referenceForces.size(); ii++ ){
                  (void) fprintf( log, "%6u r[%14.7e %14.7e %14.7e] c[%14.7e %14.7e %14.7e] p[%14.7e %14.7e %14.7e]\n", ii,
                                  referenceForces[ii][0], referenceForces[ii][1], referenceForces[ii][2],
                                  cudaForces[ii][0], cudaForces[ii][1], cudaForces[ii][2],
                                  parameterForces[ii][0], parameterForces[ii][1], parameterForces[ii][2] );
                  if( ii == maxPrint ){
                      ii = referenceForces.size() - maxPrint;
                      if( ii < maxPrint )ii = maxPrint;
                  }
               }
            } else {
               for( unsigned int ii = 0; ii < referenceForces.size(); ii++ ){
                  (void) fprintf( log, "%6u r[%14.7e %14.7e %14.7e] c[%14.7e %14.7e %14.7e]\n", ii,
                                  referenceForces[ii][0], referenceForces[ii][1], referenceForces[ii][2],
                                  cudaForces[ii][0], cudaForces[ii][1], cudaForces[ii][2] );
                  if( ii == maxPrint ){
                      ii = referenceForces.size() - maxPrint;
                      if( ii < maxPrint )ii = maxPrint;
                  }
               }
            }
         }
      
      }

      // compare reference vs cuda forces

      double maxDeltaRefCud                          = -1.0e+30;
      double maxRelativeDeltaRefCud                  = -1.0e+30;
      double maxDotRefCud                            = -1.0e+30;
      double maxDeltaPrmCud                          = -1.0e+30;
      double maxRelativeDeltaPrmCud                  = -1.0e+30;
      double maxDotPrmCud                            = -1.0e+30;
      int maxDeltaIndex;
      int maxRelativeDeltaRefCudIndex;

      std::vector<double> forceArray1Sum;
      std::vector<double> forceArray2Sum;
      std::vector<double> forceArray3Sum;

      std::vector<double> referenceForceStats;
      std::vector<double> cudaForceStats;
      std::vector<double> cudaForceStats1;
      std::vector<double> paramForceStats;

      compareForces( referenceForces, "fRef", forceArray1Sum, referenceForceStats,
                     cudaForces,      "fCud", forceArray2Sum, cudaForceStats, 
                     &maxDeltaRefCud, &maxDeltaIndex, &maxRelativeDeltaRefCud,
                     &maxRelativeDeltaRefCudIndex, &maxDotRefCud, forceTolerance, log );
      
      (void) fflush( log );

      if( compareParameterForces ){

         // compare cuda & forces retreived from parameter file

/*
         compareForces( parameterForces, "fPrm", forceArray3Sum, paramForceStats,
                        cudaForces,      "fCud", forceArray2Sum, cudaForceStats1,
                        &maxDeltaPrmCud, &maxRelativeDeltaPrmCud, &maxDotPrmCud, forceTolerance, log );
*/
      }

      // summary file info

      if( summaryFile ){

         StringVector forceStringArray;
         System system = referenceContext->getSystem();
         getForceStrings( system, forceStringArray, log );
         std::string forceString;
         if( forceStringArray.size() > 5 ){
            forceString = "All";
         } else {
            for( StringVectorCI ii = forceStringArray.begin(); ii != forceStringArray.end(); ii++ ){
               forceString += *ii;
            }
         }
         if( forceString.size() < 1 ){
            forceString = "NA";
         }
         (void) fprintf( summaryFile, "Force %s\nAtoms %u\nMaxDelta %14.7e\nMaxRelDelta %14.7e\nMaxDot %14.7e\n",
                         forceString.c_str(), referenceForces.size(), maxDeltaRefCud, maxRelativeDeltaRefCud, maxDotRefCud);

         double sum = ( fabs(forceArray1Sum[0] ) + fabs( forceArray1Sum[1] ) + fabs( forceArray1Sum[2]) )*0.33333;
         (void) fprintf( summaryFile, "SumRef %14.7e\n", sum );

                sum = ( fabs(forceArray2Sum[0] ) + fabs( forceArray2Sum[1] ) + fabs( forceArray2Sum[2]) )*0.33333;
         (void) fprintf( summaryFile, "SumCuda %14.7e\n", sum );
         double difference         = fabs( referencePotentialEnergy - cudaPotentialEnergy );
         double relativeDifference = difference/( fabs( referencePotentialEnergy ) + fabs(  cudaPotentialEnergy ) + 1.0e-10);
         (void) fprintf( summaryFile, "RefPE %14.7e\nCudaPE %14.7e\nDiffPE %14.7e\nRelDiffPE %14.7e\n",
                         referencePotentialEnergy, cudaPotentialEnergy, difference, relativeDifference );
      }

      if( log ){
         (void) fprintf( log, "max delta=%13.7e at %d maxRelDelta=%13.7e at %d maxDot=%14.7e\n",
                         maxDeltaRefCud, maxDeltaIndex, maxRelativeDeltaRefCud, maxRelativeDeltaRefCudIndex, maxDotRefCud );
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

      if( applyAssertion ){
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
         _synchContexts( *cudaContext, *referenceContext );
      }

   }

   if( log ){
      if( applyAssertion ){
         (void) fprintf( log, "\n%s tests passed\n", methodName.c_str() );
      } else {
         (void) fprintf( log, "\n%s tests off\n", methodName.c_str() );
      }
      (void) fflush( log );
   }
}

void testInputForces( std::string parameterFileName, MapStringInt& forceMap,
                      MapStringString& inputArgumentMap, FILE* inputLog, FILE* summaryFile ){

// ---------------------------------------------------------------------------------------

  static const std::string methodName      = "testInputForces";
  int PrintOn                              = 1; 
  int compareParameterForces               = 0; 

  double forceTolerance                    = 0.01; 
  double energyTolerance                   = 0.01; 
  int numberOfSteps                        = 1; 
  int steps                                = 0; 
  int applyAssertion                       = 1;
  std::string inputForceToCompare;

// ---------------------------------------------------------------------------------------

   FILE* log;
   if( PrintOn == 0 && inputLog ){
      log = NULL;
   } else {
      log = inputLog;
   } 

   setIntFromMap(       inputArgumentMap, "applyAssertion",          applyAssertion    );
   if( setStringFromMap(    inputArgumentMap, "inputForceToCompare",     inputForceToCompare ) == 0 ){
      if( log ){
         (void) fprintf( log, "%s inputForceToCompare field not set.\n", methodName.c_str() );
         (void) fflush( log );
      }   
      return;
   }

   if( log ){
      (void) fprintf( log, "%s force tolerance=%.3e energy tolerance=%.3e step=%d\n",
                      methodName.c_str(), forceTolerance, energyTolerance, numberOfSteps );
      (void) fflush( log );
   }   

   ReferencePlatform referencePlatform;
   double parameterKineticEnergy, parameterPotentialEnergy;

   std::vector<Vec3> parameterForces;
   std::vector<Vec3> parameterForces2;
   MapStringVectorOfVectors supplementary;

   Context* referenceContext       = testSetup( parameterFileName, forceMap, referencePlatform, 
                                                parameterForces, &parameterKineticEnergy, &parameterPotentialEnergy,
                                                supplementary, inputArgumentMap, log );

   MapStringVectorOfVectorsI forceVectorI = supplementary.find( inputForceToCompare );
   if(  forceVectorI == supplementary.end() ){
      if( log ){
         (void) fprintf( log, "%s inputForceToCompare=<%s> is missing.\n", methodName.c_str(), inputForceToCompare.c_str() );
         (void) fflush( log );
      }   
      return;
   }
   VectorOfVectors forceVectorToCompare = (*forceVectorI).second;

   Integrator& referenceIntegrator = referenceContext->getIntegrator();

   // Run several steps and see if relative force difference is within tolerance
   
#if 0
   for( int step = 0; step < numberOfSteps; step++ ){

      // pull info out of contexts

      int types                                       = State::Positions | State::Velocities | State::Forces | State::Energy;

      State referenceState                            = referenceContext->getState( types );

      std::vector<Vec3> referenceCoordinates          = referenceState.getPositions();
      std::vector<Vec3> referenceVelocities           = referenceState.getVelocities();
      std::vector<Vec3> referenceForces               = referenceState.getForces();
      double referenceKineticEnergy                   = referenceState.getKineticEnergy();
      double referencePotentialEnergy                 = referenceState.getPotentialEnergy();

      // diagnostics

      if( log ){
         //static const unsigned int maxPrint = MAX_PRINT;
         static const unsigned int maxPrint   = 1000000;

         // print x,y,z components separately, if formatType == 1
         // else print reference, cuda and parameter forces in blocks of 3

         static const unsigned int formatType = 1;

         (void) fprintf( log, "%s\n", methodName.c_str() );
#if 0
         if( compareParameterForces ){
            (void) fprintf( log, "Kinetic   energies: r=%14.7e c=%14.7e, p=%14.7e\n", referenceKineticEnergy, cudaKineticEnergy, parameterKineticEnergy );
            (void) fprintf( log, "Potential energies: r=%14.7e c=%14.7e, p=%14.7e\n", referencePotentialEnergy, cudaPotentialEnergy, parameterPotentialEnergy );
            (void) fprintf( log, "Sample of forces: %u (r=reference, c=cuda, p=parameter) file forces\n", referenceForces.size() );
         } else {
            (void) fprintf( log, "Kinetic   energies: r=%14.7e c=%14.7e\n", referenceKineticEnergy, cudaKineticEnergy );
            (void) fprintf( log, "Potential energies: r=%14.7e c=%14.7e\n", referencePotentialEnergy, cudaPotentialEnergy );
            (void) fprintf( log, "Sample of forces: %u (r=reference, c=cuda) file forces\n", referenceForces.size() );
         }
#endif

            for( unsigned int ii = 0; ii < referenceForces.size() && ii < maxPrint; ii++ ){
               (void) fprintf( log, "%6u 0[%14.7e %14.7e] 1[%14.7e %14.7e] 2[%14.7e %14.7e]\n", ii,
                               referenceForces[ii][0], forceVectorToCompare[ii][0],
                               referenceForces[ii][1], forceVectorToCompare[ii][1],
                               referenceForces[ii][2], forceVectorToCompare[ii][2]  );
            }
            if( referenceForces.size() > maxPrint ){
               for( unsigned int ii = referenceForces.size() - maxPrint; ii < referenceForces.size(); ii++ ){
                  (void) fprintf( log, "%6u 0[%14.7e %14.7e] 1[%14.7e %14.7e] 2[%14.7e %14.7e]\n", ii,
                                  referenceForces[ii][0], forceVectorToCompare[ii][0],
                                  referenceForces[ii][1], forceVectorToCompare[ii][1],
                                  referenceForces[ii][2], forceVectorToCompare[ii][2] );
               }
            }

         } else { 

            for( unsigned int ii = 0; ii < referenceForces.size() && ii < maxPrint; ii++ ){
               (void) fprintf( log, "%6u r[%14.7e %14.7e %14.7e] c[%14.7e %14.7e %14.7e]\n", ii,
                               referenceForces[ii][0], referenceForces[ii][1], referenceForces[ii][2],
                               forceVectorToCompare[ii][0], forceVectorToCompare[ii][1], forceVectorToCompare[ii][2] );
            }
            if( referenceForces.size() > maxPrint ){
               for( unsigned int ii = referenceForces.size() - maxPrint; ii < referenceForces.size(); ii++ ){
                  (void) fprintf( log, "%6u r[%14.7e %14.7e %14.7e] c[%14.7e %14.7e %14.7e]\n", ii,
                                  referenceForces[ii][0], referenceForces[ii][1], referenceForces[ii][2],
                                  forceVectorToCompare[ii][0], forceVectorToCompare[ii][1], forceVectorToCompare[ii][2] );
               }
            }
         }
      
      }
#endif

      // compare reference vs cuda forces

#if 0
      double maxDeltaRefCud                          = -1.0e+30;
      double maxRelativeDeltaRefCud                  = -1.0e+30;
      double maxDotRefCud                            = -1.0e+30;
      double maxDeltaPrmCud                          = -1.0e+30;
      double maxRelativeDeltaPrmCud                  = -1.0e+30;
      double maxDotPrmCud                            = -1.0e+30;

      std::vector<double> forceArray1Sum;
      std::vector<double> forceArray2Sum;
      std::vector<double> forceArray3Sum;

      std::vector<double> referenceForceStats;
      std::vector<double> cudaForceStats;
      std::vector<double> cudaForceStats1;
      std::vector<double> paramForceStats;

      compareForces( referenceForces, "fRef", forceArray1Sum, referenceForceStats,
                     forceVectorToCompare,      "fCud", forceArray2Sum, cudaForceStats, 
                     &maxDeltaRefCud, &maxRelativeDeltaRefCud, &maxDotRefCud, forceTolerance, log );
      
      (void) fflush( log );

      // summary file info

      if( summaryFile ){

         StringVector forceStringArray;
         System system = referenceContext->getSystem();
         getForceStrings( system, forceStringArray, log );
         std::string forceString;
         if( forceStringArray.size() > 5 ){
            forceString = "All";
         } else {
            for( StringVectorCI ii = forceStringArray.begin(); ii != forceStringArray.end(); ii++ ){
               forceString += *ii;
            }
         }
         if( forceString.size() < 1 ){
            forceString = "NA";
         }
         (void) fprintf( summaryFile, "Force %s\nAtoms %u\nMaxDelta %14.7e\nMaxRelDelta %14.7e\nMaxDot %14.7e\n",
                         forceString.c_str(), referenceForces.size(), maxDeltaRefCud, maxRelativeDeltaRefCud, maxDotRefCud);

         double sum = ( fabs(forceArray1Sum[0] ) + fabs( forceArray1Sum[1] ) + fabs( forceArray1Sum[2]) )*0.33333;
         (void) fprintf( summaryFile, "SumRef %14.7e\n", sum );

                sum = ( fabs(forceArray2Sum[0] ) + fabs( forceArray2Sum[1] ) + fabs( forceArray2Sum[2]) )*0.33333;
         (void) fprintf( summaryFile, "SumCuda %14.7e\n", sum );
         double difference         = fabs( referencePotentialEnergy - cudaPotentialEnergy );
         double relativeDifference = difference/( fabs( referencePotentialEnergy ) + fabs(  cudaPotentialEnergy ) + 1.0e-10);
         (void) fprintf( summaryFile, "RefPE %14.7e\nCudaPE %14.7e\nDiffPE %14.7e\nRelDiffPE %14.7e\n",
                         referencePotentialEnergy, cudaPotentialEnergy, difference, relativeDifference );
      }

      if( log ){
         (void) fprintf( log, "max delta=%14.7e maxRelDelta=%14.7e maxDot=%14.7e\n", maxDeltaRefCud, maxRelativeDeltaRefCud, maxDotRefCud);
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
         (void) fflush( log );
      }

      // check that relative force difference is small

      if( applyAssertion ){
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
         _synchContexts( *cudaContext, *referenceContext );
      }

   }

   if( log ){
      if( applyAssertion ){
         (void) fprintf( log, "\n%s tests passed\n", methodName.c_str() );
      } else {
         (void) fprintf( log, "\n%s tests off\n", methodName.c_str() );
      }
      (void) fflush( log );
   }
#endif

}

void testEnergyForcesConsistent( std::string parameterFileName, MapStringInt& forceMap, MapStringString& inputArgumentMap,
                                 FILE* inputLog, FILE* summaryFilePtr ){

// ---------------------------------------------------------------------------------------

  static const std::string methodName      = "testEnergyForcesConsistent";
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

   // get platform to test (1=cuda, 2=reference)

   std::string platformName;
   int platformInclude = -1;
   if( setStringFromMap(    inputArgumentMap, "platform",     platformName ) == 0 ){
      if( log ){
         (void) fprintf( log, "%s platform not set -- aborting.\n", methodName.c_str() );
         (void) fflush( log );
      }   
      return;
   } else if( platformName.compare( "Cuda" ) == 0 ){
      platformInclude = 1;
      if( log ){
         (void) fprintf( log, "%s Using Cuda platform.\n", methodName.c_str() );
      }
   } else if( platformName.compare( "Reference" ) == 0 ){
      platformInclude = 2;
      if( log ){
         (void) fprintf( log, "%s Using Reference platform.\n", methodName.c_str() );
      }
   } else {
      if( log ){
         (void) fprintf( log, "%s platform name not recognized: %s (valid names are Cuda & Reference).\n", methodName.c_str(), platformName.c_str() );
         (void) fflush( log );
      }   
      return;
   }

   double parameterKineticEnergy, parameterPotentialEnergy;

   std::vector<Vec3> parameterForces;
   std::vector<Vec3> parameterForces2;
   MapStringVectorOfVectors supplementary;

   if( platformInclude == 1 ){

      if( log ){
         (void) fprintf( log, "%s Testing cuda platform\n", methodName.c_str() );
         (void) fflush( log );
      }   

      Platform& cudaPlatform              = Platform::getPlatformByName( "Cuda");;

      Context* cudaContext                = testSetup( parameterFileName, forceMap,  cudaPlatform,
                                                         parameterForces2, &parameterKineticEnergy, &parameterPotentialEnergy,
                                                         supplementary, inputArgumentMap, log );

      checkEnergyForceConsistent( *cudaContext, inputArgumentMap, log, summaryFilePtr );

   } else {

      Platform& referencePlatform         =  Platform::getPlatformByName( "Reference");

      if( log ){
         (void) fprintf( log, "%s Testing reference platform\n", methodName.c_str() );
         (void) fflush( log );
      }   

      Context* referenceContext       = testSetup( parameterFileName, forceMap, referencePlatform, 
                                                   parameterForces, &parameterKineticEnergy, &parameterPotentialEnergy,
                                                   supplementary, inputArgumentMap, log );
      checkEnergyForceConsistent( *referenceContext, inputArgumentMap, log, summaryFilePtr );
   }

   return;
}

void testEnergyConservation( std::string parameterFileName, MapStringInt& forceMap, 
                             MapStringString& inputArgumentMap, FILE* inputLog, FILE* summaryFile ){

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
      (void) fprintf( log, "%s\n", methodName.c_str() );
      (void) fflush( log );
   }   

   Platform& referencePlatform         =  Platform::getPlatformByName( "Reference");

   double parameterKineticEnergy, parameterPotentialEnergy;

   std::vector<Vec3> parameterForces;
   std::vector<Vec3> parameterForces2;
   MapStringVectorOfVectors supplementary;

   Context* referenceContext       = testSetup( parameterFileName, forceMap,  referencePlatform,
                                                parameterForces2, &parameterKineticEnergy, &parameterPotentialEnergy,
                                                supplementary, inputArgumentMap, log );

   if( log ){
      (void) fprintf( log, "%s Testing cuda platform\n", methodName.c_str() );
      (void) fflush( log );
   }   

   checkEnergyConservation( *referenceContext, inputArgumentMap, log, summaryFile );
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
   (void) printf( "   -summaryFileName <summary file name> (default=no summary)\n" );
   (void) printf( "   -applyAssertion if set, apply assertion (default=apply)\n" );

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
   (void) printf( "   -checkInputForces check that cuda/reference forces agree w/ input forces\n" );
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

int getForceOffset( int argIndex, int maxArgs, char* inputArgument[], MapStringInt& forceMap ){

   // skip over '-'

   char* argument = inputArgument[argIndex];
   argument++;

   MapStringIntI forcePresent = forceMap.find( argument );
   if( forcePresent != forceMap.end() && argIndex < maxArgs ){
      (*forcePresent).second = atoi( inputArgument[argIndex+1] );
      return 1;
   }

   return 0;

}

/**---------------------------------------------------------------------------------------
 * Initialize forceMap
 *
 * @param forceMap        has w/ force name as key and int as value
 * @param initialValue    initial value
 *
 *
   --------------------------------------------------------------------------------------- */

void initializeForceMap( MapStringInt& forceMap, int initialValue ){

   forceMap[HARMONIC_BOND_FORCE]            = initialValue;
   forceMap[HARMONIC_ANGLE_FORCE]           = initialValue;
   forceMap[PERIODIC_TORSION_FORCE]         = initialValue;
   forceMap[RB_TORSION_FORCE]               = initialValue;
   forceMap[NB_FORCE]                       = initialValue;
   forceMap[NB_SOFTCORE_FORCE]              = initialValue;
   forceMap[NB_EXCEPTION_FORCE]             = initialValue;
   forceMap[NB_EXCEPTION_SOFTCORE_FORCE]    = initialValue;
   forceMap[GBSA_OBC_FORCE]                 = initialValue;
   forceMap[GBSA_OBC_SOFTCORE_FORCE]        = initialValue;
   forceMap[GBVI_FORCE]                     = initialValue;
   forceMap[GBVI_SOFTCORE_FORCE]            = initialValue;

   return;

}

// ---------------------------------------------------------------------------------------

int main( int numberOfArguments, char* argv[] ){

// ---------------------------------------------------------------------------------------

   static const std::string methodName               = "TestCudaFromFile";
   int checkForces                                   = 0;
   int checkEnergyForceConsistent                    = 0;
   int checkEnergyConservation                       = 0;
   int checkInputForces                              = 0;
   MapStringString inputArgumentMap;

   FILE* log                                         = NULL;
   FILE* summaryFile                                 = NULL;

// ---------------------------------------------------------------------------------------

   std::string defaultParameterFileName     = "OpenParameters.txt";

   if( numberOfArguments < 2 ){
      printUsage( defaultParameterFileName );
   }

   std::string parameterFileName            = defaultParameterFileName;
   std::string pluginDirectoryName          = Platform::getDefaultPluginsDirectory();
   MapStringInt forceMap;
   initializeForceMap( forceMap, 0 );
   int logFileNameIndex                     = -1;
   int summaryFileNameIndex                 = -1;

   // parse arguments

#ifdef _MSC_VER
#define STRCASECMP(X,Y)  stricmp(X,Y)
#define STRNCASECMP(X,Y,Z)  strnicmp(X,Y,Z)
#else
#define STRCASECMP(X,Y)  strcasecmp(X,Y)
#define STRNCASECMP(X,Y,Z)  strncasecmp(X,Y,Z)
#endif

   for( int ii = 1; ii < numberOfArguments; ii++ ){
      int addToMap = 0;
      if( STRCASECMP( argv[ii], "-parameterFileName" ) == 0 ){
         parameterFileName          = argv[ii+1];
         ii++;
      } else if( STRCASECMP( argv[ii], "-logFileName" ) == 0 ){
         logFileNameIndex           = ii + 1;
         ii++;
      } else if( STRCASECMP( argv[ii], "-summaryFileName" ) == 0 ){ 
         summaryFileNameIndex       = ii + 1; 
         ii++;
      } else if( STRCASECMP( argv[ii], "-checkForces" ) == 0 ){
         checkForces                = atoi( argv[ii+1] );
         ii++;
      } else if( STRCASECMP( argv[ii], "-checkInputForces" ) == 0 ){
         checkInputForces           = atoi( argv[ii+1] );
         ii++;
      } else if( STRCASECMP( argv[ii], "-checkEnergyForceConsistent" ) == 0 ){
         checkEnergyForceConsistent = atoi( argv[ii+1] );
         ii++;
      } else if( STRCASECMP( argv[ii], "-checkEnergyConservation" ) == 0 ){
         checkEnergyConservation = atoi( argv[ii+1] );;
         ii++;
      } else if( STRCASECMP( argv[ii], "-energyForceDelta" )     == 0    ||
                 STRCASECMP( argv[ii], "-energyForceTolerance" ) == 0    ||
                 STRCASECMP( argv[ii], "-cudaDeviceId" )         == 0    ||
                 STRCASECMP( argv[ii], "-platform" )             == 0    ||
                 STRCASECMP( argv[ii], "-applyAssertion" )       == 0 ){
         addToMap                   = ii;
         ii++;

      } else if( STRCASECMP( argv[ii], "-allForces" ) == 0 ){
         int flag = atoi( argv[ii+1] );
         ii++;
         initializeForceMap( forceMap, flag );
      } else if( STRCASECMP( argv[ii], "-log" ) == 0 ){
         if( atoi( argv[ii+1] ) != 0 ){
            log = stderr;
         } else {
            log = NULL;
         }
         ii++;
      } else if( STRCASECMP( argv[ii], "-help" ) == 0 ){
         printUsage( defaultParameterFileName );

      } else if( getForceOffset( ii, numberOfArguments, argv, forceMap ) ){
         ii++;
      } else if( STRNCASECMP( argv[ii], "-equilibration", 14  ) == 0 ||
                 STRNCASECMP( argv[ii], "-simulation",    11  ) == 0 ||
                 STRNCASECMP( argv[ii], "-runId",          6  ) == 0 ||
                 STRNCASECMP( argv[ii], "-nonbonded",     10  ) == 0 ||
                 STRNCASECMP( argv[ii], "-readContext",   12  ) == 0 ){
         addToMap = ii;
         ii++;
      } else {
         (void) printf( "Argument=<%s> not recognized -- aborting\n", argv[ii] );
         exit(-1);
      }
      if( addToMap && ii >= 1 ){
         char* key = argv[addToMap];
         key++;
         inputArgumentMap[key] = argv[addToMap+1];
//         (void) printf( "ArgumentMap =<%s> <%s>\n", argv[addToMap], argv[addToMap+1] );
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

   // summary file

   if( summaryFileNameIndex > -1 ){
#ifdef _MSC_VER
         fopen_s( &summaryFile, argv[summaryFileNameIndex], "w" );
#else
         summaryFile = fopen( argv[summaryFileNameIndex], "w" );
#endif
   }

   // log info

   if( log ){
      (void) fprintf( log, "Input arguments:\n" );
      for( int ii = 1; ii < numberOfArguments-1; ii += 2 ){
         (void) fprintf( log, "      %3d %30s %15s\n", ii, argv[ii], argv[ii+1] );
      }
      (void) fprintf( log, "parameter file=<%s>\n", parameterFileName.c_str() );

      if( summaryFileNameIndex > -1 ){
         (void) fprintf( log, "summary file=<%s>\n", argv[summaryFileNameIndex] );
      } else {
         (void) fprintf( log, "no summary file\n" );
      }

      (void) fprintf( log, "pluginDirectoryName         %s\n", pluginDirectoryName.c_str() );
      (void) fprintf( log, "checkEnergyForceConsistent  %d\n", checkEnergyForceConsistent );
      (void) fprintf( log, "checkEnergyConservation     %d\n", checkEnergyConservation );
      (void) fprintf( log, "checkInputForces            %d\n", checkInputForces );

      (void) fprintf( log, "ForceMap: %u\n", forceMap.size() );
      for( MapStringIntCI ii = forceMap.begin(); ii != forceMap.end(); ii++ ){
         (void) fprintf( log, "   %20s %d\n", (*ii).first.c_str(), (*ii).second );
      }
      (void) fprintf( log, "Argument map: %u\n", inputArgumentMap.size() );
      for( MapStringStringCI ii = inputArgumentMap.begin(); ii != inputArgumentMap.end(); ii++ ){
         (void) fprintf( log, "Map %s %s\n", (*ii).first.c_str(), (*ii).second.c_str() );
      }
      (void) fflush( log );
   }

   Platform::loadPluginsFromDirectory(pluginDirectoryName);

   // check forces

   if( checkForces ){
      try {
         testReferenceCudaForces( parameterFileName, forceMap, inputArgumentMap, log, summaryFile );
       } catch( const exception& e ){
         (void) fprintf( stderr, "Exception checkForces %s %s\n", methodName.c_str(),  e.what() ); (void) fflush( stderr );
         return 1;
      }   
   }

   // compare w/ input forces

   if( checkInputForces ){
      try {
         testInputForces( parameterFileName, forceMap, inputArgumentMap, log, summaryFile );
       } catch( const exception& e ){
         (void) fprintf( stderr, "Exception testInputForces %s %s\n", methodName.c_str(),  e.what() ); (void) fflush( stderr );
         return 1;
      }   
   }

   // check energy/force consistent

   if( checkEnergyForceConsistent ){
      try {
         testEnergyForcesConsistent( parameterFileName, forceMap, inputArgumentMap, log, summaryFile );
       } catch( const exception& e ){
         (void) fprintf( stderr, "Exception checkEnergyForceConsistent %s %s\n", methodName.c_str(),  e.what() ); (void) fflush( stderr );
         return 1;
      }   
   }
   

   // check energy conservation or thermal stability

   if( checkEnergyConservation ){
      try {
         testEnergyConservation( parameterFileName, forceMap, inputArgumentMap, log, summaryFile );
       } catch( const exception& e ){
         (void) fprintf( stderr, "Exception checkEnergyConservation %s %s\n", methodName.c_str(),  e.what() ); (void) fflush( stderr );
         return 1;
      }   
   }

   if( log ){
      (void) fprintf( log, "\n%s done\n", methodName.c_str() ); (void) fflush( log );
   }

   if( summaryFile ){
      (void) fclose( summaryFile );
   }

   return 0;
}
