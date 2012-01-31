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

#include "openmm/internal/AssertionUtilities.h"
#include "CudaPlatform.h"
#include "../../../platforms/opencl/include/OpenCLPlatform.h"
#include "ReferencePlatform.h"

#include "openmm/Context.h"

#include "openmm/HarmonicBondForce.h"
#include "openmm/CustomBondForce.h"

#include "openmm/HarmonicAngleForce.h"
#include "openmm/CustomAngleForce.h"

#include "openmm/PeriodicTorsionForce.h"
#include "openmm/RBTorsionForce.h"
#include "openmm/CustomTorsionForce.h"

#include "openmm/GBSAOBCForce.h"
#include "openmm/GBVIForce.h"

#include "openmm/NonbondedForce.h"
#include "openmm/CustomNonbondedForce.h"

#include "openmm/CMMotionRemover.h"
#include "openmm/System.h"
#include "openmm/LangevinIntegrator.h"
#include "openmm/VariableLangevinIntegrator.h"
#include "openmm/VerletIntegrator.h"
#include "openmm/VariableVerletIntegrator.h"
#include "openmm/BrownianIntegrator.h"
#include "sfmt/SFMT.h"

// free-energy plugin includes
//#define	INCLUDE_FREE_ENERGY_PLUGIN
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
#include <algorithm>
#include <map>

#ifdef _MSC_VER
   #define isinf !_finite
   #define isnan _isnan
#endif

// max entries to print for default output
#define MAX_PRINT 5

// force names

std::string HARMONIC_BOND_FORCE             = "HarmonicBond";
std::string CUSTOM_BOND_FORCE               = "CustomBond";

std::string HARMONIC_ANGLE_FORCE            = "HarmonicAngle"; 
std::string CUSTOM_ANGLE_FORCE              = "CustomAngle"; 

std::string PERIODIC_TORSION_FORCE          = "PeriodicTorsion";
std::string RB_TORSION_FORCE                = "RbTorsion";
std::string CUSTOM_TORSION_FORCE            = "CustomTorsion";

std::string NB_FORCE                        = "Nb";
std::string CUSTOM_NB_FORCE                 = "CustomNb";
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

static void copyMap( MapStringInt& inputMap, MapStringInt& outputMap ){

// ---------------------------------------------------------------------------------------

   static const std::string methodName             = "setIntFromMap";

// ---------------------------------------------------------------------------------------

   for( MapStringIntI ii = inputMap.begin(); ii != inputMap.end(); ii++ ){
         outputMap[(*ii).first] = (*ii).second;
   }
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

static void editMap( MapStringInt& inputMap, MapStringInt& outputMap, int newValue ){

// ---------------------------------------------------------------------------------------

   static const std::string methodName             = "editMap";

// ---------------------------------------------------------------------------------------

   copyMap( inputMap, outputMap );
   for( MapStringIntI ii = inputMap.begin(); ii != inputMap.end(); ii++ ){
      if( (*ii).second ){
         outputMap[(*ii).first] = newValue;
      }
   }
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

static int setFloatFromMap( MapStringString& argumentMap, std::string fieldToCheck, float& fieldToSet ){

// ---------------------------------------------------------------------------------------

   static const std::string methodName             = "setFloatFromMap";

// ---------------------------------------------------------------------------------------

   MapStringStringCI check = argumentMap.find( fieldToCheck );
   if( check != argumentMap.end() ){
      fieldToSet = static_cast<float>(atof( (*check).second.c_str() )); 
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

static CustomBondForce* copyToCustomBondForce(const HarmonicBondForce* bondForce ) { 

// ---------------------------------------------------------------------------------------

    //static const std::string methodName      = "copyToCustomBondForce";

// ---------------------------------------------------------------------------------------

    // bond force

    CustomBondForce* customBondForce = new CustomBondForce("scale*k*(r-r0)^2");
    customBondForce->addPerBondParameter("r0");
    customBondForce->addPerBondParameter("k");
    customBondForce->addPerBondParameter("scale");
    //customBondForce->addGlobalParameter("scale", 0.5);

    std::vector<double> parameters(3);
    unsigned int numberOfBonds = static_cast<unsigned int>(bondForce->getNumBonds());
    for( unsigned int ii = 0; ii < bondForce->getNumBonds(); ii++ ){
        int particle1, particle2;
        double length, k;
        bondForce->getBondParameters( ii, particle1, particle2, length, k ); 
        parameters[0] = length;
        parameters[1] = k; 
        parameters[2] = 0.5f; 
        customBondForce->addBond( particle1, particle2, parameters);
    }    

    return customBondForce;
}

static CustomAngleForce* copyToCustomAngleForce(const HarmonicAngleForce* angleForce, FILE* log ) { 

// ---------------------------------------------------------------------------------------

    //static const std::string methodName      = "copyToCustomAngleForce";

// ---------------------------------------------------------------------------------------

    // angle force

    CustomAngleForce* customAngleForce = new CustomAngleForce("scale*k*(theta-theta0)^2");
    customAngleForce->addPerAngleParameter("theta0");
    customAngleForce->addPerAngleParameter("k");
    //customAngleForce->addPerAngleParameter("scale");
    customAngleForce->addGlobalParameter("scale", 0.5);

    std::vector<double> parameters(2);
int targetAtom = 466;
    unsigned int numberOfAngles = static_cast<unsigned int>(angleForce->getNumAngles());
    for( unsigned int ii = 0; ii < angleForce->getNumAngles(); ii++ ){
        int particle1, particle2, particle3;
        double length, k;
        angleForce->getAngleParameters( ii, particle1, particle2, particle3, length, k ); 
if( particle1 == targetAtom || particle2 == targetAtom || particle3 == targetAtom ){
   (void) fprintf( log, "CstmQ [%5d %5d %5d] %14.6e %14.6e\n", particle1, particle2, particle3, length, k );
}
        parameters[0] = length;
        parameters[1] = k; 
//        parameters[2] = 0.5f; 
        customAngleForce->addAngle( particle1, particle2, particle3, parameters);
    }    

    return customAngleForce;
}

static CustomTorsionForce* copyToCustomPeriodicTorsionForce(const PeriodicTorsionForce* periodicTorsionForce ) { 

// ---------------------------------------------------------------------------------------

    //static const std::string methodName      = "copyToCustomTorsionForce";

// ---------------------------------------------------------------------------------------

    // periodic torsion force

    CustomTorsionForce* customTorsionForce = new CustomTorsionForce("k*(1+cos(n*theta-theta0))");
    customTorsionForce->addPerTorsionParameter("theta0");
    customTorsionForce->addPerTorsionParameter("n");
    customTorsionForce->addPerTorsionParameter("k");

    std::vector<double> parameters(3);
    unsigned int numberOfTorsions = static_cast<unsigned int>(periodicTorsionForce->getNumTorsions());
    for( unsigned int ii = 0; ii < periodicTorsionForce->getNumTorsions(); ii++ ){
        int particle1, particle2, particle3, particle4, periodicity;
        double phase, k;
        periodicTorsionForce->getTorsionParameters( ii, particle1, particle2, particle3, particle4, periodicity, phase, k ); 
        parameters[0] = phase;
        parameters[1] = static_cast<double>(periodicity); 
        parameters[2] = k;
        customTorsionForce->addTorsion( particle1, particle2, particle3, particle4, parameters);
    }    

    return customTorsionForce;
}

static CustomTorsionForce* copyToCustomRbTorsionForce(const RBTorsionForce* rbTorsionForce ) { 

// ---------------------------------------------------------------------------------------

    //static const std::string methodName      = "copyToCustomTorsionForce";

// ---------------------------------------------------------------------------------------

    // RB torsion force

    CustomTorsionForce* customTorsionForce = new CustomTorsionForce("c0+cp*(c1+cp*(c2+cp*(c3+cp*(c4+cp*c5)))); cp=-cos(theta)");
    customTorsionForce->addPerTorsionParameter("c0");
    customTorsionForce->addPerTorsionParameter("c1");
    customTorsionForce->addPerTorsionParameter("c2");
    customTorsionForce->addPerTorsionParameter("c3");
    customTorsionForce->addPerTorsionParameter("c4");
    customTorsionForce->addPerTorsionParameter("c5");

    std::vector<double> parameters(6);
    unsigned int numberOfTorsions = static_cast<unsigned int>(rbTorsionForce->getNumTorsions());
    for( unsigned int ii = 0; ii < rbTorsionForce->getNumTorsions(); ii++ ){
        int particle1, particle2, particle3, particle4;
        double c0, c1, c2, c3, c4, c5;
        rbTorsionForce->getTorsionParameters( ii, particle1, particle2, particle3, particle4, c0, c1, c2, c3, c4, c5 ); 
        parameters[0] = c0;
        parameters[1] = c1;
        parameters[2] = c2;
        parameters[3] = c3;
        parameters[4] = c4;
        parameters[5] = c5;
        customTorsionForce->addTorsion( particle1, particle2, particle3, particle4, parameters);
    }    

    return customTorsionForce;
}

static CustomNonbondedForce* copyToCustomNonbondedForce(const NonbondedForce* nonbondedForce ) { 

// ---------------------------------------------------------------------------------------

    //static const std::string methodName      = "copyToCustomTorsionForce";

// ---------------------------------------------------------------------------------------

    // nonbonded force

    CustomNonbondedForce* customNonbondedForce = new CustomNonbondedForce("4*eps*((sigma/r)^12-(sigma/r)^6)+138.935456*q/r; q=q1*q2; sigma=0.5*(sigma1+sigma2); eps=sqrt(eps1*eps2)" );
    customNonbondedForce->addPerParticleParameter("q");
    customNonbondedForce->addPerParticleParameter("sigma");
    customNonbondedForce->addPerParticleParameter("eps");

    std::vector<double> parameters(3);
    unsigned int numberOfNonbondeds = static_cast<unsigned int>(nonbondedForce->getNumParticles());
    for( unsigned int ii = 0; ii < nonbondedForce->getNumParticles(); ii++ ){
        double charge, sigma, epsilon;
        nonbondedForce->getParticleParameters(ii, charge, sigma, epsilon);
        parameters[0] = charge;
        parameters[1] = sigma;
        parameters[2] = epsilon;
        customNonbondedForce->addParticle(parameters);
    }    
    customNonbondedForce->setNonbondedMethod(CustomNonbondedForce::NoCutoff);

    // exceptions

    for( unsigned int ii = 0; ii < nonbondedForce->getNumExceptions(); ii++ ){
        int particle1, particle2;
        double chargeProd, sigma, epsilon;
        nonbondedForce->getExceptionParameters( ii, particle1, particle2, chargeProd, sigma, epsilon );
        customNonbondedForce->addExclusion( particle1, particle2 );
   }

    return customNonbondedForce;
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

   HarmonicBondForce* bondForce       = new HarmonicBondForce();
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

   // add force?

   MapStringIntI forceActive          = forceMap.find( HARMONIC_BOND_FORCE );
   if( forceActive != forceMap.end() ){
      if( (*forceActive).second == 1 ){
         if( log ){
             (void) fprintf( log, "Harmonic bond force is being included.\n" );
         }
         system.addForce( bondForce );
         return bondForce->getNumBonds();
      } else if(  (*forceActive).second == 2 ){
         CustomBondForce* customBondForce = copyToCustomBondForce( bondForce );
         if( log ){
            (void) fprintf( log, "Custom bond force is being included.\n" );
         }
         system.addForce( customBondForce );
         delete bondForce;
         return customBondForce->getNumBonds();
      } else if( log ){
         (void) fprintf( log, "force flag=%d not recognized.\n", (*forceActive).second );
      }
   } else {
      delete bondForce;
      if( log ){
         (void) fprintf( log, "Harmonic bond force is not being included.\n" );
      }
   }

   return 0;
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

   HarmonicAngleForce* angleForce = new HarmonicAngleForce();

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
         angleForce->addAngle( particle1, particle2, particle3, angle, k );
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
      unsigned int arraySize               = static_cast<unsigned int>(angleForce->getNumAngles());
      (void) fprintf( log, "%s: sample of HarmonicAngleForce parameters\n", methodName.c_str() );
      for( unsigned int ii = 0; ii < arraySize; ii++ ){
         int particle1, particle2, particle3;
         double angle, k;
         angleForce->getAngleParameters( ii, particle1, particle2, particle3, angle, k ); 
         (void) fprintf( log, "%8d %8d %8d %8d %14.7e %14.7e\n", ii, particle1, particle2, particle3, angle, k);
         if( ii == maxPrint ){
            ii = arraySize - maxPrint;
            if( ii < maxPrint )ii = maxPrint;
         }
      }
   }

   // add force?

   MapStringIntI forceActive          = forceMap.find( HARMONIC_ANGLE_FORCE );
   if( forceActive != forceMap.end() ){
      if( (*forceActive).second == 1 ){

         if( log ){
             (void) fprintf( log, "Harmonic angle force is being included.\n" );
         }

         system.addForce( angleForce );
         return angleForce->getNumAngles();

      } else if(  (*forceActive).second == 2 ){

         CustomAngleForce* customAngleForce = copyToCustomAngleForce( angleForce, log );
         if( log ){
            (void) fprintf( log, "Custom angle force is being included.\n" );
         }

         system.addForce( customAngleForce );
         delete angleForce;
         return customAngleForce->getNumAngles();

      } else if( log ){
         (void) fprintf( log, "force flag=%d not recognized.\n", (*forceActive).second );
      }
   } else {
      delete angleForce;
      if( log ){
         (void) fprintf( log, "Harmonic angle force is not being included.\n" );
      }
   }

   return 0;
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

   PeriodicTorsionForce* torsionForce = new PeriodicTorsionForce();
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
         torsionForce->addTorsion( particle1, particle2, particle3, particle4, periodicity, phase, k );
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
      unsigned int arraySize               = static_cast<unsigned int>(torsionForce->getNumTorsions());
      (void) fprintf( log, "%s: sample of PeriodicTorsionForce parameters\n", methodName.c_str() );
      for( unsigned int ii = 0; ii < arraySize; ii++ ){
         int particle1, particle2, particle3, particle4, periodicity;
         double phase, k;
         torsionForce->getTorsionParameters( ii, particle1, particle2, particle3, particle4, periodicity, phase, k );
         (void) fprintf( log, "%8d %8d %8d %8d %8d %8d %14.7e %14.7e\n", ii, particle1, particle2, particle3, particle4, periodicity, phase, k );
         if( ii == maxPrint ){
            ii = arraySize - maxPrint;
            if( ii < maxPrint )ii = maxPrint;
         }
      }
   }

   // add force?

   MapStringIntI forceActive          = forceMap.find( PERIODIC_TORSION_FORCE );
   if( forceActive != forceMap.end() ){
      if( (*forceActive).second == 1 ){
         if( log ){
             (void) fprintf( log, "Periodic torsion force is being included.\n" );
         }
         system.addForce( torsionForce );
         return torsionForce->getNumTorsions();
      } else if(  (*forceActive).second == 2 ){
         CustomTorsionForce* customTorsionForce = copyToCustomPeriodicTorsionForce( torsionForce );
         if( log ){
            (void) fprintf( log, "Custom periodic torsion force is being included.\n" );
         }
         system.addForce( customTorsionForce );
         delete torsionForce;
         return customTorsionForce->getNumTorsions();
      } else if( log ){
         (void) fprintf( log, "force flag=%d not recognized.\n", (*forceActive).second );
      }
   } else {
      delete torsionForce;
      if( log ){
         (void) fprintf( log, "Periodic torsion force is not being included.\n" );
      }
   }

   return 0;
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

   RBTorsionForce* torsionForce       = new RBTorsionForce();
   int numberOfTorsions               = atoi( tokens[1].c_str() );
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
         torsionForce->addTorsion( particle1, particle2, particle3, particle4, c0, c1, c2, c3, c4, c5 );
} 

#else


         torsionForce->addTorsion( particle1, particle2, particle3, particle4, c0, c1, c2, c3, c4, c5 );
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
      unsigned int arraySize               = static_cast<unsigned int>(torsionForce->getNumTorsions());
      (void) fprintf( log, "%s: sample of RBTorsionForce parameters\n", methodName.c_str() );
      for( unsigned int ii = 0; ii < arraySize; ii++ ){
         int particle1, particle2, particle3, particle4;
         double c0, c1, c2, c3, c4, c5;
         torsionForce->getTorsionParameters( ii, particle1, particle2, particle3, particle4, c0, c1, c2, c3, c4, c5 );
         (void) fprintf( log, "%8d %8d %8d %8d %8d %14.7e %14.7e %14.7e %14.7e %14.7e %14.7e\n",
                         ii, particle1, particle2, particle3, particle4, c0, c1, c2, c3, c4, c5 );
         if( ii == maxPrint ){
            ii = arraySize - maxPrint;
            if( ii < maxPrint )ii = maxPrint;
         }
      }
   }

   // add force?

   MapStringIntI forceActive          = forceMap.find( RB_TORSION_FORCE );
   if( forceActive != forceMap.end() ){
      if( (*forceActive).second == 1 ){
         if( log ){
             (void) fprintf( log, "RB torsion force is being included.\n" );
         }
         system.addForce( torsionForce );
         return torsionForce->getNumTorsions();
      } else if(  (*forceActive).second == 2 ){
         CustomTorsionForce* customTorsionForce = copyToCustomRbTorsionForce( torsionForce );
         if( log ){
            (void) fprintf( log, "Custom RB torsion force is being included.\n" );
         }
         system.addForce( customTorsionForce );
         delete torsionForce;
         return customTorsionForce->getNumTorsions();
      } else if( log ){
         (void) fprintf( log, "force flag=%d not recognized.\n", (*forceActive).second );
      }
   } else {
      delete torsionForce;
      if( log ){
         (void) fprintf( log, "RB torsion force is not being included.\n" );
      }
   }

   return 0;
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

   // add force?

   if( forceActive != forceMap.end() ){
      if( (*forceActive).second == 1 ){
         if( log ){
             (void) fprintf( log, "nonbonded force is being included.\n" );
         }
         system.addForce( nonbondedForce );
         return nonbondedForce->getNumParticles();
      } else if(  (*forceActive).second == 2 ){
         CustomNonbondedForce* customNonbondedForce = copyToCustomNonbondedForce( nonbondedForce );
         if( log ){
            (void) fprintf( log, "Custom nonbonded force is being included.\n" );
         }
         system.addForce( customNonbondedForce );
         delete nonbondedForce;
         return customNonbondedForce->getNumParticles();
      } else if( log ){
         (void) fprintf( log, "force flag=%d not recognized.\n", (*forceActive).second );
      }
   } else {
      delete nonbondedForce;
      if( log ){
         (void) fprintf( log, "Nonbonded force is not being included.\n" );
      }
   }

   return 0;
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
      (void) fprintf( log, "%s %s\n", methodName.c_str(), inputParameterFile.c_str() );
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
      (void) sprintf( buffer, "%s Input parameter file=<%s> could not be opened -- aborting.\n", methodName.c_str(), inputParameterFile.c_str() );
      throwException(__FILE__, __LINE__, buffer );
      (void) fflush( stderr);
      exit(-1);
   } else if( log ){
      (void) fprintf( log, "%s Input parameter file=<%s> opened.\n", methodName.c_str(), inputParameterFile.c_str() );
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

   Register forces associated w/ Reference free energy platform

   @param referencePlatform             reference platform

   --------------------------------------------------------------------------------------- */

static void registerFreeEnergyMethodsReferencePlatform( ReferencePlatform& referencePlatform ){

   // ---------------------------------------------------------------------------------------

   //static const char* methodName  = "registerFreeEnergyMethodsReferencePlatform: ";

   // ---------------------------------------------------------------------------------------

#ifdef INCLUDE_FREE_ENERGY_PLUGIN
   ReferenceFreeEnergyKernelFactory* factory  = new ReferenceFreeEnergyKernelFactory();

   referencePlatform.registerKernelFactory(CalcNonbondedSoftcoreForceKernel::Name(), factory);
   referencePlatform.registerKernelFactory(CalcGBVISoftcoreForceKernel::Name(), factory);
   referencePlatform.registerKernelFactory(CalcGBSAOBCSoftcoreForceKernel::Name(), factory);
#endif

}

/**---------------------------------------------------------------------------------------

   Register forces associated w/ Cuda free energy platform

   @param cudaPlatform             cuda platform

   --------------------------------------------------------------------------------------- */

static void registerFreeEnergyMethodsCudaPlatform( CudaPlatform& cudaPlatform ){

   // ---------------------------------------------------------------------------------------

   //static const char* methodName  = "registerFreeEnergyMethodsCudaPlatform: ";

   // ---------------------------------------------------------------------------------------

#ifdef INCLUDE_FREE_ENERGY_PLUGIN
   CudaFreeEnergyKernelFactory* factory  = new CudaFreeEnergyKernelFactory();

   cudaPlatform.registerKernelFactory(CalcNonbondedSoftcoreForceKernel::Name(), factory);
   cudaPlatform.registerKernelFactory(CalcGBVISoftcoreForceKernel::Name(), factory);
   cudaPlatform.registerKernelFactory(CalcGBSAOBCSoftcoreForceKernel::Name(), factory);
#endif

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
    registerFreeEnergyMethodsReferencePlatform( referencePlatform );

    CudaPlatform gpuPlatform;
    registerFreeEnergyMethodsCudaPlatform( gpuPlatform );

    if( platformName.compare( "ReferencePlatform" ) == 0 ){
       context = new Context( *system, *inputIntegrator, referencePlatform );
    } else {
       gpuPlatform.setPropertyDefaultValue( "CudaDevice", deviceId );
       context = new Context( *system, *inputIntegrator, gpuPlatform );
       if( log ){
          (void) fprintf( log, "OpenMM Platform: %s\n", context->getPlatform().getName().c_str() ); (void) fflush( log );
          const vector<string>& properties = gpuPlatform.getPropertyNames();
          for (unsigned int i = 0; i < properties.size(); i++) {
              fprintf( log, "%s: %s\n", properties[i].c_str(), gpuPlatform.getPropertyValue(*context, properties[i]).c_str());
          }    
       }    
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

static void getForceStrings( System& system, StringVector& forceStringArray, FILE* log ){

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
    
        if( !hit ){
            try {
               CustomBondForce& harmonicBondForce = dynamic_cast<CustomBondForce&>(force);
               forceStringArray.push_back( CUSTOM_BOND_FORCE );
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

        // custom angle

        if( !hit ){
    
            try {
               CustomAngleForce& harmonicAngleForce = dynamic_cast<CustomAngleForce&>(force);
               forceStringArray.push_back( CUSTOM_ANGLE_FORCE );
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
    
        // CustomTorsionForce
    
        if( !hit ){
            try {
               CustomTorsionForce& customTorsionForce = dynamic_cast<CustomTorsionForce&>(force);
               forceStringArray.push_back( CUSTOM_TORSION_FORCE );
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

        // nonbonded custom

        if( !hit ){
            try {
               CustomNonbondedForce& nbForce = dynamic_cast<CustomNonbondedForce&>(force);
               forceStringArray.push_back( CUSTOM_NB_FORCE );
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

    return;
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

void compareForces( const std::vector<Vec3>& forceArray1, const std::string& f1Name, std::vector<double>& forceArray1Sum, std::vector<double>& forceArray1Stats,
                   const std::vector<Vec3>& forceArray2, const std::string& f2Name, std::vector<double>& forceArray2Sum, std::vector<double>& forceArray2Stats,
                   double *averageDelta, double* averageRelativeDelta, double* maxDelta, int* maxDeltaIndex,
                   double* maxRelativeDelta, int* maxRelativeDeltaIndex, double* maxDot, double forceTolerance, FILE* inputLog ){

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
   *averageDelta                       = 0.0;
   *averageRelativeDelta               = 0.0;
   double averageRelativeDeltaCount    = 0.0;

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
      *averageDelta         += delta;
      double dotProduct      = f1[0]*f2[0] + f1[1]*f2[1] + f1[2]*f2[2];
             dotProduct     /= (normF1*normF2);
             dotProduct      = 1.0 - dotProduct;

      
      double relativeDelta;
      if( normF1 > 0.0 || normF2 > 0.0 ){
          relativeDelta               = (delta*2.0)/(normF1+normF2);
          *averageRelativeDelta      += relativeDelta;
          averageRelativeDeltaCount  += 1.0;
      } else {
          relativeDelta               = 0.0;
      }

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

   if( forceArray2.size() ){
      *averageDelta          /= (double)( forceArray2.size() );
   }

   if( averageRelativeDeltaCount ){
       *averageRelativeDelta /= averageRelativeDeltaCount;
   }

   findStatsForDouble( forceArray1Norms, forceArray1Stats );
   findStatsForDouble( forceArray2Norms, forceArray2Stats );

   return;
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

static void checkForcesDuringSimulation( int currentStep, Context& cudaContext, Context& referenceContext, FILE* log ) {
    
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
   double averageDelta;
   double averageRelativeDelta;
   int maxDeltaIndex;
   int maxRelativeDeltaRefCudIndex;
   
   std::vector<double> forceArray1Sum;
   std::vector<double> forceArray2Sum;
   std::vector<double> forceArray3Sum;
   
   std::vector<double> referenceForceStats;
   std::vector<double> cudaForceStats;
   
   compareForces( referenceForces, "fRef", forceArray1Sum, referenceForceStats,
                  cudaForces,      "fCud", forceArray2Sum, cudaForceStats, 
                  &averageDelta, &averageRelativeDelta, &maxDeltaRefCud, &maxDeltaIndex, &maxRelativeDeltaRefCud,
                  &maxRelativeDeltaRefCudIndex, &maxDotRefCud, forceTolerance, log );

   (void) fprintf( log, "MaxDelta=%13.7e at %d MaxRelativeDelta=%13.7e at %d maxDotRefCud=%14.6e averageDelta=%13.7e\n",
                   maxDeltaRefCud, maxDeltaIndex, maxRelativeDeltaRefCud, maxRelativeDeltaRefCudIndex, maxDotRefCud, averageDelta );
   (void) fprintf( log, "Reference force average=%14.7e stddev=%14.7e min=%14.7e at %6.0f max=%14.7e at %6.0f\n",
                   referenceForceStats[0], referenceForceStats[1], referenceForceStats[2], referenceForceStats[3],
                   referenceForceStats[4], referenceForceStats[5] );

   (void) fprintf( log, "     Cuda force average=%14.7e stddev=%14.7e min=%14.7e at %6.0f max=%14.7e at %6.0f\n",
                   cudaForceStats[0], cudaForceStats[1], cudaForceStats[2], cudaForceStats[3],
                   cudaForceStats[4], cudaForceStats[5] );

   (void) fflush( log );

   return;

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

   Context*       simulationContext = equilibrationContext;
   Integrator* simulationIntegrator = integrator;

   // get simulation integrator & context

/*
   Integrator* simulationIntegrator = _getIntegrator( simulationIntegratorName, simulationTimeStep,
                                                      simulationFriction, simulationTemperature,
                                                      simulationShakeTolerance, simulationErrorTolerance,
                                                      simulationSeed, log );
*/

   if( log ){
      _printIntegratorInfo( simulationIntegrator, log );
   }

   //delete equilibrationContext;
   //Context*       simulationContext = _getContext( &system, equilibrationContext, simulationIntegrator, "CudaPlatform", "SimulationContext", deviceId, log );
   //Context*       simulationContext = _getContext( &system, equilibrationContext, simulationIntegrator, "ReferencePlatform", "SimulationContext", deviceId, log );
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
fclose( summaryFile );
exit(0);
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
  int numberOfSteps                        = 2; 
  int steps                                = 0; 
  int applyAssertion                       = 1;
  int custom1                              = 0;
  int custom2                              = 0;
  int platform1                            = 0;
  int platform2                            = 1;

// ---------------------------------------------------------------------------------------

   FILE* log;
   if( PrintOn == 0 && inputLog ){
      log = NULL;
   } else {
      log = inputLog;
   } 

   setIntFromMap(    inputArgumentMap, "applyAssertion",     applyAssertion    );
   setIntFromMap(    inputArgumentMap, "custom1",            custom1           );
   setIntFromMap(    inputArgumentMap, "custom2",            custom2           );
   setIntFromMap(    inputArgumentMap, "platform1",          platform1         );
   setIntFromMap(    inputArgumentMap, "platform2",          platform2         );

   if( log ){
      (void) fprintf( log, "%s force tolerance=%.3e energy tolerance=%.3e step=%d\n",
                      methodName.c_str(), forceTolerance, energyTolerance, numberOfSteps );
      (void) fflush( log );
   }   


   MapStringInt forceMap1;
   copyMap( forceMap, forceMap1 );
   if( custom1 ){
      editMap( forceMap, forceMap1, 2 );
   } else {
      copyMap( forceMap, forceMap1 );
   }

   MapStringInt forceMap2;
   copyMap( forceMap, forceMap2 );
   if( custom2 ){
      editMap( forceMap, forceMap2, 2 );
   } else {
      copyMap( forceMap, forceMap2 );
   }

   ReferencePlatform referencePlatform1;
   registerFreeEnergyMethodsReferencePlatform( referencePlatform1 );
   ReferencePlatform referencePlatform2;
   registerFreeEnergyMethodsReferencePlatform( referencePlatform2 );

   CudaPlatform cudaPlatform1;
   registerFreeEnergyMethodsCudaPlatform( cudaPlatform1 );
   CudaPlatform cudaPlatform2;
   registerFreeEnergyMethodsCudaPlatform( cudaPlatform2 );

   double parameterKineticEnergy, parameterPotentialEnergy;

   std::vector<Vec3> parameterForces;
   std::vector<Vec3> parameterForces2;
   MapStringVectorOfVectors supplementary;

   Context* referenceContext;
   if( platform1 == 0 ){
       referenceContext      = testSetup( parameterFileName, forceMap1, referencePlatform1, 
                                          parameterForces, &parameterKineticEnergy, &parameterPotentialEnergy,
                                          supplementary, inputArgumentMap, log );
   } else {
      referenceContext       = testSetup( parameterFileName, forceMap1, cudaPlatform1, 
                                          parameterForces, &parameterKineticEnergy, &parameterPotentialEnergy,
                                          supplementary, inputArgumentMap, log );
   }

   Context* cudaContext;
   if( platform2 == 1 ){
      cudaContext            = testSetup( parameterFileName, forceMap2,  cudaPlatform2,
                                          parameterForces2, &parameterKineticEnergy, &parameterPotentialEnergy,
                                          supplementary, inputArgumentMap, log );
   } else {
      cudaContext            = testSetup( parameterFileName, forceMap2,  referencePlatform2,
                                          parameterForces2, &parameterKineticEnergy, &parameterPotentialEnergy,
                                          supplementary, inputArgumentMap, log );
   }

   (void) fprintf( log, "Platform1: %s\n", referenceContext->getPlatform().getName().c_str() ); (void) fflush( log );
   (void) fprintf( log, "Platform2: %s\n", cudaContext->getPlatform().getName().c_str() ); (void) fflush( log );
   Integrator& referenceIntegrator = referenceContext->getIntegrator();
   Integrator& cudaIntegrator      = cudaContext->getIntegrator();

   // Run several steps and see if relative force difference is within tolerance
   
   for( int step = 0; step < numberOfSteps; step++ ){

      // pull info out of contexts

      int types                                       = State::Positions | State::Velocities | State::Forces | State::Energy;

      State cudaState                                 =      cudaContext->getState( types );
      State referenceState                            = referenceContext->getState( types );

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
      double averageDelta;
      double averageRelativeDelta;
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
                     &averageDelta, &averageRelativeDelta, &maxDeltaRefCud, &maxDeltaIndex, &maxRelativeDeltaRefCud,
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
         (void) fprintf( summaryFile, "Force %s\nAtoms %u\nMaxDelta %14.7e\nMaxRelDelta %14.7e\nMaxDot %14.7e\nAverageDelta %14.7e\nAverageRelativeDelta %14.7e\n",
                         forceString.c_str(), referenceForces.size(), maxDeltaRefCud, maxRelativeDeltaRefCud, maxDotRefCud, averageDelta, averageRelativeDelta);

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

void testForces( std::string parameterFileName, MapStringInt& forceMap,
                 MapStringString& inputArgumentMap, FILE* inputLog, FILE* summaryFile ){

// ---------------------------------------------------------------------------------------

   static const std::string methodName      = "testForces";
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
 
    int customs[2] = { 0 , 0 };
    setIntFromMap(    inputArgumentMap, "applyAssertion",     applyAssertion    );
    setIntFromMap(    inputArgumentMap, "custom1",            customs[0]        );
    setIntFromMap(    inputArgumentMap, "custom2",            customs[1]        );
 
    std::string comparisonPlatform[2];
    if( setStringFromMap(    inputArgumentMap, "comparisonPlatform1", comparisonPlatform[0] ) == 0 ){
        comparisonPlatform[0] = "ReferencePlatform";
    }
    if( setStringFromMap(    inputArgumentMap, "comparisonPlatform2", comparisonPlatform[1] ) == 0 ){
        comparisonPlatform[1] = "ReferencePlatform";
    }
 
    MapStringInt forceMaps[2];
    for( int ii = 0; ii < 2; ii++ ){
         copyMap( forceMap, forceMaps[ii] );
         if( customs[ii] ){
             editMap( forceMap, forceMaps[ii], 2 );
         }
    }
 
    ReferencePlatform referencePlatform;
    CudaPlatform cudaPlatform;
    OpenCLPlatform openCLPlatform;
 
    double parameterKineticEnergy, parameterPotentialEnergy;
 
    std::vector<Vec3> parameterForces;
    std::vector<Vec3> parameterForces2;
    MapStringVectorOfVectors supplementary;
 
    Context* contexts[2];
    for( int ii = 0; ii < 2; ii++ ){
        if( comparisonPlatform[ii] == "ReferencePlatform" ){
            contexts[ii] = testSetup( parameterFileName, forceMaps[ii], referencePlatform, 
                                      parameterForces, &parameterKineticEnergy, &parameterPotentialEnergy,
                                      supplementary, inputArgumentMap, log );
        } else if( comparisonPlatform[ii] == "CudaPlatform" ){
            contexts[ii] = testSetup( parameterFileName, forceMaps[ii], cudaPlatform, 
                                      parameterForces, &parameterKineticEnergy, &parameterPotentialEnergy,
                                      supplementary, inputArgumentMap, log );
        } else if( comparisonPlatform[ii] == "OpenCLPlatform" ){
            contexts[ii] = testSetup( parameterFileName, forceMaps[ii], openCLPlatform, 
                                      parameterForces, &parameterKineticEnergy, &parameterPotentialEnergy,
                                      supplementary, inputArgumentMap, log );
        } else {
            (void) fprintf( log, "%s Platform: %d %s not recognized\n", methodName.c_str(), ii, comparisonPlatform[ii].c_str() );
            exit(0);
        }
        (void) fprintf( log, "Platform %d: %s\n", ii, contexts[ii]->getPlatform().getName().c_str() ); (void) fflush( log );
    }
 
    if( log ){
        (void) fprintf( log, "%s force tolerance=%.3e energy tolerance=%.3e step=%d\n",
                        methodName.c_str(), forceTolerance, energyTolerance, numberOfSteps );
        (void) fprintf( log, "Platform 0: %s\n", contexts[0]->getPlatform().getName().c_str() );
        (void) fprintf( log, "Platform 1: %s\n", contexts[1]->getPlatform().getName().c_str() );
        (void) fflush( log );
    }   
 
    // Run several steps and see if relative force difference is within tolerance
    
    double kineticEnergy[2];
    double potentialEnergy[2];
    std::vector<Vec3> forces[2];
    for( int step = 0; step < numberOfSteps; step++ ){
 
        // pull info out of contexts
  
        int types                                       = State::Positions | State::Velocities | State::Forces | State::Energy;
  
        std::vector<Vec3> coordinates[2];
        std::vector<Vec3> velocities[2];
        for( int ii = 0; ii < 2; ii++ ){
            State state             = contexts[ii]->getState( types );
            coordinates[ii]         = state.getPositions();
            velocities[ii]          = state.getVelocities();
            forces[ii]              = state.getForces();
            kineticEnergy[ii]       = state.getKineticEnergy();
            potentialEnergy[ii]     = state.getPotentialEnergy();
        }
  
        // diagnostics
  
        if( log ){
            //static const unsigned int maxPrint = MAX_PRINT;
            static const unsigned int maxPrint   = 1000000;
   
            (void) fprintf( log, "%s\n", methodName.c_str() );
            (void) fprintf( log, "Kinetic   energies: %14.7e %14.7e\n", kineticEnergy[0], kineticEnergy[1] );
            (void) fprintf( log, "Potential energies: %14.7e %14.7e\n", potentialEnergy[0], potentialEnergy[1] );
   
            for( unsigned int ii = 0; ii < forces[0].size(); ii++ ){
                (void) fprintf( log, "%6u [%14.7e %14.7e %14.7e] [%14.7e %14.7e %14.7e]\n", ii,
                                forces[0][ii][0], forces[0][ii][1], forces[0][ii][2],
                                forces[1][ii][0], forces[1][ii][1], forces[1][ii][2] );
                if( ii == maxPrint ){
                    ii = forces[0].size() - maxPrint;
                    if( ii < maxPrint )ii = maxPrint;
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
        double averageDelta;
        double averageRelativeDelta;
        int maxDeltaIndex;
        int maxRelativeDeltaRefCudIndex;
  
        std::vector<double> forceArraySum[2];
        std::vector<double> forceStats[2];

        compareForces( forces[0], comparisonPlatform[0], forceArraySum[0], forceStats[0],
                       forces[1], comparisonPlatform[1], forceArraySum[1], forceStats[1], 
                       &averageDelta, &averageRelativeDelta, &maxDeltaRefCud, &maxDeltaIndex, &maxRelativeDeltaRefCud,
                       &maxRelativeDeltaRefCudIndex, &maxDotRefCud, forceTolerance, log );
        
        (void) fflush( log );
 
        // summary file info
  
        if( summaryFile ){
   
            StringVector forceStringArray;
            System system = contexts[0]->getSystem();
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
            (void) fprintf( summaryFile,"Force %s\nPlatform1 %s\nPlatform2 %s\nAtoms %u\nMaxDelta %14.7e\nMaxRelDelta %14.7e\nMaxDot %14.7e\nAverageDelta %14.7e\nAverageRelativeDelta %14.7e\n",
                            forceString.c_str(),
                            contexts[0]->getPlatform().getName().c_str(),
                            contexts[1]->getPlatform().getName().c_str(),
                            forces[0].size(),
                            maxDeltaRefCud, maxRelativeDeltaRefCud, maxDotRefCud,
                            averageDelta, averageRelativeDelta);
   
            double sum = ( fabs(forceArraySum[0][0] ) + fabs( forceArraySum[0][1] ) + fabs( forceArraySum[0][2]) )*0.33333;
            (void) fprintf( summaryFile, "Sum1 %14.7e\n", sum );
   
                   sum = ( fabs(forceArraySum[1][0] ) + fabs( forceArraySum[1][1] ) + fabs( forceArraySum[1][2]) )*0.33333;
            (void) fprintf( summaryFile, "Sum2 %14.7e\n", sum );
            double difference         = fabs( potentialEnergy[0] - potentialEnergy[1] );
            double relativeDifference = difference/( fabs( potentialEnergy[0] ) + fabs(  potentialEnergy[1] ) + 1.0e-10);
            (void) fprintf( summaryFile, "PE1 %14.7e\nPE2 %14.7e\nDiffPE %14.7e\nRelDiffPE %14.7e\n",
                            potentialEnergy[0], potentialEnergy[1], difference, relativeDifference );
        }
 
        if( log ){
            (void) fprintf( log, "max delta=%13.7e at %d maxRelDelta=%13.7e at %d maxDot=%14.7e\n",
                            maxDeltaRefCud, maxDeltaIndex, maxRelativeDeltaRefCud, maxRelativeDeltaRefCudIndex, maxDotRefCud );
            for( int ii = 0; ii < 2; ii++ ){
                (void) fprintf( log, "%25s force sum [%14.7e %14.7e %14.7e]\n",
                                comparisonPlatform[ii].c_str(),
                                forceArraySum[ii][0], forceArraySum[ii][1], forceArraySum[ii][2] );
            }
            for( int ii = 0; ii < 2; ii++ ){
                (void) fprintf( log, "%25s force average=%14.7e stddev=%14.7e min=%14.7e at %6.0f max=%14.7e at %6.0f\n",
                                comparisonPlatform[ii].c_str(), 
                                forceStats[ii][0], forceStats[ii][1], forceStats[ii][2], forceStats[ii][3], forceStats[ii][4], forceStats[ii][5] ); 
            }
            (void) fflush( log );
        }
 
       // check that relative force difference is small
 
       if( applyAssertion ){
          ASSERT( maxRelativeDeltaRefCud < forceTolerance );
 
          // check energies
 
          ASSERT_EQUAL_TOL( kineticEnergy[0],    kineticEnergy[1],   energyTolerance );
          ASSERT_EQUAL_TOL( potentialEnergy[0],  potentialEnergy[1], energyTolerance );
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
          contexts[0]->getIntegrator().step( steps );
          _synchContexts( *contexts[0], *contexts[1]);
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

      CudaPlatform cudaPlatform;

      if( log ){
         (void) fprintf( log, "%s Testing cuda platform\n", methodName.c_str() );
         (void) fflush( log );
      }   

      registerFreeEnergyMethodsCudaPlatform( cudaPlatform );

      Context* cudaContext                  = testSetup( parameterFileName, forceMap,  cudaPlatform,
                                                         parameterForces2, &parameterKineticEnergy, &parameterPotentialEnergy,
                                                         supplementary, inputArgumentMap, log );

      checkEnergyForceConsistent( *cudaContext, inputArgumentMap, log, summaryFilePtr );

   } else {

      ReferencePlatform referencePlatform;
      registerFreeEnergyMethodsReferencePlatform( referencePlatform );

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

   //CudaPlatform cudaPlatform;
   ReferencePlatform referencePlatform;
   registerFreeEnergyMethodsReferencePlatform( referencePlatform );

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

// ---------------------------------------------------------------------------------------
// GB/VI test start
// ---------------------------------------------------------------------------------------

static Context* setupTwoParticle( FILE* log ){

    //ReferencePlatform platform;
    //CudaPlatform platform;
    const int numParticles         = 2;
    System* system                 = new System();
    LangevinIntegrator* integrator = new LangevinIntegrator(0, 0.1, 0.01);

    // harmonic bond

    double C_HBondDistance   = 0.05;
    HarmonicBondForce* bonds = new HarmonicBondForce();
    bonds->addBond(0, 1, C_HBondDistance, 0.0);
    system->addForce(bonds);

    double C_radius =  0.1;
    double C_gamma  =  0.0;
    double C_charge =  1.0;
    double H_radius =  1.00;
    double H_gamma  =  0.0;
    double H_charge = -1.0;

    NonbondedForce* nonbonded = new NonbondedForce();
    nonbonded->setNonbondedMethod(NonbondedForce::NoCutoff);

    (void) fprintf( log, "Applying GB/VI\n" );
    GBVIForce* forceField = new GBVIForce();
    for( int i = 0; i < numParticles; i++ ){
       system->addParticle(1.0);
       forceField->addParticle( H_charge, H_radius, H_gamma);
       nonbonded->addParticle(  H_charge, H_radius, 0.0);
    }
 
    forceField->setParticleParameters( 1, C_charge, C_radius, C_gamma);
//    nonbonded->setParticleParameters(  1, C_charge, C_radius, 0.0);
 
    forceField->addBond( 0, 1, C_HBondDistance );
    
    std::vector<pair<int, int> > bondExceptions;
    std::vector<double> bondDistances;
    
    bondExceptions.push_back(pair<int, int>(0, 1)); 
    bondDistances.push_back( C_HBondDistance );
    
    nonbonded->createExceptionsFromBonds(bondExceptions, 0.0, 0.0);
 
    system->addForce(forceField);
    system->addForce(nonbonded);

    //Context context(system, integrator, platform);
    Context* context = new Context( *system, *integrator);
    
    vector<Vec3> positions(numParticles);
    positions[0] = Vec3(0.0000,    0.0000,    0.0000);
    positions[1] = Vec3(0.500,    0.0000,    0.0000);
    context->setPositions(positions);

    State state = context->getState(State::Forces | State::Energy);
    (void) fprintf( log, "Energy %.4e\n", state.getPotentialEnergy() );

    return context;
} 

// ---------------------------------------------------------------------------------------
// GB/VI test start
// ---------------------------------------------------------------------------------------

static Context* setupEthane( FILE* log ){

    //ReferencePlatform platform;
    //CudaPlatform platform;
    const int numParticles         = 8;
    System* system                 = new System();
    LangevinIntegrator* integrator = new LangevinIntegrator(0, 0.1, 0.01);

    // harmonic bond

    double C_HBondDistance   = 0.1097;
    double C_CBondDistance   = 0.1504;
    HarmonicBondForce* bonds = new HarmonicBondForce();
    bonds->addBond(0, 1, C_HBondDistance, 0.0);
    bonds->addBond(2, 1, C_HBondDistance, 0.0);
    bonds->addBond(3, 1, C_HBondDistance, 0.0);

    bonds->addBond(1, 4, C_CBondDistance, 0.0);

    bonds->addBond(5, 4, C_HBondDistance, 0.0);
    bonds->addBond(6, 4, C_HBondDistance, 0.0);
    bonds->addBond(7, 4, C_HBondDistance, 0.0);

    system->addForce(bonds);

    double C_radius, C_gamma, C_charge, H_radius, H_gamma, H_charge;

    int AM1_BCC = 1;
    H_charge    = -0.053;
    C_charge    = -3.0*H_charge;
    if( AM1_BCC ){
       C_radius =  0.180;
       C_gamma  = -0.2863;
       H_radius =  0.125;
       H_gamma  =  0.2437;
    } else {
       C_radius =  0.215;
       C_gamma  = -1.1087;
       H_radius =  0.150;
       H_gamma  =  0.1237;
    }

    NonbondedForce* nonbonded = new NonbondedForce();
    nonbonded->setNonbondedMethod(NonbondedForce::NoCutoff);

    (void) fprintf( log, "Applying GB/VI\n" );
    GBVIForce* forceField = new GBVIForce();
    for( int i = 0; i < numParticles; i++ ){
       system->addParticle(1.0);
       forceField->addParticle( H_charge, H_radius, H_gamma);
       nonbonded->addParticle(  H_charge, H_radius, 0.0);
    }
 
    forceField->setParticleParameters( 1, C_charge, C_radius, C_gamma);
    forceField->setParticleParameters( 4, C_charge, C_radius, C_gamma);
 
    nonbonded->setParticleParameters(  1, C_charge, C_radius, 0.0);
    nonbonded->setParticleParameters(  4, C_charge, C_radius, 0.0);
 
    forceField->addBond( 0, 1, C_HBondDistance );
    forceField->addBond( 2, 1, C_HBondDistance );
    forceField->addBond( 3, 1, C_HBondDistance );
    forceField->addBond( 1, 4, C_CBondDistance );
    forceField->addBond( 5, 4, C_HBondDistance );
    forceField->addBond( 6, 4, C_HBondDistance );
    forceField->addBond( 7, 4, C_HBondDistance );
    
    std::vector<pair<int, int> > bondExceptions;
    std::vector<double> bondDistances;
    
    bondExceptions.push_back(pair<int, int>(0, 1)); 
    bondDistances.push_back( C_HBondDistance );
    
    bondExceptions.push_back(pair<int, int>(2, 1)); 
    bondDistances.push_back( C_HBondDistance );
    
    bondExceptions.push_back(pair<int, int>(3, 1)); 
    bondDistances.push_back( C_HBondDistance );
    
    bondExceptions.push_back(pair<int, int>(1, 4)); 
    bondDistances.push_back( C_CBondDistance );
    
    bondExceptions.push_back(pair<int, int>(5, 4)); 
    bondDistances.push_back( C_HBondDistance );
    
    bondExceptions.push_back(pair<int, int>(6, 4)); 
    bondDistances.push_back( C_HBondDistance );
 
    bondExceptions.push_back(pair<int, int>(7, 4));
    bondDistances.push_back( C_HBondDistance );
 
    nonbonded->createExceptionsFromBonds(bondExceptions, 0.0, 0.0);
 
    system->addForce(forceField);
    system->addForce(nonbonded);

    //Context context(system, integrator, platform);
    Context* context = new Context( *system, *integrator);
    
    vector<Vec3> positions(numParticles);
    positions[0] = Vec3(0.5480,    1.7661,    0.0000);
    positions[1] = Vec3(0.7286,    0.8978,    0.6468);
    positions[2] = Vec3(0.4974,    0.0000,    0.0588);
    positions[3] = Vec3(0.0000,    0.9459,    1.4666);
    positions[4] = Vec3(2.1421,    0.8746,    1.1615);
    positions[5] = Vec3(2.3239,    0.0050,    1.8065);
    positions[6] = Vec3(2.8705,    0.8295,    0.3416);
    positions[7] = Vec3(2.3722,    1.7711,    1.7518);
    for( int ii = 0; ii < positions.size(); ii++ ) {
        positions[ii][0] *= 0.1;
        positions[ii][1] *= 0.1;
        positions[ii][2] *= 0.1;
    }
    context->setPositions(positions);

    State state = context->getState(State::Forces | State::Energy);
    (void) fprintf( log, "Energy %.4e\n", state.getPotentialEnergy() );

    return context;
} 

// ---------------------------------------------------------------------------------------
// GB/VI small molecule input
// ---------------------------------------------------------------------------------------

static Context* setupSmallMolecule( FILE* filePtr, int* lineCount, FILE* log ){

    System* system                  = new System();
    LangevinIntegrator* integrator  = new LangevinIntegrator(0, 0.1, 0.01);

    NonbondedForce* nonbonded       = new NonbondedForce();
    nonbonded->setNonbondedMethod(NonbondedForce::NoCutoff);
    system->addForce(nonbonded);

    GBVIForce* forceField           = new GBVIForce();
    system->addForce(forceField);

    StringVector tokensName;
    char* isNotEof                  = readLine( filePtr, tokensName, lineCount, log );
    std::string moleculeName        = tokensName[0].c_str();

    StringVector tokensNumberParticles;
    isNotEof                        = readLine( filePtr, tokensNumberParticles, lineCount, log );
    const int numParticles          = atoi( tokensNumberParticles[0].c_str() );
    vector<Vec3> positions(numParticles);
 
    (void) fprintf( log, "setupSmallMolecule: %20s %5d\n", moleculeName.c_str(), numParticles );
    (void) fflush( log );

    for( unsigned int ii = 0; ii < numParticles; ii++ ){

        unsigned int tokenIndex     = 1;
        StringVector tokens;
        isNotEof                    = readLine( filePtr, tokens, lineCount, log );

       double mass                  = atof( tokens[tokenIndex++].c_str() );
       system->addParticle(mass);

       double charge                = atof( tokens[tokenIndex++].c_str() );
       double radius                = atof( tokens[tokenIndex++].c_str() );
#if 0
if( mass < 1.01 ){
    radius = .125; // H
} else if( 11.99 < mass && mass < 12.01 ){ 
    radius = 0.18; // C
} else if( 13.99 < mass && mass < 14.01 ){ 
    radius = 0.17; // N
} else if( 15.99 < mass && mass < 16.01 ){ 
    radius = 0.16; // O
} else if( 18.99 < mass && mass < 19.01 ){ 
    radius = 0.15; // F
} else if( 30.99 < mass && mass < 31.01 ){ 
    radius = 0.20; // P
} else if( 31.99 < mass && mass < 32.01 ){ 
    radius = 0.19; // S
} else if( 34.99 < mass && mass < 35.01 ){ 
    radius = 0.18; // CL
} else if( 78.99 < mass && mass < 79.01 ){ 
    radius = 0.22; // BR
} else if( 126.99 < mass && mass < 127.01 ){ 
    radius = 0.25; // I
} else {
    (void) fprintf( log, "Mass not handled %.4e   radius=-%14.6e\n", mass, radius );
}
#endif
       double gamma                 = atof( tokens[tokenIndex++].c_str() );

       double x1                    = atof( tokens[tokenIndex++].c_str() );
       double x2                    = atof( tokens[tokenIndex++].c_str() );
       double x3                    = atof( tokens[tokenIndex++].c_str() );

       forceField->addParticle( charge, radius, gamma);
       nonbonded->addParticle( charge, radius, 0.0);

       positions[ii]                = Vec3(x1,x2,x3);
    }
   
    StringVector tokensNumberBonds;
    isNotEof                        = readLine( filePtr, tokensNumberBonds, lineCount, log );
    const int numBonds              = atoi( tokensNumberBonds[0].c_str() );
    for( unsigned int ii = 0; ii < numBonds; ii++ ){

        unsigned int tokenIndex     = 1;
        StringVector tokens;
        isNotEof                    = readLine( filePtr, tokens, lineCount, log );

       int atom1                    = atoi( tokens[tokenIndex++].c_str() );
       int atom2                    = atoi( tokens[tokenIndex++].c_str() );
       double bondDistance          = atof( tokens[tokenIndex++].c_str() );

       forceField->addBond( atom1, atom2, bondDistance );
    }
    StringVector tokensJunk;
    isNotEof                    = readLine( filePtr, tokensJunk, lineCount, log );
   
    Context* context            = new Context( *system, *integrator);
    context->setPositions(positions);

    State state = context->getState(State::Forces | State::Energy);
    (void) fprintf( log, "Energy %.4e\n", state.getPotentialEnergy() );

    return context;
} 

static void getGBVIRadii( System& system, std::string forceName, std::vector<float>& radii, FILE* log ){

    for( int ii = 0; ii < system.getNumForces(); ii++ ) {

        Force& force            = system.getForce(ii);

        // GBVI
    
        if( forceName == GBVI_FORCE ){
            try {
                GBVIForce& gbviForce = dynamic_cast<GBVIForce&>(force);
    
                unsigned int numberOfParticles = static_cast<unsigned int>(gbviForce.getNumParticles());
                radii.resize( numberOfParticles );
                for( unsigned int ii = 0; ii < radii.size(); ii++ ){
                   double charge, radius, gamma;
                   gbviForce.getParticleParameters( ii, charge, radius, gamma );
                   radii[ii] = radius;
                }
                return;
            } catch( std::bad_cast ){
            }
        }
    
        // GBVI softcore
    
#ifdef INCLUDE_FREE_ENERGY_PLUGIN
        if( forceName == GBVI_SOFTCORE_FORCE ){
            try {
                GBVISoftcoreForce& gbviForce = dynamic_cast<GBVISoftcoreForce&>(force);
                unsigned int numberOfParticles = static_cast<unsigned int>(gbviForce.getNumParticles());
                radii.resize( numberOfParticles );
                for( unsigned int ii = 0; ii < radii.size(); ii++ ){
                   double charge, radius, gamma;
                   gbviForce.getParticleParameters( ii, charge, radius, gamma );
                   radii[ii] = radius;
                }
                return;
            } catch( std::bad_cast ){
            }
        }
#endif
    
    }

    // force not found

    radii.resize( 0 );

    return;
}

// get distance

static float getDistance( const float coordinates1[3], const float coordinates2[3]  ){

    float distance = (coordinates1[0] - coordinates2[0] )*(coordinates1[0] - coordinates2[0] ) +
                     (coordinates1[1] - coordinates2[1] )*(coordinates1[1] - coordinates2[1] ) +
                     (coordinates1[2] - coordinates2[2] )*(coordinates1[2] - coordinates2[2] );

    return sqrtf( distance );
}

static float getMaxRadius( const std::vector<float>& radii, FILE* log ){

    float maxRadius = -1.0e+30;
    for( int ii = 0; ii < radii.size(); ii++ ) {
       if( radii[ii] > maxRadius )maxRadius = radii[ii];
    }
    return maxRadius;
}

struct BoundingBox {
    float minCoordinates[3];
    float maxCoordinates[3];
};

// Get bounding box

void getBoundingBox( const std::vector<Vec3>& positions, struct BoundingBox& boundingBox,
                     const float maxAtomicRadius, FILE* log ){

// ---------------------------------------------------------------------------------------

    //static const std::string methodName      = "getBoundingBox";

// ---------------------------------------------------------------------------------------

    // initialize min/max bounding box  coordinates
    // and find min/max positions
    // widen box to enclose max atomic radius

    for( unsigned int ii = 0; ii < 3; ii++ ){
        boundingBox.minCoordinates[ii] =  1.0e+30;
        boundingBox.maxCoordinates[ii] = -1.0e+30;
    }

    for( unsigned int ii = 0; ii < positions.size(); ii++ ){
        for( unsigned int jj = 0; jj < 3; jj++ ){
            if( positions[ii][jj] < boundingBox.minCoordinates[jj] ){
                boundingBox.minCoordinates[jj] = positions[ii][jj];
            }
            if( positions[ii][jj] > boundingBox.maxCoordinates[jj] ){
                boundingBox.maxCoordinates[jj] = positions[ii][jj];
            }
        }
    }

    for( unsigned int ii = 0; ii < 3; ii++ ){
        boundingBox.minCoordinates[ii] -= maxAtomicRadius;
        boundingBox.maxCoordinates[ii] += maxAtomicRadius;
    }

    (void) fprintf( log, "BoundingBox min: [%12.5f %12.5f %12.5f]\n",       boundingBox.minCoordinates[0],  boundingBox.minCoordinates[1], boundingBox.minCoordinates[2]);
    (void) fprintf( log, "            max: [%12.5f %12.5f %12.5f]\n",       boundingBox.maxCoordinates[0],  boundingBox.maxCoordinates[1], boundingBox.maxCoordinates[2]);
    (void) fflush( log );

    return;
}

struct IntegrationGrid {
    float resolution;
    float origin[3];
    int numberOfGridPoints[3];
    int gridIndex[3];
    int offset[3];
    int totalNumberOfGridPoints;
    int totalNumberOfGridMaskPoints;
    int* mask;
    int* savedMask;
};

// print grid

static void printGrid( const struct IntegrationGrid& grid, FILE* log ){

// ---------------------------------------------------------------------------------------

    static const std::string methodName      = "printGrid";

// ---------------------------------------------------------------------------------------

    (void) fprintf( log, "Grid:\n" );
    (void) fprintf( log, "            : [%12d %12d %12d]\n",       grid.numberOfGridPoints[0],
                                                                   grid.numberOfGridPoints[1], grid.numberOfGridPoints[2]);
    (void) fprintf( log, "       total:  %12d %12d\n",             grid.totalNumberOfGridPoints, grid.totalNumberOfGridMaskPoints );
    (void) fprintf( log, "  resolution:  %12.5f\n",                grid.resolution );
    (void) fprintf( log, "      origin: [%12.5f %12.5f %12.5f]\n", grid.origin[0],  grid.origin[1], grid.origin[2]);
    (void) fprintf( log, "     indices: [%12d %12d %12d]\n",       grid.gridIndex[0],  grid.gridIndex[1], grid.gridIndex[2]);
    (void) fprintf( log, "      offset: [%12d %12d %12d]\n",       grid.offset[0],  grid.offset[1], grid.offset[2]);
    (void) fflush( log );

    return;
}

// Get grid

static void getGrid( const struct BoundingBox& boundingBox,  struct IntegrationGrid& grid, FILE* log ){

// ---------------------------------------------------------------------------------------

    static const std::string methodName      = "getGrid";

// ---------------------------------------------------------------------------------------

    if( grid.resolution <= 0.0f ){
        (void) fprintf( log, "%s grid resolution is invalid: %14.6e\n", methodName.c_str(), grid.resolution );
        (void) fflush( log );
        exit(-1);
    }

    // initialize min/max bounding box  coordinates
    // and find min/max positions

    grid.totalNumberOfGridPoints = 1;
    for( unsigned int ii = 0; ii < 3; ii++ ){
        float span                     = boundingBox.maxCoordinates[ii] - boundingBox.minCoordinates[ii];
		  grid.origin[ii]                = boundingBox.minCoordinates[ii];
		  grid.numberOfGridPoints[ii]    = static_cast<int>(span/grid.resolution) + 1;
        grid.totalNumberOfGridPoints  *= grid.numberOfGridPoints[ii];
        grid.offset[ii]                = (ii > 0) ? (grid.numberOfGridPoints[ii-1]*grid.offset[ii-1]) : 1;
    }
    grid.totalNumberOfGridMaskPoints = (grid.totalNumberOfGridPoints + 31)/32;
    if( grid.totalNumberOfGridMaskPoints ){
        grid.mask       = (int*) malloc( sizeof(int)*grid.totalNumberOfGridMaskPoints );
        grid.savedMask  = (int*) malloc( sizeof(int)*grid.totalNumberOfGridMaskPoints );
        memset( grid.mask, 0, sizeof(int)*grid.totalNumberOfGridMaskPoints );
    } else {
        grid.mask      = NULL;
        grid.savedMask = NULL;
    }

    return;
}

// get grid coordinates given grid index

static void getGridCoordinatesGivenIndex( int gridIndex, struct IntegrationGrid& grid,
                                          float gridCoordinates[3], int gridIndices[3], FILE* log ){

// ---------------------------------------------------------------------------------------

    static const std::string methodName      = "getGridCoordinatesGivenIndex";

// ---------------------------------------------------------------------------------------

    for( int ii = 2; ii >= 0; ii-- ){
//fprintf( log, "%s %u %d %d\n", methodName.c_str(), ii, gridIndex, grid.offset[ii] ); fflush( log );
        gridIndices[ii]    = gridIndex/grid.offset[ii]; 
        gridIndex     -= gridIndices[ii]*grid.offset[ii];
    }

    for( unsigned int ii = 0; ii < 3; ii++ ){
        gridCoordinates[ii] = gridIndices[ii]*grid.resolution + grid.origin[ii];
    } 

    return;
}

// get grid index given integer grid coordinates given grid index

static void getGridIndexGivenCoordinates( struct IntegrationGrid& grid,
                                          float gridCoordinates[3],
                                          int& gridIndex, int gridIndices[3], FILE* log ){

// ---------------------------------------------------------------------------------------

    // static const std::string methodName      = "getGridIndexGivenGridCoordinates";

// ---------------------------------------------------------------------------------------

    gridIndex = 0;
    for( unsigned int ii = 0; ii < 3; ii++ ){
        gridIndices[ii]   = static_cast<int>( (gridCoordinates[ii] - grid.origin[ii])/grid.resolution + 1.0e-04f);
        gridIndex        += gridIndices[ii]*grid.offset[ii];
    }

    return;
}

// set mask at given grid index

static std::string printMask( int maskValue, FILE* log ){

// ---------------------------------------------------------------------------------------

    //static const std::string methodName      = "maskGridPoint";

// ---------------------------------------------------------------------------------------

    std::string bitFields;
    for( int ii = 0; ii < 32; ii++ ){
       int mask = (1 << ii);
       if( mask & maskValue ){
           bitFields += "1";
       } else {
           bitFields += "0";
       }
    }

    return bitFields;
}

// set mask at given grid index

static int getMaskGridValue( int gridIndex, struct IntegrationGrid& grid, FILE* log ){

// ---------------------------------------------------------------------------------------

    //static const std::string methodName      = "maskGridPoint";

// ---------------------------------------------------------------------------------------

    int maskIndex = gridIndex/32; 
    if( maskIndex < grid.totalNumberOfGridMaskPoints ){
        return grid.mask[maskIndex];
    } else {
        (void) fprintf( log, "Grid index %d (maskIndex=%d) is out of range [0, %d] ([0, %d])\n", gridIndex, maskIndex,
                        grid.totalNumberOfGridPoints, grid.totalNumberOfGridMaskPoints );
        exit(-1);
    }

    return 0;
}


// set mask at given grid index

static void maskGridPoint( int gridIndex, struct IntegrationGrid& grid, int maskValue, FILE* log ){

// ---------------------------------------------------------------------------------------

    //static const std::string methodName      = "maskGridPoint";

// ---------------------------------------------------------------------------------------

    int maskIndex = gridIndex/32; 
    if( maskIndex < grid.totalNumberOfGridMaskPoints ){
        int offset = gridIndex - maskIndex*32;
        if( maskValue ){
            grid.mask[maskIndex] |=  (1 << offset);
        } else {
//std::string vvv = printMask( grid.mask[maskIndex], log );
            grid.mask[maskIndex] &= ~(1 << offset);
/*
std::string www = printMask( grid.mask[maskIndex], log );
(void) fprintf( log, "Grid index %8d (maskIndex=%8d) offset=%8d %s %s\n", gridIndex, maskIndex, offset, vvv.c_str(), www.c_str() );
*/

        }
    } else {
        (void) fprintf( log, "Grid index %d (maskIndex=%d) is out of range [0, %d] ([0, %d])\n", gridIndex, maskIndex,
                        grid.totalNumberOfGridPoints, grid.totalNumberOfGridMaskPoints );
    }

    return;
}

// diagnoistic

static int isGridIndexSet( const struct IntegrationGrid& grid,
                           const int gridIndex, FILE* log ){

// ---------------------------------------------------------------------------------------

   static const std::string methodName      = "isGridIndexSet";

// ---------------------------------------------------------------------------------------

    int maskIndex    = gridIndex/32;
    int offset       = gridIndex - maskIndex*32;
    if( grid.mask[maskIndex] & (1 << offset) ){
        return 1;
    } else {
        return 0;
    }
}

// diagnoistic

static void checkGridSetting( const float atomCoordinates[3], const struct IntegrationGrid& grid,
                              const int gridIndex, const std::string& idString, FILE* log ){

// ---------------------------------------------------------------------------------------

   static const std::string methodName      = "checkGridSetting";

// ---------------------------------------------------------------------------------------

    int maskIndex    = gridIndex/32;
    int offset       = gridIndex - maskIndex*32;
    std::string www  = printMask( grid.mask[maskIndex], log );
    (void) fprintf( log, "Check: %20s Grid index=%8d maskIndex=%8d offset=%2d %s ",
                    idString.c_str(), gridIndex, maskIndex, offset, www.c_str() );
    if( grid.mask[maskIndex] & (1 << offset) )fprintf( log, " hit offset" );
    fprintf( log, "\n" );
}

// get total masked count

static void getTotalMasked( struct IntegrationGrid& grid, unsigned int& totalMasked ){

// ---------------------------------------------------------------------------------------

    //static const std::string methodName      = "maskGridPoint";

// ---------------------------------------------------------------------------------------

    totalMasked  = 0;
    for( unsigned int ii = 0; ii < grid.totalNumberOfGridMaskPoints; ii++ ){
        if( grid.mask[ii] ){
            for( unsigned int jj = 0; jj < 32; jj++ ){
                if( grid.mask[ii] & (1 << jj) )totalMasked++;
            }
        }
    }

    return;
}

// used to sort atoms

typedef struct {
    int listIndex;
    unsigned int gridCount;
    unsigned int maskedCount;
    unsigned int bornMaskedCount;
    float radius;
    double bornSum;
    double bornRadius;
    double scaledRadius;
    int countCases[4];
    double minR;
    float coordinates[3];
    float range[2][3];
} sortedAtom;

static int compareAtomCoordinate( const sortedAtom& p1, const sortedAtom& p2){
    return p1.coordinates[0] < p2.coordinates[0];
}

static void sortAtomsByCoordinate( const std::vector<Vec3>& positions, 
                                   const std::vector<float>& radii,
                                   std::vector<sortedAtom>& sortedAtoms, FILE* log ){

// ---------------------------------------------------------------------------------------

  static const std::string methodName      = "sortAtomsByCoordinate";

// ---------------------------------------------------------------------------------------

    sortedAtoms.resize( positions.size() );
    for( unsigned int ii = 0; ii < positions.size(); ii++ ){
        sortedAtoms[ii].listIndex      = ii;
        sortedAtoms[ii].radius         = radii[ii];
        sortedAtoms[ii].coordinates[0] = positions[ii][0];
        sortedAtoms[ii].coordinates[1] = positions[ii][1];
        sortedAtoms[ii].coordinates[2] = positions[ii][2];
        sortedAtoms[ii].range[0][0]    = positions[ii][0] - radii[ii];
        sortedAtoms[ii].range[1][0]    = positions[ii][0] + radii[ii];
        sortedAtoms[ii].range[0][1]    = positions[ii][1] - radii[ii];
        sortedAtoms[ii].range[1][1]    = positions[ii][1] + radii[ii];
        sortedAtoms[ii].range[0][2]    = positions[ii][2] - radii[ii];
        sortedAtoms[ii].range[1][2]    = positions[ii][2] + radii[ii];
    }
    std::sort( sortedAtoms.begin(), sortedAtoms.end(), compareAtomCoordinate );
}

static int cubesOverlap( int dimension, const float cube1[2][3], const float cube2[2][3] ){

    if( cube1[0][dimension] >= cube2[0][dimension] &&
        cube1[0][dimension] <= cube2[1][dimension] )return 1;

    if( cube1[1][dimension] >= cube2[0][dimension] &&
        cube1[1][dimension] <= cube2[1][dimension] )return 1;

    return 0;
}

// get indices of atoms in cube

static void getAtomsInCube( float cube[2][3], const std::vector<sortedAtom>& sortedAtoms,
                            std::vector<unsigned int>& prospectiveAtoms, FILE* log ){

// ---------------------------------------------------------------------------------------

  static const std::string methodName      = "getAtomsInCube";

// ---------------------------------------------------------------------------------------

    prospectiveAtoms.resize(0);
    if( cube[1][0] < sortedAtoms[0].coordinates[0] ||
        cube[0][0] > sortedAtoms[sortedAtoms.size()-1].coordinates[0] ){
        return;
    }

    unsigned int lowerIndex    = 0;
    unsigned int upperIndex    = sortedAtoms.size();
    unsigned int midpoint      = (lowerIndex + upperIndex)/2;
    unsigned int done          = 0;
    while( done == 0 ){
        if( cube[0][0] > sortedAtoms[midpoint].coordinates[0] ){
            lowerIndex = midpoint; 
        } else if( cube[1][0] < sortedAtoms[midpoint].coordinates[0] ){
            upperIndex = midpoint; 
        } else {
            done = 1;
        }
        if( (upperIndex - lowerIndex) < 2 )done = 1;
        midpoint = done ? midpoint : (lowerIndex + upperIndex)/2;
    }
 
    // check if atom was found
/*
    if( cube[0][0] > sortedAtoms[midpoint].coordinates[0] ||
        cube[1][0] < sortedAtoms[midpoint].coordinates[0] ){
        if( 0 && log ){
            (void) fprintf( log, "No atom in cube [%14.6e %14.6e%14.6e] [%14.6e %14.6e%14.6e]\n",
                            cube[0][0], cube[0][1], cube[0][2],
                            cube[1][0], cube[1][1], cube[1][2] );
        }
        return;
    }
*/

    int atomIndex = (midpoint == (sortedAtoms.size()-1) ) ? midpoint : (midpoint+1);
    while( atomIndex >= 0 && cube[0][0] <= sortedAtoms[atomIndex].coordinates[0] ){
        if( cube[0][1] <= sortedAtoms[atomIndex].coordinates[1] &&
            cube[1][1] >= sortedAtoms[atomIndex].coordinates[1] &&
            cube[0][2] <= sortedAtoms[atomIndex].coordinates[2] &&
            cube[1][2] >= sortedAtoms[atomIndex].coordinates[2] &&
            cube[0][0] <= sortedAtoms[atomIndex].coordinates[0] &&
            cube[1][0] >= sortedAtoms[atomIndex].coordinates[0] ){
            prospectiveAtoms.push_back( atomIndex );
        }
        atomIndex--;
    }

    atomIndex = midpoint > 0 ? midpoint - 1 : midpoint;
    while( atomIndex < sortedAtoms.size() && cube[1][0] >= sortedAtoms[atomIndex].coordinates[0] ){
        if( cube[0][1] <= sortedAtoms[atomIndex].coordinates[1] &&
            cube[1][1] >= sortedAtoms[atomIndex].coordinates[1] &&
            cube[0][2] <= sortedAtoms[atomIndex].coordinates[2] &&
            cube[1][2] >= sortedAtoms[atomIndex].coordinates[2] &&
            cube[0][0] <= sortedAtoms[atomIndex].coordinates[0] &&
            cube[1][0] >= sortedAtoms[atomIndex].coordinates[0] &&
            (prospectiveAtoms.size() == 0 || prospectiveAtoms[0] != atomIndex) ){
            prospectiveAtoms.push_back( atomIndex );
        }
        atomIndex++;
    }

    return;

}

// get Born sum

static void getBornSum( float atomCoordinates[3], struct IntegrationGrid& grid,
                        const int exponent, double& bornSum, double& minR, FILE* log ){

// ---------------------------------------------------------------------------------------

   static const std::string methodName      = "getBornSum";

// ---------------------------------------------------------------------------------------

    bornSum                 = 0.0;
    minR                    = 1.0e+30;
    int closeIndex          = 0;
    for( unsigned int ii = 0; ii < grid.totalNumberOfGridMaskPoints; ii++ ){

  
if( ii == -380316 ){
int gridIndex = 12170113;
checkGridSetting( atomCoordinates, grid, gridIndex, "BornSum", log );
}

        if( grid.mask[ii] ){
            for( unsigned int jj = 0; jj < 32; jj++ ){
                if( grid.mask[ii] & (1 << jj) ){
                    int gridIndices[3];
                    float gridCoordinates[3];
                    int gridIndex = ii*32 + jj;
                    getGridCoordinatesGivenIndex( gridIndex, grid, gridCoordinates, gridIndices, log );
                    double r2 = static_cast<double>(
                                (gridCoordinates[0] - atomCoordinates[0])*(gridCoordinates[0] - atomCoordinates[0]) +
                                (gridCoordinates[1] - atomCoordinates[1])*(gridCoordinates[1] - atomCoordinates[1]) +
                                (gridCoordinates[2] - atomCoordinates[2])*(gridCoordinates[2] - atomCoordinates[2]) );
                    if( exponent == 4 ){
                        bornSum += 1.0/(r2*r2);
                    } else if( exponent == 6 ){
                        bornSum += 1.0/(r2*r2*r2);
                    }
                    if( r2 < minR ){
                        minR       = r2;
                        closeIndex = gridIndex;
                    }
                }
            }
        }
    }

    minR = sqrt( minR );
fprintf( log, "BornSum closest index=%8d %14.6e\n", closeIndex, minR );

    return;
}

static void readBornRadiiFile( const std::string& fileName, std::vector<sortedAtom>& bornRadii, FILE* log ){ 

// ---------------------------------------------------------------------------------------

   static const std::string methodName      = "readBornRadiiFile";
   char buffer[1024];

// ---------------------------------------------------------------------------------------

    FILE* brFile          = NULL;
#ifdef WIN32
    fopen_s( &brFile, fileName.c_str(), "r" );
#else
    brFile = fopen( fileName.c_str(), "r" );
#endif
    if( brFile == NULL ){
        (void) fprintf( log, "Born radii file=%s not opened.\n", fileName.c_str() );
        return;
    } else {
        (void) fprintf( log, "Born radii file=%s opened.\n", fileName.c_str() );
    }

    fgets( buffer, 1024, brFile );
    int numberOfAtoms = atoi( buffer );
    int lineCount     = 0;
    if( numberOfAtoms > 0 ){
        bornRadii.resize( numberOfAtoms );
        for( int ii = 0; ii < numberOfAtoms; ii++ ){
            StringVector lineTokens;
            char* isNotEof = readLine( brFile, lineTokens, &lineCount, log );
            int tokenIndex = 0;
            if( lineTokens.size() >= 4 ){
                int index                      = atoi( lineTokens[tokenIndex++].c_str() );
                double bornRadius              = atof( lineTokens[tokenIndex++].c_str() );
                double radius                  = atof( lineTokens[tokenIndex++].c_str() );
                double scaledRadius            = atof( lineTokens[tokenIndex++].c_str() );
                int countCase1                 = atoi( lineTokens[tokenIndex++].c_str() );
                int countCase2                 = atoi( lineTokens[tokenIndex++].c_str() );
                int countCase3                 = atoi( lineTokens[tokenIndex++].c_str() );
                int countCase4                 = atoi( lineTokens[tokenIndex++].c_str() );
                double coord0                  = atof( lineTokens[tokenIndex++].c_str() );
                double coord1                  = atof( lineTokens[tokenIndex++].c_str() );
                double coord2                  = atof( lineTokens[tokenIndex++].c_str() );

                bornRadii[ii].listIndex        = ii; 
                bornRadii[ii].radius           = radius; 
                bornRadii[ii].bornRadius       = bornRadius; 
                bornRadii[ii].scaledRadius     = scaledRadius; 
                bornRadii[ii].countCases[0]    = countCase1; 
                bornRadii[ii].countCases[1]    = countCase2; 
                bornRadii[ii].countCases[2]    = countCase3; 
                bornRadii[ii].countCases[3]    = countCase4; 
                bornRadii[ii].coordinates[0]   = coord0; 
                bornRadii[ii].coordinates[1]   = coord1; 
                bornRadii[ii].coordinates[2]   = coord2; 
            }
        }
    }

    std::sort( bornRadii.begin(), bornRadii.end(), compareAtomCoordinate );

    (void) fclose( brFile );
    return;
}

static void calculateBornRadiiByDirectIntegration( std::string parameterFileName, MapStringInt forceMap,
                                                   MapStringString& inputArgumentMap, FILE* inputLog,
                                                   FILE* filePtr, int* lineCount, FILE* resultsFile ){

// ---------------------------------------------------------------------------------------

   static const std::string methodName      = "calculateBornRadiiByDirectIntegration";
   int PrintOn                              = 1; 

// ---------------------------------------------------------------------------------------

    FILE* log;
    if( PrintOn == 0 && inputLog ){
       log = stderr;
    } else {
       log = inputLog;
    } 
 
    if( log ){
       (void) fprintf( log, "%s\n", methodName.c_str() );
       (void) fflush( log );
    }   
 
    //Context* context = setupEthane( log );
    Context* context = setupTwoParticle( log );
    //Context* context = setupSmallMolecule( filePtr, lineCount, log );

#if 0
    ReferencePlatform referencePlatform;
    registerFreeEnergyMethodsReferencePlatform( referencePlatform );

    if( log ){
        (void) fprintf( log, "%s Testing reference platform\n", methodName.c_str() );
        (void) fflush( log );
    }   

   double parameterKineticEnergy, parameterPotentialEnergy;

    std::vector<Vec3> parameterForces;
    std::vector<Vec3> parameterForces2;
    MapStringVectorOfVectors supplementary;

    Context* context       = testSetup( parameterFileName, forceMap, referencePlatform, 
                                        parameterForces, &parameterKineticEnergy, &parameterPotentialEnergy,
                                        supplementary, inputArgumentMap, log );
 
#endif
 	 State state            = context->getState(State::Positions | State::Energy);
 	 const std::vector<Vec3>& positions = state.getPositions();
 
    // get radii and find maximum radius

    // std::string GBVI_FORCE                      = "GBVI";
    // std::string GBVI_SOFTCORE_FORCE             = "GBVISoftcore";
    std::vector<float> radii;
    System& system = context->getSystem();
    getGBVIRadii( system, GBVI_FORCE, radii, log );
    if( radii.size() < 1 ){
       (void) fprintf( log, "%s GB/VI force not found.\n", methodName.c_str() );
       (void) fflush( log );
       exit(-1);
    }   
 
    float maxAtomicRadius = getMaxRadius( radii, log );
    
    (void) fprintf( log, "maxAtomicRadius=%14.6e\n", maxAtomicRadius ); fflush( log );

    // get bounding dimensions of system 
 
    struct BoundingBox boundingBox;
    getBoundingBox( positions, boundingBox, maxAtomicRadius, log );
 
    // allocate grid and set all entries to zero
 
    struct IntegrationGrid grid;
    memset( &grid, 0, sizeof( IntegrationGrid ) );
    int bornExponent;
    setFloatFromMap( inputArgumentMap, "gridResolution",  grid.resolution );
    setIntFromMap( inputArgumentMap, "bornExponent",  bornExponent );

    getGrid( boundingBox, grid, log );
    printGrid( grid, log );
 
    // sort atoms by x-coordinate and y-coordinate

    std::vector<sortedAtom> sortedAtoms;
    sortAtomsByCoordinate( positions, radii, sortedAtoms, log );
#if 1
    (void) fprintf( log, "\nSorted atoms:\n" );
    for( unsigned int ii = 0; ii < sortedAtoms.size(); ii++ ){
        (void) fprintf( log, "%6d [%14.6e %14.6e %14.6e] r=%8.3f\n", ii,
                        sortedAtoms[ii].coordinates[0],
                        sortedAtoms[ii].coordinates[1],
                        sortedAtoms[ii].coordinates[2],
                        sortedAtoms[ii].radius );
        for( unsigned int jj = ii+1; jj < sortedAtoms.size(); jj++ ){
            float separation = getDistance( sortedAtoms[ii].coordinates, sortedAtoms[jj].coordinates );
            if( separation <  (sortedAtoms[ii].radius +  sortedAtoms[ii].radius) ){
                (void) fprintf( log, "    %6d sep=%18.3f r+r=%8.3f r_j=%8.3f [%14.6e %14.6e %14.6e]\n", jj,
                                separation, (sortedAtoms[ii].radius + sortedAtoms[jj].radius), sortedAtoms[jj].radius,
                                sortedAtoms[jj].coordinates[0], sortedAtoms[jj].coordinates[1], sortedAtoms[jj].coordinates[2] );
            }
        }
    }
#endif

    int debug = 0;
    if( debug ){
        unsigned int misses = 0;
        for( unsigned int ii = 0; ii < grid.totalNumberOfGridPoints; ii++ ){
            float gridCoordinates[3];
            int gridIndices1[3];
            getGridCoordinatesGivenIndex( ii, grid, gridCoordinates, gridIndices1, log );
            int gridIndex;
            int gridIndices2[3];
            getGridIndexGivenCoordinates( grid, gridCoordinates, gridIndex, gridIndices2, log );
            if( ii != gridIndex && misses++ < 100 ){
               (void) fprintf( log, "Miss grid check: %6u %6d [%14.6e %14.6e %14.6e] [%4d %4d %4d] %4d %4d %4d]\n",
                               ii, gridIndex,
                               gridCoordinates[0], gridCoordinates[1], gridCoordinates[2],
                               gridIndices1[0], gridIndices1[1], gridIndices1[2],
                               gridIndices2[0], gridIndices2[1], gridIndices2[2] );
            }
            if( misses == 100 )exit(0);
        }
        (void) fprintf( log, "Total misses=%u getGridCoordinatesGivenIndex/getGridCoordinatesGivenIndex \n", misses );
        exit(0);
    }

    // set mask for grid points that are internal to the molecule to 1

    const unsigned int inside  = 1;
    const unsigned int outside = 0;
    for( unsigned int ii = 0; ii < grid.totalNumberOfGridPoints; ii++ ){

        float gridCoordinates[3];
        int gridIndices[3];
        getGridCoordinatesGivenIndex( ii, grid, gridCoordinates, gridIndices, log );

        // check if any atoms close by

        float cube[2][3];
        for( unsigned int jj = 0; jj < 3; jj++ ){
            cube[0][jj] = gridCoordinates[jj] - maxAtomicRadius;
            cube[1][jj] = gridCoordinates[jj] + maxAtomicRadius;
        }
        std::vector<unsigned int> prospectiveAtoms;
        getAtomsInCube( cube, sortedAtoms, prospectiveAtoms, log );
 
        unsigned int gridPointEnclosed = 0;
        std::vector<unsigned int> enclosingAtoms;
        std::vector<float> distanceCubeAtom;
        for( unsigned int jj = 0; jj < prospectiveAtoms.size() && gridPointEnclosed == 0; jj++ ){
            unsigned int atomIndex = prospectiveAtoms[jj];
            float distance         = sqrtf(
                                         (sortedAtoms[atomIndex].coordinates[0] - gridCoordinates[0] )*
                                         (sortedAtoms[atomIndex].coordinates[0] - gridCoordinates[0] ) + 
                                         (sortedAtoms[atomIndex].coordinates[1] - gridCoordinates[1] )*
                                         (sortedAtoms[atomIndex].coordinates[1] - gridCoordinates[1] ) +
                                         (sortedAtoms[atomIndex].coordinates[2] - gridCoordinates[2] )*
                                         (sortedAtoms[atomIndex].coordinates[2] - gridCoordinates[2] ) );
            if( distance < sortedAtoms[atomIndex].radius ){
                //gridPointEnclosed++;
                sortedAtoms[atomIndex].gridCount++;
                enclosingAtoms.push_back( atomIndex );
                distanceCubeAtom.push_back( distance );
                maskGridPoint( ii, grid, inside, log );
            }
        }

        //if( log && gridPointEnclosed )
        if( 0 && log && enclosingAtoms.size() ){
            (void) fprintf( log, "%3u atoms out of %3u prospective atoms in cube [%14.6e %14.6e%14.6e] [%14.6e %14.6e%14.6e]\n",
                            gridPointEnclosed, prospectiveAtoms.size(),
                            cube[0][0], cube[0][1], cube[0][2],
                            cube[1][0], cube[1][1], cube[1][2] );
            for( unsigned int jj = 0; jj < enclosingAtoms.size(); jj++ ){
                unsigned int atomIndex = enclosingAtoms[jj];
                (void) fprintf( log, "  %3u %3u %14.4e [%14.6e %14.6e%14.6e]\n", jj, atomIndex, distanceCubeAtom[jj],
                                sortedAtoms[atomIndex].coordinates[0],
                                sortedAtoms[atomIndex].coordinates[1], 
                                sortedAtoms[atomIndex].coordinates[2] );
            }
        }

    }

#if 0
    (void) fprintf( log, "\nSorted atoms:\n" );
    for( unsigned int ii = 0; ii < sortedAtoms.size(); ii++ ){
        float factor          = sortedAtoms[ii].radius/grid.resolution;
        float volume          = 4.1887902f*factor*factor*factor;
        int expectedGridCount = static_cast<int>(volume);
        (void) fprintf( log, "%6d [%14.6e %14.6e %14.6e] r=%14.6e %8u %8d\n", ii,
                        sortedAtoms[ii].coordinates[0],
                        sortedAtoms[ii].coordinates[1],
                        sortedAtoms[ii].coordinates[2],
                        sortedAtoms[ii].radius, sortedAtoms[ii].gridCount, expectedGridCount );
    }
    (void) fflush( log );
#endif

    // save grid result
 
    memcpy( grid.savedMask, grid.mask, sizeof(int)*grid.totalNumberOfGridMaskPoints );
//checkGridSetting( sortedAtoms[0].coordinates, grid, 12170113, "PostSave", log );

    unsigned int maskedTotal;
    getTotalMasked( grid, maskedTotal );
    double delta = static_cast<double>(bornExponent - 3 );

    // find Born radius for each atom

    for( unsigned int ii = 0; ii < sortedAtoms.size(); ii++ ){
    //for( unsigned int ii = 0; ii < 1; ii++ ){

        memcpy( grid.mask, grid.savedMask, sizeof(int)*grid.totalNumberOfGridMaskPoints );
//checkGridSetting( sortedAtoms[0].coordinates, grid, 12170113, "PostCopy", log );

        float radiusSquared = sortedAtoms[ii].radius*sortedAtoms[ii].radius;

float closest        = 1.0e+30f;
int   closestIndex   = -1;

        // set grid points enclosed by atom to 0 in temp array

        float gridCoordinates[3];
        for( int jj = 0; jj < grid.numberOfGridPoints[0]; jj++ ){
            gridCoordinates[0] = grid.origin[0] + static_cast<float>(jj)*grid.resolution;    
            if( fabsf( gridCoordinates[0] - sortedAtoms[ii].coordinates[0] ) < sortedAtoms[ii].radius ){ 
                for( int kk = 0; kk < grid.numberOfGridPoints[1]; kk++ ){

                    gridCoordinates[1] = grid.origin[1] + static_cast<float>(kk)*grid.resolution;    

                    float distance01   = (gridCoordinates[0] - sortedAtoms[ii].coordinates[0] )*(gridCoordinates[0] - sortedAtoms[ii].coordinates[0] ) +
                                         (gridCoordinates[1] - sortedAtoms[ii].coordinates[1] )*(gridCoordinates[1] - sortedAtoms[ii].coordinates[1] );

                    if( distance01 < radiusSquared ){
                        for( int mm = 0; mm < grid.numberOfGridPoints[2]; mm++ ){

                            gridCoordinates[2] = grid.origin[2] + static_cast<float>(mm)*grid.resolution;    

                            float distance     = distance01 +
                                                 (sortedAtoms[ii].coordinates[2] - gridCoordinates[2] )*(sortedAtoms[ii].coordinates[2] - gridCoordinates[2] );

                            if( distance < radiusSquared ){
                                 int gridIndices2[3];
                                 int gridIndex2;
                                 getGridIndexGivenCoordinates( grid, gridCoordinates, gridIndex2, gridIndices2, log );
                                int gridIndex = jj + kk*grid.offset[1] + mm*grid.offset[2];
                                maskGridPoint( gridIndex, grid, outside, log );
                                sortedAtoms[ii].maskedCount++;
                                if( gridIndex2 != gridIndex ){
                                    (void) fprintf( log, "Index issue: Atom=%6d [%7d %7d] [%5d %5d %5d ] [%5d %5d %5d] [%14.6e %14.6e %14.6e]\n", ii,
                                                    gridIndex, gridIndex2, jj, kk, mm, gridIndices2[0], gridIndices2[1], gridIndices2[2],
                                                    gridCoordinates[0], gridCoordinates[1], gridCoordinates[2] );
                                }
                            } else if( distance < closest ){
                                int gridIndex = jj + kk*grid.offset[1] + mm*grid.offset[2];
                                int isSet     = isGridIndexSet( grid, gridIndex, log );
                                if( isSet ){
                                    (void) fprintf( log, "Close: Atom=%6d %14.6e %7d [%5d %5d %5d ] Set=%d\n", ii,
                                                    sqrtf( distance ), gridIndex, jj, kk, mm, isSet );
                                    closest      = distance;
                                    closestIndex = gridIndex;
                                }
                            }
                        }
                    }
                }
            }
        }

        getTotalMasked( grid, sortedAtoms[ii].bornMaskedCount );

        // sum 1/r**n for points nonzero

        getBornSum( sortedAtoms[ii].coordinates, grid, bornExponent, sortedAtoms[ii].bornSum, sortedAtoms[ii].minR, log );
        sortedAtoms[ii].bornSum      *= (grid.resolution*grid.resolution*grid.resolution);
        sortedAtoms[ii].bornRadius    = std::pow( static_cast<double>(sortedAtoms[ii].radius), -delta ) -
                                        delta*0.0795774715*sortedAtoms[ii].bornSum;
        sortedAtoms[ii].bornRadius    = 1.0/std::pow( sortedAtoms[ii].bornRadius, 1.0/delta );

(void) fprintf( log, "Closest: Atom=%6d %7d %14.6e\n", ii, closestIndex, sqrtf( closest ) );
        // compute Born radius given volume

    }

    (void) fprintf( log, "\nSorted atoms: total masked=%u\n", maskedTotal );
    std::string bornRadiiFile = "br.txt"; 
    std::vector<sortedAtom> gbviRadii;
    readBornRadiiFile( bornRadiiFile, gbviRadii, log );
    for( unsigned int jj = 0; jj < sortedAtoms.size(); jj++ ){
        float factor          = sortedAtoms[jj].radius/grid.resolution;
        float volume          = 4.1887902f*factor*factor*factor;
        int expectedGridCount = static_cast<int>(volume);
        (void) fprintf( log, "%6d r=%8.3f %8u %8d masked=%8d bMask=%8u Bsum=%14.6e bR=%14.6e %14.6e minR=%14.6e",
                        sortedAtoms[jj].listIndex,
                        sortedAtoms[jj].radius, sortedAtoms[jj].gridCount, expectedGridCount, 
                        sortedAtoms[jj].maskedCount,
                        sortedAtoms[jj].bornMaskedCount,
                        sortedAtoms[jj].bornSum, 
                        sortedAtoms[jj].bornRadius, gbviRadii[jj].bornRadius,
                        sortedAtoms[jj].minR );
        if( sortedAtoms[jj].listIndex != gbviRadii[jj].listIndex ){
            (void) fprintf( log, " Indices: %5d %5d XXXX", sortedAtoms[jj].listIndex != gbviRadii[jj].listIndex );
        }
        if( 0 ){
            (void) fprintf( log, " [%14.6e %14.6e %14.6e]",
                        sortedAtoms[jj].coordinates[0],
                        sortedAtoms[jj].coordinates[1],
                        sortedAtoms[jj].coordinates[2] );
        }
        (void) fprintf( log, "\n" );
    }
    (void) fflush( log );

    FILE* brFile          = resultsFile;
    int closeBrFile       = 0;
    if( brFile == NULL ){
        closeBrFile       = 1;
        std::string fileName  = "BornResults.txt";
#ifdef WIN32
        fopen_s( &brFile, fileName.c_str(), "w" );
#else
        brFile = fopen( fileName.c_str(), "w" );
#endif
        if( brFile == NULL ){
            (void) fprintf( log, "Born radii file=%s not opened.\n", fileName.c_str() );
            return;
        } else {
            (void) fprintf( log, "Born radii file=%s opened.\n", fileName.c_str() );
        }
    }
    
    for( unsigned int jj = 0; jj < sortedAtoms.size(); jj++ ){
        (void) fprintf( brFile, "%6d %14.6e %14.6e %14.6e %14.6e %8.3f  %8.3f %5d %5d %5d %5d\n", jj,
                        sortedAtoms[jj].bornRadius, gbviRadii[jj].bornRadius,
                        1.0/sortedAtoms[jj].bornRadius, 1.0/gbviRadii[jj].bornRadius,
                        sortedAtoms[jj].radius, gbviRadii[jj].scaledRadius, 
                        gbviRadii[jj].countCases[0], gbviRadii[jj].countCases[1], gbviRadii[jj].countCases[2], gbviRadii[jj].countCases[3] );
    }
    (void) fflush( brFile );
    if( closeBrFile ){
        (void) fclose( brFile );
    }

    // compute Born radius given volume
  
    if( grid.mask ){
        free( grid.mask );
        free( grid.savedMask );
    }
    if( context ){
        delete context;
    }
}

// ---------------------------------------------------------------------------------------
// GB/VI test end
// ---------------------------------------------------------------------------------------

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
   int checkPlatformForces                           = 0;
   int checkEnergyForceConsistent                    = 0;
   int checkEnergyConservation                       = 0;
   int checkInputForces                              = 0;
   int checkBornRadii                                = 0;
   int checkSmallMoleculeBornRadii                   = 0;
   MapStringString inputArgumentMap;

   FILE* log                                         = NULL;
   FILE* summaryFile                                 = NULL;

// ---------------------------------------------------------------------------------------

   std::string defaultParameterFileName     = "OpenParameters.txt";

   if( numberOfArguments < 2 ){
      printUsage( defaultParameterFileName );
   }

   std::string parameterFileName            = defaultParameterFileName;
   std::string resultsFileName              = defaultParameterFileName;
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
      } else if( STRCASECMP( argv[ii], "-resultsFileName" ) == 0 ){
         resultsFileName          = argv[ii+1];
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
      } else if( STRCASECMP( argv[ii], "-checkPlatformForces" ) == 0 ){
         checkPlatformForces        = atoi( argv[ii+1] );
         ii++;
      } else if( STRCASECMP( argv[ii], "-checkInputForces" ) == 0 ){
         checkInputForces           = atoi( argv[ii+1] );
         ii++;
      } else if( STRCASECMP( argv[ii], "-checkEnergyForceConsistent" ) == 0 ){
         checkEnergyForceConsistent = atoi( argv[ii+1] );
         ii++;
      } else if( STRCASECMP( argv[ii], "-checkSmallMoleculeBornRadii" ) == 0 ){
         checkSmallMoleculeBornRadii = atoi( argv[ii+1] );
         ii++;
      } else if( STRCASECMP( argv[ii], "-checkBornRadii" ) == 0 ){
         checkBornRadii              = atoi( argv[ii+1] );
         ii++;
      } else if( STRCASECMP( argv[ii], "-checkEnergyConservation" ) == 0 ){
         checkEnergyConservation = atoi( argv[ii+1] );;
         ii++;
      } else if( STRCASECMP( argv[ii], "-energyForceDelta" )     == 0    ||
                 STRCASECMP( argv[ii], "-energyForceTolerance" ) == 0    ||
                 STRCASECMP( argv[ii], "-cudaDeviceId" )         == 0    ||
                 STRCASECMP( argv[ii], "-platform1" )            == 0    ||
                 STRCASECMP( argv[ii], "-platform2" )            == 0    ||
                 STRCASECMP( argv[ii], "-custom1" )              == 0    ||
                 STRCASECMP( argv[ii], "-custom2" )              == 0    ||
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
      } else if( STRNCASECMP( argv[ii], "-equilibration",           14  ) == 0 ||
                 STRNCASECMP( argv[ii], "-simulation",              11  ) == 0 ||
                 STRNCASECMP( argv[ii], "-comparisonPlatform",      19  ) == 0 ||
                 STRNCASECMP( argv[ii], "-runId",                    6  ) == 0 ||
                 STRNCASECMP( argv[ii], "-nonbonded",               10  ) == 0 ||
                 STRNCASECMP( argv[ii], "-readContext",             12  ) == 0 ||
                 STRNCASECMP( argv[ii], "-gridResolution",          15  ) == 0 ||
                 STRNCASECMP( argv[ii], "-bornExponent",            13  ) == 0 ){
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

      (void) fprintf( log, "checkEnergyForceConsistent  %d\n", checkEnergyForceConsistent );
      (void) fprintf( log, "checkEnergyConservation     %d\n", checkEnergyConservation );
      (void) fprintf( log, "checkForces                 %d\n", checkForces );
      (void) fprintf( log, "checkPlatformForces         %d\n", checkPlatformForces );
      (void) fprintf( log, "checkInputForces            %d\n", checkInputForces );
      (void) fprintf( log, "checkBornRadii              %d\n", checkBornRadii );
      (void) fprintf( log, "checkSmallMoleculeBornRadii %d\n", checkSmallMoleculeBornRadii );

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

   // check forces
   // deprecated -- replace w/ checkPlatformForces

   if( checkForces ){
      try {
         testReferenceCudaForces( parameterFileName, forceMap, inputArgumentMap, log, summaryFile );
       } catch( const exception& e ){
         (void) fprintf( stderr, "Exception checkForces %s %s\n", methodName.c_str(),  e.what() ); (void) fflush( stderr );
         return 1;
      }   
   }

   if( checkPlatformForces ){
      try {
         testForces( parameterFileName, forceMap, inputArgumentMap, log, summaryFile );
       } catch( const exception& e ){
         (void) fprintf( stderr, "Exception checkPlatformForces %s %s\n", methodName.c_str(),  e.what() ); (void) fflush( stderr );
         return 1;
      }   
   }

   // check Born radii

   if( checkBornRadii ){
      FILE* brFile          = NULL;
      int lineCount         = 0;
      try {
         calculateBornRadiiByDirectIntegration( parameterFileName, forceMap, inputArgumentMap, log, brFile, &lineCount, NULL );
       } catch( const exception& e ){
         (void) fprintf( stderr, "Exception checkBornRadii %s %s\n", methodName.c_str(),  e.what() ); (void) fflush( stderr );
         return 1;
      }   
   }

   // check Born radii

   if( checkSmallMoleculeBornRadii ){
      FILE* brFile          = NULL;
      FILE* resultsFile     = NULL;
      int lineCount         = 0;
      try {
#ifdef WIN32
         fopen_s( &brFile, parameterFileName.c_str(), "r" );
#else
         brFile = fopen( parameterFileName.c_str(), "r" );
#endif
#ifdef WIN32
         fopen_s( &resultsFile, resultsFileName.c_str(), "w" );
#else
         resultsFile = fopen( resultsFileName.c_str(), "w" );
#endif
         if( brFile == NULL ){
            (void) fprintf( log, "parameterFile=%s not opened.\n", parameterFileName.c_str() );
         } else {
           (void) fprintf( log, "parameterFile=%s opened.\n", parameterFileName.c_str() );
         }

         if( resultsFile == NULL ){
            (void) fprintf( log, "resultsFile=%s not opened.\n", resultsFileName.c_str() );
         } else {
           (void) fprintf( log, "resultsFile=%s opened.\n", resultsFileName.c_str() );
         }

         while( lineCount < 18596 ){
             calculateBornRadiiByDirectIntegration( parameterFileName, forceMap, inputArgumentMap, log, brFile, &lineCount, resultsFile );
         }
         (void) fclose( brFile );
         (void) fclose( resultsFile );
       } catch( const exception& e ){
         (void) fprintf( stderr, "Exception checkBornRadii %s %s\n", methodName.c_str(),  e.what() ); (void) fflush( stderr );
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
