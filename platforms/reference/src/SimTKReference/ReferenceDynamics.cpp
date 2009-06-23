
/* Portions copyright (c) 2006 Stanford University and Simbios.
 * Contributors: Pande Group
 *
 * Permission is hereby granted, free of charge, to any person obtaining
 * a copy of this software and associated documentation files (the
 * "Software"), to deal in the Software without restriction, including
 * without limitation the rights to use, copy, modify, merge, publish,
 * distribute, sublicense, and/or sell copies of the Software, and to
 * permit persons to whom the Software is furnished to do so, subject
 * to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included
 * in all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
 * OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
 * MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
 * IN NO EVENT SHALL THE AUTHORS, CONTRIBUTORS OR COPYRIGHT HOLDERS BE
 * LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
 * OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION
 * WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 */

#include <string.h>
#include <sstream>

#include "../SimTKUtilities/SimTKOpenMMCommon.h"
#include "../SimTKUtilities/SimTKOpenMMLog.h"
#include "../SimTKUtilities/SimTKOpenMMUtilities.h"
#include "ReferenceDynamics.h"

const int ReferenceDynamics::DefaultReturn      = 0;
const int ReferenceDynamics::ErrorReturn        = -1;


/**---------------------------------------------------------------------------------------

   ReferenceDynamics constructor

   @param numberOfAtoms  number of atoms
   @param deltaT         delta t for dynamics
   @param temperature    temperature

   --------------------------------------------------------------------------------------- */

ReferenceDynamics::ReferenceDynamics( int numberOfAtoms,  RealOpenMM deltaT, RealOpenMM temperature ) : 
                  _numberOfAtoms(numberOfAtoms), _deltaT(deltaT), _temperature(temperature) {

   // ---------------------------------------------------------------------------------------

   // static const char* methodName = "\nReferenceDynamics::ReferenceDynamics";

   static const RealOpenMM one        =  1.0;

   // ---------------------------------------------------------------------------------------

   _timeStep             = 0;

   _twoDTempArrays       = 0;
   _twoDTempArrays       = NULL;

   _oneDTempArrays       = 0;
   _oneDTempArrays       = NULL;

   _ownReferenceConstraint = false;
   _referenceConstraint    = NULL;
}

/**---------------------------------------------------------------------------------------

   ReferenceDynamics destructor

   --------------------------------------------------------------------------------------- */

ReferenceDynamics::~ReferenceDynamics( ){

   // ---------------------------------------------------------------------------------------

   // static const char* methodName = "\nReferenceDynamics::~ReferenceDynamics";

   // ---------------------------------------------------------------------------------------

   _freeTwoDArrays();
   _freeOneDArrays();

   if( _ownReferenceConstraint ){
      delete _referenceConstraint;
   }
}

/**---------------------------------------------------------------------------------------

   Free memory associated w/ 2D arrays

   @return ReferenceDynamics::DefaultReturn

   --------------------------------------------------------------------------------------- */

int ReferenceDynamics::_freeTwoDArrays( void ){

   // ---------------------------------------------------------------------------------------

   // static const char* methodName = "\nReferenceDynamics::_freeTwoDArrays";

   // ---------------------------------------------------------------------------------------

   if( _twoDTempArrays ){
      delete[] _twoDTempArrays[0][0];
      delete[] _twoDTempArrays[0];
      delete[] _twoDTempArrays;
   }

   _twoDTempArrays        = NULL;
   _numberOf2DTempArrays  = 0;

   return ReferenceDynamics::DefaultReturn;

}

/**---------------------------------------------------------------------------------------

   Free memory associated w/ 1D arrays

   @return ReferenceDynamics::DefaultReturn

   --------------------------------------------------------------------------------------- */

int ReferenceDynamics::_freeOneDArrays( void ){

   // ---------------------------------------------------------------------------------------

   // static const char* methodName = "\nReferenceDynamics::_freeOneDArrays";

   // ---------------------------------------------------------------------------------------

   if( _oneDTempArrays ){
      delete[] _oneDTempArrays[0];
      delete[] _oneDTempArrays;
   }

   _oneDTempArrays        = NULL;
   _numberOf1DTempArrays  = 0;

   return ReferenceDynamics::DefaultReturn;

}

/**---------------------------------------------------------------------------------------

   Allocate memory associated w/ 2D arrays

   @param dimension1        first dimension
   @param dimension2        second dimension
   @param numberOfArrays    number of arrays to allocate

   @return ReferenceDynamics::DefaultReturn

   --------------------------------------------------------------------------------------- */

int ReferenceDynamics::allocate2DArrays( int dimension1, int dimension2, int numberOfArrays ){

   // ---------------------------------------------------------------------------------------

   // static const char* methodName = "\nReferenceDynamics::allocate2DArrays";

   // ---------------------------------------------------------------------------------------

   _freeTwoDArrays();

   _numberOf2DTempArrays   = numberOfArrays;
   _twoDTempArrays         = new RealOpenMM**[_numberOf2DTempArrays];

   RealOpenMM** totalArray = new RealOpenMM*[dimension1*numberOfArrays];

   RealOpenMM*  totalBlock = new RealOpenMM[dimension1*dimension2*numberOfArrays];
   memset( totalBlock, 0, sizeof( RealOpenMM )*dimension1*dimension2*numberOfArrays );

   for( int ii = 0; ii < _numberOf2DTempArrays; ii++ ){
      _twoDTempArrays[ii]    = totalArray;
      totalArray            += dimension1; 

      for( int jj = 0; jj < dimension1; jj++ ){
         _twoDTempArrays[ii][jj]   = totalBlock;
         totalBlock               += dimension2;
      }
   }

   return ReferenceDynamics::DefaultReturn;

}

/**---------------------------------------------------------------------------------------

   Get array at specified index

   @param index             array index

   @return array or NULL if index invalid or arrays not allocated

   --------------------------------------------------------------------------------------- */

RealOpenMM** ReferenceDynamics::get2DArrayAtIndex( int index ) const {

   // ---------------------------------------------------------------------------------------

   static const char* methodName = "\nReferenceDynamics::get2DArrayAtIndex";

   // ---------------------------------------------------------------------------------------

   if( index < 0 || index >= _numberOf2DTempArrays ){
      std::stringstream message;
      message << methodName;
      message << " requested 2d array at index=" << index << " is unavailable.";
      SimTKOpenMMLog::printError( message );
      return NULL;
   }

   return _twoDTempArrays[index];

}

/**---------------------------------------------------------------------------------------

   Allocate memory associated w/ 1D arrays

   @param dimension1        dimension
   @param numberOfArrays    number of arrays to allocate

   @return ReferenceDynamics::DefaultReturn

   --------------------------------------------------------------------------------------- */

int ReferenceDynamics::allocate1DArrays( int dimension, int numberOfArrays ){

   // ---------------------------------------------------------------------------------------

   // static const char* methodName = "\nReferenceDynamics::allocate1DArrays";

   // ---------------------------------------------------------------------------------------

   _freeOneDArrays();

   _numberOf1DTempArrays   = numberOfArrays;
   _oneDTempArrays         = new RealOpenMM*[_numberOf1DTempArrays];

   RealOpenMM* totalArray  = new RealOpenMM[dimension*numberOfArrays];
   memset( totalArray, 0, sizeof( RealOpenMM )*dimension*numberOfArrays );

   for( int ii = 0; ii < _numberOf1DTempArrays; ii++ ){
      _oneDTempArrays[ii]    = totalArray;
      totalArray            += dimension; 
   }

   return ReferenceDynamics::DefaultReturn;

}

/**---------------------------------------------------------------------------------------

   Get array at specified index

   @param index             array index

   @return array or NULL if index invalid or arrays not allocated

   --------------------------------------------------------------------------------------- */

RealOpenMM* ReferenceDynamics::get1DArrayAtIndex( int index ) const {

   // ---------------------------------------------------------------------------------------

   static const char* methodName = "\nReferenceDynamics::get1DArrayAtIndex";

   // ---------------------------------------------------------------------------------------

   if( index < 0 || index >= _numberOf1DTempArrays ){
      std::stringstream message;
      message << methodName;
      message << " requested 1d array at index=" << index << " is unavailable.";
      SimTKOpenMMLog::printError( message );
      return NULL;
   }

   return _oneDTempArrays[index];

}

/**---------------------------------------------------------------------------------------

   Get number of atoms

   @return number of atoms

   --------------------------------------------------------------------------------------- */

int ReferenceDynamics::getNumberOfAtoms( void ) const {

   // ---------------------------------------------------------------------------------------

   // static const char* methodName  = "\nReferenceDynamics::getNumberOfAtoms";

   // ---------------------------------------------------------------------------------------

   return _numberOfAtoms;
}

/**---------------------------------------------------------------------------------------

   Get time step

   @return time step

   --------------------------------------------------------------------------------------- */

int ReferenceDynamics::getTimeStep( void ) const {

   // ---------------------------------------------------------------------------------------

   // static const char* methodName  = "\nReferenceDynamics::getTimeStep";

   // ---------------------------------------------------------------------------------------

   return _timeStep;
}

/**---------------------------------------------------------------------------------------

   Increment time step

   @return incremented time step

   --------------------------------------------------------------------------------------- */

int ReferenceDynamics::incrementTimeStep( void ){

   // ---------------------------------------------------------------------------------------

   // static const char* methodName  = "\nReferenceDynamics::getTimeStep";

   // ---------------------------------------------------------------------------------------

   return (++_timeStep);
}

/**---------------------------------------------------------------------------------------

   Get delta t

   @return deltaT

   --------------------------------------------------------------------------------------- */

RealOpenMM ReferenceDynamics::getDeltaT( void ) const {

   // ---------------------------------------------------------------------------------------

   // static const char* methodName  = "\nReferenceDynamics::getDeltaT";

   // ---------------------------------------------------------------------------------------

   return _deltaT;
}

/**---------------------------------------------------------------------------------------

   Set delta t

   --------------------------------------------------------------------------------------- */

void ReferenceDynamics::setDeltaT( RealOpenMM deltaT ) {

   // ---------------------------------------------------------------------------------------

   // static const char* methodName  = "\nReferenceDynamics::setDeltaT";

   // ---------------------------------------------------------------------------------------

   _deltaT = deltaT;
}

/**---------------------------------------------------------------------------------------

   Get temperature

   @return temperature

   --------------------------------------------------------------------------------------- */

RealOpenMM ReferenceDynamics::getTemperature( void ) const {

   // ---------------------------------------------------------------------------------------

   // static const char* methodName  = "\nReferenceDynamics::getTemperature";

   // ---------------------------------------------------------------------------------------

   return _temperature;
}

/**---------------------------------------------------------------------------------------

   Get ReferenceConstraint

   @return ReferenceConstraint  object

   --------------------------------------------------------------------------------------- */

ReferenceConstraintAlgorithm* ReferenceDynamics::getReferenceConstraintAlgorithm( void ) const {

   // ---------------------------------------------------------------------------------------

   // static const char* methodName  = "\nReferenceDynamics::getReferenceConstraint";

   // ---------------------------------------------------------------------------------------

   return _referenceConstraint;
}

/**---------------------------------------------------------------------------------------

   Set ReferenceConstraint

   @param referenceConstraint  ReferenceConstraint object

   @return ReferenceDynamics::DefaultReturn

   --------------------------------------------------------------------------------------- */

int ReferenceDynamics::setReferenceConstraintAlgorithm( ReferenceConstraintAlgorithm* referenceConstraint ){

   // ---------------------------------------------------------------------------------------

   // static const char* methodName  = "\nReferenceDynamics::setReferenceConstraint";

   // ---------------------------------------------------------------------------------------

   // delete if own

   if( _referenceConstraint && _ownReferenceConstraint ){
      delete _referenceConstraint;
   }

   _referenceConstraint = referenceConstraint;
   _ownReferenceConstraint = 0;

   return ReferenceDynamics::DefaultReturn;
}

/**---------------------------------------------------------------------------------------

   Update -- driver routine for performing stochastic dynamics update of coordinates
   and velocities

   @param numberOfAtoms       number of atoms
   @param atomCoordinates     atom coordinates
   @param velocities          velocities
   @param forces              forces
   @param masses              atom masses

   @return ReferenceDynamics::DefaultReturn

   --------------------------------------------------------------------------------------- */

int ReferenceDynamics::update( int numberOfAtoms, RealOpenMM** atomCoordinates,
                               RealOpenMM** velocities, RealOpenMM** forces, RealOpenMM* masses ){

   // ---------------------------------------------------------------------------------------

   static const char* methodName      = "\nReferenceDynamics::update";

   static const RealOpenMM zero       =  0.0;
   static const RealOpenMM one        =  1.0;

   static int debug                   =  0;  

   // ---------------------------------------------------------------------------------------

   return ReferenceDynamics::DefaultReturn;
}

/**---------------------------------------------------------------------------------------

   Remove total linear momentum

   @param numberOfAtoms      number of atoms
   @param masses             masses
   @param velocities         velocities

   @return ReferenceDynamics::DefaultReturn

   --------------------------------------------------------------------------------------- */

int ReferenceDynamics::removeTotalLinearMomentum( int numberOfAtoms, RealOpenMM* masses,
                                                  RealOpenMM** velocities ) const {

   // ---------------------------------------------------------------------------------------

   static const char* methodName  = "\nReferenceDynamics::removeTotalLinearMomentum";

   static const int debug              = 1;

   static const RealOpenMM zero        = 0.0;
   static const RealOpenMM one         = 1.0;

   // ---------------------------------------------------------------------------------------

   RealOpenMM totalMass          = zero;
   RealOpenMM linearMomentum[3]  = { zero, zero, zero };
   for( int ii = 0; ii < numberOfAtoms; ii++ ){
      totalMass          += masses[ii];
      linearMomentum[0]  += masses[ii]*velocities[ii][0];
      linearMomentum[1]  += masses[ii]*velocities[ii][1];
      linearMomentum[2]  += masses[ii]*velocities[ii][2];
/*
      if( debug ){
         std::stringstream message;
         message << ii << " TotalM=" << totalMass << " m=" << masses[ii] << " linearMomentum [";
         SimTKOpenMMUtilities::formatRealStringStream( message, linearMomentum, 3, one );
         message << "] [";
         SimTKOpenMMUtilities::formatRealStringStream( message, velocities[ii], 3, one );
         message << "]\n";
         SimTKOpenMMLog::printMessage( message );
      }
*/
   }

   if( totalMass > zero ){

      RealOpenMM totalMassI      = one/totalMass;

/*
      if( debug ){
         std::stringstream message;
         message << " Pre scale linearMomentum [";
         SimTKOpenMMUtilities::formatRealStringStream( message, linearMomentum, 3, one );
         message << "]\n";
         SimTKOpenMMLog::printMessage( message );
      }
*/

      linearMomentum[0]         *= totalMassI;
      linearMomentum[1]         *= totalMassI;
      linearMomentum[2]         *= totalMassI;

/*
      if( debug ){
         std::stringstream message;
         message << " Pre sub linearMomentum [";
         SimTKOpenMMUtilities::formatRealStringStream( message, linearMomentum, 3, one );
         message << "]\n";
         SimTKOpenMMLog::printMessage( message );
      }
*/
      for( int ii = 0; ii < numberOfAtoms; ii++ ){
         velocities[ii][0]     -= linearMomentum[0];
         velocities[ii][1]     -= linearMomentum[1];
         velocities[ii][2]     -= linearMomentum[2];
      }

      // debug

      if( debug ){
         RealOpenMM tempV[3];
         std::stringstream message;
         int maxPrint = 5;
         message << methodName;
         message << " TotalM=" << totalMass << " linearMomentum: [";
         SimTKOpenMMUtilities::formatRealStringStream( message, linearMomentum, 3, one );
         message << "]\nSample v:\n";
         for( int ii = 0; ii < maxPrint; ii++ ){
            tempV[0] = velocities[ii][0] + linearMomentum[0];
            tempV[1] = velocities[ii][1] + linearMomentum[1];
            tempV[2] = velocities[ii][2] + linearMomentum[2];
            message << "OldV[ ";
            SimTKOpenMMUtilities::formatRealStringStream( message, tempV, 3, one );
            message << "] NewV[ ";
            SimTKOpenMMUtilities::formatRealStringStream( message, velocities[ii], 3, one );
            message << "]" << std::endl;
         }
         SimTKOpenMMLog::printMessage( message );
      }
   }

   return ReferenceDynamics::DefaultReturn;
}

/**---------------------------------------------------------------------------------------

   Print parameters

   @param message             message

   @return ReferenceDynamics::DefaultReturn

   --------------------------------------------------------------------------------------- */

int ReferenceDynamics::printParameters( std::stringstream& message ) const {

   // ---------------------------------------------------------------------------------------

   static const char* methodName  = "\nReferenceDynamics::printParameters";

   // ---------------------------------------------------------------------------------------

   message << "atoms=" << getNumberOfAtoms()     << " ";
   message << "delta_t=" << getDeltaT()          << " ";
   message << "temperature=" << getTemperature() << " ";
   message << "step=" << getTimeStep()           << " ";
   message << "seed=" << SimTKOpenMMUtilities::getRandomNumberSeed()   << " ";
   message << std::endl;

   return ReferenceDynamics::DefaultReturn;
}

/**---------------------------------------------------------------------------------------

   Write state

   @param numberOfAtoms       number of atoms
   @param atomCoordinates     atom coordinates
   @param velocities          velocities
   @param forces              forces
   @param masses              atom masses
   @param state               0 if initial state; otherwise nonzero
   @param baseFileName        base file name

   @return ReferenceDynamics::DefaultReturn

   --------------------------------------------------------------------------------------- */

int ReferenceDynamics::writeState( int numberOfAtoms, RealOpenMM** atomCoordinates,
                                   RealOpenMM** velocities,
                                   RealOpenMM** forces, RealOpenMM* masses,
                                   int state, const std::string& baseFileName ) const {

   // ---------------------------------------------------------------------------------------

   static const char* methodName      = "\nReferenceStochasticDynamics::writeState";

   static const RealOpenMM zero       =  0.0;
   static const RealOpenMM one        =  1.0;
   static const int threeI            =  3;

   // ---------------------------------------------------------------------------------------

   std::stringstream stateFileName;

   stateFileName << baseFileName;
   stateFileName << "_Step" << getTimeStep();
   // stateFileName << "_State" << state;
   stateFileName << ".txt";

   // ---------------------------------------------------------------------------------------

   // open file -- return if unsuccessful

   FILE* stateFile = NULL;
#ifdef WIN32
   fopen_s( &stateFile, stateFileName.str().c_str(), "w" );
#else
   stateFile = fopen( stateFileName.str().c_str(), "w" );
#endif

   // ---------------------------------------------------------------------------------------

   // diagnostics

   if( stateFile != NULL ){
      std::stringstream message;
      message << methodName;
      message << " Opened file=<" << stateFileName.str() << ">.\n";
      SimTKOpenMMLog::printMessage( message );
   } else {
      std::stringstream message;
      message << methodName;
      message << " could not open file=<" << stateFileName.str() << "> -- abort output.\n";
      SimTKOpenMMLog::printMessage( message );
      return ReferenceDynamics::ErrorReturn;
   }   

   // ---------------------------------------------------------------------------------------

   StringVector scalarNameI;
   IntVector scalarI;

   StringVector scalarNameR;
   RealOpenMMVector scalarR;

   StringVector scalarNameR1;
   RealOpenMMPtrVector scalarR1;

   StringVector scalarNameR2;
   RealOpenMMPtrPtrVector scalarR2;

   scalarI.push_back( getNumberOfAtoms() );
   scalarNameI.push_back( "Atoms" );

   scalarI.push_back( getTimeStep() );
   scalarNameI.push_back( "Timestep" );

   scalarR.push_back( getDeltaT() );
   scalarNameR.push_back( "delta_t" );

   scalarR.push_back( static_cast<RealOpenMM>(SimTKOpenMMUtilities::getRandomNumberSeed()) );
   scalarNameR.push_back( "seed" );

   scalarR.push_back( getTemperature() );
   scalarNameR.push_back( "T" );

   if( masses ){
      scalarR1.push_back( masses );
      scalarNameR1.push_back( "mass" );
   }

   if( atomCoordinates ){
      scalarR2.push_back( atomCoordinates );
      scalarNameR2.push_back( "coord" );
   }

   if( velocities ){
      scalarR2.push_back( velocities );
      scalarNameR2.push_back( "velocities" );
   }

   if( forces ){
      scalarR2.push_back( forces );
      scalarNameR2.push_back( "forces" );
   }

   writeStateToFile( stateFile, scalarNameI, scalarI, scalarNameR, scalarR, getNumberOfAtoms(), scalarNameR1, scalarR1, threeI, scalarNameR2, scalarR2 ); 

   (void) fclose( stateFile );

   return ReferenceDynamics::DefaultReturn;

}

/**---------------------------------------------------------------------------------------

   Write state

   @param stateFile       file to write to
   @param scalarNameI     vector of scalar names for ints
   @param scalarI         vector of scalar ints
   @param scalarNameR     vector of scalar names for real
   @param scalarR         vector of scalar reals
   @param dimension1      size of first dimension for 1D & 2D real arrays
   @param scalarNameR1    vector of names for 1D real arrays
   @param scalarR1        vector of 1D real arrays
   @param dimension2      size of second dimension for 2D real arrays
   @param scalarNameR2    vector of names for 2D real arrays
   @param scalarR2        vector of 2D real arrays

   @return ReferenceDynamics::DefaultReturn

   --------------------------------------------------------------------------------------- */

int ReferenceDynamics::writeStateToFile( FILE* stateFile, 
                                         StringVector& scalarNameI, IntVector& scalarI,
                                         StringVector& scalarNameR, RealOpenMMVector& scalarR,
                                         int dimension1, StringVector& scalarNameR1,
                                         RealOpenMMPtrVector& scalarR1,
                                         int dimension2, StringVector& scalarNameR2,
                                         RealOpenMMPtrPtrVector& scalarR2 ) const {

   // ---------------------------------------------------------------------------------------

   static const char* methodName  = "\nReferenceDynamics::writeState";

   // ---------------------------------------------------------------------------------------

   // validate vector sizes

   std::stringstream message;
   message << methodName << std::endl;
   int errors  = 0;
   if( scalarNameI.size() != scalarI.size() ){
      message << "   int name size=" << scalarNameI.size() << " != " << " int vector size=" << scalarI.size() << std::endl; 
      errors++;
   }
   if( scalarNameR.size() != scalarR.size() ){
      message << "   real name size=" << scalarNameR.size() << " != " << " real vector size=" << scalarR.size() << std::endl; 
      errors++;
   }

   if( scalarNameR1.size() != scalarR1.size() ){
      message << "   real* name size=" << scalarNameR1.size() << " != " << " real* vector size=" << scalarR1.size() << std::endl; 
      errors++;
   }

   if( scalarNameR2.size() != scalarR2.size() ){
      message << "   real** name size=" << scalarNameR2.size() << " != " << " real** vector size=" << scalarR2.size() << std::endl; 
      errors++;
   }
   if( errors ){
      SimTKOpenMMLog::printError( message );
      return ReferenceDynamics::ErrorReturn;
   }

   // ---------------------------------------------------------------------------------------

   // array header

   (void) fprintf( stateFile, "# " );
   for( unsigned int ii = 0; ii < scalarNameR1.size(); ii++ ){  
      (void) fprintf( stateFile, " %s", scalarNameR1[ii].c_str() );
   }

   for( unsigned int ii = 0; ii < scalarNameR2.size(); ii++ ){  
      (void) fprintf( stateFile, " %s", scalarNameR2[ii].c_str() );
   }
   (void) fprintf( stateFile, "\n" );

   // int scalars

   for( unsigned int ii = 0; ii < scalarI.size(); ii++ ){  
      (void) fprintf( stateFile, "%6d  # %s\n", scalarI[ii], scalarNameI[ii].c_str() );
   }

   // real scalars

   for( unsigned int ii = 0; ii < scalarR.size(); ii++ ){  
      (void) fprintf( stateFile, "%14.6e  # %s\n", scalarR[ii], scalarNameR[ii].c_str() );
   }

   // arrays

   int maxPerLine = 3;
   for( unsigned int ii = 0; ii < (unsigned int) dimension1; ii++ ){  
      (void) fprintf( stateFile, "%5d ", ii );
      for( RealOpenMMPtrVectorI jj = scalarR1.begin(); jj != scalarR1.end(); jj++ ){  
         (void) fprintf( stateFile, "%14.6e ", (*jj)[ii] );
      }
      int r2Count = 0;
      for( RealOpenMMPtrPtrVectorI jj = scalarR2.begin(); jj != scalarR2.end(); jj++ ){  
         RealOpenMM** array = *jj;
         for( int kk = 0; kk < dimension2; kk++ ){
            (void) fprintf( stateFile, "%14.6e ", array[ii][kk] );
         }
         (void) fprintf( stateFile, "   " );
         r2Count++;
         if( r2Count == maxPerLine && ( (r2Count % (int) (scalarR2.size()-1)) > 2) ){
            (void) fprintf( stateFile, "\n     " );
            r2Count = 0;
         }
      }
      (void) fprintf( stateFile, "\n" );
   }
         
   return ReferenceDynamics::DefaultReturn;
}
