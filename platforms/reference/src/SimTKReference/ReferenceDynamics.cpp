
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

#include <cstring>
#include <sstream>

#include "../SimTKUtilities/SimTKOpenMMCommon.h"
#include "../SimTKUtilities/SimTKOpenMMLog.h"
#include "../SimTKUtilities/SimTKOpenMMUtilities.h"
#include "ReferenceDynamics.h"

#include <cstdio>

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
   }

   if( totalMass > zero ){

      RealOpenMM totalMassI      = one/totalMass;
      linearMomentum[0]         *= totalMassI;
      linearMomentum[1]         *= totalMassI;
      linearMomentum[2]         *= totalMassI;

      for( int ii = 0; ii < numberOfAtoms; ii++ ){
         velocities[ii][0]     -= linearMomentum[0];
         velocities[ii][1]     -= linearMomentum[1];
         velocities[ii][2]     -= linearMomentum[2];
      }
   }

   return ReferenceDynamics::DefaultReturn;
}
