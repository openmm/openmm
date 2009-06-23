
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

#ifndef __ReferenceDynamics_H__
#define __ReferenceDynamics_H__

#include "ReferenceConstraintAlgorithm.h"
#include "../SimTKUtilities/SimTKOpenMMCommon.h"
#include <cstddef>

// ---------------------------------------------------------------------------------------

/**---------------------------------------------------------------------------------------

Abstract class for dynamics

Main method (virtual) is update()

'Random' numbers are currently fixed to allow testing

--------------------------------------------------------------------------------------- */

class ReferenceDynamics {

   public:

      /**---------------------------------------------------------------------------------------
      
         Fixed static variables
      
         --------------------------------------------------------------------------------------- */

       static const int DefaultReturn;
       static const int ErrorReturn;

   private:

      int _numberOfAtoms;
      int _timeStep;

      RealOpenMM _deltaT;
      RealOpenMM _temperature;

      int _numberOf2DTempArrays;
      RealOpenMM*** _twoDTempArrays;

      int _numberOf1DTempArrays;
      RealOpenMM** _oneDTempArrays;

      int _ownReferenceConstraint;
      ReferenceConstraintAlgorithm* _referenceConstraint;

      /**---------------------------------------------------------------------------------------
      
         Free memory associated w/ 2D arrays
      
         @return ReferenceDynamics::DefaultReturn
      
         --------------------------------------------------------------------------------------- */
      
      int _freeTwoDArrays( void );
      
      /**---------------------------------------------------------------------------------------
      
         Free memory associated w/ 1D arrays
      
         @return ReferenceDynamics::DefaultReturn
      
         --------------------------------------------------------------------------------------- */
      
      int _freeOneDArrays( void );
      
   protected:

      /**---------------------------------------------------------------------------------------
      
         Allocate memory associated w/ 2D arrays
      
         @param dimension1        first dimension
         @param dimension2        second dimension
         @param numberOfArrays    number of arrays to allocate
      
         @return ReferenceDynamics::DefaultReturn
      
         --------------------------------------------------------------------------------------- */
      
      int allocate2DArrays( int dimension1, int dimension2, int numberOfArrays );
      
      /**---------------------------------------------------------------------------------------
      
         Get array at specified index
      
         @param index             array index
      
         @return array or NULL if index invalid or arrays not allocated
      
         --------------------------------------------------------------------------------------- */
      
      RealOpenMM** get2DArrayAtIndex( int index ) const;
      
      /**---------------------------------------------------------------------------------------
      
         Allocate memory associated w/ 1D arrays
      
         @param dimension1        dimension
         @param numberOfArrays    number of arrays to allocate
      
         @return DefaultReturn
      
         --------------------------------------------------------------------------------------- */
      
      int allocate1DArrays( int dimension, int numberOfArrays );
      
      /**---------------------------------------------------------------------------------------
      
         Get array at specified index
      
         @param index             array index
      
         @return array or NULL if index invalid or arrays not allocated
      
         --------------------------------------------------------------------------------------- */
      
      RealOpenMM* get1DArrayAtIndex( int index ) const;
      
   public:

      /**---------------------------------------------------------------------------------------
      
         Constructor
      
         @param numberOfAtoms  number of atoms
         @param deltaT         delta t for dynamics
         @param temperature    temperature

         --------------------------------------------------------------------------------------- */

       ReferenceDynamics( int numberOfAtoms, RealOpenMM _deltaT, RealOpenMM temperature );

      /**---------------------------------------------------------------------------------------
      
         Destructor
      
         --------------------------------------------------------------------------------------- */

       virtual ~ReferenceDynamics( );

      /**---------------------------------------------------------------------------------------
      
         Get number of atoms
      
         @return number of atoms
      
         --------------------------------------------------------------------------------------- */
      
      int getNumberOfAtoms( void ) const;
      
      /**---------------------------------------------------------------------------------------
      
         Get time step
      
         @return time step
      
         --------------------------------------------------------------------------------------- */
      
      int getTimeStep( void ) const;
      
      /**---------------------------------------------------------------------------------------
      
         Increment time step
      
         @return incremented time step
      
         --------------------------------------------------------------------------------------- */
      
      int incrementTimeStep( void );
      
      /**---------------------------------------------------------------------------------------
      
         Get delta t
      
         @return deltaT
      
         --------------------------------------------------------------------------------------- */
      
      RealOpenMM getDeltaT( void ) const;

      /**---------------------------------------------------------------------------------------

         Set delta t

         --------------------------------------------------------------------------------------- */

      void setDeltaT( RealOpenMM deltaT );

      /**---------------------------------------------------------------------------------------
      
         Get temperature
      
         @return temperature
      
         --------------------------------------------------------------------------------------- */
    
      RealOpenMM getTemperature( void ) const;
    
      /**---------------------------------------------------------------------------------------
      
         Remove total linear momentum
      
         @param numberOfAtoms      number of atoms
         @param masses             masses
         @param velocities         velocities
      
         @return ReferenceDynamics::DefaultReturn
      
         --------------------------------------------------------------------------------------- */
      
      int removeTotalLinearMomentum( int numberOfAtoms, RealOpenMM* masses, RealOpenMM** velocities ) const;
      
      /**---------------------------------------------------------------------------------------
      
         Print parameters
      
         @param message message

         @return ReferenceDynamics::DefaultReturn
      
         --------------------------------------------------------------------------------------- */
          
      int printParameters( std::stringstream& message ) const;

      /**---------------------------------------------------------------------------------------
      
         Update
      
         @param numberOfAtoms       number of atoms
         @param atomCoordinates     atom coordinates
         @param velocities          velocities
         @param forces              forces
         @param masses              atom masses
      
         @return ReferenceDynamics::DefaultReturn
      
         --------------------------------------------------------------------------------------- */
      
      virtual int update( int numberOfAtoms, RealOpenMM** atomCoordinates,
                          RealOpenMM** velocities, RealOpenMM** forces, RealOpenMM* masses );

      /**---------------------------------------------------------------------------------------
      
         Get ReferenceConstraint
      
         @return referenceConstraint  object
      
         --------------------------------------------------------------------------------------- */
      
      ReferenceConstraintAlgorithm* getReferenceConstraintAlgorithm( void ) const;
      
      /**---------------------------------------------------------------------------------------
      
         Set ReferenceConstraint
      
         @param referenceConstraint  referenceConstraint object
      
         @return ReferenceDynamics::DefaultReturn
      
         --------------------------------------------------------------------------------------- */
      
      int setReferenceConstraintAlgorithm( ReferenceConstraintAlgorithm* referenceConstraint );
      
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
    
      int writeState( int numberOfAtoms, RealOpenMM** atomCoordinates,
                      RealOpenMM** velocities, RealOpenMM** forces, RealOpenMM* masses,
                      int state, const std::string& baseFileName ) const;

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
      
      int writeStateToFile( FILE* stateFile, StringVector& scalarNameI, IntVector& scalarI,
                            StringVector& scalarNameR, RealOpenMMVector& scalarR,
                            int dimension1, StringVector& scalarNameR1, RealOpenMMPtrVector& scalarR1,
                            int dimension2, StringVector& scalarNameR2, RealOpenMMPtrPtrVector& scalarR2 ) const;
      
};

// ---------------------------------------------------------------------------------------

#endif // __ReferenceDynamics_H__
