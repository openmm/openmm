
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
#include <vector>

// ---------------------------------------------------------------------------------------

/**---------------------------------------------------------------------------------------

Abstract class for dynamics

Main method (virtual) is update()

'Random' numbers are currently fixed to allow testing

--------------------------------------------------------------------------------------- */

class OPENMM_EXPORT ReferenceDynamics {

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

      int _ownReferenceConstraint;
      ReferenceConstraintAlgorithm* _referenceConstraint;
      
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
      
      int removeTotalLinearMomentum( int numberOfAtoms, RealOpenMM* masses, std::vector<OpenMM::RealVec>& velocities ) const;

      /**---------------------------------------------------------------------------------------
      
         Update
      
         @param numberOfAtoms       number of atoms
         @param atomCoordinates     atom coordinates
         @param velocities          velocities
         @param forces              forces
         @param masses              atom masses
      
         @return ReferenceDynamics::DefaultReturn
      
         --------------------------------------------------------------------------------------- */
      
      virtual int update( int numberOfAtoms, std::vector<OpenMM::RealVec>& atomCoordinates,
                          std::vector<OpenMM::RealVec>& velocities, std::vector<OpenMM::RealVec>& forces, std::vector<RealOpenMM>& masses );

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
};

// ---------------------------------------------------------------------------------------

#endif // __ReferenceDynamics_H__
