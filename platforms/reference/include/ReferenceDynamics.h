
/* Portions copyright (c) 2006-2012 Stanford University and Simbios.
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
#include "SimTKOpenMMCommon.h"
#include "openmm/System.h"
#include <cstddef>
#include <vector>

// ---------------------------------------------------------------------------------------

/**---------------------------------------------------------------------------------------

Abstract class for dynamics

Main method (virtual) is update()

'Random' numbers are currently fixed to allow testing

--------------------------------------------------------------------------------------- */

class OPENMM_EXPORT ReferenceDynamics {

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
      
         Update
      
         @param system              the System to be integrated
         @param atomCoordinates     atom coordinates
         @param velocities          velocities
         @param forces              forces
         @param masses              atom masses
         @param tolerance           the constraint tolerance
      
         --------------------------------------------------------------------------------------- */
      
      virtual void update(const OpenMM::System& system, std::vector<OpenMM::RealVec>& atomCoordinates,
                          std::vector<OpenMM::RealVec>& velocities, std::vector<OpenMM::RealVec>& forces, std::vector<RealOpenMM>& masses, RealOpenMM tolerance);

      /**---------------------------------------------------------------------------------------
      
         Get ReferenceConstraint
      
         @return referenceConstraint  object
      
         --------------------------------------------------------------------------------------- */
      
      ReferenceConstraintAlgorithm* getReferenceConstraintAlgorithm( void ) const;
      
      /**---------------------------------------------------------------------------------------
      
         Set ReferenceConstraint
      
         @param referenceConstraint  referenceConstraint object
      
         --------------------------------------------------------------------------------------- */
      
      void setReferenceConstraintAlgorithm( ReferenceConstraintAlgorithm* referenceConstraint );
};

// ---------------------------------------------------------------------------------------

#endif // __ReferenceDynamics_H__
