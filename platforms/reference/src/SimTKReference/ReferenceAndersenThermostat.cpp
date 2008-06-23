
/* Portions copyright (c) 2008 Stanford University and Simbios.
 * Contributors: Peter Eastman
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

#include "../SimTKUtilities/SimTKOpenMMUtilities.h"
#include "ReferenceAndersenThermostat.h"

      /**---------------------------------------------------------------------------------------
      
         Constructor
      
         --------------------------------------------------------------------------------------- */

       ReferenceAndersenThermostat::ReferenceAndersenThermostat( ) {
       }

      /**---------------------------------------------------------------------------------------
      
         Destructor
      
         --------------------------------------------------------------------------------------- */

       ReferenceAndersenThermostat::~ReferenceAndersenThermostat( ) {
       }

      /**---------------------------------------------------------------------------------------
      
         Apply the thermostat at the start of a time step.
      
         @param numberOfAtoms      number of atoms
         @param atomVelocities     atom velocities
         @param temperature        thermostat temperature in Kelvin
         @param collisionFrequency collision frequency for each atom in fs^-1
         @param stepSize           integration step size in fs
                  
         --------------------------------------------------------------------------------------- */
          
      void ReferenceAndersenThermostat::applyThermostat( int numberOfAtoms, RealOpenMM** atomVelocities, RealOpenMM* atomMasses,
              RealOpenMM temperature, RealOpenMM collisionFrequency, RealOpenMM stepSize ) const {
          
          const RealOpenMM collisionProbability = collisionFrequency*stepSize;
          for (int i = 0; i < numberOfAtoms; ++i) {
              if (SimTKOpenMMUtilities::getUniformlyDistributedRandomNumber() < collisionProbability) {
                  
                  // A collision occurred, so set the velocity to a new value chosen from a Boltzmann distribution.
                  
                  const RealOpenMM velocityScale = static_cast<RealOpenMM>( sqrt(BOLTZ*temperature/atomMasses[i]) );
                  atomVelocities[i][0] = velocityScale*SimTKOpenMMUtilities::getNormallyDistributedRandomNumber();
                  atomVelocities[i][1] = velocityScale*SimTKOpenMMUtilities::getNormallyDistributedRandomNumber();
                  atomVelocities[i][2] = velocityScale*SimTKOpenMMUtilities::getNormallyDistributedRandomNumber();
              }
          }
          
      }

