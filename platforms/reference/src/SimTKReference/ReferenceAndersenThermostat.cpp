
/* Portions copyright (c) 2008-2010 Stanford University and Simbios.
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

#include <cmath>
#include <string.h>
#include <sstream>

#include "SimTKOpenMMUtilities.h"
#include "ReferenceAndersenThermostat.h"

using std::vector;
using namespace OpenMM;

      /**---------------------------------------------------------------------------------------
      
         Constructor
      
         --------------------------------------------------------------------------------------- */

       ReferenceAndersenThermostat::ReferenceAndersenThermostat() {
       }

      /**---------------------------------------------------------------------------------------
      
         Destructor
      
         --------------------------------------------------------------------------------------- */

       ReferenceAndersenThermostat::~ReferenceAndersenThermostat() {
       }

      /**---------------------------------------------------------------------------------------
      
         Apply the thermostat at the start of a time step.
      
         @param atomGroups         the groups of atoms to apply the thermostat to
         @param atomVelocities     atom velocities
         @param temperature        thermostat temperature in Kelvin
         @param collisionFrequency collision frequency for each atom in fs^-1
         @param stepSize           integration step size in fs
                  
         --------------------------------------------------------------------------------------- */
          
      void ReferenceAndersenThermostat::applyThermostat(const vector<vector<int> >& atomGroups, vector<Vec3>& atomVelocities, vector<double>& atomMasses,
              double temperature, double collisionFrequency, double stepSize) const {
          
          const double collisionProbability = 1.0f - exp(-collisionFrequency*stepSize);
          for (auto& group : atomGroups) {
              if (SimTKOpenMMUtilities::getUniformlyDistributedRandomNumber() < collisionProbability) {
                  
                  // A collision occurred, so set the velocities to new values chosen from a Boltzmann distribution.

                  for (int atom : group) {
                      if (atomMasses[atom] != 0) {
                          const double velocityScale = static_cast<double>(sqrt(BOLTZ*temperature/atomMasses[atom]));
                          atomVelocities[atom][0] = velocityScale*SimTKOpenMMUtilities::getNormallyDistributedRandomNumber();
                          atomVelocities[atom][1] = velocityScale*SimTKOpenMMUtilities::getNormallyDistributedRandomNumber();
                          atomVelocities[atom][2] = velocityScale*SimTKOpenMMUtilities::getNormallyDistributedRandomNumber();
                      }
                  }
              }
          }
          
      }

