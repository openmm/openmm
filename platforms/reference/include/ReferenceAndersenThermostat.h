
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

#ifndef __ReferenceAndersenThermostat_H__
#define __ReferenceAndersenThermostat_H__

#include "openmm/Vec3.h"
#include <vector>

namespace OpenMM {

class ReferenceAndersenThermostat {

   private:
       
   public:

      /**---------------------------------------------------------------------------------------
      
         Constructor
      
         --------------------------------------------------------------------------------------- */

       ReferenceAndersenThermostat();

      /**---------------------------------------------------------------------------------------
      
         Destructor
      
         --------------------------------------------------------------------------------------- */

       ~ReferenceAndersenThermostat();

      /**---------------------------------------------------------------------------------------
      
         Apply the thermostat at the start of a time step.
      
         @param atomGroups         the groups of atoms to apply the thermostat to
         @param atomVelocities     atom velocities
         @param atomMasses         atom masses
         @param temperature        thermostat temperature in Kelvin
         @param collisionFrequency collision frequency for each atom in fs^-1
         @param stepSize           integration step size in fs
                  
         --------------------------------------------------------------------------------------- */
          
      void applyThermostat(const std::vector<std::vector<int> >& atomGroups, std::vector<OpenMM::Vec3>& atomVelocities, std::vector<double>& atomMasses,
              double temperature, double collisionFrequency, double stepSize) const;
      
};

} // namespace OpenMM

#endif // __ReferenceAndersenThermostat_H__
