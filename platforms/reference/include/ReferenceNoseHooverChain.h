
/* Portions copyright (c) 2019 Stanford University and Simbios.
 * Contributors: Andreas Krämer and Andrew C. Simmonett
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

#ifndef __ReferenceNoseHooverChain_H__
#define __ReferenceNoseHooverChain_H__

#include "openmm/Vec3.h"
#include <vector>

namespace OpenMM {

using std::vector;

class ReferenceNoseHooverChain {

   private:
       
   public:

      /**---------------------------------------------------------------------------------------
      
         Constructor
      
         --------------------------------------------------------------------------------------- */

       ReferenceNoseHooverChain();

      /**---------------------------------------------------------------------------------------
      
         Destructor
      
         --------------------------------------------------------------------------------------- */

       ~ReferenceNoseHooverChain();

      /**---------------------------------------------------------------------------------------

         Propagate the Nose-Hoover chain a half timestep and find the appropriate velocity scaling

         @param kineticEnergy      the instantaneous kinetic energy of the particles being thermostated
         @param chainVelocities    the velocities of the chain's beads in nm / ps
         @param chainPositions     the positions of the chains's beads in nm
         @param numDOFs            the number of degrees of freedom in the system that this chain thermostats
         @param temperature        thermostat temperature in Kelvin
         @param collisionFrequency collision frequency for each atom in ps^-1
         @param timeStep           full integration step size in ps (this only propagates half way)
         @param numMTS             number of multi timestep increments used in the Trotter expansion
         @param YSWeights          vector of weights used in the Yoshida-Suzuki multi-timestepping.
         --------------------------------------------------------------------------------------- */
      double propagate(double kineticEnergy, vector<double>& chainVelocities,
                       vector<double>& chainPositions, int numDOFs,
                       double temperature, double collisionFrequency, double timeStep,
                       int numMTS, const vector<double>& YSWeights) const;
};

} // namespace OpenMM

#endif // __ReferenceNoseHooverChain_H__
