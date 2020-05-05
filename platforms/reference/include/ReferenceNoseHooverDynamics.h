
/* Portions copyright (c) 2006-2020 Stanford University and Simbios.
 * Contributors: Andy Simmonett, Peter Eastman, Pande Group
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

#ifndef __ReferenceNoseHooverDynamics_H__
#define __ReferenceNoseHooverDynamics_H__

#include "ReferenceDynamics.h"
#include <tuple>

namespace OpenMM {

class ContextImpl;

class ReferenceNoseHooverDynamics : public ReferenceDynamics {

   private:
      std::vector<OpenMM::Vec3> xPrime;
      std::vector<OpenMM::Vec3> oldx;
      std::vector<double> inverseMasses;
      int numberOfAtoms;

   public:

      /**---------------------------------------------------------------------------------------
      
         Constructor

         @param numberOfAtoms  number of atoms
         @param deltaT         delta t for dynamics
         @param friction       friction coefficient
         @param temperature    temperature
      
         --------------------------------------------------------------------------------------- */

       ReferenceNoseHooverDynamics(int numberOfAtoms, double deltaT);

      /**---------------------------------------------------------------------------------------
      
         Destructor
      
         --------------------------------------------------------------------------------------- */

       ~ReferenceNoseHooverDynamics();

      /**---------------------------------------------------------------------------------------
      
         Perform the first half of a step using the leapfrog LF-Middle scheme
      
         @param system              the System to be integrated
         @param atomCoordinates     atom coordinates
         @param velocities          velocities
         @param forces              forces
         @param masses              atom masses
         @param tolerance           the constraint tolerance
         @param forcesAreValid      whether the forces are valid (duh!)
         @param allAtoms            a list of all atoms not involved in a Drude-like pair
         @param allPairs            a list of all Drude-like pairs, and their KT values, in the system
         @param maxPairDistance     the maximum separation allowed for a Drude-like pair
      
         --------------------------------------------------------------------------------------- */
     
      void step1(OpenMM::ContextImpl &context, const OpenMM::System& system, std::vector<OpenMM::Vec3>& atomCoordinates,
                 std::vector<OpenMM::Vec3>& velocities, std::vector<OpenMM::Vec3>& forces, std::vector<double>& masses, double tolerance, bool &forcesAreValid,
                 const std::vector<int> & allAtoms, const std::vector<std::tuple<int, int, double>> & allPairs, double maxPairDistance);
      /**---------------------------------------------------------------------------------------
      
         Perform the second half of a step using the leapfrog LF-Middle scheme
      
         @param system              the System to be integrated
         @param atomCoordinates     atom coordinates
         @param velocities          velocities
         @param forces              forces
         @param masses              atom masses
         @param tolerance           the constraint tolerance
         @param forcesAreValid      whether the forces are valid (duh!)
         @param allAtoms            a list of all atoms not involved in a Drude-like pair
         @param allPairs            a list of all Drude-like pairs, and their KT values, in the system
         @param maxPairDistance     the maximum separation allowed for a Drude-like pair
      
         --------------------------------------------------------------------------------------- */
      void step2(OpenMM::ContextImpl &context, const OpenMM::System& system, std::vector<OpenMM::Vec3>& atomCoordinates,
                 std::vector<OpenMM::Vec3>& velocities, std::vector<OpenMM::Vec3>& forces, std::vector<double>& masses, double tolerance, bool &forcesAreValid,
                 const std::vector<int> & allAtoms, const std::vector<std::tuple<int, int, double>> & allPairs, double maxPairDistance);
      
};

} // namespace OpenMM

#endif // __ReferenceNoseHooverDynamics_H__
