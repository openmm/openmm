
/* Portions copyright (c) 2006-2016 Stanford University and Simbios.
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

#ifndef __ReferenceVariableStochasticDynamics_H__
#define __ReferenceVariableStochasticDynamics_H__

#include "ReferenceDynamics.h"

namespace OpenMM {

class ReferenceVariableStochasticDynamics : public ReferenceDynamics {

   private:

      std::vector<OpenMM::Vec3> xPrime;
      std::vector<double> inverseMasses;
      double friction, _accuracy;

   public:

      /**---------------------------------------------------------------------------------------

         Constructor

         @param numberOfAtoms  number of atoms
         @param friction       friction coefficient
         @param temperature    temperature
         @param accuracy       required accuracy

         --------------------------------------------------------------------------------------- */

       ReferenceVariableStochasticDynamics(int numberOfAtoms, double friction, double temperature, double accuracy);

      /**---------------------------------------------------------------------------------------

         Destructor

         --------------------------------------------------------------------------------------- */

       ~ReferenceVariableStochasticDynamics();

      /**---------------------------------------------------------------------------------------

         Get friction coefficient

         --------------------------------------------------------------------------------------- */

      double getFriction() const;
      
      /**---------------------------------------------------------------------------------------

         Get the required accuracy

         @return accuracy

         --------------------------------------------------------------------------------------- */

      double getAccuracy() const;

      /**---------------------------------------------------------------------------------------

         Set the required accuracy

         --------------------------------------------------------------------------------------- */

      void setAccuracy(double accuracy);

      /**---------------------------------------------------------------------------------------

         Update

         @param system              the System to be integrated
         @param atomCoordinates     atom coordinates
         @param velocities          velocities
         @param forces              forces
         @param masses              atom masses
         @param maxStepSize         maximum time step
         @param tolerance           the constraint tolerance

         --------------------------------------------------------------------------------------- */

      void update(const OpenMM::System& system, std::vector<OpenMM::Vec3>& atomCoordinates,
                  std::vector<OpenMM::Vec3>& velocities, std::vector<OpenMM::Vec3>& forces, std::vector<double>& masses, double maxStepSize, double tolerance);

      /**---------------------------------------------------------------------------------------

         First update; based on code in update.c do_update_sd() Gromacs 3.1.4

         @param numberOfAtoms       number of atoms
         @param atomCoordinates     atom coordinates
         @param velocities          velocities
         @param forces              forces
         @param masses              atom masses
         @param inverseMasses       inverse atom masses
         @param xPrime              xPrime
         @param maxStepSize         maximum time step

         --------------------------------------------------------------------------------------- */

      void updatePart1(int numberOfAtoms, std::vector<OpenMM::Vec3>& atomCoordinates, std::vector<OpenMM::Vec3>& velocities,
                       std::vector<OpenMM::Vec3>& forces, std::vector<double>& masses, std::vector<double>& inverseMasses,
                       std::vector<OpenMM::Vec3>& xPrime, double maxStepSize);

      /**---------------------------------------------------------------------------------------

         Second update

         @param numberOfAtoms       number of atoms
         @param atomCoordinates     atom coordinates
         @param velocities          velocities
         @param forces              forces
         @param masses              atom masses

         --------------------------------------------------------------------------------------- */

      void updatePart2(int numberOfAtoms, std::vector<OpenMM::Vec3>& atomCoordinates, std::vector<OpenMM::Vec3>& velocities,
                       std::vector<OpenMM::Vec3>& forces, std::vector<double>& inverseMasses,
                       std::vector<OpenMM::Vec3>& xPrime);
      
};

} // namespace OpenMM

#endif // __ReferenceVariableStochasticDynamics_H__
