
/* Portions copyright (c) 2006-2023 Stanford University and Simbios.
 * Contributors: Peter Eastman, Pande Group
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

#ifndef __ReferenceVariableVerletDynamics_H__
#define __ReferenceVariableVerletDynamics_H__

#include "ReferenceDynamics.h"

namespace OpenMM {

class ReferenceVariableVerletDynamics : public ReferenceDynamics {

   private:

      std::vector<OpenMM::Vec3> xPrime;
      std::vector<double> inverseMasses;
      double _accuracy;

   public:

      /**---------------------------------------------------------------------------------------

         Constructor

         @param numberOfAtoms  number of atoms
         @param accuracy       required accuracy

         --------------------------------------------------------------------------------------- */

       ReferenceVariableVerletDynamics(int numberOfAtoms, double accuracy);

      /**---------------------------------------------------------------------------------------

         Destructor

         --------------------------------------------------------------------------------------- */

       ~ReferenceVariableVerletDynamics();

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
         @param boxVectors          the current periodic box vectors

         --------------------------------------------------------------------------------------- */

      void update(const OpenMM::System& system, std::vector<OpenMM::Vec3>& atomCoordinates,
                  std::vector<OpenMM::Vec3>& velocities, std::vector<OpenMM::Vec3>& forces, std::vector<double>& masses, double maxStepSize, double tolerance, const Vec3* boxVectors);

};

} // namespace OpenMM

#endif // __ReferenceVariableVerletDynamics_H__
