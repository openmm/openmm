
/* Portions copyright (c) 2010 Stanford University and Simbios.
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

#ifndef __ReferenceMonteCarloBarostat_H__
#define __ReferenceMonteCarloBarostat_H__

#include "openmm/Vec3.h"
#include <utility>
#include <vector>

namespace OpenMM {

class ReferenceMonteCarloBarostat {

   private:

       std::vector<double> savedAtomPositions[3];
       std::vector<std::vector<int> > molecules;

   public:

      /**---------------------------------------------------------------------------------------

         Constructor

         --------------------------------------------------------------------------------------- */

       ReferenceMonteCarloBarostat(int numAtoms, const std::vector<std::vector<int> >& molecules);

      /**---------------------------------------------------------------------------------------

         Destructor

         --------------------------------------------------------------------------------------- */

       ~ReferenceMonteCarloBarostat();

      /**---------------------------------------------------------------------------------------

         Apply the barostat at the start of a time step, scaling x, y, and z coordinates independently.

         @param atomPositions      atom positions
         @param boxVectors         the periodic box vectors
         @param scaleX             the factor by which to scale atomic x coordinates
         @param scaleY             the factor by which to scale atomic y coordinates
         @param scaleZ             the factor by which to scale atomic z coordinates

         --------------------------------------------------------------------------------------- */

      void applyBarostat(std::vector<OpenMM::Vec3>& atomPositions, const OpenMM::Vec3* boxVectors, double scaleX, double scaleY, double scaleZ);

      /**---------------------------------------------------------------------------------------

         Restore atom positions to what they were before applyBarostat() was called.

         @param atomPositions      atom positions

         --------------------------------------------------------------------------------------- */

      void restorePositions(std::vector<OpenMM::Vec3>& atomPositions);

};

} // namespace OpenMM

#endif // __ReferenceMonteCarloBarostat_H__
