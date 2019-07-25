/* Portions copyright (c) 2006-2013 Stanford University and Simbios.
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

#ifndef __ReferenceForce_H__
#define __ReferenceForce_H__

#include "RealVec.h"
#include "lepton/CompiledExpression.h"
#include "openmm/internal/windowsExport.h"

namespace OpenMM {

class OPENMM_EXPORT  ReferenceForce {

   private:
       
       static RealOpenMM periodicDifference(RealOpenMM val1, RealOpenMM val2, RealOpenMM period);

   public:

      /**---------------------------------------------------------------------------------------
      
         Constructor
      
         --------------------------------------------------------------------------------------- */

       ReferenceForce();

      /**---------------------------------------------------------------------------------------
      
         Destructor
      
         --------------------------------------------------------------------------------------- */

       ~ReferenceForce();

      /**---------------------------------------------------------------------------------------
      
         Static variables
      
         --------------------------------------------------------------------------------------- */

       static const int XIndex             = 0;
       static const int YIndex             = 1;
       static const int ZIndex             = 2;
       static const int R2Index            = 3;
       static const int RIndex             = 4;
       static const int LastDeltaRIndex    = 5;

       static const int DeltaRMaxIndex     = 5;
   
      /**---------------------------------------------------------------------------------------
      
         Get deltaR and distance and distance**2 between atomI and atomJ (static method)
         deltaR: j - i
      
         @param atomCoordinatesI    atom i coordinates
         @param atomCoordinatesI    atom j coordinates
         @param deltaR              deltaX, deltaY, deltaZ, R2, R upon return
      
         --------------------------------------------------------------------------------------- */
      
      static void getDeltaR(const OpenMM::RealVec& atomCoordinatesI, const OpenMM::RealVec& atomCoordinatesJ,
                            RealOpenMM* deltaR);
      
      /**---------------------------------------------------------------------------------------
      
         Get deltaR between atomI and atomJ (static method)
         deltaR: j - i
      
         @param atomCoordinatesI    atom i coordinates
         @param atomCoordinatesI    atom j coordinates
         @return the displacement
      
         --------------------------------------------------------------------------------------- */
      
      static RealVec getDeltaR(const OpenMM::RealVec& atomCoordinatesI, const OpenMM::RealVec& atomCoordinatesJ);
      
      /**---------------------------------------------------------------------------------------

         Get deltaR and distance and distance**2 between atomI and atomJ, assuming periodic
         boundary conditions (static method); deltaR: j - i

         @param atomCoordinatesI    atom i coordinates
         @param atomCoordinatesI    atom j coordinates
         @param boxSize             X, Y, and Z sizes of the periodic box
         @param deltaR              deltaX, deltaY, deltaZ, R2, R upon return

         --------------------------------------------------------------------------------------- */

      static void getDeltaRPeriodic(const OpenMM::RealVec& atomCoordinatesI, const OpenMM::RealVec& atomCoordinatesJ,
                                    const RealOpenMM* boxSize, RealOpenMM* deltaR);

      /**---------------------------------------------------------------------------------------

         Get deltaR and distance and distance**2 between atomI and atomJ, assuming periodic
         boundary conditions (static method); deltaR: j - i

         @param atomCoordinatesI    atom i coordinates
         @param atomCoordinatesI    atom j coordinates
         @param boxVectors          the vectors defining the periodic box
         @param deltaR              deltaX, deltaY, deltaZ, R2, R upon return

         --------------------------------------------------------------------------------------- */

      static void getDeltaRPeriodic(const OpenMM::RealVec& atomCoordinatesI, const OpenMM::RealVec& atomCoordinatesJ,
                                    const OpenMM::RealVec* boxVectors, RealOpenMM* deltaR);

      /**---------------------------------------------------------------------------------------

         Get deltaR between atomI and atomJ, assuming periodic boundary conditions (static method);
         deltaR: j - i

         @param atomCoordinatesI    atom i coordinates
         @param atomCoordinatesI    atom j coordinates
         @param boxVectors          the vectors defining the periodic box
         @return the displacement

         --------------------------------------------------------------------------------------- */

      static RealVec getDeltaRPeriodic(const OpenMM::RealVec& atomCoordinatesI, const OpenMM::RealVec& atomCoordinatesJ,
                                    const OpenMM::RealVec* boxVectors);

      /**
       * Get a pointer to the memory for setting a variable in a CompiledExpression.  If the expression
       * does not use the specified variable, return NULL.
       */
      static double* getVariablePointer(Lepton::CompiledExpression& expression, const std::string& name);

      /**
       * Set the value of a variable in a CompiledExpression, using the pointer that was returned by getVariablePointer().
       */
      static void setVariable(double* pointer, double value);

};

} // namespace OpenMM

#endif // __ReferenceForce_H__
