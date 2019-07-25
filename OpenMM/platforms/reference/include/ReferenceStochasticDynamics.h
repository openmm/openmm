
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

#ifndef __ReferenceStochasticDynamics_H__
#define __ReferenceStochasticDynamics_H__

#include "ReferenceDynamics.h"
#include "openmm/internal/windowsExport.h"

namespace OpenMM {

class OPENMM_EXPORT ReferenceStochasticDynamics : public ReferenceDynamics {

   protected:

      std::vector<OpenMM::RealVec> xPrime;
      std::vector<RealOpenMM> inverseMasses;
      RealOpenMM friction;
      
   public:

      /**---------------------------------------------------------------------------------------
      
         Constructor

         @param numberOfAtoms  number of atoms
         @param deltaT         delta t for dynamics
         @param friction       friction coefficient
         @param temperature    temperature
      
         --------------------------------------------------------------------------------------- */

       ReferenceStochasticDynamics(int numberOfAtoms, RealOpenMM deltaT, RealOpenMM friction, RealOpenMM temperature);

      /**---------------------------------------------------------------------------------------
      
         Destructor
      
         --------------------------------------------------------------------------------------- */

       ~ReferenceStochasticDynamics();

      /**---------------------------------------------------------------------------------------
      
         Get friction coefficient
      
         --------------------------------------------------------------------------------------- */
      
      RealOpenMM getFriction() const;
      
      /**---------------------------------------------------------------------------------------
      
         Update
      
         @param system              the System to be integrated
         @param atomCoordinates     atom coordinates
         @param velocities          velocities
         @param forces              forces
         @param masses              atom masses
         @param tolerance           the constraint tolerance
      
         --------------------------------------------------------------------------------------- */
     
      void update(const OpenMM::System& system, std::vector<OpenMM::RealVec>& atomCoordinates,
                  std::vector<OpenMM::RealVec>& velocities, std::vector<OpenMM::RealVec>& forces, std::vector<RealOpenMM>& masses, RealOpenMM tolerance);
     
      /**---------------------------------------------------------------------------------------
      
         First update; based on code in update.c do_update_sd() Gromacs 3.1.4
      
         @param numberOfAtoms       number of atoms
         @param atomCoordinates     atom coordinates
         @param velocities          velocities
         @param forces              forces
         @param inverseMasses       inverse atom masses
         @param xPrime              xPrime
      
         --------------------------------------------------------------------------------------- */
      
      virtual void updatePart1(int numberOfAtoms, std::vector<OpenMM::RealVec>& atomCoordinates, std::vector<OpenMM::RealVec>& velocities,
                       std::vector<OpenMM::RealVec>& forces, std::vector<RealOpenMM>& inverseMasses, std::vector<OpenMM::RealVec>& xPrime);
      
      /**---------------------------------------------------------------------------------------
      
         Second update
      
         @param numberOfAtoms       number of atoms
         @param atomCoordinates     atom coordinates
         @param velocities          velocities
         @param forces              forces
         @param inverseMasses       inverse atom masses
      
         --------------------------------------------------------------------------------------- */
      
      virtual void updatePart2(int numberOfAtoms, std::vector<OpenMM::RealVec>& atomCoordinates, std::vector<OpenMM::RealVec>& velocities,
                       std::vector<OpenMM::RealVec>& forces, std::vector<RealOpenMM>& inverseMasses, std::vector<OpenMM::RealVec>& xPrime);
      
      /**---------------------------------------------------------------------------------------
      
         Third update
      
         @param numberOfAtoms       number of atoms
         @param atomCoordinates     atom coordinates
         @param velocities          velocities
         @param inverseMasses       inverse atom masses
      
         --------------------------------------------------------------------------------------- */
      
      virtual void updatePart3(int numberOfAtoms, std::vector<OpenMM::RealVec>& atomCoordinates, std::vector<OpenMM::RealVec>& velocities,
                       std::vector<RealOpenMM>& inverseMasses, std::vector<OpenMM::RealVec>& xPrime);
};

} // namespace OpenMM

#endif // __ReferenceStochasticDynamics_H__
