
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

#ifndef __ReferenceDPMECoulombIxn_H__
#define __ReferenceDPMECoulombIxn_H__

#include "ReferencePairIxn.h"
#include "ReferenceNeighborList.h"

namespace OpenMM {

class ReferenceDPMECoulombIxn {

   private:
       
      bool cutoff;
      bool useSwitch;
      bool periodic;
      bool ewald;
      bool pme_setup, dpme_setup;
      const OpenMM::NeighborList* neighborList;
      OpenMM::RealVec periodicBoxVectors[3];
      RealOpenMM cutoffDistance, dpmeShift;
      RealOpenMM krf, crf;
      RealOpenMM alphaEwald, alphaDispersionEwald;
      int numRx, numRy, numRz;
      int meshDim[3], dispersionMeshDim[3];

      // parameter indices

      static const int SigIndex = 0;
      static const int EpsIndex = 1;
      static const int   QIndex = 2;


   public:

      /**---------------------------------------------------------------------------------------
      
         Constructor
      
         --------------------------------------------------------------------------------------- */

       ReferenceDPMECoulombIxn();

      /**---------------------------------------------------------------------------------------
      
         Destructor
      
         --------------------------------------------------------------------------------------- */

       ~ReferenceDPMECoulombIxn();

      /**---------------------------------------------------------------------------------------
      
         Set the force to use a cutoff.
      
         @param distance            the cutoff distance
         @param neighbors           the neighbor list to use
         @param solventDielectric   the dielectric constant of the bulk solvent
      
         --------------------------------------------------------------------------------------- */
      
      void setUseCutoff(RealOpenMM distance, const OpenMM::NeighborList& neighbors, RealOpenMM solventDielectric);

      
      /**---------------------------------------------------------------------------------------
      
         Set the force to use periodic boundary conditions.  This requires that a cutoff has
         already been set, and the smallest side of the periodic box is at least twice the cutoff
         distance.
      
         @param vectors    the vectors defining the periodic box
      
         --------------------------------------------------------------------------------------- */
      
      void setPeriodic(OpenMM::RealVec* vectors);
       
     
      /**---------------------------------------------------------------------------------------

         Set the force to use Particle-Mesh Ewald (PME) summation for dispersion.

         @param dalpha             the dispersion Ewald separation parameter
         @param dispersionGridSize the dimensions of the disperson mesh

         --------------------------------------------------------------------------------------- */

      void setupDispersionPME(RealOpenMM dalpha, int dispersionMeshSize[3]);


      /**---------------------------------------------------------------------------------------
      
         Set the force to use Particle-Mesh Ewald (PME) summation.
      
         @param alpha    the Ewald separation parameter
         @param gridSize the dimensions of the mesh
      
         --------------------------------------------------------------------------------------- */
      
      void setupPME(RealOpenMM alpha, int meshSize[3]);
      
      /**---------------------------------------------------------------------------------------
      
         Calculate DPME Coulomb pair ixn
      
         @param numberOfAtoms    number of atoms
         @param atomCoordinates  atom coordinates
         @param atomParameters   atom parameters (charges, c6, c12, ...)     atomParameters[atomIndex][paramterIndex]
         @param exclusions       atom exclusion indices
                                 exclusions[atomIndex] contains the list of exclusions for that atom
         @param fixedParameters  non atom parameters (not currently used)
         @param forces           force array (forces added)
         @param energyByAtom     atom energy
         @param totalEnergy      total energy
         @param includeDirect      true if direct space interactions should be included
         @param includeReciprocal  true if reciprocal space interactions should be included
      
         --------------------------------------------------------------------------------------- */
          
      void calculatePairIxn(int numberOfAtoms, std::vector<OpenMM::RealVec>& atomCoordinates,
                            RealOpenMM** atomParameters, std::vector<std::set<int> >& exclusions,
                            RealOpenMM* fixedParameters, std::vector<OpenMM::RealVec>& forces,
                            RealOpenMM* energyByAtom, RealOpenMM* totalEnergy, bool includeDirect, bool includeReciprocal) const;

private:
      /**---------------------------------------------------------------------------------------
      
         Calculate Ewald ixn
      
         @param numberOfAtoms    number of atoms
         @param atomCoordinates  atom coordinates
         @param atomParameters   atom parameters (charges, c6, c12, ...)     atomParameters[atomIndex][paramterIndex]
         @param exclusions       atom exclusion indices
                                 exclusions[atomIndex] contains the list of exclusions for that atom
         @param fixedParameters  non atom parameters (not currently used)
         @param forces           force array (forces added)
         @param energyByAtom     atom energy
         @param totalEnergy      total energy
         @param includeDirect      true if direct space interactions should be included
         @param includeReciprocal  true if reciprocal space interactions should be included
            
         --------------------------------------------------------------------------------------- */
          
      void calculatePMEIxn(int numberOfAtoms, std::vector<OpenMM::RealVec>& atomCoordinates,
                            RealOpenMM** atomParameters, std::vector<std::set<int> >& exclusions,
                            RealOpenMM* fixedParameters, std::vector<OpenMM::RealVec>& forces,
                            RealOpenMM* energyByAtom, RealOpenMM* totalEnergy, bool includeDirect, bool includeReciprocal) const;
};

} // namespace OpenMM

#endif // __ReferenceDPMECoulombIxn_H__
