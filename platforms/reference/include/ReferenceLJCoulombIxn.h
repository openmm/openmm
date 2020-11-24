
/* Portions copyright (c) 2006-2020 Stanford University and Simbios.
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

#ifndef __ReferenceLJCoulombIxn_H__
#define __ReferenceLJCoulombIxn_H__

#include "ReferencePairIxn.h"
#include "ReferenceNeighborList.h"

namespace OpenMM {

class ReferenceLJCoulombIxn {

   private:
       
      bool cutoff;
      bool useSwitch;
      bool periodic, periodicExceptions;
      bool ewald;
      bool pme, ljpme;
      const OpenMM::NeighborList* neighborList;
      OpenMM::Vec3 periodicBoxVectors[3];
      double cutoffDistance, switchingDistance;
      double krf, crf;
      double alphaEwald, alphaDispersionEwald;
      int numRx, numRy, numRz;
      int meshDim[3], dispersionMeshDim[3];

      // parameter indices

      static const int SigIndex = 0;
      static const int EpsIndex = 1;
      static const int   QIndex = 2;
            
      /**---------------------------------------------------------------------------------------
      
         Calculate LJ Coulomb pair ixn between two atoms
      
         @param atom1            the index of the first atom
         @param atom2            the index of the second atom
         @param atomCoordinates  atom coordinates
         @param atomParameters   atom parameters (charges, c6, c12, ...)     atomParameters[atomIndex][paramterIndex]
         @param forces           force array (forces added)
         @param totalEnergy      total energy
            
         --------------------------------------------------------------------------------------- */
          
      void calculateOneIxn(int atom1, int atom2, std::vector<OpenMM::Vec3>& atomCoordinates,
                           std::vector<std::vector<double> >& atomParameters, std::vector<OpenMM::Vec3>& forces,
                           double* totalEnergy) const;


   public:

      /**---------------------------------------------------------------------------------------
      
         Constructor
      
         --------------------------------------------------------------------------------------- */

       ReferenceLJCoulombIxn();

      /**---------------------------------------------------------------------------------------
      
         Destructor
      
         --------------------------------------------------------------------------------------- */

       ~ReferenceLJCoulombIxn();

      /**---------------------------------------------------------------------------------------
      
         Set the force to use a cutoff.
      
         @param distance            the cutoff distance
         @param neighbors           the neighbor list to use
         @param solventDielectric   the dielectric constant of the bulk solvent
      
         --------------------------------------------------------------------------------------- */
      
      void setUseCutoff(double distance, const OpenMM::NeighborList& neighbors, double solventDielectric);

      /**---------------------------------------------------------------------------------------
      
         Set the force to use a switching function on the Lennard-Jones interaction.
      
         @param distance            the switching distance
      
         --------------------------------------------------------------------------------------- */
      
      void setUseSwitchingFunction(double distance);
      
      /**---------------------------------------------------------------------------------------
      
         Set the force to use periodic boundary conditions.  This requires that a cutoff has
         already been set, and the smallest side of the periodic box is at least twice the cutoff
         distance.
      
         @param vectors    the vectors defining the periodic box
      
         --------------------------------------------------------------------------------------- */
      
      void setPeriodic(OpenMM::Vec3* vectors);
       
      /**---------------------------------------------------------------------------------------
      
         Set the force to use Ewald summation.
      
         @param alpha  the Ewald separation parameter
         @param kmaxx  the largest wave vector in the x direction
         @param kmaxy  the largest wave vector in the y direction
         @param kmaxz  the largest wave vector in the z direction
      
         --------------------------------------------------------------------------------------- */
      
      void setUseEwald(double alpha, int kmaxx, int kmaxy, int kmaxz);

     
      /**---------------------------------------------------------------------------------------

         Set the force to use Particle-Mesh Ewald (PME) summation.

         @param alpha    the Ewald separation parameter
         @param gridSize the dimensions of the mesh

         --------------------------------------------------------------------------------------- */
      
      void setUsePME(double alpha, int meshSize[3]);
      
      /**---------------------------------------------------------------------------------------

         Set the force to use Particle-Mesh Ewald (PME) summation for dispersion.

         @param dalpha    the dispersion Ewald separation parameter
         @param dgridSize the dimensions of the dispersion mesh

         --------------------------------------------------------------------------------------- */

      void setUseLJPME(double dalpha, int dmeshSize[3]);
      
      /**---------------------------------------------------------------------------------------

         Set whether exceptions use periodic boundary conditions.

         --------------------------------------------------------------------------------------- */

      void setPeriodicExceptions(bool periodic);

      /**---------------------------------------------------------------------------------------
      
         Calculate LJ Coulomb pair ixn
      
         @param numberOfAtoms    number of atoms
         @param atomCoordinates  atom coordinates
         @param atomParameters   atom parameters (charges, c6, c12, ...)     atomParameters[atomIndex][paramterIndex]
         @param exclusions       atom exclusion indices
                                 exclusions[atomIndex] contains the list of exclusions for that atom
         @param forces           force array (forces added)
         @param totalEnergy      total energy
         @param includeDirect      true if direct space interactions should be included
         @param includeReciprocal  true if reciprocal space interactions should be included
      
         --------------------------------------------------------------------------------------- */
          
      void calculatePairIxn(int numberOfAtoms, std::vector<OpenMM::Vec3>& atomCoordinates,
                            std::vector<std::vector<double> >& atomParameters, std::vector<std::set<int> >& exclusions,
                            std::vector<OpenMM::Vec3>& forces, double* totalEnergy, bool includeDirect, bool includeReciprocal) const;

private:
      /**---------------------------------------------------------------------------------------
      
         Calculate Ewald ixn
      
         @param numberOfAtoms    number of atoms
         @param atomCoordinates  atom coordinates
         @param atomParameters   atom parameters (charges, c6, c12, ...)     atomParameters[atomIndex][paramterIndex]
         @param exclusions       atom exclusion indices
                                 exclusions[atomIndex] contains the list of exclusions for that atom
         @param forces           force array (forces added)
         @param totalEnergy      total energy
         @param includeDirect      true if direct space interactions should be included
         @param includeReciprocal  true if reciprocal space interactions should be included
            
         --------------------------------------------------------------------------------------- */
          
      void calculateEwaldIxn(int numberOfAtoms, std::vector<OpenMM::Vec3>& atomCoordinates,
                             std::vector<std::vector<double> >& atomParameters, std::vector<std::set<int> >& exclusions,
                             std::vector<OpenMM::Vec3>& forces, double* totalEnergy, bool includeDirect, bool includeReciprocal) const;
};

} // namespace OpenMM

#endif // __ReferenceLJCoulombIxn_H__
