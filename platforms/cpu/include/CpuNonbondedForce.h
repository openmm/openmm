
/* Portions copyright (c) 2006-2018 Stanford University and Simbios.
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

#ifndef OPENMM_CPU_NONBONDED_FORCE_H__
#define OPENMM_CPU_NONBONDED_FORCE_H__

#include "AlignedArray.h"
#include "CpuNeighborList.h"
#include "ReferencePairIxn.h"
#include "openmm/internal/ThreadPool.h"
#include "openmm/internal/vectorize.h"
#include <atomic>
#include <set>
#include <utility>
#include <vector>
// ---------------------------------------------------------------------------------------

namespace OpenMM {

class CpuNonbondedForce {
    public:

      /**---------------------------------------------------------------------------------------
      
         Constructor
      
         --------------------------------------------------------------------------------------- */

       CpuNonbondedForce();
       
        /**
         * Virtual destructor.
         */

        virtual ~CpuNonbondedForce();
        
      /**---------------------------------------------------------------------------------------
      
         Set the force to use a cutoff.
      
         @param distance            the cutoff distance
         @param neighbors           the neighbor list to use
         @param solventDielectric   the dielectric constant of the bulk solvent
      
         --------------------------------------------------------------------------------------- */
      
      void setUseCutoff(float distance, const CpuNeighborList& neighbors, float solventDielectric);

      /**---------------------------------------------------------------------------------------
      
         Set the force to use a switching function on the Lennard-Jones interaction.
      
         @param distance            the switching distance
      
         --------------------------------------------------------------------------------------- */
      
      void setUseSwitchingFunction(float distance);
      
      /**---------------------------------------------------------------------------------------
      
         Set the force to use periodic boundary conditions.  This requires that a cutoff has
         already been set, and the smallest side of the periodic box is at least twice the cutoff
         distance.
      
         @param periodicBoxVectors    the vectors defining the periodic box
      
         --------------------------------------------------------------------------------------- */
      
      void setPeriodic(Vec3* periodicBoxVectors);
       
      /**---------------------------------------------------------------------------------------
      
         Set the force to use Ewald summation.
      
         @param alpha  the Ewald separation parameter
         @param kmaxx  the largest wave vector in the x direction
         @param kmaxy  the largest wave vector in the y direction
         @param kmaxz  the largest wave vector in the z direction
      
         --------------------------------------------------------------------------------------- */
      
      void setUseEwald(float alpha, int kmaxx, int kmaxy, int kmaxz);

     
      /**---------------------------------------------------------------------------------------

         Set the force to use Particle-Mesh Ewald (PME) summation.

         @param alpha    the Ewald separation parameter
         @param gridSize the dimensions of the mesh

         --------------------------------------------------------------------------------------- */

      void setUsePME(float alpha, int meshSize[3]);

      /**---------------------------------------------------------------------------------------

         Set the force to use Particle-Mesh Ewald (PME) summation for dispersion.

         @param alpha    the Ewald separation parameter
         @param gridSize the dimensions of the mesh

         --------------------------------------------------------------------------------------- */

      void setUseLJPME(float alpha, int meshSize[3]);

      /**---------------------------------------------------------------------------------------
      
         Calculate Ewald ixn
      
         @param numberOfAtoms    number of atoms
         @param posq             atom coordinates and charges
         @param atomCoordinates  atom coordinates (in format needed by PME)
         @param atomParameters   atom parameters (sigma/2, 2*sqrt(epsilon))
         @param C6Paramrs        C6 parameters for multiplicative representation of dispersion
         @param exclusions       atom exclusion indices
                                 exclusions[atomIndex] contains the list of exclusions for that atom
         @param forces           force array (forces added)
         @param totalEnergy      total energy
            
         --------------------------------------------------------------------------------------- */

      void calculateReciprocalIxn(int numberOfAtoms, float* posq, const std::vector<Vec3>& atomCoordinates,
                                  const std::vector<std::pair<float, float> >& atomParameters, const std::vector<float> &C6params,
                                  const std::vector<std::set<int> >& exclusions, std::vector<Vec3>& forces, double* totalEnergy) const;
      
      /**---------------------------------------------------------------------------------------
      
         Calculate LJ Coulomb pair ixn
      
         @param numberOfAtoms    number of atoms
         @param posq             atom coordinates and charges
         @param atomCoordinates  atom coordinates (periodic boundary conditions not applied)
         @param atomParameters   atom parameters (sigma/2, 2*sqrt(epsilon))
         @param exclusions       atom exclusion indices
                                 exclusions[atomIndex] contains the list of exclusions for that atom
         @param forces           force array (forces added)
         @param totalEnergy      total energy
         @param threads          the thread pool to use
      
         --------------------------------------------------------------------------------------- */
          
      void calculateDirectIxn(int numberOfAtoms, float* posq, const std::vector<Vec3>& atomCoordinates, const std::vector<std::pair<float, float> >& atomParameters,
            const std::vector<float>& C6params, const std::vector<std::set<int> >& exclusions, std::vector<AlignedArray<float> >& threadForce, double* totalEnergy, ThreadPool& threads);

    /**
     * This routine contains the code executed by each thread.
     */
    void threadComputeDirect(ThreadPool& threads, int threadIndex);

protected:
        bool cutoff;
        bool useSwitch;
        bool periodic;
        bool triclinic;
        bool ewald;
        bool ljpme, pme;
        bool tableIsValid, expTableIsValid;
        const CpuNeighborList* neighborList;
        float recipBoxSize[3];
        Vec3 periodicBoxVectors[3];
        AlignedArray<fvec4> periodicBoxVec4;
        float cutoffDistance, switchingDistance;
        float krf, crf;
        float alphaEwald, alphaDispersionEwald;
        int numRx, numRy, numRz;
        int meshDim[3], dispersionMeshDim[3];
        std::vector<float> erfcTable, ewaldScaleTable;
        std::vector<float> exptermsTable, dExptermsTable;
        float ewaldDX, ewaldDXInv, erfcDXInv, exptermsDX, exptermsDXInv;
        std::vector<double> threadEnergy;
        // The following variables are used to make information accessible to the individual threads.
        int numberOfAtoms;
        float* posq;
        Vec3 const* atomCoordinates;
        std::pair<float, float> const* atomParameters;        
        float const *C6params;
        std::set<int> const* exclusions;
        std::vector<AlignedArray<float> >* threadForce;
        bool includeEnergy;
        float inverseRcut6;
        float inverseRcut6Expterm;
        std::atomic<int> atomicCounter;

        static const float TWO_OVER_SQRT_PI;
        static const int NUM_TABLE_POINTS;
            
      /**---------------------------------------------------------------------------------------
      
         Calculate LJ Coulomb pair ixn between two atoms
      
         @param atom1            the index of the first atom
         @param atom2            the index of the second atom
         @param forces           force array (forces added)
         @param totalEnergy      total energy
            
         --------------------------------------------------------------------------------------- */
          
      void calculateOneIxn(int atom1, int atom2, float* forces, double* totalEnergy, const fvec4& boxSize, const fvec4& invBoxSize);
            
      /**---------------------------------------------------------------------------------------
      
         Calculate all the interactions for one atom block.
      
         @param blockIndex       the index of the atom block
         @param forces           force array (forces added)
         @param totalEnergy      total energy
            
         --------------------------------------------------------------------------------------- */
          
      virtual void calculateBlockIxn(int blockIndex, float* forces, double* totalEnergy, const fvec4& boxSize, const fvec4& invBoxSize) = 0;
            
      /**---------------------------------------------------------------------------------------
      
         Calculate all the interactions for one atom block.
      
         @param blockIndex       the index of the atom block
         @param forces           force array (forces added)
         @param totalEnergy      total energy
            
         --------------------------------------------------------------------------------------- */
          
      virtual void calculateBlockEwaldIxn(int blockIndex, float* forces, double* totalEnergy, const fvec4& boxSize, const fvec4& invBoxSize) = 0;

      /**
       * Compute the displacement and squared distance between two points, optionally using
       * periodic boundary conditions.
       */
      void getDeltaR(const fvec4& posI, const fvec4& posJ, fvec4& deltaR, float& r2, bool periodic, const fvec4& boxSize, const fvec4& invBoxSize) const;

      /**
       * Create a lookup table for the scale factor used with Ewald and PME.
       */
      void tabulateEwaldScaleFactor();

      /**
       * Create a lookup table for the scale factor used with dispersion PME.
       */
      void tabulateExpTerms();

      /**
       * Compute a fast approximation to erfc(x).
       */
      float erfcApprox(float x);

      /**
       * Compute a fast approximation to (1.0 - EXP(-dar^2) * (1.0 + dar^2 + 0.5*dar^4))
       * where dar = (dispersionAlpha * R)
       * needed for LJPME energies.
       */
      float exptermsApprox(float R);

      /**
       * Compute a fast approximation to (1.0 - EXP(-dar^2) * (1.0 + dar^2 + 0.5*dar^4 + dar^6/6.0))
       * where dar = (dispersionAlpha * R)
       * needed for LJPME forces.
       */
      float dExptermsApprox(float R);
};

} // namespace OpenMM

// ---------------------------------------------------------------------------------------

#endif // OPENMM_CPU_NONBONDED_FORCE_H__
