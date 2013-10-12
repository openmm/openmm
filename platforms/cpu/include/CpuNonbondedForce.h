
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

#ifndef OPENMM_CPU_NONBONDED_FORCE_H__
#define OPENMM_CPU_NONBONDED_FORCE_H__

#include "ReferencePairIxn.h"
#include <pthread.h>
#include <set>
#include <utility>
#include <vector>
#include <smmintrin.h>
// ---------------------------------------------------------------------------------------

class CpuNonbondedForce {
    public:
        class ThreadData;

      /**---------------------------------------------------------------------------------------
      
         Constructor
      
         --------------------------------------------------------------------------------------- */

       CpuNonbondedForce();

      /**---------------------------------------------------------------------------------------
      
         Destructor
      
         --------------------------------------------------------------------------------------- */

       ~CpuNonbondedForce();

      /**---------------------------------------------------------------------------------------
      
         Set the force to use a cutoff.
      
         @param distance            the cutoff distance
         @param neighbors           the neighbor list to use
         @param solventDielectric   the dielectric constant of the bulk solvent
      
         --------------------------------------------------------------------------------------- */
      
      void setUseCutoff(float distance, const std::vector<std::pair<int, int> >& neighbors, float solventDielectric);

      /**---------------------------------------------------------------------------------------
      
         Set the force to use a switching function on the Lennard-Jones interaction.
      
         @param distance            the switching distance
      
         --------------------------------------------------------------------------------------- */
      
      void setUseSwitchingFunction(float distance);
      
      /**---------------------------------------------------------------------------------------
      
         Set the force to use periodic boundary conditions.  This requires that a cutoff has
         already been set, and the smallest side of the periodic box is at least twice the cutoff
         distance.
      
         @param boxSize             the X, Y, and Z widths of the periodic box
      
         --------------------------------------------------------------------------------------- */
      
      void setPeriodic(float* periodicBoxSize);
       
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
      
         Calculate Ewald ixn
      
         @param numberOfAtoms    number of atoms
         @param posq             atom coordinates and charges
         @param atomCoordinates  atom coordinates (in format needed by PME)
         @param atomParameters   atom parameters (sigma/2, 2*sqrt(epsilon))
         @param exclusions       atom exclusion indices
                                 exclusions[atomIndex] contains the list of exclusions for that atom
         @param forces           force array (forces added)
         @param totalEnergy      total energy
            
         --------------------------------------------------------------------------------------- */
          
      void calculateReciprocalIxn(int numberOfAtoms, float* posq, std::vector<OpenMM::RealVec>& atomCoordinates,
                            const std::vector<std::pair<float, float> >& atomParameters, const std::vector<std::set<int> >& exclusions,
                            std::vector<OpenMM::RealVec>& forces, float* totalEnergy) const;
      
      /**---------------------------------------------------------------------------------------
      
         Calculate LJ Coulomb pair ixn
      
         @param numberOfAtoms    number of atoms
         @param posq             atom coordinates and charges
         @param atomParameters   atom parameters (sigma/2, 2*sqrt(epsilon))
         @param exclusions       atom exclusion indices
                                 exclusions[atomIndex] contains the list of exclusions for that atom
         @param forces           force array (forces added)
         @param totalEnergy      total energy
      
         --------------------------------------------------------------------------------------- */
          
      void calculateDirectIxn(int numberOfAtoms, float* posq, const std::vector<std::pair<float, float> >& atomParameters,
            const std::vector<std::set<int> >& exclusions, float* forces, float* totalEnergy);

    /**
     * This routine contains the code executed by each thread.
     */
    void runThread(int index, std::vector<float>& threadForce, double& threadEnergy);

private:
        bool cutoff;
        bool useSwitch;
        bool periodic;
        bool ewald;
        bool pme;
        const std::vector<std::pair<int, int> >* neighborList;
        float periodicBoxSize[3];
        float cutoffDistance, switchingDistance;
        float krf, crf;
        float alphaEwald;
        int numRx, numRy, numRz;
        int meshDim[3];
        __m128 boxSize, invBoxSize, half;
        bool isDeleted;
        int numThreads, waitCount;
        std::vector<pthread_t> thread;
        std::vector<ThreadData*> threadData;
        pthread_cond_t startCondition, endCondition;
        pthread_mutex_t lock;
        // The following variables are used to make information accessible to the individual threads.
        float* posq;
        std::vector<std::pair<float, float> > atomParameters;        
        std::vector<std::set<int> > exclusions;
        bool includeEnergy;

        static float TWO_OVER_SQRT_PI;
            
      /**---------------------------------------------------------------------------------------
      
         Calculate LJ Coulomb pair ixn between two atoms
      
         @param atom1            the index of the first atom
         @param atom2            the index of the second atom
         @param forces           force array (forces added)
         @param totalEnergy      total energy
            
         --------------------------------------------------------------------------------------- */
          
      void calculateOneIxn(int atom1, int atom2, float* forces, double* totalEnergy);
            
      /**---------------------------------------------------------------------------------------
      
         Calculate LJ Coulomb pair ixn between two atoms
      
         @param atom1            the index of the first atom
         @param atom2            the index of the second atom
         @param forces           force array (forces added)
         @param totalEnergy      total energy
            
         --------------------------------------------------------------------------------------- */
          
      void calculateOneEwaldIxn(int atom1, int atom2, float* forces, double* totalEnergy);

      
      void getDeltaR(const __m128& posI, const __m128& posJ, __m128& deltaR, float& r2, bool periodic) const;

      static float erfcApprox(float x);
};

// ---------------------------------------------------------------------------------------

#endif // OPENMM_CPU_NONBONDED_FORCE_H__
