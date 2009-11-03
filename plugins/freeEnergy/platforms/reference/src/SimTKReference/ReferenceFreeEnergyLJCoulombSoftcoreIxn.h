#ifndef __ReferenceFreeEnergyLJCoulombSoftcoreIxn_H__
#define __ReferenceFreeEnergyLJCoulombSoftcoreIxn_H__

/* Portions copyright (c) 2006 Stanford University and Simbios.
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

#include "SimTKReference/ReferencePairIxn.h"
#include "SimTKReference/ReferenceNeighborList.h"

// ---------------------------------------------------------------------------------------

class ReferenceFreeEnergyLJCoulombSoftcoreIxn : public ReferencePairIxn {

   private:
       
      bool cutoff;
      bool periodic;
      bool ewald;
      bool pme;
      const OpenMM::NeighborList* neighborList;
      RealOpenMM periodicBoxSize[3];
      RealOpenMM cutoffDistance;
      RealOpenMM krf, crf;
      RealOpenMM softCoreLJLambda;
      int numRx, numRy, numRz;
      RealOpenMM alphaEwald;

      // parameter indices

      static const int SigIndex                   = 0;
      static const int EpsIndex                   = 1;
      static const int   QIndex                   = 2;
      static const int SoftCoreLJLambdaIndex      = 3;
            
      /**---------------------------------------------------------------------------------------
      
         Calculate LJ Coulomb pair ixn between two atoms
      
         @param atom1            the index of the first atom
         @param atom2            the index of the second atom
         @param atomCoordinates  atom coordinates
         @param atomParameters   atom parameters (charges, c6, c12, ...)     atomParameters[atomIndex][paramterIndex]
         @param forces           force array (forces added)
         @param energyByAtom     atom energy
         @param totalEnergy      total energy
      
         @return ReferenceForce::DefaultReturn
            
         --------------------------------------------------------------------------------------- */
          
      int calculateOneIxn( int atom1, int atom2, RealOpenMM** atomCoordinates,
                            RealOpenMM** atomParameters, RealOpenMM** forces,
                            RealOpenMM* energyByAtom, RealOpenMM* totalEnergy ) const;


   public:

      /**---------------------------------------------------------------------------------------
      
         Constructor
      
         --------------------------------------------------------------------------------------- */

       ReferenceFreeEnergyLJCoulombSoftcoreIxn( );

      /**---------------------------------------------------------------------------------------
      
         Destructor
      
         --------------------------------------------------------------------------------------- */

       ~ReferenceFreeEnergyLJCoulombSoftcoreIxn( );

      /**---------------------------------------------------------------------------------------
      
         Set the force to use a cutoff.
      
         @param distance            the cutoff distance
         @param neighbors           the neighbor list to use
         @param solventDielectric   the dielectric constant of the bulk solvent
      
         @return ReferenceForce::DefaultReturn
      
         --------------------------------------------------------------------------------------- */
      
      int setUseCutoff( RealOpenMM distance, const OpenMM::NeighborList& neighbors, RealOpenMM solventDielectric );
      
      /**---------------------------------------------------------------------------------------
      
         Set the force to use periodic boundary conditions.  This requires that a cutoff has
         already been set, and the smallest side of the periodic box is at least twice the cutoff
         distance.
      
         @param boxSize             the X, Y, and Z widths of the periodic box
      
         @return ReferenceForce::DefaultReturn
      
         --------------------------------------------------------------------------------------- */
      
      int setPeriodic( RealOpenMM* boxSize );      
       
      /**---------------------------------------------------------------------------------------
      
         Set the force to use Ewald summation.
      
         @param alpha  the Ewald separation parameter
         @param kmaxx  the largest wave vector in the x direction
         @param kmaxy  the largest wave vector in the y direction
         @param kmaxz  the largest wave vector in the z direction
      
         --------------------------------------------------------------------------------------- */
      
      void setUseEwald(RealOpenMM alpha, int kmaxx, int kmaxy, int kmaxz);

     
      /**---------------------------------------------------------------------------------------
      
         Set the force to use Particle-Mesh Ewald (PME) summation.
      
         @param alpha  the Ewald separation parameter
      
         --------------------------------------------------------------------------------------- */
      
      void setUsePME(RealOpenMM alpha);

      /**---------------------------------------------------------------------------------------
      
         Set the soft core LJ lambda
      
         @param lambda the soft core LJ lambda
      
         --------------------------------------------------------------------------------------- */
      
      void setSoftCoreLJLambda(RealOpenMM lambda);

      /**---------------------------------------------------------------------------------------
      
         Calculate parameters for LJ 1-4 ixn
      
         @param c6               c6
         @param c12              c12
         @param q1               q1 charge atom
         @param epsfacSqrt       epsfacSqrt (what is this?)
         @param parameters       output parameters:
                                    parameter[SigIndex]  = sqrt(c6*c6/c12)
                                    parameter[EpsIndex]  = 0.5*( (c12/c6)**1/6 )
                                    parameter[QIndex]    = epsfactorSqrt*q1
      
         @return ReferenceForce::DefaultReturn
      
         --------------------------------------------------------------------------------------- */
      
      int getDerivedParameters( RealOpenMM c6, RealOpenMM c12, RealOpenMM q1, 
                                RealOpenMM epsfacSqrt,
                                RealOpenMM* parameters ) const;
      
      /**---------------------------------------------------------------------------------------
      
         Calculate LJ Coulomb pair ixn
      
         @param numberOfAtoms    number of atoms
         @param atomCoordinates  atom coordinates
         @param atomParameters   atom parameters (charges, c6, c12, ...)     atomParameters[atomIndex][paramterIndex]
         @param exclusions       atom exclusion indices                      exclusions[atomIndex][atomToExcludeIndex]
                                 exclusions[atomIndex][0] = number of exclusions
                                 exclusions[atomIndex][1-no.] = atom indices of atoms to excluded from
                                 interacting w/ atom atomIndex
         @param fixedParameters  non atom parameters (not currently used)
         @param forces           force array (forces added)
         @param energyByAtom     atom energy
         @param totalEnergy      total energy
      
         @return ReferenceForce::DefaultReturn
            
         --------------------------------------------------------------------------------------- */
          
      int calculatePairIxn( int numberOfAtoms, RealOpenMM** atomCoordinates,
                            RealOpenMM** atomParameters, int** exclusions,
                            RealOpenMM* fixedParameters, RealOpenMM** forces,
                            RealOpenMM* energyByAtom, RealOpenMM* totalEnergy ) const;

private:
      /**---------------------------------------------------------------------------------------
      
         Calculate Ewald ixn
      
         @param numberOfAtoms    number of atoms
         @param atomCoordinates  atom coordinates
         @param atomParameters   atom parameters (charges, c6, c12, ...)     atomParameters[atomIndex][paramterIndex]
         @param exclusions       atom exclusion indices                      exclusions[atomIndex][atomToExcludeIndex]
                                 exclusions[atomIndex][0] = number of exclusions
                                 exclusions[atomIndex][1-no.] = atom indices of atoms to excluded from
                                 interacting w/ atom atomIndex
         @param fixedParameters  non atom parameters (not currently used)
         @param forces           force array (forces added)
         @param energyByAtom     atom energy
         @param totalEnergy      total energy

         @return ReferenceForce::DefaultReturn
            
         --------------------------------------------------------------------------------------- */
          
      int calculateEwaldIxn( int numberOfAtoms, RealOpenMM** atomCoordinates,
                            RealOpenMM** atomParameters, int** exclusions,
                            RealOpenMM* fixedParameters, RealOpenMM** forces,
                            RealOpenMM* energyByAtom, RealOpenMM* totalEnergy ) const;
      
      /**---------------------------------------------------------------------------------------
      
         Calculate PME ixn
      
         @param numberOfAtoms    number of atoms
         @param atomCoordinates  atom coordinates
         @param atomParameters   atom parameters (charges, c6, c12, ...)     atomParameters[atomIndex][paramterIndex]
         @param exclusions       atom exclusion indices                      exclusions[atomIndex][atomToExcludeIndex]
                                 exclusions[atomIndex][0] = number of exclusions
                                 exclusions[atomIndex][1-no.] = atom indices of atoms to excluded from
                                 interacting w/ atom atomIndex
         @param fixedParameters  non atom parameters (not currently used)
         @param forces           force array (forces added)
         @param energyByAtom     atom energy
         @param totalEnergy      total energy

         @return ReferenceForce::DefaultReturn
            
         --------------------------------------------------------------------------------------- */
          
      int calculatePMEIxn( int numberOfAtoms, RealOpenMM** atomCoordinates,
                            RealOpenMM** atomParameters, int** exclusions,
                            RealOpenMM* fixedParameters, RealOpenMM** forces,
                            RealOpenMM* energyByAtom, RealOpenMM* totalEnergy ) const;

        /**---------------------------------------------------------------------------------------
      
           Calculate LJ pair ixn between two atoms
      
           @param inverseR         1/r
           @param sig              sigma
           @param eps              epsilon
           @param dEdR             output force factor
           @param energy           LJ energy
      
           @return ReferenceForce::DefaultReturn
      
           --------------------------------------------------------------------------------------- */
      
      int calculateOneLJIxn( RealOpenMM inverseR, RealOpenMM sig, RealOpenMM eps,
                             RealOpenMM* dEdR, RealOpenMM* energy ) const;
      
        /**---------------------------------------------------------------------------------------
      
           Calculate softcore LJ pair ixn between two atoms
      
           @param r                r
           @param sig              sigma
           @param eps              epsilon
           @param lambda           lambda
           @param dEdR             output force factor
           @param energy           LJ energy
      
           @return ReferenceForce::DefaultReturn
      
           --------------------------------------------------------------------------------------- */
      
      int calculateOneSoftCoreLJIxn( RealOpenMM r, RealOpenMM sig, RealOpenMM eps,
                                     RealOpenMM lambda, RealOpenMM* dEdR, RealOpenMM* energy ) const;
      
      
};

// ---------------------------------------------------------------------------------------

#endif // __ReferenceFreeEnergyLJCoulombSoftcoreIxn_H__
