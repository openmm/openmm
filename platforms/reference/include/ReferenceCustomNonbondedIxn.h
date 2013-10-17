
/* Portions copyright (c) 2009-2013 Stanford University and Simbios.
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

#ifndef __ReferenceCustomNonbondedxIxn_H__
#define __ReferenceCustomNonbondedxIxn_H__

#include "ReferencePairIxn.h"
#include "ReferenceNeighborList.h"
#include "lepton/CompiledExpression.h"
#include <map>
#include <set>
#include <utility>
#include <vector>

// ---------------------------------------------------------------------------------------

class ReferenceCustomNonbondedIxn {

   private:

      bool cutoff;
      bool useSwitch;
      bool periodic;
      const OpenMM::NeighborList* neighborList;
      RealOpenMM periodicBoxSize[3];
      RealOpenMM cutoffDistance, switchingDistance;
      Lepton::CompiledExpression energyExpression;
      Lepton::CompiledExpression forceExpression;
      std::vector<std::string> paramNames;
      std::vector<double*> energyParticleParams;
      std::vector<double*> forceParticleParams;
      double* energyR;
      double* forceR;
      std::vector<std::pair<std::set<int>, std::set<int> > > interactionGroups;

      /**---------------------------------------------------------------------------------------

         Calculate custom pair ixn between two atoms

         @param atom1            the index of the first atom
         @param atom2            the index of the second atom
         @param atomCoordinates  atom coordinates
         @param atomParameters   atomParameters[atomIndex][parameterIndex]
         @param forces           force array (forces added)
         @param energyByAtom     atom energy
         @param totalEnergy      total energy

         --------------------------------------------------------------------------------------- */

      void calculateOneIxn( int atom1, int atom2, std::vector<OpenMM::RealVec>& atomCoordinates, std::vector<OpenMM::RealVec>& forces,
                            RealOpenMM* energyByAtom, RealOpenMM* totalEnergy );


   public:

      /**---------------------------------------------------------------------------------------

         Constructor

         --------------------------------------------------------------------------------------- */

       ReferenceCustomNonbondedIxn(const Lepton::CompiledExpression& energyExpression, const Lepton::CompiledExpression& forceExpression,
                                   const std::vector<std::string>& parameterNames);

      /**---------------------------------------------------------------------------------------

         Destructor

         --------------------------------------------------------------------------------------- */

       ~ReferenceCustomNonbondedIxn( );

      /**---------------------------------------------------------------------------------------

         Set the force to use a cutoff.

         @param distance            the cutoff distance
         @param neighbors           the neighbor list to use

         --------------------------------------------------------------------------------------- */

      void setUseCutoff( RealOpenMM distance, const OpenMM::NeighborList& neighbors );

      /**---------------------------------------------------------------------------------------

         Restrict the force to a list of interaction groups.

         @param distance            the cutoff distance
         @param neighbors           the neighbor list to use

         --------------------------------------------------------------------------------------- */

      void setInteractionGroups(const std::vector<std::pair<std::set<int>, std::set<int> > >& groups);

      /**---------------------------------------------------------------------------------------
      
         Set the force to use a switching function.
      
         @param distance            the switching distance
      
         --------------------------------------------------------------------------------------- */
      
      void setUseSwitchingFunction( RealOpenMM distance );

      /**---------------------------------------------------------------------------------------

         Set the force to use periodic boundary conditions.  This requires that a cutoff has
         already been set, and the smallest side of the periodic box is at least twice the cutoff
         distance.

         @param boxSize             the X, Y, and Z widths of the periodic box

         --------------------------------------------------------------------------------------- */

      void setPeriodic( OpenMM::RealVec& boxSize );

      /**---------------------------------------------------------------------------------------

         Calculate custom pair ixn

         @param numberOfAtoms    number of atoms
         @param atomCoordinates  atom coordinates
         @param atomParameters   atom parameters (charges, c6, c12, ...)     atomParameters[atomIndex][paramterIndex]
         @param exclusions       atom exclusion indices
                                 exclusions[atomIndex] contains the list of exclusions for that atom
         @param fixedParameters  non atom parameters (not currently used)
         @param globalParameters the values of global parameters
         @param forces           force array (forces added)
         @param energyByAtom     atom energy
         @param totalEnergy      total energy

         --------------------------------------------------------------------------------------- */

      void calculatePairIxn( int numberOfAtoms, std::vector<OpenMM::RealVec>& atomCoordinates,
                            RealOpenMM** atomParameters, std::vector<std::set<int> >& exclusions,
                            RealOpenMM* fixedParameters, const std::map<std::string, double>& globalParameters,
                            std::vector<OpenMM::RealVec>& forces, RealOpenMM* energyByAtom, RealOpenMM* totalEnergy );

// ---------------------------------------------------------------------------------------

};

#endif // __ReferenceCustomNonbondedxIxn_H__
