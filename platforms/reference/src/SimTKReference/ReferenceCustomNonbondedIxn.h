
/* Portions copyright (c) 2009 Stanford University and Simbios.
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
#include "lepton/ExpressionProgram.h"
#include <map>
#include <vector>

// ---------------------------------------------------------------------------------------

class ReferenceCustomNonbondedIxn : public ReferencePairIxn {

   private:

      bool cutoff;
      bool periodic;
      const OpenMM::NeighborList* neighborList;
      RealOpenMM periodicBoxSize[3];
      RealOpenMM cutoffDistance;
      Lepton::ExpressionProgram energyExpression;
      Lepton::ExpressionProgram forceExpression;
      std::vector<Lepton::ExpressionProgram> combiningRules;
      std::vector<std::string> paramNames;

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

      void calculateOneIxn( int atom1, int atom2, RealOpenMM** atomCoordinates,
                            std::map<std::string, double>& variables, RealOpenMM** forces,
                            RealOpenMM* energyByAtom, RealOpenMM* totalEnergy ) const;


   public:

      /**---------------------------------------------------------------------------------------

         Constructor

         --------------------------------------------------------------------------------------- */

       ReferenceCustomNonbondedIxn(const Lepton::ExpressionProgram& energyExpression, const Lepton::ExpressionProgram& forceExpression,
                                   const std::vector<Lepton::ExpressionProgram>& combiningRules);

      /**---------------------------------------------------------------------------------------

         Destructor

         --------------------------------------------------------------------------------------- */

       ~ReferenceCustomNonbondedIxn( );

      /**---------------------------------------------------------------------------------------

         Set the force to use a cutoff.

         @param distance            the cutoff distance
         @param neighbors           the neighbor list to use

         @return ReferenceForce::DefaultReturn

         --------------------------------------------------------------------------------------- */

      int setUseCutoff( RealOpenMM distance, const OpenMM::NeighborList& neighbors );

      /**---------------------------------------------------------------------------------------

         Set the force to use periodic boundary conditions.  This requires that a cutoff has
         already been set, and the smallest side of the periodic box is at least twice the cutoff
         distance.

         @param boxSize             the X, Y, and Z widths of the periodic box

         @return ReferenceForce::DefaultReturn

         --------------------------------------------------------------------------------------- */

      int setPeriodic( RealOpenMM* boxSize );

      /**---------------------------------------------------------------------------------------

         Calculate custom pair ixn

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

      /**---------------------------------------------------------------------------------------

         Calculate custom pair ixn for exceptions

         @param numberOfExceptions   number of exceptions
         @param atomIndices          indices of atoms participating in exception ixn: atomIndices[exceptionIndex][indices]
         @param atomCoordinates      atom coordinates: atomCoordinates[atomIndex][3]
         @param parameters           parameters: parameters[excptionIndex][*]; contents of array
                                     depend on ixn
         @param forces               force array (forces added to current values): forces[atomIndex][3]
         @param energyByAtom         atom energy
         @param totalEnergy          total energy

         --------------------------------------------------------------------------------------- */

      void calculateExceptionIxn( int numberOfExceptions, int** atomIndices, RealOpenMM** atomCoordinates,
                                    RealOpenMM** parameters, RealOpenMM** forces,
                                    RealOpenMM* energyByAtom, RealOpenMM* totalEnergy ) const;

// ---------------------------------------------------------------------------------------

};

#endif // __ReferenceCustomNonbondedxIxn_H__
