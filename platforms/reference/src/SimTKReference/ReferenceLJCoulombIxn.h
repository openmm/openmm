
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

#ifndef __ReferenceLJCoulombIxn_H__
#define __ReferenceLJCoulombIxn_H__

#include "ReferencePairIxn.h"

// ---------------------------------------------------------------------------------------

class ReferenceLJCoulombIxn : public ReferencePairIxn {

   private:

      // parameter indices

      static const int SigIndex = 0;
      static const int EpsIndex = 1;
      static const int   QIndex = 2;

   public:

      /**---------------------------------------------------------------------------------------
      
         Constructor
      
         --------------------------------------------------------------------------------------- */

       ReferenceLJCoulombIxn( );

      /**---------------------------------------------------------------------------------------
      
         Destructor
      
         --------------------------------------------------------------------------------------- */

       ~ReferenceLJCoulombIxn( );

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
      
};

// ---------------------------------------------------------------------------------------

#endif // __ReferenceLJCoulombIxn_H__
