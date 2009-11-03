
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

#ifndef __ReferenceFreeEnergyLJCoulomb14Softcore_H__
#define __ReferenceFreeEnergyLJCoulomb14Softcore_H__

#include "SimTKReference/ReferenceBondIxn.h"

// ---------------------------------------------------------------------------------------

class ReferenceFreeEnergyLJCoulomb14Softcore : public ReferenceBondIxn {

   private:

        bool cutoff;
        RealOpenMM cutoffDistance;
        RealOpenMM krf, crf;

   public:

      /**---------------------------------------------------------------------------------------
      
         Constructor
      
         --------------------------------------------------------------------------------------- */

       ReferenceFreeEnergyLJCoulomb14Softcore( );

      /**---------------------------------------------------------------------------------------
      
         Destructor
      
         --------------------------------------------------------------------------------------- */

       ~ReferenceFreeEnergyLJCoulomb14Softcore( );

      /**---------------------------------------------------------------------------------------
      
         Set the force to use a cutoff.
      
         @param distance            the cutoff distance
         @param solventDielectric   the dielectric constant of the bulk solvent
      
         @return ReferenceForce::DefaultReturn
      
         --------------------------------------------------------------------------------------- */
      
      int setUseCutoff( RealOpenMM distance, RealOpenMM solventDielectric );
       
      /**---------------------------------------------------------------------------------------
      
         Calculate parameters for LJ 1-4 ixn
      
         @param c6               c6
         @param c12              c12
         @param q1               q1 charge atom 1
         @param q2               q2 charge atom 2
         @param epsfac           epsfac ????????????
         @param parameters       output parameters:
                                    parameter[0]= c6*c6/c12
                                    parameter[1]= (c12/c6)**1/6
                                    parameter[2]= epsfactor*q1*q2
      
         @return ReferenceForce::DefaultReturn
      
         --------------------------------------------------------------------------------------- */
      
      int getDerivedParameters( RealOpenMM c6, RealOpenMM c12, RealOpenMM q1, 
                                               RealOpenMM q2, RealOpenMM epsfac,
                                               RealOpenMM* parameters ) const;
      
      /**---------------------------------------------------------------------------------------
      
         Calculate Ryckaert-Bellemans bond ixn
      
         @param atomIndices      atom indices of 4 atoms in bond
         @param atomCoordinates  atom coordinates
         @param parameters       six RB parameters
         @param forces           force array (forces added to current values)
         @param energiesByBond   energies by bond: energiesByBond[bondIndex]
         @param energiesByAtom   energies by atom: energiesByAtom[atomIndex]
      
         @return ReferenceForce::DefaultReturn
            
         --------------------------------------------------------------------------------------- */
      
      int calculateBondIxn( int* atomIndices, RealOpenMM** atomCoordinates,
                            RealOpenMM* parameters, RealOpenMM** forces,
                            RealOpenMM* energiesByBond, RealOpenMM* energiesByAtom ) const;
      
        /**---------------------------------------------------------------------------------------
      
           Calculate LJ pair ixn between two atoms
      
           @param inverseR         1/r
           @param sig              sigma
           @param eps              epsilon
           @param dEdR             output force factor
           @param energy           LJ energy
      
           @return ReferenceForce::DefaultReturn
      
           --------------------------------------------------------------------------------------- */
      
      int calculateOneLJ14Ixn( RealOpenMM inverseR, RealOpenMM sig, RealOpenMM eps,
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
      
      int calculateOneSoftCoreLJ14Ixn( RealOpenMM r, RealOpenMM sig, RealOpenMM eps,
                                       RealOpenMM lambda, RealOpenMM* dEdR, RealOpenMM* energy ) const;
      
};

// ---------------------------------------------------------------------------------------

#endif // __ReferenceFreeEnergyLJCoulomb14Softcore_H__
