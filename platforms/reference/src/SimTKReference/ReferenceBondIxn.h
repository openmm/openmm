
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

#ifndef __ReferenceBondIxn_H__
#define __ReferenceBondIxn_H__

// #include "ReferenceIxn.h"

// ---------------------------------------------------------------------------------------

class OPENMM_EXPORT ReferenceBondIxn {

   private:

   public:

      /**---------------------------------------------------------------------------------------
      
         Constructor
      
         --------------------------------------------------------------------------------------- */

       ReferenceBondIxn( );

      /**---------------------------------------------------------------------------------------
      
         Destructor
      
         --------------------------------------------------------------------------------------- */

       ~ReferenceBondIxn( );

      /**---------------------------------------------------------------------------------------
      
         Calculate Bond Ixn -- virtual method
      
         @param atomIndices      bond indices
         @param atomCoordinates  atom coordinates
         @param parameters       parameters
         @param forces           force array (forces added)
         @param energyByBond     bond energy 
         @param energy           atom energy
      
         --------------------------------------------------------------------------------------- */
      
      virtual int calculateBondIxn( int* atomIndices, RealOpenMM** atomCoordinates,
                                    RealOpenMM* parameters, RealOpenMM** forces,
                                    RealOpenMM* energyByBond, RealOpenMM* energyByAtom ) const;

      /**---------------------------------------------------------------------------------------
      
         Update energy
      
         @param  energy               energy value to update 
         @param  energyByBond         ptr to energyByBond accumulator (may be null)
         @param  numberOfAtomIndices  number of atoms in bond
         @param  atomIndices          array of atom indices of size 'numberOfAtomIndices'
         @param  energyByAtom         array of energies by atom (may be null)
      
         @return ReferenceForce::DefaultReturn 
      
         --------------------------------------------------------------------------------------- */
      
      int updateEnergy( RealOpenMM energy, RealOpenMM* energyByBond,
                        int numberOfAtomIndices, int* atomIndices, RealOpenMM* energyByAtom ) const;
      
      
      /**---------------------------------------------------------------------------------------
      
         Get normed dot product between two vectors
      
         @param  vector1            first vector
         @param  vector2            second vector
         @param  hasREntry          if set, then vector1[ReferenceForce::RIndex] = norm of vector;
                                    defaults to 0
      
         @return dot product
      
         --------------------------------------------------------------------------------------- */
      
      RealOpenMM getNormedDotProduct( RealOpenMM* vector1, RealOpenMM* vector2, int hasREntry ) const;
      
      /**---------------------------------------------------------------------------------------
      
         Get angle between two vectors
      
         @param  vector1            first vector
         @param  vector2            second vector
         @param  outputDotProduct   cosine of angle between two vectors (optional)
         @param  hasREntry          if set, then vector1[ReferenceForce::R2Index] = square norm of vector;
                                    defaults to 0
      
         @return cosine of angles in radians
      
         --------------------------------------------------------------------------------------- */
      
      RealOpenMM getAngleBetweenTwoVectors( RealOpenMM* vector1, RealOpenMM* vector2, 
                                            RealOpenMM* outputDotProduct, int hasREntry ) const;
      
      /**---------------------------------------------------------------------------------------
      
         Get dihedral angle between two vectors
      
         @param  vector1            first vector
         @param  vector2            second vector
         @param  vector3            third vector
         @param  outputCrossProduct output cross product vectors
         @param  cosineOfAngle      cosine of angle (output)
         @param  signVector         vector to test sign (optional)
         @param  signOfAngle        sign of angle (output) (optional)
         @param  hasREntry          if set, then vector1[ReferenceForce::R2Index] = square norm of vector
                                    defaults to 0
      
         @return cosine of angles in radians
      
         --------------------------------------------------------------------------------------- */
      
      RealOpenMM getDihedralAngleBetweenThreeVectors( RealOpenMM* vector1, RealOpenMM* vector2, 
                                                      RealOpenMM* vector3, RealOpenMM** outputCrossProduct, 
                                                      RealOpenMM* cosineOfAngle, RealOpenMM* signVector, 
                                                      RealOpenMM* signOfAngle, int hasREntry ) const;
      
};

// ---------------------------------------------------------------------------------------

#endif // __ReferenceBondIxn_H__
