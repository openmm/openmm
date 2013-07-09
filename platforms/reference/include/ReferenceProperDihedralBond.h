
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

#ifndef __ReferenceProperDihedralBond_H__
#define __ReferenceProperDihedralBond_H__

#include "ReferenceBondIxn.h"

// ---------------------------------------------------------------------------------------

class ReferenceProperDihedralBond : public ReferenceBondIxn {

   private:

   public:

      /**---------------------------------------------------------------------------------------
      
         Constructor
      
         --------------------------------------------------------------------------------------- */

       ReferenceProperDihedralBond( );

      /**---------------------------------------------------------------------------------------
      
         Destructor
      
         --------------------------------------------------------------------------------------- */

       ~ReferenceProperDihedralBond( );

      /**---------------------------------------------------------------------------------------
      
        Calculate proper dihedral bond ixn
      
         @param atomIndices      atom indices of 4 atoms in bond
         @param atomCoordinates  atom coordinates
         @param parameters       3 parameters: parameters[0] = k
                                               parameters[1] = ideal bond angle in radians
                                               parameters[2] = multiplicity
         @param forces           force array (forces added to current values)
         @param totalEnergy      if not null, the energy will be added to this
      
         --------------------------------------------------------------------------------------- */
      
      void calculateBondIxn( int* atomIndices, std::vector<OpenMM::RealVec>& atomCoordinates,
                            RealOpenMM* parameters, std::vector<OpenMM::RealVec>& forces,
                            RealOpenMM* totalEnergy ) const;

};

// ---------------------------------------------------------------------------------------

#endif // __ReferenceProperDihedralBond_H__
