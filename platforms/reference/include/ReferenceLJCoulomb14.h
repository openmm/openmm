
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

#ifndef __ReferenceLJCoulomb14_H__
#define __ReferenceLJCoulomb14_H__

#include "ReferenceBondIxn.h"

// ---------------------------------------------------------------------------------------

class ReferenceLJCoulomb14 : public ReferenceBondIxn {

   public:

      /**---------------------------------------------------------------------------------------
      
         Constructor
      
         --------------------------------------------------------------------------------------- */

       ReferenceLJCoulomb14( );

      /**---------------------------------------------------------------------------------------
      
         Destructor
      
         --------------------------------------------------------------------------------------- */

       ~ReferenceLJCoulomb14( );

      /**---------------------------------------------------------------------------------------
      
         Calculate Ryckaert-Bellemans bond ixn
      
         @param atomIndices      atom indices of 4 atoms in bond
         @param atomCoordinates  atom coordinates
         @param parameters       six RB parameters
         @param forces           force array (forces added to current values)
         @param totalEnergy      if not null, the energy will be added to this
            
         --------------------------------------------------------------------------------------- */
      
      void calculateBondIxn( int* atomIndices, std::vector<OpenMM::RealVec>& atomCoordinates,
                            RealOpenMM* parameters, std::vector<OpenMM::RealVec>& forces,
                            RealOpenMM* totalEnergy ) const;

};

// ---------------------------------------------------------------------------------------

#endif // __ReferenceLJCoulomb14_H__
