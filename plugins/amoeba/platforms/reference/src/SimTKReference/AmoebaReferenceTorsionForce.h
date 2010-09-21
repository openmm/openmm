
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

#ifndef __AmoebaReferenceTorsionForce_H__
#define __AmoebaReferenceTorsionForce_H__

#include "SimTKUtilities/SimTKOpenMMRealType.h"
#include <vector>

// ---------------------------------------------------------------------------------------

class AmoebaReferenceTorsionForce {

public:
 
    /**---------------------------------------------------------------------------------------
       
       Constructor
       
       --------------------------------------------------------------------------------------- */
 
    AmoebaReferenceTorsionForce( ){};
 
    /**---------------------------------------------------------------------------------------
       
       Destructor
       
          --------------------------------------------------------------------------------------- */
 
    ~AmoebaReferenceTorsionForce( ){};
 
     /**---------------------------------------------------------------------------------------
     
        Calculate Amoeba torsion ixns (force and energy)
     
        @param numTorsions             number of torsions
        @param posData                 particle positions
        @param particle1               particle 1 indices
        @param particle2               particle 2 indices
        @param particle3               particle 3 indices
        @param particle4               particle 4 indices
        @param torsionParameters1      first  index torsion parameters (amplitude, phase, fold)
        @param torsionParameters2      second index torsion parameters (amplitude, phase, fold)
        @param torsionParameters3      third  index torsion parameters (amplitude, phase, fold)
        @param forces                  output force vector
     
        @return total energy

     
        --------------------------------------------------------------------------------------- */

    RealOpenMM calculateForceAndEnergy( int numTorsions, RealOpenMM** posData,
                                        const std::vector<int>&  particle1,
                                        const std::vector<int>&  particle2,
                                        const std::vector<int>&  particle3,
                                        const std::vector<int>&  particle4,
                                        const std::vector< std::vector<RealOpenMM> >& torsionParameters1,
                                        const std::vector< std::vector<RealOpenMM> >& torsionParameters2,
                                        const std::vector< std::vector<RealOpenMM> >& torsionParameters3,
                                        RealOpenMM** forceData ) const;

private:

    /**---------------------------------------------------------------------------------------
    
       Calculate Amoeba harmonic angle ixn (force and energy)
    
       @param positionAtomA           Cartesian coordinates of atom A
       @param positionAtomB           Cartesian coordinates of atom B
       @param positionAtomC           Cartesian coordinates of atom C
       @param positionAtomD           Cartesian coordinates of atom D
       @param torsion1                vector of torsion params for first  index (amplitude, phase, fold)
       @param torsion2                vector of torsion params for second index (amplitude, phase, fold)
       @param torsion3                vector of torsion params for third  index (amplitude, phase, fold)
       @param forces                  force vector
    
       @return energy
    
       --------------------------------------------------------------------------------------- */
    
    RealOpenMM calculateTorsionIxn( const RealOpenMM* positionAtomA, const RealOpenMM* positionAtomB,
                                               const RealOpenMM* positionAtomC, const RealOpenMM* positionAtomD,
                                               const std::vector<RealOpenMM>& torsionParameters1,
                                               const std::vector<RealOpenMM>& torsionParameters2,
                                               const std::vector<RealOpenMM>& torsionParameters3,
                                               RealOpenMM** forces ) const;
         
};


// ---------------------------------------------------------------------------------------

#endif // _AmoebaReferenceTorsionForce___
