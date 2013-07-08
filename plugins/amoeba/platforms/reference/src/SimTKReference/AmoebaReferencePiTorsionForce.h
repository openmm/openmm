
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

#ifndef __AmoebaReferencePiTorsionForce_H__
#define __AmoebaReferencePiTorsionForce_H__

#include "RealVec.h"
#include <vector>

// ---------------------------------------------------------------------------------------

class AmoebaReferencePiTorsionForce {

public:
 
    /**---------------------------------------------------------------------------------------
       
       Constructor
       
       --------------------------------------------------------------------------------------- */
 
    AmoebaReferencePiTorsionForce( ){};
 
    /**---------------------------------------------------------------------------------------
       
       Destructor
       
          --------------------------------------------------------------------------------------- */
 
    ~AmoebaReferencePiTorsionForce( ){};
 
     /**---------------------------------------------------------------------------------------
     
        Calculate Amoeba torsion ixns (force and energy)
     
        @param numTorsions             number of torsions
        @param posData                 particle positions
        @param particle1               particle 1 indices
        @param particle2               particle 2 indices
        @param particle3               particle 3 indices
        @param particle4               particle 4 indices
        @param particle5               particle 5 indices
        @param particle6               particle 6 indices
        @param torsionParameters1      first  index torsion parameters (amplitude, phase, fold)
        @param torsionParameters2      second index torsion parameters (amplitude, phase, fold)
        @param torsionParameters3      third  index torsion parameters (amplitude, phase, fold)
        @param forces                  output force vector
     
        @return total energy

     
        --------------------------------------------------------------------------------------- */

    RealOpenMM calculateForceAndEnergy( int numPiTorsions, std::vector<OpenMM::RealVec>& posData,
                                        const std::vector<int>&  particle1,
                                        const std::vector<int>&  particle2,
                                        const std::vector<int>&  particle3,
                                        const std::vector<int>&  particle4,
                                        const std::vector<int>&  particle5,
                                        const std::vector<int>&  particle6,
                                        const std::vector<RealOpenMM>& kTorsion,
                                        std::vector<OpenMM::RealVec>& forceData ) const;


private:

    /**---------------------------------------------------------------------------------------
    
       Calculate Amoeba pi-torsion ixn (force and energy)
    
       @param positionAtomA           Cartesian coordinates of atom A
       @param positionAtomB           Cartesian coordinates of atom B
       @param positionAtomC           Cartesian coordinates of atom C
       @param positionAtomD           Cartesian coordinates of atom D
       @param positionAtomE           Cartesian coordinates of atom E
       @param positionAtomF           Cartesian coordinates of atom F
       @param kTorsion                k-torsion parameter
       @param forces                  force vector
    
       @return energy
    
       --------------------------------------------------------------------------------------- */
    
    RealOpenMM calculatePiTorsionIxn( const OpenMM::RealVec& positionAtomA, const OpenMM::RealVec& positionAtomB,
                                      const OpenMM::RealVec& positionAtomC, const OpenMM::RealVec& positionAtomD,
                                      const OpenMM::RealVec& positionAtomE, const OpenMM::RealVec& positionAtomF,
                                      RealOpenMM kTorsion, OpenMM::RealVec* forces ) const;
         
};

// ---------------------------------------------------------------------------------------

#endif // _AmoebaReferencePiTorsionForce___
