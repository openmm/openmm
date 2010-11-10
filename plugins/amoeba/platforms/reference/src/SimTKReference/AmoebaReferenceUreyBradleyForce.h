
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

#ifndef __AmoebaReferenceUreyBradleyForce_H__
#define __AmoebaReferenceUreyBradleyForce_H__

#include "SimTKUtilities/SimTKOpenMMRealType.h"
#include <vector>

// ---------------------------------------------------------------------------------------

class AmoebaReferenceUreyBradleyForce {

public:
 
    /**---------------------------------------------------------------------------------------
       
        Constructor
       
        --------------------------------------------------------------------------------------- */
 
    AmoebaReferenceUreyBradleyForce( ){};
 
    /**---------------------------------------------------------------------------------------
       
        Destructor
       
        --------------------------------------------------------------------------------------- */
 
    ~AmoebaReferenceUreyBradleyForce( ){};
 
     /**---------------------------------------------------------------------------------------
     
        Calculate Amoeba Urey-Bradley ixns (force and energy)
     
        @param numIxns                 number of ixns
        @param posData                 particle positions
        @param particle1               particle 1 indices
        @param particle2               particle 2 indices
        @param length                  ideal length
        @param kForce                  force constant
        @param cubic                   cubic force parameter
        @param quartic                 quartic force parameter
        @param forceDate               output force vector
     
        @return total energy
     
        --------------------------------------------------------------------------------------- */
     
    RealOpenMM calculateForceAndEnergy( int numIxns, RealOpenMM** posData,
                                        const std::vector<int>& particle1,
                                        const std::vector<int>&  particle2,
                                        const std::vector<RealOpenMM>& length,
                                        const std::vector<RealOpenMM>& kForce,
                                        RealOpenMM cubic, RealOpenMM quartic,
                                        RealOpenMM** forceData ) const;

private:

     /**---------------------------------------------------------------------------------------
     
        Calculate Amoeba Urey-Bradley ixns (force and energy)
     
        @param positionAtomA           Cartesian coordinates of atom A
        @param positionAtomB           Cartesian coordinates of atom B
        @param length                  ideal length
        @param kForce                  force constant
        @param cubic                   cubic force parameter
        @param quartic                 quartic force parameter
        @param forces                  force vector
     
        @return energy
     
        --------------------------------------------------------------------------------------- */
     
    RealOpenMM calculateUreyBradleyIxn( const RealOpenMM* positionAtomA, const RealOpenMM* positionAtomB,
                                        RealOpenMM idealLength, RealOpenMM kForce,
                                        RealOpenMM cubic, RealOpenMM quartic,
                                        RealOpenMM** forces ) const;
     
};

// ---------------------------------------------------------------------------------------

#endif // _AmoebaReferenceUreyBradleyForce___
