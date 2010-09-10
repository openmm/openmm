
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

#ifndef __AmoebaReferenceStretchBendForce_H__
#define __AmoebaReferenceStretchBendForce_H__

#include "SimTKUtilities/SimTKOpenMMRealType.h"

// ---------------------------------------------------------------------------------------

class AmoebaReferenceStretchBendForce {

public:
 
    /**---------------------------------------------------------------------------------------
       
       Constructor
       
       --------------------------------------------------------------------------------------- */
 
    AmoebaReferenceStretchBendForce( );
 
    /**---------------------------------------------------------------------------------------
       
       Destructor
       
          --------------------------------------------------------------------------------------- */
 
    ~AmoebaReferenceStretchBendForce( );

    /**---------------------------------------------------------------------------------------
    
       Calculate Amoeba harmonic angle ixn (force and energy)
    
       @param positionAtomA           Cartesian coordinates of atom A
       @param positionAtomB           Cartesian coordinates of atom B
       @param positionAtomC           Cartesian coordinates of atom C
       @param lengthAB                ideal AB bondlength
       @param lengthCB                ideal CB bondlength
       @param idealAngle              ideal angle
       @param kParameter              k
       @param forces                  force vector
    
       @return energy
    
       --------------------------------------------------------------------------------------- */
    static RealOpenMM calculateForceAndEnergy( const RealOpenMM* positionAtomA, const RealOpenMM* positionAtomB,
                                               const RealOpenMM* positionAtomC,
                                               RealOpenMM lengthAB,      RealOpenMM lengthCB,
                                               RealOpenMM idealAngle,    RealOpenMM kParameter,
                                               RealOpenMM** forces );
 
};

// ---------------------------------------------------------------------------------------

#endif // _AmoebaReferenceStretchBendForce___
