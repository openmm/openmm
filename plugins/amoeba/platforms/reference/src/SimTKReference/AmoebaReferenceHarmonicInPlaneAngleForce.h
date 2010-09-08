
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

#ifndef __AmoebaReferenceHarmonicInPlaneAngleForce_H__
#define __AmoebaReferenceHarmonicInPlaneAngleForce_H__

#include "SimTKUtilities/SimTKOpenMMRealType.h"

// ---------------------------------------------------------------------------------------

class AmoebaReferenceHarmonicInPlaneAngleForce {

public:
 
    /**---------------------------------------------------------------------------------------
       
       Constructor
       
       --------------------------------------------------------------------------------------- */
 
    AmoebaReferenceHarmonicInPlaneAngleForce( );
 
    /**---------------------------------------------------------------------------------------
       
       Destructor
       
          --------------------------------------------------------------------------------------- */
 
    ~AmoebaReferenceHarmonicInPlaneAngleForce( );
 
    /**---------------------------------------------------------------------------------------
    
       Get dEdT and energy prefactor given cosine of angle :: the calculation for different
       the force types is identical 
    
       @param  cosine               cosine of angle
       @param  idealAngle           ideal angle
       @param  angleK               angle k (quadratic prefactor)
       @param  angleCubic           cubic prefactor
       @param  angleQuartic         quartiic prefactor
       @param  anglePentic          pentic prefactor
       @param  angleSextic          sextic prefactor
       @param  dEdR                 dEdR
    
       @return energy
    
       --------------------------------------------------------------------------------------- */
    
    static RealOpenMM getPrefactorsGivenAngleCosine( RealOpenMM cosine, RealOpenMM idealAngle, RealOpenMM angleK,
                                                     RealOpenMM angleCubic,     RealOpenMM angleQuartic,
                                                     RealOpenMM anglePentic,    RealOpenMM angleSextic,
                                                     RealOpenMM* dEdR );
    
    /**---------------------------------------------------------------------------------------
    
       Calculate Amoeba harmonic angle ixn (force and energy)
    
       @param positionAtomA           Cartesian coordinates of atom A
       @param positionAtomB           Cartesian coordinates of atom B
       @param positionAtomC           Cartesian coordinates of atom C
       @param positionAtomD           Cartesian coordinates of atom D
       @param angleLength             angle
       @param angleK                  quadratic angle force
       @param angleCubic              cubic angle force parameter
       @param angleQuartic            quartic angle force parameter
       @param anglePentic             pentic angle force parameter
       @param angleSextic             sextic angle force parameter
       @param forces                  force vector
    
       @return energy
    
       --------------------------------------------------------------------------------------- */
    
    static RealOpenMM calculateForceAndEnergy( RealOpenMM* positionAtomA, RealOpenMM* positionAtomB,
                                               RealOpenMM* positionAtomC, RealOpenMM* positionAtomD,
                                               RealOpenMM angle,          RealOpenMM angleK,
                                               RealOpenMM angleCubic,     RealOpenMM angleQuartic,
                                               RealOpenMM anglePentic,    RealOpenMM angleSextic,
                                               RealOpenMM** forces );
         
};

// ---------------------------------------------------------------------------------------

#endif // _AmoebaReferenceHarmonicInPlaneAngleForce___
