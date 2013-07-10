
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

#ifndef __AmoebaReferenceWcaDispersionForce_H__
#define __AmoebaReferenceWcaDispersionForce_H__

#include "RealVec.h"
#include "openmm/Vec3.h"
#include <string>
#include <vector>

using namespace OpenMM;

// ---------------------------------------------------------------------------------------

class AmoebaReferenceWcaDispersionForce {

public:
 
    /**---------------------------------------------------------------------------------------
       
       Constructor

        @param epso      water oxygen eps for implicit dispersion term
        @param epsh      water hydrogen eps for implicit dispersion term
        @param rmino     water oxygen Rmin for implicit dispersion term
        @param rminh     water hydrogen Rmin for implicit dispersion term
        @param awater    water number density at standard temp & pressure
        @param shctd     overlap scale factor for HCT descreening method
        @param dispoff   dispersion radius offsets
        @param slevy     enthalpy-to-free energy scale factor for dispersion

       --------------------------------------------------------------------------------------- */
 
    AmoebaReferenceWcaDispersionForce( RealOpenMM epso, RealOpenMM epsh, RealOpenMM rmino, RealOpenMM rminh, 
                                       RealOpenMM awater, RealOpenMM shctd, RealOpenMM dispoff, RealOpenMM slevy );
 
    /**---------------------------------------------------------------------------------------
       
       Destructor
       
       --------------------------------------------------------------------------------------- */
 
    ~AmoebaReferenceWcaDispersionForce( ){};
 
    /**---------------------------------------------------------------------------------------
    
       Calculate WcaDispersion ixns
    
       @param numParticles                 number of particles
       @param particlePositions            Cartesian coordinates of particles
       @param radii                        particle radii
       @param epsilons                     particle epsilons
       @param totalMaximumDispersionEnergy total of maximum dispersion energy
       @param forces                       add forces to this vector
    
       @return energy
    
       --------------------------------------------------------------------------------------- */
    
    RealOpenMM calculateForceAndEnergy( int numParticles, const std::vector<OpenMM::RealVec>& particlePositions, 
                                        const std::vector<RealOpenMM>& radii, 
                                        const std::vector<RealOpenMM>& epsilons,
                                        RealOpenMM totalMaximumDispersionEnergy, std::vector<OpenMM::RealVec>& forces ) const;
private:

    RealOpenMM _epso; 
    RealOpenMM _epsh; 
    RealOpenMM _rmino; 
    RealOpenMM _rminh; 
    RealOpenMM _awater; 
    RealOpenMM _shctd; 
    RealOpenMM _dispoff;
    RealOpenMM _slevy;

    enum { EMIXO, RMIXO, RMIXO7, AO, EMIXH, RMIXH, RMIXH7, AH, LastIntermediateValueIndex }; 

    /**---------------------------------------------------------------------------------------
    
       Calculate pair ixn
    
       @param  radiusI              radius of particle I
       @param  radiusJ              radius of particle J
       @param  particleIPosition    particle I position 
       @param  particleJPosition    particle J position 
       @param  intermediateValues   intermediate values dependent on particle I
       @param  force                output force
    
       @return energy for ixn

       --------------------------------------------------------------------------------------- */
    
    RealOpenMM calculatePairIxn( RealOpenMM radiusI, RealOpenMM radiusJ,
                                 const OpenMM::RealVec& particleIPosition, const OpenMM::RealVec& particleJPosition,
                                 const RealOpenMM* const intermediateValues,
                                 Vec3& force ) const;

};

// ---------------------------------------------------------------------------------------

#endif // _AmoebaReferenceWcaDispersionForce___
