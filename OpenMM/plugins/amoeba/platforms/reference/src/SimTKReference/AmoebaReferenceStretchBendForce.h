
/* Portions copyright (c) 2006-2016 Stanford University and Simbios.
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

#include "RealVec.h"
#include <vector>

namespace OpenMM {

class AmoebaReferenceStretchBendForce {

public:
 
    /**---------------------------------------------------------------------------------------
       
       Constructor
       
       --------------------------------------------------------------------------------------- */
 
    AmoebaReferenceStretchBendForce() : usePeriodic(false) {};
 
    /**---------------------------------------------------------------------------------------
       
       Destructor
       
          --------------------------------------------------------------------------------------- */
 
    ~AmoebaReferenceStretchBendForce() {};

    /**---------------------------------------------------------------------------------------

       Set the force to use periodic boundary conditions.
      
       @param vectors    the vectors defining the periodic box
      
       --------------------------------------------------------------------------------------- */
      
    void setPeriodic(OpenMM::RealVec* vectors);

     /**---------------------------------------------------------------------------------------
     
        Calculate Amoeba stretch bend ixns (force and energy)
     
        @param numBonds                number of angles
        @param posData                 particle positions
        @param particle1               particle 1 indices
        @param particle2               particle 2 indices
        @param particle3               particle 3 indices
        @param lengthABParameters      ideal AB bond length 
        @param lengthCBParameters      ideal CB bond length 
        @param angle                   ideal angle 
        @param kQuadratic              force constant
        @param forces                  output force vector
     
        @return total energy

     
        --------------------------------------------------------------------------------------- */

    RealOpenMM calculateForceAndEnergy(int numAngles, std::vector<OpenMM::RealVec>& posData,
                                       const std::vector<int>& particle1,
                                       const std::vector<int>&  particle2,
                                       const std::vector<int>&  particle3,
                                       const std::vector<RealOpenMM>& lengthABParameters,
                                       const std::vector<RealOpenMM>& lengthCBParameters,
                                       const std::vector<RealOpenMM>&  angle,
                                       const std::vector<RealOpenMM>&  k1Quadratic,
                                       const std::vector<RealOpenMM>&  k2Quadratic,
                                       std::vector<OpenMM::RealVec>& forceData) const;


private:

    bool usePeriodic;
    RealVec boxVectors[3];

    /**---------------------------------------------------------------------------------------
    
       Calculate Amoeba stretch bend angle ixn (force and energy)
    
       @param positionAtomA           Cartesian coordinates of atom A
       @param positionAtomB           Cartesian coordinates of atom B
       @param positionAtomC           Cartesian coordinates of atom C
       @param lengthAB                ideal AB bondlength
       @param lengthCB                ideal CB bondlength
       @param idealAngle              ideal angle
       @param k1Parameter             k for distance A-B * angle A-B-C
       @param k2Parameter             k for distance B-C * angle A-B-C
       @param forces                  force vector
    
       @return energy
    
       --------------------------------------------------------------------------------------- */

    RealOpenMM calculateStretchBendIxn(const OpenMM::RealVec& positionAtomA, const OpenMM::RealVec& positionAtomB,
                                       const OpenMM::RealVec& positionAtomC,
                                       RealOpenMM lengthAB,      RealOpenMM lengthCB,
                                       RealOpenMM idealAngle,    RealOpenMM k1Parameter,
                                       RealOpenMM k2Parameter,   OpenMM::RealVec* forces) const;
 
};

} // namespace OpenMM

#endif // _AmoebaReferenceStretchBendForce___
