
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

#ifndef __AmoebaReferenceBondForce_H__
#define __AmoebaReferenceBondForce_H__

#include "openmm/Vec3.h"
#include <vector>

namespace OpenMM {

class AmoebaReferenceBondForce {

public:
 
    /**---------------------------------------------------------------------------------------
       
        Constructor
       
        --------------------------------------------------------------------------------------- */
 
    AmoebaReferenceBondForce() : usePeriodic(false) {};
 
    /**---------------------------------------------------------------------------------------
       
        Destructor
       
        --------------------------------------------------------------------------------------- */
 
    ~AmoebaReferenceBondForce() {};

    /**---------------------------------------------------------------------------------------

       Set the force to use periodic boundary conditions.
      
       @param vectors    the vectors defining the periodic box
      
       --------------------------------------------------------------------------------------- */
      
    void setPeriodic(OpenMM::Vec3* vectors);
 
     /**---------------------------------------------------------------------------------------
     
        Calculate Amoeba bond ixns (force and energy)
     
        @param numBonds                number of bonds
        @param posData                 particle positions
        @param particle1               particle 1 indices
        @param particle2               particle 2 indices
        @param bondLength              bond length
        @param bondK                   bond force
        @param bondCubic               cubic bond force parameter
        @param bondQuartic             quartic bond force parameter
        @param forces                  output force vector
     
        @return total energy
     
        --------------------------------------------------------------------------------------- */
     
    double calculateForceAndEnergy(int numBonds, std::vector<OpenMM::Vec3>& posData,
                                   const std::vector<int>& particle1,
                                   const std::vector<int>&  particle2,
                                   const std::vector<double>& bondLength,
                                   const std::vector<double>& bondK,
                                   double bondCubic, double bondQuartic,
                                   std::vector<OpenMM::Vec3>& forceData) const;

private:

    bool usePeriodic;
    Vec3 boxVectors[3];

     /**---------------------------------------------------------------------------------------
     
        Calculate Amoeba bond ixns (force and energy)
     
        @param positionAtomA           Cartesian coordinates of atom A
        @param positionAtomB           Cartesian coordinates of atom B
        @param bondLength              bond length
        @param bondK                   bond force
        @param bondCubic               cubic bond force parameter
        @param bondQuartic             quartic bond force parameter
        @param forces                  force vector
     
        @return energy
     
        --------------------------------------------------------------------------------------- */
     
    double calculateBondIxn(const OpenMM::Vec3& positionAtomA, const OpenMM::Vec3& positionAtomB,
                            double bondLength, double bondK,
                            double bondCubic, double bondQuartic,
                            OpenMM::Vec3* forces) const;
     
};

} // namespace OpenMM

#endif // _AmoebaReferenceBondForce___
