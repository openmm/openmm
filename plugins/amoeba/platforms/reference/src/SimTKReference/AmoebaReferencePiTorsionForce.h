
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

#ifndef __AmoebaReferencePiTorsionForce_H__
#define __AmoebaReferencePiTorsionForce_H__

#include "openmm/Vec3.h"
#include <vector>

namespace OpenMM {

class AmoebaReferencePiTorsionForce {

public:
 
    /**---------------------------------------------------------------------------------------
       
       Constructor
       
       --------------------------------------------------------------------------------------- */
 
    AmoebaReferencePiTorsionForce() : usePeriodic(false) {};
 
    /**---------------------------------------------------------------------------------------
       
       Destructor
       
          --------------------------------------------------------------------------------------- */
 
    ~AmoebaReferencePiTorsionForce() {};
 
    /**---------------------------------------------------------------------------------------

       Set the force to use periodic boundary conditions.
      
       @param vectors    the vectors defining the periodic box
      
       --------------------------------------------------------------------------------------- */
      
    void setPeriodic(OpenMM::Vec3* vectors);

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

    double calculateForceAndEnergy(int numPiTorsions, std::vector<OpenMM::Vec3>& posData,
                                   const std::vector<int>&  particle1,
                                   const std::vector<int>&  particle2,
                                   const std::vector<int>&  particle3,
                                   const std::vector<int>&  particle4,
                                   const std::vector<int>&  particle5,
                                   const std::vector<int>&  particle6,
                                   const std::vector<double>& kTorsion,
                                   std::vector<OpenMM::Vec3>& forceData) const;


private:

    bool usePeriodic;
    Vec3 boxVectors[3];

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
    
    double calculatePiTorsionIxn(const OpenMM::Vec3& positionAtomA, const OpenMM::Vec3& positionAtomB,
                                 const OpenMM::Vec3& positionAtomC, const OpenMM::Vec3& positionAtomD,
                                 const OpenMM::Vec3& positionAtomE, const OpenMM::Vec3& positionAtomF,
                                 double kTorsion, OpenMM::Vec3* forces) const;
         
};

} // namespace OpenMM

#endif // _AmoebaReferencePiTorsionForce___
