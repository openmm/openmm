
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

#ifndef __AmoebaReferenceOutOfPlaneBendForce_H__
#define __AmoebaReferenceOutOfPlaneBendForce_H__

#include "openmm/Vec3.h"
#include <vector>

namespace OpenMM {

class AmoebaReferenceOutOfPlaneBendForce {

public:
 
    /**---------------------------------------------------------------------------------------
       
       Constructor
       
       --------------------------------------------------------------------------------------- */
 
    AmoebaReferenceOutOfPlaneBendForce() : usePeriodic(false) {};
 
    /**---------------------------------------------------------------------------------------
       
       Destructor
       
          --------------------------------------------------------------------------------------- */
 
    ~AmoebaReferenceOutOfPlaneBendForce() {};
 
    /**---------------------------------------------------------------------------------------

       Set the force to use periodic boundary conditions.
      
       @param vectors    the vectors defining the periodic box
      
       --------------------------------------------------------------------------------------- */
      
    void setPeriodic(OpenMM::Vec3* vectors);

    /**---------------------------------------------------------------------------------------
     
        Calculate Amoeba out-of-plane-bend angle (force and energy)
     
        @param numOutOfPlaneBends      number of angles
        @param posData                 particle positions
        @param particle1               particle 1 indices
        @param particle2               particle 2 indices
        @param particle3               particle 3 indices
        @param particle4               particle 4 indices
        @param kAngle                  angle force constant
        @param angleCubic              cubic force parameter
        @param angleQuartic            quartic force parameter
        @param anglePentic             pentic force parameter
        @param angleSexic              sextic force parameter
        @param forces                  output force vector
     
        @return total energy

     
        --------------------------------------------------------------------------------------- */

    double calculateForceAndEnergy(int numOutOfPlaneBends, std::vector<OpenMM::Vec3>& posData,
                                   const std::vector<int>&  particle1,
                                   const std::vector<int>&  particle2,
                                   const std::vector<int>&  particle3,
                                   const std::vector<int>&  particle4,
                                   const std::vector<double>&  kAngle,
                                   double angleCubic,
                                   double angleQuartic,
                                   double anglePentic,
                                   double angleSextic,
                                   std::vector<OpenMM::Vec3>& forceData) const;

private:

    bool usePeriodic;
    Vec3 boxVectors[3];

    /**---------------------------------------------------------------------------------------
    
       Calculate Amoeba Out-Of-Plane-Bend ixn (force and energy)
    
       @param positionAtomA           Cartesian coordinates of atom A
       @param positionAtomB           Cartesian coordinates of atom B
       @param positionAtomC           Cartesian coordinates of atom C
       @param positionAtomD           Cartesian coordinates of atom D
       @param angleK                  quadratic angle force parameter
       @param angleCubic              cubic     angle force parameter
       @param angleQuartic            quartic   angle force parameter
       @param anglePentic             pentic    angle force parameter
       @param angleSextic             sextic    angle force parameter
       @param forces                  force vector
    
       @return energy
    
       --------------------------------------------------------------------------------------- */
    
    double calculateOutOfPlaneBendIxn(const OpenMM::Vec3& positionAtomA, const OpenMM::Vec3& positionAtomB,
                                     const OpenMM::Vec3& positionAtomC, const OpenMM::Vec3& positionAtomD,
                                     double angleK,
                                     double angleCubic,     double angleQuartic,
                                     double anglePentic,    double angleSextic,
                                     OpenMM::Vec3* forces) const;
};

} // namespace OpenMM

#endif // _AmoebaReferenceOutOfPlaneBendForce___
