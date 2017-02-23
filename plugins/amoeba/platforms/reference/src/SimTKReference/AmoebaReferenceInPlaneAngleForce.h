
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

#ifndef __AmoebaReferenceInPlaneAngleForce_H__
#define __AmoebaReferenceInPlaneAngleForce_H__

#include "openmm/Vec3.h"
#include <vector>

namespace OpenMM {

class AmoebaReferenceInPlaneAngleForce {

public:
 
    /**---------------------------------------------------------------------------------------
       
       Constructor
       
       --------------------------------------------------------------------------------------- */
 
    AmoebaReferenceInPlaneAngleForce() : usePeriodic(false) {};
 
    /**---------------------------------------------------------------------------------------
       
       Destructor
       
          --------------------------------------------------------------------------------------- */
 
    ~AmoebaReferenceInPlaneAngleForce() {};
 
    /**---------------------------------------------------------------------------------------

       Set the force to use periodic boundary conditions.
      
       @param vectors    the vectors defining the periodic box
      
       --------------------------------------------------------------------------------------- */
      
    void setPeriodic(OpenMM::Vec3* vectors);

     /**---------------------------------------------------------------------------------------
     
        Calculate Amoeba in-plane angle ixns (force and energy)
     
        @param numBonds                number of angles
        @param posData                 particle positions
        @param particle1               particle 1 indices
        @param particle2               particle 2 indices
        @param particle3               particle 3 indices
        @param particle4               particle 4 indices
        @param angle                   ideal angle 
        @param angleK                  angle force constant
        @param angleCubic              cubic force parameter
        @param angleQuartic            quartic force parameter
        @param anglePentic             pentic force parameter
        @param angleSexic              sextic force parameter
        @param forces                  output force vector
     
        @return total energy

     
        --------------------------------------------------------------------------------------- */

    double calculateForceAndEnergy(int numAngles, std::vector<OpenMM::Vec3>& posData,
                                   const std::vector<int>& particle1,
                                   const std::vector<int>&  particle2,
                                   const std::vector<int>&  particle3,
                                   const std::vector<int>&  particle4,
                                   const std::vector<double>& angle,
                                   const std::vector<double>& kQuadratic,
                                   double globalAngleCubic,
                                   double globalAngleQuartic,
                                   double globalAnglePentic,
                                   double globalAngleSextic,
                                   std::vector<OpenMM::Vec3>& forceData) const;

private:

    bool usePeriodic;
    Vec3 boxVectors[3];

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
    
    double getPrefactorsGivenAngleCosine(double cosine, double idealAngle, double angleK,
                                         double angleCubic,     double angleQuartic,
                                         double anglePentic,    double angleSextic,
                                         double* dEdR) const;

    /**---------------------------------------------------------------------------------------
    
       Calculate Amoeba angle ixn (force and energy)
    
       @param positionAtomA           Cartesian coordinates of atom A
       @param positionAtomB           Cartesian coordinates of atom B
       @param positionAtomC           Cartesian coordinates of atom C
       @param positionAtomD           Cartesian coordinates of atom D
       @param angle                   angle
       @param angleK                  quadratic angle force parameter
       @param angleCubic              cubic     angle force parameter
       @param angleQuartic            quartic   angle force parameter
       @param anglePentic             pentic    angle force parameter
       @param angleSextic             sextic    angle force parameter
       @param forces                  force vector
    
       @return energy
    
       --------------------------------------------------------------------------------------- */
    
    double calculateAngleIxn(const OpenMM::Vec3& positionAtomA, const OpenMM::Vec3& positionAtomB,
                             const OpenMM::Vec3& positionAtomC, const OpenMM::Vec3& positionAtomD,
                             double angle,          double angleK,
                             double angleCubic,     double angleQuartic,
                             double anglePentic,    double angleSextic,
                             OpenMM::Vec3* forces) const;
         
};

} // namespace OpenMM

#endif // _AmoebaReferenceInPlaneAngleForce___
