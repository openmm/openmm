
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

#ifndef __ReferenceBondIxn_H__
#define __ReferenceBondIxn_H__

#include "openmm/Vec3.h"
#include "openmm/internal/windowsExport.h"
#include <vector>

namespace OpenMM {

class OPENMM_EXPORT ReferenceBondIxn {

   private:

   public:

      /**---------------------------------------------------------------------------------------
      
         Constructor
      
         --------------------------------------------------------------------------------------- */

       ReferenceBondIxn();

      /**---------------------------------------------------------------------------------------
      
         Destructor
      
         --------------------------------------------------------------------------------------- */

       ~ReferenceBondIxn();

      /**---------------------------------------------------------------------------------------
      
         Calculate Bond Ixn -- virtual method
      
         @param atomIndices      bond indices
         @param atomCoordinates  atom coordinates
         @param parameters       parameters
         @param forces           force array (forces added)
         @param totalEnergy      if not null, the energy will be added to this
      
         --------------------------------------------------------------------------------------- */
      
      virtual void calculateBondIxn(std::vector<int>& atomIndices, std::vector<OpenMM::Vec3>& atomCoordinates,
                                    std::vector<double>& parameters, std::vector<OpenMM::Vec3>& forces,
                                    double* totalEnergy, double* energyParamDerivs);
      
      /**---------------------------------------------------------------------------------------
      
         Get normed dot product between two vectors
      
         @param  vector1            first vector
         @param  vector2            second vector
         @param  hasREntry          if set, then vector1[ReferenceForce::RIndex] = norm of vector;
                                    defaults to 0
      
         @return dot product
      
         --------------------------------------------------------------------------------------- */
      
      static double getNormedDotProduct(double* vector1, double* vector2, int hasREntry);
      
      /**---------------------------------------------------------------------------------------
      
         Get angle between two vectors
      
         @param  vector1            first vector
         @param  vector2            second vector
         @param  outputDotProduct   cosine of angle between two vectors (optional)
         @param  hasREntry          if set, then vector1[ReferenceForce::R2Index] = square norm of vector;
                                    defaults to 0
      
         @return cosine of angles in radians
      
         --------------------------------------------------------------------------------------- */
      
      static double getAngleBetweenTwoVectors(double* vector1, double* vector2, 
                                              double* outputDotProduct, int hasREntry);
      
      /**---------------------------------------------------------------------------------------
      
         Get dihedral angle between two vectors
      
         @param  vector1            first vector
         @param  vector2            second vector
         @param  vector3            third vector
         @param  outputCrossProduct output cross product vectors
         @param  cosineOfAngle      cosine of angle (output)
         @param  signVector         vector to test sign (optional)
         @param  signOfAngle        sign of angle (output) (optional)
         @param  hasREntry          if set, then vector1[ReferenceForce::R2Index] = square norm of vector
                                    defaults to 0
      
         @return cosine of angles in radians
      
         --------------------------------------------------------------------------------------- */
      
      static double getDihedralAngleBetweenThreeVectors(double* vector1, double* vector2, 
                                                        double* vector3, double** outputCrossProduct, 
                                                        double* cosineOfAngle, double* signVector, 
                                                        double* signOfAngle, int hasREntry);
      
};

} // namespace OpenMM

#endif // __ReferenceBondIxn_H__
