
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

#include <string.h>
#include <sstream>

#include "SimTKOpenMMUtilities.h"
#include "ReferenceForce.h"
#include "ReferenceBondIxn.h"

using std::vector;
using namespace OpenMM;

/**---------------------------------------------------------------------------------------

   ReferenceBondIxn constructor

   --------------------------------------------------------------------------------------- */

ReferenceBondIxn::ReferenceBondIxn() {
}

/**---------------------------------------------------------------------------------------

   ReferenceBondIxn destructor

   --------------------------------------------------------------------------------------- */

ReferenceBondIxn::~ReferenceBondIxn() {
}

/**---------------------------------------------------------------------------------------
      
   Calculate Bond Ixn -- virtual method -- does nothing
      
   @param atomIndices      bond indices
   @param atomCoordinates  atom coordinates
   @param parameters       parameters
   @param forces           force array (forces added)
   @param totalEnergy      if not null, the energy will be added to this
      
   --------------------------------------------------------------------------------------- */
     
   void ReferenceBondIxn::calculateBondIxn(vector<int>& atomIndices, vector<Vec3>& atomCoordinates,
                                           vector<double>& parameters, vector<Vec3>& forces,
                                           double* totalEnergy, double* energyParamDerivs) {
}
 
/**---------------------------------------------------------------------------------------

   Get normed dot product between two vectors

   Do computation in double?

   @param  vector1            first vector
   @param  vector2            second vector
   @param  hasREntry          if set, then vector1[ReferenceForce::RIndex] = norm of vector
                              defaults to 0 (i.e., R unavailable)

   @return dot product

   --------------------------------------------------------------------------------------- */

double ReferenceBondIxn::getNormedDotProduct(double* vector1, double* vector2,
                                             int hasREntry = 0) {

   double dotProduct = DOT3(vector1, vector2);
   if (dotProduct != 0.0) {
      if (hasREntry) {
         dotProduct       /= (vector1[ReferenceForce::RIndex]*vector2[ReferenceForce::RIndex]);
      } else {
         double norm1  = DOT3(vector1, vector1);
         double norm2  = DOT3(vector2, vector2);
         dotProduct   /= sqrt(norm1*norm2);
      }
   }
      
   // clamp dot product to [-1,1]

   if (dotProduct > 1.0) {
      dotProduct = 1.0;
   } else if (dotProduct < -1.0) {
      dotProduct = -1.0;
   }

   return dotProduct;

}

/**---------------------------------------------------------------------------------------

   Get angle between two vectors

   @param  vector1            first vector
   @param  vector2            second vector
   @param  outputDotProduct   output cosine of angle between two vectors (optional)
   @param  hasREntry          if set, then vector1[ReferenceForce::RIndex] = norm of vector
                              defaults to 0 -> R unavailable

   @return cosine of angles in radians

   --------------------------------------------------------------------------------------- */

double ReferenceBondIxn::getAngleBetweenTwoVectors(double* vector1, double* vector2, 
                                                   double* outputDotProduct = NULL,
                                                   int hasREntry = 0) {

    // get dot product betweenn vectors and then angle

   double dotProduct = getNormedDotProduct(vector1, vector2, hasREntry);

   double angle;
   if (dotProduct > 0.99 || dotProduct < -0.99) {
       // We're close to the singularity in acos(), so take the cross product and use asin() instead.

       double cross[3];
       SimTKOpenMMUtilities::crossProductVector3(vector1, vector2, cross);
       double scale = DOT3(vector1, vector1)*DOT3(vector2, vector2);
       angle = asin(sqrt(DOT3(cross, cross)/scale));
       if (dotProduct < 0.0)
           angle = M_PI-angle;
   } else {
      angle = acos(dotProduct);
   }

   if (outputDotProduct) {
      *outputDotProduct = dotProduct;
   }

   return angle;

}

/**---------------------------------------------------------------------------------------

   Get dihedral angle between three vectors

   @param  vector1            first vector
   @param  vector2            second vector
   @param  vector3            third vector
   @param  outputCrossProduct output cross product vectors
   @param  cosineOfAngle      cosine of angle (output)
   @param  signVector         vector to test sign (optional)
   @param  signOfAngle        sign of angle (output) (optional)
   @param  hasREntry          if set, then vector1[ReferenceForce::RIndex] = norm of vector
                              defaults to 0

   @return cosine of dihedral angle in radians

   --------------------------------------------------------------------------------------- */

double ReferenceBondIxn::getDihedralAngleBetweenThreeVectors(double*  vector1,
                                                             double*  vector2, 
                                                             double*  vector3, 
                                                             double** outputCrossProduct  = NULL, 
                                                             double*  cosineOfAngle       = NULL, 
                                                             double*  signVector          = NULL, 
                                                             double*  signOfAngle         = NULL, 
                                                             int          hasREntry = 0) {

   double   tempVectors[6]         = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };

   // get cross products between vectors and then angle between cross product vectors

   double* crossProduct[2];
   if (outputCrossProduct) {
      crossProduct[0] = outputCrossProduct[0];
      crossProduct[1] = outputCrossProduct[1];
   } else {
      crossProduct[0] = tempVectors;
      crossProduct[1] = tempVectors + 3;
   }
   
   SimTKOpenMMUtilities::crossProductVector3(vector1, vector2, crossProduct[0]);
   SimTKOpenMMUtilities::crossProductVector3(vector2, vector3, crossProduct[1]);

   double angle = getAngleBetweenTwoVectors(crossProduct[0], crossProduct[1], cosineOfAngle, 0);

   // take care of sign of angle

   if (signVector) {
      double dotProduct = DOT3(signVector, crossProduct[1]);
      double sign       = dotProduct < 0.0 ? -1.0 : 1.0; 
      if (signOfAngle) {
         *signOfAngle = sign;
      }
      angle *= sign;
   }

   return angle;

}
