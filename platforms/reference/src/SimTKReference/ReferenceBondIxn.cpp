
/* Portions copyright (c) 2006-2009 Stanford University and Simbios.
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

#include "SimTKOpenMMCommon.h"
#include "SimTKOpenMMLog.h"
#include "SimTKOpenMMUtilities.h"
#include "ReferenceForce.h"
#include "ReferenceBondIxn.h"

using std::vector;
using OpenMM::RealVec;

/**---------------------------------------------------------------------------------------

   ReferenceBondIxn constructor

   --------------------------------------------------------------------------------------- */

ReferenceBondIxn::ReferenceBondIxn( ){

   // ---------------------------------------------------------------------------------------

   // static const char* methodName = "\nReferenceBondIxn::ReferenceBondIxn";

   // ---------------------------------------------------------------------------------------

}

/**---------------------------------------------------------------------------------------

   ReferenceBondIxn destructor

   --------------------------------------------------------------------------------------- */

ReferenceBondIxn::~ReferenceBondIxn( ){

   // ---------------------------------------------------------------------------------------

   // static const char* methodName = "\nReferenceBondIxn::~ReferenceBondIxn";

   // ---------------------------------------------------------------------------------------

}

/**---------------------------------------------------------------------------------------
      
   Calculate Bond Ixn -- virtual method -- does nothing
      
   @param atomIndices      bond indices
   @param atomCoordinates  atom coordinates
   @param parameters       parameters
   @param forces           force array (forces added)
   @param totalEnergy      if not null, the energy will be added to this
      
   --------------------------------------------------------------------------------------- */
     
   void ReferenceBondIxn::calculateBondIxn( int* atomIndices, vector<RealVec>& atomCoordinates,
                                           RealOpenMM* parameters, vector<RealVec>& forces,
                                           RealOpenMM* totalEnergy ) const {
   // ---------------------------------------------------------------------------------------

   // static const std::string methodName = "\nReferenceBondIxn::calculateBondIxn";

   // ---------------------------------------------------------------------------------------

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

RealOpenMM ReferenceBondIxn::getNormedDotProduct( RealOpenMM* vector1, RealOpenMM* vector2,
                                                  int hasREntry = 0 ) {

   // ---------------------------------------------------------------------------------------

   // static const std::string methodName = "\nReferenceBondIxn::getNormedDotProduct";

   static const RealOpenMM zero        = 0.0;
   static const RealOpenMM one         = 1.0;

   // ---------------------------------------------------------------------------------------

// for angles near pi, double is required due to the 'steepness' of acos()
// in this regime.
  
//#define USE_DOUBLE_FOR_NORMED_DOT_PRODUCT

#if defined USE_DOUBLE_FOR_NORMED_DOT_PRODUCT
   double v1D[3];
   double v2D[3];
   v1D[0]                = static_cast<double>( vector1[0] );
   v1D[1]                = static_cast<double>( vector1[1] );
   v1D[2]                = static_cast<double>( vector1[2] );

   v2D[0]                = static_cast<double>( vector2[0] );
   v2D[1]                = static_cast<double>( vector2[1] );
   v2D[2]                = static_cast<double>( vector2[2] );
   double dotProductD    = DOT3( v1D, v2D );
   if( dotProductD != 0.0 ){
      if( hasREntry ){
         dotProductD    /= ( static_cast<double>(vector1[ReferenceForce::RIndex])*static_cast<double>(vector2[ReferenceForce::RIndex]) );
      } else {
         double norm1    = DOT3( v1D, v1D );
         double norm2    = DOT3( v2D, v2D);
         dotProductD    /= sqrt( norm1*norm2 );
      }
   }
   RealOpenMM dotProduct = static_cast<RealOpenMM>(dotProductD);

#else      

   RealOpenMM dotProduct = DOT3( vector1, vector2 );
   if( dotProduct != zero ){
      if( hasREntry ){
         dotProduct       /= ( vector1[ReferenceForce::RIndex]*vector2[ReferenceForce::RIndex] );
      } else {
         RealOpenMM norm1  = DOT3( vector1, vector1 );
         RealOpenMM norm2  = DOT3( vector2, vector2 );
         dotProduct       /= SQRT( norm1*norm2 );
      }
   }

#endif
#undef USE_DOUBLE_FOR_NORMED_DOT_PRODUCT
      
   // clamp dot product to [-1,1]

   if( dotProduct > one ){
      dotProduct = one;
   } else if( dotProduct < -one ){
      dotProduct = -one;
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

RealOpenMM ReferenceBondIxn::getAngleBetweenTwoVectors( RealOpenMM* vector1, RealOpenMM* vector2, 
                                                        RealOpenMM* outputDotProduct = NULL,
                                                        int hasREntry = 0 ) {

   // ---------------------------------------------------------------------------------------

   // static const std::string methodName = "\nReferenceBondIxn::getAngle";

   static const RealOpenMM zero        = 0.0;
   static const RealOpenMM one         = 1.0;

   // ---------------------------------------------------------------------------------------

   // get dot product betweenn vectors and then angle

   RealOpenMM dotProduct = getNormedDotProduct( vector1, vector2, hasREntry );

   RealOpenMM angle;
   if (dotProduct > (RealOpenMM) 0.99 || dotProduct < (RealOpenMM) -0.99) {
       // We're close to the singularity in acos(), so take the cross product and use asin() instead.

       RealOpenMM cross[3];
       SimTKOpenMMUtilities::crossProductVector3(vector1, vector2, cross);
       RealOpenMM scale = DOT3(vector1, vector1)*DOT3(vector2, vector2);
       angle = ASIN(SQRT(DOT3(cross, cross)/scale));
       if (dotProduct < zero)
           angle = (RealOpenMM) (M_PI-angle);
   } else {
      angle = ACOS(dotProduct);
   }

   if( outputDotProduct ){
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

RealOpenMM ReferenceBondIxn::getDihedralAngleBetweenThreeVectors( RealOpenMM*  vector1,
                                                                  RealOpenMM*  vector2, 
                                                                  RealOpenMM*  vector3, 
                                                                  RealOpenMM** outputCrossProduct  = NULL, 
                                                                  RealOpenMM*  cosineOfAngle       = NULL, 
                                                                  RealOpenMM*  signVector          = NULL, 
                                                                  RealOpenMM*  signOfAngle         = NULL, 
                                                                   int          hasREntry = 0 ) {

   // ---------------------------------------------------------------------------------------

   // static const std::string methodName = "\nReferenceBondIxn::getDihedralAngleBetweenThreeVectors";

   static const RealOpenMM zero        = 0.0;
   static const RealOpenMM one         = 1.0;

   RealOpenMM   tempVectors[6]         = { zero, zero, zero, zero, zero, zero };

   // ---------------------------------------------------------------------------------------

   // get cross products between vectors and then angle between cross product vectors

   RealOpenMM* crossProduct[2];
   if( outputCrossProduct ){
      crossProduct[0] = outputCrossProduct[0];
      crossProduct[1] = outputCrossProduct[1];
   } else {
      crossProduct[0] = tempVectors;
      crossProduct[1] = tempVectors + 3;
   }
   
   SimTKOpenMMUtilities::crossProductVector3( vector1, vector2, crossProduct[0] );
   SimTKOpenMMUtilities::crossProductVector3( vector2, vector3, crossProduct[1] );

   RealOpenMM angle         = getAngleBetweenTwoVectors( crossProduct[0], crossProduct[1], cosineOfAngle, 0 );

   // take care of sign of angle

   if( signVector ){
      RealOpenMM dotProduct = DOT3( signVector, crossProduct[1] );
      RealOpenMM sign       = dotProduct < zero ? -one : one; 
      if( signOfAngle ){
         *signOfAngle = sign;
      }
      angle *= sign;
   }

   return angle;

}
