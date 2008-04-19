
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

#include <string.h>
#include <sstream>

#include "../SimTKUtilities/SimTKOpenMMCommon.h"
#include "../SimTKUtilities/SimTKOpenMMLog.h"
#include "../SimTKUtilities/SimTKOpenMMUtilities.h"
#include "ReferenceForce.h"
#include "ReferenceBondIxn.h"

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

   Update energy

   @param  energy               energy value to update 
   @param  energyByBond         ptr to energyByBond accumulator (may be null)
   @param  numberOfAtomIndices  number of atoms in bond
   @param  atomIndices          array of atom indices of size 'numberOfAtomIndices'
   @param  energyByAtom         array of energies by atom (may be null)

   @return ReferenceForce::DefaultReturn 

   --------------------------------------------------------------------------------------- */

int ReferenceBondIxn::updateEnergy( RealOpenMM energy, RealOpenMM* energyByBond,
                                    int numberOfAtomIndices, int* atomIndices, RealOpenMM* energyByAtom ) const {

   // ---------------------------------------------------------------------------------------

   // static const std::string methodName = "\nReferenceBondIxn::updateEnergy";

   // ---------------------------------------------------------------------------------------

   if( energyByBond ){
      *energyByBond += energy;
   }
   if( energyByAtom ){
      for( int ii = 0; ii < numberOfAtomIndices; ii++ ){
         energyByAtom[atomIndices[ii]] += energy;
      }
   }

   return ReferenceForce::DefaultReturn;

}

/**---------------------------------------------------------------------------------------
      
   Calculate Bond Ixn -- virtual method -- does nothing
      
   @param atomIndices      bond indices
   @param atomCoordinates  atom coordinates
   @param parameters       parameters
   @param forces           force array (forces added)
   @param energyByBond     bond energy 
   @param energy           atom energy
      
   --------------------------------------------------------------------------------------- */
     
   int ReferenceBondIxn::calculateBondIxn( int* atomIndices, RealOpenMM** atomCoordinates,
                                           RealOpenMM* parameters, RealOpenMM** forces,
                                           RealOpenMM* energyByBond, RealOpenMM* energyByAtom ) const {
   // ---------------------------------------------------------------------------------------

   // static const std::string methodName = "\nReferenceBondIxn::calculateBondIxn";

   // ---------------------------------------------------------------------------------------

   return ReferenceForce::DefaultReturn;

}
 
/**---------------------------------------------------------------------------------------

   Get normed dot product between two vectors

   Do computation in double?

   @param  vector1				first vector
   @param  vector2			   second vector
   @param  hasREntry          if set, then vector1[ReferenceForce::RIndex] = norm of vector
                              defaults to 0 (i.e., R unavailable)

   @return dot product

   --------------------------------------------------------------------------------------- */

RealOpenMM ReferenceBondIxn::getNormedDotProduct( RealOpenMM* vector1, RealOpenMM* vector2,
                                                  int hasREntry = 0 ) const {

   // ---------------------------------------------------------------------------------------

   // static const std::string methodName = "\nReferenceBondIxn::getNormedDotProduct";

   static const RealOpenMM zero        = 0.0;
   static const RealOpenMM one         = 1.0;

   // ---------------------------------------------------------------------------------------

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

   @param  vector1				first vector
   @param  vector2			   second vector
   @param  outputDotProduct   output cosine of angle between two vectors (optional)
   @param  hasREntry          if set, then vector1[ReferenceForce::RIndex] = norm of vector
                              defaults to 0 -> R unavailable

   @return cosine of angles in radians

   --------------------------------------------------------------------------------------- */

RealOpenMM ReferenceBondIxn::getAngleBetweenTwoVectors( RealOpenMM* vector1, RealOpenMM* vector2, 
                                                        RealOpenMM* outputDotProduct = NULL,
                                                        int hasREntry = 0 ) const {

   // ---------------------------------------------------------------------------------------

   // static const std::string methodName = "\nReferenceBondIxn::getAngle";

   static const RealOpenMM zero        = 0.0;
   static const RealOpenMM one         = 1.0;

   // ---------------------------------------------------------------------------------------

   // get dot product betweenn vectors and then angle

   RealOpenMM dotProduct = getNormedDotProduct( vector1, vector2, hasREntry );

   RealOpenMM angle;
   if( dotProduct >= one ){
      angle = zero;
   } else if( dotProduct <= -one ){
      angle = PI_M;
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

   @param  vector1				first vector
   @param  vector2			   second vector
   @param  vector3			   third vector
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
                                                                   int          hasREntry = 0 ) const {

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
