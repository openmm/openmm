
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

#include "AmoebaReferenceForce.h"
#include <vector>

/**---------------------------------------------------------------------------------------

   Load delta of two vectors

   @param xVector      first vector
   @param yVector      second vector
   @param deltaR       output vector: y - x

   --------------------------------------------------------------------------------------- */

void AmoebaReferenceForce::loadDeltaR( const RealOpenMM* xVector, const RealOpenMM* yVector,
                                       std::vector<RealOpenMM>& deltaR ){


   // ---------------------------------------------------------------------------------------

   //static const std::string methodName = "AmoebaReferenceForce::loadDeltaR";

   // ---------------------------------------------------------------------------------------

    deltaR.resize(0);
    deltaR.push_back( yVector[0] - xVector[0] );
    deltaR.push_back( yVector[1] - xVector[1] );
    deltaR.push_back( yVector[2] - xVector[2] );
}

/**---------------------------------------------------------------------------------------

   Calculate norm squared of 3d vector

   @param inputVector    vector whose norm squared is to be computed

   @return norm squared

   --------------------------------------------------------------------------------------- */

RealOpenMM AmoebaReferenceForce::getNormSquared3( const std::vector<RealOpenMM>& inputVector ){

   // ---------------------------------------------------------------------------------------

   //static const std::string methodName = "AmoebaReferenceForce::getNorm3";

   // ---------------------------------------------------------------------------------------

   // get 3 norm

   return ( inputVector[0]*inputVector[0] + inputVector[1]*inputVector[1] + inputVector[2]*inputVector[2] );
}

/**---------------------------------------------------------------------------------------

   Calculate norm of 3d vector

   @param inputVector            vector whose norm is to be computed

   @return norm

   --------------------------------------------------------------------------------------- */

RealOpenMM AmoebaReferenceForce::getNorm3( const std::vector<RealOpenMM>& inputVector ){

   // ---------------------------------------------------------------------------------------

   //static const std::string methodName = "AmoebaReferenceForce::getNorm3";

   // ---------------------------------------------------------------------------------------

   // get 3 norm

   return SQRT( inputVector[0]*inputVector[0] + inputVector[1]*inputVector[1] + inputVector[2]*inputVector[2] );
}

/**---------------------------------------------------------------------------------------

   Calculate dot product of 3d vectors

   @param xVector   first vector
   @param yVector   second vector

   @return dot product

   --------------------------------------------------------------------------------------- */

RealOpenMM AmoebaReferenceForce::getDotProduct3( const std::vector<RealOpenMM>& xVector, const std::vector<RealOpenMM>& yVector ){

   // ---------------------------------------------------------------------------------------

   //static const std::string methodName = "AmoebaReferenceForce::getDotProduct3";

   // ---------------------------------------------------------------------------------------

   // get dot product

   return xVector[0]*yVector[0] + xVector[1]*yVector[1] + xVector[2]*yVector[2];
}

/**---------------------------------------------------------------------------------------

   Calculate z = x X y

   @param xVector      input vector
   @param yVector      input vector
   @param zVector      output vector: z = x X y

   --------------------------------------------------------------------------------------- */

void AmoebaReferenceForce::getCrossProduct( const std::vector<RealOpenMM>& xVector,
                                            const std::vector<RealOpenMM>& yVector,
                                            std::vector<RealOpenMM>& zVector ){

   // ---------------------------------------------------------------------------------------

   //static const std::string methodName = "AmoebaReferenceForce::getCrossProduct";

   // ---------------------------------------------------------------------------------------

   zVector[0]  = xVector[1]*yVector[2] - xVector[2]*yVector[1];
   zVector[1]  = xVector[2]*yVector[0] - xVector[0]*yVector[2];
   zVector[2]  = xVector[0]*yVector[1] - xVector[1]*yVector[0];

   return;
}

