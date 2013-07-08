
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

#include "SimTKOpenMMCommon.h"
#include "SimTKOpenMMLog.h"
#include "SimTKOpenMMUtilities.h"
#include "ReferenceAngleBondIxn.h"
#include "ReferenceForce.h"

using std::vector;
using OpenMM::RealVec;

/**---------------------------------------------------------------------------------------

   ReferenceAngleBondIxn constructor

   --------------------------------------------------------------------------------------- */

ReferenceAngleBondIxn::ReferenceAngleBondIxn( ){

   // ---------------------------------------------------------------------------------------

   // static const char* methodName = "\nReferenceAngleBondIxn::ReferenceAngleBondIxn";

   // ---------------------------------------------------------------------------------------

}

/**---------------------------------------------------------------------------------------

   ReferenceAngleBondIxn destructor

   --------------------------------------------------------------------------------------- */

ReferenceAngleBondIxn::~ReferenceAngleBondIxn( ){

   // ---------------------------------------------------------------------------------------

   // static const char* methodName = "\nReferenceAngleBondIxn::~ReferenceAngleBondIxn";

   // ---------------------------------------------------------------------------------------

}

/**---------------------------------------------------------------------------------------

   Get dEdR and energy term for angle bond

   @param  cosine               cosine of angle
   @param  angleParameters      angleParameters: angleParameters[0] = angle in radians
                                                 angleParameters[1] = k (force constant)
   @param  dEdR                 output dEdR
   @param  energyTerm           output energyTerm

   --------------------------------------------------------------------------------------- */

void ReferenceAngleBondIxn::getPrefactorsGivenAngleCosine( RealOpenMM cosine, RealOpenMM* angleParameters,
                                                          RealOpenMM* dEdR, RealOpenMM* energyTerm ) const {

   // ---------------------------------------------------------------------------------------

   // static const std::string methodName = "\nReferenceAngleBondIxn::getPrefactorsGivenAngleCosine";

   static const RealOpenMM zero        = 0.0;
   static const RealOpenMM one         = 1.0;
   static const RealOpenMM half        = 0.5;

   // ---------------------------------------------------------------------------------------

   RealOpenMM angle;
   if( cosine >= one ){
      angle = zero;
   } else if( cosine <= -one ){
      angle = PI_M;
   } else {
      angle = ACOS(cosine);
   }
   RealOpenMM deltaIdeal         = angle - angleParameters[0];
   RealOpenMM deltaIdeal2        = deltaIdeal*deltaIdeal;

  *dEdR                          = angleParameters[1]*deltaIdeal;
  *energyTerm                    = half*angleParameters[1]*deltaIdeal2;

}

/**---------------------------------------------------------------------------------------

   Calculate Angle Bond ixn

   @param atomIndices      two bond indices
   @param atomCoordinates  atom coordinates
   @param parameters       parameters: parameters[0] = ideal bond length
                                       parameters[1] = bond k (includes factor of 2)
   @param forces           force array (forces added)
   @param totalEnergy      if not null, the energy will be added to this

   --------------------------------------------------------------------------------------- */

void ReferenceAngleBondIxn::calculateBondIxn( int* atomIndices,
                                             vector<RealVec>& atomCoordinates,
                                             RealOpenMM* parameters,
                                             vector<RealVec>& forces,
                                             RealOpenMM* totalEnergy ) const {

   // constants -- reduce Visual Studio warnings regarding conversions between float & double

   static const RealOpenMM zero        =  0.0;
   static const RealOpenMM one         =  1.0;
   static const RealOpenMM two         =  2.0;
   static const RealOpenMM three       =  3.0;
   static const RealOpenMM oneM        = -1.0;

   static const int threeI             = 3;

   static const int LastAtomIndex      = 3;

   RealOpenMM deltaR[2][ReferenceForce::LastDeltaRIndex];

   // ---------------------------------------------------------------------------------------

   // get deltaR, R2, and R between 2 atoms

   int atomAIndex = atomIndices[0];
   int atomBIndex = atomIndices[1];
   int atomCIndex = atomIndices[2];
   ReferenceForce::getDeltaR( atomCoordinates[atomAIndex], atomCoordinates[atomBIndex], deltaR[0] );  
   ReferenceForce::getDeltaR( atomCoordinates[atomCIndex], atomCoordinates[atomBIndex], deltaR[1] );  

   RealOpenMM pVector[threeI];
   SimTKOpenMMUtilities::crossProductVector3( deltaR[0], deltaR[1], pVector );
   RealOpenMM rp              = DOT3( pVector, pVector );
   rp                         = SQRT( rp );
   if( rp < 1.0e-06 ){
      rp = (RealOpenMM) 1.0e-06;
   }   
   RealOpenMM dot             = DOT3( deltaR[0], deltaR[1] );
   RealOpenMM cosine          = dot/SQRT( (deltaR[0][ReferenceForce::R2Index]*deltaR[1][ReferenceForce::R2Index]) );

   RealOpenMM dEdR;
   RealOpenMM energy;
   getPrefactorsGivenAngleCosine( cosine, parameters, &dEdR, &energy );

   RealOpenMM termA           =  dEdR/(deltaR[0][ReferenceForce::R2Index]*rp);
   RealOpenMM termC           = -dEdR/(deltaR[1][ReferenceForce::R2Index]*rp);

   RealOpenMM deltaCrossP[LastAtomIndex][threeI];
   SimTKOpenMMUtilities::crossProductVector3( deltaR[0], pVector, deltaCrossP[0] );
   SimTKOpenMMUtilities::crossProductVector3( deltaR[1], pVector, deltaCrossP[2] );

   for( int ii = 0; ii < threeI; ii++ ){
      deltaCrossP[0][ii] *= termA;
      deltaCrossP[2][ii] *= termC;
      deltaCrossP[1][ii]  = oneM*(deltaCrossP[0][ii] + deltaCrossP[2][ii]);
   }   

   // accumulate forces
 
   for( int jj = 0; jj < LastAtomIndex; jj++ ){
      for( int ii = 0; ii < threeI; ii++ ){
         forces[atomIndices[jj]][ii] += deltaCrossP[jj][ii];
      }   
   }   

   // accumulate energies

   if (totalEnergy != NULL)
       *totalEnergy += energy;
}
