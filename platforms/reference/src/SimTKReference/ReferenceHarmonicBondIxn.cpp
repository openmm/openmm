
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
#include "ReferenceHarmonicBondIxn.h"
#include "ReferenceForce.h"

using std::vector;
using OpenMM::RealVec;

/**---------------------------------------------------------------------------------------

   ReferenceHarmonicBondIxn constructor

   --------------------------------------------------------------------------------------- */

ReferenceHarmonicBondIxn::ReferenceHarmonicBondIxn( ){

   // ---------------------------------------------------------------------------------------

   // static const char* methodName = "\nReferenceHarmonicBondIxn::ReferenceHarmonicBondIxn";

   // ---------------------------------------------------------------------------------------

}

/**---------------------------------------------------------------------------------------

   ReferenceHarmonicBondIxn destructor

   --------------------------------------------------------------------------------------- */

ReferenceHarmonicBondIxn::~ReferenceHarmonicBondIxn( ){

   // ---------------------------------------------------------------------------------------

   // static const char* methodName = "\nReferenceHarmonicBondIxn::~ReferenceHarmonicBondIxn";

   // ---------------------------------------------------------------------------------------

}

/**---------------------------------------------------------------------------------------

   Calculate Harmonic Bond Ixn

   @param atomIndices      atom indices of atom participating in bond
   @param atomCoordinates  atom coordinates
   @param parameters       parameters: parameters[0] = ideal bond length
                                       parameters[1] = bond k
   @param forces           force array (forces added to input values)
   @param totalEnergy      if not null, the energy will be added to this

   --------------------------------------------------------------------------------------- */

void ReferenceHarmonicBondIxn::calculateBondIxn( int* atomIndices,
                                                vector<RealVec>& atomCoordinates,
                                                RealOpenMM* parameters,
                                                vector<RealVec>& forces,
                                                RealOpenMM* totalEnergy ) const {

   static const std::string methodName = "\nReferenceHarmonicBondIxn::calculateBondIxn";

   static const int twoI               = 2;

   static const RealOpenMM zero        = 0.0;
   static const RealOpenMM two         = 2.0;
   static const RealOpenMM half        = 0.5;

   RealOpenMM deltaR[ReferenceForce::LastDeltaRIndex];

   // ---------------------------------------------------------------------------------------

   // get deltaR, R2, and R between 2 atoms

   int atomAIndex = atomIndices[0];
   int atomBIndex = atomIndices[1];
   ReferenceForce::getDeltaR( atomCoordinates[atomAIndex], atomCoordinates[atomBIndex], deltaR );  

   // deltaIdeal = r - r_0

   RealOpenMM deltaIdeal      = deltaR[ReferenceForce::RIndex] - parameters[0];
   RealOpenMM deltaIdeal2     = deltaIdeal*deltaIdeal;

   RealOpenMM dEdR            = parameters[1]*deltaIdeal;

   // chain rule

   dEdR                       = deltaR[ReferenceForce::RIndex] > zero ? (dEdR/deltaR[ReferenceForce::RIndex]) : zero;

   forces[atomAIndex][0]     += dEdR*deltaR[ReferenceForce::XIndex];
   forces[atomAIndex][1]     += dEdR*deltaR[ReferenceForce::YIndex];
   forces[atomAIndex][2]     += dEdR*deltaR[ReferenceForce::ZIndex];

   forces[atomBIndex][0]     -= dEdR*deltaR[ReferenceForce::XIndex];
   forces[atomBIndex][1]     -= dEdR*deltaR[ReferenceForce::YIndex];
   forces[atomBIndex][2]     -= dEdR*deltaR[ReferenceForce::ZIndex];

   if (totalEnergy != NULL)
       *totalEnergy += half*parameters[1]*deltaIdeal2;
}
