
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
#include "ReferenceLJCoulomb14.h"
#include "ReferenceForce.h"

using std::vector;
using OpenMM::RealVec;

/**---------------------------------------------------------------------------------------

   ReferenceLJCoulomb14 constructor

   --------------------------------------------------------------------------------------- */

ReferenceLJCoulomb14::ReferenceLJCoulomb14( ) {

   // ---------------------------------------------------------------------------------------

   // static const char* methodName = "\nReferenceLJCoulomb14::ReferenceLJCoulomb14";

   // ---------------------------------------------------------------------------------------

}

/**---------------------------------------------------------------------------------------

   ReferenceLJCoulomb14 destructor

   --------------------------------------------------------------------------------------- */

ReferenceLJCoulomb14::~ReferenceLJCoulomb14( ){

   // ---------------------------------------------------------------------------------------

   // static const char* methodName = "\nReferenceLJCoulomb14::~ReferenceLJCoulomb14";

   // ---------------------------------------------------------------------------------------

}

/**---------------------------------------------------------------------------------------

   Calculate LJ 1-4 ixn

   @param atomIndices      atom indices of 4 atoms in bond
   @param atomCoordinates  atom coordinates
   @param parameters       three parameters:
                                        parameters[0]= (c12/c6)**1/6  (sigma)
                                        parameters[1]= c6*c6/c12      (4*epsilon)
                                        parameters[2]= epsfac*q1*q2
   @param forces           force array (forces added to current values)
   @param totalEnergy      if not null, the energy will be added to this

   --------------------------------------------------------------------------------------- */

void ReferenceLJCoulomb14::calculateBondIxn( int* atomIndices, vector<RealVec>& atomCoordinates,
                                     RealOpenMM* parameters, vector<RealVec>& forces,
                                     RealOpenMM* totalEnergy ) const {

   static const std::string methodName = "\nReferenceLJCoulomb14::calculateBondIxn";

   // constants -- reduce Visual Studio warnings regarding conversions between float & double

   static const RealOpenMM zero        =  0.0;
   static const RealOpenMM one         =  1.0;
   static const RealOpenMM two         =  2.0;
   static const RealOpenMM three       =  3.0;
   static const RealOpenMM six         =  6.0;
   static const RealOpenMM twelve      = 12.0;
   static const RealOpenMM oneM        = -1.0;

   static const int threeI             = 3;

   // number of parameters

   static const int numberOfParameters = 3;

   static const int LastAtomIndex      = 2;

   RealOpenMM deltaR[2][ReferenceForce::LastDeltaRIndex];

   // ---------------------------------------------------------------------------------------

   // get deltaR, R2, and R between 2 atoms

   int atomAIndex = atomIndices[0];
   int atomBIndex = atomIndices[1];
   ReferenceForce::getDeltaR( atomCoordinates[atomBIndex], atomCoordinates[atomAIndex], deltaR[0] );  

   RealOpenMM r2        = deltaR[0][ReferenceForce::R2Index];
   RealOpenMM inverseR  = one/(deltaR[0][ReferenceForce::RIndex]);
   RealOpenMM sig2      = inverseR*parameters[0];
              sig2     *= sig2;
   RealOpenMM sig6      = sig2*sig2*sig2;

   RealOpenMM dEdR      = parameters[1]*( twelve*sig6 - six )*sig6;
              dEdR     += (RealOpenMM) (ONE_4PI_EPS0*parameters[2]*inverseR);
              dEdR     *= inverseR*inverseR;

   // accumulate forces

   for( int ii = 0; ii < 3; ii++ ){
      RealOpenMM force        = dEdR*deltaR[0][ii];
      forces[atomAIndex][ii] += force;
      forces[atomBIndex][ii] -= force;
   }

   // accumulate energies

   if (totalEnergy != NULL)
       *totalEnergy += parameters[1]*( sig6 - one )*sig6 + (ONE_4PI_EPS0*parameters[2]*inverseR);
}
