
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
#include "ReferenceProperDihedralBond.h"
#include "ReferenceForce.h"

using std::vector;
using OpenMM::RealVec;

/**---------------------------------------------------------------------------------------

   ReferenceProperDihedralBond constructor

   --------------------------------------------------------------------------------------- */

ReferenceProperDihedralBond::ReferenceProperDihedralBond( ){

   // ---------------------------------------------------------------------------------------

   // static const char* methodName = "\nReferenceProperDihedralBond::ReferenceProperDihedralBond";

   // ---------------------------------------------------------------------------------------

}

/**---------------------------------------------------------------------------------------

   ReferenceProperDihedralBond destructor

   --------------------------------------------------------------------------------------- */

ReferenceProperDihedralBond::~ReferenceProperDihedralBond( ){

   // ---------------------------------------------------------------------------------------

   // static const char* methodName = "\nReferenceProperDihedralBond::~ReferenceProperDihedralBond";

   // ---------------------------------------------------------------------------------------

}

/**---------------------------------------------------------------------------------------

   Calculate proper dihedral bond ixn

   @param atomIndices      atom indices of 4 atoms in bond
   @param atomCoordinates  atom coordinates
   @param parameters       3 parameters: parameters[0] = k
                                         parameters[1] = ideal bond angle in radians
                                         parameters[2] = multiplicity
   @param forces           force array (forces added to current values)
   @param totalEnergy      if not null, the energy will be added to this

   --------------------------------------------------------------------------------------- */

void ReferenceProperDihedralBond::calculateBondIxn( int* atomIndices,
                                                   vector<RealVec>& atomCoordinates,
                                                   RealOpenMM* parameters,
                                                   vector<RealVec>& forces,
                                                   RealOpenMM* totalEnergy ) const {

   static const std::string methodName = "\nReferenceProperDihedralBond::calculateBondIxn";

   // constants -- reduce Visual Studio warnings regarding conversions between float & double

   static const RealOpenMM zero        =  0.0;
   static const RealOpenMM one         =  1.0;
   static const RealOpenMM two         =  2.0;
   static const RealOpenMM three       =  3.0;
   static const RealOpenMM oneM        = -1.0;

   static const int threeI             = 3;

   static const int LastAtomIndex      = 4;

   RealOpenMM deltaR[3][ReferenceForce::LastDeltaRIndex];

   RealOpenMM crossProductMemory[6];

   // ---------------------------------------------------------------------------------------

   // get deltaR, R2, and R between three pairs of atoms: [j,i], [j,k], [l,k]

   int atomAIndex = atomIndices[0];
   int atomBIndex = atomIndices[1];
   int atomCIndex = atomIndices[2];
   int atomDIndex = atomIndices[3];
   ReferenceForce::getDeltaR( atomCoordinates[atomBIndex], atomCoordinates[atomAIndex], deltaR[0] );  
   ReferenceForce::getDeltaR( atomCoordinates[atomBIndex], atomCoordinates[atomCIndex], deltaR[1] );  
   ReferenceForce::getDeltaR( atomCoordinates[atomDIndex], atomCoordinates[atomCIndex], deltaR[2] );  

   RealOpenMM dotDihedral;
   RealOpenMM signOfAngle;
   int hasREntry             = 1;

   // Visual Studio complains if crossProduct declared as 'crossProduct[2][3]'

   RealOpenMM* crossProduct[2];
   crossProduct[0]           = crossProductMemory;
   crossProduct[1]           = crossProductMemory + 3;

   // get dihedral angle

   RealOpenMM dihedralAngle  =  getDihedralAngleBetweenThreeVectors( deltaR[0], deltaR[1], deltaR[2],
                                                                     crossProduct, &dotDihedral, deltaR[0], 
                                                                     &signOfAngle, hasREntry );

   // evaluate delta angle, dE/d(angle) 

   RealOpenMM deltaAngle     = parameters[2]*dihedralAngle - parameters[1]; 
   RealOpenMM sinDeltaAngle  = SIN( deltaAngle );
   RealOpenMM dEdAngle       = -parameters[0]*parameters[2]*sinDeltaAngle;
   RealOpenMM energy         =  parameters[0]*(one + COS( deltaAngle ) );
   
   // compute force

   RealOpenMM internalF[4][3];
   RealOpenMM forceFactors[4];
   RealOpenMM normCross1         = DOT3( crossProduct[0], crossProduct[0] );
   RealOpenMM normBC             = deltaR[1][ReferenceForce::RIndex];
              forceFactors[0]    = (-dEdAngle*normBC)/normCross1;

   RealOpenMM normCross2         = DOT3( crossProduct[1], crossProduct[1] );
              forceFactors[3]    = (dEdAngle*normBC)/normCross2;
  
              forceFactors[1]    = DOT3( deltaR[0], deltaR[1] );
              forceFactors[1]   /= deltaR[1][ReferenceForce::R2Index];

              forceFactors[2]    = DOT3( deltaR[2], deltaR[1] );
              forceFactors[2]   /= deltaR[1][ReferenceForce::R2Index];

   for( int ii = 0; ii < 3; ii++ ){

      internalF[0][ii]  = forceFactors[0]*crossProduct[0][ii];
      internalF[3][ii]  = forceFactors[3]*crossProduct[1][ii];

      RealOpenMM s      = forceFactors[1]*internalF[0][ii] - forceFactors[2]*internalF[3][ii];

      internalF[1][ii]  = internalF[0][ii] - s;
      internalF[2][ii]  = internalF[3][ii] + s;
   }
  
   // accumulate forces

   for( int ii = 0; ii < 3; ii++ ){
      forces[atomAIndex][ii] += internalF[0][ii];
      forces[atomBIndex][ii] -= internalF[1][ii];
      forces[atomCIndex][ii] -= internalF[2][ii];
      forces[atomDIndex][ii] += internalF[3][ii];
   }

   // accumulate energies

   if (totalEnergy != NULL)
       *totalEnergy += energy;
}
