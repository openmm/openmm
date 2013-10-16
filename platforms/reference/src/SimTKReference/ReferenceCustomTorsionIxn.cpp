/* Portions copyright (c) 2010-2013 Stanford University and Simbios.
 * Contributors: Peter Eastman
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
#include "ReferenceCustomTorsionIxn.h"
#include "ReferenceForce.h"

using namespace std;
using namespace OpenMM;

/**---------------------------------------------------------------------------------------

   ReferenceCustomTorsionIxn constructor

   --------------------------------------------------------------------------------------- */

ReferenceCustomTorsionIxn::ReferenceCustomTorsionIxn(const Lepton::CompiledExpression& energyExpression,
        const Lepton::CompiledExpression& forceExpression, const vector<string>& parameterNames, map<string, double> globalParameters) :
        energyExpression(energyExpression), forceExpression(forceExpression) {

    energyTheta = ReferenceForce::getVariablePointer(this->energyExpression, "theta");
    forceTheta = ReferenceForce::getVariablePointer(this->forceExpression, "theta");
    numParameters = parameterNames.size();
    for (int i = 0; i < (int) numParameters; i++) {
        energyParams.push_back(ReferenceForce::getVariablePointer(this->energyExpression, parameterNames[i]));
        forceParams.push_back(ReferenceForce::getVariablePointer(this->forceExpression, parameterNames[i]));
    }
    for (map<string, double>::const_iterator iter = globalParameters.begin(); iter != globalParameters.end(); ++iter) {
        ReferenceForce::setVariable(ReferenceForce::getVariablePointer(this->energyExpression, iter->first), iter->second);
        ReferenceForce::setVariable(ReferenceForce::getVariablePointer(this->forceExpression, iter->first), iter->second);
    }
}

/**---------------------------------------------------------------------------------------

   ReferenceCustomTorsionIxn destructor

   --------------------------------------------------------------------------------------- */

ReferenceCustomTorsionIxn::~ReferenceCustomTorsionIxn( ){

   // ---------------------------------------------------------------------------------------

   // static const char* methodName = "\nReferenceCustomTorsionIxn::~ReferenceCustomTorsionIxn";

   // ---------------------------------------------------------------------------------------

}

/**---------------------------------------------------------------------------------------

   Calculate Custom Torsion Ixn

   @param atomIndices      atom indices of atom participating in bond
   @param atomCoordinates  atom coordinates
   @param parameters       parameters values
   @param forces           force array (forces added to input values)
   @param totalEnergy      if not null, the energy will be added to this

   --------------------------------------------------------------------------------------- */

void ReferenceCustomTorsionIxn::calculateBondIxn( int* atomIndices,
                                                vector<RealVec>& atomCoordinates,
                                                RealOpenMM* parameters,
                                                vector<RealVec>& forces,
                                                RealOpenMM* totalEnergy ) const {

   static const std::string methodName = "\nReferenceCustomTorsionIxn::calculateTorsionIxn";

   static const RealOpenMM zero        = 0.0;
   static const RealOpenMM one         = 1.0;

   RealOpenMM deltaR[3][ReferenceForce::LastDeltaRIndex];
   for (int i = 0; i < numParameters; i++) {
       ReferenceForce::setVariable(energyParams[i], parameters[i]);
       ReferenceForce::setVariable(forceParams[i], parameters[i]);
   }

   // ---------------------------------------------------------------------------------------

   // get deltaR, R2, and R between three pairs of atoms: [j,i], [j,k], [l,k]

   int atomAIndex = atomIndices[0];
   int atomBIndex = atomIndices[1];
   int atomCIndex = atomIndices[2];
   int atomDIndex = atomIndices[3];
   ReferenceForce::getDeltaR(atomCoordinates[atomBIndex], atomCoordinates[atomAIndex], deltaR[0]);
   ReferenceForce::getDeltaR(atomCoordinates[atomBIndex], atomCoordinates[atomCIndex], deltaR[1]);
   ReferenceForce::getDeltaR(atomCoordinates[atomDIndex], atomCoordinates[atomCIndex], deltaR[2]);

   // Visual Studio complains if crossProduct declared as 'crossProduct[2][3]'

   RealOpenMM crossProductMemory[6];
   RealOpenMM* crossProduct[2];
   crossProduct[0] = crossProductMemory;
   crossProduct[1] = crossProductMemory + 3;

   // get dihedral angle

   RealOpenMM dotDihedral;
   RealOpenMM signOfAngle;
   RealOpenMM angle = getDihedralAngleBetweenThreeVectors(deltaR[0], deltaR[1], deltaR[2], crossProduct, &dotDihedral, deltaR[0], &signOfAngle, 1);
   ReferenceForce::setVariable(energyTheta, angle);
   ReferenceForce::setVariable(forceTheta, angle);

   // evaluate delta angle, dE/d(angle)

   RealOpenMM dEdAngle = (RealOpenMM) forceExpression.evaluate();

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
       *totalEnergy += (RealOpenMM) energyExpression.evaluate();
}

