
/* Portions copyright (c) 2009-2013 Stanford University and Simbios.
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
#include "ReferenceCustomBondIxn.h"
#include "ReferenceForce.h"

using namespace std;
using namespace OpenMM;

/**---------------------------------------------------------------------------------------

   ReferenceCustomBondIxn constructor

   --------------------------------------------------------------------------------------- */

ReferenceCustomBondIxn::ReferenceCustomBondIxn(const Lepton::CompiledExpression& energyExpression,
        const Lepton::CompiledExpression& forceExpression, const vector<string>& parameterNames, map<string, double> globalParameters) :
        energyExpression(energyExpression), forceExpression(forceExpression) {
    energyR = ReferenceForce::getVariablePointer(this->energyExpression, "r");
    forceR = ReferenceForce::getVariablePointer(this->forceExpression, "r");
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

   ReferenceCustomBondIxn destructor

   --------------------------------------------------------------------------------------- */

ReferenceCustomBondIxn::~ReferenceCustomBondIxn( ){

   // ---------------------------------------------------------------------------------------

   // static const char* methodName = "\nReferenceCustomBondIxn::~ReferenceCustomBondIxn";

   // ---------------------------------------------------------------------------------------

}

/**---------------------------------------------------------------------------------------

   Calculate Custom Bond Ixn

   @param atomIndices      atom indices of atom participating in bond
   @param atomCoordinates  atom coordinates
   @param parameters       parameters values
   @param forces           force array (forces added to input values)
   @param totalEnergy      if not null, the energy will be added to this

   --------------------------------------------------------------------------------------- */

void ReferenceCustomBondIxn::calculateBondIxn( int* atomIndices,
                                                vector<RealVec>& atomCoordinates,
                                                RealOpenMM* parameters,
                                                vector<RealVec>& forces,
                                                RealOpenMM* totalEnergy ) const {

   static const std::string methodName = "\nReferenceCustomBondIxn::calculateBondIxn";

   static const int twoI               = 2;

   static const RealOpenMM zero        = 0.0;
   static const RealOpenMM two         = 2.0;
   static const RealOpenMM half        = 0.5;

   RealOpenMM deltaR[ReferenceForce::LastDeltaRIndex];
   for (int i = 0; i < numParameters; i++) {
       ReferenceForce::setVariable(energyParams[i], parameters[i]);
       ReferenceForce::setVariable(forceParams[i], parameters[i]);
   }

   // ---------------------------------------------------------------------------------------

   // get deltaR, R2, and R between 2 atoms

   int atomAIndex = atomIndices[0];
   int atomBIndex = atomIndices[1];
   ReferenceForce::getDeltaR( atomCoordinates[atomAIndex], atomCoordinates[atomBIndex], deltaR );
   
   ReferenceForce::setVariable(energyR, deltaR[ReferenceForce::RIndex]);
   ReferenceForce::setVariable(forceR, deltaR[ReferenceForce::RIndex]);
   RealOpenMM dEdR            = (RealOpenMM) forceExpression.evaluate();
   dEdR                       = deltaR[ReferenceForce::RIndex] > zero ? (dEdR/deltaR[ReferenceForce::RIndex]) : zero;

   forces[atomAIndex][0]     += dEdR*deltaR[ReferenceForce::XIndex];
   forces[atomAIndex][1]     += dEdR*deltaR[ReferenceForce::YIndex];
   forces[atomAIndex][2]     += dEdR*deltaR[ReferenceForce::ZIndex];

   forces[atomBIndex][0]     -= dEdR*deltaR[ReferenceForce::XIndex];
   forces[atomBIndex][1]     -= dEdR*deltaR[ReferenceForce::YIndex];
   forces[atomBIndex][2]     -= dEdR*deltaR[ReferenceForce::ZIndex];

   if (totalEnergy != NULL)
       *totalEnergy += (RealOpenMM) energyExpression.evaluate();
}
