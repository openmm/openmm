
/* Portions copyright (c) 2009 Stanford University and Simbios.
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

#include "../SimTKUtilities/SimTKOpenMMCommon.h"
#include "../SimTKUtilities/SimTKOpenMMLog.h"
#include "../SimTKUtilities/SimTKOpenMMUtilities.h"
#include "ReferenceCustomBondIxn.h"
#include "ReferenceForce.h"

using namespace std;

/**---------------------------------------------------------------------------------------

   ReferenceCustomBondIxn constructor

   --------------------------------------------------------------------------------------- */

ReferenceCustomBondIxn::ReferenceCustomBondIxn(const Lepton::ExpressionProgram& energyExpression,
        const Lepton::ExpressionProgram& forceExpression, const vector<string>& parameterNames, map<string, double> globalParameters) :
        energyExpression(energyExpression), forceExpression(forceExpression), paramNames(parameterNames), globalParameters(globalParameters) {

   // ---------------------------------------------------------------------------------------

   // static const char* methodName = "\nReferenceCustomBondIxn::ReferenceCustomBondIxn";

   // ---------------------------------------------------------------------------------------

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
   @param energiesByBond   energies by bond: energiesByBond[bondIndex]
   @param energiesByAtom   energies by atom: energiesByAtom[atomIndex]

   @return ReferenceForce::DefaultReturn;

   --------------------------------------------------------------------------------------- */

int ReferenceCustomBondIxn::calculateBondIxn( int* atomIndices,
                                                RealOpenMM** atomCoordinates,
                                                RealOpenMM* parameters,
                                                RealOpenMM** forces,
                                                RealOpenMM* energiesByBond,
                                                RealOpenMM* energiesByAtom ) const {

   // ---------------------------------------------------------------------------------------

   // static const char* methodName = "\nReferenceCustomBondIxn::calculateBondIxn";

   // ---------------------------------------------------------------------------------------

   static const std::string methodName = "\nReferenceCustomBondIxn::calculateBondIxn";

   static const int twoI               = 2;

   static const RealOpenMM zero        = 0.0;
   static const RealOpenMM two         = 2.0;
   static const RealOpenMM half        = 0.5;

   RealOpenMM deltaR[ReferenceForce::LastDeltaRIndex];
   map<string, double> variables = globalParameters;
   for (int i = 0; i < (int) paramNames.size(); ++i)
       variables[paramNames[i]] = parameters[i];

   // ---------------------------------------------------------------------------------------

   // get deltaR, R2, and R between 2 atoms

   int atomAIndex = atomIndices[0];
   int atomBIndex = atomIndices[1];
   ReferenceForce::getDeltaR( atomCoordinates[atomAIndex], atomCoordinates[atomBIndex], deltaR );
   variables["r"]             = deltaR[ReferenceForce::RIndex];
   RealOpenMM dEdR            = forceExpression.evaluate(variables);
   dEdR                       = deltaR[ReferenceForce::RIndex] > zero ? (dEdR/deltaR[ReferenceForce::RIndex]) : zero;

   forces[atomAIndex][0]     += dEdR*deltaR[ReferenceForce::XIndex];
   forces[atomAIndex][1]     += dEdR*deltaR[ReferenceForce::YIndex];
   forces[atomAIndex][2]     += dEdR*deltaR[ReferenceForce::ZIndex];

   forces[atomBIndex][0]     -= dEdR*deltaR[ReferenceForce::XIndex];
   forces[atomBIndex][1]     -= dEdR*deltaR[ReferenceForce::YIndex];
   forces[atomBIndex][2]     -= dEdR*deltaR[ReferenceForce::ZIndex];

   RealOpenMM energy          = energyExpression.evaluate(variables);
   updateEnergy( energy, energiesByBond, twoI, atomIndices, energiesByAtom );

   return ReferenceForce::DefaultReturn;
}
