
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

#include "SimTKOpenMMCommon.h"
#include "SimTKOpenMMLog.h"
#include "SimTKOpenMMUtilities.h"
#include "ReferenceCustomExternalIxn.h"
#include "ReferenceForce.h"

using namespace std;
using namespace OpenMM;

/**---------------------------------------------------------------------------------------

   ReferenceCustomExternalIxn constructor

   --------------------------------------------------------------------------------------- */

ReferenceCustomExternalIxn::ReferenceCustomExternalIxn(const Lepton::ExpressionProgram& energyExpression,
        const Lepton::ExpressionProgram& forceExpressionX, const Lepton::ExpressionProgram& forceExpressionY,
        const Lepton::ExpressionProgram& forceExpressionZ, const vector<string>& parameterNames, map<string, double> globalParameters) :
        energyExpression(energyExpression), forceExpressionX(forceExpressionX), forceExpressionY(forceExpressionY),
        forceExpressionZ(forceExpressionZ), paramNames(parameterNames), globalParameters(globalParameters) {

   // ---------------------------------------------------------------------------------------

   // static const char* methodName = "\nReferenceCustomExternalIxn::ReferenceCustomExternalIxn";

   // ---------------------------------------------------------------------------------------

}

/**---------------------------------------------------------------------------------------

   ReferenceCustomExternalIxn destructor

   --------------------------------------------------------------------------------------- */

ReferenceCustomExternalIxn::~ReferenceCustomExternalIxn( ){

   // ---------------------------------------------------------------------------------------

   // static const char* methodName = "\nReferenceCustomExternalIxn::~ReferenceCustomExternalIxn";

   // ---------------------------------------------------------------------------------------

}

/**---------------------------------------------------------------------------------------

   Calculate Custom External Ixn

   @param atomIndex        the index of the atom to apply the force to
   @param atomCoordinates  atom coordinates
   @param parameters       parameters values
   @param forces           force array (forces added to input values)
   @param energy           energy is added to this

   --------------------------------------------------------------------------------------- */

void ReferenceCustomExternalIxn::calculateForce( int atomIndex,
                                                vector<RealVec>& atomCoordinates,
                                                RealOpenMM* parameters,
                                                vector<RealVec>& forces,
                                                RealOpenMM* energy ) const {

   static const std::string methodName = "\nReferenceCustomExternalIxn::calculateBondIxn";

   map<string, double> variables = globalParameters;
   for (int i = 0; i < (int) paramNames.size(); ++i)
       variables[paramNames[i]] = parameters[i];
   variables["x"] = atomCoordinates[atomIndex][0];
   variables["y"] = atomCoordinates[atomIndex][1];
   variables["z"] = atomCoordinates[atomIndex][2];

   // ---------------------------------------------------------------------------------------

   forces[atomIndex][0] -= (RealOpenMM) forceExpressionX.evaluate(variables);
   forces[atomIndex][1] -= (RealOpenMM) forceExpressionY.evaluate(variables);
   forces[atomIndex][2] -= (RealOpenMM) forceExpressionZ.evaluate(variables);
   if (energy != NULL)
       *energy += (RealOpenMM) energyExpression.evaluate(variables);
}
