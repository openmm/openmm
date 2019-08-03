
/* Portions copyright (c) 2009-2018 Stanford University and Simbios.
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

#include "SimTKOpenMMUtilities.h"
#include "ReferenceCustomExternalIxn.h"
#include "ReferenceForce.h"

using namespace std;
using namespace OpenMM;

/**---------------------------------------------------------------------------------------

   ReferenceCustomExternalIxn constructor

   --------------------------------------------------------------------------------------- */

ReferenceCustomExternalIxn::ReferenceCustomExternalIxn(const Lepton::CompiledExpression& energyExpression,
        const Lepton::CompiledExpression& forceExpressionX, const Lepton::CompiledExpression& forceExpressionY,
        const Lepton::CompiledExpression& forceExpressionZ, const vector<string>& parameterNames) :
        energyExpression(energyExpression), forceExpressionX(forceExpressionX), forceExpressionY(forceExpressionY),
        forceExpressionZ(forceExpressionZ) {

    energyX = ReferenceForce::getVariablePointer(this->energyExpression, "x");
    energyY = ReferenceForce::getVariablePointer(this->energyExpression, "y");
    energyZ = ReferenceForce::getVariablePointer(this->energyExpression, "z");
    forceXX = ReferenceForce::getVariablePointer(this->forceExpressionX, "x");
    forceXY = ReferenceForce::getVariablePointer(this->forceExpressionX, "y");
    forceXZ = ReferenceForce::getVariablePointer(this->forceExpressionX, "z");
    forceYX = ReferenceForce::getVariablePointer(this->forceExpressionY, "x");
    forceYY = ReferenceForce::getVariablePointer(this->forceExpressionY, "y");
    forceYZ = ReferenceForce::getVariablePointer(this->forceExpressionY, "z");
    forceZX = ReferenceForce::getVariablePointer(this->forceExpressionZ, "x");
    forceZY = ReferenceForce::getVariablePointer(this->forceExpressionZ, "y");
    forceZZ = ReferenceForce::getVariablePointer(this->forceExpressionZ, "z");
    numParameters = parameterNames.size();
    for (auto& param : parameterNames) {
        energyParams.push_back(ReferenceForce::getVariablePointer(this->energyExpression, param));
        forceXParams.push_back(ReferenceForce::getVariablePointer(this->forceExpressionX, param));
        forceYParams.push_back(ReferenceForce::getVariablePointer(this->forceExpressionY, param));
        forceZParams.push_back(ReferenceForce::getVariablePointer(this->forceExpressionZ, param));
    }
}

/**---------------------------------------------------------------------------------------

   ReferenceCustomExternalIxn destructor

   --------------------------------------------------------------------------------------- */

ReferenceCustomExternalIxn::~ReferenceCustomExternalIxn() {
}

void ReferenceCustomExternalIxn::setGlobalParameters(std::map<std::string, double> parameters) {
    for (auto& param : parameters) {
        ReferenceForce::setVariable(ReferenceForce::getVariablePointer(this->energyExpression, param.first), param.second);
        ReferenceForce::setVariable(ReferenceForce::getVariablePointer(this->forceExpressionX, param.first), param.second);
        ReferenceForce::setVariable(ReferenceForce::getVariablePointer(this->forceExpressionY, param.first), param.second);
        ReferenceForce::setVariable(ReferenceForce::getVariablePointer(this->forceExpressionZ, param.first), param.second);
    }
}

/**---------------------------------------------------------------------------------------

   Calculate Custom External Ixn

   @param atomIndex        the index of the atom to apply the force to
   @param atomCoordinates  atom coordinates
   @param parameters       parameters values
   @param forces           force array (forces added to input values)
   @param energy           energy is added to this

   --------------------------------------------------------------------------------------- */

void ReferenceCustomExternalIxn::calculateForce(int atomIndex,
                                                vector<Vec3>& atomCoordinates,
                                                vector<double>& parameters,
                                                vector<Vec3>& forces,
                                                double* energy) const {

   for (int i = 0; i < numParameters; i++) {
       ReferenceForce::setVariable(energyParams[i], parameters[i]);
       ReferenceForce::setVariable(forceXParams[i], parameters[i]);
       ReferenceForce::setVariable(forceYParams[i], parameters[i]);
       ReferenceForce::setVariable(forceZParams[i], parameters[i]);
   }
   ReferenceForce::setVariable(energyX, atomCoordinates[atomIndex][0]);
   ReferenceForce::setVariable(energyY, atomCoordinates[atomIndex][1]);
   ReferenceForce::setVariable(energyZ, atomCoordinates[atomIndex][2]);
   ReferenceForce::setVariable(forceXX, atomCoordinates[atomIndex][0]);
   ReferenceForce::setVariable(forceXY, atomCoordinates[atomIndex][1]);
   ReferenceForce::setVariable(forceXZ, atomCoordinates[atomIndex][2]);
   ReferenceForce::setVariable(forceYX, atomCoordinates[atomIndex][0]);
   ReferenceForce::setVariable(forceYY, atomCoordinates[atomIndex][1]);
   ReferenceForce::setVariable(forceYZ, atomCoordinates[atomIndex][2]);
   ReferenceForce::setVariable(forceZX, atomCoordinates[atomIndex][0]);
   ReferenceForce::setVariable(forceZY, atomCoordinates[atomIndex][1]);
   ReferenceForce::setVariable(forceZZ, atomCoordinates[atomIndex][2]);

   // ---------------------------------------------------------------------------------------

   forces[atomIndex][0] -= forceExpressionX.evaluate();
   forces[atomIndex][1] -= forceExpressionY.evaluate();
   forces[atomIndex][2] -= forceExpressionZ.evaluate();
   if (energy != NULL)
       *energy += energyExpression.evaluate();
}
