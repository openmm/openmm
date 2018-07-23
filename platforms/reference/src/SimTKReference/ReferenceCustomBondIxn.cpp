
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
#include "ReferenceCustomBondIxn.h"
#include "ReferenceForce.h"

using namespace std;
using namespace OpenMM;

/**---------------------------------------------------------------------------------------

   ReferenceCustomBondIxn constructor

   --------------------------------------------------------------------------------------- */

ReferenceCustomBondIxn::ReferenceCustomBondIxn(const Lepton::CompiledExpression& energyExpression,
        const Lepton::CompiledExpression& forceExpression, const vector<string>& parameterNames,
        const vector<Lepton::CompiledExpression> energyParamDerivExpressions) :
        energyExpression(energyExpression), forceExpression(forceExpression), usePeriodic(false), energyParamDerivExpressions(energyParamDerivExpressions) {
    expressionSet.registerExpression(this->energyExpression);
    expressionSet.registerExpression(this->forceExpression);
    for (int i = 0; i < this->energyParamDerivExpressions.size(); i++)
        expressionSet.registerExpression(this->energyParamDerivExpressions[i]);
    rIndex = expressionSet.getVariableIndex("r");
    numParameters = parameterNames.size();
    for (auto& param : parameterNames)
        bondParamIndex.push_back(expressionSet.getVariableIndex(param));
}

/**---------------------------------------------------------------------------------------

   ReferenceCustomBondIxn destructor

   --------------------------------------------------------------------------------------- */

ReferenceCustomBondIxn::~ReferenceCustomBondIxn() {
}

void ReferenceCustomBondIxn::setPeriodic(OpenMM::Vec3* vectors) {
    usePeriodic = true;
    boxVectors[0] = vectors[0];
    boxVectors[1] = vectors[1];
    boxVectors[2] = vectors[2];
}

void ReferenceCustomBondIxn::setGlobalParameters(std::map<std::string, double> parameters) {
    for (auto& param : parameters)
        expressionSet.setVariable(expressionSet.getVariableIndex(param.first), param.second);
}

/**---------------------------------------------------------------------------------------

   Calculate Custom Bond Ixn

   @param atomIndices      atom indices of atom participating in bond
   @param atomCoordinates  atom coordinates
   @param parameters       parameters values
   @param forces           force array (forces added to input values)
   @param totalEnergy      if not null, the energy will be added to this

   --------------------------------------------------------------------------------------- */

void ReferenceCustomBondIxn::calculateBondIxn(vector<int>& atomIndices,
                                              vector<Vec3>& atomCoordinates,
                                              vector<double>& parameters,
                                              vector<Vec3>& forces,
                                              double* totalEnergy, double* energyParamDerivs) {
   double deltaR[ReferenceForce::LastDeltaRIndex];
   for (int i = 0; i < numParameters; i++)
       expressionSet.setVariable(bondParamIndex[i], parameters[i]);

   // ---------------------------------------------------------------------------------------

   // get deltaR, R2, and R between 2 atoms

   int atomAIndex = atomIndices[0];
   int atomBIndex = atomIndices[1];
   if (usePeriodic)
       ReferenceForce::getDeltaRPeriodic(atomCoordinates[atomAIndex], atomCoordinates[atomBIndex], boxVectors, deltaR);
   else
       ReferenceForce::getDeltaR(atomCoordinates[atomAIndex], atomCoordinates[atomBIndex], deltaR);
   
   expressionSet.setVariable(rIndex, deltaR[ReferenceForce::RIndex]);
   double dEdR            = forceExpression.evaluate();
   dEdR                   = deltaR[ReferenceForce::RIndex] > 0 ? (dEdR/deltaR[ReferenceForce::RIndex]) : 0;

   forces[atomAIndex][0] += dEdR*deltaR[ReferenceForce::XIndex];
   forces[atomAIndex][1] += dEdR*deltaR[ReferenceForce::YIndex];
   forces[atomAIndex][2] += dEdR*deltaR[ReferenceForce::ZIndex];

   forces[atomBIndex][0] -= dEdR*deltaR[ReferenceForce::XIndex];
   forces[atomBIndex][1] -= dEdR*deltaR[ReferenceForce::YIndex];
   forces[atomBIndex][2] -= dEdR*deltaR[ReferenceForce::ZIndex];

   for (int i = 0; i < energyParamDerivExpressions.size(); i++)
       energyParamDerivs[i] += energyParamDerivExpressions[i].evaluate();
   if (totalEnergy != NULL)
       *totalEnergy += energyExpression.evaluate();
}
