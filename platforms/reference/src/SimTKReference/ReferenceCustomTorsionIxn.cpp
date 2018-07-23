/* Portions copyright (c) 2010-2018 Stanford University and Simbios.
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
#include "ReferenceCustomTorsionIxn.h"
#include "ReferenceForce.h"

using namespace std;
using namespace OpenMM;

/**---------------------------------------------------------------------------------------

   ReferenceCustomTorsionIxn constructor

   --------------------------------------------------------------------------------------- */

ReferenceCustomTorsionIxn::ReferenceCustomTorsionIxn(const Lepton::CompiledExpression& energyExpression,
        const Lepton::CompiledExpression& forceExpression, const vector<string>& parameterNames,
        const vector<Lepton::CompiledExpression> energyParamDerivExpressions) :
        energyExpression(energyExpression), forceExpression(forceExpression), usePeriodic(false), energyParamDerivExpressions(energyParamDerivExpressions) {
    expressionSet.registerExpression(this->energyExpression);
    expressionSet.registerExpression(this->forceExpression);
    for (int i = 0; i < this->energyParamDerivExpressions.size(); i++)
        expressionSet.registerExpression(this->energyParamDerivExpressions[i]);
    thetaIndex = expressionSet.getVariableIndex("theta");
    numParameters = parameterNames.size();
    for (auto& param : parameterNames)
        torsionParamIndex.push_back(expressionSet.getVariableIndex(param));
}

/**---------------------------------------------------------------------------------------

   ReferenceCustomTorsionIxn destructor

   --------------------------------------------------------------------------------------- */

ReferenceCustomTorsionIxn::~ReferenceCustomTorsionIxn() {
}

void ReferenceCustomTorsionIxn::setPeriodic(OpenMM::Vec3* vectors) {
    usePeriodic = true;
    boxVectors[0] = vectors[0];
    boxVectors[1] = vectors[1];
    boxVectors[2] = vectors[2];
}

void ReferenceCustomTorsionIxn::setGlobalParameters(std::map<std::string, double> parameters) {
    for (auto& param : parameters)
        expressionSet.setVariable(expressionSet.getVariableIndex(param.first), param.second);
}

/**---------------------------------------------------------------------------------------

   Calculate Custom Torsion Ixn

   @param atomIndices      atom indices of atom participating in bond
   @param atomCoordinates  atom coordinates
   @param parameters       parameters values
   @param forces           force array (forces added to input values)
   @param totalEnergy      if not null, the energy will be added to this

   --------------------------------------------------------------------------------------- */

void ReferenceCustomTorsionIxn::calculateBondIxn(vector<int>& atomIndices,
                                                vector<Vec3>& atomCoordinates,
                                                vector<double>& parameters,
                                                vector<Vec3>& forces,
                                                double* totalEnergy, double* energyParamDerivs) {
   double deltaR[3][ReferenceForce::LastDeltaRIndex];
   for (int i = 0; i < numParameters; i++)
       expressionSet.setVariable(torsionParamIndex[i], parameters[i]);

   // ---------------------------------------------------------------------------------------

   // get deltaR, R2, and R between three pairs of atoms: [j,i], [j,k], [l,k]

   int atomAIndex = atomIndices[0];
   int atomBIndex = atomIndices[1];
   int atomCIndex = atomIndices[2];
   int atomDIndex = atomIndices[3];
   if (usePeriodic) {
      ReferenceForce::getDeltaRPeriodic(atomCoordinates[atomBIndex], atomCoordinates[atomAIndex], boxVectors, deltaR[0]);  
      ReferenceForce::getDeltaRPeriodic(atomCoordinates[atomBIndex], atomCoordinates[atomCIndex], boxVectors, deltaR[1]);  
      ReferenceForce::getDeltaRPeriodic(atomCoordinates[atomDIndex], atomCoordinates[atomCIndex], boxVectors, deltaR[2]);  
   }
   else {
      ReferenceForce::getDeltaR(atomCoordinates[atomBIndex], atomCoordinates[atomAIndex], deltaR[0]);  
      ReferenceForce::getDeltaR(atomCoordinates[atomBIndex], atomCoordinates[atomCIndex], deltaR[1]);  
      ReferenceForce::getDeltaR(atomCoordinates[atomDIndex], atomCoordinates[atomCIndex], deltaR[2]);  
   }

   // Visual Studio complains if crossProduct declared as 'crossProduct[2][3]'

   double crossProductMemory[6];
   double* crossProduct[2];
   crossProduct[0] = crossProductMemory;
   crossProduct[1] = crossProductMemory + 3;

   // get dihedral angle

   double dotDihedral;
   double signOfAngle;
   double angle = getDihedralAngleBetweenThreeVectors(deltaR[0], deltaR[1], deltaR[2], crossProduct, &dotDihedral, deltaR[0], &signOfAngle, 1);
   expressionSet.setVariable(thetaIndex, angle);

   // evaluate delta angle, dE/d(angle)

   double dEdAngle = forceExpression.evaluate();

   // compute force

   double internalF[4][3];
   double forceFactors[4];
   double normCross1         = DOT3(crossProduct[0], crossProduct[0]);
   double normBC             = deltaR[1][ReferenceForce::RIndex];
          forceFactors[0]    = (-dEdAngle*normBC)/normCross1;

   double normCross2         = DOT3(crossProduct[1], crossProduct[1]);
          forceFactors[3]    = (dEdAngle*normBC)/normCross2;

          forceFactors[1]    = DOT3(deltaR[0], deltaR[1]);
          forceFactors[1]   /= deltaR[1][ReferenceForce::R2Index];

          forceFactors[2]    = DOT3(deltaR[2], deltaR[1]);
          forceFactors[2]   /= deltaR[1][ReferenceForce::R2Index];

   for (int ii = 0; ii < 3; ii++) {

      internalF[0][ii]  = forceFactors[0]*crossProduct[0][ii];
      internalF[3][ii]  = forceFactors[3]*crossProduct[1][ii];

      double s          = forceFactors[1]*internalF[0][ii] - forceFactors[2]*internalF[3][ii];

      internalF[1][ii]  = internalF[0][ii] - s;
      internalF[2][ii]  = internalF[3][ii] + s;
   }

   // accumulate forces

   for (int ii = 0; ii < 3; ii++) {
      forces[atomAIndex][ii] += internalF[0][ii];
      forces[atomBIndex][ii] -= internalF[1][ii];
      forces[atomCIndex][ii] -= internalF[2][ii];
      forces[atomDIndex][ii] += internalF[3][ii];
   }

   // Record parameter derivatives.

   for (int i = 0; i < energyParamDerivExpressions.size(); i++)
       energyParamDerivs[i] += energyParamDerivExpressions[i].evaluate();

   // accumulate energies

   if (totalEnergy != NULL)
       *totalEnergy += energyExpression.evaluate();
}

