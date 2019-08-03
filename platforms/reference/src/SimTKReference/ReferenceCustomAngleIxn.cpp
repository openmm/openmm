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
#include "ReferenceCustomAngleIxn.h"
#include "ReferenceForce.h"

using namespace OpenMM;
using namespace std;

/**---------------------------------------------------------------------------------------

   ReferenceCustomAngleIxn constructor

   --------------------------------------------------------------------------------------- */

ReferenceCustomAngleIxn::ReferenceCustomAngleIxn(const Lepton::CompiledExpression& energyExpression,
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
        angleParamIndex.push_back(expressionSet.getVariableIndex(param));
}

/**---------------------------------------------------------------------------------------

   ReferenceCustomAngleIxn destructor

   --------------------------------------------------------------------------------------- */

ReferenceCustomAngleIxn::~ReferenceCustomAngleIxn() {
}

void ReferenceCustomAngleIxn::setPeriodic(OpenMM::Vec3* vectors) {
    usePeriodic = true;
    boxVectors[0] = vectors[0];
    boxVectors[1] = vectors[1];
    boxVectors[2] = vectors[2];
}

void ReferenceCustomAngleIxn::setGlobalParameters(std::map<std::string, double> parameters) {
    for (auto& param : parameters)
        expressionSet.setVariable(expressionSet.getVariableIndex(param.first), param.second);
}

/**---------------------------------------------------------------------------------------

   Calculate Custom Angle Ixn

   @param atomIndices      atom indices of atom participating in bond
   @param atomCoordinates  atom coordinates
   @param parameters       parameters values
   @param forces           force array (forces added to input values)
   @param totalEnergy      if not null, the energy will be added to this

   --------------------------------------------------------------------------------------- */

void ReferenceCustomAngleIxn::calculateBondIxn(vector<int>& atomIndices,
                                               vector<Vec3>& atomCoordinates,
                                               vector<double>& parameters,
                                               vector<Vec3>& forces,
                                               double* totalEnergy, double* energyParamDerivs) {
   double deltaR[2][ReferenceForce::LastDeltaRIndex];
   for (int i = 0; i < numParameters; i++)
       expressionSet.setVariable(angleParamIndex[i], parameters[i]);

   // ---------------------------------------------------------------------------------------

   // Compute the angle between the three atoms

   int atomAIndex = atomIndices[0];
   int atomBIndex = atomIndices[1];
   int atomCIndex = atomIndices[2];
   if (usePeriodic) {
      ReferenceForce::getDeltaRPeriodic(atomCoordinates[atomAIndex], atomCoordinates[atomBIndex], boxVectors, deltaR[0]);
      ReferenceForce::getDeltaRPeriodic(atomCoordinates[atomCIndex], atomCoordinates[atomBIndex], boxVectors, deltaR[1]);
   }
   else {
      ReferenceForce::getDeltaR(atomCoordinates[atomAIndex], atomCoordinates[atomBIndex], deltaR[0]);
      ReferenceForce::getDeltaR(atomCoordinates[atomCIndex], atomCoordinates[atomBIndex], deltaR[1]);
   }
   double pVector[3];
   SimTKOpenMMUtilities::crossProductVector3(deltaR[0], deltaR[1], pVector);
   double rp = sqrt(DOT3(pVector, pVector));
   if (rp < 1.0e-06)
      rp = 1.0e-06;
   double dot = DOT3(deltaR[0], deltaR[1]);
   double cosine = dot/sqrt((deltaR[0][ReferenceForce::R2Index]*deltaR[1][ReferenceForce::R2Index]));
   double angle;
   if (cosine >= 1.0)
      angle = 0.0;
   else if (cosine <= -1.0)
      angle = PI_M;
   else
      angle = acos(cosine);
   expressionSet.setVariable(thetaIndex, angle);

   // Compute the force and energy, and apply them to the atoms.
   
   double energy = energyExpression.evaluate();
   double dEdR = forceExpression.evaluate();
   double termA =  dEdR/(deltaR[0][ReferenceForce::R2Index]*rp);
   double termC = -dEdR/(deltaR[1][ReferenceForce::R2Index]*rp);

   double deltaCrossP[3][3];
   SimTKOpenMMUtilities::crossProductVector3(deltaR[0], pVector, deltaCrossP[0]);
   SimTKOpenMMUtilities::crossProductVector3(deltaR[1], pVector, deltaCrossP[2]);

   for (int ii = 0; ii < 3; ii++) {
      deltaCrossP[0][ii] *= termA;
      deltaCrossP[2][ii] *= termC;
      deltaCrossP[1][ii]  = -(deltaCrossP[0][ii]+deltaCrossP[2][ii]);
   }

   // accumulate forces

   for (int jj = 0; jj < 3; jj++) {
      for (int ii = 0; ii < 3; ii++) {
         forces[atomIndices[jj]][ii] += deltaCrossP[jj][ii];
      }
   }

   // Record parameter derivatives.

   for (int i = 0; i < energyParamDerivExpressions.size(); i++)
       energyParamDerivs[i] += energyParamDerivExpressions[i].evaluate();
   
   // accumulate energies

   if (totalEnergy != NULL)
       *totalEnergy += energy;
}

