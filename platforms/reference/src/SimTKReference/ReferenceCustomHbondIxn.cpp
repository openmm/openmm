
/* Portions copyright (c) 2009-2010 Stanford University and Simbios.
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
#include <utility>

#include "../SimTKUtilities/SimTKOpenMMCommon.h"
#include "../SimTKUtilities/SimTKOpenMMLog.h"
#include "../SimTKUtilities/SimTKOpenMMUtilities.h"
#include "ReferenceForce.h"
#include "ReferenceCustomHbondIxn.h"

using std::map;
using std::pair;
using std::string;
using std::stringstream;
using std::vector;

/**---------------------------------------------------------------------------------------

   ReferenceCustomHbondIxn constructor

   --------------------------------------------------------------------------------------- */

ReferenceCustomHbondIxn::ReferenceCustomHbondIxn(const vector<pair<int, int> >& donorAtoms, const vector<pair<int, int> >& acceptorAtoms,
            const Lepton::ExpressionProgram& energyExpression,
            const Lepton::ExpressionProgram& rForceExpression, const Lepton::ExpressionProgram& thetaForceExpression,
            const Lepton::ExpressionProgram& psiForceExpression, const Lepton::ExpressionProgram& chiForceExpression,
            const vector<string>& donorParameterNames, const vector<string>& acceptorParameterNames) :
            cutoff(false), periodic(false), donorAtoms(donorAtoms), acceptorAtoms(acceptorAtoms), energyExpression(energyExpression),
            rForceExpression(rForceExpression), thetaForceExpression(thetaForceExpression), psiForceExpression(psiForceExpression),
            chiForceExpression(chiForceExpression), donorParamNames(donorParameterNames), acceptorParamNames(acceptorParameterNames) {
}

/**---------------------------------------------------------------------------------------

   ReferenceCustomHbondIxn destructor

   --------------------------------------------------------------------------------------- */

ReferenceCustomHbondIxn::~ReferenceCustomHbondIxn( ){
}

  /**---------------------------------------------------------------------------------------

     Set the force to use a cutoff.

     @param distance            the cutoff distance

     --------------------------------------------------------------------------------------- */

void ReferenceCustomHbondIxn::setUseCutoff(RealOpenMM distance) {
    cutoff = true;
    cutoffDistance = distance;
}

  /**---------------------------------------------------------------------------------------

     Set the force to use periodic boundary conditions.  This requires that a cutoff has
     also been set, and the smallest side of the periodic box is at least twice the cutoff
     distance.

     @param boxSize             the X, Y, and Z widths of the periodic box

     --------------------------------------------------------------------------------------- */

void ReferenceCustomHbondIxn::setPeriodic(RealOpenMM* boxSize) {
    assert(cutoff);
    assert(boxSize[0] >= 2.0*cutoffDistance);
    assert(boxSize[1] >= 2.0*cutoffDistance);
    assert(boxSize[2] >= 2.0*cutoffDistance);
    periodic = true;
    periodicBoxSize[0] = boxSize[0];
    periodicBoxSize[1] = boxSize[1];
    periodicBoxSize[2] = boxSize[2];
  }


/**---------------------------------------------------------------------------------------

   Calculate custom hbond interaction

   @param atomCoordinates    atom coordinates
   @param donorParameters    donor parameters values       donorParameters[donorIndex][parameterIndex]
   @param acceptorParameters acceptor parameters values    acceptorParameters[acceptorIndex][parameterIndex]
   @param exclusions         exclusion indices             exclusions[donorIndex][acceptorToExcludeIndex]
                             exclusions[donorIndex][0] = number of exclusions
                             exclusions[donorIndex][no.-1] = indices of acceptors to excluded from
                             interacting w/ donor donorIndex
   @param globalParameters   the values of global parameters
   @param forces             force array (forces added)
   @param totalEnergy        total energy

   --------------------------------------------------------------------------------------- */

void ReferenceCustomHbondIxn::calculatePairIxn(RealOpenMM** atomCoordinates, RealOpenMM** donorParameters, RealOpenMM** acceptorParameters,
                                             int** exclusions, const map<string, double>& globalParameters, RealOpenMM** forces,
                                             RealOpenMM* totalEnergy) const {

   map<string, double> variables = globalParameters;

   // allocate and initialize exclusion array

   int numDonors = donorAtoms.size();
   int numAcceptors = acceptorAtoms.size();
   int* exclusionIndices = new int[numAcceptors];
   for( int ii = 0; ii < numAcceptors; ii++ ){
      exclusionIndices[ii] = -1;
   }

   for( int donor = 0; donor < numDonors; donor++ ){

      // set exclusions

      for (int j = 1; j <= exclusions[donor][0]; j++)
         exclusionIndices[exclusions[donor][j]] = donor;

      // Initialize per-donor parameters.

      for (int j = 0; j < (int) donorParamNames.size(); j++)
          variables[donorParamNames[j]] = donorParameters[donor][j];

      // loop over atom pairs

      for( int acceptor = 0; acceptor < numAcceptors; acceptor++ ){

         if( exclusionIndices[acceptor] != donor ){
             for (int j = 0; j < (int) acceptorParamNames.size(); j++)
                 variables[acceptorParamNames[j]] = acceptorParameters[acceptor][j];
             calculateOneIxn(donor, acceptor, atomCoordinates, variables, forces, totalEnergy);
         }
      }
   }

   delete[] exclusionIndices;
}

  /**---------------------------------------------------------------------------------------

     Calculate custom interaction between a donor and an acceptor

     @param donor            the index of the donor
     @param acceptor         the index of the acceptor
     @param atomCoordinates  atom coordinates
     @param variables        the values of variables that may appear in expressions
     @param forces           force array (forces added)
     @param energyByAtom     atom energy
     @param totalEnergy      total energy

     --------------------------------------------------------------------------------------- */

void ReferenceCustomHbondIxn::calculateOneIxn(int donor, int acceptor, RealOpenMM** atomCoordinates,
                        map<string, double>& variables, RealOpenMM** forces, RealOpenMM* totalEnergy) const {

    // ---------------------------------------------------------------------------------------

    static const std::string methodName = "\nReferenceCustomHbondIxn::calculateOneIxn";

    // ---------------------------------------------------------------------------------------

    // Visual Studio complains if crossProduct declared as 'crossProduct[2][3]'

    RealOpenMM crossProductMemory[6];
    RealOpenMM* crossProduct[2];
    crossProduct[0] = crossProductMemory;
    crossProduct[1] = crossProductMemory + 3;

    // Compute the distance between the primary donor and acceptor atoms, and compare to the cutoff.

    int d1 = donorAtoms[donor].first;
    int a1 = acceptorAtoms[acceptor].first;
    RealOpenMM deltaD1A1[ReferenceForce::LastDeltaRIndex];
    if (periodic)
        ReferenceForce::getDeltaRPeriodic(atomCoordinates[a1], atomCoordinates[d1], periodicBoxSize, deltaD1A1);
    else
        ReferenceForce::getDeltaR(atomCoordinates[a1], atomCoordinates[d1], deltaD1A1);
    if (cutoff && deltaD1A1[ReferenceForce::RIndex] >= cutoffDistance)
        return;

    // Compute all of the variables the energy can depend on.
    
    int d2 = donorAtoms[donor].second;
    int a2 = acceptorAtoms[acceptor].second;
    RealOpenMM deltaD1D2[ReferenceForce::LastDeltaRIndex];
    RealOpenMM deltaA2A1[ReferenceForce::LastDeltaRIndex];
    if (periodic) {
        ReferenceForce::getDeltaRPeriodic(atomCoordinates[d2], atomCoordinates[d1], periodicBoxSize, deltaD1D2);
        ReferenceForce::getDeltaRPeriodic(atomCoordinates[a1], atomCoordinates[a2], periodicBoxSize, deltaA2A1);
    }
    else {
        ReferenceForce::getDeltaR(atomCoordinates[d2], atomCoordinates[d1], deltaD1D2);
        ReferenceForce::getDeltaR(atomCoordinates[a1], atomCoordinates[a2], deltaA2A1);
    }
    variables["r"] = deltaD1A1[ReferenceForce::RIndex];
    variables["theta"] = computeAngle(deltaD1A1, deltaD1D2, 1);
    variables["psi"] = computeAngle(deltaD1A1, deltaA2A1, 1);
    RealOpenMM dotDihedral, signOfDihedral;
    variables["chi"] =  getDihedralAngleBetweenThreeVectors(deltaA2A1, deltaD1A1, deltaD1D2,
                                  crossProduct, &dotDihedral, deltaA2A1, &signOfDihedral, 1);

    // Apply forces based on r.

    RealOpenMM dEdR = (RealOpenMM) (rForceExpression.evaluate(variables)/(deltaD1A1[ReferenceForce::RIndex]));
    if (dEdR != 0) {
        for (int i = 0; i < 3; i++) {
           RealOpenMM force  = -dEdR*deltaD1A1[i];
           forces[d1][i] += force;
           forces[a1][i] -= force;
        }
    }

    // Apply forces based on theta.

    RealOpenMM dEdTheta = (RealOpenMM) thetaForceExpression.evaluate(variables);
    if (dEdTheta != 0) {
        RealOpenMM thetaCross[ReferenceForce::LastDeltaRIndex];
        SimTKOpenMMUtilities::crossProductVector3(deltaD1D2, deltaD1A1, thetaCross);
        RealOpenMM lengthThetaCross = SQRT(DOT3(thetaCross, thetaCross));
        if (lengthThetaCross < 1.0e-06)
            lengthThetaCross = (RealOpenMM) 1.0e-06;
        RealOpenMM termA = dEdTheta/(deltaD1D2[ReferenceForce::R2Index]*lengthThetaCross);
        RealOpenMM termC = -dEdTheta/(deltaD1A1[ReferenceForce::R2Index]*lengthThetaCross);
        RealOpenMM deltaCrossP[3][3];
        SimTKOpenMMUtilities::crossProductVector3(deltaD1D2, thetaCross, deltaCrossP[0]);
        SimTKOpenMMUtilities::crossProductVector3(deltaD1A1, thetaCross, deltaCrossP[2]);
        for (int i = 0; i < 3; i++) {
            deltaCrossP[0][i] *= termA;
            deltaCrossP[2][i] *= termC;
            deltaCrossP[1][i] = -(deltaCrossP[0][i]+deltaCrossP[2][i]);
        }
        for (int i = 0; i < 3; i++) {
            forces[d2][i] += deltaCrossP[0][i];
            forces[d1][i] += deltaCrossP[1][i];
            forces[a1][i] += deltaCrossP[2][i];
        }
    }

    // Apply forces based on psi.

    RealOpenMM dEdPsi = (RealOpenMM) psiForceExpression.evaluate(variables);
    if (dEdPsi != 0) {
        RealOpenMM psiCross[ReferenceForce::LastDeltaRIndex];
        SimTKOpenMMUtilities::crossProductVector3(deltaA2A1, deltaD1A1, psiCross);
        RealOpenMM lengthPsiCross = SQRT(DOT3(psiCross, psiCross));
        if (lengthPsiCross < 1.0e-06)
            lengthPsiCross = (RealOpenMM) 1.0e-06;
        RealOpenMM termA =  dEdPsi/(deltaD1A1[ReferenceForce::R2Index]*lengthPsiCross);
        RealOpenMM termC = -dEdPsi/(deltaA2A1[ReferenceForce::R2Index]*lengthPsiCross);
        RealOpenMM deltaCrossP[3][3];
        SimTKOpenMMUtilities::crossProductVector3(deltaD1A1, psiCross, deltaCrossP[0]);
        SimTKOpenMMUtilities::crossProductVector3(deltaA2A1, psiCross, deltaCrossP[2]);
        for (int i = 0; i < 3; i++) {
            deltaCrossP[0][i] *= termA;
            deltaCrossP[2][i] *= termC;
            deltaCrossP[1][i] = -(deltaCrossP[0][i]+deltaCrossP[2][i]);
        }
        for (int i = 0; i < 3; i++) {
            forces[d1][i] += deltaCrossP[0][i];
            forces[a1][i] += deltaCrossP[1][i];
            forces[a2][i] += deltaCrossP[2][i];
        }
    }

    // Apply forces based on chi.

    RealOpenMM dEdChi = (RealOpenMM) chiForceExpression.evaluate(variables);
    if (dEdChi != 0) {
        RealOpenMM internalF[4][3];
        RealOpenMM forceFactors[4];
        RealOpenMM normCross1 = DOT3(crossProduct[0], crossProduct[0]);
        RealOpenMM normBC = deltaD1A1[ReferenceForce::RIndex];
        forceFactors[0] = (-dEdChi*normBC)/normCross1;
        RealOpenMM normCross2 = DOT3(crossProduct[1], crossProduct[1]);
                   forceFactors[3] = (dEdChi*normBC)/normCross2;
                   forceFactors[1] = DOT3(deltaA2A1, deltaD1A1);
                   forceFactors[1] /= deltaD1A1[ReferenceForce::R2Index];
                   forceFactors[2] = DOT3(deltaD1D2, deltaD1A1);
                   forceFactors[2] /= deltaD1A1[ReferenceForce::R2Index];
        for (int i = 0; i < 3; i++) {
            internalF[0][i] = forceFactors[0]*crossProduct[0][i];
            internalF[3][i] = forceFactors[3]*crossProduct[1][i];
            RealOpenMM s = forceFactors[1]*internalF[0][i] - forceFactors[2]*internalF[3][i];
            internalF[1][i] = internalF[0][i] - s;
            internalF[2][i] = internalF[3][i] + s;
        }
        for (int i = 0; i < 3; i++) {
            forces[a2][i] += internalF[0][i];
            forces[a1][i] -= internalF[1][i];
            forces[d1][i] -= internalF[2][i];
            forces[d2][i] += internalF[3][i];
        }
    }

    // Add the energy

    if (totalEnergy)
        *totalEnergy += (RealOpenMM) energyExpression.evaluate(variables);
}

RealOpenMM ReferenceCustomHbondIxn::computeAngle(RealOpenMM* vec1, RealOpenMM* vec2, RealOpenMM sign) {
    RealOpenMM dot = sign*DOT3(vec1, vec2);
    RealOpenMM cosine = dot/SQRT((vec1[ReferenceForce::R2Index]*vec2[ReferenceForce::R2Index]));
    RealOpenMM angle;
    if (cosine >= 1)
        angle = 0;
    else if (cosine <= -1)
        angle = PI_M;
    else
        angle = ACOS(cosine);
    return angle;
}
