
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

#include "SimTKOpenMMCommon.h"
#include "SimTKOpenMMLog.h"
#include "SimTKOpenMMUtilities.h"
#include "ReferenceForce.h"
#include "ReferenceCustomHbondIxn.h"

using std::map;
using std::pair;
using std::set;
using std::string;
using std::stringstream;
using std::vector;
using OpenMM::RealVec;

/**---------------------------------------------------------------------------------------

   ReferenceCustomHbondIxn constructor

   --------------------------------------------------------------------------------------- */

ReferenceCustomHbondIxn::ReferenceCustomHbondIxn(const vector<vector<int> >& donorAtoms, const vector<vector<int> >& acceptorAtoms,
            const Lepton::ParsedExpression& energyExpression, const vector<string>& donorParameterNames, const vector<string>& acceptorParameterNames,
            const map<string, vector<int> >& distances, const map<string, vector<int> >& angles, const map<string, vector<int> >& dihedrals) :
            cutoff(false), periodic(false), donorAtoms(donorAtoms), acceptorAtoms(acceptorAtoms), energyExpression(energyExpression.createProgram()),
            donorParamNames(donorParameterNames), acceptorParamNames(acceptorParameterNames) {
    for (map<string, vector<int> >::const_iterator iter = distances.begin(); iter != distances.end(); ++iter)
        distanceTerms.push_back(ReferenceCustomHbondIxn::DistanceTermInfo(iter->first, iter->second, energyExpression.differentiate(iter->first).optimize().createProgram()));
    for (map<string, vector<int> >::const_iterator iter = angles.begin(); iter != angles.end(); ++iter)
        angleTerms.push_back(ReferenceCustomHbondIxn::AngleTermInfo(iter->first, iter->second, energyExpression.differentiate(iter->first).optimize().createProgram()));
    for (map<string, vector<int> >::const_iterator iter = dihedrals.begin(); iter != dihedrals.end(); ++iter)
        dihedralTerms.push_back(ReferenceCustomHbondIxn::DihedralTermInfo(iter->first, iter->second, energyExpression.differentiate(iter->first).optimize().createProgram()));
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

void ReferenceCustomHbondIxn::setPeriodic(RealVec& boxSize) {
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
   @param exclusions         exclusion indices
                             exclusions[donorIndex] contains the list of excluded acceptors for that donor
   @param globalParameters   the values of global parameters
   @param forces             force array (forces added)
   @param totalEnergy        total energy

   --------------------------------------------------------------------------------------- */

void ReferenceCustomHbondIxn::calculatePairIxn(vector<RealVec>& atomCoordinates, RealOpenMM** donorParameters, RealOpenMM** acceptorParameters,
                                             vector<set<int> >& exclusions, const map<string, double>& globalParameters, vector<RealVec>& forces,
                                             RealOpenMM* totalEnergy) const {

   map<string, double> variables = globalParameters;

   // allocate and initialize exclusion array

   int numDonors = donorAtoms.size();
   int numAcceptors = acceptorAtoms.size();

   for( int donor = 0; donor < numDonors; donor++ ){
      // Initialize per-donor parameters.

      for (int j = 0; j < (int) donorParamNames.size(); j++)
          variables[donorParamNames[j]] = donorParameters[donor][j];

      // loop over atom pairs

      for( int acceptor = 0; acceptor < numAcceptors; acceptor++ ){
         if (exclusions[donor].find(acceptor) == exclusions[donor].end()) {
             for (int j = 0; j < (int) acceptorParamNames.size(); j++)
                 variables[acceptorParamNames[j]] = acceptorParameters[acceptor][j];
             calculateOneIxn(donor, acceptor, atomCoordinates, variables, forces, totalEnergy);
         }
      }
   }
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

void ReferenceCustomHbondIxn::calculateOneIxn(int donor, int acceptor, vector<RealVec>& atomCoordinates,
                        map<string, double>& variables, vector<RealVec>& forces, RealOpenMM* totalEnergy) const {

    // ---------------------------------------------------------------------------------------

    static const std::string methodName = "\nReferenceCustomHbondIxn::calculateOneIxn";

    // ---------------------------------------------------------------------------------------

    int atoms[6];
    atoms[0] = acceptorAtoms[acceptor][0];
    atoms[1] = acceptorAtoms[acceptor][1];
    atoms[2] = acceptorAtoms[acceptor][2];
    atoms[3] = donorAtoms[donor][0];
    atoms[4] = donorAtoms[donor][1];
    atoms[5] = donorAtoms[donor][2];

    // Compute the distance between the primary donor and acceptor atoms, and compare to the cutoff.

    if (cutoff) {
        RealOpenMM delta[ReferenceForce::LastDeltaRIndex];
        computeDelta(atoms[0], atoms[3], delta, atomCoordinates);
        if (delta[ReferenceForce::RIndex] >= cutoffDistance)
            return;
    }

    // Compute all of the variables the energy can depend on.

    for (int i = 0; i < (int) distanceTerms.size(); i++) {
        const DistanceTermInfo& term = distanceTerms[i];
        computeDelta(atoms[term.p1], atoms[term.p2], term.delta, atomCoordinates);
        variables[term.name] = term.delta[ReferenceForce::RIndex];
    }
    for (int i = 0; i < (int) angleTerms.size(); i++) {
        const AngleTermInfo& term = angleTerms[i];
        computeDelta(atoms[term.p1], atoms[term.p2], term.delta1, atomCoordinates);
        computeDelta(atoms[term.p3], atoms[term.p2], term.delta2, atomCoordinates);
        variables[term.name] = computeAngle(term.delta1, term.delta2);
    }
    for (int i = 0; i < (int) dihedralTerms.size(); i++) {
        const DihedralTermInfo& term = dihedralTerms[i];
        computeDelta(atoms[term.p2], atoms[term.p1], term.delta1, atomCoordinates);
        computeDelta(atoms[term.p2], atoms[term.p3], term.delta2, atomCoordinates);
        computeDelta(atoms[term.p4], atoms[term.p3], term.delta3, atomCoordinates);
        RealOpenMM dotDihedral, signOfDihedral;
        RealOpenMM* crossProduct[] = {term.cross1, term.cross2};
        variables[term.name] = getDihedralAngleBetweenThreeVectors(term.delta1, term.delta2, term.delta3, crossProduct, &dotDihedral, term.delta1, &signOfDihedral, 1);
    }

    // Apply forces based on distances.

    for (int i = 0; i < (int) distanceTerms.size(); i++) {
        const DistanceTermInfo& term = distanceTerms[i];
        RealOpenMM dEdR = (RealOpenMM) (term.forceExpression.evaluate(variables)/(term.delta[ReferenceForce::RIndex]));
        for (int i = 0; i < 3; i++) {
           RealOpenMM force  = -dEdR*term.delta[i];
           forces[atoms[term.p1]][i] -= force;
           forces[atoms[term.p2]][i] += force;
        }
    }

    // Apply forces based on angles.

    for (int i = 0; i < (int) angleTerms.size(); i++) {
        const AngleTermInfo& term = angleTerms[i];
        RealOpenMM dEdTheta = (RealOpenMM) term.forceExpression.evaluate(variables);
        RealOpenMM thetaCross[ReferenceForce::LastDeltaRIndex];
        SimTKOpenMMUtilities::crossProductVector3(term.delta1, term.delta2, thetaCross);
        RealOpenMM lengthThetaCross = SQRT(DOT3(thetaCross, thetaCross));
        if (lengthThetaCross < 1.0e-06)
            lengthThetaCross = (RealOpenMM) 1.0e-06;
        RealOpenMM termA = dEdTheta/(term.delta1[ReferenceForce::R2Index]*lengthThetaCross);
        RealOpenMM termC = -dEdTheta/(term.delta2[ReferenceForce::R2Index]*lengthThetaCross);
        RealOpenMM deltaCrossP[3][3];
        SimTKOpenMMUtilities::crossProductVector3(term.delta1, thetaCross, deltaCrossP[0]);
        SimTKOpenMMUtilities::crossProductVector3(term.delta2, thetaCross, deltaCrossP[2]);
        for (int i = 0; i < 3; i++) {
            deltaCrossP[0][i] *= termA;
            deltaCrossP[2][i] *= termC;
            deltaCrossP[1][i] = -(deltaCrossP[0][i]+deltaCrossP[2][i]);
        }
        for (int i = 0; i < 3; i++) {
            forces[atoms[term.p1]][i] += deltaCrossP[0][i];
            forces[atoms[term.p2]][i] += deltaCrossP[1][i];
            forces[atoms[term.p3]][i] += deltaCrossP[2][i];
        }
    }

    // Apply forces based on dihedrals.

    for (int i = 0; i < (int) dihedralTerms.size(); i++) {
        const DihedralTermInfo& term = dihedralTerms[i];
        RealOpenMM dEdTheta = (RealOpenMM) term.forceExpression.evaluate(variables);
        RealOpenMM internalF[4][3];
        RealOpenMM forceFactors[4];
        RealOpenMM normCross1 = DOT3(term.cross1, term.cross1);
        RealOpenMM normBC = term.delta2[ReferenceForce::RIndex];
        forceFactors[0] = (-dEdTheta*normBC)/normCross1;
        RealOpenMM normCross2 = DOT3(term.cross2, term.cross2);
                   forceFactors[3] = (dEdTheta*normBC)/normCross2;
                   forceFactors[1] = DOT3(term.delta1, term.delta2);
                   forceFactors[1] /= term.delta2[ReferenceForce::R2Index];
                   forceFactors[2] = DOT3(term.delta3, term.delta2);
                   forceFactors[2] /= term.delta2[ReferenceForce::R2Index];
        for (int i = 0; i < 3; i++) {
            internalF[0][i] = forceFactors[0]*term.cross1[i];
            internalF[3][i] = forceFactors[3]*term.cross2[i];
            RealOpenMM s = forceFactors[1]*internalF[0][i] - forceFactors[2]*internalF[3][i];
            internalF[1][i] = internalF[0][i] - s;
            internalF[2][i] = internalF[3][i] + s;
        }
        for (int i = 0; i < 3; i++) {
            forces[atoms[term.p1]][i] += internalF[0][i];
            forces[atoms[term.p2]][i] -= internalF[1][i];
            forces[atoms[term.p3]][i] -= internalF[2][i];
            forces[atoms[term.p4]][i] += internalF[3][i];
        }
    }

    // Add the energy

    if (totalEnergy)
        *totalEnergy += (RealOpenMM) energyExpression.evaluate(variables);
}

void ReferenceCustomHbondIxn::computeDelta(int atom1, int atom2, RealOpenMM* delta, vector<RealVec>& atomCoordinates) const {
    if (periodic)
        ReferenceForce::getDeltaRPeriodic(atomCoordinates[atom1], atomCoordinates[atom2], periodicBoxSize, delta);
    else
        ReferenceForce::getDeltaR(atomCoordinates[atom1], atomCoordinates[atom2], delta);
}

RealOpenMM ReferenceCustomHbondIxn::computeAngle(RealOpenMM* vec1, RealOpenMM* vec2) {
    RealOpenMM dot = DOT3(vec1, vec2);
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
