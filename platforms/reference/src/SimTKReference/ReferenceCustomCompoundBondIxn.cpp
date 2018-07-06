
/* Portions copyright (c) 2009-2016 Stanford University and Simbios.
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

#include "SimTKOpenMMUtilities.h"
#include "ReferenceForce.h"
#include "ReferenceCustomCompoundBondIxn.h"

using std::map;
using std::pair;
using std::string;
using std::stringstream;
using std::vector;
using namespace OpenMM;

/**---------------------------------------------------------------------------------------

   ReferenceCustomCompoundBondIxn constructor

   --------------------------------------------------------------------------------------- */

ReferenceCustomCompoundBondIxn::ReferenceCustomCompoundBondIxn(int numParticlesPerBond, const vector<vector<int> >& bondAtoms,
            const Lepton::ParsedExpression& energyExpression, const vector<string>& bondParameterNames,
            const map<string, vector<int> >& distances, const map<string, vector<int> >& angles, const map<string, vector<int> >& dihedrals,
            const std::vector<Lepton::CompiledExpression> energyParamDerivExpressions) :
            bondAtoms(bondAtoms), energyExpression(energyExpression.createCompiledExpression()), usePeriodic(false),
            energyParamDerivExpressions(energyParamDerivExpressions) {
    expressionSet.registerExpression(this->energyExpression);
    for (int i = 0; i < this->energyParamDerivExpressions.size(); i++)
        expressionSet.registerExpression(this->energyParamDerivExpressions[i]);
    for (int i = 0; i < numParticlesPerBond; i++) {
        stringstream xname, yname, zname;
        xname << 'x' << (i+1);
        yname << 'y' << (i+1);
        zname << 'z' << (i+1);
        particleTerms.push_back(ReferenceCustomCompoundBondIxn::ParticleTermInfo(xname.str(), i, 0, energyExpression.differentiate(xname.str()).createCompiledExpression()));
        particleTerms.push_back(ReferenceCustomCompoundBondIxn::ParticleTermInfo(yname.str(), i, 1, energyExpression.differentiate(yname.str()).createCompiledExpression()));
        particleTerms.push_back(ReferenceCustomCompoundBondIxn::ParticleTermInfo(zname.str(), i, 2, energyExpression.differentiate(zname.str()).createCompiledExpression()));
    }
    for (auto& term : distances)
        distanceTerms.push_back(ReferenceCustomCompoundBondIxn::DistanceTermInfo(term.first, term.second, energyExpression.differentiate(term.first).createCompiledExpression()));
    for (auto& term : angles)
        angleTerms.push_back(ReferenceCustomCompoundBondIxn::AngleTermInfo(term.first, term.second, energyExpression.differentiate(term.first).createCompiledExpression()));
    for (auto& term : dihedrals)
        dihedralTerms.push_back(ReferenceCustomCompoundBondIxn::DihedralTermInfo(term.first, term.second, energyExpression.differentiate(term.first).createCompiledExpression()));
    for (int i = 0; i < particleTerms.size(); i++) {
        expressionSet.registerExpression(particleTerms[i].forceExpression);
        particleTerms[i].index = expressionSet.getVariableIndex(particleTerms[i].name);
    }
    for (int i = 0; i < distanceTerms.size(); i++) {
        expressionSet.registerExpression(distanceTerms[i].forceExpression);
        distanceTerms[i].index = expressionSet.getVariableIndex(distanceTerms[i].name);
    }
    for (int i = 0; i < angleTerms.size(); i++) {
        expressionSet.registerExpression(angleTerms[i].forceExpression);
        angleTerms[i].index = expressionSet.getVariableIndex(angleTerms[i].name);
    }
    for (int i = 0; i < dihedralTerms.size(); i++) {
        expressionSet.registerExpression(dihedralTerms[i].forceExpression);
        dihedralTerms[i].index = expressionSet.getVariableIndex(dihedralTerms[i].name);
    }
    numParameters = bondParameterNames.size();
    for (int i = 0; i < numParameters; i++)
        bondParamIndex.push_back(expressionSet.getVariableIndex(bondParameterNames[i]));
}

/**---------------------------------------------------------------------------------------

   ReferenceCustomCompoundBondIxn destructor

   --------------------------------------------------------------------------------------- */

ReferenceCustomCompoundBondIxn::~ReferenceCustomCompoundBondIxn() {
}

void ReferenceCustomCompoundBondIxn::setPeriodic(OpenMM::Vec3* vectors) {
    usePeriodic = true;
    boxVectors[0] = vectors[0];
    boxVectors[1] = vectors[1];
    boxVectors[2] = vectors[2];
}

/**---------------------------------------------------------------------------------------

   Calculate custom hbond interaction

   @param atomCoordinates    atom coordinates
   @param bondParameters     bond parameters values       bondParameters[bondIndex][parameterIndex]
   @param globalParameters   the values of global parameters
   @param forces             force array (forces added)
   @param totalEnergy        total energy

   --------------------------------------------------------------------------------------- */

void ReferenceCustomCompoundBondIxn::calculatePairIxn(vector<Vec3>& atomCoordinates, vector<vector<double> >& bondParameters,
                                             const map<string, double>& globalParameters, vector<Vec3>& forces,
                                             double* totalEnergy, double* energyParamDerivs) {
    for (auto& param : globalParameters)
        expressionSet.setVariable(expressionSet.getVariableIndex(param.first), param.second);
    int numBonds = bondAtoms.size();
    for (int bond = 0; bond < numBonds; bond++) {
        for (int i = 0; i < numParameters; i++)
            expressionSet.setVariable(bondParamIndex[i], bondParameters[bond][i]);
        calculateOneIxn(bond, atomCoordinates, forces, totalEnergy, energyParamDerivs);
    }
}

  /**---------------------------------------------------------------------------------------

     Calculate interaction for one bond

     @param bond             the index of the bond
     @param atomCoordinates  atom coordinates
     @param forces           force array (forces added)
     @param energyByAtom     atom energy
     @param totalEnergy      total energy

     --------------------------------------------------------------------------------------- */

void ReferenceCustomCompoundBondIxn::calculateOneIxn(int bond, vector<Vec3>& atomCoordinates,
                        vector<Vec3>& forces, double* totalEnergy, double* energyParamDerivs) {
    // Compute all of the variables the energy can depend on.

    const vector<int>& atoms = bondAtoms[bond];
    for (auto& term : particleTerms)
        expressionSet.setVariable(term.index, atomCoordinates[atoms[term.atom]][term.component]);
    for (auto& term : distanceTerms) {
        computeDelta(atoms[term.p1], atoms[term.p2], term.delta, atomCoordinates);
        expressionSet.setVariable(term.index, term.delta[ReferenceForce::RIndex]);
    }
    for (auto& term : angleTerms) {
        computeDelta(atoms[term.p1], atoms[term.p2], term.delta1, atomCoordinates);
        computeDelta(atoms[term.p3], atoms[term.p2], term.delta2, atomCoordinates);
        expressionSet.setVariable(term.index, computeAngle(term.delta1, term.delta2));
    }
    for (auto& term : dihedralTerms) {
        computeDelta(atoms[term.p2], atoms[term.p1], term.delta1, atomCoordinates);
        computeDelta(atoms[term.p2], atoms[term.p3], term.delta2, atomCoordinates);
        computeDelta(atoms[term.p4], atoms[term.p3], term.delta3, atomCoordinates);
        double dotDihedral, signOfDihedral;
        double* crossProduct[] = {term.cross1, term.cross2};
        expressionSet.setVariable(term.index,getDihedralAngleBetweenThreeVectors(term.delta1, term.delta2, term.delta3, crossProduct, &dotDihedral, term.delta1, &signOfDihedral, 1));
    }
    
    // Apply forces based on individual particle coordinates.
    
    for (auto& term : particleTerms)
        forces[atoms[term.atom]][term.component] -= term.forceExpression.evaluate();

    // Apply forces based on distances.

    for (auto& term : distanceTerms) {
        double dEdR = term.forceExpression.evaluate()/(term.delta[ReferenceForce::RIndex]);
        for (int i = 0; i < 3; i++) {
           double force  = -dEdR*term.delta[i];
           forces[atoms[term.p1]][i] -= force;
           forces[atoms[term.p2]][i] += force;
        }
    }

    // Apply forces based on angles.

    for (auto& term : angleTerms) {
        double dEdTheta = term.forceExpression.evaluate();
        double thetaCross[ReferenceForce::LastDeltaRIndex];
        SimTKOpenMMUtilities::crossProductVector3(term.delta1, term.delta2, thetaCross);
        double lengthThetaCross = sqrt(DOT3(thetaCross, thetaCross));
        if (lengthThetaCross < 1.0e-06)
            lengthThetaCross = 1.0e-06;
        double termA = dEdTheta/(term.delta1[ReferenceForce::R2Index]*lengthThetaCross);
        double termC = -dEdTheta/(term.delta2[ReferenceForce::R2Index]*lengthThetaCross);
        double deltaCrossP[3][3];
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

    for (auto& term : dihedralTerms) {
        double dEdTheta = term.forceExpression.evaluate();
        double internalF[4][3];
        double forceFactors[4];
        double normCross1 = DOT3(term.cross1, term.cross1);
        double normBC = term.delta2[ReferenceForce::RIndex];
        forceFactors[0] = (-dEdTheta*normBC)/normCross1;
        double normCross2 = DOT3(term.cross2, term.cross2);
        forceFactors[3] = (dEdTheta*normBC)/normCross2;
        forceFactors[1] = DOT3(term.delta1, term.delta2);
        forceFactors[1] /= term.delta2[ReferenceForce::R2Index];
        forceFactors[2] = DOT3(term.delta3, term.delta2);
        forceFactors[2] /= term.delta2[ReferenceForce::R2Index];
        for (int i = 0; i < 3; i++) {
            internalF[0][i] = forceFactors[0]*term.cross1[i];
            internalF[3][i] = forceFactors[3]*term.cross2[i];
            double s = forceFactors[1]*internalF[0][i] - forceFactors[2]*internalF[3][i];
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
        *totalEnergy += energyExpression.evaluate();
    
    // Compute derivatives of the energy.
    
    for (int i = 0; i < energyParamDerivExpressions.size(); i++)
        energyParamDerivs[i] += energyParamDerivExpressions[i].evaluate();
}

void ReferenceCustomCompoundBondIxn::computeDelta(int atom1, int atom2, double* delta, vector<Vec3>& atomCoordinates) const {
    if (usePeriodic)
        ReferenceForce::getDeltaRPeriodic(atomCoordinates[atom1], atomCoordinates[atom2], boxVectors, delta);
    else
        ReferenceForce::getDeltaR(atomCoordinates[atom1], atomCoordinates[atom2], delta);
}

double ReferenceCustomCompoundBondIxn::computeAngle(double* vec1, double* vec2) {
    double dot = DOT3(vec1, vec2);
    double cosine = dot/sqrt((vec1[ReferenceForce::R2Index]*vec2[ReferenceForce::R2Index]));
    double angle;
    if (cosine >= 1)
        angle = 0;
    else if (cosine <= -1)
        angle = PI_M;
    else
        angle = acos(cosine);
    return angle;
}
