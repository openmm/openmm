
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
    for (int i = 0; i < particleTerms.size(); i++) {
        expressionSet.registerExpression(particleTerms[i].forceExpression);
        particleTerms[i].index = expressionSet.getVariableIndex(particleTerms[i].name);
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
    
    // Apply forces based on particle coordinates.
    
    for (auto& term : particleTerms)
        forces[atoms[term.atom]][term.component] -= term.forceExpression.evaluate();

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
