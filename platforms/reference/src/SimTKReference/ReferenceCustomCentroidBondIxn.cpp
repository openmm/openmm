
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
#include <utility>

#include "SimTKOpenMMUtilities.h"
#include "ReferenceForce.h"
#include "ReferenceCustomCentroidBondIxn.h"

using std::map;
using std::pair;
using std::string;
using std::stringstream;
using std::vector;
using namespace OpenMM;

ReferenceCustomCentroidBondIxn::ReferenceCustomCentroidBondIxn(int numGroupsPerBond, const vector<vector<int> >& groupAtoms,
            const vector<vector<double> >& normalizedWeights, const vector<vector<int> >& bondGroups,
            const Lepton::ParsedExpression& energyExpression, const vector<string>& bondParameterNames,
            const std::vector<Lepton::CompiledExpression> energyParamDerivExpressions) :
            groupAtoms(groupAtoms), normalizedWeights(normalizedWeights), bondGroups(bondGroups), energyExpression(energyExpression.createCompiledExpression()),
            usePeriodic(false), energyParamDerivExpressions(energyParamDerivExpressions) {
    expressionSet.registerExpression(this->energyExpression);
    for (int i = 0; i < this->energyParamDerivExpressions.size(); i++)
        expressionSet.registerExpression(this->energyParamDerivExpressions[i]);
    for (int i = 0; i < numGroupsPerBond; i++) {
        stringstream xname, yname, zname;
        xname << 'x' << (i+1);
        yname << 'y' << (i+1);
        zname << 'z' << (i+1);
        positionTerms.push_back(ReferenceCustomCentroidBondIxn::PositionTermInfo(xname.str(), i, 0, energyExpression.differentiate(xname.str()).createCompiledExpression()));
        positionTerms.push_back(ReferenceCustomCentroidBondIxn::PositionTermInfo(yname.str(), i, 1, energyExpression.differentiate(yname.str()).createCompiledExpression()));
        positionTerms.push_back(ReferenceCustomCentroidBondIxn::PositionTermInfo(zname.str(), i, 2, energyExpression.differentiate(zname.str()).createCompiledExpression()));
    }
    for (int i = 0; i < positionTerms.size(); i++) {
        expressionSet.registerExpression(positionTerms[i].forceExpression);
        positionTerms[i].index = expressionSet.getVariableIndex(positionTerms[i].name);
    }
    numParameters = bondParameterNames.size();
    for (int i = 0; i < numParameters; i++)
        bondParamIndex.push_back(expressionSet.getVariableIndex(bondParameterNames[i]));
}

ReferenceCustomCentroidBondIxn::~ReferenceCustomCentroidBondIxn() {
}

void ReferenceCustomCentroidBondIxn::setPeriodic(OpenMM::Vec3* vectors) {
    usePeriodic = true;
    boxVectors[0] = vectors[0];
    boxVectors[1] = vectors[1];
    boxVectors[2] = vectors[2];
}

void ReferenceCustomCentroidBondIxn::calculatePairIxn(vector<Vec3>& atomCoordinates, vector<vector<double> >& bondParameters,
                                             const map<string, double>& globalParameters, vector<Vec3>& forces,
                                             double* totalEnergy, double* energyParamDerivs) {

    // First compute the center of each group.

    int numGroups = groupAtoms.size();
    vector<Vec3> groupCenters(numGroups);
    for (int group = 0; group < numGroups; group++) {
        for (int i = 0; i < groupAtoms[group].size(); i++)
            groupCenters[group] += atomCoordinates[groupAtoms[group][i]]*normalizedWeights[group][i];
    }

    // Compute the forces on groups.

    for (auto& param : globalParameters)
        expressionSet.setVariable(expressionSet.getVariableIndex(param.first), param.second);
    vector<Vec3> groupForces(numGroups);
    int numBonds = bondGroups.size();
    for (int bond = 0; bond < numBonds; bond++) {
        for (int i = 0; i < numParameters; i++)
            expressionSet.setVariable(bondParamIndex[i], bondParameters[bond][i]);
        calculateOneIxn(bond, groupCenters, groupForces, totalEnergy, energyParamDerivs);
    }

    // Apply the forces to the individual atoms.

    for (int group = 0; group < numGroups; group++) {
        for (int i = 0; i < groupAtoms[group].size(); i++)
            forces[groupAtoms[group][i]] += groupForces[group]*normalizedWeights[group][i];
    }
}

void ReferenceCustomCentroidBondIxn::calculateOneIxn(int bond, vector<Vec3>& groupCenters,
                        vector<Vec3>& forces, double* totalEnergy, double* energyParamDerivs) {
    // Compute all of the variables the energy can depend on.

    const vector<int>& groups = bondGroups[bond];
    for (auto& term : positionTerms)
        expressionSet.setVariable(term.index, groupCenters[groups[term.group]][term.component]);

    // Apply forces based on particle coordinates.

    for (auto& term : positionTerms)
        forces[groups[term.group]][term.component] -= term.forceExpression.evaluate();

    // Add the energy

    if (totalEnergy)
        *totalEnergy += energyExpression.evaluate();
    
    // Compute derivatives of the energy.
    
    for (int i = 0; i < energyParamDerivExpressions.size(); i++)
        energyParamDerivs[i] += energyParamDerivExpressions[i].evaluate();
}

void ReferenceCustomCentroidBondIxn::computeDelta(int group1, int group2, double* delta, vector<Vec3>& groupCenters) const {
    if (usePeriodic)
        ReferenceForce::getDeltaRPeriodic(groupCenters[group1], groupCenters[group2], boxVectors, delta);
    else
        ReferenceForce::getDeltaR(groupCenters[group1], groupCenters[group2], delta);
}

double ReferenceCustomCentroidBondIxn::computeAngle(double* vec1, double* vec2) {
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
