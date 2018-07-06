
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
#include "ReferenceForce.h"
#include "ReferenceCustomNonbondedIxn.h"

using std::map;
using std::pair;
using std::string;
using std::stringstream;
using std::set;
using std::vector;
using namespace OpenMM;

/**---------------------------------------------------------------------------------------

   ReferenceCustomNonbondedIxn constructor

   --------------------------------------------------------------------------------------- */

ReferenceCustomNonbondedIxn::ReferenceCustomNonbondedIxn(const Lepton::CompiledExpression& energyExpression,
        const Lepton::CompiledExpression& forceExpression, const vector<string>& parameterNames,
        const vector<Lepton::CompiledExpression> energyParamDerivExpressions) :
            cutoff(false), useSwitch(false), periodic(false), energyExpression(energyExpression), forceExpression(forceExpression),
            paramNames(parameterNames), energyParamDerivExpressions(energyParamDerivExpressions) {
    expressionSet.registerExpression(this->energyExpression);
    expressionSet.registerExpression(this->forceExpression);
    for (int i = 0; i < this->energyParamDerivExpressions.size(); i++)
        expressionSet.registerExpression(this->energyParamDerivExpressions[i]);
    rIndex = expressionSet.getVariableIndex("r");
    for (auto& param : paramNames) {
        for (int j = 1; j < 3; j++) {
            stringstream name;
            name << param << j;
            particleParamIndex.push_back(expressionSet.getVariableIndex(name.str()));
        }
    }
}

/**---------------------------------------------------------------------------------------

   ReferenceCustomNonbondedIxn destructor

   --------------------------------------------------------------------------------------- */

ReferenceCustomNonbondedIxn::~ReferenceCustomNonbondedIxn() {
}

  /**---------------------------------------------------------------------------------------

     Set the force to use a cutoff.

     @param distance            the cutoff distance
     @param neighbors           the neighbor list to use

     --------------------------------------------------------------------------------------- */

  void ReferenceCustomNonbondedIxn::setUseCutoff(double distance, const OpenMM::NeighborList& neighbors) {

    cutoff = true;
    cutoffDistance = distance;
    neighborList = &neighbors;
  }

/**---------------------------------------------------------------------------------------

   Restrict the force to a list of interaction groups.

   @param distance            the cutoff distance
   @param neighbors           the neighbor list to use

   --------------------------------------------------------------------------------------- */

void ReferenceCustomNonbondedIxn::setInteractionGroups(const vector<pair<set<int>, set<int> > >& groups) {
    interactionGroups = groups;
}

/**---------------------------------------------------------------------------------------

   Set the force to use a switching function.

   @param distance            the switching distance

   --------------------------------------------------------------------------------------- */

void ReferenceCustomNonbondedIxn::setUseSwitchingFunction(double distance) {
    useSwitch = true;
    switchingDistance = distance;
}

  /**---------------------------------------------------------------------------------------

     Set the force to use periodic boundary conditions.  This requires that a cutoff has
     also been set, and the smallest side of the periodic box is at least twice the cutoff
     distance.

     @param vectors    the vectors defining the periodic box

     --------------------------------------------------------------------------------------- */

  void ReferenceCustomNonbondedIxn::setPeriodic(OpenMM::Vec3* vectors) {

    assert(cutoff);
    assert(vectors[0][0] >= 2.0*cutoffDistance);
    assert(vectors[1][1] >= 2.0*cutoffDistance);
    assert(vectors[2][2] >= 2.0*cutoffDistance);
    periodic = true;
    periodicBoxVectors[0] = vectors[0];
    periodicBoxVectors[1] = vectors[1];
    periodicBoxVectors[2] = vectors[2];

  }


/**---------------------------------------------------------------------------------------

   Calculate the custom pair ixn

   @param numberOfAtoms    number of atoms
   @param atomCoordinates  atom coordinates
   @param atomParameters   atom parameters                             atomParameters[atomIndex][paramterIndex]
   @param exclusions       atom exclusion indices
                           exclusions[atomIndex] contains the list of exclusions for that atom
   @param globalParameters the values of global parameters
   @param forces           force array (forces added)
   @param totalEnergy      total energy

   --------------------------------------------------------------------------------------- */

void ReferenceCustomNonbondedIxn::calculatePairIxn(int numberOfAtoms, vector<Vec3>& atomCoordinates,
                                             vector<vector<double> >& atomParameters, vector<set<int> >& exclusions,
                                             const map<string, double>& globalParameters, vector<Vec3>& forces,
                                             double* totalEnergy, double* energyParamDerivs) {

    for (auto& param : globalParameters)
        expressionSet.setVariable(expressionSet.getVariableIndex(param.first), param.second);
    if (interactionGroups.size() > 0) { 
        // The user has specified interaction groups, so compute only the requested interactions.
        
        for (auto& group : interactionGroups) {
            const set<int>& set1 = group.first;
            const set<int>& set2 = group.second;
            for (set<int>::const_iterator atom1 = set1.begin(); atom1 != set1.end(); ++atom1) {
                for (set<int>::const_iterator atom2 = set2.begin(); atom2 != set2.end(); ++atom2) {
                    if (*atom1 == *atom2 || exclusions[*atom1].find(*atom2) != exclusions[*atom1].end())
                        continue; // This is an excluded interaction.
                    if (*atom1 > *atom2 && set1.find(*atom2) != set1.end() && set2.find(*atom1) != set2.end())
                        continue; // Both atoms are in both sets, so skip duplicate interactions.
                    for (int j = 0; j < (int) paramNames.size(); j++) {
                        expressionSet.setVariable(particleParamIndex[j*2], atomParameters[*atom1][j]);
                        expressionSet.setVariable(particleParamIndex[j*2+1], atomParameters[*atom2][j]);
                    }
                    calculateOneIxn(*atom1, *atom2, atomCoordinates, forces, totalEnergy, energyParamDerivs);
                }
            }
        }
    }
    else if (cutoff) {
        // We are using a cutoff, so get the interactions from the neighbor list.
        
        for (auto& pair : *neighborList) {
            for (int j = 0; j < (int) paramNames.size(); j++) {
                expressionSet.setVariable(particleParamIndex[j*2], atomParameters[pair.first][j]);
                expressionSet.setVariable(particleParamIndex[j*2+1], atomParameters[pair.second][j]);
            }
            calculateOneIxn(pair.first, pair.second, atomCoordinates, forces, totalEnergy, energyParamDerivs);
        }
    }
    else {
        // Every particle interacts with every other one.
        
        for (int ii = 0; ii < numberOfAtoms; ii++) {
            for (int jj = ii+1; jj < numberOfAtoms; jj++) {
                if (exclusions[jj].find(ii) == exclusions[jj].end()) {
                    for (int j = 0; j < (int) paramNames.size(); j++) {
                        expressionSet.setVariable(particleParamIndex[j*2], atomParameters[ii][j]);
                        expressionSet.setVariable(particleParamIndex[j*2+1], atomParameters[jj][j]);
                    }
                    calculateOneIxn(ii, jj, atomCoordinates, forces, totalEnergy, energyParamDerivs);
                }
            }
        }
    }
}

  /**---------------------------------------------------------------------------------------

     Calculate one pair ixn between two atoms

     @param ii               the index of the first atom
     @param jj               the index of the second atom
     @param atomCoordinates  atom coordinates
     @param forces           force array (forces added)
     @param totalEnergy      total energy

     --------------------------------------------------------------------------------------- */

void ReferenceCustomNonbondedIxn::calculateOneIxn(int ii, int jj, vector<Vec3>& atomCoordinates, vector<Vec3>& forces,
                        double* totalEnergy, double* energyParamDerivs) {
    // get deltaR, R2, and R between 2 atoms

    double deltaR[ReferenceForce::LastDeltaRIndex];
    if (periodic)
        ReferenceForce::getDeltaRPeriodic(atomCoordinates[jj], atomCoordinates[ii], periodicBoxVectors, deltaR);
    else
        ReferenceForce::getDeltaR(atomCoordinates[jj], atomCoordinates[ii], deltaR);
    double r = deltaR[ReferenceForce::RIndex];
    if (cutoff && r >= cutoffDistance)
        return;

    // accumulate forces

    expressionSet.setVariable(rIndex, r);
    double dEdR = forceExpression.evaluate()/(deltaR[ReferenceForce::RIndex]);
    double energy = energyExpression.evaluate();
    double switchValue = 1.0;
    if (useSwitch) {
        if (r > switchingDistance) {
            double t = (r-switchingDistance)/(cutoffDistance-switchingDistance);
            switchValue = 1+t*t*t*(-10+t*(15-t*6));
            double switchDeriv = t*t*(-30+t*(60-t*30))/(cutoffDistance-switchingDistance);
            dEdR = switchValue*dEdR + energy*switchDeriv/r;
            energy *= switchValue;
        }
    }
    for (int kk = 0; kk < 3; kk++) {
       double force = -dEdR*deltaR[kk];
       forces[ii][kk] += force;
       forces[jj][kk] -= force;
    }
    for (int i = 0; i < energyParamDerivExpressions.size(); i++)
        energyParamDerivs[i] += switchValue*energyParamDerivExpressions[i].evaluate();

    // accumulate energies

    if (totalEnergy)
        *totalEnergy += energy;
  }


