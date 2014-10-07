
/* Portions copyright (c) 2009-2014 Stanford University and Simbios.
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

#include "SimTKOpenMMCommon.h"
#include "SimTKOpenMMLog.h"
#include "SimTKOpenMMUtilities.h"
#include "ReferenceForce.h"
#include "CpuCustomGBForce.h"

using namespace OpenMM;
using namespace std;

CpuCustomGBForce::CpuCustomGBForce(int numAtoms, const vector<Lepton::CompiledExpression>& valueExpressions,
                     const vector<vector<Lepton::CompiledExpression> >& valueDerivExpressions,
                     const vector<vector<Lepton::CompiledExpression> >& valueGradientExpressions,
                     const vector<string>& valueNames,
                     const vector<CustomGBForce::ComputationType>& valueTypes,
                     const vector<Lepton::CompiledExpression>& energyExpressions,
                     const vector<vector<Lepton::CompiledExpression> >& energyDerivExpressions,
                     const vector<vector<Lepton::CompiledExpression> >& energyGradientExpressions,
                     const vector<CustomGBForce::ComputationType>& energyTypes,
                     const vector<string>& parameterNames) :
            cutoff(false), periodic(false), valueExpressions(valueExpressions), valueDerivExpressions(valueDerivExpressions), valueGradientExpressions(valueGradientExpressions),
            valueNames(valueNames), valueTypes(valueTypes), energyExpressions(energyExpressions), energyDerivExpressions(energyDerivExpressions), energyGradientExpressions(energyGradientExpressions),
            energyTypes(energyTypes), paramNames(parameterNames) {
    for (int i = 0; i < (int) valueExpressions.size(); i++)
        expressionSet.registerExpression(this->valueExpressions[i]);
    for (int i = 0; i < (int) valueDerivExpressions.size(); i++)
        for (int j = 0; j < (int) valueDerivExpressions[i].size(); j++)
            expressionSet.registerExpression(this->valueDerivExpressions[i][j]);
    for (int i = 0; i < (int) valueGradientExpressions.size(); i++)
        for (int j = 0; j < (int) valueGradientExpressions[i].size(); j++)
            expressionSet.registerExpression(this->valueGradientExpressions[i][j]);
    for (int i = 0; i < (int) energyExpressions.size(); i++)
        expressionSet.registerExpression(this->energyExpressions[i]);
    for (int i = 0; i < (int) energyDerivExpressions.size(); i++)
        for (int j = 0; j < (int) energyDerivExpressions[i].size(); j++)
            expressionSet.registerExpression(this->energyDerivExpressions[i][j]);
    for (int i = 0; i < (int) energyGradientExpressions.size(); i++)
        for (int j = 0; j < (int) energyGradientExpressions[i].size(); j++)
            expressionSet.registerExpression(this->energyGradientExpressions[i][j]);
    xindex = expressionSet.getVariableIndex("x");
    yindex = expressionSet.getVariableIndex("y");
    zindex = expressionSet.getVariableIndex("z");
    rindex = expressionSet.getVariableIndex("r");
    for (int i = 0; i < (int) paramNames.size(); i++) {
        paramIndex.push_back(expressionSet.getVariableIndex(paramNames[i]));
        for (int j = 1; j < 3; j++) {
            stringstream name;
            name << paramNames[i] << j;
            particleParamIndex.push_back(expressionSet.getVariableIndex(name.str()));
        }
    }
    for (int i = 0; i < (int) valueNames.size(); i++) {
        valueIndex.push_back(expressionSet.getVariableIndex(valueNames[i]));
        for (int j = 1; j < 3; j++) {
            stringstream name;
            name << valueNames[i] << j;
            particleValueIndex.push_back(expressionSet.getVariableIndex(name.str()));
        }
    }
    values.resize(valueTypes.size());
    dEdV.resize(valueTypes.size());
    for (int i = 0; i < (int) values.size(); i++) {
        values[i].resize(numAtoms);
        dEdV[i].resize(numAtoms);
    }
    dVdX.resize(valueDerivExpressions.size());
    dVdY.resize(valueDerivExpressions.size());
    dVdZ.resize(valueDerivExpressions.size());
    dVdR1.resize(valueDerivExpressions.size());
    dVdR2.resize(valueDerivExpressions.size());
}

CpuCustomGBForce::~CpuCustomGBForce() {
}

void CpuCustomGBForce::setUseCutoff(RealOpenMM distance, const CpuNeighborList& neighbors) {
    cutoff = true;
    cutoffDistance = distance;
    neighborList = &neighbors;
  }

void CpuCustomGBForce::setPeriodic(RealVec& boxSize) {

    if (cutoff) {
        assert(boxSize[0] >= 2.0*cutoffDistance);
        assert(boxSize[1] >= 2.0*cutoffDistance);
        assert(boxSize[2] >= 2.0*cutoffDistance);
    }
    periodic = true;
    periodicBoxSize[0] = boxSize[0];
    periodicBoxSize[1] = boxSize[1];
    periodicBoxSize[2] = boxSize[2];
  }

void CpuCustomGBForce::calculateIxn(int numberOfAtoms, vector<RealVec>& atomCoordinates, RealOpenMM** atomParameters,
                                           const vector<set<int> >& exclusions, map<string, double>& globalParameters, vector<RealVec>& forces,
                                           RealOpenMM* totalEnergy) {
    for (map<string, double>::const_iterator iter = globalParameters.begin(); iter != globalParameters.end(); ++iter)
        expressionSet.setVariable(expressionSet.getVariableIndex(iter->first), iter->second);

    // First calculate the computed values.

    int numValues = valueTypes.size();
    for (int valueIndex = 0; valueIndex < numValues; valueIndex++) {
        if (valueTypes[valueIndex] == CustomGBForce::SingleParticle)
            calculateSingleParticleValue(valueIndex, numberOfAtoms, atomCoordinates, values, atomParameters);
        else if (valueTypes[valueIndex] == CustomGBForce::ParticlePair)
            calculateParticlePairValue(valueIndex, numberOfAtoms, atomCoordinates, atomParameters, values, exclusions, true);
        else
            calculateParticlePairValue(valueIndex, numberOfAtoms, atomCoordinates, atomParameters, values, exclusions, false);
    }

    // Now calculate the energy and its derivatives.

    for (int i = 0; i < (int) dEdV.size(); i++)
        for (int j = 0; j < (int) dEdV[i].size(); j++)
            dEdV[i][j] = 0.0;
    for (int termIndex = 0; termIndex < (int) energyExpressions.size(); termIndex++) {
        if (energyTypes[termIndex] == CustomGBForce::SingleParticle)
            calculateSingleParticleEnergyTerm(termIndex, numberOfAtoms, atomCoordinates, values, atomParameters, forces, totalEnergy, dEdV);
        else if (energyTypes[termIndex] == CustomGBForce::ParticlePair)
            calculateParticlePairEnergyTerm(termIndex, numberOfAtoms, atomCoordinates, atomParameters, values, exclusions, true, forces, totalEnergy, dEdV);
        else
            calculateParticlePairEnergyTerm(termIndex, numberOfAtoms, atomCoordinates, atomParameters, values, exclusions, false, forces, totalEnergy, dEdV);
    }

    // Apply the chain rule to evaluate forces.

    calculateChainRuleForces(numberOfAtoms, atomCoordinates, atomParameters, values, exclusions, forces, dEdV);
}

void CpuCustomGBForce::calculateSingleParticleValue(int index, int numAtoms, vector<RealVec>& atomCoordinates, vector<vector<RealOpenMM> >& values,
        RealOpenMM** atomParameters) {
    values[index].resize(numAtoms);
    for (int i = 0; i < numAtoms; i++) {
        expressionSet.setVariable(xindex, atomCoordinates[i][0]);
        expressionSet.setVariable(yindex, atomCoordinates[i][1]);
        expressionSet.setVariable(zindex, atomCoordinates[i][2]);
        for (int j = 0; j < (int) paramNames.size(); j++)
            expressionSet.setVariable(paramIndex[j], atomParameters[i][j]);
        for (int j = 0; j < index; j++)
            expressionSet.setVariable(valueIndex[j], values[j][i]);
        values[index][i] = (RealOpenMM) valueExpressions[index].evaluate();
    }
}

void CpuCustomGBForce::calculateParticlePairValue(int index, int numAtoms, vector<RealVec>& atomCoordinates, RealOpenMM** atomParameters,
        vector<vector<RealOpenMM> >& values, const vector<set<int> >& exclusions, bool useExclusions) {
    values[index].resize(numAtoms);
    for (int i = 0; i < numAtoms; i++)
        values[index][i] = (RealOpenMM) 0.0;
    if (cutoff) {
        // Loop over all pairs in the neighbor list.

        for (int blockIndex = 0; blockIndex < neighborList->getNumBlocks(); blockIndex++) {
            const int* blockAtom = &neighborList->getSortedAtoms()[4*blockIndex];
            const vector<int>& neighbors = neighborList->getBlockNeighbors(blockIndex);
            const vector<char>& blockExclusions = neighborList->getBlockExclusions(blockIndex);
            for (int i = 0; i < (int) neighbors.size(); i++) {
                int first = neighbors[i];
                for (int k = 0; k < 4; k++) {
                    if ((blockExclusions[i] & (1<<k)) == 0) {
                        int second = blockAtom[k];
                        if (useExclusions && exclusions[first].find(second) != exclusions[first].end())
                            continue;
                        calculateOnePairValue(index, first, second, atomCoordinates, atomParameters, values);
                        calculateOnePairValue(index, second, first, atomCoordinates, atomParameters, values);
                    }
                }
            }
        }
    }
    else {
        // Perform an O(N^2) loop over all atom pairs.

        for (int i = 0; i < numAtoms; i++) {
            for (int j = i+1; j < numAtoms; j++) {
                if (useExclusions && exclusions[i].find(j) != exclusions[i].end())
                    continue;
                calculateOnePairValue(index, i, j, atomCoordinates, atomParameters, values);
                calculateOnePairValue(index, j, i, atomCoordinates, atomParameters, values);
           }
        }
    }
}

void CpuCustomGBForce::calculateOnePairValue(int index, int atom1, int atom2, vector<RealVec>& atomCoordinates, RealOpenMM** atomParameters,
        vector<vector<RealOpenMM> >& values) {
    RealOpenMM deltaR[ReferenceForce::LastDeltaRIndex];
    if (periodic)
        ReferenceForce::getDeltaRPeriodic(atomCoordinates[atom2], atomCoordinates[atom1], periodicBoxSize, deltaR);
    else
        ReferenceForce::getDeltaR(atomCoordinates[atom2], atomCoordinates[atom1], deltaR);
    RealOpenMM r = deltaR[ReferenceForce::RIndex];
    if (cutoff && r >= cutoffDistance)
        return;
    for (int i = 0; i < (int) paramNames.size(); i++) {
        expressionSet.setVariable(particleParamIndex[i*2], atomParameters[atom1][i]);
        expressionSet.setVariable(particleParamIndex[i*2+1], atomParameters[atom2][i]);
    }
    expressionSet.setVariable(rindex, r);
    for (int i = 0; i < index; i++) {
        expressionSet.setVariable(particleValueIndex[i*2], values[i][atom1]);
        expressionSet.setVariable(particleValueIndex[i*2+1], values[i][atom2]);
    }
    values[index][atom1] += (RealOpenMM) valueExpressions[index].evaluate();
}

void CpuCustomGBForce::calculateSingleParticleEnergyTerm(int index, int numAtoms, vector<RealVec>& atomCoordinates, const vector<vector<RealOpenMM> >& values,
        RealOpenMM** atomParameters, vector<RealVec>& forces, RealOpenMM* totalEnergy,
        vector<vector<RealOpenMM> >& dEdV) {
    for (int i = 0; i < numAtoms; i++) {
        expressionSet.setVariable(xindex, atomCoordinates[i][0]);
        expressionSet.setVariable(yindex, atomCoordinates[i][1]);
        expressionSet.setVariable(zindex, atomCoordinates[i][2]);
        for (int j = 0; j < (int) paramNames.size(); j++)
            expressionSet.setVariable(paramIndex[j], atomParameters[i][j]);
        for (int j = 0; j < (int) valueNames.size(); j++)
            expressionSet.setVariable(valueIndex[j], values[j][i]);
        if (totalEnergy != NULL)
            *totalEnergy += (RealOpenMM) energyExpressions[index].evaluate();
        for (int j = 0; j < (int) valueNames.size(); j++)
            dEdV[j][i] += (RealOpenMM) energyDerivExpressions[index][j].evaluate();
        forces[i][0] -= (RealOpenMM) energyGradientExpressions[index][0].evaluate();
        forces[i][1] -= (RealOpenMM) energyGradientExpressions[index][1].evaluate();
        forces[i][2] -= (RealOpenMM) energyGradientExpressions[index][2].evaluate();
    }
}

void CpuCustomGBForce::calculateParticlePairEnergyTerm(int index, int numAtoms, vector<RealVec>& atomCoordinates, RealOpenMM** atomParameters,
        const vector<vector<RealOpenMM> >& values, const vector<set<int> >& exclusions, bool useExclusions,
        vector<RealVec>& forces, RealOpenMM* totalEnergy, vector<vector<RealOpenMM> >& dEdV) {
    if (cutoff) {
        // Loop over all pairs in the neighbor list.

        for (int blockIndex = 0; blockIndex < neighborList->getNumBlocks(); blockIndex++) {
            const int* blockAtom = &neighborList->getSortedAtoms()[4*blockIndex];
            const vector<int>& neighbors = neighborList->getBlockNeighbors(blockIndex);
            const vector<char>& blockExclusions = neighborList->getBlockExclusions(blockIndex);
            for (int i = 0; i < (int) neighbors.size(); i++) {
                int first = neighbors[i];
                for (int k = 0; k < 4; k++) {
                    if ((blockExclusions[i] & (1<<k)) == 0) {
                        int second = blockAtom[k];
                        if (useExclusions && exclusions[first].find(second) != exclusions[first].end())
                            continue;
                        calculateOnePairEnergyTerm(index, first, second, atomCoordinates, atomParameters, values, forces, totalEnergy, dEdV);
                    }
                }
            }
        }
    }
    else {
        // Perform an O(N^2) loop over all atom pairs.

        for (int i = 0; i < numAtoms; i++) {
            for (int j = i+1; j < numAtoms; j++) {
                if (useExclusions && exclusions[i].find(j) != exclusions[i].end())
                    continue;
                calculateOnePairEnergyTerm(index, i, j, atomCoordinates, atomParameters, values, forces, totalEnergy, dEdV);
           }
        }
    }
}

void CpuCustomGBForce::calculateOnePairEnergyTerm(int index, int atom1, int atom2, vector<RealVec>& atomCoordinates, RealOpenMM** atomParameters,
        const vector<vector<RealOpenMM> >& values, vector<RealVec>& forces, RealOpenMM* totalEnergy,
        vector<vector<RealOpenMM> >& dEdV) {
    // Compute the displacement.

    RealOpenMM deltaR[ReferenceForce::LastDeltaRIndex];
    if (periodic)
        ReferenceForce::getDeltaRPeriodic(atomCoordinates[atom2], atomCoordinates[atom1], periodicBoxSize, deltaR);
    else
        ReferenceForce::getDeltaR(atomCoordinates[atom2], atomCoordinates[atom1], deltaR);
    RealOpenMM r = deltaR[ReferenceForce::RIndex];
    if (cutoff && r >= cutoffDistance)
        return;

    // Record variables for evaluating expressions.

    for (int i = 0; i < (int) paramNames.size(); i++) {
        expressionSet.setVariable(particleParamIndex[i*2], atomParameters[atom1][i]);
        expressionSet.setVariable(particleParamIndex[i*2+1], atomParameters[atom2][i]);
    }
    expressionSet.setVariable(rindex, r);
    for (int i = 0; i < (int) valueNames.size(); i++) {
        expressionSet.setVariable(particleValueIndex[i*2], values[i][atom1]);
        expressionSet.setVariable(particleValueIndex[i*2+1], values[i][atom2]);
    }

    // Evaluate the energy and its derivatives.

    if (totalEnergy != NULL)
        *totalEnergy += (RealOpenMM) energyExpressions[index].evaluate();
    RealOpenMM dEdR = (RealOpenMM) energyDerivExpressions[index][0].evaluate();
    dEdR *= 1/r;
    for (int i = 0; i < 3; i++) {
       forces[atom1][i] -= dEdR*deltaR[i];
       forces[atom2][i] += dEdR*deltaR[i];
    }
    for (int i = 0; i < (int) valueNames.size(); i++) {
        dEdV[i][atom1] += (RealOpenMM) energyDerivExpressions[index][2*i+1].evaluate();
        dEdV[i][atom2] += (RealOpenMM) energyDerivExpressions[index][2*i+2].evaluate();
    }
}

void CpuCustomGBForce::calculateChainRuleForces(int numAtoms, vector<RealVec>& atomCoordinates, RealOpenMM** atomParameters,
        const vector<vector<RealOpenMM> >& values, const vector<set<int> >& exclusions, vector<RealVec>& forces, vector<vector<RealOpenMM> >& dEdV) {
    if (cutoff) {
        // Loop over all pairs in the neighbor list.

        for (int blockIndex = 0; blockIndex < neighborList->getNumBlocks(); blockIndex++) {
            const int* blockAtom = &neighborList->getSortedAtoms()[4*blockIndex];
            const vector<int>& neighbors = neighborList->getBlockNeighbors(blockIndex);
            const vector<char>& blockExclusions = neighborList->getBlockExclusions(blockIndex);
            for (int i = 0; i < (int) neighbors.size(); i++) {
                int first = neighbors[i];
                for (int k = 0; k < 4; k++) {
                    if ((blockExclusions[i] & (1<<k)) == 0) {
                        int second = blockAtom[k];
                        bool isExcluded = (exclusions[first].find(second) != exclusions[first].end());
                        calculateOnePairChainRule(first, second, atomCoordinates, atomParameters, values, forces, dEdV, isExcluded);
                        calculateOnePairChainRule(second, first, atomCoordinates, atomParameters, values, forces, dEdV, isExcluded);
                    }
                }
            }
        }
    }
    else {
        // Perform an O(N^2) loop over all atom pairs.

        for (int i = 0; i < numAtoms; i++) {
            for (int j = i+1; j < numAtoms; j++) {
                bool isExcluded = (exclusions[i].find(j) != exclusions[i].end());
                calculateOnePairChainRule(i, j, atomCoordinates, atomParameters, values, forces, dEdV, isExcluded);
                calculateOnePairChainRule(j, i, atomCoordinates, atomParameters, values, forces, dEdV, isExcluded);
           }
        }
    }

    // Compute chain rule terms for computed values that depend explicitly on particle coordinates.

    for (int i = 0; i < numAtoms; i++) {
        expressionSet.setVariable(xindex, atomCoordinates[i][0]);
        expressionSet.setVariable(yindex, atomCoordinates[i][1]);
        expressionSet.setVariable(zindex, atomCoordinates[i][2]);
        for (int j = 0; j < (int) paramNames.size(); j++)
            expressionSet.setVariable(paramIndex[j], atomParameters[i][j]);
        for (int j = 1; j < (int) valueNames.size(); j++) {
            expressionSet.setVariable(valueIndex[j-1], values[j-1][i]);
            dVdX[j] = 0.0;
            dVdY[j] = 0.0;
            dVdZ[j] = 0.0;
            for (int k = 1; k < j; k++) {
                RealOpenMM dVdV = (RealOpenMM) valueDerivExpressions[j][k].evaluate();
                dVdX[j] += dVdV*dVdX[k];
                dVdY[j] += dVdV*dVdY[k];
                dVdZ[j] += dVdV*dVdZ[k];
            }
            dVdX[j] += (RealOpenMM) valueGradientExpressions[j][0].evaluate();
            dVdY[j] += (RealOpenMM) valueGradientExpressions[j][1].evaluate();
            dVdZ[j] += (RealOpenMM) valueGradientExpressions[j][2].evaluate();
            forces[i][0] -= dEdV[j][i]*dVdX[j];
            forces[i][1] -= dEdV[j][i]*dVdY[j];
            forces[i][2] -= dEdV[j][i]*dVdZ[j];
        }
    }
}

void CpuCustomGBForce::calculateOnePairChainRule(int atom1, int atom2, vector<RealVec>& atomCoordinates, RealOpenMM** atomParameters,
        const vector<vector<RealOpenMM> >& values, vector<RealVec>& forces, vector<vector<RealOpenMM> >& dEdV, bool isExcluded) {
    // Compute the displacement.

    RealOpenMM deltaR[ReferenceForce::LastDeltaRIndex];
    if (periodic)
        ReferenceForce::getDeltaRPeriodic(atomCoordinates[atom2], atomCoordinates[atom1], periodicBoxSize, deltaR);
    else
        ReferenceForce::getDeltaR(atomCoordinates[atom2], atomCoordinates[atom1], deltaR);
    RealOpenMM r = deltaR[ReferenceForce::RIndex];
    if (cutoff && r >= cutoffDistance)
        return;

    // Record variables for evaluating expressions.

    for (int i = 0; i < (int) paramNames.size(); i++) {
        expressionSet.setVariable(particleParamIndex[i*2], atomParameters[atom1][i]);
        expressionSet.setVariable(particleParamIndex[i*2+1], atomParameters[atom2][i]);
    }
    expressionSet.setVariable(rindex, r);
    expressionSet.setVariable(particleValueIndex[0], values[0][atom1]);
    expressionSet.setVariable(particleValueIndex[1], values[0][atom2]);

    // Evaluate the derivative of each parameter with respect to position and apply forces.

    RealOpenMM rinv = 1/r;
    deltaR[0] *= rinv;
    deltaR[1] *= rinv;
    deltaR[2] *= rinv;
    if (!isExcluded || valueTypes[0] != CustomGBForce::ParticlePair) {
        dVdR1[0] = (RealOpenMM) valueDerivExpressions[0][0].evaluate();
        dVdR2[0] = -dVdR1[0];
        for (int i = 0; i < 3; i++) {
            forces[atom1][i] -= dEdV[0][atom1]*dVdR1[0]*deltaR[i];
            forces[atom2][i] -= dEdV[0][atom1]*dVdR2[0]*deltaR[i];
        }
    }
    for (int i = 0; i < (int) paramNames.size(); i++)
        expressionSet.setVariable(paramIndex[i], atomParameters[atom1][i]);
    expressionSet.setVariable(valueIndex[0], values[0][atom1]);
    for (int i = 1; i < (int) valueNames.size(); i++) {
        expressionSet.setVariable(valueIndex[i], values[i][atom1]);
        expressionSet.setVariable(xindex, atomCoordinates[atom1][0]);
        expressionSet.setVariable(yindex, atomCoordinates[atom1][1]);
        expressionSet.setVariable(zindex, atomCoordinates[atom1][2]);
        dVdR1[i] = 0.0;
        dVdR2[i] = 0.0;
        for (int j = 0; j < i; j++) {
            RealOpenMM dVdV = (RealOpenMM) valueDerivExpressions[i][j].evaluate();
            dVdR1[i] += dVdV*dVdR1[j];
            dVdR2[i] += dVdV*dVdR2[j];
        }
        for (int k = 0; k < 3; k++) {
            forces[atom1][k] -= dEdV[i][atom1]*dVdR1[i]*deltaR[k];
            forces[atom2][k] -= dEdV[i][atom1]*dVdR2[i]*deltaR[k];
        }
    }
}
