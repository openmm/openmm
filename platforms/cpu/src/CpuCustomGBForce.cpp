
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

void CpuCustomGBForce::setUseCutoff(float distance, const CpuNeighborList& neighbors) {
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

void CpuCustomGBForce::calculateIxn(int numberOfAtoms, float* posq, RealOpenMM** atomParameters,
                                           const vector<set<int> >& exclusions, map<string, double>& globalParameters, float* forces,
                                           double* totalEnergy) {
    fvec4 boxSize(periodicBoxSize[0], periodicBoxSize[1], periodicBoxSize[2], 0);
    fvec4 invBoxSize((1/periodicBoxSize[0]), (1/periodicBoxSize[1]), (1/periodicBoxSize[2]), 0);
    for (map<string, double>::const_iterator iter = globalParameters.begin(); iter != globalParameters.end(); ++iter)
        expressionSet.setVariable(expressionSet.getVariableIndex(iter->first), iter->second);

    // First calculate the computed values.

    int numValues = valueTypes.size();
    for (int valueIndex = 0; valueIndex < numValues; valueIndex++) {
        if (valueTypes[valueIndex] == CustomGBForce::SingleParticle)
            calculateSingleParticleValue(valueIndex, numberOfAtoms, posq, values, atomParameters);
        else if (valueTypes[valueIndex] == CustomGBForce::ParticlePair)
            calculateParticlePairValue(valueIndex, numberOfAtoms, posq, atomParameters, values, exclusions, true, boxSize, invBoxSize);
        else
            calculateParticlePairValue(valueIndex, numberOfAtoms, posq, atomParameters, values, exclusions, false, boxSize, invBoxSize);
    }

    // Now calculate the energy and its derivatives.

    for (int i = 0; i < (int) dEdV.size(); i++)
        for (int j = 0; j < (int) dEdV[i].size(); j++)
            dEdV[i][j] = 0.0;
    for (int termIndex = 0; termIndex < (int) energyExpressions.size(); termIndex++) {
        if (energyTypes[termIndex] == CustomGBForce::SingleParticle)
            calculateSingleParticleEnergyTerm(termIndex, numberOfAtoms, posq, values, atomParameters, forces, totalEnergy, dEdV);
        else if (energyTypes[termIndex] == CustomGBForce::ParticlePair)
            calculateParticlePairEnergyTerm(termIndex, numberOfAtoms, posq, atomParameters, values, exclusions, true, forces, totalEnergy, dEdV, boxSize, invBoxSize);
        else
            calculateParticlePairEnergyTerm(termIndex, numberOfAtoms, posq, atomParameters, values, exclusions, false, forces, totalEnergy, dEdV, boxSize, invBoxSize);
    }

    // Apply the chain rule to evaluate forces.

    calculateChainRuleForces(numberOfAtoms, posq, atomParameters, values, exclusions, forces, dEdV, boxSize, invBoxSize);
}

void CpuCustomGBForce::calculateSingleParticleValue(int index, int numAtoms, float* posq, vector<vector<float> >& values,
        RealOpenMM** atomParameters) {
    values[index].resize(numAtoms);
    for (int i = 0; i < numAtoms; i++) {
        expressionSet.setVariable(xindex, posq[4*i]);
        expressionSet.setVariable(yindex, posq[4*i+1]);
        expressionSet.setVariable(zindex, posq[4*i+2]);
        for (int j = 0; j < (int) paramNames.size(); j++)
            expressionSet.setVariable(paramIndex[j], atomParameters[i][j]);
        for (int j = 0; j < index; j++)
            expressionSet.setVariable(valueIndex[j], values[j][i]);
        values[index][i] = (float) valueExpressions[index].evaluate();
    }
}

void CpuCustomGBForce::calculateParticlePairValue(int index, int numAtoms, float* posq, RealOpenMM** atomParameters,
        vector<vector<float> >& values, const vector<set<int> >& exclusions, bool useExclusions, const fvec4& boxSize, const fvec4& invBoxSize) {
    values[index].resize(numAtoms);
    for (int i = 0; i < numAtoms; i++)
        values[index][i] = 0.0f;
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
                        calculateOnePairValue(index, first, second, posq, atomParameters, values, boxSize, invBoxSize);
                        calculateOnePairValue(index, second, first, posq, atomParameters, values, boxSize, invBoxSize);
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
                calculateOnePairValue(index, i, j, posq, atomParameters, values, boxSize, invBoxSize);
                calculateOnePairValue(index, j, i, posq, atomParameters, values, boxSize, invBoxSize);
           }
        }
    }
}

void CpuCustomGBForce::calculateOnePairValue(int index, int atom1, int atom2, float* posq, RealOpenMM** atomParameters,
        vector<vector<float> >& values, const fvec4& boxSize, const fvec4& invBoxSize) {
    fvec4 deltaR;
    fvec4 pos1(posq+4*atom1);
    fvec4 pos2(posq+4*atom2);
    float r2;
    getDeltaR(pos2, pos1, deltaR, r2, periodic, boxSize, invBoxSize);
    float r = sqrtf(r2);
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
    values[index][atom1] += (float) valueExpressions[index].evaluate();
}

void CpuCustomGBForce::calculateSingleParticleEnergyTerm(int index, int numAtoms, float* posq, const vector<vector<float> >& values,
        RealOpenMM** atomParameters, float* forces, double* totalEnergy,
        vector<vector<float> >& dEdV) {
    for (int i = 0; i < numAtoms; i++) {
        expressionSet.setVariable(xindex, posq[4*i]);
        expressionSet.setVariable(yindex, posq[4*i+1]);
        expressionSet.setVariable(zindex, posq[4*i+2]);
        for (int j = 0; j < (int) paramNames.size(); j++)
            expressionSet.setVariable(paramIndex[j], atomParameters[i][j]);
        for (int j = 0; j < (int) valueNames.size(); j++)
            expressionSet.setVariable(valueIndex[j], values[j][i]);
        if (totalEnergy != NULL)
            *totalEnergy += (float) energyExpressions[index].evaluate();
        for (int j = 0; j < (int) valueNames.size(); j++)
            dEdV[j][i] += (float) energyDerivExpressions[index][j].evaluate();
        forces[4*i+0] -= (float) energyGradientExpressions[index][0].evaluate();
        forces[4*i+1] -= (float) energyGradientExpressions[index][1].evaluate();
        forces[4*i+2] -= (float) energyGradientExpressions[index][2].evaluate();
    }
}

void CpuCustomGBForce::calculateParticlePairEnergyTerm(int index, int numAtoms, float* posq, RealOpenMM** atomParameters,
        const vector<vector<float> >& values, const vector<set<int> >& exclusions, bool useExclusions,
        float* forces, double* totalEnergy, vector<vector<float> >& dEdV, const fvec4& boxSize, const fvec4& invBoxSize) {
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
                        calculateOnePairEnergyTerm(index, first, second, posq, atomParameters, values, forces, totalEnergy, dEdV, boxSize, invBoxSize);
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
                calculateOnePairEnergyTerm(index, i, j, posq, atomParameters, values, forces, totalEnergy, dEdV, boxSize, invBoxSize);
           }
        }
    }
}

void CpuCustomGBForce::calculateOnePairEnergyTerm(int index, int atom1, int atom2, float* posq, RealOpenMM** atomParameters,
        const vector<vector<float> >& values, float* forces, double* totalEnergy,
        vector<vector<float> >& dEdV, const fvec4& boxSize, const fvec4& invBoxSize) {
    // Compute the displacement.

    fvec4 deltaR;
    fvec4 pos1(posq+4*atom1);
    fvec4 pos2(posq+4*atom2);
    float r2;
    getDeltaR(pos2, pos1, deltaR, r2, periodic, boxSize, invBoxSize);
    float r = sqrtf(r2);
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
        *totalEnergy += (float) energyExpressions[index].evaluate();
    float dEdR = (float) energyDerivExpressions[index][0].evaluate();
    dEdR *= 1/r;
    fvec4 result = deltaR*dEdR;
    (fvec4(forces+4*atom1)-result).store(forces+4*atom1);
    (fvec4(forces+4*atom2)+result).store(forces+4*atom2);
    for (int i = 0; i < (int) valueNames.size(); i++) {
        dEdV[i][atom1] += (float) energyDerivExpressions[index][2*i+1].evaluate();
        dEdV[i][atom2] += (float) energyDerivExpressions[index][2*i+2].evaluate();
    }
}

void CpuCustomGBForce::calculateChainRuleForces(int numAtoms, float* posq, RealOpenMM** atomParameters,
        const vector<vector<float> >& values, const vector<set<int> >& exclusions, float* forces, vector<vector<float> >& dEdV,
        const fvec4& boxSize, const fvec4& invBoxSize) {
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
                        calculateOnePairChainRule(first, second, posq, atomParameters, values, forces, dEdV, isExcluded, boxSize, invBoxSize);
                        calculateOnePairChainRule(second, first, posq, atomParameters, values, forces, dEdV, isExcluded, boxSize, invBoxSize);
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
                calculateOnePairChainRule(i, j, posq, atomParameters, values, forces, dEdV, isExcluded, boxSize, invBoxSize);
                calculateOnePairChainRule(j, i, posq, atomParameters, values, forces, dEdV, isExcluded, boxSize, invBoxSize);
           }
        }
    }

    // Compute chain rule terms for computed values that depend explicitly on particle coordinates.

    for (int i = 0; i < numAtoms; i++) {
        expressionSet.setVariable(xindex, posq[4*i]);
        expressionSet.setVariable(yindex, posq[4*i+1]);
        expressionSet.setVariable(zindex, posq[4*i+2]);
        for (int j = 0; j < (int) paramNames.size(); j++)
            expressionSet.setVariable(paramIndex[j], atomParameters[i][j]);
        for (int j = 1; j < (int) valueNames.size(); j++) {
            expressionSet.setVariable(valueIndex[j-1], values[j-1][i]);
            dVdX[j] = 0.0;
            dVdY[j] = 0.0;
            dVdZ[j] = 0.0;
            for (int k = 1; k < j; k++) {
                float dVdV = (float) valueDerivExpressions[j][k].evaluate();
                dVdX[j] += dVdV*dVdX[k];
                dVdY[j] += dVdV*dVdY[k];
                dVdZ[j] += dVdV*dVdZ[k];
            }
            dVdX[j] += (float) valueGradientExpressions[j][0].evaluate();
            dVdY[j] += (float) valueGradientExpressions[j][1].evaluate();
            dVdZ[j] += (float) valueGradientExpressions[j][2].evaluate();
            forces[4*i+0] -= dEdV[j][i]*dVdX[j];
            forces[4*i+1] -= dEdV[j][i]*dVdY[j];
            forces[4*i+2] -= dEdV[j][i]*dVdZ[j];
        }
    }
}

void CpuCustomGBForce::calculateOnePairChainRule(int atom1, int atom2, float* posq, RealOpenMM** atomParameters,
        const vector<vector<float> >& values, float* forces, vector<vector<float> >& dEdV, bool isExcluded, const fvec4& boxSize, const fvec4& invBoxSize) {
    // Compute the displacement.

    fvec4 deltaR;
    fvec4 pos1(posq+4*atom1);
    fvec4 pos2(posq+4*atom2);
    float r2;
    getDeltaR(pos2, pos1, deltaR, r2, periodic, boxSize, invBoxSize);
    float r = sqrtf(r2);
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

    float rinv = 1/r;
    deltaR *= rinv;
    if (!isExcluded || valueTypes[0] != CustomGBForce::ParticlePair) {
        dVdR1[0] = (float) valueDerivExpressions[0][0].evaluate();
        dVdR2[0] = -dVdR1[0];
        (fvec4(forces+4*atom1)-deltaR*(dEdV[0][atom1]*dVdR1[0])).store(forces+4*atom1);
        (fvec4(forces+4*atom2)-deltaR*(dEdV[0][atom1]*dVdR2[0])).store(forces+4*atom2);
    }
    for (int i = 0; i < (int) paramNames.size(); i++)
        expressionSet.setVariable(paramIndex[i], atomParameters[atom1][i]);
    expressionSet.setVariable(valueIndex[0], values[0][atom1]);
    for (int i = 1; i < (int) valueNames.size(); i++) {
        expressionSet.setVariable(valueIndex[i], values[i][atom1]);
        expressionSet.setVariable(xindex, posq[4*atom1]);
        expressionSet.setVariable(yindex, posq[4*atom1+1]);
        expressionSet.setVariable(zindex, posq[4*atom1+2]);
        dVdR1[i] = 0.0;
        dVdR2[i] = 0.0;
        for (int j = 0; j < i; j++) {
            float dVdV = (float) valueDerivExpressions[i][j].evaluate();
            dVdR1[i] += dVdV*dVdR1[j];
            dVdR2[i] += dVdV*dVdR2[j];
        }
        (fvec4(forces+4*atom1)-deltaR*(dEdV[i][atom1]*dVdR1[i])).store(forces+4*atom1);
        (fvec4(forces+4*atom2)-deltaR*(dEdV[i][atom1]*dVdR2[i])).store(forces+4*atom2);
    }
}

void CpuCustomGBForce::getDeltaR(const fvec4& posI, const fvec4& posJ, fvec4& deltaR, float& r2, bool periodic, const fvec4& boxSize, const fvec4& invBoxSize) const {
    deltaR = posJ-posI;
    if (periodic) {
        fvec4 base = round(deltaR*invBoxSize)*boxSize;
        deltaR = deltaR-base;
    }
    r2 = dot3(deltaR, deltaR);
}
