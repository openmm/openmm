
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
#include "gmx_atomic.h"

using namespace OpenMM;
using namespace std;

class CpuCustomGBForce::ComputeForceTask : public ThreadPool::Task {
public:
    ComputeForceTask(CpuCustomGBForce& owner) : owner(owner) {
    }
    void execute(ThreadPool& threads, int threadIndex) {
        owner.threadComputeForce(threads, threadIndex);
    }
    CpuCustomGBForce& owner;
};

CpuCustomGBForce::ThreadData::ThreadData(int numAtoms, int numThreads, int threadIndex,
                      const vector<Lepton::CompiledExpression>& valueExpressions,
                      const vector<vector<Lepton::CompiledExpression> >& valueDerivExpressions,
                      const vector<vector<Lepton::CompiledExpression> >& valueGradientExpressions,
                      const vector<string>& valueNames,
                      const vector<Lepton::CompiledExpression>& energyExpressions,
                      const vector<vector<Lepton::CompiledExpression> >& energyDerivExpressions,
                      const vector<vector<Lepton::CompiledExpression> >& energyGradientExpressions,
                      const vector<string>& parameterNames) :
            valueExpressions(valueExpressions), valueDerivExpressions(valueDerivExpressions), valueGradientExpressions(valueGradientExpressions),
            energyExpressions(energyExpressions), energyDerivExpressions(energyDerivExpressions), energyGradientExpressions(energyGradientExpressions) {
    firstAtom = (threadIndex*(long long) numAtoms)/numThreads;
    lastAtom = ((threadIndex+1)*(long long) numAtoms)/numThreads;
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
    for (int i = 0; i < (int) parameterNames.size(); i++) {
        paramIndex.push_back(expressionSet.getVariableIndex(parameterNames[i]));
        for (int j = 1; j < 3; j++) {
            stringstream name;
            name << parameterNames[i] << j;
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
    value0.resize(numAtoms);
    dEdV.resize(valueNames.size());
    for (int i = 0; i < (int) dEdV.size(); i++)
        dEdV[i].resize(numAtoms);
    dVdX.resize(valueDerivExpressions.size());
    dVdY.resize(valueDerivExpressions.size());
    dVdZ.resize(valueDerivExpressions.size());
    dVdR1.resize(valueDerivExpressions.size());
    dVdR2.resize(valueDerivExpressions.size());
}

CpuCustomGBForce::CpuCustomGBForce(int numAtoms, const std::vector<std::set<int> >& exclusions,
                     const vector<Lepton::CompiledExpression>& valueExpressions,
                     const vector<vector<Lepton::CompiledExpression> >& valueDerivExpressions,
                     const vector<vector<Lepton::CompiledExpression> >& valueGradientExpressions,
                     const vector<string>& valueNames,
                     const vector<CustomGBForce::ComputationType>& valueTypes,
                     const vector<Lepton::CompiledExpression>& energyExpressions,
                     const vector<vector<Lepton::CompiledExpression> >& energyDerivExpressions,
                     const vector<vector<Lepton::CompiledExpression> >& energyGradientExpressions,
                     const vector<CustomGBForce::ComputationType>& energyTypes,
                     const vector<string>& parameterNames, ThreadPool& threads) :
            exclusions(exclusions), cutoff(false), periodic(false), valueNames(valueNames), valueTypes(valueTypes),
            energyTypes(energyTypes), paramNames(parameterNames), threads(threads) {
    for (int i = 0; i < threads.getNumThreads(); i++)
        threadData.push_back(new ThreadData(numAtoms, threads.getNumThreads(), i, valueExpressions, valueDerivExpressions, valueGradientExpressions, valueNames,
                      energyExpressions, energyDerivExpressions, energyGradientExpressions, parameterNames));
    values.resize(valueNames.size());
    dEdV.resize(valueNames.size());
    for (int i = 0; i < (int) values.size(); i++) {
        values[i].resize(numAtoms);
        dEdV[i].resize(numAtoms);
    }
}

CpuCustomGBForce::~CpuCustomGBForce() {
    for (int i = 0; i < (int) threadData.size(); i++)
        delete threadData[i];
}

void CpuCustomGBForce::setUseCutoff(float distance, const CpuNeighborList& neighbors) {
    cutoff = true;
    cutoffDistance = distance;
    cutoffDistance2 = distance*distance;
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
                                           map<string, double>& globalParameters, vector<AlignedArray<float> >& threadForce,
                                           bool includeForce, bool includeEnergy, double& totalEnergy) {
    // Record the parameters for the threads.
    
    this->numberOfAtoms = numberOfAtoms;
    this->posq = posq;
    this->atomParameters = atomParameters;
    this->globalParameters = &globalParameters;
    this->threadForce = &threadForce;
    this->includeForce = includeForce;
    this->includeEnergy = includeEnergy;
    threadEnergy.resize(threads.getNumThreads());
    gmx_atomic_t counter;
    this->atomicCounter = &counter;

    // Calculate the first computed value.

    ComputeForceTask task(*this);
    gmx_atomic_set(&counter, 0);
    threads.execute(task);
    threads.waitForThreads();
    
    // Calculate the remaining computed values.
    
    threads.resumeThreads();
    threads.waitForThreads();

    // Calculate the energy terms.

    for (int i = 0; i < (int) threadData[0]->energyExpressions.size(); i++) {
        gmx_atomic_set(&counter, 0);
        threads.execute(task);
        threads.waitForThreads();
    }

    // Sum the energy derivatives.

    threads.resumeThreads();
    threads.waitForThreads();
    
    // Apply the chain rule to evaluate forces.

    gmx_atomic_set(&counter, 0);
    threads.resumeThreads();
    threads.waitForThreads();

    // Combine the energies from all the threads.
    
    if (includeEnergy) {
        int numThreads = threads.getNumThreads();
        for (int i = 0; i < numThreads; i++)
            totalEnergy += threadEnergy[i];
    }
}

void CpuCustomGBForce::threadComputeForce(ThreadPool& threads, int threadIndex) {
    // Compute this thread's subset of interactions.

    int numThreads = threads.getNumThreads();
    threadEnergy[threadIndex] = 0;
    double& energy = threadEnergy[threadIndex];
    float* forces = &(*threadForce)[threadIndex][0];
    ThreadData& data = *threadData[threadIndex];
    fvec4 boxSize(periodicBoxSize[0], periodicBoxSize[1], periodicBoxSize[2], 0);
    fvec4 invBoxSize((1/periodicBoxSize[0]), (1/periodicBoxSize[1]), (1/periodicBoxSize[2]), 0);
    for (map<string, double>::const_iterator iter = globalParameters->begin(); iter != globalParameters->end(); ++iter)
        data.expressionSet.setVariable(data.expressionSet.getVariableIndex(iter->first), iter->second);

    // Calculate the first computed value.

    for (int i = 0; i < (int) data.value0.size(); i++)
        data.value0[i] = 0.0f;
    if (valueTypes[0] == CustomGBForce::ParticlePair)
        calculateParticlePairValue(0, data, numberOfAtoms, posq, atomParameters, true, boxSize, invBoxSize);
    else
        calculateParticlePairValue(0, data, numberOfAtoms, posq, atomParameters, false, boxSize, invBoxSize);
    threads.syncThreads();

    // Sum the first computed value and calculate the remaining ones.

    int numValues = valueTypes.size();
    for (int atom = data.firstAtom; atom < data.lastAtom; atom++) {
        float sum = 0.0f;
        for (int j = 0; j < (int) threadData.size(); j++)
            sum += threadData[j]->value0[atom];
        values[0][atom] = sum;
        data.expressionSet.setVariable(data.xindex, posq[4*atom]);
        data.expressionSet.setVariable(data.yindex, posq[4*atom+1]);
        data.expressionSet.setVariable(data.zindex, posq[4*atom+2]);
        for (int j = 0; j < (int) paramNames.size(); j++)
            data.expressionSet.setVariable(data.paramIndex[j], atomParameters[atom][j]);
        for (int i = 1; i < numValues; i++) {
            data.expressionSet.setVariable(data.valueIndex[i-1], values[i-1][atom]);
            values[i][atom] = (float) data.valueExpressions[i].evaluate();
        }
    }
    threads.syncThreads();

    // Now calculate the energy and its derivatives.

    for (int i = 0; i < (int) data.dEdV.size(); i++)
        for (int j = 0; j < (int) data.dEdV[i].size(); j++)
            data.dEdV[i][j] = 0.0;
    for (int termIndex = 0; termIndex < (int) data.energyExpressions.size(); termIndex++) {
        if (energyTypes[termIndex] == CustomGBForce::SingleParticle)
            calculateSingleParticleEnergyTerm(termIndex, data, numberOfAtoms, posq, atomParameters, forces, energy);
        else if (energyTypes[termIndex] == CustomGBForce::ParticlePair)
            calculateParticlePairEnergyTerm(termIndex, data, numberOfAtoms, posq, atomParameters, true, forces, energy, boxSize, invBoxSize);
        else
            calculateParticlePairEnergyTerm(termIndex, data, numberOfAtoms, posq, atomParameters, false, forces, energy, boxSize, invBoxSize);
        threads.syncThreads();
    }

    // Sum the energy derivatives.

    for (int atom = data.firstAtom; atom < data.lastAtom; atom++) {
        for (int i = 0; i < (int) dEdV.size(); i++) {
            float sum = 0.0f;
            for (int j = 0; j < (int) threadData.size(); j++)
                sum += threadData[j]->dEdV[i][atom];
            dEdV[i][atom] = sum;
        }
    }
    threads.syncThreads();

    // Apply the chain rule to evaluate forces.

    calculateChainRuleForces(data, numberOfAtoms, posq, atomParameters, forces, boxSize, invBoxSize);
}

void CpuCustomGBForce::calculateParticlePairValue(int index, ThreadData& data, int numAtoms, float* posq, RealOpenMM** atomParameters,
        bool useExclusions, const fvec4& boxSize, const fvec4& invBoxSize) {
    for (int i = 0; i < numAtoms; i++)
        values[index][i] = 0.0f;
    vector<float>& valueArray = (index == 0 ? data.value0 : values[index]);
    if (cutoff) {
        // Loop over all pairs in the neighbor list.

        while (true) {
            int blockIndex = gmx_atomic_fetch_add(reinterpret_cast<gmx_atomic_t*>(atomicCounter), 1);
            if (blockIndex >= neighborList->getNumBlocks())
                break;
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
                        calculateOnePairValue(index, first, second, data, posq, atomParameters, valueArray, boxSize, invBoxSize);
                        calculateOnePairValue(index, second, first, data, posq, atomParameters, valueArray, boxSize, invBoxSize);
                    }
                }
            }
        }
    }
    else {
        // Perform an O(N^2) loop over all atom pairs.

        while (true) {
            int i = gmx_atomic_fetch_add(reinterpret_cast<gmx_atomic_t*>(atomicCounter), 1);
            if (i >= numAtoms)
                break;
            for (int j = i+1; j < numAtoms; j++) {
                if (useExclusions && exclusions[i].find(j) != exclusions[i].end())
                    continue;
                calculateOnePairValue(index, i, j, data, posq, atomParameters, valueArray, boxSize, invBoxSize);
                calculateOnePairValue(index, j, i, data, posq, atomParameters, valueArray, boxSize, invBoxSize);
           }
        }
    }
}

void CpuCustomGBForce::calculateOnePairValue(int index, int atom1, int atom2, ThreadData& data, float* posq, RealOpenMM** atomParameters,
        vector<float>& valueArray, const fvec4& boxSize, const fvec4& invBoxSize) {
    fvec4 deltaR;
    fvec4 pos1(posq+4*atom1);
    fvec4 pos2(posq+4*atom2);
    float r2;
    getDeltaR(pos2, pos1, deltaR, r2, periodic, boxSize, invBoxSize);
    if (cutoff && r2 >= cutoffDistance2)
        return;
    float r = sqrtf(r2);
    for (int i = 0; i < (int) paramNames.size(); i++) {
        data.expressionSet.setVariable(data.particleParamIndex[i*2], atomParameters[atom1][i]);
        data.expressionSet.setVariable(data.particleParamIndex[i*2+1], atomParameters[atom2][i]);
    }
    data.expressionSet.setVariable(data.rindex, r);
    for (int i = 0; i < index; i++) {
        data.expressionSet.setVariable(data.particleValueIndex[i*2], values[i][atom1]);
        data.expressionSet.setVariable(data.particleValueIndex[i*2+1], values[i][atom2]);
    }
    valueArray[atom1] += (float) data.valueExpressions[index].evaluate();
}

void CpuCustomGBForce::calculateSingleParticleEnergyTerm(int index, ThreadData& data, int numAtoms, float* posq,
        RealOpenMM** atomParameters, float* forces, double& totalEnergy) {
    for (int i = data.firstAtom; i < data.lastAtom; i++) {
        data.expressionSet.setVariable(data.xindex, posq[4*i]);
        data.expressionSet.setVariable(data.yindex, posq[4*i+1]);
        data.expressionSet.setVariable(data.zindex, posq[4*i+2]);
        for (int j = 0; j < (int) paramNames.size(); j++)
            data.expressionSet.setVariable(data.paramIndex[j], atomParameters[i][j]);
        for (int j = 0; j < (int) valueNames.size(); j++)
            data.expressionSet.setVariable(data.valueIndex[j], values[j][i]);
        if (includeEnergy)
            totalEnergy += (float) data.energyExpressions[index].evaluate();
        for (int j = 0; j < (int) valueNames.size(); j++)
            data.dEdV[j][i] += (float) data.energyDerivExpressions[index][j].evaluate();
        forces[4*i+0] -= (float) data.energyGradientExpressions[index][0].evaluate();
        forces[4*i+1] -= (float) data.energyGradientExpressions[index][1].evaluate();
        forces[4*i+2] -= (float) data.energyGradientExpressions[index][2].evaluate();
    }
}

void CpuCustomGBForce::calculateParticlePairEnergyTerm(int index, ThreadData& data, int numAtoms, float* posq, RealOpenMM** atomParameters,
        bool useExclusions, float* forces, double& totalEnergy, const fvec4& boxSize, const fvec4& invBoxSize) {
    if (cutoff) {
        // Loop over all pairs in the neighbor list.

        while (true) {
            int blockIndex = gmx_atomic_fetch_add(reinterpret_cast<gmx_atomic_t*>(atomicCounter), 1);
            if (blockIndex >= neighborList->getNumBlocks())
                break;
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
                        calculateOnePairEnergyTerm(index, first, second, data, posq, atomParameters, forces, totalEnergy, boxSize, invBoxSize);
                    }
                }
            }
        }
    }
    else {
        // Perform an O(N^2) loop over all atom pairs.

        while (true) {
            int i = gmx_atomic_fetch_add(reinterpret_cast<gmx_atomic_t*>(atomicCounter), 1);
            if (i >= numAtoms)
                break;
            for (int j = i+1; j < numAtoms; j++) {
                if (useExclusions && exclusions[i].find(j) != exclusions[i].end())
                    continue;
                calculateOnePairEnergyTerm(index, i, j, data, posq, atomParameters, forces, totalEnergy, boxSize, invBoxSize);
           }
        }
    }
}

void CpuCustomGBForce::calculateOnePairEnergyTerm(int index, int atom1, int atom2, ThreadData& data, float* posq, RealOpenMM** atomParameters,
        float* forces, double& totalEnergy, const fvec4& boxSize, const fvec4& invBoxSize) {
    // Compute the displacement.

    fvec4 deltaR;
    fvec4 pos1(posq+4*atom1);
    fvec4 pos2(posq+4*atom2);
    float r2;
    getDeltaR(pos2, pos1, deltaR, r2, periodic, boxSize, invBoxSize);
    if (cutoff && r2 >= cutoffDistance2)
        return;
    float r = sqrtf(r2);

    // Record variables for evaluating expressions.

    for (int i = 0; i < (int) paramNames.size(); i++) {
        data.expressionSet.setVariable(data.particleParamIndex[i*2], atomParameters[atom1][i]);
        data.expressionSet.setVariable(data.particleParamIndex[i*2+1], atomParameters[atom2][i]);
    }
    data.expressionSet.setVariable(data.rindex, r);
    for (int i = 0; i < (int) valueNames.size(); i++) {
        data.expressionSet.setVariable(data.particleValueIndex[i*2], values[i][atom1]);
        data.expressionSet.setVariable(data.particleValueIndex[i*2+1], values[i][atom2]);
    }

    // Evaluate the energy and its derivatives.

    if (includeEnergy)
        totalEnergy += (float) data.energyExpressions[index].evaluate();
    float dEdR = (float) data.energyDerivExpressions[index][0].evaluate();
    dEdR *= 1/r;
    fvec4 result = deltaR*dEdR;
    (fvec4(forces+4*atom1)-result).store(forces+4*atom1);
    (fvec4(forces+4*atom2)+result).store(forces+4*atom2);
    for (int i = 0; i < (int) valueNames.size(); i++) {
        data.dEdV[i][atom1] += (float) data.energyDerivExpressions[index][2*i+1].evaluate();
        data.dEdV[i][atom2] += (float) data.energyDerivExpressions[index][2*i+2].evaluate();
    }
}

void CpuCustomGBForce::calculateChainRuleForces(ThreadData& data, int numAtoms, float* posq, RealOpenMM** atomParameters,
        float* forces, const fvec4& boxSize, const fvec4& invBoxSize) {
    if (cutoff) {
        // Loop over all pairs in the neighbor list.

        while (true) {
            int blockIndex = gmx_atomic_fetch_add(reinterpret_cast<gmx_atomic_t*>(atomicCounter), 1);
            if (blockIndex >= neighborList->getNumBlocks())
                break;
            const int* blockAtom = &neighborList->getSortedAtoms()[4*blockIndex];
            const vector<int>& neighbors = neighborList->getBlockNeighbors(blockIndex);
            const vector<char>& blockExclusions = neighborList->getBlockExclusions(blockIndex);
            for (int i = 0; i < (int) neighbors.size(); i++) {
                int first = neighbors[i];
                for (int k = 0; k < 4; k++) {
                    if ((blockExclusions[i] & (1<<k)) == 0) {
                        int second = blockAtom[k];
                        bool isExcluded = (exclusions[first].find(second) != exclusions[first].end());
                        calculateOnePairChainRule(first, second, data, posq, atomParameters, forces, isExcluded, boxSize, invBoxSize);
                        calculateOnePairChainRule(second, first, data, posq, atomParameters, forces, isExcluded, boxSize, invBoxSize);
                    }
                }
            }
        }
    }
    else {
        // Perform an O(N^2) loop over all atom pairs.

        while (true) {
            int i = gmx_atomic_fetch_add(reinterpret_cast<gmx_atomic_t*>(atomicCounter), 1);
            if (i >= numAtoms)
                break;
            for (int j = i+1; j < numAtoms; j++) {
                bool isExcluded = (exclusions[i].find(j) != exclusions[i].end());
                calculateOnePairChainRule(i, j, data, posq, atomParameters, forces, isExcluded, boxSize, invBoxSize);
                calculateOnePairChainRule(j, i, data, posq, atomParameters, forces, isExcluded, boxSize, invBoxSize);
           }
        }
    }

    // Compute chain rule terms for computed values that depend explicitly on particle coordinates.

    for (int i = data.firstAtom; i < data.lastAtom; i++) {
        data.expressionSet.setVariable(data.xindex, posq[4*i]);
        data.expressionSet.setVariable(data.yindex, posq[4*i+1]);
        data.expressionSet.setVariable(data.zindex, posq[4*i+2]);
        for (int j = 0; j < (int) paramNames.size(); j++)
            data.expressionSet.setVariable(data.paramIndex[j], atomParameters[i][j]);
        for (int j = 1; j < (int) valueNames.size(); j++) {
            data.expressionSet.setVariable(data.valueIndex[j-1], values[j-1][i]);
            data.dVdX[j] = 0.0;
            data.dVdY[j] = 0.0;
            data.dVdZ[j] = 0.0;
            for (int k = 1; k < j; k++) {
                float dVdV = (float) data.valueDerivExpressions[j][k].evaluate();
                data.dVdX[j] += dVdV*data.dVdX[k];
                data.dVdY[j] += dVdV*data.dVdY[k];
                data.dVdZ[j] += dVdV*data.dVdZ[k];
            }
            data.dVdX[j] += (float) data.valueGradientExpressions[j][0].evaluate();
            data.dVdY[j] += (float) data.valueGradientExpressions[j][1].evaluate();
            data.dVdZ[j] += (float) data.valueGradientExpressions[j][2].evaluate();
            forces[4*i+0] -= dEdV[j][i]*data.dVdX[j];
            forces[4*i+1] -= dEdV[j][i]*data.dVdY[j];
            forces[4*i+2] -= dEdV[j][i]*data.dVdZ[j];
        }
    }
}

void CpuCustomGBForce::calculateOnePairChainRule(int atom1, int atom2, ThreadData& data, float* posq, RealOpenMM** atomParameters,
        float* forces, bool isExcluded, const fvec4& boxSize, const fvec4& invBoxSize) {
    // Compute the displacement.

    fvec4 deltaR;
    fvec4 pos1(posq+4*atom1);
    fvec4 pos2(posq+4*atom2);
    float r2;
    getDeltaR(pos2, pos1, deltaR, r2, periodic, boxSize, invBoxSize);
    if (cutoff && r2 >= cutoffDistance2)
        return;
    float r = sqrtf(r2);

    // Record variables for evaluating expressions.

    for (int i = 0; i < (int) paramNames.size(); i++) {
        data.expressionSet.setVariable(data.particleParamIndex[i*2], atomParameters[atom1][i]);
        data.expressionSet.setVariable(data.particleParamIndex[i*2+1], atomParameters[atom2][i]);
        data.expressionSet.setVariable(data.paramIndex[i], atomParameters[atom1][i]);
    }
    data.expressionSet.setVariable(data.valueIndex[0], values[0][atom1]);
    data.expressionSet.setVariable(data.xindex, posq[4*atom1]);
    data.expressionSet.setVariable(data.yindex, posq[4*atom1+1]);
    data.expressionSet.setVariable(data.zindex, posq[4*atom1+2]);
    data.expressionSet.setVariable(data.rindex, r);
    data.expressionSet.setVariable(data.particleValueIndex[0], values[0][atom1]);
    data.expressionSet.setVariable(data.particleValueIndex[1], values[0][atom2]);

    // Evaluate the derivative of each parameter with respect to position and apply forces.

    float rinv = 1/r;
    deltaR *= rinv;
    fvec4 f1(0.0f), f2(0.0f);
    if (!isExcluded || valueTypes[0] != CustomGBForce::ParticlePair) {
        data.dVdR1[0] = (float) data.valueDerivExpressions[0][0].evaluate();
        data.dVdR2[0] = -data.dVdR1[0];
        f1 -= deltaR*(dEdV[0][atom1]*data.dVdR1[0]);
        f2 -= deltaR*(dEdV[0][atom1]*data.dVdR2[0]);
    }
    for (int i = 1; i < (int) valueNames.size(); i++) {
        data.expressionSet.setVariable(data.valueIndex[i], values[i][atom1]);
        data.dVdR1[i] = 0.0;
        data.dVdR2[i] = 0.0;
        for (int j = 0; j < i; j++) {
            float dVdV = (float) data.valueDerivExpressions[i][j].evaluate();
            data.dVdR1[i] += dVdV*data.dVdR1[j];
            data.dVdR2[i] += dVdV*data.dVdR2[j];
        }
        f1 -= deltaR*(dEdV[i][atom1]*data.dVdR1[i]);
        f2 -= deltaR*(dEdV[i][atom1]*data.dVdR2[i]);
    }
    (fvec4(forces+4*atom1)+f1).store(forces+4*atom1);
    (fvec4(forces+4*atom2)+f2).store(forces+4*atom2);
}

void CpuCustomGBForce::getDeltaR(const fvec4& posI, const fvec4& posJ, fvec4& deltaR, float& r2, bool periodic, const fvec4& boxSize, const fvec4& invBoxSize) const {
    deltaR = posJ-posI;
    if (periodic) {
        fvec4 base = round(deltaR*invBoxSize)*boxSize;
        deltaR = deltaR-base;
    }
    r2 = dot3(deltaR, deltaR);
}
