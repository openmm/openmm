
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
#include "CpuCustomGBForce.h"

using namespace OpenMM;
using namespace std;

CpuCustomGBForce::ThreadData::ThreadData(int numAtoms, int numThreads, int threadIndex,
                      const vector<Lepton::CompiledExpression>& valueExpressions,
                      const vector<vector<Lepton::CompiledExpression> >& valueDerivExpressions,
                      const vector<vector<Lepton::CompiledExpression> >& valueGradientExpressions,
                      const vector<vector<Lepton::CompiledExpression> >& valueParamDerivExpressions,
                      const vector<string>& valueNames,
                      const vector<Lepton::CompiledExpression>& energyExpressions,
                      const vector<vector<Lepton::CompiledExpression> >& energyDerivExpressions,
                      const vector<vector<Lepton::CompiledExpression> >& energyGradientExpressions,
                      const vector<vector<Lepton::CompiledExpression> >& energyParamDerivExpressions,
                      const vector<string>& parameterNames) :
            valueExpressions(valueExpressions), valueDerivExpressions(valueDerivExpressions), valueGradientExpressions(valueGradientExpressions),
            valueParamDerivExpressions(valueParamDerivExpressions), energyExpressions(energyExpressions), energyDerivExpressions(energyDerivExpressions),
            energyGradientExpressions(energyGradientExpressions), energyParamDerivExpressions(energyParamDerivExpressions) {
    firstAtom = (threadIndex*(long long) numAtoms)/numThreads;
    lastAtom = ((threadIndex+1)*(long long) numAtoms)/numThreads;
    map<string, double*> variableLocations;
    variableLocations["x"] = &x;
    variableLocations["y"] = &y;
    variableLocations["z"] = &z;
    variableLocations["r"] = &r;
    param.resize(parameterNames.size());
    particleParam.resize(parameterNames.size()*2);
    for (int i = 0; i < (int) parameterNames.size(); i++) {
        variableLocations[parameterNames[i]] = &param[i];
        for (int j = 0; j < 2; j++) {
            stringstream name;
            name << parameterNames[i] << (j+1);
            variableLocations[name.str()] = &particleParam[2*i+j];
        }
    }
    value.resize(valueNames.size());
    particleValue.resize(valueNames.size()*2);
    for (int i = 0; i < (int) valueNames.size(); i++) {
        variableLocations[valueNames[i]] = &value[i];
        for (int j = 0; j < 2; j++) {
            stringstream name;
            name << valueNames[i] << (j+1);
            variableLocations[name.str()] = &particleValue[2*i+j];
        }
    }
    for (auto& expression : this->valueExpressions) {
        expression.setVariableLocations(variableLocations);
        expressionSet.registerExpression(expression);
    }
    for (auto& expressions : this->valueDerivExpressions)
        for (auto& expression : expressions) {
            expression.setVariableLocations(variableLocations);
            expressionSet.registerExpression(expression);
        }
    for (auto& expressions : this->valueGradientExpressions)
        for (auto& expression : expressions) {
            expression.setVariableLocations(variableLocations);
            expressionSet.registerExpression(expression);
        }
    for (auto& expressions : this->valueParamDerivExpressions)
        for (auto& expression : expressions) {
            expression.setVariableLocations(variableLocations);
            expressionSet.registerExpression(expression);
        }
    for (auto& expression : this->energyExpressions) {
        expression.setVariableLocations(variableLocations);
        expressionSet.registerExpression(expression);
    }
    for (auto& expressions : this->energyDerivExpressions)
        for (auto& expression : expressions) {
            expression.setVariableLocations(variableLocations);
            expressionSet.registerExpression(expression);
        }
    for (auto& expressions : this->energyGradientExpressions)
        for (auto& expression : expressions) {
            expression.setVariableLocations(variableLocations);
            expressionSet.registerExpression(expression);
        }
    for (auto& expressions : this->energyParamDerivExpressions)
        for (auto& expression : expressions) {
            expression.setVariableLocations(variableLocations);
            expressionSet.registerExpression(expression);
        }
    value0.resize(numAtoms);
    dEdV.resize(valueNames.size());
    for (auto& v : dEdV)
        v.resize(numAtoms);
    dVdX.resize(valueDerivExpressions.size());
    dVdY.resize(valueDerivExpressions.size());
    dVdZ.resize(valueDerivExpressions.size());
    dVdR1.resize(valueDerivExpressions.size());
    dVdR2.resize(valueDerivExpressions.size());
    dValue0dParam.resize(valueParamDerivExpressions[0].size(), vector<float>(numAtoms));
    energyParamDerivs.resize(valueParamDerivExpressions[0].size());
}

CpuCustomGBForce::CpuCustomGBForce(int numAtoms, const std::vector<std::set<int> >& exclusions,
                     const vector<Lepton::CompiledExpression>& valueExpressions,
                     const vector<vector<Lepton::CompiledExpression> >& valueDerivExpressions,
                     const vector<vector<Lepton::CompiledExpression> >& valueGradientExpressions,
                     const vector<vector<Lepton::CompiledExpression> >& valueParamDerivExpressions,
                     const vector<string>& valueNames,
                     const vector<CustomGBForce::ComputationType>& valueTypes,
                     const vector<Lepton::CompiledExpression>& energyExpressions,
                     const vector<vector<Lepton::CompiledExpression> >& energyDerivExpressions,
                     const vector<vector<Lepton::CompiledExpression> >& energyGradientExpressions,
                     const vector<vector<Lepton::CompiledExpression> >& energyParamDerivExpressions,
                     const vector<CustomGBForce::ComputationType>& energyTypes,
                     const vector<string>& parameterNames, ThreadPool& threads) :
            exclusions(exclusions), cutoff(false), periodic(false), valueTypes(valueTypes), energyTypes(energyTypes), numValues(valueNames.size()),
            numParams(parameterNames.size()), threads(threads) {
    for (int i = 0; i < threads.getNumThreads(); i++)
        threadData.push_back(new ThreadData(numAtoms, threads.getNumThreads(), i, valueExpressions, valueDerivExpressions, valueGradientExpressions,
                valueParamDerivExpressions, valueNames, energyExpressions, energyDerivExpressions, energyGradientExpressions, energyParamDerivExpressions, parameterNames));
    values.resize(numValues);
    dEdV.resize(numValues);
    for (int i = 0; i < (int) values.size(); i++) {
        values[i].resize(numAtoms);
        dEdV[i].resize(numAtoms);
    }
    dValuedParam.resize(numValues);
    for (int i = 0; i < numValues; i++)
        dValuedParam[i].resize(valueParamDerivExpressions[0].size(), vector<float>(numAtoms));
}

CpuCustomGBForce::~CpuCustomGBForce() {
    for (auto data : threadData)
        delete data;
}

void CpuCustomGBForce::setUseCutoff(float distance, const CpuNeighborList& neighbors) {
    cutoff = true;
    cutoffDistance = distance;
    cutoffDistance2 = distance*distance;
    neighborList = &neighbors;
  }

void CpuCustomGBForce::setPeriodic(Vec3& boxSize) {
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

void CpuCustomGBForce::calculateIxn(int numberOfAtoms, float* posq, vector<vector<double> >& atomParameters,
                                           map<string, double>& globalParameters, vector<AlignedArray<float> >& threadForce,
                                           bool includeForce, bool includeEnergy, double& totalEnergy, double* energyParamDerivs) {
    // Record the parameters for the threads.
    
    this->numberOfAtoms = numberOfAtoms;
    this->posq = posq;
    this->atomParameters = &atomParameters[0];
    this->globalParameters = &globalParameters;
    this->threadForce = &threadForce;
    this->includeForce = includeForce;
    this->includeEnergy = includeEnergy;
    threadEnergy.resize(threads.getNumThreads());

    // Calculate the first computed value.

    auto task = [&] (ThreadPool& threads, int threadIndex) { threadComputeForce(threads, threadIndex); };
    atomicCounter = 0;
    threads.execute(task);
    threads.waitForThreads();

    // Sum derivatives of the first computed value with respect to global parameters.

    bool hasParamDerivs = (threadData[0]->dValue0dParam.size() > 0);
    if (hasParamDerivs) {
        threads.resumeThreads();
        threads.waitForThreads();
    }
    
    // Calculate the remaining computed values.
    
    threads.resumeThreads();
    threads.waitForThreads();
    
    // Calculate the energy terms.

    for (int i = 0; i < (int) threadData[0]->energyExpressions.size(); i++) {
        atomicCounter = 0;
        threads.execute(task);
        threads.waitForThreads();
    }

    // Sum the energy derivatives.

    threads.resumeThreads();
    threads.waitForThreads();
    
    // Apply the chain rule to evaluate forces.

    atomicCounter = 0;
    threads.resumeThreads();
    threads.waitForThreads();

    // Combine the energies from all the threads.
    
    if (includeEnergy) {
        int numThreads = threads.getNumThreads();
        for (int i = 0; i < numThreads; i++)
            totalEnergy += threadEnergy[i];
    }
    if (hasParamDerivs)
        for (int i = 0; i < threads.getNumThreads(); i++)
            for (int j = 0; j < threadData[i]->energyParamDerivs.size(); j++)
                energyParamDerivs[j] += threadData[i]->energyParamDerivs[j];
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
    for (auto& param : *globalParameters)
        data.expressionSet.setVariable(data.expressionSet.getVariableIndex(param.first), param.second);

    // Calculate the first computed value.

    for (auto& v : data.value0)
        v = 0.0f;
    for (auto& vals : data.dValue0dParam)
        for (auto& v : vals)
            v = 0.0f;
    if (valueTypes[0] == CustomGBForce::ParticlePair)
        calculateParticlePairValue(0, data, numberOfAtoms, posq, atomParameters, true, boxSize, invBoxSize);
    else
        calculateParticlePairValue(0, data, numberOfAtoms, posq, atomParameters, false, boxSize, invBoxSize);
    threads.syncThreads();
    
    // Sum derivatives of the first computed value with respect to global parameters.
    
    bool hasParamDerivs = (data.dValue0dParam.size() > 0);
    if (hasParamDerivs) {
        for (int j = 0; j < data.dValue0dParam.size(); j++)
            for (int k = data.firstAtom; k < data.lastAtom; k++) {
                float sum = 0.0f;
                for (int m = 0; m < threadData.size(); m++)
                    sum += threadData[m]->dValue0dParam[j][k];
                dValuedParam[0][j][k] = sum;
            }
        threads.syncThreads();
    }

    // Sum the first computed value and calculate the remaining ones.

    int numValues = valueTypes.size();
    for (int atom = data.firstAtom; atom < data.lastAtom; atom++) {
        float sum = 0.0f;
        for (auto& data : threadData)
            sum += data->value0[atom];
        values[0][atom] = sum;
        data.x = posq[4*atom];
        data.y = posq[4*atom+1];
        data.z = posq[4*atom+2];
        for (int j = 0; j < numParams; j++)
            data.param[j] = atomParameters[atom][j];
        for (int i = 1; i < numValues; i++) {
            data.value[i-1] = values[i-1][atom];
            values[i][atom] = (float) data.valueExpressions[i].evaluate();

            // Calculate derivatives with respect to parameters.

            if (hasParamDerivs) {
                for (int j = 0; j < data.valueParamDerivExpressions[i].size(); j++)
                    dValuedParam[i][j][atom] = data.valueParamDerivExpressions[i][j].evaluate();
                for (int j = 0; j < i; j++) {
                    float dVdV = data.valueDerivExpressions[i][j].evaluate();
                    for (int k = 0; k < data.valueParamDerivExpressions[i].size(); k++)
                        dValuedParam[i][k][atom] += dVdV*dValuedParam[j][k][atom];
                }
            }
        }
    }
    threads.syncThreads();

    // Now calculate the energy and its derivatives.

    for (auto& vals : data.dEdV)
        for (auto& v : vals)
            v = 0.0f;
    for (auto& v : data.energyParamDerivs)
        v = 0.0f;
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
            for (auto& data : threadData)
                sum += data->dEdV[i][atom];
            dEdV[i][atom] = sum;
        }
    }
    threads.syncThreads();

    // Apply the chain rule to evaluate forces.

    calculateChainRuleForces(data, numberOfAtoms, posq, atomParameters, forces, boxSize, invBoxSize);
}

void CpuCustomGBForce::calculateParticlePairValue(int index, ThreadData& data, int numAtoms, float* posq, vector<double>* atomParameters,
        bool useExclusions, const fvec4& boxSize, const fvec4& invBoxSize) {
    for (int i = 0; i < numAtoms; i++)
        values[index][i] = 0.0f;
    vector<float>& valueArray = (index == 0 ? data.value0 : values[index]);
    if (cutoff) {
        // Loop over all pairs in the neighbor list.

        while (true) {
            int blockIndex = atomicCounter++;
            if (blockIndex >= neighborList->getNumBlocks())
                break;
            const int blockSize = neighborList->getBlockSize();
            const int* blockAtom = &neighborList->getSortedAtoms()[blockSize*blockIndex];
            const vector<int>& neighbors = neighborList->getBlockNeighbors(blockIndex);
            const vector<char>& blockExclusions = neighborList->getBlockExclusions(blockIndex);
            for (int i = 0; i < (int) neighbors.size(); i++) {
                int first = neighbors[i];
                for (int k = 0; k < blockSize; k++) {
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
            int i = atomicCounter++;
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

void CpuCustomGBForce::calculateOnePairValue(int index, int atom1, int atom2, ThreadData& data, float* posq, vector<double>* atomParameters,
        vector<float>& valueArray, const fvec4& boxSize, const fvec4& invBoxSize) {
    fvec4 deltaR;
    fvec4 pos1(posq+4*atom1);
    fvec4 pos2(posq+4*atom2);
    float r2;
    getDeltaR(pos2, pos1, deltaR, r2, periodic, boxSize, invBoxSize);
    if (cutoff && r2 >= cutoffDistance2)
        return;
    data.r = sqrtf(r2);
    for (int i = 0; i < numParams; i++) {
        data.particleParam[i*2] = atomParameters[atom1][i];
        data.particleParam[i*2+1] = atomParameters[atom2][i];
    }
    for (int i = 0; i < index; i++) {
        data.particleValue[i*2] = values[i][atom1];
        data.particleValue[i*2+1] = values[i][atom2];
    }
    valueArray[atom1] += (float) data.valueExpressions[index].evaluate();
    
    // Calculate derivatives with respect to parameters.
    
    for (int i = 0; i < data.valueParamDerivExpressions[index].size(); i++)
        data.dValue0dParam[i][atom1] += data.valueParamDerivExpressions[index][i].evaluate();
}

void CpuCustomGBForce::calculateSingleParticleEnergyTerm(int index, ThreadData& data, int numAtoms, float* posq,
        vector<double>* atomParameters, float* forces, double& totalEnergy) {
    for (int i = data.firstAtom; i < data.lastAtom; i++) {
        data.x = posq[4*i];
        data.y = posq[4*i+1];
        data.z = posq[4*i+2];
        for (int j = 0; j < numParams; j++)
            data.param[j] = atomParameters[i][j];
        for (int j = 0; j < (int) values.size(); j++)
            data.value[j] = values[j][i];
        if (includeEnergy)
            totalEnergy += (float) data.energyExpressions[index].evaluate();
        for (int j = 0; j < (int) values.size(); j++)
            data.dEdV[j][i] += (float) data.energyDerivExpressions[index][j].evaluate();
        forces[4*i+0] -= (float) data.energyGradientExpressions[index][0].evaluate();
        forces[4*i+1] -= (float) data.energyGradientExpressions[index][1].evaluate();
        forces[4*i+2] -= (float) data.energyGradientExpressions[index][2].evaluate();
        
        // Compute derivatives with respect to parameters.
        
        for (int k = 0; k < data.energyParamDerivExpressions[index].size(); k++)
            data.energyParamDerivs[k] += data.energyParamDerivExpressions[index][k].evaluate();
    }
}

void CpuCustomGBForce::calculateParticlePairEnergyTerm(int index, ThreadData& data, int numAtoms, float* posq, vector<double>* atomParameters,
        bool useExclusions, float* forces, double& totalEnergy, const fvec4& boxSize, const fvec4& invBoxSize) {
    if (cutoff) {
        // Loop over all pairs in the neighbor list.

        while (true) {
            int blockIndex = atomicCounter++;
            if (blockIndex >= neighborList->getNumBlocks())
                break;
            const int blockSize = neighborList->getBlockSize();
            const int* blockAtom = &neighborList->getSortedAtoms()[blockSize*blockIndex];
            const vector<int>& neighbors = neighborList->getBlockNeighbors(blockIndex);
            const vector<char>& blockExclusions = neighborList->getBlockExclusions(blockIndex);
            for (int i = 0; i < (int) neighbors.size(); i++) {
                int first = neighbors[i];
                for (int k = 0; k < blockSize; k++) {
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
            int i = atomicCounter++;
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

void CpuCustomGBForce::calculateOnePairEnergyTerm(int index, int atom1, int atom2, ThreadData& data, float* posq, vector<double>* atomParameters,
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
    data.r = r;

    // Record variables for evaluating expressions.

    for (int i = 0; i < numParams; i++) {
        data.particleParam[i*2] = atomParameters[atom1][i];
        data.particleParam[i*2+1] = atomParameters[atom2][i];
    }
    for (int i = 0; i < (int) values.size(); i++) {
        data.particleValue[i*2] = values[i][atom1];
        data.particleValue[i*2+1] = values[i][atom2];
    }

    // Evaluate the energy and its derivatives.

    if (includeEnergy)
        totalEnergy += (float) data.energyExpressions[index].evaluate();
    float dEdR = (float) data.energyDerivExpressions[index][0].evaluate();
    dEdR *= 1/r;
    fvec4 result = deltaR*dEdR;
    (fvec4(forces+4*atom1)-result).store(forces+4*atom1);
    (fvec4(forces+4*atom2)+result).store(forces+4*atom2);
    for (int i = 0; i < (int) values.size(); i++) {
        data.dEdV[i][atom1] += (float) data.energyDerivExpressions[index][2*i+1].evaluate();
        data.dEdV[i][atom2] += (float) data.energyDerivExpressions[index][2*i+2].evaluate();
    }
        
    // Compute derivatives with respect to parameters.

    for (int i = 0; i < data.energyParamDerivExpressions[index].size(); i++)
        data.energyParamDerivs[i] += data.energyParamDerivExpressions[index][i].evaluate();
}

void CpuCustomGBForce::calculateChainRuleForces(ThreadData& data, int numAtoms, float* posq, vector<double>* atomParameters,
        float* forces, const fvec4& boxSize, const fvec4& invBoxSize) {
    if (cutoff) {
        // Loop over all pairs in the neighbor list.

        while (true) {
            int blockIndex = atomicCounter++;
            if (blockIndex >= neighborList->getNumBlocks())
                break;
            const int blockSize = neighborList->getBlockSize();
            const int* blockAtom = &neighborList->getSortedAtoms()[blockSize*blockIndex];
            const vector<int>& neighbors = neighborList->getBlockNeighbors(blockIndex);
            const vector<char>& blockExclusions = neighborList->getBlockExclusions(blockIndex);
            for (int i = 0; i < (int) neighbors.size(); i++) {
                int first = neighbors[i];
                for (int k = 0; k < blockSize; k++) {
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
            int i = atomicCounter++;
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
        data.x = posq[4*i];
        data.y = posq[4*i+1];
        data.z = posq[4*i+2];
        for (int j = 0; j < numParams; j++)
            data.param[j] = atomParameters[i][j];
        for (int j = 1; j < (int) values.size(); j++) {
            data.value[j-1] = values[j-1][i];
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
        
    // Compute chain rule terms for derivatives with respect to parameters.

    for (int i = data.firstAtom; i < data.lastAtom; i++)
        for (int j = 0; j < data.value.size(); j++)
            for (int k = 0; k < dValuedParam[j].size(); k++)
                data.energyParamDerivs[k] += dEdV[j][i]*dValuedParam[j][k][i];
}

void CpuCustomGBForce::calculateOnePairChainRule(int atom1, int atom2, ThreadData& data, float* posq, vector<double>* atomParameters,
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
    data.r = r;

    // Record variables for evaluating expressions.

    for (int i = 0; i < numParams; i++) {
        data.particleParam[i*2] = atomParameters[atom1][i];
        data.particleParam[i*2+1] = atomParameters[atom2][i];
        data.param[i] = atomParameters[atom1][i];
    }
    data.value[0] = values[0][atom1];
    data.x = posq[4*atom1];
    data.y = posq[4*atom1+1];
    data.z = posq[4*atom1+2];
    data.particleValue[0] = values[0][atom1];
    data.particleValue[1] = values[0][atom2];

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
    for (int i = 1; i < (int) values.size(); i++) {
        data.value[i] = values[i][atom1];
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
