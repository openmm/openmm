
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
#include "CpuCustomNonbondedForce.h"

using namespace OpenMM;
using namespace std;

CpuCustomNonbondedForce::ThreadData::ThreadData(const Lepton::CompiledExpression& energyExpression, const Lepton::CompiledExpression& forceExpression,
            const vector<string>& parameterNames, const std::vector<Lepton::CompiledExpression> energyParamDerivExpressions) :
            energyExpression(energyExpression), forceExpression(forceExpression), energyParamDerivExpressions(energyParamDerivExpressions) {
    map<string, double*> variableLocations;
    variableLocations["r"] = &r;
    particleParam.resize(2*parameterNames.size());
    for (int i = 0; i < (int) parameterNames.size(); i++) {
        for (int j = 0; j < 2; j++) {
            stringstream name;
            name << parameterNames[i] << (j+1);
            variableLocations[name.str()] = &particleParam[i*2+j];
        }
    }
    energyParamDerivs.resize(energyParamDerivExpressions.size());
    this->energyExpression.setVariableLocations(variableLocations);
    this->forceExpression.setVariableLocations(variableLocations);
    expressionSet.registerExpression(this->energyExpression);
    expressionSet.registerExpression(this->forceExpression);
    for (auto& expression : this->energyParamDerivExpressions) {
        expression.setVariableLocations(variableLocations);
        expressionSet.registerExpression(expression);
    }
}

CpuCustomNonbondedForce::CpuCustomNonbondedForce(const Lepton::CompiledExpression& energyExpression,
            const Lepton::CompiledExpression& forceExpression, const vector<string>& parameterNames, const vector<set<int> >& exclusions,
            const std::vector<Lepton::CompiledExpression> energyParamDerivExpressions, ThreadPool& threads) :
            cutoff(false), useSwitch(false), periodic(false), useInteractionGroups(false), paramNames(parameterNames), exclusions(exclusions), threads(threads) {
    for (int i = 0; i < threads.getNumThreads(); i++)
        threadData.push_back(new ThreadData(energyExpression, forceExpression, parameterNames, energyParamDerivExpressions));
}

CpuCustomNonbondedForce::~CpuCustomNonbondedForce() {
    for (auto data : threadData)
        delete data;
}

void CpuCustomNonbondedForce::setUseCutoff(double distance, const CpuNeighborList& neighbors) {
    cutoff = true;
    cutoffDistance = distance;
    neighborList = &neighbors;
  }

void CpuCustomNonbondedForce::setInteractionGroups(const vector<pair<set<int>, set<int> > >& groups) {
    useInteractionGroups = true;
    for (auto& group : groups) {
        const set<int>& set1 = group.first;
        const set<int>& set2 = group.second;
        for (set<int>::const_iterator atom1 = set1.begin(); atom1 != set1.end(); ++atom1) {
            for (set<int>::const_iterator atom2 = set2.begin(); atom2 != set2.end(); ++atom2) {
                if (*atom1 == *atom2 || exclusions[*atom1].find(*atom2) != exclusions[*atom1].end())
                    continue; // This is an excluded interaction.
                if (*atom1 > *atom2 && set1.find(*atom2) != set1.end() && set2.find(*atom1) != set2.end())
                    continue; // Both atoms are in both sets, so skip duplicate interactions.
                groupInteractions.push_back(make_pair(*atom1, *atom2));
            }
        }
    }
}

void CpuCustomNonbondedForce::setUseSwitchingFunction(double distance) {
    useSwitch = true;
    switchingDistance = distance;
}

void CpuCustomNonbondedForce::setPeriodic(Vec3* periodicBoxVectors) {
    assert(cutoff);
    assert(periodicBoxVectors[0][0] >= 2.0*cutoffDistance);
    assert(periodicBoxVectors[1][1] >= 2.0*cutoffDistance);
    assert(periodicBoxVectors[2][2] >= 2.0*cutoffDistance);
    periodic = true;
    this->periodicBoxVectors[0] = periodicBoxVectors[0];
    this->periodicBoxVectors[1] = periodicBoxVectors[1];
    this->periodicBoxVectors[2] = periodicBoxVectors[2];
    recipBoxSize[0] = (float) (1.0/periodicBoxVectors[0][0]);
    recipBoxSize[1] = (float) (1.0/periodicBoxVectors[1][1]);
    recipBoxSize[2] = (float) (1.0/periodicBoxVectors[2][2]);
    periodicBoxVec4.resize(3);
    periodicBoxVec4[0] = fvec4(periodicBoxVectors[0][0], periodicBoxVectors[0][1], periodicBoxVectors[0][2], 0);
    periodicBoxVec4[1] = fvec4(periodicBoxVectors[1][0], periodicBoxVectors[1][1], periodicBoxVectors[1][2], 0);
    periodicBoxVec4[2] = fvec4(periodicBoxVectors[2][0], periodicBoxVectors[2][1], periodicBoxVectors[2][2], 0);
    triclinic = (periodicBoxVectors[0][1] != 0.0 || periodicBoxVectors[0][2] != 0.0 ||
                 periodicBoxVectors[1][0] != 0.0 || periodicBoxVectors[1][2] != 0.0 ||
                 periodicBoxVectors[2][0] != 0.0 || periodicBoxVectors[2][1] != 0.0);
}


void CpuCustomNonbondedForce::calculatePairIxn(int numberOfAtoms, float* posq, vector<Vec3>& atomCoordinates, vector<vector<double> >& atomParameters,
                                               const map<string, double>& globalParameters, vector<AlignedArray<float> >& threadForce,
                                               bool includeForce, bool includeEnergy, double& totalEnergy, double* energyParamDerivs) {
    // Record the parameters for the threads.
    
    this->numberOfAtoms = numberOfAtoms;
    this->posq = posq;
    this->atomCoordinates = &atomCoordinates[0];
    this->atomParameters = &atomParameters[0];
    this->globalParameters = &globalParameters;
    this->threadForce = &threadForce;
    this->includeForce = includeForce;
    this->includeEnergy = includeEnergy;
    threadEnergy.resize(threads.getNumThreads());
    atomicCounter = 0;
    
    // Signal the threads to start running and wait for them to finish.
    
    threads.execute([&] (ThreadPool& threads, int threadIndex) { threadComputeForce(threads, threadIndex); });
    threads.waitForThreads();
    
    // Combine the energies from all the threads.
    
    int numThreads = threads.getNumThreads();
    if (includeEnergy) {
        for (int i = 0; i < numThreads; i++)
            totalEnergy += threadEnergy[i];
    }

    // Combine the energy derivatives from all threads.
    
    int numDerivs = threadData[0]->energyParamDerivs.size();
    for (int i = 0; i < numThreads; i++)
        for (int j = 0; j < numDerivs; j++)
            energyParamDerivs[j] += threadData[i]->energyParamDerivs[j];
}

void CpuCustomNonbondedForce::threadComputeForce(ThreadPool& threads, int threadIndex) {
    // Compute this thread's subset of interactions.

    int numThreads = threads.getNumThreads();
    threadEnergy[threadIndex] = 0;
    double& energy = threadEnergy[threadIndex];
    float* forces = &(*threadForce)[threadIndex][0];
    ThreadData& data = *threadData[threadIndex];
    for (auto& param : *globalParameters)
        data.expressionSet.setVariable(data.expressionSet.getVariableIndex(param.first), param.second);
    for (auto& deriv : data.energyParamDerivs)
        deriv = 0.0;
    fvec4 boxSize(periodicBoxVectors[0][0], periodicBoxVectors[1][1], periodicBoxVectors[2][2], 0);
    fvec4 invBoxSize(recipBoxSize[0], recipBoxSize[1], recipBoxSize[2], 0);
    if (useInteractionGroups) {
        // The user has specified interaction groups, so compute only the requested interactions.
        
        int start = threadIndex*groupInteractions.size()/numThreads;
        int end = (threadIndex+1)*groupInteractions.size()/numThreads;
        for (int i = start; i < end; i++) {
            int atom1 = groupInteractions[i].first;
            int atom2 = groupInteractions[i].second;
            for (int j = 0; j < (int) paramNames.size(); j++) {
                data.particleParam[j*2] = atomParameters[atom1][j];
                data.particleParam[j*2+1] = atomParameters[atom2][j];
            }
            calculateOneIxn(atom1, atom2, data, forces, energy, boxSize, invBoxSize);
        }
    }
    else if (cutoff) {
        // We are using a cutoff, so get the interactions from the neighbor list.

        while (true) {
            int blockIndex = atomicCounter++;
            if (blockIndex >= neighborList->getNumBlocks())
                break;
            const int blockSize = neighborList->getBlockSize();
            const int* blockAtom = &neighborList->getSortedAtoms()[blockSize*blockIndex];
            const vector<int>& neighbors = neighborList->getBlockNeighbors(blockIndex);
            const vector<char>& exclusions = neighborList->getBlockExclusions(blockIndex);
            for (int i = 0; i < (int) neighbors.size(); i++) {
                int first = neighbors[i];
                for (int j = 0; j < (int) paramNames.size(); j++)
                    data.particleParam[j*2] = atomParameters[first][j];
                for (int k = 0; k < blockSize; k++) {
                    if ((exclusions[i] & (1<<k)) == 0) {
                        int second = blockAtom[k];
                        for (int j = 0; j < (int) paramNames.size(); j++)
                            data.particleParam[j*2+1] = atomParameters[second][j];
                        calculateOneIxn(first, second, data, forces, energy, boxSize, invBoxSize);
                    }
                }
            }
        }
    }
    else {
        // Every particle interacts with every other one.
        
        while (true) {
            int ii = atomicCounter++;
            if (ii >= numberOfAtoms)
                break;
            for (int jj = ii+1; jj < numberOfAtoms; jj++) {
                if (exclusions[jj].find(ii) == exclusions[jj].end()) {
                    for (int j = 0; j < (int) paramNames.size(); j++) {
                        data.particleParam[j*2] = atomParameters[ii][j];
                        data.particleParam[j*2+1] = atomParameters[jj][j];
                    }
                    calculateOneIxn(ii, jj, data, forces, energy, boxSize, invBoxSize);
                }
            }
        }
    }
}

void CpuCustomNonbondedForce::calculateOneIxn(int ii, int jj, ThreadData& data, 
        float* forces, double& totalEnergy, const fvec4& boxSize, const fvec4& invBoxSize) {
    // Get deltaR, R2, and R between 2 atoms

    fvec4 deltaR;
    fvec4 posI(posq+4*ii);
    fvec4 posJ(posq+4*jj);
    float r2;
    getDeltaR(posI, posJ, deltaR, r2, boxSize, invBoxSize);
    if (cutoff && r2 >= cutoffDistance*cutoffDistance)
        return;
    float r = sqrtf(r2);
    data.r = r;

    // accumulate forces

    double dEdR = (includeForce ? data.forceExpression.evaluate()/r : 0.0);
    double energy = 0.0;
    if (includeEnergy || (useSwitch && r > switchingDistance))
        energy = data.energyExpression.evaluate();
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
    fvec4 result = deltaR*dEdR;
    (fvec4(forces+4*ii)+result).store(forces+4*ii);
    (fvec4(forces+4*jj)-result).store(forces+4*jj);

    // accumulate energies

    totalEnergy += energy;
    
    // Accumulate energy derivatives.

    for (int i = 0; i < data.energyParamDerivExpressions.size(); i++)
        data.energyParamDerivs[i] += switchValue*data.energyParamDerivExpressions[i].evaluate();
}

void CpuCustomNonbondedForce::getDeltaR(const fvec4& posI, const fvec4& posJ, fvec4& deltaR, float& r2, const fvec4& boxSize, const fvec4& invBoxSize) const {
    deltaR = posJ-posI;
    if (periodic) {
        if (triclinic) {
            deltaR -= periodicBoxVec4[2]*floorf(deltaR[2]*recipBoxSize[2]+0.5f);
            deltaR -= periodicBoxVec4[1]*floorf(deltaR[1]*recipBoxSize[1]+0.5f);
            deltaR -= periodicBoxVec4[0]*floorf(deltaR[0]*recipBoxSize[0]+0.5f);
        }
        else {
            fvec4 base = round(deltaR*invBoxSize)*boxSize;
            deltaR = deltaR-base;
        }
    }
    r2 = dot3(deltaR, deltaR);
}
