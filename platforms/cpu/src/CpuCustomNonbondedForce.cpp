
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
#include "CpuCustomNonbondedForce.h"
#include "gmx_atomic.h"

using namespace OpenMM;
using namespace std;

class CpuCustomNonbondedForce::ComputeForceTask : public ThreadPool::Task {
public:
    ComputeForceTask(CpuCustomNonbondedForce& owner) : owner(owner) {
    }
    void execute(ThreadPool& threads, int threadIndex) {
        owner.threadComputeForce(threads, threadIndex);
    }
    CpuCustomNonbondedForce& owner;
};

CpuCustomNonbondedForce::ThreadData::ThreadData(const Lepton::CompiledExpression& energyExpression, const Lepton::CompiledExpression& forceExpression, const vector<string>& parameterNames) :
            energyExpression(energyExpression), forceExpression(forceExpression) {
    energyR = ReferenceForce::getVariablePointer(this->energyExpression, "r");
    forceR = ReferenceForce::getVariablePointer(this->forceExpression, "r");
    for (int i = 0; i < (int) parameterNames.size(); i++) {
        for (int j = 1; j < 3; j++) {
            stringstream name;
            name << parameterNames[i] << j;
            energyParticleParams.push_back(ReferenceForce::getVariablePointer(this->energyExpression, name.str()));
            forceParticleParams.push_back(ReferenceForce::getVariablePointer(this->forceExpression, name.str()));
        }
    }
}

CpuCustomNonbondedForce::CpuCustomNonbondedForce(const Lepton::CompiledExpression& energyExpression,
            const Lepton::CompiledExpression& forceExpression, const vector<string>& parameterNames, const vector<set<int> >& exclusions,ThreadPool& threads) :
            cutoff(false), useSwitch(false), periodic(false), paramNames(parameterNames), exclusions(exclusions), threads(threads) {
    for (int i = 0; i < threads.getNumThreads(); i++)
        threadData.push_back(new ThreadData(energyExpression, forceExpression, parameterNames));
}

CpuCustomNonbondedForce::~CpuCustomNonbondedForce() {
    for (int i = 0; i < (int) threadData.size(); i++)
        delete threadData[i];
}

void CpuCustomNonbondedForce::setUseCutoff(RealOpenMM distance, const CpuNeighborList& neighbors) {
    cutoff = true;
    cutoffDistance = distance;
    neighborList = &neighbors;
  }

void CpuCustomNonbondedForce::setInteractionGroups(const vector<pair<set<int>, set<int> > >& groups) {
    for (int group = 0; group < (int) groups.size(); group++) {
        const set<int>& set1 = groups[group].first;
        const set<int>& set2 = groups[group].second;
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

void CpuCustomNonbondedForce::setUseSwitchingFunction(RealOpenMM distance) {
    useSwitch = true;
    switchingDistance = distance;
}

void CpuCustomNonbondedForce::setPeriodic(RealVec& boxSize) {
    assert(cutoff);
    assert(boxSize[0] >= 2.0*cutoffDistance);
    assert(boxSize[1] >= 2.0*cutoffDistance);
    assert(boxSize[2] >= 2.0*cutoffDistance);
    periodic = true;
    periodicBoxSize[0] = boxSize[0];
    periodicBoxSize[1] = boxSize[1];
    periodicBoxSize[2] = boxSize[2];
}


void CpuCustomNonbondedForce::calculatePairIxn(int numberOfAtoms, float* posq, vector<RealVec>& atomCoordinates, RealOpenMM** atomParameters,
                                             RealOpenMM* fixedParameters, const map<string, double>& globalParameters,
                                             vector<AlignedArray<float> >& threadForce, bool includeForce, bool includeEnergy, double& totalEnergy) {
    // Record the parameters for the threads.
    
    this->numberOfAtoms = numberOfAtoms;
    this->posq = posq;
    this->atomCoordinates = &atomCoordinates[0];
    this->atomParameters = atomParameters;
    this->globalParameters = &globalParameters;
    this->threadForce = &threadForce;
    this->includeForce = includeForce;
    this->includeEnergy = includeEnergy;
    threadEnergy.resize(threads.getNumThreads());
    gmx_atomic_t counter;
    gmx_atomic_set(&counter, 0);
    this->atomicCounter = &counter;
    
    // Signal the threads to start running and wait for them to finish.
    
    ComputeForceTask task(*this);
    threads.execute(task);
    threads.waitForThreads();
    
    // Combine the energies from all the threads.
    
    if (includeEnergy) {
        int numThreads = threads.getNumThreads();
        for (int i = 0; i < numThreads; i++)
            totalEnergy += threadEnergy[i];
    }
}

void CpuCustomNonbondedForce::threadComputeForce(ThreadPool& threads, int threadIndex) {
    // Compute this thread's subset of interactions.

    int numThreads = threads.getNumThreads();
    threadEnergy[threadIndex] = 0;
    double& energy = threadEnergy[threadIndex];
    float* forces = &(*threadForce)[threadIndex][0];
    ThreadData& data = *threadData[threadIndex];
    for (map<string, double>::const_iterator iter = globalParameters->begin(); iter != globalParameters->end(); ++iter) {
        ReferenceForce::setVariable(ReferenceForce::getVariablePointer(data.energyExpression, iter->first), iter->second);
        ReferenceForce::setVariable(ReferenceForce::getVariablePointer(data.forceExpression, iter->first), iter->second);
    }
    fvec4 boxSize(periodicBoxSize[0], periodicBoxSize[1], periodicBoxSize[2], 0);
    fvec4 invBoxSize((1/periodicBoxSize[0]), (1/periodicBoxSize[1]), (1/periodicBoxSize[2]), 0);
    if (groupInteractions.size() > 0) {
        // The user has specified interaction groups, so compute only the requested interactions.
        
        while (true) {
            int i = gmx_atomic_fetch_add(reinterpret_cast<gmx_atomic_t*>(atomicCounter), 1);
            if (i >= groupInteractions.size())
                break;
            int atom1 = groupInteractions[i].first;
            int atom2 = groupInteractions[i].second;
            for (int j = 0; j < (int) paramNames.size(); j++) {
                ReferenceForce::setVariable(data.energyParticleParams[j*2], atomParameters[atom1][j]);
                ReferenceForce::setVariable(data.energyParticleParams[j*2+1], atomParameters[atom2][j]);
                ReferenceForce::setVariable(data.forceParticleParams[j*2], atomParameters[atom1][j]);
                ReferenceForce::setVariable(data.forceParticleParams[j*2+1], atomParameters[atom2][j]);
            }
            calculateOneIxn(atom1, atom2, data, forces, energy, boxSize, invBoxSize);
        }
    }
    else if (cutoff) {
        // We are using a cutoff, so get the interactions from the neighbor list.

        while (true) {
            int blockIndex = gmx_atomic_fetch_add(reinterpret_cast<gmx_atomic_t*>(atomicCounter), 1);
            if (blockIndex >= neighborList->getNumBlocks())
                break;
            const int* blockAtom = &neighborList->getSortedAtoms()[4*blockIndex];
            const vector<int>& neighbors = neighborList->getBlockNeighbors(blockIndex);
            const vector<char>& exclusions = neighborList->getBlockExclusions(blockIndex);
            for (int i = 0; i < (int) neighbors.size(); i++) {
                int first = neighbors[i];
                for (int j = 0; j < (int) paramNames.size(); j++) {
                    ReferenceForce::setVariable(data.energyParticleParams[j*2], atomParameters[first][j]);
                    ReferenceForce::setVariable(data.forceParticleParams[j*2], atomParameters[first][j]);
                }
                for (int k = 0; k < 4; k++) {
                    if ((exclusions[i] & (1<<k)) == 0) {
                        int second = blockAtom[k];
                        for (int j = 0; j < (int) paramNames.size(); j++) {
                            ReferenceForce::setVariable(data.energyParticleParams[j*2+1], atomParameters[second][j]);
                            ReferenceForce::setVariable(data.forceParticleParams[j*2+1], atomParameters[second][j]);
                        }
                        calculateOneIxn(first, second, data, forces, energy, boxSize, invBoxSize);
                    }
                }
            }
        }
    }
    else {
        // Every particle interacts with every other one.
        
        while (true) {
            int ii = gmx_atomic_fetch_add(reinterpret_cast<gmx_atomic_t*>(atomicCounter), 1);
            if (ii >= numberOfAtoms)
                break;
            for (int jj = ii+1; jj < numberOfAtoms; jj++) {
                if (exclusions[jj].find(ii) == exclusions[jj].end()) {
                    for (int j = 0; j < (int) paramNames.size(); j++) {
                        ReferenceForce::setVariable(data.energyParticleParams[j*2], atomParameters[ii][j]);
                        ReferenceForce::setVariable(data.energyParticleParams[j*2+1], atomParameters[jj][j]);
                        ReferenceForce::setVariable(data.forceParticleParams[j*2], atomParameters[ii][j]);
                        ReferenceForce::setVariable(data.forceParticleParams[j*2+1], atomParameters[jj][j]);
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

    // accumulate forces

    ReferenceForce::setVariable(data.energyR, r);
    ReferenceForce::setVariable(data.forceR, r);
    double dEdR = (includeForce ? data.forceExpression.evaluate()/r : 0.0);
    double energy = (includeEnergy ? data.energyExpression.evaluate() : 0.0);
    if (useSwitch) {
        if (r > switchingDistance) {
            RealOpenMM t = (r-switchingDistance)/(cutoffDistance-switchingDistance);
            RealOpenMM switchValue = 1+t*t*t*(-10+t*(15-t*6));
            RealOpenMM switchDeriv = t*t*(-30+t*(60-t*30))/(cutoffDistance-switchingDistance);
            dEdR = switchValue*dEdR + energy*switchDeriv/r;
            energy *= switchValue;
        }
    }
    fvec4 result = deltaR*dEdR;
    (fvec4(forces+4*ii)+result).store(forces+4*ii);
    (fvec4(forces+4*jj)-result).store(forces+4*jj);

    // accumulate energies

    totalEnergy += energy;
}

void CpuCustomNonbondedForce::getDeltaR(const fvec4& posI, const fvec4& posJ, fvec4& deltaR, float& r2, const fvec4& boxSize, const fvec4& invBoxSize) const {
    deltaR = posJ-posI;
    if (periodic) {
        fvec4 base = round(deltaR*invBoxSize)*boxSize;
        deltaR = deltaR-base;
    }
    r2 = dot3(deltaR, deltaR);
}
