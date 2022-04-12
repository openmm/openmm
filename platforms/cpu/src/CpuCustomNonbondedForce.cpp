
/* Portions copyright (c) 2009-2022 Stanford University and Simbios.
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

#include "SimTKOpenMMUtilities.h"
#include "ReferenceForce.h"
#include "CpuCustomNonbondedForce.h"

using namespace OpenMM;
using namespace Lepton;
using namespace std;

int getVecBlockSize();

CpuCustomNonbondedForce::ThreadData::ThreadData(const CompiledExpression& energyExpression, const CompiledVectorExpression& energyVecExpression,
            const CompiledExpression& forceExpression, const CompiledVectorExpression& forceVecExpression,
            const vector<string>& parameterNames, const std::vector<CompiledExpression> energyParamDerivExpressions,
            const vector<string>& computedValueNames, const vector<CompiledExpression> computedValueExpressions,
            vector<vector<double> >& atomComputedValues) :
            energyExpression(energyExpression), energyVecExpression(energyVecExpression), forceExpression(forceExpression),
            forceVecExpression(forceVecExpression), energyParamDerivExpressions(energyParamDerivExpressions),
            computedValueExpressions(computedValueExpressions), atomComputedValues(atomComputedValues) {
    // Prepare for passing variables to expressions.

    map<string, double*> variableLocations;
    variableLocations["r"] = &r;
    particleParam.resize(2*parameterNames.size());
    computedValues.resize(2*computedValueNames.size());
    for (int i = 0; i < parameterNames.size(); i++) {
        variableLocations[parameterNames[i]+"1"] = &particleParam[i*2];
        variableLocations[parameterNames[i]+"2"] = &particleParam[i*2+1];
    }
    for (int i = 0; i < computedValueNames.size(); i++) {
        variableLocations[computedValueNames[i]+"1"] = &computedValues[i*2];
        variableLocations[computedValueNames[i]+"2"] = &computedValues[i*2+1];
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

    // Prepare for passing variables to vectorized expressions.

    map<string, float*> vecVariableLocations;
    rvec.resize(blockSize);
    vecParticle1Params.resize(blockSize*parameterNames.size());
    vecParticle2Params.resize(blockSize*parameterNames.size());
    vecParticle1Values.resize(blockSize*computedValueNames.size());
    vecParticle2Values.resize(blockSize*computedValueNames.size());
    vecVariableLocations["r"] = rvec.data();
    for (int i = 0; i < parameterNames.size(); i++) {
        vecVariableLocations[parameterNames[i]+"1"] = &vecParticle1Params[i*blockSize];
        vecVariableLocations[parameterNames[i]+"2"] = &vecParticle2Params[i*blockSize];
    }
    for (int i = 0; i < computedValueNames.size(); i++) {
        vecVariableLocations[computedValueNames[i]+"1"] = &vecParticle1Values[i*blockSize];
        vecVariableLocations[computedValueNames[i]+"2"] = &vecParticle2Values[i*blockSize];
    }
    this->energyVecExpression.setVariableLocations(vecVariableLocations);
    this->forceVecExpression.setVariableLocations(vecVariableLocations);

    // Prepare for passing variables to the computed value expressions.

    map<string, double*> valueVariableLocations;
    for (int i = 0; i < parameterNames.size(); i++)
        valueVariableLocations[parameterNames[i]] = &particleParam[i];
    for (auto& expression : this->computedValueExpressions) {
        expression.setVariableLocations(valueVariableLocations);
        expressionSet.registerExpression(expression);
    }
}

CpuCustomNonbondedForce::CpuCustomNonbondedForce(const ParsedExpression& energyExpression,
            const ParsedExpression& forceExpression, const vector<string>& parameterNames, const vector<set<int> >& exclusions,
            const vector<ParsedExpression> energyParamDerivExpressions, const vector<string>& computedValueNames,
            const vector<ParsedExpression> computedValueExpressions, ThreadPool& threads) :
            cutoff(false), useSwitch(false), periodic(false), useInteractionGroups(false), paramNames(parameterNames), exclusions(exclusions),
            computedValueNames(computedValueNames), threads(threads) {
    CompiledExpression compiledEnergyExpression = energyExpression.createCompiledExpression();
    CompiledExpression compiledForceExpression = forceExpression.createCompiledExpression();
    CompiledVectorExpression energyVecExpression = energyExpression.createCompiledVectorExpression(blockSize);
    CompiledVectorExpression forceVecExpression = forceExpression.createCompiledVectorExpression(blockSize);
    vector<CompiledExpression> compiledDerivExpressions, compiledValueExpressions;
    for (auto& exp : energyParamDerivExpressions)
        compiledDerivExpressions.push_back(exp.createCompiledExpression());
    for (auto& exp : computedValueExpressions)
        compiledValueExpressions.push_back(exp.createCompiledExpression());
    for (int i = 0; i < threads.getNumThreads(); i++)
        threadData.push_back(new ThreadData(compiledEnergyExpression, energyVecExpression, compiledForceExpression, forceVecExpression, parameterNames,
                compiledDerivExpressions, computedValueNames, compiledValueExpressions, atomComputedValues));
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
    atomComputedValues.resize(computedValueNames.size(), vector<double>(numberOfAtoms));
    atomicCounter = 0;
    
    // Signal the threads to start running and wait for them to finish.
    
    threads.execute([&] (ThreadPool& threads, int threadIndex) { threadComputeForce(threads, threadIndex); });
    threads.waitForThreads(); // Computed values
    threads.resumeThreads();
    threads.waitForThreads(); // Interactions
    
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
    int numThreads = threads.getNumThreads();
    ThreadData& data = *threadData[threadIndex];
    for (auto& param : *globalParameters) {
        data.expressionSet.setVariable(data.expressionSet.getVariableIndex(param.first), param.second);
        try {
            float* p = data.energyVecExpression.getVariablePointer(param.first);
            for (int i = 0; i < blockSize; i++)
                p[i] = param.second;
        }
        catch (...) {
            // The expression doesn't use this parameter.
        }
        try {
            float* p = data.forceVecExpression.getVariablePointer(param.first);
            for (int i = 0; i < blockSize; i++)
                p[i] = param.second;
        }
        catch (...) {
            // The expression doesn't use this parameter.
        }
    }

    // Process computed values for this thread's subset of interactions.

    int numComputedValues = atomComputedValues.size();
    if (numComputedValues > 0) {
        int start = threadIndex*numberOfAtoms/numThreads;
        int end = (threadIndex+1)*numberOfAtoms/numThreads;
        for (int i = start; i < end; i++) {
            for (int j = 0; j < paramNames.size(); j++)
                data.particleParam[j] = atomParameters[i][j];
            for (int j = 0; j < numComputedValues; j++)
                atomComputedValues[j][i] = data.computedValueExpressions[j].evaluate();
        }
    }
    threads.syncThreads();

    // Compute this thread's subset of interactions.

    threadEnergy[threadIndex] = 0;
    double& energy = threadEnergy[threadIndex];
    float* forces = &(*threadForce)[threadIndex][0];
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
            for (int j = 0; j < paramNames.size(); j++) {
                data.particleParam[j*2] = atomParameters[atom1][j];
                data.particleParam[j*2+1] = atomParameters[atom2][j];
            }
            for (int j = 0; j < computedValueNames.size(); j++) {
                data.computedValues[j*2] = atomComputedValues[j][atom1];
                data.computedValues[j*2+1] = atomComputedValues[j][atom2];
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
            const int32_t* blockAtom = &neighborList->getSortedAtoms()[blockSize*blockIndex];
            const vector<int>& neighbors = neighborList->getBlockNeighbors(blockIndex);
            const auto& exclusions = neighborList->getBlockExclusions(blockIndex);
            if (data.energyParamDerivs.size() == 0)
                calculateBlockIxn(data, blockIndex, forces, energy, boxSize, invBoxSize);
            else {
                for (int i = 0; i < (int) neighbors.size(); i++) {
                    int first = neighbors[i];
                    for (int j = 0; j < paramNames.size(); j++)
                        data.particleParam[j*2] = atomParameters[first][j];
                    for (int j = 0; j < computedValueNames.size(); j++)
                        data.computedValues[j*2] = atomComputedValues[j][first];
                    for (int k = 0; k < blockSize; k++) {
                        if ((exclusions[i] & (1<<k)) == 0) {
                            int second = blockAtom[k];
                            for (int j = 0; j < paramNames.size(); j++)
                                data.particleParam[j*2+1] = atomParameters[second][j];
                            for (int j = 0; j < computedValueNames.size(); j++)
                                data.computedValues[j*2+1] = atomComputedValues[j][second];
                            calculateOneIxn(first, second, data, forces, energy, boxSize, invBoxSize);
                        }
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
                    for (int j = 0; j < paramNames.size(); j++) {
                        data.particleParam[j*2] = atomParameters[ii][j];
                        data.particleParam[j*2+1] = atomParameters[jj][j];
                    }
                    for (int j = 0; j < computedValueNames.size(); j++) {
                        data.computedValues[j*2] = atomComputedValues[j][ii];
                        data.computedValues[j*2+1] = atomComputedValues[j][jj];
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

void CpuCustomNonbondedForce::calculateBlockIxn(ThreadData& data, int blockIndex, float* forces, double& totalEnergy, const fvec4& boxSize, const fvec4& invBoxSize) {
    // Determine whether we need to apply periodic boundary conditions.

    PeriodicType periodicType;
    fvec4 blockCenter;
    if (!periodic) {
        periodicType = NoPeriodic;
        blockCenter = 0.0f;
    }
    else {
        const int32_t* blockAtom = &neighborList->getSortedAtoms()[blockSize*blockIndex];
        float minx, maxx, miny, maxy, minz, maxz;
        minx = maxx = posq[4*blockAtom[0]];
        miny = maxy = posq[4*blockAtom[0]+1];
        minz = maxz = posq[4*blockAtom[0]+2];
        for (int i = 1; i < blockSize; i++) {
            minx = min(minx, posq[4*blockAtom[i]]);
            maxx = max(maxx, posq[4*blockAtom[i]]);
            miny = min(miny, posq[4*blockAtom[i]+1]);
            maxy = max(maxy, posq[4*blockAtom[i]+1]);
            minz = min(minz, posq[4*blockAtom[i]+2]);
            maxz = max(maxz, posq[4*blockAtom[i]+2]);
        }
        blockCenter = fvec4(0.5f*(minx+maxx), 0.5f*(miny+maxy), 0.5f*(minz+maxz), 0.0f);
        if (!(minx < cutoffDistance || miny < cutoffDistance || minz < cutoffDistance ||
                maxx > boxSize[0]-cutoffDistance || maxy > boxSize[1]-cutoffDistance || maxz > boxSize[2]-cutoffDistance))
            periodicType = NoPeriodic;
        else if (triclinic)
            periodicType = PeriodicTriclinic;
        else if (0.5f*(boxSize[0]-(maxx-minx)) >= cutoffDistance &&
                 0.5f*(boxSize[1]-(maxy-miny)) >= cutoffDistance &&
                 0.5f*(boxSize[2]-(maxz-minz)) >= cutoffDistance)
            periodicType = PeriodicPerAtom;
        else
            periodicType = PeriodicPerInteraction;
    }

    // Call the appropriate version depending on what calculation is required for periodic boundary conditions.

    if (periodicType == NoPeriodic)
        calculateBlockIxnImpl<fvec4, NoPeriodic>(data, blockIndex, forces, totalEnergy, boxSize, invBoxSize, blockCenter);
    else if (periodicType == PeriodicPerAtom)
        calculateBlockIxnImpl<fvec4, PeriodicPerAtom>(data, blockIndex, forces, totalEnergy, boxSize, invBoxSize, blockCenter);
    else if (periodicType == PeriodicPerInteraction)
        calculateBlockIxnImpl<fvec4, PeriodicPerInteraction>(data, blockIndex, forces, totalEnergy, boxSize, invBoxSize, blockCenter);
    else if (periodicType == PeriodicTriclinic)
        calculateBlockIxnImpl<fvec4, PeriodicTriclinic>(data, blockIndex, forces, totalEnergy, boxSize, invBoxSize, blockCenter);
}

template <typename FVEC, int PERIODIC_TYPE>
void CpuCustomNonbondedForce::calculateBlockIxnImpl(ThreadData& data, int blockIndex, float* forces, double& totalEnergy, const fvec4& boxSize, const fvec4& invBoxSize, const fvec4& blockCenter) {
    // Load the positions and parameters of the atoms in the block.

    const int32_t* blockAtom = &neighborList->getSortedAtoms()[blockSize*blockIndex];
    fvec4 blockAtomPosq[blockSize];
    FVEC blockAtomForceX(0.0f), blockAtomForceY(0.0f), blockAtomForceZ(0.0f);
    FVEC blockAtomX, blockAtomY, blockAtomZ, blockAtomCharge;
    int numParams = paramNames.size();
    int numComputed = computedValueNames.size();
    for (int i = 0; i < blockSize; i++) {
        blockAtomPosq[i] = fvec4(posq+4*blockAtom[i]);
        if (PERIODIC_TYPE == PeriodicPerAtom)
            blockAtomPosq[i] -= floor((blockAtomPosq[i]-blockCenter)*invBoxSize+0.5f)*boxSize;
        for (int j = 0; j < numParams; j++)
            data.vecParticle1Params[j*blockSize+i] = atomParameters[blockAtom[i]][j];
        for (int j = 0; j < numComputed; j++)
            data.vecParticle1Values[j*blockSize+i] = atomComputedValues[j][blockAtom[i]];
    }
    transpose(blockAtomPosq, blockAtomX, blockAtomY, blockAtomZ, blockAtomCharge);
    const bool needPeriodic = (PERIODIC_TYPE == PeriodicPerInteraction || PERIODIC_TYPE == PeriodicTriclinic);
    const float invSwitchingInterval = 1/(cutoffDistance-switchingDistance);
    const FVEC cutoffDistanceSquared = cutoffDistance * cutoffDistance;

    // Loop over neighbors for this block.

    const auto& neighbors = neighborList->getBlockNeighbors(blockIndex);
    const auto& exclusions = neighborList->getBlockExclusions(blockIndex);
    FVEC partialEnergy = {};
    for (int i = 0; i < neighbors.size(); i++) {
        // Load the next neighbor.

        int atom = neighbors[i];
        for (int j = 0; j < numParams; j++)
            for (int k = 0; k < blockSize; k++)
                data.vecParticle2Params[j*blockSize+k] = atomParameters[atom][j];
        for (int j = 0; j < numComputed; j++)
            for (int k = 0; k < blockSize; k++)
                data.vecParticle2Values[j*blockSize+k] = atomComputedValues[j][atom];

        // Compute the distances to the block atoms.

        FVEC dx, dy, dz, r2;
        fvec4 atomPos(posq+4*atom);
        if (PERIODIC_TYPE == PeriodicPerAtom)
            atomPos -= floor((atomPos-blockCenter)*invBoxSize+0.5f)*boxSize;
        getDeltaR<FVEC, PERIODIC_TYPE>(atomPos, blockAtomX, blockAtomY, blockAtomZ, dx, dy, dz, r2, boxSize, invBoxSize);
        const auto exclNotMask = FVEC::expandBitsToMask(~exclusions[i]);
        const auto include = blendZero(r2 < cutoffDistanceSquared, exclNotMask);
        if (!any(include))
            continue; // No interactions to compute.

        // Compute the interactions.

        const auto inverseR = rsqrt(r2);
        const auto r = r2*inverseR;
        r.store(data.rvec.data());
        FVEC dEdR(data.forceVecExpression.evaluate());
        FVEC energy;
        if (includeEnergy)
            energy = FVEC(data.energyVecExpression.evaluate());
        if (useSwitch) {
            const auto t = blendZero((r-switchingDistance)*invSwitchingInterval, r>switchingDistance);
            const auto switchValue = 1+t*t*t*(-10.0f+t*(15.0f-t*6.0f));
            const auto switchDeriv = t*t*(-30.0f+t*(60.0f-t*30.0f))*invSwitchingInterval;
            dEdR = switchValue*dEdR + energy*switchDeriv;
            energy *= switchValue;
        }
        dEdR *= inverseR;

        // Accumulate forces and energies.

        if (includeEnergy) {
            energy = blendZero(energy, include);
            partialEnergy += energy;
        }
        dEdR = blendZero(dEdR, include);
        const auto fx = dx*dEdR;
        const auto fy = dy*dEdR;
        const auto fz = dz*dEdR;
        blockAtomForceX -= fx;
        blockAtomForceY -= fy;
        blockAtomForceZ -= fz;
        float* const atomForce = forces+4*atom;
        const fvec4 newAtomForce = fvec4(atomForce) + reduceToVec3(fx, fy, fz);
        newAtomForce.store(atomForce);
    }
    if (includeEnergy)
        totalEnergy += reduceAdd(partialEnergy);

    // Record the forces on the block atoms.

    fvec4 f[blockSize];
    transpose(blockAtomForceX, blockAtomForceY, blockAtomForceZ, 0.0f, f);
    for (int j = 0; j < blockSize; j++)
        (fvec4(forces+4*blockAtom[j])+f[j]).store(forces+4*blockAtom[j]);
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

template <typename FVEC, int PERIODIC_TYPE>
void CpuCustomNonbondedForce::getDeltaR(const fvec4& posI, const FVEC& x, const FVEC& y, const FVEC& z, FVEC& dx, FVEC& dy, FVEC& dz, FVEC& r2, const fvec4& boxSize, const fvec4& invBoxSize) const {
    dx = x-posI[0];
    dy = y-posI[1];
    dz = z-posI[2];
    if (PERIODIC_TYPE == PeriodicTriclinic) {
        const auto scale3 = floor(dz*recipBoxSize[2]+0.5f);
        dx -= scale3*periodicBoxVectors[2][0];
        dy -= scale3*periodicBoxVectors[2][1];
        dz -= scale3*periodicBoxVectors[2][2];
        const auto scale2 = floor(dy*recipBoxSize[1]+0.5f);
        dx -= scale2*periodicBoxVectors[1][0];
        dy -= scale2*periodicBoxVectors[1][1];
        const auto scale1 = floor(dx*recipBoxSize[0]+0.5f);
        dx -= scale1*periodicBoxVectors[0][0];
    }
    else if (PERIODIC_TYPE == PeriodicPerInteraction) {
        dx -= round(dx*invBoxSize[0])*boxSize[0];
        dy -= round(dy*invBoxSize[1])*boxSize[1];
        dz -= round(dz*invBoxSize[2])*boxSize[2];
    }
    r2 = dx*dx + dy*dy + dz*dz;
}
