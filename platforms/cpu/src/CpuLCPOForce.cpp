/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit.                   *
 * See https://openmm.org/development.                                        *
 *                                                                            *
 * Portions copyright (c) 2025 Stanford University and the Authors.           *
 * Authors: Evan Pretti                                                       *
 * Contributors:                                                              *
 *                                                                            *
 * Permission is hereby granted, free of charge, to any person obtaining a    *
 * copy of this software and associated documentation files (the "Software"), *
 * to deal in the Software without restriction, including without limitation  *
 * the rights to use, copy, modify, merge, publish, distribute, sublicense,   *
 * and/or sell copies of the Software, and to permit persons to whom the      *
 * Software is furnished to do so, subject to the following conditions:       *
 *                                                                            *
 * The above copyright notice and this permission notice shall be included in *
 * all copies or substantial portions of the Software.                        *
 *                                                                            *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR *
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,   *
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL    *
 * THE AUTHORS, CONTRIBUTORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,    *
 * DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR      *
 * OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE  *
 * USE OR OTHER DEALINGS IN THE SOFTWARE.                                     *
 * -------------------------------------------------------------------------- */

#include "CpuLCPOForce.h"
#include "SimTKOpenMMRealType.h"
#include "openmm/OpenMMException.h"

using namespace std;
using namespace OpenMM;

CpuLCPOForce::Neighbors::Neighbors(int numParticles, int maxNumNeighbors) :
        numParticles(numParticles), maxNumNeighbors(maxNumNeighbors), maxNumNeighborsFound(0),
        numNeighbors(numParticles), indices(numParticles * maxNumNeighbors), data(numParticles * maxNumNeighbors) {
}

void CpuLCPOForce::Neighbors::clear() {
    maxNumNeighborsFound = 0;
    numNeighbors.assign(numParticles, 0);
}

void CpuLCPOForce::Neighbors::insert(int i, int j, fvec4 ijData, fvec4 jiData) {
    insert(i, j, ijData);
    insert(j, i, jiData);
}

void CpuLCPOForce::Neighbors::getNeighbors(int i, int& iNumNeighbors, const int*& iIndices, const fvec4*& iData) const {
    iNumNeighbors = numNeighbors[i];
    iIndices = &indices[i * maxNumNeighbors];
    iData = &data[i * maxNumNeighbors];
}

void CpuLCPOForce::Neighbors::insert(int i, int j, fvec4 ijData) {
    int k = numNeighbors[i]++;
    if (k < maxNumNeighbors) {
        k += i * maxNumNeighbors;
        indices[k] = j;
        data[k] = ijData;
    }
    maxNumNeighborsFound = max(maxNumNeighborsFound, numNeighbors[i]);
}

CpuLCPOForce::CpuLCPOForce(ThreadPool& threads, const vector<int>& activeParticles, const vector<int>& activeParticlesInv, const vector<fvec4>& parameters, bool usePeriodic) :
        numParticles(activeParticlesInv.size()), numActiveParticles(activeParticles.size()), neighbors(numActiveParticles),
        threads(threads), activeParticles(activeParticles), activeParticlesInv(activeParticlesInv), parameters(parameters), usePeriodic(usePeriodic), exclusions(numActiveParticles) {
    // Prepare a neighbor list with only self-interactions excluded.
    neighborList = new CpuNeighborList(4);
    for (int i = 0; i < numActiveParticles; i++) {
        exclusions[i].insert(i);
    }
    updateRadii();
}

CpuLCPOForce::~CpuLCPOForce() {
    delete neighborList;
}

void CpuLCPOForce::updateRadii() {
    cutoff = 0.0f;
    for (int i = 0; i < numActiveParticles; i++) {
        cutoff = max(cutoff, parameters[i][RadiusIndex]);
    }
    cutoff *= 2.0f;
}

void CpuLCPOForce::execute(Vec3* boxVectors, AlignedArray<float>& posq, vector<AlignedArray<float> >& threadForce, bool includeForces, bool includeEnergy, double& energy) {
    if (!numActiveParticles) {
        return;
    }

    int numThreads = threads.getNumThreads();
    this->posq = &posq[0];
    threadEnergy.resize(numThreads);
    this->threadForce = &threadForce;
    threadNeighbors.resize(numThreads);

    bool isTriclinic = false;
    if (usePeriodic) {
        double minAllowedSize = 1.999999 * cutoff;
        if (boxVectors[0][0] < minAllowedSize || boxVectors[1][1] < minAllowedSize || boxVectors[2][2] < minAllowedSize) {
            throw OpenMMException("The periodic box size is less than twice the required cutoff for LCPO.");
        }

        this->boxVectors = boxVectors;
        boxSize = fvec4(boxVectors[0][0], boxVectors[1][1], boxVectors[2][2], 0.0f);
        recipBoxSize = fvec4(1.0 / boxVectors[0][0], 1.0 / boxVectors[1][1], 1.0 / boxVectors[2][2], 0.0f);
        boxVec4[0] = fvec4(boxVectors[0][0], boxVectors[0][1], boxVectors[0][2], 0.0f);
        boxVec4[1] = fvec4(boxVectors[1][0], boxVectors[1][1], boxVectors[1][2], 0.0f);
        boxVec4[2] = fvec4(boxVectors[2][0], boxVectors[2][1], boxVectors[2][2], 0.0f);
        isTriclinic = (boxVectors[0][1] != 0.0 || boxVectors[0][2] != 0.0 ||
                boxVectors[1][0] != 0.0 || boxVectors[1][2] != 0.0 ||
                boxVectors[2][0] != 0.0 || boxVectors[2][1] != 0.0);
    }

    neighborList->computeNeighborList(numActiveParticles, posq, exclusions, boxVectors, usePeriodic, cutoff, threads, &activeParticles);

    // Process neighbors from the neighbor list.

    atomicCounter = 0;
    if (usePeriodic) {
        if (isTriclinic) {
            threads.execute([&] (ThreadPool& threads, int threadIndex) { threadExecute<true, true>(threads, threadIndex); });
        }
        else {
            threads.execute([&] (ThreadPool& threads, int threadIndex) { threadExecute<true, false>(threads, threadIndex); });
        }
    }
    else {
        threads.execute([&] (ThreadPool& threads, int threadIndex) { threadExecute<false, false>(threads, threadIndex); });
    }
    threads.waitForThreads();

    // Compile all of the neighbor information from the threads into a single collection.

    neighbors.clear();
    while (true) {
        for (auto threadNeighborInfo : threadNeighbors) {
            for (NeighborInfo neighborInfo : threadNeighborInfo) {
                neighbors.insert(neighborInfo.i, neighborInfo.j, neighborInfo.ijData, neighborInfo.jiData);
            }
        }
        int maxNumNeighborsNeeded;
        if (neighbors.didOverflow(maxNumNeighborsNeeded)) {
            neighbors = Neighbors(numActiveParticles, maxNumNeighborsNeeded);
        }
        else {
            break;
        }
    }

    // Accumulate energies and forces.

    atomicCounter = 0;
    threads.resumeThreads();
    threads.waitForThreads();

    if (includeEnergy) {
        for (int i = 0; i < numThreads; i++) {
            energy += threadEnergy[i];
        }
    }
}

template<bool USE_PERIODIC, bool IS_TRICLINIC>
void CpuLCPOForce::threadExecute(ThreadPool& threads, int threadIndex) {
    int numThreads = threads.getNumThreads();

    // Process neighbors from the neighbor list.

    vector<NeighborInfo>& threadNeighborInfo = threadNeighbors[threadIndex];
    threadNeighborInfo.clear();
    while (true) {
        int blockIndex = atomicCounter++;
        if (blockIndex >= neighborList->getNumBlocks()) {
            break;
        }
        processNeighborListBlock<USE_PERIODIC, IS_TRICLINIC>(blockIndex, threadNeighborInfo);
    }

    threads.syncThreads();

    // Accumulate energies and forces.

    int groupSize = max(1, numActiveParticles / (10 * numThreads));
    threadEnergy[threadIndex] = 0.0;
    float* forces = &(*threadForce)[threadIndex][0];
    while (true) {
        int start = atomicCounter.fetch_add(groupSize);
        if (start >= numActiveParticles) {
            break;
        }
        int end = min(start + groupSize, numActiveParticles);

        for (int i = start; i < end; i++) {
            fvec4 iPos(posq + 4 * activeParticles[i]);
            float iRadius = parameters[i][RadiusIndex];
            float p2 = parameters[i][P2Index];
            float p3 = parameters[i][P3Index];
            float p4 = parameters[i][P4Index];
            float energy = 0.0f;
            fvec4 iForce(0.0f);

            int iNumNeighbors;
            const int* iIndices;
            const fvec4* iData;
            neighbors.getNeighbors(i, iNumNeighbors, iIndices, iData);
            for (int jNeighbor = 0; jNeighbor < iNumNeighbors; jNeighbor++) {
                int j = iIndices[jNeighbor];
                fvec4 Aij = iData[jNeighbor];
                fvec4 jForce(0.0f);

                // Two-body term.

                fvec4 ijForce2Body = p2 * Aij;
                iForce += ijForce2Body;
                jForce -= ijForce2Body;

                // Three-body term: includes all pairs (j, k) of neighbors of i
                // that are also neighbors of each other.

                int jNumNeighbors;
                const int* jIndices;
                const fvec4* jData;
                neighbors.getNeighbors(j, jNumNeighbors, jIndices, jData);
                for (int kNeighbor = 0; kNeighbor < jNumNeighbors; kNeighbor++) {

                    // We could check to see if k is in the list of neighbors of
                    // i but it is faster to simply recompute the distance.

                    int k = jIndices[kNeighbor];
                    if (k == i) {
                        continue;
                    }
                    fvec4 kPos(posq + 4 * activeParticles[k]);
                    fvec4 delta = kPos - iPos;
                    if (USE_PERIODIC) {
                        if (IS_TRICLINIC) {
                            delta -= boxVec4[2] * floorf(delta[2] * recipBoxSize[2] + 0.5f);
                            delta -= boxVec4[1] * floorf(delta[1] * recipBoxSize[1] + 0.5f);
                            delta -= boxVec4[0] * floorf(delta[0] * recipBoxSize[0] + 0.5f);
                        }
                        else {
                            delta -= boxSize * round(delta * recipBoxSize);
                        }
                    }
                    float rCut = iRadius + parameters[k][RadiusIndex];
                    if (dot3(delta, delta) >= rCut * rCut) {
                        continue;
                    }
                    fvec4 Ajk = jData[kNeighbor];

                    fvec4 jkForce3Body = (p3 + p4 * Aij[3]) * Ajk;
                    fvec4 ijForce3Body = p4 * Aij * Ajk[3];
                    iForce += ijForce3Body;
                    jForce += jkForce3Body - ijForce3Body;

                    energy += jkForce3Body[3];
                    float * kForcePointer = &forces[4 * activeParticles[k]];
                    (fvec4(kForcePointer) - jkForce3Body).store(kForcePointer);
                }

                energy += ijForce2Body[3];
                float * jForcePointer = &forces[4 * activeParticles[j]];
                (fvec4(jForcePointer) + jForce).store(jForcePointer);
            }

            threadEnergy[threadIndex] += energy;
            float * iForcePointer = &forces[4 * activeParticles[i]];
            (fvec4(iForcePointer) + iForce).store(iForcePointer);
        }
    }
}

template<bool USE_PERIODIC, bool IS_TRICLINIC>
void CpuLCPOForce::processNeighborListBlock(int blockIndex, vector<NeighborInfo>& threadNeighborInfo) {
    const int* iBlock = &neighborList->getSortedAtoms()[4 * blockIndex];
    fvec4 iBlockPosq[4] {
        fvec4(posq + 4 * activeParticles[iBlock[0]]),
        fvec4(posq + 4 * activeParticles[iBlock[1]]),
        fvec4(posq + 4 * activeParticles[iBlock[2]]),
        fvec4(posq + 4 * activeParticles[iBlock[3]])
    };
    fvec4 iBlockX, iBlockY, iBlockZ, iBlockQ;
    transpose(iBlockPosq, iBlockX, iBlockY, iBlockZ, iBlockQ);

    fvec4 iBlockRadius(parameters[iBlock[0]][RadiusIndex], parameters[iBlock[1]][RadiusIndex], parameters[iBlock[2]][RadiusIndex], parameters[iBlock[3]][RadiusIndex]);

    CpuNeighborList::NeighborIterator iterator = neighborList->getNeighborIterator(blockIndex);
    while (iterator.next()) {
        int j = iterator.getNeighbor();
        fvec4 jPos(posq + 4 * activeParticles[j]);

        // Only include particles close enough to each other.

        fvec4 deltaX = jPos[0] - iBlockX;
        fvec4 deltaY = jPos[1] - iBlockY;
        fvec4 deltaZ = jPos[2] - iBlockZ;
        if (USE_PERIODIC) {
            if (IS_TRICLINIC) {
                fvec4 scale3 = round(deltaZ * recipBoxSize[2]);
                deltaX -= scale3 * boxVec4[2][0];
                deltaY -= scale3 * boxVec4[2][1];
                deltaZ -= scale3 * boxVec4[2][2];
                fvec4 scale2 = round(deltaY * recipBoxSize[1]);
                deltaX -= scale2 * boxVec4[1][0];
                deltaY -= scale2 * boxVec4[1][1];
                fvec4 scale1 = round(deltaX * recipBoxSize[0]);
                deltaX -= scale1 * boxVec4[0][0];
            }
            else {
                deltaX -= boxSize[0] * round(deltaX * recipBoxSize[0]);
                deltaY -= boxSize[1] * round(deltaY * recipBoxSize[1]);
                deltaZ -= boxSize[2] * round(deltaZ * recipBoxSize[2]);
            }
        }
        fvec4 r2 = deltaX * deltaX + deltaY * deltaY + deltaZ * deltaZ;

        float jRadius = parameters[j][RadiusIndex];
        fvec4 rCut = iBlockRadius + jRadius;
        ivec4 include = blendZero(r2 < rCut * rCut, fvec4::expandBitsToMask(~iterator.getExclusions()));
        if (!any(include)) {
            continue;
        }
        fvec4 r = sqrt(r2);
        fvec4 rRecip = 1.0f / r;
        deltaX *= rRecip;
        deltaY *= rRecip;
        deltaZ *= rRecip;

        // Precompute and store the buried areas and their derivatives.

        fvec4 iBlockRadiusPi = PI_M * iBlockRadius;
        float jRadiusPi = PI_M * jRadius;
        fvec4 deltaRadiusR = (iBlockRadius * iBlockRadius - jRadius * jRadius) * rRecip;
        fvec4 deltaRadiusRSq = deltaRadiusR * rRecip;
        fvec4 iBlockForce = iBlockRadiusPi * (deltaRadiusRSq - 1.0f);
        fvec4 jForce = jRadiusPi * (deltaRadiusRSq + 1.0f);

        fvec4 ijData0 = iBlockForce * deltaX;
        fvec4 ijData1 = iBlockForce * deltaY;
        fvec4 ijData2 = iBlockForce * deltaZ;
        fvec4 ijData3 = iBlockRadiusPi * (2.0f * iBlockRadius - r - deltaRadiusR);
        fvec4 jiData0 = jForce * deltaX;
        fvec4 jiData1 = jForce * deltaY;
        fvec4 jiData2 = jForce * deltaZ;
        fvec4 jiData3 = jRadiusPi * (2.0f * jRadius - r + deltaRadiusR);

        transpose(ijData0, ijData1, ijData2, ijData3);
        transpose(jiData0, jiData1, jiData2, jiData3);

        if (include[0]) threadNeighborInfo.push_back({iBlock[0], j, ijData0, jiData0});
        if (include[1]) threadNeighborInfo.push_back({iBlock[1], j, ijData1, jiData1});
        if (include[2]) threadNeighborInfo.push_back({iBlock[2], j, ijData2, jiData2});
        if (include[3]) threadNeighborInfo.push_back({iBlock[3], j, ijData3, jiData3});
    }
}
