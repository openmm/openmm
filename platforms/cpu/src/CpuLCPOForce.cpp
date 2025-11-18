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

CpuLCPOForce::Neighbors::Neighbors(int numParticles, int maxNumNeighborsGuess) :
        numParticles(numParticles), maxNumNeighborsFound(0), numNeighbors(numParticles) {
    // Round up to the nearest multiple of 4 so that neighbor searches can be vectorized.
    maxNumNeighbors = (maxNumNeighborsGuess + 3) & -4;
    int numItems = numParticles * maxNumNeighbors;
    indices.assign(numItems, -1);
    data.resize(numItems);
}

void CpuLCPOForce::Neighbors::clear() {
    maxNumNeighborsFound = 0;
    numNeighbors.assign(numParticles, 0);
    indices.assign(numParticles * maxNumNeighbors, -1);
}

void CpuLCPOForce::Neighbors::insert(const NeighborInfo& info) {
    int& iNumNeighbors = numNeighbors[info.i];
    int k = iNumNeighbors++;
    if (k < maxNumNeighbors) {
        k += info.i * maxNumNeighbors;
        indices[k] = info.j;
        this->data[k] = info.data;
    }
    maxNumNeighborsFound = max(maxNumNeighborsFound, iNumNeighbors);
}

void CpuLCPOForce::Neighbors::getNeighbors(int i, int& iNumNeighbors, const int*& iIndices, const NeighborData*& iData) const {
    iNumNeighbors = numNeighbors[i];
    iIndices = &indices[i * maxNumNeighbors];
    iData = &data[i * maxNumNeighbors];
}

bool CpuLCPOForce::Neighbors::isNeighbor(int i, int j, NeighborData& data) {
    int start = i * maxNumNeighbors;
    int stop = start + ((numNeighbors[i] + 3) & -4);
    ivec4 j4 = ivec4(j);
    for (int k = start; k < stop; k += 4) {
        ivec4 mask = (ivec4(&indices[k]) == j4);
        if (any(mask)) {
            // Exactly one mask element should be set to -1 and the others to 0.
            // Use this to extract the appropriate index and retrieve the data.
            mask *= ivec4(0, 1, 2, 3);
            data = this->data[k - (mask[1] + mask[2] + mask[3])];
            return true;
        }
    }
    return false;
}

CpuLCPOForce::CpuLCPOForce(ThreadPool& threads, const vector<int>& activeParticles, const vector<int>& activeParticlesInv, const vector<fvec4>& parameters, bool usePeriodic) :
        numParticles(activeParticlesInv.size()), numActiveParticles(activeParticles.size()),
        threads(threads), activeParticles(activeParticles), activeParticlesInv(activeParticlesInv), parameters(parameters), usePeriodic(usePeriodic),
        neighborList(4), exclusions(numActiveParticles), neighbors(numActiveParticles) {
    for (int i = 0; i < numActiveParticles; i++) {
        exclusions[i].insert(i);
    }
    updateRadii();
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

    isTriclinic = false;
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

    neighborList.computeNeighborList(numActiveParticles, posq, exclusions, boxVectors, usePeriodic, cutoff, threads, &activeParticles);

    // Process neighbors from the neighbor list.

    atomicCounter = 0;
    threads.execute([&] (ThreadPool& threads, int threadIndex) { threadExecute(threads, threadIndex); });
    threads.waitForThreads();

    // Compile all of the neighbor information from the threads into a single collection.

    neighbors.clear();
    while (true) {
        for (auto threadNeighborInfo : threadNeighbors) {
            for (NeighborInfo neighborInfo : threadNeighborInfo) {
                neighbors.insert(neighborInfo);
            }
        }
        int maxNumNeighborsNeeded;
        if (neighbors.didOverflow(maxNumNeighborsNeeded)) {
            neighbors = Neighbors(numActiveParticles, (int) (maxNumNeighborsNeeded * 1.1));
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

void CpuLCPOForce::threadExecute(ThreadPool& threads, int threadIndex) {
    int numThreads = threads.getNumThreads();

    // Process neighbors from the neighbor list.

    vector<NeighborInfo>& threadNeighborInfo = threadNeighbors[threadIndex];
    threadNeighborInfo.clear();
    while (true) {
        int blockIndex = atomicCounter++;
        if (blockIndex >= neighborList.getNumBlocks()) {
            break;
        }
        if (usePeriodic) {
            if (isTriclinic) {
                processNeighborListBlock<true, true>(blockIndex, threadNeighborInfo);
            }
            else {
                processNeighborListBlock<true, false>(blockIndex, threadNeighborInfo);
            }
        }
        else {
            processNeighborListBlock<false, false>(blockIndex, threadNeighborInfo);
        }
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
            float p2i = parameters[i][P2Index];
            float p3i = parameters[i][P3Index];
            float p4i = parameters[i][P4Index];
            float energy = 0.0f;
            fvec4 iForce(0.0f);

            int iNumNeighbors;
            const int* iIndices;
            const NeighborData* iData;
            neighbors.getNeighbors(i, iNumNeighbors, iIndices, iData);
            for (int jNeighbor = 0; jNeighbor < iNumNeighbors; jNeighbor++) {
                int j = iIndices[jNeighbor];
                fvec4 Aij = iData[jNeighbor].ijData;
                fvec4 Aji = iData[jNeighbor].jiData;
                float p2j = parameters[j][P2Index];
                float p3j = parameters[j][P3Index];
                float p4j = parameters[j][P4Index];
                fvec4 jForce(0.0f);

                // Two-body term.

                fvec4 ijForce2 = p2i * Aij + p2j * Aji;
                iForce += ijForce2;
                jForce -= ijForce2;
                energy += ijForce2[3];

                // Three-body term: includes all pairs (j, k) of neighbors of i that are also neighbors of each other.

                int jNumNeighbors;
                const int* jIndices;
                const NeighborData* jData;
                neighbors.getNeighbors(j, jNumNeighbors, jIndices, jData);
                for (int kNeighbor = 0; kNeighbor < jNumNeighbors; kNeighbor++) {
                    int k = jIndices[kNeighbor];
                    NeighborData ikNeighborData;
                    if (!neighbors.isNeighbor(i, k, ikNeighborData)) {
                        continue;
                    }
                    fvec4 Aik = ikNeighborData.ijData;
                    fvec4 Aki = ikNeighborData.jiData;
                    fvec4 Ajk = jData[kNeighbor].ijData;
                    fvec4 Akj = jData[kNeighbor].jiData;
                    float p3k = parameters[k][P3Index];
                    float p4k = parameters[k][P4Index];

                    fvec4 ijkForce34 = (p3i + p4i * Aij[3]) * Ajk + (p3i + p4i * Aik[3]) * Akj;
                    fvec4 jikForce34 = (p3j + p4j * Aji[3]) * Aik + (p3j + p4j * Ajk[3]) * Aki;
                    fvec4 kijForce34 = (p3k + p4k * Aki[3]) * Aij + (p3k + p4k * Akj[3]) * Aji;

                    fvec4 ijkForce4 = p4i * Aij * Ajk[3];
                    fvec4 ikjForce4 = p4i * Aik * Akj[3];
                    fvec4 jikForce4 = p4j * Aji * Aik[3];
                    fvec4 jkiForce4 = p4j * Ajk * Aki[3];
                    fvec4 kijForce4 = p4k * Aki * Aij[3];
                    fvec4 kjiForce4 = p4k * Akj * Aji[3];

                    energy += ijkForce34[3] + jikForce34[3] + kijForce34[3];
                    iForce += jikForce34 + kijForce34 + ijkForce4 + ikjForce4 + jikForce4 + kijForce4;
                    jForce += ijkForce34 - kijForce34 + jkiForce4 - jikForce4 - ijkForce4 + kjiForce4;
                    fvec4 kForce = ijkForce34 + jikForce34 + kijForce4 + kjiForce4 + ikjForce4 + jkiForce4;

                    float * kForcePointer = &forces[4 * activeParticles[k]];
                    (fvec4(kForcePointer) - kForce).store(kForcePointer);
                }
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
    const int* iBlock = &neighborList.getSortedAtoms()[4 * blockIndex];
    fvec4 iBlockPosq[4] {
        fvec4(posq + 4 * activeParticles[iBlock[0]]),
        fvec4(posq + 4 * activeParticles[iBlock[1]]),
        fvec4(posq + 4 * activeParticles[iBlock[2]]),
        fvec4(posq + 4 * activeParticles[iBlock[3]])
    };
    fvec4 iBlockX, iBlockY, iBlockZ, iBlockQ;
    transpose(iBlockPosq, iBlockX, iBlockY, iBlockZ, iBlockQ);

    fvec4 iBlockRadius(parameters[iBlock[0]][RadiusIndex], parameters[iBlock[1]][RadiusIndex], parameters[iBlock[2]][RadiusIndex], parameters[iBlock[3]][RadiusIndex]);

    CpuNeighborList::NeighborIterator iterator = neighborList.getNeighborIterator(blockIndex);
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
        fvec4 iBlockForce = iBlockRadiusPi * (-1.0f + deltaRadiusRSq);
        fvec4 jForce = jRadiusPi * (-1.0f - deltaRadiusRSq);

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
