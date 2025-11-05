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

using namespace std;
using namespace OpenMM;

CpuLCPOForce::CpuLCPOForce(ThreadPool& threads, const std::vector<int>& indices, const std::vector<int>& particles, const std::vector<fvec4>& parameters, bool usePeriodic) :
        numParticles(indices.size()), numActiveParticles(particles.size()),
        threads(threads), indices(indices), particles(particles), parameters(parameters), usePeriodic(usePeriodic), exclusions(numActiveParticles) {

    // Prepare a neighbor list with only self-interactions excluded.
    neighborList = new CpuNeighborList(4);
    for (int i = 0; i < numActiveParticles; i++) {
        exclusions[i].insert(i);
    }
    neighbors.resize(numActiveParticles);
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
    threadEnergy.resize(numThreads);
    this->threadForce = &threadForce;

    neighborList->computeNeighborList(numActiveParticles, posq, exclusions, boxVectors, usePeriodic, cutoff, threads, &particles);

    if (usePeriodic) {
        recipBoxSize[0] = (float)(1.0 / boxVectors[0][0]);
        recipBoxSize[1] = (float)(1.0 / boxVectors[1][1]);
        recipBoxSize[2] = (float)(1.0 / boxVectors[2][2]);
    }

    // TODO: vectorize / parallelize

    for (int i = 0; i < numActiveParticles; i++) {
        neighbors[i].clear();
    }

    for (int blockIndex = 0; blockIndex < neighborList->getNumBlocks(); blockIndex++) {
        const vector<int>& blockNeighbors = neighborList->getBlockNeighbors(blockIndex);
        const auto& exclusions = neighborList->getBlockExclusions(blockIndex);
        int numNeighbors = blockNeighbors.size();
        for (int iBlock = 0; iBlock < 4; iBlock++) {
            int i = neighborList->getSortedAtoms()[4 * blockIndex + iBlock];
            for (int jNeighbor = 0; jNeighbor < numNeighbors; jNeighbor++) {
                if ((exclusions[jNeighbor] & (1 << iBlock)) == 0) {
                    int j = blockNeighbors[jNeighbor];

                    float iRadius = parameters[i][RadiusIndex];
                    float jRadius = parameters[j][RadiusIndex];

                    // Only include particles close enough to each other.

                    int iPos = 4 * particles[i];
                    int jPos = 4 * particles[j];
                    float deltaX = posq[jPos] - posq[iPos];
                    float deltaY = posq[jPos + 1] - posq[iPos + 1];
                    float deltaZ = posq[jPos + 2] - posq[iPos + 2];

                    if (usePeriodic) {
                        float scale3 = floor(deltaZ * recipBoxSize[2] + 0.5f);
                        deltaX -= scale3 * boxVectors[2][0];
                        deltaY -= scale3 * boxVectors[2][1];
                        deltaZ -= scale3 * boxVectors[2][2];
                        float scale2 = floor(deltaY * recipBoxSize[1] + 0.5f);
                        deltaX -= scale2 * boxVectors[1][0];
                        deltaY -= scale2 * boxVectors[1][1];
                        float scale1 = floor(deltaX * recipBoxSize[0] + 0.5f);
                        deltaX -= scale1 * boxVectors[0][0];
                    }

                    float r = sqrt(deltaX * deltaX + deltaY * deltaY + deltaZ * deltaZ);
                    if (r >= iRadius + jRadius) {
                        continue;
                    }
                    float rRecip = 1.0f / r;
                    deltaX *= rRecip;
                    deltaY *= rRecip;
                    deltaZ *= rRecip;

                    // Precompute and store the buried areas and their derivatives.

                    float iRadiusPi = PI_M * iRadius;
                    float jRadiusPi = PI_M * jRadius;
                    float deltaRadiusR = (iRadius * iRadius - jRadius * jRadius) * rRecip;
                    float deltaRadiusRSq = deltaRadiusR * rRecip;
                    float iForce = iRadiusPi * (deltaRadiusRSq - 1.0);
                    float jForce = jRadiusPi * (deltaRadiusRSq + 1.0);
                    neighbors[i][j] = fvec4(iRadiusPi * (2.0 * iRadius - r - deltaRadiusR), iForce * deltaX, iForce * deltaY, iForce * deltaZ);
                    neighbors[j][i] = fvec4(jRadiusPi * (2.0 * jRadius - r + deltaRadiusR), jForce * deltaX, jForce * deltaY, jForce * deltaZ);
                }
            }
        }
    }

    atomicCounter = 0;
    threads.execute([&] (ThreadPool& threads, int threadIndex) { threadExecute(threads, threadIndex); });
    threads.waitForThreads();

    if (includeEnergy) {
        for (int i = 0; i < numThreads; i++) {
            energy += threadEnergy[i];
        }
    }
}

void CpuLCPOForce::threadExecute(ThreadPool& threads, int threadIndex) {
    int numThreads = threads.getNumThreads();
    int groupSize = max(1, numActiveParticles / (10 * numThreads));
    float* forces = &(*threadForce)[threadIndex][0];
    threadEnergy[threadIndex] = 0.0;

    // TODO: vectorize

    while (true) {
        int start = atomicCounter.fetch_add(groupSize);
        if (start >= numActiveParticles) {
            break;
        }
        int end = min(start + groupSize, numActiveParticles);

        for (int i = start; i < end; i++) {
            float p2 = parameters[i][P2Index];
            float p3 = parameters[i][P3Index];
            float p4 = parameters[i][P4Index];
            float term2 = 0.0f;
            float term3 = 0.0f;
            float term4 = 0.0f;
            float iForce[3] = {0.0f, 0.0f, 0.0f};

            for (auto jNeighbor : neighbors[i]) {
                int j = jNeighbor.first;
                float Aij = jNeighbor.second[0];
                float dAijx = jNeighbor.second[1];
                float dAijy = jNeighbor.second[2];
                float dAijz = jNeighbor.second[3];
                float jForce[3] = {0.0f, 0.0f, 0.0f};

                // Two-body term.

                term2 += Aij;
                float ijForce2Bodyx = p2 * dAijx;
                float ijForce2Bodyy = p2 * dAijy;
                float ijForce2Bodyz = p2 * dAijz;
                iForce[0] += ijForce2Bodyx;
                iForce[1] += ijForce2Bodyy;
                iForce[2] += ijForce2Bodyz;
                jForce[0] -= ijForce2Bodyx;
                jForce[1] -= ijForce2Bodyy;
                jForce[2] -= ijForce2Bodyz;

                // Three-body term: includes all pairs (j, k) of neighbors of i that
                // are also neighbors of each other.

                for (auto kNeighbor : neighbors[j]) {
                    int k = kNeighbor.first;
                    float Ajk = kNeighbor.second[0];
                    float dAjkx = kNeighbor.second[1];
                    float dAjky = kNeighbor.second[2];
                    float dAjkz = kNeighbor.second[3];

                    if (neighbors[i].find(k) == neighbors[i].end()) {
                        continue;
                    }

                    term3 += Ajk;
                    term4 += Aij * Ajk;
                    float jkForce3Bodyx = (p3 + p4 * Aij) * dAjkx;
                    float jkForce3Bodyy = (p3 + p4 * Aij) * dAjky;
                    float jkForce3Bodyz = (p3 + p4 * Aij) * dAjkz;
                    float ijForce3Bodyx = p4 * dAijx * Ajk;
                    float ijForce3Bodyy = p4 * dAijy * Ajk;
                    float ijForce3Bodyz = p4 * dAijz * Ajk;
                    iForce[0] += ijForce3Bodyx;
                    iForce[1] += ijForce3Bodyy;
                    iForce[2] += ijForce3Bodyz;
                    jForce[0] += jkForce3Bodyx - ijForce3Bodyx;
                    jForce[1] += jkForce3Bodyy - ijForce3Bodyy;
                    jForce[2] += jkForce3Bodyz - ijForce3Bodyz;
                    forces[4 * particles[k]] -= jkForce3Bodyx;
                    forces[4 * particles[k] + 1] -= jkForce3Bodyy;
                    forces[4 * particles[k] + 2] -= jkForce3Bodyz;
                }

                forces[4 * particles[j]] += jForce[0];
                forces[4 * particles[j] + 1] += jForce[1];
                forces[4 * particles[j] + 2] += jForce[2];
            }

            threadEnergy[threadIndex] += p2 * term2 + p3 * term3 + p4 * term4;
            forces[4 * particles[i]] += iForce[0];
            forces[4 * particles[i] + 1] += iForce[1];
            forces[4 * particles[i] + 2] += iForce[2];
        }
    }
}
