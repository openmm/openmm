typedef struct {
    real4 data12;
    real4 data21;
} NeighborData;

/**
 * Load the position of an atom from global memory.
 */
inline DEVICE real3 loadCondensedPos(GLOBAL const real* RESTRICT condensedPos, int atom) {
    int offset = 3 * atom;
    return make_real3(condensedPos[offset], condensedPos[offset + 1], condensedPos[offset + 2]);
}

/**
 * Record the force on an atom to global memory.
 */
inline DEVICE void storeForce(int atom, real3 force, GLOBAL mm_ulong* RESTRICT forceBuffers) {
    ATOMIC_ADD(&forceBuffers[atom], (mm_ulong) realToFixedPoint(force.x));
    ATOMIC_ADD(&forceBuffers[atom + PADDED_NUM_ATOMS], (mm_ulong) realToFixedPoint(force.y));
    ATOMIC_ADD(&forceBuffers[atom + 2 * PADDED_NUM_ATOMS], (mm_ulong) realToFixedPoint(force.z));
}

/**
 * Compute the difference between two vectors, taking periodic boundary
 * conditions into account and setting the fourth component to the squared
 * magnitude.
 */
inline DEVICE real4 delta(real3 vec1, real3 vec2, real4 periodicBoxSize, real4 invPeriodicBoxSize,
        real4 periodicBoxVecX, real4 periodicBoxVecY, real4 periodicBoxVecZ) {

    real4 result = make_real4(vec1.x - vec2.x, vec1.y - vec2.y, vec1.z - vec2.z, 0);
#ifdef USE_PERIODIC
    APPLY_PERIODIC_TO_DELTA(result)
#endif
    result.w = result.x * result.x + result.y * result.y + result.z * result.z;
    return result;
}

/**
 * Copies positions of active particles to a separate array for performance.
 */
KERNEL void condensePos(GLOBAL const int* RESTRICT activeParticles, GLOBAL const real4* RESTRICT posq, GLOBAL real* RESTRICT condensedPos) {
    for (int i = GLOBAL_ID; i < NUM_ACTIVE; i += GLOBAL_SIZE) {
        real3 pos = trimTo3(posq[activeParticles[i]]);
        int iOffset = 3 * i;
        condensedPos[iOffset] = pos.x;
        condensedPos[iOffset + 1] = pos.y;
        condensedPos[iOffset + 2] = pos.z;
    }
}

/**
 * Find a bounding box for the atoms in each block.
 */
KERNEL void findBlockBounds(real4 periodicBoxSize, real4 invPeriodicBoxSize,
        real4 periodicBoxVecX, real4 periodicBoxVecY, real4 periodicBoxVecZ,
        GLOBAL const real* RESTRICT condensedPos, GLOBAL real4* RESTRICT blockCenter,
        GLOBAL real4* RESTRICT blockBoundingBox, GLOBAL int* RESTRICT numNeighborPairs) {

    // Each thread (index) processes one block of WARP_SIZE atoms (from base
    // through min(base + WARP_SIZE, NUM_ACTIVE) - 1) at a time.

    int index = GLOBAL_ID;
    int base = index * WARP_SIZE;

    while (base < NUM_ACTIVE) {
        real3 pos3 = loadCondensedPos(condensedPos, base);
        real4 pos = make_real4(pos3.x, pos3.y, pos3.z, 0);
#ifdef USE_PERIODIC
        APPLY_PERIODIC_TO_POS(pos)
#endif
        real4 minPos = pos;
        real4 maxPos = pos;
        int last = min(base + WARP_SIZE, NUM_ACTIVE);
        for (int i = base + 1; i < last; i++) {
            pos3 = loadCondensedPos(condensedPos, i);
            pos = make_real4(pos3.x, pos3.y, pos3.z, 0);
#ifdef USE_PERIODIC
            real4 center = 0.5f * (maxPos + minPos);
            APPLY_PERIODIC_TO_POS_WITH_CENTER(pos, center)
#endif
            minPos = make_real4(min(minPos.x, pos.x), min(minPos.y, pos.y), min(minPos.z, pos.z), 0);
            maxPos = make_real4(max(maxPos.x, pos.x), max(maxPos.y, pos.y), max(maxPos.z, pos.z), 0);
        }

        real4 blockSize = 0.5f * (maxPos - minPos);
        blockBoundingBox[index] = blockSize;
        blockCenter[index] = 0.5f * (maxPos + minPos);

        // All threads can process GLOBAL_SIZE blocks together; if there are
        // more blocks, advance to the next block for this thread.

        index += GLOBAL_SIZE;
        base = index * WARP_SIZE;
    }

    // Reset numNeighborPairs for the findNeighbors kernel.

    if (GLOBAL_ID == 0) {
        *numNeighborPairs = 0;
    }
}

/**
 * Find a list of neighbors for each atom.
 */
KERNEL void findNeighbors(real4 periodicBoxSize, real4 invPeriodicBoxSize,
        real4 periodicBoxVecX, real4 periodicBoxVecY, real4 periodicBoxVecZ,
        GLOBAL const real* RESTRICT condensedPos, GLOBAL const real4* RESTRICT parameters,
        GLOBAL const real4* RESTRICT blockCenter, GLOBAL const real4* RESTRICT blockBoundingBox,
        GLOBAL int* RESTRICT numNeighborPairs, GLOBAL int* RESTRICT numNeighborsForAtom,
        GLOBAL int2* RESTRICT neighborPairs, real cutoffSquared, int maxNeighborPairs) {

    // parameters: each real4 holds [radius, p2, p3, p4].
    // posrCache:  each real4 holds [x, y, z, radius].

    LOCAL real4 posrCache[FIND_NEIGHBORS_THREAD_BLOCK_SIZE];
    int indexInWarp = LOCAL_ID % WARP_SIZE;
    int warpStart = LOCAL_ID - indexInWarp;
#if !(defined(__CUDA_ARCH__) || defined(USE_HIP))
    LOCAL bool includeBlockFlags[FIND_NEIGHBORS_THREAD_BLOCK_SIZE];
#endif

    // Each thread will search for the neighbors of one atom at a time.  Note
    // that all threads in a warp are processing atoms from the same block.

    for (int active1 = GLOBAL_ID; active1 < PADDED_NUM_ACTIVE; active1 += GLOBAL_SIZE) {

        int numNeighborsForAtom1 = 0;

        real3 pos1 = loadCondensedPos(condensedPos, active1);
        real radius1 = parameters[active1].x;

        int block1 = active1 / WARP_SIZE;
        real4 blockCenter1 = blockCenter[block1];
        real4 blockSize1 = blockBoundingBox[block1];

        // Loop over atom blocks to search for neighbors.  The threads in a warp
        // compare block1 against 32 other blocks in parallel.

        for (int block2Base = block1; block2Base < NUM_BLOCKS; block2Base += WARP_SIZE) {
            int block2 = block2Base + indexInWarp;
            bool includeBlock2 = (block2 < NUM_BLOCKS);
            if (includeBlock2) {
                real4 blockCenter2 = blockCenter[block2];
                real4 blockSize2 = blockBoundingBox[block2];
                real4 blockDelta = blockCenter1 - blockCenter2;
#ifdef USE_PERIODIC
                APPLY_PERIODIC_TO_DELTA(blockDelta)
#endif
                blockDelta.x = max((real) 0, fabs(blockDelta.x) - blockSize1.x - blockSize2.x);
                blockDelta.y = max((real) 0, fabs(blockDelta.y) - blockSize1.y - blockSize2.y);
                blockDelta.z = max((real) 0, fabs(blockDelta.z) - blockSize1.z - blockSize2.z);
                includeBlock2 &= (blockDelta.x * blockDelta.x + blockDelta.y * blockDelta.y + blockDelta.z * blockDelta.z < cutoffSquared);
            }

            // Each thread now loops over any blocks we identified as
            // potentially containing neighbors of atoms in block1, and checks
            // to see if they are actually neighbors of atom1.

#if defined(__CUDA_ARCH__) || defined(USE_HIP)
            int includeBlockFlags = BALLOT(includeBlock2);
            while (includeBlockFlags != 0) {
                int i = __ffs(includeBlockFlags) - 1;
                includeBlockFlags &= includeBlockFlags - 1;
                {
#else
            SYNC_WARPS;
            includeBlockFlags[LOCAL_ID] = includeBlock2;
            SYNC_WARPS;
            for (int i = 0; i < WARP_SIZE; i++) {
                if (includeBlockFlags[warpStart + i]) {
#endif
                    int block2 = block2Base + i;
                    int start = block2 * WARP_SIZE;
                    int active = start + indexInWarp;

                    int included[WARP_SIZE];
                    int numIncluded = 0;

                    // All threads in a warp will compare atoms with atoms in
                    // block2, so load their parameters into shared memory.

                    SYNC_WARPS;
                    real3 pos = loadCondensedPos(condensedPos, active);
                    posrCache[LOCAL_ID] = make_real4(pos.x, pos.y, pos.z, parameters[active].x);
                    SYNC_WARPS;

                    if (active1 < NUM_ACTIVE) {
                        for (int j = 0; j < WARP_SIZE; j++) {
                            int active2 = start + j;
                            real4 posr2 = posrCache[warpStart + j];
                            real3 pos2 = trimTo3(posr2);
                            real radius2 = posr2.w;

                            // Decide whether to include this atom pair in the neighbor list.

                            real4 delta12 = delta(pos2, pos1, periodicBoxSize, invPeriodicBoxSize, periodicBoxVecX, periodicBoxVecY, periodicBoxVecZ);
                            real pairCutoff = radius1 + radius2;
                            if (active2 > active1 && active2 < NUM_ACTIVE && delta12.w < pairCutoff * pairCutoff) {
                                int includedIndex = numIncluded++;
                                included[includedIndex] = active2;
                            }
                        }
                    }

                    // Store any neighbors of this atom that were found.

                    if (numIncluded) {
                        int baseIndex = ATOMIC_ADD(numNeighborPairs, numIncluded);
                        if (baseIndex + numIncluded <= maxNeighborPairs) {
                            for (int j = 0; j < numIncluded; j++) {
                                neighborPairs[baseIndex + j] = make_int2(active1, included[j]);
                            }
                        }
                        numNeighborsForAtom1 += numIncluded;
                    }
                }
            }
        }

        numNeighborsForAtom[active1] = numNeighborsForAtom1;
        SYNC_WARPS;
    }
}

/**
 * Sum the neighbor counts to compute the start position of each atom.  This
 * kernel is executed as a single thread block group.
 */
KERNEL void computeNeighborStartIndices(GLOBAL int* RESTRICT numNeighborPairsPointer,
        GLOBAL int* RESTRICT numNeighborsForAtom, GLOBAL int* RESTRICT neighborStartIndex,
        int maxNeighborPairs) {

    LOCAL unsigned int posBuffer[THREAD_BLOCK_SIZE];

    if (*numNeighborPairsPointer > maxNeighborPairs) {
        return; // There wasn't enough memory for the neighbor list.
    }

    unsigned int globalOffset = 0;
    for (unsigned int startAtom = 0; startAtom < NUM_ACTIVE; startAtom += LOCAL_SIZE) {
        // Load the neighbor counts into local memory.

        unsigned int globalIndex = startAtom + LOCAL_ID;
        posBuffer[LOCAL_ID] = (globalIndex < NUM_ACTIVE ? numNeighborsForAtom[globalIndex] : 0);
        SYNC_THREADS;

        // Perform a parallel prefix sum.

        for (unsigned int step = 1; step < LOCAL_SIZE; step *= 2) {
            unsigned int add = (LOCAL_ID >= step ? posBuffer[LOCAL_ID - step] : 0);
            SYNC_THREADS;
            posBuffer[LOCAL_ID] += add;
            SYNC_THREADS;
        }

        // Write the results back to global memory.

        if (globalIndex < NUM_ACTIVE) {
            neighborStartIndex[globalIndex + 1] = posBuffer[LOCAL_ID] + globalOffset;
            // Clear this so the next kernel can use it as a counter.
            numNeighborsForAtom[globalIndex] = 0;
        }
        globalOffset += posBuffer[LOCAL_SIZE - 1];
        SYNC_THREADS;
    }

    if (LOCAL_ID == 0) {
        neighborStartIndex[0] = 0;
    }
}

/**
 * Assemble the final neighbor list.
 */
KERNEL void copyPairsToNeighborList(real4 periodicBoxSize, real4 invPeriodicBoxSize,
        real4 periodicBoxVecX, real4 periodicBoxVecY, real4 periodicBoxVecZ,
        GLOBAL const real* RESTRICT condensedPos, GLOBAL const real4* RESTRICT parameters,
        GLOBAL int* RESTRICT numNeighborPairsPointer, GLOBAL int* RESTRICT numNeighborsForAtom,
        GLOBAL const int* RESTRICT neighborStartIndex, GLOBAL const int2* RESTRICT neighborPairs,
        GLOBAL int2* RESTRICT neighbors, GLOBAL NeighborData* RESTRICT neighborData, int maxNeighborPairs) {

    int numNeighborPairs = *numNeighborPairsPointer;
    if (numNeighborPairs > maxNeighborPairs) {
        return; // There wasn't enough memory for the neighbor list.
    }

    for (unsigned int pairIndex = GLOBAL_ID; pairIndex < numNeighborPairs; pairIndex += GLOBAL_SIZE) {
        int2 pair = neighborPairs[pairIndex];

        int offset = ATOMIC_ADD(numNeighborsForAtom + pair.x, 1);
        int storeIndex = neighborStartIndex[pair.x] + offset;

        // Store the original index so we can get back to (i, j) pairs.

        neighbors[storeIndex] = make_int2(pair.y, pairIndex);

        // Precompute data for this pair to be looked up later.

        real3 iPos = loadCondensedPos(condensedPos, pair.x);
        real3 jPos = loadCondensedPos(condensedPos, pair.y);
        real4 ijDelta = delta(jPos, iPos, periodicBoxSize, invPeriodicBoxSize, periodicBoxVecX, periodicBoxVecY, periodicBoxVecZ);
        real iRadius = parameters[pair.x].x;
        real jRadius = parameters[pair.y].x;
        real iRadiusPi = PI * iRadius;
        real jRadiusPi = PI * jRadius;
        real r = SQRT(ijDelta.w);
        real rRecip = (real) 1 / r;
        real deltaRadiusR = (iRadius * iRadius - jRadius * jRadius) * rRecip;
        real deltaRadiusRSq = deltaRadiusR * rRecip;
        real3 forceDir = trimTo3(ijDelta) * rRecip;
        real3 iForceDir = iRadiusPi * ((real) -1 + deltaRadiusRSq) * forceDir;
        real3 jForceDir = jRadiusPi * ((real) -1 - deltaRadiusRSq) * forceDir;

        neighborData[storeIndex].data12 = make_real4(iForceDir.x, iForceDir.y, iForceDir.z, iRadiusPi * ((real) 2 * iRadius - r - deltaRadiusR));
        neighborData[storeIndex].data21 = make_real4(jForceDir.x, jForceDir.y, jForceDir.z, jRadiusPi * ((real) 2 * jRadius - r + deltaRadiusR));
    }
}

/**
 * Compute the LCPO interaction.
 */
KERNEL void computeInteraction(
        real4 periodicBoxSize, real4 invPeriodicBoxSize,
        real4 periodicBoxVecX, real4 periodicBoxVecY, real4 periodicBoxVecZ,
        GLOBAL const real* RESTRICT condensedPos, GLOBAL mm_ulong* RESTRICT forceBuffers,
        GLOBAL mixed* RESTRICT energyBuffer, GLOBAL const int* RESTRICT activeParticles,
        GLOBAL const real4* RESTRICT parameters, GLOBAL int* RESTRICT numNeighborPairsPointer,
        GLOBAL const int* RESTRICT neighborStartIndex, GLOBAL const int2* RESTRICT neighborPairs,
        GLOBAL int2* RESTRICT neighbors, GLOBAL NeighborData* RESTRICT neighborData, int maxNeighborPairs) {

    int numNeighborPairs = *numNeighborPairsPointer;
    if (numNeighborPairs > maxNeighborPairs) {
        return; // There wasn't enough memory for the neighbor list.
    }

    mixed energy = 0;
    for (unsigned int pairIndex = GLOBAL_ID; pairIndex < numNeighborPairs; pairIndex += GLOBAL_SIZE) {
        int2 pair = neighbors[pairIndex];
        int i = neighborPairs[pair.y].x;
        int j = pair.x;
        real3 iPos = loadCondensedPos(condensedPos, i);
        real4 iParams = parameters[i];
        real4 jParams = parameters[j];
        real3 iForce = make_real3(0);
        real3 jForce = make_real3(0);
        NeighborData ijData = neighborData[pairIndex];
        real4 Aij = ijData.data12;
        real4 Aji = ijData.data21;

        real4 ijForce2 = iParams.y * Aij + jParams.y * Aji;
        real3 ijForce2_3 = trimTo3(ijForce2);
        iForce += ijForce2_3;
        jForce -= ijForce2_3;
        energy += ijForce2.w;

        real iRadius = iParams.x;
        real iRadiusPi = PI * iRadius;

        int jStart = neighborStartIndex[j];
        int jEnd = neighborStartIndex[j + 1];
        for (int jNeighbor = jStart; jNeighbor < jEnd; jNeighbor++) {
            int k = neighbors[jNeighbor].x;

            // It is faster to recompute the i-k interaction distance and
            // parameters to decide whether or not to include k than to try to
            // look up its index in the neighbor list of i.

            real3 kPos = loadCondensedPos(condensedPos, k);
            real4 kParams = parameters[k];
            real kRadius = kParams.x;
            real kRadiusPi = PI * kRadius;

            real4 ikDelta = delta(kPos, iPos, periodicBoxSize, invPeriodicBoxSize, periodicBoxVecX, periodicBoxVecY, periodicBoxVecZ);
            real pairCutoff = iRadius + kRadius;
            if (ikDelta.w >= pairCutoff * pairCutoff) {
                continue;
            }

            real r = SQRT(ikDelta.w);
            real rRecip = (real) 1 / r;
            real deltaRadiusR = (iRadius * iRadius - kRadius * kRadius) * rRecip;
            real deltaRadiusRSq = deltaRadiusR * rRecip;
            real3 forceDir = trimTo3(ikDelta) * rRecip;
            real3 iForceDir = iRadiusPi * ((real) -1 + deltaRadiusRSq) * forceDir;
            real3 jForceDir = kRadiusPi * ((real) -1 - deltaRadiusRSq) * forceDir;

            real4 Aik = make_real4(iForceDir.x, iForceDir.y, iForceDir.z, iRadiusPi * ((real) 2 * iRadius - r - deltaRadiusR));
            real4 Aki = make_real4(jForceDir.x, jForceDir.y, jForceDir.z, kRadiusPi * ((real) 2 * kRadius - r + deltaRadiusR));

            NeighborData jkData = neighborData[jNeighbor];
            real4 Ajk = jkData.data12;
            real4 Akj = jkData.data21;

            real4 ijkForce34 = (iParams.z + iParams.w * Aij.w) * Ajk + (iParams.z + iParams.w * Aik.w) * Akj;
            real4 jikForce34 = (jParams.z + jParams.w * Aji.w) * Aik + (jParams.z + jParams.w * Ajk.w) * Aki;
            real4 kijForce34 = (kParams.z + kParams.w * Aki.w) * Aij + (kParams.z + kParams.w * Akj.w) * Aji;
            real3 ijkForce34_3 = trimTo3(ijkForce34);
            real3 jikForce34_3 = trimTo3(jikForce34);
            real3 kijForce34_3 = trimTo3(kijForce34);

            real3 ijkForce4 = iParams.w * trimTo3(Aij) * Ajk.w;
            real3 ikjForce4 = iParams.w * trimTo3(Aik) * Akj.w;
            real3 jikForce4 = jParams.w * trimTo3(Aji) * Aik.w;
            real3 jkiForce4 = jParams.w * trimTo3(Ajk) * Aki.w;
            real3 kijForce4 = kParams.w * trimTo3(Aki) * Aij.w;
            real3 kjiForce4 = kParams.w * trimTo3(Akj) * Aji.w;

            energy += ijkForce34.w + jikForce34.w + kijForce34.w;
            iForce += jikForce34_3 + kijForce34_3 + ijkForce4 + ikjForce4 + jikForce4 + kijForce4;
            jForce += ijkForce34_3 - kijForce34_3 + jkiForce4 - jikForce4 - ijkForce4 + kjiForce4;
            real3 kForce = -(ijkForce34_3 + jikForce34_3 + kijForce4 + kjiForce4 + ikjForce4 + jkiForce4);

            storeForce(activeParticles[k], kForce, forceBuffers);
        }

        storeForce(activeParticles[i], iForce, forceBuffers);
        storeForce(activeParticles[j], jForce, forceBuffers);
    }
    energyBuffer[GLOBAL_ID] += energy;
}
