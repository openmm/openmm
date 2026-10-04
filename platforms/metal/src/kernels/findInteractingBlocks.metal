#define GROUP_SIZE 256
#define BUFFER_SIZE 256

/**
 * Find a bounding box for the atoms in each block.
 */
KERNEL void findBlockBounds(int numAtoms, real4 periodicBoxSize, real4 invPeriodicBoxSize, real4 periodicBoxVecX, real4 periodicBoxVecY, real4 periodicBoxVecZ,
        GLOBAL const real4* RESTRICT posq, GLOBAL real4* RESTRICT blockCenter, GLOBAL real4* RESTRICT blockBoundingBox, GLOBAL int* RESTRICT rebuildNeighborList,
        GLOBAL real2* RESTRICT blockSizeRange) {
    int index = GLOBAL_ID;
    int base = index*TILE_SIZE;
    real minSize = 1e38, maxSize = 0;
    while (base < numAtoms) {
        real4 pos = posq[base];
#ifdef USE_PERIODIC
        APPLY_PERIODIC_TO_POS(pos)
#endif
        real4 minPos = pos;
        real4 maxPos = pos;
        int last = min(base+TILE_SIZE, numAtoms);
        for (int i = base+1; i < last; i++) {
            pos = posq[i];
#ifdef USE_PERIODIC
            real4 center = 0.5f*(maxPos+minPos);
            APPLY_PERIODIC_TO_POS_WITH_CENTER(pos, center)
#endif
            minPos = make_real4(min(minPos.x,pos.x), min(minPos.y,pos.y), min(minPos.z,pos.z), 0);
            maxPos = make_real4(max(maxPos.x,pos.x), max(maxPos.y,pos.y), max(maxPos.z,pos.z), 0);
        }
        real4 blockSize = 0.5f*(maxPos-minPos);
        real4 center = 0.5f*(maxPos+minPos);
        center.w = 0;
        for (int i = base; i < last; i++) {
            pos = posq[i];
            real4 delta = posq[i]-center;
#ifdef USE_PERIODIC
            APPLY_PERIODIC_TO_DELTA(delta)
#endif
            center.w = max(center.w, delta.x*delta.x+delta.y*delta.y+delta.z*delta.z);
        }
        center.w = sqrt(center.w);
        blockBoundingBox[index] = blockSize;
        blockCenter[index] = center;
        real totalSize = blockSize.x+blockSize.y+blockSize.z;
        minSize = min(minSize, totalSize);
        maxSize = max(maxSize, totalSize);
        index += GLOBAL_SIZE;
        base = index*TILE_SIZE;
    }

    // Record the range of sizes seen by threads in this block.

    LOCAL real minBuffer[64], maxBuffer[64];
    minBuffer[LOCAL_ID] = minSize;
    maxBuffer[LOCAL_ID] = maxSize;
    SYNC_THREADS
    for (int step = 1; step < 64; step *= 2) {
        if (LOCAL_ID+step < 64 && LOCAL_ID%(2*step) == 0) {
            minBuffer[LOCAL_ID] = min(minBuffer[LOCAL_ID], minBuffer[LOCAL_ID+step]);
            maxBuffer[LOCAL_ID] = max(maxBuffer[LOCAL_ID], maxBuffer[LOCAL_ID+step]);
        }
        SYNC_THREADS
    }
    if (LOCAL_ID == 0)
        blockSizeRange[GROUP_ID] = make_real2(minBuffer[0], maxBuffer[0]);
    if (GROUP_ID == 0 && LOCAL_ID == 0)
        rebuildNeighborList[0] = 0;
}

KERNEL void computeSortKeys(GLOBAL const real4* RESTRICT blockBoundingBox, GLOBAL unsigned int* RESTRICT sortedBlocks, GLOBAL real2* RESTRICT blockSizeRange, int numSizes) {
    // Find the total range of sizes recorded by all blocks.

    LOCAL real2 sizeRange;
    if (LOCAL_ID == 0) {
        sizeRange = blockSizeRange[0];
        for (int i = 1; i < numSizes; i++) {
            real2 size = blockSizeRange[i];
            if (size.x > 0)
                sizeRange.x = min(sizeRange.x, size.x);
            sizeRange.y = max(sizeRange.y, size.y);
        }
        sizeRange.x = LOG(sizeRange.x);
        sizeRange.y = LOG(sizeRange.y);
    }
    SYNC_THREADS

    // Sort keys store the bin in the high order part and the block in the low
    // order part.

    int numSizeBins = 20;
    real scale = numSizeBins/(sizeRange.y-sizeRange.x);
    for (unsigned int i = GLOBAL_ID; i < NUM_BLOCKS; i += GLOBAL_SIZE) {
        real4 box = blockBoundingBox[i];
        real size = LOG(box.x+box.y+box.z);
        int bin = (size-sizeRange.x)*scale;
        bin = max(0, min(bin, numSizeBins-1));
        sortedBlocks[i] = (((unsigned int) bin)<<BIN_SHIFT) + i;
    }
}

/**
 * Sort the data about bounding boxes so it can be accessed more efficiently in the next kernel.
 */
KERNEL void sortBoxData(GLOBAL const unsigned int* RESTRICT sortedBlocks, GLOBAL const real4* RESTRICT blockCenter,
        GLOBAL const real4* RESTRICT blockBoundingBox, GLOBAL real4* RESTRICT sortedBlockCenter, GLOBAL real3* RESTRICT sortedBlockBoundingBox,
#ifdef USE_LARGE_BLOCKS
        GLOBAL real4* RESTRICT largeBlockCenter, GLOBAL real3* RESTRICT largeBlockBoundingBox, real4 periodicBoxSize,
        real4 invPeriodicBoxSize, real4 periodicBoxVecX, real4 periodicBoxVecY, real4 periodicBoxVecZ,
#endif
        GLOBAL const real4* RESTRICT posq, GLOBAL const real4* RESTRICT oldPositions,
        GLOBAL unsigned int* RESTRICT interactionCount, GLOBAL int* RESTRICT rebuildNeighborList, bool forceRebuild) {
    for (int i = GLOBAL_ID; i < NUM_BLOCKS; i += GLOBAL_SIZE) {
        unsigned int index = sortedBlocks[i] & BLOCK_INDEX_MASK;
        sortedBlockCenter[i] = blockCenter[index];
        sortedBlockBoundingBox[i] = real3(trimTo3(blockBoundingBox[index]));

#ifdef USE_LARGE_BLOCKS
        // Compute the sizes of large blocks (composed of 32 regular blocks) starting from each block.

        real4 minPos = blockCenter[index]-blockBoundingBox[index];
        real4 maxPos = blockCenter[index]+blockBoundingBox[index];
        int last = min(i+32, NUM_BLOCKS);
        for (int j = i+1; j < last; j++) {
            unsigned int index2 = sortedBlocks[j] & BLOCK_INDEX_MASK;
            real4 blockPos = blockCenter[index2];
            real4 width = blockBoundingBox[index2];
#ifdef USE_PERIODIC
            real4 center = 0.5f*(maxPos+minPos);
            APPLY_PERIODIC_TO_POS_WITH_CENTER(blockPos, center)
#endif
            minPos = make_real4(min(minPos.x, blockPos.x-width.x), min(minPos.y, blockPos.y-width.y), min(minPos.z, blockPos.z-width.z), 0);
            maxPos = make_real4(max(maxPos.x, blockPos.x+width.x), max(maxPos.y, blockPos.y+width.y), max(maxPos.z, blockPos.z+width.z), 0);
        }
        largeBlockCenter[i] = 0.5f*(maxPos+minPos);
        largeBlockBoundingBox[i] = real3(trimTo3(0.5f*(maxPos-minPos)));
#endif
    }

    // Also check whether any atom has moved enough so that we really need to rebuild the neighbor list.

    bool rebuild = forceRebuild;
    for (int i = GLOBAL_ID; i < NUM_ATOMS; i += GLOBAL_SIZE) {
        real4 delta = oldPositions[i]-posq[i];
        if (delta.x*delta.x + delta.y*delta.y + delta.z*delta.z > 0.25f*PADDING*PADDING)
            rebuild = true;
    }
    if (rebuild) {
        rebuildNeighborList[0] = 1;
        interactionCount[0] = 0;
        interactionCount[1] = 0;
    }
}

DEVICE int saveSinglePairs(int x, LOCAL int* atoms, LOCAL int* flags, int length, unsigned int maxSinglePairs,
        GLOBAL unsigned int* singlePairCount, GLOBAL int2* singlePairs, LOCAL int* sumBuffer, LOCAL volatile unsigned int& pairStartIndex) {
    // Record interactions that should be computed as single pairs rather than in blocks.

    const int indexInWarp = LOCAL_ID%32;
    int sum = 0;
    #pragma unroll 8 // (GROUP_SIZE / TILE_SIZE)
    for (int i = indexInWarp; i < length; i += 32) {
        int count = popcount(flags[i]);
        sum += (count <= MAX_BITS_FOR_PAIRS ? count : 0);
    }
    for (int i = 1; i < 32; i *= 2) {
        int n = simd_shuffle_up(sum, i);
        if (indexInWarp >= i)
            sum += n;
    }
    if (indexInWarp == 31)
        pairStartIndex = ATOMIC_ADD(singlePairCount,(unsigned int) sum);
    SYNC_WARPS
    int prevSum = simd_shuffle_up(sum, 1);
    unsigned int pairIndex = pairStartIndex + (indexInWarp > 0 ? prevSum : 0);
    for (int i = indexInWarp; i < length; i += 32) {
        int count = popcount(flags[i]);
        if (count <= MAX_BITS_FOR_PAIRS && pairIndex+count <= maxSinglePairs) {
            int f = flags[i];
            while (f != 0) {
                singlePairs[pairIndex] = make_int2(atoms[i], x*TILE_SIZE+ctz(f)-1);
                f &= f-1;
                pairIndex++;
            }
        }
    }

    // Compact the remaining interactions.

    const int warpMask = (1<<indexInWarp)-1;
    int numCompacted = 0;
    for (int start = 0; start < length; start += 32) {
        int i = start+indexInWarp;
        int atom = atoms[i];
        int flag = flags[i];
        bool include = (i < length && popcount(flags[i]) > MAX_BITS_FOR_PAIRS);
        int includeFlags = (uint64_t) BALLOT(include);
        if (include) {
            int index = numCompacted+popcount(includeFlags&warpMask);
            atoms[index] = atom;
            flags[index] = flag;
        }
        numCompacted += popcount(includeFlags);
    }
    return numCompacted;
}

/**
 * Compare the bounding boxes for each pair of atom blocks (comprised of 32 atoms each), forming a tile. If the two
 * atom blocks are sufficiently far apart, mark them as non-interacting. There are two stages in the algorithm.
 *
 * STAGE 1:
 *
 * A coarse grained atom block against interacting atom block neighbour list is constructed.
 *
 * Each warp first loads in some block X of interest. Each thread within the warp then loads
 * in a different atom block Y. If Y has exclusions with X, then Y is not processed.  If the bounding boxes
 * of the two atom blocks are within the cutoff distance, then the two atom blocks are considered to be
 * interacting and Y is added to the buffer for X.
 *
 * STAGE 2:
 *
 * A fine grained atom block against interacting atoms neighbour list is constructed.
 *
 * The warp loops over atom blocks Y that were found to (possibly) interact with atom block X.  Each thread
 * in the warp loops over the 32 atoms in X and compares their positions to one particular atom from block Y.
 * If it finds one closer than the cutoff distance, the atom is added to the list of atoms interacting with block X.
 * This continues until the buffer fills up, at which point the results are written to global memory.
 *
 * [in] periodicBoxSize        - size of the rectangular periodic box
 * [in] invPeriodicBoxSize     - inverse of the periodic box
 * [in] blockCenter            - the center of each bounding box
 * [in] blockBoundingBox       - bounding box of each atom block
 * [out] interactionCount      - total number of tiles that have interactions
 * [out] interactingTiles      - set of blocks that have interactions
 * [out] interactingAtoms      - a list of atoms that interact with each atom block
 * [in] posq                   - x,y,z coordinates of each atom and charge q
 * [in] maxTiles               - maximum number of tiles to process, used for multi-GPUs
 * [in] startBlockIndex        - first block to process, used for multi-GPUs,
 * [in] numBlocks              - total number of atom blocks
 * [in] sortedBlocks           - a sorted list of atom blocks based on volume
 * [in] sortedBlockCenter      - sorted centers, duplicated for fast access to avoid indexing
 * [in] sortedBlockBoundingBox - sorted bounding boxes, duplicated for fast access
 * [in] exclusionIndices       - maps into exclusionRowIndices with the starting position for a given atom
 * [in] exclusionRowIndices    - stores the a continuous list of exclusions
 *           eg: block 0 is excluded from atom 3,5,6
 *               block 1 is excluded from atom 3,4
 *               block 2 is excluded from atom 1,3,5,6
 *              exclusionIndices[0][3][5][8]
 *           exclusionRowIndices[3][5][6][3][4][1][3][5][6]
 *                         index 0  1  2  3  4  5  6  7  8
 * [out] oldPos                - stores the positions of the atoms in which this neighbourlist was built on
 *                             - this is used to decide when to rebuild a neighbourlist
 * [in] rebuildNeighbourList   - whether or not to execute this kernel
 *
 */
KERNEL void findBlocksWithInteractions(real4 periodicBoxSize, real4 invPeriodicBoxSize, real4 periodicBoxVecX, real4 periodicBoxVecY, real4 periodicBoxVecZ,
        GLOBAL unsigned int* RESTRICT interactionCount, GLOBAL int* RESTRICT interactingTiles, GLOBAL unsigned int* RESTRICT interactingAtoms,
        GLOBAL int2* RESTRICT singlePairs, GLOBAL const real4* RESTRICT posq, unsigned int maxTiles, unsigned int maxSinglePairs, unsigned int startBlockIndex,
        unsigned int numBlocks, GLOBAL unsigned int* RESTRICT sortedBlocks, GLOBAL const real4* RESTRICT sortedBlockCenter, GLOBAL const real3* RESTRICT sortedBlockBoundingBox,
#ifdef USE_LARGE_BLOCKS
        GLOBAL real4* RESTRICT largeBlockCenter, GLOBAL real3* RESTRICT largeBlockBoundingBox,
#endif
        GLOBAL const unsigned int* RESTRICT exclusionIndices, GLOBAL const unsigned int* RESTRICT exclusionRowIndices,
        GLOBAL real4* RESTRICT oldPositions, GLOBAL const int* RESTRICT rebuildNeighborList) {

    if (rebuildNeighborList[0] == 0)
        return; // The neighbor list doesn't need to be rebuilt.

    const int indexInWarp = LOCAL_ID%32;
    const int warpStart = LOCAL_ID-indexInWarp;
    const int totalWarps = GLOBAL_SIZE/32;
    const int warpIndex = (GLOBAL_ID)/32;
    const int warpMask = (1<<indexInWarp)-1;
    LOCAL int workgroupBuffer[BUFFER_SIZE*(GROUP_SIZE/32)];
    LOCAL int workgroupFlagsBuffer[BUFFER_SIZE*(GROUP_SIZE/32)];
    LOCAL int warpExclusions[MAX_EXCLUSIONS*(GROUP_SIZE/32)];
    LOCAL real4 posBuffer[GROUP_SIZE];
    LOCAL volatile unsigned int workgroupTileIndex[GROUP_SIZE/32];
    LOCAL unsigned int workgroupPairStartIndex[GROUP_SIZE/32];
    LOCAL int* sumBuffer = (LOCAL int*) posBuffer; // Reuse the same buffer to save memory
    LOCAL int* buffer = workgroupBuffer+BUFFER_SIZE*(warpStart/32);
    LOCAL int* flagsBuffer = workgroupFlagsBuffer+BUFFER_SIZE*(warpStart/32);
    LOCAL int* exclusionsForX = warpExclusions+MAX_EXCLUSIONS*(warpStart/32);
    LOCAL volatile unsigned int& tileStartIndex = workgroupTileIndex[warpStart/32];
    LOCAL volatile unsigned int& pairStartIndex = workgroupPairStartIndex[warpStart/32];

    // Loop over blocks.

    for (int block1 = startBlockIndex+warpIndex; block1 < startBlockIndex+numBlocks; block1 += totalWarps) {
        // Load data for this block.  Note that all threads in a warp are processing the same block.

        int x = sortedBlocks[block1] & BLOCK_INDEX_MASK;
        real4 blockCenterX = sortedBlockCenter[block1];
        real3 blockSizeX = sortedBlockBoundingBox[block1];
        int neighborsInBuffer = 0;
        real4 pos1 = posq[x*TILE_SIZE+indexInWarp];
#ifdef USE_PERIODIC
        const bool singlePeriodicCopy = (0.5f*periodicBoxSize.x-blockSizeX.x >= PADDED_CUTOFF &&
                                         0.5f*periodicBoxSize.y-blockSizeX.y >= PADDED_CUTOFF &&
                                         0.5f*periodicBoxSize.z-blockSizeX.z >= PADDED_CUTOFF);
        if (singlePeriodicCopy) {
            // The box is small enough that we can just translate all the atoms into a single periodic
            // box, then skip having to apply periodic boundary conditions later.

            APPLY_PERIODIC_TO_POS_WITH_CENTER(pos1, blockCenterX)
        }
#endif
        pos1.w = 0.5f * (pos1.x * pos1.x + pos1.y * pos1.y + pos1.z * pos1.z);
        posBuffer[LOCAL_ID] = pos1;

        // Load exclusion data for block x.

        const int exclusionStart = exclusionRowIndices[x];
        const int exclusionEnd = exclusionRowIndices[x+1];
        const int numExclusions = exclusionEnd-exclusionStart;
        #pragma unroll 4 // (MAX_EXCLUSIONS)
        for (int j = indexInWarp; j < numExclusions; j += 32)
            exclusionsForX[j] = exclusionIndices[exclusionStart+j];
        if (MAX_EXCLUSIONS > 32)
            SYNC_THREADS

        // Loop over atom blocks to search for neighbors.  The threads in a warp compare block1 against 32
        // other blocks in parallel.

#ifdef USE_LARGE_BLOCKS
        int largeBlockFlags = 0;
        int loadedLargeBlocks = 0;
#endif
        for (int block2Base = block1+1; block2Base < NUM_BLOCKS; block2Base += 32) {
#ifdef USE_LARGE_BLOCKS
            if (loadedLargeBlocks == 0) {
                // Check the next set of large blocks.

                int largeBlockIndex = block2Base + 32*indexInWarp;
                bool includeLargeBlock = false;
                if (largeBlockIndex < NUM_BLOCKS) {
                    real4 largeCenter = largeBlockCenter[largeBlockIndex];
                    real3 largeSize = largeBlockBoundingBox[largeBlockIndex];
                    real4 blockDelta = blockCenterX-largeCenter;
#ifdef USE_PERIODIC
                    APPLY_PERIODIC_TO_DELTA(blockDelta)
#endif
                    blockDelta.x = max(0.0f, fabs(blockDelta.x)-blockSizeX.x-largeSize.x);
                    blockDelta.y = max(0.0f, fabs(blockDelta.y)-blockSizeX.y-largeSize.y);
                    blockDelta.z = max(0.0f, fabs(blockDelta.z)-blockSizeX.z-largeSize.z);
                    includeLargeBlock = (blockDelta.x*blockDelta.x+blockDelta.y*blockDelta.y+blockDelta.z*blockDelta.z < PADDED_CUTOFF_SQUARED);
#ifdef TRICLINIC
                    // The calculation to find the nearest periodic copy is only guaranteed to work if the nearest copy is less than half a box width away.
                    // If there's any possibility we might have missed it, do a detailed check.

                    if (periodicBoxSize.z/2-blockSizeX.z-largeSize.z < PADDED_CUTOFF || periodicBoxSize.y/2-blockSizeX.y-largeSize.y < PADDED_CUTOFF)
                        includeLargeBlock = true;
#endif
                }
                largeBlockFlags = (uint64_t) BALLOT(includeLargeBlock);
                loadedLargeBlocks = 32;
            }
            loadedLargeBlocks--;
            if ((largeBlockFlags&1) == 0) {
                // None of the next 32 blocks interact with block 1.

                largeBlockFlags >>= 1;
                continue;
            }
            largeBlockFlags >>= 1;
#endif
            int block2 = block2Base+indexInWarp;
            bool includeBlock2 = (block2 < NUM_BLOCKS);
            bool forceInclude = false;
            if (includeBlock2) {
                real4 blockCenterY = sortedBlockCenter[block2];
                real3 blockSizeY = sortedBlockBoundingBox[block2];
                real4 blockDelta = blockCenterX-blockCenterY;
#ifdef USE_PERIODIC
                APPLY_PERIODIC_TO_DELTA(blockDelta)
#endif
                includeBlock2 &= (blockDelta.x*blockDelta.x+blockDelta.y*blockDelta.y+blockDelta.z*blockDelta.z < (PADDED_CUTOFF+blockCenterX.w+blockCenterY.w)*(PADDED_CUTOFF+blockCenterX.w+blockCenterY.w));
                blockDelta.x = max(0.0f, fabs(blockDelta.x)-blockSizeX.x-blockSizeY.x);
                blockDelta.y = max(0.0f, fabs(blockDelta.y)-blockSizeX.y-blockSizeY.y);
                blockDelta.z = max(0.0f, fabs(blockDelta.z)-blockSizeX.z-blockSizeY.z);
                includeBlock2 &= (blockDelta.x*blockDelta.x+blockDelta.y*blockDelta.y+blockDelta.z*blockDelta.z < PADDED_CUTOFF_SQUARED);
#ifdef TRICLINIC
                // The calculation to find the nearest periodic copy is only guaranteed to work if the nearest copy is less than half a box width away.
                // If there's any possibility we might have missed it, do a detailed check.

                if (periodicBoxSize.z/2-blockSizeX.z-blockSizeY.z < PADDED_CUTOFF || periodicBoxSize.y/2-blockSizeX.y-blockSizeY.y < PADDED_CUTOFF)
                    includeBlock2 = forceInclude = true;
#endif
                if (includeBlock2) {
                    int y = sortedBlocks[block2] & BLOCK_INDEX_MASK;
                    #pragma unroll 4 // (MAX_EXCLUSIONS)
                    for (int k = 0; k < numExclusions; k++)
                        includeBlock2 &= (exclusionsForX[k] != y);
                }
            }

            // Loop over any blocks we identified as potentially containing neighbors.

            int includeBlockFlags = (uint64_t) BALLOT(includeBlock2);
            int forceIncludeFlags = (uint64_t) BALLOT(forceInclude);
            while (includeBlockFlags != 0) {
                int i = ctz(includeBlockFlags)-1;
                includeBlockFlags &= includeBlockFlags-1;
                forceInclude = (forceIncludeFlags>>i) & 1;
                int y = sortedBlocks[block2Base+i] & BLOCK_INDEX_MASK;

                // Check each atom in block Y for interactions.

                int atom2 = y*TILE_SIZE+indexInWarp;
                real4 pos2 = posq[atom2];
#ifdef USE_PERIODIC
                if (singlePeriodicCopy) {
                    APPLY_PERIODIC_TO_POS_WITH_CENTER(pos2, blockCenterX)
                }
#endif
                pos2.w = 0.5f * (pos2.x * pos2.x + pos2.y * pos2.y + pos2.z * pos2.z);

                real4 blockCenterY = sortedBlockCenter[block2Base+i];
                real3 atomDelta = trimTo3(posBuffer[warpStart+indexInWarp])-trimTo3(blockCenterY);
#ifdef USE_PERIODIC
                APPLY_PERIODIC_TO_DELTA(atomDelta)
#endif
                int atomFlags = (uint64_t) BALLOT(forceInclude || atomDelta.x*atomDelta.x+atomDelta.y*atomDelta.y+atomDelta.z*atomDelta.z < (PADDED_CUTOFF+blockCenterY.w)*(PADDED_CUTOFF+blockCenterY.w));
                int interacts = 0;
                if (atom2 < NUM_ATOMS && atomFlags != 0) {
#ifdef USE_PERIODIC
                    if (!singlePeriodicCopy) {
                        int first = ctz(atomFlags)-1;
                        int last = 32-clz(atomFlags);
                        for (int j = first; j < last; j++) {
                            real3 delta = trimTo3(pos2)-trimTo3(posBuffer[warpStart+j]);
                            APPLY_PERIODIC_TO_DELTA(delta)
                            interacts |= (delta.x*delta.x+delta.y*delta.y+delta.z*delta.z < PADDED_CUTOFF_SQUARED ? 1<<j : 0);
                        }
                    }
                    else {
#endif
                        #pragma unroll
                        for (int j = 0; j < 32; j++) {
                            real4 posj = posBuffer[warpStart+j];
                            real halfDist2 = posj.w + pos2.w - posj.x*pos2.x - posj.y*pos2.y - posj.z*pos2.z;
                            interacts |= (halfDist2 < 0.5f * PADDED_CUTOFF_SQUARED ? 1<<j : 0);
                        }
#ifdef USE_PERIODIC
                    }
#endif
                }

                // Add any interacting atoms to the buffer.

                int includeAtomFlags = (uint64_t) BALLOT(interacts);
                if (interacts) {
                    int index = neighborsInBuffer+popcount(includeAtomFlags&warpMask);
                    buffer[index] = atom2;
                    flagsBuffer[index] = interacts;
                }
                neighborsInBuffer += popcount(includeAtomFlags);
                if (neighborsInBuffer > BUFFER_SIZE-TILE_SIZE) {
                    // Store the new tiles to memory.

#if MAX_BITS_FOR_PAIRS > 0
                    neighborsInBuffer = saveSinglePairs(x, buffer, flagsBuffer, neighborsInBuffer, maxSinglePairs, &interactionCount[1], singlePairs, sumBuffer+warpStart, pairStartIndex);
#endif
                    unsigned int tilesToStore = neighborsInBuffer/TILE_SIZE;
                    if (tilesToStore > 0) {
                        if (indexInWarp == 0)
                            tileStartIndex = ATOMIC_ADD(&interactionCount[0], tilesToStore);
                        unsigned int newTileStartIndex = tileStartIndex;
                        if (newTileStartIndex+tilesToStore <= maxTiles) {
                            if (indexInWarp < tilesToStore)
                                interactingTiles[newTileStartIndex+indexInWarp] = x;
                            #pragma unroll 8 // (GROUP_SIZE / TILE_SIZE)
                            for (int j = 0; j < tilesToStore; j++)
                                interactingAtoms[(newTileStartIndex+j)*TILE_SIZE+indexInWarp] = buffer[indexInWarp+j*TILE_SIZE];
                        }
                        if (indexInWarp+TILE_SIZE*tilesToStore < BUFFER_SIZE)
                            buffer[indexInWarp] = buffer[indexInWarp+TILE_SIZE*tilesToStore];
                        neighborsInBuffer -= TILE_SIZE*tilesToStore;
                    }
                }
            }
        }

        // If we have a partially filled buffer,  store it to memory.

#if MAX_BITS_FOR_PAIRS > 0
        if (neighborsInBuffer > 32)
            neighborsInBuffer = saveSinglePairs(x, buffer, flagsBuffer, neighborsInBuffer, maxSinglePairs, &interactionCount[1], singlePairs, sumBuffer+warpStart, pairStartIndex);
#endif
        if (neighborsInBuffer > 0) {
            unsigned int tilesToStore = (neighborsInBuffer+TILE_SIZE-1)/TILE_SIZE;
            if (indexInWarp == 0)
                tileStartIndex = ATOMIC_ADD(&interactionCount[0], tilesToStore);
            unsigned int newTileStartIndex = tileStartIndex;
            if (newTileStartIndex+tilesToStore <= maxTiles) {
                if (indexInWarp < tilesToStore)
                    interactingTiles[newTileStartIndex+indexInWarp] = x;
                #pragma unroll 8 // (GROUP_SIZE / TILE_SIZE)
                for (int j = 0; j < tilesToStore; j++)
                    interactingAtoms[(newTileStartIndex+j)*TILE_SIZE+indexInWarp] = (indexInWarp+j*TILE_SIZE < neighborsInBuffer ? buffer[indexInWarp+j*TILE_SIZE] : NUM_ATOMS);
            }
        }
    }

    // Record the positions the neighbor list is based on.

    for (int i = GLOBAL_ID; i < NUM_ATOMS; i += GLOBAL_SIZE)
        oldPositions[i] = posq[i];
}
