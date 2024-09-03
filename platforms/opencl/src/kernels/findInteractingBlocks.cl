#pragma OPENCL EXTENSION cl_khr_global_int32_base_atomics : enable
#pragma OPENCL EXTENSION cl_khr_byte_addressable_store : enable

/**
 * Find a bounding box for the atoms in each block.
 */
__kernel void findBlockBounds(int numAtoms, real4 periodicBoxSize, real4 invPeriodicBoxSize, real4 periodicBoxVecX, real4 periodicBoxVecY, real4 periodicBoxVecZ,
        __global const real4* restrict posq, __global real4* restrict blockCenter, __global real4* restrict blockBoundingBox, __global int* restrict rebuildNeighborList,
        __global real2* restrict blockSizeRange) {
    int index = get_global_id(0);
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
            minPos = min(minPos, pos);
            maxPos = max(maxPos, pos);
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
        index += get_global_size(0);
        base = index*TILE_SIZE;
    }

    // Record the range of sizes seen by threads in this block.

    __local real minBuffer[64], maxBuffer[64];
    minBuffer[get_local_id(0)] = minSize;
    maxBuffer[get_local_id(0)] = maxSize;
    barrier(CLK_LOCAL_MEM_FENCE);
    for (int step = 1; step < 64; step *= 2) {
        if (get_local_id(0)+step < 64 && get_local_id(0)%(2*step) == 0) {
            minBuffer[get_local_id(0)] = min(minBuffer[get_local_id(0)], minBuffer[get_local_id(0)+step]);
            maxBuffer[get_local_id(0)] = max(maxBuffer[get_local_id(0)], maxBuffer[get_local_id(0)+step]);
        }
        barrier(CLK_LOCAL_MEM_FENCE);
    }
    if (get_local_id(0) == 0)
        blockSizeRange[get_group_id(0)] = make_real2(minBuffer[0], maxBuffer[0]);
    if (get_global_id(0) == 0)
        rebuildNeighborList[0] = 0;
}

__kernel void computeSortKeys(__global const real4* restrict blockBoundingBox, __global unsigned int* restrict sortedBlocks, __global real2* restrict blockSizeRange, int numSizes) {
    // Find the total range of sizes recorded by all blocks.

    __local real2 sizeRange;
    if (get_local_id(0) == 0) {
        sizeRange = blockSizeRange[0];
        for (int i = 1; i < numSizes; i++) {
            real2 size = blockSizeRange[i];
            sizeRange.x = min(sizeRange.x, size.x);
            sizeRange.y = max(sizeRange.y, size.y);
        }
        sizeRange.x = LOG(sizeRange.x);
        sizeRange.y = LOG(sizeRange.y);
    }
    barrier(CLK_LOCAL_MEM_FENCE);

    // Sort keys store the bin in the high order part and the block in the low
    // order part.

    int numSizeBins = 20;
    real scale = numSizeBins/(sizeRange.y-sizeRange.x);
    for (unsigned int i = get_global_id(0); i < NUM_BLOCKS; i += get_global_size(0)) {
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
__kernel void sortBoxData(__global const unsigned int* restrict sortedBlocks, __global const real4* restrict blockCenter,
        __global const real4* restrict blockBoundingBox, __global real4* restrict sortedBlockCenter,
        __global real4* restrict sortedBlockBoundingBox, __global const real4* restrict posq, __global const real4* restrict oldPositions,
        __global unsigned int* restrict interactionCount, __global int* restrict rebuildNeighborList, int forceRebuild
#ifdef USE_LARGE_BLOCKS
        , __global real4* restrict largeBlockCenter, __global real4* restrict largeBlockBoundingBox, real4 periodicBoxSize,
        real4 invPeriodicBoxSize, real4 periodicBoxVecX, real4 periodicBoxVecY, real4 periodicBoxVecZ
#endif
        ) {
    for (int i = get_global_id(0); i < NUM_BLOCKS; i += get_global_size(0)) {
        unsigned int index = sortedBlocks[i] & BLOCK_INDEX_MASK;
        sortedBlockCenter[i] = blockCenter[index];
        sortedBlockBoundingBox[i] = blockBoundingBox[index];

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
            minPos = min(minPos, blockPos-width);
            maxPos = max(maxPos, blockPos+width);
        }
        largeBlockCenter[i] = 0.5f*(maxPos+minPos);
        largeBlockBoundingBox[i] = 0.5f*(maxPos-minPos);
#endif
    }

    // Also check whether any atom has moved enough so that we really need to rebuild the neighbor list.

    bool rebuild = forceRebuild;
    for (int i = get_global_id(0); i < NUM_ATOMS; i += get_global_size(0)) {
        real4 delta = oldPositions[i]-posq[i];
        if (delta.x*delta.x + delta.y*delta.y + delta.z*delta.z > 0.25f*PADDING*PADDING)
            rebuild = true;
    }
    if (rebuild) {
        rebuildNeighborList[0] = 1;
        interactionCount[0] = 0;
    }
}

#if SIMD_WIDTH <= 32

#define BUFFER_SIZE 256

__kernel void findBlocksWithInteractions(real4 periodicBoxSize, real4 invPeriodicBoxSize, real4 periodicBoxVecX, real4 periodicBoxVecY, real4 periodicBoxVecZ,
        __global unsigned int* restrict interactionCount, __global int* restrict interactingTiles, __global unsigned int* restrict interactingAtoms,
        __global const real4* restrict posq, unsigned int maxTiles, unsigned int startBlockIndex, unsigned int numBlocks, __global unsigned int* restrict sortedBlocks,
        __global const real4* restrict sortedBlockCenter, __global const real4* restrict sortedBlockBoundingBox,
        __global const unsigned int* restrict exclusionIndices, __global const unsigned int* restrict exclusionRowIndices, __global real4* restrict oldPositions,
        __global const int* restrict rebuildNeighborList
#ifdef USE_LARGE_BLOCKS
        , __global real4* restrict largeBlockCenter, __global real4* restrict largeBlockBoundingBox
#endif
        ) {

    if (rebuildNeighborList[0] == 0)
        return; // The neighbor list doesn't need to be rebuilt.

    const int indexInWarp = get_local_id(0)%32;
    const int warpStart = get_local_id(0)-indexInWarp;
    const int totalWarps = get_global_size(0)/32;
    const int warpIndex = get_global_id(0)/32;
    const int warpMask = (1<<indexInWarp)-1;
    __local int workgroupBuffer[BUFFER_SIZE*(GROUP_SIZE/32)];
    __local int warpExclusions[MAX_EXCLUSIONS*(GROUP_SIZE/32)];
    __local real3 posBuffer[GROUP_SIZE];
    __local volatile unsigned int workgroupTileIndex[GROUP_SIZE/32];
    __local bool includeBlockFlags[GROUP_SIZE];
    __local volatile short2 atomCountBuffer[GROUP_SIZE];
    __local int* buffer = workgroupBuffer+BUFFER_SIZE*(warpStart/32);
    __local int* exclusionsForX = warpExclusions+MAX_EXCLUSIONS*(warpStart/32);
    __local volatile unsigned int* tileStartIndex = workgroupTileIndex+(warpStart/32);
#ifdef USE_LARGE_BLOCKS
    __local bool largeBlockFlags[GROUP_SIZE];
#endif

    // Loop over blocks.

    for (int block1 = startBlockIndex+warpIndex; block1 < startBlockIndex+numBlocks; block1 += totalWarps) {
        // Load data for this block.  Note that all threads in a warp are processing the same block.
        
        int x = sortedBlocks[block1] & BLOCK_INDEX_MASK;
        real4 blockCenterX = sortedBlockCenter[block1];
        real4 blockSizeX = sortedBlockBoundingBox[block1];
        int neighborsInBuffer = 0;
        real3 pos1 = posq[x*TILE_SIZE+indexInWarp].xyz;
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
        posBuffer[get_local_id(0)] = pos1;

        // Load exclusion data for block x.
        
        const int exclusionStart = exclusionRowIndices[x];
        const int exclusionEnd = exclusionRowIndices[x+1];
        const int numExclusions = exclusionEnd-exclusionStart;
        for (int j = indexInWarp; j < numExclusions; j += 32)
            exclusionsForX[j] = exclusionIndices[exclusionStart+j];
        if (MAX_EXCLUSIONS > 32)
            barrier(CLK_LOCAL_MEM_FENCE);
        else
            SYNC_WARPS;
        
        // Loop over atom blocks to search for neighbors.  The threads in a warp compare block1 against 32
        // other blocks in parallel.

#ifdef USE_LARGE_BLOCKS
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
                    real4 largeSize = largeBlockBoundingBox[largeBlockIndex];
                    real4 blockDelta = blockCenterX-largeCenter;
#ifdef USE_PERIODIC
                    APPLY_PERIODIC_TO_DELTA(blockDelta)
#endif
                    blockDelta.x = max((real) 0, fabs(blockDelta.x)-blockSizeX.x-largeSize.x);
                    blockDelta.y = max((real) 0, fabs(blockDelta.y)-blockSizeX.y-largeSize.y);
                    blockDelta.z = max((real) 0, fabs(blockDelta.z)-blockSizeX.z-largeSize.z);
                    includeLargeBlock = (blockDelta.x*blockDelta.x+blockDelta.y*blockDelta.y+blockDelta.z*blockDelta.z < PADDED_CUTOFF_SQUARED);
#ifdef TRICLINIC
                    // The calculation to find the nearest periodic copy is only guaranteed to work if the nearest copy is less than half a box width away.
                    // If there's any possibility we might have missed it, do a detailed check.

                    if (periodicBoxSize.z/2-blockSizeX.z-largeSize.z < PADDED_CUTOFF || periodicBoxSize.y/2-blockSizeX.y-largeSize.y < PADDED_CUTOFF)
                        includeLargeBlock = true;
#endif
                }
                largeBlockFlags[get_local_id(0)] = includeLargeBlock;
                loadedLargeBlocks = 32;
                SYNC_WARPS;
            }
            if (!largeBlockFlags[warpStart+32-(loadedLargeBlocks--)]) {
                // None of the next 32 blocks interact with block 1.

                continue;
            }
#endif
            int block2 = block2Base+indexInWarp;
            bool includeBlock2 = (block2 < NUM_BLOCKS);
            if (includeBlock2) {
                real4 blockCenterY = sortedBlockCenter[block2];
                real4 blockSizeY = sortedBlockBoundingBox[block2];
                real4 blockDelta = blockCenterX-blockCenterY;
#ifdef USE_PERIODIC
                APPLY_PERIODIC_TO_DELTA(blockDelta)
#endif
                includeBlock2 &= (blockDelta.x*blockDelta.x+blockDelta.y*blockDelta.y+blockDelta.z*blockDelta.z < (PADDED_CUTOFF+blockCenterX.w+blockCenterY.w)*(PADDED_CUTOFF+blockCenterX.w+blockCenterY.w));
                blockDelta.x = max((real) 0, fabs(blockDelta.x)-blockSizeX.x-blockSizeY.x);
                blockDelta.y = max((real) 0, fabs(blockDelta.y)-blockSizeX.y-blockSizeY.y);
                blockDelta.z = max((real) 0, fabs(blockDelta.z)-blockSizeX.z-blockSizeY.z);
                includeBlock2 &= (blockDelta.x*blockDelta.x+blockDelta.y*blockDelta.y+blockDelta.z*blockDelta.z < PADDED_CUTOFF_SQUARED);
#ifdef TRICLINIC
                // The calculation to find the nearest periodic copy is only guaranteed to work if the nearest copy is less than half a box width away.
                // If there's any possibility we might have missed it, do a detailed check.

                if (periodicBoxSize.z/2-blockSizeX.z-blockSizeY.z < PADDED_CUTOFF || periodicBoxSize.y/2-blockSizeX.y-blockSizeY.y < PADDED_CUTOFF)
                    includeBlock2 = true;
#endif
                if (includeBlock2) {
                    int y = sortedBlocks[block2] & BLOCK_INDEX_MASK;
                    for (int k = 0; k < numExclusions; k++)
                        includeBlock2 &= (exclusionsForX[k] != y);
                }
            }
            
            // Loop over any blocks we identified as potentially containing neighbors.
            
            includeBlockFlags[get_local_id(0)] = includeBlock2;
            SYNC_WARPS;
            for (int i = 0; i < TILE_SIZE; i++) {
                while (i < TILE_SIZE && !includeBlockFlags[warpStart+i])
                    i++;
                if (i < TILE_SIZE) {
                    int y = sortedBlocks[block2Base+i] & BLOCK_INDEX_MASK;

                    // Check each atom in block Y for interactions.

                    int atom2 = y*TILE_SIZE+indexInWarp;
                    real3 pos2 = posq[atom2].xyz;
#ifdef USE_PERIODIC
                    if (singlePeriodicCopy)
                        APPLY_PERIODIC_TO_POS_WITH_CENTER(pos2, blockCenterX)
#endif
                    bool interacts = false;
                    if (atom2 < NUM_ATOMS) {
#ifdef USE_PERIODIC
                        if (!singlePeriodicCopy) {
                            for (int j = 0; j < TILE_SIZE; j++) {
                                real3 delta = pos2-posBuffer[warpStart+j];
                                APPLY_PERIODIC_TO_DELTA(delta)
                                interacts |= (delta.x*delta.x+delta.y*delta.y+delta.z*delta.z < PADDED_CUTOFF_SQUARED);
                            }
                        }
                        else {
#endif
                            for (int j = 0; j < TILE_SIZE; j++) {
                                real3 delta = pos2-posBuffer[warpStart+j];
                                interacts |= (delta.x*delta.x+delta.y*delta.y+delta.z*delta.z < PADDED_CUTOFF_SQUARED);
                            }
#ifdef USE_PERIODIC
                        }
#endif
                    }
                    
                    // Do a prefix sum to compact the list of atoms.

                    atomCountBuffer[get_local_id(0)].x = (interacts ? 1 : 0);
                    SYNC_WARPS;
                    int whichBuffer = 0;
                    for (int offset = 1; offset < TILE_SIZE; offset *= 2) {
                        if (whichBuffer == 0)
                            atomCountBuffer[get_local_id(0)].y = (indexInWarp < offset ? atomCountBuffer[get_local_id(0)].x : atomCountBuffer[get_local_id(0)].x+atomCountBuffer[get_local_id(0)-offset].x);
                        else
                            atomCountBuffer[get_local_id(0)].x = (indexInWarp < offset ? atomCountBuffer[get_local_id(0)].y : atomCountBuffer[get_local_id(0)].y+atomCountBuffer[get_local_id(0)-offset].y);
                        whichBuffer = 1-whichBuffer;
                        SYNC_WARPS;
                    }
                    
                    // Add any interacting atoms to the buffer.

                    if (interacts)
                        buffer[neighborsInBuffer+atomCountBuffer[get_local_id(0)].y-1] = atom2;
                    neighborsInBuffer += atomCountBuffer[warpStart+TILE_SIZE-1].y;
                    if (neighborsInBuffer > BUFFER_SIZE-TILE_SIZE) {
                        // Store the new tiles to memory.

                        unsigned int tilesToStore = neighborsInBuffer/TILE_SIZE;
                        if (indexInWarp == 0)
                            *tileStartIndex = ATOMIC_ADD(interactionCount, tilesToStore);
                        SYNC_WARPS;
                        unsigned int newTileStartIndex = *tileStartIndex;
                        if (newTileStartIndex+tilesToStore <= maxTiles) {
                            if (indexInWarp < tilesToStore)
                                interactingTiles[newTileStartIndex+indexInWarp] = x;
                            for (int j = 0; j < tilesToStore; j++)
                                interactingAtoms[(newTileStartIndex+j)*TILE_SIZE+indexInWarp] = buffer[indexInWarp+j*TILE_SIZE];
                        }
                        if (indexInWarp+TILE_SIZE*tilesToStore < BUFFER_SIZE)
                            buffer[indexInWarp] = buffer[indexInWarp+TILE_SIZE*tilesToStore];
                        neighborsInBuffer -= TILE_SIZE*tilesToStore;
                   }
                }
                else {
                    SYNC_WARPS;
                }
            }
        }
        
        // If we have a partially filled buffer,  store it to memory.
        
        if (neighborsInBuffer > 0) {
            unsigned int tilesToStore = (neighborsInBuffer+TILE_SIZE-1)/TILE_SIZE;
            if (indexInWarp == 0)
                *tileStartIndex = ATOMIC_ADD(interactionCount, tilesToStore);
            SYNC_WARPS;
            unsigned int newTileStartIndex = *tileStartIndex;
            if (newTileStartIndex+tilesToStore <= maxTiles) {
                if (indexInWarp < tilesToStore)
                    interactingTiles[newTileStartIndex+indexInWarp] = x;
                for (int j = 0; j < tilesToStore; j++)
                    interactingAtoms[(newTileStartIndex+j)*TILE_SIZE+indexInWarp] = (indexInWarp+j*TILE_SIZE < neighborsInBuffer ? buffer[indexInWarp+j*TILE_SIZE] : NUM_ATOMS);
            }
        }
    }
    
    // Record the positions the neighbor list is based on.
    
    for (int i = get_global_id(0); i < NUM_ATOMS; i += get_global_size(0))
        oldPositions[i] = posq[i];
}

#else
// This is the old implementation of finding interacting blocks.  It is quite a bit more complicated,
// and slower on most GPUs.  On AMD, however, it is faster, so we keep it around to use there.

#define BUFFER_SIZE BUFFER_GROUPS*GROUP_SIZE
#define WARP_SIZE 32
#define INVALID -1

/**
 * Perform a parallel prefix sum over an array.  The input values are all assumed to be 0 or 1.
 */
void prefixSum(__local int* sum, __local int2* temp) {
    for (int i = get_local_id(0); i < BUFFER_SIZE; i += get_local_size(0))
        temp[i].x = sum[i];
    barrier(CLK_LOCAL_MEM_FENCE);
    int whichBuffer = 0;
    for (int offset = 1; offset < BUFFER_SIZE; offset *= 2) {
        if (whichBuffer == 0)
            for (int i = get_local_id(0); i < BUFFER_SIZE; i += get_local_size(0))
                temp[i].y = (i < offset ? temp[i].x : temp[i].x+temp[i-offset].x);
        else
            for (int i = get_local_id(0); i < BUFFER_SIZE; i += get_local_size(0))
                temp[i].x = (i < offset ? temp[i].y : temp[i].y+temp[i-offset].y);
        whichBuffer = 1-whichBuffer;
        barrier(CLK_LOCAL_MEM_FENCE);
    }
    if (whichBuffer == 0)
        for (int i = get_local_id(0); i < BUFFER_SIZE; i += get_local_size(0))
            sum[i] = temp[i].x;
    else
        for (int i = get_local_id(0); i < BUFFER_SIZE; i += get_local_size(0))
            sum[i] = temp[i].y;
    barrier(CLK_LOCAL_MEM_FENCE);
}

/**
 * This is called by findBlocksWithInteractions().  It compacts the list of blocks, identifies interactions
 * in them, and writes the result to global memory.
 */
void storeInteractionData(int x, __local int* buffer, __local int* sum, __local int2* temp, __local int* atoms, __local int* numAtoms,
            __local int* baseIndex, __global unsigned int* interactionCount, __global int* interactingTiles, __global unsigned int* interactingAtoms, real4 periodicBoxSize,
            real4 invPeriodicBoxSize, real4 periodicBoxVecX, real4 periodicBoxVecY, real4 periodicBoxVecZ, __global const real4* posq, __local real4* posBuffer,
            real4 blockCenterX, real4 blockSizeX, unsigned int maxTiles, bool finish) {
    const bool singlePeriodicCopy = (0.5f*periodicBoxSize.x-blockSizeX.x >= PADDED_CUTOFF &&
                                     0.5f*periodicBoxSize.y-blockSizeX.y >= PADDED_CUTOFF &&
                                     0.5f*periodicBoxSize.z-blockSizeX.z >= PADDED_CUTOFF);
    if (get_local_id(0) < TILE_SIZE) {
        real4 pos = posq[x*TILE_SIZE+get_local_id(0)];
#ifdef USE_PERIODIC
        if (singlePeriodicCopy) {
            // The box is small enough that we can just translate all the atoms into a single periodic
            // box, then skip having to apply periodic boundary conditions later.
            
            APPLY_PERIODIC_TO_POS_WITH_CENTER(pos, blockCenterX)
        }
#endif
        posBuffer[get_local_id(0)] = pos;
    }
    
    // The buffer is full, so we need to compact it and write out results.  Start by doing a parallel prefix sum.

    barrier(CLK_LOCAL_MEM_FENCE);
    for (int i = get_local_id(0); i < BUFFER_SIZE; i += get_local_size(0))
        sum[i] = (buffer[i] == INVALID ? 0 : 1);
    barrier(CLK_LOCAL_MEM_FENCE);
    prefixSum(sum, temp);
    int numValid = sum[BUFFER_SIZE-1];

    // Compact the buffer.

    for (int i = get_local_id(0); i < BUFFER_SIZE; i += get_local_size(0))
        if (buffer[i] != INVALID)
            temp[sum[i]-1].x = buffer[i];
    barrier(CLK_LOCAL_MEM_FENCE);
    for (int i = get_local_id(0); i < BUFFER_SIZE; i += get_local_size(0))
        buffer[i] = temp[i].x;
    barrier(CLK_LOCAL_MEM_FENCE);

    // Loop over the tiles and find specific interactions in them.

    const int indexInWarp = get_local_id(0)%WARP_SIZE;
    for (int base = 0; base < numValid; base += BUFFER_SIZE/WARP_SIZE) {
        for (int i = get_local_id(0)/WARP_SIZE; i < BUFFER_SIZE/WARP_SIZE && base+i < numValid; i += GROUP_SIZE/WARP_SIZE) {
            // Check each atom in block Y for interactions.
            
            real4 pos = posq[buffer[base+i]*TILE_SIZE+indexInWarp];
#ifdef USE_PERIODIC
            if (singlePeriodicCopy)
                APPLY_PERIODIC_TO_POS_WITH_CENTER(pos, blockCenterX)
#endif
            bool interacts = false;
#ifdef USE_PERIODIC
            if (!singlePeriodicCopy) {
                for (int j = 0; j < TILE_SIZE; j++) {
                    real4 delta = pos-posBuffer[j];
                    APPLY_PERIODIC_TO_DELTA(delta)
                    interacts |= (delta.x*delta.x+delta.y*delta.y+delta.z*delta.z < PADDED_CUTOFF_SQUARED);
                }
            }
            else {
#endif
                for (int j = 0; j < TILE_SIZE; j++) {
                    real4 delta = pos-posBuffer[j];
                    interacts |= (delta.x*delta.x+delta.y*delta.y+delta.z*delta.z < PADDED_CUTOFF_SQUARED);
                }
#ifdef USE_PERIODIC
            }
#endif
            sum[i*WARP_SIZE+indexInWarp] = (interacts ? 1 : 0);
        }
        for (int i = numValid-base+get_local_id(0)/WARP_SIZE; i < BUFFER_SIZE/WARP_SIZE; i += GROUP_SIZE/WARP_SIZE)
            sum[i*WARP_SIZE+indexInWarp] = 0;

        // Compact the list of atoms.

        barrier(CLK_LOCAL_MEM_FENCE);
        prefixSum(sum, temp);
        for (int i = get_local_id(0); i < BUFFER_SIZE; i += get_local_size(0))
            if (sum[i] != (i == 0 ? 0 : sum[i-1]))
                atoms[*numAtoms+sum[i]-1] = buffer[base+i/WARP_SIZE]*TILE_SIZE+indexInWarp;

        // Store them to global memory.

        int atomsToStore = *numAtoms+sum[BUFFER_SIZE-1];
        bool storePartialTile = (finish && base >= numValid-BUFFER_SIZE/WARP_SIZE);
        int tilesToStore = (storePartialTile ? (atomsToStore+TILE_SIZE-1)/TILE_SIZE : atomsToStore/TILE_SIZE);
        if (tilesToStore > 0) {
            if (get_local_id(0) == 0)
                *baseIndex = ATOMIC_ADD(interactionCount, tilesToStore);
            barrier(CLK_LOCAL_MEM_FENCE);
            if (get_local_id(0) == 0)
                *numAtoms = atomsToStore-tilesToStore*TILE_SIZE;
            if (*baseIndex+tilesToStore <= maxTiles) {
                if (get_local_id(0) < tilesToStore)
                    interactingTiles[*baseIndex+get_local_id(0)] = x;
                for (int i = get_local_id(0); i < tilesToStore*TILE_SIZE; i += get_local_size(0))
                    interactingAtoms[*baseIndex*TILE_SIZE+i] = (i < atomsToStore ? atoms[i] : NUM_ATOMS);
            }
        }
        else {
            barrier(CLK_LOCAL_MEM_FENCE);
            if (get_local_id(0) == 0)
                *numAtoms += sum[BUFFER_SIZE-1];
        }
        barrier(CLK_LOCAL_MEM_FENCE);
        if (get_local_id(0) < *numAtoms && !storePartialTile)
            atoms[get_local_id(0)] = atoms[tilesToStore*TILE_SIZE+get_local_id(0)];
    }

    if (numValid == 0 && *numAtoms > 0 && finish) {
        // We didn't have any more tiles to process, but there were some atoms left over from a
        // previous call to this function.  Save them now.

        if (get_local_id(0) == 0)
            *baseIndex = ATOMIC_ADD(interactionCount, 1);
        barrier(CLK_LOCAL_MEM_FENCE);
        if (*baseIndex < maxTiles) {
            if (get_local_id(0) == 0)
                interactingTiles[*baseIndex] = x;
            if (get_local_id(0) < TILE_SIZE)
                interactingAtoms[*baseIndex*TILE_SIZE+get_local_id(0)] = (get_local_id(0) < *numAtoms ? atoms[get_local_id(0)] : NUM_ATOMS);
        }
    }

    // Reset the buffer for processing more tiles.

    for (int i = get_local_id(0); i < BUFFER_SIZE; i += get_local_size(0))
        buffer[i] = INVALID;
    barrier(CLK_LOCAL_MEM_FENCE);
}

/**
 * Compare the bounding boxes for each pair of blocks.  If they are sufficiently far apart,
 * mark them as non-interacting.
 */
__kernel void findBlocksWithInteractions(real4 periodicBoxSize, real4 invPeriodicBoxSize, real4 periodicBoxVecX, real4 periodicBoxVecY, real4 periodicBoxVecZ,
        __global unsigned int* restrict interactionCount, __global int* restrict interactingTiles, __global unsigned int* restrict interactingAtoms,
        __global const real4* restrict posq, unsigned int maxTiles, unsigned int startBlockIndex, unsigned int numBlocks, __global unsigned int* restrict sortedBlocks,
        __global const real4* restrict sortedBlockCenter, __global const real4* restrict sortedBlockBoundingBox,
        __global const unsigned int* restrict exclusionIndices, __global const unsigned int* restrict exclusionRowIndices, __global real4* restrict oldPositions,
        __global const int* restrict rebuildNeighborList
#ifdef USE_LARGE_BLOCKS
        , __global real4* restrict largeBlockCenter, __global real4* restrict largeBlockBoundingBox
#endif
        ) {
    __local int buffer[BUFFER_SIZE];
    __local int sum[BUFFER_SIZE];
    __local int2 temp[BUFFER_SIZE];
    __local int atoms[BUFFER_SIZE+TILE_SIZE];
    __local real4 posBuffer[TILE_SIZE];
    __local int exclusionsForX[MAX_EXCLUSIONS];
    __local int bufferFull;
    __local int globalIndex;
    __local int numAtoms;
#ifdef AMD_ATOMIC_WORK_AROUND
    // Do a byte write to force all memory accesses to interactionCount to use the complete path.
    // This avoids the atomic access from causing all word accesses to other buffers from using the slow complete path.
    // The IF actually causes the write to never be executed, its presence is all that is needed.
    // AMD APP SDK 2.4 has this problem.
    if (get_global_id(0) == get_local_id(0)+1)
        ((__global char*)interactionCount)[sizeof(unsigned int)+1] = 0;
#endif

    if (rebuildNeighborList[0] == 0)
        return; // The neighbor list doesn't need to be rebuilt.

    int valuesInBuffer = 0;
    if (get_local_id(0) == 0)
        bufferFull = false;
    for (int i = 0; i < BUFFER_GROUPS; ++i)
        buffer[i*GROUP_SIZE+get_local_id(0)] = INVALID;
    barrier(CLK_LOCAL_MEM_FENCE);
    
    // Loop over blocks sorted by size.
    
    for (int i = startBlockIndex+get_group_id(0); i < startBlockIndex+numBlocks; i += get_num_groups(0)) {
        if (get_local_id(0) == get_local_size(0)-1)
            numAtoms = 0;
        unsigned int x = sortedBlocks[i] & BLOCK_INDEX_MASK;
        real4 blockCenterX = sortedBlockCenter[i];
        real4 blockSizeX = sortedBlockBoundingBox[i];

        // Load exclusion data for block x.
        
        const int exclusionStart = exclusionRowIndices[x];
        const int exclusionEnd = exclusionRowIndices[x+1];
        const int numExclusions = exclusionEnd-exclusionStart;
        for (int j = get_local_id(0); j < numExclusions; j += get_local_size(0))
            exclusionsForX[j] = exclusionIndices[exclusionStart+j];
        barrier(CLK_LOCAL_MEM_FENCE);
        
        // Compare it to other blocks after this one in sorted order.

        for (int base = i+1; base < NUM_BLOCKS; base += get_local_size(0)) {
            int j = base+get_local_id(0);
            real4 blockCenterY = (j < NUM_BLOCKS ? sortedBlockCenter[j] : (real4) 0);
            real4 blockSizeY = (j < NUM_BLOCKS ? sortedBlockBoundingBox[j] : (real4) 0);
            unsigned int y = (j < NUM_BLOCKS ? sortedBlocks[j] & BLOCK_INDEX_MASK : 0);
            real4 delta = blockCenterX-blockCenterY;
#ifdef USE_PERIODIC
            APPLY_PERIODIC_TO_DELTA(delta)
#endif
            delta.x = max((real) 0, fabs(delta.x)-blockSizeX.x-blockSizeY.x);
            delta.y = max((real) 0, fabs(delta.y)-blockSizeX.y-blockSizeY.y);
            delta.z = max((real) 0, fabs(delta.z)-blockSizeX.z-blockSizeY.z);
            bool hasExclusions = false;
            for (int k = 0; k < numExclusions; k++)
                hasExclusions |= (exclusionsForX[k] == y);
            if (j < NUM_BLOCKS && delta.x*delta.x+delta.y*delta.y+delta.z*delta.z < PADDED_CUTOFF_SQUARED && !hasExclusions) {
                // Add this tile to the buffer.

                int bufferIndex = valuesInBuffer*GROUP_SIZE+get_local_id(0);
                buffer[bufferIndex] = y;
                valuesInBuffer++;
                if (!bufferFull && valuesInBuffer == BUFFER_GROUPS)
                    bufferFull = true;
            }
            barrier(CLK_LOCAL_MEM_FENCE);
            if (bufferFull) {
                storeInteractionData(x, buffer, sum, temp, atoms, &numAtoms, &globalIndex, interactionCount, interactingTiles, interactingAtoms, periodicBoxSize, invPeriodicBoxSize, periodicBoxVecX, periodicBoxVecY, periodicBoxVecZ, posq, posBuffer, blockCenterX, blockSizeX, maxTiles, false);
                valuesInBuffer = 0;
                if (get_local_id(0) == 0)
                    bufferFull = false;
            }
            barrier(CLK_LOCAL_MEM_FENCE);
        }
        storeInteractionData(x, buffer, sum, temp, atoms, &numAtoms, &globalIndex, interactionCount, interactingTiles, interactingAtoms, periodicBoxSize, invPeriodicBoxSize, periodicBoxVecX, periodicBoxVecY, periodicBoxVecZ, posq, posBuffer, blockCenterX, blockSizeX, maxTiles, true);
    }
    
    // Record the positions the neighbor list is based on.
    
    for (int i = get_global_id(0); i < NUM_ATOMS; i += get_global_size(0))
        oldPositions[i] = posq[i];
}

#endif

