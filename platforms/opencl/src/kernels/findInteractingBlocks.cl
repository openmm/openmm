#pragma OPENCL EXTENSION cl_khr_global_int32_base_atomics : enable
#pragma OPENCL EXTENSION cl_khr_byte_addressable_store : enable
#define BUFFER_SIZE 256

/**
 * Find a bounding box for the atoms in each block.
 */
__kernel void findBlockBounds(int numAtoms, real4 periodicBoxSize, real4 invPeriodicBoxSize, __global const real4* restrict posq,
        __global real4* restrict blockCenter, __global real4* restrict blockBoundingBox, __global int* restrict rebuildNeighborList,
        __global real2* restrict sortedBlocks) {
    int index = get_global_id(0);
    int base = index*TILE_SIZE;
    while (base < numAtoms) {
        real4 pos = posq[base];
#ifdef USE_PERIODIC
        pos.xyz -= floor(pos.xyz*invPeriodicBoxSize.xyz)*periodicBoxSize.xyz;
#endif
        real4 minPos = pos;
        real4 maxPos = pos;
        int last = min(base+TILE_SIZE, numAtoms);
        for (int i = base+1; i < last; i++) {
            pos = posq[i];
#ifdef USE_PERIODIC
            real4 center = 0.5f*(maxPos+minPos);
            pos.xyz -= floor((pos.xyz-center.xyz)*invPeriodicBoxSize.xyz+0.5f)*periodicBoxSize.xyz;
#endif
            minPos = min(minPos, pos);
            maxPos = max(maxPos, pos);
        }
        real4 blockSize = 0.5f*(maxPos-minPos);
        blockBoundingBox[index] = blockSize;
        blockCenter[index] = 0.5f*(maxPos+minPos);
        sortedBlocks[index] = (real2) (blockSize.x+blockSize.y+blockSize.z, index);
        index += get_global_size(0);
        base = index*TILE_SIZE;
    }
    if (get_global_id(0) == 0)
        rebuildNeighborList[0] = 0;
}

/**
 * Sort the data about bounding boxes so it can be accessed more efficiently in the next kernel.
 */
__kernel void sortBoxData(__global const real2* restrict sortedBlock, __global const real4* restrict blockCenter,
        __global const real4* restrict blockBoundingBox, __global real4* restrict sortedBlockCenter,
        __global real4* restrict sortedBlockBoundingBox, __global const real4* restrict posq, __global const real4* restrict oldPositions,
        __global unsigned int* restrict interactionCount, __global int* restrict rebuildNeighborList) {
    for (int i = get_global_id(0); i < NUM_BLOCKS; i += get_global_size(0)) {
        int index = (int) sortedBlock[i].y;
        sortedBlockCenter[i] = blockCenter[index];
        sortedBlockBoundingBox[i] = blockBoundingBox[index];
    }
    
    // Also check whether any atom has moved enough so that we really need to rebuild the neighbor list.

    bool rebuild = false;
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

__kernel void findBlocksWithInteractions(real4 periodicBoxSize, real4 invPeriodicBoxSize, __global unsigned int* restrict interactionCount,
        __global int* restrict interactingTiles, __global unsigned int* restrict interactingAtoms, __global const real4* restrict posq, unsigned int maxTiles, unsigned int startBlockIndex,
        unsigned int numBlocks, __global real2* restrict sortedBlocks, __global const real4* restrict sortedBlockCenter, __global const real4* restrict sortedBlockBoundingBox,
        __global const unsigned int* restrict exclusionIndices, __global const unsigned int* restrict exclusionRowIndices, __global real4* restrict oldPositions,
        __global const int* restrict rebuildNeighborList) {

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
    __local volatile int workgroupTileIndex[GROUP_SIZE/32];
    __local bool includeBlockFlags[GROUP_SIZE];
    __local short2 atomCountBuffer[GROUP_SIZE];
    __local int* buffer = workgroupBuffer+BUFFER_SIZE*(warpStart/32);
    __local int* exclusionsForX = warpExclusions+MAX_EXCLUSIONS*(warpStart/32);
    __local volatile int* tileStartIndex = workgroupTileIndex+(warpStart/32);

    // Loop over blocks.
    
    for (int block1 = startBlockIndex+warpIndex; block1 < startBlockIndex+numBlocks; block1 += totalWarps) {
        // Load data for this block.  Note that all threads in a warp are processing the same block.
        
        real2 sortedKey = sortedBlocks[block1];
        int x = (int) sortedKey.y;
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
            
            pos1.xyz -= floor((pos1.xyz-blockCenterX.xyz)*invPeriodicBoxSize.xyz+0.5f)*periodicBoxSize.xyz;
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

        for (int block2Base = block1+1; block2Base < NUM_BLOCKS; block2Base += 32) {
            int block2 = block2Base+indexInWarp;
            bool includeBlock2 = (block2 < NUM_BLOCKS);
            if (includeBlock2) {
                real4 blockCenterY = (block2 < NUM_BLOCKS ? sortedBlockCenter[block2] : (real4) (0));
                real4 blockSizeY = (block2 < NUM_BLOCKS ? sortedBlockBoundingBox[block2] : (real4) (0));
                real4 blockDelta = blockCenterX-blockCenterY;
#ifdef USE_PERIODIC
                blockDelta.xyz -= floor(blockDelta.xyz*invPeriodicBoxSize.xyz+0.5f)*periodicBoxSize.xyz;
#endif
                blockDelta.x = max((real) 0, fabs(blockDelta.x)-blockSizeX.x-blockSizeY.x);
                blockDelta.y = max((real) 0, fabs(blockDelta.y)-blockSizeX.y-blockSizeY.y);
                blockDelta.z = max((real) 0, fabs(blockDelta.z)-blockSizeX.z-blockSizeY.z);
                includeBlock2 &= (blockDelta.x*blockDelta.x+blockDelta.y*blockDelta.y+blockDelta.z*blockDelta.z < PADDED_CUTOFF_SQUARED);
                if (includeBlock2) {
                    unsigned short y = (unsigned short) sortedBlocks[block2].y;
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
                    unsigned short y = (unsigned short) sortedBlocks[block2Base+i].y;

                    // Check each atom in block Y for interactions.

                    int start = y*TILE_SIZE;
                    int atom2 = start+indexInWarp;
                    real3 pos2 = posq[atom2].xyz;
#ifdef USE_PERIODIC
                    if (singlePeriodicCopy)
                        pos2.xyz -= floor((pos2.xyz-blockCenterX.xyz)*invPeriodicBoxSize.xyz+0.5f)*periodicBoxSize.xyz;
#endif
                    bool interacts = false;
                    if (atom2 < NUM_ATOMS) {
#ifdef USE_PERIODIC
                        if (!singlePeriodicCopy) {
                            for (int j = 0; j < TILE_SIZE; j++) {
                                real3 delta = pos2-posBuffer[warpStart+j];
                                delta.xyz -= floor(delta.xyz*invPeriodicBoxSize.xyz+0.5f)*periodicBoxSize.xyz;
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

                        int tilesToStore = neighborsInBuffer/TILE_SIZE;
                        if (indexInWarp == 0)
                            *tileStartIndex = atom_add(interactionCount, tilesToStore);
                        SYNC_WARPS;
                        int newTileStartIndex = *tileStartIndex;
                        if (newTileStartIndex+tilesToStore <= maxTiles) {
                            if (indexInWarp < tilesToStore)
                                interactingTiles[newTileStartIndex+indexInWarp] = x;
                            for (int j = 0; j < tilesToStore; j++)
                                interactingAtoms[(newTileStartIndex+j)*TILE_SIZE+indexInWarp] = buffer[indexInWarp+j*TILE_SIZE];
                        }
                        buffer[indexInWarp] = buffer[indexInWarp+TILE_SIZE*tilesToStore];
                        neighborsInBuffer -= TILE_SIZE*tilesToStore;
                   }
                }
            }
        }
        
        // If we have a partially filled buffer,  store it to memory.
        
        if (neighborsInBuffer > 0) {
            int tilesToStore = (neighborsInBuffer+TILE_SIZE-1)/TILE_SIZE;
            if (indexInWarp == 0)
                *tileStartIndex = atom_add(interactionCount, tilesToStore);
            SYNC_WARPS;
            int newTileStartIndex = *tileStartIndex;
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