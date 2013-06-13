#pragma OPENCL EXTENSION cl_khr_global_int32_base_atomics : enable
#pragma OPENCL EXTENSION cl_khr_byte_addressable_store : enable
#define BUFFER_SIZE BUFFER_GROUPS*GROUP_SIZE
#define WARP_SIZE 32
#define INVALID 0xFFFF

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

/**
 * Perform a parallel prefix sum over an array.  The input values are all assumed to be 0 or 1.
 */
void prefixSum(__local short* sum, __local ushort2* temp) {
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
void storeInteractionData(int x, __local unsigned short* buffer, __local short* sum, __local ushort2* temp, __local int* atoms, __local int* numAtoms,
            __local int* baseIndex, __global unsigned int* interactionCount, __global int* interactingTiles, __global unsigned int* interactingAtoms, real4 periodicBoxSize,
            real4 invPeriodicBoxSize, __global const real4* posq, __local real4* posBuffer, real4 blockCenterX, real4 blockSizeX, unsigned int maxTiles, bool finish) {
    const bool singlePeriodicCopy = (0.5f*periodicBoxSize.x-blockSizeX.x >= PADDED_CUTOFF &&
                                     0.5f*periodicBoxSize.y-blockSizeX.y >= PADDED_CUTOFF &&
                                     0.5f*periodicBoxSize.z-blockSizeX.z >= PADDED_CUTOFF);
    if (get_local_id(0) < TILE_SIZE) {
        real4 pos = posq[x*TILE_SIZE+get_local_id(0)];
#ifdef USE_PERIODIC
        if (singlePeriodicCopy) {
            // The box is small enough that we can just translate all the atoms into a single periodic
            // box, then skip having to apply periodic boundary conditions later.
            
            pos.xyz -= floor((pos.xyz-blockCenterX.xyz)*invPeriodicBoxSize.xyz+0.5f)*periodicBoxSize.xyz;
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
                pos.xyz -= floor((pos.xyz-blockCenterX.xyz)*invPeriodicBoxSize.xyz+0.5f)*periodicBoxSize.xyz;
#endif
            bool interacts = false;
#ifdef USE_PERIODIC
            if (!singlePeriodicCopy) {
                for (int j = 0; j < TILE_SIZE; j++) {
                    real4 delta = pos-posBuffer[j];
                    delta.xyz -= floor(delta.xyz*invPeriodicBoxSize.xyz+0.5f)*periodicBoxSize.xyz;
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
                *baseIndex = atom_add(interactionCount, tilesToStore);
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
            *baseIndex = atom_add(interactionCount, 1);
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
__kernel void findBlocksWithInteractions(real4 periodicBoxSize, real4 invPeriodicBoxSize, __global unsigned int* restrict interactionCount,
        __global int* restrict interactingTiles, __global unsigned int* restrict interactingAtoms, __global const real4* restrict posq, unsigned int maxTiles, unsigned int startBlockIndex,
        unsigned int numBlocks, __global real2* restrict sortedBlocks, __global const real4* restrict sortedBlockCenter, __global const real4* restrict sortedBlockBoundingBox,
        __global const unsigned int* restrict exclusionIndices, __global const unsigned int* restrict exclusionRowIndices, __global real4* restrict oldPositions,
        __global const int* restrict rebuildNeighborList) {
    __local unsigned short buffer[BUFFER_SIZE];
    __local short sum[BUFFER_SIZE];
    __local ushort2 temp[BUFFER_SIZE];
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
        real2 sortedKey = sortedBlocks[i];
        int x = (int) sortedKey.y;
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
            real2 sortedKey2 = (j < NUM_BLOCKS ? sortedBlocks[j] : (real2) 0);
            real4 blockCenterY = (j < NUM_BLOCKS ? sortedBlockCenter[j] : (real4) 0);
            real4 blockSizeY = (j < NUM_BLOCKS ? sortedBlockBoundingBox[j] : (real4) 0);
            unsigned short y = (unsigned short) sortedKey2.y;
            real4 delta = blockCenterX-blockCenterY;
#ifdef USE_PERIODIC
            delta.x -= floor(delta.x*invPeriodicBoxSize.x+0.5f)*periodicBoxSize.x;
            delta.y -= floor(delta.y*invPeriodicBoxSize.y+0.5f)*periodicBoxSize.y;
            delta.z -= floor(delta.z*invPeriodicBoxSize.z+0.5f)*periodicBoxSize.z;
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
                storeInteractionData(x, buffer, sum, temp, atoms, &numAtoms, &globalIndex, interactionCount, interactingTiles, interactingAtoms, periodicBoxSize, invPeriodicBoxSize, posq, posBuffer, blockCenterX, blockSizeX, maxTiles, false);
                valuesInBuffer = 0;
                if (get_local_id(0) == 0)
                    bufferFull = false;
            }
            barrier(CLK_LOCAL_MEM_FENCE);
        }
        storeInteractionData(x, buffer, sum, temp, atoms, &numAtoms, &globalIndex, interactionCount, interactingTiles, interactingAtoms, periodicBoxSize, invPeriodicBoxSize, posq, posBuffer, blockCenterX, blockSizeX, maxTiles, true);
    }
    
    // Record the positions the neighbor list is based on.
    
    for (int i = get_global_id(0); i < NUM_ATOMS; i += get_global_size(0))
        oldPositions[i] = posq[i];
}
