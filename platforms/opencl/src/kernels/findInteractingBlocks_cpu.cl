#pragma OPENCL EXTENSION cl_khr_global_int32_base_atomics : enable
#pragma OPENCL EXTENSION cl_khr_byte_addressable_store : enable
#define BUFFER_SIZE BUFFER_GROUPS*GROUP_SIZE

/**
 * Find a bounding box for the atoms in each block.
 */
__kernel void findBlockBounds(int numAtoms, real4 periodicBoxSize, real4 invPeriodicBoxSize, real4 periodicBoxVecX, real4 periodicBoxVecY, real4 periodicBoxVecZ,
        __global const real4* restrict posq, __global real4* restrict blockCenter, __global real4* restrict blockBoundingBox, __global int* restrict rebuildNeighborList,
        __global real2* restrict sortedBlocks) {
    int index = get_global_id(0);
    int base = index*TILE_SIZE;
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
        __global unsigned int* restrict interactionCount, __global int* restrict rebuildNeighborList, int forceRebuild) {
    for (int i = get_global_id(0); i < NUM_BLOCKS; i += get_global_size(0)) {
        int index = (int) sortedBlock[i].y;
        sortedBlockCenter[i] = blockCenter[index];
        sortedBlockBoundingBox[i] = blockBoundingBox[index];
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

/**
 * This is called by findBlocksWithInteractions().  It compacts the list of blocks and writes them
 * to global memory.
 */
void storeInteractionData(unsigned short x, unsigned short* buffer, int* atoms, int* numAtoms, int numValid, __global unsigned int* interactionCount,
            __global int* interactingTiles, __global unsigned int* interactingAtoms, real4 periodicBoxSize, real4 invPeriodicBoxSize, real4 periodicBoxVecX,
            real4 periodicBoxVecY, real4 periodicBoxVecZ, __global const real4* posq, real4 blockCenterX, real4 blockSizeX, unsigned int maxTiles, bool finish) {
    real4 posBuffer[TILE_SIZE];
    const bool singlePeriodicCopy = (0.5f*periodicBoxSize.x-blockSizeX.x >= PADDED_CUTOFF &&
                                     0.5f*periodicBoxSize.y-blockSizeX.y >= PADDED_CUTOFF &&
                                     0.5f*periodicBoxSize.z-blockSizeX.z >= PADDED_CUTOFF);
    for (int i = 0; i < TILE_SIZE; i++) {
        real4 pos = posq[x*TILE_SIZE+i];
#ifdef USE_PERIODIC
        if (singlePeriodicCopy) {
            // The box is small enough that we can just translate all the atoms into a single periodic
            // box, then skip having to apply periodic boundary conditions later.
            
            APPLY_PERIODIC_TO_POS_WITH_CENTER(pos, blockCenterX)
        }
#endif
        posBuffer[i] = pos;
    }

    // Loop over the tiles and find specific interactions in them.

    for (int tile = 0; tile < numValid; tile++) {
        for (int indexInTile = 0; indexInTile < TILE_SIZE; indexInTile++) {
            // Check each atom in block Y for interactions.
            
            int atom = buffer[tile]*TILE_SIZE+indexInTile;
            real4 pos = posq[atom];
#ifdef USE_PERIODIC
            if (singlePeriodicCopy)
		APPLY_PERIODIC_TO_POS_WITH_CENTER(pos, blockCenterX)
#endif
            bool interacts = false;
#ifdef USE_PERIODIC
            if (!singlePeriodicCopy) {
                for (int j = 0; j < TILE_SIZE && !interacts; j++) {
                    real4 delta = pos-posBuffer[j];
                    APPLY_PERIODIC_TO_DELTA(delta)
                    interacts = (delta.x*delta.x+delta.y*delta.y+delta.z*delta.z < PADDED_CUTOFF_SQUARED);
                }
            }
            else {
#endif
                for (int j = 0; j < TILE_SIZE && !interacts; j++) {
                    real4 delta = pos-posBuffer[j];
                    interacts = (delta.x*delta.x+delta.y*delta.y+delta.z*delta.z < PADDED_CUTOFF_SQUARED);
                }
#ifdef USE_PERIODIC
            }
#endif
            if (interacts)
                atoms[(*numAtoms)++] = atom;
            if (*numAtoms == BUFFER_SIZE) {
                // The atoms buffer is full, so store it to global memory.
                
                int tilesToStore = BUFFER_SIZE/TILE_SIZE;
                int baseIndex = atom_add(interactionCount, tilesToStore);
                if (baseIndex+tilesToStore <= maxTiles) {
                    for (int i = 0; i < tilesToStore; i++) {
                        interactingTiles[baseIndex+i] = x;
                        for (int j = 0; j < TILE_SIZE; j++)
                            interactingAtoms[(baseIndex+i)*TILE_SIZE+j] = atoms[i*TILE_SIZE+j];
                    }
                }
                *numAtoms = 0;
            }
        }
    }
    
    if (*numAtoms > 0 && finish) {
        // There are some leftover atoms, so save them now.
        
        int tilesToStore = (*numAtoms+TILE_SIZE-1)/TILE_SIZE;
        int baseIndex = atom_add(interactionCount, tilesToStore);
        if (baseIndex+tilesToStore <= maxTiles) {
            for (int i = 0; i < tilesToStore; i++) {
                interactingTiles[baseIndex+i] = x;
                for (int j = 0; j < TILE_SIZE; j++) {
                    int index = i*TILE_SIZE+j;
                    interactingAtoms[(baseIndex+i)*TILE_SIZE+j] = (index < *numAtoms ? atoms[index] : NUM_ATOMS);
                }
            }
        }
    }
}

/**
 * Compare the bounding boxes for each pair of blocks.  If they are sufficiently far apart,
 * mark them as non-interacting.
 */
__kernel void findBlocksWithInteractions(real4 periodicBoxSize, real4 invPeriodicBoxSize, real4 periodicBoxVecX, real4 periodicBoxVecY, real4 periodicBoxVecZ,
        __global unsigned int* restrict interactionCount, __global int* restrict interactingTiles, __global unsigned int* restrict interactingAtoms,
        __global const real4* restrict posq, unsigned int maxTiles, unsigned int startBlockIndex, unsigned int numBlocks, __global real2* restrict sortedBlocks,
        __global const real4* restrict sortedBlockCenter, __global const real4* restrict sortedBlockBoundingBox,
        __global const unsigned int* restrict exclusionIndices, __global const unsigned int* restrict exclusionRowIndices, __global real4* restrict oldPositions,
        __global const int* restrict rebuildNeighborList) {
    if (rebuildNeighborList[0] == 0)
        return; // The neighbor list doesn't need to be rebuilt.
    unsigned short buffer[BUFFER_SIZE];
    int atoms[BUFFER_SIZE];
    int exclusionsForX[MAX_EXCLUSIONS];
    int valuesInBuffer;
    int numAtoms;
    
    // Loop over blocks sorted by size.
    
    for (int i = startBlockIndex+get_group_id(0); i < startBlockIndex+numBlocks; i += get_num_groups(0)) {
        valuesInBuffer = 0;
        numAtoms = 0;
        real2 sortedKey = sortedBlocks[i];
        unsigned short x = (unsigned short) sortedKey.y;
        real4 blockCenterX = sortedBlockCenter[i];
        real4 blockSizeX = sortedBlockBoundingBox[i];

        // Load exclusion data for block x.
        
        const int exclusionStart = exclusionRowIndices[x];
        const int exclusionEnd = exclusionRowIndices[x+1];
        const int numExclusions = exclusionEnd-exclusionStart;
        for (int j = 0; j < numExclusions; j++)
            exclusionsForX[j] = exclusionIndices[exclusionStart+j];
        
        // Compare it to other blocks after this one in sorted order.
        
        for (int j = i+1; j < NUM_BLOCKS; j++) {
            real2 sortedKey2 = sortedBlocks[j];
            unsigned short y = (unsigned short) sortedKey2.y;
            bool hasExclusions = false;
            for (int k = 0; k < numExclusions; k++)
                hasExclusions |= (exclusionsForX[k] == y);
            if (hasExclusions)
                continue;
            real4 blockCenterY = sortedBlockCenter[j];
            real4 blockSizeY = sortedBlockBoundingBox[j];
            real4 delta = blockCenterX-blockCenterY;
#ifdef USE_PERIODIC
            APPLY_PERIODIC_TO_DELTA(delta)
#endif
            delta.x = max((real) 0, fabs(delta.x)-blockSizeX.x-blockSizeY.x);
            delta.y = max((real) 0, fabs(delta.y)-blockSizeX.y-blockSizeY.y);
            delta.z = max((real) 0, fabs(delta.z)-blockSizeX.z-blockSizeY.z);
            if (delta.x*delta.x+delta.y*delta.y+delta.z*delta.z < PADDED_CUTOFF_SQUARED) {
                // Add this tile to the buffer.

                buffer[valuesInBuffer++] = y;
                if (valuesInBuffer == BUFFER_SIZE) {
                    storeInteractionData(x, buffer, atoms, &numAtoms, valuesInBuffer, interactionCount, interactingTiles, interactingAtoms, periodicBoxSize, invPeriodicBoxSize, periodicBoxVecX, periodicBoxVecY, periodicBoxVecZ, posq, blockCenterX, blockSizeX, maxTiles, false);
                    valuesInBuffer = 0;
                }
            }
        }
        storeInteractionData(x, buffer, atoms, &numAtoms, valuesInBuffer, interactionCount, interactingTiles, interactingAtoms, periodicBoxSize, invPeriodicBoxSize, periodicBoxVecX, periodicBoxVecY, periodicBoxVecZ, posq, blockCenterX, blockSizeX, maxTiles, true);
    }
    
    // Record the positions the neighbor list is based on.
    
    for (int i = get_global_id(0); i < NUM_ATOMS; i += get_global_size(0))
        oldPositions[i] = posq[i];
}
