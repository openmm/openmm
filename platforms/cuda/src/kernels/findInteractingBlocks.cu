#define GROUP_SIZE 256
#define BUFFER_GROUPS 2
#define BUFFER_SIZE BUFFER_GROUPS*GROUP_SIZE
#define WARP_SIZE 32
#define INVALID 0xFFFF

/**
 * Find a bounding box for the atoms in each block.
 */
extern "C" __global__ void findBlockBounds(int numAtoms, real4 periodicBoxSize, real4 invPeriodicBoxSize, const real4* __restrict__ posq,
        real4* __restrict__ blockCenter, real4* __restrict__ blockBoundingBox, int* __restrict__ rebuildNeighborList, real2* __restrict__ sortedBlocks) {
    int index = blockIdx.x*blockDim.x+threadIdx.x;
    int base = index*TILE_SIZE;
    while (base < numAtoms) {
        real4 pos = posq[base];
#ifdef USE_PERIODIC
        pos.x -= floor(pos.x*invPeriodicBoxSize.x)*periodicBoxSize.x;
        pos.y -= floor(pos.y*invPeriodicBoxSize.y)*periodicBoxSize.y;
        pos.z -= floor(pos.z*invPeriodicBoxSize.z)*periodicBoxSize.z;
#endif
        real4 minPos = pos;
        real4 maxPos = pos;
        int last = min(base+TILE_SIZE, numAtoms);
        for (int i = base+1; i < last; i++) {
            pos = posq[i];
#ifdef USE_PERIODIC
            real4 center = 0.5f*(maxPos+minPos);
            pos.x -= floor((pos.x-center.x)*invPeriodicBoxSize.x+0.5f)*periodicBoxSize.x;
            pos.y -= floor((pos.y-center.y)*invPeriodicBoxSize.y+0.5f)*periodicBoxSize.y;
            pos.z -= floor((pos.z-center.z)*invPeriodicBoxSize.z+0.5f)*periodicBoxSize.z;
#endif
            minPos = make_real4(min(minPos.x,pos.x), min(minPos.y,pos.y), min(minPos.z,pos.z), 0);
            maxPos = make_real4(max(maxPos.x,pos.x), max(maxPos.y,pos.y), max(maxPos.z,pos.z), 0);
        }
        real4 blockSize = 0.5f*(maxPos-minPos);
        blockBoundingBox[index] = blockSize;
        blockCenter[index] = 0.5f*(maxPos+minPos);
        sortedBlocks[index] = make_real2(blockSize.x+blockSize.y+blockSize.z, index);
        index += blockDim.x*gridDim.x;
        base = index*TILE_SIZE;
    }
    if (blockIdx.x == 0 && threadIdx.x == 0)
        rebuildNeighborList[0] = 0;
}

/**
 * Sort the data about bounding boxes so it can be accessed more efficiently in the next kernel.
 */
extern "C" __global__ void sortBoxData(const real2* __restrict__ sortedBlock, const real4* __restrict__ blockCenter,
        const real4* __restrict__ blockBoundingBox, real4* __restrict__ sortedBlockCenter,
        real4* __restrict__ sortedBlockBoundingBox, const real4* __restrict__ posq, const real4* __restrict__ oldPositions,
        unsigned int* __restrict__ interactionCount, int* __restrict__ rebuildNeighborList) {
    for (int i = threadIdx.x+blockIdx.x*blockDim.x; i < NUM_BLOCKS; i += blockDim.x*gridDim.x) {
        int index = (int) sortedBlock[i].y;
        sortedBlockCenter[i] = blockCenter[index];
        sortedBlockBoundingBox[i] = blockBoundingBox[index];
    }
    
    // Also check whether any atom has moved enough so that we really need to rebuild the neighbor list.

    bool rebuild = false;
    for (int i = threadIdx.x+blockIdx.x*blockDim.x; i < NUM_ATOMS; i += blockDim.x*gridDim.x) {
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
__device__ void prefixSum(short* sum, ushort2* temp) {
#if __CUDA_ARCH__ >= 300
    const int indexInWarp = threadIdx.x%WARP_SIZE;
    const int warpMask = (2<<indexInWarp)-1;
    for (int base = 0; base < BUFFER_SIZE; base += blockDim.x)
        temp[base+threadIdx.x].x = __popc(__ballot(sum[base+threadIdx.x])&warpMask);
    __syncthreads();
    if (threadIdx.x < BUFFER_SIZE/WARP_SIZE) {
        int multiWarpSum = temp[(threadIdx.x+1)*WARP_SIZE-1].x;
        for (int offset = 1; offset < BUFFER_SIZE/WARP_SIZE; offset *= 2) {
            short n = __shfl_up(multiWarpSum, offset, WARP_SIZE);
            if (indexInWarp >= offset)
                multiWarpSum += n;
        }
        temp[threadIdx.x].y = multiWarpSum;
    }
    __syncthreads();
    for (int i = threadIdx.x; i < BUFFER_SIZE; i += blockDim.x)
        sum[i] = temp[i].x+(i < WARP_SIZE ? 0 : temp[i/WARP_SIZE-1].y);
    __syncthreads();
#else
    for (int i = threadIdx.x; i < BUFFER_SIZE; i += blockDim.x)
        temp[i].x = sum[i];
    __syncthreads();
    int whichBuffer = 0;
    for (int offset = 1; offset < BUFFER_SIZE; offset *= 2) {
        if (whichBuffer == 0)
            for (int i = threadIdx.x; i < BUFFER_SIZE; i += blockDim.x)
                temp[i].y = (i < offset ? temp[i].x : temp[i].x+temp[i-offset].x);
        else
            for (int i = threadIdx.x; i < BUFFER_SIZE; i += blockDim.x)
                temp[i].x = (i < offset ? temp[i].y : temp[i].y+temp[i-offset].y);
        whichBuffer = 1-whichBuffer;
        __syncthreads();
    }
    if (whichBuffer == 0)
        for (int i = threadIdx.x; i < BUFFER_SIZE; i += blockDim.x)
            sum[i] = temp[i].x;
    else
        for (int i = threadIdx.x; i < BUFFER_SIZE; i += blockDim.x)
            sum[i] = temp[i].y;
    __syncthreads();
#endif
}

/**
 * This is called by findBlocksWithInteractions(). It compacts the list of blocks, identifies interactions
 * in them, and writes the result to global memory.
 */
__device__ void storeInteractionData(int x, unsigned short* buffer, short* sum, ushort2* temp, int* atoms, int& numAtoms,
            int& baseIndex, unsigned int* interactionCount, int* interactingTiles, unsigned int* interactingAtoms, real4 periodicBoxSize,
            real4 invPeriodicBoxSize, const real4* posq, real3* posBuffer, real4 blockCenterX, real4 blockSizeX, unsigned int maxTiles, bool finish) {
    const bool singlePeriodicCopy = (0.5f*periodicBoxSize.x-blockSizeX.x >= PADDED_CUTOFF &&
                                     0.5f*periodicBoxSize.y-blockSizeX.y >= PADDED_CUTOFF &&
                                     0.5f*periodicBoxSize.z-blockSizeX.z >= PADDED_CUTOFF);
    if (threadIdx.x < TILE_SIZE) {
        real3 pos = trimTo3(posq[x*TILE_SIZE+threadIdx.x]);
        posBuffer[threadIdx.x] = pos;
#ifdef USE_PERIODIC
        if (singlePeriodicCopy) {
            // The box is small enough that we can just translate all the atoms into a single periodic
            // box, then skip having to apply periodic boundary conditions later.
            
            pos.x -= floor((pos.x-blockCenterX.x)*invPeriodicBoxSize.x+0.5f)*periodicBoxSize.x;
            pos.y -= floor((pos.y-blockCenterX.y)*invPeriodicBoxSize.y+0.5f)*periodicBoxSize.y;
            pos.z -= floor((pos.z-blockCenterX.z)*invPeriodicBoxSize.z+0.5f)*periodicBoxSize.z;
            posBuffer[threadIdx.x] = pos;
        }
#endif
    }
    
    // The buffer is full, so we need to compact it and write out results.  Start by doing a parallel prefix sum.

    for (int i = threadIdx.x; i < BUFFER_SIZE; i += blockDim.x)
        sum[i] = (buffer[i] == INVALID ? 0 : 1);
    __syncthreads();
    prefixSum(sum, temp);
    int numValid = sum[BUFFER_SIZE-1];

    // Compact the buffer.

    for (int i = threadIdx.x; i < BUFFER_SIZE; i += blockDim.x)
        if (buffer[i] != INVALID)
            temp[sum[i]-1].x = buffer[i];
    __syncthreads();
    for (int i = threadIdx.x; i < BUFFER_SIZE; i += blockDim.x)
        buffer[i] = temp[i].x;
    __syncthreads();

    // Loop over the tiles and find specific interactions in them.

    const int indexInWarp = threadIdx.x%WARP_SIZE;
    for (int base = 0; base < numValid; base += BUFFER_SIZE/WARP_SIZE) {
        for (int i = threadIdx.x/WARP_SIZE; i < BUFFER_SIZE/WARP_SIZE && base+i < numValid; i += GROUP_SIZE/WARP_SIZE) {
            // Check each atom in block Y for interactions.
            
            real3 pos = trimTo3(posq[buffer[base+i]*TILE_SIZE+indexInWarp]);
#ifdef USE_PERIODIC
            if (singlePeriodicCopy) {
                pos.x -= floor((pos.x-blockCenterX.x)*invPeriodicBoxSize.x+0.5f)*periodicBoxSize.x;
                pos.y -= floor((pos.y-blockCenterX.y)*invPeriodicBoxSize.y+0.5f)*periodicBoxSize.y;
                pos.z -= floor((pos.z-blockCenterX.z)*invPeriodicBoxSize.z+0.5f)*periodicBoxSize.z;
            }
#endif
            bool interacts = false;
#ifdef USE_PERIODIC
            if (!singlePeriodicCopy) {
                for (int j = 0; j < TILE_SIZE; j++) {
                    real3 delta = pos-posBuffer[j];
                    delta.x -= floor(delta.x*invPeriodicBoxSize.x+0.5f)*periodicBoxSize.x;
                    delta.y -= floor(delta.y*invPeriodicBoxSize.y+0.5f)*periodicBoxSize.y;
                    delta.z -= floor(delta.z*invPeriodicBoxSize.z+0.5f)*periodicBoxSize.z;
                    interacts |= (delta.x*delta.x+delta.y*delta.y+delta.z*delta.z < PADDED_CUTOFF_SQUARED);
                }
            }
            else {
#endif
                for (int j = 0; j < TILE_SIZE; j++) {
                    real3 delta = pos-posBuffer[j];
                    interacts |= (delta.x*delta.x+delta.y*delta.y+delta.z*delta.z < PADDED_CUTOFF_SQUARED);
                }
#ifdef USE_PERIODIC
            }
#endif
            sum[i*WARP_SIZE+indexInWarp] = (interacts ? 1 : 0);
        }
        for (int i = numValid-base+threadIdx.x/WARP_SIZE; i < BUFFER_SIZE/WARP_SIZE; i += GROUP_SIZE/WARP_SIZE)
            sum[i*WARP_SIZE+indexInWarp] = 0;

        // Compact the list of atoms.

        __syncthreads();
        prefixSum(sum, temp);
        for (int i = threadIdx.x; i < BUFFER_SIZE; i += blockDim.x)
            if (sum[i] != (i == 0 ? 0 : sum[i-1]))
                atoms[numAtoms+sum[i]-1] = buffer[base+i/WARP_SIZE]*TILE_SIZE+indexInWarp;

        // Store them to global memory.

        int atomsToStore = numAtoms+sum[BUFFER_SIZE-1];
        bool storePartialTile = (finish && base >= numValid-BUFFER_SIZE/WARP_SIZE);
        int tilesToStore = (storePartialTile ? (atomsToStore+TILE_SIZE-1)/TILE_SIZE : atomsToStore/TILE_SIZE);
        if (tilesToStore > 0) {
            if (threadIdx.x == 0)
                baseIndex = atomicAdd(interactionCount, tilesToStore);
            __syncthreads();
            if (threadIdx.x == 0)
                numAtoms = atomsToStore-tilesToStore*TILE_SIZE;
            if (baseIndex+tilesToStore <= maxTiles) {
                if (threadIdx.x < tilesToStore)
                    interactingTiles[baseIndex+threadIdx.x] = x;
                for (int i = threadIdx.x; i < tilesToStore*TILE_SIZE; i += blockDim.x)
                    interactingAtoms[baseIndex*TILE_SIZE+i] = (i < atomsToStore ? atoms[i] : NUM_ATOMS);
            }
        }
        else {
            __syncthreads();
            if (threadIdx.x == 0)
                numAtoms += sum[BUFFER_SIZE-1];
        }
        __syncthreads();
        if (threadIdx.x < numAtoms && !storePartialTile)
            atoms[threadIdx.x] = atoms[tilesToStore*TILE_SIZE+threadIdx.x];
    }

    if (numValid == 0 && numAtoms > 0 && finish) {
        // We didn't have any more tiles to process, but there were some atoms left over from a
        // previous call to this function.  Save them now.

        if (threadIdx.x == 0)
            baseIndex = atomicAdd(interactionCount, 1);
        __syncthreads();
        if (baseIndex < maxTiles) {
            if (threadIdx.x == 0)
                interactingTiles[baseIndex] = x;
            if (threadIdx.x < TILE_SIZE)
                interactingAtoms[baseIndex*TILE_SIZE+threadIdx.x] = (threadIdx.x < numAtoms ? atoms[threadIdx.x] : NUM_ATOMS);
        }
    }

    // Reset the buffer for processing more tiles.

    for (int i = threadIdx.x; i < BUFFER_SIZE; i += blockDim.x)
        buffer[i] = INVALID;
    __syncthreads();
}

/**
 * Compare the bounding boxes for each pair of atom blocks (comprised of 32 atoms each), forming a tile. If the two
 * atom blocks are sufficiently far apart, mark them as non-interacting. There are two stages in the algorithm.
 *
 * STAGE 1:
 *
 * A coarse grain atomblock against interacting atomblock neighbourlist is constructed. 
 *
 * Each threadblock first loads in some block X of interest. Each thread within the threadblock then loads 
 * in a different atomblock Y. If Y has exclusions with X, then Y is not processed.  If the bounding boxes 
 * of the two atomblocks are within the cutoff distance, then the two atomblocks are considered to be
 * interacting and Y is added to the buffer for X. If during any given iteration an atomblock (or thread) 
 * finds BUFFER_GROUP interacting blocks, the entire buffer is sent for compaction by storeInteractionData().
 *
 * STAGE 2:
 *
 * A fine grain atomblock against interacting atoms neighbourlist is constructed.
 *
 * The input is an atomblock list detailing the interactions with other atomblocks. The list of interacting 
 * atom blocks are initially stored in the buffer array in shared memory. buffer is then compacted using 
 * prefixSum. Afterwards, each threadblock processes one contiguous atomblock X. Each warp in a threadblock 
 * processes a block Y to find the atoms that interact with any given atom in X. Once BUFFER_SIZE/WARP_SIZE 
 * (eg. 16) atomblocks have been processed for a given X, the list of interacting atoms in these 16 blocks 
 * are subsequently compacted. The process repeats until all atomblocks that interact with X are computed.
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
extern "C" __global__ void findBlocksWithInteractions(real4 periodicBoxSize, real4 invPeriodicBoxSize, unsigned int* __restrict__ interactionCount,
        int* __restrict__ interactingTiles, unsigned int* __restrict__ interactingAtoms, const real4* __restrict__ posq, unsigned int maxTiles, unsigned int startBlockIndex,
        unsigned int numBlocks, real2* __restrict__ sortedBlocks, const real4* __restrict__ sortedBlockCenter, const real4* __restrict__ sortedBlockBoundingBox,
        const unsigned int* __restrict__ exclusionIndices, const unsigned int* __restrict__ exclusionRowIndices, real4* __restrict__ oldPositions,
        const int* __restrict__ rebuildNeighborList) {
    __shared__ unsigned short buffer[BUFFER_SIZE];
    __shared__ short sum[BUFFER_SIZE];
    __shared__ ushort2 temp[BUFFER_SIZE];
    __shared__ int atoms[BUFFER_SIZE+TILE_SIZE];
    __shared__ real3 posBuffer[TILE_SIZE];
    __shared__ int exclusionsForX[MAX_EXCLUSIONS];
    __shared__ int bufferFull;
    __shared__ int globalIndex;
    __shared__ int numAtoms;
    
    if (rebuildNeighborList[0] == 0)
        return; // The neighbor list doesn't need to be rebuilt.
    
    int valuesInBuffer = 0;
    if (threadIdx.x == 0)
        bufferFull = false;
    for (int i = 0; i < BUFFER_GROUPS; ++i)
        buffer[i*GROUP_SIZE+threadIdx.x] = INVALID;
    __syncthreads();
    
    // Loop over blocks sorted by size.
    
    for (int i = startBlockIndex+blockIdx.x; i < startBlockIndex+numBlocks; i += gridDim.x) {
        if (threadIdx.x == blockDim.x-1)
            numAtoms = 0;
        real2 sortedKey = sortedBlocks[i];
        int x = (int) sortedKey.y;
        real4 blockCenterX = sortedBlockCenter[i];
        real4 blockSizeX = sortedBlockBoundingBox[i];

        // Load exclusion data for block x.
        
        const int exclusionStart = exclusionRowIndices[x];
        const int exclusionEnd = exclusionRowIndices[x+1];
        const int numExclusions = exclusionEnd-exclusionStart;
        for (int j = threadIdx.x; j < numExclusions; j += blockDim.x)
            exclusionsForX[j] = exclusionIndices[exclusionStart+j];
        __syncthreads();
        
        // Compare it to other blocks after this one in sorted order.
        
        for (int base = i+1; base < NUM_BLOCKS; base += blockDim.x) {
            int j = base+threadIdx.x;
            real2 sortedKey2 = (j < NUM_BLOCKS ? sortedBlocks[j] : make_real2(0));
            real4 blockCenterY = (j < NUM_BLOCKS ? sortedBlockCenter[j] : make_real4(0));
            real4 blockSizeY = (j < NUM_BLOCKS ? sortedBlockBoundingBox[j] : make_real4(0));
            unsigned short y = (unsigned short) sortedKey2.y;
            real4 delta = blockCenterX-blockCenterY;
#ifdef USE_PERIODIC
            delta.x -= floor(delta.x*invPeriodicBoxSize.x+0.5f)*periodicBoxSize.x;
            delta.y -= floor(delta.y*invPeriodicBoxSize.y+0.5f)*periodicBoxSize.y;
            delta.z -= floor(delta.z*invPeriodicBoxSize.z+0.5f)*periodicBoxSize.z;
#endif
            delta.x = max(0.0f, fabs(delta.x)-blockSizeX.x-blockSizeY.x);
            delta.y = max(0.0f, fabs(delta.y)-blockSizeX.y-blockSizeY.y);
            delta.z = max(0.0f, fabs(delta.z)-blockSizeX.z-blockSizeY.z);
            bool hasExclusions = false;
            for (int k = 0; k < numExclusions; k++)
                hasExclusions |= (exclusionsForX[k] == y);
            if (j < NUM_BLOCKS && delta.x*delta.x+delta.y*delta.y+delta.z*delta.z < PADDED_CUTOFF_SQUARED && !hasExclusions) {
                // Add this tile to the buffer.

                int bufferIndex = valuesInBuffer*GROUP_SIZE+threadIdx.x;
                buffer[bufferIndex] = y;
                valuesInBuffer++;

                // cuda-memcheck --tool racecheck will throw errors about this as 
                // RAW/WAW/WAR race condition errors. But this is safe in all instances
                if (!bufferFull && valuesInBuffer == BUFFER_GROUPS)
                    bufferFull = true;
            }
            __syncthreads();
            if (bufferFull) {
                storeInteractionData(x, buffer, sum, temp, atoms, numAtoms, globalIndex, interactionCount, interactingTiles, interactingAtoms, periodicBoxSize, invPeriodicBoxSize, posq, posBuffer, blockCenterX, blockSizeX, maxTiles, false);
                valuesInBuffer = 0;
                if (threadIdx.x == 0)
                    bufferFull = false;
            }
            __syncthreads();
        }
        storeInteractionData(x, buffer, sum, temp, atoms, numAtoms, globalIndex, interactionCount, interactingTiles, interactingAtoms, periodicBoxSize, invPeriodicBoxSize, posq, posBuffer, blockCenterX, blockSizeX, maxTiles, true);
    }
    
    // Record the positions the neighbor list is based on.
    
    for (int i = threadIdx.x+blockIdx.x*blockDim.x; i < NUM_ATOMS; i += blockDim.x*gridDim.x)
        oldPositions[i] = posq[i];
}
