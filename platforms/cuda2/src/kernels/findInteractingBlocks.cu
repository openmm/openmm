#define TILE_SIZE 32
#define GROUP_SIZE 64
#define BUFFER_GROUPS 4
#define BUFFER_SIZE BUFFER_GROUPS*GROUP_SIZE

/**
 * Find a bounding box for the atoms in each block.
 */
extern "C" __global__ void findBlockBounds(int numAtoms, real4 periodicBoxSize, real4 invPeriodicBoxSize, const real4* __restrict__ posq, real4* __restrict__ blockCenter, real4* __restrict__ blockBoundingBox, unsigned int* __restrict__ interactionCount) {
    int index = blockIdx.x*blockDim.x+threadIdx.x;
    int base = index*TILE_SIZE;
    while (base < numAtoms) {
        real4 pos = posq[base];
#ifdef USE_PERIODIC
        pos.x -= floor(pos.x*invPeriodicBoxSize.x)*periodicBoxSize.x;
        pos.y -= floor(pos.y*invPeriodicBoxSize.y)*periodicBoxSize.y;
        pos.z -= floor(pos.z*invPeriodicBoxSize.z)*periodicBoxSize.z;
        real4 firstPoint = pos;
#endif
        real4 minPos = pos;
        real4 maxPos = pos;
        int last = min(base+TILE_SIZE, numAtoms);
        for (int i = base+1; i < last; i++) {
            pos = posq[i];
#ifdef USE_PERIODIC
            pos.x -= floor((pos.x-firstPoint.x)*invPeriodicBoxSize.x+0.5f)*periodicBoxSize.x;
            pos.y -= floor((pos.y-firstPoint.y)*invPeriodicBoxSize.y+0.5f)*periodicBoxSize.y;
            pos.z -= floor((pos.z-firstPoint.z)*invPeriodicBoxSize.z+0.5f)*periodicBoxSize.z;
#endif
            minPos = make_real4(min(minPos.x,pos.x), min(minPos.y,pos.y), min(minPos.z,pos.z), 0);
            maxPos = make_real4(max(maxPos.x,pos.x), max(maxPos.y,pos.y), max(maxPos.z,pos.z), 0);
        }
        blockBoundingBox[index] = 0.5f*(maxPos-minPos);
        blockCenter[index] = 0.5f*(maxPos+minPos);
        index += blockDim.x*gridDim.x;
        base = index*TILE_SIZE;
    }
    if (blockIdx.x == 0 && threadIdx.x == 0)
        interactionCount[0] = 0;
}

/**
 * This is called by findBlocksWithInteractions().  It compacts the list of blocks and writes them
 * to global memory.
 */
__device__ void storeInteractionData(ushort2* buffer, int* valid, short* sum, ushort2* temp, int* baseIndex,
            unsigned int* interactionCount, ushort2* interactingTiles, real4 periodicBoxSize,
            real4 invPeriodicBoxSize, const real4* posq, const real4* blockCenter, const real4* blockBoundingBox, unsigned int maxTiles) {
    // The buffer is full, so we need to compact it and write out results.  Start by doing a parallel prefix sum.

    for (int i = threadIdx.x; i < BUFFER_SIZE; i += GROUP_SIZE)
        temp[i].x = (valid[i] ? 1 : 0);
    __syncthreads();
    int whichBuffer = 0;
    for (int offset = 1; offset < BUFFER_SIZE; offset *= 2) {
        if (whichBuffer == 0)
            for (int i = threadIdx.x; i < BUFFER_SIZE; i += GROUP_SIZE)
                temp[i].y = (i < offset ? temp[i].x : temp[i].x+temp[i-offset].x);
        else
            for (int i = threadIdx.x; i < BUFFER_SIZE; i += GROUP_SIZE)
                temp[i].x = (i < offset ? temp[i].y : temp[i].y+temp[i-offset].y);
        whichBuffer = 1-whichBuffer;
        __syncthreads();
    }
    if (whichBuffer == 0)
        for (int i = threadIdx.x; i < BUFFER_SIZE; i += GROUP_SIZE)
            sum[i] = temp[i].x;
    else
        for (int i = threadIdx.x; i < BUFFER_SIZE; i += GROUP_SIZE)
            sum[i] = temp[i].y;
    __syncthreads();
    int numValid = sum[BUFFER_SIZE-1];
    __syncthreads();

    // Compact the buffer.

    for (int i = threadIdx.x; i < BUFFER_SIZE; i += GROUP_SIZE)
        if (valid[i]) {
            temp[sum[i]-1] = buffer[i];
            sum[i] = valid[i];
            valid[i] = false;
            buffer[i] = make_ushort2(1, 1);
        }
    __syncthreads();

    // Store it to global memory.

    if (threadIdx.x == 0)
        *baseIndex = atomicAdd(interactionCount, numValid);
    __syncthreads();
    if (*baseIndex+numValid <= maxTiles)
        for (int i = threadIdx.x; i < numValid; i += GROUP_SIZE)
            interactingTiles[*baseIndex+i] = temp[i];
    __syncthreads();
}

/**
 * Compare the bounding boxes for each pair of blocks.  If they are sufficiently far apart,
 * mark them as non-interacting.
 */
extern "C" __global__ void findBlocksWithInteractions(real4 periodicBoxSize, real4 invPeriodicBoxSize, const real4* __restrict__ blockCenter,
        const real4* __restrict__ blockBoundingBox, unsigned int* __restrict__ interactionCount, ushort2* __restrict__ interactingTiles,
        unsigned int* __restrict__ interactionFlags, const real4* __restrict__ posq, unsigned int maxTiles, unsigned int startTileIndex,
        unsigned int endTileIndex) {
    __shared__ ushort2 buffer[BUFFER_SIZE];
    __shared__ int valid[BUFFER_SIZE];
    __shared__ short sum[BUFFER_SIZE];
    __shared__ ushort2 temp[BUFFER_SIZE];
    __shared__ int bufferFull;
    __shared__ int globalIndex;
    int valuesInBuffer = 0;
    if (threadIdx.x == 0)
        bufferFull = false;
    for (int i = 0; i < BUFFER_GROUPS; ++i)
        valid[i*GROUP_SIZE+threadIdx.x] = false;
    __syncthreads();
    for (int baseIndex = startTileIndex+blockIdx.x*blockDim.x; baseIndex < endTileIndex; baseIndex += blockDim.x*gridDim.x) {
        // Identify the pair of blocks to compare.

        int index = baseIndex+threadIdx.x;
        if (index < endTileIndex) {
            unsigned int y = (unsigned int) floor(NUM_BLOCKS+0.5f-sqrt((NUM_BLOCKS+0.5f)*(NUM_BLOCKS+0.5f)-2*index));
            unsigned int x = (index-y*NUM_BLOCKS+y*(y+1)/2);
            if (x < y || x >= NUM_BLOCKS) { // Occasionally happens due to roundoff error.
                y += (x < y ? -1 : 1);
                x = (index-y*NUM_BLOCKS+y*(y+1)/2);
            }

            // Find the distance between the bounding boxes of the two cells.

            real4 delta = blockCenter[x]-blockCenter[y];
            real4 boxSizea = blockBoundingBox[x];
            real4 boxSizeb = blockBoundingBox[y];
#ifdef USE_PERIODIC
            delta.x -= floor(delta.x*invPeriodicBoxSize.x+0.5f)*periodicBoxSize.x;
            delta.y -= floor(delta.y*invPeriodicBoxSize.y+0.5f)*periodicBoxSize.y;
            delta.z -= floor(delta.z*invPeriodicBoxSize.z+0.5f)*periodicBoxSize.z;
#endif
            delta.x = max(0.0f, fabs(delta.x)-boxSizea.x-boxSizeb.x);
            delta.y = max(0.0f, fabs(delta.y)-boxSizea.y-boxSizeb.y);
            delta.z = max(0.0f, fabs(delta.z)-boxSizea.z-boxSizeb.z);
            if (delta.x*delta.x+delta.y*delta.y+delta.z*delta.z < CUTOFF_SQUARED) {
                // Add this tile to the buffer.

                int bufferIndex = valuesInBuffer*GROUP_SIZE+threadIdx.x;
                valid[bufferIndex] = true;
                buffer[bufferIndex] = make_ushort2(x, y);
                valuesInBuffer++;
                if (!bufferFull && valuesInBuffer == BUFFER_GROUPS)
                    bufferFull = true;
            }
        }
        __syncthreads();
        if (bufferFull) {
            storeInteractionData(buffer, valid, sum, temp, &globalIndex, interactionCount, interactingTiles, periodicBoxSize, invPeriodicBoxSize, posq, blockCenter, blockBoundingBox, maxTiles);
            valuesInBuffer = 0;
            if (threadIdx.x == 0)
                bufferFull = false;
            __syncthreads();
        }
    }
    storeInteractionData(buffer, valid, sum, temp, &globalIndex, interactionCount, interactingTiles, periodicBoxSize, invPeriodicBoxSize, posq, blockCenter, blockBoundingBox, maxTiles);
}

/**
 * Compare each atom in one block to the bounding box of another block, and set
 * flags for which ones are interacting.
 */
extern "C" __global__ void findInteractionsWithinBlocks(real4 periodicBoxSize, real4 invPeriodicBoxSize, const real4* __restrict__ posq, const ushort2* __restrict__ tiles, const real4* __restrict__ blockCenter,
            const real4* __restrict__ blockBoundingBox, unsigned int* __restrict__ interactionFlags, const unsigned int* __restrict__ interactionCount, unsigned int maxTiles) {
    unsigned int totalWarps = (blockDim.x*gridDim.x)/TILE_SIZE;
    unsigned int warp = (blockIdx.x*blockDim.x+threadIdx.x)/TILE_SIZE;
    unsigned int numTiles = interactionCount[0];
    unsigned int pos = warp*numTiles/totalWarps;
    unsigned int end = (warp+1)*numTiles/totalWarps;
    unsigned int index = threadIdx.x & (TILE_SIZE - 1);
#if (__CUDA_ARCH__ < 200)
    __shared__ unsigned int flags[128];
#endif

    if (numTiles > maxTiles)
        return;
    unsigned int lasty = 0xFFFFFFFF;
    real4 apos;
    while (pos < end) {
        // Extract the coordinates of this tile
        ushort2 tileIndices = tiles[pos];
        unsigned int x = tileIndices.x;
        unsigned int y = tileIndices.y;
        if (x == y) {
            if (index == 0)
                interactionFlags[pos] = 0xFFFFFFFF;
        }
        else {
            // Load the bounding box for x and the atom positions for y.

            real4 center = blockCenter[x];
            real4 boxSize = blockBoundingBox[x];
            if (y != lasty)
                apos = posq[y*TILE_SIZE+index];

            // Find the distance of the atom from the bounding box.

            real4 delta = apos-center;
#ifdef USE_PERIODIC
                delta.x -= floor(delta.x*invPeriodicBoxSize.x+0.5f)*periodicBoxSize.x;
                delta.y -= floor(delta.y*invPeriodicBoxSize.y+0.5f)*periodicBoxSize.y;
                delta.z -= floor(delta.z*invPeriodicBoxSize.z+0.5f)*periodicBoxSize.z;
#endif
            delta.x = max((real) 0, fabs(delta.x)-boxSize.x);
            delta.y = max((real) 0, fabs(delta.y)-boxSize.y);
            delta.z = max((real) 0, fabs(delta.z)-boxSize.z);
#if (__CUDA_ARCH__ < 200)
            flags[threadIdx.x] = (delta.x*delta.x+delta.y*delta.y+delta.z*delta.z > CUTOFF_SQUARED ? 0 : 1 << index);
            if (index % 4 == 0)
                flags[threadIdx.x] += flags[threadIdx.x+1]+flags[threadIdx.x+2]+flags[threadIdx.x+3];
            unsigned int allFlags = 0;
            if (index == 0)
                allFlags = flags[threadIdx.x]+flags[threadIdx.x+4]+flags[threadIdx.x+8]+flags[threadIdx.x+12]+flags[threadIdx.x+16]+flags[threadIdx.x+20]+flags[threadIdx.x+24]+flags[threadIdx.x+28];
#else
            unsigned int allFlags = __ballot(delta.x*delta.x+delta.y*delta.y+delta.z*delta.z < CUTOFF_SQUARED);
#endif

            // Sum the flags.

            if (index == 0) {
                // Count how many flags are set, and based on that decide whether to compute all interactions
                // or only a fraction of them.

                unsigned int bits = (allFlags&0x55555555) + ((allFlags>>1)&0x55555555);
                bits = (bits&0x33333333) + ((bits>>2)&0x33333333);
                bits = (bits&0x0F0F0F0F) + ((bits>>4)&0x0F0F0F0F);
                bits = (bits&0x00FF00FF) + ((bits>>8)&0x00FF00FF);
                bits = (bits&0x0000FFFF) + ((bits>>16)&0x0000FFFF);
                interactionFlags[pos] = (bits > 12 ? 0xFFFFFFFF : allFlags);
            }
            lasty = y;
        }
        pos++;
    }
}
