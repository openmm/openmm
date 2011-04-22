#pragma OPENCL EXTENSION cl_khr_global_int32_base_atomics : enable
#pragma OPENCL EXTENSION cl_khr_byte_addressable_store : enable
#define TILE_SIZE 32
#define GROUP_SIZE 64
#define BUFFER_GROUPS 4
#define BUFFER_SIZE BUFFER_GROUPS*GROUP_SIZE

/**
 * Find a bounding box for the atoms in each block.
 */
__kernel void findBlockBounds(int numAtoms, float4 periodicBoxSize, float4 invPeriodicBoxSize, __global float4* posq, __global float4* blockCenter, __global float4* blockBoundingBox, __global unsigned int* interactionCount) {
    int index = get_global_id(0);
    int base = index*TILE_SIZE;
    while (base < numAtoms) {
        float4 pos = posq[base];
#ifdef USE_PERIODIC
        pos.x -= floor(pos.x*invPeriodicBoxSize.x)*periodicBoxSize.x;
        pos.y -= floor(pos.y*invPeriodicBoxSize.y)*periodicBoxSize.y;
        pos.z -= floor(pos.z*invPeriodicBoxSize.z)*periodicBoxSize.z;
        float4 firstPoint = pos;
#endif
        float4 minPos = pos;
        float4 maxPos = pos;
        int last = min(base+TILE_SIZE, numAtoms);
        for (int i = base+1; i < last; i++) {
            pos = posq[i];
#ifdef USE_PERIODIC
            pos.x -= floor((pos.x-firstPoint.x)*invPeriodicBoxSize.x+0.5f)*periodicBoxSize.x;
            pos.y -= floor((pos.y-firstPoint.y)*invPeriodicBoxSize.y+0.5f)*periodicBoxSize.y;
            pos.z -= floor((pos.z-firstPoint.z)*invPeriodicBoxSize.z+0.5f)*periodicBoxSize.z;
#endif
            minPos = min(minPos, pos);
            maxPos = max(maxPos, pos);
        }
        blockBoundingBox[index] = 0.5f*(maxPos-minPos);
        blockCenter[index] = 0.5f*(maxPos+minPos);
        index += get_global_size(0);
        base = index*TILE_SIZE;
    }
    if (get_global_id(0) == 0)
        interactionCount[0] = 0;
}

/**
 * This is called by findBlocksWithInteractions().  It compacts the list of blocks and writes them
 * to global memory.
 */
void storeInteractionData(__local ushort2* buffer, __local int* valid, __local short* sum, __local ushort2* temp, __local int* baseIndex,
            __global unsigned int* interactionCount, __global ushort2* interactingTiles, float cutoffSquared, float4 periodicBoxSize,
            float4 invPeriodicBoxSize, __global float4* posq, __global float4* blockCenter, __global float4* blockBoundingBox, unsigned int maxTiles) {
    // The buffer is full, so we need to compact it and write out results.  Start by doing a parallel prefix sum.

    for (int i = get_local_id(0); i < BUFFER_SIZE; i += GROUP_SIZE)
        temp[i].x = (valid[i] ? 1 : 0);
    barrier(CLK_LOCAL_MEM_FENCE);
    int whichBuffer = 0;
    for (int offset = 1; offset < BUFFER_SIZE; offset *= 2) {
        if (whichBuffer == 0)
            for (int i = get_local_id(0); i < BUFFER_SIZE; i += GROUP_SIZE)
                temp[i].y = (i < offset ? temp[i].x : temp[i].x+temp[i-offset].x);
        else
            for (int i = get_local_id(0); i < BUFFER_SIZE; i += GROUP_SIZE)
                temp[i].x = (i < offset ? temp[i].y : temp[i].y+temp[i-offset].y);
        whichBuffer = 1-whichBuffer;
        barrier(CLK_LOCAL_MEM_FENCE);
    }
    if (whichBuffer == 0)
        for (int i = get_local_id(0); i < BUFFER_SIZE; i += GROUP_SIZE)
            sum[i] = temp[i].x;
    else
        for (int i = get_local_id(0); i < BUFFER_SIZE; i += GROUP_SIZE)
            sum[i] = temp[i].y;
    barrier(CLK_LOCAL_MEM_FENCE);

    // Compact the buffer.

    for (int i = get_local_id(0); i < BUFFER_SIZE; i += GROUP_SIZE)
        if (valid[i]) {
            temp[sum[i]-1] = buffer[i];
            valid[i] = false;
        }
    barrier(CLK_LOCAL_MEM_FENCE);
    int numValid = sum[BUFFER_SIZE-1];

    // Filter the list of tiles by comparing the distance from each atom to the other bounding box.

    int tile;
    int index = get_local_id(0)&(TILE_SIZE-1);
    int group = get_local_id(0)/TILE_SIZE;
    __local int* flag = sum;
    int lasty = -1;
    float4 center, boxSize, pos;
    for (tile = 0; tile < numValid; tile++) {
        int x = temp[tile].x;
        int y = temp[tile].y;
        if (x == y) {
            tile++;
            continue;
        }
        if (index == 0)
            flag[group] = true;
        barrier(CLK_LOCAL_MEM_FENCE);

        // Load an atom position and the bounding box the other block.

        if (group == 0) {
            center = blockCenter[x];
            boxSize = blockBoundingBox[x];
            if (y != lasty)
                pos = posq[y*TILE_SIZE+index];
        }
        else {
            if (y != lasty) {
                center = blockCenter[y];
                boxSize = blockBoundingBox[y];
            }
            pos = posq[x*TILE_SIZE+index];
        }
        lasty = y;

        // Find the distance of the atom from the bounding box.

        float4 delta = pos-center;
#ifdef USE_PERIODIC
        delta.x -= floor(delta.x*invPeriodicBoxSize.x+0.5f)*periodicBoxSize.x;
        delta.y -= floor(delta.y*invPeriodicBoxSize.y+0.5f)*periodicBoxSize.y;
        delta.z -= floor(delta.z*invPeriodicBoxSize.z+0.5f)*periodicBoxSize.z;
#endif
        delta = max((float4) 0.0f, fabs(delta)-boxSize);
        if (delta.x*delta.x+delta.y*delta.y+delta.z*delta.z < cutoffSquared)
            flag[group] = false;
        barrier(CLK_LOCAL_MEM_FENCE);
        if (flag[0] || flag[1]) {
            // This tile contains no interactions.

            numValid--;
            if (get_local_id(0) == 0)
                temp[tile] = temp[numValid];
        }
        else
            tile++;
        barrier(CLK_LOCAL_MEM_FENCE);
    }

    // Store it to global memory.

    if (get_local_id(0) == 0)
        *baseIndex = atom_add(interactionCount, numValid);
    barrier(CLK_LOCAL_MEM_FENCE);
    if (*baseIndex+numValid <= maxTiles)
        for (int i = get_local_id(0); i < numValid; i += GROUP_SIZE)
            interactingTiles[*baseIndex+i] = temp[i];
    barrier(CLK_LOCAL_MEM_FENCE);
}

/**
 * Compare the bounding boxes for each pair of blocks.  If they are sufficiently far apart,
 * mark them as non-interacting.
 */
__kernel void findBlocksWithInteractions(float cutoffSquared, float4 periodicBoxSize, float4 invPeriodicBoxSize, __global float4* blockCenter,
        __global float4* blockBoundingBox, __global unsigned int* interactionCount, __global ushort2* interactingTiles,
        __global unsigned int* interactionFlags, __global float4* posq, unsigned int maxTiles) {
    __local ushort2 buffer[BUFFER_SIZE];
    __local int valid[BUFFER_SIZE];
    __local short sum[BUFFER_SIZE];
    __local ushort2 temp[BUFFER_SIZE];
    __local int bufferFull;
    __local int globalIndex;
    int valuesInBuffer = 0;
    if (get_local_id(0) == 0)
        bufferFull = false;
    for (int i = 0; i < BUFFER_GROUPS; ++i)
        valid[i*GROUP_SIZE+get_local_id(0)] = false;
    barrier(CLK_LOCAL_MEM_FENCE);
    for (int baseIndex = START_TILE_INDEX+get_group_id(0)*get_local_size(0); baseIndex < END_TILE_INDEX; baseIndex += get_global_size(0)) {
        // Identify the pair of blocks to compare.

        int index = baseIndex+get_local_id(0);
        if (index < END_TILE_INDEX) {
            unsigned int y = (unsigned int) floor(NUM_BLOCKS+0.5f-sqrt((NUM_BLOCKS+0.5f)*(NUM_BLOCKS+0.5f)-2*index));
            unsigned int x = (index-y*NUM_BLOCKS+y*(y+1)/2);
            if (x < y || x >= NUM_BLOCKS) { // Occasionally happens due to roundoff error.
                y += (x < y ? -1 : 1);
                x = (index-y*NUM_BLOCKS+y*(y+1)/2);
            }

            // Find the distance between the bounding boxes of the two cells.

            float4 delta = blockCenter[x]-blockCenter[y];
#ifdef USE_PERIODIC
            delta.x -= floor(delta.x*invPeriodicBoxSize.x+0.5f)*periodicBoxSize.x;
            delta.y -= floor(delta.y*invPeriodicBoxSize.y+0.5f)*periodicBoxSize.y;
            delta.z -= floor(delta.z*invPeriodicBoxSize.z+0.5f)*periodicBoxSize.z;
#endif
            float4 boxSizea = blockBoundingBox[x];
            float4 boxSizeb = blockBoundingBox[y];
            delta.x = max(0.0f, fabs(delta.x)-boxSizea.x-boxSizeb.x);
            delta.y = max(0.0f, fabs(delta.y)-boxSizea.y-boxSizeb.y);
            delta.z = max(0.0f, fabs(delta.z)-boxSizea.z-boxSizeb.z);
            if (delta.x*delta.x+delta.y*delta.y+delta.z*delta.z < cutoffSquared) {
                // Add this tile to the buffer.

                int bufferIndex = valuesInBuffer*GROUP_SIZE+get_local_id(0);
                valid[bufferIndex] = true;
                buffer[bufferIndex] = (ushort2) (x, y);
                valuesInBuffer++;
                if (!bufferFull && valuesInBuffer == BUFFER_GROUPS)
                    bufferFull = true;
            }
        }
        barrier(CLK_LOCAL_MEM_FENCE);
        if (bufferFull) {
            storeInteractionData(buffer, valid, sum, temp, &globalIndex, interactionCount, interactingTiles, cutoffSquared, periodicBoxSize, invPeriodicBoxSize, posq, blockCenter, blockBoundingBox, maxTiles);
            valuesInBuffer = 0;
            if (get_local_id(0) == 0)
                bufferFull = false;
            barrier(CLK_LOCAL_MEM_FENCE);
        }
    }
    storeInteractionData(buffer, valid, sum, temp, &globalIndex, interactionCount, interactingTiles, cutoffSquared, periodicBoxSize, invPeriodicBoxSize, posq, blockCenter, blockBoundingBox, maxTiles);
}

/**
 * Compare each atom in one block to the bounding box of another block, and set
 * flags for which ones are interacting.
 */
__kernel void findInteractionsWithinBlocks(float cutoffSquared, float4 periodicBoxSize, float4 invPeriodicBoxSize, __global float4* posq, __global ushort2* tiles, __global float4* blockCenter,
            __global float4* blockBoundingBox, __global unsigned int* interactionFlags, __global unsigned int* interactionCount, __local unsigned int* flags, unsigned int maxTiles) {
    unsigned int totalWarps = get_global_size(0)/TILE_SIZE;
    unsigned int warp = get_global_id(0)/TILE_SIZE;
    unsigned int numTiles = interactionCount[0];
    unsigned int pos = warp*numTiles/totalWarps;
    unsigned int end = (warp+1)*numTiles/totalWarps;
    unsigned int index = get_local_id(0) & (TILE_SIZE - 1);

    if (numTiles > maxTiles)
        return;
    unsigned int lasty = 0xFFFFFFFF;
    float4 apos;
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

            float4 center = blockCenter[x];
            float4 boxSize = blockBoundingBox[x];
            if (y != lasty)
                apos = posq[y*TILE_SIZE+index];

            // Find the distance of the atom from the bounding box.

            float4 delta = apos-center;
#ifdef USE_PERIODIC
            delta.x -= floor(delta.x*invPeriodicBoxSize.x+0.5f)*periodicBoxSize.x;
            delta.y -= floor(delta.y*invPeriodicBoxSize.y+0.5f)*periodicBoxSize.y;
            delta.z -= floor(delta.z*invPeriodicBoxSize.z+0.5f)*periodicBoxSize.z;
#endif
            delta = max((float4) 0.0f, fabs(delta)-boxSize);
            int thread = get_local_id(0);
            flags[thread] = (delta.x*delta.x+delta.y*delta.y+delta.z*delta.z > cutoffSquared ? 0 : 1 << index);

            // Sum the flags.

#ifdef WARPS_ARE_ATOMIC
            if (index % 2 == 0)
                flags[thread] += flags[thread+1];
            if (index % 4 == 0)
                flags[thread] += flags[thread+2];
            if (index % 8 == 0)
                flags[thread] += flags[thread+4];
            if (index % 16 == 0)
                flags[thread] += flags[thread+8];
#else
            barrier(CLK_LOCAL_MEM_FENCE);
            if (index % 2 == 0)
                flags[thread] += flags[thread+1];
            barrier(CLK_LOCAL_MEM_FENCE);
            if (index % 4 == 0)
                flags[thread] += flags[thread+2];
            barrier(CLK_LOCAL_MEM_FENCE);
            if (index % 8 == 0)
                flags[thread] += flags[thread+4];
            barrier(CLK_LOCAL_MEM_FENCE);
            if (index % 16 == 0)
                flags[thread] += flags[thread+8];
            barrier(CLK_LOCAL_MEM_FENCE);
#endif
            if (index == 0) {
                unsigned int allFlags = flags[thread] + flags[thread+16];

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
