const int TileSize = 32;

/**
 * Find a bounding box for the atoms in each block.
 */
__kernel void findBlockBounds(int numAtoms, float4 periodicBoxSize, __global float4* posq, __global float4* blockCenter, __global float4* blockBoundingBox) {
    int index = get_global_id(0);
    int base = index*TileSize;
    while (base < numAtoms) {
        float4 pos = posq[base];
#ifdef USE_PERIODIC
        pos.x -= floor(pos.x/periodicBoxSize.x)*periodicBoxSize.x;
        pos.y -= floor(pos.y/periodicBoxSize.y)*periodicBoxSize.y;
        pos.z -= floor(pos.z/periodicBoxSize.z)*periodicBoxSize.z;
        float4 firstPoint = pos;
#endif
        float4 minPos = pos;
        float4 maxPos = pos;
        int last = min(base+TileSize, numAtoms);
        for (int i = base+1; i < last; i++) {
            pos = posq[i];
#ifdef USE_PERIODIC
            pos.x -= floor((pos.x-firstPoint.x)/periodicBoxSize.x+0.5f)*periodicBoxSize.x;
            pos.y -= floor((pos.y-firstPoint.y)/periodicBoxSize.y+0.5f)*periodicBoxSize.y;
            pos.z -= floor((pos.z-firstPoint.z)/periodicBoxSize.z+0.5f)*periodicBoxSize.z;
#endif
            minPos = min(minPos, pos);
            maxPos = max(maxPos, pos);
        }
        blockBoundingBox[index] = 0.5f*(maxPos-minPos);
        blockCenter[index] = 0.5f*(maxPos+minPos);
        index += get_global_size(0);
        base = index*TileSize;
    }
}

/**
 * Compare the bounding boxes for each pair of blocks.  If they are sufficiently far apart,
 * mark them as non-interacting.
 */
__kernel void findBlocksWithInteractions(int numTiles, float cutoffSquared, float4 periodicBoxSize, __global unsigned int* tiles, __global float4* blockCenter,
        __global float4* blockBoundingBox, __global unsigned int* interactionFlag) {
    int index = get_global_id(0);
    while (index < numTiles) {
        // Extract cell coordinates from appropriate work unit

        unsigned int x = tiles[index];
        unsigned int y = ((x >> 2) & 0x7fff);
        x = (x >> 17);

        // Find the distance between the bounding boxes of the two cells.

        float4 delta = blockCenter[x]-blockCenter[y];
#ifdef USE_PERIODIC
        delta.x -= floor(delta.x/periodicBoxSize.x+0.5f)*periodicBoxSize.x;
        delta.y -= floor(delta.y/periodicBoxSize.y+0.5f)*periodicBoxSize.y;
        delta.z -= floor(delta.z/periodicBoxSize.z+0.5f)*periodicBoxSize.z;
#endif
        float4 boxSizea = blockBoundingBox[x];
        float4 boxSizeb = blockBoundingBox[y];
        delta.x = max(0.0f, fabs(delta.x)-boxSizea.x-boxSizeb.x);
        delta.y = max(0.0f, fabs(delta.y)-boxSizea.y-boxSizeb.y);
        delta.z = max(0.0f, fabs(delta.z)-boxSizea.z-boxSizeb.z);
        interactionFlag[index] = (delta.x*delta.x+delta.y*delta.y+delta.z*delta.z > cutoffSquared ? 0 : 1);
        index += get_global_size(0);
    }
}

/**
 * Compare each atom in one block to the bounding box of another block, and set
 * flags for which ones are interacting.
 */
__kernel void findInteractionsWithinBlocks(float cutoffSquared, float4 periodicBoxSize, __global float4* posq, __global unsigned int* tiles, __global float4* blockCenter,
            __global float4* blockBoundingBox, __global unsigned int* interactionFlags, __global unsigned int* interactionCount, __local unsigned int* flags) {
    unsigned int totalWarps = get_global_size(0)/TileSize;
    unsigned int warp = get_global_id(0)/TileSize;
    unsigned int numTiles = interactionCount[0];
    unsigned int pos = warp*numTiles/totalWarps;
    unsigned int end = (warp+1)*numTiles/totalWarps;
    unsigned int index = get_local_id(0) & (TileSize - 1);

    unsigned int lasty = 0xFFFFFFFF;
    float4 apos;
    while (pos < end) {
        // Extract the coordinates of this tile
        unsigned int x = tiles[pos];
        unsigned int y = ((x >> 2) & 0x7fff);
        bool hasExclusions = (x & 0x1);
        x = (x >> 17);
        if (x == y || hasExclusions) {
            // Assume this tile will be dense.

            if (index == 0)
                interactionFlags[pos] = 0xFFFFFFFF;
        }
        else {
            // Load the bounding box for x and the atom positions for y.

            float4 center = blockCenter[x];
            float4 boxSize = blockBoundingBox[x];
            if (y != lasty)
                apos = posq[y*TileSize+index];

            // Find the distance of the atom from the bounding box.

            float4 delta = apos-center;
#ifdef USE_PERIODIC
            delta.x -= floor(delta.x/periodicBoxSize.x+0.5f)*periodicBoxSize.x;
            delta.y -= floor(delta.y/periodicBoxSize.y+0.5f)*periodicBoxSize.y;
            delta.z -= floor(delta.z/periodicBoxSize.z+0.5f)*periodicBoxSize.z;
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
