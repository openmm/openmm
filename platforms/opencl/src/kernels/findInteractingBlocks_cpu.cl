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
void storeInteractionData(ushort2* buffer, int numValid, __global unsigned int* interactionCount, __global ushort2* interactingTiles,
            __global unsigned int* interactionFlags, float cutoffSquared, float4 periodicBoxSize, float4 invPeriodicBoxSize,
            __global float4* posq, __global float4* blockCenter, __global float4* blockBoundingBox, unsigned int maxTiles) {
    // Filter the list of tiles by comparing the distance from each atom to the other bounding box.

    unsigned int flagsBuffer[2*BUFFER_SIZE];
    float4 atomPositions[TILE_SIZE];
    int lasty = -1;
    float4 centery, boxSizey;
    for (int tile = 0; tile < numValid; ) {
        int x = buffer[tile].x;
        int y = buffer[tile].y;
        if (x == y) {
            tile++;
            continue;
        }

        // Load the atom positions and bounding boxes.

        float4 centerx = blockCenter[x];
        float4 boxSizex = blockBoundingBox[x];
        if (y != lasty) {
            for (int atom = 0; atom < TILE_SIZE; atom++)
                atomPositions[atom] = posq[y*TILE_SIZE+atom];
            centery = blockCenter[y];
            boxSizey = blockBoundingBox[y];
            lasty = y;
        }

        // Find the distance of each atom from the bounding box.

        unsigned int flags1 = 0, flags2 = 0;
        for (int atom = 0; atom < TILE_SIZE; atom++) {
            float4 delta = atomPositions[atom]-centerx;
#ifdef USE_PERIODIC
            delta.xyz -= floor(delta.xyz*invPeriodicBoxSize.xyz+0.5f)*periodicBoxSize.xyz;
#endif
            delta = max((float4) 0.0f, fabs(delta)-boxSizex);
            if (dot(delta.xyz, delta.xyz) < cutoffSquared)
                flags1 += 1 << atom;
            delta = posq[x*TILE_SIZE+atom]-centery;
#ifdef USE_PERIODIC
            delta.xyz -= floor(delta.xyz*invPeriodicBoxSize.xyz+0.5f)*periodicBoxSize.xyz;
#endif
            delta = max((float4) 0.0f, fabs(delta)-boxSizey);
            if (dot(delta.xyz, delta.xyz) < cutoffSquared)
                flags2 += 1 << atom;
        }
        if (flags1 == 0 || flags2 == 0) {
            // This tile contains no interactions.

            numValid--;
            buffer[tile] = buffer[numValid];
        }
        else {
            flagsBuffer[2*tile] = flags1;
            flagsBuffer[2*tile+1] = flags2;
            tile++;
        }
    }

    // Store it to global memory.

    int baseIndex = atom_add(interactionCount, numValid);
    if (baseIndex+numValid <= maxTiles)
        for (int i = 0; i < numValid; i++) {
            interactingTiles[baseIndex+i] = buffer[i];
            interactionFlags[2*(baseIndex+i)] = flagsBuffer[2*i];
            interactionFlags[2*(baseIndex+i)+1] = flagsBuffer[2*i+1];
        }
}

/**
 * Compare the bounding boxes for each pair of blocks.  If they are sufficiently far apart,
 * mark them as non-interacting.
 */
__kernel void findBlocksWithInteractions(float cutoffSquared, float4 periodicBoxSize, float4 invPeriodicBoxSize, __global float4* blockCenter,
        __global float4* blockBoundingBox, __global unsigned int* interactionCount, __global ushort2* interactingTiles,
        __global unsigned int* interactionFlags, __global float4* posq, unsigned int maxTiles) {
    ushort2 buffer[BUFFER_SIZE];
    int valuesInBuffer = 0;
    const int numTiles = END_TILE_INDEX-START_TILE_INDEX;
    unsigned int start = START_TILE_INDEX+get_group_id(0)*numTiles/get_num_groups(0);
    unsigned int end = START_TILE_INDEX+(get_group_id(0)+1)*numTiles/get_num_groups(0);
    for (int index = start; index < end; index++) {
        // Identify the pair of blocks to compare.

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

            buffer[valuesInBuffer++] = (ushort2) (x, y);
            if (valuesInBuffer == BUFFER_SIZE) {
                storeInteractionData(buffer, valuesInBuffer, interactionCount, interactingTiles, interactionFlags, cutoffSquared, periodicBoxSize, invPeriodicBoxSize, posq, blockCenter, blockBoundingBox, maxTiles);
                valuesInBuffer = 0;
            }
        }
    }
    storeInteractionData(buffer, valuesInBuffer, interactionCount, interactingTiles, interactionFlags, cutoffSquared, periodicBoxSize, invPeriodicBoxSize, posq, blockCenter, blockBoundingBox, maxTiles);
}
