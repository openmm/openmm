#pragma OPENCL EXTENSION cl_khr_global_int32_base_atomics : enable
#define TILE_SIZE 32
#define STORE_DERIVATIVE_1(INDEX) derivBuffers##INDEX[offset] += deriv##INDEX##_1;
#define STORE_DERIVATIVE_2(INDEX) derivBuffers##INDEX[offset] += local_deriv##INDEX[get_local_id(0)];

/**
 * Compute a force based on pair interactions.
 */
__kernel void computeN2Energy(__global float4* forceBuffers, __global float* energyBuffer, __local float4* local_force,
	__global float4* posq, __local float4* local_posq, __global unsigned int* exclusions, __global unsigned int* exclusionIndices,
        __global unsigned int* exclusionRowIndices, __local float4* tempBuffer,
#ifdef USE_CUTOFF
        __global ushort2* tiles, __global unsigned int* interactionCount, float4 periodicBoxSize, float4 invPeriodicBoxSize, unsigned int maxTiles, __global unsigned int* interactionFlags
#else
        unsigned int numTiles
#endif
        PARAMETER_ARGUMENTS) {
    unsigned int totalWarps = get_global_size(0)/TILE_SIZE;
    unsigned int warp = get_global_id(0)/TILE_SIZE;
#ifdef USE_CUTOFF
    unsigned int numTiles = interactionCount[0];
    unsigned int pos = warp*(numTiles > maxTiles ? NUM_BLOCKS*(NUM_BLOCKS+1)/2 : numTiles)/totalWarps;
    unsigned int end = (warp+1)*(numTiles > maxTiles ? NUM_BLOCKS*(NUM_BLOCKS+1)/2 : numTiles)/totalWarps;
#else
    unsigned int pos = warp*numTiles/totalWarps;
    unsigned int end = (warp+1)*numTiles/totalWarps;
#endif
    float energy = 0.0f;
    unsigned int lasty = 0xFFFFFFFF;
    __local unsigned int exclusionRange[2*WARPS_PER_GROUP];
    __local int exclusionIndex[WARPS_PER_GROUP];
    __local int2* reservedBlocks = (__local int2*) exclusionRange;
    
    do {
        // Extract the coordinates of this tile
        const unsigned int tgx = get_local_id(0) & (TILE_SIZE-1);
        const unsigned int tbx = get_local_id(0) - tgx;
        const unsigned int localGroupIndex = get_local_id(0)/TILE_SIZE;
        unsigned int x, y;
        float4 force = 0.0f;
        DECLARE_ATOM1_DERIVATIVES
        if (pos < end) {
#ifdef USE_CUTOFF
            if (numTiles <= maxTiles) {
                ushort2 tileIndices = tiles[pos];
                x = tileIndices.x;
                y = tileIndices.y;
            }
            else
#endif
            {
                y = (unsigned int) floor(NUM_BLOCKS+0.5f-SQRT((NUM_BLOCKS+0.5f)*(NUM_BLOCKS+0.5f)-2*pos));
                x = (pos-y*NUM_BLOCKS+y*(y+1)/2);
                if (x < y || x >= NUM_BLOCKS) { // Occasionally happens due to roundoff error.
                    y += (x < y ? -1 : 1);
                    x = (pos-y*NUM_BLOCKS+y*(y+1)/2);
                }
            }
            unsigned int atom1 = x*TILE_SIZE + tgx;
            float4 posq1 = posq[atom1];
            LOAD_ATOM1_PARAMETERS

            // Locate the exclusion data for this tile.

#ifdef USE_EXCLUSIONS
            if (tgx < 2)
                exclusionRange[2*localGroupIndex+tgx] = exclusionRowIndices[x+tgx];
            if (tgx == 0)
                exclusionIndex[localGroupIndex] = -1;
            for (int i = exclusionRange[2*localGroupIndex]+tgx; i < exclusionRange[2*localGroupIndex+1]; i += TILE_SIZE)
                if (exclusionIndices[i] == y)
                    exclusionIndex[localGroupIndex] = i*TILE_SIZE;
            bool hasExclusions = (exclusionIndex[localGroupIndex] > -1);
#else
            bool hasExclusions = false;
#endif
            if (pos >= end)
                ; // This warp is done.
            else if (x == y) {
                // This tile is on the diagonal.

                const unsigned int localAtomIndex = get_local_id(0);
                local_posq[localAtomIndex] = posq1;
                LOAD_LOCAL_PARAMETERS_FROM_1
#ifdef USE_EXCLUSIONS
                unsigned int excl = exclusions[exclusionIndex[localGroupIndex]+tgx];
#endif
                for (unsigned int j = 0; j < TILE_SIZE; j++) {
#ifdef USE_EXCLUSIONS
                    bool isExcluded = !(excl & 0x1);
#endif
                    int atom2 = tbx+j;
                    float4 posq2 = local_posq[atom2];
                    float4 delta = (float4) (posq2.xyz - posq1.xyz, 0.0f);
#ifdef USE_PERIODIC
                    delta.x -= floor(delta.x*invPeriodicBoxSize.x+0.5f)*periodicBoxSize.x;
                    delta.y -= floor(delta.y*invPeriodicBoxSize.y+0.5f)*periodicBoxSize.y;
                    delta.z -= floor(delta.z*invPeriodicBoxSize.z+0.5f)*periodicBoxSize.z;
#endif
                    float r2 = delta.x*delta.x + delta.y*delta.y + delta.z*delta.z;
#ifdef USE_CUTOFF
                    if (r2 < CUTOFF_SQUARED) {
#endif
                    float r = SQRT(r2);
                    LOAD_ATOM2_PARAMETERS
                    atom2 = y*TILE_SIZE+j;
                    float dEdR = 0.0f;
                    float tempEnergy = 0.0f;
                    if (atom1 < NUM_ATOMS && atom2 < NUM_ATOMS && atom1 != atom2) {
                        COMPUTE_INTERACTION
                        dEdR /= -r;
                    }
                    energy += 0.5f*tempEnergy;
                    delta.xyz *= dEdR;
                    force.xyz -= delta.xyz;
#ifdef USE_CUTOFF
                    }
#endif
#ifdef USE_EXCLUSIONS
                    excl >>= 1;
#endif
                }
            }
            else {
                // This is an off-diagonal tile.

                const unsigned int localAtomIndex = get_local_id(0);
                if (lasty != y) {
                    unsigned int j = y*TILE_SIZE + tgx;
                    local_posq[localAtomIndex] = posq[j];
                    LOAD_LOCAL_PARAMETERS_FROM_GLOBAL
                }
                local_force[localAtomIndex] = 0.0f;
                CLEAR_LOCAL_DERIVATIVES
#ifdef USE_CUTOFF
                unsigned int flags = (numTiles <= maxTiles ? interactionFlags[pos] : 0xFFFFFFFF);
                if (!hasExclusions && flags == 0) {
                    // No interactions in this tile.
                }
                else
#endif
                {
                    // Compute the full set of interactions in this tile.

#ifdef USE_EXCLUSIONS
                    unsigned int excl = (hasExclusions ? exclusions[exclusionIndex[localGroupIndex]+tgx] : 0xFFFFFFFF);
                    excl = (excl >> tgx) | (excl << (TILE_SIZE - tgx));
#endif
                    unsigned int tj = tgx;
                    for (unsigned int j = 0; j < TILE_SIZE; j++) {
#ifdef USE_EXCLUSIONS
                        bool isExcluded = !(excl & 0x1);
#endif
                        int atom2 = tbx+tj;
                        float4 posq2 = local_posq[atom2];
                        float4 delta = (float4) (posq2.xyz - posq1.xyz, 0.0f);
#ifdef USE_PERIODIC
                        delta.x -= floor(delta.x*invPeriodicBoxSize.x+0.5f)*periodicBoxSize.x;
                        delta.y -= floor(delta.y*invPeriodicBoxSize.y+0.5f)*periodicBoxSize.y;
                        delta.z -= floor(delta.z*invPeriodicBoxSize.z+0.5f)*periodicBoxSize.z;
#endif
                        float r2 = delta.x*delta.x + delta.y*delta.y + delta.z*delta.z;
#ifdef USE_CUTOFF
                        if (r2 < CUTOFF_SQUARED) {
#endif
                        float r = SQRT(r2);
                        LOAD_ATOM2_PARAMETERS
                        atom2 = y*TILE_SIZE+tj;
                        float dEdR = 0.0f;
                        float tempEnergy = 0.0f;
                        if (atom1 < NUM_ATOMS && atom2 < NUM_ATOMS) {
                            COMPUTE_INTERACTION
                            dEdR /= -r;
                        }
                        energy += tempEnergy;
                        delta.xyz *= dEdR;
                        force.xyz -= delta.xyz;
                        atom2 = tbx+tj;
                        local_force[atom2].xyz += delta.xyz;
                        RECORD_DERIVATIVE_2
#ifdef USE_CUTOFF
                        }
#endif
#ifdef USE_EXCLUSIONS
                        excl >>= 1;
#endif
                        tj = (tj + 1) & (TILE_SIZE - 1);
                    }
                }
            }
        }
        lasty = y;
        
        // Write results.  We need to coordinate between warps to make sure no two of them
        // ever try to write to the same piece of memory at the same time.
        
        int writeX = (pos < end ? x : -1);
        int writeY = (pos < end && x != y ? y : -1);
        if (tgx == 0)
            reservedBlocks[localGroupIndex] = (int2)(writeX, writeY);
        bool done = false;
        int doneIndex = 0;
        int checkIndex = 0;
        while (true) {
            // See if any warp still needs to write its data.

            bool allDone = true;
            barrier(CLK_LOCAL_MEM_FENCE);
            while (doneIndex < WARPS_PER_GROUP && allDone) {
                if (reservedBlocks[doneIndex].x != -1)
                    allDone = false;
                else
                    doneIndex++;
            }
            if (allDone)
                break;
            if (!done) {
                // See whether this warp can write its data.  This requires that no previous warp
                // is trying to write to the same block of the buffer.

                bool canWrite = (writeX != -1);
                while (checkIndex < localGroupIndex && canWrite) {
                    if ((reservedBlocks[checkIndex].x == x || reservedBlocks[checkIndex].y == x) ||
                            (writeY != -1 && (reservedBlocks[checkIndex].x == y || reservedBlocks[checkIndex].y == y)))
                        canWrite = false;
                    else
                        checkIndex++;
                }
                if (canWrite) {
                    // Write the data to global memory, then mark this warp as done.

                    if (writeX > -1) {
                        const unsigned int offset = x*TILE_SIZE + tgx + get_group_id(0)*PADDED_NUM_ATOMS;
                        forceBuffers[offset].xyz += force.xyz;
                        STORE_DERIVATIVES_1
                    }
                    if (writeY > -1) {
                        const unsigned int offset = y*TILE_SIZE + tgx + get_group_id(0)*PADDED_NUM_ATOMS;
                        forceBuffers[offset].xyz += local_force[get_local_id(0)].xyz;
                        STORE_DERIVATIVES_2
                    }
                    done = true;
                    if (tgx == 0)
                        reservedBlocks[localGroupIndex] = (int2)(-1, -1);
                }
            }
        }
        pos++;
    } while (pos < end);
    energyBuffer[get_global_id(0)] += energy;
}
