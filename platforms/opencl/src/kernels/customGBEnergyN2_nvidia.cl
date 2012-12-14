#pragma OPENCL EXTENSION cl_khr_global_int32_base_atomics : enable
#ifdef SUPPORTS_64_BIT_ATOMICS
#pragma OPENCL EXTENSION cl_khr_int64_base_atomics : enable
#define STORE_DERIVATIVE_1(INDEX) atom_add(&derivBuffers[offset+(INDEX-1)*PADDED_NUM_ATOMS], (long) (deriv##INDEX##_1*0x100000000));
#define STORE_DERIVATIVE_2(INDEX) atom_add(&derivBuffers[offset+(INDEX-1)*PADDED_NUM_ATOMS], (long) (local_deriv##INDEX[get_local_id(0)]*0x100000000));
#else
#define STORE_DERIVATIVE_1(INDEX) derivBuffers##INDEX[offset] += deriv##INDEX##_1;
#define STORE_DERIVATIVE_2(INDEX) derivBuffers##INDEX[offset] += local_deriv##INDEX[get_local_id(0)];
#endif
#define TILE_SIZE 32

/**
 * Compute a force based on pair interactions.
 */
__kernel void computeN2Energy(
#ifdef SUPPORTS_64_BIT_ATOMICS
        __global long* restrict forceBuffers,
#else
        __global real4* restrict forceBuffers,
#endif
        __global real* restrict energyBuffer, __local real4* restrict local_force,
	__global const real4* restrict posq, __local real4* restrict local_posq, __global const unsigned int* restrict exclusions, __global const unsigned int* restrict exclusionIndices,
        __global const unsigned int* restrict exclusionRowIndices, __local real4* restrict tempBuffer,
#ifdef USE_CUTOFF
        __global const ushort2* restrict tiles, __global const unsigned int* restrict interactionCount, real4 periodicBoxSize, real4 invPeriodicBoxSize, unsigned int maxTiles, __global const unsigned int* restrict interactionFlags
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
    real energy = 0;
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
        real4 force = 0;
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
            real4 posq1 = posq[atom1];
            LOAD_ATOM1_PARAMETERS

            // Locate the exclusion data for this tile.

#ifdef USE_EXCLUSIONS
            if (tgx < 2)
                exclusionRange[2*localGroupIndex+tgx] = exclusionRowIndices[x+tgx];
            if (tgx == 0)
                exclusionIndex[localGroupIndex] = -1;
            for (unsigned int i = exclusionRange[2*localGroupIndex]+tgx; i < exclusionRange[2*localGroupIndex+1]; i += TILE_SIZE)
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
                    real4 posq2 = local_posq[atom2];
                    real4 delta = (real4) (posq2.xyz - posq1.xyz, 0);
#ifdef USE_PERIODIC
                    delta.x -= floor(delta.x*invPeriodicBoxSize.x+0.5f)*periodicBoxSize.x;
                    delta.y -= floor(delta.y*invPeriodicBoxSize.y+0.5f)*periodicBoxSize.y;
                    delta.z -= floor(delta.z*invPeriodicBoxSize.z+0.5f)*periodicBoxSize.z;
#endif
                    real r2 = delta.x*delta.x + delta.y*delta.y + delta.z*delta.z;
#ifdef USE_CUTOFF
                    if (r2 < CUTOFF_SQUARED) {
#endif
                    real invR = RSQRT(r2);
                    real r = RECIP(invR);
                    LOAD_ATOM2_PARAMETERS
                    atom2 = y*TILE_SIZE+j;
                    real dEdR = 0;
                    real tempEnergy = 0;
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
                local_force[localAtomIndex] = 0;
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
                        real4 posq2 = local_posq[atom2];
                        real4 delta = (real4) (posq2.xyz - posq1.xyz, 0);
#ifdef USE_PERIODIC
                        delta.x -= floor(delta.x*invPeriodicBoxSize.x+0.5f)*periodicBoxSize.x;
                        delta.y -= floor(delta.y*invPeriodicBoxSize.y+0.5f)*periodicBoxSize.y;
                        delta.z -= floor(delta.z*invPeriodicBoxSize.z+0.5f)*periodicBoxSize.z;
#endif
                        real r2 = delta.x*delta.x + delta.y*delta.y + delta.z*delta.z;
#ifdef USE_CUTOFF
                        if (r2 < CUTOFF_SQUARED) {
#endif
                        real invR = RSQRT(r2);
                        real r = RECIP(invR);
                        LOAD_ATOM2_PARAMETERS
                        atom2 = y*TILE_SIZE+tj;
                        real dEdR = 0;
                        real tempEnergy = 0;
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
        
#ifdef SUPPORTS_64_BIT_ATOMICS
        if (pos < end) {
            const unsigned int offset = x*TILE_SIZE + tgx;
            atom_add(&forceBuffers[offset], (long) (force.x*0x100000000));
            atom_add(&forceBuffers[offset+PADDED_NUM_ATOMS], (long) (force.y*0x100000000));
            atom_add(&forceBuffers[offset+2*PADDED_NUM_ATOMS], (long) (force.z*0x100000000));
            STORE_DERIVATIVES_1
        }
        if (pos < end && x != y) {
            const unsigned int offset = y*TILE_SIZE + tgx;
            atom_add(&forceBuffers[offset], (long) (local_force[get_local_id(0)].x*0x100000000));
            atom_add(&forceBuffers[offset+PADDED_NUM_ATOMS], (long) (local_force[get_local_id(0)].y*0x100000000));
            atom_add(&forceBuffers[offset+2*PADDED_NUM_ATOMS], (long) (local_force[get_local_id(0)].z*0x100000000));
            STORE_DERIVATIVES_2
        }
#else
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
#endif
        pos++;
    } while (pos < end);
    energyBuffer[get_global_id(0)] += energy;
}
