#pragma OPENCL EXTENSION cl_khr_global_int32_base_atomics : enable
#ifdef SUPPORTS_64_BIT_ATOMICS
#pragma OPENCL EXTENSION cl_khr_int64_base_atomics : enable
#endif
#define TILE_SIZE 32
#define WARPS_PER_GROUP (FORCE_WORK_GROUP_SIZE/TILE_SIZE)

typedef struct {
    real x, y, z;
    real q;
    real fx, fy, fz;
    ATOM_PARAMETER_DATA
#ifndef PARAMETER_SIZE_IS_EVEN
    real padding;
#endif
} AtomData;

/**
 * Compute nonbonded interactions.
 */
__kernel void computeNonbonded(
#ifdef SUPPORTS_64_BIT_ATOMICS
        __global long* restrict forceBuffers,
#else
        __global real4* restrict forceBuffers,
#endif
        __global real* restrict energyBuffer, __global const real4* restrict posq, __global const unsigned int* restrict exclusions,
        __global const unsigned int* restrict exclusionIndices, __global const unsigned int* restrict exclusionRowIndices,
        unsigned int startTileIndex, unsigned int endTileIndex,
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
    unsigned int pos = (numTiles > maxTiles ? startTileIndex+warp*(endTileIndex-startTileIndex)/totalWarps : warp*numTiles/totalWarps);
    unsigned int end = (numTiles > maxTiles ? startTileIndex+(warp+1)*(endTileIndex-startTileIndex)/totalWarps : (warp+1)*numTiles/totalWarps);
#else
    unsigned int pos = startTileIndex+warp*numTiles/totalWarps;
    unsigned int end = startTileIndex+(warp+1)*numTiles/totalWarps;
#endif
    real energy = 0;
    __local AtomData localData[FORCE_WORK_GROUP_SIZE];
    __local real tempBuffer[3*FORCE_WORK_GROUP_SIZE];
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
                localData[localAtomIndex].x = posq1.x;
                localData[localAtomIndex].y = posq1.y;
                localData[localAtomIndex].z = posq1.z;
                localData[localAtomIndex].q = posq1.w;
                LOAD_LOCAL_PARAMETERS_FROM_1
#ifdef USE_EXCLUSIONS
                unsigned int excl = exclusions[exclusionIndex[localGroupIndex]+tgx];
#endif
                for (unsigned int j = 0; j < TILE_SIZE; j++) {
#ifdef USE_EXCLUSIONS
                    bool isExcluded = !(excl & 0x1);
#endif
                    int atom2 = tbx+j;
                    real4 posq2 = (real4) (localData[atom2].x, localData[atom2].y, localData[atom2].z, localData[atom2].q);
                    real4 delta = (real4) (posq2.xyz - posq1.xyz, 0);
#ifdef USE_PERIODIC
                    delta.x -= floor(delta.x*invPeriodicBoxSize.x+0.5f)*periodicBoxSize.x;
                    delta.y -= floor(delta.y*invPeriodicBoxSize.y+0.5f)*periodicBoxSize.y;
                    delta.z -= floor(delta.z*invPeriodicBoxSize.z+0.5f)*periodicBoxSize.z;
#endif
                    real r2 = delta.x*delta.x + delta.y*delta.y + delta.z*delta.z;
                    real invR = RSQRT(r2);
                    real r = RECIP(invR);
                    LOAD_ATOM2_PARAMETERS
                    atom2 = y*TILE_SIZE+j;
#ifdef USE_SYMMETRIC
                    real dEdR = 0;
#else
                    real4 dEdR1 = (real4) 0;
                    real4 dEdR2 = (real4) 0;
#endif
                    real tempEnergy = 0;
                    COMPUTE_INTERACTION
                    energy += 0.5f*tempEnergy;
#ifdef USE_SYMMETRIC
                    force.xyz -= delta.xyz*dEdR;
#else
                    force.xyz -= dEdR1.xyz;
#endif
#ifdef USE_EXCLUSIONS
                    excl >>= 1;
#endif
                }
            }
            else {
                // This is an off-diagonal tile.

                const unsigned int localAtomIndex = get_local_id(0);
                unsigned int j = y*TILE_SIZE + tgx;
                real4 tempPosq = posq[j];
                localData[localAtomIndex].x = tempPosq.x;
                localData[localAtomIndex].y = tempPosq.y;
                localData[localAtomIndex].z = tempPosq.z;
                localData[localAtomIndex].q = tempPosq.w;
                LOAD_LOCAL_PARAMETERS_FROM_GLOBAL
                localData[localAtomIndex].fx = 0;
                localData[localAtomIndex].fy = 0;
                localData[localAtomIndex].fz = 0;
#ifdef USE_CUTOFF
                unsigned int flags = (numTiles <= maxTiles ? interactionFlags[pos] : 0xFFFFFFFF);
                if (!hasExclusions && flags != 0xFFFFFFFF) {
                    if (flags == 0) {
                        // No interactions in this tile.
                    }
                    else {
                        // Compute only a subset of the interactions in this tile.

                        for (j = 0; j < TILE_SIZE; j++) {
                            if ((flags&(1<<j)) != 0) {
                                bool isExcluded = false;
                                int atom2 = tbx+j;
                                int bufferIndex = 3*get_local_id(0);
#ifdef USE_SYMMETRIC
                                real dEdR = 0;
#else
                                real4 dEdR1 = (real4) 0;
                                real4 dEdR2 = (real4) 0;
#endif
                                real tempEnergy = 0;
                                real4 posq2 = (real4) (localData[atom2].x, localData[atom2].y, localData[atom2].z, localData[atom2].q);
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
                                    COMPUTE_INTERACTION
                                    energy += tempEnergy;
#ifdef USE_CUTOFF
                                }
#endif
#ifdef USE_SYMMETRIC
                                delta.xyz *= dEdR;
                                force.xyz -= delta.xyz;
                                tempBuffer[bufferIndex] = delta.x;
                                tempBuffer[bufferIndex+1] = delta.y;
                                tempBuffer[bufferIndex+2] = delta.z;
#else
                                force.xyz -= dEdR1.xyz;
                                tempBuffer[bufferIndex] = dEdR2.x;
                                tempBuffer[bufferIndex+1] = dEdR2.y;
                                tempBuffer[bufferIndex+2] = dEdR2.z;
#endif

                                // Sum the forces on atom2.

                                if (tgx % 4 == 0) {
                                    tempBuffer[bufferIndex] += tempBuffer[bufferIndex+3]+tempBuffer[bufferIndex+6]+tempBuffer[bufferIndex+9];
                                    tempBuffer[bufferIndex+1] += tempBuffer[bufferIndex+4]+tempBuffer[bufferIndex+7]+tempBuffer[bufferIndex+10];
                                    tempBuffer[bufferIndex+2] += tempBuffer[bufferIndex+5]+tempBuffer[bufferIndex+8]+tempBuffer[bufferIndex+11];
                                }
                                if (tgx == 0) {
                                    localData[tbx+j].fx += tempBuffer[bufferIndex]+tempBuffer[bufferIndex+12]+tempBuffer[bufferIndex+24]+tempBuffer[bufferIndex+36]+tempBuffer[bufferIndex+48]+tempBuffer[bufferIndex+60]+tempBuffer[bufferIndex+72]+tempBuffer[bufferIndex+84];
                                    localData[tbx+j].fy += tempBuffer[bufferIndex+1]+tempBuffer[bufferIndex+13]+tempBuffer[bufferIndex+25]+tempBuffer[bufferIndex+37]+tempBuffer[bufferIndex+49]+tempBuffer[bufferIndex+61]+tempBuffer[bufferIndex+73]+tempBuffer[bufferIndex+85];
                                    localData[tbx+j].fz += tempBuffer[bufferIndex+2]+tempBuffer[bufferIndex+14]+tempBuffer[bufferIndex+26]+tempBuffer[bufferIndex+38]+tempBuffer[bufferIndex+50]+tempBuffer[bufferIndex+62]+tempBuffer[bufferIndex+74]+tempBuffer[bufferIndex+86];
                                }
                            }
                        }
                    }
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
                    for (j = 0; j < TILE_SIZE; j++) {
#ifdef USE_EXCLUSIONS
                        bool isExcluded = !(excl & 0x1);
#endif
                        int atom2 = tbx+tj;
                        real4 posq2 = (real4) (localData[atom2].x, localData[atom2].y, localData[atom2].z, localData[atom2].q);
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
#ifdef USE_SYMMETRIC
                            real dEdR = 0;
#else
                            real4 dEdR1 = (real4) 0;
                            real4 dEdR2 = (real4) 0;
#endif
                            real tempEnergy = 0;
                            COMPUTE_INTERACTION
                            energy += tempEnergy;
#ifdef USE_SYMMETRIC
                            delta.xyz *= dEdR;
                            force.xyz -= delta.xyz;
                            localData[tbx+tj].fx += delta.x;
                            localData[tbx+tj].fy += delta.y;
                            localData[tbx+tj].fz += delta.z;
#else
                            force.xyz -= dEdR1.xyz;
                            localData[tbx+tj].fx += dEdR2.x;
                            localData[tbx+tj].fy += dEdR2.y;
                            localData[tbx+tj].fz += dEdR2.z;
#endif
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
        
        // Write results.  We need to coordinate between warps to make sure no two of them
        // ever try to write to the same piece of memory at the same time.
        
#ifdef SUPPORTS_64_BIT_ATOMICS
        if (pos < end) {
            const unsigned int offset = x*TILE_SIZE + tgx;
            atom_add(&forceBuffers[offset], (long) (force.x*0x100000000));
            atom_add(&forceBuffers[offset+PADDED_NUM_ATOMS], (long) (force.y*0x100000000));
            atom_add(&forceBuffers[offset+2*PADDED_NUM_ATOMS], (long) (force.z*0x100000000));
        }
        if (pos < end && x != y) {
            const unsigned int offset = y*TILE_SIZE + tgx;
            atom_add(&forceBuffers[offset], (long) (localData[get_local_id(0)].fx*0x100000000));
            atom_add(&forceBuffers[offset+PADDED_NUM_ATOMS], (long) (localData[get_local_id(0)].fy*0x100000000));
            atom_add(&forceBuffers[offset+2*PADDED_NUM_ATOMS], (long) (localData[get_local_id(0)].fz*0x100000000));
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
                    }
                    if (writeY > -1) {
                        const unsigned int offset = y*TILE_SIZE + tgx + get_group_id(0)*PADDED_NUM_ATOMS;
                        forceBuffers[offset] += (real4) (localData[get_local_id(0)].fx, localData[get_local_id(0)].fy, localData[get_local_id(0)].fz, 0);
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
