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
 * Find the maximum of a value across all threads in a warp, and return that to
 * every thread.
 */
DEVICE int reduceMax(int val, LOCAL_ARG int* temp) {
#if defined(__CUDA_ARCH__) && __CUDA_ARCH__ >= 700
    // CUDA lets us do this slightly more efficiently by using shuffle operations.
    for (int mask = 16; mask > 0; mask /= 2)
        val = max(val, __shfl_xor_sync(0xffffffff, val, mask));
    return val;
#elif defined(USE_HIP)
    for (int mask = 16; mask > 0; mask /= 2)
        val = max(val, __shfl_xor(val, mask, 32));
    return val;
#else
    int indexInWarp = LOCAL_ID%32;
    temp[LOCAL_ID] = val;
    SYNC_WARPS;
    for (int offset = 16; offset > 0; offset /= 2) {
        if (indexInWarp < offset)
            temp[LOCAL_ID] = max(temp[LOCAL_ID], temp[LOCAL_ID+offset]);
        SYNC_WARPS;
    }
    return temp[LOCAL_ID-indexInWarp];
#endif
}

KERNEL void computeInteractionGroups(
        GLOBAL mm_ulong* RESTRICT forceBuffers,
        GLOBAL mixed* RESTRICT energyBuffer, GLOBAL const real4* RESTRICT posq, GLOBAL const int4* RESTRICT groupData,
        GLOBAL const int* RESTRICT numGroupTiles, int useNeighborList,
        real4 periodicBoxSize, real4 invPeriodicBoxSize, real4 periodicBoxVecX, real4 periodicBoxVecY, real4 periodicBoxVecZ, int numDerivatives
        PARAMETER_ARGUMENTS) {
    const unsigned int totalWarps = GLOBAL_SIZE/TILE_SIZE;
    const unsigned int warp = GLOBAL_ID/TILE_SIZE; // global warpIndex
    const unsigned int tgx = LOCAL_ID & (TILE_SIZE-1); // index within the warp
    const unsigned int tbx = LOCAL_ID - tgx;           // block warpIndex
    mixed energy = 0;
    INIT_DERIVATIVES
    LOCAL AtomData localData[LOCAL_MEMORY_SIZE];
    LOCAL int reductionBuffer[LOCAL_MEMORY_SIZE];

    mm_ulong wl = warp;
    const unsigned int startTile = (unsigned int) (useNeighborList ? wl*numGroupTiles[0]/totalWarps : FIRST_TILE+wl*(LAST_TILE-FIRST_TILE)/totalWarps);
    const unsigned int endTile = (unsigned int) (useNeighborList ? (wl+1)*numGroupTiles[0]/totalWarps : FIRST_TILE+(wl+1)*(LAST_TILE-FIRST_TILE)/totalWarps);
    for (int tile = startTile; tile < endTile; tile++) {
        const int4 atomData = groupData[TILE_SIZE*tile+tgx];
        const int atom1 = atomData.x;
        const int atom2 = atomData.y;
        const int rangeStart = atomData.z&0xFFFF;
        const int rangeEnd = (atomData.z>>16)&0xFFFF;
        const int exclusions = atomData.w;
        real4 posq1 = posq[atom1];
        LOAD_ATOM1_PARAMETERS
        real3 force = make_real3(0);
        real4 posq2 = posq[atom2];
        localData[LOCAL_ID].x = posq2.x;
        localData[LOCAL_ID].y = posq2.y;
        localData[LOCAL_ID].z = posq2.z;
        localData[LOCAL_ID].q = posq2.w;
        LOAD_LOCAL_PARAMETERS
        localData[LOCAL_ID].fx = 0.0f;
        localData[LOCAL_ID].fy = 0.0f;
        localData[LOCAL_ID].fz = 0.0f;
        int tj = tgx;
        int rangeStop = rangeStart + reduceMax(rangeEnd-rangeStart, reductionBuffer);
        SYNC_WARPS;
        for (int j = rangeStart; j < rangeStop; j++) {
            if (j < rangeEnd) {
                bool isExcluded = (((exclusions>>tj)&1) == 0);
                int localIndex = tbx+tj;
                posq2 = make_real4(localData[localIndex].x, localData[localIndex].y, localData[localIndex].z, localData[localIndex].q);
                real3 delta = make_real3(posq2.x-posq1.x, posq2.y-posq1.y, posq2.z-posq1.z);
#ifdef USE_PERIODIC
                APPLY_PERIODIC_TO_DELTA(delta)
#endif
                real r2 = delta.x*delta.x + delta.y*delta.y + delta.z*delta.z;
#ifdef USE_CUTOFF
                if (!isExcluded && r2 < CUTOFF_SQUARED) {
#endif
                    real invR = RSQRT(r2);
                    real r = r2*invR;
                    LOAD_ATOM2_PARAMETERS
                    real dEdR = 0.0f;
                    real tempEnergy = 0.0f;
                    const real interactionScale = 1.0f;
                    COMPUTE_INTERACTION
                    energy += tempEnergy;
                    delta *= dEdR;
                    force.x -= delta.x;
                    force.y -= delta.y;
                    force.z -= delta.z;
                    localData[localIndex].fx += delta.x;
                    localData[localIndex].fy += delta.y;
                    localData[localIndex].fz += delta.z;
#ifdef USE_CUTOFF
                }
#endif
                tj = (tj == rangeEnd-1 ? rangeStart : tj+1);
            }
            SYNC_WARPS;
        }
        if (exclusions != 0) {
            ATOMIC_ADD(&forceBuffers[atom1], (mm_ulong) realToFixedPoint(force.x));
            ATOMIC_ADD(&forceBuffers[atom1+PADDED_NUM_ATOMS], (mm_ulong) realToFixedPoint(force.y));
            ATOMIC_ADD(&forceBuffers[atom1+2*PADDED_NUM_ATOMS], (mm_ulong) realToFixedPoint(force.z));
        }
        ATOMIC_ADD(&forceBuffers[atom2], (mm_ulong) realToFixedPoint(localData[LOCAL_ID].fx));
        ATOMIC_ADD(&forceBuffers[atom2+PADDED_NUM_ATOMS], (mm_ulong) realToFixedPoint(localData[LOCAL_ID].fy));
        ATOMIC_ADD(&forceBuffers[atom2+2*PADDED_NUM_ATOMS], (mm_ulong) realToFixedPoint(localData[LOCAL_ID].fz));
        SYNC_WARPS;
    }
    energyBuffer[GLOBAL_ID] += energy;
    SAVE_DERIVATIVES
}

/**
 * If the neighbor list needs to be rebuilt, reset the number of tiles to 0.  This is
 * executed by a single thread.
 */
KERNEL void prepareToBuildNeighborList(GLOBAL int* RESTRICT rebuildNeighborList, GLOBAL int* RESTRICT numGroupTiles) {
    if (rebuildNeighborList[0] == 1)
        numGroupTiles[0] = 0;
}

/**
 * Filter the list of tiles to include only ones that have interactions within the
 * padded cutoff.
 */
KERNEL void buildNeighborList(GLOBAL int* RESTRICT rebuildNeighborList, GLOBAL int* RESTRICT numGroupTiles,
        GLOBAL const real4* RESTRICT posq, GLOBAL const int4* RESTRICT groupData, GLOBAL int4* RESTRICT filteredGroupData,
        real4 periodicBoxSize, real4 invPeriodicBoxSize, real4 periodicBoxVecX, real4 periodicBoxVecY, real4 periodicBoxVecZ) {
    
    // If the neighbor list doesn't need to be rebuilt on this step, return immediately.
    
    if (rebuildNeighborList[0] == 0)
        return;

    const unsigned int totalWarps = GLOBAL_SIZE/TILE_SIZE;
    const unsigned int warp = GLOBAL_ID/TILE_SIZE; // global warpIndex
    const unsigned int local_warp = LOCAL_ID/TILE_SIZE; // local warpIndex
    const unsigned int tgx = LOCAL_ID & (TILE_SIZE-1); // index within the warp
    const unsigned int tbx = LOCAL_ID - tgx;           // block warpIndex
    LOCAL real4 localPos[LOCAL_MEMORY_SIZE];
    LOCAL volatile bool anyInteraction[WARPS_IN_BLOCK];
    LOCAL volatile int tileIndex[WARPS_IN_BLOCK];
    LOCAL int reductionBuffer[LOCAL_MEMORY_SIZE];

    const unsigned int startTile = (unsigned int) (warp*(mm_ulong)NUM_TILES/totalWarps);
    const unsigned int endTile = (unsigned int) ((warp+1)*(mm_ulong)NUM_TILES/totalWarps);
    for (int tile = startTile; tile < endTile; tile++) {
        const int4 atomData = groupData[TILE_SIZE*tile+tgx];
        const int atom1 = atomData.x;
        const int atom2 = atomData.y;
        const int rangeStart = atomData.z&0xFFFF;
        const int rangeEnd = (atomData.z>>16)&0xFFFF;
        const int exclusions = atomData.w;
        real4 posq1 = posq[atom1];
        localPos[LOCAL_ID] = posq[atom2];
        if (tgx == 0)
            anyInteraction[local_warp] = false;
        int tj = tgx;
        int rangeStop = rangeStart + reduceMax(rangeEnd-rangeStart, reductionBuffer);
        SYNC_WARPS;
        for (int j = rangeStart; j < rangeStop && !anyInteraction[local_warp]; j++) {
            SYNC_WARPS;
            if (j < rangeEnd && tj < rangeEnd) {
                bool isExcluded = (((exclusions>>tj)&1) == 0);
                int localIndex = tbx+tj;
                real3 delta = make_real3(localPos[localIndex].x-posq1.x, localPos[localIndex].y-posq1.y, localPos[localIndex].z-posq1.z);
#ifdef USE_PERIODIC
                APPLY_PERIODIC_TO_DELTA(delta)
#endif
                real r2 = delta.x*delta.x + delta.y*delta.y + delta.z*delta.z;
                if (!isExcluded && r2 < PADDED_CUTOFF_SQUARED)
                    anyInteraction[local_warp] = true;
            }
            tj = (tj == rangeEnd-1 ? rangeStart : tj+1);
            SYNC_WARPS;
        }
        if (anyInteraction[local_warp]) {
            SYNC_WARPS;
            if (tgx == 0)
                tileIndex[local_warp] = ATOMIC_ADD(numGroupTiles, 1);
            SYNC_WARPS;
            filteredGroupData[TILE_SIZE*tileIndex[local_warp]+tgx] = atomData;
        }
    }
}
