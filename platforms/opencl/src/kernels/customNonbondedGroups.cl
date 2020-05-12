#ifdef SUPPORTS_64_BIT_ATOMICS
#pragma OPENCL EXTENSION cl_khr_int64_base_atomics : enable
#endif

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
int reduceMax(int val, __local int* temp) {
    int indexInWarp = get_local_id(0)%32;
    temp[get_local_id(0)] = val;
    SYNC_WARPS;
    for (int offset = 16; offset > 0; offset /= 2) {
        if (indexInWarp < offset)
            temp[get_local_id(0)] = max(temp[get_local_id(0)], temp[get_local_id(0)+offset]);
        SYNC_WARPS;
    }
    return temp[get_local_id(0)-indexInWarp];
}

/**
 * This function is used on devices that don't support 64 bit atomics.  Multiple threads within
 * a single tile might have computed forces on the same atom.  This loops over them and makes sure
 * that only one thread updates the force on any given atom.
 */
void writeForces(__global real4* forceBuffers,__local AtomData* localData, int atomIndex) {
    localData[get_local_id(0)].x = atomIndex;
    SYNC_WARPS;
    real4 forceSum = (real4) 0;
    int start = (get_local_id(0)/TILE_SIZE)*TILE_SIZE;
    int end = start+32;
    bool isFirst = true;
    for (int i = start; i < end; i++)
        if (localData[i].x == atomIndex) {
            forceSum += (real4) (localData[i].fx, localData[i].fy, localData[i].fz, 0);
            isFirst &= (i >= get_local_id(0));
        }
    const unsigned int warp = get_global_id(0)/TILE_SIZE;
    unsigned int offset = atomIndex + warp*PADDED_NUM_ATOMS;
    if (isFirst)
        forceBuffers[offset] += forceSum;
    SYNC_WARPS;
}

__kernel void computeInteractionGroups(
#ifdef SUPPORTS_64_BIT_ATOMICS
        __global long* restrict forceBuffers,
#else
        __global real4* restrict forceBuffers,
#endif
        __global mixed* restrict energyBuffer, __global const real4* restrict posq, __global const int4* restrict groupData,
        __global int* restrict numGroupTiles, int useNeighborList,
        real4 periodicBoxSize, real4 invPeriodicBoxSize, real4 periodicBoxVecX, real4 periodicBoxVecY, real4 periodicBoxVecZ
        PARAMETER_ARGUMENTS) {
    const unsigned int totalWarps = get_global_size(0)/TILE_SIZE;
    const unsigned int warp = get_global_id(0)/TILE_SIZE; // global warpIndex
    const unsigned int tgx = get_local_id(0) & (TILE_SIZE-1); // index within the warp
    const unsigned int tbx = get_local_id(0) - tgx;           // block warpIndex
    mixed energy = 0;
    INIT_DERIVATIVES
    __local AtomData localData[LOCAL_MEMORY_SIZE];
    __local int reductionBuffer[LOCAL_MEMORY_SIZE];

    const unsigned int startTile = (useNeighborList ? warp*numGroupTiles[0]/totalWarps : FIRST_TILE+warp*(LAST_TILE-FIRST_TILE)/totalWarps);
    const unsigned int endTile = (useNeighborList ? (warp+1)*numGroupTiles[0]/totalWarps : FIRST_TILE+(warp+1)*(LAST_TILE-FIRST_TILE)/totalWarps);
    for (int tile = startTile; tile < endTile; tile++) {
        const int4 atomData = groupData[TILE_SIZE*tile+tgx];
        const int atom1 = atomData.x;
        const int atom2 = atomData.y;
        const int rangeStart = atomData.z&0xFFFF;
        const int rangeEnd = (atomData.z>>16)&0xFFFF;
        const int exclusions = atomData.w;
        real4 posq1 = posq[atom1];
        LOAD_ATOM1_PARAMETERS
        real4 force = (real4) (0);
        real4 posq2 = posq[atom2];
        localData[get_local_id(0)].x = posq2.x;
        localData[get_local_id(0)].y = posq2.y;
        localData[get_local_id(0)].z = posq2.z;
        localData[get_local_id(0)].q = posq2.w;
        LOAD_LOCAL_PARAMETERS
        localData[get_local_id(0)].fx = 0.0f;
        localData[get_local_id(0)].fy = 0.0f;
        localData[get_local_id(0)].fz = 0.0f;
        int tj = tgx;
        int rangeStop = rangeStart + reduceMax(rangeEnd-rangeStart, reductionBuffer);
        SYNC_WARPS;
        for (int j = rangeStart; j < rangeStop; j++) {
            if (j < rangeEnd) {
                bool isExcluded = (((exclusions>>tj)&1) == 0);
                int localIndex = tbx+tj;
                posq2 = (real4) (localData[localIndex].x, localData[localIndex].y, localData[localIndex].z, localData[localIndex].q);
                real4 delta = (real4) (posq2.xyz - posq1.xyz, 0);
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
                    force.xyz -= delta.xyz;
                    localData[localIndex].fx += delta.x;
                    localData[localIndex].fy += delta.y;
                    localData[localIndex].fz += delta.z;
#ifdef USE_CUTOFF
                }
#endif
            }
            tj = (tj == rangeEnd-1 ? rangeStart : tj+1);
            SYNC_WARPS;
        }
#ifdef SUPPORTS_64_BIT_ATOMICS
        if (exclusions != 0) {
            atom_add(&forceBuffers[atom1], (long) (force.x*0x100000000));
            atom_add(&forceBuffers[atom1+PADDED_NUM_ATOMS], (long) (force.y*0x100000000));
            atom_add(&forceBuffers[atom1+2*PADDED_NUM_ATOMS], (long) (force.z*0x100000000));
        }
        atom_add(&forceBuffers[atom2], (long) (localData[get_local_id(0)].fx*0x100000000));
        atom_add(&forceBuffers[atom2+PADDED_NUM_ATOMS], (long) (localData[get_local_id(0)].fy*0x100000000));
        atom_add(&forceBuffers[atom2+2*PADDED_NUM_ATOMS], (long) (localData[get_local_id(0)].fz*0x100000000));
#else
        writeForces(forceBuffers, localData, atom2);
        localData[get_local_id(0)].fx = force.x;
        localData[get_local_id(0)].fy = force.y;
        localData[get_local_id(0)].fz = force.z;
        writeForces(forceBuffers, localData, atom1);
#endif
    }
    energyBuffer[get_global_id(0)] += energy;
    SAVE_DERIVATIVES
}

/**
 * If the neighbor list needs to be rebuilt, reset the number of tiles to 0.  This is
 * executed by a single thread.
 */
__kernel void prepareToBuildNeighborList(__global int* restrict rebuildNeighborList, __global int* restrict numGroupTiles) {
    if (rebuildNeighborList[0] == 1)
        numGroupTiles[0] = 0;
}

/**
 * Filter the list of tiles to include only ones that have interactions within the
 * padded cutoff.
 */
__kernel void buildNeighborList(__global int* restrict rebuildNeighborList, __global int* restrict numGroupTiles,
        __global const real4* restrict posq, __global const int4* restrict groupData, __global int4* restrict filteredGroupData,
        real4 periodicBoxSize, real4 invPeriodicBoxSize, real4 periodicBoxVecX, real4 periodicBoxVecY, real4 periodicBoxVecZ) {
    
    // If the neighbor list doesn't need to be rebuilt on this step, return immediately.
    
    if (rebuildNeighborList[0] == 0)
        return;

    const unsigned int totalWarps = get_global_size(0)/TILE_SIZE;
    const unsigned int warp = get_global_id(0)/TILE_SIZE; // global warpIndex
    const unsigned int local_warp = get_local_id(0)/TILE_SIZE; // local warpIndex
    const unsigned int tgx = get_local_id(0) & (TILE_SIZE-1); // index within the warp
    const unsigned int tbx = get_local_id(0) - tgx;           // block warpIndex
    __local real4 localPos[LOCAL_MEMORY_SIZE];
    __local volatile bool anyInteraction[WARPS_IN_BLOCK];
    __local volatile int tileIndex[WARPS_IN_BLOCK];
    __local int reductionBuffer[LOCAL_MEMORY_SIZE];

    const unsigned int startTile = warp*NUM_TILES/totalWarps;
    const unsigned int endTile = (warp+1)*NUM_TILES/totalWarps;
    for (int tile = startTile; tile < endTile; tile++) {
        const int4 atomData = groupData[TILE_SIZE*tile+tgx];
        const int atom1 = atomData.x;
        const int atom2 = atomData.y;
        const int rangeStart = atomData.z&0xFFFF;
        const int rangeEnd = (atomData.z>>16)&0xFFFF;
        const int exclusions = atomData.w;
        real4 posq1 = posq[atom1];
        localPos[get_local_id(0)] = posq[atom2];
        if (tgx == 0)
            anyInteraction[local_warp] = false;
        int tj = tgx;
        int rangeStop = rangeStart + reduceMax(rangeEnd-rangeStart, reductionBuffer);
        SYNC_WARPS;
        for (int j = rangeStart; j < rangeStop && !anyInteraction[local_warp]; j++) {
            SYNC_WARPS;
            if (j < rangeEnd) {
                bool isExcluded = (((exclusions>>tj)&1) == 0);
                int localIndex = tbx+tj;
                real4 delta = (real4) (localPos[localIndex].xyz - posq1.xyz, 0);
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
                tileIndex[local_warp] = atomic_add(numGroupTiles, 1);
            SYNC_WARPS;
            filteredGroupData[TILE_SIZE*tileIndex[local_warp]+tgx] = atomData;
        }
    }
}
