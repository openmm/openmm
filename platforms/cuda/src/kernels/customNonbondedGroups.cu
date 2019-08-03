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
 * every thread.  This is only needed on Volta and later.  On earlier architectures,
 * we can just return the value that was passed in.
 */
__device__ int reduceMax(int val) {
#if __CUDA_ARCH__ >= 700
    for (int mask = 16; mask > 0; mask /= 2) 
        val = max(val, __shfl_xor_sync(0xffffffff, val, mask));
#endif
    return val;
}

extern "C" __global__ void computeInteractionGroups(
        unsigned long long* __restrict__ forceBuffers, mixed* __restrict__ energyBuffer, const real4* __restrict__ posq, const int4* __restrict__ groupData,
        const int* __restrict__ numGroupTiles, bool useNeighborList,
        real4 periodicBoxSize, real4 invPeriodicBoxSize, real4 periodicBoxVecX, real4 periodicBoxVecY, real4 periodicBoxVecZ
        PARAMETER_ARGUMENTS) {
    const unsigned int totalWarps = (blockDim.x*gridDim.x)/TILE_SIZE;
    const unsigned int warp = (blockIdx.x*blockDim.x+threadIdx.x)/TILE_SIZE; // global warpIndex
    const unsigned int tgx = threadIdx.x & (TILE_SIZE-1); // index within the warp
    const unsigned int tbx = threadIdx.x - tgx;           // block warpIndex
    mixed energy = 0;
    INIT_DERIVATIVES
    __shared__ AtomData localData[LOCAL_MEMORY_SIZE];

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
        real3 force = make_real3(0);
        real4 posq2 = posq[atom2];
        localData[threadIdx.x].x = posq2.x;
        localData[threadIdx.x].y = posq2.y;
        localData[threadIdx.x].z = posq2.z;
        localData[threadIdx.x].q = posq2.w;
        LOAD_LOCAL_PARAMETERS
        localData[threadIdx.x].fx = 0.0f;
        localData[threadIdx.x].fy = 0.0f;
        localData[threadIdx.x].fz = 0.0f;
        int tj = tgx;
        int rangeStop = rangeStart + reduceMax(rangeEnd-rangeStart);
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
            atomicAdd(&forceBuffers[atom1], static_cast<unsigned long long>((long long) (force.x*0x100000000)));
            atomicAdd(&forceBuffers[atom1+PADDED_NUM_ATOMS], static_cast<unsigned long long>((long long) (force.y*0x100000000)));
            atomicAdd(&forceBuffers[atom1+2*PADDED_NUM_ATOMS], static_cast<unsigned long long>((long long) (force.z*0x100000000)));
        }
        atomicAdd(&forceBuffers[atom2], static_cast<unsigned long long>((long long) (localData[threadIdx.x].fx*0x100000000)));
        atomicAdd(&forceBuffers[atom2+PADDED_NUM_ATOMS], static_cast<unsigned long long>((long long) (localData[threadIdx.x].fy*0x100000000)));
        atomicAdd(&forceBuffers[atom2+2*PADDED_NUM_ATOMS], static_cast<unsigned long long>((long long) (localData[threadIdx.x].fz*0x100000000)));
        SYNC_WARPS;
    }
    energyBuffer[blockIdx.x*blockDim.x+threadIdx.x] += energy;
    SAVE_DERIVATIVES
}

/**
 * If the neighbor list needs to be rebuilt, reset the number of tiles to 0.  This is
 * executed by a single thread.
 */
extern "C" __global__  void prepareToBuildNeighborList(int* __restrict__ rebuildNeighborList, int* __restrict__ numGroupTiles) {
    if (rebuildNeighborList[0] == 1)
        numGroupTiles[0] = 0;
}

/**
 * Filter the list of tiles to include only ones that have interactions within the
 * padded cutoff.
 */
extern "C" __global__  void buildNeighborList(int* __restrict__ rebuildNeighborList, int* __restrict__ numGroupTiles,
        const real4* __restrict__ posq, const int4* __restrict__ groupData, int4* __restrict__ filteredGroupData,
        real4 periodicBoxSize, real4 invPeriodicBoxSize, real4 periodicBoxVecX, real4 periodicBoxVecY, real4 periodicBoxVecZ) {
    
    // If the neighbor list doesn't need to be rebuilt on this step, return immediately.
    
    if (rebuildNeighborList[0] == 0)
        return;

    const unsigned int totalWarps = (blockDim.x*gridDim.x)/TILE_SIZE;
    const unsigned int warp = (blockIdx.x*blockDim.x+threadIdx.x)/TILE_SIZE; // global warpIndex
    const unsigned int local_warp = threadIdx.x/TILE_SIZE; // local warpIndex
    const unsigned int tgx = threadIdx.x & (TILE_SIZE-1); // index within the warp
    const unsigned int tbx = threadIdx.x - tgx;           // block warpIndex

    __shared__ real4 localPos[LOCAL_MEMORY_SIZE];
    __shared__ volatile bool anyInteraction[WARPS_IN_BLOCK];
    __shared__ volatile int tileIndex[WARPS_IN_BLOCK];

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
        localPos[threadIdx.x] = posq[atom2];
        if (tgx == 0)
            anyInteraction[local_warp] = false;
        int tj = tgx;
        int rangeStop = rangeStart + reduceMax(rangeEnd-rangeStart);
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
                tileIndex[local_warp] = atomicAdd(numGroupTiles, 1);
            SYNC_WARPS;
            filteredGroupData[TILE_SIZE*tileIndex[local_warp]+tgx] = atomData;
        }
    }
}
