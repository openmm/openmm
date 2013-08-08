#ifdef SUPPORTS_64_BIT_ATOMICS
#pragma OPENCL EXTENSION cl_khr_int64_base_atomics : enable
#endif

#define WARPS_PER_GROUP (THREAD_BLOCK_SIZE/TILE_SIZE)

typedef struct {
    real x, y, z;
    real q;
    real fx, fy, fz;
    ATOM_PARAMETER_DATA
#ifndef PARAMETER_SIZE_IS_EVEN
    real padding;
#endif
} AtomData;

__kernel void computeInteractionGroups(
        __global long* restrict forceBuffers, __global real* restrict energyBuffer, __global const real4* restrict posq,
        __global const int4* restrict groupData, real4 periodicBoxSize, real4 invPeriodicBoxSize
        PARAMETER_ARGUMENTS) {
    const unsigned int totalWarps = get_global_size(0)/TILE_SIZE;
    const unsigned int warp = get_global_id(0)/TILE_SIZE; // global warpIndex
    const unsigned int tgx = get_local_id(0) & (TILE_SIZE-1); // index within the warp
    const unsigned int tbx = get_local_id(0) - tgx;           // block warpIndex
    real energy = 0.0f;
    __local AtomData localData[THREAD_BLOCK_SIZE];

    const unsigned int startTile = FIRST_TILE+warp*(LAST_TILE-FIRST_TILE)/totalWarps;
    const unsigned int endTile = FIRST_TILE+(warp+1)*(LAST_TILE-FIRST_TILE)/totalWarps;
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
        SYNC_WARPS;
        for (int j = rangeStart; j < rangeEnd; j++) {
            bool isExcluded = (((exclusions>>tj)&1) == 0);
            int localIndex = tbx+tj;
            posq2 = (real4) (localData[localIndex].x, localData[localIndex].y, localData[localIndex].z, localData[localIndex].q);
            real4 delta = (real4) (posq2.xyz - posq1.xyz, 0);
#ifdef USE_PERIODIC
            delta.xyz -= floor(delta.xyz*invPeriodicBoxSize.xyz+0.5f)*periodicBoxSize.xyz;
#endif
            real r2 = delta.x*delta.x + delta.y*delta.y + delta.z*delta.z;
#ifdef USE_CUTOFF
            if (!isExcluded && r2 < CUTOFF_SQUARED) {
#endif
                real invR = RSQRT(r2);
                real r = RECIP(invR);
                LOAD_ATOM2_PARAMETERS
                real dEdR = 0.0f;
                real tempEnergy = 0.0f;
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
            tj = (tj == rangeEnd-1 ? rangeStart : tj+1);
            SYNC_WARPS;
        }
        if (exclusions != 0) {
            atom_add(&forceBuffers[atom1], (long) (force.x*0x100000000));
            atom_add(&forceBuffers[atom1+PADDED_NUM_ATOMS], (long) (force.y*0x100000000));
            atom_add(&forceBuffers[atom1+2*PADDED_NUM_ATOMS], (long) (force.z*0x100000000));
            atom_add(&forceBuffers[atom2], (long) (localData[get_local_id(0)].fx*0x100000000));
            atom_add(&forceBuffers[atom2+PADDED_NUM_ATOMS], (long) (localData[get_local_id(0)].fy*0x100000000));
            atom_add(&forceBuffers[atom2+2*PADDED_NUM_ATOMS], (long) (localData[get_local_id(0)].fz*0x100000000));
        }
    }
    energyBuffer[get_global_id(0)] += energy;
}