typedef struct {
    real x, y, z;
    real q;
    real fx, fy, fz;
    ATOM_PARAMETER_DATA
#ifndef PARAMETER_SIZE_IS_EVEN
    real padding;
#endif
} AtomData;

extern "C" __global__ void computeInteractionGroups(
        unsigned long long* __restrict__ forceBuffers, real* __restrict__ energyBuffer, const real4* __restrict__ posq, const int4* __restrict__ groupData,
        real4 periodicBoxSize, real4 invPeriodicBoxSize
        PARAMETER_ARGUMENTS) {
    const unsigned int totalWarps = (blockDim.x*gridDim.x)/TILE_SIZE;
    const unsigned int warp = (blockIdx.x*blockDim.x+threadIdx.x)/TILE_SIZE; // global warpIndex
    const unsigned int tgx = threadIdx.x & (TILE_SIZE-1); // index within the warp
    const unsigned int tbx = threadIdx.x - tgx;           // block warpIndex
    real energy = 0.0f;
    __shared__ AtomData localData[LOCAL_MEMORY_SIZE];

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
        for (int j = rangeStart; j < rangeEnd; j++) {
            bool isExcluded = (((exclusions>>tj)&1) == 0);
            int localIndex = tbx+tj;
            posq2 = make_real4(localData[localIndex].x, localData[localIndex].y, localData[localIndex].z, localData[localIndex].q);
            real3 delta = make_real3(posq2.x-posq1.x, posq2.y-posq1.y, posq2.z-posq1.z);
#ifdef USE_PERIODIC
            delta.x -= floor(delta.x*invPeriodicBoxSize.x+0.5f)*periodicBoxSize.x;
            delta.y -= floor(delta.y*invPeriodicBoxSize.y+0.5f)*periodicBoxSize.y;
            delta.z -= floor(delta.z*invPeriodicBoxSize.z+0.5f)*periodicBoxSize.z;
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
        if (exclusions != 0) {
            atomicAdd(&forceBuffers[atom1], static_cast<unsigned long long>((long long) (force.x*0x100000000)));
            atomicAdd(&forceBuffers[atom1+PADDED_NUM_ATOMS], static_cast<unsigned long long>((long long) (force.y*0x100000000)));
            atomicAdd(&forceBuffers[atom1+2*PADDED_NUM_ATOMS], static_cast<unsigned long long>((long long) (force.z*0x100000000)));
        }
        atomicAdd(&forceBuffers[atom2], static_cast<unsigned long long>((long long) (localData[threadIdx.x].fx*0x100000000)));
        atomicAdd(&forceBuffers[atom2+PADDED_NUM_ATOMS], static_cast<unsigned long long>((long long) (localData[threadIdx.x].fy*0x100000000)));
        atomicAdd(&forceBuffers[atom2+2*PADDED_NUM_ATOMS], static_cast<unsigned long long>((long long) (localData[threadIdx.x].fz*0x100000000)));
    }
    energyBuffer[blockIdx.x*blockDim.x+threadIdx.x] += energy;
}
