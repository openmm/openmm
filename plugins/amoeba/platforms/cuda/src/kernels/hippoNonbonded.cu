// This is a modified version of the standard nonbonded kernel for computing HippoNonbondedForce.
// This is needed because of two ways in which it differs from most nonbonded interactions:
// the force between two atoms doesn't always point along the line between them, and we need
// to accumulate torques as well as forces.

#define WARPS_PER_GROUP (THREAD_BLOCK_SIZE/TILE_SIZE)

#ifndef ENABLE_SHUFFLE
typedef struct {
    real x, y, z;
    real q;
    real fx, fy, fz;
    real tx, ty, tz;
    ATOM_PARAMETER_DATA
#ifndef PARAMETER_SIZE_IS_EVEN
    real padding;
#endif
} AtomData;
#endif

#ifdef ENABLE_SHUFFLE
//support for 64 bit shuffles
static __inline__ __device__ float real_shfl(float var, int srcLane) {
    return SHFL(var, srcLane);
}

static __inline__ __device__ double real_shfl(double var, int srcLane) {
    int hi, lo;
    asm volatile("mov.b64 { %0, %1 }, %2;" : "=r"(lo), "=r"(hi) : "d"(var));
    hi = SHFL(hi, srcLane);
    lo = SHFL(lo, srcLane);
    return __hiloint2double( hi, lo );
}

static __inline__ __device__ long long real_shfl(long long var, int srcLane) {
    int hi, lo;
    asm volatile("mov.b64 { %0, %1 }, %2;" : "=r"(lo), "=r"(hi) : "l"(var));
    hi = SHFL(hi, srcLane);
    lo = SHFL(lo, srcLane);
    // unforunately there isn't an __nv_hiloint2long(hi,lo) intrinsic cast
    int2 fuse; fuse.x = lo; fuse.y = hi;
    return *reinterpret_cast<long long*>(&fuse);
}
#endif

extern "C" __global__ void computeNonbonded(
        unsigned long long* __restrict__ forceBuffers, mixed* __restrict__ energyBuffer, const real4* __restrict__ posq, const tileflags* __restrict__ exclusions,
        const ushort2* __restrict__ exclusionTiles, unsigned int startTileIndex, unsigned int numTileIndices
#ifdef USE_CUTOFF
        , const int* __restrict__ tiles, const unsigned int* __restrict__ interactionCount, real4 periodicBoxSize, real4 invPeriodicBoxSize, 
        real4 periodicBoxVecX, real4 periodicBoxVecY, real4 periodicBoxVecZ, unsigned int maxTiles, const real4* __restrict__ blockCenter,
        const real4* __restrict__ blockSize, const unsigned int* __restrict__ interactingAtoms, unsigned int maxSinglePairs,
        const int2* __restrict__ singlePairs
#endif
        PARAMETER_ARGUMENTS) {
    const unsigned int totalWarps = (blockDim.x*gridDim.x)/TILE_SIZE;
    const unsigned int warp = (blockIdx.x*blockDim.x+threadIdx.x)/TILE_SIZE; // global warpIndex
    const unsigned int tgx = threadIdx.x & (TILE_SIZE-1); // index within the warp
    const unsigned int tbx = threadIdx.x - tgx;           // block warpIndex
    mixed energy = 0;
    // used shared memory if the device cannot shuffle
#ifndef ENABLE_SHUFFLE
    __shared__ AtomData localData[THREAD_BLOCK_SIZE];
#endif

    // First loop: process tiles that contain exclusions.

    const unsigned int firstExclusionTile = FIRST_EXCLUSION_TILE+warp*(LAST_EXCLUSION_TILE-FIRST_EXCLUSION_TILE)/totalWarps;
    const unsigned int lastExclusionTile = FIRST_EXCLUSION_TILE+(warp+1)*(LAST_EXCLUSION_TILE-FIRST_EXCLUSION_TILE)/totalWarps;
    for (int pos = firstExclusionTile; pos < lastExclusionTile; pos++) {
        const ushort2 tileIndices = exclusionTiles[pos];
        const unsigned int x = tileIndices.x;
        const unsigned int y = tileIndices.y;
        real3 force = make_real3(0);
        real3 torque = make_real3(0);
        unsigned int atom1 = x*TILE_SIZE + tgx;
        real4 posq1 = posq[atom1];
        LOAD_ATOM1_PARAMETERS
        tileflags excl = exclusions[pos*TILE_SIZE+tgx];
        const bool hasExclusions = true;
        if (x == y) {
            // This tile is on the diagonal.
#ifdef ENABLE_SHUFFLE
            real4 shflPosq = posq1;
#else
            localData[threadIdx.x].x = posq1.x;
            localData[threadIdx.x].y = posq1.y;
            localData[threadIdx.x].z = posq1.z;
            localData[threadIdx.x].q = posq1.w;
            LOAD_LOCAL_PARAMETERS_FROM_1
#endif

            // we do not need to fetch parameters from global since this is a symmetric tile
            // instead we can broadcast the values using shuffle
            for (unsigned int j = 0; j < TILE_SIZE; j++) {
                int atom2 = tbx+j;
                real4 posq2;
#ifdef ENABLE_SHUFFLE
                BROADCAST_WARP_DATA
#else   
                posq2 = make_real4(localData[atom2].x, localData[atom2].y, localData[atom2].z, localData[atom2].q);
#endif
                real3 delta = make_real3(posq2.x-posq1.x, posq2.y-posq1.y, posq2.z-posq1.z);
#ifdef USE_PERIODIC
                APPLY_PERIODIC_TO_DELTA(delta)
#endif
                real r2 = delta.x*delta.x + delta.y*delta.y + delta.z*delta.z;
                real rInv = RSQRT(r2);
                real r = r2*rInv;
                LOAD_ATOM2_PARAMETERS
                atom2 = y*TILE_SIZE+j;
                real3 tempForce = make_real3(0);
                real3 tempTorque1 = make_real3(0);
                real3 tempTorque2 = make_real3(0);
                bool isExcluded = (atom1 >= NUM_ATOMS || atom2 >= NUM_ATOMS || !(excl & 0x1));
                real tempEnergy = 0.0f;
                const real interactionScale = 0.5f;
                COMPUTE_INTERACTION
                energy += 0.5f*tempEnergy;
                force += tempForce;
                torque += tempTorque1;
                excl >>= 1;
            }
        }
        else {
            // This is an off-diagonal tile.
            unsigned int j = y*TILE_SIZE + tgx;
            real4 shflPosq = posq[j];
#ifdef ENABLE_SHUFFLE
            real3 shflForce = make_real3(0);
            real3 shflTorque = make_real3(0);
#else
            localData[threadIdx.x].x = shflPosq.x;
            localData[threadIdx.x].y = shflPosq.y;
            localData[threadIdx.x].z = shflPosq.z;
            localData[threadIdx.x].q = shflPosq.w;
            localData[threadIdx.x].fx = 0.0f;
            localData[threadIdx.x].fy = 0.0f;
            localData[threadIdx.x].fz = 0.0f;
            localData[threadIdx.x].tx = 0.0f;
            localData[threadIdx.x].ty = 0.0f;
            localData[threadIdx.x].tz = 0.0f;
#endif
            DECLARE_LOCAL_PARAMETERS
            LOAD_LOCAL_PARAMETERS_FROM_GLOBAL
            excl = (excl >> tgx) | (excl << (TILE_SIZE - tgx));
            unsigned int tj = tgx;
            for (j = 0; j < TILE_SIZE; j++) {
                int atom2 = tbx+tj;
#ifdef ENABLE_SHUFFLE
                real4 posq2 = shflPosq;
#else
                real4 posq2 = make_real4(localData[atom2].x, localData[atom2].y, localData[atom2].z, localData[atom2].q);
#endif
                real3 delta = make_real3(posq2.x-posq1.x, posq2.y-posq1.y, posq2.z-posq1.z);
#ifdef USE_PERIODIC
                APPLY_PERIODIC_TO_DELTA(delta)
#endif
                real r2 = delta.x*delta.x + delta.y*delta.y + delta.z*delta.z;
                real rInv = RSQRT(r2);
                real r = r2*rInv;
                LOAD_ATOM2_PARAMETERS
                atom2 = y*TILE_SIZE+tj;
                real3 tempForce = make_real3(0);
                real3 tempTorque1 = make_real3(0);
                real3 tempTorque2 = make_real3(0);
                bool isExcluded = (atom1 >= NUM_ATOMS || atom2 >= NUM_ATOMS || !(excl & 0x1));
                real tempEnergy = 0.0f;
                const real interactionScale = 1.0f;
                COMPUTE_INTERACTION
                energy += tempEnergy;
                force += tempForce;
                torque += tempTorque1;
#ifdef ENABLE_SHUFFLE
                shflForce -= tempForce;
                shflTorque += tempTorque2;
                SHUFFLE_WARP_DATA
                shflTorque.x = real_shfl(shflTorque.x, tgx+1);
                shflTorque.y = real_shfl(shflTorque.y, tgx+1);
                shflTorque.z = real_shfl(shflTorque.z, tgx+1);
#else
                localData[tbx+tj].fx -= tempForce.x;
                localData[tbx+tj].fy -= tempForce.y;
                localData[tbx+tj].fz -= tempForce.z;
                localData[tbx+tj].tx += tempTorque2.x;
                localData[tbx+tj].ty += tempTorque2.y;
                localData[tbx+tj].tz += tempTorque2.z;
#endif
                excl >>= 1;
                // cycles the indices
                // 0 1 2 3 4 5 6 7 -> 1 2 3 4 5 6 7 0
                tj = (tj + 1) & (TILE_SIZE - 1);
            }
            const unsigned int offset = y*TILE_SIZE + tgx;
            // write results for off diagonal tiles
#ifdef ENABLE_SHUFFLE
            atomicAdd(&forceBuffers[offset], static_cast<unsigned long long>((long long) (shflForce.x*0x100000000)));
            atomicAdd(&forceBuffers[offset+PADDED_NUM_ATOMS], static_cast<unsigned long long>((long long) (shflForce.y*0x100000000)));
            atomicAdd(&forceBuffers[offset+2*PADDED_NUM_ATOMS], static_cast<unsigned long long>((long long) (shflForce.z*0x100000000)));
            atomicAdd(&torqueBuffers[offset], static_cast<unsigned long long>((long long) (shflTorque.x*0x100000000)));
            atomicAdd(&torqueBuffers[offset+PADDED_NUM_ATOMS], static_cast<unsigned long long>((long long) (shflTorque.y*0x100000000)));
            atomicAdd(&torqueBuffers[offset+2*PADDED_NUM_ATOMS], static_cast<unsigned long long>((long long) (shflTorque.z*0x100000000)));
#else
            atomicAdd(&forceBuffers[offset], static_cast<unsigned long long>((long long) (localData[threadIdx.x].fx*0x100000000)));
            atomicAdd(&forceBuffers[offset+PADDED_NUM_ATOMS], static_cast<unsigned long long>((long long) (localData[threadIdx.x].fy*0x100000000)));
            atomicAdd(&forceBuffers[offset+2*PADDED_NUM_ATOMS], static_cast<unsigned long long>((long long) (localData[threadIdx.x].fz*0x100000000)));
            atomicAdd(&torqueBuffers[offset], static_cast<unsigned long long>((long long) (localData[threadIdx.x].tx*0x100000000)));
            atomicAdd(&torqueBuffers[offset+PADDED_NUM_ATOMS], static_cast<unsigned long long>((long long) (localData[threadIdx.x].ty*0x100000000)));
            atomicAdd(&torqueBuffers[offset+2*PADDED_NUM_ATOMS], static_cast<unsigned long long>((long long) (localData[threadIdx.x].tz*0x100000000)));
#endif
        }
        // Write results for on and off diagonal tiles

        const unsigned int offset = x*TILE_SIZE + tgx;
        atomicAdd(&forceBuffers[offset], static_cast<unsigned long long>((long long) (force.x*0x100000000)));
        atomicAdd(&forceBuffers[offset+PADDED_NUM_ATOMS], static_cast<unsigned long long>((long long) (force.y*0x100000000)));
        atomicAdd(&forceBuffers[offset+2*PADDED_NUM_ATOMS], static_cast<unsigned long long>((long long) (force.z*0x100000000)));
        atomicAdd(&torqueBuffers[offset], static_cast<unsigned long long>((long long) (torque.x*0x100000000)));
        atomicAdd(&torqueBuffers[offset+PADDED_NUM_ATOMS], static_cast<unsigned long long>((long long) (torque.y*0x100000000)));
        atomicAdd(&torqueBuffers[offset+2*PADDED_NUM_ATOMS], static_cast<unsigned long long>((long long) (torque.z*0x100000000)));
    }

    // Second loop: tiles without exclusions, either from the neighbor list (with cutoff) or just enumerating all
    // of them (no cutoff).

#ifdef USE_CUTOFF
    const unsigned int numTiles = interactionCount[0];
    if (numTiles > maxTiles)
        return; // There wasn't enough memory for the neighbor list.
    int pos = (int) (numTiles > maxTiles ? startTileIndex+warp*(long long)numTileIndices/totalWarps : warp*(long long)numTiles/totalWarps);
    int end = (int) (numTiles > maxTiles ? startTileIndex+(warp+1)*(long long)numTileIndices/totalWarps : (warp+1)*(long long)numTiles/totalWarps);
#else
    const unsigned int numTiles = numTileIndices;
    int pos = (int) (startTileIndex+warp*(long long)numTiles/totalWarps);
    int end = (int) (startTileIndex+(warp+1)*(long long)numTiles/totalWarps);
#endif
    int skipBase = 0;
    int currentSkipIndex = tbx;
    // atomIndices can probably be shuffled as well
    // but it probably wouldn't make things any faster
    __shared__ int atomIndices[THREAD_BLOCK_SIZE];
    __shared__ volatile int skipTiles[THREAD_BLOCK_SIZE];
    skipTiles[threadIdx.x] = -1;
    
    while (pos < end) {
        const bool hasExclusions = false;
        real3 force = make_real3(0);
        real3 torque = make_real3(0);
        bool includeTile = true;

        // Extract the coordinates of this tile.
        int x, y;
        bool singlePeriodicCopy = false;
#ifdef USE_CUTOFF
        x = tiles[pos];
        real4 blockSizeX = blockSize[x];
        singlePeriodicCopy = (0.5f*periodicBoxSize.x-blockSizeX.x >= MAX_CUTOFF &&
                              0.5f*periodicBoxSize.y-blockSizeX.y >= MAX_CUTOFF &&
                              0.5f*periodicBoxSize.z-blockSizeX.z >= MAX_CUTOFF);
#else
        y = (int) floor(NUM_BLOCKS+0.5f-SQRT((NUM_BLOCKS+0.5f)*(NUM_BLOCKS+0.5f)-2*pos));
        x = (pos-y*NUM_BLOCKS+y*(y+1)/2);
        if (x < y || x >= NUM_BLOCKS) { // Occasionally happens due to roundoff error.
            y += (x < y ? -1 : 1);
            x = (pos-y*NUM_BLOCKS+y*(y+1)/2);
        }

        // Skip over tiles that have exclusions, since they were already processed.

        while (skipTiles[tbx+TILE_SIZE-1] < pos) {
            if (skipBase+tgx < NUM_TILES_WITH_EXCLUSIONS) {
                ushort2 tile = exclusionTiles[skipBase+tgx];
                skipTiles[threadIdx.x] = tile.x + tile.y*NUM_BLOCKS - tile.y*(tile.y+1)/2;
            }
            else
                skipTiles[threadIdx.x] = end;
            skipBase += TILE_SIZE;            
            currentSkipIndex = tbx;
        }
        while (skipTiles[currentSkipIndex] < pos)
            currentSkipIndex++;
        includeTile = (skipTiles[currentSkipIndex] != pos);
#endif
        if (includeTile) {
            unsigned int atom1 = x*TILE_SIZE + tgx;
            // Load atom data for this tile.
            real4 posq1 = posq[atom1];
            LOAD_ATOM1_PARAMETERS
            //const unsigned int localAtomIndex = threadIdx.x;
#ifdef USE_CUTOFF
            unsigned int j = interactingAtoms[pos*TILE_SIZE+tgx];
#else
            unsigned int j = y*TILE_SIZE + tgx;
#endif
            atomIndices[threadIdx.x] = j;
#ifdef ENABLE_SHUFFLE
            DECLARE_LOCAL_PARAMETERS
            real4 shflPosq;
            real3 shflForce = make_real3(0);
            real3 shflTorque = make_real3(0);
#endif
            if (j < PADDED_NUM_ATOMS) {
                // Load position of atom j from from global memory
#ifdef ENABLE_SHUFFLE
                shflPosq = posq[j];
#else
                localData[threadIdx.x].x = posq[j].x;
                localData[threadIdx.x].y = posq[j].y;
                localData[threadIdx.x].z = posq[j].z;
                localData[threadIdx.x].q = posq[j].w;
                localData[threadIdx.x].fx = 0.0f;
                localData[threadIdx.x].fy = 0.0f;
                localData[threadIdx.x].fz = 0.0f;
#endif                
                LOAD_LOCAL_PARAMETERS_FROM_GLOBAL
            }
            else {
#ifdef ENABLE_SHUFFLE
                shflPosq = make_real4(0, 0, 0, 0);
#else
                localData[threadIdx.x].x = 0;
                localData[threadIdx.x].y = 0;
                localData[threadIdx.x].z = 0;
#endif
            }
#ifdef USE_PERIODIC
            if (singlePeriodicCopy) {
                // The box is small enough that we can just translate all the atoms into a single periodic
                // box, then skip having to apply periodic boundary conditions later.
                real4 blockCenterX = blockCenter[x];
                APPLY_PERIODIC_TO_POS_WITH_CENTER(posq1, blockCenterX)
#ifdef ENABLE_SHUFFLE
                APPLY_PERIODIC_TO_POS_WITH_CENTER(shflPosq, blockCenterX)
#else
                APPLY_PERIODIC_TO_POS_WITH_CENTER(localData[threadIdx.x], blockCenterX)
#endif
                unsigned int tj = tgx;
                for (j = 0; j < TILE_SIZE; j++) {
                    int atom2 = tbx+tj;
#ifdef ENABLE_SHUFFLE
                    real4 posq2 = shflPosq; 
#else
                    real4 posq2 = make_real4(localData[atom2].x, localData[atom2].y, localData[atom2].z, localData[atom2].q);
#endif
                    real3 delta = make_real3(posq2.x-posq1.x, posq2.y-posq1.y, posq2.z-posq1.z);
                    real r2 = delta.x*delta.x + delta.y*delta.y + delta.z*delta.z;
                    real rInv = RSQRT(r2);
                    real r = r2*rInv;
                    LOAD_ATOM2_PARAMETERS
                    atom2 = atomIndices[tbx+tj];
                    real3 tempForce = make_real3(0);
                    real3 tempTorque1 = make_real3(0);
                    real3 tempTorque2 = make_real3(0);
                    bool isExcluded = (atom1 >= NUM_ATOMS || atom2 >= NUM_ATOMS);
                    real tempEnergy = 0.0f;
                    const real interactionScale = 1.0f;
                    COMPUTE_INTERACTION
                    energy += tempEnergy;
                    force += tempForce;
                    torque += tempTorque1;
#ifdef ENABLE_SHUFFLE
                    shflForce -= tempForce;
                    shflTorque += tempTorque2;
                    SHUFFLE_WARP_DATA
                    shflTorque.x = real_shfl(shflTorque.x, tgx+1);
                    shflTorque.y = real_shfl(shflTorque.y, tgx+1);
                    shflTorque.z = real_shfl(shflTorque.z, tgx+1);
#else
                    localData[tbx+tj].fx -= tempForce.x;
                    localData[tbx+tj].fy -= tempForce.y;
                    localData[tbx+tj].fz -= tempForce.z;
                    localData[tbx+tj].tx += tempTorque2.x;
                    localData[tbx+tj].ty += tempTorque2.y;
                    localData[tbx+tj].tz += tempTorque2.z;
#endif
                    tj = (tj + 1) & (TILE_SIZE - 1);
                }
            }
            else
#endif
            {
                // We need to apply periodic boundary conditions separately for each interaction.
                unsigned int tj = tgx;
                for (j = 0; j < TILE_SIZE; j++) {
                    int atom2 = tbx+tj;
#ifdef ENABLE_SHUFFLE
                    real4 posq2 = shflPosq;
#else
                    real4 posq2 = make_real4(localData[atom2].x, localData[atom2].y, localData[atom2].z, localData[atom2].q);
#endif
                    real3 delta = make_real3(posq2.x-posq1.x, posq2.y-posq1.y, posq2.z-posq1.z);
#ifdef USE_PERIODIC
                    APPLY_PERIODIC_TO_DELTA(delta)
#endif
                    real r2 = delta.x*delta.x + delta.y*delta.y + delta.z*delta.z;
                    real rInv = RSQRT(r2);
                    real r = r2*rInv;
                    LOAD_ATOM2_PARAMETERS
                    atom2 = atomIndices[tbx+tj];
                    real3 tempForce = make_real3(0);
                    real3 tempTorque1 = make_real3(0);
                    real3 tempTorque2 = make_real3(0);
                    bool isExcluded = (atom1 >= NUM_ATOMS || atom2 >= NUM_ATOMS);
                    real tempEnergy = 0.0f;
                    const real interactionScale = 1.0f;
                    COMPUTE_INTERACTION
                    energy += tempEnergy;
                    force += tempForce;
                    torque += tempTorque1;
#ifdef ENABLE_SHUFFLE
                    shflForce -= tempForce;
                    shflTorque += tempTorque2;
                    SHUFFLE_WARP_DATA
                    shflTorque.x = real_shfl(shflTorque.x, tgx+1);
                    shflTorque.y = real_shfl(shflTorque.y, tgx+1);
                    shflTorque.z = real_shfl(shflTorque.z, tgx+1);
#else
                    localData[tbx+tj].fx -= tempForce.x;
                    localData[tbx+tj].fy -= tempForce.y;
                    localData[tbx+tj].fz -= tempForce.z;
                    localData[tbx+tj].tx += tempTorque.x;
                    localData[tbx+tj].ty += tempTorque.y;
                    localData[tbx+tj].tz += tempTorque.z;
#endif
                    tj = (tj + 1) & (TILE_SIZE - 1);
                }
            }

            // Write results.

            atomicAdd(&forceBuffers[atom1], static_cast<unsigned long long>((long long) (force.x*0x100000000)));
            atomicAdd(&forceBuffers[atom1+PADDED_NUM_ATOMS], static_cast<unsigned long long>((long long) (force.y*0x100000000)));
            atomicAdd(&forceBuffers[atom1+2*PADDED_NUM_ATOMS], static_cast<unsigned long long>((long long) (force.z*0x100000000)));
            atomicAdd(&torqueBuffers[atom1], static_cast<unsigned long long>((long long) (torque.x*0x100000000)));
            atomicAdd(&torqueBuffers[atom1+PADDED_NUM_ATOMS], static_cast<unsigned long long>((long long) (torque.y*0x100000000)));
            atomicAdd(&torqueBuffers[atom1+2*PADDED_NUM_ATOMS], static_cast<unsigned long long>((long long) (torque.z*0x100000000)));
#ifdef USE_CUTOFF
            unsigned int atom2 = atomIndices[threadIdx.x];
#else
            unsigned int atom2 = y*TILE_SIZE + tgx;
#endif
            if (atom2 < PADDED_NUM_ATOMS) {
#ifdef ENABLE_SHUFFLE
                atomicAdd(&forceBuffers[atom2], static_cast<unsigned long long>((long long) (shflForce.x*0x100000000)));
                atomicAdd(&forceBuffers[atom2+PADDED_NUM_ATOMS], static_cast<unsigned long long>((long long) (shflForce.y*0x100000000)));
                atomicAdd(&forceBuffers[atom2+2*PADDED_NUM_ATOMS], static_cast<unsigned long long>((long long) (shflForce.z*0x100000000)));
                atomicAdd(&torqueBuffers[atom2], static_cast<unsigned long long>((long long) (shflTorque.x*0x100000000)));
                atomicAdd(&torqueBuffers[atom2+PADDED_NUM_ATOMS], static_cast<unsigned long long>((long long) (shflTorque.y*0x100000000)));
                atomicAdd(&torqueBuffers[atom2+2*PADDED_NUM_ATOMS], static_cast<unsigned long long>((long long) (shflTorque.z*0x100000000)));
#else
                atomicAdd(&forceBuffers[atom2], static_cast<unsigned long long>((long long) (localData[threadIdx.x].fx*0x100000000)));
                atomicAdd(&forceBuffers[atom2+PADDED_NUM_ATOMS], static_cast<unsigned long long>((long long) (localData[threadIdx.x].fy*0x100000000)));
                atomicAdd(&forceBuffers[atom2+2*PADDED_NUM_ATOMS], static_cast<unsigned long long>((long long) (localData[threadIdx.x].fz*0x100000000)));
                atomicAdd(&torqueBuffers[atom2], static_cast<unsigned long long>((long long) (localData[threadIdx.x].tx*0x100000000)));
                atomicAdd(&torqueBuffers[atom2+PADDED_NUM_ATOMS], static_cast<unsigned long long>((long long) (localData[threadIdx.x].ty*0x100000000)));
                atomicAdd(&torqueBuffers[atom2+2*PADDED_NUM_ATOMS], static_cast<unsigned long long>((long long) (localData[threadIdx.x].tz*0x100000000)));
#endif
            }
        }
        pos++;
    }
#ifdef INCLUDE_ENERGY
    energyBuffer[blockIdx.x*blockDim.x+threadIdx.x] += energy;
#endif
}