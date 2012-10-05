#define TILE_SIZE 32
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

/**
 * Compute nonbonded interactions.
 */
extern "C" __global__ void computeNonbonded(
        unsigned long long* __restrict__ forceBuffers, real* __restrict__ energyBuffer, const real4* __restrict__ posq, const unsigned int* __restrict__ exclusions,
        const unsigned int* __restrict__ exclusionIndices, const unsigned int* __restrict__ exclusionRowIndices,
        unsigned int startTileIndex, unsigned int numTileIndices
#ifdef USE_CUTOFF
        , const ushort2* __restrict__ tiles, const unsigned int* __restrict__ interactionCount, real4 periodicBoxSize, real4 invPeriodicBoxSize, unsigned int maxTiles, const unsigned int* __restrict__ interactionFlags
#endif
        PARAMETER_ARGUMENTS) {
    unsigned int totalWarps = (blockDim.x*gridDim.x)/TILE_SIZE;
    unsigned int warp = (blockIdx.x*blockDim.x+threadIdx.x)/TILE_SIZE;
#ifdef USE_CUTOFF
    const unsigned int numTiles = interactionCount[0];
    unsigned int pos = (numTiles > maxTiles ? startTileIndex+warp*numTileIndices/totalWarps : warp*numTiles/totalWarps);
    unsigned int end = (numTiles > maxTiles ? startTileIndex+(warp+1)*numTileIndices/totalWarps : (warp+1)*numTiles/totalWarps);
#else
    const unsigned int numTiles = numTileIndices;
    unsigned int pos = startTileIndex+warp*numTiles/totalWarps;
    unsigned int end = startTileIndex+(warp+1)*numTiles/totalWarps;
#endif
    real energy = 0.0f;
    __shared__ AtomData localData[THREAD_BLOCK_SIZE];
    __shared__ unsigned int exclusionRange[2*WARPS_PER_GROUP];
    __shared__ int exclusionIndex[WARPS_PER_GROUP];
#ifndef ENABLE_SHUFFLE
    __shared__ real tempBuffer[3*THREAD_BLOCK_SIZE];
#endif
    
    do {
        // Extract the coordinates of this tile
        const unsigned int tgx = threadIdx.x & (TILE_SIZE-1);
        const unsigned int tbx = threadIdx.x - tgx;
        const unsigned int localGroupIndex = threadIdx.x/TILE_SIZE;
        unsigned int x, y;
        real3 force = make_real3(0);
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

                const unsigned int localAtomIndex = threadIdx.x;
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
                    real4 posq2 = make_real4(localData[atom2].x, localData[atom2].y, localData[atom2].z, localData[atom2].q);
                    real3 delta = make_real3(posq2.x-posq1.x, posq2.y-posq1.y, posq2.z-posq1.z);
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
                    real dEdR = 0.0f;
#else
                    real3 dEdR1 = make_real3(0);
                    real3 dEdR2 = make_real3(0);
#endif
                    real tempEnergy = 0.0f;
                    COMPUTE_INTERACTION
                    energy += 0.5f*tempEnergy;
#ifdef USE_SYMMETRIC
                    force.x -= delta.x*dEdR;
                    force.y -= delta.y*dEdR;
                    force.z -= delta.z*dEdR;
#else
                    force.x -= dEdR1.x;
                    force.y -= dEdR1.y;
                    force.z -= dEdR1.z;
#endif
#ifdef USE_EXCLUSIONS
                    excl >>= 1;
#endif
                }
            }
            else {
                // This is an off-diagonal tile.

                const unsigned int localAtomIndex = threadIdx.x;
                unsigned int j = y*TILE_SIZE + tgx;
                real4 tempPosq = posq[j];
                localData[localAtomIndex].x = tempPosq.x;
                localData[localAtomIndex].y = tempPosq.y;
                localData[localAtomIndex].z = tempPosq.z;
                localData[localAtomIndex].q = tempPosq.w;
                LOAD_LOCAL_PARAMETERS_FROM_GLOBAL
                localData[localAtomIndex].fx = 0.0f;
                localData[localAtomIndex].fy = 0.0f;
                localData[localAtomIndex].fz = 0.0f;
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
                                int bufferIndex = 3*threadIdx.x;
#ifdef USE_SYMMETRIC
                                real dEdR = 0;
#else
                                real3 dEdR1 = make_real3(0);
                                real3 dEdR2 = make_real3(0);
#endif
                                real tempEnergy = 0.0f;
                                real4 posq2 = make_real4(localData[atom2].x, localData[atom2].y, localData[atom2].z, localData[atom2].q);
                                real3 delta = make_real3(posq2.x-posq1.x, posq2.y-posq1.y, posq2.z-posq1.z);
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
#ifdef ENABLE_SHUFFLE
    #ifdef USE_SYMMETRIC
                                delta *= dEdR;
                                force.x -= delta.x;
                                force.y -= delta.y;
                                force.z -= delta.z;
                                for (int i = 16; i >= 1; i /= 2) {
                                    delta.x += __shfl_xor(delta.x, i, 32);
                                    delta.y += __shfl_xor(delta.y, i, 32);
                                    delta.z += __shfl_xor(delta.z, i, 32);
                                }
                                if (tgx == 0) {
                                    localData[tbx+j].fx += delta.x;
                                    localData[tbx+j].fy += delta.y;
                                    localData[tbx+j].fz += delta.z;
                                }
    #else
                                force.x -= dEdR1.x;
                                force.y -= dEdR1.y;
                                force.z -= dEdR1.z;
                                for (int i = 16; i >= 1; i /= 2) {
                                    dEdR2.x += __shfl_xor(dEdR2.x, i, 32);
                                    dEdR2.y += __shfl_xor(dEdR2.y, i, 32);
                                    dEdR2.z += __shfl_xor(dEdR2.z, i, 32);
                                }
                                if (tgx == 0) {
                                    localData[tbx+j].fx += dEdR2.x;
                                    localData[tbx+j].fy += dEdR2.y;
                                    localData[tbx+j].fz += dEdR2.z;
                                }
    #endif
#else
    #ifdef USE_SYMMETRIC
                                delta *= dEdR;
                                force.x -= delta.x;
                                force.y -= delta.y;
                                force.z -= delta.z;
                                tempBuffer[bufferIndex] = delta.x;
                                tempBuffer[bufferIndex+1] = delta.y;
                                tempBuffer[bufferIndex+2] = delta.z;
    #else
                                force.x -= dEdR1.x;
                                force.y -= dEdR1.y;
                                force.z -= dEdR1.z;
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
#endif
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
                        real4 posq2 = make_real4(localData[atom2].x, localData[atom2].y, localData[atom2].z, localData[atom2].q);
                        real3 delta = make_real3(posq2.x-posq1.x, posq2.y-posq1.y, posq2.z-posq1.z);
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
                            real dEdR = 0.0f;
#else
                            real3 dEdR1 = make_real3(0);
                            real3 dEdR2 = make_real3(0);
#endif
                            real tempEnergy = 0.0f;
                            COMPUTE_INTERACTION
                            energy += tempEnergy;
#ifdef USE_SYMMETRIC
                            delta *= dEdR;
                            force.x -= delta.x;
                            force.y -= delta.y;
                            force.z -= delta.z;
                            localData[tbx+tj].fx += delta.x;
                            localData[tbx+tj].fy += delta.y;
                            localData[tbx+tj].fz += delta.z;
#else
                            force.x -= dEdR1.x;
                            force.y -= dEdR1.y;
                            force.z -= dEdR1.z;
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
        
        // Write results.
        
        if (pos < end) {
            const unsigned int offset = x*TILE_SIZE + tgx;
            atomicAdd(&forceBuffers[offset], static_cast<unsigned long long>((long long) (force.x*0xFFFFFFFF)));
            atomicAdd(&forceBuffers[offset+PADDED_NUM_ATOMS], static_cast<unsigned long long>((long long) (force.y*0xFFFFFFFF)));
            atomicAdd(&forceBuffers[offset+2*PADDED_NUM_ATOMS], static_cast<unsigned long long>((long long) (force.z*0xFFFFFFFF)));
        }
        if (pos < end && x != y) {
            const unsigned int offset = y*TILE_SIZE + tgx;
            atomicAdd(&forceBuffers[offset], static_cast<unsigned long long>((long long) (localData[threadIdx.x].fx*0xFFFFFFFF)));
            atomicAdd(&forceBuffers[offset+PADDED_NUM_ATOMS], static_cast<unsigned long long>((long long) (localData[threadIdx.x].fy*0xFFFFFFFF)));
            atomicAdd(&forceBuffers[offset+2*PADDED_NUM_ATOMS], static_cast<unsigned long long>((long long) (localData[threadIdx.x].fz*0xFFFFFFFF)));
        }
        pos++;
    } while (pos < end);
    energyBuffer[blockIdx.x*blockDim.x+threadIdx.x] += energy;
}
