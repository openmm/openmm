#define TILE_SIZE 32
#define WARPS_PER_GROUP (THREAD_BLOCK_SIZE/TILE_SIZE)

typedef struct {
    real4 posq;
    real3 force, dipole, inducedDipole, inducedDipolePolar;
    real quadrupoleXX, quadrupoleXY, quadrupoleXZ;
    real quadrupoleYY, quadrupoleYZ;
    float thole, damp;
} AtomData;

__device__ void computeOneInteractionF1(AtomData& atom1, volatile AtomData& atom2, float dScale, float pScale, float mScale, real& energy, real3& outputForce);
__device__ void computeOneInteractionT1(AtomData& atom1, volatile AtomData& atom2, float dScale, float pScale, float mScale, real3& outputForce);
__device__ void computeOneInteractionT3(AtomData& atom1, volatile AtomData& atom2, float dScale, float pScale, float mScale, real3& outputForce);

inline __device__ void loadAtomData(AtomData& data, int atom, const real4* __restrict__ posq, const real* __restrict__ labFrameDipole,
        const real* __restrict__ labFrameQuadrupole, const real* __restrict__ inducedDipole, const real* __restrict__ inducedDipolePolar, const float2* __restrict__ dampingAndThole) {
    data.posq = posq[atom];
    data.dipole.x = labFrameDipole[atom*3];
    data.dipole.y = labFrameDipole[atom*3+1];
    data.dipole.z = labFrameDipole[atom*3+2];
    data.quadrupoleXX = labFrameQuadrupole[atom*5];
    data.quadrupoleXY = labFrameQuadrupole[atom*5+1];
    data.quadrupoleXZ = labFrameQuadrupole[atom*5+2];
    data.quadrupoleYY = labFrameQuadrupole[atom*5+3];
    data.quadrupoleYZ = labFrameQuadrupole[atom*5+4];
    data.inducedDipole.x = inducedDipole[atom*3];
    data.inducedDipole.y = inducedDipole[atom*3+1];
    data.inducedDipole.z = inducedDipole[atom*3+2];
    data.inducedDipolePolar.x = inducedDipolePolar[atom*3];
    data.inducedDipolePolar.y = inducedDipolePolar[atom*3+1];
    data.inducedDipolePolar.z = inducedDipolePolar[atom*3+2];
    float2 temp = dampingAndThole[atom];
    data.damp = temp.x;
    data.thole = temp.y;
}

__device__ real computeDScaleFactor(unsigned int polarizationGroup, int index) {
    return (polarizationGroup & 1<<index ? 0 : 1);
}

__device__ float computeMScaleFactor(uint2 covalent, int index) {
    int mask = 1<<index;
    bool x = (covalent.x & mask);
    bool y = (covalent.y & mask);
    return (x ? (y ? 0.0f : 0.4f) : (y ? 0.8f : 1.0f));
}

__device__ float computePScaleFactor(uint2 covalent, unsigned int polarizationGroup, int index) {
    int mask = 1<<index;
    bool x = (covalent.x & mask);
    bool y = (covalent.y & mask);
    bool p = (polarizationGroup & mask);
    return (x && y ? 0.0f : (x && p ? 0.5f : 1.0f));
}

/**
 * Compute electrostatic interactions.
 */
extern "C" __global__ void computeElectrostatics(
        unsigned long long* __restrict__ forceBuffers, unsigned long long* __restrict__ torqueBuffers, real* __restrict__ energyBuffer,
        const real4* __restrict__ posq, const unsigned int* __restrict__ exclusionIndices, const unsigned int* __restrict__ exclusionRowIndices,
        const uint2* __restrict__ covalentFlags, const unsigned int* __restrict__ polarizationGroupFlags, unsigned int startTileIndex, unsigned int numTileIndices,
#ifdef USE_CUTOFF
        const ushort2* __restrict__ tiles, const unsigned int* __restrict__ interactionCount, real4 periodicBoxSize, real4 invPeriodicBoxSize, unsigned int maxTiles, const unsigned int* __restrict__ interactionFlags,
#endif
        const real* __restrict__ labFrameDipole, const real* __restrict__ labFrameQuadrupole, const real* __restrict__ inducedDipole,
        const real* __restrict__ inducedDipolePolar, const float2* __restrict__ dampingAndThole) {
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
    real energy = 0;
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
        AtomData data;
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
            loadAtomData(data, atom1, posq, labFrameDipole, labFrameQuadrupole, inducedDipole, inducedDipolePolar, dampingAndThole);
            data.force = make_real3(0);
            
            // Locate the exclusion data for this tile.

            if (tgx < 2)
                exclusionRange[2*localGroupIndex+tgx] = exclusionRowIndices[x+tgx];
            if (tgx == 0)
                exclusionIndex[localGroupIndex] = -1;
            for (unsigned int i = exclusionRange[2*localGroupIndex]+tgx; i < exclusionRange[2*localGroupIndex+1]; i += TILE_SIZE)
                if (exclusionIndices[i] == y)
                    exclusionIndex[localGroupIndex] = i*TILE_SIZE;
            bool hasExclusions = (exclusionIndex[localGroupIndex] > -1);
            if (pos >= end)
                ; // This warp is done.
            else if (x == y) {
                // This tile is on the diagonal.

                localData[threadIdx.x].posq = data.posq;
                localData[threadIdx.x].dipole = data.dipole;
                localData[threadIdx.x].quadrupoleXX = data.quadrupoleXX;
                localData[threadIdx.x].quadrupoleXY = data.quadrupoleXY;
                localData[threadIdx.x].quadrupoleXZ = data.quadrupoleXZ;
                localData[threadIdx.x].quadrupoleYY = data.quadrupoleYY;
                localData[threadIdx.x].quadrupoleYZ = data.quadrupoleYZ;
                localData[threadIdx.x].inducedDipole = data.inducedDipole;
                localData[threadIdx.x].inducedDipolePolar = data.inducedDipolePolar;
                localData[threadIdx.x].thole = data.thole;
                localData[threadIdx.x].damp = data.damp;
                uint2 covalent = covalentFlags[exclusionIndex[localGroupIndex]+tgx];
                unsigned int polarizationGroup = polarizationGroupFlags[exclusionIndex[localGroupIndex]+tgx];
                
                // Compute forces.
                
                for (unsigned int j = 0; j < TILE_SIZE; j++) {
                    int atom2 = y*TILE_SIZE+j;
                    if (atom1 != atom2 && atom1 < NUM_ATOMS && atom2 < NUM_ATOMS) {
                        real3 tempForce;
                        real tempEnergy;
                        float d = computeDScaleFactor(polarizationGroup, j);
                        float p = computePScaleFactor(covalent, polarizationGroup, j);
                        float m = computeMScaleFactor(covalent, j);
                        computeOneInteractionF1(data, localData[tbx+j], d, p, m, tempEnergy, tempForce);
                        data.force += tempForce;
                        energy += 0.5f*tempEnergy;
                    }
                }
                data.force *= ENERGY_SCALE_FACTOR;
                atomicAdd(&forceBuffers[atom1], static_cast<unsigned long long>((long long) (data.force.x*0x100000000)));
                atomicAdd(&forceBuffers[atom1+PADDED_NUM_ATOMS], static_cast<unsigned long long>((long long) (data.force.y*0x100000000)));
                atomicAdd(&forceBuffers[atom1+2*PADDED_NUM_ATOMS], static_cast<unsigned long long>((long long) (data.force.z*0x100000000)));
                
                // Compute torques.
                
                data.force = make_real3(0);
                for (unsigned int j = 0; j < TILE_SIZE; j++) {
                    int atom2 = y*TILE_SIZE+j;
                    if (atom1 != atom2 && atom1 < NUM_ATOMS && atom2 < NUM_ATOMS) {
                        real3 tempForce;
                        float d = computeDScaleFactor(polarizationGroup, j);
                        float p = computePScaleFactor(covalent, polarizationGroup, j);
                        float m = computeMScaleFactor(covalent, j);
                        computeOneInteractionT1(data, localData[tbx+j], d, p, m, tempForce);
                        data.force += tempForce;
                    }
                }
                data.force *= ENERGY_SCALE_FACTOR;
                atomicAdd(&torqueBuffers[atom1], static_cast<unsigned long long>((long long) (data.force.x*0x100000000)));
                atomicAdd(&torqueBuffers[atom1+PADDED_NUM_ATOMS], static_cast<unsigned long long>((long long) (data.force.y*0x100000000)));
                atomicAdd(&torqueBuffers[atom1+2*PADDED_NUM_ATOMS], static_cast<unsigned long long>((long long) (data.force.z*0x100000000)));
            }
            else {
                // This is an off-diagonal tile.

                unsigned int j = y*TILE_SIZE + tgx;
                loadAtomData(localData[threadIdx.x], j, posq, labFrameDipole, labFrameQuadrupole, inducedDipole, inducedDipolePolar, dampingAndThole);
                localData[threadIdx.x].force = make_real3(0);
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
                                int atom2 = tbx+j;
                                real3 delta = make_real3(localData[atom2].posq.x-data.posq.x, localData[atom2].posq.y-data.posq.y, localData[atom2].posq.z-data.posq.z);
#ifdef USE_PERIODIC
                                delta.x -= floor(delta.x*invPeriodicBoxSize.x+0.5f)*periodicBoxSize.x;
                                delta.y -= floor(delta.y*invPeriodicBoxSize.y+0.5f)*periodicBoxSize.y;
                                delta.z -= floor(delta.z*invPeriodicBoxSize.z+0.5f)*periodicBoxSize.z;
#endif
                                real3 tempForce;
                                real tempEnergy;
                                computeOneInteractionF1(data, localData[atom2], 1, 1, 1, tempEnergy, tempForce);
                                data.force += tempForce;
                                localData[atom2].force -= tempForce;
                                energy += tempEnergy;
                                if (atom1 < NUM_ATOMS && atom2 < NUM_ATOMS) {
#ifdef ENABLE_SHUFFLE
                                    for (int i = 16; i >= 1; i /= 2) {
                                        tempForce.x += __shfl_xor(tempForce.x, i, 32);
                                        tempForce.y += __shfl_xor(tempForce.y, i, 32);
                                        tempForce.z += __shfl_xor(tempForce.z, i, 32);
                                    }
                                    if (tgx == 0)
                                        localData[atom2].force -= tempForce;
#else
                                    int bufferIndex = 3*threadIdx.x;
                                    tempBuffer[bufferIndex] = tempForce.x;
                                    tempBuffer[bufferIndex+1] = tempForce.y;
                                    tempBuffer[bufferIndex+2] = tempForce.z;
                                    if (tgx % 4 == 0) {
                                        tempBuffer[bufferIndex] += tempBuffer[bufferIndex+3]+tempBuffer[bufferIndex+6]+tempBuffer[bufferIndex+9];
                                        tempBuffer[bufferIndex+1] += tempBuffer[bufferIndex+4]+tempBuffer[bufferIndex+7]+tempBuffer[bufferIndex+10];
                                        tempBuffer[bufferIndex+2] += tempBuffer[bufferIndex+5]+tempBuffer[bufferIndex+8]+tempBuffer[bufferIndex+11];
                                    }
                                    if (tgx == 0) {
                                        localData[atom2].force.x -= tempBuffer[bufferIndex]+tempBuffer[bufferIndex+12]+tempBuffer[bufferIndex+24]+tempBuffer[bufferIndex+36]+tempBuffer[bufferIndex+48]+tempBuffer[bufferIndex+60]+tempBuffer[bufferIndex+72]+tempBuffer[bufferIndex+84];
                                        localData[atom2].force.y -= tempBuffer[bufferIndex+1]+tempBuffer[bufferIndex+13]+tempBuffer[bufferIndex+25]+tempBuffer[bufferIndex+37]+tempBuffer[bufferIndex+49]+tempBuffer[bufferIndex+61]+tempBuffer[bufferIndex+73]+tempBuffer[bufferIndex+85];
                                        localData[atom2].force.z -= tempBuffer[bufferIndex+2]+tempBuffer[bufferIndex+14]+tempBuffer[bufferIndex+26]+tempBuffer[bufferIndex+38]+tempBuffer[bufferIndex+50]+tempBuffer[bufferIndex+62]+tempBuffer[bufferIndex+74]+tempBuffer[bufferIndex+86];
                                    }
#endif
                                }
                            }
                        }
                        data.force *= ENERGY_SCALE_FACTOR;
                        localData[threadIdx.x].force *= ENERGY_SCALE_FACTOR;
                        if (pos < end) {
                            unsigned int offset = x*TILE_SIZE + tgx;
                            atomicAdd(&forceBuffers[offset], static_cast<unsigned long long>((long long) (data.force.x*0x100000000)));
                            atomicAdd(&forceBuffers[offset+PADDED_NUM_ATOMS], static_cast<unsigned long long>((long long) (data.force.y*0x100000000)));
                            atomicAdd(&forceBuffers[offset+2*PADDED_NUM_ATOMS], static_cast<unsigned long long>((long long) (data.force.z*0x100000000)));
                            offset = y*TILE_SIZE + tgx;
                            atomicAdd(&forceBuffers[offset], static_cast<unsigned long long>((long long) (localData[threadIdx.x].force.x*0x100000000)));
                            atomicAdd(&forceBuffers[offset+PADDED_NUM_ATOMS], static_cast<unsigned long long>((long long) (localData[threadIdx.x].force.y*0x100000000)));
                            atomicAdd(&forceBuffers[offset+2*PADDED_NUM_ATOMS], static_cast<unsigned long long>((long long) (localData[threadIdx.x].force.z*0x100000000)));
                        }
                        
                        // Compute torques.

                        data.force = make_real3(0);
                        localData[threadIdx.x].force = make_real3(0);
                        for (j = 0; j < TILE_SIZE; j++) {
                            if ((flags&(1<<j)) != 0) {
                                int atom2 = tbx+j;
                                real3 delta = make_real3(localData[atom2].posq.x-data.posq.x, localData[atom2].posq.y-data.posq.y, localData[atom2].posq.z-data.posq.z);
#ifdef USE_PERIODIC
                                delta.x -= floor(delta.x*invPeriodicBoxSize.x+0.5f)*periodicBoxSize.x;
                                delta.y -= floor(delta.y*invPeriodicBoxSize.y+0.5f)*periodicBoxSize.y;
                                delta.z -= floor(delta.z*invPeriodicBoxSize.z+0.5f)*periodicBoxSize.z;
#endif
                                real3 tempForce;
                                computeOneInteractionT1(data, localData[atom2], 1, 1, 1, tempForce);
                                data.force += tempForce;
                                computeOneInteractionT3(data, localData[atom2], 1, 1, 1, tempForce);
                                if (atom1 < NUM_ATOMS && atom2 < NUM_ATOMS) {
#ifdef ENABLE_SHUFFLE
                                    for (int i = 16; i >= 1; i /= 2) {
                                        tempForce.x += __shfl_xor(tempForce.x, i, 32);
                                        tempForce.y += __shfl_xor(tempForce.y, i, 32);
                                        tempForce.z += __shfl_xor(tempForce.z, i, 32);
                                    }
                                    if (tgx == 0)
                                        localData[atom2].force -= tempForce;
#else
                                    int bufferIndex = 3*threadIdx.x;
                                    tempBuffer[bufferIndex] = tempForce.x;
                                    tempBuffer[bufferIndex+1] = tempForce.y;
                                    tempBuffer[bufferIndex+2] = tempForce.z;
                                    if (tgx % 4 == 0) {
                                        tempBuffer[bufferIndex] += tempBuffer[bufferIndex+3]+tempBuffer[bufferIndex+6]+tempBuffer[bufferIndex+9];
                                        tempBuffer[bufferIndex+1] += tempBuffer[bufferIndex+4]+tempBuffer[bufferIndex+7]+tempBuffer[bufferIndex+10];
                                        tempBuffer[bufferIndex+2] += tempBuffer[bufferIndex+5]+tempBuffer[bufferIndex+8]+tempBuffer[bufferIndex+11];
                                    }
                                    if (tgx == 0) {
                                        localData[atom2].force.x += tempBuffer[bufferIndex]+tempBuffer[bufferIndex+12]+tempBuffer[bufferIndex+24]+tempBuffer[bufferIndex+36]+tempBuffer[bufferIndex+48]+tempBuffer[bufferIndex+60]+tempBuffer[bufferIndex+72]+tempBuffer[bufferIndex+84];
                                        localData[atom2].force.y += tempBuffer[bufferIndex+1]+tempBuffer[bufferIndex+13]+tempBuffer[bufferIndex+25]+tempBuffer[bufferIndex+37]+tempBuffer[bufferIndex+49]+tempBuffer[bufferIndex+61]+tempBuffer[bufferIndex+73]+tempBuffer[bufferIndex+85];
                                        localData[atom2].force.z += tempBuffer[bufferIndex+2]+tempBuffer[bufferIndex+14]+tempBuffer[bufferIndex+26]+tempBuffer[bufferIndex+38]+tempBuffer[bufferIndex+50]+tempBuffer[bufferIndex+62]+tempBuffer[bufferIndex+74]+tempBuffer[bufferIndex+86];
                                    }
#endif
                                }
                            }
                        }
                        data.force *= ENERGY_SCALE_FACTOR;
                        localData[threadIdx.x].force *= ENERGY_SCALE_FACTOR;
                        if (pos < end) {
                            unsigned int offset = x*TILE_SIZE + tgx;
                            atomicAdd(&torqueBuffers[offset], static_cast<unsigned long long>((long long) (data.force.x*0x100000000)));
                            atomicAdd(&torqueBuffers[offset+PADDED_NUM_ATOMS], static_cast<unsigned long long>((long long) (data.force.y*0x100000000)));
                            atomicAdd(&torqueBuffers[offset+2*PADDED_NUM_ATOMS], static_cast<unsigned long long>((long long) (data.force.z*0x100000000)));
                            offset = y*TILE_SIZE + tgx;
                            atomicAdd(&torqueBuffers[offset], static_cast<unsigned long long>((long long) (localData[threadIdx.x].force.x*0x100000000)));
                            atomicAdd(&torqueBuffers[offset+PADDED_NUM_ATOMS], static_cast<unsigned long long>((long long) (localData[threadIdx.x].force.y*0x100000000)));
                            atomicAdd(&torqueBuffers[offset+2*PADDED_NUM_ATOMS], static_cast<unsigned long long>((long long) (localData[threadIdx.x].force.z*0x100000000)));
                        }
                    }
                }
                else
#endif
                {
                    // Compute the full set of interactions in this tile.

                    uint2 covalent = (hasExclusions ? covalentFlags[exclusionIndex[localGroupIndex]+tgx] : make_uint2(0, 0));
                    unsigned int polarizationGroup = (hasExclusions ? polarizationGroupFlags[exclusionIndex[localGroupIndex]+tgx] : 0);
                    
                    // Compute forces.
                    
                    unsigned int tj = tgx;
                    for (j = 0; j < TILE_SIZE; j++) {
                        int atom2 = y*TILE_SIZE+tj;
                        if (atom1 < NUM_ATOMS && atom2 < NUM_ATOMS) {
                            real3 tempForce;
                            real tempEnergy;
                            float d = computeDScaleFactor(polarizationGroup, tj);
                            float p = computePScaleFactor(covalent, polarizationGroup, tj);
                            float m = computeMScaleFactor(covalent, tj);
                            computeOneInteractionF1(data, localData[tbx+tj], d, p, m, tempEnergy, tempForce);
                            data.force += tempForce;
                            localData[tbx+tj].force -= tempForce;
                            energy += tempEnergy;
                        }
                        tj = (tj + 1) & (TILE_SIZE - 1);
                    }
                    data.force *= ENERGY_SCALE_FACTOR;
                    localData[threadIdx.x].force *= ENERGY_SCALE_FACTOR;
                    if (pos < end) {
                        unsigned int offset = x*TILE_SIZE + tgx;
                        atomicAdd(&forceBuffers[offset], static_cast<unsigned long long>((long long) (data.force.x*0x100000000)));
                        atomicAdd(&forceBuffers[offset+PADDED_NUM_ATOMS], static_cast<unsigned long long>((long long) (data.force.y*0x100000000)));
                        atomicAdd(&forceBuffers[offset+2*PADDED_NUM_ATOMS], static_cast<unsigned long long>((long long) (data.force.z*0x100000000)));
                        offset = y*TILE_SIZE + tgx;
                        atomicAdd(&forceBuffers[offset], static_cast<unsigned long long>((long long) (localData[threadIdx.x].force.x*0x100000000)));
                        atomicAdd(&forceBuffers[offset+PADDED_NUM_ATOMS], static_cast<unsigned long long>((long long) (localData[threadIdx.x].force.y*0x100000000)));
                        atomicAdd(&forceBuffers[offset+2*PADDED_NUM_ATOMS], static_cast<unsigned long long>((long long) (localData[threadIdx.x].force.z*0x100000000)));
                    }
                    
                    // Compute torques.
                    
                    data.force = make_real3(0);
                    localData[threadIdx.x].force = make_real3(0);
                    for (j = 0; j < TILE_SIZE; j++) {
                        int atom2 = y*TILE_SIZE+tj;
                        if (atom1 < NUM_ATOMS && atom2 < NUM_ATOMS) {
                            real3 tempForce;
                            float d = computeDScaleFactor(polarizationGroup, tj);
                            float p = computePScaleFactor(covalent, polarizationGroup, tj);
                            float m = computeMScaleFactor(covalent, tj);
                            computeOneInteractionT1(data, localData[tbx+tj], d, p, m, tempForce);
                            data.force += tempForce;
                            computeOneInteractionT3(data, localData[tbx+tj], d, p, m, tempForce);
                            localData[tbx+tj].force += tempForce;
                        }
                        tj = (tj + 1) & (TILE_SIZE - 1);
                    }
                    data.force *= ENERGY_SCALE_FACTOR;
                    localData[threadIdx.x].force *= ENERGY_SCALE_FACTOR;
                    if (pos < end) {
                        unsigned int offset = x*TILE_SIZE + tgx;
                        atomicAdd(&torqueBuffers[offset], static_cast<unsigned long long>((long long) (data.force.x*0x100000000)));
                        atomicAdd(&torqueBuffers[offset+PADDED_NUM_ATOMS], static_cast<unsigned long long>((long long) (data.force.y*0x100000000)));
                        atomicAdd(&torqueBuffers[offset+2*PADDED_NUM_ATOMS], static_cast<unsigned long long>((long long) (data.force.z*0x100000000)));
                        offset = y*TILE_SIZE + tgx;
                        atomicAdd(&torqueBuffers[offset], static_cast<unsigned long long>((long long) (localData[threadIdx.x].force.x*0x100000000)));
                        atomicAdd(&torqueBuffers[offset+PADDED_NUM_ATOMS], static_cast<unsigned long long>((long long) (localData[threadIdx.x].force.y*0x100000000)));
                        atomicAdd(&torqueBuffers[offset+2*PADDED_NUM_ATOMS], static_cast<unsigned long long>((long long) (localData[threadIdx.x].force.z*0x100000000)));
                    }
                }
            }
        }
        pos++;
    } while (pos < end);
    energyBuffer[blockIdx.x*blockDim.x+threadIdx.x] += energy*ENERGY_SCALE_FACTOR;
}
