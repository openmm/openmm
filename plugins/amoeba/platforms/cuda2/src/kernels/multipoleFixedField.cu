#define TILE_SIZE 32
#define WARPS_PER_GROUP (THREAD_BLOCK_SIZE/TILE_SIZE)

typedef struct {
    real4 posq;
    real3 field, fieldPolar, dipole;
    real quadrupoleXX, quadrupoleXY, quadrupoleXZ;
    real quadrupoleYY, quadrupoleYZ, quadrupoleZZ;
    float thole, damp;
} AtomData;

inline __device__ void loadAtomData(AtomData& data, int atom, const real4* __restrict__ posq, const real* __restrict__ labFrameDipole, const real* __restrict__ labFrameQuadrupole, const float2* __restrict__ dampingAndThole) {
    data.posq = posq[atom];
    data.dipole.x = labFrameDipole[atom*3];
    data.dipole.y = labFrameDipole[atom*3+1];
    data.dipole.z = labFrameDipole[atom*3+2];
    data.quadrupoleXX = labFrameQuadrupole[atom*9];
    data.quadrupoleXY = labFrameQuadrupole[atom*9+1];
    data.quadrupoleXZ = labFrameQuadrupole[atom*9+2];
    data.quadrupoleYY = labFrameQuadrupole[atom*9+4];
    data.quadrupoleYZ = labFrameQuadrupole[atom*9+5];
    data.quadrupoleZZ = labFrameQuadrupole[atom*9+8];
    float2 temp = dampingAndThole[atom];
    data.damp = temp.x;
    data.thole = temp.y;
}

__device__ void computeOneInteraction(AtomData& atom1, AtomData& atom2, real3 deltaR, real3& field1, real3& field2) {
    real rI = RSQRT(dot(deltaR, deltaR));
    real r = RECIP(rI);
    real r2I = rI*rI;

    real rr3 = rI*r2I;
    real rr5 = 3*rr3*r2I;
    real rr7 = 5*rr5*r2I;
 
    // get scaling factors, if needed
    
    float damp = atom1.damp*atom2.damp;
    real dampExp;
    if (damp != 0 && r < SCALING_DISTANCE_CUTOFF) {

        // get scaling factors
      
        real ratio = r/damp;
        float pGamma = atom2.thole > atom1.thole ? atom1.thole : atom2.thole; 
        damp = ratio*ratio*ratio*pGamma;
        dampExp = EXP(-damp);
    }
    else
        dampExp = 0;
      
    rr3 *= 1 - dampExp;
    rr5 *= 1 - (1+damp)*dampExp;
    rr7 *= 1 - (1+damp+(0.6f*damp*damp))*dampExp;
      
    real rr5_2 = 2*rr5;
 
    real3 qDotDelta;
    qDotDelta.x = deltaR.x*atom2.quadrupoleXX + deltaR.y*atom2.quadrupoleXY + deltaR.z*atom2.quadrupoleXZ;
    qDotDelta.y = deltaR.x*atom2.quadrupoleXY + deltaR.y*atom2.quadrupoleYY + deltaR.z*atom2.quadrupoleYZ;
    qDotDelta.z = deltaR.x*atom2.quadrupoleXZ + deltaR.y*atom2.quadrupoleYZ + deltaR.z*atom2.quadrupoleZZ;
 
    real dotdd = dot(deltaR, atom2.dipole);
    real dotqd = dot(deltaR, qDotDelta);

    real factor = -rr3*atom2.posq.w + rr5*dotdd - rr7*dotqd;
 
    field1 = deltaR*factor - rr3*atom2.dipole + rr5_2*qDotDelta;
 
    qDotDelta.x = deltaR.x*atom1.quadrupoleXX + deltaR.y*atom1.quadrupoleXY + deltaR.z*atom1.quadrupoleXZ;
    qDotDelta.y = deltaR.x*atom1.quadrupoleXY + deltaR.y*atom1.quadrupoleYY + deltaR.z*atom1.quadrupoleYZ;
    qDotDelta.z = deltaR.x*atom1.quadrupoleXZ + deltaR.y*atom1.quadrupoleYZ + deltaR.z*atom1.quadrupoleZZ;
 
    dotdd = dot(deltaR, atom1.dipole);
    dotqd = dot(deltaR, qDotDelta);
    factor = rr3*atom1.posq.w + rr5*dotdd + rr7*dotqd;
 
    field2 = deltaR*factor - rr3*atom1.dipole - rr5_2*qDotDelta;
}

__device__ real computeDScaleFactor(unsigned int polarizationGroup) {
    return (polarizationGroup & 1 ? 0 : 1);
}

__device__ float computeDScaleFactor(uint2 covalent) {
//        int f1 = (value == 0 || value == 1 ? 1 : 0);
//        int f2 = (value == 0 || value == 2 ? 1 : 0);
    // 0 = 12 or 13: x and y: 0
    // 1 = 14: x: 0.4
    // 2 = 15: y: 0.8
    bool x = (covalent.x & 1);
    bool y = (covalent.y & 1);
    return (x ? (y ? 0.0f : 0.4f) : (y ? 0.8f : 1.0f));
}

__device__ float computePScaleFactor(uint2 covalent, unsigned int polarizationGroup) {
    bool x = (covalent.x & 1);
    bool y = (covalent.y & 1);
    bool p = (polarizationGroup & 1);
    return (x && y ? 0.0f : (x && p ? 0.5f : 1.0f));
}

/**
 * Compute nonbonded interactions.
 */
extern "C" __global__ void computeFixedField(
        unsigned long long* __restrict__ fieldBuffers, unsigned long long* __restrict__ fieldPolarBuffers, const real4* __restrict__ posq,
        const unsigned int* __restrict__ exclusions, const unsigned int* __restrict__ exclusionIndices, const unsigned int* __restrict__ exclusionRowIndices,
        const uint2* __restrict__ covalentFlags, const unsigned int* __restrict__ polarizationGroupFlags, unsigned int startTileIndex, unsigned int numTileIndices,
#ifdef USE_CUTOFF
        const ushort2* __restrict__ tiles, const unsigned int* __restrict__ interactionCount, real4 periodicBoxSize, real4 invPeriodicBoxSize, unsigned int maxTiles, const unsigned int* __restrict__ interactionFlags,
#endif
        const real* __restrict__ labFrameDipole, const real* __restrict__ labFrameQuadrupole, const float2* __restrict__ dampingAndThole) {
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
        data.field = make_real3(0);
        data.fieldPolar = make_real3(0);
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
            loadAtomData(data, atom1, posq, labFrameDipole, labFrameQuadrupole, dampingAndThole);
            
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

                const unsigned int localAtomIndex = threadIdx.x;
                localData[localAtomIndex].posq = data.posq;
                localData[localAtomIndex].dipole = data.dipole;
                localData[localAtomIndex].quadrupoleXX = data.quadrupoleXX;
                localData[localAtomIndex].quadrupoleXY = data.quadrupoleXY;
                localData[localAtomIndex].quadrupoleXZ = data.quadrupoleXZ;
                localData[localAtomIndex].quadrupoleYY = data.quadrupoleYY;
                localData[localAtomIndex].quadrupoleYZ = data.quadrupoleYZ;
                localData[localAtomIndex].quadrupoleZZ = data.quadrupoleZZ;
                localData[localAtomIndex].thole = data.thole; // IS THIS CORRECT?
                localData[localAtomIndex].damp = data.damp; // IS THIS CORRECT?
                unsigned int excl = exclusions[exclusionIndex[localGroupIndex]+tgx];
                uint2 covalent = covalentFlags[exclusionIndex[localGroupIndex]+tgx];
                unsigned int polarizationGroup = polarizationGroupFlags[exclusionIndex[localGroupIndex]+tgx];
                for (unsigned int j = 0; j < TILE_SIZE; j++) {
                    bool isExcluded = !(excl & 0x1);
                    int atom2 = tbx+j;
                    real3 delta = make_real3(localData[atom2].posq.x-data.posq.x, localData[atom2].posq.y-data.posq.y, localData[atom2].posq.z-data.posq.z);
#ifdef USE_PERIODIC
                    delta.x -= floor(delta.x*invPeriodicBoxSize.x+0.5f)*periodicBoxSize.x;
                    delta.y -= floor(delta.y*invPeriodicBoxSize.y+0.5f)*periodicBoxSize.y;
                    delta.z -= floor(delta.z*invPeriodicBoxSize.z+0.5f)*periodicBoxSize.z;
#endif
                    real3 field1;
                    real3 field2;
                    computeOneInteraction(data, localData[atom2], delta, field1, field2);
                    if (!isExcluded) {
                        float d = computeDScaleFactor(covalent);
                        data.field += d*field1;
                        float p = computePScaleFactor(covalent, polarizationGroup);
                        data.fieldPolar += p*field1;
                    }
                    excl >>= 1;
                    covalent.x >>= 1;
                    covalent.y >>= 1;
                    polarizationGroup >>= 1;
                }
            }
            else {
                // This is an off-diagonal tile.

                const unsigned int localAtomIndex = threadIdx.x;
                unsigned int j = y*TILE_SIZE + tgx;
                loadAtomData(localData[localAtomIndex], j, posq, labFrameDipole, labFrameQuadrupole, dampingAndThole);
                localData[localAtomIndex].field = make_real3(0);
                localData[localAtomIndex].fieldPolar = make_real3(0);
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
                                real3 dEdR1 = make_real3(0);
                                real3 dEdR2 = make_real3(0);
                                real3 delta = make_real3(localData[atom2].posq.x-data.posq.x, localData[atom2].posq.y-data.posq.y, localData[atom2].posq.z-data.posq.z);
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
#ifdef USE_CUTOFF
                                }
#endif
#ifdef ENABLE_SHUFFLE
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
#else
                                force.x -= dEdR1.x;
                                force.y -= dEdR1.y;
                                force.z -= dEdR1.z;
                                tempBuffer[bufferIndex] = dEdR2.x;
                                tempBuffer[bufferIndex+1] = dEdR2.y;
                                tempBuffer[bufferIndex+2] = dEdR2.z;

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

                    unsigned int excl = (hasExclusions ? exclusions[exclusionIndex[localGroupIndex]+tgx] : 0xFFFFFFFF);
                    uint2 covalent = (hasExclusions ? covalentFlags[exclusionIndex[localGroupIndex]+tgx] : make_uint2(0, 0));
                    unsigned int polarizationGroup = (hasExclusions ? polarizationGroupFlags[exclusionIndex[localGroupIndex]+tgx] : 0);
                    excl = (excl >> tgx) | (excl << (TILE_SIZE - tgx));
                    covalent.x = (covalent.x >> tgx) | (covalent.x << (TILE_SIZE - tgx));
                    covalent.y = (covalent.y >> tgx) | (covalent.y << (TILE_SIZE - tgx));
                    polarizationGroup = (polarizationGroup >> tgx) | (polarizationGroup << (TILE_SIZE - tgx));
                    unsigned int tj = tgx;
                    for (j = 0; j < TILE_SIZE; j++) {
                        bool isExcluded = !(excl & 0x1);
                        int atom2 = tbx+tj;
                        real3 delta = make_real3(localData[atom2].posq.x-data.posq.x, localData[atom2].posq.y-data.posq.y, localData[atom2].posq.z-data.posq.z);
#ifdef USE_PERIODIC
                        delta.x -= floor(delta.x*invPeriodicBoxSize.x+0.5f)*periodicBoxSize.x;
                        delta.y -= floor(delta.y*invPeriodicBoxSize.y+0.5f)*periodicBoxSize.y;
                        delta.z -= floor(delta.z*invPeriodicBoxSize.z+0.5f)*periodicBoxSize.z;
#endif
                        real3 field1;
                        real3 field2;
                        computeOneInteraction(data, localData[atom2], delta, field1, field2);
                        if (!isExcluded) {
                            float d = computeDScaleFactor(covalent);
                            data.field += d*field1;
                            localData[atom2].field += d*field2;
                            float p = computePScaleFactor(covalent, polarizationGroup);
                            data.fieldPolar += p*field1;
                            localData[atom2].fieldPolar += p*field2;
                        }
                        excl >>= 1;
                        covalent.x >>= 1;
                        covalent.y >>= 1;
                        polarizationGroup >>= 1;
                        tj = (tj + 1) & (TILE_SIZE - 1);
                    }
                }
            }
        }
        
        // Write results.
        
        if (pos < end) {
            const unsigned int offset = x*TILE_SIZE + tgx;
            atomicAdd(&fieldBuffers[offset], static_cast<unsigned long long>((long long) (data.field.x*0xFFFFFFFF)));
            atomicAdd(&fieldBuffers[offset+PADDED_NUM_ATOMS], static_cast<unsigned long long>((long long) (data.field.y*0xFFFFFFFF)));
            atomicAdd(&fieldBuffers[offset+2*PADDED_NUM_ATOMS], static_cast<unsigned long long>((long long) (data.field.z*0xFFFFFFFF)));
            atomicAdd(&fieldPolarBuffers[offset], static_cast<unsigned long long>((long long) (data.fieldPolar.x*0xFFFFFFFF)));
            atomicAdd(&fieldPolarBuffers[offset+PADDED_NUM_ATOMS], static_cast<unsigned long long>((long long) (data.fieldPolar.y*0xFFFFFFFF)));
            atomicAdd(&fieldPolarBuffers[offset+2*PADDED_NUM_ATOMS], static_cast<unsigned long long>((long long) (data.fieldPolar.z*0xFFFFFFFF)));
        }
        if (pos < end && x != y) {
            const unsigned int offset = y*TILE_SIZE + tgx;
            atomicAdd(&fieldBuffers[offset], static_cast<unsigned long long>((long long) (localData[threadIdx.x].field.x*0xFFFFFFFF)));
            atomicAdd(&fieldBuffers[offset+PADDED_NUM_ATOMS], static_cast<unsigned long long>((long long) (localData[threadIdx.x].field.y*0xFFFFFFFF)));
            atomicAdd(&fieldBuffers[offset+2*PADDED_NUM_ATOMS], static_cast<unsigned long long>((long long) (localData[threadIdx.x].field.z*0xFFFFFFFF)));
            atomicAdd(&fieldPolarBuffers[offset], static_cast<unsigned long long>((long long) (localData[threadIdx.x].fieldPolar.x*0xFFFFFFFF)));
            atomicAdd(&fieldPolarBuffers[offset+PADDED_NUM_ATOMS], static_cast<unsigned long long>((long long) (localData[threadIdx.x].fieldPolar.y*0xFFFFFFFF)));
            atomicAdd(&fieldPolarBuffers[offset+2*PADDED_NUM_ATOMS], static_cast<unsigned long long>((long long) (localData[threadIdx.x].fieldPolar.z*0xFFFFFFFF)));
        }
        pos++;
    } while (pos < end);
}
