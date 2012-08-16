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
    data.quadrupoleXX = labFrameQuadrupole[atom*5];
    data.quadrupoleXY = labFrameQuadrupole[atom*5+1];
    data.quadrupoleXZ = labFrameQuadrupole[atom*5+2];
    data.quadrupoleYY = labFrameQuadrupole[atom*5+3];
    data.quadrupoleYZ = labFrameQuadrupole[atom*5+4];
    data.quadrupoleZZ = -(data.quadrupoleXX+data.quadrupoleYY);
    float2 temp = dampingAndThole[atom];
    data.damp = temp.x;
    data.thole = temp.y;
}

#ifdef USE_EWALD
__device__ void computeOneInteraction(AtomData& atom1, AtomData& atom2, real3 deltaR, float dScale, float pScale, real3* fields) {
    real r2 = dot(deltaR, deltaR);
    if (r2 <= CUTOFF_SQUARED) {
        // calculate the error function damping terms

        real r = SQRT(r2);
        real ralpha = EWALD_ALPHA*r;
        real bn0 = erfc(ralpha)/r;
        real alsq2 = 2*EWALD_ALPHA*EWALD_ALPHA;
        real alsq2n = RECIP(SQRT_PI*EWALD_ALPHA);
        real exp2a = EXP(-(ralpha*ralpha));
        alsq2n *= alsq2;
        real bn1 = (bn0+alsq2n*exp2a)/r2;
        alsq2n *= alsq2;
        real bn2 = (3*bn1+alsq2n*exp2a)/r2;
        alsq2n *= alsq2;
        real bn3 = (5*bn2+alsq2n*exp2a)/r2;

        // compute the error function scaled and unscaled terms

        real scale3 = 1;
        real scale5 = 1;
        real scale7 = 1;
        real damp = atom1.damp*atom2.damp;
        if (damp != 0) {
            real ratio = (r/damp);
            ratio = ratio*ratio*ratio;
            real pgamma = (atom1.thole < atom2.thole ? atom1.thole : atom2.thole);
            damp = -pgamma*ratio;
            if (damp > -50) {
                real expdamp = EXP(damp);
                scale3 = 1 - expdamp;
                scale5 = 1 - expdamp*(1-damp);
                scale7 = 1 - expdamp*(1-damp+(0.6f*damp*damp));
            }
        }
        real dsc3 = dScale*scale3;
        real dsc5 = dScale*scale5;
        real dsc7 = dScale*scale7;

        real psc3 = pScale*scale3;
        real psc5 = pScale*scale5;
        real psc7 = pScale*scale7;

        real r3 = r*r2;
        real r5 = r3*r2;
        real r7 = r5*r2;
        real drr3 = (1-dsc3)/r3;
        real drr5 = 3*(1-dsc5)/r5;
        real drr7 = 15*(1-dsc7)/r7;

        real prr3 = (1-psc3)/r3;
        real prr5 = 3*(1-psc5)/r5;
        real prr7 = 15*(1-psc7)/r7;

        real dir = dot(atom1.dipole, deltaR);

        real3 qi;
        qi.x = atom1.quadrupoleXX*deltaR.x + atom1.quadrupoleXY*deltaR.y + atom1.quadrupoleXZ*deltaR.z;
        qi.y = atom1.quadrupoleXY*deltaR.x + atom1.quadrupoleYY*deltaR.y + atom1.quadrupoleYZ*deltaR.z;
        qi.z = atom1.quadrupoleXZ*deltaR.x + atom1.quadrupoleYZ*deltaR.y + atom1.quadrupoleZZ*deltaR.z;
        real qir = dot(qi, deltaR);

        real dkr = dot(atom2.dipole, deltaR);

        real3 qk;
        qk.x = atom2.quadrupoleXX*deltaR.x + atom2.quadrupoleXY*deltaR.y + atom2.quadrupoleXZ*deltaR.z;
        qk.y = atom2.quadrupoleXY*deltaR.x + atom2.quadrupoleYY*deltaR.y + atom2.quadrupoleYZ*deltaR.z;
        qk.z = atom2.quadrupoleXZ*deltaR.x + atom2.quadrupoleYZ*deltaR.y + atom2.quadrupoleZZ*deltaR.z;
        real qkr = dot(qk, deltaR);

        real3 fim = -deltaR*(bn1*atom2.posq.w-bn2*dkr+bn3*qkr) - bn1*atom2.dipole + 2*bn2*qk;
        real3 fkm = deltaR*(bn1*atom1.posq.w+bn2*dir+bn3*qir) - bn1*atom1.dipole - 2*bn2*qi;
        real3 fid = -deltaR*(drr3*atom2.posq.w-drr5*dkr+drr7*qkr) - drr3*atom2.dipole + 2*drr5*qk;
        real3 fkd = deltaR*(drr3*atom1.posq.w+drr5*dir+drr7*qir) - drr3*atom1.dipole - 2*drr5*qi;
        real3 fip = -deltaR*(prr3*atom2.posq.w-prr5*dkr+prr7*qkr) - prr3*atom2.dipole + 2*prr5*qk;
        real3 fkp = deltaR*(prr3*atom1.posq.w+prr5*dir+prr7*qir) - prr3*atom1.dipole - 2*prr5*qi;

        // increment the field at each site due to this interaction

        fields[0] = fim-fid;
        fields[1] = fim-fip;
        fields[2] = fkm-fkd;
        fields[3] = fkm-fkp;
    }
    else {
        fields[0] = make_real3(0);
        fields[1] = make_real3(0);
        fields[2] = make_real3(0);
        fields[3] = make_real3(0);
    }
}
#else
__device__ void computeOneInteraction(AtomData& atom1, AtomData& atom2, real3 deltaR, float dScale, float pScale, real3* fields) {
    real rI = RSQRT(dot(deltaR, deltaR));
    real r = RECIP(rI);
    real r2I = rI*rI;

    real rr3 = rI*r2I;
    real rr5 = 3*rr3*r2I;
    real rr7 = 5*rr5*r2I;
 
    // get scaling factors, if needed
    
    float damp = atom1.damp*atom2.damp;
    real dampExp;
    if (damp != 0) {

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
 
    real3 field1 = deltaR*factor - rr3*atom2.dipole + rr5_2*qDotDelta;
    fields[0] = dScale*field1;
    fields[1] = pScale*field1;
 
    qDotDelta.x = deltaR.x*atom1.quadrupoleXX + deltaR.y*atom1.quadrupoleXY + deltaR.z*atom1.quadrupoleXZ;
    qDotDelta.y = deltaR.x*atom1.quadrupoleXY + deltaR.y*atom1.quadrupoleYY + deltaR.z*atom1.quadrupoleYZ;
    qDotDelta.z = deltaR.x*atom1.quadrupoleXZ + deltaR.y*atom1.quadrupoleYZ + deltaR.z*atom1.quadrupoleZZ;
 
    dotdd = dot(deltaR, atom1.dipole);
    dotqd = dot(deltaR, qDotDelta);
    factor = rr3*atom1.posq.w + rr5*dotdd + rr7*dotqd;
 
    real3 field2 = deltaR*factor - rr3*atom1.dipole - rr5_2*qDotDelta;
    fields[2] = dScale*field2;
    fields[3] = pScale*field2;
}
#endif

__device__ real computeDScaleFactor(unsigned int polarizationGroup) {
    return (polarizationGroup & 1 ? 0 : 1);
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
        const unsigned int* __restrict__ exclusionIndices, const unsigned int* __restrict__ exclusionRowIndices,
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
                uint2 covalent = covalentFlags[exclusionIndex[localGroupIndex]+tgx];
                unsigned int polarizationGroup = polarizationGroupFlags[exclusionIndex[localGroupIndex]+tgx];
                for (unsigned int j = 0; j < TILE_SIZE; j++) {
                    int atom2 = tbx+j;
                    real3 delta = make_real3(localData[atom2].posq.x-data.posq.x, localData[atom2].posq.y-data.posq.y, localData[atom2].posq.z-data.posq.z);
#ifdef USE_PERIODIC
                    delta.x -= floor(delta.x*invPeriodicBoxSize.x+0.5f)*periodicBoxSize.x;
                    delta.y -= floor(delta.y*invPeriodicBoxSize.y+0.5f)*periodicBoxSize.y;
                    delta.z -= floor(delta.z*invPeriodicBoxSize.z+0.5f)*periodicBoxSize.z;
#endif
                    real3 fields[4];
                    float d = computeDScaleFactor(polarizationGroup);
                    float p = computePScaleFactor(covalent, polarizationGroup);
                    computeOneInteraction(data, localData[atom2], delta, d, p, fields);
                    atom2 = y*TILE_SIZE+j;
                    if (atom1 != atom2 && atom1 < NUM_ATOMS && atom2 < NUM_ATOMS) {
                        data.field += fields[0];
                        data.fieldPolar += fields[1];
                    }
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
                                int atom2 = tbx+j;
                                real3 delta = make_real3(localData[atom2].posq.x-data.posq.x, localData[atom2].posq.y-data.posq.y, localData[atom2].posq.z-data.posq.z);
#ifdef USE_PERIODIC
                                delta.x -= floor(delta.x*invPeriodicBoxSize.x+0.5f)*periodicBoxSize.x;
                                delta.y -= floor(delta.y*invPeriodicBoxSize.y+0.5f)*periodicBoxSize.y;
                                delta.z -= floor(delta.z*invPeriodicBoxSize.z+0.5f)*periodicBoxSize.z;
#endif
                                real3 fields[4];
                                computeOneInteraction(data, localData[atom2], delta, 1, 1, fields);
                                if (atom1 < NUM_ATOMS && atom2 < NUM_ATOMS) {
#ifdef ENABLE_SHUFFLE
                                    for (int i = 16; i >= 1; i /= 2) {
                                        fields[2].x += __shfl_xor(fields[2].x, i, 32);
                                        fields[2].y += __shfl_xor(fields[2].y, i, 32);
                                        fields[2].z += __shfl_xor(fields[2].z, i, 32);
                                        fields[3].x += __shfl_xor(fields[3].x, i, 32);
                                        fields[3].y += __shfl_xor(fields[3].y, i, 32);
                                        fields[3].z += __shfl_xor(fields[3].z, i, 32);
                                    }
                                    if (tgx == 0) {
                                        localData[atom2].field += fields[2];
                                        localData[atom2].fieldPolar += fields[3];
                                    }
#else
                                    int bufferIndex = 3*threadIdx.x;
                                    tempBuffer[bufferIndex] = fields[2].x;
                                    tempBuffer[bufferIndex+1] = fields[2].y;
                                    tempBuffer[bufferIndex+2] = fields[2].z;
                                    if (tgx % 4 == 0) {
                                        tempBuffer[bufferIndex] += tempBuffer[bufferIndex+3]+tempBuffer[bufferIndex+6]+tempBuffer[bufferIndex+9];
                                        tempBuffer[bufferIndex+1] += tempBuffer[bufferIndex+4]+tempBuffer[bufferIndex+7]+tempBuffer[bufferIndex+10];
                                        tempBuffer[bufferIndex+2] += tempBuffer[bufferIndex+5]+tempBuffer[bufferIndex+8]+tempBuffer[bufferIndex+11];
                                    }
                                    if (tgx == 0) {
                                        localData[atom2].field.x += tempBuffer[bufferIndex]+tempBuffer[bufferIndex+12]+tempBuffer[bufferIndex+24]+tempBuffer[bufferIndex+36]+tempBuffer[bufferIndex+48]+tempBuffer[bufferIndex+60]+tempBuffer[bufferIndex+72]+tempBuffer[bufferIndex+84];
                                        localData[atom2].field.y += tempBuffer[bufferIndex+1]+tempBuffer[bufferIndex+13]+tempBuffer[bufferIndex+25]+tempBuffer[bufferIndex+37]+tempBuffer[bufferIndex+49]+tempBuffer[bufferIndex+61]+tempBuffer[bufferIndex+73]+tempBuffer[bufferIndex+85];
                                        localData[atom2].field.z += tempBuffer[bufferIndex+2]+tempBuffer[bufferIndex+14]+tempBuffer[bufferIndex+26]+tempBuffer[bufferIndex+38]+tempBuffer[bufferIndex+50]+tempBuffer[bufferIndex+62]+tempBuffer[bufferIndex+74]+tempBuffer[bufferIndex+86];
                                    }
                                    tempBuffer[bufferIndex] = fields[3].x;
                                    tempBuffer[bufferIndex+1] = fields[3].y;
                                    tempBuffer[bufferIndex+2] = fields[3].z;
                                    if (tgx % 4 == 0) {
                                        tempBuffer[bufferIndex] += tempBuffer[bufferIndex+3]+tempBuffer[bufferIndex+6]+tempBuffer[bufferIndex+9];
                                        tempBuffer[bufferIndex+1] += tempBuffer[bufferIndex+4]+tempBuffer[bufferIndex+7]+tempBuffer[bufferIndex+10];
                                        tempBuffer[bufferIndex+2] += tempBuffer[bufferIndex+5]+tempBuffer[bufferIndex+8]+tempBuffer[bufferIndex+11];
                                    }
                                    if (tgx == 0) {
                                        localData[atom2].fieldPolar.x += tempBuffer[bufferIndex]+tempBuffer[bufferIndex+12]+tempBuffer[bufferIndex+24]+tempBuffer[bufferIndex+36]+tempBuffer[bufferIndex+48]+tempBuffer[bufferIndex+60]+tempBuffer[bufferIndex+72]+tempBuffer[bufferIndex+84];
                                        localData[atom2].fieldPolar.y += tempBuffer[bufferIndex+1]+tempBuffer[bufferIndex+13]+tempBuffer[bufferIndex+25]+tempBuffer[bufferIndex+37]+tempBuffer[bufferIndex+49]+tempBuffer[bufferIndex+61]+tempBuffer[bufferIndex+73]+tempBuffer[bufferIndex+85];
                                        localData[atom2].fieldPolar.z += tempBuffer[bufferIndex+2]+tempBuffer[bufferIndex+14]+tempBuffer[bufferIndex+26]+tempBuffer[bufferIndex+38]+tempBuffer[bufferIndex+50]+tempBuffer[bufferIndex+62]+tempBuffer[bufferIndex+74]+tempBuffer[bufferIndex+86];
                                    }
#endif
                                }
                            }
                        }
                    }
                }
                else
#endif
                {
                    // Compute the full set of interactions in this tile.

                    uint2 covalent = (hasExclusions ? covalentFlags[exclusionIndex[localGroupIndex]+tgx] : make_uint2(0, 0));
                    unsigned int polarizationGroup = (hasExclusions ? polarizationGroupFlags[exclusionIndex[localGroupIndex]+tgx] : 0);
                    covalent.x = (covalent.x >> tgx) | (covalent.x << (TILE_SIZE - tgx));
                    covalent.y = (covalent.y >> tgx) | (covalent.y << (TILE_SIZE - tgx));
                    polarizationGroup = (polarizationGroup >> tgx) | (polarizationGroup << (TILE_SIZE - tgx));
                    unsigned int tj = tgx;
                    for (j = 0; j < TILE_SIZE; j++) {
                        int atom2 = tbx+tj;
                        real3 delta = make_real3(localData[atom2].posq.x-data.posq.x, localData[atom2].posq.y-data.posq.y, localData[atom2].posq.z-data.posq.z);
#ifdef USE_PERIODIC
                        delta.x -= floor(delta.x*invPeriodicBoxSize.x+0.5f)*periodicBoxSize.x;
                        delta.y -= floor(delta.y*invPeriodicBoxSize.y+0.5f)*periodicBoxSize.y;
                        delta.z -= floor(delta.z*invPeriodicBoxSize.z+0.5f)*periodicBoxSize.z;
#endif
                        real3 fields[4];
                        float d = computeDScaleFactor(polarizationGroup);
                        float p = computePScaleFactor(covalent, polarizationGroup);
                        computeOneInteraction(data, localData[atom2], delta, d, p, fields);
                        if (atom1 < NUM_ATOMS && atom2 < NUM_ATOMS) {
                            data.field += fields[0];
                            data.fieldPolar += fields[1];
                            localData[atom2].field += fields[2];
                            localData[atom2].fieldPolar += fields[3];
                        }
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
