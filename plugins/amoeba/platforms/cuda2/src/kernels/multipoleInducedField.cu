#define TILE_SIZE 32
#define WARPS_PER_GROUP (THREAD_BLOCK_SIZE/TILE_SIZE)

typedef struct {
    real3 pos;
    real3 field, fieldPolar, inducedDipole, inducedDipolePolar;
    float thole, damp;
} AtomData;

inline __device__ void loadAtomData(AtomData& data, int atom, const real4* __restrict__ posq, const real* __restrict__ inducedDipole, const real* __restrict__ inducedDipolePolar, const float2* __restrict__ dampingAndThole) {
    real4 atomPosq = posq[atom];
    data.pos = make_real3(atomPosq.x, atomPosq.y, atomPosq.z);
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

#ifdef USE_EWALD
__device__ void computeOneInteraction(AtomData& atom1, AtomData& atom2, real3 deltaR, real3* fields) {
    real scale1, scale2;
    real r2 = dot(deltaR, deltaR);
    if (r2 < CUTOFF_SQUARED) {
        real rI = RSQRT(r2);
        real r = RECIP(rI);

        // calculate the error function damping terms

        real ralpha = EWALD_ALPHA*r;
        real bn0 = erfc(ralpha)*rI;
        real alsq2 = 2*EWALD_ALPHA*EWALD_ALPHA;
        real alsq2n = RECIP(SQRT(M_PI)*EWALD_ALPHA);
        real exp2a = expf(-(ralpha*ralpha));
        alsq2n *= alsq2;
        real bn1 = (bn0+alsq2n*exp2a)*rI*rI;

        alsq2n *= alsq2;
        real bn2 = (3*bn1+alsq2n*exp2a)*rI*rI;

        // compute the error function scaled and unscaled terms

        real scale3 = 1;
        real scale5 = 1;
        real damp = atom1.damp*atom2.damp;
        if (damp != 0) {
            real ratio = (r/damp);
            ratio = ratio*ratio*ratio;
            float pgamma = atom1.thole < atom2.thole ? atom1.thole : atom2.thole;
            damp = -pgamma*ratio;
            if (damp > -50) {
                real expdamp = EXP(damp);
                scale3 = 1 - expdamp;
                scale5 = 1 - expdamp*(1-damp);
            }
        }
        real dsc3 = scale3;
        real dsc5 = scale5;
        real r3 = (r*r2);
        real r5 = (r3*r2);
        real rr3 = (1-dsc3)/r3;
        real rr5 = 3*(1-dsc5)/r5;

        scale1 = rr3 - bn1;
        scale2 = bn2 - rr5;
    }
    else {
        scale1 = 0;
        scale2 = 0;
    }
    real dDotDelta = scale2*dot(deltaR, atom2.inducedDipole);
    fields[0] = scale1*atom2.inducedDipole + dDotDelta*deltaR;
    dDotDelta = scale2*dot(deltaR, atom2.inducedDipolePolar);
    fields[1] = scale1*atom2.inducedDipolePolar + dDotDelta*deltaR;
    dDotDelta = scale2*dot(deltaR, atom1.inducedDipole);
    fields[2] = scale1*atom1.inducedDipole + dDotDelta*deltaR;
    dDotDelta = scale2*dot(deltaR, atom1.inducedDipolePolar);
    fields[3] = scale1*atom1.inducedDipolePolar + dDotDelta*deltaR;
}
#else
__device__ void computeOneInteraction(AtomData& atom1, AtomData& atom2, real3 deltaR, real3* fields) {
    real rI = RSQRT(dot(deltaR, deltaR));
    real r = RECIP(rI);
    real r2I = rI*rI;
    real rr3 = -rI*r2I;
    real rr5 = -3*rr3*r2I;
    real dampProd = atom1.damp*atom2.damp;
    real ratio = (dampProd != 0 ? r/dampProd : 1);
    float pGamma = (atom2.thole > atom1.thole ? atom1.thole: atom2.thole);
    real damp = ratio*ratio*ratio*pGamma;
    real dampExp = (dampProd != 0 ? EXP(-damp) : 0); 
    rr3 *= 1 - dampExp;
    rr5 *= 1 - (1+damp)*dampExp;
    real dDotDelta = rr5*dot(deltaR, atom2.inducedDipole);
    fields[0] = rr3*atom2.inducedDipole + dDotDelta*deltaR;
    dDotDelta = rr5*dot(deltaR, atom2.inducedDipolePolar);
    fields[1] = rr3*atom2.inducedDipolePolar + dDotDelta*deltaR;
    dDotDelta = rr5*dot(deltaR, atom1.inducedDipole);
    fields[2] = rr3*atom1.inducedDipole + dDotDelta*deltaR;
    dDotDelta = rr5*dot(deltaR, atom1.inducedDipolePolar);
    fields[3] = rr3*atom1.inducedDipolePolar + dDotDelta*deltaR;
}
#endif

/**
 * Compute the mutual induced field.
 */
extern "C" __global__ void computeInducedField(
        unsigned long long* __restrict__ fieldBuffers, unsigned long long* __restrict__ fieldPolarBuffers, const real4* __restrict__ posq,
        const real* __restrict__ inducedDipole, const real* __restrict__ inducedDipolePolar, unsigned int startTileIndex, unsigned int numTileIndices,
#ifdef USE_CUTOFF
        const ushort2* __restrict__ tiles, const unsigned int* __restrict__ interactionCount, real4 periodicBoxSize, real4 invPeriodicBoxSize, unsigned int maxTiles, const unsigned int* __restrict__ interactionFlags,
#endif
        const float2* __restrict__ dampingAndThole) {
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
#ifndef ENABLE_SHUFFLE
    __shared__ real tempBuffer[3*THREAD_BLOCK_SIZE];
#endif
    
    do {
        // Extract the coordinates of this tile
        const unsigned int tgx = threadIdx.x & (TILE_SIZE-1);
        const unsigned int tbx = threadIdx.x - tgx;
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
            loadAtomData(data, atom1, posq, inducedDipole, inducedDipolePolar, dampingAndThole);
            if (pos >= end)
                ; // This warp is done.
            else if (x == y) {
                // This tile is on the diagonal.

                localData[threadIdx.x].pos = data.pos;
                localData[threadIdx.x].inducedDipole = data.inducedDipole;
                localData[threadIdx.x].inducedDipolePolar = data.inducedDipolePolar;
                localData[threadIdx.x].thole = data.thole;
                localData[threadIdx.x].damp = data.damp;
                for (unsigned int j = 0; j < TILE_SIZE; j++) {
                    int atom2 = tbx+j;
                    real3 delta = localData[atom2].pos-data.pos;
#ifdef USE_PERIODIC
                    delta.x -= floor(delta.x*invPeriodicBoxSize.x+0.5f)*periodicBoxSize.x;
                    delta.y -= floor(delta.y*invPeriodicBoxSize.y+0.5f)*periodicBoxSize.y;
                    delta.z -= floor(delta.z*invPeriodicBoxSize.z+0.5f)*periodicBoxSize.z;
#endif
                    real3 fields[4];
                    computeOneInteraction(data, localData[atom2], delta, fields);
                    atom2 = y*TILE_SIZE+j;
                    if (atom1 != atom2 && atom1 < NUM_ATOMS && atom2 < NUM_ATOMS) {
                        data.field += fields[0];
                        data.fieldPolar += fields[1];
                    }
                }
            }
            else {
                // This is an off-diagonal tile.

                loadAtomData(localData[threadIdx.x], y*TILE_SIZE+tgx, posq, inducedDipole, inducedDipolePolar, dampingAndThole);
                localData[threadIdx.x].field = make_real3(0);
                localData[threadIdx.x].fieldPolar = make_real3(0);
#ifdef USE_CUTOFF
                unsigned int flags = (numTiles <= maxTiles ? interactionFlags[pos] : 0xFFFFFFFF);
                if (flags != 0xFFFFFFFF) {
                    if (flags == 0) {
                        // No interactions in this tile.
                    }
                    else {
                        // Compute only a subset of the interactions in this tile.

                        for (int j = 0; j < TILE_SIZE; j++) {
                            if ((flags&(1<<j)) != 0) {
                                int atom2 = tbx+j;
                                real3 delta = localData[atom2].pos-data.pos;
#ifdef USE_PERIODIC
                                delta.x -= floor(delta.x*invPeriodicBoxSize.x+0.5f)*periodicBoxSize.x;
                                delta.y -= floor(delta.y*invPeriodicBoxSize.y+0.5f)*periodicBoxSize.y;
                                delta.z -= floor(delta.z*invPeriodicBoxSize.z+0.5f)*periodicBoxSize.z;
#endif
                                real3 fields[4];
                                computeOneInteraction(data, localData[atom2], delta, fields);
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

                    unsigned int tj = tgx;
                    for (unsigned int j = 0; j < TILE_SIZE; j++) {
                        int atom2 = tbx+tj;
                        real3 delta = localData[atom2].pos-data.pos;
#ifdef USE_PERIODIC
                        delta.x -= floor(delta.x*invPeriodicBoxSize.x+0.5f)*periodicBoxSize.x;
                        delta.y -= floor(delta.y*invPeriodicBoxSize.y+0.5f)*periodicBoxSize.y;
                        delta.z -= floor(delta.z*invPeriodicBoxSize.z+0.5f)*periodicBoxSize.z;
#endif
                        real3 fields[4];
                        computeOneInteraction(data, localData[atom2], delta, fields);
                        if (atom1 < NUM_ATOMS && atom2 < NUM_ATOMS) {
                            data.field += fields[0];
                            data.fieldPolar += fields[1];
                            localData[atom2].field += fields[2];
                            localData[atom2].fieldPolar += fields[3];
                        }
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

extern "C" __global__ void updateInducedFieldBySOR(const long long* __restrict__ fixedField, const long long* __restrict__ fixedFieldPolar,
        const long long* __restrict__ inducedField, const long long* __restrict__ inducedFieldPolar, real* __restrict__ inducedDipole,
        real* __restrict__ inducedDipolePolar, const float* __restrict__ polarizability, float2* __restrict__ errors) {
    extern __shared__ real2 buffer[];
    const float polarSOR = 0.55f;
    const real term = (4/(real) 3)*(EWALD_ALPHA*EWALD_ALPHA*EWALD_ALPHA)/SQRT(M_PI);
    const real fieldScale = 1/(real) 0xFFFFFFFF;
    real sumErrors = 0;
    real sumPolarErrors = 0;
    for (int atom = blockIdx.x*blockDim.x + threadIdx.x; atom < NUM_ATOMS; atom += blockDim.x*gridDim.x) {
        real scale = polarizability[atom];
        for (int component = 0; component < 3; component++) {
            int dipoleIndex = 3*atom+component;
            int fieldIndex = atom+component*PADDED_NUM_ATOMS;
            real previousDipole = inducedDipole[dipoleIndex];
            real previousDipolePolar = inducedDipolePolar[dipoleIndex];
            real newDipole = scale*((fixedField[fieldIndex]+inducedField[fieldIndex])*fieldScale+term*previousDipole);
            real newDipolePolar = scale*((fixedFieldPolar[fieldIndex]+inducedFieldPolar[fieldIndex])*fieldScale+term*previousDipolePolar);
            newDipole = previousDipole + polarSOR*(newDipole-previousDipole);
            newDipolePolar = previousDipolePolar + polarSOR*(newDipolePolar-previousDipolePolar);
            inducedDipole[dipoleIndex] = newDipole;
            inducedDipolePolar[dipoleIndex] = newDipolePolar;
            sumErrors += (newDipole-previousDipole)*(newDipole-previousDipole);
            sumPolarErrors += (newDipolePolar-previousDipolePolar)*(newDipolePolar-previousDipolePolar);
        }
    }
    
    // Sum the errors over threads and store the total for this block.
    
    buffer[threadIdx.x] = make_real2(sumErrors, sumPolarErrors);
    __syncthreads();
    for (int offset = 1; offset < blockDim.x; offset *= 2) {   
        if (threadIdx.x+offset < blockDim.x && (threadIdx.x&(2*offset-1)) == 0) {
            buffer[threadIdx.x].x += buffer[threadIdx.x+offset].x;
            buffer[threadIdx.x].y += buffer[threadIdx.x+offset].y;
        }
        __syncthreads();
    }
    if (threadIdx.x == 0)
        errors[blockIdx.x] = make_float2((float) buffer[0].x, (float) buffer[0].y);
}