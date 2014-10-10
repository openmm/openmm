#define WARPS_PER_GROUP (THREAD_BLOCK_SIZE/TILE_SIZE)

typedef struct {
    real3 pos;
    real3 field, fieldPolar, inducedDipole, inducedDipolePolar;
#ifdef USE_GK
    real3 fieldS, fieldPolarS, inducedDipoleS, inducedDipolePolarS;
    real bornRadius;
#endif
    float thole, damp;
} AtomData;

#ifdef USE_GK
inline __device__ void loadAtomData(AtomData& data, int atom, const real4* __restrict__ posq, const real* __restrict__ inducedDipole,
        const real* __restrict__ inducedDipolePolar, const float2* __restrict__ dampingAndThole, const real* __restrict__ inducedDipoleS,
        const real* __restrict__ inducedDipolePolarS, const real* __restrict__ bornRadii) {
#else
inline __device__ void loadAtomData(AtomData& data, int atom, const real4* __restrict__ posq, const real* __restrict__ inducedDipole,
        const real* __restrict__ inducedDipolePolar, const float2* __restrict__ dampingAndThole) {
#endif
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
#ifdef USE_GK
    data.inducedDipoleS.x = inducedDipoleS[atom*3];
    data.inducedDipoleS.y = inducedDipoleS[atom*3+1];
    data.inducedDipoleS.z = inducedDipoleS[atom*3+2];
    data.inducedDipolePolarS.x = inducedDipolePolarS[atom*3];
    data.inducedDipolePolarS.y = inducedDipolePolarS[atom*3+1];
    data.inducedDipolePolarS.z = inducedDipolePolarS[atom*3+2];
    data.bornRadius = bornRadii[atom];
#endif
}

inline __device__ void zeroAtomData(AtomData& data) {
    data.field = make_real3(0);
    data.fieldPolar = make_real3(0);
#ifdef USE_GK
    data.fieldS = make_real3(0);
    data.fieldPolarS = make_real3(0);
#endif
}

#ifdef USE_EWALD
__device__ void computeOneInteraction(AtomData& atom1, AtomData& atom2, real3 deltaR, bool isSelfInteraction) {
    if (isSelfInteraction)
        return;
    real scale1, scale2;
    real r2 = dot(deltaR, deltaR);
    if (r2 < CUTOFF_SQUARED) {
        real rI = RSQRT(r2);
        real r = RECIP(rI);

        // calculate the error function damping terms

        real ralpha = EWALD_ALPHA*r;
        real exp2a = EXP(-(ralpha*ralpha));
#ifdef USE_DOUBLE_PRECISION
        const real erfcAlphaR = erfc(ralpha);
#else
        // This approximation for erfc is from Abramowitz and Stegun (1964) p. 299.  They cite the following as
        // the original source: C. Hastings, Jr., Approximations for Digital Computers (1955).  It has a maximum
        // error of 1.5e-7.

        const real t = RECIP(1.0f+0.3275911f*ralpha);
        const real erfcAlphaR = (0.254829592f+(-0.284496736f+(1.421413741f+(-1.453152027f+1.061405429f*t)*t)*t)*t)*t*exp2a;
#endif
        real bn0 = erfcAlphaR*rI;
        real alsq2 = 2*EWALD_ALPHA*EWALD_ALPHA;
        real alsq2n = RECIP(SQRT_PI*EWALD_ALPHA);
        alsq2n *= alsq2;
        real bn1 = (bn0+alsq2n*exp2a)*rI*rI;

        alsq2n *= alsq2;
        real bn2 = (3*bn1+alsq2n*exp2a)*rI*rI;

        // compute the error function scaled and unscaled terms

        real scale3 = 1;
        real scale5 = 1;
        real damp = atom1.damp*atom2.damp;
        real ratio = (r/damp);
        ratio = ratio*ratio*ratio;
        float pgamma = atom1.thole < atom2.thole ? atom1.thole : atom2.thole;
        damp = damp == 0 ? 0 : -pgamma*ratio;
        real expdamp = EXP(damp);
        scale3 = 1 - expdamp;
        scale5 = 1 - expdamp*(1-damp);
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
    atom1.field += scale1*atom2.inducedDipole + dDotDelta*deltaR;
    dDotDelta = scale2*dot(deltaR, atom2.inducedDipolePolar);
    atom1.fieldPolar += scale1*atom2.inducedDipolePolar + dDotDelta*deltaR;
    dDotDelta = scale2*dot(deltaR, atom1.inducedDipole);
    atom2.field += scale1*atom1.inducedDipole + dDotDelta*deltaR;
    dDotDelta = scale2*dot(deltaR, atom1.inducedDipolePolar);
    atom2.fieldPolar += scale1*atom1.inducedDipolePolar + dDotDelta*deltaR;
}
#elif defined USE_GK
__device__ void computeOneInteraction(AtomData& atom1, AtomData& atom2, real3 deltaR, bool isSelfInteraction) {
    real r2 = dot(deltaR, deltaR);
    real r = SQRT(r2);
    if (!isSelfInteraction) {
        real rI = RECIP(r);
        real r2I = rI*rI;
        real rr3 = -rI*r2I;
        real rr5 = -3*rr3*r2I;

        real dampProd = atom1.damp*atom2.damp;
        real ratio = (dampProd != 0 ? r/dampProd : 1);
        float pGamma = (atom1.thole > atom2.thole ? atom2.thole: atom1.thole);
        real damp = ratio*ratio*ratio*pGamma;
        real dampExp = (dampProd != 0 ? EXP(-damp) : 0); 

        rr3 *= 1-dampExp;
        rr5 *= 1-(1+damp)*dampExp;

        real dDotDelta = rr5*dot(deltaR, atom2.inducedDipole);
        atom1.field += rr3*atom2.inducedDipole + dDotDelta*deltaR;
        dDotDelta = rr5*dot(deltaR, atom2.inducedDipolePolar);
        atom1.fieldPolar += rr3*atom2.inducedDipolePolar + dDotDelta*deltaR;
        dDotDelta = rr5*dot(deltaR, atom1.inducedDipole);
        atom2.field += rr3*atom1.inducedDipole + dDotDelta*deltaR;
        dDotDelta = rr5*dot(deltaR, atom1.inducedDipolePolar);
        atom2.fieldPolar += rr3*atom1.inducedDipolePolar + dDotDelta*deltaR;
        dDotDelta = rr5*dot(deltaR, atom2.inducedDipoleS);
        atom1.fieldS += rr3*atom2.inducedDipoleS + dDotDelta*deltaR;
        dDotDelta = rr5*dot(deltaR, atom2.inducedDipolePolarS);
        atom1.fieldPolarS += rr3*atom2.inducedDipolePolarS + dDotDelta*deltaR;
        dDotDelta = rr5*dot(deltaR, atom1.inducedDipoleS);
        atom2.fieldS += rr3*atom1.inducedDipoleS + dDotDelta*deltaR;
        dDotDelta = rr5*dot(deltaR, atom1.inducedDipolePolarS);
        atom2.fieldPolarS += rr3*atom1.inducedDipolePolarS + dDotDelta*deltaR;
    }

    real rb2 = atom1.bornRadius*atom2.bornRadius;
    real expterm = EXP(-r2/(GK_C*rb2));
    real expc = expterm/GK_C; 
    real gf2 = RECIP(r2+rb2*expterm);
    real gf = SQRT(gf2);
    real gf3 = gf2*gf;
    real gf5 = gf3*gf2;
    real a10 = -gf3;
    real expc1 = 1 - expc;
    real a11 = expc1 * 3 * gf5;
    real3 gux = GK_FD*make_real3(a10+deltaR.x*deltaR.x*a11, deltaR.x*deltaR.y*a11, deltaR.x*deltaR.z*a11);
    real3 guy = make_real3(gux.y, GK_FD*(a10+deltaR.y*deltaR.y*a11), GK_FD*deltaR.y*deltaR.z*a11);
    real3 guz = make_real3(gux.z, guy.z, GK_FD*(a10+deltaR.z*deltaR.z*a11));
 
    atom1.fieldS += atom2.inducedDipoleS.x*gux+atom2.inducedDipoleS.y*guy+atom2.inducedDipoleS.z*guz;
    atom2.fieldS += atom1.inducedDipoleS.x*gux+atom1.inducedDipoleS.y*guy+atom1.inducedDipoleS.z*guz;
    atom1.fieldPolarS += atom2.inducedDipolePolarS.x*gux+atom2.inducedDipolePolarS.y*guy+atom2.inducedDipolePolarS.z*guz;
    atom2.fieldPolarS += atom1.inducedDipolePolarS.x*gux+atom1.inducedDipolePolarS.y*guy+atom1.inducedDipolePolarS.z*guz;
}
#else
__device__ void computeOneInteraction(AtomData& atom1, AtomData& atom2, real3 deltaR, bool isSelfInteraction) {
    if (isSelfInteraction)
        return;
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
    atom1.field += rr3*atom2.inducedDipole + dDotDelta*deltaR;
    dDotDelta = rr5*dot(deltaR, atom2.inducedDipolePolar);
    atom1.fieldPolar += rr3*atom2.inducedDipolePolar + dDotDelta*deltaR;
    dDotDelta = rr5*dot(deltaR, atom1.inducedDipole);
    atom2.field += rr3*atom1.inducedDipole + dDotDelta*deltaR;
    dDotDelta = rr5*dot(deltaR, atom1.inducedDipolePolar);
    atom2.fieldPolar += rr3*atom1.inducedDipolePolar + dDotDelta*deltaR;
}
#endif

/**
 * Compute the mutual induced field.
 */
extern "C" __global__ void computeInducedField(
        unsigned long long* __restrict__ field, unsigned long long* __restrict__ fieldPolar, const real4* __restrict__ posq, const ushort2* __restrict__ exclusionTiles, 
        const real* __restrict__ inducedDipole, const real* __restrict__ inducedDipolePolar, unsigned int startTileIndex, unsigned int numTileIndices,
#ifdef USE_CUTOFF
        const int* __restrict__ tiles, const unsigned int* __restrict__ interactionCount, real4 periodicBoxSize, real4 invPeriodicBoxSize, unsigned int maxTiles, const real4* __restrict__ blockCenter, const unsigned int* __restrict__ interactingAtoms,
#elif defined USE_GK
        unsigned long long* __restrict__ fieldS, unsigned long long* __restrict__ fieldPolarS, const real* __restrict__ inducedDipoleS,
        const real* __restrict__ inducedDipolePolarS, const real* __restrict__ bornRadii,
#endif
        const float2* __restrict__ dampingAndThole) {
    const unsigned int totalWarps = (blockDim.x*gridDim.x)/TILE_SIZE;
    const unsigned int warp = (blockIdx.x*blockDim.x+threadIdx.x)/TILE_SIZE;
    const unsigned int tgx = threadIdx.x & (TILE_SIZE-1);
    const unsigned int tbx = threadIdx.x - tgx;
    __shared__ AtomData localData[THREAD_BLOCK_SIZE];

    // First loop: process tiles that contain exclusions.
    
    const unsigned int firstExclusionTile = FIRST_EXCLUSION_TILE+warp*(LAST_EXCLUSION_TILE-FIRST_EXCLUSION_TILE)/totalWarps;
    const unsigned int lastExclusionTile = FIRST_EXCLUSION_TILE+(warp+1)*(LAST_EXCLUSION_TILE-FIRST_EXCLUSION_TILE)/totalWarps;
    for (int pos = firstExclusionTile; pos < lastExclusionTile; pos++) {
        const ushort2 tileIndices = exclusionTiles[pos];
        const unsigned int x = tileIndices.x;
        const unsigned int y = tileIndices.y;
        AtomData data;
        zeroAtomData(data);
        unsigned int atom1 = x*TILE_SIZE + tgx;
#ifdef USE_GK
        loadAtomData(data, atom1, posq, inducedDipole, inducedDipolePolar, dampingAndThole, inducedDipoleS, inducedDipolePolarS, bornRadii);
#else
        loadAtomData(data, atom1, posq, inducedDipole, inducedDipolePolar, dampingAndThole);
#endif
        if (x == y) {
            // This tile is on the diagonal.

            localData[threadIdx.x].pos = data.pos;
            localData[threadIdx.x].inducedDipole = data.inducedDipole;
            localData[threadIdx.x].inducedDipolePolar = data.inducedDipolePolar;
            localData[threadIdx.x].thole = data.thole;
            localData[threadIdx.x].damp = data.damp;
#ifdef USE_GK
            localData[threadIdx.x].inducedDipoleS = data.inducedDipoleS;
            localData[threadIdx.x].inducedDipolePolarS = data.inducedDipolePolarS;
            localData[threadIdx.x].bornRadius = data.bornRadius;
#endif
            for (unsigned int j = 0; j < TILE_SIZE; j++) {
                real3 delta = localData[tbx+j].pos-data.pos;
#ifdef USE_PERIODIC
                delta.x -= floor(delta.x*invPeriodicBoxSize.x+0.5f)*periodicBoxSize.x;
                delta.y -= floor(delta.y*invPeriodicBoxSize.y+0.5f)*periodicBoxSize.y;
                delta.z -= floor(delta.z*invPeriodicBoxSize.z+0.5f)*periodicBoxSize.z;
#endif
                int atom2 = y*TILE_SIZE+j;
                if (atom1 < NUM_ATOMS && atom2 < NUM_ATOMS)
                    computeOneInteraction(data, localData[tbx+j], delta, atom1 == atom2);
            }
        }
        else {
            // This is an off-diagonal tile.

#ifdef USE_GK
            loadAtomData(localData[threadIdx.x], y*TILE_SIZE+tgx, posq, inducedDipole, inducedDipolePolar, dampingAndThole, inducedDipoleS, inducedDipolePolarS, bornRadii);
#else
            loadAtomData(localData[threadIdx.x], y*TILE_SIZE+tgx, posq, inducedDipole, inducedDipolePolar, dampingAndThole);
#endif
            zeroAtomData(localData[threadIdx.x]);
            unsigned int tj = tgx;
            for (unsigned int j = 0; j < TILE_SIZE; j++) {
                real3 delta = localData[tbx+tj].pos-data.pos;
#ifdef USE_PERIODIC
                delta.x -= floor(delta.x*invPeriodicBoxSize.x+0.5f)*periodicBoxSize.x;
                delta.y -= floor(delta.y*invPeriodicBoxSize.y+0.5f)*periodicBoxSize.y;
                delta.z -= floor(delta.z*invPeriodicBoxSize.z+0.5f)*periodicBoxSize.z;
#endif
                int atom2 = y*TILE_SIZE+j;
                if (atom1 < NUM_ATOMS && atom2 < NUM_ATOMS)
                    computeOneInteraction(data, localData[tbx+tj], delta, false);
                tj = (tj + 1) & (TILE_SIZE - 1);
            }
        }

        // Write results.

        unsigned int offset = x*TILE_SIZE + tgx;
        atomicAdd(&field[offset], static_cast<unsigned long long>((long long) (data.field.x*0x100000000)));
        atomicAdd(&field[offset+PADDED_NUM_ATOMS], static_cast<unsigned long long>((long long) (data.field.y*0x100000000)));
        atomicAdd(&field[offset+2*PADDED_NUM_ATOMS], static_cast<unsigned long long>((long long) (data.field.z*0x100000000)));
        atomicAdd(&fieldPolar[offset], static_cast<unsigned long long>((long long) (data.fieldPolar.x*0x100000000)));
        atomicAdd(&fieldPolar[offset+PADDED_NUM_ATOMS], static_cast<unsigned long long>((long long) (data.fieldPolar.y*0x100000000)));
        atomicAdd(&fieldPolar[offset+2*PADDED_NUM_ATOMS], static_cast<unsigned long long>((long long) (data.fieldPolar.z*0x100000000)));
#ifdef USE_GK
        atomicAdd(&fieldS[offset], static_cast<unsigned long long>((long long) (data.fieldS.x*0x100000000)));
        atomicAdd(&fieldS[offset+PADDED_NUM_ATOMS], static_cast<unsigned long long>((long long) (data.fieldS.y*0x100000000)));
        atomicAdd(&fieldS[offset+2*PADDED_NUM_ATOMS], static_cast<unsigned long long>((long long) (data.fieldS.z*0x100000000)));
        atomicAdd(&fieldPolarS[offset], static_cast<unsigned long long>((long long) (data.fieldPolarS.x*0x100000000)));
        atomicAdd(&fieldPolarS[offset+PADDED_NUM_ATOMS], static_cast<unsigned long long>((long long) (data.fieldPolarS.y*0x100000000)));
        atomicAdd(&fieldPolarS[offset+2*PADDED_NUM_ATOMS], static_cast<unsigned long long>((long long) (data.fieldPolarS.z*0x100000000)));
#endif
        if (x != y) {
            offset = y*TILE_SIZE + tgx;
            atomicAdd(&field[offset], static_cast<unsigned long long>((long long) (localData[threadIdx.x].field.x*0x100000000)));
            atomicAdd(&field[offset+PADDED_NUM_ATOMS], static_cast<unsigned long long>((long long) (localData[threadIdx.x].field.y*0x100000000)));
            atomicAdd(&field[offset+2*PADDED_NUM_ATOMS], static_cast<unsigned long long>((long long) (localData[threadIdx.x].field.z*0x100000000)));
            atomicAdd(&fieldPolar[offset], static_cast<unsigned long long>((long long) (localData[threadIdx.x].fieldPolar.x*0x100000000)));
            atomicAdd(&fieldPolar[offset+PADDED_NUM_ATOMS], static_cast<unsigned long long>((long long) (localData[threadIdx.x].fieldPolar.y*0x100000000)));
            atomicAdd(&fieldPolar[offset+2*PADDED_NUM_ATOMS], static_cast<unsigned long long>((long long) (localData[threadIdx.x].fieldPolar.z*0x100000000)));
#ifdef USE_GK
            atomicAdd(&fieldS[offset], static_cast<unsigned long long>((long long) (localData[threadIdx.x].fieldS.x*0x100000000)));
            atomicAdd(&fieldS[offset+PADDED_NUM_ATOMS], static_cast<unsigned long long>((long long) (localData[threadIdx.x].fieldS.y*0x100000000)));
            atomicAdd(&fieldS[offset+2*PADDED_NUM_ATOMS], static_cast<unsigned long long>((long long) (localData[threadIdx.x].fieldS.z*0x100000000)));
            atomicAdd(&fieldPolarS[offset], static_cast<unsigned long long>((long long) (localData[threadIdx.x].fieldPolarS.x*0x100000000)));
            atomicAdd(&fieldPolarS[offset+PADDED_NUM_ATOMS], static_cast<unsigned long long>((long long) (localData[threadIdx.x].fieldPolarS.y*0x100000000)));
            atomicAdd(&fieldPolarS[offset+2*PADDED_NUM_ATOMS], static_cast<unsigned long long>((long long) (localData[threadIdx.x].fieldPolarS.z*0x100000000)));
#endif
        }
    }

    // Second loop: tiles without exclusions, either from the neighbor list (with cutoff) or just enumerating all
    // of them (no cutoff).

#ifdef USE_CUTOFF
    const unsigned int numTiles = interactionCount[0];
    int pos = (int) (numTiles > maxTiles ? startTileIndex+warp*(long long)numTileIndices/totalWarps : warp*(long long)numTiles/totalWarps);
    int end = (int) (numTiles > maxTiles ? startTileIndex+(warp+1)*(long long)numTileIndices/totalWarps : (warp+1)*(long long)numTiles/totalWarps);
#else
    const unsigned int numTiles = numTileIndices;
    int pos = (int) (startTileIndex+warp*(long long)numTiles/totalWarps);
    int end = (int) (startTileIndex+(warp+1)*(long long)numTiles/totalWarps);
#endif
    int skipBase = 0;
    int currentSkipIndex = tbx;
    __shared__ int atomIndices[THREAD_BLOCK_SIZE];
    __shared__ volatile int skipTiles[THREAD_BLOCK_SIZE];
    skipTiles[threadIdx.x] = -1;
    
    while (pos < end) {
        bool includeTile = true;

        // Extract the coordinates of this tile.
        
        int x, y;
#ifdef USE_CUTOFF
        if (numTiles <= maxTiles)
            x = tiles[pos];
        else
#endif
        {
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
        }
        if (includeTile) {
            unsigned int atom1 = x*TILE_SIZE + tgx;

            // Load atom data for this tile.

            AtomData data;
            zeroAtomData(data);
#ifdef USE_GK
            loadAtomData(data, atom1, posq, inducedDipole, inducedDipolePolar, dampingAndThole, inducedDipoleS, inducedDipolePolarS, bornRadii);
#else
            loadAtomData(data, atom1, posq, inducedDipole, inducedDipolePolar, dampingAndThole);
#endif
#ifdef USE_CUTOFF
            unsigned int j = (numTiles <= maxTiles ? interactingAtoms[pos*TILE_SIZE+tgx] : y*TILE_SIZE + tgx);
#else
            unsigned int j = y*TILE_SIZE + tgx;
#endif
            atomIndices[threadIdx.x] = j;
#ifdef USE_GK
            loadAtomData(localData[threadIdx.x], j, posq, inducedDipole, inducedDipolePolar, dampingAndThole, inducedDipoleS, inducedDipolePolarS, bornRadii);
#else
            loadAtomData(localData[threadIdx.x], j, posq, inducedDipole, inducedDipolePolar, dampingAndThole);
#endif
            zeroAtomData(localData[threadIdx.x]);

            // Compute the full set of interactions in this tile.

            unsigned int tj = tgx;
            for (j = 0; j < TILE_SIZE; j++) {
                real3 delta = localData[tbx+tj].pos-data.pos;
#ifdef USE_PERIODIC
                delta.x -= floor(delta.x*invPeriodicBoxSize.x+0.5f)*periodicBoxSize.x;
                delta.y -= floor(delta.y*invPeriodicBoxSize.y+0.5f)*periodicBoxSize.y;
                delta.z -= floor(delta.z*invPeriodicBoxSize.z+0.5f)*periodicBoxSize.z;
#endif
                int atom2 = atomIndices[tbx+tj];
                if (atom1 < NUM_ATOMS && atom2 < NUM_ATOMS)
                    computeOneInteraction(data, localData[tbx+tj], delta, false);
                tj = (tj + 1) & (TILE_SIZE - 1);
            }

            // Write results.

            unsigned int offset = x*TILE_SIZE + tgx;
            atomicAdd(&field[offset], static_cast<unsigned long long>((long long) (data.field.x*0x100000000)));
            atomicAdd(&field[offset+PADDED_NUM_ATOMS], static_cast<unsigned long long>((long long) (data.field.y*0x100000000)));
            atomicAdd(&field[offset+2*PADDED_NUM_ATOMS], static_cast<unsigned long long>((long long) (data.field.z*0x100000000)));
            atomicAdd(&fieldPolar[offset], static_cast<unsigned long long>((long long) (data.fieldPolar.x*0x100000000)));
            atomicAdd(&fieldPolar[offset+PADDED_NUM_ATOMS], static_cast<unsigned long long>((long long) (data.fieldPolar.y*0x100000000)));
            atomicAdd(&fieldPolar[offset+2*PADDED_NUM_ATOMS], static_cast<unsigned long long>((long long) (data.fieldPolar.z*0x100000000)));
#ifdef USE_GK
            atomicAdd(&fieldS[offset], static_cast<unsigned long long>((long long) (data.fieldS.x*0x100000000)));
            atomicAdd(&fieldS[offset+PADDED_NUM_ATOMS], static_cast<unsigned long long>((long long) (data.fieldS.y*0x100000000)));
            atomicAdd(&fieldS[offset+2*PADDED_NUM_ATOMS], static_cast<unsigned long long>((long long) (data.fieldS.z*0x100000000)));
            atomicAdd(&fieldPolarS[offset], static_cast<unsigned long long>((long long) (data.fieldPolarS.x*0x100000000)));
            atomicAdd(&fieldPolarS[offset+PADDED_NUM_ATOMS], static_cast<unsigned long long>((long long) (data.fieldPolarS.y*0x100000000)));
            atomicAdd(&fieldPolarS[offset+2*PADDED_NUM_ATOMS], static_cast<unsigned long long>((long long) (data.fieldPolarS.z*0x100000000)));
#endif
#ifdef USE_CUTOFF
            offset = atomIndices[threadIdx.x];
#else
            offset = y*TILE_SIZE + tgx;
#endif
            atomicAdd(&field[offset], static_cast<unsigned long long>((long long) (localData[threadIdx.x].field.x*0x100000000)));
            atomicAdd(&field[offset+PADDED_NUM_ATOMS], static_cast<unsigned long long>((long long) (localData[threadIdx.x].field.y*0x100000000)));
            atomicAdd(&field[offset+2*PADDED_NUM_ATOMS], static_cast<unsigned long long>((long long) (localData[threadIdx.x].field.z*0x100000000)));
            atomicAdd(&fieldPolar[offset], static_cast<unsigned long long>((long long) (localData[threadIdx.x].fieldPolar.x*0x100000000)));
            atomicAdd(&fieldPolar[offset+PADDED_NUM_ATOMS], static_cast<unsigned long long>((long long) (localData[threadIdx.x].fieldPolar.y*0x100000000)));
            atomicAdd(&fieldPolar[offset+2*PADDED_NUM_ATOMS], static_cast<unsigned long long>((long long) (localData[threadIdx.x].fieldPolar.z*0x100000000)));
#ifdef USE_GK
            atomicAdd(&fieldS[offset], static_cast<unsigned long long>((long long) (localData[threadIdx.x].fieldS.x*0x100000000)));
            atomicAdd(&fieldS[offset+PADDED_NUM_ATOMS], static_cast<unsigned long long>((long long) (localData[threadIdx.x].fieldS.y*0x100000000)));
            atomicAdd(&fieldS[offset+2*PADDED_NUM_ATOMS], static_cast<unsigned long long>((long long) (localData[threadIdx.x].fieldS.z*0x100000000)));
            atomicAdd(&fieldPolarS[offset], static_cast<unsigned long long>((long long) (localData[threadIdx.x].fieldPolarS.x*0x100000000)));
            atomicAdd(&fieldPolarS[offset+PADDED_NUM_ATOMS], static_cast<unsigned long long>((long long) (localData[threadIdx.x].fieldPolarS.y*0x100000000)));
            atomicAdd(&fieldPolarS[offset+2*PADDED_NUM_ATOMS], static_cast<unsigned long long>((long long) (localData[threadIdx.x].fieldPolarS.z*0x100000000)));
#endif
        }
        pos++;
    }
}

extern "C" __global__ void updateInducedFieldBySOR(const long long* __restrict__ fixedField, const long long* __restrict__ fixedFieldPolar,
        const long long* __restrict__ fixedFieldS, const long long* __restrict__ inducedField, const long long* __restrict__ inducedFieldPolar,
        real* __restrict__ inducedDipole, real* __restrict__ inducedDipolePolar, const float* __restrict__ polarizability, float2* __restrict__ errors) {
    extern __shared__ real2 buffer[];
    const float polarSOR = 0.55f;
#ifdef USE_EWALD
    const real ewaldScale = (4/(real) 3)*(EWALD_ALPHA*EWALD_ALPHA*EWALD_ALPHA)/SQRT_PI;
#else
    const real ewaldScale = 0;
#endif
    const real fieldScale = 1/(real) 0x100000000;
    real sumErrors = 0;
    real sumPolarErrors = 0;
    for (int atom = blockIdx.x*blockDim.x + threadIdx.x; atom < NUM_ATOMS; atom += blockDim.x*gridDim.x) {
        real scale = polarizability[atom];
        for (int component = 0; component < 3; component++) {
            int dipoleIndex = 3*atom+component;
            int fieldIndex = atom+component*PADDED_NUM_ATOMS;
            real previousDipole = inducedDipole[dipoleIndex];
            real previousDipolePolar = inducedDipolePolar[dipoleIndex];
            long long fixedS = (fixedFieldS == NULL ? (long long) 0 : fixedFieldS[fieldIndex]);
            real newDipole = scale*((fixedField[fieldIndex]+fixedS+inducedField[fieldIndex])*fieldScale+ewaldScale*previousDipole);
            real newDipolePolar = scale*((fixedFieldPolar[fieldIndex]+fixedS+inducedFieldPolar[fieldIndex])*fieldScale+ewaldScale*previousDipolePolar);
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

extern "C" __global__ void recordInducedDipolesForDIIS(const long long* __restrict__ fixedField, const long long* __restrict__ fixedFieldPolar,
        const long long* __restrict__ fixedFieldS, const long long* __restrict__ inducedField, const long long* __restrict__ inducedFieldPolar,
        const real* __restrict__ inducedDipole, const real* __restrict__ inducedDipolePolar, const float* __restrict__ polarizability, float2* __restrict__ errors,
        real* __restrict__ prevDipoles, real* __restrict__ prevDipolesPolar, real* __restrict__ prevErrors, int iteration, bool recordPrevErrors, real* __restrict__ matrix) {
    extern __shared__ real2 buffer[];
#ifdef USE_EWALD
    const real ewaldScale = (4/(real) 3)*(EWALD_ALPHA*EWALD_ALPHA*EWALD_ALPHA)/SQRT_PI;
#else
    const real ewaldScale = 0;
#endif
    const real fieldScale = 1/(real) 0x100000000;
    real sumErrors = 0;
    real sumPolarErrors = 0;
    for (int atom = blockIdx.x*blockDim.x + threadIdx.x; atom < NUM_ATOMS; atom += blockDim.x*gridDim.x) {
        real scale = polarizability[atom];
        for (int component = 0; component < 3; component++) {
            int dipoleIndex = 3*atom+component;
            int fieldIndex = atom+component*PADDED_NUM_ATOMS;
            if (iteration >= MAX_PREV_DIIS_DIPOLES) {
                // We have filled up the buffer for previous dipoles, so shift them all over by one.
                
                for (int i = 1; i < MAX_PREV_DIIS_DIPOLES; i++) {
                    int index1 = dipoleIndex+(i-1)*NUM_ATOMS*3;
                    int index2 = dipoleIndex+i*NUM_ATOMS*3;
                    prevDipoles[index1] = prevDipoles[index2];
                    prevDipolesPolar[index1] = prevDipolesPolar[index2];
                    if (recordPrevErrors)
                        prevErrors[index1] = prevErrors[index2];
                }
            }
            
            // Compute the new dipole, and record it along with the error.
            
            real oldDipole = inducedDipole[dipoleIndex];
            real oldDipolePolar = inducedDipolePolar[dipoleIndex];
            long long fixedS = (fixedFieldS == NULL ? (long long) 0 : fixedFieldS[fieldIndex]);
            real newDipole = scale*((fixedField[fieldIndex]+fixedS+inducedField[fieldIndex])*fieldScale+ewaldScale*oldDipole);
            real newDipolePolar = scale*((fixedFieldPolar[fieldIndex]+fixedS+inducedFieldPolar[fieldIndex])*fieldScale+ewaldScale*oldDipolePolar);
            int storePrevIndex = dipoleIndex+min(iteration, MAX_PREV_DIIS_DIPOLES-1)*NUM_ATOMS*3;
            prevDipoles[storePrevIndex] = newDipole;
            prevDipolesPolar[storePrevIndex] = newDipolePolar;
            if (recordPrevErrors)
                prevErrors[storePrevIndex] = newDipole-oldDipole;
            sumErrors += (newDipole-oldDipole)*(newDipole-oldDipole);
            sumPolarErrors += (newDipolePolar-oldDipolePolar)*(newDipolePolar-oldDipolePolar);
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
    
    if (iteration >= MAX_PREV_DIIS_DIPOLES && recordPrevErrors && blockIdx.x == 0) {
        // Shift over the existing matrix elements.
        
        for (int i = 0; i < MAX_PREV_DIIS_DIPOLES-1; i++) {
            if (threadIdx.x < MAX_PREV_DIIS_DIPOLES-1)
                matrix[threadIdx.x+i*MAX_PREV_DIIS_DIPOLES] = matrix[(threadIdx.x+1)+(i+1)*MAX_PREV_DIIS_DIPOLES];
            __syncthreads();
        }
    }
}

extern "C" __global__ void computeDIISMatrix(real* __restrict__ prevErrors, int iteration, real* __restrict__ matrix) {
    extern __shared__ real sumBuffer[];
    int j = min(iteration, MAX_PREV_DIIS_DIPOLES-1);
    for (int i = blockIdx.x; i <= j; i += gridDim.x) {
        // All the threads in this thread block work together to compute a single matrix element.

        real sum = 0;
        for (int index = threadIdx.x; index < NUM_ATOMS*3; index += blockDim.x)
            sum += prevErrors[index+i*NUM_ATOMS*3]*prevErrors[index+j*NUM_ATOMS*3];
        sumBuffer[threadIdx.x] = sum;
        __syncthreads();
        for (int offset = 1; offset < blockDim.x; offset *= 2) { 
            if (threadIdx.x+offset < blockDim.x && (threadIdx.x&(2*offset-1)) == 0)
                sumBuffer[threadIdx.x] += sumBuffer[threadIdx.x+offset];
            __syncthreads();
        }
        if (threadIdx.x == 0) {
            matrix[i+MAX_PREV_DIIS_DIPOLES*j] = sumBuffer[0];
            if (i != j)
                matrix[j+MAX_PREV_DIIS_DIPOLES*i] = sumBuffer[0];
        }
    }
}

extern "C" __global__ void updateInducedFieldByDIIS(real* __restrict__ inducedDipole, real* __restrict__ inducedDipolePolar, 
        const real* __restrict__ prevDipoles, const real* __restrict__ prevDipolesPolar, const float* __restrict__ coefficients, int numPrev) {
    for (int index = blockIdx.x*blockDim.x + threadIdx.x; index < 3*NUM_ATOMS; index += blockDim.x*gridDim.x) {
        real sum = 0;
        real sumPolar = 0;
        for (int i = 0; i < numPrev; i++) {
            sum += coefficients[i]*prevDipoles[i*3*NUM_ATOMS+index];
            sumPolar += coefficients[i]*prevDipolesPolar[i*3*NUM_ATOMS+index];
        }
        inducedDipole[index] = sum;
        inducedDipolePolar[index] = sumPolar;
    }
}
