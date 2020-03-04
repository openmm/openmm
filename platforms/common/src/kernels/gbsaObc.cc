#define WARPS_PER_GROUP (FORCE_WORK_GROUP_SIZE/TILE_SIZE)

typedef struct {
    real x, y, z;
    real q;
    float radius, scaledRadius;
    real bornSum;
} AtomData1;

/**
 * Compute the Born sum.
 */
KERNEL void computeBornSum(
#ifdef SUPPORTS_64_BIT_ATOMICS
        GLOBAL mm_ulong* RESTRICT global_bornSum,
#else
        GLOBAL real* RESTRICT global_bornSum,
#endif
        GLOBAL const real4* RESTRICT posq, GLOBAL const real* RESTRICT charge, GLOBAL const float2* RESTRICT global_params,
#ifdef USE_CUTOFF
        GLOBAL const int* RESTRICT tiles, GLOBAL const unsigned int* RESTRICT interactionCount, real4 periodicBoxSize, real4 invPeriodicBoxSize,
        real4 periodicBoxVecX, real4 periodicBoxVecY, real4 periodicBoxVecZ, unsigned int maxTiles, GLOBAL const real4* RESTRICT blockCenter,
        GLOBAL const real4* RESTRICT blockSize, GLOBAL const int* RESTRICT interactingAtoms,
#else
        unsigned int numTiles,
#endif
        GLOBAL const ushort2* RESTRICT exclusionTiles) {
    const unsigned int totalWarps = GLOBAL_SIZE/TILE_SIZE;
    const unsigned int warp = GLOBAL_ID/TILE_SIZE;
    const unsigned int tgx = LOCAL_ID & (TILE_SIZE-1);
    const unsigned int tbx = LOCAL_ID - tgx;
    LOCAL AtomData1 localData[FORCE_WORK_GROUP_SIZE];

    // First loop: process tiles that contain exclusions.
    
    const unsigned int firstExclusionTile = FIRST_EXCLUSION_TILE+warp*(LAST_EXCLUSION_TILE-FIRST_EXCLUSION_TILE)/totalWarps;
    const unsigned int lastExclusionTile = FIRST_EXCLUSION_TILE+(warp+1)*(LAST_EXCLUSION_TILE-FIRST_EXCLUSION_TILE)/totalWarps;
    for (int pos = firstExclusionTile; pos < lastExclusionTile; pos++) {
        const ushort2 tileIndices = exclusionTiles[pos];
        const unsigned int x = tileIndices.x;
        const unsigned int y = tileIndices.y;
        real bornSum = 0;
        unsigned int atom1 = x*TILE_SIZE + tgx;
        real4 posq1 = posq[atom1];
        real charge1 = charge[atom1];
        float2 params1 = global_params[atom1];
        if (x == y) {
            // This tile is on the diagonal.

            localData[LOCAL_ID].x = posq1.x;
            localData[LOCAL_ID].y = posq1.y;
            localData[LOCAL_ID].z = posq1.z;
            localData[LOCAL_ID].q = charge1;
            localData[LOCAL_ID].radius = params1.x;
            localData[LOCAL_ID].scaledRadius = params1.y;
            SYNC_WARPS;
            for (unsigned int j = 0; j < TILE_SIZE; j++) {
                real3 delta = make_real3(localData[tbx+j].x-posq1.x, localData[tbx+j].y-posq1.y, localData[tbx+j].z-posq1.z);
#ifdef USE_PERIODIC
                APPLY_PERIODIC_TO_DELTA(delta)
#endif
                real r2 = delta.x*delta.x + delta.y*delta.y + delta.z*delta.z;
#ifdef USE_CUTOFF
                if (atom1 < NUM_ATOMS && y*TILE_SIZE+j < NUM_ATOMS && r2 < CUTOFF_SQUARED) {
#else
                if (atom1 < NUM_ATOMS && y*TILE_SIZE+j < NUM_ATOMS) {
#endif
                    real invR = RSQRT(r2);
                    real r = r2*invR;
                    float2 params2 = make_float2(localData[tbx+j].radius, localData[tbx+j].scaledRadius);
                    real rScaledRadiusJ = r+params2.y;
                    if ((j != tgx) && (params1.x < rScaledRadiusJ)) {
                        real l_ij = RECIP(max((real) params1.x, fabs(r-params2.y)));
                        real u_ij = RECIP(rScaledRadiusJ);
                        real l_ij2 = l_ij*l_ij;
                        real u_ij2 = u_ij*u_ij;
                        real ratio = LOG(u_ij * RECIP(l_ij));
                        bornSum += l_ij - u_ij + (0.50f*invR*ratio) + 0.25f*(r*(u_ij2-l_ij2) +
                                         (params2.y*params2.y*invR)*(l_ij2-u_ij2));
                        bornSum += (params1.x < params2.y-r ? 2.0f*(RECIP(params1.x)-l_ij) : 0);
                    }
                }
                SYNC_WARPS;
            }
        }
        else {
            // This is an off-diagonal tile.

            unsigned int j = y*TILE_SIZE + tgx;
            real4 tempPosq = posq[j];
            localData[LOCAL_ID].x = tempPosq.x;
            localData[LOCAL_ID].y = tempPosq.y;
            localData[LOCAL_ID].z = tempPosq.z;
            localData[LOCAL_ID].q = charge[j];
            float2 tempParams = global_params[j];
            localData[LOCAL_ID].radius = tempParams.x;
            localData[LOCAL_ID].scaledRadius = tempParams.y;
            localData[LOCAL_ID].bornSum = 0.0f;
            SYNC_WARPS;

            // Compute the full set of interactions in this tile.

            unsigned int tj = tgx;
            for (j = 0; j < TILE_SIZE; j++) {
                real3 delta = make_real3(localData[tbx+tj].x-posq1.x, localData[tbx+tj].y-posq1.y, localData[tbx+tj].z-posq1.z);
#ifdef USE_PERIODIC
                APPLY_PERIODIC_TO_DELTA(delta)
#endif
                real r2 = delta.x*delta.x + delta.y*delta.y + delta.z*delta.z;
#ifdef USE_CUTOFF
                if (atom1 < NUM_ATOMS && y*TILE_SIZE+tj < NUM_ATOMS && r2 < CUTOFF_SQUARED) {
#else
                if (atom1 < NUM_ATOMS && y*TILE_SIZE+tj < NUM_ATOMS) {
#endif
                    real invR = RSQRT(r2);
                    real r = r2*invR;
                    float2 params2 = make_float2(localData[tbx+tj].radius, localData[tbx+tj].scaledRadius);
                    real rScaledRadiusJ = r+params2.y;
                    if (params1.x < rScaledRadiusJ) {
                        real l_ij = RECIP(max((real) params1.x, fabs(r-params2.y)));
                        real u_ij = RECIP(rScaledRadiusJ);
                        real l_ij2 = l_ij*l_ij;
                        real u_ij2 = u_ij*u_ij;
                        real ratio = LOG(u_ij * RECIP(l_ij));
                        bornSum += l_ij - u_ij + (0.50f*invR*ratio) + 0.25f*(r*(u_ij2-l_ij2) +
                                         (params2.y*params2.y*invR)*(l_ij2-u_ij2));
                        bornSum += (params1.x < params2.y-r ? 2.0f*(RECIP(params1.x)-l_ij) : 0);
                    }
                    real rScaledRadiusI = r+params1.y;
                    if (params2.x < rScaledRadiusI) {
                        real l_ij = RECIP(max((real) params2.x, fabs(r-params1.y)));
                        real u_ij = RECIP(rScaledRadiusI);
                        real l_ij2 = l_ij*l_ij;
                        real u_ij2 = u_ij*u_ij;
                        real ratio = LOG(u_ij * RECIP(l_ij));
                        real term = l_ij - u_ij + (0.50f*invR*ratio) + 0.25f*(r*(u_ij2-l_ij2) +
                                         (params1.y*params1.y*invR)*(l_ij2-u_ij2));
                        term += (params2.x < params1.y-r ? 2.0f*(RECIP(params2.x)-l_ij) : 0);
                        localData[tbx+tj].bornSum += term;
                    }
                }
                tj = (tj + 1) & (TILE_SIZE - 1);
                SYNC_WARPS;
            }
        }

        // Write results.

#ifdef SUPPORTS_64_BIT_ATOMICS
        unsigned int offset = x*TILE_SIZE + tgx;
        ATOMIC_ADD(&global_bornSum[offset], (mm_ulong) ((mm_long) (bornSum*0x100000000)));
        if (x != y) {
            offset = y*TILE_SIZE + tgx;
            ATOMIC_ADD(&global_bornSum[offset], (mm_ulong) ((mm_long) (localData[LOCAL_ID].bornSum*0x100000000)));
        }
#else
        unsigned int offset1 = x*TILE_SIZE + tgx + warp*PADDED_NUM_ATOMS;
        unsigned int offset2 = y*TILE_SIZE + tgx + warp*PADDED_NUM_ATOMS;
        global_bornSum[offset1] += bornSum;
        if (x != y)
            global_bornSum[offset2] += localData[LOCAL_ID].bornSum;
#endif
    }

    // Second loop: tiles without exclusions, either from the neighbor list (with cutoff) or just enumerating all
    // of them (no cutoff).

#ifdef USE_CUTOFF
    unsigned int numTiles = interactionCount[0];
    if (numTiles > maxTiles)
        return; // There wasn't enough memory for the neighbor list.
    int pos = (int) (warp*(numTiles > maxTiles ? NUM_BLOCKS*((mm_long)NUM_BLOCKS+1)/2 : (mm_long)numTiles)/totalWarps);
    int end = (int) ((warp+1)*(numTiles > maxTiles ? NUM_BLOCKS*((mm_long)NUM_BLOCKS+1)/2 : (mm_long)numTiles)/totalWarps);
#else
    int pos = (int) (warp*(mm_long)numTiles/totalWarps);
    int end = (int) ((warp+1)*(mm_long)numTiles/totalWarps);
#endif
    int skipBase = 0;
    int currentSkipIndex = tbx;
    LOCAL int atomIndices[FORCE_WORK_GROUP_SIZE];
    LOCAL volatile int skipTiles[FORCE_WORK_GROUP_SIZE];
    skipTiles[LOCAL_ID] = -1;

    while (pos < end) {
        real bornSum = 0;
        bool includeTile = true;

        // Extract the coordinates of this tile.
        
        int x, y;
        bool singlePeriodicCopy = false;
#ifdef USE_CUTOFF
        x = tiles[pos];
        real4 blockSizeX = blockSize[x];
        singlePeriodicCopy = (0.5f*periodicBoxSize.x-blockSizeX.x >= CUTOFF &&
                              0.5f*periodicBoxSize.y-blockSizeX.y >= CUTOFF &&
                              0.5f*periodicBoxSize.z-blockSizeX.z >= CUTOFF);
#else
        y = (int) floor(NUM_BLOCKS+0.5f-SQRT((NUM_BLOCKS+0.5f)*(NUM_BLOCKS+0.5f)-2*pos));
        x = (pos-y*NUM_BLOCKS+y*(y+1)/2);
        if (x < y || x >= NUM_BLOCKS) { // Occasionally happens due to roundoff error.
            y += (x < y ? -1 : 1);
            x = (pos-y*NUM_BLOCKS+y*(y+1)/2);
        }

        // Skip over tiles that have exclusions, since they were already processed.

        SYNC_WARPS;
        while (skipTiles[tbx+TILE_SIZE-1] < pos) {
            SYNC_WARPS;
            if (skipBase+tgx < NUM_TILES_WITH_EXCLUSIONS) {
                ushort2 tile = exclusionTiles[skipBase+tgx];
                skipTiles[LOCAL_ID] = tile.x + tile.y*NUM_BLOCKS - tile.y*(tile.y+1)/2;
            }
            else
                skipTiles[LOCAL_ID] = end;
            skipBase += TILE_SIZE;            
            currentSkipIndex = tbx;
            SYNC_WARPS;
        }
        while (skipTiles[currentSkipIndex] < pos)
            currentSkipIndex++;
        includeTile = (skipTiles[currentSkipIndex] != pos);
#endif
        if (includeTile) {
            unsigned int atom1 = x*TILE_SIZE + tgx;

            // Load atom data for this tile.

            real4 posq1 = posq[atom1];
            real charge1 = charge[atom1];
            float2 params1 = global_params[atom1];
#ifdef USE_CUTOFF
            unsigned int j = interactingAtoms[pos*TILE_SIZE+tgx];
#else
            unsigned int j = y*TILE_SIZE + tgx;
#endif
            atomIndices[LOCAL_ID] = j;
            if (j < PADDED_NUM_ATOMS) {
                real4 tempPosq = posq[j];
                localData[LOCAL_ID].x = tempPosq.x;
                localData[LOCAL_ID].y = tempPosq.y;
                localData[LOCAL_ID].z = tempPosq.z;
                localData[LOCAL_ID].q = charge[j];
                float2 tempParams = global_params[j];
                localData[LOCAL_ID].radius = tempParams.x;
                localData[LOCAL_ID].scaledRadius = tempParams.y;
                localData[LOCAL_ID].bornSum = 0.0f;
            }
            SYNC_WARPS;
#ifdef USE_PERIODIC
            if (singlePeriodicCopy) {
                // The box is small enough that we can just translate all the atoms into a single periodic
                // box, then skip having to apply periodic boundary conditions later.

                real4 blockCenterX = blockCenter[x];
                APPLY_PERIODIC_TO_POS_WITH_CENTER(posq1, blockCenterX)
                APPLY_PERIODIC_TO_POS_WITH_CENTER(localData[LOCAL_ID], blockCenterX)
                SYNC_WARPS;
                unsigned int tj = tgx;
                for (j = 0; j < TILE_SIZE; j++) {
                    real3 delta = make_real3(localData[tbx+tj].x-posq1.x, localData[tbx+tj].y-posq1.y, localData[tbx+tj].z-posq1.z);
                    real r2 = delta.x*delta.x + delta.y*delta.y + delta.z*delta.z;
                    int atom2 = atomIndices[tbx+tj];
                    if (atom1 < NUM_ATOMS && atom2 < NUM_ATOMS && r2 < CUTOFF_SQUARED) {
                        real invR = RSQRT(r2);
                        real r = r2*invR;
                        float2 params2 = make_float2(localData[tbx+tj].radius, localData[tbx+tj].scaledRadius);
                        real rScaledRadiusJ = r+params2.y;
                        if (params1.x < rScaledRadiusJ) {
                            real l_ij = RECIP(max((real) params1.x, fabs(r-params2.y)));
                            real u_ij = RECIP(rScaledRadiusJ);
                            real l_ij2 = l_ij*l_ij;
                            real u_ij2 = u_ij*u_ij;
                            real ratio = LOG(u_ij * RECIP(l_ij));
                            bornSum += l_ij - u_ij + (0.50f*invR*ratio) + 0.25f*(r*(u_ij2-l_ij2) +
                                             (params2.y*params2.y*invR)*(l_ij2-u_ij2));
                            bornSum += (params1.x < params2.y-r ? 2.0f*(RECIP(params1.x)-l_ij) : 0);
                        }
                        real rScaledRadiusI = r+params1.y;
                        if (params2.x < rScaledRadiusI) {
                            real l_ij = RECIP(max((real) params2.x, fabs(r-params1.y)));
                            real u_ij = RECIP(rScaledRadiusI);
                            real l_ij2 = l_ij*l_ij;
                            real u_ij2 = u_ij*u_ij;
                            real ratio = LOG(u_ij * RECIP(l_ij));
                            real term = l_ij - u_ij + (0.50f*invR*ratio) + 0.25f*(r*(u_ij2-l_ij2) +
                                             (params1.y*params1.y*invR)*(l_ij2-u_ij2));
                            term += (params2.x < params1.y-r ? 2.0f*(RECIP(params2.x)-l_ij) : 0);
                            localData[tbx+tj].bornSum += term;
                        }
                    }
                    tj = (tj + 1) & (TILE_SIZE - 1);
                    SYNC_WARPS;
                }
            }
            else
#endif
            {
                // We need to apply periodic boundary conditions separately for each interaction.

                unsigned int tj = tgx;
                for (j = 0; j < TILE_SIZE; j++) {
                    real3 delta = make_real3(localData[tbx+tj].x-posq1.x, localData[tbx+tj].y-posq1.y, localData[tbx+tj].z-posq1.z);
#ifdef USE_PERIODIC
                    APPLY_PERIODIC_TO_DELTA(delta)
#endif
                    real r2 = delta.x*delta.x + delta.y*delta.y + delta.z*delta.z;
                    int atom2 = atomIndices[tbx+tj];
#ifdef USE_CUTOFF
                    if (atom1 < NUM_ATOMS && atom2 < NUM_ATOMS && r2 < CUTOFF_SQUARED) {
#else
                    if (atom1 < NUM_ATOMS && atom2 < NUM_ATOMS) {
#endif
                        real invR = RSQRT(r2);
                        real r = r2*invR;
                        float2 params2 = make_float2(localData[tbx+tj].radius, localData[tbx+tj].scaledRadius);
                        real rScaledRadiusJ = r+params2.y;
                        if (params1.x < rScaledRadiusJ) {
                            real l_ij = RECIP(max((real) params1.x, fabs(r-params2.y)));
                            real u_ij = RECIP(rScaledRadiusJ);
                            real l_ij2 = l_ij*l_ij;
                            real u_ij2 = u_ij*u_ij;
                            real ratio = LOG(u_ij * RECIP(l_ij));
                            bornSum += l_ij - u_ij + (0.50f*invR*ratio) + 0.25f*(r*(u_ij2-l_ij2) +
                                             (params2.y*params2.y*invR)*(l_ij2-u_ij2));
                            bornSum += (params1.x < params2.y-r ? 2.0f*(RECIP(params1.x)-l_ij) : 0);
                        }
                        real rScaledRadiusI = r+params1.y;
                        if (params2.x < rScaledRadiusI) {
                            real l_ij = RECIP(max((real) params2.x, fabs(r-params1.y)));
                            real u_ij = RECIP(rScaledRadiusI);
                            real l_ij2 = l_ij*l_ij;
                            real u_ij2 = u_ij*u_ij;
                            real ratio = LOG(u_ij * RECIP(l_ij));
                            real term = l_ij - u_ij + (0.50f*invR*ratio) + 0.25f*(r*(u_ij2-l_ij2) +
                                             (params1.y*params1.y*invR)*(l_ij2-u_ij2));
                            term += (params2.x < params1.y-r ? 2.0f*(RECIP(params2.x)-l_ij) : 0);
                            localData[tbx+tj].bornSum += term;
                        }
                    }
                    tj = (tj + 1) & (TILE_SIZE - 1);
                    SYNC_WARPS;
                }
            }

            // Write results.

#ifdef USE_CUTOFF
            unsigned int atom2 = atomIndices[LOCAL_ID];
#else
            unsigned int atom2 = y*TILE_SIZE + tgx;
#endif
#ifdef SUPPORTS_64_BIT_ATOMICS
            ATOMIC_ADD(&global_bornSum[atom1], (mm_ulong) ((mm_long) (bornSum*0x100000000)));
            if (atom2 < PADDED_NUM_ATOMS)
                ATOMIC_ADD(&global_bornSum[atom2], (mm_ulong) ((mm_long) (localData[LOCAL_ID].bornSum*0x100000000)));
#else
            unsigned int offset1 = atom1 + warp*PADDED_NUM_ATOMS;
            unsigned int offset2 = atom2 + warp*PADDED_NUM_ATOMS;
            global_bornSum[offset1] += bornSum;
            if (atom2 < PADDED_NUM_ATOMS)
                global_bornSum[offset2] += localData[LOCAL_ID].bornSum;
#endif
        }
        pos++;
    }
}

typedef struct {
    real x, y, z;
    real q;
    real fx, fy, fz, fw;
    real bornRadius;
} AtomData2;

/**
 * First part of computing the GBSA interaction.
 */

KERNEL void computeGBSAForce1(
#ifdef SUPPORTS_64_BIT_ATOMICS
        GLOBAL mm_ulong* RESTRICT forceBuffers, GLOBAL mm_ulong* RESTRICT global_bornForce,
#else
        GLOBAL real4* RESTRICT forceBuffers, GLOBAL real* RESTRICT global_bornForce,
#endif
        GLOBAL mixed* RESTRICT energyBuffer, GLOBAL const real4* RESTRICT posq, GLOBAL const real* RESTRICT charge,
        GLOBAL const real* RESTRICT global_bornRadii, int needEnergy,
#ifdef USE_CUTOFF
        GLOBAL const int* RESTRICT tiles, GLOBAL const unsigned int* RESTRICT interactionCount, real4 periodicBoxSize, real4 invPeriodicBoxSize, 
        real4 periodicBoxVecX, real4 periodicBoxVecY, real4 periodicBoxVecZ, unsigned int maxTiles, GLOBAL const real4* RESTRICT blockCenter,
        GLOBAL const real4* RESTRICT blockSize, GLOBAL const int* RESTRICT interactingAtoms,
#else
        unsigned int numTiles,
#endif
        GLOBAL const ushort2* RESTRICT exclusionTiles) {
    const unsigned int totalWarps = GLOBAL_SIZE/TILE_SIZE;
    const unsigned int warp = GLOBAL_ID/TILE_SIZE;
    const unsigned int tgx = LOCAL_ID & (TILE_SIZE-1);
    const unsigned int tbx = LOCAL_ID - tgx;
    mixed energy = 0;
    LOCAL AtomData2 localData[FORCE_WORK_GROUP_SIZE];

    // First loop: process tiles that contain exclusions.
    
    const unsigned int firstExclusionTile = FIRST_EXCLUSION_TILE+warp*(LAST_EXCLUSION_TILE-FIRST_EXCLUSION_TILE)/totalWarps;
    const unsigned int lastExclusionTile = FIRST_EXCLUSION_TILE+(warp+1)*(LAST_EXCLUSION_TILE-FIRST_EXCLUSION_TILE)/totalWarps;
    for (int pos = firstExclusionTile; pos < lastExclusionTile; pos++) {
        const ushort2 tileIndices = exclusionTiles[pos];
        const unsigned int x = tileIndices.x;
        const unsigned int y = tileIndices.y;
        real4 force = make_real4(0);
        unsigned int atom1 = x*TILE_SIZE + tgx;
        real4 posq1 = posq[atom1];
        real charge1 = charge[atom1];
        real bornRadius1 = global_bornRadii[atom1];
        if (x == y) {
            // This tile is on the diagonal.

            localData[LOCAL_ID].x = posq1.x;
            localData[LOCAL_ID].y = posq1.y;
            localData[LOCAL_ID].z = posq1.z;
            localData[LOCAL_ID].q = charge1;
            localData[LOCAL_ID].bornRadius = bornRadius1;
            SYNC_WARPS;
            for (unsigned int j = 0; j < TILE_SIZE; j++) {
                if (atom1 < NUM_ATOMS && y*TILE_SIZE+j < NUM_ATOMS) {
                    real3 pos2 = make_real3(localData[tbx+j].x, localData[tbx+j].y, localData[tbx+j].z);
                    real charge2 = localData[tbx+j].q;
                    real3 delta = make_real3(pos2.x-posq1.x, pos2.y-posq1.y, pos2.z-posq1.z);
#ifdef USE_PERIODIC
                    APPLY_PERIODIC_TO_DELTA(delta)
#endif
                    real r2 = delta.x*delta.x + delta.y*delta.y + delta.z*delta.z;
#ifdef USE_CUTOFF
                    if (r2 < CUTOFF_SQUARED) {
#endif
                        real invR = RSQRT(r2);
                        real r = r2*invR;
                        real bornRadius2 = localData[tbx+j].bornRadius;
                        real alpha2_ij = bornRadius1*bornRadius2;
                        real D_ij = r2*RECIP(4.0f*alpha2_ij);
                        real expTerm = EXP(-D_ij);
                        real denominator2 = r2 + alpha2_ij*expTerm;
                        real denominator = SQRT(denominator2);
                        real scaledChargeProduct = PREFACTOR*charge1*charge2;
                        real tempEnergy = scaledChargeProduct*RECIP(denominator);
                        real Gpol = tempEnergy*RECIP(denominator2);
                        real dGpol_dalpha2_ij = -0.5f*Gpol*expTerm*(1.0f+D_ij);
                        real dEdR = Gpol*(1.0f - 0.25f*expTerm);
                        force.w += dGpol_dalpha2_ij*bornRadius2;
#ifdef USE_CUTOFF
                        if (atom1 != y*TILE_SIZE+j)
                            tempEnergy -= scaledChargeProduct/CUTOFF;
#endif
                        if (needEnergy)
                            energy += 0.5f*tempEnergy;
                        delta *= dEdR;
                        force.x -= delta.x;
                        force.y -= delta.y;
                        force.z -= delta.z;
#ifdef USE_CUTOFF
                    }
#endif
                }
                SYNC_WARPS;
            }
        }
        else {
            // This is an off-diagonal tile.

            unsigned int j = y*TILE_SIZE + tgx;
            real4 tempPosq = posq[j];
            localData[LOCAL_ID].x = tempPosq.x;
            localData[LOCAL_ID].y = tempPosq.y;
            localData[LOCAL_ID].z = tempPosq.z;
            localData[LOCAL_ID].q = charge[j];
            localData[LOCAL_ID].bornRadius = global_bornRadii[j];
            localData[LOCAL_ID].fx = 0.0f;
            localData[LOCAL_ID].fy = 0.0f;
            localData[LOCAL_ID].fz = 0.0f;
            localData[LOCAL_ID].fw = 0.0f;
            SYNC_WARPS;
            unsigned int tj = tgx;
            for (j = 0; j < TILE_SIZE; j++) {
                if (atom1 < NUM_ATOMS && y*TILE_SIZE+tj < NUM_ATOMS) {
                    real3 pos2 = make_real3(localData[tbx+tj].x, localData[tbx+tj].y, localData[tbx+tj].z);
                    real charge2 = localData[tbx+tj].q;
                    real3 delta = make_real3(pos2.x-posq1.x, pos2.y-posq1.y, pos2.z-posq1.z);
#ifdef USE_PERIODIC
                    APPLY_PERIODIC_TO_DELTA(delta)
#endif
                    real r2 = delta.x*delta.x + delta.y*delta.y + delta.z*delta.z;
#ifdef USE_CUTOFF
                    if (r2 < CUTOFF_SQUARED) {
#endif
                        real invR = RSQRT(r2);
                        real r = r2*invR;
                        real bornRadius2 = localData[tbx+tj].bornRadius;
                        real alpha2_ij = bornRadius1*bornRadius2;
                        real D_ij = r2*RECIP(4.0f*alpha2_ij);
                        real expTerm = EXP(-D_ij);
                        real denominator2 = r2 + alpha2_ij*expTerm;
                        real denominator = SQRT(denominator2);
                        real scaledChargeProduct = PREFACTOR*charge1*charge2;
                        real tempEnergy = scaledChargeProduct*RECIP(denominator);
                        real Gpol = tempEnergy*RECIP(denominator2);
                        real dGpol_dalpha2_ij = -0.5f*Gpol*expTerm*(1.0f+D_ij);
                        real dEdR = Gpol*(1.0f - 0.25f*expTerm);
                        force.w += dGpol_dalpha2_ij*bornRadius2;
#ifdef USE_CUTOFF
                        tempEnergy -= scaledChargeProduct/CUTOFF;
#endif
                        if (needEnergy)
                            energy += tempEnergy;
                        delta *= dEdR;
                        force.x -= delta.x;
                        force.y -= delta.y;
                        force.z -= delta.z;
                        localData[tbx+tj].fx += delta.x;
                        localData[tbx+tj].fy += delta.y;
                        localData[tbx+tj].fz += delta.z;
                        localData[tbx+tj].fw += dGpol_dalpha2_ij*bornRadius1;
#ifdef USE_CUTOFF
                    }
#endif
                }
                tj = (tj + 1) & (TILE_SIZE - 1);
                SYNC_WARPS;
            }
        }
        
        // Write results.
        
#ifdef SUPPORTS_64_BIT_ATOMICS
        unsigned int offset = x*TILE_SIZE + tgx;
        ATOMIC_ADD(&forceBuffers[offset], (mm_ulong) ((mm_long) (force.x*0x100000000)));
        ATOMIC_ADD(&forceBuffers[offset+PADDED_NUM_ATOMS], (mm_ulong) ((mm_long) (force.y*0x100000000)));
        ATOMIC_ADD(&forceBuffers[offset+2*PADDED_NUM_ATOMS], (mm_ulong) ((mm_long) (force.z*0x100000000)));
        ATOMIC_ADD(&global_bornForce[offset], (mm_ulong) ((mm_long) (force.w*0x100000000)));
        if (x != y) {
            offset = y*TILE_SIZE + tgx;
            ATOMIC_ADD(&forceBuffers[offset], (mm_ulong) ((mm_long) (localData[LOCAL_ID].fx*0x100000000)));
            ATOMIC_ADD(&forceBuffers[offset+PADDED_NUM_ATOMS], (mm_ulong) ((mm_long) (localData[LOCAL_ID].fy*0x100000000)));
            ATOMIC_ADD(&forceBuffers[offset+2*PADDED_NUM_ATOMS], (mm_ulong) ((mm_long) (localData[LOCAL_ID].fz*0x100000000)));
            ATOMIC_ADD(&global_bornForce[offset], (mm_ulong) ((mm_long) (localData[LOCAL_ID].fw*0x100000000)));
        }
#else
        unsigned int offset1 = x*TILE_SIZE + tgx + warp*PADDED_NUM_ATOMS;
        unsigned int offset2 = y*TILE_SIZE + tgx + warp*PADDED_NUM_ATOMS;
        forceBuffers[offset1] += make_real4(force.x, force.y, force.z, 0);
        global_bornForce[offset1] += force.w;
        if (x != y) {
            forceBuffers[offset2] += (real4) (localData[LOCAL_ID].fx, localData[LOCAL_ID].fy, localData[LOCAL_ID].fz, 0.0f);
            global_bornForce[offset2] += localData[LOCAL_ID].fw;
        }
#endif
    }

    // Second loop: tiles without exclusions, either from the neighbor list (with cutoff) or just enumerating all
    // of them (no cutoff).

#ifdef USE_CUTOFF
    unsigned int numTiles = interactionCount[0];
    if (numTiles > maxTiles)
        return; // There wasn't enough memory for the neighbor list.
    int pos = (int) (warp*(numTiles > maxTiles ? NUM_BLOCKS*((mm_long)NUM_BLOCKS+1)/2 : (mm_long)numTiles)/totalWarps);
    int end = (int) ((warp+1)*(numTiles > maxTiles ? NUM_BLOCKS*((mm_long)NUM_BLOCKS+1)/2 : (mm_long)numTiles)/totalWarps);
#else
    int pos = (int) (warp*(mm_long)numTiles/totalWarps);
    int end = (int) ((warp+1)*(mm_long)numTiles/totalWarps);
#endif
    int skipBase = 0;
    int currentSkipIndex = tbx;
    LOCAL int atomIndices[FORCE_WORK_GROUP_SIZE];
    LOCAL volatile int skipTiles[FORCE_WORK_GROUP_SIZE];
    skipTiles[LOCAL_ID] = -1;

    while (pos < end) {
        real4 force = make_real4(0);
        bool includeTile = true;

        // Extract the coordinates of this tile.
        
        int x, y;
        bool singlePeriodicCopy = false;
#ifdef USE_CUTOFF
        x = tiles[pos];
        real4 blockSizeX = blockSize[x];
        singlePeriodicCopy = (0.5f*periodicBoxSize.x-blockSizeX.x >= CUTOFF &&
                              0.5f*periodicBoxSize.y-blockSizeX.y >= CUTOFF &&
                              0.5f*periodicBoxSize.z-blockSizeX.z >= CUTOFF);
#else
        y = (int) floor(NUM_BLOCKS+0.5f-SQRT((NUM_BLOCKS+0.5f)*(NUM_BLOCKS+0.5f)-2*pos));
        x = (pos-y*NUM_BLOCKS+y*(y+1)/2);
        if (x < y || x >= NUM_BLOCKS) { // Occasionally happens due to roundoff error.
            y += (x < y ? -1 : 1);
            x = (pos-y*NUM_BLOCKS+y*(y+1)/2);
        }

        // Skip over tiles that have exclusions, since they were already processed.

        SYNC_WARPS;
        while (skipTiles[tbx+TILE_SIZE-1] < pos) {
            SYNC_WARPS;
            if (skipBase+tgx < NUM_TILES_WITH_EXCLUSIONS) {
                ushort2 tile = exclusionTiles[skipBase+tgx];
                skipTiles[LOCAL_ID] = tile.x + tile.y*NUM_BLOCKS - tile.y*(tile.y+1)/2;
            }
            else
                skipTiles[LOCAL_ID] = end;
            skipBase += TILE_SIZE;            
            currentSkipIndex = tbx;
            SYNC_WARPS;
        }
        while (skipTiles[currentSkipIndex] < pos)
            currentSkipIndex++;
        includeTile = (skipTiles[currentSkipIndex] != pos);
#endif
        if (includeTile) {
            unsigned int atom1 = x*TILE_SIZE + tgx;

            // Load atom data for this tile.
            
            real4 posq1 = posq[atom1];
            real charge1 = charge[atom1];
            real bornRadius1 = global_bornRadii[atom1];
#ifdef USE_CUTOFF
            unsigned int j = interactingAtoms[pos*TILE_SIZE+tgx];
#else
            unsigned int j = y*TILE_SIZE + tgx;
#endif
            atomIndices[LOCAL_ID] = j;
            if (j < PADDED_NUM_ATOMS) {
                real4 tempPosq = posq[j];
                localData[LOCAL_ID].x = tempPosq.x;
                localData[LOCAL_ID].y = tempPosq.y;
                localData[LOCAL_ID].z = tempPosq.z;
                localData[LOCAL_ID].q = charge[j];
                localData[LOCAL_ID].bornRadius = global_bornRadii[j];
                localData[LOCAL_ID].fx = 0.0f;
                localData[LOCAL_ID].fy = 0.0f;
                localData[LOCAL_ID].fz = 0.0f;
                localData[LOCAL_ID].fw = 0.0f;
            }
            SYNC_WARPS;
#ifdef USE_PERIODIC
            if (singlePeriodicCopy) {
                // The box is small enough that we can just translate all the atoms into a single periodic
                // box, then skip having to apply periodic boundary conditions later.

                real4 blockCenterX = blockCenter[x];
                APPLY_PERIODIC_TO_POS_WITH_CENTER(posq1, blockCenterX)
                APPLY_PERIODIC_TO_POS_WITH_CENTER(localData[LOCAL_ID], blockCenterX)
                SYNC_WARPS;
                unsigned int tj = tgx;
                for (j = 0; j < TILE_SIZE; j++) {
                    int atom2 = atomIndices[tbx+tj];
                    if (atom1 < NUM_ATOMS && atom2 < NUM_ATOMS) {
                        real3 pos2 = make_real3(localData[tbx+tj].x, localData[tbx+tj].y, localData[tbx+tj].z);
                        real charge2 = localData[tbx+tj].q;
                        real3 delta = make_real3(pos2.x-posq1.x, pos2.y-posq1.y, pos2.z-posq1.z);
                        real r2 = delta.x*delta.x + delta.y*delta.y + delta.z*delta.z;
                        if (r2 < CUTOFF_SQUARED) {
                            real invR = RSQRT(r2);
                            real r = r2*invR;
                            real bornRadius2 = localData[tbx+tj].bornRadius;
                            real alpha2_ij = bornRadius1*bornRadius2;
                            real D_ij = r2*RECIP(4.0f*alpha2_ij);
                            real expTerm = EXP(-D_ij);
                            real denominator2 = r2 + alpha2_ij*expTerm;
                            real denominator = SQRT(denominator2);
                            real scaledChargeProduct = PREFACTOR*charge1*charge2;
                            real tempEnergy = scaledChargeProduct*RECIP(denominator);
                            real Gpol = tempEnergy*RECIP(denominator2);
                            real dGpol_dalpha2_ij = -0.5f*Gpol*expTerm*(1.0f+D_ij);
                            real dEdR = Gpol*(1.0f - 0.25f*expTerm);
                            force.w += dGpol_dalpha2_ij*bornRadius2;
#ifdef USE_CUTOFF
                            tempEnergy -= scaledChargeProduct/CUTOFF;
#endif
                            if (needEnergy)
                                energy += tempEnergy;
                            delta *= dEdR;
                            force.x -= delta.x;
                            force.y -= delta.y;
                            force.z -= delta.z;
                            localData[tbx+tj].fx += delta.x;
                            localData[tbx+tj].fy += delta.y;
                            localData[tbx+tj].fz += delta.z;
                            localData[tbx+tj].fw += dGpol_dalpha2_ij*bornRadius1;
                        }
                    }
                    tj = (tj + 1) & (TILE_SIZE - 1);
                    SYNC_WARPS;
                }
            }
            else
#endif
            {
                // We need to apply periodic boundary conditions separately for each interaction.

                unsigned int tj = tgx;
                for (j = 0; j < TILE_SIZE; j++) {
                    int atom2 = atomIndices[tbx+tj];
                    if (atom1 < NUM_ATOMS && atom2 < NUM_ATOMS) {
                        real3 pos2 = make_real3(localData[tbx+tj].x, localData[tbx+tj].y, localData[tbx+tj].z);
                        real charge2 = localData[tbx+tj].q;
                        real3 delta = make_real3(pos2.x-posq1.x, pos2.y-posq1.y, pos2.z-posq1.z);
#ifdef USE_PERIODIC
                        APPLY_PERIODIC_TO_DELTA(delta)
#endif
                        real r2 = delta.x*delta.x + delta.y*delta.y + delta.z*delta.z;
#ifdef USE_CUTOFF
                        if (r2 < CUTOFF_SQUARED) {
#endif
                            real invR = RSQRT(r2);
                            real r = r2*invR;
                            real bornRadius2 = localData[tbx+tj].bornRadius;
                            real alpha2_ij = bornRadius1*bornRadius2;
                            real D_ij = r2*RECIP(4.0f*alpha2_ij);
                            real expTerm = EXP(-D_ij);
                            real denominator2 = r2 + alpha2_ij*expTerm;
                            real denominator = SQRT(denominator2);
                            real scaledChargeProduct = PREFACTOR*charge1*charge2;
                            real tempEnergy = scaledChargeProduct*RECIP(denominator);
                            real Gpol = tempEnergy*RECIP(denominator2);
                            real dGpol_dalpha2_ij = -0.5f*Gpol*expTerm*(1.0f+D_ij);
                            real dEdR = Gpol*(1.0f - 0.25f*expTerm);
                            force.w += dGpol_dalpha2_ij*bornRadius2;
#ifdef USE_CUTOFF
                            tempEnergy -= scaledChargeProduct/CUTOFF;
#endif
                            if (needEnergy)
                                energy += tempEnergy;
                            delta *= dEdR;
                            force.x -= delta.x;
                            force.y -= delta.y;
                            force.z -= delta.z;
                            localData[tbx+tj].fx += delta.x;
                            localData[tbx+tj].fy += delta.y;
                            localData[tbx+tj].fz += delta.z;
                            localData[tbx+tj].fw += dGpol_dalpha2_ij*bornRadius1;
#ifdef USE_CUTOFF
                        }
#endif
                    }
                    tj = (tj + 1) & (TILE_SIZE - 1);
                    SYNC_WARPS;
                }
            }

            // Write results.

#ifdef USE_CUTOFF
            unsigned int atom2 = atomIndices[LOCAL_ID];
#else
            unsigned int atom2 = y*TILE_SIZE + tgx;
#endif
#ifdef SUPPORTS_64_BIT_ATOMICS
            ATOMIC_ADD(&forceBuffers[atom1], (mm_ulong) ((mm_long) (force.x*0x100000000)));
            ATOMIC_ADD(&forceBuffers[atom1+PADDED_NUM_ATOMS], (mm_ulong) ((mm_long) (force.y*0x100000000)));
            ATOMIC_ADD(&forceBuffers[atom1+2*PADDED_NUM_ATOMS], (mm_ulong) ((mm_long) (force.z*0x100000000)));
            ATOMIC_ADD(&global_bornForce[atom1], (mm_ulong) ((mm_long) (force.w*0x100000000)));
            if (atom2 < PADDED_NUM_ATOMS) {
                ATOMIC_ADD(&forceBuffers[atom2], (mm_ulong) ((mm_long) (localData[LOCAL_ID].fx*0x100000000)));
                ATOMIC_ADD(&forceBuffers[atom2+PADDED_NUM_ATOMS], (mm_ulong) ((mm_long) (localData[LOCAL_ID].fy*0x100000000)));
                ATOMIC_ADD(&forceBuffers[atom2+2*PADDED_NUM_ATOMS], (mm_ulong) ((mm_long) (localData[LOCAL_ID].fz*0x100000000)));
                ATOMIC_ADD(&global_bornForce[atom2], (mm_ulong) ((mm_long) (localData[LOCAL_ID].fw*0x100000000)));
            }
#else
            unsigned int offset1 = atom1 + warp*PADDED_NUM_ATOMS;
            unsigned int offset2 = atom2 + warp*PADDED_NUM_ATOMS;
            forceBuffers[offset1] += make_real4(force.x, force.y, force.z, 0);
            global_bornForce[offset1] += force.w;
            if (atom2 < PADDED_NUM_ATOMS) {
                forceBuffers[offset2] += (real4) (localData[LOCAL_ID].fx, localData[LOCAL_ID].fy, localData[LOCAL_ID].fz, 0.0f);
                global_bornForce[offset2] += localData[LOCAL_ID].fw;
            }
#endif
        }
        pos++;
    }
    energyBuffer[GLOBAL_ID] += energy;
}
