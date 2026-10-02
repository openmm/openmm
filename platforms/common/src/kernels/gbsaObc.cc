#define WARPS_PER_GROUP (FORCE_WORK_GROUP_SIZE/TILE_SIZE)

#if defined(USE_HIP)
    #define ALIGN alignas(16)
#else
    #define ALIGN
#endif

typedef struct ALIGN {
    real x, y, z;
    real q;
    float radius, scaledRadius;
    real bornSum;
} AtomData1;

/*
 * The formulas below are shared by all platforms.  These hooks select only
 * lane transport: local arrays by default, or Metal's register shuffle path.
 */
#if USE_GBSA_BORN_SHUFFLE
#if TILE_SIZE != 32
#error GBSA register transport requires TILE_SIZE=32
#endif
#define GBSA_BORN_DATA(index) localData
#define GBSA_BORN_DIAGONAL_DATA(index) broadcastData
#define GBSA_BORN_ATOM_INDEX(index) atomIndices
#define GBSA_BORN_SKIP_TILE(index) simdShuffle(skipTile, (unsigned int) ((index)-tbx))
#else
#define GBSA_BORN_DATA(index) localData[index]
#define GBSA_BORN_DIAGONAL_DATA(index) localData[tbx+index]
#define GBSA_BORN_ATOM_INDEX(index) atomIndices[index]
#define GBSA_BORN_SKIP_TILE(index) skipTiles[index]
#endif


/**
 * Compute the Born sum.
 */
KERNEL void computeBornSum(
        GLOBAL mm_ulong* RESTRICT global_bornSum,
        GLOBAL const real4* RESTRICT posq, GLOBAL const real* RESTRICT charge, GLOBAL const float2* RESTRICT global_params,
#ifdef USE_CUTOFF
        GLOBAL const int* RESTRICT tiles, GLOBAL const unsigned int* RESTRICT interactionCount, real4 periodicBoxSize, real4 invPeriodicBoxSize,
        real4 periodicBoxVecX, real4 periodicBoxVecY, real4 periodicBoxVecZ, unsigned int maxTiles, GLOBAL const real4* RESTRICT blockCenter,
        GLOBAL const real4* RESTRICT blockSize, GLOBAL const int* RESTRICT interactingAtoms,
#else
        unsigned int numTiles,
#endif
        GLOBAL const int2* RESTRICT exclusionTiles) {
    const unsigned int totalWarps = GLOBAL_SIZE/TILE_SIZE;
    const unsigned int warp = GLOBAL_ID/TILE_SIZE;
    const unsigned int tgx = LOCAL_ID & (TILE_SIZE-1);
    const unsigned int tbx = LOCAL_ID - tgx;
#if USE_GBSA_BORN_SHUFFLE
    AtomData1 localData = {};
    int atomIndices = 0;
#else
    LOCAL AtomData1 localData[FORCE_WORK_GROUP_SIZE];
#endif

    // First loop: process tiles that contain exclusions.
    
    const unsigned int firstExclusionTile = FIRST_EXCLUSION_TILE+warp*(LAST_EXCLUSION_TILE-FIRST_EXCLUSION_TILE)/totalWarps;
    const unsigned int lastExclusionTile = FIRST_EXCLUSION_TILE+(warp+1)*(LAST_EXCLUSION_TILE-FIRST_EXCLUSION_TILE)/totalWarps;
    for (int pos = firstExclusionTile; pos < lastExclusionTile; pos++) {
        const int2 tileIndices = exclusionTiles[pos];
        const unsigned int x = tileIndices.x;
        const unsigned int y = tileIndices.y;
        real bornSum = 0;
        unsigned int atom1 = x*TILE_SIZE + tgx;
        real4 posq1 = posq[atom1];
        real charge1 = charge[atom1];
        float2 params1 = global_params[atom1];
        if (x == y) {
            // This tile is on the diagonal.

            GBSA_BORN_DATA(LOCAL_ID).x = posq1.x;
            GBSA_BORN_DATA(LOCAL_ID).y = posq1.y;
            GBSA_BORN_DATA(LOCAL_ID).z = posq1.z;
            GBSA_BORN_DATA(LOCAL_ID).q = charge1;
            GBSA_BORN_DATA(LOCAL_ID).radius = params1.x;
            GBSA_BORN_DATA(LOCAL_ID).scaledRadius = params1.y;
#if !USE_GBSA_BORN_SHUFFLE
            SYNC_WARPS;
#endif
            for (unsigned int j = 0; j < TILE_SIZE; j++) {
#if USE_GBSA_BORN_SHUFFLE
                // Every lane broadcasts before any particle-validity branch.
                AtomData1 broadcastData = metalGbsaBroadcastBorn(localData, j);
#endif
                real3 delta = make_real3(GBSA_BORN_DIAGONAL_DATA(j).x-posq1.x, GBSA_BORN_DIAGONAL_DATA(j).y-posq1.y, GBSA_BORN_DIAGONAL_DATA(j).z-posq1.z);
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
                    float2 params2 = make_float2(GBSA_BORN_DIAGONAL_DATA(j).radius, GBSA_BORN_DIAGONAL_DATA(j).scaledRadius);
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
#if !USE_GBSA_BORN_SHUFFLE
                SYNC_WARPS;
#endif
            }
        }
        else {
            // This is an off-diagonal tile.

            unsigned int j = y*TILE_SIZE + tgx;
            real4 tempPosq = posq[j];
            GBSA_BORN_DATA(LOCAL_ID).x = tempPosq.x;
            GBSA_BORN_DATA(LOCAL_ID).y = tempPosq.y;
            GBSA_BORN_DATA(LOCAL_ID).z = tempPosq.z;
            GBSA_BORN_DATA(LOCAL_ID).q = charge[j];
            float2 tempParams = global_params[j];
            GBSA_BORN_DATA(LOCAL_ID).radius = tempParams.x;
            GBSA_BORN_DATA(LOCAL_ID).scaledRadius = tempParams.y;
            GBSA_BORN_DATA(LOCAL_ID).bornSum = 0.0f;
            SYNC_WARPS;

            // Compute the full set of interactions in this tile.

            unsigned int tj = tgx;
            for (j = 0; j < TILE_SIZE; j++) {
                real3 delta = make_real3(GBSA_BORN_DATA(tbx+tj).x-posq1.x, GBSA_BORN_DATA(tbx+tj).y-posq1.y, GBSA_BORN_DATA(tbx+tj).z-posq1.z);
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
                    float2 params2 = make_float2(GBSA_BORN_DATA(tbx+tj).radius, GBSA_BORN_DATA(tbx+tj).scaledRadius);
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
                        GBSA_BORN_DATA(tbx+tj).bornSum += term;
                    }
                }
                tj = (tj + 1) & (TILE_SIZE - 1);
#if USE_GBSA_BORN_SHUFFLE
                // Accumulators travel with their secondary particle for a full ring.
                localData = metalGbsaRotateBorn(localData, (tgx+1)&31);
                atomIndices = simdShuffle(atomIndices, (tgx+1)&31);
#else
                SYNC_WARPS;
#endif
            }
        }

        // Write results.

        unsigned int offset = x*TILE_SIZE + tgx;
        ATOMIC_ADD(&global_bornSum[offset], (mm_ulong) realToFixedPoint(bornSum));
        if (x != y) {
            offset = y*TILE_SIZE + tgx;
            ATOMIC_ADD(&global_bornSum[offset], (mm_ulong) realToFixedPoint(GBSA_BORN_DATA(LOCAL_ID).bornSum));
        }
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
#if !USE_GBSA_BORN_SHUFFLE
    LOCAL int atomIndices[FORCE_WORK_GROUP_SIZE];
#endif
#if USE_GBSA_BORN_SHUFFLE
    int skipTile = -1;
#else
    LOCAL volatile int skipTiles[FORCE_WORK_GROUP_SIZE];
    skipTiles[LOCAL_ID] = -1;
#endif

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

#if !USE_GBSA_BORN_SHUFFLE
        SYNC_WARPS;
#endif
        while (GBSA_BORN_SKIP_TILE(tbx+TILE_SIZE-1) < pos) {
#if !USE_GBSA_BORN_SHUFFLE
            SYNC_WARPS;
#endif
            if (skipBase+tgx < NUM_TILES_WITH_EXCLUSIONS) {
                int2 tile = exclusionTiles[skipBase+tgx];
#if USE_GBSA_BORN_SHUFFLE
                skipTile = tile.x + tile.y*NUM_BLOCKS - tile.y*(tile.y+1)/2;
#else
                skipTiles[LOCAL_ID] = tile.x + tile.y*NUM_BLOCKS - tile.y*(tile.y+1)/2;
#endif
            }
            else
#if USE_GBSA_BORN_SHUFFLE
                skipTile = end;
#else
                skipTiles[LOCAL_ID] = end;
#endif
            skipBase += TILE_SIZE;            
            currentSkipIndex = tbx;
#if !USE_GBSA_BORN_SHUFFLE
            SYNC_WARPS;
#endif
        }
        while (GBSA_BORN_SKIP_TILE(currentSkipIndex) < pos)
            currentSkipIndex++;
        includeTile = (GBSA_BORN_SKIP_TILE(currentSkipIndex) != pos);
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
            GBSA_BORN_ATOM_INDEX(LOCAL_ID) = j;
            if (j < PADDED_NUM_ATOMS) {
                real4 tempPosq = posq[j];
                GBSA_BORN_DATA(LOCAL_ID).x = tempPosq.x;
                GBSA_BORN_DATA(LOCAL_ID).y = tempPosq.y;
                GBSA_BORN_DATA(LOCAL_ID).z = tempPosq.z;
                GBSA_BORN_DATA(LOCAL_ID).q = charge[j];
                float2 tempParams = global_params[j];
                GBSA_BORN_DATA(LOCAL_ID).radius = tempParams.x;
                GBSA_BORN_DATA(LOCAL_ID).scaledRadius = tempParams.y;
                GBSA_BORN_DATA(LOCAL_ID).bornSum = 0.0f;
            }
            SYNC_WARPS;
#ifdef USE_PERIODIC
            if (singlePeriodicCopy) {
                // The box is small enough that we can just translate all the atoms into a single periodic
                // box, then skip having to apply periodic boundary conditions later.

                real4 blockCenterX = blockCenter[x];
                APPLY_PERIODIC_TO_POS_WITH_CENTER(posq1, blockCenterX)
                APPLY_PERIODIC_TO_POS_WITH_CENTER(GBSA_BORN_DATA(LOCAL_ID), blockCenterX)
                SYNC_WARPS;
                unsigned int tj = tgx;
                for (j = 0; j < TILE_SIZE; j++) {
                    real3 delta = make_real3(GBSA_BORN_DATA(tbx+tj).x-posq1.x, GBSA_BORN_DATA(tbx+tj).y-posq1.y, GBSA_BORN_DATA(tbx+tj).z-posq1.z);
                    real r2 = delta.x*delta.x + delta.y*delta.y + delta.z*delta.z;
                    int atom2 = GBSA_BORN_ATOM_INDEX(tbx+tj);
                    if (atom1 < NUM_ATOMS && atom2 < NUM_ATOMS && r2 < CUTOFF_SQUARED) {
                        real invR = RSQRT(r2);
                        real r = r2*invR;
                        float2 params2 = make_float2(GBSA_BORN_DATA(tbx+tj).radius, GBSA_BORN_DATA(tbx+tj).scaledRadius);
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
                            GBSA_BORN_DATA(tbx+tj).bornSum += term;
                        }
                    }
                    tj = (tj + 1) & (TILE_SIZE - 1);
#if USE_GBSA_BORN_SHUFFLE
                    // Accumulators travel with their secondary particle for a full ring.
                    localData = metalGbsaRotateBorn(localData, (tgx+1)&31);
                    atomIndices = simdShuffle(atomIndices, (tgx+1)&31);
#else
                    SYNC_WARPS;
#endif
                }
            }
            else
#endif
            {
                // We need to apply periodic boundary conditions separately for each interaction.

                unsigned int tj = tgx;
                for (j = 0; j < TILE_SIZE; j++) {
                    real3 delta = make_real3(GBSA_BORN_DATA(tbx+tj).x-posq1.x, GBSA_BORN_DATA(tbx+tj).y-posq1.y, GBSA_BORN_DATA(tbx+tj).z-posq1.z);
#ifdef USE_PERIODIC
                    APPLY_PERIODIC_TO_DELTA(delta)
#endif
                    real r2 = delta.x*delta.x + delta.y*delta.y + delta.z*delta.z;
                    int atom2 = GBSA_BORN_ATOM_INDEX(tbx+tj);
#ifdef USE_CUTOFF
                    if (atom1 < NUM_ATOMS && atom2 < NUM_ATOMS && r2 < CUTOFF_SQUARED) {
#else
                    if (atom1 < NUM_ATOMS && atom2 < NUM_ATOMS) {
#endif
                        real invR = RSQRT(r2);
                        real r = r2*invR;
                        float2 params2 = make_float2(GBSA_BORN_DATA(tbx+tj).radius, GBSA_BORN_DATA(tbx+tj).scaledRadius);
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
                            GBSA_BORN_DATA(tbx+tj).bornSum += term;
                        }
                    }
                    tj = (tj + 1) & (TILE_SIZE - 1);
#if USE_GBSA_BORN_SHUFFLE
                    // Accumulators travel with their secondary particle for a full ring.
                    localData = metalGbsaRotateBorn(localData, (tgx+1)&31);
                    atomIndices = simdShuffle(atomIndices, (tgx+1)&31);
#else
                    SYNC_WARPS;
#endif
                }
            }

            // Write results.

#ifdef USE_CUTOFF
            unsigned int atom2 = GBSA_BORN_ATOM_INDEX(LOCAL_ID);
#else
            unsigned int atom2 = y*TILE_SIZE + tgx;
#endif
            ATOMIC_ADD(&global_bornSum[atom1], (mm_ulong) realToFixedPoint(bornSum));
            if (atom2 < PADDED_NUM_ATOMS)
                ATOMIC_ADD(&global_bornSum[atom2], (mm_ulong) realToFixedPoint(GBSA_BORN_DATA(LOCAL_ID).bornSum));
        }
        pos++;
    }
}


#undef GBSA_BORN_DATA
#undef GBSA_BORN_DIAGONAL_DATA
#undef GBSA_BORN_ATOM_INDEX
#undef GBSA_BORN_SKIP_TILE
typedef struct ALIGN {
    real x, y, z;
    real q;
    real fx, fy, fz, fw;
    real bornRadius;
} AtomData2;

/*
 * The formulas below are shared by all platforms.  These hooks select only
 * lane transport: local arrays by default, or Metal's register shuffle path.
 */
#if USE_GBSA_FORCE_SHUFFLE
#if TILE_SIZE != 32
#error GBSA register transport requires TILE_SIZE=32
#endif
#define GBSA_FORCE_DATA(index) localData
#define GBSA_FORCE_DIAGONAL_DATA(index) broadcastData
#define GBSA_FORCE_ATOM_INDEX(index) atomIndices
#define GBSA_FORCE_SKIP_TILE(index) simdShuffle(skipTile, (unsigned int) ((index)-tbx))
#else
#define GBSA_FORCE_DATA(index) localData[index]
#define GBSA_FORCE_DIAGONAL_DATA(index) localData[tbx+index]
#define GBSA_FORCE_ATOM_INDEX(index) atomIndices[index]
#define GBSA_FORCE_SKIP_TILE(index) skipTiles[index]
#endif


/**
 * First part of computing the GBSA interaction.
 */

KERNEL void computeGBSAForce1(
        GLOBAL mm_ulong* RESTRICT forceBuffers, GLOBAL mm_ulong* RESTRICT global_bornForce,
        GLOBAL mixed* RESTRICT energyBuffer, GLOBAL const real4* RESTRICT posq, GLOBAL const real* RESTRICT charge,
        GLOBAL const real* RESTRICT global_bornRadii, int needEnergy,
#ifdef USE_CUTOFF
        GLOBAL const int* RESTRICT tiles, GLOBAL const unsigned int* RESTRICT interactionCount, real4 periodicBoxSize, real4 invPeriodicBoxSize, 
        real4 periodicBoxVecX, real4 periodicBoxVecY, real4 periodicBoxVecZ, unsigned int maxTiles, GLOBAL const real4* RESTRICT blockCenter,
        GLOBAL const real4* RESTRICT blockSize, GLOBAL const int* RESTRICT interactingAtoms,
#else
        unsigned int numTiles,
#endif
        GLOBAL const int2* RESTRICT exclusionTiles) {
    const unsigned int totalWarps = GLOBAL_SIZE/TILE_SIZE;
    const unsigned int warp = GLOBAL_ID/TILE_SIZE;
    const unsigned int tgx = LOCAL_ID & (TILE_SIZE-1);
    const unsigned int tbx = LOCAL_ID - tgx;
    mixed energy = 0;
#if USE_GBSA_FORCE_SHUFFLE
    AtomData2 localData = {};
    int atomIndices = 0;
#else
    LOCAL AtomData2 localData[FORCE_WORK_GROUP_SIZE];
#endif

    // First loop: process tiles that contain exclusions.
    
    const unsigned int firstExclusionTile = FIRST_EXCLUSION_TILE+warp*(LAST_EXCLUSION_TILE-FIRST_EXCLUSION_TILE)/totalWarps;
    const unsigned int lastExclusionTile = FIRST_EXCLUSION_TILE+(warp+1)*(LAST_EXCLUSION_TILE-FIRST_EXCLUSION_TILE)/totalWarps;
    for (int pos = firstExclusionTile; pos < lastExclusionTile; pos++) {
        const int2 tileIndices = exclusionTiles[pos];
        const unsigned int x = tileIndices.x;
        const unsigned int y = tileIndices.y;
        real4 force = make_real4(0);
        unsigned int atom1 = x*TILE_SIZE + tgx;
        real4 posq1 = posq[atom1];
        real charge1 = charge[atom1];
        real bornRadius1 = global_bornRadii[atom1];
        if (x == y) {
            // This tile is on the diagonal.

            GBSA_FORCE_DATA(LOCAL_ID).x = posq1.x;
            GBSA_FORCE_DATA(LOCAL_ID).y = posq1.y;
            GBSA_FORCE_DATA(LOCAL_ID).z = posq1.z;
            GBSA_FORCE_DATA(LOCAL_ID).q = charge1;
            GBSA_FORCE_DATA(LOCAL_ID).bornRadius = bornRadius1;
#if !USE_GBSA_FORCE_SHUFFLE
            SYNC_WARPS;
#endif
            for (unsigned int j = 0; j < TILE_SIZE; j++) {
#if USE_GBSA_FORCE_SHUFFLE
                // Every lane broadcasts before any particle-validity branch.
                AtomData2 broadcastData = metalGbsaBroadcastForce(localData, j);
#endif
                if (atom1 < NUM_ATOMS && y*TILE_SIZE+j < NUM_ATOMS) {
                    real3 pos2 = make_real3(GBSA_FORCE_DIAGONAL_DATA(j).x, GBSA_FORCE_DIAGONAL_DATA(j).y, GBSA_FORCE_DIAGONAL_DATA(j).z);
                    real charge2 = GBSA_FORCE_DIAGONAL_DATA(j).q;
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
                        real bornRadius2 = GBSA_FORCE_DIAGONAL_DATA(j).bornRadius;
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
#if !USE_GBSA_FORCE_SHUFFLE
                SYNC_WARPS;
#endif
            }
        }
        else {
            // This is an off-diagonal tile.

            unsigned int j = y*TILE_SIZE + tgx;
            real4 tempPosq = posq[j];
            GBSA_FORCE_DATA(LOCAL_ID).x = tempPosq.x;
            GBSA_FORCE_DATA(LOCAL_ID).y = tempPosq.y;
            GBSA_FORCE_DATA(LOCAL_ID).z = tempPosq.z;
            GBSA_FORCE_DATA(LOCAL_ID).q = charge[j];
            GBSA_FORCE_DATA(LOCAL_ID).bornRadius = global_bornRadii[j];
            GBSA_FORCE_DATA(LOCAL_ID).fx = 0.0f;
            GBSA_FORCE_DATA(LOCAL_ID).fy = 0.0f;
            GBSA_FORCE_DATA(LOCAL_ID).fz = 0.0f;
            GBSA_FORCE_DATA(LOCAL_ID).fw = 0.0f;
            SYNC_WARPS;
            unsigned int tj = tgx;
            for (j = 0; j < TILE_SIZE; j++) {
                if (atom1 < NUM_ATOMS && y*TILE_SIZE+tj < NUM_ATOMS) {
                    real3 pos2 = make_real3(GBSA_FORCE_DATA(tbx+tj).x, GBSA_FORCE_DATA(tbx+tj).y, GBSA_FORCE_DATA(tbx+tj).z);
                    real charge2 = GBSA_FORCE_DATA(tbx+tj).q;
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
                        real bornRadius2 = GBSA_FORCE_DATA(tbx+tj).bornRadius;
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
                        GBSA_FORCE_DATA(tbx+tj).fx += delta.x;
                        GBSA_FORCE_DATA(tbx+tj).fy += delta.y;
                        GBSA_FORCE_DATA(tbx+tj).fz += delta.z;
                        GBSA_FORCE_DATA(tbx+tj).fw += dGpol_dalpha2_ij*bornRadius1;
#ifdef USE_CUTOFF
                    }
#endif
                }
                tj = (tj + 1) & (TILE_SIZE - 1);
#if USE_GBSA_FORCE_SHUFFLE
                // Accumulators travel with their secondary particle for a full ring.
                localData = metalGbsaRotateForce(localData, (tgx+1)&31);
                atomIndices = simdShuffle(atomIndices, (tgx+1)&31);
#else
                SYNC_WARPS;
#endif
            }
        }
        
        // Write results.
        
        unsigned int offset = x*TILE_SIZE + tgx;
        ATOMIC_ADD(&forceBuffers[offset], (mm_ulong) realToFixedPoint(force.x));
        ATOMIC_ADD(&forceBuffers[offset+PADDED_NUM_ATOMS], (mm_ulong) realToFixedPoint(force.y));
        ATOMIC_ADD(&forceBuffers[offset+2*PADDED_NUM_ATOMS], (mm_ulong) realToFixedPoint(force.z));
        ATOMIC_ADD(&global_bornForce[offset], (mm_ulong) realToFixedPoint(force.w));
        if (x != y) {
            offset = y*TILE_SIZE + tgx;
            ATOMIC_ADD(&forceBuffers[offset], (mm_ulong) realToFixedPoint(GBSA_FORCE_DATA(LOCAL_ID).fx));
            ATOMIC_ADD(&forceBuffers[offset+PADDED_NUM_ATOMS], (mm_ulong) realToFixedPoint(GBSA_FORCE_DATA(LOCAL_ID).fy));
            ATOMIC_ADD(&forceBuffers[offset+2*PADDED_NUM_ATOMS], (mm_ulong) realToFixedPoint(GBSA_FORCE_DATA(LOCAL_ID).fz));
            ATOMIC_ADD(&global_bornForce[offset], (mm_ulong) realToFixedPoint(GBSA_FORCE_DATA(LOCAL_ID).fw));
        }
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
#if !USE_GBSA_FORCE_SHUFFLE
    LOCAL int atomIndices[FORCE_WORK_GROUP_SIZE];
#endif
#if USE_GBSA_FORCE_SHUFFLE
    int skipTile = -1;
#else
    LOCAL volatile int skipTiles[FORCE_WORK_GROUP_SIZE];
    skipTiles[LOCAL_ID] = -1;
#endif

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

#if !USE_GBSA_FORCE_SHUFFLE
        SYNC_WARPS;
#endif
        while (GBSA_FORCE_SKIP_TILE(tbx+TILE_SIZE-1) < pos) {
#if !USE_GBSA_FORCE_SHUFFLE
            SYNC_WARPS;
#endif
            if (skipBase+tgx < NUM_TILES_WITH_EXCLUSIONS) {
                int2 tile = exclusionTiles[skipBase+tgx];
#if USE_GBSA_FORCE_SHUFFLE
                skipTile = tile.x + tile.y*NUM_BLOCKS - tile.y*(tile.y+1)/2;
#else
                skipTiles[LOCAL_ID] = tile.x + tile.y*NUM_BLOCKS - tile.y*(tile.y+1)/2;
#endif
            }
            else
#if USE_GBSA_FORCE_SHUFFLE
                skipTile = end;
#else
                skipTiles[LOCAL_ID] = end;
#endif
            skipBase += TILE_SIZE;            
            currentSkipIndex = tbx;
#if !USE_GBSA_FORCE_SHUFFLE
            SYNC_WARPS;
#endif
        }
        while (GBSA_FORCE_SKIP_TILE(currentSkipIndex) < pos)
            currentSkipIndex++;
        includeTile = (GBSA_FORCE_SKIP_TILE(currentSkipIndex) != pos);
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
            GBSA_FORCE_ATOM_INDEX(LOCAL_ID) = j;
            if (j < PADDED_NUM_ATOMS) {
                real4 tempPosq = posq[j];
                GBSA_FORCE_DATA(LOCAL_ID).x = tempPosq.x;
                GBSA_FORCE_DATA(LOCAL_ID).y = tempPosq.y;
                GBSA_FORCE_DATA(LOCAL_ID).z = tempPosq.z;
                GBSA_FORCE_DATA(LOCAL_ID).q = charge[j];
                GBSA_FORCE_DATA(LOCAL_ID).bornRadius = global_bornRadii[j];
                GBSA_FORCE_DATA(LOCAL_ID).fx = 0.0f;
                GBSA_FORCE_DATA(LOCAL_ID).fy = 0.0f;
                GBSA_FORCE_DATA(LOCAL_ID).fz = 0.0f;
                GBSA_FORCE_DATA(LOCAL_ID).fw = 0.0f;
            }
            SYNC_WARPS;
#ifdef USE_PERIODIC
            if (singlePeriodicCopy) {
                // The box is small enough that we can just translate all the atoms into a single periodic
                // box, then skip having to apply periodic boundary conditions later.

                real4 blockCenterX = blockCenter[x];
                APPLY_PERIODIC_TO_POS_WITH_CENTER(posq1, blockCenterX)
                APPLY_PERIODIC_TO_POS_WITH_CENTER(GBSA_FORCE_DATA(LOCAL_ID), blockCenterX)
                SYNC_WARPS;
                unsigned int tj = tgx;
                for (j = 0; j < TILE_SIZE; j++) {
                    int atom2 = GBSA_FORCE_ATOM_INDEX(tbx+tj);
                    if (atom1 < NUM_ATOMS && atom2 < NUM_ATOMS) {
                        real3 pos2 = make_real3(GBSA_FORCE_DATA(tbx+tj).x, GBSA_FORCE_DATA(tbx+tj).y, GBSA_FORCE_DATA(tbx+tj).z);
                        real charge2 = GBSA_FORCE_DATA(tbx+tj).q;
                        real3 delta = make_real3(pos2.x-posq1.x, pos2.y-posq1.y, pos2.z-posq1.z);
                        real r2 = delta.x*delta.x + delta.y*delta.y + delta.z*delta.z;
                        if (r2 < CUTOFF_SQUARED) {
                            real invR = RSQRT(r2);
                            real r = r2*invR;
                            real bornRadius2 = GBSA_FORCE_DATA(tbx+tj).bornRadius;
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
                            GBSA_FORCE_DATA(tbx+tj).fx += delta.x;
                            GBSA_FORCE_DATA(tbx+tj).fy += delta.y;
                            GBSA_FORCE_DATA(tbx+tj).fz += delta.z;
                            GBSA_FORCE_DATA(tbx+tj).fw += dGpol_dalpha2_ij*bornRadius1;
                        }
                    }
                    tj = (tj + 1) & (TILE_SIZE - 1);
#if USE_GBSA_FORCE_SHUFFLE
                    // Accumulators travel with their secondary particle for a full ring.
                    localData = metalGbsaRotateForce(localData, (tgx+1)&31);
                    atomIndices = simdShuffle(atomIndices, (tgx+1)&31);
#else
                    SYNC_WARPS;
#endif
                }
            }
            else
#endif
            {
                // We need to apply periodic boundary conditions separately for each interaction.

                unsigned int tj = tgx;
                for (j = 0; j < TILE_SIZE; j++) {
                    int atom2 = GBSA_FORCE_ATOM_INDEX(tbx+tj);
                    if (atom1 < NUM_ATOMS && atom2 < NUM_ATOMS) {
                        real3 pos2 = make_real3(GBSA_FORCE_DATA(tbx+tj).x, GBSA_FORCE_DATA(tbx+tj).y, GBSA_FORCE_DATA(tbx+tj).z);
                        real charge2 = GBSA_FORCE_DATA(tbx+tj).q;
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
                            real bornRadius2 = GBSA_FORCE_DATA(tbx+tj).bornRadius;
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
                            GBSA_FORCE_DATA(tbx+tj).fx += delta.x;
                            GBSA_FORCE_DATA(tbx+tj).fy += delta.y;
                            GBSA_FORCE_DATA(tbx+tj).fz += delta.z;
                            GBSA_FORCE_DATA(tbx+tj).fw += dGpol_dalpha2_ij*bornRadius1;
#ifdef USE_CUTOFF
                        }
#endif
                    }
                    tj = (tj + 1) & (TILE_SIZE - 1);
#if USE_GBSA_FORCE_SHUFFLE
                    // Accumulators travel with their secondary particle for a full ring.
                    localData = metalGbsaRotateForce(localData, (tgx+1)&31);
                    atomIndices = simdShuffle(atomIndices, (tgx+1)&31);
#else
                    SYNC_WARPS;
#endif
                }
            }

            // Write results.

#ifdef USE_CUTOFF
            unsigned int atom2 = GBSA_FORCE_ATOM_INDEX(LOCAL_ID);
#else
            unsigned int atom2 = y*TILE_SIZE + tgx;
#endif
            ATOMIC_ADD(&forceBuffers[atom1], (mm_ulong) realToFixedPoint(force.x));
            ATOMIC_ADD(&forceBuffers[atom1+PADDED_NUM_ATOMS], (mm_ulong) realToFixedPoint(force.y));
            ATOMIC_ADD(&forceBuffers[atom1+2*PADDED_NUM_ATOMS], (mm_ulong) realToFixedPoint(force.z));
            ATOMIC_ADD(&global_bornForce[atom1], (mm_ulong) realToFixedPoint(force.w));
            if (atom2 < PADDED_NUM_ATOMS) {
                ATOMIC_ADD(&forceBuffers[atom2], (mm_ulong) realToFixedPoint(GBSA_FORCE_DATA(LOCAL_ID).fx));
                ATOMIC_ADD(&forceBuffers[atom2+PADDED_NUM_ATOMS], (mm_ulong) realToFixedPoint(GBSA_FORCE_DATA(LOCAL_ID).fy));
                ATOMIC_ADD(&forceBuffers[atom2+2*PADDED_NUM_ATOMS], (mm_ulong) realToFixedPoint(GBSA_FORCE_DATA(LOCAL_ID).fz));
                ATOMIC_ADD(&global_bornForce[atom2], (mm_ulong) realToFixedPoint(GBSA_FORCE_DATA(LOCAL_ID).fw));
            }
        }
        pos++;
    }
    energyBuffer[GLOBAL_ID] += energy;
}

#undef GBSA_FORCE_DATA
#undef GBSA_FORCE_DIAGONAL_DATA
#undef GBSA_FORCE_ATOM_INDEX
#undef GBSA_FORCE_SKIP_TILE
