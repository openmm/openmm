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
        GLOBAL mm_long* RESTRICT global_bornSum,
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
        GLOBAL const ushort2* exclusionTiles) {
    LOCAL AtomData1 localData[TILE_SIZE];

    // First loop: process tiles that contain exclusions.
    
    const unsigned int firstExclusionTile = FIRST_EXCLUSION_TILE+GROUP_ID*(LAST_EXCLUSION_TILE-FIRST_EXCLUSION_TILE)/NUM_GROUPS;
    const unsigned int lastExclusionTile = FIRST_EXCLUSION_TILE+(GROUP_ID+1)*(LAST_EXCLUSION_TILE-FIRST_EXCLUSION_TILE)/NUM_GROUPS;
    for (int pos = firstExclusionTile; pos < lastExclusionTile; pos++) {
        const ushort2 tileIndices = exclusionTiles[pos];
        const unsigned int x = tileIndices.x;
        const unsigned int y = tileIndices.y;

        // Load the data for this tile.

        for (int localAtomIndex = 0; localAtomIndex < TILE_SIZE; localAtomIndex++) {
            unsigned int j = y*TILE_SIZE + localAtomIndex;
            real4 tempPosq = posq[j];
            localData[localAtomIndex].x = tempPosq.x;
            localData[localAtomIndex].y = tempPosq.y;
            localData[localAtomIndex].z = tempPosq.z;
            localData[localAtomIndex].q = charge[j];
            float2 tempParams = global_params[j];
            localData[localAtomIndex].radius = tempParams.x;
            localData[localAtomIndex].scaledRadius = tempParams.y;
        }
        if (x == y) {
            // This tile is on the diagonal.

            for (unsigned int tgx = 0; tgx < TILE_SIZE; tgx++) {
                unsigned int atom1 = x*TILE_SIZE+tgx;
                real bornSum = 0;
                real4 posq1 = posq[atom1];
                float2 params1 = global_params[atom1];
                for (unsigned int j = 0; j < TILE_SIZE; j++) {
                    real3 pos2 = make_real3(localData[j].x, localData[j].y, localData[j].z);
                    real charge2 = localData[j].q;
                    real3 delta = make_real3(pos2.x-posq1.x, pos2.y-posq1.y, pos2.z-posq1.z);
#ifdef USE_PERIODIC
                    APPLY_PERIODIC_TO_DELTA(delta)
#endif
                    real r2 = dot(trimTo3(delta), trimTo3(delta));
#ifdef USE_CUTOFF
                    if (atom1 < NUM_ATOMS && y*TILE_SIZE+j < NUM_ATOMS && r2 < CUTOFF_SQUARED) {
#else
                    if (atom1 < NUM_ATOMS && y*TILE_SIZE+j < NUM_ATOMS) {
#endif
                        real invR = RSQRT(r2);
                        real r = r2*invR;
                        float2 params2 = make_float2(localData[j].radius, localData[j].scaledRadius);
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
                }

                // Write results.

#ifdef SUPPORTS_64_BIT_ATOMICS
                ATOMIC_ADD(&global_bornSum[atom1], (mm_long) (bornSum*0x100000000));
#else
                unsigned int offset = atom1 + GROUP_ID*PADDED_NUM_ATOMS;
                global_bornSum[offset] += bornSum;
#endif
            }
        }
        else {
            // This is an off-diagonal tile.

            for (int tgx = 0; tgx < TILE_SIZE; tgx++)
                localData[tgx].bornSum = 0;
            for (unsigned int tgx = 0; tgx < TILE_SIZE; tgx++) {
                unsigned int atom1 = x*TILE_SIZE+tgx;
                real bornSum = 0;
                real4 posq1 = posq[atom1];
                float2 params1 = global_params[atom1];
                for (unsigned int j = 0; j < TILE_SIZE; j++) {
                    real3 pos2 = make_real3(localData[j].x, localData[j].y, localData[j].z);
                    real charge2 = localData[j].q;
                    real3 delta = make_real3(pos2.x-posq1.x, pos2.y-posq1.y, pos2.z-posq1.z);
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
                        float2 params2 = make_float2(localData[j].radius, localData[j].scaledRadius);
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
                            localData[j].bornSum += term;
                        }
                    }
                }

               // Write results for atom1.

#ifdef SUPPORTS_64_BIT_ATOMICS
                ATOMIC_ADD(&global_bornSum[atom1], (mm_long) (bornSum*0x100000000));
#else
                unsigned int offset = atom1 + GROUP_ID*PADDED_NUM_ATOMS;
                global_bornSum[offset] += bornSum;
#endif
            }

            // Write results.

            for (int tgx = 0; tgx < TILE_SIZE; tgx++) {
#ifdef SUPPORTS_64_BIT_ATOMICS
                unsigned int offset = y*TILE_SIZE + tgx;
                ATOMIC_ADD(&global_bornSum[offset], (mm_long) (localData[tgx].bornSum*0x100000000));
#else
                unsigned int offset = y*TILE_SIZE+tgx + GROUP_ID*PADDED_NUM_ATOMS;
                global_bornSum[offset] += localData[tgx].bornSum;
#endif
            }
        }
    }

    // Second loop: tiles without exclusions, either from the neighbor list (with cutoff) or just enumerating all
    // of them (no cutoff).

#ifdef USE_CUTOFF
    unsigned int numTiles = interactionCount[0];
    if (numTiles > maxTiles)
        return; // There wasn't enough memory for the neighbor list.
    int pos = (int) (GROUP_ID*(numTiles > maxTiles ? NUM_BLOCKS*((mm_long)NUM_BLOCKS+1)/2 : numTiles)/NUM_GROUPS);
    int end = (int) ((GROUP_ID+1)*(numTiles > maxTiles ? NUM_BLOCKS*((mm_long)NUM_BLOCKS+1)/2 : numTiles)/NUM_GROUPS);
#else
    int pos = (int) (GROUP_ID*(mm_long)numTiles/NUM_GROUPS);
    int end = (int) ((GROUP_ID+1)*(mm_long)numTiles/NUM_GROUPS);
#endif
    int nextToSkip = -1;
    int currentSkipIndex = 0;
    LOCAL int atomIndices[TILE_SIZE];

    while (pos < end) {
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

        while (nextToSkip < pos) {
            if (currentSkipIndex < NUM_TILES_WITH_EXCLUSIONS) {
                ushort2 tile = exclusionTiles[currentSkipIndex++];
                nextToSkip = tile.x + tile.y*NUM_BLOCKS - tile.y*(tile.y+1)/2;
            }
            else
                nextToSkip = end;
        }
        includeTile = (nextToSkip != pos);
#endif
        if (includeTile) {
            // Load the data for this tile.

            for (int localAtomIndex = 0; localAtomIndex < TILE_SIZE; localAtomIndex++) {
#ifdef USE_CUTOFF
                unsigned int j = interactingAtoms[pos*TILE_SIZE+localAtomIndex];
#else
                unsigned int j = y*TILE_SIZE+localAtomIndex;
#endif
                atomIndices[localAtomIndex] = j;
                if (j < PADDED_NUM_ATOMS) {
                    real4 tempPosq = posq[j];
                    localData[localAtomIndex].x = tempPosq.x;
                    localData[localAtomIndex].y = tempPosq.y;
                    localData[localAtomIndex].z = tempPosq.z;
                    localData[localAtomIndex].q = charge[j];
                    float2 tempParams = global_params[j];
                    localData[localAtomIndex].radius = tempParams.x;
                    localData[localAtomIndex].scaledRadius = tempParams.y;
                    localData[localAtomIndex].bornSum = 0.0f;
                }
            }
#ifdef USE_PERIODIC
            if (singlePeriodicCopy) {
                // The box is small enough that we can just translate all the atoms into a single periodic
                // box, then skip having to apply periodic boundary conditions later.

                real4 blockCenterX = blockCenter[x];
                for (unsigned int tgx = 0; tgx < TILE_SIZE; tgx++) {
                    APPLY_PERIODIC_TO_POS_WITH_CENTER(localData[tgx], blockCenterX)
                }
                for (unsigned int tgx = 0; tgx < TILE_SIZE; tgx++) {
                    unsigned int atom1 = x*TILE_SIZE+tgx;
                    real bornSum = 0;
                    real4 posq1 = posq[atom1];
                    APPLY_PERIODIC_TO_POS_WITH_CENTER(posq1, blockCenterX)
                    float2 params1 = global_params[atom1];
                    for (unsigned int j = 0; j < TILE_SIZE; j++) {
                        real3 pos2 = make_real3(localData[j].x, localData[j].y, localData[j].z);
                        real charge2 = localData[j].q;
                        real3 delta = make_real3(pos2.x-posq1.x, pos2.y-posq1.y, pos2.z-posq1.z);
                        real r2 = delta.x*delta.x + delta.y*delta.y + delta.z*delta.z;
                        int atom2 = atomIndices[j];
                        if (atom1 < NUM_ATOMS && atom2 < NUM_ATOMS && r2 < CUTOFF_SQUARED) {
                            real invR = RSQRT(r2);
                            real r = r2*invR;
                            float2 params2 = make_float2(localData[j].radius, localData[j].scaledRadius);
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
                                localData[j].bornSum += term;
                            }
                        }
                    }

                    // Write results for atom1.

#ifdef SUPPORTS_64_BIT_ATOMICS
                    ATOMIC_ADD(&global_bornSum[atom1], (mm_long) (bornSum*0x100000000));
#else
                    unsigned int offset = atom1 + GROUP_ID*PADDED_NUM_ATOMS;
                    global_bornSum[offset] += bornSum;
#endif
                }
            }
            else
#endif
            {
                // We need to apply periodic boundary conditions separately for each interaction.

                for (unsigned int tgx = 0; tgx < TILE_SIZE; tgx++) {
                    unsigned int atom1 = x*TILE_SIZE+tgx;
                    real bornSum = 0;
                    real4 posq1 = posq[atom1];
                    float2 params1 = global_params[atom1];
                    for (unsigned int j = 0; j < TILE_SIZE; j++) {
                        real3 pos2 = make_real3(localData[j].x, localData[j].y, localData[j].z);
                        real charge2 = localData[j].q;
                        real3 delta = make_real3(pos2.x-posq1.x, pos2.y-posq1.y, pos2.z-posq1.z);
#ifdef USE_PERIODIC
                        APPLY_PERIODIC_TO_DELTA(delta)
#endif
                        real r2 = delta.x*delta.x + delta.y*delta.y + delta.z*delta.z;
                        int atom2 = atomIndices[j];
#ifdef USE_CUTOFF
                        if (atom1 < NUM_ATOMS && atom2 < NUM_ATOMS && r2 < CUTOFF_SQUARED) {
#else
                        if (atom1 < NUM_ATOMS && atom2 < NUM_ATOMS) {
#endif
                            real invR = RSQRT(r2);
                            real r = r2*invR;
                            float2 params2 = make_float2(localData[j].radius, localData[j].scaledRadius);
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
                                localData[j].bornSum += term;
                            }
                        }
                    }

                    // Write results for atom1.

#ifdef SUPPORTS_64_BIT_ATOMICS
                    ATOMIC_ADD(&global_bornSum[atom1], (mm_long) (bornSum*0x100000000));
#else
                    unsigned int offset = atom1 + GROUP_ID*PADDED_NUM_ATOMS;
                    global_bornSum[offset] += bornSum;
#endif
                }
            }

            // Write results.

            for (int tgx = 0; tgx < TILE_SIZE; tgx++) {
#ifdef USE_CUTOFF
                unsigned int atom2 = atomIndices[tgx];
#else
                unsigned int atom2 = y*TILE_SIZE + tgx;
#endif
                if (atom2 < PADDED_NUM_ATOMS) {
#ifdef SUPPORTS_64_BIT_ATOMICS
                    ATOMIC_ADD(&global_bornSum[atom2], (mm_long) (localData[tgx].bornSum*0x100000000));
#else
                    unsigned int offset = atom2 + GROUP_ID*PADDED_NUM_ATOMS;
                    global_bornSum[offset] += localData[tgx].bornSum;
#endif
                }
            }
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
        GLOBAL mm_long* RESTRICT forceBuffers, GLOBAL mm_long* RESTRICT global_bornForce,
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
        GLOBAL const ushort2* exclusionTiles) {
    mixed energy = 0;
    LOCAL AtomData2 localData[TILE_SIZE];

    // First loop: process tiles that contain exclusions.
    
    const unsigned int firstExclusionTile = FIRST_EXCLUSION_TILE+GROUP_ID*(LAST_EXCLUSION_TILE-FIRST_EXCLUSION_TILE)/NUM_GROUPS;
    const unsigned int lastExclusionTile = FIRST_EXCLUSION_TILE+(GROUP_ID+1)*(LAST_EXCLUSION_TILE-FIRST_EXCLUSION_TILE)/NUM_GROUPS;
    for (int pos = firstExclusionTile; pos < lastExclusionTile; pos++) {
        const ushort2 tileIndices = exclusionTiles[pos];
        const unsigned int x = tileIndices.x;
        const unsigned int y = tileIndices.y;

        // Load the data for this tile.

        for (int localAtomIndex = 0; localAtomIndex < TILE_SIZE; localAtomIndex++) {
            unsigned int j = y*TILE_SIZE + localAtomIndex;
            real4 tempPosq = posq[j];
            localData[localAtomIndex].x = tempPosq.x;
            localData[localAtomIndex].y = tempPosq.y;
            localData[localAtomIndex].z = tempPosq.z;
            localData[localAtomIndex].q = charge[j];
            localData[localAtomIndex].bornRadius = global_bornRadii[j];
        }
        if (x == y) {
            // This tile is on the diagonal.

            for (unsigned int tgx = 0; tgx < TILE_SIZE; tgx++) {
                unsigned int atom1 = x*TILE_SIZE+tgx;
                real4 force = make_real4(0);
                real4 posq1 = posq[atom1];
                real charge1 = charge[atom1];
                real bornRadius1 = global_bornRadii[atom1];
                for (unsigned int j = 0; j < TILE_SIZE; j++) {
                    real3 pos2 = make_real3(localData[j].x, localData[j].y, localData[j].z);
                    real charge2 = localData[j].q;
                    real3 delta = make_real3(pos2.x-posq1.x, pos2.y-posq1.y, pos2.z-posq1.z);
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
                        real bornRadius2 = localData[j].bornRadius;
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
                        energy += 0.5f*tempEnergy;
                        delta *= dEdR;
                        force.x -= delta.x;
                        force.y -= delta.y;
                        force.z -= delta.z;
                    }
                }

                // Write results.

#ifdef SUPPORTS_64_BIT_ATOMICS
                ATOMIC_ADD(&forceBuffers[atom1], (mm_long) (force.x*0x100000000));
                ATOMIC_ADD(&forceBuffers[atom1+PADDED_NUM_ATOMS], (mm_long) (force.y*0x100000000));
                ATOMIC_ADD(&forceBuffers[atom1+2*PADDED_NUM_ATOMS], (mm_long) (force.z*0x100000000));
                ATOMIC_ADD(&global_bornForce[atom1], (mm_long) (force.w*0x100000000));
#else
                unsigned int offset = atom1 + GROUP_ID*PADDED_NUM_ATOMS;
                forceBuffers[offset] += make_real4(force.x, force.y, force.z, 0);
                global_bornForce[offset] += force.w;
#endif
            }
        }
        else {
            // This is an off-diagonal tile.

            for (int tgx = 0; tgx < TILE_SIZE; tgx++) {
                localData[tgx].fx = 0;
                localData[tgx].fy = 0;
                localData[tgx].fz = 0;
                localData[tgx].fw = 0;
            }
            for (unsigned int tgx = 0; tgx < TILE_SIZE; tgx++) {
                unsigned int atom1 = x*TILE_SIZE+tgx;
                real4 force = make_real4(0);
                real4 posq1 = posq[atom1];
                real charge1 = charge[atom1];
                real bornRadius1 = global_bornRadii[atom1];
                for (unsigned int j = 0; j < TILE_SIZE; j++) {
                    real3 pos2 = make_real3(localData[j].x, localData[j].y, localData[j].z);
                    real charge2 = localData[j].q;
                    real3 delta = make_real3(pos2.x-posq1.x, pos2.y-posq1.y, pos2.z-posq1.z);
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
                        real bornRadius2 = localData[j].bornRadius;
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
                        energy += tempEnergy;
                        delta *= dEdR;
                        force.x -= delta.x;
                        force.y -= delta.y;
                        force.z -= delta.z;
                        localData[j].fx += delta.x;
                        localData[j].fy += delta.y;
                        localData[j].fz += delta.z;
                        localData[j].fw += dGpol_dalpha2_ij*bornRadius1;
                    }
                }

               // Write results for atom1.

#ifdef SUPPORTS_64_BIT_ATOMICS
                ATOMIC_ADD(&forceBuffers[atom1], (mm_long) (force.x*0x100000000));
                ATOMIC_ADD(&forceBuffers[atom1+PADDED_NUM_ATOMS], (mm_long) (force.y*0x100000000));
                ATOMIC_ADD(&forceBuffers[atom1+2*PADDED_NUM_ATOMS], (mm_long) (force.z*0x100000000));
                ATOMIC_ADD(&global_bornForce[atom1], (mm_long) (force.w*0x100000000));
#else
                unsigned int offset = atom1 + GROUP_ID*PADDED_NUM_ATOMS;
                forceBuffers[offset] += make_real4(force.x, force.y, force.z, 0);
                global_bornForce[offset] += force.w;
#endif
            }

            // Write results.

            for (int tgx = 0; tgx < TILE_SIZE; tgx++) {
#ifdef SUPPORTS_64_BIT_ATOMICS
                unsigned int offset = y*TILE_SIZE + tgx;
                ATOMIC_ADD(&forceBuffers[offset], (mm_long) (localData[tgx].fx*0x100000000));
                ATOMIC_ADD(&forceBuffers[offset+PADDED_NUM_ATOMS], (mm_long) (localData[tgx].fy*0x100000000));
                ATOMIC_ADD(&forceBuffers[offset+2*PADDED_NUM_ATOMS], (mm_long) (localData[tgx].fz*0x100000000));
                ATOMIC_ADD(&global_bornForce[offset], (mm_long) (localData[tgx].fw*0x100000000));
#else
                unsigned int offset = y*TILE_SIZE+tgx + GROUP_ID*PADDED_NUM_ATOMS;
                real4 f = forceBuffers[offset];
                f.x += localData[tgx].fx;
                f.y += localData[tgx].fy;
                f.z += localData[tgx].fz;
                forceBuffers[offset] = f;
                global_bornForce[offset] += localData[tgx].fw;
#endif
            }
        }
    }

    // Second loop: tiles without exclusions, either from the neighbor list (with cutoff) or just enumerating all
    // of them (no cutoff).

#ifdef USE_CUTOFF
    unsigned int numTiles = interactionCount[0];
    if (numTiles > maxTiles)
        return; // There wasn't enough memory for the neighbor list.
    int pos = (int) (GROUP_ID*(numTiles > maxTiles ? NUM_BLOCKS*((mm_long)NUM_BLOCKS+1)/2 : numTiles)/NUM_GROUPS);
    int end = (int) ((GROUP_ID+1)*(numTiles > maxTiles ? NUM_BLOCKS*((mm_long)NUM_BLOCKS+1)/2 : numTiles)/NUM_GROUPS);
#else
    int pos = (int) (GROUP_ID*(mm_long)numTiles/NUM_GROUPS);
    int end = (int) ((GROUP_ID+1)*(mm_long)numTiles/NUM_GROUPS);
#endif
    int nextToSkip = -1;
    int currentSkipIndex = 0;
    LOCAL int atomIndices[TILE_SIZE];

    while (pos < end) {
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

        while (nextToSkip < pos) {
            if (currentSkipIndex < NUM_TILES_WITH_EXCLUSIONS) {
                ushort2 tile = exclusionTiles[currentSkipIndex++];
                nextToSkip = tile.x + tile.y*NUM_BLOCKS - tile.y*(tile.y+1)/2;
            }
            else
                nextToSkip = end;
        }
        includeTile = (nextToSkip != pos);
#endif
        if (includeTile) {
            // Load the data for this tile.

            for (int localAtomIndex = 0; localAtomIndex < TILE_SIZE; localAtomIndex++) {
#ifdef USE_CUTOFF
                unsigned int j = interactingAtoms[pos*TILE_SIZE+localAtomIndex];
#else
                unsigned int j = y*TILE_SIZE+localAtomIndex;
#endif
                atomIndices[localAtomIndex] = j;
                if (j < PADDED_NUM_ATOMS) {
                    real4 tempPosq = posq[j];
                    localData[localAtomIndex].x = tempPosq.x;
                    localData[localAtomIndex].y = tempPosq.y;
                    localData[localAtomIndex].z = tempPosq.z;
                    localData[localAtomIndex].q = charge[j];
                    localData[localAtomIndex].bornRadius = global_bornRadii[j];
                    localData[localAtomIndex].fx = 0.0f;
                    localData[localAtomIndex].fy = 0.0f;
                    localData[localAtomIndex].fz = 0.0f;
                    localData[localAtomIndex].fw = 0.0f;
                }
            }
#ifdef USE_PERIODIC
            if (singlePeriodicCopy) {
                // The box is small enough that we can just translate all the atoms into a single periodic
                // box, then skip having to apply periodic boundary conditions later.

                real4 blockCenterX = blockCenter[x];
                for (unsigned int tgx = 0; tgx < TILE_SIZE; tgx++) {
                    APPLY_PERIODIC_TO_POS_WITH_CENTER(localData[tgx], blockCenterX)
                }
                for (unsigned int tgx = 0; tgx < TILE_SIZE; tgx++) {
                    unsigned int atom1 = x*TILE_SIZE+tgx;
                    real4 force = make_real4(0);
                    real4 posq1 = posq[atom1];
                    real charge1 = charge[atom1];
                    APPLY_PERIODIC_TO_POS_WITH_CENTER(posq1, blockCenterX)
                    float bornRadius1 = global_bornRadii[atom1];
                    for (unsigned int j = 0; j < TILE_SIZE; j++) {
                        real3 pos2 = make_real3(localData[j].x, localData[j].y, localData[j].z);
                        real charge2 = localData[j].q;
                        real3 delta = make_real3(pos2.x-posq1.x, pos2.y-posq1.y, pos2.z-posq1.z);
                        real r2 = delta.x*delta.x + delta.y*delta.y + delta.z*delta.z;
                        int atom2 = atomIndices[j];
                        if (atom1 < NUM_ATOMS && atom2 < NUM_ATOMS && r2 < CUTOFF_SQUARED) {
                            real invR = RSQRT(r2);
                            real r = r2*invR;
                            real bornRadius2 = localData[j].bornRadius;
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
                            energy += tempEnergy;
                            delta *= dEdR;
                            force.x -= delta.x;
                            force.y -= delta.y;
                            force.z -= delta.z;
                            localData[j].fx += delta.x;
                            localData[j].fy += delta.y;
                            localData[j].fz += delta.z;
                            localData[j].fw += dGpol_dalpha2_ij*bornRadius1;
                        }
                    }

                    // Write results for atom1.

#ifdef SUPPORTS_64_BIT_ATOMICS
                    ATOMIC_ADD(&forceBuffers[atom1], (mm_long) (force.x*0x100000000));
                    ATOMIC_ADD(&forceBuffers[atom1+PADDED_NUM_ATOMS], (mm_long) (force.y*0x100000000));
                    ATOMIC_ADD(&forceBuffers[atom1+2*PADDED_NUM_ATOMS], (mm_long) (force.z*0x100000000));
                    ATOMIC_ADD(&global_bornForce[atom1], (mm_long) (force.w*0x100000000));
#else
                    unsigned int offset = atom1 + GROUP_ID*PADDED_NUM_ATOMS;
                    forceBuffers[offset] += make_real4(force.x, force.y, force.z, 0);
                    global_bornForce[offset] += force.w;
#endif
                }
            }
            else
#endif
            {
                // We need to apply periodic boundary conditions separately for each interaction.

                for (unsigned int tgx = 0; tgx < TILE_SIZE; tgx++) {
                    unsigned int atom1 = x*TILE_SIZE+tgx;
                    real4 force = make_real4(0);
                    real4 posq1 = posq[atom1];
                    real charge1 = charge[atom1];
                    float bornRadius1 = global_bornRadii[atom1];
                    for (unsigned int j = 0; j < TILE_SIZE; j++) {
                        real3 pos2 = make_real3(localData[j].x, localData[j].y, localData[j].z);
                        real charge2 = localData[j].q;
                        real3 delta = make_real3(pos2.x-posq1.x, pos2.y-posq1.y, pos2.z-posq1.z);
#ifdef USE_PERIODIC
                        APPLY_PERIODIC_TO_DELTA(delta)
#endif
                        real r2 = delta.x*delta.x + delta.y*delta.y + delta.z*delta.z;
                        int atom2 = atomIndices[j];
#ifdef USE_CUTOFF
                        if (atom1 < NUM_ATOMS && atom2 < NUM_ATOMS && r2 < CUTOFF_SQUARED) {
#else
                        if (atom1 < NUM_ATOMS && atom2 < NUM_ATOMS) {
#endif
                            real invR = RSQRT(r2);
                            real r = r2*invR;
                            real bornRadius2 = localData[j].bornRadius;
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
                            energy += tempEnergy;
                            delta *= dEdR;
                            force.x -= delta.x;
                            force.y -= delta.y;
                            force.z -= delta.z;
                            localData[j].fx += delta.x;
                            localData[j].fy += delta.y;
                            localData[j].fz += delta.z;
                            localData[j].fw += dGpol_dalpha2_ij*bornRadius1;
                        }
                    }

                    // Write results for atom1.

#ifdef SUPPORTS_64_BIT_ATOMICS
                    ATOMIC_ADD(&forceBuffers[atom1], (mm_long) (force.x*0x100000000));
                    ATOMIC_ADD(&forceBuffers[atom1+PADDED_NUM_ATOMS], (mm_long) (force.y*0x100000000));
                    ATOMIC_ADD(&forceBuffers[atom1+2*PADDED_NUM_ATOMS], (mm_long) (force.z*0x100000000));
                    ATOMIC_ADD(&global_bornForce[atom1], (mm_long) (force.w*0x100000000));
#else
                    unsigned int offset = atom1 + GROUP_ID*PADDED_NUM_ATOMS;
                    forceBuffers[offset] += make_real4(force.x, force.y, force.z, 0);
                    global_bornForce[offset] += force.w;
#endif
                }
            }

            // Write results.

            for (int tgx = 0; tgx < TILE_SIZE; tgx++) {
#ifdef USE_CUTOFF
                unsigned int atom2 = atomIndices[tgx];
#else
                unsigned int atom2 = y*TILE_SIZE + tgx;
#endif
                if (atom2 < PADDED_NUM_ATOMS) {
#ifdef SUPPORTS_64_BIT_ATOMICS
                    ATOMIC_ADD(&forceBuffers[atom2], (mm_long) (localData[tgx].fx*0x100000000));
                    ATOMIC_ADD(&forceBuffers[atom2+PADDED_NUM_ATOMS], (mm_long) (localData[tgx].fy*0x100000000));
                    ATOMIC_ADD(&forceBuffers[atom2+2*PADDED_NUM_ATOMS], (mm_long) (localData[tgx].fz*0x100000000));
                    ATOMIC_ADD(&global_bornForce[atom2], (mm_long) (localData[tgx].fw*0x100000000));
#else
                    unsigned int offset = atom2 + GROUP_ID*PADDED_NUM_ATOMS;
                    real4 f = forceBuffers[offset];
                    f.x += localData[tgx].fx;
                    f.y += localData[tgx].fy;
                    f.z += localData[tgx].fz;
                    forceBuffers[offset] = f;
                    global_bornForce[offset] += localData[tgx].fw;
#endif
                }
            }
        }
        pos++;
    }
    energyBuffer[GLOBAL_ID] += energy;
}
