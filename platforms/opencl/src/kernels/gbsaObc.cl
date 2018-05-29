#ifdef SUPPORTS_64_BIT_ATOMICS
#pragma OPENCL EXTENSION cl_khr_int64_base_atomics : enable
#endif
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
__kernel void computeBornSum(
#ifdef SUPPORTS_64_BIT_ATOMICS
        __global long* restrict global_bornSum,
#else
        __global real* restrict global_bornSum,
#endif
        __global const real4* restrict posq, __global const real* restrict charge, __global const float2* restrict global_params,
#ifdef USE_CUTOFF
        __global const int* restrict tiles, __global const unsigned int* restrict interactionCount, real4 periodicBoxSize, real4 invPeriodicBoxSize, 
        real4 periodicBoxVecX, real4 periodicBoxVecY, real4 periodicBoxVecZ, unsigned int maxTiles, __global const real4* restrict blockCenter,
        __global const real4* restrict blockSize, __global const int* restrict interactingAtoms,
#else
        unsigned int numTiles,
#endif
        __global const ushort2* exclusionTiles) {
    const unsigned int totalWarps = get_global_size(0)/TILE_SIZE;
    const unsigned int warp = get_global_id(0)/TILE_SIZE;
    const unsigned int tgx = get_local_id(0) & (TILE_SIZE-1);
    const unsigned int tbx = get_local_id(0) - tgx;
    __local AtomData1 localData[FORCE_WORK_GROUP_SIZE];

    // First loop: process tiles that contain exclusions.
    
    const unsigned int firstExclusionTile = FIRST_EXCLUSION_TILE+warp*(LAST_EXCLUSION_TILE-FIRST_EXCLUSION_TILE)/totalWarps;
    const unsigned int lastExclusionTile = FIRST_EXCLUSION_TILE+(warp+1)*(LAST_EXCLUSION_TILE-FIRST_EXCLUSION_TILE)/totalWarps;
    for (int pos = firstExclusionTile; pos < lastExclusionTile; pos++) {
        const ushort2 tileIndices = exclusionTiles[pos];
        const unsigned int x = tileIndices.x;
        const unsigned int y = tileIndices.y;
        real bornSum = 0.0f;
        unsigned int atom1 = x*TILE_SIZE + tgx;
        real4 posq1 = posq[atom1];
        real charge1 = charge[atom1];
        float2 params1 = global_params[atom1];
        if (x == y) {
            // This tile is on the diagonal.

            localData[get_local_id(0)].x = posq1.x;
            localData[get_local_id(0)].y = posq1.y;
            localData[get_local_id(0)].z = posq1.z;
            localData[get_local_id(0)].q = charge1;
            localData[get_local_id(0)].radius = params1.x;
            localData[get_local_id(0)].scaledRadius = params1.y;
            SYNC_WARPS;
            for (unsigned int j = 0; j < TILE_SIZE; j++) {
                real4 delta = (real4) (localData[tbx+j].x-posq1.x, localData[tbx+j].y-posq1.y, localData[tbx+j].z-posq1.z, 0);
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
                    float2 params2 = (float2) (localData[tbx+j].radius, localData[tbx+j].scaledRadius);
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
            localData[get_local_id(0)].x = tempPosq.x;
            localData[get_local_id(0)].y = tempPosq.y;
            localData[get_local_id(0)].z = tempPosq.z;
            localData[get_local_id(0)].q = charge[j];
            float2 tempParams = global_params[j];
            localData[get_local_id(0)].radius = tempParams.x;
            localData[get_local_id(0)].scaledRadius = tempParams.y;
            localData[get_local_id(0)].bornSum = 0.0f;
            SYNC_WARPS;

            // Compute the full set of interactions in this tile.

            unsigned int tj = tgx;
            for (j = 0; j < TILE_SIZE; j++) {
                real4 delta = (real4) (localData[tbx+tj].x-posq1.x, localData[tbx+tj].y-posq1.y, localData[tbx+tj].z-posq1.z, 0);
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
                    float2 params2 = (float2) (localData[tbx+tj].radius, localData[tbx+tj].scaledRadius);
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
        atom_add(&global_bornSum[offset], (long) (bornSum*0x100000000));
        if (x != y) {
            offset = y*TILE_SIZE + tgx;
            atom_add(&global_bornSum[offset], (long) (localData[get_local_id(0)].bornSum*0x100000000));
        }
#else
        unsigned int offset1 = x*TILE_SIZE + tgx + warp*PADDED_NUM_ATOMS;
        unsigned int offset2 = y*TILE_SIZE + tgx + warp*PADDED_NUM_ATOMS;
        global_bornSum[offset1] += bornSum;
        if (x != y)
            global_bornSum[offset2] += localData[get_local_id(0)].bornSum;
#endif
    }

    // Second loop: tiles without exclusions, either from the neighbor list (with cutoff) or just enumerating all
    // of them (no cutoff).

#ifdef USE_CUTOFF
    unsigned int numTiles = interactionCount[0];
    if (numTiles > maxTiles)
        return; // There wasn't enough memory for the neighbor list.
    int pos = (int) (warp*(numTiles > maxTiles ? NUM_BLOCKS*((long)NUM_BLOCKS+1)/2 : (long)numTiles)/totalWarps);
    int end = (int) ((warp+1)*(numTiles > maxTiles ? NUM_BLOCKS*((long)NUM_BLOCKS+1)/2 : (long)numTiles)/totalWarps);
#else
    int pos = (int) (warp*(long)numTiles/totalWarps);
    int end = (int) ((warp+1)*(long)numTiles/totalWarps);
#endif
    int skipBase = 0;
    int currentSkipIndex = tbx;
    __local int atomIndices[FORCE_WORK_GROUP_SIZE];
    __local volatile int skipTiles[FORCE_WORK_GROUP_SIZE];
    skipTiles[get_local_id(0)] = -1;

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
                skipTiles[get_local_id(0)] = tile.x + tile.y*NUM_BLOCKS - tile.y*(tile.y+1)/2;
            }
            else
                skipTiles[get_local_id(0)] = end;
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
            atomIndices[get_local_id(0)] = j;
            if (j < PADDED_NUM_ATOMS) {
                real4 tempPosq = posq[j];
                localData[get_local_id(0)].x = tempPosq.x;
                localData[get_local_id(0)].y = tempPosq.y;
                localData[get_local_id(0)].z = tempPosq.z;
                localData[get_local_id(0)].q = charge[j];
                float2 tempParams = global_params[j];
                localData[get_local_id(0)].radius = tempParams.x;
                localData[get_local_id(0)].scaledRadius = tempParams.y;
                localData[get_local_id(0)].bornSum = 0.0f;
            }
            SYNC_WARPS;
#ifdef USE_PERIODIC
            if (singlePeriodicCopy) {
                // The box is small enough that we can just translate all the atoms into a single periodic
                // box, then skip having to apply periodic boundary conditions later.

                real4 blockCenterX = blockCenter[x];
                APPLY_PERIODIC_TO_POS_WITH_CENTER(posq1, blockCenterX)
                APPLY_PERIODIC_TO_POS_WITH_CENTER(localData[get_local_id(0)], blockCenterX)
                SYNC_WARPS;
                unsigned int tj = tgx;
                for (j = 0; j < TILE_SIZE; j++) {
                    real4 delta = (real4) (localData[tbx+tj].x-posq1.x, localData[tbx+tj].y-posq1.y, localData[tbx+tj].z-posq1.z, 0);
                    real r2 = delta.x*delta.x + delta.y*delta.y + delta.z*delta.z;
                    int atom2 = atomIndices[tbx+tj];
                    if (atom1 < NUM_ATOMS && atom2 < NUM_ATOMS && r2 < CUTOFF_SQUARED) {
                        real invR = RSQRT(r2);
                        real r = r2*invR;
                        float2 params2 = (float2) (localData[tbx+tj].radius, localData[tbx+tj].scaledRadius);
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
                    real4 delta = (real4) (localData[tbx+tj].x-posq1.x, localData[tbx+tj].y-posq1.y, localData[tbx+tj].z-posq1.z, 0);
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
                        float2 params2 = (float2) (localData[tbx+tj].radius, localData[tbx+tj].scaledRadius);
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
            unsigned int atom2 = atomIndices[get_local_id(0)];
#else
            unsigned int atom2 = y*TILE_SIZE + tgx;
#endif
#ifdef SUPPORTS_64_BIT_ATOMICS
            atom_add(&global_bornSum[atom1], (long) (bornSum*0x100000000));
            if (atom2 < PADDED_NUM_ATOMS)
                atom_add(&global_bornSum[atom2], (long) (localData[get_local_id(0)].bornSum*0x100000000));
#else
            unsigned int offset1 = atom1 + warp*PADDED_NUM_ATOMS;
            unsigned int offset2 = atom2 + warp*PADDED_NUM_ATOMS;
            global_bornSum[offset1] += bornSum;
            if (atom2 < PADDED_NUM_ATOMS)
                global_bornSum[offset2] += localData[get_local_id(0)].bornSum;
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

__kernel void computeGBSAForce1(
#ifdef SUPPORTS_64_BIT_ATOMICS
        __global long* restrict forceBuffers, __global long* restrict global_bornForce,
#else
        __global real4* restrict forceBuffers, __global real* restrict global_bornForce,
#endif
        __global mixed* restrict energyBuffer, __global const real4* restrict posq, __global const real* restrict charge,
        __global const real* restrict global_bornRadii, int needEnergy,
#ifdef USE_CUTOFF
        __global const int* restrict tiles, __global const unsigned int* restrict interactionCount, real4 periodicBoxSize, real4 invPeriodicBoxSize, 
        real4 periodicBoxVecX, real4 periodicBoxVecY, real4 periodicBoxVecZ, unsigned int maxTiles, __global const real4* restrict blockCenter,
        __global const real4* restrict blockSize, __global const int* restrict interactingAtoms,
#else
        unsigned int numTiles,
#endif
        __global const ushort2* exclusionTiles) {
    const unsigned int totalWarps = get_global_size(0)/TILE_SIZE;
    const unsigned int warp = get_global_id(0)/TILE_SIZE;
    const unsigned int tgx = get_local_id(0) & (TILE_SIZE-1);
    const unsigned int tbx = get_local_id(0) - tgx;
    mixed energy = 0;
    __local AtomData2 localData[FORCE_WORK_GROUP_SIZE];

    // First loop: process tiles that contain exclusions.
    
    const unsigned int firstExclusionTile = FIRST_EXCLUSION_TILE+warp*(LAST_EXCLUSION_TILE-FIRST_EXCLUSION_TILE)/totalWarps;
    const unsigned int lastExclusionTile = FIRST_EXCLUSION_TILE+(warp+1)*(LAST_EXCLUSION_TILE-FIRST_EXCLUSION_TILE)/totalWarps;
    for (int pos = firstExclusionTile; pos < lastExclusionTile; pos++) {
        const ushort2 tileIndices = exclusionTiles[pos];
        const unsigned int x = tileIndices.x;
        const unsigned int y = tileIndices.y;
        real4 force = 0.0f;
        unsigned int atom1 = x*TILE_SIZE + tgx;
        real4 posq1 = posq[atom1];
        real charge1 = charge[atom1];
        real bornRadius1 = global_bornRadii[atom1];
        if (x == y) {
            // This tile is on the diagonal.

            const unsigned int localAtomIndex = get_local_id(0);
            localData[localAtomIndex].x = posq1.x;
            localData[localAtomIndex].y = posq1.y;
            localData[localAtomIndex].z = posq1.z;
            localData[localAtomIndex].q = charge1;
            localData[get_local_id(0)].bornRadius = bornRadius1;
            SYNC_WARPS;
            for (unsigned int j = 0; j < TILE_SIZE; j++) {
                if (atom1 < NUM_ATOMS && y*TILE_SIZE+j < NUM_ATOMS) {
                    real3 pos2 = (real3) (localData[tbx+j].x, localData[tbx+j].y, localData[tbx+j].z);
                    real charge2 = localData[tbx+j].q;
                    real4 delta = (real4) (pos2 - posq1.xyz, 0);
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
                        delta.xyz *= dEdR;
                        force.xyz -= delta.xyz;
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
            localData[get_local_id(0)].x = tempPosq.x;
            localData[get_local_id(0)].y = tempPosq.y;
            localData[get_local_id(0)].z = tempPosq.z;
            localData[get_local_id(0)].q = charge[j];
            localData[get_local_id(0)].bornRadius = global_bornRadii[j];
            localData[get_local_id(0)].fx = 0.0f;
            localData[get_local_id(0)].fy = 0.0f;
            localData[get_local_id(0)].fz = 0.0f;
            localData[get_local_id(0)].fw = 0.0f;
            SYNC_WARPS;
            unsigned int tj = tgx;
            for (j = 0; j < TILE_SIZE; j++) {
                if (atom1 < NUM_ATOMS && y*TILE_SIZE+tj < NUM_ATOMS) {
                    real3 pos2 = (real3) (localData[tbx+tj].x, localData[tbx+tj].y, localData[tbx+tj].z);
                    real charge2 = localData[tbx+tj].q;
                    real4 delta = (real4) (pos2 - posq1.xyz, 0);
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
                        delta.xyz *= dEdR;
                        force.xyz -= delta.xyz;
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
        atom_add(&forceBuffers[offset], (long) (force.x*0x100000000));
        atom_add(&forceBuffers[offset+PADDED_NUM_ATOMS], (long) (force.y*0x100000000));
        atom_add(&forceBuffers[offset+2*PADDED_NUM_ATOMS], (long) (force.z*0x100000000));
        atom_add(&global_bornForce[offset], (long) (force.w*0x100000000));
        if (x != y) {
            offset = y*TILE_SIZE + tgx;
            atom_add(&forceBuffers[offset], (long) (localData[get_local_id(0)].fx*0x100000000));
            atom_add(&forceBuffers[offset+PADDED_NUM_ATOMS], (long) (localData[get_local_id(0)].fy*0x100000000));
            atom_add(&forceBuffers[offset+2*PADDED_NUM_ATOMS], (long) (localData[get_local_id(0)].fz*0x100000000));
            atom_add(&global_bornForce[offset], (long) (localData[get_local_id(0)].fw*0x100000000));
        }
#else
        unsigned int offset1 = x*TILE_SIZE + tgx + warp*PADDED_NUM_ATOMS;
        unsigned int offset2 = y*TILE_SIZE + tgx + warp*PADDED_NUM_ATOMS;
        forceBuffers[offset1].xyz += force.xyz;
        global_bornForce[offset1] += force.w;
        if (x != y) {
            forceBuffers[offset2] += (real4) (localData[get_local_id(0)].fx, localData[get_local_id(0)].fy, localData[get_local_id(0)].fz, 0.0f);
            global_bornForce[offset2] += localData[get_local_id(0)].fw;
        }
#endif
    }

    // Second loop: tiles without exclusions, either from the neighbor list (with cutoff) or just enumerating all
    // of them (no cutoff).

#ifdef USE_CUTOFF
    unsigned int numTiles = interactionCount[0];
    if (numTiles > maxTiles)
        return; // There wasn't enough memory for the neighbor list.
    int pos = (int) (warp*(numTiles > maxTiles ? NUM_BLOCKS*((long)NUM_BLOCKS+1)/2 : (long)numTiles)/totalWarps);
    int end = (int) ((warp+1)*(numTiles > maxTiles ? NUM_BLOCKS*((long)NUM_BLOCKS+1)/2 : (long)numTiles)/totalWarps);
#else
    int pos = (int) (warp*(long)numTiles/totalWarps);
    int end = (int) ((warp+1)*(long)numTiles/totalWarps);
#endif
    int skipBase = 0;
    int currentSkipIndex = tbx;
    __local int atomIndices[FORCE_WORK_GROUP_SIZE];
    __local volatile int skipTiles[FORCE_WORK_GROUP_SIZE];
    skipTiles[get_local_id(0)] = -1;

    while (pos < end) {
        real4 force = 0;
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
                skipTiles[get_local_id(0)] = tile.x + tile.y*NUM_BLOCKS - tile.y*(tile.y+1)/2;
            }
            else
                skipTiles[get_local_id(0)] = end;
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
            atomIndices[get_local_id(0)] = j;
            if (j < PADDED_NUM_ATOMS) {
                real4 tempPosq = posq[j];
                localData[get_local_id(0)].x = tempPosq.x;
                localData[get_local_id(0)].y = tempPosq.y;
                localData[get_local_id(0)].z = tempPosq.z;
                localData[get_local_id(0)].q = charge[j];
                localData[get_local_id(0)].bornRadius = global_bornRadii[j];
                localData[get_local_id(0)].fx = 0.0f;
                localData[get_local_id(0)].fy = 0.0f;
                localData[get_local_id(0)].fz = 0.0f;
                localData[get_local_id(0)].fw = 0.0f;
            }
            SYNC_WARPS;
#ifdef USE_PERIODIC
            if (singlePeriodicCopy) {
                // The box is small enough that we can just translate all the atoms into a single periodic
                // box, then skip having to apply periodic boundary conditions later.

                real4 blockCenterX = blockCenter[x];
                APPLY_PERIODIC_TO_POS_WITH_CENTER(posq1, blockCenterX)
                APPLY_PERIODIC_TO_POS_WITH_CENTER(localData[get_local_id(0)], blockCenterX)
                SYNC_WARPS;
                unsigned int tj = tgx;
                for (j = 0; j < TILE_SIZE; j++) {
                    int atom2 = atomIndices[tbx+tj];
                    if (atom1 < NUM_ATOMS && atom2 < NUM_ATOMS) {
                        real3 pos2 = (real3) (localData[tbx+tj].x, localData[tbx+tj].y, localData[tbx+tj].z);
                        real charge2 = localData[tbx+tj].q;
                        real4 delta = (real4) (pos2 - posq1.xyz, 0);
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
                            delta.xyz *= dEdR;
                            force.xyz -= delta.xyz;
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
                        real3 pos2 = (real3) (localData[tbx+tj].x, localData[tbx+tj].y, localData[tbx+tj].z);
                        real charge2 = localData[tbx+tj].q;
                        real4 delta = (real4) (pos2 - posq1.xyz, 0);
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
                            delta.xyz *= dEdR;
                            force.xyz -= delta.xyz;
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
            unsigned int atom2 = atomIndices[get_local_id(0)];
#else
            unsigned int atom2 = y*TILE_SIZE + tgx;
#endif
#ifdef SUPPORTS_64_BIT_ATOMICS
            atom_add(&forceBuffers[atom1], (long) (force.x*0x100000000));
            atom_add(&forceBuffers[atom1+PADDED_NUM_ATOMS], (long) (force.y*0x100000000));
            atom_add(&forceBuffers[atom1+2*PADDED_NUM_ATOMS], (long) (force.z*0x100000000));
            atom_add(&global_bornForce[atom1], (long) (force.w*0x100000000));
            if (atom2 < PADDED_NUM_ATOMS) {
                atom_add(&forceBuffers[atom2], (long) (localData[get_local_id(0)].fx*0x100000000));
                atom_add(&forceBuffers[atom2+PADDED_NUM_ATOMS], (long) (localData[get_local_id(0)].fy*0x100000000));
                atom_add(&forceBuffers[atom2+2*PADDED_NUM_ATOMS], (long) (localData[get_local_id(0)].fz*0x100000000));
                atom_add(&global_bornForce[atom2], (long) (localData[get_local_id(0)].fw*0x100000000));
            }
#else
            unsigned int offset1 = atom1 + warp*PADDED_NUM_ATOMS;
            unsigned int offset2 = atom2 + warp*PADDED_NUM_ATOMS;
            forceBuffers[offset1].xyz += force.xyz;
            global_bornForce[offset1] += force.w;
            if (atom2 < PADDED_NUM_ATOMS) {
                forceBuffers[offset2] += (real4) (localData[get_local_id(0)].fx, localData[get_local_id(0)].fy, localData[get_local_id(0)].fz, 0.0f);
                global_bornForce[offset2] += localData[get_local_id(0)].fw;
            }
#endif
        }
        pos++;
    }
    energyBuffer[get_global_id(0)] += energy;
}
