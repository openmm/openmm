DEVICE void computeDirectFieldDampingFactors(real alpha, real r, real* fdamp3, real* fdamp5, real* fdamp7) {
    real ar = alpha*r;
    real ar2 = ar*ar;
    real ar3 = ar2*ar;
    real ar4 = ar2*ar2;
    real expAR = EXP(-ar);
    real one = 1;
    *fdamp3 = 1 - (1 + ar + ar2*(one/2))*expAR;
    *fdamp5 = 1 - (1 + ar + ar2*(one/2) + ar3*(one/6))*expAR;
    *fdamp7 = 1 - (1 + ar + ar2*(one/2) + ar3*(one/6) + ar4*(one/30))*expAR;
}

DEVICE void computeMutualFieldDampingFactors(real alphaI, real alphaJ, real r, real* fdamp3, real* fdamp5) {
    real arI = alphaI*r;
    real arI2 = arI*arI;
    real arI3 = arI2*arI;
    real expARI = EXP(-arI);
    real one = 1;
    real seven = 7;
    if (alphaI == alphaJ) {
        real arI4 = arI3*arI;
        real arI5 = arI4*arI;
        *fdamp3 = 1 - (1 + arI + arI2*(one/2) + arI3*(seven/48) + arI4*(one/48))*expARI;
        *fdamp5 = 1 - (1 + arI + arI2*(one/2) + arI3*(one/6) + arI4*(one/24) + arI5*(one/144))*expARI;
    }
    else {
        real arJ = alphaJ*r;
        real arJ2 = arJ*arJ;
        real arJ3 = arJ2*arJ;
        real expARJ = EXP(-arJ);
        real aI2 = alphaI*alphaI;
        real aJ2 = alphaJ*alphaJ;
        real A = aJ2/(aJ2-aI2);
        real B = aI2/(aI2-aJ2);
        real A2expARI = A*A*expARI;
        real B2expARJ = B*B*expARJ;
        *fdamp3 = 1 - (1 + arI + arI2*(one/2))*A2expARI -
                      (1 + arJ + arJ2*(one/2))*B2expARJ -
                      (1 + arI)*2*B*A2expARI -
                      (1 + arJ)*2*A*B2expARJ;
        *fdamp5 = 1 - (1 + arI + arI2*(one/2) + arI3*(one/6))*A2expARI -
                      (1 + arJ + arJ2*(one/2) + arJ3*(one/6))*B2expARJ -
                      (1 + arI + arI2*(one/3))*2*B*A2expARI -
                      (1 + arJ + arJ2*(one/3))*2*A*B2expARJ;
    }
}

typedef struct {
    real x, y, z;
    real fx, fy, fz;
    ATOM_PARAMETER_DATA
} AtomData;

/**
 * Compute the electrostatic field.
 */
KERNEL void computeField(GLOBAL const real4* RESTRICT posq, GLOBAL const unsigned int* RESTRICT exclusions,
        GLOBAL const int2* RESTRICT exclusionTiles, GLOBAL mm_ulong* RESTRICT fieldBuffers,
#ifdef USE_CUTOFF
        GLOBAL const int* RESTRICT tiles, GLOBAL const unsigned int* RESTRICT interactionCount, real4 periodicBoxSize, real4 invPeriodicBoxSize,
        real4 periodicBoxVecX, real4 periodicBoxVecY, real4 periodicBoxVecZ, unsigned int maxTiles, GLOBAL const real4* RESTRICT blockCenter,
        GLOBAL const real4* RESTRICT blockSize, GLOBAL const unsigned int* RESTRICT interactingAtoms
#else
        unsigned int numTiles
#endif
        PARAMETER_ARGUMENTS) {
    const unsigned int totalWarps = (GLOBAL_SIZE)/TILE_SIZE;
    const unsigned int warp = (GLOBAL_ID)/TILE_SIZE;
    const unsigned int tgx = LOCAL_ID & (TILE_SIZE-1);
    const unsigned int tbx = LOCAL_ID - tgx;
    LOCAL AtomData localData[THREAD_BLOCK_SIZE];

    // First loop: process tiles that contain exclusions.
    
    const unsigned int firstExclusionTile = warp*NUM_TILES_WITH_EXCLUSIONS/totalWarps;
    const unsigned int lastExclusionTile = (warp+1)*NUM_TILES_WITH_EXCLUSIONS/totalWarps;
    for (int tile = firstExclusionTile; tile < lastExclusionTile; tile++) {
        const int2 tileIndices = exclusionTiles[tile];
        const unsigned int x = tileIndices.x;
        const unsigned int y = tileIndices.y;
        real3 field = make_real3(0);
        unsigned int atom1 = x*TILE_SIZE + tgx;
        real4 pos1 = posq[atom1];
        LOAD_ATOM1_PARAMETERS
        unsigned int excl = exclusions[tile*TILE_SIZE+tgx];
        if (x == y) {
            // This tile is on the diagonal.

            const unsigned int localAtomIndex = LOCAL_ID;
            localData[LOCAL_ID].x = pos1.x;
            localData[LOCAL_ID].y = pos1.y;
            localData[LOCAL_ID].z = pos1.z;
            LOAD_LOCAL_PARAMETERS_FROM_1
            SYNC_WARPS;
            for (unsigned int j = 0; j < TILE_SIZE; j++) {
                int atom2 = tbx+j;
                real3 pos2 = make_real3(localData[atom2].x, localData[atom2].y, localData[atom2].z);
                real3 delta = make_real3(pos2.x-pos1.x, pos2.y-pos1.y, pos2.z-pos1.z);
#ifdef USE_PERIODIC
                APPLY_PERIODIC_TO_DELTA(delta)
#endif
                real r2 = delta.x*delta.x + delta.y*delta.y + delta.z*delta.z;
#ifdef USE_CUTOFF
                if (r2 < CUTOFF_SQUARED) {
#endif
                    real invR = RSQRT(r2);
                    real r = r2*invR;
                    LOAD_ATOM2_PARAMETERS
                    atom2 = y*TILE_SIZE+j;
                    real3 tempField1 = make_real3(0);
                    real3 tempField2 = make_real3(0);
                    bool isExcluded = (atom1 >= NUM_ATOMS || atom2 >= NUM_ATOMS || !(excl & 0x1));
                    if (!isExcluded && atom1 != atom2) {
                        COMPUTE_FIELD
                    }
                    field += tempField1;
#ifdef USE_CUTOFF
                }
#endif
                excl >>= 1;
                SYNC_WARPS;
            }
        }
        else {
            // This is an off-diagonal tile.

            const unsigned int localAtomIndex = LOCAL_ID;
            unsigned int j = y*TILE_SIZE + tgx;
            real4 tempPosq = posq[j];
            localData[LOCAL_ID].x = tempPosq.x;
            localData[LOCAL_ID].y = tempPosq.y;
            localData[LOCAL_ID].z = tempPosq.z;
            LOAD_LOCAL_PARAMETERS_FROM_GLOBAL
            localData[localAtomIndex].fx = 0;
            localData[localAtomIndex].fy = 0;
            localData[localAtomIndex].fz = 0;
            SYNC_WARPS;
            excl = (excl >> tgx) | (excl << (TILE_SIZE - tgx));
            unsigned int tj = tgx;
            for (j = 0; j < TILE_SIZE; j++) {
                int atom2 = tbx+tj;
                real3 pos2 = make_real3(localData[atom2].x, localData[atom2].y, localData[atom2].z);
                real3 delta = make_real3(pos2.x-pos1.x, pos2.y-pos1.y, pos2.z-pos1.z);
#ifdef USE_PERIODIC
                APPLY_PERIODIC_TO_DELTA(delta)
#endif
                real r2 = delta.x*delta.x + delta.y*delta.y + delta.z*delta.z;
#ifdef USE_CUTOFF
                if (r2 < CUTOFF_SQUARED) {
#endif
                    real invR = RSQRT(r2);
                    real r = r2*invR;
                    LOAD_ATOM2_PARAMETERS
                    atom2 = y*TILE_SIZE+tj;
                    real3 tempField1 = make_real3(0);
                    real3 tempField2 = make_real3(0);
                    bool isExcluded = (atom1 >= NUM_ATOMS || atom2 >= NUM_ATOMS || !(excl & 0x1));
                    if (!isExcluded) {
                        COMPUTE_FIELD
                    }
                    field += tempField1;
                    localData[tbx+tj].fx += tempField2.x;
                    localData[tbx+tj].fy += tempField2.y;
                    localData[tbx+tj].fz += tempField2.z;
#ifdef USE_CUTOFF
                }
#endif
                excl >>= 1;
                tj = (tj + 1) & (TILE_SIZE - 1);
                SYNC_WARPS;
            }
        }

        // Write results.

        unsigned int offset1 = x*TILE_SIZE + tgx;
        ATOMIC_ADD(&fieldBuffers[offset1], (mm_ulong) ((mm_long) (field.x*0x100000000)));
        ATOMIC_ADD(&fieldBuffers[offset1+PADDED_NUM_ATOMS], (mm_ulong) ((mm_long) (field.y*0x100000000)));
        ATOMIC_ADD(&fieldBuffers[offset1+2*PADDED_NUM_ATOMS], (mm_ulong) ((mm_long) (field.z*0x100000000)));
        if (x != y) {
            unsigned int offset2 = y*TILE_SIZE + tgx;
            ATOMIC_ADD(&fieldBuffers[offset2], (mm_ulong) ((mm_long) (localData[LOCAL_ID].fx*0x100000000)));
            ATOMIC_ADD(&fieldBuffers[offset2+PADDED_NUM_ATOMS], (mm_ulong) ((mm_long) (localData[LOCAL_ID].fy*0x100000000)));
            ATOMIC_ADD(&fieldBuffers[offset2+2*PADDED_NUM_ATOMS], (mm_ulong) ((mm_long) (localData[LOCAL_ID].fz*0x100000000)));
        }
    }

    // Second loop: tiles without exclusions, either from the neighbor list (with cutoff) or just enumerating all
    // of them (no cutoff).

#ifdef USE_CUTOFF
    unsigned int numTiles = interactionCount[0];
    if (numTiles > maxTiles)
        return; // There wasn't enough memory for the neighbor list.
    int tile = (int) (warp*(numTiles > maxTiles ? NUM_BLOCKS*((mm_long)NUM_BLOCKS+1)/2 : (long)numTiles)/totalWarps);
    int end = (int) ((warp+1)*(numTiles > maxTiles ? NUM_BLOCKS*((mm_long)NUM_BLOCKS+1)/2 : (long)numTiles)/totalWarps);
#else
    int tile = (int) (warp*(mm_long)numTiles/totalWarps);
    int end = (int) ((warp+1)*(mm_long)numTiles/totalWarps);
#endif
    int skipBase = 0;
    int currentSkipIndex = tbx;
    LOCAL int atomIndices[THREAD_BLOCK_SIZE];
    LOCAL volatile int skipTiles[THREAD_BLOCK_SIZE];
    skipTiles[LOCAL_ID] = -1;
    
    while (tile < end) {
        real3 field = make_real3(0);
        bool includeTile = true;
        
        // Extract the coordinates of this tile.
        
        int x, y;
        bool singlePeriodicCopy = false;
#ifdef USE_CUTOFF
        x = tiles[tile];
        real4 blockSizeX = blockSize[x];
        singlePeriodicCopy = (0.5f*periodicBoxSize.x-blockSizeX.x >= CUTOFF &&
                              0.5f*periodicBoxSize.y-blockSizeX.y >= CUTOFF &&
                              0.5f*periodicBoxSize.z-blockSizeX.z >= CUTOFF);
#else
        y = (int) floor(NUM_BLOCKS+0.5f-SQRT((NUM_BLOCKS+0.5f)*(NUM_BLOCKS+0.5f)-2*tile));
        x = (tile-y*NUM_BLOCKS+y*(y+1)/2);
        if (x < y || x >= NUM_BLOCKS) { // Occasionally happens due to roundoff error.
            y += (x < y ? -1 : 1);
            x = (tile-y*NUM_BLOCKS+y*(y+1)/2);
        }

        // Skip over tiles that have exclusions, since they were already processed.

        SYNC_WARPS;
        while (skipTiles[tbx+TILE_SIZE-1] < tile) {
            SYNC_WARPS;
            if (skipBase+tgx < NUM_TILES_WITH_EXCLUSIONS) {
                int2 tile = exclusionTiles[skipBase+tgx];
                skipTiles[LOCAL_ID] = tile.x + tile.y*NUM_BLOCKS - tile.y*(tile.y+1)/2;
            }
            else
                skipTiles[LOCAL_ID] = end;
            skipBase += TILE_SIZE;            
            currentSkipIndex = tbx;
            SYNC_WARPS;
        }
        while (skipTiles[currentSkipIndex] < tile)
            currentSkipIndex++;
        includeTile = (skipTiles[currentSkipIndex] != tile);
#endif
        if (includeTile) {
            unsigned int atom1 = x*TILE_SIZE + tgx;

            // Load atom data for this tile.
            
            real4 pos1 = posq[atom1];
            LOAD_ATOM1_PARAMETERS
            const unsigned int localAtomIndex = LOCAL_ID;
#ifdef USE_CUTOFF
            unsigned int j = interactingAtoms[tile*TILE_SIZE+tgx];
#else
            unsigned int j = y*TILE_SIZE + tgx;
#endif
            atomIndices[LOCAL_ID] = j;
            if (j < PADDED_NUM_ATOMS) {
                real4 tempPosq = posq[j];
                localData[localAtomIndex].x = tempPosq.x;
                localData[localAtomIndex].y = tempPosq.y;
                localData[localAtomIndex].z = tempPosq.z;
                LOAD_LOCAL_PARAMETERS_FROM_GLOBAL
                localData[localAtomIndex].fx = 0;
                localData[localAtomIndex].fy = 0;
                localData[localAtomIndex].fz = 0;
            }
            SYNC_WARPS;
#ifdef USE_PERIODIC
            if (singlePeriodicCopy) {
                // The box is small enough that we can just translate all the atoms into a single periodic
                // box, then skip having to apply periodic boundary conditions later.

                real4 blockCenterX = blockCenter[x];
                APPLY_PERIODIC_TO_POS_WITH_CENTER(pos1, blockCenterX)
                APPLY_PERIODIC_TO_POS_WITH_CENTER(localData[LOCAL_ID], blockCenterX)
                SYNC_WARPS;
                unsigned int tj = tgx;
                for (unsigned int j = 0; j < TILE_SIZE; j++) {
                    int atom2 = tbx+tj;
                    real3 pos2 = make_real3(localData[atom2].x, localData[atom2].y, localData[atom2].z);
                    real3 delta = make_real3(pos2.x-pos1.x, pos2.y-pos1.y, pos2.z-pos1.z);
                    real r2 = delta.x*delta.x + delta.y*delta.y + delta.z*delta.z;
                    if (r2 < CUTOFF_SQUARED) {
                        real invR = RSQRT(r2);
                        real r = r2*invR;
                        LOAD_ATOM2_PARAMETERS
                        atom2 = atomIndices[tbx+tj];
                        real3 tempField1 = make_real3(0);
                        real3 tempField2 = make_real3(0);
                        if (atom1 < NUM_ATOMS && atom2 < NUM_ATOMS) {
                            COMPUTE_FIELD
                        }
                        field += tempField1;
                        localData[tbx+tj].fx += tempField2.x;
                        localData[tbx+tj].fy += tempField2.y;
                        localData[tbx+tj].fz += tempField2.z;
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
                for (unsigned int j = 0; j < TILE_SIZE; j++) {
                    int atom2 = tbx+tj;
                    real3 pos2 = make_real3(localData[atom2].x, localData[atom2].y, localData[atom2].z);
                    real3 delta = make_real3(pos2.x-pos1.x, pos2.y-pos1.y, pos2.z-pos1.z);
#ifdef USE_PERIODIC
                    APPLY_PERIODIC_TO_DELTA(delta)
#endif
                    real r2 = delta.x*delta.x + delta.y*delta.y + delta.z*delta.z;
#ifdef USE_CUTOFF
                    if (r2 < CUTOFF_SQUARED) {
#endif
                        real invR = RSQRT(r2);
                        real r = r2*invR;
                        LOAD_ATOM2_PARAMETERS
                        atom2 = atomIndices[tbx+tj];
                        real3 tempField1 = make_real3(0);
                        real3 tempField2 = make_real3(0);
                        if (atom1 < NUM_ATOMS && atom2 < NUM_ATOMS) {
                            COMPUTE_FIELD
                        }
                        field += tempField1;
                        localData[tbx+tj].fx += tempField2.x;
                        localData[tbx+tj].fy += tempField2.y;
                        localData[tbx+tj].fz += tempField2.z;
#ifdef USE_CUTOFF
                    }
#endif
                    tj = (tj + 1) & (TILE_SIZE - 1);
                    SYNC_WARPS;
                }
            }
        
            // Write results.

            ATOMIC_ADD(&fieldBuffers[atom1], (mm_ulong) ((mm_long) (field.x*0x100000000)));
            ATOMIC_ADD(&fieldBuffers[atom1+PADDED_NUM_ATOMS], (mm_ulong) ((mm_long) (field.y*0x100000000)));
            ATOMIC_ADD(&fieldBuffers[atom1+2*PADDED_NUM_ATOMS], (mm_ulong) ((mm_long) (field.z*0x100000000)));
#ifdef USE_CUTOFF
            unsigned int atom2 = atomIndices[LOCAL_ID];
#else
            unsigned int atom2 = y*TILE_SIZE + tgx;
#endif
            if (atom2 < PADDED_NUM_ATOMS) {
                ATOMIC_ADD(&fieldBuffers[atom2], (mm_ulong) ((mm_long) (localData[LOCAL_ID].fx*0x100000000)));
                ATOMIC_ADD(&fieldBuffers[atom2+PADDED_NUM_ATOMS], (mm_ulong) ((mm_long) (localData[LOCAL_ID].fy*0x100000000)));
                ATOMIC_ADD(&fieldBuffers[atom2+2*PADDED_NUM_ATOMS], (mm_ulong) ((mm_long) (localData[LOCAL_ID].fz*0x100000000)));
            }
        }
        tile++;
    }
}

#define COMPUTING_EXCEPTIONS

/**
 * Compute the electrostatic field from nonbonded exceptions.
 */
KERNEL void computeFieldExceptions(GLOBAL const real4* RESTRICT posq, GLOBAL mm_ulong* RESTRICT fieldBuffers,
        GLOBAL const int2* RESTRICT exceptionAtoms, GLOBAL const real* RESTRICT exceptionScale
#ifdef USE_CUTOFF
        , real4 periodicBoxSize, real4 invPeriodicBoxSize, real4 periodicBoxVecX, real4 periodicBoxVecY, real4 periodicBoxVecZ
#endif
        PARAMETER_ARGUMENTS) {
    for (int index = GLOBAL_ID; index < NUM_EXCEPTIONS; index += GLOBAL_SIZE) {
        int2 atoms = exceptionAtoms[index];
        int atom1 = atoms.x;
        int atom2 = atoms.y;
        real4 pos1 = posq[atom1];
        real4 pos2 = posq[atom2];
        LOAD_ATOM1_PARAMETERS
        LOAD_ATOM2_PARAMETERS_FROM_GLOBAL
        real scale = exceptionScale[index];
        real3 delta = make_real3(pos2.x-pos1.x, pos2.y-pos1.y, pos2.z-pos1.z);
#ifdef USE_PERIODIC
        APPLY_PERIODIC_TO_DELTA(delta)
#endif
        real r2 = delta.x*delta.x + delta.y*delta.y + delta.z*delta.z;
#ifdef USE_CUTOFF
        if (r2 < CUTOFF_SQUARED) {
#endif
            real invR = RSQRT(r2);
            real r = r2*invR;
            real3 tempField1 = make_real3(0);
            real3 tempField2 = make_real3(0);
            COMPUTE_FIELD
            ATOMIC_ADD(&fieldBuffers[atom1], (mm_ulong) ((mm_long) (tempField1.x*0x100000000)));
            ATOMIC_ADD(&fieldBuffers[atom1+PADDED_NUM_ATOMS], (mm_ulong) ((mm_long) (tempField1.y*0x100000000)));
            ATOMIC_ADD(&fieldBuffers[atom1+2*PADDED_NUM_ATOMS], (mm_ulong) ((mm_long) (tempField1.z*0x100000000)));
            ATOMIC_ADD(&fieldBuffers[atom2], (mm_ulong) ((mm_long) (tempField2.x*0x100000000)));
            ATOMIC_ADD(&fieldBuffers[atom2+PADDED_NUM_ATOMS], (mm_ulong) ((mm_long) (tempField2.y*0x100000000)));
            ATOMIC_ADD(&fieldBuffers[atom2+2*PADDED_NUM_ATOMS], (mm_ulong) ((mm_long) (tempField2.z*0x100000000)));
#ifdef USE_CUTOFF
        }
#endif
    }
}