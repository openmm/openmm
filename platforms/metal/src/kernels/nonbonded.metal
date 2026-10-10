#define WARPS_PER_GROUP (THREAD_BLOCK_SIZE/TILE_SIZE)

/**
 * Save the force on a single atom.
 */
DEVICE void saveSingleForce(int atom, real3 force, GLOBAL mm_ulong* forceBuffers) {
    if (force.x != 0)
        ATOMIC_ADD(&forceBuffers[atom], static_cast<mm_ulong>(realToFixedPoint(force.x)));
    if (force.y != 0)
        ATOMIC_ADD(&forceBuffers[atom+PADDED_NUM_ATOMS], static_cast<mm_ulong>(realToFixedPoint(force.y)));
    if (force.z != 0)
        ATOMIC_ADD(&forceBuffers[atom+2*PADDED_NUM_ATOMS], static_cast<mm_ulong>(realToFixedPoint(force.z)));
}

/**
 * Compute nonbonded interactions. The kernel is separated into two parts,
 * tiles with exclusions and tiles without exclusions. It relies heavily on
 * implicit warp-level synchronization. A tile is defined by two atom blocks
 * each of warpsize. Each warp computes a range of tiles.
 *
 * Tiles with exclusions compute the entire set of interactions across
 * atom blocks, equal to warpsize*warpsize. In order to avoid access conflicts
 * the forces are computed and accumulated diagonally in the manner shown below
 * where, suppose
 *
 * [a-h] comprise atom block 1, [i-p] comprise atom block 2
 *
 * 1 denotes the first set of calculations within the warp
 * 2 denotes the second set of calculations within the warp
 * ... etc.
 *
 *        threads
 *     0 1 2 3 4 5 6 7
 *         atom1
 * L    a b c d e f g h
 * o  i 1 2 3 4 5 6 7 8
 * c  j 8 1 2 3 4 5 6 7
 * a  k 7 8 1 2 3 4 5 6
 * l  l 6 7 8 1 2 3 4 5
 * D  m 5 6 7 8 1 2 3 4
 * a  n 4 5 6 7 8 1 2 3
 * t  o 3 4 5 6 7 8 1 2
 * a  p 2 3 4 5 6 7 8 1
 *
 * Tiles without exclusions read off directly from the neighbor list interactingAtoms
 * and follows the same force accumulation method. If more there are more interactingTiles
 * than the size of the neighbor list initially allocated, the neighbor list is rebuilt
 * and the full tileset is computed. This should happen on the first step, and very rarely
 * afterwards.
 *
 * [out]forceBuffers    - forces on each atom to eventually be accumulated
 * [out]energyBuffer    - energyBuffer to eventually be accumulated
 * [in]posq             - x,y,z,charge
 * [in]exclusions       - 1024-bit flags denoting atom-atom exclusions for each tile
 * [in]exclusionTiles   - x,y denotes the indices of tiles that have an exclusion
 * [in]startTileIndex   - index into first tile to be processed
 * [in]numTileIndices   - number of tiles this context is responsible for processing
 * [in]int tiles        - the atom block for each tile
 * [in]interactionCount - total number of tiles that have an interaction
 * [in]maxTiles         - stores the size of the neighbourlist in case it needs
 *                      - to be expanded
 * [in]periodicBoxSize  - size of the Periodic Box, last dimension (w) not used
 * [in]invPeriodicBox   - inverse of the periodicBoxSize, pre-computed for speed
 * [in]blockCenter      - the center of each block in euclidean coordinates
 * [in]blockSize        - size of the each block, radiating from the center
 *                      - x is half the distance of total length
 *                      - y is half the distance of total width
 *                      - z is half the distance of total height
 *                      - w is not used
 * [in]interactingAtoms - a list of interactions within a given tile
 *
 */
KERNEL void computeNonbonded(
        GLOBAL mm_ulong* RESTRICT forceBuffers, GLOBAL mixed* RESTRICT energyBuffer, GLOBAL const real4* RESTRICT posq,
        GLOBAL const tileflags* RESTRICT exclusions, GLOBAL const int2* RESTRICT exclusionTiles, unsigned int startTileIndex, mm_ulong numTileIndices
#ifdef USE_CUTOFF
        , GLOBAL const int* RESTRICT tiles, GLOBAL const unsigned int* RESTRICT interactionCount, real4 periodicBoxSize, real4 invPeriodicBoxSize,
        real4 periodicBoxVecX, real4 periodicBoxVecY, real4 periodicBoxVecZ, unsigned int maxTiles, GLOBAL const real4* RESTRICT blockCenter,
        GLOBAL const real4* RESTRICT blockSize, GLOBAL const unsigned int* RESTRICT interactingAtoms, unsigned int maxSinglePairs,
        GLOBAL const int2* RESTRICT singlePairs
#endif
        PARAMETER_ARGUMENTS) {
    const unsigned int totalWarps = GLOBAL_SIZE/TILE_SIZE;
    const unsigned int warp = GLOBAL_ID/TILE_SIZE; // global warpIndex
    const unsigned int tgx = LOCAL_ID & (TILE_SIZE-1); // index within the warp
    const unsigned int tbx = LOCAL_ID - tgx;           // block warpIndex
    mixed energy = 0;
    INIT_DERIVATIVES
    LOCAL real4 localPosq[THREAD_BLOCK_SIZE];
    LOCAL real3 localForce[THREAD_BLOCK_SIZE];

    // First loop: process tiles that contain exclusions.

    const unsigned int firstExclusionTile = FIRST_EXCLUSION_TILE+warp*(LAST_EXCLUSION_TILE-FIRST_EXCLUSION_TILE)/totalWarps;
    const unsigned int lastExclusionTile = FIRST_EXCLUSION_TILE+(warp+1)*(LAST_EXCLUSION_TILE-FIRST_EXCLUSION_TILE)/totalWarps;
    for (int pos = firstExclusionTile; pos < lastExclusionTile; pos++) {
        const int2 tileIndices = exclusionTiles[pos];
        const unsigned int x = tileIndices.x;
        const unsigned int y = tileIndices.y;
        real3 force = make_real3(0);
        unsigned int atom1 = x*TILE_SIZE + tgx;
        real4 posq1 = posq[atom1];
        LOAD_ATOM1_PARAMETERS
#ifdef USE_EXCLUSIONS
        tileflags excl = exclusions[pos*TILE_SIZE+tgx];
#endif
        const bool hasExclusions = true;
        if (x == y) {
            // This tile is on the diagonal.
            localPosq[LOCAL_ID] = posq1;
            SYNC_WARPS;
            for (unsigned int j = 0; j < TILE_SIZE; j++) {
                int atom2 = tbx+j;
                real4 posq2 = localPosq[atom2];
                real3 delta = make_real3(posq2.x-posq1.x, posq2.y-posq1.y, posq2.z-posq1.z);
#ifdef USE_PERIODIC
                APPLY_PERIODIC_TO_DELTA(delta)
#endif
                real r2 = delta.x*delta.x + delta.y*delta.y + delta.z*delta.z;
                real invR = RSQRT(r2);
                real r = r2*invR;
                atom2 = y*TILE_SIZE+j;
                LOAD_ATOM2_PARAMETERS_FROM_GLOBAL
#ifdef USE_SYMMETRIC
                real dEdR = 0.0f;
#else
                real3 dEdR1 = make_real3(0);
                real3 dEdR2 = make_real3(0);
#endif
#ifdef USE_EXCLUSIONS
                bool isExcluded = (atom1 >= NUM_ATOMS || atom2 >= NUM_ATOMS || !(excl & 0x1));
#endif
                real tempEnergy = 0.0f;
                const real interactionScale = 0.5f;
                COMPUTE_INTERACTION
                energy += 0.5f*tempEnergy;
#ifdef INCLUDE_FORCES
#ifdef USE_SYMMETRIC
                force.x -= delta.x*dEdR;
                force.y -= delta.y*dEdR;
                force.z -= delta.z*dEdR;
#else
                force.x -= dEdR1.x;
                force.y -= dEdR1.y;
                force.z -= dEdR1.z;
#endif
#endif
#ifdef USE_EXCLUSIONS
                excl >>= 1;
#endif
                SYNC_WARPS;
            }
        }
        else {
            // This is an off-diagonal tile.
            unsigned int j = y*TILE_SIZE + tgx;
            localPosq[LOCAL_ID] = posq[j];
            localForce[LOCAL_ID] = 0.0f;
            SYNC_WARPS;
#ifdef USE_EXCLUSIONS
            excl = (excl >> tgx) | (excl << (TILE_SIZE - tgx));
#endif
            unsigned int tj = tgx;
            for (j = 0; j < TILE_SIZE; j++) {
                int atom2 = tbx+tj;
                real4 posq2 = localPosq[atom2];
                real3 delta = make_real3(posq2.x-posq1.x, posq2.y-posq1.y, posq2.z-posq1.z);
#ifdef USE_PERIODIC
                APPLY_PERIODIC_TO_DELTA(delta)
#endif
                real r2 = delta.x*delta.x + delta.y*delta.y + delta.z*delta.z;
                real invR = RSQRT(r2);
                real r = r2*invR;
                atom2 = y*TILE_SIZE+tj;
                LOAD_ATOM2_PARAMETERS_FROM_GLOBAL
#ifdef USE_SYMMETRIC
                real dEdR = 0.0f;
#else
                real3 dEdR1 = make_real3(0);
                real3 dEdR2 = make_real3(0);
#endif
#ifdef USE_EXCLUSIONS
                bool isExcluded = (atom1 >= NUM_ATOMS || atom2 >= NUM_ATOMS || !(excl & 0x1));
#endif
                real tempEnergy = 0.0f;
                const real interactionScale = 1.0f;
                COMPUTE_INTERACTION
                energy += tempEnergy;
#ifdef INCLUDE_FORCES
#ifdef USE_SYMMETRIC
                delta *= dEdR;
                force -= delta.xyz;
                localForce[tbx+tj] += delta.xyz;
#else // !USE_SYMMETRIC
                force -= dEdR1;
                localForce[tbx+tj] += dEdR2;
#endif // end USE_SYMMETRIC
#endif
#ifdef USE_EXCLUSIONS
                excl >>= 1;
#endif
                // cycles the indices
                // 0 1 2 3 4 5 6 7 -> 1 2 3 4 5 6 7 0
                tj = (tj + 1) & (TILE_SIZE - 1);
                SYNC_WARPS;
            }
            const unsigned int offset = y*TILE_SIZE + tgx;
            // write results for off diagonal tiles
#ifdef INCLUDE_FORCES
            ATOMIC_ADD(&forceBuffers[offset], static_cast<mm_ulong>(realToFixedPoint(localForce[LOCAL_ID].x)));
            ATOMIC_ADD(&forceBuffers[offset+PADDED_NUM_ATOMS], static_cast<mm_ulong>(realToFixedPoint(localForce[LOCAL_ID].y)));
            ATOMIC_ADD(&forceBuffers[offset+2*PADDED_NUM_ATOMS], static_cast<mm_ulong>(realToFixedPoint(localForce[LOCAL_ID].z)));
#endif
        }
        // Write results for on and off diagonal tiles
#ifdef INCLUDE_FORCES
        const unsigned int offset = x*TILE_SIZE + tgx;
        ATOMIC_ADD(&forceBuffers[offset], static_cast<mm_ulong>(realToFixedPoint(force.x)));
        ATOMIC_ADD(&forceBuffers[offset+PADDED_NUM_ATOMS], static_cast<mm_ulong>(realToFixedPoint(force.y)));
        ATOMIC_ADD(&forceBuffers[offset+2*PADDED_NUM_ATOMS], static_cast<mm_ulong>(realToFixedPoint(force.z)));
#endif
    }

    // Second loop: tiles without exclusions, either from the neighbor list (with cutoff) or just enumerating all
    // of them (no cutoff).

#ifdef USE_NEIGHBOR_LIST
    const unsigned int numTiles = interactionCount[0];
    if (numTiles > maxTiles)
        return; // There wasn't enough memory for the neighbor list.
    int pos = (int) (warp*(mm_long)numTiles/totalWarps);
    int end = (int) ((warp+1)*(mm_long)numTiles/totalWarps);
#else
    int pos = (int) (startTileIndex+warp*numTileIndices/totalWarps);
    int end = (int) (startTileIndex+(warp+1)*numTileIndices/totalWarps);
    int skipBase = 0;
    int currentSkipIndex = tbx;
    LOCAL volatile int skipTiles[THREAD_BLOCK_SIZE];
    skipTiles[LOCAL_ID] = -1;
#endif
    // atomIndices can probably be shuffled as well
    // but it probably wouldn't make things any faster
    LOCAL int atomIndices[THREAD_BLOCK_SIZE];

    while (pos < end) {
        const bool hasExclusions = false;
        real3 force = make_real3(0);
        bool includeTile = true;

        // Extract the coordinates of this tile.
        int x, y;
        bool singlePeriodicCopy = false;
#ifdef USE_NEIGHBOR_LIST
        x = tiles[pos];
        real4 blockSizeX = blockSize[x];
        singlePeriodicCopy = (0.5f*periodicBoxSize.x-blockSizeX.x >= MAX_CUTOFF &&
                              0.5f*periodicBoxSize.y-blockSizeX.y >= MAX_CUTOFF &&
                              0.5f*periodicBoxSize.z-blockSizeX.z >= MAX_CUTOFF);
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
                int2 tile = exclusionTiles[skipBase+tgx];
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
            LOAD_ATOM1_PARAMETERS
#ifdef USE_NEIGHBOR_LIST
            unsigned int j = interactingAtoms[pos*TILE_SIZE+tgx];
#else
            unsigned int j = y*TILE_SIZE + tgx;
#endif
            atomIndices[LOCAL_ID] = j;
            localForce[LOCAL_ID] = 0.0f;
            if (j < PADDED_NUM_ATOMS) {
                // Load position of atom j from global memory
                localPosq[LOCAL_ID] = posq[j];
            }
            else {
                localPosq[LOCAL_ID] = make_real4(0, 0, 0, 0);
            }
            SYNC_WARPS;
#ifdef USE_PERIODIC
            if (singlePeriodicCopy) {
                // The box is small enough that we can just translate all the atoms into a single periodic
                // box, then skip having to apply periodic boundary conditions later.
                real4 blockCenterX = blockCenter[x];
                APPLY_PERIODIC_TO_POS_WITH_CENTER(posq1, blockCenterX)
                APPLY_PERIODIC_TO_POS_WITH_CENTER(localPosq[LOCAL_ID], blockCenterX)
                SYNC_WARPS;
                unsigned int tj = tgx;
                for (j = 0; j < TILE_SIZE; j++) {
                    int atom2 = tbx+tj;
                    real4 posq2 = localPosq[atom2];
                    real3 delta = make_real3(posq2.x-posq1.x, posq2.y-posq1.y, posq2.z-posq1.z);
                    real r2 = delta.x*delta.x + delta.y*delta.y + delta.z*delta.z;
                    real invR = RSQRT(r2);
                    real r = r2*invR;
                    atom2 = atomIndices[tbx+tj];
                    LOAD_ATOM2_PARAMETERS_FROM_GLOBAL
#ifdef USE_SYMMETRIC
                    real dEdR = 0.0f;
#else
                    real3 dEdR1 = make_real3(0);
                    real3 dEdR2 = make_real3(0);
#endif
#ifdef USE_EXCLUSIONS
                    bool isExcluded = (atom1 >= NUM_ATOMS || atom2 >= NUM_ATOMS);
#endif
                    real tempEnergy = 0.0f;
                    const real interactionScale = 1.0f;
                    COMPUTE_INTERACTION
                    energy += tempEnergy;
#ifdef INCLUDE_FORCES
#ifdef USE_SYMMETRIC
                    delta *= dEdR;
                    force -= delta.xyz;
                    localForce[tbx+tj] += delta.xyz;
#else // !USE_SYMMETRIC
                    force -= dEdR1;
                    localForce[tbx+tj] += dEdR2;
#endif // end USE_SYMMETRIC
#endif
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
                    int atom2 = tbx+tj;
                    real4 posq2 = localPosq[atom2];
                    real3 delta = make_real3(posq2.x-posq1.x, posq2.y-posq1.y, posq2.z-posq1.z);
#ifdef USE_PERIODIC
                    APPLY_PERIODIC_TO_DELTA(delta)
#endif
                    real r2 = delta.x*delta.x + delta.y*delta.y + delta.z*delta.z;
                    real invR = RSQRT(r2);
                    real r = r2*invR;
                    atom2 = atomIndices[tbx+tj];
                    LOAD_ATOM2_PARAMETERS_FROM_GLOBAL
#ifdef USE_SYMMETRIC
                    real dEdR = 0.0f;
#else
                    real3 dEdR1 = make_real3(0);
                    real3 dEdR2 = make_real3(0);
#endif
#ifdef USE_EXCLUSIONS
                    bool isExcluded = (atom1 >= NUM_ATOMS || atom2 >= NUM_ATOMS);
#endif
                    real tempEnergy = 0.0f;
                    const real interactionScale = 1.0f;
                    COMPUTE_INTERACTION
                    energy += tempEnergy;
#ifdef INCLUDE_FORCES
#ifdef USE_SYMMETRIC
                    delta *= dEdR;
                    force -= delta.xyz;
                    localForce[tbx+tj] += delta.xyz;
#else // !USE_SYMMETRIC
                    force -= dEdR1;
                    localForce[tbx+tj] += dEdR2;
#endif // end USE_SYMMETRIC
#endif
                    tj = (tj + 1) & (TILE_SIZE - 1);
                    SYNC_WARPS;
                }
            }

            // Write results.
#ifdef INCLUDE_FORCES
            ATOMIC_ADD(&forceBuffers[atom1], static_cast<mm_ulong>(realToFixedPoint(force.x)));
            ATOMIC_ADD(&forceBuffers[atom1+PADDED_NUM_ATOMS], static_cast<mm_ulong>(realToFixedPoint(force.y)));
            ATOMIC_ADD(&forceBuffers[atom1+2*PADDED_NUM_ATOMS], static_cast<mm_ulong>(realToFixedPoint(force.z)));
#ifdef USE_NEIGHBOR_LIST
            unsigned int atom2 = atomIndices[LOCAL_ID];
#else
            unsigned int atom2 = y*TILE_SIZE + tgx;
#endif
            if (atom2 < PADDED_NUM_ATOMS) {
                ATOMIC_ADD(&forceBuffers[atom2], static_cast<mm_ulong>(realToFixedPoint(localForce[LOCAL_ID].x)));
                ATOMIC_ADD(&forceBuffers[atom2+PADDED_NUM_ATOMS], static_cast<mm_ulong>(realToFixedPoint(localForce[LOCAL_ID].y)));
                ATOMIC_ADD(&forceBuffers[atom2+2*PADDED_NUM_ATOMS], static_cast<mm_ulong>(realToFixedPoint(localForce[LOCAL_ID].z)));
            }
#endif
        }
        pos++;
    }

    // Third loop: single pairs that aren't part of a tile.

#if USE_NEIGHBOR_LIST
    const unsigned int numPairs = interactionCount[1];
    if (numPairs > maxSinglePairs)
        return; // There wasn't enough memory for the neighbor list.
    for (int i = GLOBAL_ID; i < numPairs; i += GLOBAL_SIZE) {
        int2 pair = singlePairs[i];
        int atom1 = pair.x;
        int atom2 = pair.y;
        real4 posq1 = posq[atom1];
        real4 posq2 = posq[atom2];
        LOAD_ATOM1_PARAMETERS
        LOAD_ATOM2_PARAMETERS_FROM_GLOBAL
        real3 delta = make_real3(posq2.x-posq1.x, posq2.y-posq1.y, posq2.z-posq1.z);
#ifdef USE_PERIODIC
        APPLY_PERIODIC_TO_DELTA(delta)
#endif
        real r2 = delta.x*delta.x + delta.y*delta.y + delta.z*delta.z;
        real invR = RSQRT(r2);
        real r = r2*invR;
#ifdef USE_SYMMETRIC
        real dEdR = 0.0f;
#else
        real3 dEdR1 = make_real3(0);
        real3 dEdR2 = make_real3(0);
#endif
        bool hasExclusions = false;
        bool isExcluded = false;
        real tempEnergy = 0.0f;
        const real interactionScale = 1.0f;
        COMPUTE_INTERACTION
        energy += tempEnergy;
#ifdef INCLUDE_FORCES
#ifdef USE_SYMMETRIC
        real3 dEdR1 = delta*dEdR;
        real3 dEdR2 = -dEdR1;
#endif
        saveSingleForce(atom1, -dEdR1, forceBuffers);
        saveSingleForce(atom2, -dEdR2, forceBuffers);
#endif
    }
#endif
#ifdef INCLUDE_ENERGY
    energyBuffer[GLOBAL_ID] += energy;
#endif
    SAVE_DERIVATIVES
}
