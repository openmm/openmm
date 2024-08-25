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
 * Tiles without exclusions read off directly from the neighbourlist interactingAtoms
 * and follows the same force accumulation method. If more there are more interactingTiles
 * than the size of the neighbourlist initially allocated, the neighbourlist is rebuilt
 * and the full tileset is computed. This should happen on the first step, and very rarely
 * afterwards.
 *
 * On diagonal exclusion tiles use __shfl to broadcast. For all other types of tiles __shfl
 * is used to pass around the forces, positions, and parameters when computing the forces.
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
extern "C" __launch_bounds__(THREAD_BLOCK_SIZE) __global__ void computeNonbonded(
        unsigned long long* __restrict__ forceBuffers, mixed* __restrict__ energyBuffer, const real4* __restrict__ posq, const tileflags* __restrict__ exclusions,
        const int2* __restrict__ exclusionTiles, unsigned int startTileIndex, unsigned long long numTileIndices
#ifdef USE_CUTOFF
        , const int* __restrict__ tiles, const unsigned int* __restrict__ interactionCount, real4 periodicBoxSize, real4 invPeriodicBoxSize,
        real4 periodicBoxVecX, real4 periodicBoxVecY, real4 periodicBoxVecZ, unsigned int maxTiles, const real4* __restrict__ blockCenter,
        const real4* __restrict__ blockSize, const unsigned int* __restrict__ interactingAtoms, unsigned int maxSinglePairs,
        const int2* __restrict__ singlePairs
#endif
        PARAMETER_ARGUMENTS) {
    const unsigned int totalWarps = (blockDim.x*gridDim.x)/TILE_SIZE;
    const unsigned int warp = blockIdx.x*(THREAD_BLOCK_SIZE/TILE_SIZE) + threadIdx.x/TILE_SIZE; // global warpIndex
    const unsigned int tgx = threadIdx.x & (TILE_SIZE-1); // index within the warp
    const unsigned int tbx = threadIdx.x - tgx;           // block warpIndex
    mixed energy = 0;
    INIT_DERIVATIVES

    // First loop: process tiles that contain exclusions.

    for (int pos = FIRST_EXCLUSION_TILE+warp; pos < LAST_EXCLUSION_TILE; pos+=totalWarps) {
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
        if (x == y) {
            // This tile is on the diagonal.
            real4 shflPosq = posq1;

            // we do not need to fetch parameters from global since this is a symmetric tile
            // instead we can broadcast the values using shuffle
            #pragma unroll 2
            for (unsigned int j = 0; j < TILE_SIZE; j++) {
                real4 posq2;
                BROADCAST_WARP_DATA
                real3 delta = make_real3(posq2.x-posq1.x, posq2.y-posq1.y, posq2.z-posq1.z);
#ifdef USE_PERIODIC
                APPLY_PERIODIC_TO_DELTA(delta)
#endif
                real r2 = delta.x*delta.x + delta.y*delta.y + delta.z*delta.z;
#ifdef USE_CUTOFF
                if (r2 < MAX_CUTOFF_SQUARED) {
#endif
                    real invR = RSQRT(r2);
                    real r = r2*invR;
                    LOAD_ATOM2_PARAMETERS
                    int atom2 = y*TILE_SIZE+j;
#ifdef USE_SYMMETRIC
                    real dEdR = 0.0f;
#else
                    real3 dEdR1 = make_real3(0);
                    real3 dEdR2 = make_real3(0);
#endif
#ifdef USE_EXCLUSIONS
                    bool isExcluded = !(atom1 < NUM_ATOMS && atom2 < NUM_ATOMS && (excl & 1));
#endif
                    real tempEnergy = 0.0f;
                    const real interactionScale = 0.5f;
                    COMPUTE_INTERACTION
#ifdef INCLUDE_ENERGY
                    energy += 0.5f*tempEnergy;
#endif
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
#ifdef USE_CUTOFF
                }
#endif
#ifdef USE_EXCLUSIONS
                excl >>= 1;
#endif
            }
        }
        else {
            // This is an off-diagonal tile.
            unsigned int j = y*TILE_SIZE + tgx;
            int atomIndex = j;
            real4 shflPosq = posq[j];
            real3 shflForce;
            shflForce.x = 0.0f;
            shflForce.y = 0.0f;
            shflForce.z = 0.0f;
            DECLARE_LOCAL_PARAMETERS
            LOAD_LOCAL_PARAMETERS_FROM_GLOBAL
#ifdef USE_EXCLUSIONS
            excl = (excl >> tgx) | (excl << (TILE_SIZE - tgx));
#endif
            #pragma unroll 2
            for (j = 0; j < TILE_SIZE; j++) {
                real4 posq2 = shflPosq;
                real3 delta = make_real3(posq2.x-posq1.x, posq2.y-posq1.y, posq2.z-posq1.z);
#ifdef USE_PERIODIC
                APPLY_PERIODIC_TO_DELTA(delta)
#endif
                real r2 = delta.x*delta.x + delta.y*delta.y + delta.z*delta.z;
#ifdef USE_CUTOFF
                if (r2 < MAX_CUTOFF_SQUARED) {
#endif
                    real invR = RSQRT(r2);
                    real r = r2*invR;
                    LOAD_ATOM2_PARAMETERS
                    int atom2 = atomIndex;
#ifdef USE_SYMMETRIC
                    real dEdR = 0.0f;
#else
                    real3 dEdR1 = make_real3(0);
                    real3 dEdR2 = make_real3(0);
#endif
#ifdef USE_EXCLUSIONS
                    bool isExcluded = !(atom1 < NUM_ATOMS && atom2 < NUM_ATOMS && (excl & 1));
#endif
                    real tempEnergy = 0.0f;
                    const real interactionScale = 1.0f;
                    COMPUTE_INTERACTION
#ifdef INCLUDE_ENERGY
                    energy += tempEnergy;
#endif
#ifdef INCLUDE_FORCES
#ifdef USE_SYMMETRIC
                    delta *= dEdR;
                    force.x -= delta.x;
                    force.y -= delta.y;
                    force.z -= delta.z;
                    shflForce.x += delta.x;
                    shflForce.y += delta.y;
                    shflForce.z += delta.z;
#else // !USE_SYMMETRIC
                    force.x -= dEdR1.x;
                    force.y -= dEdR1.y;
                    force.z -= dEdR1.z;
                    shflForce.x += dEdR2.x;
                    shflForce.y += dEdR2.y;
                    shflForce.z += dEdR2.z;
#endif // end USE_SYMMETRIC
#endif
#ifdef USE_CUTOFF
                }
#endif
#ifdef USE_EXCLUSIONS
                excl >>= 1;
#endif
                SHUFFLE_WARP_DATA
                atomIndex = warpRotateLeft<TILE_SIZE>(atomIndex);
            }
            const unsigned int offset = atomIndex;
            // write results for off diagonal tiles
#ifdef INCLUDE_FORCES
            atomicAdd(&forceBuffers[offset], static_cast<unsigned long long>(realToFixedPoint(shflForce.x)));
            atomicAdd(&forceBuffers[offset+PADDED_NUM_ATOMS], static_cast<unsigned long long>(realToFixedPoint(shflForce.y)));
            atomicAdd(&forceBuffers[offset+2*PADDED_NUM_ATOMS], static_cast<unsigned long long>(realToFixedPoint(shflForce.z)));
#endif
        }
        // Write results for on and off diagonal tiles
#ifdef INCLUDE_FORCES
        const unsigned int offset = x*TILE_SIZE + tgx;
        atomicAdd(&forceBuffers[offset], static_cast<unsigned long long>(realToFixedPoint(force.x)));
        atomicAdd(&forceBuffers[offset+PADDED_NUM_ATOMS], static_cast<unsigned long long>(realToFixedPoint(force.y)));
        atomicAdd(&forceBuffers[offset+2*PADDED_NUM_ATOMS], static_cast<unsigned long long>(realToFixedPoint(force.z)));
#endif
    }

    // Second loop: tiles without exclusions, either from the neighbor list (with cutoff) or just enumerating all
    // of them (no cutoff).

#ifdef USE_NEIGHBOR_LIST
    const unsigned int numTiles = interactionCount[0];
    if (numTiles > maxTiles)
        return; // There wasn't enough memory for the neighbor list.
    for (unsigned int pos0 = warp; pos0 < LAST_EXCLUSION_TILE+numTiles; pos0+=totalWarps) {
        // Skip warps that may be still busy in the first loop
        if (pos0 < LAST_EXCLUSION_TILE) {
            continue;
        }
        const unsigned int pos = pos0-LAST_EXCLUSION_TILE;
#else
    int pos = (int) (startTileIndex+warp*numTileIndices/totalWarps);
    int end = (int) (startTileIndex+(warp+1)*numTileIndices/totalWarps);
    int skipBase = 0;
    int currentSkipIndex = tbx;
    int skipTiles;
    skipTiles = -1;
    for (; pos < end; pos++) {
#endif
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

        while (SHFL(skipTiles, TILE_SIZE-1) < pos) {
            if (skipBase+tgx < NUM_TILES_WITH_EXCLUSIONS) {
                int2 tile = exclusionTiles[skipBase+tgx];
                skipTiles = tile.x + tile.y*NUM_BLOCKS - tile.y*(tile.y+1)/2;
            }
            else
                skipTiles = end;
            skipBase += TILE_SIZE;
            currentSkipIndex = 0;
        }
        while (SHFL(skipTiles, currentSkipIndex) < pos)
            currentSkipIndex++;
        includeTile = (SHFL(skipTiles, currentSkipIndex) != pos);
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
            int atomIndex = j;
            DECLARE_LOCAL_PARAMETERS
            real4 shflPosq;
            real3 shflForce;
            shflForce.x = 0.0f;
            shflForce.y = 0.0f;
            shflForce.z = 0.0f;
            if (j < PADDED_NUM_ATOMS) {
                // Load position of atom j from from global memory
                shflPosq = posq[j];
                LOAD_LOCAL_PARAMETERS_FROM_GLOBAL
            }
            else {
                shflPosq = make_real4(0, 0, 0, 0);
                CLEAR_LOCAL_PARAMETERS
            }
#ifdef USE_PERIODIC
            if (singlePeriodicCopy) {
                // The box is small enough that we can just translate all the atoms into a single periodic
                // box, then skip having to apply periodic boundary conditions later.
                real4 blockCenterX = blockCenter[x];
                APPLY_PERIODIC_TO_POS_WITH_CENTER(posq1, blockCenterX)
                APPLY_PERIODIC_TO_POS_WITH_CENTER(shflPosq, blockCenterX)
                #pragma unroll 2
                for (j = 0; j < TILE_SIZE; j++) {
                    real4 posq2 = shflPosq;
                    real3 delta = make_real3(posq2.x-posq1.x, posq2.y-posq1.y, posq2.z-posq1.z);
                    real r2 = delta.x*delta.x + delta.y*delta.y + delta.z*delta.z;
#ifdef USE_CUTOFF
                    if (r2 < MAX_CUTOFF_SQUARED) {
#endif
                        real invR = RSQRT(r2);
                        real r = r2*invR;
                        LOAD_ATOM2_PARAMETERS
                        int atom2 = atomIndex;
#ifdef USE_SYMMETRIC
                        real dEdR = 0.0f;
#else
                        real3 dEdR1 = make_real3(0);
                        real3 dEdR2 = make_real3(0);
#endif
#ifdef USE_EXCLUSIONS
                        bool isExcluded = !(atom1 < NUM_ATOMS && atom2 < NUM_ATOMS);
#endif
                        real tempEnergy = 0.0f;
                        const real interactionScale = 1.0f;
                        COMPUTE_INTERACTION
#ifdef INCLUDE_ENERGY
                        energy += tempEnergy;
#endif
#ifdef INCLUDE_FORCES
#ifdef USE_SYMMETRIC
                        delta *= dEdR;
                        force.x -= delta.x;
                        force.y -= delta.y;
                        force.z -= delta.z;
                        shflForce.x += delta.x;
                        shflForce.y += delta.y;
                        shflForce.z += delta.z;
#else // !USE_SYMMETRIC
                        force.x -= dEdR1.x;
                        force.y -= dEdR1.y;
                        force.z -= dEdR1.z;
                        shflForce.x += dEdR2.x;
                        shflForce.y += dEdR2.y;
                        shflForce.z += dEdR2.z;
#endif // end USE_SYMMETRIC
#endif
#ifdef USE_CUTOFF
                    }
#endif
                    SHUFFLE_WARP_DATA
                    atomIndex = warpRotateLeft<TILE_SIZE>(atomIndex);
                }
            }
            else
#endif
            {
                // We need to apply periodic boundary conditions separately for each interaction.
                #pragma unroll 2
                for (j = 0; j < TILE_SIZE; j++) {
                    real4 posq2 = shflPosq;
                    real3 delta = make_real3(posq2.x-posq1.x, posq2.y-posq1.y, posq2.z-posq1.z);
#ifdef USE_PERIODIC
                    APPLY_PERIODIC_TO_DELTA(delta)
#endif
                    real r2 = delta.x*delta.x + delta.y*delta.y + delta.z*delta.z;
#ifdef USE_CUTOFF
                    if (r2 < MAX_CUTOFF_SQUARED) {
#endif
                        real invR = RSQRT(r2);
                        real r = r2*invR;
                        LOAD_ATOM2_PARAMETERS
                        int atom2 = atomIndex;
#ifdef USE_SYMMETRIC
                        real dEdR = 0.0f;
#else
                        real3 dEdR1 = make_real3(0);
                        real3 dEdR2 = make_real3(0);
#endif
#ifdef USE_EXCLUSIONS
                        bool isExcluded = !(atom1 < NUM_ATOMS && atom2 < NUM_ATOMS);
#endif
                        real tempEnergy = 0.0f;
                        const real interactionScale = 1.0f;
                        COMPUTE_INTERACTION
#ifdef INCLUDE_ENERGY
                        energy += tempEnergy;
#endif
#ifdef INCLUDE_FORCES
#ifdef USE_SYMMETRIC
                        delta *= dEdR;
                        force.x -= delta.x;
                        force.y -= delta.y;
                        force.z -= delta.z;
                        shflForce.x += delta.x;
                        shflForce.y += delta.y;
                        shflForce.z += delta.z;
#else // !USE_SYMMETRIC
                        force.x -= dEdR1.x;
                        force.y -= dEdR1.y;
                        force.z -= dEdR1.z;
                        shflForce.x += dEdR2.x;
                        shflForce.y += dEdR2.y;
                        shflForce.z += dEdR2.z;
#endif // end USE_SYMMETRIC
#endif
#ifdef USE_CUTOFF
                    }
#endif
                    SHUFFLE_WARP_DATA
                    atomIndex = warpRotateLeft<TILE_SIZE>(atomIndex);
                }
            }

            // Write results.
#ifdef INCLUDE_FORCES
            atomicAdd(&forceBuffers[atom1], static_cast<unsigned long long>(realToFixedPoint(force.x)));
            atomicAdd(&forceBuffers[atom1+PADDED_NUM_ATOMS], static_cast<unsigned long long>(realToFixedPoint(force.y)));
            atomicAdd(&forceBuffers[atom1+2*PADDED_NUM_ATOMS], static_cast<unsigned long long>(realToFixedPoint(force.z)));
            unsigned int atom2 = atomIndex;
            if (atom2 < PADDED_NUM_ATOMS) {
                atomicAdd(&forceBuffers[atom2], static_cast<unsigned long long>(realToFixedPoint(shflForce.x)));
                atomicAdd(&forceBuffers[atom2+PADDED_NUM_ATOMS], static_cast<unsigned long long>(realToFixedPoint(shflForce.y)));
                atomicAdd(&forceBuffers[atom2+2*PADDED_NUM_ATOMS], static_cast<unsigned long long>(realToFixedPoint(shflForce.z)));
            }
#endif
        }
    }

    // Third loop: single pairs that aren't part of a tile.

#if USE_NEIGHBOR_LIST
    const unsigned int numPairs = interactionCount[1];
    if (numPairs > maxSinglePairs)
        return; // There wasn't enough memory for the neighbor list.
    for (int i = blockIdx.x*blockDim.x+threadIdx.x; i < numPairs; i += blockDim.x*gridDim.x) {
        int2 pair = singlePairs[i];
        int atom1 = pair.x;
        int atom2 = pair.y;
        real4 posq1 = posq[atom1];
        real4 posq2 = posq[atom2];
        real3 delta = make_real3(posq2.x-posq1.x, posq2.y-posq1.y, posq2.z-posq1.z);
#ifdef USE_PERIODIC
        APPLY_PERIODIC_TO_DELTA(delta)
#endif
        real r2 = delta.x*delta.x + delta.y*delta.y + delta.z*delta.z;
        if (r2 < MAX_CUTOFF_SQUARED) {
            LOAD_ATOM1_PARAMETERS
            int j = atom2;
            DECLARE_LOCAL_PARAMETERS
            LOAD_LOCAL_PARAMETERS_FROM_GLOBAL
            LOAD_ATOM2_PARAMETERS
            real invR = RSQRT(r2);
            real r = r2*invR;
#ifdef USE_SYMMETRIC
            real dEdR = 0.0f;
#else
            real3 dEdR1 = make_real3(0);
            real3 dEdR2 = make_real3(0);
#endif
            bool isExcluded = false;
            real tempEnergy = 0.0f;
            const real interactionScale = 1.0f;
            COMPUTE_INTERACTION
#ifdef INCLUDE_ENERGY
            energy += tempEnergy;
#endif
#ifdef INCLUDE_FORCES
#ifdef USE_SYMMETRIC
            real3 dEdR1 = delta*dEdR;
            real3 dEdR2 = -dEdR1;
#endif
            atomicAdd(&forceBuffers[atom1], static_cast<unsigned long long>(realToFixedPoint(-dEdR1.x)));
            atomicAdd(&forceBuffers[atom1+PADDED_NUM_ATOMS], static_cast<unsigned long long>(realToFixedPoint(-dEdR1.y)));
            atomicAdd(&forceBuffers[atom1+2*PADDED_NUM_ATOMS], static_cast<unsigned long long>(realToFixedPoint(-dEdR1.z)));
            atomicAdd(&forceBuffers[atom2], static_cast<unsigned long long>(realToFixedPoint(-dEdR2.x)));
            atomicAdd(&forceBuffers[atom2+PADDED_NUM_ATOMS], static_cast<unsigned long long>(realToFixedPoint(-dEdR2.y)));
            atomicAdd(&forceBuffers[atom2+2*PADDED_NUM_ATOMS], static_cast<unsigned long long>(realToFixedPoint(-dEdR2.z)));
#endif
        }
    }
#endif
#ifdef INCLUDE_ENERGY
    energyBuffer[blockIdx.x*blockDim.x+threadIdx.x] += energy;
#endif
    SAVE_DERIVATIVES
}
