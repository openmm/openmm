#define WARPS_PER_GROUP (THREAD_BLOCK_SIZE/TILE_SIZE)

//support for 64 bit shuffles
static __inline__ __device__ float real_shfl(float var, int srcLane) {
    return SHFL(var, srcLane);
}

static __inline__ __device__ float real_shfl(int var, int srcLane) {
    return SHFL(var, srcLane);
}

static __inline__ __device__ double real_shfl(double var, int srcLane) {
    int hi, lo;
    asm volatile("mov.b64 { %0, %1 }, %2;" : "=r"(lo), "=r"(hi) : "d"(var));
    hi = SHFL(hi, srcLane);
    lo = SHFL(lo, srcLane);
    return __hiloint2double( hi, lo );
}

static __inline__ __device__ long long real_shfl(long long var, int srcLane) {
    int hi, lo;
    asm volatile("mov.b64 { %0, %1 }, %2;" : "=r"(lo), "=r"(hi) : "l"(var));
    hi = SHFL(hi, srcLane);
    lo = SHFL(lo, srcLane);
    // unforunately there isn't an __nv_hiloint2long(hi,lo) intrinsic cast
    int2 fuse; fuse.x = lo; fuse.y = hi;
    return *reinterpret_cast<long long*>(&fuse);
}

/**
 * Shuffle a 3-component vector across lanes in a warp.
 */
inline __device__ real3 real_shfl3(const real3& var, int srcLane) {
    real3 result;
    result.x = real_shfl(var.x, srcLane);
    result.y = real_shfl(var.y, srcLane);
    result.z = real_shfl(var.z, srcLane);
    return result;
}

/**
 * Shuffle a 4-component vector across lanes in a warp.
 */
inline __device__ real4 real_shfl4(const real4& var, int srcLane) {
    real4 result;
    result.x = real_shfl(var.x, srcLane);
    result.y = real_shfl(var.y, srcLane);
    result.z = real_shfl(var.z, srcLane);
    result.w = real_shfl(var.w, srcLane);
    return result;
}

/**
 * Save the force on a single atom.
 */
__device__ void saveSingleForce(int atom, real3 force, unsigned long long* forceBuffers) {
    if (force.x != 0)
        atomicAdd(&forceBuffers[atom], static_cast<unsigned long long>(realToFixedPoint(force.x)));
    if (force.y != 0)
        atomicAdd(&forceBuffers[atom+PADDED_NUM_ATOMS], static_cast<unsigned long long>(realToFixedPoint(force.y)));
    if (force.z != 0)
        atomicAdd(&forceBuffers[atom+2*PADDED_NUM_ATOMS], static_cast<unsigned long long>(realToFixedPoint(force.z)));
}

/**
 * Rotate a lane's 32-bit diagonal-tile exclusion mask right by tgx + 1 bits so bit 0
 * corresponds to initial partner lane (tgx + 1) & (TILE_SIZE - 1).
 */
inline __device__ tileflags alignDiagonalExclusions(tileflags excl, unsigned int tgx) {
    return __funnelshift_r(excl, excl, tgx + 1);
}

/**
 * Return whether the current diagonal-tile shift evaluates a non-duplicate pair within the cutoff.
 * On the final shift (j == TILE_SIZE/2 - 1), lane tgx is paired with (tgx + TILE_SIZE/2) & (TILE_SIZE - 1),
 * so only the lower half-warp evaluates the pair.
 */
inline __device__ bool isActiveDiagonalPair(unsigned int j, unsigned int tgx, real r2) {
    bool includePair = !(j == TILE_SIZE/2 - 1 && tgx >= TILE_SIZE/2);
#ifdef USE_CUTOFF
    includePair &= (r2 < MAX_CUTOFF_SQUARED);
#endif
    return includePair;
}

/**
 * Accumulate equal and opposite pair forces delta * dEdR onto atom1 (force) and atom2 (shflForce).
 */
inline __device__ void accumulateSymmetricPairForce(real3 delta, real dEdR, real3& force, real3& shflForce) {
    delta *= dEdR;
    force.x -= delta.x;
    force.y -= delta.y;
    force.z -= delta.z;
    shflForce.x += delta.x;
    shflForce.y += delta.y;
    shflForce.z += delta.z;
}

/**
 * Fold the partner reaction forces back onto force after the initial shift plus TILE_SIZE/2 loop
 * shifts (1 + TILE_SIZE/2 shifts total, which leaves lane tgx's accumulator in lane tgx + TILE_SIZE/2 - 1).
 */
inline __device__ real3 addDiagonalPartnerForce(const real3& force, const real3& shflForce, unsigned int tgx) {
    return force + real_shfl3(shflForce, tgx + TILE_SIZE/2 - 1);
}

/**
 * Accumulate symmetric pair forces onto atom1 (force) and atom2 (shflForce).
 */
inline __device__ void accumulatePairForces(real3 delta, real dEdR, real3& force, real3& shflForce) {
    accumulateSymmetricPairForce(delta, dEdR, force, shflForce);
}

/**
 * Accumulate asymmetric pair derivatives onto atom1 (force) and atom2 (shflForce).
 */
inline __device__ void accumulatePairForces(const real3& dEdR1, const real3& dEdR2, real3& force, real3& shflForce) {
    force -= dEdR1;
    shflForce += dEdR2;
}

#ifdef USE_NEIGHBOR_LIST
/**
 * Schedules neighbor-list tiles across warps, assigning 2 initial tiles statically to each warp
 * before claiming dynamic 4-tile batches from workUnitCount.
 */
struct TileWorkQueue {
    static const unsigned int INITIAL_TILES_PER_WARP = 2;
    static const unsigned int TILES_PER_WORK_UNIT = 4;
    int numTiles;
    int dynamicStart;
    int start;
    int end;
    bool firstBatch;

    inline __device__ TileWorkQueue(unsigned int warp, unsigned int totalWarps, unsigned int totalTiles) :
            numTiles((int) totalTiles),
            dynamicStart((int) (totalWarps * INITIAL_TILES_PER_WARP)),
            start((int) (warp * INITIAL_TILES_PER_WARP)),
            end(min((int) (warp * INITIAL_TILES_PER_WARP + INITIAL_TILES_PER_WARP), (int) totalTiles)),
            firstBatch(true) {
    }

    inline __device__ bool claimNextBatch(unsigned int tgx, unsigned int* workUnitCount) {
        if (firstBatch) {
            firstBatch = false;
            return (start < numTiles);
        }
        if (dynamicStart >= numTiles)
            return false;
        if (tgx == 0)
            start = dynamicStart + (int) atomicAdd(workUnitCount, TILES_PER_WORK_UNIT);
        start = SHFL(start, 0);
        end = min(start + (int) TILES_PER_WORK_UNIT, numTiles);
        return (start < numTiles);
    }
};

/**
 * Return whether block x fits within a single periodic image with margin MAX_CUTOFF so minimum-image
 * wrapping can be applied once per atom rather than once per pair.
 */
inline __device__ bool canUseSinglePeriodicCopy(const real4& blockSizeX, const real4& periodicBoxSize) {
#ifdef USE_PERIODIC
    return (0.5f*periodicBoxSize.x-blockSizeX.x >= MAX_CUTOFF &&
            0.5f*periodicBoxSize.y-blockSizeX.y >= MAX_CUTOFF &&
            0.5f*periodicBoxSize.z-blockSizeX.z >= MAX_CUTOFF);
#else
    return true;
#endif
}

/**
 * Return whether pos lies within MAX_CUTOFF (with a 1.0001 tolerance factor on squared distance)
 * of block x's axis-aligned bounding box.
 */
inline __device__ bool isAtomNearBlock(const real4& pos, const real4& blockCenterX, const real4& blockSizeX) {
    real3 blockDelta = make_real3(
        max((real) 0, FABS(pos.x-blockCenterX.x)-blockSizeX.x),
        max((real) 0, FABS(pos.y-blockCenterX.y)-blockSizeX.y),
        max((real) 0, FABS(pos.z-blockCenterX.z)-blockSizeX.z));
    return (blockDelta.x*blockDelta.x + blockDelta.y*blockDelta.y + blockDelta.z*blockDelta.z <= (real) 1.0001f*MAX_CUTOFF_SQUARED);
}

/**
 * Round the number of in-range neighbor atoms up to the smallest supported power-of-two sub-tile size
 * (4, 8, 16, or TILE_SIZE).
 */
inline __device__ unsigned int selectSubTileSize(unsigned int numAtomsInCutoff) {
    return (numAtomsInCutoff <= 4u ? 4u : (numAtomsInCutoff <= 8u ? 8u : (numAtomsInCutoff <= 16u ? 16u : (unsigned int) TILE_SIZE)));
}

/**
 * Map lane tgx to the lane holding the (tgx % subTileSize)-th set bit in atomsInCutoffMask,
 * replicating the compacted sub-tile across all TILE_SIZE / subTileSize segments of the warp.
 */
inline __device__ unsigned int findSubTileSourceLane(unsigned int tgx, unsigned int subTileSize,
        unsigned int atomsInCutoffMask, unsigned int numAtomsInCutoff) {
    const unsigned int indexInSubTile = tgx & (subTileSize - 1u);
    unsigned int sourceLane;
#if __CUDA_ARCH__ >= 700
    asm("fns.b32 %0, %1, 0, %2;" : "=r"(sourceLane) : "r"(atomsInCutoffMask), "r"(indexInSubTile + 1u));
#else
    if (indexInSubTile < numAtomsInCutoff) {
        unsigned int mask = atomsInCutoffMask;
        for (unsigned int k = 0; k < indexInSubTile; k++)
            mask &= mask - 1u;
        sourceLane = __ffs(mask) - 1;
    }
    else
        sourceLane = 0u;
#endif
    return sourceLane & (TILE_SIZE - 1u);
}

/**
 * Sum partial atom2 forces across replicated sub-tiles of width subTileSize using butterfly shuffles.
 */
inline __device__ real3 reduceSubTileForce(real3 shflForce, unsigned int tgx, unsigned int subTileSize) {
    #pragma unroll
    for (unsigned int stride = 4; stride < TILE_SIZE; stride *= 2) {
        if (subTileSize <= stride)
            shflForce += real_shfl3(shflForce, tgx ^ stride);
    }
    return shflForce;
}
#endif

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
 * On CUDA devices that support the shuffle intrinsic, on diagonal exclusion tiles use
 * __shfl to broadcast. For all other types of tiles __shfl is used to pass around the 
 * forces, positions, and parameters when computing the forces. 
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
extern "C" __global__ COMPUTE_NONBONDED_LAUNCH_BOUNDS void computeNonbonded(
        unsigned long long* __restrict__ forceBuffers, mixed* __restrict__ energyBuffer, const real4* __restrict__ posq, const tileflags* __restrict__ exclusions,
        const int2* __restrict__ exclusionTiles, unsigned int startTileIndex, unsigned long long numTileIndices
#ifdef USE_CUTOFF
        , const int* __restrict__ tiles, const unsigned int* __restrict__ interactionCount, real4 periodicBoxSize, real4 invPeriodicBoxSize, 
        real4 periodicBoxVecX, real4 periodicBoxVecY, real4 periodicBoxVecZ, unsigned int maxTiles, const real4* __restrict__ blockCenter,
        const real4* __restrict__ blockSize, const unsigned int* __restrict__ interactingAtoms, unsigned int maxSinglePairs,
        const int2* __restrict__ singlePairs, unsigned int* __restrict__ workUnitCount
#endif
        PARAMETER_ARGUMENTS) {
    const unsigned int totalWarps = (blockDim.x*gridDim.x)/TILE_SIZE;
    const unsigned int warp = (blockIdx.x*blockDim.x+threadIdx.x)/TILE_SIZE; // global warpIndex
    const unsigned int tgx = threadIdx.x & (TILE_SIZE-1); // index within the warp
    const unsigned int tbx = threadIdx.x - tgx;           // block warpIndex
    mixed energy = 0;
    INIT_DERIVATIVES

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
#if defined(USE_SYMMETRIC) && defined(USE_EXCLUSIONS)
            real4 shflPosq = posq1;
            real3 shflForce = make_real3(0);
            DECLARE_LOCAL_PARAMETERS
            COPY_ATOM1_TO_LOCAL_PARAMETERS
            excl = alignDiagonalExclusions(excl, tgx);
            SHUFFLE_WARP_DATA

            // Evaluate all TILE_SIZE*(TILE_SIZE - 1)/2 = 496 unordered pairs across TILE_SIZE/2 = 16 shifts:
            // shifts 0..14 evaluate 15*32 = 480 distinct pairs across all 32 lanes, and shift 15 pairs lane tgx
            // with (tgx + 16) & 31 so restricting evaluation to tgx < 16 covers the remaining 16 pairs once.

            #pragma unroll 4
            for (unsigned int j = 0; j < TILE_SIZE/2; j++) {
                int atom2 = y*TILE_SIZE + ((tgx+j+1) & (TILE_SIZE-1));
                real4 posq2 = shflPosq;
                real3 delta = make_real3(posq2.x-posq1.x, posq2.y-posq1.y, posq2.z-posq1.z);
#ifdef USE_PERIODIC
                APPLY_PERIODIC_TO_DELTA(delta)
#endif
                real r2 = delta.x*delta.x + delta.y*delta.y + delta.z*delta.z;
                if (atom1 < NUM_ATOMS && atom2 < NUM_ATOMS && isActiveDiagonalPair(j, tgx, r2)) {
                    real invR = RSQRT(r2);
                    real r = r2*invR;
                    LOAD_ATOM2_PARAMETERS
                    real dEdR = 0.0f;
                    bool isExcluded = !(excl & 0x1);
                    real tempEnergy = 0.0f;
                    const real interactionScale = 1.0f;
                    COMPUTE_INTERACTION
                    energy += tempEnergy;
#ifdef INCLUDE_FORCES
                    accumulateSymmetricPairForce(delta, dEdR, force, shflForce);
#endif
                }
                SHUFFLE_WARP_DATA
                excl >>= 1;
            }
#ifdef INCLUDE_FORCES
            force = addDiagonalPartnerForce(force, shflForce, tgx);
#endif
#else
            real4 shflPosq = posq1;

            // we do not need to fetch parameters from global since this is a symmetric tile
            // instead we can broadcast the values using shuffle
            for (unsigned int j = 0; j < TILE_SIZE; j++) {
                int atom2 = tbx+j;
                real4 posq2;
                BROADCAST_WARP_DATA
                real3 delta = make_real3(posq2.x-posq1.x, posq2.y-posq1.y, posq2.z-posq1.z);
#ifdef USE_PERIODIC
                APPLY_PERIODIC_TO_DELTA(delta)
#endif
                real r2 = delta.x*delta.x + delta.y*delta.y + delta.z*delta.z;
                real invR = RSQRT(r2);
                real r = r2*invR;
                LOAD_ATOM2_PARAMETERS
                atom2 = y*TILE_SIZE+j;
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
            }
#endif
        }
        else {
            // This is an off-diagonal tile.
            unsigned int j = y*TILE_SIZE + tgx;
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
            unsigned int tj = tgx;
            for (j = 0; j < TILE_SIZE; j++) {
                int atom2 = tbx+tj;
                real4 posq2 = shflPosq;
                real3 delta = make_real3(posq2.x-posq1.x, posq2.y-posq1.y, posq2.z-posq1.z);
#ifdef USE_PERIODIC
                APPLY_PERIODIC_TO_DELTA(delta)
#endif
                real r2 = delta.x*delta.x + delta.y*delta.y + delta.z*delta.z;
                real invR = RSQRT(r2);
                real r = r2*invR;
                LOAD_ATOM2_PARAMETERS
                atom2 = y*TILE_SIZE+tj;
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
                SHUFFLE_WARP_DATA
#ifdef USE_EXCLUSIONS
                excl >>= 1;
#endif
                // cycles the indices
                // 0 1 2 3 4 5 6 7 -> 1 2 3 4 5 6 7 0
                tj = (tj + 1) & (TILE_SIZE - 1);
            }
            const unsigned int offset = y*TILE_SIZE + tgx;
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

    TileWorkQueue workQueue(warp, totalWarps, numTiles);
    while (workQueue.claimNextBatch(tgx, workUnitCount)) {
        int pos = workQueue.start;
        int nextTileX = tiles[pos];
        while (pos < workQueue.end) {
            // Load atom1 data for block x and accumulate its forces across consecutive tiles sharing x.

            const bool hasExclusions = false;
            int x = nextTileX;
            unsigned int atom1 = x*TILE_SIZE + tgx;
            real3 force = make_real3(0);
            real4 posq1 = posq[atom1];
            LOAD_ATOM1_PARAMETERS
            const real4 blockCenterX = blockCenter[x];
            const real4 blockSizeX = blockSize[x];
            const bool singlePeriodicCopy = canUseSinglePeriodicCopy(blockSizeX, periodicBoxSize);
#ifdef USE_PERIODIC
            if (singlePeriodicCopy) {
                APPLY_PERIODIC_TO_POS_WITH_CENTER(posq1, blockCenterX)
            }
#endif
            do {
                unsigned int localAtom2 = interactingAtoms[pos*TILE_SIZE+tgx];
                unsigned int shuffledAtom2 = localAtom2;
                DECLARE_LOCAL_PARAMETERS
                real4 shflPosq;
                real3 shflForce = make_real3(0);
                if (localAtom2 < NUM_ATOMS) {
                    unsigned int j = localAtom2;
                    shflPosq = posq[j];
                    LOAD_LOCAL_PARAMETERS_FROM_GLOBAL
                }
                else {
                    shflPosq = make_real4(0);
                    CLEAR_LOCAL_PARAMETERS
                }

                // When block x fits in a single periodic image, pre-filter atom2 against block x's
                // cutoff box and compact active lanes into a smaller power-of-two sub-tile.

                unsigned int subTileSize = TILE_SIZE;
                if (singlePeriodicCopy) {
                    bool atomInCutoff = false;
                    if (localAtom2 < NUM_ATOMS) {
#ifdef USE_PERIODIC
                        APPLY_PERIODIC_TO_POS_WITH_CENTER(shflPosq, blockCenterX)
#endif
                        atomInCutoff = isAtomNearBlock(shflPosq, blockCenterX, blockSizeX);
                    }
                    const unsigned int atomsInCutoffMask = BALLOT(atomInCutoff);
                    if (atomsInCutoffMask == 0u) {
                        subTileSize = 0u;
                        shuffledAtom2 = NUM_ATOMS;
                    }
                    else {
                        const unsigned int numAtomsInCutoff = __popc(atomsInCutoffMask);
                        subTileSize = selectSubTileSize(numAtomsInCutoff);
                        if (subTileSize < TILE_SIZE) {
                            const unsigned int sourceLane = findSubTileSourceLane(tgx, subTileSize, atomsInCutoffMask, numAtomsInCutoff);
                            unsigned int selectedAtom2 = SHFL((int) shuffledAtom2, sourceLane);
                            shflPosq = real_shfl4(shflPosq, sourceLane);
                            SHUFFLE_LOCAL_PARAMETERS
                            shuffledAtom2 = (tgx < numAtomsInCutoff ? selectedAtom2 : NUM_ATOMS);
                        }
                        else if (!atomInCutoff) {
                            shuffledAtom2 = NUM_ATOMS;
                        }
                    }
                }

                // Evaluate interactions across the compacted sub-tile (or the full tile when periodic wrapping is per-pair).

                #pragma unroll 4
                for (unsigned int j = 0; j < subTileSize; j++) {
                    int atom2 = SHFL((int) shuffledAtom2, (tgx+j) & (subTileSize-1u));
                    real3 delta = make_real3(shflPosq.x-posq1.x, shflPosq.y-posq1.y, shflPosq.z-posq1.z);
#ifdef USE_PERIODIC
                    if (!singlePeriodicCopy) {
                        APPLY_PERIODIC_TO_DELTA(delta)
                    }
#endif
                    real r2 = delta.x*delta.x + delta.y*delta.y + delta.z*delta.z;
                    if (atom1 < NUM_ATOMS && atom2 < NUM_ATOMS && r2 < MAX_CUTOFF_SQUARED) {
                        real4 posq2 = shflPosq;
                        real invR = RSQRT(r2);
                        real r = r2*invR;
                        LOAD_ATOM2_PARAMETERS
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
                        energy += tempEnergy;
#ifdef INCLUDE_FORCES
#ifdef USE_SYMMETRIC
                        accumulatePairForces(delta, dEdR, force, shflForce);
#else
                        accumulatePairForces(dEdR1, dEdR2, force, shflForce);
#endif
#endif
                    }
                    SHUFFLE_WARP_DATA
                }
#ifdef INCLUDE_FORCES
                if (subTileSize > 0u) {
                    shflForce = reduceSubTileForce(shflForce, tgx, subTileSize);
                    if (shuffledAtom2 < NUM_ATOMS)
                        saveSingleForce(shuffledAtom2, shflForce, forceBuffers);
                }
#endif
                pos++;
                if (pos >= workQueue.end)
                    break;
                nextTileX = tiles[pos];
            } while (nextTileX == x);
#ifdef INCLUDE_FORCES
            if (atom1 < NUM_ATOMS)
                saveSingleForce(atom1, force, forceBuffers);
#endif
        }
    }
#else
    int pos = (int) (startTileIndex+warp*numTileIndices/totalWarps);
    int end = (int) (startTileIndex+(warp+1)*numTileIndices/totalWarps);
    int skipBase = 0;
    int currentSkipIndex = tbx;
    __shared__ volatile int skipTiles[THREAD_BLOCK_SIZE];
    skipTiles[threadIdx.x] = -1;
    // atomIndices can probably be shuffled as well
    // but it probably wouldn't make things any faster
    __shared__ int atomIndices[THREAD_BLOCK_SIZE];
    
    while (pos < end) {
        const bool hasExclusions = false;
        real3 force = make_real3(0);
        bool includeTile = true;

        // Extract the coordinates of this tile.
        int x, y;

        y = (int) floor(NUM_BLOCKS+0.5f-SQRT((NUM_BLOCKS+0.5f)*(NUM_BLOCKS+0.5f)-2*pos));
        x = (pos-y*NUM_BLOCKS+y*(y+1)/2);
        if (x < y || x >= NUM_BLOCKS) { // Occasionally happens due to roundoff error.
            y += (x < y ? -1 : 1);
            x = (pos-y*NUM_BLOCKS+y*(y+1)/2);
        }

        // Skip over tiles that have exclusions, since they were already processed.

        while (skipTiles[tbx+TILE_SIZE-1] < pos) {
            if (skipBase+tgx < NUM_TILES_WITH_EXCLUSIONS) {
                int2 tile = exclusionTiles[skipBase+tgx];
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
        if (includeTile) {
            unsigned int atom1 = x*TILE_SIZE + tgx;
            // Load atom data for this tile.
            real4 posq1 = posq[atom1];
            LOAD_ATOM1_PARAMETERS
            unsigned int j = y*TILE_SIZE + tgx;
            atomIndices[threadIdx.x] = j;
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
            {
                // We need to apply periodic boundary conditions separately for each interaction.
                unsigned int tj = tgx;
                for (j = 0; j < TILE_SIZE; j++) {
                    int atom2 = tbx+tj;
                    real4 posq2 = shflPosq;
                    real3 delta = make_real3(posq2.x-posq1.x, posq2.y-posq1.y, posq2.z-posq1.z);
#ifdef USE_PERIODIC
                    APPLY_PERIODIC_TO_DELTA(delta)
#endif
                    real r2 = delta.x*delta.x + delta.y*delta.y + delta.z*delta.z;
                    real invR = RSQRT(r2);
                    real r = r2*invR;
                    LOAD_ATOM2_PARAMETERS
                    atom2 = atomIndices[tbx+tj];
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
                    SHUFFLE_WARP_DATA
                    tj = (tj + 1) & (TILE_SIZE - 1);
                }
            }

            // Write results.
#ifdef INCLUDE_FORCES
            atomicAdd(&forceBuffers[atom1], static_cast<unsigned long long>(realToFixedPoint(force.x)));
            atomicAdd(&forceBuffers[atom1+PADDED_NUM_ATOMS], static_cast<unsigned long long>(realToFixedPoint(force.y)));
            atomicAdd(&forceBuffers[atom1+2*PADDED_NUM_ATOMS], static_cast<unsigned long long>(realToFixedPoint(force.z)));
            unsigned int atom2 = y*TILE_SIZE + tgx;
            if (atom2 < PADDED_NUM_ATOMS) {
                atomicAdd(&forceBuffers[atom2], static_cast<unsigned long long>(realToFixedPoint(shflForce.x)));
                atomicAdd(&forceBuffers[atom2+PADDED_NUM_ATOMS], static_cast<unsigned long long>(realToFixedPoint(shflForce.y)));
                atomicAdd(&forceBuffers[atom2+2*PADDED_NUM_ATOMS], static_cast<unsigned long long>(realToFixedPoint(shflForce.z)));
            }
#endif
        }
        pos++;
    }
#endif
    
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
            LOAD_ATOM2_PARAMETERS_FROM_GLOBAL
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
    }
#endif
#ifdef INCLUDE_ENERGY
    energyBuffer[blockIdx.x*blockDim.x+threadIdx.x] += energy;
#endif
    SAVE_DERIVATIVES
}
