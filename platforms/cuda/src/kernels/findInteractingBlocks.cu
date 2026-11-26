#define GROUP_SIZE 256
#define BUFFER_SIZE 256

/**
 * To use half precision, we're supposed to include cuda_fp16.h.  Unfortunately,
 * it isn't included in the search path automatically, and there's no reliable
 * way to find where it's located on disk.  Instead we provide our own definitions
 * for the few symbols we need.
 */
struct __align__(2) __half {
    unsigned short x;
};
__device__ __half __float2half_ru(const float f) {
    __half h;
    asm("{cvt.rp.f16.f32 %0, %1;}" : "=h"(*reinterpret_cast<unsigned short *>(&h)) : "f"(f));
    return h;
}
__device__ float __half2float(const __half h) {
    float f;
    asm("{cvt.f32.f16 %0, %1;}" : "=f"(f) : "h"(*reinterpret_cast<const unsigned short *>(&h)));
    return f;
}
struct half3 {
    __device__ half3(real3 f) {
        // Round up so we'll err on the side of making the box a little too large.
        // This ensures interactions will never be missed.
        v[0] = __float2half_ru((float) f.x);
        v[1] = __float2half_ru((float) f.y);
        v[2] = __float2half_ru((float) f.z);
    }
    __device__ real3 toReal3() const {
        return make_real3(__half2float(v[0]), __half2float(v[1]), __half2float(v[2]));
    }
private:
    __half v[3];
};

#ifdef USE_LARGE_BLOCKS
enum { REBUILD_FLAG_INDEX = 0, LAST_REBUILD_STEP_INDEX = 1, NUM_DISPERSED_BLOCKS_INDEX = 2 };

/**
 * Broadcast the x, y, and z coordinates of a position from one lane across the warp.
 */
static __inline__ __device__ real4 broadcastWarpPos(real4 pos, int sourceLane) {
    return make_real4(
        __shfl_sync(0xffffffff, pos.x, sourceLane),
        __shfl_sync(0xffffffff, pos.y, sourceLane),
        __shfl_sync(0xffffffff, pos.z, sourceLane), 0);
}

/**
 * Reduce a scalar value across all 32 lanes of a warp to compute max_i(val_i).
 */
static __inline__ __device__ real reduceWarpMax(real val) {
    for (int offset = 16; offset > 0; offset >>= 1)
        val = max(val, __shfl_xor_sync(0xffffffff, val, offset));
    return val;
}

/**
 * Reduce per-lane half-extents across a warp to compute componentwise max_i(val_i).
 */
static __inline__ __device__ real4 reduceWarpMax3(real4 val) {
    for (int offset = 16; offset > 0; offset >>= 1) {
        val.x = max(val.x, __shfl_xor_sync(0xffffffff, val.x, offset));
        val.y = max(val.y, __shfl_xor_sync(0xffffffff, val.y, offset));
        val.z = max(val.z, __shfl_xor_sync(0xffffffff, val.z, offset));
    }
    return val;
}

/**
 * Reduce per-lane positions across a warp to compute minPos = min_i(pos_i) and maxPos = max_i(pos_i).
 */
static __inline__ __device__ void reduceWarpBoundingBox(real4& minPos, real4& maxPos) {
    for (int offset = 16; offset > 0; offset >>= 1) {
        minPos.x = min(minPos.x, __shfl_xor_sync(0xffffffff, minPos.x, offset));
        minPos.y = min(minPos.y, __shfl_xor_sync(0xffffffff, minPos.y, offset));
        minPos.z = min(minPos.z, __shfl_xor_sync(0xffffffff, minPos.z, offset));
        maxPos.x = max(maxPos.x, __shfl_xor_sync(0xffffffff, maxPos.x, offset));
        maxPos.y = max(maxPos.y, __shfl_xor_sync(0xffffffff, maxPos.y, offset));
        maxPos.z = max(maxPos.z, __shfl_xor_sync(0xffffffff, maxPos.z, offset));
    }
}

/**
 * Check whether any atom has moved more than 0.5*PADDING since the last rebuild
 * (|oldPositions - pos|^2 > 0.25*PADDING^2), and update each block's half-extents
 * blockSize = max_i |periodic(pos_i - blockCenter)| on non-rebuild steps. When a rebuild is
 * needed, rebuildNeighborList[LAST_REBUILD_STEP_INDEX] is set to neighborListStep to trigger
 * sortBoxData and findBlocksWithInteractions.
 */
extern "C" __global__ void findBlockBounds(int numAtoms, real4 periodicBoxSize, real4 invPeriodicBoxSize, real4 periodicBoxVecX, real4 periodicBoxVecY, real4 periodicBoxVecZ,
        const real4* __restrict__ posq, real4* __restrict__ blockCenter, real4* __restrict__ blockBoundingBox, int* __restrict__ rebuildNeighborList,
        const real4* __restrict__ oldPositions, unsigned int* __restrict__ interactionCount, bool forceRebuild, unsigned int neighborListStep) {
    const int indexInWarp = threadIdx.x%32;
    bool resetCounts = (forceRebuild && blockIdx.x == 0 && threadIdx.x == 0);
    if (!forceRebuild) {
        // Check atom displacements since the last rebuild and update each block's bounding box.

        const int warpIndex = (blockIdx.x*blockDim.x+threadIdx.x)/32;
        const int totalWarps = (blockDim.x*gridDim.x)/32;
        bool rebuild = false;
        for (int blockIndex = warpIndex; blockIndex < NUM_BLOCKS; blockIndex += totalWarps) {
            int atomIndex = blockIndex*TILE_SIZE+indexInWarp;
            bool isValidAtom = (atomIndex < numAtoms);
            real4 pos = (isValidAtom ? posq[atomIndex] : make_real4(0, 0, 0, 0));
            if (isValidAtom) {
                real4 delta = oldPositions[atomIndex]-pos;
                rebuild |= (delta.x*delta.x + delta.y*delta.y + delta.z*delta.z > 0.25f*PADDING*PADDING);
            }
            real4 center = broadcastWarpPos(indexInWarp == 0 ? blockCenter[blockIndex] : make_real4(0, 0, 0, 0), 0);
            real4 blockSize = make_real4(0, 0, 0, 0);
            if (isValidAtom) {
#ifdef USE_PERIODIC
                APPLY_PERIODIC_TO_POS_WITH_CENTER(pos, center)
#endif
                blockSize = make_real4(FABS(pos.x-center.x), FABS(pos.y-center.y), FABS(pos.z-center.z), 0);
            }
            blockSize = reduceWarpMax3(blockSize);
            if (indexInWarp == 0)
                blockBoundingBox[blockIndex] = blockSize;
        }
        int rebuildFlags = BALLOT(rebuild);
        resetCounts = (rebuildFlags != 0 && indexInWarp == 0);
    }

    // Record that the neighbor list needs to be rebuilt on this step.

    if (resetCounts) {
        rebuildNeighborList[LAST_REBUILD_STEP_INDEX] = (int) neighborListStep;
        rebuildNeighborList[NUM_DISPERSED_BLOCKS_INDEX] = 0;
        interactionCount[0] = 0;
        interactionCount[1] = 0;
    }
}

/**
 * Compute the axis-aligned bounding box and bounding-sphere radius for a 32-atom block,
 * recording the reference positions in oldPositions. Under periodic boundary conditions,
 * a second pass re-wraps atoms around the initial midpoint center = 0.5*(maxPos + minPos)
 * so blocks straddling a cell boundary obtain a tight half-extent blockSize = 0.5*(maxPos - minPos)
 * and sphere radius center.w = max_i |pos_i - center|.
 */
static __inline__ __device__ void computeAtomBlockBounds(int blockIndex, int indexInWarp,
        real4 periodicBoxSize, real4 invPeriodicBoxSize, real4 periodicBoxVecX, real4 periodicBoxVecY, real4 periodicBoxVecZ,
        const real4* __restrict__ posq, real4* __restrict__ oldPositions, real4& center, real4& blockSize) {
    int atomIndex = blockIndex*TILE_SIZE + indexInWarp;
    bool isValidAtom = (atomIndex < NUM_ATOMS);
    real4 pos = (isValidAtom ? posq[atomIndex] : make_real4(0, 0, 0, 0));
    if (isValidAtom)
        oldPositions[atomIndex] = pos;
#ifdef USE_PERIODIC
    if (indexInWarp == 0) {
        APPLY_PERIODIC_TO_POS(pos)
    }
#endif
    real4 firstAtomPos = broadcastWarpPos(pos, 0);
    if (!isValidAtom)
        pos = firstAtomPos;
#ifdef USE_PERIODIC
    else if (indexInWarp > 0) {
        APPLY_PERIODIC_TO_POS_WITH_CENTER(pos, firstAtomPos)
    }
#endif
    real4 minPos = pos, maxPos = pos;
    reduceWarpBoundingBox(minPos, maxPos);
    center = 0.5f*(maxPos+minPos);
#ifdef USE_PERIODIC
    if (isValidAtom) {
        APPLY_PERIODIC_TO_POS_WITH_CENTER(pos, center)
    }
    else
        pos = center;
    minPos = maxPos = pos;
    reduceWarpBoundingBox(minPos, maxPos);
    center = 0.5f*(maxPos+minPos);
#endif
    blockSize = 0.5f*(maxPos-minPos);
    real4 delta = pos-center;
    real maxRadiusSquared = (isValidAtom ? delta.x*delta.x + delta.y*delta.y + delta.z*delta.z : (real) 0);
    center.w = SQRT(reduceWarpMax(maxRadiusSquared));
}

/**
 * Return true if a block's half-extent exceeds 1.65*PADDED_CUTOFF or cannot fit within a single
 * periodic image (0.5*periodicBoxSize - roundedSize < PADDED_CUTOFF), so it is excluded from its
 * superblock bounding box and tested individually via dispersedBlockIndices.
 */
static __inline__ __device__ bool isDispersedAtomBlock(real4 roundedSize, real4 periodicBoxSize) {
    const real compactBlockLimit = 1.65f*PADDED_CUTOFF;
    if (roundedSize.x > compactBlockLimit || roundedSize.y > compactBlockLimit || roundedSize.z > compactBlockLimit)
        return true;
#ifdef USE_PERIODIC
    if (0.5f*periodicBoxSize.x-roundedSize.x < PADDED_CUTOFF ||
            0.5f*periodicBoxSize.y-roundedSize.y < PADDED_CUTOFF ||
            0.5f*periodicBoxSize.z-roundedSize.z < PADDED_CUTOFF)
        return true;
#endif
    return false;
}

/**
 * On rebuild steps, compute bounding boxes and bounding spheres for each 32-atom block and
 * aggregate groups of 32 compact blocks into superblock bounding boxes (largeBlockCenter and
 * largeBlockBoundingBox). Blocks classified as dispersed by isDispersedAtomBlock are omitted from
 * compactBlockMasks and appended to dispersedBlockIndices so wide blocks do not inflate superblock volumes.
 */
extern "C" __global__ void sortBoxData(unsigned int* __restrict__ sortedBlocks, real4* __restrict__ blockCenter,
        real4* __restrict__ blockBoundingBox, half3* __restrict__ sortedBlockBoundingBox,
        real4* __restrict__ largeBlockCenter, half3* __restrict__ largeBlockBoundingBox, real4 periodicBoxSize,
        real4 invPeriodicBoxSize, real4 periodicBoxVecX, real4 periodicBoxVecY, real4 periodicBoxVecZ,
        const real4* __restrict__ posq, real4* __restrict__ oldPositions,
        int* __restrict__ rebuildNeighborList, unsigned int neighborListStep) {
    if ((unsigned int) rebuildNeighborList[LAST_REBUILD_STEP_INDEX] != neighborListStep) {
        if (blockIdx.x == 0 && threadIdx.x == 0)
            rebuildNeighborList[REBUILD_FLAG_INDEX] = 0;
        return;
    }
    if (blockIdx.x == 0 && threadIdx.x == 0)
        rebuildNeighborList[REBUILD_FLAG_INDEX] = 1;
    const int indexInWarp = threadIdx.x%32;
    const int warpInBlock = threadIdx.x/32;
    const int numSuperblocks = (NUM_BLOCKS+31)/32;
    unsigned int* compactBlockMasks = sortedBlocks;
    unsigned int* dispersedBlockIndices = compactBlockMasks+numSuperblocks;
    __shared__ real4 sharedBlockCenter[32];
    __shared__ real4 sharedBlockSize[32];
    for (int superblockIndex = blockIdx.x; superblockIndex < numSuperblocks; superblockIndex += gridDim.x) {
        // Compute the bounding box and bounding sphere radius for each block in this superblock.

        for (int pass = 0; pass < 4; pass++) {
            int blockInSuperblock = pass*8 + warpInBlock;
            int blockIndex = superblockIndex*32 + blockInSuperblock;
            if (blockIndex < NUM_BLOCKS) {
                real4 center, blockSize;
                computeAtomBlockBounds(blockIndex, indexInWarp, periodicBoxSize, invPeriodicBoxSize,
                        periodicBoxVecX, periodicBoxVecY, periodicBoxVecZ, posq, oldPositions, center, blockSize);
                if (indexInWarp == 0) {
                    half3 roundedBox(trimTo3(blockSize));
                    real3 roundedSize = roundedBox.toReal3();
                    blockBoundingBox[blockIndex] = blockSize;
                    blockCenter[blockIndex] = center;
                    sortedBlockBoundingBox[blockIndex] = roundedBox;
                    sharedBlockCenter[blockInSuperblock] = center;
                    sharedBlockSize[blockInSuperblock] = make_real4(roundedSize.x, roundedSize.y, roundedSize.z, 0);
                }
            }
        }
        __syncthreads();

        // Combine compact blocks into a superblock bounding box and append wide blocks to the dispersed list.

        if (warpInBlock == 0) {
            int blockIndex = superblockIndex*32 + indexInWarp;
            bool isValidBlock = (blockIndex < NUM_BLOCKS);
            real4 center = (isValidBlock ? sharedBlockCenter[indexInWarp] : make_real4(0, 0, 0, 0));
            real4 roundedSize = (isValidBlock ? sharedBlockSize[indexInWarp] : make_real4(0, 0, 0, 0));
            bool isDispersedBlock = (isValidBlock && isDispersedAtomBlock(roundedSize, periodicBoxSize));
            bool isCompactBlock = (isValidBlock && !isDispersedBlock);
            unsigned int compactMask = BALLOT(isCompactBlock);
            unsigned int dispersedMask = BALLOT(isDispersedBlock);
            int firstCompactLane = (compactMask != 0 ? __ffs(compactMask)-1 : 0);
            real4 firstCompactCenter = broadcastWarpPos(center, firstCompactLane);
            if (!isCompactBlock) {
                center = firstCompactCenter;
                roundedSize = make_real4(0, 0, 0, 0);
            }
#ifdef USE_PERIODIC
            else if (indexInWarp != firstCompactLane) {
                APPLY_PERIODIC_TO_POS_WITH_CENTER(center, firstCompactCenter)
            }
#endif
            real4 minPos = center-roundedSize;
            real4 maxPos = center+roundedSize;
            reduceWarpBoundingBox(minPos, maxPos);
            real4 superblockCenter = 0.5f*(maxPos+minPos);
#ifdef USE_PERIODIC
            if (isCompactBlock) {
                APPLY_PERIODIC_TO_POS_WITH_CENTER(center, superblockCenter)
            }
            else
                center = superblockCenter;
            minPos = center-roundedSize;
            maxPos = center+roundedSize;
            reduceWarpBoundingBox(minPos, maxPos);
            superblockCenter = 0.5f*(maxPos+minPos);
#endif
            if (indexInWarp == 0) {
                largeBlockCenter[superblockIndex] = superblockCenter;
                largeBlockBoundingBox[superblockIndex] = half3(trimTo3(0.5f*(maxPos-minPos)));
                compactBlockMasks[superblockIndex] = compactMask;
            }
            if (dispersedMask != 0) {
                unsigned int dispersedOffset = (indexInWarp == 0 ? atomicAdd((unsigned int*) &rebuildNeighborList[NUM_DISPERSED_BLOCKS_INDEX], (unsigned int) __popc(dispersedMask)) : 0);
                dispersedOffset = __shfl_sync(0xffffffff, dispersedOffset, 0);
                if (isDispersedBlock)
                    dispersedBlockIndices[dispersedOffset + __popc(dispersedMask & ((1u<<indexInWarp)-1u))] = blockIndex;
            }
        }
        __syncthreads();
    }
}

/**
 * Compute the squared minimum distance sum_a max(0, |delta_a| - sizeX_a - sizeY_a)^2 between
 * two axis-aligned bounding boxes separated by center-to-center vector delta with half-extents
 * sizeX and sizeY.
 */
static __inline__ __device__ real boxDistanceSquared(real3 delta, real3 sizeX, real3 sizeY) {
    real dx = max((real) 0, FABS(delta.x)-sizeX.x-sizeY.x);
    real dy = max((real) 0, FABS(delta.y)-sizeX.y-sizeY.y);
    real dz = max((real) 0, FABS(delta.z)-sizeX.z-sizeY.z);
    return dx*dx + dy*dy + dz*dz;
}

/**
 * Return true if two bounding boxes in a triclinic cell are wide enough along y or z
 * (0.5*periodicBoxSize - sizeX - sizeY < PADDED_CUTOFF) that center-to-center minimum-image
 * wrapping could miss a nearer periodic copy.
 */
static __inline__ __device__ bool triclinicBoxesMayWrap(real4 periodicBoxSize, real3 sizeX, real3 sizeY) {
#ifdef TRICLINIC
    return (0.5f*periodicBoxSize.z-sizeX.z-sizeY.z < PADDED_CUTOFF ||
            0.5f*periodicBoxSize.y-sizeX.y-sizeY.y < PADDED_CUTOFF);
#else
    return false;
#endif
}

/**
 * Test whether an atom separated from a block center by delta lies within PADDED_CUTOFF of
 * both the block's bounding sphere (|delta|^2 < (PADDED_CUTOFF + center.w)^2) and its
 * axis-aligned bounding box (sum_a max(0, |delta_a| - size_a)^2 < PADDED_CUTOFF_SQUARED).
 */
static __inline__ __device__ bool atomNearBlockBounds(real3 delta, real4 center, real3 size) {
    real sphereCutoff = PADDED_CUTOFF + center.w;
    real dx = max((real) 0, FABS(delta.x)-size.x);
    real dy = max((real) 0, FABS(delta.y)-size.y);
    real dz = max((real) 0, FABS(delta.z)-size.z);
    return (delta.x*delta.x + delta.y*delta.y + delta.z*delta.z < sphereCutoff*sphereCutoff &&
            dx*dx + dy*dy + dz*dz < PADDED_CUTOFF_SQUARED);
}

/**
 * Return the 32-bit mask of atoms j in block X (stored in posBuffer[warpStart..warpStart+31]
 * with pos.w = 0.5*(x^2 + y^2 + z^2)) located within PADDED_CUTOFF of pos2. When all atoms
 * fit in a single periodic image around blockCenterX, distances are evaluated via
 * 0.5*|posj - pos2|^2 = posj.w + pos2.w - dot(posj, pos2) < 0.5*PADDED_CUTOFF_SQUARED.
 */
static __inline__ __device__ int findInteractingAtom1Mask(real4 pos2, const real4* posBuffer, int warpStart,
        int atom1InRangeFlags, unsigned int validAtom1Mask, real4 blockCenterX, real3 blockSizeX,
        bool singlePeriodicCopy, real4 periodicBoxSize, real4 invPeriodicBoxSize,
        real4 periodicBoxVecX, real4 periodicBoxVecY, real4 periodicBoxVecZ) {
    int interacts = 0;
#ifdef USE_PERIODIC
    if (!singlePeriodicCopy) {
        int first = __ffs(atom1InRangeFlags)-1;
        int last = 32-__clz(atom1InRangeFlags);
        for (int j = first; j < last; j++) {
            real3 delta = trimTo3(pos2)-trimTo3(posBuffer[warpStart+j]);
            APPLY_PERIODIC_TO_DELTA(delta)
            interacts |= (delta.x*delta.x+delta.y*delta.y+delta.z*delta.z < PADDED_CUTOFF_SQUARED ? 1<<j : 0);
        }
        return interacts & (int) validAtom1Mask;
    }
#endif
    if (atomNearBlockBounds(trimTo3(pos2)-trimTo3(blockCenterX), blockCenterX, blockSizeX)) {
        for (int half = 0; half < 32; half += 16) {
            if ((atom1InRangeFlags & (0xffff << half)) != 0) {
                for (int j = half; j < half+16; j++) {
                    real4 posj = posBuffer[warpStart+j];
                    real halfDist2 = posj.w + pos2.w - posj.x*pos2.x - posj.y*pos2.y - posj.z*pos2.z;
                    interacts |= (halfDist2 < 0.5f*PADDED_CUTOFF_SQUARED ? 1<<j : 0);
                }
            }
        }
    }
    return interacts & (int) validAtom1Mask;
}

/**
 * Extract buffered neighbor atoms with at most MAX_BITS_FOR_PAIRS interacting partners in block x
 * into singlePairs, and compact the remaining dense neighbor entries in atoms and flags.
 */
static __inline__ __device__ int saveSinglePairs(int x, int* atoms, int* flags, int length, unsigned int maxSinglePairs, unsigned int* singlePairCount, int2* singlePairs) {
    // Record interactions that should be computed as single pairs rather than in blocks.

    const int indexInWarp = threadIdx.x%32;
    int sum = 0;
    #pragma unroll 8 // (GROUP_SIZE / TILE_SIZE)
    for (int i = indexInWarp; i < length; i += 32) {
        int count = __popc(flags[i]);
        sum += (count <= MAX_BITS_FOR_PAIRS ? count : 0);
    }
    for (int i = 1; i < 32; i *= 2) {
        int n = __shfl_up_sync(0xffffffff, sum, i);
        if (indexInWarp >= i)
            sum += n;
    }
    unsigned int pairStartIndex = (indexInWarp == 31 ? atomicAdd(singlePairCount, (unsigned int) sum) : 0);
    pairStartIndex = __shfl_sync(0xffffffff, pairStartIndex, 31);
    int prevSum = __shfl_up_sync(0xffffffff, sum, 1);
    unsigned int pairIndex = pairStartIndex + (indexInWarp > 0 ? prevSum : 0);
    for (int i = indexInWarp; i < length; i += 32) {
        int count = __popc(flags[i]);
        if (count <= MAX_BITS_FOR_PAIRS && pairIndex+count <= maxSinglePairs) {
            int f = flags[i];
            while (f != 0) {
                singlePairs[pairIndex] = make_int2(atoms[i], x*TILE_SIZE+__ffs(f)-1);
                f &= f-1;
                pairIndex++;
            }
        }
    }

    // Compact the remaining interactions.

    const int warpMask = (1<<indexInWarp)-1;
    int numCompacted = 0;
    for (int start = 0; start < length; start += 32) {
        int i = start+indexInWarp;
        int atom = (i < length ? atoms[i] : 0);
        int flag = (i < length ? flags[i] : 0);
        bool include = (i < length && __popc(flag) > MAX_BITS_FOR_PAIRS);
        int includeFlags = BALLOT(include);
        if (include) {
            int index = numCompacted+__popc(includeFlags&warpMask);
            atoms[index] = atom;
            flags[index] = flag;
        }
        numCompacted += __popc(includeFlags);
    }
    __syncwarp();
    return numCompacted;
}

/**
 * Reserve tilesToStore output slots in interactionCount[0] and write block x's buffered
 * neighbor atom indices to interactingTiles and interactingAtoms, padding partial tiles with NUM_ATOMS.
 */
static __inline__ __device__ void storeInteractionTiles(int x, const int* buffer, int neighborsInBuffer,
        unsigned int tilesToStore, bool padPartialTile, int indexInWarp, unsigned int maxTiles,
        unsigned int* __restrict__ interactionCount, int* __restrict__ interactingTiles,
        unsigned int* __restrict__ interactingAtoms) {
    unsigned int newTileStartIndex = (indexInWarp == 0 ? atomicAdd(&interactionCount[0], tilesToStore) : 0);
    newTileStartIndex = __shfl_sync(0xffffffff, newTileStartIndex, 0);
    if (newTileStartIndex+tilesToStore <= maxTiles) {
        if (indexInWarp < tilesToStore)
            interactingTiles[newTileStartIndex+indexInWarp] = x;
        for (int j = 0; j < tilesToStore; j++) {
            int offset = indexInWarp + j*TILE_SIZE;
            interactingAtoms[(newTileStartIndex+j)*TILE_SIZE+indexInWarp] = (!padPartialTile || offset < neighborsInBuffer ? buffer[offset] : NUM_ATOMS);
        }
    }
}

/**
 * Compare bounding boxes of 32-block superblocks and 32-atom blocks to construct the neighbor list
 * in two stages.
 *
 * STAGE 1:
 *
 * Coarse-grained superblock and atom-block filtering. Each warp loads a 32-atom block X (block1).
 * If X is compact, the warp tests 32 superblocks in parallel using boxDistanceSquared < PADDED_CUTOFF_SQUARED,
 * expands each interacting superblock into its 32 constituent compact blocks Y (block2 > block1),
 * and then tests all dispersed blocks Y; if X is dispersed, the warp tests only dispersed blocks Y
 * with block2 > block1 so every unordered block pair {block1, block2} is visited once. Candidate
 * blocks Y are filtered against X's bounding sphere, axis-aligned bounding box, and exclusion list.
 *
 * STAGE 2:
 *
 * Fine-grained atom-pair filtering. For each surviving block Y, the warp builds atom1InRangeFlags
 * for atoms in X within PADDED_CUTOFF of Y's bounding sphere and bounding box. Each lane then loads
 * one atom2 from Y and tests it against candidate atoms in X via findInteractingAtom1Mask. Interacting
 * atom2 indices are buffered in shared memory, sparse entries (at most MAX_BITS_FOR_PAIRS partners in X)
 * are routed to singlePairs, and tiles of 32 interacting atoms are written to global memory.
 *
 * [in] periodicBoxSize        - dimensions of the periodic box
 * [in] invPeriodicBoxSize     - inverse dimensions of the periodic box
 * [in] periodicBoxVecX        - first periodic box vector
 * [in] periodicBoxVecY        - second periodic box vector
 * [in] periodicBoxVecZ        - third periodic box vector
 * [out] interactionCount      - interactionCount[0] stores tile count; interactionCount[1] stores single-pair count
 * [out] interactingTiles      - block X index for each output tile
 * [out] interactingAtoms      - 32 interacting atom2 indices per output tile (padded with NUM_ATOMS)
 * [out] singlePairs           - sparse interacting (atom2, atom1) pairs
 * [in] posq                   - atom positions and charges
 * [in] maxTiles               - allocated capacity of interactingTiles
 * [in] maxSinglePairs         - allocated capacity of singlePairs
 * [in] startBlockIndex        - first block1 index assigned to this device
 * [in] numBlocks              - number of block1 entries assigned to this device
 * [in] sortedBlocks           - compactBlockMasks[0..numSuperblocks-1] followed by dispersedBlockIndices
 * [in] blockCenter            - center (x, y, z) and bounding-sphere radius (w) of each 32-atom block
 * [in] sortedBlockBoundingBox - half-precision half-extents of each 32-atom block
 * [in] largeBlockCenter       - center of each 32-block superblock
 * [in] largeBlockBoundingBox  - half-precision half-extents of each 32-block superblock
 * [in] exclusionIndices       - excluded block indices packed by block X
 * [in] exclusionRowIndices    - CSR row offsets into exclusionIndices for each block X
 * [in] oldPositions           - unused here (recorded in sortBoxData; retained for launch-signature compatibility)
 * [in] rebuildNeighborList    - rebuild flag and dispersed-block count written by findBlockBounds and sortBoxData
 */
extern "C" __global__ __launch_bounds__(GROUP_SIZE,3) void findBlocksWithInteractions(real4 periodicBoxSize, real4 invPeriodicBoxSize, real4 periodicBoxVecX, real4 periodicBoxVecY, real4 periodicBoxVecZ,
        unsigned int* __restrict__ interactionCount, int* __restrict__ interactingTiles, unsigned int* __restrict__ interactingAtoms,
        int2* __restrict__ singlePairs, const real4* __restrict__ posq, unsigned int maxTiles, unsigned int maxSinglePairs, unsigned int startBlockIndex,
        unsigned int numBlocks, const unsigned int* __restrict__ sortedBlocks, const real4* __restrict__ blockCenter, const half3* __restrict__ sortedBlockBoundingBox,
        const real4* __restrict__ largeBlockCenter, const half3* __restrict__ largeBlockBoundingBox,
        const unsigned int* __restrict__ exclusionIndices, const unsigned int* __restrict__ exclusionRowIndices,
        real4* __restrict__ oldPositions, const int* __restrict__ rebuildNeighborList) {

    if (rebuildNeighborList[REBUILD_FLAG_INDEX] == 0)
        return; // The neighbor list doesn't need to be rebuilt.

    const int indexInWarp = threadIdx.x%32;
    const int warpStart = threadIdx.x-indexInWarp;
    const int totalWarps = blockDim.x*gridDim.x/32;
    const int warpIndex = (blockIdx.x*blockDim.x+threadIdx.x)/32;
    const int warpMask = (1<<indexInWarp)-1;
    __shared__ int workgroupBuffer[BUFFER_SIZE*(GROUP_SIZE/32)];
    __shared__ int workgroupFlagsBuffer[BUFFER_SIZE*(GROUP_SIZE/32)];
    __shared__ int warpExclusions[MAX_EXCLUSIONS*(GROUP_SIZE/32)];
    __shared__ real4 posBuffer[GROUP_SIZE];
    int* buffer = workgroupBuffer+BUFFER_SIZE*(warpStart/32);
    int* flagsBuffer = workgroupFlagsBuffer+BUFFER_SIZE*(warpStart/32);
    int* exclusionsForX = warpExclusions+MAX_EXCLUSIONS*(warpStart/32);

    const int numSuperblocks = (NUM_BLOCKS+31)/32;
    const unsigned int* compactBlockMasks = sortedBlocks;
    const unsigned int* dispersedBlockIndices = compactBlockMasks+numSuperblocks;
    const int numDispersedBlocks = rebuildNeighborList[NUM_DISPERSED_BLOCKS_INDEX];

    // Loop over blocks.

    for (int block1 = startBlockIndex+warpIndex; block1 < startBlockIndex+numBlocks; block1 += totalWarps) {
        // Load data for this block.  Note that all threads in a warp are processing the same block.

        int x = block1;
        real4 blockCenterX = blockCenter[block1];
        real3 blockSizeX = sortedBlockBoundingBox[block1].toReal3();
        int neighborsInBuffer = 0;
        const unsigned int validAtom1Mask = BALLOT(x*TILE_SIZE+indexInWarp < NUM_ATOMS);
        real4 pos1 = (x*TILE_SIZE+indexInWarp < NUM_ATOMS ? posq[x*TILE_SIZE+indexInWarp] : make_real4(0, 0, 0, 0));
#ifdef USE_PERIODIC
        const bool singlePeriodicCopy = (0.5f*periodicBoxSize.x-blockSizeX.x >= PADDED_CUTOFF &&
                                         0.5f*periodicBoxSize.y-blockSizeX.y >= PADDED_CUTOFF &&
                                         0.5f*periodicBoxSize.z-blockSizeX.z >= PADDED_CUTOFF);
        if (singlePeriodicCopy) {
            APPLY_PERIODIC_TO_POS_WITH_CENTER(pos1, blockCenterX)
        }
#else
        const bool singlePeriodicCopy = true;
#endif
        pos1.w = 0.5f * (pos1.x * pos1.x + pos1.y * pos1.y + pos1.z * pos1.z);
        posBuffer[threadIdx.x] = pos1;

        // Load exclusion data for block x.

        const int exclusionStart = exclusionRowIndices[x];
        const int exclusionEnd = exclusionRowIndices[x+1];
        const int numExclusions = exclusionEnd-exclusionStart;
        #pragma unroll 4 // (MAX_EXCLUSIONS)
        for (int j = indexInWarp; j < numExclusions; j += 32)
            exclusionsForX[j] = exclusionIndices[exclusionStart+j];
        __syncwarp();

        // Compact blocks test compact superblocks (block2 > block1) followed by all dispersed blocks;
        // dispersed blocks test only dispersed blocks with block2 > block1, visiting each unordered pair once.

        const bool isCompactBlock1 = (((compactBlockMasks[block1/32] >> (block1%32)) & 1u) != 0);
        const int startSuperblock = (block1+1)/32;
        const int numCompactBatches = (isCompactBlock1 ? (numSuperblocks-startSuperblock+31)/32 : 0);
        const int totalBatches = numCompactBatches + (numDispersedBlocks+31)/32;
        int batch = -1;
        unsigned int laneCompactMask = 0;
        unsigned int superblockFlags = 0;
        while (true) {
            // Find the next batch of superblocks or dispersed blocks that may interact with block x.

            while (superblockFlags == 0 && ++batch < totalBatches) {
                if (batch < numCompactBatches) {
                    int superblockIndex = startSuperblock + batch*32 + indexInWarp;
                    bool includeSuperblock = false;
                    if (superblockIndex < numSuperblocks) {
                        laneCompactMask = compactBlockMasks[superblockIndex];
                        unsigned int afterBlock1Mask = (superblockIndex == block1/32 ? (0xfffffffeu << (block1%32)) : 0xffffffffu);
                        if ((laneCompactMask & afterBlock1Mask) != 0) {
                            real4 superblockCenter = largeBlockCenter[superblockIndex];
                            real3 superblockSize = largeBlockBoundingBox[superblockIndex].toReal3();
                            real4 blockDelta = blockCenterX-superblockCenter;
#ifdef USE_PERIODIC
                            APPLY_PERIODIC_TO_DELTA(blockDelta)
#endif
                            includeSuperblock = (boxDistanceSquared(trimTo3(blockDelta), blockSizeX, superblockSize) < PADDED_CUTOFF_SQUARED ||
                                                 triclinicBoxesMayWrap(periodicBoxSize, blockSizeX, superblockSize));
                        }
                    }
                    superblockFlags = BALLOT(includeSuperblock);
                }
                else {
                    superblockFlags = 1u;
                }
            }
            if (superblockFlags == 0)
                break;

            // Test the 32 candidate blocks in the selected superblock or dispersed batch.

            int superblockBit = __ffs(superblockFlags)-1;
            superblockFlags &= superblockFlags-1;
            int block2;
            bool includeBlock2;
            if (batch < numCompactBatches) {
                block2 = (startSuperblock + batch*32 + superblockBit)*32 + indexInWarp;
                unsigned int superblockCompactMask = __shfl_sync(0xffffffff, laneCompactMask, superblockBit);
                includeBlock2 = (block2 > block1 && block2 < NUM_BLOCKS && ((superblockCompactMask >> indexInWarp) & 1u) != 0);
            }
            else {
                int dispersedIndex = (batch-numCompactBatches)*32 + indexInWarp;
                block2 = (dispersedIndex < numDispersedBlocks ? (int) dispersedBlockIndices[dispersedIndex] : -1);
                includeBlock2 = (isCompactBlock1 ? (block2 >= 0) : (block2 > block1));
            }
            bool forceInclude = false;
            real4 candidateCenterY = make_real4(0, 0, 0, 0);
            real3 candidateSizeY = make_real3(0, 0, 0);
            if (includeBlock2) {
                candidateCenterY = blockCenter[block2];
                candidateSizeY = sortedBlockBoundingBox[block2].toReal3();
                real4 blockDelta = blockCenterX-candidateCenterY;
#ifdef USE_PERIODIC
                APPLY_PERIODIC_TO_DELTA(blockDelta)
#endif
                real sphereCutoff = PADDED_CUTOFF + blockCenterX.w + candidateCenterY.w;
                includeBlock2 &= (blockDelta.x*blockDelta.x+blockDelta.y*blockDelta.y+blockDelta.z*blockDelta.z < sphereCutoff*sphereCutoff);
                includeBlock2 &= (boxDistanceSquared(trimTo3(blockDelta), blockSizeX, candidateSizeY) < PADDED_CUTOFF_SQUARED);
                if (triclinicBoxesMayWrap(periodicBoxSize, blockSizeX, candidateSizeY))
                    includeBlock2 = forceInclude = true;
                if (includeBlock2) {
                    #pragma unroll 4 // (MAX_EXCLUSIONS)
                    for (int k = 0; k < numExclusions; k++)
                        includeBlock2 &= (exclusionsForX[k] != block2);
                }
            }

            // Loop over any blocks we identified as potentially containing neighbors.

            int includeBlockFlags = BALLOT(includeBlock2);
            int forceIncludeFlags = BALLOT(forceInclude);
            while (includeBlockFlags != 0) {
                int bitIndex = __ffs(includeBlockFlags)-1;
                includeBlockFlags &= includeBlockFlags-1;
                forceInclude = (forceIncludeFlags>>bitIndex) & 1;
                int y = __shfl_sync(0xffffffff, block2, bitIndex);
                real4 blockCenterY = make_real4(
                    __shfl_sync(0xffffffff, candidateCenterY.x, bitIndex),
                    __shfl_sync(0xffffffff, candidateCenterY.y, bitIndex),
                    __shfl_sync(0xffffffff, candidateCenterY.z, bitIndex),
                    __shfl_sync(0xffffffff, candidateCenterY.w, bitIndex));
                real3 blockSizeY = make_real3(
                    __shfl_sync(0xffffffff, candidateSizeY.x, bitIndex),
                    __shfl_sync(0xffffffff, candidateSizeY.y, bitIndex),
                    __shfl_sync(0xffffffff, candidateSizeY.z, bitIndex));
                real3 atomDelta = trimTo3(posBuffer[threadIdx.x])-trimTo3(blockCenterY);
#ifdef USE_PERIODIC
                APPLY_PERIODIC_TO_DELTA(atomDelta)
#endif
                int atom1InRangeFlags = BALLOT(forceInclude || atomNearBlockBounds(atomDelta, blockCenterY, blockSizeY));
                atom1InRangeFlags &= (int) validAtom1Mask;
                if (atom1InRangeFlags == 0)
                    continue;

                // Check each atom in block Y for interactions.

                int atom2 = y*TILE_SIZE+indexInWarp;
                real4 pos2 = (atom2 < NUM_ATOMS ? posq[atom2] : make_real4(0, 0, 0, 0));
#ifdef USE_PERIODIC
                if (singlePeriodicCopy) {
                    APPLY_PERIODIC_TO_POS_WITH_CENTER(pos2, blockCenterX)
                }
#endif
                pos2.w = 0.5f * (pos2.x * pos2.x + pos2.y * pos2.y + pos2.z * pos2.z);
                int interacts = (atom2 < NUM_ATOMS ? findInteractingAtom1Mask(pos2, posBuffer, warpStart, atom1InRangeFlags,
                        validAtom1Mask, blockCenterX, blockSizeX, singlePeriodicCopy, periodicBoxSize, invPeriodicBoxSize,
                        periodicBoxVecX, periodicBoxVecY, periodicBoxVecZ) : 0);

                // Add any interacting atoms to the buffer.

                int includeAtomFlags = BALLOT(interacts);
                if (interacts) {
                    int index = neighborsInBuffer+__popc(includeAtomFlags&warpMask);
                    buffer[index] = atom2;
                    flagsBuffer[index] = interacts;
                }
                __syncwarp();
                neighborsInBuffer += __popc(includeAtomFlags);
                if (neighborsInBuffer > BUFFER_SIZE-TILE_SIZE) {
                    // Store the new tiles to memory.

#if MAX_BITS_FOR_PAIRS > 0
                    neighborsInBuffer = saveSinglePairs(x, buffer, flagsBuffer, neighborsInBuffer, maxSinglePairs, &interactionCount[1], singlePairs);
#endif
                    unsigned int tilesToStore = neighborsInBuffer/TILE_SIZE;
                    if (tilesToStore > 0) {
                        storeInteractionTiles(x, buffer, neighborsInBuffer, tilesToStore, false, indexInWarp, maxTiles,
                                interactionCount, interactingTiles, interactingAtoms);
                        __syncwarp();
                        if (indexInWarp+TILE_SIZE*tilesToStore < BUFFER_SIZE) {
                            buffer[indexInWarp] = buffer[indexInWarp+TILE_SIZE*tilesToStore];
                            flagsBuffer[indexInWarp] = flagsBuffer[indexInWarp+TILE_SIZE*tilesToStore];
                        }
                        __syncwarp();
                        neighborsInBuffer -= TILE_SIZE*tilesToStore;
                    }
                }
            }
        }

        // If we have a partially filled buffer, store it to memory.

#if MAX_BITS_FOR_PAIRS > 0
        if (neighborsInBuffer > 32)
            neighborsInBuffer = saveSinglePairs(x, buffer, flagsBuffer, neighborsInBuffer, maxSinglePairs, &interactionCount[1], singlePairs);
#endif
        if (neighborsInBuffer > 0) {
            unsigned int tilesToStore = (neighborsInBuffer+TILE_SIZE-1)/TILE_SIZE;
            storeInteractionTiles(x, buffer, neighborsInBuffer, tilesToStore, true, indexInWarp, maxTiles,
                    interactionCount, interactingTiles, interactingAtoms);
        }
    }
}
#else
/**
 * Find a bounding box for the atoms in each block.
 */
extern "C" __global__ void findBlockBounds(int numAtoms, real4 periodicBoxSize, real4 invPeriodicBoxSize, real4 periodicBoxVecX, real4 periodicBoxVecY, real4 periodicBoxVecZ,
        const real4* __restrict__ posq, real4* __restrict__ blockCenter, real4* __restrict__ blockBoundingBox, int* __restrict__ rebuildNeighborList,
        real2* __restrict__ blockSizeRange) {
    int index = blockIdx.x*blockDim.x+threadIdx.x;
    int base = index*TILE_SIZE;
    real minSize = 1e38, maxSize = 0;
    while (base < numAtoms) {
        real4 pos = posq[base];
#ifdef USE_PERIODIC
        APPLY_PERIODIC_TO_POS(pos)
#endif
        real4 minPos = pos;
        real4 maxPos = pos;
        int last = min(base+TILE_SIZE, numAtoms);
        for (int i = base+1; i < last; i++) {
            pos = posq[i];
#ifdef USE_PERIODIC
            real4 center = 0.5f*(maxPos+minPos);
            APPLY_PERIODIC_TO_POS_WITH_CENTER(pos, center)
#endif
            minPos = make_real4(min(minPos.x,pos.x), min(minPos.y,pos.y), min(minPos.z,pos.z), 0);
            maxPos = make_real4(max(maxPos.x,pos.x), max(maxPos.y,pos.y), max(maxPos.z,pos.z), 0);
        }
        real4 blockSize = 0.5f*(maxPos-minPos);
        real4 center = 0.5f*(maxPos+minPos);
        center.w = 0;
        for (int i = base; i < last; i++) {
            pos = posq[i];
            real4 delta = posq[i]-center;
#ifdef USE_PERIODIC
            APPLY_PERIODIC_TO_DELTA(delta)
#endif
            center.w = max(center.w, delta.x*delta.x+delta.y*delta.y+delta.z*delta.z);
        }
        center.w = sqrt(center.w);
        blockBoundingBox[index] = blockSize;
        blockCenter[index] = center;
        real totalSize = blockSize.x+blockSize.y+blockSize.z;
        minSize = min(minSize, totalSize);
        maxSize = max(maxSize, totalSize);
        index += blockDim.x*gridDim.x;
        base = index*TILE_SIZE;
    }
    
    // Record the range of sizes seen by threads in this block.

    __shared__ real minBuffer[64], maxBuffer[64];
    minBuffer[threadIdx.x] = minSize;
    maxBuffer[threadIdx.x] = maxSize;
    __syncthreads();
    for (int step = 1; step < 64; step *= 2) {
        if (threadIdx.x+step < 64 && threadIdx.x%(2*step) == 0) {
            minBuffer[threadIdx.x] = min(minBuffer[threadIdx.x], minBuffer[threadIdx.x+step]);
            maxBuffer[threadIdx.x] = max(maxBuffer[threadIdx.x], maxBuffer[threadIdx.x+step]);
        }
        __syncthreads();
    }
    if (threadIdx.x == 0)
        blockSizeRange[blockIdx.x] = make_real2(minBuffer[0], maxBuffer[0]);
    if (blockIdx.x == 0 && threadIdx.x == 0)
        rebuildNeighborList[0] = 0;
}

extern "C" __global__ void computeSortKeys(const real4* __restrict__ blockBoundingBox, unsigned int* __restrict__ sortedBlocks, real2* __restrict__ blockSizeRange, int numSizes) {
    // Find the total range of sizes recorded by all blocks.

    __shared__ real2 sizeRange;
    if (threadIdx.x == 0) {
        sizeRange = blockSizeRange[0];
        for (int i = 1; i < numSizes; i++) {
            real2 size = blockSizeRange[i];
            if (size.x > 0)
                sizeRange.x = min(sizeRange.x, size.x);
            sizeRange.y = max(sizeRange.y, size.y);
        }
        sizeRange.x = LOG(sizeRange.x);
        sizeRange.y = LOG(sizeRange.y);
    }
    __syncthreads();

    // Sort keys store the bin in the high order part and the block in the low
    // order part.

    int numSizeBins = 20;
    real scale = numSizeBins/(sizeRange.y-sizeRange.x);
    for (unsigned int i = threadIdx.x+blockIdx.x*blockDim.x; i < NUM_BLOCKS; i += blockDim.x*gridDim.x) {
        real4 box = blockBoundingBox[i];
        real size = LOG(box.x+box.y+box.z);
        int bin = (size-sizeRange.x)*scale;
        bin = max(0, min(bin, numSizeBins-1));
        sortedBlocks[i] = (((unsigned int) bin)<<BIN_SHIFT) + i;
    }
}

/**
 * Sort the data about bounding boxes so it can be accessed more efficiently in the next kernel.
 */
extern "C" __global__ void sortBoxData(const unsigned int* __restrict__ sortedBlocks, const real4* __restrict__ blockCenter,
        const real4* __restrict__ blockBoundingBox, real4* __restrict__ sortedBlockCenter, half3* __restrict__ sortedBlockBoundingBox,
        const real4* __restrict__ posq, const real4* __restrict__ oldPositions,
        unsigned int* __restrict__ interactionCount, int* __restrict__ rebuildNeighborList, bool forceRebuild) {
    for (int i = threadIdx.x+blockIdx.x*blockDim.x; i < NUM_BLOCKS; i += blockDim.x*gridDim.x) {
        unsigned int index = sortedBlocks[i] & BLOCK_INDEX_MASK;
        sortedBlockCenter[i] = blockCenter[index];
        sortedBlockBoundingBox[i] = half3(trimTo3(blockBoundingBox[index]));

    }

    // Also check whether any atom has moved enough so that we really need to rebuild the neighbor list.

    bool rebuild = forceRebuild;
    for (int i = threadIdx.x+blockIdx.x*blockDim.x; i < NUM_ATOMS; i += blockDim.x*gridDim.x) {
        real4 delta = oldPositions[i]-posq[i];
        if (delta.x*delta.x + delta.y*delta.y + delta.z*delta.z > 0.25f*PADDING*PADDING)
            rebuild = true;
    }
    if (rebuild) {
        rebuildNeighborList[0] = 1;
        interactionCount[0] = 0;
        interactionCount[1] = 0;
    }
}

__device__ int saveSinglePairs(int x, int* atoms, int* flags, int length, unsigned int maxSinglePairs, unsigned int* singlePairCount, int2* singlePairs, int* sumBuffer, volatile unsigned int& pairStartIndex) {
    // Record interactions that should be computed as single pairs rather than in blocks.
    
    const int indexInWarp = threadIdx.x%32;
    int sum = 0;
    #pragma unroll 8 // (GROUP_SIZE / TILE_SIZE)
    for (int i = indexInWarp; i < length; i += 32) {
        int count = __popc(flags[i]);
        sum += (count <= MAX_BITS_FOR_PAIRS ? count : 0);
    }
    for (int i = 1; i < 32; i *= 2) {
        int n = __shfl_up_sync(0xffffffff, sum, i);
        if (indexInWarp >= i)
            sum += n;
    }
    if (indexInWarp == 31)
        pairStartIndex = atomicAdd(singlePairCount,(unsigned int) sum);
    __syncwarp();
    int prevSum = __shfl_up_sync(0xffffffff, sum, 1);
    unsigned int pairIndex = pairStartIndex + (indexInWarp > 0 ? prevSum : 0);
    for (int i = indexInWarp; i < length; i += 32) {
        int count = __popc(flags[i]);
        if (count <= MAX_BITS_FOR_PAIRS && pairIndex+count <= maxSinglePairs) {
            int f = flags[i];
            while (f != 0) {
                singlePairs[pairIndex] = make_int2(atoms[i], x*TILE_SIZE+__ffs(f)-1);
                f &= f-1;
                pairIndex++;
            }
        }
    }
    
    // Compact the remaining interactions.
    
    const int warpMask = (1<<indexInWarp)-1;
    int numCompacted = 0;
    for (int start = 0; start < length; start += 32) {
        int i = start+indexInWarp;
        int atom = atoms[i];
        int flag = flags[i];
        bool include = (i < length && __popc(flags[i]) > MAX_BITS_FOR_PAIRS);
        int includeFlags = BALLOT(include);
        if (include) {
            int index = numCompacted+__popc(includeFlags&warpMask);
            atoms[index] = atom;
            flags[index] = flag;
        }
        numCompacted += __popc(includeFlags);
    }
    return numCompacted;
}

/**
 * Compare the bounding boxes for each pair of atom blocks (comprised of 32 atoms each), forming a tile. If the two
 * atom blocks are sufficiently far apart, mark them as non-interacting. There are two stages in the algorithm.
 *
 * STAGE 1:
 *
 * A coarse grained atom block against interacting atom block neighbour list is constructed. 
 *
 * Each warp first loads in some block X of interest. Each thread within the warp then loads 
 * in a different atom block Y. If Y has exclusions with X, then Y is not processed.  If the bounding boxes 
 * of the two atom blocks are within the cutoff distance, then the two atom blocks are considered to be
 * interacting and Y is added to the buffer for X.
 *
 * STAGE 2:
 *
 * A fine grained atom block against interacting atoms neighbour list is constructed.
 *
 * The warp loops over atom blocks Y that were found to (possibly) interact with atom block X.  Each thread
 * in the warp loops over the 32 atoms in X and compares their positions to one particular atom from block Y.
 * If it finds one closer than the cutoff distance, the atom is added to the list of atoms interacting with block X.
 * This continues until the buffer fills up, at which point the results are written to global memory.
 *
 * [in] periodicBoxSize        - size of the rectangular periodic box
 * [in] invPeriodicBoxSize     - inverse of the periodic box
 * [in] blockCenter            - the center of each bounding box
 * [in] blockBoundingBox       - bounding box of each atom block
 * [out] interactionCount      - total number of tiles that have interactions
 * [out] interactingTiles      - set of blocks that have interactions
 * [out] interactingAtoms      - a list of atoms that interact with each atom block
 * [in] posq                   - x,y,z coordinates of each atom and charge q
 * [in] maxTiles               - maximum number of tiles to process, used for multi-GPUs
 * [in] startBlockIndex        - first block to process, used for multi-GPUs,
 * [in] numBlocks              - total number of atom blocks
 * [in] sortedBlocks           - a sorted list of atom blocks based on volume
 * [in] sortedBlockCenter      - sorted centers, duplicated for fast access to avoid indexing
 * [in] sortedBlockBoundingBox - sorted bounding boxes, duplicated for fast access
 * [in] exclusionIndices       - maps into exclusionRowIndices with the starting position for a given atom
 * [in] exclusionRowIndices    - stores the a continuous list of exclusions
 *           eg: block 0 is excluded from atom 3,5,6
 *               block 1 is excluded from atom 3,4
 *               block 2 is excluded from atom 1,3,5,6
 *              exclusionIndices[0][3][5][8]
 *           exclusionRowIndices[3][5][6][3][4][1][3][5][6]
 *                         index 0  1  2  3  4  5  6  7  8 
 * [out] oldPos                - stores the positions of the atoms in which this neighbourlist was built on
 *                             - this is used to decide when to rebuild a neighbourlist
 * [in] rebuildNeighbourList   - whether or not to execute this kernel
 *
 */
extern "C" __global__ __launch_bounds__(GROUP_SIZE,3) void findBlocksWithInteractions(real4 periodicBoxSize, real4 invPeriodicBoxSize, real4 periodicBoxVecX, real4 periodicBoxVecY, real4 periodicBoxVecZ,
        unsigned int* __restrict__ interactionCount, int* __restrict__ interactingTiles, unsigned int* __restrict__ interactingAtoms,
        int2* __restrict__ singlePairs, const real4* __restrict__ posq, unsigned int maxTiles, unsigned int maxSinglePairs, unsigned int startBlockIndex,
        unsigned int numBlocks, unsigned int* __restrict__ sortedBlocks, const real4* __restrict__ sortedBlockCenter, const half3* __restrict__ sortedBlockBoundingBox,
        const unsigned int* __restrict__ exclusionIndices, const unsigned int* __restrict__ exclusionRowIndices,
        real4* __restrict__ oldPositions, const int* __restrict__ rebuildNeighborList) {

    if (rebuildNeighborList[0] == 0)
        return; // The neighbor list doesn't need to be rebuilt.

    const int indexInWarp = threadIdx.x%32;
    const int warpStart = threadIdx.x-indexInWarp;
    const int totalWarps = blockDim.x*gridDim.x/32;
    const int warpIndex = (blockIdx.x*blockDim.x+threadIdx.x)/32;
    const int warpMask = (1<<indexInWarp)-1;
    __shared__ int workgroupBuffer[BUFFER_SIZE*(GROUP_SIZE/32)];
    __shared__ int workgroupFlagsBuffer[BUFFER_SIZE*(GROUP_SIZE/32)];
    __shared__ int warpExclusions[MAX_EXCLUSIONS*(GROUP_SIZE/32)];
    __shared__ real4 posBuffer[GROUP_SIZE];
    __shared__ volatile unsigned int workgroupTileIndex[GROUP_SIZE/32];
    __shared__ unsigned int workgroupPairStartIndex[GROUP_SIZE/32];
    int* sumBuffer = (int*) posBuffer; // Reuse the same buffer to save memory
    int* buffer = workgroupBuffer+BUFFER_SIZE*(warpStart/32);
    int* flagsBuffer = workgroupFlagsBuffer+BUFFER_SIZE*(warpStart/32);
    int* exclusionsForX = warpExclusions+MAX_EXCLUSIONS*(warpStart/32);
    volatile unsigned int& tileStartIndex = workgroupTileIndex[warpStart/32];
    volatile unsigned int& pairStartIndex = workgroupPairStartIndex[warpStart/32];

    // Loop over blocks.
    
    for (int block1 = startBlockIndex+warpIndex; block1 < startBlockIndex+numBlocks; block1 += totalWarps) {
        // Load data for this block.  Note that all threads in a warp are processing the same block.
        
        int x = sortedBlocks[block1] & BLOCK_INDEX_MASK;
        real4 blockCenterX = sortedBlockCenter[block1];
        real3 blockSizeX = sortedBlockBoundingBox[block1].toReal3();
        int neighborsInBuffer = 0;
        real4 pos1 = posq[x*TILE_SIZE+indexInWarp];
#ifdef USE_PERIODIC
        const bool singlePeriodicCopy = (0.5f*periodicBoxSize.x-blockSizeX.x >= PADDED_CUTOFF &&
                                         0.5f*periodicBoxSize.y-blockSizeX.y >= PADDED_CUTOFF &&
                                         0.5f*periodicBoxSize.z-blockSizeX.z >= PADDED_CUTOFF);
        if (singlePeriodicCopy) {
            // The box is small enough that we can just translate all the atoms into a single periodic
            // box, then skip having to apply periodic boundary conditions later.
            
            APPLY_PERIODIC_TO_POS_WITH_CENTER(pos1, blockCenterX)
        }
#endif
        pos1.w = 0.5f * (pos1.x * pos1.x + pos1.y * pos1.y + pos1.z * pos1.z);
        posBuffer[threadIdx.x] = pos1;

        // Load exclusion data for block x.
        
        const int exclusionStart = exclusionRowIndices[x];
        const int exclusionEnd = exclusionRowIndices[x+1];
        const int numExclusions = exclusionEnd-exclusionStart;
        #pragma unroll 4 // (MAX_EXCLUSIONS)
        for (int j = indexInWarp; j < numExclusions; j += 32)
            exclusionsForX[j] = exclusionIndices[exclusionStart+j];
        if (MAX_EXCLUSIONS > 32)
            __syncthreads();
        
        // Loop over atom blocks to search for neighbors.  The threads in a warp compare block1 against 32
        // other blocks in parallel.

        for (int block2Base = block1+1; block2Base < NUM_BLOCKS; block2Base += 32) {
            int block2 = block2Base+indexInWarp;
            bool includeBlock2 = (block2 < NUM_BLOCKS);
            bool forceInclude = false;
            if (includeBlock2) {
                real4 blockCenterY = sortedBlockCenter[block2];
                real3 blockSizeY = sortedBlockBoundingBox[block2].toReal3();
                real4 blockDelta = blockCenterX-blockCenterY;
#ifdef USE_PERIODIC
                APPLY_PERIODIC_TO_DELTA(blockDelta)
#endif
                includeBlock2 &= (blockDelta.x*blockDelta.x+blockDelta.y*blockDelta.y+blockDelta.z*blockDelta.z < (PADDED_CUTOFF+blockCenterX.w+blockCenterY.w)*(PADDED_CUTOFF+blockCenterX.w+blockCenterY.w));
                blockDelta.x = max(0.0f, fabs(blockDelta.x)-blockSizeX.x-blockSizeY.x);
                blockDelta.y = max(0.0f, fabs(blockDelta.y)-blockSizeX.y-blockSizeY.y);
                blockDelta.z = max(0.0f, fabs(blockDelta.z)-blockSizeX.z-blockSizeY.z);
                includeBlock2 &= (blockDelta.x*blockDelta.x+blockDelta.y*blockDelta.y+blockDelta.z*blockDelta.z < PADDED_CUTOFF_SQUARED);
#ifdef TRICLINIC
                // The calculation to find the nearest periodic copy is only guaranteed to work if the nearest copy is less than half a box width away.
                // If there's any possibility we might have missed it, do a detailed check.

                if (periodicBoxSize.z/2-blockSizeX.z-blockSizeY.z < PADDED_CUTOFF || periodicBoxSize.y/2-blockSizeX.y-blockSizeY.y < PADDED_CUTOFF)
                    includeBlock2 = forceInclude = true;
#endif
                if (includeBlock2) {
                    int y = sortedBlocks[block2] & BLOCK_INDEX_MASK;
                    #pragma unroll 4 // (MAX_EXCLUSIONS)
                    for (int k = 0; k < numExclusions; k++)
                        includeBlock2 &= (exclusionsForX[k] != y);
                }
            }
            
            // Loop over any blocks we identified as potentially containing neighbors.
            
            int includeBlockFlags = BALLOT(includeBlock2);
            int forceIncludeFlags = BALLOT(forceInclude);
            while (includeBlockFlags != 0) {
                int i = __ffs(includeBlockFlags)-1;
                includeBlockFlags &= includeBlockFlags-1;
                forceInclude = (forceIncludeFlags>>i) & 1;
                int y = sortedBlocks[block2Base+i] & BLOCK_INDEX_MASK;

                // Check each atom in block Y for interactions.

                int atom2 = y*TILE_SIZE+indexInWarp;
                real4 pos2 = posq[atom2];
#ifdef USE_PERIODIC
                if (singlePeriodicCopy) {
                    APPLY_PERIODIC_TO_POS_WITH_CENTER(pos2, blockCenterX)
                }
#endif
                pos2.w = 0.5f * (pos2.x * pos2.x + pos2.y * pos2.y + pos2.z * pos2.z);

                real4 blockCenterY = sortedBlockCenter[block2Base+i];
                real3 atomDelta = trimTo3(posBuffer[warpStart+indexInWarp])-trimTo3(blockCenterY);
#ifdef USE_PERIODIC
                APPLY_PERIODIC_TO_DELTA(atomDelta)
#endif
                int atomFlags = BALLOT(forceInclude || atomDelta.x*atomDelta.x+atomDelta.y*atomDelta.y+atomDelta.z*atomDelta.z < (PADDED_CUTOFF+blockCenterY.w)*(PADDED_CUTOFF+blockCenterY.w));
                int interacts = 0;
                if (atom2 < NUM_ATOMS && atomFlags != 0) {
#ifdef USE_PERIODIC
                    if (!singlePeriodicCopy) {
                        int first = __ffs(atomFlags)-1;
                        int last = 32-__clz(atomFlags);
                        for (int j = first; j < last; j++) {
                            real3 delta = trimTo3(pos2)-trimTo3(posBuffer[warpStart+j]);
                            APPLY_PERIODIC_TO_DELTA(delta)
                            interacts |= (delta.x*delta.x+delta.y*delta.y+delta.z*delta.z < PADDED_CUTOFF_SQUARED ? 1<<j : 0);
                        }
                    }
                    else {
#endif
                        #pragma unroll
                        for (int j = 0; j < 32; j++) {
                            real4 posj = posBuffer[warpStart+j];
                            real halfDist2 = posj.w + pos2.w - posj.x*pos2.x - posj.y*pos2.y - posj.z*pos2.z;
                            interacts |= (halfDist2 < 0.5f * PADDED_CUTOFF_SQUARED ? 1<<j : 0);
                        }
#ifdef USE_PERIODIC
                    }
#endif
                }
                
                // Add any interacting atoms to the buffer.
                
                int includeAtomFlags = BALLOT(interacts);
                if (interacts) {
                    int index = neighborsInBuffer+__popc(includeAtomFlags&warpMask);
                    buffer[index] = atom2;
                    flagsBuffer[index] = interacts;
                }
                neighborsInBuffer += __popc(includeAtomFlags);
                if (neighborsInBuffer > BUFFER_SIZE-TILE_SIZE) {
                    // Store the new tiles to memory.
                    
#if MAX_BITS_FOR_PAIRS > 0
                    neighborsInBuffer = saveSinglePairs(x, buffer, flagsBuffer, neighborsInBuffer, maxSinglePairs, &interactionCount[1], singlePairs, sumBuffer+warpStart, pairStartIndex);
#endif
                    unsigned int tilesToStore = neighborsInBuffer/TILE_SIZE;
                    if (tilesToStore > 0) {
                        if (indexInWarp == 0)
                            tileStartIndex = atomicAdd(&interactionCount[0], tilesToStore);
                        unsigned int newTileStartIndex = tileStartIndex;
                        if (newTileStartIndex+tilesToStore <= maxTiles) {
                            if (indexInWarp < tilesToStore)
                                interactingTiles[newTileStartIndex+indexInWarp] = x;
                            #pragma unroll 8 // (GROUP_SIZE / TILE_SIZE)
                            for (int j = 0; j < tilesToStore; j++)
                                interactingAtoms[(newTileStartIndex+j)*TILE_SIZE+indexInWarp] = buffer[indexInWarp+j*TILE_SIZE];
                        }
                        if (indexInWarp+TILE_SIZE*tilesToStore < BUFFER_SIZE)
                            buffer[indexInWarp] = buffer[indexInWarp+TILE_SIZE*tilesToStore];
                        neighborsInBuffer -= TILE_SIZE*tilesToStore;
                    }
                }
            }
        }
        
        // If we have a partially filled buffer,  store it to memory.
        
#if MAX_BITS_FOR_PAIRS > 0
        if (neighborsInBuffer > 32)
            neighborsInBuffer = saveSinglePairs(x, buffer, flagsBuffer, neighborsInBuffer, maxSinglePairs, &interactionCount[1], singlePairs, sumBuffer+warpStart, pairStartIndex);
#endif
        if (neighborsInBuffer > 0) {
            unsigned int tilesToStore = (neighborsInBuffer+TILE_SIZE-1)/TILE_SIZE;
            if (indexInWarp == 0)
                tileStartIndex = atomicAdd(&interactionCount[0], tilesToStore);
            unsigned int newTileStartIndex = tileStartIndex;
            if (newTileStartIndex+tilesToStore <= maxTiles) {
                if (indexInWarp < tilesToStore)
                    interactingTiles[newTileStartIndex+indexInWarp] = x;
                #pragma unroll 8 // (GROUP_SIZE / TILE_SIZE)
                for (int j = 0; j < tilesToStore; j++)
                    interactingAtoms[(newTileStartIndex+j)*TILE_SIZE+indexInWarp] = (indexInWarp+j*TILE_SIZE < neighborsInBuffer ? buffer[indexInWarp+j*TILE_SIZE] : NUM_ATOMS);
            }
        }
    }
    
    // Record the positions the neighbor list is based on.
    
    for (int i = threadIdx.x+blockIdx.x*blockDim.x; i < NUM_ATOMS; i += blockDim.x*gridDim.x)
        oldPositions[i] = posq[i];
}
#endif
