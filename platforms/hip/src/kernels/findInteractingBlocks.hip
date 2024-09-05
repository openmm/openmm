#define BUFFER_SIZE 256

#if defined(AMD_RDNA)

typedef unsigned int warpflags;

__device__ inline int warpPopc(warpflags x) {
    return __popc(x);
}

#else

typedef unsigned long long warpflags;

__device__ inline int warpPopc(warpflags x) {
    return __popcll(x);
}

#endif

struct alignas(sizeof(__half) * 4) BoundingBox {
    __device__ BoundingBox(real3 f) {
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

/**
 * Find a bounding box for the atoms in each block.
 */
extern "C" __global__ void findBlockBounds(int numAtoms, real4 periodicBoxSize, real4 invPeriodicBoxSize, real4 periodicBoxVecX, real4 periodicBoxVecY, real4 periodicBoxVecZ,
        const real4* __restrict__ posq, real4* __restrict__ blockCenter, real4* __restrict__ blockBoundingBox, int* __restrict__ rebuildNeighborList,
        real2* __restrict__ blockSizeRange) {
    const int indexInTile = threadIdx.x%TILE_SIZE;
    const int index = warpSize == TILE_SIZE ? blockIdx.x : (blockIdx.x*(warpSize/TILE_SIZE) + threadIdx.x/TILE_SIZE);
    const int base = index * TILE_SIZE;
    if (base >= numAtoms)
        return;

    real4 tPos = posq[base + indexInTile < numAtoms ? base + indexInTile : 0];
#ifdef USE_PERIODIC
    real4 pos;
    pos.x = SHFL(tPos.x, 0);
    pos.y = SHFL(tPos.y, 0);
    pos.z = SHFL(tPos.z, 0);
    APPLY_PERIODIC_TO_POS(pos)

    real4 minPos = pos;
    real4 maxPos = pos;

    for (int i = 1; i < TILE_SIZE; i++) {
        pos.x = SHFL(tPos.x, i);
        pos.y = SHFL(tPos.y, i);
        pos.z = SHFL(tPos.z, i);
        real4 center = 0.5f*(maxPos+minPos);
        APPLY_PERIODIC_TO_POS_WITH_CENTER(pos, center)
        minPos = make_real4(min(minPos.x,pos.x), min(minPos.y,pos.y), min(minPos.z,pos.z), 0);
        maxPos = make_real4(max(maxPos.x,pos.x), max(maxPos.y,pos.y), max(maxPos.z,pos.z), 0);
    }
#else
    real4 minPos = tPos;
    real4 maxPos = tPos;

    for (int i = TILE_SIZE >> 1; i > 0; i >>= 1) {
        real4 tpos1, tpos2;
        tpos1.x = __shfl_down(minPos.x, i, TILE_SIZE);
        tpos1.y = __shfl_down(minPos.y, i, TILE_SIZE);
        tpos1.z = __shfl_down(minPos.z, i, TILE_SIZE);
        tpos2.x = __shfl_down(maxPos.x, i, TILE_SIZE);
        tpos2.y = __shfl_down(maxPos.y, i, TILE_SIZE);
        tpos2.z = __shfl_down(maxPos.z, i, TILE_SIZE);

        minPos.x = min(minPos.x, tpos1.x);
        minPos.y = min(minPos.y, tpos1.y);
        minPos.z = min(minPos.z, tpos1.z);
        maxPos.x = max(maxPos.x, tpos2.x);
        maxPos.y = max(maxPos.y, tpos2.y);
        maxPos.z = max(maxPos.z, tpos2.z);
    }

    minPos.x = SHFL(minPos.x, 0);
    minPos.y = SHFL(minPos.y, 0);
    minPos.z = SHFL(minPos.z, 0);
    maxPos.x = SHFL(maxPos.x, 0);
    maxPos.y = SHFL(maxPos.y, 0);
    maxPos.z = SHFL(maxPos.z, 0);
#endif

    real4 blockSize = 0.5f*(maxPos-minPos);
    real4 center = 0.5f*(maxPos+minPos);
    center.w = 0;
    real4 delta = tPos - center;
#ifdef USE_PERIODIC
    APPLY_PERIODIC_TO_DELTA(delta)
#endif
    real tdelta = delta.x*delta.x+delta.y*delta.y+delta.z*delta.z;
    real tcenter = max(center.w, tdelta);
    for (int i = TILE_SIZE >> 1; i > 0; i >>= 1) {
        real t = __shfl_down(tcenter, i, TILE_SIZE);
        tcenter = max(tcenter, t);
    }

    if (indexInTile == 0) {
        center.w = SQRT(tcenter);
        blockBoundingBox[index] = blockSize;
        blockCenter[index] = center;

        // Record the range of sizes.
        real totalSize = blockSize.x+blockSize.y+blockSize.z;
        atomicMin(&blockSizeRange->x, totalSize);
        atomicMax(&blockSizeRange->y, totalSize);
    }
    if (blockIdx.x == 0 && threadIdx.x == 0)
        rebuildNeighborList[0] = 0;
}

extern "C" __global__ void computeSortKeys(const real4* __restrict__ blockBoundingBox, unsigned int* __restrict__ sortedBlocks, const real2* __restrict__ blockSizeRange/*, int numSizes*/) {
    // Sort keys store the bin in the high order part and the block in the low
    // order part.

    real2 sizeRange = make_real2(LOG(blockSizeRange->x), LOG(blockSizeRange->y));
    int numSizeBins = 20;
    real scale = numSizeBins/(sizeRange.y-sizeRange.x);
    unsigned int i = threadIdx.x+blockIdx.x*blockDim.x;
    if (i < NUM_BLOCKS) {
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
        const real4* __restrict__ blockBoundingBox, real4* __restrict__ sortedBlockCenter, BoundingBox* __restrict__ sortedBlockBoundingBox,
#ifdef USE_LARGE_BLOCKS
        real4* __restrict__ largeBlockCenter, BoundingBox* __restrict__ largeBlockBoundingBox, real4 periodicBoxSize,
        real4 invPeriodicBoxSize, real4 periodicBoxVecX, real4 periodicBoxVecY, real4 periodicBoxVecZ,
#endif
        const real4* __restrict__ posq, const real4* __restrict__ oldPositions,
        unsigned int* __restrict__ interactionCount, int* __restrict__ rebuildNeighborList, bool forceRebuild,
        real2* __restrict__ blockSizeRange) {
    int i = threadIdx.x+blockIdx.x*blockDim.x;
    if (i < NUM_BLOCKS) {
        unsigned int index = sortedBlocks[i] & BLOCK_INDEX_MASK;
        sortedBlockCenter[i] = blockCenter[index];
        sortedBlockBoundingBox[i] = BoundingBox(trimTo3(blockBoundingBox[index]));

#ifdef USE_LARGE_BLOCKS
        // Compute the sizes of large blocks (composed of warpSize regular blocks) starting from each block.

        real4 minPos = blockCenter[index]-blockBoundingBox[index];
        real4 maxPos = blockCenter[index]+blockBoundingBox[index];
        int last = min(i+warpSize, NUM_BLOCKS);
        for (int j = i+1; j < last; j++) {
            unsigned int index2 = sortedBlocks[j] & BLOCK_INDEX_MASK;
            real4 blockPos = blockCenter[index2];
            real4 width = blockBoundingBox[index2];
#ifdef USE_PERIODIC
            real4 center = 0.5f*(maxPos+minPos);
            APPLY_PERIODIC_TO_POS_WITH_CENTER(blockPos, center)
#endif
            minPos = make_real4(min(minPos.x, blockPos.x-width.x), min(minPos.y, blockPos.y-width.y), min(minPos.z, blockPos.z-width.z), 0);
            maxPos = make_real4(max(maxPos.x, blockPos.x+width.x), max(maxPos.y, blockPos.y+width.y), max(maxPos.z, blockPos.z+width.z), 0);
        }
        largeBlockCenter[i] = 0.5f*(maxPos+minPos);
        largeBlockBoundingBox[i] = BoundingBox(trimTo3(0.5f*(maxPos-minPos)));
#endif
    }

    // Initialize min and max for the atomic reduction in findBlockBounds (for the next step)
    if (i == 0) {
        blockSizeRange->x = 1e38;
        blockSizeRange->y = 0;
    }

    // Also check whether any atom has moved enough so that we really need to rebuild the neighbor list.

    bool rebuild = forceRebuild;
    if (i < NUM_ATOMS) {
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

#if MAX_BITS_FOR_PAIRS > 0

__device__ inline
void collectInteractions(unsigned int& interacts, float d, int bit) {
    interacts |= (__float_as_uint(d) >> 31) << bit;
}

__device__ inline
void collectInteractions(unsigned int& interacts, double d, int bit) {
    interacts |= ((unsigned int)__double2hiint(d) >> 31) << bit;
}

#else

// Simplified version that does not collect individual bits (they are not needed without single pairs),
// only sets the flag that there are any interactions.

__device__ inline
void collectInteractions(unsigned int& interacts, float d, int) {
    interacts |= __float_as_uint(d) & (1 << 31);
}

__device__ inline
void collectInteractions(unsigned int& interacts, double d, int) {
    interacts |= (unsigned int)__double2hiint(d) & (1 << 31);
}

#endif

#if !defined(USE_DOUBLE_PRECISION) && __has_builtin(__builtin_amdgcn_mfma_f32_4x4x1f32)

#define USE_MFMA

using vfloat = __attribute__((__vector_size__(4 * sizeof(float)))) float;

template<int BlockId>
inline __device__
void mfma4x4(const float4& pos1, const float4& pos2, const vfloat& c, unsigned int& interacts) {
    vfloat d;
    d = __builtin_amdgcn_mfma_f32_4x4x1f32(pos1.x, -pos2.x, c, 4, BlockId, 0);
    d = __builtin_amdgcn_mfma_f32_4x4x1f32(pos1.y, -pos2.y, d, 4, BlockId, 0);
    d = __builtin_amdgcn_mfma_f32_4x4x1f32(pos1.z, -pos2.z, d, 4, BlockId, 0);
    d = __builtin_amdgcn_mfma_f32_4x4x1f32(pos1.w,    1.0f, d, 4, BlockId, 0);
    #pragma unroll
    for (int i = 0; i < 4; i++) {
        collectInteractions(interacts, d[i], BlockId * 4 + i);
    }
}

#endif

/**
 * Compare the bounding boxes for each pair of atom blocks (comprised of TILE_SIZE atoms each), forming a tile. If the two
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
 * in the warp loops over the TILE_SIZE atoms in X and compares their positions to one particular atom from block Y.
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
extern "C" __global__ __launch_bounds__(GROUP_SIZE) void findBlocksWithInteractions(real4 periodicBoxSize, real4 invPeriodicBoxSize, real4 periodicBoxVecX, real4 periodicBoxVecY, real4 periodicBoxVecZ,
        unsigned int* __restrict__ interactionCount, int* __restrict__ interactingTiles, unsigned int* __restrict__ interactingAtoms,
        int2* __restrict__ singlePairs, const real4* __restrict__ posq, unsigned int maxTiles, unsigned int maxSinglePairs,
        unsigned int startBlockIndex, unsigned int numBlocks, const unsigned int* __restrict__ sortedBlocks, const real4* __restrict__ sortedBlockCenter,
        const BoundingBox* __restrict__ sortedBlockBoundingBox,
#ifdef USE_LARGE_BLOCKS
        const real4* __restrict__ largeBlockCenter, const BoundingBox* __restrict__ largeBlockBoundingBox,
#endif
        const unsigned int* __restrict__ exclusionIndices, const unsigned int* __restrict__ exclusionRowIndices,
        real4* __restrict__ oldPositions, const int* __restrict__ rebuildNeighborList) {

    if (rebuildNeighborList[0] == 0)
        return; // The neighbor list doesn't need to be rebuilt.

    constexpr int tilesPerWarp = warpSize/TILE_SIZE;
    constexpr int warpsPerBlock = GROUP_SIZE/warpSize;
    const int indexInWarp = threadIdx.x%warpSize;
    const int indexInTile = threadIdx.x%TILE_SIZE;
    const int tileInWarp = tilesPerWarp == 1 ? 0 : (threadIdx.x/TILE_SIZE)%tilesPerWarp;
    const int warpInBlock = warpsPerBlock == 1 ? 0 : threadIdx.x/warpSize;
    const int warpIndex = blockIdx.x*warpsPerBlock + (warpsPerBlock == 1 ? 0 : warpInBlock);
    const warpflags warpMask = (static_cast<warpflags>(1)<<indexInWarp)-1;

    __shared__ int workgroupBuffer[BUFFER_SIZE*warpsPerBlock];
    __shared__ real4 workgroupPosBuffer[TILE_SIZE*warpsPerBlock];
    __shared__ int workgroupExclusions[MAX_EXCLUSIONS*warpsPerBlock];
    __shared__ int workgroupBlock2Buffer[(warpSize+1)*warpsPerBlock];

    int* buffer = workgroupBuffer+BUFFER_SIZE*warpInBlock;
    real4* posBuffer = workgroupPosBuffer+TILE_SIZE*warpInBlock;
    int* exclusionsForX = workgroupExclusions+MAX_EXCLUSIONS*warpInBlock;
    int* block2Buffer = workgroupBlock2Buffer+(warpSize+1)*warpInBlock;

    // Loop over blocks.

    int block1 = startBlockIndex+warpIndex/NUM_TILES_IN_BATCH;
    if (block1 < startBlockIndex+numBlocks) {
        // Load data for this block.  Note that all threads in a warp are processing the same block.

        int x = sortedBlocks[block1] & BLOCK_INDEX_MASK;
        real4 blockCenterX = sortedBlockCenter[block1];
        real3 blockSizeX = sortedBlockBoundingBox[block1].toReal3();
        int neighborsInBuffer = 0;
        real4 pos1 = posq[x*TILE_SIZE+indexInTile];
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
        if (tileInWarp == 0) {
            posBuffer[indexInTile] = pos1;
        }

        // Load exclusion data for block x.

        const int exclusionStart = exclusionRowIndices[x];
        const int exclusionEnd = exclusionRowIndices[x+1];
        const int numExclusions = exclusionEnd-exclusionStart;
        #pragma unroll 1
        for (int j = indexInWarp; j < numExclusions; j += warpSize)
            exclusionsForX[j] = exclusionIndices[exclusionStart+j];

        // Loop over atom blocks to search for neighbors.  The threads in a warp compare block1 against warpSize
        // other blocks in parallel.
        // For small systems multiple warps (NUM_TILES_IN_BATCH = 4, 2...) process one block1 reducing the overall
        // duration of the kernel because first blocks block1 have to process more block2 blocks so most of compute
        // units are idle at the end of the kernel (the kernel works on the upper triangle of
        // the NUM_BLOCKS x NUM_BLOCKS matrix).

#ifdef USE_LARGE_BLOCKS
        warpflags largeBlockFlags = 0;
        int loadedLargeBlocks = 0;
#endif
        int block2Count = 0;
        // Load blocks from addresses aligned by warpSize for faster loading from sortedBlockCenter and sortedBlockBoundingBox.
        for (int block2Base = ((block1+1)/warpSize + warpIndex%NUM_TILES_IN_BATCH)*warpSize; block2Base < NUM_BLOCKS; block2Base += warpSize*NUM_TILES_IN_BATCH) {
            // Last iteration cannot be skipped (on CDNA where tilesPerWarp == 2)
            const bool lastIteration = block2Base + warpSize*NUM_TILES_IN_BATCH >= NUM_BLOCKS;
#ifdef USE_LARGE_BLOCKS
            if (loadedLargeBlocks == 0) {
                // Check the next set of large blocks.

                int largeBlockIndex = block2Base + warpSize*NUM_TILES_IN_BATCH*indexInWarp;
                bool includeLargeBlock = false;
                if (largeBlockIndex < NUM_BLOCKS) {
                    real4 largeCenter = largeBlockCenter[largeBlockIndex];
                    real3 largeSize = largeBlockBoundingBox[largeBlockIndex].toReal3();
                    real4 blockDelta = blockCenterX-largeCenter;
#ifdef USE_PERIODIC
                    APPLY_PERIODIC_TO_DELTA(blockDelta)
#endif
                    blockDelta.x = max(0.0f, fabs(blockDelta.x)-blockSizeX.x-largeSize.x);
                    blockDelta.y = max(0.0f, fabs(blockDelta.y)-blockSizeX.y-largeSize.y);
                    blockDelta.z = max(0.0f, fabs(blockDelta.z)-blockSizeX.z-largeSize.z);
                    includeLargeBlock = (blockDelta.x*blockDelta.x+blockDelta.y*blockDelta.y+blockDelta.z*blockDelta.z < PADDED_CUTOFF_SQUARED);
#ifdef TRICLINIC
                    // The calculation to find the nearest periodic copy is only guaranteed to work if the nearest copy is less than half a box width away.
                    // If there's any possibility we might have missed it, do a detailed check.

                    if (periodicBoxSize.z/2-blockSizeX.z-largeSize.z < PADDED_CUTOFF || periodicBoxSize.y/2-blockSizeX.y-largeSize.y < PADDED_CUTOFF)
                        includeLargeBlock = true;
#endif
                }
                largeBlockFlags = __ballot(includeLargeBlock);
                loadedLargeBlocks = warpSize;
            }
            loadedLargeBlocks--;
            if ((largeBlockFlags&1) == 0 && !lastIteration) {
                // None of the next warpSize blocks interact with block 1.

                largeBlockFlags >>= 1;
                continue;
            }
            largeBlockFlags >>= 1;
#endif
            int block2 = block2Base+indexInWarp;
            bool includeBlock2 = (block1 < block2 && block2 < NUM_BLOCKS);
            block2 = includeBlock2 ? block2 : block1;
            bool forceInclude = false;
            real4 blockCenterY = sortedBlockCenter[block2];
            real4 blockDelta = blockCenterX-blockCenterY;
#ifdef USE_PERIODIC
            APPLY_PERIODIC_TO_DELTA(blockDelta)
#endif
            includeBlock2 &= (blockDelta.x*blockDelta.x+blockDelta.y*blockDelta.y+blockDelta.z*blockDelta.z < (PADDED_CUTOFF+blockCenterX.w+blockCenterY.w)*(PADDED_CUTOFF+blockCenterX.w+blockCenterY.w));
#ifndef TRICLINIC
            if (!lastIteration && __ballot(includeBlock2) == 0)
                continue;
#endif
            real3 blockSizeY = sortedBlockBoundingBox[block2].toReal3();
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

            // Collect any blocks we identified as potentially containing neighbors.

            warpflags includeBlockFlags = __ballot(includeBlock2);
            if (includeBlock2) {
                int index = block2Count + warpPopc(includeBlockFlags&warpMask);
                block2Buffer[index] = (block2 << 1) | (forceInclude ? 1 : 0);
            }
            block2Count += warpPopc(includeBlockFlags);

            // Loop over the collected candidates (each warp processes 2 blocks on CDNA or 1 block on RDNA).
            // Process even number of blocks on CDNA so both half-warps have work to do (except for
            // the last iteration of the for-block2Base loop when the left-over must be processed).

            const int block2ToProcess = lastIteration ? block2Count : block2Count/tilesPerWarp*tilesPerWarp;
            for (int block2Index = 0; block2Index < block2ToProcess; block2Index += tilesPerWarp) {
                bool includeBlock2 = block2Index + tileInWarp < block2Count;
                const int b = block2Buffer[min(block2Index + tileInWarp, block2Count - 1)];
                const bool forceInclude = b & 1;
                const int block2 = b >> 1;
                int y = sortedBlocks[block2] & BLOCK_INDEX_MASK;

                #pragma unroll 1
                for (int k = indexInTile; k < numExclusions; k += TILE_SIZE)
                    includeBlock2 &= (exclusionsForX[k] != y);
                includeBlock2 = BALLOT(!includeBlock2) == 0;

                // Check each atom in block Y for interactions.

                int atom2 = y*TILE_SIZE+indexInTile;
                real4 pos2 = posq[atom2];
#ifdef USE_PERIODIC
                if (singlePeriodicCopy) {
                    APPLY_PERIODIC_TO_POS_WITH_CENTER(pos2, blockCenterX)
                }
#endif
                pos2.w = 0.5f * (pos2.x * pos2.x + pos2.y * pos2.y + pos2.z * pos2.z);
                real4 blockCenterY = sortedBlockCenter[block2];
                real3 atomDelta = trimTo3(pos1)-trimTo3(blockCenterY);
#ifdef USE_PERIODIC
                APPLY_PERIODIC_TO_DELTA(atomDelta)
#endif
                tileflags atomFlags = BALLOT(forceInclude || atomDelta.x*atomDelta.x+atomDelta.y*atomDelta.y+atomDelta.z*atomDelta.z < (PADDED_CUTOFF+blockCenterY.w)*(PADDED_CUTOFF+blockCenterY.w));
                tileflags interacts = 0;
                // The condition `posj.w + pos2.w - posj.x*pos2.x - posj.y*pos2.y - posj.z*pos2.z < 0.5f * PADDED_CUTOFF_SQUARED` is expressed as
                // `posj.x*pos2.x - posj.y*pos2.y - posj.z*pos2.z - posj.w - 0.5f * PADDED_CUTOFF_SQUARED - pos2.w` and computed using fma
                // (it saves 1 instruction).
                // Sign bit is used directly instead of `halfDist2 < 0.5f * PADDED_CUTOFF_SQUARED ? 1<<j : 0`.
#ifdef USE_PERIODIC
                if (!singlePeriodicCopy) {
                    while (atomFlags) {
                        int j = __ffs(atomFlags)-1;
                        atomFlags = atomFlags ^ (static_cast<tileflags>(1) << j);
                        real3 delta = trimTo3(pos2)-trimTo3(posBuffer[j]);
                        APPLY_PERIODIC_TO_DELTA(delta)
                        real d = delta.x*delta.x+delta.y*delta.y+delta.z*delta.z - PADDED_CUTOFF_SQUARED;
                        collectInteractions(interacts, d, j);
                    }
                }
                else {
#endif
                    const real lim = 0.5f * PADDED_CUTOFF_SQUARED - pos2.w;
#if defined(USE_MFMA)
                    const vfloat c = { -lim, -lim, -lim, -lim };
                    mfma4x4<0>(pos1, pos2, c, interacts);
                    mfma4x4<1>(pos1, pos2, c, interacts);
                    mfma4x4<2>(pos1, pos2, c, interacts);
                    mfma4x4<3>(pos1, pos2, c, interacts);
                    mfma4x4<4>(pos1, pos2, c, interacts);
                    mfma4x4<5>(pos1, pos2, c, interacts);
                    mfma4x4<6>(pos1, pos2, c, interacts);
                    mfma4x4<7>(pos1, pos2, c, interacts);
#else
                    while (atomFlags) {
                        int j = __ffs(atomFlags)-1;
                        atomFlags = atomFlags ^ (static_cast<tileflags>(1) << j);
                        real4 posj = posBuffer[j];
                        real d = fma(-posj.x, pos2.x, fma(-posj.y, pos2.y, fma(-posj.z, pos2.z, posj.w - lim)));
                        collectInteractions(interacts, d, j);
                    }
#endif
#ifdef USE_PERIODIC
                }
#endif
                if (atom2 >= NUM_ATOMS || !includeBlock2) {
                    interacts = 0;
                }

#if MAX_BITS_FOR_PAIRS > 0
                const unsigned int interactCount = __popc(interacts);

                // Record interactions that should be computed as single pairs rather than in blocks.
                const bool storeAsSinglePair = interactCount > 0 && interactCount <= MAX_BITS_FOR_PAIRS;
                if (__ballot(storeAsSinglePair)) {
                    unsigned int sum = 0;
                    unsigned int prevSum = 0;
                    for (int i = 1; i <= MAX_BITS_FOR_PAIRS; i++) {
                        warpflags b = __ballot(interactCount == i);
                        sum += warpPopc(b) * i;
                        prevSum += warpPopc(b&warpMask) * i;
                    }
                    unsigned int pairStartIndex = 0;
                    if (indexInWarp == 0)
                        pairStartIndex = atomicAdd(&interactionCount[1], sum);
                    pairStartIndex = __shfl(pairStartIndex, 0);
                    unsigned int pairIndex = pairStartIndex + prevSum;
                    if (storeAsSinglePair && pairIndex+interactCount <= maxSinglePairs) {
                        while (interacts != 0) {
                            int j = __ffs(interacts)-1;
                            singlePairs[pairIndex] = make_int2(atom2, x*TILE_SIZE+j);
                            interacts = interacts ^ (static_cast<tileflags>(1) << j);
                            pairIndex++;
                        }
                    }
                }
#else
                const unsigned int interactCount = interacts;
#endif

                // Add any interacting atoms to the buffer.

                warpflags includeAtomFlags = __ballot(interactCount > MAX_BITS_FOR_PAIRS);
                if (interactCount > MAX_BITS_FOR_PAIRS) {
                    int index = neighborsInBuffer+warpPopc(includeAtomFlags&warpMask);
                    buffer[index] = atom2;
                }
                neighborsInBuffer += warpPopc(includeAtomFlags);
                if (neighborsInBuffer > BUFFER_SIZE-warpSize) {
                    // Store the new tiles to memory.

                    unsigned int tilesToStore = neighborsInBuffer/warpSize*tilesPerWarp;
                    unsigned int tileStartIndex = 0;
                    if (indexInWarp == 0)
                        tileStartIndex = atomicAdd(&interactionCount[0], tilesToStore);
                    unsigned int newTileStartIndex = __shfl(tileStartIndex, 0);
                    if (newTileStartIndex+tilesToStore <= maxTiles) {
                        if (indexInWarp < tilesToStore)
                            interactingTiles[newTileStartIndex+indexInWarp] = x;
                        for (int j = 0; j < tilesToStore/tilesPerWarp; j++)
                            interactingAtoms[newTileStartIndex*TILE_SIZE+j*warpSize+indexInWarp] = buffer[j*warpSize+indexInWarp];
                    }
                    if (indexInWarp+TILE_SIZE*tilesToStore < BUFFER_SIZE)
                        buffer[indexInWarp] = buffer[indexInWarp+TILE_SIZE*tilesToStore];
                    neighborsInBuffer -= TILE_SIZE*tilesToStore;
                }
            }

            // Move not processed blocks to the head of block2Buffer.
            if (indexInWarp < block2Count - block2ToProcess)
                block2Buffer[indexInWarp] = block2Buffer[block2ToProcess + indexInWarp];
            block2Count = block2Count - block2ToProcess;
        }

        // If we have a partially filled buffer,  store it to memory.

        if (neighborsInBuffer > 0) {
            unsigned int tilesToStore = (neighborsInBuffer+TILE_SIZE-1)/TILE_SIZE;
            unsigned int tileStartIndex = 0;
            if (indexInWarp == 0)
                tileStartIndex = atomicAdd(&interactionCount[0], tilesToStore);
            unsigned int newTileStartIndex = __shfl(tileStartIndex, 0);
            if (newTileStartIndex+tilesToStore <= maxTiles) {
                if (indexInWarp < tilesToStore)
                    interactingTiles[newTileStartIndex+indexInWarp] = x;
                for (int j = 0; j <= tilesToStore/tilesPerWarp; j++) {
                    if (j*warpSize+indexInWarp < tilesToStore*TILE_SIZE)
                        interactingAtoms[newTileStartIndex*TILE_SIZE+j*warpSize+indexInWarp] = (j*warpSize+indexInWarp < neighborsInBuffer ? buffer[j*warpSize+indexInWarp] : PADDED_NUM_ATOMS);
                }
            }
        }
    }

    // Record the positions the neighbor list is based on.

    int i = threadIdx.x+blockIdx.x*GROUP_SIZE;
    if (i < NUM_ATOMS) {
        oldPositions[i] = posq[i];
    }
}

extern "C" __global__ void copyInteractionCounts(const unsigned int* __restrict__ interactionCount,
        unsigned int* __restrict__ pinnedInteractionCount) {
    pinnedInteractionCount[0] = interactionCount[0];
    pinnedInteractionCount[1] = interactionCount[1];
}
