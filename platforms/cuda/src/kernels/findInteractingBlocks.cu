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

__device__ int saveSinglePairs(int x, int* atoms, int* flags, int length, unsigned int maxSinglePairs, unsigned int* singlePairCount, int2* singlePairs) {
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

#ifdef USE_LARGE_BLOCKS

#define WARP_MAX3(v) \
    _Pragma("unroll") \
    for (int offset = 16; offset > 0; offset >>= 1) { \
        (v).x = max((v).x, __shfl_xor_sync(0xffffffff, (v).x, offset)); \
        (v).y = max((v).y, __shfl_xor_sync(0xffffffff, (v).y, offset)); \
        (v).z = max((v).z, __shfl_xor_sync(0xffffffff, (v).z, offset)); }

#define WARP_MINMAX3(mn, mx) \
    _Pragma("unroll") \
    for (int offset = 16; offset > 0; offset >>= 1) { \
        (mn).x = min((mn).x, __shfl_xor_sync(0xffffffff, (mn).x, offset)); \
        (mn).y = min((mn).y, __shfl_xor_sync(0xffffffff, (mn).y, offset)); \
        (mn).z = min((mn).z, __shfl_xor_sync(0xffffffff, (mn).z, offset)); \
        (mx).x = max((mx).x, __shfl_xor_sync(0xffffffff, (mx).x, offset)); \
        (mx).y = max((mx).y, __shfl_xor_sync(0xffffffff, (mx).y, offset)); \
        (mx).z = max((mx).z, __shfl_xor_sync(0xffffffff, (mx).z, offset)); }

#ifdef USE_DOUBLE_PRECISION
#define PACK_UINT_TO_REAL(u) __longlong_as_double((long long) (u))
#define UNPACK_REAL_TO_UINT(r) ((unsigned int) __double_as_longlong(r))
#else
#define PACK_UINT_TO_REAL(u) __uint_as_float(u)
#define UNPACK_REAL_TO_UINT(r) __float_as_uint(r)
#endif

__device__ __forceinline__ real4 loadSortedBlockReal4(const real4* __restrict__ ptr, int idx) {
#ifdef USE_DOUBLE_PRECISION
    return ptr[idx];
#else
    return __ldg(&ptr[idx]);
#endif
}

/**
 * Check atom displacements since the last neighbor-list rebuild and refresh block extents on non-rebuild steps.
 */
extern "C" __global__ void findBlockBounds(int numAtoms, real4 periodicBoxSize, real4 invPeriodicBoxSize, real4 periodicBoxVecX, real4 periodicBoxVecY, real4 periodicBoxVecZ,
        const real4* __restrict__ posq, real4* __restrict__ blockCenter, real4* __restrict__ blockBoundingBox, int* __restrict__ rebuildNeighborList,
        const real4* __restrict__ oldPositions, unsigned int* __restrict__ interactionCount, bool forceRebuild, unsigned int rebuildStep) {
    if (blockIdx.x == 0 && threadIdx.x == 0) {
        if (forceRebuild) {
            rebuildNeighborList[1] = (int) rebuildStep;
            rebuildNeighborList[2] = 0;
            interactionCount[0] = interactionCount[1] = 0;
        }
    }
    if (forceRebuild)
        return;
    const int lane = threadIdx.x & 31;
    for (int base = ((blockIdx.x * blockDim.x + threadIdx.x) & ~31); base < NUM_BLOCKS * TILE_SIZE; base += blockDim.x * gridDim.x) {
        int atom = base + lane;
        int index = base >> 5;
        real4 oldCenter = blockCenter[index];
        real3 ext = make_real3(0, 0, 0);
        if (atom < numAtoms) {
            real4 pos = posq[atom];
            real4 deltaMove = oldPositions[atom] - pos;
            if (deltaMove.x*deltaMove.x + deltaMove.y*deltaMove.y + deltaMove.z*deltaMove.z > 0.25f*PADDING*PADDING) {
                rebuildNeighborList[1] = (int) rebuildStep;
                rebuildNeighborList[2] = 0;
                interactionCount[0] = interactionCount[1] = 0;
            }
#ifdef USE_PERIODIC
            APPLY_PERIODIC_TO_POS_WITH_CENTER(pos, oldCenter)
#endif
            ext = make_real3(fabs(pos.x - oldCenter.x), fabs(pos.y - oldCenter.y), fabs(pos.z - oldCenter.z));
        }
        WARP_MAX3(ext)
        if (lane == 0)
            blockBoundingBox[index] = make_real4(ext.x, ext.y, ext.z, 0);
    }
}

/**
 * Compute block bounding boxes and 32-block spatial superblocks on neighbor-list rebuild steps.
 */
extern "C" __global__ void sortBoxData(real4* __restrict__ sortedBlocks, real4* __restrict__ blockCenter,
        real4* __restrict__ blockBoundingBox, real4* __restrict__ sortedBlockCenter,
        half3* __restrict__ sortedBlockBoundingBox, const real4* __restrict__ posq, real4* __restrict__ oldPositions,
        unsigned int* __restrict__ interactionCount, int* __restrict__ rebuildNeighborList, unsigned int rebuildStep,
        int numAtoms, real4 periodicBoxSize, real4 invPeriodicBoxSize, real4 periodicBoxVecX, real4 periodicBoxVecY, real4 periodicBoxVecZ) {
    if ((unsigned int) rebuildNeighborList[1] != rebuildStep) {
        if (blockIdx.x == 0 && threadIdx.x == 0)
            rebuildNeighborList[0] = 0;
        return;
    }
    if (blockIdx.x == 0 && threadIdx.x == 0)
        rebuildNeighborList[0] = 1;
    const int lane = threadIdx.x & 31;
    const int warpId = threadIdx.x >> 5;
    const int numSuperblocks = (NUM_BLOCKS + 31) / 32;
    const real compactLimit = 1.65f * PADDED_CUTOFF;
    __shared__ real4 shCenter[32];
    __shared__ real4 shSize[32];
    for (int sb = blockIdx.x; sb < numSuperblocks; sb += gridDim.x) {
        for (int pass = 0; pass < 4; pass++) {
            int sub = pass * 8 + warpId;
            int index = sb * 32 + sub;
            if (index < NUM_BLOCKS) {
                int atom = index * TILE_SIZE + lane;
                bool valid = (atom < numAtoms);
                real4 pos = (valid ? posq[atom] : make_real4(0, 0, 0, 0));
                if (valid)
                    oldPositions[atom] = pos;
#ifdef USE_PERIODIC
                if (lane == 0) {
                    APPLY_PERIODIC_TO_POS(pos)
                }
#endif
                real4 firstPos = make_real4(__shfl_sync(0xffffffff, pos.x, 0), __shfl_sync(0xffffffff, pos.y, 0), __shfl_sync(0xffffffff, pos.z, 0), 0);
#ifdef USE_PERIODIC
                if (valid && lane > 0) {
                    APPLY_PERIODIC_TO_POS_WITH_CENTER(pos, firstPos)
                }
#endif
                if (!valid)
                    pos = firstPos;
                real4 minPos = pos, maxPos = pos;
                WARP_MINMAX3(minPos, maxPos)
                real4 center = 0.5f*(maxPos+minPos);
#ifdef USE_PERIODIC
                if (valid) {
                    APPLY_PERIODIC_TO_POS_WITH_CENTER(pos, center)
                }
                else
                    pos = center;
                minPos = maxPos = pos;
                WARP_MINMAX3(minPos, maxPos)
                center = 0.5f*(maxPos+minPos);
#endif
                real4 blockSize = 0.5f*(maxPos-minPos);
                real4 delta = pos - center;
#ifdef USE_PERIODIC
                APPLY_PERIODIC_TO_DELTA(delta)
#endif
                real radius2 = (valid ? delta.x*delta.x+delta.y*delta.y+delta.z*delta.z : (real) 0);
                #pragma unroll
                for (int offset = 16; offset > 0; offset >>= 1)
                    radius2 = max(radius2, __shfl_xor_sync(0xffffffff, radius2, offset));
                center.w = sqrt(radius2);
                if (lane == 0) {
                    half3 hSize(trimTo3(blockSize));
                    real3 rSize = hSize.toReal3();
                    blockBoundingBox[index] = blockSize;
                    sortedBlockBoundingBox[index] = hSize;
                    sortedBlockCenter[index] = center;
                    shCenter[sub] = center;
                    bool isCompact = (rSize.x <= compactLimit && rSize.y <= compactLimit && rSize.z <= compactLimit
#ifdef USE_PERIODIC
                                      && 0.5f*periodicBoxSize.x-rSize.x >= PADDED_CUTOFF
                                      && 0.5f*periodicBoxSize.y-rSize.y >= PADDED_CUTOFF
                                      && 0.5f*periodicBoxSize.z-rSize.z >= PADDED_CUTOFF
#endif
                    );
                    shSize[sub] = make_real4(rSize.x, rSize.y, rSize.z, isCompact ? 1.0f : 0.0f);
                    blockCenter[index] = make_real4(center.x, center.y, center.z, 0);
                }
            }
        }
        __syncthreads();
        if (warpId == 0) {
            int b = sb * 32 + lane;
            bool validB = (b < NUM_BLOCKS);
            real4 c = (validB ? shCenter[lane] : make_real4(0, 0, 0, 0));
            real4 s = (validB ? shSize[lane] : make_real4(0, 0, 0, 1.0f));
            bool isLargeB = (validB && s.w == 0.0f), validCompact = (validB && s.w != 0.0f);
            unsigned int largeMask = BALLOT(isLargeB);
            unsigned int compactMask = BALLOT(validCompact);
            int refLane = (compactMask != 0 ? __ffs(compactMask) - 1 : (largeMask != 0 ? __ffs(largeMask) - 1 : 0));
            real4 c0 = make_real4(__shfl_sync(0xffffffff, c.x, refLane), __shfl_sync(0xffffffff, c.y, refLane), __shfl_sync(0xffffffff, c.z, refLane), 0);
            if (!validCompact) {
                c = c0;
                s = make_real4(0, 0, 0, 0);
            }
            real4 dc = c - c0;
#ifdef USE_PERIODIC
            if (validCompact) {
                APPLY_PERIODIC_TO_DELTA(dc)
            }
#endif
            real3 minP = make_real3(dc.x - s.x, dc.y - s.y, dc.z - s.z);
            real3 maxP = make_real3(dc.x + s.x, dc.y + s.y, dc.z + s.z);
            WARP_MINMAX3(minP, maxP)
            real4 mid = make_real4(c0.x + 0.5f*(minP.x + maxP.x), c0.y + 0.5f*(minP.y + maxP.y), c0.z + 0.5f*(minP.z + maxP.z), 0);
            real4 diff = (validCompact ? c - mid : make_real4(0, 0, 0, 0));
#ifdef USE_PERIODIC
            if (validCompact) {
                APPLY_PERIODIC_TO_DELTA(diff)
            }
#endif
            real3 ext = make_real3(fabs(diff.x) + s.x, fabs(diff.y) + s.y, fabs(diff.z) + s.z);
            WARP_MAX3(ext)
            if (lane == 0) {
                if (compactMask == 0)
                    mid = make_real4(1e10f, 1e10f, 1e10f, 0);
                sortedBlocks[2 * sb] = mid;
                sortedBlocks[2 * sb + 1] = make_real4(ext.x, ext.y, ext.z, PACK_UINT_TO_REAL(largeMask));
            }
            if (largeMask != 0) {
                int baseIdx = (lane == 0 ? (int) atomicAdd((unsigned int*) &rebuildNeighborList[2], (unsigned int) __popc(largeMask)) : 0);
                baseIdx = __shfl_sync(0xffffffff, baseIdx, 0);
                if (isLargeB)
                    ((int*)(sortedBlocks + 2 * numSuperblocks))[baseIdx + __popc(largeMask & ((1u << lane) - 1u))] = b;
            }
        }
        __syncthreads();
    }
}

extern "C" __global__ __launch_bounds__(GROUP_SIZE, 3) void findBlocksWithInteractions(real4 periodicBoxSize, real4 invPeriodicBoxSize, real4 periodicBoxVecX, real4 periodicBoxVecY, real4 periodicBoxVecZ,
        unsigned int* __restrict__ interactionCount, int* __restrict__ interactingTiles, unsigned int* __restrict__ interactingAtoms,
        int2* __restrict__ singlePairs, const real4* __restrict__ posq, unsigned int maxTiles, unsigned int maxSinglePairs,
        unsigned int startBlockIndex, unsigned int numBlocks, const real4* __restrict__ sortedBlocks, const real4* __restrict__ sortedBlockCenter,
        const half3* __restrict__ sortedBlockBoundingBox, const unsigned int* __restrict__ exclusionIndices, const unsigned int* __restrict__ exclusionRowIndices,
        const int* __restrict__ rebuildNeighborList, unsigned int rebuildStep) {

    if ((unsigned int) rebuildNeighborList[1] != rebuildStep)
        return; // The neighbor list doesn't need to be rebuilt.

    const int indexInWarp = threadIdx.x%32;
    const int warpStart = threadIdx.x-indexInWarp;
    const int warpMask = (1<<indexInWarp)-1;
    const int warpIndex = blockIdx.x*(GROUP_SIZE/32)+warpStart/32;
    const int totalWarps = gridDim.x*(GROUP_SIZE/32);
    const int numSuperblocks = (NUM_BLOCKS+31)/32;
    const int numLargeBlocks = rebuildNeighborList[2];
    const int* largeBlocksGlobal = (const int*)(sortedBlocks + 2 * numSuperblocks);
    const real compactLimit = 1.65f * PADDED_CUTOFF;
    __shared__ int workgroupBuffer[BUFFER_SIZE*(GROUP_SIZE/32)];
    __shared__ int workgroupFlagsBuffer[BUFFER_SIZE*(GROUP_SIZE/32)];
    __shared__ int warpExclusions[MAX_EXCLUSIONS*(GROUP_SIZE/32)];
    __shared__ real4 posBuffer[GROUP_SIZE];
    int* buffer = workgroupBuffer+BUFFER_SIZE*(warpStart/32);
    int* flagsBuffer = workgroupFlagsBuffer+BUFFER_SIZE*(warpStart/32);
    int* exclusionsForX = warpExclusions+MAX_EXCLUSIONS*(warpStart/32);

    for (int block1 = startBlockIndex + warpIndex; block1 < startBlockIndex + numBlocks; block1 += totalWarps) {
        int x = block1;
        int atom1 = x*TILE_SIZE+indexInWarp;
        real4 blockCenterX = sortedBlockCenter[block1];
        real3 blockSizeX = sortedBlockBoundingBox[block1].toReal3();
        int neighborsInBuffer = 0;
        real4 pos1 = posq[atom1];
#ifdef USE_PERIODIC
        const bool singlePeriodicCopy = (0.5f*periodicBoxSize.x-blockSizeX.x >= PADDED_CUTOFF &&
                                         0.5f*periodicBoxSize.y-blockSizeX.y >= PADDED_CUTOFF &&
                                         0.5f*periodicBoxSize.z-blockSizeX.z >= PADDED_CUTOFF);
#endif
        const bool isLargeX = (blockSizeX.x > compactLimit || blockSizeX.y > compactLimit || blockSizeX.z > compactLimit
#ifdef USE_PERIODIC
                               || !singlePeriodicCopy
#endif
        );
        if (!isLargeX) {
#ifdef USE_PERIODIC
            APPLY_PERIODIC_TO_POS_WITH_CENTER(pos1, blockCenterX)
#endif
            real3 pos1Local = trimTo3(pos1) - trimTo3(blockCenterX);
            posBuffer[threadIdx.x] = make_real4(pos1Local.x, pos1Local.y, pos1Local.z, 0.5f * (pos1Local.x * pos1Local.x + pos1Local.y * pos1Local.y + pos1Local.z * pos1Local.z));
        }
        else {
            posBuffer[threadIdx.x] = pos1;
        }

        const int exclusionStart = exclusionRowIndices[x];
        const int exclusionEnd = exclusionRowIndices[x+1];
        const int numExclusions = exclusionEnd-exclusionStart;
        #pragma unroll 4 // (MAX_EXCLUSIONS)
        for (int j = indexInWarp; j < numExclusions; j += 32)
            exclusionsForX[j] = exclusionIndices[exclusionStart+j];
        __syncwarp();

        const int startSb = (block1 + 1) >> 5;
        const int compactBatches = (isLargeX ? 0 : ((numSuperblocks - startSb + 31) >> 5));
        const int largeBatches = (numLargeBlocks + 31) >> 5;
        const int totalBatches = compactBatches + largeBatches;

        for (int batch = 0; batch < totalBatches; batch++) {
            int sb = startSb + (batch << 5) + indexInWarp;
            unsigned int laneLm = 0, sbFlags = 1u;
            if (batch < compactBatches) {
                bool includeSb = (sb < numSuperblocks);
                if (includeSb) {
                    real4 sbC = loadSortedBlockReal4(sortedBlocks, 2 * sb);
                    real4 sbS = loadSortedBlockReal4(sortedBlocks, 2 * sb + 1);
                    laneLm = UNPACK_REAL_TO_UINT(sbS.w);
                    unsigned int validMask = (sb == (block1 >> 5) ? ~((1u << ((block1 & 31) + 1)) - 1u) : 0xffffffffu);
                    includeSb = (((~laneLm) & validMask) != 0);
                    if (includeSb) {
                        real4 sbDelta = blockCenterX - sbC;
#ifdef USE_PERIODIC
                        APPLY_PERIODIC_TO_DELTA(sbDelta)
#endif
                        real3 cDelta = make_real3(
                            max(0.0f, fabs(sbDelta.x) - blockSizeX.x - sbS.x),
                            max(0.0f, fabs(sbDelta.y) - blockSizeX.y - sbS.y),
                            max(0.0f, fabs(sbDelta.z) - blockSizeX.z - sbS.z));
                        includeSb = (cDelta.x*cDelta.x + cDelta.y*cDelta.y + cDelta.z*cDelta.z < PADDED_CUTOFF_SQUARED);
#ifdef TRICLINIC
                        if (periodicBoxSize.z*0.5f - blockSizeX.z - sbS.z < PADDED_CUTOFF || periodicBoxSize.y*0.5f - blockSizeX.y - sbS.y < PADDED_CUTOFF)
                            includeSb = true;
#endif
                    }
                }
                sbFlags = BALLOT(includeSb);
            }
            while (sbFlags != 0) {
                int sbBit = __ffs(sbFlags) - 1;
                sbFlags &= sbFlags - 1;
                int block2;
                bool includeBlock2;
                if (batch < compactBatches) {
                    block2 = __shfl_sync(0xffffffff, sb, sbBit) * 32 + indexInWarp;
                    unsigned int sbLm = __shfl_sync(0xffffffff, laneLm, sbBit);
                    includeBlock2 = (block2 > block1 && block2 < NUM_BLOCKS) && (((sbLm >> indexInWarp) & 1u) == 0);
                }
                else {
                    int k = ((batch - compactBatches) << 5) + indexInWarp;
                    block2 = (k < numLargeBlocks ? __ldg(&largeBlocksGlobal[k]) : -1);
                    includeBlock2 = (isLargeX ? (block2 > block1) : (block2 >= 0));
                }
                bool forceInclude = false;
                real4 laneCenterY = make_real4(0, 0, 0, 0);
                real3 laneSizeY = make_real3(0, 0, 0);
                if (includeBlock2) {
                    laneCenterY = sortedBlockCenter[block2];
                    real4 blockDelta = blockCenterX-laneCenterY;
#ifdef USE_PERIODIC
                    APPLY_PERIODIC_TO_DELTA(blockDelta)
#endif
                    includeBlock2 &= (blockDelta.x*blockDelta.x+blockDelta.y*blockDelta.y+blockDelta.z*blockDelta.z < (PADDED_CUTOFF+blockCenterX.w+laneCenterY.w)*(PADDED_CUTOFF+blockCenterX.w+laneCenterY.w));
#ifdef TRICLINIC
                    laneSizeY = sortedBlockBoundingBox[block2].toReal3();
                    if (periodicBoxSize.z*0.5f-blockSizeX.z-laneSizeY.z < PADDED_CUTOFF || periodicBoxSize.y*0.5f-blockSizeX.y-laneSizeY.y < PADDED_CUTOFF)
                        includeBlock2 = forceInclude = true;
                    if (includeBlock2 && !forceInclude) {
#else
                    if (includeBlock2) {
                        laneSizeY = sortedBlockBoundingBox[block2].toReal3();
#endif
                        blockDelta.x = max(0.0f, fabs(blockDelta.x)-blockSizeX.x-laneSizeY.x);
                        blockDelta.y = max(0.0f, fabs(blockDelta.y)-blockSizeX.y-laneSizeY.y);
                        blockDelta.z = max(0.0f, fabs(blockDelta.z)-blockSizeX.z-laneSizeY.z);
                        includeBlock2 &= (blockDelta.x*blockDelta.x+blockDelta.y*blockDelta.y+blockDelta.z*blockDelta.z < PADDED_CUTOFF_SQUARED);
                    }
                    if (includeBlock2) {
                        #pragma unroll 4 // (MAX_EXCLUSIONS)
                        for (int k = 0; k < numExclusions; k++)
                            includeBlock2 &= (exclusionsForX[k] != block2);
                    }
                }

                int includeBlockFlags = BALLOT(includeBlock2);
                int forceIncludeFlags = BALLOT(forceInclude);
                while (includeBlockFlags != 0) {
                    int i = __ffs(includeBlockFlags)-1;
                    includeBlockFlags &= includeBlockFlags-1;
                    forceInclude = (forceIncludeFlags>>i) & 1;
                    int y = __shfl_sync(0xffffffff, block2, i);

                    real4 blockCenterY = make_real4(__shfl_sync(0xffffffff, laneCenterY.x, i), __shfl_sync(0xffffffff, laneCenterY.y, i), __shfl_sync(0xffffffff, laneCenterY.z, i), __shfl_sync(0xffffffff, laneCenterY.w, i));
                    real3 blockSizeY = make_real3(__shfl_sync(0xffffffff, laneSizeY.x, i), __shfl_sync(0xffffffff, laneSizeY.y, i), __shfl_sync(0xffffffff, laneSizeY.z, i));
                    real3 atomDelta = trimTo3(pos1)-trimTo3(blockCenterY);
#ifdef USE_PERIODIC
                    APPLY_PERIODIC_TO_DELTA(atomDelta)
#endif
                    real3 bdy = make_real3(max(0.0f, fabs(atomDelta.x)-blockSizeY.x), max(0.0f, fabs(atomDelta.y)-blockSizeY.y), max(0.0f, fabs(atomDelta.z)-blockSizeY.z));
                    unsigned int atomFlags = BALLOT(atom1 < NUM_ATOMS && (forceInclude || (atomDelta.x*atomDelta.x+atomDelta.y*atomDelta.y+atomDelta.z*atomDelta.z < (PADDED_CUTOFF+blockCenterY.w)*(PADDED_CUTOFF+blockCenterY.w) && bdy.x*bdy.x+bdy.y*bdy.y+bdy.z*bdy.z < PADDED_CUTOFF_SQUARED)));
                    if (atomFlags == 0)
                        continue;

                    int atom2 = y*TILE_SIZE+indexInWarp;
                    real4 pos2 = posq[atom2];
#ifdef USE_PERIODIC
                    if (!isLargeX) {
                        APPLY_PERIODIC_TO_POS_WITH_CENTER(pos2, blockCenterX)
                    }
#endif
                    real3 atom2Delta = trimTo3(pos2)-trimTo3(blockCenterX);
#ifdef USE_PERIODIC
                    if (isLargeX) {
                        APPLY_PERIODIC_TO_DELTA(atom2Delta)
                    }
#endif
                    real3 bdx = make_real3(max(0.0f, fabs(atom2Delta.x)-blockSizeX.x), max(0.0f, fabs(atom2Delta.y)-blockSizeX.y), max(0.0f, fabs(atom2Delta.z)-blockSizeX.z));
                    bool atom2InRange = (atom2 < NUM_ATOMS && (forceInclude || (atom2Delta.x*atom2Delta.x+atom2Delta.y*atom2Delta.y+atom2Delta.z*atom2Delta.z < (PADDED_CUTOFF+blockCenterX.w)*(PADDED_CUTOFF+blockCenterX.w) && bdx.x*bdx.x+bdx.y*bdx.y+bdx.z*bdx.z < PADDED_CUTOFF_SQUARED)));
                    unsigned int atom2Flags = BALLOT(atom2InRange);
                    int interacts = 0;
                    if (atom2Flags != 0 && atom2InRange) {
                        if (isLargeX) {
                            for (unsigned int f = atomFlags; f != 0; f &= f - 1) {
                                int j = __ffs(f) - 1;
                                real3 delta = trimTo3(pos2)-trimTo3(posBuffer[warpStart+j]);
#ifdef USE_PERIODIC
                                APPLY_PERIODIC_TO_DELTA(delta)
#endif
                                interacts |= (delta.x*delta.x+delta.y*delta.y+delta.z*delta.z < PADDED_CUTOFF_SQUARED ? 1<<j : 0);
                            }
                        }
                        else {
                            real atom2HalfMag2 = 0.5f * (atom2Delta.x * atom2Delta.x + atom2Delta.y * atom2Delta.y + atom2Delta.z * atom2Delta.z);
                            #define CHECK_SUBTILE(J0, J1) \
                                _Pragma("unroll") \
                                for (int j = (J0); j < (J1); j++) { \
                                    real4 posj = posBuffer[warpStart+j]; \
                                    if (posj.w + atom2HalfMag2 - posj.x*atom2Delta.x - posj.y*atom2Delta.y - posj.z*atom2Delta.z < 0.5f * PADDED_CUTOFF_SQUARED) \
                                        interacts |= 1<<j; }
                            if (atomFlags & 0x0000ffffu) { CHECK_SUBTILE(0, 16) }
                            if (atomFlags & 0xffff0000u) { CHECK_SUBTILE(16, TILE_SIZE) }
                            #undef CHECK_SUBTILE
                            interacts &= atomFlags;
                        }
                    }

                    int includeAtomFlags = BALLOT(interacts);
                    if (interacts) {
                        int index = neighborsInBuffer+__popc(includeAtomFlags&warpMask);
                        buffer[index] = atom2;
                        flagsBuffer[index] = interacts;
                    }
                    __syncwarp();
                    neighborsInBuffer += __popc(includeAtomFlags);
                    if (neighborsInBuffer > BUFFER_SIZE-TILE_SIZE) {
#if MAX_BITS_FOR_PAIRS > 0
                        neighborsInBuffer = saveSinglePairs(x, buffer, flagsBuffer, neighborsInBuffer, maxSinglePairs, &interactionCount[1], singlePairs);
#endif
                        unsigned int tilesToStore = neighborsInBuffer/TILE_SIZE;
                        if (tilesToStore > 0) {
                            unsigned int newTileStartIndex = (indexInWarp == 0 ? atomicAdd(&interactionCount[0], tilesToStore) : 0);
                            newTileStartIndex = __shfl_sync(0xffffffff, newTileStartIndex, 0);
                            if (newTileStartIndex+tilesToStore <= maxTiles) {
                                if (indexInWarp < tilesToStore)
                                    interactingTiles[newTileStartIndex+indexInWarp] = x;
                                #pragma unroll 8 // (GROUP_SIZE / TILE_SIZE)
                                for (int j = 0; j < tilesToStore; j++)
                                    interactingAtoms[(newTileStartIndex+j)*TILE_SIZE+indexInWarp] = buffer[indexInWarp+j*TILE_SIZE];
                            }
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
        }

#if MAX_BITS_FOR_PAIRS > 0
        if (neighborsInBuffer > 32)
            neighborsInBuffer = saveSinglePairs(x, buffer, flagsBuffer, neighborsInBuffer, maxSinglePairs, &interactionCount[1], singlePairs);
#endif
        if (neighborsInBuffer > 0) {
            unsigned int tilesToStore = (neighborsInBuffer+TILE_SIZE-1)/TILE_SIZE;
            unsigned int newTileStartIndex = (indexInWarp == 0 ? atomicAdd(&interactionCount[0], tilesToStore) : 0);
            newTileStartIndex = __shfl_sync(0xffffffff, newTileStartIndex, 0);
            if (newTileStartIndex+tilesToStore <= maxTiles) {
                if (indexInWarp < tilesToStore)
                    interactingTiles[newTileStartIndex+indexInWarp] = x;
                #pragma unroll 8 // (GROUP_SIZE / TILE_SIZE)
                for (int j = 0; j < tilesToStore; j++) {
                    int index = indexInWarp+j*TILE_SIZE;
                    interactingAtoms[(newTileStartIndex+j)*TILE_SIZE+indexInWarp] = (index < neighborsInBuffer ? buffer[index] : NUM_ATOMS);
                }
            }
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
    int* buffer = workgroupBuffer+BUFFER_SIZE*(warpStart/32);
    int* flagsBuffer = workgroupFlagsBuffer+BUFFER_SIZE*(warpStart/32);
    int* exclusionsForX = warpExclusions+MAX_EXCLUSIONS*(warpStart/32);

    // Loop over blocks.
    
    for (int block1 = startBlockIndex+warpIndex; block1 < startBlockIndex+numBlocks; block1 += totalWarps) {
        // Load data for this block.  Note that all threads in a warp are processing the same block.
        
        int x = sortedBlocks[block1] & BLOCK_INDEX_MASK;
        real4 blockCenterX = sortedBlockCenter[block1];
        real3 blockSizeX = sortedBlockBoundingBox[block1].toReal3();
        int neighborsInBuffer = 0;
        const unsigned int validAtomMask = BALLOT(x*TILE_SIZE+indexInWarp < NUM_ATOMS);
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
        __syncwarp();
        
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
                    interacts &= validAtomMask;
                }
                
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
                        unsigned int newTileStartIndex = (indexInWarp == 0 ? atomicAdd(&interactionCount[0], tilesToStore) : 0);
                        newTileStartIndex = __shfl_sync(0xffffffff, newTileStartIndex, 0);
                        if (newTileStartIndex+tilesToStore <= maxTiles) {
                            if (indexInWarp < tilesToStore)
                                interactingTiles[newTileStartIndex+indexInWarp] = x;
                            #pragma unroll 8 // (GROUP_SIZE / TILE_SIZE)
                            for (int j = 0; j < tilesToStore; j++)
                                interactingAtoms[(newTileStartIndex+j)*TILE_SIZE+indexInWarp] = buffer[indexInWarp+j*TILE_SIZE];
                        }
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
        
        // If we have a partially filled buffer,  store it to memory.
        
#if MAX_BITS_FOR_PAIRS > 0
        if (neighborsInBuffer > 32)
            neighborsInBuffer = saveSinglePairs(x, buffer, flagsBuffer, neighborsInBuffer, maxSinglePairs, &interactionCount[1], singlePairs);
#endif
        if (neighborsInBuffer > 0) {
            unsigned int tilesToStore = (neighborsInBuffer+TILE_SIZE-1)/TILE_SIZE;
            unsigned int newTileStartIndex = (indexInWarp == 0 ? atomicAdd(&interactionCount[0], tilesToStore) : 0);
            newTileStartIndex = __shfl_sync(0xffffffff, newTileStartIndex, 0);
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
