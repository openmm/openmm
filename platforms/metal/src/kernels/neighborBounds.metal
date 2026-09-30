/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * Ported from the OpenMM OpenCL and HIP Platforms:
 * platforms/opencl/src/kernels/findInteractingBlocks.cl
 * platforms/hip/src/kernels/findInteractingBlocks.hip
 * Original code: Portions copyright (c) 2009-2025 Stanford University and
 * the Authors. Authors: Peter Eastman and the OpenMM contributors.
 * Metal code: Portions copyright (c) 2026 Chun-Chi Hung.
 * Author: Chun-Chi Hung.
 * This program is free software under the GNU Lesser General Public License,
 * version 3 or (at your option) any later version. It is distributed without
 * any warranty; see <http://www.gnu.org/licenses/> for the license.
 * -------------------------------------------------------------------------- */

/**
 * @brief Compute each block cooperatively, following HIP's warp bounds.
 * Periodic wrapping retains OpenCL's ordered bounding-box expansion because
 * each nearest-image choice depends on the preceding partial box. Uniform
 * SIMD loops include inactive tail lanes by broadcasting the last valid atom.
 * The final 64-thread reduction preserves the OpenCL blockSizeRange ABI.
 */
KERNEL void findBlockBounds(int numAtoms, real4 periodicBoxSize, real4 invPeriodicBoxSize,
        real4 periodicBoxVecX, real4 periodicBoxVecY, real4 periodicBoxVecZ,
        GLOBAL const real4* posq, GLOBAL real4* blockCenter,
        GLOBAL real4* blockBoundingBox, GLOBAL int* rebuildNeighborList,
        GLOBAL real2* blockSizeRange) {
    const int lane = LOCAL_ID%32;
    real minSize = 1e38f, maxSize = 0;
    for (int index = GLOBAL_ID/32; index*32 < numAtoms; index += GLOBAL_SIZE/32) {
        int base = index*32;
        int count = min(32, numAtoms-base);
        real4 position = posq[base+min(lane, count-1)];
        real4 minPos, maxPos;
#ifdef USE_PERIODIC
        real4 pos = simdShuffle(position, 0);
        APPLY_PERIODIC_TO_POS(pos)
        minPos = maxPos = pos;
        for (int i = 1; i < count; i++) {
            pos = simdShuffle(position, i);
            real4 center = 0.5f*(minPos+maxPos);
            APPLY_PERIODIC_TO_POS_WITH_CENTER(pos, center)
            minPos = min(minPos, pos);
            maxPos = max(maxPos, pos);
        }
#else
        minPos = make_real4(simdReduceMin(position.x), simdReduceMin(position.y), simdReduceMin(position.z), simdReduceMin(position.w));
        maxPos = make_real4(simdReduceMax(position.x), simdReduceMax(position.y), simdReduceMax(position.z), simdReduceMax(position.w));
#endif
        real4 center = 0.5f*(minPos+maxPos);
        real4 width = 0.5f*(maxPos-minPos);
        real4 delta = position-center;
#ifdef USE_PERIODIC
        APPLY_PERIODIC_TO_DELTA(delta)
#endif
        real radius = dot(delta.xyz, delta.xyz);
        center.w = sqrt(simdReduceMax(radius));
        if (lane == 0) {
            blockCenter[index] = center;
            blockBoundingBox[index] = width;
        }
        real size = width.x+width.y+width.z;
        minSize = min(minSize, size);
        maxSize = max(maxSize, size);
    }
    minSize = simdReduceMin(minSize);
    maxSize = simdReduceMax(maxSize);
    LOCAL real2 ranges[2];
    if (lane == 0)
        ranges[LOCAL_ID/32] = make_real2(minSize, maxSize);
    SYNC_THREADS;
    if (LOCAL_ID == 0)
        blockSizeRange[GROUP_ID] = make_real2(min(ranges[0].x, ranges[1].x), max(ranges[0].y, ranges[1].y));
    if (GLOBAL_ID == 0)
        rebuildNeighborList[0] = 0;
}
