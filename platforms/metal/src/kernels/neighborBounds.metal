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
 * @brief Use SIMD bounds for nonperiodic blocks and ordered OpenCL periodic bounds.
 * Periodic nearest-image choices depend on the preceding partial box, so keep
 * one independent block per lane instead of repeating that serial expansion
 * in every lane of a SIMD group. Nonperiodic blocks use one SIMD group each;
 * tail lanes repeat the last valid atom and participate in every reduction.
 * Both paths reduce size ranges with SIMD operations and one threadgroup
 * barrier, preserving the OpenCL blockSizeRange ABI and 64-thread groups.
 */
KERNEL void findBlockBounds(int numAtoms, real4 periodicBoxSize, real4 invPeriodicBoxSize,
        real4 periodicBoxVecX, real4 periodicBoxVecY, real4 periodicBoxVecZ,
        GLOBAL const real4* posq, GLOBAL real4* blockCenter,
        GLOBAL real4* blockBoundingBox, GLOBAL int* rebuildNeighborList,
        GLOBAL real2* blockSizeRange) {
    const int lane = LOCAL_ID%32;
    real minSize = 1e38f, maxSize = 0;
#ifdef USE_PERIODIC
    // Match OpenCL's atom order, image convention, and radius expression.
    for (int index = GLOBAL_ID; index*32 < numAtoms; index += GLOBAL_SIZE) {
        int base = index*32;
        int last = min(base+32, numAtoms);
        real4 pos = posq[base];
        APPLY_PERIODIC_TO_POS(pos)
        real4 minPos = pos, maxPos = pos;
        for (int i = base+1; i < last; i++) {
            pos = posq[i];
            real4 center = 0.5f*(maxPos+minPos);
            APPLY_PERIODIC_TO_POS_WITH_CENTER(pos, center)
            minPos = min(minPos, pos);
            maxPos = max(maxPos, pos);
        }
        real4 width = 0.5f*(maxPos-minPos);
        real4 center = 0.5f*(maxPos+minPos);
        center.w = 0;
        for (int i = base; i < last; i++) {
            real4 delta = posq[i]-center;
            APPLY_PERIODIC_TO_DELTA(delta)
            center.w = max(center.w, delta.x*delta.x+delta.y*delta.y+delta.z*delta.z);
        }
        center.w = sqrt(center.w);
        blockCenter[index] = center;
        blockBoundingBox[index] = width;
        real size = width.x+width.y+width.z;
        minSize = min(minSize, size);
        maxSize = max(maxSize, size);
    }
#else
    for (int index = GLOBAL_ID/32; index*32 < numAtoms; index += GLOBAL_SIZE/32) {
        int base = index*32;
        int count = min(32, numAtoms-base);
        real4 position = posq[base+min(lane, count-1)];
        real4 minPos = make_real4(simdReduceMin(position.x), simdReduceMin(position.y), simdReduceMin(position.z), simdReduceMin(position.w));
        real4 maxPos = make_real4(simdReduceMax(position.x), simdReduceMax(position.y), simdReduceMax(position.z), simdReduceMax(position.w));
        real4 center = 0.5f*(minPos+maxPos);
        real4 width = 0.5f*(maxPos-minPos);
        real4 delta = position-center;
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
#endif
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
