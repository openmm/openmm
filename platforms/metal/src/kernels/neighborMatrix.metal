/* -------------------------------------------------------------------------- *
 * OpenMM — Metal Platform                                                    *
 * Portions copyright (c) 2026 Chun-Chi Hung.                                  *
 * Authors: Chun-Chi Hung                                                      *
 *                                                                            *
 * This is part of the OpenMM molecular simulation toolkit.                    *
 * See https://openmm.org/development.                                         *
 * This program is free software: you can redistribute it and/or modify        *
 * it under the terms of the GNU Lesser General Public License as published    *
 * by the Free Software Foundation, either version 3 of the License, or         *
 * (at your option) any later version.                                         *
 * This program is distributed WITHOUT ANY WARRANTY; without even the          *
 * implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.   *
 * See the GNU Lesser General Public License for more details.                 *
 * You should have received a copy of the GNU Lesser General Public License    *
 * along with this program. If not, see <http://www.gnu.org/licenses/>.          *
 * -------------------------------------------------------------------------- */

/**
 * @file
 * @brief Test-only 32-by-32 matrix distance experiment, not a neighbor list.
 *
 * MSL 3.0, one complete 32-thread SIMD group per threadgroup. Matrix loads,
 * stores, and all sixteen 8-by-8 multiply-accumulates have uniform control flow.
 * No mapping from matrix elements to hardware lanes is assumed.
 *
 * Squared distances formed from norms can lose precision. Independently
 * wrapping positions is also not a pairwise minimum-image operation. Therefore
 * candidate masks are diagnostic only: EVERY pair receives the scalar check,
 * including matrix-rejected pairs. This source is not used by any Force kernel.
 */

#include <metal_stdlib>
#include <metal_simdgroup_matrix>
using namespace metal;

#ifndef OPENMM_METAL_EXPERIMENTAL_MATRIX_SCREEN
#define OPENMM_METAL_EXPERIMENTAL_MATRIX_SCREEN 0
#endif
#ifndef METAL_MATRIX_DIAGNOSTICS
#define METAL_MATRIX_DIAGNOSTICS 0
#endif

#if OPENMM_METAL_EXPERIMENTAL_MATRIX_SCREEN

/** @brief Explicit 64-byte host/shader ABI; all box vectors use float4 storage. */
struct MatrixScreenParameters {
    float4 boxX, boxY, boxZ;
    float cutoffSquared;
    uint periodic, tileCount, reserved;
};

/** @brief The same reduced-triclinic sequential minimum image used by Common. */
inline float3 matrixScreenDelta(float3 delta, constant MatrixScreenParameters& params) {
    if (params.periodic) {
        delta -= floor(delta.z/params.boxZ.z+0.5f)*params.boxZ.xyz;
        delta -= floor(delta.y/params.boxY.y+0.5f)*params.boxY.xyz;
        delta -= floor(delta.x/params.boxX.x+0.5f)*params.boxX.xyz;
    }
    return delta;
}

/** @brief Scalar subtract-first predicate, authoritative in BOTH experimental paths. */
inline float matrixScreenDistance(float3 a, float3 b, constant MatrixScreenParameters& params) {
    float3 delta = matrixScreenDelta(b-a, params);
    return delta.x*delta.x+delta.y*delta.y+delta.z*delta.z;
}

/** @brief Resident-data scalar reference, including the same output mask writes. */
kernel void scalarScreen(device const float4* first [[buffer(0)]],
        device const float4* second [[buffer(1)]], device uint* exactMasks [[buffer(2)]],
        device uint* candidateMasks [[buffer(3)]], device float* distances [[buffer(4)]],
        constant MatrixScreenParameters& params [[buffer(5)]], device uint* shape [[buffer(6)]],
        uint lane [[thread_index_in_threadgroup]], uint group [[threadgroup_position_in_grid]],
        uint groups [[threadgroups_per_grid]], uint simdWidth [[threads_per_simdgroup]],
        uint groupWidth [[threads_per_threadgroup]]) {
    if (group == 0 && lane == 0) {
        shape[0] = simdWidth;
        shape[1] = groupWidth;
    }
    if (simdWidth != 32 || groupWidth != 32)
        return;
    for (uint tile = group; tile < params.tileCount; tile += groups) {
        uint mask = 0;
        float3 a = first[32*tile+lane].xyz;
        for (uint j = 0; j < 32; j++) {
            float distance = matrixScreenDistance(a, second[32*tile+j].xyz, params);
            if (distance < params.cutoffSquared)
                mask |= 1u<<j;
#if METAL_MATRIX_DIAGNOSTICS
            distances[1024*tile+32*lane+j] = distance;
#endif
        }
        exactMasks[32*tile+lane] = mask;
        candidateMasks[32*tile+lane] = mask;
    }
}

/**
 * @brief Pack positions and norm terms, perform sixteen MMAs, then recheck ALL pairs.
 * @note The benchmark includes this packing, shared storage, masks, and exact recheck.
 */
kernel void matrixScreen(device const float4* first [[buffer(0)]],
        device const float4* second [[buffer(1)]], device uint* exactMasks [[buffer(2)]],
        device uint* candidateMasks [[buffer(3)]], device float* distances [[buffer(4)]],
        constant MatrixScreenParameters& params [[buffer(5)]], device uint* shape [[buffer(6)]],
        uint lane [[thread_index_in_threadgroup]], uint group [[threadgroup_position_in_grid]],
        uint groups [[threadgroups_per_grid]], uint simdWidth [[threads_per_simdgroup]],
        uint groupWidth [[threads_per_threadgroup]]) {
    if (group == 0 && lane == 0) {
        shape[0] = simdWidth;
        shape[1] = groupWidth;
    }
    // This guard is uniform and executes before any matrix operation/barrier.
    if (simdWidth != 32 || groupWidth != 32)
        return;
    threadgroup float left[32*8], right[8*32], squared[32*32];
    for (uint tile = group; tile < params.tileCount; tile += groups) {
        // A shared origin reduces, but does not bound, cancellation. Independently
        // wrapped coordinates still need the pairwise periodic scalar recheck.
        float3 origin = first[32*tile].xyz;
        float3 a = matrixScreenDelta(first[32*tile+lane].xyz-origin, params);
        float3 b = matrixScreenDelta(second[32*tile+lane].xyz-origin, params);
        float normA = a.x*a.x+a.y*a.y+a.z*a.z;
        float normB = b.x*b.x+b.y*b.y+b.z*b.z;
        for (uint k = 0; k < 8; k++) {
            left[8*lane+k] = k < 3 ? -2.0f*a[k] : k == 3 ? normA : k == 4 ? 1.0f : 0.0f;
            right[32*k+lane] = k < 3 ? b[k] : k == 3 ? 1.0f : k == 4 ? normB : 0.0f;
        }
        threadgroup_barrier(mem_flags::mem_threadgroup);
        for (uint row = 0; row < 32; row += 8) {
            for (uint column = 0; column < 32; column += 8) {
                simdgroup_matrix<float, 8, 8> l, r, result;
                simdgroup_load(l, left+8*row, 8);
                simdgroup_load(r, right+column, 32);
                simdgroup_matrix<float, 8, 8> zero = make_filled_simdgroup_matrix<float, 8>(0.0f);
                simdgroup_multiply_accumulate(result, l, r, zero);
                simdgroup_store(result, squared+32*row+column, 32);
            }
        }
        threadgroup_barrier(mem_flags::mem_threadgroup);
        uint candidate = 0, exact = 0;
        for (uint j = 0; j < 32; j++) {
            float distance = squared[32*lane+j];
            if (distance < params.cutoffSquared)
                candidate |= 1u<<j;
            // Do not gate this on candidate: false negatives are not yet bounded.
            if (matrixScreenDistance(first[32*tile+lane].xyz, second[32*tile+j].xyz, params) < params.cutoffSquared)
                exact |= 1u<<j;
#if METAL_MATRIX_DIAGNOSTICS
            distances[1024*tile+32*lane+j] = distance;
#endif
        }
        candidateMasks[32*tile+lane] = candidate;
        exactMasks[32*tile+lane] = exact;
        threadgroup_barrier(mem_flags::mem_threadgroup);
    }
}
#endif
