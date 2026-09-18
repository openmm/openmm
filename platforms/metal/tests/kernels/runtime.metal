/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 *                                                                            *
 * This is part of the OpenMM molecular simulation toolkit.                   *
 * See https://openmm.org/development.                                        *
 *                                                                            *
 * Metal Platform code:                                                       *
 * Portions copyright (c) 2026 Chun-Chi Hung.                                 *
 * Authors: Chun-Chi Hung                                                     *
 *                                                                            *
 * This program is free software: you can redistribute it and/or modify       *
 * it under the terms of the GNU Lesser General Public License as published   *
 * by the Free Software Foundation, either version 3 of the License, or       *
 * (at your option) any later version.                                        *
 *                                                                            *
 * This program is distributed in the hope that it will be useful,            *
 * but WITHOUT ANY WARRANTY; without even the implied warranty of             *
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the               *
 * GNU Lesser General Public License for more details.                        *
 *                                                                            *
 * You should have received a copy of the GNU Lesser General Public License   *
 * along with this program. If not, see <http://www.gnu.org/licenses/>.       *
 * -------------------------------------------------------------------------- */

/**
 * @file
 * @brief Test-only native MSL kernels with buffer slots in Common argument order.
 */

#include <metal_stdlib>
using namespace metal;

/**
 * @brief Apply an integer transform while exercising array and scalar bindings.
 * @param input Source array at buffer 0, containing at least count integers.
 * @param output Separate destination array at buffer 1, containing at least count integers.
 * @param count Nonnegative element count at buffer 2.
 * @param offset Additive constant at buffer 3.
 * @param globalID Metal-provided thread position in the one-dimensional grid.
 * @param numGroups Metal-provided number of threadgroups in the grid.
 * @param groupWidth Metal-provided threads per threadgroup.
 * @note VALUE_SCALE is supplied at compilation. The grid-stride loop supports
 *       both padded and capped launches; stage builtins consume no buffer slots.
 */
kernel void transform(device const int* input [[buffer(0)]],
        device int* output [[buffer(1)]],
        constant int& count [[buffer(2)]],
        constant int& offset [[buffer(3)]],
        uint globalID [[thread_position_in_grid]],
        uint numGroups [[threadgroups_per_grid]],
        uint groupWidth [[threads_per_threadgroup]]) {
    for (uint i = globalID; i < uint(count); i += numGroups*groupWidth)
        output[i] = VALUE_SCALE*input[i]+offset;
}

/**
 * @brief Count participating threads to test complete threadgroup dispatch.
 * @param output Buffer 0, with one uint result per dispatched threadgroup.
 * @param localID Metal-provided thread index within the threadgroup.
 * @param groupID Metal-provided threadgroup position in the one-dimensional grid.
 * @pre Launch with exactly 64 threads per group.
 * @note All threads participate in both barriers. Zeroing all 64 entries makes
 *       an accidental partial final group detectable without reading uninitialized memory.
 */
kernel void recordGroupWidth(device uint* output [[buffer(0)]],
        uint localID [[thread_index_in_threadgroup]],
        uint groupID [[threadgroup_position_in_grid]]) {
    threadgroup uint lanes[64];
    if (localID == 0)
        for (uint i = 0; i < 64; i++)
            lanes[i] = 0;
    threadgroup_barrier(mem_flags::mem_threadgroup);
    lanes[localID] = 1;
    threadgroup_barrier(mem_flags::mem_threadgroup);
    if (localID == 0) {
        uint sum = 0;
        for (uint i = 0; i < 64; i++)
            sum += lanes[i];
        output[groupID] = sum;
    }
}
