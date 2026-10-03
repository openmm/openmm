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

#include <metal_stdlib>
using namespace metal;

// A unit-length real-to-complex axis has no redundant complex entries.
// Use a complex FFT for the remaining axes, with GPU-only layout conversion.
/**
 * @brief Expand contiguous real values into a distinct complex buffer.
 * @param input Read-only real grid at buffer 0.
 * @param output Complex grid at buffer 1, with zero imaginary components.
 * @param count Number of logical grid values at buffer 2.
 * @param id Global position within the one-dimensional dispatch.
 * @param size Total dispatched threads; the loop also handles capped grids.
 */
kernel void packRealFFT(device const float* input [[buffer(0)]],
        device float2* output [[buffer(1)]], constant int& count [[buffer(2)]],
        uint id [[thread_position_in_grid]], uint size [[threads_per_grid]]) {
    for (uint i = id; i < uint(count); i += size)
        output[i] = float2(input[i], 0.0f);
}

/**
 * @brief Extract real components into a distinct contiguous real buffer.
 * @param input Read-only complex grid at buffer 0.
 * @param output Real grid at buffer 1.
 * @param count Number of logical grid values at buffer 2.
 * @param id Global position within the one-dimensional dispatch.
 * @param size Total dispatched threads; the loop also handles capped grids.
 */
kernel void unpackRealFFT(device const float2* input [[buffer(0)]],
        device float* output [[buffer(1)]], constant int& count [[buffer(2)]],
        uint id [[thread_position_in_grid]], uint size [[threads_per_grid]]) {
    for (uint i = id; i < uint(count); i += size)
        output[i] = input[i].x;
}
