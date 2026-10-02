/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * This is part of the OpenMM molecular simulation toolkit.                   *
 * See https://openmm.org/development.                                        *
 *                                                                            *
 * Ported from platforms/opencl/src/kernels/utilities.cl                      *
 * (determineNativeAccuracy).                                                *
 * Original OpenCL Platform code:                                            *
 * Portions copyright (c) Stanford University and the Authors.                *
 * Authors: Peter Eastman                                                     *
 * Metal Platform code:                                                       *
 * Portions copyright (c) 2026 Chun-Chi Hung.                                 *
 * Authors: Chun-Chi Hung                                                     *
 *                                                                            *
 * This program is free software: you can redistribute it and/or modify       *
 * it under the terms of the GNU Lesser General Public License as published   *
 * by the Free Software Foundation, either version 3 of the License, or       *
 * (at your option) any later version.                                        *
 * This program is distributed WITHOUT ANY WARRANTY; without even the        *
 * implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. *
 * See the GNU Lesser General Public License for more details.                *
 * You should have received a copy of the GNU Lesser General Public License  *
 * along with this program. If not, see <http://www.gnu.org/licenses/>.       *
 * -------------------------------------------------------------------------- */

#include <metal_stdlib>
using namespace metal;

/** @brief Evaluate OpenCL's native-function probe using two float4s per sample. */
kernel void determineNativeAccuracy(device float4* values [[buffer(0)]],
        uint index [[thread_position_in_grid]]) {
    if (index >= 20)
        return;
    float v = values[2*index].x;
    values[2*index] = float4(v, fast::sqrt(v), fast::rsqrt(v), fast::divide(1.0f, v));
    values[2*index+1] = float4(fast::exp(v), fast::log(v), 0.0f, 0.0f);
}
