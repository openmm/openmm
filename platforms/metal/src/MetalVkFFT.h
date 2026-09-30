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

#ifndef OPENMM_METALVKFFT_H_
#define OPENMM_METALVKFFT_H_

#include <memory>

namespace OpenMM {

/**
 * @brief Private C++11-compatible boundary around VkFFT and Metal-cpp.
 *
 * Native objects are borrowed opaque handles. Only the implementation file
 * includes Metal-cpp and requires C++17; neither appears in public headers.
 */
class MetalVkFFT {
public:
    /** @brief Initialize an unnormalized, single-precision GPU transform plan. */
    MetalVkFFT(void* device, void* queue, int xsize, int ysize, int zsize, bool realToComplex);
    /** @brief Free VkFFT-owned resources without releasing borrowed runtime objects. */
    ~MetalVkFFT();
    /**
     * @brief Encode the complete transform without ending, submitting, or waiting.
     * @param commandBuffer Borrowed native command buffer owning the encoder.
     * @param encoder Borrowed serial compute encoder; caller ends it.
     * @param input Borrowed native MTLBuffer containing the transform input.
     * @param output Borrowed distinct MTLBuffer for the transform result.
     * @param forward True for forward, false for inverse.
     */
    void append(void* commandBuffer, void* encoder, void* input, void* output, bool forward);
private:
    class Impl;
    std::unique_ptr<Impl> impl;
};

} // namespace OpenMM

#endif // OPENMM_METALVKFFT_H_
