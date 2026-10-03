/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 *                                                                            *
 * This is part of the OpenMM molecular simulation toolkit.                   *
 * See https://openmm.org/development.                                        *
 *                                                                            *
 * Ported from the OpenMM OpenCL Platform's VkFFT integration.                *
 * Source: platforms/opencl/src/OpenCLFFT3D.cpp                               *
 *                                                                            *
 * Original OpenCL Platform code:                                             *
 * Portions copyright (c) 2009-2025 Stanford University and the Authors.      *
 * Authors: Peter Eastman                                                     *
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

#include "MetalVkFFT.h"
#include "openmm/OpenMMException.h"
// VkFFT defines Metal-cpp's implementation macros: include it in this TU only.
#define VKFFT_BACKEND 5
// Keep FFT compilation under the same build-time policy as other Metal kernels.
#define VKFFT_METAL_FAST_MATH OPENMM_METAL_FAST_MATH
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wdeprecated-declarations"
#include "vkFFT.h"
#pragma clang diagnostic pop

using namespace OpenMM;

class MetalVkFFT::Impl {
public:
    VkFFTApplication application = {};
    // VkFFT keeps pointers to these slots after append() returns.
    MTL::Buffer* inputBuffer = nullptr;
    MTL::Buffer* buffer = nullptr;
    ~Impl() {
        deleteVkFFT(&application);
    }
};

MetalVkFFT::MetalVkFFT(void* device, void* queue, int xsize, int ysize, int zsize, bool realToComplex)
        : impl(new Impl()) {
    VkFFTConfiguration config = {};
    // Unit axes do not change storage order, but VkFFT 1.2.33 cannot generate
    // their radix-free kernels. Keep the real-to-complex axis in its first slot.
    const int dimensions[] = {zsize, ysize, xsize};
    uint64_t stride = 1;
    for (int i = 0; i < 3; i++) {
        if (dimensions[i] > 1 || (i == 0 && realToComplex)) {
            config.size[config.FFTdim] = dimensions[i];
            stride *= dimensions[i];
            config.inputBufferStride[config.FFTdim] = stride;
            config.FFTdim++;
        }
    }
    config.performR2C = realToComplex;
    config.device = static_cast<MTL::Device*>(device);
    config.queue = static_cast<MTL::CommandQueue*>(queue);
    config.inverseReturnToInputBuffer = true;
    config.isInputFormatted = 1;
    VkFFTResult result = initializeVkFFT(&impl->application, config);
    if (result != VKFFT_SUCCESS)
        throw OpenMMException(std::string("Error initializing Metal VkFFT: ")+getVkFFTErrorString(result));
}

MetalVkFFT::~MetalVkFFT() {
}

void MetalVkFFT::append(void* commandBuffer, void* encoder, void* input, void* output, bool forward) {
    impl->inputBuffer = static_cast<MTL::Buffer*>(forward ? input : output);
    impl->buffer = static_cast<MTL::Buffer*>(forward ? output : input);
    VkFFTLaunchParams params = {};
    params.inputBuffer = &impl->inputBuffer;
    params.buffer = &impl->buffer;
    params.commandBuffer = static_cast<MTL::CommandBuffer*>(commandBuffer);
    params.commandEncoder = static_cast<MTL::ComputeCommandEncoder*>(encoder);
    VkFFTResult result = VkFFTAppend(&impl->application, forward ? -1 : 1, &params);
    if (result != VKFFT_SUCCESS)
        throw OpenMMException(std::string("Error executing Metal VkFFT: ")+getVkFFTErrorString(result));
}
