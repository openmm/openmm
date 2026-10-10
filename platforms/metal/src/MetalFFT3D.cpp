/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit.                   *
 * See https://openmm.org/development.                                        *
 *                                                                            *
 * Portions copyright (c) 2009-2026 Stanford University and the Authors.      *
 * Portions copyright (c) 2021 Advanced Micro Devices, Inc.                   *
 * Authors:                                                                   *
 * Contributors:                                                              *
 *                                                                            *
 * This program is free software: you can redistribute it and/or modify       *
 * it under the terms of the GNU Lesser General Public License as published   *
 * by the Free Software Foundation, either version 3 of the License, or       *
 * (at your option) any later version.                                        *
 *                                                                            *
 * This program is distributed in the hope that it will be useful,            *
 * but WITHOUT ANY WARRANTY; without even the implied warranty of             *
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the              *
 * GNU Lesser General Public License for more details.                        *
 *                                                                            *
 * You should have received a copy of the GNU Lesser General Public License   *
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.      *
 * -------------------------------------------------------------------------- */

#include "MetalFFT3D.h"
#include "MetalContext.h"
#include "MetalQueue.h"
#include <string>

using namespace OpenMM;
using namespace std;

MetalFFT3D::MetalFFT3D(MetalContext& context, int xsize, int ysize, int zsize, bool realToComplex) : context(context) {
    VkFFTConfiguration configuration = {};
    configuration.performR2C = realToComplex;
    configuration.device = &context.getDevice();
    configuration.queue = &dynamic_cast<MetalQueue*>(context.getCurrentQueue().get())->getQueue();
    configuration.doublePrecision = context.getUseDoublePrecision();

    configuration.FFTdim = 3;
    configuration.size[0] = zsize;
    configuration.size[1] = ysize;
    configuration.size[2] = xsize;

    configuration.inverseReturnToInputBuffer = true;
    configuration.isInputFormatted = true;
    configuration.inputBufferStride[0] = zsize;
    configuration.inputBufferStride[1] = configuration.inputBufferStride[0] * ysize;
    configuration.inputBufferStride[2] = configuration.inputBufferStride[1] * xsize;

    configuration.bufferStride[0] = realToComplex ? (zsize/2 + 1) : zsize;
    configuration.bufferStride[1] = configuration.bufferStride[0] * ysize;
    configuration.bufferStride[2] = configuration.bufferStride[1] * xsize;

    app = new VkFFTApplication();
    VkFFTResult fftResult = initializeVkFFT(app, configuration);
    if (fftResult != VKFFT_SUCCESS) {
        throw OpenMMException(string("Error executing initializeVkFFT: ")+getVkFFTErrorString(fftResult));
    }
}

MetalFFT3D::~MetalFFT3D() {
    deleteVkFFT(app);
    delete app;
}

void MetalFFT3D::execFFT(ArrayInterface& in, ArrayInterface& out, bool forward) {
    VkFFTLaunchParams params = {};
    MTL::Buffer* inputBuffer = context.unwrap(in).getBuffer();
    MTL::Buffer* outputBuffer = context.unwrap(out).getBuffer();
    if (forward) {
        params.inputBuffer = &inputBuffer;
        params.buffer = &outputBuffer;
    }
    else {
        params.inputBuffer = &outputBuffer;
        params.buffer = &inputBuffer;
    }
    MetalQueue& queue = *dynamic_cast<MetalQueue*>(context.getCurrentQueue().get());
    params.commandBuffer = &queue.getCommandBuffer();
    params.commandEncoder = &queue.getEncoder();
    VkFFTResult fftResult = VkFFTAppend(app, forward ? -1 : 1, &params);
    if (fftResult != VKFFT_SUCCESS) {
        throw OpenMMException(string("Error executing VkFFTAppend: ")+getVkFFTErrorString(fftResult));
    }
}

int MetalFFT3D::findLegalDimension(int minimum) {
    if (minimum < 1)
        return 1;
    while (true) {
        // Attempt to factor the current value.

        int unfactored = minimum;
        // VkFFT supports prime factors up to 13
        for (int factor = 2; factor <= 13; factor++) {
            while (unfactored > 1 && unfactored%factor == 0)
                unfactored /= factor;
        }
        if (unfactored == 1)
            return minimum;
        minimum++;
    }
}
