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

#include "MetalFFT3D.h"
#include "MetalContext.h"
#include "MetalKernelSources.h"
#include "MetalQueue.h"
#include "MetalVkFFT.h"
#include "openmm/OpenMMException.h"
#import <Metal/Metal.h>
#include <algorithm>
#include <limits>

using namespace OpenMM;
using namespace std;

MetalFFT3D::MetalFFT3D(MetalContext& context, int xsize, int ysize, int zsize, bool realToComplex)
        : context(context), requiredBytes(0), realToComplex(realToComplex),
          packRealAsComplex(realToComplex && zsize == 1) {
    if (xsize <= 0 || ysize <= 0 || zsize <= 0)
        throw OpenMMException("Metal FFT dimensions must be positive");
    size_t count = xsize;
    if (count > size_t(numeric_limits<int>::max())/ysize)
        throw OpenMMException("Metal FFT dimensions are too large");
    count *= ysize;
    if (count > size_t(numeric_limits<int>::max())/zsize)
        throw OpenMMException("Metal FFT dimensions are too large");
    requiredBytes = count*zsize*2*sizeof(float);
    // VkFFT has no radix stage to generate for the one-point identity transform.
    if (requiredBytes == 2*sizeof(float))
        return;
    if (packRealAsComplex) {
        complexGrid.initialize<mm_float2>(context, count*zsize, "fftComplexGrid");
        ComputeProgram program = context.compileProgram(MetalKernelSources::fft);
        packKernel = program->createKernel("packRealFFT");
        unpackKernel = program->createKernel("unpackRealFFT");
        for (ComputeKernel kernel : {packKernel, unpackKernel}) {
            kernel->addArg();
            kernel->addArg();
            kernel->addArg(int(count*zsize));
        }
    }
    @autoreleasepool {
        plan.reset(new MetalVkFFT(context.getDevice(), context.getCurrentMetalQueue().getQueue(),
                xsize, ysize, zsize, realToComplex && !packRealAsComplex));
    }
}

MetalFFT3D::~MetalFFT3D() {
}

void MetalFFT3D::execFFT(ArrayInterface& in, ArrayInterface& out, bool forward) {
    MetalArray& input = context.unwrap(in);
    MetalArray& output = context.unwrap(out);
    if (input.getBuffer() == output.getBuffer())
        throw OpenMMException("Metal FFT input and output arrays must be different");
    if (input.getSize()*input.getElementSize() < requiredBytes ||
            output.getSize()*output.getElementSize() < requiredBytes)
        throw OpenMMException("Metal FFT arrays must hold the full complex grid");
    MetalArray* transformInput = &input;
    MetalArray* transformOutput = &output;
    if (packRealAsComplex && plan) {
        if (forward) {
            packKernel->setArg(0, input);
            packKernel->setArg(1, complexGrid);
            packKernel->execute(complexGrid.getSize());
            transformInput = &context.unwrap(complexGrid);
        }
        else
            transformOutput = &context.unwrap(complexGrid);
    }
    @autoreleasepool {
        MetalQueue& queue = context.getCurrentMetalQueue();
        auto queueLock = queue.lock();
        id<MTLCommandBuffer> command = (__bridge id<MTLCommandBuffer>) queue.getCommandBuffer();
        if (!plan) {
            id<MTLBlitCommandEncoder> encoder = [command blitCommandEncoder];
            if (encoder == nil)
                throw OpenMMException("Error creating Metal FFT identity encoder");
            id<MTLBuffer> inputBuffer = (__bridge id<MTLBuffer>) input.getBuffer();
            id<MTLBuffer> outputBuffer = (__bridge id<MTLBuffer>) output.getBuffer();
            [encoder copyFromBuffer:inputBuffer sourceOffset:0 toBuffer:outputBuffer
                    destinationOffset:0 size:(realToComplex ? sizeof(float) : 2*sizeof(float))];
            if (realToComplex && forward)
                [encoder fillBuffer:outputBuffer range:NSMakeRange(sizeof(float), sizeof(float)) value:0];
            [encoder endEncoding];
            queue.submit((__bridge void*) command);
            return;
        }
        // VkFFT appends dependent FFT passes to one serial compute encoder.
        id<MTLComputeCommandEncoder> encoder = [command computeCommandEncoder];
        if (encoder == nil)
            throw OpenMMException("Error creating Metal FFT compute encoder");
        try {
            plan->append((__bridge void*) command, (__bridge void*) encoder,
                    transformInput->getBuffer(), transformOutput->getBuffer(), forward);
        }
        catch (...) {
            [encoder endEncoding];
            throw;
        }
        [encoder endEncoding];
        queue.submit((__bridge void*) command);
    }
    if (packRealAsComplex && !forward) {
        unpackKernel->setArg(0, complexGrid);
        unpackKernel->setArg(1, output);
        unpackKernel->execute(complexGrid.getSize());
    }
}

int MetalFFT3D::findLegalDimension(int minimum) {
    minimum = max(1, minimum);
    while (true) {
        int unfactored = minimum;
        for (int factor = 2; factor <= 13; factor++)
            while (unfactored > 1 && unfactored%factor == 0)
                unfactored /= factor;
        if (unfactored == 1)
            return minimum;
        if (minimum == numeric_limits<int>::max())
            throw OpenMMException("No supported Metal FFT dimension fits in an integer");
        minimum++;
    }
}
