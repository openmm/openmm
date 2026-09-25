/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 *                                                                            *
 * This is part of the OpenMM molecular simulation toolkit.                   *
 * See https://openmm.org/development.                                        *
 *                                                                            *
 * Ported from the OpenMM CUDA Platform.                                      *
 * Source: platforms/cuda/src/CudaKernel.cpp                                  *
 *                                                                            *
 * Original CUDA Platform code:                                               *
 * Portions copyright (c) 2019 Stanford University and the Authors.           *
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

#include "MetalKernel.h"
#include "MetalContext.h"
#include "MetalQueue.h"
#include "openmm/internal/AssertionUtilities.h"
#import <Metal/Metal.h>
#include <algorithm>
#include <cstring>

using namespace OpenMM;
using namespace std;

struct MetalKernel::Impl {
    id<MTLComputePipelineState> pipeline;
};

MetalKernel::MetalKernel(MetalContext& context, void* pipeline, const string& name) :
        impl(new Impl()), context(context), name(name) {
    impl->pipeline = (__bridge id<MTLComputePipelineState>) pipeline;
}

MetalKernel::~MetalKernel() {
}

int MetalKernel::getMaxBlockSize() const {
    return impl->pipeline.maxTotalThreadsPerThreadgroup;
}

void MetalKernel::execute(int threads, int blockSize) {
    if (blockSize == -1)
        blockSize = ComputeContext::ThreadBlockSize;
    if (threads < 0 || blockSize <= 0 || blockSize > getMaxBlockSize())
        throw OpenMMException("Invalid thread block size or thread count for Metal kernel "+name);
    if (threads == 0)
        return;
    if (arrayArgs.size() > 31)
        throw OpenMMException("Metal kernel exceeds the direct buffer argument limit");
    for (int i = 0; i < arrayArgs.size(); i++)
        if (arrayArgs[i] == nullptr && primitiveArgSizes[i] == 0)
            throw OpenMMException("Unbound argument for Metal kernel "+name);
    @autoreleasepool {
        MetalQueue& queue = context.getCurrentMetalQueue();
        id<MTLCommandQueue> commandQueue = (__bridge id<MTLCommandQueue>) queue.getQueue();
        id<MTLCommandBuffer> command = [commandQueue commandBuffer];
        id<MTLComputeCommandEncoder> encoder = [command computeCommandEncoder];
        if (encoder == nil)
            throw OpenMMException("Error creating Metal kernel command");
        [encoder setComputePipelineState:impl->pipeline];
        // Like CudaKernel, resolve arrays at each launch so resize/rebinding is visible.
        for (int i = 0; i < arrayArgs.size(); i++) {
            if (arrayArgs[i] != nullptr)
                [encoder setBuffer:(__bridge id<MTLBuffer>) arrayArgs[i]->getBuffer() offset:0 atIndex:i];
            else
                [encoder setBytes:&primitiveArgs[i] length:primitiveArgSizes[i] atIndex:i];
        }
        int gridSize = min(1+(threads-1)/blockSize, context.getNumThreadBlocks());
        [encoder dispatchThreadgroups:MTLSizeMake(gridSize, 1, 1)
                threadsPerThreadgroup:MTLSizeMake(blockSize, 1, 1)];
        [encoder endEncoding];
        queue.submit((__bridge void*) command);
    }
}

void MetalKernel::addArrayArg(ArrayInterface& value) {
    int index = arrayArgs.size();
    addEmptyArg();
    setArrayArg(index, value);
}

void MetalKernel::addPrimitiveArg(const void* value, int size) {
    int index = arrayArgs.size();
    addEmptyArg();
    setPrimitiveArg(index, value, size);
}

void MetalKernel::addEmptyArg() {
    primitiveArgs.push_back(mm_double4(0, 0, 0, 0));
    primitiveArgSizes.push_back(0);
    arrayArgs.push_back(nullptr);
}

void MetalKernel::setArrayArg(int index, ArrayInterface& value) {
    ASSERT_VALID_INDEX(index, arrayArgs);
    arrayArgs[index] = &context.unwrap(value);
}

void MetalKernel::setPrimitiveArg(int index, const void* value, int size) {
    ASSERT_VALID_INDEX(index, primitiveArgs);
    if (value == nullptr || size <= 0 || size > sizeof(mm_double4))
        throw OpenMMException("Unsupported value type for kernel argument");
    memcpy(&primitiveArgs[index], value, size);
    primitiveArgSizes[index] = size;
    arrayArgs[index] = nullptr;
}
