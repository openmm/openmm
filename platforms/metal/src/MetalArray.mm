/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 *                                                                            *
 * This is part of the OpenMM molecular simulation toolkit.                   *
 * See https://openmm.org/development.                                        *
 *                                                                            *
 * Ported from the OpenMM CUDA Platform.                                      *
 * Source: platforms/cuda/src/CudaArray.cpp                                   *
 *                                                                            *
 * Original CUDA Platform code:                                               *
 * Portions copyright (c) 2012-2022 Stanford University and the Authors.      *
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

#include "MetalArray.h"
#include "MetalContext.h"
#include "MetalQueue.h"
#include "openmm/OpenMMException.h"
#import <Metal/Metal.h>
#include <cstring>
#include <cstdint>
#include <limits>

using namespace OpenMM;
using namespace std;

struct MetalArray::Impl {
    id<MTLBuffer> buffer = nil;
};

static size_t getPinnedOffset(id<MTLBuffer> buffer, const void* data, size_t bytes) {
    uintptr_t start = reinterpret_cast<uintptr_t>(buffer.contents);
    uintptr_t address = reinterpret_cast<uintptr_t>(data);
    if (address < start || address-start > buffer.length || bytes > buffer.length-(address-start))
        throw OpenMMException("Nonblocking Metal transfers require memory from getPinnedBuffer()");
    return address-start;
}

MetalArray::MetalArray() : impl(new Impl()), context(nullptr), size(0), elementSize(0) {
}

MetalArray::MetalArray(MetalContext& context, size_t size, int elementSize, const string& name) : MetalArray() {
    initialize(context, size, elementSize, name);
}

MetalArray::~MetalArray() {
}

void MetalArray::initialize(ComputeContext& context, size_t size, int elementSize, const string& name) {
    if (isInitialized())
        throw OpenMMException("MetalArray has already been initialized");
    MetalContext& metal = dynamic_cast<MetalContext&>(context);
    if (elementSize <= 0 || size > numeric_limits<int>::max() ||
            size > numeric_limits<size_t>::max()/elementSize)
        throw OpenMMException("Invalid MetalArray size: "+name);
    @autoreleasepool {
        id<MTLDevice> device = (__bridge id<MTLDevice>) metal.getDevice();
        // A zero-length logical array still needs a legal Metal resource.
        impl->buffer = [device newBufferWithLength:max(size*elementSize, size_t(1)) options:MTLResourceStorageModePrivate];
        if (impl->buffer == nil)
            throw OpenMMException("Error creating Metal array "+name);
        impl->buffer.label = [NSString stringWithUTF8String:name.c_str()];
    }
    this->context = &metal;
    this->size = size;
    this->elementSize = elementSize;
    this->name = name;
}

void MetalArray::resize(size_t size) {
    if (!isInitialized())
        throw OpenMMException("MetalArray has not been initialized");
    MetalArray replacement(*context, size, elementSize, name);
    impl.swap(replacement.impl);
    this->size = size;
}

bool MetalArray::isInitialized() const {
    return impl->buffer != nil;
}

ComputeContext& MetalArray::getContext() {
    if (!isInitialized())
        throw OpenMMException("MetalArray has not been initialized");
    return *context;
}

void* MetalArray::getBuffer() const {
    if (!isInitialized())
        throw OpenMMException("MetalArray has not been initialized");
    return (__bridge void*) impl->buffer;
}

void MetalArray::uploadSubArray(const void* data, int offset, int elements, bool blocking) {
    if (!isInitialized())
        throw OpenMMException("MetalArray has not been initialized");
    if (offset < 0 || elements < 0 || size_t(offset) > size || size_t(elements) > size-offset)
        throw OpenMMException("uploadSubArray: data exceeds range of array");
    if (elements == 0)
        return;
    if (data == nullptr)
        throw OpenMMException("Null source for Metal array upload");
    @autoreleasepool {
        MetalQueue& queue = context->getCurrentMetalQueue();
        id<MTLDevice> device = (__bridge id<MTLDevice>) context->getDevice();
        size_t bytes = size_t(elements)*elementSize;
        id<MTLBuffer> staging;
        size_t sourceOffset = 0;
        if (blocking)
            staging = [device newBufferWithBytes:data length:bytes options:MTLResourceStorageModeShared];
        else {
            staging = (__bridge id<MTLBuffer>) context->getPinnedBufferHandle();
            sourceOffset = getPinnedOffset(staging, data, bytes);
        }
        if (staging == nil)
            throw OpenMMException("Error creating upload buffer for "+name);
        id<MTLCommandQueue> commandQueue = (__bridge id<MTLCommandQueue>) queue.getQueue();
        id<MTLCommandBuffer> command = [commandQueue commandBuffer];
        id<MTLBlitCommandEncoder> encoder = [command blitCommandEncoder];
        if (encoder == nil)
            throw OpenMMException("Error creating Metal upload command");
        [encoder copyFromBuffer:staging sourceOffset:sourceOffset toBuffer:impl->buffer
                destinationOffset:size_t(offset)*elementSize size:bytes];
        [encoder endEncoding];
        queue.submit((__bridge void*) command);
        if (blocking)
            queue.wait((__bridge void*) command);
    }
}

void MetalArray::download(void* data, bool blocking) const {
    if (!isInitialized())
        throw OpenMMException("MetalArray has not been initialized");
    if (size == 0)
        return;
    if (data == nullptr)
        throw OpenMMException("Null destination for Metal array download");
    @autoreleasepool {
        MetalQueue& queue = context->getCurrentMetalQueue();
        id<MTLDevice> device = (__bridge id<MTLDevice>) context->getDevice();
        size_t bytes = size*elementSize;
        id<MTLBuffer> staging;
        size_t destinationOffset = 0;
        if (blocking)
            staging = [device newBufferWithLength:bytes options:MTLResourceStorageModeShared];
        else {
            staging = (__bridge id<MTLBuffer>) context->getPinnedBufferHandle();
            destinationOffset = getPinnedOffset(staging, data, bytes);
        }
        if (staging == nil)
            throw OpenMMException("Error creating download buffer for "+name);
        id<MTLCommandQueue> commandQueue = (__bridge id<MTLCommandQueue>) queue.getQueue();
        id<MTLCommandBuffer> command = [commandQueue commandBuffer];
        id<MTLBlitCommandEncoder> encoder = [command blitCommandEncoder];
        if (encoder == nil)
            throw OpenMMException("Error creating Metal download command");
        [encoder copyFromBuffer:impl->buffer sourceOffset:0 toBuffer:staging destinationOffset:destinationOffset size:bytes];
        [encoder endEncoding];
        queue.submit((__bridge void*) command);
        if (blocking) {
            queue.wait((__bridge void*) command);
            memcpy(data, staging.contents, bytes);
        }
    }
}

void MetalArray::copyTo(ArrayInterface& dest) const {
    if (!isInitialized())
        throw OpenMMException("MetalArray has not been initialized");
    if (dest.getSize() != size || dest.getElementSize() != elementSize)
        throw OpenMMException("Error copying array "+name+": destination size does not match");
    MetalArray& metalDest = context->unwrap(dest);
    if (size == 0 || &metalDest == this)
        return;
    @autoreleasepool {
        MetalQueue& queue = context->getCurrentMetalQueue();
        id<MTLCommandQueue> commandQueue = (__bridge id<MTLCommandQueue>) queue.getQueue();
        id<MTLCommandBuffer> command = [commandQueue commandBuffer];
        id<MTLBlitCommandEncoder> encoder = [command blitCommandEncoder];
        if (encoder == nil)
            throw OpenMMException("Error creating Metal copy command");
        [encoder copyFromBuffer:impl->buffer sourceOffset:0 toBuffer:metalDest.impl->buffer
                destinationOffset:0 size:size*elementSize];
        [encoder endEncoding];
        queue.submit((__bridge void*) command);
    }
}
