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
    id<MTLArgumentEncoder> arguments;
    struct Binding {
        bool array, local;
        size_t size;
        MTLResourceUsage usage;
    };
    vector<Binding> bindings;
};

/** Return primitive ABI size; pointers and dynamic local storage are separate. */
static size_t primitiveSize(MTLDataType type) {
    switch (type) {
        case MTLDataTypeBool: case MTLDataTypeChar: case MTLDataTypeUChar: return 1;
        case MTLDataTypeShort: case MTLDataTypeUShort: case MTLDataTypeHalf: return 2;
        case MTLDataTypeInt: case MTLDataTypeUInt: case MTLDataTypeFloat: return 4;
        case MTLDataTypeLong: case MTLDataTypeULong: return 8;
        case MTLDataTypeShort2: case MTLDataTypeUShort2: case MTLDataTypeHalf2: return 4;
        case MTLDataTypeInt2: case MTLDataTypeUInt2: case MTLDataTypeFloat2: return 8;
        case MTLDataTypeShort3: case MTLDataTypeShort4: case MTLDataTypeUShort3:
        case MTLDataTypeUShort4: case MTLDataTypeHalf3: case MTLDataTypeHalf4: return 8;
        case MTLDataTypeInt3: case MTLDataTypeInt4: case MTLDataTypeUInt3:
        case MTLDataTypeUInt4: case MTLDataTypeFloat3: case MTLDataTypeFloat4: return 16;
        case MTLDataTypeLong2: case MTLDataTypeULong2: return 16;
        case MTLDataTypeLong3: case MTLDataTypeLong4: case MTLDataTypeULong3: case MTLDataTypeULong4: return 32;
        default: throw OpenMMException("Unsupported primitive type in a Metal Common kernel");
    }
}

MetalKernel::MetalKernel(MetalContext& context, const string& name, bool commonSource,
        const function<void*(bool)>& libraryLookup) : context(context), name(name), commonSource(commonSource), libraryLookup(libraryLookup) {
    getActiveImpl();
}

MetalKernel::Impl& MetalKernel::getActiveImpl() const {
    lock_guard<mutex> lock(variantMutex);
    bool floating = commonSource && context.getUseFloatingPointAccumulators();
    int index = floating ? 1 : 0;
    if (variants[index]) return *variants[index];
    @autoreleasepool {
        unique_ptr<Impl> impl(new Impl());
        id<MTLLibrary> library = (__bridge id<MTLLibrary>) libraryLookup(floating);
        id<MTLFunction> function = [library newFunctionWithName:[NSString stringWithUTF8String:name.c_str()]];
        if (function == nil)
            throw OpenMMException("Unknown Metal kernel: "+name);
        id<MTLDevice> device = (__bridge id<MTLDevice>) context.getDevice();
        NSError* error = nil;
        impl->pipeline = [device newComputePipelineStateWithFunction:function error:&error];
        if (impl->pipeline == nil)
            throw OpenMMException("Error creating Metal pipeline "+name+": "+
                    (error == nil ? string("unknown error") : string(error.localizedDescription.UTF8String)));
        if (commonSource) {
            if (impl->pipeline.threadExecutionWidth != context.getSIMDWidth())
                throw OpenMMException("Metal Common kernels require a 32-lane SIMD group: "+name);
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wdeprecated-declarations"
            MTLArgument* reflection = nil;
            impl->arguments = [function newArgumentEncoderWithBufferIndex:0 reflection:&reflection];
#pragma clang diagnostic pop
            if (impl->arguments == nil || reflection.bufferStructType == nil)
                throw OpenMMException("Missing Common argument-buffer reflection for Metal kernel "+name);
            for (MTLStructMember* member in reflection.bufferStructType.members) {
                string memberName(member.name.UTF8String);
                if (memberName == "_metal_unused") continue;
                size_t index = member.argumentIndex;
                if (index >= impl->bindings.size()) impl->bindings.resize(index+1);
                Impl::Binding& binding = impl->bindings[index];
                binding.local = memberName.compare(0, 13, "_metal_local_") == 0;
                binding.array = member.dataType == MTLDataTypePointer;
                binding.size = binding.array ? 0 : primitiveSize(member.dataType);
                binding.usage = MTLResourceUsageRead;
                if (binding.array && member.pointerType.access != MTLArgumentAccessReadOnly)
                    binding.usage |= MTLResourceUsageWrite;
                if (binding.local && index >= 31)
                    throw OpenMMException("Metal dynamic threadgroup argument index exceeds 30 in "+name);
            }
        }
        variants[index] = std::move(impl);
    }
    return *variants[index];
}

MetalKernel::~MetalKernel() {
}

int MetalKernel::getMaxBlockSize() const {
    return getActiveImpl().pipeline.maxTotalThreadsPerThreadgroup;
}

void MetalKernel::execute(int threads, int blockSize) {
    Impl* impl = &getActiveImpl();
    if (blockSize == -1)
        blockSize = ComputeContext::ThreadBlockSize;
    if (threads < 0 || blockSize <= 0 || blockSize > impl->pipeline.maxTotalThreadsPerThreadgroup)
        throw OpenMMException("Invalid thread block size or thread count for Metal kernel "+name);
    if (threads == 0)
        return;
    if (impl->arguments == nil && arrayArgs.size() > 31)
        throw OpenMMException("Metal kernel exceeds the direct buffer argument limit");
    if (impl->arguments != nil && arrayArgs.size() != impl->bindings.size())
        throw OpenMMException("Wrong argument count for Metal Common kernel "+name+": expected "+
                to_string(impl->bindings.size())+", found "+to_string(arrayArgs.size()));
    for (int i = 0; i < arrayArgs.size(); i++) {
        // An independent context may have changed modes since this array was
        // bound. Never launch a cached pipeline against an incompatible ABI.
        if (arrayArgs[i] != nullptr)
            context.unwrap(*arrayArgs[i]);
        else if (primitiveArgSizes[i] == 0 && localArgSizes[i] == 0)
            throw OpenMMException("Unbound argument for Metal kernel "+name);
    }
    @autoreleasepool {
        id<MTLBuffer> argumentBuffer = nil;
        if (impl->arguments != nil) {
            id<MTLDevice> device = (__bridge id<MTLDevice>) context.getDevice();
            size_t localMemory = impl->pipeline.staticThreadgroupMemoryLength;
            const size_t localMemoryLimit = device.maxThreadgroupMemoryLength;
            for (size_t i = 0; i < impl->bindings.size(); i++) {
                if (!impl->bindings[i].local) continue;
                if (localArgSizes[i] == 0 || localArgSizes[i] > localMemoryLimit)
                    throw OpenMMException("Invalid dynamic threadgroup memory for Metal kernel "+name);
                size_t alignedSize = ((localArgSizes[i]+15)/16)*16;
                if (localMemory > localMemoryLimit || alignedSize > localMemoryLimit-localMemory)
                    throw OpenMMException("Total threadgroup memory exceeds the device limit for Metal kernel "+name);
                localMemory += alignedSize;
            }
            argumentBuffer = [device newBufferWithLength:impl->arguments.encodedLength options:MTLResourceStorageModeShared];
            if (argumentBuffer == nil)
                throw OpenMMException("Error allocating arguments for Metal kernel "+name);
            [impl->arguments setArgumentBuffer:argumentBuffer offset:0];
            for (size_t i = 0; i < impl->bindings.size(); i++) {
                const Impl::Binding& binding = impl->bindings[i];
                if (binding.local) {
                    // A reflected placeholder keeps the original Common argument
                    // ordinal; actual local storage is bound on the encoder.
                    continue;
                }
                else if (binding.array) {
                    id<MTLBuffer> buffer = nil;
                    if (arrayArgs[i] != nullptr)
                        buffer = (__bridge id<MTLBuffer>) arrayArgs[i]->getBuffer();
                    else {
                        // Common supplies an explicit nullptr for arrays that
                        // are unused in the selected precision/feature branch.
                        void* nullValue = nullptr;
                        if (primitiveArgSizes[i] != sizeof(nullValue) ||
                                memcmp(&primitiveArgs[i], &nullValue, sizeof(nullValue)) != 0)
                            throw OpenMMException("Expected array argument "+to_string(i)+" for Metal kernel "+name);
                    }
                    [impl->arguments setBuffer:buffer offset:0 atIndex:i];
                }
                else {
                    if (arrayArgs[i] != nullptr || primitiveArgSizes[i] != binding.size)
                        throw OpenMMException("Wrong primitive size at argument "+to_string(i)+" for Metal kernel "+name);
                    memcpy([impl->arguments constantDataAtIndex:i], &primitiveArgs[i], binding.size);
                }
            }
        }
        MetalQueue& queue = context.getCurrentMetalQueue();
        auto queueLock = queue.lock();
        id<MTLCommandBuffer> command = (__bridge id<MTLCommandBuffer>) queue.getCommandBuffer();
        id<MTLComputeCommandEncoder> encoder = [command computeCommandEncoder];
        if (encoder == nil)
            throw OpenMMException("Error creating Metal kernel command");
        [encoder setComputePipelineState:impl->pipeline];
        // Like CudaKernel, resolve arrays at each launch so resize/rebinding is visible.
        if (impl->arguments != nil) {
            [encoder setBuffer:argumentBuffer offset:0 atIndex:0];
#if !OPENMM_METAL_MINIMIZE_FLOAT_ACCUMULATORS
            [encoder setBuffer:(__bridge id<MTLBuffer>) context.getFixedPointRangeBuffer().getBuffer()
                    offset:0 atIndex:1];
#endif
            for (size_t i = 0; i < impl->bindings.size(); i++) {
                if (impl->bindings[i].local)
                    [encoder setThreadgroupMemoryLength:((localArgSizes[i]+15)/16)*16 atIndex:i];
                else if (arrayArgs[i] != nullptr)
                    [encoder useResource:(__bridge id<MTLBuffer>) arrayArgs[i]->getBuffer() usage:impl->bindings[i].usage];
            }
        }
        else {
            for (int i = 0; i < arrayArgs.size(); i++) {
                if (arrayArgs[i] != nullptr)
                    [encoder setBuffer:(__bridge id<MTLBuffer>) arrayArgs[i]->getBuffer() offset:0 atIndex:i];
                else
                    [encoder setBytes:&primitiveArgs[i] length:primitiveArgSizes[i] atIndex:i];
            }
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
    localArgSizes.push_back(0);
}

void MetalKernel::setArrayArg(int index, ArrayInterface& value) {
    ASSERT_VALID_INDEX(index, arrayArgs);
    arrayArgs[index] = &context.unwrap(value);
    localArgSizes[index] = 0;
}

void MetalKernel::setPrimitiveArg(int index, const void* value, int size) {
    ASSERT_VALID_INDEX(index, primitiveArgs);
    if (value == nullptr || size <= 0 || size > sizeof(mm_double4))
        throw OpenMMException("Unsupported value type for kernel argument");
    memcpy(&primitiveArgs[index], value, size);
    primitiveArgSizes[index] = size;
    arrayArgs[index] = nullptr;
    localArgSizes[index] = 0;
}

void MetalKernel::addLocalArg(size_t bytes) {
    int index = arrayArgs.size();
    addEmptyArg();
    setLocalArg(index, bytes);
}

void MetalKernel::setLocalArg(int index, size_t bytes) {
    Impl* impl = &getActiveImpl();
    ASSERT_VALID_INDEX(index, arrayArgs);
    if (impl->arguments == nil || size_t(index) >= impl->bindings.size() || !impl->bindings[index].local || bytes == 0)
        throw OpenMMException("Invalid local-memory argument for Metal kernel "+name);
    localArgSizes[index] = bytes;
    primitiveArgSizes[index] = 0;
    arrayArgs[index] = nullptr;
}
