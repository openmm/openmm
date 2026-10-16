/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit.                   *
 * See https://openmm.org/development.                                        *
 *                                                                            *
 * Portions copyright (c) 2019-2026 Stanford University and the Authors.      *
 * Authors: Peter Eastman                                                     *
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

#include "MetalKernel.h"
#include "MetalQueue.h"
#include "openmm/common/ComputeArray.h"
#include "openmm/internal/AssertionUtilities.h"
#include <cstring>
#include <vector>

using namespace OpenMM;
using namespace std;

MetalKernel::MetalKernel(MetalContext& context, MTL::ComputePipelineState* pipeline, const string& name, bool packArgs) :
        context(context), pipeline(pipeline), name(name), dynamicLocal(0), packArgs(packArgs) {
}

MetalKernel::~MetalKernel() {
    pipeline->release();
}

string MetalKernel::getName() const {
    return name;
}

int MetalKernel::getMaxBlockSize() const {
    return pipeline->maxTotalThreadsPerThreadgroup();
}

void MetalKernel::execute(int threads, int blockSize) {
    MetalQueue& queue = dynamic_cast<MetalQueue&>(*context.getCurrentQueue());
    MTL::ComputeCommandEncoder& encoder = queue.getEncoder();
    encoder.setComputePipelineState(pipeline);
    int numArgs = arrayArgs.size();
    if (packArgs) {
        // There are too many arguments to pass individually, so we need to pack all of them into a single argument.

        vector<mm_double4> args = primitiveArgs;
        for (int i = 0; i < numArgs; i++) {
            if (arrayArgs[i] != NULL) {
                MTL::GPUAddress address = arrayArgs[i]->getBuffer()->gpuAddress();
                memcpy(&args[i], &address, sizeof(address));
                encoder.useResource(arrayArgs[i]->getBuffer(), MTL::ResourceUsageRead | MTL::ResourceUsageWrite);
            }
        }
        encoder.setBytes(args.data(), numArgs*sizeof(mm_double4), 0);
    }
    else {
        for (int i = 0; i < numArgs; i++) {
            if (arrayArgs[i] != NULL)
                encoder.setBuffer(arrayArgs[i]->getBuffer(), 0, i);
            else
                encoder.setBytes(&primitiveArgs[i], primitiveArgSizes[i], i);
        }
    }
    if (dynamicLocal > 0)
        encoder.setThreadgroupMemoryLength(dynamicLocal, 0);
    if (blockSize == -1)
        blockSize = MetalContext::ThreadBlockSize;
    int gridSize = min((threads+blockSize-1)/blockSize, context.getNumThreadBlocks());
    encoder.dispatchThreadgroups(MTL::Size(gridSize, 1, 1), MTL::Size(blockSize, 1, 1));
}

void MetalKernel::setDynamicLocalMemory(int bytes) {
    // Round up to the next multiple of 16.
    dynamicLocal = 16*((bytes+15)/16);
}

void MetalKernel::addArrayArg(ArrayInterface& value) {
    int index = arrayArgs.size();
    addEmptyArg();
    setArrayArg(index, value);
}

void MetalKernel::addPrimitiveArg(const void* value, int size, const type_info& type) {
    int index = arrayArgs.size();
    addEmptyArg();
    setPrimitiveArg(index, value, size, type);
}

void MetalKernel::addEmptyArg() {
    primitiveArgs.push_back(mm_double4(0, 0, 0, 0));
    primitiveArgSizes.push_back(0);
    arrayArgs.push_back(NULL);
}

void MetalKernel::setArrayArg(int index, ArrayInterface& value) {
    ASSERT_VALID_INDEX(index, arrayArgs);
    arrayArgs[index] = &context.unwrap(value);
}

void MetalKernel::setPrimitiveArg(int index, const void* value, int size, const type_info& type) {
    ASSERT_VALID_INDEX(index, primitiveArgs);
    if (size > sizeof(mm_double4))
        throw OpenMMException("Unsupported value type for kernel argument");
    if (type == typeid(double) || type == typeid(mm_double2) || type == typeid(mm_double4)) {
        // Convert double to df64.

        const double* src = (const double*) value;
        float* dest = (float*) &primitiveArgs[index];
        int elements = size/sizeof(double);
        for (int i = 0; i < elements; i++) {
            float f = (float) src[i];
            dest[2*i] = f;
            dest[2*i+1] = (float) (src[i]-f);
        }
    }
    else {
        memcpy(&primitiveArgs[index], value, size);
    }
    primitiveArgSizes[index] = size;
    arrayArgs[index] = NULL;
}
