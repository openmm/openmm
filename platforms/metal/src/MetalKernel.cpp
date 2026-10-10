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

MetalKernel::MetalKernel(MetalContext& context, MTL::ComputePipelineState* pipeline, const string& name) :
        context(context), pipeline(pipeline), name(name), dynamicLocal(0) {
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
    argPointers.resize(numArgs);
    for (int i = 0; i < numArgs; i++) {
        if (arrayArgs[i] != NULL)
            encoder.setBuffer(arrayArgs[i]->getBuffer(), 0, i);
        else
            encoder.setBytes(&primitiveArgs[i], primitiveArgSizes[i], i);
    }
    if (dynamicLocal > 0)
        encoder.setThreadgroupMemoryLength(dynamicLocal, 0);
    if (blockSize == -1)
        blockSize = MetalContext::ThreadBlockSize;
    int gridSize = min((threads+blockSize-1)/blockSize, context.getNumThreadBlocks());
    encoder.dispatchThreadgroups(MTL::Size(gridSize, 1, 1), MTL::Size(blockSize, 1, 1));
}

void MetalKernel::setDynamicLocalMemory(int bytes) {
    dynamicLocal = bytes;
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
    arrayArgs.push_back(NULL);
}

void MetalKernel::setArrayArg(int index, ArrayInterface& value) {
    ASSERT_VALID_INDEX(index, arrayArgs);
    arrayArgs[index] = &context.unwrap(value);
}

void MetalKernel::setPrimitiveArg(int index, const void* value, int size) {
    ASSERT_VALID_INDEX(index, primitiveArgs);
    if (size > sizeof(mm_double4))
        throw OpenMMException("Unsupported value type for kernel argument");
    memcpy(&primitiveArgs[index], value, size);
    primitiveArgSizes[index] = size;
    arrayArgs[index] = NULL;
}
