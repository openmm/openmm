/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2019 Stanford University and the Authors.           *
 * Portions copyright (c) 2020 Advanced Micro Devices, Inc.                   *
 * Authors: Peter Eastman, Nicholas Curtis                                    *
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

#include "HipKernel.h"
#include "openmm/common/ComputeArray.h"
#include "openmm/internal/AssertionUtilities.h"
#include <cstring>
#include <vector>

using namespace OpenMM;
using namespace std;

HipKernel::HipKernel(HipContext& context, hipFunction_t kernel, const string& name) : context(context), kernel(kernel), name(name) {
}

string HipKernel::getName() const {
    return name;
}

int HipKernel::getMaxBlockSize() const {
    int size;
    hipError_t result = hipFuncGetAttribute(&size, HIP_FUNC_ATTRIBUTE_MAX_THREADS_PER_BLOCK, kernel);
    if (result != hipSuccess)
        throw OpenMMException("Error querying max thread block size: "+context.getErrorString(result));
    return size;
}

void HipKernel::execute(int threads, int blockSize) {
    int numArgs = arrayArgs.size();
    argPointers.resize(numArgs);
    for (int i = 0; i < numArgs; i++) {
        if (arrayArgs[i] != NULL)
            argPointers[i] = &arrayArgs[i]->getDevicePointer();
        else
            argPointers[i] = &primitiveArgs[i];
    }
    context.executeKernel(kernel, argPointers.data(), threads, blockSize);
}

void HipKernel::addArrayArg(ArrayInterface& value) {
    int index = arrayArgs.size();
    addEmptyArg();
    setArrayArg(index, value);
}

void HipKernel::addPrimitiveArg(const void* value, int size) {
    int index = arrayArgs.size();
    addEmptyArg();
    setPrimitiveArg(index, value, size);
}

void HipKernel::addEmptyArg() {
    primitiveArgs.push_back(make_double4(0, 0, 0, 0));
    arrayArgs.push_back(NULL);
}

void HipKernel::setArrayArg(int index, ArrayInterface& value) {
    ASSERT_VALID_INDEX(index, arrayArgs);
    arrayArgs[index] = &context.unwrap(value);
}

void HipKernel::setPrimitiveArg(int index, const void* value, int size) {
    ASSERT_VALID_INDEX(index, primitiveArgs);
    if (size > sizeof(double4))
        throw OpenMMException("Unsupported value type for kernel argument");
    memcpy(&primitiveArgs[index], value, size);
    arrayArgs[index] = NULL;
}
