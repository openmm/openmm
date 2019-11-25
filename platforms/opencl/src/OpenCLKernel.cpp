/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2019 Stanford University and the Authors.           *
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

#include "OpenCLKernel.h"
#include "ComputeArray.h"

using namespace OpenMM;
using namespace std;

OpenCLKernel::OpenCLKernel(OpenCLContext& context, cl::Kernel kernel) : context(context), kernel(kernel) {
}

const string& OpenCLKernel::getName() const {
    return kernel.getInfo<CL_KERNEL_FUNCTION_NAME>();
}

void OpenCLKernel::execute(int threads, int blockSize) {
    // Set args that are specified by OpenCLArrays.  We can't do this earlier, because it's
    // possible resize() will get called on an array, causing its internal storage to be
    // recreated.
    
    for (int i = 0; i < arrayArgs.size(); i++)
        if (arrayArgs[i] != NULL)
            kernel.setArg<cl::Buffer>(i, arrayArgs[i]->getDeviceBuffer());
    context.executeKernel(kernel, threads, blockSize);
}

void OpenCLKernel::addArrayArg(ArrayInterface& value) {
    int index = arrayArgs.size();
    arrayArgs.push_back(NULL);
    setArrayArg(index, value);
}

void OpenCLKernel::addPrimitiveArg(void* value, int size) {
    int index = arrayArgs.size();
    arrayArgs.push_back(NULL);
    setPrimitiveArg(index, value, size);
}

void OpenCLKernel::setArrayArg(int index, ArrayInterface& value) {
    arrayArgs[index] = &context.unwrap(value);
}

void OpenCLKernel::setPrimitiveArg(int index, void* value, int size) {
    kernel.setArg(index, size, value);
}
