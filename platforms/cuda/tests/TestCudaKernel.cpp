/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit.                   *
 * See https://openmm.org/development.                                        *
 *                                                                            *
 * Portions copyright (c) 2008-2025 Stanford University and the Authors.      *
 * Contributors:                                                              *
 *                                                                            *
 * Permission is hereby granted, free of charge, to any person obtaining a    *
 * copy of this software and associated documentation files (the "Software"), *
 * to deal in the Software without restriction, including without limitation  *
 * the rights to use, copy, modify, merge, publish, distribute, sublicense,   *
 * and/or sell copies of the Software, and to permit persons to whom the      *
 * Software is furnished to do so, subject to the following conditions:       *
 *                                                                            *
 * The above copyright notice and this permission notice shall be included in *
 * all copies or substantial portions of the Software.                        *
 *                                                                            *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR *
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,   *
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL    *
 * THE AUTHORS, CONTRIBUTORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,    *
 * DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR      *
 * OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE  *
 * USE OR OTHER DEALINGS IN THE SOFTWARE.                                     *
 * -------------------------------------------------------------------------- */

/**
 * This tests CUDA kernel argument updates between executions.
 */

#include "openmm/internal/AssertionUtilities.h"
#include "CudaArray.h"
#include "CudaContext.h"
#include "openmm/System.h"
#include <iostream>
#include <utility>
#include <vector>

using namespace OpenMM;
using namespace std;

CudaPlatform platform;

void verifyResult(ComputeKernel kernel, CudaArray& output, int expected) {
    kernel->execute(1);
    vector<int> result;
    output.download(result);
    ASSERT_EQUAL(expected, result[0]);
}

void testArgumentUpdates() {
    System system;
    system.addParticle(0.0);
    CudaPlatform::PlatformData platformData(NULL, system, "", "true", platform.getPropertyDefaultValue("CudaPrecision"), "false",
            platform.getPropertyDefaultValue(CudaPlatform::CudaTempDirectory()),
            platform.getPropertyDefaultValue(CudaPlatform::CudaDisablePmeStream()), "false", 1, NULL);
    CudaContext& context = *platformData.contexts[0];
    context.initialize();
    context.setAsCurrent();

    ComputeProgram program = context.compileProgram(
            "extern \"C\" __global__ void computeValue(const int* input, int index, int offset, int* output) {\n"
            "    if (blockIdx.x == 0 && threadIdx.x == 0)\n"
            "        output[0] = input[index]+offset;\n"
            "}\n");
    ComputeKernel kernel = program->createKernel("computeValue");
    CudaArray first(context, 3, sizeof(int), "first");
    CudaArray second(context, 3, sizeof(int), "second");
    CudaArray output(context, 1, sizeof(int), "output");
    vector<int> firstValues = {7, 19, 41};
    vector<int> secondValues = {-13, 29, 53};
    first.upload(firstValues);
    second.upload(secondValues);
    kernel->addArg(first);
    kernel->addArg();
    kernel->addArg(0);
    kernel->addArg(output);
    kernel->setArg(1, 0);
    verifyResult(kernel, output, firstValues[0]);

    // Updating primitive values must take effect on every execution.

    for (int i = 0; i < 6; i++) {
        kernel->setArg(1, i%3);
        kernel->setArg(2, 10*i);
        verifyResult(kernel, output, firstValues[i%3]+10*i);
    }

    // Rebinding either the same array or a different one must preserve updates.

    firstValues[2] = 71;
    first.upload(firstValues);
    kernel->setArg(0, first);
    verifyResult(kernel, output, firstValues[2]+50);
    kernel->setArg(0, second);
    verifyResult(kernel, output, secondValues[2]+50);

    // A device pointer supplied as a primitive has the same kernel ABI as an
    // array argument, but its argument storage is different.  Exercise both
    // transitions, as well as updates while it remains a primitive argument.

    kernel->setArg(0, first.getDevicePointer());
    verifyResult(kernel, output, firstValues[2]+50);
    kernel->setArg(0, second.getDevicePointer());
    verifyResult(kernel, output, secondValues[2]+50);
    kernel->setArg(0, first);
    verifyResult(kernel, output, firstValues[2]+50);

    // Resizing changes the storage of an already bound array.  No setArg()
    // call should be needed for the kernel to use its current device pointer.

    first.resize(257);
    firstValues.assign(257, 101);
    first.upload(firstValues);
    verifyResult(kernel, output, firstValues[2]+50);

    // The allocator might reuse an address during resize().  Swapping the
    // pointers of two equally sized arrays guarantees a changed device address
    // while leaving the array objects, sizes, and memory ownership valid.

    second.resize(257);
    secondValues.assign(257, -37);
    second.upload(secondValues);
    swap(first.getDevicePointer(), second.getDevicePointer());
    verifyResult(kernel, output, secondValues[2]+50);
}

int main(int argc, char* argv[]) {
    try {
        if (argc > 1)
            platform.setPropertyDefaultValue("CudaPrecision", string(argv[1]));
        testArgumentUpdates();
    }
    catch(const exception& e) {
        cout << "exception: " << e.what() << endl;
        return 1;
    }
    cout << "Done" << endl;
    return 0;
}
