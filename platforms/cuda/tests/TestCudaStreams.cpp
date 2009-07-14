/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2008 Stanford University and the Authors.           *
 * Authors: Peter Eastman                                                     *
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
 * This tests the CUDA implementation of streams.
 */

#include "../../../tests/AssertionUtilities.h"
#include "CudaPlatform.h"
#include "openmm/Context.h"
#include "openmm/Stream.h"
#include "openmm/VerletIntegrator.h"
#include <iostream>

using namespace OpenMM;
using namespace std;

const double TOL = 1e-6;

template <class T, int WIDTH>
void testStream(Stream::DataType type, T scale) {
    const int size = 100;
    CudaPlatform platform;
    System system;
    for (int i = 0; i < size; i++)
        system.addParticle(1.0);
    VerletIntegrator integrator(0.01);
    Context context(system, integrator, platform);
    ContextImpl* impl = *reinterpret_cast<ContextImpl**>(&context);
    Stream stream = platform.createStream("", size, type, *impl);
    const int length = size*WIDTH;
    T array[size*WIDTH+1];
    array[length] = 0;
    T value = 3;
    
    // Test fillWithValue.
    
    stream.fillWithValue(&value);
    stream.saveToArray(array);
    for (int i = 0; i < length; ++i)
        ASSERT_EQUAL(3, array[i]);
    ASSERT_EQUAL(0, array[length]);
    
    // Test loadFromArray.
    
    for (int i = 0; i < length; ++i)
        array[i] = i*scale;
    stream.loadFromArray(array);
    for (int i = 0; i < length; ++i)
        array[i] = 0;
    stream.saveToArray(array);
    for (int i = 0; i < length; ++i)
        ASSERT_EQUAL_TOL((double) (i*scale), array[i], TOL);
    ASSERT_EQUAL_TOL(0.0, (double) array[length], TOL);
}

int main() {
    try {
        testStream<float, 1>(Stream::Float, 0.1f);
        testStream<float, 2>(Stream::Float2, 0.1f);
        testStream<float, 3>(Stream::Float3, 0.1f);
        testStream<float, 4>(Stream::Float4, 0.1f);
        testStream<double, 1>(Stream::Double, 0.1);
        testStream<double, 2>(Stream::Double2, 0.1);
        testStream<double, 3>(Stream::Double3, 0.1);
        testStream<double, 4>(Stream::Double4, 0.1);
        testStream<int, 1>(Stream::Integer, 1);
        testStream<int, 2>(Stream::Integer2, 1);
        testStream<int, 3>(Stream::Integer3, 1);
        testStream<int, 4>(Stream::Integer4, 1);
    }
    catch(const exception& e) {
        cout << "exception: " << e.what() << endl;
        return 1;
    }
    cout << "Done" << endl;
    return 0;
}
