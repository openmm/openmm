/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2014-2015 Stanford University and the Authors.      *
 * Authors: Robert T. McGibbon                                                *
 * Contributors: Daniel Towner                                                *
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
 * This tests vectorized operations.
 */

#include "openmm/internal/AssertionUtilities.h"
#include "openmm/internal/vectorizeAvx.h"
#include <iostream>

#include "TestVectorizeGeneric.h"

#ifndef __AVX__
bool isVec8Supported() {
    return false;
}
#else
/**
 * Check whether 8 component vectors are supported with the current CPU.
 */
bool isVec8Supported() {
    // Make sure the CPU supports AVX.

    int cpuInfo[4];
    cpuid(cpuInfo, 0);
    if (cpuInfo[0] >= 1) {
        cpuid(cpuInfo, 1);
        return ((cpuInfo[2] & ((int) 1 << 28)) != 0);
    }
    return false;
}
#endif

using namespace OpenMM;
using namespace std;

#define ASSERT_VEC4_EQUAL(found, expected0, expected1, expected2, expected3) {if (std::abs((found)[0]-(expected0))>1e-6 || std::abs((found)[1]-(expected1))>1e-6 || std::abs((found)[2]-(expected2))>1e-6 || std::abs((found)[3]-(expected3))>1e-6) {std::stringstream details; details << " Expected ("<<(expected0)<<","<<(expected1)<<","<<(expected2)<<","<<(expected3)<<"), found ("<<(found)[0]<<","<<(found)[1]<<","<<(found)[2]<<","<<(found)[3]<<")"; throwException(__FILE__, __LINE__, details.str());}};
#define ASSERT_VEC8_EQUAL(found, expected0, expected1, expected2, expected3, expected4, expected5, expected6, expected7) {if (std::abs((found).lowerVec()[0]-(expected0))>1e-6 || std::abs((found).lowerVec()[1]-(expected1))>1e-6 || std::abs((found).lowerVec()[2]-(expected2))>1e-6 || std::abs((found).lowerVec()[3]-(expected3))>1e-6 || std::abs((found).upperVec()[0]-(expected4))>1e-6 || std::abs((found).upperVec()[1]-(expected5))>1e-6 || std::abs((found).upperVec()[2]-(expected6))>1e-6 || std::abs((found).upperVec()[3]-(expected7))>1e-6) {std::stringstream details; details << " Expected ("<<(expected0)<<","<<(expected1)<<","<<(expected2)<<","<<(expected3)<<","<<(expected4)<<","<<(expected5)<<","<<(expected6)<<","<<(expected7)<<"), found ("<<(found).lowerVec()[0]<<","<<(found).lowerVec()[1]<<","<<(found).lowerVec()[2]<<","<<(found).lowerVec()[3]<<","<<(found).upperVec()[0]<<","<<(found).upperVec()[1]<<","<<(found).upperVec()[2]<<","<<(found).upperVec()[3]<<")"; throwException(__FILE__, __LINE__, details.str());}};
#define ASSERT_VEC8_EQUAL_INT(found, expected0, expected1, expected2, expected3, expected4, expected5, expected6, expected7) {if ((found).lowerVec()[0] != (expected0) || (found).lowerVec()[1] != (expected1) || (found).lowerVec()[2] != (expected2) || (found).lowerVec()[3] != (expected3) || (found).upperVec()[0] != (expected4) || (found).upperVec()[1] != (expected5) ||(found).upperVec()[2] != (expected6) || (found).upperVec()[3] != (expected7)) {std::stringstream details; details << " Expected ("<<(expected0)<<","<<(expected1)<<","<<(expected2)<<","<<(expected3)<<","<<(expected4)<<","<<(expected5)<<","<<(expected6)<<","<<(expected7)<<"), found ("<<(found).lowerVec()[0]<<","<<(found).lowerVec()[1]<<","<<(found).lowerVec()[2]<<","<<(found).lowerVec()[3]<<","<<(found).upperVec()[0]<<","<<(found).upperVec()[1]<<","<<(found).upperVec()[2]<<","<<(found).upperVec()[3]<<")"; throwException(__FILE__, __LINE__, details.str());}};

void testLoadStore() {
    ivec8 i1(3);
    ASSERT_VEC8_EQUAL_INT(i1, 3, 3, 3, 3, 3, 3, 3, 3);
    ivec8 i2(2, 3, 4, 5, 6, 7, 8, 9);
    ASSERT_VEC8_EQUAL_INT(i2, 2, 3, 4, 5, 6, 7, 8, 9);
    int iarray[8];
    i2.store(iarray);
    ivec8 i3(iarray);
    ASSERT_VEC8_EQUAL_INT(i3, 2, 3, 4, 5, 6, 7, 8, 9);
    ASSERT_EQUAL(i3.lowerVec()[0], 2);
    ASSERT_EQUAL(i3.lowerVec()[1], 3);
    ASSERT_EQUAL(i3.lowerVec()[2], 4);
    ASSERT_EQUAL(i3.lowerVec()[3], 5);
    ASSERT_EQUAL(i3.upperVec()[0], 6);
    ASSERT_EQUAL(i3.upperVec()[1], 7);
    ASSERT_EQUAL(i3.upperVec()[2], 8);
    ASSERT_EQUAL(i3.upperVec()[3], 9);

    // Partial store of vec3 should not overwrite beyond the 3 elements.
    // Note that this is a fvec4 method, but is conditionally compiled for AVX so needs to be
    // tested here too.
    float overwriteTest[4] = {9, 9, 9, 9};
    fvec4(1, 2, 3, 7777).storeVec3(overwriteTest);
    ASSERT_EQUAL(overwriteTest[0], 1);
    ASSERT_EQUAL(overwriteTest[1], 2);
    ASSERT_EQUAL(overwriteTest[2], 3);
    ASSERT_EQUAL(overwriteTest[3], 9);
}

void testLogic() {
    int allBits = -1;
    float allBitsf = *((float*) &allBits);
    ivec8 mask(0, allBits, allBits, 0, 0, allBits, allBits, 0);
    fvec8 fmask(0, allBitsf, allBitsf, 0, 0, allBitsf, allBitsf, 0);
    fvec8 f1(0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0);
    ivec8 i1(1, 2, 3, 4, 5, 6, 7, 8);
    ASSERT_VEC8_EQUAL(f1&fmask, 0, 1.0, 1.5, 0, 0, 3.0, 3.5, 0.0);
    fvec8 temp = f1|fmask;
    ASSERT_EQUAL(0.5, temp.lowerVec()[0]);
    ASSERT(temp.lowerVec()[1]!= temp.lowerVec()[1]); // All bits set, which is nan
    ASSERT(temp.lowerVec()[2] != temp.lowerVec()[2]); // All bits set, which is nan
    ASSERT_EQUAL(2.0, temp.lowerVec()[3]);
    ASSERT_EQUAL(2.5, temp.upperVec()[0]);
    ASSERT(temp.upperVec()[1] != temp.upperVec()[1]); // All bits set, which is nan
    ASSERT(temp.upperVec()[2] != temp.upperVec()[2]); // All bits set, which is nan
    ASSERT_EQUAL(4.0, temp.upperVec()[3]);
    ASSERT_VEC8_EQUAL_INT(i1&mask, 0, 2, 3, 0, 0, 6, 7, 0);
    ASSERT_VEC8_EQUAL_INT(i1|mask, 1, allBits, allBits, 4, 5, allBits, allBits, 8);
}

int main(int argc, char* argv[]) {
    try {
        if (!isVec8Supported()) {
            cout << "CPU is not supported.  Exiting." << endl;
            return 0;
        }
        testLoadStore();
        testLogic();

        TestFvec<fvec8>::testAll();

    }
    catch(const exception& e) {
        cout << "exception: " << e.what() << endl;
        return 1;
    }
    cout << "Done" << endl;
    return 0;
}
