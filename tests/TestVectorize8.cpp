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
 * This tests vectorized operations.
 */

#include "openmm/internal/AssertionUtilities.h"
#include "openmm/internal/vectorize8.h"
#include <iostream>


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
    fvec8 f1(2.0);
    ivec8 i1(3);
    ASSERT_VEC8_EQUAL(f1, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0);
    ASSERT_VEC8_EQUAL_INT(i1, 3, 3, 3, 3, 3, 3, 3, 3);
    fvec8 f2(2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 5.5, 6.0);
    ivec8 i2(2, 3, 4, 5, 6, 7, 8, 9);
    ASSERT_VEC8_EQUAL(f2, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 5.5, 6.0);
    ASSERT_VEC8_EQUAL_INT(i2, 2, 3, 4, 5, 6, 7, 8, 9);
    float farray[8];
    int iarray[8];
    f2.store(farray);
    i2.store(iarray);
    fvec8 f3(farray);
    ivec8 i3(iarray);
    ASSERT_VEC8_EQUAL(f3, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 5.5, 6.0);
    ASSERT_VEC8_EQUAL_INT(i3, 2, 3, 4, 5, 6, 7, 8, 9);
    ASSERT_EQUAL(f3.lowerVec()[0], 2.5);
    ASSERT_EQUAL(f3.lowerVec()[1], 3.0);
    ASSERT_EQUAL(f3.lowerVec()[2], 3.5);
    ASSERT_EQUAL(f3.lowerVec()[3], 4.0);
    ASSERT_EQUAL(f3.upperVec()[0], 4.5);
    ASSERT_EQUAL(f3.upperVec()[1], 5.0);
    ASSERT_EQUAL(f3.upperVec()[2], 5.5);
    ASSERT_EQUAL(f3.upperVec()[3], 6.0);
    ASSERT_EQUAL(i3.lowerVec()[0], 2);
    ASSERT_EQUAL(i3.lowerVec()[1], 3);
    ASSERT_EQUAL(i3.lowerVec()[2], 4);
    ASSERT_EQUAL(i3.lowerVec()[3], 5);
    ASSERT_EQUAL(i3.upperVec()[0], 6);
    ASSERT_EQUAL(i3.upperVec()[1], 7);
    ASSERT_EQUAL(i3.upperVec()[2], 8);
    ASSERT_EQUAL(i3.upperVec()[3], 9);
}

void testArithmetic() {
    fvec8 f1(0.5, 1.0, 1.5, 2.0,   2.5, 3.0, 3.5, 4.0);
    ASSERT_VEC8_EQUAL(f1+fvec8(1, 2, 3, 4, 5, 6, 7, 8), 1.5,   3. ,   4.5,   6. ,   7.5,   9. ,  10.5,  12.);
    ASSERT_VEC8_EQUAL(f1-fvec8(1, 2, 3, 4, 5, 6, 7, 8), -0.5, -1. , -1.5, -2. , -2.5, -3. , -3.5, -4.);
    ASSERT_VEC8_EQUAL(f1*fvec8(1, 2, 3, 4, 5, 6, 7, 8), 0.5,   2. ,   4.5,   8. ,  12.5,  18. ,  24.5,  32.);
    ASSERT_VEC8_EQUAL(f1/fvec8(1, 2, 3, 4, 5, 6, 7, 8), 0.5,  0.5,  0.5,  0.5,  0.5,  0.5,  0.5,  0.5);

    f1 = fvec8(0.5, 1.0, 1.5, 2.0,   2.5, 3.0, 3.5, 4.0);
    f1 += fvec8(1, 2, 3, 4, 5, 6, 7, 8);
    ASSERT_VEC8_EQUAL(f1, 1.5,   3. ,   4.5,   6. ,   7.5,   9. ,  10.5,  12.);
    f1 = fvec8(0.5, 1.0, 1.5, 2.0,   2.5, 3.0, 3.5, 4.0);
    f1 -= fvec8(1, 2, 3, 4, 5, 6, 7, 8);
    ASSERT_VEC8_EQUAL(f1, -0.5, -1. , -1.5, -2. , -2.5, -3. , -3.5, -4.);
    f1 = fvec8(0.5, 1.0, 1.5, 2.0,   2.5, 3.0, 3.5, 4.0);
    f1 *= fvec8(1, 2, 3, 4, 5, 6, 7, 8);
    ASSERT_VEC8_EQUAL(f1, 0.5,   2. ,   4.5,   8. ,  12.5,  18. ,  24.5,  32.);
    f1 = fvec8(0.5, 1.0, 1.5, 2.0,   2.5, 3.0, 3.5, 4.0);
    f1 /= fvec8(1, 2, 3, 4, 5, 6, 7, 8);
    ASSERT_VEC8_EQUAL(f1, 0.5,  0.5,  0.5,  0.5,  0.5,  0.5,  0.5,  0.5);
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

void testComparisons() {
    fvec8 v1(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
    fvec8 v2(1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5);
    ASSERT_VEC8_EQUAL(blend(v1, v2,
        fvec8(1.0, 1.5, 3.0, 2.2, 10.0, 10.5, 13.0, 12.2)==fvec8(1.1, 1.5, 3.0, 2.1, 10.1, 10.5, 13.0, 12.1)),
        0.0, 1.5, 1.5, 0.0, 0.0, 1.5, 1.5, 0.0);
    ASSERT_VEC8_EQUAL(blend(v1, v2,
        fvec8(1.0, 1.5, 3.0, 2.2, 10.0, 10.5, 13.0, 12.2)!=fvec8(1.1, 1.5, 3.0, 2.1, 10.1, 10.5, 13.0, 12.1)),
        1.5, 0.0, 0.0, 1.5, 1.5, 0.0, 0.0, 1.5);
    ASSERT_VEC8_EQUAL(blend(v1, v2,
        fvec8(1.0, 1.5, 3.0, 2.2, 10.0, 10.5, 13.0, 12.2)<fvec8(1.1, 1.5, 3.0, 2.1, 10.1, 10.5, 13.0, 12.1)),
        1.5, 0.0, 0.0, 0.0, 1.5, 0.0, 0.0, 0.0);
    ASSERT_VEC8_EQUAL(blend(v1, v2,
        fvec8(1.0, 1.5, 3.0, 2.2, 10.0, 10.5, 13.0, 12.2)>fvec8(1.1, 1.5, 3.0, 2.1, 10.1, 10.5, 13.0, 12.1)),
        0.0, 0.0, 0.0, 1.5, 0.0, 0.0, 0.0, 1.5);
    ASSERT_VEC8_EQUAL(blend(v1, v2,
        fvec8(1.0, 1.5, 3.0, 2.2, 10.0, 10.5, 13.0, 12.2)<=fvec8(1.1, 1.5, 3.0, 2.1, 10.1, 10.5, 13.0, 12.1)),
        1.5, 1.5, 1.5, 0.0, 1.5, 1.5, 1.5, 0.0);
    ASSERT_VEC8_EQUAL(blend(v1, v2,
        fvec8(1.0, 1.5, 3.0, 2.2, 10.0, 10.5, 13.0, 12.2)>=fvec8(1.1, 1.5, 3.0, 2.1, 10.1, 10.5, 13.0, 12.1)),
        0.0, 1.5, 1.5, 1.5, 0.0, 1.5, 1.5, 1.5);
}

void testMathFunctions() {
    fvec8 f1(0.4, 1.9, -1.2, -3.8, 0.4, 1.9, -1.2, -3.8);
    fvec8 f2(1.1, 1.2, 1.3, -5.0, 1.1, 1.2, 1.3, -5.0);
    ASSERT_VEC8_EQUAL(floor(f1), 0.0, 1.0, -2.0, -4.0, 0.0, 1.0, -2.0, -4.0);
    ASSERT_VEC8_EQUAL(ceil(f1), 1.0, 2.0, -1.0, -3.0, 1.0, 2.0, -1.0, -3.0);
    ASSERT_VEC8_EQUAL(round(f1), 0.0, 2.0, -1.0, -4.0, 0.0, 2.0, -1.0, -4.0);
    ASSERT_VEC8_EQUAL(abs(f1), 0.4, 1.9, 1.2, 3.8, 0.4, 1.9, 1.2, 3.8);
    ASSERT_VEC8_EQUAL(min(f1, f2), 0.4, 1.2, -1.2, -5.0, 0.4, 1.2, -1.2, -5.0);
    ASSERT_VEC8_EQUAL(max(f1, f2), 1.1, 1.9, 1.3, -3.8, 1.1, 1.9, 1.3, -3.8);
    ASSERT_VEC8_EQUAL(sqrt(fvec8(1.5, 3.1, 4.0, 15.0, 1.5, 3.1, 4.0, 15.0)), sqrt(1.5), sqrt(3.1), sqrt(4.0), sqrt(15.0), sqrt(1.5), sqrt(3.1), sqrt(4.0), sqrt(15.0));
    ASSERT_VEC8_EQUAL(rsqrt(fvec8(1.5, 3.1, 4.0, 15.0, 1.5, 3.1, 4.0, 15.0)), 1.0/sqrt(1.5), 1.0/sqrt(3.1), 1.0/sqrt(4.0), 1.0/sqrt(15.0), 1.0/sqrt(1.5), 1.0/sqrt(3.1), 1.0/sqrt(4.0), 1.0/sqrt(15.0));
    ASSERT_EQUAL_TOL(f1.lowerVec()[0]*f2.lowerVec()[0]+f1.lowerVec()[1]*f2.lowerVec()[1]+f1.lowerVec()[2]*f2.lowerVec()[2]+f1.lowerVec()[3]*f2.lowerVec()[3]+f1.upperVec()[0]*f2.upperVec()[0]+f1.upperVec()[1]*f2.upperVec()[1]+f1.upperVec()[2]*f2.upperVec()[2]+f1.upperVec()[3]*f2.upperVec()[3], dot8(f1, f2), 1e-6);
    ASSERT(any(f1 > 0.5));
    ASSERT(!any(f1 > 2.0));
    ASSERT_VEC8_EQUAL(blend(f1, f2, ivec8(-1, 0, -1, 0, -1, 0, -1, 0)), 1.1, 1.9, 1.3, -3.8, 1.1, 1.9, 1.3, -3.8);
}

void testTranspose() {
    fvec4 f1(0.0,   1.0,  2.0,  3.0);
    fvec4 f2(10.0, 11.0, 12.0, 13.0);
    fvec4 f3(20.0, 21.0, 22.0, 23.0);
    fvec4 f4(30.0, 31.0, 32.0, 33.0);
    fvec4 f5(40.0, 41.0, 42.0, 43.0);
    fvec4 f6(50.0, 51.0, 52.0, 53.0);
    fvec4 f7(60.0, 61.0, 62.0, 63.0);
    fvec4 f8(70.0, 71.0, 72.0, 73.0);
    fvec8 o1, o2, o3, o4;

    transpose(f1, f2, f3, f4, f5, f6, f7, f8, o1, o2, o3, o4);
    ASSERT_VEC8_EQUAL(o1, 0.0, 10.0, 20.0, 30.0, 40.0, 50.0, 60.0, 70.0);
    ASSERT_VEC8_EQUAL(o2, 1.0, 11.0, 21.0, 31.0, 41.0, 51.0, 61.0, 71.0);
    ASSERT_VEC8_EQUAL(o3, 2.0, 12.0, 22.0, 32.0, 42.0, 52.0, 62.0, 72.0);
    ASSERT_VEC8_EQUAL(o4, 3.0, 13.0, 23.0, 33.0, 43.0, 53.0, 63.0, 73.0);

    fvec4 g1, g2, g3, g4, g5, g6, g7, g8;
    transpose(o1, o2, o3, o4, g1, g2, g3, g4, g5, g6, g7, g8);
    ASSERT_VEC4_EQUAL(g1, 0.0,   1.0,  2.0,  3.0);
    ASSERT_VEC4_EQUAL(g2, 10.0, 11.0, 12.0, 13.0);
    ASSERT_VEC4_EQUAL(g3, 20.0, 21.0, 22.0, 23.0);
    ASSERT_VEC4_EQUAL(g4, 30.0, 31.0, 32.0, 33.0);
    ASSERT_VEC4_EQUAL(g5, 40.0, 41.0, 42.0, 43.0);
    ASSERT_VEC4_EQUAL(g6, 50.0, 51.0, 52.0, 53.0);
    ASSERT_VEC4_EQUAL(g7, 60.0, 61.0, 62.0, 63.0);
    ASSERT_VEC4_EQUAL(g8, 70.0, 71.0, 72.0, 73.0);

}

int main(int argc, char* argv[]) {
    try {
        if (!isVec8Supported()) {
            cout << "CPU is not supported.  Exiting." << endl;
            return 0;
        }
        testLoadStore();
        testArithmetic();
        testLogic();
        testComparisons();
        testMathFunctions();
        testTranspose();
    }
    catch(const exception& e) {
        cout << "exception: " << e.what() << endl;
        return 1;
    }
    cout << "Done" << endl;
    return 0;
}
