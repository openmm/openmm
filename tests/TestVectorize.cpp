/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2014-2015 Stanford University and the Authors.      *
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
 * This tests vectorized operations.
 */

#include "openmm/internal/AssertionUtilities.h"
#include "openmm/internal/vectorize.h"
#include <iostream>

using namespace OpenMM;
using namespace std;

#define ASSERT_VEC4_EQUAL(found, expected0, expected1, expected2, expected3) {if (std::abs((found)[0]-(expected0))>1e-6 || std::abs((found)[1]-(expected1))>1e-6 || std::abs((found)[2]-(expected2))>1e-6 || std::abs((found)[3]-(expected3))>1e-6) {std::stringstream details; details << " Expected ("<<(expected0)<<","<<(expected1)<<","<<(expected2)<<","<<(expected3)<<"), found ("<<(found)[0]<<","<<(found)[1]<<","<<(found)[2]<<","<<(found)[3]<<")"; throwException(__FILE__, __LINE__, details.str());}};
#define ASSERT_VEC4_EQUAL_INT(found, expected0, expected1, expected2, expected3) {if ((found)[0] != (expected0) || (found)[1] != (expected1) || (found)[2] != (expected2) || (found)[3] != (expected3)) {std::stringstream details; details << " Expected ("<<(expected0)<<","<<(expected1)<<","<<(expected2)<<","<<(expected3)<<"), found ("<<(found)[0]<<","<<(found)[1]<<","<<(found)[2]<<","<<(found)[3]<<")"; throwException(__FILE__, __LINE__, details.str());}};

void testLoadStore() {
    fvec4 f1(2.0);
    ivec4 i1(3);
    ASSERT_VEC4_EQUAL(f1, 2.0, 2.0, 2.0, 2.0);
    ASSERT_VEC4_EQUAL_INT(i1, 3, 3, 3, 3);
    fvec4 f2(2.5, 3.0, 3.5, 4.0);
    ivec4 i2(2, 3, 4, 5);
    ASSERT_VEC4_EQUAL(f2, 2.5, 3.0, 3.5, 4.0);
    ASSERT_VEC4_EQUAL_INT(i2, 2, 3, 4, 5);
    float farray[4];
    int iarray[4];
    f2.store(farray);
    i2.store(iarray);
    fvec4 f3(farray);
    ivec4 i3(iarray);
    ASSERT_VEC4_EQUAL(f3, 2.5, 3.0, 3.5, 4.0);
    ASSERT_VEC4_EQUAL_INT(i3, 2, 3, 4, 5);
    ASSERT_EQUAL(f3[0], 2.5);
    ASSERT_EQUAL(f3[1], 3.0);
    ASSERT_EQUAL(f3[2], 3.5);
    ASSERT_EQUAL(f3[3], 4.0);
    ASSERT_EQUAL(i3[0], 2);
    ASSERT_EQUAL(i3[1], 3);
    ASSERT_EQUAL(i3[2], 4);
    ASSERT_EQUAL(i3[3], 5);
}

void testArithmetic() {
    fvec4 f1(0.5, 1.0, 1.5, 2.0);
    ASSERT_VEC4_EQUAL(f1+fvec4(1, 2, 3, 4), 1.5, 3, 4.5, 6);
    ASSERT_VEC4_EQUAL(f1-fvec4(1, 2, 3, 4), -0.5, -1.0, -1.5, -2.0);
    ASSERT_VEC4_EQUAL(f1*fvec4(1, 2, 3, 4), 0.5, 2.0, 4.5, 8.0);
    ASSERT_VEC4_EQUAL(f1/fvec4(1, 2, 3, 4), 0.5, 0.5, 0.5, 0.5);
    ivec4 i1(1, 2, 3, 4);
    ASSERT_VEC4_EQUAL_INT(i1+ivec4(5, 2, 1, 3), 6, 4, 4, 7);
    ASSERT_VEC4_EQUAL_INT(i1-ivec4(5, 2, 1, 3), -4, 0, 2, 1);
    ASSERT_VEC4_EQUAL_INT(i1*ivec4(5, 2, 1, 3), 5, 4, 3, 12);
    f1 = fvec4(0.5, 1.0, 1.5, 2.0);
    f1 += fvec4(1, 2, 3, 4);
    ASSERT_VEC4_EQUAL(f1, 1.5, 3, 4.5, 6);
    f1 = fvec4(0.5, 1.0, 1.5, 2.0);
    f1 -= fvec4(1, 2, 3, 4);
    ASSERT_VEC4_EQUAL(f1, -0.5, -1.0, -1.5, -2.0);
    f1 = fvec4(0.5, 1.0, 1.5, 2.0);
    f1 *= fvec4(1, 2, 3, 4);
    ASSERT_VEC4_EQUAL(f1, 0.5, 2.0, 4.5, 8.0);
    f1 = fvec4(0.5, 1.0, 1.5, 2.0);
    f1 /= fvec4(1, 2, 3, 4);
    ASSERT_VEC4_EQUAL(f1, 0.5, 0.5, 0.5, 0.5);
    i1 = ivec4(1, 2, 3, 4);
    i1 += ivec4(5, 2, 1, 3);
    ASSERT_VEC4_EQUAL_INT(i1, 6, 4, 4, 7);
    i1 = ivec4(1, 2, 3, 4);
    i1 -= ivec4(5, 2, 1, 3);
    ASSERT_VEC4_EQUAL_INT(i1, -4, 0, 2, 1);
    i1 = ivec4(1, 2, 3, 4);
    i1 *= ivec4(5, 2, 1, 3);
    ASSERT_VEC4_EQUAL_INT(i1, 5, 4, 3, 12);
}

void testLogic() {
    int allBits = -1;
    float allBitsf = *((float*) &allBits);
    ivec4 mask(0, allBits, allBits, 0);
    fvec4 fmask(0, allBitsf, allBitsf, 0);;
    fvec4 f1(0.5, 1.0, 1.5, 2.0);
    ivec4 i1(1, 2, 3, 4);
    ASSERT_VEC4_EQUAL(f1&fmask, 0, 1.0, 1.5, 0);
    fvec4 temp = f1|fmask;
    ASSERT_EQUAL(0.5, temp[0]);
    ASSERT(temp[1] != temp[1]); // All bits set, which is nan
    ASSERT(temp[2] != temp[2]); // All bits set, which is nan
    ASSERT_EQUAL(2.0, temp[3]);
    ASSERT_VEC4_EQUAL_INT(i1&mask, 0, 2, 3, 0);
    ASSERT_VEC4_EQUAL_INT(i1|mask, 1, allBits, allBits, 4);
}

void testComparisons() {
    fvec4 v1(0.0, 0.0, 0.0, 0.0);
    fvec4 v2(1.5, 1.5, 1.5, 1.5);
    ASSERT_VEC4_EQUAL(blend(v1, v2, fvec4(1.0, 1.5, 3.0, 2.2)==fvec4(1.1, 1.5, 3.0, 2.1)), 0.0, 1.5, 1.5, 0.0);
    ASSERT_VEC4_EQUAL(blend(v1, v2, fvec4(1.0, 1.5, 3.0, 2.2)!=fvec4(1.1, 1.5, 3.0, 2.1)), 1.5, 0.0, 0.0, 1.5);
    ASSERT_VEC4_EQUAL(blend(v1, v2, fvec4(1.0, 1.5, 3.0, 2.2)<fvec4(1.1, 1.5, 3.0, 2.1)), 1.5, 0.0, 0.0, 0.0);
    ASSERT_VEC4_EQUAL(blend(v1, v2, fvec4(1.0, 1.5, 3.0, 2.2)>fvec4(1.1, 1.5, 3.0, 2.1)), 0.0, 0.0, 0.0, 1.5);
    ASSERT_VEC4_EQUAL(blend(v1, v2, fvec4(1.0, 1.5, 3.0, 2.2)<=fvec4(1.1, 1.5, 3.0, 2.1)), 1.5, 1.5, 1.5, 0.0);
    ASSERT_VEC4_EQUAL(blend(v1, v2, fvec4(1.0, 1.5, 3.0, 2.2)>=fvec4(1.1, 1.5, 3.0, 2.1)), 0.0, 1.5, 1.5, 1.5);
    fvec4 imask(3, 3, 3, 3);
    ASSERT_VEC4_EQUAL_INT((ivec4(1, 3, 7, 5)==ivec4(2, 3, 7, 4))&imask, 0, 3, 3, 0);
    ASSERT_VEC4_EQUAL_INT((ivec4(1, 3, 7, 5)!=ivec4(2, 3, 7, 4))&imask, 3, 0, 0, 3);
    ASSERT_VEC4_EQUAL_INT((ivec4(1, 3, 7, 5)<ivec4(2, 3, 7, 4))&imask, 3, 0, 0, 0);
    ASSERT_VEC4_EQUAL_INT((ivec4(1, 3, 7, 5)>ivec4(2, 3, 7, 4))&imask, 0, 0, 0, 3);
    ASSERT_VEC4_EQUAL_INT((ivec4(1, 3, 7, 5)<=ivec4(2, 3, 7, 4))&imask, 3, 3, 3, 0);
    ASSERT_VEC4_EQUAL_INT((ivec4(1, 3, 7, 5)>=ivec4(2, 3, 7, 4))&imask, 0, 3, 3, 3);
}

void testMathFunctions() {
    fvec4 f1(0.4, 1.9, -1.2, -3.8);
    fvec4 f2(1.1, 1.2, 1.3, -5.0);
    ASSERT_VEC4_EQUAL(floor(f1), 0.0, 1.0, -2.0, -4.0);
    ASSERT_VEC4_EQUAL(ceil(f1), 1.0, 2.0, -1.0, -3.0);
    ASSERT_VEC4_EQUAL(round(f1), 0.0, 2.0, -1.0, -4.0);
    ASSERT_VEC4_EQUAL(abs(f1), 0.4, 1.9, 1.2, 3.8);
    ASSERT_VEC4_EQUAL(min(f1, f2), 0.4, 1.2, -1.2, -5.0);
    ASSERT_VEC4_EQUAL(max(f1, f2), 1.1, 1.9, 1.3, -3.8);
    ASSERT_VEC4_EQUAL(sqrt(fvec4(1.5, 3.1, 4.0, 15.0)), sqrt(1.5), sqrt(3.1), sqrt(4.0), sqrt(15.0));
    ASSERT_VEC4_EQUAL(rsqrt(fvec4(1.5, 3.1, 4.0, 15.0)), 1.0/sqrt(1.5), 1.0/sqrt(3.1), 1.0/sqrt(4.0), 1.0/sqrt(15.0));
    ASSERT_VEC4_EQUAL(exp(fvec4(-2.1, -0.5, 1.5, 3.1)), expf(-2.1), expf(-0.5), expf(1.5), expf(3.1));
    ASSERT_VEC4_EQUAL(log(fvec4(1.5, 3.1, 4.0, 15.0)), logf(1.5), logf(3.1), logf(4.0), logf(15.0));
    ASSERT_EQUAL_TOL(f1[0]*f2[0]+f1[1]*f2[1]+f1[2]*f2[2], dot3(f1, f2), 1e-6);
    ASSERT_EQUAL_TOL(f1[0]*f2[0]+f1[1]*f2[1]+f1[2]*f2[2]+f1[3]*f2[3], dot4(f1, f2), 1e-6);
    ASSERT(any(f1 > 0.5));
    ASSERT(!any(f1 > 2.0));
    ASSERT_VEC4_EQUAL(blend(f1, f2, ivec4(-1, 0, -1, 0)), 1.1, 1.9, 1.3, -3.8);
    ASSERT_VEC4_EQUAL(cross(f1, f2), 3.91, -1.84, -1.61, 0.0);
}

void testTranspose() {
    fvec4 f1(1.0, 2.0, 3.0, 4.0);
    fvec4 f2(5.0, 6.0, 7.0, 8.0);
    fvec4 f3(9.0, 10.0, 11.0, 12.0);
    fvec4 f4(13.0, 14.0, 15.0, 16.0);
    transpose(f1, f2, f3, f4);
    ASSERT_VEC4_EQUAL(f1, 1.0, 5.0, 9.0, 13.0);
    ASSERT_VEC4_EQUAL(f2, 2.0, 6.0, 10.0, 14.0);
    ASSERT_VEC4_EQUAL(f3, 3.0, 7.0, 11.0, 15.0);
    ASSERT_VEC4_EQUAL(f4, 4.0, 8.0, 12.0, 16.0);
}

int main(int argc, char* argv[]) {
    try {
        if (!isVec4Supported()) {
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
