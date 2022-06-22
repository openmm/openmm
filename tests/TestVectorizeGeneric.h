/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2014-2020 Stanford University and the Authors.      *
 * Authors: Daniel Towner                                                     *
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

#pragma once

/**
 * This tests all sizes of vectorized operations using templated test code.
 */

#include <array>
#include <cstring>
#include <functional>
#include <iostream>
#include <iterator>
#include <numeric>
#include <memory.h>
#include <sstream>
#include <typeinfo>

/**
 * Return the 32-bit integer bit pattern from the given floating-point value.
 */
static int32_t floatAsIntBits(float f) {
    int32_t i;
    memcpy(&i, &f, 4);
    return i;
}

/**
 * Compare two floating-point values using units-in-last-place (ULP) as a measure of equality. Two values
 * which are only a few representable values apart can be considered to be equal. Note that IEEE
 * operations (add, mul, etc.) will always be exact, but sequences of operations might be more than
 * a few ULP apart, but still close enough to be considered equal. ULP comparisons work at any scale of
 * number, unlike an epsilon-based approach.
 */
static bool almostEqual(float a, float b) {
    // Maybe they really are equal.
    if (a == b)
        return true;

    // Infinities and NANs are never equal to anything, even other nans and infinities.
    if (std::isnan(a) || std::isinf(a) ||
        std::isnan(b) || std::isinf(b))
        return false;

    // If they are different signs then they can't be equal. For two very small denormal values they might
    // be very close to each other but either side of 0, but denormals are a corner case which don't deserve
    // to be equal.
    if (std::signbit(a) != std::signbit(b))
        return false;

    // The two numbers must be valid values with the same sign, so treat then as basic integers to
    // get at their ULP values. If they are only a few ULP apart, then they are essentially equal.
    int32_t intDiff = std::abs(floatAsIntBits(a) - floatAsIntBits(b));
    return intDiff < 4;
}

static bool exactlyEqual(float a, float b) { return a == b; }

/**
 * Write the contents of the given array-like object to a stream. No formatting is applied.
 */
template<typename FVEC>
void VecToStream(std::ostream& stream, const FVEC& vec)
{
    constexpr int numElements = sizeof(FVEC) / sizeof(float);
    const float* vptr = (const float*)&vec;
    for (int i=0; i<numElements; ++i)
        stream << vptr[i] << ", ";
}

/**
 * Given two vector-like objects compared each of their elements for equality. The vector objects can be
 * anything which in memory is a list of 32-bit floating-point values, so SIMD vectors, C arrays or
 * C++ arrays would all be valid.
 */
template<typename S, typename T>
static void checkElementsEqual(const S& computed, const T& expected,
                               std::function<bool(float, float)> equal_fn,
                               const char* file, int line) {
    // Both S and T should be arrays of floats of the same length.
    static_assert(sizeof(T) == sizeof(S), "Array-like elements must have the same size");

    constexpr int numElements = sizeof(S) / sizeof(float);

    const float* computedPtr = (const float*)&computed;
    const float* expectedPtr = (const float*)&expected;

    std::ostringstream details;
    details << "Error during test for type " << typeid(S).name() << '\n';

    bool passed = true;
    for (int i=0; i<numElements; ++i)
    {
        if (!equal_fn(computedPtr[i], expectedPtr[i]))
            passed = false;
    }

    if (!passed)
    {
        details << "Values differ. ";
        VecToStream(details, computed);
        details << " and ";
        VecToStream(details, expected);
        OpenMM::throwException(file, line, details.str());
    }

}

#define ASSERT_VEC_EQUAL(computed, expected) {checkElementsEqual(computed, expected, exactlyEqual, __FILE__, __LINE__);}
#define ASSERT_VEC_ALMOST_EQUAL(computed, expected) {checkElementsEqual(computed, expected, almostEqual, __FILE__, __LINE__);}

static float getRandomFloat () {
    // Between -50 and 50.
    return float(rand()) / float(RAND_MAX/100.0f) - 50.0f;
}

/**
 * Given an array-like memory object containing floats, apply the given function to every element.
 */
template<typename FVEC>
FVEC applyUnaryFn(const FVEC& v, std::function<float(float)> fn) {
    constexpr int numElements = sizeof(FVEC) / sizeof(float);

    FVEC result;

    float* rp = (float*)&result;
    const float* vp = (const float*)&v;

    for (int i=0; i<numElements; ++i)
        rp[i] = fn(vp[i]);

    return result;
}

/**
 * Given an array-like memory object containing floats, apply the given function to every element.
 */
template<typename FVEC>
FVEC applyBinaryFn(const FVEC& a, const FVEC& b, std::function<float(float, float)> fn) {
    constexpr int numElements = sizeof(FVEC) / sizeof(float);

    FVEC result;

    float* rp = (float*)&result;
    const float* ap = (const float*)&a;
    const float* bp = (const float*)&b;

    for (int i=0; i<numElements; ++i)
        rp[i] = fn(ap[i], bp[i]);

    return result;
}

/**
 * Provide a test fixture class which underpins all verification for a given
 * type of vector SIMD implementation, as well as providing common utility functions
 */
template<typename FVEC>
class TestFvec {
public:

    static constexpr int numElements = sizeof(FVEC) / sizeof(float);

    void testInitializers() const;
    void testUnaryOps() const;
    void testBinaryOps() const;
    void testUtilities() const;
    void testBlendAndCompare() const;
    void testTranspose() const;

    static void testAll() {
        TestFvec<FVEC> testUnit;
        testUnit.testInitializers();
        testUnit.testUnaryOps();
        testUnit.testBinaryOps();
        testUnit.testUtilities();
        testUnit.testBlendAndCompare();
        testUnit.testTranspose();
    }

    FVEC getRandomFvec() const {
        union {
            FVEC v;
            float f[numElements];
        };

        for (auto& e : f)
            e = getRandomFloat();

        return v;
    }

};

template<typename FVEC>
void TestFvec<FVEC>::testInitializers() const {
    FVEC computedZero = {};
    float expectedZero[numElements] = {};
    ASSERT_VEC_EQUAL(computedZero, expectedZero);

    FVEC computedBroadcast(14.5f);
    float expectedBroadcast[numElements];
    std::fill_n(expectedBroadcast, numElements, 14.5f);
    ASSERT_VEC_EQUAL(computedBroadcast, expectedBroadcast);

    float expectedArray[numElements];
    std::iota(expectedArray, expectedArray + numElements, 23);
    FVEC computedFromLoad(expectedArray);
    ASSERT_VEC_EQUAL(computedFromLoad, expectedArray);

    // Gather values from a table. Variants for both one vector and two vector gathers are provided.
    // The indexes to gather (multiples of 7) are also generated, along with the expected answers.
    float gatherTable[2048];
    for (int i=0; i<2048;++i)
        gatherTable[i] = -i; // Same index to make it easy to debug, but negative to avoid copying idx.

    int gatherIndexes[numElements];
    float gatherIndexesAsFloat[numElements]; // Same as above, but in float format.
    float expectedGather0[numElements];
    float expectedGather1[numElements];
    for (int i=0; i<numElements; ++i)
    {
        gatherIndexes[i] = i * 7;
        gatherIndexesAsFloat[i] = float(gatherIndexes[i]);
        expectedGather0[i] = -(i * 7);
        expectedGather1[i] = -(i * 7) - 1; // Each value is one less than previous.
    }

    // Single value gather
    FVEC computedFromGather(gatherTable, gatherIndexes);
    ASSERT_VEC_EQUAL(computedFromGather, expectedGather0);

    // Pair-wise vector gather. The first values should be the same as a normal gather, and the
    // second are just increments from the first. Note that there musty be some suitable conversion
    // from a floating-point index (i.e., an integer value in float format), and the type required
    // for the second operand of gatherVecPair. gatherVecPair can then take either an actual
    // float vector, or some suitable format like ivec4 or ivec8.
    FVEC findex(gatherIndexesAsFloat);
    FVEC p0, p1;
    gatherVecPair(gatherTable, findex, p0, p1);
    ASSERT_VEC_EQUAL(p0, expectedGather0);
    ASSERT_VEC_EQUAL(p1, expectedGather1);
}

template<typename FVEC>
void TestFvec<FVEC>::testUnaryOps() const {
    const auto v = getRandomFvec();

    // Note that these are exact comparisons because all these SIMD operators are
    // just applying the scalar operator, so there should be no loss of precision.

    ASSERT_VEC_EQUAL(abs(v), applyUnaryFn(v, [](float x) { return std::abs(x);} ));

    ASSERT_VEC_EQUAL(-v, applyUnaryFn(v, [](float x) { return 0 - x;} ));

    ASSERT_VEC_EQUAL(floor(v), applyUnaryFn(v, [](float x) { return std::floor(x);} ));
    ASSERT_VEC_EQUAL(ceil(v), applyUnaryFn(v, [](float x) { return std::ceil(x);} ));
    ASSERT_VEC_EQUAL(round(v), applyUnaryFn(v, [](float x) { return std::round(x);} ));

    // Borrow a few other functions to test sqrt neatly.
    const auto positiveValue = abs(v) + 1;
    ASSERT_VEC_ALMOST_EQUAL(sqrt(positiveValue * positiveValue), positiveValue);
    ASSERT_VEC_ALMOST_EQUAL(rsqrt(positiveValue * positiveValue), 1.0f / abs(positiveValue));
}

template<typename FVEC>
void TestFvec<FVEC>::testBinaryOps() const {
    const auto v0 = getRandomFvec();
    const auto v1 = getRandomFvec();

    // Note that most of these are exact comparisons because all these SIMD operators are
    // just applying the scalar operator, so there should be no loss of precision. The one
    // exception is division, which does often do something slightly different
    // since division is an expensive operation (e.g., multiply by reciprocal).

    // Binary operators.
    ASSERT_VEC_EQUAL(v0 + v1, applyBinaryFn(v0, v1, std::plus<float>()));
    ASSERT_VEC_EQUAL(v0 - v1, applyBinaryFn(v0, v1, std::minus<float>()));
    ASSERT_VEC_EQUAL(v0 * v1, applyBinaryFn(v0, v1, std::multiplies<float>()));
    ASSERT_VEC_ALMOST_EQUAL(v0 / v1, applyBinaryFn(v0, v1, std::divides<float>()));

    // Assignment operators.
    auto addAssign = v0;
    addAssign += v1;
    ASSERT_VEC_EQUAL(addAssign, applyBinaryFn(v0, v1, std::plus<float>()));

    auto subAssign = v0;
    subAssign -= v1;
    ASSERT_VEC_EQUAL(subAssign, applyBinaryFn(v0, v1, std::minus<float>()));

    auto mulAssign = v0;
    mulAssign *= v1;
    ASSERT_VEC_EQUAL(mulAssign, applyBinaryFn(v0, v1, std::multiplies<float>()));

    auto divAssign = v0;
    divAssign /= v1;
    ASSERT_VEC_ALMOST_EQUAL(divAssign, applyBinaryFn(v0, v1, std::divides<float>()));

    // Binary ops between SIMD and scalar.
    const float f = getRandomFloat();
    const FVEC fdup(f);

    ASSERT_VEC_EQUAL(v0 + f, applyBinaryFn(v0, fdup, std::plus<float>()));
    ASSERT_VEC_EQUAL(f + v0, applyBinaryFn(fdup, v0, std::plus<float>()));
    ASSERT_VEC_EQUAL(v0 - f, applyBinaryFn(v0, fdup, std::minus<float>()));
    ASSERT_VEC_EQUAL(f - v0, applyBinaryFn(fdup, v0, std::minus<float>()));
    ASSERT_VEC_EQUAL(v0 * f, applyBinaryFn(v0, fdup, std::multiplies<float>()));
    ASSERT_VEC_EQUAL(f * v0, applyBinaryFn(fdup, v0, std::multiplies<float>()));
    ASSERT_VEC_ALMOST_EQUAL(v0 / f, applyBinaryFn(v0, fdup, std::divides<float>()));
    ASSERT_VEC_ALMOST_EQUAL(f / v0, applyBinaryFn(fdup, v0, std::divides<float>()));

    // Binary functions.
    ASSERT_VEC_EQUAL(min(v0, v1),
                     applyBinaryFn(v0, v1, [](float x, float y) { return std::min<float>(x, y); }));
    ASSERT_VEC_EQUAL(max(v0, v1),
                     applyBinaryFn(v0, v1, [](float x, float y) { return std::max<float>(x, y); }));
}

template<typename FVEC>
void TestFvec<FVEC>::testTranspose() const {

    // A table of random data to transpose.
    float table[numElements * 4];
    for (auto& e : table) e = std::round(getRandomFloat());

    // Load the table row data into vectors.
    const auto i0 = FVEC(table + 0 * numElements);
    const auto i1 = FVEC(table + 1 * numElements);
    const auto i2 = FVEC(table + 2 * numElements);
    const auto i3 = FVEC(table + 3 * numElements);

    // Manually transpose the data.
    std::array<float, numElements * 4> expectedTranspose;
    for (auto r=0; r<4; ++r)
    {
        for (auto c=0; c<numElements; ++c)
        {
            expectedTranspose[c * 4 + r] = table[r * numElements + c];
        }
    }

    fvec4 computedTranspose[numElements];
    transpose(i0, i1, i2, i3, computedTranspose);

    ASSERT_VEC_EQUAL(computedTranspose, expectedTranspose);

    FVEC o0, o1, o2, o3;
    transpose(computedTranspose, o0, o1, o2, o3);

    ASSERT_VEC_EQUAL(i0, o0);
    ASSERT_VEC_EQUAL(i1, o1);
    ASSERT_VEC_EQUAL(i2, o2);
    ASSERT_VEC_EQUAL(i3, o3);
}

template<typename FVEC>
void TestFvec<FVEC>::testBlendAndCompare() const {
    const FVEC zero = {};
    const FVEC allOne(1.0f);
    const FVEC allTwo(2.0f);

    // Note that different targets use different types of mask, so rather than checking
    // the mask directly, instead check the output of using the mask as a blend to provide
    // an indirect test.

    const auto maskNone = FVEC::expandBitsToMask(0);
    ASSERT_VEC_EQUAL(blend(allOne, allTwo, maskNone), allOne);
    ASSERT_VEC_EQUAL(blendZero(allOne, maskNone), zero);

    const auto maskAll = FVEC::expandBitsToMask(-1);
    ASSERT_VEC_EQUAL(blend(allOne, allTwo, maskAll), allTwo);
    ASSERT_VEC_EQUAL(blendZero(allOne, maskAll), allOne);

    // Repeating pattern big enough to do most SIMD lengths.
    const int bitmask = 0b1100001101101001;
    const auto maskSome = FVEC::expandBitsToMask(bitmask);
    float expectedMaskSome[numElements];
    float expectedZeroMaskSome[numElements];
    for (int i=0; i<numElements; ++i)
    {
        expectedMaskSome[i] = (bitmask & (1 << i)) ? 2.0f : 1.0f;
        expectedZeroMaskSome[i] = (bitmask & (1 << i)) ? 2.0f : 0.0f;
    }
    ASSERT_VEC_EQUAL(blend(allOne, allTwo, maskSome), expectedMaskSome);
    ASSERT_VEC_EQUAL(blendZero(allTwo, maskSome), expectedZeroMaskSome);

    // Test comparisons too, using random numbers, and then blending in either 0 or 1.
    const auto v0 = getRandomFvec();
    const auto v1 = getRandomFvec();
    ASSERT_VEC_EQUAL(blend(allOne, allTwo, v0 < v1),
                     applyBinaryFn(v0, v1, [](float x, float y) { return x < y ? 2.0f : 1.0f; }));
    ASSERT_VEC_EQUAL(blend(allOne, allTwo, v0 <= v1),
                     applyBinaryFn(v0, v1, [](float x, float y) { return x <= y ? 2.0f : 1.0f; }));
    ASSERT_VEC_EQUAL(blend(allOne, allTwo, v0 <= v0), allTwo);
    ASSERT_VEC_EQUAL(blend(allOne, allTwo, v0 > v1),
                     applyBinaryFn(v0, v1, [](float x, float y) { return x > y ? 2.0f : 1.0f; }));
    ASSERT_VEC_EQUAL(blend(allOne, allTwo, v0 >= v1),
                     applyBinaryFn(v0, v1, [](float x, float y) { return x >= y ? 2.0f : 1.0f; }));
    ASSERT_VEC_EQUAL(blend(allOne, allTwo, v0 >= v0), allTwo);

}

template<typename FVEC>
void TestFvec<FVEC>::testUtilities() const {
    /** Use rounded (i.e., integer) values for the reductions. Reduction operations are very sensitive
     * to ordering. The correct result is found by sorting values into ascending order to ensure that
     * similar sized numbers are accumulated earlier than less similar numbers. If completely random
     * numbers were used, this effect would show up here, making it more a test of what random numbers
     * you got, than of the code itself. By rounding to integers, the numbers will behave sanely for the
     * reduction, meaning it is a test of the reduction, and not of the format.
     */
    const auto v0 = round(getRandomFvec());
    const auto v1 = round(getRandomFvec());
    const auto v2 = round(getRandomFvec());

    const float* v0p = (const float*)&v0;
    const float* v1p = (const float*)&v1;
    const float* v2p = (const float*)&v2;

    const auto expectedRedAddV0 = std::accumulate(v0p, v0p + numElements, 0.0f);
    const auto expectedRedAddV1 = std::accumulate(v1p, v1p + numElements, 0.0f);
    const auto expectedRedAddV2 = std::accumulate(v2p, v2p + numElements, 0.0f);

    ASSERT_VEC_EQUAL(reduceAdd(v0), expectedRedAddV0);

    // Reduction of three vectors by addition into a single 3-element vector. Note that the final element
    // of the reduction is undefined, so the expected value copies over whatever that undefined value is.
    const auto computedRed3 = reduceToVec3(v0, v1, v2);
    const auto expectedRed3 = fvec4(expectedRedAddV0, expectedRedAddV1, expectedRedAddV2, computedRed3[3]);
    ASSERT_VEC_EQUAL(computedRed3, expectedRed3);

}