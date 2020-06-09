#ifndef OPENMM_VECTORIZEAVX_H_
#define OPENMM_VECTORIZEAVX_H_

/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2013-2014 Stanford University and the Authors.      *
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

#include "vectorize.h"
#include <immintrin.h>

// This file defines classes and functions to simplify vectorizing code with AVX.

class ivec8;

/**
 * An eight element vector of floats.
 */
class fvec8 {
public:
    __m256 val;

    fvec8() = default;
    fvec8(float v) : val(_mm256_set1_ps(v)) {}
    fvec8(float v1, float v2, float v3, float v4, float v5, float v6, float v7, float v8) : val(_mm256_set_ps(v8, v7, v6, v5, v4, v3, v2, v1)) {}
    fvec8(__m256 v) : val(v) {}
    fvec8(const float* v) : val(_mm256_loadu_ps(v)) {}

    /** Create a vector by gathering individual indexes of data from a table. Element i of the vector will
     * be loaded from table[idx[i]].
     * @param table The table from which to do a lookup.
     * @param indexes The indexes to gather.
     */
    fvec8(const float* table, const int32_t idx[8]) {
        val = _mm256_setr_ps(table[idx[0]], table[idx[1]], table[idx[2]], table[idx[3]], table[idx[4]], table[idx[5]], table[idx[6]], table[idx[7]]);
    }

    operator __m256() const {
        return val;
    }
    fvec4 lowerVec() const {
        return _mm256_castps256_ps128(val);
    }
    fvec4 upperVec() const {
        return _mm256_extractf128_ps(val, 1);
    }
    void store(float* v) const {
        _mm256_storeu_ps(v, val);
    }
    fvec8 operator+(fvec8 other) const {
        return _mm256_add_ps(val, other);
    }
    fvec8 operator-(fvec8 other) const {
        return _mm256_sub_ps(val, other);
    }
    fvec8 operator*(fvec8 other) const {
        return _mm256_mul_ps(val, other);
    }
    fvec8 operator/(fvec8 other) const {
        return _mm256_div_ps(val, other);
    }
    void operator+=(fvec8 other) {
        val = _mm256_add_ps(val, other);
    }
    void operator-=(fvec8 other) {
        val = _mm256_sub_ps(val, other);
    }
    void operator*=(fvec8 other) {
        val = _mm256_mul_ps(val, other);
    }
    void operator/=(fvec8 other) {
        val = _mm256_div_ps(val, other);
    }
    fvec8 operator-() const {
        return _mm256_sub_ps(_mm256_set1_ps(0.0f), val);
    }
    fvec8 operator&(fvec8 other) const {
        return _mm256_and_ps(val, other);
    }
    fvec8 operator|(fvec8& other) const {
        return _mm256_or_ps(val, other);
    }
    fvec8 operator==(fvec8 other) const {
        return _mm256_cmp_ps(val, other, _CMP_EQ_OQ);
    }
    fvec8 operator!=(fvec8 other) const {
        return _mm256_cmp_ps(val, other, _CMP_NEQ_OQ);
    }
    fvec8 operator>(fvec8 other) const {
        return _mm256_cmp_ps(val, other, _CMP_GT_OQ);
    }
    fvec8 operator<(fvec8 other) const {
        return _mm256_cmp_ps(val, other, _CMP_LT_OQ);
    }
    fvec8 operator>=(fvec8 other) const {
        return _mm256_cmp_ps(val, other, _CMP_GE_OQ);
    }
    fvec8 operator<=(fvec8 other) const {
        return _mm256_cmp_ps(val, other, _CMP_LE_OQ);
    }
    operator ivec8() const;

    /**
     * Convert an integer bitmask into a full vector of elements which can be used
     * by the blend function.
     */
    static fvec8 expandBitsToMask(int bitmask);
};

/**
 * An eight element vector of ints.
 */
class ivec8 {
public:
    __m256i val;

    ivec8() {}
    ivec8(int v) : val(_mm256_set1_epi32(v)) {}
    ivec8(int v1, int v2, int v3, int v4, int v5, int v6, int v7, int v8) : val(_mm256_set_epi32(v8, v7, v6, v5, v4, v3, v2, v1)) {}
    ivec8(__m256i v) : val(v) {}
    ivec8(const int* v) : val(_mm256_loadu_si256((const __m256i*) v)) {}
    operator __m256i() const {
        return val;
    }
    ivec4 lowerVec() const {
        return _mm256_castsi256_si128(val);
    }
    ivec4 upperVec() const {
        return _mm256_extractf128_si256(val, 1);
    }
    void store(int* v) const {
        _mm256_storeu_si256((__m256i*) v, val);
    }
    ivec8 operator&(ivec8 other) const {
        return _mm256_castps_si256(_mm256_and_ps(_mm256_castsi256_ps(val), _mm256_castsi256_ps(other.val)));
    }
    ivec8 operator|(ivec8 other) const {
        return _mm256_castps_si256(_mm256_or_ps(_mm256_castsi256_ps(val), _mm256_castsi256_ps(other.val)));
    }
    operator fvec8() const;
};

// Conversion operators.

inline fvec8::operator ivec8() const {
    return _mm256_cvttps_epi32(val);
}

inline ivec8::operator fvec8() const {
    return _mm256_cvtepi32_ps(val);
}

inline fvec8 fvec8::expandBitsToMask(int bitmask) {
    // Put a copy of bit 0 in the first element, bit 1 in the second, and so on.
    const auto expandedBits =
      _mm256_and_ps(_mm256_castsi256_ps(_mm256_set1_epi8(bitmask)),
                    _mm256_castsi256_ps(_mm256_setr_epi32(1, 2, 4, 8, 16, 32, 64, 128)));

    // The individual bits are essentially extremely small floating-point values. By comparing against zero
    //  (even a floating-point zero), the individual bits are turned into a complete element mask.
    const auto elementMask = _mm256_cmp_ps(expandedBits, __m256(), _CMP_NEQ_OQ);

    return elementMask;
}

// Functions that operate on fvec8s.

static inline fvec8 floor(fvec8 v) {
    return fvec8(_mm256_round_ps(v.val, 0x09));
}

static inline fvec8 ceil(fvec8 v) {
    return fvec8(_mm256_round_ps(v.val, 0x0A));
}

static inline fvec8 round(fvec8 v) {
    return fvec8(_mm256_round_ps(v.val, _MM_FROUND_TO_NEAREST_INT));
}

static inline fvec8 min(fvec8 v1, fvec8 v2) {
    return fvec8(_mm256_min_ps(v1.val, v2.val));
}

static inline fvec8 max(fvec8 v1, fvec8 v2) {
    return fvec8(_mm256_max_ps(v1.val, v2.val));
}

static inline fvec8 abs(fvec8 v) {
    static const __m256 mask = _mm256_castsi256_ps(_mm256_set1_epi32(0x7FFFFFFF));
    return fvec8(_mm256_and_ps(v.val, mask));
}

static inline fvec8 sqrt(fvec8 v) {
    return fvec8(_mm256_sqrt_ps(v.val));
}

static inline fvec8 rsqrt(fvec8 v) {
    // Initial estimate of rsqrt().

    fvec8 y(_mm256_rsqrt_ps(v.val));

    // Perform an iteration of Newton refinement.

    fvec8 x2 = v*0.5f;
    y *= fvec8(1.5f)-x2*y*y;
    return y;
}

static inline float dot8(fvec8 v1, fvec8 v2) {
    fvec8 result = _mm256_dp_ps(v1, v2, 0xF1);
    return _mm_cvtss_f32(result.lowerVec())+_mm_cvtss_f32(result.upperVec());
}

static inline float reduceAdd(fvec8 v) {
    // :TODO: There are more efficient ways to do this.
    return dot8(v, fvec8(1.0f));
}

static inline void transpose(fvec4 in1, fvec4 in2, fvec4 in3, fvec4 in4, fvec4 in5, fvec4 in6, fvec4 in7, fvec4 in8, fvec8& out1, fvec8& out2, fvec8& out3, fvec8& out4) {
    fvec4 i1 = in1, i2 = in2, i3 = in3, i4 = in4;
    fvec4 i5 = in5, i6 = in6, i7 = in7, i8 = in8;
    _MM_TRANSPOSE4_PS(i1, i2, i3, i4);
    _MM_TRANSPOSE4_PS(i5, i6, i7, i8);
#ifdef _MSC_VER
    // Visual Studio has a bug in _mm256_castps128_ps256, so we have to use the more expensive _mm256_insertf128_ps.
    out1 = _mm256_insertf128_ps(out1, i1, 0);
    out2 = _mm256_insertf128_ps(out2, i2, 0);
    out3 = _mm256_insertf128_ps(out3, i3, 0);
    out4 = _mm256_insertf128_ps(out4, i4, 0);
#else
    out1 = _mm256_castps128_ps256(i1);
    out2 = _mm256_castps128_ps256(i2);
    out3 = _mm256_castps128_ps256(i3);
    out4 = _mm256_castps128_ps256(i4);
#endif
    out1 = _mm256_insertf128_ps(out1, i5, 1);
    out2 = _mm256_insertf128_ps(out2, i6, 1);
    out3 = _mm256_insertf128_ps(out3, i7, 1);
    out4 = _mm256_insertf128_ps(out4, i8, 1);
}

/** Given a vec4[8] input array, generate 4 vec8 outputs. The first output contains all the first elements
 * the second output the second elements, and so on. Note that the prototype is essentially differing only
 * in output type so it can be overloaded in other SIMD fvec types.
 */
static inline void transpose(const fvec4 in[8], fvec8& out1, fvec8& out2, fvec8& out3, fvec8& out4) {
    transpose(in[0], in[1], in[2], in[3], in[4], in[5], in[6], in[7], out1, out2, out3, out4);
}

static inline void transpose(fvec8 in1, fvec8 in2, fvec8 in3, fvec8 in4, fvec4& out1, fvec4& out2, fvec4& out3, fvec4& out4, fvec4& out5, fvec4& out6, fvec4& out7, fvec4& out8) {
    out1 = in1.lowerVec();
    out2 = in2.lowerVec();
    out3 = in3.lowerVec();
    out4 = in4.lowerVec();
    _MM_TRANSPOSE4_PS(out1, out2, out3, out4);
    out5 = in1.upperVec();
    out6 = in2.upperVec();
    out7 = in3.upperVec();
    out8 = in4.upperVec();
    _MM_TRANSPOSE4_PS(out5, out6, out7, out8);
}

/**
 * Given 4 input vectors of 8 elements, transpose them to form 8 output vectors of 4 elements.
 */
static inline void transpose(fvec8 in1, fvec8 in2, fvec8 in3, fvec8 in4, fvec4 out[8]) {
    transpose(in1, in2, in3, in4, out[0], out[1], out[2], out[3], out[4], out[5], out[6], out[7]);
}

// Functions that operate on ivec8s.

static inline bool any(ivec8 v) {
    return !_mm256_testz_si256(v, _mm256_set1_epi32(0xFFFFFFFF));
}

// Mathematical operators involving a scalar and a vector.

static inline fvec8 operator+(float v1, fvec8 v2) {
    return fvec8(v1)+v2;
}

static inline fvec8 operator-(float v1, fvec8 v2) {
    return fvec8(v1)-v2;
}

static inline fvec8 operator*(float v1, fvec8 v2) {
    return fvec8(v1)*v2;
}

static inline fvec8 operator/(float v1, fvec8 v2) {
    return fvec8(v1)/v2;
}

// Operation for blending fvec8 from a full bitmask.
static inline fvec8 blend(fvec8 v1, fvec8 v2, fvec8 mask) {
    return fvec8(_mm256_blendv_ps(v1.val, v2.val, mask.val));
}

static inline fvec8 blendZero(fvec8 v, fvec8 mask) {
    return blend(0.0f, v, mask);
}

/**
 * Given a table of floating-point values and a set of indexes, perform a gather read into a pair
 * of vectors. The first result vector contains the values at the given indexes, and the second
 * result vector contains the values from each respective index+1.
 */
static inline void gatherVecPair(const float* table, ivec8 index, fvec8& out0, fvec8& out1) {

    const auto lower = index.lowerVec();
    const auto upper = index.upperVec();

    // Gather all the separate memory data together. Each vector will have two values
    // which get used, and two which are ultimately discarded.
    fvec4 t0(table + lower[0]);
    fvec4 t1(table + lower[1]);
    fvec4 t2(table + lower[2]);
    fvec4 t3(table + lower[3]);
    fvec4 t4(table + upper[0]);
    fvec4 t5(table + upper[1]);
    fvec4 t6(table + upper[2]);
    fvec4 t7(table + upper[3]);

    // Tranposing the 8 vectors above will put all the first elements into one output
    // vector, all the second elements into the next vector and so on.
    fvec8 discard0, discard1;
    transpose(t0, t1, t2, t3, t4, t5, t6, t7, out0, out1, discard0, discard1);
}

/**
 * Given 3 vectors of floating-point data, reduce them to a single 3-element position
 * value by adding all the elements in each vector. Given inputs of:
 *   X0 X1 X2 X3 X4 X5 X6 X7
 *   Y0 Y1 Y2 Y3 Y4 Y5 Y6 Y7
 *   Z0 Z1 Z2 Z3 Z4 Z5 Z6 Z7
 * Each vector of values needs to be summed into a single value, and then stored into
 * the output vector:
 *   output[0] = (X0 + X1 + X2 + ...)
 *   output[1] = (Y0 + Y1 + Y2 + ...)
 *   output[2] = (Z0 + Z1 + Z2 + ...)
 *   output[3] = undefined
 */
static inline fvec4 reduceToVec3(fvec8 x, fvec8 y, fvec8 z) {
    // The general strategy for a vector reduce-add operation is to take values from
    // different parts of the vector and overlap them a different part of the vector and then
    // add together. Repeat this several times until all values have been summed. Initially 8
    // values can be reduced to 4, 4 to 2, and 2 to 1. The following code essentially does this
    // but exploits two things:
    //   - having multiple inputs means that some vectors can be combined together to amortise the
    //     cost of shuffling.
    //   - the output destinations are part of anther vector, so accumulate into the correct
    //     offsets to start with, instead of reducing to position 0 and re-inserting to the correct
    //     output location.
    //
    // As far as possible, accumulate x, y and z into their output positions in both the top and
    // bottom 128-bits to exploit in-lane permutes as much as possible early on.

    // Shuffle X and Z together to form one reduced vector.
    //   X2 X3 Z0 Z1 X6 X7 Z4 Z5
    const auto xzshuf = _mm256_shuffle_ps(x, z, 0b01001110);
    // Blend X and Z together to form another reduced vector, overlapping the previous.
    //   X0 X1 Z2 Z3 X4 X5 Z6 Z7
    const auto xzblend = _mm256_blend_ps(x, z, 0b11001100);
    // Add them together to form:
    // (X0 + X2) (X1 + X3) (Z0 + Z2) (Z1 + Z3) etc.
    const auto xz0 = _mm256_add_ps(xzshuf, xzblend);

    // Now there's only one vector containing all values. Shuffle again to form another overlap,
    // and then add.
    const auto xz1 = _mm256_permute_ps(xz0, 0b00110001);
    const auto xz2 = _mm256_add_ps(xz0, xz1);

    // Work on Z on its own as there's nothing else to work with. Start by permuting it to
    // form some overlaps, and then add:
    //   (Y0 + Y2) (Y1 + Y3) - - (Y4 + Y6) (Y5 + Y7) - -
    const auto yshuf = _mm256_permute_ps(y, 0b11101110);
    const auto y0 = _mm256_add_ps(yshuf, y);

    // Shift the bottom float of each pair to the right, into the correct Y location.
    const auto y1 = _mm256_permute_ps(y0, 0b00000000);
    const auto y2 = _mm256_add_ps(y0, y1);

    // Blend the results together to give a complete set of XYZ in the correct respective positions
    // of both top and bottom 128-bit lanes.
    const auto laneResult = fvec8(_mm256_blend_ps(xz2, y2, 0b00100010));

    return laneResult.lowerVec() + laneResult.upperVec();
}

#endif /*OPENMM_VECTORIZEAVX_H_*/
