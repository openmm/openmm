#ifndef OPENMM_VECTORIZE_SSE_H_
#define OPENMM_VECTORIZE_SSE_H_

/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2013 Stanford University and the Authors.           *
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

#ifdef __AVX__
#include <immintrin.h>
#else
#include <smmintrin.h>
#endif

#include "hardware.h"

// This file defines classes and functions to simplify vectorizing code with SSE.

// These two functions are defined in the vecmath library, which is linked into OpenMM.
__m128 exp_ps(__m128);
__m128 log_ps(__m128);

/**
 * Determine whether ivec4 and fvec4 are supported on this processor.
 */
static bool isVec4Supported() {
    int cpuInfo[4];
    cpuid(cpuInfo, 0);
    if (cpuInfo[0] >= 1) {
        cpuid(cpuInfo, 1);
        return ((cpuInfo[2] & ((int) 1 << 19)) != 0);
    }
    return false;
}

class ivec4;

/**
 * A four element vector of floats.
 */
class fvec4 {
public:
    __m128 val;
    
    fvec4() = default;
    fvec4(float v) : val(_mm_set1_ps(v)) {}
    fvec4(float v1, float v2, float v3, float v4) : val(_mm_set_ps(v4, v3, v2, v1)) {}
    fvec4(__m128 v) : val(v) {}
    fvec4(const float* v) : val(_mm_loadu_ps(v)) {}

    /**
     * Create a vector by gathering individual indexes of data from a table. Element i of the vector will
     * be loaded from table[idx[i]].
     * @param table The table from which to do a lookup.
     * @param indexes The indexes to gather.
     */
    fvec4(const float* table, const int32_t idx[4])
        : fvec4(table[idx[0]], table[idx[1]], table[idx[2]], table[idx[3]]) { }

    operator __m128() const {
        return val;
    }
    float operator[](int i) const {
        float result[4];
        store(result);
        return result[i];
    }
    void store(float* v) const {
        _mm_storeu_ps(v, val);
    }

    /**
     * Store only the lower three elements of the vector.
     */
    void storeVec3(float* v) const {
        // This code could be called from objects compiled for better SIMD domains (e.g., AVX) so conditionally
        // compile in the most efficient variant of the instruction.
#ifdef  __AVX__
        _mm_maskstore_ps(v, _mm_setr_epi32(-1, -1, -1, 0), val);
#else
        _mm_maskmoveu_si128 (_mm_castps_si128(val), _mm_setr_epi32(-1, -1, -1, 0), (char*)v);
#endif
    }

    fvec4 operator+(fvec4 other) const {
        return _mm_add_ps(val, other);
    }
    fvec4 operator-(fvec4 other) const {
        return _mm_sub_ps(val, other);
    }
    fvec4 operator*(fvec4 other) const {
        return _mm_mul_ps(val, other);
    }
    fvec4 operator/(fvec4 other) const {
        return _mm_div_ps(val, other);
    }
    void operator+=(fvec4 other) {
        val = _mm_add_ps(val, other);
    }
    void operator-=(fvec4 other) {
        val = _mm_sub_ps(val, other);
    }
    void operator*=(fvec4 other) {
        val = _mm_mul_ps(val, other);
    }
    void operator/=(fvec4 other) {
        val = _mm_div_ps(val, other);
    }
    fvec4 operator-() const {
        return _mm_sub_ps(_mm_set1_ps(0.0f), val);
    }
    fvec4 operator&(fvec4 other) const {
        return _mm_and_ps(val, other);
    }
    fvec4 operator|(fvec4 other) const {
        return _mm_or_ps(val, other);
    }
    fvec4 operator==(fvec4 other) const {
        return _mm_cmpeq_ps(val, other);
    }
    fvec4 operator!=(fvec4 other) const {
        return _mm_cmpneq_ps(val, other);
    }
    fvec4 operator>(fvec4 other) const {
        return _mm_cmpgt_ps(val, other);
    }
    fvec4 operator<(fvec4 other) const {
        return _mm_cmplt_ps(val, other);
    }
    fvec4 operator>=(fvec4 other) const {
        return _mm_cmpge_ps(val, other);
    }
    fvec4 operator<=(fvec4 other) const {
        return _mm_cmple_ps(val, other);
    }
    operator ivec4() const;

    /**
     * Convert an integer bitmask into a full vector of elements which can be used
     * by the blend function.
     */
    static fvec4 expandBitsToMask(int bitmask);
};

/**
 * A four element vector of ints.
 */
class ivec4 {
public:
    __m128i val;
    
    ivec4() {}
    ivec4(int v) : val(_mm_set1_epi32(v)) {}
    ivec4(int v1, int v2, int v3, int v4) : val(_mm_set_epi32(v4, v3, v2, v1)) {}
    ivec4(__m128i v) : val(v) {}
    ivec4(const int* v) : val(_mm_loadu_si128((const __m128i*) v)) {}
    operator __m128i() const {
        return val;
    }
    int operator[](int i) const {
        int result[4];
        store(result);
        return result[i];
    }
    void store(int* v) const {
        _mm_storeu_si128((__m128i*) v, val);
    }
    ivec4 operator+(ivec4 other) const {
        return _mm_add_epi32(val, other);
    }
    ivec4 operator-(ivec4 other) const {
        return _mm_sub_epi32(val, other);
    }
    ivec4 operator*(ivec4 other) const {
        return _mm_mullo_epi32(val, other);
    }
    void operator+=(ivec4 other) {
        val = _mm_add_epi32(val, other);
    }
    void operator-=(ivec4 other) {
        val = _mm_sub_epi32(val, other);
    }
    void operator*=(ivec4 other) {
        val = _mm_mullo_epi32(val, other);
    }
    ivec4 operator-() const {
        return _mm_sub_epi32(_mm_set1_epi32(0), val);
    }
    ivec4 operator&(ivec4 other) const {
        return _mm_and_si128(val, other);
    }
    ivec4 operator|(ivec4 other) const {
        return _mm_or_si128(val, other);
    }
    ivec4 operator==(ivec4 other) const {
        return _mm_cmpeq_epi32(val, other);
    }
    ivec4 operator!=(ivec4 other) const {
        return _mm_xor_si128(*this==other, _mm_set1_epi32(0xFFFFFFFF));
    }
    ivec4 operator>(ivec4 other) const {
        return _mm_cmpgt_epi32(val, other);
    }
    ivec4 operator<(ivec4 other) const {
        return _mm_cmplt_epi32(val, other);
    }
    ivec4 operator>=(ivec4 other) const {
        return _mm_xor_si128(_mm_cmplt_epi32(val, other), _mm_set1_epi32(0xFFFFFFFF));
    }
    ivec4 operator<=(ivec4 other) const {
        return _mm_xor_si128(_mm_cmpgt_epi32(val, other), _mm_set1_epi32(0xFFFFFFFF));
    }
    operator fvec4() const;
};

// Conversion operators.

inline fvec4::operator ivec4() const {
    return _mm_cvttps_epi32(val);
}

inline ivec4::operator fvec4() const {
    return _mm_cvtepi32_ps(val);
}

inline fvec4 fvec4::expandBitsToMask(int bitmask) {
    // Not optimal for SSE (see AVX implementation for better version)
    // but useful as an example for other SIMD architectures.
    const auto values = fvec4(bitmask & 1, bitmask & 2, bitmask & 4, bitmask & 8);
    return values != fvec4(0.0f);
}

// Functions that operate on fvec4s.

static inline fvec4 floor(fvec4 v) {
    return fvec4(_mm_floor_ps(v.val));
}

static inline fvec4 ceil(fvec4 v) {
    return fvec4(_mm_ceil_ps(v.val));
}

static inline fvec4 round(fvec4 v) {
    return fvec4(_mm_round_ps(v.val, _MM_FROUND_TO_NEAREST_INT));
}

static inline fvec4 min(fvec4 v1, fvec4 v2) {
    return fvec4(_mm_min_ps(v1.val, v2.val));
}

static inline fvec4 max(fvec4 v1, fvec4 v2) {
    return fvec4(_mm_max_ps(v1.val, v2.val));
}

static inline fvec4 abs(fvec4 v) {
    static const __m128 mask = _mm_castsi128_ps(_mm_set1_epi32(0x7FFFFFFF));
    return fvec4(_mm_and_ps(v.val, mask));
}

static inline fvec4 sqrt(fvec4 v) {
    return fvec4(_mm_sqrt_ps(v.val));
}

static inline fvec4 rsqrt(fvec4 v) {
    // Initial estimate of rsqrt().

    fvec4 y(_mm_rsqrt_ps(v.val));

    // Perform an iteration of Newton refinement.

    fvec4 x2 = v*0.5f;
    y *= fvec4(1.5f)-x2*y*y;
    return y;
}

static inline fvec4 exp(fvec4 v) {
    return fvec4(exp_ps(v.val));
}

static inline fvec4 log(fvec4 v) {
    return fvec4(log_ps(v.val));
}

static inline float dot3(fvec4 v1, fvec4 v2) {
    return _mm_cvtss_f32(_mm_dp_ps(v1, v2, 0x71));
}

static inline float dot4(fvec4 v1, fvec4 v2) {
    return _mm_cvtss_f32(_mm_dp_ps(v1, v2, 0xF1));
}

static inline float reduceAdd(fvec4 v) {
    return dot4(v, fvec4(1.0f));
}

static inline fvec4 cross(fvec4 v1, fvec4 v2) {
    fvec4 temp = fvec4(_mm_mul_ps(v1, _mm_shuffle_ps(v2, v2, _MM_SHUFFLE(3, 0, 2, 1)))) -
                 fvec4(_mm_mul_ps(v2, _mm_shuffle_ps(v1, v1, _MM_SHUFFLE(3, 0, 2, 1))));
    return _mm_shuffle_ps(temp, temp, _MM_SHUFFLE(3, 0, 2, 1));
}

static inline void transpose(fvec4& v1, fvec4& v2, fvec4& v3, fvec4& v4) {
    _MM_TRANSPOSE4_PS(v1, v2, v3, v4);
}

/**
 * Out-of-place transpose from an array into named variables.
 */
static inline void transpose(const fvec4 in[4], fvec4& v0, fvec4& v1, fvec4& v2, fvec4& v3) {
    v0 = in[0]; v1 = in[1]; v2 = in[2]; v3 = in[3];
    transpose(v0, v1, v2, v3);
}

/**
 * Out-of-place transpose from named variables into an array.
 */
static inline void transpose(fvec4 v0, fvec4 v1, fvec4 v2, fvec4 v3, fvec4 out[4]) {
    out[0] = v0; out[1] = v1; out[2] = v2; out[3] = v3;
    transpose(out[0], out[1], out[2], out[3]);
}

// Functions that operate on ivec4s.

static inline ivec4 min(ivec4 v1, ivec4 v2) {
    return ivec4(_mm_min_epi32(v1.val, v2.val));
}

static inline ivec4 max(ivec4 v1, ivec4 v2) {
    return ivec4(_mm_max_epi32(v1.val, v2.val));
}

static inline ivec4 abs(ivec4 v) {
    return ivec4(_mm_abs_epi32(v.val));
}

static inline bool any(ivec4 v) {
    return !_mm_test_all_zeros(v, _mm_set1_epi32(0xFFFFFFFF));
}

// Mathematical operators involving a scalar and a vector.

static inline fvec4 operator+(float v1, fvec4 v2) {
    return fvec4(v1)+v2;
}

static inline fvec4 operator-(float v1, fvec4 v2) {
    return fvec4(v1)-v2;
}

static inline fvec4 operator*(float v1, fvec4 v2) {
    return fvec4(v1)*v2;
}

static inline fvec4 operator/(float v1, fvec4 v2) {
    return fvec4(v1)/v2;
}

// Operations for blending fvec4
static inline fvec4 blend(fvec4 v1, fvec4 v2, fvec4 mask) {
    return fvec4(_mm_blendv_ps(v1.val, v2.val, mask.val));
}

static inline fvec4 blendZero(fvec4 v, fvec4 mask) {
    return blend(0.0f, v, mask);
}

/* Given a table of floating-point values and a set of indexes, perform a gather read into a pair
 * of vectors. The first result vector contains the values at the given indexes, and the second
 * result vector contains the values from each respective index+1.
 */
static inline void gatherVecPair(const float* table, ivec4 index, fvec4& out0, fvec4& out1) {
    fvec4 t0(table + index[0]);
    fvec4 t1(table + index[1]);
    fvec4 t2(table + index[2]);
    fvec4 t3(table + index[3]);
    transpose(t0, t1, t2, t3);
    out0 = t0;
    out1 = t1;
}

/**
 * Given 3 vectors of floating-point data, reduce them to a single 3-element position
 * value by adding all the elements in each vector. Given inputs of:
 *   X0 X1 X2 X3
 *   Y0 Y1 Y2 Y3
 *   Z0 Z1 Z2 Z3
 * Each vector of values needs to be summed into a single value, and then stored into
 * the output vector:
 *   output[0] = (X0 + X1 + X2 + X3)
 *   output[1] = (Y0 + Y1 + Y2 + Y3)
 *   output[2] = (Z0 + Z1 + Z2 + Z3)
 *   output[3] = undefined
 */
static inline fvec4 reduceToVec3(fvec4 x, fvec4 y, fvec4 z) {
    // :TODO: Could be made more efficient.
    const auto nx = reduceAdd(x);
    const auto ny = reduceAdd(y);
    const auto nz = reduceAdd(z);
    return fvec4(nx, ny, nz, 0.0);
}

#endif /*OPENMM_VECTORIZE_SSE_H_*/

