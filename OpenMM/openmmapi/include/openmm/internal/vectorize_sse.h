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

#include <smmintrin.h>
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
    
    fvec4() {}
    fvec4(float v) : val(_mm_set1_ps(v)) {}
    fvec4(float v1, float v2, float v3, float v4) : val(_mm_set_ps(v4, v3, v2, v1)) {}
    fvec4(__m128 v) : val(v) {}
    fvec4(const float* v) : val(_mm_loadu_ps(v)) {}
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
    fvec4 operator+(const fvec4& other) const {
        return _mm_add_ps(val, other);
    }
    fvec4 operator-(const fvec4& other) const {
        return _mm_sub_ps(val, other);
    }
    fvec4 operator*(const fvec4& other) const {
        return _mm_mul_ps(val, other);
    }
    fvec4 operator/(const fvec4& other) const {
        return _mm_div_ps(val, other);
    }
    void operator+=(const fvec4& other) {
        val = _mm_add_ps(val, other);
    }
    void operator-=(const fvec4& other) {
        val = _mm_sub_ps(val, other);
    }
    void operator*=(const fvec4& other) {
        val = _mm_mul_ps(val, other);
    }
    void operator/=(const fvec4& other) {
        val = _mm_div_ps(val, other);
    }
    fvec4 operator-() const {
        return _mm_sub_ps(_mm_set1_ps(0.0f), val);
    }
    fvec4 operator&(const fvec4& other) const {
        return _mm_and_ps(val, other);
    }
    fvec4 operator|(const fvec4& other) const {
        return _mm_or_ps(val, other);
    }
    fvec4 operator==(const fvec4& other) const {
        return _mm_cmpeq_ps(val, other);
    }
    fvec4 operator!=(const fvec4& other) const {
        return _mm_cmpneq_ps(val, other);
    }
    fvec4 operator>(const fvec4& other) const {
        return _mm_cmpgt_ps(val, other);
    }
    fvec4 operator<(const fvec4& other) const {
        return _mm_cmplt_ps(val, other);
    }
    fvec4 operator>=(const fvec4& other) const {
        return _mm_cmpge_ps(val, other);
    }
    fvec4 operator<=(const fvec4& other) const {
        return _mm_cmple_ps(val, other);
    }
    operator ivec4() const;
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
    ivec4 operator+(const ivec4& other) const {
        return _mm_add_epi32(val, other);
    }
    ivec4 operator-(const ivec4& other) const {
        return _mm_sub_epi32(val, other);
    }
    ivec4 operator*(const ivec4& other) const {
        return _mm_mullo_epi32(val, other);
    }
    void operator+=(const ivec4& other) {
        val = _mm_add_epi32(val, other);
    }
    void operator-=(const ivec4& other) {
        val = _mm_sub_epi32(val, other);
    }
    void operator*=(const ivec4& other) {
        val = _mm_mullo_epi32(val, other);
    }
    ivec4 operator-() const {
        return _mm_sub_epi32(_mm_set1_epi32(0), val);
    }
    ivec4 operator&(const ivec4& other) const {
        return _mm_and_si128(val, other);
    }
    ivec4 operator|(const ivec4& other) const {
        return _mm_or_si128(val, other);
    }
    ivec4 operator==(const ivec4& other) const {
        return _mm_cmpeq_epi32(val, other);
    }
    ivec4 operator!=(const ivec4& other) const {
        return _mm_xor_si128(*this==other, _mm_set1_epi32(0xFFFFFFFF));
    }
    ivec4 operator>(const ivec4& other) const {
        return _mm_cmpgt_epi32(val, other);
    }
    ivec4 operator<(const ivec4& other) const {
        return _mm_cmplt_epi32(val, other);
    }
    ivec4 operator>=(const ivec4& other) const {
        return _mm_xor_si128(_mm_cmplt_epi32(val, other), _mm_set1_epi32(0xFFFFFFFF));
    }
    ivec4 operator<=(const ivec4& other) const {
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

// Functions that operate on fvec4s.

static inline fvec4 floor(const fvec4& v) {
    return fvec4(_mm_floor_ps(v.val));
}

static inline fvec4 ceil(const fvec4& v) {
    return fvec4(_mm_ceil_ps(v.val));
}

static inline fvec4 round(const fvec4& v) {
    return fvec4(_mm_round_ps(v.val, _MM_FROUND_TO_NEAREST_INT));
}

static inline fvec4 min(const fvec4& v1, const fvec4& v2) {
    return fvec4(_mm_min_ps(v1.val, v2.val));
}

static inline fvec4 max(const fvec4& v1, const fvec4& v2) {
    return fvec4(_mm_max_ps(v1.val, v2.val));
}

static inline fvec4 abs(const fvec4& v) {
    static const __m128 mask = _mm_castsi128_ps(_mm_set1_epi32(0x7FFFFFFF));
    return fvec4(_mm_and_ps(v.val, mask));
}

static inline fvec4 sqrt(const fvec4& v) {
    return fvec4(_mm_sqrt_ps(v.val));
}

static inline fvec4 rsqrt(const fvec4& v) {
    // Initial estimate of rsqrt().

    fvec4 y(_mm_rsqrt_ps(v.val));

    // Perform an iteration of Newton refinement.

    fvec4 x2 = v*0.5f;
    y *= fvec4(1.5f)-x2*y*y;
    return y;
}

static inline fvec4 exp(const fvec4& v) {
    return fvec4(exp_ps(v.val));
}

static inline fvec4 log(const fvec4& v) {
    return fvec4(log_ps(v.val));
}

static inline float dot3(const fvec4& v1, const fvec4& v2) {
    return _mm_cvtss_f32(_mm_dp_ps(v1, v2, 0x71));
}

static inline float dot4(const fvec4& v1, const fvec4& v2) {
    return _mm_cvtss_f32(_mm_dp_ps(v1, v2, 0xF1));
}

static inline fvec4 cross(const fvec4& v1, const fvec4& v2) {
    fvec4 temp = fvec4(_mm_mul_ps(v1, _mm_shuffle_ps(v2, v2, _MM_SHUFFLE(3, 0, 2, 1)))) -
                 fvec4(_mm_mul_ps(v2, _mm_shuffle_ps(v1, v1, _MM_SHUFFLE(3, 0, 2, 1))));
    return _mm_shuffle_ps(temp, temp, _MM_SHUFFLE(3, 0, 2, 1));
}

static inline void transpose(fvec4& v1, fvec4& v2, fvec4& v3, fvec4& v4) {
    _MM_TRANSPOSE4_PS(v1, v2, v3, v4);
}

// Functions that operate on ivec4s.

static inline ivec4 min(const ivec4& v1, const ivec4& v2) {
    return ivec4(_mm_min_epi32(v1.val, v2.val));
}

static inline ivec4 max(const ivec4& v1, const ivec4& v2) {
    return ivec4(_mm_max_epi32(v1.val, v2.val));
}

static inline ivec4 abs(const ivec4& v) {
    return ivec4(_mm_abs_epi32(v.val));
}

static inline bool any(const ivec4& v) {
    return !_mm_test_all_zeros(v, _mm_set1_epi32(0xFFFFFFFF));
}

// Mathematical operators involving a scalar and a vector.

static inline fvec4 operator+(float v1, const fvec4& v2) {
    return fvec4(v1)+v2;
}

static inline fvec4 operator-(float v1, const fvec4& v2) {
    return fvec4(v1)-v2;
}

static inline fvec4 operator*(float v1, const fvec4& v2) {
    return fvec4(v1)*v2;
}

static inline fvec4 operator/(float v1, const fvec4& v2) {
    return fvec4(v1)/v2;
}

// Operations for blending fvec4s based on an ivec4.

static inline fvec4 blend(const fvec4& v1, const fvec4& v2, const ivec4& mask) {
    return fvec4(_mm_blendv_ps(v1.val, v2.val, _mm_castsi128_ps(mask.val)));
}

#endif /*OPENMM_VECTORIZE_SSE_H_*/

