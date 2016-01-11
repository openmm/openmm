#ifndef OPENMM_VECTORIZE_PNACL_H_
#define OPENMM_VECTORIZE_PNACL_H_

/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2013-2015 Stanford University and the Authors.      *
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

#include <cmath>
#include <cstdlib>

// This file defines classes and functions to simplify vectorizing code with portable SIMD vectors.

/**
 * Determine whether ivec4 and fvec4 are supported on this processor.
 */
static bool isVec4Supported() {
    return true;
}

typedef float __m128 __attribute__((vector_size(16), aligned(4)));
typedef int __m128i __attribute__((vector_size(16), aligned(4)));

class ivec4;

/**
 * A four element vector of floats.
 */
class fvec4 {
public:
    __m128 val;
    
    fvec4() {}
    fvec4(float v) {
        val = {v, v, v, v};
    }
    fvec4(float v1, float v2, float v3, float v4) {
        val = {v1, v2, v3, v4};
    }
    fvec4(__m128 v) : val(v) {}
    fvec4(const float* v) {
        val = *((__m128*) v);
    }
    operator __m128() const {
        return val;
    }
    float operator[](int i) const {
        return val[i];
    }
    void store(float* v) const {
        *((__m128*) v) = val;
    }
    fvec4 operator+(const fvec4& other) const {
        return val+other;
    }
    fvec4 operator-(const fvec4& other) const {
        return val-other;
    }
    fvec4 operator*(const fvec4& other) const {
        return val*other;
    }
    fvec4 operator/(const fvec4& other) const {
        return val/other;
    }
    void operator+=(const fvec4& other) {
        val = val+other;
    }
    void operator-=(const fvec4& other) {
        val = val-other;
    }
    void operator*=(const fvec4& other) {
        val = val*other;
    }
    void operator/=(const fvec4& other) {
        val = val/other;
    }
    fvec4 operator-() const {
        return -val;
    }
    fvec4 operator&(const fvec4& other) const {
        return (fvec4) (((__m128i)val)&((__m128i)other.val));
    }
    fvec4 operator|(const fvec4& other) const {
        return (fvec4) (((__m128i)val)|((__m128i)other.val));
    }
    ivec4 operator==(const fvec4& other) const;
    ivec4 operator!=(const fvec4& other) const;
    ivec4 operator>(const fvec4& other) const;
    ivec4 operator<(const fvec4& other) const;
    ivec4 operator>=(const fvec4& other) const;
    ivec4 operator<=(const fvec4& other) const;
    operator ivec4() const;
};

/**
 * A four element vector of ints.
 */
class ivec4 {
public:
    __m128i val;
    
    ivec4() {}
    ivec4(int v) {
        val = {v, v, v, v};
    }
    ivec4(int v1, int v2, int v3, int v4) {
        val = {v1, v2, v3, v4};
    }
    ivec4(__m128i v) : val(v) {}
    ivec4(const int* v) {
        val = *((__m128*) v);
    }
    operator __m128i() const {
        return val;
    }
    int operator[](int i) const {
        return val[i];
    }
    void store(int* v) const {
        *((__m128*) v) = val;
    }
    ivec4 operator+(const ivec4& other) const {
        return val+other;
    }
    ivec4 operator-(const ivec4& other) const {
        return val-other;
    }
    ivec4 operator*(const ivec4& other) const {
        return val*other;
    }
    void operator+=(const ivec4& other) {
        val = val+other;
    }
    void operator-=(const ivec4& other) {
        val = val-other;
    }
    void operator*=(const ivec4& other) {
        val = val*other;
    }
    ivec4 operator-() const {
        return -val;
    }
    ivec4 operator&(const ivec4& other) const {
        return val&other.val;
    }
    ivec4 operator|(const ivec4& other) const {
        return val|other.val;
    }
    ivec4 operator==(const ivec4& other) const {
        return (val==other.val);
    }
    ivec4 operator!=(const ivec4& other) const {
        return (val!=other.val);
    }
    ivec4 operator>(const ivec4& other) const {
        return (val>other.val);
    }
    ivec4 operator<(const ivec4& other) const {
        return (val<other.val);
    }
    ivec4 operator>=(const ivec4& other) const {
        return (val>=other.val);
    }
    ivec4 operator<=(const ivec4& other) const {
        return (val<=other.val);
    }
    operator fvec4() const;
};

// Conversion operators.

inline ivec4 fvec4::operator==(const fvec4& other) const {
    return (__m128i) (val==other.val);
}

inline ivec4 fvec4::operator!=(const fvec4& other) const {
    return (__m128i) (val!=other.val);
}

inline ivec4 fvec4::operator>(const fvec4& other) const {
    return (__m128i) (val>other.val);
}

inline ivec4 fvec4::operator<(const fvec4& other) const {
    return (__m128i) (val<other.val);
}

inline ivec4 fvec4::operator>=(const fvec4& other) const {
    return (__m128i) (val>=other.val);
}

inline ivec4 fvec4::operator<=(const fvec4& other) const {
    return (__m128i) (val<=other.val);
}

inline fvec4::operator ivec4() const {
    return __builtin_convertvector(val, __m128i);
}

inline ivec4::operator fvec4() const {
    return __builtin_convertvector(val, __m128);
}

// Functions that operate on fvec4s.

static inline fvec4 abs(const fvec4& v) {
    return v&(__m128) ivec4(0x7FFFFFFF);
}

static inline fvec4 exp(const fvec4& v) {
    return fvec4(expf(v[0]), expf(v[1]), expf(v[2]), expf(v[3]));
}

static inline fvec4 log(const fvec4& v) {
    return fvec4(logf(v[0]), logf(v[1]), logf(v[2]), logf(v[3]));
}

static inline float dot3(const fvec4& v1, const fvec4& v2) {
    fvec4 r = v1*v2;
    return r[0]+r[1]+r[2];
}

static inline float dot4(const fvec4& v1, const fvec4& v2) {
    fvec4 r = v1*v2;
    fvec4 temp = __builtin_shufflevector(r.val, r.val, 0, 1, -1, -1)+__builtin_shufflevector(r.val, r.val, 2, 3, -1, -1);
    return temp[0]+temp[1];
}

static inline fvec4 cross(const fvec4& v1, const fvec4& v2) {
    __m128 temp = v2.val*__builtin_shufflevector(v1.val, v1.val, 2, 0, 1, 3) -
                  v1.val*__builtin_shufflevector(v2.val, v2.val, 2, 0, 1, 3);
    return __builtin_shufflevector(temp, temp, 2, 0, 1, 3);
}

static inline void transpose(fvec4& v1, fvec4& v2, fvec4& v3, fvec4& v4) {
    __m128 a1 = __builtin_shufflevector(v1.val, v2.val, 0, 4, 2, 6);
    __m128 a2 = __builtin_shufflevector(v1.val, v2.val, 1, 5, 3, 7);
    __m128 a3 = __builtin_shufflevector(v3.val, v4.val, 0, 4, 2, 6);
    __m128 a4 = __builtin_shufflevector(v3.val, v4.val, 1, 5, 3, 7);
    v1 = __builtin_shufflevector(a1, a3, 0, 1, 4, 5);
    v2 = __builtin_shufflevector(a2, a4, 0, 1, 4, 5);
    v3 = __builtin_shufflevector(a1, a3, 2, 3, 6, 7);
    v4 = __builtin_shufflevector(a2, a4, 2, 3, 6, 7);
}

// Functions that operate on ivec4s.

static inline ivec4 min(const ivec4& v1, const ivec4& v2) {
    return ivec4(std::min(v1[0], v2[0]), std::min(v1[1], v2[1]), std::min(v1[2], v2[2]), std::min(v1[3], v2[3]));
}

static inline ivec4 max(const ivec4& v1, const ivec4& v2) {
    return ivec4(std::max(v1[0], v2[0]), std::max(v1[1], v2[1]), std::max(v1[2], v2[2]), std::max(v1[3], v2[3]));
}

static inline ivec4 abs(const ivec4& v) {
    return ivec4(abs(v[0]), abs(v[1]), abs(v[2]), abs(v[3]));
}

static inline bool any(const __m128i& v) {
    ivec4 temp = __builtin_shufflevector(v, v, 0, 1, -1, -1) | __builtin_shufflevector(v, v, 2, 3, -1, -1);
    return (temp[0] || temp[1]);
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

static inline fvec4 blend(const fvec4& v1, const fvec4& v2, const __m128i& mask) {
    return (__m128) ((mask&(__m128i)v2) + ((ivec4(0xFFFFFFFF)-ivec4(mask))&(__m128i)v1));
}

// These are at the end since they involve other functions defined above.

static inline fvec4 min(const fvec4& v1, const fvec4& v2) {
    return blend(v1, v2, v1 > v2);
}

static inline fvec4 max(const fvec4& v1, const fvec4& v2) {
    return blend(v1, v2, v1 < v2);
}

static inline fvec4 round(const fvec4& v) {
    fvec4 shift(0x1.0p23f);
    fvec4 absResult = (abs(v)+shift)-shift;
    return (__m128) ((ivec4(0x80000000)&(__m128i)v) + (ivec4(0x7FFFFFFF)&(__m128i)absResult));
}

static inline fvec4 floor(const fvec4& v) {
    fvec4 truncated = __builtin_convertvector(__builtin_convertvector(v.val, __m128i), __m128);
    return truncated + blend(0.0f, -1.0f, truncated>v);
}

static inline fvec4 ceil(const fvec4& v) {
    fvec4 truncated = __builtin_convertvector(__builtin_convertvector(v.val, __m128i), __m128);
    return truncated + blend(0.0f, 1.0f, truncated<v);
}

static inline fvec4 rsqrt(const fvec4& v) {
    // Initial estimate of rsqrt().

    ivec4 i = (__m128i) v;
    i = ivec4(0x5f375a86)-ivec4(i.val>>ivec4(1).val);
    fvec4 y = (__m128) i;

    // Perform three iterations of Newton refinement.

    fvec4 x2 = 0.5f*v;
    y *= 1.5f-x2*y*y;
    y *= 1.5f-x2*y*y;
    y *= 1.5f-x2*y*y;
    return y;
}

static inline fvec4 sqrt(const fvec4& v) {
    return rsqrt(v)*v;
}

#endif /*OPENMM_VECTORIZE_PNACL_H_*/

