#ifndef OPENMM_VECTORIZE_PPC_H_
#define OPENMM_VECTORIZE_PPC_H_

/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2013-2018 Stanford University and the Authors.      *
 * Authors: Peter Eastman, Heng Ma                                            *
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
#include <altivec.h>

// This file defines classes and functions to simplify vectorizing code with AltiVec on PPC.

/**
 * Determine whether ivec4 and fvec4 are supported on this processor.
 */
static bool isVec4Supported() {
    return true;
}

typedef vector float __m128;
typedef vector int __m128i;

class ivec4;

/**
 * A four element vector of floats.
 */
class fvec4 {
public:
    __m128 val;
    
    fvec4() {}
    fvec4(float v) {
        val = (__m128) {v, v, v, v};
    }
    fvec4(float v1, float v2, float v3, float v4) {
        val = (__m128) {v1, v2, v3, v4};
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
        return vec_add(val, other.val);
    }
    fvec4 operator-(const fvec4& other) const {
        return vec_sub(val, other.val);
    }
    fvec4 operator*(const fvec4& other) const {
        return vec_mul(val, other.val);
    }
    fvec4 operator/(const fvec4& other) const {
        return vec_div(val, other.val);
    }
    void operator+=(const fvec4& other) {
        val = vec_add(val, other.val); 
    }
    void operator-=(const fvec4& other) {
        val = vec_sub(val, other.val);
    }
    void operator*=(const fvec4& other) {
        val = vec_mul(val, other.val); 
    }
    void operator/=(const fvec4& other) {
        val = vec_div(val, other.val); 
    }
    fvec4 operator-() const {
        return -val;
    }
    fvec4 operator&(const fvec4& other) const {
        return vec_and(val, other.val);
    }
    fvec4 operator|(const fvec4& other) const {
        return vec_or(val, other.val); 
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
        val = (__m128i) {v, v, v, v};
    }
    ivec4(int v1, int v2, int v3, int v4) {
        val = (__m128i) {v1, v2, v3, v4};
    }
    ivec4(__m128i v) : val(v) {}
    ivec4(const int* v) {
        val = *((__m128i*) v);
    }
    operator __m128i() const {
        return val;
    }
    int operator[](int i) const {
        return val[i];
    }
    void store(int* v) const {
        *((__m128i*) v) = val;
    }
    ivec4 operator+(const ivec4& other) const {
        return vec_add(val, other.val);
    }
    ivec4 operator-(const ivec4& other) const {
        return vec_sub(val, other.val);
    }
    ivec4 operator*(const ivec4& other) const {
        return val*other.val;
    }
    void operator+=(const ivec4& other) {
        val = vec_add(val, other.val);
    }
    void operator-=(const ivec4& other) {
        val = vec_sub(val, other.val);
    }
    void operator*=(const ivec4& other) {
        val = val*other.val;
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
    return  (val==other.val);
}

inline ivec4 fvec4::operator!=(const fvec4& other) const {
    return  (val!=other.val);
}

inline ivec4 fvec4::operator>(const fvec4& other) const {
    return  (val>other.val);
}

inline ivec4 fvec4::operator<(const fvec4& other) const {
    return  (val<other.val);
}

inline ivec4 fvec4::operator>=(const fvec4& other) const {
    return  (val>=other.val);
}

inline ivec4 fvec4::operator<=(const fvec4& other) const {
    return  (val<=other.val);
}

inline fvec4::operator ivec4() const {
    return (__m128i) {(int)val[0], (int)val[1], (int)val[2], (int)val[3]};
}

inline ivec4::operator fvec4() const {
    return (__m128) {(float)val[0], (float)val[1], (float)val[2], (float)val[3]};
}

// Functions that operate on fvec4s.

static inline fvec4 abs(const fvec4& v) {
    return vec_abs(v.val);
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
    fvec4 temp = r + vec_sld(r.val, r.val, 8);
    return temp[0]+temp[1];
}

static inline fvec4 cross(const fvec4& v1, const fvec4& v2) {
    vector unsigned char perm = (vector unsigned char) {8, 9, 10, 11, 0, 1, 2, 3, 4, 5, 6, 7, 12, 13, 14, 15};
    __m128 temp = v2.val*vec_perm(v1.val, v1.val, perm) -
                  v1.val*vec_perm(v2.val, v2.val, perm);
    return vec_perm(temp, temp, perm);
}

static inline void transpose(fvec4& v1, fvec4& v2, fvec4& v3, fvec4& v4) {
    vector unsigned char perm1 = (vector unsigned char) {0, 1, 2, 3, 16, 17, 18, 19, 8, 9, 10, 11, 24, 25, 26, 27};
    vector unsigned char perm2 = (vector unsigned char) {4, 5, 6, 7, 20, 21, 22, 23, 12, 13, 14, 15, 28, 29, 30, 31};
    __m128 a1 = vec_perm(v1.val, v2.val, perm1);
    __m128 a2 = vec_perm(v1.val, v2.val, perm2);
    __m128 a3 = vec_perm(v3.val, v4.val, perm1);
    __m128 a4 = vec_perm(v3.val, v4.val, perm2);
    vector unsigned char perm3 = (vector unsigned char) {0, 1, 2, 3, 4, 5, 6, 7, 16, 17, 18, 19, 20, 21, 22, 23};
    vector unsigned char perm4 = (vector unsigned char) {8, 9, 10, 11, 12, 13, 14, 15, 24, 25, 26, 27, 28, 29, 30, 31};
    v1 = vec_perm(a1, a3, perm3);
    v2 = vec_perm(a2, a4, perm3);
    v3 = vec_perm(a1, a3, perm4);
    v4 = vec_perm(a2, a4, perm4);
}

// Functions that operate on ivec4s.

static inline ivec4 min(const ivec4& v1, const ivec4& v2) {
    return vec_min(v1.val, v2.val);
}

static inline ivec4 max(const ivec4& v1, const ivec4& v2) {
    return vec_max(v1.val, v2.val);
}

static inline ivec4 abs(const ivec4& v) {
    return vec_abs(v.val);
}

static inline bool any(const __m128i& v) {
    return !vec_all_eq(v, ivec4(0).val);
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
    return (__m128) ((mask&(__m128i)v2.val) + ((ivec4(0xFFFFFFFF)-ivec4(mask))&(__m128i)v1.val).val);
}

// These are at the end since they involve other functions defined above.

static inline fvec4 min(const fvec4& v1, const fvec4& v2) {
    return vec_min(v1.val, v2.val);
}

static inline fvec4 max(const fvec4& v1, const fvec4& v2) {
    return vec_max(v1.val, v2.val);
}

static inline fvec4 round(const fvec4& v) {
    return vec_round(v.val);
}

static inline fvec4 floor(const fvec4& v) {
    return vec_floor(v.val);
}

static inline fvec4 ceil(const fvec4& v) {
    return vec_ceil(v.val);
}

static inline fvec4 rsqrt(const fvec4& v) {
    // Initial estimate of rsqrt().

    fvec4 y(vec_rsqrte(v.val));

    // Perform an iteration of Newton refinement.

    fvec4 x2 = v*0.5f;
    y *= fvec4(1.5f)-x2*y*y;
    return y;
}

static inline fvec4 sqrt(const fvec4& v) {
    return vec_sqrt(v.val);
}

#endif /*OPENMM_VECTORIZE_PPC_H_*/

