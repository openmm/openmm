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
        return val * other.val; //(__m128i) {val[0]*other[0], val[1]*other[1], val[2]*other[2], val[3]*other[3]}; 
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
    return fvec4(fabs(v[0]), fabs(v[1]), fabs(v[2]), fabs(v[3]));
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
    fvec4 temp = __builtin_shuffle(r.val, r.val, (__m128i) {0, 1, -1, -1})+__builtin_shuffle(r.val, r.val, (__m128i) {2, 3, -1, -1});
    return temp[0]+temp[1];
}

static inline fvec4 cross(const fvec4& v1, const fvec4& v2) {
    __m128 temp = v2.val*__builtin_shuffle(v1.val, v1.val, (__m128i) {2, 0, 1, 3}) -
                  v1.val*__builtin_shuffle(v2.val, v2.val, (__m128i) {2, 0, 1, 3});
    return __builtin_shuffle(temp, temp, (__m128i) {2, 0, 1, 3});
}

static inline void transpose(fvec4& v1, fvec4& v2, fvec4& v3, fvec4& v4) {
    __m128 a1 = __builtin_shuffle(v1.val, v2.val, (__m128i) {0, 4, 2, 6});
    __m128 a2 = __builtin_shuffle(v1.val, v2.val, (__m128i) {1, 5, 3, 7});
    __m128 a3 = __builtin_shuffle(v3.val, v4.val, (__m128i) {0, 4, 2, 6});
    __m128 a4 = __builtin_shuffle(v3.val, v4.val, (__m128i) {1, 5, 3, 7});
    v1 = __builtin_shuffle(a1, a3, (__m128i) {0, 1, 4, 5});
    v2 = __builtin_shuffle(a2, a4, (__m128i) {0, 1, 4, 5});
    v3 = __builtin_shuffle(a1, a3, (__m128i) {2, 3, 6, 7});
    v4 = __builtin_shuffle(a2, a4, (__m128i) {2, 3, 6, 7});
}

// Functions that operate on ivec4s.

static inline ivec4 min(const ivec4& v1, const ivec4& v2) {
    return vec_min(v1.val, v2.val);
}

static inline ivec4 max(const ivec4& v1, const ivec4& v2) {
    return vec_max(v1.val, v2.val);
}

static inline ivec4 abs(const ivec4& v) {
    return ivec4(abs(v[0]), abs(v[1]), abs(v[2]), abs(v[3]));
}

static inline bool any(const __m128i& v) {
    ivec4 temp = __builtin_shuffle(v, v, (__m128i) {0, 1, -1, -1}) | __builtin_shuffle(v, v, (__m128i) {2, 3, -1, -1});
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
    return fvec4(1.0/sqrt(v[0]), 1.0/sqrt(v[1]), 1.0/sqrt(v[2]), 1.0/sqrt(v[3]));
}

static inline fvec4 sqrt(const fvec4& v) {
    return vec_sqrt(v.val);
}

#endif /*OPENMM_VECTORIZE_PPC_H_*/

