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
 * Portions copyright (c) 2013-2020 Stanford University and the Authors.      *
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
    
    fvec4() = default;
    fvec4(float v) {
        val = (__m128) {v, v, v, v};
    }
    fvec4(float v1, float v2, float v3, float v4) {
        val = (__m128) {v1, v2, v3, v4};
    }
    fvec4(__m128 v) : val(v) {}
    fvec4(const float* v) {
        val = (__m128) {v[0], v[1], v[2], v[3]};
    }

    /**
      * Create a vector by gathering individual indexes of data from a table. Element i of the vector will
      * be loaded from table[idx[i]].
      * @param table The table from which to do a lookup.
      * @param indexes The indexes to gather.
      */
    fvec4(const float* table, const int idx[4])
        : fvec4(table[idx[0]], table[idx[1]], table[idx[2]], table[idx[3]]) { }

    operator __m128() const {
        return val;
    }
    float operator[](int i) const {
        return val[i];
    }
    void store(float* v) const {
        v[0] = val[0];
        v[1] = val[1];
        v[2] = val[2];
        v[3] = val[3];
    }

    /**
      * Store only the lower three elements of the vector.
      */
    void storeVec3(float* v) const {
        v[0] = val[0];
        v[1] = val[1];
        v[2] = val[2];
    }

    fvec4 operator+(fvec4 other) const {
        return vec_add(val, other.val);
    }
    fvec4 operator-(fvec4 other) const {
        return vec_sub(val, other.val);
    }
    fvec4 operator*(fvec4 other) const {
        return vec_mul(val, other.val);
    }
    fvec4 operator/(fvec4 other) const {
        return vec_div(val, other.val);
    }
    void operator+=(fvec4 other) {
        val = vec_add(val, other.val);
    }
    void operator-=(fvec4 other) {
        val = vec_sub(val, other.val);
    }
    void operator*=(fvec4 other) {
        val = vec_mul(val, other.val);
    }
    void operator/=(fvec4 other) {
        val = vec_div(val, other.val);
    }
    fvec4 operator-() const {
        return -val;
    }
    fvec4 operator&(fvec4 other) const {
        return vec_and(val, other.val);
    }
    fvec4 operator|(fvec4 other) const {
        return vec_or(val, other.val);
    }
    ivec4 operator==(fvec4 other) const;
    ivec4 operator!=(fvec4 other) const;
    ivec4 operator>(fvec4 other) const;
    ivec4 operator<(fvec4 other) const;
    ivec4 operator>=(fvec4 other) const;
    ivec4 operator<=(fvec4 other) const;
    operator ivec4() const;

    /***
      * Convert an integer bitmask into a full vector of elements which can be used
      * by the blend function.
      */
    static ivec4 expandBitsToMask(int bitmask);

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
        val = (__m128i) {v[0], v[1], v[2], v[3]};
    }
    operator __m128i() const {
        return val;
    }
    int operator[](int i) const {
        return val[i];
    }
    void store(int* v) const {
        v[0] = val[0];
        v[1] = val[1];
        v[2] = val[2];
        v[3] = val[3];
    }
    ivec4 operator+(ivec4 other) const {
        return vec_add(val, other.val);
    }
    ivec4 operator-(ivec4 other) const {
        return vec_sub(val, other.val);
    }
    ivec4 operator*(ivec4 other) const {
        return val*other.val;
    }
    void operator+=(ivec4 other) {
        val = vec_add(val, other.val);
    }
    void operator-=(ivec4 other) {
        val = vec_sub(val, other.val);
    }
    void operator*=(ivec4 other) {
        val = val*other.val;
    }
    ivec4 operator-() const {
        return -val;
    }
    ivec4 operator&(ivec4 other) const {
        return val&other.val;
    }
    ivec4 operator|(ivec4 other) const {
        return val|other.val;
    }
    ivec4 operator==(ivec4 other) const {
        return (val==other.val);
    }
    ivec4 operator!=(ivec4 other) const {
        return (val!=other.val);
    }
    ivec4 operator>(ivec4 other) const {
        return (val>other.val);
    }
    ivec4 operator<(ivec4 other) const {
        return (val<other.val);
    }
    ivec4 operator>=(ivec4 other) const {
        return (val>=other.val);
    }
    ivec4 operator<=(ivec4 other) const {
        return (val<=other.val);
    }
    operator fvec4() const;
};

// Conversion operators.

inline ivec4 fvec4::operator==(fvec4 other) const {
    return  (val==other.val);
}

inline ivec4 fvec4::operator!=(fvec4 other) const {
    return  (val!=other.val);
}

inline ivec4 fvec4::operator>(fvec4 other) const {
    return  (val>other.val);
}

inline ivec4 fvec4::operator<(fvec4 other) const {
    return  (val<other.val);
}

inline ivec4 fvec4::operator>=(fvec4 other) const {
    return  (val>=other.val);
}

inline ivec4 fvec4::operator<=(fvec4 other) const {
    return  (val<=other.val);
}

inline fvec4::operator ivec4() const {
    return (__m128i) {(int)val[0], (int)val[1], (int)val[2], (int)val[3]};
}

inline ivec4::operator fvec4() const {
    return (__m128) {(float)val[0], (float)val[1], (float)val[2], (float)val[3]};
}

inline ivec4 fvec4::expandBitsToMask(int bitmask) {
    return ivec4(bitmask & 1 ? -1 : 0,
                 bitmask & 2 ? -1 : 0,
                 bitmask & 4 ? -1 : 0,
                 bitmask & 8 ? -1 : 0);
}

// Functions that operate on fvec4s.

static inline fvec4 abs(fvec4 v) {
    return vec_abs(v.val);
}

static inline fvec4 exp(fvec4 v) {
    return fvec4(expf(v[0]), expf(v[1]), expf(v[2]), expf(v[3]));
}

static inline fvec4 log(fvec4 v) {
    return fvec4(logf(v[0]), logf(v[1]), logf(v[2]), logf(v[3]));
}

static inline float dot3(fvec4 v1, fvec4 v2) {
    fvec4 r = v1*v2;
    return r[0]+r[1]+r[2];
}

static inline float dot4(fvec4 v1, fvec4 v2) {
    fvec4 r = v1*v2;
    fvec4 temp = r + vec_sld(r.val, r.val, 8);
    return temp[0]+temp[1];
}

static inline float reduceAdd(fvec4 v) {
    return dot4(v, fvec4(1.0f));
}

static inline fvec4 cross(fvec4 v1, fvec4 v2) {
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
    return vec_min(v1.val, v2.val);
}

static inline ivec4 max(ivec4 v1, ivec4 v2) {
    return vec_max(v1.val, v2.val);
}

static inline ivec4 abs(ivec4 v) {
    return vec_abs(v.val);
}

static inline bool any(ivec4 v) {
    return !vec_all_eq(v.val, ivec4(0).val);
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

// Operations for blending fvec4s based on an ivec4.

static inline fvec4 blend(fvec4 v1, fvec4 v2, __m128i mask) {
    return (__m128) ((mask&(__m128i)v2.val) + ((ivec4(0xFFFFFFFF)-ivec4(mask))&(__m128i)v1.val).val);
}

static inline fvec4 blendZero(fvec4 v, ivec4 mask) {
    return blend(0.0f, v, mask);
}

static inline ivec4 blendZero(ivec4 v, ivec4 mask) {
    return v & mask;
}

// These are at the end since they involve other functions defined above.

static inline fvec4 min(fvec4 v1, fvec4 v2) {
    return vec_min(v1.val, v2.val);
}

static inline fvec4 max(fvec4 v1, fvec4 v2) {
    return vec_max(v1.val, v2.val);
}

static inline fvec4 round(fvec4 v) {
    return vec_round(v.val);
}

static inline fvec4 floor(fvec4 v) {
    return vec_floor(v.val);
}

static inline fvec4 ceil(fvec4 v) {
    return vec_ceil(v.val);
}

static inline fvec4 rsqrt(fvec4 v) {
    // Initial estimate of rsqrt().

    fvec4 y(vec_rsqrte(v.val));

    // Perform an iteration of Newton refinement.

    fvec4 x2 = v*0.5f;
    y *= fvec4(1.5f)-x2*y*y;
    return y;
}

static inline fvec4 sqrt(fvec4 v) {
    return vec_sqrt(v.val);
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
    const auto nx = reduceAdd(x);
    const auto ny = reduceAdd(y);
    const auto nz = reduceAdd(z);
    return fvec4(nx, ny, nz, 0.0);
}

#endif /*OPENMM_VECTORIZE_PPC_H_*/

