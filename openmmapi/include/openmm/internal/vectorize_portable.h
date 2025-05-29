#ifndef OPENMM_VECTORIZE_PORTABLE_H_
#define OPENMM_VECTORIZE_PORTABLE_H_

/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2013-2024 Stanford University and the Authors.      *
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
    
    fvec4() = default;
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
        return val[i];
    }
    void store(float* v) const {
        *((__m128*) v) = val;
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
        return val+other.val;
    }
    fvec4 operator-(fvec4 other) const {
        return val-other.val;
    }
    fvec4 operator*(fvec4 other) const {
        return val*other.val;
    }
    fvec4 operator/(fvec4 other) const {
        return val/other.val;
    }
    void operator+=(fvec4 other) {
        val = val+other.val;
    }
    void operator-=(fvec4 other) {
        val = val-other.val;
    }
    void operator*=(fvec4 other) {
        val = val*other.val;
    }
    void operator/=(fvec4 other) {
        val = val/other.val;
    }
    fvec4 operator-() const {
        return -val;
    }
    fvec4 operator&(fvec4 other) const;
    fvec4 operator|(fvec4 other) const;
    ivec4 operator==(fvec4 other) const;
    ivec4 operator!=(fvec4 other) const;
    ivec4 operator>(fvec4 other) const;
    ivec4 operator<(fvec4 other) const;
    ivec4 operator>=(fvec4 other) const;
    ivec4 operator<=(fvec4 other) const;
    operator ivec4() const;

    /**
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
    ivec4 operator+(ivec4 other) const {
        return val+other.val;
    }
    ivec4 operator-(ivec4 other) const {
        return val-other.val;
    }
    ivec4 operator*(ivec4 other) const {
        return val*other.val;
    }
    void operator+=(ivec4 other) {
        val = val+other.val;
    }
    void operator-=(ivec4 other) {
        val = val-other.val;
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
    return (__m128i) (val==other.val);
}

inline ivec4 fvec4::operator!=(fvec4 other) const {
    return (__m128i) (val!=other.val);
}

inline ivec4 fvec4::operator>(fvec4 other) const {
    return (__m128i) (val>other.val);
}

inline ivec4 fvec4::operator<(fvec4 other) const {
    return (__m128i) (val<other.val);
}

inline ivec4 fvec4::operator>=(fvec4 other) const {
    return (__m128i) (val>=other.val);
}

inline ivec4 fvec4::operator<=(fvec4 other) const {
    return (__m128i) (val<=other.val);
}

inline fvec4::operator ivec4() const {
    return __builtin_convertvector(val, __m128i);
}

inline ivec4::operator fvec4() const {
    return __builtin_convertvector(val, __m128);
}

inline fvec4 fvec4::operator&(fvec4 other) const {
    return fvec4((__m128) (((__m128i)val)&((__m128i)other.val)));
}

inline fvec4 fvec4::operator|(fvec4 other) const {
    return fvec4((__m128) (((__m128i)val)|((__m128i)other.val)));
}

inline ivec4 fvec4::expandBitsToMask(int bitmask) {
    return ivec4(bitmask & 1 ? -1 : 0,
                 bitmask & 2 ? -1 : 0,
                 bitmask & 4 ? -1 : 0,
                 bitmask & 8 ? -1 : 0);
}

// Functions that operate on fvec4s.

static inline fvec4 abs(fvec4 v) {
    return v&(__m128) ivec4(0x7FFFFFFF).val;
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
    fvec4 temp = __builtin_shufflevector(r.val, r.val, 0, 1, -1, -1)+__builtin_shufflevector(r.val, r.val, 2, 3, -1, -1);
    return temp[0]+temp[1];
}

static inline float reduceAdd(fvec4 v) {
    return dot4(v, fvec4(1.0f));
}

static inline fvec4 cross(fvec4 v1, fvec4 v2) {
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
    return ivec4(std::min(v1[0], v2[0]), std::min(v1[1], v2[1]), std::min(v1[2], v2[2]), std::min(v1[3], v2[3]));
}

static inline ivec4 max(ivec4 v1, ivec4 v2) {
    return ivec4(std::max(v1[0], v2[0]), std::max(v1[1], v2[1]), std::max(v1[2], v2[2]), std::max(v1[3], v2[3]));
}

static inline ivec4 abs(ivec4 v) {
    return ivec4(abs(v[0]), abs(v[1]), abs(v[2]), abs(v[3]));
}

static inline bool any(__m128i v) {
    ivec4 temp = __builtin_shufflevector(v, v, 0, 1, -1, -1) | __builtin_shufflevector(v, v, 2, 3, -1, -1);
    return (temp[0] || temp[1]);
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
    return mask ? v2.val : v1.val;
}

static inline fvec4 blendZero(fvec4 v, ivec4 mask) {
    return blend(0.0f, v, mask);
}

static inline ivec4 blendZero(ivec4 v, ivec4 mask) {
    return v & mask;
}

// These are at the end since they involve other functions defined above.

static inline fvec4 min(fvec4 v1, fvec4 v2) {
    return blend(v1, v2, v1 > v2);
}

static inline fvec4 max(fvec4 v1, fvec4 v2) {
    return blend(v1, v2, v1 < v2);
}

static inline fvec4 round(fvec4 v) {
    fvec4 shift(0x1.0p23f);
    fvec4 absResult = (abs(v)+shift)-shift;
    return (__m128) ((ivec4(0x80000000).val&(__m128i)v.val) + (ivec4(0x7FFFFFFF).val&(__m128i)absResult.val));
}

static inline fvec4 floor(fvec4 v) {
    fvec4 truncated = __builtin_convertvector(__builtin_convertvector(v.val, __m128i), __m128);
    return truncated + blend(0.0f, -1.0f, truncated>v);
}

static inline fvec4 ceil(fvec4 v) {
    fvec4 truncated = __builtin_convertvector(__builtin_convertvector(v.val, __m128i), __m128);
    return truncated + blend(0.0f, 1.0f, truncated<v);
}

static inline fvec4 rsqrt(fvec4 v) {
    // Initial estimate of rsqrt().

    ivec4 i = (__m128i) v.val;
    i = ivec4(0x5f375a86)-ivec4(i.val>>ivec4(1).val);
    fvec4 y = (__m128) i.val;

    // Perform three iterations of Newton refinement.

    fvec4 x2 = 0.5f*v;
    y *= 1.5f-x2*y*y;
    y *= 1.5f-x2*y*y;
    y *= 1.5f-x2*y*y;
    return y;
}

static inline fvec4 sqrt(fvec4 v) {
    return rsqrt(v)*v;
}

/**
 * Given a table of floating-point values and a set of indexes, perform a gather read into a pair
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

#endif /*OPENMM_VECTORIZE_PORTABLE_H_*/
