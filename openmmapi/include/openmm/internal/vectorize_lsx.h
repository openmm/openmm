#ifndef OPENMM_VECTORIZE_LSX_H_
#define OPENMM_VECTORIZE_LSX_H_

/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2024 Stanford University and the Authors.           *
 * Authors: Zang Ruochen                                                      *
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
#include <lsxintrin.h>
#include <iostream>

// This file defines classes and functions to simplify vectorizing code with AltiVec on PPC.

/**
 * Determine whether ivec4 and fvec4 are supported on this processor.
 */
static bool isVec4Supported() {
    return true;
}

class ivec4;

/**
 * A four element vector of floats.
 */
class fvec4 {
public:
    v4f32 val;
    
    fvec4() = default;
    fvec4(float v) {
        val = (v4f32) {v, v, v, v};
    }
    fvec4(float v1, float v2, float v3, float v4) {
        val = (v4f32) {v1, v2, v3, v4};
    }
    fvec4(v4f32 v) : val(v) {}
    fvec4(const float* v) {
        val = (v4f32) {v[0], v[1], v[2], v[3]};
    }

    /**
      * Create a vector by gathering individual indexes of data from a table. Element i of the vector will
      * be loaded from table[idx[i]].
      * @param table The table from which to do a lookup.
      * @param indexes The indexes to gather.
      */
    fvec4(const float* table, const int idx[4])
        : fvec4(table[idx[0]], table[idx[1]], table[idx[2]], table[idx[3]]) { }

    operator v4f32() const {
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
        return __lsx_vfadd_s(val, other.val);
    }
    fvec4 operator-(fvec4 other) const {
        return __lsx_vfsub_s(val, other.val);
    }
    fvec4 operator*(fvec4 other) const {
        return __lsx_vfmul_s(val, other.val);
    }
    fvec4 operator/(fvec4 other) const {
        return __lsx_vfdiv_s(val, other.val);
    }
    void operator+=(fvec4 other) {
        val = __lsx_vfadd_s(val, other.val);
    }
    void operator-=(fvec4 other) {
        val = __lsx_vfsub_s(val, other.val);
    }
    void operator*=(fvec4 other) {
        val = __lsx_vfmul_s(val, other.val);
    }
    void operator/=(fvec4 other) {
        val = __lsx_vfdiv_s(val, other.val);
    }
    fvec4 operator-() const {
        return -val;
    }
    fvec4 operator&(fvec4 other) const {
        return (__m128)__lsx_vand_v((__m128i)val, (__m128i)other.val);
    }
    fvec4 operator|(fvec4 other) const {
        return (__m128)__lsx_vor_v((__m128i)val, (__m128i)other.val);
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
    v4i32 val;
    
    ivec4() {}
    ivec4(int v) {
        val = (v4i32) {v, v, v, v};
    }
    ivec4(int v1, int v2, int v3, int v4) {
        val = (v4i32) {v1, v2, v3, v4};
    }
    ivec4(v4i32 v) : val(v) {}
    ivec4(const int* v) {
        val = (v4i32) {v[0], v[1], v[2], v[3]};
    }
    operator v4i32() const {
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
        return (v4i32)__lsx_vadd_w((__m128i)val, (__m128i)other.val);
    }
    ivec4 operator-(ivec4 other) const {
        return (v4i32)__lsx_vsub_w((__m128i)val, (__m128i)other.val);
    }
    ivec4 operator*(ivec4 other) const {
        return (v4i32)__lsx_vmul_w((__m128i)val, (__m128i)other.val);
    }
    void operator+=(ivec4 other) {
        val = (v4i32)__lsx_vadd_w((__m128i)val, (__m128i)other.val);
    }
    void operator-=(ivec4 other) {
        val = (v4i32)__lsx_vsub_w((__m128i)val, (__m128i)other.val);
    }
    void operator*=(ivec4 other) {
        val = (v4i32)__lsx_vmul_w((__m128i)val, (__m128i)other.val);
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
    return (v4i32) {(int)val[0], (int)val[1], (int)val[2], (int)val[3]};
}

inline ivec4::operator fvec4() const {
    return (v4f32) {(float)val[0], (float)val[1], (float)val[2], (float)val[3]};
}

inline ivec4 fvec4::expandBitsToMask(int bitmask) {
    return ivec4(bitmask & 1 ? -1 : 0,
                 bitmask & 2 ? -1 : 0,
                 bitmask & 4 ? -1 : 0,
                 bitmask & 8 ? -1 : 0);
}

// Functions that operate on fvec4s.

static inline fvec4 abs(fvec4 v) {
    return (v4f32)__lsx_vbitclri_w((__m128i)v.val, 31);
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
    return r[0] + r[1] + r[2] + r[3];
}

static inline float reduceAdd(fvec4 v) {
    return dot4(v, fvec4(1.0f));
}

static inline fvec4 cross(fvec4 v1, fvec4 v2) {
    fvec4 temp = fvec4(__lsx_vfmul_s(v1, (__m128)__lsx_vshuf4i_w((__m128i)v2.val, 0xc9))) -
		 fvec4(__lsx_vfmul_s(v2, (__m128)__lsx_vshuf4i_w((__m128i)v1.val, 0xc9)));
    return (__m128)__lsx_vshuf4i_w((__m128i)temp.val, 0xc9);
}

static inline void transpose(fvec4& v1, fvec4& v2, fvec4& v3, fvec4& v4) {
    fvec4 T0 = (__m128)__lsx_vilvl_w((__m128i)v2.val, (__m128i)v1.val);
    fvec4 T1 = (__m128)__lsx_vilvh_w((__m128i)v2.val, (__m128i)v1.val);
    fvec4 T2 = (__m128)__lsx_vilvl_w((__m128i)v4.val, (__m128i)v3.val);
    fvec4 T3 = (__m128)__lsx_vilvh_w((__m128i)v4.val, (__m128i)v3.val);

    v1 = (__m128)__lsx_vilvl_d((__m128i)T2.val, (__m128i)T0.val);
    v2 = (__m128)__lsx_vilvh_d((__m128i)T2.val, (__m128i)T0.val);
    v3 = (__m128)__lsx_vilvl_d((__m128i)T3.val, (__m128i)T1.val);
    v4 = (__m128)__lsx_vilvh_d((__m128i)T3.val, (__m128i)T1.val);
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
    return (v4i32)__lsx_vmin_w((__m128i)v1.val, (__m128i)v2.val);
}

static inline ivec4 max(ivec4 v1, ivec4 v2) {
    return (v4i32)__lsx_vmax_w((__m128i)v1.val, (__m128i)v2.val);
}

static inline ivec4 abs(ivec4 v) {
    ivec4 zero = ivec4(0);
    return (v4i32)__lsx_vabsd_w((__m128i)v.val, (__m128i)zero.val);
}

static inline bool any(ivec4 v) {
    if (v[0] == 0 && v[1] == 0 && v[2] == 0 && v[3] == 0)
        return 0;
    else
        return 1;
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

static inline fvec4 blend(fvec4 v1, fvec4 v2, ivec4 mask) {
    return (__m128)__lsx_vbitsel_v((__m128i)v1.val, (__m128i)v2.val, (__m128i)mask.val);
}

static inline fvec4 blendZero(fvec4 v, ivec4 mask) {
    return blend(0.0f, v, mask.val);
}

static inline ivec4 blendZero(ivec4 v, ivec4 mask) {
    return v & mask;
}

// These are at the end since they involve other functions defined above.
static inline fvec4 min(fvec4 v1, fvec4 v2) {
    __m128i aNaN = __lsx_vfcmp_cun_s((__m128)v1.val, (__m128)v1.val);
    __m128i aMinOrNaN = __lsx_vor_v(__lsx_vfcmp_clt_s(v1, v2), aNaN);
    return (__m128)__lsx_vbitsel_v((__m128i)v2.val, (__m128i)v1.val, aMinOrNaN);
}

static inline fvec4 max(fvec4 v1, fvec4 v2) {
    __m128i aNaN = __lsx_vfcmp_cun_s((__m128)v1.val, (__m128)v1.val);
    __m128i aMaxOrNaN = __lsx_vor_v(__lsx_vfcmp_clt_s(v2, v1), aNaN);
    return (__m128)__lsx_vbitsel_v((__m128i)v2.val, (__m128i)v1.val, aMaxOrNaN);
}

static inline fvec4 round(fvec4 v) {
    return (__m128)__lsx_vfrint_s(v.val);
}

static inline fvec4 floor(fvec4 v) {
    return (__m128)__lsx_vfrintrm_s(v.val);
}

static inline fvec4 ceil(fvec4 v) {
    return (__m128)__lsx_vfrintrp_s(v.val); 
}

static inline fvec4 rsqrt(fvec4 v) {
    // Initial estimate of rsqrt().

    fvec4 y(__lsx_vfrsqrt_s(v.val));

    // Perform an iteration of Newton refinement.

    fvec4 x2 = v*0.5f;
    y *= fvec4(1.5f)-x2*y*y;
    return y;
}

static inline fvec4 sqrt(fvec4 v) {
    return __lsx_vfsqrt_s(v.val);
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

