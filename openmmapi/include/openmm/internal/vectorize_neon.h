#ifndef OPENMM_VECTORIZE_NEON_H_
#define OPENMM_VECTORIZE_NEON_H_

/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2013-2014 Stanford University and the Authors.      *
 * Authors: Mateus Lima, Peter Eastman                                        *
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

#ifdef __ANDROID__
#include <cpu-features.h>
#elif !defined(__APPLE__)
#include <sys/auxv.h>
#include <asm/hwcap.h>
#endif
#include <arm_neon.h>
#include <cmath>

typedef int int32_t;

// This file defines classes and functions to simplify vectorizing code with NEON.

// These two functions are defined in the vecmath library, which is linked into OpenMM.
float32x4_t exp_ps(float32x4_t);
float32x4_t log_ps(float32x4_t);

/**
 * Determine whether ivec4 and fvec4 are supported on this processor.
 */
static bool isVec4Supported() {
#ifdef __APPLE__
    return true;
#elif defined(__ANDROID__)
    uint64_t features = android_getCpuFeatures();
    return (features & ANDROID_CPU_ARM_FEATURE_NEON) != 0;
#elif defined(__ARM__)
    unsigned long features = getauxval(AT_HWCAP);
    return (features & HWCAP_NEON) != 0;
#else
    unsigned long features = getauxval(AT_HWCAP);
    return (features & HWCAP_ASIMD) != 0;
#endif
}

class ivec4;

/**
 * A four element vector of floats.
 */
class fvec4 {
public:
    float32x4_t val;

    fvec4() = default;
    fvec4(float v) : val(vdupq_n_f32(v)) {}
    fvec4(float v1, float v2, float v3, float v4) : val {v1, v2, v3, v4} {}
    fvec4(float32x4_t v) : val(v) {}
    fvec4(const float* v) : val(vld1q_f32(v)) {}
    operator float32x4_t() const {
        return val;
    }

    /**
      * Create a vector by gathering individual indexes of data from a table. Element i of the vector will
      * be loaded from table[idx[i]].
      * @param table The table from which to do a lookup.
      * @param indexes The indexes to gather.
      */
    fvec4(const float* table, const int32_t idx[4])
        : fvec4(table[idx[0]], table[idx[1]], table[idx[2]], table[idx[3]]) { }

    float operator[](int i) const {
        switch (i) {
            case 0:
                return vgetq_lane_f32(val, 0);
            case 1:
                return vgetq_lane_f32(val, 1);
            case 2:
                return vgetq_lane_f32(val, 2);
            case 3:
                return vgetq_lane_f32(val, 3);
        }
        return 0.0f;
    }
    void store(float* v) const {
        vst1q_f32(v, val);
    }

    /**
     * Store only the lower three elements of the vector.
     */
    void storeVec3(float* v) const {
        v[0] = vgetq_lane_f32(val, 0);
        v[1] = vgetq_lane_f32(val, 1);
        v[2] = vgetq_lane_f32(val, 2);
    }

    fvec4 operator+(fvec4 other) const {
        return vaddq_f32(val, other);
    }
    fvec4 operator-(fvec4 other) const {
        return vsubq_f32(val, other);
    }
    fvec4 operator*(fvec4 other) const {
        return vmulq_f32(val, other);
    }
    fvec4 operator/(fvec4 other) const {
        // NEON does not have a divide float-point operator, so we get the reciprocal and multiply.

        float32x4_t reciprocal = vrecpeq_f32(other);
        reciprocal = vmulq_f32(vrecpsq_f32(other, reciprocal), reciprocal);
        reciprocal = vmulq_f32(vrecpsq_f32(other, reciprocal), reciprocal);
        fvec4 result = vmulq_f32(val,reciprocal);
        return result;
    }
    void operator+=(fvec4 other) {
        val = vaddq_f32(val, other);
    }
    void operator-=(fvec4 other) {
        val = vsubq_f32(val, other);
    }
    void operator*=(fvec4 other) {
        val = vmulq_f32(val, other);
    }
    void operator/=(fvec4 other) {
        val = *this/other;
    }
    fvec4 operator-() const {
        return vnegq_f32(val);
    }
    fvec4 operator&(fvec4 other) const {
        return vreinterpretq_f32_u32(vandq_u32(vreinterpretq_u32_f32(val), vreinterpretq_u32_f32(other)));
    }
    fvec4 operator|(fvec4 other) const {
        return vreinterpretq_f32_u32(vorrq_u32(vreinterpretq_u32_f32(val), vreinterpretq_u32_f32(other)));
    }

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
    
    int32x4_t val;

    ivec4() {}
    ivec4(int v) : val(vdupq_n_s32(v)) {}
    ivec4(int v1, int v2, int v3, int v4) : val {v1, v2, v3, v4} {}
    ivec4(int32x4_t v) : val(v) {}
    ivec4(const int* v) : val(vld1q_s32(v)) {}
    operator int32x4_t() const {
        return val;
    }
    int operator[](int i) const {
        switch (i) {
            case 0:
                return vgetq_lane_s32(val, 0);
            case 1:
                return vgetq_lane_s32(val, 1);
            case 2:
                return vgetq_lane_s32(val, 2);
            case 3:
                return vgetq_lane_s32(val, 3);
        }
        return 0;
    }
    void store(int* v) const {
        vst1q_s32(v, val);
    }
    ivec4 operator+(ivec4 other) const {
        return vaddq_s32(val, other);
    }
    ivec4 operator-(ivec4 other) const {
        return vsubq_s32(val, other);
    }
    ivec4 operator*(ivec4 other) const {
        return vmulq_s32(val, other);
    }
    void operator+=(ivec4 other) {
        val = vaddq_s32(val, other);
    }
    void operator-=(ivec4 other) {
        val = vsubq_s32(val, other);
    }
    void operator*=(ivec4 other) {
        val = vmulq_s32(val, other);
    }
    ivec4 operator-() const {
        return vnegq_s32(val);
    }
    ivec4 operator&(ivec4 other) const {
        return vandq_s32(val, other);
    }
    ivec4 operator|(ivec4 other) const {
        return vorrq_s32(val, other);
    }
    ivec4 operator==(ivec4 other) const {
        return vreinterpretq_s32_u32(vceqq_s32(val, other));
    }
    ivec4 operator!=(ivec4 other) const {
        return vreinterpretq_s32_u32(vmvnq_u32(vceqq_s32(val, other))); // not(equal(val, other))
    }
    ivec4 operator>(ivec4 other) const {
        return vreinterpretq_s32_u32(vcgtq_s32(val, other));
    }
    ivec4 operator<(ivec4 other) const {
        return vreinterpretq_s32_u32(vcltq_s32(val, other));
    }
    ivec4 operator>=(ivec4 other) const {
        return vreinterpretq_s32_u32(vcgeq_s32(val, other));
    }
    ivec4 operator<=(ivec4 other) const {
        return vreinterpretq_s32_u32(vcleq_s32(val, other));
    }
    operator fvec4() const;
};

// Conversion operators.

inline fvec4::operator ivec4() const {
    return ivec4(vcvtq_s32_f32(val));
}

inline ivec4::operator fvec4() const {
    return fvec4(vcvtq_f32_s32(val));
}

inline ivec4 fvec4::expandBitsToMask(int bitmask) {
    return ivec4(bitmask & 1 ? -1 : 0,
                 bitmask & 2 ? -1 : 0,
                 bitmask & 4 ? -1 : 0,
                 bitmask & 8 ? -1 : 0);
}

// Comparison operators
inline ivec4 fvec4::operator==(fvec4 other) const {
    return vreinterpretq_s32_u32(vceqq_f32(val, other));
}

inline ivec4 fvec4::operator!=(fvec4 other) const {
    return vreinterpretq_s32_u32(vmvnq_u32(vceqq_f32(val, other))); // not(equals(val, other))
}

inline ivec4 fvec4::operator>(fvec4 other) const {
    return vreinterpretq_s32_u32(vcgtq_f32(val, other));
}

inline ivec4 fvec4::operator<(fvec4 other) const {
    return vreinterpretq_s32_u32(vcltq_f32(val, other));
}

inline ivec4 fvec4::operator>=(fvec4 other) const {
    return vreinterpretq_s32_u32(vcgeq_f32(val, other));
}

inline ivec4 fvec4::operator<=(fvec4 other) const {
    return vreinterpretq_s32_u32(vcleq_f32(val, other));
}

// Functions that operate on fvec4s.

static inline fvec4 min(fvec4 v1, fvec4 v2) {
    return vminq_f32(v1, v2);
}

static inline fvec4 max(fvec4 v1, fvec4 v2) {
    return vmaxq_f32(v1, v2);
}

static inline fvec4 abs(fvec4 v) {
    return vabsq_f32(v);
}

static inline fvec4 rsqrt(fvec4 v) {
    float32x4_t recipSqrt = vrsqrteq_f32(v);
    recipSqrt = vmulq_f32(recipSqrt, vrsqrtsq_f32(vmulq_f32(recipSqrt, v), recipSqrt));
    recipSqrt = vmulq_f32(recipSqrt, vrsqrtsq_f32(vmulq_f32(recipSqrt, v), recipSqrt));
    return recipSqrt;
}

static inline fvec4 sqrt(fvec4 v) {
    return rsqrt(v)*v;
}

static inline fvec4 exp(fvec4 v) {
    return fvec4(exp_ps(v.val));
}

static inline fvec4 log(fvec4 v) {
    return fvec4(log_ps(v.val));
}

static inline float dot3(fvec4 v1, fvec4 v2) {
    fvec4 result = v1*v2;
    return vgetq_lane_f32(result, 0) + vgetq_lane_f32(result, 1) + vgetq_lane_f32(result, 2);
}

static inline float dot4(fvec4 v1, fvec4 v2) {
    fvec4 result = v1*v2;
    return vgetq_lane_f32(result, 0) + vgetq_lane_f32(result, 1) + vgetq_lane_f32(result, 2) + vgetq_lane_f32(result,3);
}

static inline float reduceAdd(fvec4 v) {
#ifdef __ARM64__
    return vaddvq_f32(v);
#else
    return dot4(v, fvec4(1.0f));
#endif
}

static inline fvec4 cross(fvec4 v1, fvec4 v2) {
    return fvec4(v1[1]*v2[2] - v1[2]*v2[1],
                 v1[2]*v2[0] - v1[0]*v2[2],
                 v1[0]*v2[1] - v1[1]*v2[0], 0);
}

static inline void transpose(fvec4& v1, fvec4& v2, fvec4& v3, fvec4& v4) {
    float32x4x2_t t1 = vuzpq_f32(v1, v3);
    float32x4x2_t t2 = vuzpq_f32(v2, v4);
    float32x4x2_t t3 = vtrnq_f32(t1.val[0], t2.val[0]);
    float32x4x2_t t4 = vtrnq_f32(t1.val[1], t2.val[1]);
    v1 = t3.val[0];
    v2 = t4.val[0];
    v3 = t3.val[1];
    v4 = t4.val[1];
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
    return vminq_s32(v1, v2);
}

static inline ivec4 max(ivec4 v1, ivec4 v2) {
    return vmaxq_s32(v1, v2);
}

static inline ivec4 abs(ivec4 v) {
    return vabdq_s32(v, ivec4(0));
}

static inline bool any(ivec4 v) {
#ifdef __ARM64__
    return (vmaxvq_u32(vreinterpretq_u32_s32(v)) != 0);
#else
    return (vgetq_lane_s32(v, 0) != 0 || vgetq_lane_s32(v, 1) != 0 || vgetq_lane_s32(v, 2) != 0 || vgetq_lane_s32(v, 3) != 0);
#endif
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
    return vbslq_f32(vreinterpretq_u32_s32(mask), v2, v1);
}

static inline fvec4 blendZero(fvec4 v, ivec4 mask) {
    return vreinterpretq_f32_s32(vandq_s32(vreinterpretq_s32_f32(v), mask));
}

static inline ivec4 blendZero(ivec4 v, ivec4 mask) {
    return v & mask;
}

// These are at the end since they involve other functions defined above.

static inline fvec4 round(fvec4 v) {
    fvec4 shift(0x1.0p23f);
    fvec4 absResult = (abs(v)+shift)-shift;
    return blend(v, absResult, ivec4(0x7FFFFFFF));
}

static inline fvec4 floor(fvec4 v) {
    fvec4 rounded = round(v);
    return rounded + blend(0.0f, -1.0f, rounded>v);
}

static inline fvec4 ceil(fvec4 v) {
    fvec4 rounded = round(v);
    return rounded + blend(0.0f, 1.0f, rounded<v);
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

#endif /*OPENMM_VECTORIZE_NEON_H_*/
