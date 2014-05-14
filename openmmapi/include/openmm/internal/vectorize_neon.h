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

#include <cpu-features.h>
#include <arm_neon.h>
#include <cmath>

typedef int int32_t;

// This file defines classes and functions to simplify vectorizing code with NEON.

class ivec4;

/**
 * A four element vector of floats.
 */
class fvec4 {
public:
    float32x4_t val;

    fvec4() {}
    fvec4(float v) : val(vdupq_n_f32(v)) {}
    fvec4(float v1, float v2, float v3, float v4) {
        float v[] = {v1, v2, v3, v4};
        val = vld1q_f32(v);
    }
    fvec4(float32x4_t v) : val(v) {}
    fvec4(const float* v) : val(vld1q_f32(v)) {}
    operator float32x4_t() const {
        return val;
    }
    float operator[](int i) const {
        float result[4];
        store(result);
        return result[i];
    }
    void store(float* v) const {
        vst1q_f32(v, val);
    }
    fvec4 operator+(const fvec4& other) const { // Tested OK
        return vaddq_f32(val, other);
    }
    fvec4 operator-(const fvec4& other) const { // Tested OK
        return vsubq_f32(val, other);
    }
    fvec4 operator*(const fvec4& other) const { // Tested OK
        return vmulq_f32(val, other);
    }
    fvec4 operator/(const fvec4& other) const { // Tested OK
        // NEON does not have a divide float-point operator, so we get the reciprocal and multiply.

        float32x4_t reciprocal = vrecpeq_f32(other);
        reciprocal = vmulq_f32(vrecpsq_f32(other, reciprocal), reciprocal);
        reciprocal = vmulq_f32(vrecpsq_f32(other, reciprocal), reciprocal);
        fvec4 result = vmulq_f32(val,reciprocal);
        return result;
    }
    void operator+=(const fvec4& other) {
        val = vaddq_f32(val, other);
    }
    void operator-=(const fvec4& other) {
        val = vsubq_f32(val, other);
    }
    void operator*=(const fvec4& other) {
        val = vmulq_f32(val, other);
    }
    void operator/=(const fvec4& other) {
        val = val / other.val;
    }
    fvec4 operator-() const {
        return vnegq_f32(val);
    }
    fvec4 operator&(const fvec4& other) const {
        return vreinterpretq_f32_u32(vandq_u32(vreinterpretq_u32_f32(val), vreinterpretq_u32_f32(other)));
    }
    fvec4 operator|(const fvec4& other) const {
        return vcvtq_f32_s32(vreinterpretq_s32_u32(vorrq_u32(vcvtq_u32_f32(val), vcvtq_u32_f32(other))));
    }
    fvec4 operator==(const fvec4& other) const {
        return vcvtq_f32_s32(vreinterpretq_s32_u32(vceqq_f32(val, other)));
    }
    fvec4 operator!=(const fvec4& other) const {
        return vcvtq_f32_s32(vreinterpretq_s32_u32(vmvnq_u32(vceqq_f32(val, other)))); // not(equals(val, other))
    }
    fvec4 operator>(const fvec4& other) const {
        return vcvtq_f32_s32(vreinterpretq_s32_u32(vcgtq_f32(val, other)));
    }
    fvec4 operator<(const fvec4& other) const {
        return vcvtq_f32_s32(vreinterpretq_s32_u32(vcltq_f32(val, other)));
    }
    fvec4 operator>=(const fvec4& other) const {
        return vcvtq_f32_s32(vreinterpretq_s32_u32(vcgeq_f32(val, other)));
    }
    fvec4 operator<=(const fvec4& other) const {
        return vcvtq_f32_s32(vreinterpretq_s32_u32(vcleq_f32(val, other)));
    }
    operator ivec4() const;
};

/**
 * A four element vector of ints.
 */
class ivec4 {
public:
    
    int32x4_t val;

    ivec4() {}
    ivec4(int v) : val(vdupq_n_s32(v)) {}
    ivec4(int v1, int v2, int v3, int v4) {
        int v[] = {v1, v2, v3, v4};
        val = vld1q_s32(v);
    }
    ivec4(int32x4_t v) : val(v) {}
    ivec4(const int* v) : val(vld1q_s32(v)) {}
    operator int32x4_t() const {
        return val;
    }
    int operator[](int i) const {
        int result[4];
        store(result);
        return result[i];
    }
    void store(int* v) const {
        vst1q_s32(v, val);
    }
    ivec4 operator+(const ivec4& other) const {
        return vaddq_s32(val, other);
    }
    ivec4 operator-(const ivec4& other) const {
        return vsubq_s32(val, other);
    }
    ivec4 operator*(const ivec4& other) const {
        return vmulq_s32(val, other);
    }
    void operator+=(const ivec4& other) {
        val = vaddq_s32(val, other);
    }
    void operator-=(const ivec4& other) {
        val = vsubq_s32(val, other);
    }
    void operator*=(const ivec4& other) {
        val = vmulq_s32(val, other);
    }
    ivec4 operator-() const {
        return vnegq_s32(val);
    }
    ivec4 operator&(const ivec4& other) const { // Tested OK
        return ivec4(vandq_s32(val, other));
    }
    ivec4 operator|(const ivec4& other) const {
        return ivec4(vorrq_s32(val, other));
    }
    ivec4 operator==(const ivec4& other) const {
        return ivec4(vreinterpretq_s32_u32(vceqq_s32(val, other)));
    }
    ivec4 operator!=(const ivec4& other) const { // OK
        return ivec4(vreinterpretq_s32_u32(vmvnq_u32(vceqq_s32(val, other)))); // not(equal(val, other))
    }
    ivec4 operator>(const ivec4& other) const {
        return ivec4(vreinterpretq_s32_u32(vcgtq_s32(val, other)));
    }
    ivec4 operator<(const ivec4& other) const {
        return ivec4(vreinterpretq_s32_u32(vcltq_s32(val, other)));
    }
    ivec4 operator>=(const ivec4& other) const {
        return ivec4(vreinterpretq_s32_u32(vcgeq_s32(val, other)));
    }
    ivec4 operator<=(const ivec4& other) const { // OK
        return ivec4(vreinterpretq_s32_u32(vcleq_s32(val, other)));
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

// Functions that operate on fvec4s.

static inline fvec4 floor(const fvec4& v) { // Tested: OK
    fvec4 result = v + fvec4(0.5f);
    result = (fvec4) ((ivec4) result);
    return result;
}

static inline float roundToNearest(float num) {
    return (num > 0.0f) ? std::floor(num + 0.5f) : std::ceil(num - 0.5f);
}

static inline fvec4 round(const fvec4& v) { // Tested: OK - Needs optimization
    float aux[4];
    vst1q_f32(aux, v);
    return fvec4(roundToNearest(aux[0]), roundToNearest(aux[1]), roundToNearest(aux[2]), roundToNearest(aux[3]));
}

static inline fvec4 min(const fvec4& v1, const fvec4& v2) { // Tested OK
    return fvec4(vminq_f32(v1.val, v2.val));
}

static inline fvec4 max(const fvec4& v1, const fvec4& v2) { // Tested OK
    return fvec4(vmaxq_f32(v1.val, v2.val));
}

static inline fvec4 abs(const fvec4& v) { // Tested OK
    return fvec4(vabdq_f32(v.val, fvec4(0.0)));
}

static inline fvec4 ceil(const fvec4& v) { // Tested OK
    ivec4 intVersion = (ivec4) v;
    fvec4 result = min((fvec4) (v > intVersion), fvec4(1.0f));
    result += intVersion;
    return result;
}

static inline fvec4 sqrt(const fvec4& v) {
    float32x4_t recipSqrt = vrsqrteq_f32(v);
    recipSqrt = vmulq_f32(recipSqrt, vrsqrtsq_f32(vmulq_f32(recipSqrt, v), recipSqrt));
    recipSqrt = vmulq_f32(recipSqrt, vrsqrtsq_f32(vmulq_f32(recipSqrt, v), recipSqrt));
    return vmulq_f32(v, recipSqrt);
}

static inline float dot3(const fvec4& v1, const fvec4& v2) { // Tested: OK
    fvec4 result = v1 * v2;
    float aux[4];
    vst1q_f32(aux, result);
    return aux[0] + aux[1] + aux[2]; // Ignore w component
}

static inline float dot4(const fvec4& v1, const fvec4& v2) { // Tested: OK
    fvec4 result = v1 * v2;
    float aux[4];
    vst1q_f32(aux, result);
    return aux[0] + aux[1] + aux[2] + aux[3];
}

static inline void transpose(fvec4& v1, fvec4& v2, fvec4& v3, fvec4& v4) { // Tested: OK
    float aux1[4];
    float aux2[4];
    float aux3[4];
    float aux4[4];
    vst1q_f32(aux1, v1);
    vst1q_f32(aux2, v2);
    vst1q_f32(aux3, v3);
    vst1q_f32(aux4, v4);

    v1 = fvec4(aux1[0], aux2[0], aux3[0], aux4[0]);
    v2 = fvec4(aux1[1], aux2[1], aux3[1], aux4[1]);
    v3 = fvec4(aux1[2], aux2[2], aux3[2], aux4[2]);
    v4 = fvec4(aux1[3], aux2[3], aux3[3], aux4[3]);
}

// Functions that operate on ivec4s.

static inline ivec4 min(const ivec4& v1, const ivec4& v2) { // Tested: not tested
    ivec4 res = ivec4(vminq_s32(v1.val, v2.val));
    return res;
}

static inline ivec4 max(const ivec4& v1, const ivec4& v2) { // Tested: not tested
    ivec4 res = ivec4(vmaxq_s32(v1.val, v2.val));
    return res;
}

static inline ivec4 abs(const ivec4& v) { // Tested: Not tested
    ivec4 res = ivec4(vabdq_s32(v.val, ivec4(0)));
    return res;
}

static inline bool any(const ivec4& v) { // Tested: OK
    int result[4];
    vst1q_s32(result, v);
    return result[0] != 0 || result[1] != 0 || result[2] != 0 || result[3] != 0;
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

static inline fvec4 blend(const fvec4& v1, const fvec4& v2, const ivec4& mask) { //  Tested OK
    return fvec4(vbslq_f32(vreinterpretq_u32_s32(mask.val), v2, v1));
}

#endif /*OPENMM_VECTORIZE_NEON_H_*/
