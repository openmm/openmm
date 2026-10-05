/**
 * This file contains Metal definitions for the macros and functions needed for the
 * common compute framework.
 */

#include <metal_stdlib>

using namespace metal;

uint LOCAL_ID [[thread_position_in_threadgroup]];
uint LOCAL_SIZE [[threads_per_threadgroup]];
uint GLOBAL_ID [[thread_position_in_grid]];
uint GLOBAL_SIZE [[threads_per_grid]];
uint GROUP_ID [[threadgroup_position_in_grid]];
uint NUM_GROUPS [[threadgroups_per_grid]];

#define KERNEL kernel
#define DEVICE
#define LOCAL threadgroup
#define LOCAL_ARG threadgroup
#define GLOBAL device
#define RESTRICT
#define SYNC_THREADS threadgroup_barrier(mem_flags::mem_device | mem_flags::mem_threadgroup);
#define SYNC_WARPS simdgroup_barrier(mem_flags::mem_threadgroup);
#define MEM_FENCE
#define SHFL(var, srcLane) simd_shuffle(var, srcLane)
#define BALLOT(var) ((int) (uint64_t) simd_ballot(var))

inline int ATOMIC_ADD(device int* dest, int value) {
    return atomic_fetch_add_explicit((device atomic_int*) dest, value, memory_order_relaxed);
}

inline unsigned int ATOMIC_ADD(device unsigned int* dest, unsigned int value) {
    return atomic_fetch_add_explicit((device atomic_uint*) dest, value, memory_order_relaxed);
}

inline float ATOMIC_ADD(device float* dest, float value) {
    return atomic_fetch_add_explicit((device atomic_float*) dest, value, memory_order_relaxed);
}

inline unsigned long ATOMIC_ADD(device unsigned long* dest, unsigned long value) {
    device unsigned int* word = (device unsigned int*) dest;
    int lowIndex = 0;
    unsigned int lower = value;
    unsigned int upper = value >> 32;
    unsigned int result = ATOMIC_ADD(&word[lowIndex], lower);
    int carry = (lower + (unsigned long) result >= 0x100000000 ? 1 : 0);
    upper += carry;
    if (upper != 0)
        ATOMIC_ADD(&word[1-lowIndex], upper);
    return 0;
}

typedef long mm_long;
typedef unsigned long mm_ulong;

#define make_short2(x...) short2(x)
#define make_short3(x...) short3(x)
#define make_short4(x...) short4(x)
#define make_int2(x...) int2(x)
#define make_int3(x...) int3(x)
#define make_int4(x...) int4(x)
#define make_float2(x...) float2(x)
#define make_float3(x...) float3(x)
#define make_float4(x...) float4(x)

#define trimTo3(v) (v).xyz

// Metal does not support the standard names for single precision math functions.

#define sqrtf(x) sqrt(x)
#define rsqrtf(x) rsqrt(x)
#define expf(x) exp(x)
#define logf(x) log(x)
#define powf(x) pow(x)
#define cosf(x) cos(x)
#define sinf(x) sin(x)
#define tanf(x) tan(x)
#define acosf(x) acos(x)
#define asinf(x) asin(x)
#define atanf(x) atan(x)
#define atan2f(x, y) atan2(x, y)

inline float4 cross(float4 a, float4 b) {
    return float4(cross(a.xyz, b.xyz), 0.0);
}

inline long realToFixedPoint(real x) {
    return static_cast<long>(x * 0x100000000);
}

inline float erfc(float x) {
    // This approximation for erfc is from Abramowitz and Stegun (1964) p. 299.  They cite the following as
    // the original source: C. Hastings, Jr., Approximations for Digital Computers (1955).  It has a maximum
    // error of 1.5e-7.

    float t = 1.0f/(1.0f+0.3275911f*x);
    return (0.254829592f+(-0.284496736f+(1.421413741f+(-1.453152027f+1.061405429f*t)*t)*t)*t)*t*exp(-x*x);
}

inline float erf(float x) {
    return 1.0-erfc(x);
}