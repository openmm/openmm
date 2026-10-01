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
#define BALLOT(var) simd_ballot(var)

inline int ATOMIC_ADD(device int* dest, int value) {
    return atomic_fetch_add_explicit((device atomic_int*) dest, value, memory_order_relaxed);
}

inline unsigned int ATOMIC_ADD(device unsigned int* dest, unsigned int value) {
    return atomic_fetch_add_explicit((device atomic_uint*) dest, value, memory_order_relaxed);
}

inline float ATOMIC_ADD(device float* dest, float value) {
    return atomic_fetch_add_explicit((device atomic_float*) dest, value, memory_order_relaxed);
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

inline long realToFixedPoint(real x) {
    return static_cast<long>(x * 0x100000000);
}
