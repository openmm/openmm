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
#define PRIVATE thread
#define RESTRICT
#define SYNC_THREADS threadgroup_barrier(mem_flags::mem_device | mem_flags::mem_threadgroup);
#define SYNC_WARPS simdgroup_barrier(mem_flags::mem_threadgroup);
#define MEM_FENCE
#define BALLOT(var) ((int) (uint64_t) simd_ballot(var))

template <class T>
inline T SHFL(T var, ushort srcLane) {
    return simd_shuffle(var, srcLane);
}

template <>
inline long SHFL(long var, ushort srcLane) {
    uint low = SHFL((uint) (var & 0xFFFFFFFF), srcLane);
    int high = SHFL((int) (var >> 32), srcLane);
    return (((long) high) << 32) + (long) low;
}

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
typedef long2 mm_long2;
typedef long3 mm_long3;
typedef long4 mm_long4;
typedef unsigned long mm_ulong;
typedef ulong2 mm_ulong2;
typedef ulong3 mm_ulong3;
typedef ulong4 mm_ulong4;

short2 make_short2(short x, short y) {
    return short2(x, y);
}
short2 make_short2(short x) {
    return short2(x, x);
}
short3 make_short3(short x, short y, short z) {
    return short3(x, y, z);
}
short3 make_short3(short x) {
    return short3(x, x, x);
}
short4 make_short4(short x, short y, short z, short w) {
    return short4(x, y, z, w);
}
short4 make_short4(short x) {
    return short4(x, x, x, x);
}
int2 make_int2(int x, int y) {
    return int2(x, y);
}
int2 make_int2(int x) {
    return int2(x, x);
}
int3 make_int3(int x, int y, int z) {
    return int3(x, y, z);
}
int3 make_int3(int x) {
    return int3(x, x, x);
}
int4 make_int4(int x, int y, int z, int w) {
    return int4(x, y, z, w);
}
int4 make_int4(int x) {
    return int4(x, x, x, x);
}
float2 make_float2(float x, float y) {
    return float2(x, y);
}
float2 make_float2(float x) {
    return float2(x, x);
}
float3 make_float3(float x, float y, float z) {
    return float3(x, y, z);
}
float3 make_float3(float x) {
    return float3(x, x, x);
}
float4 make_float4(float x, float y, float z, float w) {
    return float4(x, y, z, w);
}
float4 make_float4(float x) {
    return float4(x, x, x, x);
}

inline short3 trimTo3(short4 v) {
    return make_short3(v.x, v.y, v.z);
}

inline int3 trimTo3(int4 v) {
    return make_int3(v.x, v.y, v.z);
}

inline float3 trimTo3(float4 v) {
    return make_float3(v.x, v.y, v.z);
}

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
    return static_cast<long>(clamp(x, (real) INT_MIN, (real) INT_MAX)*0x100000000);
}

inline float erfc(float x) {
    // This approximation for erfc is from Abramowitz and Stegun (1964) p. 299.  They cite the following as
    // the original source: C. Hastings, Jr., Approximations for Digital Computers (1955).  It has a maximum
    // error of 1.5e-7.

    float t = 1.0f/(1.0f+0.3275911f*abs(x));
    float y = (0.254829592f+(-0.284496736f+(1.421413741f+(-1.453152027f+1.061405429f*t)*t)*t)*t)*t*exp(-x*x);
    return (x > 0.0f ? y : 2.0f-y);
}

inline float erf(float x) {
    if (abs(x) < 0.5f) {
        // For small x, the Taylor series is more accurate than subtracting two numbers both close to 1.

        float x2 = x*x;
        return M_2_SQRTPI_F*x*(1.0f+x2*(-1.0f/3.0f+x2*(1.0f/10.0f+x2*(-1.0f/42.0f+x2*(1.0f/216.0f+x2*(-1.0f/1320.0f+x2*(1.0f/9360.0f)))))));
    }
    return 1.0-erfc(x);
}

inline float remainder(float x, float y) {
    return x - y*round(x/y);
}