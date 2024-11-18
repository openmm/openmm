/**
 * This file contains HIP definitions for the macros and functions needed for the
 * common compute framework.
 */

#define KERNEL extern "C" __global__
#define DEVICE __device__
#define LOCAL __shared__
#define LOCAL_ARG
#define GLOBAL
#define RESTRICT __restrict__
#define LOCAL_ID threadIdx.x
#define LOCAL_SIZE blockDim.x
#define GLOBAL_ID (blockIdx.x*blockDim.x+threadIdx.x)
#define GLOBAL_SIZE (blockDim.x*gridDim.x)
#define GROUP_ID blockIdx.x
#define NUM_GROUPS gridDim.x
#define SYNC_THREADS __syncthreads();
#define SYNC_WARPS {__builtin_amdgcn_wave_barrier(); __builtin_amdgcn_fence(__ATOMIC_ACQ_REL, "wavefront");}
#define MEM_FENCE __threadfence_block();
#define ATOMIC_ADD(dest, value) atomicAdd(dest, value)

typedef long long mm_long;
typedef unsigned long long mm_ulong;

#define SUPPORTS_DOUBLE_PRECISION 1

#ifdef USE_DOUBLE_PRECISION

__device__ inline long long realToFixedPoint(double x) {
    return static_cast<long long>(x * 0x100000000);
}

#else

__device__ inline long long realToFixedPoint(float x) {
    // Faster way to calculate static_cast<long long>(x * 0x100000000) with exactly the same
    // results but less instructions.
    float integral = truncf(x);
    float fractional = (x - integral) * 0x100000000;
    unsigned int integral_u32 = static_cast<int>(integral);
    unsigned int fractional_u32 = static_cast<unsigned int>(fabsf(fractional));
    // A negative real number (with non-zero fractional) needs rounding-down x for integral and
    // changing fractional's sign. However, -1 is used as a threshold instead of 0 because, when
    // fractional is in (-1; 0], fractional_u32 is 0 and the number is considered an integer.
    bool isNegReal = fractional <= -1.0f;
    return (static_cast<unsigned long long>(isNegReal ? integral_u32 - 1 : integral_u32) << 32) |
            static_cast<unsigned long long>(isNegReal ? 0 - fractional_u32 : fractional_u32);
}

#endif
