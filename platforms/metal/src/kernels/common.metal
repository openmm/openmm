/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * This is part of the OpenMM molecular simulation toolkit.                   *
 * See https://openmm.org/development.                                        *
 *                                                                            *
 * Ported from platforms/opencl/src/kernels/common.cl.                        *
 * Original OpenCL Platform code:                                            *
 * Portions copyright (c) Stanford University and the Authors.                *
 * Authors: Peter Eastman                                                     *
 * Metal Platform code:                                                       *
 * Portions copyright (c) 2026 Chun-Chi Hung.                                 *
 * Authors: Chun-Chi Hung                                                     *
 *                                                                            *
 * This program is free software: you can redistribute it and/or modify       *
 * it under the terms of the GNU Lesser General Public License as published   *
 * by the Free Software Foundation, either version 3 of the License, or       *
 * (at your option) any later version.                                        *
 * This program is distributed WITHOUT ANY WARRANTY; without even the        *
 * implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. *
 * See the GNU Lesser General Public License for more details.                *
 * You should have received a copy of the GNU Lesser General Public License  *
 * along with this program. If not, see <http://www.gnu.org/licenses/>.       *
 * -------------------------------------------------------------------------- */

#include <metal_stdlib>
using namespace metal;

typedef float real;
typedef float2 real2;
typedef float3 real3;
typedef float4 real4;
typedef float mixed;
typedef float2 mixed2;
typedef float3 mixed3;
typedef float4 mixed4;
typedef long mm_long;
typedef ulong mm_ulong;

/** Per-thread execution values explicitly passed to helpers in Metal 3.0. */
struct MetalExecutionContext {
    uint globalId, localId, groupId, localSize, numGroups;
#if OPENMM_METAL_CHECK_FIXED_POINT_RANGE
    // [0] enables checking; [1] is a sticky per-contribution overflow flag.
    device atomic_uint* fixedPointRange;
#endif
};

#define DEVICE
#define GLOBAL device
#define LOCAL threadgroup
#define LOCAL_ARG threadgroup
#define RESTRICT
#define __global device
#define __local threadgroup
#define __private thread
#define __constant constant
#define restrict
#define GLOBAL_ID (_metal.globalId)
#define LOCAL_ID (_metal.localId)
#define GROUP_ID (_metal.groupId)
#define GLOBAL_SIZE (_metal.localSize*_metal.numGroups)
#define LOCAL_SIZE (_metal.localSize)
#define NUM_GROUPS (_metal.numGroups)
#define get_global_id(dim) GLOBAL_ID
#define get_local_id(dim) LOCAL_ID
#define get_group_id(dim) GROUP_ID
#define get_global_size(dim) GLOBAL_SIZE
#define get_local_size(dim) LOCAL_SIZE
#define get_num_groups(dim) NUM_GROUPS
#define CLK_LOCAL_MEM_FENCE mem_flags::mem_threadgroup
#define CLK_GLOBAL_MEM_FENCE mem_flags::mem_device
#define barrier(flags) threadgroup_barrier(flags)
#define SYNC_THREADS threadgroup_barrier(mem_flags::mem_threadgroup | mem_flags::mem_device);
#ifndef SYNC_WARPS
// OpenCL's SYNC_WARPS orders only the local tile exchange. Keep the MSL
// execution rendezvous; device-wide phase dependencies use later dispatches.
#define SYNC_WARPS simdgroup_barrier(mem_flags::mem_threadgroup);
#endif
// Common uses this after independent atomic reductions, never to publish a
// payload. Consumers run in later dispatches. It is not a general OpenCL fence.
#define MEM_FENCE

#define make_short2(...) short2(__VA_ARGS__)
#define make_short3(...) short3(__VA_ARGS__)
#define make_short4(...) short4(__VA_ARGS__)
#define make_int2(...) int2(__VA_ARGS__)
#define make_int3(...) int3(__VA_ARGS__)
#define make_int4(...) int4(__VA_ARGS__)
#define make_float2(...) float2(__VA_ARGS__)
#define make_float3(...) float3(__VA_ARGS__)
#define make_float4(...) float4(__VA_ARGS__)
#ifndef make_real2
#define make_real2(...) float2(__VA_ARGS__)
#define make_real3(...) float3(__VA_ARGS__)
#define make_real4(...) float4(__VA_ARGS__)
#endif
#ifndef make_mixed2
#define make_mixed2(...) float2(__VA_ARGS__)
#define make_mixed3(...) float3(__VA_ARGS__)
#define make_mixed4(...) float4(__VA_ARGS__)
#endif
#define convert_float4(v) float4(v)
#define convert_float3(v) float3(v)
#define convert_int4(v) int4(v)
#ifndef convert_real4
#define convert_real4(v) float4(v)
#endif
#ifndef convert_mixed4
#define convert_mixed4(v) float4(v)
#endif
#define trimTo3(v) (v).xyz
#define sqrtf sqrt
#define rsqrtf rsqrt
#define expf exp
#define logf log
#define powf pow
#define cosf cos
#define sinf sin
#define tanf tan
#define acosf acos
#define asinf asin
#define atanf atan
#define atan2f atan2
#define as_uint(v) as_type<uint>(v)
#define as_float(v) as_type<float>(v)

/** OpenCL cross(float4,float4) ignores w and returns zero in its fourth lane. */
inline float4 cross(float4 first, float4 second) {
    return float4(metal::cross(first.xyz, second.xyz), 0.0f);
}

/** Scoped minimization mode accumulates directly in floating buffers. */
#ifdef OPENMM_METAL_FLOAT_ACCUMULATORS
inline float realToFixedPoint(real value) {
    return value;
}
#else
#if OPENMM_METAL_CHECK_FIXED_POINT_RANGE
/** Reject unrepresentable Q32.32 contributions only during checked minimization. */
inline long metalRealToFixedPoint(MetalExecutionContext context, real value) {
    // Bitwise magnitude also catches NaN/Inf without relying on fast-math
    // floating comparisons. The conservative supported range is |value| < 2^31.
    uint magnitude = as_type<uint>(value)&0x7fffffffu;
    if (atomic_load_explicit(context.fixedPointRange, memory_order_relaxed) != 0 && magnitude >= 0x4f000000u) {
        atomic_store_explicit(context.fixedPointRange+1, 1u, memory_order_relaxed);
        return 0;
    }
    return long(value*4294967296.0f);
}
#define realToFixedPoint(value) metalRealToFixedPoint(_metal, value)
#else
/** Convert a float contribution into Common's signed Q32.32 integer. */
inline long realToFixedPoint(real value) {
    return long(value*4294967296.0f);
}
#endif
#endif

/**
 * OpenCL's two-word reduction for devices without 64-bit atomic addition.
 * Storage keeps Common's ordinary 64-bit element/component-plane layout.
 * No thread may read the combined value until the accumulation dispatch ends.
 * Initialization and later 64-bit reads are separate GPU/transfer phases.
 */
inline ulong metalAtomicAdd(device ulong* address, ulong value) {
    device atomic_uint* words = reinterpret_cast<device atomic_uint*>(address);
    uint lower = uint(value);
    uint previous = atomic_fetch_add_explicit(words, lower, memory_order_relaxed);
    uint upper = uint(value >> 32)+uint(previous > 0xffffffffu-lower);
    if (upper != 0)
        atomic_fetch_add_explicit(words+1, upper, memory_order_relaxed);
    return 0;
}

inline uint metalAtomicAdd(device uint* address, uint value) {
    return atomic_fetch_add_explicit(reinterpret_cast<device atomic_uint*>(address), value, memory_order_relaxed);
}
inline int metalAtomicAdd(device int* address, int value) {
    return atomic_fetch_add_explicit(reinterpret_cast<device atomic_int*>(address), value, memory_order_relaxed);
}
inline uint metalAtomicAdd(threadgroup uint* address, uint value) {
    return atomic_fetch_add_explicit(reinterpret_cast<threadgroup atomic_uint*>(address), value, memory_order_relaxed);
}
inline int metalAtomicAdd(threadgroup int* address, int value) {
    return atomic_fetch_add_explicit(reinterpret_cast<threadgroup atomic_int*>(address), value, memory_order_relaxed);
}
/**
 * OpenCL-style float addition through 32-bit bitwise compare/exchange.
 * Every concurrent access is atomic, including the initial read. Failed weak
 * CAS updates expected, and comparing integer bits also permits NaN payloads.
 * The optional native path changes only the atomic primitive, not storage.
 */
inline float metalAtomicAdd(device float* address, float value) {
#if OPENMM_METAL_NATIVE_FLOAT_ATOMICS
    return atomic_fetch_add_explicit(reinterpret_cast<device atomic_float*>(address), value, memory_order_relaxed);
#else
    device atomic_uint* bits = reinterpret_cast<device atomic_uint*>(address);
    uint expected = atomic_load_explicit(bits, memory_order_relaxed);
    while (true) {
        uint desired = as_type<uint>(value+as_type<float>(expected));
        if (atomic_compare_exchange_weak_explicit(bits, &expected, desired,
                memory_order_relaxed, memory_order_relaxed))
            return as_type<float>(expected);
    }
#endif
}
#define ATOMIC_ADD(address, value) metalAtomicAdd(address, value)
#define atom_add(address, value) metalAtomicAdd(address, value)
#define atomic_add(address, value) metalAtomicAdd(address, value)
#define atom_inc(address) metalAtomicAdd(address, 1u)

/**
 * Metal-local Common Compute primitives. Collective operations require every
 * participating lane to execute the same call. Broadcast/shuffle source lanes
 * must be active; rotate additionally requires the full 32-lane SIMD group
 * enforced by MetalKernel. These helpers do not select any optional algorithm.
 */
inline uint simdBallot(bool predicate) {
    return uint(simd_vote::vote_t(simd_ballot(predicate)));
}
inline bool simdAny(bool predicate) { return simd_any(predicate); }
inline bool simdAll(bool predicate) { return simd_all(predicate); }
template <typename T> inline T simdBroadcast(T value, uint lane) { return simd_broadcast(value, lane); }
template <typename T> inline T simdShuffle(T value, uint lane) { return simd_shuffle(value, lane); }
/** @brief MSL 3.0 has no native 64-bit shuffle; transport both integer words exactly. */
inline ulong simdShuffle(ulong value, uint lane) {
    uint lower = simd_shuffle(uint(value), lane);
    uint upper = simd_shuffle(uint(value>>32), lane);
    return ulong(lower)|(ulong(upper)<<32);
}
inline long simdShuffle(long value, uint lane) { return as_type<long>(simdShuffle(as_type<ulong>(value), lane)); }
template <typename T> inline T simdShuffleDown(T value, ushort delta) { return simd_shuffle_down(value, delta); }
template <typename T> inline T simdShuffleXor(T value, ushort mask) { return simd_shuffle_xor(value, mask); }
template <typename T> inline T simdRotate(T value, uint lane, uint delta) {
    return simd_shuffle(value, (lane+delta)&31u);
}
template <typename T> inline T simdReduceAdd(T value) { return simd_sum(value); }
template <typename T> inline T simdReduceMin(T value) { return simd_min(value); }
template <typename T> inline T simdReduceMax(T value) { return simd_max(value); }
template <typename T> inline T simdPrefixInclusiveAdd(T value) { return simd_prefix_inclusive_sum(value); }
template <typename T> inline T simdPrefixExclusiveAdd(T value) { return simd_prefix_exclusive_sum(value); }

/** @brief Fetch-add a signed 32-bit integer in device or threadgroup memory. */
inline int atomicAddInt32(device int* address, int value) { return metalAtomicAdd(address, value); }
inline int atomicAddInt32(threadgroup int* address, int value) { return metalAtomicAdd(address, value); }
/** @brief Fetch-add a device float using the separately selected CAS/native primitive. */
inline float atomicAddFloat32(device float* address, float value) { return metalAtomicAdd(address, value); }
#if OPENMM_METAL_HAS_UINT64_MIN_MAX
/** @brief Atomic uint64 minimum; MSL 3.0 does not return the previous value. */
inline void atomicMinUInt64(device ulong* address, ulong value) {
    atomic_min_explicit(reinterpret_cast<device atomic_ulong*>(address), value, memory_order_relaxed);
}
/** @brief Atomic uint64 maximum; MSL 3.0 does not return the previous value. */
inline void atomicMaxUInt64(device ulong* address, ulong value) {
    atomic_max_explicit(reinterpret_cast<device atomic_ulong*>(address), value, memory_order_relaxed);
}
#endif
// No native uint64 fetch-add exists in the MSL 3.0 contract. The separate
// metalAtomicAdd(ulong*) reduction above is deliberately not a fetch-add API.
inline ulong atomicAddUInt64(device ulong*, ulong) = delete;

/**
 * @brief Bounded adjacent-target reduction for the sparse-pair force helper.
 * Every active lane must pass the same buffer base. Indices may differ and
 * inactive lanes are excluded by the ballot. Q32.32 conversion happens BEFORE
 * this call: unsigned modular addition therefore preserves every integer bit.
 * Only contiguous runs of equal targets are combined. Disconnected runs use
 * separate atomic writes, avoiding a loop over every distinct target. The
 * segmented scan takes at most five rounds; unique targets skip it entirely.
 * Shuffles use the current lane whenever a predecessor is inactive or lies
 * outside the run, so no result depends on an inactive source lane.
 * This is not a fetch-add and no result may be read until the dispatch ends.
 */
inline void metalAccumulateSparseForce(MetalExecutionContext context, device ulong* buffer, uint index, ulong value) {
#if OPENMM_METAL_FAST_SPARSE_FORCE_AGGREGATION
    uint lane = context.localId&31u;
    uint active = simdBallot(true);
    uint previous = lane == 0 ? lane : lane-1;
    bool hasPrevious = lane != 0 && (active&(1u<<previous)) != 0;
    uint previousIndex = simdShuffle(index, hasPrevious ? previous : lane);
    uint heads = simdBallot(!hasPrevious || index != previousIndex);
    if (heads == active) {
        metalAtomicAdd(buffer+index, value);
        return;
    }
    // The head at or before this lane is always active, so clz never sees zero.
    uint head = 31u-clz(heads&(0xffffffffu>>(31u-lane)));
    uint distance = lane-head;
    uint maxDistance = simdReduceMax(distance);
    for (uint offset = 1; offset <= maxDistance; offset <<= 1) {
        bool inRun = distance >= offset;
        ulong preceding = simdShuffle(value, inRun ? lane-offset : lane);
        if (inRun)
            value += preceding;
    }
    // A run ends immediately before an inactive lane or a new head (or lane 32).
    uint continued = (active&~heads)>>1;
    if ((continued&(1u<<lane)) == 0)
        metalAtomicAdd(buffer+index, value);
#else
    metalAtomicAdd(buffer+index, value);
#endif
}
/** @brief Floating minimization retains its original atomic ordering and primitive. */
inline void metalAccumulateSparseForce(MetalExecutionContext context, device float* buffer, uint index, float value) {
    metalAtomicAdd(buffer+index, value);
}
#define METAL_ACCUMULATE_SPARSE_FORCE(buffer, index, value) metalAccumulateSparseForce(_metal, buffer, index, value)

/** Float complementary error function; evaluate the positive tail directly. */
inline float metal_erfc(float value) {
    float x = fabs(value);
    float t = 1.0f/(1.0f+0.5f*x);
    float tail = t*exp(-x*x-1.26551223f+t*(1.00002368f+t*(0.37409196f+
        t*(0.09678418f+t*(-0.18628806f+t*(0.27886807f+t*(-1.13520398f+
        t*(1.48851587f+t*(-0.82215223f+t*0.17087277f)))))))));
    return value < 0 ? 2.0f-tail : tail;
}

/** Avoid cancellation at zero in the corresponding error function. */
inline float metal_erf(float value) {
    if (fabs(value) < 0.5f) {
        float x2 = value*value;
        return 1.1283791670955126f*value*(1.0f+x2*(-1.0f/3.0f+
            x2*(1.0f/10.0f+x2*(-1.0f/42.0f+x2*(1.0f/216.0f+
            x2*(-1.0f/1320.0f+x2/9360.0f))))));
    }
    return value < 0 ? metal_erfc(-value)-1.0f : 1.0f-metal_erfc(value);
}
#define erf metal_erf
#define erfc metal_erfc
#define erff metal_erf
#define erfcf metal_erfc

// Match OpenCL's independently accuracy-tested native functions, not broad
// fast-math assumptions. Robust minimization and strict solvers stay precise.
#if !defined(OPENMM_METAL_FLOAT_ACCUMULATORS) && !defined(OPENMM_METAL_REQUIRE_SAFE_MATH)
#if OPENMM_METAL_USE_NATIVE_SQRT && !defined(SQRT)
#define SQRT metal::fast::sqrt
#endif
#if OPENMM_METAL_USE_NATIVE_RSQRT && !defined(RSQRT)
#define RSQRT metal::fast::rsqrt
#endif
#if OPENMM_METAL_USE_NATIVE_RECIP && !defined(RECIP)
#define RECIP(x) metal::fast::divide(1.0f, (x))
#endif
#if OPENMM_METAL_USE_NATIVE_EXP && !defined(EXP)
#define EXP metal::fast::exp
#endif
#if OPENMM_METAL_USE_NATIVE_LOG && !defined(LOG)
#define LOG metal::fast::log
#endif
#endif
#ifndef SQRT
#define SQRT sqrt
#endif
#ifndef RSQRT
#define RSQRT rsqrt
#endif
#ifndef RECIP
#define RECIP(x) (1.0f/(x))
#endif
#ifndef EXP
#define EXP exp
#endif
#ifndef LOG
#define LOG log
#endif
#ifndef POW
#define POW pow
#endif
#ifndef COS
#define COS cos
#endif
#ifndef SIN
#define SIN sin
#endif
#ifndef TAN
#define TAN tan
#endif
#ifndef ACOS
#define ACOS acos
#endif
#ifndef ASIN
#define ASIN asin
#endif
#ifndef ATAN
#define ATAN atan
#endif
#ifndef ERF
#define ERF metal_erf
#endif
#ifndef ERFC
#define ERFC metal_erfc
#endif
#ifndef FMA
#define FMA fma
#endif
#ifndef FABS
#define FABS fabs
#endif
