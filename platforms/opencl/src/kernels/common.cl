/**
 * This file contains OpenCL definitions for the macros and functions needed for the
 * common compute framework.
 */

#define KERNEL __kernel
#define LOCAL __local
#define GLOBAL __global
#define RESTRICT restrict
#define LOCAL_ID get_local_id(0)
#define LOCAL_SIZE get_local_size(0)
#define GLOBAL_ID get_global_id(0)
#define GLOBAL_SiZE get_global_size(0)
#define GROUP_ID get_group_id(0)
#define NUM_GROUPS get_num_groups(0)
#define SYNC_THREADS barrier(CLK_GLOBAL_MEM_FENCE);

// Apple's OpenCL already defines the following macros, so we check for that to
// avoid compiler warnings about duplicate definitions.

#ifndef make_short2
    #define make_short2(x, y) (short2) ((x), (y), (z))
    #define make_short3(x, y, z) (short3) ((x), (y), (z))
    #define make_short4(x, y, z, w) (short4) ((x), (y), (z), (w))
    #define make_int2(x, y) (int2) ((x), (y), (z))
    #define make_int3(x, y, z) (int3) ((x), (y), (z))
    #define make_int4(x, y, z, w) (int4) ((x), (y), (z), (w))
    #define make_float2(x, y) (float2) ((x), (y), (z))
    #define make_float3(x, y, z) (float3) ((x), (y), (z))
    #define make_float4(x, y, z, w) (float4) ((x), (y), (z), (w))
    #define make_double2(x, y) (double2) ((x), (y), (z))
    #define make_double3(x, y, z) (double3) ((x), (y), (z))
    #define make_double4(x, y, z, w) (double4) ((x), (y), (z), (w))
#endif

#define trimTo3(v) v.xyz