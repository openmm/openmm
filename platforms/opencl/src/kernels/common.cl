/**
 * This file contains OpenCL definitions for the macros and functions needed for the
 * common compute framework.
 */

#define KERNEL __kernel
#define GLOBAL __global
#define RESTRICT restrict
#define LOCAL_ID get_local_id(0)
#define LOCAL_SIZE get_local_size(0)
#define GLOBAL_ID get_global_id(0)
#define GLOBAL_SiZE get_global_size(0)

short2 make_short2(short x, short y) {
    return (short2) (x, y);
}

short3 make_short3(short x, short y, short z) {
    return (short3) (x, y, z);
}

short4 make_short4(short x, short y, short z, short w) {
    return (short4) (x, y, z, w);
}

int2 make_int2(int x, int y) {
    return (int2) (x, y);
}

int3 make_int3(int x, int y, int z) {
    return (int3) (x, y, z);
}

int4 make_int4(int x, int y, int z, int w) {
    return (int4) (x, y, z, w);
}

float2 make_float2(float x, float y) {
    return (float2) (x, y);
}

float3 make_float3(float x, float y, float z) {
    return (float3) (x, y, z);
}

float4 make_float4(float x, float y, float z, float w) {
    return (float4) (x, y, z, w);
}

double2 make_double2(double x, double y) {
    return (double2) (x, y);
}

double3 make_double3(double x, double y, double z) {
    return (double3) (x, y, z);
}

double4 make_double4(double x, double y, double z, double w) {
    return (double4) (x, y, z, w);
}

short2 make_short2(short x) {
    return (short2) (x, y);
}

short3 make_short3(short x) {
    return (short3) (x, y, z);
}

short4 make_short4(short x) {
    return (short4) (x, y, z, w);
}

int2 make_int2(int x) {
    return (int2) x;
}

int3 make_int3(int x) {
    return (int3) x;
}

int4 make_int4(int x) {
    return (int4) x;
}

float2 make_float2(float x) {
    return (float2) x;
}

float3 make_float3(float x) {
    return (float3) x;
}

float4 make_float4(float x) {
    return (float4) x;
}

double2 make_double2(double x) {
    return (double2) x;
}

double3 make_double3(double x) {
    return (double3) x;
}

double4 make_double4(double x) {
    return (double4) x;
}
