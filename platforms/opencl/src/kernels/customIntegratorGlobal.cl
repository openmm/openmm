__kernel void computeGlobal(__global float2* restrict dt, __global float* restrict globals, __global float* restrict params,
        float uniform, float gaussian, __global const float* restrict energy) {
    COMPUTE_STEP
}
