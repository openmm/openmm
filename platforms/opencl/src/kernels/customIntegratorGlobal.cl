__kernel void computeGlobal(__global mixed2* restrict dt, __global mixed* restrict globals, __global mixed* restrict params,
        float uniform, float gaussian, __global const real* restrict energy) {
    COMPUTE_STEP
}
