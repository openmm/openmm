extern "C" __global__ void computeGlobal(mixed2* __restrict__ dt, mixed* __restrict__ globals, mixed* __restrict__ params,
        float uniform, float gaussian, const real energy) {
    COMPUTE_STEP
}
