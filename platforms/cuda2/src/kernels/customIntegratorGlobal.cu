extern "C" __global__ void computeGlobal(real2* __restrict__ dt, real* __restrict__ globals, real* __restrict__ params,
        float uniform, float gaussian, const real* __restrict__ energy) {
    COMPUTE_STEP
}
