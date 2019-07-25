/**
 * Compute chain rule terms for computed values that depend explicitly on particle coordinates.
 */

extern "C" __global__ void computeGradientChainRuleTerms(long long* __restrict__ forceBuffers, const real4* __restrict__ posq
        PARAMETER_ARGUMENTS) {
    INIT_PARAM_DERIVS
    const real scale = RECIP((real) 0x100000000);
    for (unsigned int index = blockIdx.x*blockDim.x+threadIdx.x; index < NUM_ATOMS; index += blockDim.x*gridDim.x) {
        real4 pos = posq[index];
        real3 force = make_real3(scale*forceBuffers[index], scale*forceBuffers[index+PADDED_NUM_ATOMS], scale*forceBuffers[index+PADDED_NUM_ATOMS*2]);
        COMPUTE_FORCES
        forceBuffers[index] = (long long) (force.x*0x100000000);
        forceBuffers[index+PADDED_NUM_ATOMS] = (long long) (force.y*0x100000000);
        forceBuffers[index+PADDED_NUM_ATOMS*2] = (long long) (force.z*0x100000000);
    }
    SAVE_PARAM_DERIVS
}
