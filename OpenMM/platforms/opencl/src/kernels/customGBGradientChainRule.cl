/**
 * Compute chain rule terms for computed values that depend explicitly on particle coordinates.
 */

__kernel void computeGradientChainRuleTerms(__global real4* restrict forceBuffers, __global const real4* restrict posq
        PARAMETER_ARGUMENTS) {
    INIT_PARAM_DERIVS
    unsigned int index = get_global_id(0);
    while (index < NUM_ATOMS) {
        real4 pos = posq[index];
        real4 force = forceBuffers[index];
        COMPUTE_FORCES
        forceBuffers[index] = force;
        index += get_global_size(0);
    }
    SAVE_PARAM_DERIVS
}
