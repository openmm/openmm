/**
 * Compute chain rule terms for computed values that depend explicitly on particle coordinates.
 */

KERNEL void computeGradientChainRuleTerms(GLOBAL const real4* RESTRICT posq,
    GLOBAL mm_long* RESTRICT forceBuffers
        PARAMETER_ARGUMENTS) {
    INIT_PARAM_DERIVS
    const real scale = RECIP((real) 0x100000000);
    for (int index = GLOBAL_ID; index < NUM_ATOMS; index += GLOBAL_SIZE) {
        real4 pos = posq[index];
        real3 force = make_real3(scale*forceBuffers[index], scale*forceBuffers[index+PADDED_NUM_ATOMS], scale*forceBuffers[index+PADDED_NUM_ATOMS*2]);
        COMPUTE_FORCES
        forceBuffers[index] = realToFixedPoint(force.x);
        forceBuffers[index+PADDED_NUM_ATOMS] = realToFixedPoint(force.y);
        forceBuffers[index+PADDED_NUM_ATOMS*2] = realToFixedPoint(force.z);
    }
    SAVE_PARAM_DERIVS
}
