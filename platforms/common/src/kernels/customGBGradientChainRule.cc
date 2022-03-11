/**
 * Compute chain rule terms for computed values that depend explicitly on particle coordinates.
 */

KERNEL void computeGradientChainRuleTerms(GLOBAL const real4* RESTRICT posq,
#ifdef SUPPORTS_64_BIT_ATOMICS
    GLOBAL mm_long* RESTRICT forceBuffers
#else
    GLOBAL real4* RESTRICT forceBuffers
#endif
        PARAMETER_ARGUMENTS) {
    INIT_PARAM_DERIVS
    const real scale = RECIP((real) 0x100000000);
    for (int index = GLOBAL_ID; index < NUM_ATOMS; index += GLOBAL_SIZE) {
        real4 pos = posq[index];
#ifdef SUPPORTS_64_BIT_ATOMICS
        real3 force = make_real3(scale*forceBuffers[index], scale*forceBuffers[index+PADDED_NUM_ATOMS], scale*forceBuffers[index+PADDED_NUM_ATOMS*2]);
#else
        real3 force = trimTo3(forceBuffers[index]);
#endif
        COMPUTE_FORCES
#ifdef SUPPORTS_64_BIT_ATOMICS
        forceBuffers[index] = realToFixedPoint(force.x);
        forceBuffers[index+PADDED_NUM_ATOMS] = realToFixedPoint(force.y);
        forceBuffers[index+PADDED_NUM_ATOMS*2] = realToFixedPoint(force.z);
#else
        forceBuffers[index] = make_real4(force.x, force.y, force.z, 0);
#endif
    }
    SAVE_PARAM_DERIVS
}
