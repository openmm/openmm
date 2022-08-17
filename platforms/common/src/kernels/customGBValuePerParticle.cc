/**
 * Reduce a pairwise computed value, and compute per-particle values.
 */

KERNEL void computePerParticleValues(GLOBAL real4* posq,
        GLOBAL mm_long* valueBuffers
        PARAMETER_ARGUMENTS) {
    for (int index = GLOBAL_ID; index < NUM_ATOMS; index += GLOBAL_SIZE) {
        // Reduce the pairwise value

        real sum = valueBuffers[index]/(real) 0x100000000;
        REDUCE_PARAM0_DERIV
        
        // Now calculate other values

        real4 pos = posq[index];
        COMPUTE_VALUES
    }
}
