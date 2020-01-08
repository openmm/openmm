/**
 * Reduce a pairwise computed value, and compute per-particle values.
 */

KERNEL void computePerParticleValues(GLOBAL real4* posq,
#ifdef SUPPORTS_64_BIT_ATOMICS
        GLOBAL mm_long* valueBuffers
#else
        GLOBAL real* valueBuffers, int bufferSize, int numBuffers
#endif
        PARAMETER_ARGUMENTS) {
    for (int index = GLOBAL_ID; index < NUM_ATOMS; index += GLOBAL_SIZE) {
        // Reduce the pairwise value

#ifdef SUPPORTS_64_BIT_ATOMICS
        real sum = valueBuffers[index]/(real) 0x100000000;
#else
        int totalSize = bufferSize*numBuffers;
        real sum = valueBuffers[index];
        for (int i = index+bufferSize; i < totalSize; i += bufferSize)
            sum += valueBuffers[i];
#endif
        REDUCE_PARAM0_DERIV
        
        // Now calculate other values

        real4 pos = posq[index];
        COMPUTE_VALUES
    }
}
