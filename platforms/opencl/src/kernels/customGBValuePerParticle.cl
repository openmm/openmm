/**
 * Reduce a pairwise computed value, and compute per-particle values.
 */

__kernel void computePerParticleValues(int bufferSize, int numBuffers, __global real4* posq,
#ifdef SUPPORTS_64_BIT_ATOMICS
        __global long* valueBuffers
#else
        __global real* valueBuffers
#endif
        PARAMETER_ARGUMENTS) {
    unsigned int index = get_global_id(0);
    while (index < NUM_ATOMS) {
        // Reduce the pairwise value

#ifdef SUPPORTS_64_BIT_ATOMICS
        real sum = (1.0f/0x100000000)*valueBuffers[index];
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
        index += get_global_size(0);
    }
}
