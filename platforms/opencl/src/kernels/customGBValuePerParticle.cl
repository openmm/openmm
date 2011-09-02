/**
 * Reduce a pairwise computed value, and compute per-particle values.
 */

__kernel void computePerParticleValues(int bufferSize, int numBuffers, __global float4* posq,
#ifdef SUPPORTS_64_BIT_ATOMICS
        __global long* valueBuffers
#else
        __global float* valueBuffers
#endif
        PARAMETER_ARGUMENTS) {
    unsigned int index = get_global_id(0);
    while (index < NUM_ATOMS) {
        // Reduce the pairwise value

#ifdef SUPPORTS_64_BIT_ATOMICS
        float sum = (1.0f/0xFFFFFFFF)*valueBuffers[index];
#else
        int totalSize = bufferSize*numBuffers;
        float sum = valueBuffers[index];
        for (int i = index+bufferSize; i < totalSize; i += bufferSize)
            sum += valueBuffers[i];
#endif
        
        // Now calculate other values

        float4 pos = posq[index];
        COMPUTE_VALUES
        index += get_global_size(0);
    }
}
