/**
 * Reduce a pairwise computed value, and compute per-particle values.
 */

__kernel void computePerParticleValues(int bufferSize, int numBuffers, __global float* valueBuffers, __global float4* posq
        PARAMETER_ARGUMENTS) {
    unsigned int index = get_global_id(0);
    while (index < NUM_ATOMS) {
        // Reduce the pairwise value

        int totalSize = bufferSize*numBuffers;
        float sum = valueBuffers[index];
        for (int i = index+bufferSize; i < totalSize; i += bufferSize)
            sum += valueBuffers[i];

        // Now calculate other values

        float4 pos = posq[index];
        COMPUTE_VALUES
        index += get_global_size(0);
    }
}
