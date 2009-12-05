/**
 * Reduce a pairwise computed value, and compute per-particle values.
 */

__kernel void reduceGBValue(int bufferSize, int numBuffers, __global float* valueBuffers
        PARAMETER_ARGUMENTS) {
    unsigned int index = get_global_id(0);
    while (index < NUM_ATOMS) {
        // Reduce the pairwise value

        int totalSize = bufferSize*numBuffers;
        float sum = valueBuffers[index];
        for (int i = index+bufferSize; i < totalSize; i += bufferSize)
            sum += valueBuffers[i];

        // Now calculate other values

        COMPUTE_VALUES
        index += get_global_size(0);
    }
}
