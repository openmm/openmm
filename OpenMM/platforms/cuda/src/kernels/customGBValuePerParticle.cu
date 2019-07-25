/**
 * Reduce a pairwise computed value, and compute per-particle values.
 */

extern "C" __global__ void computePerParticleValues(real4* posq, long long* valueBuffers
        PARAMETER_ARGUMENTS) {
    for (unsigned int index = blockIdx.x*blockDim.x+threadIdx.x; index < NUM_ATOMS; index += blockDim.x*gridDim.x) {
        // Load the pairwise value

        real sum = valueBuffers[index]/(real) 0x100000000;
        REDUCE_PARAM0_DERIV
        
        // Now calculate other values

        real4 pos = posq[index];
        COMPUTE_VALUES
    }
}
