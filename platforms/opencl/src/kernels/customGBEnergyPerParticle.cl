#define REDUCE_VALUE(NAME, TYPE) {\
    TYPE sum = NAME[index]; \
    for (int i = index+bufferSize; i < totalSize; i += bufferSize) \
        sum += NAME[i]; \
    NAME[index] = sum; \
}

/**
 * Reduce the derivatives computed in the N^2 energy kernel, and compute all per-particle energy terms.
 */

__kernel void computePerParticleEnergy(int bufferSize, int numBuffers, __global real4* restrict forceBuffers, __global mixed* restrict energyBuffer, __global const real4* restrict posq
        PARAMETER_ARGUMENTS) {
    mixed energy = 0;
    INIT_PARAM_DERIVS
    unsigned int index = get_global_id(0);
    while (index < NUM_ATOMS) {
        // Reduce the derivatives

        int totalSize = bufferSize*numBuffers;
        REDUCE_DERIVATIVES

        // Now calculate the per-particle energy terms.

        real4 pos = posq[index];
        real4 force = (real4) 0;
        COMPUTE_ENERGY
        index += get_global_size(0);
    }
    energyBuffer[get_global_id(0)] += energy;
    SAVE_PARAM_DERIVS
}
