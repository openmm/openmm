#define REDUCE_VALUE(NAME, TYPE) {\
    TYPE sum = NAME[index]; \
    for (int i = index+bufferSize; i < totalSize; i += bufferSize) \
        sum += NAME[i]; \
    NAME[index] = sum; \
}

/**
 * Reduce the derivatives computed in the N^2 energy kernel, and compute all per-particle energy terms.
 */

KERNEL void computePerParticleEnergy(GLOBAL mixed* RESTRICT energyBuffer, GLOBAL const real4* RESTRICT posq,
        GLOBAL mm_long* RESTRICT forceBuffers
        PARAMETER_ARGUMENTS) {
    mixed energy = 0;
    INIT_PARAM_DERIVS
    for (int index = GLOBAL_ID; index < NUM_ATOMS; index += GLOBAL_SIZE) {
        // Reduce the derivatives

        REDUCE_DERIVATIVES

        // Now calculate the per-particle energy terms.

        real4 pos = posq[index];
        real3 force = make_real3(0, 0, 0);
        COMPUTE_ENERGY
    }
    energyBuffer[GLOBAL_ID] += energy;
    SAVE_PARAM_DERIVS
}
