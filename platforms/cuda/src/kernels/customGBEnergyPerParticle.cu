/**
 * Reduce the derivatives computed in the N^2 energy kernel, and compute all per-particle energy terms.
 */

extern "C" __global__ void computePerParticleEnergy(long long* __restrict__ forceBuffers, mixed* __restrict__ energyBuffer, const real4* __restrict__ posq
        PARAMETER_ARGUMENTS) {
    mixed energy = 0;
    INIT_PARAM_DERIVS
    for (unsigned int index = blockIdx.x*blockDim.x+threadIdx.x; index < NUM_ATOMS; index += blockDim.x*gridDim.x) {
        // Load the derivatives

        LOAD_DERIVATIVES

        // Now calculate the per-particle energy terms.

        real4 pos = posq[index];
        real3 force = make_real3(0, 0, 0);
        COMPUTE_ENERGY
    }
    energyBuffer[blockIdx.x*blockDim.x+threadIdx.x] += energy;
    SAVE_PARAM_DERIVS
}
