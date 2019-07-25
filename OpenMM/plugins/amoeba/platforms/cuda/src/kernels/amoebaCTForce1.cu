/**
 * Clear the forces, and compute the position to use for each atom
 */
extern "C" __global__ void prepareToComputeForce(unsigned long long* __restrict__ forceBuffers, real4* __restrict__ posq, const real4* __restrict__ tempPosq) {
    for (unsigned int atom = blockIdx.x*blockDim.x+threadIdx.x; atom < PADDED_NUM_ATOMS; atom += blockDim.x*gridDim.x) {
        forceBuffers[atom] = 0;
        forceBuffers[atom+PADDED_NUM_ATOMS] = 0;
        forceBuffers[atom+PADDED_NUM_ATOMS*2] = 0;
        real4 pos1 = tempPosq[atom];
        posq[atom] = make_real4(pos1.x,
                                pos1.y,
                                pos1.z, pos1.w);
    }
}

/**
 * Spread the forces between atoms
 */
extern "C" __global__ void spreadForces(const unsigned long long* __restrict__ forceBuffers, unsigned long long* __restrict__ tempForceBuffers) {
    for (unsigned int atom1 = blockIdx.x*blockDim.x+threadIdx.x; atom1 < PADDED_NUM_ATOMS; atom1 += blockDim.x*gridDim.x) {
        long long fx1 = forceBuffers[atom1];
        long long fy1 = forceBuffers[atom1+PADDED_NUM_ATOMS];
        long long fz1 = forceBuffers[atom1+PADDED_NUM_ATOMS*2];
        atomicAdd(&tempForceBuffers[atom1], static_cast<unsigned long long>(fx1));
        atomicAdd(&tempForceBuffers[atom1+PADDED_NUM_ATOMS], static_cast<unsigned long long>(fy1));
        atomicAdd(&tempForceBuffers[atom1+PADDED_NUM_ATOMS*2], static_cast<unsigned long long>(fz1));
    }
}
