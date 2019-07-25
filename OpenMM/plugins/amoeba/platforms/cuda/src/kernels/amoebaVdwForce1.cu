/**
 * Clear the forces, and compute the position to use for each atom based on the bond reduction factors.
 */
extern "C" __global__ void prepareToComputeForce(unsigned long long* __restrict__ forceBuffers, real4* __restrict__ posq, const real4* __restrict__ tempPosq,
        const int* __restrict__ bondReductionAtoms, const float* __restrict__ bondReductionFactors) {
    for (unsigned int atom = blockIdx.x*blockDim.x+threadIdx.x; atom < PADDED_NUM_ATOMS; atom += blockDim.x*gridDim.x) {
        forceBuffers[atom] = 0;
        forceBuffers[atom+PADDED_NUM_ATOMS] = 0;
        forceBuffers[atom+PADDED_NUM_ATOMS*2] = 0;
        real4 pos1 = tempPosq[atom];
        real4 pos2 = tempPosq[bondReductionAtoms[atom]];
        real factor = (real) bondReductionFactors[atom];
        posq[atom] = make_real4(factor*pos1.x + (1-factor)*pos2.x,
                                    factor*pos1.y + (1-factor)*pos2.y,
                                    factor*pos1.z + (1-factor)*pos2.z, pos1.w);
    }
}

/**
 * Spread the forces between atoms based on the bond reduction factors.
 */
extern "C" __global__ void spreadForces(const unsigned long long* __restrict__ forceBuffers, unsigned long long* __restrict__ tempForceBuffers,
        const int* __restrict__ bondReductionAtoms, const float* __restrict__ bondReductionFactors) {
    for (unsigned int atom1 = blockIdx.x*blockDim.x+threadIdx.x; atom1 < PADDED_NUM_ATOMS; atom1 += blockDim.x*gridDim.x) {
        int atom2 = bondReductionAtoms[atom1];
        long long fx1 = forceBuffers[atom1];
        long long fy1 = forceBuffers[atom1+PADDED_NUM_ATOMS];
        long long fz1 = forceBuffers[atom1+PADDED_NUM_ATOMS*2];
        if (atom1 != atom2) {
            double factor = (double) bondReductionFactors[atom1];
            long long fx2 = (long long) ((1-factor)*fx1);
            long long fy2 = (long long) ((1-factor)*fy1);
            long long fz2 = (long long) ((1-factor)*fz1);
            atomicAdd(&tempForceBuffers[atom2], static_cast<unsigned long long>(fx2));
            atomicAdd(&tempForceBuffers[atom2+PADDED_NUM_ATOMS], static_cast<unsigned long long>(fy2));
            atomicAdd(&tempForceBuffers[atom2+PADDED_NUM_ATOMS*2], static_cast<unsigned long long>(fz2));
            fx1 = (long long) (factor*fx1);
            fy1 = (long long) (factor*fy1);
            fz1 = (long long) (factor*fz1);
        }
        atomicAdd(&tempForceBuffers[atom1], static_cast<unsigned long long>(fx1));
        atomicAdd(&tempForceBuffers[atom1+PADDED_NUM_ATOMS], static_cast<unsigned long long>(fy1));
        atomicAdd(&tempForceBuffers[atom1+PADDED_NUM_ATOMS*2], static_cast<unsigned long long>(fz1));
    }
}
