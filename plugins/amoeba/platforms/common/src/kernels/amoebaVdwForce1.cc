/**
 * Clear the forces, and compute the position to use for each atom based on the bond reduction factors.
 */
KERNEL void prepareToComputeForce(GLOBAL mm_ulong* RESTRICT forceBuffers, GLOBAL real4* RESTRICT posq, GLOBAL const real4* RESTRICT tempPosq,
        GLOBAL const int* RESTRICT bondReductionAtoms, GLOBAL const float* RESTRICT bondReductionFactors) {
    for (unsigned int atom = GLOBAL_ID; atom < PADDED_NUM_ATOMS; atom += GLOBAL_SIZE) {
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
KERNEL void spreadForces(GLOBAL const mm_ulong* RESTRICT forceBuffers, GLOBAL mm_ulong* RESTRICT tempForceBuffers,
        GLOBAL const int* RESTRICT bondReductionAtoms, GLOBAL const float* RESTRICT bondReductionFactors) {
    for (unsigned int atom1 = GLOBAL_ID; atom1 < PADDED_NUM_ATOMS; atom1 += GLOBAL_SIZE) {
        int atom2 = bondReductionAtoms[atom1];
        mm_long fx1 = forceBuffers[atom1];
        mm_long fy1 = forceBuffers[atom1+PADDED_NUM_ATOMS];
        mm_long fz1 = forceBuffers[atom1+PADDED_NUM_ATOMS*2];
        if (atom1 != atom2) {
            double factor = (double) bondReductionFactors[atom1];
            mm_long fx2 = (mm_long) ((1-factor)*fx1);
            mm_long fy2 = (mm_long) ((1-factor)*fy1);
            mm_long fz2 = (mm_long) ((1-factor)*fz1);
            ATOMIC_ADD(&tempForceBuffers[atom2], (mm_ulong) fx2);
            ATOMIC_ADD(&tempForceBuffers[atom2+PADDED_NUM_ATOMS], (mm_ulong) fy2);
            ATOMIC_ADD(&tempForceBuffers[atom2+PADDED_NUM_ATOMS*2], (mm_ulong) fz2);
            fx1 = (mm_long) (factor*fx1);
            fy1 = (mm_long) (factor*fy1);
            fz1 = (mm_long) (factor*fz1);
        }
        ATOMIC_ADD(&tempForceBuffers[atom1], (mm_ulong) fx1);
        ATOMIC_ADD(&tempForceBuffers[atom1+PADDED_NUM_ATOMS], (mm_ulong) fy1);
        ATOMIC_ADD(&tempForceBuffers[atom1+PADDED_NUM_ATOMS*2], (mm_ulong) fz1);
    }
}
