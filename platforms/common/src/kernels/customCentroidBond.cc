/**
 * Compute the center of each group.
 */
KERNEL void computeGroupCenters(int numParticleGroups, GLOBAL const real4* RESTRICT posq, GLOBAL const int* RESTRICT groupParticles,
        GLOBAL const real* RESTRICT groupWeights, GLOBAL const int* RESTRICT groupOffsets, GLOBAL real4* RESTRICT centerPositions) {
    LOCAL volatile real3 temp[64];
    for (int group = GROUP_ID; group < numParticleGroups; group += NUM_GROUPS) {
        // The threads in this block work together to compute the center one group.

        int firstIndex = groupOffsets[group];
        int lastIndex = groupOffsets[group+1];
        real3 center = make_real3(0);
        for (int index = LOCAL_ID; index < lastIndex-firstIndex; index += LOCAL_SIZE) {
            int atom = groupParticles[firstIndex+index];
            real weight = groupWeights[firstIndex+index];
            real4 pos = posq[atom];
            center.x += weight*pos.x;
            center.y += weight*pos.y;
            center.z += weight*pos.z;
        }

        // Sum the values.

        int thread = LOCAL_ID;
        temp[thread].x = center.x;
        temp[thread].y = center.y;
        temp[thread].z = center.z;
        SYNC_THREADS;
        if (thread < 32) {
            temp[thread].x += temp[thread+32].x;
            temp[thread].y += temp[thread+32].y;
            temp[thread].z += temp[thread+32].z;
        }
        SYNC_WARPS;
        if (thread < 16) {
            temp[thread].x += temp[thread+16].x;
            temp[thread].y += temp[thread+16].y;
            temp[thread].z += temp[thread+16].z;
        }
        SYNC_WARPS;
        if (thread < 8) {
            temp[thread].x += temp[thread+8].x;
            temp[thread].y += temp[thread+8].y;
            temp[thread].z += temp[thread+8].z;
        }
        SYNC_WARPS;
        if (thread < 4) {
            temp[thread].x += temp[thread+4].x;
            temp[thread].y += temp[thread+4].y;
            temp[thread].z += temp[thread+4].z;
        }
        SYNC_WARPS;
        if (thread < 2) {
            temp[thread].x += temp[thread+2].x;
            temp[thread].y += temp[thread+2].y;
            temp[thread].z += temp[thread+2].z;
        }
        SYNC_WARPS;
        if (thread == 0)
            centerPositions[group] = make_real4(temp[0].x+temp[1].x, temp[0].y+temp[1].y, temp[0].z+temp[1].z, 0);
    }
}

/**
 * Compute the forces on groups based on the bonds.
 */
KERNEL void computeGroupForces(int numParticleGroups, GLOBAL mm_ulong* RESTRICT groupForce, GLOBAL mixed* RESTRICT energyBuffer, GLOBAL const real4* RESTRICT centerPositions,
        GLOBAL const int* RESTRICT bondGroups, real4 periodicBoxSize, real4 invPeriodicBoxSize, real4 periodicBoxVecX, real4 periodicBoxVecY, real4 periodicBoxVecZ
        EXTRA_ARGS) {
    mixed energy = 0;
    INIT_PARAM_DERIVS
    for (int index = GLOBAL_ID; index < NUM_BONDS; index += GLOBAL_SIZE) {
        COMPUTE_FORCE
    }
    energyBuffer[GLOBAL_ID] += energy;
    SAVE_PARAM_DERIVS
}

/**
 * Apply the forces from the group centers to the individual atoms.
 */
KERNEL void applyForcesToAtoms(int numParticleGroups, GLOBAL const int* RESTRICT groupParticles, GLOBAL const real* RESTRICT groupWeights, GLOBAL const int* RESTRICT groupOffsets,
        GLOBAL const mm_long* RESTRICT groupForce, GLOBAL mm_ulong* RESTRICT atomForce) {
    for (int group = GROUP_ID; group < numParticleGroups; group += NUM_GROUPS) {
        mm_long fx = groupForce[group];
        mm_long fy = groupForce[group+numParticleGroups];
        mm_long fz = groupForce[group+numParticleGroups*2];
        int firstIndex = groupOffsets[group];
        int lastIndex = groupOffsets[group+1];
        for (int index = LOCAL_ID; index < lastIndex-firstIndex; index += LOCAL_SIZE) {
            int atom = groupParticles[firstIndex+index];
            real weight = groupWeights[firstIndex+index];
            ATOMIC_ADD(&atomForce[atom], (mm_ulong) ((mm_long) (fx*weight)));
            ATOMIC_ADD(&atomForce[atom+PADDED_NUM_ATOMS], (mm_ulong) ((mm_long) (fy*weight)));
            ATOMIC_ADD(&atomForce[atom+2*PADDED_NUM_ATOMS], (mm_ulong) ((mm_long) (fz*weight)));
        }
    }
}
