/**
 * Copy the positions and velocities to the inner context.
 */
KERNEL void copyState(GLOBAL real4* RESTRICT posq, GLOBAL real4* RESTRICT innerPosq,
#ifdef USE_MIXED_PRECISION
        GLOBAL real4* RESTRICT posqCorrection, GLOBAL real4* RESTRICT innerPosqCorrection,
#endif
        GLOBAL mixed4* RESTRICT velm, GLOBAL mixed4* RESTRICT innerVelm, GLOBAL int* RESTRICT atomOrder, GLOBAL int* RESTRICT innerInvAtomOrder, int numAtoms) {
    for (int i = GLOBAL_ID; i < numAtoms; i += GLOBAL_SIZE) {
        int index = innerInvAtomOrder[atomOrder[i]];
        innerPosq[index] = posq[i];
        innerVelm[index] = velm[i];
#ifdef USE_MIXED_PRECISION
        innerPosqCorrection[index] = posqCorrection[i];
#endif
    }
}

/**
 * Copy the forces back to the main context.
 */
KERNEL void copyForces(GLOBAL mm_long* RESTRICT forces, GLOBAL int* RESTRICT invAtomOrder, GLOBAL mm_long* RESTRICT innerForces,
        GLOBAL int* RESTRICT innerAtomOrder, int numAtoms, int paddedNumAtoms) {
    for (int i = GLOBAL_ID; i < numAtoms; i += GLOBAL_SIZE) {
        int index = invAtomOrder[innerAtomOrder[i]];
        forces[index] = innerForces[i];
        forces[index+paddedNumAtoms] = innerForces[i+paddedNumAtoms];
        forces[index+paddedNumAtoms*2] = innerForces[i+paddedNumAtoms*2];
    }
}

/**
 * Add all the forces from the CVs.
 */
KERNEL void addForces(GLOBAL mm_long* RESTRICT forces, int bufferSize
    PARAMETER_ARGUMENTS) {
    for (int i = GLOBAL_ID; i < bufferSize; i += GLOBAL_SIZE) {
        ADD_FORCES
    }
}
