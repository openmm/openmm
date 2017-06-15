/**
 * Copy the positions and velocities to the inner context.
 */
__kernel void copyState(__global real4* posq, __global real4* posqCorrection, __global mixed4* velm, __global int* restrict atomOrder,
        __global real4* innerPosq, __global real4* innerPosqCorrection, __global mixed4* innerVelm, __global int* restrict innerInvAtomOrder,
        int numAtoms) {
    for (int i = get_global_id(0); i < numAtoms; i += get_global_size(0)) {
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
__kernel void copyForces(__global real4* forces, __global int* restrict invAtomOrder, __global real4* innerForces,
        __global int* restrict innerAtomOrder, int numAtoms) {
    for (int i = get_global_id(0); i < numAtoms; i += get_global_size(0)) {
        int index = invAtomOrder[innerAtomOrder[i]];
        forces[index] = innerForces[i];
    }
}

/**
 * Add all the forces from the CVs.
 */
__kernel void addForces(__global real4* forces, int numAtoms
    PARAMETER_ARGUMENTS) {
    for (int i = get_global_id(0); i < numAtoms; i += get_global_size(0)) {
        real4 f = forces[i];
        ADD_FORCES
        forces[i] = f;
    }
}