KERNEL void addForces(GLOBAL const real* RESTRICT forces, GLOBAL mm_long* RESTRICT forceBuffers, GLOBAL int* RESTRICT atomIndex) {
    for (int atom = GLOBAL_ID; atom < NUM_ATOMS; atom += GLOBAL_SIZE) {
        int index = atomIndex[atom];
        forceBuffers[atom] += realToFixedPoint(forces[3*index]);
        forceBuffers[atom+PADDED_NUM_ATOMS] += realToFixedPoint(forces[3*index+1]);
        forceBuffers[atom+2*PADDED_NUM_ATOMS] += realToFixedPoint(forces[3*index+2]);
    }
}
