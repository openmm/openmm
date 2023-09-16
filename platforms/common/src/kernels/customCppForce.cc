KERNEL void addForces(GLOBAL const real* RESTRICT forces, GLOBAL mm_long* RESTRICT forceBuffers, GLOBAL int* RESTRICT atomIndex) {
    for (int atom = GLOBAL_ID; atom < NUM_ATOMS; atom += GLOBAL_SIZE) {
        int index = atomIndex[atom];
        forceBuffers[atom] += (mm_long) (forces[3*index]*0x100000000);
        forceBuffers[atom+PADDED_NUM_ATOMS] += (mm_long) (forces[3*index+1]*0x100000000);
        forceBuffers[atom+2*PADDED_NUM_ATOMS] += (mm_long) (forces[3*index+2]*0x100000000);
    }
}
