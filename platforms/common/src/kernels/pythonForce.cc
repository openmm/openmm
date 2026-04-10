KERNEL void copyPositions(GLOBAL const real4* RESTRICT posq, GLOBAL real* RESTRICT positions, GLOBAL int* RESTRICT particles, int numParticles) {
    for (int i = GLOBAL_ID; i < numParticles; i += GLOBAL_SIZE) {
        real4 pos = posq[particles[i]];
        positions[3*i] = pos.x;
        positions[3*i+1] = pos.y;
        positions[3*i+2] = pos.z;
    }
}

KERNEL void addForcesAll(GLOBAL const real* RESTRICT forces, GLOBAL mm_long* RESTRICT forceBuffers, GLOBAL int* RESTRICT atomIndex) {
    for (int atom = GLOBAL_ID; atom < NUM_ATOMS; atom += GLOBAL_SIZE) {
        int index = atomIndex[atom];
        forceBuffers[atom] += (mm_long) (forces[3*index]*0x100000000);
        forceBuffers[atom+PADDED_NUM_ATOMS] += (mm_long) (forces[3*index+1]*0x100000000);
        forceBuffers[atom+2*PADDED_NUM_ATOMS] += (mm_long) (forces[3*index+2]*0x100000000);
    }
}

KERNEL void addForcesSubset(GLOBAL const real* RESTRICT forces, GLOBAL mm_long* RESTRICT forceBuffers, GLOBAL int* RESTRICT atomIndex, GLOBAL int* RESTRICT particles, int numParticles) {
    for (int i = GLOBAL_ID; i < numParticles; i += GLOBAL_SIZE) {
        int index = particles[i];
        forceBuffers[index] += (mm_long) (forces[3*i]*0x100000000);
        forceBuffers[index+PADDED_NUM_ATOMS] += (mm_long) (forces[3*i+1]*0x100000000);
        forceBuffers[index+2*PADDED_NUM_ATOMS] += (mm_long) (forces[3*i+2]*0x100000000);
    }
}
