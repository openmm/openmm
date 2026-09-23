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
        forceBuffers[atom] += realToFixedPoint(forces[3*index]);
        forceBuffers[atom+PADDED_NUM_ATOMS] += realToFixedPoint(forces[3*index+1]);
        forceBuffers[atom+2*PADDED_NUM_ATOMS] += realToFixedPoint(forces[3*index+2]);
    }
}

KERNEL void addForcesSubset(GLOBAL const real* RESTRICT forces, GLOBAL mm_long* RESTRICT forceBuffers, GLOBAL int* RESTRICT atomIndex, GLOBAL int* RESTRICT particles, int numParticles) {
    for (int i = GLOBAL_ID; i < numParticles; i += GLOBAL_SIZE) {
        int index = particles[i];
        forceBuffers[index] += realToFixedPoint(forces[3*i]);
        forceBuffers[index+PADDED_NUM_ATOMS] += realToFixedPoint(forces[3*i+1]);
        forceBuffers[index+2*PADDED_NUM_ATOMS] += realToFixedPoint(forces[3*i+2]);
    }
}
