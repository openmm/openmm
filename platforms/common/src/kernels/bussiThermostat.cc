/**
 * Sum kinetic energy of particles in the Bussi thermostat subset.
 */
KERNEL void sumBussiKineticEnergy(int numParticles, GLOBAL const mixed4* RESTRICT velm,
        GLOBAL const int* RESTRICT particleIndices, GLOBAL const float* RESTRICT masses,
        GLOBAL float* RESTRICT kineticEnergy) {
    LOCAL float localKE[WORK_GROUP_SIZE];
    int localId = LOCAL_ID;
    
    float ke = 0.0f;
    for (int i = GLOBAL_ID; i < numParticles; i += GLOBAL_SIZE) {
        int idx = particleIndices[i];
        mixed4 v = velm[idx];
        float m = masses[i];
        if (m > 0) {
            ke += 0.5f * m * (v.x*v.x + v.y*v.y + v.z*v.z);
        }
    }
    
    localKE[localId] = ke;
    SYNC_THREADS;
    
    // Tree reduction
    for (int stride = WORK_GROUP_SIZE/2; stride > 0; stride >>= 1) {
        if (localId < stride)
            localKE[localId] += localKE[localId + stride];
        SYNC_THREADS;
    }
    
    if (localId == 0)
        kineticEnergy[GROUP_ID] = localKE[0];
}

/**
 * Apply the Bussi thermostat by rescaling velocities.
 */
KERNEL void applyBussiThermostat(int numParticles, GLOBAL mixed4* RESTRICT velm,
        GLOBAL const int* RESTRICT particleIndices, float alpha) {
    for (int i = GLOBAL_ID; i < numParticles; i += GLOBAL_SIZE) {
        int idx = particleIndices[i];
        mixed4 v = velm[idx];
        v.x *= alpha;
        v.y *= alpha;
        v.z *= alpha;
        velm[idx] = v;
    }
}

/**
 * Scale position deltas by alpha (when Bussi is applied after Verlet part1).
 */
KERNEL void scaleBussiPosDelta(int numAtoms, GLOBAL mixed4* RESTRICT posDelta, float alpha) {
    for (int i = GLOBAL_ID; i < numAtoms; i += GLOBAL_SIZE) {
        mixed4 d = posDelta[i];
        d.x *= alpha;
        d.y *= alpha;
        d.z *= alpha;
        posDelta[i] = d;
    }
}
