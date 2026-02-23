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
 * Sum kinetic energy using on-step velocities v(t) reconstructed from the
 * half-step velocities v(t+dt/2) = v_stored + F(t)*dt/(2m) after Verlet Part1.
 *
 * v(t) = v(t+dt/2) - F(t) * dt / (2m)
 *
 * This matches HOOMD's TwoStepConstantVolume, where the Bussi thermostat
 * computes alpha from KE(v(t)) (the on-step velocity), not from KE(v(t+dt/2)).
 */
KERNEL void sumBussiKineticEnergyPreKick(int numParticles, GLOBAL const mixed4* RESTRICT velm,
        GLOBAL const int* RESTRICT particleIndices, GLOBAL const float* RESTRICT masses,
        GLOBAL float* RESTRICT kineticEnergy, GLOBAL const mm_long* RESTRICT force,
        int paddedNumAtoms, mixed halfDt) {
    LOCAL float localKE[WORK_GROUP_SIZE];
    int localId = LOCAL_ID;
    const mixed forceScale = halfDt / (mixed) 0x100000000;

    float ke = 0.0f;
    for (int i = GLOBAL_ID; i < numParticles; i += GLOBAL_SIZE) {
        int idx = particleIndices[i];
        mixed4 v = velm[idx];
        float m = masses[i];
        if (m > 0) {
            mixed vx = v.x - forceScale * force[idx] * v.w;
            mixed vy = v.y - forceScale * force[idx + paddedNumAtoms] * v.w;
            mixed vz = v.z - forceScale * force[idx + paddedNumAtoms * 2] * v.w;
            ke += (float)(0.5f * m * (vx*vx + vy*vy + vz*vz));
        }
    }

    localKE[localId] = ke;
    SYNC_THREADS;

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
 * Scale position deltas by alpha for thermostatted particles only.
 */
KERNEL void scaleBussiPosDelta(int numParticles, GLOBAL mixed4* RESTRICT posDelta,
        GLOBAL const int* RESTRICT particleIndices, float alpha) {
    for (int i = GLOBAL_ID; i < numParticles; i += GLOBAL_SIZE) {
        int idx = particleIndices[i];
        mixed4 d = posDelta[idx];
        d.x *= alpha;
        d.y *= alpha;
        d.z *= alpha;
        posDelta[idx] = d;
    }
}
