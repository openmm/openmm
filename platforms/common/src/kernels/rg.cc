/**
 * Sum a value over all threads.
 */
DEVICE real reduceValue(real value, LOCAL_ARG volatile real* temp) {
    const int thread = LOCAL_ID;
    SYNC_THREADS;
    temp[thread] = value;
    SYNC_THREADS;
    for (int step = 1; step < 32; step *= 2) {
        if (thread+step < LOCAL_SIZE && thread%(2*step) == 0)
            temp[thread] = temp[thread] + temp[thread+step];
        SYNC_WARPS;
    }
    for (int step = 32; step < LOCAL_SIZE; step *= 2) {
        if (thread+step < LOCAL_SIZE && thread%(2*step) == 0)
            temp[thread] = temp[thread] + temp[thread+step];
        SYNC_THREADS;
    }
    return temp[0];
}

/**
 * First step of computing the radius of gyration: each thread block computes a contribution
 * to the center position.
 */
KERNEL void computeCenterPosition(int numParticles, GLOBAL const real4* RESTRICT posq, GLOBAL const int* RESTRICT particles,
        GLOBAL real* centerBuffer) {
    LOCAL volatile real temp[THREAD_BLOCK_SIZE];
    real3 center = make_real3(0);
    for (int i = GLOBAL_ID; i < numParticles; i += GLOBAL_SIZE)
        center += trimTo3(posq[particles[i]]);
    center.x = reduceValue(center.x, temp);
    center.y = reduceValue(center.y, temp);
    center.z = reduceValue(center.z, temp);
    if (LOCAL_ID == 0) {
        centerBuffer[3*GROUP_ID] = center.x;
        centerBuffer[3*GROUP_ID+1] = center.y;
        centerBuffer[3*GROUP_ID+2] = center.z;
    }
}

/**
 * Second step of computing the radius of gyration: each thread block computes a contribution to Rg^2.
 */
KERNEL void computeRg(int numParticles, GLOBAL const real4* RESTRICT posq, GLOBAL const int* RESTRICT particles,
        GLOBAL real* centerBuffer, GLOBAL real* rgBuffer) {
    LOCAL volatile real temp[THREAD_BLOCK_SIZE];

    // Sum the contributions computed in the previous kernel to find the center position.

    real3 center = make_real3(0);
    for (int i = LOCAL_ID; i < NUM_GROUPS; i += LOCAL_SIZE) {
        center.x += centerBuffer[3*i];
        center.y += centerBuffer[3*i+1];
        center.z += centerBuffer[3*i+2];
    }
    center.x = reduceValue(center.x, temp)/numParticles;
    center.y = reduceValue(center.y, temp)/numParticles;
    center.z = reduceValue(center.z, temp)/numParticles;
    if (GLOBAL_ID == 0) {
        // Save the result so we don't need to compute it again in the next kernel.

        centerBuffer[3*NUM_GROUPS] = center.x;
        centerBuffer[3*NUM_GROUPS+1] = center.y;
        centerBuffer[3*NUM_GROUPS+2] = center.z;
    }

    // Compute the contributions to Rg^2.

    real rg2 = 0;
    for (int i = GLOBAL_ID; i < numParticles; i += GLOBAL_SIZE) {
        real3 delta = trimTo3(posq[particles[i]]) - center;
        rg2 += dot(delta, delta);
    }
    rg2 = reduceValue(rg2, temp);
    if (LOCAL_ID == 0)
        rgBuffer[GROUP_ID] = rg2;
}

/**
 * Third step of computing the radius of gyration: compute the forces.
 */
KERNEL void computeForces(int numParticles, GLOBAL const real4* RESTRICT posq, GLOBAL const int* RESTRICT particles,
        GLOBAL const real* centerBuffer, GLOBAL real* rgBuffer, GLOBAL mm_ulong* RESTRICT forceBuffers, GLOBAL mixed* RESTRICT energyBuffer) {
    LOCAL volatile real temp[THREAD_BLOCK_SIZE];
    real3 center = make_real3(centerBuffer[3*NUM_GROUPS], centerBuffer[3*NUM_GROUPS+1], centerBuffer[3*NUM_GROUPS+2]);

    // Sum the contributions computed in the previous kernel to find Rg.

    real rg2 = 0;
    for (int i = LOCAL_ID; i < NUM_GROUPS; i += LOCAL_SIZE)
        rg2 += rgBuffer[i];
    real rg = SQRT(reduceValue(rg2, temp)/numParticles);
    if (GLOBAL_ID == 0)
        energyBuffer[0] += rg;

    // Compute the forces.

    real scale = 1/(rg*numParticles);
    for (int i = GLOBAL_ID; i < numParticles; i += GLOBAL_SIZE) {
        int particle = particles[i];
        real3 delta = trimTo3(posq[particle]) - center;
        ATOMIC_ADD(&forceBuffers[particle], (mm_ulong) realToFixedPoint(-scale*delta.x));
        ATOMIC_ADD(&forceBuffers[particle+PADDED_NUM_ATOMS], (mm_ulong) realToFixedPoint(-scale*delta.y));
        ATOMIC_ADD(&forceBuffers[particle+2*PADDED_NUM_ATOMS], (mm_ulong) realToFixedPoint(-scale*delta.z));
    }
}
