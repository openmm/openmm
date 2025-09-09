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
 * Compute the correlation matrix, used in finding the optimal rotation.  This is executed as a single work group.
 */
KERNEL void computeCorrelationMatrix(int numParticles, GLOBAL const real4* RESTRICT posq, GLOBAL const real4* RESTRICT referencePos,
        GLOBAL const int* RESTRICT particles, GLOBAL real* buffer) {
    LOCAL volatile real temp[THREAD_BLOCK_SIZE];

    // Compute the center of the particle positions.

    real3 center = make_real3(0);
    for (int i = LOCAL_ID; i < numParticles; i += LOCAL_SIZE)
        center += trimTo3(posq[particles[i]]);
    center.x = reduceValue(center.x, temp)/numParticles;
    center.y = reduceValue(center.y, temp)/numParticles;
    center.z = reduceValue(center.z, temp)/numParticles;

    // Compute the correlation matrix.

    real R[3][3] = {{0, 0, 0}, {0, 0, 0}, {0, 0, 0}};
    for (int i = LOCAL_ID; i < numParticles; i += LOCAL_SIZE) {
        int index = particles[i];
        real3 pos = trimTo3(posq[index]) - center;
        real3 refPos = trimTo3(referencePos[index]);
        R[0][0] += refPos.x*pos.x;
        R[0][1] += refPos.x*pos.y;
        R[0][2] += refPos.x*pos.z;
        R[1][0] += refPos.y*pos.x;
        R[1][1] += refPos.y*pos.y;
        R[1][2] += refPos.y*pos.z;
        R[2][0] += refPos.z*pos.x;
        R[2][1] += refPos.z*pos.y;
        R[2][2] += refPos.z*pos.z;
    }
    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++)
            R[i][j] = reduceValue(R[i][j], temp);

    // Copy everything into the output buffer to send back to the host.

    if (LOCAL_ID == 0)
        for (int i = 0; i < 3; i++)
            for (int j = 0; j < 3; j++)
                buffer[3*i+j] = R[i][j];
}

/**
 * Apply forces based on the orientation.
 */
KERNEL void computeOrientationForces(int numParticles, int paddedNumAtoms, GLOBAL const real4* RESTRICT referencePos,
        GLOBAL const int* RESTRICT particles, const real dxdq, const real4 eigenvalues, GLOBAL const real4* RESTRICT eigenvectors,
        GLOBAL mm_long* RESTRICT forceBuffers) {
    const int3 dsIndex[4][4] = {
        {make_int3(0, 1, 2), make_int3(0, 2, 1), make_int3(2, 0, 0), make_int3(1, 0, 0)},
        {make_int3(0, 2, 1), make_int3(0, 1, 2), make_int3(1, 0, 0), make_int3(2, 0, 0)},
        {make_int3(2, 0, 0), make_int3(1, 0, 0), make_int3(0, 1, 2), make_int3(0, 2, 1)},
        {make_int3(1, 0, 0), make_int3(2, 0, 0), make_int3(0, 2, 1), make_int3(0, 1, 2)}
    };
    const int3 dsScale[4][4] = {
        {make_int3(1, 1, 1), make_int3(0, -1, 1), make_int3(1, 0, -1), make_int3(-1, 1, 0)},
        {make_int3(0, -1, 1), make_int3(1, -1, -1), make_int3(1, 1, 0), make_int3(1, 0, 1)},
        {make_int3(1, 0, -1), make_int3(1, 1, 0), make_int3(-1, 1, -1), make_int3(0, 1, 1)},
        {make_int3(-1, 1, 0), make_int3(1, 0, 1), make_int3(0, 1, 1), make_int3(-1, -1, 1)}
    };
    for (int index = GLOBAL_ID; index < numParticles; index += GLOBAL_SIZE) {
        int particle = particles[index];
        real3 refPos = trimTo3(referencePos[particle]);
        real p[3] = {refPos.x, refPos.y, refPos.z};
        real3 dq = make_real3(0);
        for (int i = 0; i < 4; i++)
            for (int j = 0; j < 4; j++) {
                real scale = ((eigenvectors[i].z*eigenvectors[0].z) / (eigenvalues.w-eigenvalues.z) +
                              (eigenvectors[i].y*eigenvectors[0].y) / (eigenvalues.w-eigenvalues.y) +
                              (eigenvectors[i].x*eigenvectors[0].x) / (eigenvalues.w-eigenvalues.x)) * eigenvectors[j].w;
                real3 ds = make_real3(dsScale[i][j].x*p[dsIndex[i][j].x],
                                      dsScale[i][j].y*p[dsIndex[i][j].y],
                                      dsScale[i][j].z*p[dsIndex[i][j].z]);
                dq += scale*ds;
            }
        real3 force = -dxdq*dq;
        forceBuffers[particle] += realToFixedPoint(force.x);
        forceBuffers[particle+paddedNumAtoms] += realToFixedPoint(force.y);
        forceBuffers[particle+2*paddedNumAtoms] += realToFixedPoint(force.z);
    }
}
