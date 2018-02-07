// This file contains kernels to compute the RMSD and its gradient using the algorithm described
// in Coutsias et al, "Using quaternions to calculate RMSD" (doi: 10.1002/jcc.20110).

/**
 * Sum a value over all threads.
 */
real reduceValue(real value, __local volatile real* temp) {
    const int thread = get_local_id(0);
    barrier(CLK_LOCAL_MEM_FENCE);
    temp[thread] = value;
    barrier(CLK_LOCAL_MEM_FENCE);
    for (uint step = 1; step < 32; step *= 2) {
        if (thread+step < get_local_size(0) && thread%(2*step) == 0)
            temp[thread] = temp[thread] + temp[thread+step];
        SYNC_WARPS;
    }
    for (uint step = 32; step < get_local_size(0); step *= 2) {
        if (thread+step < get_local_size(0) && thread%(2*step) == 0)
            temp[thread] = temp[thread] + temp[thread+step];
        barrier(CLK_LOCAL_MEM_FENCE);
    }
    return temp[0];
}

/**
 * Perform the first step of computing the RMSD.  This is executed as a single work group.
 */
__kernel void computeRMSDPart1(int numParticles, __global const real4* restrict posq, __global const real4* restrict referencePos,
        __global const int* restrict particles, __global real* buffer, __local volatile real* restrict temp) {
    // Compute the center of the particle positions.
    
    real3 center = (real3) 0;
    for (int i = get_local_id(0); i < numParticles; i += get_local_size(0))
        center += posq[particles[i]].xyz;
    center.x = reduceValue(center.x, temp)/numParticles;
    center.y = reduceValue(center.y, temp)/numParticles;
    center.z = reduceValue(center.z, temp)/numParticles;
    
    // Compute the correlation matrix.
    
    real R[3][3] = {{0, 0, 0}, {0, 0, 0}, {0, 0, 0}};
    real sum = 0;
    for (int i = get_local_id(0); i < numParticles; i += get_local_size(0)) {
        int index = particles[i];
        real3 pos = posq[index].xyz - center;
        real3 refPos = referencePos[index].xyz;
        R[0][0] += pos.x*refPos.x;
        R[0][1] += pos.x*refPos.y;
        R[0][2] += pos.x*refPos.z;
        R[1][0] += pos.y*refPos.x;
        R[1][1] += pos.y*refPos.y;
        R[1][2] += pos.y*refPos.z;
        R[2][0] += pos.z*refPos.x;
        R[2][1] += pos.z*refPos.y;
        R[2][2] += pos.z*refPos.z;
        sum += dot(pos, pos);
    }
    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++)
            R[i][j] = reduceValue(R[i][j], temp);
    sum = reduceValue(sum, temp);

    // Copy everything into the output buffer to send back to the host.
    
    if (get_local_id(0) == 0) {
        for (int i = 0; i < 3; i++)
            for (int j = 0; j < 3; j++)
                buffer[3*i+j] = R[i][j];
        buffer[9] = sum;
        buffer[10] = center.x;
        buffer[11] = center.y;
        buffer[12] = center.z;
    }
}

/**
 * Apply forces based on the RMSD.
 */
__kernel void computeRMSDForces(int numParticles, __global const real4* restrict posq, __global const real4* restrict referencePos,
        __global const int* restrict particles, __global const real* buffer, __global real4* restrict forceBuffers) {
    real3 center = (real3) (buffer[10], buffer[11], buffer[12]);
    real scale = 1 / (real) (buffer[9]*numParticles);
    for (int i = get_global_id(0); i < numParticles; i += get_global_size(0)) {
        int index = particles[i];
        real3 pos = posq[index].xyz - center;
        real3 refPos = referencePos[index].xyz;
        real3 rotatedRef = (real3) (buffer[0]*refPos.x + buffer[3]*refPos.y + buffer[6]*refPos.z,
                                    buffer[1]*refPos.x + buffer[4]*refPos.y + buffer[7]*refPos.z,
                                    buffer[2]*refPos.x + buffer[5]*refPos.y + buffer[8]*refPos.z);
        forceBuffers[index].xyz -= (pos-rotatedRef)*scale;
    }
}
