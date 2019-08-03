// This file contains kernels to compute the RMSD and its gradient using the algorithm described
// in Coutsias et al, "Using quaternions to calculate RMSD" (doi: 10.1002/jcc.20110).

/**
 * Sum a value over all threads.
 */
__device__ real reduceValue(real value, volatile real* temp) {
    const int thread = threadIdx.x;
    __syncthreads();
    temp[thread] = value;
    __syncthreads();
    for (unsigned int step = 1; step < 32; step *= 2) {
        if (thread+step < blockDim.x && thread%(2*step) == 0)
            temp[thread] = temp[thread] + temp[thread+step];
        SYNC_WARPS
    }
    for (unsigned int step = 32; step < blockDim.x; step *= 2) {
        if (thread+step < blockDim.x && thread%(2*step) == 0)
            temp[thread] = temp[thread] + temp[thread+step];
        __syncthreads();
    }
    return temp[0];
}

/**
 * Perform the first step of computing the RMSD.  This is executed as a single work group.
 */
extern "C" __global__ void computeRMSDPart1(int numParticles, const real4* __restrict__ posq, const real4* __restrict__ referencePos,
         const int* __restrict__ particles, real* buffer) {
    extern __shared__ volatile real temp[];

    // Compute the center of the particle positions.
    
    real3 center = make_real3(0);
    for (int i = threadIdx.x; i < numParticles; i += blockDim.x)
        center += trimTo3(posq[particles[i]]);
    center.x = reduceValue(center.x, temp)/numParticles;
    center.y = reduceValue(center.y, temp)/numParticles;
    center.z = reduceValue(center.z, temp)/numParticles;
    
    // Compute the correlation matrix.
    
    real R[3][3] = {{0, 0, 0}, {0, 0, 0}, {0, 0, 0}};
    real sum = 0;
    for (int i = threadIdx.x; i < numParticles; i += blockDim.x) {
        int index = particles[i];
        real3 pos = trimTo3(posq[index]) - center;
        real3 refPos = trimTo3(referencePos[index]);
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
    
    if (threadIdx.x == 0) {
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
extern "C" __global__ void computeRMSDForces(int numParticles, int paddedNumAtoms, const real4* __restrict__ posq, const real4* __restrict__ referencePos,
         const int* __restrict__ particles, const real* buffer, unsigned long long* __restrict__ forceBuffers) {
    real3 center = make_real3(buffer[10], buffer[11], buffer[12]);
    real scale = 1 / (real) (buffer[9]*numParticles);
    for (int i = blockDim.x*blockIdx.x+threadIdx.x; i < numParticles; i += blockDim.x*gridDim.x) {
        int index = particles[i];
        real3 pos = trimTo3(posq[index]) - center;
        real3 refPos = trimTo3(referencePos[index]);
        real3 rotatedRef = make_real3(buffer[0]*refPos.x + buffer[3]*refPos.y + buffer[6]*refPos.z,
                                      buffer[1]*refPos.x + buffer[4]*refPos.y + buffer[7]*refPos.z,
                                      buffer[2]*refPos.x + buffer[5]*refPos.y + buffer[8]*refPos.z);
        real3 force = (rotatedRef-pos)*scale;
        atomicAdd(&forceBuffers[index], static_cast<unsigned long long>((long long) (force.x*0x100000000)));
        atomicAdd(&forceBuffers[index+paddedNumAtoms], static_cast<unsigned long long>((long long) (force.y*0x100000000)));
        atomicAdd(&forceBuffers[index+2*paddedNumAtoms], static_cast<unsigned long long>((long long) (force.z*0x100000000)));
    }
}
