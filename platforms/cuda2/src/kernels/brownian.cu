/**
 * Perform the first step of Brownian integration.
 */

extern "C" __global__ void integrateBrownianPart1(real tauDeltaT, real noiseAmplitude, const long long* __restrict__ force,
        real4* __restrict__ posDelta, const real4* __restrict__ velm, const float4* __restrict__ random, unsigned int randomIndex) {
    randomIndex += blockIdx.x*blockDim.x+threadIdx.x;
    const real fscale = tauDeltaT/(real) 0xFFFFFFFF;
    for (int index = blockIdx.x*blockDim.x+threadIdx.x; index < NUM_ATOMS; index += blockDim.x*gridDim.x) {
        real invMass = velm[index].w;
        if (invMass != 0) {
            posDelta[index].x = fscale*invMass*force[index] + noiseAmplitude*SQRT(invMass)*random[randomIndex].x;
            posDelta[index].y = fscale*invMass*force[index+PADDED_NUM_ATOMS] + noiseAmplitude*SQRT(invMass)*random[randomIndex].y;
            posDelta[index].z = fscale*invMass*force[index+PADDED_NUM_ATOMS*2] + noiseAmplitude*SQRT(invMass)*random[randomIndex].z;
        }
        randomIndex += blockDim.x*gridDim.x;
    }
}

/**
 * Perform the second step of Brownian integration.
 */

extern "C" __global__ void integrateBrownianPart2(real deltaT, real4* posq, real4* velm, const real4* __restrict__ posDelta) {
    const real oneOverDeltaT = RECIP(deltaT);
    for (int index = blockIdx.x*blockDim.x+threadIdx.x; index < NUM_ATOMS; index += blockDim.x*gridDim.x) {
        if (velm[index].w != 0) {
            real4 delta = posDelta[index];
            velm[index].x = oneOverDeltaT*delta.x;
            velm[index].y = oneOverDeltaT*delta.y;
            velm[index].z = oneOverDeltaT*delta.z;
            posq[index].x = posq[index].x + delta.x;
            posq[index].y = posq[index].y + delta.y;
            posq[index].z = posq[index].z + delta.z;
        }
    }
}
