/**
 * Perform the first step of Brownian integration.
 */

extern "C" __global__ void integrateBrownianPart1(int numAtoms, int paddedNumAtoms, mixed tauDeltaT, mixed noiseAmplitude, const long long* __restrict__ force,
        mixed4* __restrict__ posDelta, const mixed4* __restrict__ velm, const float4* __restrict__ random, unsigned int randomIndex) {
    randomIndex += blockIdx.x*blockDim.x+threadIdx.x;
    const mixed fscale = tauDeltaT/(mixed) 0x100000000;
    for (int index = blockIdx.x*blockDim.x+threadIdx.x; index < numAtoms; index += blockDim.x*gridDim.x) {
        mixed invMass = velm[index].w;
        if (invMass != 0) {
            posDelta[index].x = fscale*invMass*force[index] + noiseAmplitude*SQRT(invMass)*random[randomIndex].x;
            posDelta[index].y = fscale*invMass*force[index+paddedNumAtoms] + noiseAmplitude*SQRT(invMass)*random[randomIndex].y;
            posDelta[index].z = fscale*invMass*force[index+paddedNumAtoms*2] + noiseAmplitude*SQRT(invMass)*random[randomIndex].z;
        }
        randomIndex += blockDim.x*gridDim.x;
    }
}

/**
 * Perform the second step of Brownian integration.
 */

extern "C" __global__ void integrateBrownianPart2(int numAtoms, mixed deltaT, real4* posq, real4* __restrict__ posqCorrection, mixed4* velm, const mixed4* __restrict__ posDelta) {
    const mixed oneOverDeltaT = RECIP(deltaT);
    for (int index = blockIdx.x*blockDim.x+threadIdx.x; index < numAtoms; index += blockDim.x*gridDim.x) {
        if (velm[index].w != 0) {
            mixed4 delta = posDelta[index];
            velm[index].x = oneOverDeltaT*delta.x;
            velm[index].y = oneOverDeltaT*delta.y;
            velm[index].z = oneOverDeltaT*delta.z;
#ifdef USE_MIXED_PRECISION
            real4 pos1 = posq[index];
            real4 pos2 = posqCorrection[index];
            mixed4 pos = make_mixed4(pos1.x+(mixed)pos2.x, pos1.y+(mixed)pos2.y, pos1.z+(mixed)pos2.z, pos1.w);
#else
            real4 pos = posq[index];
#endif
            pos.x += delta.x;
            pos.y += delta.y;
            pos.z += delta.z;
#ifdef USE_MIXED_PRECISION
            posq[index] = make_real4((real) pos.x, (real) pos.y, (real) pos.z, (real) pos.w);
            posqCorrection[index] = make_real4(pos.x-(real) pos.x, pos.y-(real) pos.y, pos.z-(real) pos.z, 0);
#else
            posq[index] = pos;
#endif
        }
    }
}
