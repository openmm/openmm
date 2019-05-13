/**
 * Perform the first step of Velocity Verlet integration.
 */

extern "C" __global__ void integrateVelocityVerletPart1(int numAtoms, int paddedNumAtoms, const mixed2* __restrict__ dt, const real4* __restrict__ posq,
        const real4* __restrict__ posqCorrection, mixed4* __restrict__ velm, const long long* __restrict__ force, mixed4* __restrict__ posDelta) {
    const mixed2 stepSize = dt[0];
    const mixed dtPos = stepSize.y;
    const mixed dtVel = 0.5f*(stepSize.x+stepSize.y);
    const mixed scale = 0.5f*dtVel/(mixed) 0x100000000;
    for (int index = blockIdx.x*blockDim.x+threadIdx.x; index < numAtoms; index += blockDim.x*gridDim.x) {
        mixed4 velocity = velm[index];
        if (velocity.w != 0.0) {
#ifdef USE_MIXED_PRECISION
            real4 pos1 = posq[index];
            real4 pos2 = posqCorrection[index];
            mixed4 pos = make_mixed4(pos1.x+(mixed)pos2.x, pos1.y+(mixed)pos2.y, pos1.z+(mixed)pos2.z, pos1.w);
#else
            real4 pos = posq[index];
#endif
            velocity.x += scale*force[index]*velocity.w;
            velocity.y += scale*force[index+paddedNumAtoms]*velocity.w;
            velocity.z += scale*force[index+paddedNumAtoms*2]*velocity.w;
            pos.x = velocity.x*dtPos;
            pos.y = velocity.y*dtPos;
            pos.z = velocity.z*dtPos;
            posDelta[index] = pos;
            velm[index] = velocity;
        }
    }
}

/**
 * Perform the second step of Velocity Verlet integration.
 */

extern "C" __global__ void integrateVelocityVerletPart2(int numAtoms, mixed2* __restrict__ dt, real4* __restrict__ posq,
        real4* __restrict__ posqCorrection, mixed4* __restrict__ velm, const mixed4* __restrict__ posDelta) {
    mixed2 stepSize = dt[0];
    int index = blockIdx.x*blockDim.x+threadIdx.x;
    if (index == 0)
        dt[0].x = stepSize.y;
    for (; index < numAtoms; index += blockDim.x*gridDim.x) {
        mixed4 velocity = velm[index];
        if (velocity.w != 0.0) {
#ifdef USE_MIXED_PRECISION
            real4 pos1 = posq[index];
            real4 pos2 = posqCorrection[index];
            mixed4 pos = make_mixed4(pos1.x+(mixed)pos2.x, pos1.y+(mixed)pos2.y, pos1.z+(mixed)pos2.z, pos1.w);
#else
            real4 pos = posq[index];
#endif
            mixed4 delta = posDelta[index];
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

/**
 * Perform the third step of Velocity Verlet integration.
 */

extern "C" __global__ void integrateVelocityVerletPart3(int numAtoms, int paddedNumAtoms, mixed2* __restrict__ dt, real4* __restrict__ posq,
        real4* __restrict__ posqCorrection, mixed4* __restrict__ velm,  const long long* __restrict__ force, const mixed4* __restrict__ posDelta) {
    mixed2 stepSize = dt[0];
#if __CUDA_ARCH__ >= 130
    double oneOverDt = 1.0/stepSize.y;
#else
    float oneOverDt = 1.0f/stepSize.y;
    float correction = (1.0f-oneOverDt*stepSize.y)/stepSize.y;
#endif
    const mixed dtVel = 0.5f*(stepSize.x+stepSize.y);
    const mixed scale = 0.5f*dtVel/(mixed) 0x100000000;
    int index = blockIdx.x*blockDim.x+threadIdx.x;
    if (index == 0)
        dt[0].x = stepSize.y;
    for (; index < numAtoms; index += blockDim.x*gridDim.x) {
        mixed4 velocity = velm[index];
        if (velocity.w != 0.0) {
            mixed4 deltaXconstrained = posDelta[index];
            velocity.x += scale*force[index]*velocity.w + (deltaXconstrained.x - velocity.x*stepSize.y)*oneOverDt;
            velocity.y += scale*force[index+paddedNumAtoms]*velocity.w + (deltaXconstrained.y - velocity.y*stepSize.y)*oneOverDt;
            velocity.z += scale*force[index+paddedNumAtoms*2]*velocity.w + (deltaXconstrained.z - velocity.z*stepSize.y)*oneOverDt;
#if __CUDA_ARCH__ < 130
            velocity.x += (deltaXconstrained.x - velocity.x*stepSize.y)*correction;
            velocity.y += (deltaXconstrained.y - velocity.y*stepSize.y)*correction;
            velocity.z += (deltaXconstrained.z - velocity.z*stepSize.y)*correction;
#endif
            velm[index] = velocity;
        }
    }
}

