enum {VelScale, ForceScale, NoiseScale, MaxParams};

/**
 * Perform the first step of Langevin integration.
 */

extern "C" __global__ void integrateLangevinPart1(int numAtoms, int paddedNumAtoms, mixed4* __restrict__ velm, const long long* __restrict__ force, mixed4* __restrict__ posDelta,
        const mixed* __restrict__ paramBuffer, const mixed2* __restrict__ dt, const float4* __restrict__ random, unsigned int randomIndex) {
    mixed vscale = paramBuffer[VelScale];
    mixed fscale = paramBuffer[ForceScale]/(mixed) 0x100000000;
    mixed noisescale = paramBuffer[NoiseScale];
    mixed stepSize = dt[0].y;
    int index = blockIdx.x*blockDim.x+threadIdx.x;
    randomIndex += index;
    while (index < numAtoms) {
        mixed4 velocity = velm[index];
        if (velocity.w != 0) {
            mixed sqrtInvMass = SQRT(velocity.w);
            velocity.x = vscale*velocity.x + fscale*velocity.w*force[index] + noisescale*sqrtInvMass*random[randomIndex].x;
            velocity.y = vscale*velocity.y + fscale*velocity.w*force[index+paddedNumAtoms] + noisescale*sqrtInvMass*random[randomIndex].y;
            velocity.z = vscale*velocity.z + fscale*velocity.w*force[index+paddedNumAtoms*2] + noisescale*sqrtInvMass*random[randomIndex].z;
            velm[index] = velocity;
            posDelta[index] = make_mixed4(stepSize*velocity.x, stepSize*velocity.y, stepSize*velocity.z, 0);
        }
        randomIndex += blockDim.x*gridDim.x;
        index += blockDim.x*gridDim.x;
    }
}

/**
 * Perform the second step of Langevin integration.
 */

extern "C" __global__ void integrateLangevinPart2(int numAtoms, real4* __restrict__ posq, real4* __restrict__ posqCorrection, const mixed4* __restrict__ posDelta, mixed4* __restrict__ velm, const mixed2* __restrict__ dt) {
#if __CUDA_ARCH__ >= 130
    double invStepSize = 1.0/dt[0].y;
#else
    float invStepSize = 1.0f/dt[0].y;
    float correction = (1.0f-invStepSize*dt[0].y)/dt[0].y;
#endif
    int index = blockIdx.x*blockDim.x+threadIdx.x;
    while (index < numAtoms) {
        mixed4 vel = velm[index];
        if (vel.w != 0) {
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
#if __CUDA_ARCH__ >= 130
            vel.x = (mixed) (invStepSize*delta.x);
            vel.y = (mixed) (invStepSize*delta.y);
            vel.z = (mixed) (invStepSize*delta.z);
#else
            vel.x = invStepSize*delta.x + correction*delta.x;
            vel.y = invStepSize*delta.y + correction*delta.x;
            vel.z = invStepSize*delta.z + correction*delta.x;
#endif
#ifdef USE_MIXED_PRECISION
            posq[index] = make_real4((real) pos.x, (real) pos.y, (real) pos.z, (real) pos.w);
            posqCorrection[index] = make_real4(pos.x-(real) pos.x, pos.y-(real) pos.y, pos.z-(real) pos.z, 0);
#else
            posq[index] = pos;
#endif
            velm[index] = vel;
        }
        index += blockDim.x*gridDim.x;
    }
}

/**
 * Select the step size to use for the next step.
 */

extern "C" __global__ void selectLangevinStepSize(int numAtoms, int paddedNumAtoms, mixed maxStepSize, mixed errorTol, mixed friction, mixed kT, mixed2* __restrict__ dt,
        const mixed4* __restrict__ velm, const long long* __restrict__ force, mixed* __restrict__ paramBuffer) {
    // Calculate the error.

    extern __shared__ mixed params[];
    mixed* error = &params[MaxParams];
    mixed err = 0;
    unsigned int index = threadIdx.x;
    const mixed scale = RECIP((mixed) 0x100000000);
    while (index < numAtoms) {
        mixed3 f = make_mixed3(scale*force[index], scale*force[index+paddedNumAtoms], scale*force[index+paddedNumAtoms*2]);
        mixed invMass = velm[index].w;
        err += (f.x*f.x + f.y*f.y + f.z*f.z)*invMass*invMass;
        index += blockDim.x*gridDim.x;
    }
    error[threadIdx.x] = err;
    __syncthreads();

    // Sum the errors from all threads.

    for (unsigned int offset = 1; offset < blockDim.x; offset *= 2) {
        if (threadIdx.x+offset < blockDim.x && (threadIdx.x&(2*offset-1)) == 0)
            error[threadIdx.x] += error[threadIdx.x+offset];
        __syncthreads();
    }
    if (blockIdx.x*blockDim.x+threadIdx.x == 0) {
        // Select the new step size.

        mixed totalError = SQRT(error[0]/(numAtoms*3));
        mixed newStepSize = SQRT(errorTol/totalError);
        mixed oldStepSize = dt[0].y;
        if (oldStepSize > 0.0f)
            newStepSize = min(newStepSize, oldStepSize*2.0f); // For safety, limit how quickly dt can increase.
        if (newStepSize > oldStepSize && newStepSize < 1.1f*oldStepSize)
            newStepSize = oldStepSize; // Keeping dt constant between steps improves the behavior of the integrator.
        if (newStepSize > maxStepSize)
            newStepSize = maxStepSize;
        dt[0].y = newStepSize;

        // Recalculate the integration parameters.

        mixed vscale = exp(-newStepSize*friction);
        mixed fscale = (friction == 0 ? newStepSize : (1-vscale)/friction);
        mixed noisescale = sqrt(kT*(1-vscale*vscale));
        params[VelScale] = vscale;
        params[ForceScale] = fscale;
        params[NoiseScale] = noisescale;
    }
    __syncthreads();
    if (threadIdx.x < MaxParams)
        paramBuffer[threadIdx.x] = params[threadIdx.x];
}
