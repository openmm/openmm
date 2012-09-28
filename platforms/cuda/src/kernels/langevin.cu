enum {VelScale, ForceScale, NoiseScale, MaxParams};

/**
 * Perform the first step of Langevin integration.
 */

extern "C" __global__ void integrateLangevinPart1(real4* __restrict__ velm, const long long* __restrict__ force, real4* __restrict__ posDelta,
        const real* __restrict__ paramBuffer, const real2* __restrict__ dt, const float4* __restrict__ random, unsigned int randomIndex) {
    real vscale = paramBuffer[VelScale];
    real fscale = paramBuffer[ForceScale]/(real) 0xFFFFFFFF;
    real noisescale = paramBuffer[NoiseScale];
    real stepSize = dt[0].y;
    int index = blockIdx.x*blockDim.x+threadIdx.x;
    randomIndex += index;
    while (index < NUM_ATOMS) {
        real4 velocity = velm[index];
        if (velocity.w != 0) {
            real sqrtInvMass = SQRT(velocity.w);
            velocity.x = vscale*velocity.x + fscale*velocity.w*force[index] + noisescale*sqrtInvMass*random[randomIndex].x;
            velocity.y = vscale*velocity.y + fscale*velocity.w*force[index+PADDED_NUM_ATOMS] + noisescale*sqrtInvMass*random[randomIndex].y;
            velocity.z = vscale*velocity.z + fscale*velocity.w*force[index+PADDED_NUM_ATOMS*2] + noisescale*sqrtInvMass*random[randomIndex].z;
            velm[index] = velocity;
            posDelta[index] = make_real4(stepSize*velocity.x, stepSize*velocity.y, stepSize*velocity.z, 0);
        }
        randomIndex += blockDim.x*gridDim.x;
        index += blockDim.x*gridDim.x;
    }
}

/**
 * Perform the second step of Langevin integration.
 */

extern "C" __global__ void integrateLangevinPart2(real4* __restrict__ posq, const real4* __restrict__ posDelta, real4* __restrict__ velm, const real2* __restrict__ dt) {
    double invStepSize = 1.0/dt[0].y;
    int index = blockIdx.x*blockDim.x+threadIdx.x;
    while (index < NUM_ATOMS) {
        real4 vel = velm[index];
        if (vel.w != 0) {
            real4 pos = posq[index];
            real4 delta = posDelta[index];
            pos.x += delta.x;
            pos.y += delta.y;
            pos.z += delta.z;
            vel.x = (real) invStepSize*delta.x;
            vel.y = (real) invStepSize*delta.y;
            vel.z = (real) invStepSize*delta.z;
            posq[index] = pos;
            velm[index] = vel;
        }
        index += blockDim.x*gridDim.x;
    }
}

/**
 * Select the step size to use for the next step.
 */

extern "C" __global__ void selectLangevinStepSize(real maxStepSize, real errorTol, real tau, real kT, real2* __restrict__ dt,
        const real4* __restrict__ velm, const long long* __restrict__ force, real* __restrict__ paramBuffer) {
    // Calculate the error.

    extern __shared__ real params[];
    real* error = &params[MaxParams];
    real err = 0;
    unsigned int index = threadIdx.x;
    const real scale = RECIP((real) 0xFFFFFFFF);
    while (index < NUM_ATOMS) {
        real3 f = make_real3(scale*force[index], scale*force[index+PADDED_NUM_ATOMS], scale*force[index+PADDED_NUM_ATOMS*2]);
        real invMass = velm[index].w;
        err += (f.x*f.x + f.y*f.y + f.z*f.z)*invMass;
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

        real totalError = sqrt(error[0]/(NUM_ATOMS*3));
        real newStepSize = sqrt(errorTol/totalError);
        real oldStepSize = dt[0].y;
        if (oldStepSize > 0.0f)
            newStepSize = min(newStepSize, oldStepSize*2.0f); // For safety, limit how quickly dt can increase.
        if (newStepSize > oldStepSize && newStepSize < 1.1f*oldStepSize)
            newStepSize = oldStepSize; // Keeping dt constant between steps improves the behavior of the integrator.
        if (newStepSize > maxStepSize)
            newStepSize = maxStepSize;
        dt[0].y = newStepSize;

        // Recalculate the integration parameters.

        real vscale = exp(-newStepSize/tau);
        real fscale = (1-vscale)*tau;
        real noisescale = sqrt(2*kT/tau)*sqrt(0.5f*(1-vscale*vscale)*tau);
        params[VelScale] = vscale;
        params[ForceScale] = fscale;
        params[NoiseScale] = noisescale;
    }
    __syncthreads();
    if (threadIdx.x < MaxParams)
        paramBuffer[threadIdx.x] = params[threadIdx.x];
}
