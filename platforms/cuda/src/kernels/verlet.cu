/**
 * Perform the first step of Verlet integration.
 */

extern "C" __global__ void integrateVerletPart1(const real2* __restrict__ dt, const real4* __restrict__ posq, real4* __restrict__ velm, const long long* __restrict__ force, real4* __restrict__ posDelta) {
    const real2 stepSize = dt[0];
    const real dtPos = stepSize.y;
    const real dtVel = 0.5f*(stepSize.x+stepSize.y);
    const real scale = dtVel/(real) 0xFFFFFFFF;
    for (int index = blockIdx.x*blockDim.x+threadIdx.x; index < NUM_ATOMS; index += blockDim.x*gridDim.x) {
        real4 velocity = velm[index];
        if (velocity.w != 0.0) {
            real4 pos = posq[index];
            velocity.x += scale*force[index]*velocity.w;
            velocity.y += scale*force[index+PADDED_NUM_ATOMS]*velocity.w;
            velocity.z += scale*force[index+PADDED_NUM_ATOMS*2]*velocity.w;
            pos.x = velocity.x*dtPos;
            pos.y = velocity.y*dtPos;
            pos.z = velocity.z*dtPos;
            posDelta[index] = pos;
            velm[index] = velocity;
        }
    }
}

/**
 * Perform the second step of Verlet integration.
 */

extern "C" __global__ void integrateVerletPart2(real2* __restrict__ dt, real4* __restrict__ posq, real4* __restrict__ velm, const real4* __restrict__ posDelta) {
    real2 stepSize = dt[0];
    double oneOverDt = 1.0/stepSize.y;
    int index = blockIdx.x*blockDim.x+threadIdx.x;
    if (index == 0)
        dt[0].x = stepSize.y;
    for (; index < NUM_ATOMS; index += blockDim.x*gridDim.x) {
        real4 velocity = velm[index];
        if (velocity.w != 0.0) {
            real4 pos = posq[index];
            real4 delta = posDelta[index];
            pos.x += delta.x;
            pos.y += delta.y;
            pos.z += delta.z;
            velocity = make_real4((real) (delta.x*oneOverDt), (real) (delta.y*oneOverDt), (real) (delta.z*oneOverDt), velocity.w);
            posq[index] = pos;
            velm[index] = velocity;
        }
    }
}

/**
 * Select the step size to use for the next step.
 */

extern "C" __global__ void selectVerletStepSize(real maxStepSize, real errorTol, real2* __restrict__ dt, const real4* __restrict__ velm, const long long* __restrict__ force) {
    // Calculate the error.

    extern __shared__ real error[];
    real err = 0.0f;
    const real scale = RECIP((real) 0xFFFFFFFF);
    for (int index = threadIdx.x; index < NUM_ATOMS; index += blockDim.x*gridDim.x) {
        real3 f = make_real3(scale*force[index], scale*force[index+PADDED_NUM_ATOMS], scale*force[index+PADDED_NUM_ATOMS*2]);
        real invMass = velm[index].w;
        err += (f.x*f.x + f.y*f.y + f.z*f.z)*invMass;
    }
    error[threadIdx.x] = err;
    __syncthreads();

    // Sum the errors from all threads.

    for (unsigned int offset = 1; offset < blockDim.x; offset *= 2) {
        if (threadIdx.x+offset < blockDim.x && (threadIdx.x&(2*offset-1)) == 0)
            error[threadIdx.x] += error[threadIdx.x+offset];
        __syncthreads();
    }
    if (threadIdx.x == 0) {
        real totalError = SQRT(error[0]/(NUM_ATOMS*3));
        real newStepSize = SQRT(errorTol/totalError);
        real oldStepSize = dt[0].y;
        if (oldStepSize > 0.0f)
            newStepSize = min(newStepSize, oldStepSize*2.0f); // For safety, limit how quickly dt can increase.
        if (newStepSize > oldStepSize && newStepSize < 1.1f*oldStepSize)
            newStepSize = oldStepSize; // Keeping dt constant between steps improves the behavior of the integrator.
        if (newStepSize > maxStepSize)
            newStepSize = maxStepSize;
        dt[0].y = newStepSize;
    }
}
