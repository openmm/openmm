__device__ float3 multiplyComplexRealPart(float2 c1, float3 c2r, float3 c2i) {
    return c1.x*c2r-c1.y*c2i;
}

__device__ float3 multiplyComplexImagPart(float2 c1, float3 c2r, float3 c2i) {
    return c1.x*c2i+c1.y*c2r;
}

__device__ float3 multiplyComplexRealPartConj(float2 c1, float3 c2r, float3 c2i) {
    return c1.x*c2r+c1.y*c2i;
}

__device__ float3 multiplyComplexImagPartConj(float2 c1, float3 c2r, float3 c2i) {
    return c1.x*c2i-c1.y*c2r;
}

/**
 * Apply the PILE-L thermostat.
 */
extern "C" __global__ void applyPileThermostat(float4* velm, float4* random, unsigned int randomIndex,
        float dt, float kT, float friction) {
    const int numBlocks = blockDim.x*gridDim.x/NUM_COPIES;
    const int blockStart = NUM_COPIES*(threadIdx.x/NUM_COPIES);
    const int indexInBlock = threadIdx.x-blockStart;
    const float nkT = NUM_COPIES*kT;
    const float twown = 2.0f*nkT/HBAR;
    const float c1_0 = EXP(-0.5f*dt*friction);
    const float c2_0 = SQRT(1.0f-c1_0*c1_0);
    __shared__ float3 v[2*THREAD_BLOCK_SIZE];
    __shared__ float3 temp[2*THREAD_BLOCK_SIZE];
    __shared__ float2 w[NUM_COPIES];
    float3* vreal = &v[blockStart];
    float3* vimag = &v[blockStart+blockDim.x];
    if (threadIdx.x < NUM_COPIES)
        w[indexInBlock] = make_float2(cos(-indexInBlock*2*M_PI/NUM_COPIES), sin(-indexInBlock*2*M_PI/NUM_COPIES));
    __syncthreads();
    randomIndex += NUM_COPIES*((blockIdx.x*blockDim.x+threadIdx.x)/NUM_COPIES);
    for (int particle = (blockIdx.x*blockDim.x+threadIdx.x)/NUM_COPIES; particle < NUM_ATOMS; particle += numBlocks) {
        float4 particleVelm = velm[particle+indexInBlock*PADDED_NUM_ATOMS];
        float invMass = particleVelm.w;
        float c3_0 = c2_0*SQRT(nkT*invMass);
        
        // Forward FFT.
        
        vreal[indexInBlock] = SCALE*make_float3(particleVelm.x, particleVelm.y, particleVelm.z);
        vimag[indexInBlock] = make_float3(0);
        __syncthreads();
        FFT_V_FORWARD

        // Apply the thermostat.

        if (indexInBlock == 0) {
            // Apply a local Langevin thermostat to the centroid mode.

            float4 rand = random[randomIndex];
            vreal[0] = vreal[0]*c1_0 + c3_0*make_float3(rand.x, rand.y, rand.z);
        }
        else {
            // Use critical damping white noise for the remaining modes.

            int k = (indexInBlock <= NUM_COPIES/2 ? indexInBlock : NUM_COPIES-indexInBlock);
            const bool isCenter = (NUM_COPIES%2 == 0 && k == NUM_COPIES/2);
            const float wk = twown*sin(k*M_PI/NUM_COPIES);
            const float c1 = EXP(-wk*dt);
            const float c2 = SQRT((1.0f-c1*c1)/2.0f) * (isCenter ? sqrt(2.0f) : 1.0f);
            const float c3 = c2*SQRT(nkT*invMass);
            float4 rand1 = c3*random[randomIndex+k];
            float4 rand2 = (isCenter ? make_float4(0) : c3*random[randomIndex+NUM_COPIES-k]);
            vreal[indexInBlock] = c1*vreal[indexInBlock] + make_float3(rand1.x, rand1.y, rand1.z);
            vimag[indexInBlock] = c1*vimag[indexInBlock] + (indexInBlock < NUM_COPIES/2 ? make_float3(rand2.x, rand2.y, rand2.z) : make_float3(-rand2.x, -rand2.y, -rand2.z));
        }
        __syncthreads();
        
        // Inverse FFT.
        
        FFT_V_BACKWARD
        velm[particle+indexInBlock*PADDED_NUM_ATOMS] = make_float4(SCALE*vreal[indexInBlock].x, SCALE*vreal[indexInBlock].y, SCALE*vreal[indexInBlock].z, particleVelm.w);
        randomIndex += blockDim.x*gridDim.x;
    }
}

/**
 * Advance the positions and velocities.
 */
extern "C" __global__ void integrateStep(float4* posq, float4* velm, long long* force, float dt, float kT) {
    const int numBlocks = (blockDim.x*gridDim.x)/NUM_COPIES;
    const int blockStart = NUM_COPIES*(threadIdx.x/NUM_COPIES);
    const int indexInBlock = threadIdx.x-blockStart;
    const float nkT = NUM_COPIES*kT;
    const float twown = 2.0f*nkT/HBAR;
    const float forceScale = 1/(float) 0xFFFFFFFF;
    __shared__ float3 q[2*THREAD_BLOCK_SIZE];
    __shared__ float3 v[2*THREAD_BLOCK_SIZE];
    __shared__ float3 temp[2*THREAD_BLOCK_SIZE];
    __shared__ float2 w[NUM_COPIES];

    // Update velocities.
    
    for (int particle = (blockIdx.x*blockDim.x+threadIdx.x)/NUM_COPIES; particle < NUM_ATOMS; particle += numBlocks) {
        int index = particle+indexInBlock*PADDED_NUM_ATOMS;
        int forceIndex = particle+indexInBlock*PADDED_NUM_ATOMS*3;
        float4 particleVelm = velm[index];
        particleVelm.x += forceScale*force[forceIndex]*(0.5f*dt*particleVelm.w);
        particleVelm.y += forceScale*force[forceIndex+PADDED_NUM_ATOMS]*(0.5f*dt*particleVelm.w);
        particleVelm.z += forceScale*force[forceIndex+PADDED_NUM_ATOMS*2]*(0.5f*dt*particleVelm.w);
        velm[index] = particleVelm;
    }
    
    // Evolve the free ring polymer by transforming to the frequency domain.

    float3* qreal = &q[blockStart];
    float3* qimag = &q[blockStart+blockDim.x];
    float3* vreal = &v[blockStart];
    float3* vimag = &v[blockStart+blockDim.x];
    if (threadIdx.x < NUM_COPIES)
        w[indexInBlock] = make_float2(cos(-indexInBlock*2*M_PI/NUM_COPIES), sin(-indexInBlock*2*M_PI/NUM_COPIES));
    __syncthreads();
    for (int particle = (blockIdx.x*blockDim.x+threadIdx.x)/NUM_COPIES; particle < NUM_ATOMS; particle += numBlocks) {
        float4 particlePosq = posq[particle+indexInBlock*PADDED_NUM_ATOMS];
        float4 particleVelm = velm[particle+indexInBlock*PADDED_NUM_ATOMS];
        
        // Forward FFT.
        
        qreal[indexInBlock] = SCALE*make_float3(particlePosq.x, particlePosq.y, particlePosq.z);
        qimag[indexInBlock] = make_float3(0);
        vreal[indexInBlock] = SCALE*make_float3(particleVelm.x, particleVelm.y, particleVelm.z);
        vimag[indexInBlock] = make_float3(0);
        __syncthreads();
        FFT_Q_FORWARD
        FFT_V_FORWARD

        // Apply the thermostat.

        if (indexInBlock == 0) {
            qreal[0] += vreal[0]*dt;
            qimag[0] += vimag[0]*dt;
        }
        else {
            const float wk = twown*sin(indexInBlock*M_PI/NUM_COPIES);
            const float wt = wk*dt;
            const float coswt = cos(wt);
            const float sinwt = sin(wt);
            const float3 vprimereal = vreal[indexInBlock]*coswt - qreal[indexInBlock]*(wk*sinwt); // Advance velocity from t to t+dt
            const float3 vprimeimag = vimag[indexInBlock]*coswt - qimag[indexInBlock]*(wk*sinwt);
            qreal[indexInBlock] = vreal[indexInBlock]*(sinwt/wk) + qreal[indexInBlock]*coswt; // Advance position from t to t+dt
            qimag[indexInBlock] = vimag[indexInBlock]*(sinwt/wk) + qimag[indexInBlock]*coswt;
            vreal[indexInBlock] = vprimereal;
            vimag[indexInBlock] = vprimeimag;
        }
        __syncthreads();
        
        // Inverse FFT.
        
        FFT_Q_BACKWARD
        FFT_V_BACKWARD
        posq[particle+indexInBlock*PADDED_NUM_ATOMS] = make_float4(SCALE*qreal[indexInBlock].x, SCALE*qreal[indexInBlock].y, SCALE*qreal[indexInBlock].z, particlePosq.w);
        velm[particle+indexInBlock*PADDED_NUM_ATOMS] = make_float4(SCALE*vreal[indexInBlock].x, SCALE*vreal[indexInBlock].y, SCALE*vreal[indexInBlock].z, particleVelm.w);
    }
}

/**
 * Advance the velocities by a half step.
 */
extern "C" __global__ void advanceVelocities(float4* velm, long long* force, float dt) {
    const int numBlocks = (blockDim.x*gridDim.x)/NUM_COPIES;
    const int blockStart = NUM_COPIES*(threadIdx.x/NUM_COPIES);
    const int indexInBlock = threadIdx.x-blockStart;
    const float forceScale = 1/(float) 0xFFFFFFFF;

    // Update velocities.
    
    for (int particle = (blockIdx.x*blockDim.x+threadIdx.x)/NUM_COPIES; particle < NUM_ATOMS; particle += numBlocks) {
        int index = particle+indexInBlock*PADDED_NUM_ATOMS;
        int forceIndex = particle+indexInBlock*PADDED_NUM_ATOMS*3;
        float4 particleVelm = velm[index];
        particleVelm.x += forceScale*force[forceIndex]*(0.5f*dt*particleVelm.w);
        particleVelm.y += forceScale*force[forceIndex+PADDED_NUM_ATOMS]*(0.5f*dt*particleVelm.w);
        particleVelm.z += forceScale*force[forceIndex+PADDED_NUM_ATOMS*2]*(0.5f*dt*particleVelm.w);
        velm[index] = particleVelm;
    }
}

/**
 * Copy a set of per-atom values from the integrator's arrays to the context.
 */
extern "C" __global__ void copyToContext(float4* src, float4* dst, int* order, int copy) {
    const int base = copy*PADDED_NUM_ATOMS;
    for (int particle = blockIdx.x*blockDim.x+threadIdx.x; particle < NUM_ATOMS; particle += blockDim.x*gridDim.x) {
        dst[particle] = src[base+order[particle]];
    }
}

/**
 * Copy a set of per-atom force values from the context to the integrator's arrays.
 */
extern "C" __global__ void copyFromContext(long long* src, long long* dst, int* order, int copy) {
    const int base = copy*PADDED_NUM_ATOMS*3;
    for (int particle = blockIdx.x*blockDim.x+threadIdx.x; particle < NUM_ATOMS; particle += blockDim.x*gridDim.x) {
        dst[base+order[particle]] = src[particle];
        dst[base+order[particle]+PADDED_NUM_ATOMS] = src[particle+PADDED_NUM_ATOMS];
        dst[base+order[particle]+PADDED_NUM_ATOMS*2] = src[particle+PADDED_NUM_ATOMS*2];
    }
}

/**
 * Update atom positions so all copies are offset by the same number of periodic box widths.
 */
extern "C" __global__ void applyCellTranslations(float4* posq, float4* movedPos, int* order, int movedCopy) {
    for (int particle = blockIdx.x*blockDim.x+threadIdx.x; particle < NUM_ATOMS; particle += blockDim.x*gridDim.x) {
        int index = order[particle];
        float4 delta = movedPos[particle]-posq[movedCopy*PADDED_NUM_ATOMS+index];
        for (int copy = 0; copy < NUM_COPIES; copy++)
            posq[copy*PADDED_NUM_ATOMS+index] += delta;
    }
}
