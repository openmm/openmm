__device__ mixed3 multiplyComplexRealPart(mixed2 c1, mixed3 c2r, mixed3 c2i) {
    return c1.x*c2r-c1.y*c2i;
}

__device__ mixed3 multiplyComplexImagPart(mixed2 c1, mixed3 c2r, mixed3 c2i) {
    return c1.x*c2i+c1.y*c2r;
}

__device__ mixed3 multiplyComplexRealPartConj(mixed2 c1, mixed3 c2r, mixed3 c2i) {
    return c1.x*c2r+c1.y*c2i;
}

__device__ mixed3 multiplyComplexImagPartConj(mixed2 c1, mixed3 c2r, mixed3 c2i) {
    return c1.x*c2i-c1.y*c2r;
}

/**
 * Apply the PILE-L thermostat.
 */
extern "C" __global__ void applyPileThermostat(mixed4* velm, float4* random, unsigned int randomIndex,
        mixed dt, mixed kT, mixed friction) {
    const int numBlocks = blockDim.x*gridDim.x/NUM_COPIES;
    const int blockStart = NUM_COPIES*(threadIdx.x/NUM_COPIES);
    const int indexInBlock = threadIdx.x-blockStart;
    const mixed nkT = NUM_COPIES*kT;
    const mixed twown = 2.0f*nkT/HBAR;
    const mixed c1_0 = EXP(-0.5f*dt*friction);
    const mixed c2_0 = SQRT(1.0f-c1_0*c1_0);
    __shared__ mixed3 v[2*THREAD_BLOCK_SIZE];
    __shared__ mixed3 temp[2*THREAD_BLOCK_SIZE];
    __shared__ mixed2 w[NUM_COPIES];
    mixed3* vreal = &v[blockStart];
    mixed3* vimag = &v[blockStart+blockDim.x];
    if (threadIdx.x < NUM_COPIES)
        w[indexInBlock] = make_mixed2(cos(-indexInBlock*2*M_PI/NUM_COPIES), sin(-indexInBlock*2*M_PI/NUM_COPIES));
    __syncthreads();
    randomIndex += NUM_COPIES*((blockIdx.x*blockDim.x+threadIdx.x)/NUM_COPIES);
    for (int particle = (blockIdx.x*blockDim.x+threadIdx.x)/NUM_COPIES; particle < NUM_ATOMS; particle += numBlocks) {
        mixed4 particleVelm = velm[particle+indexInBlock*PADDED_NUM_ATOMS];
        mixed invMass = particleVelm.w;
        mixed c3_0 = c2_0*SQRT(nkT*invMass);
        
        // Forward FFT.
        
        vreal[indexInBlock] = SCALE*make_mixed3(particleVelm.x, particleVelm.y, particleVelm.z);
        vimag[indexInBlock] = make_mixed3(0);
        __syncthreads();
        FFT_V_FORWARD

        // Apply the thermostat.

        if (indexInBlock == 0) {
            // Apply a local Langevin thermostat to the centroid mode.

            float4 rand = random[randomIndex];
            vreal[0] = vreal[0]*c1_0 + c3_0*make_mixed3(rand.x, rand.y, rand.z);
        }
        else {
            // Use critical damping white noise for the remaining modes.

            int k = (indexInBlock <= NUM_COPIES/2 ? indexInBlock : NUM_COPIES-indexInBlock);
            const bool isCenter = (NUM_COPIES%2 == 0 && k == NUM_COPIES/2);
            const mixed wk = twown*sin(k*M_PI/NUM_COPIES);
            const mixed c1 = EXP(-wk*dt);
            const mixed c2 = SQRT((1.0f-c1*c1)/2.0f) * (isCenter ? sqrt(2.0f) : 1.0f);
            const mixed c3 = c2*SQRT(nkT*invMass);
            float4 rand1 = random[randomIndex+k];
            float4 rand2 = (isCenter ? make_float4(0) : random[randomIndex+NUM_COPIES-k]);
            vreal[indexInBlock] = c1*vreal[indexInBlock] + c3*make_mixed3(rand1.x, rand1.y, rand1.z);
            vimag[indexInBlock] = c1*vimag[indexInBlock] + c3*(indexInBlock < NUM_COPIES/2 ? make_mixed3(rand2.x, rand2.y, rand2.z) : make_mixed3(-rand2.x, -rand2.y, -rand2.z));
        }
        __syncthreads();
        
        // Inverse FFT.
        
        FFT_V_BACKWARD
        if (invMass != 0)
            velm[particle+indexInBlock*PADDED_NUM_ATOMS] = make_mixed4(SCALE*vreal[indexInBlock].x, SCALE*vreal[indexInBlock].y, SCALE*vreal[indexInBlock].z, particleVelm.w);
        randomIndex += blockDim.x*gridDim.x;
    }
}

/**
 * Advance the positions and velocities.
 */
extern "C" __global__ void integrateStep(mixed4* posq, mixed4* velm, long long* force, mixed dt, mixed kT) {
    const int numBlocks = (blockDim.x*gridDim.x)/NUM_COPIES;
    const int blockStart = NUM_COPIES*(threadIdx.x/NUM_COPIES);
    const int indexInBlock = threadIdx.x-blockStart;
    const mixed nkT = NUM_COPIES*kT;
    const mixed twown = 2.0f*nkT/HBAR;
    const mixed forceScale = 1/(mixed) 0x100000000;
    __shared__ mixed3 q[2*THREAD_BLOCK_SIZE];
    __shared__ mixed3 v[2*THREAD_BLOCK_SIZE];
    __shared__ mixed3 temp[2*THREAD_BLOCK_SIZE];
    __shared__ mixed2 w[NUM_COPIES];

    // Update velocities.
    
    for (int particle = (blockIdx.x*blockDim.x+threadIdx.x)/NUM_COPIES; particle < NUM_ATOMS; particle += numBlocks) {
        int index = particle+indexInBlock*PADDED_NUM_ATOMS;
        int forceIndex = particle+indexInBlock*PADDED_NUM_ATOMS*3;
        mixed4 particleVelm = velm[index];
        particleVelm.x += forceScale*force[forceIndex]*(0.5f*dt*particleVelm.w);
        particleVelm.y += forceScale*force[forceIndex+PADDED_NUM_ATOMS]*(0.5f*dt*particleVelm.w);
        particleVelm.z += forceScale*force[forceIndex+PADDED_NUM_ATOMS*2]*(0.5f*dt*particleVelm.w);
        if (particleVelm.w != 0)
            velm[index] = particleVelm;
    }
    
    // Evolve the free ring polymer by transforming to the frequency domain.

    mixed3* qreal = &q[blockStart];
    mixed3* qimag = &q[blockStart+blockDim.x];
    mixed3* vreal = &v[blockStart];
    mixed3* vimag = &v[blockStart+blockDim.x];
    if (threadIdx.x < NUM_COPIES)
        w[indexInBlock] = make_mixed2(cos(-indexInBlock*2*M_PI/NUM_COPIES), sin(-indexInBlock*2*M_PI/NUM_COPIES));
    __syncthreads();
    for (int particle = (blockIdx.x*blockDim.x+threadIdx.x)/NUM_COPIES; particle < NUM_ATOMS; particle += numBlocks) {
        mixed4 particlePosq = posq[particle+indexInBlock*PADDED_NUM_ATOMS];
        mixed4 particleVelm = velm[particle+indexInBlock*PADDED_NUM_ATOMS];
        
        // Forward FFT.
        
        qreal[indexInBlock] = SCALE*make_mixed3(particlePosq.x, particlePosq.y, particlePosq.z);
        qimag[indexInBlock] = make_mixed3(0);
        vreal[indexInBlock] = SCALE*make_mixed3(particleVelm.x, particleVelm.y, particleVelm.z);
        vimag[indexInBlock] = make_mixed3(0);
        __syncthreads();
        FFT_Q_FORWARD
        FFT_V_FORWARD

        // Apply the thermostat.

        if (indexInBlock == 0) {
            qreal[0] += vreal[0]*dt;
            qimag[0] += vimag[0]*dt;
        }
        else {
            const mixed wk = twown*sin(indexInBlock*M_PI/NUM_COPIES);
            const mixed wt = wk*dt;
            const mixed coswt = cos(wt);
            const mixed sinwt = sin(wt);
            const mixed3 vprimereal = vreal[indexInBlock]*coswt - qreal[indexInBlock]*(wk*sinwt); // Advance velocity from t to t+dt
            const mixed3 vprimeimag = vimag[indexInBlock]*coswt - qimag[indexInBlock]*(wk*sinwt);
            qreal[indexInBlock] = vreal[indexInBlock]*(sinwt/wk) + qreal[indexInBlock]*coswt; // Advance position from t to t+dt
            qimag[indexInBlock] = vimag[indexInBlock]*(sinwt/wk) + qimag[indexInBlock]*coswt;
            vreal[indexInBlock] = vprimereal;
            vimag[indexInBlock] = vprimeimag;
        }
        __syncthreads();
        
        // Inverse FFT.
        
        FFT_Q_BACKWARD
        FFT_V_BACKWARD
        if (particleVelm.w != 0) {
            posq[particle+indexInBlock*PADDED_NUM_ATOMS] = make_mixed4(SCALE*qreal[indexInBlock].x, SCALE*qreal[indexInBlock].y, SCALE*qreal[indexInBlock].z, particlePosq.w);
            velm[particle+indexInBlock*PADDED_NUM_ATOMS] = make_mixed4(SCALE*vreal[indexInBlock].x, SCALE*vreal[indexInBlock].y, SCALE*vreal[indexInBlock].z, particleVelm.w);
        }
    }
}

/**
 * Advance the velocities by a half step.
 */
extern "C" __global__ void advanceVelocities(mixed4* velm, long long* force, mixed dt) {
    const int numBlocks = (blockDim.x*gridDim.x)/NUM_COPIES;
    const int blockStart = NUM_COPIES*(threadIdx.x/NUM_COPIES);
    const int indexInBlock = threadIdx.x-blockStart;
    const mixed forceScale = 1/(mixed) 0x100000000;

    // Update velocities.
    
    for (int particle = (blockIdx.x*blockDim.x+threadIdx.x)/NUM_COPIES; particle < NUM_ATOMS; particle += numBlocks) {
        int index = particle+indexInBlock*PADDED_NUM_ATOMS;
        int forceIndex = particle+indexInBlock*PADDED_NUM_ATOMS*3;
        mixed4 particleVelm = velm[index];
        particleVelm.x += forceScale*force[forceIndex]*(0.5f*dt*particleVelm.w);
        particleVelm.y += forceScale*force[forceIndex+PADDED_NUM_ATOMS]*(0.5f*dt*particleVelm.w);
        particleVelm.z += forceScale*force[forceIndex+PADDED_NUM_ATOMS*2]*(0.5f*dt*particleVelm.w);
        if (particleVelm.w != 0)
            velm[index] = particleVelm;
    }
}

/**
 * Copy a set of positions and velocities from the integrator's arrays to the context.
 */
extern "C" __global__ void copyDataToContext(mixed4* srcVel, mixed4* dstVel, mixed4* srcPos, real4* dstPos, int* order, int copy) {
    const int base = copy*PADDED_NUM_ATOMS;
    for (int particle = blockIdx.x*blockDim.x+threadIdx.x; particle < NUM_ATOMS; particle += blockDim.x*gridDim.x) {
        int index = base+order[particle];
        dstVel[particle] = srcVel[index];
        mixed4 posq = srcPos[index];
        dstPos[particle].x = posq.x;
        dstPos[particle].y = posq.y;
        dstPos[particle].z = posq.z;
    }
}

/**
 * Copy a set of positions, velocities, and forces from the context to the integrator's arrays.
 */
extern "C" __global__ void copyDataFromContext(long long* srcForce, long long* dstForce, mixed4* srcVel, mixed4* dstVel,
        real4* srcPos, mixed4* dstPos, int* order, int copy) {
    const int base = copy*PADDED_NUM_ATOMS;
    for (int particle = blockIdx.x*blockDim.x+threadIdx.x; particle < NUM_ATOMS; particle += blockDim.x*gridDim.x) {
        int index = order[particle];
        dstForce[base*3+index] = srcForce[particle];
        dstForce[base*3+index+PADDED_NUM_ATOMS] = srcForce[particle+PADDED_NUM_ATOMS];
        dstForce[base*3+index+PADDED_NUM_ATOMS*2] = srcForce[particle+PADDED_NUM_ATOMS*2];
        dstVel[base+index] = srcVel[particle];
        real4 posq = srcPos[particle];
        dstPos[base+index].x = posq.x;
        dstPos[base+index].y = posq.y;
        dstPos[base+index].z = posq.z;

    }
}

/**
 * Atom positions in one copy have been modified.  Apply the same offsets to all the other copies.
 */
extern "C" __global__ void applyCellTranslations(mixed4* posq, real4* movedPos, int* order, int movedCopy) {
    for (int particle = blockIdx.x*blockDim.x+threadIdx.x; particle < NUM_ATOMS; particle += blockDim.x*gridDim.x) {
        int index = order[particle];
        real4 p = movedPos[particle];
        mixed4 delta = make_mixed4(p.x, p.y, p.z, p.w)-posq[movedCopy*PADDED_NUM_ATOMS+index];
        for (int copy = 0; copy < NUM_COPIES; copy++) {
            posq[copy*PADDED_NUM_ATOMS+index].x += delta.x;
            posq[copy*PADDED_NUM_ATOMS+index].y += delta.y;
            posq[copy*PADDED_NUM_ATOMS+index].z += delta.z;
        }
    }
}
