#ifdef cl_khr_fp64
#pragma OPENCL EXTENSION cl_khr_fp64 : enable
#endif

float4 multiplyComplexRealPart(float2 c1, float4 c2r, float4 c2i) {
    return c1.x*c2r-c1.y*c2i;
}

float4 multiplyComplexImagPart(float2 c1, float4 c2r, float4 c2i) {
    return c1.x*c2i+c1.y*c2r;
}

float4 multiplyComplexRealPartConj(float2 c1, float4 c2r, float4 c2i) {
    return c1.x*c2r+c1.y*c2i;
}

float4 multiplyComplexImagPartConj(float2 c1, float4 c2r, float4 c2i) {
    return c1.x*c2i-c1.y*c2r;
}

/**
 * Apply the PILE-L thermostat.
 */
__kernel void applyPileThermostat(__global float4* velm, __local float4* v, __local float4* temp, __local float2* w, __global float4* random, unsigned int randomIndex,
        float dt, float kT, float friction) {
    const int numBlocks = get_global_size(0)/NUM_COPIES;
    const int blockStart = NUM_COPIES*(get_local_id(0)/NUM_COPIES);
    const int indexInBlock = get_local_id(0)-blockStart;
    const float nkT = NUM_COPIES*kT;
    const float twown = 2.0f*nkT/HBAR;
    const float c1_0 = EXP(-0.5f*dt*friction);
    const float c2_0 = SQRT(1.0f-c1_0*c1_0);
    __local float4* vreal = &v[blockStart];
    __local float4* vimag = &v[blockStart+get_local_size(0)];
    if (get_local_id(0) < NUM_COPIES)
        w[indexInBlock] = (float2) (cos(-indexInBlock*2*M_PI/NUM_COPIES), sin(-indexInBlock*2*M_PI/NUM_COPIES));
    barrier(CLK_LOCAL_MEM_FENCE);
    randomIndex += blockStart;
    for (int particle = get_global_id(0)/NUM_COPIES; particle < NUM_ATOMS; particle += numBlocks) {
        float4 particleVelm = velm[particle+indexInBlock*PADDED_NUM_ATOMS];
        float invMass = particleVelm.w;
        float c3_0 = c2_0*SQRT(nkT*invMass);
        
        // Forward FFT.
        
        vreal[indexInBlock] = SCALE*particleVelm;
        vimag[indexInBlock] = (float4) (0.0f, 0.0f, 0.0f, 0.0f);
        barrier(CLK_GLOBAL_MEM_FENCE);
        FFT_V_FORWARD

        // Apply the thermostat.

        if (indexInBlock == 0) {
            // Apply a local Langevin thermostat to the centroid mode.

            vreal[0].xyz = vreal[0].xyz*c1_0 + c3_0*random[randomIndex].xyz;
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
            float4 rand2 = (isCenter ? 0.0f : c3*random[randomIndex+NUM_COPIES-k]);
            vreal[indexInBlock].xyz = c1*vreal[indexInBlock].xyz + rand1.xyz;
            vimag[indexInBlock].xyz = c1*vimag[indexInBlock].xyz + (indexInBlock < NUM_COPIES/2 ? rand2.xyz : -rand2.xyz);
        }
        barrier(CLK_GLOBAL_MEM_FENCE);
        
        // Inverse FFT.
        
        FFT_V_BACKWARD
        velm[particle+indexInBlock*PADDED_NUM_ATOMS].xyz = SCALE*vreal[indexInBlock].xyz;
        randomIndex += get_global_size(0);
    }
}

/**
 * Advance the positions and velocities.
 */
__kernel void integrateStep(__global float4* posq, __global float4* velm, __global float4* force,
        __local float4* q, __local float4* v, __local float4* temp, __local float2* w, float dt, float kT) {
    const int numBlocks = get_global_size(0)/NUM_COPIES;
    const int blockStart = NUM_COPIES*(get_local_id(0)/NUM_COPIES);
    const int indexInBlock = get_local_id(0)-blockStart;
    const float nkT = NUM_COPIES*kT;
    const float twown = 2.0f*nkT/HBAR;

    // Update velocities.
    
    for (int particle = get_global_id(0)/NUM_COPIES; particle < NUM_ATOMS; particle += numBlocks) {
        int index = particle+indexInBlock*PADDED_NUM_ATOMS;
        float4 particleVelm = velm[index];
        particleVelm.xyz += force[index].xyz*(0.5f*dt*particleVelm.w);
        velm[index] = particleVelm;
    }
    
    // Evolve the free ring polymer by transforming to the frequency domain.

    __local float4* qreal = &q[blockStart];
    __local float4* qimag = &q[blockStart+get_local_size(0)];
    __local float4* vreal = &v[blockStart];
    __local float4* vimag = &v[blockStart+get_local_size(0)];
    if (get_local_id(0) < NUM_COPIES)
        w[indexInBlock] = (float2) (cos(-indexInBlock*2*M_PI/NUM_COPIES), sin(-indexInBlock*2*M_PI/NUM_COPIES));
    barrier(CLK_LOCAL_MEM_FENCE);
    for (int particle = get_global_id(0)/NUM_COPIES; particle < NUM_ATOMS; particle += numBlocks) {
        float4 particlePosq = posq[particle+indexInBlock*PADDED_NUM_ATOMS];
        float4 particleVelm = velm[particle+indexInBlock*PADDED_NUM_ATOMS];
        float invMass = particleVelm.w;
        
        // Forward FFT.
        
        qreal[indexInBlock] = SCALE*particlePosq;
        qimag[indexInBlock] = (float4) (0.0f, 0.0f, 0.0f, 0.0f);
        vreal[indexInBlock] = SCALE*particleVelm;
        vimag[indexInBlock] = (float4) (0.0f, 0.0f, 0.0f, 0.0f);
        barrier(CLK_GLOBAL_MEM_FENCE);
        FFT_Q_FORWARD
        FFT_V_FORWARD

        // Apply the thermostat.

        if (indexInBlock == 0) {
            qreal[0].xyz += vreal[0].xyz*dt;
            qimag[0].xyz += vimag[0].xyz*dt;
        }
        else {
            const float wk = twown*sin(indexInBlock*M_PI/NUM_COPIES);
            const float wt = wk*dt;
            const float coswt = cos(wt);
            const float sinwt = sin(wt);
            const float wm = wk/particleVelm.w;
            const float4 vprimereal = vreal[indexInBlock]*coswt - qreal[indexInBlock]*(wk*sinwt); // Advance velocity from t to t+dt
            const float4 vprimeimag = vimag[indexInBlock]*coswt - qimag[indexInBlock]*(wk*sinwt);
            qreal[indexInBlock] = vreal[indexInBlock]*(sinwt/wk) + qreal[indexInBlock]*coswt; // Advance position from t to t+dt
            qimag[indexInBlock] = vimag[indexInBlock]*(sinwt/wk) + qimag[indexInBlock]*coswt;
            vreal[indexInBlock] = vprimereal;
            vimag[indexInBlock] = vprimeimag;
        }
        barrier(CLK_GLOBAL_MEM_FENCE);
        
        // Inverse FFT.
        
        FFT_Q_BACKWARD
        FFT_V_BACKWARD
        posq[particle+indexInBlock*PADDED_NUM_ATOMS].xyz = SCALE*qreal[indexInBlock].xyz;
        velm[particle+indexInBlock*PADDED_NUM_ATOMS].xyz = SCALE*vreal[indexInBlock].xyz;
    }
}

/**
 * Advance the velocities by a half step.
 */
__kernel void advanceVelocities(__global float4* velm, __global float4* force, float dt) {
    const int numBlocks = get_global_size(0)/NUM_COPIES;
    const int blockStart = NUM_COPIES*(get_local_id(0)/NUM_COPIES);
    const int indexInBlock = get_local_id(0)-blockStart;

    // Update velocities.
    
    for (int particle = get_global_id(0)/NUM_COPIES; particle < NUM_ATOMS; particle += numBlocks) {
        int index = particle+indexInBlock*PADDED_NUM_ATOMS;
        float4 particleVelm = velm[index];
        particleVelm.xyz += force[index].xyz*(0.5f*dt*particleVelm.w);
        velm[index] = particleVelm;
    }
}

/**
 * Copy a set of per-atom values from the integrator's arrays to the context.
 */
__kernel void copyToContext(__global float4* src, __global float4* dst, __global int* order, int copy) {
    const int base = copy*PADDED_NUM_ATOMS;
    for (int particle = get_global_id(0); particle < NUM_ATOMS; particle += get_global_size(0)) {
        dst[particle] = src[base+order[particle]];
    }
}

/**
 * Copy a set of per-atom values from the context to the integrator's arrays.
 */
__kernel void copyFromContext(__global float4* src, __global float4* dst, __global int* order, int copy) {
    const int base = copy*PADDED_NUM_ATOMS;
    for (int particle = get_global_id(0); particle < NUM_ATOMS; particle += get_global_size(0)) {
        dst[base+order[particle]] = src[particle];
    }
}
