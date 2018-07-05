mixed4 multiplyComplexRealPart(mixed2 c1, mixed4 c2r, mixed4 c2i) {
    return c1.x*c2r-c1.y*c2i;
}

mixed4 multiplyComplexImagPart(mixed2 c1, mixed4 c2r, mixed4 c2i) {
    return c1.x*c2i+c1.y*c2r;
}

mixed4 multiplyComplexRealPartConj(mixed2 c1, mixed4 c2r, mixed4 c2i) {
    return c1.x*c2r+c1.y*c2i;
}

mixed4 multiplyComplexImagPartConj(mixed2 c1, mixed4 c2r, mixed4 c2i) {
    return c1.x*c2i-c1.y*c2r;
}

/**
 * Apply the PILE-L thermostat.
 */
__kernel void applyPileThermostat(__global mixed4* velm, __global float4* random, unsigned int randomIndex,
        mixed dt, mixed kT, mixed friction) {
    const int numBlocks = get_global_size(0)/NUM_COPIES;
    const int blockStart = NUM_COPIES*(get_local_id(0)/NUM_COPIES);
    const int indexInBlock = get_local_id(0)-blockStart;
    const mixed nkT = NUM_COPIES*kT;
    const mixed twown = 2.0f*nkT/HBAR;
    const mixed c1_0 = exp(-0.5f*dt*friction);
    const mixed c2_0 = sqrt(1.0f-c1_0*c1_0);
    __local mixed4 v[2*THREAD_BLOCK_SIZE];
    __local mixed4 temp[2*THREAD_BLOCK_SIZE];
    __local mixed2 w[NUM_COPIES];
    __local mixed4* vreal = &v[blockStart];
    __local mixed4* vimag = &v[blockStart+get_local_size(0)];
    if (get_local_id(0) < NUM_COPIES)
        w[indexInBlock] = (mixed2) (cos(-indexInBlock*2*M_PI/NUM_COPIES), sin(-indexInBlock*2*M_PI/NUM_COPIES));
    barrier(CLK_LOCAL_MEM_FENCE);
    randomIndex += NUM_COPIES*(get_global_id(0)/NUM_COPIES);
    for (int particle = get_global_id(0)/NUM_COPIES; particle < NUM_ATOMS; particle += numBlocks) {
        mixed4 particleVelm = velm[particle+indexInBlock*PADDED_NUM_ATOMS];
        mixed invMass = particleVelm.w;
        mixed c3_0 = c2_0*sqrt(nkT*invMass);
        
        // Forward FFT.
        
        vreal[indexInBlock] = SCALE*particleVelm;
        vimag[indexInBlock] = (mixed4) (0.0f, 0.0f, 0.0f, 0.0f);
        barrier(CLK_GLOBAL_MEM_FENCE);
        FFT_V_FORWARD

        // Apply the thermostat.

        if (indexInBlock == 0) {
            // Apply a local Langevin thermostat to the centroid mode.

            vreal[0].xyz = vreal[0].xyz*c1_0 + c3_0*convert_mixed4(random[randomIndex]).xyz;
        }
        else {
            // Use critical damping white noise for the remaining modes.

            int k = (indexInBlock <= NUM_COPIES/2 ? indexInBlock : NUM_COPIES-indexInBlock);
            const bool isCenter = (NUM_COPIES%2 == 0 && k == NUM_COPIES/2);
            const mixed wk = twown*sin(k*M_PI/NUM_COPIES);
            const mixed c1 = exp(-wk*dt);
            const mixed c2 = sqrt((1.0f-c1*c1)/2.0f) * (isCenter ? sqrt(2.0f) : 1.0f);
            const mixed c3 = c2*sqrt(nkT*invMass);
            mixed4 rand1 = c3*convert_mixed4(random[randomIndex+k]);
            mixed4 rand2 = (isCenter ? 0.0f : c3*convert_mixed4(random[randomIndex+NUM_COPIES-k]));
            vreal[indexInBlock].xyz = c1*vreal[indexInBlock].xyz + rand1.xyz;
            vimag[indexInBlock].xyz = c1*vimag[indexInBlock].xyz + (indexInBlock < NUM_COPIES/2 ? rand2.xyz : -rand2.xyz);
        }
        barrier(CLK_GLOBAL_MEM_FENCE);
        
        // Inverse FFT.
        
        FFT_V_BACKWARD
        if (invMass != 0)
            velm[particle+indexInBlock*PADDED_NUM_ATOMS].xyz = SCALE*vreal[indexInBlock].xyz;
        randomIndex += get_global_size(0);
    }
}

/**
 * Advance the positions and velocities.
 */
__kernel void integrateStep(__global mixed4* posq, __global mixed4* velm, __global real4* force, mixed dt, mixed kT) {
    const int numBlocks = get_global_size(0)/NUM_COPIES;
    const int blockStart = NUM_COPIES*(get_local_id(0)/NUM_COPIES);
    const int indexInBlock = get_local_id(0)-blockStart;
    const mixed nkT = NUM_COPIES*kT;
    const mixed twown = 2.0f*nkT/HBAR;
    __local mixed4 q[2*THREAD_BLOCK_SIZE];
    __local mixed4 v[2*THREAD_BLOCK_SIZE];
    __local mixed4 temp[2*THREAD_BLOCK_SIZE];
    __local mixed2 w[NUM_COPIES];

    // Update velocities.
    
    for (int particle = get_global_id(0)/NUM_COPIES; particle < NUM_ATOMS; particle += numBlocks) {
        int index = particle+indexInBlock*PADDED_NUM_ATOMS;
        mixed4 particleVelm = velm[index];
        particleVelm.xyz += convert_mixed4(force[index]).xyz*(0.5f*dt*particleVelm.w);
        if (particleVelm.w != 0)
            velm[index] = particleVelm;
    }
    
    // Evolve the free ring polymer by transforming to the frequency domain.

    __local mixed4* qreal = &q[blockStart];
    __local mixed4* qimag = &q[blockStart+get_local_size(0)];
    __local mixed4* vreal = &v[blockStart];
    __local mixed4* vimag = &v[blockStart+get_local_size(0)];
    if (get_local_id(0) < NUM_COPIES)
        w[indexInBlock] = (mixed2) (cos(-indexInBlock*2*M_PI/NUM_COPIES), sin(-indexInBlock*2*M_PI/NUM_COPIES));
    barrier(CLK_LOCAL_MEM_FENCE);
    for (int particle = get_global_id(0)/NUM_COPIES; particle < NUM_ATOMS; particle += numBlocks) {
        mixed4 particlePosq = posq[particle+indexInBlock*PADDED_NUM_ATOMS];
        mixed4 particleVelm = velm[particle+indexInBlock*PADDED_NUM_ATOMS];
        
        // Forward FFT.
        
        qreal[indexInBlock] = SCALE*particlePosq;
        qimag[indexInBlock] = (mixed4) (0.0f, 0.0f, 0.0f, 0.0f);
        vreal[indexInBlock] = SCALE*particleVelm;
        vimag[indexInBlock] = (mixed4) (0.0f, 0.0f, 0.0f, 0.0f);
        barrier(CLK_GLOBAL_MEM_FENCE);
        FFT_Q_FORWARD
        FFT_V_FORWARD

        // Apply the thermostat.

        if (indexInBlock == 0) {
            qreal[0].xyz += vreal[0].xyz*dt;
            qimag[0].xyz += vimag[0].xyz*dt;
        }
        else {
            const mixed wk = twown*sin(indexInBlock*M_PI/NUM_COPIES);
            const mixed wt = wk*dt;
            const mixed coswt = cos(wt);
            const mixed sinwt = sin(wt);
            const mixed4 vprimereal = vreal[indexInBlock]*coswt - qreal[indexInBlock]*(wk*sinwt); // Advance velocity from t to t+dt
            const mixed4 vprimeimag = vimag[indexInBlock]*coswt - qimag[indexInBlock]*(wk*sinwt);
            qreal[indexInBlock] = vreal[indexInBlock]*(sinwt/wk) + qreal[indexInBlock]*coswt; // Advance position from t to t+dt
            qimag[indexInBlock] = vimag[indexInBlock]*(sinwt/wk) + qimag[indexInBlock]*coswt;
            vreal[indexInBlock] = vprimereal;
            vimag[indexInBlock] = vprimeimag;
        }
        barrier(CLK_GLOBAL_MEM_FENCE);
        
        // Inverse FFT.
        
        FFT_Q_BACKWARD
        FFT_V_BACKWARD
        if (particleVelm.w != 0) {
            posq[particle+indexInBlock*PADDED_NUM_ATOMS].xyz = SCALE*qreal[indexInBlock].xyz;
            velm[particle+indexInBlock*PADDED_NUM_ATOMS].xyz = SCALE*vreal[indexInBlock].xyz;
        }
    }
}

/**
 * Advance the velocities by a half step.
 */
__kernel void advanceVelocities(__global mixed4* velm, __global real4* force, mixed dt) {
    const int numBlocks = get_global_size(0)/NUM_COPIES;
    const int blockStart = NUM_COPIES*(get_local_id(0)/NUM_COPIES);
    const int indexInBlock = get_local_id(0)-blockStart;

    // Update velocities.
    
    for (int particle = get_global_id(0)/NUM_COPIES; particle < NUM_ATOMS; particle += numBlocks) {
        int index = particle+indexInBlock*PADDED_NUM_ATOMS;
        mixed4 particleVelm = velm[index];
        particleVelm.xyz += convert_mixed4(force[index]).xyz*(0.5f*dt*particleVelm.w);
        if (particleVelm.w != 0)
            velm[index] = particleVelm;
    }
}

/**
 * Copy a set of positions and velocities from the integrator's arrays to the context.
 */
__kernel void copyDataToContext(__global mixed4* srcVel, __global mixed4* dstVel, __global mixed4* srcPos,
        __global real4* dstPos, __global int* order, int copy) {
    const int base = copy*PADDED_NUM_ATOMS;
    for (int particle = get_global_id(0); particle < NUM_ATOMS; particle += get_global_size(0)) {
        int index = base+order[particle];
        dstVel[particle] = srcVel[index];
        dstPos[particle].xyz = convert_real4(srcPos[index]).xyz;
    }
}

/**
 * Copy a set of positions, velocities, and forces from the context to the integrator's arrays.
 */
__kernel void copyDataFromContext(__global real4* srcForce, __global real4* dstForce, __global mixed4* srcVel,
        __global mixed4* dstVel, __global real4* srcPos, __global mixed4* dstPos, __global int* order, int copy) {
    const int base = copy*PADDED_NUM_ATOMS;
    for (int particle = get_global_id(0); particle < NUM_ATOMS; particle += get_global_size(0)) {
        int index = base+order[particle];
        dstForce[index] = srcForce[particle];
        dstVel[index] = srcVel[particle];
        dstPos[index].xyz = convert_mixed4(srcPos[particle]).xyz;
    }
}

/**
 * Atom positions in one copy have been modified.  Apply the same offsets to all the other copies.
 */
__kernel void applyCellTranslations(__global mixed4* posq, __global real4* movedPos, __global int* order, int movedCopy) {
    for (int particle = get_global_id(0); particle < NUM_ATOMS; particle += get_global_size(0)) {
        int index = order[particle];
        mixed4 delta = convert_mixed4(movedPos[particle])-posq[movedCopy*PADDED_NUM_ATOMS+index];
        for (int copy = 0; copy < NUM_COPIES; copy++)
            posq[copy*PADDED_NUM_ATOMS+index].xyz += delta.xyz;
    }
}
