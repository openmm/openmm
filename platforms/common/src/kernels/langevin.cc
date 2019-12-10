enum {VelScale, ForceScale, NoiseScale, MaxParams};

/**
 * Perform the first step of Langevin integration.
 */

KERNEL void integrateLangevinPart1(int numAtoms, int paddedNumAtoms, GLOBAL mixed4* RESTRICT velm, GLOBAL const mm_long* RESTRICT force, GLOBAL mixed4* RESTRICT posDelta,
        GLOBAL const mixed* RESTRICT paramBuffer, GLOBAL const mixed2* RESTRICT dt, GLOBAL const float4* RESTRICT random, unsigned int randomIndex) {
    mixed vscale = paramBuffer[VelScale];
    mixed fscale = paramBuffer[ForceScale]/(mixed) 0x100000000;
    mixed noisescale = paramBuffer[NoiseScale];
    mixed stepSize = dt[0].y;
    int index = GLOBAL_ID;
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
        randomIndex += GLOBAL_SIZE;
        index += GLOBAL_SIZE;
    }
}

/**
 * Perform the second step of Langevin integration.
 */

KERNEL void integrateLangevinPart2(int numAtoms, GLOBAL real4* RESTRICT posq, GLOBAL const mixed4* RESTRICT posDelta, GLOBAL mixed4* RESTRICT velm, GLOBAL const mixed2* RESTRICT dt
#ifdef USE_MIXED_PRECISION
        , GLOBAL real4* RESTRICT posqCorrection
#endif
        ) {
#ifdef SUPPORTS_DOUBLE_PRECISION
    double invStepSize = 1.0/dt[0].y;
#else
    float invStepSize = 1.0f/dt[0].y;
    float correction = (1.0f-invStepSize*dt[0].y)/dt[0].y;
#endif
    int index = GLOBAL_ID;
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
#ifdef SUPPORTS_DOUBLE_PRECISION
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
        index += GLOBAL_SIZE;
    }
}

/**
 * Select the step size to use for the next step.
 */

KERNEL void selectLangevinStepSize(int numAtoms, int paddedNumAtoms, mixed maxStepSize, mixed errorTol, mixed friction, mixed kT, GLOBAL mixed2* RESTRICT dt,
        GLOBAL const mixed4* RESTRICT velm, GLOBAL const mm_long* RESTRICT force, GLOBAL mixed* RESTRICT paramBuffer) {
    // Calculate the error.

    LOCAL mixed error[256];
    LOCAL mixed params[MaxParams];
    mixed err = 0;
    const mixed scale = RECIP((mixed) 0x100000000);
    for (int index = LOCAL_ID; index < numAtoms; index += LOCAL_SIZE) {
        mixed3 f = make_mixed3(scale*force[index], scale*force[index+paddedNumAtoms], scale*force[index+paddedNumAtoms*2]);
        mixed invMass = velm[index].w;
        err += (f.x*f.x + f.y*f.y + f.z*f.z)*invMass*invMass;
    }
    error[LOCAL_ID] = err;
    SYNC_THREADS;

    // Sum the errors from all threads.

    for (unsigned int offset = 1; offset < LOCAL_SIZE; offset *= 2) {
        if (LOCAL_ID+offset < LOCAL_SIZE && (LOCAL_ID&(2*offset-1)) == 0)
            error[LOCAL_ID] += error[LOCAL_ID+offset];
        SYNC_THREADS;
    }
    if (GLOBAL_ID == 0) {
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
    SYNC_THREADS;
    if (LOCAL_ID < MaxParams)
        paramBuffer[LOCAL_ID] = params[LOCAL_ID];
}
