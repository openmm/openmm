/**
 * Perform the first step of Verlet integration.
 */

KERNEL void integrateVerletPart1(int numAtoms, int paddedNumAtoms, GLOBAL const mixed2* RESTRICT dt, GLOBAL const real4* RESTRICT posq,
        GLOBAL mixed4* RESTRICT velm, GLOBAL const mm_long* RESTRICT force, GLOBAL mixed4* RESTRICT posDelta
#ifdef USE_MIXED_PRECISION
        , GLOBAL const real4* RESTRICT posqCorrection
#endif
    ) {
    const mixed2 stepSize = dt[0];
    const mixed dtPos = stepSize.y;
    const mixed dtVel = 0.5f*(stepSize.x+stepSize.y);
    const mixed scale = dtVel/(mixed) 0x100000000;
    for (int index = GLOBAL_ID; index < numAtoms; index += GLOBAL_SIZE) {
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
 * Perform the second step of Verlet integration.
 */

KERNEL void integrateVerletPart2(int numAtoms, GLOBAL mixed2* RESTRICT dt, GLOBAL real4* RESTRICT posq,
        GLOBAL mixed4* RESTRICT velm, GLOBAL const mixed4* RESTRICT posDelta
#ifdef USE_MIXED_PRECISION
        , GLOBAL real4* RESTRICT posqCorrection
#endif
    ) {
    mixed2 stepSize = dt[0];
#ifdef SUPPORTS_DOUBLE_PRECISION
    double oneOverDt = 1.0/stepSize.y;
#else
    float oneOverDt = 1.0f/stepSize.y;
    float correction = (1.0f-oneOverDt*stepSize.y)/stepSize.y;
#endif
    if (GLOBAL_ID == 0)
        dt[0].x = stepSize.y;
    SYNC_THREADS;
    int index = GLOBAL_ID;
    for (; index < numAtoms; index += GLOBAL_SIZE) {
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
#ifdef SUPPORTS_DOUBLE_PRECISION
            velocity = make_mixed4((mixed) (delta.x*oneOverDt), (mixed) (delta.y*oneOverDt), (mixed) (delta.z*oneOverDt), velocity.w);
#else
            velocity = make_mixed4((mixed) (delta.x*oneOverDt+delta.x*correction), (mixed) (delta.y*oneOverDt+delta.y*correction), (mixed) (delta.z*oneOverDt+delta.z*correction), velocity.w);
#endif
#ifdef USE_MIXED_PRECISION
            posq[index] = make_real4((real) pos.x, (real) pos.y, (real) pos.z, (real) pos.w);
            posqCorrection[index] = make_real4(pos.x-(real) pos.x, pos.y-(real) pos.y, pos.z-(real) pos.z, 0);
#else
            posq[index] = pos;
#endif
            velm[index] = velocity;
        }
    }
}

/**
 * Select the step size to use for the next step.
 */

KERNEL void selectVerletStepSize(int numAtoms, int paddedNumAtoms, mixed maxStepSize, mixed errorTol, GLOBAL mixed2* RESTRICT dt, GLOBAL const mixed4* RESTRICT velm, GLOBAL const mm_long* RESTRICT force) {
    // Calculate the error.

    LOCAL mixed error[256];
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
    if (LOCAL_ID == 0) {
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
    }
}
