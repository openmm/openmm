enum {VelScale, NoiseScale, MaxParams};

/**
 * Perform the first part of integration: velocity step.
 */

KERNEL void integrateLangevinMiddlePart1(int numAtoms, int paddedNumAtoms, GLOBAL mixed4* RESTRICT velm, GLOBAL const mm_long* RESTRICT force,
        GLOBAL const mixed2* RESTRICT dt) {
    mixed fscale = dt[0].y/(mixed) 0x100000000;
    for (int index = GLOBAL_ID; index < numAtoms; index += GLOBAL_SIZE) {
        mixed4 velocity = velm[index];
        if (velocity.w != 0.0) {
            velocity.x += fscale*velocity.w*force[index];
            velocity.y += fscale*velocity.w*force[index+paddedNumAtoms];
            velocity.z += fscale*velocity.w*force[index+paddedNumAtoms*2];
            velm[index] = velocity;
        }
    }
}

/**
 * Perform the second part of integration: position half step, then interact with heat bath,
 * then another position half step.
 */

KERNEL void integrateLangevinMiddlePart2(int numAtoms, GLOBAL mixed4* RESTRICT velm, GLOBAL mixed4* RESTRICT posDelta,
        GLOBAL mixed4* RESTRICT oldDelta, GLOBAL const mixed* RESTRICT paramBuffer, GLOBAL const mixed2* RESTRICT dt, GLOBAL const float4* RESTRICT random, unsigned int startRandomIndex
        ) {
    mixed vscale = paramBuffer[VelScale];
    mixed noisescale = paramBuffer[NoiseScale];
    mixed halfdt = 0.5f*dt[0].y;
    int index = GLOBAL_ID;
    unsigned int randomIndex = startRandomIndex+index;
    while (index < numAtoms) {
        mixed4 velocity = velm[index];
        if (velocity.w != 0.0) {
            mixed4 delta = make_mixed4(halfdt*velocity.x, halfdt*velocity.y, halfdt*velocity.z, 0);
            mixed sqrtInvMass = SQRT(velocity.w);
            velocity.x = vscale*velocity.x + noisescale*sqrtInvMass*random[randomIndex].x;
            velocity.y = vscale*velocity.y + noisescale*sqrtInvMass*random[randomIndex].y;
            velocity.z = vscale*velocity.z + noisescale*sqrtInvMass*random[randomIndex].z;
            velm[index] = velocity;
            delta += make_mixed4(halfdt*velocity.x, halfdt*velocity.y, halfdt*velocity.z, 0);
            posDelta[index] = delta;
            oldDelta[index] = delta;
        }
        randomIndex += GLOBAL_SIZE;
        index += GLOBAL_SIZE;
    }
}

/**
 * Perform the third part of integration: apply constraint forces to velocities, then record
 * the constrained positions.
 */

KERNEL void integrateLangevinMiddlePart3(int numAtoms, GLOBAL real4* RESTRICT posq, GLOBAL mixed4* RESTRICT velm,
         GLOBAL mixed4* RESTRICT posDelta, GLOBAL mixed4* RESTRICT oldDelta, GLOBAL const mixed2* RESTRICT dt
#ifdef USE_MIXED_PRECISION
        , GLOBAL real4* RESTRICT posqCorrection
#endif
        ) {
    mixed invDt = 1/dt[0].y;
    for (int index = GLOBAL_ID; index < numAtoms; index += GLOBAL_SIZE) {
        mixed4 velocity = velm[index];
        if (velocity.w != 0.0) {
            mixed4 delta = posDelta[index];
            velocity.x += (delta.x-oldDelta[index].x)*invDt;
            velocity.y += (delta.y-oldDelta[index].y)*invDt;
            velocity.z += (delta.z-oldDelta[index].z)*invDt;
            velm[index] = velocity;
#ifdef USE_MIXED_PRECISION
            real4 pos1 = posq[index];
            real4 pos2 = posqCorrection[index];
            mixed4 pos = make_mixed4(pos1.x+(mixed)pos2.x, pos1.y+(mixed)pos2.y, pos1.z+(mixed)pos2.z, pos1.w);
#else
            real4 pos = posq[index];
#endif
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
        mixed noisescale = sqrt(kT*(1-vscale*vscale));
        params[VelScale] = vscale;
        params[NoiseScale] = noisescale;
    }
    SYNC_THREADS;
    if (LOCAL_ID < MaxParams)
        paramBuffer[LOCAL_ID] = params[LOCAL_ID];
}

#ifdef ENABLE_FUSED_INTEGRATION

inline DEVICE mixed4 integrateForceStep(mixed4 velocity, int atom, int paddedNumAtoms, mixed fscale, GLOBAL const mm_long* RESTRICT force) {
    if (velocity.w != 0.0) {
        velocity.x += fscale*velocity.w*force[atom];
        velocity.y += fscale*velocity.w*force[atom+paddedNumAtoms];
        velocity.z += fscale*velocity.w*force[atom+paddedNumAtoms*2];
    }
    return velocity;
}

inline DEVICE mixed4 integrateMiddleStep(PRIVATE mixed4* velocityPtr, float4 rand, mixed halfdt, mixed vscale, mixed noisescale) {
    mixed4 velocity = *velocityPtr;
    mixed4 delta = make_mixed4(halfdt*velocity.x, halfdt*velocity.y, halfdt*velocity.z, 0);
    mixed sqrtInvMass = SQRT(velocity.w);
    velocity.x = vscale*velocity.x + noisescale*sqrtInvMass*rand.x;
    velocity.y = vscale*velocity.y + noisescale*sqrtInvMass*rand.y;
    velocity.z = vscale*velocity.z + noisescale*sqrtInvMass*rand.z;
    *velocityPtr = velocity;
    return delta + make_mixed4(halfdt*velocity.x, halfdt*velocity.y, halfdt*velocity.z, 0);
}

inline DEVICE void commitConstrainedAtom(int atom, mixed4 velocity, mixed4 pos, mixed4 delta, mixed4 oldDelta, mixed invDt,
        GLOBAL mixed4* RESTRICT velm, GLOBAL real4* RESTRICT posq, GLOBAL real4* RESTRICT posqCorrection) {
    velocity.x += (delta.x-oldDelta.x)*invDt;
    velocity.y += (delta.y-oldDelta.y)*invDt;
    velocity.z += (delta.z-oldDelta.z)*invDt;
    velm[atom] = velocity;
    pos.x += delta.x;
    pos.y += delta.y;
    pos.z += delta.z;
    storePos(posq, posqCorrection, atom, pos);
}

/**
 * Perform a complete BAOAB step for SETTLE clusters, SHAKE clusters, and unconstrained atoms in a single kernel launch.
 */

KERNEL void integrateLangevinMiddleFused(int numUnconstrained, int paddedNumAtoms, int numSettle, int numShake, mixed tol,
        GLOBAL mixed4* RESTRICT velm, GLOBAL const mm_long* RESTRICT force, GLOBAL real4* RESTRICT posq,
        GLOBAL const mixed* RESTRICT paramBuffer, GLOBAL const mixed2* RESTRICT dt,
        GLOBAL const float4* RESTRICT random, unsigned int startRandomIndex,
        GLOBAL const int* RESTRICT unconstrainedAtoms, GLOBAL const int4* RESTRICT settleAtoms,
        GLOBAL const float2* RESTRICT settleParams, GLOBAL const int4* RESTRICT shakeAtoms,
        GLOBAL const float4* RESTRICT shakeParams
#ifdef USE_MIXED_PRECISION
        , GLOBAL real4* RESTRICT posqCorrection
#endif
        ) {
#ifndef USE_MIXED_PRECISION
    GLOBAL real4* posqCorrection = 0;
#endif
    mixed dtVal = dt[0].y;
    mixed fscale = dtVal/(mixed) 0x100000000;
    mixed vscale = paramBuffer[VelScale];
    mixed noisescale = paramBuffer[NoiseScale];
    mixed halfdt = 0.5f*dtVal;
    mixed invDt = 1/dtVal;
    int totalClusters = numSettle + numShake;

    for (int index = GLOBAL_ID; index < totalClusters; index += GLOBAL_SIZE) {
        if (index < numSettle) {
            int4 atoms = settleAtoms[index];
            float2 params = settleParams[index];
            mixed4 apos0 = loadPos(posq, posqCorrection, atoms.x);
            mixed4 apos1 = loadPos(posq, posqCorrection, atoms.y);
            mixed4 apos2 = loadPos(posq, posqCorrection, atoms.z);
            mixed4 v0 = integrateForceStep(velm[atoms.x], atoms.x, paddedNumAtoms, fscale, force);
            mixed4 v1 = integrateForceStep(velm[atoms.y], atoms.y, paddedNumAtoms, fscale, force);
            mixed4 v2 = integrateForceStep(velm[atoms.z], atoms.z, paddedNumAtoms, fscale, force);
            applySettleToVelocitiesCluster(apos0, apos1, apos2, &v0, &v1, &v2);
            mixed4 xp0 = integrateMiddleStep(&v0, random[startRandomIndex+atoms.x], halfdt, vscale, noisescale);
            mixed4 xp1 = integrateMiddleStep(&v1, random[startRandomIndex+atoms.y], halfdt, vscale, noisescale);
            mixed4 xp2 = integrateMiddleStep(&v2, random[startRandomIndex+atoms.z], halfdt, vscale, noisescale);
            mixed4 oldXp0 = xp0, oldXp1 = xp1, oldXp2 = xp2;
            applySettleToPositionsCluster(params, apos0, apos1, apos2, 1/v0.w, 1/v1.w, 1/v2.w, &xp0, &xp1, &xp2);
            commitConstrainedAtom(atoms.x, v0, apos0, xp0, oldXp0, invDt, velm, posq, posqCorrection);
            commitConstrainedAtom(atoms.y, v1, apos1, xp1, oldXp1, invDt, velm, posq, posqCorrection);
            commitConstrainedAtom(atoms.z, v2, apos2, xp2, oldXp2, invDt, velm, posq, posqCorrection);
        }
#ifdef HAS_SHAKE_CLUSTERS
        else {
            int shakeIndex = index - numSettle;
            int4 atoms = shakeAtoms[shakeIndex];
            float4 params = shakeParams[shakeIndex];
            mixed4 apos0 = loadPos(posq, posqCorrection, atoms.x);
            mixed4 apos1 = loadPos(posq, posqCorrection, atoms.y);
            mixed4 apos2 = (atoms.z != -1 ? loadPos(posq, posqCorrection, atoms.z) : make_mixed4(0));
            mixed4 apos3 = (atoms.w != -1 ? loadPos(posq, posqCorrection, atoms.w) : make_mixed4(0));
            mixed4 v0 = integrateForceStep(velm[atoms.x], atoms.x, paddedNumAtoms, fscale, force);
            mixed4 v1 = integrateForceStep(velm[atoms.y], atoms.y, paddedNumAtoms, fscale, force);
            mixed4 v2 = (atoms.z != -1 ? integrateForceStep(velm[atoms.z], atoms.z, paddedNumAtoms, fscale, force) : make_mixed4(0));
            mixed4 v3 = (atoms.w != -1 ? integrateForceStep(velm[atoms.w], atoms.w, paddedNumAtoms, fscale, force) : make_mixed4(0));
            applyShakeToVelocitiesCluster(atoms, params, tol, apos0, apos1, apos2, apos3, &v0, &v1, &v2, &v3);
            mixed4 xp0 = integrateMiddleStep(&v0, random[startRandomIndex+atoms.x], halfdt, vscale, noisescale);
            mixed4 xp1 = integrateMiddleStep(&v1, random[startRandomIndex+atoms.y], halfdt, vscale, noisescale);
            mixed4 xp2 = (atoms.z != -1 ? integrateMiddleStep(&v2, random[startRandomIndex+atoms.z], halfdt, vscale, noisescale) : make_mixed4(0));
            mixed4 xp3 = (atoms.w != -1 ? integrateMiddleStep(&v3, random[startRandomIndex+atoms.w], halfdt, vscale, noisescale) : make_mixed4(0));
            mixed4 oldXp0 = xp0, oldXp1 = xp1, oldXp2 = xp2, oldXp3 = xp3;
            applyShakeToPositionsCluster(atoms, params, tol, apos0, apos1, apos2, apos3, &xp0, &xp1, &xp2, &xp3);
            commitConstrainedAtom(atoms.x, v0, apos0, xp0, oldXp0, invDt, velm, posq, posqCorrection);
            commitConstrainedAtom(atoms.y, v1, apos1, xp1, oldXp1, invDt, velm, posq, posqCorrection);
            if (atoms.z != -1)
                commitConstrainedAtom(atoms.z, v2, apos2, xp2, oldXp2, invDt, velm, posq, posqCorrection);
            if (atoms.w != -1)
                commitConstrainedAtom(atoms.w, v3, apos3, xp3, oldXp3, invDt, velm, posq, posqCorrection);
        }
#endif
    }
    for (int i = GLOBAL_ID; i < numUnconstrained; i += GLOBAL_SIZE) {
        int atom = (totalClusters == 0 ? i : unconstrainedAtoms[i]);
        mixed4 velocity = velm[atom];
        if (velocity.w != 0.0) {
            velocity = integrateForceStep(velocity, atom, paddedNumAtoms, fscale, force);
            mixed4 delta = integrateMiddleStep(&velocity, random[startRandomIndex+atom], halfdt, vscale, noisescale);
            velm[atom] = velocity;
            mixed4 pos = loadPos(posq, posqCorrection, atom);
            pos.x += delta.x;
            pos.y += delta.y;
            pos.z += delta.z;
            storePos(posq, posqCorrection, atom, pos);
        }
    }
}

#endif
