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

/**
 * Perform the first part of integration for a single atom: full-step velocity kick
 *
 *   vel += forceScale * vel.w * force
 *
 * where vel.w = 1/m (0 for massless particles) and forceScale = stepSize / 0x100000000 converts
 * 64-bit fixed-point force components at force[atom + k*paddedNumAtoms].
 */

inline DEVICE mixed4 integrateForceStep(mixed4 vel, int atom, int paddedNumAtoms, mixed forceScale, GLOBAL const mm_long* RESTRICT force) {
    vel.x += forceScale*vel.w*force[atom];
    vel.y += forceScale*vel.w*force[atom+paddedNumAtoms];
    vel.z += forceScale*vel.w*force[atom+paddedNumAtoms*2];
    return vel;
}

/**
 * Perform the second part of integration for a single atom: first position half-step, stochastic
 * heat-bath velocity update, and second position half-step, returning the total unconstrained
 * displacement:
 *
 *   deltaPos = halfStepSize * vel
 *   vel = velScale * vel + noiseScale * sqrt(vel.w) * randomValues
 *   return deltaPos + halfStepSize * vel
 *
 * where halfStepSize = 0.5*stepSize, velScale = exp(-stepSize*friction), and
 * noiseScale = sqrt(kT*(1 - velScale*velScale)).
 */

inline DEVICE mixed4 integrateMiddleStep(PRIVATE mixed4* velPtr, float4 randomValues, mixed halfStepSize, mixed velScale, mixed noiseScale) {
    mixed4 vel = *velPtr;
    mixed4 deltaPos = make_mixed4(halfStepSize*vel.x, halfStepSize*vel.y, halfStepSize*vel.z, 0);
    mixed sqrtInvMass = SQRT(vel.w);
    vel.x = velScale*vel.x + noiseScale*sqrtInvMass*randomValues.x;
    vel.y = velScale*vel.y + noiseScale*sqrtInvMass*randomValues.y;
    vel.z = velScale*vel.z + noiseScale*sqrtInvMass*randomValues.z;
    *velPtr = vel;
    return deltaPos + make_mixed4(halfStepSize*vel.x, halfStepSize*vel.y, halfStepSize*vel.z, 0);
}

/**
 * Perform the third part of integration for a constrained atom: apply the constraint displacement
 * correction to the velocity and commit the updated velocity and position:
 *
 *   vel += (deltaPos - unconstrainedDelta) * invStepSize
 *   pos += deltaPos
 *
 * where invStepSize = 1/stepSize, storing vel to velm[atom] and pos to posq (plus posqCorrection
 * when USE_MIXED_PRECISION is defined).
 */

inline DEVICE void commitConstrainedAtom(int atom, mixed4 vel, mixed4 pos, mixed4 deltaPos, mixed4 unconstrainedDelta, mixed invStepSize,
        GLOBAL mixed4* RESTRICT velm, GLOBAL real4* RESTRICT posq, GLOBAL real4* RESTRICT posqCorrection) {
    vel.x += (deltaPos.x-unconstrainedDelta.x)*invStepSize;
    vel.y += (deltaPos.y-unconstrainedDelta.y)*invStepSize;
    vel.z += (deltaPos.z-unconstrainedDelta.z)*invStepSize;
    velm[atom] = vel;
    pos.x += deltaPos.x;
    pos.y += deltaPos.y;
    pos.z += deltaPos.z;
    storePos(posq, posqCorrection, atom, pos);
}

/**
 * Perform a complete LangevinMiddle step in a single kernel when all constraints belong to
 * disjoint SETTLE 3-atom clusters and SHAKE central-atom clusters (no CCMA constraints):
 *
 *   1. Constrained clusters (SETTLE for index < numSettle, SHAKE otherwise):
 *      integrateForceStep -> velocity constraint projection -> integrateMiddleStep ->
 *      position constraint projection -> commitConstrainedAtom.
 *   2. Unconstrained atoms (vel.w != 0):
 *      integrateForceStep -> integrateMiddleStep -> pos += deltaPos.
 */

KERNEL void integrateLangevinMiddleFused(int numUnconstrainedAtoms, int paddedNumAtoms, int numSettle, int numShake, mixed tol,
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
    mixed stepSize = dt[0].y;
    mixed forceScale = stepSize/(mixed) 0x100000000;
    mixed velScale = paramBuffer[VelScale];
    mixed noiseScale = paramBuffer[NoiseScale];
    mixed halfStepSize = 0.5f*stepSize;
    mixed invStepSize = 1/stepSize;
    int numClusters = numSettle + numShake;

    // Integrate SETTLE and SHAKE constraint clusters.

    for (int index = GLOBAL_ID; index < numClusters; index += GLOBAL_SIZE) {
        if (index < numSettle) {
            int4 atoms = settleAtoms[index];
            float2 params = settleParams[index];
            mixed4 pos0 = loadPos(posq, posqCorrection, atoms.x);
            mixed4 pos1 = loadPos(posq, posqCorrection, atoms.y);
            mixed4 pos2 = loadPos(posq, posqCorrection, atoms.z);
            mixed4 vel0 = integrateForceStep(velm[atoms.x], atoms.x, paddedNumAtoms, forceScale, force);
            mixed4 vel1 = integrateForceStep(velm[atoms.y], atoms.y, paddedNumAtoms, forceScale, force);
            mixed4 vel2 = integrateForceStep(velm[atoms.z], atoms.z, paddedNumAtoms, forceScale, force);
            applySettleToVelocitiesCluster(pos0, pos1, pos2, &vel0, &vel1, &vel2);
            mixed4 deltaPos0 = integrateMiddleStep(&vel0, random[startRandomIndex+atoms.x], halfStepSize, velScale, noiseScale);
            mixed4 deltaPos1 = integrateMiddleStep(&vel1, random[startRandomIndex+atoms.y], halfStepSize, velScale, noiseScale);
            mixed4 deltaPos2 = integrateMiddleStep(&vel2, random[startRandomIndex+atoms.z], halfStepSize, velScale, noiseScale);
            mixed4 unconstrainedDelta0 = deltaPos0, unconstrainedDelta1 = deltaPos1, unconstrainedDelta2 = deltaPos2;
            mixed mass0 = 1/vel0.w, mass1 = 1/vel1.w, mass2 = 1/vel2.w;
            applySettleToPositionsCluster(params, pos0, pos1, pos2, mass0, mass1, mass2, &deltaPos0, &deltaPos1, &deltaPos2);
            commitConstrainedAtom(atoms.x, vel0, pos0, deltaPos0, unconstrainedDelta0, invStepSize, velm, posq, posqCorrection);
            commitConstrainedAtom(atoms.y, vel1, pos1, deltaPos1, unconstrainedDelta1, invStepSize, velm, posq, posqCorrection);
            commitConstrainedAtom(atoms.z, vel2, pos2, deltaPos2, unconstrainedDelta2, invStepSize, velm, posq, posqCorrection);
        }
#ifdef HAS_SHAKE_CLUSTERS
        else {
            int shakeIndex = index - numSettle;
            int4 atoms = shakeAtoms[shakeIndex];
            float4 params = shakeParams[shakeIndex];
            mixed4 pos0 = loadPos(posq, posqCorrection, atoms.x);
            mixed4 pos1 = loadPos(posq, posqCorrection, atoms.y);
            mixed4 pos2 = (atoms.z != -1 ? loadPos(posq, posqCorrection, atoms.z) : make_mixed4(0));
            mixed4 pos3 = (atoms.w != -1 ? loadPos(posq, posqCorrection, atoms.w) : make_mixed4(0));
            mixed4 vel0 = integrateForceStep(velm[atoms.x], atoms.x, paddedNumAtoms, forceScale, force);
            mixed4 vel1 = integrateForceStep(velm[atoms.y], atoms.y, paddedNumAtoms, forceScale, force);
            mixed4 vel2 = (atoms.z != -1 ? integrateForceStep(velm[atoms.z], atoms.z, paddedNumAtoms, forceScale, force) : make_mixed4(0));
            mixed4 vel3 = (atoms.w != -1 ? integrateForceStep(velm[atoms.w], atoms.w, paddedNumAtoms, forceScale, force) : make_mixed4(0));
            applyShakeToVelocitiesCluster(atoms, params, tol, pos0, pos1, pos2, pos3, &vel0, &vel1, &vel2, &vel3);
            mixed4 deltaPos0 = integrateMiddleStep(&vel0, random[startRandomIndex+atoms.x], halfStepSize, velScale, noiseScale);
            mixed4 deltaPos1 = integrateMiddleStep(&vel1, random[startRandomIndex+atoms.y], halfStepSize, velScale, noiseScale);
            mixed4 deltaPos2 = (atoms.z != -1 ? integrateMiddleStep(&vel2, random[startRandomIndex+atoms.z], halfStepSize, velScale, noiseScale) : make_mixed4(0));
            mixed4 deltaPos3 = (atoms.w != -1 ? integrateMiddleStep(&vel3, random[startRandomIndex+atoms.w], halfStepSize, velScale, noiseScale) : make_mixed4(0));
            mixed4 unconstrainedDelta0 = deltaPos0, unconstrainedDelta1 = deltaPos1, unconstrainedDelta2 = deltaPos2, unconstrainedDelta3 = deltaPos3;
            applyShakeToPositionsCluster(atoms, params, tol, pos0, pos1, pos2, pos3, &deltaPos0, &deltaPos1, &deltaPos2, &deltaPos3);
            commitConstrainedAtom(atoms.x, vel0, pos0, deltaPos0, unconstrainedDelta0, invStepSize, velm, posq, posqCorrection);
            commitConstrainedAtom(atoms.y, vel1, pos1, deltaPos1, unconstrainedDelta1, invStepSize, velm, posq, posqCorrection);
            if (atoms.z != -1)
                commitConstrainedAtom(atoms.z, vel2, pos2, deltaPos2, unconstrainedDelta2, invStepSize, velm, posq, posqCorrection);
            if (atoms.w != -1)
                commitConstrainedAtom(atoms.w, vel3, pos3, deltaPos3, unconstrainedDelta3, invStepSize, velm, posq, posqCorrection);
        }
#endif
    }

    // Integrate unconstrained atoms.

    for (int index = GLOBAL_ID; index < numUnconstrainedAtoms; index += GLOBAL_SIZE) {
        int atom = (numClusters == 0 ? index : unconstrainedAtoms[index]);
        mixed4 vel = velm[atom];
        if (vel.w != 0.0) {
            vel = integrateForceStep(vel, atom, paddedNumAtoms, forceScale, force);
            mixed4 deltaPos = integrateMiddleStep(&vel, random[startRandomIndex+atom], halfStepSize, velScale, noiseScale);
            velm[atom] = vel;
            mixed4 pos = loadPos(posq, posqCorrection, atom);
            pos.x += deltaPos.x;
            pos.y += deltaPos.y;
            pos.z += deltaPos.z;
            storePos(posq, posqCorrection, atom, pos);
        }
    }
}

#endif
