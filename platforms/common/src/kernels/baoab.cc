enum {VelScale, NoiseScale};

/**
 * Perform the first part of BAOAB integration: velocity half step, then position half step.
 */

KERNEL void integrateBAOABPart1(int numAtoms, int paddedNumAtoms, GLOBAL mixed4* RESTRICT velm, GLOBAL const mm_long* RESTRICT force, GLOBAL mixed4* RESTRICT posDelta,
        GLOBAL mixed4* RESTRICT oldDelta, GLOBAL const mixed2* RESTRICT dt) {
    mixed halfdt = 0.5*dt[0].y;
    mixed fscale = halfdt/(mixed) 0x100000000;
    for (int index = GLOBAL_ID; index < numAtoms; index += GLOBAL_SIZE) {
        mixed4 velocity = velm[index];
        if (velocity.w != 0.0) {
            velocity.x += fscale*velocity.w*force[index];
            velocity.y += fscale*velocity.w*force[index+paddedNumAtoms];
            velocity.z += fscale*velocity.w*force[index+paddedNumAtoms*2];
            velm[index] = velocity;
            mixed4 delta = make_mixed4(halfdt*velocity.x, halfdt*velocity.y, halfdt*velocity.z, 0);
            posDelta[index] = delta;
            oldDelta[index] = delta;
        }
    }
}

/**
 * Perform the second part of BAOAB integration: apply constraint forces to velocities, then interact with heat bath,
 * then position half step.
 */

KERNEL void integrateBAOABPart2(int numAtoms, GLOBAL real4* RESTRICT posq, GLOBAL mixed4* RESTRICT velm, GLOBAL mixed4* RESTRICT posDelta,
        GLOBAL mixed4* RESTRICT oldDelta, GLOBAL const mixed* RESTRICT paramBuffer, GLOBAL const mixed2* RESTRICT dt, GLOBAL const float4* RESTRICT random, unsigned int randomIndex
#ifdef USE_MIXED_PRECISION
        , GLOBAL real4* RESTRICT posqCorrection
#endif
        ) {
    mixed vscale = paramBuffer[VelScale];
    mixed noisescale = paramBuffer[NoiseScale];
    mixed halfdt = 0.5*dt[0].y;
    mixed invHalfdt = 1/halfdt;
    int index = GLOBAL_ID;
    randomIndex += index;
    while (index < numAtoms) {
        mixed4 velocity = velm[index];
        if (velocity.w != 0.0) {
            mixed4 delta = posDelta[index];
            mixed sqrtInvMass = SQRT(velocity.w);
            velocity.x += (delta.x-oldDelta[index].x)*invHalfdt;
            velocity.y += (delta.y-oldDelta[index].y)*invHalfdt;
            velocity.z += (delta.z-oldDelta[index].z)*invHalfdt;
            velocity.x = vscale*velocity.x + noisescale*sqrtInvMass*random[randomIndex].x;
            velocity.y = vscale*velocity.y + noisescale*sqrtInvMass*random[randomIndex].y;
            velocity.z = vscale*velocity.z + noisescale*sqrtInvMass*random[randomIndex].z;
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
            delta = make_mixed4(halfdt*velocity.x, halfdt*velocity.y, halfdt*velocity.z, 0);
            posDelta[index] = delta;
            oldDelta[index] = delta;
        }
        randomIndex += GLOBAL_SIZE;
        index += GLOBAL_SIZE;
    }
}

/**
 * Perform the third part of BAOAB integration: apply constraint forces to velocities, then record
 * the constrained positions in preparation for computing forces.
 */

KERNEL void integrateBAOABPart3(int numAtoms, GLOBAL real4* RESTRICT posq, GLOBAL mixed4* RESTRICT velm,
         GLOBAL mixed4* RESTRICT posDelta, GLOBAL mixed4* RESTRICT oldDelta, GLOBAL const mixed2* RESTRICT dt
#ifdef USE_MIXED_PRECISION
        , GLOBAL real4* RESTRICT posqCorrection
#endif
        ) {
    mixed halfdt = 0.5*dt[0].y;
    mixed invHalfdt = 1/halfdt;
    for (int index = GLOBAL_ID; index < numAtoms; index += GLOBAL_SIZE) {
        mixed4 velocity = velm[index];
        if (velocity.w != 0.0) {
            mixed4 delta = posDelta[index];
            velocity.x += (delta.x-oldDelta[index].x)*invHalfdt;
            velocity.y += (delta.y-oldDelta[index].y)*invHalfdt;
            velocity.z += (delta.z-oldDelta[index].z)*invHalfdt;
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
 * Perform the fourth part of BAOAB integration: velocity half step.
 */

KERNEL void integrateBAOABPart4(int numAtoms, int paddedNumAtoms, GLOBAL mixed4* RESTRICT velm,
        GLOBAL const mm_long* RESTRICT force, GLOBAL const mixed2* RESTRICT dt) {
    mixed halfdt = 0.5*dt[0].y;
    mixed fscale = halfdt/(mixed) 0x100000000;
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
