typedef struct {
    mm_ulong state, increment;
    float next;
    bool nextIsValid;
} RandomState;

/**
 * Generate a random number with PCH-XSH-RR.
 */
unsigned int getRandomInt(RandomState* random) {
    unsigned int xs = ((random->state>>18)^random->state)>>27;
    unsigned int rot = random->state>>59;
    random->state = random->state*6364136223846793005ULL + random->increment;
    return (xs>>rot) | (xs<<((-rot)&31));
}

/**
 * Generate a normally distributed random number with Box-Muller.
 */
float getRandomNormal(RandomState* random) {
    if (random->nextIsValid) {
        random->nextIsValid = false;
        return random->next;
    }
    float scale = 1/(float) 0x100000000;
    float x = scale*getRandomInt(random);
    float y = scale*getRandomInt(random);
    float multiplier = SQRT(-2*LOG(x));
    float angle = 2*M_PI*y;
    random->next = multiplier*COS(angle);
    random->nextIsValid = true;
    return multiplier*SIN(angle);
}

/**
 * Perform the first part of integration: velocity step.
 */

KERNEL void integrateDPDPart1(int numAtoms, int paddedNumAtoms, GLOBAL mixed4* RESTRICT velm, GLOBAL const mm_long* RESTRICT force,
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
 * Load the position of a particle.
 */
inline DEVICE mixed4 loadPos(GLOBAL const real4* RESTRICT posq, GLOBAL const real4* RESTRICT posqCorrection, int index) {
#ifdef USE_MIXED_PRECISION
    real4 pos1 = posq[index];
    real4 pos2 = posqCorrection[index];
    return make_mixed4(pos1.x+(mixed)pos2.x, pos1.y+(mixed)pos2.y, pos1.z+(mixed)pos2.z, pos1.w);
#else
    return posq[index];
#endif
}
/**
 * Compute the friction and noise for a pair of particles.
 */
DEVICE void processPair(int i, int j, int paddedNumAtoms, GLOBAL real4* RESTRICT posq, GLOBAL mixed4* RESTRICT velm, GLOBAL mm_ulong* RESTRICT velDelta, mixed dt, float kT,
        GLOBAL const int* particleType, int numTypes, GLOBAL const float2* RESTRICT params, GLOBAL real4* RESTRICT posqCorrection, RandomState* random,
        real4 periodicBoxSize, real4 invPeriodicBoxSize, real4 periodicBoxVecX, real4 periodicBoxVecY, real4 periodicBoxVecZ) {
    mixed4 vel1 = velm[i];
    mixed4 vel2 = velm[j];
    if (vel1.w == 0 && vel2.w == 0)
        return;
    int type1 = particleType[i];
    int type2 = particleType[j];
    float2 pairParams = params[type1+type2*numTypes];
    float friction = pairParams.x;
    float cutoff = pairParams.y;
    mixed4 pos1 = loadPos(posq, posqCorrection, i);
    mixed4 pos2 = loadPos(posq, posqCorrection, j);
    mixed3 delta = make_mixed3(pos2.x-pos1.x, pos2.y-pos1.y, pos2.z-pos1.z);
#ifdef USE_PERIODIC
    APPLY_PERIODIC_TO_DELTA(delta)
#endif
    real r2 = delta.x*delta.x + delta.y*delta.y + delta.z*delta.z;
    real invR = RSQRT(r2);
    real r = r2*invR;
    if (r >= cutoff)
        return;
    mixed mass1 = (vel1.w == 0 ? 0 : 1/vel1.w);
    mixed mass2 = (vel2.w == 0 ? 0 : 1/vel2.w);
    mixed m = mass1*mass2/(mass1+mass2);
    mixed omega = 1-(r/cutoff);
    mixed vscale = EXP(-dt*friction*omega*omega);
    mixed noisescale = SQRT(1-vscale*vscale);
    mixed3 dir = delta*invR;
    mixed3 v = trimTo3(vel2)-trimTo3(vel1);
    mixed dv = (1-vscale)*dot(v, dir) + noisescale*SQRT(kT/m)*getRandomNormal(random);
    ATOMIC_ADD(&velDelta[i], (mm_ulong) realToFixedPoint(m*vel1.w*dv*dir.x));
    ATOMIC_ADD(&velDelta[i+paddedNumAtoms], (mm_ulong) realToFixedPoint(m*vel1.w*dv*dir.y));
    ATOMIC_ADD(&velDelta[i+2*paddedNumAtoms], (mm_ulong) realToFixedPoint(m*vel1.w*dv*dir.z));
    ATOMIC_ADD(&velDelta[j], (mm_ulong) realToFixedPoint(-m*vel2.w*dv*dir.x));
    ATOMIC_ADD(&velDelta[j+paddedNumAtoms], (mm_ulong) realToFixedPoint(-m*vel2.w*dv*dir.y));
    ATOMIC_ADD(&velDelta[j+2*paddedNumAtoms], (mm_ulong) realToFixedPoint(-m*vel2.w*dv*dir.z));
}

/**
 * Perform the second part of integration: position half step, then interact with heat bath.
 */
KERNEL void integrateDPDPart2(int numAtoms, int paddedNumAtoms, GLOBAL const real4* RESTRICT posq, GLOBAL const mixed4* RESTRICT velm, GLOBAL mm_long* RESTRICT velDelta,
        GLOBAL const mixed2* RESTRICT dt, GLOBAL const int* particleType, int numTypes, GLOBAL const float2* RESTRICT params, mm_long seed,
        float kT, real4 periodicBoxSize, real4 invPeriodicBoxSize, real4 periodicBoxVecX, real4 periodicBoxVecY, real4 periodicBoxVecZ
#ifdef USE_MIXED_PRECISION
        , GLOBAL const real4* RESTRICT posqCorrection
#endif
        ) {
#ifndef USE_MIXED_PRECISION
    GLOBAL real4* posqCorrection = 0;
#endif

    // Initialize the random generator for this thread.  The seed is incremented each step,
    // and the sequence ID is the global thread index.

    RandomState random;
    random.state = 0;
    random.increment = (GLOBAL_ID<<1) | 1;
    random.nextIsValid = false;
    getRandomInt(&random);
    random.state += seed;
    getRandomInt(&random);

    // Loop over atom pairs to compute the changes in velocity.

    int numPairs = numAtoms*(numAtoms+1)/2;
    for (int index = GLOBAL_ID; index < numPairs; index += GLOBAL_SIZE) {
        int j = (int) floor(numAtoms+0.5f-SQRT((numAtoms+0.5f)*(numAtoms+0.5f)-2*index));
        int i = (index-j*numAtoms+j*(j+1)/2);
        if (i < j || i >= numAtoms) { // Occasionally happens due to roundoff error.
            j += (i < j ? -1 : 1);
            i = (index-j*numAtoms+j*(j+1)/2);
        }
        if (i != j)
            processPair(i, j, paddedNumAtoms, posq, velm, (GLOBAL mm_ulong*) velDelta, dt[0].y, kT, particleType, numTypes, params, posqCorrection, &random,
                        periodicBoxSize, invPeriodicBoxSize, periodicBoxVecX, periodicBoxVecY, periodicBoxVecZ);
    }
}

/**
 * Perform the third part of integration: update the velocities and the second position half step.
 */
KERNEL void integrateDPDPart3(int numAtoms, int paddedNumAtoms, GLOBAL mixed4* RESTRICT velm, GLOBAL mm_long* RESTRICT velDelta, GLOBAL mixed4* RESTRICT posDelta,
        GLOBAL mixed4* RESTRICT oldDelta, GLOBAL const mixed2* RESTRICT dt) {
    mixed halfdt = 0.5f*dt[0].y;
    mixed deltascale = 1/(mixed) 0x100000000;
    for (int index = GLOBAL_ID; index < numAtoms; index += GLOBAL_SIZE) {
        mixed4 velocity = velm[index];
        if (velocity.w != 0.0) {
            mixed4 delta = make_mixed4(halfdt*velocity.x, halfdt*velocity.y, halfdt*velocity.z, 0);
            velocity.x += deltascale*velDelta[index];
            velocity.y += deltascale*velDelta[index+paddedNumAtoms];
            velocity.z += deltascale*velDelta[index+paddedNumAtoms*2];
            velm[index] = velocity;
            delta += make_mixed4(halfdt*velocity.x, halfdt*velocity.y, halfdt*velocity.z, 0);
            posDelta[index] = delta;
            oldDelta[index] = delta;
        }
    }
}

/**
 * Perform the third part of integration: apply constraint forces to velocities, then record
 * the constrained positions.
 */
KERNEL void integrateDPDPart4(int numAtoms, GLOBAL real4* RESTRICT posq, GLOBAL mixed4* RESTRICT velm,
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
