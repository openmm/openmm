typedef struct {
    mm_ulong state, increment;
    float next;
    bool nextIsValid;
} RandomState;

/**
 * Generate a random number with PCH-XSH-RR.
 */
DEVICE unsigned int getRandomInt(RandomState* random) {
    unsigned int xs = ((random->state>>18)^random->state)>>27;
    unsigned int rot = random->state>>59;
    random->state = random->state*6364136223846793005ULL + random->increment;
    return (xs>>rot) | (xs<<((-rot)&31));
}

/**
 * Generate a normally distributed random number with Box-Muller.
 */
DEVICE float getRandomNormal(RandomState* random) {
    if (random->nextIsValid) {
        random->nextIsValid = false;
        return random->next;
    }
    float scale = 1/(float) 0x100000000;
    float x = scale*max(getRandomInt(random), 1u);
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
        GLOBAL const mixed2* RESTRICT dt, GLOBAL int* RESTRICT tileCounter) {
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
    if (GLOBAL_ID == 0)
        *tileCounter = 0;
}

/**
 * Load the position of a particle.
 */
inline DEVICE mixed3 loadPos(GLOBAL const real4* RESTRICT posq, GLOBAL const real4* RESTRICT posqCorrection, int index) {
#ifdef USE_MIXED_PRECISION
    real4 pos1 = posq[index];
    real4 pos2 = posqCorrection[index];
    return make_mixed3(pos1.x+(mixed)pos2.x, pos1.y+(mixed)pos2.y, pos1.z+(mixed)pos2.z);
#else
    return trimTo3(posq[index]);
#endif
}
/**
 * Compute the friction and noise for a pair of particles.
 */
DEVICE void processPair(int i, int j, real3 delta, mixed4 vel1, mixed4 vel2, int type1, int type2, int paddedNumAtoms,
        GLOBAL mm_ulong* RESTRICT velDelta, mixed dt, float kT, int numTypes, GLOBAL const float2* RESTRICT params, RandomState* random,
        real4 periodicBoxSize, real4 invPeriodicBoxSize, real4 periodicBoxVecX, real4 periodicBoxVecY, real4 periodicBoxVecZ) {
    if (vel1.w == 0 && vel2.w == 0)
        return;
    float2 pairParams = params[type1+type2*numTypes];
    float friction = pairParams.x;
    float cutoff = pairParams.y;
    real r2 = delta.x*delta.x + delta.y*delta.y + delta.z*delta.z;
    real invR = RSQRT(r2);
    real r = r2*invR;
    if (r >= cutoff)
        return;
    real mass1 = (vel1.w == 0 ? 0 : 1/vel1.w);
    real mass2 = (vel2.w == 0 ? 0 : 1/vel2.w);
    real m = mass1*mass2/(mass1+mass2);
    real omega = 1-(r/cutoff);
    real vscale = EXP(-dt*2*friction*omega*omega);
    real noisescale = SQRT(1-vscale*vscale);
    real3 dir = delta*invR;
    real3 v = make_real3(vel2.x-vel1.x, vel2.y-vel1.y, vel2.z-vel1.z);
    real dv = (1-vscale)*dot(v, dir) + noisescale*SQRT(kT/m)*getRandomNormal(random);
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
        float kT, real4 periodicBoxSize, real4 invPeriodicBoxSize, real4 periodicBoxVecX, real4 periodicBoxVecY, real4 periodicBoxVecZ,
        GLOBAL const int2* RESTRICT exclusionTiles, int numExclusionTiles, GLOBAL const int* RESTRICT tiles, GLOBAL const unsigned int* RESTRICT interactionCount,
        GLOBAL const real4* RESTRICT blockCenter, GLOBAL const real4* RESTRICT blockSize, GLOBAL const int* RESTRICT interactingAtoms, GLOBAL int* RESTRICT tileCounter
#ifdef USE_MIXED_PRECISION
        , GLOBAL const real4* RESTRICT posqCorrection
#endif
        ) {
    const int totalWarps = GLOBAL_SIZE/TILE_SIZE;
    const int warp = GLOBAL_ID/TILE_SIZE;
    const int tgx = LOCAL_ID & (TILE_SIZE-1);
    const int tbx = LOCAL_ID - tgx;
    mixed halfdt = 0.5f*dt[0].y;
    LOCAL mixed3 localPos[WORK_GROUP_SIZE];
    LOCAL mixed4 localVel[WORK_GROUP_SIZE];
    LOCAL volatile int localType[WORK_GROUP_SIZE];
#ifndef USE_MIXED_PRECISION
    GLOBAL real4* posqCorrection = 0;
#endif

    // Initialize the random generator for this thread.  The seed is incremented each step,
    // and the stream ID is the global thread index.  Skipping a variable number of values
    // also seems to be necessary to decorrelate the streams.

    RandomState random;
    random.state = 0;
    random.increment = (GLOBAL_ID<<1) | 1;
    random.nextIsValid = false;
    getRandomInt(&random);
    random.state += seed;
    getRandomInt(&random);
    for (int i = 0; i < LOCAL_ID%16; i++)
        getRandomInt(&random);

    // First loop: process fixed tiles (ones that contain exclusions).

    for (int tile = warp; tile < numExclusionTiles; tile += totalWarps) {
        const int2 tileIndices = exclusionTiles[tile];
        const unsigned int x = tileIndices.x;
        const unsigned int y = tileIndices.y;
        int atom1 = x*TILE_SIZE + tgx;
        mixed4 vel1 = velm[atom1];
        mixed3 pos1 = loadPos(posq, posqCorrection, atom1);
        if (vel1.w != 0)
            pos1 += halfdt*trimTo3(vel1);
        int type1 = particleType[atom1];
        if (x == y) {
            localVel[LOCAL_ID] = vel1;
            localPos[LOCAL_ID] = pos1;
            localType[LOCAL_ID] = type1;
        }
        else {
            int atom2 = y*TILE_SIZE + tgx;
            mixed4 vel2 = velm[atom2];
            mixed3 pos2 = loadPos(posq, posqCorrection, atom2);
            if (vel2.w != 0)
                pos2 += halfdt*trimTo3(vel2);
            localVel[LOCAL_ID] = vel2;
            localPos[LOCAL_ID] = pos2;
            localType[LOCAL_ID] = particleType[atom2];
        }
        SYNC_WARPS;
        if (atom1 < numAtoms) {
            for (int i = 0; i < TILE_SIZE; i++) {
#ifdef INTEL_WORKAROUND
                // Workaround for bug in Intel's OpenCL for CPUs.
                MEM_FENCE;
#endif
                int atom2 = y*TILE_SIZE+i;
                if ((x != y || atom1 < atom2) && atom2 < numAtoms) {
                    mixed3 pos2 = localPos[tbx+i];
                    real3 delta = make_real3(pos2.x-pos1.x, pos2.y-pos1.y, pos2.z-pos1.z);
#ifdef USE_PERIODIC
                    APPLY_PERIODIC_TO_DELTA(delta)
#endif
                    processPair(atom1, atom2, delta, vel1, localVel[tbx+i], type1, localType[tbx+i],
                                paddedNumAtoms, (GLOBAL mm_ulong*) velDelta, dt[0].y, kT, numTypes, params, &random,
                                periodicBoxSize, invPeriodicBoxSize, periodicBoxVecX, periodicBoxVecY, periodicBoxVecZ);
                }
            }
        }
#ifdef NVIDIA_WORKAROUND
        SYNC_THREADS;
#else
        SYNC_WARPS;
#endif
    }

    // Second loop: process tiles from the neighbor list.

    unsigned int numTiles = interactionCount[0];
    LOCAL int atomIndices[WORK_GROUP_SIZE];
    LOCAL int nextTile[WORK_GROUP_SIZE/TILE_SIZE];
    for (int tile = warp; tile < numTiles; tile += totalWarps) {
        if (tgx == 0)
            nextTile[tbx/TILE_SIZE] = ATOMIC_ADD(tileCounter, 1);
        SYNC_WARPS;
        int tileIndex = nextTile[tbx/TILE_SIZE];
        int x = tiles[tileIndex];
        real4 blockSizeX = blockSize[x];
        bool singlePeriodicCopy = (0.5f*periodicBoxSize.x-blockSizeX.x >= MAX_CUTOFF &&
                                   0.5f*periodicBoxSize.y-blockSizeX.y >= MAX_CUTOFF &&
                                   0.5f*periodicBoxSize.z-blockSizeX.z >= MAX_CUTOFF);
        int atom1 = x*TILE_SIZE + tgx;
        mixed4 vel1 = velm[atom1];
        mixed3 pos1 = loadPos(posq, posqCorrection, atom1);
        if (vel1.w != 0)
            pos1 += halfdt*trimTo3(vel1);
        int type1 = particleType[atom1];
        int atom2 = interactingAtoms[tileIndex*TILE_SIZE+tgx];
        atomIndices[LOCAL_ID] = atom2;
        mixed4 vel2 = velm[atom2];
        mixed3 pos2 = loadPos(posq, posqCorrection, atom2);
        if (vel2.w != 0)
            pos2 += halfdt*trimTo3(vel2);
        localVel[LOCAL_ID] = vel2;
        localPos[LOCAL_ID] = pos2;
        localType[LOCAL_ID] = particleType[atom2];
#ifdef USE_PERIODIC
        if (singlePeriodicCopy) {
            // The box is small enough that we can just translate all the atoms into a single periodic
            // box, then skip having to apply periodic boundary conditions later.

            real4 blockCenterX = blockCenter[x];
            APPLY_PERIODIC_TO_POS_WITH_CENTER(pos1, blockCenterX)
            APPLY_PERIODIC_TO_POS_WITH_CENTER(localPos[LOCAL_ID], blockCenterX)
            SYNC_WARPS;
            if (atom1 < numAtoms) {
                for (int i = 0; i < TILE_SIZE; i++) {
#ifdef INTEL_WORKAROUND
                    // Workaround for bug in Intel's OpenCL for CPUs.
                    MEM_FENCE;
#endif
                    int atom2 = atomIndices[tbx+i];
                    if (atom2 < numAtoms) {
                        pos2 = localPos[tbx+i];
                        real3 delta = make_real3(pos2.x-pos1.x, pos2.y-pos1.y, pos2.z-pos1.z);
                        processPair(atom1, atom2, delta, vel1, localVel[tbx+i], type1, localType[tbx+i],
                                    paddedNumAtoms, (GLOBAL mm_ulong*) velDelta, dt[0].y, kT, numTypes, params, &random,
                                    periodicBoxSize, invPeriodicBoxSize, periodicBoxVecX, periodicBoxVecY, periodicBoxVecZ);
                    }
                }
            }
        }
        else
#endif
        {
            // We need to apply periodic boundary conditions separately for each interaction.

            SYNC_WARPS;
            if (atom1 < numAtoms) {
                for (int i = 0; i < TILE_SIZE; i++) {
#ifdef INTEL_WORKAROUND
                    // Workaround for bug in Intel's OpenCL for CPUs.
                    MEM_FENCE;
#endif
                    int atom2 = atomIndices[tbx+i];
                    if (atom2 < numAtoms) {
                        pos2 = localPos[tbx+i];
                        real3 delta = make_real3(pos2.x-pos1.x, pos2.y-pos1.y, pos2.z-pos1.z);
#ifdef USE_PERIODIC
                        APPLY_PERIODIC_TO_DELTA(delta)
#endif
                        processPair(atom1, atom2, delta, vel1, localVel[tbx+i], type1, localType[tbx+i],
                                    paddedNumAtoms, (GLOBAL mm_ulong*) velDelta, dt[0].y, kT, numTypes, params, &random,
                                    periodicBoxSize, invPeriodicBoxSize, periodicBoxVecX, periodicBoxVecY, periodicBoxVecZ);
                    }
                }
            }
        }
        SYNC_WARPS;
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
 * Perform the fourth part of integration: apply constraint forces to velocities, then record
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
