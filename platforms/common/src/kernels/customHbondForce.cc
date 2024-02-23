DEVICE void findBoundingBox(GLOBAL const real4* RESTRICT posq, GLOBAL const int4* RESTRICT atoms, int numGroups,
        real4 periodicBoxSize, real4 invPeriodicBoxSize, real4 periodicBoxVecX, real4 periodicBoxVecY, real4 periodicBoxVecZ,
        GLOBAL real4* center, GLOBAL real4* blockSize) {
    real4 pos = posq[atoms[0].x];
#ifdef USE_PERIODIC
    APPLY_PERIODIC_TO_POS(pos)
#endif
    real4 minPos = pos;
    real4 maxPos = pos;
    for (int i = 1; i < numGroups; i++) {
        pos = posq[atoms[i].x];
#ifdef USE_PERIODIC
        real4 center = 0.5f*(maxPos+minPos);
        APPLY_PERIODIC_TO_POS_WITH_CENTER(pos, center)
#endif
        minPos = make_real4(min(minPos.x,pos.x), min(minPos.y,pos.y), min(minPos.z,pos.z), 0);
        maxPos = make_real4(max(maxPos.x,pos.x), max(maxPos.y,pos.y), max(maxPos.z,pos.z), 0);
    }
    *blockSize = 0.5f*(maxPos-minPos);
    *center = 0.5f*(maxPos+minPos);
}

KERNEL void findBlockBounds(GLOBAL const int4* RESTRICT donorAtoms, GLOBAL const int4* RESTRICT acceptorAtoms,
        real4 periodicBoxSize, real4 invPeriodicBoxSize, real4 periodicBoxVecX, real4 periodicBoxVecY, real4 periodicBoxVecZ,
        GLOBAL const real4* RESTRICT posq, GLOBAL real4* RESTRICT donorBlockCenter, GLOBAL real4* RESTRICT donorBlockSize,
        GLOBAL real4* RESTRICT acceptorBlockCenter, GLOBAL real4* RESTRICT acceptorBlockSize) {
    for (int index = GLOBAL_ID; index < NUM_DONOR_BLOCKS; index += GLOBAL_SIZE) {
        findBoundingBox(posq, donorAtoms+index*32, min(32, NUM_DONORS-index*32), periodicBoxSize,
                invPeriodicBoxSize, periodicBoxVecX, periodicBoxVecY, periodicBoxVecZ, donorBlockCenter+index,
                donorBlockSize+index);
    }
    for (int index = GLOBAL_ID; index < NUM_ACCEPTOR_BLOCKS; index += GLOBAL_SIZE) {
        findBoundingBox(posq, acceptorAtoms+index*32, min(32, NUM_ACCEPTORS-index*32), periodicBoxSize,
                invPeriodicBoxSize, periodicBoxVecX, periodicBoxVecY, periodicBoxVecZ, acceptorBlockCenter+index,
                acceptorBlockSize+index);
    }
}

/**
 * Compute the difference between two vectors, optionally taking periodic boundary conditions into account
 * and setting the fourth component to the squared magnitude.
 */
inline DEVICE real4 delta(real3 vec1, real3 vec2, real4 periodicBoxSize, real4 invPeriodicBoxSize, real4 periodicBoxVecX, real4 periodicBoxVecY, real4 periodicBoxVecZ) {
    real4 result = make_real4(vec1.x-vec2.x, vec1.y-vec2.y, vec1.z-vec2.z, 0);
#ifdef USE_PERIODIC
    APPLY_PERIODIC_TO_DELTA(result)
#endif
    result.w = result.x*result.x + result.y*result.y + result.z*result.z;
    return result;
}

/**
 * Compute the angle between two vectors.  The w component of each vector should contain the squared magnitude.
 */
inline DEVICE real computeAngle(real4 vec1, real4 vec2) {
    real dotProduct = vec1.x*vec2.x + vec1.y*vec2.y + vec1.z*vec2.z;
    real cosine = dotProduct*RSQRT(vec1.w*vec2.w);
    real angle;
    if (cosine > 0.99f || cosine < -0.99f) {
        // We're close to the singularity in acos(), so take the cross product and use asin() instead.

        real3 crossProduct = cross(trimTo3(vec1), trimTo3(vec2));
        real scale = vec1.w*vec2.w;
        angle = ASIN(SQRT(dot(crossProduct, crossProduct)/scale));
        if (cosine < 0)
            angle = M_PI-angle;
    }
    else
       angle = ACOS(cosine);
    return angle;
}

/**
 * Compute the cross product of two vectors, setting the fourth component to the squared magnitude.
 */
inline DEVICE real4 computeCross(real4 vec1, real4 vec2) {
    real3 cp = cross(trimTo3(vec1), trimTo3(vec2));
    return make_real4(cp.x, cp.y, cp.z, cp.x*cp.x+cp.y*cp.y+cp.z*cp.z);
}

/**
 * Write the force on an atom to global memory.
 */
inline DEVICE void applyForce(int atom, real3 f, GLOBAL mm_ulong* force) {
    if (atom > -1) {
        if (f.x != 0)
            ATOMIC_ADD(&force[atom], (mm_ulong) realToFixedPoint(f.x));
        if (f.y != 0)
            ATOMIC_ADD(&force[atom+PADDED_NUM_ATOMS], (mm_ulong) realToFixedPoint(f.y));
        if (f.z != 0)
            ATOMIC_ADD(&force[atom+2*PADDED_NUM_ATOMS], (mm_ulong) realToFixedPoint(f.z));
        MEM_FENCE;
    }
}

typedef struct {
    real3 pos1, pos2, pos3;
    real3 f1, f2, f3;
} AcceptorData;

/**
 * Compute forces on donors and acceptors.
 */
KERNEL void computeHbondForces(
	GLOBAL mm_ulong* RESTRICT force,
	GLOBAL mixed* RESTRICT energyBuffer, GLOBAL const real4* RESTRICT posq, GLOBAL const int4* RESTRICT exclusions,
        GLOBAL const int4* RESTRICT donorAtoms, GLOBAL const int4* RESTRICT acceptorAtoms, real4 periodicBoxSize, real4 invPeriodicBoxSize,
        real4 periodicBoxVecX, real4 periodicBoxVecY, real4 periodicBoxVecZ
#ifdef USE_BOUNDING_BOXES
        , GLOBAL real4* RESTRICT donorBlockCenter, GLOBAL real4* RESTRICT donorBlockSize,
        GLOBAL real4* RESTRICT acceptorBlockCenter, GLOBAL real4* RESTRICT acceptorBlockSize
#endif
        PARAMETER_ARGUMENTS) {
    const unsigned int totalWarps = GLOBAL_SIZE/32;
    const unsigned int warp = GLOBAL_ID/32;
    const int indexInWarp = GLOBAL_ID%32;
    const int tbx = LOCAL_ID-indexInWarp;
    LOCAL AcceptorData localData[THREAD_BLOCK_SIZE];
    mixed energy = 0;
    for (int tile = warp; tile < NUM_DONOR_BLOCKS*NUM_ACCEPTOR_BLOCKS; tile += totalWarps) {
        int donorBlock = tile/NUM_ACCEPTOR_BLOCKS;
        int acceptorBlock = tile%NUM_ACCEPTOR_BLOCKS;
#ifdef USE_BOUNDING_BOXES
        real4 blockDelta = donorBlockCenter[donorBlock]-acceptorBlockCenter[acceptorBlock];
#ifdef USE_PERIODIC
        APPLY_PERIODIC_TO_DELTA(blockDelta)
#endif
        real4 donorSize = donorBlockSize[donorBlock];
        real4 acceptorSize = acceptorBlockSize[acceptorBlock];
        blockDelta.x = max((real) 0, fabs(blockDelta.x)-donorSize.x-acceptorSize.x);
        blockDelta.y = max((real) 0, fabs(blockDelta.y)-donorSize.y-acceptorSize.y);
        blockDelta.z = max((real) 0, fabs(blockDelta.z)-donorSize.z-acceptorSize.z);
        if (blockDelta.x*blockDelta.x+blockDelta.y*blockDelta.y+blockDelta.z*blockDelta.z >= CUTOFF_SQUARED)
            continue;
#endif

        // Load information about the donor this thread will compute forces on.

        real3 f1 = make_real3(0);
        real3 f2 = make_real3(0);
        real3 f3 = make_real3(0);
        int donorIndex = donorBlock*32 + indexInWarp;
        int4 atoms, exclusionIndices;
        real3 d1, d2, d3;
        if (donorIndex < NUM_DONORS) {
            atoms = donorAtoms[donorIndex];
            d1 = (atoms.x > -1 ? trimTo3(posq[atoms.x]) : make_real3(0));
            d2 = (atoms.y > -1 ? trimTo3(posq[atoms.y]) : make_real3(0));
            d3 = (atoms.z > -1 ? trimTo3(posq[atoms.z]) : make_real3(0));
#ifdef USE_EXCLUSIONS
            exclusionIndices = exclusions[donorIndex];
#endif
        }
        else
            atoms = make_int4(-1, -1, -1, -1);

        // Load information about the acceptors into local memory.

        SYNC_WARPS;
        localData[LOCAL_ID].f1 = make_real3(0);
        localData[LOCAL_ID].f2 = make_real3(0);
        localData[LOCAL_ID].f3 = make_real3(0);
        int acceptorStart = acceptorBlock*32;
        int blockSize = min(32, NUM_ACCEPTORS-acceptorStart);
        int4 atoms2 = (indexInWarp < blockSize ? acceptorAtoms[acceptorStart+indexInWarp] : make_int4(-1));
        if (indexInWarp < blockSize) {
            localData[LOCAL_ID].pos1 = (atoms2.x > -1 ? trimTo3(posq[atoms2.x]) : make_real3(0));
            localData[LOCAL_ID].pos2 = (atoms2.y > -1 ? trimTo3(posq[atoms2.y]) : make_real3(0));
            localData[LOCAL_ID].pos3 = (atoms2.z > -1 ? trimTo3(posq[atoms2.z]) : make_real3(0));
        }
        SYNC_WARPS;
        if (donorIndex < NUM_DONORS) {
            int index = indexInWarp;
            for (int j = 0; j < 32; j++) {
                int acceptorIndex = acceptorStart+index;
#ifdef USE_EXCLUSIONS
                if (acceptorIndex < NUM_ACCEPTORS && acceptorIndex != exclusionIndices.x && acceptorIndex != exclusionIndices.y && acceptorIndex != exclusionIndices.z && acceptorIndex != exclusionIndices.w) {
#else
                if (acceptorIndex < NUM_ACCEPTORS) {
#endif
                    // Compute the interaction between a donor and an acceptor.

                    real3 a1 = localData[tbx+index].pos1;
                    real3 a2 = localData[tbx+index].pos2;
                    real3 a3 = localData[tbx+index].pos3;
                    real4 deltaD1A1 = delta(d1, a1, periodicBoxSize, invPeriodicBoxSize, periodicBoxVecX, periodicBoxVecY, periodicBoxVecZ);
#ifdef USE_CUTOFF
                    if (deltaD1A1.w < CUTOFF_SQUARED) {
#endif
                        COMPUTE_FORCE
#ifdef USE_CUTOFF
                    }
#endif
                }
                index = (index+1)%32;
            }
        }

        // Write results

        if (donorIndex < NUM_DONORS) {
            applyForce(atoms.x, f1, force);
            applyForce(atoms.y, f2, force);
            applyForce(atoms.z, f3, force);
        }
        SYNC_WARPS;
        applyForce(atoms2.x, localData[LOCAL_ID].f1, force);
        applyForce(atoms2.y, localData[LOCAL_ID].f2, force);
        applyForce(atoms2.z, localData[LOCAL_ID].f3, force);
    }
    energyBuffer[GLOBAL_ID] += energy;
}
