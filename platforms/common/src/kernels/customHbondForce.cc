/**
 * Compute the difference between two vectors, optionally taking periodic boundary conditions into account
 * and setting the fourth component to the squared magnitude.
 */
inline DEVICE real4 delta(real4 vec1, real4 vec2, real4 periodicBoxSize, real4 invPeriodicBoxSize, real4 periodicBoxVecX, real4 periodicBoxVecY, real4 periodicBoxVecZ) {
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
 * Compute forces on donors.
 */
KERNEL void computeDonorForces(
#ifdef SUPPORTS_64_BIT_ATOMICS
	GLOBAL mm_ulong* RESTRICT force,
#else
	GLOBAL real4* RESTRICT forceBuffers, GLOBAL const int4* RESTRICT donorBufferIndices,
#endif
	GLOBAL mixed* RESTRICT energyBuffer, GLOBAL const real4* RESTRICT posq, GLOBAL const int4* RESTRICT exclusions,
        GLOBAL const int4* RESTRICT donorAtoms, GLOBAL const int4* RESTRICT acceptorAtoms, real4 periodicBoxSize, real4 invPeriodicBoxSize,
        real4 periodicBoxVecX, real4 periodicBoxVecY, real4 periodicBoxVecZ
        PARAMETER_ARGUMENTS) {
    LOCAL real4 posBuffer[3*THREAD_BLOCK_SIZE];
    mixed energy = 0;
    real3 f1 = make_real3(0);
    real3 f2 = make_real3(0);
    real3 f3 = make_real3(0);
    for (int donorStart = 0; donorStart < NUM_DONORS; donorStart += GLOBAL_SIZE) {
        // Load information about the donor this thread will compute forces on.

        int donorIndex = donorStart+GLOBAL_ID;
        int4 atoms, exclusionIndices;
        real4 d1, d2, d3;
        if (donorIndex < NUM_DONORS) {
            atoms = donorAtoms[donorIndex];
            d1 = (atoms.x > -1 ? posq[atoms.x] : make_real4(0));
            d2 = (atoms.y > -1 ? posq[atoms.y] : make_real4(0));
            d3 = (atoms.z > -1 ? posq[atoms.z] : make_real4(0));
#ifdef USE_EXCLUSIONS
            exclusionIndices = exclusions[donorIndex];
#endif
        }
        else
            atoms = make_int4(-1, -1, -1, -1);
        for (int acceptorStart = 0; acceptorStart < NUM_ACCEPTORS; acceptorStart += LOCAL_SIZE) {
            // Load the next block of acceptors into local memory.

            SYNC_THREADS;
            int blockSize = min((int) LOCAL_SIZE, NUM_ACCEPTORS-acceptorStart);
            if (LOCAL_ID < blockSize) {
                int4 atoms2 = acceptorAtoms[acceptorStart+LOCAL_ID];
                posBuffer[3*LOCAL_ID] = (atoms2.x > -1 ? posq[atoms2.x] : make_real4(0));
                posBuffer[3*LOCAL_ID+1] = (atoms2.y > -1 ? posq[atoms2.y] : make_real4(0));
                posBuffer[3*LOCAL_ID+2] = (atoms2.z > -1 ? posq[atoms2.z] : make_real4(0));
            }
            SYNC_THREADS;
            if (donorIndex < NUM_DONORS) {
                for (int index = 0; index < blockSize; index++) {
                    int acceptorIndex = acceptorStart+index;
#ifdef USE_EXCLUSIONS
                    if (acceptorIndex == exclusionIndices.x || acceptorIndex == exclusionIndices.y || acceptorIndex == exclusionIndices.z || acceptorIndex == exclusionIndices.w)
                        continue;
#endif
                    // Compute the interaction between a donor and an acceptor.

                    real4 a1 = posBuffer[3*index];
                    real4 a2 = posBuffer[3*index+1];
                    real4 a3 = posBuffer[3*index+2];
                    real4 deltaD1A1 = delta(d1, a1, periodicBoxSize, invPeriodicBoxSize, periodicBoxVecX, periodicBoxVecY, periodicBoxVecZ);
#ifdef USE_CUTOFF
                    if (deltaD1A1.w < CUTOFF_SQUARED) {
#endif
                        COMPUTE_DONOR_FORCE
#ifdef USE_CUTOFF
                    }
#endif
                }
            }
        }

        // Write results

        if (donorIndex < NUM_DONORS) {
#ifdef SUPPORTS_64_BIT_ATOMICS
            if (atoms.x > -1) {
                ATOMIC_ADD(&force[atoms.x], (mm_ulong) ((mm_long) (f1.x*0x100000000)));
                ATOMIC_ADD(&force[atoms.x+PADDED_NUM_ATOMS], (mm_ulong) ((mm_long) (f1.y*0x100000000)));
                ATOMIC_ADD(&force[atoms.x+2*PADDED_NUM_ATOMS], (mm_ulong) ((mm_long) (f1.z*0x100000000)));
                MEM_FENCE;
            }
            if (atoms.y > -1) {
                ATOMIC_ADD(&force[atoms.y], (mm_ulong) ((mm_long) (f2.x*0x100000000)));
                ATOMIC_ADD(&force[atoms.y+PADDED_NUM_ATOMS], (mm_ulong) ((mm_long) (f2.y*0x100000000)));
                ATOMIC_ADD(&force[atoms.y+2*PADDED_NUM_ATOMS], (mm_ulong) ((mm_long) (f2.z*0x100000000)));
                MEM_FENCE;
            }
            if (atoms.z > -1) {
                ATOMIC_ADD(&force[atoms.z], (mm_ulong) ((mm_long) (f3.x*0x100000000)));
                ATOMIC_ADD(&force[atoms.z+PADDED_NUM_ATOMS], (mm_ulong) ((mm_long) (f3.y*0x100000000)));
                ATOMIC_ADD(&force[atoms.z+2*PADDED_NUM_ATOMS], (mm_ulong) ((mm_long) (f3.z*0x100000000)));
                MEM_FENCE;
            }
#else
            int4 bufferIndices = donorBufferIndices[donorIndex];
            if (atoms.x > -1) {
                unsigned int offset = atoms.x+bufferIndices.x*PADDED_NUM_ATOMS;
                real4 force = forceBuffers[offset];
                force.xyz += f1.xyz;
                forceBuffers[offset] = force;
            }
            if (atoms.y > -1) {
                unsigned int offset = atoms.y+bufferIndices.y*PADDED_NUM_ATOMS;
                real4 force = forceBuffers[offset];
                force.xyz += f2.xyz;
                forceBuffers[offset] = force;
            }
            if (atoms.z > -1) {
                unsigned int offset = atoms.z+bufferIndices.z*PADDED_NUM_ATOMS;
                real4 force = forceBuffers[offset];
                force.xyz += f3.xyz;
                forceBuffers[offset] = force;
            }
#endif
        }
    }
    energyBuffer[GLOBAL_ID] += energy;
}
/**
 * Compute forces on acceptors.
 */
KERNEL void computeAcceptorForces(
#ifdef SUPPORTS_64_BIT_ATOMICS
	GLOBAL mm_ulong* RESTRICT force,
#else
	GLOBAL real4* RESTRICT forceBuffers, GLOBAL const int4* RESTRICT acceptorBufferIndices,
#endif
        GLOBAL mixed* RESTRICT energyBuffer, GLOBAL const real4* RESTRICT posq, GLOBAL const int4* RESTRICT exclusions,
        GLOBAL const int4* RESTRICT donorAtoms, GLOBAL const int4* RESTRICT acceptorAtoms, real4 periodicBoxSize, real4 invPeriodicBoxSize,
        real4 periodicBoxVecX, real4 periodicBoxVecY, real4 periodicBoxVecZ
        PARAMETER_ARGUMENTS) {
    LOCAL real4 posBuffer[3*THREAD_BLOCK_SIZE];
    real3 f1 = make_real3(0);
    real3 f2 = make_real3(0);
    real3 f3 = make_real3(0);
    for (int acceptorStart = 0; acceptorStart < NUM_ACCEPTORS; acceptorStart += GLOBAL_SIZE) {
        // Load information about the acceptor this thread will compute forces on.

        int acceptorIndex = acceptorStart+GLOBAL_ID;
        int4 atoms, exclusionIndices;
        real4 a1, a2, a3;
        if (acceptorIndex < NUM_ACCEPTORS) {
            atoms = acceptorAtoms[acceptorIndex];
            a1 = (atoms.x > -1 ? posq[atoms.x] : make_real4(0));
            a2 = (atoms.y > -1 ? posq[atoms.y] : make_real4(0));
            a3 = (atoms.z > -1 ? posq[atoms.z] : make_real4(0));
#ifdef USE_EXCLUSIONS
            exclusionIndices = exclusions[acceptorIndex];
#endif
        }
        else
            atoms = make_int4(-1, -1, -1, -1);
        for (int donorStart = 0; donorStart < NUM_DONORS; donorStart += LOCAL_SIZE) {
            // Load the next block of donors into local memory.

            SYNC_THREADS;
            int blockSize = min((int) LOCAL_SIZE, NUM_DONORS-donorStart);
            if (LOCAL_ID < blockSize) {
                int4 atoms2 = donorAtoms[donorStart+LOCAL_ID];
                posBuffer[3*LOCAL_ID] = (atoms2.x > -1 ? posq[atoms2.x] : make_real4(0));
                posBuffer[3*LOCAL_ID+1] = (atoms2.y > -1 ? posq[atoms2.y] : make_real4(0));
                posBuffer[3*LOCAL_ID+2] = (atoms2.z > -1 ? posq[atoms2.z] : make_real4(0));
            }
            SYNC_THREADS;
            if (acceptorIndex < NUM_ACCEPTORS) {
                for (int index = 0; index < blockSize; index++) {
                    int donorIndex = donorStart+index;
#ifdef USE_EXCLUSIONS
                    if (donorIndex == exclusionIndices.x || donorIndex == exclusionIndices.y || donorIndex == exclusionIndices.z || donorIndex == exclusionIndices.w)
                        continue;
#endif
                    // Compute the interaction between a donor and an acceptor.

                    real4 d1 = posBuffer[3*index];
                    real4 d2 = posBuffer[3*index+1];
                    real4 d3 = posBuffer[3*index+2];
                    real4 deltaD1A1 = delta(d1, a1, periodicBoxSize, invPeriodicBoxSize, periodicBoxVecX, periodicBoxVecY, periodicBoxVecZ);
#ifdef USE_CUTOFF
                    if (deltaD1A1.w < CUTOFF_SQUARED) {
#endif
                        COMPUTE_ACCEPTOR_FORCE
#ifdef USE_CUTOFF
                    }
#endif
                }
            }
        }

        // Write results

        if (acceptorIndex < NUM_ACCEPTORS) {
#ifdef SUPPORTS_64_BIT_ATOMICS
            if (atoms.x > -1) {
                ATOMIC_ADD(&force[atoms.x], (mm_ulong) ((mm_long) (f1.x*0x100000000)));
                ATOMIC_ADD(&force[atoms.x+PADDED_NUM_ATOMS], (mm_ulong) ((mm_long) (f1.y*0x100000000)));
                ATOMIC_ADD(&force[atoms.x+2*PADDED_NUM_ATOMS], (mm_ulong) ((mm_long) (f1.z*0x100000000)));
                MEM_FENCE;
            }
            if (atoms.y > -1) {
                ATOMIC_ADD(&force[atoms.y], (mm_ulong) ((mm_long) (f2.x*0x100000000)));
                ATOMIC_ADD(&force[atoms.y+PADDED_NUM_ATOMS], (mm_ulong) ((mm_long) (f2.y*0x100000000)));
                ATOMIC_ADD(&force[atoms.y+2*PADDED_NUM_ATOMS], (mm_ulong) ((mm_long) (f2.z*0x100000000)));
                MEM_FENCE;
            }
            if (atoms.z > -1) {
                ATOMIC_ADD(&force[atoms.z], (mm_ulong) ((mm_long) (f3.x*0x100000000)));
                ATOMIC_ADD(&force[atoms.z+PADDED_NUM_ATOMS], (mm_ulong) ((mm_long) (f3.y*0x100000000)));
                ATOMIC_ADD(&force[atoms.z+2*PADDED_NUM_ATOMS], (mm_ulong) ((mm_long) (f3.z*0x100000000)));
                MEM_FENCE;
            }
#else
            int4 bufferIndices = acceptorBufferIndices[acceptorIndex];
            if (atoms.x > -1) {
                unsigned int offset = atoms.x+bufferIndices.x*PADDED_NUM_ATOMS;
                real4 force = forceBuffers[offset];
                force.xyz += f1.xyz;
                forceBuffers[offset] = force;
            }
            if (atoms.y > -1) {
                unsigned int offset = atoms.y+bufferIndices.y*PADDED_NUM_ATOMS;
                real4 force = forceBuffers[offset];
                force.xyz += f2.xyz;
                forceBuffers[offset] = force;
            }
            if (atoms.z > -1) {
                unsigned int offset = atoms.z+bufferIndices.z*PADDED_NUM_ATOMS;
                real4 force = forceBuffers[offset];
                force.xyz += f3.xyz;
                forceBuffers[offset] = force;
            }
#endif
        }
    }
}
