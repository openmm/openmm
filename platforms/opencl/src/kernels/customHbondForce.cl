/**
 * Compute the difference between two vectors, optionally taking periodic boundary conditions into account
 * and setting the fourth component to the squared magnitude.
 */
real4 delta(real4 vec1, real4 vec2, real4 periodicBoxSize, real4 invPeriodicBoxSize, real4 periodicBoxVecX, real4 periodicBoxVecY, real4 periodicBoxVecZ) {
    real4 result = (real4) (vec1.x-vec2.x, vec1.y-vec2.y, vec1.z-vec2.z, 0);
#ifdef USE_PERIODIC
    APPLY_PERIODIC_TO_DELTA(result)
#endif
    result.w = result.x*result.x + result.y*result.y + result.z*result.z;
    return result;
}

/**
 * Compute the angle between two vectors.  The w component of each vector should contain the squared magnitude.
 */
real computeAngle(real4 vec1, real4 vec2) {
    real dotProduct = vec1.x*vec2.x + vec1.y*vec2.y + vec1.z*vec2.z;
    real cosine = dotProduct*RSQRT(vec1.w*vec2.w);
    real angle;
    if (cosine > 0.99f || cosine < -0.99f) {
        // We're close to the singularity in acos(), so take the cross product and use asin() instead.

        real4 crossProduct = cross(vec1, vec2);
        real scale = vec1.w*vec2.w;
        angle = asin(SQRT(dot(crossProduct, crossProduct)/scale));
        if (cosine < 0.0f)
            angle = PI-angle;
    }
    else
       angle = acos(cosine);
    return angle;
}

/**
 * Compute the cross product of two vectors, setting the fourth component to the squared magnitude.
 */
real4 computeCross(real4 vec1, real4 vec2) {
    real4 result = cross(vec1, vec2);
    result.w = result.x*result.x + result.y*result.y + result.z*result.z;
    return result;
}

/**
 * Compute forces on donors.
 */
__kernel void computeDonorForces(__global real4* restrict forceBuffers, __global mixed* restrict energyBuffer, __global const real4* restrict posq, __global const int4* restrict exclusions,
        __global const int4* restrict donorAtoms, __global const int4* restrict acceptorAtoms, __global const int4* restrict donorBufferIndices, __local real4* posBuffer, real4 periodicBoxSize, real4 invPeriodicBoxSize,
        real4 periodicBoxVecX, real4 periodicBoxVecY, real4 periodicBoxVecZ
        PARAMETER_ARGUMENTS) {
    mixed energy = 0;
    real4 f1 = (real4) 0;
    real4 f2 = (real4) 0;
    real4 f3 = (real4) 0;
    for (int donorStart = 0; donorStart < NUM_DONORS; donorStart += get_global_size(0)) {
        // Load information about the donor this thread will compute forces on.

        int donorIndex = donorStart+get_global_id(0);
        int4 atoms, exclusionIndices;
        real4 d1, d2, d3;
        if (donorIndex < NUM_DONORS) {
            atoms = donorAtoms[donorIndex];
            d1 = (atoms.x > -1 ? posq[atoms.x] : (real4) 0);
            d2 = (atoms.y > -1 ? posq[atoms.y] : (real4) 0);
            d3 = (atoms.z > -1 ? posq[atoms.z] : (real4) 0);
#ifdef USE_EXCLUSIONS
            exclusionIndices = exclusions[donorIndex];
#endif
        }
        else
            atoms = (int4) (-1, -1, -1, -1);
        for (int acceptorStart = 0; acceptorStart < NUM_ACCEPTORS; acceptorStart += get_local_size(0)) {
            // Load the next block of acceptors into local memory.

            barrier(CLK_LOCAL_MEM_FENCE);
            int blockSize = min((int) get_local_size(0), NUM_ACCEPTORS-acceptorStart);
            if (get_local_id(0) < blockSize) {
                int4 atoms2 = acceptorAtoms[acceptorStart+get_local_id(0)];
                posBuffer[3*get_local_id(0)] = (atoms2.x > -1 ? posq[atoms2.x] : (real4) 0);
                posBuffer[3*get_local_id(0)+1] = (atoms2.y > -1 ? posq[atoms2.y] : (real4) 0);
                posBuffer[3*get_local_id(0)+2] = (atoms2.z > -1 ? posq[atoms2.z] : (real4) 0);
            }
            barrier(CLK_LOCAL_MEM_FENCE);
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
        }
    }
    energyBuffer[get_global_id(0)] += energy;
}
/**
 * Compute forces on acceptors.
 */
__kernel void computeAcceptorForces(__global real4* restrict forceBuffers, __global mixed* restrict energyBuffer, __global const real4* restrict posq, __global const int4* restrict exclusions,
        __global const int4* restrict donorAtoms, __global const int4* restrict acceptorAtoms, __global const int4* restrict acceptorBufferIndices, __local real4* restrict posBuffer, real4 periodicBoxSize, real4 invPeriodicBoxSize,
        real4 periodicBoxVecX, real4 periodicBoxVecY, real4 periodicBoxVecZ
        PARAMETER_ARGUMENTS) {
    real4 f1 = (real4) 0;
    real4 f2 = (real4) 0;
    real4 f3 = (real4) 0;
    for (int acceptorStart = 0; acceptorStart < NUM_ACCEPTORS; acceptorStart += get_global_size(0)) {
        // Load information about the acceptor this thread will compute forces on.

        int acceptorIndex = acceptorStart+get_global_id(0);
        int4 atoms, exclusionIndices;
        real4 a1, a2, a3;
        if (acceptorIndex < NUM_ACCEPTORS) {
            atoms = acceptorAtoms[acceptorIndex];
            a1 = (atoms.x > -1 ? posq[atoms.x] : (real4) 0);
            a2 = (atoms.y > -1 ? posq[atoms.y] : (real4) 0);
            a3 = (atoms.z > -1 ? posq[atoms.z] : (real4) 0);
#ifdef USE_EXCLUSIONS
            exclusionIndices = exclusions[acceptorIndex];
#endif
        }
        else
            atoms = (int4) (-1, -1, -1, -1);
        for (int donorStart = 0; donorStart < NUM_DONORS; donorStart += get_local_size(0)) {
            // Load the next block of donors into local memory.

            barrier(CLK_LOCAL_MEM_FENCE);
            int blockSize = min((int) get_local_size(0), NUM_DONORS-donorStart);
            if (get_local_id(0) < blockSize) {
                int4 atoms2 = donorAtoms[donorStart+get_local_id(0)];
                posBuffer[3*get_local_id(0)] = (atoms2.x > -1 ? posq[atoms2.x] : (real4) 0);
                posBuffer[3*get_local_id(0)+1] = (atoms2.y > -1 ? posq[atoms2.y] : (real4) 0);
                posBuffer[3*get_local_id(0)+2] = (atoms2.z > -1 ? posq[atoms2.z] : (real4) 0);
            }
            barrier(CLK_LOCAL_MEM_FENCE);
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
        }
    }
}
