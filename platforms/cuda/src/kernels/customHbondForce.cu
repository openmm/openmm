/**
 * Convert a real4 to a real3 by removing its last element.
 */
inline __device__ real3 trim(real4 v) {
    return make_real3(v.x, v.y, v.z);
}

/**
 * This does nothing, and just exists to simplify the code generation.
 */
inline __device__ real3 trim(real3 v) {
    return v;
}

/**
 * Compute the difference between two vectors, optionally taking periodic boundary conditions into account
 * and setting the fourth component to the squared magnitude.
 */
inline __device__ real4 delta(real4 vec1, real4 vec2, real4 periodicBoxSize, real4 invPeriodicBoxSize, real4 periodicBoxVecX, real4 periodicBoxVecY, real4 periodicBoxVecZ) {
    real4 result = make_real4(vec1.x-vec2.x, vec1.y-vec2.y, vec1.z-vec2.z, 0.0f);
#ifdef USE_PERIODIC
    APPLY_PERIODIC_TO_DELTA(result)
#endif
    result.w = result.x*result.x + result.y*result.y + result.z*result.z;
    return result;
}

/**
 * Compute the angle between two vectors.  The w component of each vector should contain the squared magnitude.
 */
inline __device__ real computeAngle(real4 vec1, real4 vec2) {
    real dotProduct = vec1.x*vec2.x + vec1.y*vec2.y + vec1.z*vec2.z;
    real cosine = dotProduct*RSQRT(vec1.w*vec2.w);
    real angle;
    if (cosine > 0.99f || cosine < -0.99f) {
        // We're close to the singularity in acos(), so take the cross product and use asin() instead.

        real3 crossProduct = cross(vec1, vec2);
        real scale = vec1.w*vec2.w;
        angle = ASIN(SQRT(dot(crossProduct, crossProduct)/scale));
        if (cosine < 0.0f)
            angle = M_PI-angle;
    }
    else
       angle = ACOS(cosine);
    return angle;
}

/**
 * Compute the cross product of two vectors, setting the fourth component to the squared magnitude.
 */
inline __device__ real4 computeCross(real4 vec1, real4 vec2) {
    real3 result = cross(vec1, vec2);
    return make_real4(result.x, result.y, result.z, result.x*result.x + result.y*result.y + result.z*result.z);
}

/**
 * Compute forces on donors.
 */
extern "C" __global__ void computeDonorForces(unsigned long long* __restrict__ force, mixed* __restrict__ energyBuffer, const real4* __restrict__ posq,
        const int4* __restrict__ exclusions, const int4* __restrict__ donorAtoms, const int4* __restrict__ acceptorAtoms, real4 periodicBoxSize, real4 invPeriodicBoxSize,
        real4 periodicBoxVecX, real4 periodicBoxVecY, real4 periodicBoxVecZ
        PARAMETER_ARGUMENTS) {
    extern __shared__ real4 posBuffer[];
    mixed energy = 0;
    real3 f1 = make_real3(0);
    real3 f2 = make_real3(0);
    real3 f3 = make_real3(0);
    for (int donorStart = 0; donorStart < NUM_DONORS; donorStart += blockDim.x*gridDim.x) {
        // Load information about the donor this thread will compute forces on.

        int donorIndex = donorStart+blockIdx.x*blockDim.x+threadIdx.x;
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
        for (int acceptorStart = 0; acceptorStart < NUM_ACCEPTORS; acceptorStart += blockDim.x) {
            // Load the next block of acceptors into local memory.

            __syncthreads();
            int blockSize = min((int) blockDim.x, NUM_ACCEPTORS-acceptorStart);
            if (threadIdx.x < blockSize) {
                int4 atoms2 = acceptorAtoms[acceptorStart+threadIdx.x];
                posBuffer[3*threadIdx.x] = (atoms2.x > -1 ? posq[atoms2.x] : make_real4(0));
                posBuffer[3*threadIdx.x+1] = (atoms2.y > -1 ? posq[atoms2.y] : make_real4(0));
                posBuffer[3*threadIdx.x+2] = (atoms2.z > -1 ? posq[atoms2.z] : make_real4(0));
            }
            __syncthreads();
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
            if (atoms.x > -1) {
                atomicAdd(&force[atoms.x], static_cast<unsigned long long>((long long) (f1.x*0x100000000)));
                atomicAdd(&force[atoms.x+PADDED_NUM_ATOMS], static_cast<unsigned long long>((long long) (f1.y*0x100000000)));
                atomicAdd(&force[atoms.x+2*PADDED_NUM_ATOMS], static_cast<unsigned long long>((long long) (f1.z*0x100000000)));
                __threadfence_block();
            }
            if (atoms.y > -1) {
                atomicAdd(&force[atoms.y], static_cast<unsigned long long>((long long) (f2.x*0x100000000)));
                atomicAdd(&force[atoms.y+PADDED_NUM_ATOMS], static_cast<unsigned long long>((long long) (f2.y*0x100000000)));
                atomicAdd(&force[atoms.y+2*PADDED_NUM_ATOMS], static_cast<unsigned long long>((long long) (f2.z*0x100000000)));
                __threadfence_block();
            }
            if (atoms.z > -1) {
                atomicAdd(&force[atoms.z], static_cast<unsigned long long>((long long) (f3.x*0x100000000)));
                atomicAdd(&force[atoms.z+PADDED_NUM_ATOMS], static_cast<unsigned long long>((long long) (f3.y*0x100000000)));
                atomicAdd(&force[atoms.z+2*PADDED_NUM_ATOMS], static_cast<unsigned long long>((long long) (f3.z*0x100000000)));
                __threadfence_block();
            }
        }
    }
    energyBuffer[blockIdx.x*blockDim.x+threadIdx.x] += energy;
}
/**
 * Compute forces on acceptors.
 */
extern "C" __global__ void computeAcceptorForces(unsigned long long* __restrict__ force, mixed* __restrict__ energyBuffer, const real4* __restrict__ posq,
        const int4* __restrict__ exclusions, const int4* __restrict__ donorAtoms, const int4* __restrict__ acceptorAtoms, real4 periodicBoxSize, real4 invPeriodicBoxSize,
        real4 periodicBoxVecX, real4 periodicBoxVecY, real4 periodicBoxVecZ
        PARAMETER_ARGUMENTS) {
    extern __shared__ real4 posBuffer[];
    real3 f1 = make_real3(0);
    real3 f2 = make_real3(0);
    real3 f3 = make_real3(0);
    for (int acceptorStart = 0; acceptorStart < NUM_ACCEPTORS; acceptorStart += blockDim.x*gridDim.x) {
        // Load information about the acceptor this thread will compute forces on.

        int acceptorIndex = acceptorStart+blockIdx.x*blockDim.x+threadIdx.x;
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
        for (int donorStart = 0; donorStart < NUM_DONORS; donorStart += blockDim.x) {
            // Load the next block of donors into local memory.

            __syncthreads();
            int blockSize = min((int) blockDim.x, NUM_DONORS-donorStart);
            if (threadIdx.x < blockSize) {
                int4 atoms2 = donorAtoms[donorStart+threadIdx.x];
                posBuffer[3*threadIdx.x] = (atoms2.x > -1 ? posq[atoms2.x] : make_real4(0));
                posBuffer[3*threadIdx.x+1] = (atoms2.y > -1 ? posq[atoms2.y] : make_real4(0));
                posBuffer[3*threadIdx.x+2] = (atoms2.z > -1 ? posq[atoms2.z] : make_real4(0));
            }
            __syncthreads();
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
            if (atoms.x > -1) {
                atomicAdd(&force[atoms.x], static_cast<unsigned long long>((long long) (f1.x*0x100000000)));
                atomicAdd(&force[atoms.x+PADDED_NUM_ATOMS], static_cast<unsigned long long>((long long) (f1.y*0x100000000)));
                atomicAdd(&force[atoms.x+2*PADDED_NUM_ATOMS], static_cast<unsigned long long>((long long) (f1.z*0x100000000)));
                __threadfence_block();
            }
            if (atoms.y > -1) {
                atomicAdd(&force[atoms.y], static_cast<unsigned long long>((long long) (f2.x*0x100000000)));
                atomicAdd(&force[atoms.y+PADDED_NUM_ATOMS], static_cast<unsigned long long>((long long) (f2.y*0x100000000)));
                atomicAdd(&force[atoms.y+2*PADDED_NUM_ATOMS], static_cast<unsigned long long>((long long) (f2.z*0x100000000)));
                __threadfence_block();
            }
            if (atoms.z > -1) {
                atomicAdd(&force[atoms.z], static_cast<unsigned long long>((long long) (f3.x*0x100000000)));
                atomicAdd(&force[atoms.z+PADDED_NUM_ATOMS], static_cast<unsigned long long>((long long) (f3.y*0x100000000)));
                atomicAdd(&force[atoms.z+2*PADDED_NUM_ATOMS], static_cast<unsigned long long>((long long) (f3.z*0x100000000)));
                __threadfence_block();
            }
        }
    }
}
