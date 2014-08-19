/**
 * Record the force on an atom to global memory.
 */
inline __device__ void storeForce(int atom, real3 force, unsigned long long* __restrict__ forceBuffers) {
    atomicAdd(&forceBuffers[atom], static_cast<unsigned long long>((long long) (force.x*0x100000000)));
    atomicAdd(&forceBuffers[atom+PADDED_NUM_ATOMS], static_cast<unsigned long long>((long long) (force.y*0x100000000)));
    atomicAdd(&forceBuffers[atom+2*PADDED_NUM_ATOMS], static_cast<unsigned long long>((long long) (force.z*0x100000000)));
}

/**
 * Convert a real4 to a real3 by removing its last element.
 */
inline __device__ real3 trim(real4 v) {
    return make_real3(v.x, v.y, v.z);
}

/**
 * Compute the difference between two vectors, taking periodic boundary conditions into account
 * and setting the fourth component to the squared magnitude.
 */
inline __device__ real4 delta(real4 vec1, real4 vec2, real4 periodicBoxSize, real4 invPeriodicBoxSize) {
    real4 result = make_real4(vec1.x-vec2.x, vec1.y-vec2.y, vec1.z-vec2.z, 0.0f);
#ifdef USE_PERIODIC
    result.x -= floor(result.x*invPeriodicBoxSize.x+0.5f)*periodicBoxSize.x;
    result.y -= floor(result.y*invPeriodicBoxSize.y+0.5f)*periodicBoxSize.y;
    result.z -= floor(result.z*invPeriodicBoxSize.z+0.5f)*periodicBoxSize.z;
#endif
    result.w = result.x*result.x + result.y*result.y + result.z*result.z;
    return result;
}

/**
 * Compute the angle between two vectors.  The w component of each vector should contain the squared magnitude.
 */
__device__ real computeAngle(real4 vec1, real4 vec2) {
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
    real3 cp = cross(vec1, vec2);
    return make_real4(cp.x, cp.y, cp.z, cp.x*cp.x+cp.y*cp.y+cp.z*cp.z);
}

/**
 * Compute the interaction.
 */
extern "C" __global__ void computeInteraction(
        unsigned long long* __restrict__ forceBuffers, real* __restrict__ energyBuffer, const real4* __restrict__ posq,
        real4 periodicBoxSize, real4 invPeriodicBoxSize
        PARAMETER_ARGUMENTS) {
    real energy = 0.0f;
    
    // Loop over particles to be the first one in the set.
    
    for (int p1 = blockIdx.x; p1 < NUM_ATOMS; p1 += gridDim.x) {
        int numNeighbors = NUM_ATOMS-p1-1;
        int numCombinations = NUM_CANDIDATE_COMBINATIONS;
        for (int index = threadIdx.x; index < numCombinations; index += blockDim.x) {
            FIND_ATOMS_FOR_COMBINATION_INDEX;
            bool includeInteraction = IS_VALID_COMBINATION;
#ifdef USE_CUTOFF
            if (includeInteraction) {
                VERIFY_CUTOFF;
            }
#endif
            if (includeInteraction) {
                PERMUTE_ATOMS;
                LOAD_PARTICLE_DATA;
                COMPUTE_INTERACTION;
            }
        }
    }
    energyBuffer[blockIdx.x*blockDim.x+threadIdx.x] += energy;
}