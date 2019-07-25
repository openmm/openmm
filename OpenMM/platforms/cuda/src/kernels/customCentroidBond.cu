/**
 * Compute the center of each group.
 */
extern "C" __global__ void computeGroupCenters(const real4* __restrict__ posq, const int* __restrict__ groupParticles,
        const real* __restrict__ groupWeights, const int* __restrict__ groupOffsets, real4* __restrict__ centerPositions) {
    __shared__ volatile real3 temp[64];
    for (int group = blockIdx.x; group < NUM_GROUPS; group += gridDim.x) {
        // The threads in this block work together to compute the center one group.
        
        int firstIndex = groupOffsets[group];
        int lastIndex = groupOffsets[group+1];
        real3 center = make_real3(0, 0, 0);
        for (int index = threadIdx.x; index < lastIndex-firstIndex; index += blockDim.x) {
            int atom = groupParticles[firstIndex+index];
            real weight = groupWeights[firstIndex+index];
            real4 pos = posq[atom];
            center.x += weight*pos.x;
            center.y += weight*pos.y;
            center.z += weight*pos.z;
        }
        
        // Sum the values.
        
        int thread = threadIdx.x;
        temp[thread].x = center.x;
        temp[thread].y = center.y;
        temp[thread].z = center.z;
        __syncthreads();
        if (thread < 32) {
            temp[thread].x += temp[thread+32].x;
            temp[thread].y += temp[thread+32].y;
            temp[thread].z += temp[thread+32].z;
            if (thread < 16) {
                temp[thread].x += temp[thread+16].x;
                temp[thread].y += temp[thread+16].y;
                temp[thread].z += temp[thread+16].z;
            }
            if (thread < 8) {
                temp[thread].x += temp[thread+8].x;
                temp[thread].y += temp[thread+8].y;
                temp[thread].z += temp[thread+8].z;
            }
            if (thread < 4) {
                temp[thread].x += temp[thread+4].x;
                temp[thread].y += temp[thread+4].y;
                temp[thread].z += temp[thread+4].z;
            }
            if (thread < 2) {
                temp[thread].x += temp[thread+2].x;
                temp[thread].y += temp[thread+2].y;
                temp[thread].z += temp[thread+2].z;
            }
        }
        if (thread == 0)
            centerPositions[group] = make_real4(temp[0].x+temp[1].x, temp[0].y+temp[1].y, temp[0].z+temp[1].z, 0);
    }
}

/**
 * Convert a real4 to a real3 by removing its last element.
 */
inline __device__ real3 trim(real4 v) {
    return make_real3(v.x, v.y, v.z);
}

/**
 * Compute the difference between two vectors, setting the fourth component to the squared magnitude.
 */
inline __device__ real4 delta(real4 vec1, real4 vec2, bool periodic, real4 periodicBoxSize, real4 invPeriodicBoxSize, 
        real4 periodicBoxVecX, real4 periodicBoxVecY, real4 periodicBoxVecZ) {
    real4 result = make_real4(vec1.x-vec2.x, vec1.y-vec2.y, vec1.z-vec2.z, 0);
    if (periodic)
        APPLY_PERIODIC_TO_DELTA(result);
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
 * Compute the forces on groups based on the bonds.
 */
extern "C" __global__ void computeGroupForces(unsigned long long* __restrict__ groupForce, mixed* __restrict__ energyBuffer, const real4* __restrict__ centerPositions,
        const int* __restrict__ bondGroups, real4 periodicBoxSize, real4 invPeriodicBoxSize, real4 periodicBoxVecX, real4 periodicBoxVecY, real4 periodicBoxVecZ
        EXTRA_ARGS) {
    mixed energy = 0;
    INIT_PARAM_DERIVS
    for (int index = blockIdx.x*blockDim.x+threadIdx.x; index < NUM_BONDS; index += blockDim.x*gridDim.x) {
        COMPUTE_FORCE
    }
    energyBuffer[blockIdx.x*blockDim.x+threadIdx.x] += energy;
    SAVE_PARAM_DERIVS
}

/**
 * Apply the forces from the group centers to the individual atoms.
 */
extern "C" __global__ void applyForcesToAtoms(const int* __restrict__ groupParticles, const real* __restrict__ groupWeights, const int* __restrict__ groupOffsets,
        const long long* __restrict__ groupForce, unsigned long long* __restrict__ atomForce) {
    for (int group = blockIdx.x; group < NUM_GROUPS; group += gridDim.x) {
        long long fx = groupForce[group];
        long long fy = groupForce[group+NUM_GROUPS];
        long long fz = groupForce[group+NUM_GROUPS*2];
        int firstIndex = groupOffsets[group];
        int lastIndex = groupOffsets[group+1];
        for (int index = threadIdx.x; index < lastIndex-firstIndex; index += blockDim.x) {
            int atom = groupParticles[firstIndex+index];
            real weight = groupWeights[firstIndex+index];
            atomicAdd(&atomForce[atom], static_cast<unsigned long long>((long long) (fx*weight)));
            atomicAdd(&atomForce[atom+PADDED_NUM_ATOMS], static_cast<unsigned long long>((long long) (fy*weight)));
            atomicAdd(&atomForce[atom+2*PADDED_NUM_ATOMS], static_cast<unsigned long long>((long long) (fz*weight)));
        }
    }
}
