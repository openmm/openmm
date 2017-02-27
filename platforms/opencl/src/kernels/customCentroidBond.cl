#pragma OPENCL EXTENSION cl_khr_int64_base_atomics : enable

/**
 * Compute the center of each group.
 */
__kernel void computeGroupCenters(__global const real4* restrict posq, __global const int* restrict groupParticles,
        __global const real* restrict groupWeights, __global const int* restrict groupOffsets, __global real4* restrict centerPositions) {
    __local volatile real3 temp[64];
    for (int group = get_group_id(0); group < NUM_GROUPS; group += get_num_groups(0)) {
        // The threads in this block work together to compute the center one group.

        int firstIndex = groupOffsets[group];
        int lastIndex = groupOffsets[group+1];
        real3 center = (real3) 0;
        for (int index = get_local_id(0); index < lastIndex-firstIndex; index += get_local_size(0)) {
            int atom = groupParticles[firstIndex+index];
            real weight = groupWeights[firstIndex+index];
            real4 pos = posq[atom];
            center.x += weight*pos.x;
            center.y += weight*pos.y;
            center.z += weight*pos.z;
        }

        // Sum the values.

        int thread = get_local_id(0);
        temp[thread].x = center.x;
        temp[thread].y = center.y;
        temp[thread].z = center.z;

        barrier(CLK_LOCAL_MEM_FENCE);
        if (thread < 32) {
            temp[thread].x += temp[thread+32].x;
            temp[thread].y += temp[thread+32].y;
            temp[thread].z += temp[thread+32].z;
        }

        SYNC_WARPS;
        if (thread < 16) {
            temp[thread].x += temp[thread+16].x;
            temp[thread].y += temp[thread+16].y;
            temp[thread].z += temp[thread+16].z;
        }
        SYNC_WARPS;
        if (thread < 8) {
            temp[thread].x += temp[thread+8].x;
            temp[thread].y += temp[thread+8].y;
            temp[thread].z += temp[thread+8].z;
        }

        SYNC_WARPS;
        if (thread < 4) {
            temp[thread].x += temp[thread+4].x;
            temp[thread].y += temp[thread+4].y;
            temp[thread].z += temp[thread+4].z;
        }
        SYNC_WARPS;
        if (thread < 2) {
            temp[thread].x += temp[thread+2].x;
            temp[thread].y += temp[thread+2].y;
            temp[thread].z += temp[thread+2].z;
        }

        SYNC_WARPS;
        if (thread == 0)
            centerPositions[group] = (real4) (temp[0].x+temp[1].x, temp[0].y+temp[1].y, temp[0].z+temp[1].z, 0);
    }
}

/**
 * Compute the difference between two vectors, setting the fourth component to the squared magnitude.
 */
real4 delta(real4 vec1, real4 vec2, bool periodic, real4 periodicBoxSize, real4 invPeriodicBoxSize, 
        real4 periodicBoxVecX, real4 periodicBoxVecY, real4 periodicBoxVecZ) {
    real4 result = (real4) (vec1.x-vec2.x, vec1.y-vec2.y, vec1.z-vec2.z, 0);
    if (periodic)
        APPLY_PERIODIC_TO_DELTA(result);
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
        if (cosine < 0)
            angle = M_PI-angle;
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
 * Compute the forces on groups based on the bonds.
 */
__kernel void computeGroupForces(__global long* restrict groupForce, __global mixed* restrict energyBuffer, __global const real4* restrict centerPositions,
        __global const int* restrict bondGroups, real4 periodicBoxSize, real4 invPeriodicBoxSize, real4 periodicBoxVecX, real4 periodicBoxVecY, real4 periodicBoxVecZ
        EXTRA_ARGS) {
    mixed energy = 0;
    INIT_PARAM_DERIVS
    for (int index = get_global_id(0); index < NUM_BONDS; index += get_global_size(0)) {
        COMPUTE_FORCE
    }
    energyBuffer[get_global_id(0)] += energy;
    SAVE_PARAM_DERIVS
}

/**
 * Apply the forces from the group centers to the individual atoms.
 */
__kernel void applyForcesToAtoms(__global const int* restrict groupParticles, __global const real* restrict groupWeights, __global const int* restrict groupOffsets,
        __global const long* restrict groupForce, __global long* restrict atomForce) {
    for (int group = get_group_id(0); group < NUM_GROUPS; group += get_num_groups(0)) {
        long fx = groupForce[group];
        long fy = groupForce[group+NUM_GROUPS];
        long fz = groupForce[group+NUM_GROUPS*2];
        int firstIndex = groupOffsets[group];
        int lastIndex = groupOffsets[group+1];
        for (int index = get_local_id(0); index < lastIndex-firstIndex; index += get_local_size(0)) {
            int atom = groupParticles[firstIndex+index];
            real weight = groupWeights[firstIndex+index];
            atom_add(&atomForce[atom], (long) (fx*weight));
            atom_add(&atomForce[atom+PADDED_NUM_ATOMS], (long) (fy*weight));
            atom_add(&atomForce[atom+2*PADDED_NUM_ATOMS], (long) (fz*weight));
        }
    }
}
