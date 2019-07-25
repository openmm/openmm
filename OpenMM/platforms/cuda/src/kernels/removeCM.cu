/**
 * Calculate the center of mass momentum.
 */

extern "C" __global__ void calcCenterOfMassMomentum(int numAtoms, const mixed4* __restrict__ velm, float4* __restrict__ cmMomentum) {
    extern __shared__ volatile float3 temp[];
    float3 cm = make_float3(0, 0, 0);
    for (unsigned int index = blockIdx.x*blockDim.x+threadIdx.x; index < numAtoms; index += blockDim.x*gridDim.x) {
        mixed4 velocity = velm[index];
        if (velocity.w != 0) {
            mixed mass = RECIP(velocity.w);
            cm.x += (float) (velocity.x*mass);
            cm.y += (float) (velocity.y*mass);
            cm.z += (float) (velocity.z*mass);
        }
    }

    // Sum the threads in this group.

    int thread = threadIdx.x;
    temp[thread].x = cm.x;
    temp[thread].y = cm.y;
    temp[thread].z = cm.z;
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
    if (thread == 0) {
        float3 sum = make_float3(temp[thread].x+temp[thread+1].x, temp[thread].y+temp[thread+1].y, temp[thread].z+temp[thread+1].z);
        cmMomentum[blockIdx.x] = make_float4(sum.x, sum.y, sum.z, 0.0f);
    }
}

/**
 * Remove center of mass motion.
 */

extern "C" __global__ void removeCenterOfMassMomentum(unsigned int numAtoms, mixed4* __restrict__ velm, const float4* __restrict__ cmMomentum) {
    // First sum all of the momenta that were calculated by individual groups.

    extern volatile float3 temp[];
    float3 cm = make_float3(0, 0, 0);
    for (unsigned int index = threadIdx.x; index < gridDim.x; index += blockDim.x) {
        float4 momentum = cmMomentum[index];
        cm.x += momentum.x;
        cm.y += momentum.y;
        cm.z += momentum.z;
    }
    int thread = threadIdx.x;
    temp[thread].x = cm.x;
    temp[thread].y = cm.y;
    temp[thread].z = cm.z;
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
    __syncthreads();
    cm = make_float3(INVERSE_TOTAL_MASS*(temp[0].x+temp[1].x), INVERSE_TOTAL_MASS*(temp[0].y+temp[1].y), INVERSE_TOTAL_MASS*(temp[0].z+temp[1].z));

    // Now remove the center of mass velocity from each atom.

    for (unsigned int index = blockIdx.x*blockDim.x+threadIdx.x; index < numAtoms; index += blockDim.x*gridDim.x) {
        mixed4 velocity = velm[index];
        velocity.x -= cm.x;
        velocity.y -= cm.y;
        velocity.z -= cm.z;
        velm[index] = velocity;
    }
}
