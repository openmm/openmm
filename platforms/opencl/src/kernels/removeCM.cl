/**
 * Calculate the center of mass momentum.
 */

__kernel void calcCenterOfMassMomentum(int numAtoms, __global const mixed4* restrict velm, __global float4* restrict cmMomentum, __local volatile float4* restrict temp) {
    int index = get_global_id(0);
    float4 cm = 0.0f;
    while (index < numAtoms) {
        mixed4 velocity = velm[index];
        if (velocity.w != 0) {
            cm.x += velocity.x/velocity.w;
            cm.y += velocity.y/velocity.w;
            cm.z += velocity.z/velocity.w;
        }
        index += get_global_size(0);
    }

    // Sum the threads in this group.

    int thread = get_local_id(0);
    temp[thread] = cm;
    barrier(CLK_LOCAL_MEM_FENCE);
#ifdef WARPS_ARE_ATOMIC
    if (thread < 32) {
        temp[thread] += temp[thread+32];
        if (thread < 16)
            temp[thread] += temp[thread+16];
        if (thread < 8)
            temp[thread] += temp[thread+8];
        if (thread < 4)
            temp[thread] += temp[thread+4];
        if (thread < 2)
            temp[thread] += temp[thread+2];
    }
#else
    if (thread < 32)
        temp[thread] += temp[thread+32];
    barrier(CLK_LOCAL_MEM_FENCE);
    if (thread < 16)
        temp[thread] += temp[thread+16];
    barrier(CLK_LOCAL_MEM_FENCE);
    if (thread < 8)
        temp[thread] += temp[thread+8];
    barrier(CLK_LOCAL_MEM_FENCE);
    if (thread < 4)
        temp[thread] += temp[thread+4];
    barrier(CLK_LOCAL_MEM_FENCE);
    if (thread < 2)
        temp[thread] += temp[thread+2];
    barrier(CLK_LOCAL_MEM_FENCE);
#endif
    if (thread == 0)
        cmMomentum[get_group_id(0)] = temp[thread]+temp[thread+1];
}

/**
 * Remove center of mass motion.
 */

__kernel void removeCenterOfMassMomentum(unsigned int numAtoms, __global mixed4* restrict velm, __global const float4* restrict cmMomentum, __local volatile float4* restrict temp) {
    // First sum all of the momenta that were calculated by individual groups.

    unsigned int index = get_local_id(0);
    float4 cm = 0.0f;
    while (index < get_num_groups(0)) {
        cm += cmMomentum[index];
        index += get_local_size(0);
    }
    int thread = get_local_id(0);
    temp[thread] = cm;
    barrier(CLK_LOCAL_MEM_FENCE);
#ifdef WARPS_ARE_ATOMIC
    if (thread < 32) {
        temp[thread] += temp[thread+32];
        if (thread < 16)
            temp[thread] += temp[thread+16];
        if (thread < 8)
            temp[thread] += temp[thread+8];
        if (thread < 4)
            temp[thread] += temp[thread+4];
        if (thread < 2)
            temp[thread] += temp[thread+2];
    }
#else
    if (thread < 32)
        temp[thread] += temp[thread+32];
    barrier(CLK_LOCAL_MEM_FENCE);
    if (thread < 16)
        temp[thread] += temp[thread+16];
    barrier(CLK_LOCAL_MEM_FENCE);
    if (thread < 8)
        temp[thread] += temp[thread+8];
    barrier(CLK_LOCAL_MEM_FENCE);
    if (thread < 4)
        temp[thread] += temp[thread+4];
    barrier(CLK_LOCAL_MEM_FENCE);
    if (thread < 2)
        temp[thread] += temp[thread+2];
#endif
    barrier(CLK_LOCAL_MEM_FENCE);
    cm = (float4) (INVERSE_TOTAL_MASS*(temp[0].x+temp[1].x), INVERSE_TOTAL_MASS*(temp[0].y+temp[1].y), INVERSE_TOTAL_MASS*(temp[0].z+temp[1].z), 0);

    // Now remove the center of mass velocity from each atom.

    index = get_global_id(0);
    while (index < numAtoms) {
        velm[index].x -= cm.x;
        velm[index].y -= cm.y;
        velm[index].z -= cm.z;
        index += get_global_size(0);
    }
}
