/**
 * Calculate the center of mass momentum.
 */

__kernel void calcCenterOfMassMomentum(int numAtoms, __global float4* velm, __global float4* cmMomentum, __local float4* temp) {
    int index = get_global_id(0);
    float4 cm = 0.0f;
    while (index < numAtoms) {
        float4 velocity = velm[index];
        cm.xyz += velocity.xyz/velocity.w;
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

__kernel void removeCenterOfMassMomentum(int numAtoms, __global float4* velm, __global float4* cmMomentum, __local float4* temp) {
    // First sum all of the momenta that were calculated by individual groups.

    int index = get_local_id(0);
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
    cm = (temp[0]+temp[1])*INVERSE_TOTAL_MASS;

    // Now remove the center of mass velocity from each atom.

    index = get_global_id(0);
    while (index < numAtoms) {
        velm[index].xyz -= cm.xyz;
        index += get_global_size(0);
    }
}
