/**
 * Calculate the center of mass momentum.
 */

KERNEL void calcCenterOfMassMomentum(int numAtoms, GLOBAL const mixed4* RESTRICT velm, GLOBAL float4* RESTRICT cmMomentum) {
    LOCAL float4 temp[64];
    float4 cm = make_float4(0);
    for (int index = GLOBAL_ID; index < numAtoms; index += GLOBAL_SIZE) {
        mixed4 velocity = velm[index];
        if (velocity.w != 0) {
            mixed mass = RECIP(velocity.w);
            cm.x += (float) (velocity.x*mass);
            cm.y += (float) (velocity.y*mass);
            cm.z += (float) (velocity.z*mass);
        }
    }

    // Sum the threads in this group.

    int thread = LOCAL_ID;
    temp[thread] = cm;
    SYNC_THREADS;
    if (thread < 32)
        temp[thread] += temp[thread+32];
    SYNC_THREADS;
    if (thread < 16)
        temp[thread] += temp[thread+16];
    SYNC_THREADS;
    if (thread < 8)
        temp[thread] += temp[thread+8];
    SYNC_THREADS;
    if (thread < 4)
        temp[thread] += temp[thread+4];
    SYNC_THREADS;
    if (thread < 2)
        temp[thread] += temp[thread+2];
    SYNC_THREADS;
    if (thread == 0)
        cmMomentum[GROUP_ID] = temp[thread]+temp[thread+1];
}

/**
 * Remove center of mass motion.
 */

KERNEL void removeCenterOfMassMomentum(int numAtoms, GLOBAL mixed4* RESTRICT velm, GLOBAL const float4* RESTRICT cmMomentum) {
    // First sum all of the momenta that were calculated by individual groups.

    LOCAL float4 temp[64];
    float4 cm = make_float4(0);
    for (int index = LOCAL_ID; index < NUM_GROUPS; index += LOCAL_SIZE)
        cm += cmMomentum[index];
    int thread = LOCAL_ID;
    temp[thread] = cm;
    SYNC_THREADS;
    if (thread < 32)
        temp[thread] += temp[thread+32];
    SYNC_THREADS;
    if (thread < 16)
        temp[thread] += temp[thread+16];
    SYNC_THREADS;
    if (thread < 8)
        temp[thread] += temp[thread+8];
    SYNC_THREADS;
    if (thread < 4)
        temp[thread] += temp[thread+4];
    SYNC_THREADS;
    if (thread < 2)
        temp[thread] += temp[thread+2];
    SYNC_THREADS;
    cm = make_float4(INVERSE_TOTAL_MASS*(temp[0].x+temp[1].x), INVERSE_TOTAL_MASS*(temp[0].y+temp[1].y), INVERSE_TOTAL_MASS*(temp[0].z+temp[1].z), 0);

    // Now remove the center of mass velocity from each atom.

    for (int index = GLOBAL_ID; index < numAtoms; index += GLOBAL_SIZE) {
        velm[index].x -= cm.x;
        velm[index].y -= cm.y;
        velm[index].z -= cm.z;
    }
}
