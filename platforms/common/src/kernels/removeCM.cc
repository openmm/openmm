// Optional CUDA/Volta+ reduction.  Both callers launch exactly 64 threads.
// One CTA barrier transfers lanes 32..63, then warp 0 retains the exact original
// 32,16,8,4,2,1 binary32 addition tree.  Only lane 0's return value is consumed.
#if defined(CMM_WARP_REDUCTION) && defined(CMM_REDUCE_ONCE) && defined(__CUDA_ARCH__) && __CUDA_ARCH__ >= 700
#define CMM_USE_WARP_REDUCTION 1
inline DEVICE float3 reduceCMMomentumWarp64(float3 value, LOCAL_ARG float4* temp) {
    int thread = LOCAL_ID;
    temp[thread] = make_float4(value.x, value.y, value.z, 0);
    SYNC_THREADS;
    if (thread < 32) {
        value.x += temp[thread+32].x;
        value.y += temp[thread+32].y;
        value.z += temp[thread+32].z;
        for (int offset = 16; offset >= 1; offset >>= 1) {
            // Every lane in warp 0 reaches every shuffle with the full mask.
            float x = __shfl_down_sync(0xffffffff, value.x, offset);
            float y = __shfl_down_sync(0xffffffff, value.y, offset);
            float z = __shfl_down_sync(0xffffffff, value.z, offset);
            if (thread < offset) {
                value.x += x;
                value.y += y;
                value.z += z;
            }
        }
    }
    return value;
}
#endif

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

#ifdef CMM_USE_WARP_REDUCTION
    float3 total = reduceCMMomentumWarp64(make_float3(cm.x, cm.y, cm.z), temp);
    if (LOCAL_ID == 0)
        cmMomentum[GROUP_ID] = make_float4(total.x, total.y, total.z, 0);
#else
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
#endif
}

/**
 * Remove center of mass motion.
 */

#ifdef CMM_REDUCE_ONCE
KERNEL void removeCenterOfMassMomentum(int numAtoms, GLOBAL mixed4* RESTRICT velm, GLOBAL float4* RESTRICT cmMomentum) {
#else
KERNEL void removeCenterOfMassMomentum(int numAtoms, GLOBAL mixed4* RESTRICT velm, GLOBAL const float4* RESTRICT cmMomentum) {
#endif
    // First sum all of the momenta that were calculated by individual groups.

    LOCAL float4 temp[64];
    float4 cm = make_float4(0);
#ifdef CMM_REDUCE_ONCE
    for (int index = LOCAL_ID; index < CMM_ORIGINAL_NUM_GROUPS; index += LOCAL_SIZE)
#else
    for (int index = LOCAL_ID; index < NUM_GROUPS; index += LOCAL_SIZE)
#endif
        cm += cmMomentum[index];
#ifdef CMM_USE_WARP_REDUCTION
    int thread = LOCAL_ID;
    float3 total = reduceCMMomentumWarp64(make_float3(cm.x, cm.y, cm.z), temp);
    cm = make_float4(INVERSE_TOTAL_MASS*total.x, INVERSE_TOTAL_MASS*total.y, INVERSE_TOTAL_MASS*total.z, 0);
#else
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
#endif

#ifdef CMM_REDUCE_ONCE
    // This kernel runs as one 64-thread block.  Preserve the complete reduction
    // order above; do not cache across calls or alter when velocities change.
    if (thread == 0)
        cmMomentum[CMM_ORIGINAL_NUM_GROUPS] = cm;
#else
    // Now remove the center of mass velocity from each atom.

    for (int index = GLOBAL_ID; index < numAtoms; index += GLOBAL_SIZE) {
        velm[index].x -= cm.x;
        velm[index].y -= cm.y;
        velm[index].z -= cm.z;
    }
#endif
}

#ifdef CMM_REDUCE_ONCE
KERNEL void applyCenterOfMassVelocity(int numAtoms, GLOBAL mixed4* RESTRICT velm, GLOBAL const float4* RESTRICT cmMomentum) {
    float4 cm = cmMomentum[CMM_ORIGINAL_NUM_GROUPS];
    for (int index = GLOBAL_ID; index < numAtoms; index += GLOBAL_SIZE) {
        velm[index].x -= cm.x;
        velm[index].y -= cm.y;
        velm[index].z -= cm.z;
    }
}
#endif
