KERNEL void copyFloatBuffer(GLOBAL float* RESTRICT source, GLOBAL float4* RESTRICT dest, int numAtoms) {
    for (int i = GLOBAL_ID; i < numAtoms; i += GLOBAL_SIZE) {
        dest[i].x = source[3*i];
        dest[i].y = source[3*i+1];
        dest[i].z = source[3*i+2];
    }
}

// This kernel is actually used for copying double precision buffers, but we list the type as long (which has the same
// size as double) so it will compile and work on platforms that don't support double precision.  Those platforms may
// emulate double precision with a pair of floats (which also has the same size).
KERNEL void copyDoubleBuffer(GLOBAL mm_long* RESTRICT source, GLOBAL mm_long4* RESTRICT dest, int numAtoms) {
    for (int i = GLOBAL_ID; i < numAtoms; i += GLOBAL_SIZE) {
        dest[i].x = source[3*i];
        dest[i].y = source[3*i+1];
        dest[i].z = source[3*i+2];
    }
}
