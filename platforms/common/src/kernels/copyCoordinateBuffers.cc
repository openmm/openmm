KERNEL void copyFloatBuffer(GLOBAL float* RESTRICT source, GLOBAL float4* RESTRICT dest, int numAtoms) {
    for (int i = GLOBAL_ID; i < numAtoms; i += GLOBAL_SIZE) {
        dest[i].x = source[3*i];
        dest[i].y = source[3*i+1];
        dest[i].z = source[3*i+2];
    }
}

#ifdef SUPPORTS_DOUBLE_PRECISION
KERNEL void copyDoubleBuffer(GLOBAL double* RESTRICT source, GLOBAL double4* RESTRICT dest, int numAtoms) {
    for (int i = GLOBAL_ID; i < numAtoms; i += GLOBAL_SIZE) {
        dest[i].x = source[3*i];
        dest[i].y = source[3*i+1];
        dest[i].z = source[3*i+2];
    }
}
#endif