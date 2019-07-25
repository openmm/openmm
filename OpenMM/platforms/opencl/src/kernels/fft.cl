real2 multiplyComplex(real2 c1, real2 c2) {
    return (real2) (c1.x*c2.x-c1.y*c2.y, c1.x*c2.y+c1.y*c2.x);
}

/**
 * Load a value from the half-complex grid produces by a real-to-complex transform.
 */
real2 loadComplexValue(__global const real2* restrict in, int x, int y, int z) {
    const int inputZSize = ZSIZE/2+1;
    if (z < inputZSize)
        return in[x*YSIZE*inputZSize+y*inputZSize+z];
    int xp = (x == 0 ? 0 : XSIZE-x);
    int yp = (y == 0 ? 0 : YSIZE-y);
    real2 value = in[xp*YSIZE*inputZSize+yp*inputZSize+(ZSIZE-z)];
    return (real2) (value.x, -value.y);
}

/**
 * Perform a 1D FFT on each row along one axis.
 */

__kernel void execFFT(__global const INPUT_TYPE* restrict in, __global OUTPUT_TYPE* restrict out, __local real2* restrict w,
        __local real2* restrict data0, __local real2* restrict data1) {
    for (int i = get_local_id(0); i < ZSIZE; i += get_local_size(0))
        w[i] = (real2) (cos(-(SIGN)*i*2*M_PI/ZSIZE), sin(-(SIGN)*i*2*M_PI/ZSIZE));
    barrier(CLK_LOCAL_MEM_FENCE);
    
    for (int baseIndex = get_group_id(0)*BLOCKS_PER_GROUP; baseIndex < XSIZE*YSIZE; baseIndex += get_num_groups(0)*BLOCKS_PER_GROUP) {
        int index = baseIndex+get_local_id(0)/ZSIZE;
        int x = index/YSIZE;
        int y = index-x*YSIZE;
#if OUTPUT_IS_PACKED
        if (x < XSIZE/2+1) {
#endif
#if LOOP_REQUIRED
        for (int z = get_local_id(0); z < ZSIZE; z += get_local_size(0))
    #if INPUT_IS_REAL
            data0[z] = (real2) (in[x*(YSIZE*ZSIZE)+y*ZSIZE+z], 0);
    #elif INPUT_IS_PACKED
            data0[z] = loadComplexValue(in, x, y, z);
    #else
            data0[z] = in[x*(YSIZE*ZSIZE)+y*ZSIZE+z];
    #endif
#else
        if (index < XSIZE*YSIZE)
    #if INPUT_IS_REAL
            data0[get_local_id(0)] = (real2) (in[x*(YSIZE*ZSIZE)+y*ZSIZE+get_local_id(0)%ZSIZE], 0);
    #elif INPUT_IS_PACKED
            data0[get_local_id(0)] = loadComplexValue(in, x, y, get_local_id(0)%ZSIZE);
    #else
            data0[get_local_id(0)] = in[x*(YSIZE*ZSIZE)+y*ZSIZE+get_local_id(0)%ZSIZE];
    #endif
#endif
#if OUTPUT_IS_PACKED
        }
#endif
        barrier(CLK_LOCAL_MEM_FENCE);
        COMPUTE_FFT
    }
}
