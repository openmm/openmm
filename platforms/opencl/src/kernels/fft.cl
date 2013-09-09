real2 multiplyComplex(real2 c1, real2 c2) {
    return (real2) (c1.x*c2.x-c1.y*c2.y, c1.x*c2.y+c1.y*c2.x);
}

/**
 * Perform a 1D FFT on each row along one axis.
 */

__kernel void execFFT(__global const real2* restrict in, __global real2* restrict out, int sign, __local real2* restrict w,
        __local real2* restrict data0, __local real2* restrict data1) {
    for (int i = get_local_id(0); i < ZSIZE; i += get_local_size(0))
        w[i] = (real2) (cos(-sign*i*2*M_PI/ZSIZE), sin(-sign*i*2*M_PI/ZSIZE));
    barrier(CLK_LOCAL_MEM_FENCE);
    
    for (int baseIndex = get_group_id(0)*BLOCKS_PER_GROUP; baseIndex < XSIZE*YSIZE; baseIndex += get_num_groups(0)*BLOCKS_PER_GROUP) {
        int index = baseIndex+get_local_id(0)/ZSIZE;
        int x = index/YSIZE;
        int y = index-x*YSIZE;
#if LOOP_REQUIRED
        for (int z = get_local_id(0); z < ZSIZE; z += get_local_size(0))
            data0[z] = in[x*(YSIZE*ZSIZE)+y*ZSIZE+z];
#else
        if (index < XSIZE*YSIZE)
            data0[get_local_id(0)] = in[x*(YSIZE*ZSIZE)+y*ZSIZE+get_local_id(0)%ZSIZE];
#endif
        barrier(CLK_LOCAL_MEM_FENCE);
        COMPUTE_FFT
    }
}
