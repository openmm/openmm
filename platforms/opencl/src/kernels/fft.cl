float2 multiplyComplex(float2 c1, float2 c2) {
    return (float2) (c1.x*c2.x-c1.y*c2.y, c1.x*c2.y+c1.y*c2.x);
}

/**
 * Perform a 1D FFT on each row along one axis.
 */

__kernel void execFFT(__global const float2* restrict in, __global float2* restrict out, float sign, __local float2* restrict w,
        __local float2* restrict data0, __local float2* restrict data1) {
    for (int i = get_local_id(0); i < ZSIZE; i += get_local_size(0))
        w[i] = (float2) (cos(-sign*i*2*M_PI/ZSIZE), sin(-sign*i*2*M_PI/ZSIZE));
    barrier(CLK_LOCAL_MEM_FENCE);
    for (int index = get_group_id(0); index < XSIZE*YSIZE; index += get_num_groups(0)) {
        int x = index/YSIZE;
        int y = index-x*YSIZE;
#ifdef LOOP_REQUIRED
        for (int z = get_local_id(0); z < ZSIZE; z += get_local_size(0))
            data0[z] = in[x*(YSIZE*ZSIZE)+y*ZSIZE+z];
#else
        data0[get_local_id(0)] = in[x*(YSIZE*ZSIZE)+y*ZSIZE+get_local_id(0)];
#endif
        barrier(CLK_LOCAL_MEM_FENCE);
        COMPUTE_FFT
    }
}
