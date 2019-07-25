mixed4 multiplyComplexRealPart(mixed2 c1, mixed4 c2r, mixed4 c2i) {
    return c1.x*c2r-c1.y*c2i;
}

mixed4 multiplyComplexImagPart(mixed2 c1, mixed4 c2r, mixed4 c2i) {
    return c1.x*c2i+c1.y*c2r;
}

mixed4 multiplyComplexRealPartConj(mixed2 c1, mixed4 c2r, mixed4 c2i) {
    return c1.x*c2r+c1.y*c2i;
}

mixed4 multiplyComplexImagPartConj(mixed2 c1, mixed4 c2r, mixed4 c2i) {
    return c1.x*c2i-c1.y*c2r;
}

/**
 * Compute the contracted positions
 */
__kernel void contractPositions(__global mixed4* posq, __global mixed4* contracted) {
    const int numBlocks = get_global_size(0)/NUM_COPIES;
    const int blockStart = NUM_COPIES*(get_local_id(0)/NUM_COPIES);
    const int indexInBlock = get_local_id(0)-blockStart;
    __local mixed4 q[2*THREAD_BLOCK_SIZE];
    __local mixed4 temp[2*THREAD_BLOCK_SIZE];
    __local mixed2 w1[NUM_COPIES];
    __local mixed2 w2[NUM_CONTRACTED_COPIES];
    __local mixed4* qreal = &q[blockStart];
    __local mixed4* qimag = &q[blockStart+get_local_size(0)];
    __local mixed4* tempreal = &temp[blockStart];
    __local mixed4* tempimag = &temp[blockStart+get_local_size(0)];
    if (get_local_id(0) < NUM_COPIES)
        w1[indexInBlock] = (mixed2) (cos(-indexInBlock*2*M_PI/NUM_COPIES), sin(-indexInBlock*2*M_PI/NUM_COPIES));
    if (get_local_id(0) < NUM_CONTRACTED_COPIES)
        w2[indexInBlock] = (mixed2) (cos(-indexInBlock*2*M_PI/NUM_CONTRACTED_COPIES), sin(-indexInBlock*2*M_PI/NUM_CONTRACTED_COPIES));
    barrier(CLK_LOCAL_MEM_FENCE);
    for (int particle = get_global_id(0)/NUM_COPIES; particle < NUM_ATOMS; particle += numBlocks) {
        // Load the particle position.
        
        mixed4 particlePosq = convert_mixed4(posq[particle+indexInBlock*PADDED_NUM_ATOMS]);
        qreal[indexInBlock] = particlePosq;
        qimag[indexInBlock] = (mixed4) (0.0f, 0.0f, 0.0f, 0.0f);
        
        // Forward FFT.
        
        barrier(CLK_LOCAL_MEM_FENCE);
        __local mixed2* w = w1;
        FFT_Q_FORWARD
        if (NUM_CONTRACTED_COPIES > 1) {
            // Compress the data to remove high frequencies.
            
            int start = (NUM_CONTRACTED_COPIES+1)/2;
            tempreal[indexInBlock] = qreal[indexInBlock];
            tempimag[indexInBlock] = qimag[indexInBlock];
            barrier(CLK_LOCAL_MEM_FENCE);
            if (indexInBlock < NUM_CONTRACTED_COPIES) {
                qreal[indexInBlock] = tempreal[indexInBlock < start ? indexInBlock : indexInBlock+(NUM_COPIES-NUM_CONTRACTED_COPIES)];
                qimag[indexInBlock] = tempimag[indexInBlock < start ? indexInBlock : indexInBlock+(NUM_COPIES-NUM_CONTRACTED_COPIES)];
            }
            barrier(CLK_LOCAL_MEM_FENCE);
            w = w2;
            FFT_Q_BACKWARD
        }
        
        // Store results.
        
        if (indexInBlock < NUM_CONTRACTED_COPIES)
            contracted[particle+indexInBlock*PADDED_NUM_ATOMS] = (mixed4) (POS_SCALE*qreal[indexInBlock].x, POS_SCALE*qreal[indexInBlock].y, POS_SCALE*qreal[indexInBlock].z, particlePosq.w);
    }
}

/**
 * Apply the contracted forces to all copies.
 */
__kernel void contractForces(__global real4* force, __global real4* contracted) {
    const int numBlocks = get_global_size(0)/NUM_COPIES;
    const int blockStart = NUM_COPIES*(get_local_id(0)/NUM_COPIES);
    const int indexInBlock = get_local_id(0)-blockStart;
    __local mixed4 f[2*THREAD_BLOCK_SIZE];
    __local mixed4 temp[2*THREAD_BLOCK_SIZE];
    __local mixed2 w1[NUM_COPIES];
    __local mixed2 w2[NUM_CONTRACTED_COPIES];
    __local mixed4* freal = &f[blockStart];
    __local mixed4* fimag = &f[blockStart+get_local_size(0)];
    __local mixed4* tempreal = &temp[blockStart];
    __local mixed4* tempimag = &temp[blockStart+get_local_size(0)];
    if (get_local_id(0) < NUM_COPIES)
        w1[indexInBlock] = (mixed2) (cos(-indexInBlock*2*M_PI/NUM_COPIES), sin(-indexInBlock*2*M_PI/NUM_COPIES));
    if (get_local_id(0) < NUM_CONTRACTED_COPIES)
        w2[indexInBlock] = (mixed2) (cos(-indexInBlock*2*M_PI/NUM_CONTRACTED_COPIES), sin(-indexInBlock*2*M_PI/NUM_CONTRACTED_COPIES));
    barrier(CLK_LOCAL_MEM_FENCE);
    for (int particle = get_global_id(0)/NUM_COPIES; particle < NUM_ATOMS; particle += numBlocks) {
        // Load the force.
        
        int index = particle+indexInBlock*PADDED_NUM_ATOMS;
        if (indexInBlock < NUM_CONTRACTED_COPIES) {
            freal[indexInBlock] = convert_mixed4(contracted[index]);
            fimag[indexInBlock] = (mixed4) (0.0f, 0.0f, 0.0f, 0.0f);
        }
        barrier(CLK_LOCAL_MEM_FENCE);

        // Forward FFT.
        
        __local mixed2* w = w2;
        if (NUM_CONTRACTED_COPIES > 1) {
            FFT_F_FORWARD
        }
        
        // Set the high frequency components to 0.
        
        int start = (NUM_CONTRACTED_COPIES+1)/2;
        int end = NUM_COPIES-NUM_CONTRACTED_COPIES+start;
        tempreal[indexInBlock] = freal[indexInBlock];
        tempimag[indexInBlock] = fimag[indexInBlock];
        barrier(CLK_LOCAL_MEM_FENCE);
        if (indexInBlock >= start) {
            freal[indexInBlock] = (indexInBlock < end ? (mixed4) (0.0f, 0.0f, 0.0f, 0.0f) : tempreal[indexInBlock-(NUM_COPIES-NUM_CONTRACTED_COPIES)]);
            fimag[indexInBlock] = (indexInBlock < end ? (mixed4) (0.0f, 0.0f, 0.0f, 0.0f) : tempimag[indexInBlock-(NUM_COPIES-NUM_CONTRACTED_COPIES)]);
        }
        barrier(CLK_LOCAL_MEM_FENCE);
        w = w1;
        FFT_F_BACKWARD
        
        // Store results.

        force[index] += convert_real4(FORCE_SCALE*freal[indexInBlock]);
    }
}
