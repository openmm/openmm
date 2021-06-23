DEVICE mixed3 multiplyComplexRealPart(mixed2 c1, mixed3 c2r, mixed3 c2i) {
    return c1.x*c2r-c1.y*c2i;
}

DEVICE mixed3 multiplyComplexImagPart(mixed2 c1, mixed3 c2r, mixed3 c2i) {
    return c1.x*c2i+c1.y*c2r;
}

DEVICE mixed3 multiplyComplexRealPartConj(mixed2 c1, mixed3 c2r, mixed3 c2i) {
    return c1.x*c2r+c1.y*c2i;
}

DEVICE mixed3 multiplyComplexImagPartConj(mixed2 c1, mixed3 c2r, mixed3 c2i) {
    return c1.x*c2i-c1.y*c2r;
}

/**
 * Compute the contracted positions
 */
KERNEL void contractPositions(GLOBAL mixed4* posq, GLOBAL mixed4* contracted) {
    const int numBlocks = (GLOBAL_SIZE)/NUM_COPIES;
    const int blockStart = NUM_COPIES*(LOCAL_ID/NUM_COPIES);
    const int indexInBlock = LOCAL_ID-blockStart;
    LOCAL mixed3 q[2*THREAD_BLOCK_SIZE];
    LOCAL mixed3 temp[2*THREAD_BLOCK_SIZE];
    LOCAL mixed2 w1[NUM_COPIES];
    LOCAL mixed2 w2[NUM_CONTRACTED_COPIES];
    LOCAL_ARG mixed3* qreal = &q[blockStart];
    LOCAL_ARG mixed3* qimag = &q[blockStart+LOCAL_SIZE];
    LOCAL_ARG mixed3* tempreal = &temp[blockStart];
    LOCAL_ARG mixed3* tempimag = &temp[blockStart+LOCAL_SIZE];
    if (LOCAL_ID < NUM_COPIES)
        w1[indexInBlock] = make_mixed2(cos(-indexInBlock*2*M_PI/NUM_COPIES), sin(-indexInBlock*2*M_PI/NUM_COPIES));
    if (LOCAL_ID < NUM_CONTRACTED_COPIES)
        w2[indexInBlock] = make_mixed2(cos(-indexInBlock*2*M_PI/NUM_CONTRACTED_COPIES), sin(-indexInBlock*2*M_PI/NUM_CONTRACTED_COPIES));
    SYNC_THREADS;
    for (int particle = (GLOBAL_ID)/NUM_COPIES; particle < NUM_ATOMS; particle += numBlocks) {
        // Load the particle position.
        
        mixed4 particlePosq = posq[particle+indexInBlock*PADDED_NUM_ATOMS];
        qreal[indexInBlock] = make_mixed3(particlePosq.x, particlePosq.y, particlePosq.z);
        qimag[indexInBlock] = make_mixed3(0);
        
        // Forward FFT.
        
        SYNC_THREADS;
        LOCAL_ARG mixed2* w = w1;
        FFT_Q_FORWARD
        if (NUM_CONTRACTED_COPIES > 1) {
            // Compress the data to remove high frequencies.
            
            int start = (NUM_CONTRACTED_COPIES+1)/2;
            tempreal[indexInBlock] = qreal[indexInBlock];
            tempimag[indexInBlock] = qimag[indexInBlock];
            SYNC_THREADS;
            if (indexInBlock < NUM_CONTRACTED_COPIES) {
                qreal[indexInBlock] = tempreal[indexInBlock < start ? indexInBlock : indexInBlock+(NUM_COPIES-NUM_CONTRACTED_COPIES)];
                qimag[indexInBlock] = tempimag[indexInBlock < start ? indexInBlock : indexInBlock+(NUM_COPIES-NUM_CONTRACTED_COPIES)];
            }
            SYNC_THREADS;
            w = w2;
            FFT_Q_BACKWARD
        }
        
        // Store results.
        
        if (indexInBlock < NUM_CONTRACTED_COPIES)
            contracted[particle+indexInBlock*PADDED_NUM_ATOMS] = make_mixed4(POS_SCALE*qreal[indexInBlock].x, POS_SCALE*qreal[indexInBlock].y, POS_SCALE*qreal[indexInBlock].z, particlePosq.w);
    }
}

/**
 * Apply the contracted forces to all copies.
 */
KERNEL void contractForces(GLOBAL mm_long* force, GLOBAL mm_long* contracted) {
    const int numBlocks = (GLOBAL_SIZE)/NUM_COPIES;
    const int blockStart = NUM_COPIES*(LOCAL_ID/NUM_COPIES);
    const int indexInBlock = LOCAL_ID-blockStart;
    const mixed forceScale = 1/(mixed) 0x100000000;
    LOCAL mixed3 f[2*THREAD_BLOCK_SIZE];
    LOCAL mixed3 temp[2*THREAD_BLOCK_SIZE];
    LOCAL mixed2 w1[NUM_COPIES];
    LOCAL mixed2 w2[NUM_CONTRACTED_COPIES];
    LOCAL_ARG mixed3* freal = &f[blockStart];
    LOCAL_ARG mixed3* fimag = &f[blockStart+LOCAL_SIZE];
    LOCAL_ARG mixed3* tempreal = &temp[blockStart];
    LOCAL_ARG mixed3* tempimag = &temp[blockStart+LOCAL_SIZE];
    if (LOCAL_ID < NUM_COPIES)
        w1[indexInBlock] = make_mixed2(cos(-indexInBlock*2*M_PI/NUM_COPIES), sin(-indexInBlock*2*M_PI/NUM_COPIES));
    if (LOCAL_ID < NUM_CONTRACTED_COPIES)
        w2[indexInBlock] = make_mixed2(cos(-indexInBlock*2*M_PI/NUM_CONTRACTED_COPIES), sin(-indexInBlock*2*M_PI/NUM_CONTRACTED_COPIES));
    SYNC_THREADS;
    for (int particle = (GLOBAL_ID)/NUM_COPIES; particle < NUM_ATOMS; particle += numBlocks) {
        // Load the force.
        
        int forceIndex = particle+indexInBlock*PADDED_NUM_ATOMS*3;
        if (indexInBlock < NUM_CONTRACTED_COPIES) {
            freal[indexInBlock] = make_mixed3(contracted[forceIndex]*forceScale, contracted[forceIndex+PADDED_NUM_ATOMS]*forceScale, contracted[forceIndex+PADDED_NUM_ATOMS*2]*forceScale);
            fimag[indexInBlock] = make_mixed3(0);
        }
        SYNC_THREADS;

        // Forward FFT.
        
        LOCAL_ARG mixed2* w = w2;
        if (NUM_CONTRACTED_COPIES > 1) {
            FFT_F_FORWARD
        }
        
        // Set the high frequency components to 0.
        
        int start = (NUM_CONTRACTED_COPIES+1)/2;
        int end = NUM_COPIES-NUM_CONTRACTED_COPIES+start;
        tempreal[indexInBlock] = freal[indexInBlock];
        tempimag[indexInBlock] = fimag[indexInBlock];
        SYNC_THREADS;
        if (indexInBlock >= start) {
            freal[indexInBlock] = (indexInBlock < end ? make_mixed3(0) : tempreal[indexInBlock-(NUM_COPIES-NUM_CONTRACTED_COPIES)]);
            fimag[indexInBlock] = (indexInBlock < end ? make_mixed3(0) : tempimag[indexInBlock-(NUM_COPIES-NUM_CONTRACTED_COPIES)]);
        }
        SYNC_THREADS;
        w = w1;
        FFT_F_BACKWARD
        
        // Store results.
        
        force[forceIndex] += (mm_long) (FORCE_SCALE*freal[indexInBlock].x);
        force[forceIndex+PADDED_NUM_ATOMS] += (mm_long) (FORCE_SCALE*freal[indexInBlock].y);
        force[forceIndex+PADDED_NUM_ATOMS*2] += (mm_long) (FORCE_SCALE*freal[indexInBlock].z);
    }
}
