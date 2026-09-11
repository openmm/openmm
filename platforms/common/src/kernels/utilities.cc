/**
 * This is called by the various kernels below to clear a buffer.
 */
DEVICE void clearSingleBuffer(GLOBAL int* RESTRICT buffer, int size) {
    int index = GLOBAL_ID;
    GLOBAL int4* buffer4 = (GLOBAL int4*) buffer;
    int sizeDiv4 = size/4;
    while (index < sizeDiv4) {
        buffer4[index] = make_int4(0);
        index += GLOBAL_SIZE;
    }
    if (GLOBAL_ID == 0)
        for (int i = sizeDiv4*4; i < size; i++)
            buffer[i] = 0;
}

/**
 * Fill a buffer with 0.
 */
KERNEL void clearBuffer(GLOBAL int* RESTRICT buffer, int size) {
    clearSingleBuffer(buffer, size);
}

/**
 * Fill two buffers with 0.
 */
KERNEL void clearTwoBuffers(GLOBAL int* RESTRICT buffer1, int size1, GLOBAL int* RESTRICT buffer2, int size2) {
    clearSingleBuffer(buffer1, size1);
    clearSingleBuffer(buffer2, size2);
}

/**
 * Fill three buffers with 0.
 */
KERNEL void clearThreeBuffers(GLOBAL int* RESTRICT buffer1, int size1, GLOBAL int* RESTRICT buffer2, int size2, GLOBAL int* RESTRICT buffer3, int size3) {
    clearSingleBuffer(buffer1, size1);
    clearSingleBuffer(buffer2, size2);
    clearSingleBuffer(buffer3, size3);
}

/**
 * Fill four buffers with 0.
 */
KERNEL void clearFourBuffers(GLOBAL int* RESTRICT buffer1, int size1, GLOBAL int* RESTRICT buffer2, int size2, GLOBAL int* RESTRICT buffer3, int size3, GLOBAL int* RESTRICT buffer4, int size4) {
    clearSingleBuffer(buffer1, size1);
    clearSingleBuffer(buffer2, size2);
    clearSingleBuffer(buffer3, size3);
    clearSingleBuffer(buffer4, size4);
}

/**
 * Fill five buffers with 0.
 */
KERNEL void clearFiveBuffers(GLOBAL int* RESTRICT buffer1, int size1, GLOBAL int* RESTRICT buffer2, int size2, GLOBAL int* RESTRICT buffer3, int size3, GLOBAL int* RESTRICT buffer4, int size4, GLOBAL int* RESTRICT buffer5, int size5) {
    clearSingleBuffer(buffer1, size1);
    clearSingleBuffer(buffer2, size2);
    clearSingleBuffer(buffer3, size3);
    clearSingleBuffer(buffer4, size4);
    clearSingleBuffer(buffer5, size5);
}

/**
 * Fill six buffers with 0.
 */
KERNEL void clearSixBuffers(GLOBAL int* RESTRICT buffer1, int size1, GLOBAL int* RESTRICT buffer2, int size2, GLOBAL int* RESTRICT buffer3, int size3, GLOBAL int* RESTRICT buffer4, int size4, GLOBAL int* RESTRICT buffer5, int size5, GLOBAL int* RESTRICT buffer6, int size6) {
    clearSingleBuffer(buffer1, size1);
    clearSingleBuffer(buffer2, size2);
    clearSingleBuffer(buffer3, size3);
    clearSingleBuffer(buffer4, size4);
    clearSingleBuffer(buffer5, size5);
    clearSingleBuffer(buffer6, size6);
}

/**
 * Sum the energy buffer.
 */
KERNEL void reduceEnergy(GLOBAL const mixed* RESTRICT energyBuffer, GLOBAL mixed* RESTRICT result, int bufferSize, int workGroupSize) {
    LOCAL mixed tempBuffer[512];
    const unsigned int thread = LOCAL_ID;
    mixed sum = 0;
    for (unsigned int index = GLOBAL_ID; index < bufferSize; index += GLOBAL_SIZE)
        sum += energyBuffer[index];
    tempBuffer[thread] = sum;
    for (int i = 1; i < workGroupSize; i *= 2) {
        SYNC_THREADS;
        if (thread%(i*2) == 0 && thread+i < workGroupSize)
            tempBuffer[thread] += tempBuffer[thread+i];
    }
    if (thread == 0)
        result[GROUP_ID] = tempBuffer[0];
}

/**
 * Record the atomic charges into the posq array.
 */
KERNEL void setCharges(GLOBAL real* RESTRICT charges, GLOBAL real4* RESTRICT posq, GLOBAL int* RESTRICT atomOrder, int numAtoms) {
    for (int i = GLOBAL_ID; i < numAtoms; i += GLOBAL_SIZE)
        posq[i].w = charges[atomOrder[i]];
}
