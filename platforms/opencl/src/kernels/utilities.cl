/**
 * Fill a buffer with 0.
 */

__kernel void clearBuffer(__global int* restrict buffer, int size) {
    int index = get_global_id(0);
    __global int4* buffer4 = (__global int4*) buffer;
    int sizeDiv4 = size/4;
    while (index < sizeDiv4) {
        buffer4[index] = (int4) 0;
        index += get_global_size(0);
    }
    if (get_global_id(0) == 0)
        for (int i = sizeDiv4*4; i < size; i++)
            buffer[i] = 0;
}

/**
 * Fill two buffers with 0.
 */
__kernel void clearTwoBuffers(__global int* restrict buffer1, int size1, __global int* restrict buffer2, int size2) {
    clearBuffer(buffer1, size1);
    clearBuffer(buffer2, size2);
}

/**
 * Fill three buffers with 0.
 */
__kernel void clearThreeBuffers(__global int* restrict buffer1, int size1, __global int* restrict buffer2, int size2, __global int* restrict buffer3, int size3) {
    clearBuffer(buffer1, size1);
    clearBuffer(buffer2, size2);
    clearBuffer(buffer3, size3);
}

/**
 * Fill four buffers with 0.
 */
__kernel void clearFourBuffers(__global int* restrict buffer1, int size1, __global int* restrict buffer2, int size2, __global int* restrict buffer3, int size3, __global int* restrict buffer4, int size4) {
    clearBuffer(buffer1, size1);
    clearBuffer(buffer2, size2);
    clearBuffer(buffer3, size3);
    clearBuffer(buffer4, size4);
}

/**
 * Fill five buffers with 0.
 */
__kernel void clearFiveBuffers(__global int* restrict buffer1, int size1, __global int* restrict buffer2, int size2, __global int* restrict buffer3, int size3, __global int* restrict buffer4, int size4, __global int* restrict buffer5, int size5) {
    clearBuffer(buffer1, size1);
    clearBuffer(buffer2, size2);
    clearBuffer(buffer3, size3);
    clearBuffer(buffer4, size4);
    clearBuffer(buffer5, size5);
}

/**
 * Fill six buffers with 0.
 */
__kernel void clearSixBuffers(__global int* restrict buffer1, int size1, __global int* restrict buffer2, int size2, __global int* restrict buffer3, int size3, __global int* restrict buffer4, int size4, __global int* restrict buffer5, int size5, __global int* restrict buffer6, int size6) {
    clearBuffer(buffer1, size1);
    clearBuffer(buffer2, size2);
    clearBuffer(buffer3, size3);
    clearBuffer(buffer4, size4);
    clearBuffer(buffer5, size5);
    clearBuffer(buffer6, size6);
}

/**
 * Sum a collection of buffers into the first one.
 */

__kernel void reduceReal4Buffer(__global real4* restrict buffer, int bufferSize, int numBuffers) {
    int index = get_global_id(0);
    int totalSize = bufferSize*numBuffers;
    while (index < bufferSize) {
        real4 sum = buffer[index];
        for (int i = index+bufferSize; i < totalSize; i += bufferSize)
            sum += buffer[i];
        buffer[index] = sum;
        index += get_global_size(0);
    }
}

#ifdef SUPPORTS_64_BIT_ATOMICS
/**
 * Sum the various buffers containing forces.
 */
__kernel void reduceForces(__global const long* restrict longBuffer, __global real4* restrict buffer, int bufferSize, int numBuffers) {
    int totalSize = bufferSize*numBuffers;
    real scale = 1/(real) 0x100000000;
    for (int index = get_global_id(0); index < bufferSize; index += get_global_size(0)) {
        real4 sum = (real4) (scale*longBuffer[index], scale*longBuffer[index+bufferSize], scale*longBuffer[index+2*bufferSize], 0);
        for (int i = index; i < totalSize; i += bufferSize)
            sum += buffer[i];
        buffer[index] = sum;
    }
}
#endif

/**
 * Sum the energy buffer.
 */
__kernel void reduceEnergy(__global const mixed* restrict energyBuffer, __global mixed* restrict result, int bufferSize, int workGroupSize, __local mixed* tempBuffer) {
    const unsigned int thread = get_local_id(0);
    mixed sum = 0;
    for (unsigned int index = thread; index < bufferSize; index += get_local_size(0))
        sum += energyBuffer[index];
    tempBuffer[thread] = sum;
    for (int i = 1; i < workGroupSize; i *= 2) {
        barrier(CLK_LOCAL_MEM_FENCE);
        if (thread%(i*2) == 0 && thread+i < workGroupSize)
            tempBuffer[thread] += tempBuffer[thread+i];
    }
    if (thread == 0)
        *result = tempBuffer[0];
}

/**
 * This is called to determine the accuracy of various native functions.
 */

__kernel void determineNativeAccuracy(__global float8* restrict values, int numValues) {
    for (int i = get_global_id(0); i < numValues; i += get_global_size(0)) {
        float v = values[i].s0;
        values[i] = (float8) (v, native_sqrt(v), native_rsqrt(v), native_recip(v), native_exp(v), native_log(v), 0.0f, 0.0f);
    }
}

/**
 * Record the atomic charges into the posq array.
 */
__kernel void setCharges(__global real* restrict charges, __global real4* restrict posq, __global int* restrict atomOrder, int numAtoms) {
    for (int i = get_global_id(0); i < numAtoms; i += get_global_size(0))
        posq[i].w = charges[atomOrder[i]];
}