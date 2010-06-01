/**
 * Fill a buffer with 0.
 */

__kernel void clearBuffer(__global float* buffer, int size) {
    int index = get_global_id(0);
    __global float4* buffer4 = (__global float4*) buffer;
    int sizeDiv4 = size/4;
    while (index < sizeDiv4) {
        buffer4[index] = (float4) (0.0f);
        index += get_global_size(0);
    }
    if (get_global_id(0) == 0)
        for (int i = sizeDiv4*4; i < size; i++)
            buffer[i] = 0.0f;
}

/**
 * Fill two buffers with 0.
 */
__kernel void clearTwoBuffers(__global float* buffer1, int size1, __global float* buffer2, int size2) {
    clearBuffer(buffer1, size1);
    clearBuffer(buffer2, size2);
}

/**
 * Fill three buffers with 0.
 */
__kernel void clearThreeBuffers(__global float* buffer1, int size1, __global float* buffer2, int size2, __global float* buffer3, int size3) {
    clearBuffer(buffer1, size1);
    clearBuffer(buffer2, size2);
    clearBuffer(buffer3, size3);
}

/**
 * Fill four buffers with 0.
 */
__kernel void clearFourBuffers(__global float* buffer1, int size1, __global float* buffer2, int size2, __global float* buffer3, int size3, __global float* buffer4, int size4) {
    clearBuffer(buffer1, size1);
    clearBuffer(buffer2, size2);
    clearBuffer(buffer3, size3);
    clearBuffer(buffer4, size4);
}

/**
 * Sum a collection of buffers into the first one.
 */

__kernel void reduceFloat4Buffer(__global float4* buffer, int bufferSize, int numBuffers) {
    int index = get_global_id(0);
    int totalSize = bufferSize*numBuffers;
    while (index < bufferSize) {
        float4 sum = buffer[index];
        for (int i = index+bufferSize; i < totalSize; i += bufferSize)
            sum += buffer[i];
        buffer[index] = sum;
        index += get_global_size(0);
    }
}

/**
 * This is called to determine the accuracy of native_sqrt(), native_rsqrt() and native_recip().
 */

__kernel void determineNativeAccuracy(__global float4* values, int numValues) {
    for (int i = 0; i < numValues; ++i) {
        float v = values[i].x;
        values[i] = (float4) (v, native_sqrt(v), native_rsqrt(v), native_recip(v));
    }
}