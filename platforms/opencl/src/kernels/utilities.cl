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
