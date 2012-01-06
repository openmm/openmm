__kernel void computeSum(__global const float* restrict sumBuffer, __global float* result, unsigned int outputIndex) {
    __local float tempBuffer[WORK_GROUP_SIZE];
    const unsigned int thread = get_local_id(0);
    float sum = 0.0f;
    for (unsigned int index = thread; index < 3*NUM_ATOMS; index += get_local_size(0))
        sum += sumBuffer[index];
    tempBuffer[thread] = sum;
    for (int i = 1; i < WORK_GROUP_SIZE; i *= 2) {
        barrier(CLK_LOCAL_MEM_FENCE);
        if (thread%(i*2) == 0 && thread+i < WORK_GROUP_SIZE)
            tempBuffer[thread] += tempBuffer[thread+i];
    }
    if (thread == 0)
        result[outputIndex] = tempBuffer[0];
}
