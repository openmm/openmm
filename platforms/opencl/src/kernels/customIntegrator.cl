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

__kernel void applyPositionDeltas(__global float4* restrict posq, __global float4* restrict posDelta) {
    for (unsigned int index = get_global_id(0); index < NUM_ATOMS; index += get_global_size(0)) {
        float4 position = posq[index];
        position.xyz += posDelta[index].xyz;
        posq[index] = position;
        posDelta[index] = (float4) 0.0f;
    }
}

__kernel void generateRandomNumbers(__global float4* restrict random, __global uint4* restrict seed) {
    uint4 state = seed[get_global_id(0)];
    unsigned int carry = 0;
    for (int index = get_global_id(0); index < NUM_ATOMS; index += get_global_size(0)) {
        // Generate three uniform random numbers.

        state.x = state.x * 69069 + 1;
        state.y ^= state.y << 13;
        state.y ^= state.y >> 17;
        state.y ^= state.y << 5;
        unsigned int k = (state.z >> 2) + (state.w >> 3) + (carry >> 2);
        unsigned int m = state.w + state.w + state.z + carry;
        state.z = state.w;
        state.w = m;
        carry = k >> 30;
        float x1 = (float)max(state.x + state.y + state.w, 0x00000001u) / (float)0xffffffff;
        state.x = state.x * 69069 + 1;
        state.y ^= state.y << 13;
        state.y ^= state.y >> 17;
        state.y ^= state.y << 5;
        k = (state.z >> 2) + (state.w >> 3) + (carry >> 2);
        m = state.w + state.w + state.z + carry;
        state.z = state.w;
        state.w = m;
        carry = k >> 30;
        float x2 = (float)max(state.x + state.y + state.w, 0x00000001u) / (float)0xffffffff;
        state.x = state.x * 69069 + 1;
        state.y ^= state.y << 13;
        state.y ^= state.y >> 17;
        state.y ^= state.y << 5;
        k = (state.z >> 2) + (state.w >> 3) + (carry >> 2);
        m = state.w + state.w + state.z + carry;
        state.z = state.w;
        state.w = m;
        carry = k >> 30;
        float x3 = (float)max(state.x + state.y + state.w, 0x00000001u) / (float)0xffffffff;

        // Record the values.

        random[index] = (float4) (x1, x2, x3, 0.0f);
    }
    seed[get_global_id(0)] = state;
}
