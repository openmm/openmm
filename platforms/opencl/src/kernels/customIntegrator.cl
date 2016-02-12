__kernel void computeFloatSum(__global const float* restrict sumBuffer, __global float* result, int bufferSize) {
    __local float tempBuffer[WORK_GROUP_SIZE];
    const unsigned int thread = get_local_id(0);
    float sum = 0;
    for (unsigned int index = thread; index < bufferSize; index += get_local_size(0))
        sum += sumBuffer[index];
    tempBuffer[thread] = sum;
    for (int i = 1; i < WORK_GROUP_SIZE; i *= 2) {
        barrier(CLK_LOCAL_MEM_FENCE);
        if (thread%(i*2) == 0 && thread+i < WORK_GROUP_SIZE)
            tempBuffer[thread] += tempBuffer[thread+i];
    }
    if (thread == 0)
        *result = tempBuffer[0];
}

#ifdef SUPPORTS_DOUBLE_PRECISION
__kernel void computeDoubleSum(__global const double* restrict sumBuffer, __global double* result, int bufferSize) {
    __local double tempBuffer[WORK_GROUP_SIZE];
    const unsigned int thread = get_local_id(0);
    double sum = 0;
    for (unsigned int index = thread; index < bufferSize; index += get_local_size(0))
        sum += sumBuffer[index];
    tempBuffer[thread] = sum;
    for (int i = 1; i < WORK_GROUP_SIZE; i *= 2) {
        barrier(CLK_LOCAL_MEM_FENCE);
        if (thread%(i*2) == 0 && thread+i < WORK_GROUP_SIZE)
            tempBuffer[thread] += tempBuffer[thread+i];
    }
    if (thread == 0)
        *result = tempBuffer[0];
}
#endif

__kernel void applyPositionDeltas(__global real4* restrict posq, __global real4* restrict posqCorrection, __global mixed4* restrict posDelta) {
    for (unsigned int index = get_global_id(0); index < NUM_ATOMS; index += get_global_size(0)) {
#ifdef USE_MIXED_PRECISION
        real4 pos1 = posq[index];
        real4 pos2 = posqCorrection[index];
        mixed4 pos = (mixed4) (pos1.x+(mixed)pos2.x, pos1.y+(mixed)pos2.y, pos1.z+(mixed)pos2.z, pos1.w);
#else
        real4 pos = posq[index];
#endif
        pos.xyz += posDelta[index].xyz;
#ifdef USE_MIXED_PRECISION
        posq[index] = (real4) ((real) pos.x, (real) pos.y, (real) pos.z, (real) pos.w);
        posqCorrection[index] = (real4) (pos.x-(real) pos.x, pos.y-(real) pos.y, pos.z-(real) pos.z, 0);
#else
        posq[index] = pos;
#endif
        posDelta[index] = (mixed4) 0;
    }
}

__kernel void generateRandomNumbers(int numValues, __global float4* restrict random, __global uint4* restrict seed) {
    uint4 state = seed[get_global_id(0)];
    unsigned int carry = 0;
    for (int index = get_global_id(0); index < numValues; index += get_global_size(0)) {
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
