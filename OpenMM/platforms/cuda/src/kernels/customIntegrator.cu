extern "C" __global__ void computeFloatSum(const float* __restrict__ sumBuffer, float* result) {
    __shared__ float tempBuffer[WORK_GROUP_SIZE];
    const unsigned int thread = threadIdx.x;
    float sum = 0;
    for (unsigned int index = thread; index < SUM_BUFFER_SIZE; index += blockDim.x)
        sum += sumBuffer[index];
    tempBuffer[thread] = sum;
    for (int i = 1; i < WORK_GROUP_SIZE; i *= 2) {
        __syncthreads();
        if (thread%(i*2) == 0 && thread+i < WORK_GROUP_SIZE)
            tempBuffer[thread] += tempBuffer[thread+i];
    }
    if (thread == 0)
        *result = tempBuffer[0];
}

extern "C" __global__ void computeDoubleSum(const double* __restrict__ sumBuffer, double* result) {
    __shared__ double tempBuffer[WORK_GROUP_SIZE];
    const unsigned int thread = threadIdx.x;
    double sum = 0;
    for (unsigned int index = thread; index < SUM_BUFFER_SIZE; index += blockDim.x)
        sum += sumBuffer[index];
    tempBuffer[thread] = sum;
    for (int i = 1; i < WORK_GROUP_SIZE; i *= 2) {
        __syncthreads();
        if (thread%(i*2) == 0 && thread+i < WORK_GROUP_SIZE)
            tempBuffer[thread] += tempBuffer[thread+i];
    }
    if (thread == 0)
        *result = tempBuffer[0];
}

extern "C" __global__ void applyPositionDeltas(real4* __restrict__ posq, real4* __restrict__ posqCorrection, mixed4* __restrict__ posDelta) {
    for (unsigned int index = blockIdx.x*blockDim.x+threadIdx.x; index < NUM_ATOMS; index += blockDim.x*gridDim.x) {
#ifdef USE_MIXED_PRECISION
        real4 pos1 = posq[index];
        real4 pos2 = posqCorrection[index];
        mixed4 pos = make_mixed4(pos1.x+(mixed)pos2.x, pos1.y+(mixed)pos2.y, pos1.z+(mixed)pos2.z, pos1.w);
#else
        real4 pos = posq[index];
#endif
        pos.x += posDelta[index].x;
        pos.y += posDelta[index].y;
        pos.z += posDelta[index].z;
#ifdef USE_MIXED_PRECISION
        posq[index] = make_real4((real) pos.x, (real) pos.y, (real) pos.z, (real) pos.w);
        posqCorrection[index] = make_real4(pos.x-(real) pos.x, pos.y-(real) pos.y, pos.z-(real) pos.z, 0);
#else
        posq[index] = pos;
#endif
        posDelta[index] = make_mixed4(0, 0, 0, 0);
    }
}

extern "C" __global__ void generateRandomNumbers(int numValues, float4* __restrict__ random, uint4* __restrict__ seed) {
    uint4 state = seed[blockIdx.x*blockDim.x+threadIdx.x];
    unsigned int carry = 0;
    for (int index = blockIdx.x*blockDim.x+threadIdx.x; index < numValues; index += blockDim.x*gridDim.x) {
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

        random[index] = make_float4(x1, x2, x3, 0.0f);
    }
    seed[blockIdx.x*blockDim.x+threadIdx.x] = state;
}
