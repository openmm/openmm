extern "C" {

/**
 * This is called by the various functions below to clear a buffer.
 */
__device__ void clearSingleBuffer(int* __restrict__ buffer, int size) {
    int index = blockDim.x*blockIdx.x+threadIdx.x;
    int4* buffer4 = (int4*) buffer;
    int sizeDiv4 = size/4;
    while (index < sizeDiv4) {
        buffer4[index] = make_int4(0);
        index += blockDim.x*gridDim.x;
    }
    if (blockDim.x*blockIdx.x+threadIdx.x == 0)
        for (int i = sizeDiv4*4; i < size; i++)
            buffer[i] = 0;
}

/**
 * Fill a buffer with 0.
 */
__global__ void clearBuffer(int* __restrict__ buffer, int size) {
    clearSingleBuffer(buffer, size);
}

/**
 * Fill two buffers with 0.
 */
__global__ void clearTwoBuffers(int* __restrict__ buffer1, int size1, int* __restrict__ buffer2, int size2) {
    clearSingleBuffer(buffer1, size1);
    clearSingleBuffer(buffer2, size2);
}

/**
 * Fill three buffers with 0.
 */
__global__ void clearThreeBuffers(int* __restrict__ buffer1, int size1, int* __restrict__ buffer2, int size2, int* __restrict__ buffer3, int size3) {
    clearSingleBuffer(buffer1, size1);
    clearSingleBuffer(buffer2, size2);
    clearSingleBuffer(buffer3, size3);
}

/**
 * Fill four buffers with 0.
 */
__global__ void clearFourBuffers(int* __restrict__ buffer1, int size1, int* __restrict__ buffer2, int size2, int* __restrict__ buffer3, int size3, int* __restrict__ buffer4, int size4) {
    clearSingleBuffer(buffer1, size1);
    clearSingleBuffer(buffer2, size2);
    clearSingleBuffer(buffer3, size3);
    clearSingleBuffer(buffer4, size4);
}

/**
 * Fill five buffers with 0.
 */
__global__ void clearFiveBuffers(int* __restrict__ buffer1, int size1, int* __restrict__ buffer2, int size2, int* __restrict__ buffer3, int size3, int* __restrict__ buffer4, int size4, int* __restrict__ buffer5, int size5) {
    clearSingleBuffer(buffer1, size1);
    clearSingleBuffer(buffer2, size2);
    clearSingleBuffer(buffer3, size3);
    clearSingleBuffer(buffer4, size4);
    clearSingleBuffer(buffer5, size5);
}

/**
 * Fill six buffers with 0.
 */
__global__ void clearSixBuffers(int* __restrict__ buffer1, int size1, int* __restrict__ buffer2, int size2, int* __restrict__ buffer3, int size3, int* __restrict__ buffer4, int size4, int* __restrict__ buffer5, int size5, int* __restrict__ buffer6, int size6) {
    clearSingleBuffer(buffer1, size1);
    clearSingleBuffer(buffer2, size2);
    clearSingleBuffer(buffer3, size3);
    clearSingleBuffer(buffer4, size4);
    clearSingleBuffer(buffer5, size5);
    clearSingleBuffer(buffer6, size6);
}

}