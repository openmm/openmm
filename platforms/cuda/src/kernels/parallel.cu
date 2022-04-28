/**
 * Sum the forces computed by different contexts.
 */

extern "C" __global__ void sumForces(long long* __restrict__ force, long long* __restrict__ buffer, int bufferSize, int numBuffers) {
    long long totalSize = bufferSize*numBuffers;
    for (int index = blockDim.x*blockIdx.x+threadIdx.x; index < bufferSize; index += blockDim.x*gridDim.x) {
        long long sum = force[index];
        for (long long i = index; i < totalSize; i += bufferSize)
            sum += buffer[i];
        force[index] = sum;
    }
}
