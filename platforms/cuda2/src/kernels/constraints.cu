extern "C" __global__ void applyPositionDeltas(real4* __restrict__ posq, real4* __restrict__ posDelta) {
    for (unsigned int index = blockIdx.x*blockDim.x+threadIdx.x; index < NUM_ATOMS; index += blockDim.x*gridDim.x) {
        real4 position = posq[index];
        position.x += posDelta[index].x;
        position.y += posDelta[index].y;
        position.z += posDelta[index].z;
        posq[index] = position;
    }
}
