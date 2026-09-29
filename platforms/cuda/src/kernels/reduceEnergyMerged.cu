/**
 * Mixed-precision-only experiment: combine the existing PME elementwise
 * addition with the unchanged per-block energy reduction schedule.
 *
 * Materialize each merged value back to energyBuffer, preserving its contents.
 * Explicit round-to-nearest additions prevent reassociation across the original
 * addEnergy materialization and reduceEnergy accumulation boundary.
 */
extern "C" __global__ void reduceEnergyMerged(double* __restrict__ energyBuffer,
        const double* __restrict__ pmeEnergyBuffer, double* __restrict__ result,
        int bufferSize, int pmeSize, int workGroupSize) {
    extern __shared__ double tempBuffer[];
    const unsigned int thread = threadIdx.x;
    double sum = 0;
    for (unsigned int index = blockDim.x*blockIdx.x+threadIdx.x; index < bufferSize; index += blockDim.x*gridDim.x) {
        double value = energyBuffer[index];
        if (index < pmeSize) {
            value = __dadd_rn(value, pmeEnergyBuffer[index]);
            energyBuffer[index] = value;
        }
        sum = __dadd_rn(sum, value);
    }
    tempBuffer[thread] = sum;
    for (int i = 1; i < workGroupSize; i *= 2) {
        __syncthreads();
        if (thread%(i*2) == 0 && thread+i < workGroupSize)
            tempBuffer[thread] += tempBuffer[thread+i];
    }
    if (thread == 0)
        result[blockIdx.x] = tempBuffer[0];
}
