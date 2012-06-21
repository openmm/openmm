#define DIELECTRIC_OFFSET 0.009f
#define PROBE_RADIUS 0.14f
#define SURFACE_AREA_FACTOR -170.351730667551f //-6.0f*3.14159265358979323846f*0.0216f*1000.0f*0.4184f;

/**
 * Reduce the Born sums to compute the Born radii.
 */

extern "C" __global__ void reduceBornSum(float alpha, float beta, float gamma, const long long* __restrict__ bornSum,
            const float2* __restrict__ params, real* __restrict__ bornRadii, real* __restrict__ obcChain) {
    for (unsigned int index = blockIdx.x*blockDim.x+threadIdx.x; index < NUM_ATOMS; index += blockDim.x*gridDim.x) {
        // Get summed Born data

        real sum = RECIP(0xFFFFFFFF)*bornSum[index];

        // Now calculate Born radius and OBC term.

        float offsetRadius = params[index].x;
        sum *= 0.5f*offsetRadius;
        real sum2 = sum*sum;
        real sum3 = sum*sum2;
        real tanhSum = tanh(alpha*sum - beta*sum2 + gamma*sum3);
        real nonOffsetRadius = offsetRadius + DIELECTRIC_OFFSET;
        real radius = RECIP(RECIP(offsetRadius) - tanhSum/nonOffsetRadius);
        real chain = offsetRadius*(alpha - 2.0f*beta*sum + 3.0f*gamma*sum2);
        chain = (1-tanhSum*tanhSum)*chain / nonOffsetRadius;
        bornRadii[index] = radius;
        obcChain[index] = chain;
    }
}

/**
 * Reduce the Born force.
 */

extern "C" __global__ void reduceBornForce(long long* __restrict__ bornForce, real* __restrict__ energyBuffer,
        const float2* __restrict__ params, const real* __restrict__ bornRadii, const real* __restrict__ obcChain) {
    real energy = 0;
    for (unsigned int index = blockIdx.x*blockDim.x+threadIdx.x; index < NUM_ATOMS; index += blockDim.x*gridDim.x) {
        // Get summed Born force

        real force = RECIP(0xFFFFFFFF)*bornForce[index];

        // Now calculate the actual force

        float offsetRadius = params[index].x;
        real bornRadius = bornRadii[index];
        real r = offsetRadius+DIELECTRIC_OFFSET+PROBE_RADIUS;
        real ratio6 = POW((offsetRadius+DIELECTRIC_OFFSET)/bornRadius, 6);
        real saTerm = SURFACE_AREA_FACTOR*r*r*ratio6;
        force += saTerm/bornRadius;
        energy += saTerm;
        force *= bornRadius*bornRadius*obcChain[index];
        bornForce[index] = (long long) (force*0xFFFFFFFF);
    }
    energyBuffer[blockIdx.x*blockDim.x+threadIdx.x] += energy/-6;
}