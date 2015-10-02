#define DIELECTRIC_OFFSET 0.009f
#define PROBE_RADIUS 0.14f

/**
 * Reduce the Born sums to compute the Born radii.
 */

__kernel void reduceBornSum(int bufferSize, int numBuffers, float alpha, float beta, float gamma,
#ifdef SUPPORTS_64_BIT_ATOMICS
            __global const long* restrict bornSum,
#else
            __global const real* restrict bornSum,
#endif
            __global const float2* restrict params, __global real* restrict bornRadii, __global real* restrict obcChain) {
    unsigned int index = get_global_id(0);
    while (index < NUM_ATOMS) {
        // Get summed Born data

        int totalSize = bufferSize*numBuffers;
#ifdef SUPPORTS_64_BIT_ATOMICS
        real sum = (1/(real) 0x100000000)*bornSum[index];
#else
        real sum = bornSum[index];
        for (int i = index+bufferSize; i < totalSize; i += bufferSize)
            sum += bornSum[i];
#endif

        // Now calculate Born radius and OBC term.

        float offsetRadius = params[index].x;
        sum *= 0.5f*offsetRadius;
        real sum2 = sum*sum;
        real sum3 = sum*sum2;
        real tanhSum = tanh(alpha*sum - beta*sum2 + gamma*sum3);
        real nonOffsetRadius = offsetRadius + DIELECTRIC_OFFSET;
        real radius = 1/(1/offsetRadius - tanhSum/nonOffsetRadius);
        real chain = offsetRadius*(alpha - 2*beta*sum + 3*gamma*sum2);
        chain = (1-tanhSum*tanhSum)*chain / nonOffsetRadius;
        bornRadii[index] = radius;
        obcChain[index] = chain;
        index += get_global_size(0);
    }
}

/**
 * Reduce the Born force.
 */

__kernel void reduceBornForce(int bufferSize, int numBuffers, __global real* bornForce,
#ifdef SUPPORTS_64_BIT_ATOMICS
            __global const long* restrict bornForceIn,
#endif
            __global mixed* restrict energyBuffer, __global const float2* restrict params, __global const real* restrict bornRadii, __global const real* restrict obcChain) {
    mixed energy = 0;
    unsigned int index = get_global_id(0);
    while (index < NUM_ATOMS) {
        // Sum the Born force

        int totalSize = bufferSize*numBuffers;
#ifdef SUPPORTS_64_BIT_ATOMICS
        real force = (1/(real) 0x100000000)*bornForceIn[index];
#else
        real force = bornForce[index];
        for (int i = index+bufferSize; i < totalSize; i += bufferSize)
            force += bornForce[i];
#endif
        // Now calculate the actual force

        float offsetRadius = params[index].x;
        real bornRadius = bornRadii[index];
        real r = offsetRadius+DIELECTRIC_OFFSET+PROBE_RADIUS;
        real ratio6 = pow((offsetRadius+DIELECTRIC_OFFSET)/bornRadius, (real) 6);
        real saTerm = SURFACE_AREA_FACTOR*r*r*ratio6;
        force += saTerm/bornRadius;
        energy += saTerm;
        force *= bornRadius*bornRadius*obcChain[index];
        bornForce[index] = force;
        index += get_global_size(0);
    }
    energyBuffer[get_global_id(0)] += energy/-6.0f;
}