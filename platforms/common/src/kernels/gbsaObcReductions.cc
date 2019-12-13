#define DIELECTRIC_OFFSET 0.009f
#define PROBE_RADIUS 0.14f

/**
 * Reduce the Born sums to compute the Born radii.
 */

KERNEL void reduceBornSum(float alpha, float beta, float gamma,
#ifdef SUPPORTS_64_BIT_ATOMICS
            GLOBAL const mm_long* RESTRICT bornSum,
#else
            GLOBAL const real* RESTRICT bornSum, int bufferSize, int numBuffers,
#endif
            GLOBAL const float2* RESTRICT params, GLOBAL real* RESTRICT bornRadii, GLOBAL real* RESTRICT obcChain) {
    for (unsigned int index = GLOBAL_ID; index < NUM_ATOMS; index += GLOBAL_SIZE) {
        // Get summed Born data

#ifdef SUPPORTS_64_BIT_ATOMICS
        real sum = RECIP((real) 0x100000000)*bornSum[index];
#else
        real sum = bornSum[index];
        int totalSize = bufferSize*numBuffers;
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
        real radius = RECIP(RECIP(offsetRadius) - tanhSum/nonOffsetRadius);
        real chain = offsetRadius*(alpha - 2*beta*sum + 3*gamma*sum2);
        chain = (1-tanhSum*tanhSum)*chain / nonOffsetRadius;
        bornRadii[index] = radius;
        obcChain[index] = chain;
    }
}

/**
 * Reduce the Born force.
 */

KERNEL void reduceBornForce(
#ifdef SUPPORTS_64_BIT_ATOMICS
            GLOBAL mm_long* RESTRICT bornForce,
#else
            GLOBAL real* bornForce, int bufferSize, int numBuffers,
#endif
            GLOBAL mixed* RESTRICT energyBuffer, GLOBAL const float2* RESTRICT params, GLOBAL const real* RESTRICT bornRadii, GLOBAL const real* RESTRICT obcChain) {
    mixed energy = 0;
    for (unsigned int index = GLOBAL_ID; index < NUM_ATOMS; index += GLOBAL_SIZE) {
        // Get summed Born force

#ifdef SUPPORTS_64_BIT_ATOMICS
        real force = RECIP((real) 0x100000000)*bornForce[index];
#else
        real force = bornForce[index];
        int totalSize = bufferSize*numBuffers;
        for (int i = index+bufferSize; i < totalSize; i += bufferSize)
            force += bornForce[i];
#endif
        // Now calculate the actual force

        float offsetRadius = params[index].x;
        real bornRadius = bornRadii[index];
        real r = offsetRadius+DIELECTRIC_OFFSET+PROBE_RADIUS;
        real ratio6 = POW((offsetRadius+DIELECTRIC_OFFSET)/bornRadius, (real) 6);
        real saTerm = SURFACE_AREA_FACTOR*r*r*ratio6;
        force += saTerm/bornRadius;
        energy += saTerm;
        force *= bornRadius*bornRadius*obcChain[index];
#ifdef SUPPORTS_64_BIT_ATOMICS
        bornForce[index] = (mm_long) (force*0x100000000);
#else
        bornForce[index] = force;
#endif
    }
    energyBuffer[GLOBAL_ID] += energy/-6;
}