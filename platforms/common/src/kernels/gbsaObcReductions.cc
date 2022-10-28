#define DIELECTRIC_OFFSET 0.009f
#define PROBE_RADIUS 0.14f

/**
 * Reduce the Born sums to compute the Born radii.
 */

KERNEL void reduceBornSum(float alpha, float beta, float gamma,
            GLOBAL const mm_long* RESTRICT bornSum,
            GLOBAL const float2* RESTRICT params, GLOBAL real* RESTRICT bornRadii, GLOBAL real* RESTRICT obcChain) {
    for (unsigned int index = GLOBAL_ID; index < NUM_ATOMS; index += GLOBAL_SIZE) {
        // Get summed Born data

        real sum = RECIP((real) 0x100000000)*bornSum[index];

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
            GLOBAL mm_long* RESTRICT bornForce,
            GLOBAL mixed* RESTRICT energyBuffer, GLOBAL const float2* RESTRICT params, GLOBAL const real* RESTRICT bornRadii, GLOBAL const real* RESTRICT obcChain) {
    mixed energy = 0;
    for (unsigned int index = GLOBAL_ID; index < NUM_ATOMS; index += GLOBAL_SIZE) {
        // Get summed Born force

        real force = RECIP((real) 0x100000000)*bornForce[index];

        // Now calculate the actual force

        float offsetRadius = params[index].x;
        real bornRadius = bornRadii[index];
        real r = offsetRadius+DIELECTRIC_OFFSET+PROBE_RADIUS;
        real ratio6 = POW((offsetRadius+DIELECTRIC_OFFSET)/bornRadius, (real) 6);
        real saTerm = SURFACE_AREA_FACTOR*r*r*ratio6;
        force += saTerm/bornRadius;
        energy += saTerm;
        force *= bornRadius*bornRadius*obcChain[index];
        bornForce[index] = realToFixedPoint(force);
    }
    energyBuffer[GLOBAL_ID] += energy/-6;
}
