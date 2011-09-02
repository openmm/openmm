#define DIELECTRIC_OFFSET 0.009f
#define PROBE_RADIUS 0.14f
#define SURFACE_AREA_FACTOR -170.351730667551f //-6.0f*3.14159265358979323846f*0.0216f*1000.0f*0.4184f;

/**
 * Reduce the Born sums to compute the Born radii.
 */

__kernel void reduceBornSum(int bufferSize, int numBuffers, float alpha, float beta, float gamma,
#ifdef SUPPORTS_64_BIT_ATOMICS
            __global long* bornSum,
#else
            __global float* bornSum,
#endif
            __global float2* params, __global float* bornRadii, __global float* obcChain) {
    unsigned int index = get_global_id(0);
    while (index < NUM_ATOMS) {
        // Get summed Born data

        int totalSize = bufferSize*numBuffers;
#ifdef SUPPORTS_64_BIT_ATOMICS
        float sum = (1.0f/0xFFFFFFFF)*bornSum[index];
#else
        float sum = bornSum[index];
        for (int i = index+bufferSize; i < totalSize; i += bufferSize)
            sum += bornSum[i];
#endif

        // Now calculate Born radius and OBC term.

        float offsetRadius = params[index].x;
        sum *= 0.5f*offsetRadius;
        float sum2 = sum*sum;
        float sum3 = sum*sum2;
        float tanhSum = tanh(alpha*sum - beta*sum2 + gamma*sum3);
        float nonOffsetRadius = offsetRadius + DIELECTRIC_OFFSET;
        float radius = 1.0f/(1.0f/offsetRadius - tanhSum/nonOffsetRadius);
        float chain = offsetRadius*(alpha - 2.0f*beta*sum + 3.0f*gamma*sum2);
        chain = (1.0f-tanhSum*tanhSum)*chain / nonOffsetRadius;
        bornRadii[index] = radius;
        obcChain[index] = chain;
        index += get_global_size(0);
    }
}

/**
 * Reduce the Born force.
 */

__kernel void reduceBornForce(int bufferSize, int numBuffers, __global float* bornForce,
#ifdef SUPPORTS_64_BIT_ATOMICS
            __global long* bornForceIn,
#endif
            __global float* energyBuffer, __global float2* params, __global float* bornRadii, __global float* obcChain) {
    float energy = 0.0f;
    unsigned int index = get_global_id(0);
    while (index < NUM_ATOMS) {
        // Sum the Born force

        int totalSize = bufferSize*numBuffers;
#ifdef SUPPORTS_64_BIT_ATOMICS
        float force = (1.0f/0xFFFFFFFF)*bornForceIn[index];
#else
        float force = bornForce[index];
        for (int i = index+bufferSize; i < totalSize; i += bufferSize)
            force += bornForce[i];
#endif
        // Now calculate the actual force

        float offsetRadius = params[index].x;
        float bornRadius = bornRadii[index];
        float r = offsetRadius+DIELECTRIC_OFFSET+PROBE_RADIUS;
        float ratio6 = pow((offsetRadius+DIELECTRIC_OFFSET)/bornRadius, 6.0f);
        float saTerm = SURFACE_AREA_FACTOR*r*r*ratio6;
        force += saTerm/bornRadius;
        energy += saTerm;
        force *= bornRadius*bornRadius*obcChain[index];
        bornForce[index] = force;
        index += get_global_size(0);
    }
    energyBuffer[get_global_id(0)] += energy/-6.0f;
}