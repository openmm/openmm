/**
 * Perform the first step of Brownian integration.
 */

__kernel void integrateBrownianPart1(float tauDeltaT, float noiseAmplitude, __global float4* force,
        __global float4* posDelta, __global float4* random, unsigned int randomIndex) {
    randomIndex += get_global_id(0);
    for (int index = get_global_id(0); index < NUM_ATOMS; index += get_global_size(0)) {
        posDelta[index] = (float4) (tauDeltaT*force[index].xyz + noiseAmplitude*random[randomIndex].xyz, 0.0f);
        randomIndex += get_global_size(0);
    }
}

/**
 * Perform the second step of Brownian integration.
 */

__kernel void integrateBrownianPart2(float oneOverDeltaT, __global float4* posq, __global float4* velm, __global float4* posDelta) {
    for (int index = get_global_id(0); index < NUM_ATOMS; index += get_global_size(0)) {
        float4 delta = posDelta[index];
        velm[index].xyz = oneOverDeltaT*delta.xyz;
        posq[index].xyz = posq[index].xyz + delta.xyz;
    }
}
