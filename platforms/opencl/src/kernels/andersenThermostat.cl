/**
 * Apply the Andersen thermostat to adjust particle velocities.
 */

__kernel void applyAndersenThermostat(float collisionFrequency, float kT, __global float4* velm, __global const float2* restrict stepSize, __global const float4* restrict random,
        unsigned int randomIndex, __global const int* restrict atomGroups) {
    float collisionProbability = 1.0f-exp(-collisionFrequency*stepSize[0].y);
    float randomRange = erf(collisionProbability/sqrt(2.0f));
    for (int index = get_global_id(0); index < NUM_ATOMS; index += get_global_size(0)) {
        float4 velocity = velm[index];
        float4 selectRand = random[randomIndex+atomGroups[index]];
        float4 velRand = random[randomIndex+index];
        float scale = (selectRand.w > -randomRange && selectRand.w < randomRange ? 0.0f : 1.0f);
        float add = (1.0f-scale)*sqrt(kT*velocity.w);
        velocity.xyz = scale*velocity.xyz + add*velRand.xyz;
        velm[index] = velocity;
    }
}
