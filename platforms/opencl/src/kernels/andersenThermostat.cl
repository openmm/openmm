/**
 * Apply the Andersen thermostat to adjust particle velocities.
 */

__kernel void applyAndersenThermostat(float collisionFrequency, float kT, __global float4* velm, __global float2* stepSize, __global float4* random, unsigned int randomIndex) {
    randomIndex += get_global_id(0);
    float collisionProbability = 1.0f-exp(-collisionFrequency*stepSize[0].y);
    float randomRange = erf(collisionProbability/sqrt(2.0f));
    for (int index = get_global_id(0); index < NUM_ATOMS; index += get_global_size(0)) {
        float4 velocity = velm[index];
        float4 rand = random[randomIndex];
        float scale = (rand.w > -randomRange && rand.w < randomRange ? 0.0f : 1.0f);
        float add = (1.0f-scale)*sqrt(kT*velocity.w);
        velocity.xyz = scale*velocity.xyz + add*rand.xyz;
        velm[index] = velocity;
        randomIndex += get_global_size(0);
    }
}
