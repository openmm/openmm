/**
 * Apply the Andersen thermostat to adjust particle velocities.
 */

__kernel void applyAndersenThermostat(float collisionFrequency, float kT, __global mixed4* velm, __global const mixed2* restrict stepSize, __global const float4* restrict random,
        unsigned int randomIndex, __global const int* restrict atomGroups) {
    float collisionProbability = 1.0f-EXP(-collisionFrequency*stepSize[0].y);
    float randomRange = erf(collisionProbability/SQRT(2.0f));
    for (int index = get_global_id(0); index < NUM_ATOMS; index += get_global_size(0)) {
        mixed4 velocity = velm[index];
        float4 selectRand = random[randomIndex+atomGroups[index]];
        float4 velRand = random[randomIndex+index];
        real scale = (selectRand.w > -randomRange && selectRand.w < randomRange ? 0 : 1);
        real add = (1-scale)*SQRT(kT*velocity.w);
        velocity.x = scale*velocity.x + add*velRand.x;
        velocity.y = scale*velocity.y + add*velRand.y;
        velocity.z = scale*velocity.z + add*velRand.z;
        velm[index] = velocity;
    }
}
