/**
 * Apply the Andersen thermostat to adjust particle velocities.
 */

extern "C" __global__ void applyAndersenThermostat(int numAtoms, float collisionFrequency, float kT, mixed4* velm, const mixed4* __restrict__ stepSize, const float4* __restrict__ random,
        unsigned int randomIndex, const int* __restrict__ atomGroups) {
    float collisionProbability = 1.0f-expf(-(float) (collisionFrequency*stepSize[0].y));
    float randomRange = erff(collisionProbability/sqrtf(2.0f));
    for (int index = blockIdx.x*blockDim.x+threadIdx.x; index < numAtoms; index += blockDim.x*gridDim.x) {
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
