__kernel void computePerDof(__global float4* restrict posq, __global float4* restrict posDelta, __global float4* restrict velm,
        __global const float4* restrict force, __global const float2* restrict dt, __global const float* restrict globals,
        __global const float* restrict params, __global float* restrict sum, __global const float4* restrict random,
        unsigned int randomIndex, float energy
        PARAMETER_ARGUMENTS) {
    float stepSize = dt[0].y;
    int index = get_global_id(0);
    randomIndex += index;
    while (index < NUM_ATOMS) {
        float4 position = posq[index];
        float4 velocity = velm[index];
        float4 f = force[index];
        float4 gaussian = random[randomIndex];
        float mass = 1.0f/velocity.w;
        COMPUTE_STEP
        randomIndex += get_global_size(0);
        index += get_global_size(0);
    }
}
