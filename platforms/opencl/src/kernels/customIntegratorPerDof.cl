#ifdef SUPPORTS_DOUBLE_PRECISION
#pragma OPENCL EXTENSION cl_khr_fp64 : enable
#endif

__kernel void computePerDof(__global float4* restrict posq, __global float4* restrict posDelta, __global float4* restrict velm,
        __global const float4* restrict force, __global const float2* restrict dt, __global const float* restrict globals,
        __global const float* restrict params, __global float* restrict sum, __global const float4* restrict gaussianValues,
        unsigned int randomIndex, __global const float4* restrict uniformValues, __global const float* restrict energy
        PARAMETER_ARGUMENTS) {
    float stepSize = dt[0].y;
    int index = get_global_id(0);
    randomIndex += index;
    while (index < NUM_ATOMS) {
#ifdef SUPPORTS_DOUBLE_PRECISION
#ifdef LOAD_POS_AS_DELTA
        double4 position = convert_double4(posq[index]+posDelta[index]);
#else
        double4 position = convert_double4(posq[index]);
#endif
        double4 velocity = convert_double4(velm[index]);
        double4 f = convert_double4(force[index]);
        double mass = 1.0/velocity.w;
#else
#ifdef LOAD_POS_AS_DELTA
        float4 position = posq[index]+posDelta[index];
#else
        float4 position = posq[index];
#endif
        float4 velocity = velm[index];
        float4 f = force[index];
        float mass = 1.0f/velocity.w;
#endif
        if (velocity.w != 0.0) {
            float4 gaussian = gaussianValues[randomIndex];
            float4 uniform = uniformValues[index];
            COMPUTE_STEP
        }
        randomIndex += get_global_size(0);
        index += get_global_size(0);
    }
}
