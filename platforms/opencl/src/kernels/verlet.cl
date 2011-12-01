#ifdef SUPPORTS_DOUBLE_PRECISION
#pragma OPENCL EXTENSION cl_khr_fp64 : enable
#endif

/**
 * Perform the first step of verlet integration.
 */

__kernel void integrateVerletPart1(int numAtoms, __global const float2* restrict dt, __global const float4* restrict posq, __global float4* restrict velm, __global const float4* restrict force, __global float4* restrict posDelta) {
    float2 stepSize = dt[0];
    float dtPos = stepSize.y;
    float dtVel = 0.5f*(stepSize.x+stepSize.y);
    int index = get_global_id(0);
    while (index < numAtoms) {
        float4 pos = posq[index];
        float4 velocity = velm[index];
        velocity.xyz += force[index].xyz*dtVel*velocity.w;
        pos.xyz = velocity.xyz*dtPos;
        posDelta[index] = pos;
        velm[index] = velocity;
        index += get_global_size(0);
    }
}

/**
 * Perform the second step of verlet integration.
 */

__kernel void integrateVerletPart2(int numAtoms, __global float2* restrict dt, __global float4* restrict posq, __global float4* restrict velm, __global const float4* restrict posDelta) {
    float2 stepSize = dt[0];
#ifdef SUPPORTS_DOUBLE_PRECISION
    double oneOverDt = 1.0/stepSize.y;
#else
    float oneOverDt = 1.0f/stepSize.y;
#endif
    if (get_global_id(0) == 0)
        dt[0].x = stepSize.y;
    barrier(CLK_LOCAL_MEM_FENCE);
    int index = get_global_id(0);
    while (index < numAtoms) {
        float4 pos = posq[index];
        float4 delta = posDelta[index];
        float4 velocity = velm[index];
        pos.xyz += delta.xyz;
#ifdef SUPPORTS_DOUBLE_PRECISION
        velocity.xyz = convert_float4(convert_double4(delta)*oneOverDt).xyz;
#else
        velocity.xyz = delta.xyz*oneOverDt;
#endif
        posq[index] = pos;
        velm[index] = velocity;
        index += get_global_size(0);
    }
}

/**
 * Select the step size to use for the next step.
 */

__kernel void selectVerletStepSize(int numAtoms, float maxStepSize, float errorTol, __global float2* restrict dt, __global const float4* restrict velm, __global const float4* restrict force, __local float* restrict error) {
    // Calculate the error.

    float err = 0.0f;
    int index = get_local_id(0);
    while (index < numAtoms) {
        float4 f = force[index];
        float invMass = velm[index].w;
        err += (f.x*f.x + f.y*f.y + f.z*f.z)*invMass;
        index += get_global_size(0);
    }
    error[get_local_id(0)] = err;
    barrier(CLK_LOCAL_MEM_FENCE);

    // Sum the errors from all threads.

    for (unsigned int offset = 1; offset < get_local_size(0); offset *= 2) {
        if (get_local_id(0)+offset < get_local_size(0) && (get_local_id(0)&(2*offset-1)) == 0)
            error[get_local_id(0)] += error[get_local_id(0)+offset];
        barrier(CLK_LOCAL_MEM_FENCE);
    }
    if (get_local_id(0) == 0) {
        float totalError = sqrt(error[0]/(numAtoms*3));
        float newStepSize = sqrt(errorTol/totalError);
        float oldStepSize = dt[0].y;
        if (oldStepSize > 0.0f)
            newStepSize = min(newStepSize, oldStepSize*2.0f); // For safety, limit how quickly dt can increase.
        if (newStepSize > oldStepSize && newStepSize < 1.1f*oldStepSize)
            newStepSize = oldStepSize; // Keeping dt constant between steps improves the behavior of the integrator.
        if (newStepSize > maxStepSize)
            newStepSize = maxStepSize;
        dt[0].y = newStepSize;
    }
}
