enum {VelScale, ForceScale, NoiseScale, MaxParams};

/**
 * Perform the first step of Langevin integration.
 */

__kernel void integrateLangevinPart1(__global float4* velm, __global float4* force, __global float4* posDelta,
        __global float* paramBuffer, __global float2* dt, __global float4* random, unsigned int randomIndex) {
    float vscale = paramBuffer[VelScale];
    float fscale = paramBuffer[ForceScale];
    float noisescale = paramBuffer[NoiseScale];
    float stepSize = dt[0].y;
    int index = get_global_id(0);
    randomIndex += index;
    while (index < NUM_ATOMS) {
        float4 velocity = velm[index];
        float sqrtInvMass = sqrt(velocity.w);
        velocity.xyz = vscale*velocity.xyz + fscale*velocity.w*force[index].xyz + noisescale*sqrtInvMass*random[randomIndex].xyz;
        velm[index] = velocity;
        posDelta[index] = stepSize*velocity;
        randomIndex += get_global_size(0);
        index += get_global_size(0);
    }
}

/**
 * Perform the second step of Langevin integration.
 */

__kernel void integrateLangevinPart2(__global float4* posq, __global float4* posDelta, __global float4* velm, __global float2* dt) {
    float invStepSize = 1.0f/dt[0].y;
    int index = get_global_id(0);
    while (index < NUM_ATOMS) {
        float4 pos = posq[index];
        float4 delta = posDelta[index];
        float4 vel = velm[index];
        pos.xyz += delta.xyz;
        vel.xyz = invStepSize*delta.xyz;
        posq[index] = pos;
        velm[index] = vel;
        index += get_global_size(0);
    }
}

/**
 * Select the step size to use for the next step.
 */

__kernel void selectLangevinStepSize(float maxStepSize, float errorTol, float tau, float kT, __global float2* dt,
        __global float4* velm, __global float4* force, __global float* paramBuffer, __local float* params, __local float* error) {
    // Calculate the error.

    float err = 0.0f;
    unsigned int index = get_local_id(0);
    while (index < NUM_ATOMS) {
        float4 f = force[index];
        float invMass = velm[index].w;
        err += (f.x*f.x + f.y*f.y + f.z*f.z)*invMass;
        index += get_global_size(0);
    }
    error[get_local_id(0)] = err;
    barrier(CLK_LOCAL_MEM_FENCE);

    // Sum the errors from all threads.

    for (int offset = 1; offset < get_local_size(0); offset *= 2) {
        if (get_local_id(0)+offset < get_local_size(0) && (get_local_id(0)&(2*offset-1)) == 0)
            error[get_local_id(0)] += error[get_local_id(0)+offset];
        barrier(CLK_LOCAL_MEM_FENCE);
    }
    if (get_global_id(0) == 0) {
        // Select the new step size.

        float totalError = sqrt(error[0]/(NUM_ATOMS*3));
        float newStepSize = sqrt(errorTol/totalError);
        float oldStepSize = dt[0].y;
        if (oldStepSize > 0.0f)
            newStepSize = min(newStepSize, oldStepSize*2.0f); // For safety, limit how quickly dt can increase.
        if (newStepSize > oldStepSize && newStepSize < 1.1f*oldStepSize)
            newStepSize = oldStepSize; // Keeping dt constant between steps improves the behavior of the integrator.
        if (newStepSize > maxStepSize)
            newStepSize = maxStepSize;
        dt[0].y = newStepSize;

        // Recalculate the integration parameters.

        float vscale = exp(-newStepSize/tau);
        float fscale = (1-vscale)*tau;
        float noisescale = sqrt(2*kT/tau)*sqrt(0.5f*(1-vscale*vscale)*tau);
        params[VelScale] = vscale;
        params[ForceScale] = fscale;
        params[NoiseScale] = noisescale;
    }
    barrier(CLK_LOCAL_MEM_FENCE);
    if (get_local_id(0) < MaxParams)
        paramBuffer[get_local_id(0)] = params[get_local_id(0)];
}
