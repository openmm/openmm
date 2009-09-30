enum {EM, EM_V, DOverTauC, TauOneMinusEM_V, TauDOverEMMinusOne, V, X, Yv, Yx, Fix1, OneOverFix1, MaxParams};

/**
 * Perform the first step of Langevin integration.
 */

__kernel void integrateLangevinPart1(int numAtoms, __global float4* velm, __global float4* force, __global float4* posDelta,
        __global float* paramBuffer, __local float* params, __global float4* xVector, __global float4* vVector,
        __global float4* random, unsigned int randomIndex) {

    // Load the parameters into local memory for faster access.

    int index = get_global_id(0);
    if (index < MaxParams)
        params[index] = paramBuffer[index];
    barrier(CLK_LOCAL_MEM_FENCE);
    randomIndex += index;
    while (index < numAtoms) {
        float4 velocity = velm[index];
        float sqrtInvMass = sqrt(velocity.w);
        float4 vmh = (float4) (xVector[index].xyz*params[DOverTauC] + sqrtInvMass*params[Yv]*random[randomIndex].xyz, 0.0f);
        randomIndex += get_global_size(0);
        float4 vVec = (float4) (sqrtInvMass*params[V]*random[randomIndex].xyz, 0.0f);
        randomIndex += get_global_size(0);
        vVector[index] = vVec;
        velocity.xyz = velocity.xyz*params[EM_V] + velocity.w*force[index].xyz*params[TauOneMinusEM_V] + vVec.xyz - params[EM]*vmh.xyz;
        posDelta[index] = velocity*params[Fix1];
        velm[index] = velocity;
        index += get_global_size(0);
    }
}

/**
 * Perform the second step of Langevin integration.
 */

__kernel void integrateLangevinPart2(int numAtoms, __global float4* velm, __global float4* posDelta, __global float* paramBuffer,
        __local float* params, __global float4* xVector, __global float4* vVector, __global float4* random, unsigned int randomIndex) {

    // Load the parameters into local memory for faster access.

    int index = get_global_id(0);
    if (index < MaxParams)
        params[index] = paramBuffer[index];
    barrier(CLK_LOCAL_MEM_FENCE);
    randomIndex += index;
    while (index < numAtoms) {
        float4 delta = posDelta[index];
        float4 velocity = velm[index];
        float sqrtInvMass = sqrt(velocity.w);
        velocity.xyz = delta.xyz*params[OneOverFix1];
        float4 xmh = (float4) (vVector[index].xyz*params[TauDOverEMMinusOne] + sqrtInvMass*params[Yx]*random[randomIndex].xyz, 0.0f);
        randomIndex += get_global_size(0);
        float4 xVec = (float4) (sqrtInvMass*params[X]*random[randomIndex].xyz, 0.0f);
        randomIndex += get_global_size(0);
        delta.xyz += xVec.xyz - xmh.xyz;
        posDelta[index] = delta;
        velm[index] = velocity;
        xVector[index] = xVec;
        index += get_global_size(0);
    }
}

/**
 * Perform the third step of Langevin integration.
 */

__kernel void integrateLangevinPart3(int numAtoms, __global float4* posq, __global float4* posDelta) {
    int index = get_global_id(0);
    while (index < numAtoms) {
        float4 pos = posq[index];
        float4 delta = posDelta[index];
        pos.xyz += delta.xyz;
        posq[index] = pos;
        index += get_global_size(0);
    }
}
