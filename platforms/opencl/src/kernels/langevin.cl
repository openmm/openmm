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

/**
 * Select the step size to use for the next step.
 */

__kernel void selectLangevinStepSize(int numAtoms, float maxStepSize, float errorTol, float tau, float kT, __global float2* dt,
        __global float4* velm, __global float4* force, __global float* paramBuffer, __local float* params, __local float* error) {
    // Calculate the error.

    float err = 0.0f;
    unsigned int index = get_local_id(0);
    while (index < numAtoms) {
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

        // Recalculate the integration parameters.

        float gdt                  = newStepSize/tau;
        float eph                  = exp(0.5f*gdt);
        float emh                  = exp(-0.5f*gdt);
        float ep                   = exp(gdt);
        float em                   = exp(-gdt);
        float em_v                 = exp(-0.5f*(oldStepSize+newStepSize)/tau);
        float b, c, d;
        if (gdt >= 0.1f) {
            float term1 = eph-1.0f;
            term1                 *= term1;
            b                      = gdt*(ep-1.0f) - 4.0f*term1;
            c                      = gdt - 3.0f + 4.0f*emh - em;
            d                      = 2.0f - eph - emh;
        }
        else {
            float term1            = 0.5f*gdt;
            float term2            = term1*term1;
            float term4            = term2*term2;

            float third            = 1.0f/3.0f;
            float o7_9             = 7.0f/9.0f;
            float o1_12            = 1.0f/12.0f;
            float o17_90           = 17.0f/90.0f;
            float o7_30            = 7.0f/30.0f;
            float o31_1260         = 31.0f/1260.0f;
            float o_360            = 1.0f/360.0f;

            b                      = term4*(third + term1*(third + term1*(o17_90 + term1*o7_9)));
            c                      = term2*term1*(2.0f*third + term1*(-0.5f + term1*(o7_30 + term1*(-o1_12 + term1*o31_1260))));
            d                      = term2*(-1.0f + term2*(-o1_12 - term2*o_360));
        }
        float fix1                 = tau*(eph - emh);
        if (fix1 == 0.0f)
            fix1 = newStepSize;
        params[EM]                 = em;
        params[EM_V]               = em_v;
        params[DOverTauC]          = d/(tau*c);
        params[TauOneMinusEM_V]    = tau*(1.0f-em_v);
        params[TauDOverEMMinusOne] = tau*d/(em - 1.0f);
        params[Fix1]               = fix1;
        params[OneOverFix1]        = 1.0f/fix1;
        params[V]                  = sqrt(kT*(1.0f - em));
        params[X]                  = tau*sqrt(kT*c);
        params[Yv]                 = sqrt(kT*b/c);
        params[Yx]                 = tau*sqrt(kT*b/(1.0f - em));
    }
    barrier(CLK_LOCAL_MEM_FENCE);
    if (get_local_id(0) < MaxParams)
        paramBuffer[get_local_id(0)] = params[get_local_id(0)];
}
