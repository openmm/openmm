/**
 * Compute the nonbonded parameters for particles and exceptions.
 */
__kernel void computeParameters(__global mixed* restrict energyBuffer, int includeSelfEnergy, __global real* restrict globalParams,
        int numAtoms, __global const float4* restrict baseParticleParams, __global real4* restrict posq, __global real* restrict charge,
        __global float2* restrict sigmaEpsilon, __global float4* restrict particleParamOffsets, __global int* restrict particleOffsetIndices
#ifdef HAS_EXCEPTIONS
        , int numExceptions, __global const float4* restrict baseExceptionParams, __global float4* restrict exceptionParams,
        __global float4* restrict exceptionParamOffsets, __global int* restrict exceptionOffsetIndices
#endif
        ) {
    mixed energy = 0;

    // Compute particle parameters.
    
    for (int i = get_global_id(0); i < numAtoms; i += get_global_size(0)) {
        float4 params = baseParticleParams[i];
#ifdef HAS_OFFSETS
        int start = particleOffsetIndices[i], end = particleOffsetIndices[i+1];
        for (int j = start; j < end; j++) {
            float4 offset = particleParamOffsets[j];
            real value = globalParams[(int) offset.w];
            params.x += value*offset.x;
            params.y += value*offset.y;
            params.z += value*offset.z;
        }
#endif
#ifdef USE_POSQ_CHARGES
        posq[i].w = params[0];
#else
        charge[i] = params[0];
#endif
        sigmaEpsilon[i] = (float2) (0.5f*params[1], 2*SQRT(params[2]));
#ifdef INCLUDE_EWALD
        energy -= EWALD_SELF_ENERGY_SCALE*params[0]*params[0];
#endif
#ifdef INCLUDE_LJPME
        real sig3 = params[1]*params[1]*params[1];
        energy += LJPME_SELF_ENERGY_SCALE*sig3*sig3*params[2];
#endif
    }

    // Compute exception parameters.
    
#ifdef HAS_EXCEPTIONS
    for (int i = get_global_id(0); i < numExceptions; i += get_global_size(0)) {
        float4 params = baseExceptionParams[i];
#ifdef HAS_OFFSETS
        int start = exceptionOffsetIndices[i], end = exceptionOffsetIndices[i+1];
        for (int j = start; j < end; j++) {
            float4 offset = exceptionParamOffsets[j];
            real value = globalParams[(int) offset.w];
            params.x += value*offset.x;
            params.y += value*offset.y;
            params.z += value*offset.z;
        }
#endif
        exceptionParams[i] = (float4) ((float) (138.935456f*params[0]), (float) params[1], (float) (4*params[2]), 0);
    }
#endif
    if (includeSelfEnergy)
        energyBuffer[get_global_id(0)] += energy;
}
