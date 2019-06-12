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
        posq[i].w = params.x;
#else
        charge[i] = params.x;
#endif
        sigmaEpsilon[i] = (float2) (0.5f*params.y, 2*SQRT(params.z));
#ifdef HAS_OFFSETS
    #ifdef INCLUDE_EWALD
        energy -= EWALD_SELF_ENERGY_SCALE*params.x*params.x;
    #endif
    #ifdef INCLUDE_LJPME
        real sig3 = params.y*params.y*params.y;
        energy += LJPME_SELF_ENERGY_SCALE*sig3*sig3*params.z;
    #endif
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
        exceptionParams[i] = (float4) ((float) (138.935456f*params.x), (float) params.y, (float) (4*params.z), 0);
    }
#endif
    if (includeSelfEnergy)
        energyBuffer[get_global_id(0)] += energy;
}

/**
 * Compute parameters for subtracting the reciprocal part of excluded interactions.
 */
__kernel void computeExclusionParameters(__global real4* restrict posq, __global real* restrict charge, __global float2* restrict sigmaEpsilon,
        int numExclusions, __global const int2* restrict exclusionAtoms, __global float4* restrict exclusionParams) {
    for (int i = get_global_id(0); i < numExclusions; i += get_global_size(0)) {
        int2 atoms = exclusionAtoms[i];
#ifdef USE_POSQ_CHARGES
        real chargeProd = posq[atoms.x].w*posq[atoms.y].w;
#else
        real chargeProd = charge[atoms.x]*charge[atoms.y];
#endif
#ifdef INCLUDE_LJPME
        float2 sigEps1 = sigmaEpsilon[atoms.x];
        float2 sigEps2 = sigmaEpsilon[atoms.y];
        float sigma = sigEps1.x*sigEps2.x;
        float epsilon = sigEps1.y*sigEps2.y;
#else
        float sigma = 0;
        float epsilon = 0;
#endif
        exclusionParams[i] = (float4) ((float) (138.935456f*chargeProd), sigma, epsilon, 0);
    }
}
