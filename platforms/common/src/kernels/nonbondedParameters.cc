/**
 * Compute the nonbonded parameters for particles and exceptions.
 */
KERNEL void computeParameters(GLOBAL mixed* RESTRICT energyBuffer, int includeSelfEnergy, GLOBAL real* RESTRICT globalParams,
        int numAtoms, GLOBAL const float4* RESTRICT baseParticleParams, GLOBAL real4* RESTRICT posq, GLOBAL real* RESTRICT charge,
        GLOBAL float2* RESTRICT sigmaEpsilon, GLOBAL float4* RESTRICT particleParamOffsets, GLOBAL int* RESTRICT particleOffsetIndices,
        GLOBAL real* RESTRICT chargeBuffer
#ifdef HAS_EXCEPTIONS
        , int numExceptions, GLOBAL const float4* RESTRICT baseExceptionParams, GLOBAL float4* RESTRICT exceptionParams,
        GLOBAL float4* RESTRICT exceptionParamOffsets, GLOBAL int* RESTRICT exceptionOffsetIndices
#endif
        ) {
    mixed energy = 0;
    real totalCharge = 0;

    // Compute particle parameters.
    
    for (int i = GLOBAL_ID; i < numAtoms; i += GLOBAL_SIZE) {
        float4 params = baseParticleParams[i];
#ifdef HAS_PARTICLE_OFFSETS
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
        sigmaEpsilon[i] = make_float2(0.5f*params.y, 2*SQRT(params.z));
        totalCharge += params.x;
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
    for (int i = GLOBAL_ID; i < numExceptions; i += GLOBAL_SIZE) {
        float4 params = baseExceptionParams[i];
#ifdef HAS_EXCEPTION_OFFSETS
        int start = exceptionOffsetIndices[i], end = exceptionOffsetIndices[i+1];
        for (int j = start; j < end; j++) {
            float4 offset = exceptionParamOffsets[j];
            real value = globalParams[(int) offset.w];
            params.x += value*offset.x;
            params.y += value*offset.y;
            params.z += value*offset.z;
        }
#endif
        exceptionParams[i] = make_float4((float) (ONE_4PI_EPS0*params.x), (float) params.y, (float) (4*params.z), 0);
    }
#endif
    if (includeSelfEnergy)
        energyBuffer[GLOBAL_ID] += energy;

    // Record the total charge from particles processed by this block.

#if defined(HAS_OFFSETS) && defined(INCLUDE_EWALD)
    LOCAL real temp[WORK_GROUP_SIZE];
    temp[LOCAL_ID] = totalCharge;
    for (int i = 1; i < WORK_GROUP_SIZE; i *= 2) {
        SYNC_THREADS;
        if (LOCAL_ID%(i*2) == 0 && LOCAL_ID+i < WORK_GROUP_SIZE)
            temp[LOCAL_ID] += temp[LOCAL_ID+i];
    }
    if (LOCAL_ID == 0)
        chargeBuffer[GROUP_ID] = temp[0];
#endif
}

/**
 * Compute parameters for subtracting the reciprocal part of excluded interactions.
 */
KERNEL void computeExclusionParameters(GLOBAL real4* RESTRICT posq, GLOBAL real* RESTRICT charge, GLOBAL float2* RESTRICT sigmaEpsilon,
        int numExclusions, GLOBAL const int2* RESTRICT exclusionAtoms, GLOBAL float4* RESTRICT exclusionParams) {
    for (int i = GLOBAL_ID; i < numExclusions; i += GLOBAL_SIZE) {
        int2 atoms = exclusionAtoms[i];
#ifdef USE_POSQ_CHARGES
        real chargeProd = posq[atoms.x].w*posq[atoms.y].w;
#else
        real chargeProd = charge[atoms.x]*charge[atoms.y];
#endif
#ifdef INCLUDE_LJPME_EXCEPTIONS
        float2 sigEps1 = sigmaEpsilon[atoms.x];
        float2 sigEps2 = sigmaEpsilon[atoms.y];
        float sigma = sigEps1.x*sigEps2.x;
        float epsilon = sigEps1.y*sigEps2.y;
#else
        float sigma = 0;
        float epsilon = 0;
#endif
        exclusionParams[i] = make_float4((float) (ONE_4PI_EPS0*chargeProd), sigma, epsilon, 0);
    }
}

/**
 * When using Ewald or PME with parameter offsets, the total charge can change each step.
 * We therefore need to compute the correction for the neutralizing plasma on the GPU.
 * This kernel is executed by a single thread block.
 */
KERNEL void computePlasmaCorrection(GLOBAL real* RESTRICT chargeBuffer, GLOBAL mixed* RESTRICT energyBuffer,
        real alpha, real volume) {
    LOCAL real temp[WORK_GROUP_SIZE];
    real sum = 0;
    for (unsigned int index = LOCAL_ID; index < NUM_GROUPS; index += LOCAL_SIZE)
        sum += chargeBuffer[index];
    temp[LOCAL_ID] = sum;
    for (int i = 1; i < WORK_GROUP_SIZE; i *= 2) {
        SYNC_THREADS;
        if (LOCAL_ID%(i*2) == 0 && LOCAL_ID+i < WORK_GROUP_SIZE)
            temp[LOCAL_ID] += temp[LOCAL_ID+i];
    }
    if (LOCAL_ID == 0) {
        real totalCharge = temp[0];
        energyBuffer[0] -= totalCharge*totalCharge/(8*EPSILON0*volume*alpha*alpha);
    }
}