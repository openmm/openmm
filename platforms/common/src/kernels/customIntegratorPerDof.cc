#ifdef SUPPORTS_DOUBLE_PRECISION
typedef double TempType;
typedef double3 TempType3;
typedef double4 TempType4;

#define make_TempType3(a...) make_double3(a)
#define make_TempType4(a...) make_double4(a)
#define convertToTempType3(a) make_double3((a).x, (a).y, (a).z)
#define convertToTempType4(a) make_double4((a).x, (a).y, (a).z, (a).w)

inline DEVICE mixed4 convertFromDouble4(double4 a) {
    return make_mixed4(a.x, a.y, a.z, a.w);
}
#else
typedef float TempType;
typedef float3 TempType3;
typedef float4 TempType4;

#define make_TempType3(a...) make_float3(a)
#define make_TempType4(a...) make_float4(a)
#define convertToTempType3(a) make_float3((a).x, (a).y, (a).z)
#define convertToTempType4(a) make_float4((a).x, (a).y, (a).z, (a).w)
#endif

/**
 * Load the position of a particle.
 */
inline DEVICE TempType4 loadPos(GLOBAL const real4* RESTRICT posq, GLOBAL const real4* RESTRICT posqCorrection, int index) {
#ifdef USE_MIXED_PRECISION
    real4 pos1 = posq[index];
    real4 pos2 = posqCorrection[index];
    return make_TempType4(pos1.x+(mixed)pos2.x, pos1.y+(mixed)pos2.y, pos1.z+(mixed)pos2.z, pos1.w);
#else
    return convertToTempType4(posq[index]);
#endif
}

/**
 * Store the position of a particle.
 */
inline DEVICE void storePos(GLOBAL real4* RESTRICT posq, GLOBAL real4* RESTRICT posqCorrection, int index, TempType4 pos) {
#ifdef USE_MIXED_PRECISION
    posq[index] = make_real4((real) pos.x, (real) pos.y, (real) pos.z, (real) pos.w);
    posqCorrection[index] = make_real4(pos.x-(real) pos.x, pos.y-(real) pos.y, pos.z-(real) pos.z, 0);
#else
    posq[index] = make_real4(pos.x, pos.y, pos.z, pos.w);
#endif
}

KERNEL void computePerDof(GLOBAL real4* RESTRICT posq, GLOBAL real4* RESTRICT posqCorrection, GLOBAL mixed4* RESTRICT posDelta,
        GLOBAL mixed4* RESTRICT velm, GLOBAL const mm_long* RESTRICT force, GLOBAL const mixed2* RESTRICT dt, GLOBAL const mixed* RESTRICT globals,
        GLOBAL mixed* RESTRICT sum, GLOBAL const float4* RESTRICT gaussianValues, unsigned int gaussianBaseIndex, GLOBAL const float4* RESTRICT uniformValues,
        const mixed energy, GLOBAL mixed* RESTRICT energyParamDerivs
        PARAMETER_ARGUMENTS) {
    TempType3 stepSize = make_TempType3(dt[0].y);
    int index = GLOBAL_ID;
    const TempType forceScale = ((TempType) 1)/0xFFFFFFFF;
    while (index < NUM_ATOMS) {
#ifdef LOAD_POS_AS_DELTA
        TempType4 position = loadPos(posq, posqCorrection, index) + convertToTempType4(posDelta[index]);
#else
        TempType4 position = loadPos(posq, posqCorrection, index);
#endif
        TempType4 velocity = convertToTempType4(velm[index]);
        TempType3 f = make_TempType3(forceScale*force[index], forceScale*force[index+PADDED_NUM_ATOMS], forceScale*force[index+PADDED_NUM_ATOMS*2]);
        TempType3 mass = make_TempType3(RECIP(velocity.w));
        if (velocity.w != 0.0) {
            int gaussianIndex = gaussianBaseIndex;
            int uniformIndex = 0;
            COMPUTE_STEP
        }
        index += GLOBAL_SIZE;
    }
}
