/**
 * Load the position of a particle.
 */
inline __device__ mixed4 loadPos(const real4* __restrict__ posq, const real4* __restrict__ posqCorrection, int index) {
#ifdef USE_MIXED_PRECISION
    real4 pos1 = posq[index];
    real4 pos2 = posqCorrection[index];
    return make_mixed4(pos1.x+(mixed)pos2.x, pos1.y+(mixed)pos2.y, pos1.z+(mixed)pos2.z, pos1.w);
#else
    return posq[index];
#endif
}

/**
 * Store the position of a particle.
 */
inline __device__ void storePos(real4* __restrict__ posq, real4* __restrict__ posqCorrection, int index, mixed4 pos) {
#ifdef USE_MIXED_PRECISION
    posq[index] = make_real4((real) pos.x, (real) pos.y, (real) pos.z, (real) pos.w);
    posqCorrection[index] = make_real4(pos.x-(real) pos.x, pos.y-(real) pos.y, pos.z-(real) pos.z, 0);
#else
    posq[index] = pos;
#endif
}

inline __device__ double4 convertToDouble4(float4 a) {
    return make_double4(a.x, a.y, a.z, a.w);
}

inline __device__ double4 convertToDouble4(double4 a) {
    return a;
}

inline __device__ mixed4 convertFromDouble4(double4 a) {
    return make_mixed4(a.x, a.y, a.z, a.w);
}

extern "C" __global__ void computePerDof(real4* __restrict__ posq, real4* __restrict__ posqCorrection, mixed4* __restrict__ posDelta,
        mixed4* __restrict__ velm, const long long* __restrict__ force, const mixed2* __restrict__ dt, const mixed* __restrict__ globals,
        mixed* __restrict__ sum, const float4* __restrict__ gaussianValues, unsigned int gaussianBaseIndex, const float4* __restrict__ uniformValues,
        const mixed energy, mixed* __restrict__ energyParamDerivs
        PARAMETER_ARGUMENTS) {
    double3 stepSize = make_double3(dt[0].y);
    int index = blockIdx.x*blockDim.x+threadIdx.x;
    const double forceScale = 1.0/0xFFFFFFFF;
    while (index < NUM_ATOMS) {
#ifdef LOAD_POS_AS_DELTA
        double4 position = convertToDouble4(loadPos(posq, posqCorrection, index)+posDelta[index]);
#else
        double4 position = convertToDouble4(loadPos(posq, posqCorrection, index));
#endif
        double4 velocity = convertToDouble4(velm[index]);
        double4 f = make_double4(forceScale*force[index], forceScale*force[index+PADDED_NUM_ATOMS], forceScale*force[index+PADDED_NUM_ATOMS*2], 0.0);
        double3 mass = make_double3(1.0/velocity.w);
        if (velocity.w != 0.0) {
            int gaussianIndex = gaussianBaseIndex;
            int uniformIndex = 0;
            COMPUTE_STEP
        }
        index += blockDim.x*gridDim.x;
    }
}
