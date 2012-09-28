inline __device__ double4 convertToDouble4(real4 a) {
    return make_double4(a.x, a.y, a.z, a.w);
}

inline __device__ real4 convertFromDouble4(double4 a) {
    return make_real4(a.x, a.y, a.z, a.w);
}

extern "C" __global__ void computePerDof(real4* __restrict__ posq, real4* __restrict__ posDelta, real4* __restrict__ velm,
        const long long* __restrict__ force, const real2* __restrict__ dt, const real* __restrict__ globals,
        const real* __restrict__ params, real* __restrict__ sum, const float4* __restrict__ gaussianValues,
        unsigned int randomIndex, const float4* __restrict__ uniformValues, const real* __restrict__ energy
        PARAMETER_ARGUMENTS) {
    real stepSize = dt[0].y;
    int index = blockIdx.x*blockDim.x+threadIdx.x;
    randomIndex += index;
    const double forceScale = 1.0/0xFFFFFFFF;
    while (index < NUM_ATOMS) {
#ifdef LOAD_POS_AS_DELTA
        double4 position = convertToDouble4(posq[index]+posDelta[index]);
#else
        double4 position = convertToDouble4(posq[index]);
#endif
        double4 velocity = convertToDouble4(velm[index]);
        double4 f = make_double4(forceScale*force[index], forceScale*force[index+PADDED_NUM_ATOMS], forceScale*force[index+PADDED_NUM_ATOMS*2], 0.0);
        double mass = 1.0/velocity.w;
        if (velocity.w != 0.0) {
            float4 gaussian = gaussianValues[randomIndex];
            float4 uniform = uniformValues[index];
            COMPUTE_STEP
        }
        randomIndex += blockDim.x*gridDim.x;
        index += blockDim.x*gridDim.x;
    }
}
