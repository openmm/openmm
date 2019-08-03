/**
 * Load the position of a particle.
 */
mixed4 loadPos(__global const real4* restrict posq, __global const real4* restrict posqCorrection, int index) {
#ifdef USE_MIXED_PRECISION
    real4 pos1 = posq[index];
    real4 pos2 = posqCorrection[index];
    return (mixed4) (pos1.x+(mixed)pos2.x, pos1.y+(mixed)pos2.y, pos1.z+(mixed)pos2.z, pos1.w);
#else
    return posq[index];
#endif
}

/**
 * Store the position of a particle.
 */
void storePos(__global real4* restrict posq, __global real4* restrict posqCorrection, int index, mixed4 pos) {
#ifdef USE_MIXED_PRECISION
    posq[index] = (real4) ((real) pos.x, (real) pos.y, (real) pos.z, (real) pos.w);
    posqCorrection[index] = (real4) (pos.x-(real) pos.x, pos.y-(real) pos.y, pos.z-(real) pos.z, 0);
#else
    posq[index] = pos;
#endif
}

__kernel void computePerDof(__global real4* restrict posq, __global real4* restrict posqCorrection, __global mixed4* restrict posDelta,
        __global mixed4* restrict velm, __global const real4* restrict force, __global const mixed2* restrict dt, __global const mixed* restrict globals,
        __global mixed* restrict sum, __global const float4* restrict gaussianValues, unsigned int gaussianBaseIndex, __global const float4* restrict uniformValues,
        const mixed energy, __global mixed* restrict energyParamDerivs
        PARAMETER_ARGUMENTS) {
    mixed stepSize = dt[0].y;
    int index = get_global_id(0);
    while (index < NUM_ATOMS) {
#ifdef LOAD_POS_AS_DELTA
        mixed4 position = loadPos(posq, posqCorrection, index)+posDelta[index];
#else
        mixed4 position = loadPos(posq, posqCorrection, index);
#endif
        mixed4 velocity = velm[index];
        mixed4 f = convert_mixed4(force[index]);
        mixed mass = 1/velocity.w;
        if (velocity.w != 0.0) {
            int gaussianIndex = gaussianBaseIndex;
            int uniformIndex = 0;
            COMPUTE_STEP
        }
        index += get_global_size(0);
    }
}
