__kernel void applyPositionDeltas(__global real4* restrict posq, __global real4* restrict posqCorrection, __global mixed4* restrict posDelta) {
    for (unsigned int index = get_global_id(0); index < NUM_ATOMS; index += get_global_size(0)) {
#ifdef USE_MIXED_PRECISION
        real4 pos1 = posq[index];
        real4 pos2 = posqCorrection[index];
        mixed4 pos = (mixed4) (pos1.x+(mixed)pos2.x, pos1.y+(mixed)pos2.y, pos1.z+(mixed)pos2.z, pos1.w);
#else
        mixed4 pos = posq[index];
#endif
        pos.xyz += posDelta[index].xyz;
#ifdef USE_MIXED_PRECISION
        posq[index] = (real4) ((real) pos.x, (real) pos.y, (real) pos.z, (real) pos.w);
        posqCorrection[index] = (real4) (pos.x-(real) pos.x, pos.y-(real) pos.y, pos.z-(real) pos.z, 0);
#else
        posq[index] = pos;
#endif
    }
}
