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

/**
 * Compute the positions of virtual sites
 */
__kernel void computeVirtualSites(__global real4* restrict posq,
#ifdef USE_MIXED_PRECISION
        __global real4* restrict posqCorrection,
#endif
        __global const int4* restrict avg2Atoms, __global const real2* restrict avg2Weights,
        __global const int4* restrict avg3Atoms, __global const real4* restrict avg3Weights,
        __global const int4* restrict outOfPlaneAtoms, __global const real4* restrict outOfPlaneWeights) {
#ifndef USE_MIXED_PRECISION
        __global real4* posqCorrection = 0;
#endif
    
    // Two particle average sites.
    
    for (int index = get_global_id(0); index < NUM_2_AVERAGE; index += get_global_size(0)) {
        int4 atoms = avg2Atoms[index];
        real2 weights = avg2Weights[index];
        mixed4 pos = loadPos(posq, posqCorrection, atoms.x);
        mixed4 pos1 = loadPos(posq, posqCorrection, atoms.y);
        mixed4 pos2 = loadPos(posq, posqCorrection, atoms.z);
        pos.xyz = pos1.xyz*weights.x + pos2.xyz*weights.y;
        storePos(posq, posqCorrection, atoms.x, pos);
    }
    
    // Three particle average sites.
    
    for (int index = get_global_id(0); index < NUM_3_AVERAGE; index += get_global_size(0)) {
        int4 atoms = avg3Atoms[index];
        real4 weights = avg3Weights[index];
        mixed4 pos = loadPos(posq, posqCorrection, atoms.x);
        mixed4 pos1 = loadPos(posq, posqCorrection, atoms.y);
        mixed4 pos2 = loadPos(posq, posqCorrection, atoms.z);
        mixed4 pos3 = loadPos(posq, posqCorrection, atoms.w);
        pos.xyz = pos1.xyz*weights.x + pos2.xyz*weights.y + pos3.xyz*weights.z;
        storePos(posq, posqCorrection, atoms.x, pos);
    }
    
    // Out of plane sites.
    
    for (int index = get_global_id(0); index < NUM_OUT_OF_PLANE; index += get_global_size(0)) {
        int4 atoms = outOfPlaneAtoms[index];
        real4 weights = outOfPlaneWeights[index];
        mixed4 pos = loadPos(posq, posqCorrection, atoms.x);
        mixed4 pos1 = loadPos(posq, posqCorrection, atoms.y);
        mixed4 pos2 = loadPos(posq, posqCorrection, atoms.z);
        mixed4 pos3 = loadPos(posq, posqCorrection, atoms.w);
        mixed4 v12 = pos2-pos1;
        mixed4 v13 = pos3-pos1;
        pos.xyz = pos1.xyz + v12.xyz*weights.x + v13.xyz*weights.y + cross(v12, v13).xyz*weights.z;
        storePos(posq, posqCorrection, atoms.x, pos);
    }
}

/**
 * Distribute forces from virtual sites to the atoms they are based on.
 */
__kernel void distributeForces(__global const real4* restrict posq, __global real4* restrict force,
#ifdef USE_MIXED_PRECISION
        __global real4* restrict posqCorrection,
#endif
        __global const int4* restrict avg2Atoms, __global const real2* restrict avg2Weights,
        __global const int4* restrict avg3Atoms, __global const real4* restrict avg3Weights,
        __global const int4* restrict outOfPlaneAtoms, __global const real4* restrict outOfPlaneWeights) {
#ifndef USE_MIXED_PRECISION
        __global real4* posqCorrection = 0;
#endif

    // Two particle average sites.
    
    for (int index = get_global_id(0); index < NUM_2_AVERAGE; index += get_global_size(0)) {
        int4 atoms = avg2Atoms[index];
        real2 weights = avg2Weights[index];
        real4 f = force[atoms.x];
        real4 f1 = force[atoms.y];
        real4 f2 = force[atoms.z];
        f1.xyz += f.xyz*weights.x;
        f2.xyz += f.xyz*weights.y;
        force[atoms.y] = f1;
        force[atoms.z] = f2;
    }
    
    // Three particle average sites.
    
    for (int index = get_global_id(0); index < NUM_3_AVERAGE; index += get_global_size(0)) {
        int4 atoms = avg3Atoms[index];
        real4 weights = avg3Weights[index];
        real4 f = force[atoms.x];
        real4 f1 = force[atoms.y];
        real4 f2 = force[atoms.z];
        real4 f3 = force[atoms.w];
        f1.xyz += f.xyz*weights.x;
        f2.xyz += f.xyz*weights.y;
        f3.xyz += f.xyz*weights.z;
        force[atoms.y] = f1;
        force[atoms.z] = f2;
        force[atoms.w] = f3;
    }
    
    // Out of plane sites.
    
    for (int index = get_global_id(0); index < NUM_OUT_OF_PLANE; index += get_global_size(0)) {
        int4 atoms = outOfPlaneAtoms[index];
        real4 weights = outOfPlaneWeights[index];
        mixed4 pos1 = loadPos(posq, posqCorrection, atoms.y);
        mixed4 pos2 = loadPos(posq, posqCorrection, atoms.z);
        mixed4 pos3 = loadPos(posq, posqCorrection, atoms.w);
        mixed4 v12 = pos2-pos1;
        mixed4 v13 = pos3-pos1;
        real4 f = force[atoms.x];
        real4 f1 = force[atoms.y];
        real4 f2 = force[atoms.z];
        real4 f3 = force[atoms.w];
        real4 fp2 = (real4) (weights.x*f.x - weights.z*v13.z*f.y + weights.z*v13.y*f.z,
                   weights.z*v13.z*f.x + weights.x*f.y - weights.z*v13.x*f.z,
                  -weights.z*v13.y*f.x + weights.z*v13.x*f.y + weights.x*f.z, 0.0f);
        real4 fp3 = (real4) (weights.y*f.x + weights.z*v12.z*f.y - weights.z*v12.y*f.z,
                  -weights.z*v12.z*f.x + weights.y*f.y + weights.z*v12.x*f.z,
                   weights.z*v12.y*f.x - weights.z*v12.x*f.y + weights.y*f.z, 0.0f);
        f1.xyz += f.xyz-fp2.xyz-fp3.xyz;
        f2.xyz += fp2.xyz;
        f3.xyz += fp3.xyz;
        force[atoms.y] = f1;
        force[atoms.z] = f2;
        force[atoms.w] = f3;
    }
}
