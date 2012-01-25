/**
 * Compute the positions of virtual sites
 */
__kernel void computeVirtualSites(__global float4* restrict posq, __global const int4* restrict avg2Atoms, __global const float2* restrict avg2Weights,
        __global const int4* restrict avg3Atoms, __global const float4* restrict avg3Weights,
        __global const int4* restrict outOfPlaneAtoms, __global const float4* restrict outOfPlaneWeights) {
    
    // Two particle average sites.
    
    for (int index = get_global_id(0); index < NUM_2_AVERAGE; index += get_global_size(0)) {
        int4 atoms = avg2Atoms[index];
        float2 weights = avg2Weights[index];
        float4 pos = posq[atoms.x];
        float4 pos1 = posq[atoms.y];
        float4 pos2 = posq[atoms.z];
        pos.xyz = pos1.xyz*weights.x + pos2.xyz*weights.y;
        posq[atoms.x] = pos;
    }
    
    // Three particle average sites.
    
    for (int index = get_global_id(0); index < NUM_3_AVERAGE; index += get_global_size(0)) {
        int4 atoms = avg3Atoms[index];
        float4 weights = avg3Weights[index];
        float4 pos = posq[atoms.x];
        float4 pos1 = posq[atoms.y];
        float4 pos2 = posq[atoms.z];
        float4 pos3 = posq[atoms.w];
        pos.xyz = pos1.xyz*weights.x + pos2.xyz*weights.y + pos3.xyz*weights.z;
        posq[atoms.x] = pos;
    }
    
    // Out of plane sites.
    
    for (int index = get_global_id(0); index < NUM_OUT_OF_PLANE; index += get_global_size(0)) {
        int4 atoms = outOfPlaneAtoms[index];
        float4 weights = outOfPlaneWeights[index];
        float4 pos = posq[atoms.x];
        float4 pos1 = posq[atoms.y];
        float4 pos2 = posq[atoms.z];
        float4 pos3 = posq[atoms.w];
        float4 v12 = pos2-pos1;
        float4 v13 = pos3-pos1;
        pos.xyz = pos1.xyz + v12.xyz*weights.x + v13.xyz*weights.y + cross(v12, v13).xyz*weights.z;
        posq[atoms.x] = pos;
    }
}

/**
 * Distribute forces from virtual sites to the atoms they are based on.
 */
__kernel void distributeForces(__global const float4* restrict posq, __global float4* restrict force,
        __global const int4* restrict avg2Atoms, __global const float2* restrict avg2Weights,
        __global const int4* restrict avg3Atoms, __global const float4* restrict avg3Weights,
        __global const int4* restrict outOfPlaneAtoms, __global const float4* restrict outOfPlaneWeights) {
    
    // Two particle average sites.
    
    for (int index = get_global_id(0); index < NUM_2_AVERAGE; index += get_global_size(0)) {
        int4 atoms = avg2Atoms[index];
        float2 weights = avg2Weights[index];
        float4 f = force[atoms.x];
        float4 f1 = force[atoms.y];
        float4 f2 = force[atoms.z];
        f1.xyz += f.xyz*weights.x;
        f2.xyz += f.xyz*weights.y;
        force[atoms.y] = f1;
        force[atoms.z] = f2;
    }
    
    // Three particle average sites.
    
    for (int index = get_global_id(0); index < NUM_3_AVERAGE; index += get_global_size(0)) {
        int4 atoms = avg3Atoms[index];
        float4 weights = avg3Weights[index];
        float4 f = force[atoms.x];
        float4 f1 = force[atoms.y];
        float4 f2 = force[atoms.z];
        float4 f3 = force[atoms.w];
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
        float4 weights = outOfPlaneWeights[index];
        float4 pos1 = posq[atoms.y];
        float4 pos2 = posq[atoms.z];
        float4 pos3 = posq[atoms.w];
        float4 v12 = pos2-pos1;
        float4 v13 = pos3-pos1;
        float4 f = force[atoms.x];
        float4 f1 = force[atoms.y];
        float4 f2 = force[atoms.z];
        float4 f3 = force[atoms.w];
        float4 fp2 = (float4) (weights.x*f.x - weights.z*v13.z*f.y + weights.z*v13.y*f[2],
                   weights.z*v13.z*f.x + weights.x*f.y - weights.z*v13.x*f[2],
                  -weights.z*v13.y*f.x + weights.z*v13.x*f.y + weights.x*f.z, 0.0f);
        float4 fp3 = (float4) (weights.y*f.x + weights.z*v12[2]*f.y - weights.z*v12.y*f.z,
                  -weights.z*v12[2]*f.x + weights.y*f.y + weights.z*v12.x*f.z,
                   weights.z*v12.y*f.x - weights.z*v12.x*f.y + weights.y*f[2], 0.0f);
        f1.xyz += f.xyz-fp2.xyz-fp3.xyz;
        f2.xyz += fp2.xyz;
        f3.xyz += fp3.xyz;
        force[atoms.y] = f1;
        force[atoms.z] = f2;
        force[atoms.w] = f3;
    }
}
