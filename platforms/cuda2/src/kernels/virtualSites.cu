/**
 * Compute the positions of virtual sites
 */
extern "C" __global__ void computeVirtualSites(real4* __restrict__ posq, const int4* __restrict__ avg2Atoms, const real2* __restrict__ avg2Weights,
        const int4* __restrict__ avg3Atoms, const real4* __restrict__ avg3Weights,
        const int4* __restrict__ outOfPlaneAtoms, const real4* __restrict__ outOfPlaneWeights) {
    
    // Two particle average sites.
    
    for (int index = blockIdx.x*blockDim.x+threadIdx.x; index < NUM_2_AVERAGE; index += blockDim.x*gridDim.x) {
        int4 atoms = avg2Atoms[index];
        real2 weights = avg2Weights[index];
        real4 pos = posq[atoms.x];
        real4 pos1 = posq[atoms.y];
        real4 pos2 = posq[atoms.z];
        pos.x = pos1.x*weights.x + pos2.x*weights.y;
        pos.y = pos1.y*weights.x + pos2.y*weights.y;
        pos.z = pos1.z*weights.x + pos2.z*weights.y;
        posq[atoms.x] = pos;
    }
    
    // Three particle average sites.
    
    for (int index = blockIdx.x*blockDim.x+threadIdx.x; index < NUM_3_AVERAGE; index += blockDim.x*gridDim.x) {
        int4 atoms = avg3Atoms[index];
        real4 weights = avg3Weights[index];
        real4 pos = posq[atoms.x];
        real4 pos1 = posq[atoms.y];
        real4 pos2 = posq[atoms.z];
        real4 pos3 = posq[atoms.w];
        pos.x = pos1.x*weights.x + pos2.x*weights.y + pos3.x*weights.z;
        pos.y = pos1.y*weights.x + pos2.y*weights.y + pos3.y*weights.z;
        pos.z = pos1.z*weights.x + pos2.z*weights.y + pos3.z*weights.z;
        posq[atoms.x] = pos;
    }
    
    // Out of plane sites.
    
    for (int index = blockIdx.x*blockDim.x+threadIdx.x; index < NUM_OUT_OF_PLANE; index += blockDim.x*gridDim.x) {
        int4 atoms = outOfPlaneAtoms[index];
        real4 weights = outOfPlaneWeights[index];
        real4 pos = posq[atoms.x];
        real4 pos1 = posq[atoms.y];
        real4 pos2 = posq[atoms.z];
        real4 pos3 = posq[atoms.w];
        real4 v12 = pos2-pos1;
        real4 v13 = pos3-pos1;
        real3 cr = cross(v12, v13);
        pos.x = pos1.x + v12.x*weights.x + v13.x*weights.y + cr.x*weights.z;
        pos.y = pos1.y + v12.y*weights.x + v13.y*weights.y + cr.y*weights.z;
        pos.z = pos1.z + v12.z*weights.x + v13.z*weights.y + cr.z*weights.z;
        posq[atoms.x] = pos;
    }
}

inline __device__ real3 loadForce(int index, long long* __restrict__ force) {
    real scale = RECIP((real) 0xFFFFFFFF);
    return make_real3(scale*force[index], scale*force[index+PADDED_NUM_ATOMS], scale*force[index+PADDED_NUM_ATOMS*2]);
}

inline __device__ void storeForce(int index, long long* __restrict__ force, real3 value) {
    force[index] = (long long) (value.x*0xFFFFFFFF);
    force[index+PADDED_NUM_ATOMS] = (long long) (value.y*0xFFFFFFFF);
    force[index+2*PADDED_NUM_ATOMS] = (long long) (value.z*0xFFFFFFFF);
}

/**
 * Distribute forces from virtual sites to the atoms they are based on.
 */
extern "C" __global__ void distributeForces(const real4* __restrict__ posq, long long* __restrict__ force,
        const int4* __restrict__ avg2Atoms, const real2* __restrict__ avg2Weights,
        const int4* __restrict__ avg3Atoms, const real4* __restrict__ avg3Weights,
        const int4* __restrict__ outOfPlaneAtoms, const real4* __restrict__ outOfPlaneWeights) {
    
    // Two particle average sites.
    
    for (int index = blockIdx.x*blockDim.x+threadIdx.x; index < NUM_2_AVERAGE; index += blockDim.x*gridDim.x) {
        int4 atoms = avg2Atoms[index];
        real2 weights = avg2Weights[index];
        real3 f = loadForce(atoms.x, force);
        real3 f1 = loadForce(atoms.y, force);
        real3 f2 = loadForce(atoms.z, force);
        f1 += f*weights.x;
        f2 += f*weights.y;
        storeForce(atoms.y, force, f1);
        storeForce(atoms.z, force, f2);
    }
    
    // Three particle average sites.
    
    for (int index = blockIdx.x*blockDim.x+threadIdx.x; index < NUM_3_AVERAGE; index += blockDim.x*gridDim.x) {
        int4 atoms = avg3Atoms[index];
        real4 weights = avg3Weights[index];
        real3 f = loadForce(atoms.x, force);
        real3 f1 = loadForce(atoms.y, force);
        real3 f2 = loadForce(atoms.z, force);
        real3 f3 = loadForce(atoms.w, force);
        f1 += f*weights.x;
        f2 += f*weights.y;
        f3 += f*weights.z;
        storeForce(atoms.y, force, f1);
        storeForce(atoms.z, force, f2);
        storeForce(atoms.w, force, f3);
    }
    
    // Out of plane sites.
    
    for (int index = blockIdx.x*blockDim.x+threadIdx.x; index < NUM_OUT_OF_PLANE; index += blockDim.x*gridDim.x) {
        int4 atoms = outOfPlaneAtoms[index];
        real4 weights = outOfPlaneWeights[index];
        real4 pos1 = posq[atoms.y];
        real4 pos2 = posq[atoms.z];
        real4 pos3 = posq[atoms.w];
        real4 v12 = pos2-pos1;
        real4 v13 = pos3-pos1;
        real3 f = loadForce(atoms.x, force);
        real3 f1 = loadForce(atoms.y, force);
        real3 f2 = loadForce(atoms.z, force);
        real3 f3 = loadForce(atoms.w, force);
        real3 fp2 = make_real3(weights.x*f.x - weights.z*v13.z*f.y + weights.z*v13.y*f.z,
                   weights.z*v13.z*f.x + weights.x*f.y - weights.z*v13.x*f.z,
                  -weights.z*v13.y*f.x + weights.z*v13.x*f.y + weights.x*f.z);
        real3 fp3 = make_real3(weights.y*f.x + weights.z*v12.z*f.y - weights.z*v12.y*f.z,
                  -weights.z*v12.z*f.x + weights.y*f.y + weights.z*v12.x*f.z,
                   weights.z*v12.y*f.x - weights.z*v12.x*f.y + weights.y*f.z);
        f1 += f-fp2-fp3;
        f2 += fp2;
        f3 += fp3;
        storeForce(atoms.y, force, f1);
        storeForce(atoms.z, force, f2);
        storeForce(atoms.w, force, f3);
    }
}
