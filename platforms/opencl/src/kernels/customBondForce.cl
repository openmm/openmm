/**
 * Compute custom bond forces.
 */

__kernel void computeCustomBondForces(int numAtoms, int numBonds, __global float4* forceBuffers, __global float* energyBuffer,
        __global float4* posq, __global int4* indices
        EXTRA_ARGUMENTS) {
    float energy = 0.0f;
    for (int index = get_global_id(0); index < numBonds; index += get_global_size(0)) {
        // Look up the data for this bond.

        int4 atoms = indices[index];
        float4 delta = posq[atoms.y]-posq[atoms.x];

        // Compute the force.

        float r = native_sqrt(delta.x*delta.x + delta.y*delta.y + delta.z*delta.z);
        COMPUTE_FORCE
        delta.xyz *= -dEdR/r;

        // Record the force on each of the two atoms.

        unsigned int offsetA = atoms.x+atoms.z*numAtoms;
        unsigned int offsetB = atoms.y+atoms.w*numAtoms;
        float4 forceA = forceBuffers[offsetA];
        float4 forceB = forceBuffers[offsetB];
        forceA.xyz -= delta.xyz;
        forceB.xyz += delta.xyz;
        forceBuffers[offsetA] = forceA;
        forceBuffers[offsetB] = forceB;
    }
    energyBuffer[get_global_id(0)] += energy;
}
