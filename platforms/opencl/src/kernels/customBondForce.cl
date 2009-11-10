/**
 * Compute custom bond forces.
 */

__kernel void computeCustomBondForces(int numAtoms, int numBonds, __global float4* forceBuffers, __global float* energyBuffer,
        __global float4* posq, __global float4* params, __global int4* indices
        EXTRA_ARGUMENTS) {
    int index = get_global_id(0);
    float energy = 0.0f;
    while (index < numBonds) {
        // Look up the data for this exception.

        int4 atoms = indices[index];
        float4 exceptionParams = params[index];
        float4 delta = posq[atoms.y]-posq[atoms.x];

        // Compute the force.

        float r = sqrt(delta.x*delta.x + delta.y*delta.y + delta.z*delta.z);
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
        index += get_global_size(0);
    }
    energyBuffer[get_global_id(0)] += energy;
}
