/**
 * Evaluate the forces due to harmonic bonds.
 */

__kernel void calcHarmonicBondForce(int numAtoms, int numBonds, __global float4* forceBuffers, __global float* energyBuffer, __global float4* posq, __global float2* params, __global int4* indices) {
    int index = get_global_id(0);
    float energy = 0.0f;
    while (index < numBonds) {
        // Look up the data for this bonds.

        int4 atoms = indices[index];
        float4 delta = posq[atoms.y]-posq[atoms.x];
        float2 bondParams = params[index];

        // Compute the force.

        float r = SQRT(delta.x*delta.x + delta.y*delta.y + delta.z*delta.z);
        float deltaIdeal = r-bondParams.x;
        energy += 0.5f * bondParams.y*deltaIdeal*deltaIdeal;
        float dEdR = bondParams.y * deltaIdeal;
        dEdR = (r > 0.0f) ? (dEdR / r) : 0.0f;
        delta.xyz *= dEdR;

        // Record the force on each of the two atoms.

        unsigned int offsetA = atoms.x+atoms.z*numAtoms;
        unsigned int offsetB = atoms.y+atoms.w*numAtoms;
        float4 forceA = forceBuffers[offsetA];
        float4 forceB = forceBuffers[offsetB];
        forceA.xyz += delta.xyz;
        forceB.xyz -= delta.xyz;
        forceBuffers[offsetA] = forceA;
        forceBuffers[offsetB] = forceB;
        index += get_global_size(0);
    }
    energyBuffer[get_global_id(0)] += energy;
}
