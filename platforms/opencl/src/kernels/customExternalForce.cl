/**
 * Compute custom external forces.
 */

__kernel void computeCustomExternalForces(int numTerms, __global float4* forceBuffers, __global float* energyBuffer,
        __global float4* posq, __global float4* params, __global int* indices
        EXTRA_ARGUMENTS) {
    float energy = 0.0f;
    for (int index = get_global_id(0); index < numTerms; index += get_global_size(0)) {
        // Look up the data for this particle.

        int atom = indices[index];
        float4 particleParams = params[index];
        float4 pos = posq[atom];

        // Compute the force.

        COMPUTE_FORCE

        // Record the force on the atom.

        float4 force = forceBuffers[atom];
        force.x -= dEdX;
        force.y -= dEdY;
        force.z -= dEdZ;
        forceBuffers[atom] = force;
    }
    energyBuffer[get_global_id(0)] += energy;
}
