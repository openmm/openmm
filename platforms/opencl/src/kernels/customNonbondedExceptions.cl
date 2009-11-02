/**
 * Compute custom nonbonded exceptions.
 */

__kernel void computeCustomNonbondedExceptions(int numAtoms, int numExceptions, __global float4* forceBuffers, __global float* energyBuffer,
        __global float4* posq, __global float4* params, __global int4* indices
        EXTRA_ARGUMENTS) {
    int index = get_global_id(0);
    float energy = 0.0f;
    while (index < numExceptions) {
        // Look up the data for this exception.

        int4 atoms = indices[index];
        float4 exceptionParams = params[index];
        float4 delta = posq[atoms.y]-posq[atoms.x];
#ifdef USE_PERIODIC
        delta.x -= floor(delta.x/PERIODIC_BOX_SIZE_X+0.5f)*PERIODIC_BOX_SIZE_X;
        delta.y -= floor(delta.y/PERIODIC_BOX_SIZE_Y+0.5f)*PERIODIC_BOX_SIZE_Y;
        delta.z -= floor(delta.z/PERIODIC_BOX_SIZE_Z+0.5f)*PERIODIC_BOX_SIZE_Z;
#endif
        // Compute the force.

        float r2 = delta.x*delta.x + delta.y*delta.y + delta.z*delta.z;
#ifdef USE_CUTOFF
        if (r2 > CUTOFF_SQUARED) {
#else
        {
#endif
            float r = sqrt(r2);
            float dEdR;
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
        index += get_global_size(0);
    }
    energyBuffer[get_global_id(0)] += energy;
}
