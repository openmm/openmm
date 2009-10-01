/**
 * Compute nonbonded exceptions.
 */

__kernel void computeNonbondedExceptions(int numAtoms, int numExceptions, float cutoffSquared, float4 periodicBoxSize, __global float4* forceBuffers, __global float* energyBuffer,
        __global float4* posq, __global float4* params, __global int4* indices) {
    int index = get_global_id(0);
    float energy = 0.0f;
    while (index < numExceptions) {
        // Look up the data for this bonds.

        int4 atoms = indices[index];
        float4 exceptionParams = params[index];
        float4 delta = posq[atoms.y]-posq[atoms.x];
#ifdef USE_PERIODIC
        delta.x -= floor(delta.x/periodicBoxSize.x+0.5f)*periodicBoxSize.x;
        delta.y -= floor(delta.y/periodicBoxSize.y+0.5f)*periodicBoxSize.y;
        delta.z -= floor(delta.z/periodicBoxSize.z+0.5f)*periodicBoxSize.z;
#endif

        // Compute the force.

        float r2 = delta.x*delta.x + delta.y*delta.y + delta.z*delta.z;
        float invR = 1.0f/sqrt(r2);
        float sig2 = invR*exceptionParams.y;
        sig2 *= sig2;
        float sig6 = sig2*sig2*sig2;
        float dEdR = exceptionParams.z*(12.0f*sig6-6.0f)*sig6;
        float tempEnergy = exceptionParams.z*(sig6-1.0f)*sig6;
#ifdef USE_CUTOFF
        dEdR += exceptionParams.x*(invR-2.0f*cSim.reactionFieldK*r2);
        tempEnergy += exceptionParams.x*(invR+cSim.reactionFieldK*r2-cSim.reactionFieldC);
#else
        dEdR += exceptionParams.x*invR;
        tempEnergy += exceptionParams.x*invR;
#endif
        dEdR *= invR*invR;
#ifdef USE_CUTOFF
        if (r2 > cutoffSquared) {
            dEdR = 0.0f;
            tempEnergy  = 0.0f;
        }
#endif
        energy += tempEnergy;
        delta.xyz *= dEdR;

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
