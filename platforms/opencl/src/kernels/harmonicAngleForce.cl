/**
 * Evaluate the forces due to harmonic angles.
 */

__kernel void calcHarmonicAngleForce(int numAtoms, int numAngles, __global float4* forceBuffers, __global float* energyBuffer, __global float4* posq, __global float2* params, __global int8* indices) {
    int index = get_global_id(0);
    float energy = 0.0f;
    while (index < numAngles) {
        // Look up the data for this angle.

        int8 atoms = indices[index];
        float2 angleParams = params[index];
        float4 a1 = posq[atoms.x];
        float4 a2 = posq[atoms.y];
        float4 a3 = posq[atoms.z];

        // Compute the force.

        float4 v0 = a2-a1;
        float4 v1 = a2-a3;
        float4 cp = cross(v0, v1);
        float rp = dot(cp.xyz, cp.xyz);
        rp = max(sqrt(rp), 1.0e-06f);
        float r21 = dot(v0.xyz, v0.xyz);
        float r23 = dot(v1.xyz, v1.xyz);
        float dot = dot(v0.xyz, v1.xyz);
        float cosine = dot/sqrt(r21*r23);
        float deltaIdeal = acos(cosine)-angleParams.x;
        energy += 0.5f*angleParams.y*deltaIdeal*deltaIdeal;
        float dEdR = angleParams.y*deltaIdeal;
        float4 c21 = cross(v0, cp)*(dEdR/(r21*rp));
        float4 c23 = cross(cp, v1)*(dEdR/(r23*rp));

        // Record the force on each of the three atoms.

        unsigned int offsetA = atoms.s0+atoms.s3*numAtoms;
        unsigned int offsetB = atoms.s1+atoms.s4*numAtoms;
        unsigned int offsetC = atoms.s2+atoms.s5*numAtoms;
        float4 forceA = forceBuffers[offsetA];
        float4 forceB = forceBuffers[offsetB];
        float4 forceC = forceBuffers[offsetC];
        forceA.xyz += c21.xyz;
        forceB.xyz -= c21.xyz+c23.xyz;
        forceC.xyz += c23.xyz;
        forceBuffers[offsetA] = forceA;
        forceBuffers[offsetB] = forceB;
        forceBuffers[offsetC] = forceC;
        index += get_global_size(0);
    }
    energyBuffer[get_global_id(0)] += energy;
}
