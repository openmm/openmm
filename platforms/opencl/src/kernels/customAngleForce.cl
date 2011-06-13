/**
 * Compute custom angle forces.
 */

__kernel void computeCustomAngleForces(int numAtoms, int numAngles, __global float4* forceBuffers, __global float* energyBuffer,
        __global float4* posq, __global int8* indices
        EXTRA_ARGUMENTS) {
    float energy = 0.0f;
    for (int index = get_global_id(0); index < numAngles; index += get_global_size(0)) {
        int8 atoms = indices[index];
        float4 a1 = posq[atoms.s0];
        float4 a2 = posq[atoms.s1];
        float4 a3 = posq[atoms.s2];

        // Compute the force.

        float4 v0 = a2-a1;
        float4 v1 = a2-a3;
        float4 cp = cross(v0, v1);
        float rp = cp.x*cp.x + cp.y*cp.y + cp.z*cp.z;
        rp = max(sqrt(rp), 1.0e-06f);
        float r21 = v0.x*v0.x + v0.y*v0.y + v0.z*v0.z;
        float r23 = v1.x*v1.x + v1.y*v1.y + v1.z*v1.z;
        float dot = v0.x*v1.x + v0.y*v1.y + v0.z*v1.z;
        float cosine = clamp(dot/sqrt(r21*r23), -1.0f, 1.0f);
        float theta = acos(cosine);
        COMPUTE_FORCE
        float4 c21 = cross(v0, cp)*(dEdAngle/(r21*rp));
        float4 c23 = cross(cp, v1)*(dEdAngle/(r23*rp));

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
    }
    energyBuffer[get_global_id(0)] += energy;
}
