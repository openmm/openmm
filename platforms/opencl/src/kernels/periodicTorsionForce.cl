/**
 * Evaluate the forces due to periodic torsions.
 */

__kernel void calcPeriodicTorsionForce(int numAtoms, int numTorsions, __global float4* forceBuffers, __global float* energyBuffer, __global float4* posq, __global float4* params, __global int8* indices) {
    int index = get_global_id(0);
    float energy = 0.0f;
    while (index < numTorsions) {
        // Look up the data for this torsion.

        int8 atoms = indices[index];
        float4 torsionParams = params[index];
        float4 a1 = posq[atoms.x];
        float4 a2 = posq[atoms.y];
        float4 a3 = posq[atoms.z];
        float4 a4 = posq[atoms.w];

        // Compute the force.

        float4 v0 = (float4) (a1.xyz-a2.xyz, 0.0f);
        float4 v1 = (float4) (a3.xyz-a2.xyz, 0.0f);
        float4 v2 = (float4) (a3.xyz-a4.xyz, 0.0f);
        float4 cp0 = cross(v0, v1);
        float4 cp1 = cross(v1, v2);
        float cosangle = dot(normalize(cp0), normalize(cp1));
        cosangle = clamp(cosangle, -1.0f, 1.0f);
        float dihedralAngle = acos(cosangle);
        float deltaAngle = torsionParams.z*dihedralAngle-torsionParams.y;
        energy += torsionParams.x*(1.0f+cos(deltaAngle));
        float sinDeltaAngle = sin(deltaAngle);
        float dEdAngle = -torsionParams.x*torsionParams.z*sinDeltaAngle;
        float normCross1 = dot(cp0, cp0);
        float normBC = sqrt(dot(v1, v1));
        float normCross2 = dot(cp1, cp1);
        float dp = 1.0f/normBC;
        float4 ff = (float4) ((-dEdAngle*normBC)/normCross1, dot(v0, v1)*dp, dot(v2, v1)*dp, (dEdAngle*normBC)/normCross2);
        float4 internalF0 = ff.x*cp0;
        float4 internalF3 = ff.w*cp1;
        float4 s = ff.y*internalF0 - ff.z*internalF3;

        // Record the force on each of the four atoms.

        unsigned int offsetA = atoms.s0+atoms.s4*numAtoms;
        unsigned int offsetB = atoms.s1+atoms.s5*numAtoms;
        unsigned int offsetC = atoms.s2+atoms.s6*numAtoms;
        unsigned int offsetD = atoms.s3+atoms.s7*numAtoms;
        float4 forceA = forceBuffers[offsetA];
        float4 forceB = forceBuffers[offsetB];
        float4 forceC = forceBuffers[offsetC];
        float4 forceD = forceBuffers[offsetD];
        forceA.xyz += internalF0.xyz;
        forceB.xyz += s.xyz-internalF0.xyz;
        forceC.xyz += -s.xyz-internalF3.xyz;
        forceD.xyz += internalF3.xyz;
        forceBuffers[offsetA] = forceA;
        forceBuffers[offsetB] = forceB;
        forceBuffers[offsetC] = forceC;
        forceBuffers[offsetD] = forceD;
        index += get_global_size(0);
    }
    energyBuffer[get_global_id(0)] += energy;
}
