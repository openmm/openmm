/**
 * Compute CNAP torsion forces.
 */

__kernel void computeCMAPTorsionForces(int numAtoms, int numTorsions, __global float4* forceBuffers, __global float* energyBuffer,
        __global float4* posq, __global float4* coeff, __global int2* mapPositions, __global int16* indices, __global int* maps) {
    const float PI = 3.14159265358979323846f;
    float energy = 0.0f;
    for (int index = get_global_id(0); index < numTorsions; index += get_global_size(0)) {
        int16 atoms = indices[index];
        float4 a1 = posq[atoms.s0];
        float4 a2 = posq[atoms.s1];
        float4 a3 = posq[atoms.s2];
        float4 a4 = posq[atoms.s3];
        float4 b1 = posq[atoms.s4];
        float4 b2 = posq[atoms.s5];
        float4 b3 = posq[atoms.s6];
        float4 b4 = posq[atoms.s7];

        // Compute the first angle.

        float4 v0a = (float4) (a1.xyz-a2.xyz, 0.0f);
        float4 v1a = (float4) (a3.xyz-a2.xyz, 0.0f);
        float4 v2a = (float4) (a3.xyz-a4.xyz, 0.0f);
        float4 cp0a = cross(v0a, v1a);
        float4 cp1a = cross(v1a, v2a);
        float cosangle = dot(normalize(cp0a), normalize(cp1a));
        float angleA;
        if (cosangle > 0.99f || cosangle < -0.99f) {
            // We're close to the singularity in acos(), so take the cross product and use asin() instead.

            float4 cross_prod = cross(cp0a, cp1a);
            float scale = dot(cp0a, cp0a)*dot(cp1a, cp1a);
            angleA = asin(SQRT(dot(cross_prod, cross_prod)/scale));
            if (cosangle < 0.0f)
                angleA = PI-angleA;
        }
        else
           angleA = acos(cosangle);
        angleA = (dot(v0a, cp1a) >= 0 ? angleA : -angleA);
        angleA = fmod(angleA+2.0f*PI, 2.0f*PI);

        // Compute the second angle.

        float4 v0b = (float4) (b1.xyz-b2.xyz, 0.0f);
        float4 v1b = (float4) (b3.xyz-b2.xyz, 0.0f);
        float4 v2b = (float4) (b3.xyz-b4.xyz, 0.0f);
        float4 cp0b = cross(v0b, v1b);
        float4 cp1b = cross(v1b, v2b);
        cosangle = dot(normalize(cp0b), normalize(cp1b));
        float angleB;
        if (cosangle > 0.99f || cosangle < -0.99f) {
            // We're close to the singularity in acos(), so take the cross product and use asin() instead.

            float4 cross_prod = cross(cp0b, cp1b);
            float scale = dot(cp0b, cp0b)*dot(cp1b, cp1b);
            angleB = asin(SQRT(dot(cross_prod, cross_prod)/scale));
            if (cosangle < 0.0f)
                angleB = PI-angleB;
        }
        else
           angleB = acos(cosangle);
        angleB = (dot(v0b, cp1b) >= 0 ? angleB : -angleB);
        angleB = fmod(angleB+2.0f*PI, 2.0f*PI);

        // Identify which patch this is in.

        int2 pos = mapPositions[maps[index]];
        int size = pos.y;
        float delta = 2*PI/size;
        int s = (int) (angleA/delta);
        int t = (int) (angleB/delta);
        float4 c[4];
        int coeffIndex = pos.x+4*(s+size*t);
        c[0] = coeff[coeffIndex];
        c[1] = coeff[coeffIndex+1];
        c[2] = coeff[coeffIndex+2];
        c[3] = coeff[coeffIndex+3];
        float da = angleA/delta-s;
        float db = angleB/delta-t;

        // Evaluate the spline to determine the energy and gradients.

        float torsionEnergy = 0.0f;
        float dEdA = 0.0f;
        float dEdB = 0.0f;
        torsionEnergy = da*torsionEnergy + ((c[3].w*db + c[3].z)*db + c[3].y)*db + c[3].x;
        dEdA = db*dEdA + (3.0f*c[3].w*da + 2.0f*c[2].w)*da + c[1].w;
        dEdB = da*dEdB + (3.0f*c[3].w*db + 2.0f*c[3].z)*db + c[3].y;
        torsionEnergy = da*torsionEnergy + ((c[2].w*db + c[2].z)*db + c[2].y)*db + c[2].x;
        dEdA = db*dEdA + (3.0f*c[3].z*da + 2.0f*c[2].z)*da + c[1].z;
        dEdB = da*dEdB + (3.0f*c[2].w*db + 2.0f*c[2].z)*db + c[2].y;
        torsionEnergy = da*torsionEnergy + ((c[1].w*db + c[1].z)*db + c[1].y)*db + c[1].x;
        dEdA = db*dEdA + (3.0f*c[3].y*da + 2.0f*c[2].y)*da + c[1].y;
        dEdB = da*dEdB + (3.0f*c[1].w*db + 2.0f*c[1].z)*db + c[1].y;
        torsionEnergy = da*torsionEnergy + ((c[0].w*db + c[0].z)*db + c[0].y)*db + c[0].x;
        dEdA = db*dEdA + (3.0f*c[3].x*da + 2.0f*c[2].x)*da + c[1].x;
        dEdB = da*dEdB + (3.0f*c[0].w*db + 2.0f*c[0].z)*db + c[0].y;
        dEdA /= delta;
        dEdB /= delta;
        energy += torsionEnergy;

        // Apply the force to the first torsion.

        float normCross1 = dot(cp0a, cp0a);
        float normSqrBC = dot(v1a, v1a);
        float normBC = SQRT(normSqrBC);
        float normCross2 = dot(cp1a, cp1a);
        float dp = 1.0f/normSqrBC;
        float4 ff = (float4) ((-dEdA*normBC)/normCross1, dot(v0a, v1a)*dp, dot(v2a, v1a)*dp, (dEdA*normBC)/normCross2);
        float4 internalF0 = ff.x*cp0a;
        float4 internalF3 = ff.w*cp1a;
        float4 d = ff.y*internalF0 - ff.z*internalF3;
        int4 offset = atoms.lo.lo + numAtoms*atoms.hi.lo;
        float4 forceA = forceBuffers[offset.x];
        float4 forceB = forceBuffers[offset.y];
        float4 forceC = forceBuffers[offset.z];
        float4 forceD = forceBuffers[offset.w];
        forceA.xyz += internalF0.xyz;
        forceB.xyz += d.xyz-internalF0.xyz;
        forceC.xyz += -d.xyz-internalF3.xyz;
        forceD.xyz += internalF3.xyz;
        forceBuffers[offset.x] = forceA;
        forceBuffers[offset.y] = forceB;
        forceBuffers[offset.z] = forceC;
        forceBuffers[offset.w] = forceD;

        // Apply the force to the second torsion.

        normCross1 = dot(cp0b, cp0b);
        normSqrBC = dot(v1b, v1b);
        normBC = SQRT(normSqrBC);
        normCross2 = dot(cp1b, cp1b);
        dp = 1.0f/normSqrBC;
        ff = (float4) ((-dEdB*normBC)/normCross1, dot(v0b, v1b)*dp, dot(v2b, v1b)*dp, (dEdB*normBC)/normCross2);
        internalF0 = ff.x*cp0b;
        internalF3 = ff.w*cp1b;
        d = ff.y*internalF0 - ff.z*internalF3;
        offset = atoms.lo.hi + numAtoms*atoms.hi.hi;
        forceA = forceBuffers[offset.x];
        forceB = forceBuffers[offset.y];
        forceC = forceBuffers[offset.z];
        forceD = forceBuffers[offset.w];
        forceA.xyz += internalF0.xyz;
        forceB.xyz += d.xyz-internalF0.xyz;
        forceC.xyz += -d.xyz-internalF3.xyz;
        forceD.xyz += internalF3.xyz;
        forceBuffers[offset.x] = forceA;
        forceBuffers[offset.y] = forceB;
        forceBuffers[offset.z] = forceC;
        forceBuffers[offset.w] = forceD;
    }
    energyBuffer[get_global_id(0)] += energy;
}
