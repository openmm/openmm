/**
 * Compute the difference between two vectors, setting the fourth component to the squared magnitude.
 */
float4 delta(float4 vec1, float4 vec2) {
    float4 result = (float4) (vec1.x-vec2.x, vec1.y-vec2.y, vec1.z-vec2.z, 0.0f);
    result.w = result.x*result.x + result.y*result.y + result.z*result.z;
    return result;
}

/**
 * Compute the difference between two vectors, taking periodic boundary conditions into account
 * and setting the fourth component to the squared magnitude.
 */
float4 deltaPeriodic(float4 vec1, float4 vec2) {
    float4 result = (float4) (vec1.x-vec2.x, vec1.y-vec2.y, vec1.z-vec2.z, 0.0f);
#ifdef USE_PERIODIC
    result.x -= floor(result.x/PERIODIC_BOX_SIZE_X+0.5f)*PERIODIC_BOX_SIZE_X;
    result.y -= floor(result.y/PERIODIC_BOX_SIZE_Y+0.5f)*PERIODIC_BOX_SIZE_Y;
    result.z -= floor(result.z/PERIODIC_BOX_SIZE_Z+0.5f)*PERIODIC_BOX_SIZE_Z;
#endif
    result.w = result.x*result.x + result.y*result.y + result.z*result.z;
    return result;
}

/**
 * Compute the angle between two vectors.  The w component of each vector should contain the squared magnitude.
 */

float computeAngle(float4 vec1, float4 vec2) {
    float dot = vec1.x*vec2.x + vec1.y*vec2.y + vec1.z*vec2.z;
    float cosine = dot/sqrt(vec1.w*vec2.w);
    float angle;
    if (cosine > 0.99f || cosine < -0.99f) {
        // We're close to the singularity in acos(), so take the cross product and use asin() instead.

        float4 cross_prod = cross(vec1, vec2);
        float scale = vec1.w*vec2.w;
        angle = asin(sqrt(dot(cross_prod, cross_prod)/scale));
        if (cosine < 0.0f)
            angle = M_PI-angle;
    }
    else
       angle = acos(cosine);
    return angle;
}

/**
 * Compute hbond interactions.
 */

__kernel void computeHbonds(__global float4* forceBuffers, __global float* energyBuffer, __global float4* posq, /*__global unsigned int* exclusions,
        __global unsigned int* exclusionIndices, */__global int4* donorAtoms, __global int4* acceptorAtoms, __global int4* donorBufferIndices, __global int4* acceptorBufferIndices, __local float4* posBuffer, __local float4* deltaBuffer
        PARAMETER_ARGUMENTS) {
    float energy = 0.0f;
    unsigned int tgx = get_local_id(0) & (get_local_size(0)-1);
    unsigned int tbx = get_local_id(0) - tgx;
    float4 f1 = 0;
    float4 f2 = 0;
    for (int donorIndex = get_global_id(0); donorIndex < NUM_DONORS; donorIndex += get_global_size(0)) {
        // Load information about the donor this thread will compute forces on.

        int4 atoms = donorAtoms[donorIndex];
        float4 d1 = posq[atoms.x];
        float4 d2 = posq[atoms.y];
        float4 d3 = posq[atoms.z];
        float4 deltaD1D2 = delta(d1, d2);
        for (int acceptorStart = 0; acceptorStart < NUM_ACCEPTORS; acceptorStart += get_local_size(0)) {
            // Load the next block of acceptors into local memory.

            int blockSize = min((int) get_local_size(0), NUM_ACCEPTORS-acceptorStart);
            if (tgx < blockSize) {
                int4 atoms2 = acceptorAtoms[acceptorStart+tgx];
                float4 pos1 = posq[atoms2.x];
                float4 pos2 = posq[atoms2.y];
                float4 pos3 = posq[atoms2.z];
                posBuffer[get_local_id(0)] = pos1;
                deltaBuffer[get_local_id(0)] = delta(pos2, pos1);
            }
            barrier(CLK_LOCAL_MEM_FENCE);
            for (int index = 0; index < blockSize; index++) {
                // Compute the interaction between a donor and an acceptor.

                float4 a1 = posBuffer[index];
                float4 deltaD1A1 = deltaPeriodic(d1, a1);
#ifdef USE_CUTOFF
                if (deltaD1A1.w < CUTOFF_SQUARED) {
#endif
                    // Compute variables the force can depend on.

                    float r = sqrt(deltaD1A1.w);
                    float4 deltaA2A1 = deltaBuffer[index];
                    float theta = computeAngle(deltaD1A1, deltaD1D2);
                    float psi = computeAngle(deltaD1A1, deltaA2A1);
                    float4 cross1 = cross(deltaA2A1, deltaD1A1);
                    float4 cross2 = cross(deltaD1A1, deltaD1D2);
                    cross1.w = cross1.x*cross1.x + cross1.y*cross1.y + cross1.z*cross1.z;
                    cross2.w = cross2.x*cross2.x + cross2.y*cross2.y + cross2.z*cross2.z;
                    float chi = computeAngle(cross1, cross2);
                    chi = (dot(deltaA2A1, cross2) < 0 ? -chi : chi);
                    COMPUTE_FORCE

#ifdef INCLUDE_R
                    // Apply forces based on r.

                    f1.xyz -= (dEdR/r)*deltaD1A1.xyz;
#endif

#ifdef INCLUDE_THETA
                    // Apply forces based on theta.

                    float4 thetaCross = cross(deltaD1D2, deltaD1A1);
                    float lengthThetaCross = max(length(thetaCross), 1e-6f);
                    float4 deltaCross0 = cross(deltaD1D2, thetaCross)*dEdTheta/(deltaD1D2.w*lengthThetaCross);
                    float4 deltaCross2 = -cross(deltaD1A1, thetaCross)*dEdTheta/(deltaD1A1.w*lengthThetaCross);
                    float4 deltaCross1 = -(deltaCross0+deltaCross2);
                    f1.xyz += deltaCross1.xyz;
                    f2.xyz += deltaCross0.xyz;
#endif

#ifdef INCLUDE_PSI
                    // Apply forces based on psi.

                    float4 psiCross = cross(deltaA2A1, deltaD1A1);
                    float lengthPsiCross = max(length(psiCross), 1e-6f);
                    deltaCross0 = cross(deltaD1A1, psiCross)*dEdPsi/(deltaD1A1.w*lengthPsiCross);
//                    float4 deltaCross2 = -cross(deltaA2A1, psiCross)*dEdPsi/(deltaA2A1.w*lengthPsiCross);
//                    float4 deltaCross1 = -(deltaCross0+deltaCross2);
                    f1.xyz += deltaCross0.xyz;
#endif

#ifdef INCLUDE_CHI
                    // Apply forces based on chi.

                    float4 ff;
                    ff.x = (-dEdChi*r)/cross1.w;
                    ff.y = (deltaA2A1.x*deltaD1A1.x + deltaA2A1.y*deltaD1A1.y + deltaA2A1.z*deltaD1A1.z)/deltaD1A1.w;
                    ff.z = (deltaD1D2.x*deltaD1A1.x + deltaD1D2.y*deltaD1A1.y + deltaD1D2.z*deltaD1A1.z)/deltaD1A1.w;
                    ff.w = (dEdChi*r)/cross2.w;
                    float4 internalF0 = ff.x*cross1;
                    float4 internalF3 = ff.w*cross2;
                    float4 s = ff.y*internalF0 - ff.z*internalF3;
                    f1.xyz -= s.xyz+internalF3.xyz;
                    f2.xyz += internalF3.xyz;
#endif
#ifdef USE_CUTOFF
                }
#endif
            }
        }

        // Write results

        int4 bufferIndices = donorBufferIndices[donorIndex];
        unsigned int offset1 = atoms.x+bufferIndices.x*PADDED_NUM_ATOMS;
        unsigned int offset2 = atoms.y+bufferIndices.y*PADDED_NUM_ATOMS;
        float4 force1 = forceBuffers[offset1];
        float4 force2 = forceBuffers[offset2];
        force1.xyz += f1.xyz;
        force2.xyz += f2.xyz;
        forceBuffers[offset1] = force1;
        forceBuffers[offset2] = force2;
    }
    energyBuffer[get_global_id(0)] += energy;
}
