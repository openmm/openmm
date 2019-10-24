enum {VelScale, NoiseScale};

/**
 * Perform the first part of BAOAB integration: velocity half step, then position half step.
 */

__kernel void integrateBAOABPart1(__global mixed4* restrict velm, __global const real4* restrict force, __global mixed4* restrict posDelta,
        __global mixed4* restrict oldDelta, __global const mixed2* restrict dt) {
    mixed halfdt = 0.5*dt[0].y;
    for (int index = get_global_id(0); index < NUM_ATOMS; index += get_global_size(0)) {
        mixed4 velocity = velm[index];
        if (velocity.w != 0.0) {
            velocity.x += halfdt*velocity.w*force[index].x;
            velocity.y += halfdt*velocity.w*force[index].y;
            velocity.z += halfdt*velocity.w*force[index].z;
            velm[index] = velocity;
            mixed4 delta = halfdt*velocity;
            posDelta[index] = delta;
            oldDelta[index] = delta;
        }
    }
}

/**
 * Perform the second part of BAOAB integration: apply constraint forces to velocities, then interact with heat bath,
 * then position half step.
 */

__kernel void integrateBAOABPart2(__global real4* restrict posq, __global real4* restrict posqCorrection, __global mixed4* restrict velm, __global mixed4* restrict posDelta,
        __global mixed4* restrict oldDelta, __global const mixed* restrict paramBuffer, __global const mixed2* restrict dt, __global const float4* restrict random, unsigned int randomIndex) {
    mixed vscale = paramBuffer[VelScale];
    mixed noisescale = paramBuffer[NoiseScale];
    mixed halfdt = 0.5*dt[0].y;
    mixed invHalfdt = 1/halfdt;
    int index = get_global_id(0);
    randomIndex += index;
    while (index < NUM_ATOMS) {
        mixed4 velocity = velm[index];
        if (velocity.w != 0.0) {
            mixed4 delta = posDelta[index];
            mixed sqrtInvMass = SQRT(velocity.w);
            velocity.xyz += (delta.xyz-oldDelta[index].xyz)*invHalfdt;
            velocity.x = vscale*velocity.x + noisescale*sqrtInvMass*random[randomIndex].x;
            velocity.y = vscale*velocity.y + noisescale*sqrtInvMass*random[randomIndex].y;
            velocity.z = vscale*velocity.z + noisescale*sqrtInvMass*random[randomIndex].z;
            velm[index] = velocity;
#ifdef USE_MIXED_PRECISION
            real4 pos1 = posq[index];
            real4 pos2 = posqCorrection[index];
            mixed4 pos = (mixed4) (pos1.x+(mixed)pos2.x, pos1.y+(mixed)pos2.y, pos1.z+(mixed)pos2.z, pos1.w);
#else
            real4 pos = posq[index];
#endif
            pos.xyz += delta.xyz;
#ifdef USE_MIXED_PRECISION
            posq[index] = convert_real4(pos);
            posqCorrection[index] = (real4) (pos.x-(real) pos.x, pos.y-(real) pos.y, pos.z-(real) pos.z, 0);
#else
            posq[index] = pos;
#endif
            delta = halfdt*velocity;
            posDelta[index] = delta;
            oldDelta[index] = delta;
        }
        randomIndex += get_global_size(0);
        index += get_global_size(0);
    }
}

/**
 * Perform the third part of BAOAB integration: apply constraint forces to velocities, then record
 * the constrained positions in preparation for computing forces.
 */

__kernel void integrateBAOABPart3(__global real4* restrict posq, __global real4* restrict posqCorrection, __global mixed4* restrict velm,
         __global mixed4* restrict posDelta, __global mixed4* restrict oldDelta, __global const mixed2* restrict dt) {
    mixed halfdt = 0.5*dt[0].y;
    mixed invHalfdt = 1/halfdt;
    for (int index = get_global_id(0); index < NUM_ATOMS; index += get_global_size(0)) {
        mixed4 velocity = velm[index];
        if (velocity.w != 0.0) {
            mixed4 delta = posDelta[index];
            velocity.x += (delta.x-oldDelta[index].x)*invHalfdt;
            velocity.y += (delta.y-oldDelta[index].y)*invHalfdt;
            velocity.z += (delta.z-oldDelta[index].z)*invHalfdt;
            velm[index] = velocity;
#ifdef USE_MIXED_PRECISION
            real4 pos1 = posq[index];
            real4 pos2 = posqCorrection[index];
            mixed4 pos = (mixed4) (pos1.x+(mixed)pos2.x, pos1.y+(mixed)pos2.y, pos1.z+(mixed)pos2.z, pos1.w);
#else
            real4 pos = posq[index];
#endif
            pos.xyz += delta.xyz;
#ifdef USE_MIXED_PRECISION
            posq[index] = convert_real4(pos);
            posqCorrection[index] = (real4) (pos.x-(real) pos.x, pos.y-(real) pos.y, pos.z-(real) pos.z, 0);
#else
            posq[index] = pos;
#endif
        }
    }
}

/**
 * Perform the fourth part of BAOAB integration: velocity half step.
 */

__kernel void integrateBAOABPart4(__global mixed4* restrict velm, __global const real4* restrict force, __global const mixed2* restrict dt) {
    mixed halfdt = 0.5*dt[0].y;
    for (int index = get_global_id(0); index < NUM_ATOMS; index += get_global_size(0)) {
        mixed4 velocity = velm[index];
        if (velocity.w != 0.0) {
            velocity.x += halfdt*velocity.w*force[index].x;
            velocity.y += halfdt*velocity.w*force[index].y;
            velocity.z += halfdt*velocity.w*force[index].z;
            velm[index] = velocity;
        }
    }
}
