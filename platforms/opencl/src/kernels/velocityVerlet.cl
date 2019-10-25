/**
 * Perform the first step of Velocity Verlet integration.
 */

__kernel void integrateVelocityVerletPart1(int numAtoms, int paddedNumAtoms, __global const mixed2* restrict dt, __global const real4* restrict posq,
        __global const real4* restrict posqCorrection, __global mixed4* restrict velm, __global const real4* restrict force, __global mixed4* restrict posDelta) {
    const mixed2 stepSize = dt[0];
    const mixed dtPos = stepSize.y;
    const mixed dtVel = 0.5f*(stepSize.x+stepSize.y);
    int index = get_global_id(0);
    while (index < numAtoms) { 
        mixed4 velocity = velm[index];
        if (velocity.w != 0.0) {
#ifdef USE_MIXED_PRECISION
            real4 pos1 = posq[index];
            real4 pos2 = posqCorrection[index];
            mixed4 pos = (mixed4) (pos1.x+(mixed)pos2.x, pos1.y+(mixed)pos2.y, pos1.z+(mixed)pos2.z, pos1.w);
#else
            real4 pos = posq[index];
#endif
            velocity.x += 0.5f * dtVel*force[index].x*velocity.w;
            velocity.y += 0.5f * dtVel*force[index].y*velocity.w;
            velocity.z += 0.5f * dtVel*force[index].z*velocity.w;
            pos.x = velocity.x*dtPos;
            pos.y = velocity.y*dtPos;
            pos.z = velocity.z*dtPos;
            posDelta[index] = pos;
            velm[index] = velocity;
        }
        index += get_global_size(0);
    }
}

/**
 * Perform the second step of Velocity Verlet integration.
 */

__kernel void integrateVelocityVerletPart2(int numAtoms, __global mixed2* restrict dt, __global real4* restrict posq,
        __global real4* restrict posqCorrection, __global mixed4* restrict velm, __global const mixed4* restrict posDelta) {
    mixed2 stepSize = dt[0];
    int index = get_global_id(0);
    if (index == 0)
        dt[0].x = stepSize.y;
    while(index < numAtoms) {
        mixed4 velocity = velm[index];
        if (velocity.w != 0.0) {
#ifdef USE_MIXED_PRECISION
            real4 pos1 = posq[index];
            real4 pos2 = posqCorrection[index];
            mixed4 pos = (mixed4) (pos1.x+(mixed)pos2.x, pos1.y+(mixed)pos2.y, pos1.z+(mixed)pos2.z, pos1.w);
#else
            real4 pos = posq[index];
#endif
            mixed4 delta = posDelta[index];
            pos.xyz += delta.xyz;
#ifdef USE_MIXED_PRECISION
            posq[index] = (real4) ((real) pos.x, (real) pos.y, (real) pos.z, (real) pos.w);
            posqCorrection[index] = (real4) (pos.x-(real) pos.x, pos.y-(real) pos.y, pos.z-(real) pos.z, 0);
#else
            posq[index] = pos;
#endif
        }
        index += get_global_size(0);
    }
}

/**
 * Perform the third step of Velocity Verlet integration.
 */

__kernel void integrateVelocityVerletPart3(int numAtoms, int paddedNumAtoms, __global mixed2* restrict dt, __global real4* restrict posq,
        __global real4* restrict posqCorrection, __global mixed4* restrict velm, __global const real4* restrict force, __global const mixed4* restrict posDelta) {
    mixed2 stepSize = dt[0];
#ifndef SUPPORTS_DOUBLE_PRECISION
    double oneOverDt = 1.0/stepSize.y;
#else
    float oneOverDt = 1.0f/stepSize.y;
    float correction = (1.0f-oneOverDt*stepSize.y)/stepSize.y;
#endif
    const mixed dtVel = 0.5f*(stepSize.x+stepSize.y);
    int index = get_global_id(0);
    if (index == 0)
        dt[0].x = stepSize.y;
    while(index < numAtoms) {
        mixed4 velocity = velm[index];
        if (velocity.w != 0.0) {
            mixed4 deltaXconstrained = posDelta[index];
            velocity.x += 0.5f * dtVel*force[index].x*velocity.w + (deltaXconstrained.x - velocity.x*stepSize.y)*oneOverDt;
            velocity.y += 0.5f * dtVel*force[index].y*velocity.w + (deltaXconstrained.y - velocity.y*stepSize.y)*oneOverDt;
            velocity.z += 0.5f * dtVel*force[index].z*velocity.w + (deltaXconstrained.z - velocity.z*stepSize.y)*oneOverDt;
#ifdef SUPPORTS_DOUBLE_PRECISION
            velocity.x += (deltaXconstrained.x - velocity.x*stepSize.y)*correction;
            velocity.y += (deltaXconstrained.y - velocity.y*stepSize.y)*correction;
            velocity.z += (deltaXconstrained.z - velocity.z*stepSize.y)*correction;
#endif
            velm[index] = velocity;
        }
        index += get_global_size(0);
    }
}

