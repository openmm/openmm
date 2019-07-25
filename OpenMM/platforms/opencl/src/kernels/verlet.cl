/**
 * Perform the first step of verlet integration.
 */

__kernel void integrateVerletPart1(int numAtoms, __global const mixed2* restrict dt, __global const real4* restrict posq, __global const real4* restrict posqCorrection, __global mixed4* restrict velm, __global const real4* restrict force, __global mixed4* restrict posDelta) {
    mixed2 stepSize = dt[0];
    mixed dtPos = stepSize.y;
    mixed dtVel = 0.5f*(stepSize.x+stepSize.y);
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
            velocity.x += force[index].x*dtVel*velocity.w;
            velocity.y += force[index].y*dtVel*velocity.w;
            velocity.z += force[index].z*dtVel*velocity.w;
            pos.xyz = velocity.xyz*dtPos;
            posDelta[index] = pos;
            velm[index] = velocity;
        }
        index += get_global_size(0);
    }
}

/**
 * Perform the second step of verlet integration.
 */

__kernel void integrateVerletPart2(int numAtoms, __global mixed2* restrict dt, __global real4* restrict posq, __global real4* restrict posqCorrection, __global mixed4* restrict velm, __global const mixed4* restrict posDelta) {
    mixed2 stepSize = dt[0];
#ifdef SUPPORTS_DOUBLE_PRECISION
    double oneOverDt = 1.0/stepSize.y;
#else
    float oneOverDt = 1.0f/stepSize.y;
    float correction = (1.0f-oneOverDt*stepSize.y)/stepSize.y;
#endif
    if (get_global_id(0) == 0)
        dt[0].x = stepSize.y;
    barrier(CLK_LOCAL_MEM_FENCE);
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
            mixed4 delta = posDelta[index];
            pos.xyz += delta.xyz;
#ifdef SUPPORTS_DOUBLE_PRECISION
            velocity.xyz = convert_mixed4(convert_double4(delta)*oneOverDt).xyz;
#else
            velocity.xyz = delta.xyz*oneOverDt + delta.xyz*correction;
#endif
#ifdef USE_MIXED_PRECISION
            posq[index] = convert_real4(pos);
            posqCorrection[index] = (real4) (pos.x-(real) pos.x, pos.y-(real) pos.y, pos.z-(real) pos.z, 0);
#else
            posq[index] = pos;
#endif
            velm[index] = velocity;
        }
        index += get_global_size(0);
    }
}

/**
 * Select the step size to use for the next step.
 */

__kernel void selectVerletStepSize(int numAtoms, mixed maxStepSize, mixed errorTol, __global mixed2* restrict dt, __global const mixed4* restrict velm, __global const real4* restrict force, __local mixed* restrict error) {
    // Calculate the error.

    mixed err = 0;
    int index = get_local_id(0);
    while (index < numAtoms) {
        real4 f = force[index];
        mixed invMass = velm[index].w;
        err += (f.x*f.x + f.y*f.y + f.z*f.z)*invMass*invMass;
        index += get_global_size(0);
    }
    error[get_local_id(0)] = err;
    barrier(CLK_LOCAL_MEM_FENCE);

    // Sum the errors from all threads.

    for (unsigned int offset = 1; offset < get_local_size(0); offset *= 2) {
        if (get_local_id(0)+offset < get_local_size(0) && (get_local_id(0)&(2*offset-1)) == 0)
            error[get_local_id(0)] += error[get_local_id(0)+offset];
        barrier(CLK_LOCAL_MEM_FENCE);
    }
    if (get_local_id(0) == 0) {
        mixed totalError = sqrt(error[0]/(numAtoms*3));
        mixed newStepSize = sqrt(errorTol/totalError);
        mixed oldStepSize = dt[0].y;
        if (oldStepSize > 0.0f)
            newStepSize = min(newStepSize, oldStepSize*2.0f); // For safety, limit how quickly dt can increase.
        if (newStepSize > oldStepSize && newStepSize < 1.1f*oldStepSize)
            newStepSize = oldStepSize; // Keeping dt constant between steps improves the behavior of the integrator.
        if (newStepSize > maxStepSize)
            newStepSize = maxStepSize;
        dt[0].y = newStepSize;
    }
}
