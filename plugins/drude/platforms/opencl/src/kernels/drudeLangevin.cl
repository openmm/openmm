/**
 * Perform the first step of Langevin integration.
 */

__kernel void integrateDrudeLangevinPart1(__global mixed4* restrict velm, __global const real4* restrict force, __global mixed4* restrict posDelta,
        __global const int* restrict normalParticles, __global const int2* restrict pairParticles, __global const mixed2* restrict dt, mixed vscale, mixed fscale,
        mixed noisescale, mixed vscaleDrude, mixed fscaleDrude, mixed noisescaleDrude, __global const float4* restrict random, unsigned int randomIndex) {
    mixed stepSize = dt[0].y;
    
    // Update normal particles.

    for (int i = get_global_id(0); i < NUM_NORMAL_PARTICLES; i += get_global_size(0)) {
        int index = normalParticles[i];
        mixed4 velocity = velm[index];
        if (velocity.w != 0) {
            mixed sqrtInvMass = SQRT(velocity.w);
            float4 rand = random[randomIndex+index];
            real4 f = force[index];
            velocity.xyz = vscale*velocity.xyz + fscale*velocity.w*f.xyz + noisescale*sqrtInvMass*rand.xyz;
            velm[index] = velocity;
            posDelta[index] = (mixed4) (stepSize*velocity.x, stepSize*velocity.y, stepSize*velocity.z, 0);
        }
    }
    
    // Update Drude particle pairs.
    
    randomIndex += NUM_NORMAL_PARTICLES;
    for (int i = get_global_id(0); i < NUM_PAIRS; i += get_global_size(0)) {
        int2 particles = pairParticles[i];
        mixed4 velocity1 = velm[particles.x];
        mixed4 velocity2 = velm[particles.y];
        mixed mass1 = RECIP(velocity1.w);
        mixed mass2 = RECIP(velocity2.w);
        mixed invTotalMass = RECIP(mass1+mass2);
        mixed invReducedMass = (mass1+mass2)*velocity1.w*velocity2.w;
        mixed mass1fract = invTotalMass*mass1;
        mixed mass2fract = invTotalMass*mass2;
        mixed sqrtInvTotalMass = SQRT(invTotalMass);
        mixed sqrtInvReducedMass = SQRT(invReducedMass);
        mixed4 cmVel = velocity1*mass1fract+velocity2*mass2fract;
        mixed4 relVel = velocity2-velocity1;
        mixed4 force1 = force[particles.x];
        mixed4 force2 = force[particles.y];
        mixed4 cmForce = force1+force2;
        mixed4 relForce = force2*mass1fract - force1*mass2fract;
        float4 rand1 = random[randomIndex+2*i];
        float4 rand2 = random[randomIndex+2*i+1];
        cmVel.xyz = vscale*cmVel.xyz + fscale*invTotalMass*cmForce.xyz + noisescale*sqrtInvTotalMass*rand1.xyz;
        relVel.xyz = vscaleDrude*relVel.xyz + fscaleDrude*invReducedMass*relForce.xyz + noisescaleDrude*sqrtInvReducedMass*rand2.xyz;
        velocity1.xyz = cmVel.xyz-relVel.xyz*mass2fract;
        velocity2.xyz = cmVel.xyz+relVel.xyz*mass1fract;
        velm[particles.x] = velocity1;
        velm[particles.y] = velocity2;
        posDelta[particles.x] = (mixed4) (stepSize*velocity1.x, stepSize*velocity1.y, stepSize*velocity1.z, 0);
        posDelta[particles.y] = (mixed4) (stepSize*velocity2.x, stepSize*velocity2.y, stepSize*velocity2.z, 0);
    }
}

/**
 * Perform the second step of Langevin integration.
 */

__kernel void integrateDrudeLangevinPart2(__global real4* restrict posq, __global real4* restrict posqCorrection, __global const mixed4* restrict posDelta, __global mixed4* restrict velm, __global const mixed2* restrict dt) {
#ifdef SUPPORTS_DOUBLE_PRECISION
    double invStepSize = 1.0/dt[0].y;
#else
    float invStepSize = 1.0f/dt[0].y;
#endif
    int index = get_global_id(0);
    while (index < NUM_ATOMS) {
        mixed4 vel = velm[index];
        if (vel.w != 0.0) {
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
            vel.xyz = convert_mixed4(invStepSize*convert_double4(delta)).xyz;
#else
            vel.xyz = invStepSize*delta.xyz;
#endif
#ifdef USE_MIXED_PRECISION
            posq[index] = convert_real4(pos);
            posqCorrection[index] = (real4) (pos.x-(real) pos.x, pos.y-(real) pos.y, pos.z-(real) pos.z, 0);
#else
            posq[index] = pos;
#endif
            velm[index] = vel;
        }
        index += get_global_size(0);
    }
}
