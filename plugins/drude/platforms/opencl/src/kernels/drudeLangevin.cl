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
            mixed sqrtInvMass = sqrt(velocity.w);
            float4 rand = random[randomIndex+index];
            real4 f = force[index];
            velocity.x = vscale*velocity.x + fscale*velocity.w*f.x + noisescale*sqrtInvMass*rand.x;
            velocity.y = vscale*velocity.y + fscale*velocity.w*f.y + noisescale*sqrtInvMass*rand.y;
            velocity.z = vscale*velocity.z + fscale*velocity.w*f.z + noisescale*sqrtInvMass*rand.z;
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
        mixed mass1 = 1/velocity1.w;
        mixed mass2 = 1/velocity2.w;
        mixed invTotalMass = 1/(mass1+mass2);
        mixed invReducedMass = (mass1+mass2)*velocity1.w*velocity2.w;
        mixed mass1fract = invTotalMass*mass1;
        mixed mass2fract = invTotalMass*mass2;
        mixed sqrtInvTotalMass = sqrt(invTotalMass);
        mixed sqrtInvReducedMass = sqrt(invReducedMass);
        mixed4 cmVel = velocity1*mass1fract+velocity2*mass2fract;
        mixed4 relVel = velocity2-velocity1;
        mixed4 force1 = convert_mixed4(force[particles.x]);
        mixed4 force2 = convert_mixed4(force[particles.y]);
        mixed4 cmForce = force1+force2;
        mixed4 relForce = force2*mass1fract - force1*mass2fract;
        float4 rand1 = random[randomIndex+2*i];
        float4 rand2 = random[randomIndex+2*i+1];
        cmVel.x = vscale*cmVel.x + fscale*invTotalMass*cmForce.x + noisescale*sqrtInvTotalMass*rand1.x;
        cmVel.y = vscale*cmVel.y + fscale*invTotalMass*cmForce.y + noisescale*sqrtInvTotalMass*rand1.y;
        cmVel.z = vscale*cmVel.z + fscale*invTotalMass*cmForce.z + noisescale*sqrtInvTotalMass*rand1.z;
        relVel.x = vscaleDrude*relVel.x + fscaleDrude*invReducedMass*relForce.x + noisescaleDrude*sqrtInvReducedMass*rand2.x;
        relVel.y = vscaleDrude*relVel.y + fscaleDrude*invReducedMass*relForce.y + noisescaleDrude*sqrtInvReducedMass*rand2.y;
        relVel.z = vscaleDrude*relVel.z + fscaleDrude*invReducedMass*relForce.z + noisescaleDrude*sqrtInvReducedMass*rand2.z;
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

/**
 * Apply hard wall constraints
 */
__kernel void applyHardWallConstraints(__global real4* restrict posq, __global real4* restrict posqCorrection, __global mixed4* restrict velm,
        __global const int2* restrict pairParticles, __global const mixed2* restrict dt, mixed maxDrudeDistance, mixed hardwallscaleDrude) {
    mixed stepSize = dt[0].y;
    for (int i = get_global_id(0); i < NUM_PAIRS; i += get_global_size(0)) {
        int2 particles = pairParticles[i];
#ifdef USE_MIXED_PRECISION
        real4 posReal1 = posq[particles.x];
        real4 posReal2 = posq[particles.y];
        real4 posCorr1 = posqCorrection[particles.x];
        real4 posCorr2 = posqCorrection[particles.y];
        mixed4 pos1 = (mixed4) (posReal1.x+(mixed)posCorr1.x, posReal1.y+(mixed)posCorr1.y, posReal1.z+(mixed)posCorr1.z, posReal1.w);
        mixed4 pos2 = (mixed4) (posReal2.x+(mixed)posCorr2.x, posReal2.y+(mixed)posCorr2.y, posReal2.z+(mixed)posCorr2.z, posReal2.w);
#else
        mixed4 pos1 = posq[particles.x];
        mixed4 pos2 = posq[particles.y];
#endif
        mixed4 delta = pos1-pos2;
        mixed r = sqrt(delta.x*delta.x + delta.y*delta.y + delta.z*delta.z);
        mixed rInv = 1/r;
        if (rInv*maxDrudeDistance < 1) {
            // The constraint has been violated, so make the inter-particle distance "bounce"
            // off the hard wall.

            mixed4 bondDir = delta*rInv;
            mixed4 vel1 = velm[particles.x];
            mixed4 vel2 = velm[particles.y];
            mixed mass1 = 1/vel1.w;
            mixed mass2 = 1/vel2.w;
            mixed deltaR = r-maxDrudeDistance;
            mixed deltaT = stepSize;
            mixed dotvr1 = vel1.x*bondDir.x + vel1.y*bondDir.y + vel1.z*bondDir.z;
            mixed4 vb1 = bondDir*dotvr1;
            mixed4 vp1 = vel1-vb1;
            if (vel2.w == 0) {
                // The parent particle is massless, so move only the Drude particle.

                if (dotvr1 != 0)
                    deltaT = deltaR/fabs(dotvr1);
                if (deltaT > stepSize)
                    deltaT = stepSize;
                dotvr1 = -dotvr1*hardwallscaleDrude/(fabs(dotvr1)*sqrt(mass1));
                mixed dr = -deltaR + deltaT*dotvr1;
                pos1.xyz += bondDir.xyz*dr;
#ifdef USE_MIXED_PRECISION
                posq[particles.x] = (real4) ((real) pos1.x, (real) pos1.y, (real) pos1.z, (real) pos1.w);
                posqCorrection[particles.x] = (real4) (pos1.x-(real) pos1.x, pos1.y-(real) pos1.y, pos1.z-(real) pos1.z, 0);
#else
                posq[particles.x] = pos1;
#endif
                vel1.xyz = vp1.xyz + bondDir.xyz*dotvr1;
                velm[particles.x] = vel1;
            }
            else {
                // Move both particles.

                mixed invTotalMass = 1/(mass1+mass2);
                mixed dotvr2 = vel2.x*bondDir.x + vel2.y*bondDir.y + vel2.z*bondDir.z;
                mixed4 vb2 = bondDir*dotvr2;
                mixed4 vp2 = vel2-vb2;
                mixed vbCMass = (mass1*dotvr1 + mass2*dotvr2)*invTotalMass;
                dotvr1 -= vbCMass;
                dotvr2 -= vbCMass;
                if (dotvr1 != dotvr2)
                    deltaT = deltaR/fabs(dotvr1-dotvr2);
                if (deltaT > stepSize)
                    deltaT = stepSize;
                mixed vBond = hardwallscaleDrude/sqrt(mass1);
                dotvr1 = -dotvr1*vBond*mass2*invTotalMass/fabs(dotvr1);
                dotvr2 = -dotvr2*vBond*mass1*invTotalMass/fabs(dotvr2);
                mixed dr1 = -deltaR*mass2*invTotalMass + deltaT*dotvr1;
                mixed dr2 = deltaR*mass1*invTotalMass + deltaT*dotvr2;
                dotvr1 += vbCMass;
                dotvr2 += vbCMass;
                pos1.xyz += bondDir.xyz*dr1;
                pos2.xyz += bondDir.xyz*dr2;
#ifdef USE_MIXED_PRECISION
                posq[particles.x] = (real4) ((real) pos1.x, (real) pos1.y, (real) pos1.z, (real) pos1.w);
                posq[particles.y] = (real4) ((real) pos2.x, (real) pos2.y, (real) pos2.z, (real) pos2.w);
                posqCorrection[particles.x] = (real4) (pos1.x-(real) pos1.x, pos1.y-(real) pos1.y, pos1.z-(real) pos1.z, 0);
                posqCorrection[particles.y] = (real4) (pos2.x-(real) pos2.x, pos2.y-(real) pos2.y, pos2.z-(real) pos2.z, 0);
#else
                posq[particles.x] = pos1;
                posq[particles.y] = pos2;
#endif
                vel1.xyz = vp1.xyz + bondDir.xyz*dotvr1;
                vel2.xyz = vp2.xyz + bondDir.xyz*dotvr2;
                velm[particles.x] = vel1;
                velm[particles.y] = vel2;
            }
        }
    }
}
