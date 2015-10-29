/**
 * Perform the first step of Langevin integration.
 */

extern "C" __global__ void integrateDrudeLangevinPart1(mixed4* __restrict__ velm, const long long* __restrict__ force, mixed4* __restrict__ posDelta,
        const int* __restrict__ normalParticles, const int2* __restrict__ pairParticles, const mixed2* __restrict__ dt, mixed vscale, mixed fscale,
        mixed noisescale, mixed vscaleDrude, mixed fscaleDrude, mixed noisescaleDrude, const float4* __restrict__ random, unsigned int randomIndex) {
    mixed stepSize = dt[0].y;
    
    // Update normal particles.

    for (int i = blockIdx.x*blockDim.x+threadIdx.x; i < NUM_NORMAL_PARTICLES; i += blockDim.x*gridDim.x) {
        int index = normalParticles[i];
        mixed4 velocity = velm[index];
        if (velocity.w != 0) {
            mixed sqrtInvMass = SQRT(velocity.w);
            float4 rand = random[randomIndex+index];
            velocity.x = vscale*velocity.x + fscale*velocity.w*force[index] + noisescale*sqrtInvMass*rand.x;
            velocity.y = vscale*velocity.y + fscale*velocity.w*force[index+PADDED_NUM_ATOMS] + noisescale*sqrtInvMass*rand.y;
            velocity.z = vscale*velocity.z + fscale*velocity.w*force[index+PADDED_NUM_ATOMS*2] + noisescale*sqrtInvMass*rand.z;
            velm[index] = velocity;
            posDelta[index] = make_mixed4(stepSize*velocity.x, stepSize*velocity.y, stepSize*velocity.z, 0);
        }
    }
    
    // Update Drude particle pairs.
    
    randomIndex += NUM_NORMAL_PARTICLES;
    for (int i = blockIdx.x*blockDim.x+threadIdx.x; i < NUM_PAIRS; i += blockDim.x*gridDim.x) {
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
        mixed3 force1 = make_mixed3(force[particles.x], force[particles.x+PADDED_NUM_ATOMS], force[particles.x+PADDED_NUM_ATOMS*2]);
        mixed3 force2 = make_mixed3(force[particles.y], force[particles.y+PADDED_NUM_ATOMS], force[particles.y+PADDED_NUM_ATOMS*2]);
        mixed3 cmForce = force1+force2;
        mixed3 relForce = force2*mass1fract - force1*mass2fract;
        float4 rand1 = random[randomIndex+2*i];
        float4 rand2 = random[randomIndex+2*i+1];
        cmVel.x = vscale*cmVel.x + fscale*invTotalMass*cmForce.x + noisescale*sqrtInvTotalMass*rand1.x;
        cmVel.y = vscale*cmVel.y + fscale*invTotalMass*cmForce.y + noisescale*sqrtInvTotalMass*rand1.y;
        cmVel.z = vscale*cmVel.z + fscale*invTotalMass*cmForce.z + noisescale*sqrtInvTotalMass*rand1.z;
        relVel.x = vscaleDrude*relVel.x + fscaleDrude*invReducedMass*relForce.x + noisescaleDrude*sqrtInvReducedMass*rand2.x;
        relVel.y = vscaleDrude*relVel.y + fscaleDrude*invReducedMass*relForce.y + noisescaleDrude*sqrtInvReducedMass*rand2.y;
        relVel.z = vscaleDrude*relVel.z + fscaleDrude*invReducedMass*relForce.z + noisescaleDrude*sqrtInvReducedMass*rand2.z;
        velocity1.x = cmVel.x-relVel.x*mass2fract;
        velocity1.y = cmVel.y-relVel.y*mass2fract;
        velocity1.z = cmVel.z-relVel.z*mass2fract;
        velocity2.x = cmVel.x+relVel.x*mass1fract;
        velocity2.y = cmVel.y+relVel.y*mass1fract;
        velocity2.z = cmVel.z+relVel.z*mass1fract;
        velm[particles.x] = velocity1;
        velm[particles.y] = velocity2;
        posDelta[particles.x] = make_mixed4(stepSize*velocity1.x, stepSize*velocity1.y, stepSize*velocity1.z, 0);
        posDelta[particles.y] = make_mixed4(stepSize*velocity2.x, stepSize*velocity2.y, stepSize*velocity2.z, 0);
    }
}

/**
 * Perform the second step of Langevin integration.
 */

extern "C" __global__ void integrateDrudeLangevinPart2(real4* __restrict__ posq, real4* __restrict__ posqCorrection, const mixed4* __restrict__ posDelta, mixed4* __restrict__ velm, const mixed2* __restrict__ dt) {
    double invStepSize = 1.0/dt[0].y;
    int index = blockIdx.x*blockDim.x+threadIdx.x;
    while (index < NUM_ATOMS) {
        mixed4 vel = velm[index];
        if (vel.w != 0) {
#ifdef USE_MIXED_PRECISION
 
            real4 pos1 = posq[index];
            real4 pos2 = posqCorrection[index];
            mixed4 pos = make_mixed4(pos1.x+(mixed)pos2.x, pos1.y+(mixed)pos2.y, pos1.z+(mixed)pos2.z, pos1.w);
#else
            real4 pos = posq[index];
#endif
            mixed4 delta = posDelta[index];
            pos.x += delta.x;
            pos.y += delta.y;
            pos.z += delta.z;
            vel.x = (mixed) (invStepSize*delta.x);
            vel.y = (mixed) (invStepSize*delta.y);
            vel.z = (mixed) (invStepSize*delta.z);
#ifdef USE_MIXED_PRECISION
            posq[index] = make_real4((real) pos.x, (real) pos.y, (real) pos.z, (real) pos.w);
            posqCorrection[index] = make_real4(pos.x-(real) pos.x, pos.y-(real) pos.y, pos.z-(real) pos.z, 0);
#else
            posq[index] = pos;
#endif
            velm[index] = vel;
        }
        index += blockDim.x*gridDim.x;
    }
}

/**
 * Apply hard wall constraints
 */
extern "C" __global__ void applyHardWallConstraints(real4* __restrict__ posq, real4* __restrict__ posqCorrection, mixed4* __restrict__ velm,
        const int2* __restrict__ pairParticles, const mixed2* __restrict__ dt, mixed maxDrudeDistance, mixed hardwallscaleDrude) {
    mixed stepSize = dt[0].y;
    for (int i = blockIdx.x*blockDim.x+threadIdx.x; i < NUM_PAIRS; i += blockDim.x*gridDim.x) {
        int2 particles = pairParticles[i];
#ifdef USE_MIXED_PRECISION
        real4 posReal1 = posq[particles.x];
        real4 posReal2 = posq[particles.y];
        real4 posCorr1 = posqCorrection[particles.x];
        real4 posCorr2 = posqCorrection[particles.y];
        mixed4 pos1 = make_mixed4(posReal1.x+(mixed)posCorr1.x, posReal1.y+(mixed)posCorr1.y, posReal1.z+(mixed)posCorr1.z, posReal1.w);
        mixed4 pos2 = make_mixed4(posReal2.x+(mixed)posCorr2.x, posReal2.y+(mixed)posCorr2.y, posReal2.z+(mixed)posCorr2.z, posReal2.w);
#else
        mixed4 pos1 = posq[particles.x];
        mixed4 pos2 = posq[particles.y];
#endif
        mixed4 delta = pos1-pos2;
        mixed r = SQRT(delta.x*delta.x + delta.y*delta.y + delta.z*delta.z);
        mixed rInv = RECIP(r);
        if (rInv*maxDrudeDistance < 1) {
            // The constraint has been violated, so make the inter-particle distance "bounce"
            // off the hard wall.

            mixed4 bondDir = delta*rInv;
            mixed4 vel1 = velm[particles.x];
            mixed4 vel2 = velm[particles.y];
            mixed mass1 = RECIP(vel1.w);
            mixed mass2 = RECIP(vel2.w);
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
                dotvr1 = -dotvr1*hardwallscaleDrude/(fabs(dotvr1)*SQRT(mass1));
                mixed dr = -deltaR + deltaT*dotvr1;
                pos1.x += bondDir.x*dr;
                pos1.y += bondDir.y*dr;
                pos1.z += bondDir.z*dr;
#ifdef USE_MIXED_PRECISION
                posq[particles.x] = make_real4((real) pos1.x, (real) pos1.y, (real) pos1.z, (real) pos1.w);
                posqCorrection[particles.x] = make_real4(pos1.x-(real) pos1.x, pos1.y-(real) pos1.y, pos1.z-(real) pos1.z, 0);
#else
                posq[particles.x] = pos1;
#endif
                vel1.x = vp1.x + bondDir.x*dotvr1;
                vel1.y = vp1.y + bondDir.y*dotvr1;
                vel1.z = vp1.z + bondDir.z*dotvr1;
                velm[particles.x] = vel1;
            }
            else {
                // Move both particles.

                mixed invTotalMass = RECIP(mass1+mass2);
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
                mixed vBond = hardwallscaleDrude/SQRT(mass1);
                dotvr1 = -dotvr1*vBond*mass2*invTotalMass/fabs(dotvr1);
                dotvr2 = -dotvr2*vBond*mass1*invTotalMass/fabs(dotvr2);
                mixed dr1 = -deltaR*mass2*invTotalMass + deltaT*dotvr1;
                mixed dr2 = deltaR*mass1*invTotalMass + deltaT*dotvr2;
                dotvr1 += vbCMass;
                dotvr2 += vbCMass;
                pos1.x += bondDir.x*dr1;
                pos1.y += bondDir.y*dr1;
                pos1.z += bondDir.z*dr1;
                pos2.x += bondDir.x*dr2;
                pos2.y += bondDir.y*dr2;
                pos2.z += bondDir.z*dr2;
#ifdef USE_MIXED_PRECISION
                posq[particles.x] = make_real4((real) pos1.x, (real) pos1.y, (real) pos1.z, (real) pos1.w);
                posq[particles.y] = make_real4((real) pos2.x, (real) pos2.y, (real) pos2.z, (real) pos2.w);
                posqCorrection[particles.x] = make_real4(pos1.x-(real) pos1.x, pos1.y-(real) pos1.y, pos1.z-(real) pos1.z, 0);
                posqCorrection[particles.y] = make_real4(pos2.x-(real) pos2.x, pos2.y-(real) pos2.y, pos2.z-(real) pos2.z, 0);
#else
                posq[particles.x] = pos1;
                posq[particles.y] = pos2;
#endif
                vel1.x = vp1.x + bondDir.x*dotvr1;
                vel1.y = vp1.y + bondDir.y*dotvr1;
                vel1.z = vp1.z + bondDir.z*dotvr1;
                vel2.x = vp2.x + bondDir.x*dotvr2;
                vel2.y = vp2.y + bondDir.y*dotvr2;
                vel2.z = vp2.z + bondDir.z*dotvr2;
                velm[particles.x] = vel1;
                velm[particles.y] = vel2;
            }
        }
    }
}
