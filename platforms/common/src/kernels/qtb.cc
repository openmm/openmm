/**
 * Perform the first part of integration: velocity step.
 */
KERNEL void integrateQTBPart1(int numAtoms, int paddedNumAtoms, mixed dt, GLOBAL mixed4* RESTRICT velm,
        GLOBAL const mm_long* RESTRICT force, GLOBAL const int* RESTRICT atomOrder) {
    mixed fscale = dt/(mixed) 0x100000000;
    for (int atom = GLOBAL_ID; atom < numAtoms; atom += GLOBAL_SIZE) {
        int atomIndex = atomOrder[atom];
        mixed4 velocity = velm[atom];
        if (velocity.w != 0.0) {
            velocity.x += fscale*velocity.w*force[atom];
            velocity.y += fscale*velocity.w*force[atom+paddedNumAtoms];
            velocity.z += fscale*velocity.w*force[atom+paddedNumAtoms*2];
            velm[atom] = velocity;
        }
    }
}

/**
 * Perform the second part of integration: position half step, then interact with heat bath,
 * then another position half step.
 */
KERNEL void integrateQTBPart2(int numAtoms, mixed dt, mixed friction, int stepIndex, GLOBAL mixed4* RESTRICT velm,
        GLOBAL mixed4* RESTRICT posDelta, GLOBAL mixed4* RESTRICT oldDelta, GLOBAL const mixed* RESTRICT randomForce,
        GLOBAL mixed* RESTRICT segmentVelocity, GLOBAL const int* RESTRICT atomOrder) {
    mixed halfdt = 0.5f*dt;
    mixed vscale = EXP(-dt*friction);
    mixed halfvscale = EXP(-halfdt*friction);
    for (int atom = GLOBAL_ID; atom < numAtoms; atom += GLOBAL_SIZE) {
        int atomIndex = atomOrder[atom];
        mixed4 velocity = velm[atom];
        if (velocity.w != 0.0) {
            mixed4 delta = make_mixed4(halfdt*velocity.x, halfdt*velocity.y, halfdt*velocity.z, 0);
            mixed fscale = dt*velocity.w;
            mixed3 f = make_mixed3(randomForce[3*numAtoms*stepIndex + atomIndex],
                                   randomForce[3*numAtoms*stepIndex + numAtoms + atomIndex],
                                   randomForce[3*numAtoms*stepIndex + 2*numAtoms + atomIndex]);
            segmentVelocity[3*numAtoms*stepIndex + atomIndex] = halfvscale*velocity.x + 0.5f*fscale*f.x;
            segmentVelocity[3*numAtoms*stepIndex + numAtoms + atomIndex] = halfvscale*velocity.y + 0.5f*fscale*f.y;
            segmentVelocity[3*numAtoms*stepIndex + 2*numAtoms + atomIndex] = halfvscale*velocity.z + 0.5f*fscale*f.z;
            velocity.x = vscale*velocity.x + fscale*f.x;
            velocity.y = vscale*velocity.y + fscale*f.y;
            velocity.z = vscale*velocity.z + fscale*f.z;
            velm[atom] = velocity;
            delta += make_mixed4(halfdt*velocity.x, halfdt*velocity.y, halfdt*velocity.z, 0);
            posDelta[atom] = delta;
            oldDelta[atom] = delta;
        }
    }
}

/**
 * Perform the third part of integration: apply constraint forces to velocities, then record
 * the constrained positions.
 */
KERNEL void integrateQTBPart3(int numAtoms, mixed dt, GLOBAL real4* RESTRICT posq, GLOBAL mixed4* RESTRICT velm,
         GLOBAL const mixed4* RESTRICT posDelta, GLOBAL const mixed4* RESTRICT oldDelta
#ifdef USE_MIXED_PRECISION
        , GLOBAL real4* RESTRICT posqCorrection
#endif
        ) {
    mixed invDt = 1/dt;
    for (int index = GLOBAL_ID; index < numAtoms; index += GLOBAL_SIZE) {
        mixed4 velocity = velm[index];
        if (velocity.w != 0.0) {
            mixed4 delta = posDelta[index];
            velocity.x += (delta.x-oldDelta[index].x)*invDt;
            velocity.y += (delta.y-oldDelta[index].y)*invDt;
            velocity.z += (delta.z-oldDelta[index].z)*invDt;
            velm[index] = velocity;
#ifdef USE_MIXED_PRECISION
            real4 pos1 = posq[index];
            real4 pos2 = posqCorrection[index];
            mixed4 pos = make_mixed4(pos1.x+(mixed)pos2.x, pos1.y+(mixed)pos2.y, pos1.z+(mixed)pos2.z, pos1.w);
#else
            real4 pos = posq[index];
#endif
            pos.x += delta.x;
            pos.y += delta.y;
            pos.z += delta.z;
#ifdef USE_MIXED_PRECISION
            posq[index] = make_real4((real) pos.x, (real) pos.y, (real) pos.z, (real) pos.w);
            posqCorrection[index] = make_real4(pos.x-(real) pos.x, pos.y-(real) pos.y, pos.z-(real) pos.z, 0);
#else
            posq[index] = pos;
#endif
        }
    }
}

/**
 * Update the buffer of white noise.
 */
KERNEL void generateNoise(int numAtoms, int segmentLength, GLOBAL float* RESTRICT noise, GLOBAL const float* RESTRICT random, unsigned int randomIndex) {
    randomIndex = 4*randomIndex; // Interpret it as float instead of float4
    int fftLength = 3*segmentLength;
    for (int i = GROUP_ID; i < 3*numAtoms; i += NUM_GROUPS) {
        // Copy segment 2 over to segment 1

        for (int j = LOCAL_ID; j < segmentLength; j += LOCAL_SIZE)
            noise[i*fftLength+j] = noise[i*fftLength+j+segmentLength];
        SYNC_THREADS

        // Copy segment 3 over to segment 2

        for (int j = LOCAL_ID; j < segmentLength; j += LOCAL_SIZE)
            noise[i*fftLength+j+segmentLength] = noise[i*fftLength+j+2*segmentLength];
        SYNC_THREADS
 
        // Fill segment 3 with new random values.

        for (int j = LOCAL_ID; j < segmentLength; j += LOCAL_SIZE)
            noise[i*fftLength+j+2*segmentLength] = random[randomIndex+i*segmentLength+j];
        SYNC_THREADS
    }
}

DEVICE mixed2 multiplyComplex(mixed2 c1, mixed2 c2) {
    return make_mixed2(c1.x*c2.x-c1.y*c2.y, c1.x*c2.y+c1.y*c2.x);
}

DEVICE mixed2 multiplyComplexConj(mixed2 c1, mixed2 c2) {
    return make_mixed2(c1.x*c2.x+c1.y*c2.y, c1.x*c2.y-c1.y*c2.x);
}

/**
 * Generate the random force for the next segment.
 */
KERNEL void generateRandomForce(int numAtoms, int segmentLength, mixed dt, mixed friction, GLOBAL float* RESTRICT noise,
        GLOBAL mixed* RESTRICT randomForce, GLOBAL mixed4* RESTRICT velm, GLOBAL mixed* RESTRICT thetad,
        GLOBAL mixed* RESTRICT cutoffFunction, GLOBAL int* RESTRICT particleType, GLOBAL mixed* RESTRICT adaptedFriction,
        GLOBAL mixed2* RESTRICT workspace) {
    const int fftLength = 3*segmentLength;
    const int numFreq = (fftLength+1)/2;
    GLOBAL mixed2* data0 = &workspace[GROUP_ID*3*fftLength];
    GLOBAL mixed2* data1 = &data0[fftLength];
    GLOBAL mixed2* w = &data1[fftLength];
    for (int i = LOCAL_ID; i < fftLength; i += LOCAL_SIZE)
        w[i] = make_mixed2(cos(-i*2*M_PI/fftLength), sin(-i*2*M_PI/fftLength));
    for (int i = GROUP_ID; i < 3*numAtoms; i += NUM_GROUPS) {
        int atom = i/3;
        int axis = i%3;
        mixed invMass = velm[atom].w;
        int type = particleType[atom];
        if (invMass != 0) {
            for (int j = LOCAL_ID; j < fftLength; j += LOCAL_SIZE)
                data0[j] = make_mixed2(noise[i*fftLength+j], 0);
            SYNC_THREADS
            FFT_FORWARD
            for (int j = LOCAL_ID; j < numFreq; j += LOCAL_SIZE) {
                mixed f = M_PI*j/(numFreq*dt);
                mixed gamma = adaptedFriction[type*numFreq+j];
                mixed cw = (1 - 2*EXP(-dt*friction)*COS(f*dt) + EXP(-2*friction*dt)) / ((friction*friction+f*f)*dt*dt);
                mixed2 d = RECIP_DATA[j] * SQRT(cutoffFunction[j]*thetad[j]*cw*gamma/friction);
                RECIP_DATA[j] = d;
                if (j != 0 && j*2 != fftLength)
                    RECIP_DATA[fftLength-j] = make_mixed2(d.x, -d.y);
            }
            SYNC_THREADS
            FFT_BACKWARD
            const mixed scale = SQRT(2*friction/(dt*invMass))/fftLength;
            for (int j = LOCAL_ID; j < segmentLength; j += LOCAL_SIZE)
                randomForce[numAtoms*(3*j+axis)+atom] = scale*data0[segmentLength+j].x;
            SYNC_THREADS
        }
    }
}

/**
 * Update the friction rates used for generating noise, part 1: compute the error in the
 * fluctuation dissipation theorem.
 */
KERNEL void adaptFrictionPart1(int numAtoms, int segmentLength, GLOBAL const mixed4* RESTRICT velm, GLOBAL const int* RESTRICT particleType,
        GLOBAL const mixed* RESTRICT randomForce, GLOBAL const mixed* RESTRICT segmentVelocity, GLOBAL const mixed* RESTRICT adaptedFriction,
        GLOBAL mm_ulong* RESTRICT dfdt, GLOBAL mixed2* RESTRICT workspace) {
    const int fftLength = 3*segmentLength;
    const int numFreq = (fftLength+1)/2;
    GLOBAL mixed2* data0 = &workspace[GROUP_ID*3*fftLength];
    GLOBAL mixed2* data1 = &data0[fftLength];
    GLOBAL mixed2* w = &data1[fftLength];
    for (int i = LOCAL_ID; i < fftLength; i += LOCAL_SIZE)
        w[i] = make_mixed2(cos(-i*2*M_PI/fftLength), sin(-i*2*M_PI/fftLength));
    for (int i = GROUP_ID; i < 3*numAtoms; i += NUM_GROUPS) {
        int atom = i/3;
        int axis = i%3;
        int type = particleType[atom];
        mixed invMass = velm[atom].w;
        if (invMass != 0) {
            // Pack the velocities and forces together so we can transform both at once.
            for (int j = LOCAL_ID; j < segmentLength; j += LOCAL_SIZE)
                data0[j] = make_mixed2(segmentVelocity[numAtoms*(3*j+axis)+atom], randomForce[numAtoms*(3*j+axis)+atom]);
            for (int j = segmentLength+LOCAL_ID; j < fftLength; j += LOCAL_SIZE)
                data0[j] = make_mixed2(0);
            SYNC_THREADS
            ADAPTATION_FFT
            for (int j = LOCAL_ID; j < numFreq; j += LOCAL_SIZE) {
                mixed2 d1 = ADAPTATION_RECIP[j];
                mixed2 d2 = (j == 0 ? ADAPTATION_RECIP[0] : ADAPTATION_RECIP[fftLength-j]);
                mixed2 v = 0.5f*make_mixed2(d1.x+d2.x, d1.y-d2.y);
                mixed2 f = 0.5f*make_mixed2(d1.y+d2.y, -d1.x+d2.x);
                mixed cvv = v.x*v.x + v.y*v.y;
                mixed2 cvf = multiplyComplexConj(f, v);
                mixed dfdtinc = adaptedFriction[type*numFreq+j]*cvv/invMass - cvf.x;
                ATOMIC_ADD(&dfdt[type*numFreq+j], (mm_ulong) realToFixedPoint(dfdtinc));
            }
            SYNC_THREADS
        }
    }
}

/**
 * Update the friction rates used for generating noise, part 2: update the friction based
 * on the error.
 */
KERNEL void adaptFrictionPart2(int numTypes, int segmentLength, mixed dt, GLOBAL const int* RESTRICT typeParticleCount,
        GLOBAL const float* RESTRICT typeAdaptationRate, GLOBAL mixed* RESTRICT adaptedFriction, GLOBAL const mm_long* RESTRICT dfdt) {
    int numFreq = (3*segmentLength+1)/2;
    for (int type = GROUP_ID; type < numTypes; type += NUM_GROUPS) {
        mixed scale = dt*typeAdaptationRate[type]/(3*typeParticleCount[type]*segmentLength)/(mixed) 0x100000000;
        for (int i = LOCAL_ID; i < numFreq; i += LOCAL_SIZE) {
            mixed delta = -scale*dfdt[type*numFreq+i];
            adaptedFriction[type*numFreq+i] = max((mixed) 0, adaptedFriction[type*numFreq+i]+delta);
        }
    }
}