/**
 * Perform the first part of integration: velocity step.
 */
KERNEL void integrateQTBPart1(int numAtoms, int paddedNumAtoms, mixed dt, int stepIndex, GLOBAL mixed4* RESTRICT velm,
        GLOBAL const mm_long* RESTRICT force, GLOBAL mixed* RESTRICT segmentVelocity) {
    mixed fscale = dt/(mixed) 0x100000000;
    for (int atom = GLOBAL_ID; atom < numAtoms; atom += GLOBAL_SIZE) {
        mixed4 velocity = velm[atom];
        segmentVelocity[3*numAtoms*stepIndex + atom] = velocity.x;
        segmentVelocity[3*numAtoms*stepIndex + numAtoms + atom] = velocity.y;
        segmentVelocity[3*numAtoms*stepIndex + 2*numAtoms + atom] = velocity.z;
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
        GLOBAL mixed4* RESTRICT posDelta, GLOBAL mixed4* RESTRICT oldDelta, GLOBAL const mixed* RESTRICT randomForce) {
    mixed halfdt = 0.5f*dt;
    mixed vscale = EXP(-dt*friction);
    for (int atom = GLOBAL_ID; atom < numAtoms; atom += GLOBAL_SIZE) {
        mixed4 velocity = velm[atom];
        if (velocity.w != 0.0) {
            mixed4 delta = make_mixed4(halfdt*velocity.x, halfdt*velocity.y, halfdt*velocity.z, 0);
            mixed fscale = dt*velocity.w;
            velocity.x = vscale*velocity.x + fscale*randomForce[3*numAtoms*stepIndex + atom];
            velocity.y = vscale*velocity.y + fscale*randomForce[3*numAtoms*stepIndex + numAtoms + atom];
            velocity.z = vscale*velocity.z + fscale*randomForce[3*numAtoms*stepIndex + 2*numAtoms + atom];
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
         GLOBAL mixed4* RESTRICT posDelta, GLOBAL mixed4* RESTRICT oldDelta
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
 * Update the buffer of white noise for the next segment.
 */
KERNEL void generateNoise(int numAtoms, int segmentLength, GLOBAL float* RESTRICT noise, GLOBAL const float4* RESTRICT random, unsigned int randomIndex) {
    int offset = 3*numAtoms*segmentLength;
    for (int i = GLOBAL_ID; i < 6*numAtoms*segmentLength; i += GLOBAL_SIZE)
        noise[i] = noise[i+offset];
    for (int i = GLOBAL_ID; i < numAtoms*segmentLength; i += GLOBAL_SIZE) {
        float4 r = random[randomIndex+i];
        noise[2*offset+i] = r.x;
        noise[2*offset+numAtoms+i] = r.y;
        noise[2*offset+2*numAtoms+i] = r.z;
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
        GLOBAL mixed* RESTRICT cutoffFunction, GLOBAL mixed2* RESTRICT workspace) {
    const int fftLength = 3*segmentLength;
    const int numFreq = (fftLength+1)/2;
    mixed2* data0 = &workspace[GROUP_ID*3*fftLength];
    mixed2* data1 = &data0[fftLength];
    mixed2* w = &data1[fftLength];
    for (int i = LOCAL_ID; i < fftLength; i += LOCAL_SIZE)
        w[i] = make_mixed2(cos(-i*2*M_PI/fftLength), sin(-i*2*M_PI/fftLength));
    for (int i = GROUP_ID; i < 3*numAtoms; i += NUM_GROUPS) {
        int atom = i/3;
        int axis = i%3;
        mixed invMass = velm[atom].w;
        if (invMass != 0) {
            for (int j = LOCAL_ID; j < fftLength; j += LOCAL_SIZE)
                data0[j] = make_mixed2(noise[numAtoms*(3*j+axis)+atom], 0);
            SYNC_THREADS
            FFT_FORWARD
            for (int j = LOCAL_ID; j < numFreq; j += LOCAL_SIZE) {
                mixed f = M_PI*j/(numFreq*dt);
                mixed gamma = friction;//adaptedFriction[type][i];
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
