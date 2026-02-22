DEVICE mixed3 multiplyComplexRealPart(mixed2 c1, mixed3 c2r, mixed3 c2i) {
    return c1.x*c2r-c1.y*c2i;
}

DEVICE mixed3 multiplyComplexImagPart(mixed2 c1, mixed3 c2r, mixed3 c2i) {
    return c1.x*c2i+c1.y*c2r;
}

DEVICE mixed3 multiplyComplexRealPartConj(mixed2 c1, mixed3 c2r, mixed3 c2i) {
    return c1.x*c2r+c1.y*c2i;
}

DEVICE mixed3 multiplyComplexImagPartConj(mixed2 c1, mixed3 c2r, mixed3 c2i) {
    return c1.x*c2i-c1.y*c2r;
}

/**
 * Apply the PILE thermostat.
 * 
 * @param applyToCentroid If true (PILE mode), apply Langevin thermostat to centroid mode.
 *                        If false (PILE_G mode), skip centroid (it's handled by Bussi separately).
 */
KERNEL void applyPileThermostat(GLOBAL mixed4* velm, GLOBAL float4* random, unsigned int randomIndex,
        mixed dt, mixed kT, mixed friction, int applyToCentroid) {
    const int numBlocks = GLOBAL_SIZE/NUM_COPIES;
    const int blockStart = NUM_COPIES*(LOCAL_ID/NUM_COPIES);
    const int indexInBlock = LOCAL_ID-blockStart;
    const mixed nkT = NUM_COPIES*kT;
    const mixed twown = 2.0f*nkT/HBAR;
    const mixed c1_0 = exp(-0.5f*dt*friction);
    const mixed c2_0 = sqrt(1.0f-c1_0*c1_0);
    LOCAL mixed3 v[2*THREAD_BLOCK_SIZE];
    LOCAL mixed3 temp[2*THREAD_BLOCK_SIZE];
    LOCAL mixed2 w[NUM_COPIES];
    LOCAL_ARG mixed3* vreal = &v[blockStart];
    LOCAL_ARG mixed3* vimag = &v[blockStart+LOCAL_SIZE];
    if (LOCAL_ID < NUM_COPIES)
        w[indexInBlock] = make_mixed2(cos(-indexInBlock*2*M_PI/NUM_COPIES), sin(-indexInBlock*2*M_PI/NUM_COPIES));
    SYNC_THREADS;
    randomIndex += NUM_COPIES*((GLOBAL_ID)/NUM_COPIES);
    for (int particle = (GLOBAL_ID)/NUM_COPIES; particle < NUM_ATOMS; particle += numBlocks) {
        mixed4 particleVelm = velm[particle+indexInBlock*PADDED_NUM_ATOMS];
        mixed invMass = particleVelm.w;
        mixed c3_0 = c2_0*sqrt(nkT*invMass);
        
        // Forward FFT.
        
        vreal[indexInBlock] = SCALE*make_mixed3(particleVelm.x, particleVelm.y, particleVelm.z);
        vimag[indexInBlock] = make_mixed3(0);
        SYNC_THREADS;
        FFT_V_FORWARD

        // Apply the thermostat.

        if (indexInBlock == 0) {
            // Apply a local Langevin thermostat to the centroid mode (PILE mode only).
            // For PILE_G mode, centroid is handled by Bussi thermostat separately.
            
            if (applyToCentroid) {
                float4 rand = random[randomIndex];
                vreal[0] = vreal[0]*c1_0 + c3_0*make_mixed3(rand.x, rand.y, rand.z);
            }
            // If not applying to centroid, vreal[0] and vimag[0] are unchanged
        }
        else {
            // Use critical damping white noise for the remaining (internal) modes.
            // This is applied for both PILE and PILE_G modes.

            int k = (indexInBlock <= NUM_COPIES/2 ? indexInBlock : NUM_COPIES-indexInBlock);
            const bool isCenter = (NUM_COPIES%2 == 0 && k == NUM_COPIES/2);
            const mixed wk = twown*sin(k*M_PI/NUM_COPIES);
            const mixed c1 = exp(-wk*dt);
            const mixed c2 = sqrt((1.0f-c1*c1)/2.0f) * (isCenter ? sqrt(2.0f) : 1.0f);
            const mixed c3 = c2*sqrt(nkT*invMass);
            float4 rand1 = random[randomIndex+k];
            float4 rand2 = (isCenter ? make_float4(0) : random[randomIndex+NUM_COPIES-k]);
            vreal[indexInBlock] = c1*vreal[indexInBlock] + c3*make_mixed3(rand1.x, rand1.y, rand1.z);
            vimag[indexInBlock] = c1*vimag[indexInBlock] + c3*(indexInBlock < NUM_COPIES/2 ? make_mixed3(rand2.x, rand2.y, rand2.z) : make_mixed3(-rand2.x, -rand2.y, -rand2.z));
        }
        SYNC_THREADS;
        
        // Inverse FFT.
        
        FFT_V_BACKWARD
        if (invMass != 0)
            velm[particle+indexInBlock*PADDED_NUM_ATOMS] = make_mixed4(SCALE*vreal[indexInBlock].x, SCALE*vreal[indexInBlock].y, SCALE*vreal[indexInBlock].z, particleVelm.w);
        randomIndex += GLOBAL_SIZE;
    }
}

/**
 * Advance the positions and velocities.
 */
KERNEL void integrateStep(GLOBAL mixed4* posq, GLOBAL mixed4* velm, GLOBAL mm_long* force, mixed dt, mixed kT) {
    const int numBlocks = (GLOBAL_SIZE)/NUM_COPIES;
    const int blockStart = NUM_COPIES*(LOCAL_ID/NUM_COPIES);
    const int indexInBlock = LOCAL_ID-blockStart;
    const mixed nkT = NUM_COPIES*kT;
    const mixed twown = 2.0f*nkT/HBAR;
    const mixed forceScale = 1/(mixed) 0x100000000;
    LOCAL mixed3 q[2*THREAD_BLOCK_SIZE];
    LOCAL mixed3 v[2*THREAD_BLOCK_SIZE];
    LOCAL mixed3 temp[2*THREAD_BLOCK_SIZE];
    LOCAL mixed2 w[NUM_COPIES];

    // Update velocities.
    
    for (int particle = (GLOBAL_ID)/NUM_COPIES; particle < NUM_ATOMS; particle += numBlocks) {
        int index = particle+indexInBlock*PADDED_NUM_ATOMS;
        int forceIndex = particle+indexInBlock*PADDED_NUM_ATOMS*3;
        mixed4 particleVelm = velm[index];
        particleVelm.x += forceScale*force[forceIndex]*(0.5f*dt*particleVelm.w);
        particleVelm.y += forceScale*force[forceIndex+PADDED_NUM_ATOMS]*(0.5f*dt*particleVelm.w);
        particleVelm.z += forceScale*force[forceIndex+PADDED_NUM_ATOMS*2]*(0.5f*dt*particleVelm.w);
        if (particleVelm.w != 0)
            velm[index] = particleVelm;
    }
    
    // Evolve the free ring polymer by transforming to the frequency domain.

    LOCAL_ARG mixed3* qreal = &q[blockStart];
    LOCAL_ARG mixed3* qimag = &q[blockStart+LOCAL_SIZE];
    LOCAL_ARG mixed3* vreal = &v[blockStart];
    LOCAL_ARG mixed3* vimag = &v[blockStart+LOCAL_SIZE];
    if (LOCAL_ID < NUM_COPIES)
        w[indexInBlock] = make_mixed2(cos(-indexInBlock*2*M_PI/NUM_COPIES), sin(-indexInBlock*2*M_PI/NUM_COPIES));
    SYNC_THREADS;
    for (int particle = (GLOBAL_ID)/NUM_COPIES; particle < NUM_ATOMS; particle += numBlocks) {
        mixed4 particlePosq = posq[particle+indexInBlock*PADDED_NUM_ATOMS];
        mixed4 particleVelm = velm[particle+indexInBlock*PADDED_NUM_ATOMS];
        
        // Forward FFT.
        
        qreal[indexInBlock] = SCALE*make_mixed3(particlePosq.x, particlePosq.y, particlePosq.z);
        qimag[indexInBlock] = make_mixed3(0);
        vreal[indexInBlock] = SCALE*make_mixed3(particleVelm.x, particleVelm.y, particleVelm.z);
        vimag[indexInBlock] = make_mixed3(0);
        SYNC_THREADS;
        FFT_Q_FORWARD
        FFT_V_FORWARD

        // Apply the thermostat.

        if (indexInBlock == 0) {
            qreal[0] += vreal[0]*dt;
            qimag[0] += vimag[0]*dt;
        }
        else {
            const mixed wk = twown*sin(indexInBlock*M_PI/NUM_COPIES);
            const mixed wt = wk*dt;
            const mixed coswt = cos(wt);
            const mixed sinwt = sin(wt);
            const mixed3 vprimereal = vreal[indexInBlock]*coswt - qreal[indexInBlock]*(wk*sinwt); // Advance velocity from t to t+dt
            const mixed3 vprimeimag = vimag[indexInBlock]*coswt - qimag[indexInBlock]*(wk*sinwt);
            qreal[indexInBlock] = vreal[indexInBlock]*(sinwt/wk) + qreal[indexInBlock]*coswt; // Advance position from t to t+dt
            qimag[indexInBlock] = vimag[indexInBlock]*(sinwt/wk) + qimag[indexInBlock]*coswt;
            vreal[indexInBlock] = vprimereal;
            vimag[indexInBlock] = vprimeimag;
        }
        SYNC_THREADS;
        
        // Inverse FFT.
        
        FFT_Q_BACKWARD
        FFT_V_BACKWARD
        if (particleVelm.w != 0) {
            posq[particle+indexInBlock*PADDED_NUM_ATOMS] = make_mixed4(SCALE*qreal[indexInBlock].x, SCALE*qreal[indexInBlock].y, SCALE*qreal[indexInBlock].z, particlePosq.w);
            velm[particle+indexInBlock*PADDED_NUM_ATOMS] = make_mixed4(SCALE*vreal[indexInBlock].x, SCALE*vreal[indexInBlock].y, SCALE*vreal[indexInBlock].z, particleVelm.w);
        }
    }
}

/**
 * Advance the velocities by a half step.
 */
KERNEL void advanceVelocities(GLOBAL mixed4* velm, GLOBAL mm_long* force, mixed dt) {
    const int numBlocks = (GLOBAL_SIZE)/NUM_COPIES;
    const int blockStart = NUM_COPIES*(LOCAL_ID/NUM_COPIES);
    const int indexInBlock = LOCAL_ID-blockStart;
    const mixed forceScale = 1/(mixed) 0x100000000;

    // Update velocities.
    
    for (int particle = (GLOBAL_ID)/NUM_COPIES; particle < NUM_ATOMS; particle += numBlocks) {
        int index = particle+indexInBlock*PADDED_NUM_ATOMS;
        int forceIndex = particle+indexInBlock*PADDED_NUM_ATOMS*3;
        mixed4 particleVelm = velm[index];
        particleVelm.x += forceScale*force[forceIndex]*(0.5f*dt*particleVelm.w);
        particleVelm.y += forceScale*force[forceIndex+PADDED_NUM_ATOMS]*(0.5f*dt*particleVelm.w);
        particleVelm.z += forceScale*force[forceIndex+PADDED_NUM_ATOMS*2]*(0.5f*dt*particleVelm.w);
        if (particleVelm.w != 0)
            velm[index] = particleVelm;
    }
}

/**
 * Compute centroid kinetic energy for Bussi thermostat.
 * Each work item computes the contribution from one particle across all beads.
 * Results are stored in the buffer for final reduction.
 */
KERNEL void computeCentroidKE(GLOBAL mixed4* velm, GLOBAL mixed* kineticEnergy) {
    LOCAL mixed localKE[THREAD_BLOCK_SIZE];
    mixed myKE = 0.0f;
    const mixed invNumCopies = 1.0f / NUM_COPIES;
    
    // Each thread processes one particle across all beads
    for (int particle = GLOBAL_ID; particle < NUM_ATOMS; particle += GLOBAL_SIZE) {
        // Compute centroid velocity
        mixed3 centroidVel = make_mixed3(0.0f, 0.0f, 0.0f);
        mixed mass = 0.0f;
        
        for (int copy = 0; copy < NUM_COPIES; copy++) {
            mixed4 particleVelm = velm[particle + copy*PADDED_NUM_ATOMS];
            centroidVel.x += particleVelm.x;
            centroidVel.y += particleVelm.y;
            centroidVel.z += particleVelm.z;
            mass = particleVelm.w; // Same for all copies
        }
        
        centroidVel.x *= invNumCopies;
        centroidVel.y *= invNumCopies;
        centroidVel.z *= invNumCopies;
        
        // Compute kinetic energy contribution (0.5 * m * v^2)
        if (mass != 0.0f) {
            mixed v2 = centroidVel.x*centroidVel.x + centroidVel.y*centroidVel.y + centroidVel.z*centroidVel.z;
            myKE += 0.5f * v2 / mass;  // mass stored as inverse mass in .w
        }
    }
    
    // Parallel reduction within work group
    localKE[LOCAL_ID] = myKE;
    SYNC_THREADS;
    
    for (int offset = 1; offset < THREAD_BLOCK_SIZE; offset *= 2) {
        if (LOCAL_ID + offset < THREAD_BLOCK_SIZE && LOCAL_ID % (2*offset) == 0)
            localKE[LOCAL_ID] += localKE[LOCAL_ID + offset];
        SYNC_THREADS;
    }
    
    // Write result for this work group
    if (LOCAL_ID == 0)
        kineticEnergy[GROUP_ID] = localKE[0];
}

/**
 * Apply Bussi stochastic velocity rescaling to centroid mode only.
 * This kernel scales the centroid component of each bead's velocity by alpha,
 * leaving internal mode velocities unchanged.
 *
 * Formula: new_vel[bead] = old_vel[bead] + (alpha - 1) * centroid_vel
 */
KERNEL void applyBussiScaling(GLOBAL mixed4* velm, mixed alpha) {
    const mixed invNumCopies = 1.0f / NUM_COPIES;
    const mixed deltaAlpha = alpha - 1.0f;
    
    // Each thread processes one particle
    for (int particle = GLOBAL_ID; particle < NUM_ATOMS; particle += GLOBAL_SIZE) {
        // Compute centroid velocity
        mixed3 centroidVel = make_mixed3(0.0f, 0.0f, 0.0f);
        mixed mass = 0.0f;
        
        for (int copy = 0; copy < NUM_COPIES; copy++) {
            mixed4 particleVelm = velm[particle + copy*PADDED_NUM_ATOMS];
            centroidVel.x += particleVelm.x;
            centroidVel.y += particleVelm.y;
            centroidVel.z += particleVelm.z;
            mass = particleVelm.w;
        }
        
        centroidVel.x *= invNumCopies;
        centroidVel.y *= invNumCopies;
        centroidVel.z *= invNumCopies;
        
        // Apply scaling to all beads
        if (mass != 0.0f) {
            mixed3 delta = centroidVel * deltaAlpha;
            for (int copy = 0; copy < NUM_COPIES; copy++) {
                int idx = particle + copy*PADDED_NUM_ATOMS;
                mixed4 particleVelm = velm[idx];
                particleVelm.x += delta.x;
                particleVelm.y += delta.y;
                particleVelm.z += delta.z;
                velm[idx] = particleVelm;
            }
        }
    }
}

/**
 * Copy a set of positions and velocities from the integrator's arrays to the context.
 */
KERNEL void copyDataToContext(GLOBAL mixed4* srcVel, GLOBAL mixed4* dstVel, GLOBAL mixed4* srcPos, GLOBAL real4* dstPos, GLOBAL int* order, int copy) {
    const int base = copy*PADDED_NUM_ATOMS;
    for (int particle = GLOBAL_ID; particle < NUM_ATOMS; particle += GLOBAL_SIZE) {
        int index = base+order[particle];
        dstVel[particle] = srcVel[index];
        mixed4 posq = srcPos[index];
        dstPos[particle].x = posq.x;
        dstPos[particle].y = posq.y;
        dstPos[particle].z = posq.z;
    }
}

/**
 * Copy a set of positions, velocities, and forces from the context to the integrator's arrays.
 */
KERNEL void copyDataFromContext(GLOBAL mm_long* srcForce, GLOBAL mm_long* dstForce, GLOBAL mixed4* srcVel, GLOBAL mixed4* dstVel,
        GLOBAL real4* srcPos, GLOBAL mixed4* dstPos, GLOBAL int* order, int copy) {
    const int base = copy*PADDED_NUM_ATOMS;
    for (int particle = GLOBAL_ID; particle < NUM_ATOMS; particle += GLOBAL_SIZE) {
        int index = order[particle];
        dstForce[base*3+index] = srcForce[particle];
        dstForce[base*3+index+PADDED_NUM_ATOMS] = srcForce[particle+PADDED_NUM_ATOMS];
        dstForce[base*3+index+PADDED_NUM_ATOMS*2] = srcForce[particle+PADDED_NUM_ATOMS*2];
        dstVel[base+index] = srcVel[particle];
        real4 posq = srcPos[particle];
        dstPos[base+index].x = posq.x;
        dstPos[base+index].y = posq.y;
        dstPos[base+index].z = posq.z;

    }
}

/**
 * Add forces from the context to the integrator's force arrays (for batched RPMD path).
 * Unlike copyDataFromContext, this adds to existing forces instead of overwriting.
 */
KERNEL void addForcesFromContext(GLOBAL mm_long* srcForce, GLOBAL mm_long* dstForce, GLOBAL int* order, int copy) {
    const int base = copy*PADDED_NUM_ATOMS;
    for (int particle = GLOBAL_ID; particle < NUM_ATOMS; particle += GLOBAL_SIZE) {
        int index = order[particle];
        dstForce[base*3+index] += srcForce[particle];
        dstForce[base*3+index+PADDED_NUM_ATOMS] += srcForce[particle+PADDED_NUM_ATOMS];
        dstForce[base*3+index+PADDED_NUM_ATOMS*2] += srcForce[particle+PADDED_NUM_ATOMS*2];
    }
}

/**
 * Atom positions in one copy have been modified.  Apply the same offsets to all the other copies.
 */
KERNEL void applyCellTranslations(GLOBAL mixed4* posq, GLOBAL real4* movedPos, GLOBAL int* order, int movedCopy) {
    for (int particle = GLOBAL_ID; particle < NUM_ATOMS; particle += GLOBAL_SIZE) {
        int index = order[particle];
        real4 p = movedPos[particle];
        mixed4 delta = make_mixed4(p.x, p.y, p.z, p.w)-posq[movedCopy*PADDED_NUM_ATOMS+index];
        for (int copy = 0; copy < NUM_COPIES; copy++) {
            posq[copy*PADDED_NUM_ATOMS+index].x += delta.x;
            posq[copy*PADDED_NUM_ATOMS+index].y += delta.y;
            posq[copy*PADDED_NUM_ATOMS+index].z += delta.z;
        }
    }
}

// ============================================================================
// HYBRID MODE - ALL-GPU UNIFIED KERNELS
// ============================================================================
// These kernels handle both quantum and classical particles efficiently on GPU.
// All particles stored with all beads, but classical particles have coupled beads.

/**
 * Synchronize classical particle beads: average velocities and copy position from bead 0.
 * Classical particles should have identical positions across all beads, and velocities
 * should be synchronized by averaging (centroid velocity).
 * Called after velocity updates to maintain proper classical dynamics.
 */
KERNEL void syncClassicalBeads(
    GLOBAL mixed4* posq,
    GLOBAL mixed4* velm,
    GLOBAL const int* isQuantum
) {
    // Each thread handles one particle
    for (int particle = GLOBAL_ID; particle < NUM_ATOMS; particle += GLOBAL_SIZE) {
        // Skip quantum particles - they have independent beads
        if (isQuantum[particle] != 0)
            continue;
        
        // For classical particles:
        // 1. Positions should all be identical (use bead 0)
        // 2. Velocities should be averaged to get centroid velocity
        
        mixed4 pos0 = posq[particle];  // Bead 0 position
        
        // Compute average velocity across all beads
        mixed3 vel_avg = make_mixed3(0, 0, 0);
        mixed invMass = 0;
        for (int copy = 0; copy < NUM_COPIES; copy++) {
            int idx = particle + copy * PADDED_NUM_ATOMS;
            mixed4 vel_bead = velm[idx];
            vel_avg.x += vel_bead.x;
            vel_avg.y += vel_bead.y;
            vel_avg.z += vel_bead.z;
            if (copy == 0)
                invMass = vel_bead.w;
        }
        vel_avg.x /= NUM_COPIES;
        vel_avg.y /= NUM_COPIES;
        vel_avg.z /= NUM_COPIES;
        
        // Set all beads to the same position and averaged velocity
        for (int copy = 0; copy < NUM_COPIES; copy++) {
            int idx = particle + copy * PADDED_NUM_ATOMS;
            posq[idx] = pos0;
            velm[idx] = make_mixed4(vel_avg.x, vel_avg.y, vel_avg.z, invMass);
        }
    }
}

/**
 * Hybrid PILE thermostat: applies PILE to quantum particles only.
 * Classical particles are skipped (they use a separate thermostat).
 */
KERNEL void applyPileThermostatHybrid(
    GLOBAL mixed4* velm,
    GLOBAL float4* random,
    GLOBAL const int* isQuantum,
    unsigned int randomIndex,
    mixed dt,
    mixed kT,
    mixed friction,
    int applyToCentroid
) {
    const int numBlocks = GLOBAL_SIZE/NUM_COPIES;
    const int blockStart = NUM_COPIES*(LOCAL_ID/NUM_COPIES);
    const int indexInBlock = LOCAL_ID-blockStart;
    const mixed nkT = NUM_COPIES*kT;
    const mixed twown = 2.0f*nkT/HBAR;
    const mixed c1_0 = exp(-0.5f*dt*friction);
    const mixed c2_0 = sqrt(1.0f-c1_0*c1_0);
    LOCAL mixed3 v[2*THREAD_BLOCK_SIZE];
    LOCAL mixed3 temp[2*THREAD_BLOCK_SIZE];
    LOCAL mixed2 w[NUM_COPIES];
    LOCAL_ARG mixed3* vreal = &v[blockStart];
    LOCAL_ARG mixed3* vimag = &v[blockStart+LOCAL_SIZE];
    if (LOCAL_ID < NUM_COPIES)
        w[indexInBlock] = make_mixed2(cos(-indexInBlock*2*M_PI/NUM_COPIES), sin(-indexInBlock*2*M_PI/NUM_COPIES));
    SYNC_THREADS;
    randomIndex += NUM_COPIES*((GLOBAL_ID)/NUM_COPIES);
    
    for (int particle = (GLOBAL_ID)/NUM_COPIES; particle < NUM_ATOMS; particle += numBlocks) {
        // Skip classical particles - they don't use PILE thermostat
        if (isQuantum[particle] == 0) {
            randomIndex += GLOBAL_SIZE;
            continue;
        }
        
        mixed4 particleVelm = velm[particle+indexInBlock*PADDED_NUM_ATOMS];
        mixed invMass = particleVelm.w;
        mixed c3_0 = c2_0*sqrt(nkT*invMass);
        
        // Forward FFT
        vreal[indexInBlock] = SCALE*make_mixed3(particleVelm.x, particleVelm.y, particleVelm.z);
        vimag[indexInBlock] = make_mixed3(0);
        SYNC_THREADS;
        FFT_V_FORWARD

        // Apply the thermostat
        if (indexInBlock == 0) {
            if (applyToCentroid) {
                float4 rand = random[randomIndex];
                vreal[0] = vreal[0]*c1_0 + c3_0*make_mixed3(rand.x, rand.y, rand.z);
            }
        }
        else {
            int k = (indexInBlock <= NUM_COPIES/2 ? indexInBlock : NUM_COPIES-indexInBlock);
            const bool isCenter = (NUM_COPIES%2 == 0 && k == NUM_COPIES/2);
            const mixed wk = twown*sin(k*M_PI/NUM_COPIES);
            const mixed c1 = exp(-wk*dt);
            const mixed c2 = sqrt((1.0f-c1*c1)/2.0f) * (isCenter ? sqrt(2.0f) : 1.0f);
            const mixed c3 = c2*sqrt(nkT*invMass);
            float4 rand1 = random[randomIndex+k];
            float4 rand2 = (isCenter ? make_float4(0) : random[randomIndex+NUM_COPIES-k]);
            vreal[indexInBlock] = c1*vreal[indexInBlock] + c3*make_mixed3(rand1.x, rand1.y, rand1.z);
            vimag[indexInBlock] = c1*vimag[indexInBlock] + c3*(indexInBlock < NUM_COPIES/2 ? make_mixed3(rand2.x, rand2.y, rand2.z) : make_mixed3(-rand2.x, -rand2.y, -rand2.z));
        }
        SYNC_THREADS;
        
        // Inverse FFT
        FFT_V_BACKWARD
        if (invMass != 0)
            velm[particle+indexInBlock*PADDED_NUM_ATOMS] = make_mixed4(SCALE*vreal[indexInBlock].x, SCALE*vreal[indexInBlock].y, SCALE*vreal[indexInBlock].z, particleVelm.w);
        randomIndex += GLOBAL_SIZE;
    }
}

/**
 * Hybrid integration step: quantum particles use FFT ring polymer propagation,
 * classical particles use simple velocity Verlet.
 */
KERNEL void integrateStepHybrid(
    GLOBAL mixed4* posq,
    GLOBAL mixed4* velm,
    GLOBAL mm_long* force,
    GLOBAL const int* isQuantum,
    mixed dt,
    mixed kT
) {
    const int numBlocks = (GLOBAL_SIZE)/NUM_COPIES;
    const int blockStart = NUM_COPIES*(LOCAL_ID/NUM_COPIES);
    const int indexInBlock = LOCAL_ID-blockStart;
    const mixed nkT = NUM_COPIES*kT;
    const mixed twown = 2.0f*nkT/HBAR;
    const mixed forceScale = 1/(mixed) 0x100000000;
    LOCAL mixed3 q[2*THREAD_BLOCK_SIZE];
    LOCAL mixed3 v[2*THREAD_BLOCK_SIZE];
    LOCAL mixed3 temp[2*THREAD_BLOCK_SIZE];
    LOCAL mixed2 w[NUM_COPIES];

    // Update velocities (first half-step) for all particles
    for (int particle = (GLOBAL_ID)/NUM_COPIES; particle < NUM_ATOMS; particle += numBlocks) {
        int index = particle+indexInBlock*PADDED_NUM_ATOMS;
        int forceIndex = particle+indexInBlock*PADDED_NUM_ATOMS*3;
        mixed4 particleVelm = velm[index];
        particleVelm.x += forceScale*force[forceIndex]*(0.5f*dt*particleVelm.w);
        particleVelm.y += forceScale*force[forceIndex+PADDED_NUM_ATOMS]*(0.5f*dt*particleVelm.w);
        particleVelm.z += forceScale*force[forceIndex+PADDED_NUM_ATOMS*2]*(0.5f*dt*particleVelm.w);
        if (particleVelm.w != 0)
            velm[index] = particleVelm;
    }
    
    // Ring polymer propagation - different for quantum vs classical
    LOCAL_ARG mixed3* qreal = &q[blockStart];
    LOCAL_ARG mixed3* qimag = &q[blockStart+LOCAL_SIZE];
    LOCAL_ARG mixed3* vreal = &v[blockStart];
    LOCAL_ARG mixed3* vimag = &v[blockStart+LOCAL_SIZE];
    if (LOCAL_ID < NUM_COPIES)
        w[indexInBlock] = make_mixed2(cos(-indexInBlock*2*M_PI/NUM_COPIES), sin(-indexInBlock*2*M_PI/NUM_COPIES));
    SYNC_THREADS;
    
    for (int particle = (GLOBAL_ID)/NUM_COPIES; particle < NUM_ATOMS; particle += numBlocks) {
        mixed4 particlePosq = posq[particle+indexInBlock*PADDED_NUM_ATOMS];
        mixed4 particleVelm = velm[particle+indexInBlock*PADDED_NUM_ATOMS];
        
        if (isQuantum[particle] != 0) {
            // QUANTUM PARTICLE: Full FFT-based ring polymer propagation
            
            // Forward FFT
            qreal[indexInBlock] = SCALE*make_mixed3(particlePosq.x, particlePosq.y, particlePosq.z);
            qimag[indexInBlock] = make_mixed3(0);
            vreal[indexInBlock] = SCALE*make_mixed3(particleVelm.x, particleVelm.y, particleVelm.z);
            vimag[indexInBlock] = make_mixed3(0);
            SYNC_THREADS;
            FFT_Q_FORWARD
            FFT_V_FORWARD

            // Ring polymer evolution in normal mode space
            if (indexInBlock == 0) {
                qreal[0] += vreal[0]*dt;
                qimag[0] += vimag[0]*dt;
            }
            else {
                const mixed wk = twown*sin(indexInBlock*M_PI/NUM_COPIES);
                const mixed wt = wk*dt;
                const mixed coswt = cos(wt);
                const mixed sinwt = sin(wt);
                const mixed3 vprimereal = vreal[indexInBlock]*coswt - qreal[indexInBlock]*(wk*sinwt);
                const mixed3 vprimeimag = vimag[indexInBlock]*coswt - qimag[indexInBlock]*(wk*sinwt);
                qreal[indexInBlock] = vreal[indexInBlock]*(sinwt/wk) + qreal[indexInBlock]*coswt;
                qimag[indexInBlock] = vimag[indexInBlock]*(sinwt/wk) + qimag[indexInBlock]*coswt;
                vreal[indexInBlock] = vprimereal;
                vimag[indexInBlock] = vprimeimag;
            }
            SYNC_THREADS;
            
            // Inverse FFT
            FFT_Q_BACKWARD
            FFT_V_BACKWARD
            if (particleVelm.w != 0) {
                posq[particle+indexInBlock*PADDED_NUM_ATOMS] = make_mixed4(SCALE*qreal[indexInBlock].x, SCALE*qreal[indexInBlock].y, SCALE*qreal[indexInBlock].z, particlePosq.w);
                velm[particle+indexInBlock*PADDED_NUM_ATOMS] = make_mixed4(SCALE*vreal[indexInBlock].x, SCALE*vreal[indexInBlock].y, SCALE*vreal[indexInBlock].z, particleVelm.w);
            }
        }
        else {
            // CLASSICAL PARTICLE: Use centroid-averaged velocity for position update
            // All beads must be updated with the same position (centroid position)
            
            // Store velocity in shared memory for averaging
            vreal[indexInBlock] = make_mixed3(particleVelm.x, particleVelm.y, particleVelm.z);
            SYNC_THREADS;
            
            // Compute centroid velocity (average across all beads)
            if (particleVelm.w != 0) {
                mixed3 vel_centroid = make_mixed3(0, 0, 0);
                for (int k = 0; k < NUM_COPIES; k++) {
                    vel_centroid.x += vreal[k].x;
                    vel_centroid.y += vreal[k].y;
                    vel_centroid.z += vreal[k].z;
                }
                vel_centroid.x /= NUM_COPIES;
                vel_centroid.y /= NUM_COPIES;
                vel_centroid.z /= NUM_COPIES;
                
                // Update position using centroid velocity: r_c(t+dt) = r_c(t) + v_c(t+dt/2) * dt
                // Only bead 0 needs position update; sync will copy to all beads
                if (indexInBlock == 0) {
                    particlePosq.x += vel_centroid.x * dt;
                    particlePosq.y += vel_centroid.y * dt;
                    particlePosq.z += vel_centroid.z * dt;
                    posq[particle] = particlePosq;
                }
            }
        }
        SYNC_THREADS;
    }
}

/**
 * Hybrid velocity advance: second half-step for all particles.
 */
KERNEL void advanceVelocitiesHybrid(
    GLOBAL mixed4* velm,
    GLOBAL mm_long* force,
    GLOBAL const int* isQuantum,
    mixed dt
) {
    const int numBlocks = (GLOBAL_SIZE)/NUM_COPIES;
    const int blockStart = NUM_COPIES*(LOCAL_ID/NUM_COPIES);
    const int indexInBlock = LOCAL_ID-blockStart;
    const mixed forceScale = 1/(mixed) 0x100000000;

    for (int particle = (GLOBAL_ID)/NUM_COPIES; particle < NUM_ATOMS; particle += numBlocks) {
        // Update all beads for all particles (quantum and classical)
        // Classical particles will have different forces on each bead (from centroid-dependent potential)
        int index = particle+indexInBlock*PADDED_NUM_ATOMS;
        int forceIndex = particle+indexInBlock*PADDED_NUM_ATOMS*3;
        mixed4 particleVelm = velm[index];
        particleVelm.x += forceScale*force[forceIndex]*(0.5f*dt*particleVelm.w);
        particleVelm.y += forceScale*force[forceIndex+PADDED_NUM_ATOMS]*(0.5f*dt*particleVelm.w);
        particleVelm.z += forceScale*force[forceIndex+PADDED_NUM_ATOMS*2]*(0.5f*dt*particleVelm.w);
        if (particleVelm.w != 0)
            velm[index] = particleVelm;
    }
}

/**
 * Apply classical thermostat: Langevin or Bussi for classical particles only.
 * This operates on bead 0 only; syncClassicalBeads copies to other beads.
 */
KERNEL void applyClassicalThermostat(
    GLOBAL mixed4* velm,
    GLOBAL float4* random,
    GLOBAL const int* isQuantum,
    unsigned int randomIndex,
    mixed dt,
    mixed kT,
    mixed friction,
    int thermostatType  // 0=Bussi, 1=Langevin, 2=None
) {
    if (thermostatType == 2) return;  // No thermostat
    
    const mixed c1 = exp(-0.5f*dt*friction);
    const mixed c2 = sqrt(1.0f - c1*c1);
    
    for (int particle = GLOBAL_ID; particle < NUM_ATOMS; particle += GLOBAL_SIZE) {
        // Skip quantum particles
        if (isQuantum[particle] != 0)
            continue;
        
        mixed4 velocity = velm[particle];  // Bead 0
        mixed invMass = velocity.w;
        
        if (invMass == 0.0f) continue;
        
        if (thermostatType == 1) {  // Langevin
            mixed c3 = c2 * sqrt(kT * invMass);
            float4 rand = random[randomIndex + particle];
            
            velocity.x = c1 * velocity.x + c3 * rand.x;
            velocity.y = c1 * velocity.y + c3 * rand.y;
            velocity.z = c1 * velocity.z + c3 * rand.z;
            
            velm[particle] = velocity;
        }
        // Bussi thermostat requires global kinetic energy computation, done on CPU
    }
}

/**
 * Compute kinetic energy of classical particles for Bussi thermostat (isQuantum version).
 * Each work item contributes to the total KE, with parallel reduction.
 */
KERNEL void computeClassicalKEHybrid(
    GLOBAL mixed4* velm,
    GLOBAL const int* isQuantum,
    GLOBAL mixed* kineticEnergy
) {
    LOCAL mixed localKE[THREAD_BLOCK_SIZE];
    mixed myKE = 0.0f;
    
    for (int particle = GLOBAL_ID; particle < NUM_ATOMS; particle += GLOBAL_SIZE) {
        // Skip quantum particles
        if (isQuantum[particle] != 0)
            continue;
            
        mixed4 velocity = velm[particle];  // Bead 0
        mixed invMass = velocity.w;
        
        if (invMass != 0.0f) {
            mixed v2 = velocity.x*velocity.x + velocity.y*velocity.y + velocity.z*velocity.z;
            myKE += 0.5f * v2 / invMass;
        }
    }
    
    // Parallel reduction within work group
    localKE[LOCAL_ID] = myKE;
    SYNC_THREADS;
    
    for (int offset = 1; offset < THREAD_BLOCK_SIZE; offset *= 2) {
        if (LOCAL_ID + offset < THREAD_BLOCK_SIZE && LOCAL_ID % (2*offset) == 0)
            localKE[LOCAL_ID] += localKE[LOCAL_ID + offset];
        SYNC_THREADS;
    }
    
    // Write result for this work group
    if (LOCAL_ID == 0)
        kineticEnergy[GROUP_ID] = localKE[0];
}

/**
 * Apply Bussi scaling to classical particles (isQuantum version).
 */
KERNEL void applyBussiClassicalHybrid(
    GLOBAL mixed4* velm,
    GLOBAL const int* isQuantum,
    mixed scalingFactor
) {
    for (int particle = GLOBAL_ID; particle < NUM_ATOMS; particle += GLOBAL_SIZE) {
        // Skip quantum particles
        if (isQuantum[particle] != 0)
            continue;
        
        mixed4 velocity = velm[particle];  // Bead 0
        
        if (velocity.w != 0.0f) {
            velocity.x *= scalingFactor;
            velocity.y *= scalingFactor;
            velocity.z *= scalingFactor;
            velm[particle] = velocity;
        }
    }
}

