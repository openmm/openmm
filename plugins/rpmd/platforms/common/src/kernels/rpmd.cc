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
// HYBRID CLASSICAL/QUANTUM RPMD KERNELS
// ============================================================================
//
// The following kernels support hybrid simulations where some particles
// receive quantum (RPMD) treatment while others are treated classically.
//

/**
 * Integrate classical particles using velocity Verlet algorithm (first half-step).
 * This advances velocities by dt/2 and positions by dt.
 */
KERNEL void integrateClassical(
    GLOBAL mixed4* posq,
    GLOBAL mixed4* velm,
    GLOBAL mm_long* force,
    GLOBAL const int* classicalIndices,
    int numClassicalParticles,
    mixed dt
) {
    const mixed forceScale = 1/(mixed) 0x100000000;
    
    for (int i = GLOBAL_ID; i < numClassicalParticles; i += GLOBAL_SIZE) {
        int particle = classicalIndices[i];
        
        // Load current state
        mixed4 velocity = velm[particle];
        mixed4 position = posq[particle];
        mixed invMass = velocity.w;
        
        if (invMass == 0.0f) continue;
        
        // Velocity half-step: v(t+dt/2) = v(t) + f(t)*dt/(2m)
        int forceIndex = particle * 3;
        velocity.x += forceScale * force[forceIndex] * (0.5f * dt) * invMass;
        velocity.y += forceScale * force[forceIndex + PADDED_NUM_ATOMS] * (0.5f * dt) * invMass;
        velocity.z += forceScale * force[forceIndex + PADDED_NUM_ATOMS * 2] * (0.5f * dt) * invMass;
        
        // Position update: r(t+dt) = r(t) + v(t+dt/2)*dt
        position.x += velocity.x * dt;
        position.y += velocity.y * dt;
        position.z += velocity.z * dt;
        
        // Store updated state
        posq[particle] = position;
        velm[particle] = velocity;
    }
}

/**
 * Advance classical particle velocities by half step (second half of velocity Verlet).
 * This completes the velocity update: v(t+dt) = v(t+dt/2) + f(t+dt)*dt/(2m)
 */
KERNEL void advanceClassicalVelocities(
    GLOBAL mixed4* velm,
    GLOBAL mm_long* force,
    GLOBAL const int* classicalIndices,
    int numClassicalParticles,
    mixed dt
) {
    const mixed forceScale = 1/(mixed) 0x100000000;
    
    for (int i = GLOBAL_ID; i < numClassicalParticles; i += GLOBAL_SIZE) {
        int particle = classicalIndices[i];
        
        mixed4 velocity = velm[particle];
        mixed invMass = velocity.w;
        
        if (invMass == 0.0f) continue;
        
        // Velocity advance: v(t+dt) = v(t+dt/2) + f(t+dt)*dt/(2m)
        int forceIndex = particle * 3;
        velocity.x += forceScale * force[forceIndex] * (0.5f * dt) * invMass;
        velocity.y += forceScale * force[forceIndex + PADDED_NUM_ATOMS] * (0.5f * dt) * invMass;
        velocity.z += forceScale * force[forceIndex + PADDED_NUM_ATOMS * 2] * (0.5f * dt) * invMass;
        
        velm[particle] = velocity;
    }
}

/**
 * Apply Bussi stochastic velocity rescaling to classical particles.
 * This applies a pre-computed scaling factor to all velocity components.
 */
KERNEL void applyBussiClassical(
    GLOBAL mixed4* velm,
    GLOBAL const int* classicalIndices,
    int numClassicalParticles,
    mixed scalingFactor
) {
    for (int i = GLOBAL_ID; i < numClassicalParticles; i += GLOBAL_SIZE) {
        int particle = classicalIndices[i];
        
        mixed4 velocity = velm[particle];
        
        if (velocity.w != 0.0f) {
            velocity.x *= scalingFactor;
            velocity.y *= scalingFactor;
            velocity.z *= scalingFactor;
            velm[particle] = velocity;
        }
    }
}

/**
 * Apply Langevin thermostat to classical particles.
 * Uses the Langevin equation with friction and random force.
 */
KERNEL void applyLangevinClassical(
    GLOBAL mixed4* velm,
    GLOBAL float4* random,
    GLOBAL const int* classicalIndices,
    int numClassicalParticles,
    unsigned int randomIndex,
    mixed dt,
    mixed kT,
    mixed friction
) {
    const mixed c1 = exp(-friction * dt);
    const mixed c2 = sqrt(1.0f - c1*c1);
    
    for (int i = GLOBAL_ID; i < numClassicalParticles; i += GLOBAL_SIZE) {
        int particle = classicalIndices[i];
        
        mixed4 velocity = velm[particle];
        mixed invMass = velocity.w;
        
        if (invMass == 0.0f) continue;
        
        mixed c3 = c2 * sqrt(kT * invMass);
        float4 rand = random[randomIndex + i];
        
        velocity.x = c1 * velocity.x + c3 * rand.x;
        velocity.y = c1 * velocity.y + c3 * rand.y;
        velocity.z = c1 * velocity.z + c3 * rand.z;
        
        velm[particle] = velocity;
    }
}

/**
 * Compute kinetic energy of classical particles for Bussi thermostat.
 * Each work item contributes to the total KE, with parallel reduction.
 */
KERNEL void computeClassicalKE(
    GLOBAL mixed4* velm,
    GLOBAL const int* classicalIndices,
    int numClassicalParticles,
    GLOBAL mixed* kineticEnergy
) {
    LOCAL mixed localKE[THREAD_BLOCK_SIZE];
    mixed myKE = 0.0f;
    
    for (int i = GLOBAL_ID; i < numClassicalParticles; i += GLOBAL_SIZE) {
        int particle = classicalIndices[i];
        mixed4 velocity = velm[particle];
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

// ============================================================================
// HYBRID MODE - ALL-GPU UNIFIED KERNELS
// ============================================================================
// These kernels handle both quantum and classical particles efficiently on GPU.
// All particles stored with all beads, but classical particles have coupled beads.

/**
 * Synchronize classical particle beads: copy bead 0 to all other beads.
 * This ensures classical particles have identical positions/velocities across all beads.
 * Called after integration to maintain bead coupling for classical particles.
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
        
        // For classical particles, copy bead 0 to all other beads
        mixed4 pos0 = posq[particle];  // Bead 0 position
        mixed4 vel0 = velm[particle];  // Bead 0 velocity
        
        for (int copy = 1; copy < NUM_COPIES; copy++) {
            int idx = particle + copy * PADDED_NUM_ATOMS;
            posq[idx] = pos0;
            velm[idx] = vel0;
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
            // CLASSICAL PARTICLE: Simple position update on bead 0 only
            // Other beads will be synced by syncClassicalBeads kernel
            
            if (indexInBlock == 0 && particleVelm.w != 0) {
                // r(t+dt) = r(t) + v(t+dt/2) * dt
                particlePosq.x += particleVelm.x * dt;
                particlePosq.y += particleVelm.y * dt;
                particlePosq.z += particleVelm.z * dt;
                posq[particle] = particlePosq;
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
        // For classical particles, only update bead 0
        if (isQuantum[particle] == 0 && indexInBlock != 0)
            continue;
            
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

// ============================================================================
// SPARSE STORAGE KERNELS
// ============================================================================
// Memory layout: [quantum particles × numCopies][classical particles × 1]
// This saves significant memory for systems with many classical particles.
//
// Indexing:
// - storageOffset[particle]: base index for particle's data
// - Quantum particle at copy c: storageOffset[particle] + c
// - Classical particle: storageOffset[particle] (only 1 slot)

/**
 * Apply PILE thermostat to quantum particles only (sparse storage version).
 * Only processes quantum particles; classical particles are handled separately.
 * 
 * @param quantumParticles Array of original particle indices for quantum particles
 * @param numQuantumParticles Number of quantum particles
 * @param storageOffset Maps original particle index to storage base offset
 */
KERNEL void applyPileThermostatSparse(
    GLOBAL mixed4* velm,
    GLOBAL float4* random,
    GLOBAL const int* quantumParticles,
    GLOBAL const int* storageOffset,
    int numQuantumParticles,
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
    
    // Loop over quantum particles only
    for (int qIdx = (GLOBAL_ID)/NUM_COPIES; qIdx < numQuantumParticles; qIdx += numBlocks) {
        int particle = quantumParticles[qIdx];  // Original particle index
        int baseOffset = storageOffset[particle];  // Base storage offset
        int storageIdx = baseOffset + indexInBlock;  // Storage index for this bead
        
        mixed4 particleVelm = velm[storageIdx];
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
            velm[storageIdx] = make_mixed4(SCALE*vreal[indexInBlock].x, SCALE*vreal[indexInBlock].y, SCALE*vreal[indexInBlock].z, particleVelm.w);
        randomIndex += GLOBAL_SIZE;
    }
}

/**
 * Integrate quantum particles using ring polymer dynamics (sparse storage).
 */
KERNEL void integrateStepQuantumSparse(
    GLOBAL mixed4* posq,
    GLOBAL mixed4* velm,
    GLOBAL mm_long* force,
    GLOBAL const int* quantumParticles,
    GLOBAL const int* storageOffset,
    int numQuantumParticles,
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
    LOCAL_ARG mixed3* qreal = &q[blockStart];
    LOCAL_ARG mixed3* qimag = &q[blockStart+LOCAL_SIZE];
    LOCAL_ARG mixed3* vreal = &v[blockStart];
    LOCAL_ARG mixed3* vimag = &v[blockStart+LOCAL_SIZE];
    if (LOCAL_ID < NUM_COPIES)
        w[indexInBlock] = make_mixed2(cos(-indexInBlock*2*M_PI/NUM_COPIES), sin(-indexInBlock*2*M_PI/NUM_COPIES));
    SYNC_THREADS;

    for (int qIdx = (GLOBAL_ID)/NUM_COPIES; qIdx < numQuantumParticles; qIdx += numBlocks) {
        int particle = quantumParticles[qIdx];
        int baseOffset = storageOffset[particle];
        int storageIdx = baseOffset + indexInBlock;
        
        // Force index still uses original particle layout for context forces
        int forceIndex = particle + indexInBlock*PADDED_NUM_ATOMS*3;
        
        mixed4 particlePosq = posq[storageIdx];
        mixed4 particleVelm = velm[storageIdx];
        
        // First half-step velocity update
        particleVelm.x += forceScale*force[forceIndex]*(0.5f*dt*particleVelm.w);
        particleVelm.y += forceScale*force[forceIndex+PADDED_NUM_ATOMS]*(0.5f*dt*particleVelm.w);
        particleVelm.z += forceScale*force[forceIndex+PADDED_NUM_ATOMS*2]*(0.5f*dt*particleVelm.w);
        
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
            posq[storageIdx] = make_mixed4(SCALE*qreal[indexInBlock].x, SCALE*qreal[indexInBlock].y, SCALE*qreal[indexInBlock].z, particlePosq.w);
            velm[storageIdx] = make_mixed4(SCALE*vreal[indexInBlock].x, SCALE*vreal[indexInBlock].y, SCALE*vreal[indexInBlock].z, particleVelm.w);
        }
        SYNC_THREADS;
    }
}

/**
 * Integrate classical particles using simple velocity Verlet (sparse storage).
 * Classical particles have only 1 storage slot each.
 */
KERNEL void integrateStepClassicalSparse(
    GLOBAL mixed4* posq,
    GLOBAL mixed4* velm,
    GLOBAL mm_long* force,
    GLOBAL const int* classicalParticles,
    GLOBAL const int* storageOffset,
    int numClassicalParticles,
    mixed dt
) {
    const mixed forceScale = 1/(mixed) 0x100000000;
    
    for (int cIdx = GLOBAL_ID; cIdx < numClassicalParticles; cIdx += GLOBAL_SIZE) {
        int particle = classicalParticles[cIdx];
        int storageIdx = storageOffset[particle];
        
        // Force from bead 0 (context forces use original indexing)
        int forceIndex = particle;
        
        mixed4 position = posq[storageIdx];
        mixed4 velocity = velm[storageIdx];
        mixed invMass = velocity.w;
        
        if (invMass == 0.0f) continue;
        
        // Velocity half-step
        velocity.x += forceScale * force[forceIndex] * (0.5f * dt) * invMass;
        velocity.y += forceScale * force[forceIndex + PADDED_NUM_ATOMS] * (0.5f * dt) * invMass;
        velocity.z += forceScale * force[forceIndex + PADDED_NUM_ATOMS * 2] * (0.5f * dt) * invMass;
        
        // Position update
        position.x += velocity.x * dt;
        position.y += velocity.y * dt;
        position.z += velocity.z * dt;
        
        posq[storageIdx] = position;
        velm[storageIdx] = velocity;
    }
}

/**
 * Second half-step velocity update for quantum particles (sparse storage).
 */
KERNEL void advanceVelocitiesQuantumSparse(
    GLOBAL mixed4* velm,
    GLOBAL mm_long* force,
    GLOBAL const int* quantumParticles,
    GLOBAL const int* storageOffset,
    int numQuantumParticles,
    mixed dt
) {
    const int numBlocks = (GLOBAL_SIZE)/NUM_COPIES;
    const int indexInBlock = LOCAL_ID % NUM_COPIES;
    const mixed forceScale = 1/(mixed) 0x100000000;

    for (int qIdx = (GLOBAL_ID)/NUM_COPIES; qIdx < numQuantumParticles; qIdx += numBlocks) {
        int particle = quantumParticles[qIdx];
        int storageIdx = storageOffset[particle] + indexInBlock;
        int forceIndex = particle + indexInBlock*PADDED_NUM_ATOMS*3;
        
        mixed4 particleVelm = velm[storageIdx];
        particleVelm.x += forceScale*force[forceIndex]*(0.5f*dt*particleVelm.w);
        particleVelm.y += forceScale*force[forceIndex+PADDED_NUM_ATOMS]*(0.5f*dt*particleVelm.w);
        particleVelm.z += forceScale*force[forceIndex+PADDED_NUM_ATOMS*2]*(0.5f*dt*particleVelm.w);
        if (particleVelm.w != 0)
            velm[storageIdx] = particleVelm;
    }
}

/**
 * Second half-step velocity update for classical particles (sparse storage).
 */
KERNEL void advanceVelocitiesClassicalSparse(
    GLOBAL mixed4* velm,
    GLOBAL mm_long* force,
    GLOBAL const int* classicalParticles,
    GLOBAL const int* storageOffset,
    int numClassicalParticles,
    mixed dt
) {
    const mixed forceScale = 1/(mixed) 0x100000000;
    
    for (int cIdx = GLOBAL_ID; cIdx < numClassicalParticles; cIdx += GLOBAL_SIZE) {
        int particle = classicalParticles[cIdx];
        int storageIdx = storageOffset[particle];
        int forceIndex = particle;
        
        mixed4 velocity = velm[storageIdx];
        mixed invMass = velocity.w;
        
        if (invMass == 0.0f) continue;
        
        velocity.x += forceScale * force[forceIndex] * (0.5f * dt) * invMass;
        velocity.y += forceScale * force[forceIndex + PADDED_NUM_ATOMS] * (0.5f * dt) * invMass;
        velocity.z += forceScale * force[forceIndex + PADDED_NUM_ATOMS * 2] * (0.5f * dt) * invMass;
        
        velm[storageIdx] = velocity;
    }
}

/**
 * Apply Langevin thermostat to classical particles (sparse storage).
 */
KERNEL void applyClassicalThermostatSparse(
    GLOBAL mixed4* velm,
    GLOBAL float4* random,
    GLOBAL const int* classicalParticles,
    GLOBAL const int* storageOffset,
    int numClassicalParticles,
    unsigned int randomIndex,
    mixed dt,
    mixed kT,
    mixed friction
) {
    const mixed c1 = exp(-0.5f*dt*friction);
    const mixed c2 = sqrt(1.0f - c1*c1);
    
    for (int cIdx = GLOBAL_ID; cIdx < numClassicalParticles; cIdx += GLOBAL_SIZE) {
        int particle = classicalParticles[cIdx];
        int storageIdx = storageOffset[particle];
        
        mixed4 velocity = velm[storageIdx];
        mixed invMass = velocity.w;
        
        if (invMass == 0.0f) continue;
        
        mixed c3 = c2 * sqrt(kT * invMass);
        float4 rand = random[randomIndex + cIdx];
        
        velocity.x = c1 * velocity.x + c3 * rand.x;
        velocity.y = c1 * velocity.y + c3 * rand.y;
        velocity.z = c1 * velocity.z + c3 * rand.z;
        
        velm[storageIdx] = velocity;
    }
}

/**
 * Compute kinetic energy of classical particles for Bussi thermostat (sparse storage).
 */
KERNEL void computeClassicalKESparse(
    GLOBAL mixed4* velm,
    GLOBAL const int* classicalParticles,
    GLOBAL const int* storageOffset,
    int numClassicalParticles,
    GLOBAL mixed* kineticEnergy
) {
    LOCAL mixed localKE[THREAD_BLOCK_SIZE];
    mixed myKE = 0.0f;
    
    for (int cIdx = GLOBAL_ID; cIdx < numClassicalParticles; cIdx += GLOBAL_SIZE) {
        int particle = classicalParticles[cIdx];
        int storageIdx = storageOffset[particle];
        
        mixed4 velocity = velm[storageIdx];
        mixed invMass = velocity.w;
        
        if (invMass != 0.0f) {
            mixed v2 = velocity.x*velocity.x + velocity.y*velocity.y + velocity.z*velocity.z;
            myKE += 0.5f * v2 / invMass;
        }
    }
    
    localKE[LOCAL_ID] = myKE;
    SYNC_THREADS;
    
    for (int offset = 1; offset < THREAD_BLOCK_SIZE; offset *= 2) {
        if (LOCAL_ID + offset < THREAD_BLOCK_SIZE && LOCAL_ID % (2*offset) == 0)
            localKE[LOCAL_ID] += localKE[LOCAL_ID + offset];
        SYNC_THREADS;
    }
    
    if (LOCAL_ID == 0)
        kineticEnergy[GROUP_ID] = localKE[0];
}

/**
 * Apply Bussi scaling to classical particles (sparse storage).
 */
KERNEL void applyBussiClassicalSparse(
    GLOBAL mixed4* velm,
    GLOBAL const int* classicalParticles,
    GLOBAL const int* storageOffset,
    int numClassicalParticles,
    mixed scalingFactor
) {
    for (int cIdx = GLOBAL_ID; cIdx < numClassicalParticles; cIdx += GLOBAL_SIZE) {
        int particle = classicalParticles[cIdx];
        int storageIdx = storageOffset[particle];
        
        mixed4 velocity = velm[storageIdx];
        
        if (velocity.w != 0.0f) {
            velocity.x *= scalingFactor;
            velocity.y *= scalingFactor;
            velocity.z *= scalingFactor;
            velm[storageIdx] = velocity;
        }
    }
}

/**
 * Copy data from sparse RPMD storage to context (for force computation).
 * This expands classical particles to all beads (same position).
 */
KERNEL void copyDataToContextSparse(
    GLOBAL mixed4* srcVelm,
    GLOBAL mixed4* dstVelm,
    GLOBAL const int* storageOffset,
    GLOBAL const int* isQuantum,
    int copy,
    GLOBAL const mixed4* srcPosq,
    GLOBAL const int* atomIndex,
    GLOBAL mixed4* dstPosq
) {
    for (int particle = GLOBAL_ID; particle < NUM_ATOMS; particle += GLOBAL_SIZE) {
        int baseOffset = storageOffset[particle];
        int srcIdx;
        
        if (isQuantum[particle] != 0) {
            // Quantum: read from the specified copy
            srcIdx = baseOffset + copy;
        } else {
            // Classical: always read from the single storage slot
            srcIdx = baseOffset;
        }
        
        int dstIdx = atomIndex[particle];
        dstVelm[dstIdx] = srcVelm[srcIdx];
        
        mixed4 pos = srcPosq[srcIdx];
        dstPosq[dstIdx] = make_mixed4(pos.x, pos.y, pos.z, dstPosq[dstIdx].w);
    }
}

/**
 * Copy data from context back to sparse RPMD storage (after force computation).
 * Forces go to the appropriate storage location based on particle type.
 */
KERNEL void copyDataFromContextSparse(
    GLOBAL const mm_long* srcForce,
    int copy,
    GLOBAL const mixed4* srcVelm,
    GLOBAL mixed4* dstVelm,
    GLOBAL const mixed4* srcPosq,
    GLOBAL mixed4* dstPosq,
    GLOBAL const int* atomIndex,
    GLOBAL mm_long* dstForce,
    GLOBAL const int* storageOffset,
    GLOBAL const int* isQuantum
) {
    for (int particle = GLOBAL_ID; particle < NUM_ATOMS; particle += GLOBAL_SIZE) {
        int srcIdx = atomIndex[particle];
        int baseOffset = storageOffset[particle];
        int dstIdx;
        
        if (isQuantum[particle] != 0) {
            // Quantum: write to the specified copy
            dstIdx = baseOffset + copy;
        } else {
            // Classical: only copy 0 matters
            if (copy != 0) continue;
            dstIdx = baseOffset;
        }
        
        // Copy forces (scale appropriately)
        dstForce[dstIdx*3] = srcForce[srcIdx];
        dstForce[dstIdx*3 + 1] = srcForce[srcIdx + PADDED_NUM_ATOMS];
        dstForce[dstIdx*3 + 2] = srcForce[srcIdx + PADDED_NUM_ATOMS*2];
        
        // Copy velocities
        dstVelm[dstIdx] = srcVelm[srcIdx];
        
        // Copy positions
        mixed4 pos = srcPosq[srcIdx];
        dstPosq[dstIdx] = make_mixed4(pos.x, pos.y, pos.z, dstPosq[dstIdx].w);
    }
}
