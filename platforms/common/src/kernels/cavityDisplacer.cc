/**
 * Clear the dipole buffer (faster than uploading zeros from CPU).
 */
KERNEL void clearDipoleBuffer(GLOBAL float* RESTRICT dipole) {
    if (GLOBAL_ID == 0) {
        dipole[0] = 0.0f;
        dipole[1] = 0.0f;
        dipole[2] = 0.0f;
        dipole[3] = 0.0f;
    }
}

/**
 * Compute the molecular dipole moment for displacement calculation.
 * This kernel is shared with cavityForce.cc but duplicated here for modularity.
 */
KERNEL void computeCavityDipole(GLOBAL const real4* RESTRICT posq, GLOBAL const float* RESTRICT charges,
        GLOBAL float* RESTRICT dipole, int cavityParticleIndex) {
    LOCAL float localDipoleX[WORK_GROUP_SIZE];
    LOCAL float localDipoleY[WORK_GROUP_SIZE];
    LOCAL float localDipoleZ[WORK_GROUP_SIZE];
    int localId = LOCAL_ID;
    
    float dx = 0.0f, dy = 0.0f, dz = 0.0f;
    for (int i = GLOBAL_ID; i < NUM_ATOMS; i += GLOBAL_SIZE) {
        if (i != cavityParticleIndex) {
            real4 pos = posq[i];
            float q = charges[i];
            dx += q * pos.x;
            dy += q * pos.y;
            dz += q * pos.z;
        }
    }
    
    localDipoleX[localId] = dx;
    localDipoleY[localId] = dy;
    localDipoleZ[localId] = dz;
    SYNC_THREADS;
    
    // Tree reduction
    for (int stride = WORK_GROUP_SIZE/2; stride > 0; stride >>= 1) {
        if (localId < stride) {
            localDipoleX[localId] += localDipoleX[localId + stride];
            localDipoleY[localId] += localDipoleY[localId + stride];
            localDipoleZ[localId] += localDipoleZ[localId + stride];
        }
        SYNC_THREADS;
    }
    
    if (localId == 0) {
        ATOMIC_ADD(&dipole[0], localDipoleX[0]);
        ATOMIC_ADD(&dipole[1], localDipoleY[0]);
        ATOMIC_ADD(&dipole[2], localDipoleZ[0]);
    }
}

/**
 * Displace the cavity particle to its equilibrium position.
 * 
 * The equilibrium position is: q_eq = -(lambda/omega) * d_xy
 * 
 * Only x,y components are modified; z is preserved.
 */
KERNEL void displaceCavityParticle(GLOBAL real4* RESTRICT posq, GLOBAL const float* RESTRICT dipole,
        int cavityParticleIndex, float factor) {
    if (GLOBAL_ID == 0) {
        real4 pos = posq[cavityParticleIndex];
        
        // Read dipole
        float dipoleX = dipole[0];
        float dipoleY = dipole[1];
        
        // Compute new position: q = factor * d_xy (factor = -lambda/omega)
        pos.x = factor * dipoleX;
        pos.y = factor * dipoleY;
        // Keep z unchanged
        
        posq[cavityParticleIndex] = pos;
    }
}
