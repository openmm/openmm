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
 * Compute the molecular dipole moment (excluding cavity particle).
 * Uses atomic additions to accumulate the dipole.
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
 * Compute cavity forces and energies.
 * 
 * The cavity Hamiltonian is:
 *   H = (1/2) * K * q^2 + epsilon * q_xy . d_xy + (epsilon^2 / 2K) * d_xy^2
 * 
 * where:
 *   - K = photonMass * omegac^2 is the spring constant
 *   - omegac is passed in ATOMIC UNITS (Hartree) and must be converted
 *   - epsilon = lambdaCoupling * omegac_converted is the effective coupling
 * 
 * Unit conversions:
 *   - omegac [a.u.] → K [kJ/(mol·nm²)] = photonMass * omegac² * (2625.5 / 0.0529177²)
 *   - The factor (2625.5 / 0.0529177²) converts Hartree/Bohr² to kJ/(mol·nm²)
 * 
 * Forces:
 *   - On molecular particles: F = -epsilon * charge * Dq
 *     where Dq = q + (lambda/omega) * d (displaced cavity position)
 *   - On cavity particle: F = -K * q - epsilon * d_xy
 * 
 * NOTE: qPhoton is reconstructed from wrapped positions using posCellOffsets,
 * which track periodic crossings in the reordered index space.
 */
KERNEL void computeCavityForces(GLOBAL const real4* RESTRICT posq, GLOBAL const float* RESTRICT charges,
        GLOBAL mm_ulong* RESTRICT forceBuffers, GLOBAL const float* RESTRICT dipole,
        GLOBAL float* RESTRICT energyBuffer, GLOBAL const int4* RESTRICT posCellOffsets,
        real4 periodicBoxVecX, real4 periodicBoxVecY, real4 periodicBoxVecZ,
        int reorderedCavityIndex,
        float omegac, float lambdaCoupling, float photonMass, int paddedNumAtoms) {
    
    // Unit conversion constants
    const float HARTREE_TO_KJMOL = 2625.5f;
    const float BOHR_TO_NM = 0.0529177f;
    const float AMU_TO_AU = 1822.888f;  // 1 amu = 1822.888 electron masses (atomic units)
    const float CONVERSION_FACTOR = HARTREE_TO_KJMOL / (BOHR_TO_NM * BOHR_TO_NM);
    
    // Read dipole from global memory (should be called after computeCavityDipole)
    float dipoleX = dipole[0];
    float dipoleY = dipole[1];
    
    // Read wrapped cavity position and unwrap using posCellOffsets
    real4 posWrapped = posq[reorderedCavityIndex];
    int4 offset = posCellOffsets[reorderedCavityIndex];
    real4 qPhoton;
    qPhoton.x = posWrapped.x - offset.x * periodicBoxVecX.x - offset.y * periodicBoxVecY.x - offset.z * periodicBoxVecZ.x;
    qPhoton.y = posWrapped.y - offset.x * periodicBoxVecX.y - offset.y * periodicBoxVecY.y - offset.z * periodicBoxVecZ.y;
    qPhoton.z = posWrapped.z - offset.x * periodicBoxVecX.z - offset.y * periodicBoxVecY.z - offset.z * periodicBoxVecZ.z;
    
    // Compute constants with proper unit conversion
    // CRITICAL: photonMass is in amu, but omegac is in atomic units (Hartree)
    // We must convert photonMass to atomic units (electron masses) for consistency
    float photonMass_au = photonMass * AMU_TO_AU;
    // K [kJ/(mol·nm²)] = photonMass [a.u.] * omegac² [a.u.²] * CONVERSION_FACTOR
    float K = photonMass_au * omegac * omegac * CONVERSION_FACTOR;
    
    // epsilon needs units such that epsilon * d gives FORCE (kJ/(mol·nm))
    // d has units [e·nm], so epsilon needs [kJ/(mol·nm²·e)]
    // epsilon [kJ/(mol·nm²·e)] = lambdaCoupling [dimensionless] * omegac [a.u.] * (HARTREE_TO_KJMOL / BOHR_TO_NM²)
    float epsilon = lambdaCoupling * omegac * CONVERSION_FACTOR;
    
    // Compute displaced cavity position: Dq = q + (epsilon/K) * d
    // This ensures proper unit conversion: epsilon/K has units [1/e]
    float epsilonOverK = epsilon / K;
    float DqX = qPhoton.x + epsilonOverK * dipoleX;
    float DqY = qPhoton.y + epsilonOverK * dipoleY;
    
    // Compute energy components (only first thread in first block)
    if (GLOBAL_ID == 0) {
        // Harmonic energy: (1/2) * K * q^2
        float harmonic = 0.5f * K * (qPhoton.x*qPhoton.x + qPhoton.y*qPhoton.y + qPhoton.z*qPhoton.z);
        
        // Coupling energy: epsilon * q_xy . d_xy
        float coupling = epsilon * (qPhoton.x*dipoleX + qPhoton.y*dipoleY);
        
        // Dipole self-energy: (epsilon^2 / 2K) * d_xy^2
        float dipoleSelf = 0.5f * epsilon * epsilon / K * (dipoleX*dipoleX + dipoleY*dipoleY);
        
        energyBuffer[0] = harmonic;
        energyBuffer[1] = coupling;
        energyBuffer[2] = dipoleSelf;
    }
    
    // Force on molecular particles: F = -epsilon * charge * Dq (only x,y components)
    // Note: Both posq/charges and forceBuffers are indexed by REORDERED position
    for (int i = GLOBAL_ID; i < NUM_ATOMS; i += GLOBAL_SIZE) {
        if (i != reorderedCavityIndex) {
            float q = charges[i];
            real fx = -epsilon * q * DqX;
            real fy = -epsilon * q * DqY;
            
            // Convert to fixed point and add to force buffer
            ATOMIC_ADD(&forceBuffers[i], (mm_ulong) realToFixedPoint(fx));
            ATOMIC_ADD(&forceBuffers[i+paddedNumAtoms], (mm_ulong) realToFixedPoint(fy));
            // No z component
        }
    }
    
    // Force on cavity particle: F = -K*q - epsilon*d_xy
    // Use reorderedCavityIndex for force buffer access
    if (GLOBAL_ID == 0) {
        real fxCavity = -K * qPhoton.x - epsilon * dipoleX;
        real fyCavity = -K * qPhoton.y - epsilon * dipoleY;
        real fzCavity = -K * qPhoton.z;
        
        ATOMIC_ADD(&forceBuffers[reorderedCavityIndex], (mm_ulong) realToFixedPoint(fxCavity));
        ATOMIC_ADD(&forceBuffers[reorderedCavityIndex+paddedNumAtoms], (mm_ulong) realToFixedPoint(fyCavity));
        ATOMIC_ADD(&forceBuffers[reorderedCavityIndex+2*paddedNumAtoms], (mm_ulong) realToFixedPoint(fzCavity));
    }
}
