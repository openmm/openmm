/**
 * Unit conversion: picoseconds to atomic time units
 * 1 ps = 1e-12 s, 1 a.u. time = 2.418884326509e-17 s
 * So 1 ps = 1e-12 / 2.418884326509e-17 = 4.134137333e4 a.u.
 */
const float PS_TO_AU = 4.134137333e4f;

/**
 * Envelope type constants
 */
#define ENVELOPE_CONSTANT 0
#define ENVELOPE_GAUSSIAN 1
#define ENVELOPE_SQUARE 2
#define ENVELOPE_EXPONENTIAL 3

/**
 * Compute envelope function s(t)
 * 
 * @param time_ps        time in picoseconds
 * @param envelope_type  0=constant, 1=gaussian, 2=square, 3=exponential
 * @param env_param1    For gaussian: peak time μ, for square: start_time, for exponential: decay time τ
 * @param env_param2    For gaussian: width σ, for square: stop_time, for exponential: unused
 * @return envelope value s(t)
 */
DEVICE float computeEnvelope(float time_ps, int envelope_type, 
                              float env_param1, float env_param2) {
    float envelope = 1.0f;
    switch(envelope_type) {
        case ENVELOPE_CONSTANT:
            envelope = 1.0f;
            break;
        case ENVELOPE_GAUSSIAN: {
            // s(t) = exp(-(t - μ)²/(2σ²))
            // param1 = peak time μ, param2 = width σ
            // Use multiplication instead of powf to avoid OpenCL macro conflicts
            float arg = (time_ps - env_param1) / env_param2;
            envelope = expf(-0.5f * arg * arg);
            break;
        }
        case ENVELOPE_SQUARE:
            // s(t) = 1 if t_start <= t <= t_stop, else 0
            // param1 = start_time, param2 = stop_time
            envelope = (time_ps >= env_param1 && time_ps <= env_param2) ? 1.0f : 0.0f;
            break;
        case ENVELOPE_EXPONENTIAL:
            // s(t) = exp(-t/τ) for t >= 0
            // param1 = decay time τ
            envelope = (time_ps >= 0.0f) ? expf(-time_ps / env_param1) : 0.0f;
            break;
    }
    return envelope;
}

/**
 * Compute cavity drive force: f(t) = f_0 * s(t) * cos(ω_d*t + φ)
 * 
 * @param time_ps        time in picoseconds
 * @param f0            driving force amplitude f_0 (atomic units)
 * @param omega_d_au    driving frequency ω_d (atomic units)
 * @param phase         phase offset φ
 * @param envelope_type envelope type
 * @param env_param1    envelope parameter 1
 * @param env_param2    envelope parameter 2
 * @return f(t) in atomic units
 */
DEVICE float computeCavityDrive(float time_ps, float f0, float omega_d_au, 
                                 float phase, int envelope_type,
                                 float env_param1, float env_param2) {
    // Convert time from ps to atomic units
    float time_au = time_ps * PS_TO_AU;
    
    // Envelope function s(t)
    float s_t = computeEnvelope(time_ps, envelope_type, env_param1, env_param2);
    
    // Oscillating term: cos(ω_d*t + φ)
    float cos_term = cosf(omega_d_au * time_au + phase);
    
    // Combined: f(t) = f_0 * s(t) * cos(ω_d*t + φ)
    return f0 * s_t * cos_term;
}

/**
 * Compute direct laser electric field: E_ext(t) = E_0 * s(t) * cos(ω_L*t + φ)
 * 
 * @param time_ps        time in picoseconds
 * @param E0            electric field amplitude E_0 (atomic units)
 * @param omega_L_au    laser frequency ω_L (atomic units)
 * @param phase         phase offset φ
 * @param envelope_type envelope type
 * @param env_param1    envelope parameter 1
 * @param env_param2    envelope parameter 2
 * @return E_ext(t) in atomic units
 */
DEVICE float computeDirectLaserField(float time_ps, float E0, float omega_L_au,
                                      float phase, int envelope_type,
                                      float env_param1, float env_param2) {
    // Convert time from ps to atomic units
    float time_au = time_ps * PS_TO_AU;
    
    // Envelope function s(t)
    float s_t = computeEnvelope(time_ps, envelope_type, env_param1, env_param2);
    
    // Oscillating term: cos(ω_L*t + φ)
    float cos_term = cosf(omega_L_au * time_au + phase);
    
    // Combined: E_ext(t) = E_0 * s(t) * cos(ω_L*t + φ)
    return E0 * s_t * cos_term;
}

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
 * Uses UNWRAPPED positions (posq corrected by posCellOffsets) so that the
 * dipole is continuous across periodic boundaries -- matching cav-hoomd.
 */
KERNEL void computeCavityDipole(GLOBAL const real4* RESTRICT posq, GLOBAL const float* RESTRICT charges,
        GLOBAL float* RESTRICT dipole, GLOBAL const int4* RESTRICT posCellOffsets,
        real4 periodicBoxVecX, real4 periodicBoxVecY, real4 periodicBoxVecZ,
        int cavityParticleIndex) {
    LOCAL float localDipoleX[WORK_GROUP_SIZE];
    LOCAL float localDipoleY[WORK_GROUP_SIZE];
    LOCAL float localDipoleZ[WORK_GROUP_SIZE];
    int localId = LOCAL_ID;
    
    float dx = 0.0f, dy = 0.0f, dz = 0.0f;
    for (int i = GLOBAL_ID; i < NUM_ATOMS; i += GLOBAL_SIZE) {
        if (i != cavityParticleIndex) {
            real4 pos = posq[i];
            int4 offset = posCellOffsets[i];
            float ux = pos.x - offset.x * periodicBoxVecX.x - offset.y * periodicBoxVecY.x - offset.z * periodicBoxVecZ.x;
            float uy = pos.y - offset.x * periodicBoxVecX.y - offset.y * periodicBoxVecY.y - offset.z * periodicBoxVecZ.y;
            float uz = pos.z - offset.x * periodicBoxVecX.z - offset.y * periodicBoxVecY.z - offset.z * periodicBoxVecZ.z;
            float q = charges[i];
            dx += q * ux;
            dy += q * uy;
            dz += q * uz;
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
        float omegac, float lambdaCoupling, float photonMass, int paddedNumAtoms,
        float time_ps,
        float f0, float omega_d, float phase_d, int envelope_type_d, float env_param1_d, float env_param2_d, int cavityDriveEnabled,
        float E0, float omega_L, float phase_L, int envelope_type_L, float env_param1_L, float env_param2_L, int directLaserEnabled) {
    
    // Unit conversion constants (NIST 2018 CODATA, matching Reference platform)
    const float HARTREE_TO_KJMOL = 2625.4996f;
    const float BOHR_TO_NM = 0.052917721f;
    const float AMU_TO_AU = 1822.8885f;
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
    
    // Compute laser fields
    float f_drive = 0.0f;
    if (cavityDriveEnabled) {
        f_drive = computeCavityDrive(time_ps, f0, omega_d, phase_d, 
                                      envelope_type_d, env_param1_d, env_param2_d);
    }
    
    float E_ext = 0.0f;
    if (directLaserEnabled) {
        E_ext = computeDirectLaserField(time_ps, E0, omega_L, phase_L,
                                        envelope_type_L, env_param1_L, env_param2_L);
    }
    
    // Compute energy components (only first thread in first block)
    if (GLOBAL_ID == 0) {
        // Harmonic energy: (1/2) * K * q^2
        float harmonic = 0.5f * K * (qPhoton.x*qPhoton.x + qPhoton.y*qPhoton.y + qPhoton.z*qPhoton.z);
        
        // Coupling energy: epsilon * q_xy . d_xy
        float coupling = epsilon * (qPhoton.x*dipoleX + qPhoton.y*dipoleY);
        
        // Dipole self-energy: (epsilon^2 / 2K) * d_xy^2
        float dipoleSelf = 0.5f * epsilon * epsilon / K * (dipoleX*dipoleX + dipoleY*dipoleY);
        
        // Cavity drive energy: E_drive = -f(t) * q (Case 2)
        float cavityDrive = -f_drive * (qPhoton.x + qPhoton.y);
        
        // Direct laser energy: E_laser = -E_ext(t) * μ (Case 1)
        float directLaser = -E_ext * (dipoleX + dipoleY);
        
        energyBuffer[0] = harmonic;
        energyBuffer[1] = coupling;
        energyBuffer[2] = dipoleSelf;
        energyBuffer[3] = cavityDrive;
        energyBuffer[4] = directLaser;
    }
    
    // Force on molecular particles: F = -epsilon * charge * Dq + E_ext * charge (only x,y components)
    // Note: Both posq/charges and forceBuffers are indexed by REORDERED position
    // The E_ext term is from direct laser-molecule coupling (Case 1)
    for (int i = GLOBAL_ID; i < NUM_ATOMS; i += GLOBAL_SIZE) {
        if (i != reorderedCavityIndex) {
            float q = charges[i];
            real fx = -epsilon * q * DqX;
            real fy = -epsilon * q * DqY;
            
            // Add direct laser-molecule coupling force: F = E_ext(t) * charge
            if (directLaserEnabled) {
                fx += E_ext * q;
                fy += E_ext * q;
            }
            
            // Convert to fixed point and add to force buffer
            ATOMIC_ADD(&forceBuffers[i], (mm_ulong) realToFixedPoint(fx));
            ATOMIC_ADD(&forceBuffers[i+paddedNumAtoms], (mm_ulong) realToFixedPoint(fy));
            // No z component
        }
    }
    
    // Force on cavity particle: F = -K*q - epsilon*d_xy + f(t)
    // The +f(t) term is from cavity-mode driving (Case 2)
    // Use reorderedCavityIndex for force buffer access
    if (GLOBAL_ID == 0) {
        real fxCavity = -K * qPhoton.x - epsilon * dipoleX;
        real fyCavity = -K * qPhoton.y - epsilon * dipoleY;
        real fzCavity = -K * qPhoton.z;
        
        // Add cavity drive force: F = +f(t) (positive sign from H_drive = -f(t)*q)
        if (cavityDriveEnabled) {
            fxCavity += f_drive;
            fyCavity += f_drive;
        }
        
        ATOMIC_ADD(&forceBuffers[reorderedCavityIndex], (mm_ulong) realToFixedPoint(fxCavity));
        ATOMIC_ADD(&forceBuffers[reorderedCavityIndex+paddedNumAtoms], (mm_ulong) realToFixedPoint(fyCavity));
        ATOMIC_ADD(&forceBuffers[reorderedCavityIndex+2*paddedNumAtoms], (mm_ulong) realToFixedPoint(fzCavity));
    }
}
