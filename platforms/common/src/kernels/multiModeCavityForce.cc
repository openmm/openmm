/**
 * Multi-Mode Fabry-Perot Cavity Force Kernel
 *
 * Implements the Hamiltonian:
 *   H = sum_n [ (1/2)*K_n*q_n^2 + eps_n*f_n(z0)*q_n.d ]
 *       + (1/2) * ( sum_n eps_n^2/K_n * f_n(z0)^2 ) * d^2
 *
 * where for mode n (1-indexed):
 *   omega_n = n * omega_1
 *   lambda_n = sqrt(n) * lambda_1
 *   eps_n = lambda_n * omega_n = n^(3/2) * eps_1
 *   K_n = m * omega_n^2 = n^2 * K_1
 *   f_n(z0) = sin(n * pi * z0 / L)
 *
 * The molecular dipole is computed ONCE and reused for all N modes.
 * The DSE prefactor is precomputed on the host.
 *
 * Per-mode parameters are packed into a float4 array (modeParams):
 *   modeParams[n].x = K_n       (spring constant in OpenMM units)
 *   modeParams[n].y = eps_n     (effective coupling in OpenMM units)
 *   modeParams[n].z = f_n       (spatial profile, dimensionless)
 *   modeParams[n].w = (float) reordered cavity particle index for mode n
 */

/**
 * Clear the dipole and energy buffers.
 */
KERNEL void clearMultiModeBuffers(GLOBAL float* RESTRICT dipole,
        GLOBAL float* RESTRICT energyBuffer) {
    if (GLOBAL_ID == 0) {
        dipole[0] = 0.0f;
        dipole[1] = 0.0f;
        dipole[2] = 0.0f;
        dipole[3] = 0.0f;  // padding
        energyBuffer[0] = 0.0f;  // harmonic total
        energyBuffer[1] = 0.0f;  // coupling total
        energyBuffer[2] = 0.0f;  // DSE total
    }
}

/**
 * Compute the molecular dipole moment, excluding ALL cavity particles.
 *
 * The cavity particle indices are stored in cavityIndices[0..NUM_MODES-1].
 * We use a simple check against a bitmask or loop to skip them.
 */
KERNEL void computeMultiModeDipole(GLOBAL const real4* RESTRICT posq,
        GLOBAL const float* RESTRICT charges,
        GLOBAL float* RESTRICT dipole,
        GLOBAL const int* RESTRICT cavityIndices) {
    LOCAL float localDipoleX[WORK_GROUP_SIZE];
    LOCAL float localDipoleY[WORK_GROUP_SIZE];
    int localId = LOCAL_ID;

    float dx = 0.0f, dy = 0.0f;
    for (int i = GLOBAL_ID; i < NUM_ATOMS; i += GLOBAL_SIZE) {
        // Check if this particle is a cavity particle
        bool isCavity = false;
        for (int m = 0; m < NUM_MODES; m++) {
            if (i == cavityIndices[m]) {
                isCavity = true;
                break;
            }
        }
        if (!isCavity) {
            real4 pos = posq[i];
            float q = charges[i];
            dx += q * pos.x;
            dy += q * pos.y;
        }
    }

    localDipoleX[localId] = dx;
    localDipoleY[localId] = dy;
    SYNC_THREADS;

    // Tree reduction
    for (int stride = WORK_GROUP_SIZE/2; stride > 0; stride >>= 1) {
        if (localId < stride) {
            localDipoleX[localId] += localDipoleX[localId + stride];
            localDipoleY[localId] += localDipoleY[localId + stride];
        }
        SYNC_THREADS;
    }

    if (localId == 0) {
        ATOMIC_ADD(&dipole[0], localDipoleX[0]);
        ATOMIC_ADD(&dipole[1], localDipoleY[0]);
    }
}

/**
 * Compute multi-mode cavity forces and energies.
 *
 * Each molecular particle thread loops over all N modes to accumulate forces.
 * Each cavity particle gets its own force from its respective mode.
 *
 * modeParams[n]: (K_n, eps_n, f_n, reorderedCavIdx_n) as float4
 * dsePrefactor: precomputed (1/2) * sum_n(eps_n^2/K_n * f_n^2) in OpenMM units
 */
KERNEL void computeMultiModeForces(GLOBAL const real4* RESTRICT posq,
        GLOBAL const float* RESTRICT charges,
        GLOBAL mm_ulong* RESTRICT forceBuffers,
        GLOBAL const float* RESTRICT dipole,
        GLOBAL float* RESTRICT energyBuffer,
        GLOBAL const int4* RESTRICT posCellOffsets,
        real4 periodicBoxVecX, real4 periodicBoxVecY, real4 periodicBoxVecZ,
        GLOBAL const float4* RESTRICT modeParams,
        GLOBAL const int* RESTRICT cavityIndices,
        float dsePrefactor,
        int paddedNumAtoms) {

    // Read dipole from global memory (computed by computeMultiModeDipole)
    float dipoleX = dipole[0];
    float dipoleY = dipole[1];

    // Compute energy components (only first thread)
    if (GLOBAL_ID == 0) {
        float harmonicTotal = 0.0f;
        float couplingTotal = 0.0f;

        for (int n = 0; n < NUM_MODES; n++) {
            float4 mp = modeParams[n];
            float K_n = mp.x;
            float eps_n = mp.y;
            float f_n = mp.z;
            int cavIdx = (int) mp.w;

            // Unwrap cavity particle position for this mode
            real4 posWrapped = posq[cavIdx];
            int4 offset = posCellOffsets[cavIdx];
            float qx = posWrapped.x - offset.x * periodicBoxVecX.x - offset.y * periodicBoxVecY.x - offset.z * periodicBoxVecZ.x;
            float qy = posWrapped.y - offset.x * periodicBoxVecX.y - offset.y * periodicBoxVecY.y - offset.z * periodicBoxVecZ.y;
            float qz = posWrapped.z - offset.x * periodicBoxVecX.z - offset.y * periodicBoxVecY.z - offset.z * periodicBoxVecZ.z;

            // Harmonic energy: (1/2) * K_n * q_n^2
            harmonicTotal += 0.5f * K_n * (qx*qx + qy*qy + qz*qz);

            // Coupling energy: eps_n * f_n * q_n . d (x,y only)
            couplingTotal += eps_n * f_n * (qx*dipoleX + qy*dipoleY);
        }

        // Dipole self-energy: dsePrefactor * (d_x^2 + d_y^2)
        float dipoleSelfTotal = dsePrefactor * (dipoleX*dipoleX + dipoleY*dipoleY);

        energyBuffer[0] = harmonicTotal;
        energyBuffer[1] = couplingTotal;
        energyBuffer[2] = dipoleSelfTotal;
    }

    // Force on molecular particles: sum over modes of -eps_n * f_n * charge_i * Dq_n
    // where Dq_n = q_n + (eps_n * f_n / K_n) * d
    // Plus DSE force: -2 * dsePrefactor * charge_i * d
    for (int i = GLOBAL_ID; i < NUM_ATOMS; i += GLOBAL_SIZE) {
        // Check if this particle is a cavity particle
        bool isCavity = false;
        for (int m = 0; m < NUM_MODES; m++) {
            if (i == cavityIndices[m]) {
                isCavity = true;
                break;
            }
        }

        if (!isCavity) {
            float q_i = charges[i];
            float fx = 0.0f;
            float fy = 0.0f;

            for (int n = 0; n < NUM_MODES; n++) {
                float4 mp = modeParams[n];
                float K_n = mp.x;
                float eps_n = mp.y;
                float f_n = mp.z;
                int cavIdx = (int) mp.w;
                float epsf_n = eps_n * f_n;

                // Unwrap cavity particle position for this mode
                real4 posWrapped = posq[cavIdx];
                int4 offset = posCellOffsets[cavIdx];
                float qx = posWrapped.x - offset.x * periodicBoxVecX.x - offset.y * periodicBoxVecY.x - offset.z * periodicBoxVecZ.x;
                float qy = posWrapped.y - offset.x * periodicBoxVecX.y - offset.y * periodicBoxVecY.y - offset.z * periodicBoxVecZ.y;

                // Displaced cavity position: Dq_n = q_n + (eps_n*f_n/K_n) * d
                float epsfOverK_n = epsf_n / K_n;
                float DqX = qx + epsfOverK_n * dipoleX;
                float DqY = qy + epsfOverK_n * dipoleY;

                // Force contribution from mode n: -eps_n*f_n * charge_i * Dq_n
                // The displaced coordinate Dq_n includes the DSE contribution:
                //   sum_n(epsf_n^2/K_n)*q_i*d = 2*dsePrefactor*q_i*d
                fx += -epsf_n * q_i * DqX;
                fy += -epsf_n * q_i * DqY;
            }

            // Note: DSE force on molecular particles is already included via the
            // displaced coordinate formulation (DqX, DqY) above. No explicit DSE
            // force term is needed.

            // Convert to fixed point and add to force buffer
            ATOMIC_ADD(&forceBuffers[i], (mm_ulong) realToFixedPoint((real)fx));
            ATOMIC_ADD(&forceBuffers[i+paddedNumAtoms], (mm_ulong) realToFixedPoint((real)fy));
        }
    }

    // Force on cavity particles: one per mode
    // F_n = -K_n * q_n - eps_n * f_n * d
    if (GLOBAL_ID == 0) {
        for (int n = 0; n < NUM_MODES; n++) {
            float4 mp = modeParams[n];
            float K_n = mp.x;
            float eps_n = mp.y;
            float f_n = mp.z;
            int cavIdx = (int) mp.w;
            float epsf_n = eps_n * f_n;

            // Unwrap cavity particle position for this mode
            real4 posWrapped = posq[cavIdx];
            int4 offset = posCellOffsets[cavIdx];
            float qx = posWrapped.x - offset.x * periodicBoxVecX.x - offset.y * periodicBoxVecY.x - offset.z * periodicBoxVecZ.x;
            float qy = posWrapped.y - offset.x * periodicBoxVecX.y - offset.y * periodicBoxVecY.y - offset.z * periodicBoxVecZ.y;
            float qz = posWrapped.z - offset.x * periodicBoxVecX.z - offset.y * periodicBoxVecY.z - offset.z * periodicBoxVecZ.z;

            real fxCav = -K_n * qx - epsf_n * dipoleX;
            real fyCav = -K_n * qy - epsf_n * dipoleY;
            real fzCav = -K_n * qz;

            ATOMIC_ADD(&forceBuffers[cavIdx], (mm_ulong) realToFixedPoint(fxCav));
            ATOMIC_ADD(&forceBuffers[cavIdx+paddedNumAtoms], (mm_ulong) realToFixedPoint(fyCav));
            ATOMIC_ADD(&forceBuffers[cavIdx+2*paddedNumAtoms], (mm_ulong) realToFixedPoint(fzCav));
        }
    }
}
