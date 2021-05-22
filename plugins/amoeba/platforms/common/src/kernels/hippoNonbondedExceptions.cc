#define COMPUTING_EXCEPTIONS

/**
 * Compute exceptions for HIPPO.
 */
KERNEL void computeNonbondedExceptions(
        GLOBAL mm_ulong* RESTRICT forceBuffers, GLOBAL mixed* RESTRICT energyBuffer, GLOBAL mm_ulong* RESTRICT torqueBuffers,
        GLOBAL const real4* RESTRICT posq, GLOBAL const int2* RESTRICT exceptionAtoms, GLOBAL const real* RESTRICT mmScale,
        GLOBAL const real* RESTRICT dmScale, GLOBAL const real* RESTRICT ddScale, GLOBAL const real* RESTRICT dispScale, GLOBAL const real* RESTRICT repScale, GLOBAL const real* RESTRICT ctScale,
        GLOBAL const real* RESTRICT coreCharge, GLOBAL const real* RESTRICT valenceCharge, GLOBAL const real* RESTRICT alpha, GLOBAL const real* RESTRICT epsilon,
        GLOBAL const real* RESTRICT damping, GLOBAL const real* RESTRICT c6, GLOBAL const real* RESTRICT pauliK, GLOBAL const real* RESTRICT pauliQ,
        GLOBAL const real* RESTRICT pauliAlpha, GLOBAL const real* RESTRICT dipole, GLOBAL const real* RESTRICT inducedDipole, GLOBAL const real* RESTRICT qXX,
        GLOBAL const real* RESTRICT qXY, GLOBAL const real* RESTRICT qXZ, GLOBAL const real* RESTRICT qYY, GLOBAL const real* RESTRICT qYZ,
        GLOBAL const real* RESTRICT extrapolatedDipole
#ifdef USE_CUTOFF
        , real4 periodicBoxSize, real4 invPeriodicBoxSize, real4 periodicBoxVecX, real4 periodicBoxVecY, real4 periodicBoxVec
#endif
        ) {
    mixed energy = 0;
    const bool isExcluded = false;
    const real interactionScale = 1.0f;
    for (int index = GLOBAL_ID; index < NUM_EXCEPTIONS; index += GLOBAL_SIZE) {
        int2 atoms = exceptionAtoms[index];
        int atom1 = atoms.x;
        int atom2 = atoms.y;
        real4 pos1 = posq[atom1];
        real4 pos2 = posq[atom2];
        real3 delta = make_real3(pos2.x-pos1.x, pos2.y-pos1.y, pos2.z-pos1.z);
#ifdef USE_PERIODIC
        APPLY_PERIODIC_TO_DELTA(delta)
#endif
        real r2 = delta.x*delta.x + delta.y*delta.y + delta.z*delta.z;
#ifdef USE_CUTOFF
        if (r2 < CUTOFF_SQUARED) {
#endif
            real coreCharge1 = coreCharge[atom1];
            real valenceCharge1 = valenceCharge[atom1];
            real alpha1 = alpha[atom1];
            real epsilon1 = epsilon[atom1];
            real damping1 = damping[atom1];
            real c61 = c6[atom1];
            real pauliK1 = pauliK[atom1];
            real pauliQ1 = pauliQ[atom1];
            real pauliAlpha1 = pauliAlpha[atom1];
            real3 dipole1 = make_real3(dipole[3*atom1], dipole[3*atom1+1], dipole[3*atom1+2]);
            real3 inducedDipole1 = make_real3(inducedDipole[3*atom1], inducedDipole[3*atom1+1], inducedDipole[3*atom1+2]);
            real qXX1 = qXX[atom1];
            real qXY1 = qXY[atom1];
            real qXZ1 = qXZ[atom1];
            real qYY1 = qYY[atom1];
            real qYZ1 = qYZ[atom1];
            real coreCharge2 = coreCharge[atom2];
            real valenceCharge2 = valenceCharge[atom2];
            real alpha2 = alpha[atom2];
            real epsilon2 = epsilon[atom2];
            real damping2 = damping[atom2];
            real c62 = c6[atom2];
            real pauliK2 = pauliK[atom2];
            real pauliQ2 = pauliQ[atom2];
            real pauliAlpha2 = pauliAlpha[atom2];
            real3 dipole2 = make_real3(dipole[3*atom2], dipole[3*atom2+1], dipole[3*atom2+2]);
            real3 inducedDipole2 = make_real3(inducedDipole[3*atom2], inducedDipole[3*atom2+1], inducedDipole[3*atom2+2]);
            real qXX2 = qXX[atom2];
            real qXY2 = qXY[atom2];
            real qXZ2 = qXZ[atom2];
            real qYY2 = qYY[atom2];
            real qYZ2 = qYZ[atom2];
            real multipoleMultipoleScale = mmScale[index];
            real dipoleMultipoleScale = dmScale[index];
            real dipoleDipoleScale = ddScale[index];
            real repulsionScale = repScale[index];
            real dispersionScale = dispScale[index];
            real chargeTransferScale = ctScale[index];
            real rInv = RSQRT(r2);
            real r = r2*rInv;
            real3 tempForce = make_real3(0);
            real3 tempTorque1 = make_real3(0);
            real3 tempTorque2 = make_real3(0);
            real tempEnergy = 0.0f;
            COMPUTE_INTERACTION
            energy += tempEnergy;
            ATOMIC_ADD(&forceBuffers[atom1], (mm_ulong) ((mm_long) (tempForce.x*0x100000000)));
            ATOMIC_ADD(&forceBuffers[atom1+PADDED_NUM_ATOMS], (mm_ulong) ((mm_long) (tempForce.y*0x100000000)));
            ATOMIC_ADD(&forceBuffers[atom1+2*PADDED_NUM_ATOMS], (mm_ulong) ((mm_long) (tempForce.z*0x100000000)));
            ATOMIC_ADD(&forceBuffers[atom2], (mm_ulong) ((mm_long) (-tempForce.x*0x100000000)));
            ATOMIC_ADD(&forceBuffers[atom2+PADDED_NUM_ATOMS], (mm_ulong) ((mm_long) (-tempForce.y*0x100000000)));
            ATOMIC_ADD(&forceBuffers[atom2+2*PADDED_NUM_ATOMS], (mm_ulong) ((mm_long) (-tempForce.z*0x100000000)));
            ATOMIC_ADD(&torqueBuffers[atom1], (mm_ulong) ((mm_long) (tempTorque1.x*0x100000000)));
            ATOMIC_ADD(&torqueBuffers[atom1+PADDED_NUM_ATOMS], (mm_ulong) ((mm_long) (tempTorque1.y*0x100000000)));
            ATOMIC_ADD(&torqueBuffers[atom1+2*PADDED_NUM_ATOMS], (mm_ulong) ((mm_long) (tempTorque1.z*0x100000000)));
            ATOMIC_ADD(&torqueBuffers[atom2], (mm_ulong) ((mm_long) (tempTorque2.x*0x100000000)));
            ATOMIC_ADD(&torqueBuffers[atom2+PADDED_NUM_ATOMS], (mm_ulong) ((mm_long) (tempTorque2.y*0x100000000)));
            ATOMIC_ADD(&torqueBuffers[atom2+2*PADDED_NUM_ATOMS], (mm_ulong) ((mm_long) (tempTorque2.z*0x100000000)));
#ifdef USE_CUTOFF
        }
#endif
    }
    energyBuffer[GLOBAL_ID] += energy;
}
