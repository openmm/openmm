#define COMPUTING_EXCEPTIONS

/**
 * Compute exceptions for HIPPO.
 */
extern "C" __global__ void computeNonbondedExceptions(
        unsigned long long* __restrict__ forceBuffers, mixed* __restrict__ energyBuffer, unsigned long long* __restrict__ torqueBuffers,
        const real4* __restrict__ posq, const real3* __restrict__ extDipole, const int2* __restrict__ exceptionAtoms, const real* __restrict__ mmScale,
        const real* __restrict__ dmScale, const real* __restrict__ ddScale, const real* __restrict__ dispScale, const real* __restrict__ repScale, const real* __restrict__ ctScale,
        const real* __restrict__ coreCharge, const real* __restrict__ valenceCharge, const real* __restrict__ alpha, const real* __restrict__ epsilon,
        const real* __restrict__ damping, const real* __restrict__ c6, const real* __restrict__ pauliK, const real* __restrict__ pauliQ,
        const real* __restrict__ pauliAlpha, const real3* __restrict__ dipole, const real3* __restrict__ inducedDipole, const real* __restrict__ qXX,
        const real* __restrict__ qXY, const real* __restrict__ qXZ, const real* __restrict__ qYY, const real* __restrict__ qYZ,
        const real3* __restrict__ extrapolatedDipole
#ifdef USE_CUTOFF
        , real4 periodicBoxSize, real4 invPeriodicBoxSize, real4 periodicBoxVecX, real4 periodicBoxVecY, real4 periodicBoxVec
#endif
        ) {
    mixed energy = 0;
    const bool isExcluded = false;
    const real interactionScale = 1.0f;
    for (int index = blockIdx.x*blockDim.x+threadIdx.x; index < NUM_EXCEPTIONS; index += blockDim.x*gridDim.x) {
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
            real3 dipole1 = dipole[atom1];
            real3 inducedDipole1 = inducedDipole[atom1];
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
            real3 dipole2 = dipole[atom2];
            real3 inducedDipole2 = inducedDipole[atom2];
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
            atomicAdd(&forceBuffers[atom1], static_cast<unsigned long long>((long long) (tempForce.x*0x100000000)));
            atomicAdd(&forceBuffers[atom1+PADDED_NUM_ATOMS], static_cast<unsigned long long>((long long) (tempForce.y*0x100000000)));
            atomicAdd(&forceBuffers[atom1+2*PADDED_NUM_ATOMS], static_cast<unsigned long long>((long long) (tempForce.z*0x100000000)));
            atomicAdd(&forceBuffers[atom2], static_cast<unsigned long long>((long long) (-tempForce.x*0x100000000)));
            atomicAdd(&forceBuffers[atom2+PADDED_NUM_ATOMS], static_cast<unsigned long long>((long long) (-tempForce.y*0x100000000)));
            atomicAdd(&forceBuffers[atom2+2*PADDED_NUM_ATOMS], static_cast<unsigned long long>((long long) (-tempForce.z*0x100000000)));
            atomicAdd(&torqueBuffers[atom1], static_cast<unsigned long long>((long long) (tempTorque1.x*0x100000000)));
            atomicAdd(&torqueBuffers[atom1+PADDED_NUM_ATOMS], static_cast<unsigned long long>((long long) (tempTorque1.y*0x100000000)));
            atomicAdd(&torqueBuffers[atom1+2*PADDED_NUM_ATOMS], static_cast<unsigned long long>((long long) (tempTorque1.z*0x100000000)));
            atomicAdd(&torqueBuffers[atom2], static_cast<unsigned long long>((long long) (tempTorque2.x*0x100000000)));
            atomicAdd(&torqueBuffers[atom2+PADDED_NUM_ATOMS], static_cast<unsigned long long>((long long) (tempTorque2.y*0x100000000)));
            atomicAdd(&torqueBuffers[atom2+2*PADDED_NUM_ATOMS], static_cast<unsigned long long>((long long) (tempTorque2.z*0x100000000)));
#ifdef USE_CUTOFF
        }
#endif
    }
    energyBuffer[blockIdx.x*blockDim.x+threadIdx.x] += energy;
}
