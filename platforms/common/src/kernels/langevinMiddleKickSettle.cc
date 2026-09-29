// Experimental CUDA/mixed exact Middle kick + velocity SETTLE.
// The projection body is copied verbatim from integrationUtilities.cc.
#ifndef USE_MIXED_PRECISION
#error LangevinMiddle kick SETTLE fusion requires mixed precision
#endif

inline DEVICE mixed4 loadKickSettlePos(GLOBAL const real4* RESTRICT posq,
        GLOBAL const real4* RESTRICT posqCorrection, int index) {
    real4 pos1 = posq[index];
    real4 pos2 = posqCorrection[index];
    return make_mixed4(pos1.x+(mixed)pos2.x, pos1.y+(mixed)pos2.y, pos1.z+(mixed)pos2.z, pos1.w);
}

inline DEVICE mixed kickSettleStoredDouble(mixed value) {
    // Preserve the original global double materialization boundary. Front-end
    // opacity is not proof of final device bit equivalence: inspect generated
    // instructions and compare numerical results separately before deployment.
    mixed rounded;
    asm volatile("mov.b64 %0, %1;" : "=d"(rounded) : "d"(value));
    return rounded;
}

inline DEVICE mixed4 kickSettleVelocity(int index, int paddedNumAtoms, mixed fscale,
        GLOBAL const mixed4* RESTRICT velm, GLOBAL const mm_long* RESTRICT force) {
    mixed4 velocity = velm[index];
    if (velocity.w != 0.0) {
        velocity.x += fscale*velocity.w*force[index];
        velocity.y += fscale*velocity.w*force[index+paddedNumAtoms];
        velocity.z += fscale*velocity.w*force[index+paddedNumAtoms*2];
        velocity.x = kickSettleStoredDouble(velocity.x);
        velocity.y = kickSettleStoredDouble(velocity.y);
        velocity.z = kickSettleStoredDouble(velocity.z);
    }
    return velocity;
}

KERNEL void integrateLangevinMiddleKickSettle(int numClusters, int paddedNumAtoms,
        GLOBAL const real4* RESTRICT oldPos, GLOBAL const real4* RESTRICT posqCorrection,
        GLOBAL mixed4* RESTRICT velm, GLOBAL const mm_long* RESTRICT force,
        GLOBAL const mixed2* RESTRICT dt, GLOBAL const int4* RESTRICT clusterAtoms,
        int numResidualAtoms, GLOBAL const int* RESTRICT residualAtoms) {
    mixed fscale = dt[0].y/(mixed) 0x100000000;
    for (int index = GLOBAL_ID; index < numClusters; index += GLOBAL_SIZE) {
        int4 atoms = clusterAtoms[index];
        mixed4 v0 = kickSettleVelocity(atoms.x, paddedNumAtoms, fscale, velm, force);
        mixed4 v1 = kickSettleVelocity(atoms.y, paddedNumAtoms, fscale, velm, force);
        mixed4 v2 = kickSettleVelocity(atoms.z, paddedNumAtoms, fscale, velm, force);
        mixed4 apos0 = loadKickSettlePos(oldPos, posqCorrection, atoms.x);
        mixed4 apos1 = loadKickSettlePos(oldPos, posqCorrection, atoms.y);
        mixed4 apos2 = loadKickSettlePos(oldPos, posqCorrection, atoms.z);

        // Compute intermediate quantities: the atom masses, the bond directions, the relative velocities,
        // and the angle cosines and sines.

        mixed mA = 1/v0.w;
        mixed mB = 1/v1.w;
        mixed mC = 1/v2.w;
        mixed3 eAB = make_mixed3(apos1.x-apos0.x, apos1.y-apos0.y, apos1.z-apos0.z);
        mixed3 eBC = make_mixed3(apos2.x-apos1.x, apos2.y-apos1.y, apos2.z-apos1.z);
        mixed3 eCA = make_mixed3(apos0.x-apos2.x, apos0.y-apos2.y, apos0.z-apos2.z);
        eAB *= RSQRT(eAB.x*eAB.x + eAB.y*eAB.y + eAB.z*eAB.z);
        eBC *= RSQRT(eBC.x*eBC.x + eBC.y*eBC.y + eBC.z*eBC.z);
        eCA *= RSQRT(eCA.x*eCA.x + eCA.y*eCA.y + eCA.z*eCA.z);
        mixed vAB = (v1.x-v0.x)*eAB.x + (v1.y-v0.y)*eAB.y + (v1.z-v0.z)*eAB.z;
        mixed vBC = (v2.x-v1.x)*eBC.x + (v2.y-v1.y)*eBC.y + (v2.z-v1.z)*eBC.z;
        mixed vCA = (v0.x-v2.x)*eCA.x + (v0.y-v2.y)*eCA.y + (v0.z-v2.z)*eCA.z;
        mixed cA = -(eAB.x*eCA.x + eAB.y*eCA.y + eAB.z*eCA.z);
        mixed cB = -(eAB.x*eBC.x + eAB.y*eBC.y + eAB.z*eBC.z);
        mixed cC = -(eBC.x*eCA.x + eBC.y*eCA.y + eBC.z*eCA.z);
        mixed s2A = 1-cA*cA;
        mixed s2B = 1-cB*cB;
        mixed s2C = 1-cC*cC;

        // Solve the equations.  These are different from those in the SETTLE paper (JCC 13(8), pp. 952-962, 1992), because
        // in going from equations B1 to B2, they make the assumption that mB=mC (but don't bother to mention they're
        // making that assumption).  We allow all three atoms to have different masses.

        mixed mABCinv = 1/(mA*mB*mC);
        mixed denom = (((s2A*mB+s2B*mA)*mC+(s2A*mB*mB+2*(cA*cB*cC+1)*mA*mB+s2B*mA*mA))*mC+s2C*mA*mB*(mA+mB))*mABCinv;
        mixed tab = ((cB*cC*mA-cA*mB-cA*mC)*vCA + (cA*cC*mB-cB*mC-cB*mA)*vBC + (s2C*mA*mA*mB*mB*mABCinv+(mA+mB+mC))*vAB)/denom;
        mixed tbc = ((cA*cB*mC-cC*mB-cC*mA)*vCA + (s2A*mB*mB*mC*mC*mABCinv+(mA+mB+mC))*vBC + (cA*cC*mB-cB*mA-cB*mC)*vAB)/denom;
        mixed tca = ((s2B*mA*mA*mC*mC*mABCinv+(mA+mB+mC))*vCA + (cA*cB*mC-cC*mB-cC*mA)*vBC + (cB*cC*mA-cA*mB-cA*mC)*vAB)/denom;
        v0.x += (tab*eAB.x - tca*eCA.x)*v0.w;
        v0.y += (tab*eAB.y - tca*eCA.y)*v0.w;
        v0.z += (tab*eAB.z - tca*eCA.z)*v0.w;
        v1.x += (tbc*eBC.x - tab*eAB.x)*v1.w;
        v1.y += (tbc*eBC.y - tab*eAB.y)*v1.w;
        v1.z += (tbc*eBC.z - tab*eAB.z)*v1.w;
        v2.x += (tca*eCA.x - tbc*eBC.x)*v2.w;
        v2.y += (tca*eCA.y - tbc*eBC.y)*v2.w;
        v2.z += (tca*eCA.z - tbc*eBC.z)*v2.w;
        velm[atoms.x] = v0;
        velm[atoms.y] = v1;
        velm[atoms.z] = v2;
    }
    // Independent grid-stride coverage of the already validated complement.
    // Keep the original Part1 mass-zero branch and arithmetic verbatim.
    for (int slot = GLOBAL_ID; slot < numResidualAtoms; slot += GLOBAL_SIZE) {
        int index = residualAtoms[slot];
        mixed4 velocity = velm[index];
        if (velocity.w != 0.0) {
            velocity.x += fscale*velocity.w*force[index];
            velocity.y += fscale*velocity.w*force[index+paddedNumAtoms];
            velocity.z += fscale*velocity.w*force[index+paddedNumAtoms*2];
            velm[index] = velocity;
        }
    }
}
