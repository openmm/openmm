#define WARPS_PER_GROUP (THREAD_BLOCK_SIZE/TILE_SIZE)

typedef struct {
    real3 pos, force, torque, inducedDipole, inducedDipolePolar, sphericalDipole, dipole;
    real q;
    float penalpha, thole, damp, dirdamp;
    int atomic;
    real quadrupoleXX, quadrupoleXY, quadrupoleXZ;
    real quadrupoleYY, quadrupoleYZ;
#ifdef INCLUDE_QUADRUPOLES
    real sphericalQuadrupole[5];
#endif
} AtomData;

inline __device__ void loadAtomData(AtomData& data, int atom, const real4* __restrict__ posq, const real* __restrict__ sphericalDipole,
        const real* __restrict__ sphericalQuadrupole, const real* __restrict__ inducedDipole, const real* __restrict__ inducedDipolePolar,
        const float2* __restrict__ dampingAndThole, const real* __restrict__ dirDamping,
        const real* __restrict__ labFrameDipole, const real* __restrict__ labFrameQuadrupole, 
        const real* __restrict__ penAlpha, const int* __restrict__ atomicNum) {
    real4 atomPosq = posq[atom];
    data.pos = make_real3(atomPosq.x, atomPosq.y, atomPosq.z);
    data.q = atomPosq.w;
    data.sphericalDipole.x = sphericalDipole[atom*3];
    data.sphericalDipole.y = sphericalDipole[atom*3+1];
    data.sphericalDipole.z = sphericalDipole[atom*3+2];
#ifdef INCLUDE_QUADRUPOLES
    data.sphericalQuadrupole[0] = sphericalQuadrupole[atom*5];
    data.sphericalQuadrupole[1] = sphericalQuadrupole[atom*5+1];
    data.sphericalQuadrupole[2] = sphericalQuadrupole[atom*5+2];
    data.sphericalQuadrupole[3] = sphericalQuadrupole[atom*5+3];
    data.sphericalQuadrupole[4] = sphericalQuadrupole[atom*5+4];
#endif
    data.inducedDipole.x = inducedDipole[atom*3];
    data.inducedDipole.y = inducedDipole[atom*3+1];
    data.inducedDipole.z = inducedDipole[atom*3+2];
    data.inducedDipolePolar.x = inducedDipolePolar[atom*3];
    data.inducedDipolePolar.y = inducedDipolePolar[atom*3+1];
    data.inducedDipolePolar.z = inducedDipolePolar[atom*3+2];
    float2 temp = dampingAndThole[atom];
    data.damp = temp.x;
    data.thole = temp.y;
    data.dirdamp = dirDamping[atom];
    //Charge penetration used
    data.dipole.x = labFrameDipole[atom*3];
    data.dipole.y = labFrameDipole[atom*3+1];
    data.dipole.z = labFrameDipole[atom*3+2];
    data.quadrupoleXX = labFrameQuadrupole[atom*5];
    data.quadrupoleXY = labFrameQuadrupole[atom*5+1];
    data.quadrupoleXZ = labFrameQuadrupole[atom*5+2];
    data.quadrupoleYY = labFrameQuadrupole[atom*5+3];
    data.quadrupoleYZ = labFrameQuadrupole[atom*5+4];
    data.penalpha = penAlpha[atom]; 
    data.atomic = atomicNum[atom]; 
}



__device__ real computeDScaleFactor(unsigned int polarizationGroup, int index) {
    return (polarizationGroup & 1<<index ? 0 : 1);
}

__device__ float computeMScaleFactor(uint2 covalent, int index) {
    int mask = 1<<index;
    bool x = (covalent.x & mask);
    bool y = (covalent.y & mask);
    return (x ? (y ? 0.0f : 0.4f) : (y ? 0.8f : 1.0f));
}

__device__ float computePScaleFactor(uint2 covalent, unsigned int polarizationGroup, int index) {
    int mask = 1<<index;
    bool x = (covalent.x & mask);
    bool y = (covalent.y & mask);
    bool p = (polarizationGroup & mask);
    return (x && y ? 0.0f : (x && p ? 0.5f : 1.0f));
}

__device__ void computeOneInteraction(AtomData& atom1, AtomData& atom2, bool hasExclusions, float dScale, float pScale, float mScale, float forceFactor, mixed& energy) {
    // Compute the displacement.
    
    real3 delta;
    delta.x = atom2.pos.x - atom1.pos.x;
    delta.y = atom2.pos.y - atom1.pos.y;
    delta.z = atom2.pos.z - atom1.pos.z;
    real r2 = delta.x*delta.x + delta.y*delta.y + delta.z*delta.z;
    real rInv = RSQRT(r2);
    real r = r2*rInv;

    // Rotate the various dipoles and quadrupoles.

    real qiRotationMatrix[3][3];
    buildQIRotationMatrix(delta, rInv, qiRotationMatrix);

    
    real3 qiUindI = 0.5f*make_real3(qiRotationMatrix[0][1]*atom1.inducedDipole.x + qiRotationMatrix[0][2]*atom1.inducedDipole.y + qiRotationMatrix[0][0]*atom1.inducedDipole.z,
                                    qiRotationMatrix[1][1]*atom1.inducedDipole.x + qiRotationMatrix[1][2]*atom1.inducedDipole.y + qiRotationMatrix[1][0]*atom1.inducedDipole.z,
                                    qiRotationMatrix[2][1]*atom1.inducedDipole.x + qiRotationMatrix[2][2]*atom1.inducedDipole.y + qiRotationMatrix[2][0]*atom1.inducedDipole.z);
    real3 qiUindJ = 0.5f*make_real3(qiRotationMatrix[0][1]*atom2.inducedDipole.x + qiRotationMatrix[0][2]*atom2.inducedDipole.y + qiRotationMatrix[0][0]*atom2.inducedDipole.z,
                                    qiRotationMatrix[1][1]*atom2.inducedDipole.x + qiRotationMatrix[1][2]*atom2.inducedDipole.y + qiRotationMatrix[1][0]*atom2.inducedDipole.z,
                                    qiRotationMatrix[2][1]*atom2.inducedDipole.x + qiRotationMatrix[2][2]*atom2.inducedDipole.y + qiRotationMatrix[2][0]*atom2.inducedDipole.z);
    real3 qiUinpI = 0.5f*make_real3(qiRotationMatrix[0][1]*atom1.inducedDipolePolar.x + qiRotationMatrix[0][2]*atom1.inducedDipolePolar.y + qiRotationMatrix[0][0]*atom1.inducedDipolePolar.z,
                                    qiRotationMatrix[1][1]*atom1.inducedDipolePolar.x + qiRotationMatrix[1][2]*atom1.inducedDipolePolar.y + qiRotationMatrix[1][0]*atom1.inducedDipolePolar.z,
                                    qiRotationMatrix[2][1]*atom1.inducedDipolePolar.x + qiRotationMatrix[2][2]*atom1.inducedDipolePolar.y + qiRotationMatrix[2][0]*atom1.inducedDipolePolar.z);
    real3 qiUinpJ = 0.5f*make_real3(qiRotationMatrix[0][1]*atom2.inducedDipolePolar.x + qiRotationMatrix[0][2]*atom2.inducedDipolePolar.y + qiRotationMatrix[0][0]*atom2.inducedDipolePolar.z,
                                    qiRotationMatrix[1][1]*atom2.inducedDipolePolar.x + qiRotationMatrix[1][2]*atom2.inducedDipolePolar.y + qiRotationMatrix[1][0]*atom2.inducedDipolePolar.z,
                                    qiRotationMatrix[2][1]*atom2.inducedDipolePolar.x + qiRotationMatrix[2][2]*atom2.inducedDipolePolar.y + qiRotationMatrix[2][0]*atom2.inducedDipolePolar.z);
    
    real3 rotatedDipole1 = rotateDipole(atom1.sphericalDipole, qiRotationMatrix);
    real3 rotatedDipole2 = rotateDipole(atom2.sphericalDipole, qiRotationMatrix);
    real rotatedQuadrupole1[] = {0, 0, 0, 0, 0};
    real rotatedQuadrupole2[] = {0, 0, 0, 0, 0};
#ifdef INCLUDE_QUADRUPOLES
    rotateQuadupoles(qiRotationMatrix, atom1.sphericalQuadrupole, atom2.sphericalQuadrupole, rotatedQuadrupole1, rotatedQuadrupole2);
#endif    
    
    // The field derivatives at I due to permanent and induced moments on J, and vice-versa.
    // Also, their derivatives w.r.t. R, which are needed for force calculations
    real Vij[9], Vji[9], VjiR[9], VijR[9];
    // The field derivatives at I due to only permanent moments on J, and vice-versa.
    real Vijp[3], Vijd[3], Vjip[3], Vjid[3];
    real rInvVec[7];

    // The rInvVec array is defined such that the ith element is R^-i, with the
    // dieleectric constant folded in, to avoid conversions later.
    rInvVec[1] = rInv;
    for (int i = 2; i < 7; ++i)
        rInvVec[i] = rInvVec[i-1] * rInv;

    real dmp = atom1.damp*atom2.damp;
    real a = min(atom1.dirdamp, atom2.dirdamp);
    real u = SQRT(r/dmp);
        
        // Here we use the newly derived damping function
        // [REF] Chengwen Liu, Rui Qi, Qiantao Wang, J-P. Piquemal, Pengyu Ren
        // Capturing Many-Body Interactions with Classical Dipole Induction Models
        // J.Chem.TheoryComput.13,2751,(2017)

    real au3 = fabs(dmp) > 1.0e-5f ? a*u*u*u : 0;
    real expau3 = fabs(dmp) > 1.0e-5f ? EXP(-au3) : 0;
    real a2u6 = au3*au3;
    real a3u9 = a2u6*au3;
    // MBfunc damping factors for energies
    real mbfunc_c  = 1 - expau3;                              
    real mbfunc_d0 = 1 - expau3*(1 + 0.75f*au3);              
    real mbfunc_d1 = 1 - expau3;                              
    real mbfunc_q0 = 1 - expau3*(1 + 0.75f*au3 + 0.25f*a2u6); 
    real mbfunc_q1 = 1 - expau3*(1 + 0.5f*au3);               
    // MBfunc damping factors for derivatives
    real dmbfunc_c  = 1 - expau3*(1 + 0.75f*au3);
    real dmbfunc_d0 = 1 - expau3*(1 + 0.875f*au3 + 0.375f*a2u6);
    real dmbfunc_d1 = 1 - expau3*(1 + 0.5f*au3);
    real dmbfunc_q0 = 1 - expau3*(1 + 0.84375*au3 + 0.34375f*a2u6 + 0.09375f*a3u9);
    real dmbfunc_q1 = 1 - expau3*(1 + 0.6875f*au3 + 0.1875f*a2u6);

    // Now we compute the (attenuated) Coulomb operator and its derivatives, contracted with
    // permanent moments and induced dipoles.  Note that the coefficient of the permanent force
    // terms is half of the expected value; this is because we compute the interaction of I with
    // the sum of induced and permanent moments on J, as well as the interaction of J with I's
    // permanent and induced moments; doing so double counts the permanent-permanent interaction.
    real ePermCoef, dPermCoef, eUIndCoef, dUIndCoef, eUInpCoef, dUInpCoef;

    // C-C terms (m=0)
    ePermCoef = rInvVec[1]*mScale;
    dPermCoef = -0.5f*mScale*rInvVec[2];
    Vij[0]  = ePermCoef*atom2.q;
    Vji[0]  = ePermCoef*atom1.q;
    VijR[0] = dPermCoef*atom2.q;
    VjiR[0] = dPermCoef*atom1.q;

    // C-D and C-Uind terms (m=0)
    ePermCoef = rInvVec[2]*mScale;
    eUIndCoef = rInvVec[2]*pScale*mbfunc_c;
    eUInpCoef = rInvVec[2]*dScale*mbfunc_c;
    dPermCoef = -rInvVec[3]*mScale;
    dUIndCoef = -2*rInvVec[3]*pScale*dmbfunc_c;
    dUInpCoef = -2*rInvVec[3]*dScale*dmbfunc_c;
    Vij[0]  += -(ePermCoef*rotatedDipole2.x + eUIndCoef*qiUindJ.x + eUInpCoef*qiUinpJ.x);
    Vji[1]   = -(ePermCoef*atom1.q);
    VijR[0] += -(dPermCoef*rotatedDipole2.x + dUIndCoef*qiUindJ.x + dUInpCoef*qiUinpJ.x);
    VjiR[1]  = -(dPermCoef*atom1.q);
    Vjip[0]  = -(eUInpCoef*atom1.q);
    Vjid[0]  = -(eUIndCoef*atom1.q);
    // D-C and Uind-C terms (m=0)
    Vij[1]   = ePermCoef*atom2.q;
    Vji[0]  += ePermCoef*rotatedDipole1.x + eUIndCoef*qiUindI.x + eUInpCoef*qiUinpI.x;
    VijR[1]  = dPermCoef*atom2.q;
    VjiR[0] += dPermCoef*rotatedDipole1.x + dUIndCoef*qiUindI.x + dUInpCoef*qiUinpI.x;
    Vijp[0]  = eUInpCoef*atom2.q;
    Vijd[0]  = eUIndCoef*atom2.q;

    // D-D and D-Uind terms (m=0)
    ePermCoef = -2*rInvVec[3]*mScale;
    eUIndCoef = -2*rInvVec[3]*pScale*mbfunc_d0;
    eUInpCoef = -2*rInvVec[3]*dScale*mbfunc_d0;
    dPermCoef = 3*rInvVec[4]*mScale;
    dUIndCoef = 6*rInvVec[4]*pScale*dmbfunc_d0;
    dUInpCoef = 6*rInvVec[4]*dScale*dmbfunc_d0;
    Vij[1]  += ePermCoef*rotatedDipole2.x + eUIndCoef*qiUindJ.x + eUInpCoef*qiUinpJ.x;
    Vji[1]  += ePermCoef*rotatedDipole1.x + eUIndCoef*qiUindI.x + eUInpCoef*qiUinpI.x;
    VijR[1] += dPermCoef*rotatedDipole2.x + dUIndCoef*qiUindJ.x + dUInpCoef*qiUinpJ.x;
    VjiR[1] += dPermCoef*rotatedDipole1.x + dUIndCoef*qiUindI.x + dUInpCoef*qiUinpI.x;
    Vijp[0] += eUInpCoef*rotatedDipole2.x;
    Vijd[0] += eUIndCoef*rotatedDipole2.x;
    Vjip[0] += eUInpCoef*rotatedDipole1.x;
    Vjid[0] += eUIndCoef*rotatedDipole1.x;
    // D-D and D-Uind terms (m=1)
    ePermCoef = rInvVec[3]*mScale;
    eUIndCoef = rInvVec[3]*pScale*mbfunc_d1;
    eUInpCoef = rInvVec[3]*dScale*mbfunc_d1;
    dPermCoef = -1.5f*rInvVec[4]*mScale;
    dUIndCoef = -3*rInvVec[4]*pScale*dmbfunc_d1;
    dUInpCoef = -3*rInvVec[4]*dScale*dmbfunc_d1;
    Vij[2]  = ePermCoef*rotatedDipole2.y + eUIndCoef*qiUindJ.y + eUInpCoef*qiUinpJ.y;
    Vji[2]  = ePermCoef*rotatedDipole1.y + eUIndCoef*qiUindI.y + eUInpCoef*qiUinpI.y;
    VijR[2] = dPermCoef*rotatedDipole2.y + dUIndCoef*qiUindJ.y + dUInpCoef*qiUinpJ.y;
    VjiR[2] = dPermCoef*rotatedDipole1.y + dUIndCoef*qiUindI.y + dUInpCoef*qiUinpI.y;
    Vij[3]  = ePermCoef*rotatedDipole2.z + eUIndCoef*qiUindJ.z + eUInpCoef*qiUinpJ.z;
    Vji[3]  = ePermCoef*rotatedDipole1.z + eUIndCoef*qiUindI.z + eUInpCoef*qiUinpI.z;
    VijR[3] = dPermCoef*rotatedDipole2.z + dUIndCoef*qiUindJ.z + dUInpCoef*qiUinpJ.z;
    VjiR[3] = dPermCoef*rotatedDipole1.z + dUIndCoef*qiUindI.z + dUInpCoef*qiUinpI.z;
    Vijp[1] = eUInpCoef*rotatedDipole2.y;
    Vijd[1] = eUIndCoef*rotatedDipole2.y;
    Vjip[1] = eUInpCoef*rotatedDipole1.y;
    Vjid[1] = eUIndCoef*rotatedDipole1.y;
    Vijp[2] = eUInpCoef*rotatedDipole2.z;
    Vijd[2] = eUIndCoef*rotatedDipole2.z;
    Vjip[2] = eUInpCoef*rotatedDipole1.z;
    Vjid[2] = eUIndCoef*rotatedDipole1.z;

    // C-Q terms (m=0)
    ePermCoef = mScale*rInvVec[3];
    dPermCoef = -1.5f*rInvVec[4]*mScale;
    Vij[0]  += ePermCoef*rotatedQuadrupole2[0];
    Vji[4]   = ePermCoef*atom1.q;
    VijR[0] += dPermCoef*rotatedQuadrupole2[0];
    VjiR[4]  = dPermCoef*atom1.q;
    // Q-C terms (m=0)
    Vij[4]   = ePermCoef*atom2.q;
    Vji[0]  += ePermCoef*rotatedQuadrupole1[0];
    VijR[4]  = dPermCoef*atom2.q;
    VjiR[0] += dPermCoef*rotatedQuadrupole1[0];

    // D-Q and Uind-Q terms (m=0)
    ePermCoef = rInvVec[4]*3.0*mScale;
    eUIndCoef = rInvVec[4]*3.0*pScale*mbfunc_q0;
    eUInpCoef = rInvVec[4]*3.0*dScale*mbfunc_q0;
    dPermCoef = -6*rInvVec[5]*mScale;
    dUIndCoef = -12*rInvVec[5]*pScale*dmbfunc_q0;
    dUInpCoef = -12*rInvVec[5]*dScale*dmbfunc_q0;
    Vij[1]  += ePermCoef*rotatedQuadrupole2[0];
    Vji[4]  += ePermCoef*rotatedDipole1.x + eUIndCoef*qiUindI.x + eUInpCoef*qiUinpI.x;
    VijR[1] += dPermCoef*rotatedQuadrupole2[0];
    VjiR[4] += dPermCoef*rotatedDipole1.x + dUIndCoef*qiUindI.x + dUInpCoef*qiUinpI.x;
    Vijp[0] += eUInpCoef*rotatedQuadrupole2[0];
    Vijd[0] += eUIndCoef*rotatedQuadrupole2[0];
    // Q-D and Q-Uind terms (m=0)
    Vij[4]  += -(ePermCoef*rotatedDipole2.x + eUIndCoef*qiUindJ.x + eUInpCoef*qiUinpJ.x);
    Vji[1]  += -(ePermCoef*rotatedQuadrupole1[0]);
    VijR[4] += -(dPermCoef*rotatedDipole2.x + dUIndCoef*qiUindJ.x + dUInpCoef*qiUinpJ.x);
    VjiR[1] += -(dPermCoef*rotatedQuadrupole1[0]);
    Vjip[0] += -(eUInpCoef*rotatedQuadrupole1[0]);
    Vjid[0] += -(eUIndCoef*rotatedQuadrupole1[0]);

    // D-Q and Uind-Q terms (m=1)
    const real sqrtThree = SQRT((real) 3);
    ePermCoef = -sqrtThree*rInvVec[4]*mScale;
    eUIndCoef = -sqrtThree*rInvVec[4]*pScale*mbfunc_q1;
    eUInpCoef = -sqrtThree*rInvVec[4]*dScale*mbfunc_q1;
    dPermCoef = 2*sqrtThree*rInvVec[5]*mScale;
    dUIndCoef = 4*sqrtThree*rInvVec[5]*pScale*dmbfunc_q1;
    dUInpCoef = 4*sqrtThree*rInvVec[5]*dScale*dmbfunc_q1;
    Vij[2]  += ePermCoef*rotatedQuadrupole2[1];
    Vji[5]   = ePermCoef*rotatedDipole1.y + eUIndCoef*qiUindI.y + eUInpCoef*qiUinpI.y;
    VijR[2] += dPermCoef*rotatedQuadrupole2[1];
    VjiR[5]  = dPermCoef*rotatedDipole1.y + dUIndCoef*qiUindI.y + dUInpCoef*qiUinpI.y;
    Vij[3]  += ePermCoef*rotatedQuadrupole2[2];
    Vji[6]   = ePermCoef*rotatedDipole1.z + eUIndCoef*qiUindI.z + eUInpCoef*qiUinpI.z;
    VijR[3] += dPermCoef*rotatedQuadrupole2[2];
    VjiR[6]  = dPermCoef*rotatedDipole1.z + dUIndCoef*qiUindI.z + dUInpCoef*qiUinpI.z;
    Vijp[1] += eUInpCoef*rotatedQuadrupole2[1];
    Vijd[1] += eUIndCoef*rotatedQuadrupole2[1];
    Vijp[2] += eUInpCoef*rotatedQuadrupole2[2];
    Vijd[2] += eUIndCoef*rotatedQuadrupole2[2];
    // D-Q and Uind-Q terms (m=1)
    Vij[5]   = -(ePermCoef*rotatedDipole2.y + eUIndCoef*qiUindJ.y + eUInpCoef*qiUinpJ.y);
    Vji[2]  += -(ePermCoef*rotatedQuadrupole1[1]);
    VijR[5]  = -(dPermCoef*rotatedDipole2.y + dUIndCoef*qiUindJ.y + dUInpCoef*qiUinpJ.y);
    VjiR[2] += -(dPermCoef*rotatedQuadrupole1[1]);
    Vij[6]   = -(ePermCoef*rotatedDipole2.z + eUIndCoef*qiUindJ.z + eUInpCoef*qiUinpJ.z);
    Vji[3]  += -(ePermCoef*rotatedQuadrupole1[2]);
    VijR[6]  = -(dPermCoef*rotatedDipole2.z + dUIndCoef*qiUindJ.z + dUInpCoef*qiUinpJ.z);
    VjiR[3] += -(dPermCoef*rotatedQuadrupole1[2]);
    Vjip[1] += -(eUInpCoef*rotatedQuadrupole1[1]);
    Vjid[1] += -(eUIndCoef*rotatedQuadrupole1[1]);
    Vjip[2] += -(eUInpCoef*rotatedQuadrupole1[2]);
    Vjid[2] += -(eUIndCoef*rotatedQuadrupole1[2]);

    // Q-Q terms (m=0)
    ePermCoef = 6*rInvVec[5]*mScale;
    dPermCoef = -15*rInvVec[6]*mScale;
    Vij[4]  += ePermCoef*rotatedQuadrupole2[0];
    Vji[4]  += ePermCoef*rotatedQuadrupole1[0];
    VijR[4] += dPermCoef*rotatedQuadrupole2[0];
    VjiR[4] += dPermCoef*rotatedQuadrupole1[0];
    // Q-Q terms (m=1)
    ePermCoef = -4*rInvVec[5]*mScale;
    dPermCoef = 10*rInvVec[6]*mScale;
    Vij[5]  += ePermCoef*rotatedQuadrupole2[1];
    Vji[5]  += ePermCoef*rotatedQuadrupole1[1];
    VijR[5] += dPermCoef*rotatedQuadrupole2[1];
    VjiR[5] += dPermCoef*rotatedQuadrupole1[1];
    Vij[6]  += ePermCoef*rotatedQuadrupole2[2];
    Vji[6]  += ePermCoef*rotatedQuadrupole1[2];
    VijR[6] += dPermCoef*rotatedQuadrupole2[2];
    VjiR[6] += dPermCoef*rotatedQuadrupole1[2];
    // Q-Q terms (m=2)
    ePermCoef = rInvVec[5]*mScale;
    dPermCoef = -2.5f*rInvVec[6]*mScale;
    Vij[7]  = ePermCoef*rotatedQuadrupole2[3];
    Vji[7]  = ePermCoef*rotatedQuadrupole1[3];
    VijR[7] = dPermCoef*rotatedQuadrupole2[3];
    VjiR[7] = dPermCoef*rotatedQuadrupole1[3];
    Vij[8]  = ePermCoef*rotatedQuadrupole2[4];
    Vji[8]  = ePermCoef*rotatedQuadrupole1[4];
    VijR[8] = dPermCoef*rotatedQuadrupole2[4];
    VjiR[8] = dPermCoef*rotatedQuadrupole1[4];

    // Evaluate the energies, forces and torques due to permanent+induced moments
    // interacting with just the permanent moments.
    energy += forceFactor*0.5f*(
        atom1.q*Vij[0] + rotatedDipole1.x*Vij[1] + rotatedDipole1.y*Vij[2] + rotatedDipole1.z*Vij[3] + rotatedQuadrupole1[0]*Vij[4] + rotatedQuadrupole1[1]*Vij[5] + rotatedQuadrupole1[2]*Vij[6] + rotatedQuadrupole1[3]*Vij[7] + rotatedQuadrupole1[4]*Vij[8] +
        atom2.q*Vji[0] + rotatedDipole2.x*Vji[1] + rotatedDipole2.y*Vji[2] + rotatedDipole2.z*Vji[3] + rotatedQuadrupole2[0]*Vji[4] + rotatedQuadrupole2[1]*Vji[5] + rotatedQuadrupole2[2]*Vji[6] + rotatedQuadrupole2[3]*Vji[7] + rotatedQuadrupole2[4]*Vji[8]);
    real fIZ = atom1.q*VijR[0] + rotatedDipole1.x*VijR[1] + rotatedDipole1.y*VijR[2] + rotatedDipole1.z*VijR[3] + rotatedQuadrupole1[0]*VijR[4] + rotatedQuadrupole1[1]*VijR[5] + rotatedQuadrupole1[2]*VijR[6] + rotatedQuadrupole1[3]*VijR[7] + rotatedQuadrupole1[4]*VijR[8];
    real fJZ = atom2.q*VjiR[0] + rotatedDipole2.x*VjiR[1] + rotatedDipole2.y*VjiR[2] + rotatedDipole2.z*VjiR[3] + rotatedQuadrupole2[0]*VjiR[4] + rotatedQuadrupole2[1]*VjiR[5] + rotatedQuadrupole2[2]*VjiR[6] + rotatedQuadrupole2[3]*VjiR[7] + rotatedQuadrupole2[4]*VjiR[8];
    real EIX = rotatedDipole1.z*Vij[1] - rotatedDipole1.x*Vij[3] + sqrtThree*rotatedQuadrupole1[2]*Vij[4] + rotatedQuadrupole1[4]*Vij[5] - (sqrtThree*rotatedQuadrupole1[0]+rotatedQuadrupole1[3])*Vij[6] + rotatedQuadrupole1[2]*Vij[7] - rotatedQuadrupole1[1]*Vij[8];
    real EIY = -rotatedDipole1.y*Vij[1] + rotatedDipole1.x*Vij[2] - sqrtThree*rotatedQuadrupole1[1]*Vij[4] + (sqrtThree*rotatedQuadrupole1[0]-rotatedQuadrupole1[3])*Vij[5] - rotatedQuadrupole1[4]*Vij[6] + rotatedQuadrupole1[1]*Vij[7] + rotatedQuadrupole1[2]*Vij[8];
    real EIZ = -rotatedDipole1.z*Vij[2] + rotatedDipole1.y*Vij[3] - rotatedQuadrupole1[2]*Vij[5] + rotatedQuadrupole1[1]*Vij[6] - 2*rotatedQuadrupole1[4]*Vij[7] + 2*rotatedQuadrupole1[3]*Vij[8];
    real EJX = rotatedDipole2.z*Vji[1] - rotatedDipole2.x*Vji[3] + sqrtThree*rotatedQuadrupole2[2]*Vji[4] + rotatedQuadrupole2[4]*Vji[5] - (sqrtThree*rotatedQuadrupole2[0]+rotatedQuadrupole2[3])*Vji[6] + rotatedQuadrupole2[2]*Vji[7] - rotatedQuadrupole2[1]*Vji[8];
    real EJY = -rotatedDipole2.y*Vji[1] + rotatedDipole2.x*Vji[2] - sqrtThree*rotatedQuadrupole2[1]*Vji[4] + (sqrtThree*rotatedQuadrupole2[0]-rotatedQuadrupole2[3])*Vji[5] - rotatedQuadrupole2[4]*Vji[6] + rotatedQuadrupole2[1]*Vji[7] + rotatedQuadrupole2[2]*Vji[8];
    real EJZ = -rotatedDipole2.z*Vji[2] + rotatedDipole2.y*Vji[3] - rotatedQuadrupole2[2]*Vji[5] + rotatedQuadrupole2[1]*Vji[6] - 2*rotatedQuadrupole2[4]*Vji[7] + 2*rotatedQuadrupole2[3]*Vji[8];

    // Define the torque intermediates for the induced dipoles. These are simply the induced dipole torque
    // intermediates dotted with the field due to permanent moments only, at each center. We inline the
    // induced dipole torque intermediates here, for simplicity. N.B. There are no torques on the dipoles
    // themselves, so we accumulate the torque intermediates into separate variables to allow them to be
    // used only in the force calculation.
    //
    // The torque about the x axis (needed to obtain the y force on the induced dipoles, below)
    //    qiUindIx[0] = qiQUindI[2];    qiUindIx[1] = 0;    qiUindIx[2] = -qiQUindI[0]
    real iEIX = qiUinpI.z*Vijp[0] + qiUindI.z*Vijd[0] - qiUinpI.x*Vijp[2] - qiUindI.x*Vijd[2];
    real iEJX = qiUinpJ.z*Vjip[0] + qiUindJ.z*Vjid[0] - qiUinpJ.x*Vjip[2] - qiUindJ.x*Vjid[2];
    // The torque about the y axis (needed to obtain the x force on the induced dipoles, below)
    //    qiUindIy[0] = -qiQUindI[1];   qiUindIy[1] = qiQUindI[0];    qiUindIy[2] = 0
    real iEIY = qiUinpI.x*Vijp[1] + qiUindI.x*Vijd[1] - qiUinpI.y*Vijp[0] - qiUindI.y*Vijd[0];
    real iEJY = qiUinpJ.x*Vjip[1] + qiUindJ.x*Vjid[1] - qiUinpJ.y*Vjip[0] - qiUindJ.y*Vjid[0];

#ifdef MUTUAL_POLARIZATION

    //Re-calculate thole damping factors for MUTUAL polarization
    
    real dmpt = atom1.damp*atom2.damp;
    real at = min(atom1.thole, atom2.thole);
    real ut = r/dmpt;
    real atut3 = fabs(dmpt) > 1.0e-5f ? at*ut*ut*ut : 0;
    real expatut3 = fabs(dmpt) > 1.0e-5f ? EXP(-atut3) : 0;
    real at2ut6 = atut3*atut3;
    // Thole damping factors for energies
    real thole_d0 = 1 - expatut3*(1 + 1.5f*atut3);              
    real thole_d1 = 1 - expatut3;                              
    // Thole damping factors for derivatives
    real dthole_d0 = 1 - expatut3*(1 + atut3 + 1.5f*at2ut6);
    real dthole_d1 = 1 - expatut3*(1 + atut3);
    
    // Uind-Uind terms (m=0)
    real eCoef = -4*rInvVec[3]*thole_d0;
    real dCoef = 6*rInvVec[4]*dthole_d0;
    iEIX += eCoef*(qiUinpI.z*qiUindJ.x + qiUindI.z*qiUinpJ.x);
    iEJX += eCoef*(qiUinpJ.z*qiUindI.x + qiUindJ.z*qiUinpI.x);
    iEIY -= eCoef*(qiUinpI.y*qiUindJ.x + qiUindI.y*qiUinpJ.x);
    iEJY -= eCoef*(qiUinpJ.y*qiUindI.x + qiUindJ.y*qiUinpI.x);
    fIZ += dCoef*(qiUinpI.x*qiUindJ.x + qiUindI.x*qiUinpJ.x);
    fIZ += dCoef*(qiUinpJ.x*qiUindI.x + qiUindJ.x*qiUinpI.x);
    // Uind-Uind terms (m=1)
    eCoef = 2*rInvVec[3]*thole_d1;
    dCoef = -3*rInvVec[4]*dthole_d1;
    iEIX -= eCoef*(qiUinpI.x*qiUindJ.z + qiUindI.x*qiUinpJ.z);
    iEJX -= eCoef*(qiUinpJ.x*qiUindI.z + qiUindJ.x*qiUinpI.z);
    iEIY += eCoef*(qiUinpI.x*qiUindJ.y + qiUindI.x*qiUinpJ.y);
    iEJY += eCoef*(qiUinpJ.x*qiUindI.y + qiUindJ.x*qiUinpI.y);
    fIZ += dCoef*(qiUinpI.y*qiUindJ.y + qiUindI.y*qiUinpJ.y + qiUinpI.z*qiUindJ.z + qiUindI.z*qiUinpJ.z);
    fIZ += dCoef*(qiUinpJ.y*qiUindI.y + qiUindJ.y*qiUinpI.y + qiUinpJ.z*qiUindI.z + qiUindJ.z*qiUinpI.z);
#endif

    // The quasi-internal frame forces and torques.  Note that the induced torque intermediates are
    // used in the force expression, but not in the torques; the induced dipoles are isotropic.
    real qiForce[3] = {rInv*(EIY+EJY+iEIY+iEJY), -rInv*(EIX+EJX+iEIX+iEJX), -(fJZ+fIZ)};
    real qiTorqueI[3] = {-EIX, -EIY, -EIZ};
    real qiTorqueJ[3] = {-EJX, -EJY, -EJZ};


    real3 force = make_real3(qiRotationMatrix[1][1]*qiForce[0] + qiRotationMatrix[2][1]*qiForce[1] + qiRotationMatrix[0][1]*qiForce[2],
                             qiRotationMatrix[1][2]*qiForce[0] + qiRotationMatrix[2][2]*qiForce[1] + qiRotationMatrix[0][2]*qiForce[2],
                             qiRotationMatrix[1][0]*qiForce[0] + qiRotationMatrix[2][0]*qiForce[1] + qiRotationMatrix[0][0]*qiForce[2]);
    atom1.force += force;
    atom1.torque += make_real3(qiRotationMatrix[1][1]*qiTorqueI[0] + qiRotationMatrix[2][1]*qiTorqueI[1] + qiRotationMatrix[0][1]*qiTorqueI[2],
                               qiRotationMatrix[1][2]*qiTorqueI[0] + qiRotationMatrix[2][2]*qiTorqueI[1] + qiRotationMatrix[0][2]*qiTorqueI[2],
                               qiRotationMatrix[1][0]*qiTorqueI[0] + qiRotationMatrix[2][0]*qiTorqueI[1] + qiRotationMatrix[0][0]*qiTorqueI[2]);
    if (forceFactor == 1) {
        atom2.force -= force;
        atom2.torque += make_real3(qiRotationMatrix[1][1]*qiTorqueJ[0] + qiRotationMatrix[2][1]*qiTorqueJ[1] + qiRotationMatrix[0][1]*qiTorqueJ[2],
                                   qiRotationMatrix[1][2]*qiTorqueJ[0] + qiRotationMatrix[2][2]*qiTorqueJ[1] + qiRotationMatrix[0][2]*qiTorqueJ[2],
                                   qiRotationMatrix[1][0]*qiTorqueJ[0] + qiRotationMatrix[2][0]*qiTorqueJ[1] + qiRotationMatrix[0][0]*qiTorqueJ[2]);
    }
}

__device__ void computeOneInteractionCP(AtomData& atom1, AtomData& atom2, bool hasExclusions, float mScale, float forceFactor, mixed& energy) {
    /*Chengwen Liu
      Sep.2018*/
   //Charge Penetration added as a correction 

    // Define some distances 
    real xr = atom2.pos.x - atom1.pos.x;
    real yr = atom2.pos.y - atom1.pos.y;
    real zr = atom2.pos.z - atom1.pos.z;
    real r2 = xr*xr + yr*yr +zr*zr;
    real r  = SQRT(r2);

    // get reciprocal distance terms for this interaction
    real rr1 = mScale*RECIP(r);
    real rr3 = rr1/r2;
    real rr5 = 3*rr3/r2;
    real rr7 = 5*rr5/r2;
    real rr9 = 7*rr7/r2;
    real rr11 = 9*rr9/r2;

    //intemediates involving moments and distance separation
    //dipoles
    real dikx = atom1.dipole.y*atom2.dipole.z - atom1.dipole.z*atom2.dipole.y;
    real diky = atom1.dipole.z*atom2.dipole.x - atom1.dipole.x*atom2.dipole.z;
    real dikz = atom1.dipole.x*atom2.dipole.y - atom1.dipole.y*atom2.dipole.x;
    
    real dirx = atom1.dipole.y*zr - atom1.dipole.z*yr; 
    real diry = atom1.dipole.z*xr - atom1.dipole.x*zr;
    real dirz = atom1.dipole.x*yr - atom1.dipole.y*xr;

    real dkrx = atom2.dipole.y*zr - atom2.dipole.z*yr; 
    real dkry = atom2.dipole.z*xr - atom2.dipole.x*zr;
    real dkrz = atom2.dipole.x*yr - atom2.dipole.y*xr;

    real dri  = atom1.dipole.x*xr + atom1.dipole.y*yr + atom1.dipole.z*zr;
    real drk  = atom2.dipole.x*xr + atom2.dipole.y*yr + atom2.dipole.z*zr;
    real dik  = atom1.dipole.x*atom2.dipole.x + atom1.dipole.y*atom2.dipole.y + atom1.dipole.z*atom2.dipole.z;
    //quadrupoles
    //atomI <==> atom1
    //qrix = qixx*xr + qixy*yr + qixz*zr
    //qriy = qixy*xr + qiyy*yr + qiyz*zr
    //qriz = qixz*xr + qiyz*yr + qizz*zr
    real atom1quadrupoleZZ = -(atom1.quadrupoleXX+atom1.quadrupoleYY);
    real qrix = atom1.quadrupoleXX*xr + atom1.quadrupoleXY*yr + atom1.quadrupoleXZ*zr;
    real qriy = atom1.quadrupoleXY*xr + atom1.quadrupoleYY*yr + atom1.quadrupoleYZ*zr;
    real qriz = atom1.quadrupoleXZ*xr + atom1.quadrupoleYZ*yr + atom1quadrupoleZZ*zr;

    //quadrupoles
    //atomK <==> atom2
    //qrkx = qkxx*xr + qkxy*yr + qkxz*zr
    //qrky = qkxy*xr + qkyy*yr + qkyz*zr
    //qrkz = qkxz*xr + qkyz*yr + qkzz*zr
    real atom2quadrupoleZZ = -(atom2.quadrupoleXX+atom2.quadrupoleYY);
    real qrkx = atom2.quadrupoleXX*xr + atom2.quadrupoleXY*yr + atom2.quadrupoleXZ*zr;
    real qrky = atom2.quadrupoleXY*xr + atom2.quadrupoleYY*yr + atom2.quadrupoleYZ*zr;
    real qrkz = atom2.quadrupoleXZ*xr + atom2.quadrupoleYZ*yr + atom2quadrupoleZZ*zr;

    //qrri = qrix*xr + qriy*yr + qriz*zr
    //qrrk = qrkx*xr + qrky*yr + qrkz*zr
    //qrrik = qrix*qrkx + qriy*qrky + qriz*qrkz
    //qik = 2.0d0*(qixy*qkxy+qixz*qkxz+qiyz*qkyz) + qixx*qkxx + qiyy*qkyy + qizz*qkzz
    real qrri = qrix*xr + qriy*yr + qriz*zr;
    real qrrk = qrkx*xr + qrky*yr + qrkz*zr;
    real qrrik = qrix*qrkx + qriy*qrky + qriz*qrkz;
    real qik = 2.0f*(atom1.quadrupoleXY*atom2.quadrupoleXY 
                   + atom1.quadrupoleXZ*atom2.quadrupoleXZ 
                   + atom1.quadrupoleYZ*atom2.quadrupoleYZ) 
                   + atom1.quadrupoleXX*atom2.quadrupoleXX 
                   + atom1.quadrupoleYY*atom2.quadrupoleYY 
                   + atom1quadrupoleZZ*atom2quadrupoleZZ;

    //qrixr = qriz*yr - qriy*zr
    //qriyr = qrix*zr - qriz*xr
    //qrizr = qriy*xr - qrix*yr
    //qrkxr = qrkz*yr - qrky*zr
    //qrkyr = qrkx*zr - qrkz*xr
    //qrkzr = qrky*xr - qrkx*yr
    //qrrx = qrky*qriz - qrkz*qriy
    //qrry = qrkz*qrix - qrkx*qriz
    //qrrz = qrkx*qriy - qrky*qrix
    real qrixr = qriz*yr - qriy*zr;
    real qriyr = qrix*zr - qriz*xr;
    real qrizr = qriy*xr - qrix*yr;
    real qrkxr = qrkz*yr - qrky*zr;
    real qrkyr = qrkx*zr - qrkz*xr;
    real qrkzr = qrky*xr - qrkx*yr;
    real qrrx = qrky*qriz - qrkz*qriy;
    real qrry = qrkz*qrix - qrkx*qriz;
    real qrrz = qrkx*qriy - qrky*qrix;

    //qikrx = qixx*qrkx + qixy*qrky + qixz*qrkz
    //qikry = qixy*qrkx + qiyy*qrky + qiyz*qrkz
    //qikrz = qixz*qrkx + qiyz*qrky + qizz*qrkz
    //qkirx = qkxx*qrix + qkxy*qriy + qkxz*qriz
    //qkiry = qkxy*qrix + qkyy*qriy + qkyz*qriz
    //qkirz = qkxz*qrix + qkyz*qriy + qkzz*qriz
    real qikrx = atom1.quadrupoleXX*qrkx 
               + atom1.quadrupoleXY*qrky 
               + atom1.quadrupoleXZ*qrkz;
    real qikry = atom1.quadrupoleXY*qrkx 
               + atom1.quadrupoleYY*qrky 
               + atom1.quadrupoleYZ*qrkz;
    real qikrz = atom1.quadrupoleXZ*qrkx 
               + atom1.quadrupoleYZ*qrky 
               + atom1quadrupoleZZ*qrkz;
    real qkirx = atom2.quadrupoleXX*qrix 
               + atom2.quadrupoleXY*qriy 
               + atom2.quadrupoleXZ*qriz;
    real qkiry = atom2.quadrupoleXY*qrix 
               + atom2.quadrupoleYY*qriy 
               + atom2.quadrupoleYZ*qriz;
    real qkirz = atom2.quadrupoleXZ*qrix 
               + atom2.quadrupoleYZ*qriy 
               + atom2quadrupoleZZ*qriz;

    //qikrxr = qikrz*yr - qikry*zr
    //qikryr = qikrx*zr - qikrz*xr
    //qikrzr = qikry*xr - qikrx*yr
    //qkirxr = qkirz*yr - qkiry*zr
    //qkiryr = qkirx*zr - qkirz*xr
    //qkirzr = qkiry*xr - qkirx*yr
    real qikrxr = qikrz*yr - qikry*zr;
    real qikryr = qikrx*zr - qikrz*xr;
    real qikrzr = qikry*xr - qikrx*yr;
    real qkirxr = qkirz*yr - qkiry*zr;
    real qkiryr = qkirx*zr - qkirz*xr;
    real qkirzr = qkiry*xr - qkirx*yr;

    //diqkx = dix*qkxx + diy*qkxy + diz*qkxz
    //diqky = dix*qkxy + diy*qkyy + diz*qkyz
    //diqkz = dix*qkxz + diy*qkyz + diz*qkzz
    //dkqix = dkx*qixx + dky*qixy + dkz*qixz
    //dkqiy = dkx*qixy + dky*qiyy + dkz*qiyz
    //dkqiz = dkx*qixz + dky*qiyz + dkz*qizz
    //diqrk = dix*qrkx + diy*qrky + diz*qrkz
    //dkqri = dkx*qrix + dky*qriy + dkz*qriz
    real diqkx = atom1.dipole.x*atom2.quadrupoleXX 
               + atom1.dipole.y*atom2.quadrupoleXY 
               + atom1.dipole.z*atom2.quadrupoleXZ; 
    real diqky = atom1.dipole.x*atom2.quadrupoleXY 
               + atom1.dipole.y*atom2.quadrupoleYY 
               + atom1.dipole.z*atom2.quadrupoleYZ;
    real diqkz = atom1.dipole.x*atom2.quadrupoleXZ 
               + atom1.dipole.y*atom2.quadrupoleYZ 
               + atom1.dipole.z*atom2quadrupoleZZ;
    real dkqix = atom2.dipole.x*atom1.quadrupoleXX 
               + atom2.dipole.y*atom1.quadrupoleXY 
               + atom2.dipole.z*atom1.quadrupoleXZ;
    real dkqiy = atom2.dipole.x*atom1.quadrupoleXY 
               + atom2.dipole.y*atom1.quadrupoleYY 
               + atom2.dipole.z*atom1.quadrupoleYZ;
    real dkqiz = atom2.dipole.x*atom1.quadrupoleXZ 
               + atom2.dipole.y*atom1.quadrupoleYZ 
               + atom2.dipole.z*atom1quadrupoleZZ;
    real diqrk = atom1.dipole.x*qrkx 
               + atom1.dipole.y*qrky 
               + atom1.dipole.z*qrkz;
    real dkqri = atom2.dipole.x*qrix 
               + atom2.dipole.y*qriy 
               + atom2.dipole.z*qriz;

    //diqkxr = diqkz*yr - diqky*zr
    //diqkyr = diqkx*zr - diqkz*xr
    //diqkzr = diqky*xr - diqkx*yr
    //dkqixr = dkqiz*yr - dkqiy*zr
    //dkqiyr = dkqix*zr - dkqiz*xr
    //dkqizr = dkqiy*xr - dkqix*yr
    real diqkxr = diqkz*yr - diqky*zr; 
    real diqkyr = diqkx*zr - diqkz*xr;
    real diqkzr = diqky*xr - diqkx*yr;
    real dkqixr = dkqiz*yr - dkqiy*zr;
    real dkqiyr = dkqix*zr - dkqiz*xr;
    real dkqizr = dkqiy*xr - dkqix*yr;

    real dqiqkx = atom1.dipole.y*qrkz - atom1.dipole.z*qrky 
                + atom2.dipole.y*qriz - atom2.dipole.z*qriy
                - 2.0f*(atom1.quadrupoleXY*atom2.quadrupoleXZ
                +       atom1.quadrupoleYY*atom2.quadrupoleYZ
                +       atom1.quadrupoleYZ*atom2quadrupoleZZ
                -       atom1.quadrupoleXZ*atom2.quadrupoleXY
                -       atom1.quadrupoleYZ*atom2.quadrupoleYY
                -       atom1quadrupoleZZ*atom2.quadrupoleYZ);

    real dqiqky = atom1.dipole.z*qrkx - atom1.dipole.x*qrkz 
                + atom2.dipole.z*qrix - atom2.dipole.x*qriz
                - 2.0f*(atom1.quadrupoleXZ*atom2.quadrupoleXX
                +       atom1.quadrupoleYZ*atom2.quadrupoleXY
                +       atom1quadrupoleZZ*atom2.quadrupoleXZ
                -       atom1.quadrupoleXX*atom2.quadrupoleXZ
                -       atom1.quadrupoleXY*atom2.quadrupoleYZ
                -       atom1.quadrupoleXZ*atom2quadrupoleZZ);
    real dqiqkz = atom1.dipole.x*qrky - atom1.dipole.y*qrkx 
                + atom2.dipole.x*qriy - atom2.dipole.y*qrix
                - 2.0f*(atom1.quadrupoleXX*atom2.quadrupoleXY
                +       atom1.quadrupoleXY*atom2.quadrupoleYY
                +       atom1.quadrupoleXZ*atom2.quadrupoleYZ
                -       atom1.quadrupoleXY*atom2.quadrupoleXX
                -       atom1.quadrupoleYY*atom2.quadrupoleXY
                -       atom1.quadrupoleYZ*atom2.quadrupoleXZ);

  //charge penetration damping functions
     
    const real oneThirds = (real) 1/3;
    const real oneFifteens = (real) 1/15; 
    const real threeSevens = (real) 3/7;
    const real twoTwentyFirsts = (real) 2/21;
    const real oneOneOFifths = (real) 1/105;
    const real oneNinths  = (real) 1/9;
    const real fourNinths = (real) 4/9; 
    const real oneSixtyThirds  = (real) 1/63;
    const real oneNineFourFifths  = (real) 1/945;

    real alphai = atom1.penalpha; 
    real alphak = atom2.penalpha; 
    //printf("alphai, alphak, mScale %f %f %f \n", alphai, alphak, mScale);//PASSED
    real dampi = alphai*r;
    real dampk = alphak*r;
    real expdampi = EXP(-dampi);
    real expdampk = EXP(-dampk);

    real alphai2 = alphai*alphai;
    real alphak2 = alphak*alphak;

    real dampi2 = dampi*dampi; 
    real dampk2 = dampk*dampk;
    real dampi3 = dampi2*dampi;
    real dampk3 = dampk2*dampk;
    real dampi4 = dampi3*dampi;
    real dampk4 = dampk3*dampk;
    real dampi5 = dampi4*dampi;
    real dampk5 = dampk4*dampk;
    real dampi6 = dampi5*dampi;
    
    
    real scalei_1, scalei_3, scalei_5, scalei_7;
    real scalek_1, scalek_3, scalek_5, scalek_7;
    real termi, termk;
    real scaleik_1, scaleik_3, scaleik_5, scaleik_7, scaleik_9, scaleik_11; 
    const real oneSixths = (real) 1/6;
    const real oneThirtieths = (real) 1/30;
    const real fourOneOFifths = (real) 4/105;
    const real fiveOneTwoSixths = (real) 5/126;
    const real oneTwoOneZeroths = (real) 1/210;
    const real twoThreeOneFifths = (real) 2/315;
    const real oneOneEightNineZeroths = (real) 1/1890;


    //calculate one-site scale factors
    scalei_1 = -expdampi;
    scalek_1 = -expdampk;
    scalei_3 = -expdampi*(1.0f + dampi);
    scalek_3 = -expdampk*(1.0f + dampk);
    scalei_5 = -expdampi*(1.0f + dampi + oneThirds*dampi2);
    scalek_5 = -expdampk*(1.0f + dampk + oneThirds*dampk2);
    scalei_7 = -expdampi*(1.0f + dampi + 0.4f*dampi2 + oneFifteens*dampi3);
    scalek_7 = -expdampk*(1.0f + dampk + 0.4f*dampk2 + oneFifteens*dampk3);

     //calculate two-site scale factors
    if (alphai != alphak) { 
       termi = alphak2/(alphak2 - alphai2);
       termk = alphai2/(alphai2 - alphak2);
       scaleik_1 = -termi*expdampi - termk*expdampk;
       scaleik_3 = -termi*(1.0f + dampi)*expdampi 
                   -termk*(1.0f + dampk)*expdampk;
       scaleik_5 = -termi*(1.0f + dampi + oneThirds*dampi2)*expdampi 
                   -termk*(1.0f + dampk + oneThirds*dampk2)*expdampk;
       scaleik_7 = -termi*(1.0f + dampi + 0.4f*dampi2 + oneFifteens*dampi3)*expdampi 
                   -termk*(1.0f + dampk + 0.4f*dampk2 + oneFifteens*dampk3)*expdampk;
       scaleik_9 = -termi*(1.0f + dampi + threeSevens*dampi2 + twoTwentyFirsts*dampi3 + oneOneOFifths*dampi4)*expdampi 
                   -termk*(1.0f + dampk + threeSevens*dampk2 + twoTwentyFirsts*dampk3 + oneOneOFifths*dampk4)*expdampk;
       scaleik_11 = -termi*expdampi*(1.0f + dampi + fourNinths*dampi2 + oneNinths*dampi3 + oneSixtyThirds*dampi4 + oneNineFourFifths*dampi5) 
                    -termk*expdampk*(1.0f + dampk + fourNinths*dampk2 + oneNinths*dampk3 + oneSixtyThirds*dampk4 + oneNineFourFifths*dampk5);
    }
    else { 
       scaleik_1  = -expdampi*(1.0f + 0.5f*dampi);
       scaleik_3  = -expdampi*(1.0f + dampi + 0.5f*dampi2);
       scaleik_5  = -expdampi*(1.0f + dampi + 0.5f*dampi2 + oneSixths*dampi3);
       scaleik_7  = -expdampi*(1.0f + dampi + 0.5f*dampi2 + oneSixths*dampi3 + oneThirtieths*dampi4);
       scaleik_9  = -expdampi*(1.0f + dampi + 0.5f*dampi2 + oneSixths*dampi3 + fourOneOFifths*dampi4 + oneTwoOneZeroths*dampi5);
       scaleik_11 = -expdampi*(1.0f + dampi + 0.5f*dampi2 + oneSixths*dampi3 + fiveOneTwoSixths*dampi4 + twoThreeOneFifths*dampi5 + oneOneEightNineZeroths*dampi6);
    }
    
    real nuci = (real) atom1.atomic;
    real nuck = (real) atom2.atomic;       
    //printf("nuci, nuck: %f %f \n", nuci, nuck);//PASSED 
    real qi = atom1.q - nuci;
    real qk = atom2.q - nuck;

    real term1 = nuci*qk*scalek_1 + nuck*qi*scalei_1 + qi*qk*scaleik_1; 
    real term2 = nuck*dri*scalei_3 - nuci*drk*scalek_3 + (dik + qk*dri - qi*drk)*scaleik_3; 
    real term3 = nuci*qrrk*scalek_5 + nuck*qrri*scalei_5+ (-dri*drk + 2.0f*(dkqri-diqrk+qik) + qi*qrrk + qk*qrri)*scaleik_5 ;
    real term4 = (dri*qrrk-drk*qrri-4.0f*qrrik)*scaleik_7;
    real term5 = qrri*qrrk*scaleik_9;

    //compute the energy contribution for this interaction
    energy += forceFactor*(term1*rr1 + term2*rr3 + term3*rr5 + term4*rr7 + term5*rr9);

    //calculate intermediate terms for force and torque
    real de = nuci*qk*rr3*scalek_3 + nuck*qi*rr3*scalei_3 + qi*qk*rr3*scaleik_3 + nuck*dri*rr5*scalei_5 - nuci*drk*rr5*scalek_5
            + (dik + qk*dri-qi*drk)*rr5*scaleik_5 + nuci*qrrk*rr7*scalek_7 + nuck*qrri*rr7*scalei_7 +(-dri*drk + 2.0f*(dkqri-diqrk+qik) 
            + qi*qrrk + qk*qrri)*rr7*scaleik_7 + (dri*qrrk-drk*qrri-4.0f*qrrik)*rr9*scaleik_9 + qrri*qrrk*rr11*scaleik_11;

    real dterm1 = -nuck*rr3*scalei_3 - qk*rr3*scaleik_3 + drk*rr5*scaleik_5 - qrrk*rr7*scaleik_7;
    real dterm2 =  nuci*rr3*scalek_3 + qi*rr3*scaleik_3 + dri*rr5*scaleik_5 + qrri*rr7*scaleik_7;
    real dterm3 = 2.0f * rr5 * scaleik_5;
    real dterm4 = 2.0f *(-nuck*rr5*scalei_5 -qk*rr5*scaleik_5 + drk*rr7*scaleik_7 - qrrk*rr9*scaleik_9);
    real dterm5 = 2.0f *(-nuci*rr5*scalek_5 -qi*rr5*scaleik_5 - dri*rr7*scaleik_7 - qrri*rr9*scaleik_9);
    real dterm6 = 4.0f * rr7 *scaleik_7; 

     //compute the force components for this interaction
    real frcx = de*xr + dterm1*atom1.dipole.x + dterm2*atom2.dipole.x + dterm3*(diqkx-dkqix) + dterm4*qrix + dterm5*qrkx + dterm6*(qikrx+qkirx);    
    real frcy = de*yr + dterm1*atom1.dipole.y + dterm2*atom2.dipole.y + dterm3*(diqky-dkqiy) + dterm4*qriy + dterm5*qrky + dterm6*(qikry+qkiry);
    real frcz = de*zr + dterm1*atom1.dipole.z + dterm2*atom2.dipole.z + dterm3*(diqkz-dkqiz) + dterm4*qriz + dterm5*qrkz + dterm6*(qikrz+qkirz);

    //compute the torque components for this interaction

    rr3 = rr3*scaleik_3;
    real ttmi_1 = -rr3*dikx + dterm1*dirx + dterm3*(dqiqkx+dkqixr) - dterm4*qrixr - dterm6*(qikrxr+qrrx); 
    real ttmi_2 = -rr3*diky + dterm1*diry + dterm3*(dqiqky+dkqiyr) - dterm4*qriyr - dterm6*(qikryr+qrry);
    real ttmi_3 = -rr3*dikz + dterm1*dirz + dterm3*(dqiqkz+dkqizr) - dterm4*qrizr - dterm6*(qikrzr+qrrz);
    real ttmk_1 =  rr3*dikx + dterm2*dkrx - dterm3*(dqiqkx+diqkxr) - dterm5*qrkxr - dterm6*(qkirxr-qrrx);
    real ttmk_2 =  rr3*diky + dterm2*dkry - dterm3*(dqiqky+diqkyr) - dterm5*qrkyr - dterm6*(qkiryr-qrry);
    real ttmk_3 =  rr3*dikz + dterm2*dkrz - dterm3*(dqiqkz+diqkzr) - dterm5*qrkzr - dterm6*(qkirzr-qrrz);

    real3 forceCP = make_real3(frcx, frcy, frcz);
    atom1.force += forceCP;
    atom1.torque += make_real3(ttmi_1, ttmi_2, ttmi_3);
    
    if (forceFactor == 1) {
      atom2.force -= forceCP;
      atom2.torque += make_real3(ttmk_1, ttmk_2, ttmk_3);
    }
}

/**
 * Compute electrostatic interactions.
 */

extern "C" __global__ void computeElectrostatics(
        unsigned long long* __restrict__ forceBuffers, unsigned long long* __restrict__ torqueBuffers, mixed* __restrict__ energyBuffer,
        const real4* __restrict__ posq, const uint2* __restrict__ covalentFlags,const unsigned int* __restrict__ polarizationGroupFlags, 
        const ushort2* __restrict__ exclusionTiles, unsigned int startTileIndex, unsigned int numTileIndices,
#ifdef USE_CUTOFF
        const int* __restrict__ tiles, const unsigned int* __restrict__ interactionCount, real4 periodicBoxSize, real4 invPeriodicBoxSize,
        real4 periodicBoxVecX, real4 periodicBoxVecY, real4 periodicBoxVecZ, unsigned int maxTiles, const real4* __restrict__ blockCenter,
        const unsigned int* __restrict__ interactingAtoms,
#endif
        const real* __restrict__ labFrameDipole, const real* __restrict__ sphericalDipole, const real* __restrict__ labFrameQuadrupole, const real* __restrict__ sphericalQuadrupole, const real* __restrict__ inducedDipole,
        const real* __restrict__ inducedDipolePolar, const real* __restrict__ penAlpha, const int* __restrict__ atomicNum, const float2* __restrict__ dampingAndThole, const real* __restrict__ dirDamping) {
    const unsigned int totalWarps = (blockDim.x*gridDim.x)/TILE_SIZE;
    const unsigned int warp = (blockIdx.x*blockDim.x+threadIdx.x)/TILE_SIZE;
    const unsigned int tgx = threadIdx.x & (TILE_SIZE-1);
    const unsigned int tbx = threadIdx.x - tgx;
    mixed energy = 0;
    __shared__ AtomData localData[THREAD_BLOCK_SIZE];

    // First loop: process tiles that contain exclusions.
    
    const unsigned int firstExclusionTile = FIRST_EXCLUSION_TILE+warp*(LAST_EXCLUSION_TILE-FIRST_EXCLUSION_TILE)/totalWarps;
    const unsigned int lastExclusionTile = FIRST_EXCLUSION_TILE+(warp+1)*(LAST_EXCLUSION_TILE-FIRST_EXCLUSION_TILE)/totalWarps;
    for (int pos = firstExclusionTile; pos < lastExclusionTile; pos++) {
        const ushort2 tileIndices = exclusionTiles[pos];
        const unsigned int x = tileIndices.x;
        const unsigned int y = tileIndices.y;
        AtomData data;
        unsigned int atom1 = x*TILE_SIZE + tgx;
        loadAtomData(data, atom1, posq, sphericalDipole, sphericalQuadrupole, inducedDipole, inducedDipolePolar, dampingAndThole, dirDamping, 
                     labFrameDipole, labFrameQuadrupole, penAlpha, atomicNum);
        data.force = make_real3(0);
        data.torque = make_real3(0);
        uint2 covalent = covalentFlags[pos*TILE_SIZE+tgx];
        unsigned int polarizationGroup = polarizationGroupFlags[pos*TILE_SIZE+tgx];
        if (x == y) {
            // This tile is on the diagonal.

            localData[threadIdx.x].pos = data.pos;
            localData[threadIdx.x].q = data.q;
            localData[threadIdx.x].sphericalDipole = data.sphericalDipole;
            localData[threadIdx.x].dipole = data.dipole;
            localData[threadIdx.x].quadrupoleXX = data.quadrupoleXX;
            localData[threadIdx.x].quadrupoleXY = data.quadrupoleXY;
            localData[threadIdx.x].quadrupoleXZ = data.quadrupoleXZ;
            localData[threadIdx.x].quadrupoleYY = data.quadrupoleYY;
            localData[threadIdx.x].quadrupoleYZ = data.quadrupoleYZ;
#ifdef INCLUDE_QUADRUPOLES
            localData[threadIdx.x].sphericalQuadrupole[0] = data.sphericalQuadrupole[0];
            localData[threadIdx.x].sphericalQuadrupole[1] = data.sphericalQuadrupole[1];
            localData[threadIdx.x].sphericalQuadrupole[2] = data.sphericalQuadrupole[2];
            localData[threadIdx.x].sphericalQuadrupole[3] = data.sphericalQuadrupole[3];
            localData[threadIdx.x].sphericalQuadrupole[4] = data.sphericalQuadrupole[4];
#endif
            localData[threadIdx.x].inducedDipole = data.inducedDipole;
            localData[threadIdx.x].inducedDipolePolar = data.inducedDipolePolar;
            localData[threadIdx.x].penalpha = data.penalpha;
            localData[threadIdx.x].atomic = data.atomic;
            localData[threadIdx.x].thole = data.thole;
            localData[threadIdx.x].damp = data.damp;
            localData[threadIdx.x].dirdamp = data.dirdamp;

            // Compute forces.

            for (unsigned int j = 0; j < TILE_SIZE; j++) {
                int atom2 = y*TILE_SIZE+j;
                if (atom1 != atom2 && atom1 < NUM_ATOMS && atom2 < NUM_ATOMS) {
                    float d = computeDScaleFactor(polarizationGroup, j);
                    float p = computePScaleFactor(covalent, polarizationGroup, j);
                    float m = computeMScaleFactor(covalent, j);
                    computeOneInteraction(data, localData[tbx+j], true, d, p, m, 0.5f, energy);
                    computeOneInteractionCP(data, localData[tbx+j], true, m, 0.5f, energy);
                }
            }
            data.force *= -ENERGY_SCALE_FACTOR;
            data.torque *= ENERGY_SCALE_FACTOR;
            atomicAdd(&forceBuffers[atom1], static_cast<unsigned long long>((long long) (data.force.x*0x100000000)));
            atomicAdd(&forceBuffers[atom1+PADDED_NUM_ATOMS], static_cast<unsigned long long>((long long) (data.force.y*0x100000000)));
            atomicAdd(&forceBuffers[atom1+2*PADDED_NUM_ATOMS], static_cast<unsigned long long>((long long) (data.force.z*0x100000000)));
            atomicAdd(&torqueBuffers[atom1], static_cast<unsigned long long>((long long) (data.torque.x*0x100000000)));
            atomicAdd(&torqueBuffers[atom1+PADDED_NUM_ATOMS], static_cast<unsigned long long>((long long) (data.torque.y*0x100000000)));
            atomicAdd(&torqueBuffers[atom1+2*PADDED_NUM_ATOMS], static_cast<unsigned long long>((long long) (data.torque.z*0x100000000)));
        }
        else {
            // This is an off-diagonal tile.

            unsigned int j = y*TILE_SIZE + tgx;
            loadAtomData(localData[threadIdx.x], j, posq, sphericalDipole, sphericalQuadrupole, inducedDipole, inducedDipolePolar, dampingAndThole, dirDamping,
                         labFrameDipole, labFrameQuadrupole, penAlpha, atomicNum);
            localData[threadIdx.x].force = make_real3(0);
            localData[threadIdx.x].torque = make_real3(0);
            unsigned int tj = tgx;
            for (j = 0; j < TILE_SIZE; j++) {
                int atom2 = y*TILE_SIZE+tj;
                if (atom1 < NUM_ATOMS && atom2 < NUM_ATOMS) {
                    float d = computeDScaleFactor(polarizationGroup, tj);
                    float p = computePScaleFactor(covalent, polarizationGroup, tj);
                    float m = computeMScaleFactor(covalent, tj);
                    computeOneInteraction(data, localData[tbx+tj], true, d, p, m, 1, energy);
                    computeOneInteractionCP(data, localData[tbx+tj], true, m, 1, energy);
                }
                tj = (tj + 1) & (TILE_SIZE - 1);
            }
            data.force *= -ENERGY_SCALE_FACTOR;
            data.torque *= ENERGY_SCALE_FACTOR;
            localData[threadIdx.x].force *= -ENERGY_SCALE_FACTOR;
            localData[threadIdx.x].torque *= ENERGY_SCALE_FACTOR;
            unsigned int offset = x*TILE_SIZE + tgx;
            atomicAdd(&forceBuffers[offset], static_cast<unsigned long long>((long long) (data.force.x*0x100000000)));
            atomicAdd(&forceBuffers[offset+PADDED_NUM_ATOMS], static_cast<unsigned long long>((long long) (data.force.y*0x100000000)));
            atomicAdd(&forceBuffers[offset+2*PADDED_NUM_ATOMS], static_cast<unsigned long long>((long long) (data.force.z*0x100000000)));
            atomicAdd(&torqueBuffers[offset], static_cast<unsigned long long>((long long) (data.torque.x*0x100000000)));
            atomicAdd(&torqueBuffers[offset+PADDED_NUM_ATOMS], static_cast<unsigned long long>((long long) (data.torque.y*0x100000000)));
            atomicAdd(&torqueBuffers[offset+2*PADDED_NUM_ATOMS], static_cast<unsigned long long>((long long) (data.torque.z*0x100000000)));
            offset = y*TILE_SIZE + tgx;
            atomicAdd(&forceBuffers[offset], static_cast<unsigned long long>((long long) (localData[threadIdx.x].force.x*0x100000000)));
            atomicAdd(&forceBuffers[offset+PADDED_NUM_ATOMS], static_cast<unsigned long long>((long long) (localData[threadIdx.x].force.y*0x100000000)));
            atomicAdd(&forceBuffers[offset+2*PADDED_NUM_ATOMS], static_cast<unsigned long long>((long long) (localData[threadIdx.x].force.z*0x100000000)));
            atomicAdd(&torqueBuffers[offset], static_cast<unsigned long long>((long long) (localData[threadIdx.x].torque.x*0x100000000)));
            atomicAdd(&torqueBuffers[offset+PADDED_NUM_ATOMS], static_cast<unsigned long long>((long long) (localData[threadIdx.x].torque.y*0x100000000)));
            atomicAdd(&torqueBuffers[offset+2*PADDED_NUM_ATOMS], static_cast<unsigned long long>((long long) (localData[threadIdx.x].torque.z*0x100000000)));
        }
    }

    // Second loop: tiles without exclusions, either from the neighbor list (with cutoff) or just enumerating all
    // of them (no cutoff).

#ifdef USE_CUTOFF
    const unsigned int numTiles = interactionCount[0];
    if (numTiles > maxTiles)
        return; // There wasn't enough memory for the neighbor list.
    int pos = (int) (numTiles > maxTiles ? startTileIndex+warp*(long long)numTileIndices/totalWarps : warp*(long long)numTiles/totalWarps);
    int end = (int) (numTiles > maxTiles ? startTileIndex+(warp+1)*(long long)numTileIndices/totalWarps : (warp+1)*(long long)numTiles/totalWarps);
#else
    const unsigned int numTiles = numTileIndices;
    int pos = (int) (startTileIndex+warp*(long long)numTiles/totalWarps);
    int end = (int) (startTileIndex+(warp+1)*(long long)numTiles/totalWarps);
#endif
    int skipBase = 0;
    int currentSkipIndex = tbx;
    __shared__ int atomIndices[THREAD_BLOCK_SIZE];
    __shared__ volatile int skipTiles[THREAD_BLOCK_SIZE];
    skipTiles[threadIdx.x] = -1;
    
    while (pos < end) {
        bool includeTile = true;

        // Extract the coordinates of this tile.
        
        int x, y;
#ifdef USE_CUTOFF
        x = tiles[pos];
#else
        y = (int) floor(NUM_BLOCKS+0.5f-SQRT((NUM_BLOCKS+0.5f)*(NUM_BLOCKS+0.5f)-2*pos));
        x = (pos-y*NUM_BLOCKS+y*(y+1)/2);
        if (x < y || x >= NUM_BLOCKS) { // Occasionally happens due to roundoff error.
            y += (x < y ? -1 : 1);
            x = (pos-y*NUM_BLOCKS+y*(y+1)/2);
        }

        // Skip over tiles that have exclusions, since they were already processed.

        while (skipTiles[tbx+TILE_SIZE-1] < pos) {
            if (skipBase+tgx < NUM_TILES_WITH_EXCLUSIONS) {
                ushort2 tile = exclusionTiles[skipBase+tgx];
                skipTiles[threadIdx.x] = tile.x + tile.y*NUM_BLOCKS - tile.y*(tile.y+1)/2;
            }
            else
                skipTiles[threadIdx.x] = end;
            skipBase += TILE_SIZE;            
            currentSkipIndex = tbx;
        }
        while (skipTiles[currentSkipIndex] < pos)
            currentSkipIndex++;
        includeTile = (skipTiles[currentSkipIndex] != pos);
#endif
        if (includeTile) {
            unsigned int atom1 = x*TILE_SIZE + tgx;

            // Load atom data for this tile.

            AtomData data;
            loadAtomData(data, atom1, posq, sphericalDipole, sphericalQuadrupole, inducedDipole, inducedDipolePolar, dampingAndThole, dirDamping, 
                         labFrameDipole, labFrameQuadrupole, penAlpha, atomicNum);
            data.force = make_real3(0);
            data.torque = make_real3(0);
#ifdef USE_CUTOFF
            unsigned int j = interactingAtoms[pos*TILE_SIZE+tgx];
#else
            unsigned int j = y*TILE_SIZE + tgx;
#endif
            atomIndices[threadIdx.x] = j;
            loadAtomData(localData[threadIdx.x], j, posq, sphericalDipole, sphericalQuadrupole, inducedDipole, inducedDipolePolar, dampingAndThole, dirDamping,
                         labFrameDipole, labFrameQuadrupole, penAlpha, atomicNum);
            localData[threadIdx.x].force = make_real3(0);
            localData[threadIdx.x].torque = make_real3(0);

            // Compute forces.

            unsigned int tj = tgx;
            for (j = 0; j < TILE_SIZE; j++) {
                int atom2 = atomIndices[tbx+tj];
                if (atom1 < NUM_ATOMS && atom2 < NUM_ATOMS) {
                    computeOneInteraction(data, localData[tbx+tj], false, 1, 1, 1, 1, energy);
                    computeOneInteractionCP(data, localData[tbx+tj], false, 1, 1, energy);
                }
                tj = (tj + 1) & (TILE_SIZE - 1);
            }
            data.force *= -ENERGY_SCALE_FACTOR;
            data.torque *= ENERGY_SCALE_FACTOR;
            localData[threadIdx.x].force *= -ENERGY_SCALE_FACTOR;
            localData[threadIdx.x].torque *= ENERGY_SCALE_FACTOR;

            // Write results.

            unsigned int offset = x*TILE_SIZE + tgx;
            atomicAdd(&forceBuffers[offset], static_cast<unsigned long long>((long long) (data.force.x*0x100000000)));
            atomicAdd(&forceBuffers[offset+PADDED_NUM_ATOMS], static_cast<unsigned long long>((long long) (data.force.y*0x100000000)));
            atomicAdd(&forceBuffers[offset+2*PADDED_NUM_ATOMS], static_cast<unsigned long long>((long long) (data.force.z*0x100000000)));
            atomicAdd(&torqueBuffers[offset], static_cast<unsigned long long>((long long) (data.torque.x*0x100000000)));
            atomicAdd(&torqueBuffers[offset+PADDED_NUM_ATOMS], static_cast<unsigned long long>((long long) (data.torque.y*0x100000000)));
            atomicAdd(&torqueBuffers[offset+2*PADDED_NUM_ATOMS], static_cast<unsigned long long>((long long) (data.torque.z*0x100000000)));
#ifdef USE_CUTOFF
            offset = atomIndices[threadIdx.x];
#else
            offset = y*TILE_SIZE + tgx;
#endif
            atomicAdd(&forceBuffers[offset], static_cast<unsigned long long>((long long) (localData[threadIdx.x].force.x*0x100000000)));
            atomicAdd(&forceBuffers[offset+PADDED_NUM_ATOMS], static_cast<unsigned long long>((long long) (localData[threadIdx.x].force.y*0x100000000)));
            atomicAdd(&forceBuffers[offset+2*PADDED_NUM_ATOMS], static_cast<unsigned long long>((long long) (localData[threadIdx.x].force.z*0x100000000)));
            atomicAdd(&torqueBuffers[offset], static_cast<unsigned long long>((long long) (localData[threadIdx.x].torque.x*0x100000000)));
            atomicAdd(&torqueBuffers[offset+PADDED_NUM_ATOMS], static_cast<unsigned long long>((long long) (localData[threadIdx.x].torque.y*0x100000000)));
            atomicAdd(&torqueBuffers[offset+2*PADDED_NUM_ATOMS], static_cast<unsigned long long>((long long) (localData[threadIdx.x].torque.z*0x100000000)));
        }
        pos++;
    }
    energyBuffer[blockIdx.x*blockDim.x+threadIdx.x] += energy*ENERGY_SCALE_FACTOR;
}
