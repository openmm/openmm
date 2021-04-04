#if defined F1
__device__ void computeOneEDiffInteractionF1(AtomData4& atom1, volatile AtomData4& atom2, float dScale, float pScale, real& outputEnergy, real3& outputForce) {
#elif defined T1
__device__ void computeOneEDiffInteractionT1(AtomData4& atom1, volatile AtomData4& atom2, float dScale, float pScale, real3& outputForce) {
#elif defined T3
__device__ void computeOneEDiffInteractionT3(AtomData4& atom1, volatile AtomData4& atom2, float dScale, float pScale, real3& outputForce) {
#endif
    const float uscale = 1;

    // deltaR

    real xr = atom2.pos.x - atom1.pos.x;
    real yr = atom2.pos.y - atom1.pos.y;
    real zr = atom2.pos.z - atom1.pos.z;

    real r22 = xr*xr + yr*yr + zr*zr;

    real r = SQRT(r22);
    real rr1 = RECIP(r);
    real rr2 = rr1*rr1;
    real rr3 = rr1*rr2;

    real scale3 = 1;
    real scale5 = 1;
    real scale7 = 1;

#ifdef F1
    real ddsc3_1 = 0;
    real ddsc3_2 = 0;
    real ddsc3_3 = 0;

    real ddsc5_1 = 0;
    real ddsc5_2 = 0;
    real ddsc5_3 = 0;

    real ddsc7_1 = 0;
    real ddsc7_2 = 0;
    real ddsc7_3 = 0;

    real ftm2i1 = 0;
    real ftm2i2 = 0;
    real ftm2i3 = 0;
#endif

    // apply Thole polarization damping to scale factors
 
    real damp = atom1.damp*atom2.damp;
    if (damp != 0) {
        real pgamma = atom2.thole > atom1.thole ? atom1.thole : atom2.thole;
        real ratio = (r/damp);
        damp = -pgamma*ratio*ratio*ratio;
        if (damp > -50) {
            real dampE = EXP(damp);
            real damp2 = damp*damp;
            scale3 = 1 - dampE;
            scale5 = 1 - (1 - damp)*dampE;
            scale7 = 1 - (1 - damp + 0.6f*damp2)*dampE;

#ifdef F1
            ddsc3_1 = -3*damp*EXP(damp)*xr*rr2*rr3;
            ddsc3_2 = -3*damp*EXP(damp)*yr*rr2*rr3;
            ddsc3_3 = -3*damp*EXP(damp)*zr*rr2*rr3;

            ddsc5_1 = -3*damp*ddsc3_1*rr2;
            ddsc5_2 = -3*damp*ddsc3_2*rr2;
            ddsc5_3 = -3*damp*ddsc3_3*rr2;

            ddsc7_1 = -5*(0.2f+0.6f*damp)*ddsc5_1*rr2;
            ddsc7_2 = -5*(0.2f+0.6f*damp)*ddsc5_2*rr2;
            ddsc7_3 = -5*(0.2f+0.6f*damp)*ddsc5_3*rr2;
#endif
        }
    }

    real scale3i = 3*scale3*uscale*rr3*rr2;
    real scale5i = 3*scale5*uscale*rr3*rr2;

    real dsc3 = scale3*dScale*rr3;
    real dsc5 = 3*scale5*dScale*rr3*rr2;
    real dsc7 = 15*scale7*dScale*rr3*rr3*rr1;

    real psc3 = scale3*pScale*rr3;
    real psc5 = 3*scale5*pScale*rr3*rr2;
    real psc7 = 15*scale7*pScale*rr3*rr3*rr1;
 
#ifdef T1
    real dixr1 = atom1.dipole.y*zr - atom1.dipole.z*yr;
    real dixr2 = atom1.dipole.z*xr - atom1.dipole.x*zr;
    real dixr3 = atom1.dipole.x*yr - atom1.dipole.y*xr;
#endif

#ifdef T3
    real dkxr1 = atom2.dipole.y*zr - atom2.dipole.z*yr;
    real dkxr2 = atom2.dipole.z*xr - atom2.dipole.x*zr;
    real dkxr3 = atom2.dipole.x*yr - atom2.dipole.y*xr;
#endif

    real qir1 = atom1.quadrupoleXX*xr + atom1.quadrupoleXY*yr + atom1.quadrupoleXZ*zr;
    real qir2 = atom1.quadrupoleXY*xr + atom1.quadrupoleYY*yr + atom1.quadrupoleYZ*zr;
    real qir3 = atom1.quadrupoleXZ*xr + atom1.quadrupoleYZ*yr + atom1.quadrupoleZZ*zr;

    real qkr1 = atom2.quadrupoleXX*xr + atom2.quadrupoleXY*yr + atom2.quadrupoleXZ*zr;
    real qkr2 = atom2.quadrupoleXY*xr + atom2.quadrupoleYY*yr + atom2.quadrupoleYZ*zr;
    real qkr3 = atom2.quadrupoleXZ*xr + atom2.quadrupoleYZ*yr + atom2.quadrupoleZZ*zr;

#ifdef T1
    real rxqir1 = yr*qir3 - zr*qir2;
    real rxqir2 = zr*qir1 - xr*qir3;
    real rxqir3 = xr*qir2 - yr*qir1;
#endif

#ifdef T3
    real rxqkr1 = yr*qkr3 - zr*qkr2;
    real rxqkr2 = zr*qkr1 - xr*qkr3;
    real rxqkr3 = xr*qkr2 - yr*qkr1;
#endif

    // get intermediate variables for permanent energy terms
 
    real sc3 = atom1.dipole.x*xr + atom1.dipole.y*yr + atom1.dipole.z*zr;
    real sc4 = atom2.dipole.x*xr + atom2.dipole.y*yr + atom2.dipole.z*zr;
    real sc5 = qir1*xr + qir2*yr + qir3*zr;
    real sc6 = qkr1*xr + qkr2*yr + qkr3*zr;
 
#ifdef T1
    real dixuk1 = atom1.dipole.y*atom2.inducedDipoleS.z - atom1.dipole.z*atom2.inducedDipoleS.y;
    real dixuk2 = atom1.dipole.z*atom2.inducedDipoleS.x - atom1.dipole.x*atom2.inducedDipoleS.z;
    real dixuk3 = atom1.dipole.x*atom2.inducedDipoleS.y - atom1.dipole.y*atom2.inducedDipoleS.x;

    real dixukp1 = atom1.dipole.y*atom2.inducedDipolePolarS.z - atom1.dipole.z*atom2.inducedDipolePolarS.y;
    real dixukp2 = atom1.dipole.z*atom2.inducedDipolePolarS.x - atom1.dipole.x*atom2.inducedDipolePolarS.z;
    real dixukp3 = atom1.dipole.x*atom2.inducedDipolePolarS.y - atom1.dipole.y*atom2.inducedDipolePolarS.x;
#endif

#ifdef T3
    real dkxui1 = atom2.dipole.y*atom1.inducedDipoleS.z - atom2.dipole.z*atom1.inducedDipoleS.y;
    real dkxui2 = atom2.dipole.z*atom1.inducedDipoleS.x - atom2.dipole.x*atom1.inducedDipoleS.z;
    real dkxui3 = atom2.dipole.x*atom1.inducedDipoleS.y - atom2.dipole.y*atom1.inducedDipoleS.x;

    real dkxuip1 = atom2.dipole.y*atom1.inducedDipolePolarS.z - atom2.dipole.z*atom1.inducedDipolePolarS.y;
    real dkxuip2 = atom2.dipole.z*atom1.inducedDipolePolarS.x - atom2.dipole.x*atom1.inducedDipolePolarS.z;
    real dkxuip3 = atom2.dipole.x*atom1.inducedDipolePolarS.y - atom2.dipole.y*atom1.inducedDipolePolarS.x;
#endif

#if defined F1 || defined T1
    real qiuk1 = atom1.quadrupoleXX*atom2.inducedDipoleS.x + atom1.quadrupoleXY*atom2.inducedDipoleS.y + atom1.quadrupoleXZ*atom2.inducedDipoleS.z;
    real qiuk2 = atom1.quadrupoleXY*atom2.inducedDipoleS.x + atom1.quadrupoleYY*atom2.inducedDipoleS.y + atom1.quadrupoleYZ*atom2.inducedDipoleS.z;
    real qiuk3 = atom1.quadrupoleXZ*atom2.inducedDipoleS.x + atom1.quadrupoleYZ*atom2.inducedDipoleS.y + atom1.quadrupoleZZ*atom2.inducedDipoleS.z;

    real qiukp1 = atom1.quadrupoleXX*atom2.inducedDipolePolarS.x + atom1.quadrupoleXY*atom2.inducedDipolePolarS.y + atom1.quadrupoleXZ*atom2.inducedDipolePolarS.z;
    real qiukp2 = atom1.quadrupoleXY*atom2.inducedDipolePolarS.x + atom1.quadrupoleYY*atom2.inducedDipolePolarS.y + atom1.quadrupoleYZ*atom2.inducedDipolePolarS.z;
    real qiukp3 = atom1.quadrupoleXZ*atom2.inducedDipolePolarS.x + atom1.quadrupoleYZ*atom2.inducedDipolePolarS.y + atom1.quadrupoleZZ*atom2.inducedDipolePolarS.z;
#if defined F1
    qiuk1 -= atom1.quadrupoleXX*atom2.inducedDipole.x + atom1.quadrupoleXY*atom2.inducedDipole.y + atom1.quadrupoleXZ*atom2.inducedDipole.z;
    qiuk2 -= atom1.quadrupoleXY*atom2.inducedDipole.x + atom1.quadrupoleYY*atom2.inducedDipole.y + atom1.quadrupoleYZ*atom2.inducedDipole.z;
    qiuk3 -= atom1.quadrupoleXZ*atom2.inducedDipole.x + atom1.quadrupoleYZ*atom2.inducedDipole.y + atom1.quadrupoleZZ*atom2.inducedDipole.z;

    qiukp1 -= atom1.quadrupoleXX*atom2.inducedDipolePolar.x + atom1.quadrupoleXY*atom2.inducedDipolePolar.y + atom1.quadrupoleXZ*atom2.inducedDipolePolar.z;
    qiukp2 -= atom1.quadrupoleXY*atom2.inducedDipolePolar.x + atom1.quadrupoleYY*atom2.inducedDipolePolar.y + atom1.quadrupoleYZ*atom2.inducedDipolePolar.z;
    qiukp3 -= atom1.quadrupoleXZ*atom2.inducedDipolePolar.x + atom1.quadrupoleYZ*atom2.inducedDipolePolar.y + atom1.quadrupoleZZ*atom2.inducedDipolePolar.z;

    ftm2i1 -= psc5*qiuk1 + dsc5*qiukp1;
    ftm2i2 -= psc5*qiuk2 + dsc5*qiukp2;
    ftm2i3 -= psc5*qiuk3 + dsc5*qiukp3;
#endif
#endif

#if defined F1 || defined T3
    real qkui1 = atom2.quadrupoleXX*atom1.inducedDipoleS.x + atom2.quadrupoleXY*atom1.inducedDipoleS.y + atom2.quadrupoleXZ*atom1.inducedDipoleS.z;
    real qkui2 = atom2.quadrupoleXY*atom1.inducedDipoleS.x + atom2.quadrupoleYY*atom1.inducedDipoleS.y + atom2.quadrupoleYZ*atom1.inducedDipoleS.z;
    real qkui3 = atom2.quadrupoleXZ*atom1.inducedDipoleS.x + atom2.quadrupoleYZ*atom1.inducedDipoleS.y + atom2.quadrupoleZZ*atom1.inducedDipoleS.z;

    real qkuip1 = atom2.quadrupoleXX*atom1.inducedDipolePolarS.x + atom2.quadrupoleXY*atom1.inducedDipolePolarS.y + atom2.quadrupoleXZ*atom1.inducedDipolePolarS.z;
    real qkuip2 = atom2.quadrupoleXY*atom1.inducedDipolePolarS.x + atom2.quadrupoleYY*atom1.inducedDipolePolarS.y + atom2.quadrupoleYZ*atom1.inducedDipolePolarS.z;
    real qkuip3 = atom2.quadrupoleXZ*atom1.inducedDipolePolarS.x + atom2.quadrupoleYZ*atom1.inducedDipolePolarS.y + atom2.quadrupoleZZ*atom1.inducedDipolePolarS.z;

#if defined F1
    qkui1 -= atom2.quadrupoleXX*atom1.inducedDipole.x + atom2.quadrupoleXY*atom1.inducedDipole.y + atom2.quadrupoleXZ*atom1.inducedDipole.z;
    qkui2 -= atom2.quadrupoleXY*atom1.inducedDipole.x + atom2.quadrupoleYY*atom1.inducedDipole.y + atom2.quadrupoleYZ*atom1.inducedDipole.z;
    qkui3 -= atom2.quadrupoleXZ*atom1.inducedDipole.x + atom2.quadrupoleYZ*atom1.inducedDipole.y + atom2.quadrupoleZZ*atom1.inducedDipole.z;

    qkuip1 -= atom2.quadrupoleXX*atom1.inducedDipolePolar.x + atom2.quadrupoleXY*atom1.inducedDipolePolar.y + atom2.quadrupoleXZ*atom1.inducedDipolePolar.z;
    qkuip2 -= atom2.quadrupoleXY*atom1.inducedDipolePolar.x + atom2.quadrupoleYY*atom1.inducedDipolePolar.y + atom2.quadrupoleYZ*atom1.inducedDipolePolar.z;
    qkuip3 -= atom2.quadrupoleXZ*atom1.inducedDipolePolar.x + atom2.quadrupoleYZ*atom1.inducedDipolePolar.y + atom2.quadrupoleZZ*atom1.inducedDipolePolar.z;

    ftm2i1 += psc5*qkui1 + dsc5*qkuip1;
    ftm2i2 += psc5*qkui2 + dsc5*qkuip2;
    ftm2i3 += psc5*qkui3 + dsc5*qkuip3;
#endif

#endif

#ifdef T3
    real uixqkr1 = atom1.inducedDipoleS.y*qkr3 - atom1.inducedDipoleS.z*qkr2;
    real uixqkr2 = atom1.inducedDipoleS.z*qkr1 - atom1.inducedDipoleS.x*qkr3;
    real uixqkr3 = atom1.inducedDipoleS.x*qkr2 - atom1.inducedDipoleS.y*qkr1;

    real uixqkrp1 = atom1.inducedDipolePolarS.y*qkr3 - atom1.inducedDipolePolarS.z*qkr2;
    real uixqkrp2 = atom1.inducedDipolePolarS.z*qkr1 - atom1.inducedDipolePolarS.x*qkr3;
    real uixqkrp3 = atom1.inducedDipolePolarS.x*qkr2 - atom1.inducedDipolePolarS.y*qkr1;

    real rxqkuip1 = yr*qkuip3 - zr*qkuip2;
    real rxqkuip2 = zr*qkuip1 - xr*qkuip3;
    real rxqkuip3 = xr*qkuip2 - yr*qkuip1;

    real rxqkui1 = yr*qkui3 - zr*qkui2;
    real rxqkui2 = zr*qkui1 - xr*qkui3;
    real rxqkui3 = xr*qkui2 - yr*qkui1;
#endif

#ifdef T1
    real ukxqir1 = atom2.inducedDipoleS.y*qir3 - atom2.inducedDipoleS.z*qir2;
    real ukxqir2 = atom2.inducedDipoleS.z*qir1 - atom2.inducedDipoleS.x*qir3;
    real ukxqir3 = atom2.inducedDipoleS.x*qir2 - atom2.inducedDipoleS.y*qir1;

    real ukxqirp1 = atom2.inducedDipolePolarS.y*qir3 - atom2.inducedDipolePolarS.z*qir2;
    real ukxqirp2 = atom2.inducedDipolePolarS.z*qir1 - atom2.inducedDipolePolarS.x*qir3;
    real ukxqirp3 = atom2.inducedDipolePolarS.x*qir2 - atom2.inducedDipolePolarS.y*qir1;

    real rxqiuk1 = yr*qiuk3 - zr*qiuk2;
    real rxqiuk2 = zr*qiuk1 - xr*qiuk3;
    real rxqiuk3 = xr*qiuk2 - yr*qiuk1;

    real rxqiukp1 = yr*qiukp3 - zr*qiukp2;
    real rxqiukp2 = zr*qiukp1 - xr*qiukp3;
    real rxqiukp3 = xr*qiukp2 - yr*qiukp1;
#endif

    // get intermediate variables for induction energy terms

    real sci3 = atom1.inducedDipoleS.x*xr + atom1.inducedDipoleS.y*yr + atom1.inducedDipoleS.z*zr;
    real sci4 = atom2.inducedDipoleS.x*xr + atom2.inducedDipoleS.y*yr + atom2.inducedDipoleS.z*zr;
#ifdef F1
    ftm2i1 += 0.5f*scale5i*(sci4*atom1.inducedDipolePolarS.x + sci3*atom2.inducedDipolePolarS.x);
    ftm2i2 += 0.5f*scale5i*(sci4*atom1.inducedDipolePolarS.y + sci3*atom2.inducedDipolePolarS.y);
    ftm2i3 += 0.5f*scale5i*(sci4*atom1.inducedDipolePolarS.z + sci3*atom2.inducedDipolePolarS.z);
#endif
    real sci3Y = atom1.inducedDipole.x*xr + atom1.inducedDipole.y*yr + atom1.inducedDipole.z*zr;
    real sci4Y = atom2.inducedDipole.x*xr + atom2.inducedDipole.y*yr + atom2.inducedDipole.z*zr;
#ifdef F1
    ftm2i1 -= 0.5f*scale5i*(sci3Y*atom2.inducedDipolePolar.x + sci4Y*atom1.inducedDipolePolar.x);
    ftm2i2 -= 0.5f*scale5i*(sci3Y*atom2.inducedDipolePolar.y + sci4Y*atom1.inducedDipolePolar.y);
    ftm2i3 -= 0.5f*scale5i*(sci3Y*atom2.inducedDipolePolar.z + sci4Y*atom1.inducedDipolePolar.z);
#endif

    real sci7 = qir1*atom2.inducedDipoleS.x + qir2*atom2.inducedDipoleS.y + qir3*atom2.inducedDipoleS.z;
    real sci8 = qkr1*atom1.inducedDipoleS.x + qkr2*atom1.inducedDipoleS.y + qkr3*atom1.inducedDipoleS.z;
    real scip1 = atom1.inducedDipolePolarS.x*atom2.dipole.x + atom1.inducedDipolePolarS.y*atom2.dipole.y + atom1.inducedDipolePolarS.z*atom2.dipole.z +
                              atom1.dipole.x*atom2.inducedDipolePolarS.x + atom1.dipole.y*atom2.inducedDipolePolarS.y + atom1.dipole.z*atom2.inducedDipolePolarS.z;

    real scip2 = atom1.inducedDipoleS.x*atom2.inducedDipolePolarS.x + atom1.inducedDipoleS.y*atom2.inducedDipolePolarS.y + atom1.inducedDipoleS.z*atom2.inducedDipolePolarS.z +
                              atom1.inducedDipolePolarS.x*atom2.inducedDipoleS.x + atom1.inducedDipolePolarS.y*atom2.inducedDipoleS.y + atom1.inducedDipolePolarS.z*atom2.inducedDipoleS.z;

    sci7 -= qir1*atom2.inducedDipole.x + qir2*atom2.inducedDipole.y + qir3*atom2.inducedDipole.z;
    sci8 -= qkr1*atom1.inducedDipole.x + qkr2*atom1.inducedDipole.y + qkr3*atom1.inducedDipole.z;

    scip1 -= atom1.inducedDipolePolar.x*atom2.dipole.x + atom1.inducedDipolePolar.y*atom2.dipole.y + atom1.inducedDipolePolar.z*atom2.dipole.z +
                              atom1.dipole.x*atom2.inducedDipolePolar.x + atom1.dipole.y*atom2.inducedDipolePolar.y + atom1.dipole.z*atom2.inducedDipolePolar.z;


    scip2 -= atom1.inducedDipole.x*atom2.inducedDipolePolar.x + atom1.inducedDipole.y*atom2.inducedDipolePolar.y + atom1.inducedDipole.z*atom2.inducedDipolePolar.z   +
                              atom1.inducedDipolePolar.x*atom2.inducedDipole.x + atom1.inducedDipolePolar.y*atom2.inducedDipole.y + atom1.inducedDipolePolar.z*atom2.inducedDipole.z;


    real scip3 = atom1.inducedDipolePolarS.x*xr + atom1.inducedDipolePolarS.y*yr + atom1.inducedDipolePolarS.z*zr;
    real scip4 = atom2.inducedDipolePolarS.x*xr + atom2.inducedDipolePolarS.y*yr + atom2.inducedDipolePolarS.z*zr;
    real gfi1 = -2.5f*(sci3*scip4+scip3*sci4)*scale5i;

#ifdef F1
    ftm2i1 += 0.5f*scale5i*(scip4*atom1.inducedDipoleS.x + scip3*atom2.inducedDipoleS.x);
    ftm2i2 += 0.5f*scale5i*(scip4*atom1.inducedDipoleS.y + scip3*atom2.inducedDipoleS.y);
    ftm2i3 += 0.5f*scale5i*(scip4*atom1.inducedDipoleS.z + scip3*atom2.inducedDipoleS.z);
#endif

    real scip3Y = atom1.inducedDipolePolar.x*xr + atom1.inducedDipolePolar.y*yr + atom1.inducedDipolePolar.z*zr;
    real scip4Y = atom2.inducedDipolePolar.x*xr + atom2.inducedDipolePolar.y*yr + atom2.inducedDipolePolar.z*zr;
    gfi1 += 2.5f*(sci3Y*scip4Y + scip3Y*sci4Y)*scale5i;
#ifdef F1
    ftm2i1 -= 0.5f*scale5i*(scip3Y*atom2.inducedDipole.x + scip4Y*atom1.inducedDipole.x);
    ftm2i2 -= 0.5f*scale5i*(scip3Y*atom2.inducedDipole.y + scip4Y*atom1.inducedDipole.y);
    ftm2i3 -= 0.5f*scale5i*(scip3Y*atom2.inducedDipole.z + scip4Y*atom1.inducedDipole.z);
#endif
    sci3Y = sci3 - sci3Y;
    sci4Y = sci4 - sci4Y;
    scip3Y = scip3 - scip3Y;
    scip4Y = scip4 - scip4Y;

    real scip7 = qir1*atom2.inducedDipolePolarS.x + qir2*atom2.inducedDipolePolarS.y + qir3*atom2.inducedDipolePolarS.z;
    scip7 -= qir1*atom2.inducedDipolePolar.x + qir2*atom2.inducedDipolePolar.y + qir3*atom2.inducedDipolePolar.z;

    real scip8 = qkr1*atom1.inducedDipolePolarS.x + qkr2*atom1.inducedDipolePolarS.y + qkr3*atom1.inducedDipolePolarS.z;
    scip8 -= qkr1*atom1.inducedDipolePolar.x + qkr2*atom1.inducedDipolePolar.y + qkr3*atom1.inducedDipolePolar.z;

    real sci1 = atom1.inducedDipoleS.x*atom2.dipole.x + atom1.inducedDipoleS.y*atom2.dipole.y +
                              atom1.inducedDipoleS.z*atom2.dipole.z + atom1.dipole.x*atom2.inducedDipoleS.x +
                              atom1.dipole.y*atom2.inducedDipoleS.y + atom1.dipole.z*atom2.inducedDipoleS.z;
    sci1 -= atom1.inducedDipole.x*atom2.dipole.x + atom1.inducedDipole.y*atom2.dipole.y +
                              atom1.inducedDipole.z*atom2.dipole.z + atom1.dipole.x*atom2.inducedDipole.x +
                              atom1.dipole.y*atom2.inducedDipole.y + atom1.dipole.z*atom2.inducedDipole.z;

    real gli1 = atom2.q*sci3Y - atom1.q*sci4Y + sci1;
    real gli2 = -sc3*sci4Y - sci3Y*sc4 + 2*(sci7-sci8);
    real gli3 = sci3Y*sc6 - sci4Y*sc5;
    real glip1 = atom2.q*scip3Y - atom1.q*scip4Y + scip1;
    real glip2 = -sc3*scip4Y - scip3Y*sc4 + 2*(scip7-scip8);
    real glip3 = scip3Y*sc6 - scip4Y*sc5;

#ifdef F1
    ftm2i1 -= 0.5f*((gli1*pScale + glip1*dScale)*ddsc3_1 + (gli2*pScale + glip2*dScale)*ddsc5_1 + (gli3*pScale+glip3*dScale)*ddsc7_1);
    ftm2i2 -= 0.5f*((gli1*pScale + glip1*dScale)*ddsc3_2 + (gli2*pScale + glip2*dScale)*ddsc5_2 + (gli3*pScale+glip3*dScale)*ddsc7_2);
    ftm2i3 -= 0.5f*((gli1*pScale + glip1*dScale)*ddsc3_3 + (gli2*pScale + glip2*dScale)*ddsc5_3 + (gli3*pScale+glip3*dScale)*ddsc7_3);
    outputEnergy = gli1*psc3 + gli2*psc5 + gli3*psc7;
#endif

    gfi1 += 1.5f*(gli1*psc3 + glip1*dsc3);
    gfi1 += 2.5f*(gli2*psc5 + glip2*dsc5);
    gfi1 += 3.5f*(gli3*psc7 + glip3*dsc7);
    gfi1 *= rr2;
    gfi1 += 0.5f*scip2*scale3i;

#if defined F1 || defined T1
    real gfi5 =  (sci4Y*psc7+scip4Y*dsc7);
#endif

#if defined F1 || defined T3
    real gfi6 = -(sci3Y*psc7+scip3Y*dsc7);
#endif

#ifdef F1
    ftm2i1 += gfi1*xr;

    real diff0 = atom1.inducedDipoleS.x - atom1.inducedDipole.x;               
    real diff1 = atom1.inducedDipolePolarS.x - atom1.inducedDipolePolar.x;               
    ftm2i1 += 0.5f*(-atom2.q*(diff0*psc3 + diff1*dsc3) + sc4*(diff0*psc5 + diff1*dsc5) - sc6*(diff0*psc7 + diff1*dsc7));
    
    diff0 = atom2.inducedDipoleS.x - atom2.inducedDipole.x;               
    diff1 = atom2.inducedDipolePolarS.x - atom2.inducedDipolePolar.x;               
    ftm2i1 += 0.5f*(atom1.q*(diff0*psc3 + diff1*dsc3) + sc3*(diff0*psc5 + diff1*dsc5) + sc5*(diff0*psc7 + diff1*dsc7));
    ftm2i1 += 0.5f*(sci4Y*psc5+scip4Y*dsc5)*atom1.dipole.x + 0.5f*(sci3Y*psc5+scip3Y*dsc5)*atom2.dipole.x + gfi5*qir1 + gfi6*qkr1;

    ftm2i2 += gfi1*yr;

    diff0 = atom1.inducedDipoleS.y - atom1.inducedDipole.y;               
    diff1 = atom1.inducedDipolePolarS.y - atom1.inducedDipolePolar.y;               
    ftm2i2 += 0.5f*(-atom2.q*(diff0*psc3 + diff1*dsc3) + sc4*(diff0*psc5 + diff1*dsc5) - sc6*(diff0*psc7 + diff1*dsc7));

    diff0 = atom2.inducedDipoleS.y - atom2.inducedDipole.y;               
    diff1 = atom2.inducedDipolePolarS.y - atom2.inducedDipolePolar.y;               

    ftm2i2 += 0.5f*(atom1.q*(diff0*psc3 + diff1*dsc3) + sc3*(diff0*psc5 + diff1*dsc5) + sc5*(diff0*psc7 + diff1*dsc7));
    ftm2i2 += 0.5f*(sci4Y*psc5+scip4Y*dsc5)*atom1.dipole.y + 0.5f*(sci3Y*psc5+scip3Y*dsc5)*atom2.dipole.y + gfi5*qir2 + gfi6*qkr2;

    ftm2i3 += gfi1*zr;

    diff0 = atom1.inducedDipoleS.z - atom1.inducedDipole.z;               
    diff1 = atom1.inducedDipolePolarS.z - atom1.inducedDipolePolar.z;               
    ftm2i3 += 0.5f*(-atom2.q*(diff0*psc3 + diff1*dsc3) + sc4*(diff0*psc5 + diff1*dsc5) - sc6*(diff0*psc7 + diff1*dsc7));
    diff0 = atom2.inducedDipoleS.z - atom2.inducedDipole.z;               
    diff1 = atom2.inducedDipolePolarS.z - atom2.inducedDipolePolar.z;               

    ftm2i3 += 0.5f*(atom1.q*(diff0*psc3 + diff1*dsc3) + sc3*(diff0*psc5 + diff1*dsc5) + sc5*(diff0*psc7 + diff1*dsc7));
    ftm2i3 += 0.5f*(sci4Y*psc5+scip4Y*dsc5)*atom1.dipole.z + 0.5f*(sci3Y*psc5+scip3Y*dsc5)*atom2.dipole.z + gfi5*qir3 + gfi6*qkr3;
 
    // intermediate values needed for partially excluded interactions

    // correction to convert mutual to direct polarization force

#ifdef MUTUAL_POLARIZATION
    real findmp1 = uscale*(scip2*ddsc3_1 - ddsc5_1*(sci3*scip4+scip3*sci4));
    real findmp2 = uscale*(scip2*ddsc3_2 - ddsc5_2*(sci3*scip4+scip3*sci4));
    real findmp3 = uscale*(scip2*ddsc3_3 - ddsc5_3*(sci3*scip4+scip3*sci4));
    ftm2i1 -= 0.5f*findmp1;
    ftm2i2 -= 0.5f*findmp2;
    ftm2i3 -= 0.5f*findmp3;

    real sci3X = sci3 - sci3Y;
    real sci4X = sci4 - sci4Y;
    real scip3X = scip3 - scip3Y;
    real scip4X = scip4 - scip4Y;
    ftm2i1 += 0.5f*uscale*(-ddsc5_1*(sci3X*scip4X+scip3X*sci4X));
    ftm2i2 += 0.5f*uscale*(-ddsc5_2*(sci3X*scip4X+scip3X*sci4X));
    ftm2i3 += 0.5f*uscale*(-ddsc5_3*(sci3X*scip4X+scip3X*sci4X)); 
#else
    real gfd = (scip2*scale3i - 5*rr2*(scip3*sci4+sci3*scip4)*scale5i);
    real fdir1 = gfd*xr + scale5i* (sci4*atom1.inducedDipolePolarS.x+scip4*atom1.inducedDipoleS.x + sci3*atom2.inducedDipolePolarS.x+scip3*atom2.inducedDipoleS.x);
    real fdir2 = gfd*yr + scale5i* (sci4*atom1.inducedDipolePolarS.y+scip4*atom1.inducedDipoleS.y + sci3*atom2.inducedDipolePolarS.y+scip3*atom2.inducedDipoleS.y);
    real fdir3 = gfd*zr + scale5i* (sci4*atom1.inducedDipolePolarS.z+scip4*atom1.inducedDipoleS.z + sci3*atom2.inducedDipolePolarS.z+scip3*atom2.inducedDipoleS.z);
    ftm2i1 -= 0.5f*fdir1;
    ftm2i2 -= 0.5f*fdir2;
    ftm2i3 -= 0.5f*fdir3;

    real sci3X = sci3 - sci3Y;
    real sci4X = sci4 - sci4Y;
    real scip3X = scip3 - scip3Y;
    real scip4X = scip4 - scip4Y;
    gfd = -5*rr2*(scip3X*sci4X+sci3X*scip4X)*scale5i;
    fdir1 = gfd*xr + scale5i*(sci4X*atom1.inducedDipolePolar.x + scip4X*atom1.inducedDipole.x + sci3X*atom2.inducedDipolePolar.x + scip3X*atom2.inducedDipole.x);
    fdir2 = gfd*yr + scale5i*(sci4X*atom1.inducedDipolePolar.y + scip4X*atom1.inducedDipole.y + sci3X*atom2.inducedDipolePolar.y + scip3X*atom2.inducedDipole.y);
    fdir3 = gfd*zr + scale5i*(sci4X*atom1.inducedDipolePolar.z + scip4X*atom1.inducedDipole.z + sci3X*atom2.inducedDipolePolar.z + scip3X*atom2.inducedDipole.z);
    ftm2i1 += 0.5f*fdir1;
    ftm2i2 += 0.5f*fdir2;
    ftm2i3 += 0.5f*fdir3;
#endif
#endif

#ifdef T1
    real gti2 = 0.5f*(sci4Y*psc5 + scip4Y*dsc5);
    real ttm2i1 = -(dixuk1*psc3+dixukp1*dsc3)*0.5f + gti2*dixr1 + ((ukxqir1+rxqiuk1)*psc5 +(ukxqirp1+rxqiukp1)*dsc5) - gfi5*rxqir1;
    real ttm2i2 = -(dixuk2*psc3+dixukp2*dsc3)*0.5f + gti2*dixr2 + ((ukxqir2+rxqiuk2)*psc5 +(ukxqirp2+rxqiukp2)*dsc5) - gfi5*rxqir2;
    real ttm2i3 = -(dixuk3*psc3+dixukp3*dsc3)*0.5f + gti2*dixr3 + ((ukxqir3+rxqiuk3)*psc5 +(ukxqirp3+rxqiukp3)*dsc5) - gfi5*rxqir3;
#endif

#ifdef T3
    real gti3 = 0.5f*(sci3Y*psc5 + scip3Y*dsc5);
    real ttm3i1 = -(dkxui1*psc3+dkxuip1*dsc3)*0.5f + gti3*dkxr1 - ((uixqkr1+rxqkui1)*psc5 +(uixqkrp1+rxqkuip1)*dsc5) - gfi6*rxqkr1;
    real ttm3i2 = -(dkxui2*psc3+dkxuip2*dsc3)*0.5f + gti3*dkxr2 - ((uixqkr2+rxqkui2)*psc5 +(uixqkrp2+rxqkuip2)*dsc5) - gfi6*rxqkr2;
    real ttm3i3 = -(dkxui3*psc3+dkxuip3*dsc3)*0.5f + gti3*dkxr3 - ((uixqkr3+rxqkui3)*psc5 +(uixqkrp3+rxqkuip3)*dsc5) - gfi6*rxqkr3;
#endif
 
    // update force and torque on site k
    
#ifdef F1
    outputForce.x = -ftm2i1;
    outputForce.y = -ftm2i2;
    outputForce.z = -ftm2i3;
#endif

#ifdef T1
    outputForce.x = ttm2i1;
    outputForce.y = ttm2i2;
    outputForce.z = ttm2i3;
#endif

#ifdef T3
    outputForce.x = ttm3i1;
    outputForce.y = ttm3i2;
    outputForce.z = ttm3i3;
#endif

    // construct auxiliary vectors for induced terms

#ifdef T1
    dixuk1 = atom1.dipole.y*atom2.inducedDipole.z - atom1.dipole.z*atom2.inducedDipole.y;
    dixuk2 = atom1.dipole.z*atom2.inducedDipole.x - atom1.dipole.x*atom2.inducedDipole.z;
    dixuk3 = atom1.dipole.x*atom2.inducedDipole.y - atom1.dipole.y*atom2.inducedDipole.x;

    dixukp1 = atom1.dipole.y*atom2.inducedDipolePolar.z - atom1.dipole.z*atom2.inducedDipolePolar.y;
    dixukp2 = atom1.dipole.z*atom2.inducedDipolePolar.x - atom1.dipole.x*atom2.inducedDipolePolar.z;
    dixukp3 = atom1.dipole.x*atom2.inducedDipolePolar.y - atom1.dipole.y*atom2.inducedDipolePolar.x;
#endif

#ifdef T3
    dkxui1 = atom2.dipole.y*atom1.inducedDipole.z - atom2.dipole.z*atom1.inducedDipole.y;
    dkxui2 = atom2.dipole.z*atom1.inducedDipole.x - atom2.dipole.x*atom1.inducedDipole.z;
    dkxui3 = atom2.dipole.x*atom1.inducedDipole.y - atom2.dipole.y*atom1.inducedDipole.x;

    dkxuip1 = atom2.dipole.y*atom1.inducedDipolePolar.z - atom2.dipole.z*atom1.inducedDipolePolar.y;
    dkxuip2 = atom2.dipole.z*atom1.inducedDipolePolar.x - atom2.dipole.x*atom1.inducedDipolePolar.z;
    dkxuip3 = atom2.dipole.x*atom1.inducedDipolePolar.y - atom2.dipole.y*atom1.inducedDipolePolar.x;
#endif

#if defined T1
    qiuk1 = atom1.quadrupoleXX*atom2.inducedDipole.x + atom1.quadrupoleXY*atom2.inducedDipole.y + atom1.quadrupoleXZ*atom2.inducedDipole.z;
    qiuk2 = atom1.quadrupoleXY*atom2.inducedDipole.x + atom1.quadrupoleYY*atom2.inducedDipole.y + atom1.quadrupoleYZ*atom2.inducedDipole.z;
    qiuk3 = atom1.quadrupoleXZ*atom2.inducedDipole.x + atom1.quadrupoleYZ*atom2.inducedDipole.y + atom1.quadrupoleZZ*atom2.inducedDipole.z;

    qiukp1 = atom1.quadrupoleXX*atom2.inducedDipolePolar.x + atom1.quadrupoleXY*atom2.inducedDipolePolar.y + atom1.quadrupoleXZ*atom2.inducedDipolePolar.z;
    qiukp2 = atom1.quadrupoleXY*atom2.inducedDipolePolar.x + atom1.quadrupoleYY*atom2.inducedDipolePolar.y + atom1.quadrupoleYZ*atom2.inducedDipolePolar.z;
    qiukp3 = atom1.quadrupoleXZ*atom2.inducedDipolePolar.x + atom1.quadrupoleYZ*atom2.inducedDipolePolar.y + atom1.quadrupoleZZ*atom2.inducedDipolePolar.z;
#endif

#if defined T3
    qkui1 = atom2.quadrupoleXX*atom1.inducedDipole.x + atom2.quadrupoleXY*atom1.inducedDipole.y + atom2.quadrupoleXZ*atom1.inducedDipole.z;
    qkui2 = atom2.quadrupoleXY*atom1.inducedDipole.x + atom2.quadrupoleYY*atom1.inducedDipole.y + atom2.quadrupoleYZ*atom1.inducedDipole.z;
    qkui3 = atom2.quadrupoleXZ*atom1.inducedDipole.x + atom2.quadrupoleYZ*atom1.inducedDipole.y + atom2.quadrupoleZZ*atom1.inducedDipole.z;

    qkuip1 = atom2.quadrupoleXX*atom1.inducedDipolePolar.x + atom2.quadrupoleXY*atom1.inducedDipolePolar.y + atom2.quadrupoleXZ*atom1.inducedDipolePolar.z;
    qkuip2 = atom2.quadrupoleXY*atom1.inducedDipolePolar.x + atom2.quadrupoleYY*atom1.inducedDipolePolar.y + atom2.quadrupoleYZ*atom1.inducedDipolePolar.z;
    qkuip3 = atom2.quadrupoleXZ*atom1.inducedDipolePolar.x + atom2.quadrupoleYZ*atom1.inducedDipolePolar.y + atom2.quadrupoleZZ*atom1.inducedDipolePolar.z;
#endif

#ifdef T3
    uixqkr1 = atom1.inducedDipole.y*qkr3 - atom1.inducedDipole.z*qkr2;
    uixqkr2 = atom1.inducedDipole.z*qkr1 - atom1.inducedDipole.x*qkr3;
    uixqkr3 = atom1.inducedDipole.x*qkr2 - atom1.inducedDipole.y*qkr1;

    uixqkrp1 = atom1.inducedDipolePolar.y*qkr3 - atom1.inducedDipolePolar.z*qkr2;
    uixqkrp2 = atom1.inducedDipolePolar.z*qkr1 - atom1.inducedDipolePolar.x*qkr3;
    uixqkrp3 = atom1.inducedDipolePolar.x*qkr2 - atom1.inducedDipolePolar.y*qkr1;
#endif

#ifdef T1
    ukxqir1 = atom2.inducedDipole.y*qir3 - atom2.inducedDipole.z*qir2;
    ukxqir2 = atom2.inducedDipole.z*qir1 - atom2.inducedDipole.x*qir3;
    ukxqir3 = atom2.inducedDipole.x*qir2 - atom2.inducedDipole.y*qir1;

    ukxqirp1 = atom2.inducedDipolePolar.y*qir3 - atom2.inducedDipolePolar.z*qir2;
    ukxqirp2 = atom2.inducedDipolePolar.z*qir1 - atom2.inducedDipolePolar.x*qir3;
    ukxqirp3 = atom2.inducedDipolePolar.x*qir2 - atom2.inducedDipolePolar.y*qir1;

    rxqiuk1 = yr*qiuk3 - zr*qiuk2;
    rxqiuk2 = zr*qiuk1 - xr*qiuk3;
    rxqiuk3 = xr*qiuk2 - yr*qiuk1;

    rxqiukp1 = yr*qiukp3 - zr*qiukp2;
    rxqiukp2 = zr*qiukp1 - xr*qiukp3;
    rxqiukp3 = xr*qiukp2 - yr*qiukp1;
#endif

#ifdef T3
    rxqkui1 = yr*qkui3 - zr*qkui2;
    rxqkui2 = zr*qkui1 - xr*qkui3;
    rxqkui3 = xr*qkui2 - yr*qkui1;

    rxqkuip1 = yr*qkuip3 - zr*qkuip2;
    rxqkuip2 = zr*qkuip1 - xr*qkuip3;
    rxqkuip3 = xr*qkuip2 - yr*qkuip1;
#endif

#ifdef T1
     ttm2i1 = -(dixuk1*psc3+dixukp1*dsc3)*0.5f + ((ukxqir1+rxqiuk1)*psc5 +(ukxqirp1+rxqiukp1)*dsc5);
     ttm2i2 = -(dixuk2*psc3+dixukp2*dsc3)*0.5f + ((ukxqir2+rxqiuk2)*psc5 +(ukxqirp2+rxqiukp2)*dsc5);
     ttm2i3 = -(dixuk3*psc3+dixukp3*dsc3)*0.5f + ((ukxqir3+rxqiuk3)*psc5 +(ukxqirp3+rxqiukp3)*dsc5);
#endif
 
#ifdef T3
     ttm3i1 = -(dkxui1*psc3+dkxuip1*dsc3)*0.5f - ((uixqkr1+rxqkui1)*psc5 +(uixqkrp1+rxqkuip1)*dsc5);
     ttm3i2 = -(dkxui2*psc3+dkxuip2*dsc3)*0.5f - ((uixqkr2+rxqkui2)*psc5 +(uixqkrp2+rxqkuip2)*dsc5);
     ttm3i3 = -(dkxui3*psc3+dkxuip3*dsc3)*0.5f - ((uixqkr3+rxqkui3)*psc5 +(uixqkrp3+rxqkuip3)*dsc5);
#endif

    // update force and torque on site k;

#ifdef T1
    outputForce.x -= ttm2i1;
    outputForce.y -= ttm2i2;
    outputForce.z -= ttm2i3;
#endif

#ifdef T3
    outputForce.x -= ttm3i1;
    outputForce.y -= ttm3i2;
    outputForce.z -= ttm3i3;
#endif
}
