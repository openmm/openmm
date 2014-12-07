__device__ void
#ifdef APPLY_SCALE
computeOneInteractionF1(
#else
computeOneInteractionF1NoScale(
#endif
        AtomData& atom1, volatile AtomData& atom2, real4 delta, real4 bn, real bn5, float forceFactor,
#ifdef APPLY_SCALE
        float dScale, float pScale, float mScale,
#endif
        real3& force, real& energy) {
    real xr = delta.x;
    real yr = delta.y;
    real zr = delta.z;
#ifdef APPLY_SCALE
    real rr1 = delta.w;
#endif

    // set the permanent multipole and induced dipole values;

    real ci = atom1.q;

    real di1 = atom1.dipole.x;
    real di2 = atom1.dipole.y;
    real di3 = atom1.dipole.z;

    real qi1 = atom1.quadrupoleXX;
    real qi2 = atom1.quadrupoleXY;
    real qi3 = atom1.quadrupoleXZ;
    real qi5 = atom1.quadrupoleYY;
    real qi6 = atom1.quadrupoleYZ;
    real qi9 = -(atom1.quadrupoleXX + atom1.quadrupoleYY);

    real ck = atom2.q;
    real dk1 = atom2.dipole.x;
    real dk2 = atom2.dipole.y;
    real dk3 = atom2.dipole.z;

    real qk1 = atom2.quadrupoleXX;
    real qk2 = atom2.quadrupoleXY;
    real qk3 = atom2.quadrupoleXZ;
    real qk5 = atom2.quadrupoleYY;
    real qk6 = atom2.quadrupoleYZ;
    real qk9 = -(atom2.quadrupoleXX + atom2.quadrupoleYY);

    real bn1 = bn.x;
    real bn2 = bn.y;
    real bn3 = bn.z;
    real bn4 = bn.w;

#ifdef APPLY_SCALE
    real offset = 1-mScale;
    real rr3 = rr1*rr1*rr1;
    real gf4 = 2*(bn2 - 3*offset*rr3*rr1*rr1);
#else
    real gf4 = 2*bn2;
#endif
    real qidk1 = qi1*dk1 + qi2*dk2 + qi3*dk3;
    real qkdi1 = qk1*di1 + qk2*di2 + qk3*di3;
    real ftm21 = gf4*(qkdi1-qidk1);

    real qidk2 = qi2*dk1 + qi5*dk2 + qi6*dk3;
    real qkdi2 = qk2*di1 + qk5*di2 + qk6*di3;
    real ftm22 = gf4*(qkdi2-qidk2);

    real qidk3 = qi3*dk1 + qi6*dk2 + qi9*dk3;
    real qkdi3 = qk3*di1 + qk6*di2 + qk9*di3;
    real ftm23 = gf4*(qkdi3-qidk3);

    real qir1 = qi1*xr + qi2*yr + qi3*zr;
    real qir2 = qi2*xr + qi5*yr + qi6*zr;
    real qir3 = qi3*xr + qi6*yr + qi9*zr;

    real qkr1 = qk1*xr + qk2*yr + qk3*zr;
    real qkr2 = qk2*xr + qk5*yr + qk6*zr;
    real qkr3 = qk3*xr + qk6*yr + qk9*zr;

#ifdef APPLY_SCALE
    real gf7 = 4*(bn3 - 15*offset*rr3*rr3*rr1);
#else
    real gf7 = 4*bn3;
#endif
    real qiqkr1 = qi1*qkr1 + qi2*qkr2 + qi3*qkr3;
    real qkqir1 = qk1*qir1 + qk2*qir2 + qk3*qir3;
    ftm21 += gf7*(qiqkr1+qkqir1);

    real qiqkr2 = qi2*qkr1 + qi5*qkr2 + qi6*qkr3;
    real qkqir2 = qk2*qir1 + qk5*qir2 + qk6*qir3;
    ftm22 += gf7*(qiqkr2+qkqir2);

    real qiqkr3 = qi3*qkr1 + qi6*qkr2 + qi9*qkr3;
    real qkqir3 = qk3*qir1 + qk6*qir2 + qk9*qir3;
    ftm23 += gf7*(qiqkr3+qkqir3);

    // calculate the scalar products for permanent components

    real gl6 = di1*dk1 + di2*dk2 + di3*dk3;
    real gl7 =  2*(qir1*dk1 + qir2*dk2 + qir3*dk3 - (qkr1*di1 + qkr2*di2 + qkr3*di3));
    real gl5 = -4*(qir1*qkr1 + qir2*qkr2 + qir3*qkr3);

    real gl8 =  2*(qi1*qk1 + qi2*qk2 + qi3*qk3 + qi2*qk2 + qi5*qk5 + qi6*qk6 + qi3*qk3 + qi6*qk6 + qi9*qk9);

    real sc3 = di1*xr + di2*yr + di3*zr;
    real sc5 = qir1*xr + qir2*yr + qir3*zr;
    real sc4 = dk1*xr + dk2*yr + dk3*zr;
    real sc6 = qkr1*xr + qkr2*yr + qkr3*zr;

    real gl0 = ci*ck;
    real gl1 = ck*sc3 - ci*sc4;
    real gl2 = ci*sc6 + ck*sc5 - sc3*sc4;
    real gl3 = sc3*sc6 - sc4*sc5;
    real gl4 = sc5*sc6;

#ifdef APPLY_SCALE
    energy += forceFactor*(-offset*rr1*gl0 + (bn1-offset*rr3)*(gl1+gl6) + (bn2-offset*(3*rr3*rr1*rr1))*(gl2+gl7+gl8) + (bn3-offset*(15*rr3*rr3*rr1))*(gl3+gl5) + (bn4-offset*(105*rr3*rr3*rr3))*gl4);
#else
    energy += forceFactor*(bn1*(gl1+gl6) + bn2*(gl2+gl7+gl8) + bn3*(gl3+gl5) + bn4*gl4);
    
#endif

    real gf1 = bn1*gl0 + bn2*(gl1+gl6) + bn3*(gl2+gl7+gl8) + bn4*(gl3+gl5) + bn5*gl4;
#ifdef APPLY_SCALE
    gf1 -= offset*(rr3*gl0 + (3*rr3*rr1*rr1)*(gl1+gl6) + (15*rr3*rr3*rr1)*(gl2+gl7+gl8) + (105*rr3*rr3*rr3)*(gl3+gl5) + (945*rr3*rr3*rr3*rr1*rr1)*gl4);
#endif
    ftm21 += gf1*xr;
    ftm22 += gf1*yr;
    ftm23 += gf1*zr;

#ifdef APPLY_SCALE
    real gf2 = -ck*bn1 + sc4*bn2 - sc6*bn3 - offset*(-ck*rr3 + sc4*(3*rr3*rr1*rr1) - sc6*(15*rr3*rr3*rr1));
#else
    real gf2 = -ck*bn1 + sc4*bn2 - sc6*bn3;
#endif
    ftm21 += gf2*di1;
    ftm22 += gf2*di2;
    ftm23 += gf2*di3;

#ifdef APPLY_SCALE
    real gf5 = 2*(-ck*bn2+sc4*bn3-sc6*bn4 - offset*(-ck*(3*rr3*rr1*rr1)+sc4*(15*rr3*rr3*rr1)-sc6*(105*rr3*rr3*rr3)));
#else
    real gf5 = 2*(-ck*bn2+sc4*bn3-sc6*bn4);
#endif
    ftm21 += gf5*qir1;
    ftm22 += gf5*qir2;
    ftm23 += gf5*qir3;

#ifdef APPLY_SCALE
    real gf3 = ci*bn1 + sc3*bn2 + sc5*bn3 - offset*(ci*rr3 + sc3*(3*rr3*rr1*rr1) + sc5*(15*rr3*rr3*rr1));
#else
    real gf3 = ci*bn1 + sc3*bn2 + sc5*bn3;
#endif
    ftm21 += gf3*dk1;
    ftm22 += gf3*dk2;
    ftm23 += gf3*dk3;

#ifdef APPLY_SCALE
    real gf6 = 2*(-ci*bn2-sc3*bn3-sc5*bn4 - offset*(-ci*(3*rr3*rr1*rr1)-sc3*(15*rr3*rr3*rr1)-sc5*(105*rr3*rr3*rr3)));
#else
    real gf6 = 2*(-ci*bn2-sc3*bn3-sc5*bn4);
#endif

    ftm21 += gf6*qkr1;
    ftm22 += gf6*qkr2;
    ftm23 += gf6*qkr3;

    force.x = ftm21;
    force.y = ftm22;
    force.z = ftm23;
}


__device__ void
#ifdef APPLY_SCALE
computeOneInteractionF2(
#else
computeOneInteractionF2NoScale(
#endif
        AtomData& atom1, volatile AtomData& atom2, real4 delta, real4 bn, float forceFactor,
#ifdef APPLY_SCALE
        float dScale, float pScale, float mScale,
#endif
        real3& force, real& energy) {
    const float uScale = 1;
    real xr = delta.x;
    real yr = delta.y;
    real zr = delta.z;
    real rr1 = delta.w;

    // set the permanent multipole and induced dipole values;

    real ci = atom1.q;

    real di1 = atom1.dipole.x;
    real di2 = atom1.dipole.y;
    real di3 = atom1.dipole.z;

    real qi1 = atom1.quadrupoleXX;
    real qi2 = atom1.quadrupoleXY;
    real qi3 = atom1.quadrupoleXZ;
    real qi5 = atom1.quadrupoleYY;
    real qi6 = atom1.quadrupoleYZ;
    real qi9 = -(atom1.quadrupoleXX + atom1.quadrupoleYY);

    real bn1 = bn.x;
    real bn2 = bn.y;
    real bn3 = bn.z;
    real bn4 = bn.w;

    real damp = atom1.damp*atom2.damp;
    if (damp != 0) {
        real pgamma = atom1.thole < atom2.thole ? atom1.thole : atom2.thole;
        real ratio = RECIP(rr1*damp);
        damp = -pgamma*ratio*ratio*ratio;
    }

    real scale5 = (damp == 0) ? 1 : (1 - (1-damp)*EXP(damp));
    real rr5 = rr1*rr1;
          rr5 = 3*rr1*rr5*rr5;
#ifdef APPLY_SCALE
    real psc5 = rr5*(1 - scale5*pScale);
    real dsc5 = rr5*(1 - scale5*dScale);
    real usc5 = rr5*(1 - scale5*uScale);
#else
    real psc5 = rr5*(1 - scale5);
#endif

    real qiuk1 = qi1*atom2.inducedDipole.x + qi2*atom2.inducedDipole.y + qi3*atom2.inducedDipole.z;
    real qiukp1 = qi1*atom2.inducedDipolePolar.x + qi2*atom2.inducedDipolePolar.y + qi3*atom2.inducedDipolePolar.z;
    real ftm21 = -bn2*(qiuk1+qiukp1);
#ifdef APPLY_SCALE
          ftm21 += qiuk1*psc5 + qiukp1*dsc5;
#else
          ftm21 += (qiuk1 + qiukp1)*psc5;
#endif

    real qiuk2 = qi2*atom2.inducedDipole.x + qi5*atom2.inducedDipole.y + qi6*atom2.inducedDipole.z;
    real qiukp2 = qi2*atom2.inducedDipolePolar.x + qi5*atom2.inducedDipolePolar.y + qi6*atom2.inducedDipolePolar.z;
    real ftm22 = -bn2*(qiuk2+qiukp2);
#ifdef APPLY_SCALE
          ftm22 += ((qiuk2)*psc5 + (qiukp2)*dsc5);
#else
          ftm22 += (qiuk2 + qiukp2)*psc5;
#endif

    real qiuk3 = qi3*atom2.inducedDipole.x + qi6*atom2.inducedDipole.y + qi9*atom2.inducedDipole.z;
    real qiukp3 = qi3*atom2.inducedDipolePolar.x + qi6*atom2.inducedDipolePolar.y + qi9*atom2.inducedDipolePolar.z;
    real ftm23 = -bn2*(qiuk3+qiukp3);
#ifdef APPLY_SCALE
          ftm23 += ((qiuk3)*psc5 + (qiukp3)*dsc5);
#else
          ftm23 += (qiuk3 + qiukp3)*psc5;
#endif

    real expdamp = EXP(damp);
    real scale3 = (damp == 0) ? 1 : (1 - expdamp);
    real rr3 = rr1*rr1*rr1;

#ifdef APPLY_SCALE
    real psc3 = rr3*(1 - scale3*pScale);
    real dsc3 = rr3*(1 - scale3*dScale);
    real usc3 = rr3*(1 - scale3*uScale);
#else
    real psc3 = rr3*(1 - scale3);
#endif

    real scale7 = (damp == 0) ? 1 : (1 - (1-damp+0.6f*damp*damp)*expdamp);

#ifdef APPLY_SCALE
    real psc7 = (15*rr3*rr3*rr1)*(1 - scale7*pScale);
    real dsc7 = (15*rr3*rr3*rr1)*(1 - scale7*dScale);
#else
    real psc7 = (15*rr3*rr3*rr1)*(1 - scale7);
#endif

    real qir1 = qi1*xr + qi2*yr + qi3*zr;
    real qir2 = qi2*xr + qi5*yr + qi6*zr;
    real qir3 = qi3*xr + qi6*yr + qi9*zr;

    real sc3 = di1*xr + di2*yr + di3*zr;
    real sc5 = qir1*xr + qir2*yr + qir3*zr;
    real gfi3 = ci*bn1 + sc3*bn2 + sc5*bn3;

    real prefactor1;
    prefactor1 = 0.5f*(ci*psc3 + sc3*psc5 + sc5*psc7 - gfi3);
    ftm21 -= prefactor1*atom2.inducedDipole.x;
    ftm22 -= prefactor1*atom2.inducedDipole.y;
    ftm23 -= prefactor1*atom2.inducedDipole.z;

#ifdef APPLY_SCALE
    prefactor1 = 0.5f*(ci*dsc3 + sc3*dsc5 + sc5*dsc7 - gfi3);
#endif
    ftm21 -= prefactor1*atom2.inducedDipolePolar.x;
    ftm22 -= prefactor1*atom2.inducedDipolePolar.y;
    ftm23 -= prefactor1*atom2.inducedDipolePolar.z;

    real sci4 = atom2.inducedDipole.x*xr + atom2.inducedDipole.y*yr + atom2.inducedDipole.z*zr;
    energy += forceFactor*0.5f*sci4*((psc3-bn1)*ci + (psc5-bn2)*sc3 + (psc7-bn3)*sc5);

    real scip4 = atom2.inducedDipolePolar.x*xr + atom2.inducedDipolePolar.y*yr + atom2.inducedDipolePolar.z*zr;
#ifndef DIRECT_POLARIZATION
#ifdef APPLY_SCALE
    prefactor1 = 0.5f*(bn2 - usc5);
#else
    prefactor1 = 0.5f*(bn2 - psc5);
#endif
    ftm21 += prefactor1*((sci4*atom1.inducedDipolePolar.x + scip4*atom1.inducedDipole.x));
    ftm22 += prefactor1*((sci4*atom1.inducedDipolePolar.y + scip4*atom1.inducedDipole.y));
    ftm23 += prefactor1*((sci4*atom1.inducedDipolePolar.z + scip4*atom1.inducedDipole.z));
#endif

#ifdef APPLY_SCALE
    prefactor1 = 0.5f*(bn2*(sci4+scip4) - (sci4*psc5+scip4*dsc5)); 
#else
    sci4 += scip4;
    prefactor1 = 0.5f*sci4*(bn2 - psc5); 
#endif

    ftm21 += prefactor1*di1;
    ftm22 += prefactor1*di2;
    ftm23 += prefactor1*di3;

#ifdef APPLY_SCALE
    real gfi5 = bn3*(sci4+scip4) - (sci4*psc7+scip4*dsc7);
#else
    real gfi5 = sci4*(bn3 - psc7);
#endif
    ftm21 += gfi5*qir1;
    ftm22 += gfi5*qir2;
    ftm23 += gfi5*qir3;

    real sci7 = qir1*atom2.inducedDipole.x + qir2*atom2.inducedDipole.y + qir3*atom2.inducedDipole.z;
    energy += forceFactor*(bn2-psc5)*sci7;
    real scip7 = qir1*atom2.inducedDipolePolar.x + qir2*atom2.inducedDipolePolar.y + qir3*atom2.inducedDipolePolar.z;

#ifdef APPLY_SCALE
    real gli1 = -ci*sci4;
    real gli2 = -sc3*sci4 + 2*sci7;
    real gli3 = -sci4*sc5;
    real glip1 = -ci*scip4;
    real glip2 = -sc3*scip4 + 2*scip7;
    real glip3 = -scip4*sc5;
#else
    real gli1 = -ci*sci4;
    real gli2 = -sc3*sci4 + 2*(sci7 + scip7);
    real gli3 = -sci4*sc5;
#endif

#ifdef APPLY_SCALE
    real gfi1 = (bn2*(gli1+glip1) + bn3*(gli2+glip2) + bn4*(gli3+glip3));
    gfi1 -= (rr1*rr1)*(3*(gli1*psc3 + glip1*dsc3) + 5*(gli2*psc5 + glip2*dsc5) + 7*(gli3*psc7+glip3*dsc7));
#else
    real gfi1 = bn2*gli1 + bn3*gli2 + bn4*gli3;
    gfi1 -= (rr1*rr1)*(3*gli1*psc3 + 5*gli2*psc5 + 7*gli3*psc7);
#endif
    gfi1 *= 0.5f;
    ftm21 += gfi1*xr;
    ftm22 += gfi1*yr;
    ftm23 += gfi1*zr;

    {
        real expdamp = EXP(damp);
        real temp3 = -1.5f*damp*expdamp*rr1*rr1;
        real temp5 = -damp;
        real temp7 = -0.2f - 0.6f*damp;

        real ddsc31 = temp3*xr;
        real ddsc32 = temp3*yr;
        real ddsc33 = temp3*zr;

        real ddsc51 = temp5*ddsc31;
        real ddsc52 = temp5*ddsc32;
        real ddsc53 = temp5*ddsc33;

        real ddsc71 = temp7*ddsc51;
        real ddsc72 = temp7*ddsc52;
        real ddsc73 = temp7*ddsc53;

        real rr3 = rr1*rr1*rr1;
#ifdef APPLY_SCALE
        temp3 = (gli1*pScale + glip1*dScale);
        temp5 = (3*rr1*rr1)*(gli2*pScale + glip2*dScale);
        temp7 = (15*rr3*rr1)*(gli3*pScale + glip3*dScale);
#else
        temp3 = gli1;
        temp5 = (3*rr1*rr1)*gli2;
        temp7 = (15*rr3*rr1)*gli3;
#endif
        ftm21 -= rr3*(temp3*ddsc31 + temp5*ddsc51 + temp7*ddsc71);
        ftm22 -= rr3*(temp3*ddsc32 + temp5*ddsc52 + temp7*ddsc72);
        ftm23 -= rr3*(temp3*ddsc33 + temp5*ddsc53 + temp7*ddsc73);
    }

//K
    real qk1 = atom2.quadrupoleXX;
    real qk2 = atom2.quadrupoleXY;
    real qk3 = atom2.quadrupoleXZ;
    real qk5 = atom2.quadrupoleYY;
    real qk6 = atom2.quadrupoleYZ;
    real qk9 = -(qk1 + qk5);

    real qkui1 = qk1*atom1.inducedDipole.x + qk2*atom1.inducedDipole.y + qk3*atom1.inducedDipole.z;
    real qkuip1 = qk1*atom1.inducedDipolePolar.x + qk2*atom1.inducedDipolePolar.y + qk3*atom1.inducedDipolePolar.z;
          ftm21 += bn2*(qkui1+qkuip1);
#ifdef APPLY_SCALE
          ftm21 -= (qkui1*psc5 + qkuip1*dsc5);
#else
          ftm21 -= (qkui1 + qkuip1)*psc5;
#endif

    real qkui2 = qk2*atom1.inducedDipole.x + qk5*atom1.inducedDipole.y + qk6*atom1.inducedDipole.z;
    real qkuip2 = qk2*atom1.inducedDipolePolar.x + qk5*atom1.inducedDipolePolar.y + qk6*atom1.inducedDipolePolar.z;
          ftm22 += bn2*(qkui2+qkuip2);
#ifdef APPLY_SCALE
          ftm22 -= ((qkui2)*psc5 + (qkuip2)*dsc5);
#else
          ftm22 -= (qkui2 + qkuip2)*psc5;
#endif

    real qkui3 = qk3*atom1.inducedDipole.x + qk6*atom1.inducedDipole.y + qk9*atom1.inducedDipole.z;
    real qkuip3 = qk3*atom1.inducedDipolePolar.x + qk6*atom1.inducedDipolePolar.y + qk9*atom1.inducedDipolePolar.z;
          ftm23 += bn2*(qkui3+qkuip3);
#ifdef APPLY_SCALE
          ftm23 -= ((qkui3)*psc5 + (qkuip3)*dsc5);
#else
          ftm23 -= (qkui3 + qkuip3)*psc5;
#endif


    real qkr1 = qk1*xr + qk2*yr + qk3*zr;
    real qkr2 = qk2*xr + qk5*yr + qk6*zr;
    real qkr3 = qk3*xr + qk6*yr + qk9*zr;

    real dk1 = atom2.dipole.x;
    real dk2 = atom2.dipole.y;
    real dk3 = atom2.dipole.z;

    real sc4 =  dk1*xr +  dk2*yr +  dk3*zr;
    real sc6 = qkr1*xr + qkr2*yr + qkr3*zr;

    real ck = atom2.q;
    real gfi2 = (-ck*bn1 + sc4*bn2 - sc6*bn3);

    prefactor1 = 0.5f*(ck*psc3 - sc4*psc5 + sc6*psc7 + gfi2);
    ftm21 += prefactor1*atom1.inducedDipole.x;
    ftm22 += prefactor1*atom1.inducedDipole.y;
    ftm23 += prefactor1*atom1.inducedDipole.z;

#ifdef APPLY_SCALE
    prefactor1 = 0.5f*(ck*dsc3 - sc4*dsc5 + sc6*dsc7 + gfi2);
#endif
    ftm21 += prefactor1*atom1.inducedDipolePolar.x;
    ftm22 += prefactor1*atom1.inducedDipolePolar.y;
    ftm23 += prefactor1*atom1.inducedDipolePolar.z;

    real sci3 = atom1.inducedDipole.x*xr + atom1.inducedDipole.y*yr + atom1.inducedDipole.z*zr;
    energy += forceFactor*0.5f*sci3*(ck*(bn1-psc3) - sc4*(bn2-psc5) + sc6*(bn3-psc7));
    real scip3 = atom1.inducedDipolePolar.x*xr + atom1.inducedDipolePolar.y*yr + atom1.inducedDipolePolar.z*zr;

#ifndef DIRECT_POLARIZATION
#ifdef APPLY_SCALE
    prefactor1 = 0.5f*(bn2 - usc5);
#else
    prefactor1 = 0.5f*(bn2 - psc5);
#endif

    ftm21 += prefactor1*(sci3*atom2.inducedDipolePolar.x + scip3*atom2.inducedDipole.x);
    ftm22 += prefactor1*(sci3*atom2.inducedDipolePolar.y + scip3*atom2.inducedDipole.y);
    ftm23 += prefactor1*(sci3*atom2.inducedDipolePolar.z + scip3*atom2.inducedDipole.z);
    
    real sci34;
    sci4 = atom2.inducedDipole.x*xr + atom2.inducedDipole.y*yr + atom2.inducedDipole.z*zr;
    scip4 = atom2.inducedDipolePolar.x*xr + atom2.inducedDipolePolar.y*yr + atom2.inducedDipolePolar.z*zr;
    sci34 = (sci3*scip4+scip3*sci4);

#ifdef APPLY_SCALE
    gfi1 = sci34*(usc5*(5*rr1*rr1) -bn3);
#else
    gfi1 = sci34*(psc5*(5*rr1*rr1) -bn3);
#endif
#else
    gfi1 = 0;
#endif
    
#ifdef APPLY_SCALE
    prefactor1 = 0.5f*(bn2*(sci3+scip3) - (sci3*psc5+scip3*dsc5));
#else
    sci3 += scip3;
    prefactor1 = 0.5f*sci3*(bn2 - psc5);
#endif
    ftm21 += prefactor1*dk1;
    ftm22 += prefactor1*dk2;
    ftm23 += prefactor1*dk3;

#ifdef APPLY_SCALE
    real gfi6 = -bn3*(sci3+scip3) + (sci3*psc7+scip3*dsc7);
#else
    real gfi6 = sci3*(psc7 - bn3);
#endif
    ftm21 += gfi6*qkr1;
    ftm22 += gfi6*qkr2;
    ftm23 += gfi6*qkr3;

    real sci1 = atom1.inducedDipole.x*dk1 + atom1.inducedDipole.y*dk2 + atom1.inducedDipole.z*dk3 + di1*atom2.inducedDipole.x + di2*atom2.inducedDipole.y + di3*atom2.inducedDipole.z;
    energy += forceFactor*0.5f*(sci1*(bn1-psc3));

    real sci8 = qkr1*atom1.inducedDipole.x + qkr2*atom1.inducedDipole.y + qkr3*atom1.inducedDipole.z;
    energy -= forceFactor*sci8*(bn2-psc5);
    real scip1 = atom1.inducedDipolePolar.x*dk1 + atom1.inducedDipolePolar.y*dk2 + atom1.inducedDipolePolar.z*dk3 + di1*atom2.inducedDipolePolar.x + di2*atom2.inducedDipolePolar.y + di3*atom2.inducedDipolePolar.z;
#ifndef APPLY_SCALE
        sci1 += scip1;
#endif

    real scip2 = atom1.inducedDipole.x*atom2.inducedDipolePolar.x +
                                  atom1.inducedDipole.y*atom2.inducedDipolePolar.y +
                                  atom1.inducedDipole.z*atom2.inducedDipolePolar.z +
                                  atom2.inducedDipole.x*atom1.inducedDipolePolar.x +
                                  atom2.inducedDipole.y*atom1.inducedDipolePolar.y +
                                  atom2.inducedDipole.z*atom1.inducedDipolePolar.z;

    real scip8 = qkr1*atom1.inducedDipolePolar.x + qkr2*atom1.inducedDipolePolar.y + qkr3*atom1.inducedDipolePolar.z;
#ifndef APPLY_SCALE
          sci8 += scip8;
#endif

           gli1 = ck*sci3 + sci1;
           gli2 = -(sci3*sc4 + 2*sci8);
           gli3 = sci3*sc6;
#ifdef APPLY_SCALE
          glip1 = ck*scip3 + scip1;
          glip2 = -(scip3*sc4 + 2*scip8);
          glip3 = scip3*sc6;
#endif


#ifdef APPLY_SCALE
    gfi1 += (bn2*(gli1+glip1) + bn3*(gli2+glip2) + bn4*(gli3+glip3));
    gfi1 -= (rr1*rr1)*(3*(gli1*psc3 + glip1*dsc3) + 5*(gli2*psc5 + glip2*dsc5) + 7*(gli3*psc7+glip3*dsc7));
#else
    gfi1 += (bn2*gli1 + bn3*gli2 + bn4*gli3);
    gfi1 -= (rr1*rr1)*(3*gli1*psc3 + 5*gli2*psc5 + 7*gli3*psc7);
#endif
    
#ifndef DIRECT_POLARIZATION
#ifdef APPLY_SCALE
    gfi1 += scip2*(bn2 - (3*rr1*rr1)*usc3);
#else
    gfi1 += scip2*(bn2 - (3*rr1*rr1)*psc3);
#endif
#endif
    
    gfi1 *= 0.5f;

    ftm21 += gfi1*xr;
    ftm22 += gfi1*yr;
    ftm23 += gfi1*zr;

    {
        real expdamp = EXP(damp);
        real temp3 = -1.5f*damp*expdamp*rr1*rr1;
        real temp5 = -damp;
        real temp7 = -0.2f - 0.6f*damp;

        real ddsc31 = temp3*xr;
        real ddsc32 = temp3*yr;
        real ddsc33 = temp3*zr;

        real ddsc51 = temp5*ddsc31;
        real ddsc52 = temp5*ddsc32;
        real ddsc53 = temp5*ddsc33;

        real ddsc71 = temp7*ddsc51;
        real ddsc72 = temp7*ddsc52;
        real ddsc73 = temp7*ddsc53;

        real rr3 = rr1*rr1*rr1;

#ifdef APPLY_SCALE
        temp3 = gli1*pScale + glip1*dScale;
        temp5 = (3*rr1*rr1)*(gli2*pScale + glip2*dScale);
        temp7 = (15*rr3*rr1)*(gli3*pScale + glip3*dScale);
#else
        temp3 = gli1;
        temp5 = (3*rr1*rr1)*gli2;
        temp7 = (15*rr3*rr1)*(gli3);
#endif

        ftm21 -= rr3*(temp3*ddsc31 + temp5*ddsc51 + temp7*ddsc71);
        ftm22 -= rr3*(temp3*ddsc32 + temp5*ddsc52 + temp7*ddsc72);
        ftm23 -= rr3*(temp3*ddsc33 + temp5*ddsc53 + temp7*ddsc73);

#ifndef DIRECT_POLARIZATION
#ifdef APPLY_SCALE
        temp3 =  uScale*scip2;
        temp5 = -(3*rr1*rr1)*uScale*sci34;
#else
        temp3 =  scip2;
        temp5 = -(3*rr1*rr1)*sci34;
#endif
        ftm21 -= rr3*(temp3*ddsc31 + temp5*ddsc51);
        ftm22 -= rr3*(temp3*ddsc32 + temp5*ddsc52);
        ftm23 -= rr3*(temp3*ddsc33 + temp5*ddsc53);
#endif
    }

    force.x += ftm21;
    force.y += ftm22;
    force.z += ftm23;
}


__device__ void
#ifdef APPLY_SCALE
computeOneInteractionT1(
#else
computeOneInteractionT1NoScale(
#endif
        AtomData& atom1, volatile AtomData& atom2, const real4 delta, const real4 bn
#ifdef APPLY_SCALE
        , float dScale, float pScale, float mScale
#endif
        ) {

    real xr = delta.x;
    real yr = delta.y;
    real zr = delta.z;
#ifdef APPLY_SCALE
    real rr1 = delta.w;
#endif

    // set the permanent multipole and induced dipole values;

    real di1 = atom1.dipole.x;
    real di2 = atom1.dipole.y;
    real di3 = atom1.dipole.z;

    real qi1 = atom1.quadrupoleXX;
    real qi2 = atom1.quadrupoleXY;
    real qi3 = atom1.quadrupoleXZ;
    real qi5 = atom1.quadrupoleYY;
    real qi6 = atom1.quadrupoleYZ;
    //real qi9 = atom1.labFrameQuadrupole[5];
    real qi9 = -(atom1.quadrupoleXX + atom1.quadrupoleYY);

    real ck = atom2.q;

    real dk1 = atom2.dipole.x;
    real dk2 = atom2.dipole.y;
    real dk3 = atom2.dipole.z;

    real qk1 = atom2.quadrupoleXX;
    real qk2 = atom2.quadrupoleXY;
    real qk3 = atom2.quadrupoleXZ;
    real qk5 = atom2.quadrupoleYY;
    real qk6 = atom2.quadrupoleYZ;
    //real qk9 = atom2.labFrameQuadrupole[5];
    real qk9 = -(atom2.quadrupoleXX + atom2.quadrupoleYY);

    real bn1 = bn.x;
    real bn2 = bn.y;
    real bn3 = bn.z;
    real bn4 = bn.w;

    // apply Thole polarization damping to scale factors

#ifdef APPLY_SCALE
    real rr2 = rr1*rr1;
    real rr3 = rr1*rr2;
    real rr5 = 3*rr3*rr2;
    real rr7 = 5*rr5*rr2;
    real rr9 = 7*rr7*rr2;

    real scale = 1-mScale;
    real prefactor = scale*rr3 - bn1;
#else
    real prefactor = -bn1;
#endif
    real dixdk1 = di2*dk3 - di3*dk2;
    real ttm21 = prefactor*dixdk1;

    real dixdk2 = di3*dk1 - di1*dk3;
    real ttm22 = prefactor*dixdk2;

    real dixdk3 = di1*dk2 - di2*dk1;
    real ttm23 = prefactor*dixdk3;

    real qir1 = qi1*xr + qi2*yr + qi3*zr;
    real qir2 = qi2*xr + qi5*yr + qi6*zr;
    real qir3 = qi3*xr + qi6*yr + qi9*zr;

    real qkr1 = qk1*xr + qk2*yr + qk3*zr;
    real qkr2 = qk2*xr + qk5*yr + qk6*zr;
    real qkr3 = qk3*xr + qk6*yr + qk9*zr;

    real qiqkr1 = qi1*qkr1 + qi2*qkr2 + qi3*qkr3;
    real qiqkr2 = qi2*qkr1 + qi5*qkr2 + qi6*qkr3;
    real qiqkr3 = qi3*qkr1 + qi6*qkr2 + qi9*qkr3;

    real rxqikr1 = yr*qiqkr3 - zr*qiqkr2;
    real qkrxqir1 = qkr2*qir3 - qkr3*qir2;
#ifdef APPLY_SCALE
    prefactor = 4*(bn3 - scale*rr7);
#else
    prefactor = 4*bn3;
#endif
    ttm21 -= prefactor*(rxqikr1+qkrxqir1);

    real rxqikr2 = zr*qiqkr1 - xr*qiqkr3;
    real qkrxqir2 = qkr3*qir1 - qkr1*qir3;
    ttm22 -= prefactor*(rxqikr2+qkrxqir2);

    real rxqikr3 = xr*qiqkr2 - yr*qiqkr1;
    real qkrxqir3 = qkr1*qir2 - qkr2*qir1;
    ttm23 -= prefactor*(rxqikr3+qkrxqir3);

    real qidk1 = qi1*dk1 + qi2*dk2 + qi3*dk3;
    real qidk2 = qi2*dk1 + qi5*dk2 + qi6*dk3;
    real qidk3 = qi3*dk1 + qi6*dk2 + qi9*dk3;

    real dixqkr1 = di2*qkr3 - di3*qkr2;
    real dkxqir1 = dk2*qir3 - dk3*qir2;
    real rxqidk1 = yr*qidk3 - zr*qidk2;
    real qixqk1 = qi2*qk3 + qi5*qk6 + qi6*qk9 - qi3*qk2 - qi6*qk5 - qi9*qk6;
#ifdef APPLY_SCALE
    prefactor = 2*(bn2 - scale*rr5);
#else
    prefactor = 2*bn2;
#endif
    ttm21 += prefactor*(dixqkr1+dkxqir1+rxqidk1-2*qixqk1);
 
    real dixqkr2 = di3*qkr1 - di1*qkr3;
    real dkxqir2 = dk3*qir1 - dk1*qir3;
    real rxqidk2 = zr*qidk1 - xr*qidk3;
    real qixqk2 = qi3*qk1 + qi6*qk2 + qi9*qk3 - qi1*qk3 - qi2*qk6 - qi3*qk9;
    ttm22 += prefactor*(dixqkr2+dkxqir2+rxqidk2-2*qixqk2);

    real dixqkr3 = di1*qkr2 - di2*qkr1;
    real dkxqir3 = dk1*qir2 - dk2*qir1;
    real rxqidk3 = xr*qidk2 - yr*qidk1;
    real qixqk3 = qi1*qk2 + qi2*qk5 + qi3*qk6 - qi2*qk1 - qi5*qk2 - qi6*qk3;
    ttm23 += prefactor*(dixqkr3+dkxqir3+rxqidk3-2*qixqk3);

    real sc4 = dk1*xr + dk2*yr + dk3*zr;
    real sc6 = qkr1*xr + qkr2*yr + qkr3*zr;

    real gf2 = -ck*bn1 + sc4*bn2 - sc6*bn3;
#ifdef APPLY_SCALE
    real gfr2 = -ck*rr3 + sc4*rr5 - sc6*rr7;
    prefactor = (gf2 - scale*gfr2);
#else
    prefactor = gf2;
#endif
    ttm21 += prefactor*(di2*zr - di3*yr);
    ttm22 += prefactor*(di3*xr - di1*zr);
    ttm23 += prefactor*(di1*yr - di2*xr);

    real gf5 = (-ck*bn2+sc4*bn3-sc6*bn4);
#ifdef APPLY_SCALE
    real gfr5 = (-ck*rr5+sc4*rr7-sc6*rr9); 
    prefactor = 2*(gf5 - scale*gfr5);
#else
    prefactor = 2*gf5;
#endif

    real rxqir1 = yr*qir3 - zr*qir2;
    real rxqir2 = zr*qir1 - xr*qir3;
    real rxqir3 = xr*qir2 - yr*qir1;
    ttm21 -= prefactor*rxqir1; 
    ttm22 -= prefactor*rxqir2;
    ttm23 -= prefactor*rxqir3;

    atom1.torque.x += ttm21;
    atom1.torque.y += ttm22;
    atom1.torque.z += ttm23;
}


__device__ void
#ifdef APPLY_SCALE
computeOneInteractionT2(
#else
computeOneInteractionT2NoScale(
#endif
        AtomData& atom1, volatile AtomData& atom2, const real4 delta, const real4 bn
#ifdef APPLY_SCALE
        , float dScale, float pScale, float mScale
#endif
        ) {

    real xr = delta.x;
    real yr = delta.y;
    real zr = delta.z;
    real rr1 = delta.w;

    // set the permanent multipole and induced dipole values;

    real di1 = atom1.dipole.x;
    real di2 = atom1.dipole.y;
    real di3 = atom1.dipole.z;

    real qi1 = atom1.quadrupoleXX;
    real qi2 = atom1.quadrupoleXY;
    real qi3 = atom1.quadrupoleXZ;
    real qi5 = atom1.quadrupoleYY;
    real qi6 = atom1.quadrupoleYZ;
    real qi9 = -(atom1.quadrupoleXX + atom1.quadrupoleYY);

    real bn1 = bn.x;
    real bn2 = bn.y;
    real bn3 = bn.z;

    // apply Thole polarization damping to scale factors

    real scale3 = 1;
    real scale5 = 1;
    real scale7 = 1;

    real damp = atom1.damp*atom2.damp;
    if (damp != 0) {
        real pgamma = atom1.thole < atom2.thole ? atom1.thole : atom2.thole;
        real ratio = RECIP(rr1*damp);
        damp = -pgamma*ratio*ratio*ratio;
        real expdamp = EXP(damp);
        scale3 = 1 - expdamp;
        scale5 = 1 - (1-damp)*expdamp;
        scale7 = 1 - (1-damp+0.6f*damp*damp)*expdamp;
    }

    real rr3 = rr1*rr1*rr1;
#ifdef APPLY_SCALE
    real dsc3 = rr3*(1 - scale3*dScale);
    real dsc5 = (3*rr3*rr1*rr1)* (1 - scale5*dScale);
    real dsc7 = (15*rr3*rr3*rr1)*(1 - scale7*dScale);

    real psc3 = rr3*(1 - scale3*pScale);
    real psc5 = (3*rr3*rr1*rr1)*(1 - scale5*pScale);
    real psc7 = (15*rr3*rr3*rr1)*(1 - scale7*pScale);
#else
    real psc3 = rr3*(1 - scale3);
    real psc5 = (3*rr3*rr1*rr1)*(1 - scale5);
    real psc7 = (15*rr3*rr3*rr1)*(1 - scale7);
#endif

    real prefactor1 = 0.5f*(psc3 - bn1);
#ifdef APPLY_SCALE
    real prefactor2 = 0.5f*(dsc3 - bn1);
#endif

    real dixuk1 = di2*atom2.inducedDipole.z - di3*atom2.inducedDipole.y;
    real dixukp1 = di2*atom2.inducedDipolePolar.z - di3*atom2.inducedDipolePolar.y;

#ifdef APPLY_SCALE
    real ttm2i1 = prefactor1*dixuk1 + prefactor2*dixukp1;
#else
    real ttm2i1 = prefactor1*(dixuk1 + dixukp1);
#endif

    real dixuk2 = di3*atom2.inducedDipole.x - di1*atom2.inducedDipole.z;
    real dixukp2 = di3*atom2.inducedDipolePolar.x - di1*atom2.inducedDipolePolar.z;

#ifdef APPLY_SCALE
    real ttm2i2 = prefactor1*dixuk2 + prefactor2*dixukp2;
#else
    real ttm2i2 = prefactor1*(dixuk2 + dixukp2);
#endif

    real dixuk3 = di1*atom2.inducedDipole.y - di2*atom2.inducedDipole.x;
    real dixukp3 = di1*atom2.inducedDipolePolar.y - di2*atom2.inducedDipolePolar.x;
#ifdef APPLY_SCALE
    real ttm2i3 = prefactor1*dixuk3 + prefactor2*dixukp3;
#else
    real ttm2i3 = prefactor1*(dixuk3 + dixukp3);
#endif

    real sci4 = atom2.inducedDipole.x*xr + atom2.inducedDipole.y*yr + atom2.inducedDipole.z*zr;
    real scip4 = atom2.inducedDipolePolar.x*xr + atom2.inducedDipolePolar.y*yr + atom2.inducedDipolePolar.z*zr;
    real gti2 = bn2*(sci4+scip4);
#ifdef APPLY_SCALE
    real gtri2 = (sci4*psc5+scip4*dsc5);
#else
    real gtri2 = psc5*(sci4+scip4);
#endif
    prefactor1 = 0.5f*(gti2 - gtri2);

    ttm2i1 += prefactor1*(di2*zr - di3*yr);
    ttm2i2 += prefactor1*(di3*xr - di1*zr);
    ttm2i3 += prefactor1*(di1*yr - di2*xr);

    real qir1 = qi1*xr + qi2*yr + qi3*zr;
    real qir2 = qi2*xr + qi5*yr + qi6*zr;
    real qir3 = qi3*xr + qi6*yr + qi9*zr;

#ifdef APPLY_SCALE
    prefactor1 = sci4*psc7 + scip4*dsc7 - bn3*(sci4+scip4);
#else
    prefactor1 = psc7*(sci4+scip4) - bn3*(sci4+scip4);
#endif
    ttm2i1 += prefactor1*(yr*qir3 - zr*qir2);
    ttm2i2 += prefactor1*(zr*qir1 - xr*qir3);
    ttm2i3 += prefactor1*(xr*qir2 - yr*qir1);

    real qiuk1 = qi1*atom2.inducedDipole.x + qi2*atom2.inducedDipole.y + qi3*atom2.inducedDipole.z;
    real qiuk2 = qi2*atom2.inducedDipole.x + qi5*atom2.inducedDipole.y + qi6*atom2.inducedDipole.z;
    real qiuk3 = qi3*atom2.inducedDipole.x + qi6*atom2.inducedDipole.y + qi9*atom2.inducedDipole.z;

    real qiukp1 = qi1*atom2.inducedDipolePolar.x + qi2*atom2.inducedDipolePolar.y + qi3*atom2.inducedDipolePolar.z;
    real qiukp2 = qi2*atom2.inducedDipolePolar.x + qi5*atom2.inducedDipolePolar.y + qi6*atom2.inducedDipolePolar.z;
    real qiukp3 = qi3*atom2.inducedDipolePolar.x + qi6*atom2.inducedDipolePolar.y + qi9*atom2.inducedDipolePolar.z;

    prefactor1 = (bn2 - psc5);
#ifdef APPLY_SCALE
    prefactor2 = (bn2 - dsc5);
#endif
    real ukxqir1 = atom2.inducedDipole.y*qir3 - atom2.inducedDipole.z*qir2;
    real ukxqirp1 = atom2.inducedDipolePolar.y*qir3 - atom2.inducedDipolePolar.z*qir2;
    real rxqiuk1 = yr*qiuk3 - zr*qiuk2;
    real rxqiukp1 = yr*qiukp3 - zr*qiukp2;

#ifdef APPLY_SCALE
    ttm2i1 += prefactor1*(ukxqir1 + rxqiuk1) + prefactor2*(ukxqirp1 + rxqiukp1);
#else
    ttm2i1 += prefactor1*(ukxqir1 + rxqiuk1 + ukxqirp1 + rxqiukp1);
#endif

    real ukxqir2 = atom2.inducedDipole.z*qir1 - atom2.inducedDipole.x*qir3;
    real ukxqirp2 = atom2.inducedDipolePolar.z*qir1 - atom2.inducedDipolePolar.x*qir3;
    real rxqiuk2 = zr*qiuk1 - xr*qiuk3;
    real rxqiukp2 = zr*qiukp1 - xr*qiukp3;
#ifdef APPLY_SCALE
    ttm2i2 += prefactor1*(ukxqir2 + rxqiuk2) + prefactor2*(ukxqirp2 + rxqiukp2);
#else
    ttm2i2 += prefactor1*(ukxqir2 + rxqiuk2 + ukxqirp2 + rxqiukp2);
#endif

    real ukxqir3 = atom2.inducedDipole.x*qir2 - atom2.inducedDipole.y*qir1;
    real ukxqirp3 = atom2.inducedDipolePolar.x*qir2 - atom2.inducedDipolePolar.y*qir1;
    real rxqiuk3 = xr*qiuk2 - yr*qiuk1;
    real rxqiukp3 = xr*qiukp2 - yr*qiukp1;
#ifdef APPLY_SCALE
    ttm2i3 += prefactor1*(ukxqir3 + rxqiuk3) + prefactor2*(ukxqirp3 + rxqiukp3);
#else
    ttm2i3 += prefactor1*(ukxqir3 + rxqiuk3 + ukxqirp3 + rxqiukp3);
#endif

    atom1.torque.x += ttm2i1;
    atom1.torque.y += ttm2i2;
    atom1.torque.z += ttm2i3;
}
