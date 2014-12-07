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

    real ck = atom2.q;
    real dk1 = atom2.dipole.x;
    real dk2 = atom2.dipole.y;
    real dk3 = atom2.dipole.z;

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
    real ftm21 = 0;
    real ftm22 = 0;
    real ftm23 = 0;

    // calculate the scalar products for permanent components

    real gl6 = di1*dk1 + di2*dk2 + di3*dk3;
    real sc3 = di1*xr + di2*yr + di3*zr;
    real sc4 = dk1*xr + dk2*yr + dk3*zr;

    real gl0 = ci*ck;
    real gl1 = ck*sc3 - ci*sc4;
    real gl2 = -sc3*sc4;

#ifdef APPLY_SCALE
    energy += forceFactor*(-offset*rr1*gl0 + (bn1-offset*rr3)*(gl1+gl6) + (bn2-offset*(3*rr3*rr1*rr1))*gl2);
#else
    energy += forceFactor*(bn1*(gl1+gl6) + bn2*gl2);
    
#endif

    real gf1 = bn1*gl0 + bn2*(gl1+gl6) + bn3*gl2;
#ifdef APPLY_SCALE
    gf1 -= offset*(rr3*gl0 + (3*rr3*rr1*rr1)*(gl1+gl6) + (15*rr3*rr3*rr1)*gl2);
#endif
    ftm21 += gf1*xr;
    ftm22 += gf1*yr;
    ftm23 += gf1*zr;

#ifdef APPLY_SCALE
    real gf2 = -ck*bn1 + sc4*bn2 - offset*(-ck*rr3 + sc4*(3*rr3*rr1*rr1));
#else
    real gf2 = -ck*bn1 + sc4*bn2;
#endif
    ftm21 += gf2*di1;
    ftm22 += gf2*di2;
    ftm23 += gf2*di3;

#ifdef APPLY_SCALE
    real gf3 = ci*bn1 + sc3*bn2 - offset*(ci*rr3 + sc3*(3*rr3*rr1*rr1));
#else
    real gf3 = ci*bn1 + sc3*bn2;
#endif
    ftm21 += gf3*dk1;
    ftm22 += gf3*dk2;
    ftm23 += gf3*dk3;

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

    real ftm21 = 0;
    real ftm22 = 0;
    real ftm23 = 0;

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

    real sc3 = di1*xr + di2*yr + di3*zr;
    real gfi3 = ci*bn1 + sc3*bn2;

    real prefactor1;
    prefactor1 = 0.5f*(ci*psc3 + sc3*psc5 - gfi3);
    ftm21 -= prefactor1*atom2.inducedDipole.x;
    ftm22 -= prefactor1*atom2.inducedDipole.y;
    ftm23 -= prefactor1*atom2.inducedDipole.z;

#ifdef APPLY_SCALE
    prefactor1 = 0.5f*(ci*dsc3 + sc3*dsc5 - gfi3);
#endif
    ftm21 -= prefactor1*atom2.inducedDipolePolar.x;
    ftm22 -= prefactor1*atom2.inducedDipolePolar.y;
    ftm23 -= prefactor1*atom2.inducedDipolePolar.z;

    real sci4 = atom2.inducedDipole.x*xr + atom2.inducedDipole.y*yr + atom2.inducedDipole.z*zr;
    energy += forceFactor*0.5f*sci4*((psc3-bn1)*ci + (psc5-bn2)*sc3);

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
    real gli1 = -ci*sci4;
    real gli2 = -sc3*sci4;
    real glip1 = -ci*scip4;
    real glip2 = -sc3*scip4;
#else
    real gli1 = -ci*sci4;
    real gli2 = -sc3*sci4;
#endif

#ifdef APPLY_SCALE
    real gfi1 = (bn2*(gli1+glip1) + bn3*(gli2+glip2));
    gfi1 -= (rr1*rr1)*(3*(gli1*psc3 + glip1*dsc3) + 5*(gli2*psc5 + glip2*dsc5));
#else
    real gfi1 = bn2*gli1 + bn3*gli2;
    gfi1 -= (rr1*rr1)*(3*gli1*psc3 + 5*gli2*psc5);
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
#else
        temp3 = gli1;
        temp5 = (3*rr1*rr1)*gli2;
#endif
        ftm21 -= rr3*(temp3*ddsc31 + temp5*ddsc51);
        ftm22 -= rr3*(temp3*ddsc32 + temp5*ddsc52);
        ftm23 -= rr3*(temp3*ddsc33 + temp5*ddsc53);
    }

//K
    real dk1 = atom2.dipole.x;
    real dk2 = atom2.dipole.y;
    real dk3 = atom2.dipole.z;

    real sc4 =  dk1*xr +  dk2*yr +  dk3*zr;

    real ck = atom2.q;
    real gfi2 = (-ck*bn1 + sc4*bn2);

    prefactor1 = 0.5f*(ck*psc3 - sc4*psc5 + gfi2);
    ftm21 += prefactor1*atom1.inducedDipole.x;
    ftm22 += prefactor1*atom1.inducedDipole.y;
    ftm23 += prefactor1*atom1.inducedDipole.z;

#ifdef APPLY_SCALE
    prefactor1 = 0.5f*(ck*dsc3 - sc4*dsc5 + gfi2);
#endif
    ftm21 += prefactor1*atom1.inducedDipolePolar.x;
    ftm22 += prefactor1*atom1.inducedDipolePolar.y;
    ftm23 += prefactor1*atom1.inducedDipolePolar.z;

    real sci3 = atom1.inducedDipole.x*xr + atom1.inducedDipole.y*yr + atom1.inducedDipole.z*zr;
    energy += forceFactor*0.5f*sci3*(ck*(bn1-psc3) - sc4*(bn2-psc5));
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

    real sci1 = atom1.inducedDipole.x*dk1 + atom1.inducedDipole.y*dk2 + atom1.inducedDipole.z*dk3 + di1*atom2.inducedDipole.x + di2*atom2.inducedDipole.y + di3*atom2.inducedDipole.z;
    energy += forceFactor*0.5f*(sci1*(bn1-psc3));

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

           gli1 = ck*sci3 + sci1;
           gli2 = -sci3*sc4;
#ifdef APPLY_SCALE
          glip1 = ck*scip3 + scip1;
          glip2 = -scip3*sc4;
#endif


#ifdef APPLY_SCALE
    gfi1 += (bn2*(gli1+glip1) + bn3*(gli2+glip2));
    gfi1 -= (rr1*rr1)*(3*(gli1*psc3 + glip1*dsc3) + 5*(gli2*psc5 + glip2*dsc5));
#else
    gfi1 += (bn2*gli1 + bn3*gli2);
    gfi1 -= (rr1*rr1)*(3*gli1*psc3 + 5*gli2*psc5);
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
#else
        temp3 = gli1;
        temp5 = (3*rr1*rr1)*gli2;
#endif

        ftm21 -= rr3*(temp3*ddsc31 + temp5*ddsc51);
        ftm22 -= rr3*(temp3*ddsc32 + temp5*ddsc52);
        ftm23 -= rr3*(temp3*ddsc33 + temp5*ddsc53);

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

    real ck = atom2.q;

    real dk1 = atom2.dipole.x;
    real dk2 = atom2.dipole.y;
    real dk3 = atom2.dipole.z;

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

    real sc4 = dk1*xr + dk2*yr + dk3*zr;
    real sc6 = 0;

    real gf2 = -ck*bn1 + sc4*bn2;
#ifdef APPLY_SCALE
    real gfr2 = -ck*rr3 + sc4*rr5;
    prefactor = (gf2 - scale*gfr2);
#else
    prefactor = gf2;
#endif
    ttm21 += prefactor*(di2*zr - di3*yr);
    ttm22 += prefactor*(di3*xr - di1*zr);
    ttm23 += prefactor*(di1*yr - di2*xr);

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

    atom1.torque.x += ttm2i1;
    atom1.torque.y += ttm2i2;
    atom1.torque.z += ttm2i3;
}
