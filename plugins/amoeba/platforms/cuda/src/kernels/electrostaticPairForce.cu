/**
 * This defines three different closely related functions, depending on which constant (F1, T1, or T3) is defined.
 */

#if defined F1
__device__ void computeOneInteractionF1(AtomData& atom1, volatile AtomData& atom2, float dScale, float pScale, float mScale, real& energy, real3& outputForce) {
#elif defined T1
__device__ void computeOneInteractionT1(AtomData& atom1, volatile AtomData& atom2, float dScale, float pScale, float mScale, real3& outputForce) {
#else
__device__ void computeOneInteractionT3(AtomData& atom1, volatile AtomData& atom2, float dScale, float pScale, float mScale, real3& outputForce) {
#endif
    
#ifdef F1
    const float uScale = 1;
    real ddsc3_0 = 0;
    real ddsc3_1 = 0;
    real ddsc3_2 = 0;

    real ddsc5_0 = 0;
    real ddsc5_1 = 0;
    real ddsc5_2 = 0;

    real ddsc7_0 = 0;
    real ddsc7_1 = 0;
    real ddsc7_2 = 0;
#endif

    real xr = atom2.posq.x - atom1.posq.x;
    real yr = atom2.posq.y - atom1.posq.y;
    real zr = atom2.posq.z - atom1.posq.z;
    
    real r2 = xr*xr + yr*yr + zr*zr;
    real r = SQRT(r2);
    real rr1 = RECIP(r);
    real rr2 = rr1*rr1;
    real rr3 = rr1*rr2;
    real rr5 = 3*rr3*rr2;
    real rr7 = 5*rr5*rr2;
    real rr9 = 7*rr7*rr2;
#ifdef F1
    real rr11 = 9*rr9*rr2;
#endif

    real scale3 = 1;
    real scale5 = 1;
    real scale7 = 1;

    real pdamp = atom1.damp*atom2.damp;
    if (pdamp != 0) {
   
        real ratio = r/pdamp;
        float pGamma = atom2.thole > atom1.thole ? atom1.thole : atom2.thole;

        real damp = ratio*ratio*ratio*pGamma;
        real dampExp = EXP(-damp);
        real damp1 = damp + 1;
        real damp2 = damp*damp;

        scale3 = 1 - dampExp;
        scale5 = 1 - damp1*dampExp;
        scale7 = 1 - (damp1 + 0.6f*damp2)*dampExp;

#ifdef F1
        real factor = 3*damp*dampExp*rr2;
        real factor7 = -0.2f + 0.6f*damp;
        
        ddsc3_0 = factor*xr;
        ddsc5_0 = ddsc3_0*damp;
        ddsc7_0 = ddsc5_0*factor7;

        ddsc3_1 = factor*yr;
        ddsc5_1 = ddsc3_1*damp;
        ddsc7_1 = ddsc5_1*factor7;

        ddsc3_2 = factor*zr;
        ddsc5_2 = ddsc3_2*damp;
        ddsc7_2 = ddsc5_2*factor7;
#endif

    }
      
#if defined F1
    real scale3i = rr3*scale3*uScale;
    real scale5i = rr5*scale5*uScale;
#endif
    real dsc3 = rr3*scale3*dScale;
    real psc3 = rr3*scale3*pScale;

    real dsc5 = rr5*scale5*dScale;
    real psc5 = rr5*scale5*pScale;

    real dsc7 = rr7*scale7*dScale;
    real psc7 = rr7*scale7*pScale;

    real atom2quadrupoleZZ = -(atom2.quadrupoleXX+atom2.quadrupoleYY);
    real qJr_0 = atom2.quadrupoleXX*xr + atom2.quadrupoleXY*yr + atom2.quadrupoleXZ*zr;
    real qJr_1 = atom2.quadrupoleXY*xr + atom2.quadrupoleYY*yr + atom2.quadrupoleYZ*zr;
    real qJr_2 = atom2.quadrupoleXZ*xr + atom2.quadrupoleYZ*yr + atom2quadrupoleZZ*zr;

    real atom1quadrupoleZZ = -(atom1.quadrupoleXX+atom1.quadrupoleYY);
    real qIr_0 = atom1.quadrupoleXX*xr + atom1.quadrupoleXY*yr + atom1.quadrupoleXZ*zr;
    real qIr_1 = atom1.quadrupoleXY*xr + atom1.quadrupoleYY*yr + atom1.quadrupoleYZ*zr;
    real qIr_2 = atom1.quadrupoleXZ*xr + atom1.quadrupoleYZ*yr + atom1quadrupoleZZ*zr;

#if defined F1
    real sc2 = atom1.dipole.x*atom2.dipole.x + atom1.dipole.y*atom2.dipole.y + atom1.dipole.z*atom2.dipole.z;
#endif
#if defined F1 || defined T1
    real sc4 = atom2.dipole.x*xr + atom2.dipole.y*yr + atom2.dipole.z*zr;
    real sc6 = qJr_0*xr + qJr_1*yr + qJr_2*zr;
#endif

#if defined F1 || defined T3
    real sc3 = atom1.dipole.x*xr + atom1.dipole.y*yr + atom1.dipole.z*zr;
    real sc5 = qIr_0*xr + qIr_1*yr + qIr_2*zr;
#endif
    
#if defined F1
    real sc7 = qIr_0*atom2.dipole.x + qIr_1*atom2.dipole.y + qIr_2*atom2.dipole.z;
    real sc8 = qJr_0*atom1.dipole.x + qJr_1*atom1.dipole.y + qJr_2*atom1.dipole.z;
    real sc9 = qIr_0*qJr_0 + qIr_1*qJr_1 + qIr_2*qJr_2;

    real sc10 = atom1.quadrupoleXX*atom2.quadrupoleXX + atom1.quadrupoleXY*atom2.quadrupoleXY + atom1.quadrupoleXZ*atom2.quadrupoleXZ +
                atom1.quadrupoleXY*atom2.quadrupoleXY + atom1.quadrupoleYY*atom2.quadrupoleYY + atom1.quadrupoleYZ*atom2.quadrupoleYZ +
                atom1.quadrupoleXZ*atom2.quadrupoleXZ + atom1.quadrupoleYZ*atom2.quadrupoleYZ + atom1quadrupoleZZ*atom2quadrupoleZZ;

    real sci1 = atom1.inducedDipole.x*atom2.dipole.x + atom1.inducedDipole.y*atom2.dipole.y + atom1.inducedDipole.z*atom2.dipole.z +
                atom2.inducedDipole.x*atom1.dipole.x + atom2.inducedDipole.y*atom1.dipole.y + atom2.inducedDipole.z*atom1.dipole.z;
#endif
        
#if defined F1 || defined T3
    real sci3 = atom1.inducedDipole.x*xr + atom1.inducedDipole.y*yr + atom1.inducedDipole.z*zr;
#endif
#if defined F1
    real sci7 = qIr_0*atom2.inducedDipole.x + qIr_1*atom2.inducedDipole.y + qIr_2*atom2.inducedDipole.z;
    real sci8 = qJr_0*atom1.inducedDipole.x + qJr_1*atom1.inducedDipole.y + qJr_2*atom1.inducedDipole.z;
#endif
#if defined F1 || defined T1
    real sci4 = atom2.inducedDipole.x*xr + atom2.inducedDipole.y*yr + atom2.inducedDipole.z*zr;
#endif
    
#if defined F1
    real scip1 = atom1.inducedDipolePolar.x*atom2.dipole.x + atom1.inducedDipolePolar.y*atom2.dipole.y + atom1.inducedDipolePolar.z*atom2.dipole.z +
                 atom2.inducedDipolePolar.x*atom1.dipole.x + atom2.inducedDipolePolar.y*atom1.dipole.y + atom2.inducedDipolePolar.z*atom1.dipole.z;

    real scip2 = atom1.inducedDipole.x*atom2.inducedDipolePolar.x + atom1.inducedDipole.y*atom2.inducedDipolePolar.y + atom1.inducedDipole.z*atom2.inducedDipolePolar.z +
                 atom2.inducedDipole.x*atom1.inducedDipolePolar.x + atom2.inducedDipole.y*atom1.inducedDipolePolar.y + atom2.inducedDipole.z*atom1.inducedDipolePolar.z;

#endif
    
#if defined F1 || defined T3
    real scip3 = ((atom1.inducedDipolePolar.x)*(xr) + (atom1.inducedDipolePolar.y)*(yr) + (atom1.inducedDipolePolar.z)*(zr));
#endif
#if defined F1 || defined T1
    real scip4 = ((atom2.inducedDipolePolar.x)*(xr) + (atom2.inducedDipolePolar.y)*(yr) + (atom2.inducedDipolePolar.z)*(zr));
#endif

#ifdef F1
    real scip7 = ((qIr_0)*(atom2.inducedDipolePolar.x) + (qIr_1)*(atom2.inducedDipolePolar.y) + (qIr_2)*(atom2.inducedDipolePolar.z));
    real scip8 = ((qJr_0)*(atom1.inducedDipolePolar.x) + (qJr_1)*(atom1.inducedDipolePolar.y) + (qJr_2)*(atom1.inducedDipolePolar.z));

    real gli1 = atom2.posq.w*sci3 - atom1.posq.w*sci4;
    
    real gli6 = sci1;
    real glip1 = atom2.posq.w*scip3 - atom1.posq.w*scip4;
    real glip6 = scip1;
    real gli2 = -sc3*sci4 - sci3*sc4;
    real gli3 = sci3*sc6 - sci4*sc5;
    real gli7 = 2*(sci7-sci8);
    
    real glip2 = -sc3*scip4 - scip3*sc4;
    real glip3 = scip3*sc6 - scip4*sc5;
    real glip7 = 2*(scip7-scip8);
    real factor3 = rr3*((gli1  +  gli6)*pScale + (glip1  + glip6)*dScale);
    real factor5 = rr5*((gli2  +  gli7)*pScale + (glip2  + glip7)*dScale);
    real factor7 = rr7*(gli3*pScale + glip3*dScale);
      
    real ftm2i_0 = -0.5f*(factor3*ddsc3_0 + factor5*ddsc5_0 + factor7*ddsc7_0);
    real ftm2i_1 = -0.5f*(factor3*ddsc3_1 + factor5*ddsc5_1 + factor7*ddsc7_1);
    real ftm2i_2 = -0.5f*(factor3*ddsc3_2 + factor5*ddsc5_2 + factor7*ddsc7_2);
      
    real gl0 = atom1.posq.w*atom2.posq.w;
    real gl1 = atom2.posq.w*sc3 - atom1.posq.w*sc4;
    real gl2 = atom1.posq.w*sc6 + atom2.posq.w*sc5 - sc3*sc4;
    real gl3 = sc3*sc6 - sc4*sc5;
    real gl4 = sc5*sc6;
    real gl6 = sc2;
    real gl7 = 2*(sc7-sc8);
    real gl8 = 2*sc10;
    real gl5 = -4*sc9;
    
    real gf1 = rr3*gl0 + rr5*(gl1+gl6) + rr7*(gl2+gl7+gl8) + rr9*(gl3+gl5) + rr11*gl4;
#endif
#if defined F1 || defined T1
    real gf2 = -atom2.posq.w*rr3 + sc4*rr5 - sc6*rr7;
    real gf5 = 2*(-atom2.posq.w*rr5+sc4*rr7-sc6*rr9);
#endif
#if defined F1 || defined T3
    real gf3 =  atom1.posq.w*rr3 + sc3*rr5 + sc5*rr7;
    real gf6 = 2*(-atom1.posq.w*rr5-sc3*rr7-sc5*rr9);
#endif

#ifdef F1
    real em = mScale*(rr1*gl0 + rr3*(gl1+gl6) + rr5*(gl2+gl7+gl8) + rr7*(gl3+gl5) + rr9*gl4);
    real ei = 0.5f*((gli1+gli6)*psc3 + (gli2+gli7)*psc5 + gli3*psc7);
    energy = em+ei;
#endif
    
#if defined F1 || defined T1
    real qIdJ_0 = atom1.quadrupoleXX*atom2.dipole.x + atom1.quadrupoleXY*atom2.dipole.y + atom1.quadrupoleXZ*atom2.dipole.z;
    real qIdJ_1 = atom1.quadrupoleXY*atom2.dipole.x + atom1.quadrupoleYY*atom2.dipole.y + atom1.quadrupoleYZ*atom2.dipole.z;
    real qIdJ_2 = atom1.quadrupoleXZ*atom2.dipole.x + atom1.quadrupoleYZ*atom2.dipole.y + atom1quadrupoleZZ*atom2.dipole.z;

    real qIqJr_0 = atom1.quadrupoleXX*qJr_0 + atom1.quadrupoleXY*qJr_1 + atom1.quadrupoleXZ*qJr_2;
    real qIqJr_1 = atom1.quadrupoleXY*qJr_0 + atom1.quadrupoleYY*qJr_1 + atom1.quadrupoleYZ*qJr_2;
    real qIqJr_2 = atom1.quadrupoleXZ*qJr_0 + atom1.quadrupoleYZ*qJr_1 + atom1quadrupoleZZ*qJr_2;
#endif

#ifdef F1
    real qkqir_0 = atom2.quadrupoleXX*qIr_0 + atom2.quadrupoleXY*qIr_1 + atom2.quadrupoleXZ*qIr_2;
    real qkqir_1 = atom2.quadrupoleXY*qIr_0 + atom2.quadrupoleYY*qIr_1 + atom2.quadrupoleYZ*qIr_2;
    real qkqir_2 = atom2.quadrupoleXZ*qIr_0 + atom2.quadrupoleYZ*qIr_1 + atom2quadrupoleZZ*qIr_2;

    real qkdi_0 = atom2.quadrupoleXX*atom1.dipole.x + atom2.quadrupoleXY*atom1.dipole.y + atom2.quadrupoleXZ*atom1.dipole.z;
    real qkdi_1 = atom2.quadrupoleXY*atom1.dipole.x + atom2.quadrupoleYY*atom1.dipole.y + atom2.quadrupoleYZ*atom1.dipole.z;
    real qkdi_2 = atom2.quadrupoleXZ*atom1.dipole.x + atom2.quadrupoleYZ*atom1.dipole.y + atom2quadrupoleZZ*atom1.dipole.z;

    real ftm2_0 = mScale*(gf1*xr + gf2*atom1.dipole.x + gf3*atom2.dipole.x + 2*rr5*(qkdi_0 - qIdJ_0) + gf5*qIr_0 + gf6*qJr_0 + 4*rr7*(qIqJr_0 + qkqir_0));
    real ftm2_1 = mScale*(gf1*yr + gf2*atom1.dipole.y + gf3*atom2.dipole.y + 2*rr5*(qkdi_1 - qIdJ_1) + gf5*qIr_1 + gf6*qJr_1 + 4*rr7*(qIqJr_1 + qkqir_1));
    real ftm2_2 = mScale*(gf1*zr + gf2*atom1.dipole.z + gf3*atom2.dipole.z + 2*rr5*(qkdi_2 - qIdJ_2) + gf5*qIr_2 + gf6*qJr_2 + 4*rr7*(qIqJr_2 + qkqir_2));

    real gfi1 = rr2*(1.5f*((gli1+gli6)*psc3 + (glip1+glip6)*dsc3 + scip2*scale3i) + 2.5f*((gli7+gli2)*psc5 + (glip7+glip2)*dsc5 - (sci3*scip4+scip3*sci4)*scale5i) + 3.5f*(gli3*psc7+glip3*dsc7));
    ftm2i_0 += gfi1*xr;
    ftm2i_1 += gfi1*yr;
    ftm2i_2 += gfi1*zr;
#endif

#if defined F1 || defined T1
    real gfi5 = (sci4*psc7 + scip4*dsc7);
#endif
#if defined F1 || defined T3
    real gfi6 = -(sci3*psc7 + scip3*dsc7);
#endif

#if defined F1 || defined T1
    real qIuJ_0 = atom1.quadrupoleXX*atom2.inducedDipole.x   + atom1.quadrupoleXY*atom2.inducedDipole.y  + atom1.quadrupoleXZ*atom2.inducedDipole.z;
    real qIuJ_1 = atom1.quadrupoleXY*atom2.inducedDipole.x   + atom1.quadrupoleYY*atom2.inducedDipole.y  + atom1.quadrupoleYZ*atom2.inducedDipole.z;
    real qIuJ_2 = atom1.quadrupoleXZ*atom2.inducedDipole.x   + atom1.quadrupoleYZ*atom2.inducedDipole.y  + atom1quadrupoleZZ*atom2.inducedDipole.z;

    real qIuJp_0 = atom1.quadrupoleXX*atom2.inducedDipolePolar.x + atom1.quadrupoleXY*atom2.inducedDipolePolar.y + atom1.quadrupoleXZ*atom2.inducedDipolePolar.z;
    real qIuJp_1 = atom1.quadrupoleXY*atom2.inducedDipolePolar.x + atom1.quadrupoleYY*atom2.inducedDipolePolar.y + atom1.quadrupoleYZ*atom2.inducedDipolePolar.z;
    real qIuJp_2 = atom1.quadrupoleXZ*atom2.inducedDipolePolar.x + atom1.quadrupoleYZ*atom2.inducedDipolePolar.y + atom1quadrupoleZZ*atom2.inducedDipolePolar.z;
#endif

#if defined T3
    real qJuIp_0 = atom2.quadrupoleXX*atom1.inducedDipolePolar.x + atom2.quadrupoleXY*atom1.inducedDipolePolar.y + atom2.quadrupoleXZ*atom1.inducedDipolePolar.z;
    real qJuIp_1 = atom2.quadrupoleXY*atom1.inducedDipolePolar.x + atom2.quadrupoleYY*atom1.inducedDipolePolar.y + atom2.quadrupoleYZ*atom1.inducedDipolePolar.z;
    real qJuIp_2 = atom2.quadrupoleXZ*atom1.inducedDipolePolar.x + atom2.quadrupoleYZ*atom1.inducedDipolePolar.y + atom2quadrupoleZZ*atom1.inducedDipolePolar.z;

     real qJuI_0 = atom2.quadrupoleXX*atom1.inducedDipole.x + atom2.quadrupoleXY*atom1.inducedDipole.y + atom2.quadrupoleXZ*atom1.inducedDipole.z;
     real qJuI_1 = atom2.quadrupoleXY*atom1.inducedDipole.x + atom2.quadrupoleYY*atom1.inducedDipole.y + atom2.quadrupoleYZ*atom1.inducedDipole.z;
     real qJuI_2 = atom2.quadrupoleXZ*atom1.inducedDipole.x + atom2.quadrupoleYZ*atom1.inducedDipole.y + atom2quadrupoleZZ*atom1.inducedDipole.z;
#endif

#ifdef F1
    real qkui_0 = atom2.quadrupoleXX*atom1.inducedDipole.x + atom2.quadrupoleXY*atom1.inducedDipole.y + atom2.quadrupoleXZ*atom1.inducedDipole.z;
    real qkui_1 = atom2.quadrupoleXY*atom1.inducedDipole.x + atom2.quadrupoleYY*atom1.inducedDipole.y + atom2.quadrupoleYZ*atom1.inducedDipole.z;
    real qkui_2 = atom2.quadrupoleXZ*atom1.inducedDipole.x + atom2.quadrupoleYZ*atom1.inducedDipole.y + atom2quadrupoleZZ*atom1.inducedDipole.z;

    real qkuip_0 = atom2.quadrupoleXX*atom1.inducedDipolePolar.x + atom2.quadrupoleXY*atom1.inducedDipolePolar.y + atom2.quadrupoleXZ*atom1.inducedDipolePolar.z;
    real qkuip_1 = atom2.quadrupoleXY*atom1.inducedDipolePolar.x + atom2.quadrupoleYY*atom1.inducedDipolePolar.y + atom2.quadrupoleYZ*atom1.inducedDipolePolar.z;
    real qkuip_2 = atom2.quadrupoleXZ*atom1.inducedDipolePolar.x + atom2.quadrupoleYZ*atom1.inducedDipolePolar.y + atom2quadrupoleZZ*atom1.inducedDipolePolar.z;

    ftm2i_0 += 0.5f*(-atom2.posq.w*(atom1.inducedDipole.x*psc3 + atom1.inducedDipolePolar.x*dsc3) +
               sc4*(atom1.inducedDipole.x*psc5 + atom1.inducedDipolePolar.x*dsc5) -
               sc6*(atom1.inducedDipole.x*psc7 + atom1.inducedDipolePolar.x*dsc7)) +
      
               0.5f*(atom1.posq.w*(atom2.inducedDipole.x*psc3+atom2.inducedDipolePolar.x*dsc3) +
               sc3*(atom2.inducedDipole.x*psc5 +atom2.inducedDipolePolar.x*dsc5) +
               sc5*(atom2.inducedDipole.x*psc7 +atom2.inducedDipolePolar.x*dsc7)) +

               scale5i*(sci4*atom1.inducedDipolePolar.x+scip4*atom1.inducedDipole.x +
                        sci3*atom2.inducedDipolePolar.x+scip3*atom2.inducedDipole.x)*0.5f +
      
               0.5f*(sci4*psc5+scip4*dsc5)*atom1.dipole.x +
               0.5f*(sci3*psc5+scip3*dsc5)*atom2.dipole.x +
               ((qkui_0-qIuJ_0)*psc5 + (qkuip_0-qIuJp_0)*dsc5) +
               gfi5*qIr_0 + gfi6*qJr_0;
      
    ftm2i_1 += 0.5f*(-atom2.posq.w*(atom1.inducedDipole.y*psc3 + atom1.inducedDipolePolar.y*dsc3) +
               sc4*(atom1.inducedDipole.y*psc5 + atom1.inducedDipolePolar.y*dsc5) -
               sc6*(atom1.inducedDipole.y*psc7 + atom1.inducedDipolePolar.y*dsc7)) +

               (atom1.posq.w*(atom2.inducedDipole.y*psc3+atom2.inducedDipolePolar.y*dsc3) +
                    sc3*(atom2.inducedDipole.y*psc5+atom2.inducedDipolePolar.y*dsc5) +
                    sc5*(atom2.inducedDipole.y*psc7+atom2.inducedDipolePolar.y*dsc7))*0.5f +
                    scale5i*(sci4*atom1.inducedDipolePolar.y+scip4*atom1.inducedDipole.y + sci3*atom2.inducedDipolePolar.y+scip3*atom2.inducedDipole.y)*0.5f +

               0.5f*(sci4*psc5+scip4*dsc5)*atom1.dipole.y +
               0.5f*(sci3*psc5+scip3*dsc5)*atom2.dipole.y +
               ((qkui_1-qIuJ_1)*psc5 + (qkuip_1-qIuJp_1)*dsc5) +
               gfi5*qIr_1 + gfi6*qJr_1;
      
    ftm2i_2 += 0.5f*(-atom2.posq.w*(atom1.inducedDipole.z*psc3 + atom1.inducedDipolePolar.z*dsc3) +
               sc4*(atom1.inducedDipole.z*psc5 + atom1.inducedDipolePolar.z*dsc5) -
               sc6*(atom1.inducedDipole.z*psc7 + atom1.inducedDipolePolar.z*dsc7)) +

               (atom1.posq.w*(atom2.inducedDipole.z*psc3+atom2.inducedDipolePolar.z*dsc3) +
                    sc3*(atom2.inducedDipole.z*psc5+atom2.inducedDipolePolar.z*dsc5) +
                    sc5*(atom2.inducedDipole.z*psc7+atom2.inducedDipolePolar.z*dsc7))*0.5f +
                    scale5i*(sci4*atom1.inducedDipolePolar.z+scip4*atom1.inducedDipole.z +
                    sci3*atom2.inducedDipolePolar.z+scip3*atom2.inducedDipole.z)*0.5f +

               0.5f*(sci4*psc5+scip4*dsc5)*atom1.dipole.z +
               0.5f*(sci3*psc5+scip3*dsc5)*atom2.dipole.z +
               ((qkui_2-qIuJ_2)*psc5 + (qkuip_2-qIuJp_2)*dsc5) +
               gfi5*qIr_2 + gfi6*qJr_2;

#ifdef DIRECT_POLARIZATION
    real gfd = 0.5*(3*rr2*scip2*scale3i - 5*rr2*(scip3*sci4+sci3*scip4)*scale5i);
    real temp5 = 0.5*scale5i;
    real fdir_0 = gfd*xr + temp5*(sci4*atom1.inducedDipolePolar.x + scip4*atom1.inducedDipole.x + sci3*atom2.inducedDipolePolar.x + scip3*atom2.inducedDipole.x);
    real fdir_1 = gfd*yr + temp5*(sci4*atom1.inducedDipolePolar.y + scip4*atom1.inducedDipole.y + sci3*atom2.inducedDipolePolar.y + scip3*atom2.inducedDipole.y);
    real fdir_2 = gfd*zr + temp5*(sci4*atom1.inducedDipolePolar.z + scip4*atom1.inducedDipole.z + sci3*atom2.inducedDipolePolar.z + scip3*atom2.inducedDipole.z);
    ftm2i_0 -= fdir_0;
    ftm2i_1 -= fdir_1;
    ftm2i_2 -= fdir_2;
#else
    real scaleF = 0.5f*uScale;
    real inducedFactor3 = scip2*rr3*scaleF;
    real inducedFactor5 = (sci3*scip4+scip3*sci4)*rr5*scaleF;
    real findmp_0 = inducedFactor3*ddsc3_0 - inducedFactor5*ddsc5_0;
    real findmp_1 = inducedFactor3*ddsc3_1 - inducedFactor5*ddsc5_1;
    real findmp_2 = inducedFactor3*ddsc3_2 - inducedFactor5*ddsc5_2;
    ftm2i_0 -= findmp_0;
    ftm2i_1 -= findmp_1;
    ftm2i_2 -= findmp_2;
#endif
#endif

#if defined T1
    real gti2 = 0.5f*(sci4*psc5+scip4*dsc5);
    real gti5 = gfi5;
#endif
#if defined T3
    real gti3 = 0.5f*(sci3*psc5+scip3*dsc5);
    real gti6 = gfi6;
#endif

#if defined T1 || defined T3
    real dixdk_0 = atom1.dipole.y*atom2.dipole.z - atom1.dipole.z*atom2.dipole.y;
    real dixdk_1 = atom1.dipole.z*atom2.dipole.x - atom1.dipole.x*atom2.dipole.z;
    real dixdk_2 = atom1.dipole.x*atom2.dipole.y - atom1.dipole.y*atom2.dipole.x;

#if defined T1
    real dixuk_0 = atom1.dipole.y*atom2.inducedDipole.z - atom1.dipole.z*atom2.inducedDipole.y;
    real dixuk_1 = atom1.dipole.z*atom2.inducedDipole.x - atom1.dipole.x*atom2.inducedDipole.z;
    real dixuk_2 = atom1.dipole.x*atom2.inducedDipole.y - atom1.dipole.y*atom2.inducedDipole.x;
#endif
#endif

#ifdef T1
    real dixukp_0 = atom1.dipole.y*atom2.inducedDipolePolar.z - atom1.dipole.z*atom2.inducedDipolePolar.y;
    real dixukp_1 = atom1.dipole.z*atom2.inducedDipolePolar.x - atom1.dipole.x*atom2.inducedDipolePolar.z;
    real dixukp_2 = atom1.dipole.x*atom2.inducedDipolePolar.y - atom1.dipole.y*atom2.inducedDipolePolar.x;
#endif

#ifdef T1
    real dixr_0 = atom1.dipole.y*zr - atom1.dipole.z*yr;
    real dixr_1 = atom1.dipole.z*xr - atom1.dipole.x*zr;
    real dixr_2 = atom1.dipole.x*yr - atom1.dipole.y*xr;
#endif

#ifdef T1
    real rxqiukp_0 = yr*qIuJp_2 - zr*qIuJp_1;
    real rxqiukp_1 = zr*qIuJp_0 - xr*qIuJp_2;
    real rxqiukp_2 = xr*qIuJp_1 - yr*qIuJp_0;

    real rxqir_0 = yr*qIr_2 - zr*qIr_1;
    real rxqir_1 = zr*qIr_0 - xr*qIr_2;
    real rxqir_2 = xr*qIr_1 - yr*qIr_0;

    real rxqiuk_0 = yr*qIuJ_2 - zr*qIuJ_1;
    real rxqiuk_1 = zr*qIuJ_0 - xr*qIuJ_2;
    real rxqiuk_2 = xr*qIuJ_1 - yr*qIuJ_0;

    real ukxqir_0 = atom2.inducedDipole.y*qIr_2 - atom2.inducedDipole.z*qIr_1;
    real ukxqir_1 = atom2.inducedDipole.z*qIr_0 - atom2.inducedDipole.x*qIr_2;
    real ukxqir_2 = atom2.inducedDipole.x*qIr_1 - atom2.inducedDipole.y*qIr_0;

    real ukxqirp_0 = atom2.inducedDipolePolar.y*qIr_2 - atom2.inducedDipolePolar.z*qIr_1;
    real ukxqirp_1 = atom2.inducedDipolePolar.z*qIr_0 - atom2.inducedDipolePolar.x*qIr_2;
    real ukxqirp_2 = atom2.inducedDipolePolar.x*qIr_1 - atom2.inducedDipolePolar.y*qIr_0;

    real dixqkr_0 = atom1.dipole.y*qJr_2 - atom1.dipole.z*qJr_1;
    real dixqkr_1 = atom1.dipole.z*qJr_0 - atom1.dipole.x*qJr_2;
    real dixqkr_2 = atom1.dipole.x*qJr_1 - atom1.dipole.y*qJr_0;

    real dkxqir_0 = atom2.dipole.y*qIr_2 - atom2.dipole.z*qIr_1;
    real dkxqir_1 = atom2.dipole.z*qIr_0 - atom2.dipole.x*qIr_2;
    real dkxqir_2 = atom2.dipole.x*qIr_1 - atom2.dipole.y*qIr_0;

    real rxqikr_0 = yr*qIqJr_2 - zr*qIqJr_1;
    real rxqikr_1 = zr*qIqJr_0 - xr*qIqJr_2;
    real rxqikr_2 = xr*qIqJr_1 - yr*qIqJr_0;

    real rxqidk_0 = yr*qIdJ_2 - zr*qIdJ_1;
    real rxqidk_1 = zr*qIdJ_0 - xr*qIdJ_2;
    real rxqidk_2 = xr*qIdJ_1 - yr*qIdJ_0;

    real qkrxqir_0 = qJr_1*qIr_2 - qJr_2*qIr_1;
    real qkrxqir_1 = qJr_2*qIr_0 - qJr_0*qIr_2;
    real qkrxqir_2 = qJr_0*qIr_1 - qJr_1*qIr_0;
#endif

#if defined T1 || defined T3
    real qixqk_0 = atom1.quadrupoleXY*atom2.quadrupoleXZ + atom1.quadrupoleYY*atom2.quadrupoleYZ + atom1.quadrupoleYZ*atom2quadrupoleZZ -
                   atom1.quadrupoleXZ*atom2.quadrupoleXY - atom1.quadrupoleYZ*atom2.quadrupoleYY - atom1quadrupoleZZ*atom2.quadrupoleYZ;

    real qixqk_1 = atom1.quadrupoleXZ*atom2.quadrupoleXX + atom1.quadrupoleYZ*atom2.quadrupoleXY + atom1quadrupoleZZ*atom2.quadrupoleXZ -
                   atom1.quadrupoleXX*atom2.quadrupoleXZ - atom1.quadrupoleXY*atom2.quadrupoleYZ - atom1.quadrupoleXZ*atom2quadrupoleZZ;

    real qixqk_2 = atom1.quadrupoleXX*atom2.quadrupoleXY + atom1.quadrupoleXY*atom2.quadrupoleYY + atom1.quadrupoleXZ*atom2.quadrupoleYZ -
                   atom1.quadrupoleXY*atom2.quadrupoleXX - atom1.quadrupoleYY*atom2.quadrupoleXY - atom1.quadrupoleYZ*atom2.quadrupoleXZ;
#endif

#ifdef T1
    real ttm2_0 = -rr3*dixdk_0 + gf2*dixr_0-gf5*rxqir_0 + 2*rr5*(dixqkr_0 + dkxqir_0 + rxqidk_0-2*qixqk_0) - 4*rr7*(rxqikr_0 + qkrxqir_0);
    real ttm2_1 = -rr3*dixdk_1 + gf2*dixr_1-gf5*rxqir_1 + 2*rr5*(dixqkr_1 + dkxqir_1 + rxqidk_1-2*qixqk_1) - 4*rr7*(rxqikr_1 + qkrxqir_1);
    real ttm2_2 = -rr3*dixdk_2 + gf2*dixr_2-gf5*rxqir_2 + 2*rr5*(dixqkr_2 + dkxqir_2 + rxqidk_2-2*qixqk_2) - 4*rr7*(rxqikr_2 + qkrxqir_2);

    real ttm2i_0 = -(dixuk_0*psc3+dixukp_0*dsc3)*0.5f + gti2*dixr_0 + ((ukxqir_0+ rxqiuk_0)*psc5 + (ukxqirp_0 + rxqiukp_0)*dsc5) - gti5*rxqir_0;
    real ttm2i_1 = -(dixuk_1*psc3+dixukp_1*dsc3)*0.5f + gti2*dixr_1 + ((ukxqir_1+ rxqiuk_1)*psc5 + (ukxqirp_1 + rxqiukp_1)*dsc5) - gti5*rxqir_1;
    real ttm2i_2 = -(dixuk_2*psc3+dixukp_2*dsc3)*0.5f + gti2*dixr_2 + ((ukxqir_2+ rxqiuk_2)*psc5 + (ukxqirp_2 + rxqiukp_2)*dsc5) - gti5*rxqir_2;
#endif

#ifdef T3
    real qJqIr_0 = atom2.quadrupoleXX*qIr_0 + atom2.quadrupoleXY*qIr_1 + atom2.quadrupoleXZ*qIr_2;
    real qJqIr_1 = atom2.quadrupoleXY*qIr_0 + atom2.quadrupoleYY*qIr_1 + atom2.quadrupoleYZ*qIr_2;
    real qJqIr_2 = atom2.quadrupoleXZ*qIr_0 + atom2.quadrupoleYZ*qIr_1 + atom2quadrupoleZZ*qIr_2;

    real qJdI_0 = atom2.quadrupoleXX*atom1.dipole.x + atom2.quadrupoleXY*atom1.dipole.y + atom2.quadrupoleXZ*atom1.dipole.z;
    real qJdI_1 = atom2.quadrupoleXY*atom1.dipole.x + atom2.quadrupoleYY*atom1.dipole.y + atom2.quadrupoleYZ*atom1.dipole.z;
    real qJdI_2 = atom2.quadrupoleXZ*atom1.dipole.x + atom2.quadrupoleYZ*atom1.dipole.y + atom2quadrupoleZZ*atom1.dipole.z;

    real dkxr_0 = atom2.dipole.y*zr - atom2.dipole.z*yr;
    real dkxr_1 = atom2.dipole.z*xr - atom2.dipole.x*zr;
    real dkxr_2 = atom2.dipole.x*yr - atom2.dipole.y*xr;

    real rxqkr_0 = yr*qJr_2 - zr*qJr_1;
    real rxqkr_1 = zr*qJr_0 - xr*qJr_2;
    real rxqkr_2 = xr*qJr_1 - yr*qJr_0;

    real dixqkr_0 = atom1.dipole.y*qJr_2 - atom1.dipole.z*qJr_1;
    real dixqkr_1 = atom1.dipole.z*qJr_0 - atom1.dipole.x*qJr_2;
    real dixqkr_2 = atom1.dipole.x*qJr_1 - atom1.dipole.y*qJr_0;

    real dkxqir_0 = atom2.dipole.y*qIr_2 - atom2.dipole.z*qIr_1;
    real dkxqir_1 = atom2.dipole.z*qIr_0 - atom2.dipole.x*qIr_2;
    real dkxqir_2 = atom2.dipole.x*qIr_1 - atom2.dipole.y*qIr_0;

    real rxqkdi_0 = yr*qJdI_2 - zr*qJdI_1;
    real rxqkdi_1 = zr*qJdI_0 - xr*qJdI_2;
    real rxqkdi_2 = xr*qJdI_1 - yr*qJdI_0;

    real rxqkir_0 = yr*qJqIr_2 - zr*qJqIr_1;
    real rxqkir_1 = zr*qJqIr_0 - xr*qJqIr_2;
    real rxqkir_2 = xr*qJqIr_1 - yr*qJqIr_0;

    real qkrxqir_0 = qJr_1*qIr_2 - qJr_2*qIr_1;
    real qkrxqir_1 = qJr_2*qIr_0 - qJr_0*qIr_2;
    real qkrxqir_2 = qJr_0*qIr_1 - qJr_1*qIr_0;

    real dkxui_0 = atom2.dipole.y*atom1.inducedDipole.z - atom2.dipole.z*atom1.inducedDipole.y;
    real dkxui_1 = atom2.dipole.z*atom1.inducedDipole.x - atom2.dipole.x*atom1.inducedDipole.z; 
    real dkxui_2 = atom2.dipole.x*atom1.inducedDipole.y - atom2.dipole.y*atom1.inducedDipole.x;

    real dkxuip_0 = atom2.dipole.y*atom1.inducedDipolePolar.z - atom2.dipole.z*atom1.inducedDipolePolar.y;
    real dkxuip_1 = atom2.dipole.z*atom1.inducedDipolePolar.x - atom2.dipole.x*atom1.inducedDipolePolar.z;
    real dkxuip_2 = atom2.dipole.x*atom1.inducedDipolePolar.y - atom2.dipole.y*atom1.inducedDipolePolar.x;

    real uixqkrp_0 = atom1.inducedDipolePolar.y*qJr_2 - atom1.inducedDipolePolar.z*qJr_1;
    real uixqkrp_1 = atom1.inducedDipolePolar.z*qJr_0 - atom1.inducedDipolePolar.x*qJr_2;
    real uixqkrp_2 = atom1.inducedDipolePolar.x*qJr_1 - atom1.inducedDipolePolar.y*qJr_0;

    real uixqkr_0 = atom1.inducedDipole.y*qJr_2 - atom1.inducedDipole.z*qJr_1;
    real uixqkr_1 = atom1.inducedDipole.z*qJr_0 - atom1.inducedDipole.x*qJr_2;
    real uixqkr_2 = atom1.inducedDipole.x*qJr_1 - atom1.inducedDipole.y*qJr_0;

    real rxqkuip_0 = yr*qJuIp_2 - zr*qJuIp_1;
    real rxqkuip_1 = zr*qJuIp_0 - xr*qJuIp_2;
    real rxqkuip_2 = xr*qJuIp_1 - yr*qJuIp_0;

    real rxqkui_0 = yr*qJuI_2 - zr*qJuI_1;
    real rxqkui_1 = zr*qJuI_0 - xr*qJuI_2;
    real rxqkui_2 = xr*qJuI_1 - yr*qJuI_0;

    real ttm3_0 =  rr3*dixdk_0 + gf3*dkxr_0 - gf6*rxqkr_0 - 2*rr5*(dixqkr_0 + dkxqir_0 + rxqkdi_0 - 2*qixqk_0) - 4*rr7*(rxqkir_0 - qkrxqir_0);
    real ttm3_1 =  rr3*dixdk_1 + gf3*dkxr_1 - gf6*rxqkr_1 - 2*rr5*(dixqkr_1 + dkxqir_1 + rxqkdi_1 - 2*qixqk_1) - 4*rr7*(rxqkir_1 - qkrxqir_1);
    real ttm3_2 =  rr3*dixdk_2 + gf3*dkxr_2 - gf6*rxqkr_2 - 2*rr5*(dixqkr_2 + dkxqir_2 + rxqkdi_2 - 2*qixqk_2) - 4*rr7*(rxqkir_2 - qkrxqir_2);

    real ttm3i_0 = -(dkxui_0*psc3+ dkxuip_0*dsc3)*0.5f + gti3*dkxr_0 - ((uixqkr_0 + rxqkui_0)*psc5 + (uixqkrp_0 + rxqkuip_0)*dsc5) - gti6*rxqkr_0;
    real ttm3i_1 = -(dkxui_1*psc3+ dkxuip_1*dsc3)*0.5f + gti3*dkxr_1 - ((uixqkr_1 + rxqkui_1)*psc5 + (uixqkrp_1 + rxqkuip_1)*dsc5) - gti6*rxqkr_1;
    real ttm3i_2 = -(dkxui_2*psc3+ dkxuip_2*dsc3)*0.5f + gti3*dkxr_2 - ((uixqkr_2 + rxqkui_2)*psc5 + (uixqkrp_2 + rxqkuip_2)*dsc5) - gti6*rxqkr_2;
#endif

    if (mScale < 1) {
#ifdef T1
        ttm2_0 *= mScale;
        ttm2_1 *= mScale;
        ttm2_2 *= mScale;
#endif
        
#ifdef T3
        ttm3_0 *= mScale;
        ttm3_1 *= mScale;
        ttm3_2 *= mScale;
#endif
    }

#ifdef F1
    outputForce.x = -(ftm2_0+ftm2i_0);
    outputForce.y = -(ftm2_1+ftm2i_1);
    outputForce.z = -(ftm2_2+ftm2i_2);
#endif
    
#ifdef T1
    outputForce.x = (ttm2_0 + ttm2i_0);
    outputForce.y = (ttm2_1 + ttm2i_1);
    outputForce.z = (ttm2_2 + ttm2i_2);
#endif

#ifdef T3
    outputForce.x = (ttm3_0 + ttm3i_0);
    outputForce.y = (ttm3_1 + ttm3i_1);
    outputForce.z = (ttm3_2 + ttm3i_2);
#endif
}
