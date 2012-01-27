
static __device__ void SUB_METHOD_NAME( calculatePmeDirectElectrostaticPairIxnF2, _kernel )( 
                                        PmeDirectElectrostaticParticle& atomI, PmeDirectElectrostaticParticle& atomJ,
                                        float4 delta, float4 bn, float forceFactor,
#ifdef APPLY_SCALE
                                        const float* scalingFactors,
#endif
                                        float force[3], float* energy ){

    float xr                    = delta.x;
    float yr                    = delta.y;
    float zr                    = delta.z;
    float rr1                   = delta.w;

    // set the permanent multipole and induced dipole values;

    float ci                    = atomI.q;

    float di1                   = atomI.labFrameDipole[0];
    float di2                   = atomI.labFrameDipole[1];
    float di3                   = atomI.labFrameDipole[2];

    float qi1                   = atomI.labFrameQuadrupole[0];
    float qi2                   = atomI.labFrameQuadrupole[1];
    float qi3                   = atomI.labFrameQuadrupole[2];
    float qi5                   = atomI.labFrameQuadrupole[3];
    float qi6                   = atomI.labFrameQuadrupole[4];
//    float qi9                   = atomI.labFrameQuadrupole[5];
    float qi9                   = -(atomI.labFrameQuadrupole[0] + atomI.labFrameQuadrupole[3]);

    float bn1                   = bn.x;
    float bn2                   = bn.y;
    float bn3                   = bn.z;
    float bn4                   = bn.w;

    float damp                  = atomI.damp*atomJ.damp;
    if( damp != 0.0f ){
        float pgamma = atomI.thole < atomJ.thole ? atomI.thole : atomJ.thole;
        float ratio  = 1.0f/(rr1*damp);
        damp         = -pgamma*ratio*ratio*ratio;
        damp         = damp < -50.0f ? 0.0f : damp;
    }

    float scale5                = (damp == 0.0f) ? 1.0f : (1.0f - (1.0f-damp)*expf(damp));
    float rr5                   = rr1*rr1;
          rr5                   = 3.0f*rr1*rr5*rr5;
#ifdef APPLY_SCALE
    float psc5                  = rr5*(1.0f - scale5*scalingFactors[PScaleIndex]);
    float dsc5                  = rr5*(1.0f - scale5*scalingFactors[DScaleIndex]);
    float usc5                  = rr5*(1.0f - scale5*scalingFactors[UScaleIndex]);
#else
    float psc5                  = rr5*(1.0f - scale5);
#endif

    float qiuk1                 = qi1*atomJ.inducedDipole[0]  + qi2*atomJ.inducedDipole[1]  + qi3*atomJ.inducedDipole[2];
    float qiukp1                = qi1*atomJ.inducedDipoleP[0] + qi2*atomJ.inducedDipoleP[1] + qi3*atomJ.inducedDipoleP[2];
    float ftm21                 = -bn2*(qiuk1+qiukp1);
#ifdef APPLY_SCALE
          ftm21                += qiuk1*psc5 + qiukp1*dsc5;
#else
          ftm21                += (qiuk1 + qiukp1)*psc5;
#endif

    float qiuk2                 = qi2*atomJ.inducedDipole[0]  + qi5*atomJ.inducedDipole[1]  + qi6*atomJ.inducedDipole[2];
    float qiukp2                = qi2*atomJ.inducedDipoleP[0] + qi5*atomJ.inducedDipoleP[1] + qi6*atomJ.inducedDipoleP[2];
    float ftm22                 = -bn2*(qiuk2+qiukp2);
#ifdef APPLY_SCALE
          ftm22                += ((qiuk2)*psc5 + (qiukp2)*dsc5);
#else
          ftm22                += (qiuk2 + qiukp2)*psc5;
#endif

    float qiuk3                 = qi3*atomJ.inducedDipole[0]  + qi6*atomJ.inducedDipole[1]  + qi9*atomJ.inducedDipole[2];
    float qiukp3                = qi3*atomJ.inducedDipoleP[0] + qi6*atomJ.inducedDipoleP[1] + qi9*atomJ.inducedDipoleP[2];
    float ftm23                 = -bn2*(qiuk3+qiukp3);
#ifdef APPLY_SCALE
          ftm23                += ((qiuk3)*psc5 + (qiukp3)*dsc5);
#else
          ftm23                += (qiuk3 + qiukp3)*psc5;
#endif

    float expdamp               = expf(damp);
    float scale3                = (damp == 0.0f) ? 1.0f : (1.0f - expdamp);
    float rr3                   = rr1*rr1*rr1;

#ifdef APPLY_SCALE
    float psc3                  = rr3*(1.0f - scale3*scalingFactors[PScaleIndex]);
    float dsc3                  = rr3*(1.0f - scale3*scalingFactors[DScaleIndex]);
    float usc3                  = rr3*(1.0f - scale3*scalingFactors[UScaleIndex]);
#else
    float psc3                  = rr3*(1.0f - scale3);
#endif

    float scale7                = (damp == 0.0f) ? 1.0f : (1.0f - (1.0f-damp+0.6f*damp*damp)*expdamp);

#ifdef APPLY_SCALE
    float psc7                  = (15.0f*rr3*rr3*rr1)*(1.0f - scale7*scalingFactors[PScaleIndex]);
    float dsc7                  = (15.0f*rr3*rr3*rr1)*(1.0f - scale7*scalingFactors[DScaleIndex]);
#else
    float psc7                  = (15.0f*rr3*rr3*rr1)*(1.0f - scale7);
#endif

    float qir1                  = qi1*xr + qi2*yr + qi3*zr;
    float qir2                  = qi2*xr + qi5*yr + qi6*zr;
    float qir3                  = qi3*xr + qi6*yr + qi9*zr;

    float sc3                   = di1*xr  + di2*yr  + di3*zr;
    float sc5                   = qir1*xr + qir2*yr + qir3*zr;
    float gfi3                  = ci*bn1  + sc3*bn2 + sc5*bn3;

    float prefactor1;
    prefactor1                  = 0.5f*(ci*psc3 + sc3*psc5 + sc5*psc7 - gfi3);
    ftm21                      -= prefactor1*atomJ.inducedDipole[0];
    ftm22                      -= prefactor1*atomJ.inducedDipole[1];
    ftm23                      -= prefactor1*atomJ.inducedDipole[2];

#ifdef APPLY_SCALE
    prefactor1                  = 0.5f*(ci*dsc3 + sc3*dsc5 + sc5*dsc7 - gfi3);
#endif
    ftm21                      -= prefactor1*atomJ.inducedDipoleP[0];
    ftm22                      -= prefactor1*atomJ.inducedDipoleP[1];
    ftm23                      -= prefactor1*atomJ.inducedDipoleP[2];

    float sci4                  = atomJ.inducedDipole[0]*xr  + atomJ.inducedDipole[1]*yr  + atomJ.inducedDipole[2]*zr;
    //forceTorqueEnergy->w       += 0.5f*((psc3-bn1)*(ci*sci4) + (psc5-bn2)*(sc3*sci4) + (psc7-bn3)*(sci4*sc5));
    *energy                    += forceFactor*0.5f*((psc3-bn1)*(ci*sci4) + (psc5-bn2)*(sc3*sci4) + (psc7-bn3)*(sci4*sc5));

    float scip4                 = atomJ.inducedDipoleP[0]*xr + atomJ.inducedDipoleP[1]*yr + atomJ.inducedDipoleP[2]*zr;
    if( cAmoebaSim.polarizationType == 0 ){

#ifdef APPLY_SCALE
        prefactor1              = 0.5f*( bn2 - usc5 );
#else
        prefactor1              = 0.5f*( bn2 - psc5 );
#endif
        ftm21                  += prefactor1*( (sci4*atomI.inducedDipoleP[0] + scip4*atomI.inducedDipole[0]) );
        ftm22                  += prefactor1*( (sci4*atomI.inducedDipoleP[1] + scip4*atomI.inducedDipole[1]) );
        ftm23                  += prefactor1*( (sci4*atomI.inducedDipoleP[2] + scip4*atomI.inducedDipole[2]) );
   }

#ifdef APPLY_SCALE
    prefactor1                  = 0.5f*( bn2*(sci4+scip4) - (sci4*psc5+scip4*dsc5) ); 
#else
    sci4                       += scip4;
    prefactor1                  = 0.5f*sci4*( bn2 - psc5 ); 
#endif

    ftm21                      += prefactor1*di1;
    ftm22                      += prefactor1*di2;
    ftm23                      += prefactor1*di3;

#ifdef APPLY_SCALE
    float gfi5                  = bn3*(sci4+scip4) - (sci4*psc7+scip4*dsc7);
#else
    float gfi5                  = sci4*(bn3 - psc7);
#endif
    ftm21                      += gfi5*qir1;
    ftm22                      += gfi5*qir2;
    ftm23                      += gfi5*qir3;

    float sci7                  = qir1*atomJ.inducedDipole[0]  + qir2*atomJ.inducedDipole[1]  + qir3*atomJ.inducedDipole[2];
    //forceTorqueEnergy->w       += (bn2-psc5)*sci7;
    *energy                    += forceFactor*(bn2-psc5)*sci7;
    float scip7                 = qir1*atomJ.inducedDipoleP[0] + qir2*atomJ.inducedDipoleP[1] + qir3*atomJ.inducedDipoleP[2];

#ifdef APPLY_SCALE
    float gli1                  = -ci*sci4;
    float gli2                  = -sc3*sci4 + 2.0f*sci7;
    float gli3                  = -sci4*sc5;
    float glip1                 = -ci*scip4;
    float glip2                 = -sc3*scip4 + 2.0f*scip7;
    float glip3                 = -scip4*sc5;
#else
    float gli1                  = -ci*sci4;
    float gli2                  = -sc3*sci4 + 2.0f*(sci7 + scip7);
    float gli3                  = -sci4*sc5;
#endif

#ifdef APPLY_SCALE
    float gfi1                  = (bn2*(gli1+glip1) + bn3*(gli2+glip2) + bn4*(gli3+glip3));
    gfi1                       -= (rr1*rr1)*( 3.0f*(gli1*psc3 + glip1*dsc3) + 5.0f*(gli2*psc5 + glip2*dsc5 ) + 7.0f*(gli3*psc7+glip3*dsc7) );
#else
    float gfi1                  = bn2*gli1 + bn3*gli2 + bn4*gli3;
    gfi1                       -= (rr1*rr1)*( 3.0f*gli1*psc3 + 5.0f*gli2*psc5 + 7.0f*gli3*psc7);
#endif
    gfi1                       *= 0.5f;
    ftm21                      += gfi1*xr;
    ftm22                      += gfi1*yr;
    ftm23                      += gfi1*zr;

    if( damp != 0.0f ){

        float expdamp = expf(damp);
        float temp3   = -1.5f*damp*expdamp*rr1*rr1;
        float temp5   = -damp;
        float temp7   = -0.2f - 0.6f*damp;

        float ddsc31  = temp3*xr;
        float ddsc32  = temp3*yr;
        float ddsc33  = temp3*zr;

        float ddsc51  = temp5*ddsc31;
        float ddsc52  = temp5*ddsc32;
        float ddsc53  = temp5*ddsc33;

        float ddsc71  = temp7*ddsc51;
        float ddsc72  = temp7*ddsc52;
        float ddsc73  = temp7*ddsc53;

        float rr3     = rr1*rr1*rr1;
#ifdef APPLY_SCALE
        temp3         = (gli1*scalingFactors[PScaleIndex] + glip1*scalingFactors[DScaleIndex]);
        temp5         = (3.0f*rr1*rr1)*(gli2*scalingFactors[PScaleIndex] + glip2*scalingFactors[DScaleIndex]);
        temp7         = (15.0f*rr3*rr1)*(gli3*scalingFactors[PScaleIndex] + glip3*scalingFactors[DScaleIndex]);
#else
        temp3         = gli1;
        temp5         = (3.0f*rr1*rr1)*gli2;
        temp7         = (15.0f*rr3*rr1)*gli3;
#endif
        ftm21        -= rr3*(temp3*ddsc31 + temp5*ddsc51 + temp7*ddsc71);
        ftm22        -= rr3*(temp3*ddsc32 + temp5*ddsc52 + temp7*ddsc72);
        ftm23        -= rr3*(temp3*ddsc33 + temp5*ddsc53 + temp7*ddsc73);
    }

//K
    float qk1                   = atomJ.labFrameQuadrupole[0];
    float qk2                   = atomJ.labFrameQuadrupole[1];
    float qk3                   = atomJ.labFrameQuadrupole[2];
    float qk5                   = atomJ.labFrameQuadrupole[3];
    float qk6                   = atomJ.labFrameQuadrupole[4];
    //float qk9                   = atomJ.labFrameQuadrupole[5];
    float qk9                   = -(qk1 + qk5);

    float qkui1                 = qk1*atomI.inducedDipole[0]  + qk2*atomI.inducedDipole[1]  + qk3*atomI.inducedDipole[2];
    float qkuip1                = qk1*atomI.inducedDipoleP[0] + qk2*atomI.inducedDipoleP[1] + qk3*atomI.inducedDipoleP[2];
          ftm21                += bn2*(qkui1+qkuip1);
#ifdef APPLY_SCALE
          ftm21                -= (qkui1*psc5 + qkuip1*dsc5);
#else
          ftm21                -= (qkui1 + qkuip1)*psc5;
#endif

    float qkui2                 = qk2*atomI.inducedDipole[0]  + qk5*atomI.inducedDipole[1]  + qk6*atomI.inducedDipole[2];
    float qkuip2                = qk2*atomI.inducedDipoleP[0] + qk5*atomI.inducedDipoleP[1] + qk6*atomI.inducedDipoleP[2];
          ftm22                += bn2*(qkui2+qkuip2);
#ifdef APPLY_SCALE
          ftm22                -= ((qkui2)*psc5 + (qkuip2)*dsc5);
#else
          ftm22                -= (qkui2 + qkuip2)*psc5;
#endif

    float qkui3                 = qk3*atomI.inducedDipole[0]  + qk6*atomI.inducedDipole[1]  + qk9*atomI.inducedDipole[2];
    float qkuip3                = qk3*atomI.inducedDipoleP[0] + qk6*atomI.inducedDipoleP[1] + qk9*atomI.inducedDipoleP[2];
          ftm23                += bn2*(qkui3+qkuip3);
#ifdef APPLY_SCALE
          ftm23                -= ((qkui3)*psc5 + (qkuip3)*dsc5);
#else
          ftm23                -= (qkui3 + qkuip3)*psc5;
#endif


    float qkr1                  = qk1*xr + qk2*yr + qk3*zr;
    float qkr2                  = qk2*xr + qk5*yr + qk6*zr;
    float qkr3                  = qk3*xr + qk6*yr + qk9*zr;

    float dk1                   = atomJ.labFrameDipole[0];
    float dk2                   = atomJ.labFrameDipole[1];
    float dk3                   = atomJ.labFrameDipole[2];

    float sc4                   =  dk1*xr  +  dk2*yr +  dk3*zr;
    float sc6                   = qkr1*xr  + qkr2*yr + qkr3*zr;

    float ck                    = atomJ.q;
    float gfi2                  = (-ck*bn1 + sc4*bn2 - sc6*bn3);

    prefactor1                  = 0.5f*(ck*psc3 - sc4*psc5 + sc6*psc7 + gfi2);
    ftm21                      += prefactor1*atomI.inducedDipole[0];
    ftm22                      += prefactor1*atomI.inducedDipole[1];
    ftm23                      += prefactor1*atomI.inducedDipole[2];

#ifdef APPLY_SCALE
    prefactor1                  = 0.5f*(ck*dsc3 - sc4*dsc5 + sc6*dsc7 + gfi2);
#endif
    ftm21                      += prefactor1*atomI.inducedDipoleP[0];
    ftm22                      += prefactor1*atomI.inducedDipoleP[1];
    ftm23                      += prefactor1*atomI.inducedDipoleP[2];

    float sci3                  = atomI.inducedDipole[0]*xr  + atomI.inducedDipole[1]*yr  + atomI.inducedDipole[2]*zr;
    //forceTorqueEnergy->w       += 0.5f*( (ck*sci3)*(bn1-psc3) -(sci3*sc4)*(bn2-psc5) + sci3*sc6*(bn3-psc7) );
    *energy                    += forceFactor*0.5f*( (ck*sci3)*(bn1-psc3) -(sci3*sc4)*(bn2-psc5) + sci3*sc6*(bn3-psc7) );
    float scip3                 = atomI.inducedDipoleP[0]*xr + atomI.inducedDipoleP[1]*yr + atomI.inducedDipoleP[2]*zr;

    if( cAmoebaSim.polarizationType == 0 ){
#ifdef APPLY_SCALE
        prefactor1              = 0.5f*( bn2 - usc5 );
#else
        prefactor1              = 0.5f*( bn2 - psc5 );
#endif

        ftm21                  += prefactor1*( sci3*atomJ.inducedDipoleP[0] + scip3*atomJ.inducedDipole[0] );
        ftm22                  += prefactor1*( sci3*atomJ.inducedDipoleP[1] + scip3*atomJ.inducedDipole[1] );
        ftm23                  += prefactor1*( sci3*atomJ.inducedDipoleP[2] + scip3*atomJ.inducedDipole[2] );
    }

    float sci34;
    if( cAmoebaSim.polarizationType == 0 ){
        float sci4              = atomJ.inducedDipole[0]*xr  + atomJ.inducedDipole[1]*yr  + atomJ.inducedDipole[2]*zr;
        float scip4             = atomJ.inducedDipoleP[0]*xr + atomJ.inducedDipoleP[1]*yr + atomJ.inducedDipoleP[2]*zr;
        sci34                   = (sci3*scip4+scip3*sci4);
   
#ifdef APPLY_SCALE
        gfi1                    = sci34*(usc5*(5.0f*rr1*rr1) -bn3 );
#else
        gfi1                    = sci34*(psc5*(5.0f*rr1*rr1) -bn3 );
#endif

    } else {
        gfi1                    = 0.0f;
    }

#ifdef APPLY_SCALE
    prefactor1                  = 0.5f*( bn2*(sci3+scip3) - (sci3*psc5+scip3*dsc5) );
#else
    sci3                       += scip3;
    prefactor1                  = 0.5f*sci3*( bn2 - psc5 );
#endif
    ftm21                      += prefactor1*dk1;
    ftm22                      += prefactor1*dk2;
    ftm23                      += prefactor1*dk3;

#ifdef APPLY_SCALE
    float gfi6                  = -bn3*(sci3+scip3) + (sci3*psc7+scip3*dsc7);
#else
    float gfi6                  = sci3*( psc7 - bn3);
#endif
    ftm21                      += gfi6*qkr1;
    ftm22                      += gfi6*qkr2;
    ftm23                      += gfi6*qkr3;

    float sci1                  = atomI.inducedDipole[0]*dk1 + atomI.inducedDipole[1]*dk2 + atomI.inducedDipole[2]*dk3 + di1*atomJ.inducedDipole[0] + di2*atomJ.inducedDipole[1] + di3*atomJ.inducedDipole[2];
    //forceTorqueEnergy->w       += 0.5f*( sci1*(bn1-psc3) );
    *energy                    += forceFactor*0.5f*( sci1*(bn1-psc3) );

    float sci8                  = qkr1*atomI.inducedDipole[0] + qkr2*atomI.inducedDipole[1] + qkr3*atomI.inducedDipole[2];
    //forceTorqueEnergy->w       += sci8*(bn2-psc5);
    *energy                    += forceFactor*sci8*(bn2-psc5);
    float scip1                 = atomI.inducedDipoleP[0]*dk1 + atomI.inducedDipoleP[1]*dk2 + atomI.inducedDipoleP[2]*dk3 + di1*atomJ.inducedDipoleP[0] + di2*atomJ.inducedDipoleP[1] + di3*atomJ.inducedDipoleP[2];
#ifndef APPLY_SCALE
        sci1                   += scip1;
#endif

    float scip2                 = atomI.inducedDipole[0]*atomJ.inducedDipoleP[0] +
                                  atomI.inducedDipole[1]*atomJ.inducedDipoleP[1] +
                                  atomI.inducedDipole[2]*atomJ.inducedDipoleP[2] +
                                  atomJ.inducedDipole[0]*atomI.inducedDipoleP[0] +
                                  atomJ.inducedDipole[1]*atomI.inducedDipoleP[1] +
                                  atomJ.inducedDipole[2]*atomI.inducedDipoleP[2];

    float scip8                 = qkr1*atomI.inducedDipoleP[0] + qkr2*atomI.inducedDipoleP[1] + qkr3*atomI.inducedDipoleP[2];
#ifndef APPLY_SCALE
          sci8                 += scip8;
#endif

           gli1                 = ck*sci3 + sci1;
           gli2                 = -(sci3*sc4 + 2.0f*sci8);
           gli3                 = sci3*sc6;
#ifdef APPLY_SCALE
          glip1                 = ck*scip3 + scip1;
          glip2                 = -(scip3*sc4 + 2.0f*scip8);
          glip3                 = scip3*sc6;
#endif


#ifdef APPLY_SCALE
    gfi1                       += (bn2*(gli1+glip1) + bn3*(gli2+glip2) + bn4*(gli3+glip3));
    gfi1                       -= (rr1*rr1)*( 3.0f*(gli1*psc3 + glip1*dsc3) + 5.0f*(gli2*psc5 + glip2*dsc5 ) + 7.0f*(gli3*psc7+glip3*dsc7) );
#else
    gfi1                       += (bn2*gli1 + bn3*gli2 + bn4*gli3);
    gfi1                       -= (rr1*rr1)*( 3.0f*gli1*psc3 + 5.0f*gli2*psc5 + 7.0f*gli3*psc7 );
#endif
    
    if( cAmoebaSim.polarizationType == 0 ){
#ifdef APPLY_SCALE
        gfi1                       += scip2*(bn2 - (3.0f*rr1*rr1)*usc3);
#else
        gfi1                       += scip2*(bn2 - (3.0f*rr1*rr1)*psc3);
#endif

    }
    gfi1                       *= 0.5f;

    ftm21                       += gfi1*xr;
    ftm22                       += gfi1*yr;
    ftm23                       += gfi1*zr;

    if( damp != 0.0f ){

        float expdamp = expf(damp);
        float temp3   = -1.5f*damp*expdamp*rr1*rr1;
        float temp5   = -damp;
        float temp7   = -0.2f - 0.6f*damp;

        float ddsc31  = temp3*xr;
        float ddsc32  = temp3*yr;
        float ddsc33  = temp3*zr;

        float ddsc51  = temp5*ddsc31;
        float ddsc52  = temp5*ddsc32;
        float ddsc53  = temp5*ddsc33;

        float ddsc71  = temp7*ddsc51;
        float ddsc72  = temp7*ddsc52;
        float ddsc73  = temp7*ddsc53;

        float rr3     = rr1*rr1*rr1;

#ifdef APPLY_SCALE
        temp3         =                  gli1*scalingFactors[PScaleIndex] + glip1*scalingFactors[DScaleIndex];
        temp5         = (3.0f*rr1*rr1)*( gli2*scalingFactors[PScaleIndex] + glip2*scalingFactors[DScaleIndex]);
        temp7         = (15.0f*rr3*rr1)*(gli3*scalingFactors[PScaleIndex] + glip3*scalingFactors[DScaleIndex]);
#else
        temp3         = gli1;
        temp5         = (3.0f*rr1*rr1)*gli2;
        temp7         = (15.0f*rr3*rr1)*(gli3);
#endif

        ftm21        -= rr3*(temp3*ddsc31 + temp5*ddsc51 + temp7*ddsc71);
        ftm22        -= rr3*(temp3*ddsc32 + temp5*ddsc52 + temp7*ddsc72);
        ftm23        -= rr3*(temp3*ddsc33 + temp5*ddsc53 + temp7*ddsc73);

        if( cAmoebaSim.polarizationType == 0 ){
#ifdef APPLY_SCALE
            temp3   =  scalingFactors[UScaleIndex]*scip2;
            temp5   = -(3.0f*rr1*rr1)*scalingFactors[UScaleIndex]*sci34;
#else
            temp3   =  scip2;
            temp5   = -(3.0f*rr1*rr1)*sci34;
#endif
            ftm21  -= rr3*(temp3*ddsc31 + temp5*ddsc51);
            ftm22  -= rr3*(temp3*ddsc32 + temp5*ddsc52);
            ftm23  -= rr3*(temp3*ddsc33 + temp5*ddsc53);
        }
    }

    force[0] += ftm21;
    force[1] += ftm22;
    force[2] += ftm23;
/*
    if( forceFactor == 1.0f ){
        atomJ.force[0]  -= ftm21;
        atomJ.force[1]  -= ftm22;
        atomJ.force[2]  -= ftm23;
    }   
    atomI.force[0]      += ftm21;
    atomI.force[1]      += ftm22;
    atomI.force[2]      += ftm23;
*/
/*
    forceTorqueEnergy->x += ftm21;
    forceTorqueEnergy->y += ftm22;
    forceTorqueEnergy->z += ftm23;
*/

    return;

}
