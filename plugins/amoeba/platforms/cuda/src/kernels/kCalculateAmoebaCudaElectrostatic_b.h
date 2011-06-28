
__device__ void SUB_METHOD_NAME( calculateElectrostaticPairIxn, _kernel )( ElectrostaticParticle& atomI,   ElectrostaticParticle& atomJ,
                                                                           float* scalingFactors,
#ifdef F1
                                                                           float* energy,
#endif
                                                                           float* outputForce ){
  
#ifdef F1
    float ddsc3_0            = 0.0f;
    float ddsc3_1            = 0.0f;
    float ddsc3_2            = 0.0f;

    float ddsc5_0            = 0.0f;
    float ddsc5_1            = 0.0f;
    float ddsc5_2            = 0.0f;

    float ddsc7_0            = 0.0f;
    float ddsc7_1            = 0.0f;
    float ddsc7_2            = 0.0f;
#endif

    float xr                 = atomJ.x - atomI.x;
    float yr                 = atomJ.y - atomI.y;
    float zr                 = atomJ.z - atomI.z;
   
    float r2                 = xr*xr + yr*yr + zr*zr;
    float r                  = sqrtf( r2 );
    float rr1                = 1.0f/r;
    float rr2                = rr1*rr1;
    float rr3                = rr1*rr2;
    float rr5                = 3.0f*rr3*rr2;
    float rr7                = 5.0f*rr5*rr2;
    float rr9                = 7.0f*rr7*rr2;
#ifdef F1
    float rr11               = 9.0f*rr9*rr2;
#endif

    float scale3             = 1.0f;
    float scale5             = 1.0f;
    float scale7             = 1.0f;

    float pdamp              = atomI.damp*atomJ.damp;
    if( pdamp != 0.0 && r < cAmoebaSim.scalingDistanceCutoff ){
   
        float ratio                   = r/pdamp;
        float pGamma                  = atomJ.thole > atomI.thole ? atomI.thole : atomJ.thole;

        float damp                    = ratio*ratio*ratio*pGamma;
        float dampExp                 = expf( -damp );
        float damp1                   = damp + 1.0f;
        float damp2                   = damp*damp;

        scale3                        = 1.0f - dampExp;
        scale5                        = 1.0f - damp1*dampExp;
        scale7                        = 1.0f - ( damp1 + 0.6f*damp2)*dampExp;

#ifdef F1
        float factor                  = 3.0f*damp*dampExp*rr2;
        float factor7                 = -0.2f + 0.6f*damp;
        
        ddsc3_0                       = factor*xr;
        ddsc5_0                       = ddsc3_0*damp;
        ddsc7_0                       = ddsc5_0*factor7;

        ddsc3_1                       = factor*yr;
        ddsc5_1                       = ddsc3_1*damp;
        ddsc7_1                       = ddsc5_1*factor7;

        ddsc3_2                       = factor*zr;
        ddsc5_2                       = ddsc3_2*damp;
        ddsc7_2                       = ddsc5_2*factor7;
#endif

    }
      
#if defined F1
    float scale3i            = rr3*scale3*scalingFactors[UScaleIndex];
    float scale5i            = rr5*scale5*scalingFactors[UScaleIndex];
#endif
    float dsc3               = rr3*scale3*scalingFactors[DScaleIndex];
    float psc3               = rr3*scale3*scalingFactors[PScaleIndex];

    float dsc5               = rr5*scale5*scalingFactors[DScaleIndex];
    float psc5               = rr5*scale5*scalingFactors[PScaleIndex];

    float dsc7               = rr7*scale7*scalingFactors[DScaleIndex];
    float psc7               = rr7*scale7*scalingFactors[PScaleIndex];
                       
    float qJr_0              = atomJ.labFrameQuadrupole[0]*xr + atomJ.labFrameQuadrupole[3]*yr + atomJ.labFrameQuadrupole[6]*zr;
    float qJr_1              = atomJ.labFrameQuadrupole[1]*xr + atomJ.labFrameQuadrupole[4]*yr + atomJ.labFrameQuadrupole[7]*zr;
    float qJr_2              = atomJ.labFrameQuadrupole[2]*xr + atomJ.labFrameQuadrupole[5]*yr + atomJ.labFrameQuadrupole[8]*zr;

    float qIr_0              = atomI.labFrameQuadrupole[0]*xr + atomI.labFrameQuadrupole[3]*yr + atomI.labFrameQuadrupole[6]*zr;
    float qIr_1              = atomI.labFrameQuadrupole[1]*xr + atomI.labFrameQuadrupole[4]*yr + atomI.labFrameQuadrupole[7]*zr;
    float qIr_2              = atomI.labFrameQuadrupole[2]*xr + atomI.labFrameQuadrupole[5]*yr + atomI.labFrameQuadrupole[8]*zr;

#if defined F1
    float sc2                = atomI.labFrameDipole[0]*atomJ.labFrameDipole[0] + atomI.labFrameDipole[1]*atomJ.labFrameDipole[1] + atomI.labFrameDipole[2]*atomJ.labFrameDipole[2];
#endif
#if defined F1 || defined T1
    float sc4                = atomJ.labFrameDipole[0]*xr + atomJ.labFrameDipole[1]*yr + atomJ.labFrameDipole[2]*zr;
    float sc6                = qJr_0*xr + qJr_1*yr + qJr_2*zr;
#endif

#if defined F1 || defined T3
    float sc3                = atomI.labFrameDipole[0]*xr + atomI.labFrameDipole[1]*yr + atomI.labFrameDipole[2]*zr;
    float sc5                = qIr_0*xr + qIr_1*yr + qIr_2*zr;
#endif
    
#if defined F1
    float sc7                = qIr_0*atomJ.labFrameDipole[0] + qIr_1*atomJ.labFrameDipole[1] + qIr_2*atomJ.labFrameDipole[2];
    float sc8                = qJr_0*atomI.labFrameDipole[0] + qJr_1*atomI.labFrameDipole[1] + qJr_2*atomI.labFrameDipole[2];
    float sc9                = qIr_0*qJr_0 + qIr_1*qJr_1 + qIr_2*qJr_2;

	 float sc10               = atomI.labFrameQuadrupole[0]*atomJ.labFrameQuadrupole[0] + atomI.labFrameQuadrupole[1]*atomJ.labFrameQuadrupole[1] + atomI.labFrameQuadrupole[2]*atomJ.labFrameQuadrupole[2] +
                               atomI.labFrameQuadrupole[3]*atomJ.labFrameQuadrupole[3] + atomI.labFrameQuadrupole[4]*atomJ.labFrameQuadrupole[4] + atomI.labFrameQuadrupole[5]*atomJ.labFrameQuadrupole[5] +
                               atomI.labFrameQuadrupole[6]*atomJ.labFrameQuadrupole[6] + atomI.labFrameQuadrupole[7]*atomJ.labFrameQuadrupole[7] + atomI.labFrameQuadrupole[8]*atomJ.labFrameQuadrupole[8];

    float sci1               = atomI.inducedDipole[0]*atomJ.labFrameDipole[0] + atomI.inducedDipole[1]*atomJ.labFrameDipole[1] + atomI.inducedDipole[2]*atomJ.labFrameDipole[2] +
                               atomJ.inducedDipole[0]*atomI.labFrameDipole[0] + atomJ.inducedDipole[1]*atomI.labFrameDipole[1] + atomJ.inducedDipole[2]*atomI.labFrameDipole[2];
#endif
        
#if defined F1 || defined T3
    float sci3               = atomI.inducedDipole[0]*xr + atomI.inducedDipole[1]*yr + atomI.inducedDipole[2]*zr;
#endif
#if defined F1
    float sci7               = qIr_0*atomJ.inducedDipole[0] + qIr_1*atomJ.inducedDipole[1] + qIr_2*atomJ.inducedDipole[2];
    float sci8               = qJr_0*atomI.inducedDipole[0] + qJr_1*atomI.inducedDipole[1] + qJr_2*atomI.inducedDipole[2];
#endif
#if defined F1 || defined T1
    float sci4               = atomJ.inducedDipole[0]*xr + atomJ.inducedDipole[1]*yr + atomJ.inducedDipole[2]*zr;
#endif
    
#if defined F1
    float scip1              = atomI.inducedDipoleP[0]*atomJ.labFrameDipole[0] + atomI.inducedDipoleP[1]*atomJ.labFrameDipole[1] + atomI.inducedDipoleP[2]*atomJ.labFrameDipole[2] +
                               atomJ.inducedDipoleP[0]*atomI.labFrameDipole[0] + atomJ.inducedDipoleP[1]*atomI.labFrameDipole[1] + atomJ.inducedDipoleP[2]*atomI.labFrameDipole[2];

    float scip2              = atomI.inducedDipole[0]*atomJ.inducedDipoleP[0] + atomI.inducedDipole[1]*atomJ.inducedDipoleP[1] + atomI.inducedDipole[2]*atomJ.inducedDipoleP[2] +
                               atomJ.inducedDipole[0]*atomI.inducedDipoleP[0] + atomJ.inducedDipole[1]*atomI.inducedDipoleP[1] + atomJ.inducedDipole[2]*atomI.inducedDipoleP[2];

#endif
    
#if defined F1 || defined T3
    float scip3              = ((atomI.inducedDipoleP[0])*(xr) + (atomI.inducedDipoleP[1])*(yr) + (atomI.inducedDipoleP[2])*(zr));
#endif
#if defined F1 || defined T1
    float scip4              = ((atomJ.inducedDipoleP[0])*(xr) + (atomJ.inducedDipoleP[1])*(yr) + (atomJ.inducedDipoleP[2])*(zr));
#endif

#ifdef F1
    float scip7              = ((qIr_0)*(atomJ.inducedDipoleP[0]) + (qIr_1)*(atomJ.inducedDipoleP[1]) + (qIr_2)*(atomJ.inducedDipoleP[2]));
    float scip8              = ((qJr_0)*(atomI.inducedDipoleP[0]) + (qJr_1)*(atomI.inducedDipoleP[1]) + (qJr_2)*(atomI.inducedDipoleP[2]));



    float gli1               = atomJ.q*sci3 - atomI.q*sci4;
    
    float gli6               = sci1;
    float glip1              = atomJ.q*scip3 - atomI.q*scip4;
    float glip6              = scip1;
    float gli2               = -sc3*sci4 - sci3*sc4;
    float gli3               = sci3*sc6 - sci4*sc5;
    float gli7               = 2.0f*(sci7-sci8);
    
    float glip2              = -sc3*scip4 - scip3*sc4;
    float glip3              = scip3*sc6 - scip4*sc5;
    float glip7              = 2.0f*(scip7-scip8);
    float factor3            = rr3*(( gli1  +  gli6)*scalingFactors[PScaleIndex] + (glip1  + glip6)*scalingFactors[DScaleIndex]);
    float factor5            = rr5*(( gli2  +  gli7)*scalingFactors[PScaleIndex] + (glip2  + glip7)*scalingFactors[DScaleIndex]);
    float factor7            = rr7*( gli3*scalingFactors[PScaleIndex] + glip3*scalingFactors[DScaleIndex]);
      
    float ftm2i_0            = -0.5f*(factor3*ddsc3_0 + factor5*ddsc5_0 + factor7*ddsc7_0);
    float ftm2i_1            = -0.5f*(factor3*ddsc3_1 + factor5*ddsc5_1 + factor7*ddsc7_1);
    float ftm2i_2            = -0.5f*(factor3*ddsc3_2 + factor5*ddsc5_2 + factor7*ddsc7_2);
      
    float gl0                = atomI.q*atomJ.q;
    float gl1                = atomJ.q*sc3 - atomI.q*sc4;
    float gl2                = atomI.q*sc6 + atomJ.q*sc5 - sc3*sc4;
    float gl3                = sc3*sc6 - sc4*sc5;
    float gl4                = sc5*sc6;
    float gl6                = sc2;
    float gl7                = 2.0f*(sc7-sc8);
    float gl8                = 2.0f*sc10;
    float gl5                = -4.0f*sc9;
    
    float gf1                = rr3*gl0 + rr5*(gl1+gl6) + rr7*(gl2+gl7+gl8) + rr9*(gl3+gl5) + rr11*gl4;
#endif
#if defined F1 || defined T1
    float gf2                = -atomJ.q*rr3 + sc4*rr5 - sc6*rr7;
    float gf5                = 2.0f*(-atomJ.q*rr5+sc4*rr7-sc6*rr9);
#endif
#if defined F1 || defined T3
    float gf3                =  atomI.q*rr3 + sc3*rr5 + sc5*rr7;
    float gf6                = 2.0f*(-atomI.q*rr5-sc3*rr7-sc5*rr9);
#endif

#ifdef F1
    float em                 = scalingFactors[MScaleIndex]*(rr1*gl0 + rr3*(gl1+gl6) + rr5*(gl2+gl7+gl8) + rr7*(gl3+gl5) + rr9*gl4);
    float ei                 = 0.5f*((gli1+gli6)*psc3 + (gli2+gli7)*psc5 + gli3*psc7);
    *energy                  = em+ei;
#endif
    
#if defined F1 || defined T1

    float qIdJ_0 = atomI.labFrameQuadrupole[0]*atomJ.labFrameDipole[0] + atomI.labFrameQuadrupole[3]*atomJ.labFrameDipole[1] + atomI.labFrameQuadrupole[6]*atomJ.labFrameDipole[2];
    float qIdJ_1 = atomI.labFrameQuadrupole[1]*atomJ.labFrameDipole[0] + atomI.labFrameQuadrupole[4]*atomJ.labFrameDipole[1] + atomI.labFrameQuadrupole[7]*atomJ.labFrameDipole[2];
    float qIdJ_2 = atomI.labFrameQuadrupole[2]*atomJ.labFrameDipole[0] + atomI.labFrameQuadrupole[5]*atomJ.labFrameDipole[1] + atomI.labFrameQuadrupole[8]*atomJ.labFrameDipole[2];

    float qIqJr_0 = atomI.labFrameQuadrupole[0]*qJr_0 + atomI.labFrameQuadrupole[3]*qJr_1 + atomI.labFrameQuadrupole[6]*qJr_2;
    float qIqJr_1 = atomI.labFrameQuadrupole[1]*qJr_0 + atomI.labFrameQuadrupole[4]*qJr_1 + atomI.labFrameQuadrupole[7]*qJr_2;
    float qIqJr_2 = atomI.labFrameQuadrupole[2]*qJr_0 + atomI.labFrameQuadrupole[5]*qJr_1 + atomI.labFrameQuadrupole[8]*qJr_2;
#endif

#ifdef F1
    float qkqir_0 = atomJ.labFrameQuadrupole[0]*qIr_0 + atomJ.labFrameQuadrupole[3]*qIr_1 + atomJ.labFrameQuadrupole[6]*qIr_2;
    float qkqir_1 = atomJ.labFrameQuadrupole[1]*qIr_0 + atomJ.labFrameQuadrupole[4]*qIr_1 + atomJ.labFrameQuadrupole[7]*qIr_2;
    float qkqir_2 = atomJ.labFrameQuadrupole[2]*qIr_0 + atomJ.labFrameQuadrupole[5]*qIr_1 + atomJ.labFrameQuadrupole[8]*qIr_2;

    float qkdi_0 = atomJ.labFrameQuadrupole[0]*atomI.labFrameDipole[0] + atomJ.labFrameQuadrupole[3]*atomI.labFrameDipole[1] + atomJ.labFrameQuadrupole[6]*atomI.labFrameDipole[2];
    float qkdi_1 = atomJ.labFrameQuadrupole[1]*atomI.labFrameDipole[0] + atomJ.labFrameQuadrupole[4]*atomI.labFrameDipole[1] + atomJ.labFrameQuadrupole[7]*atomI.labFrameDipole[2];
    float qkdi_2 = atomJ.labFrameQuadrupole[2]*atomI.labFrameDipole[0] + atomJ.labFrameQuadrupole[5]*atomI.labFrameDipole[1] + atomJ.labFrameQuadrupole[8]*atomI.labFrameDipole[2];

    float ftm2_0   = scalingFactors[MScaleIndex]*(gf1*xr + gf2*atomI.labFrameDipole[0] + gf3*atomJ.labFrameDipole[0] + 2.0f*rr5*(qkdi_0 - qIdJ_0) + gf5*qIr_0 + gf6*qJr_0 + 4.0f*rr7*(qIqJr_0 + qkqir_0));
    float ftm2_1   = scalingFactors[MScaleIndex]*(gf1*yr + gf2*atomI.labFrameDipole[1] + gf3*atomJ.labFrameDipole[1] + 2.0f*rr5*(qkdi_1 - qIdJ_1) + gf5*qIr_1 + gf6*qJr_1 + 4.0f*rr7*(qIqJr_1 + qkqir_1));
    float ftm2_2   = scalingFactors[MScaleIndex]*(gf1*zr + gf2*atomI.labFrameDipole[2] + gf3*atomJ.labFrameDipole[2] + 2.0f*rr5*(qkdi_2 - qIdJ_2) + gf5*qIr_2 + gf6*qJr_2 + 4.0f*rr7*(qIqJr_2 + qkqir_2));

    float gfi1 = rr2*(1.5f*((gli1+gli6)*psc3 + (glip1+glip6)*dsc3 + scip2*scale3i) + 2.5f*((gli7+gli2)*psc5 + (glip7+glip2)*dsc5 - (sci3*scip4+scip3*sci4)*scale5i) + 3.5f*(gli3*psc7+glip3*dsc7));
    ftm2i_0   += gfi1*xr;
    ftm2i_1   += gfi1*yr;
    ftm2i_2   += gfi1*zr;
#endif

#if defined F1 || defined T1
    float gfi5 =  (sci4*psc7 + scip4*dsc7);
#endif
#if defined F1 || defined T3
    float gfi6 = -(sci3*psc7 + scip3*dsc7);
#endif

#if defined F1 || defined T1
    float qIuJ_0 = atomI.labFrameQuadrupole[0]*atomJ.inducedDipole[0]   + atomI.labFrameQuadrupole[3]*atomJ.inducedDipole[1]  + atomI.labFrameQuadrupole[6]*atomJ.inducedDipole[2];
    float qIuJ_1 = atomI.labFrameQuadrupole[1]*atomJ.inducedDipole[0]   + atomI.labFrameQuadrupole[4]*atomJ.inducedDipole[1]  + atomI.labFrameQuadrupole[7]*atomJ.inducedDipole[2];
    float qIuJ_2 = atomI.labFrameQuadrupole[2]*atomJ.inducedDipole[0]   + atomI.labFrameQuadrupole[5]*atomJ.inducedDipole[1]  + atomI.labFrameQuadrupole[8]*atomJ.inducedDipole[2];

    float qIuJp_0 = atomI.labFrameQuadrupole[0]*atomJ.inducedDipoleP[0] + atomI.labFrameQuadrupole[3]*atomJ.inducedDipoleP[1] + atomI.labFrameQuadrupole[6]*atomJ.inducedDipoleP[2];
    float qIuJp_1 = atomI.labFrameQuadrupole[1]*atomJ.inducedDipoleP[0] + atomI.labFrameQuadrupole[4]*atomJ.inducedDipoleP[1] + atomI.labFrameQuadrupole[7]*atomJ.inducedDipoleP[2];
    float qIuJp_2 = atomI.labFrameQuadrupole[2]*atomJ.inducedDipoleP[0] + atomI.labFrameQuadrupole[5]*atomJ.inducedDipoleP[1] + atomI.labFrameQuadrupole[8]*atomJ.inducedDipoleP[2];
#endif

#if defined T3
    float qJuIp_0 = atomJ.labFrameQuadrupole[0]*atomI.inducedDipoleP[0] + atomJ.labFrameQuadrupole[3]*atomI.inducedDipoleP[1] + atomJ.labFrameQuadrupole[6]*atomI.inducedDipoleP[2];
    float qJuIp_1 = atomJ.labFrameQuadrupole[1]*atomI.inducedDipoleP[0] + atomJ.labFrameQuadrupole[4]*atomI.inducedDipoleP[1] + atomJ.labFrameQuadrupole[7]*atomI.inducedDipoleP[2];
    float qJuIp_2 = atomJ.labFrameQuadrupole[2]*atomI.inducedDipoleP[0] + atomJ.labFrameQuadrupole[5]*atomI.inducedDipoleP[1] + atomJ.labFrameQuadrupole[8]*atomI.inducedDipoleP[2];

     float qJuI_0 = atomJ.labFrameQuadrupole[0]*atomI.inducedDipole[0] + atomJ.labFrameQuadrupole[3]*atomI.inducedDipole[1] + atomJ.labFrameQuadrupole[6]*atomI.inducedDipole[2];
     float qJuI_1 = atomJ.labFrameQuadrupole[1]*atomI.inducedDipole[0] + atomJ.labFrameQuadrupole[4]*atomI.inducedDipole[1] + atomJ.labFrameQuadrupole[7]*atomI.inducedDipole[2];
     float qJuI_2 = atomJ.labFrameQuadrupole[2]*atomI.inducedDipole[0] + atomJ.labFrameQuadrupole[5]*atomI.inducedDipole[1] + atomJ.labFrameQuadrupole[8]*atomI.inducedDipole[2];
#endif

#ifdef F1

    float qkui_0 = atomJ.labFrameQuadrupole[0]*atomI.inducedDipole[0] + atomJ.labFrameQuadrupole[3]*atomI.inducedDipole[1] + atomJ.labFrameQuadrupole[6]*atomI.inducedDipole[2];
    float qkui_1 = atomJ.labFrameQuadrupole[1]*atomI.inducedDipole[0] + atomJ.labFrameQuadrupole[4]*atomI.inducedDipole[1] + atomJ.labFrameQuadrupole[7]*atomI.inducedDipole[2];
    float qkui_2 = atomJ.labFrameQuadrupole[2]*atomI.inducedDipole[0] + atomJ.labFrameQuadrupole[5]*atomI.inducedDipole[1] + atomJ.labFrameQuadrupole[8]*atomI.inducedDipole[2];

    float qkuip_0 = atomJ.labFrameQuadrupole[0]*atomI.inducedDipoleP[0] + atomJ.labFrameQuadrupole[3]*atomI.inducedDipoleP[1] + atomJ.labFrameQuadrupole[6]*atomI.inducedDipoleP[2];
    float qkuip_1 = atomJ.labFrameQuadrupole[1]*atomI.inducedDipoleP[0] + atomJ.labFrameQuadrupole[4]*atomI.inducedDipoleP[1] + atomJ.labFrameQuadrupole[7]*atomI.inducedDipoleP[2];
    float qkuip_2 = atomJ.labFrameQuadrupole[2]*atomI.inducedDipoleP[0] + atomJ.labFrameQuadrupole[5]*atomI.inducedDipoleP[1] + atomJ.labFrameQuadrupole[8]*atomI.inducedDipoleP[2];

    ftm2i_0 += 0.5f*(-atomJ.q*(atomI.inducedDipole[0]*psc3 + atomI.inducedDipoleP[0]*dsc3) +
                    sc4*(atomI.inducedDipole[0]*psc5 + atomI.inducedDipoleP[0]*dsc5) -
                    sc6*(atomI.inducedDipole[0]*psc7 + atomI.inducedDipoleP[0]*dsc7)) +
      
                   0.5f*(atomI.q*(atomJ.inducedDipole[0]*psc3+atomJ.inducedDipoleP[0]*dsc3) +
                     sc3*(atomJ.inducedDipole[0]*psc5 +atomJ.inducedDipoleP[0]*dsc5) +
                     sc5*(atomJ.inducedDipole[0]*psc7 +atomJ.inducedDipoleP[0]*dsc7)) +

                     scale5i*(sci4*atomI.inducedDipoleP[0]+scip4*atomI.inducedDipole[0] +
                     sci3*atomJ.inducedDipoleP[0]+scip3*atomJ.inducedDipole[0])*0.5f +
      
                    0.5f*(sci4*psc5+scip4*dsc5)*atomI.labFrameDipole[0] +
                    0.5f*(sci3*psc5+scip3*dsc5)*atomJ.labFrameDipole[0] +
                    ((qkui_0-qIuJ_0)*psc5 + (qkuip_0-qIuJp_0)*dsc5) +
                    gfi5*qIr_0 + gfi6*qJr_0;
      
    ftm2i_1 += 0.5f*(-atomJ.q*(atomI.inducedDipole[1]*psc3 + atomI.inducedDipoleP[1]*dsc3) +
                    sc4*(atomI.inducedDipole[1]*psc5 + atomI.inducedDipoleP[1]*dsc5) -
                    sc6*(atomI.inducedDipole[1]*psc7 + atomI.inducedDipoleP[1]*dsc7)) +
      
                    (atomI.q*(atomJ.inducedDipole[1]*psc3+atomJ.inducedDipoleP[1]*dsc3) +
                     sc3*(atomJ.inducedDipole[1]*psc5 +atomJ.inducedDipoleP[1]*dsc5) +
                     sc5*(atomJ.inducedDipole[1]*psc7 +atomJ.inducedDipoleP[1]*dsc7))*0.5f +
                     scale5i*(sci4*atomI.inducedDipoleP[1]+scip4*atomI.inducedDipole[1] + sci3*atomJ.inducedDipoleP[1]+scip3*atomJ.inducedDipole[1])*0.5f +
      
                    0.5f*(sci4*psc5+scip4*dsc5)*atomI.labFrameDipole[1] +
                    0.5f*(sci3*psc5+scip3*dsc5)*atomJ.labFrameDipole[1] +
                    ((qkui_1-qIuJ_1)*psc5 + (qkuip_1-qIuJp_1)*dsc5) +
                    gfi5*qIr_1 + gfi6*qJr_1;
      
    ftm2i_2 += 0.5f*(-atomJ.q*(atomI.inducedDipole[2]*psc3 + atomI.inducedDipoleP[2]*dsc3) +
                    sc4*(atomI.inducedDipole[2]*psc5 + atomI.inducedDipoleP[2]*dsc5) -
                    sc6*(atomI.inducedDipole[2]*psc7 + atomI.inducedDipoleP[2]*dsc7)) +
      
                    (atomI.q*(atomJ.inducedDipole[2]*psc3+atomJ.inducedDipoleP[2]*dsc3) +
                     sc3*(atomJ.inducedDipole[2]*psc5 +atomJ.inducedDipoleP[2]*dsc5) +
                     sc5*(atomJ.inducedDipole[2]*psc7 +atomJ.inducedDipoleP[2]*dsc7))*0.5f +
                     scale5i*(sci4*atomI.inducedDipoleP[2]+scip4*atomI.inducedDipole[2] +
                     sci3*atomJ.inducedDipoleP[2]+scip3*atomJ.inducedDipole[2])*0.5f +
      
                    0.5f*(sci4*psc5+scip4*dsc5)*atomI.labFrameDipole[2] +
                    0.5f*(sci3*psc5+scip3*dsc5)*atomJ.labFrameDipole[2] +
                    ((qkui_2-qIuJ_2)*psc5 + (qkuip_2-qIuJp_2)*dsc5) +
                    gfi5*qIr_2 + gfi6*qJr_2;

    if( cAmoebaSim.polarizationType )
    {
        float gfd                 = 0.5*(3.0*rr2*scip2*scale3i - 5.0f*rr2*(scip3*sci4+sci3*scip4)*scale5i);
        float temp5               = 0.5*scale5i;
        float fdir_0              = gfd*xr + temp5*(sci4*atomI.inducedDipoleP[0] + scip4*atomI.inducedDipole[0] + sci3*atomJ.inducedDipoleP[0] + scip3*atomJ.inducedDipole[0]);
        float fdir_1              = gfd*yr + temp5*(sci4*atomI.inducedDipoleP[1] + scip4*atomI.inducedDipole[1] + sci3*atomJ.inducedDipoleP[1] + scip3*atomJ.inducedDipole[1]);
        float fdir_2              = gfd*zr + temp5*(sci4*atomI.inducedDipoleP[2] + scip4*atomI.inducedDipole[2] + sci3*atomJ.inducedDipoleP[2] + scip3*atomJ.inducedDipole[2]);
        ftm2i_0                  -= fdir_0;
        ftm2i_1                  -= fdir_1;
        ftm2i_2                  -= fdir_2;

    } else {
        float scaleF              = 0.5f*scalingFactors[UScaleIndex];
        float inducedFactor3      = scip2*rr3*scaleF;
        float inducedFactor5      = (sci3*scip4+scip3*sci4)*rr5*scaleF;
        float findmp_0            = inducedFactor3*ddsc3_0 - inducedFactor5*ddsc5_0;
        float findmp_1            = inducedFactor3*ddsc3_1 - inducedFactor5*ddsc5_1;
        float findmp_2            = inducedFactor3*ddsc3_2 - inducedFactor5*ddsc5_2;
        ftm2i_0                  -= findmp_0;
        ftm2i_1                  -= findmp_1;
        ftm2i_2                  -= findmp_2;
    }
#endif

#if defined T1
    float gti2 = 0.5f*(sci4*psc5+scip4*dsc5);
    float gti5 = gfi5;
#endif
#if defined T3
    float gti3 = 0.5f*(sci3*psc5+scip3*dsc5);
    float gti6 = gfi6;
#endif

#if defined T1 || defined T3
    float dixdk_0 = atomI.labFrameDipole[1]*atomJ.labFrameDipole[2] - atomI.labFrameDipole[2]*atomJ.labFrameDipole[1];
    float dixdk_1 = atomI.labFrameDipole[2]*atomJ.labFrameDipole[0] - atomI.labFrameDipole[0]*atomJ.labFrameDipole[2];
    float dixdk_2 = atomI.labFrameDipole[0]*atomJ.labFrameDipole[1] - atomI.labFrameDipole[1]*atomJ.labFrameDipole[0];

#if defined T1
    float dixuk_0 = atomI.labFrameDipole[1]*atomJ.inducedDipole[2] - atomI.labFrameDipole[2]*atomJ.inducedDipole[1];
    float dixuk_1 = atomI.labFrameDipole[2]*atomJ.inducedDipole[0] - atomI.labFrameDipole[0]*atomJ.inducedDipole[2];
    float dixuk_2 = atomI.labFrameDipole[0]*atomJ.inducedDipole[1] - atomI.labFrameDipole[1]*atomJ.inducedDipole[0];
#endif
#endif

#ifdef T1
    float dixukp_0 = atomI.labFrameDipole[1]*atomJ.inducedDipoleP[2] - atomI.labFrameDipole[2]*atomJ.inducedDipoleP[1];
    float dixukp_1 = atomI.labFrameDipole[2]*atomJ.inducedDipoleP[0] - atomI.labFrameDipole[0]*atomJ.inducedDipoleP[2];
    float dixukp_2 = atomI.labFrameDipole[0]*atomJ.inducedDipoleP[1] - atomI.labFrameDipole[1]*atomJ.inducedDipoleP[0];
#endif

#ifdef T1
    float dixr_0 = atomI.labFrameDipole[1]*zr - atomI.labFrameDipole[2]*yr;
    float dixr_1 = atomI.labFrameDipole[2]*xr - atomI.labFrameDipole[0]*zr;
    float dixr_2 = atomI.labFrameDipole[0]*yr - atomI.labFrameDipole[1]*xr;
#endif

#ifdef T1
    float rxqiukp_0 = yr*qIuJp_2 - zr*qIuJp_1;
    float rxqiukp_1 = zr*qIuJp_0 - xr*qIuJp_2;
    float rxqiukp_2 = xr*qIuJp_1 - yr*qIuJp_0;

    float rxqir_0   = yr*qIr_2 - zr*qIr_1;
    float rxqir_1   = zr*qIr_0 - xr*qIr_2;
    float rxqir_2   = xr*qIr_1 - yr*qIr_0;

    float rxqiuk_0 = yr*qIuJ_2 - zr*qIuJ_1;
    float rxqiuk_1 = zr*qIuJ_0 - xr*qIuJ_2;
    float rxqiuk_2 = xr*qIuJ_1 - yr*qIuJ_0;

    float ukxqir_0 = atomJ.inducedDipole[1]*qIr_2 - atomJ.inducedDipole[2]*qIr_1;
    float ukxqir_1 = atomJ.inducedDipole[2]*qIr_0 - atomJ.inducedDipole[0]*qIr_2;
    float ukxqir_2 = atomJ.inducedDipole[0]*qIr_1 - atomJ.inducedDipole[1]*qIr_0;

    float ukxqirp_0 = atomJ.inducedDipoleP[1]*qIr_2 - atomJ.inducedDipoleP[2]*qIr_1;
    float ukxqirp_1 = atomJ.inducedDipoleP[2]*qIr_0 - atomJ.inducedDipoleP[0]*qIr_2;
    float ukxqirp_2 = atomJ.inducedDipoleP[0]*qIr_1 - atomJ.inducedDipoleP[1]*qIr_0;

    float dixqkr_0 = atomI.labFrameDipole[1]*qJr_2 - atomI.labFrameDipole[2]*qJr_1;
    float dixqkr_1 = atomI.labFrameDipole[2]*qJr_0 - atomI.labFrameDipole[0]*qJr_2;
    float dixqkr_2 = atomI.labFrameDipole[0]*qJr_1 - atomI.labFrameDipole[1]*qJr_0;

    float dkxqir_0 = atomJ.labFrameDipole[1]*qIr_2 - atomJ.labFrameDipole[2]*qIr_1;
    float dkxqir_1 = atomJ.labFrameDipole[2]*qIr_0 - atomJ.labFrameDipole[0]*qIr_2;
    float dkxqir_2 = atomJ.labFrameDipole[0]*qIr_1 - atomJ.labFrameDipole[1]*qIr_0;

    float rxqikr_0 = yr*qIqJr_2 - zr*qIqJr_1;
    float rxqikr_1 = zr*qIqJr_0 - xr*qIqJr_2;
    float rxqikr_2 = xr*qIqJr_1 - yr*qIqJr_0;

    float rxqidk_0 = yr*qIdJ_2 - zr*qIdJ_1;
    float rxqidk_1 = zr*qIdJ_0 - xr*qIdJ_2;
    float rxqidk_2 = xr*qIdJ_1 - yr*qIdJ_0;

    float qkrxqir_0 = qJr_1*qIr_2 - qJr_2*qIr_1;
    float qkrxqir_1 = qJr_2*qIr_0 - qJr_0*qIr_2;
    float qkrxqir_2 = qJr_0*qIr_1 - qJr_1*qIr_0;

#endif

#if defined T1 || defined T3

    float qixqk_0 = atomI.labFrameQuadrupole[1]*atomJ.labFrameQuadrupole[2] + atomI.labFrameQuadrupole[4]*atomJ.labFrameQuadrupole[5] + atomI.labFrameQuadrupole[7]*atomJ.labFrameQuadrupole[8] -
                    atomI.labFrameQuadrupole[2]*atomJ.labFrameQuadrupole[1] - atomI.labFrameQuadrupole[5]*atomJ.labFrameQuadrupole[4] - atomI.labFrameQuadrupole[8]*atomJ.labFrameQuadrupole[7];

    float qixqk_1 = atomI.labFrameQuadrupole[2]*atomJ.labFrameQuadrupole[0] + atomI.labFrameQuadrupole[5]*atomJ.labFrameQuadrupole[3] + atomI.labFrameQuadrupole[8]*atomJ.labFrameQuadrupole[6] -
                    atomI.labFrameQuadrupole[0]*atomJ.labFrameQuadrupole[2] - atomI.labFrameQuadrupole[3]*atomJ.labFrameQuadrupole[5] - atomI.labFrameQuadrupole[6]*atomJ.labFrameQuadrupole[8];

    float qixqk_2 = atomI.labFrameQuadrupole[0]*atomJ.labFrameQuadrupole[1] + atomI.labFrameQuadrupole[3]*atomJ.labFrameQuadrupole[4] + atomI.labFrameQuadrupole[6]*atomJ.labFrameQuadrupole[7] -
                    atomI.labFrameQuadrupole[1]*atomJ.labFrameQuadrupole[0] - atomI.labFrameQuadrupole[4]*atomJ.labFrameQuadrupole[3] - atomI.labFrameQuadrupole[7]*atomJ.labFrameQuadrupole[6];

#endif

#ifdef T1
    float ttm2_0  = -rr3*dixdk_0 + gf2*dixr_0-gf5*rxqir_0 + 2.0f*rr5*(dixqkr_0 + dkxqir_0 + rxqidk_0-2.0f*qixqk_0) - 4.0f*rr7*(rxqikr_0 + qkrxqir_0);
    float ttm2_1  = -rr3*dixdk_1 + gf2*dixr_1-gf5*rxqir_1 + 2.0f*rr5*(dixqkr_1 + dkxqir_1 + rxqidk_1-2.0f*qixqk_1) - 4.0f*rr7*(rxqikr_1 + qkrxqir_1);
    float ttm2_2  = -rr3*dixdk_2 + gf2*dixr_2-gf5*rxqir_2 + 2.0f*rr5*(dixqkr_2 + dkxqir_2 + rxqidk_2-2.0f*qixqk_2) - 4.0f*rr7*(rxqikr_2 + qkrxqir_2);

    float ttm2i_0 = -(dixuk_0*psc3+dixukp_0*dsc3)*0.5f + gti2*dixr_0 + ((ukxqir_0+ rxqiuk_0)*psc5 + (ukxqirp_0 + rxqiukp_0)*dsc5) - gti5*rxqir_0;
    float ttm2i_1 = -(dixuk_1*psc3+dixukp_1*dsc3)*0.5f + gti2*dixr_1 + ((ukxqir_1+ rxqiuk_1)*psc5 + (ukxqirp_1 + rxqiukp_1)*dsc5) - gti5*rxqir_1;
    float ttm2i_2 = -(dixuk_2*psc3+dixukp_2*dsc3)*0.5f + gti2*dixr_2 + ((ukxqir_2+ rxqiuk_2)*psc5 + (ukxqirp_2 + rxqiukp_2)*dsc5) - gti5*rxqir_2;
#endif

#ifdef T3
    float qJqIr_0 = atomJ.labFrameQuadrupole[0]*qIr_0 + atomJ.labFrameQuadrupole[3]*qIr_1 + atomJ.labFrameQuadrupole[6]*qIr_2;
    float qJqIr_1 = atomJ.labFrameQuadrupole[1]*qIr_0 + atomJ.labFrameQuadrupole[4]*qIr_1 + atomJ.labFrameQuadrupole[7]*qIr_2;
    float qJqIr_2 = atomJ.labFrameQuadrupole[2]*qIr_0 + atomJ.labFrameQuadrupole[5]*qIr_1 + atomJ.labFrameQuadrupole[8]*qIr_2;

    float qJdI_0 = atomJ.labFrameQuadrupole[0]*atomI.labFrameDipole[0] + atomJ.labFrameQuadrupole[3]*atomI.labFrameDipole[1] + atomJ.labFrameQuadrupole[6]*atomI.labFrameDipole[2];
    float qJdI_1 = atomJ.labFrameQuadrupole[1]*atomI.labFrameDipole[0] + atomJ.labFrameQuadrupole[4]*atomI.labFrameDipole[1] + atomJ.labFrameQuadrupole[7]*atomI.labFrameDipole[2];
    float qJdI_2 = atomJ.labFrameQuadrupole[2]*atomI.labFrameDipole[0] + atomJ.labFrameQuadrupole[5]*atomI.labFrameDipole[1] + atomJ.labFrameQuadrupole[8]*atomI.labFrameDipole[2];

    float dkxr_0 = atomJ.labFrameDipole[1]*zr - atomJ.labFrameDipole[2]*yr;
    float dkxr_1 = atomJ.labFrameDipole[2]*xr - atomJ.labFrameDipole[0]*zr;
    float dkxr_2 = atomJ.labFrameDipole[0]*yr - atomJ.labFrameDipole[1]*xr;

    float rxqkr_0 = yr*qJr_2 - zr*qJr_1;
    float rxqkr_1 = zr*qJr_0 - xr*qJr_2;
    float rxqkr_2 = xr*qJr_1 - yr*qJr_0;

    float dixqkr_0 = atomI.labFrameDipole[1]*qJr_2 - atomI.labFrameDipole[2]*qJr_1;
    float dixqkr_1 = atomI.labFrameDipole[2]*qJr_0 - atomI.labFrameDipole[0]*qJr_2;
    float dixqkr_2 = atomI.labFrameDipole[0]*qJr_1 - atomI.labFrameDipole[1]*qJr_0;

    float dkxqir_0 = atomJ.labFrameDipole[1]*qIr_2 - atomJ.labFrameDipole[2]*qIr_1;
    float dkxqir_1 = atomJ.labFrameDipole[2]*qIr_0 - atomJ.labFrameDipole[0]*qIr_2;
    float dkxqir_2 = atomJ.labFrameDipole[0]*qIr_1 - atomJ.labFrameDipole[1]*qIr_0;

    float rxqkdi_0 = yr*qJdI_2 - zr*qJdI_1;
    float rxqkdi_1 = zr*qJdI_0 - xr*qJdI_2;
    float rxqkdi_2 = xr*qJdI_1 - yr*qJdI_0;

    float rxqkir_0 = yr*qJqIr_2 - zr*qJqIr_1;
    float rxqkir_1 = zr*qJqIr_0 - xr*qJqIr_2;
    float rxqkir_2 = xr*qJqIr_1 - yr*qJqIr_0;

    float qkrxqir_0 = qJr_1*qIr_2 - qJr_2*qIr_1;
    float qkrxqir_1 = qJr_2*qIr_0 - qJr_0*qIr_2;
    float qkrxqir_2 = qJr_0*qIr_1 - qJr_1*qIr_0;

    float dkxui_0 = atomJ.labFrameDipole[1]*atomI.inducedDipole[2] - atomJ.labFrameDipole[2]*atomI.inducedDipole[1];
    float dkxui_1 = atomJ.labFrameDipole[2]*atomI.inducedDipole[0] - atomJ.labFrameDipole[0]*atomI.inducedDipole[2]; 
    float dkxui_2 = atomJ.labFrameDipole[0]*atomI.inducedDipole[1] - atomJ.labFrameDipole[1]*atomI.inducedDipole[0];

    float dkxuip_0 = atomJ.labFrameDipole[1]*atomI.inducedDipoleP[2] - atomJ.labFrameDipole[2]*atomI.inducedDipoleP[1];
    float dkxuip_1 = atomJ.labFrameDipole[2]*atomI.inducedDipoleP[0] - atomJ.labFrameDipole[0]*atomI.inducedDipoleP[2];
    float dkxuip_2 = atomJ.labFrameDipole[0]*atomI.inducedDipoleP[1] - atomJ.labFrameDipole[1]*atomI.inducedDipoleP[0];

    float uixqkrp_0 = atomI.inducedDipoleP[1]*qJr_2 - atomI.inducedDipoleP[2]*qJr_1;
    float uixqkrp_1 = atomI.inducedDipoleP[2]*qJr_0 - atomI.inducedDipoleP[0]*qJr_2;
    float uixqkrp_2 = atomI.inducedDipoleP[0]*qJr_1 - atomI.inducedDipoleP[1]*qJr_0;

    float uixqkr_0 = atomI.inducedDipole[1]*qJr_2 - atomI.inducedDipole[2]*qJr_1;
    float uixqkr_1 = atomI.inducedDipole[2]*qJr_0 - atomI.inducedDipole[0]*qJr_2;
    float uixqkr_2 = atomI.inducedDipole[0]*qJr_1 - atomI.inducedDipole[1]*qJr_0;

    float rxqkuip_0 = yr*qJuIp_2 - zr*qJuIp_1;
    float rxqkuip_1 = zr*qJuIp_0 - xr*qJuIp_2;
    float rxqkuip_2 = xr*qJuIp_1 - yr*qJuIp_0;

    float rxqkui_0 = yr*qJuI_2 - zr*qJuI_1;
    float rxqkui_1 = zr*qJuI_0 - xr*qJuI_2;
    float rxqkui_2 = xr*qJuI_1 - yr*qJuI_0;

    float ttm3_0   =  rr3*dixdk_0 + gf3*dkxr_0 - gf6*rxqkr_0 - 2.0f*rr5*(dixqkr_0 + dkxqir_0 + rxqkdi_0 - 2.0f*qixqk_0) - 4.0f*rr7*(rxqkir_0 - qkrxqir_0);
    float ttm3_1   =  rr3*dixdk_1 + gf3*dkxr_1 - gf6*rxqkr_1 - 2.0f*rr5*(dixqkr_1 + dkxqir_1 + rxqkdi_1 - 2.0f*qixqk_1) - 4.0f*rr7*(rxqkir_1 - qkrxqir_1);
    float ttm3_2   =  rr3*dixdk_2 + gf3*dkxr_2 - gf6*rxqkr_2 - 2.0f*rr5*(dixqkr_2 + dkxqir_2 + rxqkdi_2 - 2.0f*qixqk_2) - 4.0f*rr7*(rxqkir_2 - qkrxqir_2);

    float ttm3i_0  = -(dkxui_0*psc3+ dkxuip_0*dsc3)*0.5f + gti3*dkxr_0 - ((uixqkr_0 + rxqkui_0)*psc5 + (uixqkrp_0 + rxqkuip_0)*dsc5) - gti6*rxqkr_0;
    float ttm3i_1  = -(dkxui_1*psc3+ dkxuip_1*dsc3)*0.5f + gti3*dkxr_1 - ((uixqkr_1 + rxqkui_1)*psc5 + (uixqkrp_1 + rxqkuip_1)*dsc5) - gti6*rxqkr_1;
    float ttm3i_2  = -(dkxui_2*psc3+ dkxuip_2*dsc3)*0.5f + gti3*dkxr_2 - ((uixqkr_2 + rxqkui_2)*psc5 + (uixqkrp_2 + rxqkuip_2)*dsc5) - gti6*rxqkr_2;
#endif

    if( scalingFactors[MScaleIndex] < 1.0f ){
    
#ifdef T1
        ttm2_0 *= scalingFactors[MScaleIndex];
        ttm2_1 *= scalingFactors[MScaleIndex];
        ttm2_2 *= scalingFactors[MScaleIndex];
#endif
        
#ifdef T3
        ttm3_0 *= scalingFactors[MScaleIndex];
        ttm3_1 *= scalingFactors[MScaleIndex];
        ttm3_2 *= scalingFactors[MScaleIndex];
#endif
    
    }

#ifdef F1
    outputForce[0]       = -(ftm2_0+ftm2i_0);
    outputForce[1]       = -(ftm2_1+ftm2i_1);
    outputForce[2]       = -(ftm2_2+ftm2i_2);
#endif
    
#ifdef T1
    outputForce[0]       =  (ttm2_0 + ttm2i_0);
    outputForce[1]       =  (ttm2_1 + ttm2i_1);
    outputForce[2]       =  (ttm2_2 + ttm2i_2);
#endif

#ifdef T3
    outputForce[0]       =  (ttm3_0 + ttm3i_0);
    outputForce[1]       =  (ttm3_1 + ttm3i_1);
    outputForce[2]       =  (ttm3_2 + ttm3i_2);
#endif

    return;

}
