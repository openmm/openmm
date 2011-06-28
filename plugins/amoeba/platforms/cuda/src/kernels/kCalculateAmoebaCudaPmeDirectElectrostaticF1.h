
static __device__ void SUB_METHOD_NAME( calculatePmeDirectElectrostaticPairIxnF1, _kernel)( PmeDirectElectrostaticParticle& atomI, PmeDirectElectrostaticParticle& atomJ,
                                                                                            float4 delta, float4 bn, float bn5, float forceFactor,
#ifdef APPLY_SCALE
                                                                                            const float* scalingFactors,
#endif
                                                                                            float force[3], float* energy ){

    float xr                 = delta.x;
    float yr                 = delta.y;
    float zr                 = delta.z;
#ifdef APPLY_SCALE
    float rr1                = delta.w;
#endif

    // set the permanent multipole and induced dipole values;

    float ci                 = atomI.q;

    float di1                = atomI.labFrameDipole[0];
    float di2                = atomI.labFrameDipole[1];
    float di3                = atomI.labFrameDipole[2];

    float qi1                = atomI.labFrameQuadrupole[0];
    float qi2                = atomI.labFrameQuadrupole[1];
    float qi3                = atomI.labFrameQuadrupole[2];
    float qi5                = atomI.labFrameQuadrupole[3];
    float qi6                = atomI.labFrameQuadrupole[4];
    //float qi9                = atomI.labFrameQuadrupole[5];
    float qi9                   = -(atomI.labFrameQuadrupole[0] + atomI.labFrameQuadrupole[3]);

    float ck                 = atomJ.q;
    float dk1                = atomJ.labFrameDipole[0];
    float dk2                = atomJ.labFrameDipole[1];
    float dk3                = atomJ.labFrameDipole[2];

    float qk1                = atomJ.labFrameQuadrupole[0];
    float qk2                = atomJ.labFrameQuadrupole[1];
    float qk3                = atomJ.labFrameQuadrupole[2];
    float qk5                = atomJ.labFrameQuadrupole[3];
    float qk6                = atomJ.labFrameQuadrupole[4];
//    float qk9                = atomJ.labFrameQuadrupole[5];
    float qk9                = -(atomJ.labFrameQuadrupole[0] + atomJ.labFrameQuadrupole[3]);

    float bn1                = bn.x;
    float bn2                = bn.y;
    float bn3                = bn.z;
    float bn4                = bn.w;

#ifdef APPLY_SCALE
    float offset             = 1.0f-scalingFactors[MScaleIndex];
    float rr3                = rr1*rr1*rr1;
    float gf4                = 2.0f*(bn2 - 3.0f*offset*rr3*rr1*rr1);
#else
    float gf4                = 2.0f*bn2;
#endif
    float qidk1              = qi1*dk1 + qi2*dk2 + qi3*dk3;
    float qkdi1              = qk1*di1 + qk2*di2 + qk3*di3;
    float ftm21              = gf4*(qkdi1-qidk1);

    float qidk2              = qi2*dk1 + qi5*dk2 + qi6*dk3;
    float qkdi2              = qk2*di1 + qk5*di2 + qk6*di3;
    float ftm22              = gf4*(qkdi2-qidk2);

    float qidk3              = qi3*dk1 + qi6*dk2 + qi9*dk3;
    float qkdi3              = qk3*di1 + qk6*di2 + qk9*di3;
    float ftm23              = gf4*(qkdi3-qidk3);

    float qir1               = qi1*xr + qi2*yr + qi3*zr;
    float qir2               = qi2*xr + qi5*yr + qi6*zr;
    float qir3               = qi3*xr + qi6*yr + qi9*zr;

    float qkr1               = qk1*xr + qk2*yr + qk3*zr;
    float qkr2               = qk2*xr + qk5*yr + qk6*zr;
    float qkr3               = qk3*xr + qk6*yr + qk9*zr;

#ifdef APPLY_SCALE
    float gf7                = 4.0f*(bn3 - 15.0f*offset*rr3*rr3*rr1);
#else
    float gf7                = 4.0f*bn3;
#endif
    float qiqkr1             = qi1*qkr1 + qi2*qkr2 + qi3*qkr3;
    float qkqir1             = qk1*qir1 + qk2*qir2 + qk3*qir3;
    ftm21                   += gf7*(qiqkr1+qkqir1);

    float qiqkr2             = qi2*qkr1 + qi5*qkr2 + qi6*qkr3;
    float qkqir2             = qk2*qir1 + qk5*qir2 + qk6*qir3;
    ftm22                   += gf7*(qiqkr2+qkqir2);

    float qiqkr3             = qi3*qkr1 + qi6*qkr2 + qi9*qkr3;
    float qkqir3             = qk3*qir1 + qk6*qir2 + qk9*qir3;
    ftm23                   += gf7*(qiqkr3+qkqir3);

    // calculate the scalar products for permanent components

    float gl6                = di1*dk1   + di2*dk2   + di3*dk3;
    float gl7                =  2.0f*( qir1*dk1  + qir2*dk2  + qir3*dk3 - ( qkr1*di1  + qkr2*di2  + qkr3*di3 ) );
    float gl5                = -4.0f*(qir1*qkr1 + qir2*qkr2 + qir3*qkr3);

    float gl8                =  2.0f*(qi1*qk1 + qi2*qk2 + qi3*qk3 + qi2*qk2 + qi5*qk5 + qi6*qk6 + qi3*qk3 + qi6*qk6 + qi9*qk9 );

    float sc3                = di1*xr  + di2*yr  + di3*zr;
    float sc5                = qir1*xr + qir2*yr + qir3*zr;
    float sc4                = dk1*xr  + dk2*yr  + dk3*zr;
    float sc6                = qkr1*xr + qkr2*yr + qkr3*zr;

    float gl0                = ci*ck;
    float gl1                = ck*sc3 - ci*sc4;
    float gl2                = ci*sc6 + ck*sc5 - sc3*sc4;
    float gl3                = sc3*sc6 - sc4*sc5;
    float gl4                = sc5*sc6;

#ifdef APPLY_SCALE
    //forceTorqueEnergy->w    += forceFactor*(-offset*rr1*gl0 + (bn1-offset*rr3)*(gl1+gl6) + (bn2-offset*(3.0f*rr3*rr1*rr1))*(gl2+gl7+gl8) + (bn3-offset*(15.0f*rr3*rr3*rr1))*(gl3+gl5) + (bn4-offset*(105.0f*rr3*rr3*rr3))*gl4);
    *energy                 += forceFactor*(-offset*rr1*gl0 + (bn1-offset*rr3)*(gl1+gl6) + (bn2-offset*(3.0f*rr3*rr1*rr1))*(gl2+gl7+gl8) + (bn3-offset*(15.0f*rr3*rr3*rr1))*(gl3+gl5) + (bn4-offset*(105.0f*rr3*rr3*rr3))*gl4);
#else
    //forceTorqueEnergy->w    += bn1*(gl1+gl6) + bn2*(gl2+gl7+gl8) + bn3*(gl3+gl5) + bn4*gl4;
    *energy                 += forceFactor*(bn1*(gl1+gl6) + bn2*(gl2+gl7+gl8) + bn3*(gl3+gl5) + bn4*gl4);
    
#endif

    float gf1                = bn1*gl0 + bn2*(gl1+gl6) + bn3*(gl2+gl7+gl8) + bn4*(gl3+gl5) + bn5*gl4;
#ifdef APPLY_SCALE
          gf1               -= offset*(rr3*gl0 + (3.0f*rr3*rr1*rr1)*(gl1+gl6) + (15.0f*rr3*rr3*rr1)*(gl2+gl7+gl8) + (105.0f*rr3*rr3*rr3)*(gl3+gl5) + (945.0f*rr3*rr3*rr3*rr1*rr1)*gl4);
#endif
    ftm21                   += gf1*xr;
    ftm22                   += gf1*yr;
    ftm23                   += gf1*zr;

#ifdef APPLY_SCALE
    float gf2                = -ck*bn1 + sc4*bn2 - sc6*bn3 - offset*(-ck*rr3 + sc4*(3.0f*rr3*rr1*rr1) - sc6*(15.0f*rr3*rr3*rr1));
#else
    float gf2                = -ck*bn1 + sc4*bn2 - sc6*bn3;
#endif
    ftm21                   += gf2*di1;
    ftm22                   += gf2*di2;
    ftm23                   += gf2*di3;

#ifdef APPLY_SCALE
    float gf5                = 2.0f*(-ck*bn2+sc4*bn3-sc6*bn4 - offset*(-ck*(3.0f*rr3*rr1*rr1)+sc4*(15.0f*rr3*rr3*rr1)-sc6*(105.0f*rr3*rr3*rr3)));
#else
    float gf5                = 2.0f*(-ck*bn2+sc4*bn3-sc6*bn4);
#endif
    ftm21                   += gf5*qir1;
    ftm22                   += gf5*qir2;
    ftm23                   += gf5*qir3;

#ifdef APPLY_SCALE
    float gf3                = ci*bn1 + sc3*bn2 + sc5*bn3 - offset*(ci*rr3 + sc3*(3.0f*rr3*rr1*rr1) + sc5*(15.0f*rr3*rr3*rr1));
#else
    float gf3                = ci*bn1 + sc3*bn2 + sc5*bn3;
#endif
    ftm21                   += gf3*dk1;
    ftm22                   += gf3*dk2;
    ftm23                   += gf3*dk3;

#ifdef APPLY_SCALE
    float gf6                = 2.0f*(-ci*bn2-sc3*bn3-sc5*bn4 - offset*(-ci*(3.0f*rr3*rr1*rr1)-sc3*(15.0f*rr3*rr3*rr1)-sc5*(105.0f*rr3*rr3*rr3)));
#else
    float gf6                = 2.0f*(-ci*bn2-sc3*bn3-sc5*bn4);
#endif

    ftm21                   += gf6*qkr1;
    ftm22                   += gf6*qkr2;
    ftm23                   += gf6*qkr3;

    force[0]                 = ftm21;
    force[1]                 = ftm22;
    force[2]                 = ftm23;
/*
    if( forceFactor == 1.0f ){
        atomJ.force[0]      -= ftm21;
        atomJ.force[1]      -= ftm22;
        atomJ.force[2]      -= ftm23;
    }
    atomI.force[0]      += ftm21;
    atomI.force[1]      += ftm22;
    atomI.force[2]      += ftm23;
*/
    return;

}
