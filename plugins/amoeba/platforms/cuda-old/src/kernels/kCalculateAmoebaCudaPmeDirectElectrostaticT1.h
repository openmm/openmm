
static __device__ void SUB_METHOD_NAME( calculatePmeDirectElectrostaticPairIxnT1, _kernel )( 
                                        PmeDirectElectrostaticParticle& atomI, PmeDirectElectrostaticParticle& atomJ,
                                        const float4 delta, const float4 bn
#ifdef APPLY_SCALE
                                        , const float* scalingFactors
#endif
                                        ){

    float xr                    = delta.x;
    float yr                    = delta.y;
    float zr                    = delta.z;
#ifdef APPLY_SCALE
    float rr1                   = delta.w;
#endif

    // set the permanent multipole and induced dipole values;

    float di1                   = atomI.labFrameDipole[0];
    float di2                   = atomI.labFrameDipole[1];
    float di3                   = atomI.labFrameDipole[2];

    float qi1                   = atomI.labFrameQuadrupole[0];
    float qi2                   = atomI.labFrameQuadrupole[1];
    float qi3                   = atomI.labFrameQuadrupole[2];
    float qi5                   = atomI.labFrameQuadrupole[3];
    float qi6                   = atomI.labFrameQuadrupole[4];
    //float qi9                   = atomI.labFrameQuadrupole[5];
    float qi9                   = -(atomI.labFrameQuadrupole[0] + atomI.labFrameQuadrupole[3]);

    float ck                    = atomJ.q;

    float dk1                   = atomJ.labFrameDipole[0];
    float dk2                   = atomJ.labFrameDipole[1];
    float dk3                   = atomJ.labFrameDipole[2];

    float qk1                   = atomJ.labFrameQuadrupole[0];
    float qk2                   = atomJ.labFrameQuadrupole[1];
    float qk3                   = atomJ.labFrameQuadrupole[2];
    float qk5                   = atomJ.labFrameQuadrupole[3];
    float qk6                   = atomJ.labFrameQuadrupole[4];
    //float qk9                   = atomJ.labFrameQuadrupole[5];
    float qk9                   = -(atomJ.labFrameQuadrupole[0] + atomJ.labFrameQuadrupole[3]);

    float bn1                   = bn.x;
    float bn2                   = bn.y;
    float bn3                   = bn.z;
    float bn4                   = bn.w;

    // apply Thole polarization damping to scale factors

#ifdef APPLY_SCALE
    float rr2                   = rr1*rr1;
    float rr3                   = rr1*rr2;
    float rr5                   = 3.0f*rr3*rr2;
    float rr7                   = 5.0f*rr5*rr2;
    float rr9                   = 7.0f*rr7*rr2;

    float scale                 = 1.0f-scalingFactors[MScaleIndex];
    float prefactor             = scale*rr3 - bn1;
#else
    float prefactor             = -bn1;
#endif
    float dixdk1                = di2*dk3 - di3*dk2;
    float ttm21                 = prefactor*dixdk1;

    float dixdk2                = di3*dk1 - di1*dk3;
    float ttm22                 = prefactor*dixdk2;

    float dixdk3                = di1*dk2 - di2*dk1;
    float ttm23                 = prefactor*dixdk3;

    float qir1                  = qi1*xr + qi2*yr + qi3*zr;
    float qir2                  = qi2*xr + qi5*yr + qi6*zr;
    float qir3                  = qi3*xr + qi6*yr + qi9*zr;

    float qkr1                  = qk1*xr + qk2*yr + qk3*zr;
    float qkr2                  = qk2*xr + qk5*yr + qk6*zr;
    float qkr3                  = qk3*xr + qk6*yr + qk9*zr;

    float qiqkr1                = qi1*qkr1 + qi2*qkr2 + qi3*qkr3;
    float qiqkr2                = qi2*qkr1 + qi5*qkr2 + qi6*qkr3;
    float qiqkr3                = qi3*qkr1 + qi6*qkr2 + qi9*qkr3;

    float rxqikr1               = yr*qiqkr3 - zr*qiqkr2;
    float qkrxqir1              = qkr2*qir3 - qkr3*qir2;
#ifdef APPLY_SCALE
    prefactor                   = 4.0f*(bn3 - scale*rr7);
#else
    prefactor                   = 4.0f*bn3;
#endif
    ttm21                      -= prefactor*(rxqikr1+qkrxqir1);

    float rxqikr2               = zr*qiqkr1 - xr*qiqkr3;
    float qkrxqir2              = qkr3*qir1 - qkr1*qir3;
    ttm22                      -= prefactor*(rxqikr2+qkrxqir2);

    float rxqikr3               = xr*qiqkr2 - yr*qiqkr1;
    float qkrxqir3              = qkr1*qir2 - qkr2*qir1;
    ttm23                      -= prefactor*(rxqikr3+qkrxqir3);

    float qidk1                 = qi1*dk1 + qi2*dk2 + qi3*dk3;
    float qidk2                 = qi2*dk1 + qi5*dk2 + qi6*dk3;
    float qidk3                 = qi3*dk1 + qi6*dk2 + qi9*dk3;

    float dixqkr1               = di2*qkr3 - di3*qkr2;
    float dkxqir1               = dk2*qir3 - dk3*qir2;
    float rxqidk1               = yr*qidk3 - zr*qidk2;
    float qixqk1                = qi2*qk3 + qi5*qk6 + qi6*qk9 - qi3*qk2 - qi6*qk5 - qi9*qk6;
#ifdef APPLY_SCALE
    prefactor                   = 2.0f*(bn2 - scale*rr5);
#else
    prefactor                   = 2.0f*bn2;
#endif
    ttm21                      += prefactor*(dixqkr1+dkxqir1+rxqidk1-2.0f*qixqk1);
 
    float dixqkr2               = di3*qkr1 - di1*qkr3;
    float dkxqir2               = dk3*qir1 - dk1*qir3;
    float rxqidk2               = zr*qidk1 - xr*qidk3;
    float qixqk2                = qi3*qk1 + qi6*qk2 + qi9*qk3 - qi1*qk3 - qi2*qk6 - qi3*qk9;
    ttm22                      += prefactor*(dixqkr2+dkxqir2+rxqidk2-2.0f*qixqk2);

    float dixqkr3               = di1*qkr2 - di2*qkr1;
    float dkxqir3               = dk1*qir2 - dk2*qir1;
    float rxqidk3               = xr*qidk2 - yr*qidk1;
    float qixqk3                = qi1*qk2 + qi2*qk5 + qi3*qk6 - qi2*qk1 - qi5*qk2 - qi6*qk3;
    ttm23                      += prefactor*(dixqkr3+dkxqir3+rxqidk3-2.0f*qixqk3);

    float sc4                   = dk1*xr + dk2*yr + dk3*zr;
    float sc6                   = qkr1*xr + qkr2*yr + qkr3*zr;

    float gf2                   = -ck*bn1 + sc4*bn2 - sc6*bn3;
#ifdef APPLY_SCALE
    float gfr2                  = -ck*rr3 + sc4*rr5 - sc6*rr7;
    prefactor                   = (gf2 - scale*gfr2);
#else
    prefactor                   = gf2;
#endif
    ttm21                      += prefactor*(di2*zr - di3*yr);
    ttm22                      += prefactor*(di3*xr - di1*zr);
    ttm23                      += prefactor*(di1*yr - di2*xr);

    float gf5                   = (-ck*bn2+sc4*bn3-sc6*bn4);
#ifdef APPLY_SCALE
    float gfr5                  = (-ck*rr5+sc4*rr7-sc6*rr9); 
    prefactor                   = 2.0f*(gf5 - scale*gfr5);
#else
    prefactor                   = 2.0f*gf5;
#endif

    float rxqir1                = yr*qir3 - zr*qir2;
    float rxqir2                = zr*qir1 - xr*qir3;
    float rxqir3                = xr*qir2 - yr*qir1;
    ttm21                      -= prefactor*rxqir1; 
    ttm22                      -= prefactor*rxqir2;
    ttm23                      -= prefactor*rxqir3;

    atomI.torque[0]            += ttm21;
    atomI.torque[1]            += ttm22;
    atomI.torque[2]            += ttm23;
/*

    torque[0]             = ttm21;
    torque[1]             = ttm22;
    torque[2]             = ttm23;
*/

    return;

}
