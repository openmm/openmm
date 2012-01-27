
static __device__ void SUB_METHOD_NAME( calculatePmeDirectElectrostaticPairIxnT2, _kernel)( 
                                        PmeDirectElectrostaticParticle& atomI, PmeDirectElectrostaticParticle& atomJ,
                                        const float4 delta, const float4 bn
#ifdef APPLY_SCALE
                                        , const float* scalingFactors
#endif
                                        ){

    float xr                    = delta.x;
    float yr                    = delta.y;
    float zr                    = delta.z;
    float rr1                   = delta.w;

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

    float bn1                   = bn.x;
    float bn2                   = bn.y;
    float bn3                   = bn.z;

    // apply Thole polarization damping to scale factors

    float scale3                = 1.0f;
    float scale5                = 1.0f;
    float scale7                = 1.0f;

    float damp                  = atomI.damp*atomJ.damp;
    if( damp != 0.0f ){
        float pgamma  = atomI.thole < atomJ.thole ? atomI.thole : atomJ.thole;
        float ratio   = 1.0f/(rr1*damp);
            damp      = -pgamma*ratio*ratio*ratio;
        if( damp > -50.0f ){
            float expdamp    = expf(damp);
            scale3           = 1.0f - expdamp;
            scale5           = 1.0f - (1.0f-damp)*expdamp;
            scale7           = 1.0f - (1.0f-damp+0.6f*damp*damp)*expdamp;
        }
    }


    float rr3                   = rr1*rr1*rr1;
#ifdef APPLY_SCALE
    float dsc3                  = rr3*(1.0f - scale3*scalingFactors[DScaleIndex]);
    float dsc5                  = (3.0f*rr3*rr1*rr1)* (1.0f - scale5*scalingFactors[DScaleIndex]);
    float dsc7                  = (15.0f*rr3*rr3*rr1)*(1.0f - scale7*scalingFactors[DScaleIndex]);

    float psc3                  = rr3*(1.0f - scale3*scalingFactors[PScaleIndex]);
    float psc5                  = (3.0f*rr3*rr1*rr1)*(1.0f - scale5*scalingFactors[PScaleIndex]);
    float psc7                  = (15.0f*rr3*rr3*rr1)*(1.0f - scale7*scalingFactors[PScaleIndex]);
#else
    float psc3                  = rr3*(1.0f - scale3);
    float psc5                  = (3.0f*rr3*rr1*rr1)*(1.0f - scale5);
    float psc7                  = (15.0f*rr3*rr3*rr1)*(1.0f - scale7);
#endif

    float prefactor1            = 0.5f*(psc3 - bn1);
#ifdef APPLY_SCALE
    float prefactor2            = 0.5f*(dsc3 - bn1);
#endif

    float dixuk1                = di2*atomJ.inducedDipole[2]  - di3*atomJ.inducedDipole[1];
    float dixukp1               = di2*atomJ.inducedDipoleP[2] - di3*atomJ.inducedDipoleP[1];

#ifdef APPLY_SCALE
    float ttm2i1                = prefactor1*dixuk1 + prefactor2*dixukp1;
#else
    float ttm2i1                = prefactor1*(dixuk1 + dixukp1);
#endif

    float dixuk2                = di3*atomJ.inducedDipole[0]  - di1*atomJ.inducedDipole[2];
    float dixukp2               = di3*atomJ.inducedDipoleP[0] - di1*atomJ.inducedDipoleP[2];

#ifdef APPLY_SCALE
    float ttm2i2                = prefactor1*dixuk2 + prefactor2*dixukp2;
#else
    float ttm2i2                = prefactor1*(dixuk2 + dixukp2);
#endif

    float dixuk3                = di1*atomJ.inducedDipole[1]  - di2*atomJ.inducedDipole[0];
    float dixukp3               = di1*atomJ.inducedDipoleP[1] - di2*atomJ.inducedDipoleP[0];
#ifdef APPLY_SCALE
    float ttm2i3                = prefactor1*dixuk3 + prefactor2*dixukp3;
#else
    float ttm2i3                = prefactor1*(dixuk3 + dixukp3);
#endif

    float sci4                  = atomJ.inducedDipole[0]*xr  + atomJ.inducedDipole[1]*yr  + atomJ.inducedDipole[2]*zr;
    float scip4                 = atomJ.inducedDipoleP[0]*xr + atomJ.inducedDipoleP[1]*yr + atomJ.inducedDipoleP[2]*zr;
    float gti2                  = bn2*(sci4+scip4);
#ifdef APPLY_SCALE
    float gtri2                 = (sci4*psc5+scip4*dsc5);
#else
    float gtri2                 = psc5*(sci4+scip4);
#endif
    prefactor1                  = 0.5f*(gti2 - gtri2);

    ttm2i1                     += prefactor1*( di2*zr - di3*yr );
    ttm2i2                     += prefactor1*( di3*xr - di1*zr );
    ttm2i3                     += prefactor1*( di1*yr - di2*xr );

    float qir1                  = qi1*xr + qi2*yr + qi3*zr;
    float qir2                  = qi2*xr + qi5*yr + qi6*zr;
    float qir3                  = qi3*xr + qi6*yr + qi9*zr;

#ifdef APPLY_SCALE
    prefactor1                  = sci4*psc7 + scip4*dsc7 - bn3*(sci4+scip4);
#else
    prefactor1                  = psc7*(sci4+scip4) - bn3*(sci4+scip4);
#endif
    ttm2i1                     += prefactor1*( yr*qir3 - zr*qir2 );
    ttm2i2                     += prefactor1*( zr*qir1 - xr*qir3 );
    ttm2i3                     += prefactor1*( xr*qir2 - yr*qir1 );

    float qiuk1                 = qi1*atomJ.inducedDipole[0]  + qi2*atomJ.inducedDipole[1]  + qi3*atomJ.inducedDipole[2];
    float qiuk2                 = qi2*atomJ.inducedDipole[0]  + qi5*atomJ.inducedDipole[1]  + qi6*atomJ.inducedDipole[2];
    float qiuk3                 = qi3*atomJ.inducedDipole[0]  + qi6*atomJ.inducedDipole[1]  + qi9*atomJ.inducedDipole[2];

    float qiukp1                = qi1*atomJ.inducedDipoleP[0] + qi2*atomJ.inducedDipoleP[1] + qi3*atomJ.inducedDipoleP[2];
    float qiukp2                = qi2*atomJ.inducedDipoleP[0] + qi5*atomJ.inducedDipoleP[1] + qi6*atomJ.inducedDipoleP[2];
    float qiukp3                = qi3*atomJ.inducedDipoleP[0] + qi6*atomJ.inducedDipoleP[1] + qi9*atomJ.inducedDipoleP[2];

    prefactor1                  = (bn2 - psc5);
#ifdef APPLY_SCALE
    prefactor2                  = (bn2 - dsc5);
#endif
    float ukxqir1               = atomJ.inducedDipole[1]*qir3  - atomJ.inducedDipole[2]*qir2;
    float ukxqirp1              = atomJ.inducedDipoleP[1]*qir3 - atomJ.inducedDipoleP[2]*qir2;
    float rxqiuk1               = yr*qiuk3  - zr*qiuk2;
    float rxqiukp1              = yr*qiukp3 - zr*qiukp2;

#ifdef APPLY_SCALE
    ttm2i1                     += prefactor1*(ukxqir1 + rxqiuk1) + prefactor2*(ukxqirp1 + rxqiukp1);
#else
    ttm2i1                     += prefactor1*( ukxqir1 + rxqiuk1 + ukxqirp1 + rxqiukp1 );
#endif

    float ukxqir2               = atomJ.inducedDipole[2]*qir1  - atomJ.inducedDipole[0]*qir3;
    float ukxqirp2              = atomJ.inducedDipoleP[2]*qir1 - atomJ.inducedDipoleP[0]*qir3;
    float rxqiuk2               = zr*qiuk1  - xr*qiuk3;
    float rxqiukp2              = zr*qiukp1 - xr*qiukp3;
#ifdef APPLY_SCALE
    ttm2i2                     += prefactor1*(ukxqir2 + rxqiuk2) + prefactor2*(ukxqirp2 + rxqiukp2);
#else
    ttm2i2                     += prefactor1*( ukxqir2 + rxqiuk2 + ukxqirp2 + rxqiukp2 );
#endif

    float ukxqir3               = atomJ.inducedDipole[0]*qir2  - atomJ.inducedDipole[1]*qir1;
    float ukxqirp3              = atomJ.inducedDipoleP[0]*qir2 - atomJ.inducedDipoleP[1]*qir1;
    float rxqiuk3               = xr*qiuk2  - yr*qiuk1;
    float rxqiukp3              = xr*qiukp2 - yr*qiukp1;
#ifdef APPLY_SCALE
    ttm2i3                     += prefactor1*(ukxqir3 + rxqiuk3) + prefactor2*(ukxqirp3 + rxqiukp3);
#else
    ttm2i3                     += prefactor1*(ukxqir3 + rxqiuk3 + ukxqirp3 + rxqiukp3 );
#endif

    atomI.torque[0]            += ttm2i1;
    atomI.torque[1]            += ttm2i2;
    atomI.torque[2]            += ttm2i3;

/*
    torque[0]            += ttm2i1;
    torque[1]            += ttm2i2;
    torque[2]            += ttm2i3;
*/

    return;

}
