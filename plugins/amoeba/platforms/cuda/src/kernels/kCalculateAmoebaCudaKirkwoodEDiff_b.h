
/*****************************************************************************

       ediff1 correct vacuum to SCRF derivatives

       calculates the energy and derivatives of polarizing
       the vacuum induced dipoles to their SCRF polarized values

******************************************************************************/

__device__ void SUB_METHOD_NAME( calculateKirkwoodEDiffPairIxn, _kernel)( KirkwoodEDiffParticle& atomI,  KirkwoodEDiffParticle& atomJ,
#ifdef APPLY_SCALE
                                                      float pscale, float dscale,
#endif
#ifdef F1
                                                      float*  outputEnergy,
#endif
                                                      float*  outputForce ){

    const float uscale = 1.0f;

    // deltaR

    float xr          = atomJ.x - atomI.x;
    float yr          = atomJ.y - atomI.y;
    float zr          = atomJ.z - atomI.z;

    float r22         = xr*xr + yr*yr + zr*zr;

    float r           = sqrtf(r22);
    float rr1         = 1.0f/r;
    float rr2         = rr1*rr1;
    float rr3         = rr1*rr2;

    float scale3      = 1.0f;
    float scale5      = 1.0f;
    float scale7      = 1.0f;

#ifdef F1
    float ddsc3_1     = 0.0f;
    float ddsc3_2     = 0.0f;
    float ddsc3_3     = 0.0f;

    float ddsc5_1     = 0.0f;
    float ddsc5_2     = 0.0f;
    float ddsc5_3     = 0.0f;

    float ddsc7_1     = 0.0f;
    float ddsc7_2     = 0.0f;
    float ddsc7_3     = 0.0f;

    float ftm2i1      = 0.0f;
    float ftm2i2      = 0.0f;
    float ftm2i3      = 0.0f;
#endif

    // apply Thole polarization damping to scale factors
 
    float damp = atomI.damp*atomJ.damp;
    if( damp != 0.0f ){
        float pgamma    = atomJ.thole > atomI.thole ? atomI.thole : atomJ.thole;
        float ratio     = (r/damp);
        damp            = -pgamma*ratio*ratio*ratio;
        if( damp > -50.0f){
            float dampE  = expf( damp );
            float damp2  = damp*damp;
            scale3       = 1.0f - dampE;
            scale5       = 1.0f - (1.0f - damp)*dampE;
            scale7       = 1.0f - (1.0f - damp + 0.6f*damp2)*dampE;

#ifdef F1
            ddsc3_1     = -3.0f*damp*exp(damp)*xr*rr2*rr3;
            ddsc3_2     = -3.0f*damp*exp(damp)*yr*rr2*rr3;
            ddsc3_3     = -3.0f*damp*exp(damp)*zr*rr2*rr3;

            ddsc5_1     = -3.0f*damp*ddsc3_1*rr2;
            ddsc5_2     = -3.0f*damp*ddsc3_2*rr2;
            ddsc5_3     = -3.0f*damp*ddsc3_3*rr2;

            ddsc7_1     = -5.0f*(0.2f+0.6f*damp)*ddsc5_1*rr2;
            ddsc7_2     = -5.0f*(0.2f+0.6f*damp)*ddsc5_2*rr2;
            ddsc7_3     = -5.0f*(0.2f+0.6f*damp)*ddsc5_3*rr2;
#endif
        }
    }

    float scale3i       = 3.0f*scale3*uscale*rr3*rr2;
    float scale5i       = 3.0f*scale5*uscale*rr3*rr2;

#ifdef APPLY_SCALE
    float dsc3          = scale3*dscale*rr3;
    float dsc5          = 3.0f*scale5*dscale*rr3*rr2;
    float dsc7          = 15.0f*scale7*dscale*rr3*rr3*rr1;

    float psc3          = scale3*pscale*rr3;
    float psc5          = 3.0f*scale5*pscale*rr3*rr2;
    float psc7          = 15.0f*scale7*pscale*rr3*rr3*rr1;
#else
    float psc3          = scale3*rr3;
    float psc5          = 3.0f*scale5*rr3*rr2;
    float psc7          = 15.0f*scale7*rr3*rr3*rr1;
#endif
 
#ifdef T1
    float dixr1             = atomI.labFrameDipole[1]*zr - atomI.labFrameDipole[2]*yr;
    float dixr2             = atomI.labFrameDipole[2]*xr - atomI.labFrameDipole[0]*zr;
    float dixr3             = atomI.labFrameDipole[0]*yr - atomI.labFrameDipole[1]*xr;
#endif

#ifdef T3
    float dkxr1             = atomJ.labFrameDipole[1]*zr - atomJ.labFrameDipole[2]*yr;
    float dkxr2             = atomJ.labFrameDipole[2]*xr - atomJ.labFrameDipole[0]*zr;
    float dkxr3             = atomJ.labFrameDipole[0]*yr - atomJ.labFrameDipole[1]*xr;
#endif

    float qir1              = atomI.labFrameQuadrupole_XX*xr + atomI.labFrameQuadrupole_XY*yr + atomI.labFrameQuadrupole_XZ*zr;
    float qir2              = atomI.labFrameQuadrupole_XY*xr + atomI.labFrameQuadrupole_YY*yr + atomI.labFrameQuadrupole_YZ*zr;
    float qir3              = atomI.labFrameQuadrupole_XZ*xr + atomI.labFrameQuadrupole_YZ*yr + atomI.labFrameQuadrupole_ZZ*zr;

    float qkr1              = atomJ.labFrameQuadrupole_XX*xr + atomJ.labFrameQuadrupole_XY*yr + atomJ.labFrameQuadrupole_XZ*zr;
    float qkr2              = atomJ.labFrameQuadrupole_XY*xr + atomJ.labFrameQuadrupole_YY*yr + atomJ.labFrameQuadrupole_YZ*zr;
    float qkr3              = atomJ.labFrameQuadrupole_XZ*xr + atomJ.labFrameQuadrupole_YZ*yr + atomJ.labFrameQuadrupole_ZZ*zr;

#ifdef T1
    float rxqir1            = yr*qir3 - zr*qir2;
    float rxqir2            = zr*qir1 - xr*qir3;
    float rxqir3            = xr*qir2 - yr*qir1;
#endif

#ifdef T3
    float rxqkr1            = yr*qkr3 - zr*qkr2;
    float rxqkr2            = zr*qkr1 - xr*qkr3;
    float rxqkr3            = xr*qkr2 - yr*qkr1;
#endif

    // get intermediate variables for permanent energy terms
 
    float sc3               = atomI.labFrameDipole[0]*xr  + atomI.labFrameDipole[1]*yr  + atomI.labFrameDipole[2]*zr;
    float sc4               = atomJ.labFrameDipole[0]*xr  + atomJ.labFrameDipole[1]*yr  + atomJ.labFrameDipole[2]*zr;
    float sc5               = qir1*xr + qir2*yr + qir3*zr;
    float sc6               = qkr1*xr + qkr2*yr + qkr3*zr;
 
#ifdef T1
    float dixuk1            = atomI.labFrameDipole[1]*atomJ.inducedDipoleS[2] - atomI.labFrameDipole[2]*atomJ.inducedDipoleS[1];
    float dixuk2            = atomI.labFrameDipole[2]*atomJ.inducedDipoleS[0] - atomI.labFrameDipole[0]*atomJ.inducedDipoleS[2];
    float dixuk3            = atomI.labFrameDipole[0]*atomJ.inducedDipoleS[1] - atomI.labFrameDipole[1]*atomJ.inducedDipoleS[0];

    float dixukp1           = atomI.labFrameDipole[1]*atomJ.inducedDipolePS[2] - atomI.labFrameDipole[2]*atomJ.inducedDipolePS[1];
    float dixukp2           = atomI.labFrameDipole[2]*atomJ.inducedDipolePS[0] - atomI.labFrameDipole[0]*atomJ.inducedDipolePS[2];
    float dixukp3           = atomI.labFrameDipole[0]*atomJ.inducedDipolePS[1] - atomI.labFrameDipole[1]*atomJ.inducedDipolePS[0];
#endif

#ifdef T3
    float dkxui1            = atomJ.labFrameDipole[1]*atomI.inducedDipoleS[2] - atomJ.labFrameDipole[2]*atomI.inducedDipoleS[1];
    float dkxui2            = atomJ.labFrameDipole[2]*atomI.inducedDipoleS[0] - atomJ.labFrameDipole[0]*atomI.inducedDipoleS[2];
    float dkxui3            = atomJ.labFrameDipole[0]*atomI.inducedDipoleS[1] - atomJ.labFrameDipole[1]*atomI.inducedDipoleS[0];

    float dkxuip1           = atomJ.labFrameDipole[1]*atomI.inducedDipolePS[2] - atomJ.labFrameDipole[2]*atomI.inducedDipolePS[1];
    float dkxuip2           = atomJ.labFrameDipole[2]*atomI.inducedDipolePS[0] - atomJ.labFrameDipole[0]*atomI.inducedDipolePS[2];
    float dkxuip3           = atomJ.labFrameDipole[0]*atomI.inducedDipolePS[1] - atomJ.labFrameDipole[1]*atomI.inducedDipolePS[0];
#endif

#if defined F1 || defined T1
    float qiuk1             = atomI.labFrameQuadrupole_XX*atomJ.inducedDipoleS[0] + atomI.labFrameQuadrupole_XY*atomJ.inducedDipoleS[1] + atomI.labFrameQuadrupole_XZ*atomJ.inducedDipoleS[2];
    float qiuk2             = atomI.labFrameQuadrupole_XY*atomJ.inducedDipoleS[0] + atomI.labFrameQuadrupole_YY*atomJ.inducedDipoleS[1] + atomI.labFrameQuadrupole_YZ*atomJ.inducedDipoleS[2];
    float qiuk3             = atomI.labFrameQuadrupole_XZ*atomJ.inducedDipoleS[0] + atomI.labFrameQuadrupole_YZ*atomJ.inducedDipoleS[1] + atomI.labFrameQuadrupole_ZZ*atomJ.inducedDipoleS[2];

    float qiukp1            = atomI.labFrameQuadrupole_XX*atomJ.inducedDipolePS[0] + atomI.labFrameQuadrupole_XY*atomJ.inducedDipolePS[1] + atomI.labFrameQuadrupole_XZ*atomJ.inducedDipolePS[2];
    float qiukp2            = atomI.labFrameQuadrupole_XY*atomJ.inducedDipolePS[0] + atomI.labFrameQuadrupole_YY*atomJ.inducedDipolePS[1] + atomI.labFrameQuadrupole_YZ*atomJ.inducedDipolePS[2];
    float qiukp3            = atomI.labFrameQuadrupole_XZ*atomJ.inducedDipolePS[0] + atomI.labFrameQuadrupole_YZ*atomJ.inducedDipolePS[1] + atomI.labFrameQuadrupole_ZZ*atomJ.inducedDipolePS[2];
#if defined F1
    qiuk1                  -= atomI.labFrameQuadrupole_XX*atomJ.inducedDipole[0] + atomI.labFrameQuadrupole_XY*atomJ.inducedDipole[1] + atomI.labFrameQuadrupole_XZ*atomJ.inducedDipole[2];
    qiuk2                  -= atomI.labFrameQuadrupole_XY*atomJ.inducedDipole[0] + atomI.labFrameQuadrupole_YY*atomJ.inducedDipole[1] + atomI.labFrameQuadrupole_YZ*atomJ.inducedDipole[2];
    qiuk3                  -= atomI.labFrameQuadrupole_XZ*atomJ.inducedDipole[0] + atomI.labFrameQuadrupole_YZ*atomJ.inducedDipole[1] + atomI.labFrameQuadrupole_ZZ*atomJ.inducedDipole[2];

    qiukp1                 -= atomI.labFrameQuadrupole_XX*atomJ.inducedDipoleP[0] + atomI.labFrameQuadrupole_XY*atomJ.inducedDipoleP[1] + atomI.labFrameQuadrupole_XZ*atomJ.inducedDipoleP[2];
    qiukp2                 -= atomI.labFrameQuadrupole_XY*atomJ.inducedDipoleP[0] + atomI.labFrameQuadrupole_YY*atomJ.inducedDipoleP[1] + atomI.labFrameQuadrupole_YZ*atomJ.inducedDipoleP[2];
    qiukp3                 -= atomI.labFrameQuadrupole_XZ*atomJ.inducedDipoleP[0] + atomI.labFrameQuadrupole_YZ*atomJ.inducedDipoleP[1] + atomI.labFrameQuadrupole_ZZ*atomJ.inducedDipoleP[2];

#ifdef APPLY_SCALE
    ftm2i1                 -= psc5*qiuk1 + dsc5*qiukp1;
    ftm2i2                 -= psc5*qiuk2 + dsc5*qiukp2;
    ftm2i3                 -= psc5*qiuk3 + dsc5*qiukp3;
#else
    ftm2i1                 -= psc5*(qiuk1 + qiukp1);
    ftm2i2                 -= psc5*(qiuk2 + qiukp2);
    ftm2i3                 -= psc5*(qiuk3 + qiukp3);
#endif
#endif
#endif

#if defined F1 || defined T3
    float qkui1             = atomJ.labFrameQuadrupole_XX*atomI.inducedDipoleS[0] + atomJ.labFrameQuadrupole_XY*atomI.inducedDipoleS[1] + atomJ.labFrameQuadrupole_XZ*atomI.inducedDipoleS[2];
    float qkui2             = atomJ.labFrameQuadrupole_XY*atomI.inducedDipoleS[0] + atomJ.labFrameQuadrupole_YY*atomI.inducedDipoleS[1] + atomJ.labFrameQuadrupole_YZ*atomI.inducedDipoleS[2];
    float qkui3             = atomJ.labFrameQuadrupole_XZ*atomI.inducedDipoleS[0] + atomJ.labFrameQuadrupole_YZ*atomI.inducedDipoleS[1] + atomJ.labFrameQuadrupole_ZZ*atomI.inducedDipoleS[2];

    float qkuip1            = atomJ.labFrameQuadrupole_XX*atomI.inducedDipolePS[0] + atomJ.labFrameQuadrupole_XY*atomI.inducedDipolePS[1] + atomJ.labFrameQuadrupole_XZ*atomI.inducedDipolePS[2];
    float qkuip2            = atomJ.labFrameQuadrupole_XY*atomI.inducedDipolePS[0] + atomJ.labFrameQuadrupole_YY*atomI.inducedDipolePS[1] + atomJ.labFrameQuadrupole_YZ*atomI.inducedDipolePS[2];
    float qkuip3            = atomJ.labFrameQuadrupole_XZ*atomI.inducedDipolePS[0] + atomJ.labFrameQuadrupole_YZ*atomI.inducedDipolePS[1] + atomJ.labFrameQuadrupole_ZZ*atomI.inducedDipolePS[2];

#if defined F1
    qkui1                  -= atomJ.labFrameQuadrupole_XX*atomI.inducedDipole[0] + atomJ.labFrameQuadrupole_XY*atomI.inducedDipole[1] + atomJ.labFrameQuadrupole_XZ*atomI.inducedDipole[2];
    qkui2                  -= atomJ.labFrameQuadrupole_XY*atomI.inducedDipole[0] + atomJ.labFrameQuadrupole_YY*atomI.inducedDipole[1] + atomJ.labFrameQuadrupole_YZ*atomI.inducedDipole[2];
    qkui3                  -= atomJ.labFrameQuadrupole_XZ*atomI.inducedDipole[0] + atomJ.labFrameQuadrupole_YZ*atomI.inducedDipole[1] + atomJ.labFrameQuadrupole_ZZ*atomI.inducedDipole[2];

    qkuip1                 -= atomJ.labFrameQuadrupole_XX*atomI.inducedDipoleP[0] + atomJ.labFrameQuadrupole_XY*atomI.inducedDipoleP[1] + atomJ.labFrameQuadrupole_XZ*atomI.inducedDipoleP[2];
    qkuip2                 -= atomJ.labFrameQuadrupole_XY*atomI.inducedDipoleP[0] + atomJ.labFrameQuadrupole_YY*atomI.inducedDipoleP[1] + atomJ.labFrameQuadrupole_YZ*atomI.inducedDipoleP[2];
    qkuip3                 -= atomJ.labFrameQuadrupole_XZ*atomI.inducedDipoleP[0] + atomJ.labFrameQuadrupole_YZ*atomI.inducedDipoleP[1] + atomJ.labFrameQuadrupole_ZZ*atomI.inducedDipoleP[2];

#ifdef APPLY_SCALE
    ftm2i1                 += psc5*qkui1 + dsc5*qkuip1;
    ftm2i2                 += psc5*qkui2 + dsc5*qkuip2;
    ftm2i3                 += psc5*qkui3 + dsc5*qkuip3;
#else
    ftm2i1                 += psc5*(qkui1 + qkuip1);
    ftm2i2                 += psc5*(qkui2 + qkuip2);
    ftm2i3                 += psc5*(qkui3 + qkuip3);
#endif
#endif

#endif

#ifdef T3
    float uixqkr1           = atomI.inducedDipoleS[1]*qkr3 - atomI.inducedDipoleS[2]*qkr2;
    float uixqkr2           = atomI.inducedDipoleS[2]*qkr1 - atomI.inducedDipoleS[0]*qkr3;
    float uixqkr3           = atomI.inducedDipoleS[0]*qkr2 - atomI.inducedDipoleS[1]*qkr1;

    float uixqkrp1          = atomI.inducedDipolePS[1]*qkr3 - atomI.inducedDipolePS[2]*qkr2;
    float uixqkrp2          = atomI.inducedDipolePS[2]*qkr1 - atomI.inducedDipolePS[0]*qkr3;
    float uixqkrp3          = atomI.inducedDipolePS[0]*qkr2 - atomI.inducedDipolePS[1]*qkr1;

    float rxqkuip1          = yr*qkuip3 - zr*qkuip2;
    float rxqkuip2          = zr*qkuip1 - xr*qkuip3;
    float rxqkuip3          = xr*qkuip2 - yr*qkuip1;

    float rxqkui1           = yr*qkui3 - zr*qkui2;
    float rxqkui2           = zr*qkui1 - xr*qkui3;
    float rxqkui3           = xr*qkui2 - yr*qkui1;
#endif

#ifdef T1
    float ukxqir1           = atomJ.inducedDipoleS[1]*qir3 - atomJ.inducedDipoleS[2]*qir2;
    float ukxqir2           = atomJ.inducedDipoleS[2]*qir1 - atomJ.inducedDipoleS[0]*qir3;
    float ukxqir3           = atomJ.inducedDipoleS[0]*qir2 - atomJ.inducedDipoleS[1]*qir1;

    float ukxqirp1          = atomJ.inducedDipolePS[1]*qir3 - atomJ.inducedDipolePS[2]*qir2;
    float ukxqirp2          = atomJ.inducedDipolePS[2]*qir1 - atomJ.inducedDipolePS[0]*qir3;
    float ukxqirp3          = atomJ.inducedDipolePS[0]*qir2 - atomJ.inducedDipolePS[1]*qir1;

    float rxqiuk1           = yr*qiuk3 - zr*qiuk2;
    float rxqiuk2           = zr*qiuk1 - xr*qiuk3;
    float rxqiuk3           = xr*qiuk2 - yr*qiuk1;

    float rxqiukp1          = yr*qiukp3 - zr*qiukp2;
    float rxqiukp2          = zr*qiukp1 - xr*qiukp3;
    float rxqiukp3          = xr*qiukp2 - yr*qiukp1;
#endif

    // get intermediate variables for induction energy terms

    float sci3              = atomI.inducedDipoleS[0]*xr + atomI.inducedDipoleS[1]*yr + atomI.inducedDipoleS[2]*zr;
    float sci4              = atomJ.inducedDipoleS[0]*xr + atomJ.inducedDipoleS[1]*yr + atomJ.inducedDipoleS[2]*zr;
#ifdef F1
    ftm2i1                 += 0.5f*scale5i*(sci4*atomI.inducedDipolePS[0] + sci3*atomJ.inducedDipolePS[0]);
    ftm2i2                 += 0.5f*scale5i*(sci4*atomI.inducedDipolePS[1] + sci3*atomJ.inducedDipolePS[1]);
    ftm2i3                 += 0.5f*scale5i*(sci4*atomI.inducedDipolePS[2] + sci3*atomJ.inducedDipolePS[2]);
#endif
    float sci3Y             = atomI.inducedDipole[0]*xr + atomI.inducedDipole[1]*yr + atomI.inducedDipole[2]*zr;
    float sci4Y             = atomJ.inducedDipole[0]*xr + atomJ.inducedDipole[1]*yr + atomJ.inducedDipole[2]*zr;
#ifdef F1
    ftm2i1                 -= 0.5f*scale5i*(sci3Y*atomJ.inducedDipoleP[0] + sci4Y*atomI.inducedDipoleP[0]);
    ftm2i2                 -= 0.5f*scale5i*(sci3Y*atomJ.inducedDipoleP[1] + sci4Y*atomI.inducedDipoleP[1]);
    ftm2i3                 -= 0.5f*scale5i*(sci3Y*atomJ.inducedDipoleP[2] + sci4Y*atomI.inducedDipoleP[2]);
#endif

    float sci7              = qir1*atomJ.inducedDipoleS[0] + qir2*atomJ.inducedDipoleS[1] + qir3*atomJ.inducedDipoleS[2];
    float sci8              = qkr1*atomI.inducedDipoleS[0] + qkr2*atomI.inducedDipoleS[1] + qkr3*atomI.inducedDipoleS[2];
    float scip1             = atomI.inducedDipolePS[0]*atomJ.labFrameDipole[0] + atomI.inducedDipolePS[1]*atomJ.labFrameDipole[1] + atomI.inducedDipolePS[2]*atomJ.labFrameDipole[2] +
                              atomI.labFrameDipole[0]*atomJ.inducedDipolePS[0] + atomI.labFrameDipole[1]*atomJ.inducedDipolePS[1] + atomI.labFrameDipole[2]*atomJ.inducedDipolePS[2];

    float scip2             = atomI.inducedDipoleS[0]*atomJ.inducedDipolePS[0] + atomI.inducedDipoleS[1]*atomJ.inducedDipolePS[1] + atomI.inducedDipoleS[2]*atomJ.inducedDipolePS[2] +
                              atomI.inducedDipolePS[0]*atomJ.inducedDipoleS[0] + atomI.inducedDipolePS[1]*atomJ.inducedDipoleS[1] + atomI.inducedDipolePS[2]*atomJ.inducedDipoleS[2];

    sci7                   -= qir1*atomJ.inducedDipole[0] + qir2*atomJ.inducedDipole[1] + qir3*atomJ.inducedDipole[2];
    sci8                   -= qkr1*atomI.inducedDipole[0] + qkr2*atomI.inducedDipole[1] + qkr3*atomI.inducedDipole[2];

    scip1                  -= atomI.inducedDipoleP[0]*atomJ.labFrameDipole[0]  + atomI.inducedDipoleP[1]*atomJ.labFrameDipole[1]  + atomI.inducedDipoleP[2]*atomJ.labFrameDipole[2] +
                              atomI.labFrameDipole[0]*atomJ.inducedDipoleP[0]  + atomI.labFrameDipole[1]*atomJ.inducedDipoleP[1]  + atomI.labFrameDipole[2]*atomJ.inducedDipoleP[2];


    scip2                  -= atomI.inducedDipole[0]*atomJ.inducedDipoleP[0]   + atomI.inducedDipole[1]*atomJ.inducedDipoleP[1]   + atomI.inducedDipole[2]*atomJ.inducedDipoleP[2]   +
                              atomI.inducedDipoleP[0]*atomJ.inducedDipole[0]   + atomI.inducedDipoleP[1]*atomJ.inducedDipole[1]   + atomI.inducedDipoleP[2]*atomJ.inducedDipole[2];


    float scip3             = atomI.inducedDipolePS[0]*xr + atomI.inducedDipolePS[1]*yr + atomI.inducedDipolePS[2]*zr;
    float scip4             = atomJ.inducedDipolePS[0]*xr + atomJ.inducedDipolePS[1]*yr + atomJ.inducedDipolePS[2]*zr;
    float gfi1              = -2.5f*(sci3*scip4+scip3*sci4)*scale5i;

#ifdef F1
    ftm2i1                 += 0.5f*scale5i*(scip4*atomI.inducedDipoleS[0] + scip3*atomJ.inducedDipoleS[0]);
    ftm2i2                 += 0.5f*scale5i*(scip4*atomI.inducedDipoleS[1] + scip3*atomJ.inducedDipoleS[1]);
    ftm2i3                 += 0.5f*scale5i*(scip4*atomI.inducedDipoleS[2] + scip3*atomJ.inducedDipoleS[2]);
#endif

    float scip3Y            = atomI.inducedDipoleP[0]*xr + atomI.inducedDipoleP[1]*yr + atomI.inducedDipoleP[2]*zr;
    float scip4Y            = atomJ.inducedDipoleP[0]*xr + atomJ.inducedDipoleP[1]*yr + atomJ.inducedDipoleP[2]*zr;
    gfi1                   += 2.5f*( sci3Y*scip4Y + scip3Y*sci4Y)*scale5i;
#ifdef F1
    ftm2i1                 -= 0.5f*scale5i*( scip3Y*atomJ.inducedDipole[0] + scip4Y*atomI.inducedDipole[0]);
    ftm2i2                 -= 0.5f*scale5i*( scip3Y*atomJ.inducedDipole[1] + scip4Y*atomI.inducedDipole[1]);
    ftm2i3                 -= 0.5f*scale5i*( scip3Y*atomJ.inducedDipole[2] + scip4Y*atomI.inducedDipole[2]);
#endif
    sci3Y                   = sci3  - sci3Y;
    sci4Y                   = sci4  - sci4Y;
    scip3Y                  = scip3 - scip3Y;
    scip4Y                  = scip4 - scip4Y;

    float scip7             = qir1*atomJ.inducedDipolePS[0] + qir2*atomJ.inducedDipolePS[1] + qir3*atomJ.inducedDipolePS[2];
    scip7                  -= qir1*atomJ.inducedDipoleP[0] + qir2*atomJ.inducedDipoleP[1] + qir3*atomJ.inducedDipoleP[2];

    float scip8             = qkr1*atomI.inducedDipolePS[0] + qkr2*atomI.inducedDipolePS[1] + qkr3*atomI.inducedDipolePS[2];
    scip8                  -= qkr1*atomI.inducedDipoleP[0] + qkr2*atomI.inducedDipoleP[1] + qkr3*atomI.inducedDipoleP[2];

    float sci1              = atomI.inducedDipoleS[0]*atomJ.labFrameDipole[0] + atomI.inducedDipoleS[1]*atomJ.labFrameDipole[1] +
                              atomI.inducedDipoleS[2]*atomJ.labFrameDipole[2] + atomI.labFrameDipole[0]*atomJ.inducedDipoleS[0] +
                              atomI.labFrameDipole[1]*atomJ.inducedDipoleS[1] + atomI.labFrameDipole[2]*atomJ.inducedDipoleS[2];
    sci1                   -= atomI.inducedDipole[0]*atomJ.labFrameDipole[0] + atomI.inducedDipole[1]*atomJ.labFrameDipole[1] +
                              atomI.inducedDipole[2]*atomJ.labFrameDipole[2] + atomI.labFrameDipole[0]*atomJ.inducedDipole[0] +
                              atomI.labFrameDipole[1]*atomJ.inducedDipole[1] + atomI.labFrameDipole[2]*atomJ.inducedDipole[2];

    float gli1              = atomJ.q*sci3Y - atomI.q*sci4Y + sci1;
    float gli2              = -sc3*sci4Y - sci3Y*sc4 + 2.0f*(sci7-sci8);
    float gli3              = sci3Y*sc6 - sci4Y*sc5;
    float glip1             = atomJ.q*scip3Y - atomI.q*scip4Y + scip1;
    float glip2             = -sc3*scip4Y - scip3Y*sc4 + 2.0f*(scip7-scip8);
    float glip3             = scip3Y*sc6 - scip4Y*sc5;

#ifdef F1
#ifdef APPLY_SCALE
    ftm2i1                 -= 0.5f*((gli1*pscale + glip1*dscale)*ddsc3_1 + (gli2*pscale + glip2*dscale)*ddsc5_1 + (gli3*pscale+glip3*dscale)*ddsc7_1);
    ftm2i2                 -= 0.5f*((gli1*pscale + glip1*dscale)*ddsc3_2 + (gli2*pscale + glip2*dscale)*ddsc5_2 + (gli3*pscale+glip3*dscale)*ddsc7_2);
    ftm2i3                 -= 0.5f*((gli1*pscale + glip1*dscale)*ddsc3_3 + (gli2*pscale + glip2*dscale)*ddsc5_3 + (gli3*pscale+glip3*dscale)*ddsc7_3);
#else
    ftm2i1                 -= 0.5f*((gli1 + glip1)*ddsc3_1 + (gli2 + glip2)*ddsc5_1 + (gli3 + glip3)*ddsc7_1);
    ftm2i2                 -= 0.5f*((gli1 + glip1)*ddsc3_2 + (gli2 + glip2)*ddsc5_2 + (gli3 + glip3)*ddsc7_2);
    ftm2i3                 -= 0.5f*((gli1 + glip1)*ddsc3_3 + (gli2 + glip2)*ddsc5_3 + (gli3 + glip3)*ddsc7_3);
#endif
    *outputEnergy           = gli1*psc3 + gli2*psc5 + gli3*psc7;
#endif

#ifdef APPLY_SCALE
    gfi1                   += 1.5f*(gli1*psc3 + glip1*dsc3);
    gfi1                   += 2.5f*(gli2*psc5 + glip2*dsc5);
    gfi1                   += 3.5f*(gli3*psc7 + glip3*dsc7);
#else
    gfi1                   += 1.5f*psc3*(gli1 + glip1);
    gfi1                   += 2.5f*psc5*(gli2 + glip2);
    gfi1                   += 3.5f*psc7*(gli3 + glip3);
#endif
    gfi1                   *= rr2;
    gfi1                   += 0.5f*scip2*scale3i;

#if defined F1 || defined T1
#ifdef APPLY_SCALE
    float gfi5              =  (sci4Y*psc7+scip4Y*dsc7);
#else
    float gfi5              =  psc7*(sci4Y + scip4Y);
#endif
#endif

#if defined F1 || defined T3
#ifdef APPLY_SCALE
    float gfi6              = -(sci3Y*psc7+scip3Y*dsc7);
#else
    float gfi6              = -psc7*(sci3Y + scip3Y );
#endif
#endif

#ifdef F1
    ftm2i1                 += gfi1*xr;

    float diff0             = atomI.inducedDipoleS[0]   - atomI.inducedDipole[0];               
    float diff1             = atomI.inducedDipolePS[0]  - atomI.inducedDipoleP[0];               
#ifdef APPLY_SCALE
    ftm2i1                 += 0.5f*(-atomJ.q*( diff0*psc3  + diff1*dsc3 ) + sc4*( diff0*psc5  + diff1*dsc5 ) - sc6*( diff0*psc7  + diff1*dsc7 ));
#else
    ftm2i1                 += 0.5f*(-atomJ.q*psc3*( diff0 + diff1 ) + sc4*psc5*( diff0  + diff1 ) - sc6*psc7*( diff0  + diff1 ));
#endif
    
    diff0                   = atomJ.inducedDipoleS[0]   - atomJ.inducedDipole[0];               
    diff1                   = atomJ.inducedDipolePS[0]  - atomJ.inducedDipoleP[0];               
#ifdef APPLY_SCALE
    ftm2i1                 += 0.5f*(atomI.q*( diff0*psc3  + diff1*dsc3 ) + sc3*( diff0*psc5  + diff1*dsc5 ) + sc5*( diff0*psc7  + diff1*dsc7 ));
    ftm2i1                 += 0.5f*(sci4Y*psc5+scip4Y*dsc5)*atomI.labFrameDipole[0] + 0.5f*(sci3Y*psc5+scip3Y*dsc5)*atomJ.labFrameDipole[0] + gfi5*qir1 + gfi6*qkr1;
#else
    ftm2i1                 += 0.5f*(atomI.q*psc3*( diff0  + diff1 ) + sc3*psc5*( diff0  + diff1 ) + sc5*psc7*( diff0  + diff1 ));
    ftm2i1                 += 0.5f*psc5*(sci4Y + scip4Y )*atomI.labFrameDipole[0] + 0.5f*psc5*(sci3Y + scip3Y)*atomJ.labFrameDipole[0] + gfi5*qir1 + gfi6*qkr1;
#endif


    ftm2i2                += gfi1*yr;

    diff0                   = atomI.inducedDipoleS[1]   - atomI.inducedDipole[1];               
    diff1                   = atomI.inducedDipolePS[1]  - atomI.inducedDipoleP[1];               
#ifdef APPLY_SCALE
    ftm2i2                 += 0.5f*(-atomJ.q*( diff0*psc3  + diff1*dsc3 ) + sc4*( diff0*psc5  + diff1*dsc5 ) - sc6*( diff0*psc7  + diff1*dsc7 ));
#else
    ftm2i2                 += 0.5f*(-atomJ.q*psc3*( diff0  + diff1 ) + sc4*psc5*( diff0 + diff1 ) - sc6*psc7*( diff0  + diff1 ));
#endif

    diff0                   = atomJ.inducedDipoleS[1]   - atomJ.inducedDipole[1];               
    diff1                   = atomJ.inducedDipolePS[1]  - atomJ.inducedDipoleP[1];               

#ifdef APPLY_SCALE
    ftm2i2                 += 0.5f*(atomI.q*( diff0*psc3  + diff1*dsc3 ) + sc3*( diff0*psc5  + diff1*dsc5 ) + sc5*( diff0*psc7  + diff1*dsc7 ));
    ftm2i2                 += 0.5f*(sci4Y*psc5+scip4Y*dsc5)*atomI.labFrameDipole[1] + 0.5f*(sci3Y*psc5+scip3Y*dsc5)*atomJ.labFrameDipole[1] + gfi5*qir2 + gfi6*qkr2;
#else
    ftm2i2                 += 0.5f*(atomI.q*psc3*( diff0  + diff1 ) + sc3*psc5*( diff0  + diff1 ) + sc5*psc7*( diff0  + diff1 ));
    ftm2i2                 += 0.5f*psc5*(sci4Y +scip4Y)*atomI.labFrameDipole[1] + 0.5f*psc5*(sci3Y +scip3Y)*atomJ.labFrameDipole[1] + gfi5*qir2 + gfi6*qkr2;
#endif


    ftm2i3                 += gfi1*zr;

    diff0                   = atomI.inducedDipoleS[2]   - atomI.inducedDipole[2];               
    diff1                   = atomI.inducedDipolePS[2]  - atomI.inducedDipoleP[2];               
#ifdef APPLY_SCALE
    ftm2i3                 += 0.5f*(-atomJ.q*( diff0*psc3  + diff1*dsc3 ) + sc4*( diff0*psc5  + diff1*dsc5 ) - sc6*( diff0*psc7  + diff1*dsc7 ));
#else
    ftm2i3                 += 0.5f*(-atomJ.q*psc3*( diff0  + diff1 ) + sc4*psc5*( diff0  + diff1 ) - sc6*psc7*( diff0  + diff1 ));
#endif

    diff0                   = atomJ.inducedDipoleS[2]   - atomJ.inducedDipole[2];               
    diff1                   = atomJ.inducedDipolePS[2]  - atomJ.inducedDipoleP[2];               

#ifdef APPLY_SCALE
    ftm2i3                 += 0.5f*(atomI.q*( diff0*psc3  + diff1*dsc3 ) + sc3*( diff0*psc5  + diff1*dsc5 ) + sc5*( diff0*psc7  + diff1*dsc7 ));
    ftm2i3                 += 0.5f*(sci4Y*psc5+scip4Y*dsc5)*atomI.labFrameDipole[2] + 0.5f*(sci3Y*psc5+scip3Y*dsc5)*atomJ.labFrameDipole[2] + gfi5*qir3 + gfi6*qkr3;
#else
    ftm2i3                 += 0.5f*(atomI.q*psc3*( diff0  + diff1 ) + sc3*psc5*( diff0  + diff1 ) + sc5*psc7*( diff0  + diff1 ));
    ftm2i3                 += 0.5f*psc5*(sci4Y + scip4Y)*atomI.labFrameDipole[2] + 0.5f*psc5*(sci3Y+scip3Y)*atomJ.labFrameDipole[2] + gfi5*qir3 + gfi6*qkr3;
#endif

 
    // intermediate values needed for partially excluded interactions

    // correction to convert mutual to direct polarization force

    if ( cAmoebaSim.polarizationType ){
        float gfd      = (scip2*scale3i - 5.0f*rr2*(scip3*sci4+sci3*scip4)*scale5i);
        float fdir1    = gfd*xr + scale5i* (sci4*atomI.inducedDipolePS[0]+scip4*atomI.inducedDipoleS[0] + sci3*atomJ.inducedDipolePS[0]+scip3*atomJ.inducedDipoleS[0]);
        float fdir2    = gfd*yr + scale5i* (sci4*atomI.inducedDipolePS[1]+scip4*atomI.inducedDipoleS[1] + sci3*atomJ.inducedDipolePS[1]+scip3*atomJ.inducedDipoleS[1]);
        float fdir3    = gfd*zr + scale5i* (sci4*atomI.inducedDipolePS[2]+scip4*atomI.inducedDipoleS[2] + sci3*atomJ.inducedDipolePS[2]+scip3*atomJ.inducedDipoleS[2]);
        ftm2i1        -= 0.5f*fdir1;
        ftm2i2        -= 0.5f*fdir2;
        ftm2i3        -= 0.5f*fdir3;

        float sci3X    = sci3  - sci3Y;
        float sci4X    = sci4  - sci4Y;
        float scip3X   = scip3 - scip3Y;
        float scip4X   = scip4 - scip4Y;
        gfd            = -5.0f*rr2*(scip3X*sci4X+sci3X*scip4X)*scale5i;
        fdir1          = gfd*xr + scale5i*(sci4X*atomI.inducedDipoleP[0] + scip4X*atomI.inducedDipole[0] + sci3X*atomJ.inducedDipoleP[0] + scip3X*atomJ.inducedDipole[0]);
        fdir2          = gfd*yr + scale5i*(sci4X*atomI.inducedDipoleP[1] + scip4X*atomI.inducedDipole[1] + sci3X*atomJ.inducedDipoleP[1] + scip3X*atomJ.inducedDipole[1]);
        fdir3          = gfd*zr + scale5i*(sci4X*atomI.inducedDipoleP[2] + scip4X*atomI.inducedDipole[2] + sci3X*atomJ.inducedDipoleP[2] + scip3X*atomJ.inducedDipole[2]);
        ftm2i1        += 0.5f*fdir1;
        ftm2i2        += 0.5f*fdir2;
        ftm2i3        += 0.5f*fdir3;
    } else {
        float findmp1  = uscale*(scip2*ddsc3_1 - ddsc5_1*(sci3*scip4+scip3*sci4));
        float findmp2  = uscale*(scip2*ddsc3_2 - ddsc5_2*(sci3*scip4+scip3*sci4));
        float findmp3  = uscale*(scip2*ddsc3_3 - ddsc5_3*(sci3*scip4+scip3*sci4));
        ftm2i1        -= 0.5f*findmp1;
        ftm2i2        -= 0.5f*findmp2;
        ftm2i3        -= 0.5f*findmp3;

        float sci3X    = sci3  - sci3Y;
        float sci4X    = sci4  - sci4Y;
        float scip3X   = scip3 - scip3Y;
        float scip4X   = scip4 - scip4Y;
        ftm2i1        += 0.5f*uscale*(-ddsc5_1*(sci3X*scip4X+scip3X*sci4X));
        ftm2i2        += 0.5f*uscale*(-ddsc5_2*(sci3X*scip4X+scip3X*sci4X));
        ftm2i3        += 0.5f*uscale*(-ddsc5_3*(sci3X*scip4X+scip3X*sci4X)); 
    }
#endif

#ifdef T1
#ifdef APPLY_SCALE
    float gti2              = 0.5f*(sci4Y*psc5 + scip4Y*dsc5);
    float ttm2i1            = -(dixuk1*psc3+dixukp1*dsc3)*0.5f + gti2*dixr1 + ((ukxqir1+rxqiuk1)*psc5 +(ukxqirp1+rxqiukp1)*dsc5) - gfi5*rxqir1;
    float ttm2i2            = -(dixuk2*psc3+dixukp2*dsc3)*0.5f + gti2*dixr2 + ((ukxqir2+rxqiuk2)*psc5 +(ukxqirp2+rxqiukp2)*dsc5) - gfi5*rxqir2;
    float ttm2i3            = -(dixuk3*psc3+dixukp3*dsc3)*0.5f + gti2*dixr3 + ((ukxqir3+rxqiuk3)*psc5 +(ukxqirp3+rxqiukp3)*dsc5) - gfi5*rxqir3;
#else
    float gti2              = 0.5f*psc5*(sci4Y + scip4Y);
    float ttm2i1            = -psc3*(dixuk1 + dixukp1)*0.5f + gti2*dixr1 + psc5*((ukxqir1+rxqiuk1) + (ukxqirp1+rxqiukp1)) - gfi5*rxqir1;
    float ttm2i2            = -psc3*(dixuk2 + dixukp2)*0.5f + gti2*dixr2 + psc5*((ukxqir2+rxqiuk2) + (ukxqirp2+rxqiukp2)) - gfi5*rxqir2;
    float ttm2i3            = -psc3*(dixuk3 + dixukp3)*0.5f + gti2*dixr3 + psc5*((ukxqir3+rxqiuk3) + (ukxqirp3+rxqiukp3)) - gfi5*rxqir3;
#endif
#endif

#ifdef T3
#ifdef APPLY_SCALE
    float gti3              = 0.5f*(sci3Y*psc5 + scip3Y*dsc5);
    float ttm3i1            = -(dkxui1*psc3+dkxuip1*dsc3)*0.5f + gti3*dkxr1 - ((uixqkr1+rxqkui1)*psc5 +(uixqkrp1+rxqkuip1)*dsc5) - gfi6*rxqkr1;
    float ttm3i2            = -(dkxui2*psc3+dkxuip2*dsc3)*0.5f + gti3*dkxr2 - ((uixqkr2+rxqkui2)*psc5 +(uixqkrp2+rxqkuip2)*dsc5) - gfi6*rxqkr2;
    float ttm3i3            = -(dkxui3*psc3+dkxuip3*dsc3)*0.5f + gti3*dkxr3 - ((uixqkr3+rxqkui3)*psc5 +(uixqkrp3+rxqkuip3)*dsc5) - gfi6*rxqkr3;
#else
    float gti3              = 0.5f*psc5*(sci3Y + scip3Y);
    float ttm3i1            = -psc3*(dkxui1 + dkxuip1)*0.5f + gti3*dkxr1 - psc5*((uixqkr1+rxqkui1) + (uixqkrp1+rxqkuip1)) - gfi6*rxqkr1;
    float ttm3i2            = -psc3*(dkxui2 + dkxuip2)*0.5f + gti3*dkxr2 - psc5*((uixqkr2+rxqkui2) + (uixqkrp2+rxqkuip2)) - gfi6*rxqkr2;
    float ttm3i3            = -psc3*(dkxui3 + dkxuip3)*0.5f + gti3*dkxr3 - psc5*((uixqkr3+rxqkui3) + (uixqkrp3+rxqkuip3)) - gfi6*rxqkr3;
#endif
#endif
 
    // update force and torque on site k
    
#ifdef F1
    outputForce[0]      = -ftm2i1;
    outputForce[1]      = -ftm2i2;
    outputForce[2]      = -ftm2i3;
#endif

#ifdef T1
    outputForce[0]    = ttm2i1;
    outputForce[1]    = ttm2i2;
    outputForce[2]    = ttm2i3;
#endif

#ifdef T3
    outputForce[0]    = ttm3i1;
    outputForce[1]    = ttm3i2;
    outputForce[2]    = ttm3i3;
#endif

    // construct auxiliary vectors for induced terms

#ifdef T1
    dixuk1            = atomI.labFrameDipole[1]*atomJ.inducedDipole[2] - atomI.labFrameDipole[2]*atomJ.inducedDipole[1];
    dixuk2            = atomI.labFrameDipole[2]*atomJ.inducedDipole[0] - atomI.labFrameDipole[0]*atomJ.inducedDipole[2];
    dixuk3            = atomI.labFrameDipole[0]*atomJ.inducedDipole[1] - atomI.labFrameDipole[1]*atomJ.inducedDipole[0];

    dixukp1           = atomI.labFrameDipole[1]*atomJ.inducedDipoleP[2] - atomI.labFrameDipole[2]*atomJ.inducedDipoleP[1];
    dixukp2           = atomI.labFrameDipole[2]*atomJ.inducedDipoleP[0] - atomI.labFrameDipole[0]*atomJ.inducedDipoleP[2];
    dixukp3           = atomI.labFrameDipole[0]*atomJ.inducedDipoleP[1] - atomI.labFrameDipole[1]*atomJ.inducedDipoleP[0];
#endif

#ifdef T3
    dkxui1            = atomJ.labFrameDipole[1]*atomI.inducedDipole[2] - atomJ.labFrameDipole[2]*atomI.inducedDipole[1];
    dkxui2            = atomJ.labFrameDipole[2]*atomI.inducedDipole[0] - atomJ.labFrameDipole[0]*atomI.inducedDipole[2];
    dkxui3            = atomJ.labFrameDipole[0]*atomI.inducedDipole[1] - atomJ.labFrameDipole[1]*atomI.inducedDipole[0];

    dkxuip1           = atomJ.labFrameDipole[1]*atomI.inducedDipoleP[2] - atomJ.labFrameDipole[2]*atomI.inducedDipoleP[1];
    dkxuip2           = atomJ.labFrameDipole[2]*atomI.inducedDipoleP[0] - atomJ.labFrameDipole[0]*atomI.inducedDipoleP[2];
    dkxuip3           = atomJ.labFrameDipole[0]*atomI.inducedDipoleP[1] - atomJ.labFrameDipole[1]*atomI.inducedDipoleP[0];
#endif

#if defined T1
    qiuk1             = atomI.labFrameQuadrupole_XX*atomJ.inducedDipole[0] + atomI.labFrameQuadrupole_XY*atomJ.inducedDipole[1] + atomI.labFrameQuadrupole_XZ*atomJ.inducedDipole[2];
    qiuk2             = atomI.labFrameQuadrupole_XY*atomJ.inducedDipole[0] + atomI.labFrameQuadrupole_YY*atomJ.inducedDipole[1] + atomI.labFrameQuadrupole_YZ*atomJ.inducedDipole[2];
    qiuk3             = atomI.labFrameQuadrupole_XZ*atomJ.inducedDipole[0] + atomI.labFrameQuadrupole_YZ*atomJ.inducedDipole[1] + atomI.labFrameQuadrupole_ZZ*atomJ.inducedDipole[2];

    qiukp1            = atomI.labFrameQuadrupole_XX*atomJ.inducedDipoleP[0] + atomI.labFrameQuadrupole_XY*atomJ.inducedDipoleP[1] + atomI.labFrameQuadrupole_XZ*atomJ.inducedDipoleP[2];
    qiukp2            = atomI.labFrameQuadrupole_XY*atomJ.inducedDipoleP[0] + atomI.labFrameQuadrupole_YY*atomJ.inducedDipoleP[1] + atomI.labFrameQuadrupole_YZ*atomJ.inducedDipoleP[2];
    qiukp3            = atomI.labFrameQuadrupole_XZ*atomJ.inducedDipoleP[0] + atomI.labFrameQuadrupole_YZ*atomJ.inducedDipoleP[1] + atomI.labFrameQuadrupole_ZZ*atomJ.inducedDipoleP[2];
#endif

#if defined T3
    qkui1             = atomJ.labFrameQuadrupole_XX*atomI.inducedDipole[0] + atomJ.labFrameQuadrupole_XY*atomI.inducedDipole[1] + atomJ.labFrameQuadrupole_XZ*atomI.inducedDipole[2];
    qkui2             = atomJ.labFrameQuadrupole_XY*atomI.inducedDipole[0] + atomJ.labFrameQuadrupole_YY*atomI.inducedDipole[1] + atomJ.labFrameQuadrupole_YZ*atomI.inducedDipole[2];
    qkui3             = atomJ.labFrameQuadrupole_XZ*atomI.inducedDipole[0] + atomJ.labFrameQuadrupole_YZ*atomI.inducedDipole[1] + atomJ.labFrameQuadrupole_ZZ*atomI.inducedDipole[2];

    qkuip1            = atomJ.labFrameQuadrupole_XX*atomI.inducedDipoleP[0] + atomJ.labFrameQuadrupole_XY*atomI.inducedDipoleP[1] + atomJ.labFrameQuadrupole_XZ*atomI.inducedDipoleP[2];
    qkuip2            = atomJ.labFrameQuadrupole_XY*atomI.inducedDipoleP[0] + atomJ.labFrameQuadrupole_YY*atomI.inducedDipoleP[1] + atomJ.labFrameQuadrupole_YZ*atomI.inducedDipoleP[2];
    qkuip3            = atomJ.labFrameQuadrupole_XZ*atomI.inducedDipoleP[0] + atomJ.labFrameQuadrupole_YZ*atomI.inducedDipoleP[1] + atomJ.labFrameQuadrupole_ZZ*atomI.inducedDipoleP[2];
#endif

#ifdef T3
    uixqkr1           = atomI.inducedDipole[1]*qkr3 - atomI.inducedDipole[2]*qkr2;
    uixqkr2           = atomI.inducedDipole[2]*qkr1 - atomI.inducedDipole[0]*qkr3;
    uixqkr3           = atomI.inducedDipole[0]*qkr2 - atomI.inducedDipole[1]*qkr1;

    uixqkrp1          = atomI.inducedDipoleP[1]*qkr3 - atomI.inducedDipoleP[2]*qkr2;
    uixqkrp2          = atomI.inducedDipoleP[2]*qkr1 - atomI.inducedDipoleP[0]*qkr3;
    uixqkrp3          = atomI.inducedDipoleP[0]*qkr2 - atomI.inducedDipoleP[1]*qkr1;
#endif

#ifdef T1
    ukxqir1           = atomJ.inducedDipole[1]*qir3 - atomJ.inducedDipole[2]*qir2;
    ukxqir2           = atomJ.inducedDipole[2]*qir1 - atomJ.inducedDipole[0]*qir3;
    ukxqir3           = atomJ.inducedDipole[0]*qir2 - atomJ.inducedDipole[1]*qir1;

    ukxqirp1          = atomJ.inducedDipoleP[1]*qir3 - atomJ.inducedDipoleP[2]*qir2;
    ukxqirp2          = atomJ.inducedDipoleP[2]*qir1 - atomJ.inducedDipoleP[0]*qir3;
    ukxqirp3          = atomJ.inducedDipoleP[0]*qir2 - atomJ.inducedDipoleP[1]*qir1;

    rxqiuk1           = yr*qiuk3 - zr*qiuk2;
    rxqiuk2           = zr*qiuk1 - xr*qiuk3;
    rxqiuk3           = xr*qiuk2 - yr*qiuk1;

    rxqiukp1          = yr*qiukp3 - zr*qiukp2;
    rxqiukp2          = zr*qiukp1 - xr*qiukp3;
    rxqiukp3          = xr*qiukp2 - yr*qiukp1;
#endif

#ifdef T3
    rxqkui1           = yr*qkui3 - zr*qkui2;
    rxqkui2           = zr*qkui1 - xr*qkui3;
    rxqkui3           = xr*qkui2 - yr*qkui1;

    rxqkuip1          = yr*qkuip3 - zr*qkuip2;
    rxqkuip2          = zr*qkuip1 - xr*qkuip3;
    rxqkuip3          = xr*qkuip2 - yr*qkuip1;
#endif

#ifdef T1
#ifdef APPLY_SCALE
     ttm2i1           = -(dixuk1*psc3+dixukp1*dsc3)*0.5f + ((ukxqir1+rxqiuk1)*psc5 +(ukxqirp1+rxqiukp1)*dsc5);
     ttm2i2           = -(dixuk2*psc3+dixukp2*dsc3)*0.5f + ((ukxqir2+rxqiuk2)*psc5 +(ukxqirp2+rxqiukp2)*dsc5);
     ttm2i3           = -(dixuk3*psc3+dixukp3*dsc3)*0.5f + ((ukxqir3+rxqiuk3)*psc5 +(ukxqirp3+rxqiukp3)*dsc5);
#else
     ttm2i1           = -psc3*(dixuk1+dixukp1)*0.5f + psc5*((ukxqir1+rxqiuk1) + (ukxqirp1+rxqiukp1));
     ttm2i2           = -psc3*(dixuk2+dixukp2)*0.5f + psc5*((ukxqir2+rxqiuk2) + (ukxqirp2+rxqiukp2));
     ttm2i3           = -psc3*(dixuk3+dixukp3)*0.5f + psc5*((ukxqir3+rxqiuk3) + (ukxqirp3+rxqiukp3));
#endif
#endif
 
#ifdef T3
#ifdef APPLY_SCALE
     ttm3i1           = -(dkxui1*psc3+dkxuip1*dsc3)*0.5f - ((uixqkr1+rxqkui1)*psc5 +(uixqkrp1+rxqkuip1)*dsc5);
     ttm3i2           = -(dkxui2*psc3+dkxuip2*dsc3)*0.5f - ((uixqkr2+rxqkui2)*psc5 +(uixqkrp2+rxqkuip2)*dsc5);
     ttm3i3           = -(dkxui3*psc3+dkxuip3*dsc3)*0.5f - ((uixqkr3+rxqkui3)*psc5 +(uixqkrp3+rxqkuip3)*dsc5);
#else
     ttm3i1           = -psc3*(dkxui1 + dkxuip1)*0.5f - psc5*((uixqkr1+rxqkui1) + (uixqkrp1+rxqkuip1));
     ttm3i2           = -psc3*(dkxui2 + dkxuip2)*0.5f - psc5*((uixqkr2+rxqkui2) + (uixqkrp2+rxqkuip2));
     ttm3i3           = -psc3*(dkxui3 + dkxuip3)*0.5f - psc5*((uixqkr3+rxqkui3) + (uixqkrp3+rxqkuip3));
#endif
#endif

    // update force and torque on site k;

#ifdef T1
    outputForce[0]   -= ttm2i1;
    outputForce[1]   -= ttm2i2;
    outputForce[2]   -= ttm2i3;
#endif

#ifdef T3
    outputForce[0]   -= ttm3i1;
    outputForce[1]   -= ttm3i2;
    outputForce[2]   -= ttm3i3;
#endif

}

