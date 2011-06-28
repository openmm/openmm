static __device__ void SUB_METHOD_NAME( calculateKirkwoodPairIxn, _kernel)( KirkwoodParticle& atomI, KirkwoodParticle& atomJ
#if defined F1 || defined F2
                                                                            , float* outputEnergy, float force[3]
#elif defined T1 || defined T2
#ifndef INCLUDE_TORQUE

                                                                            , float torque[3]
#endif
#endif
 ){

    // set the bulk dielectric constant to the water value

    float fc           = cAmoebaSim.electric*cAmoebaSim.fc;
    float fd           = cAmoebaSim.electric*cAmoebaSim.fd;
    float fq           = cAmoebaSim.electric*cAmoebaSim.fq;

#if defined F2 || defined B2
    float sxi          = atomI.inducedDipole[0] + atomI.inducedDipoleP[0];
    float syi          = atomI.inducedDipole[1] + atomI.inducedDipoleP[1];
    float szi          = atomI.inducedDipole[2] + atomI.inducedDipoleP[2];
#endif

#if defined F2 || defined T2 || defined B2
    float sxk          = atomJ.inducedDipole[0] + atomJ.inducedDipoleP[0];
    float syk          = atomJ.inducedDipole[1] + atomJ.inducedDipoleP[1];
    float szk          = atomJ.inducedDipole[2] + atomJ.inducedDipoleP[2];
#endif

    // decide whether to compute the current interaction;

    float xr           = atomJ.x - atomI.x;
    float yr           = atomJ.y - atomI.y;
    float zr           = atomJ.z - atomI.z;

    float xr2          = xr*xr;
    float yr2          = yr*yr;
    float zr2          = zr*zr;
    float r2           = xr2 + yr2 + zr2;

    //if( r2 > cAmoebaSim.scalingDistanceCutoff ){
    //}

    float rb2          = atomI.bornRadius*atomJ.bornRadius;

    float expterm      = expf(-r2/(cAmoebaSim.gkc*rb2));
    float expc         = expterm/cAmoebaSim.gkc;
    float expcr        = r2*expterm/(cAmoebaSim.gkc*cAmoebaSim.gkc*rb2*rb2);
    float dexpc        = -2.0f / (cAmoebaSim.gkc*rb2);
    float dexpcr       = 2.0f / (cAmoebaSim.gkc*rb2*rb2);
    float dgfdr        = 0.5f*expterm*(1.0f + r2/(rb2*cAmoebaSim.gkc));
    float gf2          = 1.0f / (r2 + rb2*expterm);

    float gf           = sqrt(gf2);
    float gf3          = gf2*gf;
    float gf5          = gf3*gf2;
    float gf7          = gf5*gf2;
    float gf9          = gf7*gf2;
    float gf11         = gf9*gf2;

    // reaction potential auxiliary terms;

    float a00          =         gf;
    float a10          =        -gf3;
    float a20          =    3.0f*gf5;
    float a30          =  -15.0f*gf7;
    float a40          =  105.0f*gf9;
    float a50          = -945.0f*gf11;

    // Born radii derivatives of reaction potential auxiliary terms;

    float b00          = dgfdr*a10;
    float b10          = dgfdr*a20;
    float b20          = dgfdr*a30;
    float b30          = dgfdr*a40;
    float b40          = dgfdr*a50;

    // reaction potential gradient auxiliary terms;

    float expc1        = 1.0f - expc;
    float a01          = expc1*a10;
    float a11          = expc1*a20;
    float a21          = expc1*a30;
    float a31          = expc1*a40;
    float a41          = expc1*a50;

    // Born radii derivs of reaction potential gradient auxiliary terms;

    float b01          = b10 - expcr*a10 - expc*b10;
    float b11          = b20 - expcr*a20 - expc*b20;
    float b21          = b30 - expcr*a30 - expc*b30;
    float b31          = b40 - expcr*a40 - expc*b40;

    // 2nd reaction potential gradient auxiliary terms;

    float expcdexpc    = -expc*dexpc;
    float a02          = expc1*a11  + expcdexpc*a10;
    float a12          = expc1*a21  + expcdexpc*a20;
    float a22          = expc1*a31  + expcdexpc*a30;
    float a32          = expc1*a41  + expcdexpc*a40;

    // Born radii derivatives of the 2nd reaction potential
    // gradient auxiliary terms

     float b02         = b11 - (expcr*(a11  + dexpc*a10) + expc*(b11  + dexpcr*a10 + dexpc*b10));
     float b12         = b21 - (expcr*(a21  + dexpc*a20) + expc*(b21  + dexpcr*a20 + dexpc*b20));
     float b22         = b31 - (expcr*(a31  + dexpc*a30) + expc*(b31  + dexpcr*a30 + dexpc*b30));

    // 3rd reaction potential gradient auxiliary terms

    expcdexpc          = 2.0f*expcdexpc;
    float a03          = expc1*a12  + expcdexpc*a11;
    float a13          = expc1*a22  + expcdexpc*a21;
    float a23          = expc1*a32  + expcdexpc*a31;

    expcdexpc          = -expc*dexpc*dexpc;
    a03                = a03  + expcdexpc*a10;
    a13                = a13  + expcdexpc*a20;
    a23                = a23  + expcdexpc*a30;

    // multiply the auxillary terms by their dieletric functions;

    a00            *= fc;
    a01            *= fc;
    a02            *= fc;
    a03            *= fc;

    b00            *= fc;
    b01            *= fc;
    b02            *= fc;

    a10            *= fd;
    a11            *= fd;
    a12            *= fd;
    a13            *= fd;

    b10            *= fd;
    b11            *= fd;
    b12            *= fd;

    a20            *= fq;
    a21            *= fq;
    a22            *= fq;
    a23            *= fq;

    b20            *= fq;
    b21            *= fq;
    b22            *= fq;

    // unweighted reaction potential tensor

#if defined F2
    float energy     = -a10*atomJ.q*(atomI.inducedDipole[0]*xr + atomI.inducedDipole[1]*yr + atomI.inducedDipole[2]*zr);
    energy          +=  a10*atomI.q*(atomJ.inducedDipole[0]*xr + atomJ.inducedDipole[1]*yr + atomJ.inducedDipole[2]*zr);
#endif

#if defined F1
    float energy     = 2.0f*atomI.q*atomJ.q*a00;
    energy          += -a10*atomJ.q*(atomI.labFrameDipole[0]*xr + atomI.labFrameDipole[1]*yr + atomI.labFrameDipole[2]*zr);
    energy          +=  a10*atomI.q*(atomJ.labFrameDipole[0]*xr + atomJ.labFrameDipole[1]*yr + atomJ.labFrameDipole[2]*zr);
    energy          +=  a20*atomJ.q*(atomI.labFrameQuadrupole_XX*xr2 + atomI.labFrameQuadrupole_YY*yr2 + atomI.labFrameQuadrupole_ZZ*zr2 + 2.0f*(atomI.labFrameQuadrupole_XY*xr*yr + atomI.labFrameQuadrupole_XZ*xr*zr + atomI.labFrameQuadrupole_YZ*yr*zr));
    energy          +=  a20*atomI.q*(atomJ.labFrameQuadrupole_XX*xr2 + atomJ.labFrameQuadrupole_YY*yr2 + atomJ.labFrameQuadrupole_ZZ*zr2 + 2.0f*(atomJ.labFrameQuadrupole_XY*xr*yr + atomJ.labFrameQuadrupole_XZ*xr*zr + atomJ.labFrameQuadrupole_YZ*yr*zr));
#endif

    // Born radii derivs of unweighted reaction potential tensor

#if defined B1
    float dsumdrB1   = 2.0f*(atomI.q*atomJ.q*b00);
    dsumdrB1        -= b10*atomJ.q*(atomI.labFrameDipole[0]*xr + atomI.labFrameDipole[1]*yr + atomI.labFrameDipole[2]*zr);
    dsumdrB1        += b10*atomI.q*(atomJ.labFrameDipole[0]*xr + atomJ.labFrameDipole[1]*yr + atomJ.labFrameDipole[2]*zr);
#endif
#if defined B2
    float dsumdrB2   = -b10*atomJ.q*(sxi*xr + syi*yr + szi*zr);
    dsumdrB2        +=  b10*atomI.q*(sxk*xr + syk*yr + szk*zr);
#endif

#if defined B1
    float gqxx21     = xr*xr;
    float gqyy21     = yr*yr;
    float gqzz21     = zr*zr;

    float gqxy21     = xr*yr;
    float gqxz21     = xr*zr;
    float gqyz21     = yr*zr;
    dsumdrB1        += b20*atomJ.q*(atomI.labFrameQuadrupole_XX*gqxx21 + atomI.labFrameQuadrupole_YY*gqyy21 + atomI.labFrameQuadrupole_ZZ*gqzz21 + 2.0f*(atomI.labFrameQuadrupole_XY*gqxy21 + atomI.labFrameQuadrupole_XZ*gqxz21 + atomI.labFrameQuadrupole_YZ*gqyz21));
    dsumdrB1        += b20*atomI.q*(atomJ.labFrameQuadrupole_XX*gqxx21 + atomJ.labFrameQuadrupole_YY*gqyy21 + atomJ.labFrameQuadrupole_ZZ*gqzz21 + 2.0f*(atomJ.labFrameQuadrupole_XY*gqxy21 + atomJ.labFrameQuadrupole_XZ*gqxz21 + atomJ.labFrameQuadrupole_YZ*gqyz21));
#endif

#if defined F1
    energy          +=  a01*atomI.q*(atomJ.labFrameDipole[0]*xr + atomJ.labFrameDipole[1]*yr + atomJ.labFrameDipole[2]*zr);
    energy          -=  a01*atomJ.q*(atomI.labFrameDipole[0]*xr + atomI.labFrameDipole[1]*yr + atomI.labFrameDipole[2]*zr);
    float factor     = a01*2.0f*atomI.q*atomJ.q;
    float dedx       = factor*xr;
    float dedy       = factor*yr;
    float dedz       = factor*zr;
#endif
#if defined F2 
    energy          += a01*atomI.q*(atomJ.inducedDipole[0]*xr + atomJ.inducedDipole[1]*yr + atomJ.inducedDipole[2]*zr);
    energy          -= a01*atomJ.q*(atomI.inducedDipole[0]*xr + atomI.inducedDipole[1]*yr + atomI.inducedDipole[2]*zr);
#endif

#if defined F1 || defined F2 || defined T1 || defined T2
    float gux2       = a10  + xr*xr*a11;
    float gux3       = xr*yr*a11;
    float gux4       = xr*zr*a11;

    float guy3       = a10  + yr*yr*a11;
    float guy4       = yr*zr*a11;
    float guz4       = a10  + zr*zr*a11;
#if defined T1
    float guy2       = gux3;
    float guz2       = gux4;
    float guz3       = guy4;
#endif
#if defined T2
    float fid1                = sxk*gux2 + syk*gux3 + szk*gux4;
    float fid2                = sxk*gux3 + syk*guy3 + szk*guy4;
    float fid3                = sxk*gux4 + syk*guy4 + szk*guz4;

    float trqi1               = atomI.labFrameDipole[1]*fid3 - atomI.labFrameDipole[2]*fid2;
    float trqi2               = atomI.labFrameDipole[2]*fid1 - atomI.labFrameDipole[0]*fid3;
    float trqi3               = atomI.labFrameDipole[0]*fid2 - atomI.labFrameDipole[1]*fid1;
#endif

#if defined F1
/*
    float sum        = atomI.labFrameDipole[0]*(atomJ.labFrameDipole[0]*gux2 + atomJ.labFrameDipole[1]*gux3 + atomJ.labFrameDipole[2]*gux4);
    sum             += atomI.labFrameDipole[1]*(atomJ.labFrameDipole[0]*gux3 + atomJ.labFrameDipole[1]*guy3 + atomJ.labFrameDipole[2]*guy4);
    sum             += atomI.labFrameDipole[2]*(atomJ.labFrameDipole[0]*gux4 + atomJ.labFrameDipole[1]*guy4 + atomJ.labFrameDipole[2]*guz4);
    energy          += -2.0f*sum;
*/
    energy          -= 2.0f*(atomI.labFrameDipole[0]*(atomJ.labFrameDipole[0]*gux2 + atomJ.labFrameDipole[1]*gux3 + atomJ.labFrameDipole[2]*gux4) +
                             atomI.labFrameDipole[1]*(atomJ.labFrameDipole[0]*gux3 + atomJ.labFrameDipole[1]*guy3 + atomJ.labFrameDipole[2]*guy4) +
                             atomI.labFrameDipole[2]*(atomJ.labFrameDipole[0]*gux4 + atomJ.labFrameDipole[1]*guy4 + atomJ.labFrameDipole[2]*guz4) );

    dedx            -= atomJ.q*(atomI.labFrameDipole[0]*gux2 + atomI.labFrameDipole[1]*gux3 + atomI.labFrameDipole[2]*gux4);
    dedx            += atomI.q*(atomJ.labFrameDipole[0]*gux2 + atomJ.labFrameDipole[1]*gux3 + atomJ.labFrameDipole[2]*gux4);

    dedy            -= atomJ.q*(atomI.labFrameDipole[0]*gux3 + atomI.labFrameDipole[1]*guy3 + atomI.labFrameDipole[2]*guy4);
    dedy            += atomI.q*(atomJ.labFrameDipole[0]*gux3 + atomJ.labFrameDipole[1]*guy3 + atomJ.labFrameDipole[2]*guy4);

    dedz            -= atomJ.q*(atomI.labFrameDipole[0]*gux4 + atomI.labFrameDipole[1]*guy4 + atomI.labFrameDipole[2]*guz4);
    dedz            += atomI.q*(atomJ.labFrameDipole[0]*gux4 + atomJ.labFrameDipole[1]*guy4 + atomJ.labFrameDipole[2]*guz4);
#endif
#if defined F2
    energy          -= 2.0f*(
                       atomI.labFrameDipole[0]*(atomJ.inducedDipole[0]*gux2 + atomJ.inducedDipole[1]*gux3 + atomJ.inducedDipole[2]*gux4)  + 
                       atomI.labFrameDipole[1]*(atomJ.inducedDipole[0]*gux3 + atomJ.inducedDipole[1]*guy3 + atomJ.inducedDipole[2]*guy4)  + 
                       atomI.labFrameDipole[2]*(atomJ.inducedDipole[0]*gux4 + atomJ.inducedDipole[1]*guy4 + atomJ.inducedDipole[2]*guz4)  + 
                       atomJ.labFrameDipole[0]*(atomI.inducedDipole[0]*gux2 + atomI.inducedDipole[1]*gux3 + atomI.inducedDipole[2]*gux4)  + 
                       atomJ.labFrameDipole[1]*(atomI.inducedDipole[0]*gux3 + atomI.inducedDipole[1]*guy3 + atomI.inducedDipole[2]*guy4)  + 
                       atomJ.labFrameDipole[2]*(atomI.inducedDipole[0]*gux4 + atomI.inducedDipole[1]*guy4 + atomI.inducedDipole[2]*guz4) );

    float dpdx       = atomI.q*(sxk*gux2 + syk*gux3 + szk*gux4);
    dpdx            -= atomJ.q*(sxi*gux2 + syi*gux3 + szi*gux4);

    float dpdy       = atomI.q*(sxk*gux3 + syk*guy3 + szk*guy4);
    dpdy            -= atomJ.q*(sxi*gux3 + syi*guy3 + szi*guy4);

    float dpdz       = atomI.q*(sxk*gux4 + syk*guy4 + szk*guz4);
    dpdz            -= atomJ.q*(sxi*gux4 + syi*guy4 + szi*guz4);

#endif
    float gqxx2      = xr*(2.0f*a20 + xr*xr*a21);
    float gqxx3      = yr*xr*xr*a21;
    float gqxx4      = zr*xr*xr*a21;
    float gqyy2      = xr*yr*yr*a21;
    float gqyy3      = yr*(2.0f*a20 + yr*yr*a21);
    float gqyy4      = zr*yr*yr*a21;
    float gqzz2      = xr*zr*zr*a21;
    float gqzz3      = yr*zr*zr*a21;
    float gqzz4      = zr*(2.0f*a20 + zr*zr*a21);
    float gqxy2      = yr*(a20 + xr*xr*a21);
    float gqxy3      = xr*(a20 + yr*yr*a21);
    float gqxy4      = zr*xr*yr*a21;
    float gqxz2      = zr*(a20 + xr*xr*a21);
    float gqxz4      = xr*(a20 + zr*zr*a21);
    float gqyz3      = zr*(a20 + yr*yr*a21);
    float gqyz4      = yr*(a20 + zr*zr*a21);
#if defined T1 || defined T2
    float gqxz3      = gqxy4;
    float gqyz2      = gqxy4;
#endif

#if defined F1
    energy          += atomJ.labFrameDipole[0]*(atomI.labFrameQuadrupole_XX*gqxx2 + atomI.labFrameQuadrupole_YY*gqyy2 + atomI.labFrameQuadrupole_ZZ*gqzz2 + 2.0f*(atomI.labFrameQuadrupole_XY*gqxy2 + atomI.labFrameQuadrupole_XZ*gqxz2 + atomI.labFrameQuadrupole_YZ*gqxy4)) +
                       atomJ.labFrameDipole[1]*(atomI.labFrameQuadrupole_XX*gqxx3 + atomI.labFrameQuadrupole_YY*gqyy3 + atomI.labFrameQuadrupole_ZZ*gqzz3 + 2.0f*(atomI.labFrameQuadrupole_XY*gqxy3 + atomI.labFrameQuadrupole_XZ*gqxy4 + atomI.labFrameQuadrupole_YZ*gqyz3)) +
                       atomJ.labFrameDipole[2]*(atomI.labFrameQuadrupole_XX*gqxx4 + atomI.labFrameQuadrupole_YY*gqyy4 + atomI.labFrameQuadrupole_ZZ*gqzz4 + 2.0f*(atomI.labFrameQuadrupole_XY*gqxy4 + atomI.labFrameQuadrupole_XZ*gqxz4 + atomI.labFrameQuadrupole_YZ*gqyz4));
    energy          -= atomI.labFrameDipole[0]*(atomJ.labFrameQuadrupole_XX*gqxx2 + atomJ.labFrameQuadrupole_YY*gqyy2 + atomJ.labFrameQuadrupole_ZZ*gqzz2 + 2.0f*(atomJ.labFrameQuadrupole_XY*gqxy2 + atomJ.labFrameQuadrupole_XZ*gqxz2 + atomJ.labFrameQuadrupole_YZ*gqxy4)) +
                       atomI.labFrameDipole[1]*(atomJ.labFrameQuadrupole_XX*gqxx3 + atomJ.labFrameQuadrupole_YY*gqyy3 + atomJ.labFrameQuadrupole_ZZ*gqzz3 + 2.0f*(atomJ.labFrameQuadrupole_XY*gqxy3 + atomJ.labFrameQuadrupole_XZ*gqxy4 + atomJ.labFrameQuadrupole_YZ*gqyz3)) +
                       atomI.labFrameDipole[2]*(atomJ.labFrameQuadrupole_XX*gqxx4 + atomJ.labFrameQuadrupole_YY*gqyy4 + atomJ.labFrameQuadrupole_ZZ*gqzz4 + 2.0f*(atomJ.labFrameQuadrupole_XY*gqxy4 + atomJ.labFrameQuadrupole_XZ*gqxz4 + atomJ.labFrameQuadrupole_YZ*gqyz4));

    dedx            += atomJ.q*(atomI.labFrameQuadrupole_XX*gqxx2 + atomI.labFrameQuadrupole_YY*gqyy2 + atomI.labFrameQuadrupole_ZZ*gqzz2 + 2.0f*(atomI.labFrameQuadrupole_XY*gqxy2 + atomI.labFrameQuadrupole_XZ*gqxz2 + atomI.labFrameQuadrupole_YZ*gqxy4));
    dedx            += atomI.q*(atomJ.labFrameQuadrupole_XX*gqxx2 + atomJ.labFrameQuadrupole_YY*gqyy2 + atomJ.labFrameQuadrupole_ZZ*gqzz2 + 2.0f*(atomJ.labFrameQuadrupole_XY*gqxy2 + atomJ.labFrameQuadrupole_XZ*gqxz2 + atomJ.labFrameQuadrupole_YZ*gqxy4));

    dedy            += atomJ.q*(atomI.labFrameQuadrupole_XX*gqxx3 + atomI.labFrameQuadrupole_YY*gqyy3 + atomI.labFrameQuadrupole_ZZ*gqzz3 + 2.0f*(atomI.labFrameQuadrupole_XY*gqxy3 + atomI.labFrameQuadrupole_XZ*gqxy4 + atomI.labFrameQuadrupole_YZ*gqyz3));
    dedy            += atomI.q*(atomJ.labFrameQuadrupole_XX*gqxx3 + atomJ.labFrameQuadrupole_YY*gqyy3 + atomJ.labFrameQuadrupole_ZZ*gqzz3 + 2.0f*(atomJ.labFrameQuadrupole_XY*gqxy3 + atomJ.labFrameQuadrupole_XZ*gqxy4 + atomJ.labFrameQuadrupole_YZ*gqyz3));

    dedz            += atomJ.q*(atomI.labFrameQuadrupole_XX*gqxx4 + atomI.labFrameQuadrupole_YY*gqyy4 + atomI.labFrameQuadrupole_ZZ*gqzz4 + 2.0f*(atomI.labFrameQuadrupole_XY*gqxy4 + atomI.labFrameQuadrupole_XZ*gqxz4 + atomI.labFrameQuadrupole_YZ*gqyz4));
    dedz            += atomI.q*(atomJ.labFrameQuadrupole_XX*gqxx4 + atomJ.labFrameQuadrupole_YY*gqyy4 + atomJ.labFrameQuadrupole_ZZ*gqzz4 + 2.0f*(atomJ.labFrameQuadrupole_XY*gqxy4 + atomJ.labFrameQuadrupole_XZ*gqxz4 + atomJ.labFrameQuadrupole_YZ*gqyz4));
#endif

#if defined F2
    energy += atomJ.inducedDipole[0]*(atomI.labFrameQuadrupole_XX*gqxx2 + atomI.labFrameQuadrupole_YY*gqyy2 + atomI.labFrameQuadrupole_ZZ*gqzz2 + 2.0f*(atomI.labFrameQuadrupole_XY*gqxy2 + atomI.labFrameQuadrupole_XZ*gqxz2 + atomI.labFrameQuadrupole_YZ*gqxy4)) +
              atomJ.inducedDipole[1]*(atomI.labFrameQuadrupole_XX*gqxx3 + atomI.labFrameQuadrupole_YY*gqyy3 + atomI.labFrameQuadrupole_ZZ*gqzz3 + 2.0f*(atomI.labFrameQuadrupole_XY*gqxy3 + atomI.labFrameQuadrupole_XZ*gqxy4 + atomI.labFrameQuadrupole_YZ*gqyz3)) +
              atomJ.inducedDipole[2]*(atomI.labFrameQuadrupole_XX*gqxx4 + atomI.labFrameQuadrupole_YY*gqyy4 + atomI.labFrameQuadrupole_ZZ*gqzz4 + 2.0f*(atomI.labFrameQuadrupole_XY*gqxy4 + atomI.labFrameQuadrupole_XZ*gqxz4 + atomI.labFrameQuadrupole_YZ*gqyz4));

    energy -= atomI.inducedDipole[0]*(atomJ.labFrameQuadrupole_XX*gqxx2 + atomJ.labFrameQuadrupole_YY*gqyy2 + atomJ.labFrameQuadrupole_ZZ*gqzz2 + 2.0f*(atomJ.labFrameQuadrupole_XY*gqxy2 + atomJ.labFrameQuadrupole_XZ*gqxz2 + atomJ.labFrameQuadrupole_YZ*gqxy4)) +
              atomI.inducedDipole[1]*(atomJ.labFrameQuadrupole_XX*gqxx3 + atomJ.labFrameQuadrupole_YY*gqyy3 + atomJ.labFrameQuadrupole_ZZ*gqzz3 + 2.0f*(atomJ.labFrameQuadrupole_XY*gqxy3 + atomJ.labFrameQuadrupole_XZ*gqxy4 + atomJ.labFrameQuadrupole_YZ*gqyz3)) +
              atomI.inducedDipole[2]*(atomJ.labFrameQuadrupole_XX*gqxx4 + atomJ.labFrameQuadrupole_YY*gqyy4 + atomJ.labFrameQuadrupole_ZZ*gqzz4 + 2.0f*(atomJ.labFrameQuadrupole_XY*gqxy4 + atomJ.labFrameQuadrupole_XZ*gqxz4 + atomJ.labFrameQuadrupole_YZ*gqyz4));

#endif
#endif

    // Born derivs of the unweighted reaction potential gradient tensor

#if defined B1
    dsumdrB1        += b01*atomI.q*(atomJ.labFrameDipole[0]*xr + atomJ.labFrameDipole[1]*yr + atomJ.labFrameDipole[2]*zr);
    dsumdrB1        -= b01*atomJ.q*(atomI.labFrameDipole[0]*xr + atomI.labFrameDipole[1]*yr + atomI.labFrameDipole[2]*zr);
#endif
#if defined B2
    dsumdrB2        += b01*atomI.q*(sxk*xr+ syk*yr + szk*zr);
    dsumdrB2        -= b01*atomJ.q*(sxi*xr+ syi*yr + szi*zr);
#endif

#if defined B1 || defined B2
    float gux22      = b10  + xr2*b11;
    float gux23      = xr*yr*b11;
    float gux24      = xr*zr*b11;
    float guy22      = gux23;
    float guy23      = b10  + yr2*b11;
    float guy24      = yr*zr*b11;
    float guz22      = gux24;
    float guz23      = guy24;
    float guz24      = b10  + zr2*b11;
#if defined B1
    dsumdrB1        -= 2.0f*( atomI.labFrameDipole[0]*(atomJ.labFrameDipole[0]*gux22 + atomJ.labFrameDipole[1]*guy22 + atomJ.labFrameDipole[2]*guz22) +
                              atomI.labFrameDipole[1]*(atomJ.labFrameDipole[0]*gux23 + atomJ.labFrameDipole[1]*guy23 + atomJ.labFrameDipole[2]*guz23) +
                              atomI.labFrameDipole[2]*(atomJ.labFrameDipole[0]*gux24 + atomJ.labFrameDipole[1]*guy24 + atomJ.labFrameDipole[2]*guz24));
#endif
#if defined B2
    dsumdrB2       -= 2.0f*(atomI.labFrameDipole[0]*(sxk*gux22 + syk*guy22 + szk*guz22) +
                            atomI.labFrameDipole[1]*(sxk*gux23 + syk*guy23 + szk*guz23) +
                            atomI.labFrameDipole[2]*(sxk*gux24 + syk*guy24 + szk*guz24) +
                            atomJ.labFrameDipole[0]*(sxi*gux22 + syi*guy22 + szi*guz22) +
                            atomJ.labFrameDipole[1]*(sxi*gux23 + syi*guy23 + szi*guz23) +
                            atomJ.labFrameDipole[2]*(sxi*gux24 + syi*guy24 + szi*guz24) );

    if( cAmoebaSim.polarizationType == 0 ){
        dsumdrB2 -=  2.0f*(    atomI.inducedDipole[0]*(atomJ.inducedDipoleP[0]*gux22 + atomJ.inducedDipoleP[1]*gux23 + atomJ.inducedDipoleP[2]*gux24)
                             + atomI.inducedDipole[1]*(atomJ.inducedDipoleP[0]*guy22 + atomJ.inducedDipoleP[1]*guy23 + atomJ.inducedDipoleP[2]*guy24)
                             + atomI.inducedDipole[2]*(atomJ.inducedDipoleP[0]*guz22 + atomJ.inducedDipoleP[1]*guz23 + atomJ.inducedDipoleP[2]*guz24)
                             + atomJ.inducedDipole[0]*(atomI.inducedDipoleP[0]*gux22 + atomI.inducedDipoleP[1]*gux23 + atomI.inducedDipoleP[2]*gux24)
                             + atomJ.inducedDipole[1]*(atomI.inducedDipoleP[0]*guy22 + atomI.inducedDipoleP[1]*guy23 + atomI.inducedDipoleP[2]*guy24)
                             + atomJ.inducedDipole[2]*(atomI.inducedDipoleP[0]*guz22 + atomI.inducedDipoleP[1]*guz23 + atomI.inducedDipoleP[2]*guz24));

    }
#endif
    float gqxx22     = xr*(2.0f*b20 + xr2*b21);
    float gqxx23     = yr*xr2*b21;
    float gqxx24     = zr*xr2*b21;
    float gqyy22     = xr*yr2*b21;
    float gqyy23     = yr*(2.0f*b20 + yr2*b21);
    float gqyy24     = zr*yr2*b21;
    float gqzz22     = xr*zr2*b21;
    float gqzz23     = yr*zr2*b21;
    float gqzz24     = zr*(2.0f*b20  + zr2*b21);
    float gqxy22     = yr*(b20 + xr2*b21);
    float gqxy23     = xr*(b20 + yr2*b21);
    float gqxy24     = zr*xr*yr*b21;
    float gqxz22     = zr*(b20 + xr2*b21);
    float gqxz23     = gqxy24;
    float gqxz24     = xr*(b20 + zr2*b21);
    float gqyz22     = gqxy24;
    float gqyz23     = zr*(b20 + yr2*b21);
    float gqyz24     = yr*(b20 + zr2*b21);
#if defined B1
    dsumdrB1 += atomJ.labFrameDipole[0]*(atomI.labFrameQuadrupole_XX*gqxx22 + atomI.labFrameQuadrupole_YY*gqyy22 + atomI.labFrameQuadrupole_ZZ*gqzz22 + 2.0f*(atomI.labFrameQuadrupole_XY*gqxy22 + atomI.labFrameQuadrupole_XZ*gqxz22 + atomI.labFrameQuadrupole_YZ*gqyz22)) +
                atomJ.labFrameDipole[1]*(atomI.labFrameQuadrupole_XX*gqxx23 + atomI.labFrameQuadrupole_YY*gqyy23 + atomI.labFrameQuadrupole_ZZ*gqzz23 + 2.0f*(atomI.labFrameQuadrupole_XY*gqxy23 + atomI.labFrameQuadrupole_XZ*gqxz23 + atomI.labFrameQuadrupole_YZ*gqyz23)) +
                atomJ.labFrameDipole[2]*(atomI.labFrameQuadrupole_XX*gqxx24 + atomI.labFrameQuadrupole_YY*gqyy24 + atomI.labFrameQuadrupole_ZZ*gqzz24 + 2.0f*(atomI.labFrameQuadrupole_XY*gqxy24 + atomI.labFrameQuadrupole_XZ*gqxz24 + atomI.labFrameQuadrupole_YZ*gqyz24));
    dsumdrB1 -= atomI.labFrameDipole[0]*(atomJ.labFrameQuadrupole_XX*gqxx22 + atomJ.labFrameQuadrupole_YY*gqyy22 + atomJ.labFrameQuadrupole_ZZ*gqzz22 + 2.0f*(atomJ.labFrameQuadrupole_XY*gqxy22 + atomJ.labFrameQuadrupole_XZ*gqxz22 + atomJ.labFrameQuadrupole_YZ*gqyz22)) +
                atomI.labFrameDipole[1]*(atomJ.labFrameQuadrupole_XX*gqxx23 + atomJ.labFrameQuadrupole_YY*gqyy23 + atomJ.labFrameQuadrupole_ZZ*gqzz23 + 2.0f*(atomJ.labFrameQuadrupole_XY*gqxy23 + atomJ.labFrameQuadrupole_XZ*gqxz23 + atomJ.labFrameQuadrupole_YZ*gqyz23)) +
                atomI.labFrameDipole[2]*(atomJ.labFrameQuadrupole_XX*gqxx24 + atomJ.labFrameQuadrupole_YY*gqyy24 + atomJ.labFrameQuadrupole_ZZ*gqzz24 + 2.0f*(atomJ.labFrameQuadrupole_XY*gqxy24 + atomJ.labFrameQuadrupole_XZ*gqxz24 + atomJ.labFrameQuadrupole_YZ*gqyz24));
#endif
#if defined B2

    dsumdrB2 += sxk*(atomI.labFrameQuadrupole_XX*gqxx22 + atomI.labFrameQuadrupole_YY*gqyy22 + atomI.labFrameQuadrupole_ZZ*gqzz22 + 2.0f*(atomI.labFrameQuadrupole_XY*gqxy22 + atomI.labFrameQuadrupole_XZ*gqxz22 + atomI.labFrameQuadrupole_YZ*gqyz22)) +
                syk*(atomI.labFrameQuadrupole_XX*gqxx23 + atomI.labFrameQuadrupole_YY*gqyy23 + atomI.labFrameQuadrupole_ZZ*gqzz23 + 2.0f*(atomI.labFrameQuadrupole_XY*gqxy23 + atomI.labFrameQuadrupole_XZ*gqxz23 + atomI.labFrameQuadrupole_YZ*gqyz23)) +
                szk*(atomI.labFrameQuadrupole_XX*gqxx24 + atomI.labFrameQuadrupole_YY*gqyy24 + atomI.labFrameQuadrupole_ZZ*gqzz24 + 2.0f*(atomI.labFrameQuadrupole_XY*gqxy24 + atomI.labFrameQuadrupole_XZ*gqxz24 + atomI.labFrameQuadrupole_YZ*gqyz24));

    dsumdrB2 -= sxi*(atomJ.labFrameQuadrupole_XX*gqxx22 + atomJ.labFrameQuadrupole_YY*gqyy22 + atomJ.labFrameQuadrupole_ZZ*gqzz22 + 2.0f*(atomJ.labFrameQuadrupole_XY*gqxy22 + atomJ.labFrameQuadrupole_XZ*gqxz22 + atomJ.labFrameQuadrupole_YZ*gqyz22)) +
                syi*(atomJ.labFrameQuadrupole_XX*gqxx23 + atomJ.labFrameQuadrupole_YY*gqyy23 + atomJ.labFrameQuadrupole_ZZ*gqzz23 + 2.0f*(atomJ.labFrameQuadrupole_XY*gqxy23 + atomJ.labFrameQuadrupole_XZ*gqxz23 + atomJ.labFrameQuadrupole_YZ*gqyz23)) +
                szi*(atomJ.labFrameQuadrupole_XX*gqxx24 + atomJ.labFrameQuadrupole_YY*gqyy24 + atomJ.labFrameQuadrupole_ZZ*gqzz24 + 2.0f*(atomJ.labFrameQuadrupole_XY*gqxy24 + atomJ.labFrameQuadrupole_XZ*gqxz24 + atomJ.labFrameQuadrupole_YZ*gqyz24));

#endif
#endif

    // unweighted 2nd reaction potential gradient tensor;

#if defined F1 || defined F2 || defined T1
    float gc5        = a01  + xr2*a02;
    float gc6        = xr*yr*a02;
    float gc7        = xr*zr*a02;
    float gc8        = a01  + yr2*a02;
    float gc9        = yr*zr*a02;
    float gc10       = a01  + zr2*a02;
#if defined F1
    energy          += atomI.q*(atomJ.labFrameQuadrupole_XX*gc5 + atomJ.labFrameQuadrupole_YY*gc8 + atomJ.labFrameQuadrupole_ZZ*gc10 + 2.0f*(atomJ.labFrameQuadrupole_XY*gc6 + atomJ.labFrameQuadrupole_XZ*gc7 + atomJ.labFrameQuadrupole_YZ*gc9));
    energy          += atomJ.q*(atomI.labFrameQuadrupole_XX*gc5 + atomI.labFrameQuadrupole_YY*gc8 + atomI.labFrameQuadrupole_ZZ*gc10 + 2.0f*(atomI.labFrameQuadrupole_XY*gc6 + atomI.labFrameQuadrupole_XZ*gc7 + atomI.labFrameQuadrupole_YZ*gc9));

    dedx            += atomI.q*(atomJ.labFrameDipole[0]*gc5 + atomJ.labFrameDipole[1]*gc6 + atomJ.labFrameDipole[2]*gc7);
    dedx            -= atomJ.q*(atomI.labFrameDipole[0]*gc5 + atomI.labFrameDipole[1]*gc6 + atomI.labFrameDipole[2]*gc7);

    dedy            += atomI.q*(atomJ.labFrameDipole[0]*gc6 + atomJ.labFrameDipole[1]*gc8 + atomJ.labFrameDipole[2]*gc9);
    dedy            -= atomJ.q*(atomI.labFrameDipole[0]*gc6 + atomI.labFrameDipole[1]*gc8 + atomI.labFrameDipole[2]*gc9);

    dedz            += atomI.q*(atomJ.labFrameDipole[0]*gc7 + atomJ.labFrameDipole[1]*gc9 + atomJ.labFrameDipole[2]*gc10);
    dedz            -= atomJ.q*(atomI.labFrameDipole[0]*gc7 + atomI.labFrameDipole[1]*gc9 + atomI.labFrameDipole[2]*gc10);
#endif

#if defined F2
    dpdx            += atomI.q*(sxk*gc5 + syk*gc6 + szk*gc7);
    dpdx            -= atomJ.q*(sxi*gc5 + syi*gc6 + szi*gc7);
    dpdy            += atomI.q*(sxk*gc6 + syk*gc8 + szk*gc9);
    dpdy            -= atomJ.q*(sxi*gc6 + syi*gc8 + szi*gc9);
    dpdz            += atomI.q*(sxk*gc7 + syk*gc9 + szk*gc10);
    dpdz            -= atomJ.q*(sxi*gc7 + syi*gc9 + szi*gc10);
#endif

#endif

#if defined F1 || defined F2 || defined T1 || defined T2
    float gux5       = xr*(3.0f*a11 + xr2*a12);
    float gux6       = yr*(a11 + xr2*a12);
    float gux7       = zr*(a11 + xr2*a12);
    float gux8       = xr*(a11 + yr2*a12);
    float gux9       = zr*xr*yr*a12;
    float gux10      = xr*(a11 + zr2*a12);
    float guy5       = yr*(a11 + xr2*a12);
    float guy6       = xr*(a11 + yr2*a12);
    float guy7       = gux9;
    float guy8       = yr*(3.0f*a11 + yr2*a12);
    float guy9       = zr*(a11 + yr2*a12);
    float guy10      = yr*(a11 + zr2*a12);
    float guz5       = zr*(a11 + xr2*a12);
    float guz6       = gux9;
    float guz7       = xr*(a11 + zr2*a12);
    float guz8       = zr*(a11 + yr2*a12);
    float guz9       = yr*(a11 + zr2*a12);
    float guz10      = zr*(3.0f*a11 + zr2*a12);
#if defined F1
    energy          -= atomI.labFrameDipole[0]*(atomJ.labFrameQuadrupole_XX*gux5 + atomJ.labFrameQuadrupole_YY*gux8 + atomJ.labFrameQuadrupole_ZZ*gux10  + 2.0f*(atomJ.labFrameQuadrupole_XY*gux6 + atomJ.labFrameQuadrupole_XZ*gux7 + atomJ.labFrameQuadrupole_YZ*gux9)) +
                       atomI.labFrameDipole[1]*(atomJ.labFrameQuadrupole_XX*guy5 + atomJ.labFrameQuadrupole_YY*guy8 + atomJ.labFrameQuadrupole_ZZ*guy10  + 2.0f*(atomJ.labFrameQuadrupole_XY*guy6 + atomJ.labFrameQuadrupole_XZ*guy7 + atomJ.labFrameQuadrupole_YZ*guy9)) +
                       atomI.labFrameDipole[2]*(atomJ.labFrameQuadrupole_XX*guz5 + atomJ.labFrameQuadrupole_YY*guz8 + atomJ.labFrameQuadrupole_ZZ*guz10  + 2.0f*(atomJ.labFrameQuadrupole_XY*guz6 + atomJ.labFrameQuadrupole_XZ*guz7 + atomJ.labFrameQuadrupole_YZ*guz9));

    energy          += atomJ.labFrameDipole[0]*(atomI.labFrameQuadrupole_XX*gux5 + atomI.labFrameQuadrupole_YY*gux8 + atomI.labFrameQuadrupole_ZZ*gux10  + 2.0f*(atomI.labFrameQuadrupole_XY*gux6 + atomI.labFrameQuadrupole_XZ*gux7 + atomI.labFrameQuadrupole_YZ*gux9)) +
                       atomJ.labFrameDipole[1]*(atomI.labFrameQuadrupole_XX*guy5 + atomI.labFrameQuadrupole_YY*guy8 + atomI.labFrameQuadrupole_ZZ*guy10  + 2.0f*(atomI.labFrameQuadrupole_XY*guy6 + atomI.labFrameQuadrupole_XZ*guy7 + atomI.labFrameQuadrupole_YZ*guy9)) +
                       atomJ.labFrameDipole[2]*(atomI.labFrameQuadrupole_XX*guz5 + atomI.labFrameQuadrupole_YY*guz8 + atomI.labFrameQuadrupole_ZZ*guz10  + 2.0f*(atomI.labFrameQuadrupole_XY*guz6 + atomI.labFrameQuadrupole_XZ*guz7 + atomI.labFrameQuadrupole_YZ*guz9));

    dedx            -= 2.0f*( atomI.labFrameDipole[0]*(atomJ.labFrameDipole[0]*gux5 + atomJ.labFrameDipole[1]*guy5 + atomJ.labFrameDipole[2]*guz5) +
                              atomI.labFrameDipole[1]*(atomJ.labFrameDipole[0]*gux6 + atomJ.labFrameDipole[1]*guy6 + atomJ.labFrameDipole[2]*guz6) +
                              atomI.labFrameDipole[2]*(atomJ.labFrameDipole[0]*gux7 + atomJ.labFrameDipole[1]*guy7 + atomJ.labFrameDipole[2]*guz7));

    dedy            -= 2.0f*( atomI.labFrameDipole[0]*(atomJ.labFrameDipole[0]*gux6 + atomJ.labFrameDipole[1]*guy6 + atomJ.labFrameDipole[2]*guz6) +
                              atomI.labFrameDipole[1]*(atomJ.labFrameDipole[0]*gux8 + atomJ.labFrameDipole[1]*guy8 + atomJ.labFrameDipole[2]*guz8) +
                              atomI.labFrameDipole[2]*(atomJ.labFrameDipole[0]*gux9 + atomJ.labFrameDipole[1]*guy9 + atomJ.labFrameDipole[2]*guz9));

    dedz            -= 2.0f*( atomI.labFrameDipole[0]*(atomJ.labFrameDipole[0]*gux7 + atomJ.labFrameDipole[1]*guy7 + atomJ.labFrameDipole[2]*guz7) +
                              atomI.labFrameDipole[1]*(atomJ.labFrameDipole[0]*gux9 + atomJ.labFrameDipole[1]*guy9 + atomJ.labFrameDipole[2]*guz9) +
                              atomI.labFrameDipole[2]*(atomJ.labFrameDipole[0]*gux10 + atomJ.labFrameDipole[1]*guy10 + atomJ.labFrameDipole[2]*guz10));

#endif

#if defined F2
    energy -= atomI.inducedDipole[0]*(atomJ.labFrameQuadrupole_XX*gux5 + atomJ.labFrameQuadrupole_YY*gux8 + atomJ.labFrameQuadrupole_ZZ*gux10 + 2.0f*(atomJ.labFrameQuadrupole_XY*gux6 + atomJ.labFrameQuadrupole_XZ*gux7 + atomJ.labFrameQuadrupole_YZ*gux9)) +
              atomI.inducedDipole[1]*(atomJ.labFrameQuadrupole_XX*guy5 + atomJ.labFrameQuadrupole_YY*guy8 + atomJ.labFrameQuadrupole_ZZ*guy10 + 2.0f*(atomJ.labFrameQuadrupole_XY*guy6 + atomJ.labFrameQuadrupole_XZ*guy7 + atomJ.labFrameQuadrupole_YZ*guy9)) +
              atomI.inducedDipole[2]*(atomJ.labFrameQuadrupole_XX*guz5 + atomJ.labFrameQuadrupole_YY*guz8 + atomJ.labFrameQuadrupole_ZZ*guz10 + 2.0f*(atomJ.labFrameQuadrupole_XY*guz6 + atomJ.labFrameQuadrupole_XZ*guz7 + atomJ.labFrameQuadrupole_YZ*guz9));

    energy += atomJ.inducedDipole[0]*(atomI.labFrameQuadrupole_XX*gux5 + atomI.labFrameQuadrupole_YY*gux8 + atomI.labFrameQuadrupole_ZZ*gux10 + 2.0f*(atomI.labFrameQuadrupole_XY*gux6 + atomI.labFrameQuadrupole_XZ*gux7 + atomI.labFrameQuadrupole_YZ*gux9)) +
              atomJ.inducedDipole[1]*(atomI.labFrameQuadrupole_XX*guy5 + atomI.labFrameQuadrupole_YY*guy8 + atomI.labFrameQuadrupole_ZZ*guy10 + 2.0f*(atomI.labFrameQuadrupole_XY*guy6 + atomI.labFrameQuadrupole_XZ*guy7 + atomI.labFrameQuadrupole_YZ*guy9)) +
              atomJ.inducedDipole[2]*(atomI.labFrameQuadrupole_XX*guz5 + atomI.labFrameQuadrupole_YY*guz8 + atomI.labFrameQuadrupole_ZZ*guz10 + 2.0f*(atomI.labFrameQuadrupole_XY*guz6 + atomI.labFrameQuadrupole_XZ*guz7 + atomI.labFrameQuadrupole_YZ*guz9));

    dpdx   -= 2.0f*(atomI.labFrameDipole[0]*(sxk*gux5 + syk*guy5 + szk*guz5) + atomI.labFrameDipole[1]*(sxk*gux6 + syk*guy6 + szk*guz6) + atomI.labFrameDipole[2]*(sxk*gux7 + syk*guy7 + szk*guz7) +
                    atomJ.labFrameDipole[0]*(sxi*gux5 + syi*guy5 + szi*guz5) + atomJ.labFrameDipole[1]*(sxi*gux6 + syi*guy6 + szi*guz6) + atomJ.labFrameDipole[2]*(sxi*gux7 + syi*guy7 + szi*guz7));

    dpdy   -= 2.0f*(atomI.labFrameDipole[0]*(sxk*gux6 + syk*guy6 + szk*guz6) + atomI.labFrameDipole[1]*(sxk*gux8 + syk*guy8 + szk*guz8) + atomI.labFrameDipole[2]*(sxk*gux9 + syk*guy9 + szk*guz9) +
                    atomJ.labFrameDipole[0]*(sxi*gux6 + syi*guy6 + szi*guz6) + atomJ.labFrameDipole[1]*(sxi*gux8 + syi*guy8 + szi*guz8) + atomJ.labFrameDipole[2]*(sxi*gux9 + syi*guy9 + szi*guz9));

    dpdz   -= 2.0f*(atomI.labFrameDipole[0]*(sxk*gux7 + syk*guy7 + szk*guz7) + atomI.labFrameDipole[1]*(sxk*gux9 + syk*guy9 + szk*guz9) + atomI.labFrameDipole[2]*(sxk*gux10 + syk*guy10 + szk*guz10) +
                    atomJ.labFrameDipole[0]*(sxi*gux7 + syi*guy7 + szi*guz7) + atomJ.labFrameDipole[1]*(sxi*gux9 + syi*guy9 + szi*guz9) + atomJ.labFrameDipole[2]*(sxi*gux10 + syi*guy10 + szi*guz10) );

    if( cAmoebaSim.polarizationType == 0 ){

        dpdx -=         2.0f*(atomI.inducedDipole[0]*(atomJ.inducedDipoleP[0]*gux5 + atomJ.inducedDipoleP[1]*gux6 + atomJ.inducedDipoleP[2]*gux7)
                            + atomI.inducedDipole[1]*(atomJ.inducedDipoleP[0]*guy5 + atomJ.inducedDipoleP[1]*guy6 + atomJ.inducedDipoleP[2]*guy7)
                            + atomI.inducedDipole[2]*(atomJ.inducedDipoleP[0]*guz5 + atomJ.inducedDipoleP[1]*guz6 + atomJ.inducedDipoleP[2]*guz7)
                            + atomJ.inducedDipole[0]*(atomI.inducedDipoleP[0]*gux5 + atomI.inducedDipoleP[1]*gux6 + atomI.inducedDipoleP[2]*gux7)
                            + atomJ.inducedDipole[1]*(atomI.inducedDipoleP[0]*guy5 + atomI.inducedDipoleP[1]*guy6 + atomI.inducedDipoleP[2]*guy7)
                            + atomJ.inducedDipole[2]*(atomI.inducedDipoleP[0]*guz5 + atomI.inducedDipoleP[1]*guz6 + atomI.inducedDipoleP[2]*guz7));

        dpdy -=         2.0f*(atomI.inducedDipole[0]*(atomJ.inducedDipoleP[0]*gux6 + atomJ.inducedDipoleP[1]*gux8 + atomJ.inducedDipoleP[2]*gux9)
                            + atomI.inducedDipole[1]*(atomJ.inducedDipoleP[0]*guy6 + atomJ.inducedDipoleP[1]*guy8 + atomJ.inducedDipoleP[2]*guy9)
                            + atomI.inducedDipole[2]*(atomJ.inducedDipoleP[0]*guz6 + atomJ.inducedDipoleP[1]*guz8 + atomJ.inducedDipoleP[2]*guz9)
                            + atomJ.inducedDipole[0]*(atomI.inducedDipoleP[0]*gux6 + atomI.inducedDipoleP[1]*gux8 + atomI.inducedDipoleP[2]*gux9)
                            + atomJ.inducedDipole[1]*(atomI.inducedDipoleP[0]*guy6 + atomI.inducedDipoleP[1]*guy8 + atomI.inducedDipoleP[2]*guy9)
                            + atomJ.inducedDipole[2]*(atomI.inducedDipoleP[0]*guz6 + atomI.inducedDipoleP[1]*guz8 + atomI.inducedDipoleP[2]*guz9));

        dpdz -=         2.0f*(atomI.inducedDipole[0]*(atomJ.inducedDipoleP[0]*gux7 + atomJ.inducedDipoleP[1]*gux9 + atomJ.inducedDipoleP[2]*gux10)
                            + atomI.inducedDipole[1]*(atomJ.inducedDipoleP[0]*guy7 + atomJ.inducedDipoleP[1]*guy9 + atomJ.inducedDipoleP[2]*guy10)
                            + atomI.inducedDipole[2]*(atomJ.inducedDipoleP[0]*guz7 + atomJ.inducedDipoleP[1]*guz9 + atomJ.inducedDipoleP[2]*guz10)
                            + atomJ.inducedDipole[0]*(atomI.inducedDipoleP[0]*gux7 + atomI.inducedDipoleP[1]*gux9 + atomI.inducedDipoleP[2]*gux10)
                            + atomJ.inducedDipole[1]*(atomI.inducedDipoleP[0]*guy7 + atomI.inducedDipoleP[1]*guy9 + atomI.inducedDipoleP[2]*guy10)
                            + atomJ.inducedDipole[2]*(atomI.inducedDipoleP[0]*guz7 + atomI.inducedDipoleP[1]*guz9 + atomI.inducedDipoleP[2]*guz10));
    }
#endif
#endif

#if defined F1 || defined F2 || defined T1
    float gqxx5      = 2.0f*a20  + xr2*(5.0f*a21 + xr2*a22);
    float gqxx6      = yr*xr*(2.0f*a21 + xr2*a22);
    float gqxx7      = zr*xr*(2.0f*a21 + xr2*a22);
    float gqxx8      = xr2*(a21 + yr2*a22);
    float gqxx9      = zr*yr*xr2*a22;
    float gqxx10     = xr2*(a21 + zr2*a22);
    float gqyy5      = yr2*(a21 + xr2*a22);
    float gqyy6      = xr*yr*(2.0f*a21 + yr2*a22);
    float gqyy7      = xr*zr*yr2*a22;
    float gqyy8      = 2.0f*a20  + yr2*(5.0f*a21 + yr2*a22);
    float gqyy9      = yr*zr*(2.0f*a21 + yr2*a22);
    float gqyy10     = yr2*(a21 + zr2*a22);
    float gqzz5      = zr2*(a21 + xr2*a22);
    float gqzz6      = xr*yr*zr2*a22;
    float gqzz7      = xr*zr*(2.0f*a21 + zr2*a22);
    float gqzz8      = zr2*(a21 + yr2*a22);
    float gqzz9      = yr*zr*(2.0f*a21 + zr2*a22);
    float gqzz10     = 2.0f*a20  + zr2*(5.0f*a21 + zr2*a22);
    float gqxy5      = xr*yr*(3.0f*a21 + xr2*a22);
    float gqxy6      = a20  + (xr2 + yr2)*a21  + xr2*yr2*a22;
    float gqxy7      = zr*yr*(a21 + xr2*a22);
    float gqxy8      = xr*yr*(3.0f*a21 + yr2*a22);
    float gqxy9      = zr*xr*(a21 + yr2*a22);
    float gqxy10     = xr*yr*(a21 + zr2*a22);
    float gqxz5      = xr*zr*(3.0f*a21 + xr2*a22);
    float gqxz6      = yr*zr*(a21 + xr2*a22);
    float gqxz7      = a20  + (xr2 + zr2)*a21  + xr2*zr2*a22;
    float gqxz8      = xr*zr*(a21 + yr2*a22);
    float gqxz9      = xr*yr*(a21 + zr2*a22);
    float gqxz10     = xr*zr*(3.0f*a21 + zr2*a22);
    float gqyz5      = zr*yr*(a21 + xr2*a22);
    float gqyz6      = xr*zr*(a21 + yr2*a22);
    float gqyz7      = xr*yr*(a21 + zr2*a22);
    float gqyz8      = yr*zr*(3.0f*a21 + yr2*a22);
    float gqyz9      = a20  + (yr2 + zr2)*a21  + yr2*zr2*a22;
    float gqyz10     = yr*zr*(3.0f*a21 + zr2*a22);
#if defined F1
    energy +=   atomI.labFrameQuadrupole_XX*(atomJ.labFrameQuadrupole_XX*gqxx5 + atomJ.labFrameQuadrupole_YY*gqxx8 + atomJ.labFrameQuadrupole_ZZ*gqxx10 + 2.0f*(atomJ.labFrameQuadrupole_XY*gqxx6 + atomJ.labFrameQuadrupole_XZ*gqxx7 + atomJ.labFrameQuadrupole_YZ*gqxx9))
              + atomI.labFrameQuadrupole_YY*(atomJ.labFrameQuadrupole_XX*gqyy5 + atomJ.labFrameQuadrupole_YY*gqyy8 + atomJ.labFrameQuadrupole_ZZ*gqyy10 + 2.0f*(atomJ.labFrameQuadrupole_XY*gqyy6 + atomJ.labFrameQuadrupole_XZ*gqyy7 + atomJ.labFrameQuadrupole_YZ*gqyy9))
              + atomI.labFrameQuadrupole_ZZ*(atomJ.labFrameQuadrupole_XX*gqzz5 + atomJ.labFrameQuadrupole_YY*gqzz8 + atomJ.labFrameQuadrupole_ZZ*gqzz10 + 2.0f*(atomJ.labFrameQuadrupole_XY*gqzz6 + atomJ.labFrameQuadrupole_XZ*gqzz7 + atomJ.labFrameQuadrupole_YZ*gqzz9))
              + 2.0f*(atomI.labFrameQuadrupole_XY*(atomJ.labFrameQuadrupole_XX*gqxy5 + atomJ.labFrameQuadrupole_YY*gqxy8 + atomJ.labFrameQuadrupole_ZZ*gqxy10
              + 2.0f*(atomJ.labFrameQuadrupole_XY*gqxy6 + atomJ.labFrameQuadrupole_XZ*gqxy7 + atomJ.labFrameQuadrupole_YZ*gqxy9))
              + atomI.labFrameQuadrupole_XZ*(atomJ.labFrameQuadrupole_XX*gqxz5 + atomJ.labFrameQuadrupole_YY*gqxz8 + atomJ.labFrameQuadrupole_ZZ*gqxz10
              + 2.0f*(atomJ.labFrameQuadrupole_XY*gqxz6 + atomJ.labFrameQuadrupole_XZ*gqxz7 + atomJ.labFrameQuadrupole_YZ*gqxz9))
              + atomI.labFrameQuadrupole_YZ*(atomJ.labFrameQuadrupole_XX*gqyz5 + atomJ.labFrameQuadrupole_YY*gqyz8 + atomJ.labFrameQuadrupole_ZZ*gqyz10
              + 2.0f*(atomJ.labFrameQuadrupole_XY*gqyz6 + atomJ.labFrameQuadrupole_XZ*gqyz7 + atomJ.labFrameQuadrupole_YZ*gqyz9)));

    energy +=  atomI.labFrameQuadrupole_XX*(atomJ.labFrameQuadrupole_XX*gqxx5 + atomJ.labFrameQuadrupole_YY*gqyy5 + atomJ.labFrameQuadrupole_ZZ*gqzz5 + 2.0f*(atomJ.labFrameQuadrupole_XY*gqxy5 + atomJ.labFrameQuadrupole_XZ*gqxz5 + atomJ.labFrameQuadrupole_YZ*gqyz5))
             + atomI.labFrameQuadrupole_YY*(atomJ.labFrameQuadrupole_XX*gqxx8 + atomJ.labFrameQuadrupole_YY*gqyy8 + atomJ.labFrameQuadrupole_ZZ*gqzz8
             + 2.0f*(atomJ.labFrameQuadrupole_XY*gqxy8 + atomJ.labFrameQuadrupole_XZ*gqxz8 + atomJ.labFrameQuadrupole_YZ*gqyz8))
             + atomI.labFrameQuadrupole_ZZ*(atomJ.labFrameQuadrupole_XX*gqxx10 + atomJ.labFrameQuadrupole_YY*gqyy10 + atomJ.labFrameQuadrupole_ZZ*gqzz10
             + 2.0f*(atomJ.labFrameQuadrupole_XY*gqxy10 + atomJ.labFrameQuadrupole_XZ*gqxz10 + atomJ.labFrameQuadrupole_YZ*gqyz10))
             + 2.0f*(atomI.labFrameQuadrupole_XY*(atomJ.labFrameQuadrupole_XX*gqxx6 + atomJ.labFrameQuadrupole_YY*gqyy6 + atomJ.labFrameQuadrupole_ZZ*gqzz6
             + 2.0f*(atomJ.labFrameQuadrupole_XY*gqxy6 + atomJ.labFrameQuadrupole_XZ*gqxz6 + atomJ.labFrameQuadrupole_YZ*gqyz6))
             + atomI.labFrameQuadrupole_XZ*(atomJ.labFrameQuadrupole_XX*gqxx7 + atomJ.labFrameQuadrupole_YY*gqyy7 + atomJ.labFrameQuadrupole_ZZ*gqzz7
             + 2.0f*(atomJ.labFrameQuadrupole_XY*gqxy7 + atomJ.labFrameQuadrupole_XZ*gqxz7 + atomJ.labFrameQuadrupole_YZ*gqyz7))
             + atomI.labFrameQuadrupole_YZ*(atomJ.labFrameQuadrupole_XX*gqxx9 + atomJ.labFrameQuadrupole_YY*gqyy9 + atomJ.labFrameQuadrupole_ZZ*gqzz9
             + 2.0f*(atomJ.labFrameQuadrupole_XY*gqxy9 + atomJ.labFrameQuadrupole_XZ*gqxz9 + atomJ.labFrameQuadrupole_YZ*gqyz9)));

    dedx   += atomJ.labFrameDipole[0]*(atomI.labFrameQuadrupole_XX*gqxx5 + atomI.labFrameQuadrupole_YY*gqyy5 + atomI.labFrameQuadrupole_ZZ*gqzz5 + 2.0f*(atomI.labFrameQuadrupole_XY*gqxy5 + atomI.labFrameQuadrupole_XZ*gqxz5 + atomI.labFrameQuadrupole_YZ*gqyz5)) +
              atomJ.labFrameDipole[1]*(atomI.labFrameQuadrupole_XX*gqxx6 + atomI.labFrameQuadrupole_YY*gqyy6 + atomI.labFrameQuadrupole_ZZ*gqzz6 + 2.0f*(atomI.labFrameQuadrupole_XY*gqxy6 + atomI.labFrameQuadrupole_XZ*gqxz6 + atomI.labFrameQuadrupole_YZ*gqyz6)) +
              atomJ.labFrameDipole[2]*(atomI.labFrameQuadrupole_XX*gqxx7 + atomI.labFrameQuadrupole_YY*gqyy7 + atomI.labFrameQuadrupole_ZZ*gqzz7 + 2.0f*(atomI.labFrameQuadrupole_XY*gqxy7 + atomI.labFrameQuadrupole_XZ*gqxz7 + atomI.labFrameQuadrupole_YZ*gqyz7));

    dedx   -= atomI.labFrameDipole[0]*(atomJ.labFrameQuadrupole_XX*gqxx5 + atomJ.labFrameQuadrupole_YY*gqyy5 + atomJ.labFrameQuadrupole_ZZ*gqzz5 + 2.0f*(atomJ.labFrameQuadrupole_XY*gqxy5 + atomJ.labFrameQuadrupole_XZ*gqxz5 + atomJ.labFrameQuadrupole_YZ*gqyz5)) +
              atomI.labFrameDipole[1]*(atomJ.labFrameQuadrupole_XX*gqxx6 + atomJ.labFrameQuadrupole_YY*gqyy6 + atomJ.labFrameQuadrupole_ZZ*gqzz6 + 2.0f*(atomJ.labFrameQuadrupole_XY*gqxy6 + atomJ.labFrameQuadrupole_XZ*gqxz6 + atomJ.labFrameQuadrupole_YZ*gqyz6)) +
              atomI.labFrameDipole[2]*(atomJ.labFrameQuadrupole_XX*gqxx7 + atomJ.labFrameQuadrupole_YY*gqyy7 + atomJ.labFrameQuadrupole_ZZ*gqzz7 + 2.0f*(atomJ.labFrameQuadrupole_XY*gqxy7 + atomJ.labFrameQuadrupole_XZ*gqxz7 + atomJ.labFrameQuadrupole_YZ*gqyz7));

    dedy   += atomJ.labFrameDipole[0]*(atomI.labFrameQuadrupole_XX*gqxx6 + atomI.labFrameQuadrupole_YY*gqyy6 + atomI.labFrameQuadrupole_ZZ*gqzz6 + 2.0f*(atomI.labFrameQuadrupole_XY*gqxy6 + atomI.labFrameQuadrupole_XZ*gqxz6 + atomI.labFrameQuadrupole_YZ*gqyz6)) +
              atomJ.labFrameDipole[1]*(atomI.labFrameQuadrupole_XX*gqxx8 + atomI.labFrameQuadrupole_YY*gqyy8 + atomI.labFrameQuadrupole_ZZ*gqzz8 + 2.0f*(atomI.labFrameQuadrupole_XY*gqxy8 + atomI.labFrameQuadrupole_XZ*gqxz8 + atomI.labFrameQuadrupole_YZ*gqyz8)) +
              atomJ.labFrameDipole[2]*(atomI.labFrameQuadrupole_XX*gqxx9 + atomI.labFrameQuadrupole_YY*gqyy9 + atomI.labFrameQuadrupole_ZZ*gqzz9 + 2.0f*(atomI.labFrameQuadrupole_XY*gqxy9 + atomI.labFrameQuadrupole_XZ*gqxz9 + atomI.labFrameQuadrupole_YZ*gqyz9));

    dedy   -= atomI.labFrameDipole[0]*(atomJ.labFrameQuadrupole_XX*gqxx6 + atomJ.labFrameQuadrupole_YY*gqyy6 + atomJ.labFrameQuadrupole_ZZ*gqzz6 + 2.0f*(atomJ.labFrameQuadrupole_XY*gqxy6 + atomJ.labFrameQuadrupole_XZ*gqxz6 + atomJ.labFrameQuadrupole_YZ*gqyz6)) +
              atomI.labFrameDipole[1]*(atomJ.labFrameQuadrupole_XX*gqxx8 + atomJ.labFrameQuadrupole_YY*gqyy8 + atomJ.labFrameQuadrupole_ZZ*gqzz8 + 2.0f*(atomJ.labFrameQuadrupole_XY*gqxy8 + atomJ.labFrameQuadrupole_XZ*gqxz8 + atomJ.labFrameQuadrupole_YZ*gqyz8)) +
              atomI.labFrameDipole[2]*(atomJ.labFrameQuadrupole_XX*gqxx9 + atomJ.labFrameQuadrupole_YY*gqyy9 + atomJ.labFrameQuadrupole_ZZ*gqzz9 + 2.0f*(atomJ.labFrameQuadrupole_XY*gqxy9 + atomJ.labFrameQuadrupole_XZ*gqxz9 + atomJ.labFrameQuadrupole_YZ*gqyz9));

    dedz   += atomJ.labFrameDipole[0]*(atomI.labFrameQuadrupole_XX*gqxx7 + atomI.labFrameQuadrupole_YY*gqyy7 + atomI.labFrameQuadrupole_ZZ*gqzz7 + 2.0f*(atomI.labFrameQuadrupole_XY*gqxy7 + atomI.labFrameQuadrupole_XZ*gqxz7 + atomI.labFrameQuadrupole_YZ*gqyz7)) +
              atomJ.labFrameDipole[1]*(atomI.labFrameQuadrupole_XX*gqxx9 + atomI.labFrameQuadrupole_YY*gqyy9 + atomI.labFrameQuadrupole_ZZ*gqzz9 + 2.0f*(atomI.labFrameQuadrupole_XY*gqxy9 + atomI.labFrameQuadrupole_XZ*gqxz9 + atomI.labFrameQuadrupole_YZ*gqyz9)) +
              atomJ.labFrameDipole[2]*(atomI.labFrameQuadrupole_XX*gqxx10 + atomI.labFrameQuadrupole_YY*gqyy10 + atomI.labFrameQuadrupole_ZZ*gqzz10 + 2.0f*(atomI.labFrameQuadrupole_XY*gqxy10 + atomI.labFrameQuadrupole_XZ*gqxz10 + atomI.labFrameQuadrupole_YZ*gqyz10));

    dedz   -= atomI.labFrameDipole[0]*(atomJ.labFrameQuadrupole_XX*gqxx7 + atomJ.labFrameQuadrupole_YY*gqyy7 + atomJ.labFrameQuadrupole_ZZ*gqzz7 + 2.0f*(atomJ.labFrameQuadrupole_XY*gqxy7 + atomJ.labFrameQuadrupole_XZ*gqxz7 + atomJ.labFrameQuadrupole_YZ*gqyz7)) +
              atomI.labFrameDipole[1]*(atomJ.labFrameQuadrupole_XX*gqxx9 + atomJ.labFrameQuadrupole_YY*gqyy9 + atomJ.labFrameQuadrupole_ZZ*gqzz9 + 2.0f*(atomJ.labFrameQuadrupole_XY*gqxy9 + atomJ.labFrameQuadrupole_XZ*gqxz9 + atomJ.labFrameQuadrupole_YZ*gqyz9)) +
              atomI.labFrameDipole[2]*(atomJ.labFrameQuadrupole_XX*gqxx10 + atomJ.labFrameQuadrupole_YY*gqyy10 + atomJ.labFrameQuadrupole_ZZ*gqzz10 + 2.0f*(atomJ.labFrameQuadrupole_XY*gqxy10 + atomJ.labFrameQuadrupole_XZ*gqxz10 + atomJ.labFrameQuadrupole_YZ*gqyz10));

#endif
#if defined F2
    dpdx += sxk*(atomI.labFrameQuadrupole_XX*gqxx5 + atomI.labFrameQuadrupole_YY*gqyy5 + atomI.labFrameQuadrupole_ZZ*gqzz5 + 2.0f*(atomI.labFrameQuadrupole_XY*gqxy5 + atomI.labFrameQuadrupole_XZ*gqxz5 + atomI.labFrameQuadrupole_YZ*gqyz5)) +
            syk*(atomI.labFrameQuadrupole_XX*gqxx6 + atomI.labFrameQuadrupole_YY*gqyy6 + atomI.labFrameQuadrupole_ZZ*gqzz6 + 2.0f*(atomI.labFrameQuadrupole_XY*gqxy6 + atomI.labFrameQuadrupole_XZ*gqxz6 + atomI.labFrameQuadrupole_YZ*gqyz6)) +
            szk*(atomI.labFrameQuadrupole_XX*gqxx7 + atomI.labFrameQuadrupole_YY*gqyy7 + atomI.labFrameQuadrupole_ZZ*gqzz7 + 2.0f*(atomI.labFrameQuadrupole_XY*gqxy7 + atomI.labFrameQuadrupole_XZ*gqxz7 + atomI.labFrameQuadrupole_YZ*gqyz7));

    dpdx -= sxi*(atomJ.labFrameQuadrupole_XX*gqxx5 + atomJ.labFrameQuadrupole_YY*gqyy5 + atomJ.labFrameQuadrupole_ZZ*gqzz5 + 2.0f*(atomJ.labFrameQuadrupole_XY*gqxy5 + atomJ.labFrameQuadrupole_XZ*gqxz5 + atomJ.labFrameQuadrupole_YZ*gqyz5)) +
            syi*(atomJ.labFrameQuadrupole_XX*gqxx6 + atomJ.labFrameQuadrupole_YY*gqyy6 + atomJ.labFrameQuadrupole_ZZ*gqzz6 + 2.0f*(atomJ.labFrameQuadrupole_XY*gqxy6 + atomJ.labFrameQuadrupole_XZ*gqxz6 + atomJ.labFrameQuadrupole_YZ*gqyz6)) +
            szi*(atomJ.labFrameQuadrupole_XX*gqxx7 + atomJ.labFrameQuadrupole_YY*gqyy7 + atomJ.labFrameQuadrupole_ZZ*gqzz7 + 2.0f*(atomJ.labFrameQuadrupole_XY*gqxy7 + atomJ.labFrameQuadrupole_XZ*gqxz7 + atomJ.labFrameQuadrupole_YZ*gqyz7));

    dpdy += sxk*(atomI.labFrameQuadrupole_XX*gqxx6 + atomI.labFrameQuadrupole_YY*gqyy6 + atomI.labFrameQuadrupole_ZZ*gqzz6 + 2.0f*(atomI.labFrameQuadrupole_XY*gqxy6 + atomI.labFrameQuadrupole_XZ*gqxz6 + atomI.labFrameQuadrupole_YZ*gqyz6)) +
            syk*(atomI.labFrameQuadrupole_XX*gqxx8 + atomI.labFrameQuadrupole_YY*gqyy8 + atomI.labFrameQuadrupole_ZZ*gqzz8 + 2.0f*(atomI.labFrameQuadrupole_XY*gqxy8 + atomI.labFrameQuadrupole_XZ*gqxz8 + atomI.labFrameQuadrupole_YZ*gqyz8)) +
            szk*(atomI.labFrameQuadrupole_XX*gqxx9 + atomI.labFrameQuadrupole_YY*gqyy9 + atomI.labFrameQuadrupole_ZZ*gqzz9 + 2.0f*(atomI.labFrameQuadrupole_XY*gqxy9 + atomI.labFrameQuadrupole_XZ*gqxz9 + atomI.labFrameQuadrupole_YZ*gqyz9));

    dpdy -= sxi*(atomJ.labFrameQuadrupole_XX*gqxx6 + atomJ.labFrameQuadrupole_YY*gqyy6 + atomJ.labFrameQuadrupole_ZZ*gqzz6 + 2.0f*(atomJ.labFrameQuadrupole_XY*gqxy6 + atomJ.labFrameQuadrupole_XZ*gqxz6 + atomJ.labFrameQuadrupole_YZ*gqyz6)) +
            syi*(atomJ.labFrameQuadrupole_XX*gqxx8 + atomJ.labFrameQuadrupole_YY*gqyy8 + atomJ.labFrameQuadrupole_ZZ*gqzz8 + 2.0f*(atomJ.labFrameQuadrupole_XY*gqxy8 + atomJ.labFrameQuadrupole_XZ*gqxz8 + atomJ.labFrameQuadrupole_YZ*gqyz8)) +
            szi*(atomJ.labFrameQuadrupole_XX*gqxx9 + atomJ.labFrameQuadrupole_YY*gqyy9 + atomJ.labFrameQuadrupole_ZZ*gqzz9 + 2.0f*(atomJ.labFrameQuadrupole_XY*gqxy9 + atomJ.labFrameQuadrupole_XZ*gqxz9 + atomJ.labFrameQuadrupole_YZ*gqyz9));

    dpdz -= sxi*(atomJ.labFrameQuadrupole_XX*gqxx7 + atomJ.labFrameQuadrupole_YY*gqyy7 + atomJ.labFrameQuadrupole_ZZ*gqzz7 + 2.0f*(atomJ.labFrameQuadrupole_XY*gqxy7 + atomJ.labFrameQuadrupole_XZ*gqxz7 + atomJ.labFrameQuadrupole_YZ*gqyz7)) +
            syi*(atomJ.labFrameQuadrupole_XX*gqxx9 + atomJ.labFrameQuadrupole_YY*gqyy9 + atomJ.labFrameQuadrupole_ZZ*gqzz9 + 2.0f*(atomJ.labFrameQuadrupole_XY*gqxy9 + atomJ.labFrameQuadrupole_XZ*gqxz9 + atomJ.labFrameQuadrupole_YZ*gqyz9)) +
            szi*(atomJ.labFrameQuadrupole_XX*gqxx10 + atomJ.labFrameQuadrupole_YY*gqyy10 + atomJ.labFrameQuadrupole_ZZ*gqzz10 + 2.0f*(atomJ.labFrameQuadrupole_XY*gqxy10 + atomJ.labFrameQuadrupole_XZ*gqxz10 + atomJ.labFrameQuadrupole_YZ*gqyz10));

    dpdz += sxk*(atomI.labFrameQuadrupole_XX*gqxx7  + atomI.labFrameQuadrupole_YY*gqyy7  + atomI.labFrameQuadrupole_ZZ*gqzz7  + 2.0f*(atomI.labFrameQuadrupole_XY*gqxy7  + atomI.labFrameQuadrupole_XZ*gqxz7  + atomI.labFrameQuadrupole_YZ*gqyz7)) +
            syk*(atomI.labFrameQuadrupole_XX*gqxx9  + atomI.labFrameQuadrupole_YY*gqyy9  + atomI.labFrameQuadrupole_ZZ*gqzz9  + 2.0f*(atomI.labFrameQuadrupole_XY*gqxy9  + atomI.labFrameQuadrupole_XZ*gqxz9  + atomI.labFrameQuadrupole_YZ*gqyz9)) +
            szk*(atomI.labFrameQuadrupole_XX*gqxx10 + atomI.labFrameQuadrupole_YY*gqyy10 + atomI.labFrameQuadrupole_ZZ*gqzz10 + 2.0f*(atomI.labFrameQuadrupole_XY*gqxy10 + atomI.labFrameQuadrupole_XZ*gqxz10 + atomI.labFrameQuadrupole_YZ*gqyz10));

#endif
#endif

    // Born radii derivatives of the unweighted 2nd reaction;
    // potential gradient tensor;

#if defined B1
    float gc25       = b01  + xr2*b02;
    float gc26       = xr*yr*b02;
    float gc27       = xr*zr*b02;
    float gc28       = b01  + yr2*b02;
    float gc29       = yr*zr*b02;
    float gc30       = b01  + zr2*b02;
    dsumdrB1        += atomI.q*(atomJ.labFrameQuadrupole_XX*gc25 + atomJ.labFrameQuadrupole_YY*gc28 + atomJ.labFrameQuadrupole_ZZ*gc30 + 2.0f*(atomJ.labFrameQuadrupole_XY*gc26 + atomJ.labFrameQuadrupole_XZ*gc27 + atomJ.labFrameQuadrupole_YZ*gc29));
    dsumdrB1        += atomJ.q*(atomI.labFrameQuadrupole_XX*gc25 + atomI.labFrameQuadrupole_YY*gc28 + atomI.labFrameQuadrupole_ZZ*gc30 + 2.0f*(atomI.labFrameQuadrupole_XY*gc26 + atomI.labFrameQuadrupole_XZ*gc27 + atomI.labFrameQuadrupole_YZ*gc29));
#endif
#if defined B1 || defined B2
    float gux25      = xr*(3.0f*b11 + xr2*b12);
    float gux26      = yr*(b11 + xr2*b12);
    float gux27      = zr*(b11 + xr2*b12);
    float gux28      = xr*(b11 + yr2*b12);
    float gux29      = zr*xr*yr*b12;
    float gux30      = xr*(b11 + zr2*b12);
    float guy25      = yr*(b11 + xr2*b12);
    float guy26      = xr*(b11 + yr2*b12);
    float guy27      = gux29;
    float guy28      = yr*(3.0f*b11 + yr2*b12);
    float guy29      = zr*(b11 + yr2*b12);
    float guy30      = yr*(b11 + zr2*b12);
    float guz25      = zr*(b11 + xr2*b12);
    float guz26      = gux29;
    float guz27      = xr*(b11 + zr2*b12);
    float guz28      = zr*(b11 + yr2*b12);
    float guz29      = yr*(b11 + zr2*b12);
    float guz30      = zr*(3.0f*b11 + zr2*b12);
#endif
#if defined B2
    dsumdrB2 -= sxi*(atomJ.labFrameQuadrupole_XX*gux25 + atomJ.labFrameQuadrupole_YY*gux28 + atomJ.labFrameQuadrupole_ZZ*gux30 + 2.0f*(atomJ.labFrameQuadrupole_XY*gux26 + atomJ.labFrameQuadrupole_XZ*gux27 + atomJ.labFrameQuadrupole_YZ*gux29)) +
                syi*(atomJ.labFrameQuadrupole_XX*guy25 + atomJ.labFrameQuadrupole_YY*guy28 + atomJ.labFrameQuadrupole_ZZ*guy30 + 2.0f*(atomJ.labFrameQuadrupole_XY*guy26 + atomJ.labFrameQuadrupole_XZ*guy27 + atomJ.labFrameQuadrupole_YZ*guy29)) +
                szi*(atomJ.labFrameQuadrupole_XX*guz25 + atomJ.labFrameQuadrupole_YY*guz28 + atomJ.labFrameQuadrupole_ZZ*guz30 + 2.0f*(atomJ.labFrameQuadrupole_XY*guz26 + atomJ.labFrameQuadrupole_XZ*guz27 + atomJ.labFrameQuadrupole_YZ*guz29));
    dsumdrB2 += sxk*(atomI.labFrameQuadrupole_XX*gux25 + atomI.labFrameQuadrupole_YY*gux28 + atomI.labFrameQuadrupole_ZZ*gux30 + 2.0f*(atomI.labFrameQuadrupole_XY*gux26 + atomI.labFrameQuadrupole_XZ*gux27 + atomI.labFrameQuadrupole_YZ*gux29)) +
                syk*(atomI.labFrameQuadrupole_XX*guy25 + atomI.labFrameQuadrupole_YY*guy28 + atomI.labFrameQuadrupole_ZZ*guy30 + 2.0f*(atomI.labFrameQuadrupole_XY*guy26 + atomI.labFrameQuadrupole_XZ*guy27 + atomI.labFrameQuadrupole_YZ*guy29)) +
                szk*(atomI.labFrameQuadrupole_XX*guz25 + atomI.labFrameQuadrupole_YY*guz28 + atomI.labFrameQuadrupole_ZZ*guz30 + 2.0f*(atomI.labFrameQuadrupole_XY*guz26 + atomI.labFrameQuadrupole_XZ*guz27 + atomI.labFrameQuadrupole_YZ*guz29));

#endif
#if defined B1
    dsumdrB1 -= atomI.labFrameDipole[0]*(atomJ.labFrameQuadrupole_XX*gux25 + atomJ.labFrameQuadrupole_YY*gux28 + atomJ.labFrameQuadrupole_ZZ*gux30 + 2.0f*(atomJ.labFrameQuadrupole_XY*gux26 + atomJ.labFrameQuadrupole_XZ*gux27 + atomJ.labFrameQuadrupole_YZ*gux29)) +
                atomI.labFrameDipole[1]*(atomJ.labFrameQuadrupole_XX*guy25 + atomJ.labFrameQuadrupole_YY*guy28 + atomJ.labFrameQuadrupole_ZZ*guy30 + 2.0f*(atomJ.labFrameQuadrupole_XY*guy26 + atomJ.labFrameQuadrupole_XZ*guy27 + atomJ.labFrameQuadrupole_YZ*guy29)) +
                atomI.labFrameDipole[2]*(atomJ.labFrameQuadrupole_XX*guz25 + atomJ.labFrameQuadrupole_YY*guz28 + atomJ.labFrameQuadrupole_ZZ*guz30 + 2.0f*(atomJ.labFrameQuadrupole_XY*guz26 + atomJ.labFrameQuadrupole_XZ*guz27 + atomJ.labFrameQuadrupole_YZ*guz29));
    dsumdrB1 += atomJ.labFrameDipole[0]*(atomI.labFrameQuadrupole_XX*gux25 + atomI.labFrameQuadrupole_YY*gux28 + atomI.labFrameQuadrupole_ZZ*gux30 + 2.0f*(atomI.labFrameQuadrupole_XY*gux26 + atomI.labFrameQuadrupole_XZ*gux27 + atomI.labFrameQuadrupole_YZ*gux29)) +
                atomJ.labFrameDipole[1]*(atomI.labFrameQuadrupole_XX*guy25 + atomI.labFrameQuadrupole_YY*guy28 + atomI.labFrameQuadrupole_ZZ*guy30 + 2.0f*(atomI.labFrameQuadrupole_XY*guy26 + atomI.labFrameQuadrupole_XZ*guy27 + atomI.labFrameQuadrupole_YZ*guy29)) +
                atomJ.labFrameDipole[2]*(atomI.labFrameQuadrupole_XX*guz25 + atomI.labFrameQuadrupole_YY*guz28 + atomI.labFrameQuadrupole_ZZ*guz30 + 2.0f*(atomI.labFrameQuadrupole_XY*guz26 + atomI.labFrameQuadrupole_XZ*guz27 + atomI.labFrameQuadrupole_YZ*guz29));

    float gqxx25     = 2.0f*b20  + xr2*(5.0f*b21 + xr2*b22);
    float gqxx26     = yr*xr*(2.0f*b21 + xr2*b22);
    float gqxx27     = zr*xr*(2.0f*b21 + xr2*b22);
    float gqxx28     = xr2*(b21 + yr2*b22);
    float gqxx29     = zr*yr*xr2*b22;
    float gqxx30     = xr2*(b21 + zr2*b22);
    float gqyy25     = yr2*(b21 + xr2*b22);
    float gqyy26     = xr*yr*(2.0f*b21 + yr2*b22);
    float gqyy27     = xr*zr*yr2*b22;
    float gqyy28     = 2.0f*b20  + yr2*(5.0f*b21 + yr2*b22);
    float gqyy29     = yr*zr*(2.0f*b21 + yr2*b22);
    float gqyy30     = yr2*(b21 + zr2*b22);
    float gqzz25     = zr2*(b21 + xr2*b22);
    float gqzz26     = xr*yr*zr2*b22;
    float gqzz27     = xr*zr*(2.0f*b21 + zr2*b22);
    float gqzz28     = zr2*(b21 + yr2*b22);
    float gqzz29     = yr*zr*(2.0f*b21 + zr2*b22);
    float gqzz30     = 2.0f*b20  + zr2*(5.0f*b21 + zr2*b22);
    float gqxy25     = xr*yr*(3.0f*b21  + xr2*b22);
    float gqxy26     = b20  + (xr2 + yr2)*b21  + xr2*yr2*b22;
    float gqxy27     = zr*yr*(b21 + xr2*b22);
    float gqxy28     = xr*yr*(3.0f*b21 + yr2*b22);
    float gqxy29     = zr*xr*(b21 + yr2*b22);
    float gqxy30     = xr*yr*(b21 + zr2*b22);
    float gqxz25     = xr*zr*(3.0f*b21 + xr2*b22);
    float gqxz26     = yr*zr*(b21 + xr2*b22);
    float gqxz27     = b20  + (xr2 + zr2)*b21  + xr2*zr2*b22;
    float gqxz28     = xr*zr*(b21 + yr2*b22);
    float gqxz29     = xr*yr*(b21 + zr2*b22);
    float gqxz30     = xr*zr*(3.0f*b21 + zr2*b22);
    float gqyz25     = zr*yr*(b21 + xr2*b22);
    float gqyz26     = xr*zr*(b21 + yr2*b22);
    float gqyz27     = xr*yr*(b21 + zr2*b22);
    float gqyz28     = yr*zr*(3.0f*b21 + yr2*b22);
    float gqyz29     = b20  + (yr2 + zr2)*b21  + yr2*zr2*b22;
    float gqyz30     = yr*zr*(3.0f*b21 + zr2*b22);

    dsumdrB1 +=
        atomI.labFrameQuadrupole_XX*(atomJ.labFrameQuadrupole_XX*gqxx25 + atomJ.labFrameQuadrupole_YY*gqxx28 + atomJ.labFrameQuadrupole_ZZ*gqxx30 + 2.0f*(atomJ.labFrameQuadrupole_XY*gqxx26 + atomJ.labFrameQuadrupole_XZ*gqxx27 + atomJ.labFrameQuadrupole_YZ*gqxx29)) +
        atomI.labFrameQuadrupole_YY*(atomJ.labFrameQuadrupole_XX*gqyy25 + atomJ.labFrameQuadrupole_YY*gqyy28 + atomJ.labFrameQuadrupole_ZZ*gqyy30 + 2.0f*(atomJ.labFrameQuadrupole_XY*gqyy26 + atomJ.labFrameQuadrupole_XZ*gqyy27 + atomJ.labFrameQuadrupole_YZ*gqyy29)) +
        atomI.labFrameQuadrupole_ZZ*(atomJ.labFrameQuadrupole_XX*gqzz25 + atomJ.labFrameQuadrupole_YY*gqzz28 + atomJ.labFrameQuadrupole_ZZ*gqzz30 + 2.0f*(atomJ.labFrameQuadrupole_XY*gqzz26 + atomJ.labFrameQuadrupole_XZ*gqzz27 + atomJ.labFrameQuadrupole_YZ*gqzz29));

    dsumdrB1 += 2.0f*(
        atomI.labFrameQuadrupole_XY*(atomJ.labFrameQuadrupole_XX*gqxy25 + atomJ.labFrameQuadrupole_YY*gqxy28 + atomJ.labFrameQuadrupole_ZZ*gqxy30 + 2.0f*(atomJ.labFrameQuadrupole_XY*gqxy26 + atomJ.labFrameQuadrupole_XZ*gqxy27 + atomJ.labFrameQuadrupole_YZ*gqxy29)) +
        atomI.labFrameQuadrupole_XZ*(atomJ.labFrameQuadrupole_XX*gqxz25 + atomJ.labFrameQuadrupole_YY*gqxz28 + atomJ.labFrameQuadrupole_ZZ*gqxz30 + 2.0f*(atomJ.labFrameQuadrupole_XY*gqxz26 + atomJ.labFrameQuadrupole_XZ*gqxz27 + atomJ.labFrameQuadrupole_YZ*gqxz29)) +
        atomI.labFrameQuadrupole_YZ*(atomJ.labFrameQuadrupole_XX*gqyz25 + atomJ.labFrameQuadrupole_YY*gqyz28 + atomJ.labFrameQuadrupole_ZZ*gqyz30 + 2.0f*(atomJ.labFrameQuadrupole_XY*gqyz26 + atomJ.labFrameQuadrupole_XZ*gqyz27 + atomJ.labFrameQuadrupole_YZ*gqyz29)));

    dsumdrB1 +=
        atomI.labFrameQuadrupole_XX*(atomJ.labFrameQuadrupole_XX*gqxx25 + atomJ.labFrameQuadrupole_YY*gqyy25 + atomJ.labFrameQuadrupole_ZZ*gqzz25 + 2.0f*(atomJ.labFrameQuadrupole_XY*gqxy25 + atomJ.labFrameQuadrupole_XZ*gqxz25 + atomJ.labFrameQuadrupole_YZ*gqyz25)) +
        atomI.labFrameQuadrupole_YY*(atomJ.labFrameQuadrupole_XX*gqxx28 + atomJ.labFrameQuadrupole_YY*gqyy28 + atomJ.labFrameQuadrupole_ZZ*gqzz28 + 2.0f*(atomJ.labFrameQuadrupole_XY*gqxy28 + atomJ.labFrameQuadrupole_XZ*gqxz28 + atomJ.labFrameQuadrupole_YZ*gqyz28)) +
        atomI.labFrameQuadrupole_ZZ*(atomJ.labFrameQuadrupole_XX*gqxx30 + atomJ.labFrameQuadrupole_YY*gqyy30 + atomJ.labFrameQuadrupole_ZZ*gqzz30 + 2.0f*(atomJ.labFrameQuadrupole_XY*gqxy30 + atomJ.labFrameQuadrupole_XZ*gqxz30 + atomJ.labFrameQuadrupole_YZ*gqyz30));

    dsumdrB1 += 2.0f*(
        atomI.labFrameQuadrupole_XY*(atomJ.labFrameQuadrupole_XX*gqxx26 + atomJ.labFrameQuadrupole_YY*gqyy26 + atomJ.labFrameQuadrupole_ZZ*gqzz26 + 2.0f*(atomJ.labFrameQuadrupole_XY*gqxy26 + atomJ.labFrameQuadrupole_XZ*gqxz26 + atomJ.labFrameQuadrupole_YZ*gqyz26)) +
        atomI.labFrameQuadrupole_XZ*(atomJ.labFrameQuadrupole_XX*gqxx27 + atomJ.labFrameQuadrupole_YY*gqyy27 + atomJ.labFrameQuadrupole_ZZ*gqzz27 + 2.0f*(atomJ.labFrameQuadrupole_XY*gqxy27 + atomJ.labFrameQuadrupole_XZ*gqxz27 + atomJ.labFrameQuadrupole_YZ*gqyz27)) +
        atomI.labFrameQuadrupole_YZ*(atomJ.labFrameQuadrupole_XX*gqxx29 + atomJ.labFrameQuadrupole_YY*gqyy29 + atomJ.labFrameQuadrupole_ZZ*gqzz29 + 2.0f*(atomJ.labFrameQuadrupole_XY*gqxy29 + atomJ.labFrameQuadrupole_XZ*gqxz29 + atomJ.labFrameQuadrupole_YZ*gqyz29)));

    dsumdrB1          *= 0.5f;
    atomI.dBornRadius += atomJ.bornRadius*dsumdrB1;
    atomJ.dBornRadius += atomI.bornRadius*dsumdrB1;
#endif

    // unweighted 3rd reaction potential gradient tensor;

#if defined F1
    float gc11       = xr*(3.0f*a02 + xr2*a03);
    float gc12       = yr*(a02 + xr2*a03);
    float gc13       = zr*(a02 + xr2*a03);
    float gc14       = xr*(a02 + yr2*a03);
    float gc15       = xr*yr*zr*a03;
    float gc16       = xr*(a02 + zr2*a03);
    float gc17       = yr*(3.0f*a02 + yr2*a03);
    float gc18       = zr*(a02 + yr2*a03);
    float gc19       = yr*(a02 + zr2*a03);
    float gc20       = zr*(3.0f*a02 + zr2*a03);
    dedx            += atomI.q*(atomJ.labFrameQuadrupole_XX*gc11 + atomJ.labFrameQuadrupole_YY*gc14 + atomJ.labFrameQuadrupole_ZZ*gc16 + 2.0f*(atomJ.labFrameQuadrupole_XY*gc12 + atomJ.labFrameQuadrupole_XZ*gc13 + atomJ.labFrameQuadrupole_YZ*gc15));
    dedx            += atomJ.q*(atomI.labFrameQuadrupole_XX*gc11 + atomI.labFrameQuadrupole_YY*gc14 + atomI.labFrameQuadrupole_ZZ*gc16 + 2.0f*(atomI.labFrameQuadrupole_XY*gc12 + atomI.labFrameQuadrupole_XZ*gc13 + atomI.labFrameQuadrupole_YZ*gc15));
    dedy            += atomI.q*(atomJ.labFrameQuadrupole_XX*gc12 + atomJ.labFrameQuadrupole_YY*gc17 + atomJ.labFrameQuadrupole_ZZ*gc19 + 2.0f*(atomJ.labFrameQuadrupole_XY*gc14 + atomJ.labFrameQuadrupole_XZ*gc15 + atomJ.labFrameQuadrupole_YZ*gc18));
    dedy            += atomJ.q*(atomI.labFrameQuadrupole_XX*gc12 + atomI.labFrameQuadrupole_YY*gc17 + atomI.labFrameQuadrupole_ZZ*gc19 + 2.0f*(atomI.labFrameQuadrupole_XY*gc14 + atomI.labFrameQuadrupole_XZ*gc15 + atomI.labFrameQuadrupole_YZ*gc18));
    dedz            += atomI.q*(atomJ.labFrameQuadrupole_XX*gc13 + atomJ.labFrameQuadrupole_YY*gc18 + atomJ.labFrameQuadrupole_ZZ*gc20 + 2.0f*(atomJ.labFrameQuadrupole_XY*gc15 + atomJ.labFrameQuadrupole_XZ*gc16 + atomJ.labFrameQuadrupole_YZ*gc19));
    dedz            += atomJ.q*(atomI.labFrameQuadrupole_XX*gc13 + atomI.labFrameQuadrupole_YY*gc18 + atomI.labFrameQuadrupole_ZZ*gc20 + 2.0f*(atomI.labFrameQuadrupole_XY*gc15 + atomI.labFrameQuadrupole_XZ*gc16 + atomI.labFrameQuadrupole_YZ*gc19));
#endif
#if defined F1 || defined F2
    float gux11      = 3.0f*a11  + xr2*(6.0f*a12 + xr2*a13);
    float gux12      = xr*yr*(3.0f*a12 + xr2*a13);
    float gux13      = xr*zr*(3.0f*a12 + xr2*a13);
    float gux14      = a11  + (xr2 + yr2)*a12  + xr2*yr2*a13;
    float gux15      = yr*zr*(a12 + xr2*a13);
    float gux16      = a11  + (xr2 + zr2)*a12  + xr2*zr2*a13;
    float gux17      = xr*yr*(3.0f*a12 + yr2*a13);
    float gux18      = xr*zr*(a12 + yr2*a13);
    float gux19      = xr*yr*(a12 + zr2*a13);
    float gux20      = xr*zr*(3.0f*a12 + zr2*a13);
    float guy11      = gux12;
    float guy12      = gux14;
    float guy13      = gux15;
    float guy14      = gux17;
    float guy15      = gux18;
    float guy16      = gux19;
    float guy17      = 3.0f*a11  + yr2*(6.0f*a12 + yr2*a13);
    float guy18      = yr*zr*(3.0f*a12 + yr2*a13);
    float guy19      = a11  + (yr2 + zr2)*a12  + yr2*zr2*a13;
    float guy20      = yr*zr*(3.0f*a12 + zr2*a13);
    float guz11      = gux13;
    float guz12      = gux15;
    float guz13      = gux16;
    float guz14      = gux18;
    float guz15      = gux19;
    float guz16      = gux20;
    float guz17      = guy18;
    float guz18      = guy19;
    float guz19      = guy20;
    float guz20      = 3.0f*a11  + zr2*(6.0f*a12 + zr2*a13);
#if defined F1
    dedx            -= atomI.labFrameDipole[0]*(atomJ.labFrameQuadrupole_XX*gux11 + atomJ.labFrameQuadrupole_YY*gux14 + atomJ.labFrameQuadrupole_ZZ*gux16 + 2.0f*(atomJ.labFrameQuadrupole_XY*gux12 + atomJ.labFrameQuadrupole_XZ*gux13 + atomJ.labFrameQuadrupole_YZ*gux15)) +
                       atomI.labFrameDipole[1]*(atomJ.labFrameQuadrupole_XX*guy11 + atomJ.labFrameQuadrupole_YY*guy14 + atomJ.labFrameQuadrupole_ZZ*guy16 + 2.0f*(atomJ.labFrameQuadrupole_XY*guy12 + atomJ.labFrameQuadrupole_XZ*guy13 + atomJ.labFrameQuadrupole_YZ*guy15)) +
                       atomI.labFrameDipole[2]*(atomJ.labFrameQuadrupole_XX*guz11 + atomJ.labFrameQuadrupole_YY*guz14 + atomJ.labFrameQuadrupole_ZZ*guz16 + 2.0f*(atomJ.labFrameQuadrupole_XY*guz12 + atomJ.labFrameQuadrupole_XZ*guz13 + atomJ.labFrameQuadrupole_YZ*guz15));

    dedx            += atomJ.labFrameDipole[0]*(atomI.labFrameQuadrupole_XX*gux11 + atomI.labFrameQuadrupole_YY*gux14 + atomI.labFrameQuadrupole_ZZ*gux16 + 2.0f*(atomI.labFrameQuadrupole_XY*gux12 + atomI.labFrameQuadrupole_XZ*gux13 + atomI.labFrameQuadrupole_YZ*gux15)) +
                       atomJ.labFrameDipole[1]*(atomI.labFrameQuadrupole_XX*guy11 + atomI.labFrameQuadrupole_YY*guy14 + atomI.labFrameQuadrupole_ZZ*guy16 + 2.0f*(atomI.labFrameQuadrupole_XY*guy12 + atomI.labFrameQuadrupole_XZ*guy13 + atomI.labFrameQuadrupole_YZ*guy15)) +
                       atomJ.labFrameDipole[2]*(atomI.labFrameQuadrupole_XX*guz11 + atomI.labFrameQuadrupole_YY*guz14 + atomI.labFrameQuadrupole_ZZ*guz16 + 2.0f*(atomI.labFrameQuadrupole_XY*guz12 + atomI.labFrameQuadrupole_XZ*guz13 + atomI.labFrameQuadrupole_YZ*guz15));

    dedy            -= atomI.labFrameDipole[0]*(atomJ.labFrameQuadrupole_XX*gux12 + atomJ.labFrameQuadrupole_YY*gux17 + atomJ.labFrameQuadrupole_ZZ*gux19 + 2.0f*(atomJ.labFrameQuadrupole_XY*gux14 + atomJ.labFrameQuadrupole_XZ*gux15 + atomJ.labFrameQuadrupole_YZ*gux18)) +
                       atomI.labFrameDipole[1]*(atomJ.labFrameQuadrupole_XX*guy12 + atomJ.labFrameQuadrupole_YY*guy17 + atomJ.labFrameQuadrupole_ZZ*guy19 + 2.0f*(atomJ.labFrameQuadrupole_XY*guy14 + atomJ.labFrameQuadrupole_XZ*guy15 + atomJ.labFrameQuadrupole_YZ*guy18)) +
                       atomI.labFrameDipole[2]*(atomJ.labFrameQuadrupole_XX*guz12 + atomJ.labFrameQuadrupole_YY*guz17 + atomJ.labFrameQuadrupole_ZZ*guz19 + 2.0f*(atomJ.labFrameQuadrupole_XY*guz14 + atomJ.labFrameQuadrupole_XZ*guz15 + atomJ.labFrameQuadrupole_YZ*guz18));

    dedy            += atomJ.labFrameDipole[0]*(atomI.labFrameQuadrupole_XX*gux12 + atomI.labFrameQuadrupole_YY*gux17 + atomI.labFrameQuadrupole_ZZ*gux19 + 2.0f*(atomI.labFrameQuadrupole_XY*gux14 + atomI.labFrameQuadrupole_XZ*gux15 + atomI.labFrameQuadrupole_YZ*gux18)) +
                       atomJ.labFrameDipole[1]*(atomI.labFrameQuadrupole_XX*guy12 + atomI.labFrameQuadrupole_YY*guy17 + atomI.labFrameQuadrupole_ZZ*guy19 + 2.0f*(atomI.labFrameQuadrupole_XY*guy14 + atomI.labFrameQuadrupole_XZ*guy15 + atomI.labFrameQuadrupole_YZ*guy18)) +
                       atomJ.labFrameDipole[2]*(atomI.labFrameQuadrupole_XX*guz12 + atomI.labFrameQuadrupole_YY*guz17 + atomI.labFrameQuadrupole_ZZ*guz19 + 2.0f*(atomI.labFrameQuadrupole_XY*guz14 + atomI.labFrameQuadrupole_XZ*guz15 + atomI.labFrameQuadrupole_YZ*guz18));

    dedz            -= atomI.labFrameDipole[0]*(atomJ.labFrameQuadrupole_XX*gux13 + atomJ.labFrameQuadrupole_YY*gux18 + atomJ.labFrameQuadrupole_ZZ*gux20 + 2.0f*(atomJ.labFrameQuadrupole_XY*gux15 + atomJ.labFrameQuadrupole_XZ*gux16 + atomJ.labFrameQuadrupole_YZ*gux19)) +
                       atomI.labFrameDipole[1]*(atomJ.labFrameQuadrupole_XX*guy13 + atomJ.labFrameQuadrupole_YY*guy18 + atomJ.labFrameQuadrupole_ZZ*guy20 + 2.0f*(atomJ.labFrameQuadrupole_XY*guy15 + atomJ.labFrameQuadrupole_XZ*guy16 + atomJ.labFrameQuadrupole_YZ*guy19)) +
                       atomI.labFrameDipole[2]*(atomJ.labFrameQuadrupole_XX*guz13 + atomJ.labFrameQuadrupole_YY*guz18 + atomJ.labFrameQuadrupole_ZZ*guz20 + 2.0f*(atomJ.labFrameQuadrupole_XY*guz15 + atomJ.labFrameQuadrupole_XZ*guz16 + atomJ.labFrameQuadrupole_YZ*guz19));

    dedz            += atomJ.labFrameDipole[0]*(atomI.labFrameQuadrupole_XX*gux13 + atomI.labFrameQuadrupole_YY*gux18 + atomI.labFrameQuadrupole_ZZ*gux20 + 2.0f*(atomI.labFrameQuadrupole_XY*gux15 + atomI.labFrameQuadrupole_XZ*gux16 + atomI.labFrameQuadrupole_YZ*gux19)) +
                       atomJ.labFrameDipole[1]*(atomI.labFrameQuadrupole_XX*guy13 + atomI.labFrameQuadrupole_YY*guy18 + atomI.labFrameQuadrupole_ZZ*guy20 + 2.0f*(atomI.labFrameQuadrupole_XY*guy15 + atomI.labFrameQuadrupole_XZ*guy16 + atomI.labFrameQuadrupole_YZ*guy19)) +
                       atomJ.labFrameDipole[2]*(atomI.labFrameQuadrupole_XX*guz13 + atomI.labFrameQuadrupole_YY*guz18 + atomI.labFrameQuadrupole_ZZ*guz20 + 2.0f*(atomI.labFrameQuadrupole_XY*guz15 + atomI.labFrameQuadrupole_XZ*guz16 + atomI.labFrameQuadrupole_YZ*guz19));

#endif
#if defined F2
    dpdx -= sxi*(atomJ.labFrameQuadrupole_XX*gux11  + atomJ.labFrameQuadrupole_YY*gux14  + atomJ.labFrameQuadrupole_ZZ*gux16  + 2.0f*(atomJ.labFrameQuadrupole_XY*gux12  + atomJ.labFrameQuadrupole_XZ*gux13  + atomJ.labFrameQuadrupole_YZ*gux15)) +
            syi*(atomJ.labFrameQuadrupole_XX*guy11  + atomJ.labFrameQuadrupole_YY*guy14  + atomJ.labFrameQuadrupole_ZZ*guy16  + 2.0f*(atomJ.labFrameQuadrupole_XY*guy12  + atomJ.labFrameQuadrupole_XZ*guy13  + atomJ.labFrameQuadrupole_YZ*guy15)) +
            szi*(atomJ.labFrameQuadrupole_XX*guz11  + atomJ.labFrameQuadrupole_YY*guz14  + atomJ.labFrameQuadrupole_ZZ*guz16  + 2.0f*(atomJ.labFrameQuadrupole_XY*guz12  + atomJ.labFrameQuadrupole_XZ*guz13  + atomJ.labFrameQuadrupole_YZ*guz15));

    dpdx += sxk*(atomI.labFrameQuadrupole_XX*gux11  + atomI.labFrameQuadrupole_YY*gux14  + atomI.labFrameQuadrupole_ZZ*gux16  + 2.0f*(atomI.labFrameQuadrupole_XY*gux12  + atomI.labFrameQuadrupole_XZ*gux13  + atomI.labFrameQuadrupole_YZ*gux15)) +
            syk*(atomI.labFrameQuadrupole_XX*guy11  + atomI.labFrameQuadrupole_YY*guy14  + atomI.labFrameQuadrupole_ZZ*guy16  + 2.0f*(atomI.labFrameQuadrupole_XY*guy12  + atomI.labFrameQuadrupole_XZ*guy13  + atomI.labFrameQuadrupole_YZ*guy15)) +
            szk*(atomI.labFrameQuadrupole_XX*guz11  + atomI.labFrameQuadrupole_YY*guz14  + atomI.labFrameQuadrupole_ZZ*guz16  + 2.0f*(atomI.labFrameQuadrupole_XY*guz12  + atomI.labFrameQuadrupole_XZ*guz13  + atomI.labFrameQuadrupole_YZ*guz15));

    dpdy -= sxi*(atomJ.labFrameQuadrupole_XX*gux12  + atomJ.labFrameQuadrupole_YY*gux17  + atomJ.labFrameQuadrupole_ZZ*gux19  + 2.0f*(atomJ.labFrameQuadrupole_XY*gux14  + atomJ.labFrameQuadrupole_XZ*gux15  + atomJ.labFrameQuadrupole_YZ*gux18)) +
            syi*(atomJ.labFrameQuadrupole_XX*guy12  + atomJ.labFrameQuadrupole_YY*guy17  + atomJ.labFrameQuadrupole_ZZ*guy19  + 2.0f*(atomJ.labFrameQuadrupole_XY*guy14  + atomJ.labFrameQuadrupole_XZ*guy15  + atomJ.labFrameQuadrupole_YZ*guy18)) +
            szi*(atomJ.labFrameQuadrupole_XX*guz12  + atomJ.labFrameQuadrupole_YY*guz17  + atomJ.labFrameQuadrupole_ZZ*guz19  + 2.0f*(atomJ.labFrameQuadrupole_XY*guz14  + atomJ.labFrameQuadrupole_XZ*guz15  + atomJ.labFrameQuadrupole_YZ*guz18));

    dpdy += sxk*(atomI.labFrameQuadrupole_XX*gux12  + atomI.labFrameQuadrupole_YY*gux17  + atomI.labFrameQuadrupole_ZZ*gux19  + 2.0f*(atomI.labFrameQuadrupole_XY*gux14  + atomI.labFrameQuadrupole_XZ*gux15  + atomI.labFrameQuadrupole_YZ*gux18)) +
            syk*(atomI.labFrameQuadrupole_XX*guy12  + atomI.labFrameQuadrupole_YY*guy17  + atomI.labFrameQuadrupole_ZZ*guy19  + 2.0f*(atomI.labFrameQuadrupole_XY*guy14  + atomI.labFrameQuadrupole_XZ*guy15  + atomI.labFrameQuadrupole_YZ*guy18)) +
            szk*(atomI.labFrameQuadrupole_XX*guz12  + atomI.labFrameQuadrupole_YY*guz17  + atomI.labFrameQuadrupole_ZZ*guz19  + 2.0f*(atomI.labFrameQuadrupole_XY*guz14  + atomI.labFrameQuadrupole_XZ*guz15  + atomI.labFrameQuadrupole_YZ*guz18));

    dpdz -= sxi*(atomJ.labFrameQuadrupole_XX*gux13  + atomJ.labFrameQuadrupole_YY*gux18  + atomJ.labFrameQuadrupole_ZZ*gux20  + 2.0f*(atomJ.labFrameQuadrupole_XY*gux15  + atomJ.labFrameQuadrupole_XZ*gux16  + atomJ.labFrameQuadrupole_YZ*gux19)) +
            syi*(atomJ.labFrameQuadrupole_XX*guy13  + atomJ.labFrameQuadrupole_YY*guy18  + atomJ.labFrameQuadrupole_ZZ*guy20  + 2.0f*(atomJ.labFrameQuadrupole_XY*guy15  + atomJ.labFrameQuadrupole_XZ*guy16  + atomJ.labFrameQuadrupole_YZ*guy19)) +
            szi*(atomJ.labFrameQuadrupole_XX*guz13  + atomJ.labFrameQuadrupole_YY*guz18  + atomJ.labFrameQuadrupole_ZZ*guz20  + 2.0f*(atomJ.labFrameQuadrupole_XY*guz15  + atomJ.labFrameQuadrupole_XZ*guz16  + atomJ.labFrameQuadrupole_YZ*guz19));

    dpdz += sxk*(atomI.labFrameQuadrupole_XX*gux13 + atomI.labFrameQuadrupole_YY*gux18 + atomI.labFrameQuadrupole_ZZ*gux20 + 2.0f*(atomI.labFrameQuadrupole_XY*gux15 + atomI.labFrameQuadrupole_XZ*gux16 + atomI.labFrameQuadrupole_YZ*gux19)) +
            syk*(atomI.labFrameQuadrupole_XX*guy13 + atomI.labFrameQuadrupole_YY*guy18 + atomI.labFrameQuadrupole_ZZ*guy20 + 2.0f*(atomI.labFrameQuadrupole_XY*guy15 + atomI.labFrameQuadrupole_XZ*guy16 + atomI.labFrameQuadrupole_YZ*guy19)) +
            szk*(atomI.labFrameQuadrupole_XX*guz13 + atomI.labFrameQuadrupole_YY*guz18 + atomI.labFrameQuadrupole_ZZ*guz20 + 2.0f*(atomI.labFrameQuadrupole_XY*guz15 + atomI.labFrameQuadrupole_XZ*guz16 + atomI.labFrameQuadrupole_YZ*guz19));

#endif

#endif

#if defined F1
    float gqxx11     = xr*(12.0f*a21 + xr2*(9.0f*a22  + xr2*a23));
    float gqxx12     = yr*(2.0f*a21 + xr2*(5.0f*a22   + xr2*a23));
    float gqxx13     = zr*(2.0f*a21 + xr2*(5.0f*a22   + xr2*a23));
    float gqxx14     = xr*(2.0f*a21 + yr2*2.0f*a22    + xr2*(a22 + yr2*a23));
    float gqxx15     = xr*yr*zr*(2.0f*a22 + xr2*a23);
    float gqxx16     = xr*(2.0f*a21 + zr2*2.0f*a22  + xr2*(a22 + zr2*a23));
    float gqxx17     = yr*xr2*(3.0f*a22 + yr2*a23);
    float gqxx18     = zr*xr2*(a22 + yr2*a23);
    float gqxx19     = yr*xr2*(a22 + zr2*a23);
    float gqxx20     = zr*xr2*(3.0f*a22 + zr2*a23);

    float gqxy11     = yr*(3.0f*a21 + xr2*(6.0f*a22  + xr2*a23));
    float gqxy12     = xr*(3.0f*(a21 + yr2*a22)  + xr2*(a22 + yr2*a23));

    float gqxy13     = xr*yr*zr*(3.0f*a22 + xr2*a23);

    float gqxy14     = yr*(3.0f*(a21 + xr2*a22)  + yr2*(a22 + xr2*a23));
    float gqxy15     = zr*(a21 + (yr2 + xr2)*a22  + yr2*xr2*a23);
    float gqxy16     = yr*(a21 + (xr2 + zr2)*a22  + xr2*zr2*a23);
    float gqxy17     = xr*(3.0f*(a21 + yr2*a22)  + yr2*(3.0f*a22 + yr2*a23));
    float gqxy18     = xr*yr*zr*(3.0f*a22 + yr2*a23);
    float gqxy19     = xr*(a21 + (yr2 + zr2)*a22  + yr2*zr2*a23);
    float gqxy20     = xr*yr*zr*(3.0f*a22 + zr2*a23);
    float gqxz11     = zr*(3.0f*a21 + xr2*(6.0f*a22  + xr2*a23));

    float gqxz12     = xr*yr*zr*(3.0f*a22 + xr2*a23);

    float gqxz13     = xr*(3.0f*(a21 + zr2*a22)  + xr2*(a22 + zr2*a23));
    float gqxz14     = zr*(a21 + (xr2 + yr2)*a22  + xr2*yr2*a23);
    float gqxz15     = yr*(a21 + (xr2 + zr2)*a22  + zr2*xr2*a23);
    float gqxz16     = zr*(3.0f*(a21 + xr2*a22)  + zr2*(a22 + xr2*a23));
    float gqxz17     = xr*yr*zr*(3.0f*a22 + yr2*a23);
    float gqxz18     = xr*(a21 + (zr2 + yr2)*a22  + zr2*yr2*a23);
    float gqxz19     = xr*yr*zr*(3.0f*a22 + zr2*a23);
    float gqxz20     = xr*(3.0f*a21 + zr2*(6.0f*a22  + zr2*a23));
    float gqyy11     = xr*yr2*(3.0f*a22 + xr2*a23);
    float gqyy12     = yr*(2.0f*a21 + xr2*2.0f*a22  + yr2*(a22 + xr2*a23));
    float gqyy13     = zr*yr2*(a22 + xr2*a23);
    float gqyy14     = xr*(2.0f*a21 + yr2*(5.0f*a22  + yr2*a23));
    float gqyy15     = xr*yr*zr*(2.0f*a22 + yr2*a23);
    float gqyy16     = xr*yr2*(a22 + zr2*a23);
    float gqyy17     = yr*(12.0f*a21 + yr2*(9.0f*a22  + yr2*a23));
    float gqyy18     = zr*(2.0f*a21 + yr2*(5.0f*a22  + yr2*a23));
    float gqyy19     = yr*(2.0f*a21 + zr2*2.0f*a22  + yr2*(a22 + zr2*a23));
    float gqyy20     = zr*yr2*(3.0f*a22 + zr2*a23);
    float gqyz11     = xr*yr*zr*(3.0f*a22 + xr2*a23);
    float gqyz12     = zr*(a21 + (xr2 + yr2)*a22  + xr2*yr2*a23);
    float gqyz13     = yr*(a21 + (xr2 + zr2)*a22  + xr2*zr2*a23);
    float gqyz14     = xr*yr*zr*(3.0f*a22 + yr2*a23);
    float gqyz15     = xr*(a21 + (yr2 + zr2)*a22  + yr2*zr2*a23);
    float gqyz16     = xr*yr*zr*(3.0f*a22 + zr2*a23);
    float gqyz17     = zr*(3.0f*a21 + yr2*(6.0f*a22  + yr2*a23));
    float gqyz18     = yr*(3.0f*(a21 + zr2*a22)  + yr2*(a22 + zr2*a23));
    float gqyz19     = zr*(3.0f*(a21 + yr2*a22)  + zr2*(a22 + yr2*a23));
    float gqyz20     = yr*(3.0f*a21 + zr2*(6.0f*a22  + zr2*a23));
    float gqzz11     = xr*zr2*(3.0f*a22 + xr2*a23);
    float gqzz12     = yr*(zr2*a22 + xr2*(zr2*a23));
    float gqzz13     = zr*(2.0f*a21 + xr2*2.0f*a22  + zr2*(a22 + xr2*a23));
    float gqzz14     = xr*zr2*(a22 + yr2*a23);
    float gqzz15     = xr*yr*zr*(2.0f*a22 + zr2*a23);
    float gqzz16     = xr*(2.0f*a21 + zr2*(5.0f*a22  + zr2*a23));
    float gqzz17     = yr*zr2*(3.0f*a22 + yr2*a23);
    float gqzz18     = zr*(2.0f*a21 + yr2*2.0f*a22  + zr2*(a22 + yr2*a23));
    float gqzz19     = yr*(2.0f*a21 + zr2*(5.0f*a22  + zr2*a23));
    float gqzz20     = zr*(12.0f*a21 + zr2*(9.0f*a22  + zr2*a23));

    dedx += atomI.labFrameQuadrupole_XX*(atomJ.labFrameQuadrupole_XX*gqxx11 + atomJ.labFrameQuadrupole_YY*gqxx14 + atomJ.labFrameQuadrupole_ZZ*gqxx16 + 2.0f*(atomJ.labFrameQuadrupole_XY*gqxx12 + atomJ.labFrameQuadrupole_XZ*gqxx13 + atomJ.labFrameQuadrupole_YZ*gqxx15)) +
            atomI.labFrameQuadrupole_YY*(atomJ.labFrameQuadrupole_XX*gqyy11 + atomJ.labFrameQuadrupole_YY*gqyy14 + atomJ.labFrameQuadrupole_ZZ*gqyy16 + 2.0f*(atomJ.labFrameQuadrupole_XY*gqyy12 + atomJ.labFrameQuadrupole_XZ*gqyy13 + atomJ.labFrameQuadrupole_YZ*gqyy15)) +
            atomI.labFrameQuadrupole_ZZ*(atomJ.labFrameQuadrupole_XX*gqzz11 + atomJ.labFrameQuadrupole_YY*gqzz14 + atomJ.labFrameQuadrupole_ZZ*gqzz16 + 2.0f*(atomJ.labFrameQuadrupole_XY*gqzz12 + atomJ.labFrameQuadrupole_XZ*gqzz13 + atomJ.labFrameQuadrupole_YZ*gqzz15)) +
          2.0f*(
            atomI.labFrameQuadrupole_XY*(atomJ.labFrameQuadrupole_XX*gqxy11 + atomJ.labFrameQuadrupole_YY*gqxy14 + atomJ.labFrameQuadrupole_ZZ*gqxy16 + 2.0f*(atomJ.labFrameQuadrupole_XY*gqxy12 + atomJ.labFrameQuadrupole_XZ*gqxy13 + atomJ.labFrameQuadrupole_YZ*gqxy15)) +
            atomI.labFrameQuadrupole_XZ*(atomJ.labFrameQuadrupole_XX*gqxz11 + atomJ.labFrameQuadrupole_YY*gqxz14 + atomJ.labFrameQuadrupole_ZZ*gqxz16 + 2.0f*(atomJ.labFrameQuadrupole_XY*gqxz12 + atomJ.labFrameQuadrupole_XZ*gqxz13 + atomJ.labFrameQuadrupole_YZ*gqxz15)) +
            atomI.labFrameQuadrupole_YZ*(atomJ.labFrameQuadrupole_XX*gqyz11 + atomJ.labFrameQuadrupole_YY*gqyz14 + atomJ.labFrameQuadrupole_ZZ*gqyz16 + 2.0f*(atomJ.labFrameQuadrupole_XY*gqyz12 + atomJ.labFrameQuadrupole_XZ*gqyz13 + atomJ.labFrameQuadrupole_YZ*gqyz15))) +

            atomI.labFrameQuadrupole_XX*(atomJ.labFrameQuadrupole_XX*gqxx11 + atomJ.labFrameQuadrupole_YY*gqyy11 + atomJ.labFrameQuadrupole_ZZ*gqzz11 + 2.0f*(atomJ.labFrameQuadrupole_XY*gqxy11 + atomJ.labFrameQuadrupole_XZ*gqxz11 + atomJ.labFrameQuadrupole_YZ*gqyz11)) +
            atomI.labFrameQuadrupole_YY*(atomJ.labFrameQuadrupole_XX*gqxx14 + atomJ.labFrameQuadrupole_YY*gqyy14 + atomJ.labFrameQuadrupole_ZZ*gqzz14 + 2.0f*(atomJ.labFrameQuadrupole_XY*gqxy14 + atomJ.labFrameQuadrupole_XZ*gqxz14 + atomJ.labFrameQuadrupole_YZ*gqyz14)) +
            atomI.labFrameQuadrupole_ZZ*(atomJ.labFrameQuadrupole_XX*gqxx16 + atomJ.labFrameQuadrupole_YY*gqyy16 + atomJ.labFrameQuadrupole_ZZ*gqzz16 + 2.0f*(atomJ.labFrameQuadrupole_XY*gqxy16 + atomJ.labFrameQuadrupole_XZ*gqxz16 + atomJ.labFrameQuadrupole_YZ*gqyz16)) +

          2.0f*(
            atomI.labFrameQuadrupole_XY*(atomJ.labFrameQuadrupole_XX*gqxx12 + atomJ.labFrameQuadrupole_YY*gqyy12 + atomJ.labFrameQuadrupole_ZZ*gqzz12 + 2.0f*(atomJ.labFrameQuadrupole_XY*gqxy12 + atomJ.labFrameQuadrupole_XZ*gqxz12 + atomJ.labFrameQuadrupole_YZ*gqyz12)) +
            atomI.labFrameQuadrupole_XZ*(atomJ.labFrameQuadrupole_XX*gqxx13 + atomJ.labFrameQuadrupole_YY*gqyy13 + atomJ.labFrameQuadrupole_ZZ*gqzz13 + 2.0f*(atomJ.labFrameQuadrupole_XY*gqxy13 + atomJ.labFrameQuadrupole_XZ*gqxz13 + atomJ.labFrameQuadrupole_YZ*gqyz13)) +
            atomI.labFrameQuadrupole_YZ*(atomJ.labFrameQuadrupole_XX*gqxx15 + atomJ.labFrameQuadrupole_YY*gqyy15 + atomJ.labFrameQuadrupole_ZZ*gqzz15 + 2.0f*(atomJ.labFrameQuadrupole_XY*gqxy15 + atomJ.labFrameQuadrupole_XZ*gqxz15 + atomJ.labFrameQuadrupole_YZ*gqyz15)));

    dedy += atomI.labFrameQuadrupole_XX*(atomJ.labFrameQuadrupole_XX*gqxx12 + atomJ.labFrameQuadrupole_YY*gqxx17 + atomJ.labFrameQuadrupole_ZZ*gqxx19   + 2.0f*(atomJ.labFrameQuadrupole_XY*gqxx14 + atomJ.labFrameQuadrupole_XZ*gqxx15 + atomJ.labFrameQuadrupole_YZ*gqxx18)) +
            atomI.labFrameQuadrupole_YY*(atomJ.labFrameQuadrupole_XX*gqyy12 + atomJ.labFrameQuadrupole_YY*gqyy17 + atomJ.labFrameQuadrupole_ZZ*gqyy19   + 2.0f*(atomJ.labFrameQuadrupole_XY*gqyy14 + atomJ.labFrameQuadrupole_XZ*gqyy15 + atomJ.labFrameQuadrupole_YZ*gqyy18)) +
            atomI.labFrameQuadrupole_ZZ*(atomJ.labFrameQuadrupole_XX*gqzz12 + atomJ.labFrameQuadrupole_YY*gqzz17 + atomJ.labFrameQuadrupole_ZZ*gqzz19 + 2.0f*(atomJ.labFrameQuadrupole_XY*gqzz14 + atomJ.labFrameQuadrupole_XZ*gqzz15 + atomJ.labFrameQuadrupole_YZ*gqzz18)) +

          2.0f*(
            atomI.labFrameQuadrupole_XY*(atomJ.labFrameQuadrupole_XX*gqxy12 + atomJ.labFrameQuadrupole_YY*gqxy17 + atomJ.labFrameQuadrupole_ZZ*gqxy19 + 2.0f*(atomJ.labFrameQuadrupole_XY*gqxy14 + atomJ.labFrameQuadrupole_XZ*gqxy15 + atomJ.labFrameQuadrupole_YZ*gqxy18)) +
            atomI.labFrameQuadrupole_XZ*(atomJ.labFrameQuadrupole_XX*gqxz12 + atomJ.labFrameQuadrupole_YY*gqxz17 + atomJ.labFrameQuadrupole_ZZ*gqxz19 + 2.0f*(atomJ.labFrameQuadrupole_XY*gqxz14 + atomJ.labFrameQuadrupole_XZ*gqxz15 + atomJ.labFrameQuadrupole_YZ*gqxz18)) +
            atomI.labFrameQuadrupole_YZ*(atomJ.labFrameQuadrupole_XX*gqyz12 + atomJ.labFrameQuadrupole_YY*gqyz17 + atomJ.labFrameQuadrupole_ZZ*gqyz19 + 2.0f*(atomJ.labFrameQuadrupole_XY*gqyz14 + atomJ.labFrameQuadrupole_XZ*gqyz15 + atomJ.labFrameQuadrupole_YZ*gqyz18))) +

            atomI.labFrameQuadrupole_XX*(atomJ.labFrameQuadrupole_XX*gqxx12 + atomJ.labFrameQuadrupole_YY*gqyy12 + atomJ.labFrameQuadrupole_ZZ*gqzz12 + 2.0f*(atomJ.labFrameQuadrupole_XY*gqxy12 + atomJ.labFrameQuadrupole_XZ*gqxz12 + atomJ.labFrameQuadrupole_YZ*gqyz12)) +
            atomI.labFrameQuadrupole_YY*(atomJ.labFrameQuadrupole_XX*gqxx17 + atomJ.labFrameQuadrupole_YY*gqyy17 + atomJ.labFrameQuadrupole_ZZ*gqzz17 + 2.0f*(atomJ.labFrameQuadrupole_XY*gqxy17 + atomJ.labFrameQuadrupole_XZ*gqxz17 + atomJ.labFrameQuadrupole_YZ*gqyz17)) +
            atomI.labFrameQuadrupole_ZZ*(atomJ.labFrameQuadrupole_XX*gqxx19 + atomJ.labFrameQuadrupole_YY*gqyy19 + atomJ.labFrameQuadrupole_ZZ*gqzz19 + 2.0f*(atomJ.labFrameQuadrupole_XY*gqxy19 + atomJ.labFrameQuadrupole_XZ*gqxz19 + atomJ.labFrameQuadrupole_YZ*gqyz19)) +

          2.0f*(
            atomI.labFrameQuadrupole_XY*(atomJ.labFrameQuadrupole_XX*gqxx14 + atomJ.labFrameQuadrupole_YY*gqyy14 + atomJ.labFrameQuadrupole_ZZ*gqzz14 + 2.0f*(atomJ.labFrameQuadrupole_XY*gqxy14 + atomJ.labFrameQuadrupole_XZ*gqxz14 + atomJ.labFrameQuadrupole_YZ*gqyz14)) +
            atomI.labFrameQuadrupole_XZ*(atomJ.labFrameQuadrupole_XX*gqxx15 + atomJ.labFrameQuadrupole_YY*gqyy15 + atomJ.labFrameQuadrupole_ZZ*gqzz15 + 2.0f*(atomJ.labFrameQuadrupole_XY*gqxy15 + atomJ.labFrameQuadrupole_XZ*gqxz15 + atomJ.labFrameQuadrupole_YZ*gqyz15)) +
            atomI.labFrameQuadrupole_YZ*(atomJ.labFrameQuadrupole_XX*gqxx18 + atomJ.labFrameQuadrupole_YY*gqyy18 + atomJ.labFrameQuadrupole_ZZ*gqzz18 + 2.0f*(atomJ.labFrameQuadrupole_XY*gqxy18 + atomJ.labFrameQuadrupole_XZ*gqxz18 + atomJ.labFrameQuadrupole_YZ*gqyz18)));

    dedz += atomI.labFrameQuadrupole_XX*(atomJ.labFrameQuadrupole_XX*gqxx13 + atomJ.labFrameQuadrupole_YY*gqxx18 + atomJ.labFrameQuadrupole_ZZ*gqxx20 + 2.0f*(atomJ.labFrameQuadrupole_XY*gqxx15 + atomJ.labFrameQuadrupole_XZ*gqxx16 + atomJ.labFrameQuadrupole_YZ*gqxx19)) +
            atomI.labFrameQuadrupole_YY*(atomJ.labFrameQuadrupole_XX*gqyy13 + atomJ.labFrameQuadrupole_YY*gqyy18 + atomJ.labFrameQuadrupole_ZZ*gqyy20 + 2.0f*(atomJ.labFrameQuadrupole_XY*gqyy15 + atomJ.labFrameQuadrupole_XZ*gqyy16 + atomJ.labFrameQuadrupole_YZ*gqyy19)) +
            atomI.labFrameQuadrupole_ZZ*(atomJ.labFrameQuadrupole_XX*gqzz13 + atomJ.labFrameQuadrupole_YY*gqzz18 + atomJ.labFrameQuadrupole_ZZ*gqzz20 + 2.0f*(atomJ.labFrameQuadrupole_XY*gqzz15 + atomJ.labFrameQuadrupole_XZ*gqzz16 + atomJ.labFrameQuadrupole_YZ*gqzz19)) +

           2.0f*(
            atomI.labFrameQuadrupole_XY*(atomJ.labFrameQuadrupole_XX*gqxy13 + atomJ.labFrameQuadrupole_YY*gqxy18 + atomJ.labFrameQuadrupole_ZZ*gqxy20 + 2.0f*(atomJ.labFrameQuadrupole_XY*gqxy15 + atomJ.labFrameQuadrupole_XZ*gqxy16 + atomJ.labFrameQuadrupole_YZ*gqxy19)) +
            atomI.labFrameQuadrupole_XZ*(atomJ.labFrameQuadrupole_XX*gqxz13 + atomJ.labFrameQuadrupole_YY*gqxz18 + atomJ.labFrameQuadrupole_ZZ*gqxz20 + 2.0f*(atomJ.labFrameQuadrupole_XY*gqxz15 + atomJ.labFrameQuadrupole_XZ*gqxz16 + atomJ.labFrameQuadrupole_YZ*gqxz19)) +
            atomI.labFrameQuadrupole_YZ*(atomJ.labFrameQuadrupole_XX*gqyz13 + atomJ.labFrameQuadrupole_YY*gqyz18 + atomJ.labFrameQuadrupole_ZZ*gqyz20 + 2.0f*(atomJ.labFrameQuadrupole_XY*gqyz15 + atomJ.labFrameQuadrupole_XZ*gqyz16 + atomJ.labFrameQuadrupole_YZ*gqyz19))) +

            atomI.labFrameQuadrupole_XX*(atomJ.labFrameQuadrupole_XX*gqxx13 + atomJ.labFrameQuadrupole_YY*gqyy13 + atomJ.labFrameQuadrupole_ZZ*gqzz13 + 2.0f*(atomJ.labFrameQuadrupole_XY*gqxy13 + atomJ.labFrameQuadrupole_XZ*gqxz13 + atomJ.labFrameQuadrupole_YZ*gqyz13)) +
            atomI.labFrameQuadrupole_YY*(atomJ.labFrameQuadrupole_XX*gqxx18 + atomJ.labFrameQuadrupole_YY*gqyy18 + atomJ.labFrameQuadrupole_ZZ*gqzz18 + 2.0f*(atomJ.labFrameQuadrupole_XY*gqxy18 + atomJ.labFrameQuadrupole_XZ*gqxz18 + atomJ.labFrameQuadrupole_YZ*gqyz18)) +
            atomI.labFrameQuadrupole_ZZ*(atomJ.labFrameQuadrupole_XX*gqxx20 + atomJ.labFrameQuadrupole_YY*gqyy20 + atomJ.labFrameQuadrupole_ZZ*gqzz20 + 2.0f*(atomJ.labFrameQuadrupole_XY*gqxy20 + atomJ.labFrameQuadrupole_XZ*gqxz20 + atomJ.labFrameQuadrupole_YZ*gqyz20)) +

           2.0f*(
            atomI.labFrameQuadrupole_XY*(atomJ.labFrameQuadrupole_XX*gqxx15 + atomJ.labFrameQuadrupole_YY*gqyy15 + atomJ.labFrameQuadrupole_ZZ*gqzz15 + 2.0f*(atomJ.labFrameQuadrupole_XY*gqxy15 + atomJ.labFrameQuadrupole_XZ*gqxz15 + atomJ.labFrameQuadrupole_YZ*gqyz15)) +
            atomI.labFrameQuadrupole_XZ*(atomJ.labFrameQuadrupole_XX*gqxx16 + atomJ.labFrameQuadrupole_YY*gqyy16 + atomJ.labFrameQuadrupole_ZZ*gqzz16 + 2.0f*(atomJ.labFrameQuadrupole_XY*gqxy16 + atomJ.labFrameQuadrupole_XZ*gqxz16 + atomJ.labFrameQuadrupole_YZ*gqyz16)) +
            atomI.labFrameQuadrupole_YZ*(atomJ.labFrameQuadrupole_XX*gqxx19 + atomJ.labFrameQuadrupole_YY*gqyy19 + atomJ.labFrameQuadrupole_ZZ*gqzz19 + 2.0f*(atomJ.labFrameQuadrupole_XY*gqxy19 + atomJ.labFrameQuadrupole_XZ*gqxz19 + atomJ.labFrameQuadrupole_YZ*gqyz19)));

#endif


#if defined T1

    if ( xr != 0.0f || yr != 0.0f || zr != 0.0f )
    {

        float gux1       = xr*a10;
        float guy1       = yr*a10;
        float guz1       = zr*a10;

        float gc2        = xr*a01;
        float gc3        = yr*a01;
        float gc4        = zr*a01;
        float fid1       = atomJ.labFrameDipole[0]*gux2  + atomJ.labFrameDipole[1]*gux3  + atomJ.labFrameDipole[2]*gux4 + 0.5f*(atomJ.q*gux1 + atomJ.labFrameQuadrupole_XX*gux5 + atomJ.labFrameQuadrupole_YY*gux8 + atomJ.labFrameQuadrupole_ZZ*gux10 +
                           2.0f*(atomJ.labFrameQuadrupole_XY*gux6 + atomJ.labFrameQuadrupole_XZ*gux7 + atomJ.labFrameQuadrupole_YZ*gux9) +
                          atomJ.q*gc2 + atomJ.labFrameQuadrupole_XX*gqxx2 + atomJ.labFrameQuadrupole_YY*gqyy2 + atomJ.labFrameQuadrupole_ZZ*gqzz2 +
                           2.0f*(atomJ.labFrameQuadrupole_XY*gqxy2 + atomJ.labFrameQuadrupole_XZ*gqxz2 + atomJ.labFrameQuadrupole_YZ*gqyz2));

        float fid2       = atomJ.labFrameDipole[0]*guy2  + atomJ.labFrameDipole[1]*guy3  + atomJ.labFrameDipole[2]*guy4 + 0.5f*(atomJ.q*guy1 + atomJ.labFrameQuadrupole_XX*guy5 + atomJ.labFrameQuadrupole_YY*guy8 + atomJ.labFrameQuadrupole_ZZ*guy10 +
                           2.0f*(atomJ.labFrameQuadrupole_XY*guy6 + atomJ.labFrameQuadrupole_XZ*guy7 + atomJ.labFrameQuadrupole_YZ*guy9) +
                          atomJ.q*gc3 + atomJ.labFrameQuadrupole_XX*gqxx3 + atomJ.labFrameQuadrupole_YY*gqyy3 + atomJ.labFrameQuadrupole_ZZ*gqzz3 + 
                           2.0f*(atomJ.labFrameQuadrupole_XY*gqxy3 + atomJ.labFrameQuadrupole_XZ*gqxz3 + atomJ.labFrameQuadrupole_YZ*gqyz3));

        float fid3       = atomJ.labFrameDipole[0]*guz2  + atomJ.labFrameDipole[1]*guz3  + atomJ.labFrameDipole[2]*guz4 + 0.5f*(atomJ.q*guz1 + atomJ.labFrameQuadrupole_XX*guz5 + atomJ.labFrameQuadrupole_YY*guz8 + atomJ.labFrameQuadrupole_ZZ*guz10 +
                           2.0f*(atomJ.labFrameQuadrupole_XY*guz6 + atomJ.labFrameQuadrupole_XZ*guz7 + atomJ.labFrameQuadrupole_YZ*guz9) +
                           atomJ.q*gc4 + atomJ.labFrameQuadrupole_XX*gqxx4 + atomJ.labFrameQuadrupole_YY*gqyy4 + atomJ.labFrameQuadrupole_ZZ*gqzz4 +
                           2.0f*(atomJ.labFrameQuadrupole_XY*gqxy4 + atomJ.labFrameQuadrupole_XZ*gqxz4 + atomJ.labFrameQuadrupole_YZ*gqyz4));

        float trq1       = atomI.labFrameDipole[1]*fid3 - atomI.labFrameDipole[2]*fid2;
        float trq2       = atomI.labFrameDipole[2]*fid1 - atomI.labFrameDipole[0]*fid3;
        float trq3       = atomI.labFrameDipole[0]*fid2 - atomI.labFrameDipole[1]*fid1;

        // torque on quadrupoles due to permanent reaction field gradient

        float fidg11 =
                (atomJ.q*xr2*a20 + atomJ.labFrameDipole[0]*gqxx2 + atomJ.labFrameDipole[1]*gqxx3 + atomJ.labFrameDipole[2]*gqxx4
                       + atomJ.labFrameQuadrupole_XX*gqxx5 + atomJ.labFrameQuadrupole_YY*gqxx8 + atomJ.labFrameQuadrupole_ZZ*gqxx10
                       + 2.0f*(atomJ.labFrameQuadrupole_XY*gqxx6 + atomJ.labFrameQuadrupole_XZ*gqxx7 + atomJ.labFrameQuadrupole_YZ*gqxx9)
                       + atomJ.q*gc5 + atomJ.labFrameDipole[0]*gux5 + atomJ.labFrameDipole[1]*guy5 + atomJ.labFrameDipole[2]*guz5
                       + atomJ.labFrameQuadrupole_XX*gqxx5 + atomJ.labFrameQuadrupole_YY*gqyy5 + atomJ.labFrameQuadrupole_ZZ*gqzz5
                       + 2.0f*(atomJ.labFrameQuadrupole_XY*gqxy5 + atomJ.labFrameQuadrupole_XZ*gqxz5 + atomJ.labFrameQuadrupole_YZ*gqyz5));

        float fidg12 =
                (atomJ.q*xr*yr*a20 + atomJ.labFrameDipole[0]*gqxy2 + atomJ.labFrameDipole[1]*gqxy3 + atomJ.labFrameDipole[2]*gqxy4
                       + atomJ.labFrameQuadrupole_XX*gqxy5 + atomJ.labFrameQuadrupole_YY*gqxy8 + atomJ.labFrameQuadrupole_ZZ*gqxy10
                       + 2.0f*(atomJ.labFrameQuadrupole_XY*gqxy6 + atomJ.labFrameQuadrupole_XZ*gqxy7 + atomJ.labFrameQuadrupole_YZ*gqxy9)
                       + atomJ.q*gc6 + atomJ.labFrameDipole[0]*gux6 + atomJ.labFrameDipole[1]*guy6 + atomJ.labFrameDipole[2]*guz6
                       + atomJ.labFrameQuadrupole_XX*gqxx6 + atomJ.labFrameQuadrupole_YY*gqyy6 + atomJ.labFrameQuadrupole_ZZ*gqzz6
                       + 2.0f*(atomJ.labFrameQuadrupole_XY*gqxy6 + atomJ.labFrameQuadrupole_XZ*gqxz6 + atomJ.labFrameQuadrupole_YZ*gqyz6));

        float fidg13 =
                (atomJ.q*xr*zr*a20 + atomJ.labFrameDipole[0]*gqxz2 + atomJ.labFrameDipole[1]*gqxz3 + atomJ.labFrameDipole[2]*gqxz4
                       + atomJ.labFrameQuadrupole_XX*gqxz5 + atomJ.labFrameQuadrupole_YY*gqxz8 + atomJ.labFrameQuadrupole_ZZ*gqxz10
                       + 2.0f*(atomJ.labFrameQuadrupole_XY*gqxz6 + atomJ.labFrameQuadrupole_XZ*gqxz7 + atomJ.labFrameQuadrupole_YZ*gqxz9)
                       + atomJ.q*gc7 + atomJ.labFrameDipole[0]*gux7 + atomJ.labFrameDipole[1]*guy7 + atomJ.labFrameDipole[2]*guz7
                       + atomJ.labFrameQuadrupole_XX*gqxx7 + atomJ.labFrameQuadrupole_YY*gqyy7 + atomJ.labFrameQuadrupole_ZZ*gqzz7
                       + 2.0f*(atomJ.labFrameQuadrupole_XY*gqxy7 + atomJ.labFrameQuadrupole_XZ*gqxz7 + atomJ.labFrameQuadrupole_YZ*gqyz7));

        float fidg22 =
                (atomJ.q*yr2*a20 + atomJ.labFrameDipole[0]*gqyy2 + atomJ.labFrameDipole[1]*gqyy3 + atomJ.labFrameDipole[2]*gqyy4
                       + atomJ.labFrameQuadrupole_XX*gqyy5 + atomJ.labFrameQuadrupole_YY*gqyy8 + atomJ.labFrameQuadrupole_ZZ*gqyy10
                       + 2.0f*(atomJ.labFrameQuadrupole_XY*gqyy6 + atomJ.labFrameQuadrupole_XZ*gqyy7 + atomJ.labFrameQuadrupole_YZ*gqyy9)
                       + atomJ.q*gc8 + atomJ.labFrameDipole[0]*gux8 + atomJ.labFrameDipole[1]*guy8 + atomJ.labFrameDipole[2]*guz8
                       + atomJ.labFrameQuadrupole_XX*gqxx8 + atomJ.labFrameQuadrupole_YY*gqyy8 + atomJ.labFrameQuadrupole_ZZ*gqzz8
                       + 2.0f*(atomJ.labFrameQuadrupole_XY*gqxy8 + atomJ.labFrameQuadrupole_XZ*gqxz8 + atomJ.labFrameQuadrupole_YZ*gqyz8));

        float fidg23 =
                (atomJ.q*yr*zr*a20 + atomJ.labFrameDipole[0]*gqyz2 + atomJ.labFrameDipole[1]*gqyz3 + atomJ.labFrameDipole[2]*gqyz4
                       + atomJ.labFrameQuadrupole_XX*gqyz5 + atomJ.labFrameQuadrupole_YY*gqyz8 + atomJ.labFrameQuadrupole_ZZ*gqyz10
                       + 2.0f*(atomJ.labFrameQuadrupole_XY*gqyz6 + atomJ.labFrameQuadrupole_XZ*gqyz7 + atomJ.labFrameQuadrupole_YZ*gqyz9)
                       + atomJ.q*gc9 + atomJ.labFrameDipole[0]*gux9 + atomJ.labFrameDipole[1]*guy9 + atomJ.labFrameDipole[2]*guz9
                       + atomJ.labFrameQuadrupole_XX*gqxx9 + atomJ.labFrameQuadrupole_YY*gqyy9 + atomJ.labFrameQuadrupole_ZZ*gqzz9
                       + 2.0f*(atomJ.labFrameQuadrupole_XY*gqxy9 + atomJ.labFrameQuadrupole_XZ*gqxz9 + atomJ.labFrameQuadrupole_YZ*gqyz9));

        float fidg33 =
                (atomJ.q*zr2*a20 + atomJ.labFrameDipole[0]*gqzz2 + atomJ.labFrameDipole[1]*gqzz3 + atomJ.labFrameDipole[2]*gqzz4
                       + atomJ.labFrameQuadrupole_XX*gqzz5 + atomJ.labFrameQuadrupole_YY*gqzz8 + atomJ.labFrameQuadrupole_ZZ*gqzz10
                       + 2.0f*(atomJ.labFrameQuadrupole_XY*gqzz6 + atomJ.labFrameQuadrupole_XZ*gqzz7 + atomJ.labFrameQuadrupole_YZ*gqzz9)
                       + atomJ.q*gc10 + atomJ.labFrameDipole[0]*gux10 + atomJ.labFrameDipole[1]*guy10 + atomJ.labFrameDipole[2]*guz10
                       + atomJ.labFrameQuadrupole_XX*gqxx10 + atomJ.labFrameQuadrupole_YY*gqyy10 + atomJ.labFrameQuadrupole_ZZ*gqzz10
                    + 2.0f*(atomJ.labFrameQuadrupole_XY*gqxy10 + atomJ.labFrameQuadrupole_XZ*gqxz10 + atomJ.labFrameQuadrupole_YZ*gqyz10));

        trq1   -= (atomI.labFrameQuadrupole_XY*fidg13 + atomI.labFrameQuadrupole_YY*fidg23 + atomI.labFrameQuadrupole_YZ*fidg33 -atomI.labFrameQuadrupole_XZ*fidg12-atomI.labFrameQuadrupole_YZ*fidg22-atomI.labFrameQuadrupole_ZZ*fidg23);
        trq2   -= (atomI.labFrameQuadrupole_XZ*fidg11 + atomI.labFrameQuadrupole_YZ*fidg12 + atomI.labFrameQuadrupole_ZZ*fidg13 -atomI.labFrameQuadrupole_XX*fidg13-atomI.labFrameQuadrupole_XY*fidg23-atomI.labFrameQuadrupole_XZ*fidg33);
        trq3   -= (atomI.labFrameQuadrupole_XX*fidg12 + atomI.labFrameQuadrupole_XY*fidg22 + atomI.labFrameQuadrupole_XZ*fidg23 -atomI.labFrameQuadrupole_XY*fidg11-atomI.labFrameQuadrupole_YY*fidg12-atomI.labFrameQuadrupole_YZ*fidg13);

#ifdef INCLUDE_TORQUE
        atomI.torque[0]         += trq1;
        atomI.torque[1]         += trq2;
        atomI.torque[2]         += trq3;
#else
        torque[0]                = trq1;
        torque[1]                = trq2;
        torque[2]                = trq3;
#endif

#ifndef INCLUDE_TORQUE
    } else {
        torque[0]                = 0.0f;
        torque[1]                = 0.0f;
        torque[2]                = 0.0f;
#endif
    }
#endif

#if defined B2 
    dsumdrB2               *= 0.5f;
    atomI.dBornRadiusPolar += atomJ.bornRadius*dsumdrB2;
    atomJ.dBornRadiusPolar += atomI.bornRadius*dsumdrB2;
#endif

#if defined T2

    // torque due to induced reaction field gradient on quadrupoles;

    float fidg11              = sxk*gqxx2 + syk*gqxx3 + szk*gqxx4 + sxk*gux5 + syk*guy5 + szk*guz5;
    float fidg12              = sxk*gqxy2 + syk*gqxy3 + szk*gqxy4 + sxk*gux6 + syk*guy6 + szk*guz6;
    float fidg13              = sxk*gqxz2 + syk*gqxz3 + szk*gqxz4 + sxk*gux7 + syk*guy7 + szk*guz7;
    float fidg22              = sxk*gqyy2 + syk*gqyy3 + szk*gqyy4 + sxk*gux8 + syk*guy8 + szk*guz8;
    float fidg23              = sxk*gqyz2 + syk*gqyz3 + szk*gqyz4 + sxk*gux9 + syk*guy9 + szk*guz9;
    float fidg33              = sxk*gqzz2 + syk*gqzz3 + szk*gqzz4 + sxk*gux10 + syk*guy10 + szk*guz10;

    trqi1                    -=  atomI.labFrameQuadrupole_XY*fidg13 + atomI.labFrameQuadrupole_YY*fidg23 + atomI.labFrameQuadrupole_YZ*fidg33
                                -atomI.labFrameQuadrupole_XZ*fidg12 - atomI.labFrameQuadrupole_YZ*fidg22 - atomI.labFrameQuadrupole_ZZ*fidg23;

    trqi2                    -=  atomI.labFrameQuadrupole_XZ*fidg11 + atomI.labFrameQuadrupole_YZ*fidg12 + atomI.labFrameQuadrupole_ZZ*fidg13
                                -atomI.labFrameQuadrupole_XX*fidg13 - atomI.labFrameQuadrupole_XY*fidg23 - atomI.labFrameQuadrupole_XZ*fidg33;

    trqi3                    -=  atomI.labFrameQuadrupole_XX*fidg12 + atomI.labFrameQuadrupole_XY*fidg22 + atomI.labFrameQuadrupole_XZ*fidg23
                                -atomI.labFrameQuadrupole_XY*fidg11 - atomI.labFrameQuadrupole_YY*fidg12 - atomI.labFrameQuadrupole_YZ*fidg13;


#ifdef INCLUDE_TORQUE
    atomI.torque[0]          += 0.5f*trqi1;
    atomI.torque[1]          += 0.5f*trqi2;
    atomI.torque[2]          += 0.5f*trqi3;
#else
    torque[0]                += 0.5f*trqi1;
    torque[1]                += 0.5f*trqi2;
    torque[2]                += 0.5f*trqi3;
#endif

#endif

#if defined F1
    
    *outputEnergy   = energy;

    if( (xr != 0.0f || yr != 0.0f || zr != 0.0f) ){
        force[0]          = dedx;
        force[1]          = dedy;
        force[2]          = dedz;
    } else {
        //*outputEnergy = force[0] = force[1] = force[2] = 0.0f;
        force[0] = force[1] = force[2] = 0.0f;
    }

#endif

#if defined F2

    *outputEnergy  += 0.5f*energy;

    dpdx           *= 0.5f;
    dpdy           *= 0.5f;
    dpdz           *= 0.5f;

    if( (xr != 0.0f || yr != 0.0f || zr != 0.0f) ){
        force[0]         += dpdx;
        force[1]         += dpdy;
        force[2]         += dpdz;
    }
#endif


}
