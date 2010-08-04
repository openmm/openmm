//-----------------------------------------------------------------------------------------

//-----------------------------------------------------------------------------------------

#include "amoebaGpuTypes.h"
#include "amoebaCudaKernels.h"
#include "kCalculateAmoebaCudaUtilities.h"
#include "kCalculateAmoebaCudaKirkwoodParticle.h"
extern void kCalculateObcGbsaForces2(gpuContext gpu);

//#define AMOEBA_DEBUG
#undef AMOEBA_DEBUG

static __constant__ cudaGmxSimulation cSim;
static __constant__ cudaAmoebaGmxSimulation cAmoebaSim;

void SetCalculateAmoebaKirkwoodSim(amoebaGpuContext amoebaGpu)
{
    cudaError_t status;
    gpuContext gpu = amoebaGpu->gpuContext;
    status         = cudaMemcpyToSymbol(cSim, &gpu->sim, sizeof(cudaGmxSimulation));
    RTERROR(status, "SetCalculateAmoebaKirkwoodSim: cudaMemcpyToSymbol: SetSim copy to cSim failed");
    status         = cudaMemcpyToSymbol(cAmoebaSim, &amoebaGpu->amoebaSim, sizeof(cudaAmoebaGmxSimulation));
    RTERROR(status, "SetCalculateAmoebaKirkwoodSim: cudaMemcpyToSymbol: SetSim copy to cAmoebaSim failed");
}

void GetCalculateAmoebaKirkwoodSim(amoebaGpuContext amoebaGpu)
{
    cudaError_t status;
    gpuContext gpu = amoebaGpu->gpuContext;
    status = cudaMemcpyFromSymbol(&gpu->sim, cSim, sizeof(cudaGmxSimulation));
    RTERROR(status, "GetCalculateAmoebaKirkwoodSim: cudaMemcpyFromSymbol: SetSim copy from cSim failed");
    status = cudaMemcpyFromSymbol(&amoebaGpu->amoebaSim, cAmoebaSim, sizeof(cudaAmoebaGmxSimulation));
    RTERROR(status, "GetCalculateAmoebaKirkwoodSim: cudaMemcpyFromSymbol: SetSim copy from cAmoebaSim failed");
}

__device__ void loadKirkwoodShared( struct KirkwoodParticle* sA, unsigned int atomI,
                                    float4* atomCoord, float* labDipole, float* labQuadrupole,
                                    float* inducedDipole, float* inducedDipolePolar, float* bornRadii )
{
    // coordinates & charge

    sA->x                        = atomCoord[atomI].x;
    sA->y                        = atomCoord[atomI].y;
    sA->z                        = atomCoord[atomI].z;
    sA->q                        = atomCoord[atomI].w;

    // lab dipole

    sA->labFrameDipole[0]        = labDipole[atomI*3];
    sA->labFrameDipole[1]        = labDipole[atomI*3+1];
    sA->labFrameDipole[2]        = labDipole[atomI*3+2];

    // lab quadrupole

    sA->labFrameQuadrupole_XX    = labQuadrupole[atomI*9];
    sA->labFrameQuadrupole_XY    = labQuadrupole[atomI*9+1];
    sA->labFrameQuadrupole_XZ    = labQuadrupole[atomI*9+2];
    sA->labFrameQuadrupole_YY    = labQuadrupole[atomI*9+4];
    sA->labFrameQuadrupole_YZ    = labQuadrupole[atomI*9+5];
    sA->labFrameQuadrupole_ZZ    = labQuadrupole[atomI*9+8];

    // induced dipole

    sA->inducedDipole[0]         = inducedDipole[atomI*3];
    sA->inducedDipole[1]         = inducedDipole[atomI*3+1];
    sA->inducedDipole[2]         = inducedDipole[atomI*3+2];

    // induced dipole polar

    sA->inducedDipoleP[0]        = inducedDipolePolar[atomI*3];
    sA->inducedDipoleP[1]        = inducedDipolePolar[atomI*3+1];
    sA->inducedDipoleP[2]        = inducedDipolePolar[atomI*3+2];

    sA->bornRadius               = bornRadii[atomI];

}

// load struct and arrays w/ shared data in sA

__device__ void loadKirkwoodData( struct KirkwoodParticle* sA,
                                  float4* jCoord, float* jDipole, float* jQuadrupole,
                                  float* jInducedDipole, float* jInducedDipolePolar, float* jBornRadius )
{

    // load coords, charge, ...

    jCoord->x               = sA->x;
    jCoord->y               = sA->y;
    jCoord->z               = sA->z;
    jCoord->w               = sA->q;

    jDipole[0]              = sA->labFrameDipole[0];
    jDipole[1]              = sA->labFrameDipole[1];
    jDipole[2]              = sA->labFrameDipole[2];

    jQuadrupole[0]          = sA->labFrameQuadrupole_XX;
    jQuadrupole[1]          = sA->labFrameQuadrupole_XY;
    jQuadrupole[2]          = sA->labFrameQuadrupole_XZ;

    jQuadrupole[3]          = sA->labFrameQuadrupole_XY;
    jQuadrupole[4]          = sA->labFrameQuadrupole_YY;
    jQuadrupole[5]          = sA->labFrameQuadrupole_YZ;

    jQuadrupole[6]          = sA->labFrameQuadrupole_XZ;
    jQuadrupole[7]          = sA->labFrameQuadrupole_YZ;
    jQuadrupole[8]          = sA->labFrameQuadrupole_ZZ;

    jInducedDipole[0]       = sA->inducedDipole[0];
    jInducedDipole[1]       = sA->inducedDipole[1];
    jInducedDipole[2]       = sA->inducedDipole[2];

    jInducedDipolePolar[0]  = sA->inducedDipoleP[0];
    jInducedDipolePolar[1]  = sA->inducedDipoleP[1];
    jInducedDipolePolar[2]  = sA->inducedDipoleP[2];

   *jBornRadius             = sA->bornRadius;

}

__device__ void calculateKirkwoodPairIxn_kernel( unsigned int sameAtom,
                                                 float4 atomCoordinatesI,       float4 atomCoordinatesJ,
                                                 float* labFrameDipoleI,        float* labFrameDipoleJ,
                                                 float* labFrameQuadrupoleI,    float* labFrameQuadrupoleJ,
                                                 float* inducedDipoleI,         float* inducedDipoleJ,
                                                 float* inducedDipolePolarI,    float* inducedDipolePolarJ,
                                                 float  bornRadiusI,            float bornRadiusJ,
                                                 float*  outputForce,           float outputTorque[2][3],
                                                 float*  outputBorn,            float*  outputBornPolar,
                                                 float* outputEnergy
#ifdef AMOEBA_DEBUG
                                                 , float4* debugArray
#endif

 ){

    float e,ei;
    float xi,yi,zi;
    float xr,yr,zr;
    float xr2,yr2,zr2;
    float ci,ck;
    float uxi,uyi,uzi;
    float uxk,uyk,uzk;
    float qxxi,qxyi,qxzi;
    float qyyi,qyzi,qzzi;
    float qxxk,qxyk,qxzk;
    float qyyk,qyzk,qzzk;
    float dxi,dyi,dzi;
    float dxk,dyk,dzk;
    float pxi,pyi,pzi;
    float pxk,pyk,pzk;
    float sxi,syi,szi;
    float sxk,syk,szk;
    float r2,rb2;
    float dedx,dedy,dedz;
    float drbi;
    float drbk;
    float dpdx,dpdy,dpdz;
    float dpbi;
    float dpbk;
    float fc,fd,fq;
    float rbi,rbk;
    float expterm;
    float gf,gf2,gf3,gf5;
    float gf7,gf9,gf11;
    float expc,dexpc;
    float expc1,expcdexpc;
    float expcr,dexpcr;
    float dgfdr;
    float esym,ewi,ewk;
    float desymdx,dewidx,dewkdx;
    float desymdy,dewidy,dewkdy;
    float desymdz,dewidz,dewkdz;
    float dsumdr,desymdr;
    float dewidr,dewkdr;
    float dsymdr;
    float esymi,ewii,ewki;
    float dpsymdx,dpwidx,dpwkdx;
    float dpsymdy,dpwidy,dpwkdy;
    float dpsymdz,dpwidz,dpwkdz;
    float dwipdr,dwkpdr;
    float duvdr;

    float gkc;

    // set the bulk dielectric constant to the water value

    gkc          = cAmoebaSim.gkc;

    fc           = cAmoebaSim.electric * cAmoebaSim.fc;
    fd           = cAmoebaSim.electric * cAmoebaSim.fd;
    fq           = cAmoebaSim.electric * cAmoebaSim.fq;

    xi           = atomCoordinatesI.x;
    yi           = atomCoordinatesI.y;
    zi           = atomCoordinatesI.z;
    ci           = atomCoordinatesI.w;

    uxi          = labFrameDipoleI[0];
    uyi          = labFrameDipoleI[1];
    uzi          = labFrameDipoleI[2];

    dxi          = inducedDipoleI[0];
    dyi          = inducedDipoleI[1];
    dzi          = inducedDipoleI[2];

    pxi          = inducedDipolePolarI[0];
    pyi          = inducedDipolePolarI[1];
    pzi          = inducedDipolePolarI[2];

    qxxi         = labFrameQuadrupoleI[0];
    qxyi         = labFrameQuadrupoleI[1];
    qxzi         = labFrameQuadrupoleI[2];
    qyyi         = labFrameQuadrupoleI[4];
    qyzi         = labFrameQuadrupoleI[5];
    qzzi         = labFrameQuadrupoleI[8];

    sxi          = dxi + pxi;
    syi          = dyi + pyi;
    szi          = dzi + pzi;

    rbi          = bornRadiusI;

    // decide whether to compute the current interaction;

    xr           = atomCoordinatesJ.x - xi;
    yr           = atomCoordinatesJ.y - yi;
    zr           = atomCoordinatesJ.z - zi;
    ck           = atomCoordinatesJ.w;

    xr2          = xr*xr;
    yr2          = yr*yr;
    zr2          = zr*zr;
    r2           = xr2 + yr2 + zr2;

    if( r2 > cAmoebaSim.scalingDistanceCutoff ){
    }
    rbk          = bornRadiusJ;

    uxk          = labFrameDipoleJ[0];
    uyk          = labFrameDipoleJ[1];
    uzk          = labFrameDipoleJ[2];

    dxk          = inducedDipoleJ[0];
    dyk          = inducedDipoleJ[1];
    dzk          = inducedDipoleJ[2];

    pxk          = inducedDipolePolarJ[0];
    pyk          = inducedDipolePolarJ[1];
    pzk          = inducedDipolePolarJ[2];

    qxxk         = labFrameQuadrupoleJ[0];
    qxyk         = labFrameQuadrupoleJ[1];
    qxzk         = labFrameQuadrupoleJ[2];
    qyyk         = labFrameQuadrupoleJ[4];
    qyzk         = labFrameQuadrupoleJ[5];
    qzzk         = labFrameQuadrupoleJ[8];

    sxk          = dxk + pxk;
    syk          = dyk + pyk;
    szk          = dzk + pzk;
    rb2          = rbi * rbk;
    expterm      = expf(-r2/(gkc*rb2));
    expc         = expterm / gkc;
    expcr        = r2*expterm / (gkc*gkc*rb2*rb2);
    dexpc        = -2.0f / (gkc*rb2);
    dexpcr       = 2.0f / (gkc*rb2*rb2);
    dgfdr        = 0.5f * expterm * (1.0f+r2/(rb2*gkc));
    gf2          = 1.0f / (r2+rb2*expterm);
    gf           = sqrt(gf2);
    gf3          = gf2 * gf;
    gf5          = gf3 * gf2;
    gf7          = gf5 * gf2;
    gf9          = gf7 * gf2;
    gf11         = gf9 * gf2;

    // reaction potential auxiliary terms;

    float a00      =            gf;
    float a10      =          -gf3;
    float a20      =    3.0f * gf5;
    float a30      =  -15.0f * gf7;
    float a40      =  105.0f * gf9;
    float a50      = -945.0f * gf11;

    // Born radii derivatives of reaction potential auxiliary terms;

    float b00      = dgfdr * a10;
    float b10      = dgfdr * a20;
    float b20      = dgfdr * a30;
    float b30      = dgfdr * a40;
    float b40      = dgfdr * a50;

    // reaction potential gradient auxiliary terms;

    expc1        = 1.0f - expc;
    float a01      = expc1 * a10;
    float a11      = expc1 * a20;
    float a21      = expc1 * a30;
    float a31      = expc1 * a40;
    float a41      = expc1 * a50;

    // Born radii derivs of reaction potential gradient auxiliary terms;

    float b01      = b10 - expcr*a10 - expc*b10;
    float b11      = b20 - expcr*a20 - expc*b20;
    float b21      = b30 - expcr*a30 - expc*b30;
    float b31      = b40 - expcr*a40 - expc*b40;

    // 2nd reaction potential gradient auxiliary terms;

    expcdexpc    = -expc * dexpc;
    float a02      = expc1*a11 + expcdexpc*a10;
    float a12      = expc1*a21 + expcdexpc*a20;
    float a22      = expc1*a31 + expcdexpc*a30;
    float a32      = expc1*a41 + expcdexpc*a40;

    // Born radii derivatives of the 2nd reaction potential
    // gradient auxiliary terms

     float b02     = b11 - (expcr*(a11 + dexpc*a10)
                           + expc*(b11 + dexpcr*a10+dexpc*b10));

     float b12     = b21 - (expcr*(a21 + dexpc*a20)
                           +   expc*(b21 + dexpcr*a20+dexpc*b20));

     float b22     = b31 - (expcr*(a31 + dexpc*a30)
                           +   expc*(b31 + dexpcr*a30+dexpc*b30));

    // 3rd reaction potential gradient auxiliary terms;

    expcdexpc    = 2.0f * expcdexpc;
    float a03      = expc1*a12 + expcdexpc*a11;
    float a13      = expc1*a22 + expcdexpc*a21;
    float a23      = expc1*a32 + expcdexpc*a31;

    //expcdexpc    = -expc * powf( dexpc, 2.0f );
    expcdexpc    = -expc*dexpc*dexpc;
    a03      = a03 + expcdexpc*a10;
    a13      = a13 + expcdexpc*a20;
    a23      = a23 + expcdexpc*a30;

    // multiply the auxillary terms by their dieletric functions;

    a00      = fc * a00;
    a01      = fc * a01;
    a02      = fc * a02;
    a03      = fc * a03;
    b00      = fc * b00;
    b01      = fc * b01;
    b02      = fc * b02;
    a10      = fd * a10;
    a11      = fd * a11;
    a12      = fd * a12;
    a13      = fd * a13;
    b10      = fd * b10;
    b11      = fd * b11;
    b12      = fd * b12;
    a20      = fq * a20;
    a21      = fq * a21;
    a22      = fq * a22;
    a23      = fq * a23;
    b20      = fq * b20;
    b21      = fq * b21;
    b22      = fq * b22;

    // unweighted reaction potential tensor;

    float gc1        = a00;
    float gux1       = xr * a10;
    float guy1       = yr * a10;
    float guz1       = zr * a10;
    float gqxx1      = xr2 * a20;
    float gqyy1      = yr2 * a20;
    float gqzz1      = zr2 * a20;
    float gqxy1      = xr * yr * a20;
    float gqxz1      = xr * zr * a20;
    float gqyz1      = yr * zr * a20;

    // Born radii derivs of unweighted reaction potential tensor;

    float gc21       = b00;
    float gux21      = xr * b10;
    float guy21      = yr * b10;
    float guz21      = zr * b10;
    float gqxx21     = xr2 * b20;
    float gqyy21     = yr2 * b20;
    float gqzz21     = zr2 * b20;
    float gqxy21     = xr * yr * b20;
    float gqxz21     = xr * zr * b20;
    float gqyz21     = yr * zr * b20;

    // unweighted reaction potential gradient tensor;

    float gc2        = xr * a01;
    float gc3        = yr * a01;
    float gc4        = zr * a01;
    float gux2       = a10 + xr2*a11;
    float gux3       = xr * yr * a11;
    float gux4       = xr * zr * a11;

    float guy2       = gux3;
    float guy3       = a10 + yr2*a11;
    float guy4       = yr * zr * a11;
    float guz2       = gux4;
    float guz3       = guy4;
    float guz4       = a10 + zr2*a11;
    float gqxx2      = xr * (2.0f*a20+xr2*a21);
    float gqxx3      = yr * xr2 * a21;
    float gqxx4      = zr * xr2 * a21;
    float gqyy2      = xr * yr2 * a21;
    float gqyy3      = yr * (2.0f*a20+yr2*a21);
    float gqyy4      = zr * yr2 * a21;
    float gqzz2      = xr * zr2 * a21;
    float gqzz3      = yr * zr2 * a21;
    float gqzz4      = zr * (2.0f*a20+zr2*a21);
    float gqxy2      = yr * (a20+xr2*a21);
    float gqxy3      = xr * (a20+yr2*a21);
    float gqxy4      = zr * xr * yr * a21;
    float gqxz2      = zr * (a20+xr2*a21);
    float gqxz3      = gqxy4;
    float gqxz4      = xr * (a20+zr2*a21);
    float gqyz2      = gqxy4;
    float gqyz3      = zr * (a20+yr2*a21);
    float gqyz4      = yr * (a20+zr2*a21);

    // Born derivs of the unweighted reaction potential gradient tensor;

    float gc22       = xr * b01;
    float gc23       = yr * b01;
    float gc24       = zr * b01;
    float gux22      = b10 + xr2*b11;
    float gux23      = xr * yr * b11;
    float gux24      = xr * zr * b11;
    float guy22      = gux23;
    float guy23      = b10 + yr2*b11;
    float guy24      = yr * zr * b11;
    float guz22      = gux24;
    float guz23      = guy24;
    float guz24      = b10 + zr2*b11;
    float gqxx22     = xr * (2.0f*b20+xr2*b21);
    float gqxx23     = yr * xr2 * b21;
    float gqxx24     = zr * xr2 * b21;
    float gqyy22     = xr * yr2 * b21;
    float gqyy23     = yr * (2.0f*b20+yr2*b21);
    float gqyy24     = zr * yr2 * b21;
    float gqzz22     = xr * zr2 * b21;
    float gqzz23     = yr * zr2 * b21;
    float gqzz24     = zr * (2.0f*b20 + zr2*b21);
    float gqxy22     = yr * (b20+xr2*b21);
    float gqxy23     = xr * (b20+yr2*b21);
    float gqxy24     = zr * xr * yr * b21;
    float gqxz22     = zr * (b20+xr2*b21);
    float gqxz23     = gqxy24;
    float gqxz24     = xr * (b20+zr2*b21);
    float gqyz22     = gqxy24;
    float gqyz23     = zr * (b20+yr2*b21);
    float gqyz24     = yr * (b20+zr2*b21);

    // unweighted 2nd reaction potential gradient tensor;

    float gc5        = a01 + xr2*a02;
    float gc6        = xr * yr * a02;
    float gc7        = xr * zr * a02;
    float gc8        = a01 + yr2*a02;
    float gc9        = yr * zr * a02;
    float gc10       = a01 + zr2*a02;
    float gux5       = xr * (3.0f*a11+xr2*a12);
    float gux6       = yr * (a11+xr2*a12);
    float gux7       = zr * (a11+xr2*a12);
    float gux8       = xr * (a11+yr2*a12);
    float gux9       = zr * xr * yr * a12;
    float gux10      = xr * (a11+zr2*a12);
    float guy5       = yr * (a11+xr2*a12);
    float guy6       = xr * (a11+yr2*a12);
    float guy7       = gux9;
    float guy8       = yr * (3.0f*a11+yr2*a12);
    float guy9       = zr * (a11+yr2*a12);
    float guy10      = yr * (a11+zr2*a12);
    float guz5       = zr * (a11+xr2*a12);
    float guz6       = gux9;
    float guz7       = xr * (a11+zr2*a12);
    float guz8       = zr * (a11+yr2*a12);
    float guz9       = yr * (a11+zr2*a12);
    float guz10      = zr * (3.0f*a11+zr2*a12);
    float gqxx5      = 2.0f*a20 + xr2*(5.0f*a21+xr2*a22);
    float gqxx6      = yr * xr * (2.0f*a21+xr2*a22);
    float gqxx7      = zr * xr * (2.0f*a21+xr2*a22);
    float gqxx8      = xr2 * (a21+yr2*a22);
    float gqxx9      = zr * yr * xr2 * a22;
    float gqxx10     = xr2 * (a21+zr2*a22);
    float gqyy5      = yr2 * (a21+xr2*a22);
    float gqyy6      = xr * yr * (2.0f*a21+yr2*a22);
    float gqyy7      = xr * zr * yr2 * a22;
    float gqyy8      = 2.0f*a20 + yr2*(5.0f*a21+yr2*a22);
    float gqyy9      = yr * zr * (2.0f*a21+yr2*a22);
    float gqyy10     = yr2 * (a21+zr2*a22);
    float gqzz5      = zr2 * (a21+xr2*a22);
    float gqzz6      = xr * yr * zr2 * a22;
    float gqzz7      = xr * zr * (2.0f*a21+zr2*a22);
    float gqzz8      = zr2 * (a21+yr2*a22);
    float gqzz9      = yr * zr * (2.0f*a21+zr2*a22);
    float gqzz10     = 2.0f*a20 + zr2*(5.0f*a21+zr2*a22);
    float gqxy5      = xr * yr * (3.0f*a21+xr2*a22);
    float gqxy6      = a20 + (xr2+yr2)*a21 + xr2*yr2*a22;
    float gqxy7      = zr * yr * (a21+xr2*a22);
    float gqxy8      = xr * yr * (3.0f*a21+yr2*a22);
    float gqxy9      = zr * xr * (a21+yr2*a22);
    float gqxy10     = xr * yr * (a21+zr2*a22);
    float gqxz5      = xr * zr * (3.0f*a21+xr2*a22);
    float gqxz6      = yr * zr * (a21+xr2*a22);
    float gqxz7      = a20 + (xr2+zr2)*a21 + xr2*zr2*a22;
    float gqxz8      = xr * zr * (a21+yr2*a22);
    float gqxz9      = xr * yr * (a21+zr2*a22);
    float gqxz10     = xr * zr * (3.0f*a21+zr2*a22);
    float gqyz5      = zr * yr * (a21+xr2*a22);
    float gqyz6      = xr * zr * (a21+yr2*a22);
    float gqyz7      = xr * yr * (a21+zr2*a22);
    float gqyz8      = yr * zr * (3.0f*a21+yr2*a22);
    float gqyz9      = a20 + (yr2+zr2)*a21 + yr2*zr2*a22;
    float gqyz10     = yr * zr * (3.0f*a21+zr2*a22);

    // Born radii derivatives of the unweighted 2nd reaction;
    // potential gradient tensor;

    float gc25       = b01 + xr2*b02;
    float gc26       = xr * yr * b02;
    float gc27       = xr * zr * b02;
    float gc28       = b01 + yr2*b02;
    float gc29       = yr * zr * b02;
    float gc30       = b01 + zr2*b02;
    float gux25      = xr * (3.0f*b11+xr2*b12);
    float gux26      = yr * (b11+xr2*b12);
    float gux27      = zr * (b11+xr2*b12);
    float gux28      = xr * (b11+yr2*b12);
    float gux29      = zr * xr * yr * b12;
    float gux30      = xr * (b11+zr2*b12);
    float guy25      = yr * (b11+xr2*b12);
    float guy26      = xr * (b11+yr2*b12);
    float guy27      = gux29;
    float guy28      = yr * (3.0f*b11+yr2*b12);
    float guy29      = zr * (b11+yr2*b12);
    float guy30      = yr * (b11+zr2*b12);
    float guz25      = zr * (b11+xr2*b12);
    float guz26      = gux29;
    float guz27      = xr * (b11+zr2*b12);
    float guz28      = zr * (b11+yr2*b12);
    float guz29      = yr * (b11+zr2*b12);
    float guz30      = zr * (3.0f*b11+zr2*b12);
    float gqxx25     = 2.0f*b20 + xr2*(5.0f*b21+xr2*b22);
    float gqxx26     = yr * xr * (2.0f*b21+xr2*b22);
    float gqxx27     = zr * xr * (2.0f*b21+xr2*b22);
    float gqxx28     = xr2 * (b21+yr2*b22);
    float gqxx29     = zr * yr * xr2 * b22;
    float gqxx30     = xr2 * (b21+zr2*b22);
    float gqyy25     = yr2 * (b21+xr2*b22);
    float gqyy26     = xr * yr * (2.0f*b21+yr2*b22);
    float gqyy27     = xr * zr * yr2 * b22;
    float gqyy28     = 2.0f*b20 + yr2*(5.0f*b21+yr2*b22);
    float gqyy29     = yr * zr * (2.0f*b21+yr2*b22);
    float gqyy30     = yr2 * (b21+zr2*b22);
    float gqzz25     = zr2 * (b21+xr2*b22);
    float gqzz26     = xr * yr * zr2 * b22;
    float gqzz27     = xr * zr * (2.0f*b21+zr2*b22);
    float gqzz28     = zr2 * (b21+yr2*b22);
    float gqzz29     = yr * zr * (2.0f*b21+zr2*b22);
    float gqzz30     = 2.0f*b20 + zr2*(5.0f*b21+zr2*b22);
    float gqxy25     = xr * yr * (3.0f*b21 + xr2*b22);
    float gqxy26     = b20 + (xr2+yr2)*b21 + xr2*yr2*b22;
    float gqxy27     = zr * yr * (b21+xr2*b22);
    float gqxy28     = xr * yr * (3.0f*b21+yr2*b22);
    float gqxy29     = zr * xr * (b21+yr2*b22);
    float gqxy30     = xr * yr * (b21+zr2*b22);
    float gqxz25     = xr * zr * (3.0f*b21+xr2*b22);
    float gqxz26     = yr * zr * (b21+xr2*b22);
    float gqxz27     = b20 + (xr2+zr2)*b21 + xr2*zr2*b22;
    float gqxz28     = xr * zr * (b21+yr2*b22);
    float gqxz29     = xr * yr * (b21+zr2*b22);
    float gqxz30     = xr * zr * (3.0f*b21+zr2*b22);
    float gqyz25     = zr * yr * (b21+xr2*b22);
    float gqyz26     = xr * zr * (b21+yr2*b22);
    float gqyz27     = xr * yr * (b21+zr2*b22);
    float gqyz28     = yr * zr * (3.0f*b21+yr2*b22);
    float gqyz29     = b20 + (yr2+zr2)*b21 + yr2*zr2*b22;
    float gqyz30     = yr * zr * (3.0f*b21+zr2*b22);

    // unweighted 3rd reaction potential gradient tensor;

    float gc11       = xr * (3.0f*a02+xr2*a03);
    float gc12       = yr * (a02+xr2*a03);
    float gc13       = zr * (a02+xr2*a03);
    float gc14       = xr * (a02+yr2*a03);
    float gc15       = xr * yr * zr * a03;
    float gc16       = xr * (a02+zr2*a03);
    float gc17       = yr * (3.0f*a02+yr2*a03);
    float gc18       = zr * (a02+yr2*a03);
    float gc19       = yr * (a02+zr2*a03);
    float gc20       = zr * (3.0f*a02+zr2*a03);
    float gux11      = 3.0f*a11 + xr2*(6.0f*a12+xr2*a13);
    float gux12      = xr * yr * (3.0f*a12+xr2*a13);
    float gux13      = xr * zr * (3.0f*a12+xr2*a13);
    float gux14      = a11 + (xr2+yr2)*a12 + xr2*yr2*a13;
    float gux15      = yr * zr * (a12+xr2*a13);
    float gux16      = a11 + (xr2+zr2)*a12 + xr2*zr2*a13;
    float gux17      = xr * yr * (3.0f*a12+yr2*a13);
    float gux18      = xr * zr * (a12+yr2*a13);
    float gux19      = xr * yr * (a12+zr2*a13);
    float gux20      = xr * zr * (3.0f*a12+zr2*a13);
    float guy11      = gux12;
    float guy12      = gux14;
    float guy13      = gux15;
    float guy14      = gux17;
    float guy15      = gux18;
    float guy16      = gux19;
    float guy17      = 3.0f*a11 + yr2*(6.0f*a12+yr2*a13);
    float guy18      = yr * zr * (3.0f*a12+yr2*a13);
    float guy19      = a11 + (yr2+zr2)*a12 + yr2*zr2*a13;
    float guy20      = yr * zr * (3.0f*a12+zr2*a13);
    float guz11      = gux13;
    float guz12      = gux15;
    float guz13      = gux16;
    float guz14      = gux18;
    float guz15      = gux19;
    float guz16      = gux20;
    float guz17      = guy18;
    float guz18      = guy19;
    float guz19      = guy20;
    float guz20      = 3.0f*a11 + zr2*(6.0f*a12+zr2*a13);

    float gqxx11     = xr * (12.0f*a21+xr2*(9.0f*a22 + xr2*a23));
    float gqxx12     = yr * (2.0f*a21+xr2*(5.0f*a22  + xr2*a23));
    float gqxx13     = zr * (2.0f*a21+xr2*(5.0f*a22  + xr2*a23));
    float gqxx14     = xr * (2.0f*a21+yr2*2.0f*a22   +xr2*(a22+yr2*a23));
    float gqxx15     = xr * yr * zr * (2.0f*a22+xr2*a23);
    float gqxx16     = xr * (2.0f*a21+zr2*2.0f*a22 +xr2*(a22+zr2*a23));
    float gqxx17     = yr * xr2 * (3.0f*a22+yr2*a23);
    float gqxx18     = zr * xr2 * (a22+yr2*a23);
    float gqxx19     = yr * xr2 * (a22+zr2*a23);
    float gqxx20     = zr * xr2 * (3.0f*a22+zr2*a23);
    float gqxy11     = yr * (3.0f*a21+xr2*(6.0f*a22 +xr2*a23));
    float gqxy12     = xr * (3.0f*(a21+yr2*a22) +xr2*(a22+yr2*a23));
    float gqxy13     = xr * yr * zr * (3.0f*a22+xr2*a23);
    float gqxy14     = yr * (3.0f*(a21+xr2*a22) +yr2*(a22+xr2*a23));
    float gqxy15     = zr * (a21+(yr2+xr2)*a22 +yr2*xr2*a23);
    float gqxy16     = yr * (a21+(xr2+zr2)*a22 +xr2*zr2*a23);
    float gqxy17     = xr * (3.0f*(a21+yr2*a22) +yr2*(3.0f*a22+yr2*a23));
    float gqxy18     = xr * yr * zr * (3.0f*a22+yr2*a23);
    float gqxy19     = xr * (a21+(yr2+zr2)*a22 +yr2*zr2*a23);
    float gqxy20     = xr * yr * zr * (3.0f*a22+zr2*a23);
    float gqxz11     = zr * (3.0f*a21+xr2*(6.0f*a22 +xr2*a23));
    float gqxz12     = xr * yr * zr * (3.0f*a22+xr2*a23);
    float gqxz13     = xr * (3.0f*(a21+zr2*a22) +xr2*(a22+zr2*a23));
    float gqxz14     = zr * (a21+(xr2+yr2)*a22 +xr2*yr2*a23);
    float gqxz15     = yr * (a21+(xr2+zr2)*a22 +zr2*xr2*a23);
    float gqxz16     = zr * (3.0f*(a21+xr2*a22) +zr2*(a22+xr2*a23));
    float gqxz17     = xr * yr * zr * (3.0f*a22+yr2*a23);
    float gqxz18     = xr * (a21+(zr2+yr2)*a22 +zr2*yr2*a23);
    float gqxz19     = xr * yr * zr * (3.0f*a22+zr2*a23);
    float gqxz20     = xr * (3.0f*a21+zr2*(6.0f*a22 +zr2*a23));
    float gqyy11     = xr * yr2 * (3.0f*a22+xr2*a23);
    float gqyy12     = yr * (2.0f*a21+xr2*2.0f*a22 +yr2*(a22+xr2*a23));
    float gqyy13     = zr * yr2 * (a22+xr2*a23);
    float gqyy14     = xr * (2.0f*a21+yr2*(5.0f*a22 +yr2*a23));
    float gqyy15     = xr * yr * zr * (2.0f*a22+yr2*a23);
    float gqyy16     = xr * yr2 * (a22+zr2*a23);
    float gqyy17     = yr * (12.0f*a21+yr2*(9.0f*a22 +yr2*a23));
    float gqyy18     = zr * (2.0f*a21+yr2*(5.0f*a22 +yr2*a23));
    float gqyy19     = yr * (2.0f*a21+zr2*2.0f*a22 +yr2*(a22+zr2*a23));
    float gqyy20     = zr * yr2 * (3.0f*a22+zr2*a23);
    float gqyz11     = xr * yr * zr * (3.0f*a22+xr2*a23);
    float gqyz12     = zr * (a21+(xr2+yr2)*a22 +xr2*yr2*a23);
    float gqyz13     = yr * (a21+(xr2+zr2)*a22 +xr2*zr2*a23);
    float gqyz14     = xr * yr * zr * (3.0f*a22+yr2*a23);
    float gqyz15     = xr * (a21+(yr2+zr2)*a22 +yr2*zr2*a23);
    float gqyz16     = xr * yr * zr * (3.0f*a22+zr2*a23);
    float gqyz17     = zr * (3.0f*a21+yr2*(6.0f*a22 +yr2*a23));
    float gqyz18     = yr * (3.0f*(a21+zr2*a22) +yr2*(a22+zr2*a23));
    float gqyz19     = zr * (3.0f*(a21+yr2*a22) +zr2*(a22+yr2*a23));
    float gqyz20     = yr * (3.0f*a21+zr2*(6.0f*a22 +zr2*a23));
    float gqzz11     = xr * zr2 * (3.0f*a22+xr2*a23);
    float gqzz12     = yr * (zr2*a22+xr2*(zr2*a23));
    float gqzz13     = zr * (2.0f*a21+xr2*2.0f*a22 +zr2*(a22+xr2*a23));
    float gqzz14     = xr * zr2 * (a22+yr2*a23);
    float gqzz15     = xr * yr * zr * (2.0f*a22+zr2*a23);
    float gqzz16     = xr * (2.0f*a21+zr2*(5.0f*a22 +zr2*a23));
    float gqzz17     = yr * zr2 * (3.0f*a22+yr2*a23);
    float gqzz18     = zr * (2.0f*a21+yr2*2.0f*a22 +zr2*(a22+yr2*a23));
    float gqzz19     = yr * (2.0f*a21+zr2*(5.0f*a22 +zr2*a23));
    float gqzz20     = zr * (12.0f*a21+zr2*(9.0f*a22 +zr2*a23));

    // electrostatic solvation energy of the permanent multipoles
    // in their own GK reaction potential

    esym = ci * ck * gc1 - (uxi*(uxk*gux2+uyk*guy2+uzk*guz2)
                           +  uyi*(uxk*gux3+uyk*guy3+uzk*guz3)
                           +  uzi*(uxk*gux4+uyk*guy4+uzk*guz4));

    ewi =  ci*(uxk*gc2+uyk*gc3+uzk*gc4)
          -ck*(uxi*gux1+uyi*guy1+uzi*guz1)
           +ci*(qxxk*gc5+qyyk*gc8+qzzk*gc10
              +2.0f*(qxyk*gc6+qxzk*gc7+qyzk*gc9))
                 +ck*(qxxi*gqxx1+qyyi*gqyy1+qzzi*gqzz1
              +2.0f*(qxyi*gqxy1+qxzi*gqxz1+qyzi*gqyz1))
               - uxi*(qxxk*gux5+qyyk*gux8+qzzk*gux10
              +2.0f*(qxyk*gux6+qxzk*gux7+qyzk*gux9))
               - uyi*(qxxk*guy5+qyyk*guy8+qzzk*guy10
              +2.0f*(qxyk*guy6+qxzk*guy7+qyzk*guy9))
               - uzi*(qxxk*guz5+qyyk*guz8+qzzk*guz10
              +2.0f*(qxyk*guz6+qxzk*guz7+qyzk*guz9))
               + uxk*(qxxi*gqxx2+qyyi*gqyy2+qzzi*gqzz2
              +2.0f*(qxyi*gqxy2+qxzi*gqxz2+qyzi*gqyz2))
               + uyk*(qxxi*gqxx3+qyyi*gqyy3+qzzi*gqzz3
              +2.0f*(qxyi*gqxy3+qxzi*gqxz3+qyzi*gqyz3))
               + uzk*(qxxi*gqxx4+qyyi*gqyy4+qzzi*gqzz4
              +2.0f*(qxyi*gqxy4+qxzi*gqxz4+qyzi*gqyz4))
              + qxxi*(qxxk*gqxx5+qyyk*gqxx8+qzzk*gqxx10
              +2.0f*(qxyk*gqxx6+qxzk*gqxx7+qyzk*gqxx9))
              + qyyi*(qxxk*gqyy5+qyyk*gqyy8+qzzk*gqyy10
              +2.0f*(qxyk*gqyy6+qxzk*gqyy7+qyzk*gqyy9))
              + qzzi*(qxxk*gqzz5+qyyk*gqzz8+qzzk*gqzz10
              +2.0f*(qxyk*gqzz6+qxzk*gqzz7+qyzk*gqzz9))
              + 2.0f*(qxyi*(qxxk*gqxy5+qyyk*gqxy8+qzzk*gqxy10
              +2.0f*(qxyk*gqxy6+qxzk*gqxy7+qyzk*gqxy9))
              + qxzi*(qxxk*gqxz5+qyyk*gqxz8+qzzk*gqxz10
              +2.0f*(qxyk*gqxz6+qxzk*gqxz7+qyzk*gqxz9))
              + qyzi*(qxxk*gqyz5+qyyk*gqyz8+qzzk*gqyz10
              +2.0f*(qxyk*gqyz6+qxzk*gqyz7+qyzk*gqyz9)));

    ewk = ci*(uxk*gux1+uyk*guy1+uzk*guz1)
                      -ck*(uxi*gc2+uyi*gc3+uzi*gc4)
                 +ci*(qxxk*gqxx1+qyyk*gqyy1+qzzk*gqzz1
              +2.0f*(qxyk*gqxy1+qxzk*gqxz1+qyzk*gqyz1))
                 +ck*(qxxi*gc5+qyyi*gc8+qzzi*gc10
              +2.0f*(qxyi*gc6+qxzi*gc7+qyzi*gc9))
               - uxi*(qxxk*gqxx2+qyyk*gqyy2+qzzk*gqzz2
              +2.0f*(qxyk*gqxy2+qxzk*gqxz2+qyzk*gqyz2))
               - uyi*(qxxk*gqxx3+qyyk*gqyy3+qzzk*gqzz3
              +2.0f*(qxyk*gqxy3+qxzk*gqxz3+qyzk*gqyz3))
               - uzi*(qxxk*gqxx4+qyyk*gqyy4+qzzk*gqzz4
              +2.0f*(qxyk*gqxy4+qxzk*gqxz4+qyzk*gqyz4))
               + uxk*(qxxi*gux5+qyyi*gux8+qzzi*gux10
              +2.0f*(qxyi*gux6+qxzi*gux7+qyzi*gux9))
               + uyk*(qxxi*guy5+qyyi*guy8+qzzi*guy10
              +2.0f*(qxyi*guy6+qxzi*guy7+qyzi*guy9))
               + uzk*(qxxi*guz5+qyyi*guz8+qzzi*guz10
              +2.0f*(qxyi*guz6+qxzi*guz7+qyzi*guz9))
              + qxxi*(qxxk*gqxx5+qyyk*gqyy5+qzzk*gqzz5
              +2.0f*(qxyk*gqxy5+qxzk*gqxz5+qyzk*gqyz5))
              + qyyi*(qxxk*gqxx8+qyyk*gqyy8+qzzk*gqzz8
              +2.0f*(qxyk*gqxy8+qxzk*gqxz8+qyzk*gqyz8))
              + qzzi*(qxxk*gqxx10+qyyk*gqyy10+qzzk*gqzz10
              +2.0f*(qxyk*gqxy10+qxzk*gqxz10+qyzk*gqyz10))
       + 2.0f*(qxyi*(qxxk*gqxx6+qyyk*gqyy6+qzzk*gqzz6
              +2.0f*(qxyk*gqxy6+qxzk*gqxz6+qyzk*gqyz6))
              + qxzi*(qxxk*gqxx7+qyyk*gqyy7+qzzk*gqzz7
              +2.0f*(qxyk*gqxy7+qxzk*gqxz7+qyzk*gqyz7))
              + qyzi*(qxxk*gqxx9+qyyk*gqyy9+qzzk*gqzz9
              +2.0f*(qxyk*gqxy9+qxzk*gqxz9+qyzk*gqyz9)));

    desymdx = ci * ck * gc2 - (uxi*(uxk*gux5+uyk*guy5+uzk*guz5)
                              +  uyi*(uxk*gux6+uyk*guy6+uzk*guz6)
                              +  uzi*(uxk*gux7+uyk*guy7+uzk*guz7));

    dewidx = ci*(uxk*gc5+uyk*gc6+uzk*gc7)
                      -ck*(uxi*gux2+uyi*guy2+uzi*guz2)
                 +ci*(qxxk*gc11+qyyk*gc14+qzzk*gc16
              +2.0f*(qxyk*gc12+qxzk*gc13+qyzk*gc15))
                 +ck*(qxxi*gqxx2+qyyi*gqyy2+qzzi*gqzz2
              +2.0f*(qxyi*gqxy2+qxzi*gqxz2+qyzi*gqyz2))
               - uxi*(qxxk*gux11+qyyk*gux14+qzzk*gux16
              +2.0f*(qxyk*gux12+qxzk*gux13+qyzk*gux15))
               - uyi*(qxxk*guy11+qyyk*guy14+qzzk*guy16
              +2.0f*(qxyk*guy12+qxzk*guy13+qyzk*guy15))
               - uzi*(qxxk*guz11+qyyk*guz14+qzzk*guz16
              +2.0f*(qxyk*guz12+qxzk*guz13+qyzk*guz15))
               + uxk*(qxxi*gqxx5+qyyi*gqyy5+qzzi*gqzz5
              +2.0f*(qxyi*gqxy5+qxzi*gqxz5+qyzi*gqyz5))
               + uyk*(qxxi*gqxx6+qyyi*gqyy6+qzzi*gqzz6
              +2.0f*(qxyi*gqxy6+qxzi*gqxz6+qyzi*gqyz6))
               + uzk*(qxxi*gqxx7+qyyi*gqyy7+qzzi*gqzz7
              +2.0f*(qxyi*gqxy7+qxzi*gqxz7+qyzi*gqyz7))
              + qxxi*(qxxk*gqxx11+qyyk*gqxx14+qzzk*gqxx16
              +2.0f*(qxyk*gqxx12+qxzk*gqxx13+qyzk*gqxx15))
              + qyyi*(qxxk*gqyy11+qyyk*gqyy14+qzzk*gqyy16
              +2.0f*(qxyk*gqyy12+qxzk*gqyy13+qyzk*gqyy15))
              + qzzi*(qxxk*gqzz11+qyyk*gqzz14+qzzk*gqzz16
              +2.0f*(qxyk*gqzz12+qxzk*gqzz13+qyzk*gqzz15))
       + 2.0f*(qxyi*(qxxk*gqxy11+qyyk*gqxy14+qzzk*gqxy16
              +2.0f*(qxyk*gqxy12+qxzk*gqxy13+qyzk*gqxy15))
              + qxzi*(qxxk*gqxz11+qyyk*gqxz14+qzzk*gqxz16
              +2.0f*(qxyk*gqxz12+qxzk*gqxz13+qyzk*gqxz15))
              + qyzi*(qxxk*gqyz11+qyyk*gqyz14+qzzk*gqyz16
              +2.0f*(qxyk*gqyz12+qxzk*gqyz13+qyzk*gqyz15)));

    dewkdx = ci*(uxk*gux2+uyk*guy2+uzk*guz2)
                      -ck*(uxi*gc5+uyi*gc6+uzi*gc7)
                 +ci*(qxxk*gqxx2+qyyk*gqyy2+qzzk*gqzz2
              +2.0f*(qxyk*gqxy2+qxzk*gqxz2+qyzk*gqyz2))
                 +ck*(qxxi*gc11+qyyi*gc14+qzzi*gc16
              +2.0f*(qxyi*gc12+qxzi*gc13+qyzi*gc15))
               - uxi*(qxxk*gqxx5+qyyk*gqyy5+qzzk*gqzz5
              +2.0f*(qxyk*gqxy5+qxzk*gqxz5+qyzk*gqyz5))
               - uyi*(qxxk*gqxx6+qyyk*gqyy6+qzzk*gqzz6
              +2.0f*(qxyk*gqxy6+qxzk*gqxz6+qyzk*gqyz6))
               - uzi*(qxxk*gqxx7+qyyk*gqyy7+qzzk*gqzz7
              +2.0f*(qxyk*gqxy7+qxzk*gqxz7+qyzk*gqyz7))
               + uxk*(qxxi*gux11+qyyi*gux14+qzzi*gux16
              +2.0f*(qxyi*gux12+qxzi*gux13+qyzi*gux15))
               + uyk*(qxxi*guy11+qyyi*guy14+qzzi*guy16
              +2.0f*(qxyi*guy12+qxzi*guy13+qyzi*guy15))
               + uzk*(qxxi*guz11+qyyi*guz14+qzzi*guz16
              +2.0f*(qxyi*guz12+qxzi*guz13+qyzi*guz15))
              + qxxi*(qxxk*gqxx11+qyyk*gqyy11+qzzk*gqzz11
              +2.0f*(qxyk*gqxy11+qxzk*gqxz11+qyzk*gqyz11))
              + qyyi*(qxxk*gqxx14+qyyk*gqyy14+qzzk*gqzz14
              +2.0f*(qxyk*gqxy14+qxzk*gqxz14+qyzk*gqyz14))
              + qzzi*(qxxk*gqxx16+qyyk*gqyy16+qzzk*gqzz16
              +2.0f*(qxyk*gqxy16+qxzk*gqxz16+qyzk*gqyz16))
       + 2.0f*(qxyi*(qxxk*gqxx12+qyyk*gqyy12+qzzk*gqzz12
              +2.0f*(qxyk*gqxy12+qxzk*gqxz12+qyzk*gqyz12))
              + qxzi*(qxxk*gqxx13+qyyk*gqyy13+qzzk*gqzz13
              +2.0f*(qxyk*gqxy13+qxzk*gqxz13+qyzk*gqyz13))
              + qyzi*(qxxk*gqxx15+qyyk*gqyy15+qzzk*gqzz15
              +2.0f*(qxyk*gqxy15+qxzk*gqxz15+qyzk*gqyz15)));

    dedx = desymdx + 0.5f*(dewidx + dewkdx);

    desymdy = ci * ck * gc3
                           - (uxi*(uxk*gux6+uyk*guy6+uzk*guz6)
                             +uyi*(uxk*gux8+uyk*guy8+uzk*guz8)
                             +uzi*(uxk*gux9+uyk*guy9+uzk*guz9));

    dewidy = ci*(uxk*gc6+uyk*gc8+uzk*gc9)
                      -ck*(uxi*gux3+uyi*guy3+uzi*guz3)
                 +ci*(qxxk*gc12+qyyk*gc17+qzzk*gc19
              +2.0f*(qxyk*gc14+qxzk*gc15+qyzk*gc18))
                 +ck*(qxxi*gqxx3+qyyi*gqyy3+qzzi*gqzz3
              +2.0f*(qxyi*gqxy3+qxzi*gqxz3+qyzi*gqyz3))
               - uxi*(qxxk*gux12+qyyk*gux17+qzzk*gux19
              +2.0f*(qxyk*gux14+qxzk*gux15+qyzk*gux18))
               - uyi*(qxxk*guy12+qyyk*guy17+qzzk*guy19
              +2.0f*(qxyk*guy14+qxzk*guy15+qyzk*guy18))
               - uzi*(qxxk*guz12+qyyk*guz17+qzzk*guz19
              +2.0f*(qxyk*guz14+qxzk*guz15+qyzk*guz18))
               + uxk*(qxxi*gqxx6+qyyi*gqyy6+qzzi*gqzz6
              +2.0f*(qxyi*gqxy6+qxzi*gqxz6+qyzi*gqyz6))
               + uyk*(qxxi*gqxx8+qyyi*gqyy8+qzzi*gqzz8
              +2.0f*(qxyi*gqxy8+qxzi*gqxz8+qyzi*gqyz8))
               + uzk*(qxxi*gqxx9+qyyi*gqyy9+qzzi*gqzz9
              +2.0f*(qxyi*gqxy9+qxzi*gqxz9+qyzi*gqyz9))
              + qxxi*(qxxk*gqxx12+qyyk*gqxx17+qzzk*gqxx19
              +2.0f*(qxyk*gqxx14+qxzk*gqxx15+qyzk*gqxx18))
              + qyyi*(qxxk*gqyy12+qyyk*gqyy17+qzzk*gqyy19
              +2.0f*(qxyk*gqyy14+qxzk*gqyy15+qyzk*gqyy18))
              + qzzi*(qxxk*gqzz12+qyyk*gqzz17+qzzk*gqzz19
              +2.0f*(qxyk*gqzz14+qxzk*gqzz15+qyzk*gqzz18))
       + 2.0f*(qxyi*(qxxk*gqxy12+qyyk*gqxy17+qzzk*gqxy19
              +2.0f*(qxyk*gqxy14+qxzk*gqxy15+qyzk*gqxy18))
              + qxzi*(qxxk*gqxz12+qyyk*gqxz17+qzzk*gqxz19
              +2.0f*(qxyk*gqxz14+qxzk*gqxz15+qyzk*gqxz18))
              + qyzi*(qxxk*gqyz12+qyyk*gqyz17+qzzk*gqyz19
              +2.0f*(qxyk*gqyz14+qxzk*gqyz15+qyzk*gqyz18)));

    dewkdy = ci*(uxk*gux3+uyk*guy3+uzk*guz3)
                      -ck*(uxi*gc6+uyi*gc8+uzi*gc9)
                 +ci*(qxxk*gqxx3+qyyk*gqyy3+qzzk*gqzz3
              +2.0f*(qxyk*gqxy3+qxzk*gqxz3+qyzk*gqyz3))
                 +ck*(qxxi*gc12+qyyi*gc17+qzzi*gc19
              +2.0f*(qxyi*gc14+qxzi*gc15+qyzi*gc18))
               - uxi*(qxxk*gqxx6+qyyk*gqyy6+qzzk*gqzz6
              +2.0f*(qxyk*gqxy6+qxzk*gqxz6+qyzk*gqyz6))
               - uyi*(qxxk*gqxx8+qyyk*gqyy8+qzzk*gqzz8
              +2.0f*(qxyk*gqxy8+qxzk*gqxz8+qyzk*gqyz8))
               - uzi*(qxxk*gqxx9+qyyk*gqyy9+qzzk*gqzz9
              +2.0f*(qxyk*gqxy9+qxzk*gqxz9+qyzk*gqyz9))
               + uxk*(qxxi*gux12+qyyi*gux17+qzzi*gux19
              +2.0f*(qxyi*gux14+qxzi*gux15+qyzi*gux18))
               + uyk*(qxxi*guy12+qyyi*guy17+qzzi*guy19
              +2.0f*(qxyi*guy14+qxzi*guy15+qyzi*guy18))
               + uzk*(qxxi*guz12+qyyi*guz17+qzzi*guz19
              +2.0f*(qxyi*guz14+qxzi*guz15+qyzi*guz18))
              + qxxi*(qxxk*gqxx12+qyyk*gqyy12+qzzk*gqzz12
              +2.0f*(qxyk*gqxy12+qxzk*gqxz12+qyzk*gqyz12))
              + qyyi*(qxxk*gqxx17+qyyk*gqyy17+qzzk*gqzz17
              +2.0f*(qxyk*gqxy17+qxzk*gqxz17+qyzk*gqyz17))
              + qzzi*(qxxk*gqxx19+qyyk*gqyy19+qzzk*gqzz19
              +2.0f*(qxyk*gqxy19+qxzk*gqxz19+qyzk*gqyz19))
       + 2.0f*(qxyi*(qxxk*gqxx14+qyyk*gqyy14+qzzk*gqzz14
              +2.0f*(qxyk*gqxy14+qxzk*gqxz14+qyzk*gqyz14))
              + qxzi*(qxxk*gqxx15+qyyk*gqyy15+qzzk*gqzz15
              +2.0f*(qxyk*gqxy15+qxzk*gqxz15+qyzk*gqyz15))
              + qyzi*(qxxk*gqxx18+qyyk*gqyy18+qzzk*gqzz18
              +2.0f*(qxyk*gqxy18+qxzk*gqxz18+qyzk*gqyz18)));

    dedy = desymdy + 0.5f*(dewidy + dewkdy);

    desymdz = ci * ck * gc4
                           - (uxi*(uxk*gux7+uyk*guy7+uzk*guz7)
                             +uyi*(uxk*gux9+uyk*guy9+uzk*guz9)
                             +uzi*(uxk*gux10+uyk*guy10+uzk*guz10));

    dewidz = ci*(uxk*gc7+uyk*gc9+uzk*gc10)
                      -ck*(uxi*gux4+uyi*guy4+uzi*guz4)
                 +ci*(qxxk*gc13+qyyk*gc18+qzzk*gc20
              +2.0f*(qxyk*gc15+qxzk*gc16+qyzk*gc19))
                 +ck*(qxxi*gqxx4+qyyi*gqyy4+qzzi*gqzz4
              +2.0f*(qxyi*gqxy4+qxzi*gqxz4+qyzi*gqyz4))
               - uxi*(qxxk*gux13+qyyk*gux18+qzzk*gux20
              +2.0f*(qxyk*gux15+qxzk*gux16+qyzk*gux19))
               - uyi*(qxxk*guy13+qyyk*guy18+qzzk*guy20
              +2.0f*(qxyk*guy15+qxzk*guy16+qyzk*guy19))
               - uzi*(qxxk*guz13+qyyk*guz18+qzzk*guz20
              +2.0f*(qxyk*guz15+qxzk*guz16+qyzk*guz19))
               + uxk*(qxxi*gqxx7+qyyi*gqyy7+qzzi*gqzz7
              +2.0f*(qxyi*gqxy7+qxzi*gqxz7+qyzi*gqyz7))
               + uyk*(qxxi*gqxx9+qyyi*gqyy9+qzzi*gqzz9
              +2.0f*(qxyi*gqxy9+qxzi*gqxz9+qyzi*gqyz9))
               + uzk*(qxxi*gqxx10+qyyi*gqyy10+qzzi*gqzz10
              +2.0f*(qxyi*gqxy10+qxzi*gqxz10+qyzi*gqyz10))
              + qxxi*(qxxk*gqxx13+qyyk*gqxx18+qzzk*gqxx20
              +2.0f*(qxyk*gqxx15+qxzk*gqxx16+qyzk*gqxx19))
              + qyyi*(qxxk*gqyy13+qyyk*gqyy18+qzzk*gqyy20
              +2.0f*(qxyk*gqyy15+qxzk*gqyy16+qyzk*gqyy19))
              + qzzi*(qxxk*gqzz13+qyyk*gqzz18+qzzk*gqzz20
              +2.0f*(qxyk*gqzz15+qxzk*gqzz16+qyzk*gqzz19))
       + 2.0f*(qxyi*(qxxk*gqxy13+qyyk*gqxy18+qzzk*gqxy20
              +2.0f*(qxyk*gqxy15+qxzk*gqxy16+qyzk*gqxy19))
              + qxzi*(qxxk*gqxz13+qyyk*gqxz18+qzzk*gqxz20
              +2.0f*(qxyk*gqxz15+qxzk*gqxz16+qyzk*gqxz19))
              + qyzi*(qxxk*gqyz13+qyyk*gqyz18+qzzk*gqyz20
              +2.0f*(qxyk*gqyz15+qxzk*gqyz16+qyzk*gqyz19)));

    dewkdz = ci*(uxk*gux4+uyk*guy4+uzk*guz4)
                      -ck*(uxi*gc7+uyi*gc9+uzi*gc10)
                 +ci*(qxxk*gqxx4+qyyk*gqyy4+qzzk*gqzz4
              +2.0f*(qxyk*gqxy4+qxzk*gqxz4+qyzk*gqyz4))
                 +ck*(qxxi*gc13+qyyi*gc18+qzzi*gc20
              +2.0f*(qxyi*gc15+qxzi*gc16+qyzi*gc19))
               - uxi*(qxxk*gqxx7+qyyk*gqyy7+qzzk*gqzz7
              +2.0f*(qxyk*gqxy7+qxzk*gqxz7+qyzk*gqyz7))
               - uyi*(qxxk*gqxx9+qyyk*gqyy9+qzzk*gqzz9
              +2.0f*(qxyk*gqxy9+qxzk*gqxz9+qyzk*gqyz9))
               - uzi*(qxxk*gqxx10+qyyk*gqyy10+qzzk*gqzz10
              +2.0f*(qxyk*gqxy10+qxzk*gqxz10+qyzk*gqyz10))
               + uxk*(qxxi*gux13+qyyi*gux18+qzzi*gux20
              +2.0f*(qxyi*gux15+qxzi*gux16+qyzi*gux19))
               + uyk*(qxxi*guy13+qyyi*guy18+qzzi*guy20
              +2.0f*(qxyi*guy15+qxzi*guy16+qyzi*guy19))
               + uzk*(qxxi*guz13+qyyi*guz18+qzzi*guz20
              +2.0f*(qxyi*guz15+qxzi*guz16+qyzi*guz19))
              + qxxi*(qxxk*gqxx13+qyyk*gqyy13+qzzk*gqzz13
              +2.0f*(qxyk*gqxy13+qxzk*gqxz13+qyzk*gqyz13))
              + qyyi*(qxxk*gqxx18+qyyk*gqyy18+qzzk*gqzz18
              +2.0f*(qxyk*gqxy18+qxzk*gqxz18+qyzk*gqyz18))
              + qzzi*(qxxk*gqxx20+qyyk*gqyy20+qzzk*gqzz20
              +2.0f*(qxyk*gqxy20+qxzk*gqxz20+qyzk*gqyz20))
       + 2.0f*(qxyi*(qxxk*gqxx15+qyyk*gqyy15+qzzk*gqzz15
              +2.0f*(qxyk*gqxy15+qxzk*gqxz15+qyzk*gqyz15))
              + qxzi*(qxxk*gqxx16+qyyk*gqyy16+qzzk*gqzz16
              +2.0f*(qxyk*gqxy16+qxzk*gqxz16+qyzk*gqyz16))
              + qyzi*(qxxk*gqxx19+qyyk*gqyy19+qzzk*gqzz19
              +2.0f*(qxyk*gqxy19+qxzk*gqxz19+qyzk*gqyz19)));

    dedz = desymdz + 0.5f*(dewidz + dewkdz);

    desymdr = ci * ck * gc21
                           - (uxi*(uxk*gux22+uyk*guy22+uzk*guz22)
                             +uyi*(uxk*gux23+uyk*guy23+uzk*guz23)
                             +uzi*(uxk*gux24+uyk*guy24+uzk*guz24));

    dewidr = ci*(uxk*gc22+uyk*gc23+uzk*gc24)
                      -ck*(uxi*gux21+uyi*guy21+uzi*guz21)
                 +ci*(qxxk*gc25+qyyk*gc28+qzzk*gc30
              +2.0f*(qxyk*gc26+qxzk*gc27+qyzk*gc29))
                 +ck*(qxxi*gqxx21+qyyi*gqyy21+qzzi*gqzz21
              +2.0f*(qxyi*gqxy21+qxzi*gqxz21+qyzi*gqyz21))
               - uxi*(qxxk*gux25+qyyk*gux28+qzzk*gux30
              +2.0f*(qxyk*gux26+qxzk*gux27+qyzk*gux29))
               - uyi*(qxxk*guy25+qyyk*guy28+qzzk*guy30
              +2.0f*(qxyk*guy26+qxzk*guy27+qyzk*guy29))
               - uzi*(qxxk*guz25+qyyk*guz28+qzzk*guz30
              +2.0f*(qxyk*guz26+qxzk*guz27+qyzk*guz29))
               + uxk*(qxxi*gqxx22+qyyi*gqyy22+qzzi*gqzz22
              +2.0f*(qxyi*gqxy22+qxzi*gqxz22+qyzi*gqyz22))
               + uyk*(qxxi*gqxx23+qyyi*gqyy23+qzzi*gqzz23
              +2.0f*(qxyi*gqxy23+qxzi*gqxz23+qyzi*gqyz23))
               + uzk*(qxxi*gqxx24+qyyi*gqyy24+qzzi*gqzz24
              +2.0f*(qxyi*gqxy24+qxzi*gqxz24+qyzi*gqyz24))
              + qxxi*(qxxk*gqxx25+qyyk*gqxx28+qzzk*gqxx30
              +2.0f*(qxyk*gqxx26+qxzk*gqxx27+qyzk*gqxx29))
              + qyyi*(qxxk*gqyy25+qyyk*gqyy28+qzzk*gqyy30
              +2.0f*(qxyk*gqyy26+qxzk*gqyy27+qyzk*gqyy29))
              + qzzi*(qxxk*gqzz25+qyyk*gqzz28+qzzk*gqzz30
              +2.0f*(qxyk*gqzz26+qxzk*gqzz27+qyzk*gqzz29))
              + 2.0f*(qxyi*(qxxk*gqxy25+qyyk*gqxy28+qzzk*gqxy30
              +2.0f*(qxyk*gqxy26+qxzk*gqxy27+qyzk*gqxy29))
              + qxzi*(qxxk*gqxz25+qyyk*gqxz28+qzzk*gqxz30
              +2.0f*(qxyk*gqxz26+qxzk*gqxz27+qyzk*gqxz29))
              + qyzi*(qxxk*gqyz25+qyyk*gqyz28+qzzk*gqyz30
              +2.0f*(qxyk*gqyz26+qxzk*gqyz27+qyzk*gqyz29)));

    dewkdr = ci*(uxk*gux21+uyk*guy21+uzk*guz21)
                      -ck*(uxi*gc22+uyi*gc23+uzi*gc24)
                 +ci*(qxxk*gqxx21+qyyk*gqyy21+qzzk*gqzz21
              +2.0f*(qxyk*gqxy21+qxzk*gqxz21+qyzk*gqyz21))
                 +ck*(qxxi*gc25+qyyi*gc28+qzzi*gc30
              +2.0f*(qxyi*gc26+qxzi*gc27+qyzi*gc29))
               - uxi*(qxxk*gqxx22+qyyk*gqyy22+qzzk*gqzz22
              +2.0f*(qxyk*gqxy22+qxzk*gqxz22+qyzk*gqyz22))
               - uyi*(qxxk*gqxx23+qyyk*gqyy23+qzzk*gqzz23
              +2.0f*(qxyk*gqxy23+qxzk*gqxz23+qyzk*gqyz23))
               - uzi*(qxxk*gqxx24+qyyk*gqyy24+qzzk*gqzz24
              +2.0f*(qxyk*gqxy24+qxzk*gqxz24+qyzk*gqyz24))
               + uxk*(qxxi*gux25+qyyi*gux28+qzzi*gux30
              +2.0f*(qxyi*gux26+qxzi*gux27+qyzi*gux29))
               + uyk*(qxxi*guy25+qyyi*guy28+qzzi*guy30
              +2.0f*(qxyi*guy26+qxzi*guy27+qyzi*guy29))
               + uzk*(qxxi*guz25+qyyi*guz28+qzzi*guz30
              +2.0f*(qxyi*guz26+qxzi*guz27+qyzi*guz29))
              + qxxi*(qxxk*gqxx25+qyyk*gqyy25+qzzk*gqzz25
              +2.0f*(qxyk*gqxy25+qxzk*gqxz25+qyzk*gqyz25))
              + qyyi*(qxxk*gqxx28+qyyk*gqyy28+qzzk*gqzz28
              +2.0f*(qxyk*gqxy28+qxzk*gqxz28+qyzk*gqyz28))
              + qzzi*(qxxk*gqxx30+qyyk*gqyy30+qzzk*gqzz30
              +2.0f*(qxyk*gqxy30+qxzk*gqxz30+qyzk*gqyz30))
              + 2.0f*(qxyi*(qxxk*gqxx26+qyyk*gqyy26+qzzk*gqzz26
              +2.0f*(qxyk*gqxy26+qxzk*gqxz26+qyzk*gqyz26))
              + qxzi*(qxxk*gqxx27+qyyk*gqyy27+qzzk*gqzz27
              +2.0f*(qxyk*gqxy27+qxzk*gqxz27+qyzk*gqyz27))
              + qyzi*(qxxk*gqxx29+qyyk*gqyy29+qzzk*gqzz29
              +2.0f*(qxyk*gqxy29+qxzk*gqxz29+qyzk*gqyz29)));

    dsumdr = desymdr + 0.5f*(dewidr + dewkdr);
    drbi = rbk*dsumdr;
    drbk = rbi*dsumdr;

    // torque on permanent dipoles due to permanent reaction field

    float trq1   = 0.0f;
    float trq2   = 0.0f;
    float trq3   = 0.0f;

    float trq_k1 = 0.0f;
    float trq_k2 = 0.0f;
    float trq_k3 = 0.0f;

    if ( sameAtom == 0 )
    {

        float fid1 = uxk*gux2 + uyk*gux3 + uzk*gux4
                + 0.5f*(ck*gux1+qxxk*gux5+qyyk*gux8+qzzk*gux10
                      +2.0f*(qxyk*gux6+qxzk*gux7+qyzk*gux9)
                      +ck*gc2+qxxk*gqxx2+qyyk*gqyy2+qzzk*gqzz2
                      +2.0f*(qxyk*gqxy2+qxzk*gqxz2+qyzk*gqyz2));

        float fid2 = uxk*guy2 + uyk*guy3 + uzk*guy4
                + 0.5f*(ck*guy1+qxxk*guy5+qyyk*guy8+qzzk*guy10
                      +2.0f*(qxyk*guy6+qxzk*guy7+qyzk*guy9)
                      +ck*gc3+qxxk*gqxx3+qyyk*gqyy3+qzzk*gqzz3
                      +2.0f*(qxyk*gqxy3+qxzk*gqxz3+qyzk*gqyz3));

        float fid3 = uxk*guz2 + uyk*guz3 + uzk*guz4
                + 0.5f*(ck*guz1+qxxk*guz5+qyyk*guz8+qzzk*guz10
                      +2.0f*(qxyk*guz6+qxzk*guz7+qyzk*guz9)
                      +ck*gc4+qxxk*gqxx4+qyyk*gqyy4+qzzk*gqzz4
                      +2.0f*(qxyk*gqxy4+qxzk*gqxz4+qyzk*gqyz4));

        float fkd1 = uxi*gux2 + uyi*gux3 + uzi*gux4
                - 0.5f*(ci*gux1+qxxi*gux5+qyyi*gux8+qzzi*gux10
                      +2.0f*(qxyi*gux6+qxzi*gux7+qyzi*gux9)
                      +ci*gc2+qxxi*gqxx2+qyyi*gqyy2+qzzi*gqzz2
                      +2.0f*(qxyi*gqxy2+qxzi*gqxz2+qyzi*gqyz2));

        float fkd2 = uxi*guy2 + uyi*guy3 + uzi*guy4
                - 0.5f*(ci*guy1+qxxi*guy5+qyyi*guy8+qzzi*guy10
                      +2.0f*(qxyi*guy6+qxzi*guy7+qyzi*guy9)
                      +ci*gc3+qxxi*gqxx3+qyyi*gqyy3+qzzi*gqzz3
                      +2.0f*(qxyi*gqxy3+qxzi*gqxz3+qyzi*gqyz3));

        float fkd3 = uxi*guz2 + uyi*guz3 + uzi*guz4
                - 0.5f*(ci*guz1+qxxi*guz5+qyyi*guz8+qzzi*guz10
                      +2.0f*(qxyi*guz6+qxzi*guz7+qyzi*guz9)
                      +ci*gc4+qxxi*gqxx4+qyyi*gqyy4+qzzi*gqzz4
                      +2.0f*(qxyi*gqxy4+qxzi*gqxz4+qyzi*gqyz4));

        trq1    = uyi*fid3 - uzi*fid2;
        trq2    = uzi*fid1 - uxi*fid3;
        trq3    = uxi*fid2 - uyi*fid1;

        trq_k1  = uyk*fkd3 - uzk*fkd2;
        trq_k2  = uzk*fkd1 - uxk*fkd3;
        trq_k3  = uxk*fkd2 - uyk*fkd1;

        // torque on quadrupoles due to permanent reaction field gradient

        float fidg11 =
                - 0.5f*(ck*gqxx1+uxk*gqxx2+uyk*gqxx3+uzk*gqxx4
                      +qxxk*gqxx5+qyyk*gqxx8+qzzk*gqxx10
                      +2.0f*(qxyk*gqxx6+qxzk*gqxx7+qyzk*gqxx9)
                      +ck*gc5+uxk*gux5+uyk*guy5+uzk*guz5
                      +qxxk*gqxx5+qyyk*gqyy5+qzzk*gqzz5
                      +2.0f*(qxyk*gqxy5+qxzk*gqxz5+qyzk*gqyz5));

        float fidg12 =
                - 0.5f*(ck*gqxy1+uxk*gqxy2+uyk*gqxy3+uzk*gqxy4
                      +qxxk*gqxy5+qyyk*gqxy8+qzzk*gqxy10
                      +2.0f*(qxyk*gqxy6+qxzk*gqxy7+qyzk*gqxy9)
                      +ck*gc6+uxk*gux6+uyk*guy6+uzk*guz6
                      +qxxk*gqxx6+qyyk*gqyy6+qzzk*gqzz6
                      +2.0f*(qxyk*gqxy6+qxzk*gqxz6+qyzk*gqyz6));

        float fidg13 =
                - 0.5f*(ck*gqxz1+uxk*gqxz2+uyk*gqxz3+uzk*gqxz4
                      +qxxk*gqxz5+qyyk*gqxz8+qzzk*gqxz10
                      +2.0f*(qxyk*gqxz6+qxzk*gqxz7+qyzk*gqxz9)
                      +ck*gc7+uxk*gux7+uyk*guy7+uzk*guz7
                      +qxxk*gqxx7+qyyk*gqyy7+qzzk*gqzz7
                      +2.0f*(qxyk*gqxy7+qxzk*gqxz7+qyzk*gqyz7));

        float fidg22 =
                - 0.5f*(ck*gqyy1+uxk*gqyy2+uyk*gqyy3+uzk*gqyy4
                      +qxxk*gqyy5+qyyk*gqyy8+qzzk*gqyy10
                      +2.0f*(qxyk*gqyy6+qxzk*gqyy7+qyzk*gqyy9)
                      +ck*gc8+uxk*gux8+uyk*guy8+uzk*guz8
                      +qxxk*gqxx8+qyyk*gqyy8+qzzk*gqzz8
                      +2.0f*(qxyk*gqxy8+qxzk*gqxz8+qyzk*gqyz8));

        float fidg23 =
                - 0.5f*(ck*gqyz1+uxk*gqyz2+uyk*gqyz3+uzk*gqyz4
                      +qxxk*gqyz5+qyyk*gqyz8+qzzk*gqyz10
                      +2.0f*(qxyk*gqyz6+qxzk*gqyz7+qyzk*gqyz9)
                      +ck*gc9+uxk*gux9+uyk*guy9+uzk*guz9
                      +qxxk*gqxx9+qyyk*gqyy9+qzzk*gqzz9
                      +2.0f*(qxyk*gqxy9+qxzk*gqxz9+qyzk*gqyz9));

        float fidg33 =
                - 0.5f*(ck*gqzz1+uxk*gqzz2+uyk*gqzz3+uzk*gqzz4
                      +qxxk*gqzz5+qyyk*gqzz8+qzzk*gqzz10
                      +2.0f*(qxyk*gqzz6+qxzk*gqzz7+qyzk*gqzz9)
                      +ck*gc10+uxk*gux10+uyk*guy10+uzk*guz10
                      +qxxk*gqxx10+qyyk*gqyy10+qzzk*gqzz10
                   +2.0f*(qxyk*gqxy10+qxzk*gqxz10+qyzk*gqyz10));

        float fidg21 = fidg12;
        float fidg31 = fidg13;
        float fidg32 = fidg23;

        float fkdg11 =
                - 0.5f*(ci*gqxx1-uxi*gqxx2-uyi*gqxx3-uzi *gqxx4
                      +qxxi*gqxx5+qyyi*gqxx8+qzzi*gqxx10
                      +2.0f*(qxyi*gqxx6+qxzi*gqxx7+qyzi*gqxx9)
                      +ci*gc5-uxi*gux5-uyi*guy5-uzi*guz5
                      +qxxi*gqxx5+qyyi*gqyy5+qzzi*gqzz5
                      +2.0f*(qxyi*gqxy5+qxzi*gqxz5+qyzi*gqyz5));

        float fkdg12 =
                - 0.5f*(ci*gqxy1-uxi*gqxy2-uyi*gqxy3-uzi*gqxy4
                      +qxxi*gqxy5+qyyi*gqxy8+qzzi*gqxy10
                      +2.0f*(qxyi*gqxy6+qxzi*gqxy7+qyzi*gqxy9)
                      +ci*gc6-uxi*gux6-uyi*guy6-uzi*guz6
                      +qxxi*gqxx6+qyyi*gqyy6+qzzi*gqzz6
                      +2.0f*(qxyi*gqxy6+qxzi*gqxz6+qyzi*gqyz6));

        float fkdg13 =
                - 0.5f*(ci*gqxz1-uxi*gqxz2-uyi*gqxz3-uzi*gqxz4
                      +qxxi*gqxz5+qyyi*gqxz8+qzzi*gqxz10
                      +2.0f*(qxyi*gqxz6+qxzi*gqxz7+qyzi*gqxz9)
                      +ci*gc7-uxi*gux7-uyi*guy7-uzi*guz7
                      +qxxi*gqxx7+qyyi*gqyy7+qzzi*gqzz7
                      +2.0f*(qxyi*gqxy7+qxzi*gqxz7+qyzi*gqyz7));

        float fkdg22 =
                - 0.5f*(ci*gqyy1-uxi*gqyy2-uyi*gqyy3-uzi*gqyy4
                      +qxxi*gqyy5+qyyi*gqyy8+qzzi*gqyy10
                      +2.0f*(qxyi*gqyy6+qxzi*gqyy7+qyzi*gqyy9)
                      +ci*gc8-uxi*gux8-uyi*guy8-uzi*guz8
                      +qxxi*gqxx8+qyyi*gqyy8+qzzi*gqzz8
                      +2.0f*(qxyi*gqxy8+qxzi*gqxz8+qyzi*gqyz8));

        float fkdg23 =
                - 0.5f*(ci*gqyz1-uxi*gqyz2-uyi*gqyz3-uzi*gqyz4
                      +qxxi*gqyz5+qyyi*gqyz8+qzzi*gqyz10
                      +2.0f*(qxyi*gqyz6+qxzi*gqyz7+qyzi*gqyz9)
                      +ci*gc9-uxi*gux9-uyi*guy9-uzi*guz9
                      +qxxi*gqxx9+qyyi*gqyy9+qzzi*gqzz9
                      +2.0f*(qxyi*gqxy9+qxzi*gqxz9+qyzi*gqyz9));
        float fkdg33 =
                - 0.5f*(ci*gqzz1-uxi*gqzz2-uyi*gqzz3-uzi*gqzz4
                      +qxxi*gqzz5+qyyi*gqzz8+qzzi*gqzz10
                      +2.0f*(qxyi*gqzz6+qxzi*gqzz7+qyzi*gqzz9)
                      +ci*gc10-uxi*gux10-uyi*guy10-uzi*guz10
                      +qxxi*gqxx10+qyyi*gqyy10+qzzi*gqzz10
                    +2.0f*(qxyi*gqxy10+qxzi*gqxz10+qyzi*gqyz10));

        float fkdg21 = fkdg12;
        float fkdg31 = fkdg13;
        float fkdg32 = fkdg23;

        trq1   += 2.0f* (qxyi*fidg13+qyyi*fidg23+qyzi*fidg33
                           -qxzi*fidg12-qyzi*fidg22-qzzi*fidg32);

        trq2   += 2.0f*(qxzi*fidg11+qyzi*fidg21+qzzi*fidg31
                         -qxxi*fidg13-qxyi*fidg23-qxzi*fidg33);

        trq3   += 2.0f*(qxxi*fidg12+qxyi*fidg22+qxzi*fidg32
                         -qxyi*fidg11-qyyi*fidg21-qyzi*fidg31);

        trq_k1 += 2.0f*
                          (qxyk*fkdg13+qyyk*fkdg23+qyzk*fkdg33
                          -qxzk*fkdg12-qyzk*fkdg22-qzzk*fkdg32);

        trq_k2 += 2.0f*
                          (qxzk*fkdg11+qyzk*fkdg21+qzzk*fkdg31
                          -qxxk*fkdg13-qxyk*fkdg23-qxzk*fkdg33);

        trq_k3 += 2.0f*
                          (qxxk*fkdg12+qxyk*fkdg22+qxzk*fkdg32
                          -qxyk*fkdg11-qyyk*fkdg21-qyzk*fkdg31);
    }

    // electrostatic solvation energy of the permanent multipoles in
    // the GK reaction potential of the induced dipoles

    esymi =              -uxi*(dxk*gux2+dyk*guy2+dzk*guz2)
                        - uyi*(dxk*gux3+dyk*guy3+dzk*guz3)
                        - uzi*(dxk*gux4+dyk*guy4+dzk*guz4)
                        - uxk*(dxi*gux2+dyi*guy2+dzi*guz2)
                        - uyk*(dxi*gux3+dyi*guy3+dzi*guz3)
                        - uzk*(dxi*gux4+dyi*guy4+dzi*guz4);

    ewii = ci*(dxk*gc2+dyk*gc3+dzk*gc4)
                      - ck*(dxi*gux1+dyi*guy1+dzi*guz1)
                      - dxi*(qxxk*gux5+qyyk*gux8+qzzk*gux10
                     +2.0f*(qxyk*gux6+qxzk*gux7+qyzk*gux9))
                      - dyi*(qxxk*guy5+qyyk*guy8+qzzk*guy10
                     +2.0f*(qxyk*guy6+qxzk*guy7+qyzk*guy9))
                      - dzi*(qxxk*guz5+qyyk*guz8+qzzk*guz10
                     +2.0f*(qxyk*guz6+qxzk*guz7+qyzk*guz9))
                      + dxk*(qxxi*gqxx2+qyyi*gqyy2+qzzi*gqzz2
                     +2.0f*(qxyi*gqxy2+qxzi*gqxz2+qyzi*gqyz2))
                      + dyk*(qxxi*gqxx3+qyyi*gqyy3+qzzi*gqzz3
                     +2.0f*(qxyi*gqxy3+qxzi*gqxz3+qyzi*gqyz3))
                      + dzk*(qxxi*gqxx4+qyyi*gqyy4+qzzi*gqzz4
                     +2.0f*(qxyi*gqxy4+qxzi*gqxz4+qyzi*gqyz4));

    ewki = ci*(dxk*gux1+dyk*guy1+dzk*guz1)
                      - ck*(dxi*gc2+dyi*gc3+dzi*gc4)
                      - dxi*(qxxk*gqxx2+qyyk*gqyy2+qzzk*gqzz2
                     +2.0f*(qxyk*gqxy2+qxzk*gqxz2+qyzk*gqyz2))
                      - dyi*(qxxk*gqxx3+qyyk*gqyy3+qzzk*gqzz3
                     +2.0f*(qxyk*gqxy3+qxzk*gqxz3+qyzk*gqyz3))
                      - dzi*(qxxk*gqxx4+qyyk*gqyy4+qzzk*gqzz4
                     +2.0f*(qxyk*gqxy4+qxzk*gqxz4+qyzk*gqyz4))
                      + dxk*(qxxi*gux5+qyyi*gux8+qzzi*gux10
                     +2.0f*(qxyi*gux6+qxzi*gux7+qyzi*gux9))
                      + dyk*(qxxi*guy5+qyyi*guy8+qzzi*guy10
                     +2.0f*(qxyi*guy6+qxzi*guy7+qyzi*guy9))
                      + dzk*(qxxi*guz5+qyyi*guz8+qzzi*guz10
                     +2.0f*(qxyi*guz6+qxzi*guz7+qyzi*guz9));

    // electrostatic solvation free energy gradient of the permanent
    // multipoles in the reaction potential of the induced dipoles

    dpsymdx = -uxi*(sxk*gux5+syk*guy5+szk*guz5)
                          - uyi*(sxk*gux6+syk*guy6+szk*guz6)
                          - uzi*(sxk*gux7+syk*guy7+szk*guz7)
                          - uxk*(sxi*gux5+syi*guy5+szi*guz5)
                          - uyk*(sxi*gux6+syi*guy6+szi*guz6)
                          - uzk*(sxi*gux7+syi*guy7+szi*guz7);

    dpwidx = ci*(sxk*gc5+syk*gc6+szk*gc7)
                        - ck*(sxi*gux2+syi*guy2+szi*guz2)
                      - sxi*(qxxk*gux11+qyyk*gux14+qzzk*gux16
                     +2.0f*(qxyk*gux12+qxzk*gux13+qyzk*gux15))
                      - syi*(qxxk*guy11+qyyk*guy14+qzzk*guy16
                     +2.0f*(qxyk*guy12+qxzk*guy13+qyzk*guy15))
                      - szi*(qxxk*guz11+qyyk*guz14+qzzk*guz16
                     +2.0f*(qxyk*guz12+qxzk*guz13+qyzk*guz15))
                      + sxk*(qxxi*gqxx5+qyyi*gqyy5+qzzi*gqzz5
                     +2.0f*(qxyi*gqxy5+qxzi*gqxz5+qyzi*gqyz5))
                      + syk*(qxxi*gqxx6+qyyi*gqyy6+qzzi*gqzz6
                     +2.0f*(qxyi*gqxy6+qxzi*gqxz6+qyzi*gqyz6))
                      + szk*(qxxi*gqxx7+qyyi*gqyy7+qzzi*gqzz7
                     +2.0f*(qxyi*gqxy7+qxzi*gqxz7+qyzi*gqyz7));

    dpwkdx = ci*(sxk*gux2+syk*guy2+szk*guz2)
                        - ck*(sxi*gc5+syi*gc6+szi*gc7)
                      - sxi*(qxxk*gqxx5+qyyk*gqyy5+qzzk*gqzz5
                     +2.0f*(qxyk*gqxy5+qxzk*gqxz5+qyzk*gqyz5))
                      - syi*(qxxk*gqxx6+qyyk*gqyy6+qzzk*gqzz6
                     +2.0f*(qxyk*gqxy6+qxzk*gqxz6+qyzk*gqyz6))
                      - szi*(qxxk*gqxx7+qyyk*gqyy7+qzzk*gqzz7
                     +2.0f*(qxyk*gqxy7+qxzk*gqxz7+qyzk*gqyz7))
                      + sxk*(qxxi*gux11+qyyi*gux14+qzzi*gux16
                     +2.0f*(qxyi*gux12+qxzi*gux13+qyzi*gux15))
                      + syk*(qxxi*guy11+qyyi*guy14+qzzi*guy16
                     +2.0f*(qxyi*guy12+qxzi*guy13+qyzi*guy15))
                      + szk*(qxxi*guz11+qyyi*guz14+qzzi*guz16
                     +2.0f*(qxyi*guz12+qxzi*guz13+qyzi*guz15));

    dpdx = 0.5f * (dpsymdx + 0.5f*(dpwidx + dpwkdx));

    dpsymdy = -uxi*(sxk*gux6+syk*guy6+szk*guz6)
                          - uyi*(sxk*gux8+syk*guy8+szk*guz8)
                          - uzi*(sxk*gux9+syk*guy9+szk*guz9)
                          - uxk*(sxi*gux6+syi*guy6+szi*guz6)
                          - uyk*(sxi*gux8+syi*guy8+szi*guz8)
                          - uzk*(sxi*gux9+syi*guy9+szi*guz9);

    dpwidy = ci*(sxk*gc6+syk*gc8+szk*gc9)
                        - ck*(sxi*gux3+syi*guy3+szi*guz3)
                         - sxi*(qxxk*gux12+qyyk*gux17+qzzk*gux19
                        +2.0f*(qxyk*gux14+qxzk*gux15+qyzk*gux18))
                         - syi*(qxxk*guy12+qyyk*guy17+qzzk*guy19
                        +2.0f*(qxyk*guy14+qxzk*guy15+qyzk*guy18))
                         - szi*(qxxk*guz12+qyyk*guz17+qzzk*guz19
                        +2.0f*(qxyk*guz14+qxzk*guz15+qyzk*guz18))
                         + sxk*(qxxi*gqxx6+qyyi*gqyy6+qzzi*gqzz6
                        +2.0f*(qxyi*gqxy6+qxzi*gqxz6+qyzi*gqyz6))
                         + syk*(qxxi*gqxx8+qyyi*gqyy8+qzzi*gqzz8
                        +2.0f*(qxyi*gqxy8+qxzi*gqxz8+qyzi*gqyz8))
                         + szk*(qxxi*gqxx9+qyyi*gqyy9+qzzi*gqzz9
                        +2.0f*(qxyi*gqxy9+qxzi*gqxz9+qyzi*gqyz9));

    dpwkdy = ci*(sxk*gux3+syk*guy3+szk*guz3)
                        - ck*(sxi*gc6+syi*gc8+szi*gc9)
                      - sxi*(qxxk*gqxx6+qyyk*gqyy6+qzzk*gqzz6
                     +2.0f*(qxyk*gqxy6+qxzk*gqxz6+qyzk*gqyz6))
                      - syi*(qxxk*gqxx8+qyyk*gqyy8+qzzk*gqzz8
                     +2.0f*(qxyk*gqxy8+qxzk*gqxz8+qyzk*gqyz8))
                      - szi*(qxxk*gqxx9+qyyk*gqyy9+qzzk*gqzz9
                     +2.0f*(qxyk*gqxy9+qxzk*gqxz9+qyzk*gqyz9))
                      + sxk*(qxxi*gux12+qyyi*gux17+qzzi*gux19
                     +2.0f*(qxyi*gux14+qxzi*gux15+qyzi*gux18))
                      + syk*(qxxi*guy12+qyyi*guy17+qzzi*guy19
                     +2.0f*(qxyi*guy14+qxzi*guy15+qyzi*guy18))
                      + szk*(qxxi*guz12+qyyi*guz17+qzzi*guz19
                     +2.0f*(qxyi*guz14+qxzi*guz15+qyzi*guz18));

    dpdy    = 0.5f * (dpsymdy + 0.5f*(dpwidy + dpwkdy));

    dpsymdz = -uxi*(sxk*gux7+syk*guy7+szk*guz7)
                          - uyi*(sxk*gux9+syk*guy9+szk*guz9)
                          - uzi*(sxk*gux10+syk*guy10+szk*guz10)
                          - uxk*(sxi*gux7+syi*guy7+szi*guz7)
                          - uyk*(sxi*gux9+syi*guy9+szi*guz9)
                          - uzk*(sxi*gux10+syi*guy10+szi*guz10);

    dpwidz = ci*(sxk*gc7+syk*gc9+szk*gc10)
                        - ck*(sxi*gux4+syi*guy4+szi*guz4)
                      - sxi*(qxxk*gux13+qyyk*gux18+qzzk*gux20
                     +2.0f*(qxyk*gux15+qxzk*gux16+qyzk*gux19))
                      - syi*(qxxk*guy13+qyyk*guy18+qzzk*guy20
                     +2.0f*(qxyk*guy15+qxzk*guy16+qyzk*guy19))
                      - szi*(qxxk*guz13+qyyk*guz18+qzzk*guz20
                     +2.0f*(qxyk*guz15+qxzk*guz16+qyzk*guz19))
                      + sxk*(qxxi*gqxx7+qyyi*gqyy7+qzzi*gqzz7
                     +2.0f*(qxyi*gqxy7+qxzi*gqxz7+qyzi*gqyz7))
                      + syk*(qxxi*gqxx9+qyyi*gqyy9+qzzi*gqzz9
                     +2.0f*(qxyi*gqxy9+qxzi*gqxz9+qyzi*gqyz9))
                      + szk*(qxxi*gqxx10+qyyi*gqyy10+qzzi*gqzz10
                     +2.0f*(qxyi*gqxy10+qxzi*gqxz10+qyzi*gqyz10));

    dpwkdz = ci*(sxk*gux4+syk*guy4+szk*guz4)
                        - ck*(sxi*gc7+syi*gc9+szi*gc10)
                      - sxi*(qxxk*gqxx7+qyyk*gqyy7+qzzk*gqzz7
                     +2.0f*(qxyk*gqxy7+qxzk*gqxz7+qyzk*gqyz7))
                      - syi*(qxxk*gqxx9+qyyk*gqyy9+qzzk*gqzz9
                     +2.0f*(qxyk*gqxy9+qxzk*gqxz9+qyzk*gqyz9))
                      - szi*(qxxk*gqxx10+qyyk*gqyy10+qzzk*gqzz10
                     +2.0f*(qxyk*gqxy10+qxzk*gqxz10+qyzk*gqyz10))
                      + sxk*(qxxi*gux13+qyyi*gux18+qzzi*gux20
                     +2.0f*(qxyi*gux15+qxzi*gux16+qyzi*gux19))
                      + syk*(qxxi*guy13+qyyi*guy18+qzzi*guy20
                     +2.0f*(qxyi*guy15+qxzi*guy16+qyzi*guy19))
                      + szk*(qxxi*guz13+qyyi*guz18+qzzi*guz20
                     +2.0f*(qxyi*guz15+qxzi*guz16+qyzi*guz19));

    dpdz = 0.5f * (dpsymdz + 0.5f*(dpwidz + dpwkdz));

    // effective radii chain rule terms for the;
    // electrostatic solvation free energy gradient of the permanent;
    // multipoles in the reaction potential of the induced dipoles;

    dsymdr = -uxi*(sxk*gux22+syk*guy22+szk*guz22)
                          - uyi*(sxk*gux23+syk*guy23+szk*guz23)
                          - uzi*(sxk*gux24+syk*guy24+szk*guz24)
                          - uxk*(sxi*gux22+syi*guy22+szi*guz22)
                          - uyk*(sxi*gux23+syi*guy23+szi*guz23)
                          - uzk*(sxi*gux24+syi*guy24+szi*guz24);

    dwipdr = ci*(sxk*gc22+syk*gc23+szk*gc24)
                         - ck*(sxi*gux21+syi*guy21+szi*guz21)
                      - sxi*(qxxk*gux25+qyyk*gux28+qzzk*gux30
                     +2.0f*(qxyk*gux26+qxzk*gux27+qyzk*gux29))
                      - syi*(qxxk*guy25+qyyk*guy28+qzzk*guy30
                     +2.0f*(qxyk*guy26+qxzk*guy27+qyzk*guy29))
                      - szi*(qxxk*guz25+qyyk*guz28+qzzk*guz30
                     +2.0f*(qxyk*guz26+qxzk*guz27+qyzk*guz29))
                      + sxk*(qxxi*gqxx22+qyyi*gqyy22+qzzi*gqzz22
                     +2.0f*(qxyi*gqxy22+qxzi*gqxz22+qyzi*gqyz22))
                      + syk*(qxxi*gqxx23+qyyi*gqyy23+qzzi*gqzz23
                     +2.0f*(qxyi*gqxy23+qxzi*gqxz23+qyzi*gqyz23))
                      + szk*(qxxi*gqxx24+qyyi*gqyy24+qzzi*gqzz24
                     +2.0f*(qxyi*gqxy24+qxzi*gqxz24+qyzi*gqyz24));

    dwkpdr = ci*(sxk*gux21+syk*guy21+szk*guz21)
                         - ck*(sxi*gc22+syi*gc23+szi*gc24)
                      - sxi*(qxxk*gqxx22+qyyk*gqyy22+qzzk*gqzz22
                     +2.0f*(qxyk*gqxy22+qxzk*gqxz22+qyzk*gqyz22))
                      - syi*(qxxk*gqxx23+qyyk*gqyy23+qzzk*gqzz23
                     +2.0f*(qxyk*gqxy23+qxzk*gqxz23+qyzk*gqyz23))
                      - szi*(qxxk*gqxx24+qyyk*gqyy24+qzzk*gqzz24
                     +2.0f*(qxyk*gqxy24+qxzk*gqxz24+qyzk*gqyz24))
                      + sxk*(qxxi*gux25+qyyi*gux28+qzzi*gux30
                     +2.0f*(qxyi*gux26+qxzi*gux27+qyzi*gux29))
                      + syk*(qxxi*guy25+qyyi*guy28+qzzi*guy30
                     +2.0f*(qxyi*guy26+qxzi*guy27+qyzi*guy29))
                      + szk*(qxxi*guz25+qyyi*guz28+qzzi*guz30
                     +2.0f*(qxyi*guz26+qxzi*guz27+qyzi*guz29));

    dsumdr = dsymdr + 0.5f*(dwipdr + dwkpdr);
    dpbi = 0.5f*rbk*dsumdr;
    dpbk = 0.5f*rbi*dsumdr;

    // mutual polarization electrostatic solvation free energy gradient

//   if (poltyp .eq. 'MUTUAL'){

        dpdx = dpdx - 0.5f *
                           (dxi*(pxk*gux5+pyk*gux6+pzk*gux7)
                           +dyi*(pxk*guy5+pyk*guy6+pzk*guy7)
                           +dzi*(pxk*guz5+pyk*guz6+pzk*guz7)
                           +dxk*(pxi*gux5+pyi*gux6+pzi*gux7)
                           +dyk*(pxi*guy5+pyi*guy6+pzi*guy7)
                           +dzk*(pxi*guz5+pyi*guz6+pzi*guz7));

        dpdy = dpdy - 0.5f *
                           (dxi*(pxk*gux6+pyk*gux8+pzk*gux9)
                           +dyi*(pxk*guy6+pyk*guy8+pzk*guy9)
                           +dzi*(pxk*guz6+pyk*guz8+pzk*guz9)
                           +dxk*(pxi*gux6+pyi*gux8+pzi*gux9)
                           +dyk*(pxi*guy6+pyi*guy8+pzi*guy9)
                           +dzk*(pxi*guz6+pyi*guz8+pzi*guz9));

        dpdz = dpdz - 0.5f *
                           (dxi*(pxk*gux7+pyk*gux9+pzk*gux10)
                           +dyi*(pxk*guy7+pyk*guy9+pzk*guy10)
                           +dzi*(pxk*guz7+pyk*guz9+pzk*guz10)
                           +dxk*(pxi*gux7+pyi*gux9+pzi*gux10)
                           +dyk*(pxi*guy7+pyi*guy9+pzi*guy10)
                           +dzk*(pxi*guz7+pyi*guz9+pzi*guz10));

        duvdr = dxi*(pxk*gux22+pyk*gux23+pzk*gux24)
                            + dyi*(pxk*guy22+pyk*guy23+pzk*guy24)
                            + dzi*(pxk*guz22+pyk*guz23+pzk*guz24)
                            + dxk*(pxi*gux22+pyi*gux23+pzi*gux24)
                            + dyk*(pxi*guy22+pyi*guy23+pzi*guy24)
                            + dzk*(pxi*guz22+pyi*guz23+pzi*guz24);
        dpbi = dpbi - 0.5f*rbk*duvdr;
        dpbk = dpbk - 0.5f*rbi*duvdr;
//    }

    // torque due to induced reaction field on permanent dipoles

    float fid1 = 0.5f * (sxk*gux2+syk*guy2+szk*guz2);
    float fid2 = 0.5f * (sxk*gux3+syk*guy3+szk*guz3);
    float fid3 = 0.5f * (sxk*gux4+syk*guy4+szk*guz4);
    float fkd1 = 0.5f * (sxi*gux2+syi*guy2+szi*guz2);
    float fkd2 = 0.5f * (sxi*gux3+syi*guy3+szi*guz3);
    float fkd3 = 0.5f * (sxi*gux4+syi*guy4+szi*guz4);


    // the factor 0.5 appears to be included since trqi1[i] & trqi1[k]
    // are identical in the Tinker code (inner loop starts at k = i
    // factor not needed here since
/*
    if ( sameAtom )
    {
        fid1 = 0.5f * fid1;
        fid2 = 0.5f * fid2;
        fid3 = 0.5f * fid3;
        fkd1 = 0.5f * fkd1;
        fkd2 = 0.5f * fkd2;
        fkd3 = 0.5f * fkd3;
    }
*/
    float trqi1   = uyi*fid3 - uzi*fid2;
    float trqi2   = uzi*fid1 - uxi*fid3;
    float trqi3   = uxi*fid2 - uyi*fid1;

    float trqi_k1 = uyk*fkd3 - uzk*fkd2;
    float trqi_k2 = uzk*fkd1 - uxk*fkd3;
    float trqi_k3 = uxk*fkd2 - uyk*fkd1;


    // torque due to induced reaction field gradient on quadrupoles;

    float fidg11 = -0.25f *
                              ( (sxk*gqxx2+syk*gqxx3+szk*gqxx4)
                              + (sxk*gux5+syk*guy5+szk*guz5));

    float fidg12 = -0.25f *
                              ( (sxk*gqxy2+syk*gqxy3+szk*gqxy4)
                              + (sxk*gux6+syk*guy6+szk*guz6));

    float fidg13 = -0.25f *
                              ( (sxk*gqxz2+syk*gqxz3+szk*gqxz4)
                              + (sxk*gux7+syk*guy7+szk*guz7));

    float fidg22 = -0.25f *
                              ( (sxk*gqyy2+syk*gqyy3+szk*gqyy4)
                              + (sxk*gux8+syk*guy8+szk*guz8));

    float fidg23 = -0.25f *
                              ( (sxk*gqyz2+syk*gqyz3+szk*gqyz4)
                              + (sxk*gux9+syk*guy9+szk*guz9));

    float fidg33 = -0.25f *
                              ( (sxk*gqzz2+syk*gqzz3+szk*gqzz4)
                              + (sxk*gux10+syk*guy10+szk*guz10));

    float fidg21 = fidg12;
    float fidg31 = fidg13;
    float fidg32 = fidg23;

    float fkdg11 = 0.25f *
                              ( (sxi*gqxx2+syi*gqxx3+szi*gqxx4)
                              + (sxi*gux5+syi*guy5+szi*guz5));

    float fkdg12 = 0.25f *
                              ( (sxi*gqxy2+syi*gqxy3+szi*gqxy4)
                              + (sxi*gux6+syi*guy6+szi*guz6));
    float fkdg13 = 0.25f *
                              ( (sxi*gqxz2+syi*gqxz3+szi*gqxz4)
                              + (sxi*gux7+syi*guy7+szi*guz7));
    float fkdg22 = 0.25f *
                              ( (sxi*gqyy2+syi*gqyy3+szi*gqyy4)
                              + (sxi*gux8+syi*guy8+szi*guz8));
    float fkdg23 = 0.25f *
                              ( (sxi*gqyz2+syi*gqyz3+szi*gqyz4)
                              + (sxi*gux9+syi*guy9+szi*guz9));
    float fkdg33 = 0.25f *
                              ( (sxi*gqzz2+syi*gqzz3+szi*gqzz4)
                              + (sxi*gux10+syi*guy10+szi*guz10));
    float fkdg21 = fkdg12;
    float fkdg31 = fkdg13;
    float fkdg32 = fkdg23;

/*
    if ( sameAtom )
    {
        fidg11 = 0.5f * fidg11;
        fidg12 = 0.5f * fidg12;
        fidg13 = 0.5f * fidg13;
        fidg21 = 0.5f * fidg21;
        fidg22 = 0.5f * fidg22;
        fidg23 = 0.5f * fidg23;
        fidg31 = 0.5f * fidg31;
        fidg32 = 0.5f * fidg32;
        fidg33 = 0.5f * fidg33;
        fkdg11 = 0.5f * fkdg11;
        fkdg12 = 0.5f * fkdg12;
        fkdg13 = 0.5f * fkdg13;
        fkdg21 = 0.5f * fkdg21;
        fkdg22 = 0.5f * fkdg22;
        fkdg23 = 0.5f * fkdg23;
        fkdg31 = 0.5f * fkdg31;
        fkdg32 = 0.5f * fkdg32;
        fkdg33 = 0.5f * fkdg33;
     }
*/

    trqi1 += 2.0f*(qxyi*fidg13+qyyi*fidg23+qyzi*fidg33
                        -qxzi*fidg12-qyzi*fidg22-qzzi*fidg32);
    trqi2 += 2.0f*(qxzi*fidg11+qyzi*fidg21+qzzi*fidg31
                        -qxxi*fidg13-qxyi*fidg23-qxzi*fidg33);

    trqi3 += 2.0f*(qxxi*fidg12+qxyi*fidg22+qxzi*fidg32
                        -qxyi*fidg11-qyyi*fidg21-qyzi*fidg31);

    trqi_k1 += 2.0f*
                        (qxyk*fkdg13+qyyk*fkdg23+qyzk*fkdg33
                        -qxzk*fkdg12-qyzk*fkdg22-qzzk*fkdg32);

    trqi_k2 += 2.0f*
                        (qxzk*fkdg11+qyzk*fkdg21+qzzk*fkdg31
                        -qxxk*fkdg13-qxyk*fkdg23-qxzk*fkdg33);

    trqi_k3 += 2.0f*
                        (qxxk*fkdg12+qxyk*fkdg22+qxzk*fkdg32
                        -qxyk*fkdg11-qyyk*fkdg21-qyzk*fkdg31);

    // total permanent and induced energies for this interaction;

    e                        = esym + 0.5f*(ewi+ewk);
    ei                       = 0.5f * (esymi + 0.5f*(ewii+ewki));

    outputForce[0]           = (dedx + dpdx);
    outputForce[1]           = (dedy + dpdy);
    outputForce[2]           = (dedz + dpdz);

    outputTorque[0][0]       = (trq1 + trqi1);
    outputTorque[0][1]       = (trq2 + trqi2);
    outputTorque[0][2]       = (trq3 + trqi3);

    outputTorque[1][0]       = (trq_k1 + trqi_k1);
    outputTorque[1][1]       = (trq_k2 + trqi_k2);
    outputTorque[1][2]       = (trq_k3 + trqi_k3);

    outputBorn[0]            = drbi;
    outputBorn[1]            = drbk;

    outputBornPolar[0]       = dpbi;
    outputBornPolar[1]       = dpbk;

    *outputEnergy            = (e + ei);

}

__device__ void zeroKirkwoodParticleSharedField( struct KirkwoodParticle* sA )
{
    // zero shared fields

    sA->force[0]              = 0.0f;
    sA->force[1]              = 0.0f;
    sA->force[2]              = 0.0f;

    sA->torque[0]             = 0.0f;
    sA->torque[1]             = 0.0f;
    sA->torque[2]             = 0.0f;

    sA->dBornRadius           = 0.0f;
    sA->dBornRadiusPolar      = 0.0f;

}

// Include versions of the kernels for N^2 calculations.

#undef USE_OUTPUT_BUFFER_PER_WARP
#define METHOD_NAME(a, b) a##N2##b
#include "kCalculateAmoebaCudaKirkwood.h"
#define USE_OUTPUT_BUFFER_PER_WARP
#undef METHOD_NAME
#define METHOD_NAME(a, b) a##N2ByWarp##b
#include "kCalculateAmoebaCudaKirkwood.h"

// reduce psWorkArray_3_1 -> force
// reduce psWorkArray_3_2 -> torque

static void kReduceForceTorque(amoebaGpuContext amoebaGpu )
{

    kReduceFields_kernel<<<amoebaGpu->nonbondBlocks, amoebaGpu->fieldReduceThreadsPerBlock>>>(
                           amoebaGpu->paddedNumberOfAtoms*3, amoebaGpu->outputBuffers,
                           amoebaGpu->psWorkArray_3_1->_pDevStream[0], amoebaGpu->psKirkwoodForce->_pDevStream[0] );
    LAUNCHERROR("kReduceForceTorque1");

    kReduceFields_kernel<<<amoebaGpu->nonbondBlocks, amoebaGpu->fieldReduceThreadsPerBlock>>>(
                           amoebaGpu->paddedNumberOfAtoms*3, amoebaGpu->outputBuffers,
                           amoebaGpu->psWorkArray_3_2->_pDevStream[0], amoebaGpu->psTorque->_pDevStream[0] );
    LAUNCHERROR("kReduceForceTorque2");
}

// reduce psWorkArray_1_1 -> dBorn
// reduce psWorkArray_1_2 -> dBornPolar

#ifdef AMOEBA_DEBUG
static void kReduce_dBorn(amoebaGpuContext amoebaGpu )
{

/*
    kReduceFields_kernel<<<amoebaGpu->nonbondBlocks, amoebaGpu->fieldReduceThreadsPerBlock>>>(
                           amoebaGpu->paddedNumberOfAtoms, amoebaGpu->outputBuffers,
                           amoebaGpu->psWorkArray_1_1->_pDevStream[0], amoebaGpu->psBorn->_pDevStream[0] );
    LAUNCHERROR("kReduce_dBorn1");
*/

    kReduceFields_kernel<<<amoebaGpu->nonbondBlocks, amoebaGpu->fieldReduceThreadsPerBlock>>>(
                           amoebaGpu->paddedNumberOfAtoms, amoebaGpu->outputBuffers,
                           amoebaGpu->psWorkArray_1_2->_pDevStream[0], amoebaGpu->psBornPolar->_pDevStream[0] );

    LAUNCHERROR("kReduce_dBorn2");
}
#endif

__global__
#if (__CUDA_ARCH__ >= 200)
__launch_bounds__(GF1XX_THREADS_PER_BLOCK, 1)
#elif (__CUDA_ARCH__ >= 130)
__launch_bounds__(GT2XX_THREADS_PER_BLOCK, 1)
#else
__launch_bounds__(G8X_THREADS_PER_BLOCK, 1)
#endif
 void kReduceToBornForcePrefactor_kernel( unsigned int fieldComponents, unsigned int outputBuffers, float* fieldIn1, float* fieldIn2,
                                                    float* fieldOut )
{
    unsigned int pos = blockIdx.x * blockDim.x + threadIdx.x;

    // Reduce field

    while (pos < fieldComponents)
    {

        float totalField = 0.0f;

        float* pFt1      = fieldIn1 + pos;
        float* pFt2      = fieldIn2 + pos;

        float bornRadius = cSim.pBornRadii[pos];
        float obcChain   = cSim.pObcChain[pos];

        unsigned int i   = outputBuffers;

        while (i >= 4)
        {
            totalField += pFt1[0] + pFt1[fieldComponents] + pFt1[2*fieldComponents] + pFt1[3*fieldComponents];
            totalField += pFt2[0] + pFt2[fieldComponents] + pFt2[2*fieldComponents] + pFt2[3*fieldComponents];
            pFt1       += fieldComponents*4;
            pFt2       += fieldComponents*4;
            i          -= 4;
        }

        if (i >= 2)
        {
            totalField += pFt1[0] + pFt1[fieldComponents];
            totalField += pFt2[0] + pFt2[fieldComponents];
            pFt1       += fieldComponents*2;
            pFt2       += fieldComponents*2;
            i          -= 2;
        }

        if (i > 0)
        {
            totalField += pFt1[0];
            totalField += pFt2[0];
        }

        fieldOut[pos]   = totalField*bornRadius*bornRadius*obcChain;
        pos            += gridDim.x * blockDim.x;
    }
}

__global__
#if (__CUDA_ARCH__ >= 200)
__launch_bounds__(GF1XX_THREADS_PER_BLOCK, 1)
#elif (__CUDA_ARCH__ >= 130)
__launch_bounds__(GT2XX_THREADS_PER_BLOCK, 1)
#else
__launch_bounds__(G8X_THREADS_PER_BLOCK, 1)
#endif
void kReduceToBornForcePrefactorAndSASA_kernel( unsigned int fieldComponents, unsigned int outputBuffers, float* fieldIn1, float* fieldIn2,
                                                           float* fieldOut )
{
    unsigned int pos = blockIdx.x * blockDim.x + threadIdx.x;

    float energy     = 0.0f;

    // Reduce field

    while (pos < fieldComponents)
    {

        float totalForce = 0.0f;

        float* pFt1      = fieldIn1 + pos;
        float* pFt2      = fieldIn2 + pos;

        float bornRadius = cSim.pBornRadii[pos];
        float obcChain   = cSim.pObcChain[pos];
        float2 obcData   = cSim.pObcData[pos];

        unsigned int i   = outputBuffers;

        while (i >= 4)
        {
            totalForce += pFt1[0] + pFt1[fieldComponents] + pFt1[2*fieldComponents] + pFt1[3*fieldComponents];
            totalForce += pFt2[0] + pFt2[fieldComponents] + pFt2[2*fieldComponents] + pFt2[3*fieldComponents];
            pFt1       += fieldComponents*4;
            pFt2       += fieldComponents*4;
            i          -= 4;
        }

        if (i >= 2)
        {
            totalForce += pFt1[0] + pFt1[fieldComponents];
            totalForce += pFt2[0] + pFt2[fieldComponents];
            pFt1       += fieldComponents*2;
            pFt2       += fieldComponents*2;
            i          -= 2;
        }

        if (i > 0)
        {
            totalForce += pFt1[0];
            totalForce += pFt2[0];
        }

        float r        = (obcData.x + cSim.dielectricOffset + cSim.probeRadius);
        float ratio6   = ( (obcData.x + cSim.dielectricOffset) / bornRadius);
              ratio6   = ratio6*ratio6*ratio6;
              ratio6   = ratio6*ratio6;
        float saTerm   = cSim.surfaceAreaFactor * r * r * ratio6;

        totalForce     += saTerm / bornRadius;
        totalForce     *= bornRadius * bornRadius * obcChain;

        energy         += saTerm;

        fieldOut[pos]   = totalForce*bornRadius*bornRadius*obcChain;
        pos            += gridDim.x * blockDim.x;
    }

    cSim.pEnergy[blockIdx.x * blockDim.x + threadIdx.x] += energy / -6.0f;
}

/*
static void kReduceAndCombine_dBorn(amoebaGpuContext amoebaGpu )
{

    kReduceAndCombineFields_kernel<<<amoebaGpu->nonbondBlocks, amoebaGpu->fieldReduceThreadsPerBlock>>>(
                                     amoebaGpu->paddedNumberOfAtoms, amoebaGpu->outputBuffers,
                                     amoebaGpu->psWorkArray_1_1->_pDevStream[0],
                                     amoebaGpu->psWorkArray_1_2->_pDevStream[0],
                                     amoebaGpu->psBorn->_pDevStream[0] );
    LAUNCHERROR("kReduce_dBorn");
} */

static void kReduceToBornForcePrefactor( amoebaGpuContext amoebaGpu )
{

    if( amoebaGpu->includeObcCavityTerm ){
        kReduceToBornForcePrefactorAndSASA_kernel<<<amoebaGpu->nonbondBlocks, amoebaGpu->fieldReduceThreadsPerBlock>>>(
                                                    amoebaGpu->paddedNumberOfAtoms, amoebaGpu->outputBuffers,
                                                    amoebaGpu->psWorkArray_1_1->_pDevStream[0],
                                                    amoebaGpu->psWorkArray_1_2->_pDevStream[0],
                                                    amoebaGpu->gpuContext->psBornForce->_pDevStream[0] );
    } else {
        kReduceToBornForcePrefactor_kernel<<<amoebaGpu->nonbondBlocks, amoebaGpu->fieldReduceThreadsPerBlock>>>(
                                             amoebaGpu->paddedNumberOfAtoms, amoebaGpu->outputBuffers,
                                             amoebaGpu->psWorkArray_1_1->_pDevStream[0],
                                             amoebaGpu->psWorkArray_1_2->_pDevStream[0],
                                             amoebaGpu->gpuContext->psBornForce->_pDevStream[0] );
    }
    LAUNCHERROR("kReduceToBornForcePrefactor");
}

//#ifdef AMOEBA_DEBUG
#if 0
static void printKirkwoodBuffer( amoebaGpuContext amoebaGpu, unsigned int bufferIndex )
{
    (void) fprintf( amoebaGpu->log, "Kirkwood Buffer %u\n", bufferIndex );
    unsigned int start = bufferIndex*3*amoebaGpu->paddedNumberOfAtoms;
    unsigned int stop  = (bufferIndex+1)*3*amoebaGpu->paddedNumberOfAtoms;
    for( unsigned int ii = start; ii < stop; ii += 3 ){
        unsigned int ii3Index      = ii/3;
        unsigned int bufferIndex   = ii3Index/(amoebaGpu->paddedNumberOfAtoms);
        unsigned int particleIndex = ii3Index - bufferIndex*(amoebaGpu->paddedNumberOfAtoms);
        (void) fprintf( amoebaGpu->log, "   %6u %3u %6u [%14.6e %14.6e %14.6e] [%14.6e %14.6e %14.6e]\n",
                            ii/3,  bufferIndex, particleIndex,
                            amoebaGpu->psWorkArray_3_1->_pSysStream[0][ii],
                            amoebaGpu->psWorkArray_3_1->_pSysStream[0][ii+1],
                            amoebaGpu->psWorkArray_3_1->_pSysStream[0][ii+2],
                            amoebaGpu->psWorkArray_3_2->_pSysStream[0][ii],
                            amoebaGpu->psWorkArray_3_2->_pSysStream[0][ii+1],
                            amoebaGpu->psWorkArray_3_2->_pSysStream[0][ii+2] );
    }

/*
    start = 0;
    stop  = -146016;
    float maxV = -1.0e+99;
    for( unsigned int ii = start; ii < stop; ii += 3 ){
        if(  amoebaGpu->psWorkArray_3_1->_pSysStream[0][ii] > maxV ){
            unsigned int ii3Index      = ii/3;
            unsigned int bufferIndex   = ii3Index/(amoebaGpu->paddedNumberOfAtoms);
            unsigned int particleIndex = ii3Index - bufferIndex*(amoebaGpu->paddedNumberOfAtoms);
            (void) fprintf( amoebaGpu->log, "MaxQ %6u %3u %6u %14.6e\n",
                            ii/3,  bufferIndex, particleIndex,
                            amoebaGpu->psWorkArray_3_1->_pSysStream[0][ii] );
            maxV = amoebaGpu->psWorkArray_3_1->_pSysStream[0][ii];
        }
    }
*/
}

static void printKirkwoodAtomBuffers( amoebaGpuContext amoebaGpu, unsigned int targetAtom )
{
    (void) fprintf( amoebaGpu->log, "Kirkwood atom %u\n", targetAtom );
    for( unsigned int ii = 0; ii < amoebaGpu->outputBuffers; ii++ ){
        unsigned int particleIndex = 3*(targetAtom + ii*amoebaGpu->paddedNumberOfAtoms);
        (void) fprintf( amoebaGpu->log, " %2u %6u [%14.6e %14.6e %14.6e] [%14.6e %14.6e %14.6e]\n",
                        ii, particleIndex,
                        amoebaGpu->psWorkArray_3_1->_pSysStream[0][particleIndex],
                        amoebaGpu->psWorkArray_3_1->_pSysStream[0][particleIndex+1],
                        amoebaGpu->psWorkArray_3_1->_pSysStream[0][particleIndex+2],
                        amoebaGpu->psWorkArray_3_2->_pSysStream[0][particleIndex],
                        amoebaGpu->psWorkArray_3_2->_pSysStream[0][particleIndex+1],
                        amoebaGpu->psWorkArray_3_2->_pSysStream[0][particleIndex+2] );
    }
}
#endif

/**---------------------------------------------------------------------------------------

   Compute Amoeba electrostatic force & torque

   @param amoebaGpu        amoebaGpu context
   @param gpu              OpenMM gpu Cuda context

   --------------------------------------------------------------------------------------- */

void kCalculateAmoebaKirkwood( amoebaGpuContext amoebaGpu )
{

   // ---------------------------------------------------------------------------------------

    static unsigned int threadsPerBlock = 0;

#ifdef AMOEBA_DEBUG
    static const char* methodName       = "kCalculateAmoebaKirkwood";
    static int timestep = 0;
    std::vector<int> fileId;
    timestep++;
    fileId.resize( 2 );
    fileId[0] = timestep;
    fileId[1] = 1;
#endif

    // ---------------------------------------------------------------------------------------

    gpuContext gpu = amoebaGpu->gpuContext;

    // apparently debug array can take up nontrivial no. registers

#ifdef AMOEBA_DEBUG
    if( amoebaGpu->log ){
        (void) fprintf( amoebaGpu->log, "%s %d maxCovalentDegreeSz=%d"
                        " gamma=%.3e scalingDistanceCutoff=%.3f ZZZ\n",
                        methodName, gpu->natoms,
                        amoebaGpu->maxCovalentDegreeSz, amoebaGpu->pGamma,
                        amoebaGpu->scalingDistanceCutoff );
    }
   int paddedNumberOfAtoms                    = amoebaGpu->gpuContext->sim.paddedNumberOfAtoms;
    CUDAStream<float4>* debugArray            = new CUDAStream<float4>(paddedNumberOfAtoms*paddedNumberOfAtoms, 1, "DebugArray");
    memset( debugArray->_pSysStream[0],      0, sizeof( float )*4*paddedNumberOfAtoms*paddedNumberOfAtoms);
    debugArray->Upload();
    unsigned int targetAtom                   = 0;

    gpu->psBornRadii->Download();
    (void) fprintf( amoebaGpu->log, "Kirkwood input\n" ); (void) fflush( amoebaGpu->log );
    for( int ii = 0; ii < amoebaGpu->gpuContext->sim.paddedNumberOfAtoms; ii++ ){
        (void) fprintf( amoebaGpu->log,"Born %6d %16.9e\n", ii,
                        gpu->psBornRadii->_pSysStream[0][ii] );
    }
#endif

    // on first pass, set threads/block and based on that setting the energy buffer array

    if( threadsPerBlock == 0 ){
        unsigned int maxThreads;
        if (gpu->sm_version >= SM_20)
            maxThreads = 256;
        else if (gpu->sm_version >= SM_12)
            maxThreads = 128;
        else
            maxThreads = 64;
        threadsPerBlock = std::min(getThreadsPerBlock(amoebaGpu, sizeof(KirkwoodParticle)), maxThreads);
        //unsigned int eDiffhreadsPerBlock            = getThreadsPerBlock( amoebaGpu, sizeof(KirkwoodEDiffParticle));
        //unsigned int maxThreadsPerBlock             = threadsPerBlock> eDiffhreadsPerBlock ? threadsPerBlock : eDiffhreadsPerBlock;

        if( amoebaGpu->log ){

            (void) fprintf( amoebaGpu->log, "kCalculateAmoebaCudaKirkwood: blcks=%u tds=%u %u bPrWrp=%u atm=%lu shrd=%lu Ebuf=%u ixnCt=%lu workUnits=%u\n",
                            amoebaGpu->nonbondBlocks, threadsPerBlock, maxThreads, amoebaGpu->bOutputBufferPerWarp,
                            sizeof(KirkwoodParticle), sizeof(KirkwoodParticle)*threadsPerBlock,
                            amoebaGpu->energyOutputBuffers, (*gpu->psInteractionCount)[0], gpu->sim.workUnits );
            (void) fflush( amoebaGpu->log );
        }
    }

    kClearFields_1( amoebaGpu );
    kClearFields_3( amoebaGpu, 6 );

    if (gpu->bOutputBufferPerWarp){
    } else {

#ifdef AMOEBA_DEBUG
        (void) fprintf( amoebaGpu->log, "kCalculateAmoebaCudaKirkwoodN2Forces no warp:  numBlocks=%u numThreads=%u bufferPerWarp=%u atm=%u shrd=%u Ebuf=%u ixnCt=%u workUnits=%u\n",
                        amoebaGpu->nonbondBlocks, threadsPerBlock, amoebaGpu->bOutputBufferPerWarp,
                        sizeof(KirkwoodParticle), sizeof(KirkwoodParticle)*threadsPerBlock,
                        amoebaGpu->energyOutputBuffers, (*gpu->psInteractionCount)[0], gpu->sim.workUnits );
        (void) fflush( amoebaGpu->log );
#endif

        kCalculateAmoebaCudaKirkwoodN2Forces_kernel<<<amoebaGpu->nonbondBlocks, threadsPerBlock, sizeof(KirkwoodParticle)*threadsPerBlock>>>(
                                                                           amoebaGpu->psWorkUnit->_pDevStream[0]
#ifdef AMOEBA_DEBUG
                                                                           , debugArray->_pDevStream[0], targetAtom );
#else
                                                                           );
#endif
    }

    LAUNCHERROR("kCalculateAmoebaCudaKirkwoodN2Forces");

#ifdef AMOEBA_DEBUG
    if( amoebaGpu->log ){

        amoebaGpu->psWorkArray_3_1->Download();
        amoebaGpu->psWorkArray_3_2->Download();
        amoebaGpu->psWorkArray_1_1->Download();
        amoebaGpu->psWorkArray_1_2->Download();
/*
        amoebaGpu->psLabFrameDipole->Download();
        amoebaGpu->psLabFrameQuadrupole->Download();
        amoebaGpu->psInducedDipoleS->Download();
        amoebaGpu->psInducedDipolePolarS->Download();

        for( int ii = 0; ii < gpu->natoms; ii++ ){
           int indexOffset3    = ii*3;
           int indexOffset9    = ii*9;
           (void) fprintf( amoebaGpu->log, "%5d [%14.7e %14.7e %14.7e] q[%14.7e %14.7e %14.7e]\n", ii,
                           amoebaGpu->psLabFrameDipole->_pSysStream[0][indexOffset3],
                           amoebaGpu->psLabFrameDipole->_pSysStream[0][indexOffset3+1],
                           amoebaGpu->psLabFrameDipole->_pSysStream[0][indexOffset3+2],
                           amoebaGpu->psLabFrameQuadrupole->_pSysStream[0][indexOffset9],
                           amoebaGpu->psLabFrameQuadrupole->_pSysStream[0][indexOffset9+1],
                           amoebaGpu->psLabFrameQuadrupole->_pSysStream[0][indexOffset9+2] );
           (void) fprintf( amoebaGpu->log, "%5d [%14.7e %14.7e %14.7e] q[%14.7e %14.7e %14.7e]\n", ii,
                           amoebaGpu->psInducedDipoleS->_pSysStream[0][indexOffset3],
                           amoebaGpu->psInducedDipoleS->_pSysStream[0][indexOffset3+1],
                           amoebaGpu->psInducedDipoleS->_pSysStream[0][indexOffset3+2],
                           amoebaGpu->psInducedDipolePolarS->_pSysStream[0][indexOffset3],
                           amoebaGpu->psInducedDipolePolarS->_pSysStream[0][indexOffset3+1],
                           amoebaGpu->psInducedDipolePolarS->_pSysStream[0][indexOffset3+2] );
        }
*/
        debugArray->Download();

        (void) fprintf( amoebaGpu->log, "Target Info\n" );
        (void) fflush( amoebaGpu->log );

        int paddedNumberOfAtoms          = amoebaGpu->gpuContext->sim.paddedNumberOfAtoms;
        for( int jj = 0; jj < gpu->natoms; jj++ ){
            int debugIndex = jj;
            (void) fprintf( amoebaGpu->log,"%5d %5d DebugGk\n", targetAtom, jj );
            for( int kk = 0; kk < 8; kk++ ){
                (void) fprintf( amoebaGpu->log,"[%16.9e %16.9e %16.9e %16.9e]\n",
                                debugArray->_pSysStream[0][debugIndex].x, debugArray->_pSysStream[0][debugIndex].y,
                                debugArray->_pSysStream[0][debugIndex].z, debugArray->_pSysStream[0][debugIndex].w );
                debugIndex += paddedNumberOfAtoms;
            }
            (void) fprintf( amoebaGpu->log,"\n" );
       }
       //printKirkwoodAtomBuffers( amoebaGpu, (targetAtom + 0) );
       //printKirkwoodAtomBuffers( amoebaGpu, (targetAtom + 1231) );
       //printKirkwoodBuffer( amoebaGpu, 0 );
       //printKirkwoodBuffer( amoebaGpu, 38 );
    }
#endif

    kReduceForceTorque( amoebaGpu );
    kReduceToBornForcePrefactor( amoebaGpu );

#ifdef AMOEBA_DEBUG
    if( amoebaGpu->log ){

        kReduce_dBorn( amoebaGpu );
        amoebaGpu->psKirkwoodForce->Download();
        amoebaGpu->psTorque->Download();
        amoebaGpu->psBorn->Download();
        amoebaGpu->psBornPolar->Download();
        gpu->psBornRadii->Download();
        gpu->psObcChain->Download();
        debugArray->Download();

        int maxPrint        = 10;
        for( int ii = 0; ii < gpu->natoms; ii++ ){
           (void) fprintf( amoebaGpu->log, "%5d ", ii);

            int indexOffset     = ii*3;

           // force

           (void) fprintf( amoebaGpu->log,"KirkwoodF [%16.9e %16.9e %16.9e] ",
                           amoebaGpu->psKirkwoodForce->_pSysStream[0][indexOffset],
                           amoebaGpu->psKirkwoodForce->_pSysStream[0][indexOffset+1],
                           amoebaGpu->psKirkwoodForce->_pSysStream[0][indexOffset+2] );

           // torque

           (void) fprintf( amoebaGpu->log,"T [%16.9e %16.9e %16.9e] ",
                           amoebaGpu->psTorque->_pSysStream[0][indexOffset],
                           amoebaGpu->psTorque->_pSysStream[0][indexOffset+1],
                           amoebaGpu->psTorque->_pSysStream[0][indexOffset+2] );

           // d_Born

           //float bornForceValue = amoebaGpu->psBorn->_pSysStream[0][ii];
           float bornForceValue = gpu->psBornForce->_pSysStream[0][ii];
           float bornRadius     = gpu->psBornRadii->_pSysStream[0][ii];
           float obcChain       = gpu->psObcChain->_pSysStream[0][ii];
           float bornSumValue   = bornRadius*obcChain != 0.0f ? bornForceValue/(bornRadius*bornRadius*obcChain) : 0.0f;
           float bornValue      = bornSumValue - amoebaGpu->psBornPolar->_pSysStream[0][ii];
           (void) fprintf( amoebaGpu->log,"dB br=%16.9e obcC=%16.9e bSum=%16.9e  [%16.9e %16.9e]",
                           bornRadius,obcChain, bornSumValue, bornValue,
                           amoebaGpu->psBornPolar->_pSysStream[0][ii] );

           // coords

           (void) fprintf( amoebaGpu->log,"\n" );
           if( ii == maxPrint && (gpu->natoms - maxPrint) > ii ){
                ii = gpu->natoms - maxPrint;
           }
        }
/*
        //kReduceAndCombine_dBorn( amoebaGpu );
        amoebaGpu->psBorn->Download();
        for( int ii = 0; ii < gpu->natoms; ii++ ){
           (void) fprintf( amoebaGpu->log, "%5d ", ii);

           // d_Born

           (void) fprintf( amoebaGpu->log,"dBrnSum %16.9e ",
                           amoebaGpu->psBorn->_pSysStream[0][ii] );

           (void) fprintf( amoebaGpu->log,"\n" );
           if( ii == maxPrint && (gpu->natoms - maxPrint) > ii ){
                ii = gpu->natoms - maxPrint;
           }
        }
*/
        (void) fflush( amoebaGpu->log );

        if( 1 ){
            std::vector<int> fileId;
            //fileId.push_back( 0 );
            VectorOfDoubleVectors outputVector;
            cudaLoadCudaFloat4Array( gpu->natoms, 3, gpu->psPosq4,                    outputVector );
            cudaLoadCudaFloatArray( gpu->natoms,  3, amoebaGpu->psKirkwoodForce,      outputVector );
            cudaLoadCudaFloatArray( gpu->natoms,  3, amoebaGpu->psTorque,             outputVector);
            cudaWriteVectorOfDoubleVectorsToFile( "CudaForceTorque", fileId, outputVector );

         }

    }
    delete debugArray;
#endif

    // map torques to forces

    cudaComputeAmoebaMapTorquesAndAddTotalForce( amoebaGpu, amoebaGpu->psTorque, amoebaGpu->psKirkwoodForce, amoebaGpu->gpuContext->psForce4 );

#ifdef AMOEBA_DEBUG
    if( amoebaGpu->log ){

        cudaComputeAmoebaMapTorques( amoebaGpu, amoebaGpu->psTorque, amoebaGpu->psKirkwoodForce );
        amoebaGpu->psKirkwoodForce->Download();
        (void) fprintf( amoebaGpu->log, "Mapped Kirkwood torques to forces.\n" ); (void) fflush( amoebaGpu->log );

        int maxPrint        = 10;
        for( int ii = 0; ii < gpu->natoms; ii++ ){
           (void) fprintf( amoebaGpu->log, "%5d ", ii);

            int indexOffset     = ii*3;

           // force

           (void) fprintf( amoebaGpu->log,"KirkwoodF [%16.9e %16.9e %16.9e] ",
                           amoebaGpu->psKirkwoodForce->_pSysStream[0][indexOffset],
                           amoebaGpu->psKirkwoodForce->_pSysStream[0][indexOffset+1],
                           amoebaGpu->psKirkwoodForce->_pSysStream[0][indexOffset+2] );

           (void) fprintf( amoebaGpu->log,"\n" );
           if( ii == maxPrint && (gpu->natoms - maxPrint) > ii ){
                ii = gpu->natoms - maxPrint;
           }
        }
        (void) fflush( amoebaGpu->log );
        if( 1 ){
            std::vector<int> fileId;
            //fileId.push_back( 0 );
            VectorOfDoubleVectors outputVector;
            cudaLoadCudaFloat4Array( gpu->natoms, 3, gpu->psPosq4,                    outputVector );
            cudaLoadCudaFloatArray( gpu->natoms,  3, amoebaGpu->psKirkwoodForce,      outputVector );
            cudaWriteVectorOfDoubleVectorsToFile( "CudaKirkwoodForce", fileId, outputVector );
         }

    }
#endif

   // Tinker's Born1

   //kClearForces(amoebaGpu->gpuContext );
   //kCalculateAmoebaObcGbsaForces2( amoebaGpu );
   kCalculateObcGbsaForces2( amoebaGpu->gpuContext );

   // E-diff

   kCalculateAmoebaKirkwoodEDiff( amoebaGpu );

   // ---------------------------------------------------------------------------------------
}
