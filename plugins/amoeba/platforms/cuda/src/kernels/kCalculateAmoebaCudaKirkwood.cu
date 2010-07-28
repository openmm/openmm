//-----------------------------------------------------------------------------------------

//-----------------------------------------------------------------------------------------

#include "amoebaGpuTypes.h"
#include "amoebaCudaKernels.h"
#include "kCalculateAmoebaCudaUtilities.h"
#ifdef _MSC_VER
extern void kCalculateObcGbsaForces2(gpuContext gpu);
#else
#include "cudaKernels.h"
#endif
#include "kCalculateAmoebaCudaKirkwoodParticle.h"

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
    float a[6][4];
    float b[5][3];
    float fid[4],fkd[4];
    float fidg[4][4],fkdg[4][4];
    float trq[4];
    float trqi[4];
    float trq_k[4];
    float trqi_k[4];
    float gc[31];
    float gux[31],guy[31],guz[31];
    float gqxx[31],gqxy[31];
    float gqxz[31],gqyy[31];
    float gqyz[31],gqzz[31];

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

    a[0][0]      =            gf;
    a[1][0]      =          -gf3;
    a[2][0]      =    3.0f * gf5;
    a[3][0]      =  -15.0f * gf7;
    a[4][0]      =  105.0f * gf9;
    a[5][0]      = -945.0f * gf11;
 
    // Born radii derivatives of reaction potential auxiliary terms;
                 
    b[0][0]      = dgfdr * a[1][0];
    b[1][0]      = dgfdr * a[2][0];
    b[2][0]      = dgfdr * a[3][0];
    b[3][0]      = dgfdr * a[4][0];
    b[4][0]      = dgfdr * a[5][0];

    // reaction potential gradient auxiliary terms;

    expc1        = 1.0f - expc;
    a[0][1]      = expc1 * a[1][0];
    a[1][1]      = expc1 * a[2][0];
    a[2][1]      = expc1 * a[3][0];
    a[3][1]      = expc1 * a[4][0];
    a[4][1]      = expc1 * a[5][0];

    // Born radii derivs of reaction potential gradient auxiliary terms;

    b[0][1]      = b[1][0] - expcr*a[1][0] - expc*b[1][0];
    b[1][1]      = b[2][0] - expcr*a[2][0] - expc*b[2][0];
    b[2][1]      = b[3][0] - expcr*a[3][0] - expc*b[3][0];
    b[3][1]      = b[4][0] - expcr*a[4][0] - expc*b[4][0];
 
    // 2nd reaction potential gradient auxiliary terms;
 
    expcdexpc    = -expc * dexpc;
    a[0][2]      = expc1*a[1][1] + expcdexpc*a[1][0];
    a[1][2]      = expc1*a[2][1] + expcdexpc*a[2][0];
    a[2][2]      = expc1*a[3][1] + expcdexpc*a[3][0];
    a[3][2]      = expc1*a[4][1] + expcdexpc*a[4][0];
 
    // Born radii derivatives of the 2nd reaction potential
    // gradient auxiliary terms

     b[0][2]     = b[1][1] - (expcr*(a[1][1] + dexpc*a[1][0])
                           + expc*(b[1][1] + dexpcr*a[1][0]+dexpc*b[1][0]));

     b[1][2]     = b[2][1] - (expcr*(a[2][1] + dexpc*a[2][0])
                           +   expc*(b[2][1] + dexpcr*a[2][0]+dexpc*b[2][0]));

     b[2][2]     = b[3][1] - (expcr*(a[3][1] + dexpc*a[3][0])
                           +   expc*(b[3][1] + dexpcr*a[3][0]+dexpc*b[3][0]));
 
    // 3rd reaction potential gradient auxiliary terms;
 
    expcdexpc    = 2.0f * expcdexpc;
    a[0][3]      = expc1*a[1][2] + expcdexpc*a[1][1];
    a[1][3]      = expc1*a[2][2] + expcdexpc*a[2][1];
    a[2][3]      = expc1*a[3][2] + expcdexpc*a[3][1];

    //expcdexpc    = -expc * powf( dexpc, 2.0f );
    expcdexpc    = -expc*dexpc*dexpc;
    a[0][3]      = a[0][3] + expcdexpc*a[1][0];
    a[1][3]      = a[1][3] + expcdexpc*a[2][0];
    a[2][3]      = a[2][3] + expcdexpc*a[3][0];
 
    // multiply the auxillary terms by their dieletric functions;
 
    a[0][0]      = fc * a[0][0];
    a[0][1]      = fc * a[0][1];
    a[0][2]      = fc * a[0][2];
    a[0][3]      = fc * a[0][3];
    b[0][0]      = fc * b[0][0];
    b[0][1]      = fc * b[0][1];
    b[0][2]      = fc * b[0][2];
    a[1][0]      = fd * a[1][0];
    a[1][1]      = fd * a[1][1];
    a[1][2]      = fd * a[1][2];
    a[1][3]      = fd * a[1][3];
    b[1][0]      = fd * b[1][0];
    b[1][1]      = fd * b[1][1];
    b[1][2]      = fd * b[1][2];
    a[2][0]      = fq * a[2][0];
    a[2][1]      = fq * a[2][1];
    a[2][2]      = fq * a[2][2];
    a[2][3]      = fq * a[2][3];
    b[2][0]      = fq * b[2][0];
    b[2][1]      = fq * b[2][1];
    b[2][2]      = fq * b[2][2];
 
    // unweighted reaction potential tensor;

    gc[1]        = a[0][0];
    gux[1]       = xr * a[1][0];
    guy[1]       = yr * a[1][0];
    guz[1]       = zr * a[1][0];
    gqxx[1]      = xr2 * a[2][0];
    gqyy[1]      = yr2 * a[2][0];
    gqzz[1]      = zr2 * a[2][0];
    gqxy[1]      = xr * yr * a[2][0];
    gqxz[1]      = xr * zr * a[2][0];
    gqyz[1]      = yr * zr * a[2][0];
 
    // Born radii derivs of unweighted reaction potential tensor;

    gc[21]       = b[0][0];
    gux[21]      = xr * b[1][0];
    guy[21]      = yr * b[1][0];
    guz[21]      = zr * b[1][0];
    gqxx[21]     = xr2 * b[2][0];
    gqyy[21]     = yr2 * b[2][0];
    gqzz[21]     = zr2 * b[2][0];
    gqxy[21]     = xr * yr * b[2][0];
    gqxz[21]     = xr * zr * b[2][0];
    gqyz[21]     = yr * zr * b[2][0];

    // unweighted reaction potential gradient tensor;

    gc[2]        = xr * a[0][1];
    gc[3]        = yr * a[0][1];
    gc[4]        = zr * a[0][1];
    gux[2]       = a[1][0] + xr2*a[1][1];
    gux[3]       = xr * yr * a[1][1];
    gux[4]       = xr * zr * a[1][1];

    guy[2]       = gux[3];
    guy[3]       = a[1][0] + yr2*a[1][1];
    guy[4]       = yr * zr * a[1][1];
    guz[2]       = gux[4];
    guz[3]       = guy[4];
    guz[4]       = a[1][0] + zr2*a[1][1];
    gqxx[2]      = xr * (2.0f*a[2][0]+xr2*a[2][1]);
    gqxx[3]      = yr * xr2 * a[2][1];
    gqxx[4]      = zr * xr2 * a[2][1];
    gqyy[2]      = xr * yr2 * a[2][1];
    gqyy[3]      = yr * (2.0f*a[2][0]+yr2*a[2][1]);
    gqyy[4]      = zr * yr2 * a[2][1];
    gqzz[2]      = xr * zr2 * a[2][1];
    gqzz[3]      = yr * zr2 * a[2][1];
    gqzz[4]      = zr * (2.0f*a[2][0]+zr2*a[2][1]);
    gqxy[2]      = yr * (a[2][0]+xr2*a[2][1]);
    gqxy[3]      = xr * (a[2][0]+yr2*a[2][1]);
    gqxy[4]      = zr * xr * yr * a[2][1];
    gqxz[2]      = zr * (a[2][0]+xr2*a[2][1]);
    gqxz[3]      = gqxy[4];
    gqxz[4]      = xr * (a[2][0]+zr2*a[2][1]);
    gqyz[2]      = gqxy[4];
    gqyz[3]      = zr * (a[2][0]+yr2*a[2][1]);
    gqyz[4]      = yr * (a[2][0]+zr2*a[2][1]);

    // Born derivs of the unweighted reaction potential gradient tensor;
 
    gc[22]       = xr * b[0][1];
    gc[22]       = xr * b[0][1];
    gc[23]       = yr * b[0][1];
    gc[24]       = zr * b[0][1];
    gux[22]      = b[1][0] + xr2*b[1][1];
    gux[23]      = xr * yr * b[1][1];
    gux[24]      = xr * zr * b[1][1];
    guy[22]      = gux[23];
    guy[23]      = b[1][0] + yr2*b[1][1];
    guy[24]      = yr * zr * b[1][1];
    guz[22]      = gux[24];
    guz[23]      = guy[24];
    guz[24]      = b[1][0] + zr2*b[1][1];
    gqxx[22]     = xr * (2.0f*b[2][0]+xr2*b[2][1]);
    gqxx[23]     = yr * xr2 * b[2][1];
    gqxx[24]     = zr * xr2 * b[2][1];
    gqyy[22]     = xr * yr2 * b[2][1];
    gqyy[23]     = yr * (2.0f*b[2][0]+yr2*b[2][1]);
    gqyy[24]     = zr * yr2 * b[2][1];
    gqzz[22]     = xr * zr2 * b[2][1];
    gqzz[23]     = yr * zr2 * b[2][1];
    gqzz[24]     = zr * (2.0f*b[2][0] + zr2*b[2][1]);
    gqxy[22]     = yr * (b[2][0]+xr2*b[2][1]);
    gqxy[23]     = xr * (b[2][0]+yr2*b[2][1]);
    gqxy[24]     = zr * xr * yr * b[2][1];
    gqxz[22]     = zr * (b[2][0]+xr2*b[2][1]);
    gqxz[23]     = gqxy[24];
    gqxz[24]     = xr * (b[2][0]+zr2*b[2][1]);
    gqyz[22]     = gqxy[24];
    gqyz[23]     = zr * (b[2][0]+yr2*b[2][1]);
    gqyz[24]     = yr * (b[2][0]+zr2*b[2][1]);
 
    // unweighted 2nd reaction potential gradient tensor;
 
    gc[5]        = a[0][1] + xr2*a[0][2];
    gc[6]        = xr * yr * a[0][2];
    gc[7]        = xr * zr * a[0][2];
    gc[8]        = a[0][1] + yr2*a[0][2];
    gc[9]        = yr * zr * a[0][2];
    gc[10]       = a[0][1] + zr2*a[0][2];
    gux[5]       = xr * (3.0f*a[1][1]+xr2*a[1][2]);
    gux[6]       = yr * (a[1][1]+xr2*a[1][2]);
    gux[7]       = zr * (a[1][1]+xr2*a[1][2]);
    gux[8]       = xr * (a[1][1]+yr2*a[1][2]);
    gux[9]       = zr * xr * yr * a[1][2];
    gux[10]      = xr * (a[1][1]+zr2*a[1][2]);
    guy[5]       = yr * (a[1][1]+xr2*a[1][2]);
    guy[6]       = xr * (a[1][1]+yr2*a[1][2]);
    guy[7]       = gux[9];
    guy[8]       = yr * (3.0f*a[1][1]+yr2*a[1][2]);
    guy[9]       = zr * (a[1][1]+yr2*a[1][2]);
    guy[10]      = yr * (a[1][1]+zr2*a[1][2]);
    guz[5]       = zr * (a[1][1]+xr2*a[1][2]);
    guz[6]       = gux[9];
    guz[7]       = xr * (a[1][1]+zr2*a[1][2]);
    guz[8]       = zr * (a[1][1]+yr2*a[1][2]);
    guz[9]       = yr * (a[1][1]+zr2*a[1][2]);
    guz[10]      = zr * (3.0f*a[1][1]+zr2*a[1][2]);
    gqxx[5]      = 2.0f*a[2][0] + xr2*(5.0f*a[2][1]+xr2*a[2][2]);
    gqxx[6]      = yr * xr * (2.0f*a[2][1]+xr2*a[2][2]);
    gqxx[7]      = zr * xr * (2.0f*a[2][1]+xr2*a[2][2]);
    gqxx[8]      = xr2 * (a[2][1]+yr2*a[2][2]);
    gqxx[9]      = zr * yr * xr2 * a[2][2];
    gqxx[10]     = xr2 * (a[2][1]+zr2*a[2][2]);
    gqyy[5]      = yr2 * (a[2][1]+xr2*a[2][2]);
    gqyy[6]      = xr * yr * (2.0f*a[2][1]+yr2*a[2][2]);
    gqyy[7]      = xr * zr * yr2 * a[2][2];
    gqyy[8]      = 2.0f*a[2][0] + yr2*(5.0f*a[2][1]+yr2*a[2][2]);
    gqyy[9]      = yr * zr * (2.0f*a[2][1]+yr2*a[2][2]);
    gqyy[10]     = yr2 * (a[2][1]+zr2*a[2][2]);
    gqzz[5]      = zr2 * (a[2][1]+xr2*a[2][2]);
    gqzz[6]      = xr * yr * zr2 * a[2][2];
    gqzz[7]      = xr * zr * (2.0f*a[2][1]+zr2*a[2][2]);
    gqzz[8]      = zr2 * (a[2][1]+yr2*a[2][2]);
    gqzz[9]      = yr * zr * (2.0f*a[2][1]+zr2*a[2][2]);
    gqzz[10]     = 2.0f*a[2][0] + zr2*(5.0f*a[2][1]+zr2*a[2][2]);
    gqxy[5]      = xr * yr * (3.0f*a[2][1]+xr2*a[2][2]);
    gqxy[6]      = a[2][0] + (xr2+yr2)*a[2][1] + xr2*yr2*a[2][2];
    gqxy[7]      = zr * yr * (a[2][1]+xr2*a[2][2]);
    gqxy[8]      = xr * yr * (3.0f*a[2][1]+yr2*a[2][2]);
    gqxy[9]      = zr * xr * (a[2][1]+yr2*a[2][2]);
    gqxy[10]     = xr * yr * (a[2][1]+zr2*a[2][2]);
    gqxz[5]      = xr * zr * (3.0f*a[2][1]+xr2*a[2][2]);
    gqxz[6]      = yr * zr * (a[2][1]+xr2*a[2][2]);
    gqxz[7]      = a[2][0] + (xr2+zr2)*a[2][1] + xr2*zr2*a[2][2];
    gqxz[8]      = xr * zr * (a[2][1]+yr2*a[2][2]);
    gqxz[9]      = xr * yr * (a[2][1]+zr2*a[2][2]);
    gqxz[10]     = xr * zr * (3.0f*a[2][1]+zr2*a[2][2]);
    gqyz[5]      = zr * yr * (a[2][1]+xr2*a[2][2]);
    gqyz[6]      = xr * zr * (a[2][1]+yr2*a[2][2]);
    gqyz[7]      = xr * yr * (a[2][1]+zr2*a[2][2]);
    gqyz[8]      = yr * zr * (3.0f*a[2][1]+yr2*a[2][2]);
    gqyz[9]      = a[2][0] + (yr2+zr2)*a[2][1] + yr2*zr2*a[2][2];
    gqyz[10]     = yr * zr * (3.0f*a[2][1]+zr2*a[2][2]);
 
    // Born radii derivatives of the unweighted 2nd reaction;
    // potential gradient tensor;

    gc[25]       = b[0][1] + xr2*b[0][2];
    gc[26]       = xr * yr * b[0][2];
    gc[27]       = xr * zr * b[0][2];
    gc[28]       = b[0][1] + yr2*b[0][2];
    gc[29]       = yr * zr * b[0][2];
    gc[30]       = b[0][1] + zr2*b[0][2];
    gux[25]      = xr * (3.0f*b[1][1]+xr2*b[1][2]);
    gux[26]      = yr * (b[1][1]+xr2*b[1][2]);
    gux[27]      = zr * (b[1][1]+xr2*b[1][2]);
    gux[28]      = xr * (b[1][1]+yr2*b[1][2]);
    gux[29]      = zr * xr * yr * b[1][2];
    gux[30]      = xr * (b[1][1]+zr2*b[1][2]);
    guy[25]      = yr * (b[1][1]+xr2*b[1][2]);
    guy[26]      = xr * (b[1][1]+yr2*b[1][2]);
    guy[27]      = gux[29];
    guy[28]      = yr * (3.0f*b[1][1]+yr2*b[1][2]);
    guy[29]      = zr * (b[1][1]+yr2*b[1][2]);
    guy[30]      = yr * (b[1][1]+zr2*b[1][2]);
    guz[25]      = zr * (b[1][1]+xr2*b[1][2]);
    guz[26]      = gux[29];
    guz[27]      = xr * (b[1][1]+zr2*b[1][2]);
    guz[28]      = zr * (b[1][1]+yr2*b[1][2]);
    guz[29]      = yr * (b[1][1]+zr2*b[1][2]);
    guz[30]      = zr * (3.0f*b[1][1]+zr2*b[1][2]);
    gqxx[25]     = 2.0f*b[2][0] + xr2*(5.0f*b[2][1]+xr2*b[2][2]);
    gqxx[26]     = yr * xr * (2.0f*b[2][1]+xr2*b[2][2]);
    gqxx[27]     = zr * xr * (2.0f*b[2][1]+xr2*b[2][2]);
    gqxx[28]     = xr2 * (b[2][1]+yr2*b[2][2]);
    gqxx[29]     = zr * yr * xr2 * b[2][2];
    gqxx[30]     = xr2 * (b[2][1]+zr2*b[2][2]);
    gqyy[25]     = yr2 * (b[2][1]+xr2*b[2][2]);
    gqyy[26]     = xr * yr * (2.0f*b[2][1]+yr2*b[2][2]);
    gqyy[27]     = xr * zr * yr2 * b[2][2];
    gqyy[28]     = 2.0f*b[2][0] + yr2*(5.0f*b[2][1]+yr2*b[2][2]);
    gqyy[29]     = yr * zr * (2.0f*b[2][1]+yr2*b[2][2]);
    gqyy[30]     = yr2 * (b[2][1]+zr2*b[2][2]);
    gqzz[25]     = zr2 * (b[2][1]+xr2*b[2][2]);
    gqzz[26]     = xr * yr * zr2 * b[2][2];
    gqzz[27]     = xr * zr * (2.0f*b[2][1]+zr2*b[2][2]);
    gqzz[28]     = zr2 * (b[2][1]+yr2*b[2][2]);
    gqzz[29]     = yr * zr * (2.0f*b[2][1]+zr2*b[2][2]);
    gqzz[30]     = 2.0f*b[2][0] + zr2*(5.0f*b[2][1]+zr2*b[2][2]);
    gqxy[25]     = xr * yr * (3.0f*b[2][1] + xr2*b[2][2]);
    gqxy[26]     = b[2][0] + (xr2+yr2)*b[2][1] + xr2*yr2*b[2][2];
    gqxy[27]     = zr * yr * (b[2][1]+xr2*b[2][2]);
    gqxy[28]     = xr * yr * (3.0f*b[2][1]+yr2*b[2][2]);
    gqxy[29]     = zr * xr * (b[2][1]+yr2*b[2][2]);
    gqxy[30]     = xr * yr * (b[2][1]+zr2*b[2][2]);
    gqxz[25]     = xr * zr * (3.0f*b[2][1]+xr2*b[2][2]);
    gqxz[26]     = yr * zr * (b[2][1]+xr2*b[2][2]);
    gqxz[27]     = b[2][0] + (xr2+zr2)*b[2][1] + xr2*zr2*b[2][2];
    gqxz[28]     = xr * zr * (b[2][1]+yr2*b[2][2]);
    gqxz[29]     = xr * yr * (b[2][1]+zr2*b[2][2]);
    gqxz[30]     = xr * zr * (3.0f*b[2][1]+zr2*b[2][2]);
    gqyz[25]     = zr * yr * (b[2][1]+xr2*b[2][2]);
    gqyz[26]     = xr * zr * (b[2][1]+yr2*b[2][2]);
    gqyz[27]     = xr * yr * (b[2][1]+zr2*b[2][2]);
    gqyz[28]     = yr * zr * (3.0f*b[2][1]+yr2*b[2][2]);
    gqyz[29]     = b[2][0] + (yr2+zr2)*b[2][1] + yr2*zr2*b[2][2];
    gqyz[30]     = yr * zr * (3.0f*b[2][1]+zr2*b[2][2]);

    // unweighted 3rd reaction potential gradient tensor;

    gc[11]       = xr * (3.0f*a[0][2]+xr2*a[0][3]);
    gc[12]       = yr * (a[0][2]+xr2*a[0][3]);
    gc[13]       = zr * (a[0][2]+xr2*a[0][3]);
    gc[14]       = xr * (a[0][2]+yr2*a[0][3]);
    gc[15]       = xr * yr * zr * a[0][3];
    gc[16]       = xr * (a[0][2]+zr2*a[0][3]);
    gc[17]       = yr * (3.0f*a[0][2]+yr2*a[0][3]);
    gc[18]       = zr * (a[0][2]+yr2*a[0][3]);
    gc[19]       = yr * (a[0][2]+zr2*a[0][3]);
    gc[20]       = zr * (3.0f*a[0][2]+zr2*a[0][3]);
    gux[11]      = 3.0f*a[1][1] + xr2*(6.0f*a[1][2]+xr2*a[1][3]);
    gux[12]      = xr * yr * (3.0f*a[1][2]+xr2*a[1][3]);
    gux[13]      = xr * zr * (3.0f*a[1][2]+xr2*a[1][3]);
    gux[14]      = a[1][1] + (xr2+yr2)*a[1][2] + xr2*yr2*a[1][3];
    gux[15]      = yr * zr * (a[1][2]+xr2*a[1][3]);
    gux[16]      = a[1][1] + (xr2+zr2)*a[1][2] + xr2*zr2*a[1][3];
    gux[17]      = xr * yr * (3.0f*a[1][2]+yr2*a[1][3]);
    gux[18]      = xr * zr * (a[1][2]+yr2*a[1][3]);
    gux[19]      = xr * yr * (a[1][2]+zr2*a[1][3]);
    gux[20]      = xr * zr * (3.0f*a[1][2]+zr2*a[1][3]);
    guy[11]      = gux[12];
    guy[12]      = gux[14];
    guy[13]      = gux[15];
    guy[14]      = gux[17];
    guy[15]      = gux[18];
    guy[16]      = gux[19];
    guy[17]      = 3.0f*a[1][1] + yr2*(6.0f*a[1][2]+yr2*a[1][3]);
    guy[18]      = yr * zr * (3.0f*a[1][2]+yr2*a[1][3]);
    guy[19]      = a[1][1] + (yr2+zr2)*a[1][2] + yr2*zr2*a[1][3];
    guy[20]      = yr * zr * (3.0f*a[1][2]+zr2*a[1][3]);
    guz[11]      = gux[13];
    guz[12]      = gux[15];
    guz[13]      = gux[16];
    guz[14]      = gux[18];
    guz[15]      = gux[19];
    guz[16]      = gux[20];
    guz[17]      = guy[18];
    guz[18]      = guy[19];
    guz[19]      = guy[20];
    guz[20]      = 3.0f*a[1][1] + zr2*(6.0f*a[1][2]+zr2*a[1][3]);

    gqxx[11]     = xr * (12.0f*a[2][1]+xr2*(9.0f*a[2][2] + xr2*a[2][3]));
    gqxx[12]     = yr * (2.0f*a[2][1]+xr2*(5.0f*a[2][2]  + xr2*a[2][3]));
    gqxx[13]     = zr * (2.0f*a[2][1]+xr2*(5.0f*a[2][2]  + xr2*a[2][3]));
    gqxx[14]     = xr * (2.0f*a[2][1]+yr2*2.0f*a[2][2]   +xr2*(a[2][2]+yr2*a[2][3]));
    gqxx[15]     = xr * yr * zr * (2.0f*a[2][2]+xr2*a[2][3]);
    gqxx[16]     = xr * (2.0f*a[2][1]+zr2*2.0f*a[2][2] +xr2*(a[2][2]+zr2*a[2][3]));
    gqxx[17]     = yr * xr2 * (3.0f*a[2][2]+yr2*a[2][3]);
    gqxx[18]     = zr * xr2 * (a[2][2]+yr2*a[2][3]);
    gqxx[19]     = yr * xr2 * (a[2][2]+zr2*a[2][3]);
    gqxx[20]     = zr * xr2 * (3.0f*a[2][2]+zr2*a[2][3]);
    gqxy[11]     = yr * (3.0f*a[2][1]+xr2*(6.0f*a[2][2] +xr2*a[2][3]));
    gqxy[12]     = xr * (3.0f*(a[2][1]+yr2*a[2][2]) +xr2*(a[2][2]+yr2*a[2][3]));
    gqxy[13]     = xr * yr * zr * (3.0f*a[2][2]+xr2*a[2][3]);
    gqxy[14]     = yr * (3.0f*(a[2][1]+xr2*a[2][2]) +yr2*(a[2][2]+xr2*a[2][3]));
    gqxy[15]     = zr * (a[2][1]+(yr2+xr2)*a[2][2] +yr2*xr2*a[2][3]);
    gqxy[16]     = yr * (a[2][1]+(xr2+zr2)*a[2][2] +xr2*zr2*a[2][3]);
    gqxy[17]     = xr * (3.0f*(a[2][1]+yr2*a[2][2]) +yr2*(3.0f*a[2][2]+yr2*a[2][3]));
    gqxy[18]     = xr * yr * zr * (3.0f*a[2][2]+yr2*a[2][3]);
    gqxy[19]     = xr * (a[2][1]+(yr2+zr2)*a[2][2] +yr2*zr2*a[2][3]);
    gqxy[20]     = xr * yr * zr * (3.0f*a[2][2]+zr2*a[2][3]);
    gqxz[11]     = zr * (3.0f*a[2][1]+xr2*(6.0f*a[2][2] +xr2*a[2][3]));
    gqxz[12]     = xr * yr * zr * (3.0f*a[2][2]+xr2*a[2][3]);
    gqxz[13]     = xr * (3.0f*(a[2][1]+zr2*a[2][2]) +xr2*(a[2][2]+zr2*a[2][3]));
    gqxz[14]     = zr * (a[2][1]+(xr2+yr2)*a[2][2] +xr2*yr2*a[2][3]);
    gqxz[15]     = yr * (a[2][1]+(xr2+zr2)*a[2][2] +zr2*xr2*a[2][3]);
    gqxz[16]     = zr * (3.0f*(a[2][1]+xr2*a[2][2]) +zr2*(a[2][2]+xr2*a[2][3]));
    gqxz[17]     = xr * yr * zr * (3.0f*a[2][2]+yr2*a[2][3]);
    gqxz[18]     = xr * (a[2][1]+(zr2+yr2)*a[2][2] +zr2*yr2*a[2][3]);
    gqxz[19]     = xr * yr * zr * (3.0f*a[2][2]+zr2*a[2][3]);
    gqxz[20]     = xr * (3.0f*a[2][1]+zr2*(6.0f*a[2][2] +zr2*a[2][3]));
    gqyy[11]     = xr * yr2 * (3.0f*a[2][2]+xr2*a[2][3]);
    gqyy[12]     = yr * (2.0f*a[2][1]+xr2*2.0f*a[2][2] +yr2*(a[2][2]+xr2*a[2][3]));
    gqyy[13]     = zr * yr2 * (a[2][2]+xr2*a[2][3]);
    gqyy[14]     = xr * (2.0f*a[2][1]+yr2*(5.0f*a[2][2] +yr2*a[2][3]));
    gqyy[15]     = xr * yr * zr * (2.0f*a[2][2]+yr2*a[2][3]);
    gqyy[16]     = xr * yr2 * (a[2][2]+zr2*a[2][3]);
    gqyy[17]     = yr * (12.0f*a[2][1]+yr2*(9.0f*a[2][2] +yr2*a[2][3]));
    gqyy[18]     = zr * (2.0f*a[2][1]+yr2*(5.0f*a[2][2] +yr2*a[2][3]));
    gqyy[19]     = yr * (2.0f*a[2][1]+zr2*2.0f*a[2][2] +yr2*(a[2][2]+zr2*a[2][3]));
    gqyy[20]     = zr * yr2 * (3.0f*a[2][2]+zr2*a[2][3]);
    gqyz[11]     = xr * yr * zr * (3.0f*a[2][2]+xr2*a[2][3]);
    gqyz[12]     = zr * (a[2][1]+(xr2+yr2)*a[2][2] +xr2*yr2*a[2][3]);
    gqyz[13]     = yr * (a[2][1]+(xr2+zr2)*a[2][2] +xr2*zr2*a[2][3]);
    gqyz[14]     = xr * yr * zr * (3.0f*a[2][2]+yr2*a[2][3]);
    gqyz[15]     = xr * (a[2][1]+(yr2+zr2)*a[2][2] +yr2*zr2*a[2][3]);
    gqyz[16]     = xr * yr * zr * (3.0f*a[2][2]+zr2*a[2][3]);
    gqyz[17]     = zr * (3.0f*a[2][1]+yr2*(6.0f*a[2][2] +yr2*a[2][3]));
    gqyz[18]     = yr * (3.0f*(a[2][1]+zr2*a[2][2]) +yr2*(a[2][2]+zr2*a[2][3]));
    gqyz[19]     = zr * (3.0f*(a[2][1]+yr2*a[2][2]) +zr2*(a[2][2]+yr2*a[2][3]));
    gqyz[20]     = yr * (3.0f*a[2][1]+zr2*(6.0f*a[2][2] +zr2*a[2][3]));
    gqzz[11]     = xr * zr2 * (3.0f*a[2][2]+xr2*a[2][3]);
    gqzz[12]     = yr * (zr2*a[2][2]+xr2*(zr2*a[2][3]));
    gqzz[13]     = zr * (2.0f*a[2][1]+xr2*2.0f*a[2][2] +zr2*(a[2][2]+xr2*a[2][3]));
    gqzz[14]     = xr * zr2 * (a[2][2]+yr2*a[2][3]);
    gqzz[15]     = xr * yr * zr * (2.0f*a[2][2]+zr2*a[2][3]);
    gqzz[16]     = xr * (2.0f*a[2][1]+zr2*(5.0f*a[2][2] +zr2*a[2][3]));
    gqzz[17]     = yr * zr2 * (3.0f*a[2][2]+yr2*a[2][3]);
    gqzz[18]     = zr * (2.0f*a[2][1]+yr2*2.0f*a[2][2] +zr2*(a[2][2]+yr2*a[2][3]));
    gqzz[19]     = yr * (2.0f*a[2][1]+zr2*(5.0f*a[2][2] +zr2*a[2][3]));
    gqzz[20]     = zr * (12.0f*a[2][1]+zr2*(9.0f*a[2][2] +zr2*a[2][3]));
 
    // electrostatic solvation energy of the permanent multipoles
    // in their own GK reaction potential
 
    esym = ci * ck * gc[1] - (uxi*(uxk*gux[2]+uyk*guy[2]+uzk*guz[2])
                           +  uyi*(uxk*gux[3]+uyk*guy[3]+uzk*guz[3])
                           +  uzi*(uxk*gux[4]+uyk*guy[4]+uzk*guz[4]));

    ewi =  ci*(uxk*gc[2]+uyk*gc[3]+uzk*gc[4])
          -ck*(uxi*gux[1]+uyi*guy[1]+uzi*guz[1])
           +ci*(qxxk*gc[5]+qyyk*gc[8]+qzzk*gc[10]
              +2.0f*(qxyk*gc[6]+qxzk*gc[7]+qyzk*gc[9]))
                 +ck*(qxxi*gqxx[1]+qyyi*gqyy[1]+qzzi*gqzz[1]
              +2.0f*(qxyi*gqxy[1]+qxzi*gqxz[1]+qyzi*gqyz[1]))
               - uxi*(qxxk*gux[5]+qyyk*gux[8]+qzzk*gux[10]
              +2.0f*(qxyk*gux[6]+qxzk*gux[7]+qyzk*gux[9]))
               - uyi*(qxxk*guy[5]+qyyk*guy[8]+qzzk*guy[10]
              +2.0f*(qxyk*guy[6]+qxzk*guy[7]+qyzk*guy[9]))
               - uzi*(qxxk*guz[5]+qyyk*guz[8]+qzzk*guz[10]
              +2.0f*(qxyk*guz[6]+qxzk*guz[7]+qyzk*guz[9]))
               + uxk*(qxxi*gqxx[2]+qyyi*gqyy[2]+qzzi*gqzz[2]
              +2.0f*(qxyi*gqxy[2]+qxzi*gqxz[2]+qyzi*gqyz[2]))
               + uyk*(qxxi*gqxx[3]+qyyi*gqyy[3]+qzzi*gqzz[3]
              +2.0f*(qxyi*gqxy[3]+qxzi*gqxz[3]+qyzi*gqyz[3]))
               + uzk*(qxxi*gqxx[4]+qyyi*gqyy[4]+qzzi*gqzz[4]
              +2.0f*(qxyi*gqxy[4]+qxzi*gqxz[4]+qyzi*gqyz[4]))
              + qxxi*(qxxk*gqxx[5]+qyyk*gqxx[8]+qzzk*gqxx[10]
              +2.0f*(qxyk*gqxx[6]+qxzk*gqxx[7]+qyzk*gqxx[9]))
              + qyyi*(qxxk*gqyy[5]+qyyk*gqyy[8]+qzzk*gqyy[10]
              +2.0f*(qxyk*gqyy[6]+qxzk*gqyy[7]+qyzk*gqyy[9]))
              + qzzi*(qxxk*gqzz[5]+qyyk*gqzz[8]+qzzk*gqzz[10]
              +2.0f*(qxyk*gqzz[6]+qxzk*gqzz[7]+qyzk*gqzz[9]))
              + 2.0f*(qxyi*(qxxk*gqxy[5]+qyyk*gqxy[8]+qzzk*gqxy[10]
              +2.0f*(qxyk*gqxy[6]+qxzk*gqxy[7]+qyzk*gqxy[9]))
              + qxzi*(qxxk*gqxz[5]+qyyk*gqxz[8]+qzzk*gqxz[10]
              +2.0f*(qxyk*gqxz[6]+qxzk*gqxz[7]+qyzk*gqxz[9]))
              + qyzi*(qxxk*gqyz[5]+qyyk*gqyz[8]+qzzk*gqyz[10]
              +2.0f*(qxyk*gqyz[6]+qxzk*gqyz[7]+qyzk*gqyz[9])));

#ifdef AMOEBA_DEBUG
//int n2         = numberOfAtoms*numberOfAtoms;
//int debugIndex = atomI*numberOfAtoms + atomJ;
#if 0
if( atomI == -1 ){
// 4, [3,3], [5,5], [5,5],[5,5],
int debugIndex               = atomJ;

    debugArray[debugIndex].x = gkc;
    debugArray[debugIndex].y = fc;
    debugArray[debugIndex].z = gf2;

    debugIndex              += numberOfAtoms;
    debugArray[debugIndex].x = ci;
    debugArray[debugIndex].y = ck;
    debugArray[debugIndex].z = esym;
    debugArray[debugIndex].w = ewi;

    debugIndex              += numberOfAtoms;
    debugArray[debugIndex].x = uxi;
    debugArray[debugIndex].y = uyi;
    debugArray[debugIndex].z = uzi;

    debugIndex              += numberOfAtoms;
    debugArray[debugIndex].x = uxk;
    debugArray[debugIndex].y = uyk;
    debugArray[debugIndex].z = uzk;

    debugIndex              += numberOfAtoms;
    debugArray[debugIndex].x = gc[1];
    debugArray[debugIndex].y = gc[2];
    debugArray[debugIndex].z = gc[3];
    debugArray[debugIndex].w = gc[4];

    debugIndex              += numberOfAtoms;
    debugArray[debugIndex].x = gc[5];
    debugArray[debugIndex].y = gc[6];
    debugArray[debugIndex].z = gc[7];
    debugArray[debugIndex].w = gc[8];

    debugIndex              += numberOfAtoms;
    debugArray[debugIndex].x = gc[9];
    debugArray[debugIndex].y = gc[10];


    debugIndex              += numberOfAtoms;
    debugArray[debugIndex].x = gux[1];
    debugArray[debugIndex].y = gux[2];
    debugArray[debugIndex].z = gux[3];
    debugArray[debugIndex].w = gux[4];

    debugIndex              += numberOfAtoms;
    debugArray[debugIndex].x = gux[5];
    debugArray[debugIndex].y = gux[6];
    debugArray[debugIndex].z = gux[7];
    debugArray[debugIndex].w = gux[8];

    debugIndex              += numberOfAtoms;
    debugArray[debugIndex].x = gux[9];
    debugArray[debugIndex].y = gux[10];

    debugIndex              += numberOfAtoms;
    debugArray[debugIndex].x = guy[1];
    debugArray[debugIndex].y = guy[2];
    debugArray[debugIndex].z = guy[3];
    debugArray[debugIndex].w = guy[4];

    debugIndex              += numberOfAtoms;
    debugArray[debugIndex].x = guy[5];
    debugArray[debugIndex].y = guy[6];
    debugArray[debugIndex].z = guy[7];
    debugArray[debugIndex].w = guy[8];

    debugIndex              += numberOfAtoms;
    debugArray[debugIndex].x = guy[9];
    debugArray[debugIndex].y = guy[10];

    debugIndex              += numberOfAtoms;
    debugArray[debugIndex].x = guz[1];
    debugArray[debugIndex].y = guz[2];
    debugArray[debugIndex].z = guz[3];
    debugArray[debugIndex].w = guz[4];

    debugIndex              += numberOfAtoms;
    debugArray[debugIndex].x = guz[5];
    debugArray[debugIndex].y = guz[6];
    debugArray[debugIndex].z = guz[7];
    debugArray[debugIndex].w = guz[8];

    debugIndex              += numberOfAtoms;
    debugArray[debugIndex].x = guz[9];
    debugArray[debugIndex].y = guz[10];

    debugIndex              += numberOfAtoms;
    debugArray[debugIndex].x = qxxi;
    debugArray[debugIndex].y = qyyi;
    debugArray[debugIndex].z = qzzi;

    debugIndex              += numberOfAtoms;
    debugArray[debugIndex].x = qxxk;
    debugArray[debugIndex].y = qyyk;
    debugArray[debugIndex].z = qzzk;

    debugIndex              += numberOfAtoms;
    debugArray[debugIndex].x = qxyi;
    debugArray[debugIndex].y = qxzi;
    debugArray[debugIndex].z = qyzi;

    debugIndex              += numberOfAtoms;
    debugArray[debugIndex].x = qxyk;
    debugArray[debugIndex].y = qxzk;
    debugArray[debugIndex].z = qyzk;

    debugIndex              += numberOfAtoms;
    debugArray[debugIndex].x = gqxx[1];
    debugArray[debugIndex].y = gqxx[2];
    debugArray[debugIndex].z = gqxx[3];
    debugArray[debugIndex].w = gqxx[4];

    debugIndex              += numberOfAtoms;
    debugArray[debugIndex].x = gqxx[5];
    debugArray[debugIndex].y = gqxx[6];
    debugArray[debugIndex].z = gqxx[7];
    debugArray[debugIndex].w = gqxx[8];

    debugIndex              += numberOfAtoms;
    debugArray[debugIndex].x = gqxx[9];
    debugArray[debugIndex].y = gqxx[10];

    debugIndex              += numberOfAtoms;
    debugArray[debugIndex].x = gqxy[1];
    debugArray[debugIndex].y = gqxy[2];
    debugArray[debugIndex].z = gqxy[3];
    debugArray[debugIndex].w = gqxy[4];

    debugIndex              += numberOfAtoms;
    debugArray[debugIndex].x = gqxy[5];
    debugArray[debugIndex].y = gqxy[6];
    debugArray[debugIndex].z = gqxy[7];
    debugArray[debugIndex].w = gqxy[8];

    debugIndex              += numberOfAtoms;
    debugArray[debugIndex].x = gqxy[9];
    debugArray[debugIndex].y = gqxy[10];

    debugIndex              += numberOfAtoms;
    debugArray[debugIndex].x = gqxz[1];
    debugArray[debugIndex].y = gqxz[2];
    debugArray[debugIndex].z = gqxz[3];
    debugArray[debugIndex].w = gqxz[4];

    debugIndex              += numberOfAtoms;
    debugArray[debugIndex].x = gqxz[5];
    debugArray[debugIndex].y = gqxz[6];
    debugArray[debugIndex].z = gqxz[7];
    debugArray[debugIndex].w = gqxz[8];

    debugIndex              += numberOfAtoms;
    debugArray[debugIndex].x = gqxz[9];
    debugArray[debugIndex].y = gqxz[10];

    // float a[6][4];
    // float b[5][3];
    unsigned int idx1, idx2;

    idx1                     = 0;
    idx2                     = 0;
    debugIndex              += numberOfAtoms;
    debugArray[debugIndex].x = a[idx1][idx2++];
    debugArray[debugIndex].y = a[idx1][idx2++];
    debugArray[debugIndex].z = a[idx1][idx2++];
    debugArray[debugIndex].w = a[idx1][idx2++];

    idx1                     = 1;
    idx2                     = 0;
    debugIndex              += numberOfAtoms;
    debugArray[debugIndex].x = a[idx1][idx2++];
    debugArray[debugIndex].y = a[idx1][idx2++];
    debugArray[debugIndex].z = a[idx1][idx2++];
    debugArray[debugIndex].w = a[idx1][idx2++];

    idx1                     = 2;
    idx2                     = 0;
    debugIndex              += numberOfAtoms;
    debugArray[debugIndex].x = a[idx1][idx2++];
    debugArray[debugIndex].y = a[idx1][idx2++];
    debugArray[debugIndex].z = a[idx1][idx2++];
    debugArray[debugIndex].w = a[idx1][idx2++];

    idx1                     = 3;
    idx2                     = 0;
    debugIndex              += numberOfAtoms;
    debugArray[debugIndex].x = a[idx1][idx2++];
    debugArray[debugIndex].y = a[idx1][idx2++];
    debugArray[debugIndex].z = a[idx1][idx2++];
    debugArray[debugIndex].w = a[idx1][idx2++];

    idx1                     = 4;
    idx2                     = 0;
    debugIndex              += numberOfAtoms;
    debugArray[debugIndex].x = a[idx1][idx2++];
    debugArray[debugIndex].y = a[idx1][idx2++];
    debugArray[debugIndex].z = a[idx1][idx2++];
    debugArray[debugIndex].w = a[idx1][idx2++];

    idx1                     = 5;
    idx2                     = 0;
    debugIndex              += numberOfAtoms;
    debugArray[debugIndex].x = a[idx1][idx2++];
    debugArray[debugIndex].y = a[idx1][idx2++];
    debugArray[debugIndex].z = a[idx1][idx2++];
    debugArray[debugIndex].w = a[idx1][idx2++];

    idx1                     = 0;
    idx2                     = 0;
    debugIndex              += numberOfAtoms;
    debugArray[debugIndex].x = b[idx1][idx2++];
    debugArray[debugIndex].y = b[idx1][idx2++];
    debugArray[debugIndex].z = b[idx1][idx2++];

    idx1                     = 1;
    idx2                     = 0;
    debugIndex              += numberOfAtoms;
    debugArray[debugIndex].x = b[idx1][idx2++];
    debugArray[debugIndex].y = b[idx1][idx2++];
    debugArray[debugIndex].z = b[idx1][idx2++];

    idx1                     = 2;
    idx2                     = 0;
    debugIndex              += numberOfAtoms;
    debugArray[debugIndex].x = b[idx1][idx2++];
    debugArray[debugIndex].y = b[idx1][idx2++];
    debugArray[debugIndex].z = b[idx1][idx2++];

    idx1                     = 3;
    idx2                     = 0;
    debugIndex              += numberOfAtoms;
    debugArray[debugIndex].x = b[idx1][idx2++];
    debugArray[debugIndex].y = b[idx1][idx2++];
    debugArray[debugIndex].z = b[idx1][idx2++];

    idx1                     = 4;
    idx2                     = 0;
    debugIndex              += numberOfAtoms;
    debugArray[debugIndex].x = b[idx1][idx2++];
    debugArray[debugIndex].y = b[idx1][idx2++];
    debugArray[debugIndex].z = b[idx1][idx2++];

}
#endif
#endif
    ewk = ci*(uxk*gux[1]+uyk*guy[1]+uzk*guz[1])
                      -ck*(uxi*gc[2]+uyi*gc[3]+uzi*gc[4])
                 +ci*(qxxk*gqxx[1]+qyyk*gqyy[1]+qzzk*gqzz[1]
              +2.0f*(qxyk*gqxy[1]+qxzk*gqxz[1]+qyzk*gqyz[1]))
                 +ck*(qxxi*gc[5]+qyyi*gc[8]+qzzi*gc[10]
              +2.0f*(qxyi*gc[6]+qxzi*gc[7]+qyzi*gc[9]))
               - uxi*(qxxk*gqxx[2]+qyyk*gqyy[2]+qzzk*gqzz[2]
              +2.0f*(qxyk*gqxy[2]+qxzk*gqxz[2]+qyzk*gqyz[2]))
               - uyi*(qxxk*gqxx[3]+qyyk*gqyy[3]+qzzk*gqzz[3]
              +2.0f*(qxyk*gqxy[3]+qxzk*gqxz[3]+qyzk*gqyz[3]))
               - uzi*(qxxk*gqxx[4]+qyyk*gqyy[4]+qzzk*gqzz[4]
              +2.0f*(qxyk*gqxy[4]+qxzk*gqxz[4]+qyzk*gqyz[4]))
               + uxk*(qxxi*gux[5]+qyyi*gux[8]+qzzi*gux[10]
              +2.0f*(qxyi*gux[6]+qxzi*gux[7]+qyzi*gux[9]))
               + uyk*(qxxi*guy[5]+qyyi*guy[8]+qzzi*guy[10]
              +2.0f*(qxyi*guy[6]+qxzi*guy[7]+qyzi*guy[9]))
               + uzk*(qxxi*guz[5]+qyyi*guz[8]+qzzi*guz[10]
              +2.0f*(qxyi*guz[6]+qxzi*guz[7]+qyzi*guz[9]))
              + qxxi*(qxxk*gqxx[5]+qyyk*gqyy[5]+qzzk*gqzz[5]
              +2.0f*(qxyk*gqxy[5]+qxzk*gqxz[5]+qyzk*gqyz[5]))
              + qyyi*(qxxk*gqxx[8]+qyyk*gqyy[8]+qzzk*gqzz[8]
              +2.0f*(qxyk*gqxy[8]+qxzk*gqxz[8]+qyzk*gqyz[8]))
              + qzzi*(qxxk*gqxx[10]+qyyk*gqyy[10]+qzzk*gqzz[10]
              +2.0f*(qxyk*gqxy[10]+qxzk*gqxz[10]+qyzk*gqyz[10]))
       + 2.0f*(qxyi*(qxxk*gqxx[6]+qyyk*gqyy[6]+qzzk*gqzz[6]
              +2.0f*(qxyk*gqxy[6]+qxzk*gqxz[6]+qyzk*gqyz[6]))
              + qxzi*(qxxk*gqxx[7]+qyyk*gqyy[7]+qzzk*gqzz[7]
              +2.0f*(qxyk*gqxy[7]+qxzk*gqxz[7]+qyzk*gqyz[7]))
              + qyzi*(qxxk*gqxx[9]+qyyk*gqyy[9]+qzzk*gqzz[9]
              +2.0f*(qxyk*gqxy[9]+qxzk*gqxz[9]+qyzk*gqyz[9])));
 
    desymdx = ci * ck * gc[2] - (uxi*(uxk*gux[5]+uyk*guy[5]+uzk*guz[5])
                              +  uyi*(uxk*gux[6]+uyk*guy[6]+uzk*guz[6])
                              +  uzi*(uxk*gux[7]+uyk*guy[7]+uzk*guz[7]));

    dewidx = ci*(uxk*gc[5]+uyk*gc[6]+uzk*gc[7])
                      -ck*(uxi*gux[2]+uyi*guy[2]+uzi*guz[2])
                 +ci*(qxxk*gc[11]+qyyk*gc[14]+qzzk*gc[16]
              +2.0f*(qxyk*gc[12]+qxzk*gc[13]+qyzk*gc[15]))
                 +ck*(qxxi*gqxx[2]+qyyi*gqyy[2]+qzzi*gqzz[2]
              +2.0f*(qxyi*gqxy[2]+qxzi*gqxz[2]+qyzi*gqyz[2]))
               - uxi*(qxxk*gux[11]+qyyk*gux[14]+qzzk*gux[16]
              +2.0f*(qxyk*gux[12]+qxzk*gux[13]+qyzk*gux[15]))
               - uyi*(qxxk*guy[11]+qyyk*guy[14]+qzzk*guy[16]
              +2.0f*(qxyk*guy[12]+qxzk*guy[13]+qyzk*guy[15]))
               - uzi*(qxxk*guz[11]+qyyk*guz[14]+qzzk*guz[16]
              +2.0f*(qxyk*guz[12]+qxzk*guz[13]+qyzk*guz[15]))
               + uxk*(qxxi*gqxx[5]+qyyi*gqyy[5]+qzzi*gqzz[5]
              +2.0f*(qxyi*gqxy[5]+qxzi*gqxz[5]+qyzi*gqyz[5]))
               + uyk*(qxxi*gqxx[6]+qyyi*gqyy[6]+qzzi*gqzz[6]
              +2.0f*(qxyi*gqxy[6]+qxzi*gqxz[6]+qyzi*gqyz[6]))
               + uzk*(qxxi*gqxx[7]+qyyi*gqyy[7]+qzzi*gqzz[7]
              +2.0f*(qxyi*gqxy[7]+qxzi*gqxz[7]+qyzi*gqyz[7]))
              + qxxi*(qxxk*gqxx[11]+qyyk*gqxx[14]+qzzk*gqxx[16]
              +2.0f*(qxyk*gqxx[12]+qxzk*gqxx[13]+qyzk*gqxx[15]))
              + qyyi*(qxxk*gqyy[11]+qyyk*gqyy[14]+qzzk*gqyy[16]
              +2.0f*(qxyk*gqyy[12]+qxzk*gqyy[13]+qyzk*gqyy[15]))
              + qzzi*(qxxk*gqzz[11]+qyyk*gqzz[14]+qzzk*gqzz[16]
              +2.0f*(qxyk*gqzz[12]+qxzk*gqzz[13]+qyzk*gqzz[15]))
       + 2.0f*(qxyi*(qxxk*gqxy[11]+qyyk*gqxy[14]+qzzk*gqxy[16]
              +2.0f*(qxyk*gqxy[12]+qxzk*gqxy[13]+qyzk*gqxy[15]))
              + qxzi*(qxxk*gqxz[11]+qyyk*gqxz[14]+qzzk*gqxz[16]
              +2.0f*(qxyk*gqxz[12]+qxzk*gqxz[13]+qyzk*gqxz[15]))
              + qyzi*(qxxk*gqyz[11]+qyyk*gqyz[14]+qzzk*gqyz[16]
              +2.0f*(qxyk*gqyz[12]+qxzk*gqyz[13]+qyzk*gqyz[15])));

    dewkdx = ci*(uxk*gux[2]+uyk*guy[2]+uzk*guz[2])
                      -ck*(uxi*gc[5]+uyi*gc[6]+uzi*gc[7])
                 +ci*(qxxk*gqxx[2]+qyyk*gqyy[2]+qzzk*gqzz[2]
              +2.0f*(qxyk*gqxy[2]+qxzk*gqxz[2]+qyzk*gqyz[2]))
                 +ck*(qxxi*gc[11]+qyyi*gc[14]+qzzi*gc[16]
              +2.0f*(qxyi*gc[12]+qxzi*gc[13]+qyzi*gc[15]))
               - uxi*(qxxk*gqxx[5]+qyyk*gqyy[5]+qzzk*gqzz[5]
              +2.0f*(qxyk*gqxy[5]+qxzk*gqxz[5]+qyzk*gqyz[5]))
               - uyi*(qxxk*gqxx[6]+qyyk*gqyy[6]+qzzk*gqzz[6]
              +2.0f*(qxyk*gqxy[6]+qxzk*gqxz[6]+qyzk*gqyz[6]))
               - uzi*(qxxk*gqxx[7]+qyyk*gqyy[7]+qzzk*gqzz[7]
              +2.0f*(qxyk*gqxy[7]+qxzk*gqxz[7]+qyzk*gqyz[7]))
               + uxk*(qxxi*gux[11]+qyyi*gux[14]+qzzi*gux[16]
              +2.0f*(qxyi*gux[12]+qxzi*gux[13]+qyzi*gux[15]))
               + uyk*(qxxi*guy[11]+qyyi*guy[14]+qzzi*guy[16]
              +2.0f*(qxyi*guy[12]+qxzi*guy[13]+qyzi*guy[15]))
               + uzk*(qxxi*guz[11]+qyyi*guz[14]+qzzi*guz[16]
              +2.0f*(qxyi*guz[12]+qxzi*guz[13]+qyzi*guz[15]))
              + qxxi*(qxxk*gqxx[11]+qyyk*gqyy[11]+qzzk*gqzz[11]
              +2.0f*(qxyk*gqxy[11]+qxzk*gqxz[11]+qyzk*gqyz[11]))
              + qyyi*(qxxk*gqxx[14]+qyyk*gqyy[14]+qzzk*gqzz[14]
              +2.0f*(qxyk*gqxy[14]+qxzk*gqxz[14]+qyzk*gqyz[14]))
              + qzzi*(qxxk*gqxx[16]+qyyk*gqyy[16]+qzzk*gqzz[16]
              +2.0f*(qxyk*gqxy[16]+qxzk*gqxz[16]+qyzk*gqyz[16]))
       + 2.0f*(qxyi*(qxxk*gqxx[12]+qyyk*gqyy[12]+qzzk*gqzz[12]
              +2.0f*(qxyk*gqxy[12]+qxzk*gqxz[12]+qyzk*gqyz[12]))
              + qxzi*(qxxk*gqxx[13]+qyyk*gqyy[13]+qzzk*gqzz[13]
              +2.0f*(qxyk*gqxy[13]+qxzk*gqxz[13]+qyzk*gqyz[13]))
              + qyzi*(qxxk*gqxx[15]+qyyk*gqyy[15]+qzzk*gqzz[15]
              +2.0f*(qxyk*gqxy[15]+qxzk*gqxz[15]+qyzk*gqyz[15])));

    dedx = desymdx + 0.5f*(dewidx + dewkdx);
 
    desymdy = ci * ck * gc[3]
                           - (uxi*(uxk*gux[6]+uyk*guy[6]+uzk*guz[6])
                             +uyi*(uxk*gux[8]+uyk*guy[8]+uzk*guz[8])
                             +uzi*(uxk*gux[9]+uyk*guy[9]+uzk*guz[9]));

    dewidy = ci*(uxk*gc[6]+uyk*gc[8]+uzk*gc[9])
                      -ck*(uxi*gux[3]+uyi*guy[3]+uzi*guz[3])
                 +ci*(qxxk*gc[12]+qyyk*gc[17]+qzzk*gc[19]
              +2.0f*(qxyk*gc[14]+qxzk*gc[15]+qyzk*gc[18]))
                 +ck*(qxxi*gqxx[3]+qyyi*gqyy[3]+qzzi*gqzz[3]
              +2.0f*(qxyi*gqxy[3]+qxzi*gqxz[3]+qyzi*gqyz[3]))
               - uxi*(qxxk*gux[12]+qyyk*gux[17]+qzzk*gux[19]
              +2.0f*(qxyk*gux[14]+qxzk*gux[15]+qyzk*gux[18]))
               - uyi*(qxxk*guy[12]+qyyk*guy[17]+qzzk*guy[19]
              +2.0f*(qxyk*guy[14]+qxzk*guy[15]+qyzk*guy[18]))
               - uzi*(qxxk*guz[12]+qyyk*guz[17]+qzzk*guz[19]
              +2.0f*(qxyk*guz[14]+qxzk*guz[15]+qyzk*guz[18]))
               + uxk*(qxxi*gqxx[6]+qyyi*gqyy[6]+qzzi*gqzz[6]
              +2.0f*(qxyi*gqxy[6]+qxzi*gqxz[6]+qyzi*gqyz[6]))
               + uyk*(qxxi*gqxx[8]+qyyi*gqyy[8]+qzzi*gqzz[8]
              +2.0f*(qxyi*gqxy[8]+qxzi*gqxz[8]+qyzi*gqyz[8]))
               + uzk*(qxxi*gqxx[9]+qyyi*gqyy[9]+qzzi*gqzz[9]
              +2.0f*(qxyi*gqxy[9]+qxzi*gqxz[9]+qyzi*gqyz[9]))
              + qxxi*(qxxk*gqxx[12]+qyyk*gqxx[17]+qzzk*gqxx[19]
              +2.0f*(qxyk*gqxx[14]+qxzk*gqxx[15]+qyzk*gqxx[18]))
              + qyyi*(qxxk*gqyy[12]+qyyk*gqyy[17]+qzzk*gqyy[19]
              +2.0f*(qxyk*gqyy[14]+qxzk*gqyy[15]+qyzk*gqyy[18]))
              + qzzi*(qxxk*gqzz[12]+qyyk*gqzz[17]+qzzk*gqzz[19]
              +2.0f*(qxyk*gqzz[14]+qxzk*gqzz[15]+qyzk*gqzz[18]))
       + 2.0f*(qxyi*(qxxk*gqxy[12]+qyyk*gqxy[17]+qzzk*gqxy[19]
              +2.0f*(qxyk*gqxy[14]+qxzk*gqxy[15]+qyzk*gqxy[18]))
              + qxzi*(qxxk*gqxz[12]+qyyk*gqxz[17]+qzzk*gqxz[19]
              +2.0f*(qxyk*gqxz[14]+qxzk*gqxz[15]+qyzk*gqxz[18]))
              + qyzi*(qxxk*gqyz[12]+qyyk*gqyz[17]+qzzk*gqyz[19]
              +2.0f*(qxyk*gqyz[14]+qxzk*gqyz[15]+qyzk*gqyz[18])));

    dewkdy = ci*(uxk*gux[3]+uyk*guy[3]+uzk*guz[3])
                      -ck*(uxi*gc[6]+uyi*gc[8]+uzi*gc[9])
                 +ci*(qxxk*gqxx[3]+qyyk*gqyy[3]+qzzk*gqzz[3]
              +2.0f*(qxyk*gqxy[3]+qxzk*gqxz[3]+qyzk*gqyz[3]))
                 +ck*(qxxi*gc[12]+qyyi*gc[17]+qzzi*gc[19]
              +2.0f*(qxyi*gc[14]+qxzi*gc[15]+qyzi*gc[18]))
               - uxi*(qxxk*gqxx[6]+qyyk*gqyy[6]+qzzk*gqzz[6]
              +2.0f*(qxyk*gqxy[6]+qxzk*gqxz[6]+qyzk*gqyz[6]))
               - uyi*(qxxk*gqxx[8]+qyyk*gqyy[8]+qzzk*gqzz[8]
              +2.0f*(qxyk*gqxy[8]+qxzk*gqxz[8]+qyzk*gqyz[8]))
               - uzi*(qxxk*gqxx[9]+qyyk*gqyy[9]+qzzk*gqzz[9]
              +2.0f*(qxyk*gqxy[9]+qxzk*gqxz[9]+qyzk*gqyz[9]))
               + uxk*(qxxi*gux[12]+qyyi*gux[17]+qzzi*gux[19]
              +2.0f*(qxyi*gux[14]+qxzi*gux[15]+qyzi*gux[18]))
               + uyk*(qxxi*guy[12]+qyyi*guy[17]+qzzi*guy[19]
              +2.0f*(qxyi*guy[14]+qxzi*guy[15]+qyzi*guy[18]))
               + uzk*(qxxi*guz[12]+qyyi*guz[17]+qzzi*guz[19]
              +2.0f*(qxyi*guz[14]+qxzi*guz[15]+qyzi*guz[18]))
              + qxxi*(qxxk*gqxx[12]+qyyk*gqyy[12]+qzzk*gqzz[12]
              +2.0f*(qxyk*gqxy[12]+qxzk*gqxz[12]+qyzk*gqyz[12]))
              + qyyi*(qxxk*gqxx[17]+qyyk*gqyy[17]+qzzk*gqzz[17]
              +2.0f*(qxyk*gqxy[17]+qxzk*gqxz[17]+qyzk*gqyz[17]))
              + qzzi*(qxxk*gqxx[19]+qyyk*gqyy[19]+qzzk*gqzz[19]
              +2.0f*(qxyk*gqxy[19]+qxzk*gqxz[19]+qyzk*gqyz[19]))
       + 2.0f*(qxyi*(qxxk*gqxx[14]+qyyk*gqyy[14]+qzzk*gqzz[14]
              +2.0f*(qxyk*gqxy[14]+qxzk*gqxz[14]+qyzk*gqyz[14]))
              + qxzi*(qxxk*gqxx[15]+qyyk*gqyy[15]+qzzk*gqzz[15]
              +2.0f*(qxyk*gqxy[15]+qxzk*gqxz[15]+qyzk*gqyz[15]))
              + qyzi*(qxxk*gqxx[18]+qyyk*gqyy[18]+qzzk*gqzz[18]
              +2.0f*(qxyk*gqxy[18]+qxzk*gqxz[18]+qyzk*gqyz[18])));

    dedy = desymdy + 0.5f*(dewidy + dewkdy);
 
    desymdz = ci * ck * gc[4]
                           - (uxi*(uxk*gux[7]+uyk*guy[7]+uzk*guz[7])
                             +uyi*(uxk*gux[9]+uyk*guy[9]+uzk*guz[9])
                             +uzi*(uxk*gux[10]+uyk*guy[10]+uzk*guz[10]));

    dewidz = ci*(uxk*gc[7]+uyk*gc[9]+uzk*gc[10])
                      -ck*(uxi*gux[4]+uyi*guy[4]+uzi*guz[4])
                 +ci*(qxxk*gc[13]+qyyk*gc[18]+qzzk*gc[20]
              +2.0f*(qxyk*gc[15]+qxzk*gc[16]+qyzk*gc[19]))
                 +ck*(qxxi*gqxx[4]+qyyi*gqyy[4]+qzzi*gqzz[4]
              +2.0f*(qxyi*gqxy[4]+qxzi*gqxz[4]+qyzi*gqyz[4]))
               - uxi*(qxxk*gux[13]+qyyk*gux[18]+qzzk*gux[20]
              +2.0f*(qxyk*gux[15]+qxzk*gux[16]+qyzk*gux[19]))
               - uyi*(qxxk*guy[13]+qyyk*guy[18]+qzzk*guy[20]
              +2.0f*(qxyk*guy[15]+qxzk*guy[16]+qyzk*guy[19]))
               - uzi*(qxxk*guz[13]+qyyk*guz[18]+qzzk*guz[20]
              +2.0f*(qxyk*guz[15]+qxzk*guz[16]+qyzk*guz[19]))
               + uxk*(qxxi*gqxx[7]+qyyi*gqyy[7]+qzzi*gqzz[7]
              +2.0f*(qxyi*gqxy[7]+qxzi*gqxz[7]+qyzi*gqyz[7]))
               + uyk*(qxxi*gqxx[9]+qyyi*gqyy[9]+qzzi*gqzz[9]
              +2.0f*(qxyi*gqxy[9]+qxzi*gqxz[9]+qyzi*gqyz[9]))
               + uzk*(qxxi*gqxx[10]+qyyi*gqyy[10]+qzzi*gqzz[10]
              +2.0f*(qxyi*gqxy[10]+qxzi*gqxz[10]+qyzi*gqyz[10]))
              + qxxi*(qxxk*gqxx[13]+qyyk*gqxx[18]+qzzk*gqxx[20]
              +2.0f*(qxyk*gqxx[15]+qxzk*gqxx[16]+qyzk*gqxx[19]))
              + qyyi*(qxxk*gqyy[13]+qyyk*gqyy[18]+qzzk*gqyy[20]
              +2.0f*(qxyk*gqyy[15]+qxzk*gqyy[16]+qyzk*gqyy[19]))
              + qzzi*(qxxk*gqzz[13]+qyyk*gqzz[18]+qzzk*gqzz[20]
              +2.0f*(qxyk*gqzz[15]+qxzk*gqzz[16]+qyzk*gqzz[19]))
       + 2.0f*(qxyi*(qxxk*gqxy[13]+qyyk*gqxy[18]+qzzk*gqxy[20]
              +2.0f*(qxyk*gqxy[15]+qxzk*gqxy[16]+qyzk*gqxy[19]))
              + qxzi*(qxxk*gqxz[13]+qyyk*gqxz[18]+qzzk*gqxz[20]
              +2.0f*(qxyk*gqxz[15]+qxzk*gqxz[16]+qyzk*gqxz[19]))
              + qyzi*(qxxk*gqyz[13]+qyyk*gqyz[18]+qzzk*gqyz[20]
              +2.0f*(qxyk*gqyz[15]+qxzk*gqyz[16]+qyzk*gqyz[19])));

    dewkdz = ci*(uxk*gux[4]+uyk*guy[4]+uzk*guz[4])
                      -ck*(uxi*gc[7]+uyi*gc[9]+uzi*gc[10])
                 +ci*(qxxk*gqxx[4]+qyyk*gqyy[4]+qzzk*gqzz[4]
              +2.0f*(qxyk*gqxy[4]+qxzk*gqxz[4]+qyzk*gqyz[4]))
                 +ck*(qxxi*gc[13]+qyyi*gc[18]+qzzi*gc[20]
              +2.0f*(qxyi*gc[15]+qxzi*gc[16]+qyzi*gc[19]))
               - uxi*(qxxk*gqxx[7]+qyyk*gqyy[7]+qzzk*gqzz[7]
              +2.0f*(qxyk*gqxy[7]+qxzk*gqxz[7]+qyzk*gqyz[7]))
               - uyi*(qxxk*gqxx[9]+qyyk*gqyy[9]+qzzk*gqzz[9]
              +2.0f*(qxyk*gqxy[9]+qxzk*gqxz[9]+qyzk*gqyz[9]))
               - uzi*(qxxk*gqxx[10]+qyyk*gqyy[10]+qzzk*gqzz[10]
              +2.0f*(qxyk*gqxy[10]+qxzk*gqxz[10]+qyzk*gqyz[10]))
               + uxk*(qxxi*gux[13]+qyyi*gux[18]+qzzi*gux[20]
              +2.0f*(qxyi*gux[15]+qxzi*gux[16]+qyzi*gux[19]))
               + uyk*(qxxi*guy[13]+qyyi*guy[18]+qzzi*guy[20]
              +2.0f*(qxyi*guy[15]+qxzi*guy[16]+qyzi*guy[19]))
               + uzk*(qxxi*guz[13]+qyyi*guz[18]+qzzi*guz[20]
              +2.0f*(qxyi*guz[15]+qxzi*guz[16]+qyzi*guz[19]))
              + qxxi*(qxxk*gqxx[13]+qyyk*gqyy[13]+qzzk*gqzz[13]
              +2.0f*(qxyk*gqxy[13]+qxzk*gqxz[13]+qyzk*gqyz[13]))
              + qyyi*(qxxk*gqxx[18]+qyyk*gqyy[18]+qzzk*gqzz[18]
              +2.0f*(qxyk*gqxy[18]+qxzk*gqxz[18]+qyzk*gqyz[18]))
              + qzzi*(qxxk*gqxx[20]+qyyk*gqyy[20]+qzzk*gqzz[20]
              +2.0f*(qxyk*gqxy[20]+qxzk*gqxz[20]+qyzk*gqyz[20]))
       + 2.0f*(qxyi*(qxxk*gqxx[15]+qyyk*gqyy[15]+qzzk*gqzz[15]
              +2.0f*(qxyk*gqxy[15]+qxzk*gqxz[15]+qyzk*gqyz[15]))
              + qxzi*(qxxk*gqxx[16]+qyyk*gqyy[16]+qzzk*gqzz[16]
              +2.0f*(qxyk*gqxy[16]+qxzk*gqxz[16]+qyzk*gqyz[16]))
              + qyzi*(qxxk*gqxx[19]+qyyk*gqyy[19]+qzzk*gqzz[19]
              +2.0f*(qxyk*gqxy[19]+qxzk*gqxz[19]+qyzk*gqyz[19])));

    dedz = desymdz + 0.5f*(dewidz + dewkdz);
 
    desymdr = ci * ck * gc[21]
                           - (uxi*(uxk*gux[22]+uyk*guy[22]+uzk*guz[22])
                             +uyi*(uxk*gux[23]+uyk*guy[23]+uzk*guz[23])
                             +uzi*(uxk*gux[24]+uyk*guy[24]+uzk*guz[24]));

    dewidr = ci*(uxk*gc[22]+uyk*gc[23]+uzk*gc[24])
                      -ck*(uxi*gux[21]+uyi*guy[21]+uzi*guz[21])
                 +ci*(qxxk*gc[25]+qyyk*gc[28]+qzzk*gc[30]
              +2.0f*(qxyk*gc[26]+qxzk*gc[27]+qyzk*gc[29]))
                 +ck*(qxxi*gqxx[21]+qyyi*gqyy[21]+qzzi*gqzz[21]
              +2.0f*(qxyi*gqxy[21]+qxzi*gqxz[21]+qyzi*gqyz[21]))
               - uxi*(qxxk*gux[25]+qyyk*gux[28]+qzzk*gux[30]
              +2.0f*(qxyk*gux[26]+qxzk*gux[27]+qyzk*gux[29]))
               - uyi*(qxxk*guy[25]+qyyk*guy[28]+qzzk*guy[30]
              +2.0f*(qxyk*guy[26]+qxzk*guy[27]+qyzk*guy[29]))
               - uzi*(qxxk*guz[25]+qyyk*guz[28]+qzzk*guz[30]
              +2.0f*(qxyk*guz[26]+qxzk*guz[27]+qyzk*guz[29]))
               + uxk*(qxxi*gqxx[22]+qyyi*gqyy[22]+qzzi*gqzz[22]
              +2.0f*(qxyi*gqxy[22]+qxzi*gqxz[22]+qyzi*gqyz[22]))
               + uyk*(qxxi*gqxx[23]+qyyi*gqyy[23]+qzzi*gqzz[23]
              +2.0f*(qxyi*gqxy[23]+qxzi*gqxz[23]+qyzi*gqyz[23]))
               + uzk*(qxxi*gqxx[24]+qyyi*gqyy[24]+qzzi*gqzz[24]
              +2.0f*(qxyi*gqxy[24]+qxzi*gqxz[24]+qyzi*gqyz[24]))
              + qxxi*(qxxk*gqxx[25]+qyyk*gqxx[28]+qzzk*gqxx[30]
              +2.0f*(qxyk*gqxx[26]+qxzk*gqxx[27]+qyzk*gqxx[29]))
              + qyyi*(qxxk*gqyy[25]+qyyk*gqyy[28]+qzzk*gqyy[30]
              +2.0f*(qxyk*gqyy[26]+qxzk*gqyy[27]+qyzk*gqyy[29]))
              + qzzi*(qxxk*gqzz[25]+qyyk*gqzz[28]+qzzk*gqzz[30]
              +2.0f*(qxyk*gqzz[26]+qxzk*gqzz[27]+qyzk*gqzz[29]))
              + 2.0f*(qxyi*(qxxk*gqxy[25]+qyyk*gqxy[28]+qzzk*gqxy[30]
              +2.0f*(qxyk*gqxy[26]+qxzk*gqxy[27]+qyzk*gqxy[29]))
              + qxzi*(qxxk*gqxz[25]+qyyk*gqxz[28]+qzzk*gqxz[30]
              +2.0f*(qxyk*gqxz[26]+qxzk*gqxz[27]+qyzk*gqxz[29]))
              + qyzi*(qxxk*gqyz[25]+qyyk*gqyz[28]+qzzk*gqyz[30]
              +2.0f*(qxyk*gqyz[26]+qxzk*gqyz[27]+qyzk*gqyz[29])));

    dewkdr = ci*(uxk*gux[21]+uyk*guy[21]+uzk*guz[21])
                      -ck*(uxi*gc[22]+uyi*gc[23]+uzi*gc[24])
                 +ci*(qxxk*gqxx[21]+qyyk*gqyy[21]+qzzk*gqzz[21]
              +2.0f*(qxyk*gqxy[21]+qxzk*gqxz[21]+qyzk*gqyz[21]))
                 +ck*(qxxi*gc[25]+qyyi*gc[28]+qzzi*gc[30]
              +2.0f*(qxyi*gc[26]+qxzi*gc[27]+qyzi*gc[29]))
               - uxi*(qxxk*gqxx[22]+qyyk*gqyy[22]+qzzk*gqzz[22]
              +2.0f*(qxyk*gqxy[22]+qxzk*gqxz[22]+qyzk*gqyz[22]))
               - uyi*(qxxk*gqxx[23]+qyyk*gqyy[23]+qzzk*gqzz[23]
              +2.0f*(qxyk*gqxy[23]+qxzk*gqxz[23]+qyzk*gqyz[23]))
               - uzi*(qxxk*gqxx[24]+qyyk*gqyy[24]+qzzk*gqzz[24]
              +2.0f*(qxyk*gqxy[24]+qxzk*gqxz[24]+qyzk*gqyz[24]))
               + uxk*(qxxi*gux[25]+qyyi*gux[28]+qzzi*gux[30]
              +2.0f*(qxyi*gux[26]+qxzi*gux[27]+qyzi*gux[29]))
               + uyk*(qxxi*guy[25]+qyyi*guy[28]+qzzi*guy[30]
              +2.0f*(qxyi*guy[26]+qxzi*guy[27]+qyzi*guy[29]))
               + uzk*(qxxi*guz[25]+qyyi*guz[28]+qzzi*guz[30]
              +2.0f*(qxyi*guz[26]+qxzi*guz[27]+qyzi*guz[29]))
              + qxxi*(qxxk*gqxx[25]+qyyk*gqyy[25]+qzzk*gqzz[25]
              +2.0f*(qxyk*gqxy[25]+qxzk*gqxz[25]+qyzk*gqyz[25]))
              + qyyi*(qxxk*gqxx[28]+qyyk*gqyy[28]+qzzk*gqzz[28]
              +2.0f*(qxyk*gqxy[28]+qxzk*gqxz[28]+qyzk*gqyz[28]))
              + qzzi*(qxxk*gqxx[30]+qyyk*gqyy[30]+qzzk*gqzz[30]
              +2.0f*(qxyk*gqxy[30]+qxzk*gqxz[30]+qyzk*gqyz[30]))
              + 2.0f*(qxyi*(qxxk*gqxx[26]+qyyk*gqyy[26]+qzzk*gqzz[26]
              +2.0f*(qxyk*gqxy[26]+qxzk*gqxz[26]+qyzk*gqyz[26]))
              + qxzi*(qxxk*gqxx[27]+qyyk*gqyy[27]+qzzk*gqzz[27]
              +2.0f*(qxyk*gqxy[27]+qxzk*gqxz[27]+qyzk*gqyz[27]))
              + qyzi*(qxxk*gqxx[29]+qyyk*gqyy[29]+qzzk*gqzz[29]
              +2.0f*(qxyk*gqxy[29]+qxzk*gqxz[29]+qyzk*gqyz[29])));

    dsumdr = desymdr + 0.5f*(dewidr + dewkdr);
    drbi = rbk*dsumdr;
    drbk = rbi*dsumdr;
 
    // torque on permanent dipoles due to permanent reaction field
 
    trq[1]   = 0.0f;
    trq[2]   = 0.0f;
    trq[3]   = 0.0f;

    trq_k[1] = 0.0f;
    trq_k[2] = 0.0f;
    trq_k[3] = 0.0f;

    if ( sameAtom == 0 )
    {

        fid[1] = uxk*gux[2] + uyk*gux[3] + uzk*gux[4]
                + 0.5f*(ck*gux[1]+qxxk*gux[5]+qyyk*gux[8]+qzzk*gux[10]
                      +2.0f*(qxyk*gux[6]+qxzk*gux[7]+qyzk*gux[9])
                      +ck*gc[2]+qxxk*gqxx[2]+qyyk*gqyy[2]+qzzk*gqzz[2]
                      +2.0f*(qxyk*gqxy[2]+qxzk*gqxz[2]+qyzk*gqyz[2]));

        fid[2] = uxk*guy[2] + uyk*guy[3] + uzk*guy[4]
                + 0.5f*(ck*guy[1]+qxxk*guy[5]+qyyk*guy[8]+qzzk*guy[10]
                      +2.0f*(qxyk*guy[6]+qxzk*guy[7]+qyzk*guy[9])
                      +ck*gc[3]+qxxk*gqxx[3]+qyyk*gqyy[3]+qzzk*gqzz[3]
                      +2.0f*(qxyk*gqxy[3]+qxzk*gqxz[3]+qyzk*gqyz[3]));

        fid[3] = uxk*guz[2] + uyk*guz[3] + uzk*guz[4]
                + 0.5f*(ck*guz[1]+qxxk*guz[5]+qyyk*guz[8]+qzzk*guz[10]
                      +2.0f*(qxyk*guz[6]+qxzk*guz[7]+qyzk*guz[9])
                      +ck*gc[4]+qxxk*gqxx[4]+qyyk*gqyy[4]+qzzk*gqzz[4]
                      +2.0f*(qxyk*gqxy[4]+qxzk*gqxz[4]+qyzk*gqyz[4]));

        fkd[1] = uxi*gux[2] + uyi*gux[3] + uzi*gux[4]
                - 0.5f*(ci*gux[1]+qxxi*gux[5]+qyyi*gux[8]+qzzi*gux[10]
                      +2.0f*(qxyi*gux[6]+qxzi*gux[7]+qyzi*gux[9])
                      +ci*gc[2]+qxxi*gqxx[2]+qyyi*gqyy[2]+qzzi*gqzz[2]
                      +2.0f*(qxyi*gqxy[2]+qxzi*gqxz[2]+qyzi*gqyz[2]));

        fkd[2] = uxi*guy[2] + uyi*guy[3] + uzi*guy[4]
                - 0.5f*(ci*guy[1]+qxxi*guy[5]+qyyi*guy[8]+qzzi*guy[10]
                      +2.0f*(qxyi*guy[6]+qxzi*guy[7]+qyzi*guy[9])
                      +ci*gc[3]+qxxi*gqxx[3]+qyyi*gqyy[3]+qzzi*gqzz[3]
                      +2.0f*(qxyi*gqxy[3]+qxzi*gqxz[3]+qyzi*gqyz[3]));

        fkd[3] = uxi*guz[2] + uyi*guz[3] + uzi*guz[4]
                - 0.5f*(ci*guz[1]+qxxi*guz[5]+qyyi*guz[8]+qzzi*guz[10]
                      +2.0f*(qxyi*guz[6]+qxzi*guz[7]+qyzi*guz[9])
                      +ci*gc[4]+qxxi*gqxx[4]+qyyi*gqyy[4]+qzzi*gqzz[4]
                      +2.0f*(qxyi*gqxy[4]+qxzi*gqxz[4]+qyzi*gqyz[4]));

        trq[1]    = uyi*fid[3] - uzi*fid[2];
        trq[2]    = uzi*fid[1] - uxi*fid[3];
        trq[3]    = uxi*fid[2] - uyi*fid[1];

        trq_k[1]  = uyk*fkd[3] - uzk*fkd[2];
        trq_k[2]  = uzk*fkd[1] - uxk*fkd[3];
        trq_k[3]  = uxk*fkd[2] - uyk*fkd[1];
 
        // torque on quadrupoles due to permanent reaction field gradient
 
        fidg[1][1] =
                - 0.5f*(ck*gqxx[1]+uxk*gqxx[2]+uyk*gqxx[3]+uzk*gqxx[4]
                      +qxxk*gqxx[5]+qyyk*gqxx[8]+qzzk*gqxx[10]
                      +2.0f*(qxyk*gqxx[6]+qxzk*gqxx[7]+qyzk*gqxx[9])
                      +ck*gc[5]+uxk*gux[5]+uyk*guy[5]+uzk*guz[5]
                      +qxxk*gqxx[5]+qyyk*gqyy[5]+qzzk*gqzz[5]
                      +2.0f*(qxyk*gqxy[5]+qxzk*gqxz[5]+qyzk*gqyz[5]));

        fidg[1][2] =
                - 0.5f*(ck*gqxy[1]+uxk*gqxy[2]+uyk*gqxy[3]+uzk*gqxy[4]
                      +qxxk*gqxy[5]+qyyk*gqxy[8]+qzzk*gqxy[10]
                      +2.0f*(qxyk*gqxy[6]+qxzk*gqxy[7]+qyzk*gqxy[9])
                      +ck*gc[6]+uxk*gux[6]+uyk*guy[6]+uzk*guz[6]
                      +qxxk*gqxx[6]+qyyk*gqyy[6]+qzzk*gqzz[6]
                      +2.0f*(qxyk*gqxy[6]+qxzk*gqxz[6]+qyzk*gqyz[6]));

        fidg[1][3] =
                - 0.5f*(ck*gqxz[1]+uxk*gqxz[2]+uyk*gqxz[3]+uzk*gqxz[4]
                      +qxxk*gqxz[5]+qyyk*gqxz[8]+qzzk*gqxz[10]
                      +2.0f*(qxyk*gqxz[6]+qxzk*gqxz[7]+qyzk*gqxz[9])
                      +ck*gc[7]+uxk*gux[7]+uyk*guy[7]+uzk*guz[7]
                      +qxxk*gqxx[7]+qyyk*gqyy[7]+qzzk*gqzz[7]
                      +2.0f*(qxyk*gqxy[7]+qxzk*gqxz[7]+qyzk*gqyz[7]));

        fidg[2][2] =
                - 0.5f*(ck*gqyy[1]+uxk*gqyy[2]+uyk*gqyy[3]+uzk*gqyy[4]
                      +qxxk*gqyy[5]+qyyk*gqyy[8]+qzzk*gqyy[10]
                      +2.0f*(qxyk*gqyy[6]+qxzk*gqyy[7]+qyzk*gqyy[9])
                      +ck*gc[8]+uxk*gux[8]+uyk*guy[8]+uzk*guz[8]
                      +qxxk*gqxx[8]+qyyk*gqyy[8]+qzzk*gqzz[8]
                      +2.0f*(qxyk*gqxy[8]+qxzk*gqxz[8]+qyzk*gqyz[8]));

        fidg[2][3] =
                - 0.5f*(ck*gqyz[1]+uxk*gqyz[2]+uyk*gqyz[3]+uzk*gqyz[4]
                      +qxxk*gqyz[5]+qyyk*gqyz[8]+qzzk*gqyz[10]
                      +2.0f*(qxyk*gqyz[6]+qxzk*gqyz[7]+qyzk*gqyz[9])
                      +ck*gc[9]+uxk*gux[9]+uyk*guy[9]+uzk*guz[9]
                      +qxxk*gqxx[9]+qyyk*gqyy[9]+qzzk*gqzz[9]
                      +2.0f*(qxyk*gqxy[9]+qxzk*gqxz[9]+qyzk*gqyz[9]));

        fidg[3][3] =
                - 0.5f*(ck*gqzz[1]+uxk*gqzz[2]+uyk*gqzz[3]+uzk*gqzz[4]
                      +qxxk*gqzz[5]+qyyk*gqzz[8]+qzzk*gqzz[10]
                      +2.0f*(qxyk*gqzz[6]+qxzk*gqzz[7]+qyzk*gqzz[9])
                      +ck*gc[10]+uxk*gux[10]+uyk*guy[10]+uzk*guz[10]
                      +qxxk*gqxx[10]+qyyk*gqyy[10]+qzzk*gqzz[10]
                   +2.0f*(qxyk*gqxy[10]+qxzk*gqxz[10]+qyzk*gqyz[10]));

        fidg[2][1] = fidg[1][2];
        fidg[3][1] = fidg[1][3];
        fidg[3][2] = fidg[2][3];

        fkdg[1][1] =
                - 0.5f*(ci*gqxx[1]-uxi*gqxx[2]-uyi*gqxx[3]-uzi *gqxx[4]
                      +qxxi*gqxx[5]+qyyi*gqxx[8]+qzzi*gqxx[10]
                      +2.0f*(qxyi*gqxx[6]+qxzi*gqxx[7]+qyzi*gqxx[9])
                      +ci*gc[5]-uxi*gux[5]-uyi*guy[5]-uzi*guz[5]
                      +qxxi*gqxx[5]+qyyi*gqyy[5]+qzzi*gqzz[5]
                      +2.0f*(qxyi*gqxy[5]+qxzi*gqxz[5]+qyzi*gqyz[5]));

        fkdg[1][2] =
                - 0.5f*(ci*gqxy[1]-uxi*gqxy[2]-uyi*gqxy[3]-uzi*gqxy[4]
                      +qxxi*gqxy[5]+qyyi*gqxy[8]+qzzi*gqxy[10]
                      +2.0f*(qxyi*gqxy[6]+qxzi*gqxy[7]+qyzi*gqxy[9])
                      +ci*gc[6]-uxi*gux[6]-uyi*guy[6]-uzi*guz[6]
                      +qxxi*gqxx[6]+qyyi*gqyy[6]+qzzi*gqzz[6]
                      +2.0f*(qxyi*gqxy[6]+qxzi*gqxz[6]+qyzi*gqyz[6]));

        fkdg[1][3] =
                - 0.5f*(ci*gqxz[1]-uxi*gqxz[2]-uyi*gqxz[3]-uzi*gqxz[4]
                      +qxxi*gqxz[5]+qyyi*gqxz[8]+qzzi*gqxz[10]
                      +2.0f*(qxyi*gqxz[6]+qxzi*gqxz[7]+qyzi*gqxz[9])
                      +ci*gc[7]-uxi*gux[7]-uyi*guy[7]-uzi*guz[7]
                      +qxxi*gqxx[7]+qyyi*gqyy[7]+qzzi*gqzz[7]
                      +2.0f*(qxyi*gqxy[7]+qxzi*gqxz[7]+qyzi*gqyz[7]));

        fkdg[2][2] =
                - 0.5f*(ci*gqyy[1]-uxi*gqyy[2]-uyi*gqyy[3]-uzi*gqyy[4]
                      +qxxi*gqyy[5]+qyyi*gqyy[8]+qzzi*gqyy[10]
                      +2.0f*(qxyi*gqyy[6]+qxzi*gqyy[7]+qyzi*gqyy[9])
                      +ci*gc[8]-uxi*gux[8]-uyi*guy[8]-uzi*guz[8]
                      +qxxi*gqxx[8]+qyyi*gqyy[8]+qzzi*gqzz[8]
                      +2.0f*(qxyi*gqxy[8]+qxzi*gqxz[8]+qyzi*gqyz[8]));

        fkdg[2][3] =
                - 0.5f*(ci*gqyz[1]-uxi*gqyz[2]-uyi*gqyz[3]-uzi*gqyz[4]
                      +qxxi*gqyz[5]+qyyi*gqyz[8]+qzzi*gqyz[10]
                      +2.0f*(qxyi*gqyz[6]+qxzi*gqyz[7]+qyzi*gqyz[9])
                      +ci*gc[9]-uxi*gux[9]-uyi*guy[9]-uzi*guz[9]
                      +qxxi*gqxx[9]+qyyi*gqyy[9]+qzzi*gqzz[9]
                      +2.0f*(qxyi*gqxy[9]+qxzi*gqxz[9]+qyzi*gqyz[9]));
        fkdg[3][3] =
                - 0.5f*(ci*gqzz[1]-uxi*gqzz[2]-uyi*gqzz[3]-uzi*gqzz[4]
                      +qxxi*gqzz[5]+qyyi*gqzz[8]+qzzi*gqzz[10]
                      +2.0f*(qxyi*gqzz[6]+qxzi*gqzz[7]+qyzi*gqzz[9])
                      +ci*gc[10]-uxi*gux[10]-uyi*guy[10]-uzi*guz[10]
                      +qxxi*gqxx[10]+qyyi*gqyy[10]+qzzi*gqzz[10]
                    +2.0f*(qxyi*gqxy[10]+qxzi*gqxz[10]+qyzi*gqyz[10]));

        fkdg[2][1] = fkdg[1][2];
        fkdg[3][1] = fkdg[1][3];
        fkdg[3][2] = fkdg[2][3];

        trq[1]   += 2.0f* (qxyi*fidg[1][3]+qyyi*fidg[2][3]+qyzi*fidg[3][3]
                           -qxzi*fidg[1][2]-qyzi*fidg[2][2]-qzzi*fidg[3][2]);

        trq[2]   += 2.0f*(qxzi*fidg[1][1]+qyzi*fidg[2][1]+qzzi*fidg[3][1]
                         -qxxi*fidg[1][3]-qxyi*fidg[2][3]-qxzi*fidg[3][3]);

        trq[3]   += 2.0f*(qxxi*fidg[1][2]+qxyi*fidg[2][2]+qxzi*fidg[3][2]
                         -qxyi*fidg[1][1]-qyyi*fidg[2][1]-qyzi*fidg[3][1]);

        trq_k[1] += 2.0f*
                          (qxyk*fkdg[1][3]+qyyk*fkdg[2][3]+qyzk*fkdg[3][3]
                          -qxzk*fkdg[1][2]-qyzk*fkdg[2][2]-qzzk*fkdg[3][2]);

        trq_k[2] += 2.0f*
                          (qxzk*fkdg[1][1]+qyzk*fkdg[2][1]+qzzk*fkdg[3][1]
                          -qxxk*fkdg[1][3]-qxyk*fkdg[2][3]-qxzk*fkdg[3][3]);

        trq_k[3] += 2.0f*
                          (qxxk*fkdg[1][2]+qxyk*fkdg[2][2]+qxzk*fkdg[3][2]
                          -qxyk*fkdg[1][1]-qyyk*fkdg[2][1]-qyzk*fkdg[3][1]);
    }
 
    // electrostatic solvation energy of the permanent multipoles in
    // the GK reaction potential of the induced dipoles
 
    esymi =              -uxi*(dxk*gux[2]+dyk*guy[2]+dzk*guz[2])
                        - uyi*(dxk*gux[3]+dyk*guy[3]+dzk*guz[3])
                        - uzi*(dxk*gux[4]+dyk*guy[4]+dzk*guz[4])
                        - uxk*(dxi*gux[2]+dyi*guy[2]+dzi*guz[2])
                        - uyk*(dxi*gux[3]+dyi*guy[3]+dzi*guz[3])
                        - uzk*(dxi*gux[4]+dyi*guy[4]+dzi*guz[4]);

    ewii = ci*(dxk*gc[2]+dyk*gc[3]+dzk*gc[4])
                      - ck*(dxi*gux[1]+dyi*guy[1]+dzi*guz[1])
                      - dxi*(qxxk*gux[5]+qyyk*gux[8]+qzzk*gux[10]
                     +2.0f*(qxyk*gux[6]+qxzk*gux[7]+qyzk*gux[9]))
                      - dyi*(qxxk*guy[5]+qyyk*guy[8]+qzzk*guy[10]
                     +2.0f*(qxyk*guy[6]+qxzk*guy[7]+qyzk*guy[9]))
                      - dzi*(qxxk*guz[5]+qyyk*guz[8]+qzzk*guz[10]
                     +2.0f*(qxyk*guz[6]+qxzk*guz[7]+qyzk*guz[9]))
                      + dxk*(qxxi*gqxx[2]+qyyi*gqyy[2]+qzzi*gqzz[2]
                     +2.0f*(qxyi*gqxy[2]+qxzi*gqxz[2]+qyzi*gqyz[2]))
                      + dyk*(qxxi*gqxx[3]+qyyi*gqyy[3]+qzzi*gqzz[3]
                     +2.0f*(qxyi*gqxy[3]+qxzi*gqxz[3]+qyzi*gqyz[3]))
                      + dzk*(qxxi*gqxx[4]+qyyi*gqyy[4]+qzzi*gqzz[4]
                     +2.0f*(qxyi*gqxy[4]+qxzi*gqxz[4]+qyzi*gqyz[4]));

    ewki = ci*(dxk*gux[1]+dyk*guy[1]+dzk*guz[1])
                      - ck*(dxi*gc[2]+dyi*gc[3]+dzi*gc[4])
                      - dxi*(qxxk*gqxx[2]+qyyk*gqyy[2]+qzzk*gqzz[2]
                     +2.0f*(qxyk*gqxy[2]+qxzk*gqxz[2]+qyzk*gqyz[2]))
                      - dyi*(qxxk*gqxx[3]+qyyk*gqyy[3]+qzzk*gqzz[3]
                     +2.0f*(qxyk*gqxy[3]+qxzk*gqxz[3]+qyzk*gqyz[3]))
                      - dzi*(qxxk*gqxx[4]+qyyk*gqyy[4]+qzzk*gqzz[4]
                     +2.0f*(qxyk*gqxy[4]+qxzk*gqxz[4]+qyzk*gqyz[4]))
                      + dxk*(qxxi*gux[5]+qyyi*gux[8]+qzzi*gux[10]
                     +2.0f*(qxyi*gux[6]+qxzi*gux[7]+qyzi*gux[9]))
                      + dyk*(qxxi*guy[5]+qyyi*guy[8]+qzzi*guy[10]
                     +2.0f*(qxyi*guy[6]+qxzi*guy[7]+qyzi*guy[9]))
                      + dzk*(qxxi*guz[5]+qyyi*guz[8]+qzzi*guz[10]
                     +2.0f*(qxyi*guz[6]+qxzi*guz[7]+qyzi*guz[9]));
 
    // electrostatic solvation free energy gradient of the permanent
    // multipoles in the reaction potential of the induced dipoles
 
    dpsymdx = -uxi*(sxk*gux[5]+syk*guy[5]+szk*guz[5])
                          - uyi*(sxk*gux[6]+syk*guy[6]+szk*guz[6])
                          - uzi*(sxk*gux[7]+syk*guy[7]+szk*guz[7])
                          - uxk*(sxi*gux[5]+syi*guy[5]+szi*guz[5])
                          - uyk*(sxi*gux[6]+syi*guy[6]+szi*guz[6])
                          - uzk*(sxi*gux[7]+syi*guy[7]+szi*guz[7]);

    dpwidx = ci*(sxk*gc[5]+syk*gc[6]+szk*gc[7])
                        - ck*(sxi*gux[2]+syi*guy[2]+szi*guz[2])
                      - sxi*(qxxk*gux[11]+qyyk*gux[14]+qzzk*gux[16]
                     +2.0f*(qxyk*gux[12]+qxzk*gux[13]+qyzk*gux[15]))
                      - syi*(qxxk*guy[11]+qyyk*guy[14]+qzzk*guy[16]
                     +2.0f*(qxyk*guy[12]+qxzk*guy[13]+qyzk*guy[15]))
                      - szi*(qxxk*guz[11]+qyyk*guz[14]+qzzk*guz[16]
                     +2.0f*(qxyk*guz[12]+qxzk*guz[13]+qyzk*guz[15]))
                      + sxk*(qxxi*gqxx[5]+qyyi*gqyy[5]+qzzi*gqzz[5]
                     +2.0f*(qxyi*gqxy[5]+qxzi*gqxz[5]+qyzi*gqyz[5]))
                      + syk*(qxxi*gqxx[6]+qyyi*gqyy[6]+qzzi*gqzz[6]
                     +2.0f*(qxyi*gqxy[6]+qxzi*gqxz[6]+qyzi*gqyz[6]))
                      + szk*(qxxi*gqxx[7]+qyyi*gqyy[7]+qzzi*gqzz[7]
                     +2.0f*(qxyi*gqxy[7]+qxzi*gqxz[7]+qyzi*gqyz[7]));

    dpwkdx = ci*(sxk*gux[2]+syk*guy[2]+szk*guz[2])
                        - ck*(sxi*gc[5]+syi*gc[6]+szi*gc[7])
                      - sxi*(qxxk*gqxx[5]+qyyk*gqyy[5]+qzzk*gqzz[5]
                     +2.0f*(qxyk*gqxy[5]+qxzk*gqxz[5]+qyzk*gqyz[5]))
                      - syi*(qxxk*gqxx[6]+qyyk*gqyy[6]+qzzk*gqzz[6]
                     +2.0f*(qxyk*gqxy[6]+qxzk*gqxz[6]+qyzk*gqyz[6]))
                      - szi*(qxxk*gqxx[7]+qyyk*gqyy[7]+qzzk*gqzz[7]
                     +2.0f*(qxyk*gqxy[7]+qxzk*gqxz[7]+qyzk*gqyz[7]))
                      + sxk*(qxxi*gux[11]+qyyi*gux[14]+qzzi*gux[16]
                     +2.0f*(qxyi*gux[12]+qxzi*gux[13]+qyzi*gux[15]))
                      + syk*(qxxi*guy[11]+qyyi*guy[14]+qzzi*guy[16]
                     +2.0f*(qxyi*guy[12]+qxzi*guy[13]+qyzi*guy[15]))
                      + szk*(qxxi*guz[11]+qyyi*guz[14]+qzzi*guz[16]
                     +2.0f*(qxyi*guz[12]+qxzi*guz[13]+qyzi*guz[15]));

    dpdx = 0.5f * (dpsymdx + 0.5f*(dpwidx + dpwkdx));

    dpsymdy = -uxi*(sxk*gux[6]+syk*guy[6]+szk*guz[6])
                          - uyi*(sxk*gux[8]+syk*guy[8]+szk*guz[8])
                          - uzi*(sxk*gux[9]+syk*guy[9]+szk*guz[9])
                          - uxk*(sxi*gux[6]+syi*guy[6]+szi*guz[6])
                          - uyk*(sxi*gux[8]+syi*guy[8]+szi*guz[8])
                          - uzk*(sxi*gux[9]+syi*guy[9]+szi*guz[9]);

    dpwidy = ci*(sxk*gc[6]+syk*gc[8]+szk*gc[9])
                        - ck*(sxi*gux[3]+syi*guy[3]+szi*guz[3])
                         - sxi*(qxxk*gux[12]+qyyk*gux[17]+qzzk*gux[19]
                        +2.0f*(qxyk*gux[14]+qxzk*gux[15]+qyzk*gux[18]))
                         - syi*(qxxk*guy[12]+qyyk*guy[17]+qzzk*guy[19]
                        +2.0f*(qxyk*guy[14]+qxzk*guy[15]+qyzk*guy[18]))
                         - szi*(qxxk*guz[12]+qyyk*guz[17]+qzzk*guz[19]
                        +2.0f*(qxyk*guz[14]+qxzk*guz[15]+qyzk*guz[18]))
                         + sxk*(qxxi*gqxx[6]+qyyi*gqyy[6]+qzzi*gqzz[6]
                        +2.0f*(qxyi*gqxy[6]+qxzi*gqxz[6]+qyzi*gqyz[6]))
                         + syk*(qxxi*gqxx[8]+qyyi*gqyy[8]+qzzi*gqzz[8]
                        +2.0f*(qxyi*gqxy[8]+qxzi*gqxz[8]+qyzi*gqyz[8]))
                         + szk*(qxxi*gqxx[9]+qyyi*gqyy[9]+qzzi*gqzz[9]
                        +2.0f*(qxyi*gqxy[9]+qxzi*gqxz[9]+qyzi*gqyz[9]));

    dpwkdy = ci*(sxk*gux[3]+syk*guy[3]+szk*guz[3])
                        - ck*(sxi*gc[6]+syi*gc[8]+szi*gc[9])
                      - sxi*(qxxk*gqxx[6]+qyyk*gqyy[6]+qzzk*gqzz[6]
                     +2.0f*(qxyk*gqxy[6]+qxzk*gqxz[6]+qyzk*gqyz[6]))
                      - syi*(qxxk*gqxx[8]+qyyk*gqyy[8]+qzzk*gqzz[8]
                     +2.0f*(qxyk*gqxy[8]+qxzk*gqxz[8]+qyzk*gqyz[8]))
                      - szi*(qxxk*gqxx[9]+qyyk*gqyy[9]+qzzk*gqzz[9]
                     +2.0f*(qxyk*gqxy[9]+qxzk*gqxz[9]+qyzk*gqyz[9]))
                      + sxk*(qxxi*gux[12]+qyyi*gux[17]+qzzi*gux[19]
                     +2.0f*(qxyi*gux[14]+qxzi*gux[15]+qyzi*gux[18]))
                      + syk*(qxxi*guy[12]+qyyi*guy[17]+qzzi*guy[19]
                     +2.0f*(qxyi*guy[14]+qxzi*guy[15]+qyzi*guy[18]))
                      + szk*(qxxi*guz[12]+qyyi*guz[17]+qzzi*guz[19]
                     +2.0f*(qxyi*guz[14]+qxzi*guz[15]+qyzi*guz[18]));

    dpdy    = 0.5f * (dpsymdy + 0.5f*(dpwidy + dpwkdy));

    dpsymdz = -uxi*(sxk*gux[7]+syk*guy[7]+szk*guz[7])
                          - uyi*(sxk*gux[9]+syk*guy[9]+szk*guz[9])
                          - uzi*(sxk*gux[10]+syk*guy[10]+szk*guz[10])
                          - uxk*(sxi*gux[7]+syi*guy[7]+szi*guz[7])
                          - uyk*(sxi*gux[9]+syi*guy[9]+szi*guz[9])
                          - uzk*(sxi*gux[10]+syi*guy[10]+szi*guz[10]);

    dpwidz = ci*(sxk*gc[7]+syk*gc[9]+szk*gc[10])
                        - ck*(sxi*gux[4]+syi*guy[4]+szi*guz[4])
                      - sxi*(qxxk*gux[13]+qyyk*gux[18]+qzzk*gux[20]
                     +2.0f*(qxyk*gux[15]+qxzk*gux[16]+qyzk*gux[19]))
                      - syi*(qxxk*guy[13]+qyyk*guy[18]+qzzk*guy[20]
                     +2.0f*(qxyk*guy[15]+qxzk*guy[16]+qyzk*guy[19]))
                      - szi*(qxxk*guz[13]+qyyk*guz[18]+qzzk*guz[20]
                     +2.0f*(qxyk*guz[15]+qxzk*guz[16]+qyzk*guz[19]))
                      + sxk*(qxxi*gqxx[7]+qyyi*gqyy[7]+qzzi*gqzz[7]
                     +2.0f*(qxyi*gqxy[7]+qxzi*gqxz[7]+qyzi*gqyz[7]))
                      + syk*(qxxi*gqxx[9]+qyyi*gqyy[9]+qzzi*gqzz[9]
                     +2.0f*(qxyi*gqxy[9]+qxzi*gqxz[9]+qyzi*gqyz[9]))
                      + szk*(qxxi*gqxx[10]+qyyi*gqyy[10]+qzzi*gqzz[10]
                     +2.0f*(qxyi*gqxy[10]+qxzi*gqxz[10]+qyzi*gqyz[10]));

    dpwkdz = ci*(sxk*gux[4]+syk*guy[4]+szk*guz[4])
                        - ck*(sxi*gc[7]+syi*gc[9]+szi*gc[10])
                      - sxi*(qxxk*gqxx[7]+qyyk*gqyy[7]+qzzk*gqzz[7]
                     +2.0f*(qxyk*gqxy[7]+qxzk*gqxz[7]+qyzk*gqyz[7]))
                      - syi*(qxxk*gqxx[9]+qyyk*gqyy[9]+qzzk*gqzz[9]
                     +2.0f*(qxyk*gqxy[9]+qxzk*gqxz[9]+qyzk*gqyz[9]))
                      - szi*(qxxk*gqxx[10]+qyyk*gqyy[10]+qzzk*gqzz[10]
                     +2.0f*(qxyk*gqxy[10]+qxzk*gqxz[10]+qyzk*gqyz[10]))
                      + sxk*(qxxi*gux[13]+qyyi*gux[18]+qzzi*gux[20]
                     +2.0f*(qxyi*gux[15]+qxzi*gux[16]+qyzi*gux[19]))
                      + syk*(qxxi*guy[13]+qyyi*guy[18]+qzzi*guy[20]
                     +2.0f*(qxyi*guy[15]+qxzi*guy[16]+qyzi*guy[19]))
                      + szk*(qxxi*guz[13]+qyyi*guz[18]+qzzi*guz[20]
                     +2.0f*(qxyi*guz[15]+qxzi*guz[16]+qyzi*guz[19]));

    dpdz = 0.5f * (dpsymdz + 0.5f*(dpwidz + dpwkdz));
 
    // effective radii chain rule terms for the;
    // electrostatic solvation free energy gradient of the permanent;
    // multipoles in the reaction potential of the induced dipoles;
 
    dsymdr = -uxi*(sxk*gux[22]+syk*guy[22]+szk*guz[22])
                          - uyi*(sxk*gux[23]+syk*guy[23]+szk*guz[23])
                          - uzi*(sxk*gux[24]+syk*guy[24]+szk*guz[24])
                          - uxk*(sxi*gux[22]+syi*guy[22]+szi*guz[22])
                          - uyk*(sxi*gux[23]+syi*guy[23]+szi*guz[23])
                          - uzk*(sxi*gux[24]+syi*guy[24]+szi*guz[24]);

    dwipdr = ci*(sxk*gc[22]+syk*gc[23]+szk*gc[24])
                         - ck*(sxi*gux[21]+syi*guy[21]+szi*guz[21])
                      - sxi*(qxxk*gux[25]+qyyk*gux[28]+qzzk*gux[30]
                     +2.0f*(qxyk*gux[26]+qxzk*gux[27]+qyzk*gux[29]))
                      - syi*(qxxk*guy[25]+qyyk*guy[28]+qzzk*guy[30]
                     +2.0f*(qxyk*guy[26]+qxzk*guy[27]+qyzk*guy[29]))
                      - szi*(qxxk*guz[25]+qyyk*guz[28]+qzzk*guz[30]
                     +2.0f*(qxyk*guz[26]+qxzk*guz[27]+qyzk*guz[29]))
                      + sxk*(qxxi*gqxx[22]+qyyi*gqyy[22]+qzzi*gqzz[22]
                     +2.0f*(qxyi*gqxy[22]+qxzi*gqxz[22]+qyzi*gqyz[22]))
                      + syk*(qxxi*gqxx[23]+qyyi*gqyy[23]+qzzi*gqzz[23]
                     +2.0f*(qxyi*gqxy[23]+qxzi*gqxz[23]+qyzi*gqyz[23]))
                      + szk*(qxxi*gqxx[24]+qyyi*gqyy[24]+qzzi*gqzz[24]
                     +2.0f*(qxyi*gqxy[24]+qxzi*gqxz[24]+qyzi*gqyz[24]));

    dwkpdr = ci*(sxk*gux[21]+syk*guy[21]+szk*guz[21])
                         - ck*(sxi*gc[22]+syi*gc[23]+szi*gc[24])
                      - sxi*(qxxk*gqxx[22]+qyyk*gqyy[22]+qzzk*gqzz[22]
                     +2.0f*(qxyk*gqxy[22]+qxzk*gqxz[22]+qyzk*gqyz[22]))
                      - syi*(qxxk*gqxx[23]+qyyk*gqyy[23]+qzzk*gqzz[23]
                     +2.0f*(qxyk*gqxy[23]+qxzk*gqxz[23]+qyzk*gqyz[23]))
                      - szi*(qxxk*gqxx[24]+qyyk*gqyy[24]+qzzk*gqzz[24]
                     +2.0f*(qxyk*gqxy[24]+qxzk*gqxz[24]+qyzk*gqyz[24]))
                      + sxk*(qxxi*gux[25]+qyyi*gux[28]+qzzi*gux[30]
                     +2.0f*(qxyi*gux[26]+qxzi*gux[27]+qyzi*gux[29]))
                      + syk*(qxxi*guy[25]+qyyi*guy[28]+qzzi*guy[30]
                     +2.0f*(qxyi*guy[26]+qxzi*guy[27]+qyzi*guy[29]))
                      + szk*(qxxi*guz[25]+qyyi*guz[28]+qzzi*guz[30]
                     +2.0f*(qxyi*guz[26]+qxzi*guz[27]+qyzi*guz[29]));

    dsumdr = dsymdr + 0.5f*(dwipdr + dwkpdr);
    dpbi = 0.5f*rbk*dsumdr;
    dpbk = 0.5f*rbi*dsumdr;
 
    // mutual polarization electrostatic solvation free energy gradient
 
//   if (poltyp .eq. 'MUTUAL'){

        dpdx = dpdx - 0.5f *
                           (dxi*(pxk*gux[5]+pyk*gux[6]+pzk*gux[7])
                           +dyi*(pxk*guy[5]+pyk*guy[6]+pzk*guy[7])
                           +dzi*(pxk*guz[5]+pyk*guz[6]+pzk*guz[7])
                           +dxk*(pxi*gux[5]+pyi*gux[6]+pzi*gux[7])
                           +dyk*(pxi*guy[5]+pyi*guy[6]+pzi*guy[7])
                           +dzk*(pxi*guz[5]+pyi*guz[6]+pzi*guz[7]));

        dpdy = dpdy - 0.5f *
                           (dxi*(pxk*gux[6]+pyk*gux[8]+pzk*gux[9])
                           +dyi*(pxk*guy[6]+pyk*guy[8]+pzk*guy[9])
                           +dzi*(pxk*guz[6]+pyk*guz[8]+pzk*guz[9])
                           +dxk*(pxi*gux[6]+pyi*gux[8]+pzi*gux[9])
                           +dyk*(pxi*guy[6]+pyi*guy[8]+pzi*guy[9])
                           +dzk*(pxi*guz[6]+pyi*guz[8]+pzi*guz[9]));

        dpdz = dpdz - 0.5f *
                           (dxi*(pxk*gux[7]+pyk*gux[9]+pzk*gux[10])
                           +dyi*(pxk*guy[7]+pyk*guy[9]+pzk*guy[10])
                           +dzi*(pxk*guz[7]+pyk*guz[9]+pzk*guz[10])
                           +dxk*(pxi*gux[7]+pyi*gux[9]+pzi*gux[10])
                           +dyk*(pxi*guy[7]+pyi*guy[9]+pzi*guy[10])
                           +dzk*(pxi*guz[7]+pyi*guz[9]+pzi*guz[10]));

        duvdr = dxi*(pxk*gux[22]+pyk*gux[23]+pzk*gux[24])
                            + dyi*(pxk*guy[22]+pyk*guy[23]+pzk*guy[24])
                            + dzi*(pxk*guz[22]+pyk*guz[23]+pzk*guz[24])
                            + dxk*(pxi*gux[22]+pyi*gux[23]+pzi*gux[24])
                            + dyk*(pxi*guy[22]+pyi*guy[23]+pzi*guy[24])
                            + dzk*(pxi*guz[22]+pyi*guz[23]+pzi*guz[24]);
        dpbi = dpbi - 0.5f*rbk*duvdr;
        dpbk = dpbk - 0.5f*rbi*duvdr;
//    }
 
    // torque due to induced reaction field on permanent dipoles

    fid[1] = 0.5f * (sxk*gux[2]+syk*guy[2]+szk*guz[2]);
    fid[2] = 0.5f * (sxk*gux[3]+syk*guy[3]+szk*guz[3]);
    fid[3] = 0.5f * (sxk*gux[4]+syk*guy[4]+szk*guz[4]);
    fkd[1] = 0.5f * (sxi*gux[2]+syi*guy[2]+szi*guz[2]);
    fkd[2] = 0.5f * (sxi*gux[3]+syi*guy[3]+szi*guz[3]);
    fkd[3] = 0.5f * (sxi*gux[4]+syi*guy[4]+szi*guz[4]);


    // the factor 0.5 appears to be included since trqi[1][i] & trqi[1][k]
    // are identical in the Tinker code (inner loop starts at k = i
    // factor not needed here since
/*
    if ( sameAtom )
    {
        fid[1] = 0.5f * fid[1];
        fid[2] = 0.5f * fid[2];
        fid[3] = 0.5f * fid[3];
        fkd[1] = 0.5f * fkd[1];
        fkd[2] = 0.5f * fkd[2];
        fkd[3] = 0.5f * fkd[3];
    }
*/
    trqi[1]   = uyi*fid[3] - uzi*fid[2];
    trqi[2]   = uzi*fid[1] - uxi*fid[3];
    trqi[3]   = uxi*fid[2] - uyi*fid[1];

    trqi_k[1] = uyk*fkd[3] - uzk*fkd[2];
    trqi_k[2] = uzk*fkd[1] - uxk*fkd[3];
    trqi_k[3] = uxk*fkd[2] - uyk*fkd[1];

 
    // torque due to induced reaction field gradient on quadrupoles;

    fidg[1][1] = -0.25f *
                              ( (sxk*gqxx[2]+syk*gqxx[3]+szk*gqxx[4])
                              + (sxk*gux[5]+syk*guy[5]+szk*guz[5]));

    fidg[1][2] = -0.25f *
                              ( (sxk*gqxy[2]+syk*gqxy[3]+szk*gqxy[4])
                              + (sxk*gux[6]+syk*guy[6]+szk*guz[6]));

    fidg[1][3] = -0.25f *
                              ( (sxk*gqxz[2]+syk*gqxz[3]+szk*gqxz[4])
                              + (sxk*gux[7]+syk*guy[7]+szk*guz[7]));

    fidg[2][2] = -0.25f *
                              ( (sxk*gqyy[2]+syk*gqyy[3]+szk*gqyy[4])
                              + (sxk*gux[8]+syk*guy[8]+szk*guz[8]));

    fidg[2][3] = -0.25f *
                              ( (sxk*gqyz[2]+syk*gqyz[3]+szk*gqyz[4])
                              + (sxk*gux[9]+syk*guy[9]+szk*guz[9]));

    fidg[3][3] = -0.25f *
                              ( (sxk*gqzz[2]+syk*gqzz[3]+szk*gqzz[4])
                              + (sxk*gux[10]+syk*guy[10]+szk*guz[10]));

    fidg[2][1] = fidg[1][2];
    fidg[3][1] = fidg[1][3];
    fidg[3][2] = fidg[2][3];

    fkdg[1][1] = 0.25f *
                              ( (sxi*gqxx[2]+syi*gqxx[3]+szi*gqxx[4])
                              + (sxi*gux[5]+syi*guy[5]+szi*guz[5]));

    fkdg[1][2] = 0.25f *
                              ( (sxi*gqxy[2]+syi*gqxy[3]+szi*gqxy[4])
                              + (sxi*gux[6]+syi*guy[6]+szi*guz[6]));
    fkdg[1][3] = 0.25f *
                              ( (sxi*gqxz[2]+syi*gqxz[3]+szi*gqxz[4])
                              + (sxi*gux[7]+syi*guy[7]+szi*guz[7]));
    fkdg[2][2] = 0.25f *
                              ( (sxi*gqyy[2]+syi*gqyy[3]+szi*gqyy[4])
                              + (sxi*gux[8]+syi*guy[8]+szi*guz[8]));
    fkdg[2][3] = 0.25f *
                              ( (sxi*gqyz[2]+syi*gqyz[3]+szi*gqyz[4])
                              + (sxi*gux[9]+syi*guy[9]+szi*guz[9]));
    fkdg[3][3] = 0.25f *
                              ( (sxi*gqzz[2]+syi*gqzz[3]+szi*gqzz[4])
                              + (sxi*gux[10]+syi*guy[10]+szi*guz[10]));
    fkdg[2][1] = fkdg[1][2];
    fkdg[3][1] = fkdg[1][3];
    fkdg[3][2] = fkdg[2][3];

/*
    if ( sameAtom )
    {
        fidg[1][1] = 0.5f * fidg[1][1];
        fidg[1][2] = 0.5f * fidg[1][2];
        fidg[1][3] = 0.5f * fidg[1][3];
        fidg[2][1] = 0.5f * fidg[2][1];
        fidg[2][2] = 0.5f * fidg[2][2];
        fidg[2][3] = 0.5f * fidg[2][3];
        fidg[3][1] = 0.5f * fidg[3][1];
        fidg[3][2] = 0.5f * fidg[3][2];
        fidg[3][3] = 0.5f * fidg[3][3];
        fkdg[1][1] = 0.5f * fkdg[1][1];
        fkdg[1][2] = 0.5f * fkdg[1][2];
        fkdg[1][3] = 0.5f * fkdg[1][3];
        fkdg[2][1] = 0.5f * fkdg[2][1];
        fkdg[2][2] = 0.5f * fkdg[2][2];
        fkdg[2][3] = 0.5f * fkdg[2][3];
        fkdg[3][1] = 0.5f * fkdg[3][1];
        fkdg[3][2] = 0.5f * fkdg[3][2];
        fkdg[3][3] = 0.5f * fkdg[3][3];
     }
*/

    trqi[1] += 2.0f*(qxyi*fidg[1][3]+qyyi*fidg[2][3]+qyzi*fidg[3][3]
                        -qxzi*fidg[1][2]-qyzi*fidg[2][2]-qzzi*fidg[3][2]);
    trqi[2] += 2.0f*(qxzi*fidg[1][1]+qyzi*fidg[2][1]+qzzi*fidg[3][1]
                        -qxxi*fidg[1][3]-qxyi*fidg[2][3]-qxzi*fidg[3][3]);

    trqi[3] += 2.0f*(qxxi*fidg[1][2]+qxyi*fidg[2][2]+qxzi*fidg[3][2]
                        -qxyi*fidg[1][1]-qyyi*fidg[2][1]-qyzi*fidg[3][1]);

    trqi_k[1] += 2.0f*
                        (qxyk*fkdg[1][3]+qyyk*fkdg[2][3]+qyzk*fkdg[3][3]
                        -qxzk*fkdg[1][2]-qyzk*fkdg[2][2]-qzzk*fkdg[3][2]);

    trqi_k[2] += 2.0f*
                        (qxzk*fkdg[1][1]+qyzk*fkdg[2][1]+qzzk*fkdg[3][1]
                        -qxxk*fkdg[1][3]-qxyk*fkdg[2][3]-qxzk*fkdg[3][3]);

    trqi_k[3] += 2.0f*
                        (qxxk*fkdg[1][2]+qxyk*fkdg[2][2]+qxzk*fkdg[3][2]
                        -qxyk*fkdg[1][1]-qyyk*fkdg[2][1]-qyzk*fkdg[3][1]);
 
#if 1
#ifdef AMOEBA_DEBUG
if( 1 ){
/*
int debugIndex               = atomJ;
    debugArray[debugIndex].x = rbi;
    debugArray[debugIndex].y = rbk;
    debugArray[debugIndex].z = esym;
    debugArray[debugIndex].w = ewi;
*/
    debugArray[0].x          = dedx;
    debugArray[0].y          = dedy;
    debugArray[0].z          = dedz;
/*
    dewidx = ci*(uxk*gc[5]+uyk*gc[6]+uzk*gc[7])
                      -ck*(uxi*gux[2]+uyi*guy[2]+uzi*guz[2])
                 +ci*(qxxk*gc[11]+qyyk*gc[14]+qzzk*gc[16]
              +2.0f*(qxyk*gc[12]+qxzk*gc[13]+qyzk*gc[15]))
 xr * yr * zr * a[0][3];
*/

    debugArray[1].x          = 5.0f;
    debugArray[1].y          = expc*powf( dexpc, 2.0f );
    debugArray[1].z          = powf( dexpc, 2.0f );
    debugArray[1].w          = dexpc*dexpc;

    //a[0][3]      = expc1*a[1][2] + expcdexpc*a[1][1];
    //expcdexpc    = -expc * dexpc;
    debugArray[2].x          = expc;
    debugArray[2].y          = dexpc;
    debugArray[2].z          = expcdexpc;
    debugArray[2].w          = a[1][1];
/*
    debugIndex              += numberOfAtoms;
    debugArray[debugIndex].x = trq[1] + trqi[1];
    debugArray[debugIndex].y = trq[2] + trqi[2];
    debugArray[debugIndex].z = trq[3] + trqi[3];

    debugIndex              += numberOfAtoms;
    debugArray[debugIndex].x = drbi;
    debugArray[debugIndex].y = dpbi;

*/
/*
    debugIndex              += numberOfAtoms;
    debugArray[debugIndex].x = dedx;
    debugArray[debugIndex].y = dedy;
    debugArray[debugIndex].z = dedz;

    debugIndex              += numberOfAtoms;
    debugArray[debugIndex].x = dpdx;
    debugArray[debugIndex].y = dpdy;
    debugArray[debugIndex].z = dpdz;

    debugIndex              += numberOfAtoms;
    debugArray[debugIndex].x = trq[1];
    debugArray[debugIndex].y = trq[2];
    debugArray[debugIndex].z = trq[3];

    debugIndex              += numberOfAtoms;
    debugArray[debugIndex].x = trqi[1];
    debugArray[debugIndex].y = trqi[2];
    debugArray[debugIndex].z = trqi[3];

    debugIndex              += numberOfAtoms;
    debugArray[debugIndex].x = dpsymdx;
    debugArray[debugIndex].y = dpwidx;
    debugArray[debugIndex].z = dpwkdx;

    debugIndex              += numberOfAtoms;
    debugArray[debugIndex].x = dpsymdy;
    debugArray[debugIndex].y = dpwidy;
    debugArray[debugIndex].z = dpwkdy;

    debugIndex              += numberOfAtoms;
    debugArray[debugIndex].x = dxi;
    debugArray[debugIndex].y = dyi;
    debugArray[debugIndex].z = dzi;

    debugIndex              += numberOfAtoms;
    debugArray[debugIndex].x = dxk;
    debugArray[debugIndex].y = dyk;
    debugArray[debugIndex].z = dzk;

    debugIndex              += numberOfAtoms;
    debugArray[debugIndex].x = pxi;
    debugArray[debugIndex].y = pyi;
    debugArray[debugIndex].z = pzi;

    debugIndex              += numberOfAtoms;
    debugArray[debugIndex].x = pxk;
    debugArray[debugIndex].y = pyk;
    debugArray[debugIndex].z = pzk;

    debugIndex              += numberOfAtoms;
    debugArray[debugIndex].x = fid[1];
    debugArray[debugIndex].y = fid[2];
    debugArray[debugIndex].z = fid[3];

    debugIndex              += numberOfAtoms;
    debugArray[debugIndex].x = fidg[1][3];
    debugArray[debugIndex].y = fidg[2][3];
    debugArray[debugIndex].z = fidg[3][3];

*/

}
#endif
#endif

    // total permanent and induced energies for this interaction;

    e                        = esym + 0.5f*(ewi+ewk);
    ei                       = 0.5f * (esymi + 0.5f*(ewii+ewki));
 
    outputForce[0]           = (dedx + dpdx);
    outputForce[1]           = (dedy + dpdy);
    outputForce[2]           = (dedz + dpdz);
    
    outputTorque[0][0]       = (trq[1] + trqi[1]);
    outputTorque[0][1]       = (trq[2] + trqi[2]);
    outputTorque[0][2]       = (trq[3] + trqi[3]);

    outputTorque[1][0]       = (trq_k[1] + trqi_k[1]);
    outputTorque[1][1]       = (trq_k[2] + trqi_k[2]);
    outputTorque[1][2]       = (trq_k[3] + trqi_k[3]);

    outputBorn[0]            = drbi;
    outputBorn[1]            = drbk;

    outputBornPolar[0]       = dpbi;
    outputBornPolar[1]       = dpbk;

    *outputEnergy            = (e + ei);
 
}

__device__ static int debugAccumulate( int index, float4* debugArray, float* field, unsigned int addMask, float idLabel )
{
    index                             += cAmoebaSim.paddedNumberOfAtoms;
    debugArray[index].x                = addMask ? field[0] : 0.0f;
    debugArray[index].y                = addMask ? field[1] : 0.0f;
    debugArray[index].z                = addMask ? field[2] : 0.0f;
    debugArray[index].w                = idLabel;

    return index;
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
        threadsPerBlock                             = getThreadsPerBlock( amoebaGpu, sizeof(KirkwoodParticle));
threadsPerBlock = 32;
        //unsigned int eDiffhreadsPerBlock            = getThreadsPerBlock( amoebaGpu, sizeof(KirkwoodEDiffParticle));
        //unsigned int maxThreadsPerBlock             = threadsPerBlock> eDiffhreadsPerBlock ? threadsPerBlock : eDiffhreadsPerBlock;

        if( amoebaGpu->log ){
       
#if (__CUDA_ARCH__ >= 200)
            unsigned int maxThreads = GF1XX_NONBOND_THREADS_PER_BLOCK;
#elif (__CUDA_ARCH__ >= 130)
            unsigned int maxThreads = GT2XX_NONBOND_THREADS_PER_BLOCK;
#else
            unsigned int maxThreads = G8X_NONBOND_THREADS_PER_BLOCK;
#endif

            (void) fprintf( amoebaGpu->log, "kCalculateAmoebaCudaKirkwood: blcks=%u tds=%u %u bPrWrp=%u atm=%u shrd=%u Ebuf=%u ixnCt=%u workUnits=%u\n",
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
