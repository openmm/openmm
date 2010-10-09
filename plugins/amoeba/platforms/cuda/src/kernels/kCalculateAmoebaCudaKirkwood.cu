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

__device__ void calculateKirkwoodPairIxn_kernel( KirkwoodParticle& atomI,       KirkwoodParticle& atomJ,
                                                 unsigned int sameAtom,
                                                 float*  outputForce,           float outputTorque[2][3],
                                                 float*  outputBorn,            float*  outputBornPolar,
                                                 float* outputEnergy
#ifdef AMOEBA_DEBUG
                                                 , float4* debugArray
#endif

 ){

    float e,ei;
    float xr,yr,zr;
    float xr2,yr2,zr2;
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

    // set the bulk dielectric constant to the water value

    fc           = cAmoebaSim.electric * cAmoebaSim.fc;
    fd           = cAmoebaSim.electric * cAmoebaSim.fd;
    fq           = cAmoebaSim.electric * cAmoebaSim.fq;

    sxi          = atomI.inducedDipole[0] + atomI.inducedDipoleP[0];
    syi          = atomI.inducedDipole[1] + atomI.inducedDipoleP[1];
    szi          = atomI.inducedDipole[2] + atomI.inducedDipoleP[2];

    // decide whether to compute the current interaction;

    xr           = atomJ.x - atomI.x;
    yr           = atomJ.y - atomI.y;
    zr           = atomJ.z - atomI.z;

    xr2          = xr*xr;
    yr2          = yr*yr;
    zr2          = zr*zr;
    r2           = xr2 + yr2 + zr2;

    if( r2 > cAmoebaSim.scalingDistanceCutoff ){
    }

    sxk          = atomJ.inducedDipole[0] + atomJ.inducedDipoleP[0];
    syk          = atomJ.inducedDipole[1] + atomJ.inducedDipoleP[1];
    szk          = atomJ.inducedDipole[2] + atomJ.inducedDipoleP[2];
    rb2          = atomI.bornRadius * atomJ.bornRadius;
    expterm      = expf(-r2/(cAmoebaSim.gkc*rb2));
    expc         = expterm / cAmoebaSim.gkc;
    expcr        = r2*expterm / (cAmoebaSim.gkc*cAmoebaSim.gkc*rb2*rb2);
    dexpc        = -2.0f / (cAmoebaSim.gkc*rb2);
    dexpcr       = 2.0f / (cAmoebaSim.gkc*rb2*rb2);
    dgfdr        = 0.5f * expterm * (1.0f+r2/(rb2*cAmoebaSim.gkc));
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

    esym = atomI.q * atomJ.q * gc1 - (atomI.labFrameDipole[0]*(atomJ.labFrameDipole[0]*gux2+atomJ.labFrameDipole[1]*guy2+atomJ.labFrameDipole[2]*guz2)
                           +  atomI.labFrameDipole[1]*(atomJ.labFrameDipole[0]*gux3+atomJ.labFrameDipole[1]*guy3+atomJ.labFrameDipole[2]*guz3)
                           +  atomI.labFrameDipole[2]*(atomJ.labFrameDipole[0]*gux4+atomJ.labFrameDipole[1]*guy4+atomJ.labFrameDipole[2]*guz4));

    ewi =  atomI.q*(atomJ.labFrameDipole[0]*gc2+atomJ.labFrameDipole[1]*gc3+atomJ.labFrameDipole[2]*gc4)
          -atomJ.q*(atomI.labFrameDipole[0]*gux1+atomI.labFrameDipole[1]*guy1+atomI.labFrameDipole[2]*guz1)
           +atomI.q*(atomJ.labFrameQuadrupole_XX*gc5+atomJ.labFrameQuadrupole_YY*gc8+atomJ.labFrameQuadrupole_ZZ*gc10
              +2.0f*(atomJ.labFrameQuadrupole_XY*gc6+atomJ.labFrameQuadrupole_XZ*gc7+atomJ.labFrameQuadrupole_YZ*gc9))
                 +atomJ.q*(atomI.labFrameQuadrupole_XX*gqxx1+atomI.labFrameQuadrupole_YY*gqyy1+atomI.labFrameQuadrupole_ZZ*gqzz1
              +2.0f*(atomI.labFrameQuadrupole_XY*gqxy1+atomI.labFrameQuadrupole_XZ*gqxz1+atomI.labFrameQuadrupole_YZ*gqyz1))
               - atomI.labFrameDipole[0]*(atomJ.labFrameQuadrupole_XX*gux5+atomJ.labFrameQuadrupole_YY*gux8+atomJ.labFrameQuadrupole_ZZ*gux10
              +2.0f*(atomJ.labFrameQuadrupole_XY*gux6+atomJ.labFrameQuadrupole_XZ*gux7+atomJ.labFrameQuadrupole_YZ*gux9))
               - atomI.labFrameDipole[1]*(atomJ.labFrameQuadrupole_XX*guy5+atomJ.labFrameQuadrupole_YY*guy8+atomJ.labFrameQuadrupole_ZZ*guy10
              +2.0f*(atomJ.labFrameQuadrupole_XY*guy6+atomJ.labFrameQuadrupole_XZ*guy7+atomJ.labFrameQuadrupole_YZ*guy9))
               - atomI.labFrameDipole[2]*(atomJ.labFrameQuadrupole_XX*guz5+atomJ.labFrameQuadrupole_YY*guz8+atomJ.labFrameQuadrupole_ZZ*guz10
              +2.0f*(atomJ.labFrameQuadrupole_XY*guz6+atomJ.labFrameQuadrupole_XZ*guz7+atomJ.labFrameQuadrupole_YZ*guz9))
               + atomJ.labFrameDipole[0]*(atomI.labFrameQuadrupole_XX*gqxx2+atomI.labFrameQuadrupole_YY*gqyy2+atomI.labFrameQuadrupole_ZZ*gqzz2
              +2.0f*(atomI.labFrameQuadrupole_XY*gqxy2+atomI.labFrameQuadrupole_XZ*gqxz2+atomI.labFrameQuadrupole_YZ*gqyz2))
               + atomJ.labFrameDipole[1]*(atomI.labFrameQuadrupole_XX*gqxx3+atomI.labFrameQuadrupole_YY*gqyy3+atomI.labFrameQuadrupole_ZZ*gqzz3
              +2.0f*(atomI.labFrameQuadrupole_XY*gqxy3+atomI.labFrameQuadrupole_XZ*gqxz3+atomI.labFrameQuadrupole_YZ*gqyz3))
               + atomJ.labFrameDipole[2]*(atomI.labFrameQuadrupole_XX*gqxx4+atomI.labFrameQuadrupole_YY*gqyy4+atomI.labFrameQuadrupole_ZZ*gqzz4
              +2.0f*(atomI.labFrameQuadrupole_XY*gqxy4+atomI.labFrameQuadrupole_XZ*gqxz4+atomI.labFrameQuadrupole_YZ*gqyz4))
              + atomI.labFrameQuadrupole_XX*(atomJ.labFrameQuadrupole_XX*gqxx5+atomJ.labFrameQuadrupole_YY*gqxx8+atomJ.labFrameQuadrupole_ZZ*gqxx10
              +2.0f*(atomJ.labFrameQuadrupole_XY*gqxx6+atomJ.labFrameQuadrupole_XZ*gqxx7+atomJ.labFrameQuadrupole_YZ*gqxx9))
              + atomI.labFrameQuadrupole_YY*(atomJ.labFrameQuadrupole_XX*gqyy5+atomJ.labFrameQuadrupole_YY*gqyy8+atomJ.labFrameQuadrupole_ZZ*gqyy10
              +2.0f*(atomJ.labFrameQuadrupole_XY*gqyy6+atomJ.labFrameQuadrupole_XZ*gqyy7+atomJ.labFrameQuadrupole_YZ*gqyy9))
              + atomI.labFrameQuadrupole_ZZ*(atomJ.labFrameQuadrupole_XX*gqzz5+atomJ.labFrameQuadrupole_YY*gqzz8+atomJ.labFrameQuadrupole_ZZ*gqzz10
              +2.0f*(atomJ.labFrameQuadrupole_XY*gqzz6+atomJ.labFrameQuadrupole_XZ*gqzz7+atomJ.labFrameQuadrupole_YZ*gqzz9))
              + 2.0f*(atomI.labFrameQuadrupole_XY*(atomJ.labFrameQuadrupole_XX*gqxy5+atomJ.labFrameQuadrupole_YY*gqxy8+atomJ.labFrameQuadrupole_ZZ*gqxy10
              +2.0f*(atomJ.labFrameQuadrupole_XY*gqxy6+atomJ.labFrameQuadrupole_XZ*gqxy7+atomJ.labFrameQuadrupole_YZ*gqxy9))
              + atomI.labFrameQuadrupole_XZ*(atomJ.labFrameQuadrupole_XX*gqxz5+atomJ.labFrameQuadrupole_YY*gqxz8+atomJ.labFrameQuadrupole_ZZ*gqxz10
              +2.0f*(atomJ.labFrameQuadrupole_XY*gqxz6+atomJ.labFrameQuadrupole_XZ*gqxz7+atomJ.labFrameQuadrupole_YZ*gqxz9))
              + atomI.labFrameQuadrupole_YZ*(atomJ.labFrameQuadrupole_XX*gqyz5+atomJ.labFrameQuadrupole_YY*gqyz8+atomJ.labFrameQuadrupole_ZZ*gqyz10
              +2.0f*(atomJ.labFrameQuadrupole_XY*gqyz6+atomJ.labFrameQuadrupole_XZ*gqyz7+atomJ.labFrameQuadrupole_YZ*gqyz9)));

    ewk = atomI.q*(atomJ.labFrameDipole[0]*gux1+atomJ.labFrameDipole[1]*guy1+atomJ.labFrameDipole[2]*guz1)
                      -atomJ.q*(atomI.labFrameDipole[0]*gc2+atomI.labFrameDipole[1]*gc3+atomI.labFrameDipole[2]*gc4)
                 +atomI.q*(atomJ.labFrameQuadrupole_XX*gqxx1+atomJ.labFrameQuadrupole_YY*gqyy1+atomJ.labFrameQuadrupole_ZZ*gqzz1
              +2.0f*(atomJ.labFrameQuadrupole_XY*gqxy1+atomJ.labFrameQuadrupole_XZ*gqxz1+atomJ.labFrameQuadrupole_YZ*gqyz1))
                 +atomJ.q*(atomI.labFrameQuadrupole_XX*gc5+atomI.labFrameQuadrupole_YY*gc8+atomI.labFrameQuadrupole_ZZ*gc10
              +2.0f*(atomI.labFrameQuadrupole_XY*gc6+atomI.labFrameQuadrupole_XZ*gc7+atomI.labFrameQuadrupole_YZ*gc9))
               - atomI.labFrameDipole[0]*(atomJ.labFrameQuadrupole_XX*gqxx2+atomJ.labFrameQuadrupole_YY*gqyy2+atomJ.labFrameQuadrupole_ZZ*gqzz2
              +2.0f*(atomJ.labFrameQuadrupole_XY*gqxy2+atomJ.labFrameQuadrupole_XZ*gqxz2+atomJ.labFrameQuadrupole_YZ*gqyz2))
               - atomI.labFrameDipole[1]*(atomJ.labFrameQuadrupole_XX*gqxx3+atomJ.labFrameQuadrupole_YY*gqyy3+atomJ.labFrameQuadrupole_ZZ*gqzz3
              +2.0f*(atomJ.labFrameQuadrupole_XY*gqxy3+atomJ.labFrameQuadrupole_XZ*gqxz3+atomJ.labFrameQuadrupole_YZ*gqyz3))
               - atomI.labFrameDipole[2]*(atomJ.labFrameQuadrupole_XX*gqxx4+atomJ.labFrameQuadrupole_YY*gqyy4+atomJ.labFrameQuadrupole_ZZ*gqzz4
              +2.0f*(atomJ.labFrameQuadrupole_XY*gqxy4+atomJ.labFrameQuadrupole_XZ*gqxz4+atomJ.labFrameQuadrupole_YZ*gqyz4))
               + atomJ.labFrameDipole[0]*(atomI.labFrameQuadrupole_XX*gux5+atomI.labFrameQuadrupole_YY*gux8+atomI.labFrameQuadrupole_ZZ*gux10
              +2.0f*(atomI.labFrameQuadrupole_XY*gux6+atomI.labFrameQuadrupole_XZ*gux7+atomI.labFrameQuadrupole_YZ*gux9))
               + atomJ.labFrameDipole[1]*(atomI.labFrameQuadrupole_XX*guy5+atomI.labFrameQuadrupole_YY*guy8+atomI.labFrameQuadrupole_ZZ*guy10
              +2.0f*(atomI.labFrameQuadrupole_XY*guy6+atomI.labFrameQuadrupole_XZ*guy7+atomI.labFrameQuadrupole_YZ*guy9))
               + atomJ.labFrameDipole[2]*(atomI.labFrameQuadrupole_XX*guz5+atomI.labFrameQuadrupole_YY*guz8+atomI.labFrameQuadrupole_ZZ*guz10
              +2.0f*(atomI.labFrameQuadrupole_XY*guz6+atomI.labFrameQuadrupole_XZ*guz7+atomI.labFrameQuadrupole_YZ*guz9))
              + atomI.labFrameQuadrupole_XX*(atomJ.labFrameQuadrupole_XX*gqxx5+atomJ.labFrameQuadrupole_YY*gqyy5+atomJ.labFrameQuadrupole_ZZ*gqzz5
              +2.0f*(atomJ.labFrameQuadrupole_XY*gqxy5+atomJ.labFrameQuadrupole_XZ*gqxz5+atomJ.labFrameQuadrupole_YZ*gqyz5))
              + atomI.labFrameQuadrupole_YY*(atomJ.labFrameQuadrupole_XX*gqxx8+atomJ.labFrameQuadrupole_YY*gqyy8+atomJ.labFrameQuadrupole_ZZ*gqzz8
              +2.0f*(atomJ.labFrameQuadrupole_XY*gqxy8+atomJ.labFrameQuadrupole_XZ*gqxz8+atomJ.labFrameQuadrupole_YZ*gqyz8))
              + atomI.labFrameQuadrupole_ZZ*(atomJ.labFrameQuadrupole_XX*gqxx10+atomJ.labFrameQuadrupole_YY*gqyy10+atomJ.labFrameQuadrupole_ZZ*gqzz10
              +2.0f*(atomJ.labFrameQuadrupole_XY*gqxy10+atomJ.labFrameQuadrupole_XZ*gqxz10+atomJ.labFrameQuadrupole_YZ*gqyz10))
       + 2.0f*(atomI.labFrameQuadrupole_XY*(atomJ.labFrameQuadrupole_XX*gqxx6+atomJ.labFrameQuadrupole_YY*gqyy6+atomJ.labFrameQuadrupole_ZZ*gqzz6
              +2.0f*(atomJ.labFrameQuadrupole_XY*gqxy6+atomJ.labFrameQuadrupole_XZ*gqxz6+atomJ.labFrameQuadrupole_YZ*gqyz6))
              + atomI.labFrameQuadrupole_XZ*(atomJ.labFrameQuadrupole_XX*gqxx7+atomJ.labFrameQuadrupole_YY*gqyy7+atomJ.labFrameQuadrupole_ZZ*gqzz7
              +2.0f*(atomJ.labFrameQuadrupole_XY*gqxy7+atomJ.labFrameQuadrupole_XZ*gqxz7+atomJ.labFrameQuadrupole_YZ*gqyz7))
              + atomI.labFrameQuadrupole_YZ*(atomJ.labFrameQuadrupole_XX*gqxx9+atomJ.labFrameQuadrupole_YY*gqyy9+atomJ.labFrameQuadrupole_ZZ*gqzz9
              +2.0f*(atomJ.labFrameQuadrupole_XY*gqxy9+atomJ.labFrameQuadrupole_XZ*gqxz9+atomJ.labFrameQuadrupole_YZ*gqyz9)));

    desymdx = atomI.q * atomJ.q * gc2 - (atomI.labFrameDipole[0]*(atomJ.labFrameDipole[0]*gux5+atomJ.labFrameDipole[1]*guy5+atomJ.labFrameDipole[2]*guz5)
                              +  atomI.labFrameDipole[1]*(atomJ.labFrameDipole[0]*gux6+atomJ.labFrameDipole[1]*guy6+atomJ.labFrameDipole[2]*guz6)
                              +  atomI.labFrameDipole[2]*(atomJ.labFrameDipole[0]*gux7+atomJ.labFrameDipole[1]*guy7+atomJ.labFrameDipole[2]*guz7));

    dewidx = atomI.q*(atomJ.labFrameDipole[0]*gc5+atomJ.labFrameDipole[1]*gc6+atomJ.labFrameDipole[2]*gc7)
                      -atomJ.q*(atomI.labFrameDipole[0]*gux2+atomI.labFrameDipole[1]*guy2+atomI.labFrameDipole[2]*guz2)
                 +atomI.q*(atomJ.labFrameQuadrupole_XX*gc11+atomJ.labFrameQuadrupole_YY*gc14+atomJ.labFrameQuadrupole_ZZ*gc16
              +2.0f*(atomJ.labFrameQuadrupole_XY*gc12+atomJ.labFrameQuadrupole_XZ*gc13+atomJ.labFrameQuadrupole_YZ*gc15))
                 +atomJ.q*(atomI.labFrameQuadrupole_XX*gqxx2+atomI.labFrameQuadrupole_YY*gqyy2+atomI.labFrameQuadrupole_ZZ*gqzz2
              +2.0f*(atomI.labFrameQuadrupole_XY*gqxy2+atomI.labFrameQuadrupole_XZ*gqxz2+atomI.labFrameQuadrupole_YZ*gqyz2))
               - atomI.labFrameDipole[0]*(atomJ.labFrameQuadrupole_XX*gux11+atomJ.labFrameQuadrupole_YY*gux14+atomJ.labFrameQuadrupole_ZZ*gux16
              +2.0f*(atomJ.labFrameQuadrupole_XY*gux12+atomJ.labFrameQuadrupole_XZ*gux13+atomJ.labFrameQuadrupole_YZ*gux15))
               - atomI.labFrameDipole[1]*(atomJ.labFrameQuadrupole_XX*guy11+atomJ.labFrameQuadrupole_YY*guy14+atomJ.labFrameQuadrupole_ZZ*guy16
              +2.0f*(atomJ.labFrameQuadrupole_XY*guy12+atomJ.labFrameQuadrupole_XZ*guy13+atomJ.labFrameQuadrupole_YZ*guy15))
               - atomI.labFrameDipole[2]*(atomJ.labFrameQuadrupole_XX*guz11+atomJ.labFrameQuadrupole_YY*guz14+atomJ.labFrameQuadrupole_ZZ*guz16
              +2.0f*(atomJ.labFrameQuadrupole_XY*guz12+atomJ.labFrameQuadrupole_XZ*guz13+atomJ.labFrameQuadrupole_YZ*guz15))
               + atomJ.labFrameDipole[0]*(atomI.labFrameQuadrupole_XX*gqxx5+atomI.labFrameQuadrupole_YY*gqyy5+atomI.labFrameQuadrupole_ZZ*gqzz5
              +2.0f*(atomI.labFrameQuadrupole_XY*gqxy5+atomI.labFrameQuadrupole_XZ*gqxz5+atomI.labFrameQuadrupole_YZ*gqyz5))
               + atomJ.labFrameDipole[1]*(atomI.labFrameQuadrupole_XX*gqxx6+atomI.labFrameQuadrupole_YY*gqyy6+atomI.labFrameQuadrupole_ZZ*gqzz6
              +2.0f*(atomI.labFrameQuadrupole_XY*gqxy6+atomI.labFrameQuadrupole_XZ*gqxz6+atomI.labFrameQuadrupole_YZ*gqyz6))
               + atomJ.labFrameDipole[2]*(atomI.labFrameQuadrupole_XX*gqxx7+atomI.labFrameQuadrupole_YY*gqyy7+atomI.labFrameQuadrupole_ZZ*gqzz7
              +2.0f*(atomI.labFrameQuadrupole_XY*gqxy7+atomI.labFrameQuadrupole_XZ*gqxz7+atomI.labFrameQuadrupole_YZ*gqyz7))
              + atomI.labFrameQuadrupole_XX*(atomJ.labFrameQuadrupole_XX*gqxx11+atomJ.labFrameQuadrupole_YY*gqxx14+atomJ.labFrameQuadrupole_ZZ*gqxx16
              +2.0f*(atomJ.labFrameQuadrupole_XY*gqxx12+atomJ.labFrameQuadrupole_XZ*gqxx13+atomJ.labFrameQuadrupole_YZ*gqxx15))
              + atomI.labFrameQuadrupole_YY*(atomJ.labFrameQuadrupole_XX*gqyy11+atomJ.labFrameQuadrupole_YY*gqyy14+atomJ.labFrameQuadrupole_ZZ*gqyy16
              +2.0f*(atomJ.labFrameQuadrupole_XY*gqyy12+atomJ.labFrameQuadrupole_XZ*gqyy13+atomJ.labFrameQuadrupole_YZ*gqyy15))
              + atomI.labFrameQuadrupole_ZZ*(atomJ.labFrameQuadrupole_XX*gqzz11+atomJ.labFrameQuadrupole_YY*gqzz14+atomJ.labFrameQuadrupole_ZZ*gqzz16
              +2.0f*(atomJ.labFrameQuadrupole_XY*gqzz12+atomJ.labFrameQuadrupole_XZ*gqzz13+atomJ.labFrameQuadrupole_YZ*gqzz15))
       + 2.0f*(atomI.labFrameQuadrupole_XY*(atomJ.labFrameQuadrupole_XX*gqxy11+atomJ.labFrameQuadrupole_YY*gqxy14+atomJ.labFrameQuadrupole_ZZ*gqxy16
              +2.0f*(atomJ.labFrameQuadrupole_XY*gqxy12+atomJ.labFrameQuadrupole_XZ*gqxy13+atomJ.labFrameQuadrupole_YZ*gqxy15))
              + atomI.labFrameQuadrupole_XZ*(atomJ.labFrameQuadrupole_XX*gqxz11+atomJ.labFrameQuadrupole_YY*gqxz14+atomJ.labFrameQuadrupole_ZZ*gqxz16
              +2.0f*(atomJ.labFrameQuadrupole_XY*gqxz12+atomJ.labFrameQuadrupole_XZ*gqxz13+atomJ.labFrameQuadrupole_YZ*gqxz15))
              + atomI.labFrameQuadrupole_YZ*(atomJ.labFrameQuadrupole_XX*gqyz11+atomJ.labFrameQuadrupole_YY*gqyz14+atomJ.labFrameQuadrupole_ZZ*gqyz16
              +2.0f*(atomJ.labFrameQuadrupole_XY*gqyz12+atomJ.labFrameQuadrupole_XZ*gqyz13+atomJ.labFrameQuadrupole_YZ*gqyz15)));

    dewkdx = atomI.q*(atomJ.labFrameDipole[0]*gux2+atomJ.labFrameDipole[1]*guy2+atomJ.labFrameDipole[2]*guz2)
                      -atomJ.q*(atomI.labFrameDipole[0]*gc5+atomI.labFrameDipole[1]*gc6+atomI.labFrameDipole[2]*gc7)
                 +atomI.q*(atomJ.labFrameQuadrupole_XX*gqxx2+atomJ.labFrameQuadrupole_YY*gqyy2+atomJ.labFrameQuadrupole_ZZ*gqzz2
              +2.0f*(atomJ.labFrameQuadrupole_XY*gqxy2+atomJ.labFrameQuadrupole_XZ*gqxz2+atomJ.labFrameQuadrupole_YZ*gqyz2))
                 +atomJ.q*(atomI.labFrameQuadrupole_XX*gc11+atomI.labFrameQuadrupole_YY*gc14+atomI.labFrameQuadrupole_ZZ*gc16
              +2.0f*(atomI.labFrameQuadrupole_XY*gc12+atomI.labFrameQuadrupole_XZ*gc13+atomI.labFrameQuadrupole_YZ*gc15))
               - atomI.labFrameDipole[0]*(atomJ.labFrameQuadrupole_XX*gqxx5+atomJ.labFrameQuadrupole_YY*gqyy5+atomJ.labFrameQuadrupole_ZZ*gqzz5
              +2.0f*(atomJ.labFrameQuadrupole_XY*gqxy5+atomJ.labFrameQuadrupole_XZ*gqxz5+atomJ.labFrameQuadrupole_YZ*gqyz5))
               - atomI.labFrameDipole[1]*(atomJ.labFrameQuadrupole_XX*gqxx6+atomJ.labFrameQuadrupole_YY*gqyy6+atomJ.labFrameQuadrupole_ZZ*gqzz6
              +2.0f*(atomJ.labFrameQuadrupole_XY*gqxy6+atomJ.labFrameQuadrupole_XZ*gqxz6+atomJ.labFrameQuadrupole_YZ*gqyz6))
               - atomI.labFrameDipole[2]*(atomJ.labFrameQuadrupole_XX*gqxx7+atomJ.labFrameQuadrupole_YY*gqyy7+atomJ.labFrameQuadrupole_ZZ*gqzz7
              +2.0f*(atomJ.labFrameQuadrupole_XY*gqxy7+atomJ.labFrameQuadrupole_XZ*gqxz7+atomJ.labFrameQuadrupole_YZ*gqyz7))
               + atomJ.labFrameDipole[0]*(atomI.labFrameQuadrupole_XX*gux11+atomI.labFrameQuadrupole_YY*gux14+atomI.labFrameQuadrupole_ZZ*gux16
              +2.0f*(atomI.labFrameQuadrupole_XY*gux12+atomI.labFrameQuadrupole_XZ*gux13+atomI.labFrameQuadrupole_YZ*gux15))
               + atomJ.labFrameDipole[1]*(atomI.labFrameQuadrupole_XX*guy11+atomI.labFrameQuadrupole_YY*guy14+atomI.labFrameQuadrupole_ZZ*guy16
              +2.0f*(atomI.labFrameQuadrupole_XY*guy12+atomI.labFrameQuadrupole_XZ*guy13+atomI.labFrameQuadrupole_YZ*guy15))
               + atomJ.labFrameDipole[2]*(atomI.labFrameQuadrupole_XX*guz11+atomI.labFrameQuadrupole_YY*guz14+atomI.labFrameQuadrupole_ZZ*guz16
              +2.0f*(atomI.labFrameQuadrupole_XY*guz12+atomI.labFrameQuadrupole_XZ*guz13+atomI.labFrameQuadrupole_YZ*guz15))
              + atomI.labFrameQuadrupole_XX*(atomJ.labFrameQuadrupole_XX*gqxx11+atomJ.labFrameQuadrupole_YY*gqyy11+atomJ.labFrameQuadrupole_ZZ*gqzz11
              +2.0f*(atomJ.labFrameQuadrupole_XY*gqxy11+atomJ.labFrameQuadrupole_XZ*gqxz11+atomJ.labFrameQuadrupole_YZ*gqyz11))
              + atomI.labFrameQuadrupole_YY*(atomJ.labFrameQuadrupole_XX*gqxx14+atomJ.labFrameQuadrupole_YY*gqyy14+atomJ.labFrameQuadrupole_ZZ*gqzz14
              +2.0f*(atomJ.labFrameQuadrupole_XY*gqxy14+atomJ.labFrameQuadrupole_XZ*gqxz14+atomJ.labFrameQuadrupole_YZ*gqyz14))
              + atomI.labFrameQuadrupole_ZZ*(atomJ.labFrameQuadrupole_XX*gqxx16+atomJ.labFrameQuadrupole_YY*gqyy16+atomJ.labFrameQuadrupole_ZZ*gqzz16
              +2.0f*(atomJ.labFrameQuadrupole_XY*gqxy16+atomJ.labFrameQuadrupole_XZ*gqxz16+atomJ.labFrameQuadrupole_YZ*gqyz16))
       + 2.0f*(atomI.labFrameQuadrupole_XY*(atomJ.labFrameQuadrupole_XX*gqxx12+atomJ.labFrameQuadrupole_YY*gqyy12+atomJ.labFrameQuadrupole_ZZ*gqzz12
              +2.0f*(atomJ.labFrameQuadrupole_XY*gqxy12+atomJ.labFrameQuadrupole_XZ*gqxz12+atomJ.labFrameQuadrupole_YZ*gqyz12))
              + atomI.labFrameQuadrupole_XZ*(atomJ.labFrameQuadrupole_XX*gqxx13+atomJ.labFrameQuadrupole_YY*gqyy13+atomJ.labFrameQuadrupole_ZZ*gqzz13
              +2.0f*(atomJ.labFrameQuadrupole_XY*gqxy13+atomJ.labFrameQuadrupole_XZ*gqxz13+atomJ.labFrameQuadrupole_YZ*gqyz13))
              + atomI.labFrameQuadrupole_YZ*(atomJ.labFrameQuadrupole_XX*gqxx15+atomJ.labFrameQuadrupole_YY*gqyy15+atomJ.labFrameQuadrupole_ZZ*gqzz15
              +2.0f*(atomJ.labFrameQuadrupole_XY*gqxy15+atomJ.labFrameQuadrupole_XZ*gqxz15+atomJ.labFrameQuadrupole_YZ*gqyz15)));

    dedx = desymdx + 0.5f*(dewidx + dewkdx);

    desymdy = atomI.q * atomJ.q * gc3
                           - (atomI.labFrameDipole[0]*(atomJ.labFrameDipole[0]*gux6+atomJ.labFrameDipole[1]*guy6+atomJ.labFrameDipole[2]*guz6)
                             +atomI.labFrameDipole[1]*(atomJ.labFrameDipole[0]*gux8+atomJ.labFrameDipole[1]*guy8+atomJ.labFrameDipole[2]*guz8)
                             +atomI.labFrameDipole[2]*(atomJ.labFrameDipole[0]*gux9+atomJ.labFrameDipole[1]*guy9+atomJ.labFrameDipole[2]*guz9));

    dewidy = atomI.q*(atomJ.labFrameDipole[0]*gc6+atomJ.labFrameDipole[1]*gc8+atomJ.labFrameDipole[2]*gc9)
                      -atomJ.q*(atomI.labFrameDipole[0]*gux3+atomI.labFrameDipole[1]*guy3+atomI.labFrameDipole[2]*guz3)
                 +atomI.q*(atomJ.labFrameQuadrupole_XX*gc12+atomJ.labFrameQuadrupole_YY*gc17+atomJ.labFrameQuadrupole_ZZ*gc19
              +2.0f*(atomJ.labFrameQuadrupole_XY*gc14+atomJ.labFrameQuadrupole_XZ*gc15+atomJ.labFrameQuadrupole_YZ*gc18))
                 +atomJ.q*(atomI.labFrameQuadrupole_XX*gqxx3+atomI.labFrameQuadrupole_YY*gqyy3+atomI.labFrameQuadrupole_ZZ*gqzz3
              +2.0f*(atomI.labFrameQuadrupole_XY*gqxy3+atomI.labFrameQuadrupole_XZ*gqxz3+atomI.labFrameQuadrupole_YZ*gqyz3))
               - atomI.labFrameDipole[0]*(atomJ.labFrameQuadrupole_XX*gux12+atomJ.labFrameQuadrupole_YY*gux17+atomJ.labFrameQuadrupole_ZZ*gux19
              +2.0f*(atomJ.labFrameQuadrupole_XY*gux14+atomJ.labFrameQuadrupole_XZ*gux15+atomJ.labFrameQuadrupole_YZ*gux18))
               - atomI.labFrameDipole[1]*(atomJ.labFrameQuadrupole_XX*guy12+atomJ.labFrameQuadrupole_YY*guy17+atomJ.labFrameQuadrupole_ZZ*guy19
              +2.0f*(atomJ.labFrameQuadrupole_XY*guy14+atomJ.labFrameQuadrupole_XZ*guy15+atomJ.labFrameQuadrupole_YZ*guy18))
               - atomI.labFrameDipole[2]*(atomJ.labFrameQuadrupole_XX*guz12+atomJ.labFrameQuadrupole_YY*guz17+atomJ.labFrameQuadrupole_ZZ*guz19
              +2.0f*(atomJ.labFrameQuadrupole_XY*guz14+atomJ.labFrameQuadrupole_XZ*guz15+atomJ.labFrameQuadrupole_YZ*guz18))
               + atomJ.labFrameDipole[0]*(atomI.labFrameQuadrupole_XX*gqxx6+atomI.labFrameQuadrupole_YY*gqyy6+atomI.labFrameQuadrupole_ZZ*gqzz6
              +2.0f*(atomI.labFrameQuadrupole_XY*gqxy6+atomI.labFrameQuadrupole_XZ*gqxz6+atomI.labFrameQuadrupole_YZ*gqyz6))
               + atomJ.labFrameDipole[1]*(atomI.labFrameQuadrupole_XX*gqxx8+atomI.labFrameQuadrupole_YY*gqyy8+atomI.labFrameQuadrupole_ZZ*gqzz8
              +2.0f*(atomI.labFrameQuadrupole_XY*gqxy8+atomI.labFrameQuadrupole_XZ*gqxz8+atomI.labFrameQuadrupole_YZ*gqyz8))
               + atomJ.labFrameDipole[2]*(atomI.labFrameQuadrupole_XX*gqxx9+atomI.labFrameQuadrupole_YY*gqyy9+atomI.labFrameQuadrupole_ZZ*gqzz9
              +2.0f*(atomI.labFrameQuadrupole_XY*gqxy9+atomI.labFrameQuadrupole_XZ*gqxz9+atomI.labFrameQuadrupole_YZ*gqyz9))
              + atomI.labFrameQuadrupole_XX*(atomJ.labFrameQuadrupole_XX*gqxx12+atomJ.labFrameQuadrupole_YY*gqxx17+atomJ.labFrameQuadrupole_ZZ*gqxx19
              +2.0f*(atomJ.labFrameQuadrupole_XY*gqxx14+atomJ.labFrameQuadrupole_XZ*gqxx15+atomJ.labFrameQuadrupole_YZ*gqxx18))
              + atomI.labFrameQuadrupole_YY*(atomJ.labFrameQuadrupole_XX*gqyy12+atomJ.labFrameQuadrupole_YY*gqyy17+atomJ.labFrameQuadrupole_ZZ*gqyy19
              +2.0f*(atomJ.labFrameQuadrupole_XY*gqyy14+atomJ.labFrameQuadrupole_XZ*gqyy15+atomJ.labFrameQuadrupole_YZ*gqyy18))
              + atomI.labFrameQuadrupole_ZZ*(atomJ.labFrameQuadrupole_XX*gqzz12+atomJ.labFrameQuadrupole_YY*gqzz17+atomJ.labFrameQuadrupole_ZZ*gqzz19
              +2.0f*(atomJ.labFrameQuadrupole_XY*gqzz14+atomJ.labFrameQuadrupole_XZ*gqzz15+atomJ.labFrameQuadrupole_YZ*gqzz18))
       + 2.0f*(atomI.labFrameQuadrupole_XY*(atomJ.labFrameQuadrupole_XX*gqxy12+atomJ.labFrameQuadrupole_YY*gqxy17+atomJ.labFrameQuadrupole_ZZ*gqxy19
              +2.0f*(atomJ.labFrameQuadrupole_XY*gqxy14+atomJ.labFrameQuadrupole_XZ*gqxy15+atomJ.labFrameQuadrupole_YZ*gqxy18))
              + atomI.labFrameQuadrupole_XZ*(atomJ.labFrameQuadrupole_XX*gqxz12+atomJ.labFrameQuadrupole_YY*gqxz17+atomJ.labFrameQuadrupole_ZZ*gqxz19
              +2.0f*(atomJ.labFrameQuadrupole_XY*gqxz14+atomJ.labFrameQuadrupole_XZ*gqxz15+atomJ.labFrameQuadrupole_YZ*gqxz18))
              + atomI.labFrameQuadrupole_YZ*(atomJ.labFrameQuadrupole_XX*gqyz12+atomJ.labFrameQuadrupole_YY*gqyz17+atomJ.labFrameQuadrupole_ZZ*gqyz19
              +2.0f*(atomJ.labFrameQuadrupole_XY*gqyz14+atomJ.labFrameQuadrupole_XZ*gqyz15+atomJ.labFrameQuadrupole_YZ*gqyz18)));

    dewkdy = atomI.q*(atomJ.labFrameDipole[0]*gux3+atomJ.labFrameDipole[1]*guy3+atomJ.labFrameDipole[2]*guz3)
                      -atomJ.q*(atomI.labFrameDipole[0]*gc6+atomI.labFrameDipole[1]*gc8+atomI.labFrameDipole[2]*gc9)
                 +atomI.q*(atomJ.labFrameQuadrupole_XX*gqxx3+atomJ.labFrameQuadrupole_YY*gqyy3+atomJ.labFrameQuadrupole_ZZ*gqzz3
              +2.0f*(atomJ.labFrameQuadrupole_XY*gqxy3+atomJ.labFrameQuadrupole_XZ*gqxz3+atomJ.labFrameQuadrupole_YZ*gqyz3))
                 +atomJ.q*(atomI.labFrameQuadrupole_XX*gc12+atomI.labFrameQuadrupole_YY*gc17+atomI.labFrameQuadrupole_ZZ*gc19
              +2.0f*(atomI.labFrameQuadrupole_XY*gc14+atomI.labFrameQuadrupole_XZ*gc15+atomI.labFrameQuadrupole_YZ*gc18))
               - atomI.labFrameDipole[0]*(atomJ.labFrameQuadrupole_XX*gqxx6+atomJ.labFrameQuadrupole_YY*gqyy6+atomJ.labFrameQuadrupole_ZZ*gqzz6
              +2.0f*(atomJ.labFrameQuadrupole_XY*gqxy6+atomJ.labFrameQuadrupole_XZ*gqxz6+atomJ.labFrameQuadrupole_YZ*gqyz6))
               - atomI.labFrameDipole[1]*(atomJ.labFrameQuadrupole_XX*gqxx8+atomJ.labFrameQuadrupole_YY*gqyy8+atomJ.labFrameQuadrupole_ZZ*gqzz8
              +2.0f*(atomJ.labFrameQuadrupole_XY*gqxy8+atomJ.labFrameQuadrupole_XZ*gqxz8+atomJ.labFrameQuadrupole_YZ*gqyz8))
               - atomI.labFrameDipole[2]*(atomJ.labFrameQuadrupole_XX*gqxx9+atomJ.labFrameQuadrupole_YY*gqyy9+atomJ.labFrameQuadrupole_ZZ*gqzz9
              +2.0f*(atomJ.labFrameQuadrupole_XY*gqxy9+atomJ.labFrameQuadrupole_XZ*gqxz9+atomJ.labFrameQuadrupole_YZ*gqyz9))
               + atomJ.labFrameDipole[0]*(atomI.labFrameQuadrupole_XX*gux12+atomI.labFrameQuadrupole_YY*gux17+atomI.labFrameQuadrupole_ZZ*gux19
              +2.0f*(atomI.labFrameQuadrupole_XY*gux14+atomI.labFrameQuadrupole_XZ*gux15+atomI.labFrameQuadrupole_YZ*gux18))
               + atomJ.labFrameDipole[1]*(atomI.labFrameQuadrupole_XX*guy12+atomI.labFrameQuadrupole_YY*guy17+atomI.labFrameQuadrupole_ZZ*guy19
              +2.0f*(atomI.labFrameQuadrupole_XY*guy14+atomI.labFrameQuadrupole_XZ*guy15+atomI.labFrameQuadrupole_YZ*guy18))
               + atomJ.labFrameDipole[2]*(atomI.labFrameQuadrupole_XX*guz12+atomI.labFrameQuadrupole_YY*guz17+atomI.labFrameQuadrupole_ZZ*guz19
              +2.0f*(atomI.labFrameQuadrupole_XY*guz14+atomI.labFrameQuadrupole_XZ*guz15+atomI.labFrameQuadrupole_YZ*guz18))
              + atomI.labFrameQuadrupole_XX*(atomJ.labFrameQuadrupole_XX*gqxx12+atomJ.labFrameQuadrupole_YY*gqyy12+atomJ.labFrameQuadrupole_ZZ*gqzz12
              +2.0f*(atomJ.labFrameQuadrupole_XY*gqxy12+atomJ.labFrameQuadrupole_XZ*gqxz12+atomJ.labFrameQuadrupole_YZ*gqyz12))
              + atomI.labFrameQuadrupole_YY*(atomJ.labFrameQuadrupole_XX*gqxx17+atomJ.labFrameQuadrupole_YY*gqyy17+atomJ.labFrameQuadrupole_ZZ*gqzz17
              +2.0f*(atomJ.labFrameQuadrupole_XY*gqxy17+atomJ.labFrameQuadrupole_XZ*gqxz17+atomJ.labFrameQuadrupole_YZ*gqyz17))
              + atomI.labFrameQuadrupole_ZZ*(atomJ.labFrameQuadrupole_XX*gqxx19+atomJ.labFrameQuadrupole_YY*gqyy19+atomJ.labFrameQuadrupole_ZZ*gqzz19
              +2.0f*(atomJ.labFrameQuadrupole_XY*gqxy19+atomJ.labFrameQuadrupole_XZ*gqxz19+atomJ.labFrameQuadrupole_YZ*gqyz19))
       + 2.0f*(atomI.labFrameQuadrupole_XY*(atomJ.labFrameQuadrupole_XX*gqxx14+atomJ.labFrameQuadrupole_YY*gqyy14+atomJ.labFrameQuadrupole_ZZ*gqzz14
              +2.0f*(atomJ.labFrameQuadrupole_XY*gqxy14+atomJ.labFrameQuadrupole_XZ*gqxz14+atomJ.labFrameQuadrupole_YZ*gqyz14))
              + atomI.labFrameQuadrupole_XZ*(atomJ.labFrameQuadrupole_XX*gqxx15+atomJ.labFrameQuadrupole_YY*gqyy15+atomJ.labFrameQuadrupole_ZZ*gqzz15
              +2.0f*(atomJ.labFrameQuadrupole_XY*gqxy15+atomJ.labFrameQuadrupole_XZ*gqxz15+atomJ.labFrameQuadrupole_YZ*gqyz15))
              + atomI.labFrameQuadrupole_YZ*(atomJ.labFrameQuadrupole_XX*gqxx18+atomJ.labFrameQuadrupole_YY*gqyy18+atomJ.labFrameQuadrupole_ZZ*gqzz18
              +2.0f*(atomJ.labFrameQuadrupole_XY*gqxy18+atomJ.labFrameQuadrupole_XZ*gqxz18+atomJ.labFrameQuadrupole_YZ*gqyz18)));

    dedy = desymdy + 0.5f*(dewidy + dewkdy);

    desymdz = atomI.q * atomJ.q * gc4
                           - (atomI.labFrameDipole[0]*(atomJ.labFrameDipole[0]*gux7+atomJ.labFrameDipole[1]*guy7+atomJ.labFrameDipole[2]*guz7)
                             +atomI.labFrameDipole[1]*(atomJ.labFrameDipole[0]*gux9+atomJ.labFrameDipole[1]*guy9+atomJ.labFrameDipole[2]*guz9)
                             +atomI.labFrameDipole[2]*(atomJ.labFrameDipole[0]*gux10+atomJ.labFrameDipole[1]*guy10+atomJ.labFrameDipole[2]*guz10));

    dewidz = atomI.q*(atomJ.labFrameDipole[0]*gc7+atomJ.labFrameDipole[1]*gc9+atomJ.labFrameDipole[2]*gc10)
                      -atomJ.q*(atomI.labFrameDipole[0]*gux4+atomI.labFrameDipole[1]*guy4+atomI.labFrameDipole[2]*guz4)
                 +atomI.q*(atomJ.labFrameQuadrupole_XX*gc13+atomJ.labFrameQuadrupole_YY*gc18+atomJ.labFrameQuadrupole_ZZ*gc20
              +2.0f*(atomJ.labFrameQuadrupole_XY*gc15+atomJ.labFrameQuadrupole_XZ*gc16+atomJ.labFrameQuadrupole_YZ*gc19))
                 +atomJ.q*(atomI.labFrameQuadrupole_XX*gqxx4+atomI.labFrameQuadrupole_YY*gqyy4+atomI.labFrameQuadrupole_ZZ*gqzz4
              +2.0f*(atomI.labFrameQuadrupole_XY*gqxy4+atomI.labFrameQuadrupole_XZ*gqxz4+atomI.labFrameQuadrupole_YZ*gqyz4))
               - atomI.labFrameDipole[0]*(atomJ.labFrameQuadrupole_XX*gux13+atomJ.labFrameQuadrupole_YY*gux18+atomJ.labFrameQuadrupole_ZZ*gux20
              +2.0f*(atomJ.labFrameQuadrupole_XY*gux15+atomJ.labFrameQuadrupole_XZ*gux16+atomJ.labFrameQuadrupole_YZ*gux19))
               - atomI.labFrameDipole[1]*(atomJ.labFrameQuadrupole_XX*guy13+atomJ.labFrameQuadrupole_YY*guy18+atomJ.labFrameQuadrupole_ZZ*guy20
              +2.0f*(atomJ.labFrameQuadrupole_XY*guy15+atomJ.labFrameQuadrupole_XZ*guy16+atomJ.labFrameQuadrupole_YZ*guy19))
               - atomI.labFrameDipole[2]*(atomJ.labFrameQuadrupole_XX*guz13+atomJ.labFrameQuadrupole_YY*guz18+atomJ.labFrameQuadrupole_ZZ*guz20
              +2.0f*(atomJ.labFrameQuadrupole_XY*guz15+atomJ.labFrameQuadrupole_XZ*guz16+atomJ.labFrameQuadrupole_YZ*guz19))
               + atomJ.labFrameDipole[0]*(atomI.labFrameQuadrupole_XX*gqxx7+atomI.labFrameQuadrupole_YY*gqyy7+atomI.labFrameQuadrupole_ZZ*gqzz7
              +2.0f*(atomI.labFrameQuadrupole_XY*gqxy7+atomI.labFrameQuadrupole_XZ*gqxz7+atomI.labFrameQuadrupole_YZ*gqyz7))
               + atomJ.labFrameDipole[1]*(atomI.labFrameQuadrupole_XX*gqxx9+atomI.labFrameQuadrupole_YY*gqyy9+atomI.labFrameQuadrupole_ZZ*gqzz9
              +2.0f*(atomI.labFrameQuadrupole_XY*gqxy9+atomI.labFrameQuadrupole_XZ*gqxz9+atomI.labFrameQuadrupole_YZ*gqyz9))
               + atomJ.labFrameDipole[2]*(atomI.labFrameQuadrupole_XX*gqxx10+atomI.labFrameQuadrupole_YY*gqyy10+atomI.labFrameQuadrupole_ZZ*gqzz10
              +2.0f*(atomI.labFrameQuadrupole_XY*gqxy10+atomI.labFrameQuadrupole_XZ*gqxz10+atomI.labFrameQuadrupole_YZ*gqyz10))
              + atomI.labFrameQuadrupole_XX*(atomJ.labFrameQuadrupole_XX*gqxx13+atomJ.labFrameQuadrupole_YY*gqxx18+atomJ.labFrameQuadrupole_ZZ*gqxx20
              +2.0f*(atomJ.labFrameQuadrupole_XY*gqxx15+atomJ.labFrameQuadrupole_XZ*gqxx16+atomJ.labFrameQuadrupole_YZ*gqxx19))
              + atomI.labFrameQuadrupole_YY*(atomJ.labFrameQuadrupole_XX*gqyy13+atomJ.labFrameQuadrupole_YY*gqyy18+atomJ.labFrameQuadrupole_ZZ*gqyy20
              +2.0f*(atomJ.labFrameQuadrupole_XY*gqyy15+atomJ.labFrameQuadrupole_XZ*gqyy16+atomJ.labFrameQuadrupole_YZ*gqyy19))
              + atomI.labFrameQuadrupole_ZZ*(atomJ.labFrameQuadrupole_XX*gqzz13+atomJ.labFrameQuadrupole_YY*gqzz18+atomJ.labFrameQuadrupole_ZZ*gqzz20
              +2.0f*(atomJ.labFrameQuadrupole_XY*gqzz15+atomJ.labFrameQuadrupole_XZ*gqzz16+atomJ.labFrameQuadrupole_YZ*gqzz19))
       + 2.0f*(atomI.labFrameQuadrupole_XY*(atomJ.labFrameQuadrupole_XX*gqxy13+atomJ.labFrameQuadrupole_YY*gqxy18+atomJ.labFrameQuadrupole_ZZ*gqxy20
              +2.0f*(atomJ.labFrameQuadrupole_XY*gqxy15+atomJ.labFrameQuadrupole_XZ*gqxy16+atomJ.labFrameQuadrupole_YZ*gqxy19))
              + atomI.labFrameQuadrupole_XZ*(atomJ.labFrameQuadrupole_XX*gqxz13+atomJ.labFrameQuadrupole_YY*gqxz18+atomJ.labFrameQuadrupole_ZZ*gqxz20
              +2.0f*(atomJ.labFrameQuadrupole_XY*gqxz15+atomJ.labFrameQuadrupole_XZ*gqxz16+atomJ.labFrameQuadrupole_YZ*gqxz19))
              + atomI.labFrameQuadrupole_YZ*(atomJ.labFrameQuadrupole_XX*gqyz13+atomJ.labFrameQuadrupole_YY*gqyz18+atomJ.labFrameQuadrupole_ZZ*gqyz20
              +2.0f*(atomJ.labFrameQuadrupole_XY*gqyz15+atomJ.labFrameQuadrupole_XZ*gqyz16+atomJ.labFrameQuadrupole_YZ*gqyz19)));

    dewkdz = atomI.q*(atomJ.labFrameDipole[0]*gux4+atomJ.labFrameDipole[1]*guy4+atomJ.labFrameDipole[2]*guz4)
                      -atomJ.q*(atomI.labFrameDipole[0]*gc7+atomI.labFrameDipole[1]*gc9+atomI.labFrameDipole[2]*gc10)
                 +atomI.q*(atomJ.labFrameQuadrupole_XX*gqxx4+atomJ.labFrameQuadrupole_YY*gqyy4+atomJ.labFrameQuadrupole_ZZ*gqzz4
              +2.0f*(atomJ.labFrameQuadrupole_XY*gqxy4+atomJ.labFrameQuadrupole_XZ*gqxz4+atomJ.labFrameQuadrupole_YZ*gqyz4))
                 +atomJ.q*(atomI.labFrameQuadrupole_XX*gc13+atomI.labFrameQuadrupole_YY*gc18+atomI.labFrameQuadrupole_ZZ*gc20
              +2.0f*(atomI.labFrameQuadrupole_XY*gc15+atomI.labFrameQuadrupole_XZ*gc16+atomI.labFrameQuadrupole_YZ*gc19))
               - atomI.labFrameDipole[0]*(atomJ.labFrameQuadrupole_XX*gqxx7+atomJ.labFrameQuadrupole_YY*gqyy7+atomJ.labFrameQuadrupole_ZZ*gqzz7
              +2.0f*(atomJ.labFrameQuadrupole_XY*gqxy7+atomJ.labFrameQuadrupole_XZ*gqxz7+atomJ.labFrameQuadrupole_YZ*gqyz7))
               - atomI.labFrameDipole[1]*(atomJ.labFrameQuadrupole_XX*gqxx9+atomJ.labFrameQuadrupole_YY*gqyy9+atomJ.labFrameQuadrupole_ZZ*gqzz9
              +2.0f*(atomJ.labFrameQuadrupole_XY*gqxy9+atomJ.labFrameQuadrupole_XZ*gqxz9+atomJ.labFrameQuadrupole_YZ*gqyz9))
               - atomI.labFrameDipole[2]*(atomJ.labFrameQuadrupole_XX*gqxx10+atomJ.labFrameQuadrupole_YY*gqyy10+atomJ.labFrameQuadrupole_ZZ*gqzz10
              +2.0f*(atomJ.labFrameQuadrupole_XY*gqxy10+atomJ.labFrameQuadrupole_XZ*gqxz10+atomJ.labFrameQuadrupole_YZ*gqyz10))
               + atomJ.labFrameDipole[0]*(atomI.labFrameQuadrupole_XX*gux13+atomI.labFrameQuadrupole_YY*gux18+atomI.labFrameQuadrupole_ZZ*gux20
              +2.0f*(atomI.labFrameQuadrupole_XY*gux15+atomI.labFrameQuadrupole_XZ*gux16+atomI.labFrameQuadrupole_YZ*gux19))
               + atomJ.labFrameDipole[1]*(atomI.labFrameQuadrupole_XX*guy13+atomI.labFrameQuadrupole_YY*guy18+atomI.labFrameQuadrupole_ZZ*guy20
              +2.0f*(atomI.labFrameQuadrupole_XY*guy15+atomI.labFrameQuadrupole_XZ*guy16+atomI.labFrameQuadrupole_YZ*guy19))
               + atomJ.labFrameDipole[2]*(atomI.labFrameQuadrupole_XX*guz13+atomI.labFrameQuadrupole_YY*guz18+atomI.labFrameQuadrupole_ZZ*guz20
              +2.0f*(atomI.labFrameQuadrupole_XY*guz15+atomI.labFrameQuadrupole_XZ*guz16+atomI.labFrameQuadrupole_YZ*guz19))
              + atomI.labFrameQuadrupole_XX*(atomJ.labFrameQuadrupole_XX*gqxx13+atomJ.labFrameQuadrupole_YY*gqyy13+atomJ.labFrameQuadrupole_ZZ*gqzz13
              +2.0f*(atomJ.labFrameQuadrupole_XY*gqxy13+atomJ.labFrameQuadrupole_XZ*gqxz13+atomJ.labFrameQuadrupole_YZ*gqyz13))
              + atomI.labFrameQuadrupole_YY*(atomJ.labFrameQuadrupole_XX*gqxx18+atomJ.labFrameQuadrupole_YY*gqyy18+atomJ.labFrameQuadrupole_ZZ*gqzz18
              +2.0f*(atomJ.labFrameQuadrupole_XY*gqxy18+atomJ.labFrameQuadrupole_XZ*gqxz18+atomJ.labFrameQuadrupole_YZ*gqyz18))
              + atomI.labFrameQuadrupole_ZZ*(atomJ.labFrameQuadrupole_XX*gqxx20+atomJ.labFrameQuadrupole_YY*gqyy20+atomJ.labFrameQuadrupole_ZZ*gqzz20
              +2.0f*(atomJ.labFrameQuadrupole_XY*gqxy20+atomJ.labFrameQuadrupole_XZ*gqxz20+atomJ.labFrameQuadrupole_YZ*gqyz20))
       + 2.0f*(atomI.labFrameQuadrupole_XY*(atomJ.labFrameQuadrupole_XX*gqxx15+atomJ.labFrameQuadrupole_YY*gqyy15+atomJ.labFrameQuadrupole_ZZ*gqzz15
              +2.0f*(atomJ.labFrameQuadrupole_XY*gqxy15+atomJ.labFrameQuadrupole_XZ*gqxz15+atomJ.labFrameQuadrupole_YZ*gqyz15))
              + atomI.labFrameQuadrupole_XZ*(atomJ.labFrameQuadrupole_XX*gqxx16+atomJ.labFrameQuadrupole_YY*gqyy16+atomJ.labFrameQuadrupole_ZZ*gqzz16
              +2.0f*(atomJ.labFrameQuadrupole_XY*gqxy16+atomJ.labFrameQuadrupole_XZ*gqxz16+atomJ.labFrameQuadrupole_YZ*gqyz16))
              + atomI.labFrameQuadrupole_YZ*(atomJ.labFrameQuadrupole_XX*gqxx19+atomJ.labFrameQuadrupole_YY*gqyy19+atomJ.labFrameQuadrupole_ZZ*gqzz19
              +2.0f*(atomJ.labFrameQuadrupole_XY*gqxy19+atomJ.labFrameQuadrupole_XZ*gqxz19+atomJ.labFrameQuadrupole_YZ*gqyz19)));

    dedz = desymdz + 0.5f*(dewidz + dewkdz);

    desymdr = atomI.q * atomJ.q * gc21
                           - (atomI.labFrameDipole[0]*(atomJ.labFrameDipole[0]*gux22+atomJ.labFrameDipole[1]*guy22+atomJ.labFrameDipole[2]*guz22)
                             +atomI.labFrameDipole[1]*(atomJ.labFrameDipole[0]*gux23+atomJ.labFrameDipole[1]*guy23+atomJ.labFrameDipole[2]*guz23)
                             +atomI.labFrameDipole[2]*(atomJ.labFrameDipole[0]*gux24+atomJ.labFrameDipole[1]*guy24+atomJ.labFrameDipole[2]*guz24));

    dewidr = atomI.q*(atomJ.labFrameDipole[0]*gc22+atomJ.labFrameDipole[1]*gc23+atomJ.labFrameDipole[2]*gc24)
                      -atomJ.q*(atomI.labFrameDipole[0]*gux21+atomI.labFrameDipole[1]*guy21+atomI.labFrameDipole[2]*guz21)
                 +atomI.q*(atomJ.labFrameQuadrupole_XX*gc25+atomJ.labFrameQuadrupole_YY*gc28+atomJ.labFrameQuadrupole_ZZ*gc30
              +2.0f*(atomJ.labFrameQuadrupole_XY*gc26+atomJ.labFrameQuadrupole_XZ*gc27+atomJ.labFrameQuadrupole_YZ*gc29))
                 +atomJ.q*(atomI.labFrameQuadrupole_XX*gqxx21+atomI.labFrameQuadrupole_YY*gqyy21+atomI.labFrameQuadrupole_ZZ*gqzz21
              +2.0f*(atomI.labFrameQuadrupole_XY*gqxy21+atomI.labFrameQuadrupole_XZ*gqxz21+atomI.labFrameQuadrupole_YZ*gqyz21))
               - atomI.labFrameDipole[0]*(atomJ.labFrameQuadrupole_XX*gux25+atomJ.labFrameQuadrupole_YY*gux28+atomJ.labFrameQuadrupole_ZZ*gux30
              +2.0f*(atomJ.labFrameQuadrupole_XY*gux26+atomJ.labFrameQuadrupole_XZ*gux27+atomJ.labFrameQuadrupole_YZ*gux29))
               - atomI.labFrameDipole[1]*(atomJ.labFrameQuadrupole_XX*guy25+atomJ.labFrameQuadrupole_YY*guy28+atomJ.labFrameQuadrupole_ZZ*guy30
              +2.0f*(atomJ.labFrameQuadrupole_XY*guy26+atomJ.labFrameQuadrupole_XZ*guy27+atomJ.labFrameQuadrupole_YZ*guy29))
               - atomI.labFrameDipole[2]*(atomJ.labFrameQuadrupole_XX*guz25+atomJ.labFrameQuadrupole_YY*guz28+atomJ.labFrameQuadrupole_ZZ*guz30
              +2.0f*(atomJ.labFrameQuadrupole_XY*guz26+atomJ.labFrameQuadrupole_XZ*guz27+atomJ.labFrameQuadrupole_YZ*guz29))
               + atomJ.labFrameDipole[0]*(atomI.labFrameQuadrupole_XX*gqxx22+atomI.labFrameQuadrupole_YY*gqyy22+atomI.labFrameQuadrupole_ZZ*gqzz22
              +2.0f*(atomI.labFrameQuadrupole_XY*gqxy22+atomI.labFrameQuadrupole_XZ*gqxz22+atomI.labFrameQuadrupole_YZ*gqyz22))
               + atomJ.labFrameDipole[1]*(atomI.labFrameQuadrupole_XX*gqxx23+atomI.labFrameQuadrupole_YY*gqyy23+atomI.labFrameQuadrupole_ZZ*gqzz23
              +2.0f*(atomI.labFrameQuadrupole_XY*gqxy23+atomI.labFrameQuadrupole_XZ*gqxz23+atomI.labFrameQuadrupole_YZ*gqyz23))
               + atomJ.labFrameDipole[2]*(atomI.labFrameQuadrupole_XX*gqxx24+atomI.labFrameQuadrupole_YY*gqyy24+atomI.labFrameQuadrupole_ZZ*gqzz24
              +2.0f*(atomI.labFrameQuadrupole_XY*gqxy24+atomI.labFrameQuadrupole_XZ*gqxz24+atomI.labFrameQuadrupole_YZ*gqyz24))
              + atomI.labFrameQuadrupole_XX*(atomJ.labFrameQuadrupole_XX*gqxx25+atomJ.labFrameQuadrupole_YY*gqxx28+atomJ.labFrameQuadrupole_ZZ*gqxx30
              +2.0f*(atomJ.labFrameQuadrupole_XY*gqxx26+atomJ.labFrameQuadrupole_XZ*gqxx27+atomJ.labFrameQuadrupole_YZ*gqxx29))
              + atomI.labFrameQuadrupole_YY*(atomJ.labFrameQuadrupole_XX*gqyy25+atomJ.labFrameQuadrupole_YY*gqyy28+atomJ.labFrameQuadrupole_ZZ*gqyy30
              +2.0f*(atomJ.labFrameQuadrupole_XY*gqyy26+atomJ.labFrameQuadrupole_XZ*gqyy27+atomJ.labFrameQuadrupole_YZ*gqyy29))
              + atomI.labFrameQuadrupole_ZZ*(atomJ.labFrameQuadrupole_XX*gqzz25+atomJ.labFrameQuadrupole_YY*gqzz28+atomJ.labFrameQuadrupole_ZZ*gqzz30
              +2.0f*(atomJ.labFrameQuadrupole_XY*gqzz26+atomJ.labFrameQuadrupole_XZ*gqzz27+atomJ.labFrameQuadrupole_YZ*gqzz29))
              + 2.0f*(atomI.labFrameQuadrupole_XY*(atomJ.labFrameQuadrupole_XX*gqxy25+atomJ.labFrameQuadrupole_YY*gqxy28+atomJ.labFrameQuadrupole_ZZ*gqxy30
              +2.0f*(atomJ.labFrameQuadrupole_XY*gqxy26+atomJ.labFrameQuadrupole_XZ*gqxy27+atomJ.labFrameQuadrupole_YZ*gqxy29))
              + atomI.labFrameQuadrupole_XZ*(atomJ.labFrameQuadrupole_XX*gqxz25+atomJ.labFrameQuadrupole_YY*gqxz28+atomJ.labFrameQuadrupole_ZZ*gqxz30
              +2.0f*(atomJ.labFrameQuadrupole_XY*gqxz26+atomJ.labFrameQuadrupole_XZ*gqxz27+atomJ.labFrameQuadrupole_YZ*gqxz29))
              + atomI.labFrameQuadrupole_YZ*(atomJ.labFrameQuadrupole_XX*gqyz25+atomJ.labFrameQuadrupole_YY*gqyz28+atomJ.labFrameQuadrupole_ZZ*gqyz30
              +2.0f*(atomJ.labFrameQuadrupole_XY*gqyz26+atomJ.labFrameQuadrupole_XZ*gqyz27+atomJ.labFrameQuadrupole_YZ*gqyz29)));

    dewkdr = atomI.q*(atomJ.labFrameDipole[0]*gux21+atomJ.labFrameDipole[1]*guy21+atomJ.labFrameDipole[2]*guz21)
                      -atomJ.q*(atomI.labFrameDipole[0]*gc22+atomI.labFrameDipole[1]*gc23+atomI.labFrameDipole[2]*gc24)
                 +atomI.q*(atomJ.labFrameQuadrupole_XX*gqxx21+atomJ.labFrameQuadrupole_YY*gqyy21+atomJ.labFrameQuadrupole_ZZ*gqzz21
              +2.0f*(atomJ.labFrameQuadrupole_XY*gqxy21+atomJ.labFrameQuadrupole_XZ*gqxz21+atomJ.labFrameQuadrupole_YZ*gqyz21))
                 +atomJ.q*(atomI.labFrameQuadrupole_XX*gc25+atomI.labFrameQuadrupole_YY*gc28+atomI.labFrameQuadrupole_ZZ*gc30
              +2.0f*(atomI.labFrameQuadrupole_XY*gc26+atomI.labFrameQuadrupole_XZ*gc27+atomI.labFrameQuadrupole_YZ*gc29))
               - atomI.labFrameDipole[0]*(atomJ.labFrameQuadrupole_XX*gqxx22+atomJ.labFrameQuadrupole_YY*gqyy22+atomJ.labFrameQuadrupole_ZZ*gqzz22
              +2.0f*(atomJ.labFrameQuadrupole_XY*gqxy22+atomJ.labFrameQuadrupole_XZ*gqxz22+atomJ.labFrameQuadrupole_YZ*gqyz22))
               - atomI.labFrameDipole[1]*(atomJ.labFrameQuadrupole_XX*gqxx23+atomJ.labFrameQuadrupole_YY*gqyy23+atomJ.labFrameQuadrupole_ZZ*gqzz23
              +2.0f*(atomJ.labFrameQuadrupole_XY*gqxy23+atomJ.labFrameQuadrupole_XZ*gqxz23+atomJ.labFrameQuadrupole_YZ*gqyz23))
               - atomI.labFrameDipole[2]*(atomJ.labFrameQuadrupole_XX*gqxx24+atomJ.labFrameQuadrupole_YY*gqyy24+atomJ.labFrameQuadrupole_ZZ*gqzz24
              +2.0f*(atomJ.labFrameQuadrupole_XY*gqxy24+atomJ.labFrameQuadrupole_XZ*gqxz24+atomJ.labFrameQuadrupole_YZ*gqyz24))
               + atomJ.labFrameDipole[0]*(atomI.labFrameQuadrupole_XX*gux25+atomI.labFrameQuadrupole_YY*gux28+atomI.labFrameQuadrupole_ZZ*gux30
              +2.0f*(atomI.labFrameQuadrupole_XY*gux26+atomI.labFrameQuadrupole_XZ*gux27+atomI.labFrameQuadrupole_YZ*gux29))
               + atomJ.labFrameDipole[1]*(atomI.labFrameQuadrupole_XX*guy25+atomI.labFrameQuadrupole_YY*guy28+atomI.labFrameQuadrupole_ZZ*guy30
              +2.0f*(atomI.labFrameQuadrupole_XY*guy26+atomI.labFrameQuadrupole_XZ*guy27+atomI.labFrameQuadrupole_YZ*guy29))
               + atomJ.labFrameDipole[2]*(atomI.labFrameQuadrupole_XX*guz25+atomI.labFrameQuadrupole_YY*guz28+atomI.labFrameQuadrupole_ZZ*guz30
              +2.0f*(atomI.labFrameQuadrupole_XY*guz26+atomI.labFrameQuadrupole_XZ*guz27+atomI.labFrameQuadrupole_YZ*guz29))
              + atomI.labFrameQuadrupole_XX*(atomJ.labFrameQuadrupole_XX*gqxx25+atomJ.labFrameQuadrupole_YY*gqyy25+atomJ.labFrameQuadrupole_ZZ*gqzz25
              +2.0f*(atomJ.labFrameQuadrupole_XY*gqxy25+atomJ.labFrameQuadrupole_XZ*gqxz25+atomJ.labFrameQuadrupole_YZ*gqyz25))
              + atomI.labFrameQuadrupole_YY*(atomJ.labFrameQuadrupole_XX*gqxx28+atomJ.labFrameQuadrupole_YY*gqyy28+atomJ.labFrameQuadrupole_ZZ*gqzz28
              +2.0f*(atomJ.labFrameQuadrupole_XY*gqxy28+atomJ.labFrameQuadrupole_XZ*gqxz28+atomJ.labFrameQuadrupole_YZ*gqyz28))
              + atomI.labFrameQuadrupole_ZZ*(atomJ.labFrameQuadrupole_XX*gqxx30+atomJ.labFrameQuadrupole_YY*gqyy30+atomJ.labFrameQuadrupole_ZZ*gqzz30
              +2.0f*(atomJ.labFrameQuadrupole_XY*gqxy30+atomJ.labFrameQuadrupole_XZ*gqxz30+atomJ.labFrameQuadrupole_YZ*gqyz30))
              + 2.0f*(atomI.labFrameQuadrupole_XY*(atomJ.labFrameQuadrupole_XX*gqxx26+atomJ.labFrameQuadrupole_YY*gqyy26+atomJ.labFrameQuadrupole_ZZ*gqzz26
              +2.0f*(atomJ.labFrameQuadrupole_XY*gqxy26+atomJ.labFrameQuadrupole_XZ*gqxz26+atomJ.labFrameQuadrupole_YZ*gqyz26))
              + atomI.labFrameQuadrupole_XZ*(atomJ.labFrameQuadrupole_XX*gqxx27+atomJ.labFrameQuadrupole_YY*gqyy27+atomJ.labFrameQuadrupole_ZZ*gqzz27
              +2.0f*(atomJ.labFrameQuadrupole_XY*gqxy27+atomJ.labFrameQuadrupole_XZ*gqxz27+atomJ.labFrameQuadrupole_YZ*gqyz27))
              + atomI.labFrameQuadrupole_YZ*(atomJ.labFrameQuadrupole_XX*gqxx29+atomJ.labFrameQuadrupole_YY*gqyy29+atomJ.labFrameQuadrupole_ZZ*gqzz29
              +2.0f*(atomJ.labFrameQuadrupole_XY*gqxy29+atomJ.labFrameQuadrupole_XZ*gqxz29+atomJ.labFrameQuadrupole_YZ*gqyz29)));

    dsumdr = desymdr + 0.5f*(dewidr + dewkdr);
    drbi = atomJ.bornRadius*dsumdr;
    drbk = atomI.bornRadius*dsumdr;

    // torque on permanent dipoles due to permanent reaction field

    float trq1   = 0.0f;
    float trq2   = 0.0f;
    float trq3   = 0.0f;

    float trq_k1 = 0.0f;
    float trq_k2 = 0.0f;
    float trq_k3 = 0.0f;

    if ( sameAtom == 0 )
    {

        float fid1 = atomJ.labFrameDipole[0]*gux2 + atomJ.labFrameDipole[1]*gux3 + atomJ.labFrameDipole[2]*gux4
                + 0.5f*(atomJ.q*gux1+atomJ.labFrameQuadrupole_XX*gux5+atomJ.labFrameQuadrupole_YY*gux8+atomJ.labFrameQuadrupole_ZZ*gux10
                      +2.0f*(atomJ.labFrameQuadrupole_XY*gux6+atomJ.labFrameQuadrupole_XZ*gux7+atomJ.labFrameQuadrupole_YZ*gux9)
                      +atomJ.q*gc2+atomJ.labFrameQuadrupole_XX*gqxx2+atomJ.labFrameQuadrupole_YY*gqyy2+atomJ.labFrameQuadrupole_ZZ*gqzz2
                      +2.0f*(atomJ.labFrameQuadrupole_XY*gqxy2+atomJ.labFrameQuadrupole_XZ*gqxz2+atomJ.labFrameQuadrupole_YZ*gqyz2));

        float fid2 = atomJ.labFrameDipole[0]*guy2 + atomJ.labFrameDipole[1]*guy3 + atomJ.labFrameDipole[2]*guy4
                + 0.5f*(atomJ.q*guy1+atomJ.labFrameQuadrupole_XX*guy5+atomJ.labFrameQuadrupole_YY*guy8+atomJ.labFrameQuadrupole_ZZ*guy10
                      +2.0f*(atomJ.labFrameQuadrupole_XY*guy6+atomJ.labFrameQuadrupole_XZ*guy7+atomJ.labFrameQuadrupole_YZ*guy9)
                      +atomJ.q*gc3+atomJ.labFrameQuadrupole_XX*gqxx3+atomJ.labFrameQuadrupole_YY*gqyy3+atomJ.labFrameQuadrupole_ZZ*gqzz3
                      +2.0f*(atomJ.labFrameQuadrupole_XY*gqxy3+atomJ.labFrameQuadrupole_XZ*gqxz3+atomJ.labFrameQuadrupole_YZ*gqyz3));

        float fid3 = atomJ.labFrameDipole[0]*guz2 + atomJ.labFrameDipole[1]*guz3 + atomJ.labFrameDipole[2]*guz4
                + 0.5f*(atomJ.q*guz1+atomJ.labFrameQuadrupole_XX*guz5+atomJ.labFrameQuadrupole_YY*guz8+atomJ.labFrameQuadrupole_ZZ*guz10
                      +2.0f*(atomJ.labFrameQuadrupole_XY*guz6+atomJ.labFrameQuadrupole_XZ*guz7+atomJ.labFrameQuadrupole_YZ*guz9)
                      +atomJ.q*gc4+atomJ.labFrameQuadrupole_XX*gqxx4+atomJ.labFrameQuadrupole_YY*gqyy4+atomJ.labFrameQuadrupole_ZZ*gqzz4
                      +2.0f*(atomJ.labFrameQuadrupole_XY*gqxy4+atomJ.labFrameQuadrupole_XZ*gqxz4+atomJ.labFrameQuadrupole_YZ*gqyz4));

        float fkd1 = atomI.labFrameDipole[0]*gux2 + atomI.labFrameDipole[1]*gux3 + atomI.labFrameDipole[2]*gux4
                - 0.5f*(atomI.q*gux1+atomI.labFrameQuadrupole_XX*gux5+atomI.labFrameQuadrupole_YY*gux8+atomI.labFrameQuadrupole_ZZ*gux10
                      +2.0f*(atomI.labFrameQuadrupole_XY*gux6+atomI.labFrameQuadrupole_XZ*gux7+atomI.labFrameQuadrupole_YZ*gux9)
                      +atomI.q*gc2+atomI.labFrameQuadrupole_XX*gqxx2+atomI.labFrameQuadrupole_YY*gqyy2+atomI.labFrameQuadrupole_ZZ*gqzz2
                      +2.0f*(atomI.labFrameQuadrupole_XY*gqxy2+atomI.labFrameQuadrupole_XZ*gqxz2+atomI.labFrameQuadrupole_YZ*gqyz2));

        float fkd2 = atomI.labFrameDipole[0]*guy2 + atomI.labFrameDipole[1]*guy3 + atomI.labFrameDipole[2]*guy4
                - 0.5f*(atomI.q*guy1+atomI.labFrameQuadrupole_XX*guy5+atomI.labFrameQuadrupole_YY*guy8+atomI.labFrameQuadrupole_ZZ*guy10
                      +2.0f*(atomI.labFrameQuadrupole_XY*guy6+atomI.labFrameQuadrupole_XZ*guy7+atomI.labFrameQuadrupole_YZ*guy9)
                      +atomI.q*gc3+atomI.labFrameQuadrupole_XX*gqxx3+atomI.labFrameQuadrupole_YY*gqyy3+atomI.labFrameQuadrupole_ZZ*gqzz3
                      +2.0f*(atomI.labFrameQuadrupole_XY*gqxy3+atomI.labFrameQuadrupole_XZ*gqxz3+atomI.labFrameQuadrupole_YZ*gqyz3));

        float fkd3 = atomI.labFrameDipole[0]*guz2 + atomI.labFrameDipole[1]*guz3 + atomI.labFrameDipole[2]*guz4
                - 0.5f*(atomI.q*guz1+atomI.labFrameQuadrupole_XX*guz5+atomI.labFrameQuadrupole_YY*guz8+atomI.labFrameQuadrupole_ZZ*guz10
                      +2.0f*(atomI.labFrameQuadrupole_XY*guz6+atomI.labFrameQuadrupole_XZ*guz7+atomI.labFrameQuadrupole_YZ*guz9)
                      +atomI.q*gc4+atomI.labFrameQuadrupole_XX*gqxx4+atomI.labFrameQuadrupole_YY*gqyy4+atomI.labFrameQuadrupole_ZZ*gqzz4
                      +2.0f*(atomI.labFrameQuadrupole_XY*gqxy4+atomI.labFrameQuadrupole_XZ*gqxz4+atomI.labFrameQuadrupole_YZ*gqyz4));

        trq1    = atomI.labFrameDipole[1]*fid3 - atomI.labFrameDipole[2]*fid2;
        trq2    = atomI.labFrameDipole[2]*fid1 - atomI.labFrameDipole[0]*fid3;
        trq3    = atomI.labFrameDipole[0]*fid2 - atomI.labFrameDipole[1]*fid1;

        trq_k1  = atomJ.labFrameDipole[1]*fkd3 - atomJ.labFrameDipole[2]*fkd2;
        trq_k2  = atomJ.labFrameDipole[2]*fkd1 - atomJ.labFrameDipole[0]*fkd3;
        trq_k3  = atomJ.labFrameDipole[0]*fkd2 - atomJ.labFrameDipole[1]*fkd1;

        // torque on quadrupoles due to permanent reaction field gradient

        float fidg11 =
                - 0.5f*(atomJ.q*gqxx1+atomJ.labFrameDipole[0]*gqxx2+atomJ.labFrameDipole[1]*gqxx3+atomJ.labFrameDipole[2]*gqxx4
                      +atomJ.labFrameQuadrupole_XX*gqxx5+atomJ.labFrameQuadrupole_YY*gqxx8+atomJ.labFrameQuadrupole_ZZ*gqxx10
                      +2.0f*(atomJ.labFrameQuadrupole_XY*gqxx6+atomJ.labFrameQuadrupole_XZ*gqxx7+atomJ.labFrameQuadrupole_YZ*gqxx9)
                      +atomJ.q*gc5+atomJ.labFrameDipole[0]*gux5+atomJ.labFrameDipole[1]*guy5+atomJ.labFrameDipole[2]*guz5
                      +atomJ.labFrameQuadrupole_XX*gqxx5+atomJ.labFrameQuadrupole_YY*gqyy5+atomJ.labFrameQuadrupole_ZZ*gqzz5
                      +2.0f*(atomJ.labFrameQuadrupole_XY*gqxy5+atomJ.labFrameQuadrupole_XZ*gqxz5+atomJ.labFrameQuadrupole_YZ*gqyz5));

        float fidg12 =
                - 0.5f*(atomJ.q*gqxy1+atomJ.labFrameDipole[0]*gqxy2+atomJ.labFrameDipole[1]*gqxy3+atomJ.labFrameDipole[2]*gqxy4
                      +atomJ.labFrameQuadrupole_XX*gqxy5+atomJ.labFrameQuadrupole_YY*gqxy8+atomJ.labFrameQuadrupole_ZZ*gqxy10
                      +2.0f*(atomJ.labFrameQuadrupole_XY*gqxy6+atomJ.labFrameQuadrupole_XZ*gqxy7+atomJ.labFrameQuadrupole_YZ*gqxy9)
                      +atomJ.q*gc6+atomJ.labFrameDipole[0]*gux6+atomJ.labFrameDipole[1]*guy6+atomJ.labFrameDipole[2]*guz6
                      +atomJ.labFrameQuadrupole_XX*gqxx6+atomJ.labFrameQuadrupole_YY*gqyy6+atomJ.labFrameQuadrupole_ZZ*gqzz6
                      +2.0f*(atomJ.labFrameQuadrupole_XY*gqxy6+atomJ.labFrameQuadrupole_XZ*gqxz6+atomJ.labFrameQuadrupole_YZ*gqyz6));

        float fidg13 =
                - 0.5f*(atomJ.q*gqxz1+atomJ.labFrameDipole[0]*gqxz2+atomJ.labFrameDipole[1]*gqxz3+atomJ.labFrameDipole[2]*gqxz4
                      +atomJ.labFrameQuadrupole_XX*gqxz5+atomJ.labFrameQuadrupole_YY*gqxz8+atomJ.labFrameQuadrupole_ZZ*gqxz10
                      +2.0f*(atomJ.labFrameQuadrupole_XY*gqxz6+atomJ.labFrameQuadrupole_XZ*gqxz7+atomJ.labFrameQuadrupole_YZ*gqxz9)
                      +atomJ.q*gc7+atomJ.labFrameDipole[0]*gux7+atomJ.labFrameDipole[1]*guy7+atomJ.labFrameDipole[2]*guz7
                      +atomJ.labFrameQuadrupole_XX*gqxx7+atomJ.labFrameQuadrupole_YY*gqyy7+atomJ.labFrameQuadrupole_ZZ*gqzz7
                      +2.0f*(atomJ.labFrameQuadrupole_XY*gqxy7+atomJ.labFrameQuadrupole_XZ*gqxz7+atomJ.labFrameQuadrupole_YZ*gqyz7));

        float fidg22 =
                - 0.5f*(atomJ.q*gqyy1+atomJ.labFrameDipole[0]*gqyy2+atomJ.labFrameDipole[1]*gqyy3+atomJ.labFrameDipole[2]*gqyy4
                      +atomJ.labFrameQuadrupole_XX*gqyy5+atomJ.labFrameQuadrupole_YY*gqyy8+atomJ.labFrameQuadrupole_ZZ*gqyy10
                      +2.0f*(atomJ.labFrameQuadrupole_XY*gqyy6+atomJ.labFrameQuadrupole_XZ*gqyy7+atomJ.labFrameQuadrupole_YZ*gqyy9)
                      +atomJ.q*gc8+atomJ.labFrameDipole[0]*gux8+atomJ.labFrameDipole[1]*guy8+atomJ.labFrameDipole[2]*guz8
                      +atomJ.labFrameQuadrupole_XX*gqxx8+atomJ.labFrameQuadrupole_YY*gqyy8+atomJ.labFrameQuadrupole_ZZ*gqzz8
                      +2.0f*(atomJ.labFrameQuadrupole_XY*gqxy8+atomJ.labFrameQuadrupole_XZ*gqxz8+atomJ.labFrameQuadrupole_YZ*gqyz8));

        float fidg23 =
                - 0.5f*(atomJ.q*gqyz1+atomJ.labFrameDipole[0]*gqyz2+atomJ.labFrameDipole[1]*gqyz3+atomJ.labFrameDipole[2]*gqyz4
                      +atomJ.labFrameQuadrupole_XX*gqyz5+atomJ.labFrameQuadrupole_YY*gqyz8+atomJ.labFrameQuadrupole_ZZ*gqyz10
                      +2.0f*(atomJ.labFrameQuadrupole_XY*gqyz6+atomJ.labFrameQuadrupole_XZ*gqyz7+atomJ.labFrameQuadrupole_YZ*gqyz9)
                      +atomJ.q*gc9+atomJ.labFrameDipole[0]*gux9+atomJ.labFrameDipole[1]*guy9+atomJ.labFrameDipole[2]*guz9
                      +atomJ.labFrameQuadrupole_XX*gqxx9+atomJ.labFrameQuadrupole_YY*gqyy9+atomJ.labFrameQuadrupole_ZZ*gqzz9
                      +2.0f*(atomJ.labFrameQuadrupole_XY*gqxy9+atomJ.labFrameQuadrupole_XZ*gqxz9+atomJ.labFrameQuadrupole_YZ*gqyz9));

        float fidg33 =
                - 0.5f*(atomJ.q*gqzz1+atomJ.labFrameDipole[0]*gqzz2+atomJ.labFrameDipole[1]*gqzz3+atomJ.labFrameDipole[2]*gqzz4
                      +atomJ.labFrameQuadrupole_XX*gqzz5+atomJ.labFrameQuadrupole_YY*gqzz8+atomJ.labFrameQuadrupole_ZZ*gqzz10
                      +2.0f*(atomJ.labFrameQuadrupole_XY*gqzz6+atomJ.labFrameQuadrupole_XZ*gqzz7+atomJ.labFrameQuadrupole_YZ*gqzz9)
                      +atomJ.q*gc10+atomJ.labFrameDipole[0]*gux10+atomJ.labFrameDipole[1]*guy10+atomJ.labFrameDipole[2]*guz10
                      +atomJ.labFrameQuadrupole_XX*gqxx10+atomJ.labFrameQuadrupole_YY*gqyy10+atomJ.labFrameQuadrupole_ZZ*gqzz10
                   +2.0f*(atomJ.labFrameQuadrupole_XY*gqxy10+atomJ.labFrameQuadrupole_XZ*gqxz10+atomJ.labFrameQuadrupole_YZ*gqyz10));

        float fidg21 = fidg12;
        float fidg31 = fidg13;
        float fidg32 = fidg23;

        float fkdg11 =
                - 0.5f*(atomI.q*gqxx1-atomI.labFrameDipole[0]*gqxx2-atomI.labFrameDipole[1]*gqxx3-atomI.labFrameDipole[2] *gqxx4
                      +atomI.labFrameQuadrupole_XX*gqxx5+atomI.labFrameQuadrupole_YY*gqxx8+atomI.labFrameQuadrupole_ZZ*gqxx10
                      +2.0f*(atomI.labFrameQuadrupole_XY*gqxx6+atomI.labFrameQuadrupole_XZ*gqxx7+atomI.labFrameQuadrupole_YZ*gqxx9)
                      +atomI.q*gc5-atomI.labFrameDipole[0]*gux5-atomI.labFrameDipole[1]*guy5-atomI.labFrameDipole[2]*guz5
                      +atomI.labFrameQuadrupole_XX*gqxx5+atomI.labFrameQuadrupole_YY*gqyy5+atomI.labFrameQuadrupole_ZZ*gqzz5
                      +2.0f*(atomI.labFrameQuadrupole_XY*gqxy5+atomI.labFrameQuadrupole_XZ*gqxz5+atomI.labFrameQuadrupole_YZ*gqyz5));

        float fkdg12 =
                - 0.5f*(atomI.q*gqxy1-atomI.labFrameDipole[0]*gqxy2-atomI.labFrameDipole[1]*gqxy3-atomI.labFrameDipole[2]*gqxy4
                      +atomI.labFrameQuadrupole_XX*gqxy5+atomI.labFrameQuadrupole_YY*gqxy8+atomI.labFrameQuadrupole_ZZ*gqxy10
                      +2.0f*(atomI.labFrameQuadrupole_XY*gqxy6+atomI.labFrameQuadrupole_XZ*gqxy7+atomI.labFrameQuadrupole_YZ*gqxy9)
                      +atomI.q*gc6-atomI.labFrameDipole[0]*gux6-atomI.labFrameDipole[1]*guy6-atomI.labFrameDipole[2]*guz6
                      +atomI.labFrameQuadrupole_XX*gqxx6+atomI.labFrameQuadrupole_YY*gqyy6+atomI.labFrameQuadrupole_ZZ*gqzz6
                      +2.0f*(atomI.labFrameQuadrupole_XY*gqxy6+atomI.labFrameQuadrupole_XZ*gqxz6+atomI.labFrameQuadrupole_YZ*gqyz6));

        float fkdg13 =
                - 0.5f*(atomI.q*gqxz1-atomI.labFrameDipole[0]*gqxz2-atomI.labFrameDipole[1]*gqxz3-atomI.labFrameDipole[2]*gqxz4
                      +atomI.labFrameQuadrupole_XX*gqxz5+atomI.labFrameQuadrupole_YY*gqxz8+atomI.labFrameQuadrupole_ZZ*gqxz10
                      +2.0f*(atomI.labFrameQuadrupole_XY*gqxz6+atomI.labFrameQuadrupole_XZ*gqxz7+atomI.labFrameQuadrupole_YZ*gqxz9)
                      +atomI.q*gc7-atomI.labFrameDipole[0]*gux7-atomI.labFrameDipole[1]*guy7-atomI.labFrameDipole[2]*guz7
                      +atomI.labFrameQuadrupole_XX*gqxx7+atomI.labFrameQuadrupole_YY*gqyy7+atomI.labFrameQuadrupole_ZZ*gqzz7
                      +2.0f*(atomI.labFrameQuadrupole_XY*gqxy7+atomI.labFrameQuadrupole_XZ*gqxz7+atomI.labFrameQuadrupole_YZ*gqyz7));

        float fkdg22 =
                - 0.5f*(atomI.q*gqyy1-atomI.labFrameDipole[0]*gqyy2-atomI.labFrameDipole[1]*gqyy3-atomI.labFrameDipole[2]*gqyy4
                      +atomI.labFrameQuadrupole_XX*gqyy5+atomI.labFrameQuadrupole_YY*gqyy8+atomI.labFrameQuadrupole_ZZ*gqyy10
                      +2.0f*(atomI.labFrameQuadrupole_XY*gqyy6+atomI.labFrameQuadrupole_XZ*gqyy7+atomI.labFrameQuadrupole_YZ*gqyy9)
                      +atomI.q*gc8-atomI.labFrameDipole[0]*gux8-atomI.labFrameDipole[1]*guy8-atomI.labFrameDipole[2]*guz8
                      +atomI.labFrameQuadrupole_XX*gqxx8+atomI.labFrameQuadrupole_YY*gqyy8+atomI.labFrameQuadrupole_ZZ*gqzz8
                      +2.0f*(atomI.labFrameQuadrupole_XY*gqxy8+atomI.labFrameQuadrupole_XZ*gqxz8+atomI.labFrameQuadrupole_YZ*gqyz8));

        float fkdg23 =
                - 0.5f*(atomI.q*gqyz1-atomI.labFrameDipole[0]*gqyz2-atomI.labFrameDipole[1]*gqyz3-atomI.labFrameDipole[2]*gqyz4
                      +atomI.labFrameQuadrupole_XX*gqyz5+atomI.labFrameQuadrupole_YY*gqyz8+atomI.labFrameQuadrupole_ZZ*gqyz10
                      +2.0f*(atomI.labFrameQuadrupole_XY*gqyz6+atomI.labFrameQuadrupole_XZ*gqyz7+atomI.labFrameQuadrupole_YZ*gqyz9)
                      +atomI.q*gc9-atomI.labFrameDipole[0]*gux9-atomI.labFrameDipole[1]*guy9-atomI.labFrameDipole[2]*guz9
                      +atomI.labFrameQuadrupole_XX*gqxx9+atomI.labFrameQuadrupole_YY*gqyy9+atomI.labFrameQuadrupole_ZZ*gqzz9
                      +2.0f*(atomI.labFrameQuadrupole_XY*gqxy9+atomI.labFrameQuadrupole_XZ*gqxz9+atomI.labFrameQuadrupole_YZ*gqyz9));
        float fkdg33 =
                - 0.5f*(atomI.q*gqzz1-atomI.labFrameDipole[0]*gqzz2-atomI.labFrameDipole[1]*gqzz3-atomI.labFrameDipole[2]*gqzz4
                      +atomI.labFrameQuadrupole_XX*gqzz5+atomI.labFrameQuadrupole_YY*gqzz8+atomI.labFrameQuadrupole_ZZ*gqzz10
                      +2.0f*(atomI.labFrameQuadrupole_XY*gqzz6+atomI.labFrameQuadrupole_XZ*gqzz7+atomI.labFrameQuadrupole_YZ*gqzz9)
                      +atomI.q*gc10-atomI.labFrameDipole[0]*gux10-atomI.labFrameDipole[1]*guy10-atomI.labFrameDipole[2]*guz10
                      +atomI.labFrameQuadrupole_XX*gqxx10+atomI.labFrameQuadrupole_YY*gqyy10+atomI.labFrameQuadrupole_ZZ*gqzz10
                    +2.0f*(atomI.labFrameQuadrupole_XY*gqxy10+atomI.labFrameQuadrupole_XZ*gqxz10+atomI.labFrameQuadrupole_YZ*gqyz10));

        float fkdg21 = fkdg12;
        float fkdg31 = fkdg13;
        float fkdg32 = fkdg23;

        trq1   += 2.0f* (atomI.labFrameQuadrupole_XY*fidg13+atomI.labFrameQuadrupole_YY*fidg23+atomI.labFrameQuadrupole_YZ*fidg33
                           -atomI.labFrameQuadrupole_XZ*fidg12-atomI.labFrameQuadrupole_YZ*fidg22-atomI.labFrameQuadrupole_ZZ*fidg32);

        trq2   += 2.0f*(atomI.labFrameQuadrupole_XZ*fidg11+atomI.labFrameQuadrupole_YZ*fidg21+atomI.labFrameQuadrupole_ZZ*fidg31
                         -atomI.labFrameQuadrupole_XX*fidg13-atomI.labFrameQuadrupole_XY*fidg23-atomI.labFrameQuadrupole_XZ*fidg33);

        trq3   += 2.0f*(atomI.labFrameQuadrupole_XX*fidg12+atomI.labFrameQuadrupole_XY*fidg22+atomI.labFrameQuadrupole_XZ*fidg32
                         -atomI.labFrameQuadrupole_XY*fidg11-atomI.labFrameQuadrupole_YY*fidg21-atomI.labFrameQuadrupole_YZ*fidg31);

        trq_k1 += 2.0f*
                          (atomJ.labFrameQuadrupole_XY*fkdg13+atomJ.labFrameQuadrupole_YY*fkdg23+atomJ.labFrameQuadrupole_YZ*fkdg33
                          -atomJ.labFrameQuadrupole_XZ*fkdg12-atomJ.labFrameQuadrupole_YZ*fkdg22-atomJ.labFrameQuadrupole_ZZ*fkdg32);

        trq_k2 += 2.0f*
                          (atomJ.labFrameQuadrupole_XZ*fkdg11+atomJ.labFrameQuadrupole_YZ*fkdg21+atomJ.labFrameQuadrupole_ZZ*fkdg31
                          -atomJ.labFrameQuadrupole_XX*fkdg13-atomJ.labFrameQuadrupole_XY*fkdg23-atomJ.labFrameQuadrupole_XZ*fkdg33);

        trq_k3 += 2.0f*
                          (atomJ.labFrameQuadrupole_XX*fkdg12+atomJ.labFrameQuadrupole_XY*fkdg22+atomJ.labFrameQuadrupole_XZ*fkdg32
                          -atomJ.labFrameQuadrupole_XY*fkdg11-atomJ.labFrameQuadrupole_YY*fkdg21-atomJ.labFrameQuadrupole_YZ*fkdg31);
    }

    // electrostatic solvation energy of the permanent multipoles in
    // the GK reaction potential of the induced dipoles

    esymi =              -atomI.labFrameDipole[0]*(atomJ.inducedDipole[0]*gux2+atomJ.inducedDipole[1]*guy2+atomJ.inducedDipole[2]*guz2)
                        - atomI.labFrameDipole[1]*(atomJ.inducedDipole[0]*gux3+atomJ.inducedDipole[1]*guy3+atomJ.inducedDipole[2]*guz3)
                        - atomI.labFrameDipole[2]*(atomJ.inducedDipole[0]*gux4+atomJ.inducedDipole[1]*guy4+atomJ.inducedDipole[2]*guz4)
                        - atomJ.labFrameDipole[0]*(atomI.inducedDipole[0]*gux2+atomI.inducedDipole[1]*guy2+atomI.inducedDipole[2]*guz2)
                        - atomJ.labFrameDipole[1]*(atomI.inducedDipole[0]*gux3+atomI.inducedDipole[1]*guy3+atomI.inducedDipole[2]*guz3)
                        - atomJ.labFrameDipole[2]*(atomI.inducedDipole[0]*gux4+atomI.inducedDipole[1]*guy4+atomI.inducedDipole[2]*guz4);

    ewii = atomI.q*(atomJ.inducedDipole[0]*gc2+atomJ.inducedDipole[1]*gc3+atomJ.inducedDipole[2]*gc4)
                      - atomJ.q*(atomI.inducedDipole[0]*gux1+atomI.inducedDipole[1]*guy1+atomI.inducedDipole[2]*guz1)
                      - atomI.inducedDipole[0]*(atomJ.labFrameQuadrupole_XX*gux5+atomJ.labFrameQuadrupole_YY*gux8+atomJ.labFrameQuadrupole_ZZ*gux10
                     +2.0f*(atomJ.labFrameQuadrupole_XY*gux6+atomJ.labFrameQuadrupole_XZ*gux7+atomJ.labFrameQuadrupole_YZ*gux9))
                      - atomI.inducedDipole[1]*(atomJ.labFrameQuadrupole_XX*guy5+atomJ.labFrameQuadrupole_YY*guy8+atomJ.labFrameQuadrupole_ZZ*guy10
                     +2.0f*(atomJ.labFrameQuadrupole_XY*guy6+atomJ.labFrameQuadrupole_XZ*guy7+atomJ.labFrameQuadrupole_YZ*guy9))
                      - atomI.inducedDipole[2]*(atomJ.labFrameQuadrupole_XX*guz5+atomJ.labFrameQuadrupole_YY*guz8+atomJ.labFrameQuadrupole_ZZ*guz10
                     +2.0f*(atomJ.labFrameQuadrupole_XY*guz6+atomJ.labFrameQuadrupole_XZ*guz7+atomJ.labFrameQuadrupole_YZ*guz9))
                      + atomJ.inducedDipole[0]*(atomI.labFrameQuadrupole_XX*gqxx2+atomI.labFrameQuadrupole_YY*gqyy2+atomI.labFrameQuadrupole_ZZ*gqzz2
                     +2.0f*(atomI.labFrameQuadrupole_XY*gqxy2+atomI.labFrameQuadrupole_XZ*gqxz2+atomI.labFrameQuadrupole_YZ*gqyz2))
                      + atomJ.inducedDipole[1]*(atomI.labFrameQuadrupole_XX*gqxx3+atomI.labFrameQuadrupole_YY*gqyy3+atomI.labFrameQuadrupole_ZZ*gqzz3
                     +2.0f*(atomI.labFrameQuadrupole_XY*gqxy3+atomI.labFrameQuadrupole_XZ*gqxz3+atomI.labFrameQuadrupole_YZ*gqyz3))
                      + atomJ.inducedDipole[2]*(atomI.labFrameQuadrupole_XX*gqxx4+atomI.labFrameQuadrupole_YY*gqyy4+atomI.labFrameQuadrupole_ZZ*gqzz4
                     +2.0f*(atomI.labFrameQuadrupole_XY*gqxy4+atomI.labFrameQuadrupole_XZ*gqxz4+atomI.labFrameQuadrupole_YZ*gqyz4));

    ewki = atomI.q*(atomJ.inducedDipole[0]*gux1+atomJ.inducedDipole[1]*guy1+atomJ.inducedDipole[2]*guz1)
                      - atomJ.q*(atomI.inducedDipole[0]*gc2+atomI.inducedDipole[1]*gc3+atomI.inducedDipole[2]*gc4)
                      - atomI.inducedDipole[0]*(atomJ.labFrameQuadrupole_XX*gqxx2+atomJ.labFrameQuadrupole_YY*gqyy2+atomJ.labFrameQuadrupole_ZZ*gqzz2
                     +2.0f*(atomJ.labFrameQuadrupole_XY*gqxy2+atomJ.labFrameQuadrupole_XZ*gqxz2+atomJ.labFrameQuadrupole_YZ*gqyz2))
                      - atomI.inducedDipole[1]*(atomJ.labFrameQuadrupole_XX*gqxx3+atomJ.labFrameQuadrupole_YY*gqyy3+atomJ.labFrameQuadrupole_ZZ*gqzz3
                     +2.0f*(atomJ.labFrameQuadrupole_XY*gqxy3+atomJ.labFrameQuadrupole_XZ*gqxz3+atomJ.labFrameQuadrupole_YZ*gqyz3))
                      - atomI.inducedDipole[2]*(atomJ.labFrameQuadrupole_XX*gqxx4+atomJ.labFrameQuadrupole_YY*gqyy4+atomJ.labFrameQuadrupole_ZZ*gqzz4
                     +2.0f*(atomJ.labFrameQuadrupole_XY*gqxy4+atomJ.labFrameQuadrupole_XZ*gqxz4+atomJ.labFrameQuadrupole_YZ*gqyz4))
                      + atomJ.inducedDipole[0]*(atomI.labFrameQuadrupole_XX*gux5+atomI.labFrameQuadrupole_YY*gux8+atomI.labFrameQuadrupole_ZZ*gux10
                     +2.0f*(atomI.labFrameQuadrupole_XY*gux6+atomI.labFrameQuadrupole_XZ*gux7+atomI.labFrameQuadrupole_YZ*gux9))
                      + atomJ.inducedDipole[1]*(atomI.labFrameQuadrupole_XX*guy5+atomI.labFrameQuadrupole_YY*guy8+atomI.labFrameQuadrupole_ZZ*guy10
                     +2.0f*(atomI.labFrameQuadrupole_XY*guy6+atomI.labFrameQuadrupole_XZ*guy7+atomI.labFrameQuadrupole_YZ*guy9))
                      + atomJ.inducedDipole[2]*(atomI.labFrameQuadrupole_XX*guz5+atomI.labFrameQuadrupole_YY*guz8+atomI.labFrameQuadrupole_ZZ*guz10
                     +2.0f*(atomI.labFrameQuadrupole_XY*guz6+atomI.labFrameQuadrupole_XZ*guz7+atomI.labFrameQuadrupole_YZ*guz9));

    // electrostatic solvation free energy gradient of the permanent
    // multipoles in the reaction potential of the induced dipoles

    dpsymdx = -atomI.labFrameDipole[0]*(sxk*gux5+syk*guy5+szk*guz5)
                          - atomI.labFrameDipole[1]*(sxk*gux6+syk*guy6+szk*guz6)
                          - atomI.labFrameDipole[2]*(sxk*gux7+syk*guy7+szk*guz7)
                          - atomJ.labFrameDipole[0]*(sxi*gux5+syi*guy5+szi*guz5)
                          - atomJ.labFrameDipole[1]*(sxi*gux6+syi*guy6+szi*guz6)
                          - atomJ.labFrameDipole[2]*(sxi*gux7+syi*guy7+szi*guz7);

    dpwidx = atomI.q*(sxk*gc5+syk*gc6+szk*gc7)
                        - atomJ.q*(sxi*gux2+syi*guy2+szi*guz2)
                      - sxi*(atomJ.labFrameQuadrupole_XX*gux11+atomJ.labFrameQuadrupole_YY*gux14+atomJ.labFrameQuadrupole_ZZ*gux16
                     +2.0f*(atomJ.labFrameQuadrupole_XY*gux12+atomJ.labFrameQuadrupole_XZ*gux13+atomJ.labFrameQuadrupole_YZ*gux15))
                      - syi*(atomJ.labFrameQuadrupole_XX*guy11+atomJ.labFrameQuadrupole_YY*guy14+atomJ.labFrameQuadrupole_ZZ*guy16
                     +2.0f*(atomJ.labFrameQuadrupole_XY*guy12+atomJ.labFrameQuadrupole_XZ*guy13+atomJ.labFrameQuadrupole_YZ*guy15))
                      - szi*(atomJ.labFrameQuadrupole_XX*guz11+atomJ.labFrameQuadrupole_YY*guz14+atomJ.labFrameQuadrupole_ZZ*guz16
                     +2.0f*(atomJ.labFrameQuadrupole_XY*guz12+atomJ.labFrameQuadrupole_XZ*guz13+atomJ.labFrameQuadrupole_YZ*guz15))
                      + sxk*(atomI.labFrameQuadrupole_XX*gqxx5+atomI.labFrameQuadrupole_YY*gqyy5+atomI.labFrameQuadrupole_ZZ*gqzz5
                     +2.0f*(atomI.labFrameQuadrupole_XY*gqxy5+atomI.labFrameQuadrupole_XZ*gqxz5+atomI.labFrameQuadrupole_YZ*gqyz5))
                      + syk*(atomI.labFrameQuadrupole_XX*gqxx6+atomI.labFrameQuadrupole_YY*gqyy6+atomI.labFrameQuadrupole_ZZ*gqzz6
                     +2.0f*(atomI.labFrameQuadrupole_XY*gqxy6+atomI.labFrameQuadrupole_XZ*gqxz6+atomI.labFrameQuadrupole_YZ*gqyz6))
                      + szk*(atomI.labFrameQuadrupole_XX*gqxx7+atomI.labFrameQuadrupole_YY*gqyy7+atomI.labFrameQuadrupole_ZZ*gqzz7
                     +2.0f*(atomI.labFrameQuadrupole_XY*gqxy7+atomI.labFrameQuadrupole_XZ*gqxz7+atomI.labFrameQuadrupole_YZ*gqyz7));

    dpwkdx = atomI.q*(sxk*gux2+syk*guy2+szk*guz2)
                        - atomJ.q*(sxi*gc5+syi*gc6+szi*gc7)
                      - sxi*(atomJ.labFrameQuadrupole_XX*gqxx5+atomJ.labFrameQuadrupole_YY*gqyy5+atomJ.labFrameQuadrupole_ZZ*gqzz5
                     +2.0f*(atomJ.labFrameQuadrupole_XY*gqxy5+atomJ.labFrameQuadrupole_XZ*gqxz5+atomJ.labFrameQuadrupole_YZ*gqyz5))
                      - syi*(atomJ.labFrameQuadrupole_XX*gqxx6+atomJ.labFrameQuadrupole_YY*gqyy6+atomJ.labFrameQuadrupole_ZZ*gqzz6
                     +2.0f*(atomJ.labFrameQuadrupole_XY*gqxy6+atomJ.labFrameQuadrupole_XZ*gqxz6+atomJ.labFrameQuadrupole_YZ*gqyz6))
                      - szi*(atomJ.labFrameQuadrupole_XX*gqxx7+atomJ.labFrameQuadrupole_YY*gqyy7+atomJ.labFrameQuadrupole_ZZ*gqzz7
                     +2.0f*(atomJ.labFrameQuadrupole_XY*gqxy7+atomJ.labFrameQuadrupole_XZ*gqxz7+atomJ.labFrameQuadrupole_YZ*gqyz7))
                      + sxk*(atomI.labFrameQuadrupole_XX*gux11+atomI.labFrameQuadrupole_YY*gux14+atomI.labFrameQuadrupole_ZZ*gux16
                     +2.0f*(atomI.labFrameQuadrupole_XY*gux12+atomI.labFrameQuadrupole_XZ*gux13+atomI.labFrameQuadrupole_YZ*gux15))
                      + syk*(atomI.labFrameQuadrupole_XX*guy11+atomI.labFrameQuadrupole_YY*guy14+atomI.labFrameQuadrupole_ZZ*guy16
                     +2.0f*(atomI.labFrameQuadrupole_XY*guy12+atomI.labFrameQuadrupole_XZ*guy13+atomI.labFrameQuadrupole_YZ*guy15))
                      + szk*(atomI.labFrameQuadrupole_XX*guz11+atomI.labFrameQuadrupole_YY*guz14+atomI.labFrameQuadrupole_ZZ*guz16
                     +2.0f*(atomI.labFrameQuadrupole_XY*guz12+atomI.labFrameQuadrupole_XZ*guz13+atomI.labFrameQuadrupole_YZ*guz15));

    dpdx = 0.5f * (dpsymdx + 0.5f*(dpwidx + dpwkdx));

    dpsymdy = -atomI.labFrameDipole[0]*(sxk*gux6+syk*guy6+szk*guz6)
                          - atomI.labFrameDipole[1]*(sxk*gux8+syk*guy8+szk*guz8)
                          - atomI.labFrameDipole[2]*(sxk*gux9+syk*guy9+szk*guz9)
                          - atomJ.labFrameDipole[0]*(sxi*gux6+syi*guy6+szi*guz6)
                          - atomJ.labFrameDipole[1]*(sxi*gux8+syi*guy8+szi*guz8)
                          - atomJ.labFrameDipole[2]*(sxi*gux9+syi*guy9+szi*guz9);

    dpwidy = atomI.q*(sxk*gc6+syk*gc8+szk*gc9)
                        - atomJ.q*(sxi*gux3+syi*guy3+szi*guz3)
                         - sxi*(atomJ.labFrameQuadrupole_XX*gux12+atomJ.labFrameQuadrupole_YY*gux17+atomJ.labFrameQuadrupole_ZZ*gux19
                        +2.0f*(atomJ.labFrameQuadrupole_XY*gux14+atomJ.labFrameQuadrupole_XZ*gux15+atomJ.labFrameQuadrupole_YZ*gux18))
                         - syi*(atomJ.labFrameQuadrupole_XX*guy12+atomJ.labFrameQuadrupole_YY*guy17+atomJ.labFrameQuadrupole_ZZ*guy19
                        +2.0f*(atomJ.labFrameQuadrupole_XY*guy14+atomJ.labFrameQuadrupole_XZ*guy15+atomJ.labFrameQuadrupole_YZ*guy18))
                         - szi*(atomJ.labFrameQuadrupole_XX*guz12+atomJ.labFrameQuadrupole_YY*guz17+atomJ.labFrameQuadrupole_ZZ*guz19
                        +2.0f*(atomJ.labFrameQuadrupole_XY*guz14+atomJ.labFrameQuadrupole_XZ*guz15+atomJ.labFrameQuadrupole_YZ*guz18))
                         + sxk*(atomI.labFrameQuadrupole_XX*gqxx6+atomI.labFrameQuadrupole_YY*gqyy6+atomI.labFrameQuadrupole_ZZ*gqzz6
                        +2.0f*(atomI.labFrameQuadrupole_XY*gqxy6+atomI.labFrameQuadrupole_XZ*gqxz6+atomI.labFrameQuadrupole_YZ*gqyz6))
                         + syk*(atomI.labFrameQuadrupole_XX*gqxx8+atomI.labFrameQuadrupole_YY*gqyy8+atomI.labFrameQuadrupole_ZZ*gqzz8
                        +2.0f*(atomI.labFrameQuadrupole_XY*gqxy8+atomI.labFrameQuadrupole_XZ*gqxz8+atomI.labFrameQuadrupole_YZ*gqyz8))
                         + szk*(atomI.labFrameQuadrupole_XX*gqxx9+atomI.labFrameQuadrupole_YY*gqyy9+atomI.labFrameQuadrupole_ZZ*gqzz9
                        +2.0f*(atomI.labFrameQuadrupole_XY*gqxy9+atomI.labFrameQuadrupole_XZ*gqxz9+atomI.labFrameQuadrupole_YZ*gqyz9));

    dpwkdy = atomI.q*(sxk*gux3+syk*guy3+szk*guz3)
                        - atomJ.q*(sxi*gc6+syi*gc8+szi*gc9)
                      - sxi*(atomJ.labFrameQuadrupole_XX*gqxx6+atomJ.labFrameQuadrupole_YY*gqyy6+atomJ.labFrameQuadrupole_ZZ*gqzz6
                     +2.0f*(atomJ.labFrameQuadrupole_XY*gqxy6+atomJ.labFrameQuadrupole_XZ*gqxz6+atomJ.labFrameQuadrupole_YZ*gqyz6))
                      - syi*(atomJ.labFrameQuadrupole_XX*gqxx8+atomJ.labFrameQuadrupole_YY*gqyy8+atomJ.labFrameQuadrupole_ZZ*gqzz8
                     +2.0f*(atomJ.labFrameQuadrupole_XY*gqxy8+atomJ.labFrameQuadrupole_XZ*gqxz8+atomJ.labFrameQuadrupole_YZ*gqyz8))
                      - szi*(atomJ.labFrameQuadrupole_XX*gqxx9+atomJ.labFrameQuadrupole_YY*gqyy9+atomJ.labFrameQuadrupole_ZZ*gqzz9
                     +2.0f*(atomJ.labFrameQuadrupole_XY*gqxy9+atomJ.labFrameQuadrupole_XZ*gqxz9+atomJ.labFrameQuadrupole_YZ*gqyz9))
                      + sxk*(atomI.labFrameQuadrupole_XX*gux12+atomI.labFrameQuadrupole_YY*gux17+atomI.labFrameQuadrupole_ZZ*gux19
                     +2.0f*(atomI.labFrameQuadrupole_XY*gux14+atomI.labFrameQuadrupole_XZ*gux15+atomI.labFrameQuadrupole_YZ*gux18))
                      + syk*(atomI.labFrameQuadrupole_XX*guy12+atomI.labFrameQuadrupole_YY*guy17+atomI.labFrameQuadrupole_ZZ*guy19
                     +2.0f*(atomI.labFrameQuadrupole_XY*guy14+atomI.labFrameQuadrupole_XZ*guy15+atomI.labFrameQuadrupole_YZ*guy18))
                      + szk*(atomI.labFrameQuadrupole_XX*guz12+atomI.labFrameQuadrupole_YY*guz17+atomI.labFrameQuadrupole_ZZ*guz19
                     +2.0f*(atomI.labFrameQuadrupole_XY*guz14+atomI.labFrameQuadrupole_XZ*guz15+atomI.labFrameQuadrupole_YZ*guz18));

    dpdy    = 0.5f * (dpsymdy + 0.5f*(dpwidy + dpwkdy));

    dpsymdz = -atomI.labFrameDipole[0]*(sxk*gux7+syk*guy7+szk*guz7)
                          - atomI.labFrameDipole[1]*(sxk*gux9+syk*guy9+szk*guz9)
                          - atomI.labFrameDipole[2]*(sxk*gux10+syk*guy10+szk*guz10)
                          - atomJ.labFrameDipole[0]*(sxi*gux7+syi*guy7+szi*guz7)
                          - atomJ.labFrameDipole[1]*(sxi*gux9+syi*guy9+szi*guz9)
                          - atomJ.labFrameDipole[2]*(sxi*gux10+syi*guy10+szi*guz10);

    dpwidz = atomI.q*(sxk*gc7+syk*gc9+szk*gc10)
                        - atomJ.q*(sxi*gux4+syi*guy4+szi*guz4)
                      - sxi*(atomJ.labFrameQuadrupole_XX*gux13+atomJ.labFrameQuadrupole_YY*gux18+atomJ.labFrameQuadrupole_ZZ*gux20
                     +2.0f*(atomJ.labFrameQuadrupole_XY*gux15+atomJ.labFrameQuadrupole_XZ*gux16+atomJ.labFrameQuadrupole_YZ*gux19))
                      - syi*(atomJ.labFrameQuadrupole_XX*guy13+atomJ.labFrameQuadrupole_YY*guy18+atomJ.labFrameQuadrupole_ZZ*guy20
                     +2.0f*(atomJ.labFrameQuadrupole_XY*guy15+atomJ.labFrameQuadrupole_XZ*guy16+atomJ.labFrameQuadrupole_YZ*guy19))
                      - szi*(atomJ.labFrameQuadrupole_XX*guz13+atomJ.labFrameQuadrupole_YY*guz18+atomJ.labFrameQuadrupole_ZZ*guz20
                     +2.0f*(atomJ.labFrameQuadrupole_XY*guz15+atomJ.labFrameQuadrupole_XZ*guz16+atomJ.labFrameQuadrupole_YZ*guz19))
                      + sxk*(atomI.labFrameQuadrupole_XX*gqxx7+atomI.labFrameQuadrupole_YY*gqyy7+atomI.labFrameQuadrupole_ZZ*gqzz7
                     +2.0f*(atomI.labFrameQuadrupole_XY*gqxy7+atomI.labFrameQuadrupole_XZ*gqxz7+atomI.labFrameQuadrupole_YZ*gqyz7))
                      + syk*(atomI.labFrameQuadrupole_XX*gqxx9+atomI.labFrameQuadrupole_YY*gqyy9+atomI.labFrameQuadrupole_ZZ*gqzz9
                     +2.0f*(atomI.labFrameQuadrupole_XY*gqxy9+atomI.labFrameQuadrupole_XZ*gqxz9+atomI.labFrameQuadrupole_YZ*gqyz9))
                      + szk*(atomI.labFrameQuadrupole_XX*gqxx10+atomI.labFrameQuadrupole_YY*gqyy10+atomI.labFrameQuadrupole_ZZ*gqzz10
                     +2.0f*(atomI.labFrameQuadrupole_XY*gqxy10+atomI.labFrameQuadrupole_XZ*gqxz10+atomI.labFrameQuadrupole_YZ*gqyz10));

    dpwkdz = atomI.q*(sxk*gux4+syk*guy4+szk*guz4)
                        - atomJ.q*(sxi*gc7+syi*gc9+szi*gc10)
                      - sxi*(atomJ.labFrameQuadrupole_XX*gqxx7+atomJ.labFrameQuadrupole_YY*gqyy7+atomJ.labFrameQuadrupole_ZZ*gqzz7
                     +2.0f*(atomJ.labFrameQuadrupole_XY*gqxy7+atomJ.labFrameQuadrupole_XZ*gqxz7+atomJ.labFrameQuadrupole_YZ*gqyz7))
                      - syi*(atomJ.labFrameQuadrupole_XX*gqxx9+atomJ.labFrameQuadrupole_YY*gqyy9+atomJ.labFrameQuadrupole_ZZ*gqzz9
                     +2.0f*(atomJ.labFrameQuadrupole_XY*gqxy9+atomJ.labFrameQuadrupole_XZ*gqxz9+atomJ.labFrameQuadrupole_YZ*gqyz9))
                      - szi*(atomJ.labFrameQuadrupole_XX*gqxx10+atomJ.labFrameQuadrupole_YY*gqyy10+atomJ.labFrameQuadrupole_ZZ*gqzz10
                     +2.0f*(atomJ.labFrameQuadrupole_XY*gqxy10+atomJ.labFrameQuadrupole_XZ*gqxz10+atomJ.labFrameQuadrupole_YZ*gqyz10))
                      + sxk*(atomI.labFrameQuadrupole_XX*gux13+atomI.labFrameQuadrupole_YY*gux18+atomI.labFrameQuadrupole_ZZ*gux20
                     +2.0f*(atomI.labFrameQuadrupole_XY*gux15+atomI.labFrameQuadrupole_XZ*gux16+atomI.labFrameQuadrupole_YZ*gux19))
                      + syk*(atomI.labFrameQuadrupole_XX*guy13+atomI.labFrameQuadrupole_YY*guy18+atomI.labFrameQuadrupole_ZZ*guy20
                     +2.0f*(atomI.labFrameQuadrupole_XY*guy15+atomI.labFrameQuadrupole_XZ*guy16+atomI.labFrameQuadrupole_YZ*guy19))
                      + szk*(atomI.labFrameQuadrupole_XX*guz13+atomI.labFrameQuadrupole_YY*guz18+atomI.labFrameQuadrupole_ZZ*guz20
                     +2.0f*(atomI.labFrameQuadrupole_XY*guz15+atomI.labFrameQuadrupole_XZ*guz16+atomI.labFrameQuadrupole_YZ*guz19));

    dpdz = 0.5f * (dpsymdz + 0.5f*(dpwidz + dpwkdz));

    // effective radii chain rule terms for the;
    // electrostatic solvation free energy gradient of the permanent;
    // multipoles in the reaction potential of the induced dipoles;

    dsymdr = -atomI.labFrameDipole[0]*(sxk*gux22+syk*guy22+szk*guz22)
                          - atomI.labFrameDipole[1]*(sxk*gux23+syk*guy23+szk*guz23)
                          - atomI.labFrameDipole[2]*(sxk*gux24+syk*guy24+szk*guz24)
                          - atomJ.labFrameDipole[0]*(sxi*gux22+syi*guy22+szi*guz22)
                          - atomJ.labFrameDipole[1]*(sxi*gux23+syi*guy23+szi*guz23)
                          - atomJ.labFrameDipole[2]*(sxi*gux24+syi*guy24+szi*guz24);

    dwipdr = atomI.q*(sxk*gc22+syk*gc23+szk*gc24)
                         - atomJ.q*(sxi*gux21+syi*guy21+szi*guz21)
                      - sxi*(atomJ.labFrameQuadrupole_XX*gux25+atomJ.labFrameQuadrupole_YY*gux28+atomJ.labFrameQuadrupole_ZZ*gux30
                     +2.0f*(atomJ.labFrameQuadrupole_XY*gux26+atomJ.labFrameQuadrupole_XZ*gux27+atomJ.labFrameQuadrupole_YZ*gux29))
                      - syi*(atomJ.labFrameQuadrupole_XX*guy25+atomJ.labFrameQuadrupole_YY*guy28+atomJ.labFrameQuadrupole_ZZ*guy30
                     +2.0f*(atomJ.labFrameQuadrupole_XY*guy26+atomJ.labFrameQuadrupole_XZ*guy27+atomJ.labFrameQuadrupole_YZ*guy29))
                      - szi*(atomJ.labFrameQuadrupole_XX*guz25+atomJ.labFrameQuadrupole_YY*guz28+atomJ.labFrameQuadrupole_ZZ*guz30
                     +2.0f*(atomJ.labFrameQuadrupole_XY*guz26+atomJ.labFrameQuadrupole_XZ*guz27+atomJ.labFrameQuadrupole_YZ*guz29))
                      + sxk*(atomI.labFrameQuadrupole_XX*gqxx22+atomI.labFrameQuadrupole_YY*gqyy22+atomI.labFrameQuadrupole_ZZ*gqzz22
                     +2.0f*(atomI.labFrameQuadrupole_XY*gqxy22+atomI.labFrameQuadrupole_XZ*gqxz22+atomI.labFrameQuadrupole_YZ*gqyz22))
                      + syk*(atomI.labFrameQuadrupole_XX*gqxx23+atomI.labFrameQuadrupole_YY*gqyy23+atomI.labFrameQuadrupole_ZZ*gqzz23
                     +2.0f*(atomI.labFrameQuadrupole_XY*gqxy23+atomI.labFrameQuadrupole_XZ*gqxz23+atomI.labFrameQuadrupole_YZ*gqyz23))
                      + szk*(atomI.labFrameQuadrupole_XX*gqxx24+atomI.labFrameQuadrupole_YY*gqyy24+atomI.labFrameQuadrupole_ZZ*gqzz24
                     +2.0f*(atomI.labFrameQuadrupole_XY*gqxy24+atomI.labFrameQuadrupole_XZ*gqxz24+atomI.labFrameQuadrupole_YZ*gqyz24));

    dwkpdr = atomI.q*(sxk*gux21+syk*guy21+szk*guz21)
                         - atomJ.q*(sxi*gc22+syi*gc23+szi*gc24)
                      - sxi*(atomJ.labFrameQuadrupole_XX*gqxx22+atomJ.labFrameQuadrupole_YY*gqyy22+atomJ.labFrameQuadrupole_ZZ*gqzz22
                     +2.0f*(atomJ.labFrameQuadrupole_XY*gqxy22+atomJ.labFrameQuadrupole_XZ*gqxz22+atomJ.labFrameQuadrupole_YZ*gqyz22))
                      - syi*(atomJ.labFrameQuadrupole_XX*gqxx23+atomJ.labFrameQuadrupole_YY*gqyy23+atomJ.labFrameQuadrupole_ZZ*gqzz23
                     +2.0f*(atomJ.labFrameQuadrupole_XY*gqxy23+atomJ.labFrameQuadrupole_XZ*gqxz23+atomJ.labFrameQuadrupole_YZ*gqyz23))
                      - szi*(atomJ.labFrameQuadrupole_XX*gqxx24+atomJ.labFrameQuadrupole_YY*gqyy24+atomJ.labFrameQuadrupole_ZZ*gqzz24
                     +2.0f*(atomJ.labFrameQuadrupole_XY*gqxy24+atomJ.labFrameQuadrupole_XZ*gqxz24+atomJ.labFrameQuadrupole_YZ*gqyz24))
                      + sxk*(atomI.labFrameQuadrupole_XX*gux25+atomI.labFrameQuadrupole_YY*gux28+atomI.labFrameQuadrupole_ZZ*gux30
                     +2.0f*(atomI.labFrameQuadrupole_XY*gux26+atomI.labFrameQuadrupole_XZ*gux27+atomI.labFrameQuadrupole_YZ*gux29))
                      + syk*(atomI.labFrameQuadrupole_XX*guy25+atomI.labFrameQuadrupole_YY*guy28+atomI.labFrameQuadrupole_ZZ*guy30
                     +2.0f*(atomI.labFrameQuadrupole_XY*guy26+atomI.labFrameQuadrupole_XZ*guy27+atomI.labFrameQuadrupole_YZ*guy29))
                      + szk*(atomI.labFrameQuadrupole_XX*guz25+atomI.labFrameQuadrupole_YY*guz28+atomI.labFrameQuadrupole_ZZ*guz30
                     +2.0f*(atomI.labFrameQuadrupole_XY*guz26+atomI.labFrameQuadrupole_XZ*guz27+atomI.labFrameQuadrupole_YZ*guz29));

    dsumdr = dsymdr + 0.5f*(dwipdr + dwkpdr);
    dpbi = 0.5f*atomJ.bornRadius*dsumdr;
    dpbk = 0.5f*atomI.bornRadius*dsumdr;

    // mutual polarization electrostatic solvation free energy gradient

//   if (poltyp .eq. 'MUTUAL'){

        dpdx = dpdx - 0.5f *
                           (atomI.inducedDipole[0]*(atomJ.inducedDipoleP[0]*gux5+atomJ.inducedDipoleP[1]*gux6+atomJ.inducedDipoleP[2]*gux7)
                           +atomI.inducedDipole[1]*(atomJ.inducedDipoleP[0]*guy5+atomJ.inducedDipoleP[1]*guy6+atomJ.inducedDipoleP[2]*guy7)
                           +atomI.inducedDipole[2]*(atomJ.inducedDipoleP[0]*guz5+atomJ.inducedDipoleP[1]*guz6+atomJ.inducedDipoleP[2]*guz7)
                           +atomJ.inducedDipole[0]*(atomI.inducedDipoleP[0]*gux5+atomI.inducedDipoleP[1]*gux6+atomI.inducedDipoleP[2]*gux7)
                           +atomJ.inducedDipole[1]*(atomI.inducedDipoleP[0]*guy5+atomI.inducedDipoleP[1]*guy6+atomI.inducedDipoleP[2]*guy7)
                           +atomJ.inducedDipole[2]*(atomI.inducedDipoleP[0]*guz5+atomI.inducedDipoleP[1]*guz6+atomI.inducedDipoleP[2]*guz7));

        dpdy = dpdy - 0.5f *
                           (atomI.inducedDipole[0]*(atomJ.inducedDipoleP[0]*gux6+atomJ.inducedDipoleP[1]*gux8+atomJ.inducedDipoleP[2]*gux9)
                           +atomI.inducedDipole[1]*(atomJ.inducedDipoleP[0]*guy6+atomJ.inducedDipoleP[1]*guy8+atomJ.inducedDipoleP[2]*guy9)
                           +atomI.inducedDipole[2]*(atomJ.inducedDipoleP[0]*guz6+atomJ.inducedDipoleP[1]*guz8+atomJ.inducedDipoleP[2]*guz9)
                           +atomJ.inducedDipole[0]*(atomI.inducedDipoleP[0]*gux6+atomI.inducedDipoleP[1]*gux8+atomI.inducedDipoleP[2]*gux9)
                           +atomJ.inducedDipole[1]*(atomI.inducedDipoleP[0]*guy6+atomI.inducedDipoleP[1]*guy8+atomI.inducedDipoleP[2]*guy9)
                           +atomJ.inducedDipole[2]*(atomI.inducedDipoleP[0]*guz6+atomI.inducedDipoleP[1]*guz8+atomI.inducedDipoleP[2]*guz9));

        dpdz = dpdz - 0.5f *
                           (atomI.inducedDipole[0]*(atomJ.inducedDipoleP[0]*gux7+atomJ.inducedDipoleP[1]*gux9+atomJ.inducedDipoleP[2]*gux10)
                           +atomI.inducedDipole[1]*(atomJ.inducedDipoleP[0]*guy7+atomJ.inducedDipoleP[1]*guy9+atomJ.inducedDipoleP[2]*guy10)
                           +atomI.inducedDipole[2]*(atomJ.inducedDipoleP[0]*guz7+atomJ.inducedDipoleP[1]*guz9+atomJ.inducedDipoleP[2]*guz10)
                           +atomJ.inducedDipole[0]*(atomI.inducedDipoleP[0]*gux7+atomI.inducedDipoleP[1]*gux9+atomI.inducedDipoleP[2]*gux10)
                           +atomJ.inducedDipole[1]*(atomI.inducedDipoleP[0]*guy7+atomI.inducedDipoleP[1]*guy9+atomI.inducedDipoleP[2]*guy10)
                           +atomJ.inducedDipole[2]*(atomI.inducedDipoleP[0]*guz7+atomI.inducedDipoleP[1]*guz9+atomI.inducedDipoleP[2]*guz10));

        duvdr = atomI.inducedDipole[0]*(atomJ.inducedDipoleP[0]*gux22+atomJ.inducedDipoleP[1]*gux23+atomJ.inducedDipoleP[2]*gux24)
                            + atomI.inducedDipole[1]*(atomJ.inducedDipoleP[0]*guy22+atomJ.inducedDipoleP[1]*guy23+atomJ.inducedDipoleP[2]*guy24)
                            + atomI.inducedDipole[2]*(atomJ.inducedDipoleP[0]*guz22+atomJ.inducedDipoleP[1]*guz23+atomJ.inducedDipoleP[2]*guz24)
                            + atomJ.inducedDipole[0]*(atomI.inducedDipoleP[0]*gux22+atomI.inducedDipoleP[1]*gux23+atomI.inducedDipoleP[2]*gux24)
                            + atomJ.inducedDipole[1]*(atomI.inducedDipoleP[0]*guy22+atomI.inducedDipoleP[1]*guy23+atomI.inducedDipoleP[2]*guy24)
                            + atomJ.inducedDipole[2]*(atomI.inducedDipoleP[0]*guz22+atomI.inducedDipoleP[1]*guz23+atomI.inducedDipoleP[2]*guz24);
        dpbi = dpbi - 0.5f*atomJ.bornRadius*duvdr;
        dpbk = dpbk - 0.5f*atomI.bornRadius*duvdr;
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
    float trqi1   = atomI.labFrameDipole[1]*fid3 - atomI.labFrameDipole[2]*fid2;
    float trqi2   = atomI.labFrameDipole[2]*fid1 - atomI.labFrameDipole[0]*fid3;
    float trqi3   = atomI.labFrameDipole[0]*fid2 - atomI.labFrameDipole[1]*fid1;

    float trqi_k1 = atomJ.labFrameDipole[1]*fkd3 - atomJ.labFrameDipole[2]*fkd2;
    float trqi_k2 = atomJ.labFrameDipole[2]*fkd1 - atomJ.labFrameDipole[0]*fkd3;
    float trqi_k3 = atomJ.labFrameDipole[0]*fkd2 - atomJ.labFrameDipole[1]*fkd1;


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

    trqi1 += 2.0f*(atomI.labFrameQuadrupole_XY*fidg13+atomI.labFrameQuadrupole_YY*fidg23+atomI.labFrameQuadrupole_YZ*fidg33
                        -atomI.labFrameQuadrupole_XZ*fidg12-atomI.labFrameQuadrupole_YZ*fidg22-atomI.labFrameQuadrupole_ZZ*fidg32);
    trqi2 += 2.0f*(atomI.labFrameQuadrupole_XZ*fidg11+atomI.labFrameQuadrupole_YZ*fidg21+atomI.labFrameQuadrupole_ZZ*fidg31
                        -atomI.labFrameQuadrupole_XX*fidg13-atomI.labFrameQuadrupole_XY*fidg23-atomI.labFrameQuadrupole_XZ*fidg33);

    trqi3 += 2.0f*(atomI.labFrameQuadrupole_XX*fidg12+atomI.labFrameQuadrupole_XY*fidg22+atomI.labFrameQuadrupole_XZ*fidg32
                        -atomI.labFrameQuadrupole_XY*fidg11-atomI.labFrameQuadrupole_YY*fidg21-atomI.labFrameQuadrupole_YZ*fidg31);

    trqi_k1 += 2.0f*
                        (atomJ.labFrameQuadrupole_XY*fkdg13+atomJ.labFrameQuadrupole_YY*fkdg23+atomJ.labFrameQuadrupole_YZ*fkdg33
                        -atomJ.labFrameQuadrupole_XZ*fkdg12-atomJ.labFrameQuadrupole_YZ*fkdg22-atomJ.labFrameQuadrupole_ZZ*fkdg32);

    trqi_k2 += 2.0f*
                        (atomJ.labFrameQuadrupole_XZ*fkdg11+atomJ.labFrameQuadrupole_YZ*fkdg21+atomJ.labFrameQuadrupole_ZZ*fkdg31
                        -atomJ.labFrameQuadrupole_XX*fkdg13-atomJ.labFrameQuadrupole_XY*fkdg23-atomJ.labFrameQuadrupole_XZ*fkdg33);

    trqi_k3 += 2.0f*
                        (atomJ.labFrameQuadrupole_XX*fkdg12+atomJ.labFrameQuadrupole_XY*fkdg22+atomJ.labFrameQuadrupole_XZ*fkdg32
                        -atomJ.labFrameQuadrupole_XY*fkdg11-atomJ.labFrameQuadrupole_YY*fkdg21-atomJ.labFrameQuadrupole_YZ*fkdg31);

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

        totalForce    += saTerm / bornRadius;
        totalForce    *= bornRadius * bornRadius * obcChain;

        fieldOut[pos]  = totalForce;

        energy        += saTerm;
        pos           += gridDim.x * blockDim.x;
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
#ifdef AMOEBA_DEBUG
    if( amoebaGpu->log ){

        // kClearEnergy() should be called prior to kReduceToBornForcePrefactorAndSASA_kernel

        double energy = kReduceEnergy( amoebaGpu->gpuContext );
        amoebaGpu->gpuContext->psBornForce->Download();
        amoebaGpu->gpuContext->psObcData->Download();
        amoebaGpu->gpuContext->psBornRadii->Download();
        (void) fprintf( amoebaGpu->log, "Born force w/ cavity energy=%15.7e.\n", energy ); (void) fflush( amoebaGpu->log );

        for( int ii = 0; ii < amoebaGpu->gpuContext->natoms; ii++ ){
           (void) fprintf( amoebaGpu->log, "%5d ", ii);
           (void) fprintf( amoebaGpu->log,"bF %16.9e obc=%16.9e bR=%16.9e\n",
                           amoebaGpu->gpuContext->psBornForce->_pSysStream[0][ii],
                           amoebaGpu->gpuContext->psObcData->_pSysStream[0][ii].x,
                           amoebaGpu->gpuContext->psBornRadii->_pSysStream[0][ii] );
        }
        (void) fflush( amoebaGpu->log );
        if( 1 ){
            std::vector<int> fileId;
            //fileId.push_back( 0 );
            VectorOfDoubleVectors outputVector;
            cudaLoadCudaFloat4Array( amoebaGpu->gpuContext->natoms,  3, amoebaGpu->gpuContext->psPosq4,       outputVector );
            cudaLoadCudaFloatArray(  amoebaGpu->gpuContext->natoms,  1, amoebaGpu->gpuContext->psBornRadii,   outputVector );
            cudaLoadCudaFloat2Array( amoebaGpu->gpuContext->natoms,  2, amoebaGpu->gpuContext->psObcData,     outputVector );
            cudaLoadCudaFloatArray(  amoebaGpu->gpuContext->natoms,  1, amoebaGpu->gpuContext->psBornForce,   outputVector );
            cudaWriteVectorOfDoubleVectorsToFile( "CudaBornForce", fileId, outputVector );
            (void) fprintf( amoebaGpu->log, "kReduceToBornForcePrefactor: exiting.\n" );
            (void) fprintf( stderr, "kReduceToBornForcePrefactor: exiting.\n" ); (void) fflush( stderr );
            exit(0);
        }

    }
#endif

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
        (void) fprintf( amoebaGpu->log, "%s %d maxCovalentDegreeSz=%d ZZZ\n",
                        methodName, gpu->natoms, amoebaGpu->maxCovalentDegreeSz );
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
            maxThreads = 512;
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
