//-----------------------------------------------------------------------------------------

//-----------------------------------------------------------------------------------------

#include "amoebaGpuTypes.h"
#include "amoebaCudaKernels.h"
#include "kCalculateAmoebaCudaUtilities.h"
#include "kCalculateAmoebaCudaKirkwoodParticle.h"

//#define AMOEBA_DEBUG

static __constant__ cudaGmxSimulation cSim;
static __constant__ cudaAmoebaGmxSimulation cAmoebaSim;

void SetCalculateAmoebaKirkwoodEDiffSim(amoebaGpuContext amoebaGpu)
{
    cudaError_t status;
    gpuContext gpu = amoebaGpu->gpuContext;
    status         = cudaMemcpyToSymbol(cSim, &gpu->sim, sizeof(cudaGmxSimulation));    
    RTERROR(status, "SetCalculateAmoebaKirkwoodEDiffSim: cudaMemcpyToSymbol: SetSim copy to cSim failed");
    status         = cudaMemcpyToSymbol(cAmoebaSim, &amoebaGpu->amoebaSim, sizeof(cudaAmoebaGmxSimulation));    
    RTERROR(status, "SetCalculateAmoebaKirkwoodEDiffSim: cudaMemcpyToSymbol: SetSim copy to cAmoebaSim failed");
}

void GetCalculateAmoebaKirkwoodEDiffSim(amoebaGpuContext amoebaGpu)
{
    cudaError_t status;
    gpuContext gpu = amoebaGpu->gpuContext;
    status = cudaMemcpyFromSymbol(&gpu->sim, cSim, sizeof(cudaGmxSimulation));    
    RTERROR(status, "GetCalculateAmoebaKirkwoodEDiffSim: cudaMemcpyFromSymbol: SetSim copy from cSim failed");
    status = cudaMemcpyFromSymbol(&amoebaGpu->amoebaSim, cAmoebaSim, sizeof(cudaAmoebaGmxSimulation));    
    RTERROR(status, "GetCalculateAmoebaKirkwoodEDiffSim: cudaMemcpyFromSymbol: SetSim copy from cAmoebaSim failed");
}

__device__ void loadKirkwoodEDiffShared( struct KirkwoodEDiffParticle* sA, unsigned int atomI,
                                         float4* atomCoord,
                                         float* labDipole,          float* labQuadrupole,
                                         float* inducedDipole,      float* inducedDipolePolar,
                                         float* inducedDipoleS,     float* inducedDipolePolarS )
{
    // coordinates & charge

    sA->x                        = atomCoord[atomI].x;
    sA->y                        = atomCoord[atomI].y;
    sA->z                        = atomCoord[atomI].z;
    sA->q                        = atomCoord[atomI].w;

    sA->damp                     = cAmoebaSim.pDampingFactorAndThole[atomI].x;
    sA->thole                    = cAmoebaSim.pDampingFactorAndThole[atomI].y;

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

    // induced dipole

    sA->inducedDipoleS[0]        = inducedDipoleS[atomI*3];
    sA->inducedDipoleS[1]        = inducedDipoleS[atomI*3+1];
    sA->inducedDipoleS[2]        = inducedDipoleS[atomI*3+2];

    // induced dipole polar

    sA->inducedDipolePS[0]       = inducedDipolePolarS[atomI*3];
    sA->inducedDipolePS[1]       = inducedDipolePolarS[atomI*3+1];
    sA->inducedDipolePS[2]       = inducedDipolePolarS[atomI*3+2];

}

/*****************************************************************************

       ediff1 correct vacuum to SCRF derivatives

       calculates the energy and derivatives of polarizing
       the vacuum induced dipoles to their SCRF polarized values

******************************************************************************/

__device__ void calculateKirkwoodEDiffPairIxn_kernel( KirkwoodEDiffParticle& atomI,  KirkwoodEDiffParticle& atomJ,
                                                      float pscale,                  float dscale,
                                                      float*  outputEnergy,          float*  outputForce,
                                                      float*  outputTorqueI,         float* outputTorqueJ
#ifdef AMOEBA_DEBUG
                                                     , float4* debugArray 
#endif
){

    //float f;
    //float gfd;
    float scale3,scale5;
    float scale7;
    //float scale9;
    float psc3,psc5,psc7;
    //float psc9;
    float dsc3,dsc5,dsc7;
    //float dsc9;
    float scale3i,scale5i;
    //float scale7i;
    float r,rr1,rr3;
    float rr5,rr7,rr9;
    float pgamma;
    const float uscale = 1.0f;

//      float ftm1i(4,maxatm);
//      float ttm1i(4,maxatm);
//      float trqi(4,maxatm);
 
    // set conversion factor, cutoff and scaling coefficients

    //f           = cAmoebaSim.electric / cAmoebaSim.dwater;

    // deltaR

    float xr          = atomJ.x - atomI.x;
    float yr          = atomJ.y - atomI.y;
    float zr          = atomJ.z - atomI.z;

    float r2          = xr*xr + yr*yr + zr*zr;
    if( r2 > cAmoebaSim.scalingDistanceCutoff ){
    }

    r           = sqrtf(r2);
    rr1         = 1.0f / r;
    rr3         = rr1 / r2;
    rr5         = 3.0f * rr3 / r2;
    rr7         = 5.0f * rr5 / r2;
    rr9         = 7.0f * rr7 / r2;

    scale3      = 1.0f;
    scale5      = 1.0f;
    scale7      = 1.0f;
    //scale9      = 1.0f;

    float ddsc3_1    = 0.0f;
    float ddsc3_2    = 0.0f;
    float ddsc3_3    = 0.0f;

    float ddsc5_1    = 0.0f;
    float ddsc5_2    = 0.0f;
    float ddsc5_3    = 0.0f;

    float ddsc7_1    = 0.0f;
    float ddsc7_2    = 0.0f;
    float ddsc7_3    = 0.0f;

    // apply Thole polarization damping to scale factors
 
    float damp = atomI.damp * atomJ.damp;
    if( damp != 0.0f ){
        pgamma          = atomJ.thole > atomI.thole ? atomI.thole : atomJ.thole;
        float ratio     = (r/damp);
        damp            = -pgamma * ratio*ratio*ratio;
        if( damp > -50.0f){
            float dampE  = expf( damp );
            float damp2  = damp*damp;
            scale3       = 1.0f - dampE;
            scale5       = 1.0f - (1.0f - damp)*dampE;
            scale7       = 1.0f - (1.0f - damp + 0.6f*damp2)*dampE;
            //scale9       = 1.0f - (1.0f - damp + (18.0f*damp2 - (9.0f*damp*damp2))/35.0f)*dampE;

            ddsc3_1     = -3.0f*damp*exp(damp) * xr/r2;
            ddsc3_2     = -3.0f*damp*exp(damp) * yr/r2;
            ddsc3_3     = -3.0f*damp*exp(damp) * zr/r2;

            ddsc5_1     = -damp * ddsc3_1;
            ddsc5_2     = -damp * ddsc3_2;
            ddsc5_3     = -damp * ddsc3_3;

            ddsc7_1     = (-0.2f-0.6f*damp) * ddsc5_1;
            ddsc7_2     = (-0.2f-0.6f*damp) * ddsc5_2;
            ddsc7_3     = (-0.2f-0.6f*damp) * ddsc5_3;
        }
    }

    scale3i             = scale3 * uscale;
    scale5i             = scale5 * uscale;
    //scale7i             = scale7 * uscale;

    dsc3                = scale3 * dscale;
    dsc5                = scale5 * dscale;
    dsc7                = scale7 * dscale;
    //dsc9                = scale9 * dscale;

    psc3                = scale3 * pscale;
    psc5                = scale5 * pscale;
    psc7                = scale7 * pscale;
    //psc9                = scale9 * pscale;
 
    // construct auxiliary vectors for permanent terms
 
#if 0
    float dixdk1            = atomI.labFrameDipole[1]*atomJ.labFrameDipole[2] - atomI.labFrameDipole[2]*atomJ.labFrameDipole[1];
    float dixdk2            = atomI.labFrameDipole[2]*atomJ.labFrameDipole[0] - atomI.labFrameDipole[0]*atomJ.labFrameDipole[2];
    float dixdk3            = atomI.labFrameDipole[0]*atomJ.labFrameDipole[1] - atomI.labFrameDipole[1]*atomJ.labFrameDipole[0];
#endif

    float dixr1             = atomI.labFrameDipole[1]*zr - atomI.labFrameDipole[2]*yr;
    float dixr2             = atomI.labFrameDipole[2]*xr - atomI.labFrameDipole[0]*zr;
    float dixr3             = atomI.labFrameDipole[0]*yr - atomI.labFrameDipole[1]*xr;

    float dkxr1             = atomJ.labFrameDipole[1]*zr - atomJ.labFrameDipole[2]*yr;
    float dkxr2             = atomJ.labFrameDipole[2]*xr - atomJ.labFrameDipole[0]*zr;
    float dkxr3             = atomJ.labFrameDipole[0]*yr - atomJ.labFrameDipole[1]*xr;

    float qir1              = atomI.labFrameQuadrupole_XX*xr + atomI.labFrameQuadrupole_XY*yr + atomI.labFrameQuadrupole_XZ*zr;
    float qir2              = atomI.labFrameQuadrupole_XY*xr + atomI.labFrameQuadrupole_YY*yr + atomI.labFrameQuadrupole_YZ*zr;
    float qir3              = atomI.labFrameQuadrupole_XZ*xr + atomI.labFrameQuadrupole_YZ*yr + atomI.labFrameQuadrupole_ZZ*zr;

    float qkr1              = atomJ.labFrameQuadrupole_XX*xr + atomJ.labFrameQuadrupole_XY*yr + atomJ.labFrameQuadrupole_XZ*zr;
    float qkr2              = atomJ.labFrameQuadrupole_XY*xr + atomJ.labFrameQuadrupole_YY*yr + atomJ.labFrameQuadrupole_YZ*zr;
    float qkr3              = atomJ.labFrameQuadrupole_XZ*xr + atomJ.labFrameQuadrupole_YZ*yr + atomJ.labFrameQuadrupole_ZZ*zr;

#if 0
    float qiqkr1            = atomI.labFrameQuadrupole_XX*qkr1 + atomI.labFrameQuadrupole_XY*qkr2 + atomI.labFrameQuadrupole_XZ*qkr3;
    float qiqkr2            = atomI.labFrameQuadrupole_XY*qkr1 + atomI.labFrameQuadrupole_YY*qkr2 + atomI.labFrameQuadrupole_YZ*qkr3;
    float qiqkr3            = atomI.labFrameQuadrupole_XZ*qkr1 + atomI.labFrameQuadrupole_YZ*qkr2 + atomI.labFrameQuadrupole_ZZ*qkr3;

    float qkqir1            = atomJ.labFrameQuadrupole_XX*qir1 + atomJ.labFrameQuadrupole_XY*qir2 + atomJ.labFrameQuadrupole_XZ*qir3;
    float qkqir2            = atomJ.labFrameQuadrupole_XY*qir1 + atomJ.labFrameQuadrupole_YY*qir2 + atomJ.labFrameQuadrupole_YZ*qir3;
    float qkqir3            = atomJ.labFrameQuadrupole_XZ*qir1 + atomJ.labFrameQuadrupole_YZ*qir2 + atomJ.labFrameQuadrupole_ZZ*qir3;

    float qixqk1            = atomI.labFrameQuadrupole_XY*atomJ.labFrameQuadrupole_XZ + atomI.labFrameQuadrupole_YY*atomJ.labFrameQuadrupole_YZ + atomI.labFrameQuadrupole_YZ*atomJ.labFrameQuadrupole_ZZ -
                              atomI.labFrameQuadrupole_XZ*atomJ.labFrameQuadrupole_XY - atomI.labFrameQuadrupole_YZ*atomJ.labFrameQuadrupole_YY - atomI.labFrameQuadrupole_ZZ*atomJ.labFrameQuadrupole_YZ;

    float qixqk2            = atomI.labFrameQuadrupole_XZ*atomJ.labFrameQuadrupole_XX + atomI.labFrameQuadrupole_YZ*atomJ.labFrameQuadrupole_XY + atomI.labFrameQuadrupole_ZZ*atomJ.labFrameQuadrupole_XZ -
                              atomI.labFrameQuadrupole_XX*atomJ.labFrameQuadrupole_XZ - atomI.labFrameQuadrupole_XY*atomJ.labFrameQuadrupole_YZ - atomI.labFrameQuadrupole_XZ*atomJ.labFrameQuadrupole_ZZ;

    float qixqk3            = atomI.labFrameQuadrupole_XX*atomJ.labFrameQuadrupole_XY + atomI.labFrameQuadrupole_XY*atomJ.labFrameQuadrupole_YY + atomI.labFrameQuadrupole_XZ*atomJ.labFrameQuadrupole_YZ -
                              atomI.labFrameQuadrupole_XY*atomJ.labFrameQuadrupole_XX - atomI.labFrameQuadrupole_YY*atomJ.labFrameQuadrupole_XY - atomI.labFrameQuadrupole_YZ*atomJ.labFrameQuadrupole_XZ;
#endif

    float rxqir1            = yr*qir3 - zr*qir2;
    float rxqir2            = zr*qir1 - xr*qir3;
    float rxqir3            = xr*qir2 - yr*qir1;

    float rxqkr1            = yr*qkr3 - zr*qkr2;
    float rxqkr2            = zr*qkr1 - xr*qkr3;
    float rxqkr3            = xr*qkr2 - yr*qkr1;

#if 0
    float rxqikr1           = yr*qiqkr3 - zr*qiqkr2;
    float rxqikr2           = zr*qiqkr1 - xr*qiqkr3;
    float rxqikr3           = xr*qiqkr2 - yr*qiqkr1;

    float rxqkir1           = yr*qkqir3 - zr*qkqir2;
    float rxqkir2           = zr*qkqir1 - xr*qkqir3;
    float rxqkir3           = xr*qkqir2 - yr*qkqir1;

    float qkrxqir1          = qkr2*qir3 - qkr3*qir2;
    float qkrxqir2          = qkr3*qir1 - qkr1*qir3;
    float qkrxqir3          = qkr1*qir2 - qkr2*qir1;

    float qidk1             = atomI.labFrameQuadrupole_XX*atomJ.labFrameDipole[0] + atomI.labFrameQuadrupole_XY*atomJ.labFrameDipole[1] + atomI.labFrameQuadrupole_XZ*atomJ.labFrameDipole[2];
    float qidk2             = atomI.labFrameQuadrupole_XY*atomJ.labFrameDipole[0] + atomI.labFrameQuadrupole_YY*atomJ.labFrameDipole[1] + atomI.labFrameQuadrupole_YZ*atomJ.labFrameDipole[2];
    float qidk3             = atomI.labFrameQuadrupole_XZ*atomJ.labFrameDipole[0] + atomI.labFrameQuadrupole_YZ*atomJ.labFrameDipole[1] + atomI.labFrameQuadrupole_ZZ*atomJ.labFrameDipole[2];

    float qkdi1             = atomJ.labFrameQuadrupole_XX*atomI.labFrameDipole[0] + atomJ.labFrameQuadrupole_XY*atomI.labFrameDipole[1] + atomJ.labFrameQuadrupole_XZ*atomI.labFrameDipole[2];
    float qkdi2             = atomJ.labFrameQuadrupole_XY*atomI.labFrameDipole[0] + atomJ.labFrameQuadrupole_YY*atomI.labFrameDipole[1] + atomJ.labFrameQuadrupole_YZ*atomI.labFrameDipole[2];
    float qkdi3             = atomJ.labFrameQuadrupole_XZ*atomI.labFrameDipole[0] + atomJ.labFrameQuadrupole_YZ*atomI.labFrameDipole[1] + atomJ.labFrameQuadrupole_ZZ*atomI.labFrameDipole[2];

    float dixqkr1           = atomI.labFrameDipole[1]*qkr3 - atomI.labFrameDipole[2]*qkr2;
    float dixqkr2           = atomI.labFrameDipole[2]*qkr1 - atomI.labFrameDipole[0]*qkr3;
    float dixqkr3           = atomI.labFrameDipole[0]*qkr2 - atomI.labFrameDipole[1]*qkr1;

    float dkxqir1           = atomJ.labFrameDipole[1]*qir3 - atomJ.labFrameDipole[2]*qir2;
    float dkxqir2           = atomJ.labFrameDipole[2]*qir1 - atomJ.labFrameDipole[0]*qir3;
    float dkxqir3           = atomJ.labFrameDipole[0]*qir2 - atomJ.labFrameDipole[1]*qir1;

    float rxqidk1           = yr*qidk3 - zr*qidk2;
    float rxqidk2           = zr*qidk1 - xr*qidk3;
    float rxqidk3           = xr*qidk2 - yr*qidk1;

    float rxqkdi1           = yr*qkdi3 - zr*qkdi2;
    float rxqkdi2           = zr*qkdi1 - xr*qkdi3;
    float rxqkdi3           = xr*qkdi2 - yr*qkdi1;
#endif
 
    // get intermediate variables for permanent energy terms
 
    float sc3               = atomI.labFrameDipole[0]*xr  + atomI.labFrameDipole[1]*yr  + atomI.labFrameDipole[2]*zr;
    float sc4               = atomJ.labFrameDipole[0]*xr  + atomJ.labFrameDipole[1]*yr  + atomJ.labFrameDipole[2]*zr;
    float sc5               = qir1*xr + qir2*yr + qir3*zr;
    float sc6               = qkr1*xr + qkr2*yr + qkr3*zr;
 
    // construct auxiliary vectors for induced terms
 
    float dixuk1            = atomI.labFrameDipole[1]*atomJ.inducedDipoleS[2] - atomI.labFrameDipole[2]*atomJ.inducedDipoleS[1];
    float dixuk2            = atomI.labFrameDipole[2]*atomJ.inducedDipoleS[0] - atomI.labFrameDipole[0]*atomJ.inducedDipoleS[2];
    float dixuk3            = atomI.labFrameDipole[0]*atomJ.inducedDipoleS[1] - atomI.labFrameDipole[1]*atomJ.inducedDipoleS[0];

    float dkxui1            = atomJ.labFrameDipole[1]*atomI.inducedDipoleS[2] - atomJ.labFrameDipole[2]*atomI.inducedDipoleS[1];
    float dkxui2            = atomJ.labFrameDipole[2]*atomI.inducedDipoleS[0] - atomJ.labFrameDipole[0]*atomI.inducedDipoleS[2];
    float dkxui3            = atomJ.labFrameDipole[0]*atomI.inducedDipoleS[1] - atomJ.labFrameDipole[1]*atomI.inducedDipoleS[0];

    float dixukp1           = atomI.labFrameDipole[1]*atomJ.inducedDipolePS[2] - atomI.labFrameDipole[2]*atomJ.inducedDipolePS[1];
    float dixukp2           = atomI.labFrameDipole[2]*atomJ.inducedDipolePS[0] - atomI.labFrameDipole[0]*atomJ.inducedDipolePS[2];
    float dixukp3           = atomI.labFrameDipole[0]*atomJ.inducedDipolePS[1] - atomI.labFrameDipole[1]*atomJ.inducedDipolePS[0];

    float dkxuip1           = atomJ.labFrameDipole[1]*atomI.inducedDipolePS[2] - atomJ.labFrameDipole[2]*atomI.inducedDipolePS[1];
    float dkxuip2           = atomJ.labFrameDipole[2]*atomI.inducedDipolePS[0] - atomJ.labFrameDipole[0]*atomI.inducedDipolePS[2];
    float dkxuip3           = atomJ.labFrameDipole[0]*atomI.inducedDipolePS[1] - atomJ.labFrameDipole[1]*atomI.inducedDipolePS[0];

    float qiuk1             = atomI.labFrameQuadrupole_XX*atomJ.inducedDipoleS[0] + atomI.labFrameQuadrupole_XY*atomJ.inducedDipoleS[1] + atomI.labFrameQuadrupole_XZ*atomJ.inducedDipoleS[2];
    float qiuk2             = atomI.labFrameQuadrupole_XY*atomJ.inducedDipoleS[0] + atomI.labFrameQuadrupole_YY*atomJ.inducedDipoleS[1] + atomI.labFrameQuadrupole_YZ*atomJ.inducedDipoleS[2];
    float qiuk3             = atomI.labFrameQuadrupole_XZ*atomJ.inducedDipoleS[0] + atomI.labFrameQuadrupole_YZ*atomJ.inducedDipoleS[1] + atomI.labFrameQuadrupole_ZZ*atomJ.inducedDipoleS[2];

    float qkui1             = atomJ.labFrameQuadrupole_XX*atomI.inducedDipoleS[0] + atomJ.labFrameQuadrupole_XY*atomI.inducedDipoleS[1] + atomJ.labFrameQuadrupole_XZ*atomI.inducedDipoleS[2];
    float qkui2             = atomJ.labFrameQuadrupole_XY*atomI.inducedDipoleS[0] + atomJ.labFrameQuadrupole_YY*atomI.inducedDipoleS[1] + atomJ.labFrameQuadrupole_YZ*atomI.inducedDipoleS[2];
    float qkui3             = atomJ.labFrameQuadrupole_XZ*atomI.inducedDipoleS[0] + atomJ.labFrameQuadrupole_YZ*atomI.inducedDipoleS[1] + atomJ.labFrameQuadrupole_ZZ*atomI.inducedDipoleS[2];

    float qiukp1            = atomI.labFrameQuadrupole_XX*atomJ.inducedDipolePS[0] + atomI.labFrameQuadrupole_XY*atomJ.inducedDipolePS[1] + atomI.labFrameQuadrupole_XZ*atomJ.inducedDipolePS[2];
    float qiukp2            = atomI.labFrameQuadrupole_XY*atomJ.inducedDipolePS[0] + atomI.labFrameQuadrupole_YY*atomJ.inducedDipolePS[1] + atomI.labFrameQuadrupole_YZ*atomJ.inducedDipolePS[2];
    float qiukp3            = atomI.labFrameQuadrupole_XZ*atomJ.inducedDipolePS[0] + atomI.labFrameQuadrupole_YZ*atomJ.inducedDipolePS[1] + atomI.labFrameQuadrupole_ZZ*atomJ.inducedDipolePS[2];

    float qkuip1            = atomJ.labFrameQuadrupole_XX*atomI.inducedDipolePS[0] + atomJ.labFrameQuadrupole_XY*atomI.inducedDipolePS[1] + atomJ.labFrameQuadrupole_XZ*atomI.inducedDipolePS[2];
    float qkuip2            = atomJ.labFrameQuadrupole_XY*atomI.inducedDipolePS[0] + atomJ.labFrameQuadrupole_YY*atomI.inducedDipolePS[1] + atomJ.labFrameQuadrupole_YZ*atomI.inducedDipolePS[2];
    float qkuip3            = atomJ.labFrameQuadrupole_XZ*atomI.inducedDipolePS[0] + atomJ.labFrameQuadrupole_YZ*atomI.inducedDipolePS[1] + atomJ.labFrameQuadrupole_ZZ*atomI.inducedDipolePS[2];

    float uixqkr1           = atomI.inducedDipoleS[1]*qkr3 - atomI.inducedDipoleS[2]*qkr2;
    float uixqkr2           = atomI.inducedDipoleS[2]*qkr1 - atomI.inducedDipoleS[0]*qkr3;
    float uixqkr3           = atomI.inducedDipoleS[0]*qkr2 - atomI.inducedDipoleS[1]*qkr1;

    float ukxqir1           = atomJ.inducedDipoleS[1]*qir3 - atomJ.inducedDipoleS[2]*qir2;
    float ukxqir2           = atomJ.inducedDipoleS[2]*qir1 - atomJ.inducedDipoleS[0]*qir3;
    float ukxqir3           = atomJ.inducedDipoleS[0]*qir2 - atomJ.inducedDipoleS[1]*qir1;

    float uixqkrp1          = atomI.inducedDipolePS[1]*qkr3 - atomI.inducedDipolePS[2]*qkr2;
    float uixqkrp2          = atomI.inducedDipolePS[2]*qkr1 - atomI.inducedDipolePS[0]*qkr3;
    float uixqkrp3          = atomI.inducedDipolePS[0]*qkr2 - atomI.inducedDipolePS[1]*qkr1;

    float ukxqirp1          = atomJ.inducedDipolePS[1]*qir3 - atomJ.inducedDipolePS[2]*qir2;
    float ukxqirp2          = atomJ.inducedDipolePS[2]*qir1 - atomJ.inducedDipolePS[0]*qir3;
    float ukxqirp3          = atomJ.inducedDipolePS[0]*qir2 - atomJ.inducedDipolePS[1]*qir1;

    float rxqiuk1           = yr*qiuk3 - zr*qiuk2;
    float rxqiuk2           = zr*qiuk1 - xr*qiuk3;
    float rxqiuk3           = xr*qiuk2 - yr*qiuk1;

    float rxqkui1           = yr*qkui3 - zr*qkui2;
    float rxqkui2           = zr*qkui1 - xr*qkui3;
    float rxqkui3           = xr*qkui2 - yr*qkui1;

    float rxqiukp1          = yr*qiukp3 - zr*qiukp2;
    float rxqiukp2          = zr*qiukp1 - xr*qiukp3;
    float rxqiukp3          = xr*qiukp2 - yr*qiukp1;

    float rxqkuip1          = yr*qkuip3 - zr*qkuip2;
    float rxqkuip2          = zr*qkuip1 - xr*qkuip3;
    float rxqkuip3          = xr*qkuip2 - yr*qkuip1;
 
    // get intermediate variables for induction energy terms

    float sci1              = atomI.inducedDipoleS[0]*atomJ.labFrameDipole[0] + atomI.inducedDipoleS[1]*atomJ.labFrameDipole[1] +
                              atomI.inducedDipoleS[2]*atomJ.labFrameDipole[2] + atomI.labFrameDipole[0]*atomJ.inducedDipoleS[0] +
                              atomI.labFrameDipole[1]*atomJ.inducedDipoleS[1] + atomI.labFrameDipole[2]*atomJ.inducedDipoleS[2];

    float sci3              = atomI.inducedDipoleS[0]*xr + atomI.inducedDipoleS[1]*yr + atomI.inducedDipoleS[2]*zr;
    float sci4              = atomJ.inducedDipoleS[0]*xr + atomJ.inducedDipoleS[1]*yr + atomJ.inducedDipoleS[2]*zr;

    float sci7              = qir1*atomJ.inducedDipoleS[0] + qir2*atomJ.inducedDipoleS[1] + qir3*atomJ.inducedDipoleS[2];
    float sci8              = qkr1*atomI.inducedDipoleS[0] + qkr2*atomI.inducedDipoleS[1] + qkr3*atomI.inducedDipoleS[2];

    float scip1             = atomI.inducedDipolePS[0]*atomJ.labFrameDipole[0] + atomI.inducedDipolePS[1]*atomJ.labFrameDipole[1] +
                              atomI.inducedDipolePS[2]*atomJ.labFrameDipole[2] + atomI.labFrameDipole[0]*atomJ.inducedDipolePS[0] +
                              atomI.labFrameDipole[1]*atomJ.inducedDipolePS[1] + atomI.labFrameDipole[2]*atomJ.inducedDipolePS[2];

    float scip2             = atomI.inducedDipoleS[0]*atomJ.inducedDipolePS[0] + atomI.inducedDipoleS[1]*atomJ.inducedDipolePS[1] +
                              atomI.inducedDipoleS[2]*atomJ.inducedDipolePS[2] + atomI.inducedDipolePS[0]*atomJ.inducedDipoleS[0] +
                              atomI.inducedDipolePS[1]*atomJ.inducedDipoleS[1] + atomI.inducedDipolePS[2]*atomJ.inducedDipoleS[2];

    float scip3             = atomI.inducedDipolePS[0]*xr + atomI.inducedDipolePS[1]*yr + atomI.inducedDipolePS[2]*zr;
    float scip4             = atomJ.inducedDipolePS[0]*xr + atomJ.inducedDipolePS[1]*yr + atomJ.inducedDipolePS[2]*zr;

    float scip7             = qir1*atomJ.inducedDipolePS[0] + qir2*atomJ.inducedDipolePS[1] + qir3*atomJ.inducedDipolePS[2];
    float scip8             = qkr1*atomI.inducedDipolePS[0] + qkr2*atomI.inducedDipolePS[1] + qkr3*atomI.inducedDipolePS[2];

 
    // calculate the gl functions for potential energy

    float gli1              = atomJ.q*sci3 - atomI.q*sci4;
    float gli2              = -sc3*sci4 - sci3*sc4;
    float gli3              = sci3*sc6 - sci4*sc5;
    float gli6              = sci1;
    float gli7              = 2.0f * (sci7-sci8);
    float glip1             = atomJ.q*scip3 - atomI.q*scip4;
    float glip2             = -sc3*scip4 - scip3*sc4;
    float glip3             = scip3*sc6 - scip4*sc5;
    float glip6             = scip1;
    float glip7             = 2.0f * (scip7-scip8);

    // get the permanent multipole and induced energies;

    *outputEnergy       = 0.5f * (rr3*(gli1+gli6)*psc3 +
                                  rr5*(gli2+gli7)*psc5 +
                                  rr7*gli3*psc7);

    // intermediate variables for the induced-permanent terms

    float gfi1              = 0.5f*rr5*((gli1+gli6)*psc3 +
                                        (glip1+glip6)*dsc3+scip2*scale3i) +
                              0.5f*rr7*((gli7+gli2)*psc5 +
                                        (glip7+glip2)*dsc5 -
                              (sci3*scip4+scip3*sci4)*scale5i) +
                              0.5f*rr9*(gli3*psc7+glip3*dsc7);

    float gfi4              = 2.0f * rr5;
    float gfi5              = rr7 * (sci4*psc7+scip4*dsc7);
    float gfi6              = -rr7 * (sci3*psc7+scip3*dsc7);

    // get the induced force;
 
    float ftm2i1            = gfi1*xr + 0.5f*
                          (- rr3*atomJ.q*(atomI.inducedDipoleS[0]*psc3+atomI.inducedDipolePS[0]*dsc3)
                           + rr5*sc4*(atomI.inducedDipoleS[0]*psc5+atomI.inducedDipolePS[0]*dsc5)
                           - rr7*sc6*(atomI.inducedDipoleS[0]*psc7+atomI.inducedDipolePS[0]*dsc7))
                           +(rr3*atomI.q*(atomJ.inducedDipoleS[0]*psc3+atomJ.inducedDipolePS[0]*dsc3)
                           + rr5*sc3*(atomJ.inducedDipoleS[0]*psc5+atomJ.inducedDipolePS[0]*dsc5)
                           + rr7*sc5*(atomJ.inducedDipoleS[0]*psc7+atomJ.inducedDipolePS[0]*dsc7))*0.5f
                           + rr5*scale5i*(sci4*atomI.inducedDipolePS[0]+scip4*atomI.inducedDipoleS[0]
                           + sci3*atomJ.inducedDipolePS[0]+scip3*atomJ.inducedDipoleS[0])*0.5f
                           + 0.5f*(sci4*psc5+scip4*dsc5)*rr5*atomI.labFrameDipole[0]
                           + 0.5f*(sci3*psc5+scip3*dsc5)*rr5*atomJ.labFrameDipole[0]
                           + 0.5f*gfi4*((qkui1-qiuk1)*psc5
                           + (qkuip1-qiukp1)*dsc5)
                           + gfi5*qir1 + gfi6*qkr1;
 
    float ftm2i2            = gfi1*yr + 0.5f*
                          (- rr3*atomJ.q*(atomI.inducedDipoleS[1]*psc3+atomI.inducedDipolePS[1]*dsc3)
                           + rr5*sc4*(atomI.inducedDipoleS[1]*psc5+atomI.inducedDipolePS[1]*dsc5)
                           - rr7*sc6*(atomI.inducedDipoleS[1]*psc7+atomI.inducedDipolePS[1]*dsc7))
                           +(rr3*atomI.q*(atomJ.inducedDipoleS[1]*psc3+atomJ.inducedDipolePS[1]*dsc3)
                           + rr5*sc3*(atomJ.inducedDipoleS[1]*psc5+atomJ.inducedDipolePS[1]*dsc5)
                           + rr7*sc5*(atomJ.inducedDipoleS[1]*psc7+atomJ.inducedDipolePS[1]*dsc7))*0.5f
                           + rr5*scale5i*(sci4*atomI.inducedDipolePS[1]+scip4*atomI.inducedDipoleS[1]
                           + sci3*atomJ.inducedDipolePS[1]+scip3*atomJ.inducedDipoleS[1])*0.5f
                           + 0.5f*(sci4*psc5+scip4*dsc5)*rr5*atomI.labFrameDipole[1]
                           + 0.5f*(sci3*psc5+scip3*dsc5)*rr5*atomJ.labFrameDipole[1]
                           + 0.5f*gfi4*((qkui2-qiuk2)*psc5
                           + (qkuip2-qiukp2)*dsc5)
                           + gfi5*qir2 + gfi6*qkr2;

    float ftm2i3            = gfi1*zr  + 0.5f*
                          (- rr3*atomJ.q*(atomI.inducedDipoleS[2]*psc3+atomI.inducedDipolePS[2]*dsc3)
                           + rr5*sc4*(atomI.inducedDipoleS[2]*psc5+atomI.inducedDipolePS[2]*dsc5)
                           - rr7*sc6*(atomI.inducedDipoleS[2]*psc7+atomI.inducedDipolePS[2]*dsc7))
                           +(rr3*atomI.q*(atomJ.inducedDipoleS[2]*psc3+atomJ.inducedDipolePS[2]*dsc3)
                           + rr5*sc3*(atomJ.inducedDipoleS[2]*psc5+atomJ.inducedDipolePS[2]*dsc5)
                           + rr7*sc5*(atomJ.inducedDipoleS[2]*psc7+atomJ.inducedDipolePS[2]*dsc7))*0.5f
                           + rr5*scale5i*(sci4*atomI.inducedDipolePS[2]+scip4*atomI.inducedDipoleS[2]
                           + sci3*atomJ.inducedDipolePS[2]+scip3*atomJ.inducedDipoleS[2])*0.5f
                           + 0.5f*(sci4*psc5+scip4*dsc5)*rr5*atomI.labFrameDipole[2]
                           + 0.5f*(sci3*psc5+scip3*dsc5)*rr5*atomJ.labFrameDipole[2]
                           + 0.5f*gfi4*((qkui3-qiuk3)*psc5
                           + (qkuip3-qiukp3)*dsc5)
                           + gfi5*qir3 + gfi6*qkr3;
 
    // intermediate values needed for partially excluded interactions

    float fridmp1           = 0.5f * (rr3*((gli1+gli6)*pscale
                          +(glip1+glip6)*dscale)*ddsc3_1
                          + rr5*((gli2+gli7)*pscale
                          +(glip2+glip7)*dscale)*ddsc5_1
                          + rr7*(gli3*pscale+glip3*dscale)*ddsc7_1);

    float fridmp2           = 0.5f * (rr3*((gli1+gli6)*pscale
                          +(glip1+glip6)*dscale)*ddsc3_2
                          + rr5*((gli2+gli7)*pscale
                          +(glip2+glip7)*dscale)*ddsc5_2
                          + rr7*(gli3*pscale+glip3*dscale)*ddsc7_2);

    float fridmp3           = 0.5f * (rr3*((gli1+gli6)*pscale
                          +(glip1+glip6)*dscale)*ddsc3_3
                          + rr5*((gli2+gli7)*pscale
                          +(glip2+glip7)*dscale)*ddsc5_3
                          + rr7*(gli3*pscale+glip3*dscale)*ddsc7_3);

    // get the induced-induced derivative terms
 
    float findmp1           = 0.5f * uscale * (scip2*rr3*ddsc3_1
                          - rr5*ddsc5_1*(sci3*scip4+scip3*sci4));

    float findmp2           = 0.5f * uscale * (scip2*rr3*ddsc3_2
                          - rr5*ddsc5_2*(sci3*scip4+scip3*sci4));

    float findmp3           = 0.5f * uscale * (scip2*rr3*ddsc3_3
                          - rr5*ddsc5_3*(sci3*scip4+scip3*sci4));

    // handle of scaling for partially excluded interactions

    ftm2i1           -= fridmp1 + findmp1;
    ftm2i2           -= fridmp2 + findmp2;
    ftm2i3           -= fridmp3 + findmp3;

    // correction to convert mutual to direct polarization force

#if 0
               if (poltyp .eq. 'DIRECT') then;
                  gfd             = 0.5f * (rr5*scip2*scale3i;
     &                  - rr7*(scip3*sci4+sci3*scip4)*scale5i);
                  fdir1             = gfd*xr + 0.5f*rr5*scale5i;
     &                         * (sci4*atomI.inducedDipolePS[0]+scip4*atomI.inducedDipoleS[0];
     &                           +sci3*atomJ.inducedDipolePS[0]+scip3*atomJ.inducedDipoleS[0]);
                  fdir2             = gfd*yr + 0.5f*rr5*scale5i;
     &                         * (sci4*atomI.inducedDipolePS[1]+scip4*atomI.inducedDipoleS[1];
     &                           +sci3*atomJ.inducedDipolePS[1]+scip3*atomJ.inducedDipoleS[1]);
                  fdir3             = gfd*zr + 0.5f*rr5*scale5i;
     &                         * (sci4*atomI.inducedDipolePS[2]+scip4*atomI.inducedDipoleS[2];
     &                           +sci3*atomJ.inducedDipolePS[2]+scip3*atomJ.inducedDipoleS[2]);
                  ftm2i1             = ftm2i1 - fdir1 + findmp1;
                  ftm2i2             = ftm2i2 - fdir2 + findmp2;
                  ftm2i3             = ftm2i3 - fdir3 + findmp3;
               end if;
#endif

    // now perform the torque calculation
    // intermediate terms for torque between multipoles i and k
 
    float gti2              = 0.5f * (sci4*psc5+scip4*dsc5) * rr5;
    float gti3              = 0.5f * (sci3*psc5+scip3*dsc5) * rr5;
    float gti4              = gfi4;
    float gti5              = gfi5;
    float gti6              = gfi6;

    // calculate the induced torque components

    float ttm2i1            = -rr3*(dixuk1*psc3+dixukp1*dsc3)*0.5f
                          + gti2*dixr1 + gti4*((ukxqir1+rxqiuk1)*psc5
                          +(ukxqirp1+rxqiukp1)*dsc5)*0.5f - gti5*rxqir1;

    float ttm2i2            = -rr3*(dixuk2*psc3+dixukp2*dsc3)*0.5f
                          + gti2*dixr2 + gti4*((ukxqir2+rxqiuk2)*psc5
                          +(ukxqirp2+rxqiukp2)*dsc5)*0.5f - gti5*rxqir2;

    float ttm2i3            = -rr3*(dixuk3*psc3+dixukp3*dsc3)*0.5f
                          + gti2*dixr3 + gti4*((ukxqir3+rxqiuk3)*psc5
                          +(ukxqirp3+rxqiukp3)*dsc5)*0.5f - gti5*rxqir3;

    float ttm3i1            = -rr3*(dkxui1*psc3+dkxuip1*dsc3)*0.5f
                          + gti3*dkxr1 - gti4*((uixqkr1+rxqkui1)*psc5
                          +(uixqkrp1+rxqkuip1)*dsc5)*0.5f - gti6*rxqkr1;

    float ttm3i2            = -rr3*(dkxui2*psc3+dkxuip2*dsc3)*0.5f
                          + gti3*dkxr2 - gti4*((uixqkr2+rxqkui2)*psc5
                          +(uixqkrp2+rxqkuip2)*dsc5)*0.5f - gti6*rxqkr2;

    float ttm3i3            = -rr3*(dkxui3*psc3+dkxuip3*dsc3)*0.5f
                          + gti3*dkxr3 - gti4*((uixqkr3+rxqkui3)*psc5
                          +(uixqkrp3+rxqkuip3)*dsc5)*0.5f - gti6*rxqkr3;
 
    // update force and torque on site k
#if 0
    ftm1i(1,k)          = ftm1i(1,k) + ftm2i1;
    ftm1i(2,k)          = ftm1i(2,k) + ftm2i2;
    ftm1i(3,k)          = ftm1i(3,k) + ftm2i3;

    ttm1i(1,k)          = ttm1i(1,k) + ttm3i1;
    ttm1i(2,k)          = ttm1i(2,k) + ttm3i2;
    ttm1i(3,k)          = ttm1i(3,k) + ttm3i3;

    // update force and torque on site i

    ftm1i(1,i)          = ftm1i(1,i) - ftm2i1;
    ftm1i(2,i)          = ftm1i(2,i) - ftm2i2;
    ftm1i(3,i)          = ftm1i(3,i) - ftm2i3;

    ttm1i(1,i)          = ttm1i(1,i) + ttm2i1;
    ttm1i(2,i)          = ttm1i(2,i) + ttm2i2;
    ttm1i(3,i)          = ttm1i(3,i) + ttm2i3;
#endif
    outputForce[0]      = -ftm2i1;
    outputForce[1]      = -ftm2i2;
    outputForce[2]      = -ftm2i3;

    outputTorqueI[0]    = ttm2i1;
    outputTorqueI[1]    = ttm2i2;
    outputTorqueI[2]    = ttm2i3;

    outputTorqueJ[0]    = ttm3i1;
    outputTorqueJ[1]    = ttm3i2;
    outputTorqueJ[2]    = ttm3i3;

    // construct auxiliary vectors for induced terms

    dixuk1            = atomI.labFrameDipole[1]*atomJ.inducedDipole[2] - atomI.labFrameDipole[2]*atomJ.inducedDipole[1];
    dixuk2            = atomI.labFrameDipole[2]*atomJ.inducedDipole[0] - atomI.labFrameDipole[0]*atomJ.inducedDipole[2];
    dixuk3            = atomI.labFrameDipole[0]*atomJ.inducedDipole[1] - atomI.labFrameDipole[1]*atomJ.inducedDipole[0];

    dkxui1            = atomJ.labFrameDipole[1]*atomI.inducedDipole[2] - atomJ.labFrameDipole[2]*atomI.inducedDipole[1];
    dkxui2            = atomJ.labFrameDipole[2]*atomI.inducedDipole[0] - atomJ.labFrameDipole[0]*atomI.inducedDipole[2];
    dkxui3            = atomJ.labFrameDipole[0]*atomI.inducedDipole[1] - atomJ.labFrameDipole[1]*atomI.inducedDipole[0];

    dixukp1           = atomI.labFrameDipole[1]*atomJ.inducedDipoleP[2] - atomI.labFrameDipole[2]*atomJ.inducedDipoleP[1];
    dixukp2           = atomI.labFrameDipole[2]*atomJ.inducedDipoleP[0] - atomI.labFrameDipole[0]*atomJ.inducedDipoleP[2];
    dixukp3           = atomI.labFrameDipole[0]*atomJ.inducedDipoleP[1] - atomI.labFrameDipole[1]*atomJ.inducedDipoleP[0];

    dkxuip1           = atomJ.labFrameDipole[1]*atomI.inducedDipoleP[2] - atomJ.labFrameDipole[2]*atomI.inducedDipoleP[1];
    dkxuip2           = atomJ.labFrameDipole[2]*atomI.inducedDipoleP[0] - atomJ.labFrameDipole[0]*atomI.inducedDipoleP[2];
    dkxuip3           = atomJ.labFrameDipole[0]*atomI.inducedDipoleP[1] - atomJ.labFrameDipole[1]*atomI.inducedDipoleP[0];

    qiuk1             = atomI.labFrameQuadrupole_XX*atomJ.inducedDipole[0] + atomI.labFrameQuadrupole_XY*atomJ.inducedDipole[1] + atomI.labFrameQuadrupole_XZ*atomJ.inducedDipole[2];
    qiuk2             = atomI.labFrameQuadrupole_XY*atomJ.inducedDipole[0] + atomI.labFrameQuadrupole_YY*atomJ.inducedDipole[1] + atomI.labFrameQuadrupole_YZ*atomJ.inducedDipole[2];
    qiuk3             = atomI.labFrameQuadrupole_XZ*atomJ.inducedDipole[0] + atomI.labFrameQuadrupole_YZ*atomJ.inducedDipole[1] + atomI.labFrameQuadrupole_ZZ*atomJ.inducedDipole[2];

    qkui1             = atomJ.labFrameQuadrupole_XX*atomI.inducedDipole[0] + atomJ.labFrameQuadrupole_XY*atomI.inducedDipole[1] + atomJ.labFrameQuadrupole_XZ*atomI.inducedDipole[2];
    qkui2             = atomJ.labFrameQuadrupole_XY*atomI.inducedDipole[0] + atomJ.labFrameQuadrupole_YY*atomI.inducedDipole[1] + atomJ.labFrameQuadrupole_YZ*atomI.inducedDipole[2];
    qkui3             = atomJ.labFrameQuadrupole_XZ*atomI.inducedDipole[0] + atomJ.labFrameQuadrupole_YZ*atomI.inducedDipole[1] + atomJ.labFrameQuadrupole_ZZ*atomI.inducedDipole[2];

    qiukp1            = atomI.labFrameQuadrupole_XX*atomJ.inducedDipoleP[0] + atomI.labFrameQuadrupole_XY*atomJ.inducedDipoleP[1] + atomI.labFrameQuadrupole_XZ*atomJ.inducedDipoleP[2];
    qiukp2            = atomI.labFrameQuadrupole_XY*atomJ.inducedDipoleP[0] + atomI.labFrameQuadrupole_YY*atomJ.inducedDipoleP[1] + atomI.labFrameQuadrupole_YZ*atomJ.inducedDipoleP[2];
    qiukp3            = atomI.labFrameQuadrupole_XZ*atomJ.inducedDipoleP[0] + atomI.labFrameQuadrupole_YZ*atomJ.inducedDipoleP[1] + atomI.labFrameQuadrupole_ZZ*atomJ.inducedDipoleP[2];

    qkuip1            = atomJ.labFrameQuadrupole_XX*atomI.inducedDipoleP[0] + atomJ.labFrameQuadrupole_XY*atomI.inducedDipoleP[1] + atomJ.labFrameQuadrupole_XZ*atomI.inducedDipoleP[2];
    qkuip2            = atomJ.labFrameQuadrupole_XY*atomI.inducedDipoleP[0] + atomJ.labFrameQuadrupole_YY*atomI.inducedDipoleP[1] + atomJ.labFrameQuadrupole_YZ*atomI.inducedDipoleP[2];
    qkuip3            = atomJ.labFrameQuadrupole_XZ*atomI.inducedDipoleP[0] + atomJ.labFrameQuadrupole_YZ*atomI.inducedDipoleP[1] + atomJ.labFrameQuadrupole_ZZ*atomI.inducedDipoleP[2];

    uixqkr1           = atomI.inducedDipole[1]*qkr3 - atomI.inducedDipole[2]*qkr2;
    uixqkr2           = atomI.inducedDipole[2]*qkr1 - atomI.inducedDipole[0]*qkr3;
    uixqkr3           = atomI.inducedDipole[0]*qkr2 - atomI.inducedDipole[1]*qkr1;

    ukxqir1           = atomJ.inducedDipole[1]*qir3 - atomJ.inducedDipole[2]*qir2;
    ukxqir2           = atomJ.inducedDipole[2]*qir1 - atomJ.inducedDipole[0]*qir3;
    ukxqir3           = atomJ.inducedDipole[0]*qir2 - atomJ.inducedDipole[1]*qir1;

    uixqkrp1          = atomI.inducedDipoleP[1]*qkr3 - atomI.inducedDipoleP[2]*qkr2;
    uixqkrp2          = atomI.inducedDipoleP[2]*qkr1 - atomI.inducedDipoleP[0]*qkr3;
    uixqkrp3          = atomI.inducedDipoleP[0]*qkr2 - atomI.inducedDipoleP[1]*qkr1;

    ukxqirp1          = atomJ.inducedDipoleP[1]*qir3 - atomJ.inducedDipoleP[2]*qir2;
    ukxqirp2          = atomJ.inducedDipoleP[2]*qir1 - atomJ.inducedDipoleP[0]*qir3;
    ukxqirp3          = atomJ.inducedDipoleP[0]*qir2 - atomJ.inducedDipoleP[1]*qir1;

    rxqiuk1           = yr*qiuk3 - zr*qiuk2;
    rxqiuk2           = zr*qiuk1 - xr*qiuk3;
    rxqiuk3           = xr*qiuk2 - yr*qiuk1;

    rxqkui1           = yr*qkui3 - zr*qkui2;
    rxqkui2           = zr*qkui1 - xr*qkui3;
    rxqkui3           = xr*qkui2 - yr*qkui1;

    rxqiukp1          = yr*qiukp3 - zr*qiukp2;
    rxqiukp2          = zr*qiukp1 - xr*qiukp3;
    rxqiukp3          = xr*qiukp2 - yr*qiukp1;

    rxqkuip1          = yr*qkuip3 - zr*qkuip2;
    rxqkuip2          = zr*qkuip1 - xr*qkuip3;
    rxqkuip3          = xr*qkuip2 - yr*qkuip1;

    // get intermediate variables for induction energy terms

    sci1              = atomI.inducedDipole[0]*atomJ.labFrameDipole[0] + atomI.inducedDipole[1]*atomJ.labFrameDipole[1]
                          + atomI.inducedDipole[2]*atomJ.labFrameDipole[2] + atomI.labFrameDipole[0]*atomJ.inducedDipole[0]
                          + atomI.labFrameDipole[1]*atomJ.inducedDipole[1] + atomI.labFrameDipole[2]*atomJ.inducedDipole[2];

    sci3              = atomI.inducedDipole[0]*xr + atomI.inducedDipole[1]*yr + atomI.inducedDipole[2]*zr;
    sci4              = atomJ.inducedDipole[0]*xr + atomJ.inducedDipole[1]*yr + atomJ.inducedDipole[2]*zr;

    sci7              = qir1*atomJ.inducedDipole[0] + qir2*atomJ.inducedDipole[1] + qir3*atomJ.inducedDipole[2];
    sci8              = qkr1*atomI.inducedDipole[0] + qkr2*atomI.inducedDipole[1] + qkr3*atomI.inducedDipole[2];

    scip1             = atomI.inducedDipoleP[0]*atomJ.labFrameDipole[0] + atomI.inducedDipoleP[1]*atomJ.labFrameDipole[1] + atomI.inducedDipoleP[2]*atomJ.labFrameDipole[2] + atomI.labFrameDipole[0]*atomJ.inducedDipoleP[0] + atomI.labFrameDipole[1]*atomJ.inducedDipoleP[1] + atomI.labFrameDipole[2]*atomJ.inducedDipoleP[2];
    scip2             = atomI.inducedDipole[0]*atomJ.inducedDipoleP[0]+atomI.inducedDipole[1]*atomJ.inducedDipoleP[1] + atomI.inducedDipole[2]*atomJ.inducedDipoleP[2]+atomI.inducedDipoleP[0]*atomJ.inducedDipole[0] + atomI.inducedDipoleP[1]*atomJ.inducedDipole[1]+atomI.inducedDipoleP[2]*atomJ.inducedDipole[2];

    scip3             = atomI.inducedDipoleP[0]*xr + atomI.inducedDipoleP[1]*yr + atomI.inducedDipoleP[2]*zr;
    scip4             = atomJ.inducedDipoleP[0]*xr + atomJ.inducedDipoleP[1]*yr + atomJ.inducedDipoleP[2]*zr;

    scip7             = qir1*atomJ.inducedDipoleP[0] + qir2*atomJ.inducedDipoleP[1] + qir3*atomJ.inducedDipoleP[2];
    scip8             = qkr1*atomI.inducedDipoleP[0] + qkr2*atomI.inducedDipoleP[1] + qkr3*atomI.inducedDipoleP[2];

    // calculate the gl functions for potential energy

    gli1              = atomJ.q*sci3 - atomI.q*sci4;
    gli2              = -sc3*sci4 - sci3*sc4;
    gli3              = sci3*sc6 - sci4*sc5;
    gli6              = sci1;
    gli7              = 2.0f * (sci7-sci8);

    glip1             = atomJ.q*scip3 - atomI.q*scip4;
    glip2             = -sc3*scip4 - scip3*sc4;
    glip3             = scip3*sc6 - scip4*sc5;
    glip6             = scip1;
    glip7             = 2.0f * (scip7-scip8);

    // get the permanent multipole and induced energies

    *outputEnergy      += -0.5f * (rr3*(gli1+gli6)*psc3 + rr5*(gli2+gli7)*psc5 + rr7*gli3*psc7);

    // intermediate variables for the induced-permanent terms;

    gfi1             = 0.5f*rr5*((gli1+gli6)*psc3 + (glip1+glip6)*dsc3+scip2*scale3i)
                         + 0.5f*rr7*((gli7+gli2)*psc5 +(glip7+glip2)*dsc5
                         -(sci3*scip4+scip3*sci4)*scale5i)
                         + 0.5f*rr9*(gli3*psc7+glip3*dsc7);

    gfi4             = 2.0f * rr5;
    gfi5             = rr7 * (sci4*psc7+scip4*dsc7);
    gfi6             = -rr7 * (sci3*psc7+scip3*dsc7);

    // get the induced force

    ftm2i1           = gfi1*xr + 0.5f*
                        (- rr3*atomJ.q*(atomI.inducedDipole[0]*psc3+atomI.inducedDipoleP[0]*dsc3)
                         + rr5*sc4*(atomI.inducedDipole[0]*psc5+atomI.inducedDipoleP[0]*dsc5)
                         - rr7*sc6*(atomI.inducedDipole[0]*psc7+atomI.inducedDipoleP[0]*dsc7))
                         +(rr3*atomI.q*(atomJ.inducedDipole[0]*psc3+atomJ.inducedDipoleP[0]*dsc3)
                         + rr5*sc3*(atomJ.inducedDipole[0]*psc5+atomJ.inducedDipoleP[0]*dsc5)
                         + rr7*sc5*(atomJ.inducedDipole[0]*psc7+atomJ.inducedDipoleP[0]*dsc7))*0.5f
                         + rr5*scale5i*(sci4*atomI.inducedDipoleP[0]+scip4*atomI.inducedDipole[0]
                         + sci3*atomJ.inducedDipoleP[0]+scip3*atomJ.inducedDipole[0])*0.5f
                         + 0.5f*(sci4*psc5+scip4*dsc5)*rr5*atomI.labFrameDipole[0]
                         + 0.5f*(sci3*psc5+scip3*dsc5)*rr5*atomJ.labFrameDipole[0]
                         + 0.5f*gfi4*((qkui1-qiuk1)*psc5
                         + (qkuip1-qiukp1)*dsc5)
                         + gfi5*qir1 + gfi6*qkr1;

    ftm2i2           = gfi1*yr + 0.5f*
                        (- rr3*atomJ.q*(atomI.inducedDipole[1]*psc3+atomI.inducedDipoleP[1]*dsc3)
                         + rr5*sc4*(atomI.inducedDipole[1]*psc5+atomI.inducedDipoleP[1]*dsc5)
                         - rr7*sc6*(atomI.inducedDipole[1]*psc7+atomI.inducedDipoleP[1]*dsc7))
                         +(rr3*atomI.q*(atomJ.inducedDipole[1]*psc3+atomJ.inducedDipoleP[1]*dsc3)
                         + rr5*sc3*(atomJ.inducedDipole[1]*psc5+atomJ.inducedDipoleP[1]*dsc5)
                         + rr7*sc5*(atomJ.inducedDipole[1]*psc7+atomJ.inducedDipoleP[1]*dsc7))*0.5f
                         + rr5*scale5i*(sci4*atomI.inducedDipoleP[1]+scip4*atomI.inducedDipole[1]
                         + sci3*atomJ.inducedDipoleP[1]+scip3*atomJ.inducedDipole[1])*0.5f
                         + 0.5f*(sci4*psc5+scip4*dsc5)*rr5*atomI.labFrameDipole[1]
                         + 0.5f*(sci3*psc5+scip3*dsc5)*rr5*atomJ.labFrameDipole[1]
                         + 0.5f*gfi4*((qkui2-qiuk2)*psc5
                         + (qkuip2-qiukp2)*dsc5)
                         + gfi5*qir2 + gfi6*qkr2;

    ftm2i3           = gfi1*zr  + 0.5f*
                         (- rr3*atomJ.q*(atomI.inducedDipole[2]*psc3+atomI.inducedDipoleP[2]*dsc3)
                         + rr5*sc4*(atomI.inducedDipole[2]*psc5+atomI.inducedDipoleP[2]*dsc5)
                         - rr7*sc6*(atomI.inducedDipole[2]*psc7+atomI.inducedDipoleP[2]*dsc7))
                         +(rr3*atomI.q*(atomJ.inducedDipole[2]*psc3+atomJ.inducedDipoleP[2]*dsc3)
                         + rr5*sc3*(atomJ.inducedDipole[2]*psc5+atomJ.inducedDipoleP[2]*dsc5)
                         + rr7*sc5*(atomJ.inducedDipole[2]*psc7+atomJ.inducedDipoleP[2]*dsc7))*0.5f
                         + rr5*scale5i*(sci4*atomI.inducedDipoleP[2]+scip4*atomI.inducedDipole[2]
                         + sci3*atomJ.inducedDipoleP[2]+scip3*atomJ.inducedDipole[2])*0.5f
                         + 0.5f*(sci4*psc5+scip4*dsc5)*rr5*atomI.labFrameDipole[2]
                         + 0.5f*(sci3*psc5+scip3*dsc5)*rr5*atomJ.labFrameDipole[2]
                         + 0.5f*gfi4*((qkui3-qiuk3)*psc5
                         + (qkuip3-qiukp3)*dsc5)
                         + gfi5*qir3 + gfi6*qkr3;

    // intermediate values needed for partially excluded interactions

    fridmp1          = 0.5f * (rr3*((gli1+gli6)*pscale
                         +(glip1+glip6)*dscale)*ddsc3_1
                         + rr5*((gli2+gli7)*pscale
                         +(glip2+glip7)*dscale)*ddsc5_1
                         + rr7*(gli3*pscale+glip3*dscale)*ddsc7_1);

    fridmp2          = 0.5f * (rr3*((gli1+gli6)*pscale
                         +(glip1+glip6)*dscale)*ddsc3_2
                         + rr5*((gli2+gli7)*pscale
                         +(glip2+glip7)*dscale)*ddsc5_2
                         + rr7*(gli3*pscale+glip3*dscale)*ddsc7_2);

    fridmp3          = 0.5f * (rr3*((gli1+gli6)*pscale
                         +(glip1+glip6)*dscale)*ddsc3_3
                         + rr5*((gli2+gli7)*pscale
                         +(glip2+glip7)*dscale)*ddsc5_3
                         + rr7*(gli3*pscale+glip3*dscale)*ddsc7_3);

    // get the induced-induced derivative terms;

    findmp1          = 0.5f * uscale * (scip2*rr3*ddsc3_1
                         - rr5*ddsc5_1*(sci3*scip4+scip3*sci4));

    findmp2          = 0.5f * uscale * (scip2*rr3*ddsc3_2
                         - rr5*ddsc5_2*(sci3*scip4+scip3*sci4));

    findmp3          = 0.5f * uscale * (scip2*rr3*ddsc3_3
                         - rr5*ddsc5_3*(sci3*scip4+scip3*sci4));

    // handle of scaling for partially excluded interactions;

    ftm2i1           = ftm2i1 - fridmp1 - findmp1;
    ftm2i2           = ftm2i2 - fridmp2 - findmp2;
    ftm2i3           = ftm2i3 - fridmp3 - findmp3;

    // correction to convert mutual to direct polarization force;

#if 0
               if (poltyp .eq. 'DIRECT') then;
                  gfd = 0.5f * (rr5*scip2*scale3i;
     &                  - rr7*(scip3*sci4+sci3*scip4)*scale5i);
                  fdir1 = gfd*xr + 0.5f*rr5*scale5i;
     &                         * (sci4*atomI.inducedDipoleP[0]+scip4*atomI.inducedDipole[0];
     &                           +sci3*atomJ.inducedDipoleP[0]+scip3*atomJ.inducedDipole[0]);
                  fdir2 = gfd*yr + 0.5f*rr5*scale5i;
     &                         * (sci4*atomI.inducedDipoleP[1]+scip4*atomI.inducedDipole[1];
     &                           +sci3*atomJ.inducedDipoleP[1]+scip3*atomJ.inducedDipole[1]);
                  fdir3 = gfd*zr + 0.5f*rr5*scale5i;
     &                         * (sci4*atomI.inducedDipoleP[2]+scip4*atomI.inducedDipole[2];
     &                           +sci3*atomJ.inducedDipoleP[2]+scip3*atomJ.inducedDipole[2]);
                  ftm2i1 = ftm2i1 - fdir1 + findmp1;
                  ftm2i2 = ftm2i2 - fdir2 + findmp2;
                  ftm2i3 = ftm2i3 - fdir3 + findmp3;
               end if;
#endif

    // now perform the torque calculation
    // intermediate terms for torque between multipoles i and k

    gti2             = 0.5f * (sci4*psc5+scip4*dsc5) * rr5;
    gti3             = 0.5f * (sci3*psc5+scip3*dsc5) * rr5;
    gti4             = gfi4;
    gti5             = gfi5;
    gti6             = gfi6;

    // calculate the induced torque components

    ttm2i1           = -rr3*(dixuk1*psc3+dixukp1*dsc3)*0.5f
                         + gti2*dixr1 + gti4*((ukxqir1+rxqiuk1)*psc5
                         +(ukxqirp1+rxqiukp1)*dsc5)*0.5f - gti5*rxqir1;

    ttm2i2           = -rr3*(dixuk2*psc3+dixukp2*dsc3)*0.5f
                         + gti2*dixr2 + gti4*((ukxqir2+rxqiuk2)*psc5
                         +(ukxqirp2+rxqiukp2)*dsc5)*0.5f - gti5*rxqir2;

    ttm2i3           = -rr3*(dixuk3*psc3+dixukp3*dsc3)*0.5f
                         + gti2*dixr3 + gti4*((ukxqir3+rxqiuk3)*psc5
                         +(ukxqirp3+rxqiukp3)*dsc5)*0.5f - gti5*rxqir3;

    ttm3i1           = -rr3*(dkxui1*psc3+dkxuip1*dsc3)*0.5f
                         + gti3*dkxr1 - gti4*((uixqkr1+rxqkui1)*psc5
                         +(uixqkrp1+rxqkuip1)*dsc5)*0.5f - gti6*rxqkr1;

    ttm3i2           = -rr3*(dkxui2*psc3+dkxuip2*dsc3)*0.5f
                         + gti3*dkxr2 - gti4*((uixqkr2+rxqkui2)*psc5
                         +(uixqkrp2+rxqkuip2)*dsc5)*0.5f - gti6*rxqkr2;

    ttm3i3           = -rr3*(dkxui3*psc3+dkxuip3*dsc3)*0.5f
                         + gti3*dkxr3 - gti4*((uixqkr3+rxqkui3)*psc5
                         +(uixqkrp3+rxqkuip3)*dsc5)*0.5f - gti6*rxqkr3;

    // update force and torque on site k;

#if 0
    ftm1i(1,k) = ftm1i(1,k) - ftm2i1;
    ftm1i(2,k) = ftm1i(2,k) - ftm2i2;
    ftm1i(3,k) = ftm1i(3,k) - ftm2i3;
    ttm1i(1,k) = ttm1i(1,k) - ttm3i1;
    ttm1i(2,k) = ttm1i(2,k) - ttm3i2;
    ttm1i(3,k) = ttm1i(3,k) - ttm3i3;

    // update force and torque on site i

    ftm1i(1,i) = ftm1i(1,i) + ftm2i1;
    ftm1i(2,i) = ftm1i(2,i) + ftm2i2;
    ftm1i(3,i) = ftm1i(3,i) + ftm2i3;
    ttm1i(1,i) = ttm1i(1,i) - ttm2i1;
    ttm1i(2,i) = ttm1i(2,i) - ttm2i2;
    ttm1i(3,i) = ttm1i(3,i) - ttm2i3;
#endif
 
    outputForce[0]     += ftm2i1;
    outputForce[1]     += ftm2i2;
    outputForce[2]     += ftm2i3;

    outputTorqueI[0]   -= ttm2i1;
    outputTorqueI[1]   -= ttm2i2;
    outputTorqueI[2]   -= ttm2i3;

    outputTorqueJ[0]   -= ttm3i1;
    outputTorqueJ[1]   -= ttm3i2;
    outputTorqueJ[2]   -= ttm3i3;

}

#ifdef AMOEBA_DEBUG
__device__ static int debugAccumulate( unsigned int index, float4* debugArray, float* field, unsigned int addMask, float idLabel )
{
    index                             += cAmoebaSim.paddedNumberOfAtoms;
    debugArray[index].x                = addMask ? field[0] : 0.0f;
    debugArray[index].y                = addMask ? field[1] : 0.0f;
    debugArray[index].z                = addMask ? field[2] : 0.0f;
    debugArray[index].w                = idLabel;

    return index;
}
#endif

__device__ void zeroKirkwoodEDiffParticleSharedField( struct KirkwoodEDiffParticle* sA )
{
    // zero shared fields

    sA->force[0]              = 0.0f;
    sA->force[1]              = 0.0f;
    sA->force[2]              = 0.0f;

    sA->torque[0]             = 0.0f;
    sA->torque[1]             = 0.0f;
    sA->torque[2]             = 0.0f;

}

// Include versions of the kernels for N^2 calculations.

#undef USE_OUTPUT_BUFFER_PER_WARP
#define METHOD_NAME(a, b) a##N2##b
#include "kCalculateAmoebaCudaKirkwoodEDiff.h"
#define USE_OUTPUT_BUFFER_PER_WARP
#undef METHOD_NAME
#define METHOD_NAME(a, b) a##N2ByWarp##b
#include "kCalculateAmoebaCudaKirkwoodEDiff.h"

// reduce psWorkArray_3_1 -> force
// reduce psWorkArray_3_2 -> torque

static void kReduceForceTorque( amoebaGpuContext amoebaGpu )
{

    kReduceFields_kernel<<<amoebaGpu->nonbondBlocks, amoebaGpu->fieldReduceThreadsPerBlock>>>(
                           amoebaGpu->paddedNumberOfAtoms*3, amoebaGpu->outputBuffers,
                           amoebaGpu->psWorkArray_3_1->_pDevStream[0], amoebaGpu->psKirkwoodEDiffForce->_pDevStream[0] );
    LAUNCHERROR("kReduceForceTorqueKirkwoodEDiff1");

    kReduceFields_kernel<<<amoebaGpu->nonbondBlocks, amoebaGpu->fieldReduceThreadsPerBlock>>>(
                           amoebaGpu->paddedNumberOfAtoms*3, amoebaGpu->outputBuffers,
                           amoebaGpu->psWorkArray_3_2->_pDevStream[0], amoebaGpu->psTorque->_pDevStream[0] );

    LAUNCHERROR("kReduceForceTorqueKirkwoodEDiff2");
}

#ifdef AMOEBA_DEBUG
//#if 1
static void printKirkwoodEDiffBuffer( amoebaGpuContext amoebaGpu, unsigned int bufferIndex )
{
    (void) fprintf( amoebaGpu->log, "KirkwoodEDiff Buffer %u\n", bufferIndex );
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

static void printKirkwoodEDiffAtomBuffers( amoebaGpuContext amoebaGpu, unsigned int targetAtom )
{
    (void) fprintf( amoebaGpu->log, "KirkwoodEDiff atom %u\n", targetAtom );
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

void kCalculateAmoebaKirkwoodEDiff( amoebaGpuContext amoebaGpu )
{
  
   // ---------------------------------------------------------------------------------------

    static unsigned int threadsPerBlock = 0;
    static int timestep                 = 0;
    timestep++;
#ifdef AMOEBA_DEBUG
    static const char* methodName       = "kCalculateAmoebaKirkwoodEDiff";
    std::vector<int> fileId;
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
        (void) fflush( amoebaGpu->log );
    }   
    int paddedNumberOfAtoms                   = amoebaGpu->gpuContext->sim.paddedNumberOfAtoms;
    CUDAStream<float4>* debugArray            = new CUDAStream<float4>(paddedNumberOfAtoms*paddedNumberOfAtoms, 1, "DebugArray");
    memset( debugArray->_pSysStream[0],      0, sizeof( float )*4*paddedNumberOfAtoms*paddedNumberOfAtoms);
    debugArray->Upload();
    unsigned int targetAtom                   = 0;
#endif

    kClearFields_3( amoebaGpu, 6 );

    if( threadsPerBlock == 0 ){
        unsigned int maxThreads;
        if (gpu->sm_version >= SM_20)
            maxThreads = 384;
        else if (gpu->sm_version >= SM_12)
            maxThreads = 96;
        else
            maxThreads = 32;
        threadsPerBlock = std::min(getThreadsPerBlock( amoebaGpu, sizeof(KirkwoodEDiffParticle)), maxThreads);
    }   
    
    if( amoebaGpu->log && timestep == 1 ){
        (void) fprintf( amoebaGpu->log, "kCalculateAmoebaCudaKirkwoodEDiffN2Forces: blocks=%u threads=%u bffr/Warp=%u atm=%lu shrd=%lu"
                                        " Ebuf=%u ixnCt=%lu workUnits=%u sm=%d device=%d sharedMemoryPerBlock=%u\n",
                        amoebaGpu->nonbondBlocks, threadsPerBlock, amoebaGpu->bOutputBufferPerWarp,
                        sizeof(KirkwoodEDiffParticle), sizeof(KirkwoodEDiffParticle)*threadsPerBlock,
                        amoebaGpu->energyOutputBuffers, (*gpu->psInteractionCount)[0], gpu->sim.workUnits, gpu->sm_version, gpu->device, gpu->sharedMemoryPerBlock );
        //gpuPrintCudaAmoebaGmxSimulation(amoebaGpu, amoebaGpu->log );
        (void) fflush( amoebaGpu->log );
    }   

    if (gpu->bOutputBufferPerWarp){
#if 0
        (void) fprintf( amoebaGpu->log, "kCalculateAmoebaCudaKirkwoodEDiffN2Forces warp:  numBlocks=%u numThreads=%u bufferPerWarp=%u atm=%u shrd=%u Ebuf=%u ixnCt=%u workUnits=%u\n",
                        amoebaGpu->nonbondBlocks, amoebaGpu->nonbondElectrostaticThreadsPerBlock, amoebaGpu->bOutputBufferPerWarp,
                        sizeof(KirkwoodEDiffParticle), sizeof(KirkwoodEDiffParticle)*amoebaGpu->nonbondElectrostaticThreadsPerBlock, amoebaGpu->energyOutputBuffers, (*gpu->psInteractionCount)[0], gpu->sim.workUnits );
        (void) fflush( amoebaGpu->log );
        kCalculateAmoebaCudaKirkwoodEDiffN2ByWarpForces_kernel<<<amoebaGpu->nonbondBlocks, amoebaGpu->nonbondElectrostaticThreadsPerBlock,
                                                            sizeof(KirkwoodEDiffParticle)*amoebaGpu->nonbondElectrostaticThreadsPerBlock>>>(

                                                         amoebaGpu->psWorkUnit->_pDevStream[0],
                                                         gpu->psPosq4->_pDevStream[0],
                                                         amoebaGpu->psLabFrameDipole->_pDevStream[0],
                                                         amoebaGpu->psLabFrameQuadrupole->_pDevStream[0],
                                                         amoebaGpu->psInducedDipole->_pDevStream[0],
                                                         amoebaGpu->psInducedDipolePolar->_pDevStream[0],
                                                         amoebaGpu->psWorkArray_3_1->_pDevStream[0],
                                                         amoebaGpu->psWorkArray_3_2->_pDevStream[0],
#ifdef AMOEBA_DEBUG
                                                                           amoebaGpu->psWorkArray_3_2->_pDevStream[0],
                                                                           debugArray->_pDevStream[0], targetAtom );
#else
                                                                           amoebaGpu->psWorkArray_3_2->_pDevStream[0] );
#endif
#endif

    } else {

        kCalculateAmoebaCudaKirkwoodEDiffN2Forces_kernel<<<amoebaGpu->nonbondBlocks, threadsPerBlock, sizeof(KirkwoodEDiffParticle)*threadsPerBlock>>>(
                                                                           amoebaGpu->psWorkUnit->_pDevStream[0],
                                                                           gpu->psPosq4->_pDevStream[0],
                                                                           amoebaGpu->psLabFrameDipole->_pDevStream[0],
                                                                           amoebaGpu->psLabFrameQuadrupole->_pDevStream[0],
                                                                           amoebaGpu->psInducedDipole->_pDevStream[0],
                                                                           amoebaGpu->psInducedDipolePolar->_pDevStream[0],
                                                                           amoebaGpu->psInducedDipoleS->_pDevStream[0],
                                                                           amoebaGpu->psInducedDipolePolarS->_pDevStream[0],
                                                                           amoebaGpu->psWorkArray_3_1->_pDevStream[0],
#ifdef AMOEBA_DEBUG
                                                                           amoebaGpu->psWorkArray_3_2->_pDevStream[0],
                                                                           debugArray->_pDevStream[0], targetAtom );
#else
                                                                           amoebaGpu->psWorkArray_3_2->_pDevStream[0] );
#endif
    }

    LAUNCHERROR("kCalculateAmoebaCudaKirkwoodEDiffN2Forces");

    kReduceForceTorque( amoebaGpu );

#ifdef AMOEBA_DEBUG
    if( amoebaGpu->log ){

        amoebaGpu->psWorkArray_3_1->Download();
        amoebaGpu->psWorkArray_3_2->Download();

        //printKirkwoodEDiffAtomBuffers( amoebaGpu, (targetAtom + 0) );
        //printKirkwoodEDiffAtomBuffers( amoebaGpu, (targetAtom + 1231) );
        //printKirkwoodEDiffBuffer( amoebaGpu, 0 );
        //printKirkwoodEDiffBuffer( amoebaGpu, 38 );

        amoebaGpu->psKirkwoodEDiffForce->Download();
        amoebaGpu->psTorque->Download();
        debugArray->Download();

        int maxPrint        = 1400;
        for( int ii = 0; ii < gpu->natoms; ii++ ){
           (void) fprintf( amoebaGpu->log, "%5d ", ii); 

            int indexOffset     = ii*3;
    
           // force

           (void) fprintf( amoebaGpu->log,"KirkwoodEDiffF [%16.9e %16.9e %16.9e] ",
                           amoebaGpu->psKirkwoodEDiffForce->_pSysStream[0][indexOffset],
                           amoebaGpu->psKirkwoodEDiffForce->_pSysStream[0][indexOffset+1],
                           amoebaGpu->psKirkwoodEDiffForce->_pSysStream[0][indexOffset+2] );
    
           // torque

           (void) fprintf( amoebaGpu->log,"T [%16.9e %16.9e %16.9e] ",
                           amoebaGpu->psTorque->_pSysStream[0][indexOffset],
                           amoebaGpu->psTorque->_pSysStream[0][indexOffset+1],
                           amoebaGpu->psTorque->_pSysStream[0][indexOffset+2] );

           // coords

#if 0
            (void) fprintf( amoebaGpu->log,"x[%16.9e %16.9e %16.9e] ",
                            gpu->psPosq4->_pSysStream[0][ii].x,
                            gpu->psPosq4->_pSysStream[0][ii].y,
                            gpu->psPosq4->_pSysStream[0][ii].z);


           for( int jj = 0; jj < gpu->natoms && jj < 5; jj++ ){
               int debugIndex = jj*gpu->natoms + ii;
               float xx       =  gpu->psPosq4->_pSysStream[0][jj].x -  gpu->psPosq4->_pSysStream[0][ii].x;
               float yy       =  gpu->psPosq4->_pSysStream[0][jj].y -  gpu->psPosq4->_pSysStream[0][ii].y;
               float zz       =  gpu->psPosq4->_pSysStream[0][jj].z -  gpu->psPosq4->_pSysStream[0][ii].z;
               (void) fprintf( amoebaGpu->log,"\n%4d %4d delta [%16.9e %16.9e %16.9e] [%16.9e %16.9e %16.9e] ",
                               ii, jj, xx, yy, zz,
                               debugArray->_pSysStream[0][debugIndex].x, debugArray->_pSysStream[0][debugIndex].y, debugArray->_pSysStream[0][debugIndex].z );

           }
#endif
           if( ii == targetAtom ){
               (void) fprintf( amoebaGpu->log,"\n" );
               int paddedNumberOfAtoms                    = amoebaGpu->gpuContext->sim.paddedNumberOfAtoms;
               for( int jj = 0; jj < gpu->natoms; jj++ ){
                   int debugIndex = jj;
                   (void) fprintf( amoebaGpu->log,"%5d %5d ediff F%T\n", ii, jj );
                   for( int kk = 0; kk < 5; kk++ ){
                       (void) fprintf( amoebaGpu->log,"[%16.9e %16.9e %16.9e %16.9e]\n",
                                       debugArray->_pSysStream[0][debugIndex].x, debugArray->_pSysStream[0][debugIndex].y,
                                       debugArray->_pSysStream[0][debugIndex].z, debugArray->_pSysStream[0][debugIndex].w );
                       debugIndex += paddedNumberOfAtoms;
                   }
                   (void) fprintf( amoebaGpu->log,"\n" );

               }

               (void) fprintf( amoebaGpu->log,"\n" );
           }
           (void) fprintf( amoebaGpu->log,"\n" );
           if( ii == maxPrint && (gpu->natoms - maxPrint) > ii ){
                ii = gpu->natoms - maxPrint;
           }
        }
        (void) fflush( amoebaGpu->log );
        {
            (void) fprintf( amoebaGpu->log, "%s Tiled F & T\n", methodName ); fflush( amoebaGpu->log );
            int maxPrint = 12;
            for( int ii = 0; ii < gpu->natoms; ii++ ){
    
                // print cpu & gpu reductions
    
                int offset  = 3*ii;
    
                (void) fprintf( amoebaGpu->log,"%6d F[%16.7e %16.7e %16.7e] T[%16.7e %16.7e %16.7e]\n", ii,
                                amoebaGpu->psKirkwoodEDiffForce->_pSysStream[0][offset],
                                amoebaGpu->psKirkwoodEDiffForce->_pSysStream[0][offset+1],
                                amoebaGpu->psKirkwoodEDiffForce->_pSysStream[0][offset+2],
                                amoebaGpu->psTorque->_pSysStream[0][offset],
                                amoebaGpu->psTorque->_pSysStream[0][offset+1],
                                amoebaGpu->psTorque->_pSysStream[0][offset+2] );
                if( (ii == maxPrint) && (ii < (gpu->natoms - maxPrint)) )ii = gpu->natoms - maxPrint; 
            }   
        }   

        if( 1 ){
            std::vector<int> fileId;
            //fileId.push_back( 0 );
            VectorOfDoubleVectors outputVector;
            cudaLoadCudaFloat4Array( gpu->natoms, 3, gpu->psPosq4,                    outputVector );
            cudaLoadCudaFloatArray( gpu->natoms,  3, amoebaGpu->psKirkwoodEDiffForce,      outputVector );
            cudaLoadCudaFloatArray( gpu->natoms,  3, amoebaGpu->psTorque,             outputVector);
            cudaWriteVectorOfDoubleVectorsToFile( "CudaForceTorque", fileId, outputVector );
         }

    }   
    delete debugArray;

#endif

    // map torques to forces

    cudaComputeAmoebaMapTorquesAndAddTotalForce( amoebaGpu, amoebaGpu->psTorque, amoebaGpu->psKirkwoodEDiffForce, amoebaGpu->gpuContext->psForce4 );

#ifdef AMOEBA_DEBUG
    if( amoebaGpu->log ){

        cudaComputeAmoebaMapTorques( amoebaGpu, amoebaGpu->psTorque, amoebaGpu->psKirkwoodEDiffForce );
        amoebaGpu->psKirkwoodEDiffForce->Download();
        (void) fprintf( amoebaGpu->log, "Mapped KirkwoodEDiff torques to forces.\n" ); (void) fflush( amoebaGpu->log );

        int maxPrint        = 1400;
        for( int ii = 0; ii < gpu->natoms; ii++ ){
           (void) fprintf( amoebaGpu->log, "%5d ", ii); 

            int indexOffset     = ii*3;
    
           // force

           (void) fprintf( amoebaGpu->log,"KirkwoodEDiffF [%16.9e %16.9e %16.9e] ",
                           amoebaGpu->psKirkwoodEDiffForce->_pSysStream[0][indexOffset],
                           amoebaGpu->psKirkwoodEDiffForce->_pSysStream[0][indexOffset+1],
                           amoebaGpu->psKirkwoodEDiffForce->_pSysStream[0][indexOffset+2] );
    
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
            cudaLoadCudaFloatArray( gpu->natoms,  3, amoebaGpu->psKirkwoodEDiffForce,      outputVector );
            cudaWriteVectorOfDoubleVectorsToFile( "CudaKirkwoodEDiffForce", fileId, outputVector );
         }

    }   
#endif

   // ---------------------------------------------------------------------------------------
}

