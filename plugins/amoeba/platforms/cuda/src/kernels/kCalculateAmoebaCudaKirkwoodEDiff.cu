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
    float fridmp[4],findmp[4];
    float ftm2i[4];
    float ttm2i[4],ttm3i[4];
    //float dixdk[4],fdir[4];
    float dixuk[4],dkxui[4];
    float dixukp[4],dkxuip[4];
    float uixqkr[4],ukxqir[4];
    float uixqkrp[4],ukxqirp[4];
    float qiuk[4],qkui[4];
    float qiukp[4],qkuip[4];
    float rxqiuk[4],rxqkui[4];
    float rxqiukp[4],rxqkuip[4];
    //float qidk[4],qkdi[4];
    float qir[4],qkr[4];
    //float qiqkr[4],qkqir[4];
    //float qixqk[4];
    float rxqir[4];
    float dixr[4],dkxr[4];
    //float dixqkr[4],dkxqir[4];
    float rxqkr[4];
    //float qkrxqir[4];
    //float rxqikr[4];
    //float rxqkir[4];
    //float rxqidk[4],rxqkdi[4];
    float ddsc3[4],ddsc5[4];
    float ddsc7[4];
    float gli[8],glip[8];
    float sc[11];
    float sci[9],scip[9];
    float gfi[7],gti[7];
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

    ddsc3[1]    = 0.0f;
    ddsc3[2]    = 0.0f;
    ddsc3[3]    = 0.0f;

    ddsc5[1]    = 0.0f;
    ddsc5[2]    = 0.0f;
    ddsc5[3]    = 0.0f;

    ddsc7[1]    = 0.0f;
    ddsc7[2]    = 0.0f;
    ddsc7[3]    = 0.0f;

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

            ddsc3[1]     = -3.0f*damp*exp(damp) * xr/r2;
            ddsc3[2]     = -3.0f*damp*exp(damp) * yr/r2;
            ddsc3[3]     = -3.0f*damp*exp(damp) * zr/r2;

            ddsc5[1]     = -damp * ddsc3[1];
            ddsc5[2]     = -damp * ddsc3[2];
            ddsc5[3]     = -damp * ddsc3[3];

            ddsc7[1]     = (-0.2f-0.6f*damp) * ddsc5[1];
            ddsc7[2]     = (-0.2f-0.6f*damp) * ddsc5[2];
            ddsc7[3]     = (-0.2f-0.6f*damp) * ddsc5[3];
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
    dixdk[1]            = atomI.labFrameDipole[1]*atomJ.labFrameDipole[2] - atomI.labFrameDipole[2]*atomJ.labFrameDipole[1];
    dixdk[2]            = atomI.labFrameDipole[2]*atomJ.labFrameDipole[0] - atomI.labFrameDipole[0]*atomJ.labFrameDipole[2];
    dixdk[3]            = atomI.labFrameDipole[0]*atomJ.labFrameDipole[1] - atomI.labFrameDipole[1]*atomJ.labFrameDipole[0];
#endif

    dixr[1]             = atomI.labFrameDipole[1]*zr - atomI.labFrameDipole[2]*yr;
    dixr[2]             = atomI.labFrameDipole[2]*xr - atomI.labFrameDipole[0]*zr;
    dixr[3]             = atomI.labFrameDipole[0]*yr - atomI.labFrameDipole[1]*xr;

    dkxr[1]             = atomJ.labFrameDipole[1]*zr - atomJ.labFrameDipole[2]*yr;
    dkxr[2]             = atomJ.labFrameDipole[2]*xr - atomJ.labFrameDipole[0]*zr;
    dkxr[3]             = atomJ.labFrameDipole[0]*yr - atomJ.labFrameDipole[1]*xr;

    qir[1]              = atomI.labFrameQuadrupole_XX*xr + atomI.labFrameQuadrupole_XY*yr + atomI.labFrameQuadrupole_XZ*zr;
    qir[2]              = atomI.labFrameQuadrupole_XY*xr + atomI.labFrameQuadrupole_YY*yr + atomI.labFrameQuadrupole_YZ*zr;
    qir[3]              = atomI.labFrameQuadrupole_XZ*xr + atomI.labFrameQuadrupole_YZ*yr + atomI.labFrameQuadrupole_ZZ*zr;

    qkr[1]              = atomJ.labFrameQuadrupole_XX*xr + atomJ.labFrameQuadrupole_XY*yr + atomJ.labFrameQuadrupole_XZ*zr;
    qkr[2]              = atomJ.labFrameQuadrupole_XY*xr + atomJ.labFrameQuadrupole_YY*yr + atomJ.labFrameQuadrupole_YZ*zr;
    qkr[3]              = atomJ.labFrameQuadrupole_XZ*xr + atomJ.labFrameQuadrupole_YZ*yr + atomJ.labFrameQuadrupole_ZZ*zr;

#if 0
    qiqkr[1]            = atomI.labFrameQuadrupole_XX*qkr[1] + atomI.labFrameQuadrupole_XY*qkr[2] + atomI.labFrameQuadrupole_XZ*qkr[3];
    qiqkr[2]            = atomI.labFrameQuadrupole_XY*qkr[1] + atomI.labFrameQuadrupole_YY*qkr[2] + atomI.labFrameQuadrupole_YZ*qkr[3];
    qiqkr[3]            = atomI.labFrameQuadrupole_XZ*qkr[1] + atomI.labFrameQuadrupole_YZ*qkr[2] + atomI.labFrameQuadrupole_ZZ*qkr[3];

    qkqir[1]            = atomJ.labFrameQuadrupole_XX*qir[1] + atomJ.labFrameQuadrupole_XY*qir[2] + atomJ.labFrameQuadrupole_XZ*qir[3];
    qkqir[2]            = atomJ.labFrameQuadrupole_XY*qir[1] + atomJ.labFrameQuadrupole_YY*qir[2] + atomJ.labFrameQuadrupole_YZ*qir[3];
    qkqir[3]            = atomJ.labFrameQuadrupole_XZ*qir[1] + atomJ.labFrameQuadrupole_YZ*qir[2] + atomJ.labFrameQuadrupole_ZZ*qir[3];

    qixqk[1]            = atomI.labFrameQuadrupole_XY*atomJ.labFrameQuadrupole_XZ + atomI.labFrameQuadrupole_YY*atomJ.labFrameQuadrupole_YZ + atomI.labFrameQuadrupole_YZ*atomJ.labFrameQuadrupole_ZZ -
                          atomI.labFrameQuadrupole_XZ*atomJ.labFrameQuadrupole_XY - atomI.labFrameQuadrupole_YZ*atomJ.labFrameQuadrupole_YY - atomI.labFrameQuadrupole_ZZ*atomJ.labFrameQuadrupole_YZ;

    qixqk[2]            = atomI.labFrameQuadrupole_XZ*atomJ.labFrameQuadrupole_XX + atomI.labFrameQuadrupole_YZ*atomJ.labFrameQuadrupole_XY + atomI.labFrameQuadrupole_ZZ*atomJ.labFrameQuadrupole_XZ -
                          atomI.labFrameQuadrupole_XX*atomJ.labFrameQuadrupole_XZ - atomI.labFrameQuadrupole_XY*atomJ.labFrameQuadrupole_YZ - atomI.labFrameQuadrupole_XZ*atomJ.labFrameQuadrupole_ZZ;

    qixqk[3]            = atomI.labFrameQuadrupole_XX*atomJ.labFrameQuadrupole_XY + atomI.labFrameQuadrupole_XY*atomJ.labFrameQuadrupole_YY + atomI.labFrameQuadrupole_XZ*atomJ.labFrameQuadrupole_YZ -
                          atomI.labFrameQuadrupole_XY*atomJ.labFrameQuadrupole_XX - atomI.labFrameQuadrupole_YY*atomJ.labFrameQuadrupole_XY - atomI.labFrameQuadrupole_YZ*atomJ.labFrameQuadrupole_XZ;
#endif

    rxqir[1]            = yr*qir[3] - zr*qir[2];
    rxqir[2]            = zr*qir[1] - xr*qir[3];
    rxqir[3]            = xr*qir[2] - yr*qir[1];

    rxqkr[1]            = yr*qkr[3] - zr*qkr[2];
    rxqkr[2]            = zr*qkr[1] - xr*qkr[3];
    rxqkr[3]            = xr*qkr[2] - yr*qkr[1];

#if 0
    rxqikr[1]           = yr*qiqkr[3] - zr*qiqkr[2];
    rxqikr[2]           = zr*qiqkr[1] - xr*qiqkr[3];
    rxqikr[3]           = xr*qiqkr[2] - yr*qiqkr[1];

    rxqkir[1]           = yr*qkqir[3] - zr*qkqir[2];
    rxqkir[2]           = zr*qkqir[1] - xr*qkqir[3];
    rxqkir[3]           = xr*qkqir[2] - yr*qkqir[1];

    qkrxqir[1]          = qkr[2]*qir[3] - qkr[3]*qir[2];
    qkrxqir[2]          = qkr[3]*qir[1] - qkr[1]*qir[3];
    qkrxqir[3]          = qkr[1]*qir[2] - qkr[2]*qir[1];

    qidk[1]             = atomI.labFrameQuadrupole_XX*atomJ.labFrameDipole[0] + atomI.labFrameQuadrupole_XY*atomJ.labFrameDipole[1] + atomI.labFrameQuadrupole_XZ*atomJ.labFrameDipole[2];
    qidk[2]             = atomI.labFrameQuadrupole_XY*atomJ.labFrameDipole[0] + atomI.labFrameQuadrupole_YY*atomJ.labFrameDipole[1] + atomI.labFrameQuadrupole_YZ*atomJ.labFrameDipole[2];
    qidk[3]             = atomI.labFrameQuadrupole_XZ*atomJ.labFrameDipole[0] + atomI.labFrameQuadrupole_YZ*atomJ.labFrameDipole[1] + atomI.labFrameQuadrupole_ZZ*atomJ.labFrameDipole[2];

    qkdi[1]             = atomJ.labFrameQuadrupole_XX*atomI.labFrameDipole[0] + atomJ.labFrameQuadrupole_XY*atomI.labFrameDipole[1] + atomJ.labFrameQuadrupole_XZ*atomI.labFrameDipole[2];
    qkdi[2]             = atomJ.labFrameQuadrupole_XY*atomI.labFrameDipole[0] + atomJ.labFrameQuadrupole_YY*atomI.labFrameDipole[1] + atomJ.labFrameQuadrupole_YZ*atomI.labFrameDipole[2];
    qkdi[3]             = atomJ.labFrameQuadrupole_XZ*atomI.labFrameDipole[0] + atomJ.labFrameQuadrupole_YZ*atomI.labFrameDipole[1] + atomJ.labFrameQuadrupole_ZZ*atomI.labFrameDipole[2];

    dixqkr[1]           = atomI.labFrameDipole[1]*qkr[3] - atomI.labFrameDipole[2]*qkr[2];
    dixqkr[2]           = atomI.labFrameDipole[2]*qkr[1] - atomI.labFrameDipole[0]*qkr[3];
    dixqkr[3]           = atomI.labFrameDipole[0]*qkr[2] - atomI.labFrameDipole[1]*qkr[1];

    dkxqir[1]           = atomJ.labFrameDipole[1]*qir[3] - atomJ.labFrameDipole[2]*qir[2];
    dkxqir[2]           = atomJ.labFrameDipole[2]*qir[1] - atomJ.labFrameDipole[0]*qir[3];
    dkxqir[3]           = atomJ.labFrameDipole[0]*qir[2] - atomJ.labFrameDipole[1]*qir[1];

    rxqidk[1]           = yr*qidk[3] - zr*qidk[2];
    rxqidk[2]           = zr*qidk[1] - xr*qidk[3];
    rxqidk[3]           = xr*qidk[2] - yr*qidk[1];

    rxqkdi[1]           = yr*qkdi[3] - zr*qkdi[2];
    rxqkdi[2]           = zr*qkdi[1] - xr*qkdi[3];
    rxqkdi[3]           = xr*qkdi[2] - yr*qkdi[1];
#endif
 
    // get intermediate variables for permanent energy terms
 
    sc[3]               = atomI.labFrameDipole[0]*xr  + atomI.labFrameDipole[1]*yr  + atomI.labFrameDipole[2]*zr;
    sc[4]               = atomJ.labFrameDipole[0]*xr  + atomJ.labFrameDipole[1]*yr  + atomJ.labFrameDipole[2]*zr;
    sc[5]               = qir[1]*xr + qir[2]*yr + qir[3]*zr;
    sc[6]               = qkr[1]*xr + qkr[2]*yr + qkr[3]*zr;
 
    // construct auxiliary vectors for induced terms
 
    dixuk[1]            = atomI.labFrameDipole[1]*atomJ.inducedDipoleS[2] - atomI.labFrameDipole[2]*atomJ.inducedDipoleS[1];
    dixuk[2]            = atomI.labFrameDipole[2]*atomJ.inducedDipoleS[0] - atomI.labFrameDipole[0]*atomJ.inducedDipoleS[2];
    dixuk[3]            = atomI.labFrameDipole[0]*atomJ.inducedDipoleS[1] - atomI.labFrameDipole[1]*atomJ.inducedDipoleS[0];

    dkxui[1]            = atomJ.labFrameDipole[1]*atomI.inducedDipoleS[2] - atomJ.labFrameDipole[2]*atomI.inducedDipoleS[1];
    dkxui[2]            = atomJ.labFrameDipole[2]*atomI.inducedDipoleS[0] - atomJ.labFrameDipole[0]*atomI.inducedDipoleS[2];
    dkxui[3]            = atomJ.labFrameDipole[0]*atomI.inducedDipoleS[1] - atomJ.labFrameDipole[1]*atomI.inducedDipoleS[0];

    dixukp[1]           = atomI.labFrameDipole[1]*atomJ.inducedDipolePS[2] - atomI.labFrameDipole[2]*atomJ.inducedDipolePS[1];
    dixukp[2]           = atomI.labFrameDipole[2]*atomJ.inducedDipolePS[0] - atomI.labFrameDipole[0]*atomJ.inducedDipolePS[2];
    dixukp[3]           = atomI.labFrameDipole[0]*atomJ.inducedDipolePS[1] - atomI.labFrameDipole[1]*atomJ.inducedDipolePS[0];

    dkxuip[1]           = atomJ.labFrameDipole[1]*atomI.inducedDipolePS[2] - atomJ.labFrameDipole[2]*atomI.inducedDipolePS[1];
    dkxuip[2]           = atomJ.labFrameDipole[2]*atomI.inducedDipolePS[0] - atomJ.labFrameDipole[0]*atomI.inducedDipolePS[2];
    dkxuip[3]           = atomJ.labFrameDipole[0]*atomI.inducedDipolePS[1] - atomJ.labFrameDipole[1]*atomI.inducedDipolePS[0];

    qiuk[1]             = atomI.labFrameQuadrupole_XX*atomJ.inducedDipoleS[0] + atomI.labFrameQuadrupole_XY*atomJ.inducedDipoleS[1] + atomI.labFrameQuadrupole_XZ*atomJ.inducedDipoleS[2];
    qiuk[2]             = atomI.labFrameQuadrupole_XY*atomJ.inducedDipoleS[0] + atomI.labFrameQuadrupole_YY*atomJ.inducedDipoleS[1] + atomI.labFrameQuadrupole_YZ*atomJ.inducedDipoleS[2];
    qiuk[3]             = atomI.labFrameQuadrupole_XZ*atomJ.inducedDipoleS[0] + atomI.labFrameQuadrupole_YZ*atomJ.inducedDipoleS[1] + atomI.labFrameQuadrupole_ZZ*atomJ.inducedDipoleS[2];

    qkui[1]             = atomJ.labFrameQuadrupole_XX*atomI.inducedDipoleS[0] + atomJ.labFrameQuadrupole_XY*atomI.inducedDipoleS[1] + atomJ.labFrameQuadrupole_XZ*atomI.inducedDipoleS[2];
    qkui[2]             = atomJ.labFrameQuadrupole_XY*atomI.inducedDipoleS[0] + atomJ.labFrameQuadrupole_YY*atomI.inducedDipoleS[1] + atomJ.labFrameQuadrupole_YZ*atomI.inducedDipoleS[2];
    qkui[3]             = atomJ.labFrameQuadrupole_XZ*atomI.inducedDipoleS[0] + atomJ.labFrameQuadrupole_YZ*atomI.inducedDipoleS[1] + atomJ.labFrameQuadrupole_ZZ*atomI.inducedDipoleS[2];

    qiukp[1]            = atomI.labFrameQuadrupole_XX*atomJ.inducedDipolePS[0] + atomI.labFrameQuadrupole_XY*atomJ.inducedDipolePS[1] + atomI.labFrameQuadrupole_XZ*atomJ.inducedDipolePS[2];
    qiukp[2]            = atomI.labFrameQuadrupole_XY*atomJ.inducedDipolePS[0] + atomI.labFrameQuadrupole_YY*atomJ.inducedDipolePS[1] + atomI.labFrameQuadrupole_YZ*atomJ.inducedDipolePS[2];
    qiukp[3]            = atomI.labFrameQuadrupole_XZ*atomJ.inducedDipolePS[0] + atomI.labFrameQuadrupole_YZ*atomJ.inducedDipolePS[1] + atomI.labFrameQuadrupole_ZZ*atomJ.inducedDipolePS[2];

    qkuip[1]            = atomJ.labFrameQuadrupole_XX*atomI.inducedDipolePS[0] + atomJ.labFrameQuadrupole_XY*atomI.inducedDipolePS[1] + atomJ.labFrameQuadrupole_XZ*atomI.inducedDipolePS[2];
    qkuip[2]            = atomJ.labFrameQuadrupole_XY*atomI.inducedDipolePS[0] + atomJ.labFrameQuadrupole_YY*atomI.inducedDipolePS[1] + atomJ.labFrameQuadrupole_YZ*atomI.inducedDipolePS[2];
    qkuip[3]            = atomJ.labFrameQuadrupole_XZ*atomI.inducedDipolePS[0] + atomJ.labFrameQuadrupole_YZ*atomI.inducedDipolePS[1] + atomJ.labFrameQuadrupole_ZZ*atomI.inducedDipolePS[2];

    uixqkr[1]           = atomI.inducedDipoleS[1]*qkr[3] - atomI.inducedDipoleS[2]*qkr[2];
    uixqkr[2]           = atomI.inducedDipoleS[2]*qkr[1] - atomI.inducedDipoleS[0]*qkr[3];
    uixqkr[3]           = atomI.inducedDipoleS[0]*qkr[2] - atomI.inducedDipoleS[1]*qkr[1];

    ukxqir[1]           = atomJ.inducedDipoleS[1]*qir[3] - atomJ.inducedDipoleS[2]*qir[2];
    ukxqir[2]           = atomJ.inducedDipoleS[2]*qir[1] - atomJ.inducedDipoleS[0]*qir[3];
    ukxqir[3]           = atomJ.inducedDipoleS[0]*qir[2] - atomJ.inducedDipoleS[1]*qir[1];

    uixqkrp[1]          = atomI.inducedDipolePS[1]*qkr[3] - atomI.inducedDipolePS[2]*qkr[2];
    uixqkrp[2]          = atomI.inducedDipolePS[2]*qkr[1] - atomI.inducedDipolePS[0]*qkr[3];
    uixqkrp[3]          = atomI.inducedDipolePS[0]*qkr[2] - atomI.inducedDipolePS[1]*qkr[1];

    ukxqirp[1]          = atomJ.inducedDipolePS[1]*qir[3] - atomJ.inducedDipolePS[2]*qir[2];
    ukxqirp[2]          = atomJ.inducedDipolePS[2]*qir[1] - atomJ.inducedDipolePS[0]*qir[3];
    ukxqirp[3]          = atomJ.inducedDipolePS[0]*qir[2] - atomJ.inducedDipolePS[1]*qir[1];

    rxqiuk[1]           = yr*qiuk[3] - zr*qiuk[2];
    rxqiuk[2]           = zr*qiuk[1] - xr*qiuk[3];
    rxqiuk[3]           = xr*qiuk[2] - yr*qiuk[1];

    rxqkui[1]           = yr*qkui[3] - zr*qkui[2];
    rxqkui[2]           = zr*qkui[1] - xr*qkui[3];
    rxqkui[3]           = xr*qkui[2] - yr*qkui[1];

    rxqiukp[1]          = yr*qiukp[3] - zr*qiukp[2];
    rxqiukp[2]          = zr*qiukp[1] - xr*qiukp[3];
    rxqiukp[3]          = xr*qiukp[2] - yr*qiukp[1];

    rxqkuip[1]          = yr*qkuip[3] - zr*qkuip[2];
    rxqkuip[2]          = zr*qkuip[1] - xr*qkuip[3];
    rxqkuip[3]          = xr*qkuip[2] - yr*qkuip[1];
 
    // get intermediate variables for induction energy terms

    sci[1]              = atomI.inducedDipoleS[0]*atomJ.labFrameDipole[0] + atomI.inducedDipoleS[1]*atomJ.labFrameDipole[1] +
                          atomI.inducedDipoleS[2]*atomJ.labFrameDipole[2] + atomI.labFrameDipole[0]*atomJ.inducedDipoleS[0] +
                          atomI.labFrameDipole[1]*atomJ.inducedDipoleS[1] + atomI.labFrameDipole[2]*atomJ.inducedDipoleS[2];

    sci[2]              = atomI.inducedDipoleS[0]*atomJ.inducedDipoleS[0] + atomI.inducedDipoleS[1]*atomJ.inducedDipoleS[1] + atomI.inducedDipoleS[2]*atomJ.inducedDipoleS[2];

    sci[3]              = atomI.inducedDipoleS[0]*xr + atomI.inducedDipoleS[1]*yr + atomI.inducedDipoleS[2]*zr;
    sci[4]              = atomJ.inducedDipoleS[0]*xr + atomJ.inducedDipoleS[1]*yr + atomJ.inducedDipoleS[2]*zr;

    sci[7]              = qir[1]*atomJ.inducedDipoleS[0] + qir[2]*atomJ.inducedDipoleS[1] + qir[3]*atomJ.inducedDipoleS[2];
    sci[8]              = qkr[1]*atomI.inducedDipoleS[0] + qkr[2]*atomI.inducedDipoleS[1] + qkr[3]*atomI.inducedDipoleS[2];

    scip[1]             = atomI.inducedDipolePS[0]*atomJ.labFrameDipole[0] + atomI.inducedDipolePS[1]*atomJ.labFrameDipole[1] +
                          atomI.inducedDipolePS[2]*atomJ.labFrameDipole[2] + atomI.labFrameDipole[0]*atomJ.inducedDipolePS[0] +
                          atomI.labFrameDipole[1]*atomJ.inducedDipolePS[1] + atomI.labFrameDipole[2]*atomJ.inducedDipolePS[2];

    scip[2]             = atomI.inducedDipoleS[0]*atomJ.inducedDipolePS[0] + atomI.inducedDipoleS[1]*atomJ.inducedDipolePS[1] +
                          atomI.inducedDipoleS[2]*atomJ.inducedDipolePS[2] + atomI.inducedDipolePS[0]*atomJ.inducedDipoleS[0] +
                          atomI.inducedDipolePS[1]*atomJ.inducedDipoleS[1] + atomI.inducedDipolePS[2]*atomJ.inducedDipoleS[2];

    scip[3]             = atomI.inducedDipolePS[0]*xr + atomI.inducedDipolePS[1]*yr + atomI.inducedDipolePS[2]*zr;
    scip[4]             = atomJ.inducedDipolePS[0]*xr + atomJ.inducedDipolePS[1]*yr + atomJ.inducedDipolePS[2]*zr;

    scip[7]             = qir[1]*atomJ.inducedDipolePS[0] + qir[2]*atomJ.inducedDipolePS[1] + qir[3]*atomJ.inducedDipolePS[2];
    scip[8]             = qkr[1]*atomI.inducedDipolePS[0] + qkr[2]*atomI.inducedDipolePS[1] + qkr[3]*atomI.inducedDipolePS[2];

 
    // calculate the gl functions for potential energy

    gli[1]              = atomJ.q*sci[3] - atomI.q*sci[4];
    gli[2]              = -sc[3]*sci[4] - sci[3]*sc[4];
    gli[3]              = sci[3]*sc[6] - sci[4]*sc[5];
    gli[6]              = sci[1];
    gli[7]              = 2.0f * (sci[7]-sci[8]);
    glip[1]             = atomJ.q*scip[3] - atomI.q*scip[4];
    glip[2]             = -sc[3]*scip[4] - scip[3]*sc[4];
    glip[3]             = scip[3]*sc[6] - scip[4]*sc[5];
    glip[6]             = scip[1];
    glip[7]             = 2.0f * (scip[7]-scip[8]);

    // get the permanent multipole and induced energies;

    *outputEnergy       = 0.5f * (rr3*(gli[1]+gli[6])*psc3 +
                                  rr5*(gli[2]+gli[7])*psc5 +
                                  rr7*gli[3]*psc7);

    // intermediate variables for the induced-permanent terms

    gfi[1]              = 0.5f*rr5*((gli[1]+gli[6])*psc3 +
                                    (glip[1]+glip[6])*dsc3+scip[2]*scale3i) +
                          0.5f*rr7*((gli[7]+gli[2])*psc5 +
                                    (glip[7]+glip[2])*dsc5 -
                          (sci[3]*scip[4]+scip[3]*sci[4])*scale5i) +
                          0.5f*rr9*(gli[3]*psc7+glip[3]*dsc7);

    gfi[2]              = -rr3*atomJ.q + rr5*sc[4] - rr7*sc[6];
    gfi[3]              = rr3*atomI.q + rr5*sc[3] + rr7*sc[5];
    gfi[4]              = 2.0f * rr5;
    gfi[5]              = rr7 * (sci[4]*psc7+scip[4]*dsc7);
    gfi[6]              = -rr7 * (sci[3]*psc7+scip[3]*dsc7);

    // get the induced force;
 
    ftm2i[1]            = gfi[1]*xr + 0.5f*
                          (- rr3*atomJ.q*(atomI.inducedDipoleS[0]*psc3+atomI.inducedDipolePS[0]*dsc3)
                           + rr5*sc[4]*(atomI.inducedDipoleS[0]*psc5+atomI.inducedDipolePS[0]*dsc5)
                           - rr7*sc[6]*(atomI.inducedDipoleS[0]*psc7+atomI.inducedDipolePS[0]*dsc7))
                           +(rr3*atomI.q*(atomJ.inducedDipoleS[0]*psc3+atomJ.inducedDipolePS[0]*dsc3)
                           + rr5*sc[3]*(atomJ.inducedDipoleS[0]*psc5+atomJ.inducedDipolePS[0]*dsc5)
                           + rr7*sc[5]*(atomJ.inducedDipoleS[0]*psc7+atomJ.inducedDipolePS[0]*dsc7))*0.5f
                           + rr5*scale5i*(sci[4]*atomI.inducedDipolePS[0]+scip[4]*atomI.inducedDipoleS[0]
                           + sci[3]*atomJ.inducedDipolePS[0]+scip[3]*atomJ.inducedDipoleS[0])*0.5f
                           + 0.5f*(sci[4]*psc5+scip[4]*dsc5)*rr5*atomI.labFrameDipole[0]
                           + 0.5f*(sci[3]*psc5+scip[3]*dsc5)*rr5*atomJ.labFrameDipole[0]
                           + 0.5f*gfi[4]*((qkui[1]-qiuk[1])*psc5
                           + (qkuip[1]-qiukp[1])*dsc5)
                           + gfi[5]*qir[1] + gfi[6]*qkr[1];
 
    ftm2i[2]            = gfi[1]*yr + 0.5f*
                          (- rr3*atomJ.q*(atomI.inducedDipoleS[1]*psc3+atomI.inducedDipolePS[1]*dsc3)
                           + rr5*sc[4]*(atomI.inducedDipoleS[1]*psc5+atomI.inducedDipolePS[1]*dsc5)
                           - rr7*sc[6]*(atomI.inducedDipoleS[1]*psc7+atomI.inducedDipolePS[1]*dsc7))
                           +(rr3*atomI.q*(atomJ.inducedDipoleS[1]*psc3+atomJ.inducedDipolePS[1]*dsc3)
                           + rr5*sc[3]*(atomJ.inducedDipoleS[1]*psc5+atomJ.inducedDipolePS[1]*dsc5)
                           + rr7*sc[5]*(atomJ.inducedDipoleS[1]*psc7+atomJ.inducedDipolePS[1]*dsc7))*0.5f
                           + rr5*scale5i*(sci[4]*atomI.inducedDipolePS[1]+scip[4]*atomI.inducedDipoleS[1]
                           + sci[3]*atomJ.inducedDipolePS[1]+scip[3]*atomJ.inducedDipoleS[1])*0.5f
                           + 0.5f*(sci[4]*psc5+scip[4]*dsc5)*rr5*atomI.labFrameDipole[1]
                           + 0.5f*(sci[3]*psc5+scip[3]*dsc5)*rr5*atomJ.labFrameDipole[1]
                           + 0.5f*gfi[4]*((qkui[2]-qiuk[2])*psc5
                           + (qkuip[2]-qiukp[2])*dsc5)
                           + gfi[5]*qir[2] + gfi[6]*qkr[2];

    ftm2i[3]            = gfi[1]*zr  + 0.5f*
                          (- rr3*atomJ.q*(atomI.inducedDipoleS[2]*psc3+atomI.inducedDipolePS[2]*dsc3)
                           + rr5*sc[4]*(atomI.inducedDipoleS[2]*psc5+atomI.inducedDipolePS[2]*dsc5)
                           - rr7*sc[6]*(atomI.inducedDipoleS[2]*psc7+atomI.inducedDipolePS[2]*dsc7))
                           +(rr3*atomI.q*(atomJ.inducedDipoleS[2]*psc3+atomJ.inducedDipolePS[2]*dsc3)
                           + rr5*sc[3]*(atomJ.inducedDipoleS[2]*psc5+atomJ.inducedDipolePS[2]*dsc5)
                           + rr7*sc[5]*(atomJ.inducedDipoleS[2]*psc7+atomJ.inducedDipolePS[2]*dsc7))*0.5f
                           + rr5*scale5i*(sci[4]*atomI.inducedDipolePS[2]+scip[4]*atomI.inducedDipoleS[2]
                           + sci[3]*atomJ.inducedDipolePS[2]+scip[3]*atomJ.inducedDipoleS[2])*0.5f
                           + 0.5f*(sci[4]*psc5+scip[4]*dsc5)*rr5*atomI.labFrameDipole[2]
                           + 0.5f*(sci[3]*psc5+scip[3]*dsc5)*rr5*atomJ.labFrameDipole[2]
                           + 0.5f*gfi[4]*((qkui[3]-qiuk[3])*psc5
                           + (qkuip[3]-qiukp[3])*dsc5)
                           + gfi[5]*qir[3] + gfi[6]*qkr[3];
 
    // intermediate values needed for partially excluded interactions

    fridmp[1]           = 0.5f * (rr3*((gli[1]+gli[6])*pscale
                          +(glip[1]+glip[6])*dscale)*ddsc3[1]
                          + rr5*((gli[2]+gli[7])*pscale
                          +(glip[2]+glip[7])*dscale)*ddsc5[1]
                          + rr7*(gli[3]*pscale+glip[3]*dscale)*ddsc7[1]);

    fridmp[2]           = 0.5f * (rr3*((gli[1]+gli[6])*pscale
                          +(glip[1]+glip[6])*dscale)*ddsc3[2]
                          + rr5*((gli[2]+gli[7])*pscale
                          +(glip[2]+glip[7])*dscale)*ddsc5[2]
                          + rr7*(gli[3]*pscale+glip[3]*dscale)*ddsc7[2]);

    fridmp[3]           = 0.5f * (rr3*((gli[1]+gli[6])*pscale
                          +(glip[1]+glip[6])*dscale)*ddsc3[3]
                          + rr5*((gli[2]+gli[7])*pscale
                          +(glip[2]+glip[7])*dscale)*ddsc5[3]
                          + rr7*(gli[3]*pscale+glip[3]*dscale)*ddsc7[3]);

    // get the induced-induced derivative terms
 
    findmp[1]           = 0.5f * uscale * (scip[2]*rr3*ddsc3[1]
                          - rr5*ddsc5[1]*(sci[3]*scip[4]+scip[3]*sci[4]));

    findmp[2]           = 0.5f * uscale * (scip[2]*rr3*ddsc3[2]
                          - rr5*ddsc5[2]*(sci[3]*scip[4]+scip[3]*sci[4]));

    findmp[3]           = 0.5f * uscale * (scip[2]*rr3*ddsc3[3]
                          - rr5*ddsc5[3]*(sci[3]*scip[4]+scip[3]*sci[4]));

    // handle of scaling for partially excluded interactions

    ftm2i[1]           -= fridmp[1] + findmp[1];
    ftm2i[2]           -= fridmp[2] + findmp[2];
    ftm2i[3]           -= fridmp[3] + findmp[3];

    // correction to convert mutual to direct polarization force

#if 0
               if (poltyp .eq. 'DIRECT') then;
                  gfd             = 0.5f * (rr5*scip[2]*scale3i;
     &                  - rr7*(scip[3]*sci[4]+sci[3]*scip[4])*scale5i);
                  fdir[1]             = gfd*xr + 0.5f*rr5*scale5i;
     &                         * (sci[4]*atomI.inducedDipolePS[0]+scip[4]*atomI.inducedDipoleS[0];
     &                           +sci[3]*atomJ.inducedDipolePS[0]+scip[3]*atomJ.inducedDipoleS[0]);
                  fdir[2]             = gfd*yr + 0.5f*rr5*scale5i;
     &                         * (sci[4]*atomI.inducedDipolePS[1]+scip[4]*atomI.inducedDipoleS[1];
     &                           +sci[3]*atomJ.inducedDipolePS[1]+scip[3]*atomJ.inducedDipoleS[1]);
                  fdir[3]             = gfd*zr + 0.5f*rr5*scale5i;
     &                         * (sci[4]*atomI.inducedDipolePS[2]+scip[4]*atomI.inducedDipoleS[2];
     &                           +sci[3]*atomJ.inducedDipolePS[2]+scip[3]*atomJ.inducedDipoleS[2]);
                  ftm2i[1]             = ftm2i[1] - fdir[1] + findmp[1];
                  ftm2i[2]             = ftm2i[2] - fdir[2] + findmp[2];
                  ftm2i[3]             = ftm2i[3] - fdir[3] + findmp[3];
               end if;
#endif

    // now perform the torque calculation
    // intermediate terms for torque between multipoles i and k
 
    gti[2]              = 0.5f * (sci[4]*psc5+scip[4]*dsc5) * rr5;
    gti[3]              = 0.5f * (sci[3]*psc5+scip[3]*dsc5) * rr5;
    gti[4]              = gfi[4];
    gti[5]              = gfi[5];
    gti[6]              = gfi[6];

    // calculate the induced torque components

    ttm2i[1]            = -rr3*(dixuk[1]*psc3+dixukp[1]*dsc3)*0.5f
                          + gti[2]*dixr[1] + gti[4]*((ukxqir[1]+rxqiuk[1])*psc5
                          +(ukxqirp[1]+rxqiukp[1])*dsc5)*0.5f - gti[5]*rxqir[1];

    ttm2i[2]            = -rr3*(dixuk[2]*psc3+dixukp[2]*dsc3)*0.5f
                          + gti[2]*dixr[2] + gti[4]*((ukxqir[2]+rxqiuk[2])*psc5
                          +(ukxqirp[2]+rxqiukp[2])*dsc5)*0.5f - gti[5]*rxqir[2];

    ttm2i[3]            = -rr3*(dixuk[3]*psc3+dixukp[3]*dsc3)*0.5f
                          + gti[2]*dixr[3] + gti[4]*((ukxqir[3]+rxqiuk[3])*psc5
                          +(ukxqirp[3]+rxqiukp[3])*dsc5)*0.5f - gti[5]*rxqir[3];

    ttm3i[1]            = -rr3*(dkxui[1]*psc3+dkxuip[1]*dsc3)*0.5f
                          + gti[3]*dkxr[1] - gti[4]*((uixqkr[1]+rxqkui[1])*psc5
                          +(uixqkrp[1]+rxqkuip[1])*dsc5)*0.5f - gti[6]*rxqkr[1];

    ttm3i[2]            = -rr3*(dkxui[2]*psc3+dkxuip[2]*dsc3)*0.5f
                          + gti[3]*dkxr[2] - gti[4]*((uixqkr[2]+rxqkui[2])*psc5
                          +(uixqkrp[2]+rxqkuip[2])*dsc5)*0.5f - gti[6]*rxqkr[2];

    ttm3i[3]            = -rr3*(dkxui[3]*psc3+dkxuip[3]*dsc3)*0.5f
                          + gti[3]*dkxr[3] - gti[4]*((uixqkr[3]+rxqkui[3])*psc5
                          +(uixqkrp[3]+rxqkuip[3])*dsc5)*0.5f - gti[6]*rxqkr[3];
 
    // update force and torque on site k
#if 0
    ftm1i(1,k)          = ftm1i(1,k) + ftm2i[1];
    ftm1i(2,k)          = ftm1i(2,k) + ftm2i[2];
    ftm1i(3,k)          = ftm1i(3,k) + ftm2i[3];

    ttm1i(1,k)          = ttm1i(1,k) + ttm3i[1];
    ttm1i(2,k)          = ttm1i(2,k) + ttm3i[2];
    ttm1i(3,k)          = ttm1i(3,k) + ttm3i[3];

    // update force and torque on site i

    ftm1i(1,i)          = ftm1i(1,i) - ftm2i[1];
    ftm1i(2,i)          = ftm1i(2,i) - ftm2i[2];
    ftm1i(3,i)          = ftm1i(3,i) - ftm2i[3];

    ttm1i(1,i)          = ttm1i(1,i) + ttm2i[1];
    ttm1i(2,i)          = ttm1i(2,i) + ttm2i[2];
    ttm1i(3,i)          = ttm1i(3,i) + ttm2i[3];
#endif
    outputForce[0]      = -ftm2i[1];
    outputForce[1]      = -ftm2i[2];
    outputForce[2]      = -ftm2i[3];

    outputTorqueI[0]    = ttm2i[1];
    outputTorqueI[1]    = ttm2i[2];
    outputTorqueI[2]    = ttm2i[3];

    outputTorqueJ[0]    = ttm3i[1];
    outputTorqueJ[1]    = ttm3i[2];
    outputTorqueJ[2]    = ttm3i[3];

    // construct auxiliary vectors for induced terms

    dixuk[1]            = atomI.labFrameDipole[1]*atomJ.inducedDipole[2] - atomI.labFrameDipole[2]*atomJ.inducedDipole[1];
    dixuk[2]            = atomI.labFrameDipole[2]*atomJ.inducedDipole[0] - atomI.labFrameDipole[0]*atomJ.inducedDipole[2];
    dixuk[3]            = atomI.labFrameDipole[0]*atomJ.inducedDipole[1] - atomI.labFrameDipole[1]*atomJ.inducedDipole[0];

    dkxui[1]            = atomJ.labFrameDipole[1]*atomI.inducedDipole[2] - atomJ.labFrameDipole[2]*atomI.inducedDipole[1];
    dkxui[2]            = atomJ.labFrameDipole[2]*atomI.inducedDipole[0] - atomJ.labFrameDipole[0]*atomI.inducedDipole[2];
    dkxui[3]            = atomJ.labFrameDipole[0]*atomI.inducedDipole[1] - atomJ.labFrameDipole[1]*atomI.inducedDipole[0];

    dixukp[1]           = atomI.labFrameDipole[1]*atomJ.inducedDipoleP[2] - atomI.labFrameDipole[2]*atomJ.inducedDipoleP[1];
    dixukp[2]           = atomI.labFrameDipole[2]*atomJ.inducedDipoleP[0] - atomI.labFrameDipole[0]*atomJ.inducedDipoleP[2];
    dixukp[3]           = atomI.labFrameDipole[0]*atomJ.inducedDipoleP[1] - atomI.labFrameDipole[1]*atomJ.inducedDipoleP[0];

    dkxuip[1]           = atomJ.labFrameDipole[1]*atomI.inducedDipoleP[2] - atomJ.labFrameDipole[2]*atomI.inducedDipoleP[1];
    dkxuip[2]           = atomJ.labFrameDipole[2]*atomI.inducedDipoleP[0] - atomJ.labFrameDipole[0]*atomI.inducedDipoleP[2];
    dkxuip[3]           = atomJ.labFrameDipole[0]*atomI.inducedDipoleP[1] - atomJ.labFrameDipole[1]*atomI.inducedDipoleP[0];

    qiuk[1]             = atomI.labFrameQuadrupole_XX*atomJ.inducedDipole[0] + atomI.labFrameQuadrupole_XY*atomJ.inducedDipole[1] + atomI.labFrameQuadrupole_XZ*atomJ.inducedDipole[2];
    qiuk[2]             = atomI.labFrameQuadrupole_XY*atomJ.inducedDipole[0] + atomI.labFrameQuadrupole_YY*atomJ.inducedDipole[1] + atomI.labFrameQuadrupole_YZ*atomJ.inducedDipole[2];
    qiuk[3]             = atomI.labFrameQuadrupole_XZ*atomJ.inducedDipole[0] + atomI.labFrameQuadrupole_YZ*atomJ.inducedDipole[1] + atomI.labFrameQuadrupole_ZZ*atomJ.inducedDipole[2];

    qkui[1]             = atomJ.labFrameQuadrupole_XX*atomI.inducedDipole[0] + atomJ.labFrameQuadrupole_XY*atomI.inducedDipole[1] + atomJ.labFrameQuadrupole_XZ*atomI.inducedDipole[2];
    qkui[2]             = atomJ.labFrameQuadrupole_XY*atomI.inducedDipole[0] + atomJ.labFrameQuadrupole_YY*atomI.inducedDipole[1] + atomJ.labFrameQuadrupole_YZ*atomI.inducedDipole[2];
    qkui[3]             = atomJ.labFrameQuadrupole_XZ*atomI.inducedDipole[0] + atomJ.labFrameQuadrupole_YZ*atomI.inducedDipole[1] + atomJ.labFrameQuadrupole_ZZ*atomI.inducedDipole[2];

    qiukp[1]            = atomI.labFrameQuadrupole_XX*atomJ.inducedDipoleP[0] + atomI.labFrameQuadrupole_XY*atomJ.inducedDipoleP[1] + atomI.labFrameQuadrupole_XZ*atomJ.inducedDipoleP[2];
    qiukp[2]            = atomI.labFrameQuadrupole_XY*atomJ.inducedDipoleP[0] + atomI.labFrameQuadrupole_YY*atomJ.inducedDipoleP[1] + atomI.labFrameQuadrupole_YZ*atomJ.inducedDipoleP[2];
    qiukp[3]            = atomI.labFrameQuadrupole_XZ*atomJ.inducedDipoleP[0] + atomI.labFrameQuadrupole_YZ*atomJ.inducedDipoleP[1] + atomI.labFrameQuadrupole_ZZ*atomJ.inducedDipoleP[2];

    qkuip[1]            = atomJ.labFrameQuadrupole_XX*atomI.inducedDipoleP[0] + atomJ.labFrameQuadrupole_XY*atomI.inducedDipoleP[1] + atomJ.labFrameQuadrupole_XZ*atomI.inducedDipoleP[2];
    qkuip[2]            = atomJ.labFrameQuadrupole_XY*atomI.inducedDipoleP[0] + atomJ.labFrameQuadrupole_YY*atomI.inducedDipoleP[1] + atomJ.labFrameQuadrupole_YZ*atomI.inducedDipoleP[2];
    qkuip[3]            = atomJ.labFrameQuadrupole_XZ*atomI.inducedDipoleP[0] + atomJ.labFrameQuadrupole_YZ*atomI.inducedDipoleP[1] + atomJ.labFrameQuadrupole_ZZ*atomI.inducedDipoleP[2];

    uixqkr[1]           = atomI.inducedDipole[1]*qkr[3] - atomI.inducedDipole[2]*qkr[2];
    uixqkr[2]           = atomI.inducedDipole[2]*qkr[1] - atomI.inducedDipole[0]*qkr[3];
    uixqkr[3]           = atomI.inducedDipole[0]*qkr[2] - atomI.inducedDipole[1]*qkr[1];

    ukxqir[1]           = atomJ.inducedDipole[1]*qir[3] - atomJ.inducedDipole[2]*qir[2];
    ukxqir[2]           = atomJ.inducedDipole[2]*qir[1] - atomJ.inducedDipole[0]*qir[3];
    ukxqir[3]           = atomJ.inducedDipole[0]*qir[2] - atomJ.inducedDipole[1]*qir[1];

    uixqkrp[1]          = atomI.inducedDipoleP[1]*qkr[3] - atomI.inducedDipoleP[2]*qkr[2];
    uixqkrp[2]          = atomI.inducedDipoleP[2]*qkr[1] - atomI.inducedDipoleP[0]*qkr[3];
    uixqkrp[3]          = atomI.inducedDipoleP[0]*qkr[2] - atomI.inducedDipoleP[1]*qkr[1];

    ukxqirp[1]          = atomJ.inducedDipoleP[1]*qir[3] - atomJ.inducedDipoleP[2]*qir[2];
    ukxqirp[2]          = atomJ.inducedDipoleP[2]*qir[1] - atomJ.inducedDipoleP[0]*qir[3];
    ukxqirp[3]          = atomJ.inducedDipoleP[0]*qir[2] - atomJ.inducedDipoleP[1]*qir[1];

    rxqiuk[1]           = yr*qiuk[3] - zr*qiuk[2];
    rxqiuk[2]           = zr*qiuk[1] - xr*qiuk[3];
    rxqiuk[3]           = xr*qiuk[2] - yr*qiuk[1];

    rxqkui[1]           = yr*qkui[3] - zr*qkui[2];
    rxqkui[2]           = zr*qkui[1] - xr*qkui[3];
    rxqkui[3]           = xr*qkui[2] - yr*qkui[1];

    rxqiukp[1]          = yr*qiukp[3] - zr*qiukp[2];
    rxqiukp[2]          = zr*qiukp[1] - xr*qiukp[3];
    rxqiukp[3]          = xr*qiukp[2] - yr*qiukp[1];

    rxqkuip[1]          = yr*qkuip[3] - zr*qkuip[2];
    rxqkuip[2]          = zr*qkuip[1] - xr*qkuip[3];
    rxqkuip[3]          = xr*qkuip[2] - yr*qkuip[1];

    // get intermediate variables for induction energy terms

    sci[1]              = atomI.inducedDipole[0]*atomJ.labFrameDipole[0] + atomI.inducedDipole[1]*atomJ.labFrameDipole[1]
                          + atomI.inducedDipole[2]*atomJ.labFrameDipole[2] + atomI.labFrameDipole[0]*atomJ.inducedDipole[0]
                          + atomI.labFrameDipole[1]*atomJ.inducedDipole[1] + atomI.labFrameDipole[2]*atomJ.inducedDipole[2];

    sci[2]              = atomI.inducedDipole[0]*atomJ.inducedDipole[0] + atomI.inducedDipole[1]*atomJ.inducedDipole[1] + atomI.inducedDipole[2]*atomJ.inducedDipole[2];

    sci[3]              = atomI.inducedDipole[0]*xr + atomI.inducedDipole[1]*yr + atomI.inducedDipole[2]*zr;
    sci[4]              = atomJ.inducedDipole[0]*xr + atomJ.inducedDipole[1]*yr + atomJ.inducedDipole[2]*zr;

    sci[7]              = qir[1]*atomJ.inducedDipole[0] + qir[2]*atomJ.inducedDipole[1] + qir[3]*atomJ.inducedDipole[2];
    sci[8]              = qkr[1]*atomI.inducedDipole[0] + qkr[2]*atomI.inducedDipole[1] + qkr[3]*atomI.inducedDipole[2];

    scip[1]             = atomI.inducedDipoleP[0]*atomJ.labFrameDipole[0] + atomI.inducedDipoleP[1]*atomJ.labFrameDipole[1] + atomI.inducedDipoleP[2]*atomJ.labFrameDipole[2] + atomI.labFrameDipole[0]*atomJ.inducedDipoleP[0] + atomI.labFrameDipole[1]*atomJ.inducedDipoleP[1] + atomI.labFrameDipole[2]*atomJ.inducedDipoleP[2];
    scip[2]             = atomI.inducedDipole[0]*atomJ.inducedDipoleP[0]+atomI.inducedDipole[1]*atomJ.inducedDipoleP[1] + atomI.inducedDipole[2]*atomJ.inducedDipoleP[2]+atomI.inducedDipoleP[0]*atomJ.inducedDipole[0] + atomI.inducedDipoleP[1]*atomJ.inducedDipole[1]+atomI.inducedDipoleP[2]*atomJ.inducedDipole[2];

    scip[3]             = atomI.inducedDipoleP[0]*xr + atomI.inducedDipoleP[1]*yr + atomI.inducedDipoleP[2]*zr;
    scip[4]             = atomJ.inducedDipoleP[0]*xr + atomJ.inducedDipoleP[1]*yr + atomJ.inducedDipoleP[2]*zr;

    scip[7]             = qir[1]*atomJ.inducedDipoleP[0] + qir[2]*atomJ.inducedDipoleP[1] + qir[3]*atomJ.inducedDipoleP[2];
    scip[8]             = qkr[1]*atomI.inducedDipoleP[0] + qkr[2]*atomI.inducedDipoleP[1] + qkr[3]*atomI.inducedDipoleP[2];

    // calculate the gl functions for potential energy

    gli[1]              = atomJ.q*sci[3] - atomI.q*sci[4];
    gli[2]              = -sc[3]*sci[4] - sci[3]*sc[4];
    gli[3]              = sci[3]*sc[6] - sci[4]*sc[5];
    gli[6]              = sci[1];
    gli[7]              = 2.0f * (sci[7]-sci[8]);

    glip[1]             = atomJ.q*scip[3] - atomI.q*scip[4];
    glip[2]             = -sc[3]*scip[4] - scip[3]*sc[4];
    glip[3]             = scip[3]*sc[6] - scip[4]*sc[5];
    glip[6]             = scip[1];
    glip[7]             = 2.0f * (scip[7]-scip[8]);

    // get the permanent multipole and induced energies

    *outputEnergy      += -0.5f * (rr3*(gli[1]+gli[6])*psc3 + rr5*(gli[2]+gli[7])*psc5 + rr7*gli[3]*psc7);

    // intermediate variables for the induced-permanent terms;

    gfi[1]             = 0.5f*rr5*((gli[1]+gli[6])*psc3 + (glip[1]+glip[6])*dsc3+scip[2]*scale3i) 
                         + 0.5f*rr7*((gli[7]+gli[2])*psc5 +(glip[7]+glip[2])*dsc5
                         -(sci[3]*scip[4]+scip[3]*sci[4])*scale5i)
                         + 0.5f*rr9*(gli[3]*psc7+glip[3]*dsc7);

    gfi[2]             = -rr3*atomJ.q + rr5*sc[4] - rr7*sc[6];
    gfi[3]             = rr3*atomI.q + rr5*sc[3] + rr7*sc[5];
    gfi[4]             = 2.0f * rr5;
    gfi[5]             = rr7 * (sci[4]*psc7+scip[4]*dsc7);
    gfi[6]             = -rr7 * (sci[3]*psc7+scip[3]*dsc7);

    // get the induced force

    ftm2i[1]           = gfi[1]*xr + 0.5f*
                        (- rr3*atomJ.q*(atomI.inducedDipole[0]*psc3+atomI.inducedDipoleP[0]*dsc3)
                         + rr5*sc[4]*(atomI.inducedDipole[0]*psc5+atomI.inducedDipoleP[0]*dsc5)
                         - rr7*sc[6]*(atomI.inducedDipole[0]*psc7+atomI.inducedDipoleP[0]*dsc7))
                         +(rr3*atomI.q*(atomJ.inducedDipole[0]*psc3+atomJ.inducedDipoleP[0]*dsc3)
                         + rr5*sc[3]*(atomJ.inducedDipole[0]*psc5+atomJ.inducedDipoleP[0]*dsc5)
                         + rr7*sc[5]*(atomJ.inducedDipole[0]*psc7+atomJ.inducedDipoleP[0]*dsc7))*0.5f
                         + rr5*scale5i*(sci[4]*atomI.inducedDipoleP[0]+scip[4]*atomI.inducedDipole[0]
                         + sci[3]*atomJ.inducedDipoleP[0]+scip[3]*atomJ.inducedDipole[0])*0.5f
                         + 0.5f*(sci[4]*psc5+scip[4]*dsc5)*rr5*atomI.labFrameDipole[0]
                         + 0.5f*(sci[3]*psc5+scip[3]*dsc5)*rr5*atomJ.labFrameDipole[0]
                         + 0.5f*gfi[4]*((qkui[1]-qiuk[1])*psc5
                         + (qkuip[1]-qiukp[1])*dsc5)
                         + gfi[5]*qir[1] + gfi[6]*qkr[1];

    ftm2i[2]           = gfi[1]*yr + 0.5f*
                        (- rr3*atomJ.q*(atomI.inducedDipole[1]*psc3+atomI.inducedDipoleP[1]*dsc3)
                         + rr5*sc[4]*(atomI.inducedDipole[1]*psc5+atomI.inducedDipoleP[1]*dsc5)
                         - rr7*sc[6]*(atomI.inducedDipole[1]*psc7+atomI.inducedDipoleP[1]*dsc7))
                         +(rr3*atomI.q*(atomJ.inducedDipole[1]*psc3+atomJ.inducedDipoleP[1]*dsc3)
                         + rr5*sc[3]*(atomJ.inducedDipole[1]*psc5+atomJ.inducedDipoleP[1]*dsc5)
                         + rr7*sc[5]*(atomJ.inducedDipole[1]*psc7+atomJ.inducedDipoleP[1]*dsc7))*0.5f
                         + rr5*scale5i*(sci[4]*atomI.inducedDipoleP[1]+scip[4]*atomI.inducedDipole[1]
                         + sci[3]*atomJ.inducedDipoleP[1]+scip[3]*atomJ.inducedDipole[1])*0.5f
                         + 0.5f*(sci[4]*psc5+scip[4]*dsc5)*rr5*atomI.labFrameDipole[1]
                         + 0.5f*(sci[3]*psc5+scip[3]*dsc5)*rr5*atomJ.labFrameDipole[1]
                         + 0.5f*gfi[4]*((qkui[2]-qiuk[2])*psc5
                         + (qkuip[2]-qiukp[2])*dsc5)
                         + gfi[5]*qir[2] + gfi[6]*qkr[2];

    ftm2i[3]           = gfi[1]*zr  + 0.5f*
                         (- rr3*atomJ.q*(atomI.inducedDipole[2]*psc3+atomI.inducedDipoleP[2]*dsc3)
                         + rr5*sc[4]*(atomI.inducedDipole[2]*psc5+atomI.inducedDipoleP[2]*dsc5)
                         - rr7*sc[6]*(atomI.inducedDipole[2]*psc7+atomI.inducedDipoleP[2]*dsc7))
                         +(rr3*atomI.q*(atomJ.inducedDipole[2]*psc3+atomJ.inducedDipoleP[2]*dsc3)
                         + rr5*sc[3]*(atomJ.inducedDipole[2]*psc5+atomJ.inducedDipoleP[2]*dsc5)
                         + rr7*sc[5]*(atomJ.inducedDipole[2]*psc7+atomJ.inducedDipoleP[2]*dsc7))*0.5f
                         + rr5*scale5i*(sci[4]*atomI.inducedDipoleP[2]+scip[4]*atomI.inducedDipole[2]
                         + sci[3]*atomJ.inducedDipoleP[2]+scip[3]*atomJ.inducedDipole[2])*0.5f
                         + 0.5f*(sci[4]*psc5+scip[4]*dsc5)*rr5*atomI.labFrameDipole[2]
                         + 0.5f*(sci[3]*psc5+scip[3]*dsc5)*rr5*atomJ.labFrameDipole[2]
                         + 0.5f*gfi[4]*((qkui[3]-qiuk[3])*psc5
                         + (qkuip[3]-qiukp[3])*dsc5)
                         + gfi[5]*qir[3] + gfi[6]*qkr[3];

    // intermediate values needed for partially excluded interactions

    fridmp[1]          = 0.5f * (rr3*((gli[1]+gli[6])*pscale
                         +(glip[1]+glip[6])*dscale)*ddsc3[1]
                         + rr5*((gli[2]+gli[7])*pscale
                         +(glip[2]+glip[7])*dscale)*ddsc5[1]
                         + rr7*(gli[3]*pscale+glip[3]*dscale)*ddsc7[1]);

    fridmp[2]          = 0.5f * (rr3*((gli[1]+gli[6])*pscale
                         +(glip[1]+glip[6])*dscale)*ddsc3[2]
                         + rr5*((gli[2]+gli[7])*pscale
                         +(glip[2]+glip[7])*dscale)*ddsc5[2]
                         + rr7*(gli[3]*pscale+glip[3]*dscale)*ddsc7[2]);

    fridmp[3]          = 0.5f * (rr3*((gli[1]+gli[6])*pscale
                         +(glip[1]+glip[6])*dscale)*ddsc3[3]
                         + rr5*((gli[2]+gli[7])*pscale
                         +(glip[2]+glip[7])*dscale)*ddsc5[3]
                         + rr7*(gli[3]*pscale+glip[3]*dscale)*ddsc7[3]);

    // get the induced-induced derivative terms;

    findmp[1]          = 0.5f * uscale * (scip[2]*rr3*ddsc3[1]
                         - rr5*ddsc5[1]*(sci[3]*scip[4]+scip[3]*sci[4]));

    findmp[2]          = 0.5f * uscale * (scip[2]*rr3*ddsc3[2]
                         - rr5*ddsc5[2]*(sci[3]*scip[4]+scip[3]*sci[4]));

    findmp[3]          = 0.5f * uscale * (scip[2]*rr3*ddsc3[3]
                         - rr5*ddsc5[3]*(sci[3]*scip[4]+scip[3]*sci[4]));

    // handle of scaling for partially excluded interactions;

    ftm2i[1]           = ftm2i[1] - fridmp[1] - findmp[1];
    ftm2i[2]           = ftm2i[2] - fridmp[2] - findmp[2];
    ftm2i[3]           = ftm2i[3] - fridmp[3] - findmp[3];

    // correction to convert mutual to direct polarization force;

#if 0
               if (poltyp .eq. 'DIRECT') then;
                  gfd = 0.5f * (rr5*scip[2]*scale3i;
     &                  - rr7*(scip[3]*sci[4]+sci[3]*scip[4])*scale5i);
                  fdir[1] = gfd*xr + 0.5f*rr5*scale5i;
     &                         * (sci[4]*atomI.inducedDipoleP[0]+scip[4]*atomI.inducedDipole[0];
     &                           +sci[3]*atomJ.inducedDipoleP[0]+scip[3]*atomJ.inducedDipole[0]);
                  fdir[2] = gfd*yr + 0.5f*rr5*scale5i;
     &                         * (sci[4]*atomI.inducedDipoleP[1]+scip[4]*atomI.inducedDipole[1];
     &                           +sci[3]*atomJ.inducedDipoleP[1]+scip[3]*atomJ.inducedDipole[1]);
                  fdir[3] = gfd*zr + 0.5f*rr5*scale5i;
     &                         * (sci[4]*atomI.inducedDipoleP[2]+scip[4]*atomI.inducedDipole[2];
     &                           +sci[3]*atomJ.inducedDipoleP[2]+scip[3]*atomJ.inducedDipole[2]);
                  ftm2i[1] = ftm2i[1] - fdir[1] + findmp[1];
                  ftm2i[2] = ftm2i[2] - fdir[2] + findmp[2];
                  ftm2i[3] = ftm2i[3] - fdir[3] + findmp[3];
               end if;
#endif

    // now perform the torque calculation
    // intermediate terms for torque between multipoles i and k

    gti[2]             = 0.5f * (sci[4]*psc5+scip[4]*dsc5) * rr5;
    gti[3]             = 0.5f * (sci[3]*psc5+scip[3]*dsc5) * rr5;
    gti[4]             = gfi[4];
    gti[5]             = gfi[5];
    gti[6]             = gfi[6];

    // calculate the induced torque components

    ttm2i[1]           = -rr3*(dixuk[1]*psc3+dixukp[1]*dsc3)*0.5f
                         + gti[2]*dixr[1] + gti[4]*((ukxqir[1]+rxqiuk[1])*psc5
                         +(ukxqirp[1]+rxqiukp[1])*dsc5)*0.5f - gti[5]*rxqir[1];

    ttm2i[2]           = -rr3*(dixuk[2]*psc3+dixukp[2]*dsc3)*0.5f
                         + gti[2]*dixr[2] + gti[4]*((ukxqir[2]+rxqiuk[2])*psc5
                         +(ukxqirp[2]+rxqiukp[2])*dsc5)*0.5f - gti[5]*rxqir[2];

    ttm2i[3]           = -rr3*(dixuk[3]*psc3+dixukp[3]*dsc3)*0.5f
                         + gti[2]*dixr[3] + gti[4]*((ukxqir[3]+rxqiuk[3])*psc5
                         +(ukxqirp[3]+rxqiukp[3])*dsc5)*0.5f - gti[5]*rxqir[3];

    ttm3i[1]           = -rr3*(dkxui[1]*psc3+dkxuip[1]*dsc3)*0.5f
                         + gti[3]*dkxr[1] - gti[4]*((uixqkr[1]+rxqkui[1])*psc5
                         +(uixqkrp[1]+rxqkuip[1])*dsc5)*0.5f - gti[6]*rxqkr[1];

    ttm3i[2]           = -rr3*(dkxui[2]*psc3+dkxuip[2]*dsc3)*0.5f
                         + gti[3]*dkxr[2] - gti[4]*((uixqkr[2]+rxqkui[2])*psc5
                         +(uixqkrp[2]+rxqkuip[2])*dsc5)*0.5f - gti[6]*rxqkr[2];

    ttm3i[3]           = -rr3*(dkxui[3]*psc3+dkxuip[3]*dsc3)*0.5f
                         + gti[3]*dkxr[3] - gti[4]*((uixqkr[3]+rxqkui[3])*psc5
                         +(uixqkrp[3]+rxqkuip[3])*dsc5)*0.5f - gti[6]*rxqkr[3];

    // update force and torque on site k;

#if 0
    ftm1i(1,k) = ftm1i(1,k) - ftm2i[1];
    ftm1i(2,k) = ftm1i(2,k) - ftm2i[2];
    ftm1i(3,k) = ftm1i(3,k) - ftm2i[3];
    ttm1i(1,k) = ttm1i(1,k) - ttm3i[1];
    ttm1i(2,k) = ttm1i(2,k) - ttm3i[2];
    ttm1i(3,k) = ttm1i(3,k) - ttm3i[3];

    // update force and torque on site i

    ftm1i(1,i) = ftm1i(1,i) + ftm2i[1];
    ftm1i(2,i) = ftm1i(2,i) + ftm2i[2];
    ftm1i(3,i) = ftm1i(3,i) + ftm2i[3];
    ttm1i(1,i) = ttm1i(1,i) - ttm2i[1];
    ttm1i(2,i) = ttm1i(2,i) - ttm2i[2];
    ttm1i(3,i) = ttm1i(3,i) - ttm2i[3];
#endif
 
    outputForce[0]     += ftm2i[1];
    outputForce[1]     += ftm2i[2];
    outputForce[2]     += ftm2i[3];

    outputTorqueI[0]   -= ttm2i[1];
    outputTorqueI[1]   -= ttm2i[2];
    outputTorqueI[2]   -= ttm2i[3];

    outputTorqueJ[0]   -= ttm3i[1];
    outputTorqueJ[1]   -= ttm3i[2];
    outputTorqueJ[2]   -= ttm3i[3];

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
        (void) fprintf( amoebaGpu->log, "%s %d maxCovalentDegreeSz=%d"
                        " gamma=%.3e scalingDistanceCutoff=%.3f ZZZ\n",
                        methodName, gpu->natoms,
                        amoebaGpu->maxCovalentDegreeSz, amoebaGpu->pGamma,
                        amoebaGpu->scalingDistanceCutoff );
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
            maxThreads = 192;
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

