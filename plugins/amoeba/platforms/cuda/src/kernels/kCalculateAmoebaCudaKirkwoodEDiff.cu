//-----------------------------------------------------------------------------------------

//-----------------------------------------------------------------------------------------

#include "amoebaGpuTypes.h"
#include "amoebaCudaKernels.h"
#include "kCalculateAmoebaCudaUtilities.h"
//#include "kCalculateAmoebaCudaKirkwoodParticle.h"
#include "kCalculateAmoebaCudaKirkwoodEDiffParticle.h"

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

__device__ void loadKirkwoodEDiffShared( struct KirkwoodEDiffParticle* sA, unsigned int atomI )
{
    // coordinates & charge

    sA->x                        = cSim.pPosq[atomI].x;
    sA->y                        = cSim.pPosq[atomI].y;
    sA->z                        = cSim.pPosq[atomI].z;
    sA->q                        = cSim.pPosq[atomI].w;

    sA->damp                     = cAmoebaSim.pDampingFactorAndThole[atomI].x;
    sA->thole                    = cAmoebaSim.pDampingFactorAndThole[atomI].y;

    // lab dipole

    sA->labFrameDipole[0]        = cAmoebaSim.pLabFrameDipole[atomI*3];
    sA->labFrameDipole[1]        = cAmoebaSim.pLabFrameDipole[atomI*3+1];
    sA->labFrameDipole[2]        = cAmoebaSim.pLabFrameDipole[atomI*3+2];

    // lab quadrupole

    sA->labFrameQuadrupole_XX    = cAmoebaSim.pLabFrameQuadrupole[atomI*9];
    sA->labFrameQuadrupole_XY    = cAmoebaSim.pLabFrameQuadrupole[atomI*9+1];
    sA->labFrameQuadrupole_XZ    = cAmoebaSim.pLabFrameQuadrupole[atomI*9+2];
    sA->labFrameQuadrupole_YY    = cAmoebaSim.pLabFrameQuadrupole[atomI*9+4];
    sA->labFrameQuadrupole_YZ    = cAmoebaSim.pLabFrameQuadrupole[atomI*9+5];
    sA->labFrameQuadrupole_ZZ    = cAmoebaSim.pLabFrameQuadrupole[atomI*9+8];

    // induced dipole

    sA->inducedDipole[0]         = cAmoebaSim.pInducedDipole[atomI*3];
    sA->inducedDipole[1]         = cAmoebaSim.pInducedDipole[atomI*3+1];
    sA->inducedDipole[2]         = cAmoebaSim.pInducedDipole[atomI*3+2];

    // induced dipole polar

    sA->inducedDipoleP[0]        = cAmoebaSim.pInducedDipolePolar[atomI*3];
    sA->inducedDipoleP[1]        = cAmoebaSim.pInducedDipolePolar[atomI*3+1];
    sA->inducedDipoleP[2]        = cAmoebaSim.pInducedDipolePolar[atomI*3+2];

    // induced dipole

    sA->inducedDipoleS[0]        = cAmoebaSim.pInducedDipoleS[atomI*3];
    sA->inducedDipoleS[1]        = cAmoebaSim.pInducedDipoleS[atomI*3+1];
    sA->inducedDipoleS[2]        = cAmoebaSim.pInducedDipoleS[atomI*3+2];

    // induced dipole polar

    sA->inducedDipolePS[0]       = cAmoebaSim.pInducedDipolePolarS[atomI*3];
    sA->inducedDipolePS[1]       = cAmoebaSim.pInducedDipolePolarS[atomI*3+1];
    sA->inducedDipolePS[2]       = cAmoebaSim.pInducedDipolePolarS[atomI*3+2];

}

/*****************************************************************************

       ediff1 correct vacuum to SCRF derivatives

       calculates the energy and derivatives of polarizing
       the vacuum induced dipoles to their SCRF polarized values

******************************************************************************/
#ifdef INCLUDE_TORQUE
__device__ void calculateKirkwoodEDiffPairIxnOrig_kernel( KirkwoodEDiffParticle& atomI,  KirkwoodEDiffParticle& atomJ,
                                                      float pscale,                  float dscale,
                                                      float*  outputEnergy,          float*  outputForce,
                                                      float*  outputTorqueI,         float* outputTorqueJ){

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

    // set conversion factor, cutoff and scaling coefficients

    //f           = cAmoebaSim.electric / cAmoebaSim.dwater;

    // deltaR

    float xr          = atomJ.x - atomI.x;
    float yr          = atomJ.y - atomI.y;
    float zr          = atomJ.z - atomI.z;

    float r2          = xr*xr + yr*yr + zr*zr;

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

    float rxqir1            = yr*qir3 - zr*qir2;
    float rxqir2            = zr*qir1 - xr*qir3;
    float rxqir3            = xr*qir2 - yr*qir1;

    float rxqkr1            = yr*qkr3 - zr*qkr2;
    float rxqkr2            = zr*qkr1 - xr*qkr3;
    float rxqkr3            = xr*qkr2 - yr*qkr1;

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

    *outputEnergy           = 0.5f * (rr3*(gli1+gli6)*psc3 + rr5*(gli2+gli7)*psc5 + rr7*gli3*psc7);

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

    if ( cAmoebaSim.polarizationType ){
        float gfd      = 0.5f * (rr5*scip2*scale3i - rr7*(scip3*sci4+sci3*scip4)*scale5i);
        float fdir1    = gfd*xr + 0.5f*rr5*scale5i* (sci4*atomI.inducedDipolePS[0]+scip4*atomI.inducedDipoleS[0] + sci3*atomJ.inducedDipolePS[0]+scip3*atomJ.inducedDipoleS[0]);
        float fdir2    = gfd*yr + 0.5f*rr5*scale5i* (sci4*atomI.inducedDipolePS[1]+scip4*atomI.inducedDipoleS[1] + sci3*atomJ.inducedDipolePS[1]+scip3*atomJ.inducedDipoleS[1]);
        float fdir3    = gfd*zr + 0.5f*rr5*scale5i* (sci4*atomI.inducedDipolePS[2]+scip4*atomI.inducedDipoleS[2] + sci3*atomJ.inducedDipolePS[2]+scip3*atomJ.inducedDipoleS[2]);
        ftm2i1        -= fdir1 - findmp1;
        ftm2i2        -= fdir2 - findmp2;
        ftm2i3        -= fdir3 - findmp3;
     
    }

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
    //
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

    if ( cAmoebaSim.polarizationType ){

        float gfd    = 0.5f * (rr5*scip2*scale3i- rr7*(scip3*sci4+sci3*scip4)*scale5i);
        float fdir1  = gfd*xr + 0.5f*rr5*scale5i* (sci4*atomI.inducedDipoleP[0]+scip4*atomI.inducedDipole[0] + sci3*atomJ.inducedDipoleP[0]+scip3*atomJ.inducedDipole[0]);
        float fdir2  = gfd*yr + 0.5f*rr5*scale5i* (sci4*atomI.inducedDipoleP[1]+scip4*atomI.inducedDipole[1] + sci3*atomJ.inducedDipoleP[1]+scip3*atomJ.inducedDipole[1]);
        float fdir3  = gfd*zr + 0.5f*rr5*scale5i* (sci4*atomI.inducedDipoleP[2]+scip4*atomI.inducedDipole[2] + sci3*atomJ.inducedDipoleP[2]+scip3*atomJ.inducedDipole[2]);
        ftm2i1      -= fdir1 - findmp1;
        ftm2i2      -= fdir2 - findmp2;
        ftm2i3      -= fdir3 - findmp3;
    }

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
#endif

#undef SUB_METHOD_NAME
#undef F1
#define SUB_METHOD_NAME(a, b) a##F1##b
#define F1
#include "kCalculateAmoebaCudaKirkwoodEDiff_b.h"
#undef F1
#undef SUB_METHOD_NAME

#undef SUB_METHOD_NAME
#undef T1
#define SUB_METHOD_NAME(a, b) a##T1##b
#define T1
#include "kCalculateAmoebaCudaKirkwoodEDiff_b.h"
#undef T1
#undef SUB_METHOD_NAME

#undef SUB_METHOD_NAME
#undef T3
#define SUB_METHOD_NAME(a, b) a##T3##b
#define T3
#include "kCalculateAmoebaCudaKirkwoodEDiff_b.h"
#undef T3
#undef SUB_METHOD_NAME

#define APPLY_SCALE
#undef SUB_METHOD_NAME
#undef F1
#define SUB_METHOD_NAME(a, b) a##F1Scale##b
#define F1
#include "kCalculateAmoebaCudaKirkwoodEDiff_b.h"
#undef F1
#undef SUB_METHOD_NAME

#undef SUB_METHOD_NAME
#undef T1
#define SUB_METHOD_NAME(a, b) a##T1Scale##b
#define T1
#include "kCalculateAmoebaCudaKirkwoodEDiff_b.h"
#undef T1
#undef SUB_METHOD_NAME

#undef SUB_METHOD_NAME
#undef T3
#define SUB_METHOD_NAME(a, b) a##T3Scale##b
#define T3
#include "kCalculateAmoebaCudaKirkwoodEDiff_b.h"
#undef T3
#undef SUB_METHOD_NAME
#undef APPLY_SCALE

#ifdef INCLUDE_TORQUE

/*****************************************************************************
 *
 *        ediff1 correct vacuum to SCRF derivatives
 *
 *               calculates the energy and derivatives of polarizing
 *                      the vacuum induced dipoles to their SCRF polarized values
 *
*******************************************************************************/

__device__ void calculateKirkwoodEDiffPairIxn_kernel( KirkwoodEDiffParticle& atomI,  KirkwoodEDiffParticle& atomJ,
                                                      float pscale,                  float dscale,
                                                      float*  outputEnergy, float forceFactor ){

#ifdef Orig
    return calculateKirkwoodEDiffPairIxnOrig_kernel( atomI,  atomJ, pscale, dscale, outputEnergy,outputForce,
                                                     outputTorqueI, outputTorqueJ );
#else

    float force[3];
    float energy;
    calculateKirkwoodEDiffPairIxnF1_kernel( atomI,  atomJ, pscale, dscale, &energy, force);
    atomI.force[0] += force[0];
    atomI.force[1] += force[1];
    atomI.force[2] += force[2];
    if( forceFactor == 1.0f ){
        atomJ.force[0] -= force[0];
        atomJ.force[1] -= force[1];
        atomJ.force[2] -= force[2];
        energy *= 0.5f;
    } else {
        energy *= 0.25f;
    }
    *outputEnergy += energy;

    calculateKirkwoodEDiffPairIxnT1_kernel( atomI,  atomJ, pscale, dscale, outputEnergy, force);
    atomI.torque[0] += force[0];
    atomI.torque[1] += force[1];
    atomI.torque[2] += force[2];
    //calculateKirkwoodEDiffPairIxnT1_kernel( atomJ,  atomI, pscale, dscale, outputEnergy, force );
    calculateKirkwoodEDiffPairIxnT3_kernel( atomI,  atomJ, pscale, dscale, outputEnergy, force );
    atomJ.torque[0] += force[0];
    atomJ.torque[1] += force[1];
    atomJ.torque[2] += force[2];

    return;

#endif

}

#endif

#ifdef AMOEBA_DEBUG
__device__ static int debugAccumulate( unsigned int index, float4* debugArray, float* field, unsigned int addMask, float idLabel )
{
    index                             += cSim.paddedNumberOfAtoms;
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

#ifdef INCLUDE_TORQUE
    sA->torque[0]             = 0.0f;
    sA->torque[1]             = 0.0f;
    sA->torque[2]             = 0.0f;
#endif

}

// Include versions of the kernels for N^2 calculations.

#undef USE_OUTPUT_BUFFER_PER_WARP
#define METHOD_NAME(a, b) a##N2##b
#include "kCalculateAmoebaCudaKirkwoodEDiff.h"
#define USE_OUTPUT_BUFFER_PER_WARP
#undef METHOD_NAME
#define METHOD_NAME(a, b) a##N2ByWarp##b
#include "kCalculateAmoebaCudaKirkwoodEDiff.h"

// reduce psWorkArray_3_1 -> torque

static void kReduceTorque( amoebaGpuContext amoebaGpu )
{

    gpuContext gpu = amoebaGpu->gpuContext;
    kReduceFields_kernel<<<gpu->sim.nonbond_blocks, gpu->sim.bsf_reduce_threads_per_block>>>(
                           gpu->sim.paddedNumberOfAtoms*3, gpu->sim.outputBuffers,
                           amoebaGpu->psWorkArray_3_1->_pDevData, amoebaGpu->psTorque->_pDevData, 0 );

    LAUNCHERROR("kReduceForceTorqueKirkwoodEDiff");
}

/**---------------------------------------------------------------------------------------

   Compute Amoeba electrostatic force & torque

   @param amoebaGpu        amoebaGpu context
   @param gpu              OpenMM gpu Cuda context

   --------------------------------------------------------------------------------------- */

void kCalculateAmoebaKirkwoodEDiff( amoebaGpuContext amoebaGpu )
{
  
   // ---------------------------------------------------------------------------------------

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
    //
    static unsigned int threadsPerBlock = 0;
    if( threadsPerBlock == 0 ){
        unsigned int maxThreads;
        if (gpu->sm_version >= SM_20)
            //maxThreads = 384;
            maxThreads = 512;
        else if (gpu->sm_version >= SM_12)
            maxThreads = 96;
        else
            maxThreads = 32;
        threadsPerBlock = std::min(getThreadsPerBlock( amoebaGpu, sizeof(KirkwoodEDiffParticle), gpu->sharedMemoryPerBlock ), maxThreads);
    }   
    
#ifdef AMOEBA_DEBUG
    if( amoebaGpu->log ){
        (void) fprintf( amoebaGpu->log, "kCalculateAmoebaCudaKirkwoodEDiffN2Forces: blocks=%u threads=%u bffr/Warp=%u atm=%lu shrd=%lu"
                                        " ixnCt=%lu workUnits=%u sm=%d device=%d sharedMemoryPerBlock=%u step=%d\n",
                        gpu->sim.nonbond_blocks, threadsPerBlock, gpu->bOutputBufferPerWarp,
                        sizeof(KirkwoodEDiffParticle), sizeof(KirkwoodEDiffParticle)*threadsPerBlock,
                        (*gpu->psInteractionCount)[0], gpu->sim.workUnits, gpu->sm_version, gpu->device, gpu->sharedMemoryPerBlock, timestep );
        (void) fflush( amoebaGpu->log );
    }   
#endif

    if (gpu->bOutputBufferPerWarp){
        kCalculateAmoebaCudaKirkwoodEDiffN2ByWarpForces_kernel<<<gpu->sim.nonbond_blocks, threadsPerBlock, sizeof(KirkwoodEDiffParticle)*threadsPerBlock>>>(
                                                                           amoebaGpu->psWorkUnit->_pDevData, amoebaGpu->psWorkArray_3_1->_pDevData );

    } else {

        kCalculateAmoebaCudaKirkwoodEDiffN2Forces_kernel<<<gpu->sim.nonbond_blocks, threadsPerBlock, sizeof(KirkwoodEDiffParticle)*threadsPerBlock>>>(
                                                                           amoebaGpu->psWorkUnit->_pDevData, amoebaGpu->psWorkArray_3_1->_pDevData );
    }
    LAUNCHERROR("kCalculateAmoebaCudaKirkwoodEDiffN2Forces");

    // reduce and map torques to forces

    kReduceTorque( amoebaGpu );
    cudaComputeAmoebaMapTorqueAndAddToForce( amoebaGpu, amoebaGpu->psTorque );

   // ---------------------------------------------------------------------------------------
}
