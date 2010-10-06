//-----------------------------------------------------------------------------------------

//-----------------------------------------------------------------------------------------

#include "amoebaGpuTypes.h"
#include "amoebaCudaKernels.h"
#include "kCalculateAmoebaCudaUtilities.h"

//#define AMOEBA_DEBUG

static __constant__ cudaGmxSimulation cSim;
static __constant__ cudaAmoebaGmxSimulation cAmoebaSim;

void SetCalculateAmoebaPmeDirectElectrostaticSim(amoebaGpuContext amoebaGpu)
{
    cudaError_t status;
    gpuContext gpu = amoebaGpu->gpuContext;
    status         = cudaMemcpyToSymbol(cSim, &gpu->sim, sizeof(cudaGmxSimulation));    
    RTERROR(status, "SetCalculateAmoebaPmeDirectElectrostaticSim: cudaMemcpyToSymbol: SetSim copy to cSim failed");
    status         = cudaMemcpyToSymbol(cAmoebaSim, &amoebaGpu->amoebaSim, sizeof(cudaAmoebaGmxSimulation));    
    RTERROR(status, "SetCalculateAmoebaPmeDirectElectrostaticSim: cudaMemcpyToSymbol: SetSim copy to cAmoebaSim failed");
}

void GetCalculateAmoebaPmeDirectElectrostaticSim(amoebaGpuContext amoebaGpu)
{
    cudaError_t status;
    gpuContext gpu = amoebaGpu->gpuContext;
    status = cudaMemcpyFromSymbol(&gpu->sim, cSim, sizeof(cudaGmxSimulation));    
    RTERROR(status, "GetCalculateAmoebaPmeDirectElectrostaticSim: cudaMemcpyFromSymbol: SetSim copy from cSim failed");
    status = cudaMemcpyFromSymbol(&amoebaGpu->amoebaSim, cAmoebaSim, sizeof(cudaAmoebaGmxSimulation));    
    RTERROR(status, "GetCalculateAmoebaPmeDirectElectrostaticSim: cudaMemcpyFromSymbol: SetSim copy from cAmoebaSim failed");
}

static int const PScaleIndex            =  0; 
static int const DScaleIndex            =  1; 
static int const UScaleIndex            =  2; 
static int const MScaleIndex            =  3;
static int const LastScalingIndex       =  4;

struct PmeDirectElectrostaticParticle {

    // coordinates charge

    float x;
    float y;
    float z;
    float q;

    // lab frame dipole

    float labFrameDipole[3];

    // lab frame quadrupole

    float labFrameQuadrupole[9];

    // induced dipole

    float inducedDipole[3];

    // polar induced dipole

    float inducedDipoleP[3];

    // scaling factors

    float thole;
    float damp;

    float force[3];

    float torque[3];
    float padding;

    float tempForce[3];
    float tempTorque[3];
};

__device__ void sumTempBuffer( PmeDirectElectrostaticParticle& atomI, PmeDirectElectrostaticParticle& atomJ ){

    atomI.tempForce[0]  += atomJ.tempForce[0];
    atomI.tempForce[1]  += atomJ.tempForce[1];
    atomI.tempForce[2]  += atomJ.tempForce[2];

    atomI.tempTorque[0] += atomJ.tempTorque[0];
    atomI.tempTorque[1] += atomJ.tempTorque[1];
    atomI.tempTorque[2] += atomJ.tempTorque[2];
}

/*
__device__ static void debugSetup( unsigned int atomI, unsigned int atomJ,
                                   float4* debugArray, float4* pullBack )
{
    unsigned int index                 = atomI + atomJ*cAmoebaSim.paddedNumberOfAtoms;
    float blockId                      = 111.0f;

    debugArray[index].x                = (float) atomI;
    debugArray[index].y                = (float) atomJ;
    debugArray[index].z                = 0.0f;
    debugArray[index].w                = blockId;

    for( int pullIndex = 0; pullIndex < 1; pullIndex++ ){
        index                             += cAmoebaSim.paddedNumberOfAtoms;
        debugArray[index].x                = pullBack[pullIndex].x;
        debugArray[index].y                = pullBack[pullIndex].y;
        debugArray[index].z                = pullBack[pullIndex].z;
        debugArray[index].w                = pullBack[pullIndex].w;
    }   
}
*/

// self-energy for PME

__device__ static void calculatePmeSelfEnergyElectrostaticPairIxn_kernel( PmeDirectElectrostaticParticle& atomI, float* energy)
{
    float term   = 2.0f*cSim.alphaEwald*cSim.alphaEwald;
    float fterm  = -(cAmoebaSim.electric/cAmoebaSim.dielec)*cSim.alphaEwald/cAmoebaSim.sqrtPi;

    float cii    = atomI.q*atomI.q;

    float dii    = atomI.labFrameDipole[0]*atomI.labFrameDipole[0] +
                   atomI.labFrameDipole[1]*atomI.labFrameDipole[1] +
                   atomI.labFrameDipole[2]*atomI.labFrameDipole[2];

    float qii    = atomI.labFrameQuadrupole[0]*atomI.labFrameQuadrupole[0] +
                   atomI.labFrameQuadrupole[4]*atomI.labFrameQuadrupole[4] +
                   atomI.labFrameQuadrupole[8]*atomI.labFrameQuadrupole[8] + 2.0f*(
                   atomI.labFrameQuadrupole[1]*atomI.labFrameQuadrupole[1] +
                   atomI.labFrameQuadrupole[2]*atomI.labFrameQuadrupole[2] +
                   atomI.labFrameQuadrupole[5]*atomI.labFrameQuadrupole[5]);

    float uii    = atomI.labFrameDipole[0]*atomI.inducedDipole[0] + atomI.labFrameDipole[1]*atomI.inducedDipole[1] + atomI.labFrameDipole[2]*atomI.inducedDipole[2];

    *energy      = (cii + term*(dii/3.0f + 2.0f*term*qii/5.0f));
    *energy     += term*uii/3.0f;
    *energy     *= fterm;
}

// self-torque for PME

__device__ static void calculatePmeSelfTorqueElectrostaticPairIxn_kernel( PmeDirectElectrostaticParticle& atomI)
{
    float term       = (4.0f/3.0f)*(cAmoebaSim.electric/cAmoebaSim.dielec)*(cSim.alphaEwald*cSim.alphaEwald*cSim.alphaEwald)/cAmoebaSim.sqrtPi;

    float uix        = 0.5f*(atomI.inducedDipole[0] + atomI.inducedDipoleP[0]);
    float uiy        = 0.5f*(atomI.inducedDipole[1] + atomI.inducedDipoleP[1]);
    float uiz        = 0.5f*(atomI.inducedDipole[2] + atomI.inducedDipoleP[2]);

    atomI.torque[0] += term*(atomI.labFrameDipole[1]*uiz - atomI.labFrameDipole[2]*uiy);
    atomI.torque[1] += term*(atomI.labFrameDipole[2]*uix - atomI.labFrameDipole[0]*uiz);
    atomI.torque[2] += term*(atomI.labFrameDipole[0]*uiy - atomI.labFrameDipole[1]*uix);
}

__device__ void calculatePmeDirectElectrostaticPairIxn_kernel( PmeDirectElectrostaticParticle& atomI,   PmeDirectElectrostaticParticle& atomJ,
                                                               float* scalingFactors, float*  outputForce, float  outputTorque[2][3], float* energy
#ifdef AMOEBA_DEBUG
                                                               ,float4* debugArray 
#endif
 ){
  

    float e,ei;
    float erl,erli;
    float di[4],qi[10];
    float dk[4],qk[10];
    float fridmp[4],findmp[4];
    float ftm2[4],ftm2i[4];
    float ftm2r[4],ftm2ri[4];
    float ttm2[4],ttm3[4];
    float ttm2i[4],ttm3i[4];
    float ttm2r[4],ttm3r[4];
    float ttm2ri[4],ttm3ri[4];
    float dixdk[4];
    float dkxui[4],dixuk[4];
    float dixukp[4],dkxuip[4];
    float uixqkr[4],ukxqir[4];
    float uixqkrp[4],ukxqirp[4];
    float qiuk[4],qkui[4];
    float qiukp[4],qkuip[4];
    float rxqiuk[4],rxqkui[4];
    float rxqiukp[4],rxqkuip[4];
    float qidk[4],qkdi[4];
    float qir[4],qkr[4];
    float qiqkr[4],qkqir[4];
    float qixqk[4],rxqir[4];
    float dixr[4],dkxr[4];
    float dixqkr[4],dkxqir[4];
    float rxqkr[4],qkrxqir[4];
    float rxqikr[4],rxqkir[4];
    float rxqidk[4],rxqkdi[4];
    float ddsc3[4],ddsc5[4];
    float ddsc7[4];
    float bn[6];
    float sc[11],gl[9];
    float sci[9],scip[9];
    float gli[8],glip[8];
    float gf[8],gfi[7];
    float gfr[8],gfri[7];
    float gti[7],gtri[7];

    float conversionFactor   = (-cAmoebaSim.electric/cAmoebaSim.dielec);

    // set the permanent multipole and induced dipole values;

    float pdi   = atomI.damp;
    float pti   = atomI.thole;
    float ci    = atomI.q;

    di[1]       = atomI.labFrameDipole[0];
    di[2]       = atomI.labFrameDipole[1];
    di[3]       = atomI.labFrameDipole[2];

    qi[1]       = atomI.labFrameQuadrupole[0];
    qi[2]       = atomI.labFrameQuadrupole[1];
    qi[3]       = atomI.labFrameQuadrupole[2];
    qi[4]       = atomI.labFrameQuadrupole[3];
    qi[5]       = atomI.labFrameQuadrupole[4];
    qi[6]       = atomI.labFrameQuadrupole[5];
    qi[7]       = atomI.labFrameQuadrupole[6];
    qi[8]       = atomI.labFrameQuadrupole[7];
    qi[9]       = atomI.labFrameQuadrupole[8];

    float xr    = atomJ.x - atomI.x;
    float yr    = atomJ.y - atomI.y;
    float zr    = atomJ.z - atomI.z;

    // periodic box 

    xr         -= floor(xr*cSim.invPeriodicBoxSizeX+0.5f)*cSim.periodicBoxSizeX;
    yr         -= floor(yr*cSim.invPeriodicBoxSizeY+0.5f)*cSim.periodicBoxSizeY;
    zr         -= floor(zr*cSim.invPeriodicBoxSizeZ+0.5f)*cSim.periodicBoxSizeZ;

    float r2    = xr*xr + yr*yr + zr*zr;
    if( r2 <= cSim.nonbondedCutoffSqr ){

        float r      = sqrt(r2);
        float ck     = atomJ.q;
      
        dk[1]  = atomJ.labFrameDipole[0];
        dk[2]  = atomJ.labFrameDipole[1];
        dk[3]  = atomJ.labFrameDipole[2];
      
        qk[1]  = atomJ.labFrameQuadrupole[0];
        qk[2]  = atomJ.labFrameQuadrupole[1];
        qk[3]  = atomJ.labFrameQuadrupole[2];
        qk[4]  = atomJ.labFrameQuadrupole[3];
        qk[5]  = atomJ.labFrameQuadrupole[4];
        qk[6]  = atomJ.labFrameQuadrupole[5];
        qk[7]  = atomJ.labFrameQuadrupole[6];
        qk[8]  = atomJ.labFrameQuadrupole[7];
        qk[9]  = atomJ.labFrameQuadrupole[8];
      
        // calculate the real space error function terms;

        float ralpha = cSim.alphaEwald*r;

               bn[0] = erfc(ralpha)/r;

        float alsq2  = 2.0f*cSim.alphaEwald*cSim.alphaEwald;
        float alsq2n = 0.0f;
        if( cSim.alphaEwald > 0.0f){
            alsq2n = 1.0f/(cAmoebaSim.sqrtPi*cSim.alphaEwald);
        }
        float exp2a  = exp(-(ralpha*ralpha));

        alsq2n      *= alsq2;
        bn[1]        = (bn[0]+alsq2n*exp2a)/r2;

        alsq2n      *= alsq2;
        bn[2]        = (3.0f*bn[1]+alsq2n*exp2a)/r2;

        alsq2n      *= alsq2;
        bn[3]        = (5.0f*bn[2]+alsq2n*exp2a)/r2;

        alsq2n      *= alsq2;
        bn[4]        = (7.0f*bn[3]+alsq2n*exp2a)/r2;

        alsq2n      *= alsq2;
        bn[5]        = (9.0f*bn[4]+alsq2n*exp2a)/r2;

        // apply Thole polarization damping to scale factors;

        float rr1    = 1.0f/r;
        float rr3    = rr1 / r2;
        float rr5    = 3.0f * rr3 / r2;
        float rr7    = 5.0f * rr5 / r2;
        float rr9    = 7.0f * rr7 / r2;
        float rr11   = 9.0f * rr9 / r2;
        float scale3 = 1.0f;
        float scale5 = 1.0f;
        float scale7 = 1.0f;

            ddsc3[1] = 0.0f;
            ddsc3[2] = 0.0f;
            ddsc3[3] = 0.0f;

            ddsc5[1] = 0.0f;
            ddsc5[2] = 0.0f;
            ddsc5[3] = 0.0f;

            ddsc7[1] = 0.0f;
            ddsc7[2] = 0.0f;
            ddsc7[3] = 0.0f;

        float pdk    = atomJ.damp;
        float ptk    = atomJ.thole;
        float damp   = pdi*pdk;
        if( damp != 0.0f ){
            float pgamma = pti < ptk ? pti : ptk;
            float ratio  = r/damp;
                damp     = -pgamma * ratio*ratio*ratio;
            if( damp > -50.0f ){
                float expdamp  = exp(damp);
                scale3         = 1.0f - expdamp;
                   scale5      = 1.0f - (1.0f-damp)*expdamp;
                   scale7      = 1.0f - (1.0f-damp+0.6f*damp*damp)*expdamp;
                float temp3    = -3.0f * damp * expdamp / r2;
                float temp5    = -damp;
                float temp7    = -0.2f - 0.6f*damp;

                ddsc3[1]       = temp3 * xr;
                ddsc3[2]       = temp3 * yr;
                ddsc3[3]       = temp3 * zr;

                ddsc5[1]       = temp5 * ddsc3[1];
                ddsc5[2]       = temp5 * ddsc3[2];
                ddsc5[3]       = temp5 * ddsc3[3];

                ddsc7[1]       = temp7 * ddsc5[1];
                ddsc7[2]       = temp7 * ddsc5[2];
                ddsc7[3]       = temp7 * ddsc5[3];
            }
        }

        float dsc3 = 1.0f - scale3*scalingFactors[DScaleIndex];
        float dsc5 = 1.0f - scale5*scalingFactors[DScaleIndex];
        float dsc7 = 1.0f - scale7*scalingFactors[DScaleIndex];

        float psc3 = 1.0f - scale3*scalingFactors[PScaleIndex];
        float psc5 = 1.0f - scale5*scalingFactors[PScaleIndex];
        float psc7 = 1.0f - scale7*scalingFactors[PScaleIndex];

        float usc3 = 1.0f - scale3*scalingFactors[UScaleIndex];
        float usc5 = 1.0f - scale5*scalingFactors[UScaleIndex];

        // construct necessary auxiliary vectors

        dixdk[1]       = di[2]*dk[3] - di[3]*dk[2];
        dixdk[2]       = di[3]*dk[1] - di[1]*dk[3];
        dixdk[3]       = di[1]*dk[2] - di[2]*dk[1];

        dixuk[1]       = di[2]*atomJ.inducedDipole[2] - di[3]*atomJ.inducedDipole[1];
        dixuk[2]       = di[3]*atomJ.inducedDipole[0] - di[1]*atomJ.inducedDipole[2];
        dixuk[3]       = di[1]*atomJ.inducedDipole[1] - di[2]*atomJ.inducedDipole[0];
        dkxui[1]       = dk[2]*atomI.inducedDipole[2] - dk[3]*atomI.inducedDipole[1];
        dkxui[2]       = dk[3]*atomI.inducedDipole[0] - dk[1]*atomI.inducedDipole[2];
        dkxui[3]       = dk[1]*atomI.inducedDipole[1] - dk[2]*atomI.inducedDipole[0];
        dixukp[1]      = di[2]*atomJ.inducedDipoleP[2] - di[3]*atomJ.inducedDipoleP[1];
        dixukp[2]      = di[3]*atomJ.inducedDipoleP[0] - di[1]*atomJ.inducedDipoleP[2];
        dixukp[3]      = di[1]*atomJ.inducedDipoleP[1] - di[2]*atomJ.inducedDipoleP[0];
        dkxuip[1]      = dk[2]*atomI.inducedDipoleP[2] - dk[3]*atomI.inducedDipoleP[1];
        dkxuip[2]      = dk[3]*atomI.inducedDipoleP[0] - dk[1]*atomI.inducedDipoleP[2];
        dkxuip[3]      = dk[1]*atomI.inducedDipoleP[1] - dk[2]*atomI.inducedDipoleP[0];
        dixr[1]        = di[2]*zr - di[3]*yr;
        dixr[2]        = di[3]*xr - di[1]*zr;
        dixr[3]        = di[1]*yr - di[2]*xr;
        dkxr[1]        = dk[2]*zr - dk[3]*yr;
        dkxr[2]        = dk[3]*xr - dk[1]*zr;
        dkxr[3]        = dk[1]*yr - dk[2]*xr;
        qir[1]         = qi[1]*xr + qi[4]*yr + qi[7]*zr;
        qir[2]         = qi[2]*xr + qi[5]*yr + qi[8]*zr;
        qir[3]         = qi[3]*xr + qi[6]*yr + qi[9]*zr;
        qkr[1]         = qk[1]*xr + qk[4]*yr + qk[7]*zr;
        qkr[2]         = qk[2]*xr + qk[5]*yr + qk[8]*zr;
        qkr[3]         = qk[3]*xr + qk[6]*yr + qk[9]*zr;
        qiqkr[1]       = qi[1]*qkr[1] + qi[4]*qkr[2] + qi[7]*qkr[3];
        qiqkr[2]       = qi[2]*qkr[1] + qi[5]*qkr[2] + qi[8]*qkr[3];
        qiqkr[3]       = qi[3]*qkr[1] + qi[6]*qkr[2] + qi[9]*qkr[3];
        qkqir[1]       = qk[1]*qir[1] + qk[4]*qir[2] + qk[7]*qir[3];
        qkqir[2]       = qk[2]*qir[1] + qk[5]*qir[2] + qk[8]*qir[3];
        qkqir[3]       = qk[3]*qir[1] + qk[6]*qir[2] + qk[9]*qir[3];
        qixqk[1]       = qi[2]*qk[3] + qi[5]*qk[6] + qi[8]*qk[9]
                       - qi[3]*qk[2] - qi[6]*qk[5] - qi[9]*qk[8];
        qixqk[2]       = qi[3]*qk[1] + qi[6]*qk[4] + qi[9]*qk[7]
                       - qi[1]*qk[3] - qi[4]*qk[6] - qi[7]*qk[9];
        qixqk[3]       = qi[1]*qk[2] + qi[4]*qk[5] + qi[7]*qk[8]
                       - qi[2]*qk[1] - qi[5]*qk[4] - qi[8]*qk[7];
        rxqir[1]       = yr*qir[3] - zr*qir[2];
        rxqir[2]       = zr*qir[1] - xr*qir[3];
        rxqir[3]       = xr*qir[2] - yr*qir[1];
        rxqkr[1]       = yr*qkr[3] - zr*qkr[2];
        rxqkr[2]       = zr*qkr[1] - xr*qkr[3];
        rxqkr[3]       = xr*qkr[2] - yr*qkr[1];
        rxqikr[1]      = yr*qiqkr[3] - zr*qiqkr[2];
        rxqikr[2]      = zr*qiqkr[1] - xr*qiqkr[3];
        rxqikr[3]      = xr*qiqkr[2] - yr*qiqkr[1];
        rxqkir[1]      = yr*qkqir[3] - zr*qkqir[2];
        rxqkir[2]      = zr*qkqir[1] - xr*qkqir[3];
        rxqkir[3]      = xr*qkqir[2] - yr*qkqir[1];
        qkrxqir[1]     = qkr[2]*qir[3] - qkr[3]*qir[2];
        qkrxqir[2]     = qkr[3]*qir[1] - qkr[1]*qir[3];
        qkrxqir[3]     = qkr[1]*qir[2] - qkr[2]*qir[1];
        qidk[1]        = qi[1]*dk[1] + qi[4]*dk[2] + qi[7]*dk[3];
        qidk[2]        = qi[2]*dk[1] + qi[5]*dk[2] + qi[8]*dk[3];
        qidk[3]        = qi[3]*dk[1] + qi[6]*dk[2] + qi[9]*dk[3];
        qkdi[1]        = qk[1]*di[1] + qk[4]*di[2] + qk[7]*di[3];
        qkdi[2]        = qk[2]*di[1] + qk[5]*di[2] + qk[8]*di[3];
        qkdi[3]        = qk[3]*di[1] + qk[6]*di[2] + qk[9]*di[3];
        qiuk[1]        = qi[1]*atomJ.inducedDipole[0] + qi[4]*atomJ.inducedDipole[1]
                       + qi[7]*atomJ.inducedDipole[2];
        qiuk[2]        = qi[2]*atomJ.inducedDipole[0] + qi[5]*atomJ.inducedDipole[1]
                       + qi[8]*atomJ.inducedDipole[2];
        qiuk[3]        = qi[3]*atomJ.inducedDipole[0] + qi[6]*atomJ.inducedDipole[1] 
                       + qi[9]*atomJ.inducedDipole[2];
        qkui[1]        = qk[1]*atomI.inducedDipole[0] + qk[4]*atomI.inducedDipole[1]
                       + qk[7]*atomI.inducedDipole[2];
        qkui[2]        = qk[2]*atomI.inducedDipole[0] + qk[5]*atomI.inducedDipole[1]
                       + qk[8]*atomI.inducedDipole[2];
        qkui[3]        = qk[3]*atomI.inducedDipole[0] + qk[6]*atomI.inducedDipole[1]
                       + qk[9]*atomI.inducedDipole[2];
        qiukp[1]       = qi[1]*atomJ.inducedDipoleP[0] + qi[4]*atomJ.inducedDipoleP[1]
                        + qi[7]*atomJ.inducedDipoleP[2];
        qiukp[2]       = qi[2]*atomJ.inducedDipoleP[0] + qi[5]*atomJ.inducedDipoleP[1]
                        + qi[8]*atomJ.inducedDipoleP[2];
        qiukp[3]       = qi[3]*atomJ.inducedDipoleP[0] + qi[6]*atomJ.inducedDipoleP[1]
                        + qi[9]*atomJ.inducedDipoleP[2];
        qkuip[1]       = qk[1]*atomI.inducedDipoleP[0] + qk[4]*atomI.inducedDipoleP[1]
                        + qk[7]*atomI.inducedDipoleP[2];
        qkuip[2]       = qk[2]*atomI.inducedDipoleP[0] + qk[5]*atomI.inducedDipoleP[1]
                        + qk[8]*atomI.inducedDipoleP[2];
        qkuip[3]       = qk[3]*atomI.inducedDipoleP[0] + qk[6]*atomI.inducedDipoleP[1]
                        + qk[9]*atomI.inducedDipoleP[2];
        dixqkr[1]      = di[2]*qkr[3] - di[3]*qkr[2];
        dixqkr[2]      = di[3]*qkr[1] - di[1]*qkr[3];
        dixqkr[3]      = di[1]*qkr[2] - di[2]*qkr[1];
        dkxqir[1]      = dk[2]*qir[3] - dk[3]*qir[2];
        dkxqir[2]      = dk[3]*qir[1] - dk[1]*qir[3];
        dkxqir[3]      = dk[1]*qir[2] - dk[2]*qir[1];
        uixqkr[1]      = atomI.inducedDipole[1]*qkr[3] - atomI.inducedDipole[2]*qkr[2];
        uixqkr[2]      = atomI.inducedDipole[2]*qkr[1] - atomI.inducedDipole[0]*qkr[3];
        uixqkr[3]      = atomI.inducedDipole[0]*qkr[2] - atomI.inducedDipole[1]*qkr[1];
        ukxqir[1]      = atomJ.inducedDipole[1]*qir[3] - atomJ.inducedDipole[2]*qir[2];
        ukxqir[2]      = atomJ.inducedDipole[2]*qir[1] - atomJ.inducedDipole[0]*qir[3];
        ukxqir[3]      = atomJ.inducedDipole[0]*qir[2] - atomJ.inducedDipole[1]*qir[1];
        uixqkrp[1]     = atomI.inducedDipoleP[1]*qkr[3] - atomI.inducedDipoleP[2]*qkr[2];
        uixqkrp[2]     = atomI.inducedDipoleP[2]*qkr[1] - atomI.inducedDipoleP[0]*qkr[3];
        uixqkrp[3]     = atomI.inducedDipoleP[0]*qkr[2] - atomI.inducedDipoleP[1]*qkr[1];
        ukxqirp[1]     = atomJ.inducedDipoleP[1]*qir[3] - atomJ.inducedDipoleP[2]*qir[2];
        ukxqirp[2]     = atomJ.inducedDipoleP[2]*qir[1] - atomJ.inducedDipoleP[0]*qir[3];
        ukxqirp[3]     = atomJ.inducedDipoleP[0]*qir[2] - atomJ.inducedDipoleP[1]*qir[1];
        rxqidk[1]      = yr*qidk[3] - zr*qidk[2];
        rxqidk[2]      = zr*qidk[1] - xr*qidk[3];
        rxqidk[3]      = xr*qidk[2] - yr*qidk[1];
        rxqkdi[1]      = yr*qkdi[3] - zr*qkdi[2];
        rxqkdi[2]      = zr*qkdi[1] - xr*qkdi[3];
        rxqkdi[3]      = xr*qkdi[2] - yr*qkdi[1];
        rxqiuk[1]      = yr*qiuk[3] - zr*qiuk[2];
        rxqiuk[2]      = zr*qiuk[1] - xr*qiuk[3];
        rxqiuk[3]      = xr*qiuk[2] - yr*qiuk[1];
        rxqkui[1]      = yr*qkui[3] - zr*qkui[2];
        rxqkui[2]      = zr*qkui[1] - xr*qkui[3];
        rxqkui[3]      = xr*qkui[2] - yr*qkui[1];
        rxqiukp[1]     = yr*qiukp[3] - zr*qiukp[2];
        rxqiukp[2]     = zr*qiukp[1] - xr*qiukp[3];
        rxqiukp[3]     = xr*qiukp[2] - yr*qiukp[1];
        rxqkuip[1]     = yr*qkuip[3] - zr*qkuip[2];
        rxqkuip[2]     = zr*qkuip[1] - xr*qkuip[3];
        rxqkuip[3]     = xr*qkuip[2] - yr*qkuip[1];

        // calculate the scalar products for permanent components

        sc[2]          = di[1]*dk[1] + di[2]*dk[2] + di[3]*dk[3];
        sc[3]          = di[1]*xr + di[2]*yr + di[3]*zr;
        sc[4]          = dk[1]*xr + dk[2]*yr + dk[3]*zr;
        sc[5]          = qir[1]*xr + qir[2]*yr + qir[3]*zr;
        sc[6]          = qkr[1]*xr + qkr[2]*yr + qkr[3]*zr;
        sc[7]          = qir[1]*dk[1] + qir[2]*dk[2] + qir[3]*dk[3];
        sc[8]          = qkr[1]*di[1] + qkr[2]*di[2] + qkr[3]*di[3];
        sc[9]          = qir[1]*qkr[1] + qir[2]*qkr[2] + qir[3]*qkr[3];
        sc[10]         = qi[1]*qk[1] + qi[2]*qk[2] + qi[3]*qk[3]
                       + qi[4]*qk[4] + qi[5]*qk[5] + qi[6]*qk[6]
                       + qi[7]*qk[7] + qi[8]*qk[8] + qi[9]*qk[9];

        // calculate the scalar products for induced components

        sci[1]          = atomI.inducedDipole[0]*dk[1] + atomI.inducedDipole[1]*dk[2]
                      + atomI.inducedDipole[2]*dk[3] + di[1]*atomJ.inducedDipole[0]
                      + di[2]*atomJ.inducedDipole[1] + di[3]*atomJ.inducedDipole[2];

        sci[2]          = atomI.inducedDipole[0]*atomJ.inducedDipole[0] + atomI.inducedDipole[1]*atomJ.inducedDipole[1]
                        + atomI.inducedDipole[2]*atomJ.inducedDipole[2];
        sci[3]          = atomI.inducedDipole[0]*xr + atomI.inducedDipole[1]*yr + atomI.inducedDipole[2]*zr;
        sci[4]          = atomJ.inducedDipole[0]*xr + atomJ.inducedDipole[1]*yr + atomJ.inducedDipole[2]*zr;
        sci[7]          = qir[1]*atomJ.inducedDipole[0] + qir[2]*atomJ.inducedDipole[1]
                        + qir[3]*atomJ.inducedDipole[2];
        sci[8]          = qkr[1]*atomI.inducedDipole[0] + qkr[2]*atomI.inducedDipole[1]
                        + qkr[3]*atomI.inducedDipole[2];
        scip[1]         = atomI.inducedDipoleP[0]*dk[1] + atomI.inducedDipoleP[1]*dk[2]
                        + atomI.inducedDipoleP[2]*dk[3] + di[1]*atomJ.inducedDipoleP[0]
                        + di[2]*atomJ.inducedDipoleP[1] + di[3]*atomJ.inducedDipoleP[2];
        scip[2]         = atomI.inducedDipole[0]*atomJ.inducedDipoleP[0]+atomI.inducedDipole[1]*atomJ.inducedDipoleP[1]
                        + atomI.inducedDipole[2]*atomJ.inducedDipoleP[2]+atomI.inducedDipoleP[0]*atomJ.inducedDipole[0]
                        + atomI.inducedDipoleP[1]*atomJ.inducedDipole[1]+atomI.inducedDipoleP[2]*atomJ.inducedDipole[2];
        scip[3]         = atomI.inducedDipoleP[0]*xr + atomI.inducedDipoleP[1]*yr + atomI.inducedDipoleP[2]*zr;
        scip[4]         = atomJ.inducedDipoleP[0]*xr + atomJ.inducedDipoleP[1]*yr + atomJ.inducedDipoleP[2]*zr;
        scip[7]         = qir[1]*atomJ.inducedDipoleP[0] + qir[2]*atomJ.inducedDipoleP[1]
                        + qir[3]*atomJ.inducedDipoleP[2];
        scip[8]         = qkr[1]*atomI.inducedDipoleP[0] + qkr[2]*atomI.inducedDipoleP[1]
                        + qkr[3]*atomI.inducedDipoleP[2];

        // calculate the gl functions for permanent components

        gl[0]           = ci*ck;
        gl[1]           = ck*sc[3] - ci*sc[4];
        gl[2]           = ci*sc[6] + ck*sc[5] - sc[3]*sc[4];
        gl[3]           = sc[3]*sc[6] - sc[4]*sc[5];
        gl[4]           = sc[5]*sc[6];
        gl[5]           = -4.0f * sc[9];
        gl[6]           = sc[2];
        gl[7]           = 2.0f * (sc[7]-sc[8]);
        gl[8]           = 2.0f * sc[10];

        // calculate the gl functions for induced components

        gli[1]          = ck*sci[3] - ci*sci[4];
        gli[2]          = -sc[3]*sci[4] - sci[3]*sc[4];
        gli[3]          = sci[3]*sc[6] - sci[4]*sc[5];
        gli[6]          = sci[1];
        gli[7]          = 2.0f * (sci[7]-sci[8]);
        glip[1]         = ck*scip[3] - ci*scip[4];
        glip[2]         = -sc[3]*scip[4] - scip[3]*sc[4];
        glip[3]         = scip[3]*sc[6] - scip[4]*sc[5];
        glip[6]         = scip[1];
        glip[7]         = 2.0f * (scip[7]-scip[8]);

        // compute the energy contributions for this interaction

        e    = bn[0]*gl[0] + bn[1]*(gl[1]+gl[6])
                 + bn[2]*(gl[2]+gl[7]+gl[8])
                 + bn[3]*(gl[3]+gl[5]) + bn[4]*gl[4];
        ei    = 0.5f * (bn[1]*(gli[1]+gli[6])
                       + bn[2]*(gli[2]+gli[7]) + bn[3]*gli[3]);

        // get the real energy without any screening function

        erl = rr1*gl[0] + rr3*(gl[1]+gl[6])
                   + rr5*(gl[2]+gl[7]+gl[8])
                   + rr7*(gl[3]+gl[5]) + rr9*gl[4];
        erli = 0.5f*(rr3*(gli[1]+gli[6])*psc3
                    + rr5*(gli[2]+gli[7])*psc5
                    + rr7*gli[3]*psc7);
        e = e - (1.0f-scalingFactors[MScaleIndex])*erl;
        ei = ei - erli;

        *energy = -conversionFactor*(e + ei);

        // increment the total intramolecular energy; assumes;
        // intramolecular distances are less than half of cell;
        // length and less than the ewald cutoff;
/*
        if (molcule(ii) .eq. molcule(kk)) {
           eintra = eintra + mscale(kk)*erl*f;
           eintra = eintra + 0.5f*pscale(kk);
&                        * (rr3*(gli[1]+gli[6])*scale3;
&                              + rr5*(gli[2]+gli[7])*scale5;
&                              + rr7*gli[3]*scale7);
        }
*/

        // intermediate variables for permanent force terms

        gf[1] = bn[1]*gl[0] + bn[2]*(gl[1]+gl[6])
                     + bn[3]*(gl[2]+gl[7]+gl[8])
                     + bn[4]*(gl[3]+gl[5]) + bn[5]*gl[4];
        gf[2] = -ck*bn[1] + sc[4]*bn[2] - sc[6]*bn[3];
        gf[3] = ci*bn[1] + sc[3]*bn[2] + sc[5]*bn[3];
        gf[4] = 2.0f * bn[2];
        gf[5] = 2.0f * (-ck*bn[2]+sc[4]*bn[3]-sc[6]*bn[4]);
        gf[6] = 2.0f * (-ci*bn[2]-sc[3]*bn[3]-sc[5]*bn[4]);
        gf[7] = 4.0f * bn[3];
        gfr[1] = rr3*gl[0] + rr5*(gl[1]+gl[6])
                      + rr7*(gl[2]+gl[7]+gl[8])
                      + rr9*(gl[3]+gl[5]) + rr11*gl[4];
        gfr[2] = -ck*rr3 + sc[4]*rr5 - sc[6]*rr7;
        gfr[3] = ci*rr3 + sc[3]*rr5 + sc[5]*rr7;
        gfr[4] = 2.0f * rr5;
        gfr[5] = 2.0f * (-ck*rr5+sc[4]*rr7-sc[6]*rr9);
        gfr[6] = 2.0f * (-ci*rr5-sc[3]*rr7-sc[5]*rr9);
        gfr[7] = 4.0f * rr7;

        // intermediate variables for induced force terms

        gfi[1] = 0.5f*bn[2]*(gli[1]+glip[1]+gli[6]+glip[6])
                      + 0.5f*bn[2]*scip[2]
                      + 0.5f*bn[3]*(gli[2]+glip[2]+gli[7]+glip[7])
                      - 0.5f*bn[3]*(sci[3]*scip[4]+scip[3]*sci[4])
                      + 0.5f*bn[4]*(gli[3]+glip[3]);
        gfi[2] = -ck*bn[1] + sc[4]*bn[2] - sc[6]*bn[3];
        gfi[3] = ci*bn[1] + sc[3]*bn[2] + sc[5]*bn[3];
        gfi[4] = 2.0f * bn[2];
        gfi[5] = bn[3] * (sci[4]+scip[4]);
        gfi[6] = -bn[3] * (sci[3]+scip[3]);
        gfri[1] = 0.5f*rr5*((gli[1]+gli[6])*psc3
                             + (glip[1]+glip[6])*dsc3
                             + scip[2]*usc3)
                  + 0.5f*rr7*((gli[7]+gli[2])*psc5
                             + (glip[7]+glip[2])*dsc5
                      - (sci[3]*scip[4]+scip[3]*sci[4])*usc5)
                  + 0.5f*rr9*(gli[3]*psc7+glip[3]*dsc7);
        gfri[2] = -rr3*ck + rr5*sc[4] - rr7*sc[6];
        gfri[3] = rr3*ci + rr5*sc[3] + rr7*sc[5];
        gfri[4] = 2.0f * rr5;
        gfri[5] = rr7 * (sci[4]*psc7+scip[4]*dsc7);
        gfri[6] = -rr7 * (sci[3]*psc7+scip[3]*dsc7);

        // get the permanent force with screening

        ftm2[1] = gf[1]*xr + gf[2]*di[1] + gf[3]*dk[1]
                       + gf[4]*(qkdi[1]-qidk[1]) + gf[5]*qir[1]
                       + gf[6]*qkr[1] + gf[7]*(qiqkr[1]+qkqir[1]);
        ftm2[2] = gf[1]*yr + gf[2]*di[2] + gf[3]*dk[2]
                       + gf[4]*(qkdi[2]-qidk[2]) + gf[5]*qir[2]
                       + gf[6]*qkr[2] + gf[7]*(qiqkr[2]+qkqir[2]);
        ftm2[3] = gf[1]*zr + gf[2]*di[3] + gf[3]*dk[3]
                       + gf[4]*(qkdi[3]-qidk[3]) + gf[5]*qir[3]
                       + gf[6]*qkr[3] + gf[7]*(qiqkr[3]+qkqir[3]);

        // get the permanent force without screening

        ftm2r[1] = gfr[1]*xr + gfr[2]*di[1] + gfr[3]*dk[1]
                       + gfr[4]*(qkdi[1]-qidk[1]) + gfr[5]*qir[1]
                       + gfr[6]*qkr[1] + gfr[7]*(qiqkr[1]+qkqir[1]);
        ftm2r[2] = gfr[1]*yr + gfr[2]*di[2] + gfr[3]*dk[2]
                       + gfr[4]*(qkdi[2]-qidk[2]) + gfr[5]*qir[2]
                       + gfr[6]*qkr[2] + gfr[7]*(qiqkr[2]+qkqir[2]);
        ftm2r[3] = gfr[1]*zr + gfr[2]*di[3] + gfr[3]*dk[3]
                       + gfr[4]*(qkdi[3]-qidk[3]) + gfr[5]*qir[3]
                       + gfr[6]*qkr[3] + gfr[7]*(qiqkr[3]+qkqir[3]);

        // get the induced force with screening

        ftm2i[1] = gfi[1]*xr + 0.5f*
              (gfi[2]*(atomI.inducedDipole[0]+atomI.inducedDipoleP[0])
             + bn[2]*(sci[4]*atomI.inducedDipoleP[0]+scip[4]*atomI.inducedDipole[0])
             + gfi[3]*(atomJ.inducedDipole[0]+atomJ.inducedDipoleP[0])
             + bn[2]*(sci[3]*atomJ.inducedDipoleP[0]+scip[3]*atomJ.inducedDipole[0])
             + (sci[4]+scip[4])*bn[2]*di[1]
             + (sci[3]+scip[3])*bn[2]*dk[1]
             + gfi[4]*(qkui[1]+qkuip[1]-qiuk[1]-qiukp[1]))
             + gfi[5]*qir[1] + gfi[6]*qkr[1];
        ftm2i[2] = gfi[1]*yr + 0.5f*
              (gfi[2]*(atomI.inducedDipole[1]+atomI.inducedDipoleP[1])
             + bn[2]*(sci[4]*atomI.inducedDipoleP[1]+scip[4]*atomI.inducedDipole[1])
             + gfi[3]*(atomJ.inducedDipole[1]+atomJ.inducedDipoleP[1])
             + bn[2]*(sci[3]*atomJ.inducedDipoleP[1]+scip[3]*atomJ.inducedDipole[1])
             + (sci[4]+scip[4])*bn[2]*di[2]
             + (sci[3]+scip[3])*bn[2]*dk[2]
             + gfi[4]*(qkui[2]+qkuip[2]-qiuk[2]-qiukp[2]))
             + gfi[5]*qir[2] + gfi[6]*qkr[2];
        ftm2i[3] = gfi[1]*zr + 0.5f*
              (gfi[2]*(atomI.inducedDipole[2]+atomI.inducedDipoleP[2])
             + bn[2]*(sci[4]*atomI.inducedDipoleP[2]+scip[4]*atomI.inducedDipole[2])
             + gfi[3]*(atomJ.inducedDipole[2]+atomJ.inducedDipoleP[2])
             + bn[2]*(sci[3]*atomJ.inducedDipoleP[2]+scip[3]*atomJ.inducedDipole[2])
             + (sci[4]+scip[4])*bn[2]*di[3]
             + (sci[3]+scip[3])*bn[2]*dk[3]
             + gfi[4]*(qkui[3]+qkuip[3]-qiuk[3]-qiukp[3]))
             + gfi[5]*qir[3] + gfi[6]*qkr[3];

        // get the induced force without screening

        ftm2ri[1] = gfri[1]*xr + 0.5f*
            (- rr3*ck*(atomI.inducedDipole[0]*psc3+atomI.inducedDipoleP[0]*dsc3)
             + rr5*sc[4]*(atomI.inducedDipole[0]*psc5+atomI.inducedDipoleP[0]*dsc5)
             - rr7*sc[6]*(atomI.inducedDipole[0]*psc7+atomI.inducedDipoleP[0]*dsc7))
             + (rr3*ci*(atomJ.inducedDipole[0]*psc3+atomJ.inducedDipoleP[0]*dsc3)
             + rr5*sc[3]*(atomJ.inducedDipole[0]*psc5+atomJ.inducedDipoleP[0]*dsc5)
             + rr7*sc[5]*(atomJ.inducedDipole[0]*psc7+atomJ.inducedDipoleP[0]*dsc7))*0.5f
             + rr5*usc5*(sci[4]*atomI.inducedDipoleP[0]+scip[4]*atomI.inducedDipole[0]
             + sci[3]*atomJ.inducedDipoleP[0]+scip[3]*atomJ.inducedDipole[0])*0.5f
             + 0.5f*(sci[4]*psc5+scip[4]*dsc5)*rr5*di[1]
             + 0.5f*(sci[3]*psc5+scip[3]*dsc5)*rr5*dk[1]
             + 0.5f*gfri[4]*((qkui[1]-qiuk[1])*psc5
             + (qkuip[1]-qiukp[1])*dsc5)
             + gfri[5]*qir[1] + gfri[6]*qkr[1];
        ftm2ri[2] = gfri[1]*yr + 0.5f*
            (- rr3*ck*(atomI.inducedDipole[1]*psc3+atomI.inducedDipoleP[1]*dsc3)
             + rr5*sc[4]*(atomI.inducedDipole[1]*psc5+atomI.inducedDipoleP[1]*dsc5)
             - rr7*sc[6]*(atomI.inducedDipole[1]*psc7+atomI.inducedDipoleP[1]*dsc7))
             + (rr3*ci*(atomJ.inducedDipole[1]*psc3+atomJ.inducedDipoleP[1]*dsc3)
             + rr5*sc[3]*(atomJ.inducedDipole[1]*psc5+atomJ.inducedDipoleP[1]*dsc5)
             + rr7*sc[5]*(atomJ.inducedDipole[1]*psc7+atomJ.inducedDipoleP[1]*dsc7))*0.5f
             + rr5*usc5*(sci[4]*atomI.inducedDipoleP[1]+scip[4]*atomI.inducedDipole[1]
             + sci[3]*atomJ.inducedDipoleP[1]+scip[3]*atomJ.inducedDipole[1])*0.5f
             + 0.5f*(sci[4]*psc5+scip[4]*dsc5)*rr5*di[2]
             + 0.5f*(sci[3]*psc5+scip[3]*dsc5)*rr5*dk[2]
             + 0.5f*gfri[4]*((qkui[2]-qiuk[2])*psc5
             + (qkuip[2]-qiukp[2])*dsc5)
             + gfri[5]*qir[2] + gfri[6]*qkr[2];
        ftm2ri[3] = gfri[1]*zr + 0.5f*
            (- rr3*ck*(atomI.inducedDipole[2]*psc3+atomI.inducedDipoleP[2]*dsc3)
             + rr5*sc[4]*(atomI.inducedDipole[2]*psc5+atomI.inducedDipoleP[2]*dsc5)
             - rr7*sc[6]*(atomI.inducedDipole[2]*psc7+atomI.inducedDipoleP[2]*dsc7))
             + (rr3*ci*(atomJ.inducedDipole[2]*psc3+atomJ.inducedDipoleP[2]*dsc3)
             + rr5*sc[3]*(atomJ.inducedDipole[2]*psc5+atomJ.inducedDipoleP[2]*dsc5)
             + rr7*sc[5]*(atomJ.inducedDipole[2]*psc7+atomJ.inducedDipoleP[2]*dsc7))*0.5f
             + rr5*usc5*(sci[4]*atomI.inducedDipoleP[2]+scip[4]*atomI.inducedDipole[2]
             + sci[3]*atomJ.inducedDipoleP[2]+scip[3]*atomJ.inducedDipole[2])*0.5f
             + 0.5f*(sci[4]*psc5+scip[4]*dsc5)*rr5*di[3]
             + 0.5f*(sci[3]*psc5+scip[3]*dsc5)*rr5*dk[3]
             + 0.5f*gfri[4]*((qkui[3]-qiuk[3])*psc5
             + (qkuip[3]-qiukp[3])*dsc5)
             + gfri[5]*qir[3] + gfri[6]*qkr[3];

        // account for partially excluded induced interactions

        float temp3 = 0.5f * rr3 * ((gli[1]+gli[6])*scalingFactors[PScaleIndex]
                                   +(glip[1]+glip[6])*scalingFactors[DScaleIndex]);
        float temp5 = 0.5f * rr5 * ((gli[2]+gli[7])*scalingFactors[PScaleIndex]
                                   +(glip[2]+glip[7])*scalingFactors[DScaleIndex]);
        float temp7 = 0.5f * rr7 * (gli[3]*scalingFactors[PScaleIndex]
                                   +glip[3]*scalingFactors[DScaleIndex]);
        fridmp[1] = temp3*ddsc3[1] + temp5*ddsc5[1] + temp7*ddsc7[1];
        fridmp[2] = temp3*ddsc3[2] + temp5*ddsc5[2] + temp7*ddsc7[2];
        fridmp[3] = temp3*ddsc3[3] + temp5*ddsc5[3] + temp7*ddsc7[3];

        // find some scaling terms for induced-induced force

        temp3 = 0.5f * rr3 * scalingFactors[UScaleIndex] * scip[2];
        temp5 = -0.5f * rr5 * scalingFactors[UScaleIndex] * (sci[3]*scip[4]+scip[3]*sci[4]);
        findmp[1] = temp3*ddsc3[1] + temp5*ddsc5[1];
        findmp[2] = temp3*ddsc3[2] + temp5*ddsc5[2];
        findmp[3] = temp3*ddsc3[3] + temp5*ddsc5[3];

        // modify the forces for partially excluded interactions

        ftm2i[1] = ftm2i[1] - fridmp[1] - findmp[1];
        ftm2i[2] = ftm2i[2] - fridmp[2] - findmp[2];
        ftm2i[3] = ftm2i[3] - fridmp[3] - findmp[3];

        // correction to convert mutual to direct polarization force

/*
        if (poltyp .eq. 'DIRECT') {
           gfd = 0.5f * (bn[2]*scip[2];
&                     - bn[3]*(scip[3]*sci[4]+sci[3]*scip[4]));
           gfdr = 0.5f * (rr5*scip[2]*usc3;
&                     - rr7*(scip[3]*sci[4];
&                           +sci[3]*scip[4])*usc5);
           ftm2i[1] = ftm2i[1] - gfd*xr - 0.5f*bn[2]*;
&                          (sci[4]*atomI.inducedDipoleP[0]+scip[4]*atomI.inducedDipole[0];
&                          +sci[3]*atomJ.inducedDipoleP[0]+scip[3]*atomJ.inducedDipole[0]);
           ftm2i[2] = ftm2i[2] - gfd*yr - 0.5f*bn[2]*;
&                          (sci[4]*atomI.inducedDipoleP[1]+scip[4]*atomI.inducedDipole[1];
&                          +sci[3]*atomJ.inducedDipoleP[1]+scip[3]*atomJ.inducedDipole[1]);
           ftm2i[3] = ftm2i[3] - gfd*zr - 0.5f*bn[2]*;
&                          (sci[4]*atomI.inducedDipoleP[2]+scip[4]*atomI.inducedDipole[2];
&                          +sci[3]*atomJ.inducedDipoleP[2]+scip[3]*atomJ.inducedDipole[2]);
           fdir[1] = gfdr*xr + 0.5f*usc5*rr5*;
&                         (sci[4]*atomI.inducedDipoleP[0]+scip[4]*atomI.inducedDipole[0];
&                        + sci[3]*atomJ.inducedDipoleP[0]+scip[3]*atomJ.inducedDipole[0]);
           fdir[2] = gfdr*yr + 0.5f*usc5*rr5*;
&                         (sci[4]*atomI.inducedDipoleP[1]+scip[4]*atomI.inducedDipole[1];
&                        + sci[3]*atomJ.inducedDipoleP[1]+scip[3]*atomJ.inducedDipole[1]);
           fdir[3] = gfdr*zr + 0.5f*usc5*rr5*;
&                         (sci[4]*atomI.inducedDipoleP[2]+scip[4]*atomI.inducedDipole[2];
&                        + sci[3]*atomJ.inducedDipoleP[2]+scip[3]*atomJ.inducedDipole[2]);
           ftm2i[1] = ftm2i[1] + fdir[1] + findmp[1];
           ftm2i[2] = ftm2i[2] + fdir[2] + findmp[2];
           ftm2i[3] = ftm2i[3] + fdir[3] + findmp[3];
        }
*/

        // intermediate variables for induced torque terms

        gti[2] = 0.5f * bn[2] * (sci[4]+scip[4]);
        gti[3] = 0.5f * bn[2] * (sci[3]+scip[3]);
        gti[4] = gfi[4];
        gti[5] = gfi[5];
        gti[6] = gfi[6];
        gtri[2] = 0.5f * rr5 * (sci[4]*psc5+scip[4]*dsc5);
        gtri[3] = 0.5f * rr5 * (sci[3]*psc5+scip[3]*dsc5);
        gtri[4] = gfri[4];
        gtri[5] = gfri[5];
        gtri[6] = gfri[6];

        // get the permanent torque with screening

        ttm2[1] = -bn[1]*dixdk[1] + gf[2]*dixr[1]
            + gf[4]*(dixqkr[1]+dkxqir[1]+rxqidk[1]-2.0f*qixqk[1])
            - gf[5]*rxqir[1] - gf[7]*(rxqikr[1]+qkrxqir[1]);
        ttm2[2] = -bn[1]*dixdk[2] + gf[2]*dixr[2]
            + gf[4]*(dixqkr[2]+dkxqir[2]+rxqidk[2]-2.0f*qixqk[2])
            - gf[5]*rxqir[2] - gf[7]*(rxqikr[2]+qkrxqir[2]);
        ttm2[3] = -bn[1]*dixdk[3] + gf[2]*dixr[3]
            + gf[4]*(dixqkr[3]+dkxqir[3]+rxqidk[3]-2.0f*qixqk[3])
            - gf[5]*rxqir[3] - gf[7]*(rxqikr[3]+qkrxqir[3]);
        ttm3[1] = bn[1]*dixdk[1] + gf[3]*dkxr[1]
            - gf[4]*(dixqkr[1]+dkxqir[1]+rxqkdi[1]-2.0f*qixqk[1])
            - gf[6]*rxqkr[1] - gf[7]*(rxqkir[1]-qkrxqir[1]);
        ttm3[2] = bn[1]*dixdk[2] + gf[3]*dkxr[2]
            - gf[4]*(dixqkr[2]+dkxqir[2]+rxqkdi[2]-2.0f*qixqk[2])
            - gf[6]*rxqkr[2] - gf[7]*(rxqkir[2]-qkrxqir[2]);
        ttm3[3] = bn[1]*dixdk[3] + gf[3]*dkxr[3]
            - gf[4]*(dixqkr[3]+dkxqir[3]+rxqkdi[3]-2.0f*qixqk[3])
            - gf[6]*rxqkr[3] - gf[7]*(rxqkir[3]-qkrxqir[3]);

        // get the permanent torque without screening

        ttm2r[1] = -rr3*dixdk[1] + gfr[2]*dixr[1]-gfr[5]*rxqir[1]
            + gfr[4]*(dixqkr[1]+dkxqir[1]+rxqidk[1]-2.0f*qixqk[1])
            - gfr[7]*(rxqikr[1]+qkrxqir[1]);
        ttm2r[2] = -rr3*dixdk[2] + gfr[2]*dixr[2]-gfr[5]*rxqir[2]
            + gfr[4]*(dixqkr[2]+dkxqir[2]+rxqidk[2]-2.0f*qixqk[2])
            - gfr[7]*(rxqikr[2]+qkrxqir[2]);
        ttm2r[3] = -rr3*dixdk[3] + gfr[2]*dixr[3]-gfr[5]*rxqir[3]
            + gfr[4]*(dixqkr[3]+dkxqir[3]+rxqidk[3]-2.0f*qixqk[3])
            - gfr[7]*(rxqikr[3]+qkrxqir[3]);
        ttm3r[1] = rr3*dixdk[1] + gfr[3]*dkxr[1] -gfr[6]*rxqkr[1]
            - gfr[4]*(dixqkr[1]+dkxqir[1]+rxqkdi[1]-2.0f*qixqk[1])
            - gfr[7]*(rxqkir[1]-qkrxqir[1]);
        ttm3r[2] = rr3*dixdk[2] + gfr[3]*dkxr[2] -gfr[6]*rxqkr[2]
            - gfr[4]*(dixqkr[2]+dkxqir[2]+rxqkdi[2]-2.0f*qixqk[2])
            - gfr[7]*(rxqkir[2]-qkrxqir[2]);
        ttm3r[3] = rr3*dixdk[3] + gfr[3]*dkxr[3] -gfr[6]*rxqkr[3]
            - gfr[4]*(dixqkr[3]+dkxqir[3]+rxqkdi[3]-2.0f*qixqk[3])
            - gfr[7]*(rxqkir[3]-qkrxqir[3]);

        // get the induced torque with screening

        ttm2i[1] = -bn[1]*(dixuk[1]+dixukp[1])*0.5f
            + gti[2]*dixr[1] + gti[4]*(ukxqir[1]+rxqiuk[1]
            + ukxqirp[1]+rxqiukp[1])*0.5f - gti[5]*rxqir[1];
        ttm2i[2] = -bn[1]*(dixuk[2]+dixukp[2])*0.5f
            + gti[2]*dixr[2] + gti[4]*(ukxqir[2]+rxqiuk[2]
            + ukxqirp[2]+rxqiukp[2])*0.5f - gti[5]*rxqir[2];
        ttm2i[3] = -bn[1]*(dixuk[3]+dixukp[3])*0.5f
            + gti[2]*dixr[3] + gti[4]*(ukxqir[3]+rxqiuk[3]
            + ukxqirp[3]+rxqiukp[3])*0.5f - gti[5]*rxqir[3];
        ttm3i[1] = -bn[1]*(dkxui[1]+dkxuip[1])*0.5f
            + gti[3]*dkxr[1] - gti[4]*(uixqkr[1]+rxqkui[1]
            + uixqkrp[1]+rxqkuip[1])*0.5f - gti[6]*rxqkr[1];
        ttm3i[2] = -bn[1]*(dkxui[2]+dkxuip[2])*0.5f
            + gti[3]*dkxr[2] - gti[4]*(uixqkr[2]+rxqkui[2]
            + uixqkrp[2]+rxqkuip[2])*0.5f - gti[6]*rxqkr[2];
        ttm3i[3] = -bn[1]*(dkxui[3]+dkxuip[3])*0.5f
            + gti[3]*dkxr[3] - gti[4]*(uixqkr[3]+rxqkui[3]
            + uixqkrp[3]+rxqkuip[3])*0.5f - gti[6]*rxqkr[3];

        // get the induced torque without screening

        ttm2ri[1] = -rr3*(dixuk[1]*psc3+dixukp[1]*dsc3)*0.5f
            + gtri[2]*dixr[1] + gtri[4]*((ukxqir[1]+rxqiuk[1])*psc5
            +(ukxqirp[1]+rxqiukp[1])*dsc5)*0.5f - gtri[5]*rxqir[1];
        ttm2ri[2] = -rr3*(dixuk[2]*psc3+dixukp[2]*dsc3)*0.5f
            + gtri[2]*dixr[2] + gtri[4]*((ukxqir[2]+rxqiuk[2])*psc5
            +(ukxqirp[2]+rxqiukp[2])*dsc5)*0.5f - gtri[5]*rxqir[2];
        ttm2ri[3] = -rr3*(dixuk[3]*psc3+dixukp[3]*dsc3)*0.5f
            + gtri[2]*dixr[3] + gtri[4]*((ukxqir[3]+rxqiuk[3])*psc5
            +(ukxqirp[3]+rxqiukp[3])*dsc5)*0.5f - gtri[5]*rxqir[3];
        ttm3ri[1] = -rr3*(dkxui[1]*psc3+dkxuip[1]*dsc3)*0.5f
            + gtri[3]*dkxr[1] - gtri[4]*((uixqkr[1]+rxqkui[1])*psc5
            +(uixqkrp[1]+rxqkuip[1])*dsc5)*0.5f - gtri[6]*rxqkr[1];
        ttm3ri[2] = -rr3*(dkxui[2]*psc3+dkxuip[2]*dsc3)*0.5f
            + gtri[3]*dkxr[2] - gtri[4]*((uixqkr[2]+rxqkui[2])*psc5
            +(uixqkrp[2]+rxqkuip[2])*dsc5)*0.5f - gtri[6]*rxqkr[2];
        ttm3ri[3] = -rr3*(dkxui[3]*psc3+dkxuip[3]*dsc3)*0.5f
            + gtri[3]*dkxr[3] - gtri[4]*((uixqkr[3]+rxqkui[3])*psc5
            +(uixqkrp[3]+rxqkuip[3])*dsc5)*0.5f - gtri[6]*rxqkr[3];

        // handle the case where scaling is used

        for( int j = 1; j <= 3; j++ ){
           ftm2[j]  = (ftm2[j]-(1.0f-scalingFactors[MScaleIndex])*ftm2r[j]);
           ftm2i[j] = (ftm2i[j]-ftm2ri[j]);
           ttm2[j]  = (ttm2[j]-(1.0f-scalingFactors[MScaleIndex])*ttm2r[j]);
           ttm2i[j] = (ttm2i[j]-ttm2ri[j]);
           ttm3[j]  = (ttm3[j]-(1.0f-scalingFactors[MScaleIndex])*ttm3r[j]);
           ttm3i[j] = (ttm3i[j]-ttm3ri[j]);
        }

        // increment gradient due to force and torque on first site;

        outputForce[0]           = conversionFactor*(ftm2[1] + ftm2i[1]);
        outputForce[1]           = conversionFactor*(ftm2[2] + ftm2i[2]);
        outputForce[2]           = conversionFactor*(ftm2[3] + ftm2i[3]);
        
        conversionFactor        *= -1.0;
        outputTorque[0][0]       =  conversionFactor*(ttm2[1] + ttm2i[1]);
        outputTorque[0][1]       =  conversionFactor*(ttm2[2] + ttm2i[2]);
        outputTorque[0][2]       =  conversionFactor*(ttm2[3] + ttm2i[3]);
    
        outputTorque[1][0]       =  conversionFactor*(ttm3[1] + ttm3i[1]);
        outputTorque[1][1]       =  conversionFactor*(ttm3[2] + ttm3i[2]);
        outputTorque[1][2]       =  conversionFactor*(ttm3[3] + ttm3i[3]);

#ifdef AMOEBA_DEBUG
    int debugIndex               = 0;
    float idTracker              = 1.0f;
/*
    debugArray[debugIndex].x = atomI.labFrameDipole[0];
    debugArray[debugIndex].y = atomI.labFrameDipole[1];
    debugArray[debugIndex].z = atomI.labFrameDipole[2];
    debugArray[debugIndex].w = r2;

    debugIndex++;
    idTracker               += 1.0;
    debugArray[debugIndex].x = atomJ.labFrameDipole[0];
    debugArray[debugIndex].y = atomJ.labFrameDipole[1];
    debugArray[debugIndex].z = atomJ.labFrameDipole[2];
    debugArray[debugIndex].w = cSim.alphaEwald;

    debugIndex++;
    idTracker               += 1.0;
    debugArray[debugIndex].x = atomI.inducedDipole[0];
    debugArray[debugIndex].y = atomI.inducedDipole[1];
    debugArray[debugIndex].z = atomI.inducedDipole[2];
    debugArray[debugIndex].w = idTracker;

    debugIndex++;
    idTracker               += 1.0;
    debugArray[debugIndex].x = atomJ.inducedDipole[0];
    debugArray[debugIndex].y = atomJ.inducedDipole[1];
    debugArray[debugIndex].z = atomJ.inducedDipole[2];
    debugArray[debugIndex].w = idTracker;

    debugIndex++;
    idTracker               += 1.0;
    debugArray[debugIndex].x = atomI.inducedDipoleP[0];
    debugArray[debugIndex].y = atomI.inducedDipoleP[1];
    debugArray[debugIndex].z = atomI.inducedDipoleP[2];
    debugArray[debugIndex].w = idTracker;

    debugIndex++;
    idTracker               += 1.0;
    debugArray[debugIndex].x = atomJ.inducedDipoleP[0];
    debugArray[debugIndex].y = atomJ.inducedDipoleP[1];
    debugArray[debugIndex].z = atomJ.inducedDipoleP[2];
    debugArray[debugIndex].w = idTracker;

    debugIndex++;
    idTracker               += 1.0;
    debugArray[debugIndex].x = conversionFactor*ftm2[1];
    debugArray[debugIndex].y = conversionFactor*ftm2[2];
    debugArray[debugIndex].z = conversionFactor*ftm2[3];
    debugArray[debugIndex].w = idTracker;
    debugIndex++;
*/
    idTracker               += 1.0;
    debugArray[debugIndex].x = e;
    debugArray[debugIndex].y = ei;
    debugArray[debugIndex].z = erl;
    debugArray[debugIndex].w = erli;

    debugIndex++;
    idTracker               += 1.0;
    debugArray[debugIndex].x = r2;
    debugArray[debugIndex].y = cSim.alphaEwald;
    debugArray[debugIndex].z = conversionFactor*fridmp[3];
    debugArray[debugIndex].w = 115.0;

    debugIndex++;
    idTracker               += 1.0;
    debugArray[debugIndex].x = conversionFactor*findmp[1];
    debugArray[debugIndex].y = conversionFactor*findmp[2];
    debugArray[debugIndex].z = conversionFactor*findmp[3];
    debugArray[debugIndex].w = cSim.alphaEwald + 1.0f;

    debugIndex++;
    idTracker               += 1.0;
    debugArray[debugIndex].x = conversionFactor*ttm2[1];
    debugArray[debugIndex].y = conversionFactor*ttm2[2];
    debugArray[debugIndex].z = conversionFactor*ttm2[3];
    debugArray[debugIndex].w = idTracker;

    debugIndex++;
    idTracker               += 1.0;
    debugArray[debugIndex].x = conversionFactor*ttm2i[1];
    debugArray[debugIndex].y = conversionFactor*ttm2i[2];
    debugArray[debugIndex].z = conversionFactor*ttm2i[3];
    debugArray[debugIndex].w = idTracker;
#endif

    } else {

        outputForce[0]           = 0.0f;
        outputForce[1]           = 0.0f;
        outputForce[2]           = 0.0f;
        
        outputTorque[0][0]       = 0.0f;
        outputTorque[0][1]       = 0.0f;
        outputTorque[0][2]       = 0.0f;
    
        outputTorque[1][0]       = 0.0f;
        outputTorque[1][1]       = 0.0f;
        outputTorque[1][2]       = 0.0f;

        *energy                  = 0.0f;

#ifdef AMOEBA_DEBUG
for( int ii = 0; ii < 5; ii++ ){
    debugArray[ii].x = 0.0f;
    debugArray[ii].y = 0.0f;
    debugArray[ii].z = 0.0f;
    debugArray[ii].w = (float) (11*ii);
}
#endif
   
    }

    return;

}

__device__ void loadPmeDirectElectrostaticShared( struct PmeDirectElectrostaticParticle* sA, unsigned int atomI )
{
    // coordinates & charge
    sA->x                        = cSim.pPosq[atomI].x;
    sA->y                        = cSim.pPosq[atomI].y;
    sA->z                        = cSim.pPosq[atomI].z;
    sA->q                        = cSim.pPosq[atomI].w;

    // lab dipole

    sA->labFrameDipole[0]        = cAmoebaSim.pLabFrameDipole[atomI*3];
    sA->labFrameDipole[1]        = cAmoebaSim.pLabFrameDipole[atomI*3+1];
    sA->labFrameDipole[2]        = cAmoebaSim.pLabFrameDipole[atomI*3+2];


    // lab quadrupole

    sA->labFrameQuadrupole[0]    = cAmoebaSim.pLabFrameQuadrupole[atomI*9];
    sA->labFrameQuadrupole[1]    = cAmoebaSim.pLabFrameQuadrupole[atomI*9+1];
    sA->labFrameQuadrupole[2]    = cAmoebaSim.pLabFrameQuadrupole[atomI*9+2];
    sA->labFrameQuadrupole[3]    = cAmoebaSim.pLabFrameQuadrupole[atomI*9+3];
    sA->labFrameQuadrupole[4]    = cAmoebaSim.pLabFrameQuadrupole[atomI*9+4];
    sA->labFrameQuadrupole[5]    = cAmoebaSim.pLabFrameQuadrupole[atomI*9+5];
    sA->labFrameQuadrupole[6]    = cAmoebaSim.pLabFrameQuadrupole[atomI*9+6];
    sA->labFrameQuadrupole[7]    = cAmoebaSim.pLabFrameQuadrupole[atomI*9+7];
    sA->labFrameQuadrupole[8]    = cAmoebaSim.pLabFrameQuadrupole[atomI*9+8];

    // induced dipole

    sA->inducedDipole[0]         =  cAmoebaSim.pInducedDipole[atomI*3];
    sA->inducedDipole[1]         =  cAmoebaSim.pInducedDipole[atomI*3+1];
    sA->inducedDipole[2]         =  cAmoebaSim.pInducedDipole[atomI*3+2];

    // induced dipole polar

    sA->inducedDipoleP[0]        =  cAmoebaSim.pInducedDipolePolar[atomI*3];
    sA->inducedDipoleP[1]        =  cAmoebaSim.pInducedDipolePolar[atomI*3+1];
    sA->inducedDipoleP[2]        =  cAmoebaSim.pInducedDipolePolar[atomI*3+2];

    sA->damp                     = cAmoebaSim.pDampingFactorAndThole[atomI].x;
    sA->thole                    = cAmoebaSim.pDampingFactorAndThole[atomI].y;

}

// Include versions of the kernels for N^2 calculations.

#undef USE_OUTPUT_BUFFER_PER_WARP
#define METHOD_NAME(a, b) a##N2##b
#include "kCalculateAmoebaCudaPmeDirectElectrostatic.h"
#define USE_OUTPUT_BUFFER_PER_WARP
#undef METHOD_NAME
#define METHOD_NAME(a, b) a##N2ByWarp##b
#include "kCalculateAmoebaCudaPmeDirectElectrostatic.h"

// reduce psWorkArray_3_1 -> force
// reduce psWorkArray_3_2 -> torque

static void kReduceForceTorque(amoebaGpuContext amoebaGpu )
{
    kReduceFields_kernel<<<amoebaGpu->nonbondBlocks, amoebaGpu->fieldReduceThreadsPerBlock>>>(
                             amoebaGpu->paddedNumberOfAtoms*3, amoebaGpu->outputBuffers,
                             amoebaGpu->psWorkArray_3_1->_pDevStream[0], amoebaGpu->psForce->_pDevStream[0] );
    LAUNCHERROR("kReducePmeDirectElectrostaticForce");
    kReduceFields_kernel<<<amoebaGpu->nonbondBlocks, amoebaGpu->fieldReduceThreadsPerBlock>>>(
                             amoebaGpu->paddedNumberOfAtoms*3, amoebaGpu->outputBuffers,
                             amoebaGpu->psWorkArray_3_2->_pDevStream[0], amoebaGpu->psTorque->_pDevStream[0] );
    LAUNCHERROR("kReducePmeDirectElectrostaticTorque");
}

//#define GET_INDUCED_DIPOLE_FROM_FILE
#ifdef GET_INDUCED_DIPOLE_FROM_FILE
#include <stdlib.h>
#endif

/**---------------------------------------------------------------------------------------

   Compute Amoeba dirrect space portion of electrostatic force & torque

   @param amoebaGpu        amoebaGpu context

   --------------------------------------------------------------------------------------- */

void cudaComputeAmoebaPmeDirectElectrostatic( amoebaGpuContext amoebaGpu )
{
  
   // ---------------------------------------------------------------------------------------

    static unsigned int threadsPerBlock = 0;

#ifdef AMOEBA_DEBUG
    static const char* methodName = "cudaComputeAmoebaPmeDirectElectrostatic";
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
    unsigned int targetAtom                   = 10;
#endif

#ifdef GET_INDUCED_DIPOLE_FROM_FILE
    std::string fileName = "waterInducedDipole.txt";
    StringVectorVector fileContents;
    readFile( fileName, fileContents );
    unsigned int offset  = 0;
    (void) fprintf( amoebaGpu->log, "Read file: %s %u\n", fileName.c_str(), fileContents.size() ); fflush(  amoebaGpu->log );
    for( unsigned int ii = 1; ii < fileContents.size()-1; ii++ ){

        StringVector lineTokens     = fileContents[ii];
        unsigned int lineTokenIndex = 1;

        // (void) fprintf( amoebaGpu->log, "   %u %s %s\n", ii, lineTokens[0].c_str(), lineTokens[lineTokenIndex].c_str() ); fflush(  amoebaGpu->log );
        amoebaGpu->psInducedDipole->_pSysStream[0][offset++]       = static_cast<float>(atof(lineTokens[lineTokenIndex++].c_str()));
        amoebaGpu->psInducedDipole->_pSysStream[0][offset++]       = static_cast<float>(atof(lineTokens[lineTokenIndex++].c_str()));
        amoebaGpu->psInducedDipole->_pSysStream[0][offset++]       = static_cast<float>(atof(lineTokens[lineTokenIndex++].c_str())); 
        offset                                              -= 3;        
        amoebaGpu->psInducedDipolePolar->_pSysStream[0][offset++]  = static_cast<float>(atof(lineTokens[lineTokenIndex++].c_str()));
        amoebaGpu->psInducedDipolePolar->_pSysStream[0][offset++]  = static_cast<float>(atof(lineTokens[lineTokenIndex++].c_str()));
        amoebaGpu->psInducedDipolePolar->_pSysStream[0][offset++]  = static_cast<float>(atof(lineTokens[lineTokenIndex++].c_str()));
    }
    float conversion = 0.1f;
    for( int ii = 0; ii < 3*gpu->natoms; ii++ ){
        amoebaGpu->psInducedDipole->_pSysStream[0][ii]       *= conversion;
        amoebaGpu->psInducedDipolePolar->_pSysStream[0][ii]  *= conversion;
    }
    amoebaGpu->gpuContext->sim.alphaEwald = 5.4459052e+00f;
    SetCalculateAmoebaPmeDirectElectrostaticSim(amoebaGpu);
    amoebaGpu->psInducedDipole->Upload();
    amoebaGpu->psInducedDipolePolar->Upload();
#endif

    // on first pass, set threads/block

    if( threadsPerBlock == 0 ){
        unsigned int maxThreads;
        if (gpu->sm_version >= SM_20)
            maxThreads = 384;
        else if (gpu->sm_version >= SM_12)
            maxThreads = 128;
        else
            maxThreads = 64;
        threadsPerBlock = std::min(getThreadsPerBlock(amoebaGpu, sizeof(PmeDirectElectrostaticParticle)+sizeof(float3)), maxThreads);
    }

    kClearFields_3( amoebaGpu, 2 );

#ifdef AMOEBA_DEBUG
    (void) fprintf( amoebaGpu->log, "kCalculateAmoebaPmeDirectElectrostaticN2Forces:  threadsPerBlock=%u getThreadsPerBlock=%d sizeof=%u\n", 
                    threadsPerBlock, getThreadsPerBlock(amoebaGpu, sizeof(PmeDirectElectrostaticParticle)+sizeof(float3)),
                    (sizeof(PmeDirectElectrostaticParticle)+sizeof(float3)) );

      (void) fprintf( amoebaGpu->log, "kCalculateAmoebaPmeDirectElectrostaticN2Forces no warp:  numBlocks=%u numThreads=%u bufferPerWarp=%u atm=%u shrd=%u Obuf=%u ixnCt=%u workUnits=%u gpu->nonbond_threads_per_block=%u\n",
                      amoebaGpu->nonbondBlocks, threadsPerBlock, amoebaGpu->bOutputBufferPerWarp,
                      sizeof(PmeDirectElectrostaticParticle)+sizeof(float3), (sizeof(PmeDirectElectrostaticParticle)+sizeof(float3))*threadsPerBlock, amoebaGpu->energyOutputBuffers, (*gpu->psInteractionCount)[0], gpu->sim.workUnits,
                      gpu->sim.nonbond_threads_per_block );
      (void) fflush( amoebaGpu->log );
#endif

    if (gpu->bOutputBufferPerWarp){

      kCalculateAmoebaPmeDirectElectrostaticN2ByWarpForces_kernel<<<amoebaGpu->nonbondBlocks, threadsPerBlock, (sizeof(PmeDirectElectrostaticParticle)+sizeof(float3))*threadsPerBlock>>>(
                                                                         gpu->sim.pInteractingWorkUnit,
                                                                         amoebaGpu->psWorkArray_3_1->_pDevStream[0],
#ifdef AMOEBA_DEBUG
                                                                         amoebaGpu->psWorkArray_3_2->_pDevStream[0],
                                                                         debugArray->_pDevStream[0], targetAtom );
#else
                                                                         amoebaGpu->psWorkArray_3_2->_pDevStream[0] );
#endif

    } else {


//                                                                         gpu->sim.pInteractingWorkUnit,
//                                                                         amoebaGpu->psWorkUnit->_pDevStream[0],
      kCalculateAmoebaPmeDirectElectrostaticN2Forces_kernel<<<amoebaGpu->nonbondBlocks, threadsPerBlock, (sizeof(PmeDirectElectrostaticParticle)+sizeof(float3))*threadsPerBlock>>>(
                                                                         gpu->sim.pInteractingWorkUnit,
                                                                         amoebaGpu->psWorkArray_3_1->_pDevStream[0],
#ifdef AMOEBA_DEBUG
                                                                         amoebaGpu->psWorkArray_3_2->_pDevStream[0],
                                                                         debugArray->_pDevStream[0], targetAtom );
#else
                                                                         amoebaGpu->psWorkArray_3_2->_pDevStream[0] );
#endif
    }
    LAUNCHERROR("kCalculateAmoebaPmeDirectElectrostaticN2Forces");
    kReduceForceTorque( amoebaGpu );

#ifdef AMOEBA_DEBUG
    if( amoebaGpu->log ){

        amoebaGpu->psForce->Download();
        amoebaGpu->psTorque->Download();
        debugArray->Download();
  
        (void) fprintf( amoebaGpu->log, "Finished PmeDirectElectrostatic kernel execution\n" ); (void) fflush( amoebaGpu->log );
  
        int maxPrint        = 5;
        float conversion    = 1.0f/41.84f;
        float forceSum[3]   = { 0.0f, 0.0f, 0.0f};
        for( int ii = 0; ii < gpu->natoms; ii++ ){
            (void) fprintf( amoebaGpu->log, "%5d ", ii); 
  
            int indexOffset     = ii*3;
      
            // force
  
            (void) fprintf( amoebaGpu->log,"PmeDirectElectrostaticF [%16.9e %16.9e %16.9e] ",
                            conversion*amoebaGpu->psForce->_pSysStream[0][indexOffset],
                            conversion*amoebaGpu->psForce->_pSysStream[0][indexOffset+1],
                            conversion*amoebaGpu->psForce->_pSysStream[0][indexOffset+2] );

            forceSum[0]         += amoebaGpu->psForce->_pSysStream[0][indexOffset];
            forceSum[1]         += amoebaGpu->psForce->_pSysStream[0][indexOffset+1];
            forceSum[2]         += amoebaGpu->psForce->_pSysStream[0][indexOffset+2];

            // torque
  
            (void) fprintf( amoebaGpu->log,"PmeDirectElectrostaticT [%16.9e %16.9e %16.9e] ",
                            conversion*amoebaGpu->psTorque->_pSysStream[0][indexOffset],
                            conversion*amoebaGpu->psTorque->_pSysStream[0][indexOffset+1],
                            conversion*amoebaGpu->psTorque->_pSysStream[0][indexOffset+2] );
  
            (void) fprintf( amoebaGpu->log,"\n" );
            if( ii == maxPrint && (gpu->natoms - maxPrint) > ii ){
                ii = gpu->natoms - maxPrint;
            }
        }
        (void) fflush( amoebaGpu->log );
        gpu->psEnergy->Download();
        double energy = 0.0;
        for( unsigned int ii = 0; ii < gpu->sim.energyOutputBuffers; ii++ ){
            if( (*gpu->psEnergy)[ii] !=  (*gpu->psEnergy)[ii] || (*gpu->psEnergy)[ii] == std::numeric_limits<double>::infinity() || (*gpu->psEnergy)[ii] == -std::numeric_limits<double>::infinity() ){
                (void) fprintf( amoebaGpu->log,"Energy nan at index=%d\n", ii );
            } else {
               energy += (*gpu->psEnergy)[ii];
            }
        }   
        (void) fprintf( amoebaGpu->log,"Force sums: [%16.9e %16.9e %16.9e] Energy=%16.9e\n", forceSum[0], forceSum[1], forceSum[2], energy );

        if( 0 ){
            (void) fprintf( amoebaGpu->log,"DebugElecAll\n" );
            int paddedNumberOfAtoms = amoebaGpu->gpuContext->sim.paddedNumberOfAtoms;
            for( int jj = 0; jj < gpu->natoms*gpu->natoms; jj++ ){
                if( fabs( debugArray->_pSysStream[0][jj].w - 111.0 ) < 1.0e-04 ){
                    int debugIndex = jj;
                    (void) fprintf( amoebaGpu->log,"%8d [%16.9e %16.9e %16.9e %16.9e] Enr11\n", jj,
                                    debugArray->_pSysStream[0][debugIndex].x, debugArray->_pSysStream[0][debugIndex].y,
                                    debugArray->_pSysStream[0][debugIndex].z, debugArray->_pSysStream[0][debugIndex].w );
                    debugIndex += paddedNumberOfAtoms;
                    (void) fprintf( amoebaGpu->log,"%8d [%16.9e %16.9e %16.9e %16.9e] Enr12\n", jj,
                                    debugArray->_pSysStream[0][debugIndex].x, debugArray->_pSysStream[0][debugIndex].y,
                                    debugArray->_pSysStream[0][debugIndex].z, debugArray->_pSysStream[0][debugIndex].w );
                }
            }
        }
        (void) fprintf( amoebaGpu->log,"\n" );

        if( 0 ){
            (void) fprintf( amoebaGpu->log,"DebugElec\n" );
            int paddedNumberOfAtoms = amoebaGpu->gpuContext->sim.paddedNumberOfAtoms;
            for( int jj = 0; jj < gpu->natoms; jj++ ){
                int debugIndex = jj;
                for( int kk = 0; kk < 6; kk++ ){
                    (void) fprintf( amoebaGpu->log,"%5d %5d [%16.9e %16.9e %16.9e %16.9e] E11\n", targetAtom, jj,
                                    debugArray->_pSysStream[0][debugIndex].x, debugArray->_pSysStream[0][debugIndex].y,
                                    debugArray->_pSysStream[0][debugIndex].z, debugArray->_pSysStream[0][debugIndex].w );
                    debugIndex += paddedNumberOfAtoms;
                }
                (void) fprintf( amoebaGpu->log,"\n" );
            }
        }
        (void) fflush( amoebaGpu->log );
  
        if( 1 ){
            std::vector<int> fileId;
            //fileId.push_back( 0 );
            VectorOfDoubleVectors outputVector;
            cudaLoadCudaFloat4Array( gpu->natoms, 3, gpu->psPosq4,            outputVector );
            cudaLoadCudaFloatArray( gpu->natoms,  3, amoebaGpu->psForce,      outputVector );
            cudaLoadCudaFloatArray( gpu->natoms,  3, amoebaGpu->psTorque,     outputVector);
            cudaWriteVectorOfDoubleVectorsToFile( "CudaPmeDirectForceTorque", fileId, outputVector );
         }

    }   
    delete debugArray;
#endif

    cudaComputeAmoebaMapTorquesAndAddTotalForce( amoebaGpu, amoebaGpu->psTorque, amoebaGpu->psForce, gpu->psForce4 );

}

/**---------------------------------------------------------------------------------------

   Compute Amoeba electrostatic force & torque using PME

   @param amoebaGpu        amoebaGpu context

   --------------------------------------------------------------------------------------- */

void cudaComputeAmoebaPmeElectrostatic( amoebaGpuContext amoebaGpu )
{
    cudaComputeAmoebaPmeDirectElectrostatic( amoebaGpu );
    kCalculateAmoebaPMEInducedDipoleForces( amoebaGpu );
} 

