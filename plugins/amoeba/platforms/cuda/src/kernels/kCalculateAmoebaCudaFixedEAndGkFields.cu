
//-----------------------------------------------------------------------------------------

//-----------------------------------------------------------------------------------------

#include "amoebaCudaKernels.h"
#include "kCalculateAmoebaCudaUtilities.h"

//#define AMOEBA_DEBUG

static __constant__ cudaGmxSimulation cSim;
static __constant__ cudaAmoebaGmxSimulation cAmoebaSim;

void SetCalculateAmoebaCudaFixedEAndGKFieldsSim(amoebaGpuContext amoebaGpu)
{
    cudaError_t status;
    gpuContext gpu = amoebaGpu->gpuContext;
    status         = cudaMemcpyToSymbol(cSim, &gpu->sim, sizeof(cudaGmxSimulation));    
    RTERROR(status, "SetCalculateAmoebaCudaFixedEAndGKFieldSim: cudaMemcpyToSymbol: SetSim copy to cSim failed");
    status         = cudaMemcpyToSymbol(cAmoebaSim, &amoebaGpu->amoebaSim, sizeof(cudaAmoebaGmxSimulation));    
    RTERROR(status, "SetCalculateAmoebaCudaFixedEAndGKFieldSim: cudaMemcpyToSymbol: SetSim copy to cAmoebaSim failed");
}

void GetCalculateAmoebaCudaFixedEAndGKFieldSim(amoebaGpuContext amoebaGpu)
{
    cudaError_t status;
    gpuContext gpu = amoebaGpu->gpuContext;
    status = cudaMemcpyFromSymbol(&gpu->sim, cSim, sizeof(cudaGmxSimulation));    
    RTERROR(status, "GetCalculateAmoebaCudaFixedEAndGKFieldSim: cudaMemcpyFromSymbol: SetSim copy from cSim failed");
    status = cudaMemcpyFromSymbol(&amoebaGpu->amoebaSim, cAmoebaSim, sizeof(cudaAmoebaGmxSimulation));         
    RTERROR(status, "GetCalculateAmoebaCudaFixedEAndGKFieldSim: cudaMemcpyFromSymbol: SetSim copy from cAmoebaSim failed");
}

// reduce psWorkArray_3_1 -> E_Field
// reduce psWorkArray_3_2 -> E_FieldPolar
// reduce psWorkArray_3_3 -> Gk_FieldPolar

static void kReduceEAndGkFields(amoebaGpuContext amoebaGpu )
{

    kReduceFields_kernel<<<amoebaGpu->nonbondBlocks, amoebaGpu->fieldReduceThreadsPerBlock>>>(
                               amoebaGpu->paddedNumberOfAtoms*3, amoebaGpu->outputBuffers,
                               amoebaGpu->psWorkArray_3_1->_pDevData, amoebaGpu->psE_Field->_pDevData );

    LAUNCHERROR("kReduceEAndGK_Fields1");
    kReduceFields_kernel<<<amoebaGpu->nonbondBlocks, amoebaGpu->fieldReduceThreadsPerBlock>>>(
                               amoebaGpu->paddedNumberOfAtoms*3, amoebaGpu->outputBuffers,
                               amoebaGpu->psWorkArray_3_2->_pDevData, amoebaGpu->psE_FieldPolar->_pDevData );
    LAUNCHERROR("kReduceEAndGK_Fields2");

    kReduceFields_kernel<<<amoebaGpu->nonbondBlocks, amoebaGpu->fieldReduceThreadsPerBlock>>>(
                               amoebaGpu->paddedNumberOfAtoms*3, amoebaGpu->outputBuffers,
                               amoebaGpu->psWorkArray_3_3->_pDevData, amoebaGpu->psGk_Field->_pDevData );
    LAUNCHERROR("kReduceEAndGK_Fields3");
}

// file includes FixedFieldParticle struct definition/load/unload struct and kernel body for fixed E-field

#define GK
#include "kCalculateAmoebaCudaFixedFieldParticle.h"
#undef GK

__device__ void calculateFixedGkFieldPairIxn_kernel( float4 atomCoordinatesI,       float4 atomCoordinatesJ,
                                                     float* labFrameDipoleI,        float* labFrameDipoleJ,
                                                     float* labFrameQuadrupoleI,    float* labFrameQuadrupoleJ,
                                                     float  rb2,
                                                     float outputField[2][3]
#ifdef AMOEBA_DEBUG
                                          , float4* debugArray
#endif

 ){
  
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
    float r2;
    float fc,fd,fq;
    float expterm;
    float gf,gf2,gf3,gf5;
    float gf7;
    float expc,dexpc;
    float expc1,expcdexpc;
    float a[4][4];
    float gc[5];
    float gux[11],guy[11],guz[11];
    float gqxx[5],gqxy[5];
    float gqxz[5],gqyy[5];
    float gqyz[5],gqzz[5];

    float gkc;

    gkc          = cAmoebaSim.gkc;

    fc           = cAmoebaSim.fc;
    fd           = cAmoebaSim.fd;
    fq           = cAmoebaSim.fq;

    xi           = atomCoordinatesI.x;
    yi           = atomCoordinatesI.y;
    zi           = atomCoordinatesI.z;
    ci           = atomCoordinatesI.w;

    uxi          = labFrameDipoleI[0];
    uyi          = labFrameDipoleI[1];
    uzi          = labFrameDipoleI[2];

    qxxi         = labFrameQuadrupoleI[0];
    qxyi         = labFrameQuadrupoleI[1];
    qxzi         = labFrameQuadrupoleI[2];
    qyyi         = labFrameQuadrupoleI[4];
    qyzi         = labFrameQuadrupoleI[5];
    qzzi         = labFrameQuadrupoleI[8];

    xr           = atomCoordinatesJ.x - xi;
    yr           = atomCoordinatesJ.y - yi;
    zr           = atomCoordinatesJ.z - zi;
    ck           = atomCoordinatesJ.w;

    xr2          = xr*xr;
    yr2          = yr*yr;
    zr2          = zr*zr;
    r2           = xr2 + yr2 + zr2;

    uxk          = labFrameDipoleJ[0];
    uyk          = labFrameDipoleJ[1];
    uzk          = labFrameDipoleJ[2];

    qxxk         = labFrameQuadrupoleJ[0];
    qxyk         = labFrameQuadrupoleJ[1];
    qxzk         = labFrameQuadrupoleJ[2];
    qyyk         = labFrameQuadrupoleJ[4];
    qyzk         = labFrameQuadrupoleJ[5];
    qzzk         = labFrameQuadrupoleJ[8];

    expterm      = exp(-r2/(gkc*rb2));
    expc         = expterm / gkc;
    dexpc        = -2.0f / (gkc*rb2);
    gf2          = 1.0f / (r2+rb2*expterm);
    gf           = sqrtf(gf2);
    gf3          = gf2 * gf;
    gf5          = gf3 * gf2;
    gf7          = gf5 * gf2;

    // reaction potential auxiliary terms

    a[0][0]      = gf;
    a[1][0]      = -gf3;
    a[2][0]      = 3.0f * gf5;
    a[3][0]      = -15.0f * gf7;

    // reaction potential gradient auxiliary terms

    expc1        = 1.0f - expc;
    a[0][1]      = expc1 * a[1][0];
    a[1][1]      = expc1 * a[2][0];
    a[2][1]      = expc1 * a[3][0];

    // dipole second reaction potential gradient auxiliary term

    expcdexpc    = -expc * dexpc;
    a[1][2]      = expc1*a[2][1] + expcdexpc*a[2][0];

    // multiply the auxillary terms by dielectric functions;

    a[0][1]      = fc * a[0][1];
    a[1][0]      = fd * a[1][0];
    a[1][1]      = fd * a[1][1];
    a[1][2]      = fd * a[1][2];
    a[2][0]      = fq * a[2][0];
    a[2][1]      = fq * a[2][1];

    // unweighted dipole reaction potential tensor

    gux[1]       = xr * a[1][0];
    guy[1]       = yr * a[1][0];
    guz[1]       = zr * a[1][0];

    // unweighted reaction potential gradient tensor

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
    gqxx[3]      = yr * xr2*a[2][1];
    gqxx[4]      = zr * xr2*a[2][1];
    gqyy[2]      = xr * yr2*a[2][1];
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

    // unweighted dipole second reaction potential gradient tensor

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

    // generalized Kirkwood permanent reaction field

    outputField[0][0] = uxk*gux[2] + uyk*gux[3] + uzk*gux[4]
                                   + 0.5f * (ck*gux[1] + qxxk*gux[5]
                                   + qyyk*gux[8] + qzzk*gux[10]
                                   + 2.0f*(qxyk*gux[6]+qxzk*gux[7]
                                   + qyzk*gux[9]))
                                   + 0.5f * (ck*gc[2] + qxxk*gqxx[2]
                                   + qyyk*gqyy[2] + qzzk*gqzz[2]
                                   + 2.0f*(qxyk*gqxy[2]+qxzk*gqxz[2]
                                   + qyzk*gqyz[2]));

    outputField[0][1] = uxk*guy[2] + uyk*guy[3] + uzk*guy[4]
                                   + 0.5f * (ck*guy[1] + qxxk*guy[5]
                                   + qyyk*guy[8] + qzzk*guy[10]
                                   + 2.0f*(qxyk*guy[6]+qxzk*guy[7]
                                   + qyzk*guy[9]))
                                   + 0.5f * (ck*gc[3] + qxxk*gqxx[3]
                                   + qyyk*gqyy[3] + qzzk*gqzz[3]
                                   + 2.0f*(qxyk*gqxy[3]+qxzk*gqxz[3]
                                   + qyzk*gqyz[3]));

    outputField[0][2] = uxk*guz[2] + uyk*guz[3] + uzk*guz[4]
                                   + 0.5f * (ck*guz[1] + qxxk*guz[5]
                                   + qyyk*guz[8] + qzzk*guz[10]
                                   + 2.0f*(qxyk*guz[6]+qxzk*guz[7]
                                   + qyzk*guz[9]))
                                   + 0.5f * (ck*gc[4] + qxxk*gqxx[4]
                                   + qyyk*gqyy[4] + qzzk*gqzz[4]
                                   + 2.0f*(qxyk*gqxy[4]+qxzk*gqxz[4]
                                   + qyzk*gqyz[4]));

    outputField[1][0] = uxi*gux[2] + uyi*gux[3] + uzi*gux[4]
                                   - 0.5f * (ci*gux[1] + qxxi*gux[5]
                                   + qyyi*gux[8] + qzzi*gux[10]
                                   + 2.0f*(qxyi*gux[6]+qxzi*gux[7]
                                   + qyzi*gux[9]))
                                   - 0.5f * (ci*gc[2] + qxxi*gqxx[2]
                                   + qyyi*gqyy[2] + qzzi*gqzz[2]
                                   + 2.0f*(qxyi*gqxy[2]+qxzi*gqxz[2]
                                   + qyzi*gqyz[2]));

    outputField[1][1] = uxi*guy[2] + uyi*guy[3] + uzi*guy[4]
                                   - 0.5f * (ci*guy[1] + qxxi*guy[5]
                                   + qyyi*guy[8] + qzzi*guy[10]
                                   + 2.0f*(qxyi*guy[6]+qxzi*guy[7]
                                   + qyzi*guy[9]))
                                   - 0.5f * (ci*gc[3]      + qxxi*gqxx[3]
                                   + qyyi*gqyy[3] + qzzi*gqzz[3]
                                   + 2.0f*(qxyi*gqxy[3]+qxzi*gqxz[3]
                                   + qyzi*gqyz[3]));

    outputField[1][2] = uxi*guz[2] + uyi*guz[3] + uzi*guz[4]
                                   - 0.5f * (ci*guz[1] + qxxi*guz[5]
                                   + qyyi*guz[8] + qzzi*guz[10]
                                   + 2.0f*(qxyi*guz[6]+qxzi*guz[7]
                                   + qyzi*guz[9]))
                                   - 0.5f * (ci*gc[4] + qxxi*gqxx[4]
                                   + qyyi*gqyy[4] + qzzi*gqzz[4]
                                   + 2.0f*(qxyi*gqxy[4]+qxzi*gqxz[4]
                                   + qyzi*gqyz[4]));

}

// Include versions of the kernels for N^2 calculations.

#define METHOD_NAME(a, b) a##N2##b
#include "kCalculateAmoebaCudaFixedEAndGkFields.h"
#define USE_OUTPUT_BUFFER_PER_WARP
#undef METHOD_NAME
#define METHOD_NAME(a, b) a##N2ByWarp##b
#include "kCalculateAmoebaCudaFixedEAndGkFields.h"

#ifdef AMOEBA_DEBUG
#if 0
static void printEFieldBuffer( amoebaGpuContext amoebaGpu, unsigned int bufferIndex )
{
    (void) fprintf( amoebaGpu->log, "EField Buffer %u\n", bufferIndex );
    unsigned int start = bufferIndex*3*amoebaGpu->paddedNumberOfAtoms;
    unsigned int stop  = (bufferIndex+1)*3*amoebaGpu->paddedNumberOfAtoms;
    for( unsigned int ii = start; ii < stop; ii += 3 ){
        unsigned int ii3Index      = ii/3;
        unsigned int bufferIndex   = ii3Index/(amoebaGpu->paddedNumberOfAtoms);
        unsigned int particleIndex = ii3Index - bufferIndex*(amoebaGpu->paddedNumberOfAtoms);
        (void) fprintf( amoebaGpu->log, "   %6u %3u %6u [%14.6e %14.6e %14.6e] [%14.6e %14.6e %14.6e]\n", 
                            ii/3,  bufferIndex, particleIndex,
                            amoebaGpu->psWorkArray_3_1->_pSysData[ii],
                            amoebaGpu->psWorkArray_3_1->_pSysData[ii+1],
                            amoebaGpu->psWorkArray_3_1->_pSysData[ii+2],
                            amoebaGpu->psWorkArray_3_2->_pSysData[ii],
                            amoebaGpu->psWorkArray_3_2->_pSysData[ii+1],
                            amoebaGpu->psWorkArray_3_2->_pSysData[ii+2] );
    } 
}

static void printEFieldAtomBuffers( amoebaGpuContext amoebaGpu, unsigned int targetAtom )
{
    (void) fprintf( amoebaGpu->log, "EField atom %u\n", targetAtom );
    for( unsigned int ii = 0; ii < amoebaGpu->outputBuffers; ii++ ){
        unsigned int particleIndex = targetAtom + ii*3*amoebaGpu->paddedNumberOfAtoms;
        (void) fprintf( amoebaGpu->log, " %2u %6u [%14.6e %14.6e %14.6e] [%14.6e %14.6e %14.6e]\n", 
                        ii, particleIndex,
                        amoebaGpu->psWorkArray_3_1->_pSysData[particleIndex],
                        amoebaGpu->psWorkArray_3_1->_pSysData[particleIndex+1],
                        amoebaGpu->psWorkArray_3_1->_pSysData[particleIndex+2],
                        amoebaGpu->psWorkArray_3_2->_pSysData[particleIndex],
                        amoebaGpu->psWorkArray_3_2->_pSysData[particleIndex+1],
                        amoebaGpu->psWorkArray_3_2->_pSysData[particleIndex+2] );
    } 
}
#endif
#endif

/**---------------------------------------------------------------------------------------

   Compute fixed electric field

   @param amoebaGpu        amoebaGpu context
   @param gpu              OpenMM gpu Cuda context

   --------------------------------------------------------------------------------------- */

void cudaComputeAmoebaFixedEAndGkFields( amoebaGpuContext amoebaGpu )
{
  
   // ---------------------------------------------------------------------------------------

#ifdef AMOEBA_DEBUG
    static const char* methodName = "computeCudaAmoebaFixedEAndGKFields";
#endif

   // ---------------------------------------------------------------------------------------

    static unsigned int threadsPerBlock = 0;

    gpuContext gpu                             = amoebaGpu->gpuContext;

#ifdef AMOEBA_DEBUG
    if( amoebaGpu->log ){
        (void) fprintf( amoebaGpu->log, "\n%s\n", methodName ); (void) fflush( amoebaGpu->log );
    }
    int paddedNumberOfAtoms                    = amoebaGpu->gpuContext->sim.paddedNumberOfAtoms;

    // N2 debug array

    CUDAStream<float4>* debugArray             = new CUDAStream<float4>(paddedNumberOfAtoms*paddedNumberOfAtoms, 1, "DebugArray");
    memset( debugArray->_pSysData,      0, sizeof( float )*4*paddedNumberOfAtoms*paddedNumberOfAtoms);
    debugArray->Upload();

    (*gpu->psInteractionCount)[0]              = gpu->sim.workUnits;
    gpu->psInteractionCount->Upload();

    // print intermediate results for the targetAtom 

    unsigned int targetAtom  = 0;

#endif

    // on first pass, set threads/block

    if( threadsPerBlock == 0 ){
        unsigned int maxThreads;
        if (gpu->sm_version >= SM_20)
            maxThreads = 256;
        else if (gpu->sm_version >= SM_12)
            maxThreads = 128;
        else
            maxThreads = 64;
        threadsPerBlock = std::min(getThreadsPerBlock(amoebaGpu, sizeof(FixedFieldParticle)), maxThreads);
    }

    kClearFields_3( amoebaGpu, 3 );

    if (gpu->bOutputBufferPerWarp){
        (void) fprintf( amoebaGpu->log, "N2 warp\n" );
        kCalculateAmoebaFixedEAndGkFieldN2ByWarp_kernel<<<amoebaGpu->nonbondBlocks, threadsPerBlock, sizeof(FixedFieldParticle)*threadsPerBlock>>>(
                                                                           amoebaGpu->psWorkUnit->_pDevData,
                                                                           gpu->psPosq4->_pDevData,
                                                                           amoebaGpu->psLabFrameDipole->_pDevData,
                                                                           amoebaGpu->psLabFrameQuadrupole->_pDevData,
                                                                           gpu->psBornRadii->_pDevData,
                                                                           amoebaGpu->psWorkArray_3_1->_pDevData,
                                                                           amoebaGpu->psWorkArray_3_2->_pDevData,
#ifdef AMOEBA_DEBUG
                                                                           amoebaGpu->psWorkArray_3_3->_pDevData,
                                                                           debugArray->_pDevData, targetAtom );
#else
                                                                           amoebaGpu->psWorkArray_3_3->_pDevData );
#endif
    } else {

#ifdef AMOEBA_DEBUG
        (void) fprintf( amoebaGpu->log, "N2 no warp\n" );
        (void) fprintf( amoebaGpu->log, "AmoebaN2Forces_kernel numBlocks=%u numThreads=%u bufferPerWarp=%u atm=%u shrd=%u ixnCt=%u workUnits=%u\n",
                        amoebaGpu->nonbondBlocks, threadsPerBlock, amoebaGpu->bOutputBufferPerWarp,
                        sizeof(FixedFieldParticle), sizeof(FixedFieldParticle)*threadsPerBlock, (*gpu->psInteractionCount)[0], gpu->sim.workUnits );
        (void) fflush( amoebaGpu->log );
#endif

        kCalculateAmoebaFixedEAndGkFieldN2_kernel<<<amoebaGpu->nonbondBlocks, threadsPerBlock, sizeof(FixedFieldParticle)*threadsPerBlock>>>(
                                                          amoebaGpu->psWorkUnit->_pDevData,
                                                          gpu->psPosq4->_pDevData,
                                                          amoebaGpu->psLabFrameDipole->_pDevData,
                                                          amoebaGpu->psLabFrameQuadrupole->_pDevData,
                                                          gpu->psBornRadii->_pDevData,
                                                          amoebaGpu->psWorkArray_3_1->_pDevData,
                                                          amoebaGpu->psWorkArray_3_2->_pDevData,
#ifdef AMOEBA_DEBUG
                                                          amoebaGpu->psWorkArray_3_3->_pDevData,
                                                          debugArray->_pDevData, targetAtom );
#else
                                                          amoebaGpu->psWorkArray_3_3->_pDevData );
#endif
    }
    LAUNCHERROR("kCalculateAmoebaFixedE_FieldN2Forces_kernel");

#if 0
        for( unsigned int ii = 0; ii < amoebaGpu->outputBuffers; ii++ ){
            //float index = 1.0f;
            float index = (float) ii;
            for( unsigned int jj = 0; jj < 3*amoebaGpu->paddedNumberOfAtoms; jj += 3 ){
                unsigned int kk = 3*ii*amoebaGpu->paddedNumberOfAtoms + jj;
                amoebaGpu->psWorkArray_3_1->_pSysData[kk]   = index;
                amoebaGpu->psWorkArray_3_1->_pSysData[kk+1] = index;
                amoebaGpu->psWorkArray_3_1->_pSysData[kk+2] = index;
            }
        }
        amoebaGpu->psWorkArray_3_1->Upload();
#endif

    kReduceEAndGkFields( amoebaGpu );

#ifdef AMOEBA_DEBUG
    if( amoebaGpu->log ){
        gpu->psInteractionCount->Download();
        (void) fprintf( amoebaGpu->log, "AmoebaN2Forces_kernel numBlocks=%u numThreads=%u bufferPerWarp=%u atm=%u shrd=%u ixnCt=%u workUnits=%u\n",
                        amoebaGpu->nonbondBlocks, amoebaGpu->nonbondThreadsPerBlock, amoebaGpu->bOutputBufferPerWarp,
                        sizeof(FixedFieldParticle), sizeof(FixedFieldParticle)*amoebaGpu->nonbondThreadsPerBlock, (*gpu->psInteractionCount)[0], gpu->sim.workUnits );
        (void) fflush( amoebaGpu->log );

        amoebaGpu->psWorkArray_3_1->Download();
        amoebaGpu->psWorkArray_3_2->Download();
        amoebaGpu->psWorkArray_3_3->Download();

        amoebaGpu->psE_Field->Download();
        amoebaGpu->psE_FieldPolar->Download();
        amoebaGpu->psGk_Field->Download();
        gpu->psBornRadii->Download();

        (void) fprintf( amoebaGpu->log, "OutE & Gk Fields\n" );
        int maxPrint        = 32;

        for( int ii = 0; ii < gpu->natoms; ii++ ){
           (void) fprintf( amoebaGpu->log, "%5d ", ii); 

            int indexOffset     = ii*3;

           // E_Field

           (void) fprintf( amoebaGpu->log,"E[%16.9e %16.9e %16.9e] ",
                           amoebaGpu->psE_Field->_pSysData[indexOffset],
                           amoebaGpu->psE_Field->_pSysData[indexOffset+1],
                           amoebaGpu->psE_Field->_pSysData[indexOffset+2] );
   
           // E_Field polar

           (void) fprintf( amoebaGpu->log,"Epol[%16.9e %16.9e %16.9e] ",
                           amoebaGpu->psE_FieldPolar->_pSysData[indexOffset],
                           amoebaGpu->psE_FieldPolar->_pSysData[indexOffset+1],
                           amoebaGpu->psE_FieldPolar->_pSysData[indexOffset+2] );

           // Gk_Field polar

           (void) fprintf( amoebaGpu->log,"Gk[%16.9e %16.9e %16.9e] ",
                           amoebaGpu->psGk_Field->_pSysData[indexOffset],
                           amoebaGpu->psGk_Field->_pSysData[indexOffset+1],
                           amoebaGpu->psGk_Field->_pSysData[indexOffset+2] );

           (void) fprintf( amoebaGpu->log,"\n" );
           if( ii == maxPrint && (gpu->natoms - maxPrint) > ii ){
                ii = gpu->natoms - maxPrint;
           }
        }
        (void) fflush( amoebaGpu->log );

        //printEFieldAtomBuffers( amoebaGpu, (targetAtom + 0) );
        //printEFieldAtomBuffers( amoebaGpu, (targetAtom + 1) );
        //printEFieldAtomBuffers( amoebaGpu, 100 );
        //printEFieldBuffer( amoebaGpu, 0 );
        //printEFieldBuffer( amoebaGpu, 1 );
        //printEFieldBuffer( amoebaGpu, 37 );
        //printEFieldBuffer( amoebaGpu, 38 );

        (void) fprintf( amoebaGpu->log, "EFields End\n" );
        (void) fprintf( amoebaGpu->log, "DebugQ\n" );
        debugArray->Download();
        if( 0 ){
            int ii = targetAtom;
            (void) fprintf( amoebaGpu->log,"\n" );
            int paddedNumberOfAtoms                    = amoebaGpu->gpuContext->sim.paddedNumberOfAtoms;
            unsigned int count                         = 0;
            for( int jj = 0; jj < gpu->natoms; jj++ ){
                int debugIndex = jj; 
                (void) fprintf( amoebaGpu->log,"%4d %4d Qint [%16.9e %16.9e %16.9e %16.9e] %16.9e ",
                                   ii, jj, 
                                   debugArray->_pSysData[debugIndex].x, debugArray->_pSysData[debugIndex].y,
                                   debugArray->_pSysData[debugIndex].z, debugArray->_pSysData[debugIndex].w,
                                   gpu->psBornRadii->_pSysData[jj] );

                for( int kk = 0; kk < 2; kk++ ){
                    debugIndex += paddedNumberOfAtoms;
                    (void) fprintf( amoebaGpu->log,"[%16.9e %16.9e %16.9e] ",
                                    debugArray->_pSysData[debugIndex].x, debugArray->_pSysData[debugIndex].y, debugArray->_pSysData[debugIndex].z );
                }
                (void) fprintf( amoebaGpu->log,"\n" );
            }
        }

        // write results to file

        if( 1 ){
            std::vector<int> fileId;
            //fileId.push_back( 0 );
            VectorOfDoubleVectors outputVector;
            cudaLoadCudaFloat4Array( gpu->natoms, 3, gpu->psPosq4,              outputVector, NULL, 1.0f );
            cudaLoadCudaFloatArray( gpu->natoms,  3, amoebaGpu->psE_Field,      outputVector, NULL, 1.0f );
            cudaLoadCudaFloatArray( gpu->natoms,  3, amoebaGpu->psE_FieldPolar, outputVector, NULL, 1.0f);
            cudaLoadCudaFloatArray( gpu->natoms,  3, amoebaGpu->psGk_Field,     outputVector, NULL, 1.0f);
            cudaWriteVectorOfDoubleVectorsToFile( "CudaEAndGkField", fileId, outputVector );

         }
         delete debugArray;
    }
#endif
 
//exit(0);
}
