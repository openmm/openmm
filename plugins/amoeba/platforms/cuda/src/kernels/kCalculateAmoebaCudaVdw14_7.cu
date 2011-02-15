//-----------------------------------------------------------------------------------------

//-----------------------------------------------------------------------------------------

#include "amoebaGpuTypes.h"
#include "amoebaCudaKernels.h"
#include "cudaKernels.h"
#include "kCalculateAmoebaCudaUtilities.h"
#include "kCalculateAmoebaCudaVdwParticle.h"
#include "amoebaScaleFactors.h"

#include <stdio.h>
extern int isNanOrInfinity( double number );

using namespace std;

static __constant__ cudaGmxSimulation cSim;
static __constant__ cudaAmoebaGmxSimulation cAmoebaSim;

void SetCalculateAmoebaCudaVdw14_7Sim(amoebaGpuContext amoebaGpu)
{
    cudaError_t status;
    gpuContext gpu = amoebaGpu->gpuContext;
    status         = cudaMemcpyToSymbol(cSim, &gpu->sim, sizeof(cudaGmxSimulation));    
    RTERROR(status, "SetCalculateAmoebaCudaVdw14_7Sim: cudaMemcpyToSymbol: SetSim copy to cSim failed");
    status         = cudaMemcpyToSymbol(cAmoebaSim, &amoebaGpu->amoebaSim, sizeof(cudaAmoebaGmxSimulation));    
    RTERROR(status, "SetCalculateAmoebaCudaVdw14_7FieldSim: cudaMemcpyToSymbol: SetSim copy to cAmoebaSim failed");
}

void GetCalculateAmoebaCudaVdw14_7Sim(amoebaGpuContext amoebaGpu)
{
    cudaError_t status;
    gpuContext gpu = amoebaGpu->gpuContext;
    status = cudaMemcpyFromSymbol(&gpu->sim, cSim, sizeof(cudaGmxSimulation));    
    RTERROR(status, "GetCalculateAmoebaCudaVdw14_7Sim: cudaMemcpyFromSymbol: SetSim copy from cSim failed");
    status = cudaMemcpyFromSymbol(&amoebaGpu->amoebaSim, cAmoebaSim, sizeof(cudaAmoebaGmxSimulation));    
    RTERROR(status, "GetCalculateAmoebaCudaVdw14_7Sim: cudaMemcpyFromSymbol: SetSim copy from cAmoebaSim failed");
}

//#define AMOEBA_DEBUG_PRINT
#undef AMOEBA_DEBUG_PRINT

//#define AMOEBA_DEBUG
#undef AMOEBA_DEBUG

__device__ void zeroVdw14_7SharedForce( struct Vdw14_7Particle* sA ) 
{
    // zero shared fields

    sA->force[0]              = 0.0f;
    sA->force[1]              = 0.0f;
    sA->force[2]              = 0.0f;

}

__device__ void loadVdw14_7Shared( struct Vdw14_7Particle* sA, unsigned int atomI,
                                   float4* atomCoord, float2* vdwParameters )
{
    // coordinates, sigma, epsilon

    sA->x                        = atomCoord[atomI].x;
    sA->y                        = atomCoord[atomI].y;
    sA->z                        = atomCoord[atomI].z;

    sA->sigma                    = vdwParameters[atomI].x;
    sA->epsilon                  = vdwParameters[atomI].y;

}

__device__ void getVdw14_7CombindedSigmaEpsilon_kernel( int sigmaCombiningRule, float iSigma, float jSigma, float* combindedSigma,
                                                        int epsilonCombiningRule, float iEpsilon, float jEpsilon, float* combindedEpsilon )
{
    if( sigmaCombiningRule == 1 ){
        *combindedSigma      = iSigma + jSigma;
    } else if( sigmaCombiningRule == 2 ){
        *combindedSigma      = 2.0f*sqrtf( iSigma*jSigma );
    } else {
        float iSigma2        = iSigma*iSigma;
        float jSigma2        = jSigma*jSigma;
        *combindedSigma      = 2.0f*( iSigma2*iSigma + jSigma2*jSigma )/( iSigma2 + jSigma2 );
    }

    if( epsilonCombiningRule == 1 ){
        *combindedEpsilon    = iEpsilon + jEpsilon;
    } else if( epsilonCombiningRule == 2 ){
        *combindedEpsilon    = 2.0f*sqrtf( iEpsilon*jEpsilon );
    } else if( epsilonCombiningRule == 3 ){
        float iEpsilon2      = iEpsilon*iEpsilon;
        float jEpsilon2      = jEpsilon*jEpsilon;
        *combindedEpsilon    = 2.0f*( iEpsilon2*iEpsilon + jEpsilon2*jEpsilon )/( iEpsilon2 + jEpsilon2 );
    } else {
        float epsilonS       = sqrtf( iEpsilon ) + sqrtf( jEpsilon );
        *combindedEpsilon    = 4.0f*( iEpsilon*jEpsilon )/( epsilonS*epsilonS );
    }   

}

__device__ void calculateVdw14_7PairIxn_kernel( float combindedSigma,    float combindedEpsilon,
                                                float force[3], float* energy
#ifdef AMOEBA_DEBUG
               , float4* debugArray
#endif
)
{

    const float deltaHalM1 = 0.07f;
    const float deltaHal   = 1.07f;
    const float gammaHal   = 1.12f;
    const float gammaHalM1 = 0.12f;

    // ---------------------------------------------------------------------------------------
    
    // on input force[i] is assummed to contain delta[i] for coordinates of atom I and J 
    
    float r2                                     = force[0]*force[0] + force[1]*force[1] + force[2]*force[2];
    if( r2 > cAmoebaSim.vdwCutoff2 ){
        *energy = force[0] = force[1] = force[2] = 0.0f;
        return;
    }
    float rI                                     =  rsqrtf( r2 );
    float r                                      =  1.0f/rI;
    float r6                                     =  r2*r2*r2;
    float r7                                     =  r6*r;
 
    float combindedSigma7                        = combindedSigma*combindedSigma;
    combindedSigma7                              = combindedSigma7*combindedSigma7*combindedSigma7*combindedSigma;

    float rho                                    = r7 + combindedSigma7*gammaHalM1;
    float rhoInverse                             = 1.0f/rho;
 
    float tau                                    = deltaHal/(r + deltaHalM1*combindedSigma);
    float tau7                                   = tau*tau*tau;
         tau7                                    = tau7*tau7*tau;
    float dTau                                   = tau/deltaHal;
    
    float tmp                                    = combindedSigma7*rhoInverse;
    float gTau                                   = combindedEpsilon*tau7*r6*gammaHal*tmp*tmp;
 
    *energy                                      = combindedEpsilon*combindedSigma7*tau7*( (combindedSigma7*gammaHal*rhoInverse) - 2.0f);
    float deltaE                                 = (-7.0f*(dTau*(*energy) + gTau))*rI;
 
    force[0]                                    *= deltaE;
    force[1]                                    *= deltaE;
    force[2]                                    *= deltaE;

#ifdef AMOEBA_DEBUG
    debugArray[0].x                              = r;
    debugArray[0].y                              = deltaE;
    debugArray[0].z                              = combindedSigma;
    debugArray[0].w                              = combindedEpsilon;

    debugArray[1].x                              = tau;
    debugArray[1].y                              = rho;
    debugArray[1].z                              = gTau;
    debugArray[1].w                              = r;
#endif
}

// perform reduction of force on H's and add to heavy atom partner
// input force is the Vdw force
// output force is the cumulative force

__global__ 
#if (__CUDA_ARCH__ >= 200)
__launch_bounds__(GF1XX_NONBOND_THREADS_PER_BLOCK, 1)
#elif (__CUDA_ARCH__ >= 120)
__launch_bounds__(GT2XX_NONBOND_THREADS_PER_BLOCK, 1)
#else
__launch_bounds__(G8X_NONBOND_THREADS_PER_BLOCK, 1)
#endif
void kCalculateAmoebaVdw14_7Reduction_kernel( float* inputForce, float4* outputForce )
{
    unsigned int pos = blockIdx.x * blockDim.x + threadIdx.x;

    while (pos < cAmoebaSim.amoebaVdwReductions )
    {
        int4   atomIndices              = cAmoebaSim.pAmoebaVdwReductionID[pos];
        float4 forceA;
        float4 forceB;
        float4 forceC;
        float4 forceD;

        int index                       = 3*atomIndices.x;
        forceA.x                        = inputForce[index];
        forceA.y                        = inputForce[index+1];
        forceA.z                        = inputForce[index+2];

        index                           = 3*atomIndices.y;
        forceB.x                        = inputForce[index];
        forceB.y                        = inputForce[index+1];
        forceB.z                        = inputForce[index+2];

        index                           = 3*atomIndices.z;
        forceC.x                        = inputForce[index];
        forceC.y                        = inputForce[index+1];
        forceC.z                        = inputForce[index+2];

        index                           = 3*atomIndices.w;
        forceD.x                        = inputForce[index];
        forceD.y                        = inputForce[index+1];
        forceD.z                        = inputForce[index+2];

        float  reductionFactor          = cAmoebaSim.pAmoebaVdwReduction[pos];
        float  reductionFactorM1        = 1.0f - reductionFactor;
        
        float4 forceTemp1;
        forceTemp1.x                    = reductionFactor*forceB.x;
        forceTemp1.y                    = reductionFactor*forceB.y;
        forceTemp1.z                    = reductionFactor*forceB.z;

        forceA.x                       += reductionFactorM1*forceB.x;
        forceA.y                       += reductionFactorM1*forceB.y;
        forceA.z                       += reductionFactorM1*forceB.z;

        outputForce[atomIndices.y].x   += forceTemp1.x;
        outputForce[atomIndices.y].y   += forceTemp1.y;
        outputForce[atomIndices.y].z   += forceTemp1.z;

        reductionFactor                 = atomIndices.x != atomIndices.z ? reductionFactor   : 0.0f;
        reductionFactorM1               = atomIndices.x != atomIndices.z ? reductionFactorM1 : 0.0f;

        forceTemp1.x                    = reductionFactor*forceC.x;
        forceTemp1.y                    = reductionFactor*forceC.y;
        forceTemp1.z                    = reductionFactor*forceC.z;

        forceA.x                       += reductionFactorM1*forceC.x;
        forceA.y                       += reductionFactorM1*forceC.y;
        forceA.z                       += reductionFactorM1*forceC.z;

        outputForce[atomIndices.z].x   += forceTemp1.x;
        outputForce[atomIndices.z].y   += forceTemp1.y;
        outputForce[atomIndices.z].z   += forceTemp1.z;

        reductionFactor                 = atomIndices.x != atomIndices.w ? reductionFactor   : 0.0f;
        reductionFactorM1               = atomIndices.x != atomIndices.w ? reductionFactorM1 : 0.0f;

        forceTemp1.x                    = reductionFactor*forceD.x;
        forceTemp1.y                    = reductionFactor*forceD.y;
        forceTemp1.z                    = reductionFactor*forceD.z;

        forceA.x                       += reductionFactorM1*forceD.x;
        forceA.y                       += reductionFactorM1*forceD.y;
        forceA.z                       += reductionFactorM1*forceD.z;

        outputForce[atomIndices.w].x   += forceTemp1.x;
        outputForce[atomIndices.w].y   += forceTemp1.y;
        outputForce[atomIndices.w].z   += forceTemp1.z;

        outputForce[atomIndices.x].x   += forceA.x;
        outputForce[atomIndices.x].y   += forceA.y;
        outputForce[atomIndices.x].z   += forceA.z;
        
        pos                            += blockDim.x * gridDim.x;
    }
}

__device__ void sumTempBuffer( Vdw14_7Particle& atomI, Vdw14_7Particle& atomJ ){

    atomI.tempForce[0]  += atomJ.tempForce[0];
    atomI.tempForce[1]  += atomJ.tempForce[1];
    atomI.tempForce[2]  += atomJ.tempForce[2];
}

static void kCalculateAmoebaVdw14_7Reduction(amoebaGpuContext amoebaGpu, CUDAStream<float>* vdwOutputArray, CUDAStream<float4>* forceOutputArray )
{
    unsigned int threadsPerBlock;
    if (amoebaGpu->gpuContext->sm_version >= SM_20)
        threadsPerBlock = GF1XX_NONBOND_THREADS_PER_BLOCK;
    else if (amoebaGpu->gpuContext->sm_version >= SM_12)
        threadsPerBlock = GT2XX_NONBOND_THREADS_PER_BLOCK;
    else
        threadsPerBlock = G8X_NONBOND_THREADS_PER_BLOCK;
    kCalculateAmoebaVdw14_7Reduction_kernel<<<amoebaGpu->gpuContext->sim.blocks, threadsPerBlock>>>(
                               vdwOutputArray->_pDevData, forceOutputArray->_pDevData );
    LAUNCHERROR("kCalculateAmoebaVdw14_7Reduction");
}

// perform reduction of coordinate on H's and add to heavy atom partner
// input coordinate is the Vdw coordinate
// output coordinate is the cumulative coordinate

__global__ 
#if (__CUDA_ARCH__ >= 200)
__launch_bounds__(GF1XX_THREADS_PER_BLOCK, 1)
#elif (__CUDA_ARCH__ >= 120)
__launch_bounds__(GT2XX_THREADS_PER_BLOCK, 1)
#else
__launch_bounds__(G8X_THREADS_PER_BLOCK, 1)
#endif
void kCalculateAmoebaVdw14_7CoordinateReduction_kernel( float4* inputCoordinate, float4* outputCoordinate )
{
    unsigned int pos = blockIdx.x * blockDim.x + threadIdx.x;

    while (pos < cAmoebaSim.amoebaVdwReductions )
    {
        int4   atomIndices              = cAmoebaSim.pAmoebaVdwReductionID[pos];
        float4 coordinateA;
        float4 coordinateB;
        float4 coordinateC;
        float4 coordinateD;

        int index                             = atomIndices.x;
        coordinateA.x                         = inputCoordinate[index].x;
        coordinateA.y                         = inputCoordinate[index].y;
        coordinateA.z                         = inputCoordinate[index].z;

        index                                 = atomIndices.y;
        coordinateB.x                         = inputCoordinate[index].x;
        coordinateB.y                         = inputCoordinate[index].y;
        coordinateB.z                         = inputCoordinate[index].z;

        index                                 = atomIndices.z;
        coordinateC.x                         = inputCoordinate[index].x;
        coordinateC.y                         = inputCoordinate[index].y;
        coordinateC.z                         = inputCoordinate[index].z;

        index                                 = atomIndices.w;
        coordinateD.x                         = inputCoordinate[index].x;
        coordinateD.y                         = inputCoordinate[index].y;
        coordinateD.z                         = inputCoordinate[index].z;

        float  reductionFactor                = cAmoebaSim.pAmoebaVdwReduction[pos];
        float  reductionFactorM1              = 1.0f - reductionFactor;
        
        float4 coordinateTemp1;
        coordinateTemp1.x                     = reductionFactor*coordinateB.x + reductionFactorM1*coordinateA.x;
        coordinateTemp1.y                     = reductionFactor*coordinateB.y + reductionFactorM1*coordinateA.y;
        coordinateTemp1.z                     = reductionFactor*coordinateB.z + reductionFactorM1*coordinateA.z;

        outputCoordinate[atomIndices.y].x     = coordinateTemp1.x;
        outputCoordinate[atomIndices.y].y     = coordinateTemp1.y;
        outputCoordinate[atomIndices.y].z     = coordinateTemp1.z;

        reductionFactor                       = atomIndices.x != atomIndices.z ? reductionFactor   : 1.0f;
        reductionFactorM1                     = atomIndices.x != atomIndices.z ? reductionFactorM1 : 0.0f;

        coordinateTemp1.x                     = reductionFactor*coordinateC.x + reductionFactorM1*coordinateA.x;
        coordinateTemp1.y                     = reductionFactor*coordinateC.y + reductionFactorM1*coordinateA.y;
        coordinateTemp1.z                     = reductionFactor*coordinateC.z + reductionFactorM1*coordinateA.z;

        outputCoordinate[atomIndices.z].x     = coordinateTemp1.x;
        outputCoordinate[atomIndices.z].y     = coordinateTemp1.y;
        outputCoordinate[atomIndices.z].z     = coordinateTemp1.z;

        reductionFactor                       = atomIndices.x != atomIndices.w ? reductionFactor   : 1.0f;
        reductionFactorM1                     = atomIndices.x != atomIndices.w ? reductionFactorM1 : 0.0f;

        coordinateTemp1.x                     = reductionFactor*coordinateD.x + reductionFactorM1*coordinateA.x;
        coordinateTemp1.y                     = reductionFactor*coordinateD.y + reductionFactorM1*coordinateA.y;
        coordinateTemp1.z                     = reductionFactor*coordinateD.z + reductionFactorM1*coordinateA.z;

        outputCoordinate[atomIndices.w].x     = coordinateTemp1.x;
        outputCoordinate[atomIndices.w].y     = coordinateTemp1.y;
        outputCoordinate[atomIndices.w].z     = coordinateTemp1.z;

        pos                                  += blockDim.x * gridDim.x;
    }
}

static void kCalculateAmoebaVdw14_7CoordinateReduction(amoebaGpuContext amoebaGpu,
                                                       CUDAStream<float4>* coordinateArray,
                                                       CUDAStream<float4>* reducedCoordinateArray)
{
    unsigned int threadsPerBlock;
    if (amoebaGpu->gpuContext->sm_version >= SM_20)
        threadsPerBlock = GF1XX_THREADS_PER_BLOCK;
    else if (amoebaGpu->gpuContext->sm_version >= SM_12)
        threadsPerBlock = GT2XX_THREADS_PER_BLOCK;
    else
        threadsPerBlock = G8X_THREADS_PER_BLOCK;
    kCalculateAmoebaVdw14_7CoordinateReduction_kernel<<<amoebaGpu->gpuContext->sim.blocks, threadsPerBlock>>>(
                               coordinateArray->_pDevData, reducedCoordinateArray->_pDevData );
    LAUNCHERROR("kCalculateAmoebaVdw14_7CoordinateReduction");
}

// perform reduction of force on H's and add to heavy atom partner
// input force is the Vdw force
// output force is the cumulative force

__global__ 
#if (__CUDA_ARCH__ >= 200)
__launch_bounds__(GF1XX_THREADS_PER_BLOCK, 1)
#elif (__CUDA_ARCH__ >= 120)
__launch_bounds__(GT2XX_THREADS_PER_BLOCK, 1)
#else
__launch_bounds__(G8X_THREADS_PER_BLOCK, 1)
#endif
void kCalculateAmoebaVdw14_7NonReduction_kernel( float* inputForce, float4* outputForce )
{
    unsigned int pos = blockIdx.x * blockDim.x + threadIdx.x;

    while (pos < cAmoebaSim.amoebaVdwNonReductions )
    {
        int   atomIndex             = cAmoebaSim.pAmoebaVdwNonReductionID[pos];
        int index                   = 3*atomIndex;
        outputForce[atomIndex].x   += inputForce[index];
        outputForce[atomIndex].y   += inputForce[index+1];
        outputForce[atomIndex].z   += inputForce[index+2];
        
        pos                        += blockDim.x * gridDim.x;
    }
}

static void kCalculateAmoebaVdw14_7NonReduction(amoebaGpuContext amoebaGpu, CUDAStream<float>* vdwOutputArray, CUDAStream<float4>* forceOutputArray )
{
    kCalculateAmoebaVdw14_7NonReduction_kernel<<<amoebaGpu->gpuContext->sim.blocks, 384>>>(
                               vdwOutputArray->_pDevData, forceOutputArray->_pDevData );
    LAUNCHERROR("kCalculateAmoebaVdw14_7MonReduction");
}

// Include versions of the kernels for N^2 calculations.

#define METHOD_NAME(a, b) a##N2##b
#include "kCalculateAmoebaCudaVdw14_7.h"

#define USE_OUTPUT_BUFFER_PER_WARP
#undef METHOD_NAME
#define METHOD_NAME(a, b) a##N2ByWarp##b
#include "kCalculateAmoebaCudaVdw14_7.h"
#undef METHOD_NAME
#undef USE_OUTPUT_BUFFER_PER_WARP 

#define USE_CUTOFF
#define METHOD_NAME(a, b) a##Cutoff##b
#include "kCalculateAmoebaCudaVdw14_7.h"
#undef METHOD_NAME

#define USE_OUTPUT_BUFFER_PER_WARP
#undef METHOD_NAME
#define METHOD_NAME(a, b) a##CutoffByWarp##b
#include "kCalculateAmoebaCudaVdw14_7.h"
#undef METHOD_NAME
#undef USE_OUTPUT_BUFFER_PER_WARP 
#undef USE_CUTOFF

// reduce psWorkArray_3_1 -> outputArray

static void kReduceVdw14_7(amoebaGpuContext amoebaGpu, CUDAStream<float>* outputArray )
{
    kReduceFields_kernel<<<amoebaGpu->nonbondBlocks, amoebaGpu->fieldReduceThreadsPerBlock>>>(
                           amoebaGpu->paddedNumberOfAtoms*3, amoebaGpu->outputBuffers,
                           amoebaGpu->psWorkArray_3_1->_pDevData, outputArray->_pDevData );
    LAUNCHERROR("kReduceVdw14_7");
}

__global__ 
#if (__CUDA_ARCH__ >= 200)
__launch_bounds__(GF1XX_THREADS_PER_BLOCK, 1)
#elif (__CUDA_ARCH__ >= 120)
__launch_bounds__(GT2XX_THREADS_PER_BLOCK, 1)
#else
__launch_bounds__(G8X_THREADS_PER_BLOCK, 1)
#endif
void kCalculateAmoebaVdw14_7CopyCoordinates_kernel( unsigned int bufferLength, float4* toCopy, float4* copy )
{
    unsigned int pos = blockIdx.x * blockDim.x + threadIdx.x;
    while (pos < bufferLength )
    {   
        copy[pos].x       = toCopy[pos].x;
        copy[pos].y       = toCopy[pos].y;
        copy[pos].z       = toCopy[pos].z;
        copy[pos].w       = toCopy[pos].w;
        pos              += gridDim.x * blockDim.x;
    }   
}

void kCalculateAmoebaVdw14_7CopyCoordinates( amoebaGpuContext amoebaGpu, CUDAStream<float4>* toCopy, CUDAStream<float4>* copy )
{
    kCalculateAmoebaVdw14_7CopyCoordinates_kernel<<<amoebaGpu->gpuContext->blocksPerSM, 384>>>( amoebaGpu->gpuContext->sim.paddedNumberOfAtoms, 
                      toCopy->_pDevData, copy->_pDevData );
    LAUNCHERROR("kCalculateAmoebaVdw14_7CopyCoordinates");
}

/**---------------------------------------------------------------------------------------

   Compute Vdw 14-7

   @param amoebaGpu        amoebaGpu context

   --------------------------------------------------------------------------------------- */

void kCalculateAmoebaVdw14_7Forces( amoebaGpuContext amoebaGpu, int applyCutoff )
{
  
   // ---------------------------------------------------------------------------------------

    static int threadsPerBlock    = 0;

   // ---------------------------------------------------------------------------------------

     gpuContext gpu    = amoebaGpu->gpuContext;

#ifdef AMOEBA_DEBUG_PRINT
    static const char* methodName = "kCalculateAmoebaVdw14_7Forces";
    if( 1 && amoebaGpu->log ){
        (void) fprintf( amoebaGpu->log, "%s: \n", methodName );
        (void) fflush( amoebaGpu->log );
    }   
#ifdef AMOEBA_DEBUG
    int paddedNumberOfAtoms                    = amoebaGpu->gpuContext->sim.paddedNumberOfAtoms;
    CUDAStream<float4>* debugArray             = new CUDAStream<float4>(paddedNumberOfAtoms*paddedNumberOfAtoms, 1, "DebugArray");
    memset( debugArray->_pSysData,      0, sizeof( float )*4*paddedNumberOfAtoms*paddedNumberOfAtoms);
    debugArray->Upload();
    int targetAtom                             = 342;
#endif
#endif

    // set threads/block first time through

    // on first pass, set threads/block

    if( threadsPerBlock == 0 ){ 
        unsigned int maxThreads;
        if (gpu->sm_version >= SM_20)
            maxThreads = 384; 
        else if (gpu->sm_version >= SM_12)
            maxThreads = 192; 
        else
            maxThreads = 128;
        threadsPerBlock = std::min(getThreadsPerBlock(amoebaGpu, sizeof(Vdw14_7Particle)), maxThreads);
    }    

    kCalculateAmoebaVdw14_7CopyCoordinates( amoebaGpu, gpu->psPosq4, amoebaGpu->psAmoebaVdwCoordinates );
    kCalculateAmoebaVdw14_7CoordinateReduction( amoebaGpu, amoebaGpu->psAmoebaVdwCoordinates, amoebaGpu->psAmoebaVdwCoordinates );

#ifdef AMOEBA_DEBUG_PRINT
    (void) fprintf( amoebaGpu->log, "Apply cutoff=%d warp=%d\n", applyCutoff, gpu->bOutputBufferPerWarp );
    (void) fprintf( amoebaGpu->log, "numBlocks=%u numThreads=%u bufferPerWarp=%u atm=%u shrd=%u Ebuf=%u ixnCt=%u workUnits=%u\n",
                    amoebaGpu->nonbondBlocks, threadsPerBlock, amoebaGpu->bOutputBufferPerWarp,
                    sizeof(Vdw14_7Particle), sizeof(Vdw14_7Particle)*threadsPerBlock,
                    amoebaGpu->energyOutputBuffers, (*gpu->psInteractionCount)[0], gpu->sim.workUnits );
    if( 0 ){  
        gpu->psInteractionCount->Download();
        amoebaGpu->psVdwWorkUnit->Download();
        unsigned int totalWarps  = (amoebaGpu->nonbondBlocks*threadsPerBlock)/GRID;
        float        ratiof     = (float)totalWarps/(float)amoebaGpu->psVdwWorkUnit->_length;
        (void) fprintf( amoebaGpu->log, "Ixn warps=%u count=%u\n", totalWarps, gpu->psInteractionCount->_pSysData[0] );
        for( unsigned int ii = 0; ii < amoebaGpu->psVdwWorkUnit->_length; ii++ ){
    
            unsigned int x          = amoebaGpu->psVdwWorkUnit->_pSysData[ii];
            unsigned int y          = ((x >> 2) & 0x7fff) << GRIDBITS;
            unsigned int exclusions = (x & 0x1);
                         x          = (x >> 17) << GRIDBITS;
            float        warp       = (float)(ii)*ratiof; 
            (void) fprintf( amoebaGpu->log, "GpuCell %8u  [%5u %5u %1u] %10u warp=%15.6f\n", ii, x,y,exclusions, warp );
    
        }           
    }
    
    (void) fflush( amoebaGpu->log );
#endif

    // clear output arrays

    kClearFields_3( amoebaGpu, 1 );

    if( applyCutoff ){

        kFindBlockBoundsPeriodic_kernel<<<(gpu->psGridBoundingBox->_length+63)/64, 64>>>();
        LAUNCHERROR("kFindBlockBoundsVdwPeriodic");
        kFindBlocksWithInteractionsVdwPeriodic_kernel<<<gpu->sim.interaction_blocks, gpu->sim.interaction_threads_per_block>>>();
        LAUNCHERROR("kFindBlocksWithInteractionsVdwPeriodic");
        compactStream(gpu->compactPlan, gpu->sim.pInteractingWorkUnit, amoebaGpu->amoebaSim.pVdwWorkUnit, gpu->sim.pInteractionFlag, gpu->sim.workUnits, gpu->sim.pInteractionCount);
        kFindInteractionsWithinBlocksVdwPeriodic_kernel<<<gpu->sim.nonbond_blocks, gpu->sim.nonbond_threads_per_block,
                    sizeof(unsigned int)*gpu->sim.nonbond_threads_per_block>>>(gpu->sim.pInteractingWorkUnit);
        LAUNCHERROR("kFindInteractionsWithinBlocksVdwPeriodic");

        if( 0 ){  
            gpu->psInteractionCount->Download();
            gpu->psInteractingWorkUnit->Download();
            gpu->psInteractionFlag->Download();
            amoebaGpu->psVdwWorkUnit->Download();
            (void) fprintf( amoebaGpu->log, "Vdw Ixn count=%u\n", gpu->psInteractionCount->_pSysData[0] );
            for( unsigned int ii = 0; ii < gpu->psInteractingWorkUnit->_length; ii++ ){
        
                unsigned int x          = gpu->psInteractingWorkUnit->_pSysData[ii];
                unsigned int y          = ((x >> 2) & 0x7fff) << GRIDBITS;
                unsigned int exclusions = (x & 0x1);
                             x          = (x >> 17) << GRIDBITS;
                (void) fprintf( amoebaGpu->log, "GpuCell %8u  %8u [%5u %5u %1u] %10u ", ii, gpu->psInteractingWorkUnit->_pSysData[ii], x,y,exclusions, gpu->psInteractionFlag->_pSysData[ii] );
        
                             x          = amoebaGpu->psVdwWorkUnit->_pSysData[ii];
                             y          = ((x >> 2) & 0x7fff) << GRIDBITS;
                             exclusions = (x & 0x1);
                             x          = (x >> 17) << GRIDBITS;
                (void) fprintf( amoebaGpu->log, "   AmGpu %8u [%5u %5u %1u]\n", amoebaGpu->psWorkUnit->_pSysData[ii], x,y,exclusions );
            }    
                (void) fflush( amoebaGpu->log );
        }

        if (gpu->bOutputBufferPerWarp){
            kCalculateAmoebaVdw14_7CutoffByWarp_kernel<<<amoebaGpu->nonbondBlocks, threadsPerBlock, sizeof(Vdw14_7Particle)*threadsPerBlock>>>(
                                                                gpu->sim.pInteractingWorkUnit,
                                                                amoebaGpu->psAmoebaVdwCoordinates->_pDevData,
                                                                amoebaGpu->psVdwSigmaEpsilon->_pDevData,
                                                                amoebaGpu->vdwSigmaCombiningRule,
                                                                amoebaGpu->vdwEpsilonCombiningRule,
#ifdef AMOEBA_DEBUG
                                                                amoebaGpu->psWorkArray_3_1->_pDevData,
                                                                debugArray->_pDevData, targetAtom );
#else
                                                                amoebaGpu->psWorkArray_3_1->_pDevData );
#endif
        } else {

            kCalculateAmoebaVdw14_7Cutoff_kernel<<<amoebaGpu->nonbondBlocks, threadsPerBlock, sizeof(Vdw14_7Particle)*threadsPerBlock>>>(
                                                                gpu->sim.pInteractingWorkUnit,
                                                                amoebaGpu->psAmoebaVdwCoordinates->_pDevData,
                                                                amoebaGpu->psVdwSigmaEpsilon->_pDevData,
                                                                amoebaGpu->vdwSigmaCombiningRule,
                                                                amoebaGpu->vdwEpsilonCombiningRule,
#ifdef AMOEBA_DEBUG
                                                                amoebaGpu->psWorkArray_3_1->_pDevData,
                                                                debugArray->_pDevData, targetAtom );
#else
                                                                amoebaGpu->psWorkArray_3_1->_pDevData );
#endif
    
        }
        LAUNCHERROR("kCalculateAmoebaVdw14_7Cutoff");  

    } else {

        if (gpu->bOutputBufferPerWarp){

            kCalculateAmoebaVdw14_7N2ByWarp_kernel<<<amoebaGpu->nonbondBlocks, threadsPerBlock, sizeof(Vdw14_7Particle)*threadsPerBlock>>>(
                                                                amoebaGpu->psVdwWorkUnit->_pDevData,
                                                                amoebaGpu->psAmoebaVdwCoordinates->_pDevData,
                                                                amoebaGpu->psVdwSigmaEpsilon->_pDevData,
                                                                amoebaGpu->vdwSigmaCombiningRule,
                                                                amoebaGpu->vdwEpsilonCombiningRule,
#ifdef AMOEBA_DEBUG
                                                                amoebaGpu->psWorkArray_3_1->_pDevData,
                                                                debugArray->_pDevData, targetAtom );
#else
                                                                amoebaGpu->psWorkArray_3_1->_pDevData );
#endif
        } else {

            kCalculateAmoebaVdw14_7N2_kernel<<<amoebaGpu->nonbondBlocks, threadsPerBlock, sizeof(Vdw14_7Particle)*threadsPerBlock>>>(
                                                                amoebaGpu->psVdwWorkUnit->_pDevData,
                                                                amoebaGpu->psAmoebaVdwCoordinates->_pDevData,
                                                                amoebaGpu->psVdwSigmaEpsilon->_pDevData,
                                                                amoebaGpu->vdwSigmaCombiningRule,
                                                                amoebaGpu->vdwEpsilonCombiningRule,
#ifdef AMOEBA_DEBUG
                                                                amoebaGpu->psWorkArray_3_1->_pDevData,
                                                                debugArray->_pDevData, targetAtom );
#else
                                                                amoebaGpu->psWorkArray_3_1->_pDevData );
#endif
    
        }
        LAUNCHERROR("kCalculateAmoebaVdw14_7N2");  
    }

#ifdef AMOEBA_DEBUG_PRINT
    if( amoebaGpu->log ){

        (void) fprintf( amoebaGpu->log, "Finished 14-7 kernel execution\n" );
        (void) fflush( amoebaGpu->log );

#ifdef AMOEBA_DEBUG
        debugArray->Download();
        int paddedNumberOfAtoms          = amoebaGpu->gpuContext->sim.paddedNumberOfAtoms;
        double cutOff                    = 1.0e+03;
        for( int jj = 0; jj < gpu->natoms; jj++ ){
            int debugIndex = jj;
            (void) fprintf( amoebaGpu->log,"%5d %5d DebugVdw\n", targetAtom, jj );
            for( int kk = 0; kk < 5; kk++ ){
                (void) fprintf( amoebaGpu->log,"[%16.9e %16.9e %16.9e %16.9e]\n",
                                debugArray->_pSysData[debugIndex].x, debugArray->_pSysData[debugIndex].y,
                                debugArray->_pSysData[debugIndex].z, debugArray->_pSysData[debugIndex].w );
                if( kk == 4 && ( fabs(  debugArray->_pSysData[debugIndex].x ) > cutOff ||
                                 fabs(  debugArray->_pSysData[debugIndex].y ) > cutOff ||
                                 fabs(  debugArray->_pSysData[debugIndex].z ) > cutOff ) ){
                    (void) fprintf( amoebaGpu->log," XXXX\n" );
                }
                debugIndex += paddedNumberOfAtoms;
            }
            (void) fprintf( amoebaGpu->log,"\n" );
        }
#endif

        amoebaGpu->psWorkArray_3_2->Download();
        amoebaGpu->psWorkArray_3_1->Download();
        //for( int jj = 0; jj < 3*gpu->natoms; jj += 3 )
        for( int jj = 0; jj < 3*gpu->natoms; jj += 3 ){
            for( int kk = 0; kk < amoebaGpu->outputBuffers; kk++ ){
                float delta = fabs(amoebaGpu->psWorkArray_3_1->_pSysStream[kk][jj+2] + 1.0f);
                if( delta < 5.0e-06 || isNanOrInfinity( (double) amoebaGpu->psWorkArray_3_1->_pSysStream[kk][jj] ) || isNanOrInfinity( (double) amoebaGpu->psWorkArray_3_1->_pSysStream[kk][jj+2] )  )
                (void) fprintf( amoebaGpu->log,"%6d %6d [%16.9e %16.9e %16.9e] [%16.9e %16.9e %16.9e]\n", jj, kk,
                                amoebaGpu->psWorkArray_3_1->_pSysStream[kk][jj],
                                amoebaGpu->psWorkArray_3_1->_pSysStream[kk][jj+1],
                                amoebaGpu->psWorkArray_3_1->_pSysStream[kk][jj+2],
                                amoebaGpu->psWorkArray_3_2->_pSysStream[kk][jj],
                                amoebaGpu->psWorkArray_3_2->_pSysStream[kk][jj+1],
                                amoebaGpu->psWorkArray_3_2->_pSysStream[kk][jj+2] );
            }
        }

    }
#endif

    kReduceVdw14_7( amoebaGpu, amoebaGpu->psWorkArray_3_2 );
    kCalculateAmoebaVdw14_7Reduction( amoebaGpu, amoebaGpu->psWorkArray_3_2, amoebaGpu->gpuContext->psForce4 );
    kCalculateAmoebaVdw14_7NonReduction( amoebaGpu, amoebaGpu->psWorkArray_3_2, amoebaGpu->gpuContext->psForce4 );

    if( 0 ){
        int paddedNumberOfAtoms             = amoebaGpu->gpuContext->sim.paddedNumberOfAtoms;
        CUDAStream<float4>* psTempForce     = new CUDAStream<float4>(paddedNumberOfAtoms, 1, "psTempForce");
        kClearFloat4( amoebaGpu, paddedNumberOfAtoms, psTempForce );
        kCalculateAmoebaVdw14_7Reduction( amoebaGpu, amoebaGpu->psWorkArray_3_2, psTempForce );
        //kCalculateAmoebaVdw14_7NonReduction( amoebaGpu, amoebaGpu->psWorkArray_3_2, psTempForce );
        std::vector<int> fileId;
        //fileId.push_back( 0 );
        VectorOfDoubleVectors outputVector;
        cudaLoadCudaFloat4Array( gpu->natoms, 3, gpu->psPosq4, outputVector, gpu->psAtomIndex->_pSysData, 1.0f );
        cudaLoadCudaFloat4Array( gpu->natoms,  3, psTempForce, outputVector, gpu->psAtomIndex->_pSysData, 1.0f );
        cudaWriteVectorOfDoubleVectorsToFile( "CudaVdw", fileId, outputVector );
        delete psTempForce;
        //exit(0);
     }

#ifdef AMOEBA_DEBUG
    delete debugArray;
#endif

   // ---------------------------------------------------------------------------------------
}
