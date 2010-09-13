//-----------------------------------------------------------------------------------------

//-----------------------------------------------------------------------------------------

#include "amoebaGpuTypes.h"
#include "amoebaCudaKernels.h"
#include "kCalculateAmoebaCudaUtilities.h"

#include <stdio.h>

using namespace std;

static __constant__ cudaGmxSimulation cSim;
static __constant__ cudaAmoebaGmxSimulation cAmoebaSim;

void SetCalculateAmoebaCudaMutualInducedAndGkFieldsSim(amoebaGpuContext amoebaGpu)
{
    cudaError_t status;
    gpuContext gpu = amoebaGpu->gpuContext;
    status         = cudaMemcpyToSymbol(cSim, &gpu->sim, sizeof(cudaGmxSimulation));    
    RTERROR(status, "SetCalculateAmoebaCudaMutualInducedAndGkFieldSim: cudaMemcpyToSymbol: SetSim copy to cSim failed");
    status         = cudaMemcpyToSymbol(cAmoebaSim, &amoebaGpu->amoebaSim, sizeof(cudaAmoebaGmxSimulation));    
    RTERROR(status, "SetCalculateAmoebaCudaMutualInducedAndGkFieldSim: cudaMemcpyToSymbol: SetSim copy to cAmoebaSim failed");
}

void GetCalculateAmoebaCudaMutualInducedAndGkFieldsSim(amoebaGpuContext amoebaGpu)
{
    cudaError_t status;
    gpuContext gpu = amoebaGpu->gpuContext;
    status = cudaMemcpyFromSymbol(&gpu->sim, cSim, sizeof(cudaGmxSimulation));    
    RTERROR(status, "GetCalculateAmoebaCudaMutualInducedAndGkFieldSim: cudaMemcpyFromSymbol: SetSim copy from cSim failed");
    status = cudaMemcpyFromSymbol(&amoebaGpu->amoebaSim, cAmoebaSim, sizeof(cudaAmoebaGmxSimulation));    
    RTERROR(status, "GetCalculateAmoebaCudaMutualInducedAndGkFieldSim: cudaMemcpyFromSymbol: SetSim copy from cAmoebaSim failed");
}

//#define AMOEBA_DEBUG
#undef AMOEBA_DEBUG

#define GK
#include "kCalculateAmoebaCudaMutualInducedParticle.h"
#undef GK

__device__ void calculateMutualInducedAndGkFieldsPairIxn_kernel( MutualInducedParticle& atomI, MutualInducedParticle& atomJ,
                                                                 float fields[8][3]

#ifdef AMOEBA_DEBUG
               , float4* debugArray
#endif

 )
{

    float deltaR[3];
    
    // ---------------------------------------------------------------------------------------
    
    // get deltaR, and r between 2 atoms
    
    deltaR[0]                                    = atomJ.x - atomI.x;
    deltaR[1]                                    = atomJ.y - atomI.y;
    deltaR[2]                                    = atomJ.z - atomI.z;

    float r                                      =  sqrtf( deltaR[0]*deltaR[0] + deltaR[1]*deltaR[1] + deltaR[2]*deltaR[2] );
    float rI                                     =  1.0f/r;
    float r2I                                    =  rI*rI;
    float rr3                                    = -rI*r2I;
    float rr5                                    = -3.0f*rr3*r2I;
    
    float dampProd                               = atomI.damp*atomJ.damp;
    float ratio                                  = (dampProd != 0.0f) ? (r/dampProd) : 1.0f;
    float pGamma                                 = atomI.thole > atomJ.thole ? atomJ.thole: atomI.thole;
    float damp                                   = ratio*ratio*ratio*pGamma;
    float dampExp                                = ( (dampProd != 0.0f) && (r < cAmoebaSim.scalingDistanceCutoff) ) ? expf( -damp ) : 0.0f; 

    rr3                                         *= (1.0f - dampExp);
    rr5                                         *= (1.0f - ( 1.0f + damp )*dampExp);
        
    float dDotDelta                              = rr5*(deltaR[0]*atomJ.inducedDipole[0]    + deltaR[1]*atomJ.inducedDipole[1]    + deltaR[2]*atomJ.inducedDipole[2] );
    fields[0][0]                                 = rr3*atomJ.inducedDipole[0] + dDotDelta*deltaR[0];
    fields[0][1]                                 = rr3*atomJ.inducedDipole[1] + dDotDelta*deltaR[1];
    fields[0][2]                                 = rr3*atomJ.inducedDipole[2] + dDotDelta*deltaR[2];
   
    dDotDelta                                    = rr5*(deltaR[0]*atomJ.inducedDipolePolar[0]    + deltaR[1]*atomJ.inducedDipolePolar[1]    + deltaR[2]*atomJ.inducedDipolePolar[2] );
    fields[1][0]                                 = rr3*atomJ.inducedDipolePolar[0] + dDotDelta*deltaR[0];
    fields[1][1]                                 = rr3*atomJ.inducedDipolePolar[1] + dDotDelta*deltaR[1];
    fields[1][2]                                 = rr3*atomJ.inducedDipolePolar[2] + dDotDelta*deltaR[2];
  
    dDotDelta                                    = rr5*(deltaR[0]*atomI.inducedDipole[0]    + deltaR[1]*atomI.inducedDipole[1]    + deltaR[2]*atomI.inducedDipole[2] );
    fields[2][0]                                 = rr3*atomI.inducedDipole[0] + dDotDelta*deltaR[0];
    fields[2][1]                                 = rr3*atomI.inducedDipole[1] + dDotDelta*deltaR[1];
    fields[2][2]                                 = rr3*atomI.inducedDipole[2] + dDotDelta*deltaR[2];
   
    dDotDelta                                    = rr5*(deltaR[0]*atomI.inducedDipolePolar[0]    + deltaR[1]*atomI.inducedDipolePolar[1]    + deltaR[2]*atomI.inducedDipolePolar[2] );
    fields[3][0]                                 = rr3*atomI.inducedDipolePolar[0] + dDotDelta*deltaR[0];
    fields[3][1]                                 = rr3*atomI.inducedDipolePolar[1] + dDotDelta*deltaR[1];
    fields[3][2]                                 = rr3*atomI.inducedDipolePolar[2] + dDotDelta*deltaR[2];

    dDotDelta                                    = rr5*(deltaR[0]*atomJ.inducedDipoleS[0]    + deltaR[1]*atomJ.inducedDipoleS[1]    + deltaR[2]*atomJ.inducedDipoleS[2] );
    fields[4][0]                                 = rr3*atomJ.inducedDipoleS[0] + dDotDelta*deltaR[0];
    fields[4][1]                                 = rr3*atomJ.inducedDipoleS[1] + dDotDelta*deltaR[1];
    fields[4][2]                                 = rr3*atomJ.inducedDipoleS[2] + dDotDelta*deltaR[2];
   
    dDotDelta                                    = rr5*(deltaR[0]*atomJ.inducedDipolePolarS[0]    + deltaR[1]*atomJ.inducedDipolePolarS[1]    + deltaR[2]*atomJ.inducedDipolePolarS[2] );
    fields[5][0]                                 = rr3*atomJ.inducedDipolePolarS[0] + dDotDelta*deltaR[0];
    fields[5][1]                                 = rr3*atomJ.inducedDipolePolarS[1] + dDotDelta*deltaR[1];
    fields[5][2]                                 = rr3*atomJ.inducedDipolePolarS[2] + dDotDelta*deltaR[2];
  
    dDotDelta                                    = rr5*(deltaR[0]*atomI.inducedDipoleS[0]    + deltaR[1]*atomI.inducedDipoleS[1]    + deltaR[2]*atomI.inducedDipoleS[2] );
    fields[6][0]                                 = rr3*atomI.inducedDipoleS[0] + dDotDelta*deltaR[0];
    fields[6][1]                                 = rr3*atomI.inducedDipoleS[1] + dDotDelta*deltaR[1];
    fields[6][2]                                 = rr3*atomI.inducedDipoleS[2] + dDotDelta*deltaR[2];
   
    dDotDelta                                    = rr5*(deltaR[0]*atomI.inducedDipolePolarS[0]    + deltaR[1]*atomI.inducedDipolePolarS[1]    + deltaR[2]*atomI.inducedDipolePolarS[2] );
    fields[7][0]                                 = rr3*atomI.inducedDipolePolarS[0] + dDotDelta*deltaR[0];
    fields[7][1]                                 = rr3*atomI.inducedDipolePolarS[1] + dDotDelta*deltaR[1];
    fields[7][2]                                 = rr3*atomI.inducedDipolePolarS[2] + dDotDelta*deltaR[2];


}

__device__ void calculateMutualInducedAndGkFieldsGkPairIxn_kernel( MutualInducedParticle& atomI, MutualInducedParticle& atomJ,
                                                                   float gkField[8][3]

#ifdef AMOEBA_DEBUG
               , float4* debugArray
#endif

 )
{

    float gux[5];
    float guy[5];
    float guz[5];
    float a[3][3];
    
    // ---------------------------------------------------------------------------------------
    
    float xr               = atomJ.x - atomI.x;
    float yr               = atomJ.y - atomI.y;
    float zr               = atomJ.z - atomI.z;

    float xr2              = xr*xr;
    float yr2              = yr*yr;
    float zr2              = zr*zr;

    float rb2              = atomI.bornRadius*atomJ.bornRadius;

    float r2               = xr2 + yr2 + zr2;
    float expterm          = expf(-r2/(cAmoebaSim.gkc*rb2));
    float expc             = expterm /cAmoebaSim.gkc; 
    //float dexpc            = -2.0f / (cAmoebaSim.gkc*rb2);

    float gf2              = 1.0f / (r2+rb2*expterm);
    float gf               = sqrtf(gf2);
    float gf3              = gf2 * gf;
    float gf5              = gf3 * gf2;

    float duixs            = atomI.inducedDipoleS[0];
    float duiys            = atomI.inducedDipoleS[1];
    float duizs            = atomI.inducedDipoleS[2];

    float puixs            = atomI.inducedDipolePolarS[0];
    float puiys            = atomI.inducedDipolePolarS[1];
    float puizs            = atomI.inducedDipolePolarS[2];
 
    float dukxs            = atomJ.inducedDipoleS[0];
    float dukys            = atomJ.inducedDipoleS[1];
    float dukzs            = atomJ.inducedDipoleS[2];

    float pukxs            = atomJ.inducedDipolePolarS[0];
    float pukys            = atomJ.inducedDipolePolarS[1];
    float pukzs            = atomJ.inducedDipolePolarS[2];
 
    // reaction potential auxiliary terms
 
    a[1][0]                = -gf3;
    a[2][0]                = 3.0f * gf5;

    // reaction potential gradient auxiliary terms

    float expc1            = 1.0f - expc;
    a[1][1]                = expc1 * a[2][0];
 
    // unweighted dipole reaction potential gradient tensor

    gux[2]                 = cAmoebaSim.fd * (a[1][0] + xr2*a[1][1]);
    gux[3]                 = cAmoebaSim.fd * xr*yr*a[1][1];
    gux[4]                 = cAmoebaSim.fd * xr*zr*a[1][1];

    guy[2]                 = gux[3];
    guy[3]                 = cAmoebaSim.fd * (a[1][0] + yr2*a[1][1]);
    guy[4]                 = cAmoebaSim.fd * yr*zr*a[1][1];

    guz[2]                 = gux[4];
    guz[3]                 = guy[4];
    guz[4]                 = cAmoebaSim.fd * (a[1][0] + zr2*a[1][1]);
 
    gkField[0][0]          = dukxs*gux[2]+dukys*guy[2]+dukzs*guz[2];
    gkField[0][1]          = dukxs*gux[3]+dukys*guy[3]+dukzs*guz[3];
    gkField[0][2]          = dukxs*gux[4]+dukys*guy[4]+dukzs*guz[4];

    gkField[1][0]          = duixs*gux[2]+duiys*guy[2]+duizs*guz[2];
    gkField[1][1]          = duixs*gux[3]+duiys*guy[3]+duizs*guz[3];
    gkField[1][2]          = duixs*gux[4]+duiys*guy[4]+duizs*guz[4];

    gkField[2][0]          = pukxs*gux[2]+pukys*guy[2]+pukzs*guz[2];
    gkField[2][1]          = pukxs*gux[3]+pukys*guy[3]+pukzs*guz[3];
    gkField[2][2]          = pukxs*gux[4]+pukys*guy[4]+pukzs*guz[4];

    gkField[3][0]          = puixs*gux[2]+puiys*guy[2]+puizs*guz[2];
    gkField[3][1]          = puixs*gux[3]+puiys*guy[3]+puizs*guz[3];
    gkField[3][2]          = puixs*gux[4]+puiys*guy[4]+puizs*guz[4];

}

#ifdef AMOEBA_DEBUG
__device__ static int debugAccumulate( int index, float4* debugArray, float* field, unsigned int addMask, float idLabel )
{
    index                             += cAmoebaSim.paddedNumberOfAtoms;
    debugArray[index].x                = addMask ? field[0] : 0.0f;
    debugArray[index].y                = addMask ? field[1] : 0.0f;
    debugArray[index].z                = addMask ? field[2] : 0.0f;
    debugArray[index].w                = idLabel;

    return index;
}
#endif


// Include versions of the kernels for N^2 calculations.

#define METHOD_NAME(a, b) a##N2##b
#include "kCalculateAmoebaCudaMutualInducedAndGkFields.h"
#define USE_OUTPUT_BUFFER_PER_WARP
#undef METHOD_NAME
#define METHOD_NAME(a, b) a##N2ByWarp##b
#include "kCalculateAmoebaCudaMutualInducedAndGkFields.h"

__global__ 
#if (__CUDA_ARCH__ >= 200)
__launch_bounds__(GF1XX_THREADS_PER_BLOCK, 1)
#elif (__CUDA_ARCH__ >= 130)
__launch_bounds__(GT2XX_THREADS_PER_BLOCK, 1)
#else
__launch_bounds__(G8X_THREADS_PER_BLOCK, 1)
#endif
void kInitializeMutualInducedAndGkField_kernel(
                   float* fixedEField,
                   float* fixedEFieldPolar,
                   float* fixedGkField,
                   float* polarizability,
                   float* inducedDipole,
                   float* inducedDipolePolar,
                   float* inducedDipoleS,
                   float* inducedDipolePolarS )
{

    int threadId = __mul24(blockIdx.x,blockDim.x) + threadIdx.x;
    if( threadId >= 3*cAmoebaSim.numberOfAtoms )return;

    fixedEField[threadId]          *= polarizability[threadId];
    inducedDipole[threadId]         = fixedEField[threadId];

    fixedEFieldPolar[threadId]     *= polarizability[threadId];
    inducedDipolePolar[threadId]    = fixedEFieldPolar[threadId];

    fixedGkField[threadId]         *= polarizability[threadId];
    inducedDipoleS[threadId]        = fixedEField[threadId]       + fixedGkField[threadId];
    inducedDipolePolarS[threadId]   = fixedEFieldPolar[threadId]  + fixedGkField[threadId];

}

__global__ 
#if (__CUDA_ARCH__ >= 200)
__launch_bounds__(GF1XX_THREADS_PER_BLOCK, 1)
#elif (__CUDA_ARCH__ >= 130)
__launch_bounds__(GT2XX_THREADS_PER_BLOCK, 1)
#else
__launch_bounds__(G8X_THREADS_PER_BLOCK, 1)
#endif
void kReduceMutualInducedAndGkFieldDelta_kernel( float* arrayOfDeltas1, float* arrayOfDeltas2,
                                                 float* arrayOfDeltas3, float* arrayOfDeltas4, float* epsilon )
{
    extern __shared__ float4 delta[];

    delta[threadIdx.x].x    = 0.0f;
    delta[threadIdx.x].y    = 0.0f;
    delta[threadIdx.x].z    = 0.0f;
    delta[threadIdx.x].w    = 0.0f;

    unsigned int pos        = threadIdx.x;

    // load deltas

    while( pos <  3*cAmoebaSim.numberOfAtoms )
    {   
        delta[threadIdx.x].x  += arrayOfDeltas1[pos];
        delta[threadIdx.x].y  += arrayOfDeltas2[pos];
        delta[threadIdx.x].z  += arrayOfDeltas3[pos];
        delta[threadIdx.x].w  += arrayOfDeltas4[pos];
        pos                   += blockDim.x*gridDim.x;
    }   
    __syncthreads();

    // sum the deltas

    for (int offset = 1; offset < blockDim.x; offset *= 2 )
    {   
        if (threadIdx.x + offset < blockDim.x && (threadIdx.x & (2*offset-1)) == 0)
        {
            delta[threadIdx.x].x   += delta[threadIdx.x+offset].x;
            delta[threadIdx.x].y   += delta[threadIdx.x+offset].y;
            delta[threadIdx.x].z   += delta[threadIdx.x+offset].z;
            delta[threadIdx.x].w   += delta[threadIdx.x+offset].w;
        }
        __syncthreads();
    }   

    // set epsilons

    if (threadIdx.x == 0)
    {   
        epsilon[0]  = delta[0].x;
        epsilon[0]  = epsilon[0] < delta[0].y ? delta[0].y : epsilon[0];
        epsilon[0]  = epsilon[0] < delta[0].z ? delta[0].z : epsilon[0];
        epsilon[0]  = epsilon[0] < delta[0].w ? delta[0].w : epsilon[0];
        epsilon[0]  = 4.8033324f*sqrtf( epsilon[0]/( (float) cAmoebaSim.numberOfAtoms ) );
#ifdef AMOEBA_DEBUG
        epsilon[1]  = 4.8033324f*sqrtf( delta[0].x/( (float) cAmoebaSim.numberOfAtoms ) );
        epsilon[2]  = 4.8033324f*sqrtf( delta[0].y/( (float) cAmoebaSim.numberOfAtoms ) );
        epsilon[3]  = 4.8033324f*sqrtf( delta[0].z/( (float) cAmoebaSim.numberOfAtoms ) );
        epsilon[4]  = 4.8033324f*sqrtf( delta[0].w/( (float) cAmoebaSim.numberOfAtoms ) );
#endif
    }   
}

/**

   matrixProduct/matrixProductP contains epsilon**2 on output

*/
__global__ 
#if (__CUDA_ARCH__ >= 200)
__launch_bounds__(GF1XX_THREADS_PER_BLOCK, 1)
#elif (__CUDA_ARCH__ >= 130)
__launch_bounds__(GT2XX_THREADS_PER_BLOCK, 1)
#else
__launch_bounds__(G8X_THREADS_PER_BLOCK, 1)
#endif
void kSorUpdateMutualInducedAndGkField_kernel(
                   float* polarizability,
                   float* inducedDipole, float* inducedDipoleP,
                   float* fixedEField,   float* fixedEFieldP,
                   float* matrixProduct, float* matrixProductP )
{

    float polarSOR = 0.70f;
    int threadId   = __mul24(blockIdx.x,blockDim.x) + threadIdx.x;
    if( threadId  >= 3*cAmoebaSim.numberOfAtoms)return;

    float previousDipole                = inducedDipole[threadId];
    float previousDipoleP               = inducedDipoleP[threadId];

    inducedDipole[threadId]             = fixedEField[threadId]     + polarizability[threadId]*matrixProduct[threadId];
    inducedDipoleP[threadId]            = fixedEFieldP[threadId]    + polarizability[threadId]*matrixProductP[threadId];

    inducedDipole[threadId]             = previousDipole   + polarSOR*( inducedDipole[threadId]   - previousDipole  );   
    inducedDipoleP[threadId]            = previousDipoleP  + polarSOR*( inducedDipoleP[threadId]  - previousDipoleP );

    matrixProduct[threadId]             = ( inducedDipole[threadId]  - previousDipole  )*( inducedDipole[threadId]  - previousDipole  );
    matrixProductP[threadId]            = ( inducedDipoleP[threadId] - previousDipoleP )*( inducedDipoleP[threadId] - previousDipoleP );

}

__global__ 
#if (__CUDA_ARCH__ >= 200)
__launch_bounds__(GF1XX_THREADS_PER_BLOCK, 1)
#elif (__CUDA_ARCH__ >= 130)
__launch_bounds__(GT2XX_THREADS_PER_BLOCK, 1)
#else
__launch_bounds__(G8X_THREADS_PER_BLOCK, 1)
#endif
void kSorUpdateMutualInducedAndGkFieldS_kernel(
                   float* polarizability,
                   float* inducedDipole, float* inducedDipoleP,
                   float* fixedEField,   float* fixedEFieldP,
                   float* fixedGkField,
                   float* matrixProduct, float* matrixProductP )
{

    float polarSOR = 0.70f;
    int threadId   = __mul24(blockIdx.x,blockDim.x) + threadIdx.x;
    if( threadId  >= 3*cAmoebaSim.numberOfAtoms)return;

    float previousDipole                = inducedDipole[threadId];
    float previousDipoleP               = inducedDipoleP[threadId];

    inducedDipole[threadId]             = fixedGkField[threadId]    + fixedEField[threadId]     + polarizability[threadId]*matrixProduct[threadId];
    inducedDipoleP[threadId]            = fixedGkField[threadId]    + fixedEFieldP[threadId]    + polarizability[threadId]*matrixProductP[threadId];

    inducedDipole[threadId]             = previousDipole   + polarSOR*( inducedDipole[threadId]   - previousDipole  );   
    inducedDipoleP[threadId]            = previousDipoleP  + polarSOR*( inducedDipoleP[threadId]  - previousDipoleP );

    matrixProduct[threadId]             = ( inducedDipole[threadId]  - previousDipole  )*( inducedDipole[threadId]  - previousDipole  );
    matrixProductP[threadId]            = ( inducedDipoleP[threadId] - previousDipoleP )*( inducedDipoleP[threadId] - previousDipoleP );

}

// reduce psWorkArray_3_1 -> outputArray
// reduce psWorkArray_3_2 -> outputPolarArray
// reduce psWorkArray_3_3 -> outputArrayS
// reduce psWorkArray_3_4 -> outputPolarArrayS

static void kReduceMutualInducedAndGkFields(amoebaGpuContext amoebaGpu,
                                            CUDAStream<float>* outputArray,  CUDAStream<float>* outputPolarArray,
                                            CUDAStream<float>* outputArrayS, CUDAStream<float>* outputPolarArrayS )
{
    kReduceFields_kernel<<<amoebaGpu->nonbondBlocks, amoebaGpu->fieldReduceThreadsPerBlock>>>(
                               amoebaGpu->paddedNumberOfAtoms*3, amoebaGpu->outputBuffers,
                               amoebaGpu->psWorkArray_3_1->_pDevStream[0], outputArray->_pDevStream[0] );
    LAUNCHERROR("kReduceMutualInducedAndGkFields1");

    kReduceFields_kernel<<<amoebaGpu->nonbondBlocks, amoebaGpu->fieldReduceThreadsPerBlock>>>(
                               amoebaGpu->paddedNumberOfAtoms*3, amoebaGpu->outputBuffers,
                               amoebaGpu->psWorkArray_3_2->_pDevStream[0], outputPolarArray->_pDevStream[0] );
    LAUNCHERROR("kReduceMutualInducedAndGkFields2");

    kReduceFields_kernel<<<amoebaGpu->nonbondBlocks, amoebaGpu->fieldReduceThreadsPerBlock>>>(
                               amoebaGpu->paddedNumberOfAtoms*3, amoebaGpu->outputBuffers,
                               amoebaGpu->psWorkArray_3_3->_pDevStream[0], outputArrayS->_pDevStream[0] );
    LAUNCHERROR("kReduceMutualInducedAndGkFields3");

    kReduceFields_kernel<<<amoebaGpu->nonbondBlocks, amoebaGpu->fieldReduceThreadsPerBlock>>>(
                               amoebaGpu->paddedNumberOfAtoms*3, amoebaGpu->outputBuffers,
                               amoebaGpu->psWorkArray_3_4->_pDevStream[0], outputPolarArrayS->_pDevStream[0] );
    LAUNCHERROR("kReduceMutualInducedAndGkFields4");
}

#ifdef AMOEBA_DEBUG
#if 0
static void printMiFieldBuffer( amoebaGpuContext amoebaGpu, unsigned int bufferIndex )
{
    (void) fprintf( amoebaGpu->log, "MI Field Buffer %u\n", bufferIndex );
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
}

static void printMiFieldAtomBuffers( amoebaGpuContext amoebaGpu, unsigned int targetAtom )
{
    (void) fprintf( amoebaGpu->log, "MI Field atom %u\n", targetAtom );
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
#endif

/**---------------------------------------------------------------------------------------

   Compute mutual induce field

   @param amoebaGpu        amoebaGpu context

   --------------------------------------------------------------------------------------- */

static void cudaComputeAmoebaMutualInducedAndGkFieldMatrixMultiply( amoebaGpuContext amoebaGpu,
                                                                    CUDAStream<float>* outputArray,  CUDAStream<float>* outputPolarArray,
                                                                    CUDAStream<float>* outputArrayS, CUDAStream<float>* outputPolarArrayS )
{
  
   // ---------------------------------------------------------------------------------------

    static unsigned int threadsPerBlock = 0;

   // ---------------------------------------------------------------------------------------

  gpuContext gpu    = amoebaGpu->gpuContext;

#ifdef AMOEBA_DEBUG
    static int iteration = 1;
    int targetAtom    = 0;
    static const char* methodName       = "cudaComputeAmoebaMutualInducedAndGkFieldMatrixMultiply";
    if( 1 && amoebaGpu->log ){
        (void) fprintf( amoebaGpu->log, "%s: scalingDistanceCutoff=%.5f\n",
                        methodName, amoebaGpu->scalingDistanceCutoff );
        (void) fflush( amoebaGpu->log );
    }
    int paddedNumberOfAtoms                    = amoebaGpu->gpuContext->sim.paddedNumberOfAtoms;
    CUDAStream<float4>* debugArray             = new CUDAStream<float4>(paddedNumberOfAtoms*paddedNumberOfAtoms, 1, "DebugArray");
    memset( debugArray->_pSysStream[0],      0, sizeof( float )*4*paddedNumberOfAtoms*paddedNumberOfAtoms);
    debugArray->Upload();
#endif

    // clear output arrays

    kClearFields_3( amoebaGpu, 4 );

    // set threads/block first time through

    if( threadsPerBlock == 0 ){
        unsigned int maxThreads;
        if (gpu->sm_version >= SM_20)
            maxThreads = 384;
        else if (gpu->sm_version >= SM_12)
            maxThreads = 128;
        else
            maxThreads = 64;
        threadsPerBlock = std::min(getThreadsPerBlock( amoebaGpu, sizeof(MutualInducedParticle)), maxThreads);
    }
    
    if (gpu->bOutputBufferPerWarp){
        kCalculateAmoebaMutualInducedAndGkFieldsN2ByWarp_kernel<<<amoebaGpu->nonbondBlocks, threadsPerBlock, sizeof(MutualInducedParticle)*threadsPerBlock>>>(
                                                                 amoebaGpu->psWorkUnit->_pDevStream[0],
                                                                 amoebaGpu->psWorkArray_3_1->_pDevStream[0],
                                                                 amoebaGpu->psWorkArray_3_2->_pDevStream[0],
                                                                 amoebaGpu->psWorkArray_3_3->_pDevStream[0],
#ifdef AMOEBA_DEBUG
                                                                 amoebaGpu->psWorkArray_3_4->_pDevStream[0],
                                                                 debugArray->_pDevStream[0], targetAtom );
#else
                                                                 amoebaGpu->psWorkArray_3_4->_pDevStream[0] );
#endif
    } else {

#ifdef AMOEBA_DEBUG
        (void) fprintf( amoebaGpu->log, "N2 no warp\n" );
        (void) fprintf( amoebaGpu->log, "cudaComputeAmoebaMutualInducedAndGkFieldMatrixMultiply numBlocks=%u numThreads=%u bufferPerWarp=%u atm=%u shrd=%u Ebuf=%u ixnCt=%u workUnits=%u\n",
                        amoebaGpu->nonbondBlocks, threadsPerBlock, amoebaGpu->bOutputBufferPerWarp,
                        sizeof(MutualInducedParticle), sizeof(MutualInducedParticle)*threadsPerBlock,
                        amoebaGpu->energyOutputBuffers, (*gpu->psInteractionCount)[0], gpu->sim.workUnits );
        (void) fflush( amoebaGpu->log );
#endif

        kCalculateAmoebaMutualInducedAndGkFieldsN2_kernel<<<amoebaGpu->nonbondBlocks, threadsPerBlock, sizeof(MutualInducedParticle)*threadsPerBlock>>>(
                                                            amoebaGpu->psWorkUnit->_pDevStream[0],
                                                            amoebaGpu->psWorkArray_3_1->_pDevStream[0],
                                                            amoebaGpu->psWorkArray_3_2->_pDevStream[0],
                                                            amoebaGpu->psWorkArray_3_3->_pDevStream[0],
#ifdef AMOEBA_DEBUG
                                                            amoebaGpu->psWorkArray_3_4->_pDevStream[0],
                                                            debugArray->_pDevStream[0], targetAtom );
#else
                                                            amoebaGpu->psWorkArray_3_4->_pDevStream[0] );
#endif

    }
    LAUNCHERROR("kCalculateAmoebaMutualInducedAndGkFields");  

    kReduceMutualInducedAndGkFields( amoebaGpu, outputArray, outputPolarArray, outputArrayS, outputPolarArrayS );

#ifdef AMOEBA_DEBUG
        amoebaGpu->psWorkArray_3_1->Download();
        amoebaGpu->psWorkArray_3_2->Download();
        amoebaGpu->psWorkArray_3_3->Download();
        amoebaGpu->psWorkArray_3_4->Download();

        //printMiFieldAtomBuffers( amoebaGpu, (targetAtom + 0) );
        //printMiFieldAtomBuffers( amoebaGpu, (targetAtom + 1) );
        //printMiFieldAtomBuffers( amoebaGpu, 100 );
        //printMiFieldBuffer( amoebaGpu, 0 );
        //printMiFieldBuffer( amoebaGpu, 1 );
        //printMiFieldBuffer( amoebaGpu, 37 );
        //printMiFieldBuffer( amoebaGpu, 38 );

    if( amoebaGpu->log && iteration == -1 ){

        (void) fprintf( amoebaGpu->log, "Finished MI kernel execution %d\n", iteration ); (void) fflush( amoebaGpu->log );

        outputArray->Download();
        outputPolarArray->Download();

        outputArrayS->Download();
        outputPolarArrayS->Download();

        debugArray->Download();

        int maxPrint        = 33;
        for( int ii = 0; ii < gpu->natoms; ii++ ){
           (void) fprintf( amoebaGpu->log, "%5d ", ii); 

            int indexOffset     = ii*3;
    
           // Mi

           (void) fprintf( amoebaGpu->log," Mult[%16.9e %16.9e %16.9e] ",
                           outputArray->_pSysStream[0][indexOffset],
                           outputArray->_pSysStream[0][indexOffset+1],
                           outputArray->_pSysStream[0][indexOffset+2] );
    
           // Mi polar

           (void) fprintf( amoebaGpu->log," MultP[%16.9e %16.9e %16.9e]\n",
                           outputPolarArray->_pSysStream[0][indexOffset],
                           outputPolarArray->_pSysStream[0][indexOffset+1],
                           outputPolarArray->_pSysStream[0][indexOffset+2] );

           // MiS

           (void) fprintf( amoebaGpu->log,"      MultS[%16.9e %16.9e %16.9e] ",
                           outputArrayS->_pSysStream[0][indexOffset],
                           outputArrayS->_pSysStream[0][indexOffset+1],
                           outputArrayS->_pSysStream[0][indexOffset+2] );
    
           // Mi polarS

           (void) fprintf( amoebaGpu->log,"MultPS[%16.9e %16.9e %16.9e]\n",
                           outputPolarArrayS->_pSysStream[0][indexOffset],
                           outputPolarArrayS->_pSysStream[0][indexOffset+1],
                           outputPolarArrayS->_pSysStream[0][indexOffset+2] );

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
                   (void) fprintf( amoebaGpu->log,"%4d %4d Rint [%16.9e %16.9e %16.9e %16.9e]\n",
                                   ii, jj,
                                   debugArray->_pSysStream[0][debugIndex].x, debugArray->_pSysStream[0][debugIndex].y,
                                   debugArray->_pSysStream[0][debugIndex].z, debugArray->_pSysStream[0][debugIndex].w );

                   for( int kk = 0; kk < 9; kk++ ){
                       debugIndex += paddedNumberOfAtoms;
                       (void) fprintf( amoebaGpu->log,"[%16.9e %16.9e %16.9e %5.1f]\n",
                                       debugArray->_pSysStream[0][debugIndex].x, debugArray->_pSysStream[0][debugIndex].y, debugArray->_pSysStream[0][debugIndex].z, debugArray->_pSysStream[0][debugIndex].w );
                   }
               }
               (void) fprintf( amoebaGpu->log,"\n" );
           }
           (void) fprintf( amoebaGpu->log,"\n" );
           if( ii == maxPrint && (gpu->natoms - maxPrint) > ii ){
                ii = gpu->natoms - maxPrint;
           }
        }
        (void) fflush( amoebaGpu->log );
        iteration++;

     }
     delete debugArray;
#endif

}

/**---------------------------------------------------------------------------------------

   Compute mutual induce field

   @param amoebaGpu        amoebaGpu context

   --------------------------------------------------------------------------------------- */

static void cudaComputeAmoebaMutualInducedAndGkFieldBySOR( amoebaGpuContext amoebaGpu )
{
  
   // ---------------------------------------------------------------------------------------

    static int timestep                  = 0;
    timestep++;
#ifdef AMOEBA_DEBUG
    static const char* methodName        = "cudaComputeAmoebaMutualInducedAndGkFieldBySOR";
    std::vector<int> fileId;
    fileId.resize( 2 );
    fileId[0] = timestep;
    fileId[1] = 1;
#endif

   // ---------------------------------------------------------------------------------------

    int done;
    int iteration;

    gpuContext gpu    = amoebaGpu->gpuContext;
    int numOfElems     = gpu->natoms*3;
    int numThreads     = min( THREADS_PER_BLOCK, numOfElems );
    int numBlocks      = numOfElems/numThreads;

    if( (numOfElems % numThreads) != 0 )numBlocks++;

#ifdef AMOEBA_DEBUG
    if( amoebaGpu->log && timestep == 1 ){
        (void) fprintf( amoebaGpu->log, "%s %d numOfElems=%d numThreads=%d numBlocks=%d "
                        "maxIterations=%d targetEpsilon=%.3e\n", 
                        methodName, gpu->natoms, numOfElems, numThreads, numBlocks,
                        amoebaGpu->mutualInducedMaxIterations, amoebaGpu->mutualInducedTargetEpsilon);
        (void) fflush( amoebaGpu->log );
    }   
#endif

   // ---------------------------------------------------------------------------------------

    // set  E_Field & E_FieldPolar] to [ E_Field & E_FieldPolar]*Polarizability
    // initialize [ InducedDipole & InducedDipolePolar ] to [ E_Field & E_FieldPolar]*Polarizability

    kInitializeMutualInducedAndGkField_kernel<<< numBlocks, numThreads >>>(
         amoebaGpu->psE_Field->_pDevStream[0],
         amoebaGpu->psE_FieldPolar->_pDevStream[0],
         amoebaGpu->psGk_Field->_pDevStream[0],
         amoebaGpu->psPolarizability->_pDevStream[0],
         amoebaGpu->psInducedDipole->_pDevStream[0],
         amoebaGpu->psInducedDipolePolar->_pDevStream[0],
         amoebaGpu->psInducedDipoleS->_pDevStream[0],
         amoebaGpu->psInducedDipolePolarS->_pDevStream[0] );
    LAUNCHERROR("kInitializeMutualInducedAndGkField");  

#ifdef AMOEBA_DEBUG
    if( amoebaGpu->log ){

        amoebaGpu->psE_Field->Download();
        amoebaGpu->psE_FieldPolar->Download();
        amoebaGpu->psInducedDipole->Download(),
        amoebaGpu->psInducedDipolePolar->Download();
        amoebaGpu->psInducedDipoleS->Download(),
        amoebaGpu->psInducedDipolePolarS->Download();
        amoebaGpu->psPolarizability->Download();

        (void) fprintf( amoebaGpu->log, "%s Initial setup for matrix multiply\n", methodName );
        int offset   = 0;
        int maxPrint = 10;
        for( int ii = 0; ii < gpu->natoms; ii++ ){
            (void) fprintf( amoebaGpu->log, "\n%4d pol=%12.4e\n", ii, 
                            amoebaGpu->psPolarizability->_pSysStream[0][offset] );
            if( amoebaGpu->psPolarizability->_pSysStream[0][offset] != amoebaGpu->psPolarizability->_pSysStream[0][offset+1] ||
                amoebaGpu->psPolarizability->_pSysStream[0][offset] != amoebaGpu->psPolarizability->_pSysStream[0][offset+2] ){
                (void) fprintf( amoebaGpu->log, "PolX!!! %12.4e %12.4e ", amoebaGpu->psPolarizability->_pSysStream[0][offset+1], amoebaGpu->psPolarizability->_pSysStream[0][offset+2] ); 
            }

            (void) fprintf( amoebaGpu->log," E[%14.6e %14.6e %14.6e]  Mi[%14.6e %14.6e %14.6e]\n",
                            amoebaGpu->psE_Field->_pSysStream[0][offset],       amoebaGpu->psE_Field->_pSysStream[0][offset+1],       amoebaGpu->psE_Field->_pSysStream[0][offset+2],
                            amoebaGpu->psInducedDipole->_pSysStream[0][offset], amoebaGpu->psInducedDipole->_pSysStream[0][offset+1], amoebaGpu->psInducedDipole->_pSysStream[0][offset+2] );
            (void) fprintf( amoebaGpu->log,"Ep[%14.6e %14.6e %14.6e] Mip[%14.6e %14.6e %14.6e]\n",
                            amoebaGpu->psE_FieldPolar->_pSysStream[0][offset],       amoebaGpu->psE_FieldPolar->_pSysStream[0][offset+1],       amoebaGpu->psE_FieldPolar->_pSysStream[0][offset+2],
                            amoebaGpu->psInducedDipolePolar->_pSysStream[0][offset], amoebaGpu->psInducedDipolePolar->_pSysStream[0][offset+1], amoebaGpu->psInducedDipolePolar->_pSysStream[0][offset+2] );

            (void) fprintf( amoebaGpu->log,"Gk[%14.6e %14.6e %14.6e] MiS[%14.6e %14.6e %14.6e] MipS[%14.6e %14.6e %14.6e]\n",
                            amoebaGpu->psGk_Field->_pSysStream[0][offset],            amoebaGpu->psGk_Field->_pSysStream[0][offset+1],            amoebaGpu->psGk_Field->_pSysStream[0][offset+2],
                            amoebaGpu->psInducedDipoleS->_pSysStream[0][offset],      amoebaGpu->psInducedDipoleS->_pSysStream[0][offset+1],      amoebaGpu->psInducedDipoleS->_pSysStream[0][offset+2],
                            amoebaGpu->psInducedDipolePolarS->_pSysStream[0][offset], amoebaGpu->psInducedDipolePolarS->_pSysStream[0][offset+1], amoebaGpu->psInducedDipolePolarS->_pSysStream[0][offset+2] );
            offset += 3;
            if( ii == maxPrint && (ii < (gpu->natoms - maxPrint) ) )ii =  (gpu->natoms - maxPrint);
        }   
        (void) fflush( amoebaGpu->log );
    }   
#endif

    // ---------------------------------------------------------------------------------------
 
    done      = 0;
    iteration = 1;

    while( !done ){

        // matrix multiply

        cudaComputeAmoebaMutualInducedAndGkFieldMatrixMultiply( amoebaGpu,
                                                                amoebaGpu->psWorkVector[0],  amoebaGpu->psWorkVector[1],
                                                                amoebaGpu->psWorkVector[2],  amoebaGpu->psWorkVector[3] );

        LAUNCHERROR("cudaComputeAmoebaMutualInducedAndGkFieldMatrixMultiply");  

    // ---------------------------------------------------------------------------------------

        // post matrix multiply

        kSorUpdateMutualInducedAndGkField_kernel<<< numBlocks, numThreads >>>(
           amoebaGpu->psPolarizability->_pDevStream[0],
           amoebaGpu->psInducedDipole->_pDevStream[0],     amoebaGpu->psInducedDipolePolar->_pDevStream[0],
           amoebaGpu->psE_Field->_pDevStream[0],           amoebaGpu->psE_FieldPolar->_pDevStream[0],
           amoebaGpu->psWorkVector[0]->_pDevStream[0],     amoebaGpu->psWorkVector[1]->_pDevStream[0] );
        LAUNCHERROR("cudaComputeAmoebaMutualInducedAndGkFieldSorUpdate1");  

        kSorUpdateMutualInducedAndGkFieldS_kernel<<< numBlocks, numThreads >>>(
           amoebaGpu->psPolarizability->_pDevStream[0],
           amoebaGpu->psInducedDipoleS->_pDevStream[0],    amoebaGpu->psInducedDipolePolarS->_pDevStream[0],
           amoebaGpu->psE_Field->_pDevStream[0],          amoebaGpu->psE_FieldPolar->_pDevStream[0],
           amoebaGpu->psGk_Field->_pDevStream[0],
           amoebaGpu->psWorkVector[2]->_pDevStream[0],     amoebaGpu->psWorkVector[3]->_pDevStream[0] );
        LAUNCHERROR("cudaComputeAmoebaMutualInducedAndGkFieldSorUpdate2");  

        // get total epsilon -- performing sums on gpu

        kReduceMutualInducedAndGkFieldDelta_kernel<<<1, amoebaGpu->epsilonThreadsPerBlock, 4*sizeof(float)*amoebaGpu->epsilonThreadsPerBlock>>>(
           amoebaGpu->psWorkVector[0]->_pDevStream[0], amoebaGpu->psWorkVector[1]->_pDevStream[0],
           amoebaGpu->psWorkVector[2]->_pDevStream[0], amoebaGpu->psWorkVector[3]->_pDevStream[0],
           amoebaGpu->psCurrentEpsilon->_pDevStream[0] );
        LAUNCHERROR("kReduceMutualInducedAndGkFieldDelta_kernel");

#if 0
        // get total epsilon -- performing sums on cpu
{
        float sum[4];
        float currentEpsilon = -1.0e30;
        for( int ii = 0; ii < 4; ii++ ){
            sum[ii]                  = cudaGetSum( 3*gpu->natoms, amoebaGpu->psWorkVector[ii]);
            sum[ii]                  = 4.8033324f*sqrtf( sum[ii]/( (float) gpu->natoms) );
            if( sum[ii] > currentEpsilon ){
                currentEpsilon = sum[ii];
            }
        }

        amoebaGpu->mutualInducedCurrentEpsilon  = currentEpsilon;
        (void) fprintf( amoebaGpu->log, "%s iteration=%3d eps %14.6e [%14.6e %14.6e] done=%d sums=%14.6e %14.6e %14.6e %14.6e\n",
                        methodName, iteration, amoebaGpu->mutualInducedCurrentEpsilon,
                           amoebaGpu->psCurrentEpsilon->_pSysStream[0][1], 
                           amoebaGpu->psCurrentEpsilon->_pSysStream[0][2], done, sum[0], sum[1], sum[2], sum[3] );
}
#endif

        // Debye=4.8033324f

        amoebaGpu->psCurrentEpsilon->Download();
        float currentEpsilon                     = amoebaGpu->psCurrentEpsilon->_pSysStream[0][0];
        amoebaGpu->mutualInducedCurrentEpsilon   = currentEpsilon;

        if( iteration > amoebaGpu->mutualInducedMaxIterations || amoebaGpu->mutualInducedCurrentEpsilon < amoebaGpu->mutualInducedTargetEpsilon ){ 
            done = 1;
        }

#ifdef AMOEBA_DEBUG
        if( amoebaGpu->log ){
           amoebaGpu->psInducedDipole->Download();
           amoebaGpu->psInducedDipolePolar->Download();
           amoebaGpu->psInducedDipoleS->Download();
           amoebaGpu->psInducedDipolePolarS->Download();
           amoebaGpu->psWorkVector[2]->Download();
           amoebaGpu->psWorkVector[3]->Download();
           amoebaGpu->psGk_Field->Download();

#if 1
           (void) fprintf( amoebaGpu->log, "%s iteration=%3d eps %14.6e [%14.6e %14.6e] done=%d\n",
                           methodName, iteration, amoebaGpu->mutualInducedCurrentEpsilon,
                           amoebaGpu->psCurrentEpsilon->_pSysStream[0][1], 
                           amoebaGpu->psCurrentEpsilon->_pSysStream[0][2], done );
#else
           (void) fprintf( amoebaGpu->log, "%s iteration=%3d eps %14.6e %14.6e crrntEps=%14.6e %14.6e %14.6e %14.6e done=%d\n",
                           methodName, iteration, sum1, sum2, amoebaGpu->mutualInducedCurrentEpsilon,
                           amoebaGpu->psCurrentEpsilon->_pSysStream[0][0], 
                           amoebaGpu->psCurrentEpsilon->_pSysStream[0][1], 
                           amoebaGpu->psCurrentEpsilon->_pSysStream[0][2], done );
#endif
           (void) fflush( amoebaGpu->log );

            int offset   = 0;
            int maxPrint = 20;
            for( int ii = 0; ii < gpu->natoms; ii++ ){
                (void) fprintf( amoebaGpu->log, "%4d ", ii ); 
    
                (void) fprintf( amoebaGpu->log," Mi[%14.6e %14.6e %14.6e] ",
                                amoebaGpu->psInducedDipole->_pSysStream[0][offset], amoebaGpu->psInducedDipole->_pSysStream[0][offset+1], amoebaGpu->psInducedDipole->_pSysStream[0][offset+2] );
                (void) fprintf( amoebaGpu->log,"Mip[%14.6e %14.6e %14.6e]",
                                amoebaGpu->psInducedDipolePolar->_pSysStream[0][offset], amoebaGpu->psInducedDipolePolar->_pSysStream[0][offset+1], amoebaGpu->psInducedDipolePolar->_pSysStream[0][offset+2] );
                (void) fprintf( amoebaGpu->log," MiS[%14.6e %14.6e %14.6e] ",
                                amoebaGpu->psInducedDipoleS->_pSysStream[0][offset], amoebaGpu->psInducedDipoleS->_pSysStream[0][offset+1], amoebaGpu->psInducedDipoleS->_pSysStream[0][offset+2] );
                (void) fprintf( amoebaGpu->log,"MipS[%14.6e %14.6e %14.6e]\n",
                                amoebaGpu->psInducedDipolePolarS->_pSysStream[0][offset], amoebaGpu->psInducedDipolePolarS->_pSysStream[0][offset+1], amoebaGpu->psInducedDipolePolarS->_pSysStream[0][offset+2] );
/*
                (void) fprintf( amoebaGpu->log,"Gk [%14.6e %14.6e %14.6e]\n",
                                amoebaGpu->psGk_Field->_pSysStream[0][offset], amoebaGpu->psGk_Field->_pSysStream[0][offset+1], amoebaGpu->psGk_Field->_pSysStream[0][offset+2] );
                (void) fprintf( amoebaGpu->log,"W2 [%14.6e %14.6e %14.6e] ",
                                amoebaGpu->psWorkVector[2]->_pSysStream[0][offset], amoebaGpu->psWorkVector[2]->_pSysStream[0][offset+1], amoebaGpu->psWorkVector[2]->_pSysStream[0][offset+2] );
                (void) fprintf( amoebaGpu->log,"W3 [%14.6e %14.6e %14.6e]\n",
                                amoebaGpu->psWorkVector[3]->_pSysStream[0][offset], amoebaGpu->psWorkVector[3]->_pSysStream[0][offset+1], amoebaGpu->psWorkVector[3]->_pSysStream[0][offset+2] );
*/
                if( ii == maxPrint && (ii < (gpu->natoms - maxPrint) ) ){
                    ii =  (gpu->natoms - maxPrint);
                    offset = 3*(ii+1);
                } else {
                    offset += 3;
                }
            }   
            (void) fflush( amoebaGpu->log );
        }
#endif
        iteration++;
    }

    amoebaGpu->mutualInducedDone             = done;
    amoebaGpu->mutualInducedConverged        = ( !done || iteration > amoebaGpu->mutualInducedMaxIterations ) ? 0 : 1;

    if( amoebaGpu->log ){
        trackMutualInducedIterations( amoebaGpu, iteration );
    }

#ifdef AMOEBA_DEBUG
    if( 0 ){
        std::vector<int> fileId;
        //fileId.push_back( 0 );
        VectorOfDoubleVectors outputVector;
        cudaLoadCudaFloat4Array( gpu->natoms, 3, gpu->psPosq4,                    outputVector );
        cudaLoadCudaFloatArray( gpu->natoms,  3, amoebaGpu->psInducedDipole,      outputVector );
        cudaLoadCudaFloatArray( gpu->natoms,  3, amoebaGpu->psInducedDipolePolar, outputVector );
        cudaWriteVectorOfDoubleVectorsToFile( "CudaMI_GK", fileId, outputVector );
     }
#endif

   // ---------------------------------------------------------------------------------------
}

void cudaComputeAmoebaMutualInducedAndGkField( amoebaGpuContext amoebaGpu )
{
    if( amoebaGpu->mutualInducedIterativeMethod == 0 ){
        cudaComputeAmoebaMutualInducedAndGkFieldBySOR( amoebaGpu );
    }
}
