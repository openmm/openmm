//-----------------------------------------------------------------------------------------

//-----------------------------------------------------------------------------------------

#include "amoebaGpuTypes.h"
#include "amoebaCudaKernels.h"
#include "kCalculateAmoebaCudaUtilities.h"

#include <stdio.h>

using namespace std;

static __constant__ cudaGmxSimulation cSim;
static __constant__ cudaAmoebaGmxSimulation cAmoebaSim;

void SetCalculateAmoebaCudaMutualInducedFieldSim(amoebaGpuContext amoebaGpu)
{
    cudaError_t status;
    gpuContext gpu = amoebaGpu->gpuContext;
    status         = cudaMemcpyToSymbol(cSim, &gpu->sim, sizeof(cudaGmxSimulation));    
    RTERROR(status, "SetCalculateAmoebaCudaMutualInducedFieldSim: cudaMemcpyToSymbol: SetSim copy to cSim failed");
    status         = cudaMemcpyToSymbol(cAmoebaSim, &amoebaGpu->amoebaSim, sizeof(cudaAmoebaGmxSimulation));    
    RTERROR(status, "SetCalculateAmoebaCudaMutualInducedFieldSim: cudaMemcpyToSymbol: SetSim copy to cAmoebaSim failed");
}

void GetCalculateAmoebaCudaMutualInducedFieldSim(amoebaGpuContext amoebaGpu)
{
    cudaError_t status;
    gpuContext gpu = amoebaGpu->gpuContext;
    status = cudaMemcpyFromSymbol(&gpu->sim, cSim, sizeof(cudaGmxSimulation));    
    RTERROR(status, "GetCalculateAmoebaCudaMutualInducedFieldSim: cudaMemcpyFromSymbol: SetSim copy from cSim failed");
    status = cudaMemcpyFromSymbol(&amoebaGpu->amoebaSim, cAmoebaSim, sizeof(cudaAmoebaGmxSimulation));    
    RTERROR(status, "GetCalculateAmoebaCudaMutualInducedFieldSim: cudaMemcpyFromSymbol: SetSim copy from cAmoebaSim failed");
}

//#define AMOEBA_DEBUG
#undef AMOEBA_DEBUG

#include "kCalculateAmoebaCudaMutualInducedParticle.h"

__device__ void calculateMutualInducedFieldPairIxn_kernel( MutualInducedParticle& atomI, MutualInducedParticle& atomJ,
                                                           float fields[4][3]

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
    float pGamma                                 = atomJ.thole > atomI.thole ? atomI.thole: atomJ.thole;
    float damp                                   = ratio*ratio*ratio*pGamma;
    float dampExp                                = ( (dampProd != 0.0f) && (r < cAmoebaSim.scalingDistanceCutoff) ) ? expf( -damp ) : 0.0f; 

    rr3                                         *= (1.0f - dampExp);
    rr5                                         *= (1.0f - ( 1.0f + damp )*dampExp);
        
    float dDotDelta                              = rr5*(deltaR[0]*atomJ.inducedDipole[0]         + deltaR[1]*atomJ.inducedDipole[1]       + deltaR[2]*atomJ.inducedDipole[2] );
    fields[0][0]                                 = rr3*atomJ.inducedDipole[0] + dDotDelta*deltaR[0];
    fields[0][1]                                 = rr3*atomJ.inducedDipole[1] + dDotDelta*deltaR[1];
    fields[0][2]                                 = rr3*atomJ.inducedDipole[2] + dDotDelta*deltaR[2];
   
    dDotDelta                                    = rr5*(deltaR[0]*atomJ.inducedDipolePolar[0]    + deltaR[1]*atomJ.inducedDipolePolar[1]  + deltaR[2]*atomJ.inducedDipolePolar[2] );
    fields[1][0]                                 = rr3*atomJ.inducedDipolePolar[0] + dDotDelta*deltaR[0];
    fields[1][1]                                 = rr3*atomJ.inducedDipolePolar[1] + dDotDelta*deltaR[1];
    fields[1][2]                                 = rr3*atomJ.inducedDipolePolar[2] + dDotDelta*deltaR[2];
  
    dDotDelta                                    = rr5*(deltaR[0]*atomI.inducedDipole[0]         + deltaR[1]*atomI.inducedDipole[1]       + deltaR[2]*atomI.inducedDipole[2] );
    fields[2][0]                                 = rr3*atomI.inducedDipole[0] + dDotDelta*deltaR[0];
    fields[2][1]                                 = rr3*atomI.inducedDipole[1] + dDotDelta*deltaR[1];
    fields[2][2]                                 = rr3*atomI.inducedDipole[2] + dDotDelta*deltaR[2];
   
    dDotDelta                                    = rr5*(deltaR[0]*atomI.inducedDipolePolar[0]    + deltaR[1]*atomI.inducedDipolePolar[1]  + deltaR[2]*atomI.inducedDipolePolar[2] );
    fields[3][0]                                 = rr3*atomI.inducedDipolePolar[0] + dDotDelta*deltaR[0];
    fields[3][1]                                 = rr3*atomI.inducedDipolePolar[1] + dDotDelta*deltaR[1];
    fields[3][2]                                 = rr3*atomI.inducedDipolePolar[2] + dDotDelta*deltaR[2];
}

// Include versions of the kernels for N^2 calculations.

#define METHOD_NAME(a, b) a##N2##b
#include "kCalculateAmoebaCudaMutualInducedField.h"
#define USE_OUTPUT_BUFFER_PER_WARP
#undef METHOD_NAME
#define METHOD_NAME(a, b) a##N2ByWarp##b
#include "kCalculateAmoebaCudaMutualInducedField.h"

__global__ 
#if (__CUDA_ARCH__ >= 200)
__launch_bounds__(GF1XX_THREADS_PER_BLOCK, 1)
#elif (__CUDA_ARCH__ >= 120)
__launch_bounds__(GT2XX_THREADS_PER_BLOCK, 1)
#else
__launch_bounds__(G8X_THREADS_PER_BLOCK, 1)
#endif
void kInitializeMutualInducedField_kernel(
                   int numberOfAtoms,
                   float* fixedEField,
                   float* fixedEFieldPolar,
                   float* polarizability,
                   float* inducedDipole,
                   float* inducedDipolePolar )
{

    int threadId = __mul24(blockIdx.x,blockDim.x) + threadIdx.x;
    if( threadId >= 3*numberOfAtoms )return;

    fixedEField[threadId]         *= polarizability[threadId];
    inducedDipole[threadId]        = fixedEField[threadId];

    fixedEFieldPolar[threadId]    *= polarizability[threadId];
    inducedDipolePolar[threadId]   = fixedEFieldPolar[threadId];

}

__global__ 
#if (__CUDA_ARCH__ >= 200)
__launch_bounds__(GF1XX_THREADS_PER_BLOCK, 1)
#elif (__CUDA_ARCH__ >= 120)
__launch_bounds__(GT2XX_THREADS_PER_BLOCK, 1)
#else
__launch_bounds__(G8X_THREADS_PER_BLOCK, 1)
#endif
void kReduceMutualInducedFieldDelta_kernel(int numberOfEntries, float* arrayOfDeltas1, float* arrayOfDeltas2, float* epsilon )
{
    extern __shared__ float2 delta[];

    delta[threadIdx.x].x    = 0.0f;
    delta[threadIdx.x].y    = 0.0f;

    unsigned int pos = threadIdx.x;

    // load deltas

    while( pos < numberOfEntries )
    {   
        delta[threadIdx.x].x  += arrayOfDeltas1[pos];
        delta[threadIdx.x].y  += arrayOfDeltas2[pos];
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
        }
        __syncthreads();
    }   

    // set epsilons

    if (threadIdx.x == 0)
    {   
        epsilon[0]  = delta[0].x > delta[0].y ? delta[0].x : delta[0].y;
        epsilon[0]  = 48.033324f*sqrtf( epsilon[0]/( (float) (numberOfEntries/3)) );
    }   
}

__global__ 
#if (__CUDA_ARCH__ >= 200)
__launch_bounds__(GF1XX_THREADS_PER_BLOCK, 1)
#elif (__CUDA_ARCH__ >= 120)
__launch_bounds__(GT2XX_THREADS_PER_BLOCK, 1)
#else
__launch_bounds__(G8X_THREADS_PER_BLOCK, 1)
#endif
void kSorUpdateMutualInducedField_kernel(
                   int numberOfEntries,    float* polarizability,
                   float* inducedDipole, float* inducedDipoleP,
                   float* fixedEField,   float* fixedEFieldP,
                   float* matrixProduct, float* matrixProductP )
{

    float polarSOR = 0.70f;
    int threadId   = __mul24(blockIdx.x,blockDim.x) + threadIdx.x;
    if( threadId  >= 3*numberOfEntries )return;

    float previousDipole                = inducedDipole[threadId];
    float previousDipoleP               = inducedDipoleP[threadId];

    inducedDipole[threadId]             = fixedEField[threadId]     + polarizability[threadId]*matrixProduct[threadId];
    inducedDipoleP[threadId]            = fixedEFieldP[threadId]    + polarizability[threadId]*matrixProductP[threadId];

    inducedDipole[threadId]             = previousDipole   + polarSOR*( inducedDipole[threadId]   - previousDipole  );   
    inducedDipoleP[threadId]            = previousDipoleP  + polarSOR*( inducedDipoleP[threadId]  - previousDipoleP );

    matrixProduct[threadId]             = ( inducedDipole[threadId]  - previousDipole  )*( inducedDipole[threadId]  - previousDipole  );
    matrixProductP[threadId]            = ( inducedDipoleP[threadId] - previousDipoleP )*( inducedDipoleP[threadId] - previousDipoleP );

}

// reduce psWorkArray_3_1
// reduce psWorkArray_3_2

static void kReduceMutualInducedFields(amoebaGpuContext amoebaGpu, CUDAStream<float>* outputArray, CUDAStream<float>* outputPolarArray )
{
    kReduceFields_kernel<<<amoebaGpu->nonbondBlocks, amoebaGpu->fieldReduceThreadsPerBlock>>>(
                               amoebaGpu->paddedNumberOfAtoms*3, amoebaGpu->outputBuffers,
                               amoebaGpu->psWorkArray_3_1->_pDevStream[0], outputArray->_pDevStream[0] );
    LAUNCHERROR("kReduceMI_Fields1");

    kReduceFields_kernel<<<amoebaGpu->nonbondBlocks, amoebaGpu->fieldReduceThreadsPerBlock>>>(
                               amoebaGpu->paddedNumberOfAtoms*3, amoebaGpu->outputBuffers,
                               amoebaGpu->psWorkArray_3_2->_pDevStream[0], outputPolarArray->_pDevStream[0] );
    LAUNCHERROR("kReduceMI_Fields2");
}

/**---------------------------------------------------------------------------------------

   Compute mutual induce field

   @param amoebaGpu        amoebaGpu context

   --------------------------------------------------------------------------------------- */

static void cudaComputeAmoebaMutualInducedFieldMatrixMultiply( amoebaGpuContext amoebaGpu,
                                                               CUDAStream<float>* outputArray, CUDAStream<float>* outputPolarArray )
{
  
  gpuContext gpu    = amoebaGpu->gpuContext;

#ifdef AMOEBA_DEBUG
    int targetAtom    = 1231;
    static const char* methodName = "cudaComputeAmoebaMutualInducedFieldMatrixMultiply";
    static int iteration = 1;
    if( 1 && amoebaGpu->log ){
        (void) fprintf( amoebaGpu->log, "%s\n", methodName );
        (void) fflush( amoebaGpu->log );
    }
    int paddedNumberOfAtoms                    = amoebaGpu->gpuContext->sim.paddedNumberOfAtoms;
    CUDAStream<float4>* debugArray             = new CUDAStream<float4>(paddedNumberOfAtoms*paddedNumberOfAtoms, 1, "DebugArray");
    memset( debugArray->_pSysStream[0],      0, sizeof( float )*4*paddedNumberOfAtoms*paddedNumberOfAtoms);
    debugArray->Upload();
#endif

    kClearFields_3( amoebaGpu, 2 );

    if (gpu->bOutputBufferPerWarp){
        kCalculateAmoebaMutualInducedFieldN2ByWarp_kernel<<<amoebaGpu->nonbondBlocks, amoebaGpu->nonbondThreadsPerBlock, sizeof(MutualInducedParticle)*amoebaGpu->nonbondThreadsPerBlock>>>(
                                                                 amoebaGpu->psWorkUnit->_pDevStream[0],
                                                                 amoebaGpu->psWorkArray_3_1->_pDevStream[0],
#ifdef AMOEBA_DEBUG
                                                                 amoebaGpu->psWorkArray_3_2->_pDevStream[0],
                                                                 debugArray->_pDevStream[0], targetAtom );
#else
                                                                 amoebaGpu->psWorkArray_3_2->_pDevStream[0] );
#endif

    } else {

#ifdef AMOEBA_DEBUG
        (void) fprintf( amoebaGpu->log, "N2 no warp\n" );
        (void) fprintf( amoebaGpu->log, "AmoebaN2Forces_kernel numBlocks=%u numThreads=%u bufferPerWarp=%u atm=%u shrd=%u Ebuf=%u ixnCt=%u workUnits=%u\n",
                        amoebaGpu->nonbondBlocks, amoebaGpu->nonbondThreadsPerBlock, amoebaGpu->bOutputBufferPerWarp,
                        sizeof(MutualInducedParticle), sizeof(MutualInducedParticle)*amoebaGpu->nonbondThreadsPerBlock,
                        amoebaGpu->energyOutputBuffers, (*gpu->psInteractionCount)[0], gpu->sim.workUnits );
        (void) fflush( amoebaGpu->log );
#endif
        kCalculateAmoebaMutualInducedFieldN2_kernel<<<amoebaGpu->nonbondBlocks, amoebaGpu->nonbondThreadsPerBlock, sizeof(MutualInducedParticle)*amoebaGpu->nonbondThreadsPerBlock>>>(
                                                                 amoebaGpu->psWorkUnit->_pDevStream[0],
                                                                 amoebaGpu->psWorkArray_3_1->_pDevStream[0],
#ifdef AMOEBA_DEBUG
                                                                 amoebaGpu->psWorkArray_3_2->_pDevStream[0],
                                                                 debugArray->_pDevStream[0], targetAtom );
#else
                                                                 amoebaGpu->psWorkArray_3_2->_pDevStream[0] );
#endif


    }
    LAUNCHERROR("kCalculateAmoebaMutualInducedField");

    kReduceMutualInducedFields( amoebaGpu, outputArray, outputPolarArray );

#ifdef AMOEBA_DEBUG
    amoebaGpu->psWorkArray_3_1->Download();
    amoebaGpu->psWorkArray_3_2->Download();

    if( amoebaGpu->log && iteration == -1 ){
        (void) fprintf( amoebaGpu->log, "Finished MI kernel execution %d\n", iteration ); (void) fflush( amoebaGpu->log );
        outputArray->Download();
        outputPolarArray->Download();
        debugArray->Download();

        int maxPrint        = 1400;
        for( int ii = 0; ii < gpu->natoms; ii++ ){
           (void) fprintf( amoebaGpu->log, "%5d ", ii); 

            int indexOffset     = ii*3;
    
           // MI

           (void) fprintf( amoebaGpu->log,"Mult[%16.9e %16.9e %16.9e] ",
                           outputArray->_pSysStream[0][indexOffset],
                           outputArray->_pSysStream[0][indexOffset+1],
                           outputArray->_pSysStream[0][indexOffset+2] );
    
           // MI polar

           (void) fprintf( amoebaGpu->log,"MultP[%16.9e %16.9e %16.9e] ",
                           outputPolarArray->_pSysStream[0][indexOffset],
                           outputPolarArray->_pSysStream[0][indexOffset+1],
                           outputPolarArray->_pSysStream[0][indexOffset+2] );

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
               float sums[4][3] = { { 0.0f, 0.0f, 0.0f },
                                    { 0.0f, 0.0f, 0.0f },
                                    { 0.0f, 0.0f, 0.0f },
                                    { 0.0f, 0.0f, 0.0f } };
               (void) fprintf( amoebaGpu->log,"\n" );
               int paddedNumberOfAtoms                    = amoebaGpu->gpuContext->sim.paddedNumberOfAtoms;
               unsigned int count                         = 0;
               for( int jj = 0; jj < gpu->natoms; jj++ ){
                   int debugIndex = jj;
                   (void) fprintf( amoebaGpu->log,"%4d %4d Pint [%16.9e %16.9e %16.9e %16.9e] ",
                                   ii, jj,
                                   debugArray->_pSysStream[0][debugIndex].x, debugArray->_pSysStream[0][debugIndex].y,
                                   debugArray->_pSysStream[0][debugIndex].z, debugArray->_pSysStream[0][debugIndex].w );

                   //debugIndex += gpu->natoms;
                   debugIndex += paddedNumberOfAtoms;
                   (void) fprintf( amoebaGpu->log,"[%16.9e %16.9e %16.9e] ",
                                   debugArray->_pSysStream[0][debugIndex].x, debugArray->_pSysStream[0][debugIndex].y, debugArray->_pSysStream[0][debugIndex].z );

                   int index = 0;
                   sums[index][0] += debugArray->_pSysStream[0][debugIndex].x; 
                   sums[index][1] += debugArray->_pSysStream[0][debugIndex].y; 
                   sums[index][2] += debugArray->_pSysStream[0][debugIndex].z; 
                   
                   if( count && ( (count % 31) == 0) ){
                      static float saveSum[3] = { 0.0f, 0.0f, 0.0f };
                      (void) fprintf( amoebaGpu->log,"Block sum [%16.9e %16.9e %16.9e] ",
                                      sums[index][0] - saveSum[0], sums[index][1] - saveSum[1], sums[index][2] - saveSum[2] );
                      saveSum[0] = sums[index][0];
                      saveSum[1] = sums[index][1];
                      saveSum[2] = sums[index][2];
                     
                   }
                   

                   debugIndex += paddedNumberOfAtoms;
                   (void) fprintf( amoebaGpu->log,"[%16.9e %16.9e %16.9e] ",
                                   debugArray->_pSysStream[0][debugIndex].x, debugArray->_pSysStream[0][debugIndex].y, debugArray->_pSysStream[0][debugIndex].z );

                   index++;
                   sums[index][0] += debugArray->_pSysStream[0][debugIndex].x; 
                   sums[index][1] += debugArray->_pSysStream[0][debugIndex].y; 
                   sums[index][2] += debugArray->_pSysStream[0][debugIndex].z; 

                   if( count && ( (count % 31) == 0) ){
                      static float saveSum[3] = { 0.0f, 0.0f, 0.0f };
                      (void) fprintf( amoebaGpu->log,"Block sumP [%16.9e %16.9e %16.9e] ",
                                      sums[index][0] - saveSum[0], sums[index][1] - saveSum[1], sums[index][2] - saveSum[2] );
                      saveSum[0] = sums[index][0];
                      saveSum[1] = sums[index][1];
                      saveSum[2] = sums[index][2];
                   }
                   (void) fprintf( amoebaGpu->log,"\n" );
                   count++;
               }

               (void) fprintf( amoebaGpu->log,"\n" );
               int index = 0;
               (void) fprintf( amoebaGpu->log,"Sum1 [%16.9e %16.9e %16.9e]\n", sums[index][0], sums[index][1],sums[index][2] ); index++;
               (void) fprintf( amoebaGpu->log,"Sum2 [%16.9e %16.9e %16.9e]\n", sums[index][0], sums[index][1],sums[index][2] ); index++;
               (void) fprintf( amoebaGpu->log,"Sum3 [%16.9e %16.9e %16.9e]\n", sums[index][0], sums[index][1],sums[index][2] ); index++;
               (void) fprintf( amoebaGpu->log,"Sum4 [%16.9e %16.9e %16.9e]\n", sums[index][0], sums[index][1],sums[index][2] ); index++;
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

static void cudaComputeAmoebaMutualInducedFieldBySOR( amoebaGpuContext amoebaGpu )
{
  
   // ---------------------------------------------------------------------------------------

#ifdef AMOEBA_DEBUG
    static const char* methodName = "cudaComputeAmoebaMutualInducedFieldBySOR";
    static int timestep = 0;
    std::vector<int> fileId;
    timestep++;
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
    if( amoebaGpu->log ){
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

    kInitializeMutualInducedField_kernel<<< numBlocks, numThreads >>>(
         gpu->natoms,
         amoebaGpu->psE_Field->_pDevStream[0],
         amoebaGpu->psE_FieldPolar->_pDevStream[0],
         amoebaGpu->psPolarizability->_pDevStream[0],
         amoebaGpu->psInducedDipole->_pDevStream[0],
         amoebaGpu->psInducedDipolePolar->_pDevStream[0] );
    LAUNCHERROR("AmoebaMutualInducedFieldSetup");  

#ifdef AMOEBA_DEBUG
    if( amoebaGpu->log ){

        amoebaGpu->psE_Field->Download();
        amoebaGpu->psE_FieldPolar->Download();
        amoebaGpu->psInducedDipole->Download(),
        amoebaGpu->psInducedDipolePolar->Download();
        amoebaGpu->psPolarizability->Download();
        (void) fprintf( amoebaGpu->log, "%s Initial setup for matrix multiply\n", methodName );
        int offset   = 0;
        int maxPrint = 20000;
        for( int ii = 0; ii < gpu->natoms; ii++ ){
            (void) fprintf( amoebaGpu->log, "%4d pol=%12.4e ", ii, 
                            amoebaGpu->psPolarizability->_pSysStream[0][offset] );
            if( amoebaGpu->psPolarizability->_pSysStream[0][offset] != amoebaGpu->psPolarizability->_pSysStream[0][offset+1] ||
                amoebaGpu->psPolarizability->_pSysStream[0][offset] != amoebaGpu->psPolarizability->_pSysStream[0][offset+2] ){
                (void) fprintf( amoebaGpu->log, "PolX!!! %12.4e %12.4e ", amoebaGpu->psPolarizability->_pSysStream[0][offset+1], amoebaGpu->psPolarizability->_pSysStream[0][offset+2] ); 
            }

            (void) fprintf( amoebaGpu->log," E[%14.6e %14.6e %14.6e] Mi[%14.6e %14.6e %14.6e] ",
                            amoebaGpu->psE_Field->_pSysStream[0][offset],       amoebaGpu->psE_Field->_pSysStream[0][offset+1],       amoebaGpu->psE_Field->_pSysStream[0][offset+2],
                            amoebaGpu->psInducedDipole->_pSysStream[0][offset], amoebaGpu->psInducedDipole->_pSysStream[0][offset+1], amoebaGpu->psInducedDipole->_pSysStream[0][offset+2] );
            (void) fprintf( amoebaGpu->log,"Ep[%14.6e %14.6e %14.6e] Mip[%14.6e %14.6e %14.6e]\n",
                            amoebaGpu->psE_FieldPolar->_pSysStream[0][offset],       amoebaGpu->psE_FieldPolar->_pSysStream[0][offset+1],       amoebaGpu->psE_FieldPolar->_pSysStream[0][offset+2],
                            amoebaGpu->psInducedDipolePolar->_pSysStream[0][offset], amoebaGpu->psInducedDipolePolar->_pSysStream[0][offset+1], amoebaGpu->psInducedDipolePolar->_pSysStream[0][offset+2] );
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

        cudaComputeAmoebaMutualInducedFieldMatrixMultiply( amoebaGpu, amoebaGpu->psWorkVector[0],  amoebaGpu->psWorkVector[1] );
        LAUNCHERROR("cudaComputeAmoebaMutualInducedFieldMatrixMultiply Loop\n");  

        // post matrix multiply

        kSorUpdateMutualInducedField_kernel<<< numBlocks, numThreads >>>(
           gpu->natoms, amoebaGpu->psPolarizability->_pDevStream[0],
           amoebaGpu->psInducedDipole->_pDevStream[0], amoebaGpu->psInducedDipolePolar->_pDevStream[0],
           amoebaGpu->psE_Field->_pDevStream[0],       amoebaGpu->psE_FieldPolar->_pDevStream[0],
           amoebaGpu->psWorkVector[0]->_pDevStream[0],     amoebaGpu->psWorkVector[1]->_pDevStream[0] );
        LAUNCHERROR("kSorUpdateMutualInducedField");  

        // get total epsilon -- performing sums on gpu

        kReduceMutualInducedFieldDelta_kernel<<<1, amoebaGpu->epsilonThreadsPerBlock, 2*sizeof(float)*amoebaGpu->epsilonThreadsPerBlock>>>(
           3*gpu->natoms, amoebaGpu->psWorkVector[0]->_pDevStream[0], amoebaGpu->psWorkVector[1]->_pDevStream[0],
           amoebaGpu->psCurrentEpsilon->_pDevStream[0] );
        LAUNCHERROR("kReduceMutualInducedFieldDelta");

        if( amoebaGpu->log ){
            trackMutualInducedIterations( amoebaGpu, iteration);
        }

        // Debye=48.033324f
        amoebaGpu->psCurrentEpsilon->Download();
        float currentEpsilon          = amoebaGpu->psCurrentEpsilon->_pSysStream[0][0];
        amoebaGpu->mutualInducedCurrentEpsilon   = currentEpsilon;

        if( iteration > amoebaGpu->mutualInducedMaxIterations || amoebaGpu->mutualInducedCurrentEpsilon < amoebaGpu->mutualInducedTargetEpsilon ){ 
            done = 1;
        }

#ifdef AMOEBA_DEBUG
        if( amoebaGpu->log ){
           amoebaGpu->psInducedDipole->Download();
           amoebaGpu->psInducedDipolePolar->Download();
           (void) fprintf( amoebaGpu->log, "%s iteration=%3d eps %14.6e done=%d\n",
                           methodName, iteration, amoebaGpu->mutualInducedCurrentEpsilon, done );
           (void) fflush( amoebaGpu->log );

            int offset   = 0;
            int maxPrint = 20;
            for( int ii = 0; ii < gpu->natoms; ii++ ){
                (void) fprintf( amoebaGpu->log, "%4d ", ii ); 
    
                (void) fprintf( amoebaGpu->log," Mi[%14.6e %14.6e %14.6e] ",
                                amoebaGpu->psInducedDipole->_pSysStream[0][offset], amoebaGpu->psInducedDipole->_pSysStream[0][offset+1], amoebaGpu->psInducedDipole->_pSysStream[0][offset+2] );
                (void) fprintf( amoebaGpu->log,"Mip[%14.6e %14.6e %14.6e]\n",
                                amoebaGpu->psInducedDipolePolar->_pSysStream[0][offset], amoebaGpu->psInducedDipolePolar->_pSysStream[0][offset+1], amoebaGpu->psInducedDipolePolar->_pSysStream[0][offset+2] );
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

#ifdef AMOEBA_DEBUG
    if( 0 ){
        std::vector<int> fileId;
        //fileId.push_back( 0 );
        VectorOfDoubleVectors outputVector;
//        cudaLoadCudaFloat4Array( gpu->natoms, 3, gpu->psPosq4,                    outputVector );
        cudaLoadCudaFloatArray( gpu->natoms,  3, amoebaGpu->psInducedDipole,      outputVector, NULL );
        cudaLoadCudaFloatArray( gpu->natoms,  3, amoebaGpu->psInducedDipolePolar, outputVector, NULL );
        cudaWriteVectorOfDoubleVectorsToFile( "CudaMI", fileId, outputVector );
     }

#endif

   // ---------------------------------------------------------------------------------------
}

void cudaComputeAmoebaMutualInducedField( amoebaGpuContext amoebaGpu )
{
    if( amoebaGpu->mutualInducedIterativeMethod == 0 ){
        cudaComputeAmoebaMutualInducedFieldBySOR( amoebaGpu );
    }
}
