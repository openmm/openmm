//-----------------------------------------------------------------------------------------

//-----------------------------------------------------------------------------------------

#include "amoebaGpuTypes.h"
#include "amoebaCudaKernels.h"
#include "kCalculateAmoebaCudaUtilities.h"

#include <stdio.h>

using namespace std;

static __constant__ cudaGmxSimulation cSim;
static __constant__ cudaAmoebaGmxSimulation cAmoebaSim;

void SetCalculateAmoebaCudaPmeMutualInducedFieldSim(amoebaGpuContext amoebaGpu)
{
    cudaError_t status;
    gpuContext gpu = amoebaGpu->gpuContext;
    status         = cudaMemcpyToSymbol(cSim, &gpu->sim, sizeof(cudaGmxSimulation));    
    RTERROR(status, "SetCalculateAmoebaCudaPmeMutualInducedFieldSim: cudaMemcpyToSymbol: SetSim copy to cSim failed");
    status         = cudaMemcpyToSymbol(cAmoebaSim, &amoebaGpu->amoebaSim, sizeof(cudaAmoebaGmxSimulation));    
    RTERROR(status, "SetCalculateAmoebaCudaPmeMutualInducedFieldSim: cudaMemcpyToSymbol: SetSim copy to cAmoebaSim failed");
}

void GetCalculateAmoebaCudaPmeMutualInducedFieldSim(amoebaGpuContext amoebaGpu)
{
    cudaError_t status;
    gpuContext gpu = amoebaGpu->gpuContext;
    status = cudaMemcpyFromSymbol(&gpu->sim, cSim, sizeof(cudaGmxSimulation));    
    RTERROR(status, "GetCalculateAmoebaCudaPmeMutualInducedFieldSim: cudaMemcpyFromSymbol: SetSim copy from cSim failed");
    status = cudaMemcpyFromSymbol(&amoebaGpu->amoebaSim, cAmoebaSim, sizeof(cudaAmoebaGmxSimulation));    
    RTERROR(status, "GetCalculateAmoebaCudaPmeMutualInducedFieldSim: cudaMemcpyFromSymbol: SetSim copy from cAmoebaSim failed");
}

//#define AMOEBA_DEBUG
#undef AMOEBA_DEBUG

#include "kCalculateAmoebaCudaMutualInducedParticle.h"

// file includes FixedFieldParticle struct definition/load/unload struct and body kernel for fixed E-field

__device__ void calculatePmeDirectMutualInducedFieldPairIxn_kernel( MutualInducedParticle& atomI, MutualInducedParticle& atomJ,
                                                                    float uscale, float fields[4][3]
#ifdef AMOEBA_DEBUG
                                                            , float4* pullBack
#endif

 ){

    // compute the real space portion of the Ewald summation
  
    float xr          = atomJ.x - atomI.x;
    float yr          = atomJ.y - atomI.y;
    float zr          = atomJ.z - atomI.z;

    // periodic boundary conditions

    xr               -= floor(xr*cSim.invPeriodicBoxSizeX+0.5f)*cSim.periodicBoxSizeX;
    yr               -= floor(yr*cSim.invPeriodicBoxSizeY+0.5f)*cSim.periodicBoxSizeY;
    zr               -= floor(zr*cSim.invPeriodicBoxSizeZ+0.5f)*cSim.periodicBoxSizeZ;

    float r2          = xr*xr + yr* yr + zr*zr;
    float r           = sqrtf(r2);

    // calculate the error function damping terms

    float ralpha      = cSim.alphaEwald*r;
    float bn[3];

    bn[0]             = erfc(ralpha)/r;
    float alsq2       = 2.0f*cSim.alphaEwald*cSim.alphaEwald;
    float alsq2n      = 1.0f/(cAmoebaSim.sqrtPi*cSim.alphaEwald);
    float exp2a       = exp(-(ralpha*ralpha));
    alsq2n           *= alsq2;
    bn[1]             = (bn[0]+alsq2n*exp2a)/r2;

    alsq2n           *= alsq2;
    bn[2]             = (3.0f*bn[1]+alsq2n*exp2a)/r2;

    // compute the error function scaled and unscaled terms

    float scale3      = 1.0f;
    float scale5      = 1.0f;
    float damp        = atomI.damp*atomJ.damp;
    if( damp != 0.0f ){

        float ratio  = (r/damp);
              ratio  = ratio*ratio*ratio;
        float pgamma = atomI.thole < atomJ.thole ? atomI.thole : atomJ.thole;
              damp   = -pgamma*ratio;

        if( damp > -50.0f) {
            float expdamp = exp(damp);
            scale3        = 1.0f - expdamp;
            scale5        = 1.0f - expdamp*(1.0f-damp);
        }
    }
    float dsc3        = uscale*scale3;
    float dsc5        = uscale*scale5;

    float r3          = (r*r2);
    float r5          = (r3*r2);
    float rr3         = (1.0f-dsc3)/r3;
    float rr5         = 3.0f * (1.0f-dsc5)/r5;

    float duir        = atomI.inducedDipole[0]*xr      + atomI.inducedDipole[1]*yr      + atomI.inducedDipole[2]*zr;
    float dukr        = atomJ.inducedDipole[0]*xr      + atomJ.inducedDipole[1]*yr      + atomJ.inducedDipole[2]*zr;

    float puir        = atomI.inducedDipolePolar[0]*xr + atomI.inducedDipolePolar[1]*yr + atomI.inducedDipolePolar[2]*zr;
    float pukr        = atomJ.inducedDipolePolar[0]*xr + atomJ.inducedDipolePolar[1]*yr + atomJ.inducedDipolePolar[2]*zr;

    float fimd[3],fkmd[3];
    float fimp[3],fkmp[3];
    float fid[3],fkd[3];
    float fip[3],fkp[3];

    bn[1]            *= -1.0f;

    fimd[0]           = bn[1]*atomJ.inducedDipole[0]      + bn[2]*dukr*xr;
    fimd[1]           = bn[1]*atomJ.inducedDipole[1]      + bn[2]*dukr*yr;
    fimd[2]           = bn[1]*atomJ.inducedDipole[2]      + bn[2]*dukr*zr;

    fkmd[0]           = bn[1]*atomI.inducedDipole[0]      + bn[2]*duir*xr;
    fkmd[1]           = bn[1]*atomI.inducedDipole[1]      + bn[2]*duir*yr;
    fkmd[2]           = bn[1]*atomI.inducedDipole[2]      + bn[2]*duir*zr;

    fimp[0]           = bn[1]*atomJ.inducedDipolePolar[0] + bn[2]*pukr*xr;
    fimp[1]           = bn[1]*atomJ.inducedDipolePolar[1] + bn[2]*pukr*yr;
    fimp[2]           = bn[1]*atomJ.inducedDipolePolar[2] + bn[2]*pukr*zr;

    fkmp[0]           = bn[1]*atomI.inducedDipolePolar[0] + bn[2]*puir*xr;
    fkmp[1]           = bn[1]*atomI.inducedDipolePolar[1] + bn[2]*puir*yr;
    fkmp[2]           = bn[1]*atomI.inducedDipolePolar[2] + bn[2]*puir*zr;

    rr3              *= -1.0f;;
    fid[0]            = rr3*atomJ.inducedDipole[0]      + rr5*dukr*xr;
    fid[1]            = rr3*atomJ.inducedDipole[1]      + rr5*dukr*yr;
    fid[2]            = rr3*atomJ.inducedDipole[2]      + rr5*dukr*zr;

    fkd[0]            = rr3*atomI.inducedDipole[0]      + rr5*duir*xr;
    fkd[1]            = rr3*atomI.inducedDipole[1]      + rr5*duir*yr;
    fkd[2]            = rr3*atomI.inducedDipole[2]      + rr5*duir*zr;

    fip[0]            = rr3*atomJ.inducedDipolePolar[0] + rr5*pukr*xr;
    fip[1]            = rr3*atomJ.inducedDipolePolar[1] + rr5*pukr*yr;
    fip[2]            = rr3*atomJ.inducedDipolePolar[2] + rr5*pukr*zr;

    fkp[0]            = rr3*atomI.inducedDipolePolar[0] + rr5*puir*xr;
    fkp[1]            = rr3*atomI.inducedDipolePolar[1] + rr5*puir*yr;
    fkp[2]            = rr3*atomI.inducedDipolePolar[2] + rr5*puir*zr;

    // increment the field at each site due to this interaction

    if( r2 <= cAmoebaSim.cutoffDistance2 ){

        fields[0][0]       = fimd[0] - fid[0];
        fields[1][0]       = fkmd[0] - fkd[0];
        fields[2][0]       = fimp[0] - fip[0];
        fields[3][0]       = fkmp[0] - fkp[0];
    
        fields[0][1]       = fimd[1] - fid[1];
        fields[1][1]       = fkmd[1] - fkd[1];
        fields[2][1]       = fimp[1] - fip[1];
        fields[3][1]       = fkmp[1] - fkp[1];
    
        fields[0][2]       = fimd[2] - fid[2];
        fields[1][2]       = fkmd[2] - fkd[2];
        fields[2][2]       = fimp[2] - fip[2];
        fields[3][2]       = fkmp[2] - fkp[2];
 
    } else {

        fields[0][0]       = 0.0f;
        fields[1][0]       = 0.0f;
        fields[2][0]       = 0.0f;
        fields[3][0]       = 0.0f;
    
        fields[0][1]       = 0.0f;
        fields[1][1]       = 0.0f;
        fields[2][1]       = 0.0f;
        fields[3][1]       = 0.0f;
    
        fields[0][2]       = 0.0f;
        fields[1][2]       = 0.0f;
        fields[2][2]       = 0.0f;
        fields[3][2]       = 0.0f;
    }
#ifdef AMOEBA_DEBUG
    pullBack[0].x = xr;
    pullBack[0].y = yr;
    pullBack[0].z = zr;
    pullBack[0].w = r2;

    pullBack[1].x = alsq2;
    pullBack[1].y = bn[0];
    pullBack[1].z = bn[2];
    pullBack[1].w = exp2a;

/*
    pullBack[1].x = atomJ.x - atomI.x;
    pullBack[1].y = atomJ.y - atomI.y;
    pullBack[1].z = atomJ.z - atomI.z;
    pullBack[1].w = (atomJ.x - atomI.x)*(atomJ.x - atomI.x) + (atomJ.y - atomI.y)*(atomJ.y - atomI.y)+ (atomJ.z - atomI.z)*(atomJ.z - atomI.z);
    pullBack[1].x = scale3;
    pullBack[1].y = scale5;
    pullBack[1].z = scale7;
*/
#endif
}

// Include versions of the kernels for N^2 calculations.

#define METHOD_NAME(a, b) a##N2##b
#include "kCalculateAmoebaCudaPmeMutualInducedField.h"
#define USE_OUTPUT_BUFFER_PER_WARP
#undef METHOD_NAME
#define METHOD_NAME(a, b) a##N2ByWarp##b
#include "kCalculateAmoebaCudaPmeMutualInducedField.h"

__global__ 
#if (__CUDA_ARCH__ >= 200)
__launch_bounds__(GF1XX_THREADS_PER_BLOCK, 1)
#elif (__CUDA_ARCH__ >= 130)
__launch_bounds__(GT2XX_THREADS_PER_BLOCK, 1)
#else
__launch_bounds__(G8X_THREADS_PER_BLOCK, 1)
#endif
static void kInitializeMutualInducedField_kernel(
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
#elif (__CUDA_ARCH__ >= 130)
__launch_bounds__(GT2XX_THREADS_PER_BLOCK, 1)
#else
__launch_bounds__(G8X_THREADS_PER_BLOCK, 1)
#endif
static void kReduceMutualInducedFieldDelta_kernel(int numberOfEntries, float* arrayOfDeltas1, float* arrayOfDeltas2, float* epsilon )
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
        epsilon[0]  = 4.8033324f*sqrtf( epsilon[0]/( (float) (numberOfEntries/3)) );
#ifdef AMOEBA_DEBUG
        epsilon[1]  = 4.8033324f*sqrtf( delta[0].x/( (float) (numberOfEntries/3)) );
        epsilon[2]  = 4.8033324f*sqrtf( delta[0].y/( (float) (numberOfEntries/3)) );
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
static void kSorUpdateMutualInducedField_kernel(
                   int numberOfEntries,    float* polarizability,
                   float* inducedDipole, float* inducedDipoleP,
                   float* fixedEField,   float* fixedEFieldP,
                   float* matrixProduct, float* matrixProductP )
{

    int threadId         = __mul24(blockIdx.x,blockDim.x) + threadIdx.x;
    if( threadId  >= 3*numberOfEntries )return;

    float previousDipole                = inducedDipole[threadId];
    float previousDipoleP               = inducedDipoleP[threadId];

    // add self terms to fields

    const float term                    = (4.0f/3.0f)*(cSim.alphaEwald*cSim.alphaEwald*cSim.alphaEwald)/cAmoebaSim.sqrtPi;
    matrixProduct[threadId]            +=  term*previousDipole;
    matrixProductP[threadId]           +=  term*previousDipoleP;

    inducedDipole[threadId]             = fixedEField[threadId]     + polarizability[threadId]*matrixProduct[threadId];
    inducedDipoleP[threadId]            = fixedEFieldP[threadId]    + polarizability[threadId]*matrixProductP[threadId];

    const float polarSOR                = 0.70f;
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
    LAUNCHERROR("kReducePmeMI_Fields1");

    kReduceFields_kernel<<<amoebaGpu->nonbondBlocks, amoebaGpu->fieldReduceThreadsPerBlock>>>(
                               amoebaGpu->paddedNumberOfAtoms*3, amoebaGpu->outputBuffers,
                               amoebaGpu->psWorkArray_3_2->_pDevStream[0], outputPolarArray->_pDevStream[0] );
    LAUNCHERROR("kReducePmeMI_Fields2");
}

/**---------------------------------------------------------------------------------------

   Compute mutual induce field

   @param amoebaGpu        amoebaGpu context

   --------------------------------------------------------------------------------------- */

static void cudaComputeAmoebaPmeMutualInducedFieldMatrixMultiply( amoebaGpuContext amoebaGpu,
                                                                  CUDAStream<float>* outputArray, CUDAStream<float>* outputPolarArray )
{
  
  gpuContext gpu    = amoebaGpu->gpuContext;

#ifdef AMOEBA_DEBUG
    int targetAtom    = 0;
    static const char* methodName = "cudaComputeAmoebaPmeMutualInducedFieldMatrixMultiply";
    static int iteration = 1;
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

    kClearFields_3( amoebaGpu, 2 );

    if (gpu->bOutputBufferPerWarp){
        kCalculateAmoebaPmeMutualInducedFieldN2ByWarp_kernel<<<amoebaGpu->nonbondBlocks, amoebaGpu->nonbondThreadsPerBlock, sizeof(MutualInducedParticle)*amoebaGpu->nonbondThreadsPerBlock>>>(
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
        kCalculateAmoebaPmeMutualInducedFieldN2_kernel<<<amoebaGpu->nonbondBlocks, amoebaGpu->nonbondThreadsPerBlock,
                                                         sizeof(MutualInducedParticle)*amoebaGpu->nonbondThreadsPerBlock>>>(
                                                                 amoebaGpu->psWorkUnit->_pDevStream[0],
                                                                 amoebaGpu->psWorkArray_3_1->_pDevStream[0],
#ifdef AMOEBA_DEBUG
                                                                 amoebaGpu->psWorkArray_3_2->_pDevStream[0],
                                                                 debugArray->_pDevStream[0], targetAtom );
#else
                                                                 amoebaGpu->psWorkArray_3_2->_pDevStream[0] );
#endif


    }
    LAUNCHERROR("kCalculateAmoebaPmeMutualInducedField");

    kReduceMutualInducedFields( amoebaGpu, outputArray, outputPolarArray );

#ifdef AMOEBA_DEBUG
    if( amoebaGpu->log && iteration == 1 ){
        (void) fprintf( amoebaGpu->log, "Finished maxtrixMultiply kernel execution %d -- Direct only -- self added in kSorUpdateMutualInducedField_kernel\n",
                        iteration ); (void) fflush( amoebaGpu->log );
        outputArray->Download();
        outputPolarArray->Download();
        debugArray->Download();
        int maxPrint = 5;
        for( int ii = 0; ii < gpu->natoms; ii++ ){
            (void) fprintf( amoebaGpu->log, "%5d ", ii); 
 
             int indexOffset     = ii*3;
     
            // MI
 
            (void) fprintf( amoebaGpu->log,"Mult[%16.9e %16.9e %16.9e] ",
                            outputArray->_pSysStream[0][indexOffset],
                            outputArray->_pSysStream[0][indexOffset+1],
                            outputArray->_pSysStream[0][indexOffset+2] );
     
            // MI polar
 
            (void) fprintf( amoebaGpu->log,"MultP[%16.9e %16.9e %16.9e]\n",
                            outputPolarArray->_pSysStream[0][indexOffset],
                            outputPolarArray->_pSysStream[0][indexOffset+1],
                            outputPolarArray->_pSysStream[0][indexOffset+2] );
            if( ii == maxPrint && (gpu->natoms - maxPrint) > ii ){
                ii = gpu->natoms - maxPrint;
            }

        }
/*
        int paddedNumberOfAtoms = amoebaGpu->gpuContext->sim.paddedNumberOfAtoms;
        for( int jj = 0; jj < gpu->natoms; jj++ ){
            int debugIndex = jj; 
            (void) fprintf( amoebaGpu->log,"%5d PmeMIMult\n", jj );
            for( int kk = 0; kk < 7; kk++ ){
                (void) fprintf( amoebaGpu->log,"[%16.9e %16.9e %16.9e %16.9e]\n",
                                debugArray->_pSysStream[0][debugIndex].x, debugArray->_pSysStream[0][debugIndex].y,
                                debugArray->_pSysStream[0][debugIndex].z, debugArray->_pSysStream[0][debugIndex].w );
                debugIndex += paddedNumberOfAtoms;
            }
            (void) fprintf( amoebaGpu->log,"\n" );

        }
*/
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

//#define GET_EFIELD_FROM_FILE
#ifdef GET_EFIELD_FROM_FILE
#include <stdlib.h>
#endif

static void cudaComputeAmoebaPmeMutualInducedFieldBySOR( amoebaGpuContext amoebaGpu )
{
  
   // ---------------------------------------------------------------------------------------

#ifdef AMOEBA_DEBUG
    static const char* methodName = "cudaComputeAmoebaPmeMutualInducedFieldBySOR";
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

#ifdef GET_EFIELD_FROM_FILE
    std::string fileName = "waterEField.txt";
    StringVectorVector fileContents;
    readFile( fileName, fileContents );
    unsigned int offset  = 0;
    (void) fprintf( amoebaGpu->log, "Read file: %s %u\n", fileName.c_str(), fileContents.size() ); fflush(  amoebaGpu->log );
    for( unsigned int ii = 1; ii < fileContents.size()-1; ii++ ){

        StringVector lineTokens     = fileContents[ii];
        unsigned int lineTokenIndex = 1;

        // (void) fprintf( amoebaGpu->log, "   %u %s %s\n", ii, lineTokens[0].c_str(), lineTokens[lineTokenIndex].c_str() ); fflush(  amoebaGpu->log );
        amoebaGpu->psE_Field->_pSysStream[0][offset++]       = static_cast<float>(atof(lineTokens[lineTokenIndex++].c_str()));
        amoebaGpu->psE_Field->_pSysStream[0][offset++]       = static_cast<float>(atof(lineTokens[lineTokenIndex++].c_str()));
        amoebaGpu->psE_Field->_pSysStream[0][offset++]       = static_cast<float>(atof(lineTokens[lineTokenIndex++].c_str())); 
        offset                                              -= 3;        
        amoebaGpu->psE_FieldPolar->_pSysStream[0][offset++]  = static_cast<float>(atof(lineTokens[lineTokenIndex++].c_str()));
        amoebaGpu->psE_FieldPolar->_pSysStream[0][offset++]  = static_cast<float>(atof(lineTokens[lineTokenIndex++].c_str()));
        amoebaGpu->psE_FieldPolar->_pSysStream[0][offset++]  = static_cast<float>(atof(lineTokens[lineTokenIndex++].c_str()));
    }
    float conversion = 100.0f;
    for( int ii = 0; ii < 3*gpu->natoms; ii++ ){
        amoebaGpu->psE_Field->_pSysStream[0][ii]       *= conversion;
        amoebaGpu->psE_FieldPolar->_pSysStream[0][ii]  *= conversion;
    }
    amoebaGpu->gpuContext->sim.alphaEwald = 5.4459052e+00f;
    (void) fprintf( amoebaGpu->log, "cudaComputeAmoebaPmeMutualInducedFieldBySOR: forcing alphaEwald=%15.7e\n", amoebaGpu->gpuContext->sim.alphaEwald );
    SetCalculateAmoebaCudaPmeMutualInducedFieldSim(amoebaGpu);
    amoebaGpu->psE_Field->Upload();
    amoebaGpu->psE_FieldPolar->Upload();
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
    LAUNCHERROR("AmoebaPmeMutualInducedFieldSetup");  

#ifdef AMOEBA_DEBUG
    if( amoebaGpu->log ){

        amoebaGpu->psE_Field->Download();
        amoebaGpu->psE_FieldPolar->Download();
        amoebaGpu->psInducedDipole->Download(),
        amoebaGpu->psInducedDipolePolar->Download();
        amoebaGpu->psPolarizability->Download();
        (void) fprintf( amoebaGpu->log, "%s Initial setup for matrix multiply\n", methodName );
        int offset   = 0;
        int maxPrint = 10;
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

        cudaComputeAmoebaPmeMutualInducedFieldMatrixMultiply( amoebaGpu, amoebaGpu->psWorkVector[0],  amoebaGpu->psWorkVector[1] );
        kCalculateAmoebaPMEInducedDipoleField( amoebaGpu );
        LAUNCHERROR("cudaComputeAmoebaPmeMutualInducedFieldMatrixMultiply Loop\n");  


#ifdef GET_EFIELD_FROM_FILE
{
    std::string fileName = "waterInduceRecip.txt";
    StringVectorVector fileContents;
    readFile( fileName, fileContents );
    unsigned int offset  = 0;
    amoebaGpu->psWorkVector[0]->Download();
    amoebaGpu->psWorkVector[1]->Download();
    (void) fprintf( amoebaGpu->log, "Read file: %s %u\n", fileName.c_str(), fileContents.size() ); fflush(  amoebaGpu->log );
    float conversion = 100.0f;
    for( unsigned int ii = 1; ii < fileContents.size()-1; ii++ ){

        StringVector lineTokens     = fileContents[ii];
        unsigned int lineTokenIndex = 1;

        // (void) fprintf( amoebaGpu->log, "   %u %s %s\n", ii, lineTokens[0].c_str(), lineTokens[lineTokenIndex].c_str() ); fflush(  amoebaGpu->log );
        amoebaGpu->psWorkVector[0]->_pSysStream[0][offset++]       += conversion*static_cast<float>(atof(lineTokens[lineTokenIndex++].c_str()));
        amoebaGpu->psWorkVector[0]->_pSysStream[0][offset++]       += conversion*static_cast<float>(atof(lineTokens[lineTokenIndex++].c_str()));
        amoebaGpu->psWorkVector[0]->_pSysStream[0][offset++]       += conversion*static_cast<float>(atof(lineTokens[lineTokenIndex++].c_str())); 
        offset                                                     -= 3;        
        amoebaGpu->psWorkVector[1]->_pSysStream[0][offset++]       += conversion*static_cast<float>(atof(lineTokens[lineTokenIndex++].c_str()));
        amoebaGpu->psWorkVector[1]->_pSysStream[0][offset++]       += conversion*static_cast<float>(atof(lineTokens[lineTokenIndex++].c_str()));
        amoebaGpu->psWorkVector[1]->_pSysStream[0][offset++]       += conversion*static_cast<float>(atof(lineTokens[lineTokenIndex++].c_str()));
    }
    amoebaGpu->psWorkVector[0]->Upload();
    amoebaGpu->psWorkVector[1]->Upload();
}
#endif

        // post matrix multiply

        kSorUpdateMutualInducedField_kernel<<< numBlocks, numThreads >>>(
           gpu->natoms, amoebaGpu->psPolarizability->_pDevStream[0],
           amoebaGpu->psInducedDipole->_pDevStream[0], amoebaGpu->psInducedDipolePolar->_pDevStream[0],
           amoebaGpu->psE_Field->_pDevStream[0],       amoebaGpu->psE_FieldPolar->_pDevStream[0],
           amoebaGpu->psWorkVector[0]->_pDevStream[0], amoebaGpu->psWorkVector[1]->_pDevStream[0] );
        LAUNCHERROR("kSorUpdatePmeMutualInducedField");  

        // get total epsilon -- performing sums on gpu

        kReduceMutualInducedFieldDelta_kernel<<<1, amoebaGpu->epsilonThreadsPerBlock, 2*sizeof(float)*amoebaGpu->epsilonThreadsPerBlock>>>(
           3*gpu->natoms, amoebaGpu->psWorkVector[0]->_pDevStream[0], amoebaGpu->psWorkVector[1]->_pDevStream[0],
           amoebaGpu->psCurrentEpsilon->_pDevStream[0] );
        LAUNCHERROR("kReducePmeMutualInducedFieldDelta");

        if( amoebaGpu->log ){
            trackMutualInducedIterations( amoebaGpu, iteration);
        }

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

    if( 0 ){
        std::vector<int> fileId;
        //fileId.push_back( 0 );
        VectorOfDoubleVectors outputVector;
        cudaLoadCudaFloat4Array( gpu->natoms, 3, gpu->psPosq4,                    outputVector );
        cudaLoadCudaFloatArray( gpu->natoms,  3, amoebaGpu->psInducedDipole,      outputVector );
        cudaLoadCudaFloatArray( gpu->natoms,  3, amoebaGpu->psInducedDipolePolar, outputVector );
        cudaWriteVectorOfDoubleVectorsToFile( "CudaPmeMI", fileId, outputVector );
     }

   // ---------------------------------------------------------------------------------------
}

void cudaComputeAmoebaPmeMutualInducedField( amoebaGpuContext amoebaGpu )
{
    if( amoebaGpu->mutualInducedIterativeMethod == 0 ){
        cudaComputeAmoebaPmeMutualInducedFieldBySOR( amoebaGpu );
    }
}
