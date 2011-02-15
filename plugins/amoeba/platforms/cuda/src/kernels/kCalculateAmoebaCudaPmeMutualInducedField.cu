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

#undef INCLUDE_MI_FIELD_BUFFERS
#define INCLUDE_MI_FIELD_BUFFERS 
#include "kCalculateAmoebaCudaMutualInducedParticle.h"
#undef INCLUDE_MI_FIELD_BUFFERS

__device__ void sumTempBuffer( MutualInducedParticle& atomI, MutualInducedParticle& atomJ ){

    atomI.tempBuffer[0]  += atomJ.tempBuffer[0];
    atomI.tempBuffer[1]  += atomJ.tempBuffer[1];
    atomI.tempBuffer[2]  += atomJ.tempBuffer[2];

    atomI.tempBufferP[0] += atomJ.tempBufferP[0];
    atomI.tempBufferP[1] += atomJ.tempBufferP[1];
    atomI.tempBufferP[2] += atomJ.tempBufferP[2];
}

// file includes FixedFieldParticle struct definition/load/unload struct and body kernel for fixed E-field

__device__ void calculatePmeDirectMutualInducedFieldPairIxn_kernel( MutualInducedParticle& atomI, MutualInducedParticle& atomJ,
                                                                    float uscale, float4 fields[3]
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
    if( r2 <= cSim.nonbondedCutoffSqr ){
        float r           = sqrtf(r2);

        // calculate the error function damping terms

        float ralpha      = cSim.alphaEwald*r;

        float bn0             = erfc(ralpha)/r;
        float alsq2       = 2.0f*cSim.alphaEwald*cSim.alphaEwald;
        float alsq2n      = 1.0f/(cAmoebaSim.sqrtPi*cSim.alphaEwald);
        float exp2a       = exp(-(ralpha*ralpha));
        alsq2n           *= alsq2;
        float bn1             = (bn0+alsq2n*exp2a)/r2;

        alsq2n           *= alsq2;
        float bn2         = (3.0f*bn1+alsq2n*exp2a)/r2;

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

        bn1              *= -1.0f;

        float fimd0       = bn1*atomJ.inducedDipole[0]      + bn2*dukr*xr;
        float fimd1       = bn1*atomJ.inducedDipole[1]      + bn2*dukr*yr;
        float fimd2       = bn1*atomJ.inducedDipole[2]      + bn2*dukr*zr;

        float fkmd0       = bn1*atomI.inducedDipole[0]      + bn2*duir*xr;
        float fkmd1       = bn1*atomI.inducedDipole[1]      + bn2*duir*yr;
        float fkmd2       = bn1*atomI.inducedDipole[2]      + bn2*duir*zr;

        float fimp0       = bn1*atomJ.inducedDipolePolar[0] + bn2*pukr*xr;
        float fimp1       = bn1*atomJ.inducedDipolePolar[1] + bn2*pukr*yr;
        float fimp2       = bn1*atomJ.inducedDipolePolar[2] + bn2*pukr*zr;

        float fkmp0       = bn1*atomI.inducedDipolePolar[0] + bn2*puir*xr;
        float fkmp1       = bn1*atomI.inducedDipolePolar[1] + bn2*puir*yr;
        float fkmp2       = bn1*atomI.inducedDipolePolar[2] + bn2*puir*zr;

        rr3              *= -1.0f;
        float fid0        = rr3*atomJ.inducedDipole[0]      + rr5*dukr*xr;
        float fid1        = rr3*atomJ.inducedDipole[1]      + rr5*dukr*yr;
        float fid2        = rr3*atomJ.inducedDipole[2]      + rr5*dukr*zr;

        float fkd0        = rr3*atomI.inducedDipole[0]      + rr5*duir*xr;
        float fkd1        = rr3*atomI.inducedDipole[1]      + rr5*duir*yr;
        float fkd2        = rr3*atomI.inducedDipole[2]      + rr5*duir*zr;

        float fip0        = rr3*atomJ.inducedDipolePolar[0] + rr5*pukr*xr;
        float fip1        = rr3*atomJ.inducedDipolePolar[1] + rr5*pukr*yr;
        float fip2        = rr3*atomJ.inducedDipolePolar[2] + rr5*pukr*zr;

        float fkp0        = rr3*atomI.inducedDipolePolar[0] + rr5*puir*xr;
        float fkp1        = rr3*atomI.inducedDipolePolar[1] + rr5*puir*yr;
        float fkp2        = rr3*atomI.inducedDipolePolar[2] + rr5*puir*zr;

        // increment the field at each site due to this interaction

        fields[0].x       = fimd0 - fid0;
        fields[0].y       = fkmd0 - fkd0;
        fields[0].z       = fimp0 - fip0;
        fields[0].w       = fkmp0 - fkp0;
    
        fields[1].x       = fimd1 - fid1;
        fields[1].y       = fkmd1 - fkd1;
        fields[1].z       = fimp1 - fip1;
        fields[1].w       = fkmp1 - fkp1;
    
        fields[2].x       = fimd2 - fid2;
        fields[2].y       = fkmd2 - fkd2;
        fields[2].z       = fimp2 - fip2;
        fields[2].w       = fkmp2 - fkp2;
 
    } else {

        fields[0].x       = 0.0f;
        fields[0].y       = 0.0f;
        fields[0].z       = 0.0f;
        fields[0].w       = 0.0f;
    
        fields[1].x       = 0.0f;
        fields[1].y       = 0.0f;
        fields[1].z       = 0.0f;
        fields[1].w       = 0.0f;
    
        fields[2].x       = 0.0f;
        fields[2].y       = 0.0f;
        fields[2].z       = 0.0f;
        fields[2].w       = 0.0f;
    }
/*
#ifdef AMOEBA_DEBUG
    pullBack[0].x = xr;
    pullBack[0].y = yr;
    pullBack[0].z = zr;
    pullBack[0].w = r2;

    pullBack[1].x = alsq2;
    pullBack[1].y = bn0;
    pullBack[1].z = bn2;
    pullBack[1].w = exp2a;

    pullBack[1].x = atomJ.x - atomI.x;
    pullBack[1].y = atomJ.y - atomI.y;
    pullBack[1].z = atomJ.z - atomI.z;
    pullBack[1].w = (atomJ.x - atomI.x)*(atomJ.x - atomI.x) + (atomJ.y - atomI.y)*(atomJ.y - atomI.y)+ (atomJ.z - atomI.z)*(atomJ.z - atomI.z);
    pullBack[1].x = scale3;
    pullBack[1].y = scale5;
    pullBack[1].z = scale7;
#endif
*/
}

// Include versions of the kernels for N^2 calculations.

#define METHOD_NAME(a, b) a##Cutoff##b
#include "kCalculateAmoebaCudaPmeMutualInducedField.h"
#define USE_OUTPUT_BUFFER_PER_WARP
#undef METHOD_NAME
#define METHOD_NAME(a, b) a##CutoffByWarp##b
#include "kCalculateAmoebaCudaPmeMutualInducedField.h"

__global__ 
#if (__CUDA_ARCH__ >= 200)
__launch_bounds__(GF1XX_THREADS_PER_BLOCK, 1)
#elif (__CUDA_ARCH__ >= 120)
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
#elif (__CUDA_ARCH__ >= 120)
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
        epsilon[0]  = 48.033324f*sqrtf( epsilon[0]/( (float) (numberOfEntries/3)) );
#ifdef AMOEBA_DEBUG
        epsilon[1]  = 48.033324f*sqrtf( delta[0].x/( (float) (numberOfEntries/3)) );
        epsilon[2]  = 48.033324f*sqrtf( delta[0].y/( (float) (numberOfEntries/3)) );
#endif
    }   
}

/**

   matrixProduct/matrixProductP contains epsilon**2 on output

*/
__global__ 
#if (__CUDA_ARCH__ >= 200)
__launch_bounds__(GF1XX_THREADS_PER_BLOCK, 1)
#elif (__CUDA_ARCH__ >= 120)
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

    int threadId                        = __mul24(blockIdx.x,blockDim.x) + threadIdx.x;
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
                               amoebaGpu->psWorkArray_3_1->_pDevData, outputArray->_pDevData );
    LAUNCHERROR("kReducePmeMI_Fields1");

    kReduceFields_kernel<<<amoebaGpu->nonbondBlocks, amoebaGpu->fieldReduceThreadsPerBlock>>>(
                               amoebaGpu->paddedNumberOfAtoms*3, amoebaGpu->outputBuffers,
                               amoebaGpu->psWorkArray_3_2->_pDevData, outputPolarArray->_pDevData );
    LAUNCHERROR("kReducePmeMI_Fields2");
}

/**---------------------------------------------------------------------------------------

   Compute mutual induce field

   @param amoebaGpu        amoebaGpu context

   --------------------------------------------------------------------------------------- */

static void cudaComputeAmoebaPmeMutualInducedFieldMatrixMultiply( amoebaGpuContext amoebaGpu,
                                                                  CUDAStream<float>* outputArray, CUDAStream<float>* outputPolarArray )
{
  
  static unsigned int threadsPerBlock  = 0;
  gpuContext gpu                       = amoebaGpu->gpuContext;

#ifdef AMOEBA_DEBUG
    int targetAtom                = 546;
    static const char* methodName = "cudaComputeAmoebaPmeMutualInducedFieldMatrixMultiply";
    static int iteration          = 1;
    if( 1 && amoebaGpu->log ){
        (void) fprintf( amoebaGpu->log, "%s\n", methodName );
        (void) fflush( amoebaGpu->log );
    }
    int paddedNumberOfAtoms                    = amoebaGpu->gpuContext->sim.paddedNumberOfAtoms;
    int maxSlots                               = amoebaGpu->gpuContext->sim.paddedNumberOfAtoms;
    CUDAStream<float4>* debugArray             = new CUDAStream<float4>(maxSlots*paddedNumberOfAtoms, 1, "DebugArray");
    memset( debugArray->_pSysData,      0, sizeof( float )*4*maxSlots*paddedNumberOfAtoms);
    debugArray->Upload();
#endif

    kClearFields_3( amoebaGpu, 2 );

    // on first pass, set threads/block

    if( threadsPerBlock == 0 ){  
        unsigned int maxThreads;
        if (gpu->sm_version >= SM_20)
            maxThreads = 384; 
        else if (gpu->sm_version >= SM_12)
            maxThreads = 128; 
        else
            maxThreads = 64; 
        threadsPerBlock = std::min(getThreadsPerBlock(amoebaGpu, sizeof(MutualInducedParticle)), maxThreads);
    }    

    if (gpu->bOutputBufferPerWarp){

#ifdef AMOEBA_DEBUG
        (void) fprintf( amoebaGpu->log, "Cutoff -- use warp\n" );
        (void) fprintf( amoebaGpu->log, "AmoebaCutoffForces_kernel numBlocks=%u numThreads=%u bufferPerWarp=%u atm=%u shrd=%u Ebuf=%u ixnCt=%u workUnits=%u\n",
                        amoebaGpu->nonbondBlocks, threadsPerBlock, amoebaGpu->bOutputBufferPerWarp,
                        sizeof(MutualInducedParticle), sizeof(MutualInducedParticle)*threadsPerBlock,
                        amoebaGpu->energyOutputBuffers, (*gpu->psInteractionCount)[0], gpu->sim.workUnits );
        (void) fflush( amoebaGpu->log );
#endif
                                                                 //gpu->sim.pInteractingWorkUnit,
                                                                 //amoebaGpu->psWorkUnit->_pDevData,
        kCalculateAmoebaPmeMutualInducedFieldCutoffByWarp_kernel<<<amoebaGpu->nonbondBlocks, threadsPerBlock, sizeof(MutualInducedParticle)*threadsPerBlock>>>(
                                                                 gpu->sim.pInteractingWorkUnit,
                                                                 amoebaGpu->psWorkArray_3_1->_pDevData,
#ifdef AMOEBA_DEBUG
                                                                 amoebaGpu->psWorkArray_3_2->_pDevData,
                                                                 debugArray->_pDevData, targetAtom );
#else
                                                                 amoebaGpu->psWorkArray_3_2->_pDevData );
#endif

    } else {

#ifdef AMOEBA_DEBUG
        (void) fprintf( amoebaGpu->log, "Cutoff no warp\n" );
        (void) fprintf( amoebaGpu->log, "AmoebaCutoffForces_kernel numBlocks=%u numThreads=%u bufferPerWarp=%u atm=%u shrd=%u Ebuf=%u ixnCt=%u workUnits=%u\n",
                        amoebaGpu->nonbondBlocks, threadsPerBlock, amoebaGpu->bOutputBufferPerWarp,
                        sizeof(MutualInducedParticle), sizeof(MutualInducedParticle)*threadsPerBlock,
                        amoebaGpu->energyOutputBuffers, (*gpu->psInteractionCount)[0], gpu->sim.workUnits );
        (void) fflush( amoebaGpu->log );
#endif
        kCalculateAmoebaPmeMutualInducedFieldCutoff_kernel<<<amoebaGpu->nonbondBlocks, threadsPerBlock, sizeof(MutualInducedParticle)*threadsPerBlock>>>(
                                                                 gpu->sim.pInteractingWorkUnit,
                                                                 amoebaGpu->psWorkArray_3_1->_pDevData,
#ifdef AMOEBA_DEBUG
                                                                 amoebaGpu->psWorkArray_3_2->_pDevData,
                                                                 debugArray->_pDevData, targetAtom );
#else
                                                                 amoebaGpu->psWorkArray_3_2->_pDevData );
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
                            outputArray->_pSysData[indexOffset],
                            outputArray->_pSysData[indexOffset+1],
                            outputArray->_pSysData[indexOffset+2] );
     
            // MI polar
 
            (void) fprintf( amoebaGpu->log,"MultP[%16.9e %16.9e %16.9e]\n",
                            outputPolarArray->_pSysData[indexOffset],
                            outputPolarArray->_pSysData[indexOffset+1],
                            outputPolarArray->_pSysData[indexOffset+2] );
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
                                debugArray->_pSysData[debugIndex].x, debugArray->_pSysData[debugIndex].y,
                                debugArray->_pSysData[debugIndex].z, debugArray->_pSysData[debugIndex].w );
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

static void cudaComputeAmoebaPmeMutualInducedFieldBySOR( amoebaGpuContext amoebaGpu )
{
  
   // ---------------------------------------------------------------------------------------

//#define AMOEBA_DEBUG
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

   // ---------------------------------------------------------------------------------------

    // set  E_Field & E_FieldPolar] to [ E_Field & E_FieldPolar]*Polarizability
    // initialize [ InducedDipole & InducedDipolePolar ] to [ E_Field & E_FieldPolar]*Polarizability

    kInitializeMutualInducedField_kernel<<< numBlocks, numThreads >>>(
         gpu->natoms,
         amoebaGpu->psE_Field->_pDevData,
         amoebaGpu->psE_FieldPolar->_pDevData,
         amoebaGpu->psPolarizability->_pDevData,
         amoebaGpu->psInducedDipole->_pDevData,
         amoebaGpu->psInducedDipolePolar->_pDevData );
    LAUNCHERROR("AmoebaPmeMutualInducedFieldSetup");  

#ifdef AMOEBA_DEBUG
    if( amoebaGpu->log ){

        gpuContext gpu = amoebaGpu->gpuContext;
        std::vector<int> fileId;
        VectorOfDoubleVectors outputVector;
        cudaLoadCudaFloatArray( gpu->natoms,  3, amoebaGpu->psE_Field,            outputVector, gpu->psAtomIndex->_pSysData, 1.0f );
        cudaLoadCudaFloatArray( gpu->natoms,  3, amoebaGpu->psE_FieldPolar,       outputVector, gpu->psAtomIndex->_pSysData, 1.0f );
        cudaLoadCudaFloatArray( gpu->natoms,  3, amoebaGpu->psInducedDipole,      outputVector, gpu->psAtomIndex->_pSysData, 1.0f );
        cudaLoadCudaFloatArray( gpu->natoms,  3, amoebaGpu->psInducedDipolePolar, outputVector, gpu->psAtomIndex->_pSysData, 1.0f );
        cudaWriteVectorOfDoubleVectorsToFile( "CudaEFieldPolarity", fileId, outputVector );
/*
        amoebaGpu->psE_FieldPolar->Download();
        amoebaGpu->psInducedDipole->Download(),
        amoebaGpu->psInducedDipolePolar->Download();
        amoebaGpu->psPolarizability->Download();
        (void) fprintf( amoebaGpu->log, "%s Initial setup for matrix multiply\n", methodName );
        int offset   = 0;
        int maxPrint = 10;
         cudaLoadCudaFloatArray( gpu->natoms,  3, amoebaGpu->psForce,      outputVector, gpu->psAtomIndex->_pSysData, 1.0f );
        for( int ii = 0; ii < gpu->natoms; ii++ ){
            (void) fprintf( amoebaGpu->log, "%4d pol=%12.4e ", ii, 
                            amoebaGpu->psPolarizability->_pSysData[offset] );
            if( amoebaGpu->psPolarizability->_pSysData[offset] != amoebaGpu->psPolarizability->_pSysData[offset+1] ||
                amoebaGpu->psPolarizability->_pSysData[offset] != amoebaGpu->psPolarizability->_pSysData[offset+2] ){
                (void) fprintf( amoebaGpu->log, "PolX!!! %12.4e %12.4e ", amoebaGpu->psPolarizability->_pSysData[offset+1], amoebaGpu->psPolarizability->_pSysData[offset+2] ); 
            }

            (void) fprintf( amoebaGpu->log," E[%14.6e %14.6e %14.6e] Mi[%14.6e %14.6e %14.6e] ",
                            amoebaGpu->psE_Field->_pSysData[offset],       amoebaGpu->psE_Field->_pSysData[offset+1],       amoebaGpu->psE_Field->_pSysData[offset+2],
                            amoebaGpu->psInducedDipole->_pSysData[offset], amoebaGpu->psInducedDipole->_pSysData[offset+1], amoebaGpu->psInducedDipole->_pSysData[offset+2] );
            (void) fprintf( amoebaGpu->log,"Ep[%14.6e %14.6e %14.6e] Mip[%14.6e %14.6e %14.6e]\n",
                            amoebaGpu->psE_FieldPolar->_pSysData[offset],       amoebaGpu->psE_FieldPolar->_pSysData[offset+1],       amoebaGpu->psE_FieldPolar->_pSysData[offset+2],
                            amoebaGpu->psInducedDipolePolar->_pSysData[offset], amoebaGpu->psInducedDipolePolar->_pSysData[offset+1], amoebaGpu->psInducedDipolePolar->_pSysData[offset+2] );
            offset += 3;
            if( ii == maxPrint && (ii < (gpu->natoms - maxPrint) ) )ii =  (gpu->natoms - maxPrint);
        }   
        
void) fflush( amoebaGpu->log );
*/
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
        amoebaGpu->psWorkVector[0]->_pSysData[offset++]       += conversion*static_cast<float>(atof(lineTokens[lineTokenIndex++].c_str()));
        amoebaGpu->psWorkVector[0]->_pSysData[offset++]       += conversion*static_cast<float>(atof(lineTokens[lineTokenIndex++].c_str()));
        amoebaGpu->psWorkVector[0]->_pSysData[offset++]       += conversion*static_cast<float>(atof(lineTokens[lineTokenIndex++].c_str())); 
        offset                                                     -= 3;        
        amoebaGpu->psWorkVector[1]->_pSysData[offset++]       += conversion*static_cast<float>(atof(lineTokens[lineTokenIndex++].c_str()));
        amoebaGpu->psWorkVector[1]->_pSysData[offset++]       += conversion*static_cast<float>(atof(lineTokens[lineTokenIndex++].c_str()));
        amoebaGpu->psWorkVector[1]->_pSysData[offset++]       += conversion*static_cast<float>(atof(lineTokens[lineTokenIndex++].c_str()));
    }
    amoebaGpu->psWorkVector[0]->Upload();
    amoebaGpu->psWorkVector[1]->Upload();
}
#endif

        // post matrix multiply

        kSorUpdateMutualInducedField_kernel<<< numBlocks, numThreads >>>(
           gpu->natoms, amoebaGpu->psPolarizability->_pDevData,
           amoebaGpu->psInducedDipole->_pDevData, amoebaGpu->psInducedDipolePolar->_pDevData,
           amoebaGpu->psE_Field->_pDevData,       amoebaGpu->psE_FieldPolar->_pDevData,
           amoebaGpu->psWorkVector[0]->_pDevData, amoebaGpu->psWorkVector[1]->_pDevData );
        LAUNCHERROR("kSorUpdatePmeMutualInducedField");  

            if( 0 ){
                gpuContext gpu = amoebaGpu->gpuContext;
                std::vector<int> fileId;
                fileId.push_back( iteration );
                VectorOfDoubleVectors outputVector;
//                cudaLoadCudaFloatArray( gpu->natoms,  3, amoebaGpu->psE_Field, outputVector, gpu->psAtomIndex->_pSysData );
//                cudaLoadCudaFloatArray( gpu->natoms,  3, amoebaGpu->psE_FieldPolar, outputVector, gpu->psAtomIndex->_pSysData );
                cudaLoadCudaFloatArray( gpu->natoms,  3, amoebaGpu->psInducedDipole, outputVector, gpu->psAtomIndex->_pSysData, 1.0f );
                cudaLoadCudaFloatArray( gpu->natoms,  3, amoebaGpu->psInducedDipolePolar, outputVector, gpu->psAtomIndex->_pSysData, 1.0f );
                cudaWriteVectorOfDoubleVectorsToFile( "CudaPmeDirectMI", fileId, outputVector );
            }

        // get total epsilon -- performing sums on gpu

        kReduceMutualInducedFieldDelta_kernel<<<1, amoebaGpu->epsilonThreadsPerBlock, 2*sizeof(float)*amoebaGpu->epsilonThreadsPerBlock>>>(
           3*gpu->natoms, amoebaGpu->psWorkVector[0]->_pDevData, amoebaGpu->psWorkVector[1]->_pDevData,
           amoebaGpu->psCurrentEpsilon->_pDevData );
        LAUNCHERROR("kReducePmeMutualInducedFieldDelta");

        if( amoebaGpu->log ){
            trackMutualInducedIterations( amoebaGpu, iteration);
        }

        // Debye=48.033324f
        amoebaGpu->psCurrentEpsilon->Download();
        float currentEpsilon                     = amoebaGpu->psCurrentEpsilon->_pSysData[0];
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
                           amoebaGpu->psCurrentEpsilon->_pSysData[1], 
                           amoebaGpu->psCurrentEpsilon->_pSysData[2], done );
#else
           (void) fprintf( amoebaGpu->log, "%s iteration=%3d eps %14.6e %14.6e crrntEps=%14.6e %14.6e %14.6e %14.6e done=%d\n",
                           methodName, iteration, sum1, sum2, amoebaGpu->mutualInducedCurrentEpsilon,
                           amoebaGpu->psCurrentEpsilon->_pSysData[0], 
                           amoebaGpu->psCurrentEpsilon->_pSysData[1], 
                           amoebaGpu->psCurrentEpsilon->_pSysData[2], done );
#endif
           (void) fflush( amoebaGpu->log );

            if( 0 ){
                gpuContext gpu = amoebaGpu->gpuContext;
                std::vector<int> fileId;
                fileId.push_back( iteration );
                VectorOfDoubleVectors outputVector;
                cudaLoadCudaFloatArray( gpu->natoms,  3, amoebaGpu->psE_Field, outputVector, gpu->psAtomIndex->_pSysData, 1.0f );
                cudaLoadCudaFloatArray( gpu->natoms,  3, amoebaGpu->psE_FieldPolar, outputVector, gpu->psAtomIndex->_pSysData, 1.0f );
                cudaLoadCudaFloatArray( gpu->natoms,  3, amoebaGpu->psInducedDipole, outputVector, gpu->psAtomIndex->_pSysData, 1.0f );
                cudaLoadCudaFloatArray( gpu->natoms,  3, amoebaGpu->psInducedDipolePolar, outputVector, gpu->psAtomIndex->_pSysData, 1.0f );
                cudaWriteVectorOfDoubleVectorsToFile( "CudaPmeMI", fileId, outputVector );
            }
/*
            int offset   = 0;
            int maxPrint = 10;
            for( int ii = 0; ii < gpu->natoms; ii++ ){
                (void) fprintf( amoebaGpu->log, "%4d ", ii ); 
    
                (void) fprintf( amoebaGpu->log," Mi[%14.6e %14.6e %14.6e] ",
                                amoebaGpu->psInducedDipole->_pSysData[offset],
                                amoebaGpu->psInducedDipole->_pSysData[offset+1],
                                amoebaGpu->psInducedDipole->_pSysData[offset+2] );
                (void) fprintf( amoebaGpu->log,"Mip[%14.6e %14.6e %14.6e]\n",
                                amoebaGpu->psInducedDipolePolar->_pSysData[offset],
                                amoebaGpu->psInducedDipolePolar->_pSysData[offset+1],
                                amoebaGpu->psInducedDipolePolar->_pSysData[offset+2] );
                if( ii == maxPrint && (ii < (gpu->natoms - maxPrint) ) ){
                    ii =  (gpu->natoms - maxPrint);
                    offset = 3*(ii+1);
                } else {
                    offset += 3;
                }
            }   
            (void) fflush( amoebaGpu->log );
*/

            if( 0 ){
                std::vector<int> fileId;
                fileId.push_back( iteration );
                VectorOfDoubleVectors outputVector;
                cudaLoadCudaFloat4Array( gpu->natoms, 3, gpu->psPosq4,                    outputVector, gpu->psAtomIndex->_pSysData, 1.0f );
                cudaLoadCudaFloatArray( gpu->natoms,  3, amoebaGpu->psInducedDipole,      outputVector, gpu->psAtomIndex->_pSysData, 1.0f );
                cudaLoadCudaFloatArray( gpu->natoms,  3, amoebaGpu->psInducedDipolePolar, outputVector, gpu->psAtomIndex->_pSysData, 1.0f );
                cudaWriteVectorOfDoubleVectorsToFile( "CudaPmeMI", fileId, outputVector );
            }

        }

        (void) fprintf( amoebaGpu->log, "MI iteration=%3d eps %14.6e [%14.6e %14.6e] done=%d\n",
                        iteration, amoebaGpu->mutualInducedCurrentEpsilon,
                        amoebaGpu->psCurrentEpsilon->_pSysData[1], 
                        amoebaGpu->psCurrentEpsilon->_pSysData[2], done );
        (void) fflush( amoebaGpu->log );
#endif
        // exit if nan
        if( amoebaGpu->mutualInducedCurrentEpsilon != amoebaGpu->mutualInducedCurrentEpsilon ){
            (void) fprintf( amoebaGpu->log, "PME MI iteration=%3d eps is nan -- exiting.\n", iteration );
            exit(0);
        }


        iteration++;
    }

    amoebaGpu->mutualInducedDone             = done;
    amoebaGpu->mutualInducedConverged        = ( !done || iteration > amoebaGpu->mutualInducedMaxIterations ) ? 0 : 1;

    if( 0 ){
        std::vector<int> fileId;
        //fileId.push_back( 0 );
        VectorOfDoubleVectors outputVector;
        //cudaLoadCudaFloat4Array( gpu->natoms, 3, gpu->psPosq4,                    outputVector, 1.0f );
        cudaLoadCudaFloatArray( gpu->natoms,  3, amoebaGpu->psInducedDipole,      outputVector, gpu->psAtomIndex->_pSysData, 1.0f );
        cudaLoadCudaFloatArray( gpu->natoms,  3, amoebaGpu->psInducedDipolePolar, outputVector, gpu->psAtomIndex->_pSysData, 1.0f );
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
