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

        float bn0         = erfc(ralpha)/r;
        float alsq2       = 2.0f*cSim.alphaEwald*cSim.alphaEwald;
        float alsq2n      = 1.0f/(cAmoebaSim.sqrtPi*cSim.alphaEwald);
        float exp2a       = exp(-(ralpha*ralpha));
        alsq2n           *= alsq2;
        float bn1         = (bn0+alsq2n*exp2a)/r2;

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
        float rr5         = 3.0f*(1.0f-dsc5)/r5;

        float preFactor1  = rr3 - bn1;
        float preFactor2  = bn2 - rr5;

        float dukr        = atomJ.inducedDipole[0]*xr      + atomJ.inducedDipole[1]*yr      + atomJ.inducedDipole[2]*zr;
        float preFactor3  = preFactor2*dukr;

        fields[0].x       = preFactor3*xr + preFactor1*atomJ.inducedDipole[0];
        fields[1].x       = preFactor3*yr + preFactor1*atomJ.inducedDipole[1];
        fields[2].x       = preFactor3*zr + preFactor1*atomJ.inducedDipole[2];


        float duir        = atomI.inducedDipole[0]*xr      + atomI.inducedDipole[1]*yr      + atomI.inducedDipole[2]*zr;
        preFactor3        = preFactor2*duir;

        fields[0].y       = preFactor3*xr + preFactor1*atomI.inducedDipole[0];
        fields[1].y       = preFactor3*yr + preFactor1*atomI.inducedDipole[1];
        fields[2].y       = preFactor3*zr + preFactor1*atomI.inducedDipole[2];


        float pukr        = atomJ.inducedDipolePolar[0]*xr + atomJ.inducedDipolePolar[1]*yr + atomJ.inducedDipolePolar[2]*zr;
        preFactor3        = preFactor2*pukr;

        fields[0].z       = preFactor3*xr + preFactor1*atomJ.inducedDipolePolar[0];
        fields[1].z       = preFactor3*yr + preFactor1*atomJ.inducedDipolePolar[1];
        fields[2].z       = preFactor3*zr + preFactor1*atomJ.inducedDipolePolar[2];


        float puir        = atomI.inducedDipolePolar[0]*xr + atomI.inducedDipolePolar[1]*yr + atomI.inducedDipolePolar[2]*zr;
        preFactor3        = preFactor2*puir;
        fields[0].w       = preFactor3*xr + preFactor1*atomI.inducedDipolePolar[0];
        fields[1].w       = preFactor3*yr + preFactor1*atomI.inducedDipolePolar[1];
        fields[2].w       = preFactor3*zr + preFactor1*atomI.inducedDipolePolar[2];

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
                   float* polarizability )
{

    int pos = blockIdx.x*blockDim.x + threadIdx.x;
    while( pos < 3*cSim.atoms )
    {   
        fixedEField[pos]         *= polarizability[pos];
        fixedEFieldPolar[pos]    *= polarizability[pos];

        pos                      += blockDim.x*gridDim.x;
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

    int pos = blockIdx.x*blockDim.x + threadIdx.x;
    while( pos < 3*cSim.atoms )
    {   

        float previousDipole           = inducedDipole[pos];
        float previousDipoleP          = inducedDipoleP[pos];
    
        // add self terms to fields
    
        const float term               = (4.0f/3.0f)*(cSim.alphaEwald*cSim.alphaEwald*cSim.alphaEwald)/cAmoebaSim.sqrtPi;
        matrixProduct[pos]            +=  term*previousDipole;
        matrixProductP[pos]           +=  term*previousDipoleP;
    
        inducedDipole[pos]             = fixedEField[pos]     + polarizability[pos]*matrixProduct[pos];
        inducedDipoleP[pos]            = fixedEFieldP[pos]    + polarizability[pos]*matrixProductP[pos];
    
        const float polarSOR           = 0.70f;
        inducedDipole[pos]             = previousDipole   + polarSOR*( inducedDipole[pos]   - previousDipole  );   
        inducedDipoleP[pos]            = previousDipoleP  + polarSOR*( inducedDipoleP[pos]  - previousDipoleP );
    
        matrixProduct[pos]             = ( inducedDipole[pos]  - previousDipole  )*( inducedDipole[pos]  - previousDipole  );
        matrixProductP[pos]            = ( inducedDipoleP[pos] - previousDipoleP )*( inducedDipoleP[pos] - previousDipoleP );

        pos                           += blockDim.x*gridDim.x;
    }

}

// reduce psWorkArray_3_1
// reduce psWorkArray_3_2

static void kReduceMutualInducedFields(amoebaGpuContext amoebaGpu, CUDAStream<float>* outputArray, CUDAStream<float>* outputPolarArray )
{
    gpuContext gpu = amoebaGpu->gpuContext;
    kReduceFields_kernel<<<gpu->sim.nonbond_blocks, gpu->sim.bsf_reduce_threads_per_block>>>(
                               gpu->sim.paddedNumberOfAtoms*3, gpu->sim.outputBuffers,
                               amoebaGpu->psWorkArray_3_1->_pDevData, outputArray->_pDevData, 0 );
    LAUNCHERROR("kReducePmeMI_Fields1");

    kReduceFields_kernel<<<gpu->sim.nonbond_blocks, gpu->sim.bsf_reduce_threads_per_block>>>(
                               gpu->sim.paddedNumberOfAtoms*3, gpu->sim.outputBuffers,
                               amoebaGpu->psWorkArray_3_2->_pDevData, outputPolarArray->_pDevData, 0 );
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
    int maxSlots                               = 10;
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

#ifdef AMOEBA_DEBUG
    if( amoebaGpu->log ){
        (void) fprintf( amoebaGpu->log, "cudaComputeAmoebaPmeMutualInducedFieldMatrixMultiply: numBlocks=%u numThreads=%u bufferPerWarp=%u atm=%lu shrd=%lu ixnCt=%lu workUnits=%u\n",
                        gpu->sim.nonbond_blocks, threadsPerBlock, gpu->bOutputBufferPerWarp,
                        sizeof(MutualInducedParticle), sizeof(MutualInducedParticle)*threadsPerBlock,
                        (*gpu->psInteractionCount)[0], gpu->sim.workUnits );
        (void) fflush( amoebaGpu->log );
    }
#endif

    if (gpu->bOutputBufferPerWarp){

        kCalculateAmoebaPmeMutualInducedFieldCutoffByWarp_kernel<<<gpu->sim.nonbond_blocks, threadsPerBlock, sizeof(MutualInducedParticle)*threadsPerBlock>>>(
                                                                 gpu->sim.pInteractingWorkUnit,
                                                                 amoebaGpu->psWorkArray_3_1->_pDevData,
#ifdef AMOEBA_DEBUG
                                                                 amoebaGpu->psWorkArray_3_2->_pDevData,
                                                                 debugArray->_pDevData, targetAtom );
#else
                                                                 amoebaGpu->psWorkArray_3_2->_pDevData );
#endif

    } else {

        kCalculateAmoebaPmeMutualInducedFieldCutoff_kernel<<<gpu->sim.nonbond_blocks, threadsPerBlock, sizeof(MutualInducedParticle)*threadsPerBlock>>>(
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

   // ---------------------------------------------------------------------------------------

    // set  E_Field & E_FieldPolar] to [ E_Field & E_FieldPolar]*Polarizability
    // initialize [ InducedDipole & InducedDipolePolar ] to [ E_Field & E_FieldPolar]*Polarizability

    kInitializeMutualInducedField_kernel<<< gpu->sim.nonbond_blocks, gpu->sim.bsf_reduce_threads_per_block >>>(
         gpu->natoms,
         amoebaGpu->psE_Field->_pDevData,
         amoebaGpu->psE_FieldPolar->_pDevData,
         amoebaGpu->psPolarizability->_pDevData );
    LAUNCHERROR("AmoebaPmeMutualInducedFieldSetup");  

    cudaMemcpy( amoebaGpu->psInducedDipole->_pDevData,        amoebaGpu->psE_Field->_pDevData,       3*gpu->sim.paddedNumberOfAtoms*sizeof( float ), cudaMemcpyDeviceToDevice );
    cudaMemcpy( amoebaGpu->psInducedDipolePolar->_pDevData,   amoebaGpu->psE_FieldPolar->_pDevData,  3*gpu->sim.paddedNumberOfAtoms*sizeof( float ), cudaMemcpyDeviceToDevice );

#ifdef AMOEBA_DEBUG
    if( amoebaGpu->log ){

        std::vector<int> fileId;
        VectorOfDoubleVectors outputVector;
        cudaLoadCudaFloatArray( gpu->natoms,  3, amoebaGpu->psE_Field,            outputVector, gpu->psAtomIndex->_pSysData, 1.0f );
        cudaLoadCudaFloatArray( gpu->natoms,  3, amoebaGpu->psE_FieldPolar,       outputVector, gpu->psAtomIndex->_pSysData, 1.0f );
        cudaLoadCudaFloatArray( gpu->natoms,  3, amoebaGpu->psInducedDipole,      outputVector, gpu->psAtomIndex->_pSysData, 1.0f );
        cudaLoadCudaFloatArray( gpu->natoms,  3, amoebaGpu->psInducedDipolePolar, outputVector, gpu->psAtomIndex->_pSysData, 1.0f );
        cudaWriteVectorOfDoubleVectorsToFile( "CudaPmeEFieldPolarity", fileId, outputVector );
    }   
#endif

    // if polarization type is direct, set flags signalling done and return

    if( amoebaGpu->amoebaSim.polarizationType )
    {
        amoebaGpu->mutualInducedDone          = 1;
        amoebaGpu->mutualInducedConverged     = 1;
        kCalculateAmoebaPMEInducedDipoleField( amoebaGpu );
        return;
    }

    // ---------------------------------------------------------------------------------------
 
    done      = 0;
    iteration = 1;

    while( !done ){

        // matrix multiply

        cudaComputeAmoebaPmeMutualInducedFieldMatrixMultiply( amoebaGpu, amoebaGpu->psWorkVector[0],  amoebaGpu->psWorkVector[1] );
        kCalculateAmoebaPMEInducedDipoleField( amoebaGpu );

        // post matrix multiply

        kSorUpdateMutualInducedField_kernel<<< gpu->sim.nonbond_blocks, gpu->sim.bsf_reduce_threads_per_block >>>(
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

        if( 0 && amoebaGpu->log ){ // trackMutualInducedIterations
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

        if( 0 && amoebaGpu->mutualInducedCurrentEpsilon != amoebaGpu->mutualInducedCurrentEpsilon ){
            (void) fprintf( stderr, "PME MI iteration=%3d eps is nan -- exiting.\n", iteration );
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

    if( 0 ){
        static int iteration = 0;
        checkForNans( gpu->natoms,  3, amoebaGpu->psInducedDipole, gpu->psAtomIndex->_pSysData,    ++iteration, "CudaPmeMI", stderr );
        checkForNans( gpu->natoms,  3, amoebaGpu->psInducedDipolePolar, gpu->psAtomIndex->_pSysData, iteration, "CudaPmeMIPolar", stderr );
     }

   // ---------------------------------------------------------------------------------------
}

void cudaComputeAmoebaPmeMutualInducedField( amoebaGpuContext amoebaGpu )
{
    if( amoebaGpu->mutualInducedIterativeMethod == 0 ){
        cudaComputeAmoebaPmeMutualInducedFieldBySOR( amoebaGpu );
    }
}
