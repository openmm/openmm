
//-----------------------------------------------------------------------------------------

//-----------------------------------------------------------------------------------------

#include "amoebaCudaKernels.h"
#include "kCalculateAmoebaCudaUtilities.h"

//#define AMOEBA_DEBUG

static __constant__ cudaGmxSimulation cSim;
static __constant__ cudaAmoebaGmxSimulation cAmoebaSim;

void SetCalculateAmoebaCudaFixedEFieldSim(amoebaGpuContext amoebaGpu)
{
    cudaError_t status;
    gpuContext gpu = amoebaGpu->gpuContext;
    status         = cudaMemcpyToSymbol(cSim, &gpu->sim, sizeof(cudaGmxSimulation));    
    RTERROR(status, "SetCalculateAmoebaCudaFixedEFieldSim: cudaMemcpyToSymbol: SetSim copy to cSim failed");
    status         = cudaMemcpyToSymbol(cAmoebaSim, &amoebaGpu->amoebaSim, sizeof(cudaAmoebaGmxSimulation));    
    RTERROR(status, "SetCalculateAmoebaCudaFixedEFieldSim: cudaMemcpyToSymbol: SetSim copy to cAmoebaSim failed");
}

void GetCalculateAmoebaCudaFixedEFieldSim(amoebaGpuContext amoebaGpu)
{
    cudaError_t status;
    gpuContext gpu = amoebaGpu->gpuContext;
    status = cudaMemcpyFromSymbol(&gpu->sim, cSim, sizeof(cudaGmxSimulation));    
    RTERROR(status, "GetCalculateAmoebaCudaFixedEFieldSim: cudaMemcpyFromSymbol: SetSim copy from cSim failed");
    status = cudaMemcpyFromSymbol(&amoebaGpu->amoebaSim, cAmoebaSim, sizeof(cudaAmoebaGmxSimulation));         
    RTERROR(status, "GetCalculateAmoebaCudaFixedEFieldSim: cudaMemcpyFromSymbol: SetSim copy from cAmoebaSim failed");
}

// reduce psWorkArray_3_1 -> EField
// reduce psWorkArray_3_2 -> EFieldPolar

static void kReduceE_Fields_kernel(amoebaGpuContext amoebaGpu )
{
    gpuContext gpu = amoebaGpu->gpuContext;
    kReduceFields_kernel<<<gpu->sim.nonbond_blocks, gpu->sim.bsf_reduce_threads_per_block>>>(
                               gpu->sim.paddedNumberOfAtoms*3, gpu->sim.outputBuffers,
                               amoebaGpu->psWorkArray_3_1->_pDevData, amoebaGpu->psE_Field->_pDevData, 0 );
    LAUNCHERROR("kReduceE_Fields1");

    kReduceFields_kernel<<<gpu->sim.nonbond_blocks, gpu->sim.bsf_reduce_threads_per_block>>>(
                               gpu->sim.paddedNumberOfAtoms*3, gpu->sim.outputBuffers,
                               amoebaGpu->psWorkArray_3_2->_pDevData, amoebaGpu->psE_FieldPolar->_pDevData, 0 );
    LAUNCHERROR("kReduceE_Fields2");
}

// file includes FixedFieldParticle struct definition/load/unload struct and body kernel for fixed E-field

#undef GK
#include "kCalculateAmoebaCudaFixedFieldParticle.h"

// Include versions of the kernels for N^2 calculations.

#define METHOD_NAME(a, b) a##N2##b
#include "kCalculateAmoebaCudaFixedEField.h"
#define USE_OUTPUT_BUFFER_PER_WARP
#undef METHOD_NAME
#define METHOD_NAME(a, b) a##N2ByWarp##b
#include "kCalculateAmoebaCudaFixedEField.h"

/**---------------------------------------------------------------------------------------

   Compute fixed electric field

   @param amoebaGpu        amoebaGpu context
   @param gpu              OpenMM gpu Cuda context

   --------------------------------------------------------------------------------------- */

void cudaComputeAmoebaFixedEField( amoebaGpuContext amoebaGpu )
{
  
    gpuContext gpu    = amoebaGpu->gpuContext;

#ifdef AMOEBA_DEBUG
    static const char* methodName = "computeCudaAmoebaFixedEField";
    if( amoebaGpu->log ){
        (void) fprintf( amoebaGpu->log, "\n%s\n", methodName ); (void) fflush( amoebaGpu->log );
    }
    int paddedNumberOfAtoms                    = amoebaGpu->gpuContext->sim.paddedNumberOfAtoms;

    // N2 debug array

    CUDAStream<float4>* debugArray             = new CUDAStream<float4>(10*paddedNumberOfAtoms, 1, "DebugArray");
    memset( debugArray->_pSysData,      0, sizeof( float )*4*paddedNumberOfAtoms*paddedNumberOfAtoms);
    debugArray->Upload();

    // print intermediate results for the targetAtom 

    unsigned int targetAtom  = 3;
#endif

    kClearFields_3( amoebaGpu, 2 );

    static unsigned int threadsPerBlock = 0;
    if( threadsPerBlock == 0 ){ 
        unsigned int maxThreads;
        if (gpu->sm_version >= SM_20)
            maxThreads = 512; 
        else if (gpu->sm_version >= SM_12)
            maxThreads = 128; 
        else 
            maxThreads = 64;
        threadsPerBlock = std::min(getThreadsPerBlock(amoebaGpu, sizeof(FixedFieldParticle), gpu->sharedMemoryPerBlock ), maxThreads);
    }

#ifdef AMOEBA_DEBUG
    if( amoebaGpu->log ){
        (void) fprintf( amoebaGpu->log, "cudaComputeAmoebaFixedEField numBlocks=%u numThreads=%u bufferPerWarp=%u atm=%u shrd=%lu ixnCt=%lu workUnits=%lu\n",
                        gpu->sim.nonbond_blocks, threadsPerBlock, gpu->bOutputBufferPerWarp,
                        sizeof(FixedFieldParticle), sizeof(FixedFieldParticle)*threadsPerBlock, (*gpu->psInteractionCount)[0], gpu->sim.workUnits );
        (void) fflush( amoebaGpu->log );
    }
#endif

    if (gpu->bOutputBufferPerWarp){
        kCalculateAmoebaFixedE_FieldN2ByWarpForces_kernel<<<gpu->sim.nonbond_blocks, threadsPerBlock, sizeof(FixedFieldParticle)*threadsPerBlock>>>(
                                                                           gpu->psWorkUnit->_pDevData,
                                                                           amoebaGpu->psWorkArray_3_1->_pDevData,
                                                                           amoebaGpu->psWorkArray_3_2->_pDevData );
    } else {

        kCalculateAmoebaFixedE_FieldN2Forces_kernel<<<gpu->sim.nonbond_blocks, threadsPerBlock, sizeof(FixedFieldParticle)*threadsPerBlock>>>(
                                                                           gpu->psWorkUnit->_pDevData,
                                                                           amoebaGpu->psWorkArray_3_1->_pDevData,
                                                                           amoebaGpu->psWorkArray_3_2->_pDevData );
    }

    LAUNCHERROR("kCalculateAmoebaFixedE_FieldN2Forces_kernel");

#if 0
        for( unsigned int ii = 0; ii < gpu->sim.outputBuffers; ii++ ){
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

    kReduceE_Fields_kernel( amoebaGpu );

#ifdef AMOEBA_DEBUG
    if( amoebaGpu->log ){
        gpu->psInteractionCount->Download();
        (void) fprintf( amoebaGpu->log, "AmoebaN2Forces_kernel numBlocks=%u numThreads=%u bufferPerWarp=%u atm=%u shrd=%lu ixnCt=%lu workUnits=%lu\n",
                        gpu->sim.nonbond_blocks, threadsPerBlock, gpu->bOutputBufferPerWarp,
                        sizeof(FixedFieldParticle), sizeof(FixedFieldParticle)*threadsPerBlock, (*gpu->psInteractionCount)[0], gpu->sim.workUnits );
        (void) fflush( amoebaGpu->log );
        amoebaGpu->psWorkArray_3_1->Download();
        amoebaGpu->psWorkArray_3_2->Download();
        amoebaGpu->psE_Field->Download();
        amoebaGpu->psE_FieldPolar->Download();
        (void) fprintf( amoebaGpu->log, "OutEFields\n" );
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

           (void) fprintf( amoebaGpu->log,"\n" );
           if( ii == maxPrint && (gpu->natoms - maxPrint) > ii ){
                ii = gpu->natoms - maxPrint;
           }
        }
        (void) fflush( amoebaGpu->log );

        (void) fprintf( amoebaGpu->log, "EFields End\n" );

        // write results to file

        if( 0 ){
            std::vector<int> fileId;
            //fileId.push_back( 0 );
            VectorOfDoubleVectors outputVector;
            //cudaLoadCudaFloat4Array( gpu->natoms, 3, gpu->psPosq4,              outputVector, NULL, 1.0f );
            cudaLoadCudaFloatArray( gpu->natoms,  3, amoebaGpu->psE_Field,      outputVector, NULL, 1.0f );
            cudaLoadCudaFloatArray( gpu->natoms,  3, amoebaGpu->psE_FieldPolar, outputVector, NULL, 1.0f);
            cudaWriteVectorOfDoubleVectorsToFile( "CudaEField", fileId, outputVector );

         }
         //delete debugArray;
    }
#endif

}
