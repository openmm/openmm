
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
    kReduceFields_kernel<<<amoebaGpu->nonbondBlocks, amoebaGpu->fieldReduceThreadsPerBlock>>>(
                               amoebaGpu->paddedNumberOfAtoms*3, amoebaGpu->outputBuffers,
                               amoebaGpu->psWorkArray_3_1->_pDevData, amoebaGpu->psE_Field->_pDevData );
    LAUNCHERROR("kReduceE_Fields1");

    kReduceFields_kernel<<<amoebaGpu->nonbondBlocks, amoebaGpu->fieldReduceThreadsPerBlock>>>(
                               amoebaGpu->paddedNumberOfAtoms*3, amoebaGpu->outputBuffers,
                               amoebaGpu->psWorkArray_3_2->_pDevData, amoebaGpu->psE_FieldPolar->_pDevData );
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
  
   // ---------------------------------------------------------------------------------------
   // ---------------------------------------------------------------------------------------

    gpuContext gpu    = amoebaGpu->gpuContext;

#ifdef AMOEBA_DEBUG
    static const char* methodName = "computeCudaAmoebaFixedEField";
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

    unsigned int targetAtom  = 3;
#endif

    kClearFields_3( amoebaGpu, 2 );

    if (gpu->bOutputBufferPerWarp){
#ifdef AMOEBA_DEBUG
        (void) fprintf( amoebaGpu->log, "N2 warp\n" );
        (void) fprintf( amoebaGpu->log, "%s numBlocks=%u numThreads=%u bufferPerWarp=%u atm=%u shrd=%u ixnCt=%u workUnits=%u\n", methodName,
                        amoebaGpu->nonbondBlocks, amoebaGpu->nonbondThreadsPerBlock, amoebaGpu->bOutputBufferPerWarp,
                        sizeof(FixedFieldParticle), sizeof(FixedFieldParticle)*amoebaGpu->nonbondThreadsPerBlock, (*gpu->psInteractionCount)[0], gpu->sim.workUnits );
        (void) fflush( amoebaGpu->log );
#endif
        kCalculateAmoebaFixedE_FieldN2ByWarpForces_kernel<<<amoebaGpu->nonbondBlocks, amoebaGpu->nonbondThreadsPerBlock, sizeof(FixedFieldParticle)*amoebaGpu->nonbondThreadsPerBlock>>>(
                                                                           amoebaGpu->psWorkUnit->_pDevData,
                                                                           amoebaGpu->psWorkArray_3_1->_pDevData,
#ifdef AMOEBA_DEBUG
                                                                           amoebaGpu->psWorkArray_3_2->_pDevData,
                                                                           debugArray->_pDevData, targetAtom );
#else
                                                                           amoebaGpu->psWorkArray_3_2->_pDevData );
#endif
    } else {

#ifdef AMOEBA_DEBUG
        (void) fprintf( amoebaGpu->log, "N2 no warp\n" );
        (void) fprintf( amoebaGpu->log, "%s numBlocks=%u numThreads=%u bufferPerWarp=%u atm=%u shrd=%u ixnCt=%u workUnits=%u\n", methodName,
                        amoebaGpu->nonbondBlocks, amoebaGpu->nonbondThreadsPerBlock, amoebaGpu->bOutputBufferPerWarp,
                        sizeof(FixedFieldParticle), sizeof(FixedFieldParticle)*amoebaGpu->nonbondThreadsPerBlock, (*gpu->psInteractionCount)[0], gpu->sim.workUnits );
        (void) fflush( amoebaGpu->log );
#endif

        kCalculateAmoebaFixedE_FieldN2Forces_kernel<<<amoebaGpu->nonbondBlocks, amoebaGpu->nonbondThreadsPerBlock, sizeof(FixedFieldParticle)*amoebaGpu->nonbondThreadsPerBlock>>>(
                                                                           amoebaGpu->psWorkUnit->_pDevData,
                                                                           amoebaGpu->psWorkArray_3_1->_pDevData,
#ifdef AMOEBA_DEBUG
                                                                           amoebaGpu->psWorkArray_3_2->_pDevData,
                                                                           debugArray->_pDevData, targetAtom );
#else
                                                                           amoebaGpu->psWorkArray_3_2->_pDevData );
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

    kReduceE_Fields_kernel( amoebaGpu );

#ifdef AMOEBA_DEBUG
    if( amoebaGpu->log ){
        gpu->psInteractionCount->Download();
        (void) fprintf( amoebaGpu->log, "AmoebaN2Forces_kernel numBlocks=%u numThreads=%u bufferPerWarp=%u atm=%u shrd=%u ixnCt=%u workUnits=%u\n",
                        amoebaGpu->nonbondBlocks, amoebaGpu->nonbondThreadsPerBlock, amoebaGpu->bOutputBufferPerWarp,
                        sizeof(FixedFieldParticle), sizeof(FixedFieldParticle)*amoebaGpu->nonbondThreadsPerBlock, (*gpu->psInteractionCount)[0], gpu->sim.workUnits );
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
               float sum[2][3] = { { 0.0f, 0.0f, 0.0f },  { 0.0f, 0.0f, 0.0f } };

               (void) fprintf( amoebaGpu->log,"\n" );
               for( int jj = 0; jj < 1248; jj++ ){
                   int debugIndex = jj;
                   if( jj == ii )continue;
                   (void) fprintf( amoebaGpu->log,"\n\n%4d %4d rrs\n[%16.9e %16.9e %16.9e %16.9e]\n",
                                   ii, jj,
                                   debugArray->_pSysData[debugIndex].x, debugArray->_pSysData[debugIndex].y,
                                   debugArray->_pSysData[debugIndex].z, debugArray->_pSysData[debugIndex].w );
    
                   debugIndex += amoebaGpu->paddedNumberOfAtoms;
                   (void) fprintf( amoebaGpu->log,"[%16.9e %16.9e %16.9e]\n",
                                   debugArray->_pSysData[debugIndex].x, debugArray->_pSysData[debugIndex].y,
                                   debugArray->_pSysData[debugIndex].z );
    
                   debugIndex += amoebaGpu->paddedNumberOfAtoms;
                   (void) fprintf( amoebaGpu->log,"[%16.9e %16.9e %16.9e]\n",
                                   debugArray->_pSysData[debugIndex].x, debugArray->_pSysData[debugIndex].y,
                                   debugArray->_pSysData[debugIndex].z );
    
    
                   debugIndex += amoebaGpu->paddedNumberOfAtoms;
                   (void) fprintf( amoebaGpu->log,"[%16.9e %16.9e %16.9e]\n",
                                   debugArray->_pSysData[debugIndex].x, debugArray->_pSysData[debugIndex].y,
                                   debugArray->_pSysData[debugIndex].z );
    
    
                   debugIndex += amoebaGpu->paddedNumberOfAtoms;
                   (void) fprintf( amoebaGpu->log,"[%16.9e %16.9e %16.9e]\n",
                                   debugArray->_pSysData[debugIndex].x, debugArray->_pSysData[debugIndex].y,
                                   debugArray->_pSysData[debugIndex].z );

                   debugIndex += amoebaGpu->paddedNumberOfAtoms;
                   (void) fprintf( amoebaGpu->log,"Y1 %5d %16.9e %16.9e %16.9e\n", jj,
                                   debugArray->_pSysData[debugIndex].x, debugArray->_pSysData[debugIndex].y,
                                   debugArray->_pSysData[debugIndex].z );

                   sum[0][0] += debugArray->_pSysData[debugIndex].x;
                   sum[0][1] += debugArray->_pSysData[debugIndex].y;
                   sum[0][2] += debugArray->_pSysData[debugIndex].z;
    
                   debugIndex += amoebaGpu->paddedNumberOfAtoms;
                   (void) fprintf( amoebaGpu->log,"Y2 %5d %16.9e %16.9e %16.9e\n", jj,
                                   debugArray->_pSysData[debugIndex].x, debugArray->_pSysData[debugIndex].y,
                                   debugArray->_pSysData[debugIndex].z );

                   sum[1][0] += debugArray->_pSysData[debugIndex].x;
                   sum[1][1] += debugArray->_pSysData[debugIndex].y;
                   sum[1][2] += debugArray->_pSysData[debugIndex].z;

/*
                   debugIndex += amoebaGpu->paddedNumberOfAtoms;
                   (void) fprintf( amoebaGpu->log,"atmJ[%16.9e %16.9e %16.9e]\n",
                                   debugArray->_pSysData[debugIndex].x, debugArray->_pSysData[debugIndex].y,
                                   debugArray->_pSysData[debugIndex].z );

                   debugIndex += amoebaGpu->paddedNumberOfAtoms;
                   (void) fprintf( amoebaGpu->log,"atmJ[%16.9e %16.9e %16.9e]\n",
                                   debugArray->_pSysData[debugIndex].x, debugArray->_pSysData[debugIndex].y,
                                   debugArray->_pSysData[debugIndex].z );

    
                   debugIndex += gpu->natoms;
                   (void) fprintf( amoebaGpu->log,"[%16.9e %16.9e %16.9e %16.9e]\n",
                                   debugArray->_pSysData[debugIndex].x, debugArray->_pSysData[debugIndex].y,
                                   debugArray->_pSysData[debugIndex].z, debugArray->_pSysData[debugIndex].w );
 */   
               }
               (void) fprintf( amoebaGpu->log,"SumQ [%16.9e %16.9e %16.9e] [%16.9e %16.9e %16.9e]\n",
                               sum[0][0], sum[0][1], sum[0][2],
                               sum[1][0], sum[1][1], sum[1][2] );
           }
        for( unsigned int ii = 0; ii < debugArray->_stride; ii++ ){
           int print;
           if( debugArray->_pSysData[ii].x  != 0.0f || debugArray->_pSysData[ii].y != 0.0f ||
               debugArray->_pSysData[ii].y  != 0.0f || debugArray->_pSysData[ii].w != 0.0f ||
               debugArray->_pSysData[ii].x  != debugArray->_pSysData[ii].x               ||
               debugArray->_pSysData[ii].y  != debugArray->_pSysData[ii].y               ||
               debugArray->_pSysData[ii].z  != debugArray->_pSysData[ii].z               ||
               debugArray->_pSysData[ii].w  != debugArray->_pSysData[ii].w               ){
               print = 0;
           } else {
               print = 0;
           }
           if( print ){
                unsigned int atomI = ii/amoebaGpu->paddedNumberOfAtoms;
                unsigned int atomJ = ii - atomI*amoebaGpu->paddedNumberOfAtoms;
                (void) fprintf( amoebaGpu->log, "%5u [%5u %5u]  ", ii, atomI, atomJ);
                (void) fprintf( amoebaGpu->log, "%14.6e %14.6e %14.6e %14.6e\n",
                                debugArray->_pSysData[ii].x,
                                debugArray->_pSysData[ii].y,
                                debugArray->_pSysData[ii].z,
                                debugArray->_pSysData[ii].w );
            }
        }

        // write results to file

        if( 1 ){
            std::vector<int> fileId;
            //fileId.push_back( 0 );
            VectorOfDoubleVectors outputVector;
            //cudaLoadCudaFloat4Array( gpu->natoms, 3, gpu->psPosq4,              outputVector, NULL, 1.0f );
            cudaLoadCudaFloatArray( gpu->natoms,  3, amoebaGpu->psE_Field,      outputVector, NULL, 1.0f );
            cudaLoadCudaFloatArray( gpu->natoms,  3, amoebaGpu->psE_FieldPolar, outputVector, NULL, 1.0f);
            cudaWriteVectorOfDoubleVectorsToFile( "CudaEField", fileId, outputVector );

         }
         delete debugArray;
    }
#endif

}
