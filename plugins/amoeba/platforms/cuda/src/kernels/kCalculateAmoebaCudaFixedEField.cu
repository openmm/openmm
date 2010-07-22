
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
                               amoebaGpu->psWorkArray_3_1->_pDevStream[0], amoebaGpu->psE_Field->_pDevStream[0] );
    LAUNCHERROR("kReduceE_Fields1");

    kReduceFields_kernel<<<amoebaGpu->nonbondBlocks, amoebaGpu->fieldReduceThreadsPerBlock>>>(
                               amoebaGpu->paddedNumberOfAtoms*3, amoebaGpu->outputBuffers,
                               amoebaGpu->psWorkArray_3_2->_pDevStream[0], amoebaGpu->psE_FieldPolar->_pDevStream[0] );
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
                            amoebaGpu->psWorkArray_3_1->_pSysStream[0][ii],
                            amoebaGpu->psWorkArray_3_1->_pSysStream[0][ii+1],
                            amoebaGpu->psWorkArray_3_1->_pSysStream[0][ii+2],
                            amoebaGpu->psWorkArray_3_2->_pSysStream[0][ii],
                            amoebaGpu->psWorkArray_3_2->_pSysStream[0][ii+1],
                            amoebaGpu->psWorkArray_3_2->_pSysStream[0][ii+2] );
    } 
}

static void printEFieldAtomBuffers( amoebaGpuContext amoebaGpu, unsigned int targetAtom )
{
    (void) fprintf( amoebaGpu->log, "EField atom %u\n", targetAtom );
    for( unsigned int ii = 0; ii < amoebaGpu->outputBuffers; ii++ ){
        unsigned int particleIndex = targetAtom + ii*3*amoebaGpu->paddedNumberOfAtoms;
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
    memset( debugArray->_pSysStream[0],      0, sizeof( float )*4*paddedNumberOfAtoms*paddedNumberOfAtoms);
    debugArray->Upload();

    (*gpu->psInteractionCount)[0]              = gpu->sim.workUnits;
    gpu->psInteractionCount->Upload();

    // print intermediate results for the targetAtom 

    unsigned int targetAtom  = 3;
#endif

    kClearFields_3( amoebaGpu, 2 );

    if (gpu->bOutputBufferPerWarp){
        (void) fprintf( amoebaGpu->log, "N2 warp\n" );
        kCalculateAmoebaFixedE_FieldN2ByWarpForces_kernel<<<amoebaGpu->nonbondBlocks, amoebaGpu->nonbondThreadsPerBlock, sizeof(FixedFieldParticle)*amoebaGpu->nonbondThreadsPerBlock>>>(
                                                                           amoebaGpu->psWorkUnit->_pDevStream[0],
                                                                           gpu->psPosq4->_pDevStream[0],
                                                                           amoebaGpu->psLabFrameDipole->_pDevStream[0],
                                                                           amoebaGpu->psLabFrameQuadrupole->_pDevStream[0],
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
                        sizeof(FixedFieldParticle), sizeof(FixedFieldParticle)*amoebaGpu->nonbondThreadsPerBlock, amoebaGpu->energyOutputBuffers, (*gpu->psInteractionCount)[0], gpu->sim.workUnits );
        (void) fflush( amoebaGpu->log );
#endif

        kCalculateAmoebaFixedE_FieldN2Forces_kernel<<<amoebaGpu->nonbondBlocks, amoebaGpu->nonbondThreadsPerBlock, sizeof(FixedFieldParticle)*amoebaGpu->nonbondThreadsPerBlock>>>(
                                                                           amoebaGpu->psWorkUnit->_pDevStream[0],
                                                                           gpu->psPosq4->_pDevStream[0],
                                                                           amoebaGpu->psLabFrameDipole->_pDevStream[0],
                                                                           amoebaGpu->psLabFrameQuadrupole->_pDevStream[0],
                                                                           amoebaGpu->psWorkArray_3_1->_pDevStream[0],
#ifdef AMOEBA_DEBUG
                                                                           amoebaGpu->psWorkArray_3_2->_pDevStream[0],
                                                                           debugArray->_pDevStream[0], targetAtom );
#else
                                                                           amoebaGpu->psWorkArray_3_2->_pDevStream[0] );
#endif
    }

    LAUNCHERROR("kCalculateAmoebaFixedE_FieldN2Forces_kernel");

#if 0
        for( unsigned int ii = 0; ii < amoebaGpu->outputBuffers; ii++ ){
            //float index = 1.0f;
            float index = (float) ii;
            for( unsigned int jj = 0; jj < 3*amoebaGpu->paddedNumberOfAtoms; jj += 3 ){
                unsigned int kk = 3*ii*amoebaGpu->paddedNumberOfAtoms + jj;
                amoebaGpu->psWorkArray_3_1->_pSysStream[0][kk]   = index;
                amoebaGpu->psWorkArray_3_1->_pSysStream[0][kk+1] = index;
                amoebaGpu->psWorkArray_3_1->_pSysStream[0][kk+2] = index;
            }
        }
        amoebaGpu->psWorkArray_3_1->Upload();
#endif

    kReduceE_Fields_kernel( amoebaGpu );

#ifdef AMOEBA_DEBUG
    if( amoebaGpu->log ){
        gpu->psInteractionCount->Download();
        (void) fprintf( amoebaGpu->log, "AmoebaN2Forces_kernel numBlocks=%u numThreads=%u bufferPerWarp=%u atm=%u shrd=%u Ebuf=%u ixnCt=%u workUnits=%u\n",
                        amoebaGpu->nonbondBlocks, amoebaGpu->nonbondThreadsPerBlock, amoebaGpu->bOutputBufferPerWarp,
                        sizeof(FixedFieldParticle), sizeof(FixedFieldParticle)*amoebaGpu->nonbondThreadsPerBlock, amoebaGpu->energyOutputBuffers, (*gpu->psInteractionCount)[0], gpu->sim.workUnits );
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
                           amoebaGpu->psE_Field->_pSysStream[0][indexOffset],
                           amoebaGpu->psE_Field->_pSysStream[0][indexOffset+1],
                           amoebaGpu->psE_Field->_pSysStream[0][indexOffset+2] );
   
           // E_Field polar

           (void) fprintf( amoebaGpu->log,"Epol[%16.9e %16.9e %16.9e] ",
                           amoebaGpu->psE_FieldPolar->_pSysStream[0][indexOffset],
                           amoebaGpu->psE_FieldPolar->_pSysStream[0][indexOffset+1],
                           amoebaGpu->psE_FieldPolar->_pSysStream[0][indexOffset+2] );

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
                                   debugArray->_pSysStream[0][debugIndex].x, debugArray->_pSysStream[0][debugIndex].y,
                                   debugArray->_pSysStream[0][debugIndex].z, debugArray->_pSysStream[0][debugIndex].w );
    
                   debugIndex += amoebaGpu->paddedNumberOfAtoms;
                   (void) fprintf( amoebaGpu->log,"[%16.9e %16.9e %16.9e]\n",
                                   debugArray->_pSysStream[0][debugIndex].x, debugArray->_pSysStream[0][debugIndex].y,
                                   debugArray->_pSysStream[0][debugIndex].z );
    
                   debugIndex += amoebaGpu->paddedNumberOfAtoms;
                   (void) fprintf( amoebaGpu->log,"[%16.9e %16.9e %16.9e]\n",
                                   debugArray->_pSysStream[0][debugIndex].x, debugArray->_pSysStream[0][debugIndex].y,
                                   debugArray->_pSysStream[0][debugIndex].z );
    
    
                   debugIndex += amoebaGpu->paddedNumberOfAtoms;
                   (void) fprintf( amoebaGpu->log,"[%16.9e %16.9e %16.9e]\n",
                                   debugArray->_pSysStream[0][debugIndex].x, debugArray->_pSysStream[0][debugIndex].y,
                                   debugArray->_pSysStream[0][debugIndex].z );
    
    
                   debugIndex += amoebaGpu->paddedNumberOfAtoms;
                   (void) fprintf( amoebaGpu->log,"[%16.9e %16.9e %16.9e]\n",
                                   debugArray->_pSysStream[0][debugIndex].x, debugArray->_pSysStream[0][debugIndex].y,
                                   debugArray->_pSysStream[0][debugIndex].z );

                   debugIndex += amoebaGpu->paddedNumberOfAtoms;
                   (void) fprintf( amoebaGpu->log,"Y1 %5d %16.9e %16.9e %16.9e\n", jj,
                                   debugArray->_pSysStream[0][debugIndex].x, debugArray->_pSysStream[0][debugIndex].y,
                                   debugArray->_pSysStream[0][debugIndex].z );

                   sum[0][0] += debugArray->_pSysStream[0][debugIndex].x;
                   sum[0][1] += debugArray->_pSysStream[0][debugIndex].y;
                   sum[0][2] += debugArray->_pSysStream[0][debugIndex].z;
    
                   debugIndex += amoebaGpu->paddedNumberOfAtoms;
                   (void) fprintf( amoebaGpu->log,"Y2 %5d %16.9e %16.9e %16.9e\n", jj,
                                   debugArray->_pSysStream[0][debugIndex].x, debugArray->_pSysStream[0][debugIndex].y,
                                   debugArray->_pSysStream[0][debugIndex].z );

                   sum[1][0] += debugArray->_pSysStream[0][debugIndex].x;
                   sum[1][1] += debugArray->_pSysStream[0][debugIndex].y;
                   sum[1][2] += debugArray->_pSysStream[0][debugIndex].z;

/*
                   debugIndex += amoebaGpu->paddedNumberOfAtoms;
                   (void) fprintf( amoebaGpu->log,"atmJ[%16.9e %16.9e %16.9e]\n",
                                   debugArray->_pSysStream[0][debugIndex].x, debugArray->_pSysStream[0][debugIndex].y,
                                   debugArray->_pSysStream[0][debugIndex].z );

                   debugIndex += amoebaGpu->paddedNumberOfAtoms;
                   (void) fprintf( amoebaGpu->log,"atmJ[%16.9e %16.9e %16.9e]\n",
                                   debugArray->_pSysStream[0][debugIndex].x, debugArray->_pSysStream[0][debugIndex].y,
                                   debugArray->_pSysStream[0][debugIndex].z );

    
                   debugIndex += gpu->natoms;
                   (void) fprintf( amoebaGpu->log,"[%16.9e %16.9e %16.9e %16.9e]\n",
                                   debugArray->_pSysStream[0][debugIndex].x, debugArray->_pSysStream[0][debugIndex].y,
                                   debugArray->_pSysStream[0][debugIndex].z, debugArray->_pSysStream[0][debugIndex].w );
 */   
               }
               (void) fprintf( amoebaGpu->log,"SumQ [%16.9e %16.9e %16.9e] [%16.9e %16.9e %16.9e]\n",
                               sum[0][0], sum[0][1], sum[0][2],
                               sum[1][0], sum[1][1], sum[1][2] );
           }
        for( unsigned int ii = 0; ii < debugArray->_stride; ii++ ){
           int print;
           if( debugArray->_pSysStream[0][ii].x  != 0.0f || debugArray->_pSysStream[0][ii].y != 0.0f ||
               debugArray->_pSysStream[0][ii].y  != 0.0f || debugArray->_pSysStream[0][ii].w != 0.0f ||
               debugArray->_pSysStream[0][ii].x  != debugArray->_pSysStream[0][ii].x               ||
               debugArray->_pSysStream[0][ii].y  != debugArray->_pSysStream[0][ii].y               ||
               debugArray->_pSysStream[0][ii].z  != debugArray->_pSysStream[0][ii].z               ||
               debugArray->_pSysStream[0][ii].w  != debugArray->_pSysStream[0][ii].w               ){
               print = 0;
           } else {
               print = 0;
           }
           if( print ){
                unsigned int atomI = ii/amoebaGpu->paddedNumberOfAtoms;
                unsigned int atomJ = ii - atomI*amoebaGpu->paddedNumberOfAtoms;
                (void) fprintf( amoebaGpu->log, "%5u [%5u %5u]  ", ii, atomI, atomJ);
                (void) fprintf( amoebaGpu->log, "%14.6e %14.6e %14.6e %14.6e\n",
                                debugArray->_pSysStream[0][ii].x,
                                debugArray->_pSysStream[0][ii].y,
                                debugArray->_pSysStream[0][ii].z,
                                debugArray->_pSysStream[0][ii].w );
            }
        }

        // write results to file

        if( 1 ){
            std::vector<int> fileId;
            //fileId.push_back( 0 );
            VectorOfDoubleVectors outputVector;
            cudaLoadCudaFloat4Array( gpu->natoms, 3, gpu->psPosq4,              outputVector );
            cudaLoadCudaFloatArray( gpu->natoms,  3, amoebaGpu->psE_Field,      outputVector );
            cudaLoadCudaFloatArray( gpu->natoms,  3, amoebaGpu->psE_FieldPolar, outputVector);
            cudaWriteVectorOfDoubleVectorsToFile( "CudaEField", fileId, outputVector );

         }
         delete debugArray;
    }
#endif

}
