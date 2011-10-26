/* -------------------------------------------------------------------------- *
 *                               OpenMMFreeEnergy                             *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2009 Stanford University and the Authors.           *
 * Authors: Scott Le Grand, Peter Eastman, Mark Friedrichs                    *
 * Contributors:                                                              *
 *                                                                            *
 * Permission is hereby granted, free of charge, to any person obtaining a    *
 * copy of this software and associated documentation files (the "Software"), *
 * to deal in the Software without restriction, including without limitation  *
 * the rights to use, copy, modify, merge, publish, distribute, sublicense,   *
 * and/or sell copies of the Software, and to permit persons to whom the      *
 * Software is furnished to do so, subject to the following conditions:       *
 *                                                                            *
 * The above copyright notice and this permission notice shall be included in *
 * all copies or substantial portions of the Software.                        *
 *                                                                            *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR *
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,   *
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL    *
 * THE AUTHORS, CONTRIBUTORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,    *
 * DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR      *
 * OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE  *
 * USE OR OTHER DEALINGS IN THE SOFTWARE.                                     *
 * -------------------------------------------------------------------------- */

#ifdef WIN32
  #define _USE_MATH_DEFINES /* M_PI */
#endif

#define PARAMETER_PRINT 1
#define MAX_PARAMETER_PRINT 10

#include "openmm/OpenMMException.h"
#include "cudaKernels.h"
#include "GpuFreeEnergyCudaKernels.h"
#include "freeEnergyGpuTypes.h"

// for some reason, these are not being included w/ cudaKernels.h on Windows
//extern void OPENMMCUDA_EXPORT SetCalculateObcGbsaForces2Sim(gpuContext gpu);
extern void OPENMMCUDA_EXPORT SetForcesSim(gpuContext gpu);

#include <cmath>
#include <iostream>
#include <sstream>
#include <limits>
#include <cstring>
#include <vector>
#include <stdio.h>

#ifdef WIN32
#include <windows.h>
#else
#include <sys/time.h>
#endif

using std::vector;

extern "C"
freeEnergyGpuContext freeEnergyGpuInit( _gpuContext* gpu ){

    // allocate and zero block

    freeEnergyGpuContext freeEnergyGpu             = new _freeEnergyGpuContext;
    memset( freeEnergyGpu, 0, sizeof( struct _freeEnergyGpuContext ) );

    freeEnergyGpu->gpuContext                      = gpu; 

    return freeEnergyGpu;
}

extern "C"
void gpuPrintCudaStream( std::string name,
                         unsigned int length, unsigned int subStreams, unsigned int stride,
                         unsigned int memoryFootprint,
                         void*  pSysStream, void* pDevStream,
                         void*  pSysData,   void* pDevData, FILE* log)
{
   
    (void) fprintf( log, "     %-35s [%8u %5u %8u %8u] Stream[%p %p] Data[%16p %16p]\n",
                    name.c_str(), length, subStreams,
                    stride, memoryFootprint, pSysStream, pDevStream, pSysData, pDevData );
}

extern "C"
int gpuPrintCudaStreamFloat( CUDAStream<float>* cUDAStream, FILE* log )
{
   
    if( cUDAStream == NULL )return 0;
    gpuPrintCudaStream( cUDAStream->_name.c_str(),
                        cUDAStream->_length, cUDAStream->_subStreams, cUDAStream->_stride,
                        cUDAStream->_length*cUDAStream->_subStreams*sizeof( float ),
                        static_cast<void*>(cUDAStream->_pSysStream), static_cast<void*>(cUDAStream->_pDevStream),
                        static_cast<void*>(cUDAStream->_pSysData), static_cast<void*>(cUDAStream->_pDevData), log );
    return cUDAStream->_length*cUDAStream->_subStreams*sizeof( float );
}

extern "C"
int gpuPrintCudaStreamFloat2( CUDAStream<float2>* cUDAStream, FILE* log )
{
   
    if( cUDAStream == NULL )return 0;
    gpuPrintCudaStream( cUDAStream->_name.c_str(),
                        cUDAStream->_length, cUDAStream->_subStreams, cUDAStream->_stride,
                        cUDAStream->_length*cUDAStream->_subStreams*sizeof( float2 ),
                        static_cast<void*>(cUDAStream->_pSysStream), static_cast<void*>(cUDAStream->_pDevStream),
                        static_cast<void*>(cUDAStream->_pSysData), static_cast<void*>(cUDAStream->_pDevData), log );
    return cUDAStream->_length*cUDAStream->_subStreams*2*sizeof( float );
}

extern "C"
int gpuPrintCudaStreamFloat4( CUDAStream<float4>* cUDAStream, FILE* log )
{
   
    if( cUDAStream == NULL )return 0;
    gpuPrintCudaStream( cUDAStream->_name.c_str(),
                        cUDAStream->_length, cUDAStream->_subStreams, cUDAStream->_stride,
                        cUDAStream->_length*cUDAStream->_subStreams*sizeof( float4 ),
                        static_cast<void*>(cUDAStream->_pSysStream), static_cast<void*>(cUDAStream->_pDevStream),
                        static_cast<void*>(cUDAStream->_pSysData), static_cast<void*>(cUDAStream->_pDevData), log );
    return cUDAStream->_length*cUDAStream->_subStreams*4*sizeof( float );
}

extern "C"
int gpuPrintCudaStreamUnsignedInt( CUDAStream<unsigned int>* cUDAStream, FILE* log )
{
   
    if( cUDAStream == NULL )return 0;
    gpuPrintCudaStream( cUDAStream->_name.c_str(),
                        cUDAStream->_length, cUDAStream->_subStreams, cUDAStream->_stride,
                        cUDAStream->_length*cUDAStream->_subStreams*sizeof( unsigned int ),
                        static_cast<void*>(cUDAStream->_pSysStream), static_cast<void*>(cUDAStream->_pDevStream),
                        static_cast<void*>(cUDAStream->_pSysData), static_cast<void*>(cUDAStream->_pDevData), log );
    return cUDAStream->_length*cUDAStream->_subStreams*sizeof( unsigned int );
}

extern "C"
int gpuPrintCudaStreamInt( CUDAStream<int>* cUDAStream, FILE* log )
{
   
    if( cUDAStream == NULL )return 0;
    gpuPrintCudaStream( cUDAStream->_name.c_str(),
                        cUDAStream->_length, cUDAStream->_subStreams, cUDAStream->_stride,
                        cUDAStream->_length*cUDAStream->_subStreams*sizeof( int ),
                        static_cast<void*>(cUDAStream->_pSysStream), static_cast<void*>(cUDAStream->_pDevStream),
                        static_cast<void*>(cUDAStream->_pSysData), static_cast<void*>(cUDAStream->_pDevData), log );
    return cUDAStream->_length*cUDAStream->_subStreams*sizeof( int );
}

extern "C"
int gpuPrintCudaStreamInt2( CUDAStream<int2>* cUDAStream, FILE* log )
{
   
    if( cUDAStream == NULL )return 0;
    gpuPrintCudaStream( cUDAStream->_name.c_str(),
                        cUDAStream->_length, cUDAStream->_subStreams, cUDAStream->_stride,
                        cUDAStream->_length*cUDAStream->_subStreams*sizeof( int2 ),
                        static_cast<void*>(cUDAStream->_pSysStream), static_cast<void*>(cUDAStream->_pDevStream),
                        static_cast<void*>(cUDAStream->_pSysData), static_cast<void*>(cUDAStream->_pDevData), log );
    return cUDAStream->_length*cUDAStream->_subStreams*2*sizeof( int );
}

extern "C"
int gpuPrintCudaStreamInt4( CUDAStream<int4>* cUDAStream, FILE* log )
{
   
    if( cUDAStream == NULL )return 0;
    gpuPrintCudaStream( cUDAStream->_name.c_str(),
                        cUDAStream->_length, cUDAStream->_subStreams, cUDAStream->_stride,
                        cUDAStream->_length*cUDAStream->_subStreams*sizeof( int4 ),
                        static_cast<void*>(cUDAStream->_pSysStream), static_cast<void*>(cUDAStream->_pDevStream),
                        static_cast<void*>(cUDAStream->_pSysData), static_cast<void*>(cUDAStream->_pDevData), log );
    return cUDAStream->_length*cUDAStream->_subStreams*4*sizeof( int );
}

extern "C"
void gpuPrintCudaFreeEnergyGmxSimulation(freeEnergyGpuContext freeEnergyGpu, FILE* log )
{

    if( log == NULL )return;

    _gpuContext* gpu                            = freeEnergyGpu->gpuContext;
    int totalMemory                             = 0;

    (void) fprintf( log, "cudaFreeEnergyGmxSimulation:\n\n" );

    (void) fprintf( log, "\n" );
    (void) fprintf( log, "     numberOfAtoms                      %u\n",      gpu->natoms );
    (void) fprintf( log, "     paddedNumberOfAtoms                %u\n",      gpu->sim.paddedNumberOfAtoms );


    (void) fprintf( log, "\n\n" );
    (void) fprintf( log, "     gpuContext                         %p\n",      freeEnergyGpu->gpuContext );
    (void) fprintf( log, "     log                                %p %s\n",   freeEnergyGpu->log, freeEnergyGpu->log == stderr ? "is stderr" : "is not stderr");
    (void) fprintf( log, "     sm_version                         %u\n",      gpu->sm_version );
    (void) fprintf( log, "     device                             %u\n",      gpu->device );
    (void) fprintf( log, "     sharedMemoryPerBlock               %u\n",      gpu->sharedMemoryPerBlock );
    (void) fprintf( log, "     bOutputBufferPerWarp               %d\n",      gpu->bOutputBufferPerWarp );
    (void) fprintf( log, "     blocks                             %u\n",      gpu->sim.blocks );
    (void) fprintf( log, "     threads_per_block                  %u\n",      gpu->sim.threads_per_block);
    (void) fprintf( log, "     update_threads_per_block           %u\n",      gpu->sim.update_threads_per_block);
    (void) fprintf( log, "     nonbondBlocks                      %u\n",      gpu->sim.nonbond_blocks );
    (void) fprintf( log, "     nonbondThreadsPerBlock             %u\n",      gpu->sim.nonbond_threads_per_block);
    (void) fprintf( log, "     bsf_reduce_threads_per_block       %u\n",      gpu->sim.bsf_reduce_threads_per_block);
    (void) fprintf( log, "     nonbondOutputBuffers               %u\n",      gpu->sim.nonbondOutputBuffers );
    (void) fprintf( log, "     outputBuffers                      %u\n",      gpu->sim.outputBuffers );

    totalMemory += gpuPrintCudaStreamFloat(  freeEnergyGpu->gpuContext->psEnergy,    log );
    totalMemory += gpuPrintCudaStreamFloat4( freeEnergyGpu->gpuContext->psForce4,    log );

    (void) fflush( log );
}

extern "C"
void freeEnergyGpuShutDown( freeEnergyGpuContext freeEnergyGpu ){

    if( freeEnergyGpu->log ){
        (void) fprintf( freeEnergyGpu->log, "freeEnergyGpuShutDown called.\n" );
        (void) fflush( freeEnergyGpu->log );
    }

    // free free energy Cuda arrays
   
    delete freeEnergyGpu->psLJ14ID;
    delete freeEnergyGpu->psLJ14Parameter;

    delete freeEnergyGpu->psSigEps4;
    delete freeEnergyGpu->psSwitchDerivative;
    delete freeEnergyGpu->psNonPolarScalingFactors;

    delete freeEnergyGpu;

    return;
}

extern "C"
void freeEnergyGpuSetConstants( freeEnergyGpuContext freeEnergyGpu ){

    if( freeEnergyGpu->log ){
        (void) fprintf( freeEnergyGpu->log, "FreeEnergyGpuSetConstants called\n" );
        (void) fflush( freeEnergyGpu->log );
    }

    SetCalculateLocalSoftcoreGpuSim( freeEnergyGpu );
    SetCalculateCDLJSoftcoreGpuSim( freeEnergyGpu );
    SetCalculateGBVISoftcoreBornSumGpuSim( freeEnergyGpu );
    SetCalculateCDLJObcGbsaSoftcoreGpu1Sim( freeEnergyGpu );
    SetCalculateGBVISoftcoreForces2Sim( freeEnergyGpu );
    SetCalculateObcGbsaSoftcoreBornSumSim( freeEnergyGpu );
    SetCalculateObcGbsaSoftcoreForces2Sim( freeEnergyGpu );

}

extern "C"
void freeEnergyGpuSetPeriodicBoxSize( freeEnergyGpuContext freeEnergyGpu, float xsize, float ysize, float zsize)
{
    freeEnergyGpu->freeEnergySim.periodicBoxSizeX    = xsize;
    freeEnergyGpu->freeEnergySim.periodicBoxSizeY    = ysize;
    freeEnergyGpu->freeEnergySim.periodicBoxSizeZ    = zsize;

    freeEnergyGpu->freeEnergySim.invPeriodicBoxSizeX = 1.0f/xsize;
    freeEnergyGpu->freeEnergySim.invPeriodicBoxSizeY = 1.0f/ysize;
    freeEnergyGpu->freeEnergySim.invPeriodicBoxSizeZ = 1.0f/zsize;

    freeEnergyGpu->freeEnergySim.recipBoxSizeX       = 2.0f*PI/freeEnergyGpu->freeEnergySim.periodicBoxSizeX;
    freeEnergyGpu->freeEnergySim.recipBoxSizeY       = 2.0f*PI/freeEnergyGpu->freeEnergySim.periodicBoxSizeY;
    freeEnergyGpu->freeEnergySim.recipBoxSizeZ       = 2.0f*PI/freeEnergyGpu->freeEnergySim.periodicBoxSizeZ;

    freeEnergyGpu->freeEnergySim.cellVolume          = freeEnergyGpu->freeEnergySim.periodicBoxSizeX*freeEnergyGpu->freeEnergySim.periodicBoxSizeY*freeEnergyGpu->freeEnergySim.periodicBoxSizeZ;

    gpuSetPeriodicBoxSize( freeEnergyGpu->gpuContext, xsize, ysize, zsize );
}

extern "C"
void gpuSetNonbondedSoftcoreParameters( freeEnergyGpuContext freeEnergyGpu, float epsfac, const std::vector<int>& atom, const std::vector<float>& c6,
                                        const std::vector<float>& c12, const std::vector<float>& q,
                                        const std::vector<float>& softcoreLJLambdaArray, const std::vector<char>& symbol,
                                        const std::vector<std::vector<int> >& exclusions, CudaFreeEnergyNonbondedMethod method,
                                        float cutoffDistance, float solventDielectric ){

    unsigned int numberOfParticles                         = c6.size();
    gpuContext gpu                                         = freeEnergyGpu->gpuContext;
    int paddedNumberOfAtoms                                = gpu->sim.paddedNumberOfAtoms;

    // sanity checks

    if( paddedNumberOfAtoms < 1 ){
        std::stringstream msg;
        msg << "gpuSetNonbondedSoftcoreParameters: number of padded atoms=" <<  gpu->sim.paddedNumberOfAtoms << " is less than 1.";
        throw OpenMM::OpenMMException( msg.str() );
    }

    if( freeEnergyGpu->gpuContext->sim.atoms != numberOfParticles  ){
        std::stringstream msg;
        msg << "gpuSetNonbondedSoftcoreParameters: number of atoms in gpuContext does not match input count: " << freeEnergyGpu->gpuContext->sim.atoms << " " << numberOfParticles << ".";
        throw OpenMM::OpenMMException( msg.str() );
    }

    freeEnergyGpu->freeEnergySim.epsfac                    = epsfac;
    freeEnergyGpu->freeEnergySim.nonbondedMethod           = method;

    freeEnergyGpu->freeEnergySim.nonbondedCutoff           = cutoffDistance;
    freeEnergyGpu->freeEnergySim.nonbondedCutoffSqr        = cutoffDistance*cutoffDistance;

    gpu->sim.nonbondedCutoff                               = cutoffDistance;
    gpu->sim.nonbondedCutoffSqr                            = cutoffDistance*cutoffDistance;

    if( cutoffDistance > 0.0f ){
        freeEnergyGpu->freeEnergySim.reactionFieldK        = pow(cutoffDistance, -3.0f)*(solventDielectric-1.0f)/(2.0f*solventDielectric+1.0f);
        freeEnergyGpu->freeEnergySim.reactionFieldC        = (1.0f / cutoffDistance)*(3.0f*solventDielectric)/(2.0f*solventDielectric+1.0f);
        gpu->sim.reactionFieldK                            = freeEnergyGpu->freeEnergySim.reactionFieldK;
        gpu->sim.reactionFieldC                            = freeEnergyGpu->freeEnergySim.reactionFieldC;
    } else {
        freeEnergyGpu->freeEnergySim.reactionFieldK        = 0.0f;
        freeEnergyGpu->freeEnergySim.reactionFieldC        = 0.0f;
    }

    setExclusions( gpu, exclusions );
 
    // parameters
 
    freeEnergyGpu->psSigEps4                               = new CUDAStream<float4>( paddedNumberOfAtoms, 1, "freeEnergyGpuSigEps4");
    freeEnergyGpu->freeEnergySim.pSigEps4                  = freeEnergyGpu->psSigEps4->_pDevData;

    for( unsigned int ii = 0; ii < numberOfParticles; ii++ ){

        float p1 = 0.5f;
        float p2 = 0.0f;               

        if( (c6[ii] > 0.0f) && (c12[ii] > 0.0f) ){
            p1 = 0.5f * powf(c12[ii] / c6[ii], 1.0f / 6.0f);
            p2 = c6[ii] * sqrtf(1.0f / c12[ii]);
        }
/*
            if (symbol.size() > 0)
                freeEnergyGpu->pAtomSymbol[ii] = symbol[ii];
*/

        (*freeEnergyGpu->psSigEps4)[ii].x        = p1;
        (*freeEnergyGpu->psSigEps4)[ii].y        = p2;
        (*freeEnergyGpu->psSigEps4)[ii].z        = softcoreLJLambdaArray[ii];
        (*freeEnergyGpu->psSigEps4)[ii].w        = q[ii];
    }

    // Dummy out extra atom data

    for( unsigned int ii = numberOfParticles; ii < paddedNumberOfAtoms; ii++ ){

        (*freeEnergyGpu->psSigEps4)[ii].x              = 1.0f;
        (*freeEnergyGpu->psSigEps4)[ii].y              = 0.0f;
        (*freeEnergyGpu->psSigEps4)[ii].z              = 0.0f;
        (*freeEnergyGpu->psSigEps4)[ii].w              = 0.0f;

        (*gpu->psPosq4)[ii].x                          = 100000.0f + ii * 10.0f;
        (*gpu->psPosq4)[ii].y                          = 100000.0f + ii * 10.0f;
        (*gpu->psPosq4)[ii].z                          = 100000.0f + ii * 10.0f;
        (*gpu->psPosq4)[ii].w                          = 0.0f;

    }

    if( freeEnergyGpu->log ){
        (void) fprintf( freeEnergyGpu->log,"freeEnergyGpuSetNonbondedSoftcoreParameters: %5u padded=%u epsfac=%14.7e method=%d cutoffDistance=%9.2f solventDielectric=%9.2f\n",
                        numberOfParticles, freeEnergyGpu->gpuContext->sim.paddedNumberOfAtoms, epsfac, method, cutoffDistance, solventDielectric );
#ifdef PARAMETER_PRINT
        int maxPrint = MAX_PARAMETER_PRINT;
        for (unsigned int ii = 0; ii < numberOfParticles; ii++){
            (void) fprintf( freeEnergyGpu->log,"%6u sig[%14.7e %14.7e] lambda=%10.3f q=%10.3f\n",
                            ii, 
                            (*freeEnergyGpu->psSigEps4)[ii].x, (*freeEnergyGpu->psSigEps4)[ii].y, (*freeEnergyGpu->psSigEps4)[ii].z, (*freeEnergyGpu->psSigEps4)[ii].w );
            if( ii == maxPrint && ii < freeEnergyGpu->gpuContext->sim.paddedNumberOfAtoms - maxPrint ){
               ii = numberOfParticles - maxPrint;
            }
        }
        unsigned int offset = paddedNumberOfAtoms - maxPrint;
        if( offset > 0 ){
            if( offset > numberOfParticles ){
                (void) fprintf( freeEnergyGpu->log,"Dummy padded entries\n" );
                for (unsigned int ii = offset; ii < paddedNumberOfAtoms; ii++){
                    (void) fprintf( freeEnergyGpu->log,"%6u sig[%14.7e %14.7e] lambda=%10.3f q=%10.3f\n",
                                    ii, 
                                    (*freeEnergyGpu->psSigEps4)[ii].x, (*freeEnergyGpu->psSigEps4)[ii].y, (*freeEnergyGpu->psSigEps4)[ii].z, (*freeEnergyGpu->psSigEps4)[ii].w );
                }
            }
        }
#endif
        (void) fflush( freeEnergyGpu->log );
    }
 
    // upload data to board

    freeEnergyGpu->psSigEps4->Upload();
    gpu->psPosq4->Upload();

    return;
}

/**---------------------------------------------------------------------------------------

   Get threads/block

   @param amoebaGpu              amoebaGpuContext
   @param sharedMemoryPerThread  shared memory/thread
   @param sharedMemoryPerBlock   shared memory/block

   @return threadsPerBlock

   --------------------------------------------------------------------------------------- */

extern "C"
unsigned int getThreadsPerBlockFEP( freeEnergyGpuContext freeEnergyGpu, unsigned int sharedMemoryPerThread, unsigned int sharedMemoryPerBlock )
{
    unsigned int grid               = freeEnergyGpu->gpuContext->grid;
    unsigned int threadsPerBlock    = (sharedMemoryPerBlock + grid -1)/(grid*sharedMemoryPerThread);
    threadsPerBlock                 = threadsPerBlock < 1 ? 1 : threadsPerBlock;
    threadsPerBlock                *= grid;

   return threadsPerBlock;
}

static int decodeCell( int cellCode, unsigned int* x, unsigned int* y, unsigned int* exclusion ){
    *x         =  cellCode >> 17; 
    *y         = (cellCode >> 2 ) & 0x7FFF;
    *exclusion = (cellCode & 1) ? 1 : 0;  
   return 0;
}

void showWorkUnitsFreeEnergy( freeEnergyGpuContext freeEnergyGpu, int interactingWorkUnit ){

    gpuContext gpu = freeEnergyGpu->gpuContext;

    gpu->psWorkUnit->Download();
    gpu->psInteractingWorkUnit->Download();
    gpu->psInteractionFlag->Download();

    unsigned int totalWarps      = (gpu->sim.nonbond_blocks*gpu->sim.nonbond_threads_per_block)/GRID;
    //unsigned int warp            = (blockIdx.x*blockDim.x+threadIdx.x)/GRID;
    //unsigned int numWorkUnits    = cSim.pInteractionCount[0];
    unsigned int numWorkUnits    = gpu->psInteractionCount->_pSysData[0];

    (void) fprintf( stderr, "Total warps=%u blocks=%u threads=%u GRID=%u wus=%u warpOutput=%d\n",
                    totalWarps, gpu->sim.nonbond_blocks, gpu->sim.nonbond_threads_per_block, GRID, 
                    numWorkUnits, gpu->bOutputBufferPerWarp );

    unsigned int maxPrint = 3;
    std::stringstream message;
    char buffer[2048];
    unsigned int targetAtom = 18;
    for( unsigned int ii = 0; ii < gpu->sim.nonbond_blocks; ii++ )
    {
        unsigned int blockId = ii;
        for( unsigned int jj = 0; jj < gpu->sim.nonbond_threads_per_block; jj++ )
        {
            unsigned int warp    = (ii*gpu->sim.nonbond_threads_per_block+jj)/GRID;
            unsigned int pos     = warp*numWorkUnits/totalWarps;
            unsigned int end     = (warp+1)*numWorkUnits/totalWarps;
            unsigned int print   = 0;
            while( pos < end ){
                unsigned int x, y, exclusion, flags;
                int flagInt;
                if( interactingWorkUnit ){
                    decodeCell( gpu->psInteractingWorkUnit->_pSysData[pos], &x, &y, &exclusion );
                    flags = gpu->psInteractionFlag->_pSysData[pos];
                    if( flags == 0xFFFFFFFF ){
                        flagInt = -2;
                    } else {
                        flagInt = flags;
                    }
                } else {
                    decodeCell( gpu->psWorkUnit->_pSysData[pos], &x, &y, &exclusion );
                    flagInt = -1;
                }

                x                   *= GRID;
                y                   *= GRID;
                unsigned int bufferX = gpu->bOutputBufferPerWarp ? warp*gpu->sim.stride : (x >> GRIDBITS) * gpu->sim.stride;
                unsigned int bufferY = gpu->bOutputBufferPerWarp ? warp*gpu->sim.stride : (y >> GRIDBITS) * gpu->sim.stride;
                bufferX             /= gpu->sim.paddedNumberOfAtoms;
                bufferY             /= gpu->sim.paddedNumberOfAtoms;
                if( jj == 1 ){
                //if( bufferX >= 62 || bufferY >= 62 ){
                     (void) sprintf( buffer, "Block %4u thread %4u warp=%4u pos[%4u %4u] bufXY[%4u %4u]", ii, jj, warp, pos, end, bufferX, bufferY );
                     message << buffer;
                     (void) sprintf( buffer, "    x[%4u %4u] y[%4u %4u] excl=%u", x, x+32, y, y+32, exclusion );
                     message << buffer;
                     if( interactingWorkUnit ){
                         (void) sprintf( buffer, " Flg=%d (-2=all) %u", flagInt, flags );
                     }
                     message << buffer;
                     message << std::endl;
      
                }
                pos++;
            }
        }
    }
    (void) fprintf( stderr, "%s\n\n", message.str().c_str() );

}
