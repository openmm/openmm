/* -------------------------------------------------------------------------- *
 *                               AmoebaOpenMM                                 *
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

#include "cudaKernels.h"
#include "amoebaCudaKernels.h"

// for some reason, these are not being included w/ cudaKernels.h on Windows
extern void OPENMMCUDA_EXPORT SetCalculateObcGbsaForces2Sim(gpuContext gpu);
extern void OPENMMCUDA_EXPORT SetForcesSim(gpuContext gpu);

#include <sstream>
#include <limits>
#include <cstring>

#ifdef WIN32
//  #include <windows.h>
#endif

#define DUMP_PARAMETERS 0
//#define AMOEBA_DEBUG
#undef  AMOEBA_DEBUG

extern "C"
amoebaGpuContext amoebaGpuInit( _gpuContext* gpu )
{
    amoebaGpuContext amoebaGpu                 = new _amoebaGpuContext;

    // zero block

    memset( amoebaGpu, 0, sizeof( struct _amoebaGpuContext ) );

    amoebaGpu->gpuContext                      = gpu; 
#ifdef AMOEBA_DEBUG
    amoebaGpu->log                             = stderr; 
#endif
    amoebaGpu->numberOfSorWorkVectors          = 4; 
    
    amoebaGpu->paddedNumberOfAtoms             = gpu->sim.paddedNumberOfAtoms; 
    amoebaGpu->amoebaSim.numberOfAtoms         = gpu->natoms;
    amoebaGpu->amoebaSim.paddedNumberOfAtoms   = gpu->sim.paddedNumberOfAtoms;

    return amoebaGpu;
}

extern "C"
void gpuPrintCudaStream( std::string name,
                         unsigned int length, unsigned int subStreams, unsigned int stride,
                         void*  pSysStream, void* pDevStream,
                         void*  pSysData,   void* pDevData, FILE* log)
{
   
    (void) fprintf( log, "     %-35s [%8u %5u %8u] Stream[%p %p] Data[%16p %16p]\n",
                    name.c_str(), length, subStreams,
                    stride, pSysStream, pDevStream, pSysData, pDevData );
}

extern "C"
void gpuPrintCudaStreamFloat( CUDAStream<float>* cUDAStream, FILE* log )
{
   
    if( cUDAStream == NULL )return;
    gpuPrintCudaStream( cUDAStream->_name.c_str(),
                        cUDAStream->_length, cUDAStream->_subStreams, cUDAStream->_stride,
                        static_cast<void*>(cUDAStream->_pSysStream), static_cast<void*>(cUDAStream->_pDevStream),
                        static_cast<void*>(cUDAStream->_pSysData), static_cast<void*>(cUDAStream->_pDevData), log );
}

extern "C"
void gpuPrintCudaStreamFloat2( CUDAStream<float2>* cUDAStream, FILE* log )
{
   
    if( cUDAStream == NULL )return;
    gpuPrintCudaStream( cUDAStream->_name.c_str(),
                        cUDAStream->_length, cUDAStream->_subStreams, cUDAStream->_stride,
                        static_cast<void*>(cUDAStream->_pSysStream), static_cast<void*>(cUDAStream->_pDevStream),
                        static_cast<void*>(cUDAStream->_pSysData), static_cast<void*>(cUDAStream->_pDevData), log );
}

extern "C"
void gpuPrintCudaStreamFloat4( CUDAStream<float4>* cUDAStream, FILE* log )
{
   
    if( cUDAStream == NULL )return;
    gpuPrintCudaStream( cUDAStream->_name.c_str(),
                        cUDAStream->_length, cUDAStream->_subStreams, cUDAStream->_stride,
                        static_cast<void*>(cUDAStream->_pSysStream), static_cast<void*>(cUDAStream->_pDevStream),
                        static_cast<void*>(cUDAStream->_pSysData), static_cast<void*>(cUDAStream->_pDevData), log );
}

extern "C"
void gpuPrintCudaStreamUnsignedInt( CUDAStream<unsigned int>* cUDAStream, FILE* log )
{
   
    if( cUDAStream == NULL )return;
    gpuPrintCudaStream( cUDAStream->_name.c_str(),
                        cUDAStream->_length, cUDAStream->_subStreams, cUDAStream->_stride,
                        static_cast<void*>(cUDAStream->_pSysStream), static_cast<void*>(cUDAStream->_pDevStream),
                        static_cast<void*>(cUDAStream->_pSysData), static_cast<void*>(cUDAStream->_pDevData), log );
}

extern "C"
void gpuPrintCudaStreamInt( CUDAStream<int>* cUDAStream, FILE* log )
{
   
    if( cUDAStream == NULL )return;
    gpuPrintCudaStream( cUDAStream->_name.c_str(),
                        cUDAStream->_length, cUDAStream->_subStreams, cUDAStream->_stride,
                        static_cast<void*>(cUDAStream->_pSysStream), static_cast<void*>(cUDAStream->_pDevStream),
                        static_cast<void*>(cUDAStream->_pSysData), static_cast<void*>(cUDAStream->_pDevData), log );
}

extern "C"
void gpuPrintCudaStreamInt2( CUDAStream<int2>* cUDAStream, FILE* log )
{
   
    if( cUDAStream == NULL )return;
    gpuPrintCudaStream( cUDAStream->_name.c_str(),
                        cUDAStream->_length, cUDAStream->_subStreams, cUDAStream->_stride,
                        static_cast<void*>(cUDAStream->_pSysStream), static_cast<void*>(cUDAStream->_pDevStream),
                        static_cast<void*>(cUDAStream->_pSysData), static_cast<void*>(cUDAStream->_pDevData), log );
}

extern "C"
void gpuPrintCudaStreamInt4( CUDAStream<int4>* cUDAStream, FILE* log )
{
   
    if( cUDAStream == NULL )return;
    gpuPrintCudaStream( cUDAStream->_name.c_str(),
                        cUDAStream->_length, cUDAStream->_subStreams, cUDAStream->_stride,
                        static_cast<void*>(cUDAStream->_pSysStream), static_cast<void*>(cUDAStream->_pDevStream),
                        static_cast<void*>(cUDAStream->_pSysData), static_cast<void*>(cUDAStream->_pDevData), log );
}

extern "C"
void gpuPrintCudaAmoebaGmxSimulation(amoebaGpuContext amoebaGpu, FILE* log )
{
    if( log == NULL )return;

    _gpuContext* gpu                            = amoebaGpu->gpuContext;
    (void) fprintf( log, "cudaAmoebaGmxSimulation:\n\n" );

    (void) fprintf( log, "\n" );
    (void) fprintf( log, "     numberOfAtoms                      %u\n",      amoebaGpu->amoebaSim.numberOfAtoms );
    (void) fprintf( log, "     paddedNumberOfAtoms                %u\n",      amoebaGpu->amoebaSim.paddedNumberOfAtoms );


    (void) fprintf( log, "\n\n" );
    (void) fprintf( log, "     gpuContext                         %p\n",      amoebaGpu->gpuContext );
    (void) fprintf( log, "     log                                %p\n",      amoebaGpu->log );
    (void) fprintf( log, "     sm_version                         %u\n",      gpu->sm_version );
    (void) fprintf( log, "     device                             %u\n",      gpu->device );
    (void) fprintf( log, "     sharedMemoryPerBlock               %u\n",      gpu->sharedMemoryPerBlock );
    (void) fprintf( log, "     pMapArray                          %p\n",      amoebaGpu->pMapArray );
    (void) fprintf( log, "     dMapArray                          %p\n",      amoebaGpu->dMapArray );
    (void) fprintf( log, "     bOutputBufferPerWarp               %d\n",      amoebaGpu->bOutputBufferPerWarp );
    (void) fprintf( log, "     paddedNumberOfAtoms                %u\n",      amoebaGpu->paddedNumberOfAtoms );
    (void) fprintf( log, "     nonbondBlocks                      %u\n",      amoebaGpu->nonbondBlocks );
    (void) fprintf( log, "     nonbondThreadsPerBlock             %u\n",      amoebaGpu->nonbondThreadsPerBlock );
    (void) fprintf( log, "     nonbondElectrostaticThreadsPerBlock%u\n",      amoebaGpu->nonbondElectrostaticThreadsPerBlock );
    (void) fprintf( log, "     nonbondOutputBuffers               %u\n",      amoebaGpu->nonbondOutputBuffers );
    (void) fprintf( log, "     energyOutputBuffers                %u\n",      amoebaGpu->energyOutputBuffers );
    (void) fprintf( log, "     threadsPerBlock                    %u\n",      amoebaGpu->threadsPerBlock );
    (void) fprintf( log, "     fieldReduceThreadsPerBlock         %u\n",      amoebaGpu->fieldReduceThreadsPerBlock );
    (void) fprintf( log, "     outputBuffers                      %u\n",      amoebaGpu->outputBuffers );
    (void) fprintf( log, "     workUnits                          %u\n",      amoebaGpu->workUnits );

    gpuPrintCudaStreamFloat(  amoebaGpu->psWorkArray_3_1, log );
    gpuPrintCudaStreamFloat(  amoebaGpu->psWorkArray_3_2, log );
    gpuPrintCudaStreamFloat(  amoebaGpu->psWorkArray_3_3, log );
    gpuPrintCudaStreamFloat(  amoebaGpu->psWorkArray_3_4, log );
    gpuPrintCudaStreamFloat(  amoebaGpu->psWorkArray_3_5, log );
    gpuPrintCudaStreamFloat(  amoebaGpu->psWorkArray_3_6, log );

    gpuPrintCudaStreamFloat(  amoebaGpu->psWorkArray_1_1, log );
    gpuPrintCudaStreamFloat(  amoebaGpu->psWorkArray_1_2, log );
    
    (void) fprintf( log, "\n\n" );

    gpuPrintCudaStreamUnsignedInt( amoebaGpu->psWorkUnit, log );
    gpuPrintCudaStreamInt(  amoebaGpu->psScalingIndicesIndex, log );
    gpuPrintCudaStreamInt(  amoebaGpu->ps_D_ScaleIndices, log );
    gpuPrintCudaStreamInt2( amoebaGpu->ps_P_ScaleIndices, log );
    gpuPrintCudaStreamInt2( amoebaGpu->ps_M_ScaleIndices, log );

    if( amoebaGpu->psAmoebaBondParameter)(void) fprintf( log, "\n" );
    gpuPrintCudaStreamInt4( amoebaGpu->psAmoebaBondID, log );
    gpuPrintCudaStreamFloat2( amoebaGpu->psAmoebaBondParameter, log );
    (void) fprintf( log, "     amoebaBonds                       %u\n",      amoebaGpu->amoebaSim.amoebaBonds );
    (void) fprintf( log, "     amoebaBond_offset                 %u\n",      amoebaGpu->amoebaSim.amoebaBond_offset );
    (void) fprintf( log, "     cubic                             %15.7e\n",  amoebaGpu->amoebaSim.amoebaBondCubicParameter);
    (void) fprintf( log, "     quartic                           %15.7e\n",  amoebaGpu->amoebaSim.amoebaBondQuarticicParameter);
    (void) fprintf( log, "     pAmoebaBondID                     %p\n",      amoebaGpu->amoebaSim.pAmoebaBondID );
    (void) fprintf( log, "     pAmoebaBondParameter              %p\n",      amoebaGpu->amoebaSim.pAmoebaBondParameter );
    
    gpuPrintCudaStreamInt4( amoebaGpu->psAmoebaAngleID1, log );
    gpuPrintCudaStreamInt2( amoebaGpu->psAmoebaAngleID2, log );
    gpuPrintCudaStreamFloat2( amoebaGpu->psAmoebaAngleParameter, log );
    (void) fprintf( log, "\n" );
    (void) fprintf( log, "     amoebaAngles                      %u\n",      amoebaGpu->amoebaSim.amoebaAngles );
    (void) fprintf( log, "     amoebaAngle_offset                %u\n",      amoebaGpu->amoebaSim.amoebaAngle_offset );
    (void) fprintf( log, "     amoebaAngleCubicK                 %15.7e\n",  amoebaGpu->amoebaSim.amoebaAngleCubicK );
    (void) fprintf( log, "     amoebaAngleQuarticK               %15.7e\n",  amoebaGpu->amoebaSim.amoebaAngleQuarticK );
    (void) fprintf( log, "     amoebaAnglePenticK                %15.7e\n",  amoebaGpu->amoebaSim.amoebaAnglePenticK );
    (void) fprintf( log, "     amoebaAngleSexticK                %15.7e\n",  amoebaGpu->amoebaSim.amoebaAngleSexticK );
    (void) fprintf( log, "     pAmoebaAngleID1                   %p\n",      amoebaGpu->amoebaSim.pAmoebaAngleID1 );
    (void) fprintf( log, "     pAmoebaAngleID2                   %p\n",      amoebaGpu->amoebaSim.pAmoebaAngleID2 );
    (void) fprintf( log, "     pAmoebaAngleParameter             %p\n",      amoebaGpu->amoebaSim.pAmoebaAngleParameter );
    
    if( amoebaGpu->psAmoebaInPlaneAngleID1 )(void) fprintf( log, "\n" );
    gpuPrintCudaStreamInt4( amoebaGpu->psAmoebaInPlaneAngleID1, log );
    gpuPrintCudaStreamInt4( amoebaGpu->psAmoebaInPlaneAngleID2, log );
    gpuPrintCudaStreamFloat2( amoebaGpu->psAmoebaInPlaneAngleParameter, log );
    (void) fprintf( log, "\n" );
    (void) fprintf( log, "     amoebaInPlaneAngles               %u\n",      amoebaGpu->amoebaSim.amoebaInPlaneAngles );
    (void) fprintf( log, "     amoebaInPlaneAngle_offset         %u\n",      amoebaGpu->amoebaSim.amoebaInPlaneAngle_offset );
    (void) fprintf( log, "     amoebaInPlaneAngleCubicK          %15.7e\n",  amoebaGpu->amoebaSim.amoebaInPlaneAngleCubicK );
    (void) fprintf( log, "     amoebaInPlaneAngleQuarticK        %15.7e\n",  amoebaGpu->amoebaSim.amoebaInPlaneAngleQuarticK );
    (void) fprintf( log, "     amoebaInPlaneAnglePenticK         %15.7e\n",  amoebaGpu->amoebaSim.amoebaInPlaneAnglePenticK );
    (void) fprintf( log, "     amoebaInPlaneAngleSexticK         %15.7e\n",  amoebaGpu->amoebaSim.amoebaInPlaneAngleSexticK );
    (void) fprintf( log, "     pAmoebaInPlaneAngleID1            %p\n",      amoebaGpu->amoebaSim.pAmoebaInPlaneAngleID1 );
    (void) fprintf( log, "     pAmoebaInPlaneAngleID2            %p\n",      amoebaGpu->amoebaSim.pAmoebaInPlaneAngleID2 );
    (void) fprintf( log, "     pAmoebaInPlaneAngleParameter      %p\n",      amoebaGpu->amoebaSim.pAmoebaInPlaneAngleParameter );

    
    if( amoebaGpu->psAmoebaTorsionID1)(void) fprintf( log, "\n" );
    gpuPrintCudaStreamInt4( amoebaGpu->psAmoebaTorsionID1, log );
    gpuPrintCudaStreamInt4( amoebaGpu->psAmoebaTorsionID2, log );
    gpuPrintCudaStreamFloat4( amoebaGpu->psAmoebaTorsionParameter1, log );
    gpuPrintCudaStreamFloat2( amoebaGpu->psAmoebaTorsionParameter2, log );
    (void) fprintf( log, "     amoebaTorsions                    %u\n",      amoebaGpu->amoebaSim.amoebaTorsions );
    (void) fprintf( log, "     amoebaTorsion_offset              %u\n",      amoebaGpu->amoebaSim.amoebaTorsion_offset );
    (void) fprintf( log, "     pAmoebaTorsionID1                 %p\n",      amoebaGpu->amoebaSim.pAmoebaTorsionID1 );
    (void) fprintf( log, "     pAmoebaTorsionID2                 %p\n",      amoebaGpu->amoebaSim.pAmoebaTorsionID2 );
    (void) fprintf( log, "     pAmoebaTorsionParameter1          %p\n",      amoebaGpu->amoebaSim.pAmoebaTorsionParameter1 );
    (void) fprintf( log, "     pAmoebaTorsionParameter2          %p\n",      amoebaGpu->amoebaSim.pAmoebaTorsionParameter2 );
    
    if( amoebaGpu->psAmoebaPiTorsionID1)(void) fprintf( log, "\n" );
    gpuPrintCudaStreamInt4( amoebaGpu->psAmoebaPiTorsionID1, log );
    gpuPrintCudaStreamInt4( amoebaGpu->psAmoebaPiTorsionID2, log );
    gpuPrintCudaStreamInt4( amoebaGpu->psAmoebaPiTorsionID3, log );
    gpuPrintCudaStreamFloat( amoebaGpu->psAmoebaPiTorsionParameter, log );

    (void) fprintf( log, "     amoebaPiTorsions                  %u\n",      amoebaGpu->amoebaSim.amoebaPiTorsions );
    (void) fprintf( log, "     amoebaPiTorsion_offset            %u\n",      amoebaGpu->amoebaSim.amoebaPiTorsion_offset );
    (void) fprintf( log, "     pAmoebaPiTorsionID1               %p\n",      amoebaGpu->amoebaSim.pAmoebaPiTorsionID1 );
    (void) fprintf( log, "     pAmoebaPiTorsionID2               %p\n",      amoebaGpu->amoebaSim.pAmoebaPiTorsionID2 );
    (void) fprintf( log, "     pAmoebaPiTorsionID3               %p\n",      amoebaGpu->amoebaSim.pAmoebaPiTorsionID3 );
    (void) fprintf( log, "     pAmoebaPiTorsionParameter         %p\n",      amoebaGpu->amoebaSim.pAmoebaPiTorsionParameter );
    
    if( amoebaGpu->psAmoebaStretchBendID1)(void) fprintf( log, "\n" );
    gpuPrintCudaStreamInt4( amoebaGpu->psAmoebaStretchBendID1, log );
    gpuPrintCudaStreamInt2( amoebaGpu->psAmoebaStretchBendID2, log );
    gpuPrintCudaStreamFloat4( amoebaGpu->psAmoebaStretchBendParameter, log );
    (void) fprintf( log, "     amoebaStretchBend                  %u\n",      amoebaGpu->amoebaSim.amoebaStretchBends );
    (void) fprintf( log, "     amoebaStretchBend_offset           %u\n",      amoebaGpu->amoebaSim.amoebaStretchBend_offset );
    (void) fprintf( log, "     pAmoebaStretchBendID1              %p\n",      amoebaGpu->amoebaSim.pAmoebaStretchBendID1 );
    (void) fprintf( log, "     pAmoebaStretchBendID2              %p\n",      amoebaGpu->amoebaSim.pAmoebaStretchBendID2 );
    (void) fprintf( log, "     pAmoebaStretchBendParameter        %p\n",      amoebaGpu->amoebaSim.pAmoebaStretchBendParameter );

    if( amoebaGpu->psAmoebaOutOfPlaneBendID1)(void) fprintf( log, "\n" );
    gpuPrintCudaStreamInt4( amoebaGpu->psAmoebaOutOfPlaneBendID1, log );
    gpuPrintCudaStreamInt4( amoebaGpu->psAmoebaOutOfPlaneBendID2, log );
    gpuPrintCudaStreamFloat( amoebaGpu->psAmoebaOutOfPlaneBendParameter, log );
    (void) fprintf( log, "     amoebaOutOfPlaneBend               %u\n",      amoebaGpu->amoebaSim.amoebaOutOfPlaneBends );
    (void) fprintf( log, "     amoebaOutOfPlaneBend_offset        %u\n",      amoebaGpu->amoebaSim.amoebaOutOfPlaneBend_offset );
    (void) fprintf( log, "     amoebaOutOfPlaneBendCubicK         %15.7e\n",  amoebaGpu->amoebaSim.amoebaOutOfPlaneBendCubicK );
    (void) fprintf( log, "     amoebaOutOfPlaneBendQuarticK       %15.7e\n",  amoebaGpu->amoebaSim.amoebaOutOfPlaneBendQuarticK );
    (void) fprintf( log, "     amoebaOutOfPlaneBendPenticK        %15.7e\n",  amoebaGpu->amoebaSim.amoebaOutOfPlaneBendPenticK );
    (void) fprintf( log, "     amoebaOutOfPlaneBendSexticK        %15.7e\n",  amoebaGpu->amoebaSim.amoebaOutOfPlaneBendSexticK );
    (void) fprintf( log, "     pAmoebaOutOfPlaneBendID1           %p\n",      amoebaGpu->amoebaSim.pAmoebaOutOfPlaneBendID1 );
    (void) fprintf( log, "     pAmoebaOutOfPlaneBendID2           %p\n",      amoebaGpu->amoebaSim.pAmoebaOutOfPlaneBendID2 );
    (void) fprintf( log, "     pAmoebaOutOfPlaneBendParameter     %p\n",      amoebaGpu->amoebaSim.pAmoebaOutOfPlaneBendParameter );
    
    if( amoebaGpu->psAmoebaTorsionTorsionID1)(void) fprintf( log, "\n" );
    gpuPrintCudaStreamInt4( amoebaGpu->psAmoebaTorsionTorsionID1, log );
    gpuPrintCudaStreamInt4( amoebaGpu->psAmoebaTorsionTorsionID2, log );
    gpuPrintCudaStreamInt4( amoebaGpu->psAmoebaTorsionTorsionID3, log );
    gpuPrintCudaStreamFloat4( amoebaGpu->psAmoebaTorsionTorsionGrids, log );
    (void) fprintf( log, "\n" );
    (void) fprintf( log, "     amoebaTorsionTorsions              %u\n",      amoebaGpu->amoebaSim.amoebaTorsionTorsions );
    (void) fprintf( log, "     amoebaTorsionTorsion_offset        %u\n",      amoebaGpu->amoebaSim.amoebaTorsionTorsion_offset );
    (void) fprintf( log, "     pAmoebaTorsionTorsionID1           %p\n",      amoebaGpu->amoebaSim.pAmoebaTorsionTorsionID1 );
    (void) fprintf( log, "     pAmoebaTorsionTorsionID2           %p\n",      amoebaGpu->amoebaSim.pAmoebaTorsionTorsionID2 );
    (void) fprintf( log, "     pAmoebaTorsionTorsionID3           %p\n",      amoebaGpu->amoebaSim.pAmoebaTorsionTorsionID3 );

    (void) fprintf( log, "     pOutputBufferCounter               %p\n",      amoebaGpu->gpuContext->pOutputBufferCounter );

    if( amoebaGpu->psRotationMatrix)(void) fprintf( log, "\n" );
    gpuPrintCudaStreamFloat( amoebaGpu->psRotationMatrix, log );
    gpuPrintCudaStreamInt4( amoebaGpu->psMultipoleParticlesIdsAndAxisType, log );
    gpuPrintCudaStreamFloat( amoebaGpu->psMolecularDipole, log );
    gpuPrintCudaStreamFloat( amoebaGpu->psMolecularQuadrupole, log );
    
    gpuPrintCudaStreamFloat( amoebaGpu->psLabFrameDipole, log );
    gpuPrintCudaStreamFloat( amoebaGpu->psLabFrameQuadrupole, log );
    
    (void) fprintf( log, "     maxCovalentDegreeSz                %d\n",      amoebaGpu->maxCovalentDegreeSz );
    (void) fprintf( log, "     solventDielectric                  %10.3f\n",  amoebaGpu->solventDielectric);
    (void) fprintf( log, "     pGamma                             %10.3f\n",  amoebaGpu->pGamma );
    (void) fprintf( log, "     scalingDistanceCutoff              %10.3f\n",  amoebaGpu->scalingDistanceCutoff );
    (void) fprintf( log, "     scalingDistanceCutoff              %15.7e\n",  amoebaGpu->amoebaSim.scalingDistanceCutoff );
    (void) fprintf( log, "     pDampingFactorAndThole             %p\n",      amoebaGpu->amoebaSim.pDampingFactorAndThole );
    (void) fprintf( log, "     pScaleIndicesIndex                 %p\n",      amoebaGpu->amoebaSim.pScaleIndicesIndex );
    (void) fprintf( log, "     pD_ScaleIndices                    %p\n",      amoebaGpu->amoebaSim.pD_ScaleIndices );
    (void) fprintf( log, "     pP_ScaleIndices                    %p\n",      amoebaGpu->amoebaSim.pP_ScaleIndices );
    (void) fprintf( log, "     pM_ScaleIndices                    %p\n",      amoebaGpu->amoebaSim.pM_ScaleIndices );
    (void) fprintf( log, "     electric                           %15.7e\n",  amoebaGpu->amoebaSim.electric );
    (void) fprintf( log, "     gkc                                %15.7e\n",  amoebaGpu->amoebaSim.gkc );
    (void) fprintf( log, "     dielec                             %15.7e\n",  amoebaGpu->amoebaSim.dielec );
    (void) fprintf( log, "     dwater                             %15.7e\n",  amoebaGpu->amoebaSim.dwater );
    (void) fprintf( log, "     fc                                 %15.7e\n",  amoebaGpu->amoebaSim.fc );
    (void) fprintf( log, "     fd                                 %15.7e\n",  amoebaGpu->amoebaSim.fd );
    (void) fprintf( log, "     fq                                 %15.7e\n",  amoebaGpu->amoebaSim.fq );

    gpuPrintCudaStreamFloat2( amoebaGpu->psDampingFactorAndThole, log );

    gpuPrintCudaStreamInt( amoebaGpu->psCovalentDegree, log );
    gpuPrintCudaStreamInt( amoebaGpu->psPolarizationDegree, log );

    gpuPrintCudaStreamFloat( amoebaGpu->psE_Field, log );
    gpuPrintCudaStreamFloat( amoebaGpu->psE_FieldPolar, log );


    (void) fprintf( log, "     mutualInducedIterativeMethod       %d\n",  amoebaGpu->mutualInducedIterativeMethod);
    (void) fprintf( log, "     mutualInducedMaxIterations         %d\n",  amoebaGpu->mutualInducedMaxIterations);
    (void) fprintf( log, "     epsilonThreadsPerBlock             %d\n",  amoebaGpu->epsilonThreadsPerBlock);
    (void) fprintf( log, "     mutualInducedTargetEpsilon         %10.3e\n",  amoebaGpu->mutualInducedTargetEpsilon);
    (void) fprintf( log, "     mutualInducedCurrentEpsilon        %10.3e\n",  amoebaGpu->mutualInducedCurrentEpsilon );

    gpuPrintCudaStreamFloat( amoebaGpu->psInducedDipole, log );
    gpuPrintCudaStreamFloat( amoebaGpu->psInducedDipolePolar, log );
    gpuPrintCudaStreamFloat( amoebaGpu->psInducedDipolePolar, log );
    gpuPrintCudaStreamFloat( amoebaGpu->psCurrentEpsilon, log );

    (void) fprintf( log, "     numberOfSorWorkVectors             %u\n",  amoebaGpu->numberOfSorWorkVectors);
    gpuPrintCudaStreamFloat( amoebaGpu->psWorkVector[0], log );
    gpuPrintCudaStreamFloat( amoebaGpu->psForce, log );
    gpuPrintCudaStreamFloat( amoebaGpu->psTorque, log );
    gpuPrintCudaStreamFloat( amoebaGpu->psEnergy, log );
    gpuPrintCudaStreamFloat( amoebaGpu->torqueMapForce, log );

    (void) fprintf( log, "     maxMapTorqueDifference             %d\n",  amoebaGpu->maxMapTorqueDifference );
    (void) fprintf( log, "     maxMapTorqueDifferencePow2         %d\n",  amoebaGpu->maxMapTorqueDifferencePow2);

    gpuPrintCudaStreamFloat( amoebaGpu->psGk_Field, log );
    gpuPrintCudaStreamFloat( amoebaGpu->psInducedDipoleS, log );
    gpuPrintCudaStreamFloat( amoebaGpu->psInducedDipolePolarS, log );
    gpuPrintCudaStreamFloat( amoebaGpu->psBorn, log );
    gpuPrintCudaStreamFloat( amoebaGpu->psBornPolar, log );
    gpuPrintCudaStreamFloat( amoebaGpu->psKirkwoodForce, log );
    gpuPrintCudaStreamFloat( amoebaGpu->psKirkwoodEDiffForce, log );
    (void) fprintf( log, "     includeObcCavityTerm               %d\n",      amoebaGpu->includeObcCavityTerm );
    (void) fprintf( log, "     dielectricOffset                   %15.7e\n",  gpu->sim.dielectricOffset );
    (void) fprintf( log, "     probeRadius                        %15.7e\n",  gpu->sim.probeRadius );
    (void) fprintf( log, "     surfaceAreaFactor                  %15.7e\n",  gpu->sim.surfaceAreaFactor );
    (void) fprintf( log, "\n" );


    (void) fprintf( log, "     useVdwTable                        %u\n",  amoebaGpu->useVdwTable);
    (void) fprintf( log, "     vdwTableSize                       %u\n",  amoebaGpu->vdwTableSize);
    (void) fprintf( log, "     vdwSigmaCombiningRule              %d\n",  amoebaGpu->vdwSigmaCombiningRule);
    (void) fprintf( log, "     vdwEpsilonCombiningRule            %d\n",  amoebaGpu->vdwEpsilonCombiningRule);
    gpuPrintCudaStreamFloat2( amoebaGpu->psVdwSigmaEpsilon, log );
    gpuPrintCudaStreamFloat2( amoebaGpu->psVdwTable, log );
    gpuPrintCudaStreamInt( amoebaGpu->psAmoebaVdwNonReductionID, log );
    gpuPrintCudaStreamInt4( amoebaGpu->psAmoebaVdwReductionID, log );
    gpuPrintCudaStreamFloat( amoebaGpu->psAmoebaVdwReduction, log );
    gpuPrintCudaStreamFloat4( amoebaGpu->psAmoebaVdwCoordinates, log );
    gpuPrintCudaStreamUnsignedInt( amoebaGpu->psVdwWorkUnit, log );
    gpuPrintCudaStreamInt( amoebaGpu->psVdwExclusionIndicesIndex, log );
    gpuPrintCudaStreamInt( amoebaGpu->psVdwExclusionIndices, log );
    (void) fprintf( log, "     amoebaVdwNonReductions             %u\n",      amoebaGpu->amoebaSim.amoebaVdwNonReductions );
    (void) fprintf( log, "     pAmoebaVdwNonReductionID           %p\n",      amoebaGpu->amoebaSim.pAmoebaVdwNonReductionID );
    (void) fprintf( log, "     amoebaVdwReductions                %u\n",      amoebaGpu->amoebaSim.amoebaVdwReductions );
    (void) fprintf( log, "     pAmoebaVdwReductionID              %p\n",      amoebaGpu->amoebaSim.pAmoebaVdwReductionID );
    (void) fprintf( log, "     pAmoebaVdwReduction                %p\n",      amoebaGpu->amoebaSim.pAmoebaVdwReduction );
    (void) fprintf( log, "     pVdwExclusionIndicesIndex          %p\n",      amoebaGpu->amoebaSim.pVdwExclusionIndicesIndex);
    (void) fprintf( log, "     pVdwExclusionIndices               %p\n",      amoebaGpu->amoebaSim.pVdwExclusionIndices);

    gpuPrintCudaStreamFloat2( amoebaGpu->psWcaDispersionRadiusEpsilon, log );
    (void) fprintf( log, "\n" );
    (void) fprintf( log, "     epso                               %15.7e\n",  amoebaGpu->amoebaSim.epso );
    (void) fprintf( log, "     epsh                               %15.7e\n",  amoebaGpu->amoebaSim.epsh );
    (void) fprintf( log, "     rmino                              %15.7e\n",  amoebaGpu->amoebaSim.rmino );
    (void) fprintf( log, "     rminh                              %15.7e\n",  amoebaGpu->amoebaSim.rminh );
    (void) fprintf( log, "     awater                             %15.7e\n",  amoebaGpu->amoebaSim.awater );
    (void) fprintf( log, "     shctd                              %15.7e\n",  amoebaGpu->amoebaSim.shctd );
    (void) fprintf( log, "     dispoff                            %15.7e\n",  amoebaGpu->amoebaSim.dispoff );
    (void) fprintf( log, "     totalMaxWcaDispersionEnergy        %15.7e\n",  amoebaGpu->amoebaSim.totalMaxWcaDispersionEnergy );

    (void) fflush( log );

}

extern "C"
void gpuSetAmoebaBondParameters(amoebaGpuContext amoebaGpu, const std::vector<int>& particles1, const std::vector<int>& particles2, 
                                const std::vector<float>& length, const std::vector<float>& k, float cubic, float quartic)
{
    _gpuContext* gpu                                  = amoebaGpu->gpuContext;
    int bonds                                         = particles1.size();
    amoebaGpu->amoebaSim.amoebaBonds                  = bonds;

    CUDAStream<int4>* psBondID                        = new CUDAStream<int4>(bonds, 1, "AmoebaBondID");
    amoebaGpu->psAmoebaBondID                         = psBondID;
    amoebaGpu->amoebaSim.pAmoebaBondID                = psBondID->_pDevStream[0];

    CUDAStream<float2>* psBondParameter               = new CUDAStream<float2>(bonds, 1, "AmoebaBondParameter");
    amoebaGpu->psAmoebaBondParameter                  = psBondParameter;
    amoebaGpu->amoebaSim.pAmoebaBondParameter         = psBondParameter->_pDevStream[0];

    amoebaGpu->amoebaSim.amoebaBondCubicParameter     = cubic;
    amoebaGpu->amoebaSim.amoebaBondQuarticicParameter = quartic;
    for (int i = 0; i < bonds; i++)
    {
        (*psBondID)[i].x         = particles1[i];
        (*psBondID)[i].y         = particles2[i];
        (*psBondID)[i].z         = gpu->pOutputBufferCounter[(*psBondID)[i].x]++;
        (*psBondID)[i].w         = gpu->pOutputBufferCounter[(*psBondID)[i].y]++;
        (*psBondParameter)[i].x  = length[i];
        (*psBondParameter)[i].y  = k[i];

#undef DUMP_PARAMETERS
#define DUMP_PARAMETERS 5
#if (DUMP_PARAMETERS > 0 )                
if( amoebaGpu->log && (i < DUMP_PARAMETERS || i > bonds - (DUMP_PARAMETERS + 1)  ) )
        fprintf( amoebaGpu->log, "Bonds: %5d [%5d %5d %5d %5d] L=%15.7e k[%15.7e %15.7e %15.7e] [%5d %5d]\n",
            i, (*psBondID)[i].x, (*psBondID)[i].y, (*psBondID)[i].z, (*psBondID)[i].w, 
            (*psBondParameter)[i].x, (*psBondParameter)[i].y, cubic, quartic, 
            gpu->pOutputBufferCounter[(*psBondID)[i].x],
            gpu->pOutputBufferCounter[(*psBondID)[i].y] );
#endif
#undef DUMP_PARAMETERS

    }
    psBondID->Upload();
    psBondParameter->Upload();
}

extern "C"
void gpuSetAmoebaAngleParameters(amoebaGpuContext amoebaGpu, const std::vector<int>& particles1, const std::vector<int>& particles2, const std::vector<int>& particles3,
                                 const std::vector<float>& angle, const std::vector<float>& k,
                                 float cubicK, float quarticK, float penticK, float sexticK)
{

    _gpuContext* gpu                                  = amoebaGpu->gpuContext;
    int bond_angles                                   = particles1.size();
    amoebaGpu->amoebaSim.amoebaAngles                 = bond_angles;

    CUDAStream<int4>* psAngleID1                      = new CUDAStream<int4>(bond_angles, 1, "AmoebaAngleID1");
    amoebaGpu->psAmoebaAngleID1                       = psAngleID1;
    amoebaGpu->amoebaSim.pAmoebaAngleID1              = psAngleID1->_pDevStream[0];

    CUDAStream<int2>* psAngleID2                      = new CUDAStream<int2>(bond_angles, 1, "AmoebaAngleID2");
    amoebaGpu->psAmoebaAngleID2                       = psAngleID2;
    amoebaGpu->amoebaSim.pAmoebaAngleID2              = psAngleID2->_pDevStream[0];

    CUDAStream<float2>* psAngleParameter              = new CUDAStream<float2>(bond_angles, 1, "AmoebaAngleParameter");
    amoebaGpu->psAmoebaAngleParameter                 = psAngleParameter;
    amoebaGpu->amoebaSim.pAmoebaAngleParameter        = psAngleParameter->_pDevStream[0];        

    amoebaGpu->amoebaSim.amoebaAngleCubicK            = cubicK;
    amoebaGpu->amoebaSim.amoebaAngleQuarticK          = quarticK;
    amoebaGpu->amoebaSim.amoebaAnglePenticK           = penticK;
    amoebaGpu->amoebaSim.amoebaAngleSexticK           = sexticK;

    for (int i = 0; i < bond_angles; i++)
    {
        (*psAngleID1)[i].x         = particles1[i];
        (*psAngleID1)[i].y         = particles2[i];
        (*psAngleID1)[i].z         = particles3[i];
        (*psAngleParameter)[i].x   = angle[i];
        (*psAngleParameter)[i].y   = k[i];
        psAngleID1->_pSysData[i].w = gpu->pOutputBufferCounter[psAngleID1->_pSysData[i].x]++;
        psAngleID2->_pSysData[i].x = gpu->pOutputBufferCounter[psAngleID1->_pSysData[i].y]++;
        psAngleID2->_pSysData[i].y = gpu->pOutputBufferCounter[psAngleID1->_pSysData[i].z]++;

#undef DUMP_PARAMETERS
#define DUMP_PARAMETERS 5
#if (DUMP_PARAMETERS > 0 )
if( (i < DUMP_PARAMETERS || i > bond_angles - (DUMP_PARAMETERS + 1)) && amoebaGpu->log )
         fprintf( amoebaGpu->log, "Angles: %5d [%5d %5d %5d %5d] [%5d %5d] A=%15.7e k=%15.7e [%5d %5d %5d]\n", i, 
                  (*psAngleID1)[i].x, (*psAngleID1)[i].y, (*psAngleID1)[i].z, (*psAngleID1)[i].w,
                  (*psAngleID2)[i].x, (*psAngleID2)[i].y,
                  (*psAngleParameter)[i].x, (*psAngleParameter)[i].y,
                  gpu->pOutputBufferCounter[psAngleID1->_pSysData[i].x],
                  gpu->pOutputBufferCounter[psAngleID1->_pSysData[i].y],
                  gpu->pOutputBufferCounter[psAngleID1->_pSysData[i].z] );
#endif
#undef DUMP_PARAMETERS
    }
    psAngleID1->Upload();
    psAngleID2->Upload();
    psAngleParameter->Upload();
}

extern "C"
void gpuSetAmoebaInPlaneAngleParameters(amoebaGpuContext amoebaGpu, const std::vector<int>& particles1, const std::vector<int>& particles2,
                                                                    const std::vector<int>& particles3, const std::vector<int>& particles4,
                                        const std::vector<float>& angle, const std::vector<float>& k,
                                        float cubicK, float quarticK, float penticK, float sexticK)
{

    _gpuContext* gpu                                        = amoebaGpu->gpuContext;
    int bond_angles                                         = particles1.size();
    amoebaGpu->amoebaSim.amoebaInPlaneAngles                = bond_angles;

    CUDAStream<int4>* psAngleID1                            = new CUDAStream<int4>(bond_angles, 1, "AmoebaInPlaneAngleID1");
    amoebaGpu->psAmoebaInPlaneAngleID1                      = psAngleID1;
    amoebaGpu->amoebaSim.pAmoebaInPlaneAngleID1             = psAngleID1->_pDevStream[0];

    CUDAStream<int4>* psAngleID2                            = new CUDAStream<int4>(bond_angles, 1, "AmoebaInPlaneAngleID2");
    amoebaGpu->psAmoebaInPlaneAngleID2                      = psAngleID2;
    amoebaGpu->amoebaSim.pAmoebaInPlaneAngleID2             = psAngleID2->_pDevStream[0];

    CUDAStream<float2>* psAngleParameter                    = new CUDAStream<float2>(bond_angles, 1, "AmoebaInPlaneAngleParameter");
    amoebaGpu->psAmoebaInPlaneAngleParameter                = psAngleParameter;
    amoebaGpu->amoebaSim.pAmoebaInPlaneAngleParameter       = psAngleParameter->_pDevStream[0];        

    amoebaGpu->amoebaSim.amoebaInPlaneAngleCubicK           = cubicK;
    amoebaGpu->amoebaSim.amoebaInPlaneAngleQuarticK         = quarticK;
    amoebaGpu->amoebaSim.amoebaInPlaneAnglePenticK          = penticK;
    amoebaGpu->amoebaSim.amoebaInPlaneAngleSexticK          = sexticK;

    for (int i = 0; i < bond_angles; i++)
    {
        (*psAngleID1)[i].x         = particles1[i];
        (*psAngleID1)[i].y         = particles2[i];
        (*psAngleID1)[i].z         = particles3[i];
        (*psAngleID1)[i].w         = particles4[i];
        (*psAngleParameter)[i].x   = angle[i];
        (*psAngleParameter)[i].y   = k[i];
        psAngleID2->_pSysData[i].x = gpu->pOutputBufferCounter[psAngleID1->_pSysData[i].x]++;
        psAngleID2->_pSysData[i].y = gpu->pOutputBufferCounter[psAngleID1->_pSysData[i].y]++;
        psAngleID2->_pSysData[i].z = gpu->pOutputBufferCounter[psAngleID1->_pSysData[i].z]++;
        psAngleID2->_pSysData[i].w = gpu->pOutputBufferCounter[psAngleID1->_pSysData[i].w]++;

#undef DUMP_PARAMETERS
#define DUMP_PARAMETERS 5
#if (DUMP_PARAMETERS > 0 )
if( (i < DUMP_PARAMETERS || i > bond_angles - (DUMP_PARAMETERS + 1)) && amoebaGpu->log )
         fprintf( amoebaGpu->log, "InPlaneAngles: %5d [%5d %5d %5d %5d] [%5d %5d %5d %5d] A=%15.7e k=%15.7e [%5d %5d %5d %5d]\n", i, 
                  (*psAngleID1)[i].x, (*psAngleID1)[i].y, (*psAngleID1)[i].z, (*psAngleID1)[i].w,
                  (*psAngleID2)[i].x, (*psAngleID2)[i].y, (*psAngleID2)[i].z, (*psAngleID2)[i].w,
                  (*psAngleParameter)[i].x, (*psAngleParameter)[i].y,
                  gpu->pOutputBufferCounter[psAngleID1->_pSysData[i].x],
                  gpu->pOutputBufferCounter[psAngleID1->_pSysData[i].y],
                  gpu->pOutputBufferCounter[psAngleID1->_pSysData[i].z],
                  gpu->pOutputBufferCounter[psAngleID1->_pSysData[i].w] );
#endif
#undef DUMP_PARAMETERS
    }

    psAngleID1->Upload();
    psAngleID2->Upload();
    psAngleParameter->Upload();
}

extern "C"
void gpuSetAmoebaTorsionParameters(amoebaGpuContext amoebaGpu, const std::vector<int>& particles1, const std::vector<int>& particles2,
                                                               const std::vector<int>& particles3, const std::vector<int>& particles4,
                                                               const std::vector< std::vector<float> >& torsion1, 
                                                               const std::vector< std::vector<float> >& torsion2, 
                                                               const std::vector< std::vector<float> >& torsion3 ) 
{

    _gpuContext* gpu                                   = amoebaGpu->gpuContext;
    int torsions                                       = particles1.size();
    amoebaGpu->amoebaSim.amoebaTorsions                = torsions;

    CUDAStream<int4>* psTorsionID1                     = new CUDAStream<int4>(torsions, 1, "AmoebaTorsionID1");
    amoebaGpu->psAmoebaTorsionID1                      = psTorsionID1;
    amoebaGpu->amoebaSim.pAmoebaTorsionID1             = psTorsionID1->_pDevStream[0];

    CUDAStream<int4>* psTorsionID2                     = new CUDAStream<int4>(torsions, 1, "AmoebaTorsionID2");
    amoebaGpu->psAmoebaTorsionID2                      = psTorsionID2;
    amoebaGpu->amoebaSim.pAmoebaTorsionID2             = psTorsionID2->_pDevStream[0];

    CUDAStream<float4>* psTorsionParameter1            = new CUDAStream<float4>(torsions, 1, "AmoebaTorsionParameter1");
    amoebaGpu->psAmoebaTorsionParameter1               = psTorsionParameter1;
    amoebaGpu->amoebaSim.pAmoebaTorsionParameter1      = psTorsionParameter1->_pDevStream[0];        

    CUDAStream<float2>* psTorsionParameter2            = new CUDAStream<float2>(torsions, 1, "AmoebaTorsionParameter2");
    amoebaGpu->psAmoebaTorsionParameter2               = psTorsionParameter2;
    amoebaGpu->amoebaSim.pAmoebaTorsionParameter2      = psTorsionParameter2->_pDevStream[0];        

    for (int i = 0; i < torsions; i++)
    {
        (*psTorsionID1)[i].x         = particles1[i];
        (*psTorsionID1)[i].y         = particles2[i];
        (*psTorsionID1)[i].z         = particles3[i];
        (*psTorsionID1)[i].w         = particles4[i];

        (*psTorsionParameter1)[i].x  = torsion1[i][0];
        (*psTorsionParameter1)[i].y  = torsion1[i][1];
        (*psTorsionParameter1)[i].z  = torsion2[i][0];

        (*psTorsionParameter1)[i].w  = torsion2[i][1];
        (*psTorsionParameter2)[i].x  = torsion3[i][0];
        (*psTorsionParameter2)[i].y  = torsion3[i][1];

        psTorsionID2->_pSysData[i].x = gpu->pOutputBufferCounter[psTorsionID1->_pSysData[i].x]++;
        psTorsionID2->_pSysData[i].y = gpu->pOutputBufferCounter[psTorsionID1->_pSysData[i].y]++;
        psTorsionID2->_pSysData[i].z = gpu->pOutputBufferCounter[psTorsionID1->_pSysData[i].z]++;
        psTorsionID2->_pSysData[i].w = gpu->pOutputBufferCounter[psTorsionID1->_pSysData[i].w]++;

#undef DUMP_PARAMETERS
#define DUMP_PARAMETERS 5
#if (DUMP_PARAMETERS > 0 )
if( (i < DUMP_PARAMETERS || i > torsions - (DUMP_PARAMETERS + 1)) && amoebaGpu->log )
         fprintf( amoebaGpu->log, "Torsions: %5d [%5d %5d %5d %5d] [%5d %5d %5d %5d] 0[%15.7e %15.7e] 1[%15.7e %15.7e] 2[%15.7e %15.7e] [%5d %5d %5d %5d]\n", i, 
                  (*psTorsionID1)[i].x, (*psTorsionID1)[i].y, (*psTorsionID1)[i].z, (*psTorsionID1)[i].w,
                  (*psTorsionID2)[i].x, (*psTorsionID2)[i].y, (*psTorsionID2)[i].z, (*psTorsionID2)[i].w,
                  (*psTorsionParameter1)[i].x, (*psTorsionParameter1)[i].y, (*psTorsionParameter1)[i].z, (*psTorsionParameter1)[i].w,
                  (*psTorsionParameter2)[i].x, (*psTorsionParameter2)[i].y,
                  gpu->pOutputBufferCounter[psTorsionID1->_pSysData[i].x],
                  gpu->pOutputBufferCounter[psTorsionID1->_pSysData[i].y],
                  gpu->pOutputBufferCounter[psTorsionID1->_pSysData[i].z],
                  gpu->pOutputBufferCounter[psTorsionID1->_pSysData[i].w] );
#endif
#undef DUMP_PARAMETERS
    }

    psTorsionID1->Upload();
    psTorsionID2->Upload();
    psTorsionParameter1->Upload();
    psTorsionParameter2->Upload();
}

extern "C"
void gpuSetAmoebaPiTorsionParameters(amoebaGpuContext amoebaGpu, const std::vector<int>& particles1, const std::vector<int>& particles2,
                                                                 const std::vector<int>& particles3, const std::vector<int>& particles4,
                                                                 const std::vector<int>& particles5, const std::vector<int>& particles6,
                                                                 const std::vector<float>& torsionK ) 
{

    _gpuContext* gpu                                     = amoebaGpu->gpuContext;
    int piTorsions                                       = particles1.size();
    amoebaGpu->amoebaSim.amoebaPiTorsions                = piTorsions;

    CUDAStream<int4>* psPiTorsionID1                     = new CUDAStream<int4>(piTorsions, 1, "AmoebaPiTorsionID1");
    amoebaGpu->psAmoebaPiTorsionID1                      = psPiTorsionID1;
    amoebaGpu->amoebaSim.pAmoebaPiTorsionID1             = psPiTorsionID1->_pDevStream[0];

    CUDAStream<int4>* psPiTorsionID2                     = new CUDAStream<int4>(piTorsions, 1, "AmoebaPiTorsionID2");
    amoebaGpu->psAmoebaPiTorsionID2                      = psPiTorsionID2;
    amoebaGpu->amoebaSim.pAmoebaPiTorsionID2             = psPiTorsionID2->_pDevStream[0];

    CUDAStream<int4>* psPiTorsionID3                     = new CUDAStream<int4>(piTorsions, 1, "AmoebaPiTorsionID3");
    amoebaGpu->psAmoebaPiTorsionID3                      = psPiTorsionID3;
    amoebaGpu->amoebaSim.pAmoebaPiTorsionID3             = psPiTorsionID3->_pDevStream[0];

    CUDAStream<float>* psPiTorsionParameter              = new CUDAStream<float>(piTorsions, 1, "AmoebaPiTorsionParameter1");
    amoebaGpu->psAmoebaPiTorsionParameter                = psPiTorsionParameter;
    amoebaGpu->amoebaSim.pAmoebaPiTorsionParameter       = psPiTorsionParameter->_pDevStream[0];        

    for (int i = 0; i < piTorsions; i++)
    {
        (*psPiTorsionID1)[i].x         = particles1[i];
        (*psPiTorsionID1)[i].y         = particles2[i];
        (*psPiTorsionID1)[i].z         = particles3[i];
        (*psPiTorsionID1)[i].w         = particles4[i];
        (*psPiTorsionID2)[i].x         = particles5[i];
        (*psPiTorsionID2)[i].y         = particles6[i];

        (*psPiTorsionParameter)[i]     = torsionK[i];

        psPiTorsionID2->_pSysData[i].z = gpu->pOutputBufferCounter[psPiTorsionID1->_pSysData[i].x]++;
        psPiTorsionID2->_pSysData[i].w = gpu->pOutputBufferCounter[psPiTorsionID1->_pSysData[i].y]++;
        psPiTorsionID3->_pSysData[i].x = gpu->pOutputBufferCounter[psPiTorsionID1->_pSysData[i].z]++;
        psPiTorsionID3->_pSysData[i].y = gpu->pOutputBufferCounter[psPiTorsionID1->_pSysData[i].w]++;
        psPiTorsionID3->_pSysData[i].z = gpu->pOutputBufferCounter[psPiTorsionID2->_pSysData[i].x]++;
        psPiTorsionID3->_pSysData[i].w = gpu->pOutputBufferCounter[psPiTorsionID2->_pSysData[i].y]++;

#undef DUMP_PARAMETERS
#define DUMP_PARAMETERS 5
#if (DUMP_PARAMETERS > 0 )
if( (i < DUMP_PARAMETERS || i > piTorsions - (DUMP_PARAMETERS + 1)) && amoebaGpu->log )
         fprintf( amoebaGpu->log, "PiTorsions: %5d [%5d %5d %5d %5d %5d %5d [%5d %5d %5d %5d %5d %5d]  k=%15.7e [%5d %5d %5d %5d %5d %5d]\n", i, 
                  (*psPiTorsionID1)[i].x, (*psPiTorsionID1)[i].y, (*psPiTorsionID1)[i].z, (*psPiTorsionID1)[i].w,
                  (*psPiTorsionID2)[i].x, (*psPiTorsionID2)[i].y, (*psPiTorsionID2)[i].z, (*psPiTorsionID2)[i].w,
                  (*psPiTorsionID3)[i].x, (*psPiTorsionID3)[i].y, (*psPiTorsionID3)[i].z, (*psPiTorsionID3)[i].w,
                  (*psPiTorsionParameter)[i],
                  gpu->pOutputBufferCounter[psPiTorsionID1->_pSysData[i].x],
                  gpu->pOutputBufferCounter[psPiTorsionID1->_pSysData[i].y],
                  gpu->pOutputBufferCounter[psPiTorsionID1->_pSysData[i].z],
                  gpu->pOutputBufferCounter[psPiTorsionID1->_pSysData[i].w],
                  gpu->pOutputBufferCounter[psPiTorsionID2->_pSysData[i].x],
                  gpu->pOutputBufferCounter[psPiTorsionID2->_pSysData[i].y] );
#endif
#undef DUMP_PARAMETERS
    }

    psPiTorsionID1->Upload();
    psPiTorsionID2->Upload();
    psPiTorsionID3->Upload();
    psPiTorsionParameter->Upload();
}

extern "C"
void gpuSetAmoebaStretchBendParameters(amoebaGpuContext amoebaGpu, const std::vector<int>& particles1, const std::vector<int>& particles2, const std::vector<int>& particles3,
                                       const std::vector<float>& lengthAB, const std::vector<float>& lengthCB,
                                       const std::vector<float>& angle, const std::vector<float>& k )
{

    _gpuContext* gpu                                  = amoebaGpu->gpuContext;
    int stretchBends                                  = particles1.size();
    amoebaGpu->amoebaSim.amoebaStretchBends           = stretchBends;

    CUDAStream<int4>* psStretchBendID1                = new CUDAStream<int4>(stretchBends, 1, "AmoebaStretchBendID1");
    amoebaGpu->psAmoebaStretchBendID1                 = psStretchBendID1;
    amoebaGpu->amoebaSim.pAmoebaStretchBendID1        = psStretchBendID1->_pDevStream[0];

    CUDAStream<int2>* psStretchBendID2                = new CUDAStream<int2>(stretchBends, 1, "AmoebaStretchBendID2");
    amoebaGpu->psAmoebaStretchBendID2                 = psStretchBendID2;
    amoebaGpu->amoebaSim.pAmoebaStretchBendID2        = psStretchBendID2->_pDevStream[0];

    CUDAStream<float4>* psStretchBendParameter        = new CUDAStream<float4>(stretchBends, 1, "AmoebaStretchBendParameter1");
    amoebaGpu->psAmoebaStretchBendParameter           = psStretchBendParameter;
    amoebaGpu->amoebaSim.pAmoebaStretchBendParameter  = psStretchBendParameter->_pDevStream[0];        

    for (int i = 0; i < stretchBends; i++)
    {
        (*psStretchBendID1)[i].x         = particles1[i];
        (*psStretchBendID1)[i].y         = particles2[i];
        (*psStretchBendID1)[i].z         = particles3[i];
        (*psStretchBendParameter)[i].x   = lengthAB[i];
        (*psStretchBendParameter)[i].y   = lengthCB[i];
        (*psStretchBendParameter)[i].z   = angle[i];
        (*psStretchBendParameter)[i].w   = k[i];
        psStretchBendID1->_pSysData[i].w = gpu->pOutputBufferCounter[psStretchBendID1->_pSysData[i].x]++;
        psStretchBendID2->_pSysData[i].x = gpu->pOutputBufferCounter[psStretchBendID1->_pSysData[i].y]++;
        psStretchBendID2->_pSysData[i].y = gpu->pOutputBufferCounter[psStretchBendID1->_pSysData[i].z]++;

#undef DUMP_PARAMETERS
#define DUMP_PARAMETERS 5
#if (DUMP_PARAMETERS > 0 )
if( (i < DUMP_PARAMETERS || i > stretchBends - (DUMP_PARAMETERS + 1)) && amoebaGpu->log )
         fprintf( amoebaGpu->log, "StretchBends: %5d [%5d %5d %5d] [%5d %5d %5d] [%15.7e %15.7e %15.7e %15.7e [%5d %5d %5d]\n", i, 
                  (*psStretchBendID1)[i].x, (*psStretchBendID1)[i].y, (*psStretchBendID1)[i].z, (*psStretchBendID1)[i].w,
                  (*psStretchBendID2)[i].x, (*psStretchBendID2)[i].y,
                  (*psStretchBendParameter)[i].x, (*psStretchBendParameter)[i].y, (*psStretchBendParameter)[i].z, (*psStretchBendParameter)[i].w,
                  gpu->pOutputBufferCounter[psStretchBendID1->_pSysData[i].x],
                  gpu->pOutputBufferCounter[psStretchBendID1->_pSysData[i].y],
                  gpu->pOutputBufferCounter[psStretchBendID1->_pSysData[i].z] );
#endif
#undef DUMP_PARAMETERS
    }
    psStretchBendID1->Upload();
    psStretchBendID2->Upload();
    psStretchBendParameter->Upload();
}

extern "C"
void gpuSetAmoebaOutOfPlaneBendParameters(amoebaGpuContext amoebaGpu, const std::vector<int>& particles1, const std::vector<int>& particles2, const std::vector<int>& particles3,
                                          const std::vector<int>& particles4, const std::vector<float>& k,
                                          float cubicK, float quarticK, float penticK, float sexticK  )
{

    _gpuContext* gpu                                     = amoebaGpu->gpuContext;
    int outOfPlaneBends                                  = particles1.size();
    amoebaGpu->amoebaSim.amoebaOutOfPlaneBends           = outOfPlaneBends;

    CUDAStream<int4>* psOutOfPlaneBendID1                = new CUDAStream<int4>(outOfPlaneBends, 1, "AmoebaOutOfPlaneBendID1");
    amoebaGpu->psAmoebaOutOfPlaneBendID1                 = psOutOfPlaneBendID1;
    amoebaGpu->amoebaSim.pAmoebaOutOfPlaneBendID1        = psOutOfPlaneBendID1->_pDevStream[0];

    CUDAStream<int4>* psOutOfPlaneBendID2                = new CUDAStream<int4>(outOfPlaneBends, 1, "AmoebaOutOfPlaneBendID2");
    amoebaGpu->psAmoebaOutOfPlaneBendID2                 = psOutOfPlaneBendID2;
    amoebaGpu->amoebaSim.pAmoebaOutOfPlaneBendID2        = psOutOfPlaneBendID2->_pDevStream[0];

    CUDAStream<float>* psOutOfPlaneBendParameter         = new CUDAStream<float>(outOfPlaneBends, 1, "AmoebaOutOfPlaneBendParameter");
    amoebaGpu->psAmoebaOutOfPlaneBendParameter           = psOutOfPlaneBendParameter;
    amoebaGpu->amoebaSim.pAmoebaOutOfPlaneBendParameter  = psOutOfPlaneBendParameter->_pDevStream[0];        

    amoebaGpu->amoebaSim.amoebaOutOfPlaneBendCubicK      = cubicK;
    amoebaGpu->amoebaSim.amoebaOutOfPlaneBendQuarticK    = quarticK;
    amoebaGpu->amoebaSim.amoebaOutOfPlaneBendPenticK     = penticK;
    amoebaGpu->amoebaSim.amoebaOutOfPlaneBendSexticK     = sexticK;

#undef DUMP_PARAMETERS
#define DUMP_PARAMETERS 5
#if (DUMP_PARAMETERS > 0 )
    if( amoebaGpu->log )
        fprintf( amoebaGpu->log, "OutOfPlaneBends: global ks[%15.7e %15.7e %15.7e %15.7e]\n", cubicK, quarticK, penticK, sexticK );
#endif

    for (int i = 0; i < outOfPlaneBends; i++)
    {
        (*psOutOfPlaneBendID1)[i].x         = particles1[i];
        (*psOutOfPlaneBendID1)[i].y         = particles2[i];
        (*psOutOfPlaneBendID1)[i].z         = particles3[i];
        (*psOutOfPlaneBendID1)[i].w         = particles4[i];
        (*psOutOfPlaneBendParameter)[i]     = k[i];
        psOutOfPlaneBendID2->_pSysData[i].x = gpu->pOutputBufferCounter[psOutOfPlaneBendID1->_pSysData[i].x]++;
        psOutOfPlaneBendID2->_pSysData[i].y = gpu->pOutputBufferCounter[psOutOfPlaneBendID1->_pSysData[i].y]++;
        psOutOfPlaneBendID2->_pSysData[i].z = gpu->pOutputBufferCounter[psOutOfPlaneBendID1->_pSysData[i].z]++;
        psOutOfPlaneBendID2->_pSysData[i].w = gpu->pOutputBufferCounter[psOutOfPlaneBendID1->_pSysData[i].w]++;

#if (DUMP_PARAMETERS > 0 )
if( (i < DUMP_PARAMETERS || i > outOfPlaneBends - (DUMP_PARAMETERS + 1)) && amoebaGpu->log )
         fprintf( amoebaGpu->log, "OutOfPlaneBends: %5d [%5d %5d %5d %5d] [%5d %5d %5d %5d] k=%15.7e [%5d %5d %5d %5d]\n", i, 
                  (*psOutOfPlaneBendID1)[i].x, (*psOutOfPlaneBendID1)[i].y, (*psOutOfPlaneBendID1)[i].z, (*psOutOfPlaneBendID1)[i].w,
                  (*psOutOfPlaneBendID2)[i].x, (*psOutOfPlaneBendID2)[i].y, (*psOutOfPlaneBendID2)[i].z, (*psOutOfPlaneBendID2)[i].w,
                  (*psOutOfPlaneBendParameter)[i], 
                  gpu->pOutputBufferCounter[psOutOfPlaneBendID1->_pSysData[i].x],
                  gpu->pOutputBufferCounter[psOutOfPlaneBendID1->_pSysData[i].y],
                  gpu->pOutputBufferCounter[psOutOfPlaneBendID1->_pSysData[i].z],
                  gpu->pOutputBufferCounter[psOutOfPlaneBendID1->_pSysData[i].w] );
#endif
#undef DUMP_PARAMETERS
    }
    psOutOfPlaneBendID1->Upload();
    psOutOfPlaneBendID2->Upload();
    psOutOfPlaneBendParameter->Upload();
}

extern "C"
void gpuSetAmoebaTorsionTorsionParameters(amoebaGpuContext amoebaGpu, const std::vector<int>& particles1, const std::vector<int>& particles2, const std::vector<int>& particles3,
                                          const std::vector<int>& particles4, const std::vector<int>& particles5,  const std::vector<int>& chiralParticleIndex, const std::vector<int>& gridIndices ) 
{

    _gpuContext* gpu                                     = amoebaGpu->gpuContext;
    int torsionTorsions                                  = particles1.size();
    amoebaGpu->amoebaSim.amoebaTorsionTorsions           = torsionTorsions;

    CUDAStream<int4>* psTorsionTorsionID1                = new CUDAStream<int4>(torsionTorsions, 1, "AmoebaTorsionTorsionID1");
    amoebaGpu->psAmoebaTorsionTorsionID1                 = psTorsionTorsionID1;
    amoebaGpu->amoebaSim.pAmoebaTorsionTorsionID1        = psTorsionTorsionID1->_pDevStream[0];

    CUDAStream<int4>* psTorsionTorsionID2                = new CUDAStream<int4>(torsionTorsions, 1, "AmoebaTorsionTorsionID2");
    amoebaGpu->psAmoebaTorsionTorsionID2                 = psTorsionTorsionID2;
    amoebaGpu->amoebaSim.pAmoebaTorsionTorsionID2        = psTorsionTorsionID2->_pDevStream[0];

    CUDAStream<int4>* psTorsionTorsionID3                = new CUDAStream<int4>(torsionTorsions, 1, "AmoebaTorsionTorsionID3");
    amoebaGpu->psAmoebaTorsionTorsionID3                 = psTorsionTorsionID3;
    amoebaGpu->amoebaSim.pAmoebaTorsionTorsionID3        = psTorsionTorsionID3->_pDevStream[0];

    for (int i = 0; i < torsionTorsions; i++)
    {
        (*psTorsionTorsionID1)[i].x         = particles1[i];
        (*psTorsionTorsionID1)[i].y         = particles2[i];
        (*psTorsionTorsionID1)[i].z         = particles3[i];
        (*psTorsionTorsionID1)[i].w         = particles4[i];
        (*psTorsionTorsionID2)[i].x         = particles5[i];

        (*psTorsionTorsionID2)[i].y         = chiralParticleIndex[i];
        (*psTorsionTorsionID2)[i].z         = gridIndices[i];

        psTorsionTorsionID2->_pSysData[i].w = gpu->pOutputBufferCounter[psTorsionTorsionID1->_pSysData[i].x]++;
        psTorsionTorsionID3->_pSysData[i].x = gpu->pOutputBufferCounter[psTorsionTorsionID1->_pSysData[i].y]++;
        psTorsionTorsionID3->_pSysData[i].y = gpu->pOutputBufferCounter[psTorsionTorsionID1->_pSysData[i].z]++;
        psTorsionTorsionID3->_pSysData[i].z = gpu->pOutputBufferCounter[psTorsionTorsionID1->_pSysData[i].w]++;
        psTorsionTorsionID3->_pSysData[i].w = gpu->pOutputBufferCounter[psTorsionTorsionID2->_pSysData[i].x]++;

#undef DUMP_PARAMETERS
#define DUMP_PARAMETERS 5
#if (DUMP_PARAMETERS > 0 )
if( (i < DUMP_PARAMETERS || i > torsionTorsions - (DUMP_PARAMETERS + 1)) && amoebaGpu->log )
         fprintf( amoebaGpu->log, "TorsionTorsions: %5d [%5d %5d %5d %5d %5d] chiral=%5d grid=%5d [%5d %5d %5d %5d %5d] [%5d %5d %5d %5d %5d]\n", i, 
                  (*psTorsionTorsionID1)[i].x, (*psTorsionTorsionID1)[i].y, (*psTorsionTorsionID1)[i].z, (*psTorsionTorsionID1)[i].w,
                  (*psTorsionTorsionID2)[i].x, (*psTorsionTorsionID2)[i].y, (*psTorsionTorsionID2)[i].z, (*psTorsionTorsionID2)[i].w,
                  (*psTorsionTorsionID3)[i].x, (*psTorsionTorsionID3)[i].y, (*psTorsionTorsionID3)[i].z,  (*psTorsionTorsionID3)[i].w,
                  gpu->pOutputBufferCounter[psTorsionTorsionID1->_pSysData[i].x],
                  gpu->pOutputBufferCounter[psTorsionTorsionID1->_pSysData[i].y],
                  gpu->pOutputBufferCounter[psTorsionTorsionID1->_pSysData[i].z],
                  gpu->pOutputBufferCounter[psTorsionTorsionID1->_pSysData[i].w], 
                  gpu->pOutputBufferCounter[psTorsionTorsionID2->_pSysData[i].x] );
#endif
#undef DUMP_PARAMETERS
    }
    psTorsionTorsionID1->Upload();
    psTorsionTorsionID2->Upload();
    psTorsionTorsionID3->Upload();
}

static void testAmoebaTorsionTorsionGridLookup( amoebaGpuContext amoebaGpu, float* grids, int gridIndex, float angle1, float angle2, std::vector<float>& values )
{

    int index1    = static_cast<int>((angle1 - amoebaGpu->amoebaSim.amoebaTorTorGridBegin[gridIndex])/amoebaGpu->amoebaSim.amoebaTorTorGridDelta[gridIndex]);
    int index2    = static_cast<int>((angle2 - amoebaGpu->amoebaSim.amoebaTorTorGridBegin[gridIndex])/amoebaGpu->amoebaSim.amoebaTorTorGridDelta[gridIndex]);
    int index     = index2 + index1*amoebaGpu->amoebaSim.amoebaTorTorGridNy[gridIndex];
        index    *= 4;
        index    += amoebaGpu->amoebaSim.amoebaTorTorGridOffset[gridIndex];
    int i         = 0;
//(void) fprintf( amoebaGpu->log, "%d %d   %d  [%10.3f %10.3f] %d\n", index1, index2, index, angle1, angle2, gridIndex );
    values.resize( 4 );
    values[i++]   = grids[index++];
    values[i++]   = grids[index++];
    values[i++]   = grids[index++];
    values[i++]   = grids[index++];
}

extern "C"
void gpuSetAmoebaTorsionTorsionGrids(amoebaGpuContext amoebaGpu, const std::vector< std::vector< std::vector< std::vector<float> > > >& floatGrids )
{

    _gpuContext* gpu                                     = amoebaGpu->gpuContext;
    unsigned int torsionTorsionGrids                     = floatGrids.size();
    unsigned int totalGridEntries                        = 0;
    // 4 (grids) * (25 *25 grid)*(2 +4 a1, a2, f, f1,f2, f12) = 15000
    for (unsigned int ii = 0; ii < floatGrids.size(); ii++) {
        amoebaGpu->amoebaSim.amoebaTorTorGridOffset[ii] = (totalGridEntries/4);
        amoebaGpu->amoebaSim.amoebaTorTorGridBegin[ii]  = -180.0f;
        amoebaGpu->amoebaSim.amoebaTorTorGridDelta[ii]  =   15.0f;
        amoebaGpu->amoebaSim.amoebaTorTorGridNy[ii]     =   25;
        for (unsigned int jj = 0; jj < floatGrids[ii].size(); jj++) {
            for (unsigned int kk = 0; kk < floatGrids[ii][jj].size(); kk++) {
                totalGridEntries += (floatGrids[ii][jj][kk].size() - 2);
            }
        }
    }

    unsigned int totalEntries                            = totalGridEntries/4;
    CUDAStream<float4>* psTorsionTorsionGrids            = new CUDAStream<float4>(totalEntries, 1, "AmoebaTorsionTorsionGrids");
    amoebaGpu->psAmoebaTorsionTorsionGrids               = psTorsionTorsionGrids;
    amoebaGpu->amoebaSim.pAmoebaTorsionTorsionGrids      = psTorsionTorsionGrids->_pDevStream[0];

    if( amoebaGpu->log ){
        (void) fprintf( amoebaGpu->log, "totalGridEntries=%u totalFloat4 entries=%u\n", totalGridEntries, totalEntries );
    }

    unsigned int index    = 0;
    for (unsigned int ii = 0; ii < floatGrids.size(); ii++) {
        for (unsigned int jj = 0; jj < floatGrids[ii].size(); jj++) {
            for (unsigned int kk = 0; kk < floatGrids[ii][jj].size(); kk++) {
                (*psTorsionTorsionGrids)[index].x    = floatGrids[ii][jj][kk][2];
                (*psTorsionTorsionGrids)[index].y    = floatGrids[ii][jj][kk][3];
                (*psTorsionTorsionGrids)[index].z    = floatGrids[ii][jj][kk][4];
                (*psTorsionTorsionGrids)[index].w    = floatGrids[ii][jj][kk][5];

#undef DUMP_PARAMETERS
#define DUMP_PARAMETERS 5
#if (DUMP_PARAMETERS > 0 )
if( (index < DUMP_PARAMETERS || index > totalEntries - (DUMP_PARAMETERS + 1)) && amoebaGpu->log )
         fprintf( amoebaGpu->log, "TorsionTorsionGrid: %d %5d [%5d %5d ] [%10.3f %10.3f] [%15.7e %15.7e %15.7e %15.7e]\n", index, ii, jj, kk,
                  floatGrids[ii][jj][kk][0], floatGrids[ii][jj][kk][1],
                  (*psTorsionTorsionGrids)[index].x, (*psTorsionTorsionGrids)[index].y, (*psTorsionTorsionGrids)[index].z, (*psTorsionTorsionGrids)[index].w );
#endif
#undef DUMP_PARAMETERS
                index++;
            }
        }
    }
    
#if 0
    float* grids = (float*) malloc( totalGridEntries*sizeof( float ) );
    int index    = 0;
    for (unsigned int ii = 0; ii < floatGrids.size(); ii++) {
        for (unsigned int jj = 0; jj < floatGrids[ii].size(); jj++) {
            for (unsigned int kk = 0; kk < floatGrids[ii][jj].size(); kk++) {
                for (unsigned int mm = 2; mm < floatGrids[ii][jj][kk].size(); mm++) {
                    grids[index++] = floatGrids[ii][jj][kk][mm]; 
                }
            }
        }
    }
    
    unsigned int i = 0;
    for (unsigned int ii = 0; ii < ; ii++) {
        for (unsigned int jj = 0; jj < floatGrids[ii].size(); jj++) {
            for (unsigned int kk = 0; kk < floatGrids[ii][jj].size(); kk++) {
                (*psTorsionTorsionGrids)[i].x    = floatGrids[ii][jj][kk][2];
                (*psTorsionTorsionGrids)[i].y    = floatGrids[ii][jj][kk][3];
                (*psTorsionTorsionGrids)[i].z    = floatGrids[ii][jj][kk][4];
                (*psTorsionTorsionGrids)[i].w    = floatGrids[ii][jj][kk][5];
                i++;
            }
        }
    }

    float epsilon = 1.0e-06;
    int errors    = 0;
    for (unsigned int ii = 0; ii < floatGrids.size(); ii++) {
        for (unsigned int jj = 0; jj < floatGrids[ii].size(); jj++) {
            for (unsigned int kk = 0; kk < floatGrids[ii][jj].size(); kk++) {
                    std::vector<float> values;
                    testAmoebaTorsionTorsionGridLookup( amoebaGpu, grids, ii, floatGrids[ii][jj][kk][0] + 1.0f, floatGrids[ii][jj][kk][1]+ 1.0f, values );
                    if( fabsf( values[0] - floatGrids[ii][jj][kk][2] ) > epsilon ||
                        fabsf( values[1] - floatGrids[ii][jj][kk][3] ) > epsilon ||
                        fabsf( values[2] - floatGrids[ii][jj][kk][4] ) > epsilon ||
                        fabsf( values[3] - floatGrids[ii][jj][kk][5] ) > epsilon ){
                        (void) fprintf( amoebaGpu->log, "Error %u %u %u [%10.3f %10.3f]  [%15.7e %15.7e %15.7e %15.7e] lk=[%15.7e %15.7e %15.7e %15.7e]\n", ii, jj, kk,
                                        floatGrids[ii][jj][kk][0], floatGrids[ii][jj][kk][1],
                                        floatGrids[ii][jj][kk][2], floatGrids[ii][jj][kk][3], floatGrids[ii][jj][kk][4], floatGrids[ii][jj][kk][5],
                                        values[0], values[1], values[2], values[3] );
                        if( errors++ > 10 ){
                           (void) fflush( amoebaGpu->log );
                           exit(0);
                        }
                    } 
            }
        }
    }
    if( !errors ){
        (void) fprintf( amoebaGpu->log, "No errors in grid readback\n" );
    }
    free( grids );
exit(0);
#endif



    psTorsionTorsionGrids->Upload();

}

extern "C"
void gpuSetAmoebaBondOffsets(amoebaGpuContext amoebaGpu )
{
 
    // make sure only flip once

    static int flipped = 0;
    if( amoebaGpu && flipped ){
        return;
    }
    flipped = 1;
    _gpuContext* gpu                                           = amoebaGpu->gpuContext;


    amoebaGpu->amoebaSim.amoebaBond_offset                     = amoebaGpu->psAmoebaBondParameter ? amoebaGpu->psAmoebaBondParameter->_stride : 0;

    amoebaGpu->amoebaSim.amoebaAngle_offset                    = amoebaGpu->amoebaSim.amoebaBond_offset            + 
                                                                 (amoebaGpu->psAmoebaAngleParameter ? amoebaGpu->psAmoebaAngleParameter->_stride : 0);

    amoebaGpu->amoebaSim.amoebaInPlaneAngle_offset             = amoebaGpu->amoebaSim.amoebaAngle_offset           +
                                                                 (amoebaGpu->psAmoebaInPlaneAngleParameter ? amoebaGpu->psAmoebaInPlaneAngleParameter->_stride : 0);

    amoebaGpu->amoebaSim.amoebaTorsion_offset                  = amoebaGpu->amoebaSim.amoebaInPlaneAngle_offset    +
                                                                 (amoebaGpu->psAmoebaTorsionParameter1 ?  amoebaGpu->psAmoebaTorsionParameter1->_stride : 0);

    amoebaGpu->amoebaSim.amoebaPiTorsion_offset                = amoebaGpu->amoebaSim.amoebaTorsion_offset         +
                                                                 (amoebaGpu->psAmoebaPiTorsionParameter ? amoebaGpu->psAmoebaPiTorsionParameter->_stride : 0);

    amoebaGpu->amoebaSim.amoebaStretchBend_offset              = amoebaGpu->amoebaSim.amoebaPiTorsion_offset       +
                                                                 (amoebaGpu->psAmoebaStretchBendParameter ? amoebaGpu->psAmoebaStretchBendParameter->_stride : 0);

    amoebaGpu->amoebaSim.amoebaOutOfPlaneBend_offset           = amoebaGpu->amoebaSim.amoebaStretchBend_offset     +
                                                                 (amoebaGpu->psAmoebaOutOfPlaneBendParameter ? amoebaGpu->psAmoebaOutOfPlaneBendParameter->_stride : 0);

    amoebaGpu->amoebaSim.amoebaTorsionTorsion_offset           = amoebaGpu->amoebaSim.amoebaOutOfPlaneBend_offset  +
                                                                 (amoebaGpu->psAmoebaTorsionTorsionID1 ?  amoebaGpu->psAmoebaTorsionTorsionID1->_stride : 0);

    //gpu->sim.localForces_threads_per_block  = (std::max(amoebaGpu->amoebaSim.amoebaTorsionTorsion_offset, gpu->sim.customBonds) / gpu->sim.blocks + 15) & 0xfffffff0;
    unsigned int maxI                                         = (amoebaGpu->amoebaSim.amoebaTorsionTorsion_offset > gpu->sim.customBonds) ? amoebaGpu->amoebaSim.amoebaTorsionTorsion_offset : gpu->sim.customBonds;
    gpu->sim.localForces_threads_per_block                    = (maxI/gpu->sim.blocks + 15) & 0xfffffff0;
    if (gpu->sim.localForces_threads_per_block > gpu->sim.max_localForces_threads_per_block)
        gpu->sim.localForces_threads_per_block = gpu->sim.max_localForces_threads_per_block;
    if (gpu->sim.localForces_threads_per_block < 1) 
        gpu->sim.localForces_threads_per_block = 1; 

    // Flip local force output buffers

    int flip = gpu->sim.outputBuffers - 1; 
    if( flip > 0 ){
        for (unsigned int i = 0; i < amoebaGpu->amoebaSim.amoebaBonds; i++) 
        {    
            (*amoebaGpu->psAmoebaBondID)[i].z = flip - (*amoebaGpu->psAmoebaBondID)[i].z;
            (*amoebaGpu->psAmoebaBondID)[i].w = flip - (*amoebaGpu->psAmoebaBondID)[i].w;
        }    
        if( amoebaGpu->psAmoebaBondID )
            amoebaGpu->psAmoebaBondID->Upload();
    
        for (unsigned int i = 0; i < amoebaGpu->amoebaSim.amoebaAngles; i++) 
        {    
            (*amoebaGpu->psAmoebaAngleID1)[i].w = flip - (*amoebaGpu->psAmoebaAngleID1)[i].w;
            (*amoebaGpu->psAmoebaAngleID2)[i].x = flip - (*amoebaGpu->psAmoebaAngleID2)[i].x;
            (*amoebaGpu->psAmoebaAngleID2)[i].y = flip - (*amoebaGpu->psAmoebaAngleID2)[i].y;
        }    
        if( amoebaGpu->psAmoebaAngleID1 && amoebaGpu->psAmoebaAngleID2 ){
            amoebaGpu->psAmoebaAngleID1->Upload();
            amoebaGpu->psAmoebaAngleID2->Upload();
        }

        for (unsigned int i = 0; i < amoebaGpu->amoebaSim.amoebaInPlaneAngles; i++) 
        {    
            (*amoebaGpu->psAmoebaInPlaneAngleID2)[i].x = flip - (*amoebaGpu->psAmoebaInPlaneAngleID2)[i].x;
            (*amoebaGpu->psAmoebaInPlaneAngleID2)[i].y = flip - (*amoebaGpu->psAmoebaInPlaneAngleID2)[i].y;
            (*amoebaGpu->psAmoebaInPlaneAngleID2)[i].z = flip - (*amoebaGpu->psAmoebaInPlaneAngleID2)[i].z;
            (*amoebaGpu->psAmoebaInPlaneAngleID2)[i].w = flip - (*amoebaGpu->psAmoebaInPlaneAngleID2)[i].w;
        }    
        if( amoebaGpu->psAmoebaInPlaneAngleID2 ){
            amoebaGpu->psAmoebaInPlaneAngleID2->Upload();
        }

        for (unsigned int i = 0; i < amoebaGpu->amoebaSim.amoebaTorsions; i++) 
        {    
            (*amoebaGpu->psAmoebaTorsionID2)[i].x = flip - (*amoebaGpu->psAmoebaTorsionID2)[i].x;
            (*amoebaGpu->psAmoebaTorsionID2)[i].y = flip - (*amoebaGpu->psAmoebaTorsionID2)[i].y;
            (*amoebaGpu->psAmoebaTorsionID2)[i].z = flip - (*amoebaGpu->psAmoebaTorsionID2)[i].z;
            (*amoebaGpu->psAmoebaTorsionID2)[i].w = flip - (*amoebaGpu->psAmoebaTorsionID2)[i].w;
        }    
        if( amoebaGpu->psAmoebaTorsionID2 )
            amoebaGpu->psAmoebaTorsionID2->Upload();

        for (unsigned int i = 0; i < amoebaGpu->amoebaSim.amoebaPiTorsions; i++) 
        {    
            (*amoebaGpu->psAmoebaPiTorsionID2)[i].z = flip - (*amoebaGpu->psAmoebaPiTorsionID2)[i].z;
            (*amoebaGpu->psAmoebaPiTorsionID2)[i].w = flip - (*amoebaGpu->psAmoebaPiTorsionID2)[i].w;
            (*amoebaGpu->psAmoebaPiTorsionID3)[i].x = flip - (*amoebaGpu->psAmoebaPiTorsionID3)[i].x;
            (*amoebaGpu->psAmoebaPiTorsionID3)[i].y = flip - (*amoebaGpu->psAmoebaPiTorsionID3)[i].y;
            (*amoebaGpu->psAmoebaPiTorsionID3)[i].z = flip - (*amoebaGpu->psAmoebaPiTorsionID3)[i].z;
            (*amoebaGpu->psAmoebaPiTorsionID3)[i].w = flip - (*amoebaGpu->psAmoebaPiTorsionID3)[i].w;
        }    
        if( amoebaGpu->psAmoebaPiTorsionID2 && amoebaGpu->psAmoebaPiTorsionID3 ){
            amoebaGpu->psAmoebaPiTorsionID2->Upload();
            amoebaGpu->psAmoebaPiTorsionID3->Upload();
        }

        for (unsigned int i = 0; i < amoebaGpu->amoebaSim.amoebaStretchBends; i++) 
        {    
            (*amoebaGpu->psAmoebaStretchBendID1)[i].w = flip - (*amoebaGpu->psAmoebaStretchBendID1)[i].w;
            (*amoebaGpu->psAmoebaStretchBendID2)[i].x = flip - (*amoebaGpu->psAmoebaStretchBendID2)[i].x;
            (*amoebaGpu->psAmoebaStretchBendID2)[i].y = flip - (*amoebaGpu->psAmoebaStretchBendID2)[i].y;
        }    
        if( amoebaGpu->psAmoebaStretchBendID1 && amoebaGpu->psAmoebaStretchBendID2 ){
            amoebaGpu->psAmoebaStretchBendID1->Upload();
            amoebaGpu->psAmoebaStretchBendID2->Upload();
        }

        for (unsigned int i = 0; i < amoebaGpu->amoebaSim.amoebaOutOfPlaneBends; i++) 
        {    
            (*amoebaGpu->psAmoebaOutOfPlaneBendID2)[i].x = flip - (*amoebaGpu->psAmoebaOutOfPlaneBendID2)[i].x;
            (*amoebaGpu->psAmoebaOutOfPlaneBendID2)[i].y = flip - (*amoebaGpu->psAmoebaOutOfPlaneBendID2)[i].y;
            (*amoebaGpu->psAmoebaOutOfPlaneBendID2)[i].z = flip - (*amoebaGpu->psAmoebaOutOfPlaneBendID2)[i].z;
            (*amoebaGpu->psAmoebaOutOfPlaneBendID2)[i].w = flip - (*amoebaGpu->psAmoebaOutOfPlaneBendID2)[i].w;
        }    
        if( amoebaGpu->psAmoebaOutOfPlaneBendID2 ){
            amoebaGpu->psAmoebaOutOfPlaneBendID2->Upload();
        }

        for (unsigned int i = 0; i < amoebaGpu->amoebaSim.amoebaTorsionTorsions; i++) 
        {    
            (*amoebaGpu->psAmoebaTorsionTorsionID2)[i].w = flip - (*amoebaGpu->psAmoebaTorsionTorsionID2)[i].w;
            (*amoebaGpu->psAmoebaTorsionTorsionID3)[i].x = flip - (*amoebaGpu->psAmoebaTorsionTorsionID3)[i].x;
            (*amoebaGpu->psAmoebaTorsionTorsionID3)[i].y = flip - (*amoebaGpu->psAmoebaTorsionTorsionID3)[i].y;
            (*amoebaGpu->psAmoebaTorsionTorsionID3)[i].z = flip - (*amoebaGpu->psAmoebaTorsionTorsionID3)[i].z;
            (*amoebaGpu->psAmoebaTorsionTorsionID3)[i].w = flip - (*amoebaGpu->psAmoebaTorsionTorsionID3)[i].w;
        }    
        if( amoebaGpu->psAmoebaTorsionTorsionID2 && amoebaGpu->psAmoebaTorsionTorsionID3 ){
            amoebaGpu->psAmoebaTorsionTorsionID2->Upload();
            amoebaGpu->psAmoebaTorsionTorsionID3->Upload();
        }

    }    
    

}

/**---------------------------------------------------------------------------------------

   Allocate data structs associated w/ molecular -> lab frame calculation

   @param amoebaGpu        amoebaGpu context

   --------------------------------------------------------------------------------------- */

static void gpuRotationToLabFrameAllocate( amoebaGpuContext amoebaGpu )
{
    // ---------------------------------------------------------------------------------------

    static const std::string methodName = "gpuRotationToLabFrameAllocate";

    // ---------------------------------------------------------------------------------------

    if( amoebaGpu->psRotationMatrix != NULL ){
        return;
    }
   
#ifdef AMOEBA_DEBUG
    if( amoebaGpu->log ){
        (void) fprintf( amoebaGpu->log,"%s: paddedNumberOfAtoms=%d\n",
                        methodName.c_str(), amoebaGpu->paddedNumberOfAtoms ); (void) fflush( amoebaGpu->log );
    }
#endif

    // work space
 
    amoebaGpu->psRotationMatrix                      = new CUDAStream<float>(9*amoebaGpu->paddedNumberOfAtoms, 1, "RotationMatrix");
 
    // parameters
 
    amoebaGpu->psMultipoleParticlesIdsAndAxisType    = new CUDAStream<int4>(amoebaGpu->paddedNumberOfAtoms,    1, "MultipoleParticlesIdsAndAxisType");
    amoebaGpu->psMolecularDipole                     = new CUDAStream<float>(3*amoebaGpu->paddedNumberOfAtoms, 1, "MolecularDipole");
    amoebaGpu->psMolecularQuadrupole                 = new CUDAStream<float>(9*amoebaGpu->paddedNumberOfAtoms, 1, "MolecularQuadrupole");
 
    // output
 
    amoebaGpu->psLabFrameDipole                      = new CUDAStream<float>(3*amoebaGpu->paddedNumberOfAtoms, 1, "LabFrameDipole");
    amoebaGpu->amoebaSim.pLabFrameDipole             = amoebaGpu->psLabFrameDipole->_pDevStream[0];

    amoebaGpu->psLabFrameQuadrupole                  = new CUDAStream<float>(9*amoebaGpu->paddedNumberOfAtoms, 1, "LabFrameQuadrupole");
    amoebaGpu->amoebaSim.pLabFrameQuadrupole         = amoebaGpu->psLabFrameQuadrupole->_pDevStream[0];

    memset( amoebaGpu->psMultipoleParticlesIdsAndAxisType->_pSysStream[0], 0, sizeof(int)*4*amoebaGpu->paddedNumberOfAtoms );
}

static void gpuFixedEFieldAllocate( amoebaGpuContext amoebaGpu )
{
    // ---------------------------------------------------------------------------------------

    static const std::string methodName = "gpuFixedEFieldAllocate";

    // ---------------------------------------------------------------------------------------

    if( amoebaGpu->psE_Field != NULL ){
        return;
    }

    int paddedNumberOfAtoms                    = amoebaGpu->paddedNumberOfAtoms;

#ifdef AMOEBA_DEBUG
    if( amoebaGpu->log ){
        (void) fprintf( amoebaGpu->log,"%s: paddedNumberOfAtoms=%d maxCovalentDegreeSz=%d\n", methodName.c_str(),
                        paddedNumberOfAtoms, amoebaGpu->maxCovalentDegreeSz );
        (void) fflush( amoebaGpu->log );
    }    
#endif

    amoebaGpu->psE_Field                             = new CUDAStream<float>(paddedNumberOfAtoms*3, 1, "E_Field");
    amoebaGpu->psE_FieldPolar                        = new CUDAStream<float>(paddedNumberOfAtoms*3, 1, "E_FieldPolar");

    // parameters

    amoebaGpu->psDampingFactorAndThole               = new CUDAStream<float2>(paddedNumberOfAtoms, 1, "DampingFactorAndThole");
    amoebaGpu->amoebaSim.pDampingFactorAndThole      = amoebaGpu->psDampingFactorAndThole->_pDevStream[0];

    amoebaGpu->psCovalentDegree                      = new CUDAStream<int>(amoebaGpu->maxCovalentDegreeSz*paddedNumberOfAtoms, 1, "CovalentDegree");
    amoebaGpu->psPolarizationDegree                  = new CUDAStream<int>(amoebaGpu->maxCovalentDegreeSz*paddedNumberOfAtoms, 1, "PolarizationDegree");
    
    unsigned int offset                              = paddedNumberOfAtoms*sizeof( float );
    memset( amoebaGpu->psDampingFactorAndThole->_pSysStream[0],              0,2*offset );
    //memset( amoebaGpu->psE_Field->_pSysStream[0],            0, offset*3 );
    //memset( amoebaGpu->psE_FieldPolar->_pSysStream[0],       0, offset*3 );

    // should be removed XXXXX

    offset                                           = amoebaGpu->maxCovalentDegreeSz*paddedNumberOfAtoms*sizeof( int );
    memset( amoebaGpu->psCovalentDegree->_pSysStream[0],     0, offset );
    memset( amoebaGpu->psPolarizationDegree->_pSysStream[0], 0, offset );

}

/**---------------------------------------------------------------------------------------

   Allocate data structs associated w/ computing mutual induced field

   @param amoebaGpu        amoebaGpu context

   --------------------------------------------------------------------------------------- */

extern "C"
void gpuMutualInducedFieldAllocate( amoebaGpuContext amoebaGpu )
{
    // ---------------------------------------------------------------------------------------

    static const std::string methodName = "gpuMutualInducedFieldAllocate";

    // ---------------------------------------------------------------------------------------

    int paddedNumberOfAtoms                    = amoebaGpu->paddedNumberOfAtoms;

#ifdef AMOEBA_DEBUG
    if( amoebaGpu->log ){
        (void) fprintf( amoebaGpu->log,"%s: paddedNumberOfAtoms=%d\n", methodName.c_str(), paddedNumberOfAtoms );
        (void) fflush( amoebaGpu->log );
    }
#endif

    amoebaGpu->psInducedDipole                 = new CUDAStream<float>(paddedNumberOfAtoms*3, 1, "InducedDipole");
    amoebaGpu->amoebaSim.pInducedDipole        = amoebaGpu->psInducedDipole->_pDevStream[0];

    amoebaGpu->psInducedDipolePolar            = new CUDAStream<float>(paddedNumberOfAtoms*3, 1, "InducedDipolePolar");
    amoebaGpu->amoebaSim.pInducedDipolePolar   = amoebaGpu->psInducedDipolePolar->_pDevStream[0];

    amoebaGpu->psCurrentEpsilon                = new CUDAStream<float>(5, 1, "CurrentEpsilon");
    amoebaGpu->epsilonThreadsPerBlock          = 384;
 
    amoebaGpu->psPolarizability                = new CUDAStream<float>(paddedNumberOfAtoms*3, 1, "Polarizability");
    unsigned int offset                        = paddedNumberOfAtoms*3*sizeof( float );

    // currently only SOR

    if( amoebaGpu->mutualInducedIterativeMethod ==  0 ){

        for( unsigned int ii = 0; ii < amoebaGpu->numberOfSorWorkVectors; ii++ ){
            amoebaGpu->psWorkVector[ii]  = new CUDAStream<float>( paddedNumberOfAtoms*3, 1, "SorWorkVector"  );
        }
    }
 
    memset( amoebaGpu->psInducedDipole->_pSysStream[0],      0, offset );
    memset( amoebaGpu->psInducedDipolePolar->_pSysStream[0], 0, offset );
    memset( amoebaGpu->psPolarizability->_pSysStream[0],     0, offset );

}

/**---------------------------------------------------------------------------------------

   Allocate data structs associated w/ computing electrostatic  force/torque

   @param amoebaGpu        amoebaGpu context

   --------------------------------------------------------------------------------------- */

extern "C"
void gpuElectrostaticAllocate( amoebaGpuContext amoebaGpu )
{
    // ---------------------------------------------------------------------------------------

    static const std::string methodName = "gpuElectrostaticAllocate";

    // ---------------------------------------------------------------------------------------

    if( amoebaGpu->psForce != NULL ){
        return;
    }
    int paddedNumberOfAtoms                    = amoebaGpu->paddedNumberOfAtoms;

#ifdef AMOEBA_DEBUG
    if( amoebaGpu->log ){
        (void) fprintf( amoebaGpu->log,"%s: paddedNumberOfAtoms=%d\n", methodName.c_str(), paddedNumberOfAtoms );
        (void) fflush( amoebaGpu->log );
    }
#endif

    amoebaGpu->psForce                         = new CUDAStream<float>(paddedNumberOfAtoms*3, 1, "ElectrostaticForce");
    amoebaGpu->psTorque                        = new CUDAStream<float>(paddedNumberOfAtoms*3, 1, "Torque");

    unsigned int offset                        = 3*paddedNumberOfAtoms*sizeof( float );
    memset( amoebaGpu->psForce->_pSysStream[0],          0, offset );
    memset( amoebaGpu->psTorque->_pSysStream[0],         0, offset );

    amoebaGpu->psForce->Download();
    amoebaGpu->psTorque->Download();
    
}

/**---------------------------------------------------------------------------------------

   Allocate data structs associated w/ computing Kirkwood force/torque

   @param amoebaGpu        amoebaGpu context

   --------------------------------------------------------------------------------------- */

extern "C"
void gpuKirkwoodAllocate( amoebaGpuContext amoebaGpu )
{
    // ---------------------------------------------------------------------------------------

    static const std::string methodName = "gpuKirkwoodAllocate";

    // ---------------------------------------------------------------------------------------

    if( amoebaGpu->psBorn != NULL ){
        return;
    }

    int paddedNumberOfAtoms                      = amoebaGpu->paddedNumberOfAtoms;

#ifdef AMOEBA_DEBUG
    if( amoebaGpu->log ){
        (void) fprintf( amoebaGpu->log,"%s: paddedNumberOfAtoms      =%d\n", methodName.c_str(), paddedNumberOfAtoms );
        (void) fflush( amoebaGpu->log );
    }
#endif

    amoebaGpu->psBorn                            = new CUDAStream<float>(paddedNumberOfAtoms,   1, "KirkwoodBorn");
    amoebaGpu->psBornPolar                       = new CUDAStream<float>(paddedNumberOfAtoms,   1, "KirkwoodBornPolar");
    amoebaGpu->psGk_Field                        = new CUDAStream<float>(paddedNumberOfAtoms*3, 1, "Gk_Fixed_Field");

    amoebaGpu->psInducedDipoleS                  = new CUDAStream<float>(paddedNumberOfAtoms*3, 1, "InducedDipoleS");
    amoebaGpu->amoebaSim.pInducedDipoleS         = amoebaGpu->psInducedDipoleS->_pDevStream[0];

    amoebaGpu->psInducedDipolePolarS             = new CUDAStream<float>(paddedNumberOfAtoms*3, 1, "InducedDipolePolarS");
    amoebaGpu->amoebaSim.pInducedDipolePolarS    = amoebaGpu->psInducedDipolePolarS->_pDevStream[0];

    amoebaGpu->psKirkwoodForce                   = new CUDAStream<float>(paddedNumberOfAtoms*3, 1, "KirkwoodForce");
    amoebaGpu->psKirkwoodEDiffForce              = new CUDAStream<float>(paddedNumberOfAtoms*3, 1, "KirkwoodEDiffForce");

    unsigned int offset                    = paddedNumberOfAtoms*sizeof( float );
    memset( amoebaGpu->psBorn->_pSysStream[0],               0, offset );
    memset( amoebaGpu->psBornPolar->_pSysStream[0],          0, offset );

    amoebaGpu->psBorn->Download();
    amoebaGpu->psBornPolar->Download();
    
}

/**---------------------------------------------------------------------------------------

   Create/initialize data structs associated w/ molecular -> lab frame calculation

   @param amoebaGpu        amoebaGpu context

   --------------------------------------------------------------------------------------- */

extern "C" 
void gpuSetAmoebaMultipoleParameters(amoebaGpuContext amoebaGpu, const std::vector<float>& charges, const std::vector<float>& dipoles, const std::vector<float>& quadrupoles,
                                     const std::vector<int>& axisType, const std::vector<int>& multipoleParticleId1, const std::vector<int>& multipoleParticleId2,
                                     const std::vector<float>& tholes, float scalingDistanceCutoff,const std::vector<float>& dampingFactors, const std::vector<float>& polarity,
                                     const std::vector< std::vector< std::vector<int> > >& multipoleParticleCovalentInfo, const std::vector<int>& covalentDegree,
                                     const std::vector<int>& minCovalentIndices,  const std::vector<int>& minCovalentPolarizationIndices, int maxCovalentRange, 
                                     int mutualInducedIterativeMethod, int mutualInducedMaxIterations, float mutualInducedTargetEpsilon, float electricConstant ){

// ---------------------------------------------------------------------------------------

    static const std::string methodName     = "gpuSetAmoebaMultipoleParameters";
    int errorCount                          = 0;
    std::stringstream message;

// ---------------------------------------------------------------------------------------

    // initialize data structs associated w/ molecular -> lab frame calculation

    amoebaGpu->maxCovalentDegreeSz                                            = maxCovalentRange;
    amoebaGpu->paddedNumberOfAtoms                                            = amoebaGpu->gpuContext->sim.paddedNumberOfAtoms;
    amoebaGpu->scalingDistanceCutoff                                          = static_cast<float>(scalingDistanceCutoff);

    gpuRotationToLabFrameAllocate( amoebaGpu );
    gpuFixedEFieldAllocate( amoebaGpu );

    std::string minId;
    if( mutualInducedIterativeMethod == 0 ){
        minId = "SOR";
    } else if( mutualInducedIterativeMethod == 1 ){
        minId = "ConjugateGradient";
    }

#ifdef AMOEBA_DEBUG
    static unsigned int targetAtoms[2] = { 0, 700000 };

    float* pScaleCheckSum = (float*) malloc( sizeof( float )*charges.size() );
    float* dScaleCheckSum = (float*) malloc( sizeof( float )*charges.size() );
    float* mScaleCheckSum = (float*) malloc( sizeof( float )*charges.size() );

    memset( pScaleCheckSum, 0, charges.size()*sizeof( float ) );
    memset( dScaleCheckSum, 0, charges.size()*sizeof( float ) );
    memset( mScaleCheckSum, 0, charges.size()*sizeof( float ) );

    amoebaGpu->pMapArray = (MapIntFloat**) malloc( sizeof( MapIntFloat* )*amoebaGpu->paddedNumberOfAtoms );
    amoebaGpu->dMapArray = (MapIntFloat**) malloc( sizeof( MapIntFloat* )*amoebaGpu->paddedNumberOfAtoms );
    for( unsigned int ii = 0; ii < amoebaGpu->paddedNumberOfAtoms; ii++ ){
        amoebaGpu->pMapArray[ii]         =  new MapIntFloat;
        (*amoebaGpu->pMapArray[ii])[ii]  =  0.0f;

     //   amoebaGpu->mMapArray[ii]         =  new MapIntFloat;
     //   (*amoebaGpu->mMapArray[ii])[ii]  =  0.0f;

        amoebaGpu->dMapArray[ii]         =  new MapIntFloat;
        (*amoebaGpu->dMapArray[ii])[ii]  =  0.0f;
        for( unsigned int jj = amoebaGpu->gpuContext->natoms; jj < amoebaGpu->paddedNumberOfAtoms; jj++ ){
            (*amoebaGpu->pMapArray[ii])[jj]  =  0.0f;
            (*amoebaGpu->dMapArray[ii])[jj]  =  0.0f;
      //      (*amoebaGpu->mMapArray[ii])[jj]  =  0.0f;
        }
    }

    if( amoebaGpu->log ){
        (void) fprintf( amoebaGpu->log,"%s MinimizerId=%s %d maxIter=%d tgtEps=%.3e\n",
                methodName.c_str(), minId.c_str(), mutualInducedIterativeMethod, mutualInducedMaxIterations, mutualInducedTargetEpsilon );
        (void) fflush( amoebaGpu->log );
    }
#endif

    // allocate memory

    amoebaGpu->mutualInducedIterativeMethod = mutualInducedIterativeMethod;
    gpuMutualInducedFieldAllocate( amoebaGpu );
    amoebaGpu->mutualInducedMaxIterations = mutualInducedMaxIterations;
    amoebaGpu->mutualInducedTargetEpsilon = mutualInducedTargetEpsilon;

    unsigned int dipoleIndex                                                  = 0;
    unsigned int quadrupoleIndex                                              = 0;
    unsigned int maxPrint                                                     = 5;
    
    int* maxIndices = (int*) malloc( charges.size()*sizeof(int) );
    for( unsigned int ii = 0; ii < charges.size(); ii++ ){
        maxIndices[ii]                                                        = ii;
        amoebaGpu->psMultipoleParticlesIdsAndAxisType->_pSysStream[0][ii].z   = ii;
    }

    amoebaGpu->amoebaSim.electric    = electricConstant;
    if( amoebaGpu->amoebaSim.dielec < 1.0e-05 ){
        amoebaGpu->amoebaSim.dielec      = 1.0f;
    }

    for( int ii = 0; ii < static_cast<int>(charges.size()); ii++ ){

        // axis type & multipole particles ids
 
        amoebaGpu->psMultipoleParticlesIdsAndAxisType->_pSysStream[0][ii].x       = multipoleParticleId1[ii];
        amoebaGpu->psMultipoleParticlesIdsAndAxisType->_pSysStream[0][ii].y       = multipoleParticleId2[ii];
        amoebaGpu->psMultipoleParticlesIdsAndAxisType->_pSysStream[0][ii].w       = axisType[ii];
        int axisParticleIndex                                                     = multipoleParticleId1[ii];
        if( maxIndices[axisParticleIndex] < ii ){
            maxIndices[axisParticleIndex] = ii;
        }
        if( amoebaGpu->psMultipoleParticlesIdsAndAxisType->_pSysStream[0][axisParticleIndex].z > ii ){
            amoebaGpu->psMultipoleParticlesIdsAndAxisType->_pSysStream[0][axisParticleIndex].z = ii;
        }

        axisParticleIndex                                                         = multipoleParticleId2[ii];
        if( maxIndices[axisParticleIndex] < ii ){
            maxIndices[axisParticleIndex] = ii;
        }
        if( amoebaGpu->psMultipoleParticlesIdsAndAxisType->_pSysStream[0][axisParticleIndex].z > ii ){
            amoebaGpu->psMultipoleParticlesIdsAndAxisType->_pSysStream[0][axisParticleIndex].z = ii;
        }

        if( 0 && amoebaGpu->log )
            fprintf( amoebaGpu->log, "Z1 %4d [%4d %4d] %4d %4d %4d %4d   %d %d\n", ii,
                     multipoleParticleId1[ii], multipoleParticleId2[ii],
                     amoebaGpu->psMultipoleParticlesIdsAndAxisType->_pSysStream[0][multipoleParticleId1[ii]].z,
                     maxIndices[multipoleParticleId1[ii]],
                     amoebaGpu->psMultipoleParticlesIdsAndAxisType->_pSysStream[0][multipoleParticleId2[ii]].z,
                     maxIndices[multipoleParticleId2[ii]],
                     amoebaGpu->psMultipoleParticlesIdsAndAxisType->_pSysStream[0][0].z, maxIndices[0] );

        // charges

        amoebaGpu->gpuContext->psPosq4->_pSysStream[0][ii].w                  = charges[ii];
 
        // molecule dipole
 
        amoebaGpu->psMolecularDipole->_pSysStream[0][dipoleIndex]             = dipoles[dipoleIndex]; 
        amoebaGpu->psMolecularDipole->_pSysStream[0][dipoleIndex+1]           = dipoles[dipoleIndex+1]; 
        amoebaGpu->psMolecularDipole->_pSysStream[0][dipoleIndex+2]           = dipoles[dipoleIndex+2]; 
        dipoleIndex                                                          += 3;
 
        // molecule quadrupole
 
        amoebaGpu->psMolecularQuadrupole->_pSysStream[0][quadrupoleIndex]     = quadrupoles[quadrupoleIndex]; 
        amoebaGpu->psMolecularQuadrupole->_pSysStream[0][quadrupoleIndex+1]   = quadrupoles[quadrupoleIndex+1]; 
        amoebaGpu->psMolecularQuadrupole->_pSysStream[0][quadrupoleIndex+2]   = quadrupoles[quadrupoleIndex+2]; 
        amoebaGpu->psMolecularQuadrupole->_pSysStream[0][quadrupoleIndex+3]   = quadrupoles[quadrupoleIndex+3]; 
        amoebaGpu->psMolecularQuadrupole->_pSysStream[0][quadrupoleIndex+4]   = quadrupoles[quadrupoleIndex+4]; 
        amoebaGpu->psMolecularQuadrupole->_pSysStream[0][quadrupoleIndex+5]   = quadrupoles[quadrupoleIndex+5]; 
        amoebaGpu->psMolecularQuadrupole->_pSysStream[0][quadrupoleIndex+6]   = quadrupoles[quadrupoleIndex+6]; 
        amoebaGpu->psMolecularQuadrupole->_pSysStream[0][quadrupoleIndex+7]   = quadrupoles[quadrupoleIndex+7]; 
        amoebaGpu->psMolecularQuadrupole->_pSysStream[0][quadrupoleIndex+8]   = quadrupoles[quadrupoleIndex+8]; 
        quadrupoleIndex                                                      += 9;

        // damping factor 

        amoebaGpu->psDampingFactorAndThole->_pSysStream[0][ii].x              = dampingFactors[ii];
        amoebaGpu->psDampingFactorAndThole->_pSysStream[0][ii].y              = tholes[ii];

        // polarizability

        amoebaGpu->psPolarizability->_pSysStream[0][3*ii]                     = polarity[ii];
        amoebaGpu->psPolarizability->_pSysStream[0][3*ii+1]                   = polarity[ii];
        amoebaGpu->psPolarizability->_pSysStream[0][3*ii+2]                   = polarity[ii];
 
        // psCovalentDegree & psPolarizationDegree are arrays of size maxCovalentDegreeSz*paddedNumberOfAtoms

        const int particlesOffset                                             = ii*amoebaGpu->maxCovalentDegreeSz;
        const int minCovalentIndex                                            = minCovalentIndices[ii];
        amoebaGpu->psCovalentDegree->_pSysStream[0][particlesOffset]          = minCovalentIndex;

        // covalent info

        const std::vector< std::vector<int> >& covalentInfo = multipoleParticleCovalentInfo[ii];
        for( unsigned int jj = 0; jj < 4; jj++ ){
            const std::vector<int> covalentList = covalentInfo[jj];
            for( unsigned int kk = 0; kk < covalentList.size(); kk++ ){
                int covalentIndex = covalentList[kk] - minCovalentIndex + 1; 
                if( covalentIndex > amoebaGpu->maxCovalentDegreeSz || covalentIndex < 0 ){ 
                     message << "particles=" << ii << " covalent degree covalentIndex=" << covalentList[kk] << " (offset value=" << covalentIndex << " min for offset=" << minCovalentIndex <<
                               ") is out of range -- maxCovalentDegreeSz needs to be increased." << std::endl;
                    errorCount++;
                } else {
                    amoebaGpu->psCovalentDegree->_pSysStream[0][particlesOffset+covalentIndex] = covalentDegree[jj] + 1;
                }
            }
        }

        // polarization covalent info

        const int minCovalentPolarizationIndex                                = minCovalentPolarizationIndices[ii];
        amoebaGpu->psPolarizationDegree->_pSysStream[0][particlesOffset]      = minCovalentPolarizationIndex;

        for( unsigned int jj = 4; jj < covalentInfo.size(); jj++ ){
            const std::vector<int> covalentList = covalentInfo[jj];
            for( unsigned int kk = 0; kk < covalentList.size(); kk++ ){
                int covalentIndex = covalentList[kk] - minCovalentPolarizationIndex + 1; 
                if( covalentIndex > amoebaGpu->maxCovalentDegreeSz || covalentIndex < 0 ){ 
                     message << "particles=" << ii << " covalent polarization degree covalentIndex=" << covalentList[kk] << " (offset value=" << covalentIndex << " min for offset=" << minCovalentPolarizationIndex <<
                               ") is out of range -- maxCovalentDegreeSz needs to be increased." << std::endl;
                    errorCount++;
                } else {
                    amoebaGpu->psPolarizationDegree->_pSysStream[0][particlesOffset+covalentIndex] = covalentDegree[jj] + 1;
                }
            }
        }

#ifdef AMOEBA_DEBUG
        //if( amoebaGpu->log && ( ( ( ii < maxPrint ) || (ii >= (charges.size() - maxPrint) ) ) || ( ii == targetAtoms[0] || ii == targetAtoms[1]) ) ){
        if( amoebaGpu->log && ( ( ( ii < maxPrint ) || (ii >= (charges.size() - maxPrint) ) ) ) ){

            // axis particles

            (void) fprintf( amoebaGpu->log,"%u axis particles [%6d %6d %6d diff=%d %d] ", ii,
                            amoebaGpu->psMultipoleParticlesIdsAndAxisType->_pSysStream[0][ii].x,
                            amoebaGpu->psMultipoleParticlesIdsAndAxisType->_pSysStream[0][ii].y,
                            amoebaGpu->psMultipoleParticlesIdsAndAxisType->_pSysStream[0][ii].z,
                            maxIndices[ii], maxIndices[ii] - amoebaGpu->psMultipoleParticlesIdsAndAxisType->_pSysStream[0][ii].z,
                            amoebaGpu->psMultipoleParticlesIdsAndAxisType->_pSysStream[0][ii].w );

            // dipole

            int dipoleOffset = 3*ii; 
            (void) fprintf( amoebaGpu->log,"d[%16.9e %16.9e %16.9e]\n",
                            amoebaGpu->psMolecularDipole->_pSysStream[0][dipoleOffset],
                            amoebaGpu->psMolecularDipole->_pSysStream[0][dipoleOffset+1],
                            amoebaGpu->psMolecularDipole->_pSysStream[0][dipoleOffset+2] );

            // quadrupole

            int quadrupoleOffset = 9*ii;
            (void) fprintf( amoebaGpu->log,"\nq[" );
            for( int jj = 0; jj < 9; jj++ ){
               (void) fprintf( amoebaGpu->log,"%16.9e ", 
                               amoebaGpu->psMolecularQuadrupole->_pSysStream[0][quadrupoleOffset+jj] );
               if( jj == 2 || jj == 5 ){
                  (void) fprintf( amoebaGpu->log,"\n  ");
               }
            }
            (void) fprintf( amoebaGpu->log,"]\n\n" );

            // covalent/polarization degree

            if( ii < 1 ){ 
                (void) fprintf( amoebaGpu->log,"Gamma=%.5f scaledDistCutoff=%.5f\n",
                                amoebaGpu->pGamma, amoebaGpu->scalingDistanceCutoff );
            }

            (void) fprintf( amoebaGpu->log,"%3d covalent/polarization degree: minIdx[%6d %6d] Thole=%12.5f dampingFactor=%12.5f\n", ii,
                            amoebaGpu->psCovalentDegree->_pSysStream[0][particlesOffset],  amoebaGpu->psPolarizationDegree->_pSysStream[0][particlesOffset],
                            amoebaGpu->psDampingFactorAndThole->_pSysStream[0][ii].y, amoebaGpu->psDampingFactorAndThole->_pSysStream[0][ii].x );

            // covalent

            for( int kk = 1; kk < 6; kk++ ){

                const float polarScale[5]   = { 0.0f, 0.0f, 0.0f, 1.0f, 1.0f };

                // print entries w/ degree=kk

                int count = 0;
                for( int jj = 1; jj < amoebaGpu->maxCovalentDegreeSz; jj++ ){
                    if( amoebaGpu->psCovalentDegree->_pSysStream[0][particlesOffset+jj] == kk ){ 
                        if( count == 0 ){
                            (void) fprintf( amoebaGpu->log,"%d [", kk );
                        }
                        float pScale        = polarScale[kk-1];
                        int particle2Index  = amoebaGpu->psCovalentDegree->_pSysStream[0][particlesOffset] + jj - 1;
                        if( kk == 4 && particle2Index >= amoebaGpu->psPolarizationDegree->_pSysStream[0][particlesOffset] ){
                            int particle2Offset = particle2Index - amoebaGpu->psPolarizationDegree->_pSysStream[0][particlesOffset] + 1;
                            if( particle2Offset < amoebaGpu->maxCovalentDegreeSz && amoebaGpu->psPolarizationDegree->_pSysStream[0][particlesOffset+particle2Offset] == 1 ){
                                pScale *= 0.5;
                            }
                        }
                        (void) fprintf( amoebaGpu->log,"%5d %5.1f   ",
                                        amoebaGpu->psCovalentDegree->_pSysStream[0][particlesOffset] + jj - 1, pScale );
                        count++;
                   }
                }
                if( count ){
                    (void) fprintf( amoebaGpu->log,"] Sz=%5d\n", count );
                }
            }

            // polarization

            for( int kk = 1; kk < 5; kk++ ){

                // print entries w/ degree=kk

                int count = 0;
                for( int jj = 1; jj < amoebaGpu->maxCovalentDegreeSz; jj++ ){
                    if( amoebaGpu->psPolarizationDegree->_pSysStream[0][particlesOffset+jj] == kk ){ 
                        if( count == 0 ){
                            (void) fprintf( amoebaGpu->log,"%d [", kk );
                        }
                        (void) fprintf( amoebaGpu->log,"%5d ", amoebaGpu->psPolarizationDegree->_pSysStream[0][particlesOffset] + jj - 1 );
                        count++;
                    }
                }
                if( count ){
                    const float directScale[5]  = { 0.0f, 1.0f, 1.0f, 1.0f, 1.0f };
                    float dScale = directScale[kk-1];
                    (void) fprintf( amoebaGpu->log,"] Sz=%5d d=%.2f\n", count, dScale );
                }
            }
            (void) fprintf( amoebaGpu->log,"\n" );
            (void) fflush( amoebaGpu->log );
        }
        if( amoebaGpu->psDampingFactorAndThole->_pSysStream[0][ii].x != amoebaGpu->psDampingFactorAndThole->_pSysStream[0][ii].x ||
            amoebaGpu->psDampingFactorAndThole->_pSysStream[0][ii].x == std::numeric_limits<double>::infinity() || 
            amoebaGpu->psDampingFactorAndThole->_pSysStream[0][ii].x == -std::numeric_limits<double>::infinity()){
            (void) fprintf( amoebaGpu->log,"Nan detected at index=%d in psDampingFactor\n", ii );
        }
#endif

#ifdef AMOEBA_DEBUG
        if( amoebaGpu->log ){

            // covalent

            const float polarScale[5]   = { 0.0f, 0.0f, 0.0f, 1.0f, 1.0f };
            const float directScale[5]  = { 0.0f, 1.0f, 1.0f, 1.0f, 1.0f };
            const float mpoleScale[5]   = { 0.0f, 0.0f, 0.0f, 0.4f, 0.8f };

            // print entries w/ degree=kk

            for( int jj = 1; jj < amoebaGpu->maxCovalentDegreeSz; jj++ ){
                if( amoebaGpu->psCovalentDegree->_pSysStream[0][particlesOffset+jj] ){ 
                    int index           = amoebaGpu->psCovalentDegree->_pSysStream[0][particlesOffset+jj];
                    float pScale        = polarScale[index-1];
                    float mScale        = mpoleScale[index-1];
                    int particle2Index  = amoebaGpu->psCovalentDegree->_pSysStream[0][particlesOffset] + jj - 1;
                    if( index == 4 && particle2Index >= amoebaGpu->psPolarizationDegree->_pSysStream[0][particlesOffset] ){
                        int particle2Offset = particle2Index - amoebaGpu->psPolarizationDegree->_pSysStream[0][particlesOffset] + 1;
                        if( particle2Offset < amoebaGpu->maxCovalentDegreeSz && amoebaGpu->psPolarizationDegree->_pSysStream[0][particlesOffset+particle2Offset] == 1 ){
                            pScale *= 0.5;
                        }
                    }
                    pScaleCheckSum[ii] += (pScale - 1.0f);
                    int covIndex        = amoebaGpu->psCovalentDegree->_pSysStream[0][particlesOffset];
                    if( pScale != 1.0f ){
                        MapIntFloat* pMap = amoebaGpu->pMapArray[ii];
                        (*pMap)[covIndex+jj-1] = pScale;
                    }
                }
            }

            // polarization

            for( int jj = 1; jj < amoebaGpu->maxCovalentDegreeSz; jj++ ){
                if( amoebaGpu->psPolarizationDegree->_pSysStream[0][particlesOffset+jj] ){ 
                    int index    = amoebaGpu->psPolarizationDegree->_pSysStream[0][particlesOffset+jj];
                    dScaleCheckSum[ii] += (directScale[index-1] - 1.0f);
                    int covIndex        = amoebaGpu->psPolarizationDegree->_pSysStream[0][particlesOffset];
                    if( directScale[index-1] != 1.0f ){
                        MapIntFloat* dMap      = amoebaGpu->dMapArray[ii];
                        (*dMap)[covIndex+jj-1] = directScale[index-1];
                    }
                }
            }
        }
#endif

    }

#if 0
    if( amoebaGpu->log ){
        FILE* filePtr = fopen( "oldScale.txt", "w" );
        for( unsigned int kk = 0; kk < charges.size(); kk++ ){
            (void) fprintf( filePtr, "%6u %14.6e %14.6e\n", kk, pScaleCheckSum[kk], dScaleCheckSum[kk] );
        }
        (void) fclose( filePtr );
        free( pScaleCheckSum );
        free( dScaleCheckSum );
        free( mScaleCheckSum );
        filePtr = fopen( "oldScaleMap.txt", "w" );
        for( unsigned int kk = 0; kk < charges.size(); kk++ ){
            MapIntFloat* pMap = amoebaGpu->pMapArray[kk];
            
            float sum = 0.0f;
            for( MapIntFloatCI ii = pMap->begin(); ii != pMap->end(); ii++ ){
                sum += (*ii).second;
            }
            (void) fprintf( filePtr, "%6u sz=%u sum=%.3f total=%.3f %.3f\n", kk, pMap->size(),
                            sum, (float) charges.size() - sum, 1248.0f - sum );
            for( MapIntFloatCI ii = pMap->begin(); ii != pMap->end(); ii++ ){
                (void) fprintf( filePtr, "    %6d %14.6e\n", (*ii).first, (*ii).second );
            }
        }
        (void) fclose( filePtr );
    }
#endif

    amoebaGpu->maxMapTorqueDifference = 0;
    for( unsigned int ii = 0; ii < charges.size(); ii++ ){

        // axis type & multipole particles ids
        int diff = maxIndices[ii] - amoebaGpu->psMultipoleParticlesIdsAndAxisType->_pSysStream[0][ii].z;
        if( diff > amoebaGpu->maxMapTorqueDifference ){
            amoebaGpu->maxMapTorqueDifference = diff;
        }
    }

    // set maxMapTorqueDifferencePow2 to smallest power of 2 greater than maxMapTorqueDifference

    float logDiff                                = logf( static_cast<float>(amoebaGpu->maxMapTorqueDifference) )/logf(2.0f);    
    int logDiffI                                 = static_cast<int>(logDiff) + 1;
    amoebaGpu->maxMapTorqueDifferencePow2        = static_cast<int>(powf( 2.0f, static_cast<float>(logDiffI) ));

    amoebaGpu->torqueMapForce                    = new CUDAStream<float>(amoebaGpu->paddedNumberOfAtoms*3*amoebaGpu->maxMapTorqueDifference, 1, "torqueMapForce");
    memset( amoebaGpu->torqueMapForce->_pSysStream[0], 0, amoebaGpu->paddedNumberOfAtoms*3*amoebaGpu->maxMapTorqueDifference*sizeof( float ) );
    amoebaGpu->torqueMapForce->Upload();

    free( maxIndices );

    amoebaGpuBuildOutputBuffers( amoebaGpu );
    amoebaGpuBuildThreadBlockWorkList( amoebaGpu );
    amoebaGpuBuildScalingList( amoebaGpu );

    // upload
 
    amoebaGpu->amoebaSim.scalingDistanceCutoff = amoebaGpu->scalingDistanceCutoff;
    amoebaGpu->amoebaSim.numberOfAtoms         = amoebaGpu->gpuContext->natoms;
    amoebaGpu->amoebaSim.paddedNumberOfAtoms   = amoebaGpu->paddedNumberOfAtoms;

    amoebaGpu->psMultipoleParticlesIdsAndAxisType->Upload();
    amoebaGpu->psMolecularDipole->Upload();
    amoebaGpu->psMolecularQuadrupole->Upload();
    amoebaGpu->psCovalentDegree->Upload();
    amoebaGpu->psPolarizationDegree->Upload();
    amoebaGpu->psDampingFactorAndThole->Upload();
    amoebaGpu->psPolarizability->Upload();
    amoebaGpu->gpuContext->psPosq4->Upload();

    gpuElectrostaticAllocate( amoebaGpu );
}

extern "C"
void gpuSetAmoebaObcParameters( amoebaGpuContext amoebaGpu, float innerDielectric, float solventDielectric, float dielectricOffset,
                                const std::vector<float>& radius, const std::vector<float>& scale, const std::vector<float>& charge,
                                int includeCavityTerm, float probeRadius, float surfaceAreaFactor )
{

    gpuContext gpu                         = amoebaGpu->gpuContext;
    gpu->sim.dielectricOffset              = dielectricOffset;
    amoebaGpu->includeObcCavityTerm        = includeCavityTerm;
    gpu->sim.probeRadius                   = probeRadius;
    gpu->sim.surfaceAreaFactor             = surfaceAreaFactor;
    unsigned int particles                 = radius.size();

    for (unsigned int i = 0; i < particles; i++) 
    {    
            (*gpu->psObcData)[i].x = radius[i] - dielectricOffset;
            (*gpu->psObcData)[i].y = scale[i] * (*gpu->psObcData)[i].x;
            (*gpu->psPosq4)[i].w   = charge[i];
    }    

    // Dummy out extra particles data
    for (unsigned int i = particles; i < amoebaGpu->paddedNumberOfAtoms; i++) 
    {    
        (*gpu->psBornRadii)[i]     = 0.2f;
        (*gpu->psObcData)[i].x     = 0.01f;
        (*gpu->psObcData)[i].y     = 0.01f;
    }    

    gpu->psBornRadii->Upload();
    gpu->psObcData->Upload();
    gpu->psPosq4->Upload();

    amoebaGpu->amoebaSim.gkc         = 2.455f;
    amoebaGpu->amoebaSim.dwater      = solventDielectric;
    amoebaGpu->amoebaSim.dielec      = innerDielectric;

//    amoebaGpu->amoebaSim.fc          = amoebaGpu->amoebaSim.electric*1.0f*(1.0f-solventDielectric)/(0.0f+1.0f*solventDielectric);;
//    amoebaGpu->amoebaSim.fd          = amoebaGpu->amoebaSim.electric*2.0f*(1.0f-solventDielectric)/(1.0f+2.0f*solventDielectric);;
//    amoebaGpu->amoebaSim.fq          = amoebaGpu->amoebaSim.electric*3.0f*(1.0f-solventDielectric)/(2.0f+3.0f*solventDielectric);;

    amoebaGpu->amoebaSim.fc          = 1.0f*(1.0f-solventDielectric)/(0.0f+1.0f*solventDielectric);;
    amoebaGpu->amoebaSim.fd          = 2.0f*(1.0f-solventDielectric)/(1.0f+2.0f*solventDielectric);;
    amoebaGpu->amoebaSim.fq          = 3.0f*(1.0f-solventDielectric)/(2.0f+3.0f*solventDielectric);;

    gpu->sim.preFactor               = -amoebaGpu->amoebaSim.electric*((1.0f/innerDielectric)-(1.0f/solventDielectric));

    if( amoebaGpu->log ){
        (void) fprintf( amoebaGpu->log,"gpuSetAmoebaObcParameters: cavity=%d dielectricOffset=%15.7e probeRadius=%15.7e surfaceAreaFactor=%15.7e\n", 
                        includeCavityTerm, dielectricOffset, probeRadius, surfaceAreaFactor );
        (void) fprintf( amoebaGpu->log,"                           gkc=%12.3f solventDielectric=%15.7e innerDielectric=%15.7e sim.preFactor=%15.7e\n", 
                        amoebaGpu->amoebaSim.gkc, amoebaGpu->amoebaSim.dwater, amoebaGpu->amoebaSim.dielec, gpu->sim.preFactor );
        (void) fprintf( amoebaGpu->log,"                           fc=%15.7e fd=%15.7e fq=%15.7e\n",
                        amoebaGpu->amoebaSim.fc, amoebaGpu->amoebaSim.fq, amoebaGpu->amoebaSim.fq );
        (void) fprintf( amoebaGpu->log,"\nRadius (r-off) scl*(r-off) scl\n" );
        for (unsigned int i = 0; i < amoebaGpu->paddedNumberOfAtoms && i < 10; i++) 
        {
            (void) fprintf( amoebaGpu->log,"%6d %15.7e %15.7e %15.7e %15.7e\n", i,
                            radius[i] , (*gpu->psObcData)[i].x, (*gpu->psObcData)[i].y, scale[i] );
        }
    }

    gpuRotationToLabFrameAllocate( amoebaGpu );
    gpuFixedEFieldAllocate( amoebaGpu );
    gpuElectrostaticAllocate( amoebaGpu );
    gpuKirkwoodAllocate( amoebaGpu );

}

static int encodeCell( unsigned int x, unsigned int y ){
    return ( (x << 17) | (y << 2) );
}

static int encodeCellExclusion( unsigned int cellCode ){
    return cellCode |= 1;
}

static int decodeCell( int cellCode, unsigned int* x, unsigned int* y, unsigned int* exclusion ){
    *x         =  cellCode >> 17;
    *y         = (cellCode >> 2 ) & 0x7FFF;
    *exclusion = (cellCode & 1) ? 1 : 0;
	return 0;
}

extern "C"
void gpuSetAmoebaVdwParameters( amoebaGpuContext amoebaGpu,
                                const std::vector<int>& indexIVs, 
                                const std::vector<int>& indexClasses, 
                                const std::vector<float>& sigmas,
                                const std::vector<float>& epsilons,
                                const std::vector<float>& sigma4s,
                                const std::vector<float>& epsilon4s,
                                const std::vector<float>& reductions,
                                const std::string& vdwSigmaCombiningRule,
                                const std::string& vdwEpsilonCombiningRule,
                                const std::vector< std::vector< std::vector<float> > >& sigEpsTable,
                                const std::vector< std::vector<int> >& allExclusions )
{
   // ---------------------------------------------------------------------------------------

    static const char* methodName = "gpuSetAmoebaVdwParameters";

   // ---------------------------------------------------------------------------------------

    gpuContext gpu                         = amoebaGpu->gpuContext;
    amoebaGpu->paddedNumberOfAtoms         = amoebaGpu->gpuContext->sim.paddedNumberOfAtoms;
    unsigned int particles                 = sigmas.size();
    amoebaGpu->useVdwTable                 = 0;
    
    // set sigma combining rule flag

    if( vdwSigmaCombiningRule.compare( "ARITHMETIC" ) == 0 ){
        amoebaGpu->vdwSigmaCombiningRule = 1;
    } else if( vdwSigmaCombiningRule.compare( "GEOMETRIC" ) == 0 ){
        amoebaGpu->vdwSigmaCombiningRule = 2;
    } else if( vdwSigmaCombiningRule.compare( "CUBIC-MEAN" ) == 0 ){
        amoebaGpu->vdwSigmaCombiningRule = 3;
    } else {
        amoebaGpu->vdwSigmaCombiningRule = 1;
        if( amoebaGpu->log ){
            (void) fprintf( amoebaGpu->log, "%s sigma combining rule=<%s> not recognized; using ARITHMETIC\n", methodName, vdwSigmaCombiningRule.c_str() );
        } else {
            (void) fprintf( stderr, "%s sigma combining rule=<%s> not recognized; using ARITHMETIC\n", methodName, vdwSigmaCombiningRule.c_str() );
        }
    }

    // set epsilon combining rule flag

    if( vdwEpsilonCombiningRule.compare( "ARITHMETIC" ) == 0 ){
        amoebaGpu->vdwEpsilonCombiningRule = 1;
    } else if( vdwEpsilonCombiningRule.compare( "GEOMETRIC" ) == 0 ){
        amoebaGpu->vdwEpsilonCombiningRule = 2;
    } else if( vdwEpsilonCombiningRule.compare( "HARMONIC" ) == 0 ){
        amoebaGpu->vdwEpsilonCombiningRule = 3;
    } else if( vdwEpsilonCombiningRule.compare( "HHG" ) == 0 ){
        amoebaGpu->vdwEpsilonCombiningRule = 4;
    } else {
        amoebaGpu->vdwEpsilonCombiningRule = 1;
        if( amoebaGpu->log ){
            (void) fprintf( amoebaGpu->log, "%s epsilon combining rule=<%s> not recognized; using ARITHMETIC\n", methodName, vdwEpsilonCombiningRule.c_str() );
        } else {
            (void) fprintf( stderr, "%s epsilon combining rule=<%s> not recognized; using ARITHMETIC\n", methodName, vdwEpsilonCombiningRule.c_str() );
        }
    }

    if( amoebaGpu->useVdwTable ){
        amoebaGpu->vdwTableSize            = sigEpsTable.size();
        amoebaGpu->psVdwTable              = new CUDAStream<float2>( amoebaGpu->vdwTableSize*amoebaGpu->vdwTableSize,  1, "VdwTable");
        for (unsigned int ii = 0; ii < amoebaGpu->vdwTableSize; ii++) 
        {    
            for (unsigned int jj = 0; jj < amoebaGpu->vdwTableSize; jj++) 
            {    
                amoebaGpu->psVdwTable->_pSysStream[0][ii].x  = sigEpsTable[ii][jj][0];
                amoebaGpu->psVdwTable->_pSysStream[0][ii].y  = sigEpsTable[ii][jj][1];
            }    
        }    
        amoebaGpu->psVdwTable->Upload();
    } else {
        if( particles < 1 ){
            (void) fprintf( stderr, "%s number of particles\n", methodName );
            return;
        } 
    
        amoebaGpu->psVdwSigmaEpsilon           = new CUDAStream<float2>(amoebaGpu->paddedNumberOfAtoms,   1, "VdwSigmaEpsilon");
        for (unsigned int ii = 0; ii < particles; ii++) 
        {    
            amoebaGpu->psVdwSigmaEpsilon->_pSysStream[0][ii].x    = sigmas[ii];
            amoebaGpu->psVdwSigmaEpsilon->_pSysStream[0][ii].y    = epsilons[ii];
        }    
    
        // Dummy out extra particles data

        for (unsigned int ii = particles; ii < amoebaGpu->paddedNumberOfAtoms; ii++) 
        {    
            amoebaGpu->psVdwSigmaEpsilon->_pSysStream[0][ii].x     = 1.0f;
            amoebaGpu->psVdwSigmaEpsilon->_pSysStream[0][ii].y     = 0.0f;
        }    
        amoebaGpu->psVdwSigmaEpsilon->Upload();
    }

    std::vector< std::vector<unsigned int> > ivMapping;
    std::vector< unsigned int > ivNonMapping;
    std::vector<float> reductionFactors;
    ivMapping.resize( particles );
    ivNonMapping.resize( particles );
    reductionFactors.resize( particles );
    for (unsigned int ii = 0; ii < particles; ii++) 
    {    
        ivNonMapping[ii] = 1;
    }
    for (unsigned int ii = 0; ii < particles; ii++) 
    {    
        if( ii != indexIVs[ii] )
        {
            ivMapping[indexIVs[ii]].push_back( ii ); 

            reductionFactors[indexIVs[ii]] = reductions[ii];

            ivNonMapping[ii]               = 0;
            ivNonMapping[indexIVs[ii]]     = 0;
        }
    }    

   unsigned int numberOfNonReductions = 0;
   unsigned int numberOfReductions    = 0;
    for (unsigned int ii = 0; ii < particles; ii++) 
    {    
        if( ivNonMapping[ii] )
        {
            numberOfNonReductions++;
        }
        if( ivMapping[ii].size() > 0 )
        {
            numberOfReductions++;
#ifdef AMOEBA_DEBUG
            (void) fprintf( amoebaGpu->log, "Atom %u has %u reductions: [", ii, ivMapping[ii].size() );
            for (unsigned int jj = 0; jj < ivMapping[ii].size(); jj++) 
            {
                (void) fprintf( amoebaGpu->log, "%u ", ivMapping[ii][jj] );
            }
            (void) fprintf( amoebaGpu->log, "]  %12.4f\n", reductionFactors[ii] );
#endif
        }
    }    
    
    amoebaGpu->amoebaSim.amoebaVdwNonReductions        = numberOfNonReductions;
    amoebaGpu->amoebaSim.amoebaVdwReductions           = numberOfReductions;

    CUDAStream<int>* psVdwNonReductionID               = new CUDAStream<int>(numberOfNonReductions, 1, "AmoebaVdwNonReductionID");
    amoebaGpu->psAmoebaVdwNonReductionID               = psVdwNonReductionID;
    amoebaGpu->amoebaSim.pAmoebaVdwNonReductionID      = psVdwNonReductionID->_pDevStream[0];

    CUDAStream<int4>* psVdwReductionID                 = new CUDAStream<int4>(numberOfReductions, 1, "AmoebaVdwReductionID");
    amoebaGpu->psAmoebaVdwReductionID                  = psVdwReductionID;
    amoebaGpu->amoebaSim.pAmoebaVdwReductionID         = psVdwReductionID->_pDevStream[0];

    CUDAStream<float>* psAmoebaVdwReduction            = new CUDAStream<float>(numberOfReductions, 1, "AmoebaVdwReduction");
    amoebaGpu->psAmoebaVdwReduction                    = psAmoebaVdwReduction;
    amoebaGpu->amoebaSim.pAmoebaVdwReduction           = psAmoebaVdwReduction->_pDevStream[0];        
    amoebaGpu->psAmoebaVdwCoordinates                  = new CUDAStream<float4>( amoebaGpu->paddedNumberOfAtoms,  1, "VdwCoordinates");

   unsigned int count    = 0;
   unsigned int nonCount = 0;
    for (unsigned int ii = 0; ii < particles; ii++) 
    {    
        if( ivNonMapping[ii] )
        {
            psVdwNonReductionID->_pSysStream[0][nonCount++] = ii;
        }

        if( ivMapping[ii].size() > 0 )
        {
            psAmoebaVdwReduction->_pSysStream[0][count] = reductionFactors[ii];
            psVdwReductionID->_pSysStream[0][count].x   = ii;
            psVdwReductionID->_pSysStream[0][count].y   = ivMapping[ii][0];
            if( ivMapping[ii].size() > 1 ){
                psVdwReductionID->_pSysStream[0][count].z   = ivMapping[ii][1];
            } else {
                psVdwReductionID->_pSysStream[0][count].z   = ii;
            }
            if( ivMapping[ii].size() > 2 ){
                psVdwReductionID->_pSysStream[0][count].w   = ivMapping[ii][2];
            } else {
                psVdwReductionID->_pSysStream[0][count].w   = ii;
            }
            if( ivMapping[ii].size() > 3 ){
                (void) fprintf( stderr, "Atom %u has %u reductions -- invalid -- aborting", ii, static_cast<unsigned int>(ivMapping[ii].size()) );
                exit(1);
            }
            count++;
        }
    }    
    psVdwNonReductionID->Upload();
    psVdwReductionID->Upload();
    psAmoebaVdwReduction->Upload();

    amoebaGpuBuildVdwExclusionList( amoebaGpu, allExclusions );

#ifdef AMOEBA_DEBUG
    if( amoebaGpu->log ){
        unsigned int maxPrint = 32;
        (void) fprintf( amoebaGpu->log, "%s useVdwTable=%d size=%d\n", methodName, amoebaGpu->useVdwTable, 
                        (amoebaGpu->useVdwTable ? amoebaGpu->vdwTableSize : 0) ); 
        (void) fprintf( amoebaGpu->log, "%s sigma/epsilon combining rules=%d %d\n", methodName, 
                        amoebaGpu->vdwSigmaCombiningRule, amoebaGpu->vdwEpsilonCombiningRule);
        for (unsigned int ii = 0; ii < gpu->natoms; ii++) 
        {    
            (void) fprintf( amoebaGpu->log, "%5u %15.7e %15.7e\n", ii, sigmas[ii], epsilons[ii] );
            if( ii == maxPrint && ii < (amoebaGpu->paddedNumberOfAtoms - maxPrint) )
            {
                ii = (amoebaGpu->paddedNumberOfAtoms - maxPrint);
            }
        } 
        (void) fflush( amoebaGpu->log );
    }    
#endif

}

extern "C"
void amoebaGpuBuildVdwExclusionList( amoebaGpuContext amoebaGpu,  const std::vector< std::vector<int> >& exclusions )
{
    // ---------------------------------------------------------------------------------------

    static const std::string methodName = "amoebaGpuBuildVdwExclusionList";
    static const int debugOn            = 0;

    // ---------------------------------------------------------------------------------------

    amoebaGpuBuildThreadBlockWorkList( amoebaGpu );

    const unsigned int paddedAtoms     = amoebaGpu->paddedNumberOfAtoms;
    const unsigned int actualAtoms     = amoebaGpu->gpuContext->natoms;
    const unsigned int grid            = amoebaGpu->gpuContext->grid;
    const unsigned int dim             = paddedAtoms/grid;
    const unsigned int cells           = dim * (dim + 1) / 2;
    unsigned int* pWorkList            = amoebaGpu->psVdwWorkUnit->_pSysData;

    // minCellIndex & maxCellIndex track min/max atom index for each cell

    std::vector<int> minCellIndex( dim + 1 );
    std::vector<int> maxCellIndex( dim + 1 );
    for(unsigned int ii = 0; ii <= dim; ii++)
    {
        minCellIndex[ii] = paddedAtoms + 1;
        maxCellIndex[ii] = 0;
    }

    for(unsigned int atom1 = 0; atom1 < actualAtoms; atom1++)
    {
        int x                  = atom1/grid;
        for ( unsigned int jj =  0; jj < exclusions[atom1].size(); jj++ )
        {
            if( exclusions[atom1][jj] > maxCellIndex[x] )
            { 
                maxCellIndex[x] =  exclusions[atom1][jj];
            }
            if( exclusions[atom1][jj] < minCellIndex[x] )
            { 
                minCellIndex[x] =  exclusions[atom1][jj];
            }
        }
    }

    // diagnostics

    if( debugOn && amoebaGpu->log ){
        (void) fprintf( amoebaGpu->log, "%s min/max cell indices:\n", methodName.c_str() );
        for ( unsigned int ii = 0; ii < dim; ii++)
        {
            (void) fprintf( amoebaGpu->log, "%6d [%6d %6d]\n", ii, minCellIndex[ii], maxCellIndex[ii] );
        }
        (void) fflush( amoebaGpu->log );
    }

    // Build a list of indexes for the work units with exclusions

    CUDAStream<int>* psVdwExclusionIndicesIndex     = new CUDAStream<int>(cells, 1u, "VdwExclusionIndicesIndex");
    amoebaGpu->psVdwExclusionIndicesIndex           = psVdwExclusionIndicesIndex;
    amoebaGpu->amoebaSim.pVdwExclusionIndicesIndex  = psVdwExclusionIndicesIndex->_pDevStream[0];

    //memset( amoebaGpu->psVdwExclusionIndicesIndex->_pSysStream[0], 0, sizeof(cells)*sizeof(int) );
    for (unsigned int ii = 0; ii < cells; ii++)
    {
        amoebaGpu->psVdwExclusionIndicesIndex->_pSysStream[0][ii] = -1;
    }
    int numWithExclusionIndices                     = 0;
    int gridOffset                                  = grid - 1;
    int lastBlock                                   = (static_cast<int>(paddedAtoms) > amoebaGpu->gpuContext->natoms) ? (amoebaGpu->gpuContext->natoms)/grid : -1;
    for (unsigned int ii = 0; ii < cells; ii++)
    {
        unsigned int x, y, exclusion;
        decodeCell( pWorkList[ii], &x, &y, &exclusion );
        int xAtomMin = x*grid;
        int xAtomMax = xAtomMin + gridOffset;
        if( (maxCellIndex[y] >= xAtomMin && minCellIndex[y] <= xAtomMax) || (x == lastBlock || y == lastBlock) ){
            pWorkList[ii]                                   = encodeCellExclusion( pWorkList[ii] );
            psVdwExclusionIndicesIndex->_pSysStream[0][ii]  = grid*numWithExclusionIndices;
            numWithExclusionIndices++;
//(void) fprintf( amoebaGpu->log, "%5d cellIdx[%6d %6d] atomCell[%6d %6d] exclCell[%6d %6d] num=%5d last=%5d\n",
//                ii, x, y, xAtomMin, xAtomMax, minCellIndex[y], maxCellIndex[y], numWithExclusionIndices, lastBlock );
        }
    }
    for ( unsigned int ii = 0; ii < cells; ii++)
    {
        if( amoebaGpu->psVdwExclusionIndicesIndex->_pSysStream[0][ii] == -1 )
        {
            amoebaGpu->psVdwExclusionIndicesIndex->_pSysStream[0][ii] = grid*numWithExclusionIndices;
        }
    }

    // diagnostics

    if( debugOn && amoebaGpu->log ){
        (void) fprintf( amoebaGpu->log, "%s %d cells w/ exclusions\n", methodName.c_str(), numWithExclusionIndices );
        for (unsigned int ii = 0; ii < cells; ii++)
        {
            unsigned int x, y, exclusion;
            decodeCell( pWorkList[ii], &x, &y, &exclusion );
            if( exclusion ){
                (void) fprintf( amoebaGpu->log, "%6d [%6u %6u] %8u %8u [%6u %6u] has excl: indexInToIndices=%8d\n",
                                ii, x, y, exclusion, pWorkList[ii], 32*x, 32*y,
                                psVdwExclusionIndicesIndex->_pSysStream[0][ii] );
            } else {
                (void) fprintf( amoebaGpu->log, "%6d [%6u %6u] %8u %8u [%6u %6u]\n",
                                ii, x, y, exclusion, pWorkList[ii], 32*x, 32*y );
            }
        }
        (void) fflush( amoebaGpu->log );
    }

    // Record the scaling data

    CUDAStream<int>* psVdwExclusionIndices       = new CUDAStream<int>((numWithExclusionIndices+1)*grid, 1u, "VdwExclusionIndices");
    amoebaGpu->psVdwExclusionIndices             = psVdwExclusionIndices;
    amoebaGpu->amoebaSim.pVdwExclusionIndices    = psVdwExclusionIndices->_pDevStream[0];

    memset( psVdwExclusionIndices->_pSysStream[0], 0, (numWithExclusionIndices+1)*grid*sizeof( int ) );

    // load scaling indices

    // psVdwExclusionIndicesIndex[cell] gives the index into array of masks for that cell

    // for each cell, ps_X_ScaleIndices is an int4 array of GRID(32) masks for that cell
    // ps_X_ScaleIndices[scaleIndex + threadIdI] is the mask for atomI (x + threadIdI) w/ the the y-atoms 
    // the ith bits of the mask give the scale index for atomI w/ atom (y + i)
 
    for ( unsigned int ii = 0; ii < cells; ++ii)
    {
        unsigned int x, y, exclusion;
        decodeCell( pWorkList[ii], &x, &y, &exclusion );
        if( exclusion )
        {
            int scaleIndex = psVdwExclusionIndicesIndex->_pSysStream[0][ii];  
            for (unsigned int jj = 0; jj < grid; jj++)
            {
                unsigned int atomI = grid*x + jj;
                int scaleOffset    = scaleIndex + jj;
                for (unsigned int kk = 0; kk < grid; kk++)
                {
                    unsigned int atomJ = grid*y + kk;
                    unsigned int hit   = 0;
                    if( atomI < exclusions.size() ){
                        for ( unsigned int mm =  0; mm < exclusions[atomI].size() && hit == 0; mm++ )
                        {
                            if( exclusions[atomI][mm] == atomJ )hit = 1;
                        }
                    }
                    if( (hit == 1) || (atomI >= actualAtoms) || (atomJ >= actualAtoms) ){
                        psVdwExclusionIndices->_pSysStream[0][scaleOffset] |= (1 << kk);
//if( hit )
//(void) fprintf( amoebaGpu->log, "Excluding %u %u\n", atomI, atomJ );
                    }

                }
            }
        }
    }

    // diagnostics

    if( debugOn && amoebaGpu->log ){

        (void) fprintf( amoebaGpu->log, "%s Echo exclusions\n", methodName.c_str() );
        (void) fflush( amoebaGpu->log );
        std::vector< std::map<int,int> > echoExclusions;
        echoExclusions.resize( paddedAtoms );
        for (unsigned int ii = 0; ii < cells; ii++)
        {
            unsigned int x, y, exclusion;
            decodeCell( pWorkList[ii], &x, &y, &exclusion );
            if( exclusion ){
                int scaleIndex     = psVdwExclusionIndicesIndex->_pSysStream[0][ii];  
                for (unsigned int jj = 0; jj < grid; jj++)
                {
                    unsigned int atomI        = grid*x + jj;
                    unsigned int scaleOffset  = scaleIndex + jj;
                    for (unsigned int kk = 0; kk < grid; kk++)
                    {
                        unsigned int atomJ         = grid*y + kk;
                        unsigned int mask          = 1 << kk;
                        unsigned int exclude       = psVdwExclusionIndices->_pSysStream[0][scaleOffset] & mask ? 1 : 0;
                        if( exclude ){
                             if( (atomJ < actualAtoms) && (atomI < actualAtoms) ){
                                 echoExclusions[atomI][atomJ] = 1;
                                 echoExclusions[atomJ][atomI] = 1;
                             }
                             (void) fprintf( amoebaGpu->log, "%6d cell[%6u %6u] %1u atom[%6u %6u] %6u  scaleOffset=%u %d kk=%u\n",
                                             ii, x, y, exclusion, atomI, atomJ, exclude, scaleOffset,
                                             psVdwExclusionIndices->_pSysStream[0][scaleOffset],kk );
                        }
                    }
                }
            }
        }

        (void) fprintf( amoebaGpu->log, "Check exclusions w/ echo\n" );
        int totalErrors = 0;
        for (unsigned int ii = 0; ii < actualAtoms; ii++){
            bool error = false;
            if( exclusions[ii].size() != echoExclusions[ii].size() ){
                 (void) fprintf( amoebaGpu->log, "\nAtom %6d sz %6u %6u XX\n", ii, static_cast<unsigned int>(exclusions[ii].size()), static_cast<unsigned int>(echoExclusions[ii].size()) );
                 error = true;
                 totalErrors++;
            }
            for (unsigned int jj = 0; jj < exclusions[ii].size(); jj++){
                int hit        = 0;
                int searchAtom = exclusions[ii][jj];
                for ( std::map<int,int>::iterator kk = echoExclusions[ii].begin(); kk != echoExclusions[ii].end() & !hit; kk++){
                    if( kk->first == searchAtom ){
                        hit =1;
                    }
                }
                if( !hit ){
                     error = true;
                     totalErrors++;
                     (void) fprintf( amoebaGpu->log, "Atom %6d missing %6u  XX\n", ii, searchAtom );
                }
            }
            if( !error && totalErrors ){
                (void) fprintf( amoebaGpu->log, "Atom %6d Ok\n", ii );
            }
        }
        if( totalErrors == 0 ){
            (void) fprintf( amoebaGpu->log, "No exclusion errors\n" );
        }
        (void) fflush( amoebaGpu->log );
    }

    amoebaGpu->psVdwExclusionIndices->Upload();
    amoebaGpu->psVdwExclusionIndicesIndex->Upload();
    amoebaGpu->psVdwWorkUnit->Upload();
}

extern "C"
void gpuSetAmoebaWcaDispersionParameters( amoebaGpuContext amoebaGpu,
                                          const std::vector<float>& radii,
                                          const std::vector<float>& epsilons,
                                          const float totalMaxWcaDisperionEnergy,
                                          const float epso, const float epsh, const float rmino, const float rminh,
                                          const float awater, const float shctd, const float dispoff )
{
   // ---------------------------------------------------------------------------------------

    static const char* methodName = "gpuSetAmoebaWcaDispersionParameters";

   // ---------------------------------------------------------------------------------------

    gpuContext gpu                           = amoebaGpu->gpuContext;
    amoebaGpu->paddedNumberOfAtoms           = amoebaGpu->gpuContext->sim.paddedNumberOfAtoms;
    unsigned int particles                   = radii.size();
    
    amoebaGpu->psWcaDispersionRadiusEpsilon  = new CUDAStream<float2>(amoebaGpu->paddedNumberOfAtoms,   1, "WcaDispersionRadiusEpsilon");
    for (unsigned int ii = 0; ii < particles; ii++) 
    {    
        amoebaGpu->psWcaDispersionRadiusEpsilon->_pSysStream[0][ii].x    = radii[ii];
        amoebaGpu->psWcaDispersionRadiusEpsilon->_pSysStream[0][ii].y    = epsilons[ii];
    }    
    
    // Dummy out extra particles data

    for (unsigned int ii = particles; ii < amoebaGpu->paddedNumberOfAtoms; ii++) 
    {    
        amoebaGpu->psWcaDispersionRadiusEpsilon->_pSysStream[0][ii].x     = 1.0f;
        amoebaGpu->psWcaDispersionRadiusEpsilon->_pSysStream[0][ii].y     = 0.0f;
    }    
    amoebaGpu->psWcaDispersionRadiusEpsilon->Upload();
    amoebaGpu->amoebaSim.totalMaxWcaDispersionEnergy = totalMaxWcaDisperionEnergy;
    amoebaGpu->amoebaSim.epso                        = epso;
    amoebaGpu->amoebaSim.epsh                        = epsh;
    amoebaGpu->amoebaSim.rmino                       = rmino;
    amoebaGpu->amoebaSim.rminh                       = rminh;
    amoebaGpu->amoebaSim.awater                      = awater;
    amoebaGpu->amoebaSim.shctd                       = shctd;
    amoebaGpu->amoebaSim.dispoff                     = dispoff;

#ifdef AMOEBA_DEBUG
    if( amoebaGpu->log ){
        unsigned int maxPrint = 10;
        (void) fprintf( amoebaGpu->log, "%s particles=%u total max dispersion energy=%14.5e eps[%14.5e %14.5e] rmin[%14.5e %14.5e] awtr=%14.5e shctd=%14.5e dispoff=%14.5e\n",
                        methodName, radii.size(), totalMaxWcaDisperionEnergy, epso, epsh, rmino, rminh, awater, shctd, dispoff );
        for (unsigned int ii = 0; ii < gpu->natoms; ii++) 
        {    
            (void) fprintf( amoebaGpu->log, "%5u %15.7e %15.7e\n", ii, radii[ii], epsilons[ii] );
            if( ii == maxPrint && ii < (amoebaGpu->paddedNumberOfAtoms - maxPrint) )
            {
                ii = (amoebaGpu->paddedNumberOfAtoms - maxPrint);
            }
        } 
        (void) fflush( amoebaGpu->log );
    }    
#endif

}

extern "C"
void gpuSetAmoebaSASAParameters( amoebaGpuContext amoebaGpu, float probeRadius, 
                                 const std::vector<float>& radii, const std::vector<float>& weights )
{
   // ---------------------------------------------------------------------------------------

    static const char* methodName = "gpuSetAmoebaSASAParameters";
    static const int maxarc       = 500;

   // ---------------------------------------------------------------------------------------


    gpuContext gpu                         = amoebaGpu->gpuContext;
    amoebaGpu->paddedNumberOfAtoms         = amoebaGpu->gpuContext->sim.paddedNumberOfAtoms;
    unsigned int particles                 = radii.size();
    if( particles < 1 ){
        (void) fprintf( stderr, "%s no particles\n", methodName );
        return;
    } 
(void) fprintf( stderr, "%s radius converted Ang !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n", methodName );
float scaleRadius = 10.0f;

    amoebaGpu->amoebaSim.probeRadius       = probeRadius;
    amoebaGpu->amoebaSim.maxarc            = maxarc;
    amoebaGpu->psSASA_Radii                = new CUDAStream<float>(amoebaGpu->paddedNumberOfAtoms,   1, "SASARadii");
    amoebaGpu->psSASA_Weights              = new CUDAStream<float>(amoebaGpu->paddedNumberOfAtoms,   1, "SASAWeights");
    for (unsigned int ii = 0; ii < particles; ii++) 
    {    
        amoebaGpu->psSASA_Radii->_pSysStream[0][ii]     = scaleRadius*( radii[ii] + probeRadius );
        amoebaGpu->psSASA_Weights->_pSysStream[0][ii]   = weights[ii];
    }    

    // Dummy out extra particles data
    for (unsigned int ii = particles; ii < amoebaGpu->paddedNumberOfAtoms; ii++) 
    {    
        amoebaGpu->psSASA_Radii->_pSysStream[0][ii]     = 1.0f;
        amoebaGpu->psSASA_Weights->_pSysStream[0][ii]   = 0.0f;
    }    
    for (unsigned int ii = 0; ii < 4; ii++) 
    {
        amoebaGpu->psIntWorkArray[ii] = new CUDAStream<int>(amoebaGpu->amoebaSim.maxarc*amoebaGpu->paddedNumberOfAtoms,   1, "SASAIntWorkArray");
    }
    amoebaGpu->psFloatWorkArray = new CUDAStream<float>(amoebaGpu->amoebaSim.maxarc*amoebaGpu->paddedNumberOfAtoms,   1, "SASAFloatWorkArray");
    amoebaGpu->psIoListCount    = new CUDAStream<int>(amoebaGpu->paddedNumberOfAtoms,   1, "SASAIoCount");
    amoebaGpu->psDoneAtom       = new CUDAStream<int>(amoebaGpu->paddedNumberOfAtoms,   1, "SASADoneAtom");

    amoebaGpu->psSASA_Radii->Upload();
    amoebaGpu->psSASA_Weights->Upload();

#ifdef AMOEBA_DEBUG
    if( amoebaGpu->log ){
        unsigned int maxPrint = 32;
        (void) fprintf( amoebaGpu->log, "%s probeRadius=%12.3f\n", methodName, probeRadius );
        for (unsigned int ii = 0; ii < gpu->natoms; ii++) 
        {    
            (void) fprintf( amoebaGpu->log, "%5u %15.7e %15.7e\n", ii, radii[ii], weights[ii] );
            if( ii == maxPrint && ii < (amoebaGpu->paddedNumberOfAtoms - maxPrint) )
            {
                ii = (amoebaGpu->paddedNumberOfAtoms - maxPrint);
            }
        } 
        (void) fflush( amoebaGpu->log );
    }    
#endif

}

extern "C"
void amoebaGpuShutDown(amoebaGpuContext gpu)
{
    // free Cuda arrays

    delete gpu->psAmoebaBondID;
    delete gpu->psAmoebaBondParameter;

    delete gpu->psAmoebaAngleID1;
    delete gpu->psAmoebaAngleID2;
    delete gpu->psAmoebaAngleParameter;

    delete gpu->psAmoebaInPlaneAngleID1;
    delete gpu->psAmoebaInPlaneAngleID2;
    delete gpu->psAmoebaInPlaneAngleParameter;

    delete gpu->psAmoebaTorsionID1;
    delete gpu->psAmoebaTorsionID2;
    delete gpu->psAmoebaTorsionParameter1;
    delete gpu->psAmoebaTorsionParameter2;

    delete gpu->psAmoebaPiTorsionID1;
    delete gpu->psAmoebaPiTorsionID2;
    delete gpu->psAmoebaPiTorsionID3;
    delete gpu->psAmoebaPiTorsionParameter;

    delete gpu->psAmoebaStretchBendID1;
    delete gpu->psAmoebaStretchBendID2;
    delete gpu->psAmoebaStretchBendParameter;

    delete gpu->psAmoebaOutOfPlaneBendID1;
    delete gpu->psAmoebaOutOfPlaneBendID2;
    delete gpu->psAmoebaOutOfPlaneBendParameter;

    delete gpu->psAmoebaTorsionTorsionID1;
    delete gpu->psAmoebaTorsionTorsionID2;
    delete gpu->psAmoebaTorsionTorsionID3;
    delete gpu->psAmoebaTorsionTorsionGrids;

    // molecular frame multipoles

    delete gpu->psRotationMatrix;
    delete gpu->psMultipoleParticlesIdsAndAxisType;
    delete gpu->psMolecularDipole;
    delete gpu->psMolecularQuadrupole;
    delete gpu->psLabFrameDipole;
    delete gpu->psLabFrameQuadrupole;
    delete gpu->psDampingFactorAndThole;
    delete gpu->psCovalentDegree;
    delete gpu->psPolarizationDegree;
    delete gpu->psE_Field;
    delete gpu->psE_FieldPolar;
    delete gpu->psInducedDipole;
    delete gpu->psInducedDipolePolar;
    delete gpu->psPolarizability;
    delete gpu->psCurrentEpsilon;
    delete gpu->psWorkVector[0];
    delete gpu->psWorkVector[1];
    delete gpu->psWorkVector[2];
    delete gpu->psWorkVector[3];
    delete gpu->psForce;
    delete gpu->psTorque;
    delete gpu->psEnergy;
    delete gpu->torqueMapForce;

    delete gpu->psGk_Field;
    delete gpu->psInducedDipoleS;
    delete gpu->psInducedDipolePolarS;
    delete gpu->psBorn;
    delete gpu->psBornPolar;
    delete gpu->psKirkwoodForce;
    delete gpu->psKirkwoodEDiffForce;

    delete gpu->psVdwSigmaEpsilon;
    delete gpu->psVdwTable;
    delete gpu->psAmoebaVdwNonReductionID;
    delete gpu->psAmoebaVdwReductionID;
    delete gpu->psAmoebaVdwReduction;
    delete gpu->psAmoebaVdwCoordinates;
    delete gpu->psVdwWorkUnit;
    delete gpu->psVdwExclusionIndicesIndex;
    delete gpu->psVdwExclusionIndices;

    delete gpu->psWcaDispersionRadiusEpsilon;

    delete gpu->psSASA_Radii;
    delete gpu->psSASA_Weights;
    //delete gpu->psSASA_WeightIntWorkArray[0];
    //delete gpu->psSASA_WeightIntWorkArray[1];
    //delete gpu->psSASA_WeightIntWorkArray[2];
    //delete gpu->psSASA_WeightIntWorkArray[3];
    delete gpu->psDoneAtom;
    delete gpu->psIoListCount;
    delete gpu->psFloatWorkArray;

    delete gpu->psWorkArray_3_1; 
    delete gpu->psWorkArray_3_2; 
    delete gpu->psWorkArray_3_3; 
    delete gpu->psWorkArray_3_4; 
    delete gpu->psWorkArray_3_5; 
    delete gpu->psWorkArray_3_6; 

    delete gpu->psWorkArray_1_1; 
    delete gpu->psWorkArray_1_2; 

    delete gpu->psWorkUnit; 
    delete gpu->psScalingIndicesIndex; 
    delete gpu->ps_D_ScaleIndices; 
    delete gpu->ps_P_ScaleIndices; 
    delete gpu->ps_M_ScaleIndices; 

    if( gpu->pMapArray ){
        for( unsigned int ii = 0; ii < gpu->paddedNumberOfAtoms; ii++ ){
            delete gpu->pMapArray[ii];
        }
    }
    delete gpu->pMapArray;

    if( gpu->dMapArray ){
        for( unsigned int ii = 0; ii < gpu->paddedNumberOfAtoms; ii++ ){
            delete gpu->dMapArray[ii];
        }
    }
    delete gpu->dMapArray;

    // ---------------------------------------------------------------------------------------
    //delete gpu->amoebaSim;

    //gpuShutDown( gpu->gpuContext );
    delete gpu;

    return;
}

extern "C"
void amoebaGpuSetConstants(amoebaGpuContext amoebaGpu) 
{

    if( amoebaGpu->log ){
        (void) fprintf( amoebaGpu->log, "In amoebaGpuSetConstants\n" );
        (void) fflush( amoebaGpu->log );
    }

    if( amoebaGpu->amoebaSim.dielec > 0.0f && amoebaGpu->amoebaSim.dwater > 0.0f ){
        amoebaGpu->gpuContext->sim.preFactor = -amoebaGpu->amoebaSim.electric*((1.0f/amoebaGpu->amoebaSim.dielec)-(1.0f/amoebaGpu->amoebaSim.dwater));
    }

    gpuSetAmoebaBondOffsets( amoebaGpu );
    SetCalculateAmoebaLocalForcesSim( amoebaGpu );
    //SetCalculateAmoebaCudaSASAForcesSim( amoebaGpu );
    SetForcesSim( amoebaGpu->gpuContext );
    SetCalculateAmoebaMultipoleForcesSim( amoebaGpu );
    SetCalculateAmoebaCudaFixedEFieldSim( amoebaGpu );
    SetCalculateAmoebaCudaVdw14_7Sim( amoebaGpu );
    SetCalculateAmoebaCudaWcaDispersionSim( amoebaGpu );
    SetCalculateAmoebaCudaMutualInducedFieldSim( amoebaGpu );
    SetCalculateAmoebaElectrostaticSim( amoebaGpu );
    SetCalculateAmoebaCudaMapTorquesSim( amoebaGpu );
    SetCalculateAmoebaKirkwoodSim( amoebaGpu );
    SetCalculateAmoebaKirkwoodEDiffSim( amoebaGpu );
    SetCalculateAmoebaCudaFixedEAndGKFieldsSim( amoebaGpu );
    SetCalculateAmoebaCudaMutualInducedAndGkFieldsSim( amoebaGpu );
    //SetCalculateAmoebaObcGbsaForces2Sim( amoebaGpu );
    SetCalculateObcGbsaForces2Sim(  amoebaGpu->gpuContext  );
}

extern "C"
void amoebaGpuBuildOutputBuffers( amoebaGpuContext amoebaGpu )
{

    unsigned int paddedNumberOfAtoms                = amoebaGpu->paddedNumberOfAtoms;
    amoebaGpu->nonbondBlocks                        = amoebaGpu->gpuContext->sim.blocks;
    amoebaGpu->threadsPerBlock                      = amoebaGpu->gpuContext->sim.threads_per_block;

    // nonbondThreadsPerBlock & nonbondElectrostaticThreadsPerBlock need to be multiples of 32

    amoebaGpu->nonbondThreadsPerBlock               = 192;
    amoebaGpu->nonbondElectrostaticThreadsPerBlock  = 128;

    amoebaGpu->fieldReduceThreadsPerBlock           = (amoebaGpu->paddedNumberOfAtoms*3 + amoebaGpu->gpuContext->natoms + amoebaGpu->nonbondBlocks - 1) / amoebaGpu->nonbondBlocks;
    amoebaGpu->fieldReduceThreadsPerBlock           = ((amoebaGpu->fieldReduceThreadsPerBlock + (amoebaGpu->gpuContext->grid - 1)) /  amoebaGpu->gpuContext->grid) *  amoebaGpu->gpuContext->grid;

    if (amoebaGpu->fieldReduceThreadsPerBlock > amoebaGpu->nonbondThreadsPerBlock )
        amoebaGpu->fieldReduceThreadsPerBlock = amoebaGpu->nonbondThreadsPerBlock;

    if (amoebaGpu->fieldReduceThreadsPerBlock < 1) 
        amoebaGpu->fieldReduceThreadsPerBlock = 1; 

    // Select the number of output buffer to use.
    amoebaGpu->bOutputBufferPerWarp           = true;
    amoebaGpu->nonbondOutputBuffers           = (amoebaGpu->nonbondBlocks*amoebaGpu->nonbondThreadsPerBlock)/GRID;
    if (amoebaGpu->nonbondOutputBuffers >= (paddedNumberOfAtoms/GRID))
    {
        // For small systems, it is more efficient to have one output buffer per block of 32 atoms instead of one per warp.
        amoebaGpu->bOutputBufferPerWarp           = false;
        amoebaGpu->nonbondOutputBuffers           = paddedNumberOfAtoms/GRID;
    }
    amoebaGpu->outputBuffers           = amoebaGpu->nonbondOutputBuffers;
    //amoebaGpu->energyOutputBuffers     = max(amoebaGpu->nonbondThreadsPerBlock, amoebaGpu->sim.localForces_threads_per_block)*amoebaGpu->sim.blocks;
    amoebaGpu->energyOutputBuffers     = amoebaGpu->nonbondThreadsPerBlock*amoebaGpu->nonbondBlocks;
 
    if( amoebaGpu->log ){
        (void) fprintf(  amoebaGpu->log, "amoebaGpuBuildOutputBuffers: bOutputBufferPerWarp=%u nonbondBuffers=%u "
                         "outputBuffers=%u energyBuffers=%d nonbondBlocks=%u fieldReduceThreadsPerBlock=%u\n",
                         amoebaGpu->bOutputBufferPerWarp,
                         amoebaGpu->nonbondOutputBuffers,
                         amoebaGpu->outputBuffers,
                         amoebaGpu->energyOutputBuffers,
                         amoebaGpu->nonbondBlocks,
                         amoebaGpu->fieldReduceThreadsPerBlock ); 
        (void) fflush( amoebaGpu->log );
    }
    amoebaGpu->psWorkArray_3_1            = new CUDAStream<float>(3*paddedNumberOfAtoms, (amoebaGpu->outputBuffers), "AmoebaField_3_1");
    amoebaGpu->amoebaSim.pWorkArray_3_1   = amoebaGpu->psWorkArray_3_1->_pDevStream[0];

    amoebaGpu->psWorkArray_3_2            = new CUDAStream<float>(3*paddedNumberOfAtoms, (amoebaGpu->outputBuffers), "AmoebaField_3_2");
    amoebaGpu->amoebaSim.pWorkArray_3_2   = amoebaGpu->psWorkArray_3_2->_pDevStream[0];

    // used GK
    amoebaGpu->psWorkArray_3_3            = new CUDAStream<float>(3*paddedNumberOfAtoms, (amoebaGpu->outputBuffers), "AmoebaField_3_3");
    amoebaGpu->psWorkArray_3_4            = new CUDAStream<float>(3*paddedNumberOfAtoms, (amoebaGpu->outputBuffers), "AmoebaField_3_4");
    amoebaGpu->psWorkArray_3_5            = new CUDAStream<float>(3*paddedNumberOfAtoms, (amoebaGpu->outputBuffers), "AmoebaField_3_5");
    amoebaGpu->psWorkArray_3_6            = new CUDAStream<float>(3*paddedNumberOfAtoms, (amoebaGpu->outputBuffers), "AmoebaField_3_6");

    amoebaGpu->psWorkArray_1_1            = new CUDAStream<float>(  paddedNumberOfAtoms, (amoebaGpu->outputBuffers), "AmoebaField_1_1");
    amoebaGpu->amoebaSim.pWorkArray_1_1   = amoebaGpu->psWorkArray_1_1->_pDevStream[0];

    amoebaGpu->psWorkArray_1_2            = new CUDAStream<float>(  paddedNumberOfAtoms, (amoebaGpu->outputBuffers), "AmoebaField_1_2");
    amoebaGpu->amoebaSim.pWorkArray_1_2   = amoebaGpu->psWorkArray_1_2->_pDevStream[0];

    amoebaGpu->psEnergy                   = new CUDAStream<float>(amoebaGpu->energyOutputBuffers, 1, "AmoebaEnergy");

    return;
}

static int matchMaps( std::string idString, MapIntFloat* map1,  MapIntFloat* map2, FILE* log ){

    int equalSizes = 1;
    int error      = 0;
    if( map1->size() != map2->size() ){
        (void) fprintf( log, "%s sizes unequal: %u %u\n", idString.c_str(), static_cast<unsigned int>(map1->size()), static_cast<unsigned int>(map2->size()) );
        equalSizes = 0;
        error++;
    }
    for( MapIntFloatCI ii = map2->begin(); ii != map2->end(); ii++ ){
        int hit   = 0;
        int idHit = 0;
        for( MapIntFloatCI jj = map1->begin(); jj != map1->end() && hit == 0; jj++ ){
            if( (*ii).first == (*jj).first ){
                idHit = 1;
                if( (*ii).second == (*jj).second ){
                    hit = 1;
                }
            }
        }
        if( hit == 0 ){
            error++;
            (void) fprintf( log, "%s (1::2) no match for %d %.3f idHit=%d\n", idString.c_str(), (*ii).first, (*ii).second, idHit );
        }
    }
    for( MapIntFloatCI ii = map1->begin(); ii != map1->end(); ii++ ){
        int hit   = 0;
        int idHit = 0;
        for( MapIntFloatCI jj = map2->begin(); jj != map2->end() && hit == 0; jj++ ){
            if( (*ii).first == (*jj).first ){
                idHit = 1;
                if( (*ii).second == (*jj).second ){
                    hit = 1;
                }
            }
        }
        if( hit == 0 ){
            error++;
            (void) fprintf( log, "%s (2::1) no match for %d %.3f idHit=%d\n", idString.c_str(), (*ii).first, (*ii).second, idHit );
        }
    }

    if( error ){
        (void) fprintf( log, "Map1 [" );
        for( MapIntFloatCI ii = map1->begin(); ii != map1->end(); ii++ ){
            (void) fprintf( log, "%d %.3f  ", (*ii).first, (*ii).second );
        }
        (void) fprintf( log, "]\n" );
        (void) fprintf( log, "Map2 [" );
        for( MapIntFloatCI ii = map2->begin(); ii != map2->end(); ii++ ){
            (void) fprintf( log, "%d %.3f  ", (*ii).first, (*ii).second );
        }
        (void) fprintf( log, "]\n" );
    }
    return error;
}

static void getScalingDegrees( amoebaGpuContext amoebaGpu, unsigned int particleI, unsigned int particleJ, int* covalentDegree, int* polarizationDegree )
{
    int particlesOffset                        = particleI*amoebaGpu->maxCovalentDegreeSz;

    unsigned int minCovalentIndex              = static_cast<unsigned int>(amoebaGpu->psCovalentDegree->_pSysStream[0][particlesOffset]);
    unsigned int minCovalentPolarizationIndex  = static_cast<unsigned int>( amoebaGpu->psPolarizationDegree->_pSysStream[0][particlesOffset]);

    if( particleJ < minCovalentIndex || particleJ > (minCovalentIndex + amoebaGpu->maxCovalentDegreeSz) ){
        *covalentDegree     = 0;
    } else {
        *covalentDegree     = amoebaGpu->psCovalentDegree->_pSysStream[0][particlesOffset + (particleJ-minCovalentIndex) + 1];
    }

    if( particleJ < minCovalentPolarizationIndex || particleJ > (minCovalentPolarizationIndex + amoebaGpu->maxCovalentDegreeSz) ){
        *polarizationDegree = 0;
    } else {
        *polarizationDegree = amoebaGpu->psPolarizationDegree->_pSysStream[0][particlesOffset + (particleJ-minCovalentPolarizationIndex) + 1];
    }
}

extern "C"
int amoebaGpuBuildThreadBlockWorkList( amoebaGpuContext amoebaGpu )
{
    if( amoebaGpu->psWorkUnit != NULL ){
        return 0;
    }
    const unsigned int atoms = amoebaGpu->paddedNumberOfAtoms;
    const unsigned int grid  = amoebaGpu->gpuContext->grid;
    const unsigned int dim   = (atoms + (grid - 1)) / grid;
    const unsigned int cells = dim * (dim + 1) / 2;

    CUDAStream<unsigned int>* psWorkUnit       = new CUDAStream<unsigned int>(cells, 1u, "WorkUnit");
    unsigned int* pWorkList                    = psWorkUnit->_pSysData;
    amoebaGpu->psWorkUnit                      = psWorkUnit;

    CUDAStream<unsigned int>* psVdwWorkUnit    = new CUDAStream<unsigned int>(cells, 1u, "VdwWorkUnit");
    unsigned int* pVdwWorkList                 = psVdwWorkUnit->_pSysData;
    amoebaGpu->psVdwWorkUnit                   = psVdwWorkUnit;

/*
    CUDAStream<unsigned int>* psInteractingWorkUnit  = new CUDAStream<unsigned int>(cells, 1u, "InteractingWorkUnit");
    amoebaGpu->psInteractingWorkUnit                 = psInteractingWorkUnit;
    amoebaGpu->workUnits                             = cells;

    CUDAStream<unsigned int>* psInteractionFlag = new CUDAStream<unsigned int>(cells, 1u, "InteractionFlag");
    amoebaGpu->psInteractionFlag = psInteractionFlag;
    amoebaGpu->sim.pInteractionFlag = psInteractionFlag->_pDevStream[0];
    CUDAStream<size_t>* psInteractionCount = new CUDAStream<size_t>(1, 1u, "InteractionCount");
    amoebaGpu->psInteractionCount = psInteractionCount;
    amoebaGpu->sim.pInteractionCount = psInteractionCount->_pDevStream[0];
    CUDAStream<float4>* psGridBoundingBox = new CUDAStream<float4>(dim, 1u, "GridBoundingBox");
    amoebaGpu->psGridBoundingBox = psGridBoundingBox;
    amoebaGpu->sim.pGridBoundingBox = psGridBoundingBox->_pDevStream[0];
    CUDAStream<float4>* psGridCenter = new CUDAStream<float4>(dim, 1u, "GridCenter");
    amoebaGpu->psGridCenter = psGridCenter;
    amoebaGpu->sim.pGridCenter = psGridCenter->_pDevStream[0];

    amoebaGpu->sim.nonbond_workBlock      = amoebaGpu->sim.nonbondThreadsPerBlock / GRID;
    amoebaGpu->sim.bornForce2_workBlock   = amoebaGpu->sim.bornForce2_threads_per_block / GRID;
    amoebaGpu->sim.workUnits = cells;
*/

    // Initialize the plan for doing stream compaction.

    // planCompaction(amoebaGpu->compactPlan);

    // Increase block count if necessary for extra large molecules that would
    // otherwise overflow the SM workunit buffers
//    int minimumBlocks = (cells + amoebaGpu->sim.workUnitsPerSM - 1) / amoebaGpu->sim.workUnitsPerSM;
//    if ((int) amoebaGpu->sim.nonbond_blocks < minimumBlocks)
//    {
//        amoebaGpu->sim.nonbond_blocks = amoebaGpu->sim.nonbond_blocks * ((minimumBlocks + amoebaGpu->sim.nonbond_blocks - 1) / amoebaGpu->sim.nonbond_blocks);
//    }
//    if ((int) amoebaGpu->sim.bornForce2_blocks < minimumBlocks)
//    {
//        amoebaGpu->sim.bornForce2_blocks = amoebaGpu->sim.bornForce2_blocks * ((minimumBlocks + amoebaGpu->sim.bornForce2_blocks - 1) / amoebaGpu->sim.bornForce2_blocks);
//    }


/*
    amoebaGpu->sim.nbWorkUnitsPerBlock            = cells / amoebaGpu->sim.nonbond_blocks;
    amoebaGpu->sim.nbWorkUnitsPerBlockRemainder   = cells - amoebaGpu->sim.nonbond_blocks * amoebaGpu->sim.nbWorkUnitsPerBlock;
    amoebaGpu->sim.interaction_threads_per_block  = 64;

    amoebaGpu->sim.interaction_blocks = (amoebaGpu->workUnits + amoebaGpu->sim.interaction_threads_per_block - 1) / amoebaGpu->sim.interaction_threads_per_block;
    if (amoebaGpu->sim.interaction_blocks > 8*amoebaGpu->sim.blocks)
        amoebaGpu->sim.interaction_blocks = 8*amoebaGpu->sim.blocks;
    if (activeWorkUnits > (int) cells)
    {
        int balancedWorkBlock                   = (cells + amoebaGpu->sim.nonbond_blocks - 1) / amoebaGpu->sim.nonbond_blocks;
        amoebaGpu->sim.nonbondThreadsPerBlock      = balancedWorkBlock * GRID;
        amoebaGpu->sim.nonbond_workBlock              = balancedWorkBlock;
    }
    activeWorkUnits = amoebaGpu->sim.bornForce2_blocks * amoebaGpu->sim.bornForce2_workBlock;
    if (activeWorkUnits > (int) cells)
    {
        int balancedWorkBlock                   = (cells + amoebaGpu->sim.bornForce2_blocks - 1) / amoebaGpu->sim.bornForce2_blocks;
        amoebaGpu->sim.bornForce2_threads_per_block   = balancedWorkBlock * GRID;
        amoebaGpu->sim.bornForce2_workBlock           = balancedWorkBlock;
    }
*/

    unsigned int count = 0;
    for (unsigned int y = 0; y < dim; y++)
    {
        for (unsigned int x = y; x < dim; x++)
        {
            pWorkList[count]    = encodeCell( x, y );
            pVdwWorkList[count] = encodeCell( x, y );
            count++;
        }
    }
    //(*amoebaGpu->psInteractionCount)[0] = amoebaGpu->workUnits;
    //amoebaGpu->psInteractionCount->Upload();

    psWorkUnit->Upload();
    psVdwWorkUnit->Upload();

    if( amoebaGpu->log ){
        (void) fprintf( amoebaGpu->log, "amoebaGpuBuildThreadBlockWorkList %u %u dim=%u cells=%u nonbondThreadsPerBlock=%u nonbond_blocks=%u\n",
                                        atoms, grid, dim, cells, amoebaGpu->nonbondThreadsPerBlock, amoebaGpu->nonbondBlocks );
        (void) fflush( amoebaGpu->log );
    }

    return cells;
}

extern "C"
void amoebaGpuBuildScalingList( amoebaGpuContext amoebaGpu )
{
    // ---------------------------------------------------------------------------------------

    static const std::string methodName = "amoebaGpuBuildScalingList";
    static const int debugOn            = 0;

    // ---------------------------------------------------------------------------------------

    if( amoebaGpu->psCovalentDegree == NULL ){
        return;
    }    

    const unsigned int paddedAtoms     = amoebaGpu->paddedNumberOfAtoms;
    const unsigned int actualAtoms     = amoebaGpu->gpuContext->natoms;
    const unsigned int grid            = amoebaGpu->gpuContext->grid;
    const unsigned int dim             = paddedAtoms/grid;
    const unsigned int cells           = dim * (dim + 1) / 2;
    unsigned int* pWorkList            = amoebaGpu->psWorkUnit->_pSysData;

    // minCellIndex & maxCellIndex track min/max atom index for each cell

    std::vector<int> minCellIndex;
    std::vector<int> maxCellIndex;
    minCellIndex.resize( dim + 1 );
    maxCellIndex.resize( dim + 1 );

    for ( unsigned int ii = 0; ii <= dim; ii++)
    {
        minCellIndex[ii] = paddedAtoms + 1;
        maxCellIndex[ii] = 0;
    }

    for (unsigned int atom1 = 0; atom1 < paddedAtoms; atom1++)
    {
        int x                  = atom1/grid;
        int particlesOffset    = atom1*amoebaGpu->maxCovalentDegreeSz;
        int minCovalentIndex   = amoebaGpu->psCovalentDegree->_pSysStream[0][particlesOffset];
        int minPolarCovIndex   = amoebaGpu->psPolarizationDegree->_pSysStream[0][particlesOffset];
        int maxCIndex          = 0;
        int maxPIndex          = 0;
        for (int jj = amoebaGpu->maxCovalentDegreeSz - 1; jj >= 1 && (maxPIndex == 0 || maxCIndex == 0); jj-- )
        {
            if( amoebaGpu->psCovalentDegree->_pSysStream[0][particlesOffset+jj] && maxCellIndex[x] < (minCovalentIndex+jj) )
            { 
                maxCellIndex[x] =  minCovalentIndex + jj;
                maxCIndex++; 
            }
            if( amoebaGpu->psPolarizationDegree->_pSysStream[0][particlesOffset+jj] && maxCellIndex[x] < (minPolarCovIndex+jj) )
            { 
                maxCellIndex[x] =  minPolarCovIndex + jj;
                maxPIndex++; 
            }
        }
        if( minCellIndex[x] > minCovalentIndex ){ 
            minCellIndex[x] = minCovalentIndex;
        }
        if( minCellIndex[x] > minPolarCovIndex ){ 
            minCellIndex[x] = minPolarCovIndex;
        }
    }

    // diagnostics

    if( debugOn && amoebaGpu->log ){
        (void) fprintf( amoebaGpu->log, "%s min/max cell indices:\n", methodName.c_str() );
        for (unsigned int ii = 0; ii < dim; ii++)
        {
            (void) fprintf( amoebaGpu->log, "%6d [%6d %6d]\n", ii, minCellIndex[ii], maxCellIndex[ii] );
        }
        (void) fflush( amoebaGpu->log );
    }

    // Build a list of indexes for the work units with scaling values different from 1

    CUDAStream<int>* psScalingIndicesIndex          = new CUDAStream<int>(cells, 1u, "ScalingIndicesIndex");
    amoebaGpu->psScalingIndicesIndex                = psScalingIndicesIndex;
    amoebaGpu->amoebaSim.pScaleIndicesIndex         = psScalingIndicesIndex->_pDevStream[0];

    memset( amoebaGpu->psScalingIndicesIndex->_pSysStream[0], 0, sizeof( cells)*sizeof( unsigned int) );
    int numWithScalingIndices                       = 0;
    int gridOffset                                  = grid - 1;
    int lastBlock                                   = (static_cast<int>(paddedAtoms) > amoebaGpu->gpuContext->natoms) ? (amoebaGpu->gpuContext->natoms)/grid : -1;
    for (unsigned int ii = 0; ii < cells; ii++)
    {
        unsigned int x, y, exclusion;
        decodeCell( pWorkList[ii], &x, &y, &exclusion );
        int xAtomMin = x*grid;
        int xAtomMax = xAtomMin + gridOffset;
        if( (maxCellIndex[y] >= xAtomMin && minCellIndex[y] <= xAtomMax) || (x == lastBlock || y == lastBlock) ){
            pWorkList[ii]                              = encodeCellExclusion( pWorkList[ii] );
            psScalingIndicesIndex->_pSysStream[0][ii]  = numWithScalingIndices*grid;
            numWithScalingIndices++;
//(void) fprintf( amoebaGpu->log, "%5d [%6d %6d] [%6d %6d] [%6d %6d] num=%5d last=%5d\n",
//                ii, x, y, xAtomMin, xAtomMax, minCellIndex[y], maxCellIndex[y], numWithScalingIndices, lastBlock );
        }
    }

    // diagnostics

#if 0
    if( 0 && amoebaGpu->log ){
        (void) fprintf( amoebaGpu->log, "%s %d cells\n",
                                        methodName.c_str(), numWithScalingIndices );
        for (int ii = 0; ii < cells; ii++)
        {
            unsigned int x, y, exclusion;
            decodeCell( pWorkList[ii], &x, &y, &exclusion );
            if( exclusion ){
                (void) fprintf( amoebaGpu->log, "%6d [%6u %6u] %8u %8u [%6u %6u]indexInToIndices=%8d\n",
                                ii, x, y, exclusion, pWorkList[ii], 32*x, 32*y,
                                psScalingIndicesIndex->_pSysStream[0][ii] );
            } else {
                (void) fprintf( amoebaGpu->log, "%6d [%6u %6u] %8u %8u [%6u %6u]\n",
                                ii, x, y, exclusion, pWorkList[ii], 32*x, 32*y );
            }
        }
        (void) fflush( amoebaGpu->log );
    }
#else
    if( debugOn && amoebaGpu->log ){
        (void) fprintf( amoebaGpu->log, "%s %d cells w/ exclusions\n",
                                        methodName.c_str(), numWithScalingIndices );
        for (unsigned int ii = 0; ii < cells; ii++)
        {
            unsigned int x, y, exclusion;
            decodeCell( pWorkList[ii], &x, &y, &exclusion );
            if( exclusion ){
                (void) fprintf( amoebaGpu->log, "%6d [%6u %6u] %8u %8u indexInToIndices=%8d\n", ii, x, y, exclusion, pWorkList[ii],
                                psScalingIndicesIndex->_pSysStream[0][ii] );
            }
        }
        (void) fflush( amoebaGpu->log );
    }
#endif
    // Record the scaling data

    CUDAStream<int>* ps_D_ScaleIndices           = new CUDAStream<int>(numWithScalingIndices*grid, 1u, "ps_D_ScaleIndices");
    amoebaGpu->ps_D_ScaleIndices                 = ps_D_ScaleIndices;
    amoebaGpu->amoebaSim.pD_ScaleIndices         = ps_D_ScaleIndices->_pDevStream[0];

    CUDAStream<int2>* ps_P_ScaleIndices          = new CUDAStream<int2>(numWithScalingIndices*grid, 1u, "ps_P_ScaleIndices");
    amoebaGpu->ps_P_ScaleIndices                 = ps_P_ScaleIndices;
    amoebaGpu->amoebaSim.pP_ScaleIndices         = ps_P_ScaleIndices->_pDevStream[0];

    CUDAStream<int2>* ps_M_ScaleIndices          = new CUDAStream<int2>(numWithScalingIndices*grid, 1u, "ps_M_ScaleIndices");
    amoebaGpu->ps_M_ScaleIndices                 = ps_M_ScaleIndices;
    amoebaGpu->amoebaSim.pM_ScaleIndices         = ps_M_ScaleIndices->_pDevStream[0];

    memset( ps_D_ScaleIndices->_pSysStream[0], 0, numWithScalingIndices*grid*sizeof( int )   );
    memset( ps_P_ScaleIndices->_pSysStream[0], 0, numWithScalingIndices*grid*sizeof( int )*2 );
    memset( ps_M_ScaleIndices->_pSysStream[0], 0, numWithScalingIndices*grid*sizeof( int )*2 );

    // load scaling indices

    // psScalingIndicesIndex[cell] gives the index into array of scalingMasks (ps_X_ScaleIndices)
    // for that cell

    // for each cell, ps_X_ScaleIndices is an int4 array of GRID(32) masks for that cell
    // ps_X_ScaleIndices[scaleIndex + threadIdI] is the mask for atomI (x + threadIdI) w/ the the y-atoms 
    // the ith bits of the mask give the scale index for atomI w/ atom (y + i)
 
static unsigned int targetAtoms[2] = { 0, 1};
    for ( unsigned int ii = 0; ii < cells; ++ii)
    {
        unsigned int x, y, exclusion;
        decodeCell( pWorkList[ii], &x, &y, &exclusion );
        if ( exclusion )
        {
            int scaleIndex = psScalingIndicesIndex->_pSysStream[0][ii];  
            for (unsigned int jj = 0; jj < grid; jj++)
            {
                unsigned int atomI = grid*x + jj;
                int scaleOffset    = scaleIndex + jj;
                for (unsigned int kk = 0; kk < grid; kk++)
                {
                    int covalentDegree;
                    int polarizationDegree;
                    unsigned int atomJ = grid*y + kk;
                    getScalingDegrees( amoebaGpu, atomI, atomJ, &covalentDegree, &polarizationDegree ); 
       
                    int pX, pY;
                    // polarScale[5]   = { 0.0f, 0.0f, 0.0f, 1.0f, 1.0f };
                    // pScale[3]       = { 1.0f, 0.5f, 0.0f };
                    if( (atomI == atomJ) || (atomI >= actualAtoms) || (atomJ >= actualAtoms) ){ 
 
                        // 0.0

                        pX  = 0;
                        pY  = 1;

                    } else if( covalentDegree == 0 ){

                        // 1.0

                        pX  = 0;
                        pY  = 0;

                    } else if( covalentDegree <= 3 ){ 
 
                        // 0.0

                        pX  = 0;
                        pY  = 1;

                    } else if( covalentDegree == 4 && polarizationDegree == 1 ){ 
 
                        // 0.5

                        pX  = 1;
                        pY  = 0;

                    } else {

                        // 1.0

                        pX  = 0;
                        pY  = 0;
                    }

                    if( pX ){
                        ps_P_ScaleIndices->_pSysStream[0][scaleOffset].x |= (pX << kk);
                    }
                    if( pY ){
                        ps_P_ScaleIndices->_pSysStream[0][scaleOffset].y |= (pY << kk);
                    }

                    // directScale[5]  = { 0.0f, 1.0f, 1.0f, 1.0f, 1.0f };
                    // dScale[2]       = { 1.0f, 0.0f };
                    if( (polarizationDegree == 1) || (atomI >= actualAtoms) || (atomJ >= actualAtoms) ){
                        ps_D_ScaleIndices->_pSysStream[0][scaleOffset] |= (1 << kk);
                    }

                    int mX = 0;
                    int mY = 0;
                    // mpoleScale[5]   = { 0.0f, 0.0f, 0.0f, 0.4f, 0.8f };
                    // mScale[4]       = { 1.0f, 0.4f, 0.8f, 0.0f };
                    if( (atomI == atomJ) || (atomI >= actualAtoms) || (atomJ >= actualAtoms) ||  (covalentDegree > 0 && covalentDegree <= 3) ){ 
 
                        // 0.0

                        mX  = 1;
                        mY  = 1;

                    } else if( covalentDegree == 4 ){ 
 
                        // 0.4

                        mX  = 0;
                        mY  = 1;

                    } else if( covalentDegree == 5 ){ 
 
                        // 0.8

                        mX  = 1;
                        mY  = 0;

                    }

                    if( mX ){
                        ps_M_ScaleIndices->_pSysStream[0][scaleOffset].x |= (1 << kk);
                    }
                    if( mY ){
                        ps_M_ScaleIndices->_pSysStream[0][scaleOffset].y |= (1 << kk);
                    }

                    if( 0 && amoebaGpu->log && ( (atomI == targetAtoms[0]) || (atomI == targetAtoms[1]) ) ){
                        (void) fprintf( amoebaGpu->log, "XXX cell=%u [%u %u] [%d %d] p[%d %d] m[%d %d] scaleOffset=%d kk=%d\n", 
                                        ii, atomI, atomJ, covalentDegree, polarizationDegree, pX, pY, mX, mY, scaleOffset, kk );
                    }
                }
            }
        }
    }

    // diagnostics

    if( debugOn && amoebaGpu->log ){

        float* pScaleCheckSum = (float*) malloc( sizeof( float )*paddedAtoms );
        float* dScaleCheckSum = (float*) malloc( sizeof( float )*paddedAtoms );
        float* mScaleCheckSum = (float*) malloc( sizeof( float )*paddedAtoms );
    
        memset( pScaleCheckSum, 0, paddedAtoms*sizeof( float ) );
        memset( dScaleCheckSum, 0, paddedAtoms*sizeof( float ) );
        memset( mScaleCheckSum, 0, paddedAtoms*sizeof( float ) );

        MapIntFloat** pMapArray = (MapIntFloat**) malloc( sizeof( MapIntFloat* )*paddedAtoms );
        MapIntFloat** dMapArray = (MapIntFloat**) malloc( sizeof( MapIntFloat* )*paddedAtoms );
        for( unsigned int ii = 0; ii < paddedAtoms; ii++ ){
            pMapArray[ii] =  new MapIntFloat;
            dMapArray[ii] =  new MapIntFloat;
        }

        float pScale[4]  = { 1.0f, 0.5f, 0.0f, -1.0f };
        (void) fprintf( amoebaGpu->log, "%s Pscale values\n",
                                        methodName.c_str(), numWithScalingIndices );
        for (unsigned int ii = 0; ii < cells; ii++)
        {
            unsigned int x, y, exclusion;
            decodeCell( pWorkList[ii], &x, &y, &exclusion );
            if( exclusion ){
                int scaleIndex     = psScalingIndicesIndex->_pSysStream[0][ii];  
                for (unsigned int jj = 0; jj < grid; jj++)
                {
                    unsigned int atomI        = grid*x + jj;
                    unsigned int scaleOffset  = scaleIndex + jj;
                    for (unsigned int kk = 0; kk < grid; kk++)
                    {
                        int covalentDegree;
                        int polarizationDegree;
                        unsigned int atomJ         = grid*y + kk;
                        unsigned int mask          = 1 << kk;
                        unsigned int valueP_X      = ps_P_ScaleIndices->_pSysStream[0][scaleOffset].x & mask;
                        unsigned int valueP_Y      = ps_P_ScaleIndices->_pSysStream[0][scaleOffset].y & mask;
                        unsigned int valueM_X      = ps_M_ScaleIndices->_pSysStream[0][scaleOffset].x & mask;
                        unsigned int valueM_Y      = ps_M_ScaleIndices->_pSysStream[0][scaleOffset].y & mask;
                        unsigned int valueD        = ps_D_ScaleIndices->_pSysStream[0][scaleOffset]   & mask;
                        unsigned int valueI        = 0;
                        if( valueP_X ){
                            valueI++;
                        }
                        if( valueP_Y ){
                            valueI += 2;
                        }
                        float dScale = valueD ? 0.0f : 1.0f;

                        if( pScale[valueI] != 1.0f ){
                            MapIntFloat* pMap = pMapArray[atomI];
                            (*pMap)[atomJ]    = pScale[valueI];
                            pMap              = pMapArray[atomJ];
                            (*pMap)[atomI]    = pScale[valueI];
                        }
                        if( atomI < paddedAtoms ){
                            pScaleCheckSum[atomI] += (pScale[valueI] - 1.0f);
                            dScaleCheckSum[atomI] += (dScale - 1.0f);
                        }
                        if( pScale[valueI] != 1.0f || dScale != 1.0f ){
                             getScalingDegrees( amoebaGpu, atomI, atomJ, &covalentDegree, &polarizationDegree ); 
                             (void) fprintf( amoebaGpu->log, "%6d cell[%6u %6u] %1u atom[%6u %6u] deg[%6u %6u] %6u %6u %2u %.3f %.3f [0]=%u %u\n",
                                             ii, x, y, exclusion,
                                             atomI, atomJ, covalentDegree, polarizationDegree, 
                                             (valueP_X ? 1 : 0), (valueP_Y ? 1 : 0), valueI, pScale[valueI], dScale,
                                             (ps_P_ScaleIndices->_pSysStream[0][0].x & 1),
                                             (ps_P_ScaleIndices->_pSysStream[0][0].y & 1)
                                             );
                        }

                    }
                }
            }
        }
        (void) fflush( amoebaGpu->log );

#if 1
        unsigned int totalWarps      = (amoebaGpu->nonbondBlocks*amoebaGpu->nonbondThreadsPerBlock)/GRID;
        //unsigned int warp            = (blockIdx.x*blockDim.x+threadIdx.x)/GRID;
        //unsigned int numWorkUnits    = cSim.pInteractionCount[0];
        unsigned int numWorkUnits    = 780;
    
        (void) fprintf( amoebaGpu->log, "XX Total warps=%u blocks=%u threads=%u GRID=%u\n",
                        totalWarps, amoebaGpu->nonbondBlocks, amoebaGpu->nonbondThreadsPerBlock, GRID );
        unsigned int maxPrint = 3;
        std::stringstream message;
        char buffer[2048];
        unsigned int targetAtom = 0;
        for( unsigned int ii = 0; ii < amoebaGpu->nonbondBlocks; ii++ )
        {
            unsigned int blockId = ii;
            for( unsigned int jj = 0; jj < amoebaGpu->nonbondThreadsPerBlock; jj++ )
            {
					 unsigned int warp = (ii*amoebaGpu->nonbondThreadsPerBlock+jj)/GRID;
                unsigned int pos  = warp*numWorkUnits/totalWarps;
                unsigned int end  = (warp+1)*numWorkUnits/totalWarps;
                (void) sprintf( buffer, "Block %4u thread %4u warp=%4u pos[%4u %4u]\n",
                                ii, jj, warp, pos, end );
                unsigned int print = 0;
                while( pos < end ){
                    unsigned int x, y, exclusion;
                    decodeCell( pWorkList[pos], &x, &y, &exclusion );
                    x                   *= GRID;
                    y                   *= GRID;
                    unsigned int tgx     = jj & (GRID - 1); 
                    unsigned int tbx     = jj - tgx;
                    unsigned int tj      = tgx;
tgx     = 0;
                    unsigned int atomI   = x + tgx;
                    unsigned int offset1 = (x + tgx + (y >> GRIDBITS) * amoebaGpu->paddedNumberOfAtoms);
                    unsigned int offset2 = (y + tgx + (x >> GRIDBITS) * amoebaGpu->paddedNumberOfAtoms);
                    unsigned int offset3 = x + tgx;
                    unsigned int offset4 = y + tgx;
                    unsigned int offset5 = (x >> GRIDBITS);
                    unsigned int offset6 = (y >> GRIDBITS);
                    if( (x <= targetAtom && targetAtom < (x+32)) ||
                        (y <= targetAtom && targetAtom < (y+32)) ){
                        if( print == 0 ){
                            print++;
                            message << buffer;
                        }
                        (void) sprintf( buffer, "   pos=%3u atomI=%4u tgx=%4u tbx=%4u tj=%4u x/y[%4u %4u]"
                                        " scl=%u off[%u %u] [%u %u]",
                                        pos, atomI, tgx, tbx, tj, x, y, exclusion, offset3, offset4, offset5, offset6 );
                        message << buffer;
                        if( exclusion ){
                            unsigned int xi        = x >> GRIDBITS;
                            unsigned int yi        = y >> GRIDBITS;
                            unsigned int cell      = xi+yi*amoebaGpu->paddedNumberOfAtoms/GRID-yi*(yi+1)/2;
                            unsigned int cellIndex = psScalingIndicesIndex->_pSysStream[0][cell]+tgx;
                            int  dScaleMask        = ps_D_ScaleIndices->_pSysStream[0][cellIndex];
                            int2 pScaleMask        = ps_P_ScaleIndices->_pSysStream[0][cellIndex];
                            (void) sprintf( buffer, "   xi=%u yi=%u cell=%u cellIndex=%u",
                                           xi, yi, cell, cellIndex );
                            message << buffer;
         
                            if( (x <= targetAtom && targetAtom < (x+32)) ){
                                unsigned int offset    = targetAtom - x;
                                unsigned int mask      =  1 << offset;
                                (void) sprintf( buffer, "   off=%u [%u %u %u] mask=%u", offset,
                                               (pScaleMask.x & mask) ? 1 : 0,
                                               (pScaleMask.y & mask) ? 1 : 0,
                                               (dScaleMask   & mask) ? 1 : 0, mask );
                                message << buffer;
                            }
                        }
                        message << std::endl;
                    }
                    pos++;
                }
#if 0
                if( jj == maxPrint && (maxPrint < amoebaGpu->nonbondThreadsPerBlock-jj) ){
                    jj = amoebaGpu->nonbondThreadsPerBlock - maxPrint;
                    (void) fprintf( amoebaGpu->log, "\n\n" );
                }
#endif
            }
        }
        (void) fprintf( amoebaGpu->log, "%s\n\n", message.str().c_str() );
#endif

#if 0
        (void) fprintf( amoebaGpu->log, "ZXX Total warps=%u blocks=%u threads=%u GRID=%u\n",
                        amoebaGpu->nonbondBlocks, amoebaGpu->fieldReduceThreadsPerBlock );
        unsigned int maxPrint = 3;
        for( unsigned int ii = 0; ii < amoebaGpu->nonbondBlocks; ii++ )
        {
            unsigned int blockId = ii;
            for( unsigned int jj = 0; jj < amoebaGpu->fieldReduceThreadsPerBlock; jj++ )
            {
					 unsigned int warp = (ii*amoebaGpu->nonbondThreadsPerBlock+jj)/GRID;
                unsigned int pos  = warp*numWorkUnits/totalWarps;
                unsigned int end  = (warp+1)*numWorkUnits/totalWarps;
                (void) fprintf( amoebaGpu->log, "Block %4u thread %4u warp=%4u pos[%4u %4u]\n",
                                ii, jj, warp, pos, end );
                while( pos < end ){
                    unsigned int x, y, exclusion;
                    decodeCell( pWorkList[pos], &x, &y, &exclusion );
                    x *= GRID;
                    y *= GRID;
                    unsigned int tgx     = jj & (GRID - 1); 
                    unsigned int tbx     = jj - tgx;
                    unsigned int tj      = tgx;
                    unsigned int atomI   = x + tgx;
                    unsigned int offset1 = (x + tgx + (y >> GRIDBITS) * amoebaGpu->paddedNumberOfAtoms);
                    unsigned int offset2 = (y + tgx + (x >> GRIDBITS) * amoebaGpu->paddedNumberOfAtoms);
                    unsigned int offset3 = x + tgx;
                    unsigned int offset4 = y + tgx;
                    unsigned int offset5 = (y >> GRIDBITS);
                    unsigned int offset6 = (x >> GRIDBITS);
                    (void) fprintf( amoebaGpu->log, "   atomI=%4u tgx=%4u tbx=%4u tj=%4u x/y[%4u %4u]"
                                    " scl=%u off[%6u %6u] [%6u %6u]\n",
                                    atomI, tgx, tbx, tj, x, y, exclusion, offset3, offset4, offset5, offset6 );
                    pos++;
                }
                if( jj == maxPrint && (maxPrint < amoebaGpu->nonbondThreadsPerBlock-jj) ){
                    jj = amoebaGpu->nonbondThreadsPerBlock - maxPrint;
                    (void) fprintf( amoebaGpu->log, "\n\n" );
                }
            }
        }
#endif

#if 0
        FILE* filePtr = fopen( "newScale.txt", "w" );
        for( unsigned int kk = 0; kk < actualAtoms; kk++ ){
            (void) fprintf( filePtr, "%6u %14.6e %14.6e\n", kk, pScaleCheckSum[kk], dScaleCheckSum[kk] );
        }
        (void) fclose( filePtr );
        free( pScaleCheckSum );
        free( dScaleCheckSum );
        free( mScaleCheckSum );
        filePtr = fopen( "newScaleMap.txt", "w" );
        //char buffer[1024];

        for( unsigned int kk = 0; kk < actualAtoms; kk++ ){
            MapIntFloat* pMap1 = amoebaGpu->pMapArray[kk];
            MapIntFloat* pMap2 = pMapArray[kk];

            float sum = 0.0f;
            for( MapIntFloatCI ii = pMap2->begin(); ii != pMap2->end(); ii++ ){
                sum += (*ii).second > 0.0f ? 0.5f : 1.0f;
            }
            (void) fprintf( filePtr, "%6u sz=%u sum=%.3f total=%.3f %.3f\n", kk, pMap2->size(),
                            sum, (float) actualAtoms - sum, 1248.0f - sum );
            for( MapIntFloatCI ii = pMap2->begin(); ii != pMap2->end(); ii++ ){
                (void) fprintf( filePtr, "    %6d %14.6e\n", (*ii).first, (*ii).second );
            }
        }
        for( unsigned int kk = 0; kk < actualAtoms; kk++ ){
            MapIntFloat* pMap1 = amoebaGpu->pMapArray[kk];
            MapIntFloat* pMap2 = pMapArray[kk];
            (void) sprintf( buffer, "Atom %u ", kk );
            if( matchMaps( buffer, pMap1, pMap2, amoebaGpu->log ) == 0 ){
                (void) fprintf( amoebaGpu->log, "%s ok\n", buffer );
            }
        }
        (void) fclose( filePtr );
#endif
    }

    // Mark all interactions that involve a padding atom as being excluded.

#if 0
    for (int atom1 = gpu->natoms; atom1 < (int)atoms; ++atom1)
    {
        int x = atom1/grid;
        int offset1 = atom1-x*grid;
        for (int atom2 = 0; atom2 < (int)atoms; ++atom2)
        {
            int y = atom2/grid;
            int offset2 = atom2-y*grid;
            if (x >= y)
            {
                int cell = x+y*dim-y*(y+1)/2;
                pScalingIndices[pScalingIndicesIndex[cell]+offset1] &= 0xFFFFFFFF-(1<<offset2);
            }
            if (y >= x)
            {
                int cell = y+x*dim-x*(x+1)/2;
                pScalingIndices[pScalingIndicesIndex[cell]+offset2] &= 0xFFFFFFFF-(1<<offset1);
            }
        }
    }
#endif

    amoebaGpu->ps_D_ScaleIndices->Upload();
    amoebaGpu->ps_P_ScaleIndices->Upload();
    amoebaGpu->ps_M_ScaleIndices->Upload();
    amoebaGpu->psScalingIndicesIndex->Upload();
    amoebaGpu->psWorkUnit->Upload();
}

/**---------------------------------------------------------------------------------------

   Get threads/block

   @param amoebaGpu        amoebaGpuContext
   @param sharedMemoryPerThread shared memory/thread

   @return threadsPerBlock

   --------------------------------------------------------------------------------------- */

unsigned int getThreadsPerBlock( amoebaGpuContext amoebaGpu, unsigned int sharedMemoryPerThread )
{
    unsigned int grid               = amoebaGpu->gpuContext->grid;
    unsigned int threadsPerBlock    = (amoebaGpu->gpuContext->sharedMemoryPerBlock + grid -1)/(grid*sharedMemoryPerThread);
    threadsPerBlock                 = threadsPerBlock < 1 ? 1 : threadsPerBlock;
    threadsPerBlock                 *= grid;

   return threadsPerBlock;
}

/**---------------------------------------------------------------------------------------

   Open file for writing -- centralized, utility routine

   @param fname            base file name
   @param step             timestep -- appended to base file name

   --------------------------------------------------------------------------------------- */

FILE* getWriteToFilePtr( char* fname, int step )
{
   std::stringstream fileName;
   fileName << fname;
   fileName << "_" << step;
   fileName << ".txt";

#ifdef WIN32
   FILE* filePtr;
   fopen_s( &filePtr, fileName.str().c_str(), "w" );
#else
   FILE* filePtr = fopen( fileName.str().c_str(), "w" );
#endif
   if( filePtr == NULL ){
      (void) fprintf( stderr, "Could not open file=<%s> for writitng.", fileName.str().c_str() );
      exit(-1);
   }
   return filePtr;
}

/**---------------------------------------------------------------------------------------

   Open file for writing -- centralized, utility routine

   @param fname            base file name
   @param step             timestep -- appended to base file name

   --------------------------------------------------------------------------------------- */

FILE* getWriteToFilePtrV( char* fname, std::vector<int>& fileId )
{
    std::stringstream fileName;
    fileName << fname;
 
    for( std::vector<int>::const_iterator ii = fileId.begin(); ii != fileId.end(); ii++ ){
       fileName << "_" << *ii;
    }
 
    fileName << ".txt";
 #ifdef WIN32
    FILE* filePtr;
    fopen_s( &filePtr, fileName.str().c_str(), "w" );
 #else
    FILE* filePtr = fopen( fileName.str().c_str(), "w" );
 #endif
    if( filePtr == NULL ){
       (void) fprintf( stderr, "Could not open file=<%s> for writing.", fileName.str().c_str() );
       exit(-1);
    }
    return filePtr;
}

/**---------------------------------------------------------------------------------------

   Print values in array -- utility routine -- centralized -- one line per call

   @param filePtr          file ptr to write values to
   @param index            particles/line index (?)
   @param numberOfValues   number of values in array
   @param values           array of values

   --------------------------------------------------------------------------------------- */

extern "C" {
static void printValues( FILE* filePtr, int index, int numberOfValues, float* values )
{
    int ii;
    (void) fprintf( filePtr, "%5d ", index );
    for ( ii = 0; ii < numberOfValues; ii++ ) { 
       (void) fprintf( filePtr, " %18.10e", values[ii] );
    }
    (void) fprintf( filePtr, "\n" );
} 
}

/**---------------------------------------------------------------------------------------

   Write contents of two arrays to file

   @param numberOfParticles    number of particles
   @param fname            base file name
   @param timestep         tomestep -- appended to base file name
   @param entriesPerParticle1  entries to be written for first array (usually 3 for xyz)
   @param array1           first array (typically coordinate array)
   @param entriesPerParticle2  entries to be written for second array (3 for dipoles, 9 for quadrupole, ...)
   @param array2           second array

   --------------------------------------------------------------------------------------- */

//extern "C"
void cudaWriteFloat4AndFloat1ArraysToFile( int numberOfParticles, char* fname, int timestep, int entriesPerParticle1, CUDAStream<float4>* array1, 
                                           int entriesPerParticle2, CUDAStream<float>* array2 )
{
    // ---------------------------------------------------------------------------------------

    // static const std::string methodName = "cudaWrite4And1ArraysToFile";

    // ---------------------------------------------------------------------------------------

    int ii, jj;
    FILE* filePtr;
    float values[20];
 
    array1->Download();
    array2->Download();
 
    filePtr          = getWriteToFilePtr( fname, timestep );
 
    int runningIndex = 0;
    int offset       = entriesPerParticle1 == 4 ? 4 : 3;
    for ( ii = 0; ii < numberOfParticles; ii++ ){ 
 
       values[0] = array1->_pSysStream[0][ii].x;
       values[1] = array1->_pSysStream[0][ii].y;
       values[2] = array1->_pSysStream[0][ii].z;
 
       for( jj = 0; jj < entriesPerParticle2; jj++ ) { 
          values[offset+jj] = array2->_pSysStream[0][runningIndex++];
       }
       printValues( filePtr, ii, (offset+entriesPerParticle2), values ); 
    }
    (void) fflush( filePtr );
    (void) fclose( filePtr );
}

/**---------------------------------------------------------------------------------------

   Write contents of two arrays to file

   @param numberOfParticles    number of particles
   @param fname            base file name
   @param timestep         timestep -- appended to base file name
   @param entriesPerParticle1  entries to be written for first array
   @param array1           first array
   @param entriesPerParticle2  entries to be written for second array
   @param array2           second array

   --------------------------------------------------------------------------------------- */

void cudaWriteFloat1AndFloat1ArraysToFile( int numberOfParticles, char* fname, std::vector<int>& fileId, int entriesPerParticle1, CUDAStream<float>* array1, 
                                           int entriesPerParticle2, CUDAStream<float>* array2 )
{
    // ---------------------------------------------------------------------------------------

    // static const std::string methodName = "cudaWrite1And1ArraysToFile";

    // ---------------------------------------------------------------------------------------

    int ii, jj;
    FILE* filePtr;
    float values[50];
 
    array1->Download();
    if( entriesPerParticle2 > 0 && array2 ){
        array2->Download();
    }
 
    filePtr            = getWriteToFilePtrV( fname, fileId );
 
    int runningIndex1  = 0;
    int runningIndex2  = 0;
 
    float sum1         = 0.0f;
    float sum2         = 0.0f;
 
    for ( ii = 0; ii < numberOfParticles; ii++ ){ 
       for( jj = 0; jj < entriesPerParticle1; jj++ ) { 
          sum1 += array1->_pSysStream[0][runningIndex1]*array1->_pSysStream[0][runningIndex1];
          runningIndex1++;
       }
       for( jj = 0; jj < entriesPerParticle2; jj++ ) { 
          sum2 += array2->_pSysStream[0][runningIndex2]*array2->_pSysStream[0][runningIndex2];
          runningIndex2++;
       }
    }
    (void) fprintf( filePtr, "%d   # norm: %.6e %.6e\n", numberOfParticles, sqrtf( sum1 ), sqrtf( sum2 ) );
 
    runningIndex1  = 0;
    runningIndex2  = 0;
    for ( ii = 0; ii < numberOfParticles; ii++ ){ 
  
       int index = 0;
       for( jj = 0; jj < entriesPerParticle1; jj++ ) { 
          values[index++] = array1->_pSysStream[0][runningIndex1++];
       }
       for( jj = 0; jj < entriesPerParticle2; jj++ ) { 
          values[index++] = array2->_pSysStream[0][runningIndex2++];
       }
       printValues( filePtr, ii, index, values ); 
    }
 
    (void) fflush( filePtr );
    (void) fclose( filePtr );
}

/**---------------------------------------------------------------------------------------

   Write contents of two arrays to file

   @param numberOfParticles    number of particles
   @param fname            base file name
   @param timestep         timestep -- appended to base file name
   @param entriesPerParticle1  entries to be written for first array
   @param array1           first array
   @param entriesPerParticle2  entries to be written for second array
   @param array2           second array

   --------------------------------------------------------------------------------------- */

void cudaWriteVectorOfDoubleVectorsToFile( char* fname, std::vector<int>& fileId,
                                           VectorOfDoubleVectors& outputVector )
{
    // ---------------------------------------------------------------------------------------

    // static const std::string methodName = "cudaWrite1And1ArraysToFile";

    // ---------------------------------------------------------------------------------------

    FILE* filePtr            = getWriteToFilePtrV( fname, fileId );
    (void) fprintf( filePtr, "%u\n", static_cast<unsigned int>(outputVector.size()) );
 
    float values[50];
    for ( unsigned int ii = 0; ii < outputVector.size(); ii++ ){ 
        int index = 0;
        for ( unsigned int jj = 0; jj < outputVector[ii].size(); jj++ ){ 
            values[index++] = static_cast<float>(outputVector[ii][jj]);
        }
        printValues( filePtr, static_cast<int>(ii), index, values ); 
    }
    (void) fflush( filePtr );
    (void) fclose( filePtr );
}

/**---------------------------------------------------------------------------------------

   Load contents of arrays into vector

   @param numberOfParticles    number of particles
   @param entriesPerParticle   entries/particles array
   @param array                cuda array
   @param outputVector         output vector

   --------------------------------------------------------------------------------------- */

void cudaLoadCudaFloatArray( int numberOfParticles, int entriesPerParticle, CUDAStream<float>* array, VectorOfDoubleVectors& outputVector )
{
    // ---------------------------------------------------------------------------------------

    // static const std::string methodName = "cudaLoadCudaFloatArray";

    // ---------------------------------------------------------------------------------------

    array->Download();
    int runningIndex  = 0;
    
    outputVector.resize( numberOfParticles ); 
 
    for( int ii = 0; ii < numberOfParticles; ii++ ){ 
       for( int jj = 0; jj < entriesPerParticle; jj++ ) { 
          outputVector[ii].push_back( array->_pSysStream[0][runningIndex++] );
       }
    }
}

/**---------------------------------------------------------------------------------------

   Load contents of arrays into vector

   @param numberOfParticles    number of particles
   @param entriesPerParticle   entries/particles array
   @param array                cuda array
   @param outputVector         output vector

   --------------------------------------------------------------------------------------- */

void cudaLoadCudaFloat2Array( int numberOfParticles, int entriesPerParticle, CUDAStream<float2>* array, VectorOfDoubleVectors& outputVector ) 
{
    // ---------------------------------------------------------------------------------------

    // static const std::string methodName = "cudaLoadCudaFloat2Array";

    // ---------------------------------------------------------------------------------------

    array->Download();
    int runningIndex  = 0;
    
    outputVector.resize( numberOfParticles ); 
 
    for( int ii = 0; ii < numberOfParticles; ii++ ){ 
        if( entriesPerParticle > 0 ){
            outputVector[ii].push_back( array->_pSysStream[0][runningIndex].x );
        }
        if( entriesPerParticle > 1 ){
            outputVector[ii].push_back( array->_pSysStream[0][runningIndex].y );
        }
        runningIndex++;
    }
}

/**---------------------------------------------------------------------------------------

   Load contents of arrays into vector

   @param numberOfParticles    number of particles
   @param entriesPerParticle   entries/particles array
   @param array                cuda array
   @param outputVector         output vector

   --------------------------------------------------------------------------------------- */

void cudaLoadCudaFloat4Array( int numberOfParticles, int entriesPerParticle, CUDAStream<float4>* array, VectorOfDoubleVectors& outputVector ) 
{
    // ---------------------------------------------------------------------------------------

    // static const std::string methodName = "cudaLoadCudaFloat4Array";

    // ---------------------------------------------------------------------------------------

    array->Download();
    int runningIndex  = 0;
    
    outputVector.resize( numberOfParticles ); 
 
    for( int ii = 0; ii < numberOfParticles; ii++ ){ 
        if( entriesPerParticle > 0 ){
            outputVector[ii].push_back( array->_pSysStream[0][runningIndex].x );
        }
        if( entriesPerParticle > 1 ){
            outputVector[ii].push_back( array->_pSysStream[0][runningIndex].y );
        }
        if( entriesPerParticle > 2 ){
            outputVector[ii].push_back( array->_pSysStream[0][runningIndex].z );
        }
        if( entriesPerParticle > 3 ){
            outputVector[ii].push_back( array->_pSysStream[0][runningIndex].w );
        }
        runningIndex++;
    }
}

/**---------------------------------------------------------------------------------------

   Get norm squared of vector

   @param numberOfElements number of elements
   @param array            array

   @return norm**2

   --------------------------------------------------------------------------------------- */

float cudaGetSum( int numberOfElements, CUDAStream<float>* array )
{
    // ---------------------------------------------------------------------------------------

    // static const std::string methodName = "cudaGetNorm2";

    // ---------------------------------------------------------------------------------------

    int ii;
    float sum;
 
    array->Download();
 
    sum = 0.0f;
    for ( ii = 0; ii < numberOfElements; ii++ ){ 
       sum += array->_pSysStream[0][ii];
    }
    return sum;
}

/**---------------------------------------------------------------------------------------

   Get norm squared of vector

   @param numberOfElements number of elements
   @param array            array

   @return norm**2

   --------------------------------------------------------------------------------------- */

float cudaGetNorm2( int numberOfElements, CUDAStream<float>* array )
{
    // ---------------------------------------------------------------------------------------

    // static const std::string methodName = "cudaGetNorm2";

    // ---------------------------------------------------------------------------------------

    int ii;
    float sum;
 
    array->Download();
 
    sum = 0.0f;
    for ( ii = 0; ii < numberOfElements; ii++ ){ 
       sum += (array->_pSysStream[0][ii]*array->_pSysStream[0][ii]);
    }
    return sum;
}

/**---------------------------------------------------------------------------------------

   Return count of nans/infinities in array

   @param numberOfElements number of elements
   @param array            array

   @return  count of nans/infinities 

   --------------------------------------------------------------------------------------- */

int checkForNansAndInfinities( int numberOfElements, CUDAStream<float>* array )
{
    // ---------------------------------------------------------------------------------------

    // static const std::string methodName = "cudaGetNorm2";

    // ---------------------------------------------------------------------------------------

    array->Download();
    int nansDetected = 0;
    for( int ii = 0; ii < numberOfElements; ii++ ){
        if( array->_pSysStream[0][ii] !=  array->_pSysStream[0][ii]   ||  
            array->_pSysStream[0][ii] ==  std::numeric_limits<double>::infinity() ||
            array->_pSysStream[0][ii] == -std::numeric_limits<double>::infinity() ){
            nansDetected++;
        }
    }   
 
    return nansDetected;
}

/**---------------------------------------------------------------------------------------

   Replacement of sorts for strtok()
   Used to parse parameter file lines

   @param lineBuffer           string to tokenize
   @param delimiter            token delimter

   @return number of args

   --------------------------------------------------------------------------------------- */

static char* strsepLocal( char** lineBuffer, const char* delimiter ){

   // ---------------------------------------------------------------------------------------

   // static const std::string methodName = "strsepLocal";

   char *s;
   const char *spanp;
   int c, sc; 
   char *tok;
   
   // ---------------------------------------------------------------------------------------

   s = *lineBuffer;
   if( s == NULL ){
      return (NULL);
   }
   
   for( tok = s;; ){
      c     = *s++;
      spanp = delimiter;
      do {
         if( (sc = *spanp++) == c ){
            if( c == 0 ){
               s = NULL;
            } else {
               s[-1] = 0;
            }
/*
            if( *s == '\n' ){ 
               *s = NULL;
            }
*/
            *lineBuffer = s;
            return( tok );
         }
      } while( sc != 0 );
   }
}  

/**---------------------------------------------------------------------------------------

   Tokenize a string

   @param lineBuffer           string to tokenize
   @param tokenArray           upon return vector of tokens
   @param delimiter            token delimter

   @return number of tokens

   --------------------------------------------------------------------------------------- */

static int tokenizeString( char* lineBuffer, StringVector& tokenArray, const std::string delimiter ){

   // ---------------------------------------------------------------------------------------

   // static const std::string methodName = "\nSimTKOpenMMUtilities::tokenizeString";

   // ---------------------------------------------------------------------------------------

   char *ptr_c = NULL;

   for( ; (ptr_c = strsepLocal( &lineBuffer, delimiter.c_str() )) != NULL; ){
      if( *ptr_c ){
/*
         char* endOfLine = ptr_c;
         while( endOfLine ){
printf( "%c", *endOfLine ); fflush( stdout );
            if( *endOfLine == '\n' )*endOfLine = '\0';
            endOfLine++;
         }  
*/
         tokenArray.push_back( std::string( ptr_c ) );
      }
   }

   return (int) tokenArray.size();
}

/**---------------------------------------------------------------------------------------

   Read a line from a file and tokenize into an array of strings

   @param filePtr              file to read from
   @param tokens               array of token strings
   @param lineCount            line count
   @param log                  optional file ptr for logging

   @return ptr to string containing line

   --------------------------------------------------------------------------------------- */

static char* readLineFromFile( FILE* filePtr, StringVector& tokens ){

// ---------------------------------------------------------------------------------------

   //static const std::string methodName      = "readLine";

   std::string delimiter                    = " \n";
   const int bufferSize                     = 4096;
   char buffer[bufferSize];

// ---------------------------------------------------------------------------------------

   char* isNotEof = fgets( buffer, bufferSize, filePtr );
   if( isNotEof ){
      tokenizeString( buffer, tokens, delimiter );
   }
   return isNotEof;

}

/**---------------------------------------------------------------------------------------

   Read a file

   @param fileName             file name
   @param fileContents         output file contents

   --------------------------------------------------------------------------------------- */

void readFile( std::string fileName, StringVectorVector& fileContents ){

// ---------------------------------------------------------------------------------------

   //static const std::string methodName      = "readFile";

// ---------------------------------------------------------------------------------------

    fileContents.resize(0);
    FILE* filePtr = fopen( fileName.c_str(), "r" );
    StringVector firstLine;
    char* isNotEof = readLineFromFile( filePtr, firstLine);
    fileContents.push_back( firstLine );
    //int lineCount  = 0;
    while( isNotEof ){
        StringVector lineTokens;
        isNotEof = readLineFromFile( filePtr, lineTokens );
        fileContents.push_back( lineTokens );
    }
    (void) fclose( filePtr );

    return;
}

/**---------------------------------------------------------------------------------------

   Report whether a number is a nan or infinity

   @param number               number to test
   @return 1 if number is  nan or infinity; else return 0

   --------------------------------------------------------------------------------------- */

int isNanOrInfinity( double number ){
    return (number != number || number == std::numeric_limits<double>::infinity() || number == -std::numeric_limits<double>::infinity()) ? 1 : 0; 
}

/**---------------------------------------------------------------------------------------

   Track iterations for MI dipoles

   @param amoebaGpu            amoebaGpuContext reference
   @param iteration            MI iteration

   --------------------------------------------------------------------------------------- */

void trackMutualInducedIterations( amoebaGpuContext amoebaGpu, int iteration){

// ---------------------------------------------------------------------------------------

    static const std::string methodName      = "trackMutualInducedIterations";
    static int currentStep                   = 0; 
    static double iterationStat[6]           = { 0.0, 0.0, 1000.0, 0.0, 0.0, 0.0 };

// ---------------------------------------------------------------------------------------

    if( amoebaGpu->log == NULL || currentStep > 20000 )return;
    //if( amoebaGpu->log == NULL )return;

    gpuContext gpu                       = amoebaGpu->gpuContext;
    currentStep++;

    double interationD = static_cast<double>(iteration);
    iterationStat[0]  += interationD;
    iterationStat[1]  += interationD*interationD;
    iterationStat[2]   = interationD < iterationStat[2] ? interationD : iterationStat[2];
    iterationStat[3]   = interationD > iterationStat[3] ? interationD : iterationStat[3];
    iterationStat[4]  += 1.0; 
    if( iterationStat[4] >= 1000.0 ){
        double average = iterationStat[0]/iterationStat[4]; 
        double stddev  = iterationStat[1] - average*average*iterationStat[4]; 
               stddev  = sqrt( stddev )/(iterationStat[4]-1.0);
        (void) fprintf( amoebaGpu->log, "%s %8d iteration=%10.3f stddev=%10.3f min/max[%10.3f %10.3f] %10.1f eps=%14.7e\n",
                        methodName.c_str(), currentStep, average, stddev, iterationStat[2], iterationStat[3], iterationStat[4], amoebaGpu->mutualInducedCurrentEpsilon );
        (void) fflush( amoebaGpu->log );
        iterationStat[0] = iterationStat[1] = iterationStat[4] = 0.0; 
    }
    if( 0 ){ 
        std::vector<int> fileId;
        if( interationD < (amoebaGpu->mutualInducedMaxIterations-10) ) {
            int id = (currentStep % 20); 
            fileId.push_back( id );
        } else {
            fileId.push_back( currentStep );
        }    
        if( (currentStep % 20) == 0 || fileId[0] > 20 ){
             (void) fprintf( amoebaGpu->log, "step=%d fileId=%d iterations=%d\n", currentStep, fileId[0], iteration );
        }
        (void) fflush( amoebaGpu->log );
        VectorOfDoubleVectors outputVector;
        cudaLoadCudaFloat4Array( gpu->natoms, 3, gpu->psPosq4,                    outputVector );
        cudaLoadCudaFloat4Array( gpu->natoms, 3, gpu->psVelm4,                    outputVector );
/*
            cudaLoadCudaFloatArray( gpu->natoms,  3, amoebaGpu->psInducedDipole,      outputVector );
            cudaLoadCudaFloatArray( gpu->natoms,  3, amoebaGpu->psInducedDipolePolar, outputVector );
            cudaLoadCudaFloatArray( gpu->natoms,  3, amoebaGpu->psInducedDipoleS,     outputVector );
            cudaLoadCudaFloatArray( gpu->natoms,  3, amoebaGpu->psInducedDipolePolarS,outputVector );
*/
        cudaWriteVectorOfDoubleVectorsToFile( "CudaMI", fileId, outputVector );
        int nansPresent = isNanOrInfinity( amoebaGpu->mutualInducedCurrentEpsilon );
        if( nansPresent == 0 ){ 
            for( int ii = 0; ii < gpu->natoms && nansPresent == 0; ii++ ){
                if( isNanOrInfinity( gpu->psPosq4->_pSysStream[0][ii].x ) || 
                    isNanOrInfinity( gpu->psPosq4->_pSysStream[0][ii].y ) || 
                    isNanOrInfinity( gpu->psPosq4->_pSysStream[0][ii].z ) || 
                    isNanOrInfinity( gpu->psVelm4->_pSysStream[0][ii].x  ) || 
                    isNanOrInfinity( gpu->psVelm4->_pSysStream[0][ii].y  ) || 
                    isNanOrInfinity( gpu->psVelm4->_pSysStream[0][ii].z  ) ){ 
                    nansPresent = 1; 
                }
            }
        }
        if( nansPresent ){
            (void) fprintf( amoebaGpu->log, "epsilon nan - exiting\n" );
            (void) fflush( amoebaGpu->log );
            exit(-1);
        }
    }
}

#undef  AMOEBA_DEBUG
