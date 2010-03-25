/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2009 Stanford University and the Authors.           *
 * Authors: Scott Le Grand, Peter Eastman                                     *
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

#include <stdio.h>
#include <cuda.h>
#include <vector_functions.h>
#include <cstdlib>
#include <string>
#include <iostream>
#include <fstream>
using namespace std;

#include "gputypes.h"

#define UNROLLXX 0
#define UNROLLXY 0

struct Atom {
    float x;
    float y;
    float z;
    float r;
    float sr;
    float sum;
    float gamma;
};

static __constant__ cudaGmxSimulation cSim;

void SetCalculateGBVIBornSumSim(gpuContext gpu)
{
    cudaError_t status;
    status = cudaMemcpyToSymbol(cSim, &gpu->sim, sizeof(cudaGmxSimulation));     
    RTERROR(status, "cudaMemcpyToSymbol: SetSim copy to cSim failed");
}

void GetCalculateGBVIBornSumSim(gpuContext gpu)
{
    cudaError_t status;
    status = cudaMemcpyFromSymbol(&gpu->sim, cSim, sizeof(cudaGmxSimulation));     
    RTERROR(status, "cudaMemcpyFromSymbol: SetSim copy from cSim failed");
}

// Include versions of the kernels for N^2 calculations.

#define METHOD_NAME(a, b) a##N2##b
#include "kCalculateGBVIBornSum.h"
#define USE_OUTPUT_BUFFER_PER_WARP
#undef METHOD_NAME
#define METHOD_NAME(a, b) a##N2ByWarp##b
#include "kCalculateGBVIBornSum.h"

// Include versions of the kernels with cutoffs.

#undef METHOD_NAME
#undef USE_OUTPUT_BUFFER_PER_WARP
#define USE_CUTOFF
#define METHOD_NAME(a, b) a##Cutoff##b
#include "kCalculateGBVIBornSum.h"
#define USE_OUTPUT_BUFFER_PER_WARP
#undef METHOD_NAME
#define METHOD_NAME(a, b) a##CutoffByWarp##b
#include "kCalculateGBVIBornSum.h"

// Include versions of the kernels with periodic boundary conditions.

#undef METHOD_NAME
#undef USE_OUTPUT_BUFFER_PER_WARP
#define USE_PERIODIC
#define METHOD_NAME(a, b) a##Periodic##b
#include "kCalculateGBVIBornSum.h"
#define USE_OUTPUT_BUFFER_PER_WARP
#undef METHOD_NAME
#define METHOD_NAME(a, b) a##PeriodicByWarp##b
#include "kCalculateGBVIBornSum.h"

__global__ void kReduceGBVIBornSum_kernel()
{
    unsigned int pos = (blockIdx.x * blockDim.x + threadIdx.x);
    
    while (pos < cSim.atoms)
    {
        float sum = 0.0f;
        float* pSt = cSim.pBornSum + pos;
        float4 atom = cSim.pGBVIData[pos];
        
        // Get summed Born data
        for (int i = 0; i < cSim.nonbondOutputBuffers; i++)
        {
            sum += *pSt;
       //     printf("%4d %4d A: %9.4f\n", pos, i, *pSt);
            pSt += cSim.stride;
        }
        
        // Now calculate Born radius

        float Rinv           = 1.0f/atom.x;
        sum                  = Rinv*Rinv*Rinv - sum; 
        cSim.pBornRadii[pos] = pow( sum, (-1.0f/3.0f) ); 
        pos += gridDim.x * blockDim.x;
    }   
}

void kReduceGBVIBornSum(gpuContext gpu)
{
    //printf("kReduceGBVIBornSum\n");
#define GBVI_DEBUG 0
#if ( GBVI_DEBUG == 1 )
               gpu->psGBVIData->Download();
               gpu->psBornSum->Download();
               gpu->psPosq4->Download();
                (void) fprintf( stderr, "\nkReduceGBVIBornSum: Post BornSum %s Born radii & params\n", 
                               (gpu->bIncludeGBVI ? "GBVI" : "Obc") );
                for( int ii = 0; ii < gpu->natoms; ii++ ){
                   (void) fprintf( stderr, "%d bSum=%14.6e param[%14.6e %14.6e %14.6e] x[%14.6f %14.6f %14.6f %14.6f]\n",
                                   ii, 
                                   gpu->psBornSum->_pSysStream[0][ii],
                                   gpu->psGBVIData->_pSysStream[0][ii].x,
                                   gpu->psGBVIData->_pSysStream[0][ii].y,
                                   gpu->psGBVIData->_pSysStream[0][ii].z,
                                   gpu->psPosq4->_pSysStream[0][ii].x, gpu->psPosq4->_pSysStream[0][ii].y,
                                   gpu->psPosq4->_pSysStream[0][ii].z, gpu->psPosq4->_pSysStream[0][ii].w
                                 );  
                }   
#endif
#undef GBVI_DEBUG


    kReduceGBVIBornSum_kernel<<<gpu->sim.blocks, 384>>>();
    gpu->bRecalculateBornRadii = false;
    LAUNCHERROR("kReduceGBVIBornSum");
}

void kCalculateGBVIBornSum(gpuContext gpu)
{
    //printf("kCalculateGBVIBornSum\n");
    //size_t numWithInteractions;
    switch (gpu->sim.nonbondedMethod)
    {
        case NO_CUTOFF:
#define GBVI 0
#if GBVI == 1
int maxPrint = 10;
gpu->psWorkUnit->Download();
fprintf( stderr, "kCalculateGBVIBornSum: bOutputBufferPerWarp=%u blks=%u th/blk=%u wu=%u %u shrd=%u\n", gpu->bOutputBufferPerWarp,
                 gpu->sim.nonbond_blocks, gpu->sim.nonbond_threads_per_block, gpu->sim.workUnits, gpu->psWorkUnit->_pSysStream[0][0],
        sizeof(Atom)*gpu->sim.nonbond_threads_per_block );

               gpu->psGBVIData->Download();
               gpu->psBornSum->Download();
               gpu->psPosq4->Download();

                (void) fprintf( stderr, "\nkCalculateGBVIBornSum: pre BornSum %s Born radii & params\n",
                               (gpu->bIncludeGBVI ? "GBVI" : "Obc") );
                for( int ii = 0; ii < gpu->natoms; ii++ ){
                   (void) fprintf( stderr, "%d bSum=%14.6e param[%14.6e %14.6e %14.6e] x[%14.6f %14.6f %14.6f %14.6f]\n",
                                   ii, 
                                   gpu->psBornSum->_pSysStream[0][ii],
                                   gpu->psGBVIData->_pSysStream[0][ii].x,
                                   gpu->psGBVIData->_pSysStream[0][ii].y,
                                   gpu->psGBVIData->_pSysStream[0][ii].z,
                                   gpu->psPosq4->_pSysStream[0][ii].x, gpu->psPosq4->_pSysStream[0][ii].y,
                                   gpu->psPosq4->_pSysStream[0][ii].z, gpu->psPosq4->_pSysStream[0][ii].w
                                 );
                   if( (ii == maxPrint) && ( ii < (gpu->natoms - maxPrint)) ){
                      ii = gpu->natoms - maxPrint;
                   }
                }

#endif
#undef GBVI
            if (gpu->bOutputBufferPerWarp){
                kCalculateGBVIN2ByWarpBornSum_kernel<<<gpu->sim.nonbond_blocks, gpu->sim.nonbond_threads_per_block,
                        sizeof(Atom)*gpu->sim.nonbond_threads_per_block>>>(gpu->sim.pWorkUnit);
            } else {
                kCalculateGBVIN2BornSum_kernel<<<gpu->sim.nonbond_blocks, gpu->sim.nonbond_threads_per_block,
                        sizeof(Atom)*gpu->sim.nonbond_threads_per_block>>>(gpu->sim.pWorkUnit);
            }
            break;
        case CUTOFF:
            if (gpu->bOutputBufferPerWarp)
                kCalculateGBVICutoffByWarpBornSum_kernel<<<gpu->sim.nonbond_blocks, gpu->sim.nonbond_threads_per_block,
                        (sizeof(Atom)+sizeof(float))*gpu->sim.nonbond_threads_per_block>>>(gpu->sim.pInteractingWorkUnit);
            else
                kCalculateGBVICutoffBornSum_kernel<<<gpu->sim.nonbond_blocks, gpu->sim.nonbond_threads_per_block,
                        (sizeof(Atom)+sizeof(float))*gpu->sim.nonbond_threads_per_block>>>(gpu->sim.pInteractingWorkUnit );
            break;
        case PERIODIC:
            if (gpu->bOutputBufferPerWarp)
                kCalculateGBVIPeriodicByWarpBornSum_kernel<<<gpu->sim.nonbond_blocks, gpu->sim.nonbond_threads_per_block,
                        (sizeof(Atom)+sizeof(float))*gpu->sim.nonbond_threads_per_block>>>(gpu->sim.pInteractingWorkUnit );
            else
                kCalculateGBVIPeriodicBornSum_kernel<<<gpu->sim.nonbond_blocks, gpu->sim.nonbond_threads_per_block,
                        (sizeof(Atom)+sizeof(float))*gpu->sim.nonbond_threads_per_block>>>(gpu->sim.pInteractingWorkUnit );
            break;
    }
    LAUNCHERROR("kCalculateGBVIBornSum");
}
