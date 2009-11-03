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
#include "cudaKernels.h"

struct Atom {
    float x;
    float y;
    float z;
    float r;
    float sr;
    float fx;
    float fy;
    float fz;
    float fb;
};


static __constant__ cudaGmxSimulation cSim;

void SetCalculateGBVIForces2Sim(gpuContext gpu)
{
    cudaError_t status;
    status = cudaMemcpyToSymbol(cSim, &gpu->sim, sizeof(cudaGmxSimulation));     
    RTERROR(status, "cudaMemcpyToSymbol: SetSim copy to cSim failed");
}

void GetCalculateGBVIForces2Sim(gpuContext gpu)
{
    cudaError_t status;
    status = cudaMemcpyFromSymbol(&gpu->sim, cSim, sizeof(cudaGmxSimulation));     
    RTERROR(status, "cudaMemcpyFromSymbol: SetSim copy from cSim failed");
}

#include "kCalculateGBVIAux.h"

/**
 * This file contains the kernel for evalauating the second stage of GBSA.  It is included
 * several times in kCalculateGBVIForces2.cu with different #defines to generate
 * different versions of the kernels.
 */

__global__ void kCalculateGBVIForces2a_kernel()
{
    unsigned int pos  = (blockIdx.x * blockDim.x + threadIdx.x);
    if( pos >=  cSim.atoms )return;

    float4 apos                     = cSim.pPosq[pos];
    float4 ar                       = cSim.pGBVIData[pos];
    float fb                        = cSim.pBornForce[pos];
    unsigned int posJ               = 0;
    float4 force;
    force.x = force.y = force.z = force.w = 0.0f;
    while ( posJ < cSim.atoms )
    {

        float4 aposJ                = cSim.pPosq[posJ];
        float4 arJ                  = cSim.pGBVIData[posJ];
        float fbJ                   = cSim.pBornForce[posJ];

        float dx                    = aposJ.x - apos.x;
        float dy                    = aposJ.y - apos.y;
        float dz                    = aposJ.z - apos.z;

        float r2                    = dx * dx + dy * dy + dz * dz;
        float r                     = sqrt(r2);

        float dE                    = getGBVI_dE2( r, ar.x, arJ.y, fb );
        dE                          = r > 1.0e-08f ? dE : 0.0f;

//dx = dy = dz = 1.0f;
        float d                     = dx*dE;
        force.x                    -= d;
        d                           = dy*dE;
        force.y                    -= d;
        d                           = dz*dE;
        force.z                    -= d;
#if 1
        dE                          = getGBVI_dE2( r, arJ.x, ar.y, fbJ );
        dE                          = r > 1.0e-08f ? dE : 0.0f;
        d                           = dx*dE;
        force.x                    -= d;
        d                           = dy*dE;
        force.y                    -= d;
        d                           = dz*dE;
        force.z                    -= d;
#endif

        posJ                       += 1;
    }

    // Write results
    cSim.pForce4a[pos]             = force;

}

// Include versions of the kernels for N^2 calculations.

#define METHOD_NAME(a, b) a##N2##b
#include "kCalculateGBVIForces2.h"
#define USE_OUTPUT_BUFFER_PER_WARP
#undef METHOD_NAME
#define METHOD_NAME(a, b) a##N2ByWarp##b
#include "kCalculateGBVIForces2.h"

// Include versions of the kernels with cutoffs.

#undef METHOD_NAME
#undef USE_OUTPUT_BUFFER_PER_WARP
#define USE_CUTOFF
#define METHOD_NAME(a, b) a##Cutoff##b
#include "kCalculateGBVIForces2.h"
#define USE_OUTPUT_BUFFER_PER_WARP
#undef METHOD_NAME
#define METHOD_NAME(a, b) a##CutoffByWarp##b
#include "kCalculateGBVIForces2.h"

// Include versions of the kernels with periodic boundary conditions.

#undef METHOD_NAME
#undef USE_OUTPUT_BUFFER_PER_WARP
#define USE_PERIODIC
#define METHOD_NAME(a, b) a##Periodic##b
#include "kCalculateGBVIForces2.h"
#define USE_OUTPUT_BUFFER_PER_WARP
#undef METHOD_NAME
#define METHOD_NAME(a, b) a##PeriodicByWarp##b
#include "kCalculateGBVIForces2.h"

void kCalculateGBVIForces2(gpuContext gpu)
{
    //printf("kCalculateGBVIForces2\n");
    size_t numWithInteractions;

#if 0
    kClearForces(gpu);
    (void) fprintf( stderr, "\nkCalculateGBVIForces2: cleared force prior loop2\n" ); (void) fflush( stderr );
    kCalculateGBVIForces2a_kernel<<<gpu->sim.blocks, 384>>>();
    (void) fprintf( stderr, "\ncalled kCalculateGBVIForces2a\n" ); (void) fflush( stderr );
    return;
#endif

    switch (gpu->sim.nonbondedMethod)
    {
        case NO_CUTOFF:
            if (gpu->bOutputBufferPerWarp)
                kCalculateGBVIN2ByWarpForces2_kernel<<<gpu->sim.bornForce2_blocks, gpu->sim.bornForce2_threads_per_block,
                        sizeof(Atom)*gpu->sim.bornForce2_threads_per_block>>>(gpu->sim.pWorkUnit, gpu->sim.workUnits);
            else
                kCalculateGBVIN2Forces2_kernel<<<gpu->sim.bornForce2_blocks, gpu->sim.bornForce2_threads_per_block,
                        sizeof(Atom)*gpu->sim.bornForce2_threads_per_block>>>(gpu->sim.pWorkUnit, gpu->sim.workUnits);
//(void) fprintf( stderr, "\nkCalculateGBVIForces2: Born radii/force forces warp=%u\n", gpu->bOutputBufferPerWarp ); (void) fflush( stderr );
#define GBVI_DEBUG 0
#if ( GBVI_DEBUG == 1 )
                (void) fprintf( stderr, "\nkCalculateGBVIForces2: Born radii/force forces:\n" ); (void) fflush( stderr );
                gpu->psBornForce->Download();
                gpu->psForce4->Download();
                for( int ii = 0; ii < gpu->natoms; ii++ ){
                    (void) fprintf( stderr, "%d bF=%14.6e Fa[%14.6e %14.6e %14.6e] Fb[%14.6e %14.6e %14.6e]\n",
                                    ii,
                                    gpu->psBornForce->_pSysStream[0][ii],
                                    gpu->psForce4->_pSysStream[0][ii].x,
                                    gpu->psForce4->_pSysStream[0][ii].y,
                                    gpu->psForce4->_pSysStream[0][ii].z,
                                    gpu->psForce4->_pSysStream[1][ii].x,
                                    gpu->psForce4->_pSysStream[1][ii].y,
                                    gpu->psForce4->_pSysStream[1][ii].z
                                  );  
                }   
                for( int ii = 0; ii < gpu->sim.paddedNumberOfAtoms*2; ii++ ){
                    (void) fprintf( stderr, "%d bF=%14.6e Fa[%14.6e %14.6e %14.6e %14.6e]\n",
                                    ii,
                                    gpu->psBornForce->_pSysStream[0][ii],
                                    gpu->psForce4->_pSysStream[0][ii].x,
                                    gpu->psForce4->_pSysStream[0][ii].y,
                                    gpu->psForce4->_pSysStream[0][ii].z,
                                    gpu->psForce4->_pSysStream[0][ii].w
                                  );  
                }   
#endif
#undef GBVI_DEBUG

            break;
        case CUTOFF:
            numWithInteractions = gpu->psInteractionCount->_pSysData[0];
            if (gpu->bOutputBufferPerWarp)
                kCalculateGBVICutoffByWarpForces2_kernel<<<gpu->sim.bornForce2_blocks, gpu->sim.bornForce2_threads_per_block,
                        (sizeof(Atom)+sizeof(float3))*gpu->sim.bornForce2_threads_per_block>>>(gpu->sim.pInteractingWorkUnit, numWithInteractions);
            else
                kCalculateGBVICutoffForces2_kernel<<<gpu->sim.bornForce2_blocks, gpu->sim.bornForce2_threads_per_block,
                        (sizeof(Atom)+sizeof(float3))*gpu->sim.bornForce2_threads_per_block>>>(gpu->sim.pInteractingWorkUnit, numWithInteractions);
            break;
        case PERIODIC:
            numWithInteractions = gpu->psInteractionCount->_pSysData[0];
            if (gpu->bOutputBufferPerWarp)
                kCalculateGBVIPeriodicByWarpForces2_kernel<<<gpu->sim.bornForce2_blocks, gpu->sim.bornForce2_threads_per_block,
                        (sizeof(Atom)+sizeof(float3))*gpu->sim.bornForce2_threads_per_block>>>(gpu->sim.pInteractingWorkUnit, numWithInteractions);
            else
                kCalculateGBVIPeriodicForces2_kernel<<<gpu->sim.bornForce2_blocks, gpu->sim.bornForce2_threads_per_block,
                        (sizeof(Atom)+sizeof(float3))*gpu->sim.bornForce2_threads_per_block>>>(gpu->sim.pInteractingWorkUnit, numWithInteractions);
            break;
    }
    LAUNCHERROR("kCalculateGBVIForces2");
}
