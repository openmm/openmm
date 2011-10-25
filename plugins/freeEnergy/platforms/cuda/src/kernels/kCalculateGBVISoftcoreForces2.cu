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

#include "gputypes.h"
#include "GpuFreeEnergyCudaKernels.h"
#include "openmm/OpenMMException.h"

#include <stdio.h>
#include <cuda.h>
#include <vector_functions.h>
#include <cstdlib>
#include <string>
#include <iostream>
#include <fstream>

using namespace std;

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
    float bornRadiusScaleFactor;
};


static __constant__ cudaGmxSimulation cSim;

void SetCalculateGBVISoftcoreForces2Sim( freeEnergyGpuContext gpu)
{
    cudaError_t status;
    status = cudaMemcpyToSymbol(cSim, &gpu->gpuContext->sim, sizeof(cudaGmxSimulation));     
    RTERROR(status, "cudaMemcpyToSymbol: SetCalculateGBVISoftcoreForces2Sim copy to cSim failed");
}

#include "kCalculateGBVIAux.h"

/**
 * This file contains the kernel for evalauating the second stage of GBSA.  It is included
 * several times in kCalculateGBVIForces2.cu with different #defines to generate
 * different versions of the kernels.
 */

__global__ void kCalculateGBVISoftcoreForces2a_kernel()
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
    cSim.pForce4[pos]              = force;

}

// Include versions of the kernels for N^2 calculations.

#define METHOD_NAME(a, b) a##N2##b
#include "kCalculateGBVISoftcoreForces2.h"
#define USE_OUTPUT_BUFFER_PER_WARP
#undef METHOD_NAME
#define METHOD_NAME(a, b) a##N2ByWarp##b
#include "kCalculateGBVISoftcoreForces2.h"

// Include versions of the kernels with cutoffs.

#undef METHOD_NAME
#undef USE_OUTPUT_BUFFER_PER_WARP
#define USE_CUTOFF
#define METHOD_NAME(a, b) a##Cutoff##b
#include "kCalculateGBVISoftcoreForces2.h"
#define USE_OUTPUT_BUFFER_PER_WARP
#undef METHOD_NAME
#define METHOD_NAME(a, b) a##CutoffByWarp##b
#include "kCalculateGBVISoftcoreForces2.h"

// Include versions of the kernels with periodic boundary conditions.

#undef METHOD_NAME
#undef USE_OUTPUT_BUFFER_PER_WARP
#define USE_PERIODIC
#define METHOD_NAME(a, b) a##Periodic##b
#include "kCalculateGBVISoftcoreForces2.h"
#define USE_OUTPUT_BUFFER_PER_WARP
#undef METHOD_NAME
#define METHOD_NAME(a, b) a##PeriodicByWarp##b
#include "kCalculateGBVISoftcoreForces2.h"

void kCalculateGBVISoftcoreForces2( freeEnergyGpuContext freeEnergyGpu )
{
    unsigned int threadsPerBlock;
    static unsigned int threadsPerBlockPerMethod[3] = { 0, 0, 0 };
    static unsigned int natoms[3]                   = { 0, 0, 0 };
    gpuContext gpu                                  = freeEnergyGpu->gpuContext;
    unsigned int methodIndex                        = static_cast<unsigned int>(freeEnergyGpu->freeEnergySim.nonbondedMethod);

    if( methodIndex > 2 ){
        throw OpenMM::OpenMMException( "kCalculateGBVISoftcoreForces2 method index invalid." );
    }   

    if( natoms[methodIndex] != gpu->natoms ){
        unsigned int extra                    = methodIndex == 0 ? 0 : sizeof(float3);
        threadsPerBlockPerMethod[methodIndex] = std::min(getThreadsPerBlockFEP( freeEnergyGpu, (sizeof(Atom) + extra), gpu->sharedMemoryPerBlock ), gpu->sim.nonbond_threads_per_block );
        natoms[methodIndex]                   = gpu->natoms;
    }
    threadsPerBlock = threadsPerBlockPerMethod[methodIndex];

#ifdef DEBUG
    fprintf( stderr,"kCalculateGBVISoftcoreForces2 nonbondedMethod=%d bornForce2_blocks=%u threadsPerBlock=%u shMem=%u\n",
             freeEnergyGpu->freeEnergySim.nonbondedMethod,
             gpu->sim.bornForce2_blocks, threadsPerBlock, (sizeof(Atom)+sizeof(float3))*threadsPerBlock ); fflush( stderr );
    
int psize = 64;
CUDAStream<float4>* pdE1 = new CUDAStream<float4>( psize, 1, "pdE");
CUDAStream<float4>* pdE2 = new CUDAStream<float4>( psize, 1, "pdE");
for( int ii = 0; ii < 32; ii++ ){

pdE1->_pSysData[ii].x = 0.0f;
pdE1->_pSysData[ii].y = 0.0f;
pdE1->_pSysData[ii].z = 0.0f;
pdE1->_pSysData[ii].w = 0.0f;

pdE2->_pSysData[ii].x = 0.0f;
pdE2->_pSysData[ii].y = 0.0f;
pdE2->_pSysData[ii].z = 0.0f;
pdE2->_pSysData[ii].w = 0.0f;
}
pdE1->Upload();
pdE2->Upload();

            if (gpu->bOutputBufferPerWarp)
                kCalculateGBVISoftcoreN2ByWarpForces2_kernel<<<gpu->sim.bornForce2_blocks, threadsPerBlock,
                        sizeof(Atom)*threadsPerBlock>>>(gpu->sim.pWorkUnit, gpu->sim.workUnits, pdE1->_pDevData, pdE2->_pDevData);
            else
                kCalculateGBVISoftcoreN2Forces2_kernel<<<gpu->sim.bornForce2_blocks, threadsPerBlock,
                        sizeof(Atom)*threadsPerBlock>>>(gpu->sim.pWorkUnit, gpu->sim.workUnits, pdE1->_pDevData, pdE2->_pDevData);
pdE1->Download();
pdE2->Download();
fprintf( stderr, "Pde\n" );
for( int ii = 0; ii < 32; ii++ ){
fprintf( stderr, "%4d %15.7e %15.7e %15.7e %15.7e    %15.7e %15.7e %15.7e %15.7e\n", ii, 
         pdE1->_pSysData[ii].x, pdE1->_pSysData[ii].y, pdE1->_pSysData[ii].z, pdE1->_pSysData[ii].w,
         pdE2->_pSysData[ii].x, pdE2->_pSysData[ii].y, pdE2->_pSysData[ii].z, pdE2->_pSysData[ii].w );
}
            break;
#endif

    switch (freeEnergyGpu->freeEnergySim.nonbondedMethod)
    {
        case FREE_ENERGY_NO_CUTOFF:


            if (gpu->bOutputBufferPerWarp)
                kCalculateGBVISoftcoreN2ByWarpForces2_kernel<<<gpu->sim.bornForce2_blocks, threadsPerBlock,
                        sizeof(Atom)*threadsPerBlock>>>(gpu->sim.pWorkUnit);
            else
                kCalculateGBVISoftcoreN2Forces2_kernel<<<gpu->sim.bornForce2_blocks, threadsPerBlock,
                        sizeof(Atom)*threadsPerBlock>>>(gpu->sim.pWorkUnit);
            break;

        case FREE_ENERGY_CUTOFF:

            if (gpu->bOutputBufferPerWarp)
                kCalculateGBVISoftcoreCutoffByWarpForces2_kernel<<<gpu->sim.bornForce2_blocks, threadsPerBlock,
                        (sizeof(Atom)+sizeof(float3))*threadsPerBlock>>>(gpu->sim.pInteractingWorkUnit );
            else
                kCalculateGBVISoftcoreCutoffForces2_kernel<<<gpu->sim.bornForce2_blocks, threadsPerBlock,
                        (sizeof(Atom)+sizeof(float3))*threadsPerBlock>>>(gpu->sim.pInteractingWorkUnit );
            break;

        case FREE_ENERGY_PERIODIC:

            if (gpu->bOutputBufferPerWarp)
                kCalculateGBVISoftcorePeriodicByWarpForces2_kernel<<<gpu->sim.bornForce2_blocks, threadsPerBlock,
                        (sizeof(Atom)+sizeof(float3))*threadsPerBlock>>>(gpu->sim.pInteractingWorkUnit );
            else
                kCalculateGBVISoftcorePeriodicForces2_kernel<<<gpu->sim.bornForce2_blocks, threadsPerBlock,
                        (sizeof(Atom)+sizeof(float3))*threadsPerBlock>>>(gpu->sim.pInteractingWorkUnit );
            break;

    }
    LAUNCHERROR("kCalculateGBVISoftcoreForces2");
}
