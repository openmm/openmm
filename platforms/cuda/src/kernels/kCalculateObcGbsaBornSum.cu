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
 * This program is free software: you can redistribute it and/or modify       *
 * it under the terms of the GNU Lesser General Public License as published   *
 * by the Free Software Foundation, either version 3 of the License, or       *
 * (at your option) any later version.                                        *
 *                                                                            *
 * This program is distributed in the hope that it will be useful,            *
 * but WITHOUT ANY WARRANTY; without even the implied warranty of             *
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the              *
 * GNU Lesser General Public License for more details.                        *
 *                                                                            *
 * You should have received a copy of the GNU Lesser General Public License   *
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.      *
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
};

static __constant__ cudaGmxSimulation cSim;

void SetCalculateObcGbsaBornSumSim(gpuContext gpu)
{
    cudaError_t status;
    status = cudaMemcpyToSymbol(cSim, &gpu->sim, sizeof(cudaGmxSimulation));     
    RTERROR(status, "cudaMemcpyToSymbol: SetSim copy to cSim failed");
}

void GetCalculateObcGbsaBornSumSim(gpuContext gpu)
{
    cudaError_t status;
    status = cudaMemcpyFromSymbol(&gpu->sim, cSim, sizeof(cudaGmxSimulation));     
    RTERROR(status, "cudaMemcpyFromSymbol: SetSim copy from cSim failed");
}

// Include versions of the kernels for N^2 calculations.

#define METHOD_NAME(a, b) a##N2##b
#include "kCalculateObcGbsaBornSum.h"
#define USE_OUTPUT_BUFFER_PER_WARP
#undef METHOD_NAME
#define METHOD_NAME(a, b) a##N2ByWarp##b
#include "kCalculateObcGbsaBornSum.h"

// Include versions of the kernels with cutoffs.

#undef METHOD_NAME
#undef USE_OUTPUT_BUFFER_PER_WARP
#define USE_CUTOFF
#define METHOD_NAME(a, b) a##Cutoff##b
#include "kCalculateObcGbsaBornSum.h"
#define USE_OUTPUT_BUFFER_PER_WARP
#undef METHOD_NAME
#define METHOD_NAME(a, b) a##CutoffByWarp##b
#include "kCalculateObcGbsaBornSum.h"

// Include versions of the kernels with periodic boundary conditions.

#undef METHOD_NAME
#undef USE_OUTPUT_BUFFER_PER_WARP
#define USE_PERIODIC
#define METHOD_NAME(a, b) a##Periodic##b
#include "kCalculateObcGbsaBornSum.h"
#define USE_OUTPUT_BUFFER_PER_WARP
#undef METHOD_NAME
#define METHOD_NAME(a, b) a##PeriodicByWarp##b
#include "kCalculateObcGbsaBornSum.h"


__global__ void kClearObcGbsaBornSum_kernel()
{
    unsigned int pos = blockIdx.x * blockDim.x + threadIdx.x;
    while (pos < cSim.stride * cSim.nonbondOutputBuffers)
    {
        ((float*)cSim.pBornSum)[pos] = 0.0f;
        pos += gridDim.x * blockDim.x;
    }
}

__global__ void kReduceObcGbsaBornSum_kernel()
{
    unsigned int pos = (blockIdx.x * blockDim.x + threadIdx.x);
    
    while (pos < cSim.atoms)
    {
        float sum = 0.0f;
        float* pSt = cSim.pBornSum + pos;
        float2 atom = cSim.pObcData[pos];
        
        // Get summed Born data
        for (int i = 0; i < cSim.nonbondOutputBuffers; i++)
        {
            sum += *pSt;
       //     printf("%4d %4d A: %9.4f\n", pos, i, *pSt);
            pSt += cSim.stride;
        }
        
        
        // Now calculate Born radius and OBC term.
        sum                    *= 0.5f * atom.x;
        float sum2              = sum * sum;
        float sum3              = sum * sum2;
        float tanhSum           = tanh(cSim.alphaOBC * sum - cSim.betaOBC * sum2 + cSim.gammaOBC * sum3);
        float nonOffsetRadii    = atom.x + cSim.dielectricOffset;
        float bornRadius        = 1.0f / (1.0f / atom.x - tanhSum / nonOffsetRadii); 
        float obcChain          = atom.x * (cSim.alphaOBC - 2.0f * cSim.betaOBC * sum + 3.0f * cSim.gammaOBC * sum2);
        obcChain                = (1.0f - tanhSum * tanhSum) * obcChain / nonOffsetRadii;        
        cSim.pBornRadii[pos] = bornRadius;
        cSim.pObcChain[pos]  = obcChain;
        pos += gridDim.x * blockDim.x;
    }   
}

void kReduceObcGbsaBornSum(gpuContext gpu)
{
//    printf("kReduceObcGbsaBornSum\n");
    kReduceObcGbsaBornSum_kernel<<<gpu->sim.blocks, 384>>>();
    gpu->bRecalculateBornRadii = false;
    LAUNCHERROR("kReduceObcGbsaBornSum");
}

void kClearObcGbsaBornSum(gpuContext gpu)
{
  //  printf("kClearObcGbsaBornSum\n");
    kClearObcGbsaBornSum_kernel<<<gpu->sim.blocks, 384>>>();
}

void kCalculateObcGbsaBornSum(gpuContext gpu)
{
  //  printf("kCalculateObcgbsaBornSum\n");
    kClearObcGbsaBornSum(gpu);
    LAUNCHERROR("kClearBornSum");
    switch (gpu->sim.nonbondedMethod)
    {
        case NO_CUTOFF:

            if (gpu->bOutputBufferPerWarp)
                kCalculateObcGbsaN2ByWarpBornSum_kernel<<<gpu->sim.nonbond_blocks, gpu->sim.nonbond_threads_per_block,
                        sizeof(Atom)*gpu->sim.nonbond_threads_per_block>>>(gpu->sim.pWorkUnit);
            else
                kCalculateObcGbsaN2BornSum_kernel<<<gpu->sim.nonbond_blocks, gpu->sim.nonbond_threads_per_block,
                        sizeof(Atom)*gpu->sim.nonbond_threads_per_block>>>(gpu->sim.pWorkUnit);
            break;
        case CUTOFF:
            if (gpu->bOutputBufferPerWarp)
                kCalculateObcGbsaCutoffByWarpBornSum_kernel<<<gpu->sim.nonbond_blocks, gpu->sim.nonbond_threads_per_block,
                        (sizeof(Atom)+sizeof(float))*gpu->sim.nonbond_threads_per_block>>>(gpu->sim.pInteractingWorkUnit);
            else
                kCalculateObcGbsaCutoffBornSum_kernel<<<gpu->sim.nonbond_blocks, gpu->sim.nonbond_threads_per_block,
                        (sizeof(Atom)+sizeof(float))*gpu->sim.nonbond_threads_per_block>>>(gpu->sim.pInteractingWorkUnit);
            break;
        case PERIODIC:
            if (gpu->bOutputBufferPerWarp)
                kCalculateObcGbsaPeriodicByWarpBornSum_kernel<<<gpu->sim.nonbond_blocks, gpu->sim.nonbond_threads_per_block,
                        (sizeof(Atom)+sizeof(float))*gpu->sim.nonbond_threads_per_block>>>(gpu->sim.pInteractingWorkUnit);
            else
                kCalculateObcGbsaPeriodicBornSum_kernel<<<gpu->sim.nonbond_blocks, gpu->sim.nonbond_threads_per_block,
                        (sizeof(Atom)+sizeof(float))*gpu->sim.nonbond_threads_per_block>>>(gpu->sim.pInteractingWorkUnit);
            break;
    }
    LAUNCHERROR("kCalculateBornSum");
}
