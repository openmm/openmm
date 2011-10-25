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

struct Atom {
    float x;
    float y;
    float z;
    float r;
    float sr;
    float sum;
    float padding;
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

__global__ 
__launch_bounds__(384, 1)
void kReduceObcGbsaBornSum_kernel()
{
    unsigned int pos = (blockIdx.x * blockDim.x + threadIdx.x);
    
    while (pos < cSim.atoms)
    {
        float sum   = 0.0f;
        float* pSt  = cSim.pBornSum + pos;
        float2 atom = cSim.pObcData[pos];
        
        // Get summed Born data
        for (int i = 0; i < cSim.nonbondOutputBuffers; i++)
        {
            sum += *pSt;
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
        cSim.pBornRadii[pos]    = bornRadius;
        cSim.pObcChain[pos]     = obcChain;
        pos                    += gridDim.x * blockDim.x;
    }   
}

void OPENMMCUDA_EXPORT kReduceObcGbsaBornSum(gpuContext gpu)
{
    kReduceObcGbsaBornSum_kernel<<<gpu->sim.blocks, 384>>>();
    gpu->bRecalculateBornRadii = false;
    LAUNCHERROR("kReduceObcGbsaBornSum");
}

void kPrintObc( gpuContext gpu, std::string callId, int call, FILE* log)
{

    gpu->psObcData->Download();
    gpu->psBornRadii->Download();
    gpu->psObcChain->Download();
    gpu->psBornForce->Download();
    gpu->psPosq4->Download();
    gpu->psSigEps2->Download();

    (void) fprintf( log, "kPrintObc Cuda bCh bR bF prm[2]   sigeps[2]\n" );
    (void) fprintf( stderr, "bOutputWarp=%u blks=%u th/blk=%u wu=%u %u shrd=%u\n", gpu->bOutputBufferPerWarp,
                    gpu->sim.nonbond_blocks, gpu->sim.nonbond_threads_per_block, gpu->sim.workUnits, gpu->psWorkUnit->_pSysStream[0][0],
                    sizeof(Atom)*gpu->sim.nonbond_threads_per_block );
    for( int ii = 0; ii < gpu->sim.paddedNumberOfAtoms; ii++ ){
        (void) fprintf( log, "%6d %15.7e %15.7e %15.7e    %15.7e %15.7e   %15.7e %15.7e \n", ii, 
                        gpu->psObcChain->_pSysData[ii],
                        gpu->psBornRadii->_pSysData[ii],
                        gpu->psBornForce->_pSysData[ii],

                        gpu->psObcData->_pSysData[ii].x,
                        gpu->psObcData->_pSysData[ii].y,

                        gpu->psSigEps2->_pSysData[ii].x,
                        gpu->psSigEps2->_pSysData[ii].y );

    }   

}

void OPENMMCUDA_EXPORT kCalculateObcGbsaBornSum(gpuContext gpu)
{
  //  printf("kCalculateObcgbsaBornSum\n");
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
