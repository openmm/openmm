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
using namespace std;

#include "gputypes.h"
#include "cudatypes.h"
#include "cudaKernels.h"

struct Atom {
    float x;
    float y;
    float z;
    float q;
    float sig;
    float eps;
    float br;
    float fx;
    float fy;
    float fz;
    float fb;
};

static __constant__ cudaGmxSimulation cSim;

void SetCalculateCDLJObcGbsaForces1Sim(gpuContext gpu)
{
    cudaError_t status;
    status = cudaMemcpyToSymbol(cSim, &gpu->sim, sizeof(cudaGmxSimulation));     
    RTERROR(status, "cudaMemcpyToSymbol: SetSim copy to cSim failed");
}

void GetCalculateCDLJObcGbsaForces1Sim(gpuContext gpu)
{
    cudaError_t status;
    status = cudaMemcpyFromSymbol(&gpu->sim, cSim, sizeof(cudaGmxSimulation));     
    RTERROR(status, "cudaMemcpyFromSymbol: SetSim copy from cSim failed");
}

// Include versions of the kernel for N^2 calculations.

#define METHOD_NAME(a, b) a##N2##b
#include "kCalculateCDLJObcGbsaForces1.h"
#define USE_OUTPUT_BUFFER_PER_WARP
#undef METHOD_NAME
#define METHOD_NAME(a, b) a##N2ByWarp##b
#include "kCalculateCDLJObcGbsaForces1.h"

// Include versions of the kernel with cutoffs.

#undef METHOD_NAME
#undef USE_OUTPUT_BUFFER_PER_WARP
#define USE_CUTOFF
#define METHOD_NAME(a, b) a##Cutoff##b
#include "kCalculateCDLJObcGbsaForces1.h"
#define USE_OUTPUT_BUFFER_PER_WARP
#undef METHOD_NAME
#define METHOD_NAME(a, b) a##CutoffByWarp##b
#include "kCalculateCDLJObcGbsaForces1.h"

// Include versions of the kernel with periodic boundary conditions.

#undef METHOD_NAME
#undef USE_OUTPUT_BUFFER_PER_WARP
#define USE_PERIODIC
#define METHOD_NAME(a, b) a##Periodic##b
#include "kCalculateCDLJObcGbsaForces1.h"
#define USE_OUTPUT_BUFFER_PER_WARP
#undef METHOD_NAME
#define METHOD_NAME(a, b) a##PeriodicByWarp##b
#include "kCalculateCDLJObcGbsaForces1.h"

extern __global__ void kFindBlockBoundsCutoff_kernel();
extern __global__ void kFindBlockBoundsPeriodic_kernel();
extern __global__ void kFindBlocksWithInteractionsCutoff_kernel();
extern __global__ void kFindBlocksWithInteractionsPeriodic_kernel();
extern __global__ void kFindInteractionsWithinBlocksCutoff_kernel(unsigned int*, int);
extern __global__ void kFindInteractionsWithinBlocksPeriodic_kernel(unsigned int*, int);

void kCalculateCDLJObcGbsaForces1(gpuContext gpu)
{
//    printf("kCalculateCDLJObcGbsaForces1\n");

    // check if Born radii need to be calculated

    kClearBornForces(gpu);
    CUDPPResult result;
    size_t numWithInteractions;
    switch (gpu->sim.nonbondedMethod)
    {
        case NO_CUTOFF:
            if (gpu->bRecalculateBornRadii)
            {
                kCalculateObcGbsaBornSum(gpu);
                kReduceObcGbsaBornSum(gpu);
            }
            if (gpu->bOutputBufferPerWarp)
                kCalculateCDLJObcGbsaN2ByWarpForces1_kernel<<<gpu->sim.nonbond_blocks, gpu->sim.nonbond_threads_per_block,
                        sizeof(Atom)*gpu->sim.nonbond_threads_per_block>>>(gpu->sim.pWorkUnit, gpu->sim.workUnits);
            else
                kCalculateCDLJObcGbsaN2Forces1_kernel<<<gpu->sim.nonbond_blocks, gpu->sim.nonbond_threads_per_block,
                        sizeof(Atom)*gpu->sim.nonbond_threads_per_block>>>(gpu->sim.pWorkUnit, gpu->sim.workUnits);
            LAUNCHERROR("kCalculateCDLJObcGbsaN2Forces1");
            break;
        case CUTOFF:
            kFindBlockBoundsCutoff_kernel<<<(gpu->psGridBoundingBox->_length+63)/64, 64>>>();
            LAUNCHERROR("kFindBlockBoundsCutoff");
            kFindBlocksWithInteractionsCutoff_kernel<<<gpu->sim.interaction_blocks, gpu->sim.interaction_threads_per_block>>>();
            LAUNCHERROR("kFindBlocksWithInteractionsCutoff");
            result = cudppCompact(gpu->cudpp, gpu->sim.pInteractingWorkUnit, gpu->sim.pInteractionCount,
                    gpu->sim.pWorkUnit, gpu->sim.pInteractionFlag, gpu->sim.workUnits);
            if (result != CUDPP_SUCCESS)
            {
                printf("Error in cudppCompact: %d\n", result);
                exit(-1);
            }
            gpu->psInteractionCount->Download();
            numWithInteractions = gpu->psInteractionCount->_pSysData[0];
            kFindInteractionsWithinBlocksCutoff_kernel<<<gpu->sim.nonbond_blocks, gpu->sim.nonbond_threads_per_block,
                    sizeof(unsigned int)*gpu->sim.nonbond_threads_per_block>>>(gpu->sim.pInteractingWorkUnit, numWithInteractions);
            if (gpu->bRecalculateBornRadii)
            {
                kCalculateObcGbsaBornSum(gpu);
                kReduceObcGbsaBornSum(gpu);
            }
            if (gpu->bOutputBufferPerWarp)
                kCalculateCDLJObcGbsaCutoffByWarpForces1_kernel<<<gpu->sim.nonbond_blocks, gpu->sim.nonbond_threads_per_block,
                        (sizeof(Atom)+sizeof(float))*gpu->sim.nonbond_threads_per_block>>>(gpu->sim.pInteractingWorkUnit, numWithInteractions);
            else
                kCalculateCDLJObcGbsaCutoffForces1_kernel<<<gpu->sim.nonbond_blocks, gpu->sim.nonbond_threads_per_block,
                        (sizeof(Atom)+sizeof(float))*gpu->sim.nonbond_threads_per_block>>>(gpu->sim.pInteractingWorkUnit, numWithInteractions);
            LAUNCHERROR("kCalculateCDLJObcGbsaCutoffForces1");
            break;
        case PERIODIC:
            kFindBlockBoundsPeriodic_kernel<<<(gpu->psGridBoundingBox->_length+63)/64, 64>>>();
            LAUNCHERROR("kFindBlockBoundsPeriodic");
            kFindBlocksWithInteractionsPeriodic_kernel<<<gpu->sim.interaction_blocks, gpu->sim.interaction_threads_per_block>>>();
            LAUNCHERROR("kFindBlocksWithInteractionsPeriodic");
            result = cudppCompact(gpu->cudpp, gpu->sim.pInteractingWorkUnit, gpu->sim.pInteractionCount,
                    gpu->sim.pWorkUnit, gpu->sim.pInteractionFlag, gpu->sim.workUnits);
            if (result != CUDPP_SUCCESS)
            {
                printf("Error in cudppCompact: %d\n", result);
                exit(-1);
            }
            gpu->psInteractionCount->Download();
            numWithInteractions = gpu->psInteractionCount->_pSysData[0];
            kFindInteractionsWithinBlocksPeriodic_kernel<<<gpu->sim.nonbond_blocks, gpu->sim.nonbond_threads_per_block,
                    sizeof(unsigned int)*gpu->sim.nonbond_threads_per_block>>>(gpu->sim.pInteractingWorkUnit, numWithInteractions);
            if (gpu->bRecalculateBornRadii)
            {
                kCalculateObcGbsaBornSum(gpu);
                kReduceObcGbsaBornSum(gpu);
            }
            if (gpu->bOutputBufferPerWarp)
                kCalculateCDLJObcGbsaPeriodicByWarpForces1_kernel<<<gpu->sim.nonbond_blocks, gpu->sim.nonbond_threads_per_block,
                        (sizeof(Atom)+sizeof(float))*gpu->sim.nonbond_threads_per_block>>>(gpu->sim.pInteractingWorkUnit, numWithInteractions);
            else
                kCalculateCDLJObcGbsaPeriodicForces1_kernel<<<gpu->sim.nonbond_blocks, gpu->sim.nonbond_threads_per_block,
                        (sizeof(Atom)+sizeof(float))*gpu->sim.nonbond_threads_per_block>>>(gpu->sim.pInteractingWorkUnit, numWithInteractions);
            LAUNCHERROR("kCalculateCDLJObcGbsaPeriodicForces1");
            break;
    }
}
