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
extern __global__ void kFindInteractionsWithinBlocksCutoff_kernel(unsigned int*);
extern __global__ void kFindInteractionsWithinBlocksPeriodic_kernel(unsigned int*);

void kCalculateCDLJObcGbsaForces1(gpuContext gpu)
{
//    printf("kCalculateCDLJObcGbsaForces1\n");

    // check if Born radii need to be calculated

    kClearBornForces(gpu);
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
                        sizeof(Atom)*gpu->sim.nonbond_threads_per_block>>>(gpu->sim.pWorkUnit);
            else
                kCalculateCDLJObcGbsaN2Forces1_kernel<<<gpu->sim.nonbond_blocks, gpu->sim.nonbond_threads_per_block,
                        sizeof(Atom)*gpu->sim.nonbond_threads_per_block>>>(gpu->sim.pWorkUnit);
            LAUNCHERROR("kCalculateCDLJObcGbsaN2Forces1");
            break;
        case CUTOFF:
            kFindBlockBoundsCutoff_kernel<<<(gpu->psGridBoundingBox->_length+63)/64, 64>>>();
            LAUNCHERROR("kFindBlockBoundsCutoff");
            kFindBlocksWithInteractionsCutoff_kernel<<<gpu->sim.interaction_blocks, gpu->sim.interaction_threads_per_block>>>();
            LAUNCHERROR("kFindBlocksWithInteractionsCutoff");
            compactStream(gpu->compactPlan, gpu->sim.pInteractingWorkUnit, gpu->sim.pWorkUnit, gpu->sim.pInteractionFlag, gpu->sim.workUnits, gpu->sim.pInteractionCount);
            kFindInteractionsWithinBlocksCutoff_kernel<<<gpu->sim.nonbond_blocks, gpu->sim.nonbond_threads_per_block,
                    sizeof(unsigned int)*gpu->sim.nonbond_threads_per_block>>>(gpu->sim.pInteractingWorkUnit);
            if (gpu->bRecalculateBornRadii)
            {
                kCalculateObcGbsaBornSum(gpu);
                kReduceObcGbsaBornSum(gpu);
            }
            if (gpu->bOutputBufferPerWarp)
                kCalculateCDLJObcGbsaCutoffByWarpForces1_kernel<<<gpu->sim.nonbond_blocks, gpu->sim.nonbond_threads_per_block,
                        (sizeof(Atom)+sizeof(float))*gpu->sim.nonbond_threads_per_block>>>(gpu->sim.pInteractingWorkUnit);
            else
                kCalculateCDLJObcGbsaCutoffForces1_kernel<<<gpu->sim.nonbond_blocks, gpu->sim.nonbond_threads_per_block,
                        (sizeof(Atom)+sizeof(float))*gpu->sim.nonbond_threads_per_block>>>(gpu->sim.pInteractingWorkUnit);
            LAUNCHERROR("kCalculateCDLJObcGbsaCutoffForces1");
            break;
        case PERIODIC:
            kFindBlockBoundsPeriodic_kernel<<<(gpu->psGridBoundingBox->_length+63)/64, 64>>>();
            LAUNCHERROR("kFindBlockBoundsPeriodic");
            kFindBlocksWithInteractionsPeriodic_kernel<<<gpu->sim.interaction_blocks, gpu->sim.interaction_threads_per_block>>>();
            LAUNCHERROR("kFindBlocksWithInteractionsPeriodic");
            compactStream(gpu->compactPlan, gpu->sim.pInteractingWorkUnit, gpu->sim.pWorkUnit, gpu->sim.pInteractionFlag, gpu->sim.workUnits, gpu->sim.pInteractionCount);
            kFindInteractionsWithinBlocksPeriodic_kernel<<<gpu->sim.nonbond_blocks, gpu->sim.nonbond_threads_per_block,
                    sizeof(unsigned int)*gpu->sim.nonbond_threads_per_block>>>(gpu->sim.pInteractingWorkUnit);
            if (gpu->bRecalculateBornRadii)
            {
                kCalculateObcGbsaBornSum(gpu);
                kReduceObcGbsaBornSum(gpu);
            }
            if (gpu->bOutputBufferPerWarp)
                kCalculateCDLJObcGbsaPeriodicByWarpForces1_kernel<<<gpu->sim.nonbond_blocks, gpu->sim.nonbond_threads_per_block,
                        (sizeof(Atom)+sizeof(float))*gpu->sim.nonbond_threads_per_block>>>(gpu->sim.pInteractingWorkUnit);
            else
                kCalculateCDLJObcGbsaPeriodicForces1_kernel<<<gpu->sim.nonbond_blocks, gpu->sim.nonbond_threads_per_block,
                        (sizeof(Atom)+sizeof(float))*gpu->sim.nonbond_threads_per_block>>>(gpu->sim.pInteractingWorkUnit);
            LAUNCHERROR("kCalculateCDLJObcGbsaPeriodicForces1");
            break;
    }
}
