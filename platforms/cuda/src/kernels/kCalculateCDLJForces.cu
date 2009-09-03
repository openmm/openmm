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
#include "cudatypes.h"

#define UNROLLXX 0
#define UNROLLXY 0

struct Atom {
    float x;
    float y;
    float z;
    float q;
    float sig;
    float eps;
    float fx;
    float fy;
    float fz;
};

static __constant__ cudaGmxSimulation cSim;

void SetCalculateCDLJForcesSim(gpuContext gpu)
{
    cudaError_t status;
    status = cudaMemcpyToSymbol(cSim, &gpu->sim, sizeof(cudaGmxSimulation));
    RTERROR(status, "cudaMemcpyToSymbol: SetSim copy to cSim failed");
}

void GetCalculateCDLJForcesSim(gpuContext gpu)
{
    cudaError_t status;
    status = cudaMemcpyFromSymbol(&gpu->sim, cSim, sizeof(cudaGmxSimulation));
    RTERROR(status, "cudaMemcpyFromSymbol: SetSim copy from cSim failed");
}

// Include versions of the kernels for N^2 calculations.

#define METHOD_NAME(a, b) a##N2##b
#include "kCalculateCDLJForces.h"
#define USE_OUTPUT_BUFFER_PER_WARP
#undef METHOD_NAME
#define METHOD_NAME(a, b) a##N2ByWarp##b
#include "kCalculateCDLJForces.h"

// Include versions of the kernels with cutoffs.

#undef METHOD_NAME
#undef USE_OUTPUT_BUFFER_PER_WARP
#define USE_CUTOFF
#define METHOD_NAME(a, b) a##Cutoff##b
#include "kCalculateCDLJForces.h"
#include "kFindInteractingBlocks.h"
#define USE_OUTPUT_BUFFER_PER_WARP
#undef METHOD_NAME
#define METHOD_NAME(a, b) a##CutoffByWarp##b
#include "kCalculateCDLJForces.h"

// Include versions of the kernels with periodic boundary conditions.

#undef METHOD_NAME
#undef USE_OUTPUT_BUFFER_PER_WARP
#define USE_PERIODIC
#define METHOD_NAME(a, b) a##Periodic##b
#include "kCalculateCDLJForces.h"
#include "kFindInteractingBlocks.h"
#define USE_OUTPUT_BUFFER_PER_WARP
#undef METHOD_NAME
#define METHOD_NAME(a, b) a##PeriodicByWarp##b
#include "kCalculateCDLJForces.h"

// Include versions of the kernels for Ewald

#undef METHOD_NAME
#undef USE_OUTPUT_BUFFER_PER_WARP
#define USE_PERIODIC
#define USE_EWALD
#define METHOD_NAME(a, b) a##Ewald##b
#include "kCalculateCDLJForces.h"
#define USE_OUTPUT_BUFFER_PER_WARP
#undef METHOD_NAME
#define METHOD_NAME(a, b) a##EwaldByWarp##b
#include "kCalculateCDLJForces.h"

// Reciprocal Space Ewald summation is in a separate kernel
#include "kCalculateCDLJEwaldFastReciprocal.h"

void kCalculatePME(gpuContext gpu);

void kCalculateCDLJForces(gpuContext gpu)
{
//    printf("kCalculateCDLJCutoffForces\n");
    switch (gpu->sim.nonbondedMethod)
    {
        case NO_CUTOFF:
            if (gpu->bOutputBufferPerWarp)
                kCalculateCDLJN2ByWarpForces_kernel<<<gpu->sim.nonbond_blocks, gpu->sim.nonbond_threads_per_block,
                        sizeof(Atom)*gpu->sim.nonbond_threads_per_block>>>(gpu->sim.pWorkUnit);
            else
                kCalculateCDLJN2Forces_kernel<<<gpu->sim.nonbond_blocks, gpu->sim.nonbond_threads_per_block,
                        sizeof(Atom)*gpu->sim.nonbond_threads_per_block>>>(gpu->sim.pWorkUnit);
            LAUNCHERROR("kCalculateCDLJN2Forces");
            break;
        case CUTOFF:
            kFindBlockBoundsCutoff_kernel<<<(gpu->psGridBoundingBox->_length+63)/64, 64>>>();
            LAUNCHERROR("kFindBlockBoundsCutoff");
            kFindBlocksWithInteractionsCutoff_kernel<<<gpu->sim.interaction_blocks, gpu->sim.interaction_threads_per_block>>>();
            LAUNCHERROR("kFindBlocksWithInteractionsCutoff");
            compactStream(gpu->compactPlan, gpu->sim.pInteractingWorkUnit, gpu->sim.pWorkUnit, gpu->sim.pInteractionFlag, gpu->sim.workUnits, gpu->sim.pInteractionCount);
            kFindInteractionsWithinBlocksCutoff_kernel<<<gpu->sim.nonbond_blocks, gpu->sim.nonbond_threads_per_block,
                    sizeof(unsigned int)*gpu->sim.nonbond_threads_per_block>>>(gpu->sim.pInteractingWorkUnit);
            if (gpu->bOutputBufferPerWarp)
                kCalculateCDLJCutoffByWarpForces_kernel<<<gpu->sim.nonbond_blocks, gpu->sim.nonbond_threads_per_block,
                        (sizeof(Atom)+sizeof(float3))*gpu->sim.nonbond_threads_per_block>>>(gpu->sim.pInteractingWorkUnit);
            else
                kCalculateCDLJCutoffForces_kernel<<<gpu->sim.nonbond_blocks, gpu->sim.nonbond_threads_per_block,
                        (sizeof(Atom)+sizeof(float3))*gpu->sim.nonbond_threads_per_block>>>(gpu->sim.pInteractingWorkUnit);
            LAUNCHERROR("kCalculateCDLJCutoffForces");
            break;
        case PERIODIC:
            kFindBlockBoundsPeriodic_kernel<<<(gpu->psGridBoundingBox->_length+63)/64, 64>>>();
            LAUNCHERROR("kFindBlockBoundsPeriodic");
            kFindBlocksWithInteractionsPeriodic_kernel<<<gpu->sim.interaction_blocks, gpu->sim.interaction_threads_per_block>>>();
            LAUNCHERROR("kFindBlocksWithInteractionsPeriodic");
            compactStream(gpu->compactPlan, gpu->sim.pInteractingWorkUnit, gpu->sim.pWorkUnit, gpu->sim.pInteractionFlag, gpu->sim.workUnits, gpu->sim.pInteractionCount);
            kFindInteractionsWithinBlocksPeriodic_kernel<<<gpu->sim.nonbond_blocks, gpu->sim.nonbond_threads_per_block,
                    sizeof(unsigned int)*gpu->sim.nonbond_threads_per_block>>>(gpu->sim.pInteractingWorkUnit);
            if (gpu->bOutputBufferPerWarp)
                kCalculateCDLJPeriodicByWarpForces_kernel<<<gpu->sim.nonbond_blocks, gpu->sim.nonbond_threads_per_block,
                        (sizeof(Atom)+sizeof(float3))*gpu->sim.nonbond_threads_per_block>>>(gpu->sim.pInteractingWorkUnit);
            else
                kCalculateCDLJPeriodicForces_kernel<<<gpu->sim.nonbond_blocks, gpu->sim.nonbond_threads_per_block,
                        (sizeof(Atom)+sizeof(float3))*gpu->sim.nonbond_threads_per_block>>>(gpu->sim.pInteractingWorkUnit);
            LAUNCHERROR("kCalculateCDLJPeriodicForces");
            break;
        case EWALD:
            kFindBlockBoundsPeriodic_kernel<<<(gpu->psGridBoundingBox->_length+63)/64, 64>>>();
            LAUNCHERROR("kFindBlockBoundsPeriodic");
            kFindBlocksWithInteractionsPeriodic_kernel<<<gpu->sim.interaction_blocks, gpu->sim.interaction_threads_per_block>>>();
            LAUNCHERROR("kFindBlocksWithInteractionsPeriodic");
            compactStream(gpu->compactPlan, gpu->sim.pInteractingWorkUnit, gpu->sim.pWorkUnit, gpu->sim.pInteractionFlag, gpu->sim.workUnits, gpu->sim.pInteractionCount);
            kFindInteractionsWithinBlocksPeriodic_kernel<<<gpu->sim.nonbond_blocks, gpu->sim.nonbond_threads_per_block,
                    sizeof(unsigned int)*gpu->sim.nonbond_threads_per_block>>>(gpu->sim.pInteractingWorkUnit);
            LAUNCHERROR("kFindInteractionsWithinBlocksPeriodic");
            if (gpu->bOutputBufferPerWarp)
                kCalculateCDLJEwaldByWarpForces_kernel<<<gpu->sim.nonbond_blocks, gpu->sim.nonbond_threads_per_block,
                        (sizeof(Atom)+sizeof(float3))*gpu->sim.nonbond_threads_per_block>>>(gpu->sim.pInteractingWorkUnit);
            else
                kCalculateCDLJEwaldForces_kernel<<<gpu->sim.nonbond_blocks, gpu->sim.nonbond_threads_per_block,
                        (sizeof(Atom)+sizeof(float3))*gpu->sim.nonbond_threads_per_block>>>(gpu->sim.pInteractingWorkUnit);
            LAUNCHERROR("kCalculateCDLJEwaldForces");
            // Ewald summation
            kCalculateEwaldFastCosSinSums_kernel<<<gpu->sim.nonbond_blocks, gpu->sim.nonbond_threads_per_block>>>();
            LAUNCHERROR("kCalculateEwaldFastCosSinSums");
            kCalculateEwaldFastForces_kernel<<<gpu->sim.blocks, gpu->sim.update_threads_per_block>>>();
            LAUNCHERROR("kCalculateEwaldFastForces");
            break;
        case PARTICLE_MESH_EWALD:
            kFindBlockBoundsPeriodic_kernel<<<(gpu->psGridBoundingBox->_length+63)/64, 64>>>();
            LAUNCHERROR("kFindBlockBoundsPeriodic");
            kFindBlocksWithInteractionsPeriodic_kernel<<<gpu->sim.interaction_blocks, gpu->sim.interaction_threads_per_block>>>();
            LAUNCHERROR("kFindBlocksWithInteractionsPeriodic");
            compactStream(gpu->compactPlan, gpu->sim.pInteractingWorkUnit, gpu->sim.pWorkUnit, gpu->sim.pInteractionFlag, gpu->sim.workUnits, gpu->sim.pInteractionCount);
            kFindInteractionsWithinBlocksPeriodic_kernel<<<gpu->sim.nonbond_blocks, gpu->sim.nonbond_threads_per_block,
                    sizeof(unsigned int)*gpu->sim.nonbond_threads_per_block>>>(gpu->sim.pInteractingWorkUnit);
            LAUNCHERROR("kFindInteractionsWithinBlocksPeriodic");
            if (gpu->bOutputBufferPerWarp)
                kCalculateCDLJEwaldByWarpForces_kernel<<<gpu->sim.nonbond_blocks, gpu->sim.nonbond_threads_per_block,
                        (sizeof(Atom)+sizeof(float3))*gpu->sim.nonbond_threads_per_block>>>(gpu->sim.pInteractingWorkUnit);
            else
                kCalculateCDLJEwaldForces_kernel<<<gpu->sim.nonbond_blocks, gpu->sim.nonbond_threads_per_block,
                        (sizeof(Atom)+sizeof(float3))*gpu->sim.nonbond_threads_per_block>>>(gpu->sim.pInteractingWorkUnit);
            LAUNCHERROR("kCalculateCDLJEwaldForces");
            kCalculatePME(gpu);
    }
}
