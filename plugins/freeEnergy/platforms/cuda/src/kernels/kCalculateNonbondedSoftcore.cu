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

#include "freeEnergyGpuTypes.h"
#include "GpuFreeEnergyCudaKernels.h"
#include "openmm/OpenMMException.h"
#include <iostream>
#include <sstream>

#define PARAMETER_PRINT 0
#define MAX_PARAMETER_PRINT 10

// device handles

static __constant__ cudaGmxSimulation cSim;
static __constant__ cudaFreeEnergyGmxSimulation feSimDev;

// write address of structs to devices

void SetCalculateCDLJSoftcoreGpuSim( freeEnergyGpuContext freeEnergyGpu ){
    cudaError_t status;
    status = cudaMemcpyToSymbol(cSim, &freeEnergyGpu->gpuContext->sim, sizeof(cudaGmxSimulation));
    RTERROR(status, "cudaMemcpyToSymbol: SetCalculateCDLJSoftcoreGpuSim copy to cSim failed");

    status = cudaMemcpyToSymbol( feSimDev, &freeEnergyGpu->freeEnergySim, sizeof(cudaFreeEnergyGmxSimulation));
    RTERROR(status, "cudaMemcpyToSymbol: SetCalculateCDLJSoftcoreGpuSim copy to feSimDev failed");
}

extern "C"
bool gpuIsAvailableSoftcore()
{
    int deviceCount;
    cudaGetDeviceCount(&deviceCount);
    return (deviceCount > 0);
}

struct Atom {
    float x;
    float y;
    float z;
    float q;
    float sig;
    float eps;
    float softCoreLJLambda;
    float fx;
    float fy;
    float fz;
};

// Include versions of the kernels for N^2 calculations with softcore LJ.

#define USE_SOFTCORE_LJ
#ifdef USE_SOFTCORE_LJ
#include "kSoftcoreLJ.h"
#endif

#undef METHOD_NAME
#define METHOD_NAME(a, b) a##N2SoftcoreLJ##b
#undef USE_OUTPUT_BUFFER_PER_WARP
#include "kCalculateNonbondedSoftcore.h"

#define USE_OUTPUT_BUFFER_PER_WARP
#undef METHOD_NAME
#define METHOD_NAME(a, b) a##N2SoftcoreLJByWarp##b
#include "kCalculateNonbondedSoftcore.h"

// Include versions of the kernels with cutoffs.

#undef METHOD_NAME
#undef USE_OUTPUT_BUFFER_PER_WARP
#define USE_CUTOFF
#define METHOD_NAME(a, b) a##Cutoff##b
#include "kCalculateNonbondedSoftcore.h"
#include "kFindInteractingBlocks.h"
#define USE_OUTPUT_BUFFER_PER_WARP
#undef METHOD_NAME
#define METHOD_NAME(a, b) a##CutoffByWarp##b
#include "kCalculateNonbondedSoftcore.h"

// Include versions of the kernels with periodic boundary conditions.

#undef METHOD_NAME
#undef USE_OUTPUT_BUFFER_PER_WARP
#define USE_PERIODIC
#define METHOD_NAME(a, b) a##Periodic##b
#include "kCalculateNonbondedSoftcore.h"
#include "kFindInteractingBlocks.h"
#define USE_OUTPUT_BUFFER_PER_WARP
#undef METHOD_NAME
#define METHOD_NAME(a, b) a##PeriodicByWarp##b
#include "kCalculateNonbondedSoftcore.h"

void kCalculateCDLJSoftcoreForces( freeEnergyGpuContext freeEnergyGpu )
{

    gpuContext gpu = freeEnergyGpu->gpuContext;
    switch (freeEnergyGpu->freeEnergySim.nonbondedMethod)
    {
        case FREE_ENERGY_NO_CUTOFF:

           if (gpu->bOutputBufferPerWarp)
                kCalculateCDLJSoftcoreN2SoftcoreLJByWarpForces_kernel<<<gpu->sim.nonbond_blocks, gpu->sim.nonbond_threads_per_block,
                         sizeof(Atom)*gpu->sim.nonbond_threads_per_block>>>(gpu->sim.pWorkUnit);
           else
                   kCalculateCDLJSoftcoreN2SoftcoreLJForces_kernel<<<gpu->sim.nonbond_blocks, gpu->sim.nonbond_threads_per_block,
                           sizeof(Atom)*gpu->sim.nonbond_threads_per_block>>>(gpu->sim.pWorkUnit );
            LAUNCHERROR("kCalculateCDLJSoftcoreN2Forces");

            break;

        case FREE_ENERGY_CUTOFF:

            kFindBlockBoundsCutoff_kernel<<<(gpu->psGridBoundingBox->_length+63)/64, 64>>>();
            LAUNCHERROR("kFindBlockBoundsCutoff");
            kFindBlocksWithInteractionsCutoff_kernel<<<gpu->sim.interaction_blocks, gpu->sim.interaction_threads_per_block>>>();
            LAUNCHERROR("kFindBlocksWithInteractionsCutoff");
            compactStream(gpu->compactPlan, gpu->sim.pInteractingWorkUnit, gpu->sim.pWorkUnit, gpu->sim.pInteractionFlag, gpu->sim.workUnits, gpu->sim.pInteractionCount);
            kFindInteractionsWithinBlocksCutoff_kernel<<<gpu->sim.nonbond_blocks, gpu->sim.nonbond_threads_per_block,
                    sizeof(unsigned int)*gpu->sim.nonbond_threads_per_block>>>(gpu->sim.pInteractingWorkUnit);

            if (gpu->bOutputBufferPerWarp)
                kCalculateCDLJSoftcoreCutoffByWarpForces_kernel<<<gpu->sim.nonbond_blocks, gpu->sim.nonbond_threads_per_block,
                        (sizeof(Atom)+sizeof(float3))*gpu->sim.nonbond_threads_per_block>>>(gpu->sim.pInteractingWorkUnit);
            else
                kCalculateCDLJSoftcoreCutoffForces_kernel<<<gpu->sim.nonbond_blocks, gpu->sim.nonbond_threads_per_block,
                        (sizeof(Atom)+sizeof(float3))*gpu->sim.nonbond_threads_per_block>>>(gpu->sim.pInteractingWorkUnit);

            LAUNCHERROR("kCalculateCDLJSoftcoreCutoffForces");
            break;

        case FREE_ENERGY_PERIODIC:

            kFindBlockBoundsPeriodic_kernel<<<(gpu->psGridBoundingBox->_length+63)/64, 64>>>();
            LAUNCHERROR("kFindBlockBoundsPeriodic");
            kFindBlocksWithInteractionsPeriodic_kernel<<<gpu->sim.interaction_blocks, gpu->sim.interaction_threads_per_block>>>();
            LAUNCHERROR("kFindBlocksWithInteractionsPeriodic");
            compactStream(gpu->compactPlan, gpu->sim.pInteractingWorkUnit, gpu->sim.pWorkUnit, gpu->sim.pInteractionFlag, gpu->sim.workUnits, gpu->sim.pInteractionCount);
            kFindInteractionsWithinBlocksPeriodic_kernel<<<gpu->sim.nonbond_blocks, gpu->sim.nonbond_threads_per_block,
                    sizeof(unsigned int)*gpu->sim.nonbond_threads_per_block>>>(gpu->sim.pInteractingWorkUnit);
            if (gpu->bOutputBufferPerWarp)
                kCalculateCDLJSoftcorePeriodicByWarpForces_kernel<<<gpu->sim.nonbond_blocks, gpu->sim.nonbond_threads_per_block,
                        (sizeof(Atom)+sizeof(float3))*gpu->sim.nonbond_threads_per_block>>>(gpu->sim.pInteractingWorkUnit);
            else
                kCalculateCDLJSoftcorePeriodicForces_kernel<<<gpu->sim.nonbond_blocks, gpu->sim.nonbond_threads_per_block,
                        (sizeof(Atom)+sizeof(float3))*gpu->sim.nonbond_threads_per_block>>>(gpu->sim.pInteractingWorkUnit);
            LAUNCHERROR("kCalculateCDLJSoftcorePeriodicForces");
            break;

        default:
            throw OpenMM::OpenMMException( "Nonbonded softcore method not recognized." );

    }
}

