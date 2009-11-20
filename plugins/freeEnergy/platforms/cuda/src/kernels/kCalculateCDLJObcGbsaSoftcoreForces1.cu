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

#include "kernels/gputypes.h"
#include "kernels/cudatypes.h"
#include "kernels/cudaKernels.h"
#include "GpuFreeEnergyCudaKernels.h"

#include <stdio.h>
#include <cuda.h>
#include <vector_functions.h>
#include <cstdlib>

struct Atom {
    float x;
    float y;
    float z;
    float q;
    float sig;
    float eps;
    float br;
    float softCoreLJLambda;
    float fx;
    float fy;
    float fz;
    float fb;
};

struct cudaFreeEnergySimulation {
    float* pParticleSoftCoreLJLambda;
};

static __constant__ cudaGmxSimulation cSim;
static __constant__ cudaFreeEnergySimulation feSimDev;

void SetCalculateCDLJObcGbsaSoftcoreGpu1Sim( gpuContext gpu )
{
    cudaError_t status;

    //(void) fprintf( stderr, "SetCalculateCDLJObcGbsaSoftcoreGpu1Sim gpu=%p cSim=%p sizeof=%u\n", gpu, &gpu->sim, sizeof(cudaGmxSimulation) ); fflush( stderr );
    status = cudaMemcpyToSymbol(cSim, &gpu->sim, sizeof(cudaGmxSimulation));     
    RTERROR(status, "cudaMemcpyToSymbol: SetCalculateCDLJObcGbsaSoftcoreGpu1Sim copy to cSim failed");

}

void SetCalculateCDLJObcGbsaSoftcoreSupplementary1Sim( float* gpuParticleSoftCoreLJLambda)
{
    cudaError_t status;
    //(void) fprintf( stderr, "SetCalculateCDLJObcGbsaSoftcoreSupplementary1Sim\n" );

    struct cudaFreeEnergySimulation feSim;
    feSim.pParticleSoftCoreLJLambda = gpuParticleSoftCoreLJLambda;
    status = cudaMemcpyToSymbol(feSimDev, &feSim, sizeof(cudaFreeEnergySimulation));
    RTERROR(status, "cudaMemcpyToSymbol: SetCalculateCDLJObcGbsaSoftcoreSupplementary1Sim failed");

}

void GetCalculateCDLJObcGbsaSoftcoreForces1Sim( gpuContext gpu )
{
    cudaError_t status;

    status = cudaMemcpyFromSymbol(&gpu->sim, cSim, sizeof(cudaGmxSimulation));     
    RTERROR(status, "cudaMemcpyFromSymbol: SetSim copy from cSim failed");

}

#if 0
__device__ float fastErfc(float r)
{
    float normalized = cSim.tabulatedErfcScale*r;
    int index = (int) normalized;
    float fract2 = normalized-index;
    float fract1 = 1.0f-fract2;
    return fract1*tex1Dfetch(tabulatedErfcRef, index) + fract2*tex1Dfetch(tabulatedErfcRef, index+1);
}
#endif

// Include versions of the kernel for N^2 calculations.

#if 0
#define METHOD_NAME(a, b) a##N2##b
#include "kCalculateCDLJObcGbsaForces1.h"
#define USE_OUTPUT_BUFFER_PER_WARP
#undef METHOD_NAME
#define METHOD_NAME(a, b) a##N2ByWarp##b
#include "kCalculateCDLJObcGbsaForces1.h"
#endif

// Include versions of the kernel for N^2 calculations with softcore LJ.

#undef METHOD_NAME
#define METHOD_NAME(a, b) a##N2##b
#undef USE_OUTPUT_BUFFER_PER_WARP
#define USE_SOFTCORE_LJ
#include "kCalculateCDLJObcGbsaSoftcoreForces1.h"

#define USE_OUTPUT_BUFFER_PER_WARP
#undef METHOD_NAME
#define METHOD_NAME(a, b) a##N2ByWarp##b
#include "kCalculateCDLJObcGbsaSoftcoreForces1.h"
#undef USE_SOFTCORE_LJ
#undef USE_OUTPUT_BUFFER_PER_WARP

// Include versions of the kernel with cutoffs.

#if 0
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

// Include versions of the kernels for Ewald

#undef METHOD_NAME
#undef USE_OUTPUT_BUFFER_PER_WARP
#define USE_PERIODIC
#define USE_EWALD
#define METHOD_NAME(a, b) a##Ewald##b
#include "kCalculateCDLJObcGbsaForces1.h"
#define USE_OUTPUT_BUFFER_PER_WARP
#undef METHOD_NAME
#define METHOD_NAME(a, b) a##EwaldByWarp##b
#include "kCalculateCDLJObcGbsaForces1.h"

extern __global__ void kFindBlockBoundsCutoff_kernel();
extern __global__ void kFindBlockBoundsPeriodic_kernel();
extern __global__ void kFindBlocksWithInteractionsCutoff_kernel();
extern __global__ void kFindBlocksWithInteractionsPeriodic_kernel();
extern __global__ void kFindInteractionsWithinBlocksCutoff_kernel(unsigned int*);
extern __global__ void kFindInteractionsWithinBlocksPeriodic_kernel(unsigned int*);
extern __global__ void kCalculateEwaldFastCosSinSums_kernel();
extern __global__ void kCalculateEwaldFastForces_kernel();
extern void kCalculatePME(gpuContext gpu);
#endif

/**
 * 
 * Calculate Born radii and first GBSA loop forces/energy
 *
 * @param gpu     gpu contexct
 * @param gbsaObc if set, calculate Born radii for OBC
 *                otherwise calculate Born radii for GB/VI
 *
 */
void kCalculateCDLJObcGbsaSoftcoreForces1(gpuContext gpu )
{
//    printf("kCalculateCDLJObcGbsaForces1\n");
    switch (gpu->sim.nonbondedMethod)
    {
        case NO_CUTOFF:

            // use softcore LJ potential

            if (gpu->bOutputBufferPerWarp)
                   kCalculateCDLJObcGbsaSoftcoreN2ByWarpForces1_kernel<<<gpu->sim.nonbond_blocks, gpu->sim.nonbond_threads_per_block,
                           sizeof(Atom)*gpu->sim.nonbond_threads_per_block>>>(gpu->sim.pWorkUnit);
            else
                   kCalculateCDLJObcGbsaSoftcoreN2Forces1_kernel<<<gpu->sim.nonbond_blocks, gpu->sim.nonbond_threads_per_block,
                           sizeof(Atom)*gpu->sim.nonbond_threads_per_block>>>(gpu->sim.pWorkUnit);
   
            LAUNCHERROR("kCalculateCDLJObcGbsaSoftcoreForces1");
            break;
#if 0
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
        case EWALD:
        case PARTICLE_MESH_EWALD:
            kFindBlockBoundsPeriodic_kernel<<<(gpu->psGridBoundingBox->_length+63)/64, 64>>>();
            LAUNCHERROR("kFindBlockBoundsPeriodic");
            kFindBlocksWithInteractionsPeriodic_kernel<<<gpu->sim.interaction_blocks, gpu->sim.interaction_threads_per_block>>>();
            LAUNCHERROR("kFindBlocksWithInteractionsPeriodic");
            compactStream(gpu->compactPlan, gpu->sim.pInteractingWorkUnit, gpu->sim.pWorkUnit, gpu->sim.pInteractionFlag, gpu->sim.workUnits, gpu->sim.pInteractionCount);
            kFindInteractionsWithinBlocksPeriodic_kernel<<<gpu->sim.nonbond_blocks, gpu->sim.nonbond_threads_per_block,
                    sizeof(unsigned int)*gpu->sim.nonbond_threads_per_block>>>(gpu->sim.pInteractingWorkUnit);
            LAUNCHERROR("kFindInteractionsWithinBlocksPeriodic");
            if (gpu->bRecalculateBornRadii)
            {
                kCalculateObcGbsaBornSum(gpu);
                kReduceObcGbsaBornSum(gpu);
            }
            cudaChannelFormatDesc channelDesc = cudaCreateChannelDesc<float>();
            cudaBindTexture(NULL, &tabulatedErfcRef, gpu->psTabulatedErfc->_pDevData, &channelDesc, gpu->psTabulatedErfc->_length*sizeof(float));
            if (gpu->bOutputBufferPerWarp)
                kCalculateCDLJObcGbsaEwaldByWarpForces1_kernel<<<gpu->sim.nonbond_blocks, gpu->sim.nonbond_threads_per_block,
                        (sizeof(Atom)+sizeof(float3))*gpu->sim.nonbond_threads_per_block>>>(gpu->sim.pInteractingWorkUnit);
            else
                kCalculateCDLJObcGbsaEwaldForces1_kernel<<<gpu->sim.nonbond_blocks, gpu->sim.nonbond_threads_per_block,
                        (sizeof(Atom)+sizeof(float3))*gpu->sim.nonbond_threads_per_block>>>(gpu->sim.pInteractingWorkUnit);
            LAUNCHERROR("kCalculateCDLJObcGbsaEwaldForces");
            if (gpu->sim.nonbondedMethod == EWALD)
            {
                // Ewald summation
                kCalculateEwaldFastCosSinSums_kernel<<<gpu->sim.nonbond_blocks, gpu->sim.nonbond_threads_per_block>>>();
                LAUNCHERROR("kCalculateEwaldFastCosSinSums");
                kCalculateEwaldFastForces_kernel<<<gpu->sim.blocks, gpu->sim.update_threads_per_block>>>();
                LAUNCHERROR("kCalculateEwaldFastForces");
            }
            else
                kCalculatePME(gpu);
#endif
    }
}
