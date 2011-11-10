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
#include "openmm/OpenMMException.h"

#include <stdio.h>
#include <cuda.h>
#include <vector_functions.h>
#include <cstdlib>
#include <iostream>
#include <sstream>

#define USE_SOFTCORE_LJ

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

static __constant__ cudaGmxSimulation cSim;
static __constant__ cudaFreeEnergyGmxSimulation feSimDev;

void SetCalculateCDLJObcGbsaSoftcoreGpu1Sim( freeEnergyGpuContext freeEnergyGpu ){

    cudaError_t status;
    status = cudaMemcpyToSymbol(cSim, &freeEnergyGpu->gpuContext->sim, sizeof(cudaGmxSimulation));
    RTERROR(status, "cudaMemcpyToSymbol: SetCalculateCDLJObcGbsaSoftcoreGpu1Sim copy to cSim failed");

    status = cudaMemcpyToSymbol( feSimDev, &freeEnergyGpu->freeEnergySim, sizeof(cudaFreeEnergyGmxSimulation));
    RTERROR(status, "cudaMemcpyToSymbol: SetCalculateCDLJObcGbsaSoftcoreGpu1Sim copy to feSimDev failed");
}

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

#undef METHOD_NAME
#undef USE_OUTPUT_BUFFER_PER_WARP
#define USE_CUTOFF
#define METHOD_NAME(a, b) a##Cutoff##b
#include "kCalculateCDLJObcGbsaSoftcoreForces1.h"
#define USE_OUTPUT_BUFFER_PER_WARP
#undef METHOD_NAME
#define METHOD_NAME(a, b) a##CutoffByWarp##b
#include "kCalculateCDLJObcGbsaSoftcoreForces1.h"

// Include versions of the kernel with periodic boundary conditions.

#undef METHOD_NAME
#undef USE_OUTPUT_BUFFER_PER_WARP
#define USE_PERIODIC
#define METHOD_NAME(a, b) a##Periodic##b
#include "kCalculateCDLJObcGbsaSoftcoreForces1.h"
#define USE_OUTPUT_BUFFER_PER_WARP
#undef METHOD_NAME
#define METHOD_NAME(a, b) a##PeriodicByWarp##b
#include "kCalculateCDLJObcGbsaSoftcoreForces1.h"


/**
 * 
 * Calculate Born radii and first GBSA loop forces/energy
 *
 * @param gpu     gpu context
 *
 */
void kCalculateCDLJObcGbsaSoftcoreForces1( freeEnergyGpuContext freeEnergyGpu )
{
    unsigned int threadsPerBlock;
    static unsigned int threadsPerBlockPerMethod[3] = { 0, 0, 0 };
    static unsigned int natoms[3]                   = { 0, 0, 0 };
    gpuContext gpu                                  = freeEnergyGpu->gpuContext;
    unsigned int methodIndex                        = static_cast<unsigned int>(freeEnergyGpu->freeEnergySim.nonbondedMethod);

    if( methodIndex > 2 ){
        throw OpenMM::OpenMMException( "kCalculateCDLJObcGbsaSoftcoreForces1 method index invalid." );
    }   

    if( natoms[methodIndex] != gpu->natoms ){
        unsigned int extra                    = methodIndex == 0 ? 0 : sizeof(float);
        threadsPerBlockPerMethod[methodIndex] = std::min(getThreadsPerBlockFEP( freeEnergyGpu, (sizeof(Atom) + extra), gpu->sharedMemoryPerBlock ), gpu->sim.nonbond_threads_per_block );
        natoms[methodIndex]                   = gpu->natoms;
    }   
    threadsPerBlock = threadsPerBlockPerMethod[methodIndex];

    switch( freeEnergyGpu->freeEnergySim.nonbondedMethod )
    {
        case FREE_ENERGY_NO_CUTOFF:

            // use softcore LJ potential
            if (gpu->bOutputBufferPerWarp)
                   kCalculateCDLJObcGbsaSoftcoreN2ByWarpForces1_kernel<<<gpu->sim.nonbond_blocks, threadsPerBlock,
                           sizeof(Atom)*threadsPerBlock>>>(gpu->sim.pWorkUnit );
            else
                   kCalculateCDLJObcGbsaSoftcoreN2Forces1_kernel<<<gpu->sim.nonbond_blocks, threadsPerBlock,
                           sizeof(Atom)*threadsPerBlock>>>(gpu->sim.pWorkUnit );
   
            LAUNCHERROR("kCalculateCDLJObcGbsaSoftcoreForces1");
            break;

        case FREE_ENERGY_CUTOFF:

            if (gpu->bOutputBufferPerWarp)
                kCalculateCDLJObcGbsaSoftcoreCutoffByWarpForces1_kernel<<<gpu->sim.nonbond_blocks, threadsPerBlock,
                        (sizeof(Atom)+sizeof(float))*threadsPerBlock>>>(gpu->sim.pInteractingWorkUnit );
            else
                kCalculateCDLJObcGbsaSoftcoreCutoffForces1_kernel<<<gpu->sim.nonbond_blocks, threadsPerBlock,
                        (sizeof(Atom)+sizeof(float))*threadsPerBlock>>>(gpu->sim.pInteractingWorkUnit );

            LAUNCHERROR("kCalculateCDLJObcGbsaSoftcoreCutoffForces1");

            break;

        case FREE_ENERGY_PERIODIC:

            if (gpu->bOutputBufferPerWarp)
                kCalculateCDLJObcGbsaSoftcorePeriodicByWarpForces1_kernel<<<gpu->sim.nonbond_blocks, threadsPerBlock,
                        (sizeof(Atom)+sizeof(float))*threadsPerBlock>>>(gpu->sim.pInteractingWorkUnit);
            else
                kCalculateCDLJObcGbsaSoftcorePeriodicForces1_kernel<<<gpu->sim.nonbond_blocks, threadsPerBlock,
                        (sizeof(Atom)+sizeof(float))*threadsPerBlock>>>(gpu->sim.pInteractingWorkUnit);

            LAUNCHERROR("kCalculateCDLJObcGbsaSoftcorePeriodicForces1");
            break;
    }

}
