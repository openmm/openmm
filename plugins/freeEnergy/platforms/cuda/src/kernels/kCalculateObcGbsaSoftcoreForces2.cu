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
    float fx;
    float fy;
    float fz;
    float fb;
};


static __constant__ cudaGmxSimulation cSim;

extern "C"
void SetCalculateObcGbsaSoftcoreForces2Sim(gpuContext gpu)
{
    cudaError_t status;
    status = cudaMemcpyToSymbol(cSim, &gpu->sim, sizeof(cudaGmxSimulation));     
    RTERROR(status, "cudaMemcpyToSymbol: SetSim copy to cSim failed");
}

void GetCalculateObcGbsaSoftcoreForces2Sim(gpuContext gpu)
{
    cudaError_t status;
    status = cudaMemcpyFromSymbol(&gpu->sim, cSim, sizeof(cudaGmxSimulation));     
    RTERROR(status, "cudaMemcpyFromSymbol: SetSim copy from cSim failed");
}

// Include versions of the kernels for N^2 calculations.

#define METHOD_NAME(a, b) a##N2##b
#include "kCalculateObcGbsaSoftcoreForces2.h"
#define USE_OUTPUT_BUFFER_PER_WARP
#undef METHOD_NAME
#define METHOD_NAME(a, b) a##N2ByWarp##b
#include "kCalculateObcGbsaSoftcoreForces2.h"

// Include versions of the kernels with cutoffs.

#undef METHOD_NAME
#undef USE_OUTPUT_BUFFER_PER_WARP
#define USE_CUTOFF
#define METHOD_NAME(a, b) a##Cutoff##b
#include "kCalculateObcGbsaSoftcoreForces2.h"
#define USE_OUTPUT_BUFFER_PER_WARP
#undef METHOD_NAME
#define METHOD_NAME(a, b) a##CutoffByWarp##b
#include "kCalculateObcGbsaSoftcoreForces2.h"

// Include versions of the kernels with periodic boundary conditions.

#undef METHOD_NAME
#undef USE_OUTPUT_BUFFER_PER_WARP
#define USE_PERIODIC
#define METHOD_NAME(a, b) a##Periodic##b
#include "kCalculateObcGbsaSoftcoreForces2.h"
#define USE_OUTPUT_BUFFER_PER_WARP
#undef METHOD_NAME
#define METHOD_NAME(a, b) a##PeriodicByWarp##b
#include "kCalculateObcGbsaSoftcoreForces2.h"

void kCalculateObcGbsaSoftcoreForces2(gpuContext gpu)
{
    //printf("kCalculateObcGbsaSoftcoreForces2\n");
    //fprintf( stderr, "kCalculateObcGbsaSoftcoreForces2 nonbondedMethod=%d warp=%d\n", gpu->sim.nonbondedMethod, gpu->bOutputBufferPerWarp);
//fprintf( stderr, "kCalculateObcGbsaSoftcoreForces2 nonbondedMethod=%d calling kReduceForces\n", gpu->sim.nonbondedMethod);
//kReduceForces(gpu);
    switch (gpu->sim.nonbondedMethod)
    {
        case NO_CUTOFF:
            if (gpu->bOutputBufferPerWarp)
                kCalculateObcGbsaSoftcoreN2ByWarpForces2_kernel<<<gpu->sim.bornForce2_blocks, gpu->sim.bornForce2_threads_per_block,
                        sizeof(Atom)*gpu->sim.bornForce2_threads_per_block>>>(gpu->sim.pWorkUnit);
            else
                kCalculateObcGbsaSoftcoreN2Forces2_kernel<<<gpu->sim.bornForce2_blocks, gpu->sim.bornForce2_threads_per_block,
                        sizeof(Atom)*gpu->sim.bornForce2_threads_per_block>>>(gpu->sim.pWorkUnit);
            break;
        case CUTOFF:
            if (gpu->bOutputBufferPerWarp)
                kCalculateObcGbsaSoftcoreCutoffByWarpForces2_kernel<<<gpu->sim.bornForce2_blocks, gpu->sim.bornForce2_threads_per_block,
                        (sizeof(Atom)+sizeof(float3))*gpu->sim.bornForce2_threads_per_block>>>(gpu->sim.pInteractingWorkUnit);
            else
                kCalculateObcGbsaSoftcoreCutoffForces2_kernel<<<gpu->sim.bornForce2_blocks, gpu->sim.bornForce2_threads_per_block,
                        (sizeof(Atom)+sizeof(float3))*gpu->sim.bornForce2_threads_per_block>>>(gpu->sim.pInteractingWorkUnit);
            break;
        case PERIODIC:
            if (gpu->bOutputBufferPerWarp)
                kCalculateObcGbsaSoftcorePeriodicByWarpForces2_kernel<<<gpu->sim.bornForce2_blocks, gpu->sim.bornForce2_threads_per_block,
                        (sizeof(Atom)+sizeof(float3))*gpu->sim.bornForce2_threads_per_block>>>(gpu->sim.pInteractingWorkUnit);
            else
                kCalculateObcGbsaSoftcorePeriodicForces2_kernel<<<gpu->sim.bornForce2_blocks, gpu->sim.bornForce2_threads_per_block,
                        (sizeof(Atom)+sizeof(float3))*gpu->sim.bornForce2_threads_per_block>>>(gpu->sim.pInteractingWorkUnit);
            break;
    }
    LAUNCHERROR("kCalculateObcGbsaSoftcoreForces2");
}
