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
