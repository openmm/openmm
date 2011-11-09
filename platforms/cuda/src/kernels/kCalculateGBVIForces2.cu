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
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
using namespace std;

#include "gputypes.h"
#include "cudaKernels.h"

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

void SetCalculateGBVIForces2Sim(gpuContext gpu)
{
    cudaError_t status;
    status = cudaMemcpyToSymbol(cSim, &gpu->sim, sizeof(cudaGmxSimulation));     
    RTERROR(status, "cudaMemcpyToSymbol: SetSim copy to cSim failed");
}

void GetCalculateGBVIForces2Sim(gpuContext gpu)
{
    cudaError_t status;
    status = cudaMemcpyFromSymbol(&gpu->sim, cSim, sizeof(cudaGmxSimulation));     
    RTERROR(status, "cudaMemcpyFromSymbol: SetSim copy from cSim failed");
}

#include "kCalculateGBVIAux.h"

// Include versions of the kernels for N^2 calculations.

#define METHOD_NAME(a, b) a##N2##b
#include "kCalculateGBVIForces2.h"
#define USE_OUTPUT_BUFFER_PER_WARP
#undef METHOD_NAME
#define METHOD_NAME(a, b) a##N2ByWarp##b
#include "kCalculateGBVIForces2.h"

// Include versions of the kernels with cutoffs.

#undef METHOD_NAME
#undef USE_OUTPUT_BUFFER_PER_WARP
#define USE_CUTOFF
#define METHOD_NAME(a, b) a##Cutoff##b
#include "kCalculateGBVIForces2.h"
#define USE_OUTPUT_BUFFER_PER_WARP
#undef METHOD_NAME
#define METHOD_NAME(a, b) a##CutoffByWarp##b
#include "kCalculateGBVIForces2.h"

// Include versions of the kernels with periodic boundary conditions.

#undef METHOD_NAME
#undef USE_OUTPUT_BUFFER_PER_WARP
#define USE_PERIODIC
#define METHOD_NAME(a, b) a##Periodic##b
#include "kCalculateGBVIForces2.h"
#define USE_OUTPUT_BUFFER_PER_WARP
#undef METHOD_NAME
#define METHOD_NAME(a, b) a##PeriodicByWarp##b
#include "kCalculateGBVIForces2.h"

void kCalculateGBVIForces2(gpuContext gpu)
{

    switch (gpu->sim.nonbondedMethod)
    {
        case NO_CUTOFF:

            if (gpu->bOutputBufferPerWarp)
                kCalculateGBVIN2ByWarpForces2_kernel<<<gpu->sim.bornForce2_blocks, gpu->sim.bornForce2_threads_per_block,
                        sizeof(Atom)*gpu->sim.bornForce2_threads_per_block>>>(gpu->sim.pWorkUnit );
            else
                kCalculateGBVIN2Forces2_kernel<<<gpu->sim.bornForce2_blocks, gpu->sim.bornForce2_threads_per_block,
                        sizeof(Atom)*gpu->sim.bornForce2_threads_per_block>>>(gpu->sim.pWorkUnit );
            break;

        case CUTOFF:

            if (gpu->bOutputBufferPerWarp)
                kCalculateGBVICutoffByWarpForces2_kernel<<<gpu->sim.bornForce2_blocks, gpu->sim.bornForce2_threads_per_block,
                        (sizeof(Atom)+sizeof(float3))*gpu->sim.bornForce2_threads_per_block>>>(gpu->sim.pInteractingWorkUnit );
            else
                kCalculateGBVICutoffForces2_kernel<<<gpu->sim.bornForce2_blocks, gpu->sim.bornForce2_threads_per_block,
                        (sizeof(Atom)+sizeof(float3))*gpu->sim.bornForce2_threads_per_block>>>(gpu->sim.pInteractingWorkUnit);
            break;

        case PERIODIC:

            if (gpu->bOutputBufferPerWarp)
                kCalculateGBVIPeriodicByWarpForces2_kernel<<<gpu->sim.bornForce2_blocks, gpu->sim.bornForce2_threads_per_block,
                        (sizeof(Atom)+sizeof(float3))*gpu->sim.bornForce2_threads_per_block>>>(gpu->sim.pInteractingWorkUnit );
            else
                kCalculateGBVIPeriodicForces2_kernel<<<gpu->sim.bornForce2_blocks, gpu->sim.bornForce2_threads_per_block,
                        (sizeof(Atom)+sizeof(float3))*gpu->sim.bornForce2_threads_per_block>>>(gpu->sim.pInteractingWorkUnit );
            break;

    }
    LAUNCHERROR("kCalculateGBVIForces2");

    //kPrintGBVI( gpu, "kCalculateGBVIForces2", 0, stderr);
}
