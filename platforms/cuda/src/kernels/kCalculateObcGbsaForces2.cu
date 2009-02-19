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
using namespace std;

#include "gputypes.h"
#include "cudaKernels.h"

struct Atom {
    float x;
    float y;
    float z;
    float r;
    float sr;
    float sr2;
    float fx;
    float fy;
    float fz;
    float fb;
};


static __constant__ cudaGmxSimulation cSim;

void SetCalculateObcGbsaForces2Sim(gpuContext gpu)
{
    cudaError_t status;
    status = cudaMemcpyToSymbol(cSim, &gpu->sim, sizeof(cudaGmxSimulation));     
    RTERROR(status, "cudaMemcpyToSymbol: SetSim copy to cSim failed");
}

void GetCalculateObcGbsaForces2Sim(gpuContext gpu)
{
    cudaError_t status;
    status = cudaMemcpyFromSymbol(&gpu->sim, cSim, sizeof(cudaGmxSimulation));     
    RTERROR(status, "cudaMemcpyFromSymbol: SetSim copy from cSim failed");
}


// Include versions of the kernels for N^2 calculations.

#define METHOD_NAME(a, b) a##N2##b
#include "kCalculateObcGbsaForces2.h"
#define USE_OUTPUT_BUFFER_PER_WARP
#undef METHOD_NAME
#define METHOD_NAME(a, b) a##N2ByWarp##b
#include "kCalculateObcGbsaForces2.h"

// Include versions of the kernels with cutoffs.

#undef METHOD_NAME
#undef USE_OUTPUT_BUFFER_PER_WARP
#define USE_CUTOFF
#define METHOD_NAME(a, b) a##Cutoff##b
#include "kCalculateObcGbsaForces2.h"
#define USE_OUTPUT_BUFFER_PER_WARP
#undef METHOD_NAME
#define METHOD_NAME(a, b) a##CutoffByWarp##b
#include "kCalculateObcGbsaForces2.h"

// Include versions of the kernels with periodic boundary conditions.

#undef METHOD_NAME
#undef USE_OUTPUT_BUFFER_PER_WARP
#define USE_PERIODIC
#define METHOD_NAME(a, b) a##Periodic##b
#include "kCalculateObcGbsaForces2.h"
#define USE_OUTPUT_BUFFER_PER_WARP
#undef METHOD_NAME
#define METHOD_NAME(a, b) a##PeriodicByWarp##b
#include "kCalculateObcGbsaForces2.h"

void kCalculateObcGbsaForces2(gpuContext gpu)
{
    //printf("kCalculateObcGbsaForces2\n");
    size_t numWithInteractions;
    switch (gpu->sim.nonbondedMethod)
    {
        case NO_CUTOFF:
            if (gpu->bOutputBufferPerWarp)
                kCalculateObcGbsaN2ByWarpForces2_kernel<<<gpu->sim.bornForce2_blocks, gpu->sim.bornForce2_threads_per_block,
                        sizeof(Atom)*gpu->sim.bornForce2_threads_per_block>>>(gpu->sim.pWorkUnit, gpu->sim.workUnits);
            else
                kCalculateObcGbsaN2Forces2_kernel<<<gpu->sim.bornForce2_blocks, gpu->sim.bornForce2_threads_per_block,
                        sizeof(Atom)*gpu->sim.bornForce2_threads_per_block>>>(gpu->sim.pWorkUnit, gpu->sim.workUnits);
            break;
        case CUTOFF:
            numWithInteractions = gpu->psInteractionCount->_pSysData[0];
            if (gpu->bOutputBufferPerWarp)
                kCalculateObcGbsaCutoffByWarpForces2_kernel<<<gpu->sim.bornForce2_blocks, gpu->sim.bornForce2_threads_per_block,
                        sizeof(Atom)*gpu->sim.bornForce2_threads_per_block>>>(gpu->sim.pInteractingWorkUnit, numWithInteractions);
            else
                kCalculateObcGbsaCutoffForces2_kernel<<<gpu->sim.bornForce2_blocks, gpu->sim.bornForce2_threads_per_block,
                        sizeof(Atom)*gpu->sim.bornForce2_threads_per_block>>>(gpu->sim.pInteractingWorkUnit, numWithInteractions);
            break;
        case PERIODIC:
            numWithInteractions = gpu->psInteractionCount->_pSysData[0];
            if (gpu->bOutputBufferPerWarp)
                kCalculateObcGbsaPeriodicByWarpForces2_kernel<<<gpu->sim.bornForce2_blocks, gpu->sim.bornForce2_threads_per_block,
                        sizeof(Atom)*gpu->sim.bornForce2_threads_per_block>>>(gpu->sim.pInteractingWorkUnit, numWithInteractions);
            else
                kCalculateObcGbsaPeriodicForces2_kernel<<<gpu->sim.bornForce2_blocks, gpu->sim.bornForce2_threads_per_block,
                        sizeof(Atom)*gpu->sim.bornForce2_threads_per_block>>>(gpu->sim.pInteractingWorkUnit, numWithInteractions);
            break;
    }
    LAUNCHERROR("kCalculateObcGbsaForces2");
}
