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
//#include <fstream>
using namespace std;

#include "gputypes.h"

#define DeltaShake

static __constant__ cudaGmxSimulation cSim;

void SetBrownianUpdateSim(gpuContext gpu)
{
    cudaError_t status;
    status = cudaMemcpyToSymbol(cSim, &gpu->sim, sizeof(cudaGmxSimulation));     
    RTERROR(status, "cudaMemcpyToSymbol: SetSim copy to cSim failed");
}

void GetBrownianUpdateSim(gpuContext gpu)
{
    cudaError_t status;
    status = cudaMemcpyFromSymbol(&gpu->sim, cSim, sizeof(cudaGmxSimulation));     
    RTERROR(status, "cudaMemcpyFromSymbol: SetSim copy from cSim failed");
}

__global__ void kBrownianUpdatePart1_kernel()
{
    unsigned int pos    = threadIdx.x + blockIdx.x * blockDim.x;
    unsigned int rpos   = cSim.pRandomPosition[blockIdx.x];
    __syncthreads();
    
    while (pos < cSim.atoms)
    {
        float4 random4a         = cSim.pRandom4a[rpos + pos];
        float4 apos             = cSim.pPosq[pos];
        float4 force            = cSim.pForce4[pos];

        cSim.pOldPosq[pos]      = apos;
#ifndef DeltaShake
        apos.x                 += force.x*cSim.GDT + random4a.x;
        apos.y                 += force.y*cSim.GDT + random4a.y;
        apos.z                 += force.z*cSim.GDT + random4a.z;
#else
        apos.x                  = force.x*cSim.GDT + random4a.x;
        apos.y                  = force.y*cSim.GDT + random4a.y;
        apos.z                  = force.z*cSim.GDT + random4a.z;
#endif
        cSim.pPosqP[pos]        = apos;
        pos                    += blockDim.x * gridDim.x;
    }
}

void kBrownianUpdatePart1(gpuContext gpu)
{
//    printf("kBrownianUpdatePart1\n");
    kBrownianUpdatePart1_kernel<<<gpu->sim.blocks, gpu->sim.update_threads_per_block>>>();
    LAUNCHERROR("kBrownianUpdatePart1");
}

__global__ void kBrownianUpdatePart2_kernel()
{
    unsigned int pos            = threadIdx.x + blockIdx.x * blockDim.x;
    unsigned int rpos           = cSim.pRandomPosition[blockIdx.x];
    __syncthreads();
    
    while (pos < cSim.atoms)
    {
        float4 velocity         = cSim.pVelm4[pos];
        float4 apos             = cSim.pPosq[pos];
        float4 xPrime           = cSim.pPosqP[pos];

#ifndef DeltaShake
        velocity.x              = cSim.oneOverDeltaT*(xPrime.x-apos.x);
        velocity.y              = cSim.oneOverDeltaT*(xPrime.y-apos.y);
        velocity.z              = cSim.oneOverDeltaT*(xPrime.z-apos.z);
#else
        velocity.x              = cSim.oneOverDeltaT*(xPrime.x);
        velocity.y              = cSim.oneOverDeltaT*(xPrime.y);
        velocity.z              = cSim.oneOverDeltaT*(xPrime.z);

        xPrime.x               += apos.x;
        xPrime.y               += apos.y;
        xPrime.z               += apos.z;
#endif
        cSim.pPosq[pos]         = xPrime;
        cSim.pVelm4[pos]        = velocity;
         
        pos                    += blockDim.x * gridDim.x;    
    }

    // Update random position pointer
    if (threadIdx.x == 0)
    {
        rpos                   += cSim.paddedNumberOfAtoms;
        if (rpos > cSim.randoms)
            rpos               -= cSim.randoms;
        cSim.pRandomPosition[blockIdx.x] = rpos;
    }
}

extern void kGenerateRandoms(gpuContext gpu);
void kBrownianUpdatePart2(gpuContext gpu)
{
//    printf("kBrownianUpdatePart2\n");
    kBrownianUpdatePart2_kernel<<<gpu->sim.blocks, gpu->sim.update_threads_per_block>>>();
    LAUNCHERROR("kBrownianUpdatePart2");
    
    // Update randoms if necessary
    gpu->iterations++;
    if (gpu->iterations == gpu->sim.randomIterations)
    {
        kGenerateRandoms(gpu);
        gpu->iterations = 0;
    }
}

