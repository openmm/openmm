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
//#include <fstream>
using namespace std;

#include "gputypes.h"

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
        apos.x                  = force.x*cSim.tauDeltaT + cSim.noiseAmplitude*random4a.x;
        apos.y                  = force.y*cSim.tauDeltaT + cSim.noiseAmplitude*random4a.y;
        apos.z                  = force.z*cSim.tauDeltaT + cSim.noiseAmplitude*random4a.z;
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

        velocity.x              = cSim.oneOverDeltaT*(xPrime.x);
        velocity.y              = cSim.oneOverDeltaT*(xPrime.y);
        velocity.z              = cSim.oneOverDeltaT*(xPrime.z);

        xPrime.x               += apos.x;
        xPrime.y               += apos.y;
        xPrime.z               += apos.z;

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

