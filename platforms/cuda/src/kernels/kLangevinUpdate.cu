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
using namespace std;

#include "gputypes.h"

static __constant__ cudaGmxSimulation cSim;

void SetLangevinUpdateSim(gpuContext gpu)
{
    cudaError_t status;
    status = cudaMemcpyToSymbol(cSim, &gpu->sim, sizeof(cudaGmxSimulation));     
    RTERROR(status, "cudaMemcpyToSymbol: SetSim copy to cSim failed");
}

void GetLangevinUpdateSim(gpuContext gpu)
{
    cudaError_t status;
    status = cudaMemcpyFromSymbol(&gpu->sim, cSim, sizeof(cudaGmxSimulation));     
    RTERROR(status, "cudaMemcpyFromSymbol: SetSim copy from cSim failed");
}

// Include versions of the kernels with and with center of mass motion removal.

#include "kLangevinUpdate.h"
#define REMOVE_CM
#include "kLangevinUpdate.h"

void kLangevinUpdatePart1(gpuContext gpu)
{
//    printf("kLangevinUpdatePart1\n");
    if (gpu->bRemoveCM)
    {
        kLangevinUpdatePart1CM_kernel<<<gpu->sim.blocks, gpu->sim.update_threads_per_block, gpu->sim.update_threads_per_block * sizeof(float3)>>>();
        LAUNCHERROR("kLangevinUpdatePart1CM");
        gpu->bRemoveCM = false;
    }
    else
    {    
        kLangevinUpdatePart1_kernel<<<gpu->sim.blocks, gpu->sim.update_threads_per_block>>>();
        LAUNCHERROR("kLangevinUpdatePart1");
    }
}

extern void kGenerateRandoms(gpuContext gpu);
void kLangevinUpdatePart2(gpuContext gpu)
{
//    printf("kLangevinUpdatePart2\n");
    if (gpu->bCalculateCM)
    {
        kLangevinUpdatePart2CM_kernel<<<gpu->sim.blocks, gpu->sim.update_threads_per_block, gpu->sim.update_threads_per_block * sizeof(float3)>>>();
        LAUNCHERROR("kLangevinUpdatePart2CM");
        gpu->bCalculateCM = false;
        gpu->bRemoveCM = true;
    }
    else
    {
        kLangevinUpdatePart2_kernel<<<gpu->sim.blocks, gpu->sim.update_threads_per_block>>>();
        LAUNCHERROR("kLangevinUpdatePart2");
    }
    
    // Update randoms if necessary
    gpu->iterations++;
    if (gpu->iterations == gpu->sim.randomIterations)
    {
        kGenerateRandoms(gpu);
        gpu->iterations = 0;
    }
}

