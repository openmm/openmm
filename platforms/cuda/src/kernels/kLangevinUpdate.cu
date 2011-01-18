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

enum {VelScale, ForceScale, NoiseScale, MaxParams};

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


__global__ 
#if (__CUDA_ARCH__ >= 200)
__launch_bounds__(GF1XX_UPDATE_THREADS_PER_BLOCK, 1)
#elif (__CUDA_ARCH__ >= 120)
__launch_bounds__(GT2XX_UPDATE_THREADS_PER_BLOCK, 1)
#else
__launch_bounds__(G8X_UPDATE_THREADS_PER_BLOCK, 1)
#endif
void kSelectLangevinStepSize_kernel(float maxStepSize)
{
    // Calculate the error.

    extern __shared__ float error[];
    __shared__ float params[MaxParams];
    error[threadIdx.x] = 0.0f;
    unsigned int pos = threadIdx.x;
    while (pos < cSim.atoms)
    {
        float4 force  = cSim.pForce4[pos];
        float invMass = cSim.pVelm4[pos].w;
        error[threadIdx.x] += (force.x*force.x + force.y*force.y + force.z*force.z)*invMass;
        pos += blockDim.x * gridDim.x;
    }
    __syncthreads();

    // Sum the errors from all threads.

    for (int offset = 1; offset < blockDim.x; offset *= 2)
    {
        if (threadIdx.x+offset < blockDim.x && (threadIdx.x&(2*offset-1)) == 0)
            error[threadIdx.x] += error[threadIdx.x+offset];
        __syncthreads();
    }
    if (threadIdx.x == 0)
    {
        // Select the new step size.
        
        float totalError = sqrt(error[0]/(cSim.atoms*3));
        float newStepSize = sqrt(cSim.errorTol/totalError);
        float oldStepSize = cSim.pStepSize[0].y;
        if (oldStepSize > 0.0f)
            newStepSize = min(newStepSize, oldStepSize*2.0f); // For safety, limit how quickly dt can increase.
        if (newStepSize > oldStepSize && newStepSize < 1.1f*oldStepSize)
            newStepSize = oldStepSize; // Keeping dt constant between steps improves the behavior of the integrator.
        if (newStepSize > maxStepSize)
            newStepSize = maxStepSize;
        cSim.pStepSize[0].y = newStepSize;

        // Recalculate the integration parameters.

        float vscale = exp(-newStepSize/cSim.tau);
        float fscale = (1-vscale)*cSim.tau;
        float noisescale = sqrt(2*cSim.kT/cSim.tau)*sqrt(0.5f*(1-vscale*vscale)*cSim.tau);
        params[VelScale] = vscale;
        params[ForceScale] = fscale;
        params[NoiseScale] = noisescale;
    }
    __syncthreads();
    if (threadIdx.x < MaxParams)
        cSim.pLangevinParameters[threadIdx.x] = params[threadIdx.x];
}

void kSelectLangevinStepSize(gpuContext gpu, float maxTimeStep)
{
//    printf("kSelectLangevinStepSize\n");
    kSelectLangevinStepSize_kernel<<<1, gpu->sim.update_threads_per_block, sizeof(float)*gpu->sim.update_threads_per_block>>>(maxTimeStep);
    LAUNCHERROR("kSelectLangevinStepSize");
}

__global__
#if (__CUDA_ARCH__ >= 200)
__launch_bounds__(GF1XX_UPDATE_THREADS_PER_BLOCK, 1)
#elif (__CUDA_ARCH__ >= 120)
__launch_bounds__(GT2XX_UPDATE_THREADS_PER_BLOCK, 1)
#else
__launch_bounds__(G8X_UPDATE_THREADS_PER_BLOCK, 1)
#endif
void kSetVelocitiesFromPositions_kernel()
{
    float2 stepSize = cSim.pStepSize[0];
    double oneOverDt = 2.0/(stepSize.x+stepSize.y);
    unsigned int pos = threadIdx.x;
    while (pos < cSim.atoms)
    {
        float4 posq = cSim.pPosq[pos];
        float4 oldPosq = cSim.pOldPosq[pos];
        float4 velm = cSim.pVelm4[pos];
        velm.x = (float) (oneOverDt*((double)posq.x-(double)oldPosq.x));
        velm.y = (float) (oneOverDt*((double)posq.y-(double)oldPosq.y));
        velm.z = (float) (oneOverDt*((double)posq.z-(double)oldPosq.z));
        cSim.pVelm4[pos] = velm;
        pos += blockDim.x * gridDim.x;
    }
}

void kSetVelocitiesFromPositions(gpuContext gpu)
{
//    printf("kSetVelocitiesFromPositions\n");
    kSetVelocitiesFromPositions_kernel<<<1, gpu->sim.update_threads_per_block, sizeof(float)*gpu->sim.update_threads_per_block>>>();
    LAUNCHERROR("kSetVelocitiesFromPositions");
}
