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

enum {EM, EM_V, DOverTauC, TauOneMinusEM_V, TauDOverEMMinusOne, V, X, Yv, Yx, Fix1, OneOverFix1, MaxParams};

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


__global__ void kSelectLangevinStepSize_kernel(float maxStepSize)
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

        float gdt                  = newStepSize / cSim.tau;
        float eph                  = exp(0.5f * gdt);
        float emh                  = exp(-0.5f * gdt);
        float ep                   = exp(gdt);
        float em                   = exp(-gdt);
        float em_v                 = exp(-0.5f*(oldStepSize+newStepSize)/cSim.tau);
        float b, c, d;
        if (gdt >= 0.1f)
        {
            float term1 = eph - 1.0f;
            term1                 *= term1;
            b                      = gdt * (ep - 1.0f) - 4.0f * term1;
            c                      = gdt - 3.0f + 4.0f * emh - em;
            d                      = 2.0f - eph - emh;
        }
        else
        {
            float term1            = 0.5f * gdt;
            float term2            = term1 * term1;
            float term4            = term2 * term2;

            float third            = 1.0f / 3.0f;
            float o7_9             = 7.0f / 9.0f;
            float o1_12            = 1.0f / 12.0f;
            float o17_90           = 17.0f / 90.0f;
            float o7_30            = 7.0f / 30.0f;
            float o31_1260         = 31.0f / 1260.0f;
            float o_360            = 1.0f / 360.0f;

            b                      = term4 * (third + term1 * (third + term1 * (o17_90 + term1 * o7_9)));
            c                      = term2 * term1 * (2.0f * third + term1 * (-0.5f + term1 * (o7_30 + term1 * (-o1_12 + term1 * o31_1260))));
            d                      = term2 * (-1.0f + term2 * (-o1_12 - term2 * o_360));
        }
        float fix1                 = cSim.tau * (eph - emh);
        if (fix1 == 0.0f)
            fix1 = newStepSize;
        params[EM]                 = em;
        params[EM_V]               = em_v;
        params[DOverTauC]          = d / (cSim.tau * c);
        params[TauOneMinusEM_V]    = cSim.tau * (1.0f-em_v);
        params[TauDOverEMMinusOne] = cSim.tau * d / (em - 1.0f);
        params[Fix1]               = fix1;
        params[OneOverFix1]        = 1.0f / fix1;
        params[V]                  = sqrt(cSim.kT * (1.0f - em));
        params[X]                  = cSim.tau * sqrt(cSim.kT * c);
        params[Yv]                 = sqrt(cSim.kT * b / c);
        params[Yx]                 = cSim.tau * sqrt(cSim.kT * b / (1.0f - em));
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
