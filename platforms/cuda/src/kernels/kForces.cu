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

#define FABS(a) ((a) > 0.0f ? (a) : -(a))

static __constant__ cudaGmxSimulation cSim;

void OPENMMCUDA_EXPORT SetForcesSim(gpuContext gpu)
{
    cudaError_t status;
    status = cudaMemcpyToSymbol(cSim, &gpu->sim, sizeof(cudaGmxSimulation));     
    RTERROR(status, "cudaMemcpyToSymbol: SetForcesSim copy to cSim failed");
}

void GetForcesSim(gpuContext gpu)
{
    cudaError_t status;
    status = cudaMemcpyFromSymbol(&gpu->sim, cSim, sizeof(cudaGmxSimulation));     
    RTERROR(status, "cudaMemcpyFromSymbol: GetForcesSim copy from cSim failed");
}

__global__ 
__launch_bounds__(384, 1)
void kClearForces_kernel()
{
    unsigned int pos = blockIdx.x * blockDim.x + threadIdx.x;
    while (pos < cSim.stride * cSim.outputBuffers)
    {
        cSim.pForce4[pos] = make_float4(0.0f, 0.0f, 0.0f, 0.0f);
        pos += gridDim.x * blockDim.x;
    }
}

void OPENMMCUDA_EXPORT kClearForces(gpuContext gpu)
{
//    printf("kClearForces\n");
    kClearForces_kernel<<<gpu->sim.blocks, 384>>>();
    LAUNCHERROR("kClearForces");
}

__global__ 
__launch_bounds__(384, 1)
void kClearBornSumAndForces_kernel()
{
    unsigned int pos = blockIdx.x * blockDim.x + threadIdx.x;
    while (pos < cSim.stride * cSim.nonbondOutputBuffers)
    {
        cSim.pBornSum[pos] = 0.0f;
        cSim.pBornForce[pos] = 0.0f;
        cSim.pForce4[pos] = make_float4(0.0f, 0.0f, 0.0f, 0.0f);
        pos += gridDim.x * blockDim.x;
    }
    while (pos < cSim.stride * cSim.outputBuffers)
    {
        cSim.pForce4[pos] = make_float4(0.0f, 0.0f, 0.0f, 0.0f);
        pos += gridDim.x * blockDim.x;
    }
}

void kClearBornSumAndForces(gpuContext gpu)
{
  //  printf("kClearBornSumAndForces\n");
    kClearBornSumAndForces_kernel<<<gpu->sim.blocks, 384>>>();
    LAUNCHERROR("kClearBornSumAndForces");
}

__global__ 
__launch_bounds__(384, 1)
void kClearEnergy_kernel()
{
    unsigned int pos = blockIdx.x * blockDim.x + threadIdx.x;
    while (pos < cSim.energyOutputBuffers)
    {
        ((float*)cSim.pEnergy)[pos] = 0.0f;
        pos += gridDim.x * blockDim.x;
    }
}

void kClearEnergy(gpuContext gpu)
{
  //  printf("kClearEnergy\n");
    kClearEnergy_kernel<<<gpu->sim.blocks, 384>>>();
    LAUNCHERROR("kClearEnergy");
}

__global__ 
#if (__CUDA_ARCH__ >= 200)
__launch_bounds__(GF1XX_THREADS_PER_BLOCK, 1)
#elif (__CUDA_ARCH__ >= 130)
__launch_bounds__(GT2XX_THREADS_PER_BLOCK, 1)
#else
__launch_bounds__(G8X_THREADS_PER_BLOCK, 1)
#endif
void kReduceBornSumAndForces_kernel()
{
    unsigned int pos = (blockIdx.x * blockDim.x + threadIdx.x);
   
    // Reduce forces
    while (pos < cSim.stride4)
    {
        float totalForce = 0.0f;
        float* pFt = (float*)cSim.pForce4 + pos;
        int i = cSim.outputBuffers;
        while (i >= 4)
        {
            float f1    = *pFt;
            pFt        += cSim.stride4;
            float f2    = *pFt;
            pFt        += cSim.stride4;
            float f3    = *pFt;
            pFt        += cSim.stride4;
            float f4    = *pFt;
            pFt        += cSim.stride4;
            totalForce += f1 + f2 + f3 + f4;
            i -= 4;
        }
        if (i >= 2)
        {
            float f1    = *pFt;
            pFt        += cSim.stride4;
            float f2    = *pFt;
            pFt        += cSim.stride4;
            totalForce += f1 + f2;
            i -= 2;
        }
        if (i > 0)
        {
            totalForce += *pFt;
        }
        
        pFt = (float*)cSim.pForce4 + pos;
        *pFt = totalForce;
        pos += gridDim.x * blockDim.x;
    }   
    
    
    // Reduce Born Sum
    while (pos - cSim.stride4 < cSim.atoms)
    {
        float sum = 0.0f;
        float* pSt = cSim.pBornSum + pos - cSim.stride4;
        float2 atom = cSim.pObcData[pos - cSim.stride4];
        
    
        // Get summed Born data
        int i = cSim.nonbondOutputBuffers;
        while (i >= 4)
        {
            float f1    = *pSt;
            pSt        += cSim.stride;
            float f2    = *pSt;
            pSt        += cSim.stride;
            float f3    = *pSt;
            pSt        += cSim.stride;
            float f4    = *pSt;
            pSt        += cSim.stride;
            sum += f1 + f2 + f3 + f4;
            i -= 4;
        }
        if (i >= 2)
        {
            float f1    = *pSt;
            pSt        += cSim.stride;
            float f2    = *pSt;
            pSt        += cSim.stride;
            sum += f1 + f2;
            i -= 2;
        }
        if (i > 0)
        {
            sum += *pSt;
        }
       
        // Now calculate Born radius and OBC term.
        cSim.pBornSum[pos - cSim.stride4] = sum; 
        sum                    *= 0.5f * atom.x;
        float sum2              = sum * sum;
        float sum3              = sum * sum2;
        float tanhSum           = tanh(cSim.alphaOBC * sum - cSim.betaOBC * sum2 + cSim.gammaOBC * sum3);
        float nonOffsetRadii    = atom.x + cSim.dielectricOffset;
        float bornRadius        = 1.0f / (1.0f / atom.x - tanhSum / nonOffsetRadii); 
        float obcChain          = atom.x * (cSim.alphaOBC - 2.0f * cSim.betaOBC * sum + 3.0f * cSim.gammaOBC * sum2);
        obcChain                = (1.0f - tanhSum * tanhSum) * obcChain / nonOffsetRadii;              
        cSim.pBornRadii[pos - cSim.stride4] = bornRadius;
        cSim.pObcChain[pos - cSim.stride4]  = obcChain;
        pos += gridDim.x * blockDim.x;
    }
}

void kReduceBornSumAndForces(gpuContext gpu)
{
    //printf("kReduceBornSumAndForces\n");
    kReduceBornSumAndForces_kernel<<<gpu->sim.blocks, gpu->sim.bsf_reduce_threads_per_block>>>();
    LAUNCHERROR("kReduceBornSumAndForces");
}

__global__ 
#if (__CUDA_ARCH__ >= 200)
__launch_bounds__(GF1XX_THREADS_PER_BLOCK, 1)
#elif (__CUDA_ARCH__ >= 130)
__launch_bounds__(GT2XX_THREADS_PER_BLOCK, 1)
#else
__launch_bounds__(G8X_THREADS_PER_BLOCK, 1)
#endif
void kReduceForces_kernel()
{
    unsigned int pos = (blockIdx.x * blockDim.x + threadIdx.x);
   
    // Reduce forces
    while (pos < cSim.stride4)
    {
        float totalForce = 0.0f;
        float* pFt = (float*)cSim.pForce4 + pos;
        int i = cSim.outputBuffers;
        while (i >= 4)
        {
            float f1    = *pFt;
            pFt        += cSim.stride4;
            float f2    = *pFt;
            pFt        += cSim.stride4;
            float f3    = *pFt;
            pFt        += cSim.stride4;
            float f4    = *pFt;
            pFt        += cSim.stride4;
            totalForce += f1 + f2 + f3 + f4;
            i -= 4;
        }
        if (i >= 2)
        {
            float f1    = *pFt;
            pFt        += cSim.stride4;
            float f2    = *pFt;
            pFt        += cSim.stride4;
            totalForce += f1 + f2;
            i -= 2;
        }
        if (i > 0)
        {
            totalForce += *pFt;
        }

        pFt = (float*)cSim.pForce4 + pos;
        *pFt = totalForce;
        pos += gridDim.x * blockDim.x;
    }   
}

void OPENMMCUDA_EXPORT kReduceForces(gpuContext gpu)
{
 //   printf("kReduceForces\n");
    kReduceForces_kernel<<<gpu->sim.blocks, gpu->sim.bsf_reduce_threads_per_block>>>();
    LAUNCHERROR("kReduceForces");
}

double kReduceEnergy(gpuContext gpu)
{
    //printf("kReduceEnergy\n");
    gpu->psEnergy->Download();
    double sum = 0.0;
    for (int i = 0; i < gpu->sim.energyOutputBuffers; i++){
        sum += (*gpu->psEnergy)[i];
    }

    return sum;
}

__global__ 
#if (__CUDA_ARCH__ >= 200)
__launch_bounds__(GF1XX_UPDATE_THREADS_PER_BLOCK, 1)
#elif (__CUDA_ARCH__ >= 130)
__launch_bounds__(GT2XX_UPDATE_THREADS_PER_BLOCK, 1)
#else
__launch_bounds__(G8X_UPDATE_THREADS_PER_BLOCK, 1)
#endif
void kReduceObcGbsaBornForces_kernel()
{
    unsigned int pos = (blockIdx.x * blockDim.x + threadIdx.x);
    float energy = 0.0f;
    while (pos < cSim.atoms)
    {
        float bornRadius = cSim.pBornRadii[pos];
        float obcChain   = cSim.pObcChain[pos];
        float2 obcData   = cSim.pObcData[pos];
        float totalForce = 0.0f;
        float* pFt = cSim.pBornForce + pos;

        int i = cSim.nonbondOutputBuffers;
        while (i >= 4)
        {
            float f1    = *pFt;
            pFt        += cSim.stride;
            float f2    = *pFt;
            pFt        += cSim.stride;
            float f3    = *pFt;
            pFt        += cSim.stride;
            float f4    = *pFt;
            pFt        += cSim.stride;
            totalForce += f1 + f2 + f3 + f4;
            i -= 4;
        }
        if (i >= 2)
        {
            float f1    = *pFt;
            pFt        += cSim.stride;
            float f2    = *pFt;
            pFt        += cSim.stride;
            totalForce += f1 + f2;
            i -= 2;
        }
        if (i > 0)
        {
            totalForce += *pFt;
        }
        float r            = (obcData.x + cSim.dielectricOffset + cSim.probeRadius);
        float ratio6       = pow((obcData.x + cSim.dielectricOffset) / bornRadius, 6.0f);
        float saTerm       = cSim.surfaceAreaFactor * r * r * ratio6;

        totalForce        += saTerm / bornRadius;
        totalForce        *= bornRadius * bornRadius * obcChain;

        energy            += saTerm;

        pFt                = cSim.pBornForce + pos;
        *pFt               = totalForce;

        pos               += gridDim.x * blockDim.x;
    }

    // correct for surface area factor of -6
    cSim.pEnergy[blockIdx.x * blockDim.x + threadIdx.x] += energy / -6.0f;
}

__global__ 
#if (__CUDA_ARCH__ >= 200)
__launch_bounds__(GF1XX_UPDATE_THREADS_PER_BLOCK, 1)
#elif (__CUDA_ARCH__ >= 130)
__launch_bounds__(GT2XX_UPDATE_THREADS_PER_BLOCK, 1)
#else
__launch_bounds__(G8X_UPDATE_THREADS_PER_BLOCK, 1)
#endif
void kReduceGBVIBornForces_kernel()
{
    unsigned int pos = (blockIdx.x * blockDim.x + threadIdx.x);
    float energy = 0.0f;
    while (pos < cSim.atoms)
    {
        float bornRadius  = cSim.pBornRadii[pos];
        float4 gbviData   = cSim.pGBVIData[pos];
        float totalForce  = 0.0f;
        float* pFt        = cSim.pBornForce + pos;

        int i = cSim.nonbondOutputBuffers;
        while (i >= 4)
        {
            float f1    = *pFt;
            pFt        += cSim.stride;
            float f2    = *pFt;
            pFt        += cSim.stride;
            float f3    = *pFt;
            pFt        += cSim.stride;
            float f4    = *pFt;
            pFt        += cSim.stride;
            totalForce += f1 + f2 + f3 + f4;
            i -= 4;
        }
        if (i >= 2)
        {
            float f1    = *pFt;
            pFt        += cSim.stride;
            float f2    = *pFt;
            pFt        += cSim.stride;
            totalForce += f1 + f2;
            i -= 2;
        }
        if (i > 0)
        {
            totalForce += *pFt;
        }

        float ratio         = (gbviData.x/bornRadius);
        float ratio3        = ratio*ratio*ratio;
        energy             -= gbviData.z*ratio3;
        totalForce         += (3.0f*gbviData.z*ratio3)/bornRadius; // 'cavity' term
        float br2           = bornRadius*bornRadius;
        totalForce         *= (1.0f/3.0f)*br2*br2;

        pFt = cSim.pBornForce + pos;
        *pFt = totalForce;
        pos += gridDim.x * blockDim.x;
    }
    cSim.pEnergy[blockIdx.x * blockDim.x + threadIdx.x] += energy;
}

void kReduceObcGbsaBornForces(gpuContext gpu)
{
    //printf("kReduceObcGbsaBornForces\n");
    if( gpu->bIncludeGBSA ){
       kReduceObcGbsaBornForces_kernel<<<gpu->sim.blocks, gpu->sim.bf_reduce_threads_per_block>>>();
       LAUNCHERROR("kReduceObcGbsaBornForces");
    } else if( gpu->bIncludeGBVI ){
       kReduceGBVIBornForces_kernel<<<gpu->sim.blocks, gpu->sim.bf_reduce_threads_per_block>>>();
       LAUNCHERROR("kReduceGBVIBornForces");
    }   

}


