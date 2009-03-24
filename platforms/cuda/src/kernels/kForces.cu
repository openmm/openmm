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

#define FABS(a) ((a) > 0.0f ? (a) : -(a))

static __constant__ cudaGmxSimulation cSim;

void SetForcesSim(gpuContext gpu)
{
    cudaError_t status;
    status = cudaMemcpyToSymbol(cSim, &gpu->sim, sizeof(cudaGmxSimulation));     
    RTERROR(status, "cudaMemcpyToSymbol: SetSim copy to cSim failed");
}

void GetForcesSim(gpuContext gpu)
{
    cudaError_t status;
    status = cudaMemcpyFromSymbol(&gpu->sim, cSim, sizeof(cudaGmxSimulation));     
    RTERROR(status, "cudaMemcpyFromSymbol: SetSim copy from cSim failed");
}

__global__ void kClearForces_kernel()
{
    unsigned int pos = blockIdx.x * blockDim.x + threadIdx.x;
    while (pos < cSim.stride * cSim.outputBuffers)
    {
        cSim.pForce4[pos] = make_float4(0.0f, 0.0f, 0.0f, 0.0f);
        pos += gridDim.x * blockDim.x;
    }
}

void kClearForces(gpuContext gpu)
{
//    printf("kClearForces\n");
    kClearForces_kernel<<<gpu->sim.blocks, 384>>>();
    LAUNCHERROR("kClearForces");
}

__global__ void kClearBornForces_kernel()
{
    unsigned int pos = blockIdx.x * blockDim.x + threadIdx.x;
    while (pos < cSim.stride * cSim.nonbondOutputBuffers)
    {
        ((float*)cSim.pBornForce)[pos] = 0.0f;
        pos += gridDim.x * blockDim.x;
    }
}

void kClearBornForces(gpuContext gpu)
{
  //  printf("kClearBornForces\n");
    kClearBornForces_kernel<<<gpu->sim.blocks, 384>>>();
    LAUNCHERROR("kClearBornForces");
}

__global__ void kReduceBornSumAndForces_kernel()
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
    
#if 0
    //gpuDumpObcLoop1( gpu );
	 /*
    gpu->psForce4->Download();
    for (int i = 0; i < gpu->natoms; i++)
    {
        printf("%4d: %12.6f %12.6f %12.6f\n", i, 
            gpu->psForce4->_pSysStream[0][i].x,
            gpu->psForce4->_pSysStream[0][i].y,
            gpu->psForce4->_pSysStream[0][i].z
        );
    } */
#endif
}

__global__ void kReduceForces_kernel()
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

void kReduceForces(gpuContext gpu)
{
 //   printf("kReduceForces\n");
    kReduceForces_kernel<<<gpu->sim.blocks, gpu->sim.bsf_reduce_threads_per_block>>>();
    LAUNCHERROR("kReduceForces");
}


__global__ void kReduceObcGbsaBornForces_kernel()
{
    unsigned int pos = (blockIdx.x * blockDim.x + threadIdx.x);
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
        //float saTerm       = cSim.surfaceAreaFactor * r * r * ratio6;
        float saTerm       = cSim.surfaceAreaFactor * r * r * ratio6;
        totalForce        += saTerm / bornRadius; // 1.102 == Temp mysterious fudge factor, FIX FIX FIX

        totalForce *= bornRadius * bornRadius * obcChain;

        pFt = cSim.pBornForce + pos;
        *pFt = totalForce;
        pos += gridDim.x * blockDim.x;
    }
}

void kReduceObcGbsaBornForces(gpuContext gpu)
{
    //printf("kReduceObcGbsaBornForces\n");
    kReduceObcGbsaBornForces_kernel<<<gpu->sim.blocks, gpu->sim.bf_reduce_threads_per_block>>>();
    LAUNCHERROR("kReduceObcGbsaBornForces");
}


