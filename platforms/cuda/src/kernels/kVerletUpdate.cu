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

void SetVerletUpdateSim(gpuContext gpu)
{
    cudaError_t status;
    status = cudaMemcpyToSymbol(cSim, &gpu->sim, sizeof(cudaGmxSimulation));     
    RTERROR(status, "cudaMemcpyToSymbol: SetSim copy to cSim failed");
}

void GetVerletUpdateSim(gpuContext gpu)
{
    cudaError_t status;
    status = cudaMemcpyFromSymbol(&gpu->sim, cSim, sizeof(cudaGmxSimulation));     
    RTERROR(status, "cudaMemcpyFromSymbol: SetSim copy from cSim failed");
}

__global__ void kVerletUpdatePart1_kernel()
{
    unsigned int pos    = threadIdx.x + blockIdx.x * blockDim.x;
    
    while (pos < cSim.atoms)
    {
        float4 apos             = cSim.pPosq[pos];
        float4 velocity         = cSim.pVelm4[pos];
        float4 force            = cSim.pForce4[pos];
        float dtOverMass        = cSim.deltaT*velocity.w;

        cSim.pOldPosq[pos]      = apos;        
        velocity.x             += dtOverMass*force.x;
        velocity.y             += dtOverMass*force.y;
        velocity.z             += dtOverMass*force.z;

#ifndef DeltaShake
        apos.x                 += velocity.x*cSim.deltaT;
        apos.y                 += velocity.y*cSim.deltaT;
        apos.z                 += velocity.z*cSim.deltaT;
#else
        apos.x                  = velocity.x*cSim.deltaT;
        apos.y                  = velocity.y*cSim.deltaT;
        apos.z                  = velocity.z*cSim.deltaT;
#endif
        cSim.pPosqP[pos]        = apos;
        cSim.pVelm4[pos]        = velocity;        
        pos                    += blockDim.x * gridDim.x;
    }
}

__global__ void kVerletUpdatePart1CM_kernel()
{
    extern __shared__ float3 sCM[];
    unsigned int pos    = threadIdx.x + blockIdx.x * blockDim.x;
    float3 CM           = { 0.0f, 0.0f, 0.0f};
    float4 CM1          = { 0.0f, 0.0f, 0.0f, 0.0f };
    
    // Read CM outputs from previous step
    unsigned int cpos = threadIdx.x;
    while (cpos < gridDim.x)
    {
        CM1             = cSim.pLinearMomentum[cpos];
        CM.x           += CM1.x;
        CM.y           += CM1.y;
        CM.z           += CM1.z;
        cpos           += blockDim.x;
    }
    sCM[threadIdx.x].x  = CM.x;
    sCM[threadIdx.x].y  = CM.y;
    sCM[threadIdx.x].z  = CM.z;
    __syncthreads();
    
    // Reduce CM
    unsigned int offset = 1;
    unsigned int mask   = 1;
    while (offset < blockDim.x)
    {
        if (((threadIdx.x & mask) == 0) && (threadIdx.x + offset < blockDim.x))
        {
            sCM[threadIdx.x].x += sCM[threadIdx.x + offset].x;
            sCM[threadIdx.x].y += sCM[threadIdx.x + offset].y;
            sCM[threadIdx.x].z += sCM[threadIdx.x + offset].z;
        }
        mask = 2 * mask + 1;
        offset *= 2;
        __syncthreads();
    }       
    
    while (pos < cSim.atoms)
    {
        float4 apos             = cSim.pPosq[pos];
        float4 velocity         = cSim.pVelm4[pos];
        float4 force            = cSim.pForce4[pos];
        float dtOverMass        = cSim.deltaT*velocity.w;

        cSim.pOldPosq[pos]      = apos;        
        velocity.x             += dtOverMass*force.x-sCM[0].x;
        velocity.y             += dtOverMass*force.y-sCM[0].y;
        velocity.z             += dtOverMass*force.z-sCM[0].z;

#ifndef DeltaShake
        apos.x                 += velocity.x*cSim.deltaT;
        apos.y                 += velocity.y*cSim.deltaT;
        apos.z                 += velocity.z*cSim.deltaT;
#else
        apos.x                  = velocity.x*cSim.deltaT;
        apos.y                  = velocity.y*cSim.deltaT;
        apos.z                  = velocity.z*cSim.deltaT;
#endif

        cSim.pPosqP[pos]        = apos;
        cSim.pVelm4[pos]        = velocity;        
        pos                    += blockDim.x * gridDim.x;
    }
}

void kVerletUpdatePart1(gpuContext gpu)
{
//    printf("kVerletUpdatePart1\n");
    if (gpu->bRemoveCM)
    {
        kVerletUpdatePart1CM_kernel<<<gpu->sim.blocks, gpu->sim.update_threads_per_block, gpu->sim.update_threads_per_block * sizeof(float3)>>>();
        LAUNCHERROR("kVerletUpdatePart1CM");
        gpu->bRemoveCM = false;
    }
    else
    {    
        kVerletUpdatePart1_kernel<<<gpu->sim.blocks, gpu->sim.update_threads_per_block>>>();
        LAUNCHERROR("kVerletUpdatePart1");
    }
}

__global__ void kVerletUpdatePart2_kernel()
{
    unsigned int pos            = threadIdx.x + blockIdx.x * blockDim.x;
    
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
}

__global__ void kVerletUpdatePart2CM_kernel()
{
    extern __shared__ float3 sCM[];
    unsigned int pos            = threadIdx.x + blockIdx.x * blockDim.x;
    float3 CM                   = {0.0f, 0.0f, 0.0f};
    
    while (pos < cSim.atoms)
    {
        float4 velocity         = cSim.pVelm4[pos];
        float4 apos             = cSim.pPosq[pos];
        float4 xPrime           = cSim.pPosqP[pos];
        float mass              = 1.0f / velocity.w;

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

        CM.x                   += mass * velocity.x;
        CM.y                   += mass * velocity.y;
        CM.z                   += mass * velocity.z;
        cSim.pPosq[pos]         = xPrime;
        cSim.pVelm4[pos]        = velocity;
         
        pos                    += blockDim.x * gridDim.x;    
    }
    
    // Scale CM
    CM.x *= cSim.inverseTotalMass;
    CM.y *= cSim.inverseTotalMass;
    CM.z *= cSim.inverseTotalMass;
    sCM[threadIdx.x] = CM;
    __syncthreads();
    
    // Reduce CM for CTA
    unsigned int offset = 1;
    unsigned int mask   = 1;
    while (offset < blockDim.x)
    {
        if (((threadIdx.x & mask) == 0) && (threadIdx.x + offset < blockDim.x))
        {
            sCM[threadIdx.x].x += sCM[threadIdx.x + offset].x;
            sCM[threadIdx.x].y += sCM[threadIdx.x + offset].y;
            sCM[threadIdx.x].z += sCM[threadIdx.x + offset].z;
        }
        mask = 2 * mask + 1;
        offset *= 2;
        __syncthreads();
    }
    if (threadIdx.x == 0)
    {
        float4 CM;
        CM.x                                = sCM[0].x;
        CM.y                                = sCM[0].y;
        CM.z                                = sCM[0].z;
        CM.w                                = 0.0f;
        cSim.pLinearMomentum[blockIdx.x]    = CM;
    }  
}

void kVerletUpdatePart2(gpuContext gpu)
{
//    printf("kVerletUpdatePart2\n");
    if (gpu->bCalculateCM)
    {
        kVerletUpdatePart2CM_kernel<<<gpu->sim.blocks, gpu->sim.update_threads_per_block, gpu->sim.update_threads_per_block * sizeof(float3)>>>();
        LAUNCHERROR("kVerletUpdatePart2CM");
        gpu->bCalculateCM = false;
        gpu->bRemoveCM = true;
    }
    else
    {
        kVerletUpdatePart2_kernel<<<gpu->sim.blocks, gpu->sim.update_threads_per_block>>>();
        LAUNCHERROR("kVerletUpdatePart2");
    }
}

