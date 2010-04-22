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

static __constant__ cudaGmxSimulation cSim;

void SetRandomSim(gpuContext gpu)
{
    cudaError_t status;
    status = cudaMemcpyToSymbol(cSim, &gpu->sim, sizeof(cudaGmxSimulation));     
    RTERROR(status, "cudaMemcpyToSymbol: SetSim copy to cSim failed");
}

void GetRandomSim(gpuContext gpu)
{
    cudaError_t status;
    status = cudaMemcpyFromSymbol(&gpu->sim, cSim, sizeof(cudaGmxSimulation));     
    RTERROR(status, "cudaMemcpyFromSymbol: SetSim copy from cSim failed");
}

extern __shared__ float3 sRand[];


__global__ 
#if (__CUDA_ARCH__ >= 200)
__launch_bounds__(GF1XX_RANDOM_THREADS_PER_BLOCK, 1)
#elif (__CUDA_ARCH__ >= 130)
__launch_bounds__(GT2XX_RANDOM_THREADS_PER_BLOCK, 1)
#else
__launch_bounds__(G8X_RANDOM_THREADS_PER_BLOCK, 1)
#endif
void kGenerateRandoms_kernel()
{
    unsigned int pos            = blockIdx.x * blockDim.x + threadIdx.x;
    unsigned int increment      = blockDim.x * gridDim.x;
    
    // Read generator state
    uint4 state                 = cSim.pRandomSeed[pos];
    unsigned int carry          = 0;
    
    float4 random4;
    float2 random2;
    while (pos < cSim.totalRandoms)
    {
        
        // Generate 6 randoms in GRF
        unsigned int pos1       = threadIdx.x;
        for (int i = 0; i < 2; i++)
        {
            state.x             = state.x * 69069 + 1;
            state.y            ^= state.y << 13;
            state.y            ^= state.y >> 17;
            state.y            ^= state.y << 5;
            unsigned int k      = (state.z >> 2) + (state.w >> 3) + (carry >> 2);
            unsigned int m      = state.w + state.w + state.z + carry;
            state.z             = state.w;
            state.w             = m;
            carry               = k >> 30;
            float x1            = (float)max(state.x + state.y + state.w, 0x00000001) / (float)0xffffffff;
            state.x             = state.x * 69069 + 1;
            state.y            ^= state.y << 13;
            state.y            ^= state.y >> 17;
            state.y            ^= state.y << 5;
            x1                  = sqrt(-2.0f * log(x1));
            k                   = (state.z >> 2) + (state.w >> 3) + (carry >> 2);
            m                   = state.w + state.w + state.z + carry;
            state.z             = state.w;
            state.w             = m;
            carry               = k >> 30;
            float x2            = (float)(state.x + state.y + state.w) / (float)0xffffffff;
            
            state.x             = state.x * 69069 + 1;
            state.y            ^= state.y << 13;
            state.y            ^= state.y >> 17;
            state.y            ^= state.y << 5;
            sRand[pos1].x       = x1 * cos(2.0f * 3.14159265f * x2);
            k                   = (state.z >> 2) + (state.w >> 3) + (carry >> 2);
            m                   = state.w + state.w + state.z + carry;
            state.z             = state.w;
            state.w             = m;
            carry               = k >> 30;
            float x3            = (float)max(state.x + state.y + state.w, 0x00000001) / (float)0xffffffff;
            state.x             = state.x * 69069 + 1;
            state.y            ^= state.y << 13;
            state.y            ^= state.y >> 17;
            state.y            ^= state.y << 5;
            x3                  = sqrt(-2.0f * log(x3));
            k                   = (state.z >> 2) + (state.w >> 3) + (carry >> 2);
            m                   = state.w + state.w + state.z + carry;
            state.z             = state.w;
            state.w             = m;
            carry               = k >> 30;
            float x4            = (float)(state.x + state.y + state.w) / (float)0xffffffff;
            
            state.x             = state.x * 69069 + 1;
            state.y            ^= state.y << 13;
            state.y            ^= state.y >> 17;
            state.y            ^= state.y << 5;
            sRand[pos1].y       = x3 * cos(2.0f * 3.14159265f * x4);
            k                   = (state.z >> 2) + (state.w >> 3) + (carry >> 2);
            m                   = state.w + state.w + state.z + carry;
            state.z             = state.w;
            state.w             = m;
            carry               = k >> 30;
            float x5            = (float)max(state.x + state.y + state.w, 0x00000001) / (float)0xffffffff;
            state.x             = state.x * 69069 + 1;
            state.y            ^= state.y << 13;
            state.y            ^= state.y >> 17;
            state.y            ^= state.y << 5;
            x5                  = sqrt(-2.0f * log(x5));
            k                   = (state.z >> 2) + (state.w >> 3) + (carry >> 2);
            m                   = state.w + state.w + state.z + carry;
            state.z             = state.w;
            state.w             = m;
            carry               = k >> 30;
            float x6            = (float)(state.x + state.y + state.w) / (float)0xffffffff;
            sRand[pos1].z       = x5 * cos(2.0f * 3.14159265f * x6); 
            pos1               += blockDim.x;
        }
        
        // Output final randoms
        random4.x               = sRand[threadIdx.x].x;
        random4.y               = sRand[threadIdx.x].y;
        random4.z               = sRand[threadIdx.x].z;
        random4.w               = sRand[threadIdx.x + blockDim.x].x;
        cSim.pRandom4[pos]     = random4;
        random2.x               = sRand[threadIdx.x + blockDim.x].y;
        random2.y               = sRand[threadIdx.x + blockDim.x].z;
        cSim.pRandom2[pos]     = random2;
        
        pos += increment;
    }
    
    
    // Write generator state
    pos                     = blockIdx.x * blockDim.x + threadIdx.x;
    cSim.pRandomSeed[pos]   = state;
}

void kGenerateRandoms(gpuContext gpu)
{
    kGenerateRandoms_kernel<<<gpu->sim.blocks, gpu->sim.random_threads_per_block, gpu->sim.random_threads_per_block * 2 * sizeof(float3)>>>();
}

