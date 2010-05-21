/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2010 Stanford University and the Authors.           *
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
#include "cudatypes.h"

static __constant__ cudaGmxSimulation cSim;
static __constant__ Expression<256> forceExp;
static __constant__ Expression<256> energyExp;

#include "kEvaluateExpression.h"

void SetCalculateCustomAngleForcesSim(gpuContext gpu)
{
    cudaError_t status;
    status = cudaMemcpyToSymbol(cSim, &gpu->sim, sizeof(cudaGmxSimulation));
    RTERROR(status, "cudaMemcpyToSymbol: SetSim copy to cSim failed");
}

void GetCalculateCustomAngleForcesSim(gpuContext gpu)
{
    cudaError_t status;
    status = cudaMemcpyFromSymbol(&gpu->sim, cSim, sizeof(cudaGmxSimulation));
    RTERROR(status, "cudaMemcpyFromSymbol: SetSim copy from cSim failed");
}

void SetCustomAngleForceExpression(const Expression<256>& expression)
{
    cudaError_t status;
    status = cudaMemcpyToSymbol(forceExp, &expression, sizeof(forceExp));
    RTERROR(status, "SetCustomAngleForceExpression: cudaMemcpyToSymbol failed");
}

void SetCustomAngleEnergyExpression(const Expression<256>& expression)
{
    cudaError_t status;
    status = cudaMemcpyToSymbol(energyExp, &expression, sizeof(energyExp));
    RTERROR(status, "SetCustomAngleEnergyExpression: cudaMemcpyToSymbol failed");
}

void SetCustomAngleGlobalParams(const vector<float>& paramValues)
{
    cudaError_t status;
    status = cudaMemcpyToSymbol(globalParams, &paramValues[0], paramValues.size()*sizeof(float));
    RTERROR(status, "SetCustomAngleGlobalParams: cudaMemcpyToSymbol failed");
}

#define DOT3(v1, v2) (v1.x*v2.x + v1.y*v2.y + v1.z*v2.z)

#define CROSS_PRODUCT(v1, v2) make_float3(v1.y*v2.z - v1.z*v2.y, v1.z*v2.x - v1.x*v2.z, v1.x*v2.y - v1.y*v2.x)

__global__ 
#if (__CUDA_ARCH__ >= 200)
__launch_bounds__(1024, 1)
#elif (__CUDA_ARCH__ >= 130)
__launch_bounds__(512, 1)
#else
__launch_bounds__(256, 1)
#endif
void kCalculateCustomAngleForces_kernel()
{
    extern __shared__ float stack[];
    float* variables = (float*) &stack[cSim.customExpressionStackSize*blockDim.x];
    unsigned int pos = blockIdx.x * blockDim.x + threadIdx.x;
    float totalEnergy = 0.0f;

    while (pos < cSim.customAngles)
    {
        int4 atom       = cSim.pCustomAngleID1[pos];
        int2 atom2      = cSim.pCustomAngleID2[pos];
        float4 params   = cSim.pCustomAngleParams[pos];
        float4 a1       = cSim.pPosq[atom.x];
        float4 a2       = cSim.pPosq[atom.y];
        float4 a3       = cSim.pPosq[atom.z];
        float3 v0 = make_float3(a2.x-a1.x, a2.y-a1.y, a2.z-a1.z);
        float3 v1 = make_float3(a2.x-a3.x, a2.y-a3.y, a2.z-a3.z);
        float3 cp = CROSS_PRODUCT(v0, v1);
        float rp = DOT3(cp, cp);
        rp = max(sqrt(rp), 1.0e-06f);
        float r21 = DOT3(v0, v0);
        float r23 = DOT3(v1, v1);
        float dot = DOT3(v0, v1);
        float cosine = dot/sqrt(r21*r23);
        VARIABLE(0) = acos(cosine);
        VARIABLE(1) = params.x;
        VARIABLE(2) = params.y;
        VARIABLE(3) = params.z;
        VARIABLE(4) = params.w;
        float dEdR = kEvaluateExpression_kernel(&forceExp, stack, variables);
        totalEnergy += kEvaluateExpression_kernel(&energyExp, stack, variables);
        float termA =  dEdR/(r21*rp);
        float termC = -dEdR/(r23*rp);
        float3 c21 = CROSS_PRODUCT(v0, cp);
        float3 c23 = CROSS_PRODUCT(v1, cp);
        c21.x *= termA;
        c21.y *= termA;
        c21.z *= termA;
        c23.x *= termC;
        c23.y *= termC;
        c23.z *= termC;
        unsigned int offsetA  = atom.x + atom.w * cSim.stride;
        unsigned int offsetB  = atom.y + atom2.x * cSim.stride;
        unsigned int offsetC  = atom.z + atom2.y * cSim.stride;
        float4 forceA         = cSim.pForce4[offsetA];
        float4 forceB         = cSim.pForce4[offsetB];
        float4 forceC         = cSim.pForce4[offsetC];
        forceA.x             += c21.x;
        forceA.y             += c21.y;
        forceA.z             += c21.z;
        forceB.x             -= c21.x+c23.x;
        forceB.y             -= c21.y+c23.y;
        forceB.z             -= c21.z+c23.z;
        forceC.x             += c23.x;
        forceC.y             += c23.y;
        forceC.z             += c23.z;
        cSim.pForce4[offsetA] = forceA;
        cSim.pForce4[offsetB] = forceB;
        cSim.pForce4[offsetC] = forceC;
        pos                  += blockDim.x * gridDim.x;
    }
    cSim.pEnergy[blockIdx.x * blockDim.x + threadIdx.x] += totalEnergy;
}

void kCalculateCustomAngleForces(gpuContext gpu)
{
//    printf("kCalculateCustomAngleForces\n");
    int memoryPerThread = (gpu->sim.customExpressionStackSize+9)*sizeof(float);
    int maxThreads = (gpu->sharedMemoryPerBlock-16)/memoryPerThread;
    int threads = min(gpu->sim.localForces_threads_per_block, (maxThreads/64)*64);
    kCalculateCustomAngleForces_kernel<<<gpu->sim.blocks, threads, memoryPerThread*threads>>>();
    LAUNCHERROR("kCalculateCustomAngleForces");
}
