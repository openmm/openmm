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

void SetCalculateCustomTorsionForcesSim(gpuContext gpu)
{
    cudaError_t status;
    status = cudaMemcpyToSymbol(cSim, &gpu->sim, sizeof(cudaGmxSimulation));
    RTERROR(status, "cudaMemcpyToSymbol: SetSim copy to cSim failed");
}

void GetCalculateCustomTorsionForcesSim(gpuContext gpu)
{
    cudaError_t status;
    status = cudaMemcpyFromSymbol(&gpu->sim, cSim, sizeof(cudaGmxSimulation));
    RTERROR(status, "cudaMemcpyFromSymbol: SetSim copy from cSim failed");
}

void SetCustomTorsionForceExpression(const Expression<256>& expression)
{
    cudaError_t status;
    status = cudaMemcpyToSymbol(forceExp, &expression, sizeof(forceExp));
    RTERROR(status, "SetCustomTorsionForceExpression: cudaMemcpyToSymbol failed");
}

void SetCustomTorsionEnergyExpression(const Expression<256>& expression)
{
    cudaError_t status;
    status = cudaMemcpyToSymbol(energyExp, &expression, sizeof(energyExp));
    RTERROR(status, "SetCustomTorsionEnergyExpression: cudaMemcpyToSymbol failed");
}

void SetCustomTorsionGlobalParams(const vector<float>& paramValues)
{
    cudaError_t status;
    status = cudaMemcpyToSymbol(globalParams, &paramValues[0], paramValues.size()*sizeof(float));
    RTERROR(status, "SetCustomTorsionGlobalParams: cudaMemcpyToSymbol failed");
}

#define LOCAL_HACK_PI 3.1415926535897932384626433832795

#define DOT3(v1, v2) (v1.x*v2.x + v1.y*v2.y + v1.z*v2.z)

#define CROSS_PRODUCT(v1, v2) make_float3(v1.y*v2.z - v1.z*v2.y, v1.z*v2.x - v1.x*v2.z, v1.x*v2.y - v1.y*v2.x)

#define GETNORMEDDOTPRODUCT(v1, v2, dp) \
{ \
    dp          = DOT3(v1, v2); \
    float norm1 = DOT3(v1, v1); \
    float norm2 = DOT3(v2, v2); \
    dp /= sqrt(norm1 * norm2); \
    dp = min(dp, 1.0f); \
    dp = max(dp, -1.0f); \
}

#define GETANGLEBETWEENTWOVECTORS(v1, v2, angle) \
{ \
    float dp; \
    GETNORMEDDOTPRODUCT(v1, v2, dp); \
    if (dp > 0.99f || dp < -0.99f) { \
        float3 cross = CROSS_PRODUCT(v1, v2); \
        float scale = DOT3(v1, v1)*DOT3(v2, v2); \
        angle = asin(sqrt(DOT3(cross, cross)/scale)); \
        if (dp < 0.0f) \
            angle = LOCAL_HACK_PI-angle; \
    } \
    else { \
        angle = acos(dp); \
    } \
}

#define GETDIHEDRALANGLEBETWEENTHREEVECTORS(vector1, vector2, vector3, signVector, cp0, cp1, angle) \
{ \
    cp0 = CROSS_PRODUCT(vector1, vector2); \
    cp1 = CROSS_PRODUCT(vector2, vector3); \
    GETANGLEBETWEENTWOVECTORS(cp0, cp1, angle); \
    float dp = DOT3(signVector, cp1); \
    angle = (dp >= 0) ? angle : -angle; \
}

__global__ 
#if (__CUDA_ARCH__ >= 200)
__launch_bounds__(GF1XX_LOCALFORCES_THREADS_PER_BLOCK, 1)
#elif (__CUDA_ARCH__ >= 130)
__launch_bounds__(GT2XX_LOCALFORCES_THREADS_PER_BLOCK, 1)
#else
__launch_bounds__(G8X_LOCALFORCES_THREADS_PER_BLOCK, 1)
#endif
void kCalculateCustomTorsionForces_kernel()
{
    extern __shared__ float stack[];
    float* variables = (float*) &stack[cSim.customExpressionStackSize*blockDim.x];
    unsigned int pos = blockIdx.x * blockDim.x + threadIdx.x;
    float totalEnergy = 0.0f;

    while (pos < cSim.customTorsions)
    {
        int4 atom       = cSim.pCustomTorsionID1[pos];
        int4 atom2      = cSim.pCustomTorsionID2[pos];
        float4 params   = cSim.pCustomTorsionParams[pos];
        float4 a1       = cSim.pPosq[atom.x];
        float4 a2       = cSim.pPosq[atom.y];
        float4 a3       = cSim.pPosq[atom.z];
        float4 a4       = cSim.pPosq[atom.w];
        float3 v0 = make_float3(a1.x-a2.x, a1.y-a2.y, a1.z-a2.z);
        float3 v1 = make_float3(a3.x-a2.x, a3.y-a2.y, a3.z-a2.z);
        float3 v2 = make_float3(a3.x-a4.x, a3.y-a4.y, a3.z-a4.z);
        float3 cp0, cp1;
        float dihedralAngle;
        GETDIHEDRALANGLEBETWEENTHREEVECTORS(v0, v1, v2, v0, cp0, cp1, dihedralAngle);
        VARIABLE(0) = dihedralAngle;
        VARIABLE(1) = params.x;
        VARIABLE(2) = params.y;
        VARIABLE(3) = params.z;
        VARIABLE(4) = params.w;
        float dEdAngle = kEvaluateExpression_kernel(&forceExp, stack, variables);
        totalEnergy += kEvaluateExpression_kernel(&energyExp, stack, variables);
        float normBC = sqrt(DOT3(v1, v1));
        float dp = 1.0f / DOT3(v1, v1);
        float4 ff = make_float4((-dEdAngle*normBC)/DOT3(cp0, cp0), DOT3(v0, v1)*dp, DOT3(v2, v1)*dp, (dEdAngle*normBC)/DOT3(cp1, cp1));
        float3 internalF0 = make_float3(ff.x*cp0.x, ff.x*cp0.y, ff.x*cp0.z);
        float3 internalF3 = make_float3(ff.w*cp1.x, ff.w*cp1.y, ff.w*cp1.z);
        float3 s = make_float3(ff.y*internalF0.x - ff.z*internalF3.x,
                               ff.y*internalF0.y - ff.z*internalF3.y,
                               ff.y*internalF0.z - ff.z*internalF3.z);
        unsigned int offsetA  = atom.x + atom2.x * cSim.stride;
        unsigned int offsetB  = atom.y + atom2.y * cSim.stride;
        unsigned int offsetC  = atom.z + atom2.z * cSim.stride;
        unsigned int offsetD  = atom.w + atom2.w * cSim.stride;
        float4 forceA         = cSim.pForce4[offsetA];
        float4 forceB         = cSim.pForce4[offsetB];
        float4 forceC         = cSim.pForce4[offsetC];
        float4 forceD         = cSim.pForce4[offsetD];
        forceA.x             += internalF0.x;
        forceA.y             += internalF0.y;
        forceA.z             += internalF0.z;
        forceB.x             += -internalF0.x + s.x;
        forceB.y             += -internalF0.y + s.y;
        forceB.z             += -internalF0.z + s.z;
        forceC.x             += -internalF3.x - s.x;
        forceC.y             += -internalF3.y - s.y;
        forceC.z             += -internalF3.z - s.z;
        forceD.x             += internalF3.x;
        forceD.y             += internalF3.y;
        forceD.z             += internalF3.z;
        cSim.pForce4[offsetA] = forceA;
        cSim.pForce4[offsetB] = forceB;
        cSim.pForce4[offsetC] = forceC;
        cSim.pForce4[offsetD] = forceD;
        pos                  += blockDim.x * gridDim.x;
    }
    cSim.pEnergy[blockIdx.x * blockDim.x + threadIdx.x] += totalEnergy;
}

void kCalculateCustomTorsionForces(gpuContext gpu)
{
//    printf("kCalculateCustomTorsionForces\n");
    int memoryPerThread = (gpu->sim.customExpressionStackSize+9)*sizeof(float);
    int maxThreads = (gpu->sharedMemoryPerBlock-16)/memoryPerThread;
    int threads = min(gpu->sim.localForces_threads_per_block, (maxThreads/64)*64);
    kCalculateCustomTorsionForces_kernel<<<gpu->sim.blocks, threads, memoryPerThread*threads>>>();
    LAUNCHERROR("kCalculateCustomTorsionForces");
}
