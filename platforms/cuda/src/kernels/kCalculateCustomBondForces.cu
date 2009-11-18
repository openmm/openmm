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
#include "cudatypes.h"

static __constant__ cudaGmxSimulation cSim;
static __constant__ Expression<128> forceExp;
static __constant__ Expression<128> energyExp;

#include "kEvaluateExpression.h"

void SetCalculateCustomBondForcesSim(gpuContext gpu)
{
    cudaError_t status;
    status = cudaMemcpyToSymbol(cSim, &gpu->sim, sizeof(cudaGmxSimulation));
    RTERROR(status, "cudaMemcpyToSymbol: SetSim copy to cSim failed");
}

void GetCalculateCustomBondForcesSim(gpuContext gpu)
{
    cudaError_t status;
    status = cudaMemcpyFromSymbol(&gpu->sim, cSim, sizeof(cudaGmxSimulation));
    RTERROR(status, "cudaMemcpyFromSymbol: SetSim copy from cSim failed");
}

void SetCustomBondForceExpression(const Expression<128>& expression)
{
    cudaError_t status;
    status = cudaMemcpyToSymbol(forceExp, &expression, sizeof(forceExp));
    RTERROR(status, "SetCustomBondForceExpression: cudaMemcpyToSymbol failed");
}

void SetCustomBondEnergyExpression(const Expression<128>& expression)
{
    cudaError_t status;
    status = cudaMemcpyToSymbol(energyExp, &expression, sizeof(energyExp));
    RTERROR(status, "SetCustomBondEnergyExpression: cudaMemcpyToSymbol failed");
}

void SetCustomBondGlobalParams(const vector<float>& paramValues)
{
    cudaError_t status;
    status = cudaMemcpyToSymbol(globalParams, &paramValues[0], paramValues.size()*sizeof(float));
    RTERROR(status, "SetCustomBondGlobalParams: cudaMemcpyToSymbol failed");
}


__global__ void kCalculateCustomBondForces_kernel()
{
    extern __shared__ float stack[];
    float* variables = (float*) &stack[cSim.customExpressionStackSize*blockDim.x];
    unsigned int pos = blockIdx.x * blockDim.x + threadIdx.x;
    float totalEnergy = 0.0f;

    while (pos < cSim.customBonds)
    {
        int4 atom       = cSim.pCustomBondID[pos];
        float4 params   = cSim.pCustomBondParams[pos];
        float4 a1       = cSim.pPosq[atom.x];
        float4 a2       = cSim.pPosq[atom.y];
        float dx        = a1.x - a2.x;
        float dy        = a1.y - a2.y;
        float dz        = a1.z - a2.z;
        float r         = sqrt(dx*dx + dy*dy + dz*dz);
        float invR      = 1.0f/r;
        VARIABLE(0)     = r;
        VARIABLE(1)     = params.x;
        VARIABLE(2)     = params.y;
        VARIABLE(3)     = params.z;
        VARIABLE(4)     = params.w;
        float dEdR      = -kEvaluateExpression_kernel(&forceExp, stack, variables)*invR;
        float energy    = kEvaluateExpression_kernel(&energyExp, stack, variables);
        totalEnergy          += energy;
        dx                   *= dEdR;
        dy                   *= dEdR;
        dz                   *= dEdR;
        unsigned int offsetA  = atom.x + atom.z * cSim.stride;
        unsigned int offsetB  = atom.y + atom.w * cSim.stride;
        float4 forceA         = cSim.pForce4[offsetA];
        float4 forceB         = cSim.pForce4[offsetB];
        forceA.x             += dx;
        forceA.y             += dy;
        forceA.z             += dz;
        forceB.x             -= dx;
        forceB.y             -= dy;
        forceB.z             -= dz;
        cSim.pForce4[offsetA] = forceA;
        cSim.pForce4[offsetB] = forceB;
        pos                  += blockDim.x * gridDim.x;
    }
    cSim.pEnergy[blockIdx.x * blockDim.x + threadIdx.x] += totalEnergy;
}

void kCalculateCustomBondForces(gpuContext gpu)
{
//    printf("kCalculateCustomBondForces\n");
    kCalculateCustomBondForces_kernel<<<gpu->sim.blocks, gpu->sim.localForces_threads_per_block,
            gpu->sim.customExpressionStackSize*sizeof(float)*gpu->sim.localForces_threads_per_block>>>();
    LAUNCHERROR("kCalculateCustomBondForces");
}
