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
static __constant__ Expression<128> forceExpX;
static __constant__ Expression<128> forceExpY;
static __constant__ Expression<128> forceExpZ;
static __constant__ Expression<128> energyExp;

#include "kEvaluateExpression.h"

void SetCalculateCustomExternalForcesSim(gpuContext gpu)
{
    cudaError_t status;
    status = cudaMemcpyToSymbol(cSim, &gpu->sim, sizeof(cudaGmxSimulation));
    RTERROR(status, "cudaMemcpyToSymbol: SetSim copy to cSim failed");
}

void GetCalculateCustomExternalForcesSim(gpuContext gpu)
{
    cudaError_t status;
    status = cudaMemcpyFromSymbol(&gpu->sim, cSim, sizeof(cudaGmxSimulation));
    RTERROR(status, "cudaMemcpyFromSymbol: SetSim copy from cSim failed");
}

void SetCustomExternalForceExpressions(const Expression<128>& expressionX, const Expression<128>& expressionY, const Expression<128>& expressionZ)
{
    cudaError_t status;
    status = cudaMemcpyToSymbol(forceExpX, &expressionX, sizeof(forceExpX));
    status = cudaMemcpyToSymbol(forceExpY, &expressionY, sizeof(forceExpY));
    status = cudaMemcpyToSymbol(forceExpZ, &expressionZ, sizeof(forceExpZ));
    RTERROR(status, "SetCustomExternalForceExpression: cudaMemcpyToSymbol failed");
}

void SetCustomExternalEnergyExpression(const Expression<128>& expression)
{
    cudaError_t status;
    status = cudaMemcpyToSymbol(energyExp, &expression, sizeof(energyExp));
    RTERROR(status, "SetCustomExternalEnergyExpression: cudaMemcpyToSymbol failed");
}

void SetCustomExternalGlobalParams(float* paramValues)
{
    cudaError_t status;
    status = cudaMemcpyToSymbol(globalParams, paramValues, sizeof(globalParams));
    RTERROR(status, "SetCustomExternalGlobalParams: cudaMemcpyToSymbol failed");
}


__global__ void kCalculateCustomExternalForces_kernel()
{
    extern __shared__ float stack[];
    float* variables = (float*) &stack[cSim.customExpressionStackSize*blockDim.x];
    unsigned int index = blockIdx.x * blockDim.x + threadIdx.x;
    float totalEnergy = 0.0f;

    while (index < cSim.customExternals)
    {
        int atom        = cSim.pCustomExternalID[index];
        float4 params   = cSim.pCustomExternalParams[index];
        float4 pos      = cSim.pPosq[atom];
        VARIABLE(0)     = pos.x;
        VARIABLE(1)     = pos.y;
        VARIABLE(2)     = pos.z;
        VARIABLE(3)     = params.x;
        VARIABLE(4)     = params.y;
        VARIABLE(5)     = params.z;
        VARIABLE(6)     = params.w;
        totalEnergy       += kEvaluateExpression_kernel(&energyExp, stack, variables);;
        float4 force       = cSim.pForce4[atom];
        force.x           -= kEvaluateExpression_kernel(&forceExpX, stack, variables);
        force.y           -= kEvaluateExpression_kernel(&forceExpY, stack, variables);
        force.z           -= kEvaluateExpression_kernel(&forceExpZ, stack, variables);
        cSim.pForce4[atom] = force;
        index             += blockDim.x * gridDim.x;
    }
    cSim.pEnergy[blockIdx.x * blockDim.x + threadIdx.x] += totalEnergy;
}

void kCalculateCustomExternalForces(gpuContext gpu)
{
//    printf("kCalculateCustomExternalForces\n");
    kCalculateCustomExternalForces_kernel<<<gpu->sim.blocks, gpu->sim.localForces_threads_per_block,
            gpu->sim.customExpressionStackSize*sizeof(float)*gpu->sim.localForces_threads_per_block>>>();
    LAUNCHERROR("kCalculateCustomExternalForces");
}
