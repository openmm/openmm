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

#define UNROLLXX 0
#define UNROLLXY 0

struct Atom {
    float x;
    float y;
    float z;
    float4 params;
    float fx;
    float fy;
    float fz;
};

static __constant__ cudaGmxSimulation cSim;
static __constant__ Expression<128> forceExp;
static __constant__ Expression<128> energyExp;
static __constant__ Expression<64> combiningRules[4];

void SetCalculateCustomNonbondedForcesSim(gpuContext gpu)
{
    cudaError_t status;
    status = cudaMemcpyToSymbol(cSim, &gpu->sim, sizeof(cudaGmxSimulation));
    RTERROR(status, "cudaMemcpyToSymbol: SetSim copy to cSim failed");
}

void GetCalculateCustomNonbondedForcesSim(gpuContext gpu)
{
    cudaError_t status;
    status = cudaMemcpyFromSymbol(&gpu->sim, cSim, sizeof(cudaGmxSimulation));
    RTERROR(status, "cudaMemcpyFromSymbol: SetSim copy from cSim failed");
}

void SetCustomNonbondedForceExpression(const Expression<128>& expression)
{
    cudaError_t status;
    status = cudaMemcpyToSymbol(forceExp, &expression, sizeof(forceExp));
    RTERROR(status, "SetCustomNonbondedForceExpression: cudaMemcpyToSymbol failed");
}

void SetCustomNonbondedEnergyExpression(const Expression<128>& expression)
{
    cudaError_t status;
    status = cudaMemcpyToSymbol(energyExp, &expression, sizeof(energyExp));
    RTERROR(status, "SetCustomNonbondedEnergyExpression: cudaMemcpyToSymbol failed");
}

void SetCustomNonbondedCombiningRules(const Expression<64>* expressions)
{
    cudaError_t status;
    status = cudaMemcpyToSymbol(combiningRules, expressions, sizeof(combiningRules));
    RTERROR(status, "SetCustomNonbondedCombiningRules: cudaMemcpyToSymbol failed");
}

template<int SIZE>
__device__ float kEvaluateExpression_kernel(Expression<SIZE>* expression, float* stack, float var0, float4 vars1, float4 vars2)
{
    int stackPointer = -1;
    for (int i = 0; i < expression->length; i++)
    {
        switch (expression->op[i])
        {
            case CONSTANT:
                stack[++stackPointer] = expression->arg[i];
                break;
            case VARIABLE0:
                stack[++stackPointer] = var0;
                break;
            case VARIABLE1:
                stack[++stackPointer] = vars1.x;
                break;
            case VARIABLE2:
                stack[++stackPointer] = vars1.y;
                break;
            case VARIABLE3:
                stack[++stackPointer] = vars1.z;
                break;
            case VARIABLE4:
                stack[++stackPointer] = vars1.w;
                break;
            case VARIABLE5:
                stack[++stackPointer] = vars2.x;
                break;
            case VARIABLE6:
                stack[++stackPointer] = vars2.y;
                break;
            case VARIABLE7:
                stack[++stackPointer] = vars2.z;
                break;
            case VARIABLE8:
                stack[++stackPointer] = vars2.w;
                break;
            case ADD:
            {
                float temp = stack[stackPointer];
                stack[stackPointer] = temp+stack[--stackPointer];
                break;
            }
            case SUBTRACT:
            {
                float temp = stack[stackPointer];
                stack[stackPointer] = temp-stack[--stackPointer];
                break;
            }
            case MULTIPLY:
            {
                float temp = stack[stackPointer];
                stack[stackPointer] = temp*stack[--stackPointer];
                break;
            }
            case DIVIDE:
            {
                float temp = stack[stackPointer];
                stack[stackPointer] = temp/stack[--stackPointer];
                break;
            }
            case POWER:
            {
                float temp = stack[stackPointer];
                stack[stackPointer] = pow(temp, stack[--stackPointer]);
                break;
            }
            case NEGATE:
                stack[stackPointer] = -stack[stackPointer];
                break;
            case SQRT:
                stack[stackPointer] = sqrt(stack[stackPointer]);
                break;
            case EXP:
                stack[stackPointer] = exp(stack[stackPointer]);
                break;
            case LOG:
                stack[stackPointer] = log(stack[stackPointer]);
                break;
            case SIN:
                stack[stackPointer] = sin(stack[stackPointer]);
                break;
            case COS:
                stack[stackPointer] = cos(stack[stackPointer]);
                break;
            case SEC:
                stack[stackPointer] = 1.0f/cos(stack[stackPointer]);
                break;
            case CSC:
                stack[stackPointer] = 1.0f/sin(stack[stackPointer]);
                break;
            case TAN:
                stack[stackPointer] = tan(stack[stackPointer]);
                break;
            case COT:
                stack[stackPointer] = 1.0f/tan(stack[stackPointer]);
                break;
            case ASIN:
                stack[stackPointer] = asin(stack[stackPointer]);
                break;
            case ACOS:
                stack[stackPointer] = acos(stack[stackPointer]);
                break;
            case ATAN:
                stack[stackPointer] = atan(stack[stackPointer]);
                break;
            case SQUARE:
            {
                float temp = stack[stackPointer];
                stack[stackPointer] = temp*temp;
                break;
            }
            case CUBE:
            {
                float temp = stack[stackPointer];
                stack[stackPointer] = temp*temp*temp;
                break;
            }
            case RECIPROCAL:
                stack[stackPointer] = 1.0f/stack[stackPointer];
                break;
            case INCREMENT:
                stack[stackPointer] = stack[stackPointer]+1.0f;
                break;
            case DECREMENT:
                stack[stackPointer] = stack[stackPointer]-1.0f;
                break;
        }
    }
    return stack[stackPointer];
}

// Include versions of the kernels for N^2 calculations.

#define METHOD_NAME(a, b) a##N2##b
#include "kCalculateCustomNonbondedForces.h"
#define USE_OUTPUT_BUFFER_PER_WARP
#undef METHOD_NAME
#define METHOD_NAME(a, b) a##N2ByWarp##b
#include "kCalculateCustomNonbondedForces.h"

// Include versions of the kernels with cutoffs.

#undef METHOD_NAME
#undef USE_OUTPUT_BUFFER_PER_WARP
#define USE_CUTOFF
#define METHOD_NAME(a, b) a##Cutoff##b
#include "kCalculateCustomNonbondedForces.h"
#define USE_OUTPUT_BUFFER_PER_WARP
#undef METHOD_NAME
#define METHOD_NAME(a, b) a##CutoffByWarp##b
#include "kCalculateCustomNonbondedForces.h"

// Include versions of the kernels with periodic boundary conditions.

#undef METHOD_NAME
#undef USE_OUTPUT_BUFFER_PER_WARP
#define USE_PERIODIC
#define METHOD_NAME(a, b) a##Periodic##b
#include "kCalculateCustomNonbondedForces.h"
#define USE_OUTPUT_BUFFER_PER_WARP
#undef METHOD_NAME
#define METHOD_NAME(a, b) a##PeriodicByWarp##b
#include "kCalculateCustomNonbondedForces.h"

__global__ void kFindBlockBoundsCutoff_kernel();
__global__ void kFindBlocksWithInteractionsCutoff_kernel();
__global__ void kFindInteractionsWithinBlocksCutoff_kernel(unsigned int* workUnit);
__global__ void kFindBlockBoundsPeriodic_kernel();
__global__ void kFindBlocksWithInteractionsPeriodic_kernel();
__global__ void kFindInteractionsWithinBlocksPeriodic_kernel(unsigned int* workUnit);

void kCalculateCustomNonbondedForces(gpuContext gpu, bool neighborListValid)
{
//    printf("kCalculateCustomNonbondedCutoffForces\n");
    CUDPPResult result;
    switch (gpu->sim.customNonbondedMethod)
    {
        case NO_CUTOFF:
            if (gpu->bOutputBufferPerWarp)
                kCalculateCustomNonbondedN2ByWarpForces_kernel<<<gpu->sim.nonbond_blocks, gpu->sim.nonbond_threads_per_block,
                        (sizeof(Atom)+MAX_STACK_SIZE*sizeof(float))*gpu->sim.nonbond_threads_per_block>>>(gpu->sim.pWorkUnit);
            else
                kCalculateCustomNonbondedN2Forces_kernel<<<gpu->sim.nonbond_blocks, gpu->sim.nonbond_threads_per_block,
                        (sizeof(Atom)+MAX_STACK_SIZE*sizeof(float))*gpu->sim.nonbond_threads_per_block>>>(gpu->sim.pWorkUnit);
            LAUNCHERROR("kCalculateCustomNonbondedN2Forces");
            break;
        case CUTOFF:
            if (!neighborListValid)
            {
                kFindBlockBoundsCutoff_kernel<<<(gpu->psGridBoundingBox->_length+63)/64, 64>>>();
                LAUNCHERROR("kFindBlockBoundsCutoff");
                kFindBlocksWithInteractionsCutoff_kernel<<<gpu->sim.interaction_blocks, gpu->sim.interaction_threads_per_block>>>();
                LAUNCHERROR("kFindBlocksWithInteractionsCutoff");
                result = cudppCompact(gpu->cudpp, gpu->sim.pInteractingWorkUnit, gpu->sim.pInteractionCount,
                        gpu->sim.pWorkUnit, gpu->sim.pInteractionFlag, gpu->sim.workUnits);
                if (result != CUDPP_SUCCESS)
                {
                    printf("Error in cudppCompact: %d\n", result);
                    exit(-1);
                }
                kFindInteractionsWithinBlocksCutoff_kernel<<<gpu->sim.nonbond_blocks, gpu->sim.nonbond_threads_per_block,
                        sizeof(unsigned int)*gpu->sim.nonbond_threads_per_block>>>(gpu->sim.pInteractingWorkUnit);
            }
            if (gpu->bOutputBufferPerWarp)
                kCalculateCustomNonbondedCutoffByWarpForces_kernel<<<gpu->sim.nonbond_blocks, gpu->sim.nonbond_threads_per_block,
                        (sizeof(Atom)+MAX_STACK_SIZE*sizeof(float)+sizeof(float3))*gpu->sim.nonbond_threads_per_block>>>(gpu->sim.pInteractingWorkUnit);
            else
                kCalculateCustomNonbondedCutoffForces_kernel<<<gpu->sim.nonbond_blocks, gpu->sim.nonbond_threads_per_block,
                        (sizeof(Atom)+MAX_STACK_SIZE*sizeof(float)+sizeof(float3))*gpu->sim.nonbond_threads_per_block>>>(gpu->sim.pInteractingWorkUnit);
            LAUNCHERROR("kCalculateCustomNonbondedCutoffForces");
            break;
        case PERIODIC:
            if (!neighborListValid)
            {
                kFindBlockBoundsPeriodic_kernel<<<(gpu->psGridBoundingBox->_length+63)/64, 64>>>();
                LAUNCHERROR("kFindBlockBoundsPeriodic");
                kFindBlocksWithInteractionsPeriodic_kernel<<<gpu->sim.interaction_blocks, gpu->sim.interaction_threads_per_block>>>();
                LAUNCHERROR("kFindBlocksWithInteractionsPeriodic");
                result = cudppCompact(gpu->cudpp, gpu->sim.pInteractingWorkUnit, gpu->sim.pInteractionCount,
                        gpu->sim.pWorkUnit, gpu->sim.pInteractionFlag, gpu->sim.workUnits);
                if (result != CUDPP_SUCCESS)
                {
                    printf("Error in cudppCompact: %d\n", result);
                    exit(-1);
                }
                kFindInteractionsWithinBlocksPeriodic_kernel<<<gpu->sim.nonbond_blocks, gpu->sim.nonbond_threads_per_block,
                        sizeof(unsigned int)*gpu->sim.nonbond_threads_per_block>>>(gpu->sim.pInteractingWorkUnit);
            }
            if (gpu->bOutputBufferPerWarp)
                kCalculateCustomNonbondedPeriodicByWarpForces_kernel<<<gpu->sim.nonbond_blocks, gpu->sim.nonbond_threads_per_block,
                        (sizeof(Atom)+MAX_STACK_SIZE*sizeof(float)+sizeof(float3))*gpu->sim.nonbond_threads_per_block>>>(gpu->sim.pInteractingWorkUnit);
            else
                kCalculateCustomNonbondedPeriodicForces_kernel<<<gpu->sim.nonbond_blocks, gpu->sim.nonbond_threads_per_block,
                        (sizeof(Atom)+MAX_STACK_SIZE*sizeof(float)+sizeof(float3))*gpu->sim.nonbond_threads_per_block>>>(gpu->sim.pInteractingWorkUnit);
            LAUNCHERROR("kCalculateCustomNonbondedPeriodicForces");
            break;
    }
}
