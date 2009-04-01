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


static __constant__ cudaGmxSimulation cSim;

void SetLincsSim(gpuContext gpu)
{
    cudaError_t status;
    status = cudaMemcpyToSymbol(cSim, &gpu->sim, sizeof(cudaGmxSimulation));
    RTERROR(status, "cudaMemcpyToSymbol: SetSim copy to cSim failed");
}

void GetLincsSim(gpuContext gpu)
{
    cudaError_t status;
    status = cudaMemcpyFromSymbol(&gpu->sim, cSim, sizeof(cudaGmxSimulation));
    RTERROR(status, "cudaMemcpyFromSymbol: SetSim copy from cSim failed");
}

/**
 * Synchronize all threads across all blocks.
 */
__device__ void kSyncAllThreads_kernel(unsigned int* syncCounter)
{
    __syncthreads();
    if (threadIdx.x == 0)
        atomicInc(syncCounter, gridDim.x-1);
    __shared__ int counterValue;
    do
    {
        if (threadIdx.x == 0)
            counterValue = *syncCounter;
    } while (counterValue > 0);
}

__device__ void kSolveMatrix_kernel(int numTerms, unsigned int* syncCounter)
{
    for (int iteration = 0; iteration < numTerms; iteration++) {
        float* rhs1 = (iteration%2 == 0 ? cSim.pLincsRhs1 : cSim.pLincsRhs2);
        float* rhs2 = (iteration%2 == 0 ? cSim.pLincsRhs2 : cSim.pLincsRhs1);
        unsigned int pos = threadIdx.x + blockIdx.x * blockDim.x;
        while (pos < cSim.lincsConstraints)
        {
            float rhs = 0.0f;
            int start = cSim.pLincsConnectionsIndex[pos];
            int end = cSim.pLincsConnectionsIndex[pos+1];
            for (int i = start; i < end; i++)
            {
                int otherConstraint = cSim.pLincsConnections[i];
                rhs += cSim.pLincsCoupling[i]*rhs1[otherConstraint];
            }
            rhs2[pos] = rhs;
            cSim.pLincsSolution[pos] += rhs;
            pos += blockDim.x * gridDim.x;
        }
        kSyncAllThreads_kernel(&syncCounter[iteration]);
    }
}

__device__ void kUpdateAtomPositions_kernel(float4* atomPositions)
{
    unsigned int pos = threadIdx.x + blockIdx.x * blockDim.x;
    while (pos < cSim.atoms)
    {
        float4 atomPos = atomPositions[pos];
        float invMass = cSim.pVelm4[pos].w;
        int start = cSim.pLincsAtomConstraintsIndex[pos];
        int end = cSim.pLincsAtomConstraintsIndex[pos+1];
        for (int i = start; i < end; i++)
        {
            int constraint = cSim.pLincsAtomConstraints[i];
            float4 dir = cSim.pLincsDistance[constraint];
            float c = invMass*cSim.pLincsS[constraint]*cSim.pLincsSolution[constraint];
            c = (cSim.pLincsAtoms[constraint].x == pos ? -c : c);
            atomPos.x += c*dir.x;
            atomPos.y += c*dir.y;
            atomPos.z += c*dir.z;
        }
        atomPositions[pos] = atomPos;
        pos += blockDim.x * gridDim.x;
    }
}

__global__ void kApplyLincs_kernel(int numTerms, float4* atomPositions, bool addOldPosition)
{
   // Calculate the direction of each constraint, along with the initial RHS and solution vectors.

    unsigned int pos = threadIdx.x + blockIdx.x * blockDim.x;
    while (pos < cSim.lincsConstraints)
    {
        int2 atoms = cSim.pLincsAtoms[pos];
        float4 delta1 = atomPositions[atoms.x];
        float4 delta2 = atomPositions[atoms.y];
        float4 dir = cSim.pLincsDistance[pos];
        if (addOldPosition)
        {
            float4 oldPos1 = cSim.pOldPosq[atoms.x];
            float4 oldPos2 = cSim.pOldPosq[atoms.y];
            dir.x = (oldPos1.x-oldPos2.x)+(delta1.x-delta2.x);
            dir.y = (oldPos1.y-oldPos2.y)+(delta1.y-delta2.y);
            dir.z = (oldPos1.z-oldPos2.z)+(delta1.z-delta2.z);
        }
        else
        {
            dir.x = delta1.x-delta2.x;
            dir.y = delta1.y-delta2.y;
            dir.z = delta1.z-delta2.z;
        }
        float invLength = 1.0f/sqrt(dir.x*dir.x + dir.y*dir.y + dir.z*dir.z);
        dir.x *= invLength;
        dir.y *= invLength;
        dir.z *= invLength;
        cSim.pLincsDistance[pos] = dir;
        float diff = cSim.pLincsS[pos]*(1.0f/invLength-dir.w);
        cSim.pLincsRhs1[pos] = diff;
        cSim.pLincsSolution[pos] = diff;
        pos += blockDim.x * gridDim.x;
    }
    kSyncAllThreads_kernel(cSim.pSyncCounter);

    // Build the coupling matrix.

    pos = threadIdx.x + blockIdx.x * blockDim.x;
    while (pos < cSim.lincsConstraints)
    {
        float4 dir1 = cSim.pLincsDistance[pos];
        int2 atoms1 = cSim.pLincsAtoms[pos];
        int start = cSim.pLincsConnectionsIndex[pos];
        int end = cSim.pLincsConnectionsIndex[pos+1];
        float s = cSim.pLincsS[pos];
        float invMass = cSim.pVelm4[atoms1.x].w;
        for (int i = start; i < end; i++)
        {
            int otherConstraint = cSim.pLincsConnections[i];
            float4 dir2 = cSim.pLincsDistance[otherConstraint];
            int2 atoms2 = cSim.pLincsAtoms[otherConstraint];
            float sign = (atoms1.x == atoms2.x || atoms1.y == atoms2.y ? -1.0f : 1.0f);
            cSim.pLincsCoupling[i] = sign*invMass*s*(dir1.x*dir2.x+dir1.y*dir2.y+dir1.z*dir2.z)*cSim.pLincsS[otherConstraint]; // ***** Is this the correct mass? *****
        }
        pos += blockDim.x * gridDim.x;
    }

    // Solve the matrix equation and update the atom positions.

    kSolveMatrix_kernel(numTerms, cSim.pSyncCounter+1);
    kUpdateAtomPositions_kernel(atomPositions);

    // Correct for rotational lengthening.

    pos = threadIdx.x + blockIdx.x * blockDim.x;
    while (pos < cSim.lincsConstraints)
    {
        int2 atoms = cSim.pLincsAtoms[pos];
        float4 delta1 = atomPositions[atoms.x];
        float4 delta2 = atomPositions[atoms.y];
        float3 delta;
        if (addOldPosition)
        {
            float4 oldPos1 = cSim.pOldPosq[atoms.x];
            float4 oldPos2 = cSim.pOldPosq[atoms.y];
            delta = make_float3((oldPos1.x-oldPos2.x)+(delta1.x-delta2.x),
                                (oldPos1.y-oldPos2.y)+(delta1.y-delta2.y),
                                (oldPos1.z-oldPos2.z)+(delta1.z-delta2.z));
        }
        else
        {
            delta = make_float3(delta1.x-delta2.x, delta1.y-delta2.y, delta1.z-delta2.z);
        }
        float distance = cSim.pLincsDistance[pos].w;
        float p2 = 2.0f*distance*distance-(delta.x*delta.x+delta.y*delta.y+delta.z*delta.z);
        p2 = (p2 < 0.0f ? 0.0f : p2);
        float diff = cSim.pLincsS[pos]*(distance-sqrt(p2));
        cSim.pLincsRhs1[pos] = diff;
        cSim.pLincsSolution[pos] = diff;
        pos += blockDim.x * gridDim.x;
    }

    // Solve the matrix equation and update the atom positions.

    kSolveMatrix_kernel(numTerms, cSim.pSyncCounter+numTerms+1);
    kUpdateAtomPositions_kernel(atomPositions);
}

void printDist(float4 v1, float4 v2)
{
    float dx = v1.x-v2.x;
    float dy = v1.y-v2.y;
    float dz = v1.z-v2.z;
    printf("%f ", sqrt(dx*dx+dy*dy+dz*dz));
}
void kApplyFirstLincs(gpuContext gpu)
{
//    printf("kApplyFirstLincs\n");
    if (gpu->sim.lincsConstraints > 0)
    {
        kApplyLincs_kernel<<<gpu->sim.blocks, gpu->sim.lincs_threads_per_block>>>(4, gpu->sim.pPosqP, true);
        LAUNCHERROR("kApplyFirstLincs");
    }
}

void kApplySecondLincs(gpuContext gpu)
{
//    printf("kApplySecondLincs\n");
    if (gpu->sim.lincsConstraints > 0)
    {
        kApplyLincs_kernel<<<gpu->sim.blocks, gpu->sim.lincs_threads_per_block>>>(4, gpu->sim.pPosq, false);
        LAUNCHERROR("kApplySecondLincs");
    }
}
