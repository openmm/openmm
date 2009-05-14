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

#include <cuda.h>
#include <vector_functions.h>
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
__device__ void kSyncAllThreads_kernel(short* syncCounter, short newCount)
{
//    short* syncCounter = &cSim.pSyncCounter[newCount%2 == 0 ? 0 : gridDim.x];
    __syncthreads();
    if (threadIdx.x == 0)
        syncCounter[blockIdx.x] = newCount;
    if (threadIdx.x < gridDim.x)
    {
        volatile short* counter = &syncCounter[threadIdx.x];
        do
        {
        } while (*counter != newCount);
    }
    __syncthreads();
}

__global__ void kSolveLincsMatrix_kernel(float4* atomPositions)
{
    for (unsigned int iteration = 0; iteration < cSim.lincsTerms; iteration++) {
        float* rhs1 = (iteration%2 == 0 ? cSim.pLincsRhs1 : cSim.pLincsRhs2);
        float* rhs2 = (iteration%2 == 0 ? cSim.pLincsRhs2 : cSim.pLincsRhs1);
        unsigned int pos = threadIdx.x + blockIdx.x * blockDim.x;
        while (pos < cSim.lincsConstraints)
        {
            float rhs = 0.0f;
            int num = cSim.pLincsNumConnections[pos];
            for (int i = 0; i < num; i++)
            {
                int index = pos+i*cSim.lincsConstraints;
                int otherConstraint = cSim.pLincsConnections[index];
                rhs += cSim.pLincsCoupling[index]*rhs1[otherConstraint];
            }
            rhs2[pos] = rhs;
            cSim.pLincsSolution[pos] += rhs;
            pos += blockDim.x * gridDim.x;
        }
        kSyncAllThreads_kernel(&cSim.pSyncCounter[iteration%2 == 0 ? 0 : gridDim.x], iteration);
    }

    // Update the atom positions based on the solution to the matrix equations.

    unsigned int pos = threadIdx.x + blockIdx.x * blockDim.x;
    while (pos < cSim.atoms)
    {
        float4 atomPos = atomPositions[pos];
        float invMass = cSim.pVelm4[pos].w;
        int num = cSim.pLincsNumAtomConstraints[pos];
        for (int i = 0; i < num; i++)
        {
            int index = pos+i*cSim.atoms;
            int constraint = cSim.pLincsAtomConstraints[index];
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

__global__ void kApplyLincsPart1_kernel(float4* atomPositions, bool addOldPosition)
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
    kSyncAllThreads_kernel(cSim.pSyncCounter, cSim.lincsTerms+1);

    // Build the coupling matrix.

    pos = threadIdx.x + blockIdx.x * blockDim.x;
    while (pos < cSim.lincsConstraints)
    {
        float4 dir1 = cSim.pLincsDistance[pos];
        int2 atoms1 = cSim.pLincsAtoms[pos];
        int num = cSim.pLincsNumConnections[pos];
        float s = cSim.pLincsS[pos];
        float invMass = cSim.pVelm4[atoms1.x].w;
        for (int i = 0; i < num; i++)
        {
            int index = pos+i*cSim.lincsConstraints;
            int otherConstraint = cSim.pLincsConnections[index];
            float4 dir2 = cSim.pLincsDistance[otherConstraint];
            int2 atoms2 = cSim.pLincsAtoms[otherConstraint];
            float signedMass = (atoms1.x == atoms2.x || atoms1.y == atoms2.y ? -invMass : cSim.pVelm4[atoms1.y].w);
            cSim.pLincsCoupling[index] = signedMass*s*(dir1.x*dir2.x+dir1.y*dir2.y+dir1.z*dir2.z)*cSim.pLincsS[otherConstraint];
        }
        pos += blockDim.x * gridDim.x;
    }
}

__global__ void kApplyLincsPart2_kernel(float4* atomPositions, bool addOldPosition)
{
    // Correct for rotational lengthening.

    unsigned int pos = threadIdx.x + blockIdx.x * blockDim.x;
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
}

static void kApplyLincs(gpuContext gpu, float4* atomPositions, bool addOldPosition)
{
    kApplyLincsPart1_kernel<<<gpu->sim.blocks, gpu->sim.lincs_threads_per_block>>>(atomPositions, addOldPosition);
    LAUNCHERROR("kApplyLincsPart1");
    kSolveLincsMatrix_kernel<<<gpu->sim.blocks, gpu->sim.lincs_threads_per_block>>>(atomPositions);
    LAUNCHERROR("kSolveLincsMatrix_kernel");
    kApplyLincsPart2_kernel<<<gpu->sim.blocks, gpu->sim.lincs_threads_per_block>>>(atomPositions, addOldPosition);
    LAUNCHERROR("kApplyLincsPart2");
    kSolveLincsMatrix_kernel<<<gpu->sim.blocks, gpu->sim.lincs_threads_per_block>>>(atomPositions);
    LAUNCHERROR("kSolveLincsMatrix_kernel");
}

void kApplyFirstLincs(gpuContext gpu)
{
//    printf("kApplyFirstLincs\n");
    if (gpu->sim.lincsConstraints > 0)
        kApplyLincs(gpu, gpu->sim.pPosqP, true);
}

void kApplySecondLincs(gpuContext gpu)
{
//    printf("kApplySecondLincs\n");
    if (gpu->sim.lincsConstraints > 0)
        kApplyLincs(gpu, gpu->sim.pPosq, false);
}
