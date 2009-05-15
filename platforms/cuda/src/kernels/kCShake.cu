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
#include <vector>
#include "gputypes.h"

using namespace std;


static __constant__ cudaGmxSimulation cSim;

void SetCShakeSim(gpuContext gpu)
{
    cudaError_t status;
    status = cudaMemcpyToSymbol(cSim, &gpu->sim, sizeof(cudaGmxSimulation));
    RTERROR(status, "cudaMemcpyToSymbol: SetSim copy to cSim failed");
}

void GetCShakeSim(gpuContext gpu)
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

__global__ void kApplyCShake_kernel(float4* atomPositions, bool addOldPosition)
{
    extern __shared__ float temp[];

    // Initialize counters used for monitoring convergence and doing global thread synchronization.

    __shared__ unsigned int requiredIterations;
    if (threadIdx.x == 0)
    {
        requiredIterations = 0;
        cSim.pSyncCounter[gridDim.x+blockIdx.x] = -1;
        cSim.pSyncCounter[2*gridDim.x+blockIdx.x] = -1;
        cSim.pRequiredIterations[0] = 0;
    }

    // Calculate the direction of each constraint.

    unsigned int pos = threadIdx.x + blockIdx.x * blockDim.x;
    while (pos < cSim.lincsConstraints)
    {
        int2 atoms = cSim.pLincsAtoms[pos];
        float4 dir = cSim.pLincsDistance[pos];
        float4 oldPos1 = cSim.pOldPosq[atoms.x];
        float4 oldPos2 = cSim.pOldPosq[atoms.y];
        dir.x = oldPos1.x-oldPos2.x;
        dir.y = oldPos1.y-oldPos2.y;
        dir.z = oldPos1.z-oldPos2.z;
        cSim.pLincsDistance[pos] = dir;
        pos += blockDim.x*gridDim.x;
    }
    __syncthreads();

    // Iteratively update the atom positions.

    unsigned int maxIterations = 150;
    float lowerTol = 1.0f-2.0f*cSim.shakeTolerance+cSim.shakeTolerance*cSim.shakeTolerance;
    float upperTol = 1.0f+2.0f*cSim.shakeTolerance+cSim.shakeTolerance*cSim.shakeTolerance;
    for (unsigned int iteration = 0; iteration < maxIterations && iteration == requiredIterations; iteration++)
    {
        // Calculate the constraint force for each constraint.

        pos = threadIdx.x + blockIdx.x * blockDim.x;
        while (pos < cSim.lincsConstraints)
        {
            int2 atoms = cSim.pLincsAtoms[pos];
            float4 delta1 = atomPositions[atoms.x];
            float4 delta2 = atomPositions[atoms.y];
            float4 dir = cSim.pLincsDistance[pos];
            float3 rp_ij = make_float3(delta1.x-delta2.x, delta1.y-delta2.y, delta1.z-delta2.z);
            if (addOldPosition)
            {
                rp_ij.x += dir.x;
                rp_ij.y += dir.y;
                rp_ij.z += dir.z;
            }
            float rp2 = rp_ij.x*rp_ij.x + rp_ij.y*rp_ij.y + rp_ij.z*rp_ij.z;
            float dist2 = dir.w*dir.w;
            float diff = dist2 - rp2;
            float rrpr  = rp_ij.x*dir.x + rp_ij.y*dir.y + rp_ij.z*dir.z;
            float d_ij2  = dir.x*dir.x + dir.y*dir.y + dir.z*dir.z;
            float reducedMass = cSim.pShakeReducedMass[pos];
            cSim.pLincsSolution[pos] = (rrpr > d_ij2*1e-6f ? reducedMass*diff/rrpr : 0.0f);
            if (requiredIterations == iteration && (rp2 < lowerTol*dist2 || rp2 > upperTol*dist2))
                requiredIterations = iteration+1;
            pos += blockDim.x * gridDim.x;
        }
        kSyncAllThreads_kernel(cSim.pSyncCounter, iteration);
        if (threadIdx.x == 0 && requiredIterations > iteration)
            cSim.pRequiredIterations[0] = requiredIterations;

        // Multiply by the inverse constraint matrix for each rigid cluster.

        if (cSim.rigidClusters > 0)
        {
            pos = threadIdx.x + blockIdx.x * blockDim.x;
            unsigned int block = pos/cSim.clusterShakeBlockSize;
            unsigned int indexInBlock = pos-block*cSim.clusterShakeBlockSize;
            while (block < cSim.rigidClusters)
            {
                unsigned int firstConstraint = cSim.pRigidClusterConstraintIndex[block];
                unsigned int blockSize = cSim.pRigidClusterConstraintIndex[block+1]-firstConstraint;
                if (indexInBlock < blockSize)
                {
                    // Load the constraint forces and matrix.

                    temp[threadIdx.x] = cSim.pLincsSolution[firstConstraint+indexInBlock];
                    unsigned int firstMatrixIndex = cSim.pRigidClusterMatrixIndex[block];

                    // Multiply by the matrix.

                    float sum = 0.0f;
                    for (unsigned int i = 0; i < blockSize; i++)
                        sum += temp[threadIdx.x-indexInBlock+i]*cSim.pRigidClusterMatrix[firstMatrixIndex+i*blockSize+indexInBlock];
                    cSim.pLincsSolution[firstConstraint+indexInBlock] = sum;
                }
                block += (blockDim.x*gridDim.x)/cSim.clusterShakeBlockSize;
            }
            kSyncAllThreads_kernel(&cSim.pSyncCounter[gridDim.x], iteration);
        }

        // Update the position of each atom.

        pos = threadIdx.x + blockIdx.x * blockDim.x;
        float damping = (iteration < 2 ? 0.5f : 1.0f);
        while (pos < cSim.atoms)
        {
            float4 atomPos = atomPositions[pos];
            float invMass = cSim.pVelm4[pos].w;
            int num = cSim.pLincsNumAtomConstraints[pos];
            for (int i = 0; i < num; i++)
            {
                int index = pos+i*cSim.atoms;
                int constraint = cSim.pLincsAtomConstraints[index];
                bool forward = (constraint > 0);
                constraint = (forward ? constraint-1 : -constraint-1);
                float constraintForce = damping*invMass*cSim.pLincsSolution[constraint];
                constraintForce = (forward ? constraintForce : -constraintForce);
                float4 dir = cSim.pLincsDistance[constraint];
                atomPos.x += constraintForce*dir.x;
                atomPos.y += constraintForce*dir.y;
                atomPos.z += constraintForce*dir.z;
            }
            atomPositions[pos] = atomPos;
            pos += blockDim.x*gridDim.x;
        }
        kSyncAllThreads_kernel(&cSim.pSyncCounter[2*gridDim.x], iteration);
        requiredIterations = cSim.pRequiredIterations[0];
    }

    // Reset the initial sync counter to be ready for the next call.

    if (threadIdx.x == 0)
        cSim.pSyncCounter[blockIdx.x] = -1;
}

void kApplyFirstCShake(gpuContext gpu)
{
//    printf("kApplyFirstCShake\n");
    if (gpu->sim.lincsConstraints > 0)
    {
        kApplyCShake_kernel<<<gpu->sim.blocks, gpu->sim.lincs_threads_per_block, 4*gpu->sim.lincs_threads_per_block>>>(gpu->sim.pPosqP, true);
        LAUNCHERROR("kApplyCShake");
    }
}

void kApplySecondCShake(gpuContext gpu)
{
//    printf("kApplySecondCShake\n");
    if (gpu->sim.lincsConstraints > 0)
    {
        kApplyCShake_kernel<<<gpu->sim.blocks, gpu->sim.lincs_threads_per_block, 4*gpu->sim.lincs_threads_per_block>>>(gpu->sim.pPosq, false);
        LAUNCHERROR("kApplyCShake");
    }
}
