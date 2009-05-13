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

#include <cuda.h>
#include <vector_functions.h>
#include <vector>
#include "jama_svd.h"
#include "gputypes.h"

using namespace std;
using TNT::Array2D;
using JAMA::SVD;


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
                unsigned int firstIndex = cSim.pRigidClusterConstraintIndex[block];
                unsigned int blockSize = cSim.pRigidClusterConstraintIndex[block+1]-firstIndex;
                if (indexInBlock < blockSize)
                {
                    // Load the constraint forces and matrix.

                    unsigned int constraint = cSim.pRigidClusterConstraints[firstIndex+indexInBlock];
                    temp[threadIdx.x] = cSim.pLincsSolution[constraint];
                    unsigned int firstMatrixIndex = cSim.pRigidClusterMatrixIndex[block];

                    // Multiply by the matrix.

                    float sum = 0.0f;
                    for (unsigned int i = 0; i < blockSize; i++)
                        sum += temp[threadIdx.x-indexInBlock+i]*cSim.pRigidClusterMatrix[firstMatrixIndex+i*blockSize+indexInBlock];
                    cSim.pLincsSolution[constraint] = sum;
                }
                block += (blockDim.x*gridDim.x)/cSim.clusterShakeBlockSize;
            }
            kSyncAllThreads_kernel(&cSim.pSyncCounter[gridDim.x], iteration);
        }

        // Update the position of each atom.

        pos = threadIdx.x + blockIdx.x * blockDim.x;
        while (pos < cSim.atoms)
        {
            float4 atomPos = atomPositions[pos];
            float invMass = cSim.pVelm4[pos].w;
            int num = cSim.pLincsNumAtomConstraints[pos];
            for (int i = 0; i < num; i++)
            {
                int index = pos+i*cSim.atoms;
                int constraint = cSim.pLincsAtomConstraints[index];
                float constraintForce = invMass*cSim.pLincsSolution[constraint];
                constraintForce = (cSim.pLincsAtoms[constraint].x == pos ? constraintForce : -constraintForce);
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

static void initInverseMatrices(gpuContext gpu)
{
    // Build the inverse constraint matrix for each cluster.

    gpu->psPosq4->Download();
    gpu->psVelm4->Download();
    unsigned int elementIndex = 0;
    for (unsigned int i = 0; i < gpu->sim.rigidClusters; i++) {
        // Compute the constraint coupling matrix for this cluster.

        unsigned int startIndex = (*gpu->psRigidClusterConstraintIndex)[i];
        unsigned int endIndex = (*gpu->psRigidClusterConstraintIndex)[i+1];
        unsigned int size = endIndex-startIndex;
        vector<float3> r(size);
        for (unsigned int j = 0; j < size; j++) {
            int2 atoms = (*gpu->psLincsAtoms)[(*gpu->psRigidClusterConstraints)[startIndex+j]];
            float4 pos1 = (*gpu->psPosq4)[atoms.x];
            float4 pos2 = (*gpu->psPosq4)[atoms.y];
            r[j] = make_float3(pos1.x-pos2.x, pos1.y-pos2.y, pos1.z-pos2.z);
            float invLength = 1.0f/sqrt(r[j].x*r[j].x + r[j].y*r[j].y + r[j].z*r[j].z);
            r[j].x *= invLength;
            r[j].y *= invLength;
            r[j].z *= invLength;
        }
        Array2D<double> matrix(size, size);
        for (unsigned int j = 0; j < size; j++) {
            int constraintj = (*gpu->psRigidClusterConstraints)[startIndex+j];
            int2 atomsj = (*gpu->psLincsAtoms)[constraintj];
            for (unsigned int k = 0; k < size; k++) {
                int constraintk = (*gpu->psRigidClusterConstraints)[startIndex+k];
                int2 atomsk = (*gpu->psLincsAtoms)[constraintk];
                float invMassj0 = (*gpu->psVelm4)[atomsj.x].w;
                float invMassj1 = (*gpu->psVelm4)[atomsj.y].w;
                double dot = r[j].x*r[k].x + r[j].y*r[k].y + r[j].z*r[k].z;
                if (atomsj.x == atomsk.x)
                    dot *= invMassj0/(invMassj0+invMassj1);
                else if (atomsj.y == atomsk.y)
                    dot *= invMassj1/(invMassj0+invMassj1);
                else if (atomsj.x == atomsk.y)
                    dot *= -invMassj0/(invMassj0+invMassj1);
                else if (atomsj.y == atomsk.x)
                    dot *= -invMassj1/(invMassj0+invMassj1);
                else
                    dot = 0.0;
                matrix[j][k] = dot;
            }
            matrix[j][j] = 1.0;
        }

        // Invert it using SVD.

        Array2D<double> u, v;
        Array1D<double> w;
        SVD<double> svd(matrix);
        svd.getU(u);
        svd.getV(v);
        svd.getSingularValues(w);
        double singularValueCutoff = 0.01*w[0];
        for (unsigned int j = 0; j < size; j++)
            w[j] = (w[j] < singularValueCutoff ? 0.0 : 1.0/w[j]);
        for (unsigned int j = 0; j < size; j++) {
            for (unsigned int k = 0; k < size; k++) {
                matrix[j][k] = 0.0;
                for (unsigned int m = 0; m < size; m++)
                    matrix[j][k] += v[j][m]*w[m]*u[k][m];
            }
        }

        // Record the inverted matrix.

        (*gpu->psRigidClusterMatrixIndex)[i] = elementIndex;
        for (unsigned int j = 0; j < size; j++)
        {
            float distance1 = (*gpu->psLincsDistance)[(*gpu->psRigidClusterConstraints)[startIndex+j]].w;
            for (unsigned int k = 0; k < size; k++)
            {
                float distance2 = (*gpu->psLincsDistance)[(*gpu->psRigidClusterConstraints)[startIndex+k]].w;
                (*gpu->psRigidClusterMatrix)[elementIndex++] = matrix[k][j]*distance1/distance2;
            }
        }
    }
    (*gpu->psRigidClusterMatrixIndex)[gpu->sim.rigidClusters] = elementIndex;
    gpu->psRigidClusterMatrix->Upload();
    gpu->psRigidClusterMatrixIndex->Upload();
    gpu->hasInitializedRigidClusters = true;
}

void kApplyFirstCShake(gpuContext gpu)
{
//    printf("kApplyFirstCShake\n");
    if (gpu->sim.lincsConstraints > 0)
    {
        if (!gpu->hasInitializedRigidClusters)
            initInverseMatrices(gpu);
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
