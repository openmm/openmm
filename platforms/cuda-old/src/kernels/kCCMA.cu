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

void SetCCMASim(gpuContext gpu)
{
    cudaError_t status;
    status = cudaMemcpyToSymbol(cSim, &gpu->sim, sizeof(cudaGmxSimulation));
    RTERROR(status, "cudaMemcpyToSymbol: SetSim copy to cSim failed");
}

void GetCCMASim(gpuContext gpu)
{
    cudaError_t status;
    status = cudaMemcpyFromSymbol(&gpu->sim, cSim, sizeof(cudaGmxSimulation));
    RTERROR(status, "cudaMemcpyFromSymbol: SetSim copy from cSim failed");
}

__global__ void
#if (__CUDA_ARCH__ >= 200)
__launch_bounds__(1024, 1)
#elif (__CUDA_ARCH__ >= 120)
__launch_bounds__(512, 1)
#else
__launch_bounds__(256, 1)
#endif
kComputeCCMAConstraintDirections()
{
    // Calculate the direction of each constraint.

    for (unsigned int index = threadIdx.x+blockIdx.x*blockDim.x; index < cSim.ccmaConstraints; index += blockDim.x*gridDim.x)
    {
        int2 atoms = cSim.pCcmaAtoms[index];
        float4 dir = cSim.pCcmaDistance[index];
        float4 oldPos1 = cSim.pOldPosq[atoms.x];
        float4 oldPos2 = cSim.pOldPosq[atoms.y];
        dir.x = oldPos1.x-oldPos2.x;
        dir.y = oldPos1.y-oldPos2.y;
        dir.z = oldPos1.z-oldPos2.z;
        cSim.pCcmaDistance[index] = dir;
    }
}

__global__ void
#if (__CUDA_ARCH__ >= 200)
__launch_bounds__(1024, 1)
#elif (__CUDA_ARCH__ >= 120)
__launch_bounds__(512, 1)
#else
__launch_bounds__(256, 1)
#endif
kComputeCCMAConstraintForces(float4* atomPositions, bool addOldPosition)
{
    __shared__ int converged;
    float lowerTol = 1.0f-2.0f*cSim.shakeTolerance+cSim.shakeTolerance*cSim.shakeTolerance;
    float upperTol = 1.0f+2.0f*cSim.shakeTolerance+cSim.shakeTolerance*cSim.shakeTolerance;
    if (threadIdx.x == 0)
        converged = 1;
    __syncthreads();

    // Calculate the constraint force for each constraint.

    for (unsigned int index = threadIdx.x+blockIdx.x*blockDim.x; index < cSim.ccmaConstraints; index += blockDim.x*gridDim.x)
    {
        int2 atoms = cSim.pCcmaAtoms[index];
        float4 delta1 = atomPositions[atoms.x];
        float4 delta2 = atomPositions[atoms.y];
        float4 dir = cSim.pCcmaDistance[index];
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
        float reducedMass = cSim.pCcmaReducedMass[index];
        cSim.pCcmaDelta1[index] = (rrpr > d_ij2*1e-6f ? reducedMass*diff/rrpr : 0.0f);

        // See whether it has converged.

        if (converged && (rp2 < lowerTol*dist2 || rp2 > upperTol*dist2))
        {
            converged = 0;
            *cSim.ccmaConvergedDeviceMarker = 0;
        }
    }
}

__global__ void
#if (__CUDA_ARCH__ >= 200)
__launch_bounds__(1024, 1)
#elif (__CUDA_ARCH__ >= 120)
__launch_bounds__(512, 1)
#else
__launch_bounds__(256, 1)
#endif
kMultiplyByCCMAConstraintMatrix()
{
    if (*cSim.ccmaConvergedDeviceMarker)
        return; // The constraint iteration has already converged

    // Multiply by the inverse constraint matrix.

    for (unsigned int index = threadIdx.x+blockIdx.x*blockDim.x; index < cSim.ccmaConstraints; index += blockDim.x*gridDim.x)
    {
        float sum = 0.0f;
        for (unsigned int i = 0; ; i++)
        {
            unsigned int element = index+i*cSim.ccmaConstraints;
            unsigned int column = cSim.pConstraintMatrixColumn[element];
            if (column >= cSim.ccmaConstraints)
                break;
            sum += cSim.pCcmaDelta1[column]*cSim.pConstraintMatrixValue[element];
        }
        cSim.pCcmaDelta2[index] = sum;
    }
}

__global__ void
#if (__CUDA_ARCH__ >= 200)
__launch_bounds__(1024, 1)
#elif (__CUDA_ARCH__ >= 120)
__launch_bounds__(512, 1)
#else
__launch_bounds__(256, 1)
#endif
kUpdateCCMAAtomPositions(float4* atomPositions, int iteration)
{
    if (*cSim.ccmaConvergedDeviceMarker)
        return; // The constraint iteration has already converged.
    float damping = (iteration < 2 ? 0.5f : 1.0f);
    for (unsigned int index = threadIdx.x+blockIdx.x*blockDim.x; index < cSim.atoms; index += blockDim.x*gridDim.x)
    {
        float4 atomPos = atomPositions[index];
        float invMass = cSim.pVelm4[index].w;
        int num = cSim.pCcmaNumAtomConstraints[index];
        for (int i = 0; i < num; i++)
        {
            int constraint = cSim.pCcmaAtomConstraints[index+i*cSim.atoms];
            bool forward = (constraint > 0);
            constraint = (forward ? constraint-1 : -constraint-1);
            float constraintForce = damping*invMass*cSim.pCcmaDelta2[constraint];
            constraintForce = (forward ? constraintForce : -constraintForce);
            float4 dir = cSim.pCcmaDistance[constraint];
            atomPos.x += constraintForce*dir.x;
            atomPos.y += constraintForce*dir.y;
            atomPos.z += constraintForce*dir.z;
        }
        atomPositions[index] = atomPos;
    }
}

void kApplyCCMA(gpuContext gpu, float4* posq, bool addOldPosition)
{
    kComputeCCMAConstraintDirections<<<gpu->sim.blocks, gpu->sim.ccma_threads_per_block>>>();
    LAUNCHERROR("kComputeCCMAConstraintDirections");
    const int checkInterval = 3;
    for (int i = 0; i < 150; i++) {
        if ((i+1)%checkInterval == 0)
            *gpu->ccmaConvergedHostMarker = 1;
        kComputeCCMAConstraintForces<<<gpu->sim.blocks, gpu->sim.ccma_threads_per_block, gpu->sim.ccma_threads_per_block*sizeof(int)>>>(posq, addOldPosition);
        cudaEventRecord(gpu->ccmaEvent, 0);
        kMultiplyByCCMAConstraintMatrix<<<gpu->sim.blocks, gpu->sim.ccma_threads_per_block, gpu->sim.ccma_threads_per_block*sizeof(int)>>>();
        kUpdateCCMAAtomPositions<<<gpu->sim.blocks, gpu->sim.ccma_threads_per_block>>>(posq, 3*i+2);
        cudaEventSynchronize(gpu->ccmaEvent);
        if ((i+1)%checkInterval == 0 && *gpu->ccmaConvergedHostMarker)
            break;
    }
}

void kApplyCCMA(gpuContext gpu)
{
//    printf("kApplyCCMA\n");
    if (gpu->sim.ccmaConstraints > 0)
        kApplyCCMA(gpu, gpu->sim.pPosqP, true);
}
