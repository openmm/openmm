/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2010 Stanford University and the Authors.           *
 * Authors: Peter Eastman                                                     *
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
//#include <fstream>
using namespace std;

#include "gputypes.h"
#include "cudaKernels.h"

__global__ void kPrepareConstraints_kernel(int numAtoms, float4* oldPos, float4* posq, float4* posqP) {
    for (int index = threadIdx.x+blockIdx.x*blockDim.x; index < numAtoms; index += blockDim.x*gridDim.x) {
        float4 pos = posq[index];
        oldPos[index] = pos;
        posqP[index] = make_float4(0.0f, 0.0f, 0.0f, pos.w);
    }
}

__global__ void kFinishConstraints_kernel(int numAtoms, float4* posq, float4* posqP) {
    for (int index = threadIdx.x+blockIdx.x*blockDim.x; index < numAtoms; index += blockDim.x*gridDim.x) {
        float4 pos = posq[index];
        float4 delta = posqP[index];
        posq[index] = make_float4(pos.x+delta.x, pos.y+delta.y, pos.z+delta.z, pos.w);
    }
}

void kApplyConstraints(gpuContext gpu)
{
    kPrepareConstraints_kernel<<<gpu->sim.blocks, gpu->sim.update_threads_per_block>>>(gpu->natoms, gpu->sim.pOldPosq, gpu->sim.pPosq, gpu->sim.pPosqP);
    LAUNCHERROR("kPrepareConstraints");
    kApplyFirstShake(gpu);
    kApplyFirstSettle(gpu);
    kApplyFirstCCMA(gpu);
    kFinishConstraints_kernel<<<gpu->sim.blocks, gpu->sim.update_threads_per_block>>>(gpu->natoms, gpu->sim.pPosq, gpu->sim.pPosqP);
    LAUNCHERROR("kFinishConstraints");
}

