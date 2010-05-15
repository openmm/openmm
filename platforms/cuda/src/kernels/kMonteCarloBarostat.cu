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

__global__ void kScaleAtomCoordinates_kernel(float scale, int numMolecules, float3 periodicBoxSize, float4* posq, int* moleculeAtoms, int* moleculeStartIndex) {
    float3 invPeriodicBoxSize = make_float3(1.0f/periodicBoxSize.x, 1.0f/periodicBoxSize.y, 1.0f/periodicBoxSize.z);
    for (int index = threadIdx.x+blockIdx.x*blockDim.x; index < numMolecules; index += blockDim.x*gridDim.x) {
        int first = moleculeStartIndex[index];
        int last = moleculeStartIndex[index+1];
        int numAtoms = last-first;

        // Find the center of each molecule.

        float3 center = make_float3(0, 0, 0);
        for (int atom = first; atom < last; atom++) {
            float4 pos = posq[moleculeAtoms[atom]];
            center.x += pos.x;
            center.y += pos.y;
            center.z += pos.z;
        }
        center.x /= (float) numAtoms;
        center.y /= (float) numAtoms;
        center.z /= (float) numAtoms;

        // Move it into the first periodic box.

        int xcell = (int) floor(center.x*invPeriodicBoxSize.x);
        int ycell = (int) floor(center.y*invPeriodicBoxSize.y);
        int zcell = (int) floor(center.z*invPeriodicBoxSize.z);
        float3 delta = make_float3(xcell*periodicBoxSize.x, ycell*periodicBoxSize.y, zcell*periodicBoxSize.z);
        center.x -= delta.x;
        center.y -= delta.y;
        center.z -= delta.z;

        // Now scale the position of the molecule center.

        delta.x = center.x*(scale-1)-delta.x;
        delta.y = center.y*(scale-1)-delta.y;
        delta.z = center.z*(scale-1)-delta.z;
        for (int atom = first; atom < last; atom++) {
            float4 pos = posq[moleculeAtoms[atom]];
            pos.x += delta.x;
            pos.y += delta.y;
            pos.z += delta.z;
            posq[moleculeAtoms[atom]] = pos;
        }
    }
}

void kScaleAtomCoordinates(gpuContext gpu, float scale, CUDAStream<int>& moleculeAtoms, CUDAStream<int>& moleculeStartIndex)
{
//    printf("kScaleAtomCoordinates\n");
    kScaleAtomCoordinates_kernel<<<gpu->sim.blocks, gpu->sim.update_threads_per_block>>>(scale, moleculeStartIndex._length-1,
            make_float3(gpu->sim.periodicBoxSizeX, gpu->sim.periodicBoxSizeY, gpu->sim.periodicBoxSizeZ), gpu->sim.pPosq,
            moleculeAtoms._pDevData, moleculeStartIndex._pDevData);
    LAUNCHERROR("kScaleAtomCoordinates");
}

