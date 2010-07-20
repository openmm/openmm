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
#include <fstream>
using namespace std;

#include "gputypes.h"
#include "cudatypes.h"


#define LOCAL_HACK_PI 3.1415926535897932384626433832795

#define DOT3(v1, v2) (v1.x*v2.x + v1.y*v2.y + v1.z*v2.z)

#define CROSS_PRODUCT(v1, v2) make_float3(v1.y*v2.z - v1.z*v2.y, v1.z*v2.x - v1.x*v2.z, v1.x*v2.y - v1.y*v2.x)

#define GETNORMEDDOTPRODUCT(v1, v2, dp) \
{ \
    dp          = DOT3(v1, v2); \
    float norm1 = DOT3(v1, v1); \
    float norm2 = DOT3(v2, v2); \
    dp /= sqrt(norm1 * norm2); \
    dp = min(dp, 1.0f); \
    dp = max(dp, -1.0f); \
}

#define GETANGLEBETWEENTWOVECTORS(v1, v2, angle) \
{ \
    float dp; \
    GETNORMEDDOTPRODUCT(v1, v2, dp); \
    if (dp > 0.99f || dp < -0.99f) { \
        float3 cross = CROSS_PRODUCT(v1, v2); \
        float scale = DOT3(v1, v1)*DOT3(v2, v2); \
        angle = asin(sqrt(DOT3(cross, cross)/scale)); \
        if (dp < 0.0f) \
            angle = LOCAL_HACK_PI-angle; \
    } \
    else { \
        angle = acos(dp); \
    } \
}

#define GETDIHEDRALANGLEBETWEENTHREEVECTORS(vector1, vector2, vector3, signVector, cp0, cp1, angle) \
{ \
    cp0 = CROSS_PRODUCT(vector1, vector2); \
    cp1 = CROSS_PRODUCT(vector2, vector3); \
    GETANGLEBETWEENTWOVECTORS(cp0, cp1, angle); \
    float dp = DOT3(signVector, cp1); \
    angle = (dp >= 0) ? angle : -angle; \
}

__global__
#if (__CUDA_ARCH__ >= 200)
__launch_bounds__(GF1XX_LOCALFORCES_THREADS_PER_BLOCK, 1)
#elif (__CUDA_ARCH__ >= 130)
__launch_bounds__(GT2XX_LOCALFORCES_THREADS_PER_BLOCK, 1)
#else
__launch_bounds__(G8X_LOCALFORCES_THREADS_PER_BLOCK, 1)
#endif
void kCalculateCMAPTorsionForces_kernel(int numAtoms, int numTorsions, float4* forceBuffers, float* energyBuffer,
        float4* posq, float4* coeff, int2* mapPositions, int4* indices, int* maps)
{
    const float PI = 3.14159265358979323846f;
    float energy = 0.0f;
    for (int index = blockIdx.x*blockDim.x+threadIdx.x; index < numTorsions; index += blockDim.x*gridDim.x) {
        int4 atoms1 = indices[4*index];
        int4 atoms2 = indices[4*index+1];
        int4 atoms3 = indices[4*index+2];
        int4 atoms4 = indices[4*index+3];
        float4 a1 = posq[atoms1.x];
        float4 a2 = posq[atoms1.y];
        float4 a3 = posq[atoms1.z];
        float4 a4 = posq[atoms1.w];
        float4 b1 = posq[atoms2.x];
        float4 b2 = posq[atoms2.y];
        float4 b3 = posq[atoms2.z];
        float4 b4 = posq[atoms2.w];

        // Compute the first angle.

        float3 v0a = make_float3(a1.x-a2.x, a1.y-a2.y, a1.z-a2.z);
        float3 v1a = make_float3(a3.x-a2.x, a3.y-a2.y, a3.z-a2.z);
        float3 v2a = make_float3(a3.x-a4.x, a3.y-a4.y, a3.z-a4.z);
        float3 cp0a, cp1a;
        float angleA;
        GETDIHEDRALANGLEBETWEENTHREEVECTORS(v0a, v1a, v2a, v0a, cp0a, cp1a, angleA);
        angleA = fmod(angleA+2.0f*PI, 2.0f*PI);

        // Compute the second angle.

        float3 v0b = make_float3(b1.x-b2.x, b1.y-b2.y, b1.z-b2.z);
        float3 v1b = make_float3(b3.x-b2.x, b3.y-b2.y, b3.z-b2.z);
        float3 v2b = make_float3(b3.x-b4.x, b3.y-b4.y, b3.z-b4.z);
        float3 cp0b, cp1b;
        float angleB;
        GETDIHEDRALANGLEBETWEENTHREEVECTORS(v0b, v1b, v2b, v0b, cp0b, cp1b, angleB);
        angleB = fmod(angleB+2.0f*PI, 2.0f*PI);

        // Identify which patch this is in.

        int2 pos = mapPositions[maps[index]];
        int size = pos.y;
        float delta = 2.0f*PI/size;
        int s = (int) (angleA/delta);
        int t = (int) (angleB/delta);
        float4 c[4];
        int coeffIndex = 4*(pos.x+s+size*t);
        c[0] = coeff[coeffIndex];
        c[1] = coeff[coeffIndex+1];
        c[2] = coeff[coeffIndex+2];
        c[3] = coeff[coeffIndex+3];
        float da = angleA/delta-s;
        float db = angleB/delta-t;

        // Evaluate the spline to determine the energy and gradients.

        float torsionEnergy = 0.0f;
        float dEdA = 0.0f;
        float dEdB = 0.0f;
        torsionEnergy = da*torsionEnergy + ((c[3].w*db + c[3].z)*db + c[3].y)*db + c[3].x;
        dEdA = db*dEdA + (3.0f*c[3].w*da + 2.0f*c[2].w)*da + c[1].w;
        dEdB = da*dEdB + (3.0f*c[3].w*db + 2.0f*c[3].z)*db + c[3].y;
        torsionEnergy = da*torsionEnergy + ((c[2].w*db + c[2].z)*db + c[2].y)*db + c[2].x;
        dEdA = db*dEdA + (3.0f*c[3].z*da + 2.0f*c[2].z)*da + c[1].z;
        dEdB = da*dEdB + (3.0f*c[2].w*db + 2.0f*c[2].z)*db + c[2].y;
        torsionEnergy = da*torsionEnergy + ((c[1].w*db + c[1].z)*db + c[1].y)*db + c[1].x;
        dEdA = db*dEdA + (3.0f*c[3].y*da + 2.0f*c[2].y)*da + c[1].y;
        dEdB = da*dEdB + (3.0f*c[1].w*db + 2.0f*c[1].z)*db + c[1].y;
        torsionEnergy = da*torsionEnergy + ((c[0].w*db + c[0].z)*db + c[0].y)*db + c[0].x;
        dEdA = db*dEdA + (3.0f*c[3].x*da + 2.0f*c[2].x)*da + c[1].x;
        dEdB = da*dEdB + (3.0f*c[0].w*db + 2.0f*c[0].z)*db + c[0].y;
        dEdA /= delta;
        dEdB /= delta;
        energy += torsionEnergy;

        // Apply the force to the first torsion.

        float normCross1 = DOT3(cp0a, cp0a);
        float normSqrBC = DOT3(v1a, v1a);
        float normBC = sqrt(normSqrBC);
        float normCross2 = DOT3(cp1a, cp1a);
        float dp = 1.0f/normSqrBC;
        float4 ff = make_float4((-dEdA*normBC)/normCross1, DOT3(v0a, v1a)*dp, DOT3(v2a, v1a)*dp, (dEdA*normBC)/normCross2);
        float3 internalF0 = make_float3(ff.x*cp0a.x, ff.x*cp0a.y, ff.x*cp0a.z);
        float3 internalF3 = make_float3(ff.w*cp1a.x, ff.w*cp1a.y, ff.w*cp1a.z);
        float3 d = make_float3(ff.y*internalF0.x - ff.z*internalF3.x,
                               ff.y*internalF0.y - ff.z*internalF3.y,
                               ff.y*internalF0.z - ff.z*internalF3.z);
        unsigned int offsetA = atoms1.x+atoms3.x*numAtoms;
        unsigned int offsetB = atoms1.y+atoms3.y*numAtoms;
        unsigned int offsetC = atoms1.z+atoms3.z*numAtoms;
        unsigned int offsetD = atoms1.w+atoms3.w*numAtoms;
        float4 forceA = forceBuffers[offsetA];
        float4 forceB = forceBuffers[offsetB];
        float4 forceC = forceBuffers[offsetC];
        float4 forceD = forceBuffers[offsetD];
        forceA.x += internalF0.x;
        forceA.y += internalF0.y;
        forceA.z += internalF0.z;
        forceB.x += d.x-internalF0.x;
        forceB.y += d.y-internalF0.y;
        forceB.z += d.z-internalF0.z;
        forceC.x += -d.x-internalF3.x;
        forceC.y += -d.y-internalF3.y;
        forceC.z += -d.z-internalF3.z;
        forceD.x += internalF3.x;
        forceD.y += internalF3.y;
        forceD.z += internalF3.z;
        forceBuffers[offsetA] = forceA;
        forceBuffers[offsetB] = forceB;
        forceBuffers[offsetC] = forceC;
        forceBuffers[offsetD] = forceD;

        // Apply the force to the second torsion.

        normCross1 = DOT3(cp0b, cp0b);
        normSqrBC = DOT3(v1b, v1b);
        normBC = sqrt(normSqrBC);
        normCross2 = DOT3(cp1b, cp1b);
        dp = 1.0f/normSqrBC;
        ff = make_float4((-dEdB*normBC)/normCross1, DOT3(v0b, v1b)*dp, DOT3(v2b, v1b)*dp, (dEdB*normBC)/normCross2);
        internalF0 = make_float3(ff.x*cp0b.x, ff.x*cp0b.y, ff.x*cp0b.z);
        internalF3 = make_float3(ff.w*cp1b.x, ff.w*cp1b.y, ff.w*cp1b.z);
        d = make_float3(ff.y*internalF0.x - ff.z*internalF3.x,
                        ff.y*internalF0.y - ff.z*internalF3.y,
                        ff.y*internalF0.z - ff.z*internalF3.z);
        offsetA = atoms2.x+atoms4.x*numAtoms;
        offsetB = atoms2.y+atoms4.y*numAtoms;
        offsetC = atoms2.z+atoms4.z*numAtoms;
        offsetD = atoms2.w+atoms4.w*numAtoms;
        forceA = forceBuffers[offsetA];
        forceB = forceBuffers[offsetB];
        forceC = forceBuffers[offsetC];
        forceD = forceBuffers[offsetD];
        forceA.x += internalF0.x;
        forceA.y += internalF0.y;
        forceA.z += internalF0.z;
        forceB.x += d.x-internalF0.x;
        forceB.y += d.y-internalF0.y;
        forceB.z += d.z-internalF0.z;
        forceC.x += -d.x-internalF3.x;
        forceC.y += -d.y-internalF3.y;
        forceC.z += -d.z-internalF3.z;
        forceD.x += internalF3.x;
        forceD.y += internalF3.y;
        forceD.z += internalF3.z;
        forceBuffers[offsetA] = forceA;
        forceBuffers[offsetB] = forceB;
        forceBuffers[offsetC] = forceC;
        forceBuffers[offsetD] = forceD;
    }
    energyBuffer[blockIdx.x*blockDim.x+threadIdx.x] += energy;
}

void kCalculateCMAPTorsionForces(gpuContext gpu, CUDAStream<float4>& coefficients, CUDAStream<int2>& mapPositions, CUDAStream<int4>& torsionIndices, CUDAStream<int>& torsionMaps)
{
    kCalculateCMAPTorsionForces_kernel<<<gpu->sim.blocks, gpu->sim.localForces_threads_per_block>>>(gpu->sim.stride,
            torsionMaps._length, gpu->sim.pForce4, gpu->sim.pEnergy, gpu->sim.pPosq, coefficients._pDevData,
            mapPositions._pDevData, torsionIndices._pDevData, torsionMaps._pDevData);
    LAUNCHERROR("kCalculateCMAPTorsionForces");
}
