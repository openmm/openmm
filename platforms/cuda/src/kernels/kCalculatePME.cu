/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2009 Stanford University and the Authors.           *
 * Authors: Erik Lindahl, Rossen Apostolov, Szilard Pall, Peter Eastman       *
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

#include "gputypes.h"
#include <cuda.h>

using namespace std;

static __constant__ cudaGmxSimulation cSim;

void SetCalculatePMESim(gpuContext gpu)
{
    cudaError_t status;
    status = cudaMemcpyToSymbol(cSim, &gpu->sim, sizeof(cudaGmxSimulation));
    RTERROR(status, "cudaMemcpyToSymbol: SetSim copy to cSim failed");
}

void GetCalculatePMESim(gpuContext gpu)
{
    cudaError_t status;
    status = cudaMemcpyFromSymbol(&gpu->sim, cSim, sizeof(cudaGmxSimulation));
    RTERROR(status, "cudaMemcpyFromSymbol: SetSim copy from cSim failed");
}

inline __host__ __device__ int fast_mod(int a, int b)
{
    return (b & (b - 1)) ? a % b : a & (b - 1);
}
inline __host__ __device__ float4 make_float4(float s)
{
    return make_float4(s, s, s, s);
}
inline __host__ __device__ float4 operator-(float4 &a)
{
    return make_float4(-a.x, -a.y, -a.z, -a.w);
}
inline __host__ __device__ float4 operator-(float4 a, float4 b)
{
    return make_float4(a.x - b.x, a.y - b.y, a.z - b.z,  a.w - b.w);
}
inline __host__ __device__ float4 operator+(float4 a, float4 b)
{
    return make_float4(a.x + b.x, a.y + b.y, a.z + b.z,  a.w + b.w);
}
inline __host__ __device__ float4 operator+(float4 a, float b)
{
    return make_float4(a.x + b, a.y + b, a.z + b, a.w + b);
}
inline __host__ __device__ float4 operator+(float a, float4 b)
{
    return make_float4(a + b.x, a + b.y, a + b.z,  a + b.w);
}
inline __host__ __device__ float4 operator*(float s, float4 a)
{
    return make_float4(a.x * s, a.y * s, a.z * s, a.w * s);
}
inline __host__ __device__ float4 operator*(float4 a, float4 b)
{
    return make_float4(a.x * b.x, a.y * b.y, a.z * b.z, a.w + b.w);
}
inline __host__ __device__ float4 make_float4(int3 a)
{
    return make_float4(a.x, a.y, a.z, 0);
}

__global__ void kUpdateGridIndexAndFraction()
{
    unsigned int tnb = blockDim.x * gridDim.x;
    unsigned int tid = blockIdx.x * blockDim.x + threadIdx.x;

    for (int i = tid; i < cSim.atoms; i += tnb)
    {
        float4 ftmp = cSim.pPosq[i];

        __syncthreads();

        float3 t = make_float3((ftmp.x/cSim.periodicBoxSizeX+1.0f)*cSim.pmeGridSize.x,
                               (ftmp.y/cSim.periodicBoxSizeY+1.0f)*cSim.pmeGridSize.y,
                               (ftmp.z/cSim.periodicBoxSizeZ+1.0f)*cSim.pmeGridSize.z);
        float3 tix;
        ftmp.x = modff(t.x, &tix.x);
        ftmp.y = modff(t.y, &tix.y);
        ftmp.z = modff(t.z, &tix.z);

        cSim.pPmeParticleFraction[i] = ftmp;

        /* avoid costly % operations if possible that is if dc_ngrid.* is pow. of 2 */
        int4 itmp = make_int4(fast_mod(__float2int_rd(tix.x), cSim.pmeGridSize.x),
                              fast_mod(__float2int_rd(tix.y), cSim.pmeGridSize.y),
                              fast_mod(__float2int_rd(tix.z), cSim.pmeGridSize.z), 0);

        cSim.pPmeParticleIndex[i] = itmp;

        __syncthreads();
    }
}

__global__ void kUpdateBsplines()
{
    unsigned int    tnb = blockDim.x * gridDim.x;
    unsigned int    tid = blockIdx.x * blockDim.x + threadIdx.x;
    extern __shared__ float4 bsplines_cache[]; /*size = 2 * block_size * pme_order*/

    const float4 div_o   = make_float4(1.0f/(PME_ORDER - 1));

    for (int i = tid; i < cSim.atoms; i += tnb)
    {

        float4* data    = &bsplines_cache[threadIdx.x*PME_ORDER];
        float4* ddata   = &bsplines_cache[threadIdx.x*PME_ORDER + blockDim.x*PME_ORDER];

        for (int j = 0; j < PME_ORDER; j++)
        {
	    data[j] = make_float4(0.0f);
            ddata[j] = make_float4(0.0f);
        }

        /* load data */
        float4 dr = cSim.pPmeParticleFraction[i];

        __syncthreads();

        data[PME_ORDER - 1] = make_float4(0.0f);
        data[1]            = dr;
        data[0]            = make_float4(1.0f) - dr;

        for (int j = 3; j < PME_ORDER; j++)
        {
            float div = 1.0f / ((float)j - 1.0f);
            data[j - 1] = div * dr * data[j - 2];

            for (int k = 1; k < (j - 1); k++)
            {
                data[j - k - 1] =
                   div * (
                           (dr + float(k))          * data[j - k - 2] +
                           (-dr + ((float)(j - k))) * data[j - k - 1]);
            }
            data[0] = div * (- dr + 1) * data[0];
        }

        ddata[0] = -data[0];

        for (int j = 1; j < PME_ORDER; j++)
        {
            ddata[j] = data[j - 1] - data[j];
        }

        data[PME_ORDER - 1] = div_o * dr * data[PME_ORDER - 2];

        for (int j = 1; j < (PME_ORDER - 1); j++)
        {
            data[PME_ORDER - j - 1] =
                div_o * (
                    (dr + (float)j)                 * data[PME_ORDER - j - 2] +
                    (-dr + ((float)(PME_ORDER - j))) * data[PME_ORDER - j - 1]
                );
        }
        data[0] = div_o * (-dr + 1.0f) * data[0];

        __syncthreads();

        /* write back */
        for (int j = 0; j < PME_ORDER; j++)
        {
            cSim.pPmeBsplineTheta[i * PME_ORDER + j] =  data[j];
            cSim.pPmeBsplineDtheta[i * PME_ORDER + j] = ddata[j];
        }
        __syncthreads();
    }
}

void kCalculatePME(gpuContext gpu)
{
//    printf("kCalculatePME\n");
    kUpdateGridIndexAndFraction<<<gpu->sim.blocks, gpu->sim.update_threads_per_block>>>();
    LAUNCHERROR("kUpdateGridIndexAndFraction");
    kUpdateBsplines<<<gpu->sim.blocks, gpu->sim.update_threads_per_block, 2*gpu->sim.update_threads_per_block*PME_ORDER*sizeof(float4)>>>();
    LAUNCHERROR("kUpdateBsplines");
}
