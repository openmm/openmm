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

__global__ void kUpdateGridIndexAndFraction_kernel()
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

__global__ void kUpdateBsplines_kernel()
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

__global__ void kGridSpreadCharge_kernel()
{
    extern __shared__ float4 atomPos[];
    int4* atomGridIndex = (int4*) &atomPos[blockDim.x];
    const unsigned int totalWarps = cSim.nonbond_blocks*cSim.nonbond_threads_per_block/GRID;
    const unsigned int warp = (blockIdx.x*blockDim.x+threadIdx.x)/GRID;
    const int3 groupDim = make_int3(4, 4, 2);
    const int3 numGroups = make_int3((cSim.pmeGridSize.x+groupDim.x-1)/groupDim.x, (cSim.pmeGridSize.y+groupDim.y-1)/groupDim.y, (cSim.pmeGridSize.z+groupDim.z-1)/groupDim.z);
    const unsigned int totalGroups = numGroups.x*numGroups.y*numGroups.z;
    unsigned int group = warp*totalGroups/totalWarps;
    const unsigned int end = (warp+1)*totalGroups/totalWarps;
    const unsigned int index = threadIdx.x & (GRID - 1);

    while (group < end)
    {
        // Process a group of grid points of size groupDim.  First figure out the base index for the group,
        // and the index of the specific point this thread will handle.

        int3 gridBase;
        gridBase.x = group/(numGroups.y*numGroups.z);
        int remainder = group-gridBase.x*numGroups.y*numGroups.z;
        gridBase.y = remainder/numGroups.z;
        gridBase.z = remainder-gridBase.y*numGroups.z;
        gridBase.x *= groupDim.x;
        gridBase.y *= groupDim.y;
        gridBase.z *= groupDim.z;
        int3 gridPoint;
        gridPoint.x = index/(groupDim.y*groupDim.z);
        remainder = index-gridPoint.x*groupDim.y*groupDim.z;
        gridPoint.y = remainder/groupDim.z;
        gridPoint.z = remainder-gridPoint.y*groupDim.z;
        gridPoint.x += gridBase.x;
        gridPoint.y += gridBase.y;
        gridPoint.z += gridBase.z;

        // Loop over blocks of atoms.

        float result = 0.0f;
        for (int atomBlock = 0; atomBlock < cSim.paddedNumberOfAtoms>>GRIDBITS; atomBlock++)
        {
            int atomIndex = (atomBlock<<GRIDBITS)+index;
            if (atomIndex < cSim.atoms)
            {
                atomPos[threadIdx.x] = cSim.pPosq[atomIndex];
                atomGridIndex[threadIdx.x] = cSim.pPmeParticleIndex[atomIndex];
            }
            int maxAtoms = min(GRID, cSim.atoms-(atomBlock<<GRIDBITS));
            for (int i = 0; i < maxAtoms; i++)
            {
                int atomIndex = (atomBlock<<GRIDBITS)+i;
                int ix = gridPoint.x-atomGridIndex[threadIdx.x-index+i].x;
                int iy = gridPoint.y-atomGridIndex[threadIdx.x-index+i].y;
                int iz = gridPoint.z-atomGridIndex[threadIdx.x-index+i].z;
                if (ix < 0)
                    ix += cSim.pmeGridSize.x;
                if (iy < 0)
                    iy += cSim.pmeGridSize.y;
                if (iz < 0)
                    iz += cSim.pmeGridSize.z;
                if (ix < PME_ORDER && iy < PME_ORDER && iz < PME_ORDER)
                    result += atomPos[threadIdx.x-index+i].w*cSim.pPmeBsplineTheta[atomIndex*PME_ORDER+ix].x*cSim.pPmeBsplineTheta[atomIndex*PME_ORDER+iy].y*cSim.pPmeBsplineTheta[atomIndex*PME_ORDER+iz].z;
            }
        }
        unsigned int gridIndex = gridPoint.x*cSim.pmeGridSize.y*cSim.pmeGridSize.z+gridPoint.y*cSim.pmeGridSize.z+gridPoint.z;
        if (gridIndex < cSim.pmeGridSize.x*cSim.pmeGridSize.y*cSim.pmeGridSize.z)
            cSim.pPmeGrid[gridIndex] = make_cuComplex(result*sqrt(cSim.epsfac), 0.0f);
        group++;
    }
}

__global__ void kReciprocalConvolution_kernel()
{
    const unsigned int gridSize = cSim.pmeGridSize.x*cSim.pmeGridSize.y*cSim.pmeGridSize.z;
    float expFactor = PI*PI/(cSim.alphaEwald*cSim.alphaEwald);
    float scaleFactor = 1.0/(PI*cSim.periodicBoxSizeX*cSim.periodicBoxSizeY*cSim.periodicBoxSizeZ);
    float energy = 0.0f;
    for (int index = blockIdx.x*blockDim.x+threadIdx.x; index < gridSize; index += blockDim.x*gridDim.x)
    {
        int kx = index/(cSim.pmeGridSize.y*cSim.pmeGridSize.z);
        int remainder = index-kx*cSim.pmeGridSize.y*cSim.pmeGridSize.z;
        int ky = remainder/cSim.pmeGridSize.z;
        int kz = remainder-ky*cSim.pmeGridSize.z;
        if (kx == 0 && ky == 0 && kz == 0)
            continue;
        int mx = (kx < (cSim.pmeGridSize.x+1)/2) ? kx : (kx-cSim.pmeGridSize.x);
        int my = (ky < (cSim.pmeGridSize.y+1)/2) ? ky : (ky-cSim.pmeGridSize.y);
        int mz = (kz < (cSim.pmeGridSize.z+1)/2) ? kz : (kz-cSim.pmeGridSize.z);
        float mhx = mx/cSim.periodicBoxSizeX;
        float mhy = my/cSim.periodicBoxSizeY;
        float mhz = mz/cSim.periodicBoxSizeZ;
        float bx = cSim.pPmeBsplineModuli[0][kx];
        float by = cSim.pPmeBsplineModuli[1][ky];
        float bz = cSim.pPmeBsplineModuli[2][kz];
        cuComplex grid = cSim.pPmeGrid[index];
        float m2 = mhx*mhx+mhy*mhy+mhz*mhz;
        float denom = m2*bx*by*bz;
        float eterm = scaleFactor*exp(-expFactor*m2)/denom;
        cSim.pPmeGrid[index] = make_cuComplex(grid.x*eterm, grid.y*eterm);
        energy += eterm*(grid.x*grid.x + grid.y*grid.y);
    }
    cSim.pEnergy[blockIdx.x*blockDim.x+threadIdx.x] += 0.5f*energy;
}

__global__ void kGridInterpolateForce_kernel()
{
    for (int atom = blockIdx.x*blockDim.x+threadIdx.x; atom < cSim.atoms; atom += blockDim.x*gridDim.x)
    {
        float3 force = make_float3(0.0f, 0.0f, 0.0f);
        float4 posq = cSim.pPosq[atom];
        int4 gridIndex = cSim.pPmeParticleIndex[atom];
        for (int ix = 0; ix < PME_ORDER; ix++)
        {
            int xindex = (gridIndex.x + ix) % cSim.pmeGridSize.x;
            float tx = cSim.pPmeBsplineTheta[atom*PME_ORDER+ix].x;
            float dtx = cSim.pPmeBsplineDtheta[atom*PME_ORDER+ix].x;
            for (int iy = 0; iy < PME_ORDER; iy++)
            {
                int yindex = (gridIndex.y + iy) % cSim.pmeGridSize.y;
                float ty = cSim.pPmeBsplineTheta[atom*PME_ORDER+iy].y;
                float dty = cSim.pPmeBsplineDtheta[atom*PME_ORDER+iy].y;
                for (int iz = 0; iz < PME_ORDER; iz++)
                {
                    int zindex               = (gridIndex.z + iz) % cSim.pmeGridSize.z;
                    float tz = cSim.pPmeBsplineTheta[atom*PME_ORDER+iz].z;
                    float dtz = cSim.pPmeBsplineDtheta[atom*PME_ORDER+iz].z;
                    int index                = xindex*cSim.pmeGridSize.y*cSim.pmeGridSize.z + yindex*cSim.pmeGridSize.z + zindex;
                    float gridvalue            = cSim.pPmeGrid[index].x;
                    force.x                  += dtx*ty*tz*gridvalue;
                    force.y                  += tx*dty*tz*gridvalue;
                    force.z                  += tx*ty*dtz*gridvalue;
                }
            }
        }
        float4 totalForce = cSim.pForce4[atom];
        float q = posq.w*sqrt(cSim.epsfac);
        totalForce.x -= q*force.x*cSim.pmeGridSize.x/cSim.periodicBoxSizeX;
        totalForce.y -= q*force.y*cSim.pmeGridSize.y/cSim.periodicBoxSizeY;
        totalForce.z -= q*force.z*cSim.pmeGridSize.z/cSim.periodicBoxSizeZ;
        cSim.pForce4[atom] = totalForce;
    }
}

void kCalculatePME(gpuContext gpu)
{
//    printf("kCalculatePME\n");
    kUpdateGridIndexAndFraction_kernel<<<gpu->sim.blocks, gpu->sim.update_threads_per_block>>>();
    LAUNCHERROR("kUpdateGridIndexAndFraction");
    kUpdateBsplines_kernel<<<gpu->sim.blocks, gpu->sim.update_threads_per_block, 2*gpu->sim.update_threads_per_block*PME_ORDER*sizeof(float4)>>>();
    LAUNCHERROR("kUpdateBsplines");
    kGridSpreadCharge_kernel<<<gpu->sim.blocks, 64, 64*(sizeof(float4)+sizeof(int4))>>>();
    LAUNCHERROR("kGridSpreadChargeß");
    cufftExecC2C(gpu->fftplan, gpu->psPmeGrid->_pDevData, gpu->psPmeGrid->_pDevData, CUFFT_FORWARD);
    kReciprocalConvolution_kernel<<<gpu->sim.blocks, gpu->sim.nonbond_threads_per_block>>>();
    LAUNCHERROR("kReciprocalConvolution");
    cufftExecC2C(gpu->fftplan, gpu->psPmeGrid->_pDevData, gpu->psPmeGrid->_pDevData, CUFFT_INVERSE);
    kGridInterpolateForce_kernel<<<gpu->sim.blocks, gpu->sim.update_threads_per_block>>>();
    LAUNCHERROR("kGridInterpolateForce");
}
