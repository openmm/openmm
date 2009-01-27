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
#include <fstream>
using namespace std;

#include "gputypes.h"

#define UNROLLXX 0
#define UNROLLXY 0

struct Atom {
    float x;
    float y;
    float z;
    float r;
    float sr;
    float sum;
    float junk;
};

__shared__ Atom sA[GT2XX_NONBOND_THREADS_PER_BLOCK];
__shared__ unsigned int sWorkUnit[GT2XX_NONBOND_WORKUNITS_PER_SM];
__shared__ unsigned int sNext[GRID];

static __constant__ cudaGmxSimulation cSim;

void SetCalculateObcGbsaBornSumSim(gpuContext gpu)
{
    cudaError_t status;
    status = cudaMemcpyToSymbol(cSim, &gpu->sim, sizeof(cudaGmxSimulation));     
    RTERROR(status, "cudaMemcpyToSymbol: SetSim copy to cSim failed");
}

void GetCalculateObcGbsaBornSumSim(gpuContext gpu)
{
    cudaError_t status;
    status = cudaMemcpyFromSymbol(&gpu->sim, cSim, sizeof(cudaGmxSimulation));     
    RTERROR(status, "cudaMemcpyFromSymbol: SetSim copy from cSim failed");
}

__global__ void kReduceObcGbsaBornSum_kernel()
{
    unsigned int pos = (blockIdx.x * blockDim.x + threadIdx.x);
    
    while (pos < cSim.atoms)
    {
        float sum = 0.0f;
        float* pSt = cSim.pBornSum + pos;
        float2 atom = cSim.pObcData[pos];
        
        // Get summed Born data
        for (int i = 0; i < cSim.nonbondOutputBuffers; i++)
        {
            sum += *pSt;
       //     printf("%4d %4d A: %9.4f\n", pos, i, *pSt);
            pSt += cSim.stride;
        }
        
        
        // Now calculate Born radius and OBC term.
        sum                    *= 0.5f * atom.x;
        float sum2              = sum * sum;
        float sum3              = sum * sum2;
        float tanhSum           = tanh(cSim.alphaOBC * sum - cSim.betaOBC * sum2 + cSim.gammaOBC * sum3);
        float nonOffsetRadii    = atom.x + cSim.dielectricOffset;
        float bornRadius        = 1.0f / (1.0f / atom.x - tanhSum / nonOffsetRadii); 
        float obcChain          = atom.x * (cSim.alphaOBC - 2.0f * cSim.betaOBC * sum + 3.0f * cSim.gammaOBC * sum2);
        obcChain                = (1.0f - tanhSum * tanhSum) * obcChain / nonOffsetRadii;        
        cSim.pBornRadii[pos] = bornRadius;
        cSim.pObcChain[pos]  = obcChain;
        pos += gridDim.x * blockDim.x;
    }   
}

void kReduceObcGbsaBornSum(gpuContext gpu)
{
//    printf("kReduceObcGbsaBornSum\n");
    kReduceObcGbsaBornSum_kernel<<<gpu->sim.blocks, 384>>>();
    gpu->bRecalculateBornRadii = false;

if( 0 ){
   static int step = 0;
   int numPrint    = -1;
   step++;
   WriteArrayToFile1( gpu, "ObcGbsaBornBRad", step, gpu->psBornRadii, numPrint );
   WriteArrayToFile1( gpu, "ObcGbsaBornSum", step, gpu->psBornSum, numPrint );
   WriteArrayToFile2( gpu, "ObcGbsaObcData", step, gpu->psObcData, numPrint );
   WriteArrayToFile4( gpu, "ObcGbsaBornPos", step, gpu->psPosq4, numPrint );
   //gpuDumpCoordinates( gpu );
   gpuDumpObcInfo( gpu );
}
    LAUNCHERROR("kReduceObcGbsaBornSum");
}


__global__ void kCalculateObcGbsaBornSum_kernel()
{
    // Read queue of work blocks once so the remainder of
    // kernel can run asynchronously    
    int pos = (blockIdx.x * cSim.workUnits) / gridDim.x;
    int end = ((blockIdx.x + 1) * cSim.workUnits) / gridDim.x;
    if (threadIdx.x < end - pos)
    {
        sWorkUnit[threadIdx.x] = cSim.pWorkUnit[pos + threadIdx.x];
    }
    if (threadIdx.x < GRID)
    {
        sNext[threadIdx.x] = (threadIdx.x - 1) & (GRID - 1);
    }
    __syncthreads();

    // Now change pos and end to reflect work queue just read
    // into shared memory
    end = end - pos; 
    pos = end - (threadIdx.x >> GRIDBITS) - 1;
       
    while (pos >= 0)
    {  
        // Extract cell coordinates from appropriate work unit
        unsigned int x = sWorkUnit[pos];
        unsigned int y = ((x >> 2) & 0x7fff) << GRIDBITS;
        x = (x >> 17) << GRIDBITS;
        float       dx; 
        float       dy; 
        float       dz; 
        float       r2; 
        float       r;

        unsigned int tgx = threadIdx.x & (GRID - 1);
        unsigned int tbx = threadIdx.x - tgx;
        int tj = tgx; 
        Atom* psA = &sA[tbx];
     
        if (x == y) // Handle diagonals uniquely at 50% efficiency
        { 
            // Read fixed atom data into registers and GRF       
            unsigned int i = x + tgx;
            float4 apos = cSim.pPosq[i];    // Local atom x, y, z, sum
            float2 ar = cSim.pObcData[i];   // Local atom vr, sr
            sA[threadIdx.x].x           = apos.x;
            sA[threadIdx.x].y           = apos.y;
            sA[threadIdx.x].z           = apos.z;
            sA[threadIdx.x].r           = ar.x;
            sA[threadIdx.x].sr          = ar.y;
            apos.w                      = 0.0f;

            for (unsigned int j = 0; j < GRID; j++)
            {
                dx                      = psA[j].x - apos.x;
                dy                      = psA[j].y - apos.y;
                dz                      = psA[j].z - apos.z;
                r2                      = dx * dx + dy * dy + dz * dz; 
                r                       = sqrt(r2);
                float rInverse          = 1.0f / r; 
                float rScaledRadiusJ    = r + psA[j].sr;
                if ((j != tgx) && (ar.x < rScaledRadiusJ))
                {
                    float l_ij     = 1.0f / max(ar.x, fabs(r - psA[j].sr));
                    float u_ij     = 1.0f / rScaledRadiusJ;
                    float l_ij2    = l_ij * l_ij;
                    float u_ij2    = u_ij * u_ij;
                    float ratio    = log(u_ij / l_ij);
                    apos.w        += l_ij - 
                                     u_ij + 
                                     0.25f * r * (u_ij2 - l_ij2) + 
                                     (0.50f * rInverse * ratio) + 
                                     (0.25f * psA[j].sr * psA[j].sr * rInverse) *
                                     (l_ij2 - u_ij2);
                                                                                                              
                    if (ar.x < (psA[j].r - r))
                    {
                        apos.w += 2.0f * ((1.0f / ar.x) - l_ij);
                    }
                }
            }             

            // Write results
            int offset = x + tgx + (x >> GRIDBITS) * cSim.stride;
            cSim.pBornSum[offset] = apos.w;
        }         
        else        // 100% utilization
        {
            // Read fixed atom data into registers and GRF
            int j                           = y + tgx;
            unsigned int i                  = x + tgx;      
            
            float4 temp                     = cSim.pPosq[j];
            float2 temp1                    = cSim.pObcData[j];
            float4 apos                     = cSim.pPosq[i];        // Local atom x, y, z, sum
            float2 ar                       = cSim.pObcData[i];    // Local atom vr, sr
            sA[threadIdx.x].x               = temp.x;
            sA[threadIdx.x].y               = temp.y;
            sA[threadIdx.x].z               = temp.z;
            sA[threadIdx.x].r               = temp1.x;
            sA[threadIdx.x].sr              = temp1.y;
            sA[threadIdx.x].sum = apos.w    = 0.0f;

            for (unsigned int j = 0; j < GRID; j++)
            {
                dx                      = psA[tj].x - apos.x; 
                dy                      = psA[tj].y - apos.y; 
                dz                      = psA[tj].z - apos.z; 
                r2                      = dx * dx + dy * dy + dz * dz; 
                r                       = sqrt(r2);
                float rInverse          = 1.0f / r; 
                float rScaledRadiusJ    = r + psA[tj].sr;
                if (ar.x < rScaledRadiusJ)
                {
                    float l_ij     = 1.0f / max(ar.x, fabs(r - psA[tj].sr));
                    float u_ij     = 1.0f / rScaledRadiusJ;
                    float l_ij2    = l_ij * l_ij;
                    float u_ij2    = u_ij * u_ij;
                    float ratio    = log(u_ij / l_ij);
                    float term     = l_ij - 
                                     u_ij + 
                                     0.25f * r * (u_ij2 - l_ij2) + 
                                     (0.50f * rInverse * ratio) + 
                                     (0.25f * psA[tj].sr * psA[tj].sr * rInverse) *
                                     (l_ij2 - u_ij2);
                    if (ar.x < (psA[tj].sr - r))
                    {
                        term += 2.0f * ((1.0f / ar.x) - l_ij);
                    }
                    apos.w        += term;
                }
                float rScaledRadiusI    = r + ar.y;
                if (psA[tj].r < rScaledRadiusI)
                {
                    float l_ij     = 1.0f / max(psA[tj].r, fabs(r - ar.y));
                    float u_ij     = 1.0f / rScaledRadiusI;
                    float l_ij2    = l_ij * l_ij;
                    float u_ij2    = u_ij * u_ij;
                    float ratio    = log(u_ij / l_ij);
                    float term     = l_ij - 
                                     u_ij + 
                                     0.25f * r * (u_ij2 - l_ij2) + 
                                     (0.50f * rInverse * ratio) + 
                                     (0.25f * ar.y * ar.y * rInverse) *
                                     (l_ij2 - u_ij2);
 
                    if (psA[tj].r < (ar.y - r))
                    {
                        term += 2.0f * ((1.0f / psA[tj].r) - l_ij);
                    }
                    psA[tj].sum    += term;
                }      
                tj = sNext[tj];
            }    
                
            // Write results
            int offset = x + tgx + (y >> GRIDBITS) * cSim.stride;
            cSim.pBornSum[offset] = apos.w;
            offset = y + tgx + (x >> GRIDBITS) * cSim.stride;
            cSim.pBornSum[offset] = sA[threadIdx.x].sum;
        }       
       
        pos -= cSim.nonbond_workBlock;     
    }
}

void kCalculateObcGbsaBornSum(gpuContext gpu)
{
  //  printf("kCalculateObcgbsaBornSum\n");
    kCalculateObcGbsaBornSum_kernel<<<gpu->sim.nonbond_blocks, gpu->sim.nonbond_threads_per_block>>>();
    LAUNCHERROR("kCalculateBornSum");
}
