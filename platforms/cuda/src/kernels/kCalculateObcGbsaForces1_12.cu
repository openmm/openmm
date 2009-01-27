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

struct Atom {
    float x;
    float y;
    float z;
    float q;
    float br;
    float fx;
    float fy;
    float fz;
    float fb;
};

__shared__ Atom sA[GT2XX_NONBOND_THREADS_PER_BLOCK];
__shared__ unsigned int sWorkUnit[GT2XX_NONBOND_WORKUNITS_PER_SM];
__shared__ unsigned int sNext[GRID];

static __constant__ cudaGmxSimulation cSim;

void SetCalculateObcGbsaForces1_12Sim(gpuContext gpu)
{
    cudaError_t status;
    status = cudaMemcpyToSymbol(cSim, &gpu->sim, sizeof(cudaGmxSimulation));     
    RTERROR(status, "cudaMemcpyToSymbol: SetSim copy to cSim failed");
}

void GetCalculateObcGbsaForces1_12Sim(gpuContext gpu)
{
    cudaError_t status;
    status = cudaMemcpyFromSymbol(&gpu->sim, cSim, sizeof(cudaGmxSimulation));     
    RTERROR(status, "cudaMemcpyFromSymbol: SetSim copy from cSim failed");
}

__global__ void kCalculateObcGbsaForces1_12_kernel()
{
    // Read queue of work blocks once so the remainder of
    // kernel can run asynchronously    
    int pos = cSim.nbWorkUnitsPerBlock * blockIdx.x + min(blockIdx.x, cSim.nbWorkUnitsPerBlockRemainder);
    int end = cSim.nbWorkUnitsPerBlock * (blockIdx.x + 1) + min((blockIdx.x + 1), cSim.nbWorkUnitsPerBlockRemainder);    
    if (threadIdx.x < end - pos)
    {
        sWorkUnit[threadIdx.x] = cSim.pWorkUnit[pos + threadIdx.x];
    }
    if (threadIdx.x < GRID)
    {
        sNext[threadIdx.x] = (threadIdx.x + 1) & (GRID - 1);
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
        float4      apos;   // Local atom x, y, z, q
        float4      af;     // Local atom fx, fy, fz, fb
        unsigned int tgx = threadIdx.x & (GRID - 1);
        unsigned int tbx = threadIdx.x - tgx;
        int tj = tgx; 
        Atom* psA = &sA[tbx];

        if (x == y) // Handle diagonals uniquely at 50% efficiency
        { 
            // Read fixed atom data into registers and GRF
            unsigned int i          = x + tgx;
            apos                    = cSim.pPosq[i];
            float br                = cSim.pBornRadii[i];
            sA[threadIdx.x].x       = apos.x;
            sA[threadIdx.x].y       = apos.y;
            sA[threadIdx.x].z       = apos.z;
            sA[threadIdx.x].q       = apos.w;
            sA[threadIdx.x].br      = br;
            af.x                    = 0.0f;
            af.y                    = 0.0f;
            af.z                    = 0.0f;
            af.w                    = 0.0f;
            apos.w                 *= cSim.preFactor;
            
            for (unsigned int j = 0; j < GRID; j++)
            {
                float dx                = psA[j].x - apos.x; 
                float dy                = psA[j].y - apos.y; 
                float dz                = psA[j].z - apos.z; 
                float r2                = dx * dx + dy * dy + dz * dz; 
                float alpha2_ij         = br * psA[j].br; 
                float D_ij              = r2 / (4.0f * alpha2_ij); 
                float expTerm           = exp(-D_ij); 
                float denominator2      = r2 + alpha2_ij * expTerm; 
                float denominator       = sqrt(denominator2); 
                float Gpol              = (apos.w * psA[j].q) / (denominator * denominator2); 
                float dGpol_dr          = Gpol * (1.0f - 0.25f * expTerm); 
                float dGpol_dalpha2_ij  = -0.5f * Gpol * expTerm * (1.0f + D_ij); 
                dx                     *= dGpol_dr; 
                dy                     *= dGpol_dr; 
                dz                     *= dGpol_dr; 
                af.x                   -= dx; 
                af.y                   -= dy; 
                af.z                   -= dz; 
                af.w                   += dGpol_dalpha2_ij * psA[j].br;      
            }
            
            // Write results
            int offset                          = x + tgx + (x >> GRIDBITS) * cSim.stride;
            cSim.pForce4a[offset]               = af;
            cSim.pBornForce[offset]             = af.w;
        }         
        else        // 100% utilization
        {
            // Read fixed atom data into registers and GRF
            int j                   = y + tgx;
            unsigned int i          = x + tgx;
            float4 temp             = cSim.pPosq[j];
            float temp1             = cSim.pBornRadii[j];
            apos                    = cSim.pPosq[i];
            float br                = cSim.pBornRadii[i];
            sA[threadIdx.x].x       = temp.x;
            sA[threadIdx.x].y       = temp.y;
            sA[threadIdx.x].z       = temp.z;
            sA[threadIdx.x].q       = temp.w;
            sA[threadIdx.x].br      = temp1;
            sA[threadIdx.x].fx      = af.x = 0.0f;
            sA[threadIdx.x].fy      = af.y = 0.0f;
            sA[threadIdx.x].fz      = af.z = 0.0f;
            sA[threadIdx.x].fb      = af.w = 0.0f;
            apos.w                 *= cSim.preFactor;

            for (j = 0; j < GRID; j++)
            {   
                float dx                = psA[tj].x - apos.x; 
                float dy                = psA[tj].y - apos.y; 
                float dz                = psA[tj].z - apos.z; 
                float r2                = dx * dx + dy * dy + dz * dz; 
                float alpha2_ij         = br * psA[tj].br; 
                float D_ij              = r2 / (4.0f * alpha2_ij); 
                float expTerm           = exp(-D_ij); 
                float denominator2      = r2 + alpha2_ij * expTerm; 
                float denominator       = sqrt(denominator2); 
                float Gpol              = (apos.w * psA[tj].q) / (denominator * denominator2); 
                float dGpol_dr          = Gpol * (1.0f - 0.25f * expTerm); 
                float dGpol_dalpha2_ij  = -0.5f * Gpol * expTerm * (1.0f + D_ij); 
                dx                     *= dGpol_dr; 
                dy                     *= dGpol_dr; 
                dz                     *= dGpol_dr; 
                af.x                   -= dx; 
                af.y                   -= dy; 
                af.z                   -= dz; 
                psA[tj].fx             += dx; 
                psA[tj].fy             += dy; 
                psA[tj].fz             += dz; 
                af.w                   += dGpol_dalpha2_ij * psA[tj].br; 
                psA[tj].fb             += dGpol_dalpha2_ij * br;        
                tj                      = sNext[tj]; 
            }
           
            // Write results
            int offset                          = x + tgx + (y >> GRIDBITS) * cSim.stride;
            cSim.pForce4a[offset]               = af;
            cSim.pBornForce[offset]             = af.w;
            offset                              = y + tgx + (x >> GRIDBITS) * cSim.stride;
            af.x                                = sA[threadIdx.x].fx;
            af.y                                = sA[threadIdx.x].fy;
            af.z                                = sA[threadIdx.x].fz;
            af.w                                = sA[threadIdx.x].fb;
            offset                              = y + tgx + (x >> GRIDBITS) * cSim.stride;
            cSim.pForce4a[offset]               = af;
            cSim.pBornForce[offset]             = af.w;
        }
        pos -= cSim.nonbond_workBlock;     
    }
}

void kCalculateObcGbsaForces1_12(gpuContext gpu)
{
  //  printf("kCalculateObcGbsaForces1_12\n");
    kCalculateObcGbsaForces1_12_kernel<<<gpu->sim.nonbond_blocks, gpu->sim.nonbond_threads_per_block>>>();
    LAUNCHERROR("kCalculateObcGbsaForce1_12");
}
