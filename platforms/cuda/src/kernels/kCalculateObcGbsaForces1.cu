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

__shared__ Atom sA[G8X_NONBOND_THREADS_PER_BLOCK];
__shared__ unsigned int sWorkUnit[G8X_NONBOND_WORKUNITS_PER_SM];
__shared__ unsigned int sNext[GRID];

static __constant__ cudaGmxSimulation cSim;

void SetCalculateObcGbsaForces1Sim(gpuContext gpu)
{
    cudaError_t status;
    status = cudaMemcpyToSymbol(cSim, &gpu->sim, sizeof(cudaGmxSimulation));     
    RTERROR(status, "cudaMemcpyToSymbol: SetSim copy to cSim failed");
}

void GetCalculateObcGbsaForces1Sim(gpuContext gpu)
{
    cudaError_t status;
    status = cudaMemcpyFromSymbol(&gpu->sim, cSim, sizeof(cudaGmxSimulation));     
    RTERROR(status, "cudaMemcpyFromSymbol: SetSim copy from cSim failed");
}

__global__ void kReduceObcGbsaBornForces_kernel()
{



    unsigned int pos = (blockIdx.x * blockDim.x + threadIdx.x);
    while (pos < cSim.atoms)
    {
        float bornRadius = cSim.pBornRadii[pos];
        float obcChain   = cSim.pObcChain[pos];
        float2 obcData   = cSim.pObcData[pos];
        float totalForce = 0.0f;
        float* pFt = cSim.pBornForce + pos;
        
        int i = cSim.nonbondOutputBuffers;
        while (i >= 4)
        {
            float f1    = *pFt;
            pFt        += cSim.stride;
            float f2    = *pFt;
            pFt        += cSim.stride;
            float f3    = *pFt;
            pFt        += cSim.stride;
            float f4    = *pFt;
            pFt        += cSim.stride;
            totalForce += f1 + f2 + f3 + f4;
            i -= 4;
        }
        if (i >= 2)
        {
            float f1    = *pFt;
            pFt        += cSim.stride;
            float f2    = *pFt;
            pFt        += cSim.stride;
            totalForce += f1 + f2;
            i -= 2;
        }
        if (i > 0)
        {
            totalForce += *pFt;
        }
        
// __syncthreads();       
        //printf("%4d: %9.4f %9.4f %9.4f\n", pos, totalForce, bornRadius, obcChain);
//totalForce = 0.0f;        

//        if (bornRadius > 0.0f)
//        {
            float r            = (obcData.x + cSim.dielectricOffset + cSim.probeRadius);
            float ratio6       = pow((obcData.x + cSim.dielectricOffset) / bornRadius, 6.0f);
            //float saTerm       = cSim.surfaceAreaFactor * r * r * ratio6;
            float saTerm       = cSim.surfaceAreaFactor * r * r * ratio6;
            totalForce        += saTerm / bornRadius; // 1.102 == Temp mysterious fudge factor, FIX FIX FIX
//        }

        totalForce *= bornRadius * bornRadius * obcChain;
        
        pFt = cSim.pBornForce + pos;
        *pFt = totalForce;
        pos += gridDim.x * blockDim.x;
    }   
}

__global__ void kReduceObcGbsaBornForces1_kernel()
{


    unsigned int pos = (blockIdx.x * blockDim.x + threadIdx.x);
    while (pos < cSim.atoms)
    {
        float bornRadius = cSim.pBornRadii[pos];
        float obcChain   = cSim.pObcChain[pos];
        //float2 obcData   = cSim.pObcData[pos];
        float totalForce = 0.0f;
        float* pFt = cSim.pBornForce + pos;
        
        int i = cSim.nonbondOutputBuffers;
        while (i >= 4)
        {
            float f1    = *pFt;
            pFt        += cSim.stride;
            float f2    = *pFt;
            pFt        += cSim.stride;
            float f3    = *pFt;
            pFt        += cSim.stride;
            float f4    = *pFt;
            pFt        += cSim.stride;
            totalForce += f1 + f2 + f3 + f4;
            i -= 4;
        }
        if (i >= 2)
        {
            float f1    = *pFt;
            pFt        += cSim.stride;
            float f2    = *pFt;
            pFt        += cSim.stride;
            totalForce += f1 + f2;
            i -= 2;
        }
        if (i > 0)
        {
            totalForce += *pFt;
        }
        
// __syncthreads();       
        //printf("%4d: %9.4f %9.4f %9.4f\n", pos, totalForce, bornRadius, obcChain);
//totalForce = 0.0f;        

/*
//        if (bornRadius > 0.0f)
//        {
            float r            = (obcData.x + cSim.dielectricOffset + cSim.probeRadius);
            float ratio6       = pow((obcData.x + cSim.dielectricOffset) / bornRadius, 6.0f);
            float saTerm       = cSim.surfaceAreaFactor * r * r * ratio6;
            totalForce        += saTerm / bornRadius; // 1.102 == Temp mysterious fudge factor, FIX FIX FIX
//        }

		  */

        totalForce *= bornRadius * bornRadius * obcChain;
        
        cSim.pBornForce[pos] = totalForce;
        pos += gridDim.x * blockDim.x;
    }   
}


__global__ void kAceGbsa_kernel()
{

    unsigned int pos = (blockIdx.x * blockDim.x + threadIdx.x);
    while (pos < cSim.atoms)
    {
        float bornRadius = cSim.pBornRadii[pos];
        float obcChain   = cSim.pObcChain[pos];
        float2 obcData   = cSim.pObcData[pos];
        float totalForce = cSim.pBornForce[pos];
        //float totalForce = 0.0f;
        
        float r            = (obcData.x + cSim.dielectricOffset + cSim.probeRadius);
		  
        float ratio6       = pow((obcData.x + cSim.dielectricOffset) / bornRadius, 6.0f);
		  /*
        float ratio6       = (obcData.x + cSim.dielectricOffset) / bornRadius;
		        ratio6       = ratio6*ratio6;
		        ratio6       = ratio6*ratio6*ratio6;
	*/	
	
        //float saTerm       = 41.84f*cSim.surfaceAreaFactor * r * r * ratio6;
        float saTerm       = cSim.surfaceAreaFactor * r * r * ratio6;
        totalForce        += saTerm / bornRadius; // 1.102 == Temp mysterious fudge factor, FIX FIX FIX
        totalForce        *= bornRadius * bornRadius * obcChain;

        cSim.pBornForce[pos] = totalForce;
        pos += gridDim.x * blockDim.x;
    }   
}

void kReduceObcGbsaBornForces(gpuContext gpu)
{
    //printf("kReduceObcGbsaBornForces QQ\n");
    kReduceObcGbsaBornForces_kernel<<<gpu->sim.blocks, gpu->sim.bf_reduce_threads_per_block>>>();
    //kReduceObcGbsaBornForces1_kernel<<<gpu->sim.blocks, gpu->sim.bf_reduce_threads_per_block>>>();
    //kAceGbsa_kernel<<<gpu->sim.blocks, gpu->sim.bf_reduce_threads_per_block>>>();
    //printf("kReduceObcGbsaBornForces calling gpuDumpObcLoop1 QQ\n");
	 //gpuDumpObcLoop1(gpu);
}


__global__ void kCalculateObcGbsaForces1_kernel()
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

__global__ extern void kCalculateObcGbsaForces1_12_kernel();

void kCalculateObcGbsaForces1(gpuContext gpu)
{
    //printf("kCalculateObcGbsaForces1 version=%d sm_12=%d QQ\n", gpu->sm_version, SM_12);
    if (gpu->sm_version < SM_12)
        kCalculateObcGbsaForces1_kernel<<<gpu->sim.nonbond_blocks, gpu->sim.nonbond_threads_per_block>>>();
    else
        kCalculateObcGbsaForces1_12_kernel<<<gpu->sim.nonbond_blocks, gpu->sim.nonbond_threads_per_block>>>();
    LAUNCHERROR("kCalculateObcGbsaForce1");
}
