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

/**
 * This file contains the kernel for calculating Born sums.  It is included
 * several times in kCalculateGBVIBornSum.cu with different #defines to generate
 * different versions of the kernels.
 */

#include "kCalculateGBVIAux.h"

__global__ void METHOD_NAME(kCalculateGBVI, BornSum_kernel)(unsigned int* workUnit)
{
    extern __shared__ Atom sA[];

    unsigned int totalWarps   = cSim.nonbond_blocks*cSim.nonbond_threads_per_block/GRID;
    unsigned int warp         = (blockIdx.x*blockDim.x+threadIdx.x)/GRID;
    unsigned int numWorkUnits = cSim.pInteractionCount[0];
    unsigned int pos          = warp*numWorkUnits/totalWarps;
    unsigned int end          = (warp+1)*numWorkUnits/totalWarps;

//    int end = workUnits / gridDim.x;
//    int pos = end - (threadIdx.x >> GRIDBITS) - 1;
#ifdef USE_CUTOFF
    float* tempBuffer = (float*) &sA[cSim.nonbond_threads_per_block];
#endif

    while ( pos < end )
    {
        // Extract cell coordinates from appropriate work unit
        unsigned int x = workUnit[pos];
        unsigned int y = ((x >> 2) & 0x7fff) << GRIDBITS;
        x              = (x >> 17)           << GRIDBITS;

        float       dx;
        float       dy;
        float       dz;
        float       r2;
        float       r;

        // forces tgx into interval [0,31] 
        // forces tbx 0 
        unsigned int tgx = threadIdx.x & (GRID - 1);
        unsigned int tbx = threadIdx.x - tgx;
        unsigned int tj  = tgx;
        Atom* psA        = &sA[tbx];

        if (x == y) // Handle diagonals uniquely at 50% efficiency
        {
            // Read fixed atom data into registers and GRF
            unsigned int i                          = x + tgx;
            float4 apos                             = cSim.pPosq[i];    // Local atom x, y, z, sum
            float4 ar                               = cSim.pGBVIData[i];  // Local atom vr, sr
            sA[threadIdx.x].x                       = apos.x;
            sA[threadIdx.x].y                       = apos.y;
            sA[threadIdx.x].z                       = apos.z;
            sA[threadIdx.x].r                       = ar.x;
            sA[threadIdx.x].sr                      = ar.y;
            apos.w                                  = 0.0f;

            for (unsigned int j             = 0; j < GRID; j++)
            {
                dx                                  = psA[j].x - apos.x;
                dy                                  = psA[j].y - apos.y;
                dz                                  = psA[j].z - apos.z;
#ifdef USE_PERIODIC
                dx -= floor(dx/cSim.periodicBoxSizeX+0.5f)*cSim.periodicBoxSizeX;
                dy -= floor(dy/cSim.periodicBoxSizeY+0.5f)*cSim.periodicBoxSizeY;
                dz -= floor(dz/cSim.periodicBoxSizeZ+0.5f)*cSim.periodicBoxSizeZ;
#endif
                r2                      = dx * dx + dy * dy + dz * dz;
#if defined USE_PERIODIC
                if (i < cSim.atoms && x+j < cSim.atoms && r2 < cSim.nonbondedCutoffSqr)
#elif defined USE_CUTOFF
                if (r2 < cSim.nonbondedCutoffSqr)
#endif
                {
                    r                       = sqrt(r2);
                    if ((j != tgx) )
                    {
                        apos.w             += getGBVI_Volume( r, ar.x, psA[j].sr );
                    }
                }
            }

            // Write results
#ifdef USE_OUTPUT_BUFFER_PER_WARP
            unsigned int offset = x + tgx + warp*cSim.stride;
            cSim.pBornSum[offset] += apos.w;
#else
            unsigned int offset = x + tgx + (x >> GRIDBITS) * cSim.stride;
            cSim.pBornSum[offset] = apos.w;
#endif
        }
        else        // 100% utilization
        {
            // Read fixed atom data into registers and GRF
            unsigned int j                              = y + tgx;
            unsigned int i                              = x + tgx;

            float4 temp                                 = cSim.pPosq[j];
            float4 temp1                                = cSim.pGBVIData[j];
            float4 apos                                 = cSim.pPosq[i];        // Local atom x, y, z, sum
            float4 ar                                   = cSim.pGBVIData[i];    // Local atom vr, sr
            sA[threadIdx.x].x                           = temp.x;
            sA[threadIdx.x].y                           = temp.y;
            sA[threadIdx.x].z                           = temp.z;
            sA[threadIdx.x].r                           = temp1.x;
            sA[threadIdx.x].sr                          = temp1.y;
            sA[threadIdx.x].sum             = apos.w    = 0.0f;

#ifdef USE_CUTOFF
            //unsigned int flags = cSim.pInteractionFlag[pos + (blockIdx.x*workUnits)/gridDim.x];
            unsigned int flags = cSim.pInteractionFlag[pos];
            if (flags == 0)
            {
                // No interactions in this block.
            }
            else if (flags == 0xFFFFFFFF)
#endif
            {
                // Compute all interactions within this block.

                for (unsigned int j = 0; j < GRID; j++)
                {
                    dx                      = psA[tj].x - apos.x;
                    dy                      = psA[tj].y - apos.y;
                    dz                      = psA[tj].z - apos.z;
#ifdef USE_PERIODIC
                    dx -= floor(dx/cSim.periodicBoxSizeX+0.5f)*cSim.periodicBoxSizeX;
                    dy -= floor(dy/cSim.periodicBoxSizeY+0.5f)*cSim.periodicBoxSizeY;
                    dz -= floor(dz/cSim.periodicBoxSizeZ+0.5f)*cSim.periodicBoxSizeZ;
#endif
                    r2                      = dx * dx + dy * dy + dz * dz;
#ifdef USE_PERIODIC
                    if (i < cSim.atoms && y+tj < cSim.atoms && r2 < cSim.nonbondedCutoffSqr)
#elif defined USE_CUTOFF
                    if (r2 < cSim.nonbondedCutoffSqr)
#endif
                    {
                        r                       = sqrt(r2);

                        // psA[tj].sr = Sj
                        // ar.x       = Ri

                        apos.w                 += getGBVI_Volume( r, ar.x,      psA[tj].sr );
                        psA[tj].sum            += getGBVI_Volume( r, psA[tj].r, ar.y );
                    }
                    tj = (tj - 1) & (GRID - 1);
                }
            }
#ifdef USE_CUTOFF
            else
            {
                // Compute only a subset of the interactions in this block.

                for (unsigned int j = 0; j < GRID; j++)
                {
                    if ((flags&(1<<j)) != 0)
                    {
                        tempBuffer[threadIdx.x] = 0.0f;
                        dx                      = psA[j].x - apos.x;
                        dy                      = psA[j].y - apos.y;
                        dz                      = psA[j].z - apos.z;
#ifdef USE_PERIODIC
                        dx -= floor(dx/cSim.periodicBoxSizeX+0.5f)*cSim.periodicBoxSizeX;
                        dy -= floor(dy/cSim.periodicBoxSizeY+0.5f)*cSim.periodicBoxSizeY;
                        dz -= floor(dz/cSim.periodicBoxSizeZ+0.5f)*cSim.periodicBoxSizeZ;
#endif
                        r2                      = dx * dx + dy * dy + dz * dz;
#ifdef USE_PERIODIC
                        if (i < cSim.atoms && y+j < cSim.atoms && r2 < cSim.nonbondedCutoffSqr)
#elif defined USE_CUTOFF
                        if (r2 < cSim.nonbondedCutoffSqr)
#endif
                        {
                            r                       = sqrt(r2);
                            tempBuffer[threadIdx.x] = getGBVI_Volume( r, psA[tj].r, ar.y );
                        }

                        // Sum the terms.

                        if (tgx % 2 == 0)
                            tempBuffer[threadIdx.x] += tempBuffer[threadIdx.x+1];
                        if (tgx % 4 == 0)
                            tempBuffer[threadIdx.x] += tempBuffer[threadIdx.x+2];
                        if (tgx % 8 == 0)
                            tempBuffer[threadIdx.x] += tempBuffer[threadIdx.x+4];
                        if (tgx % 16 == 0)
                            tempBuffer[threadIdx.x] += tempBuffer[threadIdx.x+8];
                        if (tgx == 0)
                            psA[j].sum += tempBuffer[threadIdx.x] + tempBuffer[threadIdx.x+16];
                    }
                }
            }
#endif

            // Write results
#ifdef USE_OUTPUT_BUFFER_PER_WARP
            unsigned int offset = x + tgx + warp*cSim.stride;
            cSim.pBornSum[offset] += apos.w;
            offset = y + tgx + warp*cSim.stride;
            cSim.pBornSum[offset] += sA[threadIdx.x].sum;
#else
            unsigned int offset = x + tgx + (y >> GRIDBITS) * cSim.stride;
            cSim.pBornSum[offset] = apos.w;
            offset = y + tgx + (x >> GRIDBITS) * cSim.stride;
            cSim.pBornSum[offset] = sA[threadIdx.x].sum;
#endif
        }

        pos++;
    }
}
