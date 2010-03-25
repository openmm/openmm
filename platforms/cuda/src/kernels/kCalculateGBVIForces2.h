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

#include "kCalculateGBVIAux.h"

/**
 * This file contains the kernel for evalauating the second stage of GBSA.  It is included
 * several times in kCalculateGBVIForces2.cu with different #defines to generate
 * different versions of the kernels.
 */

__global__ void 
#if (__CUDA_ARCH__ >= 200)
__launch_bounds__(GF1XX_BORNFORCE2_THREADS_PER_BLOCK, 1)
#elif (__CUDA_ARCH__ >= 130)
__launch_bounds__(GT2XX_BORNFORCE2_THREADS_PER_BLOCK, 1)
#else
__launch_bounds__(G8X_BORNFORCE2_THREADS_PER_BLOCK, 1)
#endif
METHOD_NAME(kCalculateGBVI, Forces2_kernel)(unsigned int* workUnit, unsigned int numWorkUnits)
{
    extern __shared__ Atom sA[];
    unsigned int totalWarps  = cSim.bornForce2_blocks*cSim.bornForce2_threads_per_block/GRID;
    unsigned int warp        = (blockIdx.x*blockDim.x+threadIdx.x)/GRID;
    unsigned int pos         = warp*numWorkUnits/totalWarps;
    unsigned int end         = (warp+1)*numWorkUnits/totalWarps;
#ifdef USE_CUTOFF
    float3* tempBuffer       = (float3*) &sA[cSim.bornForce2_threads_per_block];
#endif

    unsigned int lasty = -0xFFFFFFFF;
    while (pos < end)
    {

        // Extract cell coordinates from appropriate work unit
        unsigned int x                  = workUnit[pos];
        unsigned int y                  = ((x >> 2) & 0x7fff) << GRIDBITS;
        x                               = (x >> 17) << GRIDBITS;
        unsigned int tgx                = threadIdx.x & (GRID - 1);
        unsigned int i                  = x + tgx;
        float4 apos                     = cSim.pPosq[i];
        float4 ar                       = cSim.pGBVIData[i];
        float fb                        = cSim.pBornForce[i];
        unsigned int tbx                = threadIdx.x - tgx;
        unsigned int tj                 = tgx;
        Atom* psA                       = &sA[tbx];
        float3 af;
        sA[threadIdx.x].fx = af.x   = 0.0f;
        sA[threadIdx.x].fy = af.y   = 0.0f;
        sA[threadIdx.x].fz = af.z   = 0.0f;
        if (x == y) // Handle diagonals uniquely at 50% efficiency
        {
            // Read fixed atom data into registers and GRF

            sA[threadIdx.x].x           = apos.x;
            sA[threadIdx.x].y           = apos.y;
            sA[threadIdx.x].z           = apos.z;
            sA[threadIdx.x].r           = ar.x;
            sA[threadIdx.x].sr          = ar.y;
            sA[threadIdx.x].fb          = fb;

            for (unsigned int j = (tgx+1)&(GRID-1); j != tgx; j = (j+1)&(GRID-1))
            {
                float dx                = psA[j].x - apos.x;
                float dy                = psA[j].y - apos.y;
                float dz                = psA[j].z - apos.z;
#ifdef USE_PERIODIC
                dx -= floor(dx/cSim.periodicBoxSizeX+0.5f)*cSim.periodicBoxSizeX;
                dy -= floor(dy/cSim.periodicBoxSizeY+0.5f)*cSim.periodicBoxSizeY;
                dz -= floor(dz/cSim.periodicBoxSizeZ+0.5f)*cSim.periodicBoxSizeZ;
#endif
                float r2                = dx * dx + dy * dy + dz * dz;
                float r                 = sqrt(r2);

                // Atom I Born forces and sum
                float dE                = getGBVI_dE2( r, ar.x, psA[j].sr, fb );
               
#if defined USE_PERIODIC
                if (i >= cSim.atoms || x+j >= cSim.atoms || r2 > cSim.nonbondedCutoffSqr)
                {
                    dE              = 0.0f;
                }
#endif
#if defined USE_CUTOFF
                if (r2 > cSim.nonbondedCutoffSqr)
                {
                    dE              = 0.0f;
                }
#endif
                float d             = dx * dE;
                af.x               -= d;
                psA[j].fx          += d;
                d                   = dy * dE;
                af.y               -= d;
                psA[j].fy          += d;
                d                   = dz * dE;
                af.z               -= d;
                psA[j].fz          += d;
            }

            // Write results
            float4 of;
#ifdef USE_OUTPUT_BUFFER_PER_WARP
            unsigned int offset         = x + tgx + warp*cSim.stride;
#else
            unsigned int offset         = x + tgx + (x >> GRIDBITS) * cSim.stride;
#endif
            of                          = cSim.pForce4[offset];
            of.x                       += af.x + sA[threadIdx.x].fx;
            of.y                       += af.y + sA[threadIdx.x].fy;
            of.z                       += af.z + sA[threadIdx.x].fz;
            cSim.pForce4[offset]       = of;
        }
        else
        {
            // Read fixed atom data into registers and GRF
            if (lasty != y)
            {
                unsigned int j              = y + tgx;
                float4 temp                 = cSim.pPosq[j];
                float4 temp1                = cSim.pGBVIData[j];
                float fb                    = cSim.pBornForce[j];
                sA[threadIdx.x].fb          = fb;
                sA[threadIdx.x].x           = temp.x;
                sA[threadIdx.x].y           = temp.y;
                sA[threadIdx.x].z           = temp.z;
                sA[threadIdx.x].r           = temp1.x;
                sA[threadIdx.x].sr          = temp1.y;
            }
#ifdef USE_CUTOFF
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
                    float dx                = psA[tj].x - apos.x;
                    float dy                = psA[tj].y - apos.y;
                    float dz                = psA[tj].z - apos.z;
#ifdef USE_PERIODIC
                    dx -= floor(dx/cSim.periodicBoxSizeX+0.5f)*cSim.periodicBoxSizeX;
                    dy -= floor(dy/cSim.periodicBoxSizeY+0.5f)*cSim.periodicBoxSizeY;
                    dz -= floor(dz/cSim.periodicBoxSizeZ+0.5f)*cSim.periodicBoxSizeZ;
#endif
                    float r2                = dx * dx + dy * dy + dz * dz;
                    float r                 = sqrt(r2);

                    float dE                = getGBVI_dE2( r, ar.x, psA[tj].sr, fb );

#if defined USE_PERIODIC
                    if (i >= cSim.atoms || y+tj >= cSim.atoms || r2 > cSim.nonbondedCutoffSqr)
                    {
                        dE                  = 0.0f;
                    }
#endif
#if defined USE_CUTOFF
                    if (r2 > cSim.nonbondedCutoffSqr)
                    {
                        dE                  = 0.0f;
                    }
#endif

                    float d                 = dx * dE;
                    af.x                   -= d;
                    psA[tj].fx             += d;
                    d                       = dy * dE;
                    af.y                   -= d;
                    psA[tj].fy             += d;
                    d                       = dz * dE;
                    af.z                   -= d;
                    psA[tj].fz             += d;

                    // Atom J Born sum term
                    dE                      = getGBVI_dE2( r, psA[tj].r, ar.y, psA[tj].fb );

#ifdef USE_PERIODIC
                    if (i >= cSim.atoms || y+tj >= cSim.atoms || r2 > cSim.nonbondedCutoffSqr)
                    {
                        dE                  = 0.0f;
                    }
#endif
#if defined USE_CUTOFF
                    if (r2 > cSim.nonbondedCutoffSqr)
                    {
                        dE                  = 0.0f;
                    }
#endif
                    dx                     *= dE;
                    dy                     *= dE;
                    dz                     *= dE;
                    psA[tj].fx             += dx;
                    psA[tj].fy             += dy;
                    psA[tj].fz             += dz;
                    af.x                   -= dx;
                    af.y                   -= dy;
                    af.z                   -= dz;
                    tj                      = (tj + 1) & (GRID - 1);
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
                        float dx                = psA[j].x - apos.x;
                        float dy                = psA[j].y - apos.y;
                        float dz                = psA[j].z - apos.z;
#ifdef USE_PERIODIC
                        dx -= floor(dx/cSim.periodicBoxSizeX+0.5f)*cSim.periodicBoxSizeX;
                        dy -= floor(dy/cSim.periodicBoxSizeY+0.5f)*cSim.periodicBoxSizeY;
                        dz -= floor(dz/cSim.periodicBoxSizeZ+0.5f)*cSim.periodicBoxSizeZ;
#endif
                        float r2                = dx * dx + dy * dy + dz * dz;
                        float r                 = sqrt(r2);

                        // Interleaved Atom I and J Born Forces and sum components
                        float dE                = getGBVI_dE2( r, ar.x, psA[j].sr, fb );

#if defined USE_PERIODIC
                        if (i >= cSim.atoms || y+j >= cSim.atoms || r2 > cSim.nonbondedCutoffSqr)
                        {
                            dE                  = 0.0f;
                        }
#endif
#if defined USE_CUTOFF
                        if (r2 > cSim.nonbondedCutoffSqr)
                        {
                            dE                  = 0.0f;
                        }
#endif

                        float d                 = dx * dE;
                        af.x                   -= d;
                        tempBuffer[threadIdx.x].x = d;
                        d                       = dy * dE;
                        af.y                   -= d;
                        tempBuffer[threadIdx.x].y = d;
                        d                       = dz * dE;
                        af.z                   -= d;
                        tempBuffer[threadIdx.x].z = d;

                        // Atom J Born sum term
                        dE                      = getGBVI_dE2( r, psA[j].r, ar.y, psA[j].fb );

#ifdef USE_PERIODIC
                        if (i >= cSim.atoms || y+j >= cSim.atoms || r2 > cSim.nonbondedCutoffSqr)
                        {
                            dE                  = 0.0f;
                        }
#endif
#if defined USE_CUTOFF
                        if (r2 > cSim.nonbondedCutoffSqr)
                        {
                            dE                  = 0.0f;
                        }
#endif
                        dx                     *= dE;
                        dy                     *= dE;
                        dz                     *= dE;
                        tempBuffer[threadIdx.x].x += dx;
                        tempBuffer[threadIdx.x].y += dy;
                        tempBuffer[threadIdx.x].z += dz;
                        af.x                   -= dx;
                        af.y                   -= dy;
                        af.z                   -= dz;

                        // Sum the forces on atom j.

                        if (tgx % 2 == 0)
                        {
                            tempBuffer[threadIdx.x].x += tempBuffer[threadIdx.x+1].x;
                            tempBuffer[threadIdx.x].y += tempBuffer[threadIdx.x+1].y;
                            tempBuffer[threadIdx.x].z += tempBuffer[threadIdx.x+1].z;
                        }
                        if (tgx % 4 == 0)
                        {
                            tempBuffer[threadIdx.x].x += tempBuffer[threadIdx.x+2].x;
                            tempBuffer[threadIdx.x].y += tempBuffer[threadIdx.x+2].y;
                            tempBuffer[threadIdx.x].z += tempBuffer[threadIdx.x+2].z;
                        }
                        if (tgx % 8 == 0)
                        {
                            tempBuffer[threadIdx.x].x += tempBuffer[threadIdx.x+4].x;
                            tempBuffer[threadIdx.x].y += tempBuffer[threadIdx.x+4].y;
                            tempBuffer[threadIdx.x].z += tempBuffer[threadIdx.x+4].z;
                        }
                        if (tgx % 16 == 0)
                        {
                            tempBuffer[threadIdx.x].x += tempBuffer[threadIdx.x+8].x;
                            tempBuffer[threadIdx.x].y += tempBuffer[threadIdx.x+8].y;
                            tempBuffer[threadIdx.x].z += tempBuffer[threadIdx.x+8].z;
                        }
                        if (tgx == 0)
                        {
                            psA[j].fx += tempBuffer[threadIdx.x].x + tempBuffer[threadIdx.x+16].x;
                            psA[j].fy += tempBuffer[threadIdx.x].y + tempBuffer[threadIdx.x+16].y;
                            psA[j].fz += tempBuffer[threadIdx.x].z + tempBuffer[threadIdx.x+16].z;
                        }
                    }
                }
            }
#endif

            // Write results
            float4 of;
#ifdef USE_OUTPUT_BUFFER_PER_WARP
            unsigned int offset         = x + tgx + warp*cSim.stride;
#else
            unsigned int offset         = x + tgx + (y >> GRIDBITS) * cSim.stride;
#endif
            of                          = cSim.pForce4[offset];
            of.x                       += af.x;
            of.y                       += af.y;
            of.z                       += af.z;
            cSim.pForce4[offset]       = of;
#ifdef USE_OUTPUT_BUFFER_PER_WARP
            offset                      = y + tgx + warp*cSim.stride;
#else
            offset                      = y + tgx + (x >> GRIDBITS) * cSim.stride;
#endif
            of                          = cSim.pForce4[offset];
            of.x                       += sA[threadIdx.x].fx;
            of.y                       += sA[threadIdx.x].fy;
            of.z                       += sA[threadIdx.x].fz;
            cSim.pForce4[offset]       = of;
        }
        lasty = y;
        pos++;
    }
}
