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
 * several times in kCalculateObcGbsaBornSum.cu with different #defines to generate
 * different versions of the kernels.
 */

__global__ void METHOD_NAME(kCalculateObcGbsa, BornSum_kernel)(unsigned int* workUnit, int workUnits)
{
    extern __shared__ Atom sA[];
    int end = workUnits / gridDim.x;
    int pos = end - (threadIdx.x >> GRIDBITS) - 1;
#ifdef USE_OUTPUT_BUFFER_PER_WARP
    unsigned int warp = (blockIdx.x*blockDim.x+threadIdx.x)/GRID;
#endif

    while (pos >= 0)
    {
        // Extract cell coordinates from appropriate work unit
        unsigned int x = workUnit[pos + (blockIdx.x*workUnits)/gridDim.x];
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
            }

            // Write results
#ifdef USE_OUTPUT_BUFFER_PER_WARP
            int offset = x + tgx + warp*cSim.stride;
            cSim.pBornSum[offset] += apos.w;
#else
            int offset = x + tgx + (x >> GRIDBITS) * cSim.stride;
            cSim.pBornSum[offset] = apos.w;
#endif
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
                }
                tj = (tj - 1) & (GRID - 1);
            }

            // Write results
#ifdef USE_OUTPUT_BUFFER_PER_WARP
            int offset = x + tgx + warp*cSim.stride;
            cSim.pBornSum[offset] += apos.w;
            offset = y + tgx + warp*cSim.stride;
            cSim.pBornSum[offset] += sA[threadIdx.x].sum;
#else
            int offset = x + tgx + (y >> GRIDBITS) * cSim.stride;
            cSim.pBornSum[offset] = apos.w;
            offset = y + tgx + (x >> GRIDBITS) * cSim.stride;
            cSim.pBornSum[offset] = sA[threadIdx.x].sum;
#endif
        }

        pos -= cSim.nonbond_workBlock;
    }
}
