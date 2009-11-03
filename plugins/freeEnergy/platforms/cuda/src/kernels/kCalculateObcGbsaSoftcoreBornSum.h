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

/**
 * This file contains the kernel for calculating Born sums.  It is included
 * several times in kCalculateObcGbsaBornSum.cu with different #defines to generate
 * different versions of the kernels.
 */

__global__ void METHOD_NAME(kCalculateObcGbsa, BornSum_kernel)(unsigned int* workUnit)
{
    extern __shared__ Atom sA[];
    unsigned int totalWarps   = cSim.nonbond_blocks*cSim.nonbond_threads_per_block/GRID;
    unsigned int warp         = (blockIdx.x*blockDim.x+threadIdx.x)/GRID;
    unsigned int numWorkUnits = cSim.pInteractionCount[0];
    unsigned int pos          = warp*numWorkUnits/totalWarps;
    unsigned int end          = (warp+1)*numWorkUnits/totalWarps;

#ifdef USE_CUTOFF
    float* tempBuffer = (float*) &sA[cSim.nonbond_threads_per_block];
#endif

    while (pos < end)
    {
        // Extract cell coordinates from appropriate work unit
        
        //unsigned int x = workUnit[pos + (blockIdx.x*numWorkUnits)/gridDim.x];
        unsigned int x = workUnit[pos];
        unsigned int y = ((x >> 2) & 0x7fff) << GRIDBITS;
        x = (x >> 17) << GRIDBITS;
        float       dx;
        float       dy;
        float       dz;
        float       r2;
        float       r;

        unsigned int tgx = threadIdx.x & (GRID - 1);
        unsigned int tbx = threadIdx.x - tgx;
        unsigned int tj = tgx;
        Atom* psA = &sA[tbx];

        if (x == y) // Handle diagonals uniquely at 50% efficiency
        {
            // Read fixed atom data into registers and GRF
            unsigned int i = x + tgx;
            float4 apos                             = cSim.pPosq[i];    // Local atom x, y, z, sum
            float2 ar                               = cSim.pObcData[i];   // Local atom vr, sr
            float polarScaleData                    = gbsaSimDev.pNonPolarScalingFactors[i];  // scale contribution
            sA[threadIdx.x].x                       = apos.x;
            sA[threadIdx.x].y                       = apos.y;
            sA[threadIdx.x].z                       = apos.z;
            sA[threadIdx.x].r                       = ar.x;
            sA[threadIdx.x].sr                      = ar.y;
            sA[threadIdx.x].polarScaleData          = polarScaleData;
            apos.w                                  = 0.0f;

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
                        float sum      = l_ij -
                                         u_ij +
                                         0.25f * r * (u_ij2 - l_ij2) +
                                         (0.50f * rInverse * ratio) +
                                         (0.25f * psA[j].sr * psA[j].sr * rInverse) *
                                         (l_ij2 - u_ij2);
                        float rj = psA[j].r;
                        if (ar.x < (rj - r))
                        {
                            sum += 2.0f * ((1.0f / ar.x) - l_ij);
                        }
                        apos.w +=  psA[j].polarScaleData*sum;
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
            unsigned int j                  = y + tgx;
            unsigned int i                  = x + tgx;

            float4 temp                     = cSim.pPosq[j];
            float2 temp1                    = cSim.pObcData[j];
            float polarScaleDataJ           = gbsaSimDev.pNonPolarScalingFactors[j];  // scale contribution
            float4 apos                     = cSim.pPosq[i];        // Local atom x, y, z, sum
            float2 ar                       = cSim.pObcData[i];    // Local atom vr, sr
            float polarScaleDataI           = gbsaSimDev.pNonPolarScalingFactors[i];  // scale contribution
            sA[threadIdx.x].x               = temp.x;
            sA[threadIdx.x].y               = temp.y;
            sA[threadIdx.x].z               = temp.z;
            sA[threadIdx.x].r               = temp1.x;
            sA[threadIdx.x].sr              = temp1.y;
            sA[threadIdx.x].polarScaleData  = polarScaleDataJ;
            sA[threadIdx.x].sum = apos.w    = 0.0f;

#ifdef USE_CUTOFF
            //unsigned int flags = cSim.pInteractionFlag[pos + (blockIdx.x*numWorkUnits)/gridDim.x];
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
                            float srj = psA[tj].sr;
                            float scale = psA[tj].polarScaleData;
                            if (ar.x < (srj - r))
                            {
                                term += 2.0f * ((1.0f / ar.x) - l_ij);
                            }
                            //apos.w        += term;
                            apos.w        += (scale*term);
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
                            float rj = psA[tj].r;
                            if (rj < (ar.y - r))
                            {
                                term += 2.0f * ((1.0f / psA[tj].r) - l_ij);
                            }
                            psA[tj].sum    += polarScaleDataI*term;
                        }
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
                            float rInverse          = 1.0f / r;
                            float rScaledRadiusJ    = r + psA[j].sr;
                            if (ar.x < rScaledRadiusJ)
                            {
                                float l_ij     = 1.0f / max(ar.x, fabs(r - psA[j].sr));
                                float u_ij     = 1.0f / rScaledRadiusJ;
                                float l_ij2    = l_ij * l_ij;
                                float u_ij2    = u_ij * u_ij;
                                float ratio    = log(u_ij / l_ij);
                                float term     = l_ij -
                                                 u_ij +
                                                 0.25f * r * (u_ij2 - l_ij2) +
                                                 (0.50f * rInverse * ratio) +
                                                 (0.25f * psA[j].sr * psA[j].sr * rInverse) *
                                                 (l_ij2 - u_ij2);
                                float srj = psA[j].sr;
                                if (ar.x < (srj - r))
                                {
                                    term += 2.0f * ((1.0f / ar.x) - l_ij);
                                }
                                apos.w        += psA[j].polarScaleData*term;
                            }
                            float rScaledRadiusI    = r + ar.y;
                            if (psA[j].r < rScaledRadiusI)
                            {
                                float l_ij     = 1.0f / max(psA[j].r, fabs(r - ar.y));
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
                                float rj = psA[j].r;
                                if (rj < (ar.y - r))
                                {
                                    term += 2.0f * ((1.0f / psA[j].r) - l_ij);
                                }
                                tempBuffer[threadIdx.x] = polarScaleDataI*term;
                            }
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
