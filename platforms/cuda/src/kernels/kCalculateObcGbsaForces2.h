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
 * This file contains the kernel for evalauating the second stage of GBSA.  It is included
 * several times in kCalculateObcGbsaForces2.cu with different #defines to generate
 * different versions of the kernels.
 */

__global__ void METHOD_NAME(kCalculateObcGbsa, Forces2_kernel)(unsigned int* workUnit, int numWorkUnits)
{
    extern __shared__ Atom sA[];
    unsigned int totalWarps = cSim.nonbond_blocks*cSim.nonbond_threads_per_block/GRID;
    unsigned int warp = (blockIdx.x*blockDim.x+threadIdx.x)/GRID;
    int pos = warp*numWorkUnits/totalWarps;
    int end = (warp+1)*numWorkUnits/totalWarps;

    int lasty = -1;
    while (pos < end)
    {

        // Extract cell coordinates from appropriate work unit
        unsigned int x                  = workUnit[pos];
        unsigned int y                  = ((x >> 2) & 0x7fff) << GRIDBITS;
        x                               = (x >> 17) << GRIDBITS;
        unsigned int tgx                = threadIdx.x & (GRID - 1);
        unsigned int i                  = x + tgx;
        float4 apos                     = cSim.pPosq[i];
        float2 a                        = cSim.pObcData[i];
        float fb                        = cSim.pBornForce[i];
        unsigned int tbx                = threadIdx.x - tgx;
        int tj                          = tgx;
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
            sA[threadIdx.x].r           = a.x;
            sA[threadIdx.x].sr          = a.y;
            sA[threadIdx.x].sr2         = a.y * a.y;
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
                float rScaledRadiusJ    = r + psA[j].sr;

                float l_ij          = 1.0f / max(a.x, fabs(r - psA[j].sr));
                float u_ij          = 1.0f / rScaledRadiusJ;
                float rInverse      = 1.0f / r;
                float l_ij2         = l_ij * l_ij;
                float u_ij2         = u_ij * u_ij;
                float r2Inverse     = rInverse * rInverse;
                float t1            = log (u_ij / l_ij);
                float t2            = (l_ij2 - u_ij2);
                float t3            = t2 * rInverse;
                t1                 *= rInverse;

                // Born Forces term
                float term          =  0.125f *
                                      (1.000f + psA[j].sr2 * r2Inverse) * t3 +
                                       0.250f * t1 * r2Inverse;
                float dE            = fb * term;

#if defined USE_PERIODIC
                if (a.x >= rScaledRadiusJ || i >= cSim.atoms || x+j >= cSim.atoms || r2 > cSim.nonbondedCutoffSqr)
#elif defined USE_CUTOFF
                if (a.x >= rScaledRadiusJ || r2 > cSim.nonbondedCutoffSqr)
#else
                if (a.x >= rScaledRadiusJ)
#endif
                {
                    dE              = 0.0f;
                }
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
            int offset                  = x + tgx + warp*cSim.stride;
            of                          = cSim.pForce4b[offset];
            of.x                       += af.x + sA[threadIdx.x].fx;
            of.y                       += af.y + sA[threadIdx.x].fy;
            of.z                       += af.z + sA[threadIdx.x].fz;
            cSim.pForce4b[offset]       = of;
#else
            int offset                  = x + tgx + (x >> GRIDBITS) * cSim.stride;
            of.x                        = af.x + sA[threadIdx.x].fx;
            of.y                        = af.y + sA[threadIdx.x].fy;
            of.z                        = af.z + sA[threadIdx.x].fz;
            of.w                        = 0.0f;
            cSim.pForce4b[offset]       = of;
#endif
        }
        else
        {
            // Read fixed atom data into registers and GRF
            if (lasty != y)
            {
                int j                       = y + tgx;
                float4 temp                 = cSim.pPosq[j];
                float2 temp1                = cSim.pObcData[j];
                sA[threadIdx.x].fb          = cSim.pBornForce[j];
                sA[threadIdx.x].x           = temp.x;
                sA[threadIdx.x].y           = temp.y;
                sA[threadIdx.x].z           = temp.z;
                sA[threadIdx.x].r           = temp1.x;
                sA[threadIdx.x].sr          = temp1.y;
                sA[threadIdx.x].sr2         = temp1.y * temp1.y;
            }
            float sr2                   = a.y * a.y;
            for (int j = 0; j < GRID; j++)
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

                // Interleaved Atom I and J Born Forces and sum components
                float r2Inverse         = 1.0f / r2;
                float rScaledRadiusJ    = r + psA[tj].sr;
                float rScaledRadiusI    = r + a.y;
                float rInverse          = 1.0f / r;
                float l_ijJ             = 1.0f / max(a.x, fabs(r - psA[tj].sr));
                float l_ijI             = 1.0f / max(psA[tj].r, fabs(r - a.y));
                float u_ijJ             = 1.0f / rScaledRadiusJ;
                float u_ijI             = 1.0f / rScaledRadiusI;
                float l_ij2J            = l_ijJ * l_ijJ;
                float l_ij2I            = l_ijI * l_ijI;
                float u_ij2J            = u_ijJ * u_ijJ;
                float u_ij2I            = u_ijI * u_ijI;
                float t1J               = log (u_ijJ / l_ijJ);
                float t1I               = log (u_ijI / l_ijI);
                float t2J               = (l_ij2J - u_ij2J);
                float t2I               = (l_ij2I - u_ij2I);
                float t3J               = t2J * rInverse;
                float t3I               = t2I * rInverse;
                t1J                    *= rInverse;
                t1I                    *= rInverse;

                // Born Forces term
                float term              =  0.125f *
                                          (1.000f + psA[tj].sr2 * r2Inverse) * t3J +
                                           0.250f * t1J * r2Inverse;
                float dE                = fb * term;

#if defined USE_PERIODIC
                if (a.x >= rScaledRadiusJ || i >= cSim.atoms || y+tj >= cSim.atoms || r2 > cSim.nonbondedCutoffSqr)
#elif defined USE_CUTOFF
                if (a.x >= rScaledRadiusJ || r2 > cSim.nonbondedCutoffSqr)
#else
                if (a.x >= rScaledRadiusJ)
#endif
                {
                    dE                  = 0.0f;
                }

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
                term                    =  0.125f *
                                          (1.000f + sr2 * r2Inverse) * t3I +
                                           0.250f * t1I * r2Inverse;
                dE                      = psA[tj].fb * term;

#ifdef USE_PERIODIC
                if (psA[tj].r >= rScaledRadiusI || i >= cSim.atoms || y+tj >= cSim.atoms || r2 > cSim.nonbondedCutoffSqr)
#elif defined USE_CUTOFF
                if (psA[tj].r >= rScaledRadiusI || r2 > cSim.nonbondedCutoffSqr)
#else
                if (psA[tj].r >= rScaledRadiusI)
#endif
                {
                    dE                  = 0.0f;
                }
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

            // Write results
            float4 of;
#ifdef USE_OUTPUT_BUFFER_PER_WARP
            int offset                  = x + tgx + warp*cSim.stride;
            of                          = cSim.pForce4b[offset];
            of.x                       += af.x;
            of.y                       += af.y;
            of.z                       += af.z;
            cSim.pForce4b[offset]       = of;
            offset                      = y + tgx + warp*cSim.stride;
            of                          = cSim.pForce4b[offset];
            of.x                       += sA[threadIdx.x].fx;
            of.y                       += sA[threadIdx.x].fy;
            of.z                       += sA[threadIdx.x].fz;
            cSim.pForce4b[offset]       = of;
#else
            int offset                  = x + tgx + (y >> GRIDBITS) * cSim.stride;
            of.x                        = af.x;
            of.y                        = af.y;
            of.z                        = af.z;
            of.w                        = 0.0f;
            cSim.pForce4b[offset]       = of;
            offset                      = y + tgx + (x >> GRIDBITS) * cSim.stride;
            of.x                        = sA[threadIdx.x].fx;
            of.y                        = sA[threadIdx.x].fy;
            of.z                        = sA[threadIdx.x].fz;
            cSim.pForce4b[offset]       = of;
#endif
        }
        lasty = y;
        pos++;
    }
}
