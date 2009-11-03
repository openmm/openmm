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
 * This file contains the kernel for evalauating the second stage of GBSA.  It is included
 * several times in kCalculateObcGbsaSoftcoreForces2.cu with different #defines to generate
 * different versions of the kernels.
 */

__global__ void METHOD_NAME(kCalculateObcGbsaSoftcore, Forces2_kernel)(unsigned int* workUnit)
{
    extern __shared__ Atom sA[];
    unsigned int totalWarps = cSim.bornForce2_blocks*cSim.bornForce2_threads_per_block/GRID;
    unsigned int warp = (blockIdx.x*blockDim.x+threadIdx.x)/GRID;
    unsigned int numWorkUnits = cSim.pInteractionCount[0];
    unsigned int pos = warp*numWorkUnits/totalWarps;
    unsigned int end = (warp+1)*numWorkUnits/totalWarps;
#ifdef USE_CUTOFF
    float3* tempBuffer = (float3*) &sA[cSim.bornForce2_threads_per_block];
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
        float2 a                        = cSim.pObcData[i];
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
            sA[threadIdx.x].r           = a.x;
            sA[threadIdx.x].sr          = a.y;
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
                                      (1.000f + psA[j].sr * psA[j].sr * r2Inverse) * t3 +
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
            unsigned int offset         = x + tgx + warp*cSim.stride;
#else
            unsigned int offset         = x + tgx + (x >> GRIDBITS) * cSim.stride;
#endif
            of                          = cSim.pForce4b[offset];
            of.x                       += af.x + sA[threadIdx.x].fx;
            of.y                       += af.y + sA[threadIdx.x].fy;
            of.z                       += af.z + sA[threadIdx.x].fz;
            of.w                        = 0.0f;
            cSim.pForce4b[offset]       = of;
        }
        else
        {
            // Read fixed atom data into registers and GRF
            if (lasty != y)
            {
                unsigned int j              = y + tgx;
                float4 temp                 = cSim.pPosq[j];
                float2 temp1                = cSim.pObcData[j];
                sA[threadIdx.x].fb          = cSim.pBornForce[j];
                sA[threadIdx.x].x           = temp.x;
                sA[threadIdx.x].y           = temp.y;
                sA[threadIdx.x].z           = temp.z;
                sA[threadIdx.x].r           = temp1.x;
                sA[threadIdx.x].sr          = temp1.y;
            }
            float sr2                   = a.y * a.y;
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
                                              (1.000f + psA[tj].sr * psA[tj].sr * r2Inverse) * t3J +
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

                    float rj = psA[tj].r;
#ifdef USE_PERIODIC
                    if (rj >= rScaledRadiusI || i >= cSim.atoms || y+tj >= cSim.atoms || r2 > cSim.nonbondedCutoffSqr)
#elif defined USE_CUTOFF
                    if (rj >= rScaledRadiusI || r2 > cSim.nonbondedCutoffSqr)
#else
                    if (rj >= rScaledRadiusI)
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
                        float r2Inverse         = 1.0f / r2;
                        float rScaledRadiusJ    = r + psA[j].sr;
                        float rScaledRadiusI    = r + a.y;
                        float rInverse          = 1.0f / r;
                        float l_ijJ             = 1.0f / max(a.x, fabs(r - psA[j].sr));
                        float l_ijI             = 1.0f / max(psA[j].r, fabs(r - a.y));
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
                                                  (1.000f + psA[j].sr * psA[j].sr * r2Inverse) * t3J +
                                                   0.250f * t1J * r2Inverse;
                        float dE                = fb * term;

    #if defined USE_PERIODIC
                        if (a.x >= rScaledRadiusJ || i >= cSim.atoms || y+j >= cSim.atoms || r2 > cSim.nonbondedCutoffSqr)
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
                        tempBuffer[threadIdx.x].x = d;
                        d                       = dy * dE;
                        af.y                   -= d;
                        tempBuffer[threadIdx.x].y = d;
                        d                       = dz * dE;
                        af.z                   -= d;
                        tempBuffer[threadIdx.x].z = d;

                        // Atom J Born sum term
                        term                    =  0.125f *
                                                  (1.000f + sr2 * r2Inverse) * t3I +
                                                   0.250f * t1I * r2Inverse;
                        dE                      = psA[j].fb * term;

                        float rj = psA[j].r;
    #ifdef USE_PERIODIC
                        if (rj >= rScaledRadiusI || i >= cSim.atoms || y+j >= cSim.atoms || r2 > cSim.nonbondedCutoffSqr)
    #elif defined USE_CUTOFF
                        if (rj >= rScaledRadiusI || r2 > cSim.nonbondedCutoffSqr)
    #else
                        if (rj >= rScaledRadiusI)
    #endif
                        {
                            dE                  = 0.0f;
                        }
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
            of                          = cSim.pForce4b[offset];
            of.x                       += af.x;
            of.y                       += af.y;
            of.z                       += af.z;
            cSim.pForce4b[offset]       = of;

#ifdef USE_OUTPUT_BUFFER_PER_WARP
            offset                      = y + tgx + warp*cSim.stride;
#else
            offset                      = y + tgx + (x >> GRIDBITS) * cSim.stride;
#endif
            of                          = cSim.pForce4b[offset];
            of.x                       += sA[threadIdx.x].fx;
            of.y                       += sA[threadIdx.x].fy;
            of.z                       += sA[threadIdx.x].fz;
            cSim.pForce4b[offset]       = of;
        }
        lasty = y;
        pos++;
    }
}
