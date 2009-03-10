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
 * This file contains the kernel for evalauating nonbonded forces and the first stage of GBSA.
 * It is included several times in kCalculateCDLJObcGbsaForces1.cu with different #defines to generate
 * different versions of the kernels.
 */

__global__ void METHOD_NAME(kCalculateCDLJObcGbsa, Forces1_kernel)(unsigned int* workUnit, int numWorkUnits)
{
    extern __shared__ Atom sA[];
    unsigned int totalWarps = cSim.nonbond_blocks*cSim.nonbond_threads_per_block/GRID;
    unsigned int warp = (blockIdx.x*blockDim.x+threadIdx.x)/GRID;
    int pos = warp*numWorkUnits/totalWarps;
    int end = (warp+1)*numWorkUnits/totalWarps;
#ifdef USE_CUTOFF
    float* tempBuffer = (float*) &sA[cSim.nonbond_threads_per_block];
#endif

    int lasty = -1;
    while (pos < end)
    {

        // Extract cell coordinates from appropriate work unit
        unsigned int x                      = workUnit[pos];
        unsigned int y                      = ((x >> 2) & 0x7fff) << GRIDBITS;
        bool bExclusionFlag                 = (x & 0x1);
        x                                   = (x >> 17) << GRIDBITS;
        unsigned int tgx                    = threadIdx.x & (GRID - 1);
        unsigned int i                      = x + tgx;
        float4 apos                         = cSim.pPosq[i];
        float2 a                            = cSim.pAttr[i];
        float br                            = cSim.pBornRadii[i];
        unsigned int tbx                    = threadIdx.x - tgx;
        int tj                              = tgx;
        Atom* psA                           = &sA[tbx];
        float4 af;
        af.x                        = 0.0f;
        af.y                        = 0.0f;
        af.z                        = 0.0f;
        af.w                        = 0.0f;
        if (x == y) // Handle diagonals uniquely at 50% efficiency
        {
            // Read fixed atom data into registers and GRF
            sA[threadIdx.x].x           = apos.x;
            sA[threadIdx.x].y           = apos.y;
            sA[threadIdx.x].z           = apos.z;
            sA[threadIdx.x].q           = apos.w;
            float q2                    = cSim.preFactor * apos.w;
            apos.w                     *= cSim.epsfac;
            sA[threadIdx.x].sig         = a.x;
            sA[threadIdx.x].eps         = a.y;
            sA[threadIdx.x].br          = br;
            if (!bExclusionFlag)
            {
                for (unsigned int j = 0; j < GRID; j++)
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

                    // CDLJ part
                    float invR              = 1.0f / sqrt(r2);
                    float sig               = a.x + psA[j].sig;
                    float sig2              = invR * sig;
                    sig2                   *= sig2;
                    float sig6              = sig2 * sig2 * sig2;
                    float eps               = a.y * psA[j].eps;
                    float dEdR              = eps * (12.0f * sig6 - 6.0f) * sig6;
#ifdef USE_CUTOFF
                    dEdR                   += apos.w * psA[j].q * (invR - 2.0f * cSim.reactionFieldK * r2);
#else
                    dEdR                   += apos.w * psA[j].q * invR;
#endif
                    dEdR                   *= invR * invR;

                    // ObcGbsaForce1 part
                    float alpha2_ij         = br * psA[j].br;
                    float D_ij              = r2 / (4.0f * alpha2_ij);
                    float expTerm           = exp(-D_ij);
                    float denominator2      = r2 + alpha2_ij * expTerm;
                    float denominator       = sqrt(denominator2);
                    float Gpol              = (q2 * psA[j].q) / (denominator * denominator2);
                    float dGpol_dalpha2_ij  = -0.5f * Gpol * expTerm * (1.0f + D_ij);
                    af.w                   += dGpol_dalpha2_ij * psA[j].br;
                    dEdR                   += Gpol * (1.0f - 0.25f * expTerm);
#ifdef USE_CUTOFF
                    if (r2 > cSim.nonbondedCutoffSqr)
                    {
                        dEdR = 0.0f;
                    }
#endif

                    // Add Forces
                    dx                     *= dEdR;
                    dy                     *= dEdR;
                    dz                     *= dEdR;
                    af.x                   -= dx;
                    af.y                   -= dy;
                    af.z                   -= dz;
                }
            }
            else  // bExclusion
            {
                unsigned int xi   = x>>GRIDBITS;
                int cell          = xi+xi*cSim.paddedNumberOfAtoms/GRID-xi*(xi+1)/2;
                unsigned int excl = cSim.pExclusion[cSim.pExclusionIndex[cell]+tgx];
                for (unsigned int j = 0; j < GRID; j++)
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

                    // CDLJ part
                    float invR              = 1.0f / sqrt(r2);
                    float sig               = a.x + psA[j].sig;
                    float sig2              = invR * sig;
                    sig2                   *= sig2;
                    float sig6              = sig2 * sig2 * sig2;
                    float eps               = a.y * psA[j].eps;
                    float dEdR              = eps * (12.0f * sig6 - 6.0f) * sig6;
#ifdef USE_CUTOFF
                    dEdR                   += apos.w * psA[j].q * (invR - 2.0f * cSim.reactionFieldK * r2);
#else
                    dEdR                   += apos.w * psA[j].q * invR;
#endif
                    dEdR                   *= invR * invR;
                    if (!(excl & 0x1))
                    {
                        dEdR = 0.0f;
                    }

                    // ObcGbsaForce1 part
                    float alpha2_ij         = br * psA[j].br;
                    float D_ij              = r2 / (4.0f * alpha2_ij);
                    float expTerm           = exp(-D_ij);
                    float denominator2      = r2 + alpha2_ij * expTerm;
                    float denominator       = sqrt(denominator2);
                    float Gpol              = (q2 * psA[j].q) / (denominator * denominator2);
                    float dGpol_dalpha2_ij  = -0.5f * Gpol * expTerm * (1.0f + D_ij);
                    af.w                   += dGpol_dalpha2_ij * psA[j].br;
                    dEdR                   += Gpol * (1.0f - 0.25f * expTerm);
#if defined USE_PERIODIC
                    if (i >= cSim.atoms || x+j >= cSim.atoms || r2 > cSim.nonbondedCutoffSqr)
                    {
                        dEdR = 0.0f;
                    }
#elif defined USE_CUTOFF
                    if (r2 > cSim.nonbondedCutoffSqr)
                    {
                        dEdR = 0.0f;
                    }
#endif

                    // Add Forces
                    dx                     *= dEdR;
                    dy                     *= dEdR;
                    dz                     *= dEdR;
                    af.x                   -= dx;
                    af.y                   -= dy;
                    af.z                   -= dz;
                    excl                  >>= 1;
                }
            }

            // Write results
#ifdef USE_OUTPUT_BUFFER_PER_WARP
            int offset                  = x + tgx + warp*cSim.stride;
            float4 of                   = cSim.pForce4a[offset];
            of.x                       += af.x;
            of.y                       += af.y;
            of.z                       += af.z;
            of.w                       += af.w;
            cSim.pForce4a[offset]       = of;
            cSim.pBornForce[offset]     = af.w;
#else
            int offset                  = x + tgx + (x >> GRIDBITS) * cSim.stride;
            cSim.pForce4a[offset]       = af;
            cSim.pBornForce[offset]     = af.w;
#endif
        }
        else        // 100% utilization
        {
            // Read fixed atom data into registers and GRF
            if (lasty != y)
            {
                int j                       = y + tgx;
                float4 temp                 = cSim.pPosq[j];
                float2 temp1                = cSim.pAttr[j];
                sA[threadIdx.x].br          = cSim.pBornRadii[j];
                sA[threadIdx.x].x           = temp.x;
                sA[threadIdx.x].y           = temp.y;
                sA[threadIdx.x].z           = temp.z;
                sA[threadIdx.x].q           = temp.w;
                sA[threadIdx.x].sig         = temp1.x;
                sA[threadIdx.x].eps         = temp1.y;
            }
            sA[threadIdx.x].fx          = 0.0f;
            sA[threadIdx.x].fy          = 0.0f;
            sA[threadIdx.x].fz          = 0.0f;
            sA[threadIdx.x].fb          = 0.0f;
            float q2                    = apos.w * cSim.preFactor;
            apos.w                     *= cSim.epsfac;
            if (!bExclusionFlag)
            {
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

                        // CDLJ part
                        float invR              = 1.0f / sqrt(r2);
                        float sig               = a.x + psA[tj].sig;
                        float sig2              = invR * sig;
                        sig2                   *= sig2;
                        float sig6              = sig2 * sig2 * sig2;
                        float eps               = a.y * psA[tj].eps;
                        float dEdR              = eps * (12.0f * sig6 - 6.0f) * sig6;
#ifdef USE_CUTOFF
                        dEdR                   += apos.w * psA[tj].q * (invR - 2.0f * cSim.reactionFieldK * r2);
#else
                        dEdR                   += apos.w * psA[tj].q * invR;
#endif
                        dEdR                   *= invR * invR;

                        // ObcGbsaForce1 part
                        float alpha2_ij         = br * psA[tj].br;
                        float D_ij              = r2 / (4.0f * alpha2_ij);
                        float expTerm           = exp(-D_ij);
                        float denominator2      = r2 + alpha2_ij * expTerm;
                        float denominator       = sqrt(denominator2);
                        float Gpol              = (q2 * psA[tj].q) / (denominator * denominator2);
                        float dGpol_dalpha2_ij  = -0.5f * Gpol * expTerm * (1.0f + D_ij);
                        af.w                   += dGpol_dalpha2_ij * psA[tj].br;
                        psA[tj].fb             += dGpol_dalpha2_ij * br;
                        dEdR                   += Gpol * (1.0f - 0.25f * expTerm);
#ifdef USE_CUTOFF
                        if (r2 > cSim.nonbondedCutoffSqr)
                        {
                            dEdR = 0.0f;
                        }
#endif

                        // Add forces
                        dx                     *= dEdR;
                        dy                     *= dEdR;
                        dz                     *= dEdR;
                        af.x                   -= dx;
                        af.y                   -= dy;
                        af.z                   -= dz;
                        psA[tj].fx             += dx;
                        psA[tj].fy             += dy;
                        psA[tj].fz             += dz;
                        tj                      = (tj + 1) & (GRID - 1);
                    }
                }
#ifdef USE_CUTOFF
                else
                {
                    // Compute only a subset of the interactions in this block.

                    for (int j = 0; j < GRID; j++)
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

                            // CDLJ part
                            float invR              = 1.0f / sqrt(r2);
                            float sig               = a.x + psA[j].sig;
                            float sig2              = invR * sig;
                            sig2                   *= sig2;
                            float sig6              = sig2 * sig2 * sig2;
                            float eps               = a.y * psA[j].eps;
                            float dEdR              = eps * (12.0f * sig6 - 6.0f) * sig6;
#ifdef USE_CUTOFF
                            dEdR                   += apos.w * psA[j].q * (invR - 2.0f * cSim.reactionFieldK * r2);
#else
                            dEdR                   += apos.w * psA[j].q * invR;
#endif
                            dEdR                   *= invR * invR;

                            // ObcGbsaForce1 part
                            float alpha2_ij         = br * psA[j].br;
                            float D_ij              = r2 / (4.0f * alpha2_ij);
                            float expTerm           = exp(-D_ij);
                            float denominator2      = r2 + alpha2_ij * expTerm;
                            float denominator       = sqrt(denominator2);
                            float Gpol              = (q2 * psA[j].q) / (denominator * denominator2);
                            float dGpol_dalpha2_ij  = -0.5f * Gpol * expTerm * (1.0f + D_ij);
                            af.w                   += dGpol_dalpha2_ij * psA[j].br;
                            dEdR                   += Gpol * (1.0f - 0.25f * expTerm);

                            // Sum the Born forces.

                            tempBuffer[threadIdx.x] = dGpol_dalpha2_ij * br;
                            if (tgx % 2 == 0)
                                tempBuffer[threadIdx.x] += tempBuffer[threadIdx.x+1];
                            if (tgx % 4 == 0)
                                tempBuffer[threadIdx.x] += tempBuffer[threadIdx.x+2];
                            if (tgx % 8 == 0)
                                tempBuffer[threadIdx.x] += tempBuffer[threadIdx.x+4];
                            if (tgx % 16 == 0)
                                tempBuffer[threadIdx.x] += tempBuffer[threadIdx.x+8];
                            if (tgx == 0)
                                psA[j].fb += tempBuffer[threadIdx.x] + tempBuffer[threadIdx.x+16];
#ifdef USE_CUTOFF
                            if (r2 > cSim.nonbondedCutoffSqr)
                            {
                                dEdR = 0.0f;
                            }
#endif

                            // Add forces
                            dx                     *= dEdR;
                            dy                     *= dEdR;
                            dz                     *= dEdR;
                            af.x                   -= dx;
                            af.y                   -= dy;
                            af.z                   -= dz;
                            tempBuffer[threadIdx.x] = dx;
                            if (tgx % 2 == 0)
                                tempBuffer[threadIdx.x] += tempBuffer[threadIdx.x+1];
                            if (tgx % 4 == 0)
                                tempBuffer[threadIdx.x] += tempBuffer[threadIdx.x+2];
                            if (tgx % 8 == 0)
                                tempBuffer[threadIdx.x] += tempBuffer[threadIdx.x+4];
                            if (tgx % 16 == 0)
                                tempBuffer[threadIdx.x] += tempBuffer[threadIdx.x+8];
                            if (tgx == 0)
                                psA[j].fx += tempBuffer[threadIdx.x] + tempBuffer[threadIdx.x+16];
                            tempBuffer[threadIdx.x] = dy;
                            if (tgx % 2 == 0)
                                tempBuffer[threadIdx.x] += tempBuffer[threadIdx.x+1];
                            if (tgx % 4 == 0)
                                tempBuffer[threadIdx.x] += tempBuffer[threadIdx.x+2];
                            if (tgx % 8 == 0)
                                tempBuffer[threadIdx.x] += tempBuffer[threadIdx.x+4];
                            if (tgx % 16 == 0)
                                tempBuffer[threadIdx.x] += tempBuffer[threadIdx.x+8];
                            if (tgx == 0)
                                psA[j].fy += tempBuffer[threadIdx.x] + tempBuffer[threadIdx.x+16];
                            tempBuffer[threadIdx.x] = dz;
                            if (tgx % 2 == 0)
                                tempBuffer[threadIdx.x] += tempBuffer[threadIdx.x+1];
                            if (tgx % 4 == 0)
                                tempBuffer[threadIdx.x] += tempBuffer[threadIdx.x+2];
                            if (tgx % 8 == 0)
                                tempBuffer[threadIdx.x] += tempBuffer[threadIdx.x+4];
                            if (tgx % 16 == 0)
                                tempBuffer[threadIdx.x] += tempBuffer[threadIdx.x+8];
                            if (tgx == 0)
                                psA[j].fz += tempBuffer[threadIdx.x] + tempBuffer[threadIdx.x+16];
                        }
                    }
                }
#endif
            }
            else  // bExclusion
            {
                unsigned int xi   = x>>GRIDBITS;
                unsigned int yi   = y>>GRIDBITS;
                int cell          = xi+yi*cSim.paddedNumberOfAtoms/GRID-yi*(yi+1)/2;
                unsigned int excl = cSim.pExclusion[cSim.pExclusionIndex[cell]+tgx];
                excl              = (excl >> tgx) | (excl << (GRID - tgx));
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

                    // CDLJ part
                    float invR              = 1.0f / sqrt(r2);
                    float sig               = a.x + psA[tj].sig;
                    float sig2              = invR * sig;
                    sig2                   *= sig2;
                    float sig6              = sig2 * sig2 * sig2;
                    float eps               = a.y * psA[tj].eps;
                    float dEdR              = eps * (12.0f * sig6 - 6.0f) * sig6;
#ifdef USE_CUTOFF
                    dEdR                   += apos.w * psA[tj].q * (invR - 2.0f * cSim.reactionFieldK * r2);
#else
                    dEdR                   += apos.w * psA[tj].q * invR;
#endif
                    dEdR                   *= invR * invR;
                    if (!(excl & 0x1))
                    {
                        dEdR = 0.0f;
                    }

                    // ObcGbsaForce1 part
                    float alpha2_ij         = br * psA[tj].br;
                    float D_ij              = r2 / (4.0f * alpha2_ij);
                    float expTerm           = exp(-D_ij);
                    float denominator2      = r2 + alpha2_ij * expTerm;
                    float denominator       = sqrt(denominator2);
                    float Gpol              = (q2 * psA[tj].q) / (denominator * denominator2);
                    float dGpol_dalpha2_ij  = -0.5f * Gpol * expTerm * (1.0f + D_ij);
                    af.w                   += dGpol_dalpha2_ij * psA[tj].br;
                    psA[tj].fb             += dGpol_dalpha2_ij * br;
                    dEdR                   += Gpol * (1.0f - 0.25f * expTerm);
#if defined USE_PERIODIC
                    if (i >= cSim.atoms || y+tj >= cSim.atoms || r2 > cSim.nonbondedCutoffSqr)
                    {
                        dEdR = 0.0f;
                    }
#elif defined USE_CUTOFF
                    if (r2 > cSim.nonbondedCutoffSqr)
                    {
                        dEdR = 0.0f;
                    }
#endif

                    // Add forces
                    dx                     *= dEdR;
                    dy                     *= dEdR;
                    dz                     *= dEdR;
                    af.x                   -= dx;
                    af.y                   -= dy;
                    af.z                   -= dz;
                    psA[tj].fx             += dx;
                    psA[tj].fy             += dy;
                    psA[tj].fz             += dz;
                    excl                  >>= 1;
                    tj                      = (tj + 1) & (GRID - 1);
                }
            }

            // Write results
#ifdef USE_OUTPUT_BUFFER_PER_WARP
            int offset                  = x + tgx + warp*cSim.stride;
            float4 of                   = cSim.pForce4a[offset];
            of.x                       += af.x;
            of.y                       += af.y;
            of.z                       += af.z;
            of.w                       += af.w;
            cSim.pForce4a[offset]       = of;
            cSim.pBornForce[offset]     = af.w;
            offset                      = y + tgx + warp*cSim.stride;
            of                          = cSim.pForce4a[offset];
            of.x                       += sA[threadIdx.x].fx;
            of.y                       += sA[threadIdx.x].fy;
            of.z                       += sA[threadIdx.x].fz;
            of.w                       += sA[threadIdx.x].fb;
            cSim.pForce4a[offset]       = of;
            cSim.pBornForce[offset]     = af.w;
#else
            int offset                  = x + tgx + (y >> GRIDBITS) * cSim.stride;
            cSim.pForce4a[offset]       = af;
            cSim.pBornForce[offset]     = af.w;
            af.x                        = sA[threadIdx.x].fx;
            af.y                        = sA[threadIdx.x].fy;
            af.z                        = sA[threadIdx.x].fz;
            af.w                        = sA[threadIdx.x].fb;
            offset                      = y + tgx + (x >> GRIDBITS) * cSim.stride;
            cSim.pForce4a[offset]       = af;
            cSim.pBornForce[offset]     = af.w;
#endif
            lasty = y;
        }
        pos++;
    }
}
