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
 * This file contains the kernels for evalauating nonbonded forces.  It is included
 * several times in kCalculateCDLJForces.cu with different #defines to generate
 * different versions of the kernels.
 */

__global__ void METHOD_NAME(kCalculateCDLJ, Forces_kernel)(unsigned int* workUnit, int numWorkUnits)
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
        unsigned int x = workUnit[pos];
        unsigned int y = ((x >> 2) & 0x7fff) << GRIDBITS;
        bool bExclusionFlag = (x & 0x1);
        x = (x >> 17) << GRIDBITS;
        float4      apos;   // Local atom x, y, z, q
        float3      af;     // Local atom fx, fy, fz
        float dx;
        float dy;
        float dz;
        float r2;
        float invR;
        float sig;
        float sig2;
        float sig6;
        float eps;
        float dEdR;
        unsigned int tgx = threadIdx.x & (GRID - 1);
        unsigned int tbx = threadIdx.x - tgx;
        int tj = tgx;
        Atom* psA = &sA[tbx];
        unsigned int i      = x + tgx;
        apos                = cSim.pPosq[i];
        float2 a            = cSim.pAttr[i];
        af.x                = 0.0f;
        af.y                = 0.0f;
        af.z                = 0.0f;
        if (x == y) // Handle diagonals uniquely at 50% efficiency
        {
            // Read fixed atom data into registers and GRF
            sA[threadIdx.x].x   = apos.x;
            sA[threadIdx.x].y   = apos.y;
            sA[threadIdx.x].z   = apos.z;
            sA[threadIdx.x].q   = apos.w;
            sA[threadIdx.x].sig = a.x;
            sA[threadIdx.x].eps = a.y;
            apos.w             *= cSim.epsfac;
            if (!bExclusionFlag)
            {
                for (unsigned int j = 0; j < GRID; j++)
                {
                    dx              = psA[j].x - apos.x;
                    dy              = psA[j].y - apos.y;
                    dz              = psA[j].z - apos.z;
#ifdef USE_PERIODIC
                    dx -= floor(dx/cSim.periodicBoxSizeX+0.5f)*cSim.periodicBoxSizeX;
                    dy -= floor(dy/cSim.periodicBoxSizeY+0.5f)*cSim.periodicBoxSizeY;
                    dz -= floor(dz/cSim.periodicBoxSizeZ+0.5f)*cSim.periodicBoxSizeZ;
#endif
                    r2              = dx * dx + dy * dy + dz * dz;
                    invR            = 1.0f / sqrt(r2);
                    sig             = a.x + psA[j].sig;
                    sig2            = invR * sig;
                    sig2           *= sig2;
                    sig6            = sig2 * sig2 * sig2;
                    eps             = a.y * psA[j].eps;
                    dEdR            = eps * (12.0f * sig6 - 6.0f) * sig6;
#ifdef USE_CUTOFF
                    dEdR           += apos.w * psA[j].q * (invR - 2.0f * cSim.reactionFieldK * r2);
#else
                    dEdR           += apos.w * psA[j].q * invR;
#endif
                    dEdR           *= invR * invR;
#ifdef USE_CUTOFF
                    if (r2 > cSim.nonbondedCutoffSqr)
                    {
                        dEdR = 0.0f;
                    }
#endif
                    dx             *= dEdR;
                    dy             *= dEdR;
                    dz             *= dEdR;
                    af.x           -= dx;
                    af.y           -= dy;
                    af.z           -= dz;
                }
            }
            else  // bExclusion
            {
                unsigned int excl = cSim.pExclusion[x * cSim.exclusionStride + y + tgx];
                for (unsigned int j = 0; j < GRID; j++)
                {
                    dx              = psA[j].x - apos.x;
                    dy              = psA[j].y - apos.y;
                    dz              = psA[j].z - apos.z;
#ifdef USE_PERIODIC
                    dx -= floor(dx/cSim.periodicBoxSizeX+0.5f)*cSim.periodicBoxSizeX;
                    dy -= floor(dy/cSim.periodicBoxSizeY+0.5f)*cSim.periodicBoxSizeY;
                    dz -= floor(dz/cSim.periodicBoxSizeZ+0.5f)*cSim.periodicBoxSizeZ;
#endif
                    r2              = dx * dx + dy * dy + dz * dz;
                    invR            = 1.0f / sqrt(r2);
                    sig             = a.x + psA[j].sig;
                    sig2            = invR * sig;
                    sig2           *= sig2;
                    sig6            = sig2 * sig2 * sig2;
                    eps             = a.y * psA[j].eps;
                    dEdR            = eps * (12.0f * sig6 - 6.0f) * sig6;
#ifdef USE_CUTOFF
                    dEdR           += apos.w * psA[j].q * (invR - 2.0f * cSim.reactionFieldK * r2);
#else
                    dEdR           += apos.w * psA[j].q * invR;
#endif
                    dEdR           *= invR * invR;
#ifdef USE_CUTOFF
                    if (!(excl & 0x1) || r2 > cSim.nonbondedCutoffSqr)
#else
                    if (!(excl & 0x1))
#endif
                    {
                        dEdR = 0.0f;
                    }
                    dx             *= dEdR;
                    dy             *= dEdR;
                    dz             *= dEdR;
                    af.x           -= dx;
                    af.y           -= dy;
                    af.z           -= dz;
                    excl          >>= 1;
                }
            }

            // Write results
            float4 of;
#ifdef USE_OUTPUT_BUFFER_PER_WARP
            int offset                          = x + tgx + warp*cSim.stride;
            of                                  = cSim.pForce4a[offset];
            of.x                               += af.x;
            of.y                               += af.y;
            of.z                               += af.z;
            cSim.pForce4a[offset]               = of;
#else
            of.x                                = af.x;
            of.y                                = af.y;
            of.z                                = af.z;
            of.w                                = 0.0f;
            int offset                          = x + tgx + (x >> GRIDBITS) * cSim.stride;
            cSim.pForce4a[offset]               = of;
#endif
        }
        else        // 100% utilization
        {
            // Read fixed atom data into registers and GRF
            if (lasty != y)
            {
                int j                   = y + tgx;
                float4 temp             = cSim.pPosq[j];
                float2 temp1            = cSim.pAttr[j];
                sA[threadIdx.x].x       = temp.x;
                sA[threadIdx.x].y       = temp.y;
                sA[threadIdx.x].z       = temp.z;
                sA[threadIdx.x].q       = temp.w;
                sA[threadIdx.x].sig     = temp1.x;
                sA[threadIdx.x].eps     = temp1.y;
            }
            sA[threadIdx.x].fx      = 0.0f;
            sA[threadIdx.x].fy      = 0.0f;
            sA[threadIdx.x].fz      = 0.0f;
            apos.w                 *= cSim.epsfac;
            if (!bExclusionFlag)
            {
                for (unsigned int j = 0; j < GRID; j++)
                {
                    dx              = psA[tj].x - apos.x;
                    dy              = psA[tj].y - apos.y;
                    dz              = psA[tj].z - apos.z;
#ifdef USE_PERIODIC
                    dx -= floor(dx/cSim.periodicBoxSizeX+0.5f)*cSim.periodicBoxSizeX;
                    dy -= floor(dy/cSim.periodicBoxSizeY+0.5f)*cSim.periodicBoxSizeY;
                    dz -= floor(dz/cSim.periodicBoxSizeZ+0.5f)*cSim.periodicBoxSizeZ;
#endif
                    r2              = dx * dx + dy * dy + dz * dz;
                    invR            = 1.0f / sqrt(r2);
                    sig             = a.x + psA[tj].sig;
                    sig2            = invR * sig;
                    sig2           *= sig2;
                    sig6            = sig2 * sig2 * sig2;
                    eps             = a.y * psA[tj].eps;
                    dEdR            = eps * (12.0f * sig6 - 6.0f) * sig6;
#ifdef USE_CUTOFF
                    dEdR           += apos.w * psA[tj].q * (invR - 2.0f * cSim.reactionFieldK * r2);
#else
                    dEdR           += apos.w * psA[tj].q * invR;
#endif
                    dEdR           *= invR * invR;
#ifdef USE_CUTOFF
                    if (r2 > cSim.nonbondedCutoffSqr)
                    {
                        dEdR = 0.0f;
                    }
#endif
                    dx             *= dEdR;
                    dy             *= dEdR;
                    dz             *= dEdR;
                    af.x           -= dx;
                    af.y           -= dy;
                    af.z           -= dz;
                    psA[tj].fx     += dx;
                    psA[tj].fy     += dy;
                    psA[tj].fz     += dz;
                    tj              = (tj + 1) & (GRID - 1);
                }
            }
            else  // bExclusion
            {
                // Read fixed atom data into registers and GRF
                unsigned int excl       = cSim.pExclusion[x * cSim.exclusionStride + y + tgx];
                excl                    = (excl >> tgx) | (excl << (GRID - tgx));
                for (unsigned int j = 0; j < GRID; j++)
                {
                    dx              = psA[tj].x - apos.x;
                    dy              = psA[tj].y - apos.y;
                    dz              = psA[tj].z - apos.z;
#ifdef USE_PERIODIC
                    dx -= floor(dx/cSim.periodicBoxSizeX+0.5f)*cSim.periodicBoxSizeX;
                    dy -= floor(dy/cSim.periodicBoxSizeY+0.5f)*cSim.periodicBoxSizeY;
                    dz -= floor(dz/cSim.periodicBoxSizeZ+0.5f)*cSim.periodicBoxSizeZ;
#endif
                    r2              = dx * dx + dy * dy + dz * dz;
                    invR            = 1.0f / sqrt(r2);
                    sig             = a.x + psA[tj].sig;
                    sig2            = invR * sig;
                    sig2           *= sig2;
                    sig6            = sig2 * sig2 * sig2;
                    eps             = a.y * psA[tj].eps;
                    dEdR            = eps * (12.0f * sig6 - 6.0f) * sig6;
#ifdef USE_CUTOFF
                    dEdR           += apos.w * psA[tj].q * (invR - 2.0f * cSim.reactionFieldK * r2);
#else
                    dEdR           += apos.w * psA[tj].q * invR;
#endif
                    dEdR           *= invR * invR;
#ifdef USE_CUTOFF
                    if (!(excl & 0x1) || r2 > cSim.nonbondedCutoffSqr)
#else
                    if (!(excl & 0x1))
#endif
                    {
                        dEdR = 0.0f;
                    }
                    dx             *= dEdR;
                    dy             *= dEdR;
                    dz             *= dEdR;
                    af.x           -= dx;
                    af.y           -= dy;
                    af.z           -= dz;
                    psA[tj].fx     += dx;
                    psA[tj].fy     += dy;
                    psA[tj].fz     += dz;
                    excl          >>= 1;
                    tj              = (tj + 1) & (GRID - 1);
                }
            }

            // Write results
            float4 of;
#ifdef USE_OUTPUT_BUFFER_PER_WARP
            int offset                          = x + tgx + warp*cSim.stride;
            of                                  = cSim.pForce4a[offset];
            of.x                               += af.x;
            of.y                               += af.y;
            of.z                               += af.z;
            cSim.pForce4a[offset]               = of;
            offset                              = y + tgx + warp*cSim.stride;
            of                                  = cSim.pForce4a[offset];
            of.x                               += sA[threadIdx.x].fx;
            of.y                               += sA[threadIdx.x].fy;
            of.z                               += sA[threadIdx.x].fz;
            cSim.pForce4a[offset]               = of;
#else
            of.x                                = af.x;
            of.y                                = af.y;
            of.z                                = af.z;
            of.w                                = 0.0f;
            int offset                          = x + tgx + (y >> GRIDBITS) * cSim.stride;
            cSim.pForce4a[offset]               = of;
            of.x                                = sA[threadIdx.x].fx;
            of.y                                = sA[threadIdx.x].fy;
            of.z                                = sA[threadIdx.x].fz;
            offset                              = y + tgx + (x >> GRIDBITS) * cSim.stride;
            cSim.pForce4a[offset]               = of;
#endif
            lasty = y;
        }

        pos++;
    }
}
