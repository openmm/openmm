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
 * This file contains the kernel for evalauating nonbonded forces and the first stage of GBSA.
 * It is included several times in kCalculateCDLJObcGbsaForces1.cu with different #defines to generate
 * different versions of the kernels.
 */

#define USE_SOFTCORE_LJ

#ifdef USE_SOFTCORE_LJ
#include "kSoftcoreLJ.h"
#endif

__global__ 
#if (__CUDA_ARCH__ >= 200)
__launch_bounds__(GF1XX_NONBOND_THREADS_PER_BLOCK, 1)
#elif (__CUDA_ARCH__ >= 120)
__launch_bounds__(GT2XX_NONBOND_THREADS_PER_BLOCK, 1)
#else
__launch_bounds__(G8X_NONBOND_THREADS_PER_BLOCK, 1)
#endif
#ifdef DEBUG
void METHOD_NAME(kCalculateCDLJObcGbsaSoftcore, Forces1_kernel)(unsigned int* workUnit, float4* pdE1, float4* pdE2 )
#else
void METHOD_NAME(kCalculateCDLJObcGbsaSoftcore, Forces1_kernel)(unsigned int* workUnit )
#endif
{
    extern __shared__ Atom sA[];
    unsigned int totalWarps        = gridDim.x*blockDim.x/GRID;
    unsigned int warp              = (blockIdx.x*blockDim.x+threadIdx.x)/GRID;

    unsigned int numWorkUnits      = cSim.pInteractionCount[0];
    unsigned int pos               = warp*numWorkUnits/totalWarps;
    unsigned int end               = (warp+1)*numWorkUnits/totalWarps;
    float CDLJObcGbsa_energy;
    float energy                   = 0.0f;
#ifdef USE_CUTOFF
    float* tempBuffer              = (float*) &sA[blockDim.x];
#endif

    unsigned int lasty             = -0xFFFFFFFF;
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
        float4 a                            = feSimDev.pSigEps4[i];
        float  softCoreLJLambda             = a.z;
        float br                            = cSim.pBornRadii[i];
        unsigned int tbx                    = threadIdx.x - tgx;
        unsigned int tj                     = tgx;
        Atom* psA                           = &sA[tbx];

        float4 af;
        af.x                                = 0.0f;
        af.y                                = 0.0f;
        af.z                                = 0.0f;
        af.w                                = 0.0f;

        if (x == y) // Handle diagonals uniquely at 50% efficiency
        {
            // Read fixed atom data into registers and GRF

            sA[threadIdx.x].x                    = apos.x;
            sA[threadIdx.x].y                    = apos.y;
            sA[threadIdx.x].z                    = apos.z;

            sA[threadIdx.x].q                    = a.w;
            sA[threadIdx.x].sig                  = a.x;
            sA[threadIdx.x].eps                  = a.y;
            sA[threadIdx.x].br                   = br;
            sA[threadIdx.x].softCoreLJLambda     = softCoreLJLambda;

            float q2                             = cSim.preFactor*a.w;
            a.w                                 *= cSim.epsfac;

            if (!bExclusionFlag)
            {
                for (unsigned int j = 0; j < GRID; j++)
                {
                    float dx                = psA[j].x - apos.x;
                    float dy                = psA[j].y - apos.y;
                    float dz                = psA[j].z - apos.z;
#ifdef USE_PERIODIC
                    dx                     -= floor(dx/cSim.periodicBoxSizeX+0.5f)*cSim.periodicBoxSizeX;
                    dy                     -= floor(dy/cSim.periodicBoxSizeY+0.5f)*cSim.periodicBoxSizeY;
                    dz                     -= floor(dz/cSim.periodicBoxSizeZ+0.5f)*cSim.periodicBoxSizeZ;
#endif
                    float r2                = dx * dx + dy * dy + dz * dz;
                    float invR              = 1.0f / sqrt(r2);
                    float sig               = a.x + psA[j].sig;
                    float eps               = a.y * psA[j].eps;
#ifdef USE_SOFTCORE_LJ
                    float dEdR              = getSoftCoreLJ( r2, sig, eps, softCoreLJLambda, psA[j].softCoreLJLambda, &CDLJObcGbsa_energy );
#else

                    // CDLJ part
                    float sig2              = invR * sig;
                    sig2                   *= sig2;
                    float sig6              = sig2 * sig2 * sig2;
                    float dEdR              = eps * (12.0f * sig6 - 6.0f) * sig6;
		              CDLJObcGbsa_energy      = eps * (sig6 - 1.0f) * sig6;
#endif
#ifdef USE_CUTOFF
                    dEdR                   += a.w * psA[j].q * (invR - 2.0f * feSimDev.reactionFieldK * r2);
                    CDLJObcGbsa_energy     += a.w * psA[j].q * (invR + feSimDev.reactionFieldK * r2 - feSimDev.reactionFieldC);
#else

                    float factorX           = a.w * psA[j].q * invR;
                    dEdR                   += factorX;
                    CDLJObcGbsa_energy     += factorX;

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
                    dEdR                   += Gpol * (1.0f - 0.25f * expTerm);
		              CDLJObcGbsa_energy     += (q2 * psA[j].q) / denominator;
#ifdef USE_CUTOFF
                    if ( i >= cSim.atoms || (x+j) >= cSim.atoms || r2 > cSim.nonbondedCutoffSqr)
#else
                    if ( i >= cSim.atoms || (x+j) >= cSim.atoms)
#endif
                    {
                        dEdR                = 0.0f;
			               CDLJObcGbsa_energy  = 0.0f;
                        dGpol_dalpha2_ij    = 0.0f;
                    }
                    af.w                   += dGpol_dalpha2_ij * psA[j].br;
                    energy                 += 0.5f*CDLJObcGbsa_energy;

                    // Add Forces

                    dx                     *= dEdR;
                    dy                     *= dEdR;
                    dz                     *= dEdR;

                    af.x                   -= dx;
                    af.y                   -= dy;
                    af.z                   -= dz;
                }

            } else {

                unsigned int xi   = x>>GRIDBITS;
                unsigned int cell = xi+xi*cSim.paddedNumberOfAtoms/GRID-xi*(xi+1)/2;
                unsigned int excl = cSim.pExclusion[cSim.pExclusionIndex[cell]+tgx];

                for (unsigned int j = 0; j < GRID; j++)
                {
                    float dx                = psA[j].x - apos.x;
                    float dy                = psA[j].y - apos.y;
                    float dz                = psA[j].z - apos.z;
#ifdef USE_PERIODIC
                    dx                     -= floor(dx/cSim.periodicBoxSizeX+0.5f)*cSim.periodicBoxSizeX;
                    dy                     -= floor(dy/cSim.periodicBoxSizeY+0.5f)*cSim.periodicBoxSizeY;
                    dz                     -= floor(dz/cSim.periodicBoxSizeZ+0.5f)*cSim.periodicBoxSizeZ;
#endif
                    float r2                = dx * dx + dy * dy + dz * dz;
                    float invR              = 1.0f / sqrt(r2);
                    float sig               = a.x + psA[j].sig;
                    float eps               = a.y * psA[j].eps;

#ifdef USE_SOFTCORE_LJ
                    float dEdR              = getSoftCoreLJ( r2, sig, eps, softCoreLJLambda, psA[j].softCoreLJLambda, &CDLJObcGbsa_energy );
                    //float dEdR              = getSoftCoreLJMod( (invR*sig), eps, softCoreLJLambda, psA[j].softCoreLJLambda, &CDLJObcGbsa_energy );
#else

                    // CDLJ part

                    float sig2              = invR * sig;
                    sig2                   *= sig2;
                    float sig6              = sig2 * sig2 * sig2;
                    float dEdR              = eps * (12.0f * sig6 - 6.0f) * sig6;
		              CDLJObcGbsa_energy      = eps * (sig6 - 1.0f) * sig6;
#endif

#ifdef USE_CUTOFF
                    dEdR                   += a.w * psA[j].q * (invR - 2.0f * feSimDev.reactionFieldK * r2);
                    CDLJObcGbsa_energy     += a.w * psA[j].q * (invR + feSimDev.reactionFieldK * r2 - feSimDev.reactionFieldC);
#else

                    float factorX           = a.w * psA[j].q * invR;
                    dEdR                   += factorX;
                    CDLJObcGbsa_energy     += factorX;

#endif
                    dEdR                   *= invR * invR;

                    if (!(excl & 0x1))
                    {
                        dEdR                = 0.0f;
                        CDLJObcGbsa_energy  = 0.0f;
                    }

                    // ObcGbsaForce1 part

                    float alpha2_ij         = br * psA[j].br;
                    float D_ij              = r2 / (4.0f * alpha2_ij);
                    float expTerm           = exp(-D_ij);
                    float denominator2      = r2 + alpha2_ij * expTerm;
                    float denominator       = sqrt(denominator2);
                    float Gpol              = (q2 * psA[j].q) / (denominator * denominator2);
                    float dGpol_dalpha2_ij  = -0.5f * Gpol * expTerm * (1.0f + D_ij);
                    dEdR                   += Gpol * (1.0f - 0.25f * expTerm);
                    CDLJObcGbsa_energy     += (q2 * psA[j].q) / denominator;

#if defined USE_CUTOFF
                    if (i >= cSim.atoms || x+j >= cSim.atoms || r2 > cSim.nonbondedCutoffSqr)
#else
                    if (i >= cSim.atoms || x+j >= cSim.atoms )
#endif
                    {
                        dEdR               = 0.0f;
		                  CDLJObcGbsa_energy = 0.0f;
                        dGpol_dalpha2_ij   = 0.0f;
                    }

                    af.w                  += dGpol_dalpha2_ij * psA[j].br;
                    energy                += 0.5f*CDLJObcGbsa_energy;
                     
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
            unsigned int offset         = x + tgx + warp*cSim.stride;
#else
            unsigned int offset         = x + tgx + (x >> GRIDBITS) * cSim.stride;
#endif

            float4 of                   = cSim.pForce4[offset];
            float  bf                   = cSim.pBornForce[offset];
            of.x                       += af.x;
            of.y                       += af.y;
            of.z                       += af.z;
            bf                         += af.w;
            cSim.pForce4[offset]        = of;
            cSim.pBornForce[offset]     = bf;

        } else { 

            // Read fixed atom data into registers and GRF

            if (lasty != y)
            {
                unsigned int j                       = y + tgx;
                float4 temp                          = cSim.pPosq[j];
                float4 temp1                         = feSimDev.pSigEps4[j];
                sA[threadIdx.x].br                   = cSim.pBornRadii[j];
                sA[threadIdx.x].x                    = temp.x;
                sA[threadIdx.x].y                    = temp.y;
                sA[threadIdx.x].z                    = temp.z;
                sA[threadIdx.x].q                    = temp1.w;
                sA[threadIdx.x].sig                  = temp1.x;
                sA[threadIdx.x].eps                  = temp1.y;
                sA[threadIdx.x].softCoreLJLambda     = temp1.z;
            }

            sA[threadIdx.x].fx          = 0.0f;
            sA[threadIdx.x].fy          = 0.0f;
            sA[threadIdx.x].fz          = 0.0f;
            sA[threadIdx.x].fb          = 0.0f;

            float q2                    = a.w * cSim.preFactor;
            a.w                        *= cSim.epsfac;
            if (!bExclusionFlag)
            {
#ifdef USE_CUTOFF
                unsigned int flags = cSim.pInteractionFlag[pos];
                if (flags == 0)
                {
                    // No interactions in this block.
                }
                //else if (flags == 0xFFFFFFFF)
                else if (flags)
#endif
                {
                    // Compute all interactions within this block.

                    for (unsigned int j = 0; j < GRID; j++)
                    {
                        float dx                = psA[tj].x - apos.x;
                        float dy                = psA[tj].y - apos.y;
                        float dz                = psA[tj].z - apos.z;
#ifdef USE_PERIODIC
                        dx                     -= floor(dx/cSim.periodicBoxSizeX+0.5f)*cSim.periodicBoxSizeX;
                        dy                     -= floor(dy/cSim.periodicBoxSizeY+0.5f)*cSim.periodicBoxSizeY;
                        dz                     -= floor(dz/cSim.periodicBoxSizeZ+0.5f)*cSim.periodicBoxSizeZ;
#endif
                        float r2                = dx * dx + dy * dy + dz * dz;
                        float invR              = 1.0f / sqrt(r2);

                        float sig               = a.x + psA[tj].sig;
                        float eps               = a.y * psA[tj].eps;
#ifdef USE_SOFTCORE_LJ
                        float dEdR              = getSoftCoreLJ( r2, sig, eps, softCoreLJLambda, psA[tj].softCoreLJLambda, &CDLJObcGbsa_energy );
#else
                        // CDLJ part
                        float sig2              = invR * sig;
                        sig2                   *= sig2;
                        float sig6              = sig2 * sig2 * sig2;
                        float dEdR              = eps * (12.0f * sig6 - 6.0f) * sig6;
                        CDLJObcGbsa_energy      = eps * (sig6 - 1.0f) * sig6;
#endif
#ifdef USE_CUTOFF
                        dEdR                   += a.w * psA[tj].q * (invR - 2.0f * feSimDev.reactionFieldK * r2);
                        CDLJObcGbsa_energy     += a.w * psA[tj].q * (invR + feSimDev.reactionFieldK * r2 - feSimDev.reactionFieldC);
#else

                        float factorX           = a.w * psA[tj].q * invR;
                        dEdR                   += factorX;
                        CDLJObcGbsa_energy     += factorX;
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
                        dEdR                   += Gpol * (1.0f - 0.25f * expTerm);
                        CDLJObcGbsa_energy     += (q2 * psA[tj].q) / denominator;
#ifdef USE_CUTOFF
                        if ( i >= cSim.atoms || (y+tj) >= cSim.atoms || r2 > cSim.nonbondedCutoffSqr)
#else
                        if ( i >= cSim.atoms || (y+tj) >= cSim.atoms)
#endif
                        {
                            dEdR               = 0.0f;
       			             CDLJObcGbsa_energy = 0.0f;
                            dGpol_dalpha2_ij   = 0.0f;
                        }
                        psA[tj].fb             += dGpol_dalpha2_ij * br;
                        af.w                   += dGpol_dalpha2_ij * psA[tj].br;
                        energy                 += CDLJObcGbsa_energy;

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
                else {

                    // Compute only a subset of the interactions in this block.

                    for (unsigned int j = 0; j < GRID; j++)
                    {
                        if ((flags&(1<<j)) != 0)
                        {
                            float dx                = psA[j].x - apos.x;
                            float dy                = psA[j].y - apos.y;
                            float dz                = psA[j].z - apos.z;
#ifdef USE_PERIODIC
                            dx                     -= floor(dx/cSim.periodicBoxSizeX+0.5f)*cSim.periodicBoxSizeX;
                            dy                     -= floor(dy/cSim.periodicBoxSizeY+0.5f)*cSim.periodicBoxSizeY;
                            dz                     -= floor(dz/cSim.periodicBoxSizeZ+0.5f)*cSim.periodicBoxSizeZ;
#endif
                            float r2                = dx * dx + dy * dy + dz * dz;
                            float invR              = 1.0f / sqrt(r2);
                            float sig               = a.x + psA[j].sig;
                            float eps               = a.y * psA[j].eps;
#ifdef USE_SOFTCORE_LJ
                            float dEdR              = getSoftCoreLJ( r2, sig, eps, softCoreLJLambda, psA[j].softCoreLJLambda, &CDLJObcGbsa_energy );
#else

                            // CDLJ part
                            float sig2              = invR * sig;
                            sig2                   *= sig2;
                            float sig6              = sig2 * sig2 * sig2;
                            float dEdR              = eps * (12.0f * sig6 - 6.0f) * sig6;
                            CDLJObcGbsa_energy      = eps * (sig6 - 1.0f) * sig6;
#endif
#ifdef USE_CUTOFF
                            dEdR                   += a.w * psA[j].q * (invR - 2.0f * feSimDev.reactionFieldK * r2);
                            CDLJObcGbsa_energy     += a.w * psA[j].q * (invR + feSimDev.reactionFieldK * r2 - feSimDev.reactionFieldC);
#else

                            float factorX           = a.w * psA[j].q * invR;
                            dEdR                   += factorX;
                            CDLJObcGbsa_energy     += factorX;
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
                            dEdR                   += Gpol * (1.0f - 0.25f * expTerm);
                            CDLJObcGbsa_energy     += (q2 * psA[j].q) / denominator;

#ifdef USE_CUTOFF
                            if ( i >= cSim.atoms || (y+j) >= cSim.atoms || r2 > cSim.nonbondedCutoffSqr)
#else
                            if ( i >= cSim.atoms || (y+j) >= cSim.atoms)
#endif
                            {
                                dEdR                = 0.0f;
				                    CDLJObcGbsa_energy  = 0.0f;
                                dGpol_dalpha2_ij    = 0.0f;
                            }
                            af.w                   += dGpol_dalpha2_ij * psA[j].br;

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

                            energy                 += CDLJObcGbsa_energy;

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
            } else {

                unsigned int xi   = x>>GRIDBITS;
                unsigned int yi   = y>>GRIDBITS;
                unsigned int cell = xi+yi*cSim.paddedNumberOfAtoms/GRID-yi*(yi+1)/2;
                unsigned int excl = cSim.pExclusion[cSim.pExclusionIndex[cell]+tgx];
                excl              = (excl >> tgx) | (excl << (GRID - tgx));
                for (unsigned int j = 0; j < GRID; j++)
                {
                    float dx                = psA[tj].x - apos.x;
                    float dy                = psA[tj].y - apos.y;
                    float dz                = psA[tj].z - apos.z;
#ifdef USE_PERIODIC
                    dx                     -= floor(dx/cSim.periodicBoxSizeX+0.5f)*cSim.periodicBoxSizeX;
                    dy                     -= floor(dy/cSim.periodicBoxSizeY+0.5f)*cSim.periodicBoxSizeY;
                    dz                     -= floor(dz/cSim.periodicBoxSizeZ+0.5f)*cSim.periodicBoxSizeZ;
#endif
                    float r2                = dx * dx + dy * dy + dz * dz;
                    float invR              = 1.0f / sqrt(r2);
                    float sig               = a.x + psA[tj].sig;
                    float eps               = a.y * psA[tj].eps;
#ifdef USE_SOFTCORE_LJ
                    float dEdR              = getSoftCoreLJ( r2, sig, eps, softCoreLJLambda, psA[tj].softCoreLJLambda, &CDLJObcGbsa_energy );
#else

                    // CDLJ part
                    float sig2              = invR * sig;
                    sig2                   *= sig2;
                    float sig6              = sig2 * sig2 * sig2;
                    float dEdR              = eps * (12.0f * sig6 - 6.0f) * sig6;
		              CDLJObcGbsa_energy      = eps * (sig6 - 1.0f) * sig6;
#endif
#ifdef USE_CUTOFF
                    dEdR                   += a.w * psA[tj].q * (invR - 2.0f * feSimDev.reactionFieldK * r2);
                    CDLJObcGbsa_energy     += a.w * psA[tj].q * (invR + feSimDev.reactionFieldK * r2 - feSimDev.reactionFieldC);
#else
                    float factorX           = a.w * psA[tj].q * invR;
                    dEdR                   += factorX;
                    CDLJObcGbsa_energy     += factorX;
#endif
                    dEdR                   *= invR * invR;
                    if (!(excl & 0x1))
                    {
                        dEdR               = 0.0f;
			               CDLJObcGbsa_energy = 0.0f;
                    }

                    // ObcGbsaForce1 part

                    float alpha2_ij         = br * psA[tj].br;
                    float D_ij              = r2 / (4.0f * alpha2_ij);
                    float expTerm           = exp(-D_ij);
                    float denominator2      = r2 + alpha2_ij * expTerm;
                    float denominator       = sqrt(denominator2);
                    float Gpol              = (q2 * psA[tj].q) / (denominator * denominator2);
                    float dGpol_dalpha2_ij  = -0.5f * Gpol * expTerm * (1.0f + D_ij);
                    dEdR                   += Gpol * (1.0f - 0.25f * expTerm);
		              CDLJObcGbsa_energy     += (q2 * psA[tj].q) / denominator;
#if defined USE_CUTOFF
                    if (i >= cSim.atoms || y+tj >= cSim.atoms || r2 > cSim.nonbondedCutoffSqr)
#else
                    if (i >= cSim.atoms || y+tj >= cSim.atoms)
#endif
                    {
                        dEdR               = 0.0f;
			               CDLJObcGbsa_energy = 0.0f;
                        dGpol_dalpha2_ij   = 0.0f;
                    }

                    af.w                   += dGpol_dalpha2_ij * psA[tj].br;
                    psA[tj].fb             += dGpol_dalpha2_ij * br;
                    energy                 += CDLJObcGbsa_energy;

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
            unsigned int offset         = x + tgx + warp*cSim.stride;
#else
            unsigned int offset         = x + tgx + (y >> GRIDBITS) * cSim.stride;
#endif
            float4 of                   = cSim.pForce4[offset];
            float  bf                   = cSim.pBornForce[offset];
            of.x                       += af.x;
            of.y                       += af.y;
            of.z                       += af.z;
            bf                         += af.w;
            cSim.pForce4[offset]        = of;
            cSim.pBornForce[offset]     = bf;

#ifdef USE_OUTPUT_BUFFER_PER_WARP
            offset                      = y + tgx + warp*cSim.stride;
#else
            offset                      = y + tgx + (x >> GRIDBITS) * cSim.stride;
#endif
            of                          = cSim.pForce4[offset];
            bf                          = cSim.pBornForce[offset];
            of.x                       += sA[threadIdx.x].fx;
            of.y                       += sA[threadIdx.x].fy;
            of.z                       += sA[threadIdx.x].fz;
            bf                         += sA[threadIdx.x].fb;
            cSim.pForce4[offset]        = of;
            cSim.pBornForce[offset]     = bf;

            lasty                       = y;
        }
        pos++;
    }
    cSim.pEnergy[blockIdx.x*blockDim.x+threadIdx.x] += energy;
}
