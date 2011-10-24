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
 * This file contains the kernels for evaluating nonbonded softcore forces.  It is included
 * several times in kCalculateNonbondedSoftcore.cu with different #defines to generate
 * different versions of the kernels.
 */

__global__ 
#if (__CUDA_ARCH__ >= 200)
__launch_bounds__(GF1XX_NONBOND_THREADS_PER_BLOCK, 1)
#elif (__CUDA_ARCH__ >= 120)
__launch_bounds__(GT2XX_NONBOND_THREADS_PER_BLOCK, 1)
#else
__launch_bounds__(G8X_NONBOND_THREADS_PER_BLOCK, 1)
#endif
void METHOD_NAME(kCalculateCDLJSoftcore, Forces_kernel)(unsigned int* workUnit )
{
    extern __shared__ Atom sA[];
    unsigned int totalWarps   = cSim.nonbond_blocks*cSim.nonbond_threads_per_block/GRID;
    unsigned int warp         = (blockIdx.x*blockDim.x+threadIdx.x)/GRID;
    unsigned int numWorkUnits = cSim.pInteractionCount[0];
    unsigned int pos          = warp*numWorkUnits/totalWarps;
    unsigned int end          = (warp+1)*numWorkUnits/totalWarps;
    float CDLJ_energy;
    float energy              = 0.0f;
#ifdef USE_CUTOFF
    float3* tempBuffer        = (float3*) &sA[cSim.nonbond_threads_per_block];
#endif

    unsigned int lasty = 0xFFFFFFFF;
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
        float eps;
        float dEdR;

        unsigned int tgx             = threadIdx.x & (GRID - 1);
        unsigned int tbx             = threadIdx.x - tgx;
        unsigned int tj              = tgx;

        Atom* psA                    = &sA[tbx];
        unsigned int i               = x + tgx;

        apos                         = cSim.pPosq[i];
        float4 a                     = feSimDev.pSigEps4[i];
        float softCoreLJLambda       = a.z;
        af.x                         = 0.0f;
        af.y                         = 0.0f;
        af.z                         = 0.0f;
        if (x == y) // Handle diagonals uniquely at 50% efficiency
        {
            // Read fixed atom data into registers and GRF
            sA[threadIdx.x].x                     = apos.x;
            sA[threadIdx.x].y                     = apos.y;
            sA[threadIdx.x].z                     = apos.z;
            sA[threadIdx.x].q                     = a.w;
            sA[threadIdx.x].sig                   = a.x;
            sA[threadIdx.x].eps                   = a.y;
            sA[threadIdx.x].softCoreLJLambda      = a.z;
            a.w                                  *= cSim.epsfac;
            if (!bExclusionFlag)
            {
                for (unsigned int j = 0; j < GRID; j++)
                {
                    dx              = psA[j].x - apos.x;
                    dy              = psA[j].y - apos.y;
                    dz              = psA[j].z - apos.z;
#ifdef USE_PERIODIC
                    dx             -= floor(dx/cSim.periodicBoxSizeX+0.5f)*cSim.periodicBoxSizeX;
                    dy             -= floor(dy/cSim.periodicBoxSizeY+0.5f)*cSim.periodicBoxSizeY;
                    dz             -= floor(dz/cSim.periodicBoxSizeZ+0.5f)*cSim.periodicBoxSizeZ;
#endif
                    r2              = dx * dx + dy * dy + dz * dz;
                    invR            = 1.0f / sqrt(r2);
                    sig             = a.x + psA[j].sig;
                    eps             = a.y * psA[j].eps;
#ifdef USE_SOFTCORE_LJ
                    dEdR            = getSoftCoreLJ( r2, sig, eps, softCoreLJLambda, psA[j].softCoreLJLambda, &CDLJ_energy );
#else
                    float sig2      = invR * sig;
                    sig2           *= sig2;
                    float sig6      = sig2 * sig2 * sig2;
                    dEdR            = eps * (12.0f * sig6 - 6.0f) * sig6;
	          	     CDLJ_energy     = eps * (sig6 - 1.0f) * sig6;
#endif

#ifdef USE_CUTOFF
                    dEdR           += a.w * psA[j].q * (invR - 2.0f * feSimDev.reactionFieldK * r2);
		              CDLJ_energy    += a.w * psA[j].q * (invR + feSimDev.reactionFieldK * r2 - feSimDev.reactionFieldC);
#else
                    dEdR           += a.w * psA[j].q * invR;
		              CDLJ_energy    += a.w * psA[j].q * invR;
#endif
                    dEdR           *= invR * invR;
#ifdef USE_CUTOFF
                    if (r2 > cSim.nonbondedCutoffSqr)
                    {
                        dEdR        = 0.0f;
                        CDLJ_energy = 0.0f;
                    }
#endif
		              energy         += 0.5f*CDLJ_energy;
                    dx             *= dEdR;
                    dy             *= dEdR;
                    dz             *= dEdR;
                    af.x           -= dx;
                    af.y           -= dy;
                    af.z           -= dz;
                }

            } else  {

                unsigned int xi   = x>>GRIDBITS;
                unsigned int cell = xi+xi*cSim.paddedNumberOfAtoms/GRID-xi*(xi+1)/2;
                unsigned int excl = cSim.pExclusion[cSim.pExclusionIndex[cell]+tgx];

                for (unsigned int j = 0; j < GRID; j++)
                {
                    dx              = psA[j].x - apos.x;
                    dy              = psA[j].y - apos.y;
                    dz              = psA[j].z - apos.z;
#ifdef USE_PERIODIC
                    dx             -= floor(dx/cSim.periodicBoxSizeX+0.5f)*cSim.periodicBoxSizeX;
                    dy             -= floor(dy/cSim.periodicBoxSizeY+0.5f)*cSim.periodicBoxSizeY;
                    dz             -= floor(dz/cSim.periodicBoxSizeZ+0.5f)*cSim.periodicBoxSizeZ;
#endif
                    r2              = dx * dx + dy * dy + dz * dz;
                    invR            = 1.0f / sqrt(r2);
                    sig             = a.x + psA[j].sig;
                    eps             = a.y * psA[j].eps;
#ifdef USE_SOFTCORE_LJ
                    dEdR            = getSoftCoreLJ( r2, sig, eps, softCoreLJLambda, psA[j].softCoreLJLambda, &CDLJ_energy );
#else
                    float sig2      = invR * sig;
                    sig2           *= sig2;
                    float sig6      = sig2 * sig2 * sig2;
                    dEdR            = eps * (12.0f * sig6 - 6.0f) * sig6;
                    CDLJ_energy     = eps * (sig6 - 1.0f) * sig6;
#endif

#ifdef USE_CUTOFF
                    dEdR           += a.w * psA[j].q * (invR - 2.0f * feSimDev.reactionFieldK * r2);
		              CDLJ_energy    += a.w * psA[j].q * (invR + feSimDev.reactionFieldK * r2 - feSimDev.reactionFieldC);
#else
                    dEdR           += a.w * psA[j].q * invR;
		              CDLJ_energy    += a.w * psA[j].q * invR;
#endif
                    dEdR           *= invR * invR;
#ifdef USE_CUTOFF
                    if (!(excl & 0x1) || r2 > cSim.nonbondedCutoffSqr)
#else
                    if (!(excl & 0x1))
#endif
                    {
                        dEdR = 0.0f;
		                  CDLJ_energy  = 0.0f;
                    }

                    energy         += 0.5f*CDLJ_energy;
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
            unsigned int offset                 = x + tgx + warp*cSim.stride;
            of                                  = cSim.pForce4[offset];
            of.x                               += af.x;
            of.y                               += af.y;
            of.z                               += af.z;
            cSim.pForce4[offset]                = of;
#else
            of.x                                = af.x;
            of.y                                = af.y;
            of.z                                = af.z;
            of.w                                = 0.0f;
            unsigned int offset                 = x + tgx + (x >> GRIDBITS) * cSim.stride;
            cSim.pForce4[offset]                = of;
#endif

        } else {

            // Read fixed atom data into registers and GRF
            if (lasty != y)
            {
                unsigned int j                   = y + tgx;
                float4 temp                      = cSim.pPosq[j];
                //float2 temp1                     = cSim.pAttr[j];
                float4 temp1                     = feSimDev.pSigEps4[j];
                //float  temp3                     = cSim.pSoftCoreLJLambda[j];
                //float  temp3                     = softCoreLJLambdaArray[j];
                float temp3                      = temp1.z;
                sA[threadIdx.x].x                = temp.x;
                sA[threadIdx.x].y                = temp.y;
                sA[threadIdx.x].z                = temp.z;
                sA[threadIdx.x].q                = temp1.w;
                sA[threadIdx.x].sig              = temp1.x;
                sA[threadIdx.x].eps              = temp1.y;
                sA[threadIdx.x].softCoreLJLambda = temp3;
            }
            sA[threadIdx.x].fx      = 0.0f;
            sA[threadIdx.x].fy      = 0.0f;
            sA[threadIdx.x].fz      = 0.0f;
            a.w                *= cSim.epsfac;
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
                        eps             = a.y * psA[tj].eps;
#ifdef USE_SOFTCORE_LJ
                        dEdR            = getSoftCoreLJ( r2, sig, eps, softCoreLJLambda, psA[tj].softCoreLJLambda, &CDLJ_energy );
#else
                        float sig2      = invR * sig;
                        sig2           *= sig2;
                        float sig6      = sig2 * sig2 * sig2;
                        dEdR            = eps * (12.0f * sig6 - 6.0f) * sig6;
			               CDLJ_energy     = eps * (sig6 - 1.0f) * sig6;
#endif
#ifdef USE_CUTOFF
                        dEdR           += a.w * psA[tj].q * (invR - 2.0f * feSimDev.reactionFieldK * r2);
                        CDLJ_energy    += a.w * psA[tj].q * (invR + feSimDev.reactionFieldK * r2 - feSimDev.reactionFieldC);
#else
                        dEdR           += a.w * psA[tj].q * invR;
                        CDLJ_energy    += a.w * psA[tj].q * invR;
#endif
                        dEdR           *= invR * invR;
#ifdef USE_CUTOFF
                        if (r2 > cSim.nonbondedCutoffSqr)
                        {
                            dEdR = 0.0f;
       			             CDLJ_energy = 0.0f;
                        }
#endif

			               energy         += CDLJ_energy;
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
#ifdef USE_CUTOFF
                else
                {
                    // Compute only a subset of the interactions in this block.

                    for (unsigned int j = 0; j < GRID; j++)
                    {
                        if ((flags&(1<<j)) != 0)
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
                            eps             = a.y * psA[j].eps;
#ifdef USE_SOFTCORE_LJ
                            dEdR            = getSoftCoreLJ( r2, sig, eps, softCoreLJLambda, psA[j].softCoreLJLambda, &CDLJ_energy );
#else
                            float sig2      = invR * sig;
                            sig2           *= sig2;
                            float sig6      = sig2 * sig2 * sig2;
                            dEdR            = eps * (12.0f * sig6 - 6.0f) * sig6;
                			    CDLJ_energy     = eps * (sig6 - 1.0f) * sig6;
#endif
#ifdef USE_CUTOFF
                            dEdR           += a.w * psA[j].q * (invR - 2.0f * feSimDev.reactionFieldK * r2);
                            CDLJ_energy    += a.w * psA[j].q * (invR + feSimDev.reactionFieldK * r2 - feSimDev.reactionFieldC);
#else
                            dEdR           += a.w * psA[j].q * invR;
                            CDLJ_energy    += a.w * psA[j].q * invR;
#endif
                            dEdR           *= invR * invR;
#ifdef USE_CUTOFF
                            if (r2 > cSim.nonbondedCutoffSqr)
                            {
                                dEdR = 0.0f;
				                    CDLJ_energy = 0.0f;
                            }
#endif
			                   energy         += CDLJ_energy;
                            dx             *= dEdR;
                            dy             *= dEdR;
                            dz             *= dEdR;
                            af.x           -= dx;
                            af.y           -= dy;
                            af.z           -= dz;
                            tempBuffer[threadIdx.x].x = dx;
                            tempBuffer[threadIdx.x].y = dy;
                            tempBuffer[threadIdx.x].z = dz;

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
            }
            else  // bExclusion
            {
                // Read fixed atom data into registers and GRF
                unsigned int xi   = x>>GRIDBITS;
                unsigned int yi   = y>>GRIDBITS;
                unsigned int cell          = xi+yi*cSim.paddedNumberOfAtoms/GRID-yi*(yi+1)/2;
                unsigned int excl = cSim.pExclusion[cSim.pExclusionIndex[cell]+tgx];
                excl              = (excl >> tgx) | (excl << (GRID - tgx));
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
                    eps             = a.y * psA[tj].eps;
#ifdef USE_SOFTCORE_LJ
                    dEdR            = getSoftCoreLJ( r2, sig, eps, softCoreLJLambda, psA[tj].softCoreLJLambda, &CDLJ_energy );
#else
                    float sig2      = invR * sig;
                    sig2           *= sig2;
                    float sig6      = sig2 * sig2 * sig2;
                    dEdR            = eps * (12.0f * sig6 - 6.0f) * sig6;
		              CDLJ_energy     = eps * (sig6 - 1.0f) * sig6;
#endif

#ifdef USE_CUTOFF
                    dEdR           += a.w * psA[tj].q * (invR - 2.0f * feSimDev.reactionFieldK * r2);
	                 CDLJ_energy    += a.w * psA[tj].q * (invR + feSimDev.reactionFieldK * r2 - feSimDev.reactionFieldC);
#else
                    dEdR           += a.w * psA[tj].q * invR;
                    CDLJ_energy    += a.w * psA[tj].q * invR;
#endif
                    dEdR           *= invR * invR;
#ifdef USE_CUTOFF
                    if (!(excl & 0x1) || r2 > cSim.nonbondedCutoffSqr)
#else
                    if (!(excl & 0x1))
#endif
                    {
                        dEdR = 0.0f;			
	                     CDLJ_energy  = 0.0f;
                    }

		              energy         += CDLJ_energy;
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
            unsigned int offset                          = x + tgx + warp*cSim.stride;
            of                                  = cSim.pForce4[offset];
            of.x                               += af.x;
            of.y                               += af.y;
            of.z                               += af.z;
            cSim.pForce4[offset]                = of;
            offset                              = y + tgx + warp*cSim.stride;
            of                                  = cSim.pForce4[offset];
            of.x                               += sA[threadIdx.x].fx;
            of.y                               += sA[threadIdx.x].fy;
            of.z                               += sA[threadIdx.x].fz;
            cSim.pForce4[offset]                = of;
#else
            of.x                                = af.x;
            of.y                                = af.y;
            of.z                                = af.z;
            of.w                                = 0.0f;
            unsigned int offset                          = x + tgx + (y >> GRIDBITS) * cSim.stride;
            cSim.pForce4[offset]                = of;
            of.x                                = sA[threadIdx.x].fx;
            of.y                                = sA[threadIdx.x].fy;
            of.z                                = sA[threadIdx.x].fz;
            offset                              = y + tgx + (x >> GRIDBITS) * cSim.stride;
            cSim.pForce4[offset]                = of;
#endif
            lasty = y;
        }

        pos++;
    }
    cSim.pEnergy[blockIdx.x*blockDim.x+threadIdx.x] += energy;
}
