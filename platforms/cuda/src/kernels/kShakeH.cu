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

#include <stdio.h>
#include <cuda.h>
#include <vector_functions.h>
#include <cstdlib>
#include <string>
#include <iostream>
using namespace std;

#include "gputypes.h"

struct Atom
{
    float3 rij1;
    float3 rij2;
    float3 rij3;
    float  M;
    float  d2;
    float  InvMassI;
    float  rij1sq;
    float  rij2sq;
    float  rij3sq;
};


static __constant__ cudaGmxSimulation cSim;

void SetShakeHSim(gpuContext gpu)
{
    cudaError_t status;
    status = cudaMemcpyToSymbol(cSim, &gpu->sim, sizeof(cudaGmxSimulation));
    RTERROR(status, "cudaMemcpyToSymbol: SetSim copy to cSim failed");
}

void GetShakeHSim(gpuContext gpu)
{
    cudaError_t status;
    status = cudaMemcpyFromSymbol(&gpu->sim, cSim, sizeof(cudaGmxSimulation));
    RTERROR(status, "cudaMemcpyFromSymbol: SetSim copy from cSim failed");
}

__global__ void kApplyFirstShake_kernel()
{
    __shared__ Atom sA[G8X_THREADS_PER_BLOCK];
    Atom* psA = &sA[threadIdx.x];
    unsigned int pos = threadIdx.x + blockIdx.x * blockDim.x;
    while (pos < cSim.ShakeConstraints)
    {
        int4 atomID         = cSim.pShakeID[pos];
        float4 params       = cSim.pShakeParameter[pos];
        float4 apos         = cSim.pOldPosq[atomID.x];
        float4 xpi          = cSim.pPosqP[atomID.x];
        float4 apos1        = cSim.pOldPosq[atomID.y];
        float4 xpj1         = cSim.pPosqP[atomID.y];
        float4 apos2        = {0.0f, 0.0f, 0.0f, 0.0f};
        float4 xpj2         = {0.0f, 0.0f, 0.0f, 0.0f};
        psA->InvMassI       = params.x;
        psA->M              = params.y;
        psA->d2             = params.z;
        float invMassJ      = params.w;
        if (atomID.z != -1)
        {
            apos2           = cSim.pOldPosq[atomID.z];
            xpj2            = cSim.pPosqP[atomID.z];
        }
        float4 apos3        = {0.0f, 0.0f, 0.0f, 0.0f};
        float4 xpj3         = {0.0f, 0.0f, 0.0f, 0.0f};
        if (atomID.w != -1)
        {
            apos3           = cSim.pOldPosq[atomID.w];
            xpj3            = cSim.pPosqP[atomID.w];
        }

        float3 xi, xj1, xj2, xj3;
        xi.x                = apos.x;
        xi.y                = apos.y;
        xi.z                = apos.z;
        xj1.x               = apos1.x;
        xj1.y               = apos1.y;
        xj1.z               = apos1.z;
        xj2.x               = apos2.x;
        xj2.y               = apos2.y;
        xj2.z               = apos2.z;
        xj3.x               = apos3.x;
        xj3.y               = apos3.y;
        xj3.z               = apos3.z;
        psA->rij1.x         = xi.x - xj1.x;
        psA->rij1.y         = xi.y - xj1.y;
        psA->rij1.z         = xi.z - xj1.z;
        psA->rij2.x         = xi.x - xj2.x;
        psA->rij2.y         = xi.y - xj2.y;
        psA->rij2.z         = xi.z - xj2.z;
        psA->rij3.x         = xi.x - xj3.x;
        psA->rij3.y         = xi.y - xj3.y;
        psA->rij3.z         = xi.z - xj3.z;
        psA->rij1sq         = psA->rij1.x * psA->rij1.x + psA->rij1.y * psA->rij1.y + psA->rij1.z * psA->rij1.z;
        psA->rij2sq         = psA->rij2.x * psA->rij2.x + psA->rij2.y * psA->rij2.y + psA->rij2.z * psA->rij2.z;
        psA->rij3sq         = psA->rij3.x * psA->rij3.x + psA->rij3.y * psA->rij3.y + psA->rij3.z * psA->rij3.z;
        float ld1           = psA->d2 - psA->rij1sq;
        float ld2           = psA->d2 - psA->rij2sq;
        float ld3           = psA->d2 - psA->rij3sq;


        bool converged = false;
        int iteration = 0;
        while (iteration < 15 && !converged)
        {
            converged = true;
            float3 rpij;
            rpij.x          = xpi.x - xpj1.x;
            rpij.y          = xpi.y - xpj1.y;
            rpij.z          = xpi.z - xpj1.z;
		    float rpsqij    = rpij.x * rpij.x + rpij.y * rpij.y + rpij.z * rpij.z;
		    float rrpr      = psA->rij1.x * rpij.x + psA->rij1.y * rpij.y + psA->rij1.z * rpij.z;
		    float diff      = fabs(ld1 - 2.0f * rrpr - rpsqij) / (psA->d2 * cSim.shakeTolerance);
            if (diff >= 1.0f)
            {
                float acor  = (ld1 - 2.0f * rrpr - rpsqij) * psA->M / (rrpr + psA->rij1sq);
                float3 dr;
                dr.x    = psA->rij1.x * acor;
                dr.y    = psA->rij1.y * acor;
                dr.z    = psA->rij1.z * acor;
		        xpi.x  += dr.x * psA->InvMassI;
		        xpi.y  += dr.y * psA->InvMassI;
		        xpi.z  += dr.z * psA->InvMassI;
		        xpj1.x -= dr.x * invMassJ;
		        xpj1.y -= dr.y * invMassJ;
		        xpj1.z -= dr.z * invMassJ;
                converged = false;
            }

            if (atomID.z != -1)
            {
                rpij.x          = xpi.x - xpj2.x;
                rpij.y          = xpi.y - xpj2.y;
                rpij.z          = xpi.z - xpj2.z;
		        rpsqij          = rpij.x * rpij.x + rpij.y * rpij.y + rpij.z * rpij.z;
		        rrpr            = psA->rij2.x * rpij.x + psA->rij2.y * rpij.y + psA->rij2.z * rpij.z;
		        diff            = fabs(ld2 - 2.0f * rrpr - rpsqij) / (psA->d2 * cSim.shakeTolerance);
                if (diff >= 1.0f)
                {
                    float acor  = (ld2 - 2.0f * rrpr - rpsqij) * psA->M / (rrpr + psA->rij2sq);
                    float3 dr;
                    dr.x    = psA->rij2.x * acor;
                    dr.y    = psA->rij2.y * acor;
                    dr.z    = psA->rij2.z * acor;
		            xpi.x  += dr.x * psA->InvMassI;
		            xpi.y  += dr.y * psA->InvMassI;
		            xpi.z  += dr.z * psA->InvMassI;
		            xpj2.x -= dr.x * invMassJ;
		            xpj2.y -= dr.y * invMassJ;
		            xpj2.z -= dr.z * invMassJ;
                    converged = false;
                }
            }

            if (atomID.w != -1)
            {
                rpij.x          = xpi.x - xpj3.x;
                rpij.y          = xpi.y - xpj3.y;
                rpij.z          = xpi.z - xpj3.z;
		        rpsqij          = rpij.x * rpij.x + rpij.y * rpij.y + rpij.z * rpij.z;
		        rrpr            = psA->rij3.x * rpij.x + psA->rij3.y * rpij.y + psA->rij3.z * rpij.z;
		        diff            = fabs(ld3 - 2.0f * rrpr - rpsqij) / (psA->d2 * cSim.shakeTolerance);
                if (diff >= 1.0f)
                {
                    float acor  = (ld3 - 2.0f * rrpr - rpsqij) * psA->M / (rrpr + psA->rij3sq);
                    float3 dr;
                    dr.x    = psA->rij3.x * acor;
                    dr.y    = psA->rij3.y * acor;
                    dr.z    = psA->rij3.z * acor;
		            xpi.x  += dr.x * psA->InvMassI;
		            xpi.y  += dr.y * psA->InvMassI;
		            xpi.z  += dr.z * psA->InvMassI;
		            xpj3.x -= dr.x * invMassJ;
		            xpj3.y -= dr.y * invMassJ;
		            xpj3.z -= dr.z * invMassJ;
                    converged = false;
                }
            }
            iteration++;
        }
        cSim.pPosqP[atomID.x] = xpi;
        cSim.pPosqP[atomID.y] = xpj1;
        if (atomID.z != -1)
            cSim.pPosqP[atomID.z] = xpj2;
        if (atomID.w != -1)
            cSim.pPosqP[atomID.w] = xpj3;

        pos += blockDim.x * gridDim.x;
    }
}

void kApplyFirstShake(gpuContext gpu)
{
//    printf("kApplyFirstShake\n");
    if (gpu->sim.ShakeConstraints > 0)
    {
        kApplyFirstShake_kernel<<<gpu->sim.blocks, gpu->sim.shake_threads_per_block>>>();
        LAUNCHERROR("kApplyFirstShake");
    }
}

__global__ void kApplySecondShake_kernel()
{
    __shared__ Atom sA[G8X_THREADS_PER_BLOCK];
    Atom* psA = &sA[threadIdx.x];
    unsigned int pos = threadIdx.x + blockIdx.x * blockDim.x;
    while (pos < cSim.ShakeConstraints)
    {
        int4 atomID         = cSim.pShakeID[pos];
        float4 params       = cSim.pShakeParameter[pos];
        float4 apos         = cSim.pOldPosq[atomID.x];
        float4 xpi          = cSim.pPosq[atomID.x];
        float4 apos1        = cSim.pOldPosq[atomID.y];
        float4 xpj1         = cSim.pPosq[atomID.y];
        float4 apos2        = {0.0f, 0.0f, 0.0f, 0.0f};
        float4 xpj2         = {0.0f, 0.0f, 0.0f, 0.0f};
        psA->InvMassI       = params.x;
        psA->M              = params.y;
        psA->d2             = params.z;
        float invMassJ      = params.w;
        if (atomID.z != -1)
        {
            apos2           = cSim.pOldPosq[atomID.z];
            xpj2            = cSim.pPosq[atomID.z];
        }
        float4 apos3        = {0.0f, 0.0f, 0.0f, 0.0f};
        float4 xpj3         = {0.0f, 0.0f, 0.0f, 0.0f};
        if (atomID.w != -1)
        {
            apos3           = cSim.pOldPosq[atomID.w];
            xpj3            = cSim.pPosq[atomID.w];
        }

        float3 xi, xj1, xj2, xj3;
        xi.x                = apos.x;
        xi.y                = apos.y;
        xi.z                = apos.z;
        xj1.x               = apos1.x;
        xj1.y               = apos1.y;
        xj1.z               = apos1.z;
        xj2.x               = apos2.x;
        xj2.y               = apos2.y;
        xj2.z               = apos2.z;
        xj3.x               = apos3.x;
        xj3.y               = apos3.y;
        xj3.z               = apos3.z;
        psA->rij1.x         = xi.x - xj1.x;
        psA->rij1.y         = xi.y - xj1.y;
        psA->rij1.z         = xi.z - xj1.z;
        psA->rij2.x         = xi.x - xj2.x;
        psA->rij2.y         = xi.y - xj2.y;
        psA->rij2.z         = xi.z - xj2.z;
        psA->rij3.x         = xi.x - xj3.x;
        psA->rij3.y         = xi.y - xj3.y;
        psA->rij3.z         = xi.z - xj3.z;
        psA->rij1sq         = psA->rij1.x * psA->rij1.x + psA->rij1.y * psA->rij1.y + psA->rij1.z * psA->rij1.z;
        psA->rij2sq         = psA->rij2.x * psA->rij2.x + psA->rij2.y * psA->rij2.y + psA->rij2.z * psA->rij2.z;
        psA->rij3sq         = psA->rij3.x * psA->rij3.x + psA->rij3.y * psA->rij3.y + psA->rij3.z * psA->rij3.z;
        float ld1           = psA->d2 - psA->rij1sq;
        float ld2           = psA->d2 - psA->rij2sq;
        float ld3           = psA->d2 - psA->rij3sq;


        bool converged = false;
        int iteration = 0;
        while (iteration < 15 && !converged)
        {
            converged = true;
            float3 rpij;
            rpij.x          = xpi.x - xpj1.x;
            rpij.y          = xpi.y - xpj1.y;
            rpij.z          = xpi.z - xpj1.z;
		    float rpsqij    = rpij.x * rpij.x + rpij.y * rpij.y + rpij.z * rpij.z;
		    float rrpr      = psA->rij1.x * rpij.x + psA->rij1.y * rpij.y + psA->rij1.z * rpij.z;
		    float diff      = fabs(ld1 - 2.0f * rrpr - rpsqij) / (psA->d2 * cSim.shakeTolerance );
            if (diff >= 1.0f)
            {
                float acor  = (ld1 - 2.0f * rrpr - rpsqij) * psA->M / (rrpr + psA->rij1sq);
                float3 dr;
                dr.x    = psA->rij1.x * acor;
                dr.y    = psA->rij1.y * acor;
                dr.z    = psA->rij1.z * acor;
		        xpi.x  += dr.x * psA->InvMassI;
		        xpi.y  += dr.y * psA->InvMassI;
		        xpi.z  += dr.z * psA->InvMassI;
		        xpj1.x -= dr.x * invMassJ;
		        xpj1.y -= dr.y * invMassJ;
		        xpj1.z -= dr.z * invMassJ;
                converged = false;
            }

            if (atomID.z != -1)
            {
                rpij.x          = xpi.x - xpj2.x;
                rpij.y          = xpi.y - xpj2.y;
                rpij.z          = xpi.z - xpj2.z;
		        rpsqij          = rpij.x * rpij.x + rpij.y * rpij.y + rpij.z * rpij.z;
		        rrpr            = psA->rij2.x * rpij.x + psA->rij2.y * rpij.y + psA->rij2.z * rpij.z;
		        diff            = fabs(ld2 - 2.0f * rrpr - rpsqij) / (psA->d2 * cSim.shakeTolerance );
                if (diff >= 1.0f)
                {
                    float acor  = (ld2 - 2.0f * rrpr - rpsqij) * psA->M / (rrpr + psA->rij2sq);
                    float3 dr;
                    dr.x    = psA->rij2.x * acor;
                    dr.y    = psA->rij2.y * acor;
                    dr.z    = psA->rij2.z * acor;
		            xpi.x  += dr.x * psA->InvMassI;
		            xpi.y  += dr.y * psA->InvMassI;
		            xpi.z  += dr.z * psA->InvMassI;
		            xpj2.x -= dr.x * invMassJ;
		            xpj2.y -= dr.y * invMassJ;
		            xpj2.z -= dr.z * invMassJ;
                    converged = false;
                }
            }

            if (atomID.w != -1)
            {
                rpij.x          = xpi.x - xpj3.x;
                rpij.y          = xpi.y - xpj3.y;
                rpij.z          = xpi.z - xpj3.z;
		        rpsqij          = rpij.x * rpij.x + rpij.y * rpij.y + rpij.z * rpij.z;
		        rrpr            = psA->rij3.x * rpij.x + psA->rij3.y * rpij.y + psA->rij3.z * rpij.z;
		        diff            = fabs(ld3 - 2.0f * rrpr - rpsqij) / (psA->d2 * cSim.shakeTolerance );
                if (diff >= 1.0f)
                {
                    float acor  = (ld3 - 2.0f * rrpr - rpsqij) * psA->M / (rrpr + psA->rij3sq);
                    float3 dr;
                    dr.x    = psA->rij3.x * acor;
                    dr.y    = psA->rij3.y * acor;
                    dr.z    = psA->rij3.z * acor;
		            xpi.x  += dr.x * psA->InvMassI;
		            xpi.y  += dr.y * psA->InvMassI;
		            xpi.z  += dr.z * psA->InvMassI;
		            xpj3.x -= dr.x * invMassJ;
		            xpj3.y -= dr.y * invMassJ;
		            xpj3.z -= dr.z * invMassJ;
                    converged = false;
                }
            }
            iteration++;
        }

        xpi.x += xi.x;
        xpi.y += xi.y;
        xpi.z += xi.z;
        xpj1.x += xj1.x;
        xpj1.y += xj1.y;
        xpj1.z += xj1.z;
        xpj2.x += xj2.x;
        xpj2.y += xj2.y;
        xpj2.z += xj2.z;
        xpj3.x += xj3.x;
        xpj3.y += xj3.y;
        xpj3.z += xj3.z;

        cSim.pPosq[atomID.x] = xpi;
        cSim.pPosq[atomID.y] = xpj1;

        if (atomID.z != -1)
            cSim.pPosq[atomID.z] = xpj2;

        if (atomID.w != -1)
            cSim.pPosq[atomID.w] = xpj3;

        pos += blockDim.x * gridDim.x;
    }
}

__global__ void kApplyNoShake_kernel()
{
    unsigned int pos = threadIdx.x + blockIdx.x * blockDim.x;
    while (pos < cSim.NonShakeConstraints)
    {
        int  atomID          = cSim.pNonShakeID[pos];
        float4 apos          = cSim.pOldPosq[atomID];
        float4 xpi           = cSim.pPosq[atomID];
        xpi.x               += apos.x;
        xpi.y               += apos.y;
        xpi.z               += apos.z;
        cSim.pPosq[atomID]   = xpi;

        pos += blockDim.x * gridDim.x;
    }
}


void kApplySecondShake(gpuContext gpu)
{
  //  printf("kApplySecondShake\n");
    if (gpu->sim.ShakeConstraints > 0)
    {
        kApplySecondShake_kernel<<<gpu->sim.blocks, gpu->sim.shake_threads_per_block>>>();
        LAUNCHERROR("kApplySecondShake");
    }

    // handle non-Shake atoms

    if (gpu->sim.NonShakeConstraints > 0)
    {
        kApplyNoShake_kernel<<<gpu->sim.blocks, gpu->sim.nonshake_threads_per_block>>>();
        LAUNCHERROR("kApplyNoShake");
    }
}

