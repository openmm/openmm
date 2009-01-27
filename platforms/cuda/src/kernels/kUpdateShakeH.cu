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

#include <stdio.h>
#include <cuda.h>
#include <vector_functions.h>
#include <cstdlib>
#include <string>
#include <iostream>
//#include <fstream>
using namespace std;

#define DeltaShake

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

void SetUpdateShakeHSim(gpuContext gpu)
{
    cudaError_t status;
    status = cudaMemcpyToSymbol(cSim, &gpu->sim, sizeof(cudaGmxSimulation));     
    RTERROR(status, "cudaMemcpyToSymbol: SetSim copy to cSim failed");
}

void GetUpdateShakeHSim(gpuContext gpu)
{
    cudaError_t status;
    status = cudaMemcpyFromSymbol(&gpu->sim, cSim, sizeof(cudaGmxSimulation));     
    RTERROR(status, "cudaMemcpyFromSymbol: SetSim copy from cSim failed");
}

__global__ void kUpdatePart1_kernel()
{
    unsigned int pos    = threadIdx.x + blockIdx.x * blockDim.x;
    unsigned int rpos   = cSim.pRandomPosition[blockIdx.x];
    __syncthreads();
    
    while (pos < cSim.atoms)
    {
        float4 velocity         = cSim.pVelm4[pos];
        float4 xVector          = cSim.pxVector4[pos];
        float4 random4a         = cSim.pRandom4a[rpos + pos];
        float2 random2a         = cSim.pRandom2a[rpos + pos];
        float4 apos             = cSim.pPosq[pos];
        float4 force            = cSim.pForce4[pos];
        
        float3 Vmh;
        float sqrtInvMass       = sqrt(velocity.w);
        Vmh.x                   = xVector.x * cSim.DOverTauC + sqrtInvMass * random4a.x;
        Vmh.y                   = xVector.y * cSim.DOverTauC + sqrtInvMass * random4a.y;
        Vmh.z                   = xVector.z * cSim.DOverTauC + sqrtInvMass * random4a.z;
        float4 vVector;
        vVector.x               = sqrtInvMass * random4a.w;
        vVector.y               = sqrtInvMass * random2a.x;
        vVector.z               = sqrtInvMass * random2a.y;
        vVector.w               = 0.0f;
        cSim.pvVector4[pos]     = vVector;
        velocity.x              = velocity.x * cSim.EM + 
                                  velocity.w * force.x * cSim.TauOneMinusEM +
                                  vVector.x -
                                  cSim.EM * Vmh.x;
        velocity.y              = velocity.y * cSim.EM + 
                                  velocity.w * force.y * cSim.TauOneMinusEM +
                                  vVector.y -
                                  cSim.EM * Vmh.y;
        velocity.z              = velocity.z * cSim.EM + 
                                  velocity.w * force.z * cSim.TauOneMinusEM +
                                  vVector.z -
                                  cSim.EM * Vmh.z;
        cSim.pOldPosq[pos]      = apos;
#ifndef DeltaShake
        apos.x                 += velocity.x * cSim.fix1;
        apos.y                 += velocity.y * cSim.fix1;
        apos.z                 += velocity.z * cSim.fix1;
#else
        apos.x                  = velocity.x * cSim.fix1;
        apos.y                  = velocity.y * cSim.fix1;
        apos.z                  = velocity.z * cSim.fix1;
#endif
        cSim.pPosqP[pos]        = apos;
        cSim.pVelm4[pos]        = velocity;        
        pos                    += blockDim.x * gridDim.x;
    }
}

__global__ void kUpdatePart1CM_kernel()
{
    extern __shared__ float3 sCM[];
    unsigned int pos    = threadIdx.x + blockIdx.x * blockDim.x;
    unsigned int rpos   = cSim.pRandomPosition[blockIdx.x];
    float3 CM           = { 0.0f, 0.0f, 0.0f};
    float4 CM1          = { 0.0f, 0.0f, 0.0f, 0.0f };
    
    // Read CM outputs from previous step
    unsigned int cpos = threadIdx.x;
#if 0
    float4 CM2          = { 0.0f, 0.0f, 0.0f, 0.0f };
    float4 CM3          = { 0.0f, 0.0f, 0.0f, 0.0f };
    float4 CM4          = { 0.0f, 0.0f, 0.0f, 0.0f };
    if (cpos < gridDim.x)
        CM1             = cSim.pLinearMomentum[cpos];
    cpos               += gridDim.x;
    if (cpos < gridDim.x)
        CM2             = cSim.pLinearMomentum[cpos];
    cpos               += gridDim.x;
    if (cpos < gridDim.x)
        CM3             = cSim.pLinearMomentum[cpos];
    cpos               += gridDim.x;
    if (cpos < gridDim.x)
        CM4             = cSim.pLinearMomentum[cpos];
    sCM[threadIdx.x].x  = CM1.x + CM2.x + CM3.x + CM4.x;
    sCM[threadIdx.x].y  = CM1.y + CM2.y + CM3.y + CM4.y;
    sCM[threadIdx.x].z  = CM1.z + CM2.z + CM3.z + CM4.z;
#else
    while (cpos < gridDim.x)
    {
        CM1             = cSim.pLinearMomentum[cpos];
        CM.x           += CM1.x;
        CM.y           += CM1.y;
        CM.z           += CM1.z;
        cpos           += blockDim.x;
    }
    sCM[threadIdx.x].x  = CM.x;
    sCM[threadIdx.x].y  = CM.y;
    sCM[threadIdx.x].z  = CM.z;
#endif
    __syncthreads();
    
    // Reduce CM
    unsigned int offset = 1;
    unsigned int mask   = 1;
    while (offset < blockDim.x)
    {
        if (((threadIdx.x & mask) == 0) && (threadIdx.x + offset < blockDim.x))
        {
            sCM[threadIdx.x].x += sCM[threadIdx.x + offset].x;
            sCM[threadIdx.x].y += sCM[threadIdx.x + offset].y;
            sCM[threadIdx.x].z += sCM[threadIdx.x + offset].z;
        }
        mask = 2 * mask + 1;
        offset *= 2;
        __syncthreads();
    }       
    
    while (pos < cSim.atoms)
    {
        float4 velocity         = cSim.pVelm4[pos];
        float4 xVector          = cSim.pxVector4[pos];
        float4 random4a         = cSim.pRandom4a[rpos + pos];
        float2 random2a         = cSim.pRandom2a[rpos + pos];
        float4 apos             = cSim.pPosq[pos];
        float4 force            = cSim.pForce4[pos];
        
        float3 Vmh;
        float sqrtInvMass       = sqrt(velocity.w);
        Vmh.x                   = xVector.x * cSim.DOverTauC + sqrtInvMass * random4a.x;
        Vmh.y                   = xVector.y * cSim.DOverTauC + sqrtInvMass * random4a.y;
        Vmh.z                   = xVector.z * cSim.DOverTauC + sqrtInvMass * random4a.z;
        float4 vVector;
        vVector.x               = sqrtInvMass * random4a.w;
        vVector.y               = sqrtInvMass * random2a.x;
        vVector.z               = sqrtInvMass * random2a.y;
        vVector.w               = 0.0f;
        cSim.pvVector4[pos]     = vVector;
        velocity.x              = velocity.x * cSim.EM + 
                                  velocity.w * force.x * cSim.TauOneMinusEM +
                                  vVector.x -
                                  cSim.EM * Vmh.x -
                                  sCM[0].x;
        velocity.y              = velocity.y * cSim.EM + 
                                  velocity.w * force.y * cSim.TauOneMinusEM +
                                  vVector.y -
                                  cSim.EM * Vmh.y -
                                  sCM[0].y;
        velocity.z              = velocity.z * cSim.EM + 
                                  velocity.w * force.z * cSim.TauOneMinusEM +
                                  vVector.z -
                                  cSim.EM * Vmh.z -
                                  sCM[0].z;
        cSim.pOldPosq[pos]      = apos;
#ifndef DeltaShake
        apos.x                 += velocity.x * cSim.fix1;
        apos.y                 += velocity.y * cSim.fix1;
        apos.z                 += velocity.z * cSim.fix1;
#else
        apos.x                  = velocity.x * cSim.fix1;
        apos.y                  = velocity.y * cSim.fix1;
        apos.z                  = velocity.z * cSim.fix1;
#endif
        cSim.pPosqP[pos]        = apos;
        cSim.pVelm4[pos]        = velocity;        
        pos                    += blockDim.x * gridDim.x;
    }
}



void kUpdatePart1(gpuContext gpu)
{
//    printf("kUpdatePart1\n");
#if 0
    static int iteration = 0;
    if (iteration == 0)
    {
        gpu->psPosq4->Download();
        gpu->psVelm4->Download();
        printf("# %d atoms\n", gpu->natoms);
        for (int i = 0; i < gpu->natoms; i++)
        {
            printf("%5d %11.5f %11.5f %11.5f %11.5f %11.5f %11.5f %11.5f %11.5f\n", i,
                gpu->psPosq4->_pSysStream[0][i].x, gpu->psPosq4->_pSysStream[0][i].y,
                gpu->psPosq4->_pSysStream[0][i].z, gpu->psPosq4->_pSysStream[0][i].w,
                gpu->psVelm4->_pSysStream[0][i].x, gpu->psVelm4->_pSysStream[0][i].y,
                gpu->psVelm4->_pSysStream[0][i].z, gpu->psVelm4->_pSysStream[0][i].w
            );       
        }
    }
    iteration++;
#endif
#if 0
    static const float KILO 		        =    1e3;              		// Thousand
    static const float BOLTZMANN	        =    1.380658e-23f;            // (J/K)	
    static const float AVOGADRO	            =    6.0221367e23f;		    // ()		
    static const float RGAS                 =    BOLTZMANN * AVOGADRO;     // (J/(mol K))
    static const float BOLTZ                =    (RGAS / KILO);            // (kJ/(mol K)) 
    static int iteration = 0;

    // Check T
    if (iteration % 1000 == 0)
    {
        gpu->psVelm4->Download();
        float ke = 0.0f;
        for (int i = 0; i < gpu->natoms; i++)
        {
            float vx = gpu->psVelm4->_pSysStream[0][i].x;
            float vy = gpu->psVelm4->_pSysStream[0][i].y;
            float vz = gpu->psVelm4->_pSysStream[0][i].z;
            float m = 1.0f / gpu->psVelm4->_pSysStream[0][i].w;
            ke += m * (vx * vx + vy * vy + vz * vz);
        }
        float T = ke / (BOLTZ  * gpu->sim.degreesOfFreedom);
        printf("Iteration %d, Temperature is %f\n", iteration, T);
    }
    iteration++;
#endif    
    if (gpu->bRemoveCM)
    {
        kUpdatePart1CM_kernel<<<gpu->sim.blocks, gpu->sim.update_threads_per_block, gpu->sim.update_threads_per_block * sizeof(float3)>>>();
        LAUNCHERROR("kUpdatePart1CM");
        gpu->bRemoveCM = false;

#if 0
        gpu->psLinearMomentum->Download();
        gpu->psVelm4->Download();
        float3 mv = {0.0f, 0.0f, 0.0f};
        for (int i = 0; i < gpu->natoms; i++)
        {
            float mass = 1.0f / gpu->psVelm4->_pSysStream[0][i].w;
            mv.x += mass * gpu->psVelm4->_pSysStream[0][i].x;
            mv.y += mass * gpu->psVelm4->_pSysStream[0][i].y;
            mv.z += mass * gpu->psVelm4->_pSysStream[0][i].z;
        }
        mv.x *= gpu->sim.inverseTotalMass;
        mv.y *= gpu->sim.inverseTotalMass;
        mv.z *= gpu->sim.inverseTotalMass;
        
        float3 mv1 = {0.0f, 0.0f, 0.0f};
        for (int i = 0; i < gpu->sim.blocks; i++)
        {
            mv1.x += gpu->psLinearMomentum->_pSysStream[0][i].x;
            mv1.y += gpu->psLinearMomentum->_pSysStream[0][i].y;
            mv1.z += gpu->psLinearMomentum->_pSysStream[0][i].z;
        }
        printf("%11.5f %11.5f %11.5f | %11.5f %11.5f %11.5f\n", mv.x, mv.y, mv.z, mv1.x, mv1.y, mv1.z);
#endif
    }
    else
    {    
        kUpdatePart1_kernel<<<gpu->sim.blocks, gpu->sim.update_threads_per_block>>>();
        LAUNCHERROR("kUpdatePart1");
    }
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
#ifndef DeltaShake
        xpi.x              -= xi.x;
        xpi.y              -= xi.y;
        xpi.z              -= xi.z;
        xpj1.x             -= xj1.x;
        xpj1.y             -= xj1.y;
        xpj1.z             -= xj1.z;
        xpj2.x             -= xj2.x;
        xpj2.y             -= xj2.y;
        xpj2.z             -= xj2.z;
        xpj3.x             -= xj3.x;
        xpj3.y             -= xj3.y;
        xpj3.z             -= xj3.z;
#endif
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
        
#ifndef DeltaShake
        xpi.x  += xi.x;
        xpi.y  += xi.y;
        xpi.z  += xi.z;

        xpj1.x += xj1.x;
        xpj1.y += xj1.y;
        xpj1.z += xj1.z;
        
        xpj2.x += xj2.x;
        xpj2.y += xj2.y;
        xpj2.z += xj2.z;
        
        xpj3.x += xj3.x;
        xpj3.y += xj3.y;
        xpj3.z += xj3.z;
#endif
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


__global__ void kUpdatePart2_kernel()
{
    unsigned int pos            = threadIdx.x + blockIdx.x * blockDim.x;
    unsigned int rpos           = cSim.pRandomPosition[blockIdx.x];
    __syncthreads();
    
    while (pos < cSim.atoms)
    {
        float4 velocity         = cSim.pVelm4[pos];
#ifndef DeltaShake
        float4 apos             = cSim.pPosq[pos];
#endif
        float4 xPrime           = cSim.pPosqP[pos];
        float4 vVector          = cSim.pvVector4[pos];
        float4 xVector;
        float4 random4b         = cSim.pRandom4b[rpos + pos];
        float2 random2b         = cSim.pRandom2b[rpos + pos];
        float3 Xmh;
        
        float sqrtInvMass       = sqrt(velocity.w);
#ifdef DeltaShake
        velocity.x              = xPrime.x * cSim.oneOverFix1;
        velocity.y              = xPrime.y * cSim.oneOverFix1;
        velocity.z              = xPrime.z * cSim.oneOverFix1;
#else
        velocity.x              = (xPrime.x - apos.x) * cSim.oneOverFix1;
        velocity.y              = (xPrime.y - apos.y) * cSim.oneOverFix1;
        velocity.z              = (xPrime.z - apos.z) * cSim.oneOverFix1;
#endif
        Xmh.x                   = vVector.x * cSim.TauDOverEMMinusOne +
                                  sqrtInvMass * random4b.x;
        Xmh.y                   = vVector.y * cSim.TauDOverEMMinusOne +
                                  sqrtInvMass * random4b.y;
        Xmh.z                   = vVector.z * cSim.TauDOverEMMinusOne +
                                  sqrtInvMass * random4b.z;
        xVector.x               = sqrtInvMass * random4b.w;
        xVector.y               = sqrtInvMass * random2b.x;
        xVector.z               = sqrtInvMass * random2b.y;                    
        xPrime.x               += xVector.x - Xmh.x;
        xPrime.y               += xVector.y - Xmh.y;
        xPrime.z               += xVector.z - Xmh.z;
        
    
        cSim.pPosq[pos]         = xPrime;
        cSim.pVelm4[pos]        = velocity;
        cSim.pxVector4[pos]     = xVector;
         
        pos                    += blockDim.x * gridDim.x;    
    }

    // Update random position pointer
    if (threadIdx.x == 0)
    {
        rpos                   += cSim.paddedNumberOfAtoms;
        if (rpos > cSim.randoms)
            rpos               -= cSim.randoms;
        cSim.pRandomPosition[blockIdx.x] = rpos;
    }
}

__global__ void kUpdatePart2CM_kernel()
{
    extern __shared__ float3 sCM[];
    unsigned int pos            = threadIdx.x + blockIdx.x * blockDim.x;
    unsigned int rpos           = cSim.pRandomPosition[blockIdx.x];
    float3 CM                   = {0.0f, 0.0f, 0.0f};
    __syncthreads();
    
    while (pos < cSim.atoms)
    {
        float4 velocity         = cSim.pVelm4[pos];
#ifndef DeltaShake
        float4 apos             = cSim.pPosq[pos];
#endif
        float4 xPrime           = cSim.pPosqP[pos];
        float4 vVector          = cSim.pvVector4[pos];
        float4 xVector;
        float4 random4b         = cSim.pRandom4b[rpos + pos];
        float2 random2b         = cSim.pRandom2b[rpos + pos];
        float3 Xmh;
        float mass              = 1.0f / velocity.w;
        float sqrtInvMass       = sqrt(velocity.w);
#ifdef DeltaShake
        velocity.x              = xPrime.x * cSim.oneOverFix1;
        velocity.y              = xPrime.y * cSim.oneOverFix1;
        velocity.z              = xPrime.z * cSim.oneOverFix1;
#else
        velocity.x              = (xPrime.x - apos.x) * cSim.oneOverFix1;
        velocity.y              = (xPrime.y - apos.y) * cSim.oneOverFix1;
        velocity.z              = (xPrime.z - apos.z) * cSim.oneOverFix1;
#endif
        CM.x                   += mass * velocity.x;
        CM.y                   += mass * velocity.y;
        CM.z                   += mass * velocity.z;
        
        Xmh.x                   = vVector.x * cSim.TauDOverEMMinusOne +
                                  sqrtInvMass * random4b.x;
        Xmh.y                   = vVector.y * cSim.TauDOverEMMinusOne +
                                  sqrtInvMass * random4b.y;
        Xmh.z                   = vVector.z * cSim.TauDOverEMMinusOne +
                                  sqrtInvMass * random4b.z;
        xVector.x               = sqrtInvMass * random4b.w;
        xVector.y               = sqrtInvMass * random2b.x;
        xVector.z               = sqrtInvMass * random2b.y;                    
        xPrime.x               += xVector.x - Xmh.x;
        xPrime.y               += xVector.y - Xmh.y;
        xPrime.z               += xVector.z - Xmh.z;
        
    
        cSim.pPosq[pos]         = xPrime;
        cSim.pVelm4[pos]        = velocity;
        cSim.pxVector4[pos]     = xVector;
        
        pos                    += blockDim.x * gridDim.x;    
    }

    // Update random position pointer
    if (threadIdx.x == 0)
    {
        rpos                   += cSim.paddedNumberOfAtoms;
        if (rpos > cSim.randoms)
            rpos               -= cSim.randoms;
        cSim.pRandomPosition[blockIdx.x] = rpos;
    }
    
    // Scale CM
    CM.x *= cSim.inverseTotalMass;
    CM.y *= cSim.inverseTotalMass;
    CM.z *= cSim.inverseTotalMass;
    sCM[threadIdx.x] = CM;
    __syncthreads();
    
    // Reduce CM for CTA
    unsigned int offset = 1;
    unsigned int mask   = 1;
    while (offset < blockDim.x)
    {
        if (((threadIdx.x & mask) == 0) && (threadIdx.x + offset < blockDim.x))
        {
            sCM[threadIdx.x].x += sCM[threadIdx.x + offset].x;
            sCM[threadIdx.x].y += sCM[threadIdx.x + offset].y;
            sCM[threadIdx.x].z += sCM[threadIdx.x + offset].z;
        }
        mask = 2 * mask + 1;
        offset *= 2;
        __syncthreads();
    }
    if (threadIdx.x == 0)
    {
        float4 CM;
        CM.x                                = sCM[0].x;
        CM.y                                = sCM[0].y;
        CM.z                                = sCM[0].z;
        CM.w                                = 0.0f;
        cSim.pLinearMomentum[blockIdx.x]    = CM;
    }  
}

extern void kGenerateRandoms(gpuContext gpu);
void kUpdatePart2(gpuContext gpu)
{
//    printf("kUpdatePart2\n");
    if (gpu->bCalculateCM)
    {
        kUpdatePart2CM_kernel<<<gpu->sim.blocks, gpu->sim.update_threads_per_block, gpu->sim.update_threads_per_block * sizeof(float3)>>>();
        LAUNCHERROR("kUpdatePart2CM");
        gpu->bCalculateCM = false;
        gpu->bRemoveCM = true;
       
#if 0
        gpu->psLinearMomentum->Download();
        gpu->psVelm4->Download();
        float3 mv = {0.0f, 0.0f, 0.0f};
        for (int i = 0; i < gpu->natoms; i++)
        {
            float mass = 1.0f / gpu->psVelm4->_pSysStream[0][i].w;
            mv.x += mass * gpu->psVelm4->_pSysStream[0][i].x;
            mv.y += mass * gpu->psVelm4->_pSysStream[0][i].y;
            mv.z += mass * gpu->psVelm4->_pSysStream[0][i].z;
        }
        mv.x *= gpu->sim.inverseTotalMass;
        mv.y *= gpu->sim.inverseTotalMass;
        mv.z *= gpu->sim.inverseTotalMass;
        
        float3 mv1 = {0.0f, 0.0f, 0.0f};
        for (int i = 0; i < gpu->sim.blocks; i++)
        {
            mv1.x += gpu->psLinearMomentum->_pSysStream[0][i].x;
            mv1.y += gpu->psLinearMomentum->_pSysStream[0][i].y;
            mv1.z += gpu->psLinearMomentum->_pSysStream[0][i].z;
        }
        printf("%11.5f %11.5f %11.5f | %11.5f %11.5f %11.5f\n", mv.x, mv.y, mv.z, mv1.x, mv1.y, mv1.z);
#endif
    }
    else
    {
        kUpdatePart2_kernel<<<gpu->sim.blocks, gpu->sim.update_threads_per_block>>>();
        LAUNCHERROR("kUpdatePart2");
    }
    
    // Update randoms if necessary
    static int iteration = 0;
    iteration++;
    if (iteration == gpu->sim.randomIterations)
    {
        kGenerateRandoms(gpu);
        iteration = 0;
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
#ifndef DeltaShake
        xpi.x              -= xi.x;
        xpi.y              -= xi.y;
        xpi.z              -= xi.z;
        xpj1.x             -= xj1.x;
        xpj1.y             -= xj1.y;
        xpj1.z             -= xj1.z;
        xpj2.x             -= xj2.x;
        xpj2.y             -= xj2.y;
        xpj2.z             -= xj2.z;
        xpj3.x             -= xj3.x;
        xpj3.y             -= xj3.y;
        xpj3.z             -= xj3.z;
#endif
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

void kCPUShake2(gpuContext gpu)
{

}

void kApplySecondShake(gpuContext gpu)
{
  //  printf("kApplySecondShake\n");
  //  kCPUShake2(gpu);
    if (gpu->sim.ShakeConstraints > 0)
    {
        kApplySecondShake_kernel<<<gpu->sim.blocks, gpu->sim.shake_threads_per_block>>>();
        LAUNCHERROR("kApplySecondShake");
    }

    // handle non-Shake atoms

#ifdef DeltaShake
    if (gpu->sim.NonShakeConstraints > 0)
    {
        //fprintf( gpu->log, "kApplyNoShake_kernel %d %d \n", gpu->sim.blocks, gpu->sim.nonshake_threads_per_block); fflush( gpu->log );
        kApplyNoShake_kernel<<<gpu->sim.blocks, gpu->sim.nonshake_threads_per_block>>>();
        LAUNCHERROR("kApplyNoShake");
    }
#endif

}

