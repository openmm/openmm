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
#include <fstream>
using namespace std;

#include "gputypes.h"

#define UNROLLXX 0
#define UNROLLXY 0

struct Atom {
    float x;
    float y;
    float z;
    float r;
    float sr;
    float sum;
    float gamma;
};

static __constant__ cudaGmxSimulation cSim;

void SetCalculateGBVIBornSumSim(gpuContext gpu)
{
    cudaError_t status;
    status = cudaMemcpyToSymbol(cSim, &gpu->sim, sizeof(cudaGmxSimulation));     
    RTERROR(status, "cudaMemcpyToSymbol: SetSim copy to cSim failed");
}

void GetCalculateGBVIBornSumSim(gpuContext gpu)
{
    cudaError_t status;
    status = cudaMemcpyFromSymbol(&gpu->sim, cSim, sizeof(cudaGmxSimulation));     
    RTERROR(status, "cudaMemcpyFromSymbol: SetSim copy from cSim failed");
}

// Include versions of the kernels for N^2 calculations.

#define METHOD_NAME(a, b) a##N2##b
#include "kCalculateGBVIBornSum.h"
#define USE_OUTPUT_BUFFER_PER_WARP
#undef METHOD_NAME
#define METHOD_NAME(a, b) a##N2ByWarp##b
#include "kCalculateGBVIBornSum.h"

// Include versions of the kernels with cutoffs.

#undef METHOD_NAME
#undef USE_OUTPUT_BUFFER_PER_WARP
#define USE_CUTOFF
#define METHOD_NAME(a, b) a##Cutoff##b
#include "kCalculateGBVIBornSum.h"
#define USE_OUTPUT_BUFFER_PER_WARP
#undef METHOD_NAME
#define METHOD_NAME(a, b) a##CutoffByWarp##b
#include "kCalculateGBVIBornSum.h"

// Include versions of the kernels with periodic boundary conditions.

#undef METHOD_NAME
#undef USE_OUTPUT_BUFFER_PER_WARP
#define USE_PERIODIC
#define METHOD_NAME(a, b) a##Periodic##b
#include "kCalculateGBVIBornSum.h"
#define USE_OUTPUT_BUFFER_PER_WARP
#undef METHOD_NAME
#define METHOD_NAME(a, b) a##PeriodicByWarp##b
#include "kCalculateGBVIBornSum.h"

/**---------------------------------------------------------------------------------------

   Compute quintic spline value and associated derviative

   @param x                   value to compute spline at
   @param rl                  lower cutoff value
   @param ru                  upper cutoff value
   @param outValue            value of spline at x
   @param outDerivative       value of derivative of spline at x

   --------------------------------------------------------------------------------------- */

static __device__ void quinticSpline_kernel( float x, float rl, float ru,
                                             float* outValue, float* outDerivative ){

   // ---------------------------------------------------------------------------------------

   const float one           =    1.0f;
   const float minusSix      =   -6.0f;
   const float minusTen      =  -10.0f;
   const float minusThirty   =  -30.0f;
   const float fifteen       =   15.0f;
   const float sixty         =   60.0f;

   // ---------------------------------------------------------------------------------------

   float numerator    = x  - rl;
   float denominator  = ru - rl;
   float ratio        = numerator/denominator;
   float ratio2       = ratio*ratio;
   float ratio3       = ratio2*ratio;

   *outValue               = one + ratio3*(minusTen + fifteen*ratio + minusSix*ratio2);
   *outDerivative          = ratio2*(minusThirty + sixty*ratio + minusThirty*ratio2)/denominator;
}

/**---------------------------------------------------------------------------------------

   Compute Born radii based on Eq. 3 of Labute paper [JCC 29 p. 1693-1698 2008])
   and quintic splice switching function

   @param atomicRadius3       atomic radius cubed
   @param bornSum             Born sum (volume integral)
   @param bornRadius          output Born radius
   @param switchDeriviative   output switching function deriviative

   --------------------------------------------------------------------------------------- */

__device__ void computeBornRadiiUsingQuinticSpline( float atomicRadius3, float bornSum,
                                                    float* bornRadius, float* switchDeriviative ){

   // ---------------------------------------------------------------------------------------

   const float zero          =   0.0f;
   const float one           =   1.0f;
   const float minusOneThird =  (-1.0f/3.0f);

   // ---------------------------------------------------------------------------------------

   // R                = [ S(V)*(A - V) ]**(-1/3)

   // S(V)             = 1                                 V < L
   // S(V)             = qSpline + U/(A-V)                 L < V < A
   // S(V)             = U/(A-V)                           U < V 

   // dR/dr            = (-1/3)*[ S(V)*(A - V) ]**(-4/3)*[ d{ S(V)*(A-V) }/dr

   // d{ S(V)*(A-V) }/dr   = (dV/dr)*[ (A-V)*dS/dV - S(V) ]

   //  (A - V)*dS/dV - S(V)  = 0 - 1                             V < L

   //  (A - V)*dS/dV - S(V)  = (A-V)*d(qSpline) + (A-V)*U/(A-V)**2 - qSpline - U/(A-V) 

	//                        = (A-V)*d(qSpline) - qSpline        L < V < A**(-3)

   //  (A - V)*dS/dV - S(V)  = (A-V)*U*/(A-V)**2 - U/(A-V) = 0   U < V

   float splineL          = cSim.gbviQuinticLowerLimitFactor*atomicRadius3;
   float sum;
   if( bornSum > splineL ){
      if( bornSum < atomicRadius3 ){
         float splineValue, splineDerivative;
         quinticSpline_kernel( bornSum, splineL, atomicRadius3, &splineValue, &splineDerivative ); 
         sum                 = (atomicRadius3 - bornSum)*splineValue + cSim.gbviQuinticUpperBornRadiusLimit;
         *switchDeriviative  = splineValue - (atomicRadius3 - bornSum)*splineDerivative;
      } else {   
         sum                 = cSim.gbviQuinticUpperBornRadiusLimit;
         *switchDeriviative  = zero;
      }
   } else {
      sum                = atomicRadius3 - bornSum; 
      *switchDeriviative = one;
   }
   *bornRadius = pow( sum, minusOneThird );
}

__global__ void kReduceGBVIBornSum_kernel()
{
    unsigned int pos = (blockIdx.x * blockDim.x + threadIdx.x);
    
    while (pos < cSim.atoms)
    {
        float sum = 0.0f;
        float* pSt = cSim.pBornSum + pos;
        float4 atom = cSim.pGBVIData[pos];
        
        // Get summed Born data
        for (int i = 0; i < cSim.nonbondOutputBuffers; i++)
        {
            sum += *pSt;
            pSt += cSim.stride;
        }
        
        // Now calculate Born radius

        float Rinv           = 1.0f/atom.x;
        Rinv                 = Rinv*Rinv*Rinv;
        if( cSim.gbviBornRadiusScalingMethod == 0 ){
            sum                             = Rinv - sum; 
            cSim.pBornRadii[pos]            = pow( sum, (-1.0f/3.0f) ); 
            cSim.pGBVISwitchDerivative[pos] = 1.0f; 
        } else {
            float bornRadius;
            float switchDeriviative;
            computeBornRadiiUsingQuinticSpline( Rinv, sum, &bornRadius, &switchDeriviative );
            cSim.pBornRadii[pos]             = bornRadius; 
            cSim.pGBVISwitchDerivative[pos]  = switchDeriviative; 
        }
        pos += gridDim.x * blockDim.x;
    }   
}

void kReduceGBVIBornSum(gpuContext gpu)
{
    //printf("kReduceGBVIBornSum\n");
    kReduceGBVIBornSum_kernel<<<gpu->sim.blocks, 384>>>();
    gpu->bRecalculateBornRadii = false;
    LAUNCHERROR("kReduceGBVIBornSum");
}

void kPrintGBVI( gpuContext gpu, std::string callId, int call, FILE* log)
{

    gpu->psGBVIData->Download();
    gpu->psBornRadii->Download();
    gpu->psBornForce->Download();
    gpu->psPosq4->Download();
    gpu->psSigEps2->Download();

    (void) fprintf( log, "kPrintGBVI Cuda comp bR bF prm    sigeps2\n" );
    (void) fprintf( stderr, "kCalculateGBVIBornSum: bOutputBufferPerWarp=%u blks=%u th/blk=%u wu=%u %u shrd=%u\n", gpu->bOutputBufferPerWarp,
                    gpu->sim.nonbond_blocks, gpu->sim.nonbond_threads_per_block, gpu->sim.workUnits, gpu->psWorkUnit->_pSysStream[0][0],
                    sizeof(Atom)*gpu->sim.nonbond_threads_per_block );
    for( int ii = 0; ii < gpu->sim.paddedNumberOfAtoms; ii++ ){
        (void) fprintf( log, "%6d %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e \n", ii,
                        gpu->psBornRadii->_pSysData[ii],
                        gpu->psBornForce->_pSysData[ii],

                        gpu->psGBVIData->_pSysData[ii].x,
                        gpu->psGBVIData->_pSysData[ii].y,
                        gpu->psGBVIData->_pSysData[ii].z,
                        gpu->psGBVIData->_pSysData[ii].w,

                        gpu->psSigEps2->_pSysData[ii].x,
                        gpu->psSigEps2->_pSysData[ii].y );

    }   

}

void kCalculateGBVIBornSum(gpuContext gpu)
{
    //printf("kCalculateGBVIBornSum\n");

    switch (gpu->sim.nonbondedMethod)
    {
        case NO_CUTOFF:
            if (gpu->bOutputBufferPerWarp){
                kCalculateGBVIN2ByWarpBornSum_kernel<<<gpu->sim.nonbond_blocks, gpu->sim.nonbond_threads_per_block,
                        sizeof(Atom)*gpu->sim.nonbond_threads_per_block>>>(gpu->sim.pWorkUnit);
            } else {
                kCalculateGBVIN2BornSum_kernel<<<gpu->sim.nonbond_blocks, gpu->sim.nonbond_threads_per_block,
                        sizeof(Atom)*gpu->sim.nonbond_threads_per_block>>>(gpu->sim.pWorkUnit);
            }
            break;

        case CUTOFF:
            if (gpu->bOutputBufferPerWarp)
                kCalculateGBVICutoffByWarpBornSum_kernel<<<gpu->sim.nonbond_blocks, gpu->sim.nonbond_threads_per_block,
                        (sizeof(Atom)+sizeof(float))*gpu->sim.nonbond_threads_per_block>>>(gpu->sim.pInteractingWorkUnit);
            else
                kCalculateGBVICutoffBornSum_kernel<<<gpu->sim.nonbond_blocks, gpu->sim.nonbond_threads_per_block,
                        (sizeof(Atom)+sizeof(float))*gpu->sim.nonbond_threads_per_block>>>(gpu->sim.pInteractingWorkUnit );
            break;

        case PERIODIC:
            if (gpu->bOutputBufferPerWarp)
                kCalculateGBVIPeriodicByWarpBornSum_kernel<<<gpu->sim.nonbond_blocks, gpu->sim.nonbond_threads_per_block,
                        (sizeof(Atom)+sizeof(float))*gpu->sim.nonbond_threads_per_block>>>(gpu->sim.pInteractingWorkUnit );
            else
                kCalculateGBVIPeriodicBornSum_kernel<<<gpu->sim.nonbond_blocks, gpu->sim.nonbond_threads_per_block,
                        (sizeof(Atom)+sizeof(float))*gpu->sim.nonbond_threads_per_block>>>(gpu->sim.pInteractingWorkUnit );
            break;

    }
    LAUNCHERROR("kCalculateGBVIBornSum");
}
