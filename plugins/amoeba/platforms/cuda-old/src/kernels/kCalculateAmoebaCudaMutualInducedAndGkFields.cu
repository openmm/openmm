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

#include "amoebaGpuTypes.h"
#include "amoebaCudaKernels.h"
#include "kCalculateAmoebaCudaUtilities.h"
#include "openmm/OpenMMException.h"

#include <stdio.h>
#include <sstream>

using namespace std;

static __constant__ cudaGmxSimulation cSim;
static __constant__ cudaAmoebaGmxSimulation cAmoebaSim;

void SetCalculateAmoebaCudaMutualInducedAndGkFieldsSim(amoebaGpuContext amoebaGpu)
{
    cudaError_t status;
    gpuContext gpu = amoebaGpu->gpuContext;
    status         = cudaMemcpyToSymbol(cSim, &gpu->sim, sizeof(cudaGmxSimulation));    
    RTERROR(status, "SetCalculateAmoebaCudaMutualInducedAndGkFieldSim: cudaMemcpyToSymbol: SetSim copy to cSim failed");
    status         = cudaMemcpyToSymbol(cAmoebaSim, &amoebaGpu->amoebaSim, sizeof(cudaAmoebaGmxSimulation));    
    RTERROR(status, "SetCalculateAmoebaCudaMutualInducedAndGkFieldSim: cudaMemcpyToSymbol: SetSim copy to cAmoebaSim failed");
}

void GetCalculateAmoebaCudaMutualInducedAndGkFieldsSim(amoebaGpuContext amoebaGpu)
{
    cudaError_t status;
    gpuContext gpu = amoebaGpu->gpuContext;
    status = cudaMemcpyFromSymbol(&gpu->sim, cSim, sizeof(cudaGmxSimulation));    
    RTERROR(status, "GetCalculateAmoebaCudaMutualInducedAndGkFieldSim: cudaMemcpyFromSymbol: SetSim copy from cSim failed");
    status = cudaMemcpyFromSymbol(&amoebaGpu->amoebaSim, cAmoebaSim, sizeof(cudaAmoebaGmxSimulation));    
    RTERROR(status, "GetCalculateAmoebaCudaMutualInducedAndGkFieldSim: cudaMemcpyFromSymbol: SetSim copy from cAmoebaSim failed");
}

#define GK
#include "kCalculateAmoebaCudaMutualInducedParticle.h"
#undef GK

__device__ void calculateMutualInducedAndGkFieldsPairIxn_kernel( MutualInducedParticle& atomI, MutualInducedParticle& atomJ,
                                                                 float fields[8][3] )
{

    float deltaR[3];
    
    // ---------------------------------------------------------------------------------------
    
    // get deltaR, and r between 2 atoms
    
    deltaR[0]                                    = atomJ.x - atomI.x;
    deltaR[1]                                    = atomJ.y - atomI.y;
    deltaR[2]                                    = atomJ.z - atomI.z;

    float r                                      =  sqrtf( deltaR[0]*deltaR[0] + deltaR[1]*deltaR[1] + deltaR[2]*deltaR[2] );
    float rI                                     =  1.0f/r;
    float r2I                                    =  rI*rI;
    float rr3                                    = -rI*r2I;
    float rr5                                    = -3.0f*rr3*r2I;
    
    float dampProd                               = atomI.damp*atomJ.damp;
    float ratio                                  = (dampProd != 0.0f) ? (r/dampProd) : 1.0f;
    float pGamma                                 = atomI.thole > atomJ.thole ? atomJ.thole: atomI.thole;
    float damp                                   = ratio*ratio*ratio*pGamma;
    float dampExp                                = ( (dampProd != 0.0f) && (r < cAmoebaSim.scalingDistanceCutoff) ) ? expf( -damp ) : 0.0f; 

    rr3                                         *= (1.0f - dampExp);
    rr5                                         *= (1.0f - ( 1.0f + damp )*dampExp);
        
    float dDotDelta                              = rr5*(deltaR[0]*atomJ.inducedDipole[0]    + deltaR[1]*atomJ.inducedDipole[1]    + deltaR[2]*atomJ.inducedDipole[2] );
    fields[0][0]                                 = rr3*atomJ.inducedDipole[0] + dDotDelta*deltaR[0];
    fields[0][1]                                 = rr3*atomJ.inducedDipole[1] + dDotDelta*deltaR[1];
    fields[0][2]                                 = rr3*atomJ.inducedDipole[2] + dDotDelta*deltaR[2];
   
    dDotDelta                                    = rr5*(deltaR[0]*atomJ.inducedDipolePolar[0]    + deltaR[1]*atomJ.inducedDipolePolar[1]    + deltaR[2]*atomJ.inducedDipolePolar[2] );
    fields[1][0]                                 = rr3*atomJ.inducedDipolePolar[0] + dDotDelta*deltaR[0];
    fields[1][1]                                 = rr3*atomJ.inducedDipolePolar[1] + dDotDelta*deltaR[1];
    fields[1][2]                                 = rr3*atomJ.inducedDipolePolar[2] + dDotDelta*deltaR[2];
  
    dDotDelta                                    = rr5*(deltaR[0]*atomI.inducedDipole[0]    + deltaR[1]*atomI.inducedDipole[1]    + deltaR[2]*atomI.inducedDipole[2] );
    fields[2][0]                                 = rr3*atomI.inducedDipole[0] + dDotDelta*deltaR[0];
    fields[2][1]                                 = rr3*atomI.inducedDipole[1] + dDotDelta*deltaR[1];
    fields[2][2]                                 = rr3*atomI.inducedDipole[2] + dDotDelta*deltaR[2];
   
    dDotDelta                                    = rr5*(deltaR[0]*atomI.inducedDipolePolar[0]    + deltaR[1]*atomI.inducedDipolePolar[1]    + deltaR[2]*atomI.inducedDipolePolar[2] );
    fields[3][0]                                 = rr3*atomI.inducedDipolePolar[0] + dDotDelta*deltaR[0];
    fields[3][1]                                 = rr3*atomI.inducedDipolePolar[1] + dDotDelta*deltaR[1];
    fields[3][2]                                 = rr3*atomI.inducedDipolePolar[2] + dDotDelta*deltaR[2];

    dDotDelta                                    = rr5*(deltaR[0]*atomJ.inducedDipoleS[0]    + deltaR[1]*atomJ.inducedDipoleS[1]    + deltaR[2]*atomJ.inducedDipoleS[2] );
    fields[4][0]                                 = rr3*atomJ.inducedDipoleS[0] + dDotDelta*deltaR[0];
    fields[4][1]                                 = rr3*atomJ.inducedDipoleS[1] + dDotDelta*deltaR[1];
    fields[4][2]                                 = rr3*atomJ.inducedDipoleS[2] + dDotDelta*deltaR[2];
   
    dDotDelta                                    = rr5*(deltaR[0]*atomJ.inducedDipolePolarS[0]    + deltaR[1]*atomJ.inducedDipolePolarS[1]    + deltaR[2]*atomJ.inducedDipolePolarS[2] );
    fields[5][0]                                 = rr3*atomJ.inducedDipolePolarS[0] + dDotDelta*deltaR[0];
    fields[5][1]                                 = rr3*atomJ.inducedDipolePolarS[1] + dDotDelta*deltaR[1];
    fields[5][2]                                 = rr3*atomJ.inducedDipolePolarS[2] + dDotDelta*deltaR[2];
  
    dDotDelta                                    = rr5*(deltaR[0]*atomI.inducedDipoleS[0]    + deltaR[1]*atomI.inducedDipoleS[1]    + deltaR[2]*atomI.inducedDipoleS[2] );
    fields[6][0]                                 = rr3*atomI.inducedDipoleS[0] + dDotDelta*deltaR[0];
    fields[6][1]                                 = rr3*atomI.inducedDipoleS[1] + dDotDelta*deltaR[1];
    fields[6][2]                                 = rr3*atomI.inducedDipoleS[2] + dDotDelta*deltaR[2];
   
    dDotDelta                                    = rr5*(deltaR[0]*atomI.inducedDipolePolarS[0]    + deltaR[1]*atomI.inducedDipolePolarS[1]    + deltaR[2]*atomI.inducedDipolePolarS[2] );
    fields[7][0]                                 = rr3*atomI.inducedDipolePolarS[0] + dDotDelta*deltaR[0];
    fields[7][1]                                 = rr3*atomI.inducedDipolePolarS[1] + dDotDelta*deltaR[1];
    fields[7][2]                                 = rr3*atomI.inducedDipolePolarS[2] + dDotDelta*deltaR[2];


}

__device__ void calculateMutualInducedAndGkFieldsGkPairIxn_kernel( MutualInducedParticle& atomI, MutualInducedParticle& atomJ,
                                                                   float gkField[8][3] )
{

    float gux[5];
    float guy[5];
    float guz[5];
    float a[3][3];
    
    // ---------------------------------------------------------------------------------------
    
    float xr               = atomJ.x - atomI.x;
    float yr               = atomJ.y - atomI.y;
    float zr               = atomJ.z - atomI.z;

    float xr2              = xr*xr;
    float yr2              = yr*yr;
    float zr2              = zr*zr;

    float rb2              = atomI.bornRadius*atomJ.bornRadius;

    float r2               = xr2 + yr2 + zr2;
    float expterm          = expf(-r2/(cAmoebaSim.gkc*rb2));
    float expc             = expterm /cAmoebaSim.gkc; 
    //float dexpc            = -2.0f / (cAmoebaSim.gkc*rb2);

    float gf2              = 1.0f / (r2+rb2*expterm);
    float gf               = sqrtf(gf2);
    float gf3              = gf2 * gf;
    float gf5              = gf3 * gf2;

    float duixs            = atomI.inducedDipoleS[0];
    float duiys            = atomI.inducedDipoleS[1];
    float duizs            = atomI.inducedDipoleS[2];

    float puixs            = atomI.inducedDipolePolarS[0];
    float puiys            = atomI.inducedDipolePolarS[1];
    float puizs            = atomI.inducedDipolePolarS[2];
 
    float dukxs            = atomJ.inducedDipoleS[0];
    float dukys            = atomJ.inducedDipoleS[1];
    float dukzs            = atomJ.inducedDipoleS[2];

    float pukxs            = atomJ.inducedDipolePolarS[0];
    float pukys            = atomJ.inducedDipolePolarS[1];
    float pukzs            = atomJ.inducedDipolePolarS[2];
 
    // reaction potential auxiliary terms
 
    a[1][0]                = -gf3;
    a[2][0]                = 3.0f * gf5;

    // reaction potential gradient auxiliary terms

    float expc1            = 1.0f - expc;
    a[1][1]                = expc1 * a[2][0];
 
    // unweighted dipole reaction potential gradient tensor

    gux[2]                 = cAmoebaSim.fd * (a[1][0] + xr2*a[1][1]);
    gux[3]                 = cAmoebaSim.fd * xr*yr*a[1][1];
    gux[4]                 = cAmoebaSim.fd * xr*zr*a[1][1];

    guy[2]                 = gux[3];
    guy[3]                 = cAmoebaSim.fd * (a[1][0] + yr2*a[1][1]);
    guy[4]                 = cAmoebaSim.fd * yr*zr*a[1][1];

    guz[2]                 = gux[4];
    guz[3]                 = guy[4];
    guz[4]                 = cAmoebaSim.fd * (a[1][0] + zr2*a[1][1]);
 
    gkField[0][0]          = dukxs*gux[2]+dukys*guy[2]+dukzs*guz[2];
    gkField[0][1]          = dukxs*gux[3]+dukys*guy[3]+dukzs*guz[3];
    gkField[0][2]          = dukxs*gux[4]+dukys*guy[4]+dukzs*guz[4];

    gkField[1][0]          = duixs*gux[2]+duiys*guy[2]+duizs*guz[2];
    gkField[1][1]          = duixs*gux[3]+duiys*guy[3]+duizs*guz[3];
    gkField[1][2]          = duixs*gux[4]+duiys*guy[4]+duizs*guz[4];

    gkField[2][0]          = pukxs*gux[2]+pukys*guy[2]+pukzs*guz[2];
    gkField[2][1]          = pukxs*gux[3]+pukys*guy[3]+pukzs*guz[3];
    gkField[2][2]          = pukxs*gux[4]+pukys*guy[4]+pukzs*guz[4];

    gkField[3][0]          = puixs*gux[2]+puiys*guy[2]+puizs*guz[2];
    gkField[3][1]          = puixs*gux[3]+puiys*guy[3]+puizs*guz[3];
    gkField[3][2]          = puixs*gux[4]+puiys*guy[4]+puizs*guz[4];

}

// Include versions of the kernels for N^2 calculations.

#define METHOD_NAME(a, b) a##N2##b
#include "kCalculateAmoebaCudaMutualInducedAndGkFields.h"
#define USE_OUTPUT_BUFFER_PER_WARP
#undef METHOD_NAME
#define METHOD_NAME(a, b) a##N2ByWarp##b
#include "kCalculateAmoebaCudaMutualInducedAndGkFields.h"

__global__ 
#if (__CUDA_ARCH__ >= 200)
__launch_bounds__(GF1XX_THREADS_PER_BLOCK, 1)
#elif (__CUDA_ARCH__ >= 120)
__launch_bounds__(GT2XX_THREADS_PER_BLOCK, 1)
#else
__launch_bounds__(G8X_THREADS_PER_BLOCK, 1)
#endif
void kInitializeMutualInducedAndGkField_kernel(
                   float* fixedEField,
                   float* fixedEFieldPolar,
                   float* fixedGkField,
                   float* polarizability,
                   float* inducedDipoleS,
                   float* inducedDipolePolarS )
{

    int pos = blockIdx.x*blockDim.x + threadIdx.x;
    while( pos < 3*cSim.atoms )
    {

        fixedEField[pos]          *= polarizability[pos];
        fixedEFieldPolar[pos]     *= polarizability[pos];
        fixedGkField[pos]         *= polarizability[pos];

        inducedDipoleS[pos]        = fixedEField[pos]       + fixedGkField[pos];
        inducedDipolePolarS[pos]   = fixedEFieldPolar[pos]  + fixedGkField[pos];
   
        pos                       += blockDim.x*gridDim.x;
    }

}

__global__ 
#if (__CUDA_ARCH__ >= 200)
__launch_bounds__(GF1XX_THREADS_PER_BLOCK, 1)
#elif (__CUDA_ARCH__ >= 120)
__launch_bounds__(GT2XX_THREADS_PER_BLOCK, 1)
#else
__launch_bounds__(G8X_THREADS_PER_BLOCK, 1)
#endif
void kReduceMutualInducedAndGkFieldDelta_kernel( float* arrayOfDeltas1, float* arrayOfDeltas2,
                                                 float* arrayOfDeltas3, float* arrayOfDeltas4, float* epsilon )
{
    extern __shared__ float4 delta[];

    delta[threadIdx.x].x    = 0.0f;
    delta[threadIdx.x].y    = 0.0f;
    delta[threadIdx.x].z    = 0.0f;
    delta[threadIdx.x].w    = 0.0f;

    unsigned int pos        = threadIdx.x;

    // load deltas

    while( pos < 3*cSim.atoms )
    {   
        delta[threadIdx.x].x  += arrayOfDeltas1[pos];
        delta[threadIdx.x].y  += arrayOfDeltas2[pos];
        delta[threadIdx.x].z  += arrayOfDeltas3[pos];
        delta[threadIdx.x].w  += arrayOfDeltas4[pos];
        pos                   += blockDim.x*gridDim.x;
    }   
    __syncthreads();

    // sum the deltas

    for (int offset = 1; offset < blockDim.x; offset *= 2 )
    {   
        if (threadIdx.x + offset < blockDim.x && (threadIdx.x & (2*offset-1)) == 0)
        {
            delta[threadIdx.x].x   += delta[threadIdx.x+offset].x;
            delta[threadIdx.x].y   += delta[threadIdx.x+offset].y;
            delta[threadIdx.x].z   += delta[threadIdx.x+offset].z;
            delta[threadIdx.x].w   += delta[threadIdx.x+offset].w;
        }
        __syncthreads();
    }   

    // set epsilons

    if (threadIdx.x == 0)
    {   
        epsilon[0]  = delta[0].x;
        epsilon[0]  = epsilon[0] < delta[0].y ? delta[0].y : epsilon[0];
        epsilon[0]  = epsilon[0] < delta[0].z ? delta[0].z : epsilon[0];
        epsilon[0]  = epsilon[0] < delta[0].w ? delta[0].w : epsilon[0];
        epsilon[0]  = 48.033324f*sqrtf( epsilon[0]/( (float) cSim.atoms ) );
    }   
}

/**

   matrixProduct/matrixProductP contains epsilon**2 on output

*/
__global__ 
#if (__CUDA_ARCH__ >= 200)
__launch_bounds__(GF1XX_THREADS_PER_BLOCK, 1)
#elif (__CUDA_ARCH__ >= 120)
__launch_bounds__(GT2XX_THREADS_PER_BLOCK, 1)
#else
__launch_bounds__(G8X_THREADS_PER_BLOCK, 1)
#endif
void kSorUpdateMutualInducedAndGkField_kernel(
                   float* polarizability,
                   float* inducedDipole, float* inducedDipoleP,
                   float* fixedEField,   float* fixedEFieldP,
                   float* matrixProduct, float* matrixProductP )
{

    float polarSOR = 0.55f;
    int pos        = blockIdx.x*blockDim.x + threadIdx.x;
    while( pos < 3*cSim.atoms )
    {

        float previousDipole           = inducedDipole[pos];
        float previousDipoleP          = inducedDipoleP[pos];
    
        inducedDipole[pos]             = fixedEField[pos]     + polarizability[pos]*matrixProduct[pos];
        inducedDipoleP[pos]            = fixedEFieldP[pos]    + polarizability[pos]*matrixProductP[pos];
    
        inducedDipole[pos]             = previousDipole   + polarSOR*( inducedDipole[pos]   - previousDipole  );   
        inducedDipoleP[pos]            = previousDipoleP  + polarSOR*( inducedDipoleP[pos]  - previousDipoleP );
    
        matrixProduct[pos]             = ( inducedDipole[pos]  - previousDipole  )*( inducedDipole[pos]  - previousDipole  );
        matrixProductP[pos]            = ( inducedDipoleP[pos] - previousDipoleP )*( inducedDipoleP[pos] - previousDipoleP );

        pos                           += blockDim.x*gridDim.x;
    }
}

__global__ 
#if (__CUDA_ARCH__ >= 200)
__launch_bounds__(GF1XX_THREADS_PER_BLOCK, 1)
#elif (__CUDA_ARCH__ >= 120)
__launch_bounds__(GT2XX_THREADS_PER_BLOCK, 1)
#else
__launch_bounds__(G8X_THREADS_PER_BLOCK, 1)
#endif
void kSorUpdateMutualInducedAndGkFieldS_kernel(
                   float* polarizability,
                   float* inducedDipole, float* inducedDipoleP,
                   float* fixedEField,   float* fixedEFieldP,
                   float* fixedGkField,
                   float* matrixProduct, float* matrixProductP )
{

    float polarSOR = 0.55f;
    int pos        = blockIdx.x*blockDim.x + threadIdx.x;
    while( pos < 3*cSim.atoms )
    {
        float previousDipole      = inducedDipole[pos];
        float previousDipoleP     = inducedDipoleP[pos];
    
        inducedDipole[pos]        = fixedGkField[pos]    + fixedEField[pos]     + polarizability[pos]*matrixProduct[pos];
        inducedDipoleP[pos]       = fixedGkField[pos]    + fixedEFieldP[pos]    + polarizability[pos]*matrixProductP[pos];
    
        inducedDipole[pos]        = previousDipole   + polarSOR*( inducedDipole[pos]   - previousDipole  );   
        inducedDipoleP[pos]       = previousDipoleP  + polarSOR*( inducedDipoleP[pos]  - previousDipoleP );
    
        matrixProduct[pos]        = ( inducedDipole[pos]  - previousDipole  )*( inducedDipole[pos]  - previousDipole  );
        matrixProductP[pos]       = ( inducedDipoleP[pos] - previousDipoleP )*( inducedDipoleP[pos] - previousDipoleP );
    
        pos                      += blockDim.x*gridDim.x;
    }
}

// reduce psWorkArray_3_1 -> outputArray
// reduce psWorkArray_3_2 -> outputPolarArray
// reduce psWorkArray_3_3 -> outputArrayS
// reduce psWorkArray_3_4 -> outputPolarArrayS

static void kReduceMutualInducedAndGkFields(amoebaGpuContext amoebaGpu,
                                            CUDAStream<float>* outputArray,  CUDAStream<float>* outputPolarArray,
                                            CUDAStream<float>* outputArrayS, CUDAStream<float>* outputPolarArrayS )
{
    gpuContext gpu = amoebaGpu->gpuContext;
    kReduceFields_kernel<<<gpu->sim.nonbond_blocks, gpu->sim.bsf_reduce_threads_per_block>>>(
                               gpu->sim.paddedNumberOfAtoms*3, gpu->sim.outputBuffers,
                               amoebaGpu->psWorkArray_3_1->_pDevData, outputArray->_pDevData, 0 );
    LAUNCHERROR("kReduceMutualInducedAndGkFields1");

    kReduceFields_kernel<<<gpu->sim.nonbond_blocks, gpu->sim.bsf_reduce_threads_per_block>>>(
                               gpu->sim.paddedNumberOfAtoms*3, gpu->sim.outputBuffers,
                               amoebaGpu->psWorkArray_3_2->_pDevData, outputPolarArray->_pDevData, 0 );
    LAUNCHERROR("kReduceMutualInducedAndGkFields2");

    kReduceFields_kernel<<<gpu->sim.nonbond_blocks, gpu->sim.bsf_reduce_threads_per_block>>>(
                               gpu->sim.paddedNumberOfAtoms*3, gpu->sim.outputBuffers,
                               amoebaGpu->psWorkArray_3_3->_pDevData, outputArrayS->_pDevData, 0 );
    LAUNCHERROR("kReduceMutualInducedAndGkFields3");

    kReduceFields_kernel<<<gpu->sim.nonbond_blocks, gpu->sim.bsf_reduce_threads_per_block>>>(
                               gpu->sim.paddedNumberOfAtoms*3, gpu->sim.outputBuffers,
                               amoebaGpu->psWorkArray_3_4->_pDevData, outputPolarArrayS->_pDevData, 0 );
    LAUNCHERROR("kReduceMutualInducedAndGkFields4");
}

/**---------------------------------------------------------------------------------------

   Compute mutual induce field

   @param amoebaGpu        amoebaGpu context

   --------------------------------------------------------------------------------------- */

static void cudaComputeAmoebaMutualInducedAndGkFieldMatrixMultiply( amoebaGpuContext amoebaGpu,
                                                                    CUDAStream<float>* outputArray,  CUDAStream<float>* outputPolarArray,
                                                                    CUDAStream<float>* outputArrayS, CUDAStream<float>* outputPolarArrayS )
{
  
   // ---------------------------------------------------------------------------------------

    static unsigned int threadsPerBlock = 0;

   // ---------------------------------------------------------------------------------------

    gpuContext gpu    = amoebaGpu->gpuContext;

    // clear output arrays

    kClearFields_3( amoebaGpu, 4 );

    // set threads/block first time through

    if( threadsPerBlock == 0 ){
        unsigned int maxThreads;
        if (gpu->sm_version >= SM_20)
            maxThreads = 384;
        else if (gpu->sm_version >= SM_12)
            maxThreads = 128;
        else
            maxThreads = 64;
        threadsPerBlock = std::min(getThreadsPerBlock( amoebaGpu, sizeof(MutualInducedParticle), gpu->sharedMemoryPerBlock ), maxThreads);
    }
    
    if (gpu->bOutputBufferPerWarp){
        kCalculateAmoebaMutualInducedAndGkFieldsN2ByWarp_kernel<<<gpu->sim.nonbond_blocks, threadsPerBlock, sizeof(MutualInducedParticle)*threadsPerBlock>>>(
                                                                 gpu->psWorkUnit->_pDevData,
                                                                 amoebaGpu->psWorkArray_3_1->_pDevData,
                                                                 amoebaGpu->psWorkArray_3_2->_pDevData,
                                                                 amoebaGpu->psWorkArray_3_3->_pDevData,
                                                                 amoebaGpu->psWorkArray_3_4->_pDevData );
    } else {

        kCalculateAmoebaMutualInducedAndGkFieldsN2_kernel<<<gpu->sim.nonbond_blocks, threadsPerBlock, sizeof(MutualInducedParticle)*threadsPerBlock>>>(
                                                            gpu->psWorkUnit->_pDevData,
                                                            amoebaGpu->psWorkArray_3_1->_pDevData,
                                                            amoebaGpu->psWorkArray_3_2->_pDevData,
                                                            amoebaGpu->psWorkArray_3_3->_pDevData,
                                                            amoebaGpu->psWorkArray_3_4->_pDevData );

    }
    LAUNCHERROR("kCalculateAmoebaMutualInducedAndGkFields");  

    kReduceMutualInducedAndGkFields( amoebaGpu, outputArray, outputPolarArray, outputArrayS, outputPolarArrayS );

}

/**---------------------------------------------------------------------------------------

   Compute mutual induce field

   @param amoebaGpu        amoebaGpu context

   --------------------------------------------------------------------------------------- */

static void cudaComputeAmoebaMutualInducedAndGkFieldBySOR( amoebaGpuContext amoebaGpu )
{
  
   // ---------------------------------------------------------------------------------------

    int done;
    int iteration;
    static int timestep = 0;
    timestep++;

    gpuContext gpu     = amoebaGpu->gpuContext;

   // ---------------------------------------------------------------------------------------

    // set  E_Field & E_FieldPolar] to [ E_Field & E_FieldPolar]*Polarizability
    // initialize [ InducedDipole & InducedDipolePolar ] to [ E_Field & E_FieldPolar]*Polarizability

    kInitializeMutualInducedAndGkField_kernel<<< gpu->sim.nonbond_blocks, gpu->sim.bsf_reduce_threads_per_block >>>(
         amoebaGpu->psE_Field->_pDevData,
         amoebaGpu->psE_FieldPolar->_pDevData,
         amoebaGpu->psGk_Field->_pDevData,
         amoebaGpu->psPolarizability->_pDevData,
         amoebaGpu->psInducedDipoleS->_pDevData,
         amoebaGpu->psInducedDipolePolarS->_pDevData );
    LAUNCHERROR("kInitializeMutualInducedAndGkField");  

    cudaMemcpy( amoebaGpu->psInducedDipole->_pDevData,        amoebaGpu->psE_Field->_pDevData,       3*gpu->sim.paddedNumberOfAtoms*sizeof( float ), cudaMemcpyDeviceToDevice );
    cudaMemcpy( amoebaGpu->psInducedDipolePolar->_pDevData,   amoebaGpu->psE_FieldPolar->_pDevData,  3*gpu->sim.paddedNumberOfAtoms*sizeof( float ), cudaMemcpyDeviceToDevice );

    // if polarization type is direct, set flags signalling done and return

    if( amoebaGpu->amoebaSim.polarizationType )
    {   
        amoebaGpu->mutualInducedDone          = 1;
        amoebaGpu->mutualInducedConverged     = 1;
        return;
    }   

    // ---------------------------------------------------------------------------------------
 
    done      = 0;
    iteration = 1;

    while( !done ){

        // matrix multiply

        cudaComputeAmoebaMutualInducedAndGkFieldMatrixMultiply( amoebaGpu,
                                                                amoebaGpu->psWorkVector[0],  amoebaGpu->psWorkVector[1],
                                                                amoebaGpu->psWorkVector[2],  amoebaGpu->psWorkVector[3] );

        LAUNCHERROR("cudaComputeAmoebaMutualInducedAndGkFieldMatrixMultiply");  

    // ---------------------------------------------------------------------------------------

        // post matrix multiply

        kSorUpdateMutualInducedAndGkField_kernel<<< gpu->sim.nonbond_blocks, gpu->sim.bsf_reduce_threads_per_block >>>(
           amoebaGpu->psPolarizability->_pDevData,
           amoebaGpu->psInducedDipole->_pDevData,     amoebaGpu->psInducedDipolePolar->_pDevData,
           amoebaGpu->psE_Field->_pDevData,           amoebaGpu->psE_FieldPolar->_pDevData,
           amoebaGpu->psWorkVector[0]->_pDevData,     amoebaGpu->psWorkVector[1]->_pDevData );
        LAUNCHERROR("cudaComputeAmoebaMutualInducedAndGkFieldSorUpdate1");  

        kSorUpdateMutualInducedAndGkFieldS_kernel<<< gpu->sim.nonbond_blocks, gpu->sim.bsf_reduce_threads_per_block >>>(
           amoebaGpu->psPolarizability->_pDevData,
           amoebaGpu->psInducedDipoleS->_pDevData,    amoebaGpu->psInducedDipolePolarS->_pDevData,
           amoebaGpu->psE_Field->_pDevData,          amoebaGpu->psE_FieldPolar->_pDevData,
           amoebaGpu->psGk_Field->_pDevData,
           amoebaGpu->psWorkVector[2]->_pDevData,     amoebaGpu->psWorkVector[3]->_pDevData );
        LAUNCHERROR("cudaComputeAmoebaMutualInducedAndGkFieldSorUpdate2");  

        // get total epsilon -- performing sums on gpu

        kReduceMutualInducedAndGkFieldDelta_kernel<<<1, amoebaGpu->epsilonThreadsPerBlock, 4*sizeof(float)*amoebaGpu->epsilonThreadsPerBlock>>>(
           amoebaGpu->psWorkVector[0]->_pDevData, amoebaGpu->psWorkVector[1]->_pDevData,
           amoebaGpu->psWorkVector[2]->_pDevData, amoebaGpu->psWorkVector[3]->_pDevData,
           amoebaGpu->psCurrentEpsilon->_pDevData );
        LAUNCHERROR("kReduceMutualInducedAndGkFieldDelta_kernel");

        // Debye=48.033324f

        amoebaGpu->psCurrentEpsilon->Download();
        float currentEpsilon                     = amoebaGpu->psCurrentEpsilon->_pSysData[0];
        amoebaGpu->mutualInducedCurrentEpsilon   = currentEpsilon;
 
        // check for nans

        if( currentEpsilon != currentEpsilon ){
             std::stringstream msg;            
             msg << "GkFieldBySOR: Nans detected in induced dipole calculation at iteration=" << iteration << " call=" << timestep << ".";
             throw OpenMM::OpenMMException( msg.str() );
        }

        // converged?

        if( iteration > amoebaGpu->mutualInducedMaxIterations || amoebaGpu->mutualInducedCurrentEpsilon < amoebaGpu->mutualInducedTargetEpsilon ){ 
            done = 1;
        }

        iteration++;
    }

    amoebaGpu->mutualInducedDone             = done;
    amoebaGpu->mutualInducedConverged        = ( !done || iteration > amoebaGpu->mutualInducedMaxIterations ) ? 0 : 1;
}

void cudaComputeAmoebaMutualInducedAndGkField( amoebaGpuContext amoebaGpu )
{
    if( amoebaGpu->mutualInducedIterativeMethod == 0 ){
        cudaComputeAmoebaMutualInducedAndGkFieldBySOR( amoebaGpu );
    }
}
