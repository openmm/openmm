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

#include "cudaKernels.h"
#include "amoebaCudaKernels.h"
#include "kCalculateAmoebaCudaUtilities.h"
#include "openmm/OpenMMException.h"

#include <stdio.h>
#include <cuda.h>
#include <cstdlib>
using namespace std; 

#define SQRT sqrtf

static __constant__ cudaGmxSimulation cSim;
static __constant__ cudaAmoebaGmxSimulation cAmoebaSim;
extern __global__ void kFindInteractionsWithinBlocksPeriodic_kernel(unsigned int*);

void SetCalculateAmoebaMultipolePotentialSim(amoebaGpuContext amoebaGpu)
{
    cudaError_t status;
    gpuContext gpu = amoebaGpu->gpuContext;
    status         = cudaMemcpyToSymbol(cSim, &gpu->sim, sizeof(cudaGmxSimulation));     
    RTERROR(status, "SetCalculateAmoebaMultipolePotentialSim: cudaMemcpyToSymbol: SetSim copy to cSim failed");
    status         = cudaMemcpyToSymbol(cAmoebaSim, &amoebaGpu->amoebaSim, sizeof(cudaAmoebaGmxSimulation));     
    RTERROR(status, "SetCalculateAmoebaMultipolePotentialSim: cudaMemcpyToSymbol: SetSim copy to cAmoebaSim failed");
}

void GetCalculateAmoebaMultipolePotentialSim(amoebaGpuContext amoebaGpu)
{
    cudaError_t status;
    gpuContext gpu = amoebaGpu->gpuContext;
    status = cudaMemcpyFromSymbol(&gpu->sim, cSim, sizeof(cudaGmxSimulation));     
    RTERROR(status, "GetCalculateAmoebaMultipolePotentialSim: cudaMemcpyFromSymbol: SetSim copy from cSim failed");
    status = cudaMemcpyFromSymbol(&amoebaGpu->amoebaSim, cAmoebaSim, sizeof(cudaAmoebaGmxSimulation));     
    RTERROR(status, "GetCalculateAmoebaMultipolePotentialSim: cudaMemcpyFromSymbol: SetSim copy from cAmoebaSim failed");
}

struct ElectrostaticPotentialParticle {

    // coordinates charge

    float x;
    float y;
    float z;
    float q;

    // lab frame dipole

    float labFrameDipole[3];

    // lab frame quadrupole

    float labFrameQuadrupole[9];

    // induced dipole

    float inducedDipole[3];

};

/**---------------------------------------------------------------------------------------

   Load data for particle w/ index=atomI

   @param sa        address to store atomI's coordinates and multipole moments
   @param atomI     index of atom whose data is to be stored

   --------------------------------------------------------------------------------------- */

static __device__ void loadElectrostaticPotentialParticle( volatile struct ElectrostaticPotentialParticle* sA, unsigned int atomI ){

    // coordinates & charge

    sA->x                        = cSim.pPosq[atomI].x;
    sA->y                        = cSim.pPosq[atomI].y;
    sA->z                        = cSim.pPosq[atomI].z;
    sA->q                        = cSim.pPosq[atomI].w;

    // lab dipole

    sA->labFrameDipole[0]        = cAmoebaSim.pLabFrameDipole[atomI*3];
    sA->labFrameDipole[1]        = cAmoebaSim.pLabFrameDipole[atomI*3+1];
    sA->labFrameDipole[2]        = cAmoebaSim.pLabFrameDipole[atomI*3+2];

    // lab quadrupole

    sA->labFrameQuadrupole[0]    = cAmoebaSim.pLabFrameQuadrupole[atomI*9];
    sA->labFrameQuadrupole[1]    = cAmoebaSim.pLabFrameQuadrupole[atomI*9+1];
    sA->labFrameQuadrupole[2]    = cAmoebaSim.pLabFrameQuadrupole[atomI*9+2];
    sA->labFrameQuadrupole[3]    = cAmoebaSim.pLabFrameQuadrupole[atomI*9+3];
    sA->labFrameQuadrupole[4]    = cAmoebaSim.pLabFrameQuadrupole[atomI*9+4];
    sA->labFrameQuadrupole[5]    = cAmoebaSim.pLabFrameQuadrupole[atomI*9+5];
    sA->labFrameQuadrupole[6]    = cAmoebaSim.pLabFrameQuadrupole[atomI*9+6];
    sA->labFrameQuadrupole[7]    = cAmoebaSim.pLabFrameQuadrupole[atomI*9+7];
    sA->labFrameQuadrupole[8]    = cAmoebaSim.pLabFrameQuadrupole[atomI*9+8];

    // induced dipole

    sA->inducedDipole[0]         = cAmoebaSim.pInducedDipole[atomI*3];
    sA->inducedDipole[1]         = cAmoebaSim.pInducedDipole[atomI*3+1];
    sA->inducedDipole[2]         = cAmoebaSim.pInducedDipole[atomI*3+2];

}

/**---------------------------------------------------------------------------------------

   Calculate potential at grid point due atomI
   Code adapted from TINKER routine potpoint in potpoint.f

   @param atomI     atomI's coordinates and multipole moments
   @param gridPoint grid coordinates
   @param potential output potential

   --------------------------------------------------------------------------------------- */

__device__ void calculateElectrostaticPotentialForAtomGridPoint_kernel( volatile ElectrostaticPotentialParticle& atomI, volatile float4& gridPoint, float* potential ){
  
    float xr                 = atomI.x - gridPoint.x;
    float yr                 = atomI.y - gridPoint.y;
    float zr                 = atomI.z - gridPoint.z;
   
    float r2                 = xr*xr + yr*yr + zr*zr;
    float r                  = sqrtf( r2 );

    float rr1                = 1.0f/r;
    *potential               = atomI.q*rr1;
    float rr2                = rr1*rr1;
    float rr3                = rr1*rr2;

    float scd                = atomI.labFrameDipole[0]*xr     +  atomI.labFrameDipole[1]*yr    + atomI.labFrameDipole[2]*zr;
    float scu                =  atomI.inducedDipole[0]*xr     +   atomI.inducedDipole[1]*yr    +  atomI.inducedDipole[2]*zr;
    *potential              -= (scd + scu)*rr3;

    float rr5                = 3.0f*rr3*rr2;
    float scq                = xr*(atomI.labFrameQuadrupole[0]*xr + atomI.labFrameQuadrupole[1]*yr + atomI.labFrameQuadrupole[2]*zr);
          scq               += yr*(atomI.labFrameQuadrupole[1]*xr + atomI.labFrameQuadrupole[4]*yr + atomI.labFrameQuadrupole[5]*zr);
          scq               += zr*(atomI.labFrameQuadrupole[2]*xr + atomI.labFrameQuadrupole[5]*yr + atomI.labFrameQuadrupole[8]*zr);
    *potential              += scq*rr5;

    return;

}

// Include versions of the kernels for N x PotentialGridSize calculations.

#undef USE_OUTPUT_BUFFER_PER_WARP
#undef METHOD_NAME 
#define METHOD_NAME(a, b) a##NxG##b
#include "kCalculateAmoebaCudaElectrostaticPotential.h"

#define USE_OUTPUT_BUFFER_PER_WARP
#undef METHOD_NAME
#define METHOD_NAME(a, b) a##NxGByWarp##b
#include "kCalculateAmoebaCudaElectrostaticPotential.h"

// Kernel to reduce potential

__global__ 
#if (__CUDA_ARCH__ >= 200)
__launch_bounds__(GF1XX_THREADS_PER_BLOCK, 1)
#elif (__CUDA_ARCH__ >= 120)
__launch_bounds__(GT2XX_THREADS_PER_BLOCK, 1)
#else
__launch_bounds__(G8X_THREADS_PER_BLOCK, 1)
#endif
void kReducePotential_kernel()
{
    unsigned int pos             = (blockIdx.x * blockDim.x + threadIdx.x);
    float conversionFactor       = (cAmoebaSim.electric/cAmoebaSim.dielec);
   
    // Reduce potential
    while (pos < cAmoebaSim.paddedPotentialGridSize)
    {
        float totalPotential         = 0.0f;
        float* pFt                   = cAmoebaSim.pPotential + pos;
        int i                        = cSim.outputBuffers;
        while (i >= 4)
        {
            float f1             = *pFt;
            pFt                 += cAmoebaSim.paddedPotentialGridSize;
            float f2             = *pFt;
            pFt                 += cAmoebaSim.paddedPotentialGridSize;
            float f3             = *pFt;
            pFt                 += cAmoebaSim.paddedPotentialGridSize;
            float f4             = *pFt;
            pFt                 += cAmoebaSim.paddedPotentialGridSize;
            totalPotential      += f1 + f2 + f3 + f4;
            i                   -= 4;
        }
        if (i >= 2)
        {
            float f1             = *pFt;
            pFt                 += cAmoebaSim.paddedPotentialGridSize;
            float f2             = *pFt;
            pFt                 += cAmoebaSim.paddedPotentialGridSize;
            totalPotential      += f1 + f2;
            i                   -= 2;
        }
        if (i > 0)
        {
            totalPotential += *pFt;
        }
        totalPotential *= conversionFactor;
        pFt             = cAmoebaSim.pPotential + pos;
        *pFt            = totalPotential;
        pos            += gridDim.x*blockDim.x;
    }   
}

/**---------------------------------------------------------------------------------------

   Reduce Amoeba electrostatic potential

   @param gpu        gpu context

   --------------------------------------------------------------------------------------- */

void kReducePotential(gpuContext gpu)
{
    kReducePotential_kernel<<<gpu->sim.blocks, gpu->sim.bsf_reduce_threads_per_block>>>();
    LAUNCHERROR("kReducePotential");
}

/**---------------------------------------------------------------------------------------

   Compute Amoeba electrostatic potential

   @param amoebaGpu        amoebaGpu context

   --------------------------------------------------------------------------------------- */

void cudaComputeAmoebaElectrostaticPotential( amoebaGpuContext amoebaGpu ){
  
   // ---------------------------------------------------------------------------------------

    gpuContext gpu = amoebaGpu->gpuContext;

    // on first pass, set threads/block

    static unsigned int threadsPerBlock = 0;
    if( threadsPerBlock == 0 ){
        unsigned int maxThreads;
        if (gpu->sm_version >= SM_20)
            //maxThreads = 384;
            maxThreads = 512;
        else if (gpu->sm_version >= SM_12)
            maxThreads = 128;
        else
            maxThreads = 64;
        threadsPerBlock = std::min(getThreadsPerBlock(amoebaGpu, sizeof(ElectrostaticPotentialParticle), gpu->sharedMemoryPerBlock), maxThreads);
    }

    if (gpu->bOutputBufferPerWarp){
        kCalculateAmoebaCudaElectrostaticPotentialNxGByWarp_kernel<<<gpu->sim.nonbond_blocks, threadsPerBlock, sizeof(ElectrostaticPotentialParticle)*threadsPerBlock>>>( );
    } else {
        kCalculateAmoebaCudaElectrostaticPotentialNxG_kernel<<<gpu->sim.nonbond_blocks, threadsPerBlock, sizeof(ElectrostaticPotentialParticle)*threadsPerBlock>>>( );
    }
    LAUNCHERROR("kCalculateAmoebaCudaElectrostaticPotential");

    kReducePotential( amoebaGpu->gpuContext );

   // ---------------------------------------------------------------------------------------
}

void kCalculateAmoebaMultipolePotential(amoebaGpuContext amoebaGpu ) 
{

    // setup

    kSetupAmoebaMultipoleForces(amoebaGpu, false ); 

    // calculate electrostatic potential

    cudaComputeAmoebaElectrostaticPotential( amoebaGpu );

}
