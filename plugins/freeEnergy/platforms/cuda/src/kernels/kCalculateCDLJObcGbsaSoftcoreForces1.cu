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

#include "kernels/gputypes.h"
#include "kernels/cudatypes.h"
#include "kernels/cudaKernels.h"
#include "GpuFreeEnergyCudaKernels.h"

#include <stdio.h>
#include <cuda.h>
#include <vector_functions.h>
#include <cstdlib>
#include <sstream>

#define USE_SOFTCORE_LJ

struct Atom {
    float x;
    float y;
    float z;
    float q;
    float sig;
    float eps;
    float br;
    float softCoreLJLambda;
    float fx;
    float fy;
    float fz;
    float fb;
};

static __constant__ cudaGmxSimulation cSim;
static __constant__ cudaFreeEnergyGmxSimulation feSimDev;

void SetCalculateCDLJObcGbsaSoftcoreGpu1Sim( freeEnergyGpuContext freeEnergyGpu ){

    cudaError_t status;
    status = cudaMemcpyToSymbol(cSim, &freeEnergyGpu->gpuContext->sim, sizeof(cudaGmxSimulation));
    RTERROR(status, "cudaMemcpyToSymbol: SetCalculateCDLJObcGbsaSoftcoreGpu1Sim copy to cSim failed");

    status = cudaMemcpyToSymbol( feSimDev, &freeEnergyGpu->freeEnergySim, sizeof(cudaFreeEnergyGmxSimulation));
    RTERROR(status, "cudaMemcpyToSymbol: SetCalculateCDLJObcGbsaSoftcoreGpu1Sim copy to feSimDev failed");
}

// Include versions of the kernel for N^2 calculations with softcore LJ.

#undef METHOD_NAME
#define METHOD_NAME(a, b) a##N2##b
#undef USE_OUTPUT_BUFFER_PER_WARP
#define USE_SOFTCORE_LJ
#include "kCalculateCDLJObcGbsaSoftcoreForces1.h"

#define USE_OUTPUT_BUFFER_PER_WARP
#undef METHOD_NAME
#define METHOD_NAME(a, b) a##N2ByWarp##b
#include "kCalculateCDLJObcGbsaSoftcoreForces1.h"
#undef USE_SOFTCORE_LJ
#undef USE_OUTPUT_BUFFER_PER_WARP

// Include versions of the kernel with cutoffs.

#undef METHOD_NAME
#undef USE_OUTPUT_BUFFER_PER_WARP
#define USE_CUTOFF
#define METHOD_NAME(a, b) a##Cutoff##b
#include "kCalculateCDLJObcGbsaSoftcoreForces1.h"
#define USE_OUTPUT_BUFFER_PER_WARP
#undef METHOD_NAME
#define METHOD_NAME(a, b) a##CutoffByWarp##b
#include "kCalculateCDLJObcGbsaSoftcoreForces1.h"

// Include versions of the kernel with periodic boundary conditions.

#undef METHOD_NAME
#undef USE_OUTPUT_BUFFER_PER_WARP
#define USE_PERIODIC
#define METHOD_NAME(a, b) a##Periodic##b
#include "kCalculateCDLJObcGbsaSoftcoreForces1.h"
#define USE_OUTPUT_BUFFER_PER_WARP
#undef METHOD_NAME
#define METHOD_NAME(a, b) a##PeriodicByWarp##b
#include "kCalculateCDLJObcGbsaSoftcoreForces1.h"

/**
 * 
 * Calculate Born radii and first GBSA loop forces/energy
 *
 * @param gpu     gpu context
 *
 */
void kCalculateCDLJObcGbsaSoftcoreForces1( freeEnergyGpuContext freeEnergyGpu )
{
//    printf("kCalculateCDLJObcGbsaForces1\n");
    gpuContext gpu = freeEnergyGpu->gpuContext;

//fprintf( stderr, "kCalculateCDLJObcGbsaSoftcoreForces1 cutoff=%15.7e blks=%u thread/.block=%u nbMethod==%d warp=%u\n",
//         gpu->sim.nonbondedCutoffSqr, gpu->sim.nonbond_blocks, gpu->sim.nonbond_threads_per_block, 
//         freeEnergyGpu->freeEnergySim.nonbondedMethod, gpu->bOutputBufferPerWarp);

//#define DEBUG
#ifdef DEBUG 
fprintf( stderr, "kCalculateCDLJObcGbsaSoftcoreForces1 cutoff=%15.7e\n", gpu->sim.nonbondedCutoffSqr );
int psize = gpu->sim.paddedNumberOfAtoms;
CUDAStream<float4>* pdE1 = new CUDAStream<float4>( psize, 1, "pdE");
CUDAStream<float4>* pdE2 = new CUDAStream<float4>( psize, 1, "pdE");
float bF,bR;
float bF1,b2;
float ratio;
float atomicRadii;
showWorkUnitsFreeEnergy( freeEnergyGpu, 1 );
for( int ii = 0; ii < psize; ii++ ){

pdE1->_pSysData[ii].x = 0.001f;
pdE1->_pSysData[ii].y = 0.001f;
pdE1->_pSysData[ii].z = 0.001f;
pdE1->_pSysData[ii].w = 0.001f;

pdE2->_pSysData[ii].x = 0.001f;
pdE2->_pSysData[ii].y = 0.001f;
pdE2->_pSysData[ii].z = 0.001f;
pdE2->_pSysData[ii].w = 0.001f;
}
pdE1->Upload();
pdE2->Upload();
#endif

    switch( freeEnergyGpu->freeEnergySim.nonbondedMethod )
    {
        case FREE_ENERGY_NO_CUTOFF:

            // use softcore LJ potential
#ifdef DEBUG
            if (gpu->bOutputBufferPerWarp)
                   kCalculateCDLJObcGbsaSoftcoreN2ByWarpForces1_kernel<<<gpu->sim.nonbond_blocks, gpu->sim.nonbond_threads_per_block,
                           sizeof(Atom)*gpu->sim.nonbond_threads_per_block>>>(gpu->sim.pWorkUnit,pdE1->_pDevData, pdE2->_pDevData);
            else
                   kCalculateCDLJObcGbsaSoftcoreN2Forces1_kernel<<<gpu->sim.nonbond_blocks, gpu->sim.nonbond_threads_per_block,
                           sizeof(Atom)*gpu->sim.nonbond_threads_per_block>>>(gpu->sim.pWorkUnit, pdE1->_pDevData, pdE2->_pDevData);
#else   
            if (gpu->bOutputBufferPerWarp)
                   kCalculateCDLJObcGbsaSoftcoreN2ByWarpForces1_kernel<<<gpu->sim.nonbond_blocks, gpu->sim.nonbond_threads_per_block,
                           sizeof(Atom)*gpu->sim.nonbond_threads_per_block>>>(gpu->sim.pWorkUnit );
            else
                   kCalculateCDLJObcGbsaSoftcoreN2Forces1_kernel<<<gpu->sim.nonbond_blocks, gpu->sim.nonbond_threads_per_block,
                           sizeof(Atom)*gpu->sim.nonbond_threads_per_block>>>(gpu->sim.pWorkUnit );
   
#endif
            LAUNCHERROR("kCalculateCDLJObcGbsaSoftcoreForces1");
            break;

        case FREE_ENERGY_CUTOFF:

#ifdef DEBUG

            if (gpu->bOutputBufferPerWarp)
                kCalculateCDLJObcGbsaSoftcoreCutoffByWarpForces1_kernel<<<gpu->sim.nonbond_blocks, gpu->sim.nonbond_threads_per_block,
                        (sizeof(Atom)+sizeof(float))*gpu->sim.nonbond_threads_per_block>>>(gpu->sim.pInteractingWorkUnit, pdE1->_pDevData, pdE2->_pDevData);
            else
                kCalculateCDLJObcGbsaSoftcoreCutoffForces1_kernel<<<gpu->sim.nonbond_blocks, gpu->sim.nonbond_threads_per_block,
                        (sizeof(Atom)+sizeof(float))*gpu->sim.nonbond_threads_per_block>>>(gpu->sim.pInteractingWorkUnit, pdE1->_pDevData, pdE2->_pDevData);
#else
            if (gpu->bOutputBufferPerWarp)
                kCalculateCDLJObcGbsaSoftcoreCutoffByWarpForces1_kernel<<<gpu->sim.nonbond_blocks, gpu->sim.nonbond_threads_per_block,
                        (sizeof(Atom)+sizeof(float))*gpu->sim.nonbond_threads_per_block>>>(gpu->sim.pInteractingWorkUnit );
            else
                kCalculateCDLJObcGbsaSoftcoreCutoffForces1_kernel<<<gpu->sim.nonbond_blocks, gpu->sim.nonbond_threads_per_block,
                        (sizeof(Atom)+sizeof(float))*gpu->sim.nonbond_threads_per_block>>>(gpu->sim.pInteractingWorkUnit );
#endif

            break;

        case FREE_ENERGY_PERIODIC:

            if (gpu->bOutputBufferPerWarp)
                kCalculateCDLJObcGbsaSoftcorePeriodicByWarpForces1_kernel<<<gpu->sim.nonbond_blocks, gpu->sim.nonbond_threads_per_block,
                        (sizeof(Atom)+sizeof(float))*gpu->sim.nonbond_threads_per_block>>>(gpu->sim.pInteractingWorkUnit);
            else
                kCalculateCDLJObcGbsaSoftcorePeriodicForces1_kernel<<<gpu->sim.nonbond_blocks, gpu->sim.nonbond_threads_per_block,
                        (sizeof(Atom)+sizeof(float))*gpu->sim.nonbond_threads_per_block>>>(gpu->sim.pInteractingWorkUnit);

            LAUNCHERROR("kCalculateCDLJObcGbsaSoftcorePeriodicForces1");
            break;
    }

#ifdef DEBUG 
pdE1->Download();
pdE2->Download();
gpu->psPosq4->Download();
gpu->psGBVIData->Download();
gpu->psBornRadii->Download();
freeEnergyGpu->psSwitchDerivative->Download();
fprintf( stderr, "PdeCud %d\n", TARGET );
bF = 0.0;
for( int ii = 0; ii < gpu->natoms; ii++ ){
bF += pdE1->_pSysData[ii].x;
if( fabs( pdE1->_pSysData[ii].w ) > 1.0e-03 ){
fprintf( stderr, "%4d %15.7e %15.7e %15.7e %15.7e    %15.7e %15.7e %15.7e %15.7e\n", ii, 
         pdE1->_pSysData[ii].x, pdE1->_pSysData[ii].y, pdE1->_pSysData[ii].z, pdE1->_pSysData[ii].w,
         pdE2->_pSysData[ii].x, pdE2->_pSysData[ii].y, pdE2->_pSysData[ii].z, pdE2->_pSysData[ii].w );
}
}
bR      = gpu->psBornRadii->_pSysData[TARGET];
atomicRadii = gpu->psGBVIData->_pSysData[TARGET].x; 
ratio   = (atomicRadii/bR);
bF1     = bF + (3.0f*gpu->psGBVIData->_pSysData[TARGET].z*ratio*ratio*ratio)/bR; 
b2      = bR*bR;
bF1     *= (1.0f/3.0f)*b2*b2;
fprintf( stderr, "sumbF Cud %6d %15.7e %15.7e %15.7e\n", TARGET, bF, bF1, bR);
#endif

}
