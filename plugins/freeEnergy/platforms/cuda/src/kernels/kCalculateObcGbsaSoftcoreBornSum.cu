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

#include "gputypes.h"
#include "freeEnergyGpuTypes.h"
#include "GpuFreeEnergyCudaKernels.h"
#include "kernels/cudaKernels.h"
#include "openmm/OpenMMException.h"

#include <cuda.h>

#define PARAMETER_PRINT 0
#define MAX_PARAMETER_PRINT 10

struct Atom {
    float x;
    float y;
    float z;
    float r;
    float sr;
    float sum;
    float polarScaleData;
};

static __constant__ cudaGmxSimulation cSim;
static __constant__ cudaFreeEnergyGmxSimulation gbsaSimDev;

extern "C" void SetCalculateObcGbsaSoftcoreBornSumSim( freeEnergyGpuContext freeEnergyGpu)
{
    cudaError_t status;
    status = cudaMemcpyToSymbol( cSim, &freeEnergyGpu->gpuContext->sim, sizeof(cudaGmxSimulation));
    RTERROR(status, "cudaMemcpyToSymbol: SetCalculateObcGbsaSoftcoreBornSumSim copy to cSim failed.");

    status = cudaMemcpyToSymbol( gbsaSimDev, &freeEnergyGpu->freeEnergySim, sizeof(cudaFreeEnergyGmxSimulation));
    RTERROR(status, "cudaMemcpyToSymbol: SetCalculateObcGbsaSoftcoreBornSumSim copy to gbsaSimDev failed.");
}

__global__ void kClearObcGbsaSoftcoreBornSum_kernel()
{
    unsigned int pos = blockIdx.x * blockDim.x + threadIdx.x;
    while (pos < cSim.stride * cSim.nonbondOutputBuffers)
    {
        ((float*)cSim.pBornSum)[pos]   = 0.0f;
        pos                           += gridDim.x * blockDim.x;
    }
}

__global__ void kClearSoftcoreBornForces_kernel()
{
    unsigned int pos = blockIdx.x * blockDim.x + threadIdx.x;
    while (pos < cSim.stride * cSim.nonbondOutputBuffers)
    {   
        ((float*)cSim.pBornForce)[pos]  = 0.0f;
        pos                            += gridDim.x * blockDim.x;
    }   
}

void kClearSoftcoreBornForces(gpuContext gpu)
{
  //  printf("kClearSoftcoreBornForces\n");
    kClearSoftcoreBornForces_kernel<<<gpu->sim.blocks, 384>>>();
    LAUNCHERROR("kClearSoftcoreBornForces");
}

void kClearObcGbsaSoftcoreBornSum(gpuContext gpu)
{
  //  printf("kClearObcGbsaBornSum\n");
    kClearObcGbsaSoftcoreBornSum_kernel<<<gpu->sim.blocks, 384>>>();
}

__global__ void kReduceObcGbsaSoftcoreBornForces_kernel()
{
    unsigned int pos = (blockIdx.x * blockDim.x + threadIdx.x);
    float energy     = 0.0f;

    while (pos < cSim.atoms)
    {
        float bornRadius         = cSim.pBornRadii[pos];
        float obcChain           = cSim.pObcChain[pos];
        float2 obcData           = cSim.pObcData[pos];
        float  nonPolarScaleData = gbsaSimDev.pNonPolarScalingFactors[pos];
        float totalForce         = 0.0f;
        float* pFt               = cSim.pBornForce + pos;

        int i                    = cSim.nonbondOutputBuffers;
        while (i >= 4)
        {
            float f1    = *pFt;
            pFt        += cSim.stride;
            float f2    = *pFt;
            pFt        += cSim.stride;
            float f3    = *pFt;
            pFt        += cSim.stride;
            float f4    = *pFt;
            pFt        += cSim.stride;
            totalForce += f1 + f2 + f3 + f4;
            i -= 4;
        }
        if (i >= 2)
        {
            float f1    = *pFt;
            pFt        += cSim.stride;
            float f2    = *pFt;
            pFt        += cSim.stride;
            totalForce += f1 + f2;
            i -= 2;
        }
        if (i > 0)
        {
            totalForce += *pFt;
        }
        float r            = (obcData.x + cSim.dielectricOffset + cSim.probeRadius);
        float ratio6       = pow((obcData.x + cSim.dielectricOffset) / bornRadius, 6.0f);
        float saTerm       = nonPolarScaleData*cSim.surfaceAreaFactor * r * r * ratio6;
        totalForce        += saTerm / bornRadius; // 1.102 == Temp mysterious fudge factor, FIX FIX FIX

        energy            += saTerm;

        totalForce        *= bornRadius * bornRadius * obcChain;

        pFt                = cSim.pBornForce + pos;
        *pFt               = totalForce;
        pos               += gridDim.x * blockDim.x;
    }

    // correct for surface area factor of -6
    cSim.pEnergy[blockIdx.x * blockDim.x + threadIdx.x] += energy / -6.0f;
}


void kReduceObcGbsaSoftcoreBornForces( gpuContext gpu ){

    kReduceObcGbsaSoftcoreBornForces_kernel<<<gpu->sim.blocks, gpu->sim.bsf_reduce_threads_per_block>>>();
    LAUNCHERROR("kReduceObcGbsaSoftcoreBornForces");
}


// Include versions of the kernels for N^2 calculations.

#define METHOD_NAME(a, b) a##N2##b
#include "kCalculateObcGbsaSoftcoreBornSum.h"
#define USE_OUTPUT_BUFFER_PER_WARP
#undef METHOD_NAME
#define METHOD_NAME(a, b) a##N2ByWarp##b
#include "kCalculateObcGbsaSoftcoreBornSum.h"

// Include versions of the kernels with cutoffs.

#undef METHOD_NAME
#undef USE_OUTPUT_BUFFER_PER_WARP
#define USE_CUTOFF
#define METHOD_NAME(a, b) a##Cutoff##b
#include "kCalculateObcGbsaSoftcoreBornSum.h"
#define USE_OUTPUT_BUFFER_PER_WARP
#undef METHOD_NAME
#define METHOD_NAME(a, b) a##CutoffByWarp##b
#include "kCalculateObcGbsaSoftcoreBornSum.h"

// Include versions of the kernels with periodic boundary conditions.

#undef METHOD_NAME
#undef USE_OUTPUT_BUFFER_PER_WARP
#define USE_PERIODIC
#define METHOD_NAME(a, b) a##Periodic##b
#include "kCalculateObcGbsaSoftcoreBornSum.h"
#define USE_OUTPUT_BUFFER_PER_WARP
#undef METHOD_NAME
#define METHOD_NAME(a, b) a##PeriodicByWarp##b
#include "kCalculateObcGbsaSoftcoreBornSum.h"

__global__ void kReduceObcGbsaSoftcoreBornSum_kernel()
{
    unsigned int pos = (blockIdx.x * blockDim.x + threadIdx.x);
    
    while (pos < cSim.atoms)
    {
        float sum   = 0.0f;
        float* pSt  = cSim.pBornSum + pos;
        float2 atom = cSim.pObcData[pos];
        
        // Get summed Born data
        for( int i = 0; i < cSim.nonbondOutputBuffers; i++ ){
            sum += *pSt;
            pSt += cSim.stride;
        }
        
        
        // Now calculate Born radius and OBC term.
        sum                    *= 0.5f * atom.x;
        float sum2              = sum * sum;
        float sum3              = sum * sum2;
        float tanhSum           = tanh(cSim.alphaOBC * sum - cSim.betaOBC * sum2 + cSim.gammaOBC * sum3);
        float nonOffsetRadii    = atom.x + cSim.dielectricOffset;
        float bornRadius        = 1.0f / (1.0f / atom.x - tanhSum / nonOffsetRadii); 
        float obcChain          = atom.x * (cSim.alphaOBC - 2.0f * cSim.betaOBC * sum + 3.0f * cSim.gammaOBC * sum2);
        obcChain                = (1.0f - tanhSum * tanhSum) * obcChain / nonOffsetRadii;        
        cSim.pBornRadii[pos]    = bornRadius;
        cSim.pObcChain[pos]     = obcChain;
        pos                    += gridDim.x * blockDim.x;
    }   
}

void kReduceObcGbsaSoftcoreBornSum(gpuContext gpu)
{
//    printf("kReduceObcGbsaSoftcoreBornSum\n");
    kReduceObcGbsaSoftcoreBornSum_kernel<<<gpu->sim.blocks, 384>>>();
    gpu->bRecalculateBornRadii = false;
    LAUNCHERROR("kReduceObcGbsaSoftcoreBornSum");
}

/** 
 * Initialize parameters for Cuda Obc softcore
 * 
 * @param freeEnergyGpu        freeEnergyGpu context
 * @param innerDielectric      solute dielectric
 * @param solventDielectric    solvent dielectric
 * @param radius               intrinsic Born radii
 * @param scale                Obc scaling factors
 * @param charge               atomic charges (possibly overwritten by other methods?)
 * @param nonPolarScalingFactors non-polar scaling factors
 *
 */

extern "C"
void  gpuSetObcSoftcoreParameters( freeEnergyGpuContext freeEnergyGpu, float innerDielectric, float solventDielectric, float nonPolarPrefactor,
                                   const std::vector<float>& radius, const std::vector<float>& scale,
                                   const std::vector<float>& charge, const std::vector<float>& nonPolarScalingFactors)
{

// ---------------------------------------------------------------------------------------

   static const float dielectricOffset    =    0.009f;
   static const float electricConstant    = -166.02691f;
   static const std::string methodName    = "gpuSetObcSoftcoreParameters";

// ---------------------------------------------------------------------------------------

    unsigned int numberOfParticles                       = radius.size();
    gpuContext gpu                                       = freeEnergyGpu->gpuContext;

    // initialize parameters

    freeEnergyGpu->psNonPolarScalingFactors              = new CUDAStream<float>( gpu->sim.paddedNumberOfAtoms, 1, "ObcSoftcoreNonPolarScaling");
    freeEnergyGpu->freeEnergySim.pNonPolarScalingFactors = freeEnergyGpu->psNonPolarScalingFactors->_pDevData;

    gpu->sim.surfaceAreaFactor                           =  -6.0f*PI*4.0f*nonPolarPrefactor;
    gpu->sim.preFactor                                   = 2.0f*electricConstant*((1.0f/innerDielectric)-(1.0f/solventDielectric))*gpu->sim.forceConversionFactor;

    for( unsigned int ii = 0; ii < numberOfParticles; ii++ ){
        (*gpu->psObcData)[ii].x                        = radius[ii] - dielectricOffset;
        (*gpu->psObcData)[ii].y                        = scale[ii] * (*gpu->psObcData)[ii].x;
        (*gpu->psPosq4)[ii].w                          = charge[ii];
        (*gpu->psBornRadii)[ii]                        = 0.0f;
        (*freeEnergyGpu->psNonPolarScalingFactors)[ii] = nonPolarScalingFactors[ii];
    }

    // diagnostics

    if( freeEnergyGpu->log ){
        (void) fprintf( freeEnergyGpu->log, "%s %u %u\n", methodName.c_str(), gpu->natoms, gpu->sim.paddedNumberOfAtoms );
        (void) fprintf( freeEnergyGpu->log, "surfaceAreaFactor=%15.7e preFactor=%15.7e\n", gpu->sim.surfaceAreaFactor, gpu->sim.preFactor);
#ifdef PARAMETER_PRINT
        int maxPrint = MAX_PARAMETER_PRINT;
        for( unsigned int ii = 0; ii < numberOfParticles; ii++ ){
            (void) fprintf( freeEnergyGpu->log, "%6u %13.6e %13.6e %8.3f %8.3f\n", ii, 
                            (*gpu->psObcData)[ii].x, (*gpu->psObcData)[ii].y, (*gpu->psPosq4)[ii].w, (*freeEnergyGpu->psNonPolarScalingFactors)[ii] );
             if( ii == maxPrint ){
                ii = numberOfParticles - maxPrint;
                if( ii < maxPrint )ii = maxPrint;
            }
        }
#endif
    }

    // dummy out extra atom data

    for (unsigned int ii = gpu->natoms; ii < gpu->sim.paddedNumberOfAtoms; ii++ ){
        (*gpu->psObcData)[ii].x                         = 0.01f;
        (*gpu->psObcData)[ii].y                         = 0.01f;
        (*freeEnergyGpu->psNonPolarScalingFactors)[ii]  = 0.0f;
        (*gpu->psBornRadii)[ii]                         = 0.0f;
    }

    // load data to board

    gpu->psObcData->Upload();
    gpu->psPosq4->Upload();
    gpu->psBornRadii->Upload();
    freeEnergyGpu->psNonPolarScalingFactors->Upload();

    return;
}

void kPrintObcGbsaSoftcore( freeEnergyGpuContext freeEnergyGpu, std::string callId, int call, FILE* log){

    gpuContext gpu = freeEnergyGpu->gpuContext;
    int maxPrint   = gpu->natoms;

    (void) fprintf( log, "kPrintObcGbsaSoftcore %s %d\n", callId.c_str(), call );

    gpu->psObcData->Download();
    gpu->psBornRadii->Download();
    gpu->psBornForce->Download();
    gpu->psPosq4->Download();
    freeEnergyGpu->psNonPolarScalingFactors->Download();

    CUDAStream<float4>* sigEps4          = freeEnergyGpu->psSigEps4;
    sigEps4->Download();

    (void) fprintf( log, "BornSum Born radii & params\n" );
    for( int ii = 0; ii < gpu->sim.paddedNumberOfAtoms; ii++ ){
        //(void) fprintf( log, "%6d prm[%15.7e %15.7e %15.7e %15.7e] [%15.7e %15.7e %15.7e %15.7e] bR=%15.7e bF=%15.7e swDrv=%3.1f x[%8.3f %8.3f %8.3f %15.7f]\n",
        (void) fprintf( log, "%6d prm[%15.7e %15.7e %15.7e] sig/eps4[%15.7e %15.7e %15.7e %15.7e] bR=%15.7e bF=%15.7e\n",
                        ii,

                        gpu->psObcData->_pSysData[ii].x,
                        gpu->psObcData->_pSysData[ii].y,
                        freeEnergyGpu->psNonPolarScalingFactors->_pSysData[ii],

                        sigEps4->_pSysData[ii].x,
                        sigEps4->_pSysData[ii].y,
                        sigEps4->_pSysData[ii].z,
                        sigEps4->_pSysData[ii].w,

                        gpu->psBornRadii->_pSysData[ii],
                        gpu->psBornForce->_pSysData[ii]
/*
                        gpu->psPosq4->_pSysData[ii].x,
                        gpu->psPosq4->_pSysData[ii].y,
                        gpu->psPosq4->_pSysData[ii].z,
                        gpu->psPosq4->_pSysData[ii].w );
*/
                        );

        if( (ii == maxPrint) && ( ii < (gpu->natoms - maxPrint)) ){
            ii = gpu->natoms - maxPrint;
        }
    }

}
extern __global__ void kFindBlockBoundsCutoff_kernel();
extern __global__ void kFindBlockBoundsPeriodic_kernel();

extern __global__ void kFindBlocksWithInteractionsCutoff_kernel();
extern __global__ void kFindBlocksWithInteractionsPeriodic_kernel();

extern __global__ void kFindInteractionsWithinBlocksCutoff_kernel(unsigned int*);
extern __global__ void kFindInteractionsWithinBlocksPeriodic_kernel(unsigned int*);

void kCalculateObcGbsaSoftcoreBornSum( freeEnergyGpuContext freeEnergyGpu )
{
  //  printf("kCalculateObcGbsaSoftcoreBornSum\n");
    gpuContext gpu = freeEnergyGpu->gpuContext;

    kClearObcGbsaSoftcoreBornSum(gpu);
    LAUNCHERROR("kClearBornSum from kCalculateObcGbsaSoftcoreBornSum");

    switch ( freeEnergyGpu->freeEnergySim.nonbondedMethod )
    {
        case FREE_ENERGY_NO_CUTOFF:

            if (gpu->bOutputBufferPerWarp)
                kCalculateObcGbsaSoftcoreN2ByWarpBornSum_kernel<<<gpu->sim.nonbond_blocks, gpu->sim.nonbond_threads_per_block,
                        sizeof(Atom)*gpu->sim.nonbond_threads_per_block>>>(gpu->sim.pWorkUnit);
            else
                kCalculateObcGbsaSoftcoreN2BornSum_kernel<<<gpu->sim.nonbond_blocks, gpu->sim.nonbond_threads_per_block,
                        sizeof(Atom)*gpu->sim.nonbond_threads_per_block>>>(gpu->sim.pWorkUnit);

            break;

        case FREE_ENERGY_CUTOFF:

            kFindBlockBoundsCutoff_kernel<<<(gpu->psGridBoundingBox->_length+63)/64, 64>>>();
            LAUNCHERROR("kFindBlockBoundsCutoff");
            kFindBlocksWithInteractionsCutoff_kernel<<<gpu->sim.interaction_blocks, gpu->sim.interaction_threads_per_block>>>();
            LAUNCHERROR("kFindBlocksWithInteractionsCutoff");
            compactStream(gpu->compactPlan, gpu->sim.pInteractingWorkUnit, gpu->sim.pWorkUnit, gpu->sim.pInteractionFlag, gpu->sim.workUnits, gpu->sim.pInteractionCount);
            kFindInteractionsWithinBlocksCutoff_kernel<<<gpu->sim.nonbond_blocks, gpu->sim.nonbond_threads_per_block,
                    sizeof(unsigned int)*gpu->sim.nonbond_threads_per_block>>>(gpu->sim.pInteractingWorkUnit);


            if (gpu->bOutputBufferPerWarp)
                kCalculateObcGbsaSoftcoreCutoffByWarpBornSum_kernel<<<gpu->sim.nonbond_blocks, gpu->sim.nonbond_threads_per_block,
                        (sizeof(Atom)+sizeof(float))*gpu->sim.nonbond_threads_per_block>>>(gpu->sim.pInteractingWorkUnit);
            else
                kCalculateObcGbsaSoftcoreCutoffBornSum_kernel<<<gpu->sim.nonbond_blocks, gpu->sim.nonbond_threads_per_block,
                        (sizeof(Atom)+sizeof(float))*gpu->sim.nonbond_threads_per_block>>>(gpu->sim.pInteractingWorkUnit);

            break;

        case FREE_ENERGY_PERIODIC:

            kFindBlockBoundsPeriodic_kernel<<<(gpu->psGridBoundingBox->_length+63)/64, 64>>>();
            LAUNCHERROR("kFindBlockBoundsPeriodic");
            kFindBlocksWithInteractionsPeriodic_kernel<<<gpu->sim.interaction_blocks, gpu->sim.interaction_threads_per_block>>>();
            LAUNCHERROR("kFindBlocksWithInteractionsPeriodic");
            compactStream(gpu->compactPlan, gpu->sim.pInteractingWorkUnit, gpu->sim.pWorkUnit, gpu->sim.pInteractionFlag, gpu->sim.workUnits, gpu->sim.pInteractionCount);
            kFindInteractionsWithinBlocksPeriodic_kernel<<<gpu->sim.nonbond_blocks, gpu->sim.nonbond_threads_per_block,
                    sizeof(unsigned int)*gpu->sim.nonbond_threads_per_block>>>(gpu->sim.pInteractingWorkUnit);

            if (gpu->bOutputBufferPerWarp)
                kCalculateObcGbsaSoftcorePeriodicByWarpBornSum_kernel<<<gpu->sim.nonbond_blocks, gpu->sim.nonbond_threads_per_block,
                        (sizeof(Atom)+sizeof(float))*gpu->sim.nonbond_threads_per_block>>>(gpu->sim.pInteractingWorkUnit);
            else
                kCalculateObcGbsaSoftcorePeriodicBornSum_kernel<<<gpu->sim.nonbond_blocks, gpu->sim.nonbond_threads_per_block,
                        (sizeof(Atom)+sizeof(float))*gpu->sim.nonbond_threads_per_block>>>(gpu->sim.pInteractingWorkUnit);


            break;

        default:
            throw OpenMM::OpenMMException( "Nonbonded softcore method not recognized." );

    }
    LAUNCHERROR("kCalculateObcGbsaSoftcoreBornSum");

}
