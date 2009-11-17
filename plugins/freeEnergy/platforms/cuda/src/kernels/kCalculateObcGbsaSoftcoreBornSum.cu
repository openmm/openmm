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
#include "kernels/cudaKernels.h"
#include "GpuObcGbsaSoftcore.h"

#include <cuda.h>
#include <string>

struct Atom {
    float x;
    float y;
    float z;
    float r;
    float sr;
    float sum;
    float polarScaleData;
};

struct cudaFreeEnergySimulationObcGbsaSoftcore {
    float* pNonPolarScalingFactors;
};
struct cudaFreeEnergySimulationObcGbsaSoftcore gbsaSim;

static __constant__ cudaGmxSimulation cSim;
static __constant__ cudaFreeEnergySimulationObcGbsaSoftcore gbsaSimDev;

extern "C"
void SetCalculateObcGbsaSoftcoreBornSumSim(gpuContext gpu)
{
    cudaError_t status;
    status = cudaMemcpyToSymbol(cSim, &gpu->sim, sizeof(cudaGmxSimulation));     
    RTERROR(status, "SetCalculateObcGbsaSoftcoreBornSumSim: cudaMemcpyToSymbol: SetSim copy to cSim failed");
}

extern "C"
void SetCalculateObcGbsaSoftcoreNonPolarScalingFactorsSim( float* nonPolarScalingFactors )
{
    cudaError_t status;
    gbsaSim.pNonPolarScalingFactors = nonPolarScalingFactors;
    status                          = cudaMemcpyToSymbol(gbsaSimDev, &gbsaSim, sizeof(cudaFreeEnergySimulationObcGbsaSoftcore));
    RTERROR(status, "cudaMemcpyToSymbol: SetCalculateObcGbsaSoftcoreNonPolarScalingFactorsSim");

    //(void) fprintf( stderr, "In SetCalculateObcGbsaSoftcoreNonPolarScalingFactorsSim\n" );
}

void GetCalculateObcGbsaSoftcoreBornSumSim(gpuContext gpu)
{
    cudaError_t status;
    status = cudaMemcpyFromSymbol(&gpu->sim, cSim, sizeof(cudaGmxSimulation));     
    RTERROR(status, "GetCalculateObcGbsaSoftcoreBornSumSim: cudaMemcpyFromSymbol: SetSim copy from cSim failed");
}

__global__ void kReduceObcGbsaSoftcoreBornForces_kernel()
{
    unsigned int pos = (blockIdx.x * blockDim.x + threadIdx.x);
    float energy = 0.0f;
    while (pos < cSim.atoms)
    {
        float bornRadius         = cSim.pBornRadii[pos];
        float obcChain           = cSim.pObcChain[pos];
        float2 obcData           = cSim.pObcData[pos];
        float  nonPolarScaleData = gbsaSimDev.pNonPolarScalingFactors[pos];
        float totalForce = 0.0f;
        float* pFt = cSim.pBornForce + pos;

        int i = cSim.nonbondOutputBuffers;
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

        /* E */
        energy += saTerm;

        totalForce *= bornRadius * bornRadius * obcChain;

        pFt = cSim.pBornForce + pos;
        *pFt = totalForce;
        pos += gridDim.x * blockDim.x;
    }
    /* E */
    // correct for surface area factor of -6
    cSim.pEnergy[blockIdx.x * blockDim.x + threadIdx.x] += energy / -6.0f;
}


void kReduceObcGbsaSoftcoreBornForces(gpuContext gpu)
{

    kReduceObcGbsaSoftcoreBornForces_kernel<<<gpu->sim.blocks, gpu->sim.bf_reduce_threads_per_block>>>();
    LAUNCHERROR("kReduceObcGbsaBornForces");
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

#if 0
__global__ void kClearObcGbsaBornSum_kernel()
{
    unsigned int pos = blockIdx.x * blockDim.x + threadIdx.x;
    while (pos < cSim.stride * cSim.nonbondOutputBuffers)
    {
        ((float*)cSim.pBornSum)[pos] = 0.0f;
        pos += gridDim.x * blockDim.x;
    }
}

__global__ void kReduceObcGbsaBornSum_kernel()
{
    unsigned int pos = (blockIdx.x * blockDim.x + threadIdx.x);
    
    while (pos < cSim.atoms)
    {
        float sum = 0.0f;
        float* pSt = cSim.pBornSum + pos;
        float2 atom = cSim.pObcData[pos];
        
        // Get summed Born data
        for (int i = 0; i < cSim.nonbondOutputBuffers; i++)
        {
            sum += *pSt;
       //     printf("%4d %4d A: %9.4f\n", pos, i, *pSt);
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
        cSim.pBornRadii[pos] = bornRadius;
        cSim.pObcChain[pos]  = obcChain;
        pos += gridDim.x * blockDim.x;
    }   
}

void kReduceObcGbsaBornSum(gpuContext gpu)
{
//    printf("kReduceObcGbsaBornSum\n");
    kReduceObcGbsaBornSum_kernel<<<gpu->sim.blocks, 384>>>();
    gpu->bRecalculateBornRadii = false;
    LAUNCHERROR("kReduceObcGbsaBornSum");
}
#endif

/** 
 * Initialize parameters for Cuda Obc softcore
 * 
 * @param gpu                  gpu context
 * @param innerDielectric      solute dielectric
 * @param solventDielectric    solvent dielectric
 * @param radius               intrinsic Born radii
 * @param scale                Obc scaling factors
 * @param charge               atomic charges (possibly overwritten by other methods?)
 * @param nonPolarScalingFactors non-polar scaling factors
 *
 */

extern "C"
GpuObcGbsaSoftcore* gpuSetObcSoftcoreParameters(gpuContext gpu, float innerDielectric, float solventDielectric, float nonPolarPrefactor,
                                                const std::vector<float>& radius, const std::vector<float>& scale,
                                                const std::vector<float>& charge, const std::vector<float>& nonPolarScalingFactors)
{

// ---------------------------------------------------------------------------------------

   static const float dielectricOffset    =    0.009f;
   static const float electricConstant    = -166.02691f;
   static const std::string methodName    = "gpuSetObcSoftcoreParameters";

// ---------------------------------------------------------------------------------------

    unsigned int atoms                     = radius.size();

    // initialize parameters

//    gpu->bIncludeGBSA = true;
    GpuObcGbsaSoftcore* gpuObcGbsaSoftcore = new GpuObcGbsaSoftcore();
    gpuObcGbsaSoftcore->initializeNonPolarScalingFactors( gpu->sim.paddedNumberOfAtoms );

    gpu->sim.surfaceAreaFactor             =  -6.0f*PI*4.0f*nonPolarPrefactor;
    for (unsigned int i = 0; i < atoms; i++)
    {
            (*gpu->psObcData)[i].x = radius[i] - dielectricOffset;
            (*gpu->psObcData)[i].y = scale[i] * (*gpu->psObcData)[i].x;
            (*gpu->psPosq4)[i].w   = charge[i];
            gpuObcGbsaSoftcore->setNonPolarScalingFactors( i, nonPolarScalingFactors[i] );

    }

    // diagnostics
#define DUMP_PARAMETERS 0
#if (DUMP_PARAMETERS == 1)
    (void) fprintf( stderr, "%s %u %u\n", methodName.c_str(), gpu->natoms, gpu->sim.paddedNumberOfAtoms );
    for (unsigned int i = 0; i < atoms; i++)
    {
       (void) fprintf( stderr, "%6u %13.6e %13.6e %8.3f %8.3f\n", i,  (*gpu->psObcData)[i].x, (*gpu->psObcData)[i].y, (*gpu->psPosq4)[i].w , nonPolarScalingFactors[i] );
    }
#endif

    // dummy out extra atom data

    for (unsigned int i = gpu->natoms; i < gpu->sim.paddedNumberOfAtoms; i++)
    {
        (*gpu->psBornRadii)[i]     = 0.2f;
        (*gpu->psObcData)[i].x     = 0.01f;
        (*gpu->psObcData)[i].y     = 0.01f;
    }

    // load data to board

    gpuObcGbsaSoftcore->upload( gpu );
    gpu->psBornRadii->Upload();
    gpu->psObcData->Upload();
    gpu->psPosq4->Upload();

    gpu->sim.preFactor = 2.0f*electricConstant*((1.0f/innerDielectric)-(1.0f/solventDielectric))*gpu->sim.forceConversionFactor;

    return gpuObcGbsaSoftcore;
}

void kCalculateObcGbsaSoftcoreBornSum(gpuContext gpu)
{
  //  printf("kCalculateObcGbsaSoftcoreBornSum\n");
    kClearObcGbsaBornSum( gpu );
    LAUNCHERROR("kClearBornSum from kCalculateObcGbsaSoftcoreBornSum");

    switch (gpu->sim.nonbondedMethod)
    {
        case NO_CUTOFF:
#define GBSA 0
#if GBSA == 1
gpu->psWorkUnit->Download();
fprintf( stderr, "kCalculateObcGbsaSoftcoreBornSum: bOutputBufferPerWarp=%u blks=%u th/blk=%u wu=%u %u\n", gpu->bOutputBufferPerWarp,
                 gpu->sim.nonbond_blocks, gpu->sim.nonbond_threads_per_block, gpu->sim.workUnits, gpu->psWorkUnit->_pSysData[0] );
#endif
#undef GBSA

            if (gpu->bOutputBufferPerWarp)
                kCalculateObcGbsaN2ByWarpBornSum_kernel<<<gpu->sim.nonbond_blocks, gpu->sim.nonbond_threads_per_block,
                        sizeof(Atom)*gpu->sim.nonbond_threads_per_block>>>(gpu->sim.pWorkUnit);
            else
                kCalculateObcGbsaN2BornSum_kernel<<<gpu->sim.nonbond_blocks, gpu->sim.nonbond_threads_per_block,
                        sizeof(Atom)*gpu->sim.nonbond_threads_per_block>>>(gpu->sim.pWorkUnit);
            break;
        case CUTOFF:
            if (gpu->bOutputBufferPerWarp)
                kCalculateObcGbsaCutoffByWarpBornSum_kernel<<<gpu->sim.nonbond_blocks, gpu->sim.nonbond_threads_per_block,
                        (sizeof(Atom)+sizeof(float))*gpu->sim.nonbond_threads_per_block>>>(gpu->sim.pInteractingWorkUnit);
            else
                kCalculateObcGbsaCutoffBornSum_kernel<<<gpu->sim.nonbond_blocks, gpu->sim.nonbond_threads_per_block,
                        (sizeof(Atom)+sizeof(float))*gpu->sim.nonbond_threads_per_block>>>(gpu->sim.pInteractingWorkUnit);
            break;
        case PERIODIC:
            if (gpu->bOutputBufferPerWarp)
                kCalculateObcGbsaPeriodicByWarpBornSum_kernel<<<gpu->sim.nonbond_blocks, gpu->sim.nonbond_threads_per_block,
                        (sizeof(Atom)+sizeof(float))*gpu->sim.nonbond_threads_per_block>>>(gpu->sim.pInteractingWorkUnit);
            else
                kCalculateObcGbsaPeriodicBornSum_kernel<<<gpu->sim.nonbond_blocks, gpu->sim.nonbond_threads_per_block,
                        (sizeof(Atom)+sizeof(float))*gpu->sim.nonbond_threads_per_block>>>(gpu->sim.pInteractingWorkUnit);
            break;
    }
    LAUNCHERROR("kCalculateObcGbsaSoftcoreBornSum");
}
