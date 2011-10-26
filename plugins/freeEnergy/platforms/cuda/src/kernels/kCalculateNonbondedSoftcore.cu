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

#include "freeEnergyGpuTypes.h"
#include "GpuFreeEnergyCudaKernels.h"
#include "openmm/OpenMMException.h"
#include <iostream>
#include <sstream>

#define PARAMETER_PRINT 0
#define MAX_PARAMETER_PRINT 10

// device handles

static __constant__ cudaGmxSimulation cSim;
static __constant__ cudaFreeEnergyGmxSimulation feSimDev;

// write address of structs to devices

void SetCalculateCDLJSoftcoreGpuSim( freeEnergyGpuContext freeEnergyGpu ){
    cudaError_t status;
    status = cudaMemcpyToSymbol(cSim, &freeEnergyGpu->gpuContext->sim, sizeof(cudaGmxSimulation));
    RTERROR(status, "cudaMemcpyToSymbol: SetCalculateCDLJSoftcoreGpuSim copy to cSim failed");

    status = cudaMemcpyToSymbol( feSimDev, &freeEnergyGpu->freeEnergySim, sizeof(cudaFreeEnergyGmxSimulation));
    RTERROR(status, "cudaMemcpyToSymbol: SetCalculateCDLJSoftcoreGpuSim copy to feSimDev failed");
}

extern "C"
void freeEnergyGpuSetPeriodicBoxSize( freeEnergyGpuContext freeEnergyGpu, float xsize, float ysize, float zsize)
{
    freeEnergyGpu->freeEnergySim.periodicBoxSizeX    = xsize;
    freeEnergyGpu->freeEnergySim.periodicBoxSizeY    = ysize;
    freeEnergyGpu->freeEnergySim.periodicBoxSizeZ    = zsize;

    freeEnergyGpu->freeEnergySim.invPeriodicBoxSizeX = 1.0f/xsize;
    freeEnergyGpu->freeEnergySim.invPeriodicBoxSizeY = 1.0f/ysize;
    freeEnergyGpu->freeEnergySim.invPeriodicBoxSizeZ = 1.0f/zsize;

    freeEnergyGpu->freeEnergySim.recipBoxSizeX       = 2.0f*PI/freeEnergyGpu->freeEnergySim.periodicBoxSizeX;
    freeEnergyGpu->freeEnergySim.recipBoxSizeY       = 2.0f*PI/freeEnergyGpu->freeEnergySim.periodicBoxSizeY;
    freeEnergyGpu->freeEnergySim.recipBoxSizeZ       = 2.0f*PI/freeEnergyGpu->freeEnergySim.periodicBoxSizeZ;

    freeEnergyGpu->freeEnergySim.cellVolume          = freeEnergyGpu->freeEnergySim.periodicBoxSizeX*freeEnergyGpu->freeEnergySim.periodicBoxSizeY*freeEnergyGpu->freeEnergySim.periodicBoxSizeZ;

    gpuSetPeriodicBoxSize( freeEnergyGpu->gpuContext, xsize, ysize, zsize );
}

extern "C"
void gpuSetNonbondedSoftcoreParameters( freeEnergyGpuContext freeEnergyGpu, float epsfac, const std::vector<int>& atom, const std::vector<float>& c6,
                                        const std::vector<float>& c12, const std::vector<float>& q,
                                        const std::vector<float>& softcoreLJLambdaArray, const std::vector<char>& symbol,
                                        const std::vector<std::vector<int> >& exclusions, CudaFreeEnergyNonbondedMethod method,
                                        float cutoffDistance, float solventDielectric ){

    unsigned int numberOfParticles                         = c6.size();
    gpuContext gpu                                         = freeEnergyGpu->gpuContext;
    int paddedNumberOfAtoms                                = gpu->sim.paddedNumberOfAtoms;

    // sanity checks

    if( paddedNumberOfAtoms < 1 ){
        std::stringstream msg;
        msg << "gpuSetNonbondedSoftcoreParameters: number of padded atoms=" <<  gpu->sim.paddedNumberOfAtoms << " is less than 1.";
        throw OpenMM::OpenMMException( msg.str() );
    }

    if( freeEnergyGpu->gpuContext->sim.atoms != numberOfParticles  ){
        std::stringstream msg;
        msg << "gpuSetNonbondedSoftcoreParameters: number of atoms in gpuContext does not match input count: " << freeEnergyGpu->gpuContext->sim.atoms << " " << numberOfParticles << ".";
        throw OpenMM::OpenMMException( msg.str() );
    }

    freeEnergyGpu->freeEnergySim.epsfac                    = epsfac;
    freeEnergyGpu->freeEnergySim.nonbondedMethod           = method;

    freeEnergyGpu->freeEnergySim.nonbondedCutoff           = cutoffDistance;
    freeEnergyGpu->freeEnergySim.nonbondedCutoffSqr        = cutoffDistance*cutoffDistance;

    gpu->sim.nonbondedCutoff                               = cutoffDistance;
    gpu->sim.nonbondedCutoffSqr                            = cutoffDistance*cutoffDistance;

    if( cutoffDistance > 0.0f ){
        freeEnergyGpu->freeEnergySim.reactionFieldK        = pow(cutoffDistance, -3.0f)*(solventDielectric-1.0f)/(2.0f*solventDielectric+1.0f);
        freeEnergyGpu->freeEnergySim.reactionFieldC        = (1.0f / cutoffDistance)*(3.0f*solventDielectric)/(2.0f*solventDielectric+1.0f);
        gpu->sim.reactionFieldK                            = freeEnergyGpu->freeEnergySim.reactionFieldK;
        gpu->sim.reactionFieldC                            = freeEnergyGpu->freeEnergySim.reactionFieldC;
    } else {
        freeEnergyGpu->freeEnergySim.reactionFieldK        = 0.0f;
        freeEnergyGpu->freeEnergySim.reactionFieldC        = 0.0f;
    }

    setExclusions( gpu, exclusions );
 
    // parameters
 
    freeEnergyGpu->psSigEps4                               = new CUDAStream<float4>( paddedNumberOfAtoms, 1, "freeEnergyGpuSigEps4");
    freeEnergyGpu->freeEnergySim.pSigEps4                  = freeEnergyGpu->psSigEps4->_pDevData;

    for( unsigned int ii = 0; ii < numberOfParticles; ii++ ){

        float p1 = 0.5f;
        float p2 = 0.0f;               

        if( (c6[ii] > 0.0f) && (c12[ii] > 0.0f) ){
            p1 = 0.5f * powf(c12[ii] / c6[ii], 1.0f / 6.0f);
            p2 = c6[ii] * sqrtf(1.0f / c12[ii]);
        }
/*
            if (symbol.size() > 0)
                freeEnergyGpu->pAtomSymbol[ii] = symbol[ii];
*/

        (*freeEnergyGpu->psSigEps4)[ii].x        = p1;
        (*freeEnergyGpu->psSigEps4)[ii].y        = p2;
        (*freeEnergyGpu->psSigEps4)[ii].z        = softcoreLJLambdaArray[ii];
        (*freeEnergyGpu->psSigEps4)[ii].w        = q[ii];
    }

    // Dummy out extra atom data

    for( unsigned int ii = numberOfParticles; ii < paddedNumberOfAtoms; ii++ ){

        (*freeEnergyGpu->psSigEps4)[ii].x              = 1.0f;
        (*freeEnergyGpu->psSigEps4)[ii].y              = 0.0f;
        (*freeEnergyGpu->psSigEps4)[ii].z              = 0.0f;
        (*freeEnergyGpu->psSigEps4)[ii].w              = 0.0f;

        (*gpu->psPosq4)[ii].x                          = 100000.0f + ii * 10.0f;
        (*gpu->psPosq4)[ii].y                          = 100000.0f + ii * 10.0f;
        (*gpu->psPosq4)[ii].z                          = 100000.0f + ii * 10.0f;
        (*gpu->psPosq4)[ii].w                          = 0.0f;

    }

    if( freeEnergyGpu->log ){
        (void) fprintf( freeEnergyGpu->log,"freeEnergyGpuSetNonbondedSoftcoreParameters: %5u padded=%u epsfac=%14.7e method=%d cutoffDistance=%9.2f solventDielectric=%9.2f\n",
                        numberOfParticles, freeEnergyGpu->gpuContext->sim.paddedNumberOfAtoms, epsfac, method, cutoffDistance, solventDielectric );
#ifdef PARAMETER_PRINT
        int maxPrint = MAX_PARAMETER_PRINT;
        for (unsigned int ii = 0; ii < numberOfParticles; ii++){
            (void) fprintf( freeEnergyGpu->log,"%6u sig[%14.7e %14.7e] lambda=%10.3f q=%10.3f\n",
                            ii, 
                            (*freeEnergyGpu->psSigEps4)[ii].x, (*freeEnergyGpu->psSigEps4)[ii].y, (*freeEnergyGpu->psSigEps4)[ii].z, (*freeEnergyGpu->psSigEps4)[ii].w );
            if( ii == maxPrint && ii < freeEnergyGpu->gpuContext->sim.paddedNumberOfAtoms - maxPrint ){
               ii = numberOfParticles - maxPrint;
            }
        }
        unsigned int offset = paddedNumberOfAtoms - maxPrint;
        if( offset > 0 ){
            if( offset > numberOfParticles ){
                (void) fprintf( freeEnergyGpu->log,"Dummy padded entries\n" );
                for (unsigned int ii = offset; ii < paddedNumberOfAtoms; ii++){
                    (void) fprintf( freeEnergyGpu->log,"%6u sig[%14.7e %14.7e] lambda=%10.3f q=%10.3f\n",
                                    ii, 
                                    (*freeEnergyGpu->psSigEps4)[ii].x, (*freeEnergyGpu->psSigEps4)[ii].y, (*freeEnergyGpu->psSigEps4)[ii].z, (*freeEnergyGpu->psSigEps4)[ii].w );
                }
            }
        }
#endif
        (void) fflush( freeEnergyGpu->log );
    }
 
    // upload data to board

    freeEnergyGpu->psSigEps4->Upload();
    gpu->psPosq4->Upload();

    return;
}

extern "C"
bool gpuIsAvailableSoftcore()
{
    int deviceCount;
    cudaGetDeviceCount(&deviceCount);
    return (deviceCount > 0);
}

struct Atom {
    float x;
    float y;
    float z;
    float q;
    float sig;
    float eps;
    float softCoreLJLambda;
    float fx;
    float fy;
    float fz;
};

// Include versions of the kernels for N^2 calculations with softcore LJ.

#define USE_SOFTCORE_LJ
#ifdef USE_SOFTCORE_LJ
#include "kSoftcoreLJ.h"
#endif

#undef METHOD_NAME
#define METHOD_NAME(a, b) a##N2SoftcoreLJ##b
#undef USE_OUTPUT_BUFFER_PER_WARP
#include "kCalculateNonbondedSoftcore.h"

#define USE_OUTPUT_BUFFER_PER_WARP
#undef METHOD_NAME
#define METHOD_NAME(a, b) a##N2SoftcoreLJByWarp##b
#include "kCalculateNonbondedSoftcore.h"

// Include versions of the kernels with cutoffs.

#undef METHOD_NAME
#undef USE_OUTPUT_BUFFER_PER_WARP
#define USE_CUTOFF
#define METHOD_NAME(a, b) a##Cutoff##b
#include "kCalculateNonbondedSoftcore.h"
#include "kFindInteractingBlocks.h"
#define USE_OUTPUT_BUFFER_PER_WARP
#undef METHOD_NAME
#define METHOD_NAME(a, b) a##CutoffByWarp##b
#include "kCalculateNonbondedSoftcore.h"

// Include versions of the kernels with periodic boundary conditions.

#undef METHOD_NAME
#undef USE_OUTPUT_BUFFER_PER_WARP
#define USE_PERIODIC
#define METHOD_NAME(a, b) a##Periodic##b
#include "kCalculateNonbondedSoftcore.h"
#include "kFindInteractingBlocks.h"
#define USE_OUTPUT_BUFFER_PER_WARP
#undef METHOD_NAME
#define METHOD_NAME(a, b) a##PeriodicByWarp##b
#include "kCalculateNonbondedSoftcore.h"

void kCalculateCDLJSoftcoreForces( freeEnergyGpuContext freeEnergyGpu )
{

    gpuContext gpu = freeEnergyGpu->gpuContext;
    // (void) fprintf( stderr,"kCalculateCDLJCutoffForces %d warp=%u nonbond_blocks=%u nonbond_threads_per_block=%u rfK=%15.7e rfC=%15.7e\n", freeEnergyGpu->freeEnergySim.nonbondedMethod,
    //                 gpu->bOutputBufferPerWarp, gpu->sim.nonbond_blocks, gpu->sim.nonbond_threads_per_block, gpu->sim.reactionFieldK, gpu->sim.reactionFieldC); fflush( stderr );

    switch (freeEnergyGpu->freeEnergySim.nonbondedMethod)
    {
        case FREE_ENERGY_NO_CUTOFF:

           if (gpu->bOutputBufferPerWarp)
                kCalculateCDLJSoftcoreN2SoftcoreLJByWarpForces_kernel<<<gpu->sim.nonbond_blocks, gpu->sim.nonbond_threads_per_block,
                         sizeof(Atom)*gpu->sim.nonbond_threads_per_block>>>(gpu->sim.pWorkUnit);
           else
                   kCalculateCDLJSoftcoreN2SoftcoreLJForces_kernel<<<gpu->sim.nonbond_blocks, gpu->sim.nonbond_threads_per_block,
                           sizeof(Atom)*gpu->sim.nonbond_threads_per_block>>>(gpu->sim.pWorkUnit );
            LAUNCHERROR("kCalculateCDLJSoftcoreN2Forces");

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
                kCalculateCDLJSoftcoreCutoffByWarpForces_kernel<<<gpu->sim.nonbond_blocks, gpu->sim.nonbond_threads_per_block,
                        (sizeof(Atom)+sizeof(float3))*gpu->sim.nonbond_threads_per_block>>>(gpu->sim.pInteractingWorkUnit);
            else
                kCalculateCDLJSoftcoreCutoffForces_kernel<<<gpu->sim.nonbond_blocks, gpu->sim.nonbond_threads_per_block,
                        (sizeof(Atom)+sizeof(float3))*gpu->sim.nonbond_threads_per_block>>>(gpu->sim.pInteractingWorkUnit);

            LAUNCHERROR("kCalculateCDLJSoftcoreCutoffForces");
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
                kCalculateCDLJSoftcorePeriodicByWarpForces_kernel<<<gpu->sim.nonbond_blocks, gpu->sim.nonbond_threads_per_block,
                        (sizeof(Atom)+sizeof(float3))*gpu->sim.nonbond_threads_per_block>>>(gpu->sim.pInteractingWorkUnit);
            else
                kCalculateCDLJSoftcorePeriodicForces_kernel<<<gpu->sim.nonbond_blocks, gpu->sim.nonbond_threads_per_block,
                        (sizeof(Atom)+sizeof(float3))*gpu->sim.nonbond_threads_per_block>>>(gpu->sim.pInteractingWorkUnit);
            LAUNCHERROR("kCalculateCDLJSoftcorePeriodicForces");
            break;

        default:
            throw OpenMM::OpenMMException( "Nonbonded softcore method not recognized." );

    }
}

