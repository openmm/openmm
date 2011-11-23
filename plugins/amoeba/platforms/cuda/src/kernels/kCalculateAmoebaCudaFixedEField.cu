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

#include "amoebaCudaKernels.h"
#include "kCalculateAmoebaCudaUtilities.h"

static __constant__ cudaGmxSimulation cSim;
static __constant__ cudaAmoebaGmxSimulation cAmoebaSim;

void SetCalculateAmoebaCudaFixedEFieldSim(amoebaGpuContext amoebaGpu)
{
    cudaError_t status;
    gpuContext gpu = amoebaGpu->gpuContext;
    status         = cudaMemcpyToSymbol(cSim, &gpu->sim, sizeof(cudaGmxSimulation));    
    RTERROR(status, "SetCalculateAmoebaCudaFixedEFieldSim: cudaMemcpyToSymbol: SetSim copy to cSim failed");
    status         = cudaMemcpyToSymbol(cAmoebaSim, &amoebaGpu->amoebaSim, sizeof(cudaAmoebaGmxSimulation));    
    RTERROR(status, "SetCalculateAmoebaCudaFixedEFieldSim: cudaMemcpyToSymbol: SetSim copy to cAmoebaSim failed");
}

void GetCalculateAmoebaCudaFixedEFieldSim(amoebaGpuContext amoebaGpu)
{
    cudaError_t status;
    gpuContext gpu = amoebaGpu->gpuContext;
    status = cudaMemcpyFromSymbol(&gpu->sim, cSim, sizeof(cudaGmxSimulation));    
    RTERROR(status, "GetCalculateAmoebaCudaFixedEFieldSim: cudaMemcpyFromSymbol: SetSim copy from cSim failed");
    status = cudaMemcpyFromSymbol(&amoebaGpu->amoebaSim, cAmoebaSim, sizeof(cudaAmoebaGmxSimulation));         
    RTERROR(status, "GetCalculateAmoebaCudaFixedEFieldSim: cudaMemcpyFromSymbol: SetSim copy from cAmoebaSim failed");
}

// reduce psWorkArray_3_1 -> EField
// reduce psWorkArray_3_2 -> EFieldPolar

static void kReduceE_Fields_kernel(amoebaGpuContext amoebaGpu )
{
    gpuContext gpu = amoebaGpu->gpuContext;
    kReduceFields_kernel<<<gpu->sim.nonbond_blocks, gpu->sim.bsf_reduce_threads_per_block>>>(
                               gpu->sim.paddedNumberOfAtoms*3, gpu->sim.outputBuffers,
                               amoebaGpu->psWorkArray_3_1->_pDevData, amoebaGpu->psE_Field->_pDevData, 0 );
    LAUNCHERROR("kReduceE_Fields1");

    kReduceFields_kernel<<<gpu->sim.nonbond_blocks, gpu->sim.bsf_reduce_threads_per_block>>>(
                               gpu->sim.paddedNumberOfAtoms*3, gpu->sim.outputBuffers,
                               amoebaGpu->psWorkArray_3_2->_pDevData, amoebaGpu->psE_FieldPolar->_pDevData, 0 );
    LAUNCHERROR("kReduceE_Fields2");
}

// file includes FixedFieldParticle struct definition/load/unload struct and body kernel for fixed E-field

#undef GK
#include "kCalculateAmoebaCudaFixedFieldParticle.h"

// Include versions of the kernels for N^2 calculations.

#define METHOD_NAME(a, b) a##N2##b
#include "kCalculateAmoebaCudaFixedEField.h"
#define USE_OUTPUT_BUFFER_PER_WARP
#undef METHOD_NAME
#define METHOD_NAME(a, b) a##N2ByWarp##b
#include "kCalculateAmoebaCudaFixedEField.h"

/**---------------------------------------------------------------------------------------

   Compute fixed electric field

   @param amoebaGpu        amoebaGpu context
   @param gpu              OpenMM gpu Cuda context

   --------------------------------------------------------------------------------------- */

void cudaComputeAmoebaFixedEField( amoebaGpuContext amoebaGpu )
{
  
    gpuContext gpu    = amoebaGpu->gpuContext;

    kClearFields_3( amoebaGpu, 2 );

    static unsigned int threadsPerBlock = 0;
    if( threadsPerBlock == 0 ){ 
        unsigned int maxThreads;
        if (gpu->sm_version >= SM_20)
            maxThreads = 512; 
        else if (gpu->sm_version >= SM_12)
            maxThreads = 128; 
        else 
            maxThreads = 64;
        threadsPerBlock = std::min(getThreadsPerBlock(amoebaGpu, sizeof(FixedFieldParticle), gpu->sharedMemoryPerBlock ), maxThreads);
    }

    if (gpu->bOutputBufferPerWarp){
        kCalculateAmoebaFixedE_FieldN2ByWarpForces_kernel<<<gpu->sim.nonbond_blocks, threadsPerBlock, sizeof(FixedFieldParticle)*threadsPerBlock>>>(
                                                                           gpu->psWorkUnit->_pDevData,
                                                                           amoebaGpu->psWorkArray_3_1->_pDevData,
                                                                           amoebaGpu->psWorkArray_3_2->_pDevData );
    } else {

        kCalculateAmoebaFixedE_FieldN2Forces_kernel<<<gpu->sim.nonbond_blocks, threadsPerBlock, sizeof(FixedFieldParticle)*threadsPerBlock>>>(
                                                                           gpu->psWorkUnit->_pDevData,
                                                                           amoebaGpu->psWorkArray_3_1->_pDevData,
                                                                           amoebaGpu->psWorkArray_3_2->_pDevData );
    }

    LAUNCHERROR("kCalculateAmoebaFixedE_FieldN2Forces_kernel");
    kReduceE_Fields_kernel( amoebaGpu );
}
