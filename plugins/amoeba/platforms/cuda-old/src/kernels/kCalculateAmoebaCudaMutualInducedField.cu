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

#include <stdio.h>

using namespace std;

static __constant__ cudaGmxSimulation cSim;
static __constant__ cudaAmoebaGmxSimulation cAmoebaSim;

void SetCalculateAmoebaCudaMutualInducedFieldSim(amoebaGpuContext amoebaGpu)
{
    cudaError_t status;
    gpuContext gpu = amoebaGpu->gpuContext;
    status         = cudaMemcpyToSymbol(cSim, &gpu->sim, sizeof(cudaGmxSimulation));    
    RTERROR(status, "SetCalculateAmoebaCudaMutualInducedFieldSim: cudaMemcpyToSymbol: SetSim copy to cSim failed");
    status         = cudaMemcpyToSymbol(cAmoebaSim, &amoebaGpu->amoebaSim, sizeof(cudaAmoebaGmxSimulation));    
    RTERROR(status, "SetCalculateAmoebaCudaMutualInducedFieldSim: cudaMemcpyToSymbol: SetSim copy to cAmoebaSim failed");
}

void GetCalculateAmoebaCudaMutualInducedFieldSim(amoebaGpuContext amoebaGpu)
{
    cudaError_t status;
    gpuContext gpu = amoebaGpu->gpuContext;
    status = cudaMemcpyFromSymbol(&gpu->sim, cSim, sizeof(cudaGmxSimulation));    
    RTERROR(status, "GetCalculateAmoebaCudaMutualInducedFieldSim: cudaMemcpyFromSymbol: SetSim copy from cSim failed");
    status = cudaMemcpyFromSymbol(&amoebaGpu->amoebaSim, cAmoebaSim, sizeof(cudaAmoebaGmxSimulation));    
    RTERROR(status, "GetCalculateAmoebaCudaMutualInducedFieldSim: cudaMemcpyFromSymbol: SetSim copy from cAmoebaSim failed");
}

#include "kCalculateAmoebaCudaMutualInducedParticle.h"

__device__ void calculateMutualInducedFieldPairIxn_kernel( MutualInducedParticle& atomI, MutualInducedParticle& atomJ,
                                                           float fields[4][3] )
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
    float pGamma                                 = atomJ.thole > atomI.thole ? atomI.thole: atomJ.thole;
    float damp                                   = ratio*ratio*ratio*pGamma;
    float dampExp                                = ( (dampProd != 0.0f) && (r < cAmoebaSim.scalingDistanceCutoff) ) ? expf( -damp ) : 0.0f; 

    rr3                                         *= (1.0f - dampExp);
    rr5                                         *= (1.0f - ( 1.0f + damp )*dampExp);
        
    float dDotDelta                              = rr5*(deltaR[0]*atomJ.inducedDipole[0]         + deltaR[1]*atomJ.inducedDipole[1]       + deltaR[2]*atomJ.inducedDipole[2] );
    fields[0][0]                                 = rr3*atomJ.inducedDipole[0] + dDotDelta*deltaR[0];
    fields[0][1]                                 = rr3*atomJ.inducedDipole[1] + dDotDelta*deltaR[1];
    fields[0][2]                                 = rr3*atomJ.inducedDipole[2] + dDotDelta*deltaR[2];
   
    dDotDelta                                    = rr5*(deltaR[0]*atomJ.inducedDipolePolar[0]    + deltaR[1]*atomJ.inducedDipolePolar[1]  + deltaR[2]*atomJ.inducedDipolePolar[2] );
    fields[1][0]                                 = rr3*atomJ.inducedDipolePolar[0] + dDotDelta*deltaR[0];
    fields[1][1]                                 = rr3*atomJ.inducedDipolePolar[1] + dDotDelta*deltaR[1];
    fields[1][2]                                 = rr3*atomJ.inducedDipolePolar[2] + dDotDelta*deltaR[2];
  
    dDotDelta                                    = rr5*(deltaR[0]*atomI.inducedDipole[0]         + deltaR[1]*atomI.inducedDipole[1]       + deltaR[2]*atomI.inducedDipole[2] );
    fields[2][0]                                 = rr3*atomI.inducedDipole[0] + dDotDelta*deltaR[0];
    fields[2][1]                                 = rr3*atomI.inducedDipole[1] + dDotDelta*deltaR[1];
    fields[2][2]                                 = rr3*atomI.inducedDipole[2] + dDotDelta*deltaR[2];
   
    dDotDelta                                    = rr5*(deltaR[0]*atomI.inducedDipolePolar[0]    + deltaR[1]*atomI.inducedDipolePolar[1]  + deltaR[2]*atomI.inducedDipolePolar[2] );
    fields[3][0]                                 = rr3*atomI.inducedDipolePolar[0] + dDotDelta*deltaR[0];
    fields[3][1]                                 = rr3*atomI.inducedDipolePolar[1] + dDotDelta*deltaR[1];
    fields[3][2]                                 = rr3*atomI.inducedDipolePolar[2] + dDotDelta*deltaR[2];
}

// Include versions of the kernels for N^2 calculations.

#define METHOD_NAME(a, b) a##N2##b
#include "kCalculateAmoebaCudaMutualInducedField.h"
#define USE_OUTPUT_BUFFER_PER_WARP
#undef METHOD_NAME
#define METHOD_NAME(a, b) a##N2ByWarp##b
#include "kCalculateAmoebaCudaMutualInducedField.h"

__global__ 
#if (__CUDA_ARCH__ >= 200)
__launch_bounds__(GF1XX_THREADS_PER_BLOCK, 1)
#elif (__CUDA_ARCH__ >= 120)
__launch_bounds__(GT2XX_THREADS_PER_BLOCK, 1)
#else
__launch_bounds__(G8X_THREADS_PER_BLOCK, 1)
#endif
void kInitializeMutualInducedField_kernel(
                   int numberOfAtoms,
                   float* fixedEField,
                   float* fixedEFieldPolar,
                   float* polarizability )
{

    int pos = blockIdx.x*blockDim.x + threadIdx.x;
    while( pos < 3*cSim.atoms )
    {

        fixedEField[pos]         *= polarizability[pos];
        fixedEFieldPolar[pos]    *= polarizability[pos];
       
        pos                      += blockDim.x*gridDim.x;
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
void kReduceMutualInducedFieldDelta_kernel(int numberOfEntries, float* arrayOfDeltas1, float* arrayOfDeltas2, float* epsilon )
{
    extern __shared__ float2 delta[];

    delta[threadIdx.x].x    = 0.0f;
    delta[threadIdx.x].y    = 0.0f;

    unsigned int pos = threadIdx.x;

    // load deltas

    while( pos < numberOfEntries )
    {   
        delta[threadIdx.x].x  += arrayOfDeltas1[pos];
        delta[threadIdx.x].y  += arrayOfDeltas2[pos];
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
        }
        __syncthreads();
    }   

    // set epsilons

    if (threadIdx.x == 0)
    {   
        epsilon[0]  = delta[0].x > delta[0].y ? delta[0].x : delta[0].y;
        epsilon[0]  = 48.033324f*sqrtf( epsilon[0]/( (float) (numberOfEntries/3)) );
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
void kSorUpdateMutualInducedField_kernel(
                   int numberOfEntries,    float* polarizability,
                   float* inducedDipole, float* inducedDipoleP,
                   float* fixedEField,   float* fixedEFieldP,
                   float* matrixProduct, float* matrixProductP )
{

    float polarSOR = 0.55f;
    int pos        = blockIdx.x*blockDim.x + threadIdx.x;
    while( pos < 3*cSim.atoms )
    {

        float previousDipole       = inducedDipole[pos];
        float previousDipoleP      = inducedDipoleP[pos];
    
        inducedDipole[pos]         = fixedEField[pos]     + polarizability[pos]*matrixProduct[pos];
        inducedDipoleP[pos]        = fixedEFieldP[pos]    + polarizability[pos]*matrixProductP[pos];
    
        inducedDipole[pos]         = previousDipole   + polarSOR*( inducedDipole[pos]   - previousDipole  );   
        inducedDipoleP[pos]        = previousDipoleP  + polarSOR*( inducedDipoleP[pos]  - previousDipoleP );
    
        matrixProduct[pos]         = ( inducedDipole[pos]  - previousDipole  )*( inducedDipole[pos]  - previousDipole  );
        matrixProductP[pos]        = ( inducedDipoleP[pos] - previousDipoleP )*( inducedDipoleP[pos] - previousDipoleP );
 
        pos                       += blockDim.x*gridDim.x;
    }

}

// reduce psWorkArray_3_1
// reduce psWorkArray_3_2

static void kReduceMutualInducedFields(amoebaGpuContext amoebaGpu, CUDAStream<float>* outputArray, CUDAStream<float>* outputPolarArray )
{
    gpuContext gpu = amoebaGpu->gpuContext;
    kReduceFields_kernel<<<gpu->sim.nonbond_blocks, gpu->sim.bsf_reduce_threads_per_block>>>(
                               gpu->sim.paddedNumberOfAtoms*3, gpu->sim.outputBuffers,
                               amoebaGpu->psWorkArray_3_1->_pDevData, outputArray->_pDevData, 0 );
    LAUNCHERROR("kReduceMI_Fields1");

    kReduceFields_kernel<<<gpu->sim.nonbond_blocks, gpu->sim.bsf_reduce_threads_per_block>>>(
                               gpu->sim.paddedNumberOfAtoms*3, gpu->sim.outputBuffers,
                               amoebaGpu->psWorkArray_3_2->_pDevData, outputPolarArray->_pDevData, 0 );
    LAUNCHERROR("kReduceMI_Fields2");
}

/**---------------------------------------------------------------------------------------

   Compute mutual induce field

   @param amoebaGpu        amoebaGpu context

   --------------------------------------------------------------------------------------- */

static void cudaComputeAmoebaMutualInducedFieldMatrixMultiply( amoebaGpuContext amoebaGpu,
                                                               CUDAStream<float>* outputArray, CUDAStream<float>* outputPolarArray )
{
  
   // ---------------------------------------------------------------------------------------

    static unsigned int threadsPerBlock = 0;

   // ---------------------------------------------------------------------------------------
  
    gpuContext gpu    = amoebaGpu->gpuContext;

    kClearFields_3( amoebaGpu, 2 );

    if( threadsPerBlock == 0 ){  
        unsigned int maxThreads;
        if (gpu->sm_version >= SM_20)
            maxThreads = 512; 
        else if (gpu->sm_version >= SM_12)
            maxThreads = 128; 
        else 
            maxThreads = 64; 
        threadsPerBlock = std::min(getThreadsPerBlock(amoebaGpu, sizeof(MutualInducedParticle), gpu->sharedMemoryPerBlock ), maxThreads);
    }   

    if (gpu->bOutputBufferPerWarp){
        kCalculateAmoebaMutualInducedFieldN2ByWarp_kernel<<<gpu->sim.nonbond_blocks, threadsPerBlock, sizeof(MutualInducedParticle)*threadsPerBlock>>>(
                                                                 gpu->psWorkUnit->_pDevData,
                                                                 amoebaGpu->psWorkArray_3_1->_pDevData,
                                                                 amoebaGpu->psWorkArray_3_2->_pDevData );

    } else {

        kCalculateAmoebaMutualInducedFieldN2_kernel<<<gpu->sim.nonbond_blocks, threadsPerBlock, sizeof(MutualInducedParticle)*threadsPerBlock>>>(
                                                                 gpu->psWorkUnit->_pDevData,
                                                                 amoebaGpu->psWorkArray_3_1->_pDevData,
                                                                 amoebaGpu->psWorkArray_3_2->_pDevData );
    }
    LAUNCHERROR("kCalculateAmoebaMutualInducedField");

    kReduceMutualInducedFields( amoebaGpu, outputArray, outputPolarArray );

}

/**---------------------------------------------------------------------------------------

   Compute mutual induce field

   @param amoebaGpu        amoebaGpu context

   --------------------------------------------------------------------------------------- */

static void cudaComputeAmoebaMutualInducedFieldBySOR( amoebaGpuContext amoebaGpu )
{
  
   // ---------------------------------------------------------------------------------------

    int done;
    int iteration;

    gpuContext gpu    = amoebaGpu->gpuContext;

   // ---------------------------------------------------------------------------------------

    // set  E_Field & E_FieldPolar] to [ E_Field & E_FieldPolar]*Polarizability
    // initialize [ InducedDipole & InducedDipolePolar ] to [ E_Field & E_FieldPolar]*Polarizability

    kInitializeMutualInducedField_kernel<<< gpu->sim.nonbond_blocks, gpu->sim.bsf_reduce_threads_per_block >>>(
         gpu->natoms,
         amoebaGpu->psE_Field->_pDevData,
         amoebaGpu->psE_FieldPolar->_pDevData,
         amoebaGpu->psPolarizability->_pDevData );
    LAUNCHERROR("AmoebaMutualInducedFieldSetup");  

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

        cudaComputeAmoebaMutualInducedFieldMatrixMultiply( amoebaGpu, amoebaGpu->psWorkVector[0],  amoebaGpu->psWorkVector[1] );
        LAUNCHERROR("cudaComputeAmoebaMutualInducedFieldMatrixMultiply Loop\n");  

        // post matrix multiply

        kSorUpdateMutualInducedField_kernel<<< gpu->sim.nonbond_blocks, gpu->sim.bsf_reduce_threads_per_block >>>(
           gpu->natoms, amoebaGpu->psPolarizability->_pDevData,
           amoebaGpu->psInducedDipole->_pDevData, amoebaGpu->psInducedDipolePolar->_pDevData,
           amoebaGpu->psE_Field->_pDevData,       amoebaGpu->psE_FieldPolar->_pDevData,
           amoebaGpu->psWorkVector[0]->_pDevData,     amoebaGpu->psWorkVector[1]->_pDevData );
        LAUNCHERROR("kSorUpdateMutualInducedField");  

        // get total epsilon -- performing sums on gpu

        kReduceMutualInducedFieldDelta_kernel<<<1, amoebaGpu->epsilonThreadsPerBlock, 2*sizeof(float)*amoebaGpu->epsilonThreadsPerBlock>>>(
           3*gpu->natoms, amoebaGpu->psWorkVector[0]->_pDevData, amoebaGpu->psWorkVector[1]->_pDevData,
           amoebaGpu->psCurrentEpsilon->_pDevData );
        LAUNCHERROR("kReduceMutualInducedFieldDelta");

        // Debye=48.033324f
        amoebaGpu->psCurrentEpsilon->Download();
        float currentEpsilon          = amoebaGpu->psCurrentEpsilon->_pSysData[0];
        amoebaGpu->mutualInducedCurrentEpsilon   = currentEpsilon;

        if( iteration > amoebaGpu->mutualInducedMaxIterations || amoebaGpu->mutualInducedCurrentEpsilon < amoebaGpu->mutualInducedTargetEpsilon ){ 
            done = 1;
        }
        iteration++;
    }

    amoebaGpu->mutualInducedDone             = done;
    amoebaGpu->mutualInducedConverged        = ( !done || iteration > amoebaGpu->mutualInducedMaxIterations ) ? 0 : 1;

}

void cudaComputeAmoebaMutualInducedField( amoebaGpuContext amoebaGpu )
{
    if( amoebaGpu->mutualInducedIterativeMethod == 0 ){
        cudaComputeAmoebaMutualInducedFieldBySOR( amoebaGpu );
    }
}
