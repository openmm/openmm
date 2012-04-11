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

using namespace std;

static __constant__ cudaGmxSimulation cSim;
static __constant__ cudaAmoebaGmxSimulation cAmoebaSim;

void SetCalculateAmoebaCudaPmeMutualInducedFieldSim(amoebaGpuContext amoebaGpu)
{
    cudaError_t status;
    gpuContext gpu = amoebaGpu->gpuContext;
    status         = cudaMemcpyToSymbol(cSim, &gpu->sim, sizeof(cudaGmxSimulation));    
    RTERROR(status, "SetCalculateAmoebaCudaPmeMutualInducedFieldSim: cudaMemcpyToSymbol: SetSim copy to cSim failed");
    status         = cudaMemcpyToSymbol(cAmoebaSim, &amoebaGpu->amoebaSim, sizeof(cudaAmoebaGmxSimulation));    
    RTERROR(status, "SetCalculateAmoebaCudaPmeMutualInducedFieldSim: cudaMemcpyToSymbol: SetSim copy to cAmoebaSim failed");
}

void GetCalculateAmoebaCudaPmeMutualInducedFieldSim(amoebaGpuContext amoebaGpu)
{
    cudaError_t status;
    gpuContext gpu = amoebaGpu->gpuContext;
    status = cudaMemcpyFromSymbol(&gpu->sim, cSim, sizeof(cudaGmxSimulation));    
    RTERROR(status, "GetCalculateAmoebaCudaPmeMutualInducedFieldSim: cudaMemcpyFromSymbol: SetSim copy from cSim failed");
    status = cudaMemcpyFromSymbol(&amoebaGpu->amoebaSim, cAmoebaSim, sizeof(cudaAmoebaGmxSimulation));    
    RTERROR(status, "GetCalculateAmoebaCudaPmeMutualInducedFieldSim: cudaMemcpyFromSymbol: SetSim copy from cAmoebaSim failed");
}

#undef INCLUDE_MI_FIELD_BUFFERS
#define INCLUDE_MI_FIELD_BUFFERS 
#include "kCalculateAmoebaCudaMutualInducedParticle.h"
#ifdef INCLUDE_MI_FIELD_BUFFERS
__device__ void sumTempBuffer( MutualInducedParticle& atomI, MutualInducedParticle& atomJ ){

    atomI.tempBuffer[0]  += atomJ.tempBuffer[0];
    atomI.tempBuffer[1]  += atomJ.tempBuffer[1];
    atomI.tempBuffer[2]  += atomJ.tempBuffer[2];

    atomI.tempBufferP[0] += atomJ.tempBufferP[0];
    atomI.tempBufferP[1] += atomJ.tempBufferP[1];
    atomI.tempBufferP[2] += atomJ.tempBufferP[2];
}
#endif

// file includes FixedFieldParticle struct definition/load/unload struct and body kernel for fixed E-field

__device__ void setupMutualInducedFieldPairIxn_kernel( const MutualInducedParticle& atomI, const MutualInducedParticle& atomJ,
                                                       const float uscale, float4* delta, float* preFactor2 ) {

    // compute thedelta->xeal space portion of the Ewald summation
  
    delta->x                = atomJ.x - atomI.x;
    delta->y                = atomJ.y - atomI.y;
    delta->z                = atomJ.z - atomI.z;

    // pdelta->xiodic boundary conditions

    delta->x               -= floorf(delta->x*cSim.invPeriodicBoxSizeX+0.5f)*cSim.periodicBoxSizeX;
    delta->y               -= floorf(delta->y*cSim.invPeriodicBoxSizeY+0.5f)*cSim.periodicBoxSizeY;
    delta->z               -= floorf(delta->z*cSim.invPeriodicBoxSizeZ+0.5f)*cSim.periodicBoxSizeZ;

    float r2                = (delta->x*delta->x) + (delta->y*delta->y) + (delta->z*delta->z); 
    if( r2 <= cSim.nonbondedCutoffSqr ){
        float r           = sqrtf(r2);

        // calculate the error function damping terms

        float ralpha      = cSim.alphaEwald*r;

        float bn0         = erfcf(ralpha)/r;
        float alsq2       = 2.0f*cSim.alphaEwald*cSim.alphaEwald;
        float alsq2n      = 1.0f/(cAmoebaSim.sqrtPi*cSim.alphaEwald);
        float exp2a       = expf(-(ralpha*ralpha));
        alsq2n           *= alsq2;
        float bn1         = (bn0+alsq2n*exp2a)/r2;

        alsq2n           *= alsq2;
        float bn2         = (3.0f*bn1+alsq2n*exp2a)/r2;

        // compute the error function scaled and unscaled terms

        float scale3      = 1.0f;
        float scale5      = 1.0f;
        float damp        = atomI.damp*atomJ.damp;
        if( damp != 0.0f ){

            float ratio  = (r/damp);
                  ratio  = ratio*ratio*ratio;
            float pgamma = atomI.thole < atomJ.thole ? atomI.thole : atomJ.thole;
                  damp   = -pgamma*ratio;

            if( damp > -50.0f) {
                float expdamp = expf(damp);
                scale3        = 1.0f - expdamp;
                scale5        = 1.0f - expdamp*(1.0f-damp);
            }
        }
        float dsc3        = uscale*scale3;
        float dsc5        = uscale*scale5;

        float r3          = (r*r2);
        float r5          = (r3*r2);
        float rr3         = (1.0f-dsc3)/r3;
        float rr5         = 3.0f*(1.0f-dsc5)/r5;

        delta->w          = rr3 - bn1;
        *preFactor2       = bn2 - rr5;
    } else {
        delta->w = *preFactor2 = 0.0f;
    }
}

__device__ void calculateMutualInducedFieldPairIxn_kernel( const float inducedDipole[3], const float4 delta, const float preFactor2, float fieldSum[3] ) {

    float preFactor3  = preFactor2*(inducedDipole[0]*delta.x   + inducedDipole[1]*delta.y  + inducedDipole[2]*delta.z);

    fieldSum[0]      += preFactor3*delta.x + delta.w*inducedDipole[0];
    fieldSum[1]      += preFactor3*delta.y + delta.w*inducedDipole[1];
    fieldSum[2]      += preFactor3*delta.z + delta.w*inducedDipole[2];
}

__device__ void calculateMutualInducedFieldPairIxnNoAdd_kernel( const float inducedDipole[3], const float4 delta, const float preFactor2, float fieldSum[3] ) {

    float preFactor3  = preFactor2*(inducedDipole[0]*delta.x   + inducedDipole[1]*delta.y  + inducedDipole[2]*delta.z);

    fieldSum[0]       = preFactor3*delta.x + delta.w*inducedDipole[0];
    fieldSum[1]       = preFactor3*delta.y + delta.w*inducedDipole[1];
    fieldSum[2]       = preFactor3*delta.z + delta.w*inducedDipole[2];
}

// file includes FixedFieldParticle struct definition/load/unload struct and body kernel for fixed E-field

__device__ void calculatePmeDirectMutualInducedFieldPairIxn_kernel( MutualInducedParticle& atomI, MutualInducedParticle& atomJ,
                                                                    float uscale, float4 fields[3] ){

    // compute the real space portion of the Ewald summation
  
    float xr          = atomJ.x - atomI.x;
    float yr          = atomJ.y - atomI.y;
    float zr          = atomJ.z - atomI.z;

    // periodic boundary conditions

    xr               -= floorf(xr*cSim.invPeriodicBoxSizeX+0.5f)*cSim.periodicBoxSizeX;
    yr               -= floorf(yr*cSim.invPeriodicBoxSizeY+0.5f)*cSim.periodicBoxSizeY;
    zr               -= floorf(zr*cSim.invPeriodicBoxSizeZ+0.5f)*cSim.periodicBoxSizeZ;

    float r2          = xr*xr + yr* yr + zr*zr;
    if( r2 <= cSim.nonbondedCutoffSqr ){
        float r           = sqrtf(r2);

        // calculate the error function damping terms

        float ralpha      = cSim.alphaEwald*r;

        float bn0         = erfcf(ralpha)/r;
        float alsq2       = 2.0f*cSim.alphaEwald*cSim.alphaEwald;
        float alsq2n      = 1.0f/(cAmoebaSim.sqrtPi*cSim.alphaEwald);
        float exp2a       = expf(-(ralpha*ralpha));
        alsq2n           *= alsq2;
        float bn1         = (bn0+alsq2n*exp2a)/r2;

        alsq2n           *= alsq2;
        float bn2         = (3.0f*bn1+alsq2n*exp2a)/r2;

        // compute the error function scaled and unscaled terms

        float scale3      = 1.0f;
        float scale5      = 1.0f;
        float damp        = atomI.damp*atomJ.damp;
        if( damp != 0.0f ){

            float ratio  = (r/damp);
                  ratio  = ratio*ratio*ratio;
            float pgamma = atomI.thole < atomJ.thole ? atomI.thole : atomJ.thole;
                  damp   = -pgamma*ratio;

            if( damp > -50.0f) {
                float expdamp = expf(damp);
                scale3        = 1.0f - expdamp;
                scale5        = 1.0f - expdamp*(1.0f-damp);
            }
        }
        float dsc3        = uscale*scale3;
        float dsc5        = uscale*scale5;

        float r3          = (r*r2);
        float r5          = (r3*r2);
        float rr3         = (1.0f-dsc3)/r3;
        float rr5         = 3.0f*(1.0f-dsc5)/r5;

        float preFactor1  = rr3 - bn1;
        float preFactor2  = bn2 - rr5;

        float dukr        = atomJ.inducedDipole[0]*xr      + atomJ.inducedDipole[1]*yr      + atomJ.inducedDipole[2]*zr;
        float preFactor3  = preFactor2*dukr;

        fields[0].x       = preFactor3*xr + preFactor1*atomJ.inducedDipole[0];
        fields[1].x       = preFactor3*yr + preFactor1*atomJ.inducedDipole[1];
        fields[2].x       = preFactor3*zr + preFactor1*atomJ.inducedDipole[2];


        float duir        = atomI.inducedDipole[0]*xr      + atomI.inducedDipole[1]*yr      + atomI.inducedDipole[2]*zr;
        preFactor3        = preFactor2*duir;

        fields[0].y       = preFactor3*xr + preFactor1*atomI.inducedDipole[0];
        fields[1].y       = preFactor3*yr + preFactor1*atomI.inducedDipole[1];
        fields[2].y       = preFactor3*zr + preFactor1*atomI.inducedDipole[2];


        float pukr        = atomJ.inducedDipolePolar[0]*xr + atomJ.inducedDipolePolar[1]*yr + atomJ.inducedDipolePolar[2]*zr;
        preFactor3        = preFactor2*pukr;

        fields[0].z       = preFactor3*xr + preFactor1*atomJ.inducedDipolePolar[0];
        fields[1].z       = preFactor3*yr + preFactor1*atomJ.inducedDipolePolar[1];
        fields[2].z       = preFactor3*zr + preFactor1*atomJ.inducedDipolePolar[2];


        float puir        = atomI.inducedDipolePolar[0]*xr + atomI.inducedDipolePolar[1]*yr + atomI.inducedDipolePolar[2]*zr;
        preFactor3        = preFactor2*puir;
        fields[0].w       = preFactor3*xr + preFactor1*atomI.inducedDipolePolar[0];
        fields[1].w       = preFactor3*yr + preFactor1*atomI.inducedDipolePolar[1];
        fields[2].w       = preFactor3*zr + preFactor1*atomI.inducedDipolePolar[2];

    } else {

        fields[0].x       = 0.0f;
        fields[0].y       = 0.0f;
        fields[0].z       = 0.0f;
        fields[0].w       = 0.0f;
    
        fields[1].x       = 0.0f;
        fields[1].y       = 0.0f;
        fields[1].z       = 0.0f;
        fields[1].w       = 0.0f;
    
        fields[2].x       = 0.0f;
        fields[2].y       = 0.0f;
        fields[2].z       = 0.0f;
        fields[2].w       = 0.0f;
    }
}

// Include versions of the kernels for N^2 calculations.

#define METHOD_NAME(a, b) a##Cutoff##b
#include "kCalculateAmoebaCudaPmeMutualInducedField.h"
#define USE_OUTPUT_BUFFER_PER_WARP
#undef METHOD_NAME
#define METHOD_NAME(a, b) a##CutoffByWarp##b
#include "kCalculateAmoebaCudaPmeMutualInducedField.h"

__global__ 
#if (__CUDA_ARCH__ >= 200)
__launch_bounds__(GF1XX_THREADS_PER_BLOCK, 1)
#elif (__CUDA_ARCH__ >= 120)
__launch_bounds__(GT2XX_THREADS_PER_BLOCK, 1)
#else
__launch_bounds__(G8X_THREADS_PER_BLOCK, 1)
#endif
static void kInitializeMutualInducedField_kernel(
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
static void kReduceMutualInducedFieldDelta_kernel(int numberOfEntries, float* arrayOfDeltas1, float* arrayOfDeltas2, float* epsilon )
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
static void kSorUpdateMutualInducedField_kernel(
                   float* polarizability,
                   float* inducedDipole, float* inducedDipoleP,
                   float* fixedEField,   float* fixedEFieldP,
                   float* matrixProduct, float* matrixProductP )
{

    int pos                        = blockIdx.x*blockDim.x + threadIdx.x;
    const float term               = (4.0f/3.0f)*(cSim.alphaEwald*cSim.alphaEwald*cSim.alphaEwald)/cAmoebaSim.sqrtPi;
    const float polarSOR           = 0.55f;

    while( pos < 3*cSim.atoms )
    {   

        float previousDipole           = inducedDipole[pos];
        float previousDipoleP          = inducedDipoleP[pos];
    
        // add self terms to fields
    
        float mProd                    = matrixProduct[pos];
        float mProdP                   = matrixProductP[pos];

        mProd                         += term*previousDipole;
        mProdP                        += term*previousDipoleP;
    
        float inducedDipoleI           = fixedEField[pos]     + polarizability[pos]*mProd;
        float inducedDipoleIP          = fixedEFieldP[pos]    + polarizability[pos]*mProdP;
    
        inducedDipole[pos]             = previousDipole   + polarSOR*( inducedDipoleI   - previousDipole  );   
        inducedDipoleP[pos]            = previousDipoleP  + polarSOR*( inducedDipoleIP  - previousDipoleP );
    
        matrixProduct[pos]             = ( inducedDipole[pos]  - previousDipole  )*( inducedDipole[pos]  - previousDipole  );
        matrixProductP[pos]            = ( inducedDipoleP[pos] - previousDipoleP )*( inducedDipoleP[pos] - previousDipoleP );

        pos                           += blockDim.x*gridDim.x;
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
    LAUNCHERROR("kReducePmeMI_Fields1");

    kReduceFields_kernel<<<gpu->sim.nonbond_blocks, gpu->sim.bsf_reduce_threads_per_block>>>(
                               gpu->sim.paddedNumberOfAtoms*3, gpu->sim.outputBuffers,
                               amoebaGpu->psWorkArray_3_2->_pDevData, outputPolarArray->_pDevData, 0 );
    LAUNCHERROR("kReducePmeMI_Fields2");
}

/**---------------------------------------------------------------------------------------

   Compute mutual induce field

   @param amoebaGpu        amoebaGpu context

   --------------------------------------------------------------------------------------- */

static void cudaComputeAmoebaPmeMutualInducedFieldMatrixMultiply( amoebaGpuContext amoebaGpu,
                                                                  CUDAStream<float>* outputArray, CUDAStream<float>* outputPolarArray )
{
  
    static unsigned int threadsPerBlock  = 0;
    gpuContext gpu                       = amoebaGpu->gpuContext;

    kClearFields_3( amoebaGpu, 2 );

    // on first pass, set threads/block

    if( threadsPerBlock == 0 ){  
        unsigned int maxThreads;
        if (gpu->sm_version >= SM_20)
            maxThreads = 384; 
        else if (gpu->sm_version >= SM_12)
            maxThreads = 128; 
        else
            maxThreads = 64; 
        threadsPerBlock = std::min(getThreadsPerBlock(amoebaGpu, sizeof(MutualInducedParticle), gpu->sharedMemoryPerBlock ), maxThreads);
    }    

    if (gpu->bOutputBufferPerWarp){

        kCalculateAmoebaPmeMutualInducedFieldCutoffByWarp_kernel<<<gpu->sim.nonbond_blocks, threadsPerBlock, sizeof(MutualInducedParticle)*threadsPerBlock>>>(
                                                                 gpu->sim.pInteractingWorkUnit,
                                                                 amoebaGpu->psWorkArray_3_1->_pDevData,
                                                                 amoebaGpu->psWorkArray_3_2->_pDevData );

    } else {

        kCalculateAmoebaPmeMutualInducedFieldCutoff_kernel<<<gpu->sim.nonbond_blocks, threadsPerBlock, sizeof(MutualInducedParticle)*threadsPerBlock>>>(
                                                                 gpu->sim.pInteractingWorkUnit,
                                                                 amoebaGpu->psWorkArray_3_1->_pDevData,
                                                                 amoebaGpu->psWorkArray_3_2->_pDevData );

    }
    LAUNCHERROR("kCalculateAmoebaPmeMutualInducedField");

    kReduceMutualInducedFields( amoebaGpu, outputArray, outputPolarArray );

}

/**---------------------------------------------------------------------------------------

   Compute mutual induce field

   @param amoebaGpu        amoebaGpu context

   --------------------------------------------------------------------------------------- */

static void cudaComputeAmoebaPmeMutualInducedFieldBySOR( amoebaGpuContext amoebaGpu )
{
  
   // ---------------------------------------------------------------------------------------

    int done;
    int iteration;

    gpuContext gpu    = amoebaGpu->gpuContext;

    // ---------------------------------------------------------------------------------------

    // set  E_Field & E_FieldPolar] to [ E_Field & E_FieldPolar]*Polarizability
    // initialize [ InducedDipole & InducedDipolePolar ] to [ E_Field & E_FieldPolar]*Polarizability

    kInitializeMutualInducedField_kernel<<< gpu->sim.blocks, gpu->sim.threads_per_block >>>(
         gpu->natoms,
         amoebaGpu->psE_Field->_pDevData,
         amoebaGpu->psE_FieldPolar->_pDevData,
         amoebaGpu->psPolarizability->_pDevData );
    LAUNCHERROR("AmoebaPmeMutualInducedFieldSetup");  

    cudaMemcpy( amoebaGpu->psInducedDipole->_pDevData,        amoebaGpu->psE_Field->_pDevData,       3*gpu->sim.paddedNumberOfAtoms*sizeof( float ), cudaMemcpyDeviceToDevice );
    cudaMemcpy( amoebaGpu->psInducedDipolePolar->_pDevData,   amoebaGpu->psE_FieldPolar->_pDevData,  3*gpu->sim.paddedNumberOfAtoms*sizeof( float ), cudaMemcpyDeviceToDevice );

    // if polarization type is direct, set flags signalling done and return

    if( amoebaGpu->amoebaSim.polarizationType )
    {
        amoebaGpu->mutualInducedDone          = 1;
        amoebaGpu->mutualInducedConverged     = 1;
        kCalculateAmoebaPMEInducedDipoleField( amoebaGpu );
        return;
    }

    // ---------------------------------------------------------------------------------------
 
    done      = 0;
    iteration = 1;

    while( !done ){

        //  apply SOR

        cudaComputeAmoebaPmeMutualInducedFieldMatrixMultiply( amoebaGpu, amoebaGpu->psWorkVector[0],  amoebaGpu->psWorkVector[1] );
        kCalculateAmoebaPMEInducedDipoleField( amoebaGpu );

        // post matrix multiply

        kSorUpdateMutualInducedField_kernel<<< gpu->sim.blocks, gpu->sim.threads_per_block >>>(
           amoebaGpu->psPolarizability->_pDevData,
           amoebaGpu->psInducedDipole->_pDevData, amoebaGpu->psInducedDipolePolar->_pDevData,
           amoebaGpu->psE_Field->_pDevData,       amoebaGpu->psE_FieldPolar->_pDevData,
           amoebaGpu->psWorkVector[0]->_pDevData, amoebaGpu->psWorkVector[1]->_pDevData );
        LAUNCHERROR("kSorUpdatePmeMutualInducedField");  

        // get total epsilon -- performing sums on gpu

        kReduceMutualInducedFieldDelta_kernel<<<1, amoebaGpu->epsilonThreadsPerBlock, 2*sizeof(float)*amoebaGpu->epsilonThreadsPerBlock>>>(
           3*gpu->natoms, amoebaGpu->psWorkVector[0]->_pDevData, amoebaGpu->psWorkVector[1]->_pDevData,
           amoebaGpu->psCurrentEpsilon->_pDevData );
        LAUNCHERROR("kReducePmeMutualInducedFieldDelta");

        // Debye=48.033324f
        amoebaGpu->psCurrentEpsilon->Download();
        float currentEpsilon                     = amoebaGpu->psCurrentEpsilon->_pSysData[0];
        amoebaGpu->mutualInducedCurrentEpsilon   = currentEpsilon;

        if( iteration > amoebaGpu->mutualInducedMaxIterations || amoebaGpu->mutualInducedCurrentEpsilon < amoebaGpu->mutualInducedTargetEpsilon ){ 
            done = 1;
        }

        // throw exception if nan detected

        if( amoebaGpu->mutualInducedCurrentEpsilon != amoebaGpu->mutualInducedCurrentEpsilon ){
            throw OpenMM::OpenMMException("PME induced dipole calculation detected nans." );
        }

        iteration++;
    }

    amoebaGpu->mutualInducedDone             = done;
    amoebaGpu->mutualInducedConverged        = ( !done || iteration > amoebaGpu->mutualInducedMaxIterations ) ? 0 : 1;

}

void cudaComputeAmoebaPmeMutualInducedField( amoebaGpuContext amoebaGpu )
{
    if( amoebaGpu->mutualInducedIterativeMethod == 0 ){
        cudaComputeAmoebaPmeMutualInducedFieldBySOR( amoebaGpu );
    }
}
