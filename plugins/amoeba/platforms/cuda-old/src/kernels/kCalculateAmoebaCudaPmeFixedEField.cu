
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

static __constant__ cudaGmxSimulation cSim;
static __constant__ cudaAmoebaGmxSimulation cAmoebaSim;

void SetCalculateAmoebaCudaPmeFixedEFieldSim(amoebaGpuContext amoebaGpu)
{
    cudaError_t status;
    gpuContext gpu = amoebaGpu->gpuContext;
    status         = cudaMemcpyToSymbol(cSim, &gpu->sim, sizeof(cudaGmxSimulation));    
    RTERROR(status, "SetCalculateAmoebaCudaPmeFixedEFieldSim: cudaMemcpyToSymbol: SetSim copy to cSim failed");
    status         = cudaMemcpyToSymbol(cAmoebaSim, &amoebaGpu->amoebaSim, sizeof(cudaAmoebaGmxSimulation));    
    RTERROR(status, "SetCalculateAmoebaCudaPmeFixedEFieldSim: cudaMemcpyToSymbol: SetSim copy to cAmoebaSim failed");
}

void GetCalculateAmoebaCudaPmeFixedEFieldSim(amoebaGpuContext amoebaGpu)
{
    cudaError_t status;
    gpuContext gpu = amoebaGpu->gpuContext;
    status = cudaMemcpyFromSymbol(&gpu->sim, cSim, sizeof(cudaGmxSimulation));    
    RTERROR(status, "GetCalculateAmoebaCudaPmeFixedEFieldSim: cudaMemcpyFromSymbol: SetSim copy from cSim failed");
    status = cudaMemcpyFromSymbol(&amoebaGpu->amoebaSim, cAmoebaSim, sizeof(cudaAmoebaGmxSimulation));         
    RTERROR(status, "GetCalculateAmoebaCudaPmeFixedEFieldSim: cudaMemcpyFromSymbol: SetSim copy from cAmoebaSim failed");
}

__global__ 
#if (__CUDA_ARCH__ >= 200)
__launch_bounds__(GF1XX_THREADS_PER_BLOCK, 1)
#elif (__CUDA_ARCH__ >= 120)
__launch_bounds__(GT2XX_THREADS_PER_BLOCK, 1)
#else
__launch_bounds__(G8X_THREADS_PER_BLOCK, 1)
#endif
static void kReducePmeEFieldPolar_kernel( unsigned int fieldComponents, unsigned int outputBuffers, float* EFieldReciprocal,  float* fieldIn, float* fieldOut )
{
    unsigned int pos = blockIdx.x * blockDim.x + threadIdx.x;

    // Reduce field

    const float term = (4.0f/3.0f)*(cSim.alphaEwald*cSim.alphaEwald*cSim.alphaEwald)/cAmoebaSim.sqrtPi;
    //const float term = 0.0f;
    while (pos < fieldComponents)
    {   

        // self-term included here

        float totalField = EFieldReciprocal[pos] + term*cAmoebaSim.pLabFrameDipole[pos];

        float* pFt       = fieldIn + pos;
        unsigned int i   = outputBuffers;
        while (i >= 4)
        {   
            totalField += pFt[0] + pFt[fieldComponents] + pFt[2*fieldComponents] + pFt[3*fieldComponents];
            pFt        += fieldComponents*4;
            i          -= 4;
        }   

        if (i >= 2)
        {   
            totalField += pFt[0] + pFt[fieldComponents];
            pFt        += fieldComponents*2;
            i          -= 2;
        }   

        if (i > 0)
        {   
            totalField += pFt[0];
        }   

        fieldOut[pos]   = totalField;
        pos            += gridDim.x * blockDim.x;
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
static void kReducePmeEField_kernel( unsigned int fieldComponents, unsigned int outputBuffers,  float* fieldIn, float* fieldOut )
{
    unsigned int pos = blockIdx.x * blockDim.x + threadIdx.x;

    // Reduce field

    const float term = (4.0f/3.0f)*(cSim.alphaEwald*cSim.alphaEwald*cSim.alphaEwald)/cAmoebaSim.sqrtPi;
    //const float term = 0.0;
    while (pos < fieldComponents)
    {   

        // self-term included here

        float totalField = term*cAmoebaSim.pLabFrameDipole[pos];

        float* pFt       = fieldIn + pos;
        unsigned int i   = outputBuffers;
        while (i >= 4)
        {   
            totalField += pFt[0] + pFt[fieldComponents] + pFt[2*fieldComponents] + pFt[3*fieldComponents];
            pFt        += fieldComponents*4;
            i          -= 4;
        }   

        if (i >= 2)
        {   
            totalField += pFt[0] + pFt[fieldComponents];
            pFt        += fieldComponents*2;
            i          -= 2;
        }   

        if (i > 0)
        {   
            totalField += pFt[0];
        }   

        fieldOut[pos]  += totalField;
        pos            += gridDim.x * blockDim.x;
    }   
}

// reduce psWorkArray_3_1 -> EField
// reduce psWorkArray_3_2 -> EFieldPolar

static void kReducePmeDirectE_Fields(amoebaGpuContext amoebaGpu )
{

    gpuContext gpu = amoebaGpu->gpuContext;

    // E_FieldPolar = E_Field (reciprocal) + E_FieldPolar (direct) + self

    kReducePmeEFieldPolar_kernel<<<gpu->sim.nonbond_blocks, gpu->sim.bsf_reduce_threads_per_block>>>(
                                   gpu->sim.paddedNumberOfAtoms*3, gpu->sim.outputBuffers,
                                   amoebaGpu->psE_Field->_pDevData, amoebaGpu->psWorkArray_3_2->_pDevData, amoebaGpu->psE_FieldPolar->_pDevData );
    LAUNCHERROR("kReducePmeE_Fields1");

    // E_Field = E_Field (reciprocal) + E_Field (direct) + self

    kReducePmeEField_kernel<<<gpu->sim.nonbond_blocks, gpu->sim.bsf_reduce_threads_per_block>>>(
                              gpu->sim.paddedNumberOfAtoms*3, gpu->sim.outputBuffers,
                              amoebaGpu->psWorkArray_3_1->_pDevData, amoebaGpu->psE_Field->_pDevData );
    LAUNCHERROR("kReducePmeE_Fields2");
}

// file includes FixedFieldParticle struct definition/load/unload struct and body kernel for fixed E-field

#undef GK
#undef INCLUDE_FIXED_FIELD_BUFFERS
#define INCLUDE_FIXED_FIELD_BUFFERS
#include "kCalculateAmoebaCudaFixedFieldParticle.h"
#undef INCLUDE_FIXED_FIELD_BUFFERS
__device__ void sumTempBuffer( FixedFieldParticle& atomI, FixedFieldParticle& atomJ ){
    atomI.tempBuffer[0]  += atomJ.tempBuffer[0];
    atomI.tempBuffer[1]  += atomJ.tempBuffer[1];
    atomI.tempBuffer[2]  += atomJ.tempBuffer[2];

    atomI.tempBufferP[0] += atomJ.tempBufferP[0];
    atomI.tempBufferP[1] += atomJ.tempBufferP[1];
    atomI.tempBufferP[2] += atomJ.tempBufferP[2];
}

__device__ void calculateFixedFieldRealSpacePairIxn_kernel( FixedFieldParticle& atomI, FixedFieldParticle& atomJ,
                                                            float dscale, float pscale, float4 fields[3]){

    // compute the real space portion of the Ewald summation
  
    float xr          = atomJ.x - atomI.x;
    float yr          = atomJ.y - atomI.y;
    float zr          = atomJ.z - atomI.z;

    // periodic boundary conditions

    xr               -= floorf(xr*cSim.invPeriodicBoxSizeX+0.5f)*cSim.periodicBoxSizeX;
    yr               -= floorf(yr*cSim.invPeriodicBoxSizeY+0.5f)*cSim.periodicBoxSizeY;
    zr               -= floorf(zr*cSim.invPeriodicBoxSizeZ+0.5f)*cSim.periodicBoxSizeZ;

    float r2          = xr*xr + yr*yr + zr*zr;
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

        alsq2n           *= alsq2;
        float bn3         = (5.0f*bn2+alsq2n*exp2a)/r2;

        // compute the error function scaled and unscaled terms

        float scale3      = 1.0f;
        float scale5      = 1.0f;
        float scale7      = 1.0f;
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
                scale7        = 1.0f - expdamp*(1.0f-damp+(0.6f*damp*damp));
            }
        }
        float dsc3        = dscale*scale3;
        float dsc5        = dscale*scale5;
        float dsc7        = dscale*scale7;

        float psc3        = pscale*scale3;
        float psc5        = pscale*scale5;
        float psc7        = pscale*scale7;

        float r3          = (r*r2);
        float r5          = (r3*r2);
        float r7          = (r5*r2);
        float drr3        = (1.0f-dsc3)/r3;
        float drr5        = 3.0f * (1.0f-dsc5)/r5;
        float drr7        = 15.0f * (1.0f-dsc7)/r7;

        float prr3        = (1.0f-psc3) / r3;
        float prr5        = 3.0f *(1.0f-psc5)/r5;
        float prr7        = 15.0f*(1.0f-psc7)/r7;

        float dir         = atomI.labFrameDipole_X*xr      + atomI.labFrameDipole_Y*yr      + atomI.labFrameDipole_Z*zr;

        float qix         = atomI.labFrameQuadrupole_XX*xr + atomI.labFrameQuadrupole_XY*yr + atomI.labFrameQuadrupole_XZ*zr;
        float qiy         = atomI.labFrameQuadrupole_XY*xr + atomI.labFrameQuadrupole_YY*yr + atomI.labFrameQuadrupole_YZ*zr;
        float qiz         = atomI.labFrameQuadrupole_XZ*xr + atomI.labFrameQuadrupole_YZ*yr + atomI.labFrameQuadrupole_ZZ*zr;

        float qir         = qix*xr + qiy*yr + qiz*zr;

        float dkr         = atomJ.labFrameDipole_X*xr      + atomJ.labFrameDipole_Y*yr      + atomJ.labFrameDipole_Z*zr;

        float qkx         = atomJ.labFrameQuadrupole_XX*xr + atomJ.labFrameQuadrupole_XY*yr + atomJ.labFrameQuadrupole_XZ*zr;
        float qky         = atomJ.labFrameQuadrupole_XY*xr + atomJ.labFrameQuadrupole_YY*yr + atomJ.labFrameQuadrupole_YZ*zr;
        float qkz         = atomJ.labFrameQuadrupole_XZ*xr + atomJ.labFrameQuadrupole_YZ*yr + atomJ.labFrameQuadrupole_ZZ*zr;

        float qkr         = qkx*xr + qky*yr + qkz*zr;

        float fim0        = -xr*(bn1*atomJ.q-bn2*dkr+bn3*qkr)    - bn1*atomJ.labFrameDipole_X  + 2.0f*bn2*qkx;
        float fim1        = -yr*(bn1*atomJ.q-bn2*dkr+bn3*qkr)    - bn1*atomJ.labFrameDipole_Y  + 2.0f*bn2*qky;
        float fim2        = -zr*(bn1*atomJ.q-bn2*dkr+bn3*qkr)    - bn1*atomJ.labFrameDipole_Z  + 2.0f*bn2*qkz;

        float fkm0        = xr*(bn1*atomI.q+bn2*dir+bn3*qir)     - bn1*atomI.labFrameDipole_X  - 2.0f*bn2*qix;
        float fkm1        = yr*(bn1*atomI.q+bn2*dir+bn3*qir)     - bn1*atomI.labFrameDipole_Y  - 2.0f*bn2*qiy;
        float fkm2        = zr*(bn1*atomI.q+bn2*dir+bn3*qir)     - bn1*atomI.labFrameDipole_Z  - 2.0f*bn2*qiz;

        float fid0        = -xr*(drr3*atomJ.q-drr5*dkr+drr7*qkr) - drr3*atomJ.labFrameDipole_X + 2.0f*drr5*qkx;
        float fid1        = -yr*(drr3*atomJ.q-drr5*dkr+drr7*qkr) - drr3*atomJ.labFrameDipole_Y + 2.0f*drr5*qky;
        float fid2        = -zr*(drr3*atomJ.q-drr5*dkr+drr7*qkr) - drr3*atomJ.labFrameDipole_Z + 2.0f*drr5*qkz;

        float fkd0        = xr*(drr3*atomI.q+drr5*dir+drr7*qir)  - drr3*atomI.labFrameDipole_X - 2.0f*drr5*qix;
        float fkd1        = yr*(drr3*atomI.q+drr5*dir+drr7*qir)  - drr3*atomI.labFrameDipole_Y - 2.0f*drr5*qiy;
        float fkd2        = zr*(drr3*atomI.q+drr5*dir+drr7*qir)  - drr3*atomI.labFrameDipole_Z - 2.0f*drr5*qiz;

        float fip0        = -xr*(prr3*atomJ.q-prr5*dkr+prr7*qkr) - prr3*atomJ.labFrameDipole_X + 2.0f*prr5*qkx;
        float fip1        = -yr*(prr3*atomJ.q-prr5*dkr+prr7*qkr) - prr3*atomJ.labFrameDipole_Y + 2.0f*prr5*qky;
        float fip2        = -zr*(prr3*atomJ.q-prr5*dkr+prr7*qkr) - prr3*atomJ.labFrameDipole_Z + 2.0f*prr5*qkz;

        float fkp0        = xr*(prr3*atomI.q+prr5*dir+prr7*qir)  - prr3*atomI.labFrameDipole_X - 2.0f*prr5*qix;
        float fkp1        = yr*(prr3*atomI.q+prr5*dir+prr7*qir)  - prr3*atomI.labFrameDipole_Y - 2.0f*prr5*qiy;
        float fkp2        = zr*(prr3*atomI.q+prr5*dir+prr7*qir)  - prr3*atomI.labFrameDipole_Z - 2.0f*prr5*qiz;

        // increment the field at each site due to this interaction

        fields[0].x       = fim0 - fid0;
        fields[1].x       = fim1 - fid1;
        fields[2].x       = fim2 - fid2;

        fields[0].y       = fkm0 - fkd0;
        fields[1].y       = fkm1 - fkd1;
        fields[2].y       = fkm2 - fkd2;

        fields[0].z       = fim0 - fip0;
        fields[1].z       = fim1 - fip1;
        fields[2].z       = fim2 - fip2;

        fields[0].w       = fkm0 - fkp0;
        fields[1].w       = fkm1 - fkp1;
        fields[2].w       = fkm2 - fkp2;
 
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
#include "kCalculateAmoebaCudaPmeFixedEField.h"
#define USE_OUTPUT_BUFFER_PER_WARP
#undef METHOD_NAME
#define METHOD_NAME(a, b) a##CutoffByWarp##b
#include "kCalculateAmoebaCudaPmeFixedEField.h"

/**---------------------------------------------------------------------------------------

   Report whether a number is a nan or infinity

   @param number               number to test
   @return 1 if number is  nan or infinity; else return 0

   --------------------------------------------------------------------------------------- */

/**---------------------------------------------------------------------------------------

   Compute fixed electric field using PME

   @param amoebaGpu        amoebaGpu context

   --------------------------------------------------------------------------------------- */

static void cudaComputeAmoebaPmeDirectFixedEField( amoebaGpuContext amoebaGpu )
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
            maxThreads = 192;
        else
            maxThreads = 64;
        threadsPerBlock = std::min(getThreadsPerBlock(amoebaGpu, sizeof(FixedFieldParticle), gpu->sharedMemoryPerBlock ), maxThreads);
    }    

    if (gpu->bOutputBufferPerWarp){
        kCalculateAmoebaPmeDirectFixedE_FieldCutoffByWarp_kernel<<<gpu->sim.nonbond_blocks, threadsPerBlock, sizeof(FixedFieldParticle)*threadsPerBlock>>>(
                                                                           gpu->sim.pInteractingWorkUnit,
                                                                           amoebaGpu->psWorkArray_3_1->_pDevData,
                                                                           amoebaGpu->psWorkArray_3_2->_pDevData );
    } else {
        kCalculateAmoebaPmeDirectFixedE_FieldCutoff_kernel<<<gpu->sim.nonbond_blocks, threadsPerBlock, sizeof(FixedFieldParticle)*threadsPerBlock>>>(
                                                                           gpu->sim.pInteractingWorkUnit,
                                                                           amoebaGpu->psWorkArray_3_1->_pDevData,
                                                                           amoebaGpu->psWorkArray_3_2->_pDevData );
    }
    LAUNCHERROR("kCalculateAmoebaPmeDirectFixedE_Field_kernel");

    kReducePmeDirectE_Fields( amoebaGpu );

}

void cudaComputeAmoebaPmeFixedEField( amoebaGpuContext amoebaGpu )
{

    kCalculateAmoebaPMEFixedMultipoles( amoebaGpu );
    cudaComputeAmoebaPmeDirectFixedEField( amoebaGpu );

}
