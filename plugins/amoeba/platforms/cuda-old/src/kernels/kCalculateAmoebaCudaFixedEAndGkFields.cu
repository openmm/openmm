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

void SetCalculateAmoebaCudaFixedEAndGKFieldsSim(amoebaGpuContext amoebaGpu)
{
    cudaError_t status;
    gpuContext gpu = amoebaGpu->gpuContext;
    status         = cudaMemcpyToSymbol(cSim, &gpu->sim, sizeof(cudaGmxSimulation));    
    RTERROR(status, "SetCalculateAmoebaCudaFixedEAndGKFieldSim: cudaMemcpyToSymbol: SetSim copy to cSim failed");
    status         = cudaMemcpyToSymbol(cAmoebaSim, &amoebaGpu->amoebaSim, sizeof(cudaAmoebaGmxSimulation));    
    RTERROR(status, "SetCalculateAmoebaCudaFixedEAndGKFieldSim: cudaMemcpyToSymbol: SetSim copy to cAmoebaSim failed");
}

void GetCalculateAmoebaCudaFixedEAndGKFieldSim(amoebaGpuContext amoebaGpu)
{
    cudaError_t status;
    gpuContext gpu = amoebaGpu->gpuContext;
    status = cudaMemcpyFromSymbol(&gpu->sim, cSim, sizeof(cudaGmxSimulation));    
    RTERROR(status, "GetCalculateAmoebaCudaFixedEAndGKFieldSim: cudaMemcpyFromSymbol: SetSim copy from cSim failed");
    status = cudaMemcpyFromSymbol(&amoebaGpu->amoebaSim, cAmoebaSim, sizeof(cudaAmoebaGmxSimulation));         
    RTERROR(status, "GetCalculateAmoebaCudaFixedEAndGKFieldSim: cudaMemcpyFromSymbol: SetSim copy from cAmoebaSim failed");
}

// reduce psWorkArray_3_1 -> E_Field
// reduce psWorkArray_3_2 -> E_FieldPolar
// reduce psWorkArray_3_3 -> Gk_FieldPolar

static void kReduceEAndGkFields(amoebaGpuContext amoebaGpu )
{

    gpuContext gpu = amoebaGpu->gpuContext;
    kReduceFields_kernel<<<gpu->sim.nonbond_blocks, gpu->sim.bsf_reduce_threads_per_block>>>(
                               gpu->sim.paddedNumberOfAtoms*3, gpu->sim.outputBuffers,
                               amoebaGpu->psWorkArray_3_1->_pDevData, amoebaGpu->psE_Field->_pDevData, 0 );

    LAUNCHERROR("kReduceEAndGK_Fields1");
    kReduceFields_kernel<<<gpu->sim.nonbond_blocks, gpu->sim.bsf_reduce_threads_per_block>>>(
                               gpu->sim.paddedNumberOfAtoms*3, gpu->sim.outputBuffers,
                               amoebaGpu->psWorkArray_3_2->_pDevData, amoebaGpu->psE_FieldPolar->_pDevData, 0 );
    LAUNCHERROR("kReduceEAndGK_Fields2");

    kReduceFields_kernel<<<gpu->sim.nonbond_blocks, gpu->sim.bsf_reduce_threads_per_block>>>(
                               gpu->sim.paddedNumberOfAtoms*3, gpu->sim.outputBuffers,
                               amoebaGpu->psWorkArray_3_3->_pDevData, amoebaGpu->psGk_Field->_pDevData, 0 );
    LAUNCHERROR("kReduceEAndGK_Fields3");
}

// file includes FixedFieldParticle struct definition/load/unload struct and kernel body for fixed E-field

#define GK
#include "kCalculateAmoebaCudaFixedFieldParticle.h"
#undef GK

__device__ void calculateFixedGkFieldPairIxn_kernel( float4 atomCoordinatesI,       float4 atomCoordinatesJ,
                                                     float* labFrameDipoleI,        float* labFrameDipoleJ,
                                                     float* labFrameQuadrupoleI,    float* labFrameQuadrupoleJ,
                                                     float  rb2,
                                                     float outputField[2][3]
 ){
  
    float xi,yi,zi;
    float xr,yr,zr;
    float xr2,yr2,zr2;
    float ci,ck;
    float uxi,uyi,uzi;
    float uxk,uyk,uzk;
    float qxxi,qxyi,qxzi;
    float qyyi,qyzi,qzzi;
    float qxxk,qxyk,qxzk;
    float qyyk,qyzk,qzzk;
    float r2;
    float fc,fd,fq;
    float expterm;
    float gf,gf2,gf3,gf5;
    float gf7;
    float expc,dexpc;
    float expc1,expcdexpc;
    float a[4][4];
    float gc[5];
    float gux[11],guy[11],guz[11];
    float gqxx[5],gqxy[5];
    float gqxz[5],gqyy[5];
    float gqyz[5],gqzz[5];

    float gkc;

    gkc          = cAmoebaSim.gkc;

    fc           = cAmoebaSim.fc;
    fd           = cAmoebaSim.fd;
    fq           = cAmoebaSim.fq;

    xi           = atomCoordinatesI.x;
    yi           = atomCoordinatesI.y;
    zi           = atomCoordinatesI.z;
    ci           = atomCoordinatesI.w;

    uxi          = labFrameDipoleI[0];
    uyi          = labFrameDipoleI[1];
    uzi          = labFrameDipoleI[2];

    qxxi         = labFrameQuadrupoleI[0];
    qxyi         = labFrameQuadrupoleI[1];
    qxzi         = labFrameQuadrupoleI[2];
    qyyi         = labFrameQuadrupoleI[4];
    qyzi         = labFrameQuadrupoleI[5];
    qzzi         = labFrameQuadrupoleI[8];

    xr           = atomCoordinatesJ.x - xi;
    yr           = atomCoordinatesJ.y - yi;
    zr           = atomCoordinatesJ.z - zi;
    ck           = atomCoordinatesJ.w;

    xr2          = xr*xr;
    yr2          = yr*yr;
    zr2          = zr*zr;
    r2           = xr2 + yr2 + zr2;

    uxk          = labFrameDipoleJ[0];
    uyk          = labFrameDipoleJ[1];
    uzk          = labFrameDipoleJ[2];

    qxxk         = labFrameQuadrupoleJ[0];
    qxyk         = labFrameQuadrupoleJ[1];
    qxzk         = labFrameQuadrupoleJ[2];
    qyyk         = labFrameQuadrupoleJ[4];
    qyzk         = labFrameQuadrupoleJ[5];
    qzzk         = labFrameQuadrupoleJ[8];

    expterm      = expf(-r2/(gkc*rb2));
    expc         = expterm / gkc;
    dexpc        = -2.0f / (gkc*rb2);
    gf2          = 1.0f / (r2+rb2*expterm);
    gf           = sqrtf(gf2);
    gf3          = gf2 * gf;
    gf5          = gf3 * gf2;
    gf7          = gf5 * gf2;

    // reaction potential auxiliary terms

    a[0][0]      = gf;
    a[1][0]      = -gf3;
    a[2][0]      = 3.0f * gf5;
    a[3][0]      = -15.0f * gf7;

    // reaction potential gradient auxiliary terms

    expc1        = 1.0f - expc;
    a[0][1]      = expc1 * a[1][0];
    a[1][1]      = expc1 * a[2][0];
    a[2][1]      = expc1 * a[3][0];

    // dipole second reaction potential gradient auxiliary term

    expcdexpc    = -expc * dexpc;
    a[1][2]      = expc1*a[2][1] + expcdexpc*a[2][0];

    // multiply the auxillary terms by dielectric functions;

    a[0][1]      = fc * a[0][1];
    a[1][0]      = fd * a[1][0];
    a[1][1]      = fd * a[1][1];
    a[1][2]      = fd * a[1][2];
    a[2][0]      = fq * a[2][0];
    a[2][1]      = fq * a[2][1];

    // unweighted dipole reaction potential tensor

    gux[1]       = xr * a[1][0];
    guy[1]       = yr * a[1][0];
    guz[1]       = zr * a[1][0];

    // unweighted reaction potential gradient tensor

    gc[2]        = xr * a[0][1];
    gc[3]        = yr * a[0][1];
    gc[4]        = zr * a[0][1];
    gux[2]       = a[1][0] + xr2*a[1][1];
    gux[3]       = xr * yr * a[1][1];
    gux[4]       = xr * zr * a[1][1];
    guy[2]       = gux[3];
    guy[3]       = a[1][0] + yr2*a[1][1];
    guy[4]       = yr * zr * a[1][1];
    guz[2]       = gux[4];
    guz[3]       = guy[4];
    guz[4]       = a[1][0] + zr2*a[1][1];
    gqxx[2]      = xr * (2.0f*a[2][0]+xr2*a[2][1]);
    gqxx[3]      = yr * xr2*a[2][1];
    gqxx[4]      = zr * xr2*a[2][1];
    gqyy[2]      = xr * yr2*a[2][1];
    gqyy[3]      = yr * (2.0f*a[2][0]+yr2*a[2][1]);
    gqyy[4]      = zr * yr2 * a[2][1];
    gqzz[2]      = xr * zr2 * a[2][1];
    gqzz[3]      = yr * zr2 * a[2][1];
    gqzz[4]      = zr * (2.0f*a[2][0]+zr2*a[2][1]);
    gqxy[2]      = yr * (a[2][0]+xr2*a[2][1]);
    gqxy[3]      = xr * (a[2][0]+yr2*a[2][1]);
    gqxy[4]      = zr * xr * yr * a[2][1];
    gqxz[2]      = zr * (a[2][0]+xr2*a[2][1]);
    gqxz[3]      = gqxy[4];
    gqxz[4]      = xr * (a[2][0]+zr2*a[2][1]);
    gqyz[2]      = gqxy[4];
    gqyz[3]      = zr * (a[2][0]+yr2*a[2][1]);
    gqyz[4]      = yr * (a[2][0]+zr2*a[2][1]);

    // unweighted dipole second reaction potential gradient tensor

    gux[5]       = xr * (3.0f*a[1][1]+xr2*a[1][2]);
    gux[6]       = yr * (a[1][1]+xr2*a[1][2]);
    gux[7]       = zr * (a[1][1]+xr2*a[1][2]);
    gux[8]       = xr * (a[1][1]+yr2*a[1][2]);
    gux[9]       = zr * xr * yr * a[1][2];
    gux[10]      = xr * (a[1][1]+zr2*a[1][2]);
    guy[5]       = yr * (a[1][1]+xr2*a[1][2]);
    guy[6]       = xr * (a[1][1]+yr2*a[1][2]);
    guy[7]       = gux[9];
    guy[8]       = yr * (3.0f*a[1][1]+yr2*a[1][2]);
    guy[9]       = zr * (a[1][1]+yr2*a[1][2]);
    guy[10]      = yr * (a[1][1]+zr2*a[1][2]);
    guz[5]       = zr * (a[1][1]+xr2*a[1][2]);
    guz[6]       = gux[9];
    guz[7]       = xr * (a[1][1]+zr2*a[1][2]);
    guz[8]       = zr * (a[1][1]+yr2*a[1][2]);
    guz[9]       = yr * (a[1][1]+zr2*a[1][2]);
    guz[10]      = zr * (3.0f*a[1][1]+zr2*a[1][2]);

    // generalized Kirkwood permanent reaction field

    outputField[0][0] = uxk*gux[2] + uyk*gux[3] + uzk*gux[4]
                                   + 0.5f * (ck*gux[1] + qxxk*gux[5]
                                   + qyyk*gux[8] + qzzk*gux[10]
                                   + 2.0f*(qxyk*gux[6]+qxzk*gux[7]
                                   + qyzk*gux[9]))
                                   + 0.5f * (ck*gc[2] + qxxk*gqxx[2]
                                   + qyyk*gqyy[2] + qzzk*gqzz[2]
                                   + 2.0f*(qxyk*gqxy[2]+qxzk*gqxz[2]
                                   + qyzk*gqyz[2]));

    outputField[0][1] = uxk*guy[2] + uyk*guy[3] + uzk*guy[4]
                                   + 0.5f * (ck*guy[1] + qxxk*guy[5]
                                   + qyyk*guy[8] + qzzk*guy[10]
                                   + 2.0f*(qxyk*guy[6]+qxzk*guy[7]
                                   + qyzk*guy[9]))
                                   + 0.5f * (ck*gc[3] + qxxk*gqxx[3]
                                   + qyyk*gqyy[3] + qzzk*gqzz[3]
                                   + 2.0f*(qxyk*gqxy[3]+qxzk*gqxz[3]
                                   + qyzk*gqyz[3]));

    outputField[0][2] = uxk*guz[2] + uyk*guz[3] + uzk*guz[4]
                                   + 0.5f * (ck*guz[1] + qxxk*guz[5]
                                   + qyyk*guz[8] + qzzk*guz[10]
                                   + 2.0f*(qxyk*guz[6]+qxzk*guz[7]
                                   + qyzk*guz[9]))
                                   + 0.5f * (ck*gc[4] + qxxk*gqxx[4]
                                   + qyyk*gqyy[4] + qzzk*gqzz[4]
                                   + 2.0f*(qxyk*gqxy[4]+qxzk*gqxz[4]
                                   + qyzk*gqyz[4]));

    outputField[1][0] = uxi*gux[2] + uyi*gux[3] + uzi*gux[4]
                                   - 0.5f * (ci*gux[1] + qxxi*gux[5]
                                   + qyyi*gux[8] + qzzi*gux[10]
                                   + 2.0f*(qxyi*gux[6]+qxzi*gux[7]
                                   + qyzi*gux[9]))
                                   - 0.5f * (ci*gc[2] + qxxi*gqxx[2]
                                   + qyyi*gqyy[2] + qzzi*gqzz[2]
                                   + 2.0f*(qxyi*gqxy[2]+qxzi*gqxz[2]
                                   + qyzi*gqyz[2]));

    outputField[1][1] = uxi*guy[2] + uyi*guy[3] + uzi*guy[4]
                                   - 0.5f * (ci*guy[1] + qxxi*guy[5]
                                   + qyyi*guy[8] + qzzi*guy[10]
                                   + 2.0f*(qxyi*guy[6]+qxzi*guy[7]
                                   + qyzi*guy[9]))
                                   - 0.5f * (ci*gc[3]      + qxxi*gqxx[3]
                                   + qyyi*gqyy[3] + qzzi*gqzz[3]
                                   + 2.0f*(qxyi*gqxy[3]+qxzi*gqxz[3]
                                   + qyzi*gqyz[3]));

    outputField[1][2] = uxi*guz[2] + uyi*guz[3] + uzi*guz[4]
                                   - 0.5f * (ci*guz[1] + qxxi*guz[5]
                                   + qyyi*guz[8] + qzzi*guz[10]
                                   + 2.0f*(qxyi*guz[6]+qxzi*guz[7]
                                   + qyzi*guz[9]))
                                   - 0.5f * (ci*gc[4] + qxxi*gqxx[4]
                                   + qyyi*gqyy[4] + qzzi*gqzz[4]
                                   + 2.0f*(qxyi*gqxy[4]+qxzi*gqxz[4]
                                   + qyzi*gqyz[4]));

}

// Include versions of the kernels for N^2 calculations.

#define METHOD_NAME(a, b) a##N2##b
#include "kCalculateAmoebaCudaFixedEAndGkFields.h"
#define USE_OUTPUT_BUFFER_PER_WARP
#undef METHOD_NAME
#define METHOD_NAME(a, b) a##N2ByWarp##b
#include "kCalculateAmoebaCudaFixedEAndGkFields.h"

/**---------------------------------------------------------------------------------------

   Compute fixed electric field

   @param amoebaGpu        amoebaGpu context
   @param gpu              OpenMM gpu Cuda context

   --------------------------------------------------------------------------------------- */

void cudaComputeAmoebaFixedEAndGkFields( amoebaGpuContext amoebaGpu )
{
  
   // ---------------------------------------------------------------------------------------

   // ---------------------------------------------------------------------------------------

    gpuContext gpu                             = amoebaGpu->gpuContext;

    // on first pass, set threads/block

    static unsigned int threadsPerBlock        = 0;
    if( threadsPerBlock == 0 ){
        unsigned int maxThreads;
        if (gpu->sm_version >= SM_20)
            maxThreads = 256;
        else if (gpu->sm_version >= SM_12)
            maxThreads = 128;
        else
            maxThreads = 64;
        threadsPerBlock = std::min(getThreadsPerBlock(amoebaGpu, sizeof(FixedFieldParticle), gpu->sharedMemoryPerBlock ), maxThreads);
    }

    kClearFields_3( amoebaGpu, 3 );

    if (gpu->bOutputBufferPerWarp){
        kCalculateAmoebaFixedEAndGkFieldN2ByWarp_kernel<<<gpu->sim.nonbond_blocks, threadsPerBlock, sizeof(FixedFieldParticle)*threadsPerBlock>>>(
                                                                           gpu->psWorkUnit->_pDevData,
                                                                           amoebaGpu->psWorkArray_3_1->_pDevData, amoebaGpu->psWorkArray_3_2->_pDevData,
                                                                           amoebaGpu->psWorkArray_3_3->_pDevData );

    } else {

        kCalculateAmoebaFixedEAndGkFieldN2_kernel<<<gpu->sim.nonbond_blocks, threadsPerBlock, sizeof(FixedFieldParticle)*threadsPerBlock>>>(
                                                          gpu->psWorkUnit->_pDevData,
                                                          amoebaGpu->psWorkArray_3_1->_pDevData, amoebaGpu->psWorkArray_3_2->_pDevData,
                                                          amoebaGpu->psWorkArray_3_3->_pDevData );
    }
    LAUNCHERROR("kCalculateAmoebaFixedEAndGkFieldN2_kernel");

    kReduceEAndGkFields( amoebaGpu );
 
}
