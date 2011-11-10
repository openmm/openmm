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

#include "GpuFreeEnergyCudaKernels.h"
#include "freeEnergyGpuTypes.h"
#include <cudatypes.h>
#include "kSoftcoreLJ.h"

#define PARAMETER_PRINT 0
#define MAX_PARAMETER_PRINT 10

static __constant__ cudaGmxSimulation cSim;
static __constant__ cudaFreeEnergyGmxSimulation feSim;

void SetCalculateLocalSoftcoreGpuSim( freeEnergyGpuContext gpu)
{
    cudaError_t status;
    status = cudaMemcpyToSymbol(cSim, &gpu->gpuContext->sim, sizeof(cudaGmxSimulation));     
    RTERROR(status, "cudaMemcpyToSymbol: SetCalculateLocalSoftcoreGpuSim copy to cSim failed");

    status = cudaMemcpyToSymbol(feSim, &gpu->freeEnergySim, sizeof(cudaFreeEnergyGmxSimulation));     
    RTERROR(status, "cudaMemcpyToSymbol: SetCalculateLocalSoftcoreGpuSim copy to feSim failed");

}

extern "C"
void gpuSetLJ14SoftcoreParameters( freeEnergyGpuContext gpu, float epsfac, const std::vector<int>& atom1, const std::vector<int>& atom2,
                                   const std::vector<float>& c6, const std::vector<float>& c12, const std::vector<float>& qProd,
                                   const std::vector<float>& softcoreLJLambdaArray ){

    unsigned int LJ14s                                  = atom1.size();
    gpu->freeEnergySim.LJ14_count                       = LJ14s;

    gpu->psLJ14ID                                       = new CUDAStream<int4>(LJ14s, 1, "LJ14SoftcoreID");
    CUDAStream<int4>* psLJ14ID                          = gpu->psLJ14ID;
    gpu->freeEnergySim.pLJ14ID                          = psLJ14ID->_pDevData;

    gpu->psLJ14Parameter                                = new CUDAStream<float4>(LJ14s, 1, "LJ14SoftcoreParameter");
    CUDAStream<float4>* psLJ14Parameter                 = gpu->psLJ14Parameter;
    gpu->freeEnergySim.pLJ14Parameter  = psLJ14Parameter->_pDevData;

    std::vector<int> outputBufferCounter( gpu->gpuContext->sim.atoms, 0 );

    for( int ii = 0; ii < LJ14s; ii++ ){

        (*psLJ14ID)[ii].x          = atom1[ii];
        (*psLJ14ID)[ii].y          = atom2[ii];
        (*psLJ14ID)[ii].z          = outputBufferCounter[atom1[ii]]++;
        (*psLJ14ID)[ii].w          = outputBufferCounter[atom2[ii]]++;

        float p0, p1, p2, p3;
        if( c12[ii] == 0.0f ){
            p0 = 0.0f;
            p1 = 1.0f;

        } else {
            p0 = c6[ii] * c6[ii] / c12[ii];
            p1 = pow(c12[ii] / c6[ii], 1.0f / 6.0f);
        }

        p2                       = epsfac*qProd[ii];
        p3                       = softcoreLJLambdaArray[ii];

        (*psLJ14Parameter)[ii].x = p0;
        (*psLJ14Parameter)[ii].y = p1;
        (*psLJ14Parameter)[ii].z = p2;
        (*psLJ14Parameter)[ii].w = p3;
    }

    // logging info

    if( gpu->log ){
        (void) fprintf( gpu->log, "gpuSetLJ14SoftcoreParameters: number of 1-4 bonds=%5u\n", LJ14s );
#ifdef PARAMETER_PRINT
        unsigned int maxPrint = MAX_PARAMETER_PRINT;
        for( unsigned int ii = 0; ii < LJ14s; ii++ ){
            (void) fprintf( gpu->log, "    %5d [%5d %5d %5d %5d] %15.7e %15.7e %15.7e %15.7e\n",
                            ii, (*psLJ14ID)[ii].x, (*psLJ14ID)[ii].y, (*psLJ14ID)[ii].z, (*psLJ14ID)[ii].w, 
                            (*psLJ14Parameter)[ii].x, (*psLJ14Parameter)[ii].y,
                            (*psLJ14Parameter)[ii].z/epsfac, (*psLJ14Parameter)[ii].w );
            if( ii == maxPrint ){
                (void) fprintf( gpu->log, "\n" );
                ii = LJ14s - maxPrint;
                if( ii < maxPrint )ii = maxPrint;
            }    
        }    
        (void) fprintf( gpu->log, "\n" );
#endif
        (void) fflush( gpu->log );
    }    

    psLJ14ID->Upload();
    psLJ14Parameter->Upload();

    return;
}

#define DOT3(v1, v2) (v1.x * v2.x + v1.y * v2.y + v1.z * v2.z)

__global__ void kCalculateLocalSoftcoreForces_kernel()
{
    unsigned int pos = blockIdx.x * blockDim.x + threadIdx.x;

    float energy     = 0.0f;

    if (feSim.nonbondedMethod == NO_CUTOFF)
    {
        while (pos < feSim.LJ14_count)
        {
            int4 atom               = feSim.pLJ14ID[pos];
            float4 LJ14             = feSim.pLJ14Parameter[pos];

            float4 a1               = cSim.pPosq[atom.x];
            float4 a2               = cSim.pPosq[atom.y];

            float3 d;
            d.x                     = a1.x - a2.x;
            d.y                     = a1.y - a2.y;
            d.z                     = a1.z - a2.z;

            float r2                = DOT3(d, d);
            float inverseR          = 1.0f / sqrt(r2);
#ifdef USE_SOFTCORE_LJ
            float CDLJ_energy       = 0.0f;
            float dEdR              = getSoftCoreLJ( r2, LJ14.y, LJ14.x, LJ14.w, LJ14.w, &CDLJ_energy );
            energy                 += CDLJ_energy;
#else
            float sig2              = inverseR * LJ14.y;
            sig2                   *= sig2;
            float sig6              = sig2 * sig2 * sig2;
            float dEdR              = LJ14.x * (12.0f * sig6 - 6.0f) * sig6;
            energy                 += LJ14.x * (sig6 - 1.0f) * sig6;
#endif
            energy                 += LJ14.z * inverseR;
            dEdR                   += LJ14.z * inverseR;
            dEdR                   *= inverseR * inverseR;
            unsigned int offsetA    = atom.x + atom.z * cSim.stride;
            unsigned int offsetB    = atom.y + atom.w * cSim.stride;
            float4 forceA           = cSim.pForce4[offsetA];
            float4 forceB           = cSim.pForce4[offsetB];
            d.x                    *= dEdR;
            d.y                    *= dEdR;
            d.z                    *= dEdR;
            forceA.x               += d.x;
            forceA.y               += d.y;
            forceA.z               += d.z;
            forceB.x               -= d.x;
            forceB.y               -= d.y;
            forceB.z               -= d.z;
            cSim.pForce4[offsetA]   = forceA;
            cSim.pForce4[offsetB]   = forceB;
            pos                    += blockDim.x * gridDim.x;
        }

    } else if (feSim.nonbondedMethod == CUTOFF) {
        float LJ14_energy;
        while (pos < feSim.LJ14_count ){
            int4 atom               = feSim.pLJ14ID[pos];
            float4 LJ14             = feSim.pLJ14Parameter[pos];
            float4 a1               = cSim.pPosq[atom.x];
            float4 a2               = cSim.pPosq[atom.y];
            float3 d;
            d.x                     = a1.x - a2.x;
            d.y                     = a1.y - a2.y;
            d.z                     = a1.z - a2.z;
            float r2                = DOT3(d, d);
            float inverseR          = 1.0f / sqrt(r2);
#ifdef USE_SOFTCORE_LJ
            float dEdR              = getSoftCoreLJ( r2, LJ14.y, LJ14.x, LJ14.w, LJ14.w, &LJ14_energy);
#else
            float sig2              = inverseR * LJ14.y;
            sig2                   *= sig2;
            float sig6              = sig2 * sig2 * sig2;
            float dEdR              = LJ14.x * (12.0f * sig6 - 6.0f) * sig6;                
            LJ14_energy             = LJ14.x * (sig6 - 1.0f) * sig6;
#endif
            LJ14_energy            += LJ14.z * (inverseR + cSim.reactionFieldK * r2 - cSim.reactionFieldC);
            dEdR                   += LJ14.z * (inverseR - 2.0f * cSim.reactionFieldK * r2);
            dEdR                   *= inverseR * inverseR;
            if (r2 > feSim.nonbondedCutoffSqr)
            {                   
                dEdR = 0.0f;
                LJ14_energy = 0.0f;
            }
            energy                 += LJ14_energy;
 
            unsigned int offsetA    = atom.x + atom.z * cSim.stride;
            unsigned int offsetB    = atom.y + atom.w * cSim.stride;
            float4 forceA           = cSim.pForce4[offsetA];
            float4 forceB           = cSim.pForce4[offsetB];
            d.x                    *= dEdR;
            d.y                    *= dEdR;
            d.z                    *= dEdR;
            forceA.x               += d.x;
            forceA.y               += d.y;
            forceA.z               += d.z;
            forceB.x               -= d.x;
            forceB.y               -= d.y;
            forceB.z               -= d.z;
            cSim.pForce4[offsetA]   = forceA;
            cSim.pForce4[offsetB]   = forceB;
            pos                    += blockDim.x * gridDim.x;
        }
    } else if (feSim.nonbondedMethod == PERIODIC ){
        float LJ14_energy;
        while (pos < feSim.LJ14_count ){
            int4 atom               = feSim.pLJ14ID[pos];
            float4 LJ14             = feSim.pLJ14Parameter[pos];
            float4 a1               = cSim.pPosq[atom.x];
            float4 a2               = cSim.pPosq[atom.y];
            float3 d;
            d.x                     = a1.x - a2.x;
            d.y                     = a1.y - a2.y;
            d.z                     = a1.z - a2.z;
            d.x                    -= floor(d.x/cSim.periodicBoxSizeX+0.5f)*cSim.periodicBoxSizeX;
            d.y                    -= floor(d.y/cSim.periodicBoxSizeY+0.5f)*cSim.periodicBoxSizeY;
            d.z                    -= floor(d.z/cSim.periodicBoxSizeZ+0.5f)*cSim.periodicBoxSizeZ;
            float r2                = DOT3(d, d);
            float inverseR          = 1.0f / sqrt(r2);
#ifdef USE_SOFTCORE_LJ
            float dEdR              = getSoftCoreLJ( r2, LJ14.y, LJ14.x, LJ14.w, LJ14.w, &LJ14_energy);
#else
            float sig2              = inverseR * LJ14.y;
            sig2                   *= sig2;
            float sig6              = sig2 * sig2 * sig2;
            float dEdR              = LJ14.x * (12.0f * sig6 - 6.0f) * sig6;
            LJ14_energy             = LJ14.x * (sig6 - 1.0f) * sig6;
#endif
            LJ14_energy            += LJ14.z * (inverseR + cSim.reactionFieldK * r2 - cSim.reactionFieldC);

            dEdR                   += LJ14.z * (inverseR - 2.0f * cSim.reactionFieldK * r2);
            dEdR                   *= inverseR * inverseR;
            if (r2 > feSim.nonbondedCutoffSqr)
            {
                dEdR = 0.0f;
                LJ14_energy = 0.0f;
            }
            energy                 += LJ14_energy;

            unsigned int offsetA    = atom.x + atom.z * cSim.stride;
            unsigned int offsetB    = atom.y + atom.w * cSim.stride;
            float4 forceA           = cSim.pForce4[offsetA];
            float4 forceB           = cSim.pForce4[offsetB];
            d.x                    *= dEdR;
            d.y                    *= dEdR;
            d.z                    *= dEdR;
            forceA.x               += d.x;
            forceA.y               += d.y;
            forceA.z               += d.z;
            forceB.x               -= d.x;
            forceB.y               -= d.y;
            forceB.z               -= d.z;
            cSim.pForce4[offsetA]   = forceA;
            cSim.pForce4[offsetB]   = forceB;
            pos                    += blockDim.x * gridDim.x;
        }
    }
    cSim.pEnergy[blockIdx.x * blockDim.x + threadIdx.x] += energy;
}

void kCalculateLocalSoftcoreForces( freeEnergyGpuContext freeEnergyGpuContext )
{
    gpuContext gpu = freeEnergyGpuContext->gpuContext;
    kCalculateLocalSoftcoreForces_kernel<<<gpu->sim.blocks, gpu->sim.localForces_threads_per_block, gpu->sim.localForces_threads_per_block * sizeof(Vectors)>>>();
    LAUNCHERROR("kCalculateLocalSoftcoreForces");
}

