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

#include "amoebaScaleFactors.h"

__global__ 
#if (__CUDA_ARCH__ >= 200)
__launch_bounds__(512, 1)
#elif (__CUDA_ARCH__ >= 120)
__launch_bounds__(128, 1)
#else
__launch_bounds__(64, 1)
#endif
void METHOD_NAME(kCalculateAmoebaCudaElectrostatic, Forces_kernel)(
                            unsigned int* workUnit, float* outputTorque

#ifdef AMOEBA_DEBUG
                           , float4* debugArray, unsigned int targetAtom
#endif
){

    extern __shared__ ElectrostaticParticle sA[];

    unsigned int totalWarps      = gridDim.x*blockDim.x/GRID;
    unsigned int warp            = (blockIdx.x*blockDim.x+threadIdx.x)/GRID;
    unsigned int numWorkUnits    = cSim.pInteractionCount[0];
    unsigned int pos             = warp*numWorkUnits/totalWarps;
    unsigned int end             = (warp+1)*numWorkUnits/totalWarps;
    unsigned int lasty           = 0xFFFFFFFF;
    float totalEnergy            = 0.0f;     
    float conversionFactor       = (cAmoebaSim.electric/cAmoebaSim.dielec);
    float scalingFactors[LastScalingIndex];

    while (pos < end)
    {

        unsigned int x;
        unsigned int y;
        bool bExclusionFlag;

        // Extract cell coordinates

        decodeCell( workUnit[pos], &x, &y, &bExclusionFlag );

        unsigned int tgx              = threadIdx.x & (GRID - 1);
        unsigned int tbx              = threadIdx.x - tgx;
        unsigned int tj               = tgx;

        ElectrostaticParticle* psA    = &sA[tbx];
        unsigned int atomI            = x + tgx;
        ElectrostaticParticle localParticle;
        loadElectrostaticParticle( &localParticle, atomI );
        zeroElectrostaticParticle( &localParticle );

        scalingFactors[PScaleIndex]   = 1.0f;
        scalingFactors[DScaleIndex]   = 1.0f;
        scalingFactors[UScaleIndex]   = 1.0f;
        scalingFactors[MScaleIndex]   = 1.0f;

        if (x == y) // Handle diagonals uniquely at 50% efficiency
        {

            // load shared data

            loadElectrostaticParticle( &(sA[threadIdx.x]), atomI );
            unsigned int xi       = x >> GRIDBITS;
            unsigned int cell     = xi + xi*cSim.paddedNumberOfAtoms/GRID-xi*(xi+1)/2;
            int  dScaleMask       = cAmoebaSim.pD_ScaleIndices[cAmoebaSim.pScaleIndicesIndex[cell]+tgx];
            int2 pScaleMask       = cAmoebaSim.pP_ScaleIndices[cAmoebaSim.pScaleIndicesIndex[cell]+tgx];
            int2 mScaleMask       = cAmoebaSim.pM_ScaleIndices[cAmoebaSim.pScaleIndicesIndex[cell]+tgx];

            for (unsigned int j = 0; j < GRID; j++)
            {

                unsigned int atomJ = y + j;
                if( (atomI != atomJ) && (atomI < cSim.atoms) && (atomJ < cSim.atoms) ){

                    getMaskedDScaleFactor( j, dScaleMask, scalingFactors + DScaleIndex );
                    getMaskedPScaleFactor( j, pScaleMask, scalingFactors + PScaleIndex );
                    getMaskedMScaleFactor( j, mScaleMask, scalingFactors + MScaleIndex );
    
                    float force[3];
                    float energy;
                    calculateElectrostaticPairIxnF1_kernel( localParticle, psA[j], scalingFactors, &energy, force);

                    localParticle.force[0]            += force[0];
                    localParticle.force[1]            += force[1];
                    localParticle.force[2]            += force[2];
                    totalEnergy                       += 0.5f*energy;

                }
            }

            // Write results

            localParticle.force[0]  *= conversionFactor;
            localParticle.force[1]  *= conversionFactor;
            localParticle.force[2]  *= conversionFactor;

#ifdef USE_OUTPUT_BUFFER_PER_WARP
            unsigned int offset      = (x + tgx + warp*cSim.paddedNumberOfAtoms);
#else
            unsigned int offset      = (x + tgx + (x >> GRIDBITS) * cSim.paddedNumberOfAtoms);
#endif
            add3dArrayToFloat4( offset, localParticle.force, cSim.pForce4 );

            zeroElectrostaticParticle( &localParticle );
            for (unsigned int j = 0; j < GRID; j++)
            {

                unsigned int atomJ = y + j;
                if( (atomI != atomJ) && (atomI < cSim.atoms) && (atomJ < cSim.atoms) ){

                    getMaskedDScaleFactor( j, dScaleMask, scalingFactors + DScaleIndex );
                    getMaskedPScaleFactor( j, pScaleMask, scalingFactors + PScaleIndex );
                    getMaskedMScaleFactor( j, mScaleMask, scalingFactors + MScaleIndex );
    
                    float force[3];
                    calculateElectrostaticPairIxnT1_kernel( localParticle, psA[j], scalingFactors, force);
                    localParticle.force[0]  += force[0];
                    localParticle.force[1]  += force[1];
                    localParticle.force[2]  += force[2];

                }
            }

            localParticle.force[0] *= conversionFactor;
            localParticle.force[1] *= conversionFactor;
            localParticle.force[2] *= conversionFactor;

#ifdef USE_OUTPUT_BUFFER_PER_WARP
            offset                 = (x + tgx + warp*cSim.paddedNumberOfAtoms);
            add3dArray( 3*offset, localParticle.force, outputTorque );
#else
            offset                 = (x + tgx + (x >> GRIDBITS) * cSim.paddedNumberOfAtoms);
            load3dArray( 3*offset, localParticle.force, outputTorque );
#endif

        } else {

            // Read fixed atom data into registers and GRF

            if( lasty != y ){
               loadElectrostaticParticle( &(sA[threadIdx.x]), (y+tgx) );
            }

            zeroElectrostaticParticle( &(sA[threadIdx.x]) );
       
            int  dScaleMask;
            int2 pScaleMask;
            int2 mScaleMask;

            if( bExclusionFlag ){
                unsigned int xi   = x >> GRIDBITS;
                unsigned int yi   = y >> GRIDBITS;
                unsigned int cell = xi+yi*cSim.paddedNumberOfAtoms/GRID-yi*(yi+1)/2;
                dScaleMask        = cAmoebaSim.pD_ScaleIndices[cAmoebaSim.pScaleIndicesIndex[cell]+tgx];
                pScaleMask        = cAmoebaSim.pP_ScaleIndices[cAmoebaSim.pScaleIndicesIndex[cell]+tgx];
                mScaleMask        = cAmoebaSim.pM_ScaleIndices[cAmoebaSim.pScaleIndicesIndex[cell]+tgx];
            }
       
            for (unsigned int j = 0; j < GRID; j++){
       
                unsigned int atomJ = y + tj;
                if( (atomI < cSim.atoms) && (atomJ < cSim.atoms) ){

                    if( bExclusionFlag ){
                        getMaskedDScaleFactor( tj, dScaleMask, scalingFactors + DScaleIndex );
                        getMaskedPScaleFactor( tj, pScaleMask, scalingFactors + PScaleIndex );
                        getMaskedMScaleFactor( tj, mScaleMask, scalingFactors + MScaleIndex );
                    }
           
                    float force[3];
                    float energy;
                    calculateElectrostaticPairIxnF1_kernel( localParticle, psA[tj], scalingFactors, &energy, force);
    
                    totalEnergy                       += energy;
    
                    localParticle.force[0]            += force[0];
                    localParticle.force[1]            += force[1];
                    localParticle.force[2]            += force[2];
    
                    psA[tj].force[0]                  -= force[0];
                    psA[tj].force[1]                  -= force[1];
                    psA[tj].force[2]                  -= force[2];

                }

                tj = (tj + 1) & (GRID - 1);
            }

            // Write results

            localParticle.force[0]     *= conversionFactor;
            localParticle.force[1]     *= conversionFactor;
            localParticle.force[2]     *= conversionFactor;
 
            sA[threadIdx.x].force[0]   *= conversionFactor;
            sA[threadIdx.x].force[1]   *= conversionFactor;
            sA[threadIdx.x].force[2]   *= conversionFactor;

#ifdef USE_OUTPUT_BUFFER_PER_WARP

            unsigned int offset                 = (x + tgx + warp*cSim.paddedNumberOfAtoms);
            add3dArrayToFloat4(   offset, localParticle.force,   cSim.pForce4 );

            offset                              = (y + tgx + warp*cSim.paddedNumberOfAtoms);
            add3dArrayToFloat4(   offset, sA[threadIdx.x].force,   cSim.pForce4 );

#else
            unsigned int offset                 = (x + tgx + (y >> GRIDBITS) * cSim.paddedNumberOfAtoms);
            add3dArrayToFloat4(   offset, localParticle.force,   cSim.pForce4 );

            offset                              = (y + tgx + (x >> GRIDBITS) * cSim.paddedNumberOfAtoms);
            add3dArrayToFloat4( offset, sA[threadIdx.x].force,  cSim.pForce4 );

#endif
            zeroElectrostaticParticle( &(sA[threadIdx.x]) );
            zeroElectrostaticParticle( &localParticle );
            tj = tgx;
            for (unsigned int j = 0; j < GRID; j++){
       
                unsigned int atomJ = y + tj;
                if( (atomI < cSim.atoms) && (atomJ < cSim.atoms) ){

                    if( bExclusionFlag ){
                        getMaskedDScaleFactor( tj, dScaleMask, scalingFactors + DScaleIndex );
                        getMaskedPScaleFactor( tj, pScaleMask, scalingFactors + PScaleIndex );
                        getMaskedMScaleFactor( tj, mScaleMask, scalingFactors + MScaleIndex );
                    }
           
                    float force[3];
                    calculateElectrostaticPairIxnT1_kernel( localParticle, psA[tj], scalingFactors, force);
                    localParticle.force[0]           += force[0];
                    localParticle.force[1]           += force[1];
                    localParticle.force[2]           += force[2];
    
                    calculateElectrostaticPairIxnT3_kernel( localParticle, psA[tj], scalingFactors, force);
                    psA[tj].force[0]                 += force[0];
                    psA[tj].force[1]                 += force[1];
                    psA[tj].force[2]                 += force[2];

                }

                tj = (tj + 1) & (GRID - 1);
            }

            localParticle.force[0]    *= conversionFactor;
            localParticle.force[1]    *= conversionFactor;
            localParticle.force[2]    *= conversionFactor;
 
            sA[threadIdx.x].force[0]  *= conversionFactor;
            sA[threadIdx.x].force[1]  *= conversionFactor;
            sA[threadIdx.x].force[2]  *= conversionFactor;

#ifdef USE_OUTPUT_BUFFER_PER_WARP

            offset                 = (x + tgx + warp*cSim.paddedNumberOfAtoms);
            add3dArray( 3*offset, localParticle.force,  outputTorque );

            offset                              = (y + tgx + warp*cSim.paddedNumberOfAtoms);
            add3dArray( 3*offset, sA[threadIdx.x].force,  outputTorque );

#else
            offset                 = (x + tgx + (y >> GRIDBITS) * cSim.paddedNumberOfAtoms);
            load3dArray(         3*offset, localParticle.force, outputTorque );

            offset                              = (y + tgx + (x >> GRIDBITS) * cSim.paddedNumberOfAtoms);
            load3dArray(       3*offset, sA[threadIdx.x].force, outputTorque );

#endif
            lasty = y;
        }

        pos++;
    }

    cSim.pEnergy[blockIdx.x * blockDim.x + threadIdx.x] += (conversionFactor*totalEnergy);
}
