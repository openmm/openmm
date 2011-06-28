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
__launch_bounds__(384, 1)
#elif (__CUDA_ARCH__ >= 120)
__launch_bounds__(128, 1)
#else
__launch_bounds__(64, 1)
#endif
void METHOD_NAME(kCalculateAmoebaPmeDirectElectrostatic, Forces_kernel)( unsigned int* workUnit, float* outputTorque ){

    extern __shared__ PmeDirectElectrostaticParticle sA[];

    unsigned int totalWarps      = gridDim.x*blockDim.x/GRID;
    unsigned int warp            = (blockIdx.x*blockDim.x+threadIdx.x)/GRID;
    unsigned int numWorkUnits    = cSim.pInteractionCount[0];
    unsigned int pos             = warp*numWorkUnits/totalWarps;
    unsigned int end             = (warp+1)*numWorkUnits/totalWarps;
    unsigned int lasty           = 0xFFFFFFFF;
    float totalEnergy            = 0.0f;     

    float scalingFactors[LastScalingIndex];
    float conversionFactor       = (-cAmoebaSim.electric/cAmoebaSim.dielec);

    while (pos < end)
    {

        unsigned int x;
        unsigned int y;
        bool bExclusionFlag;
        int  dScaleMask;
        int2 pScaleMask;
        int2 mScaleMask;

        // Extract cell coordinates

        decodeCell( workUnit[pos], &x, &y, &bExclusionFlag );

        unsigned int tgx                                = threadIdx.x & (GRID - 1);
        unsigned int tbx                                = threadIdx.x - tgx;
        unsigned int tj                                 = tgx;

        PmeDirectElectrostaticParticle* psA             = &sA[tbx];
        unsigned int atomI                              = x + tgx;
        PmeDirectElectrostaticParticle localParticleI;
        loadPmeDirectElectrostaticParticle( atomI, &localParticleI );

        zeroPmeDirectElectrostaticParticle( &localParticleI );
        scalingFactors[UScaleIndex]                     = 1.0f;

        if (x == y) // Handle diagonals uniquely at 50% efficiency
        {

            // load shared data

            loadPmeDirectElectrostaticParticle( atomI, &(sA[threadIdx.x]) );

            if( atomI < cSim.atoms ){
                if (bExclusionFlag)
                {
                    unsigned int xi       = x >> GRIDBITS;
                    unsigned int cell     = xi + xi*cSim.paddedNumberOfAtoms/GRID-xi*(xi+1)/2;
                    dScaleMask            = cAmoebaSim.pD_ScaleIndices[cAmoebaSim.pScaleIndicesIndex[cell]+tgx];
                    pScaleMask            = cAmoebaSim.pP_ScaleIndices[cAmoebaSim.pScaleIndicesIndex[cell]+tgx];
                    mScaleMask            = cAmoebaSim.pM_ScaleIndices[cAmoebaSim.pScaleIndicesIndex[cell]+tgx];
                } else {
                    scalingFactors[DScaleIndex] = scalingFactors[PScaleIndex] = scalingFactors[MScaleIndex] = 1.0f;
                }
    
                for (unsigned int j = 0; j < GRID && (y+j) < cSim.atoms; j++)
                {
    
                    if( atomI != (y+j) )
                    {
                        if (bExclusionFlag)
                        {
                            getMaskedDScaleFactor( j, dScaleMask, scalingFactors + DScaleIndex );
                            getMaskedPScaleFactor( j, pScaleMask, scalingFactors + PScaleIndex );
                            getMaskedMScaleFactor( j, mScaleMask, scalingFactors + MScaleIndex );
                        }
                        calculatePmeDirectElectrostaticPairIxn_kernel( localParticleI, psA[j], bExclusionFlag, scalingFactors, 0.5f, &totalEnergy);
                    }
    
                } // end of j-loop
    
                // include self energy and self torque
    
                calculatePmeSelfTorqueElectrostaticPairIxn_kernel( localParticleI );
                calculatePmeSelfEnergyElectrostaticPairIxn_kernel( localParticleI, &totalEnergy );
    
                localParticleI.force[0]  *= conversionFactor;
                localParticleI.force[1]  *= conversionFactor;
                localParticleI.force[2]  *= conversionFactor;
    
                localParticleI.torque[0] *= -conversionFactor;
                localParticleI.torque[1] *= -conversionFactor;
                localParticleI.torque[2] *= -conversionFactor;
    
                // Write results
    
#ifdef USE_OUTPUT_BUFFER_PER_WARP
                unsigned int offset                 = (x + tgx + warp*cSim.paddedNumberOfAtoms);
                add3dArrayToFloat4( offset, localParticleI.force, cSim.pForce4 );
                add3dArray( 3*offset, localParticleI.torque, outputTorque );
#else
                unsigned int offset                 = (x + tgx + (x >> GRIDBITS) * cSim.paddedNumberOfAtoms);
                add3dArrayToFloat4( offset, localParticleI.force, cSim.pForce4 );
                load3dArray( 3*offset, localParticleI.torque, outputTorque );
#endif
            }

        } else {

            if (lasty != y) {
               loadPmeDirectElectrostaticParticle( (y+tgx), &(sA[threadIdx.x]) );
            }

            if (cSim.pInteractionFlag[pos] != 0 ) {

                zeroPmeDirectElectrostaticParticle( &(sA[threadIdx.x]) );
/*
                sA[threadIdx.x].force[0]     = 0.0f;
                sA[threadIdx.x].force[1]     = 0.0f;
                sA[threadIdx.x].force[2]     = 0.0f;
    
                sA[threadIdx.x].torque[0]    = 0.0f;
                sA[threadIdx.x].torque[1]    = 0.0f;
                sA[threadIdx.x].torque[2]    = 0.0f;
 */   
                if( bExclusionFlag )
                {
                    unsigned int xi   = x >> GRIDBITS;
                    unsigned int yi   = y >> GRIDBITS;
                    unsigned int cell = xi+yi*cSim.paddedNumberOfAtoms/GRID-yi*(yi+1)/2;
                    dScaleMask        = cAmoebaSim.pD_ScaleIndices[cAmoebaSim.pScaleIndicesIndex[cell]+tgx];
                    pScaleMask        = cAmoebaSim.pP_ScaleIndices[cAmoebaSim.pScaleIndicesIndex[cell]+tgx];
                    mScaleMask        = cAmoebaSim.pM_ScaleIndices[cAmoebaSim.pScaleIndicesIndex[cell]+tgx];
                } else {
                    scalingFactors[DScaleIndex] = scalingFactors[PScaleIndex] = scalingFactors[MScaleIndex] = 1.0f;
                }
       
                for (unsigned int j = 0; j < GRID; j++)
                {

                    // set scale factors and calculate force

                    if( bExclusionFlag ){
                        getMaskedDScaleFactor( tj, dScaleMask, scalingFactors + DScaleIndex );
                        getMaskedPScaleFactor( tj, pScaleMask, scalingFactors + PScaleIndex );
                        getMaskedMScaleFactor( tj, mScaleMask, scalingFactors + MScaleIndex );
                    }

                    // check if atoms out-of-bounds

                    if( (atomI < cSim.atoms) && ((y+tj) < cSim.atoms) )
                    {
                        calculatePmeDirectElectrostaticPairIxn_kernel( localParticleI, psA[tj], bExclusionFlag, scalingFactors, 1.0f, &totalEnergy);
                    } 
 
                    tj = (tj + 1) & (GRID - 1);

                } // end of j-loop

                localParticleI.force[0]    *=  conversionFactor;
                localParticleI.force[1]    *=  conversionFactor;
                localParticleI.force[2]    *=  conversionFactor;
    
                localParticleI.torque[0]   *= -conversionFactor;
                localParticleI.torque[1]   *= -conversionFactor;
                localParticleI.torque[2]   *= -conversionFactor;
    
                sA[threadIdx.x].force[0]   *=  conversionFactor;
                sA[threadIdx.x].force[1]   *=  conversionFactor;
                sA[threadIdx.x].force[2]   *=  conversionFactor;
    
                sA[threadIdx.x].torque[0]  *= -conversionFactor;
                sA[threadIdx.x].torque[1]  *= -conversionFactor;
                sA[threadIdx.x].torque[2]  *= -conversionFactor;
    
                // Write results
    
#ifdef USE_OUTPUT_BUFFER_PER_WARP
    
                unsigned int offset                 = (x + tgx + warp*cSim.paddedNumberOfAtoms);
                add3dArrayToFloat4( offset, localParticleI.force,  cSim.pForce4 );
                add3dArray(       3*offset, localParticleI.torque, outputTorque );
    
                offset                              = (y + tgx + warp*cSim.paddedNumberOfAtoms);
                add3dArrayToFloat4( offset, sA[threadIdx.x].force,  cSim.pForce4 );
                add3dArray(       3*offset, sA[threadIdx.x].torque, outputTorque );
    
#else
                unsigned int offset                 = (x + tgx + (y >> GRIDBITS) * cSim.paddedNumberOfAtoms);
                add3dArrayToFloat4( offset,  localParticleI.force,  cSim.pForce4 );
                load3dArray(       3*offset, localParticleI.torque, outputTorque );
    
                offset                              = (y + tgx + (x >> GRIDBITS) * cSim.paddedNumberOfAtoms);
    
                add3dArrayToFloat4( offset,  sA[threadIdx.x].force,  cSim.pForce4 );
                load3dArray(       3*offset, sA[threadIdx.x].torque, outputTorque );
    
#endif

            } // end of pInteractionFlag block
            lasty = y;
        }
        pos++;
    }
    cSim.pEnergy[blockIdx.x * blockDim.x + threadIdx.x] -= conversionFactor*totalEnergy;
}
