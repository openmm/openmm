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
__launch_bounds__(192, 1)
#else
__launch_bounds__(64, 1)
#endif
void METHOD_NAME(kCalculateAmoebaPmeDirectFixedE_Field, _kernel)(
                            unsigned int* workUnit,
                            float* outputEField,
                            float* outputEFieldPolar){ 

    extern __shared__ FixedFieldParticle sA[];

    unsigned int totalWarps      = gridDim.x*blockDim.x/GRID;
    unsigned int warp            = (blockIdx.x*blockDim.x+threadIdx.x)/GRID;
    unsigned int numWorkUnits    = cSim.pInteractionCount[0];
    unsigned int pos             = warp*numWorkUnits/totalWarps;
    unsigned int end             = (warp+1)*numWorkUnits/totalWarps;
    unsigned int lasty           = 0xFFFFFFFF;

    while (pos < end)
    {

        unsigned int x;
        unsigned int y;
        bool bExclusionFlag;
        float dScaleValue;
        float pScaleValue;
        int  dScaleMask;
        int2 pScaleMask;

        // extract cell coordinates

        decodeCell( workUnit[pos], &x, &y, &bExclusionFlag );

        unsigned int tgx           = threadIdx.x & (GRID - 1);
        unsigned int tbx           = threadIdx.x - tgx;
        unsigned int tj            = tgx;

        FixedFieldParticle* psA    = &sA[tbx];
        unsigned int atomI         = x + tgx;
        FixedFieldParticle localParticle;
        loadFixedFieldShared( &localParticle, atomI );

        float fieldSum[3];
        float fieldPolarSum[3];

        fieldSum[0]                = 0.0f;
        fieldSum[1]                = 0.0f;
        fieldSum[2]                = 0.0f;

        fieldPolarSum[0]           = 0.0f;
        fieldPolarSum[1]           = 0.0f;
        fieldPolarSum[2]           = 0.0f;

        if (x == y)
        {

            // load coordinates, charge, ...

            loadFixedFieldShared( &(sA[threadIdx.x]), atomI );

            if( bExclusionFlag ){
                unsigned int xi       = x >> GRIDBITS;
                unsigned int cell     = xi + xi*cSim.paddedNumberOfAtoms/GRID-xi*(xi+1)/2;
                dScaleMask            = cAmoebaSim.pD_ScaleIndices[cAmoebaSim.pScaleIndicesIndex[cell]+tgx];
                pScaleMask            = cAmoebaSim.pP_ScaleIndices[cAmoebaSim.pScaleIndicesIndex[cell]+tgx];
            } else {
                dScaleValue = pScaleValue = 1.0f;

            }

            for (unsigned int j = 0; j < GRID; j++)
            {

                if( bExclusionFlag ){
                    getMaskedDScaleFactor( j, dScaleMask, &dScaleValue );
                    getMaskedPScaleFactor( j, pScaleMask, &pScaleValue );
                }

                float4 ijField[3];
                calculateFixedFieldRealSpacePairIxn_kernel( localParticle, psA[j], dScaleValue, pScaleValue, ijField);

                // nan*0.0 = nan not 0.0, so explicitly exclude (atomI == atomJ) contribution
                // by setting match flag

                unsigned int match      = ( (atomI == (y + j)) || (atomI >= cSim.atoms) || ((y+j) >= cSim.atoms) ) ? 1 : 0;

                // add to field at atomI the field due atomJ's charge/dipole/quadrupole

                fieldSum[0]            += match ? 0.0f : ijField[0].x;
                fieldSum[1]            += match ? 0.0f : ijField[1].x;
                fieldSum[2]            += match ? 0.0f : ijField[2].x;

                fieldPolarSum[0]       += match ? 0.0f : ijField[0].z;
                fieldPolarSum[1]       += match ? 0.0f : ijField[1].z;
                fieldPolarSum[2]       += match ? 0.0f : ijField[2].z;

            }

            // Write results

#ifdef USE_OUTPUT_BUFFER_PER_WARP
            unsigned int offset                 = 3*(x + tgx + warp*cSim.paddedNumberOfAtoms);
            load3dArrayBufferPerWarp( offset, fieldSum,       outputEField );
            load3dArrayBufferPerWarp( offset, fieldPolarSum,  outputEFieldPolar );
#else
            unsigned int offset                 = 3*(x + tgx + (x >> GRIDBITS) * cSim.paddedNumberOfAtoms);
            load3dArray( offset, fieldSum,       outputEField );
            load3dArray( offset, fieldPolarSum,  outputEFieldPolar );
#endif

        } else {

            if (lasty != y ) {
    
                // load coordinates, charge, ...
    
                loadFixedFieldShared( &(sA[threadIdx.x]), (y+tgx) );
    
            }

            unsigned int flags = cSim.pInteractionFlag[pos];
            if (flags == 0) {
                // No interactions in this block.
            } else {

                // zero shared fields

                zeroFixedFieldParticleSharedField( &(sA[threadIdx.x]) );

                if( bExclusionFlag ) {
                    unsigned int xi   = x >> GRIDBITS;
                    unsigned int yi   = y >> GRIDBITS;
                    unsigned int cell = xi+yi*cSim.paddedNumberOfAtoms/GRID-yi*(yi+1)/2;
                    dScaleMask        = cAmoebaSim.pD_ScaleIndices[cAmoebaSim.pScaleIndicesIndex[cell]+tgx];
                    pScaleMask        = cAmoebaSim.pP_ScaleIndices[cAmoebaSim.pScaleIndicesIndex[cell]+tgx];
                } else {
                    dScaleValue = pScaleValue  = 1.0f;
                }

                for (unsigned int j = 0; j < GRID; j++){

                    if ((flags&(1<<j)) != 0) {
                        unsigned int jIdx = (flags == 0xFFFFFFFF) ? tj : j;
                        if( bExclusionFlag ){
                            getMaskedDScaleFactor( jIdx, dScaleMask, &dScaleValue );
                            getMaskedPScaleFactor( jIdx, pScaleMask, &pScaleValue );
                        }

                        float4 ijField[3];
                        calculateFixedFieldRealSpacePairIxn_kernel( localParticle, psA[jIdx], dScaleValue, pScaleValue, ijField);

                        unsigned int outOfBounds     = ( (atomI >= cSim.atoms) || ((y+jIdx) >= cSim.atoms) ) ? 1 : 0;

                        // add to field at atomI the field due atomJ's charge/dipole/quadrupole

                        fieldSum[0]                 += outOfBounds ? 0.0f : ijField[0].x;
                        fieldSum[1]                 += outOfBounds ? 0.0f : ijField[1].x;
                        fieldSum[2]                 += outOfBounds ? 0.0f : ijField[2].x;

                        fieldPolarSum[0]            += outOfBounds ? 0.0f : ijField[0].z;
                        fieldPolarSum[1]            += outOfBounds ? 0.0f : ijField[1].z;
                        fieldPolarSum[2]            += outOfBounds ? 0.0f : ijField[2].z;

                        if( flags == 0xFFFFFFFF ){

                            // add to field at atomJ the field due atomI's charge/dipole/quadrupole

                            psA[jIdx].eField[0]        += outOfBounds ? 0.0f : ijField[0].y;
                            psA[jIdx].eField[1]        += outOfBounds ? 0.0f : ijField[1].y;
                            psA[jIdx].eField[2]        += outOfBounds ? 0.0f : ijField[2].y;

                            psA[jIdx].eFieldP[0]       += outOfBounds ? 0.0f : ijField[0].w;
                            psA[jIdx].eFieldP[1]       += outOfBounds ? 0.0f : ijField[1].w;
                            psA[jIdx].eFieldP[2]       += outOfBounds ? 0.0f : ijField[2].w;

                        } else {

                            sA[threadIdx.x].tempBuffer[0]  = outOfBounds ? 0.0f : ijField[0].y;
                            sA[threadIdx.x].tempBuffer[1]  = outOfBounds ? 0.0f : ijField[1].y;
                            sA[threadIdx.x].tempBuffer[2]  = outOfBounds ? 0.0f : ijField[2].y;

                            sA[threadIdx.x].tempBufferP[0] = outOfBounds ? 0.0f : ijField[0].w;
                            sA[threadIdx.x].tempBufferP[1] = outOfBounds ? 0.0f : ijField[1].w;
                            sA[threadIdx.x].tempBufferP[2] = outOfBounds ? 0.0f : ijField[2].w;

                            if( tgx % 2 == 0 ){
                                sumTempBuffer( sA[threadIdx.x], sA[threadIdx.x+1] );
                            }
                            if( tgx % 4 == 0 ){
                                sumTempBuffer( sA[threadIdx.x], sA[threadIdx.x+2] );
                            }
                            if( tgx % 8 == 0 ){
                                sumTempBuffer( sA[threadIdx.x], sA[threadIdx.x+4] );
                            }
                            if( tgx % 16 == 0 ){
                                sumTempBuffer( sA[threadIdx.x], sA[threadIdx.x+8] );
                            }

                            if (tgx == 0)
                            {
                                psA[jIdx].eField[0]  += sA[threadIdx.x].tempBuffer[0]  + sA[threadIdx.x+16].tempBuffer[0];
                                psA[jIdx].eField[1]  += sA[threadIdx.x].tempBuffer[1]  + sA[threadIdx.x+16].tempBuffer[1];
                                psA[jIdx].eField[2]  += sA[threadIdx.x].tempBuffer[2]  + sA[threadIdx.x+16].tempBuffer[2];

                                psA[jIdx].eFieldP[0] += sA[threadIdx.x].tempBufferP[0] + sA[threadIdx.x+16].tempBufferP[0];
                                psA[jIdx].eFieldP[1] += sA[threadIdx.x].tempBufferP[1] + sA[threadIdx.x+16].tempBufferP[1];
                                psA[jIdx].eFieldP[2] += sA[threadIdx.x].tempBufferP[2] + sA[threadIdx.x+16].tempBufferP[2];
                            }
                        }

                    }
                    tj = (tj + 1) & (GRID - 1);

                } // j-loop block
    
                // Write results
    
#ifdef USE_OUTPUT_BUFFER_PER_WARP
                unsigned int offset                 = 3*(x + tgx + warp*cSim.paddedNumberOfAtoms);
                load3dArrayBufferPerWarp( offset, fieldSum,       outputEField );
                load3dArrayBufferPerWarp( offset, fieldPolarSum,  outputEFieldPolar );
    
                offset                              = 3*(y + tgx + warp*cSim.paddedNumberOfAtoms);
                load3dArrayBufferPerWarp( offset, sA[threadIdx.x].eField,  outputEField );
                load3dArrayBufferPerWarp( offset, sA[threadIdx.x].eFieldP, outputEFieldPolar );
    
#else
                unsigned int offset                 = 3*(x + tgx + (y >> GRIDBITS) * cSim.paddedNumberOfAtoms);
                load3dArray( offset, fieldSum,       outputEField );
                load3dArray( offset, fieldPolarSum,  outputEFieldPolar );
    
                offset                              = 3*(y + tgx + (x >> GRIDBITS) * cSim.paddedNumberOfAtoms);
                load3dArray( offset, sA[threadIdx.x].eField,  outputEField );
                load3dArray( offset, sA[threadIdx.x].eFieldP, outputEFieldPolar );
     
#endif
            } // end of pInteractionFlag block 
            lasty = y;
        } // x == y block

        pos++;
    }
}
