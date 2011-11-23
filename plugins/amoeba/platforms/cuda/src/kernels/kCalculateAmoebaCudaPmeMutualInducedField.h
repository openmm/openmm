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
__launch_bounds__(GF1XX_NONBOND_THREADS_PER_BLOCK, 1)
#elif (__CUDA_ARCH__ >= 120)
__launch_bounds__(GT2XX_NONBOND_THREADS_PER_BLOCK, 1)
#else
__launch_bounds__(G8X_NONBOND_THREADS_PER_BLOCK, 1)
#endif
void METHOD_NAME(kCalculateAmoebaPmeMutualInducedField, _kernel)(
                            unsigned int* workUnit,
                            float* outputField, float* outputFieldPolar
){

    extern __shared__ MutualInducedParticle sA[];

    unsigned int totalWarps      = gridDim.x*blockDim.x/GRID;
    unsigned int warp            = (blockIdx.x*blockDim.x+threadIdx.x)/GRID;
    unsigned int numWorkUnits    = cSim.pInteractionCount[0];
    unsigned int pos             = warp*numWorkUnits/totalWarps;
    unsigned int end             = (warp+1)*numWorkUnits/totalWarps;
    unsigned int lasty           = 0xFFFFFFFF;
    const float uscale           = 1.0f;

    while (pos < end)
    {

        unsigned int x;
        unsigned int y;
        bool bExclusionFlag;

        // Extract cell coordinates

        decodeCell( workUnit[pos], &x, &y, &bExclusionFlag );

        unsigned int tgx                 = threadIdx.x & (GRID - 1);
        unsigned int tbx                 = threadIdx.x - tgx;
        unsigned int tj                  = tgx;

        MutualInducedParticle*  psA      = &sA[tbx];
        unsigned int atomI               = x + tgx;
        MutualInducedParticle localParticle;
        loadMutualInducedShared( &localParticle, atomI );

        float fieldSum[3];
        float fieldPolarSum[3];

        // 0: field at i due to j
        // 1: field at i due to j polar

        fieldSum[0]                      = 0.0f;
        fieldSum[1]                      = 0.0f;
        fieldSum[2]                      = 0.0f;

        fieldPolarSum[0]                 = 0.0f;
        fieldPolarSum[1]                 = 0.0f;
        fieldPolarSum[2]                 = 0.0f;

        if (x == y ){

            // load shared data

            loadMutualInducedShared( &(sA[threadIdx.x]), atomI );

            for (unsigned int j = 0; j < GRID; j++) {
                if(  ( (atomI != (y + j)) && (atomI < cSim.atoms) && ((y+j) < cSim.atoms) ) ){
                    float4 delta;
                    float prefactor2;
                    setupMutualInducedFieldPairIxn_kernel( localParticle, psA[j], uscale, &delta, &prefactor2 );
                    calculateMutualInducedFieldPairIxn_kernel(  psA[j].inducedDipole,      delta, prefactor2, fieldSum );
                    calculateMutualInducedFieldPairIxn_kernel(  psA[j].inducedDipolePolar, delta, prefactor2, fieldPolarSum );
                }

            }

            // Write results

#ifdef USE_OUTPUT_BUFFER_PER_WARP
            unsigned int offset            = 3*(x + tgx + warp*cSim.paddedNumberOfAtoms);
#else
            unsigned int offset            = 3*(x + tgx + (x >> GRIDBITS) * cSim.paddedNumberOfAtoms);
#endif
            load3dArray( offset, fieldSum,      outputField );
            load3dArray( offset, fieldPolarSum, outputFieldPolar);

        } else {

            if( lasty != y ){
                unsigned int atomJ        = y + tgx;
                loadMutualInducedShared( &(sA[threadIdx.x]), atomJ );
            }
    
            unsigned int flags = cSim.pInteractionFlag[pos];
            if( flags != 0 ){

#ifndef INCLUDE_MI_FIELD_BUFFERS
                flags = 0xFFFFFFFF;
#endif

               // zero shared fields
    
                zeroMutualInducedParticleSharedField(  &(sA[threadIdx.x]) );
    
                for (unsigned int j = 0; j < GRID; j++){
                    if ((flags&(1<<j)) != 0) {
                        unsigned int jIdx = (flags == 0xFFFFFFFF) ? tj : j;
                        float4 delta;
                        float prefactor2;
                        if( (atomI < cSim.atoms) && ((y+jIdx) < cSim.atoms) ){
                            setupMutualInducedFieldPairIxn_kernel( localParticle, psA[jIdx], uscale, &delta, &prefactor2 );
                            calculateMutualInducedFieldPairIxn_kernel(  psA[jIdx].inducedDipole,          delta, prefactor2, fieldSum );
                            calculateMutualInducedFieldPairIxn_kernel(  psA[jIdx].inducedDipolePolar,     delta, prefactor2, fieldPolarSum );
#ifndef INCLUDE_MI_FIELD_BUFFERS
                            calculateMutualInducedFieldPairIxn_kernel(  localParticle.inducedDipole,      delta, prefactor2, psA[jIdx].field );
                            calculateMutualInducedFieldPairIxn_kernel(  localParticle.inducedDipolePolar, delta, prefactor2, psA[jIdx].fieldPolar );
#else
                            if( flags == 0xFFFFFFFF ){
                                calculateMutualInducedFieldPairIxn_kernel(  localParticle.inducedDipole,      delta, prefactor2, psA[jIdx].field );
                                calculateMutualInducedFieldPairIxn_kernel(  localParticle.inducedDipolePolar, delta, prefactor2, psA[jIdx].fieldPolar );
                            } else {
                                calculateMutualInducedFieldPairIxnNoAdd_kernel(  localParticle.inducedDipole,      delta, prefactor2,  sA[threadIdx.x].tempBuffer );
                                calculateMutualInducedFieldPairIxnNoAdd_kernel(  localParticle.inducedDipolePolar, delta, prefactor2,  sA[threadIdx.x].tempBufferP );
    
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
                                    psA[jIdx].field[0]         += sA[threadIdx.x].tempBuffer[0]  + sA[threadIdx.x+16].tempBuffer[0];
                                    psA[jIdx].field[1]         += sA[threadIdx.x].tempBuffer[1]  + sA[threadIdx.x+16].tempBuffer[1];
                                    psA[jIdx].field[2]         += sA[threadIdx.x].tempBuffer[2]  + sA[threadIdx.x+16].tempBuffer[2];
    
                                    psA[jIdx].fieldPolar[0]    += sA[threadIdx.x].tempBufferP[0] + sA[threadIdx.x+16].tempBufferP[0];
                                    psA[jIdx].fieldPolar[1]    += sA[threadIdx.x].tempBufferP[1] + sA[threadIdx.x+16].tempBufferP[1];
                                    psA[jIdx].fieldPolar[2]    += sA[threadIdx.x].tempBufferP[2] + sA[threadIdx.x+16].tempBufferP[2];
                                }
    
                            }
#endif
                        }
                    }
    
                    tj = (tj + 1) & (GRID - 1);
    
                } // end of j-loop
    
                // Write results
    
#ifdef USE_OUTPUT_BUFFER_PER_WARP
                unsigned int offset     = 3*(x + tgx + warp*cSim.paddedNumberOfAtoms);
                load3dArrayBufferPerWarp( offset, fieldSum,      outputField );
                load3dArrayBufferPerWarp( offset, fieldPolarSum, outputFieldPolar);
    
                offset                  = 3*(y + tgx + warp*cSim.paddedNumberOfAtoms);
    
                load3dArrayBufferPerWarp( offset, sA[threadIdx.x].field,      outputField );
                load3dArrayBufferPerWarp( offset, sA[threadIdx.x].fieldPolar, outputFieldPolar);
    
#else
                unsigned int offset     = 3*(x + tgx + (y >> GRIDBITS) * cSim.paddedNumberOfAtoms);
                load3dArray( offset, fieldSum,      outputField );
                load3dArray( offset, fieldPolarSum, outputFieldPolar);
    
                offset                  = 3*(y + tgx + (x >> GRIDBITS) * cSim.paddedNumberOfAtoms);
                load3dArray( offset, sA[threadIdx.x].field,      outputField );
                load3dArray( offset, sA[threadIdx.x].fieldPolar, outputFieldPolar);
    
#endif
                lasty = y;
    
            } // end of pInteractionFlag block

        } // end of x == y block
        pos++;
    }
}
