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
#elif (__CUDA_ARCH__ >= 130)
__launch_bounds__(GT2XX_NONBOND_THREADS_PER_BLOCK, 1)
#else
__launch_bounds__(G8X_NONBOND_THREADS_PER_BLOCK, 1)
#endif
void METHOD_NAME(kCalculateAmoebaPmeMutualInducedField, _kernel)(
                            unsigned int* workUnit,
                            float* outputField, float* outputFieldPolar
#ifdef AMOEBA_DEBUG
                           , float4* debugArray, unsigned int targetAtom
#endif
){

    extern __shared__ MutualInducedParticle sA[];

    unsigned int totalWarps      = gridDim.x*blockDim.x/GRID;
    unsigned int warp            = (blockIdx.x*blockDim.x+threadIdx.x)/GRID;
    unsigned int numWorkUnits    = cSim.pInteractionCount[0];
    unsigned int pos             = warp*numWorkUnits/totalWarps;
    unsigned int end             = (warp+1)*numWorkUnits/totalWarps;
    unsigned int lasty           = 0xFFFFFFFF;
    const float uscale           = 1.0f;

#ifdef AMOEBA_DEBUG
    float4 pullBack[4];
#endif

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

        if (x == y) // Handle diagonals uniquely at 50% efficiency
        {

            // load shared data

            loadMutualInducedShared( &(sA[threadIdx.x]), atomI );

            for (unsigned int j = 0; j < GRID; j++)
            {

                float ijField[4][3];

                // load coords, charge, ...

                calculatePmeDirectMutualInducedFieldPairIxn_kernel( localParticle, psA[j], uscale, ijField
#ifdef AMOEBA_DEBUG
, pullBack 
#endif
);

                unsigned int mask       =  ( (atomI == (y + j)) || (atomI >= cAmoebaSim.numberOfAtoms) || ((y+j) >= cAmoebaSim.numberOfAtoms) ) ? 0 : 1;

                // add to field at atomI the field due atomJ's dipole

                fieldSum[0]            += mask ? ijField[0][0] : 0.0f;
                fieldSum[1]            += mask ? ijField[0][1] : 0.0f;
                fieldSum[2]            += mask ? ijField[0][2] : 0.0f;

                fieldPolarSum[0]       += mask ? ijField[2][0] : 0.0f;
                fieldPolarSum[1]       += mask ? ijField[2][1] : 0.0f;
                fieldPolarSum[2]       += mask ? ijField[2][2] : 0.0f;

#ifdef AMOEBA_DEBUG
if( atomI == targetAtom || (y+j) == targetAtom ){
            unsigned int index                 = atomI == targetAtom ? (y+j) : atomI;
            unsigned int pullBackIndex         = 0;
            unsigned int indexI                = 0;
            unsigned int indexJ                = indexI ? 0 : 2;

            debugArray[index].x                = (float) atomI;
            debugArray[index].y                = (float) (y + j);
            debugArray[index].z                = cSim.nonbondedCutoffSqr;
            debugArray[index].w                = 6.0f;


            index                             += cAmoebaSim.paddedNumberOfAtoms;
            debugArray[index].x                = pullBack[pullBackIndex].x;
            debugArray[index].y                = pullBack[pullBackIndex].y;
            debugArray[index].z                = pullBack[pullBackIndex].z;
            debugArray[index].w                = pullBack[pullBackIndex].w;

            pullBackIndex++;
            index                             += cAmoebaSim.paddedNumberOfAtoms;
            debugArray[index].x                = pullBack[pullBackIndex].x;
            debugArray[index].y                = pullBack[pullBackIndex].y;
            debugArray[index].z                = pullBack[pullBackIndex].z;
            debugArray[index].w                = pullBack[pullBackIndex].w;

            index                             += cAmoebaSim.paddedNumberOfAtoms;
            float flag                         = 6.0f;
            debugArray[index].x                = ijField[indexI][0];
            debugArray[index].y                = ijField[indexI][1];
            debugArray[index].z                = ijField[indexI][2];
            debugArray[index].w                = flag;

            index                             += cAmoebaSim.paddedNumberOfAtoms;
            debugArray[index].x                = ijField[indexJ][0];
            debugArray[index].y                = ijField[indexJ][1];
            debugArray[index].z                = ijField[indexJ][2];
            debugArray[index].w                = flag;

            index                             += cAmoebaSim.paddedNumberOfAtoms;
            debugArray[index].x                = ijField[indexI+1][0];
            debugArray[index].y                = ijField[indexI+1][1];
            debugArray[index].z                = ijField[indexI+1][2];
            debugArray[index].w                = flag;

            index                             += cAmoebaSim.paddedNumberOfAtoms;
            debugArray[index].x                = ijField[indexJ+1][0];
            debugArray[index].y                = ijField[indexJ+1][1];
            debugArray[index].z                = ijField[indexJ+1][2];
            debugArray[index].w                = flag;

/*
            index                             += cAmoebaSim.paddedNumberOfAtoms;

            index                             += cAmoebaSim.paddedNumberOfAtoms;
            debugArray[index].x                = match ? 0.0f : ijField[indexI][0];
            debugArray[index].y                = match ? 0.0f : ijField[indexI][1];
            debugArray[index].z                = match ? 0.0f : ijField[indexI][2];
            index                             += cAmoebaSim.paddedNumberOfAtoms;
            unsigned int mask                  = 1 << j;
            unsigned int pScaleIndex           = (scaleMask.x & mask) ? 1 : 0;
            pScaleIndex                       += (scaleMask.y & mask) ? 2 : 0;
            debugArray[index].x                = (float) pScaleIndex;

            debugArray[index].y                = scaleMask.x & mask ? 1.0f : -1.0f;
            debugArray[index].z                = scaleMask.y & mask ? 1.0f : -1.0f;
            debugArray[index].w                = + 10.0f;
*/

}
#endif
            }

            // Write results

#ifdef USE_OUTPUT_BUFFER_PER_WARP
            unsigned int offset            = 3*(x + tgx + warp*cAmoebaSim.paddedNumberOfAtoms);
            load3dArrayBufferPerWarp( offset, fieldSum,      outputField );
            load3dArrayBufferPerWarp( offset, fieldPolarSum, outputFieldPolar);

#else
            unsigned int offset            = 3*(x + tgx + (x >> GRIDBITS) * cAmoebaSim.paddedNumberOfAtoms);
            load3dArray( offset, fieldSum,      outputField );
            load3dArray( offset, fieldPolarSum, outputFieldPolar);

#endif

        } else {

            if (lasty != y)
            {
                unsigned int atomJ        = y + tgx;

                // load coordinates, charge, ...

                loadMutualInducedShared( &(sA[threadIdx.x]), atomJ );
            }
    
            unsigned int flags = cSim.pInteractionFlag[pos];
            if (flags == 0) {
                // No interactions in this block.
            } else {

               // zero shared fields
    
                zeroMutualInducedParticleSharedField(  &(sA[threadIdx.x]) );
    
                for (unsigned int j = 0; j < GRID; j++)
                {
    
                    unsigned int jIdx = (flags == 0xFFFFFFFF) ? tj : j;
                    float ijField[4][3];
    
                    // load coords, charge, ...
    
                    calculatePmeDirectMutualInducedFieldPairIxn_kernel( localParticle, psA[jIdx], uscale, ijField
#ifdef AMOEBA_DEBUG
    , pullBack 
#endif
       );
    
                    unsigned int mask   =  ( (atomI >= cAmoebaSim.numberOfAtoms) || ((y+jIdx) >= cAmoebaSim.numberOfAtoms) ) ? 0 : 1;
               
                    // add to field at atomI the field due atomJ's dipole
    
                    fieldSum[0]              += mask ? ijField[0][0] : 0.0f;
                    fieldSum[1]              += mask ? ijField[0][1] : 0.0f;
                    fieldSum[2]              += mask ? ijField[0][2] : 0.0f;
        
                    // add to polar field at atomI the field due atomJ's dipole
    
                    fieldPolarSum[0]         += mask ? ijField[2][0] : 0.0f;
                    fieldPolarSum[1]         += mask ? ijField[2][1] : 0.0f;
                    fieldPolarSum[2]         += mask ? ijField[2][2] : 0.0f;
    
                    // add to field at atomJ the field due atomI's dipole
    
                    if( flags == 0xFFFFFFFF ){

                        psA[jIdx].field[0]             += mask ? ijField[1][0] : 0.0f;
                        psA[jIdx].field[1]             += mask ? ijField[1][1] : 0.0f;
                        psA[jIdx].field[2]             += mask ? ijField[1][2] : 0.0f;
        
                        // add to polar field at atomJ the field due atomI's dipole
        
                        psA[jIdx].fieldPolar[0]        += mask ? ijField[3][0] : 0.0f;
                        psA[jIdx].fieldPolar[1]        += mask ? ijField[3][1] : 0.0f;
                        psA[jIdx].fieldPolar[2]        += mask ? ijField[3][2] : 0.0f;

                    } else {

                        sA[threadIdx.x].tempBuffer[0]  = mask ? 0.0f : ijField[1][0];
                        sA[threadIdx.x].tempBuffer[1]  = mask ? 0.0f : ijField[1][1];
                        sA[threadIdx.x].tempBuffer[2]  = mask ? 0.0f : ijField[1][2];

                        sA[threadIdx.x].tempBufferP[0] = mask ? 0.0f : ijField[3][0];
                        sA[threadIdx.x].tempBufferP[1] = mask ? 0.0f : ijField[3][1];
                        sA[threadIdx.x].tempBufferP[2] = mask ? 0.0f : ijField[3][2];

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
    
#ifdef AMOEBA_DEBUG
if( atomI == targetAtom || (y+jIdx) == targetAtom ){
            unsigned int index                 = atomI == targetAtom ? (y+jIdx) : atomI;
            unsigned int pullBackIndex         = 0;
            unsigned int indexI                = 0;
            unsigned int indexJ                = indexI ? 0 : 2;

            debugArray[index].x                = (float) atomI;
            debugArray[index].y                = (float) (y + jIdx);
            debugArray[index].z                = cSim.nonbondedCutoffSqr;
            debugArray[index].w                = 7.0f;


            index                             += cAmoebaSim.paddedNumberOfAtoms;
            debugArray[index].x                = pullBack[pullBackIndex].x;
            debugArray[index].y                = pullBack[pullBackIndex].y;
            debugArray[index].z                = pullBack[pullBackIndex].z;
            debugArray[index].w                = pullBack[pullBackIndex].w;

            pullBackIndex++;
            index                             += cAmoebaSim.paddedNumberOfAtoms;
            debugArray[index].x                = pullBack[pullBackIndex].x;
            debugArray[index].y                = pullBack[pullBackIndex].y;
            debugArray[index].z                = pullBack[pullBackIndex].z;
            debugArray[index].w                = pullBack[pullBackIndex].w;

            index                             += cAmoebaSim.paddedNumberOfAtoms;
            float flag                         = 7.0f;
            debugArray[index].x                = ijField[indexI][0];
            debugArray[index].y                = ijField[indexI][1];
            debugArray[index].z                = ijField[indexI][2];
            debugArray[index].w                = flag;

            index                             += cAmoebaSim.paddedNumberOfAtoms;
            debugArray[index].x                = ijField[indexJ][0];
            debugArray[index].y                = ijField[indexJ][1];
            debugArray[index].z                = ijField[indexJ][2];
            debugArray[index].w                = flag;

            index                             += cAmoebaSim.paddedNumberOfAtoms;
            debugArray[index].x                = ijField[indexI+1][0];
            debugArray[index].y                = ijField[indexI+1][1];
            debugArray[index].z                = ijField[indexI+1][2];
            debugArray[index].w                = flag;

            index                             += cAmoebaSim.paddedNumberOfAtoms;
            debugArray[index].x                = ijField[indexJ+1][0];
            debugArray[index].y                = ijField[indexJ+1][1];
            debugArray[index].z                = ijField[indexJ+1][2];
            debugArray[index].w                = flag;
}
#endif
    
                    tj                  = (tj + 1) & (GRID - 1);
    
                } // end of j-loop
    
                // Write results
    
#ifdef USE_OUTPUT_BUFFER_PER_WARP
                unsigned int offset     = 3*(x + tgx + warp*cAmoebaSim.paddedNumberOfAtoms);
                load3dArrayBufferPerWarp( offset, fieldSum,      outputField );
                load3dArrayBufferPerWarp( offset, fieldPolarSum, outputFieldPolar);
    
                offset                  = 3*(y + tgx + warp*cAmoebaSim.paddedNumberOfAtoms);
    
                load3dArrayBufferPerWarp( offset, sA[threadIdx.x].field,      outputField );
                load3dArrayBufferPerWarp( offset, sA[threadIdx.x].fieldPolar, outputFieldPolar);
    
#else
                unsigned int offset     = 3*(x + tgx + (y >> GRIDBITS) * cAmoebaSim.paddedNumberOfAtoms);
                load3dArray( offset, fieldSum,      outputField );
                load3dArray( offset, fieldPolarSum, outputFieldPolar);
    
                offset                  = 3*(y + tgx + (x >> GRIDBITS) * cAmoebaSim.paddedNumberOfAtoms);
                load3dArray( offset, sA[threadIdx.x].field,      outputField );
                load3dArray( offset, sA[threadIdx.x].fieldPolar, outputFieldPolar);
    
#endif
                lasty = y;
    
            } // end of pInteractionFlag block

        } // end of x == y block
        pos++;
    }
}
