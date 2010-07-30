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
__launch_bounds__(256, 1)
#elif (__CUDA_ARCH__ >= 130)
__launch_bounds__(128, 1)
#else
__launch_bounds__(64, 1)
#endif
void METHOD_NAME(kCalculateAmoebaMutualInducedAndGkFields, _kernel)(
                            unsigned int* workUnit,
                            float4* atomCoord,
                            float*  bornRadii,
                            float* inducedDipole,
                            float* inducedDipolePolar,
                            float* inducedDipoleS,
                            float* inducedDipolePolarS,
                            float* outputField,
                            float* outputFieldPolar,
                            float* outputFieldS,
                            float* outputFieldPolarS
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

    float4 jCoord;
    float  jBornRadius;
    float  jDipole[3];     
    float  jDipolePolar[3];     
    float  jDipoleS[3];     
    float  jDipolePolarS[3];     

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
        float4 iCoord                    = atomCoord[atomI];

        float fieldSum[3];
        float fieldPolarSum[3];
        float fieldSumS[3];
        float fieldPolarSumS[3];

        // fieldSum:      field at i due to j
        // fieldPolarSum: field at i due to j polar

        fieldSum[0]                      = 0.0f;
        fieldSum[1]                      = 0.0f;
        fieldSum[2]                      = 0.0f;

        fieldPolarSum[0]                 = 0.0f;
        fieldPolarSum[1]                 = 0.0f;
        fieldPolarSum[2]                 = 0.0f;

        fieldSumS[0]                     = 0.0f;
        fieldSumS[1]                     = 0.0f;
        fieldSumS[2]                     = 0.0f;

        fieldPolarSumS[0]                = 0.0f;
        fieldPolarSumS[1]                = 0.0f;
        fieldPolarSumS[2]                = 0.0f;

        if (x == y) 
        {

            // load shared data

            loadMutualInducedShared( &(sA[threadIdx.x]), atomI,
                                     atomCoord, inducedDipole, inducedDipolePolar, cAmoebaSim.pDampingFactorAndThole,
                                     bornRadii, inducedDipoleS, inducedDipolePolarS );

            for (unsigned int j = 0; j < GRID; j++)
            {

                float ijField[8][3];

                // load coords, charge, ...

                loadMutualInducedData( &(psA[j]), &jCoord, jDipole, jDipolePolar, &jBornRadius, jDipoleS, jDipolePolarS ); 

                calculateMutualInducedAndGkFieldsPairIxn_kernel( iCoord,                                          jCoord,
                                         cAmoebaSim.pDampingFactorAndThole[atomI].x,      psA[j].damp,
                                         cAmoebaSim.pDampingFactorAndThole[atomI].y,      psA[j].thole,
                                         &(inducedDipole[atomI*3]),                       &(inducedDipolePolar[atomI*3]),
                                         jDipole,                                         jDipolePolar,
                                         &(inducedDipoleS[atomI*3]),                      &(inducedDipolePolarS[atomI*3]),
                                         jDipoleS,                                        jDipolePolarS,
                                         cAmoebaSim.scalingDistanceCutoff,                ijField
#ifdef AMOEBA_DEBUG
,  debugArray
#endif
);

                unsigned int mask       =  ( (atomI == (y + j)) || (atomI >= cAmoebaSim.numberOfAtoms) || ((y+j) >= cAmoebaSim.numberOfAtoms) ) ? 0 : 1;

                // add to field at atomI the field due atomJ's dipole

                fieldSum[0]            += mask ? ijField[0][0] : 0.0f;
                fieldSum[1]            += mask ? ijField[0][1] : 0.0f;
                fieldSum[2]            += mask ? ijField[0][2] : 0.0f;

                fieldPolarSum[0]       += mask ? ijField[1][0] : 0.0f;
                fieldPolarSum[1]       += mask ? ijField[1][1] : 0.0f;
                fieldPolarSum[2]       += mask ? ijField[1][2] : 0.0f;

                fieldSumS[0]           += mask ? ijField[4][0] : 0.0f;
                fieldSumS[1]           += mask ? ijField[4][1] : 0.0f;
                fieldSumS[2]           += mask ? ijField[4][2] : 0.0f;

                fieldPolarSumS[0]      += mask ? ijField[5][0] : 0.0f;
                fieldPolarSumS[1]      += mask ? ijField[5][1] : 0.0f;
                fieldPolarSumS[2]      += mask ? ijField[5][2] : 0.0f;

//#ifdef AMOEBA_DEBUG
#if 0
unsigned int index                 = y + j;
if( atomI == targetAtom ){

        debugArray[index].x                = (float) atomI;
        debugArray[index].y                = (float) (y + j);
        debugArray[index].z                = bornRadii[atomI];
        debugArray[index].w                = jBornRadius;

        index                              = debugAccumulate( index, debugArray, ijField[0],   mask, 1.0f );
        index                              = debugAccumulate( index, debugArray, ijField[1],   mask, 2.0f );
        index                              = debugAccumulate( index, debugArray, ijField[4],   mask, 3.0f );
        index                              = debugAccumulate( index, debugArray, ijField[5],   mask, 4.0f );
}
#endif
                calculateMutualInducedAndGkFieldsGkPairIxn_kernel( iCoord,                                          jCoord, bornRadii[atomI]*jBornRadius,
                                                                   &(inducedDipoleS[atomI*3]),                      &(inducedDipolePolarS[atomI*3]),
                                                                   jDipoleS,                                        jDipolePolarS,
                                                                   ijField
#ifdef AMOEBA_DEBUG
                                                                   , debugArray
#endif
);

                // atomI == atomJ contribution included

                mask                    =  ( (atomI >= cAmoebaSim.numberOfAtoms) || ((y+j) >= cAmoebaSim.numberOfAtoms) ) ? 0 : 1;
                fieldSumS[0]           += mask ? ijField[0][0] : 0.0f;
                fieldSumS[1]           += mask ? ijField[0][1] : 0.0f;
                fieldSumS[2]           += mask ? ijField[0][2] : 0.0f;

                fieldPolarSumS[0]      += mask ? ijField[2][0] : 0.0f;
                fieldPolarSumS[1]      += mask ? ijField[2][1] : 0.0f;
                fieldPolarSumS[2]      += mask ? ijField[2][2] : 0.0f;

//#ifdef AMOEBA_DEBUG
#if 0
if( atomI == targetAtom ){

        index                              = debugAccumulate( index, debugArray, ijField[0],   mask, 5.0f );
        index                              = debugAccumulate( index, debugArray, ijField[2],   mask, 6.0f );

        index                              = debugAccumulate( index, debugArray, jDipoleS, 1, 7.0f );
        index                              = debugAccumulate( index, debugArray, jDipolePolarS, 1, 8.0f );

        index                             += cAmoebaSim.paddedNumberOfAtoms;
        debugArray[index].x                = bornRadii[atomI];
        debugArray[index].y                = jBornRadius;
        debugArray[index].w                = 9.0f;


}
#endif
            }

            // Write results

#ifdef USE_OUTPUT_BUFFER_PER_WARP

            unsigned int offset                 = 3*(x + tgx + warp*cAmoebaSim.paddedNumberOfAtoms);

            load3dArrayBufferPerWarp( offset, fieldSum,       outputField );
            load3dArrayBufferPerWarp( offset, fieldPolarSum,  outputFieldPolar );

            load3dArrayBufferPerWarp( offset, fieldSumS,      outputFieldS );
            load3dArrayBufferPerWarp( offset, fieldPolarSumS, outputFieldPolarS );

#else
            unsigned int offset                   = 3*(x + tgx + (x >> GRIDBITS) * cAmoebaSim.paddedNumberOfAtoms);

            load3dArray( offset, fieldSum,        outputField );
            load3dArray( offset, fieldPolarSum,   outputFieldPolar);

            load3dArray( offset, fieldSumS,       outputFieldS );
            load3dArray( offset, fieldPolarSumS,  outputFieldPolarS );
#endif

        }
        else        // 100% utilization
        {
            // Read fixed atom data into registers and GRF
            if (lasty != y)
            {
                // load coordinates, charge, ...

                loadMutualInducedShared( &(sA[threadIdx.x]), (y+tgx),
                                         atomCoord, inducedDipole,
                                         inducedDipolePolar, cAmoebaSim.pDampingFactorAndThole, bornRadii, inducedDipoleS, inducedDipolePolarS );
            }

           // zero shared fields

            zeroMutualInducedParticleSharedField(  &(sA[threadIdx.x]) );

            for (unsigned int j = 0; j < GRID; j++)
            {

                float ijField[8][3];

                // load coords, charge, ...

                loadMutualInducedData( &(psA[tj]), &jCoord, jDipole, jDipolePolar, &jBornRadius, jDipoleS, jDipolePolarS ); 

                calculateMutualInducedAndGkFieldsPairIxn_kernel( iCoord,                  jCoord,
                                         cAmoebaSim.pDampingFactorAndThole[atomI].x,      psA[tj].damp,
                                         cAmoebaSim.pDampingFactorAndThole[atomI].y,      psA[tj].thole,
                                         &(inducedDipole[atomI*3]),                       &(inducedDipolePolar[atomI*3]),
                                         jDipole,                                         jDipolePolar,
                                         &(inducedDipoleS[atomI*3]),                      &(inducedDipolePolarS[atomI*3]),
                                         jDipoleS,                                        jDipolePolarS,
                                         cAmoebaSim.scalingDistanceCutoff, ijField
#ifdef AMOEBA_DEBUG
,  debugArray
#endif
   );

                if( (atomI < cAmoebaSim.numberOfAtoms) && ((y+tj) < cAmoebaSim.numberOfAtoms) ){
           
                    // add to field at atomI the field due atomJ's dipole
    
                    fieldSum[0]              += ijField[0][0];
                    fieldSum[1]              += ijField[0][1];
                    fieldSum[2]              += ijField[0][2];
        
                    // add to polar field at atomI the field due atomJ's dipole
    
                    fieldPolarSum[0]         += ijField[1][0];
                    fieldPolarSum[1]         += ijField[1][1];
                    fieldPolarSum[2]         += ijField[1][2];
    
                    fieldSumS[0]             += ijField[4][0];
                    fieldSumS[1]             += ijField[4][1];
                    fieldSumS[2]             += ijField[4][2];
    
                    fieldPolarSumS[0]        += ijField[5][0];
                    fieldPolarSumS[1]        += ijField[5][1];
                    fieldPolarSumS[2]        += ijField[5][2];
    
                    // add to field at atomJ the field due atomI's dipole
    
                    psA[tj].field[0]         += ijField[2][0];
                    psA[tj].field[1]         += ijField[2][1];
                    psA[tj].field[2]         += ijField[2][2];
    
                    // add to polar field at atomJ the field due atomI's dipole
    
                    psA[tj].fieldPolar[0]    += ijField[3][0];
                    psA[tj].fieldPolar[1]    += ijField[3][1];
                    psA[tj].fieldPolar[2]    += ijField[3][2];
    
                    // add to field at atomJ the field due atomI's dipole
    
                    psA[tj].fieldS[0]        += ijField[6][0];
                    psA[tj].fieldS[1]        += ijField[6][1];
                    psA[tj].fieldS[2]        += ijField[6][2];
    
                    // add to polar field at atomJ the field due atomI's dipole
    
                    psA[tj].fieldPolarS[0]   += ijField[7][0];
                    psA[tj].fieldPolarS[1]   += ijField[7][1];
                    psA[tj].fieldPolarS[2]   += ijField[7][2];
    
                }

//#ifdef AMOEBA_DEBUG
#if 0
unsigned int index                 = (atomI == targetAtom) ? (y + tj) : atomI;
if( atomI == targetAtom  || (y + tj) == targetAtom ){
        unsigned int indexI                = (atomI == targetAtom) ? 0 : 2;
        unsigned int maskD                 = (atomI < cAmoebaSim.numberOfAtoms) && ((y+tj) < cAmoebaSim.numberOfAtoms);

        debugArray[index].x                = (float) atomI;
        debugArray[index].y                = (float) (y + tj);
        debugArray[index].z                = bornRadii[atomI];
        debugArray[index].w                = jBornRadius;

        index                              = debugAccumulate( index, debugArray, ijField[indexI],   maskD, -1.0f );
        index                              = debugAccumulate( index, debugArray, ijField[indexI+1], maskD, -2.0f );
        index                              = debugAccumulate( index, debugArray, ijField[indexI+4], maskD, -3.0f );
        index                              = debugAccumulate( index, debugArray, ijField[indexI+5], maskD, -4.0f );
}
#endif
                calculateMutualInducedAndGkFieldsGkPairIxn_kernel( iCoord,                                          jCoord, bornRadii[atomI]*jBornRadius,
                                                                   &(inducedDipoleS[atomI*3]),                      &(inducedDipolePolarS[atomI*3]),
                                                                   jDipoleS,                                        jDipolePolarS,
                                                                   ijField
#ifdef AMOEBA_DEBUG
,  debugArray
#endif
);


                if( (atomI < cAmoebaSim.numberOfAtoms) && ((y+tj) < cAmoebaSim.numberOfAtoms) ){

                    fieldSumS[0]           += ijField[0][0];
                    fieldSumS[1]           += ijField[0][1];
                    fieldSumS[2]           += ijField[0][2];
    
                    fieldPolarSumS[0]      += ijField[2][0];
                    fieldPolarSumS[1]      += ijField[2][1];
                    fieldPolarSumS[2]      += ijField[2][2];

                    // add to field at atomJ the field due atomI's dipole
    
                    psA[tj].fieldS[0]      += ijField[1][0];
                    psA[tj].fieldS[1]      += ijField[1][1];
                    psA[tj].fieldS[2]      += ijField[1][2];
    
                    // add to polar field at atomJ the field due atomI's dipole
    
                    psA[tj].fieldPolarS[0] += ijField[3][0];
                    psA[tj].fieldPolarS[1] += ijField[3][1];
                    psA[tj].fieldPolarS[2] += ijField[3][2];
                }
   
//#ifdef AMOEBA_DEBUG
#if 0
if( atomI == targetAtom  || (y + tj) == targetAtom ){
        unsigned int indexI                = (atomI == targetAtom) ? 0 : 1;
        unsigned int maskD                 = (atomI < cAmoebaSim.numberOfAtoms) && ((y+tj) < cAmoebaSim.numberOfAtoms);

        index                              = debugAccumulate( index, debugArray, ijField[indexI],   maskD, -5.0f );
        index                              = debugAccumulate( index, debugArray, ijField[indexI+2], maskD, -6.0f );
        index                              = debugAccumulate( index, debugArray, jDipoleS, 1, -7.0f );
        index                              = debugAccumulate( index, debugArray, jDipolePolarS, 1, -8.0f );

        index                             += cAmoebaSim.paddedNumberOfAtoms;
        debugArray[index].x                = bornRadii[atomI];
        debugArray[index].y                = jBornRadius;
        debugArray[index].w                = -9.0f;
}
#endif

                tj                  = (tj + 1) & (GRID - 1);

            }

            // Write results

#ifdef USE_OUTPUT_BUFFER_PER_WARP
            unsigned int offset                 = 3*(x + tgx + warp*cAmoebaSim.paddedNumberOfAtoms);
            load3dArrayBufferPerWarp( offset, fieldSum,       outputField );
            load3dArrayBufferPerWarp( offset, fieldPolarSum,  outputFieldPolar);
            load3dArrayBufferPerWarp( offset, fieldSumS,      outputFieldS );
            load3dArrayBufferPerWarp( offset, fieldPolarSumS, outputFieldPolarS );

            offset                              = 3*(y + tgx + warp*cAmoebaSim.paddedNumberOfAtoms);

            load3dArrayBufferPerWarp( offset, sA[threadIdx.x].field,       outputField );
            load3dArrayBufferPerWarp( offset, sA[threadIdx.x].fieldPolar,  outputFieldPolar);
            load3dArrayBufferPerWarp( offset, sA[threadIdx.x].fieldS,      outputFieldS );
            load3dArrayBufferPerWarp( offset, sA[threadIdx.x].fieldPolarS, outputFieldPolarS);

#else
            unsigned int offset                 = 3*(x + tgx + (y >> GRIDBITS) * cAmoebaSim.paddedNumberOfAtoms);
            load3dArray( offset, fieldSum,       outputField );
            load3dArray( offset, fieldPolarSum,  outputFieldPolar);
            load3dArray( offset, fieldSumS,      outputFieldS );
            load3dArray( offset, fieldPolarSumS, outputFieldPolarS);

            offset                              = 3*(y + tgx + (x >> GRIDBITS) * cAmoebaSim.paddedNumberOfAtoms);
            load3dArray( offset, sA[threadIdx.x].field,       outputField );
            load3dArray( offset, sA[threadIdx.x].fieldPolar,  outputFieldPolar);
            load3dArray( offset, sA[threadIdx.x].fieldS,      outputFieldS );
            load3dArray( offset, sA[threadIdx.x].fieldPolarS, outputFieldPolarS);

#endif
            lasty = y;
        }

        pos++;
    }
}
