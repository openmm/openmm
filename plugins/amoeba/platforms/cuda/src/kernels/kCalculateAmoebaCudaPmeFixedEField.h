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
void METHOD_NAME(kCalculateAmoebaPmeDirectFixedE_Field, _kernel)(
                            unsigned int* workUnit,
                            float* outputEField,
                            float* outputEFieldPolar
#ifdef AMOEBA_DEBUG
                           , float4* debugArray, unsigned int targetAtom
#endif
){

#ifdef AMOEBA_DEBUG
    int pullIndexMax = 12;
    float4 pullBack[12];
    float dScaleVal;
    float pScaleVal;
#endif

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

            if (!bExclusionFlag)
            {

                // this branch is never exercised since it includes the
                // interaction between atomI and itself which is always excluded

                for (unsigned int j = 0; j < GRID; j++)
                {

                    float ijField[4][3];

                    // load coords, charge, ...

#ifdef AMOEBA_DEBUG
dScaleVal = 1.0f;
pScaleVal = 1.0f;
#endif
                    calculateFixedFieldRealSpacePairIxn_kernel( localParticle, psA[j], 1.0f, 1.0f, ijField
#ifdef AMOEBA_DEBUG
                                                , pullBack
#endif
                    );

                    unsigned int match      = (atomI == (y + j)) ? 1 : 0;
match = 1;

                    // add to field at atomI the field due atomJ's charge/dipole/quadrupole

                    fieldSum[0]            += match ? 0.0f : ijField[0][0];
                    fieldSum[1]            += match ? 0.0f : ijField[0][1];
                    fieldSum[2]            += match ? 0.0f : ijField[0][2];

                    fieldPolarSum[0]       += match ? 0.0f : ijField[2][0];
                    fieldPolarSum[1]       += match ? 0.0f : ijField[2][1];
                    fieldPolarSum[2]       += match ? 0.0f : ijField[2][2];

                }

            }
            else  // bExclusion
            {
                unsigned int xi       = x >> GRIDBITS;
                unsigned int cell     = xi + xi*cAmoebaSim.paddedNumberOfAtoms/GRID-xi*(xi+1)/2;
                int  dScaleMask       = cAmoebaSim.pD_ScaleIndices[cAmoebaSim.pScaleIndicesIndex[cell]+tgx];
                int2 pScaleMask       = cAmoebaSim.pP_ScaleIndices[cAmoebaSim.pScaleIndicesIndex[cell]+tgx];

                for (unsigned int j = 0; j < GRID; j++)
                {

                    // load coords, charge, ...

                    float ijField[4][3];
                    float dScaleValue;
                    float pScaleValue;

                    getMaskedDScaleFactor( j, dScaleMask, &dScaleValue );
                    getMaskedPScaleFactor( j, pScaleMask, &pScaleValue );

#ifdef AMOEBA_DEBUG
dScaleVal = dScaleValue;
pScaleVal = pScaleValue;
#endif
                    calculateFixedFieldRealSpacePairIxn_kernel( localParticle, psA[j], dScaleValue, pScaleValue, ijField
#ifdef AMOEBA_DEBUG
                                                , pullBack
#endif
                    );


                    // nan*0.0 = nan not 0.0, so explicitly exclude (atomI == atomJ) contribution
                    // by setting match flag

                    unsigned int match      = ( (atomI == (y + j)) || (atomI >= cAmoebaSim.numberOfAtoms) || ((y+j) >= cAmoebaSim.numberOfAtoms) ) ? 1 : 0;

                    // add to field at atomI the field due atomJ's charge/dipole/quadrupole

                    fieldSum[0]            += match ? 0.0f : ijField[0][0];
                    fieldSum[1]            += match ? 0.0f : ijField[0][1];
                    fieldSum[2]            += match ? 0.0f : ijField[0][2];

                    fieldPolarSum[0]       += match ? 0.0f : ijField[2][0];
                    fieldPolarSum[1]       += match ? 0.0f : ijField[2][1];
                    fieldPolarSum[2]       += match ? 0.0f : ijField[2][2];

#ifdef AMOEBA_DEBUG
if( atomI == targetAtom ){
        unsigned int index                 = atomI == targetAtom ? (y + j) : atomI;
        unsigned int indexI                = 0;
        unsigned int indexJ                = indexI ? 0 : 2;
        unsigned int indices[4]            = { indexI, indexJ, indexI+1, indexJ+1 };

        debugArray[index].x                = (float) atomI;
        debugArray[index].y                = (float) (y + j);
        debugArray[index].z                = dScaleValue;
        debugArray[index].w                = pScaleValue;

        index                             += cAmoebaSim.paddedNumberOfAtoms;
        unsigned int off                   = 3*(x + tgx + (x >> GRIDBITS) * cAmoebaSim.paddedNumberOfAtoms);
        debugArray[index].x                = (float) x;
        debugArray[index].y                = (float) tgx;
        debugArray[index].z                = -2;
        debugArray[index].w                = (float) off;

        float flag                         = 7.0f;
        for( int ii = 0; ii < 4; ii++ ){
            index                             += cAmoebaSim.paddedNumberOfAtoms;
            debugArray[index].x                = match ? 0.0f : ijField[indices[ii]][0];
            debugArray[index].y                = match ? 0.0f : ijField[indices[ii]][1];
            debugArray[index].z                = match ? 0.0f : ijField[indices[ii]][2];
            debugArray[index].w                = flag;
        }

        for( int pullIndex = 0; pullIndex < pullIndexMax; pullIndex++ ){
            index                             += cAmoebaSim.paddedNumberOfAtoms;
            debugArray[index].x                = pullBack[pullIndex].x;
            debugArray[index].y                = pullBack[pullIndex].y;
            debugArray[index].z                = pullBack[pullIndex].z;
            debugArray[index].w                = pullBack[pullIndex].w;
        }   


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
            }

            // Write results

#ifdef USE_OUTPUT_BUFFER_PER_WARP
            unsigned int offset                 = 3*(x + tgx + warp*cAmoebaSim.paddedNumberOfAtoms);
            load3dArrayBufferPerWarp( offset, fieldSum,       outputEField );
            load3dArrayBufferPerWarp( offset, fieldPolarSum,  outputEFieldPolar );
#else
            unsigned int offset                 = 3*(x + tgx + (x >> GRIDBITS) * cAmoebaSim.paddedNumberOfAtoms);
            load3dArray( offset, fieldSum,       outputEField );
            load3dArray( offset, fieldPolarSum,  outputEFieldPolar );
#endif

        }
        else        // 100% utilization
        {
            // Read fixed atom data into registers and GRF
            if (lasty != y)
            {

                // load coordinates, charge, ...

                loadFixedFieldShared( &(sA[threadIdx.x]), (y+tgx) );

            }

            // zero shared fields

            zeroFixedFieldParticleSharedField( &(sA[threadIdx.x]) );

            if (!bExclusionFlag)
            {
                for (unsigned int j = 0; j < GRID; j++)
                {

                    float ijField[4][3];
    
                    // load coords, charge, ...
    
#ifdef AMOEBA_DEBUG
dScaleVal = 1.0f;
pScaleVal = 1.0f;
#endif
                    calculateFixedFieldRealSpacePairIxn_kernel( localParticle, psA[tj], 1.0f, 1.0f, ijField
#ifdef AMOEBA_DEBUG
                                                 , pullBack
#endif
                    );
    
                    unsigned int outOfBounds     = ( (atomI >= cAmoebaSim.numberOfAtoms) || ((y+tj) >= cAmoebaSim.numberOfAtoms) ) ? 1 : 0;

                    // add to field at atomI the field due atomJ's charge/dipole/quadrupole

                    fieldSum[0]                 += outOfBounds ? 0.0f : ijField[0][0];
                    fieldSum[1]                 += outOfBounds ? 0.0f : ijField[0][1];
                    fieldSum[2]                 += outOfBounds ? 0.0f : ijField[0][2];
        
                    fieldPolarSum[0]            += outOfBounds ? 0.0f : ijField[2][0];
                    fieldPolarSum[1]            += outOfBounds ? 0.0f : ijField[2][1];
                    fieldPolarSum[2]            += outOfBounds ? 0.0f : ijField[2][2];

                    // add to field at atomJ the field due atomI's charge/dipole/quadrupole

                    psA[tj].eField[0]           += outOfBounds ? 0.0f : ijField[1][0];
                    psA[tj].eField[1]           += outOfBounds ? 0.0f : ijField[1][1];
                    psA[tj].eField[2]           += outOfBounds ? 0.0f : ijField[1][2];

                    psA[tj].eFieldP[0]          += outOfBounds ? 0.0f : ijField[3][0];
                    psA[tj].eFieldP[1]          += outOfBounds ? 0.0f : ijField[3][1];
                    psA[tj].eFieldP[2]          += outOfBounds ? 0.0f : ijField[3][2];


#ifdef AMOEBA_DEBUG
if( (atomI == targetAtom  || (y + tj) == targetAtom) ){
            unsigned int index                 = (atomI == targetAtom) ? (y + tj) : atomI;
            unsigned int indexI                = (atomI == targetAtom) ? 0 : 2;
            unsigned int indexJ                = (atomI == targetAtom) ? 2 : 0;

            debugArray[index].x                = (float) atomI;
            debugArray[index].y                = (float) (y + tj);
            debugArray[index].z                = dScaleVal;
            debugArray[index].w                = pScaleVal;

            unsigned int pullBackIndex         = 0;
            index                             += cAmoebaSim.paddedNumberOfAtoms;
            debugArray[index].x                = pullBack[pullBackIndex].x;
            debugArray[index].y                = pullBack[pullBackIndex].y;
            debugArray[index].z                = pullBack[pullBackIndex].z;
            debugArray[index].w                = pullBack[pullBackIndex].w;;

            pullBackIndex++;
            index                             += cAmoebaSim.paddedNumberOfAtoms;
            debugArray[index].x                = pullBack[pullBackIndex].x;
            debugArray[index].y                = pullBack[pullBackIndex].y;
            debugArray[index].z                = pullBack[pullBackIndex].z;
            debugArray[index].w                = pullBack[pullBackIndex].w;;


            float flag                         = 8.0f;
            index                             += cAmoebaSim.paddedNumberOfAtoms;
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

#if 0

            index                             += cAmoebaSim.paddedNumberOfAtoms;
            unsigned int mask                  = 1 << j;
            unsigned int pScaleIndex           = (scaleMask.x & mask) ? 1 : 0;
            pScaleIndex                       += (scaleMask.y & mask) ? 2 : 0;
            debugArray[index].x                = (float) pScaleIndex;

            debugArray[index].y                = scaleMask.x & mask ? 1.0f : -1.0f;
            debugArray[index].z                = scaleMask.y & mask ? 1.0f : -1.0f;
            debugArray[index].w                = pScaleValue + 10.0f;
#endif
}
#endif

                    tj                  = (tj + 1) & (GRID - 1);

                }
            }
            else  // bExclusion
            {
                // Read fixed atom data into registers and GRF

                unsigned int xi   = x >> GRIDBITS;
                unsigned int yi   = y >> GRIDBITS;
                unsigned int cell = xi+yi*cAmoebaSim.paddedNumberOfAtoms/GRID-yi*(yi+1)/2;
                int  dScaleMask   = cAmoebaSim.pD_ScaleIndices[cAmoebaSim.pScaleIndicesIndex[cell]+tgx];
                int2 pScaleMask   = cAmoebaSim.pP_ScaleIndices[cAmoebaSim.pScaleIndicesIndex[cell]+tgx];

                for (unsigned int j = 0; j < GRID; j++)
                {
                    // load coords, charge, ...

                    float ijField[4][3];

                    float dScaleValue;
                    float pScaleValue;
                    getMaskedDScaleFactor( tj, dScaleMask, &dScaleValue );
                    getMaskedPScaleFactor( tj, pScaleMask, &pScaleValue );

#ifdef AMOEBA_DEBUG
dScaleVal = dScaleValue;
pScaleVal = pScaleValue;
#endif
                    calculateFixedFieldRealSpacePairIxn_kernel( localParticle, psA[tj], dScaleValue, pScaleValue, ijField
#ifdef AMOEBA_DEBUG
                                                , pullBack
#endif
                    );

                    unsigned int outOfBounds     = ( (atomI >= cAmoebaSim.numberOfAtoms) || ((y+tj) >= cAmoebaSim.numberOfAtoms) ) ? 1 : 0;

                    // add to field at atomI the field due atomJ's charge/dipole/quadrupole

                    fieldSum[0]                 += outOfBounds ? 0.0f : ijField[0][0];
                    fieldSum[1]                 += outOfBounds ? 0.0f : ijField[0][1];
                    fieldSum[2]                 += outOfBounds ? 0.0f : ijField[0][2];

                    fieldPolarSum[0]            += outOfBounds ? 0.0f : ijField[2][0];
                    fieldPolarSum[1]            += outOfBounds ? 0.0f : ijField[2][1];
                    fieldPolarSum[2]            += outOfBounds ? 0.0f : ijField[2][2];

                    // add to field at atomJ the field due atomI's charge/dipole/quadrupole

                    psA[tj].eField[0]           += outOfBounds ? 0.0f : ijField[1][0];
                    psA[tj].eField[1]           += outOfBounds ? 0.0f : ijField[1][1];
                    psA[tj].eField[2]           += outOfBounds ? 0.0f : ijField[1][2];

                    psA[tj].eFieldP[0]          += outOfBounds ? 0.0f : ijField[3][0];
                    psA[tj].eFieldP[1]          += outOfBounds ? 0.0f : ijField[3][1];
                    psA[tj].eFieldP[2]          += outOfBounds ? 0.0f : ijField[3][2];


#ifdef AMOEBA_DEBUG
if( (atomI == targetAtom || (y + tj) == targetAtom) ){

            unsigned int index                 = (atomI == targetAtom) ? (y + tj) : atomI;
            unsigned int indexI                = (atomI == targetAtom) ? 0 : 2;
            unsigned int indexJ                = (atomI == targetAtom) ? 2 : 0;

            debugArray[index].x                = (float) atomI;
            debugArray[index].y                = (float) (y + tj);
            debugArray[index].z                = dScaleVal;
            debugArray[index].w                = pScaleVal;

            unsigned int pullBackIndex         = 0;
            index                             += cAmoebaSim.paddedNumberOfAtoms;
            debugArray[index].x                = pullBack[pullBackIndex].x;
            debugArray[index].y                = pullBack[pullBackIndex].y;
            debugArray[index].z                = pullBack[pullBackIndex].z;
            debugArray[index].w                = pullBack[pullBackIndex].w;;

            pullBackIndex++;
            index                             += cAmoebaSim.paddedNumberOfAtoms;
            debugArray[index].x                = pullBack[pullBackIndex].x;
            debugArray[index].y                = pullBack[pullBackIndex].y;
            debugArray[index].z                = pullBack[pullBackIndex].z;
            debugArray[index].w                = pullBack[pullBackIndex].w;;


            float flag                         = 9.0f;
            index                             += cAmoebaSim.paddedNumberOfAtoms;
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
                }
            }

            // Write results


#ifdef USE_OUTPUT_BUFFER_PER_WARP
            unsigned int offset                 = 3*(x + tgx + warp*cAmoebaSim.paddedNumberOfAtoms);
            load3dArrayBufferPerWarp( offset, fieldSum,       outputEField );
            load3dArrayBufferPerWarp( offset, fieldPolarSum,  outputEFieldPolar );

            offset                              = 3*(y + tgx + warp*cAmoebaSim.paddedNumberOfAtoms);
            load3dArrayBufferPerWarp( offset, sA[threadIdx.x].eField,  outputEField );
            load3dArrayBufferPerWarp( offset, sA[threadIdx.x].eFieldP, outputEFieldPolar );

#else
            unsigned int offset                 = 3*(x + tgx + (y >> GRIDBITS) * cAmoebaSim.paddedNumberOfAtoms);
            load3dArray( offset, fieldSum,       outputEField );
            load3dArray( offset, fieldPolarSum,  outputEFieldPolar );

            offset                              = 3*(y + tgx + (x >> GRIDBITS) * cAmoebaSim.paddedNumberOfAtoms);
            load3dArray( offset, sA[threadIdx.x].eField,  outputEField );
            load3dArray( offset, sA[threadIdx.x].eFieldP, outputEFieldPolar );
 
#endif
            lasty = y;
        }

        pos++;
    }
}
