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
void METHOD_NAME(kCalculateAmoebaFixedEAndGkField, _kernel)(
                            unsigned int* workUnit,
                            float4* atomCoord,
                            float* labFrameDipole,
                            float* labFrameQuadrupole,
                            float* bornRadii,
                            float* outputEField,
                            float* outputEFieldPolar,
                            float* outputGkField
#ifdef AMOEBA_DEBUG
                           , float4* debugArray, unsigned int targetAtom
#endif
){

#ifdef AMOEBA_DEBUG
    float4 pullBack[12];
#endif

    extern __shared__ FixedFieldParticle sA[];

    unsigned int totalWarps      = gridDim.x*blockDim.x/GRID;
    unsigned int warp            = (blockIdx.x*blockDim.x+threadIdx.x)/GRID;
    unsigned int numWorkUnits    = cSim.pInteractionCount[0];
    unsigned int pos             = warp*numWorkUnits/totalWarps;
    unsigned int end             = (warp+1)*numWorkUnits/totalWarps;
    unsigned int lasty           = 0xFFFFFFFF;

    float4 jCoord;
    float  jBornRadius;
    float  jDipole[3];    
    float  jQuadrupole[9];    
 
    while (pos < end)
    {

        unsigned int x;
        unsigned int y;
        bool bExclusionFlag;

        // Extract cell coordinates

        decodeCell( workUnit[pos], &x, &y, &bExclusionFlag );

        unsigned int tgx           = threadIdx.x & (GRID - 1);
        unsigned int tbx           = threadIdx.x - tgx;
        unsigned int tj            = tgx;

        FixedFieldParticle* psA    = &sA[tbx];
        unsigned int atomI         = x + tgx;
        float4 iCoord              = atomCoord[atomI];

        float eFieldSum[3];
        float eFieldPolarSum[3];
        float gkFieldSum[3];

        eFieldSum[0]               = 0.0f;
        eFieldSum[1]               = 0.0f;
        eFieldSum[2]               = 0.0f;

        eFieldPolarSum[0]          = 0.0f;
        eFieldPolarSum[1]          = 0.0f;
        eFieldPolarSum[2]          = 0.0f;

        gkFieldSum[0]              = 0.0f;
        gkFieldSum[1]              = 0.0f;
        gkFieldSum[2]              = 0.0f;

        if (x == y) // Handle diagonals uniquely at 50% efficiency
        {

            // load coordinates, charge, ...

            loadFixedFieldShared( &(sA[threadIdx.x]), atomI,
                                  atomCoord, labFrameDipole, labFrameQuadrupole,
                                  cAmoebaSim.pDampingFactorAndThole, bornRadii );

            if (!bExclusionFlag)
            {

                // this branch is never exercised since it includes the
                // interaction between atomI and itself which is always excluded

                for (unsigned int j = 0; j < GRID; j++)
                {

                    float ijField[2][3];

                    // load coords, charge, ...

                    loadFixedFieldParticleData( &(psA[j]), &jCoord, jDipole, jQuadrupole, &jBornRadius );

                    calculateFixedEFieldPairIxn_kernel( iCoord,                                          jCoord,
                                                        cAmoebaSim.pDampingFactorAndThole[atomI].x,      psA[j].damp,
                                                        cAmoebaSim.pDampingFactorAndThole[atomI].y,      psA[j].thole,
                                                        &(labFrameDipole[atomI*3]),                      jDipole,
                                                        &(labFrameQuadrupole[atomI*9]),                  jQuadrupole,
                                                        cAmoebaSim.scalingDistanceCutoff,                ijField
#ifdef AMOEBA_DEBUG
                                                , pullBack
#endif
                    );

                    unsigned int match      = (atomI == (y + j)) ? 1 : 0;

                    // add to field at atomI the field due atomJ's charge/dipole/quadrupole

                    eFieldSum[0]           += match ? 0.0f : ijField[0][0];
                    eFieldSum[1]           += match ? 0.0f : ijField[0][1];
                    eFieldSum[2]           += match ? 0.0f : ijField[0][2];

                    eFieldPolarSum[0]      += match ? 0.0f : ijField[0][0];
                    eFieldPolarSum[1]      += match ? 0.0f : ijField[0][1];
                    eFieldPolarSum[2]      += match ? 0.0f : ijField[0][2];

                    // GK field

                    calculateFixedGkFieldPairIxn_kernel( iCoord,                             jCoord,
                                                         &(labFrameDipole[atomI*3]),         jDipole,
                                                         &(labFrameQuadrupole[atomI*9]),     jQuadrupole,
                                                         bornRadii[atomI]*jBornRadius,       ijField
#ifdef AMOEBA_DEBUG
                                                  , pullBack
#endif
                    );

                    gkFieldSum[0]          += match ? 0.0f : ijField[0][0];
                    gkFieldSum[1]          += match ? 0.0f : ijField[0][1];
                    gkFieldSum[2]          += match ? 0.0f : ijField[0][2];

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

                    float ijField[2][3];

                    loadFixedFieldParticleData( &(psA[j]), &jCoord, jDipole, jQuadrupole, &jBornRadius );

                    calculateFixedEFieldPairIxn_kernel( iCoord,                                          jCoord,
                                                        cAmoebaSim.pDampingFactorAndThole[atomI].x,      psA[j].damp,
                                                        cAmoebaSim.pDampingFactorAndThole[atomI].y,      psA[j].thole,
                                                        &(labFrameDipole[atomI*3]),                      jDipole,
                                                        &(labFrameQuadrupole[atomI*9]),                  jQuadrupole,
                                                        cAmoebaSim.scalingDistanceCutoff,                ijField
#ifdef AMOEBA_DEBUG
                                                , pullBack
#endif
                    );


                    float dScaleVal;
                    float pScaleVal;
                    getMaskedDScaleFactor( j, dScaleMask, &dScaleVal );
                    getMaskedPScaleFactor( j, pScaleMask, &pScaleVal );

                    // nan*0.0 = nan not 0.0, so explicitly exclude (atomI == atomJ) contribution
                    // by setting match flag

                    unsigned int match      = (atomI == (y + j)) ? 1 : 0;

                    // add to field at atomI the field due atomJ's charge/dipole/quadrupole

                    eFieldSum[0]           += match ? 0.0f : dScaleVal*ijField[0][0];
                    eFieldSum[1]           += match ? 0.0f : dScaleVal*ijField[0][1];
                    eFieldSum[2]           += match ? 0.0f : dScaleVal*ijField[0][2];

                    eFieldPolarSum[0]      += match ? 0.0f : pScaleVal*ijField[0][0];
                    eFieldPolarSum[1]      += match ? 0.0f : pScaleVal*ijField[0][1];
                    eFieldPolarSum[2]      += match ? 0.0f : pScaleVal*ijField[0][2];

//#ifdef AMOEBA_DEBUG
#if 0
if( 0 && atomI == targetAtom ){
            unsigned int index                 = atomI == targetAtom ? (y + j) : atomI;
            unsigned int pullBackIndex         = 0;
            unsigned int indexI                = 0;
            unsigned int indexJ                = indexI ? 0 : 1;

            debugArray[index].x                = (float) atomI;
            debugArray[index].y                = (float) (y + j);
            debugArray[index].z                = dScaleVal;
            debugArray[index].w                = pScaleVal;
            pullBackIndex                     += 2;

            index                             += cAmoebaSim.paddedNumberOfAtoms;
            debugArray[index].x                = pullBack[pullBackIndex].x;
            debugArray[index].y                = pullBack[pullBackIndex].y;
            debugArray[index].z                = pullBack[pullBackIndex].z;
            pullBackIndex++;

            index                             += cAmoebaSim.paddedNumberOfAtoms;
            debugArray[index].x                = pullBack[pullBackIndex].x;
            debugArray[index].y                = pullBack[pullBackIndex].y;
            debugArray[index].z                = pullBack[pullBackIndex].z;
            pullBackIndex++;

            index                             += cAmoebaSim.paddedNumberOfAtoms;
            debugArray[index].x                = ijField[indexI][0];
            debugArray[index].y                = ijField[indexI][1];
            debugArray[index].z                = ijField[indexI][2];

            index                             += cAmoebaSim.paddedNumberOfAtoms;
            debugArray[index].x                = ijField[indexJ][0];
            debugArray[index].y                = ijField[indexJ][1];
            debugArray[index].z                = ijField[indexJ][2];

            index                             += cAmoebaSim.paddedNumberOfAtoms;
            debugArray[index].x                = match ? 0.0f : dScaleVal*ijField[indexI][0];
            debugArray[index].y                = match ? 0.0f : dScaleVal*ijField[indexI][1];
            debugArray[index].z                = match ? 0.0f : dScaleVal*ijField[indexI][2];

            index                             += cAmoebaSim.paddedNumberOfAtoms;
            debugArray[index].x                = match ? 0.0f : pScaleVal*ijField[indexI][0];
            debugArray[index].y                = match ? 0.0f : pScaleVal*ijField[indexI][1];
            debugArray[index].z                = match ? 0.0f : pScaleVal*ijField[indexI][2];

/*
            index                             += cAmoebaSim.paddedNumberOfAtoms;
            unsigned int mask                  = 1 << j;
            unsigned int pScaleIndex           = (scaleMask.x & mask) ? 1 : 0;
            pScaleIndex                       += (scaleMask.y & mask) ? 2 : 0;
            debugArray[index].x                = (float) pScaleIndex;

            debugArray[index].y                = scaleMask.x & mask ? 1.0f : -1.0f;
            debugArray[index].z                = scaleMask.y & mask ? 1.0f : -1.0f;
            debugArray[index].w                = pScaleVal + 10.0f;
            index                             += cAmoebaSim.paddedNumberOfAtoms;
            debugArray[index].x                = jCoord.x;
            debugArray[index].y                = jCoord.y;
            debugArray[index].z                = jCoord.z;

            index                             += cAmoebaSim.paddedNumberOfAtoms;
            debugArray[index].x                = iCoord.x;
            debugArray[index].y                = iCoord.y;
            debugArray[index].z                = iCoord.z;
*/

}
#endif
                    // GK field

                    calculateFixedGkFieldPairIxn_kernel( iCoord,                                 jCoord,
                                                         &(labFrameDipole[atomI*3]),             jDipole,
                                                         &(labFrameQuadrupole[atomI*9]),         jQuadrupole,
                                                         bornRadii[atomI]*jBornRadius,           ijField
#ifdef AMOEBA_DEBUG
                                                , pullBack
#endif
                    );

                    match                   = (atomI >= cAmoebaSim.numberOfAtoms) || ((y+tj) >= cAmoebaSim.numberOfAtoms) ? 1 : 0;
                    gkFieldSum[0]          += match ? 0.0f : ijField[0][0];
                    gkFieldSum[1]          += match ? 0.0f : ijField[0][1];
                    gkFieldSum[2]          += match ? 0.0f : ijField[0][2];


//#ifdef AMOEBA_DEBUG
#if 0
if( atomI == targetAtom ){
            unsigned int index                 = atomI == targetAtom ? (y + j) : atomI;
            ///unsigned int pullBackIndex         = 0;
            unsigned int indexI                = 0;
            unsigned int indexJ                = indexI ? 0 : 1;

            debugArray[index].x                = (float) atomI;
            debugArray[index].y                = (float) (y + j);
            debugArray[index].z                = bornRadii[atomI];
            debugArray[index].w                = jBornRadius;

/*
            pullBackIndex                     += 2;
            index                             += cAmoebaSim.paddedNumberOfAtoms;
            debugArray[index].x                = pullBack[pullBackIndex].x;
            debugArray[index].y                = pullBack[pullBackIndex].y;
            debugArray[index].z                = pullBack[pullBackIndex].z;
            pullBackIndex++;

            index                             += cAmoebaSim.paddedNumberOfAtoms;
            debugArray[index].x                = pullBack[pullBackIndex].x;
            debugArray[index].y                = pullBack[pullBackIndex].y;
            debugArray[index].z                = pullBack[pullBackIndex].z;
            pullBackIndex++;
*/
            index                             += cAmoebaSim.paddedNumberOfAtoms;
            debugArray[index].x                =  match ? 0.0f :  ijField[indexI][0];
            debugArray[index].y                =  match ? 0.0f :  ijField[indexI][1];
            debugArray[index].z                =  match ? 0.0f :  ijField[indexI][2];

            index                             +=  match ? 0.0f :  cAmoebaSim.paddedNumberOfAtoms;
            debugArray[index].x                =  match ? 0.0f :  ijField[indexJ][0];
            debugArray[index].y                =  match ? 0.0f :  ijField[indexJ][1];
            debugArray[index].z                =  match ? 0.0f :  ijField[indexJ][2];

}
#endif
                }
            }

            // Write results

#ifdef USE_OUTPUT_BUFFER_PER_WARP
            unsigned int offset                 = 3*(x + tgx + warp*cAmoebaSim.paddedNumberOfAtoms);
            load3dArrayBufferPerWarp( offset, eFieldSum,       outputEField );
            load3dArrayBufferPerWarp( offset, eFieldPolarSum,  outputEFieldPolar );
            load3dArrayBufferPerWarp( offset, gkFieldSum,      outputGkField );

#else
            unsigned int offset                 = 3*(x + tgx + (x >> GRIDBITS) * cAmoebaSim.paddedNumberOfAtoms);
            load3dArray( offset, eFieldSum,       outputEField );
            load3dArray( offset, eFieldPolarSum,  outputEFieldPolar );
            load3dArray( offset, gkFieldSum,      outputGkField );

#endif

        }
        else        // 100% utilization
        {
            // Read fixed atom data into registers and GRF
            if (lasty != y)
            {
                // load coordinates, charge, ...

                loadFixedFieldShared( &(sA[threadIdx.x]), (y+tgx),
                                      atomCoord, labFrameDipole, labFrameQuadrupole,
                                      cAmoebaSim.pDampingFactorAndThole, bornRadii );

            }

            // zero shared fields

            zeroFixedFieldParticleSharedField( &(sA[threadIdx.x]) );

            if (!bExclusionFlag)
            {
                for (unsigned int j = 0; j < GRID; j++)
                {

                    float ijField[2][3];
    
                    // load coords, charge, ...
    
                    loadFixedFieldParticleData( &(psA[tj]),  &jCoord, jDipole, jQuadrupole, &jBornRadius );

                    calculateFixedEFieldPairIxn_kernel( iCoord,                                          jCoord,
                                                        cAmoebaSim.pDampingFactorAndThole[atomI].x,      psA[tj].damp,
                                                        cAmoebaSim.pDampingFactorAndThole[atomI].y,      psA[tj].thole,
                                                        &(labFrameDipole[atomI*3]),                      jDipole,
                                                        &(labFrameQuadrupole[atomI*9]),                  jQuadrupole,
                                                        cAmoebaSim.scalingDistanceCutoff,                ijField
#ifdef AMOEBA_DEBUG
                                                 , pullBack
#endif
                    );
    
                    // add to field at atomI the field due atomJ's charge/dipole/quadrupole

                    eFieldSum[0]       += ijField[0][0];
                    eFieldSum[1]       += ijField[0][1];
                    eFieldSum[2]       += ijField[0][2];
        
                    eFieldPolarSum[0]  += ijField[0][0];
                    eFieldPolarSum[1]  += ijField[0][1];
                    eFieldPolarSum[2]  += ijField[0][2];

                    // add to field at atomJ the field due atomI's charge/dipole/quadrupole

                    psA[tj].eField[0]  += ijField[1][0];
                    psA[tj].eField[1]  += ijField[1][1];
                    psA[tj].eField[2]  += ijField[1][2];

                    psA[tj].eFieldP[0] += ijField[1][0];
                    psA[tj].eFieldP[1] += ijField[1][1];
                    psA[tj].eFieldP[2] += ijField[1][2];

// #ifdef AMOEBA_DEBUG
#if 0
if( 0 && (atomI == targetAtom  || (y + tj) == targetAtom) ){
            unsigned int index                 = (atomI == targetAtom) ? (y + tj) : atomI;
            unsigned int indexI                = (atomI == targetAtom) ? 0 : 1;
            unsigned int indexJ                = (atomI == targetAtom) ? 1 : 0;
            //unsigned int pullBackIndex         = 0;

            debugArray[index].x                = (float) atomI;
            debugArray[index].y                = (float) (y + tj);
/*
            debugArray[index].z                = pullBack[pullBackIndex++];
            debugArray[index].w                = pullBack[pullBackIndex++];

            index                             += cAmoebaSim.paddedNumberOfAtoms;
            debugArray[index].x                = pullBack[pullBackIndex++];
            debugArray[index].y                = pullBack[pullBackIndex++];
            debugArray[index].z                = pullBack[pullBackIndex++];

            index                             += cAmoebaSim.paddedNumberOfAtoms;
            debugArray[index].x                = pullBack[pullBackIndex++];
            debugArray[index].y                = pullBack[pullBackIndex++];
            debugArray[index].z                = pullBack[pullBackIndex++];
*/

            index                             += cAmoebaSim.paddedNumberOfAtoms;
            debugArray[index].x                = ijField[indexI][0];
            debugArray[index].y                = ijField[indexI][1];
            debugArray[index].z                = ijField[indexI][2];

            index                             += cAmoebaSim.paddedNumberOfAtoms;
            debugArray[index].x                = ijField[indexJ][0];
            debugArray[index].y                = ijField[indexJ][1];
            debugArray[index].z                = ijField[indexJ][2];

            index                             += cAmoebaSim.paddedNumberOfAtoms;
            debugArray[index].x                = ijField[indexI][0];
            debugArray[index].y                = ijField[indexI][1];
            debugArray[index].z                = ijField[indexI][2];

            index                             += cAmoebaSim.paddedNumberOfAtoms;
            debugArray[index].x                = ijField[indexI][0];
            debugArray[index].y                = ijField[indexI][1];
            debugArray[index].z                = ijField[indexI][2];

#if 0
            index                             += cAmoebaSim.paddedNumberOfAtoms;
            debugArray[index].x                = jCoord.x;
            debugArray[index].y                = jCoord.y;
            debugArray[index].z                = jCoord.z;

            index                             += cAmoebaSim.paddedNumberOfAtoms;
            debugArray[index].x                = iCoord.x;
            debugArray[index].y                = iCoord.y;
            debugArray[index].z                = iCoord.z;

            index                             += cAmoebaSim.paddedNumberOfAtoms;
            unsigned int mask                  = 1 << j;
            unsigned int pScaleIndex           = (scaleMask.x & mask) ? 1 : 0;
            pScaleIndex                       += (scaleMask.y & mask) ? 2 : 0;
            debugArray[index].x                = (float) pScaleIndex;

            debugArray[index].y                = scaleMask.x & mask ? 1.0f : -1.0f;
            debugArray[index].z                = scaleMask.y & mask ? 1.0f : -1.0f;
            debugArray[index].w                = pScaleVal + 10.0f;
#endif
}
#endif
                    // Gk field

                    calculateFixedGkFieldPairIxn_kernel( iCoord,                                          jCoord,
                                                         &(labFrameDipole[atomI*3]),                      jDipole,
                                                         &(labFrameQuadrupole[atomI*9]),                  jQuadrupole,
                                                         bornRadii[atomI]*jBornRadius,                    ijField
#ifdef AMOEBA_DEBUG
                                                , pullBack
#endif
                    );

                    gkFieldSum[0]              += ijField[0][0];
                    gkFieldSum[1]              += ijField[0][1];
                    gkFieldSum[2]              += ijField[0][2];

                    psA[tj].gkField[0]         += ijField[1][0];
                    psA[tj].gkField[1]         += ijField[1][1];
                    psA[tj].gkField[2]         += ijField[1][2];


//#ifdef AMOEBA_DEBUG
#if 0
if( (atomI == targetAtom  || (y + tj) == targetAtom) ){
            unsigned int index                 = (atomI == targetAtom) ? (y + tj) : atomI;
            unsigned int indexI                = (atomI == targetAtom) ? 0 : 1;
            unsigned int indexJ                = (atomI == targetAtom) ? 1 : 0;
            //unsigned int pullBackIndex         = 0;

            debugArray[index].x                = (float) atomI;
            debugArray[index].y                = (float) (y + tj);
            debugArray[index].z                = bornRadii[atomI];
            debugArray[index].w                = jBornRadius;
/*
            debugArray[index].z                = pullBack[pullBackIndex++];
            debugArray[index].w                = pullBack[pullBackIndex++];

            index                             += cAmoebaSim.paddedNumberOfAtoms;
            debugArray[index].x                = pullBack[pullBackIndex++];
            debugArray[index].y                = pullBack[pullBackIndex++];
            debugArray[index].z                = pullBack[pullBackIndex++];

*/

            index                             += cAmoebaSim.paddedNumberOfAtoms;
            debugArray[index].x                = ijField[indexI][0];
            debugArray[index].y                = ijField[indexI][1];
            debugArray[index].z                = ijField[indexI][2];

            index                             += cAmoebaSim.paddedNumberOfAtoms;
            debugArray[index].x                = ijField[indexJ][0];
            debugArray[index].y                = ijField[indexJ][1];
            debugArray[index].z                = ijField[indexJ][2];

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

                    float ijField[2][3];

                    loadFixedFieldParticleData( &(psA[tj]),  &jCoord, jDipole, jQuadrupole, &jBornRadius );
 
                    calculateFixedEFieldPairIxn_kernel( iCoord,                                          jCoord,
                                                        cAmoebaSim.pDampingFactorAndThole[atomI].x,      psA[tj].damp,
                                                        cAmoebaSim.pDampingFactorAndThole[atomI].y,      psA[tj].thole,
                                                        &(labFrameDipole[atomI*3]),                      jDipole,
                                                        &(labFrameQuadrupole[atomI*9]),                  jQuadrupole,
                                                        cAmoebaSim.scalingDistanceCutoff,                ijField
#ifdef AMOEBA_DEBUG
                                                , pullBack
#endif
                    );

                    float dScaleVal;
                    float pScaleVal;
                    getMaskedDScaleFactor( tj, dScaleMask, &dScaleVal );
                    getMaskedPScaleFactor( tj, pScaleMask, &pScaleVal );

                    // add to field at atomI the field due atomJ's charge/dipole/quadrupole

                    eFieldSum[0]         += dScaleVal*ijField[0][0];
                    eFieldSum[1]         += dScaleVal*ijField[0][1];
                    eFieldSum[2]         += dScaleVal*ijField[0][2];

                    eFieldPolarSum[0]    += pScaleVal*ijField[0][0];
                    eFieldPolarSum[1]    += pScaleVal*ijField[0][1];
                    eFieldPolarSum[2]    += pScaleVal*ijField[0][2];

                    // add to field at atomJ the field due atomI's charge/dipole/quadrupole

                    psA[tj].eField[0]    += dScaleVal*ijField[1][0];
                    psA[tj].eField[1]    += dScaleVal*ijField[1][1];
                    psA[tj].eField[2]    += dScaleVal*ijField[1][2];

                    psA[tj].eFieldP[0]   += pScaleVal*ijField[1][0];
                    psA[tj].eFieldP[1]   += pScaleVal*ijField[1][1];
                    psA[tj].eFieldP[2]   += pScaleVal*ijField[1][2];

//#ifdef AMOEBA_DEBUG
#if 0
if( 0 && (atomI == targetAtom || (y + tj) == targetAtom) ){
            unsigned int index                 = (atomI == targetAtom) ? (y + tj) : atomI;
//            unsigned int indexI                = (atomI == targetAtom) ? 0 : 1;
//            unsigned int indexJ                = (atomI == targetAtom) ? 1 : 0;
//            unsigned int pullBackIndex         = 0;

            debugArray[index].x                = (float) atomI;
            debugArray[index].y                = (float) (y + tj);
            debugArray[index].z                = dScaleVal;
            debugArray[index].w                = pScaleVal;
//            pullBackIndex                     += 2;

/*
            index                             += cAmoebaSim.paddedNumberOfAtoms;
            debugArray[index].x                = pullBack[pullBackIndex++];
            debugArray[index].y                = pullBack[pullBackIndex++];
            debugArray[index].z                = pullBack[pullBackIndex++];

            index                             += cAmoebaSim.paddedNumberOfAtoms;
            debugArray[index].x                = pullBack[pullBackIndex++];
            debugArray[index].y                = pullBack[pullBackIndex++];
            debugArray[index].z                = pullBack[pullBackIndex++];

            index                             += cAmoebaSim.paddedNumberOfAtoms;
            debugArray[index].x                = ijField[indexI][0];
            debugArray[index].y                = ijField[indexI][1];
            debugArray[index].z                = ijField[indexI][2];

            index                             += cAmoebaSim.paddedNumberOfAtoms;
            debugArray[index].x                = ijField[indexJ][0];
            debugArray[index].y                = ijField[indexJ][1];
            debugArray[index].z                = ijField[indexJ][2];

            index                             += cAmoebaSim.paddedNumberOfAtoms;
            debugArray[index].x                = dScaleVal*ijField[indexI][0];
            debugArray[index].y                = dScaleVal*ijField[indexI][1];
            debugArray[index].z                = dScaleVal*ijField[indexI][2];

            index                             += cAmoebaSim.paddedNumberOfAtoms;
            debugArray[index].x                = pScaleVal*ijField[indexI][0];
            debugArray[index].y                = pScaleVal*ijField[indexI][1];
            debugArray[index].z                = pScaleVal*ijField[indexI][2];
*/
/*
            index                             += cAmoebaSim.paddedNumberOfAtoms;
            unsigned int mask                  = 1 << j;
            unsigned int pScaleIndex           = (scaleMask.x & mask) ? 1 : 0;
            pScaleIndex                       += (scaleMask.y & mask) ? 2 : 0;
            debugArray[index].x                = (float) pScaleIndex;

            debugArray[index].y                = scaleMask.x & mask ? 1.0f : -1.0f;
            debugArray[index].z                = scaleMask.y & mask ? 1.0f : -1.0f;
            debugArray[index].w                = pScaleVal + 10.0f;

            index                             += cAmoebaSim.paddedNumberOfAtoms;
            debugArray[index].x                = jCoord.x;
            debugArray[index].y                = jCoord.y;
            debugArray[index].z                = jCoord.z;

            index                             += cAmoebaSim.paddedNumberOfAtoms;
            debugArray[index].x                = iCoord.x;
            debugArray[index].y                = iCoord.y;
            debugArray[index].z                = iCoord.z;
*/

}
#endif
                    // GK field

                    calculateFixedGkFieldPairIxn_kernel( iCoord,                             jCoord,
                                                         &(labFrameDipole[atomI*3]),         jDipole,
                                                         &(labFrameQuadrupole[atomI*9]),     jQuadrupole,
                                                         bornRadii[atomI]*jBornRadius,       ijField
#ifdef AMOEBA_DEBUG
                                                , pullBack
#endif
                    );

                    if( (atomI < cAmoebaSim.numberOfAtoms) && ((y+tj) < cAmoebaSim.numberOfAtoms) ){
                        gkFieldSum[0]        += ijField[0][0];
                        gkFieldSum[1]        += ijField[0][1];
                        gkFieldSum[2]        += ijField[0][2];

                        psA[tj].gkField[0]   += ijField[1][0];
                        psA[tj].gkField[1]   += ijField[1][1];
                        psA[tj].gkField[2]   += ijField[1][2];
                    }


//#ifdef AMOEBA_DEBUG
#if 0
if( (atomI == targetAtom || (y + tj) == targetAtom) ){
            unsigned int index                 = (atomI == targetAtom) ? (y + tj) : atomI;
            unsigned int indexI                = (atomI == targetAtom) ? 0 : 1;
            unsigned int indexJ                = (atomI == targetAtom) ? 1 : 0;
//            unsigned int pullBackIndex         = 0;

            debugArray[index].x                = (float) atomI;
            debugArray[index].y                = (float) (y + tj);
            debugArray[index].z                = bornRadii[atomI];
            debugArray[index].w                = jBornRadius;
//            pullBackIndex                     += 2;

/*
            index                             += cAmoebaSim.paddedNumberOfAtoms;
            debugArray[index].x                = pullBack[pullBackIndex++];
            debugArray[index].y                = pullBack[pullBackIndex++];
            debugArray[index].z                = pullBack[pullBackIndex++];

            index                             += cAmoebaSim.paddedNumberOfAtoms;
            debugArray[index].x                = pullBack[pullBackIndex++];
            debugArray[index].y                = pullBack[pullBackIndex++];
            debugArray[index].z                = pullBack[pullBackIndex++];

*/
            index                             += cAmoebaSim.paddedNumberOfAtoms;
            debugArray[index].x                = (atomI == targetAtom || (y + tj) == targetAtom) ?  ijField[indexI][0] : 0.0f;
            debugArray[index].y                = (atomI == targetAtom || (y + tj) == targetAtom) ?  ijField[indexI][1] : 0.0f;
            debugArray[index].z                = (atomI == targetAtom || (y + tj) == targetAtom) ?  ijField[indexI][2] : 0.0f;

            index                             += cAmoebaSim.paddedNumberOfAtoms;
            debugArray[index].x                = (atomI == targetAtom || (y + tj) == targetAtom) ?  ijField[indexJ][0] : 0.0f;
            debugArray[index].y                = (atomI == targetAtom || (y + tj) == targetAtom) ?  ijField[indexJ][1] : 0.0f;
            debugArray[index].z                = (atomI == targetAtom || (y + tj) == targetAtom) ?  ijField[indexJ][2] : 0.0f;

/*
            index                             += cAmoebaSim.paddedNumberOfAtoms;
            debugArray[index].x                = jCoord.x;
            debugArray[index].y                = jCoord.y;
            debugArray[index].z                = jCoord.z;

            index                             += cAmoebaSim.paddedNumberOfAtoms;
            debugArray[index].x                = iCoord.x;
            debugArray[index].y                = iCoord.y;
            debugArray[index].z                = iCoord.z;
*/

}
#endif



                    tj                  = (tj + 1) & (GRID - 1);
                }
            }

            // Write results

#ifdef USE_OUTPUT_BUFFER_PER_WARP

            unsigned int offset                 = 3*(x + tgx + warp * cAmoebaSim.paddedNumberOfAtoms);
            load3dArrayBufferPerWarp( offset, eFieldSum,       outputEField );
            load3dArrayBufferPerWarp( offset, eFieldPolarSum,  outputEFieldPolar );
            load3dArrayBufferPerWarp( offset, gkFieldSum,      outputGkField );


            offset                              = 3*(y + tgx + warp*cAmoebaSim.paddedNumberOfAtoms);
            load3dArrayBufferPerWarp( offset, sA[threadIdx.x].eField,  outputEField );
            load3dArrayBufferPerWarp( offset, sA[threadIdx.x].eFieldP, outputEFieldPolar );
            load3dArrayBufferPerWarp( offset, sA[threadIdx.x].gkField, outputGkField );

#else
            unsigned int offset                 = 3*(x + tgx + (y >> GRIDBITS) * cAmoebaSim.paddedNumberOfAtoms);
            load3dArray( offset, eFieldSum,       outputEField );
            load3dArray( offset, eFieldPolarSum,  outputEFieldPolar );
            load3dArray( offset, gkFieldSum,      outputGkField );

            offset                              = 3*(y + tgx + (x >> GRIDBITS) * cAmoebaSim.paddedNumberOfAtoms);
            load3dArray( offset, sA[threadIdx.x].eField,  outputEField );
            load3dArray( offset, sA[threadIdx.x].eFieldP, outputEFieldPolar );
            load3dArray( offset, sA[threadIdx.x].gkField, outputGkField );
#endif
            lasty = y;
        }

        pos++;
    }
}
