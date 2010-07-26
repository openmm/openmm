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
/*
#if (__CUDA_ARCH__ >= 200)
__launch_bounds__(GF1XX_NONBOND_THREADS_PER_BLOCK, 1)
#elif (__CUDA_ARCH__ >= 130)
__launch_bounds__(GT2XX_NONBOND_THREADS_PER_BLOCK, 1)
#else
__launch_bounds__(G8X_NONBOND_THREADS_PER_BLOCK, 1)
#endif
*/
void METHOD_NAME(kCalculateAmoebaCudaKirkwood, Forces_kernel)(
                            unsigned int* workUnit
#ifdef AMOEBA_DEBUG
                           , float4* debugArray, unsigned int targetAtom
#endif
){

#ifdef AMOEBA_DEBUG
    float4 pullBack[20];
#endif

    extern __shared__ KirkwoodParticle sA[];

    unsigned int totalWarps      = gridDim.x*blockDim.x/GRID;
    unsigned int warp            = (blockIdx.x*blockDim.x+threadIdx.x)/GRID;
    unsigned int numWorkUnits    = cSim.pInteractionCount[0];
    unsigned int pos             = warp*numWorkUnits/totalWarps;
    unsigned int end             = (warp+1)*numWorkUnits/totalWarps;
    unsigned int lasty           = 0xFFFFFFFF;

    // pWorkArray_3_1 == force
    // pWorkArray_3_2 == torque

    // pWorkArray_1_1 == dBorn
    // pWorkArray_1_2 == dBornPolar

    float4 jCoord;
    float jDipole[3];     
    float jQuadrupole[9];     
    float jInducedDipole[3];     
    float jInducedDipolePolar[3];     
    float jBornRadius;     

    float energySum = 0.0f;     

    while (pos < end)
    {

        unsigned int x;
        unsigned int y;
        bool bExclusionFlag;

        // extract cell coordinates

        decodeCell( workUnit[pos], &x, &y, &bExclusionFlag );

        unsigned int tgx              = threadIdx.x & (GRID - 1);
        unsigned int tbx              = threadIdx.x - tgx;
        unsigned int tj               = tgx;

        KirkwoodParticle* psA         = &sA[tbx];

        unsigned int atomI            = x + tgx;
        float4 iCoord                 =  cSim.pPosq[atomI];

        float forceSum[3];
        float torqueSum[3];
        float dBornSum;
        float dBornPolarSum;

        forceSum[0]                   = 0.0f;
        forceSum[1]                   = 0.0f;
        forceSum[2]                   = 0.0f;

        torqueSum[0]                  = 0.0f;
        torqueSum[1]                  = 0.0f;
        torqueSum[2]                  = 0.0f;

        dBornSum                      = 0.0f;
        dBornPolarSum                 = 0.0f;

        // handle diagonals uniquely at 50% efficiency

        if (x == y) 
        {

            // load shared data

            loadKirkwoodShared( &(sA[threadIdx.x]), atomI,
                                cSim.pPosq, cAmoebaSim.pLabFrameDipole, cAmoebaSim.pLabFrameQuadrupole,
                                cAmoebaSim.pInducedDipoleS,  cAmoebaSim.pInducedDipolePolarS, cSim.pBornRadii );

            // this branch is never exercised since it includes the
            // interaction between atomI and itself which is always excluded

            for (unsigned int j = 0; j < GRID; j++)
            {

                float force[3];
                float torque[2][3];
                float dBorn[2];
                float dBornPolar[2];
                float energy;

                unsigned int atomJ            = y + j;
                unsigned int sameAtom         = atomI == atomJ ? 1 : 0;

                // load coords, charge, ...

                loadKirkwoodData( &(psA[j]), &jCoord, jDipole, jQuadrupole,
                                  jInducedDipole, jInducedDipolePolar, &jBornRadius );

                calculateKirkwoodPairIxn_kernel( sameAtom,
                                                 iCoord,                                     jCoord,
                                                 &(cAmoebaSim.pLabFrameDipole[3*atomI]),     jDipole,
                                                 &(cAmoebaSim.pLabFrameQuadrupole[9*atomI]), jQuadrupole,
                                                 &(cAmoebaSim.pInducedDipoleS[3*atomI]),     jInducedDipole,
                                                 &(cAmoebaSim.pInducedDipolePolarS[3*atomI]),jInducedDipolePolar,
                                                 cSim.pBornRadii[atomI],                     jBornRadius,
                                                 force, torque, dBorn, dBornPolar, &energy
#ifdef AMOEBA_DEBUG
                                          , pullBack
#endif
                                        );

                unsigned int mask       =  ( (atomI >= cAmoebaSim.numberOfAtoms) || (atomJ >= cAmoebaSim.numberOfAtoms) ) ? 0 : 1;

                // torques include i == j contribution

                torqueSum[0]           += mask ? torque[0][0]  : 0.0f;
                torqueSum[1]           += mask ? torque[0][1]  : 0.0f;
                torqueSum[2]           += mask ? torque[0][2]  : 0.0f;

                dBornSum               += mask ? dBorn[0]      : 0.0f;
                dBornPolarSum          += mask ? dBornPolar[0] : 0.0f;
                energySum              += mask ? 0.5f*energy   : 0.0f;

                // add to field at atomI the field due atomJ's charge/dipole/quadrupole

                mask                    =  (atomI == atomJ) ? 0 : mask;

                forceSum[0]            += mask ? force[0]      : 0.0f;
                forceSum[1]            += mask ? force[1]      : 0.0f;
                forceSum[2]            += mask ? force[2]      : 0.0f;


#ifdef AMOEBA_DEBUG
if( atomI == targetAtom  || atomJ == targetAtom ){

        unsigned int index                 = (atomI == targetAtom) ? atomJ : atomI;
        unsigned int indexI                = (atomI == targetAtom) ? 0 : 1;
        unsigned int indexJ                = (atomI == targetAtom) ? 1 : 0;
        //float forceSign                    = (atomI == targetAtom) ? 1.0f : -1.0f;

        debugArray[index].x                = (float) atomI;
        debugArray[index].y                = (float) atomJ;
        debugArray[index].z                = (float) (mask + 1);
        //debugArray[index].z                = cSim.pBornRadii[atomI];
        //debugArray[index].z                = energy;
        //debugArray[index].w                = (float) (blockIdx.x*blockDim.x+threadIdx.x);
        debugArray[index].w                = jBornRadius;

        index = debugAccumulate( index, debugArray, force,          mask, 1.0f );

        mask  =  ( (atomI >= cAmoebaSim.numberOfAtoms) || (atomJ >= cAmoebaSim.numberOfAtoms) ) ? 0 : 1;
        index = debugAccumulate( index, debugArray, torque[indexI], mask, 2.0f );
        index = debugAccumulate( index, debugArray, torque[indexJ], mask, 3.0f );

        index                             += cAmoebaSim.paddedNumberOfAtoms;
        debugArray[index].x                = pullBack[0].x;
        debugArray[index].y                = pullBack[0].y;
        debugArray[index].z                = pullBack[0].z;
        debugArray[index].w                = pullBack[0].w;

        index                             += cAmoebaSim.paddedNumberOfAtoms;
        debugArray[index].x                = pullBack[1].x;
        debugArray[index].y                = pullBack[1].y;
        debugArray[index].z                = pullBack[1].z;
        debugArray[index].w                = pullBack[1].w;

        index                             += cAmoebaSim.paddedNumberOfAtoms;
        debugArray[index].x                = pullBack[2].x;
        debugArray[index].y                = pullBack[2].y;
        debugArray[index].z                = pullBack[2].z;
        debugArray[index].w                = pullBack[2].w;
/*
        index                             += cAmoebaSim.paddedNumberOfAtoms;
        debugArray[index].x                = dBorn[0];
        debugArray[index].y                = dBornPolar[0];
        debugArray[index].z                = dBorn[1];
        debugArray[index].w                = dBornPolar[1];
*/
}
#endif

            } // end of j-loop

            // write results

#ifdef USE_OUTPUT_BUFFER_PER_WARP
            float  of;
            unsigned int offset                 = x + tgx + warp*cAmoebaSim.paddedNumberOfAtoms;

            of                                  = cAmoebaSim.pWorkArray_1_1[offset];
            of                                 += dBornSum;
            cAmoebaSim.pWorkArray_1_1[offset]   = of;

            of                                  = cAmoebaSim.pWorkArray_1_2[offset];
            of                                 += dBornPolarSum;
            cAmoebaSim.pWorkArray_1_2[offset]   = of;


            offset                             *= 3;

            load3dArrayBufferPerWarp( offset, forceSum,  cAmoebaSim.pWorkArray_3_1 );
            load3dArrayBufferPerWarp( offset, torqueSum, cAmoebaSim.pWorkArray_3_2 );

#else
            unsigned int offset                 = x + tgx + (x >> GRIDBITS) * cAmoebaSim.paddedNumberOfAtoms;

            cAmoebaSim.pWorkArray_1_1[offset]   = dBornSum;
            cAmoebaSim.pWorkArray_1_2[offset]   = dBornPolarSum;

            offset                             *= 3;

            load3dArray( offset, forceSum,  cAmoebaSim.pWorkArray_3_1 );
            load3dArray( offset, torqueSum, cAmoebaSim.pWorkArray_3_2 );

#endif

        }
        else        // 100% utilization
        {
            // Read fixed atom data into registers and GRF

            if (lasty != y)
            {
                // load shared data

                loadKirkwoodShared( &(sA[threadIdx.x]), (y+tgx),
                                    cSim.pPosq, cAmoebaSim.pLabFrameDipole, cAmoebaSim.pLabFrameQuadrupole,
                                    cAmoebaSim.pInducedDipoleS, cAmoebaSim.pInducedDipolePolarS, cSim.pBornRadii);

            }

            // zero j-atom output fields

            zeroKirkwoodParticleSharedField( &(sA[threadIdx.x]) );

            for (unsigned int j = 0; j < GRID; j++)
            {

                float force[3];
                float torque[2][3];

                float dBorn[2];
                float dBornPolar[2];
                float energy;

                unsigned int atomJ    = y + tj;
                unsigned int sameAtom = 0;

                // load coords, charge, ...

                loadKirkwoodData( &(psA[tj]), &jCoord, jDipole, jQuadrupole,
                                  jInducedDipole, jInducedDipolePolar, &jBornRadius );

                calculateKirkwoodPairIxn_kernel( sameAtom, 
                                                 iCoord,                                       jCoord,
                                                 &(cAmoebaSim.pLabFrameDipole[3*atomI]),       jDipole,
                                                 &(cAmoebaSim.pLabFrameQuadrupole[9*atomI]),   jQuadrupole,
                                                 &(cAmoebaSim.pInducedDipoleS[3*atomI]),       jInducedDipole,
                                                 &(cAmoebaSim.pInducedDipolePolarS[3*atomI]),  jInducedDipolePolar,
                                                 cSim.pBornRadii[atomI],                       jBornRadius,
                                                 force, torque, dBorn, dBornPolar, &energy
#ifdef AMOEBA_DEBUG
                                          , pullBack
#endif
                                        );

                unsigned int mask          =  ( (atomI >= cAmoebaSim.numberOfAtoms) || ( atomJ >= cAmoebaSim.numberOfAtoms) ) ? 0 : 1;

                // add force and torque to atom I due atom J

                forceSum[0]               += mask ? force[0]      : 0.0f;
                forceSum[1]               += mask ? force[1]      : 0.0f;
                forceSum[2]               += mask ? force[2]      : 0.0f;

                torqueSum[0]              += mask ? torque[0][0]  : 0.0f;
                torqueSum[1]              += mask ? torque[0][1]  : 0.0f;
                torqueSum[2]              += mask ? torque[0][2]  : 0.0f;

                dBornSum                  += mask ? dBorn[0]      : 0.0f;
                dBornPolarSum             += mask ? dBornPolar[0] : 0.0f;
                energySum                 += mask ? energy        : 0.0f;

                // add force and torque to atom J due atom I

                psA[tj].force[0]          -= mask ? force[0]      : 0.0f;
                psA[tj].force[1]          -= mask ? force[1]      : 0.0f;
                psA[tj].force[2]          -= mask ? force[2]      : 0.0f;

                psA[tj].torque[0]         += mask ? torque[1][0]  : 0.0f;
                psA[tj].torque[1]         += mask ? torque[1][1]  : 0.0f;
                psA[tj].torque[2]         += mask ? torque[1][2]  : 0.0f;

                psA[tj].dBornRadius       += mask ? dBorn[1]      : 0.0f;
                psA[tj].dBornRadiusPolar  += mask ? dBornPolar[1] : 0.0f;

#ifdef AMOEBA_DEBUG
if( atomI == targetAtom || atomJ == targetAtom ){
        unsigned int index                 = (atomI == targetAtom) ? atomJ : atomI;
        unsigned int indexI                = (atomI == targetAtom) ? 0 : 1;
        unsigned int indexJ                = (atomI == targetAtom) ? 1 : 0;
        // float forceSign                    = (atomI == targetAtom) ? 1.0f : -1.0f;

        debugArray[index].x                = (float) atomI;
        debugArray[index].y                = (float) atomJ;
        debugArray[index].z                = (float) (mask+1);
        //debugArray[index].z                = cSim.pBornRadii[atomI];
        //debugArray[index].z                = energy;
        debugArray[index].w                = (float) (blockIdx.x*blockDim.x+threadIdx.x);

        index = debugAccumulate( index, debugArray, force,          mask, -1.0f );
        index = debugAccumulate( index, debugArray, torque[indexI], mask, -2.0f );
        index = debugAccumulate( index, debugArray, torque[indexJ], mask, -3.0f );

        index                             += cAmoebaSim.paddedNumberOfAtoms;
        debugArray[index].x                = pullBack[0].x;
        debugArray[index].y                = pullBack[0].y;
        debugArray[index].z                = pullBack[0].z;
        debugArray[index].w                = -1.0f;

        index                             += cAmoebaSim.paddedNumberOfAtoms;
        debugArray[index].x                = pullBack[1].x;
        debugArray[index].y                = pullBack[1].y;
        debugArray[index].z                = pullBack[1].z;
        debugArray[index].w                = pullBack[1].w;

        index                             += cAmoebaSim.paddedNumberOfAtoms;
        debugArray[index].x                = pullBack[2].x;
        debugArray[index].y                = pullBack[2].y;
        debugArray[index].z                = pullBack[2].z;
        debugArray[index].w                = pullBack[2].w;
/*
        index                             += cAmoebaSim.paddedNumberOfAtoms;
        debugArray[index].x                = dBorn[0];
        debugArray[index].y                = dBornPolar[0];
        debugArray[index].z                = dBorn[1];
        debugArray[index].w                = dBornPolar[1];
*/

}
#endif
//#ifdef AMOEBA_DEBUG
#if 0
if( mask || !mask ){
        unsigned int index                 = atomJ + atomI*cAmoebaSim.paddedNumberOfAtoms;
        debugArray[index].x                = (float) atomI;
        debugArray[index].y                = (float) atomJ;
        debugArray[index].z                = energy;
        debugArray[index].w                = jBornRadius;
}
#endif

                tj                  = (tj + 1) & (GRID - 1);

            }

            // Write results

#ifdef USE_OUTPUT_BUFFER_PER_WARP

            float of;
            unsigned int offset                 = x + tgx + warp*cAmoebaSim.paddedNumberOfAtoms;

            of                                  = cAmoebaSim.pWorkArray_1_1[offset];
            of                                 += dBornSum;
            cAmoebaSim.pWorkArray_1_1[offset]   = of;

            of                                  = cAmoebaSim.pWorkArray_1_2[offset];
            of                                 += dBornPolarSum;
            cAmoebaSim.pWorkArray_1_2[offset]   = of;

            offset                             *= 3;

            load3dArrayBufferPerWarp( offset, forceSum,  cAmoebaSim.pWorkArray_3_1 );
            load3dArrayBufferPerWarp( offset, torqueSum, cAmoebaSim.pWorkArray_3_2 );

            offset                              = y + tgx + warp*cAmoebaSim.paddedNumberOfAtoms;

            of                                  = cAmoebaSim.pWorkArray_1_1[offset];
            of                                 += dBornSum;
            cAmoebaSim.pWorkArray_1_1[offset]   = of;

            of                                  = cAmoebaSim.pWorkArray_1_2[offset];
            of                                 += dBornPolarSum;
            cAmoebaSim.pWorkArray_1_2[offset]   = of;

            offset                             *= 3;

            load3dArrayBufferPerWarp( offset, sA[threadIdx.x].force,  cAmoebaSim.pWorkArray_3_1 );
            load3dArrayBufferPerWarp( offset, sA[threadIdx.x].torque, cAmoebaSim.pWorkArray_3_2 );
#else
            unsigned int offset                 = x + tgx + (y >> GRIDBITS) * cAmoebaSim.paddedNumberOfAtoms;

            cAmoebaSim.pWorkArray_1_1[offset]   = dBornSum;
            cAmoebaSim.pWorkArray_1_2[offset]   = dBornPolarSum;

            offset                             *= 3;

            load3dArray( offset, forceSum,  cAmoebaSim.pWorkArray_3_1 );
            load3dArray( offset, torqueSum, cAmoebaSim.pWorkArray_3_2 );


            offset                              = y + tgx + (x >> GRIDBITS) * cAmoebaSim.paddedNumberOfAtoms;

            cAmoebaSim.pWorkArray_1_1[offset]   = sA[threadIdx.x].dBornRadius;
            cAmoebaSim.pWorkArray_1_2[offset]   = sA[threadIdx.x].dBornRadiusPolar;

            offset                             *= 3;

            load3dArray( offset, sA[threadIdx.x].force,  cAmoebaSim.pWorkArray_3_1 );
            load3dArray( offset, sA[threadIdx.x].torque, cAmoebaSim.pWorkArray_3_2 );

#endif
            lasty = y;
        }
        pos++;
    }
    cSim.pEnergy[blockIdx.x*blockDim.x+threadIdx.x] += energySum;
}
