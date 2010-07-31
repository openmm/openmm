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
void METHOD_NAME(kCalculateAmoebaCudaElectrostatic, Forces_kernel)(
                            unsigned int* workUnit,
                            float4* atomCoord,
                            float* labFrameDipole,
                            float* labFrameQuadrupole,
                            float* inducedDipole,
                            float* inducedDipolePolar,
                            float* outputForce,
                            float* outputTorque

#ifdef AMOEBA_DEBUG
                           , float4* debugArray, unsigned int targetAtom
#endif
){

#ifdef AMOEBA_DEBUG
    float4 pullBack[20];
#endif

    extern __shared__ ElectrostaticParticle sA[];

    unsigned int totalWarps      = gridDim.x*blockDim.x/GRID;
    unsigned int warp            = (blockIdx.x*blockDim.x+threadIdx.x)/GRID;
    unsigned int numWorkUnits    = cSim.pInteractionCount[0];
    unsigned int pos             = warp*numWorkUnits/totalWarps;
    unsigned int end             = (warp+1)*numWorkUnits/totalWarps;
    unsigned int lasty           = 0xFFFFFFFF;
    float totalEnergy            = 0.0f;     

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
        loadElectrostaticShared(&localParticle, atomI,
                                atomCoord, labFrameDipole, labFrameQuadrupole,
                                inducedDipole, inducedDipolePolar, cAmoebaSim.pDampingFactorAndThole );

        localParticle.force[0]                   = 0.0f;
        localParticle.force[1]                   = 0.0f;
        localParticle.force[2]                   = 0.0f;

        localParticle.torque[0]                  = 0.0f;
        localParticle.torque[1]                  = 0.0f;
        localParticle.torque[2]                  = 0.0f;

        scalingFactors[PScaleIndex]   = 1.0f;
        scalingFactors[DScaleIndex]   = 1.0f;
        scalingFactors[UScaleIndex]   = 1.0f;
        scalingFactors[MScaleIndex]   = 1.0f;

        if (x == y) // Handle diagonals uniquely at 50% efficiency
        {

            // load shared data

           loadElectrostaticShared( &(sA[threadIdx.x]), atomI,
                                    atomCoord, labFrameDipole, labFrameQuadrupole,
                                    inducedDipole, inducedDipolePolar, cAmoebaSim.pDampingFactorAndThole );

            if (!bExclusionFlag)
            {

                // this branch is never exercised since it includes the
                // interaction between atomI and itself which is always excluded

                for (unsigned int j = 0; j < GRID; j++)
                {

                    float force[3];
                    float torque[2][3];

                    float energy;
                    calculateElectrostaticPairIxn_kernel( localParticle,                              psA[j],
                                                          cAmoebaSim.scalingDistanceCutoff,           scalingFactors,
                                                          force, torque, &energy
#ifdef AMOEBA_DEBUG
, pullBack
#endif
 );

                    unsigned int mask       =  ( (atomI == (y + j)) || (atomI >= cAmoebaSim.numberOfAtoms) || ((y+j) >= cAmoebaSim.numberOfAtoms) ) ? 0 : 1;

                    // add to field at atomI the field due atomJ's charge/dipole/quadrupole

                    localParticle.force[0]            += mask ? force[0]     : 0.0f;
                    localParticle.force[1]            += mask ? force[1]     : 0.0f;
                    localParticle.force[2]            += mask ? force[2]     : 0.0f;

                    localParticle.torque[0]           += mask ? torque[0][0] : 0.0f;
                    localParticle.torque[1]           += mask ? torque[0][1] : 0.0f;
                    localParticle.torque[2]           += mask ? torque[0][2] : 0.0f;
 
                    totalEnergy            += mask ? 0.5*energy   : 0.0f;
                }

            }
            else  // bExclusion
            {
                unsigned int xi       = x >> GRIDBITS;
                unsigned int cell     = xi + xi*cAmoebaSim.paddedNumberOfAtoms/GRID-xi*(xi+1)/2;
                int  dScaleMask       = cAmoebaSim.pD_ScaleIndices[cAmoebaSim.pScaleIndicesIndex[cell]+tgx];
                int2 pScaleMask       = cAmoebaSim.pP_ScaleIndices[cAmoebaSim.pScaleIndicesIndex[cell]+tgx];
                int2 mScaleMask       = cAmoebaSim.pM_ScaleIndices[cAmoebaSim.pScaleIndicesIndex[cell]+tgx];

                for (unsigned int j = 0; j < GRID; j++)
                {

                    float force[3];
                    float torque[2][3];

                    unsigned int atomJ = y + j;

                    // set scale factors

                    getMaskedDScaleFactor( j, dScaleMask, scalingFactors + DScaleIndex );
                    getMaskedPScaleFactor( j, pScaleMask, scalingFactors + PScaleIndex );
                    getMaskedMScaleFactor( j, mScaleMask, scalingFactors + MScaleIndex );

                    // force

                    float energy;
                    calculateElectrostaticPairIxn_kernel( localParticle,                              psA[j],
                                                          cAmoebaSim.scalingDistanceCutoff,           scalingFactors,
                                                          force, torque, &energy
#ifdef AMOEBA_DEBUG
, pullBack
#endif
 );

                    // nan*0.0 = nan not 0.0, so explicitly exclude (atomI == atomJ) contribution
                    // by setting match flag

                    unsigned int mask       =  ( (atomI == atomJ) || (atomI >= cAmoebaSim.numberOfAtoms) || (atomJ >= cAmoebaSim.numberOfAtoms) ) ? 0 : 1;

                    // add to field at atomI the field due atomJ's charge/dipole/quadrupole

                    localParticle.force[0]            += mask ? force[0]     : 0.0f;
                    localParticle.force[1]            += mask ? force[1]     : 0.0f;
                    localParticle.force[2]            += mask ? force[2]     : 0.0f;

                    localParticle.torque[0]           += mask ? torque[0][0] : 0.0f;
                    localParticle.torque[1]           += mask ? torque[0][1] : 0.0f;
                    localParticle.torque[2]           += mask ? torque[0][2] : 0.0f;

                    totalEnergy            += mask ? 0.5*energy   : 0.0f;

#ifdef AMOEBA_DEBUG
if( atomI == targetAtom ){
            unsigned int index                 = (atomI == targetAtom) ? atomJ : atomI;

            debugArray[index].x                = (float) atomI;
            debugArray[index].y                = (float) atomJ;
            debugArray[index].z                = 1.0f;
            debugArray[index].w                = (float) y;
/*
            index                             += cAmoebaSim.paddedNumberOfAtoms;
            debugArray[index].x                = mask ? force[0]     : 0.0f;
            debugArray[index].y                = mask ? force[1]     : 0.0f;
            debugArray[index].z                = mask ? force[2]     : 0.0f;


            index                             += cAmoebaSim.paddedNumberOfAtoms;
            debugArray[index].x                = mask ? torque[0][0] : 0.0f;
            debugArray[index].y                = mask ? torque[0][1] : 0.0f;
            debugArray[index].z                = mask ? torque[0][2] : 0.0f;
*/

            index                             += cAmoebaSim.paddedNumberOfAtoms;
            int pullIndex                      = 0;
            debugArray[index].x                = pullBack[pullIndex].x;
            debugArray[index].y                = pullBack[pullIndex].y;
            debugArray[index].z                = pullBack[pullIndex].z;
            debugArray[index].w                = pullBack[pullIndex].w;

            index                             += cAmoebaSim.paddedNumberOfAtoms;
            pullIndex++;
            debugArray[index].x                = pullBack[pullIndex].x;
            debugArray[index].y                = pullBack[pullIndex].y;
            debugArray[index].z                = pullBack[pullIndex].z;
            debugArray[index].w                = pullBack[pullIndex].w;

            index                             += cAmoebaSim.paddedNumberOfAtoms;
            pullIndex++;
            debugArray[index].x                = pullBack[pullIndex].x;
            debugArray[index].y                = pullBack[pullIndex].y;
            debugArray[index].z                = pullBack[pullIndex].z;
            debugArray[index].w                = pullBack[pullIndex].w;

            index                             += cAmoebaSim.paddedNumberOfAtoms;
            pullIndex++;
            debugArray[index].x                = pullBack[pullIndex].x;
            debugArray[index].y                = pullBack[pullIndex].y;
            debugArray[index].z                = pullBack[pullIndex].z;
            debugArray[index].w                = pullBack[pullIndex].w;

}
#endif

                }
            }

            // Write results

#ifdef USE_OUTPUT_BUFFER_PER_WARP
            float  of;
            unsigned int offset                 = 3*(x + tgx + warp*cAmoebaSim.paddedNumberOfAtoms);
            of                                  = outputForce[offset];
            of                                 += localParticle.force[0];
            outputForce[offset]                 = of;

            of                                  = outputForce[offset+1];
            of                                 += localParticle.force[1];
            outputForce[offset+1]               = of;

            of                                  = outputForce[offset+2];
            of                                 += localParticle.force[2];
            outputForce[offset+2]               = of;

            of                                  = outputTorque[offset];
            of                                 += localParticle.torque[0];
            outputTorque[offset]                = of;

            of                                  = outputTorque[offset+1];
            of                                 += localParticle.torque[1];
            outputTorque[offset+1]              = of;

            of                                  = outputTorque[offset+2];
            of                                 += localParticle.torque[2];
            outputTorque[offset+2]              = of;

#else
            unsigned int offset                 = 3*(x + tgx + (x >> GRIDBITS) * cAmoebaSim.paddedNumberOfAtoms);
            outputForce[offset]                 = localParticle.force[0];
            outputForce[offset+1]               = localParticle.force[1];
            outputForce[offset+2]               = localParticle.force[2];

            outputTorque[offset]                = localParticle.torque[0];
            outputTorque[offset+1]              = localParticle.torque[1];
            outputTorque[offset+2]              = localParticle.torque[2];
#endif

        }
        else        // 100% utilization
        {
            // Read fixed atom data into registers and GRF

            if (lasty != y)
            {
                // load shared data

               loadElectrostaticShared( &(sA[threadIdx.x]), (y+tgx),
                                        atomCoord, labFrameDipole, labFrameQuadrupole,
                                        inducedDipole, inducedDipolePolar, cAmoebaSim.pDampingFactorAndThole );

            }

            sA[threadIdx.x].force[0]     = 0.0f;
            sA[threadIdx.x].force[1]     = 0.0f;
            sA[threadIdx.x].force[2]     = 0.0f;

            sA[threadIdx.x].torque[0]    = 0.0f;
            sA[threadIdx.x].torque[1]    = 0.0f;
            sA[threadIdx.x].torque[2]    = 0.0f;

            if (!bExclusionFlag)
            {
                for (unsigned int j = 0; j < GRID; j++)
                {

                    float force[3];
                    float torque[2][3];

                    unsigned int atomJ = y + tj;

                    float energy;
                    calculateElectrostaticPairIxn_kernel( localParticle,                              psA[tj],
                                                                 cAmoebaSim.scalingDistanceCutoff,           scalingFactors,
                                                                 force, torque, &energy
#ifdef AMOEBA_DEBUG
, pullBack
#endif
 );
       
                           unsigned int mask       =  ( (atomI >= cAmoebaSim.numberOfAtoms) || ( atomJ >= cAmoebaSim.numberOfAtoms) ) ? 0 : 1;
       
                           // add force and torque to atom I due atom J
       
                           localParticle.force[0]         += mask ? force[0]      : 0.0f;
                           localParticle.force[1]         += mask ? force[1]      : 0.0f;
                           localParticle.force[2]         += mask ? force[2]      : 0.0f;
       
                           localParticle.torque[0]        += mask ? torque[0][0]  : 0.0f;
                           localParticle.torque[1]        += mask ? torque[0][1]  : 0.0f;
                           localParticle.torque[2]        += mask ? torque[0][2]  : 0.0f;

                           totalEnergy         += mask ? energy        : 0.0f;
           
                           // add force and torque to atom J due atom I
       
                           psA[tj].force[0]     -= mask ?  force[0]     : 0.0f;
                           psA[tj].force[1]     -= mask ?  force[1]     : 0.0f;
                           psA[tj].force[2]     -= mask ?  force[2]     : 0.0f;
       
                           psA[tj].torque[0]    += mask ?  torque[1][0] : 0.0f;
                           psA[tj].torque[1]    += mask ?  torque[1][1] : 0.0f;
                           psA[tj].torque[2]    += mask ?  torque[1][2] : 0.0f;
       
       
#ifdef AMOEBA_DEBUG
if( atomI == targetAtom  || atomJ == targetAtom ){
        unsigned int index                 = (atomI == targetAtom) ? atomJ : atomI;
        unsigned int indexI                = (atomI == targetAtom) ? 0 : 1;
        unsigned int indexJ                = (atomI == targetAtom) ? 1 : 0;
        float forceSign                    = (atomI == targetAtom) ? 1.0f : -1.0f;

        debugArray[index].x                = (float) atomI;
        debugArray[index].y                = (float) atomJ;
        debugArray[index].z                = 2.0f;
        debugArray[index].w                = (float) y;
/*
        index                             += cAmoebaSim.paddedNumberOfAtoms;
        debugArray[index].x                = mask ? forceSign*force[0]     : 0.0f;
        debugArray[index].y                = mask ? forceSign*force[1]     : 0.0f;
        debugArray[index].z                = mask ? forceSign*force[2]     : 0.0f;


        index                             += cAmoebaSim.paddedNumberOfAtoms;
        debugArray[index].x                = mask ? torque[indexI][0] : 0.0f;
        debugArray[index].y                = mask ? torque[indexI][1] : 0.0f;
        debugArray[index].z                = mask ? torque[indexI][2] : 0.0f;

        index                             += cAmoebaSim.paddedNumberOfAtoms;
        debugArray[index].x                = mask ? torque[indexJ][0] : 0.0f;
        debugArray[index].y                = mask ? torque[indexJ][1] : 0.0f;
        debugArray[index].z                = mask ? torque[indexJ][2] : 0.0f;
*/
        index                             += cAmoebaSim.paddedNumberOfAtoms;
        int pullIndex                      = 0;
        debugArray[index].x                = pullBack[pullIndex].x;
        debugArray[index].y                = pullBack[pullIndex].y;
        debugArray[index].z                = pullBack[pullIndex].z;
        debugArray[index].w                = pullBack[pullIndex].w;

        index                             += cAmoebaSim.paddedNumberOfAtoms;
        pullIndex++;
        debugArray[index].x                = pullBack[pullIndex].x;
        debugArray[index].y                = pullBack[pullIndex].y;
        debugArray[index].z                = pullBack[pullIndex].z;
        debugArray[index].w                = pullBack[pullIndex].w;

        index                             += cAmoebaSim.paddedNumberOfAtoms;
        pullIndex++;
        debugArray[index].x                = pullBack[pullIndex].x;
        debugArray[index].y                = pullBack[pullIndex].y;
        debugArray[index].z                = pullBack[pullIndex].z;
        debugArray[index].w                = pullBack[pullIndex].w;
       
        index                             += cAmoebaSim.paddedNumberOfAtoms;
        pullIndex++;
        debugArray[index].x                = pullBack[pullIndex].x;
        debugArray[index].y                = pullBack[pullIndex].y;
        debugArray[index].z                = pullBack[pullIndex].z;
        debugArray[index].w                = pullBack[pullIndex].w;
       
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
                       int2 mScaleMask   = cAmoebaSim.pM_ScaleIndices[cAmoebaSim.pScaleIndicesIndex[cell]+tgx];
       
                       for (unsigned int j = 0; j < GRID; j++)
                       {
       
                           float force[3];
                           float torque[2][3];
       
                           unsigned int atomJ = y + tj;
       
                           // set scale factors
       
                           getMaskedDScaleFactor( tj, dScaleMask, scalingFactors + DScaleIndex );
                           getMaskedPScaleFactor( tj, pScaleMask, scalingFactors + PScaleIndex );
                           getMaskedMScaleFactor( tj, mScaleMask, scalingFactors + MScaleIndex );
       
                           // force
       
                    float energy;
                    calculateElectrostaticPairIxn_kernel( localParticle,                              psA[tj],
                                                          cAmoebaSim.scalingDistanceCutoff,           scalingFactors,
                                                          force, torque, &energy
#ifdef AMOEBA_DEBUG
, pullBack
#endif
 );

                    // check if atoms out-of-bounds

                    unsigned int mask    =  ( (atomI >= cAmoebaSim.numberOfAtoms) || (atomJ >= cAmoebaSim.numberOfAtoms) ) ? 0 : 1;

                    // add force and torque to atom I due atom J

                    localParticle.force[0]         += mask ? force[0]      : 0.0f;
                    localParticle.force[1]         += mask ? force[1]      : 0.0f;
                    localParticle.force[2]         += mask ? force[2]      : 0.0f;

                    localParticle.torque[0]        += mask ? torque[0][0]  : 0.0f;
                    localParticle.torque[1]        += mask ? torque[0][1]  : 0.0f;
                    localParticle.torque[2]        += mask ? torque[0][2]  : 0.0f;
    
                    totalEnergy         += mask ? energy        : 0.0f;

                    // add force and torque to atom J due atom I

                    psA[tj].force[0]     -= mask ?  force[0]     : 0.0f;
                    psA[tj].force[1]     -= mask ?  force[1]     : 0.0f;
                    psA[tj].force[2]     -= mask ?  force[2]     : 0.0f;

                    psA[tj].torque[0]    += mask ?  torque[1][0] : 0.0f;
                    psA[tj].torque[1]    += mask ?  torque[1][1] : 0.0f;
                    psA[tj].torque[2]    += mask ?  torque[1][2] : 0.0f;


#ifdef AMOEBA_DEBUG
if( atomI == targetAtom  || atomJ == targetAtom ){
         unsigned int index                 = (atomI == targetAtom) ? atomJ : atomI;
         unsigned int indexI                = (atomI == targetAtom) ? 0 : 1;
         unsigned int indexJ                = (atomI == targetAtom) ? 1 : 0;
         float forceSign                    = (atomI == targetAtom) ? 1.0f : -1.0f;

         debugArray[index].x                = (float) atomI;
         debugArray[index].y                = (float) atomJ;
         debugArray[index].z                = 3.0f;
         debugArray[index].w                = (float) y;
/*
         index                             += cAmoebaSim.paddedNumberOfAtoms;
         debugArray[index].x                = mask ? forceSign*force[0]     : 0.0f;
         debugArray[index].y                = mask ? forceSign*force[1]     : 0.0f;
         debugArray[index].z                = mask ? forceSign*force[2]     : 0.0f;


         index                             += cAmoebaSim.paddedNumberOfAtoms;
         debugArray[index].x                = mask ? torque[indexI][0] : 0.0f;
         debugArray[index].y                = mask ? torque[indexI][1] : 0.0f;
         debugArray[index].z                = mask ? torque[indexI][2] : 0.0f;

         index                             += cAmoebaSim.paddedNumberOfAtoms;
         debugArray[index].x                = mask ? torque[indexJ][0] : 0.0f;
         debugArray[index].y                = mask ? torque[indexJ][1] : 0.0f;
         debugArray[index].z                = mask ? torque[indexJ][2] : 0.0f;
*/
        index                             += cAmoebaSim.paddedNumberOfAtoms;
        int pullIndex                      = 0;
        debugArray[index].x                = pullBack[pullIndex].x;
        debugArray[index].y                = pullBack[pullIndex].y;
        debugArray[index].z                = pullBack[pullIndex].z;
        debugArray[index].w                = pullBack[pullIndex].w;

        index                             += cAmoebaSim.paddedNumberOfAtoms;
        pullIndex++;
        debugArray[index].x                = pullBack[pullIndex].x;
        debugArray[index].y                = pullBack[pullIndex].y;
        debugArray[index].z                = pullBack[pullIndex].z;
        debugArray[index].w                = pullBack[pullIndex].w;

        index                             += cAmoebaSim.paddedNumberOfAtoms;
        pullIndex++;
        debugArray[index].x                = pullBack[pullIndex].x;
        debugArray[index].y                = pullBack[pullIndex].y;
        debugArray[index].z                = pullBack[pullIndex].z;
        debugArray[index].w                = pullBack[pullIndex].w;
       
        index                             += cAmoebaSim.paddedNumberOfAtoms;
        pullIndex++;
        debugArray[index].x                = pullBack[pullIndex].x;
        debugArray[index].y                = pullBack[pullIndex].y;
        debugArray[index].z                = pullBack[pullIndex].z;
        debugArray[index].w                = pullBack[pullIndex].w;
       
#if 0
            index                             += cAmoebaSim.paddedNumberOfAtoms;
            debugArray[index].x                = scalingFactors[DScaleIndex];
            debugArray[index].y                = dScaleVal;
            debugArray[index].z                = scalingFactors[PScaleIndex];
            debugArray[index].w                = pScaleVal;

            index                             += cAmoebaSim.paddedNumberOfAtoms;
            debugArray[index].x                = scalingFactors[MScaleIndex];
            debugArray[index].y                = mScaleVal;
            for( int pIndex = 0; pIndex < 14; pIndex++ ){
                index                             += cAmoebaSim.paddedNumberOfAtoms;
                debugArray[index].x                = pullBack[pIndex].x;
                debugArray[index].y                = pullBack[pIndex].y;
                debugArray[index].z                = pullBack[pIndex].z;
                debugArray[index].w                = pullBack[pIndex].w;
            }

            index                             += cAmoebaSim.paddedNumberOfAtoms;
            debugArray[index].x                = labFrameDipole[3*atomI];
            debugArray[index].y                = labFrameDipole[3*atomI+1];
            debugArray[index].z                = labFrameDipole[3*atomI+2];
            debugArray[index].w                = 25.0f;

            index                             += cAmoebaSim.paddedNumberOfAtoms;
            debugArray[index].x                = labFrameDipole[3*atomJ];
            debugArray[index].y                = labFrameDipole[3*atomJ+1];
            debugArray[index].z                = labFrameDipole[3*atomJ+2];
            debugArray[index].w                = 26.0f;

            index                             += cAmoebaSim.paddedNumberOfAtoms;
            debugArray[index].x                = jDipole[0];
            debugArray[index].y                = jDipole[1];
            debugArray[index].z                = jDipole[2];
            debugArray[index].w                = 27.0f;
#endif


}
#endif


                    tj                  = (tj + 1) & (GRID - 1);
                }
            }

            // Write results

#ifdef USE_OUTPUT_BUFFER_PER_WARP

            float of;
            unsigned int offset                 = 3*(x + tgx + warp*cAmoebaSim.paddedNumberOfAtoms);
            of                                  = outputForce[offset];
            of                                 += localParticle.force[0];
            outputForce[offset]                 = of;

            of                                  = outputForce[offset+1];
            of                                 += localParticle.force[1];
            outputForce[offset+1]               = of;

            of                                  = outputForce[offset+2];
            of                                 += localParticle.force[2];
            outputForce[offset+2]               = of;

            of                                  = outputTorque[offset];
            of                                 += localParticle.torque[0];
            outputTorque[offset]                 = of;

            of                                  = outputTorque[offset+1];
            of                                 += localParticle.torque[1];
            outputTorque[offset+1]               = of;

            of                                  = outputTorque[offset+2];
            of                                 += localParticle.torque[2];
            outputTorque[offset+2]               = of;

            offset                              = 3*(y + tgx + warp*cAmoebaSim.paddedNumberOfAtoms);

            of                                  = outputForce[offset];
            of                                 += sA[threadIdx.x].force[0];
            outputForce[offset]                 = of;

            of                                  = outputForce[offset+1];
            of                                 += sA[threadIdx.x].force[1];
            outputForce[offset+1]               = of;

            of                                  = outputForce[offset+2];
            of                                 += sA[threadIdx.x].force[2];
            outputForce[offset+2]               = of;

            of                                  = outputTorque[offset];
            of                                 += sA[threadIdx.x].torque[0];
            outputTorque[offset]                = of;

            of                                  = outputTorque[offset+1];
            of                                 += sA[threadIdx.x].torque[1];
            outputTorque[offset+1]              = of;

            of                                  = outputTorque[offset+2];
            of                                 += sA[threadIdx.x].torque[2];
            outputTorque[offset+2]              = of;

#else
            unsigned int offset                 = 3*(x + tgx + (y >> GRIDBITS) * cAmoebaSim.paddedNumberOfAtoms);

            outputForce[offset]                 = localParticle.force[0];
            outputForce[offset+1]               = localParticle.force[1];
            outputForce[offset+2]               = localParticle.force[2];

            outputTorque[offset]                = localParticle.torque[0];
            outputTorque[offset+1]              = localParticle.torque[1];
            outputTorque[offset+2]              = localParticle.torque[2];

            offset                              = 3*(y + tgx + (x >> GRIDBITS) * cAmoebaSim.paddedNumberOfAtoms);

            outputForce[offset]                 = sA[threadIdx.x].force[0];
            outputForce[offset+1]               = sA[threadIdx.x].force[1];
            outputForce[offset+2]               = sA[threadIdx.x].force[2];

            outputTorque[offset]                = sA[threadIdx.x].torque[0];
            outputTorque[offset+1]              = sA[threadIdx.x].torque[1];
            outputTorque[offset+2]              = sA[threadIdx.x].torque[2];


#endif
            lasty = y;
        }

        pos++;
    }
    cSim.pEnergy[blockIdx.x * blockDim.x + threadIdx.x] += totalEnergy;
}
