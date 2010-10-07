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
void METHOD_NAME(kCalculateAmoebaPmeDirectElectrostatic, Forces_kernel)(
                            unsigned int* workUnit, float* outputForce, float* outputTorque

#ifdef AMOEBA_DEBUG
                           , float4* debugArray, unsigned int targetAtom
#endif
){

#ifdef AMOEBA_DEBUG
    int maxPullIndex = 2;
    float4 pullBack[12];
#endif

    extern __shared__ PmeDirectElectrostaticParticle sA[];

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
        int  dScaleMask;
        int2 pScaleMask;
        int2 mScaleMask;

        // Extract cell coordinates

        decodeCell( workUnit[pos], &x, &y, &bExclusionFlag );

        unsigned int tgx              = threadIdx.x & (GRID - 1);
        unsigned int tbx              = threadIdx.x - tgx;
        unsigned int tj               = tgx;

        PmeDirectElectrostaticParticle* psA    = &sA[tbx];
        unsigned int atomI            = x + tgx;
        PmeDirectElectrostaticParticle localParticle;
        loadPmeDirectElectrostaticShared(&localParticle, atomI );

        localParticle.force[0]        = 0.0f;
        localParticle.force[1]        = 0.0f;
        localParticle.force[2]        = 0.0f;

        localParticle.torque[0]       = 0.0f;
        localParticle.torque[1]       = 0.0f;
        localParticle.torque[2]       = 0.0f;

        scalingFactors[PScaleIndex]   = 1.0f;
        scalingFactors[DScaleIndex]   = 1.0f;
        scalingFactors[UScaleIndex]   = 1.0f;
        scalingFactors[MScaleIndex]   = 1.0f;

        if (x == y) // Handle diagonals uniquely at 50% efficiency
        {

            // load shared data

            loadPmeDirectElectrostaticShared( &(sA[threadIdx.x]), atomI );

            if (bExclusionFlag)
            {
                unsigned int xi       = x >> GRIDBITS;
                unsigned int cell     = xi + xi*cAmoebaSim.paddedNumberOfAtoms/GRID-xi*(xi+1)/2;
                dScaleMask            = cAmoebaSim.pD_ScaleIndices[cAmoebaSim.pScaleIndicesIndex[cell]+tgx];
                pScaleMask            = cAmoebaSim.pP_ScaleIndices[cAmoebaSim.pScaleIndicesIndex[cell]+tgx];
                mScaleMask            = cAmoebaSim.pM_ScaleIndices[cAmoebaSim.pScaleIndicesIndex[cell]+tgx];
            } else {
                scalingFactors[DScaleIndex] = scalingFactors[PScaleIndex] = scalingFactors[MScaleIndex] = 1.0f;
            }

            for (unsigned int j = 0; j < GRID; j++)
            {

                float force[3];
                float torque[2][3];

                unsigned int atomJ = y + j;

                // set scale factors

                if (bExclusionFlag)
                {
                    getMaskedDScaleFactor( j, dScaleMask, scalingFactors + DScaleIndex );
                    getMaskedPScaleFactor( j, pScaleMask, scalingFactors + PScaleIndex );
                    getMaskedMScaleFactor( j, mScaleMask, scalingFactors + MScaleIndex );
                }

                // force

                float energy;
                calculatePmeDirectElectrostaticPairIxn_kernel( localParticle,   psA[j],
                                                               scalingFactors, force, torque, &energy
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

                totalEnergy                       += mask ? 0.5*energy   : 0.0f;

#ifdef AMOEBA_DEBUG
/*
energy =  mask ? 0.5*energy   : 0.0f;
if( atomI < 200 && (fabs( energy ) > 1.0e+8 || energy != energy) ){
    debugSetup( atomI, atomJ, debugArray, pullBack );
} */

if( atomI == targetAtom ){
    unsigned int index                 = (atomI == targetAtom) ? atomJ : atomI;
    float blockId                      = 1.0f;

    debugArray[index].x                = (float) atomI;
    debugArray[index].y                = (float) atomJ;
    debugArray[index].z                = (float) y;
    debugArray[index].w                = blockId;

    index                             += cAmoebaSim.paddedNumberOfAtoms;
    debugArray[index].x                = mask ? force[0]     : 0.0f;
    debugArray[index].y                = mask ? force[1]     : 0.0f;
    debugArray[index].z                = mask ? force[2]     : 0.0f;
    debugArray[index].w                = blockId;


    index                             += cAmoebaSim.paddedNumberOfAtoms;
    debugArray[index].x                = mask ? torque[0][0] : 0.0f;
    debugArray[index].y                = mask ? torque[0][1] : 0.0f;
    debugArray[index].z                = mask ? torque[0][2] : 0.0f;
    debugArray[index].w                = mask ? energy       : 0.0f;

    index                             += cAmoebaSim.paddedNumberOfAtoms;
    debugArray[index].x                = mask ? torque[0][0] : 0.0f;
    debugArray[index].y                = mask ? torque[0][1] : 0.0f;
    debugArray[index].z                = mask ? torque[0][2] : 0.0f;
    debugArray[index].w                = (float) (blockIdx.x * blockDim.x + threadIdx.x);

    for( int pullIndex = 0; pullIndex < maxPullIndex; pullIndex++ ){
        index                             += cAmoebaSim.paddedNumberOfAtoms;
        debugArray[index].x                = pullBack[pullIndex].x;
        debugArray[index].y                = pullBack[pullIndex].y;
        debugArray[index].z                = pullBack[pullIndex].z;
        debugArray[index].w                = pullBack[pullIndex].w;
    }
}
#endif

            } // end of j-loop

            // include self energy and self torque

            if( atomI < cAmoebaSim.numberOfAtoms ){
                calculatePmeSelfTorqueElectrostaticPairIxn_kernel( localParticle );
                float energy;
                calculatePmeSelfEnergyElectrostaticPairIxn_kernel( localParticle, &energy );
                totalEnergy += energy;
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

        } else {

            if (lasty != y) {

                // load shared data

               loadPmeDirectElectrostaticShared( &(sA[threadIdx.x]), (y+tgx) );

            }

            unsigned int flags           = cSim.pInteractionFlag[pos];
            if (flags == 0) {
                // No interactions in this block.
            } else {

                if (lasty != y) {
    
                    // load shared data
    
                   loadPmeDirectElectrostaticShared( &(sA[threadIdx.x]), (y+tgx) );
    
                }
    
                sA[threadIdx.x].force[0]     = 0.0f;
                sA[threadIdx.x].force[1]     = 0.0f;
                sA[threadIdx.x].force[2]     = 0.0f;
    
                sA[threadIdx.x].torque[0]    = 0.0f;
                sA[threadIdx.x].torque[1]    = 0.0f;
                sA[threadIdx.x].torque[2]    = 0.0f;
    
                if( bExclusionFlag )
                {
                    unsigned int xi   = x >> GRIDBITS;
                    unsigned int yi   = y >> GRIDBITS;
                    unsigned int cell = xi+yi*cAmoebaSim.paddedNumberOfAtoms/GRID-yi*(yi+1)/2;
                    dScaleMask        = cAmoebaSim.pD_ScaleIndices[cAmoebaSim.pScaleIndicesIndex[cell]+tgx];
                    pScaleMask        = cAmoebaSim.pP_ScaleIndices[cAmoebaSim.pScaleIndicesIndex[cell]+tgx];
                    mScaleMask        = cAmoebaSim.pM_ScaleIndices[cAmoebaSim.pScaleIndicesIndex[cell]+tgx];
                } else {
                    scalingFactors[DScaleIndex] = scalingFactors[PScaleIndex] = scalingFactors[MScaleIndex] = 1.0f;
                }
       
                for (unsigned int j = 0; j < GRID; j++)
                {
           
                    unsigned int jIdx  = (flags == 0xFFFFFFFF) ? tj : j;
                    unsigned int atomJ = y + jIdx;
           
                    float force[3];
                    float torque[2][3];
           
                    // set scale factors
           
                    if( bExclusionFlag )
                    {
                        getMaskedDScaleFactor( jIdx, dScaleMask, scalingFactors + DScaleIndex );
                        getMaskedPScaleFactor( jIdx, pScaleMask, scalingFactors + PScaleIndex );
                        getMaskedMScaleFactor( jIdx, mScaleMask, scalingFactors + MScaleIndex );
                    }
           
                    // force
           
                    float energy;
                    calculatePmeDirectElectrostaticPairIxn_kernel( localParticle,   psA[jIdx],
                                                                   scalingFactors, force, torque, &energy
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
    
                    totalEnergy                    += mask ? energy        : 0.0f;

                    // add force and torque to atom J due atom I

                    if( flags == 0xFFFFFFFF ){

                        psA[jIdx].force[0]               -= mask ?  force[0]     : 0.0f;
                        psA[jIdx].force[1]               -= mask ?  force[1]     : 0.0f;
                        psA[jIdx].force[2]               -= mask ?  force[2]     : 0.0f;
    
                        psA[jIdx].torque[0]              += mask ?  torque[1][0] : 0.0f;
                        psA[jIdx].torque[1]              += mask ?  torque[1][1] : 0.0f;
                        psA[jIdx].torque[2]              += mask ?  torque[1][2] : 0.0f;

                    } else {

                        sA[threadIdx.x].tempForce[0]     = mask ? 0.0f : force[0];
                        sA[threadIdx.x].tempForce[1]     = mask ? 0.0f : force[1];
                        sA[threadIdx.x].tempForce[2]     = mask ? 0.0f : force[2];
   
                        sA[threadIdx.x].tempTorque[0]    = mask ? 0.0f : torque[1][0];
                        sA[threadIdx.x].tempTorque[1]    = mask ? 0.0f : torque[1][1];
                        sA[threadIdx.x].tempTorque[2]    = mask ? 0.0f : torque[1][2];

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
                            psA[jIdx].force[0]  -= sA[threadIdx.x].tempForce[0]  + sA[threadIdx.x+16].tempForce[0];
                            psA[jIdx].force[1]  -= sA[threadIdx.x].tempForce[1]  + sA[threadIdx.x+16].tempForce[1];
                            psA[jIdx].force[2]  -= sA[threadIdx.x].tempForce[2]  + sA[threadIdx.x+16].tempForce[2];

                            psA[jIdx].torque[0] += sA[threadIdx.x].tempTorque[0] + sA[threadIdx.x+16].tempTorque[0];
                            psA[jIdx].torque[1] += sA[threadIdx.x].tempTorque[1] + sA[threadIdx.x+16].tempTorque[1];
                            psA[jIdx].torque[2] += sA[threadIdx.x].tempTorque[2] + sA[threadIdx.x+16].tempTorque[2];
                        }
                    }
 
                    tj = (tj + 1) & (GRID - 1);

                } // end of j-loop

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

            } // end of pInteractionFlag block
        }
        pos++;
    }
    cSim.pEnergy[blockIdx.x * blockDim.x + threadIdx.x] += totalEnergy;
}
