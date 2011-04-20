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
                            unsigned int* workUnit, float* outputTorque

#ifdef AMOEBA_DEBUG
                           , float4* debugArray, unsigned int targetAtom
#endif
){
#ifdef AMOEBA_DEBUG
    int maxPullIndex = 7;
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
    float4 forceTorqueEnergy[3];

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

        unsigned int tgx                       = threadIdx.x & (GRID - 1);
        unsigned int tbx                       = threadIdx.x - tgx;
        unsigned int tj                        = tgx;

        PmeDirectElectrostaticParticle* psA    = &sA[tbx];
        unsigned int atomI                     = x + tgx;
        PmeDirectElectrostaticParticle localParticle;
        loadPmeDirectElectrostaticShared(&localParticle, atomI );

        localParticle.force[0]                 = 0.0f;
        localParticle.force[1]                 = 0.0f;
        localParticle.force[2]                 = 0.0f;

        localParticle.torque[0]                = 0.0f;
        localParticle.torque[1]                = 0.0f;
        localParticle.torque[2]                = 0.0f;

        scalingFactors[UScaleIndex]            = 1.0f;

        if (x == y) // Handle diagonals uniquely at 50% efficiency
        {

            // load shared data

            loadPmeDirectElectrostaticShared( &(sA[threadIdx.x]), atomI );

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

            for (unsigned int j = 0; j < GRID; j++)
            {

                unsigned int atomJ = y + j;

                // set scale factors

                if (bExclusionFlag)
                {
                    getMaskedDScaleFactor( j, dScaleMask, scalingFactors + DScaleIndex );
                    getMaskedPScaleFactor( j, pScaleMask, scalingFactors + PScaleIndex );
                    getMaskedMScaleFactor( j, mScaleMask, scalingFactors + MScaleIndex );
                }

                // force

                calculatePmeDirectElectrostaticPairIxn_kernel( localParticle, psA[j], scalingFactors, forceTorqueEnergy
#ifdef AMOEBA_DEBUG
, pullBack
#endif
 );

                // nan*0.0 = nan not 0.0, so explicitly exclude (atomI == atomJ) contribution
                // by setting match flag

                if( (atomI != atomJ) && (atomI < cSim.atoms) && (atomJ < cSim.atoms) )
                {
                    localParticle.force[0]      += forceTorqueEnergy[0].x;
                    localParticle.force[1]      += forceTorqueEnergy[0].y;
                    localParticle.force[2]      += forceTorqueEnergy[0].z;
    
                    localParticle.torque[0]     += forceTorqueEnergy[1].x;
                    localParticle.torque[1]     += forceTorqueEnergy[1].y;
                    localParticle.torque[2]     += forceTorqueEnergy[1].z;
    
                    // energy for each diagonal-block ixn included twice, hence factor of 0.5

                    totalEnergy                 += 0.5*forceTorqueEnergy[0].w;
                }

#ifdef AMOEBA_DEBUG
if( atomI == targetAtom || atomJ == targetAtom ){

    unsigned int mask       =  ( (atomI == atomJ) || (atomI >= cSim.atoms) || (atomJ >= cSim.atoms) ) ? 0 : 1;
    unsigned int index                 = (atomI == targetAtom) ? atomJ : atomI;
    float blockId                      = 1.0f;

    debugArray[index].x                = (float) atomI;
    debugArray[index].y                = (float) atomJ;
    debugArray[index].z                = (float) y;
    debugArray[index].w                = blockId;

    index                             += cSim.paddedNumberOfAtoms;
    debugArray[index].x                = mask ? forceTorqueEnergy[0].x  : 0.0f;
    debugArray[index].y                = mask ? forceTorqueEnergy[0].y  : 0.0f;
    debugArray[index].z                = mask ? forceTorqueEnergy[0].z  : 0.0f;
    debugArray[index].w                = mask ? forceTorqueEnergy[0].w  : 0.0f;


    index                             += cSim.paddedNumberOfAtoms;
    debugArray[index].x                = mask ? forceTorqueEnergy[1].x : 0.0f;
    debugArray[index].y                = mask ? forceTorqueEnergy[1].y : 0.0f;
    debugArray[index].z                = mask ? forceTorqueEnergy[1].z : 0.0f;
    float offsetF                      = (float)(3*(x + tgx + (x >> GRIDBITS) * cSim.paddedNumberOfAtoms));
    debugArray[index].w                = offsetF;

    index                             += cSim.paddedNumberOfAtoms;
    debugArray[index].x                = mask ? forceTorqueEnergy[2].x : 0.0f;
    debugArray[index].y                = mask ? forceTorqueEnergy[2].y : 0.0f;
    debugArray[index].z                = mask ? forceTorqueEnergy[2].z : 0.0f;
    debugArray[index].w                = offsetF;

    for( int pullIndex = 0; pullIndex < maxPullIndex; pullIndex++ ){
        index                             += cSim.paddedNumberOfAtoms;
        debugArray[index].x                = pullBack[pullIndex].x;
        debugArray[index].y                = pullBack[pullIndex].y;
        debugArray[index].z                = pullBack[pullIndex].z;
        debugArray[index].w                = pullBack[pullIndex].w;
    }
}
#endif

            } // end of j-loop

            // include self energy and self torque

            if( atomI < cSim.atoms ){
                calculatePmeSelfTorqueElectrostaticPairIxn_kernel( localParticle );
                float energy;
                calculatePmeSelfEnergyElectrostaticPairIxn_kernel( localParticle, &energy );
                totalEnergy += energy;
            }

            localParticle.force[0]  *= conversionFactor;
            localParticle.force[1]  *= conversionFactor;
            localParticle.force[2]  *= conversionFactor;

            localParticle.torque[0] *= -conversionFactor;
            localParticle.torque[1] *= -conversionFactor;
            localParticle.torque[2] *= -conversionFactor;

            // Write results

#ifdef USE_OUTPUT_BUFFER_PER_WARP
            unsigned int offset                 = (x + tgx + warp*cSim.paddedNumberOfAtoms);
            add3dArrayToFloat4( offset, localParticle.force, cSim.pForce4 );
            add3dArray( 3*offset, localParticle.torque, outputTorque );
#else
            unsigned int offset                 = (x + tgx + (x >> GRIDBITS) * cSim.paddedNumberOfAtoms);
            add3dArrayToFloat4( offset, localParticle.force, cSim.pForce4 );
            load3dArray( 3*offset, localParticle.torque, outputTorque );
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

#ifdef CALCULATE_FULL_TILE
                flags = 0xFFFFFFFF;
#endif
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
                    unsigned int cell = xi+yi*cSim.paddedNumberOfAtoms/GRID-yi*(yi+1)/2;
                    dScaleMask        = cAmoebaSim.pD_ScaleIndices[cAmoebaSim.pScaleIndicesIndex[cell]+tgx];
                    pScaleMask        = cAmoebaSim.pP_ScaleIndices[cAmoebaSim.pScaleIndicesIndex[cell]+tgx];
                    mScaleMask        = cAmoebaSim.pM_ScaleIndices[cAmoebaSim.pScaleIndicesIndex[cell]+tgx];
                } else {
                    scalingFactors[DScaleIndex] = scalingFactors[PScaleIndex] = scalingFactors[MScaleIndex] = 1.0f;
                }
       
                for (unsigned int j = 0; j < GRID; j++)
                {
                    if( (flags & (1<<j) ) != 0)
                    {
                        unsigned int jIdx  = (flags == 0xFFFFFFFF) ? tj : j;
                        unsigned int atomJ = y + jIdx;

                        // set scale factors

                        if( bExclusionFlag )
                        {
                            getMaskedDScaleFactor( jIdx, dScaleMask, scalingFactors + DScaleIndex );
                            getMaskedPScaleFactor( jIdx, pScaleMask, scalingFactors + PScaleIndex );
                            getMaskedMScaleFactor( jIdx, mScaleMask, scalingFactors + MScaleIndex );
                        }

                        // force

                        calculatePmeDirectElectrostaticPairIxn_kernel( localParticle, psA[jIdx],
                                                                       scalingFactors, forceTorqueEnergy
#ifdef AMOEBA_DEBUG
    , pullBack
#endif
         );

                        // check if atoms out-of-bounds

                        if( (atomI < cSim.atoms) && (atomJ < cSim.atoms) )
                        {
                            // add force and torque to atom I due atom J
    
                            localParticle.force[0]         += forceTorqueEnergy[0].x;
                            localParticle.force[1]         += forceTorqueEnergy[0].y;
                            localParticle.force[2]         += forceTorqueEnergy[0].z;
    
                            totalEnergy                    += forceTorqueEnergy[0].w;
    
                            localParticle.torque[0]        += forceTorqueEnergy[1].x;
                            localParticle.torque[1]        += forceTorqueEnergy[1].y;
                            localParticle.torque[2]        += forceTorqueEnergy[1].z;
    
                            // add force and torque to atom J due atom I
    
                            if( flags == 0xFFFFFFFF ){
    
                                psA[jIdx].force[0]         -= forceTorqueEnergy[0].x;
                                psA[jIdx].force[1]         -= forceTorqueEnergy[0].y;
                                psA[jIdx].force[2]         -= forceTorqueEnergy[0].z;
    
                                psA[jIdx].torque[0]        += forceTorqueEnergy[2].x;
                                psA[jIdx].torque[1]        += forceTorqueEnergy[2].y;
                                psA[jIdx].torque[2]        += forceTorqueEnergy[2].z;

#ifndef CALCULATE_FULL_TILE
                            } else {
    
                                sA[threadIdx.x].tempForce[0]  = forceTorqueEnergy[0].x;
                                sA[threadIdx.x].tempForce[1]  = forceTorqueEnergy[0].y;
                                sA[threadIdx.x].tempForce[2]  = forceTorqueEnergy[0].z;
    
                                sA[threadIdx.x].tempTorque[0] = forceTorqueEnergy[2].x;
                                sA[threadIdx.x].tempTorque[1] = forceTorqueEnergy[2].y;
                                sA[threadIdx.x].tempTorque[2] = forceTorqueEnergy[2].z;
    
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
#endif
                            }
                        } // end of atoms out-of-bounds
                    } // end of flags&(1<<j block
 
#ifdef AMOEBA_DEBUG
unsigned int jIdx  = (flags == 0xFFFFFFFF) ? tj : j;
unsigned int atomJ = y + jIdx;
unsigned int mask  =  ( (atomI >= cSim.atoms) || (atomJ >= cSim.atoms) ) ? 0 : 1;
if( atomI == targetAtom || atomJ == targetAtom ){
    unsigned int index                 = (atomI == targetAtom) ? atomJ : atomI;

    debugArray[index].x                = (float) atomI;
    debugArray[index].y                = (float) atomJ;
    debugArray[index].z                = (float) y;
    debugArray[index].w                = (flags == 0xFFFFFFFF) ? (float) -141.0f : -151.0f;

    index                             += cSim.paddedNumberOfAtoms;
    debugArray[index].x                = mask ? forceTorqueEnergy[0].x  : 0.0f;
    debugArray[index].y                = mask ? forceTorqueEnergy[0].y  : 0.0f;
    debugArray[index].z                = mask ? forceTorqueEnergy[0].z  : 0.0f;
    debugArray[index].w                = mask ? forceTorqueEnergy[0].w  : 0.0f;


    index                             += cSim.paddedNumberOfAtoms;
    debugArray[index].x                = mask ? forceTorqueEnergy[1].x : 0.0f;
    debugArray[index].y                = mask ? forceTorqueEnergy[1].y : 0.0f;
    debugArray[index].z                = mask ? forceTorqueEnergy[1].z : 0.0f;
    float offsetF                      = (float)(3*(y + tgx + (x >> GRIDBITS) * cSim.paddedNumberOfAtoms));
    debugArray[index].w                = offsetF;

    index                             += cSim.paddedNumberOfAtoms;
    debugArray[index].x                = mask ? forceTorqueEnergy[2].x : 0.0f;
    debugArray[index].y                = mask ? forceTorqueEnergy[2].y : 0.0f;
    debugArray[index].z                = mask ? forceTorqueEnergy[2].z : 0.0f;
    offsetF                            = (float) (3*(x + tgx + (y >> GRIDBITS) * cSim.paddedNumberOfAtoms));
    debugArray[index].w                = offsetF;

    for( int pullIndex = 0; pullIndex < maxPullIndex; pullIndex++ ){
        index                             += cSim.paddedNumberOfAtoms;
        debugArray[index].x                = pullBack[pullIndex].x;
        debugArray[index].y                = pullBack[pullIndex].y;
        debugArray[index].z                = pullBack[pullIndex].z;
        debugArray[index].w                = pullBack[pullIndex].w;
    }
}
#endif
                    tj = (tj + 1) & (GRID - 1);

                } // end of j-loop

                localParticle.force[0]    *=  conversionFactor;
                localParticle.force[1]    *=  conversionFactor;
                localParticle.force[2]    *=  conversionFactor;
    
                localParticle.torque[0]   *= -conversionFactor;
                localParticle.torque[1]   *= -conversionFactor;
                localParticle.torque[2]   *= -conversionFactor;
    
                sA[threadIdx.x].force[0]  *=  conversionFactor;
                sA[threadIdx.x].force[1]  *=  conversionFactor;
                sA[threadIdx.x].force[2]  *=  conversionFactor;
    
                sA[threadIdx.x].torque[0] *= -conversionFactor;
                sA[threadIdx.x].torque[1] *= -conversionFactor;
                sA[threadIdx.x].torque[2] *= -conversionFactor;
    
                // Write results
    
#ifdef USE_OUTPUT_BUFFER_PER_WARP
    
                unsigned int offset                 = (x + tgx + warp*cSim.paddedNumberOfAtoms);
                add3dArrayToFloat4( offset, localParticle.force,  cSim.pForce4 );
                add3dArray(       3*offset, localParticle.torque, outputTorque );
    
                offset                              = (y + tgx + warp*cSim.paddedNumberOfAtoms);
                add3dArrayToFloat4( offset, sA[threadIdx.x].force,  cSim.pForce4 );
                add3dArray(       3*offset, sA[threadIdx.x].torque, outputTorque );
    
#else
                unsigned int offset                 = (x + tgx + (y >> GRIDBITS) * cSim.paddedNumberOfAtoms);
                add3dArrayToFloat4( offset, localParticle.force,  cSim.pForce4 );
                load3dArray(       3*offset, localParticle.torque, outputTorque );
    
                offset                              = (y + tgx + (x >> GRIDBITS) * cSim.paddedNumberOfAtoms);
    
                add3dArrayToFloat4( offset, sA[threadIdx.x].force,  cSim.pForce4 );
                load3dArray(       3*offset, sA[threadIdx.x].torque, outputTorque );
    
#endif
                lasty = y;

            } // end of pInteractionFlag block
        }
        pos++;
    }
    cSim.pEnergy[blockIdx.x * blockDim.x + threadIdx.x] -= conversionFactor*totalEnergy;
}
