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
__launch_bounds__(96, 1)
#else
__launch_bounds__(32, 1)
#endif
void METHOD_NAME(kCalculateAmoebaCudaKirkwoodEDiff, Forces_kernel)(
                            unsigned int* workUnit,
                            float4* atomCoord,
                            float* labFrameDipole,
                            float* labFrameQuadrupole,
                            float* inducedDipole,
                            float* inducedDipolePolar,
                            float* inducedDipoleS,
                            float* inducedDipolePolarS,
                            float* outputTorque
#ifdef AMOEBA_DEBUG
                           , float4* debugArray, unsigned int targetAtom
#endif
){

#ifdef AMOEBA_DEBUG
    float4 pullBack[20];
#endif

    extern __shared__ KirkwoodEDiffParticle sA[];

    unsigned int totalWarps      = gridDim.x*blockDim.x/GRID;
    unsigned int warp            = (blockIdx.x*blockDim.x+threadIdx.x)/GRID;
    unsigned int numWorkUnits    = cSim.pInteractionCount[0];
    unsigned int pos             = warp*numWorkUnits/totalWarps;
    unsigned int end             = (warp+1)*numWorkUnits/totalWarps;
    unsigned int lasty           = 0xFFFFFFFF;

    float totalEnergy     = 0.0f;
    float tinker_f        = (cAmoebaSim.electric/cAmoebaSim.dielec);

    while (pos < end)
    {

        unsigned int x;
        unsigned int y;
        bool bExclusionFlag;

        float force[3];
        float torqueI[3];
        float torqueJ[3];
        float energy;

        // Extract cell coordinates

        decodeCell( workUnit[pos], &x, &y, &bExclusionFlag );

        unsigned int tgx              = threadIdx.x & (GRID - 1);
        unsigned int tbx              = threadIdx.x - tgx;
        unsigned int tj               = tgx;

        KirkwoodEDiffParticle* psA    = &sA[tbx];

        unsigned int atomI            = x + tgx;
        KirkwoodEDiffParticle localParticle;
        loadKirkwoodEDiffShared(&localParticle, atomI,
                                atomCoord,
                                labFrameDipole, labFrameQuadrupole,
                                inducedDipole,  inducedDipolePolar,
                                inducedDipoleS, inducedDipolePolarS );

        localParticle.force[0]                   = 0.0f;
        localParticle.force[1]                   = 0.0f;
        localParticle.force[2]                   = 0.0f;

        localParticle.torque[0]                  = 0.0f;
        localParticle.torque[1]                  = 0.0f;
        localParticle.torque[2]                  = 0.0f;

        if (x == y) // Handle diagonals uniquely at 50% efficiency
        {

            // load shared data

             loadKirkwoodEDiffShared( &(sA[threadIdx.x]), atomI,
                                      atomCoord,
                                      labFrameDipole, labFrameQuadrupole,
                                      inducedDipole,  inducedDipolePolar,
                                      inducedDipoleS, inducedDipolePolarS );

            if (!bExclusionFlag)
            {

                float pScale = 1.0f;
                float dScale = 1.0f;

                // this branch is never exercised since it includes the
                // interaction between atomI and itself which is always excluded

                for (unsigned int j = 0; j < GRID; j++)
                {

                    unsigned int atomJ      = (y + j);

                    
                    calculateKirkwoodEDiffPairIxn_kernel( localParticle,                              psA[j],
                                                          pScale,                                     dScale,
                                                          &energy, force, torqueI, torqueJ
#ifdef AMOEBA_DEBUG
                                              , pullBack
#endif
                                            );

                    unsigned int mask       =  ( (atomI >= cSim.atoms) || (atomJ >= cSim.atoms) ) ? 0 : 1;

                    // torques include i == j contribution

                    localParticle.torque[0]           += mask ? torqueI[0]  : 0.0f;
                    localParticle.torque[1]           += mask ? torqueI[1]  : 0.0f;
                    localParticle.torque[2]           += mask ? torqueI[2]  : 0.0f;
 
                    totalEnergy            += mask ? 0.5f*energy : 0.0f;

                    // add to field at atomI the field due atomJ's charge/dipole/quadrupole

                    mask                    =  (atomI == atomJ) ? 0 : mask;

                    localParticle.force[0]            += mask ? force[0]    : 0.0f;
                    localParticle.force[1]            += mask ? force[1]    : 0.0f;
                    localParticle.force[2]            += mask ? force[2]    : 0.0f;


#ifdef AMOEBA_DEBUG
if( atomI == targetAtom  || atomJ == targetAtom ){

            unsigned int index                 = (atomI == targetAtom) ? atomJ : atomI;
            float* torqueIPtr                  = (atomI == targetAtom) ? torqueI : torqueJ;
            float* torqueJPtr                  = (atomI == targetAtom) ? torqueJ : torqueI;

            debugArray[index].x                = (float) atomI;
            debugArray[index].y                = (float) atomJ;
            mask                               =  ( (atomI >= cSim.atoms) || (atomJ >= cSim.atoms) ) ? 0 : 1;
            debugArray[index].z                = mask ? tinker_f*energy : 0.0f;

            index = debugAccumulate( index, debugArray, force,          mask, 1.0f );

            mask  =  ( (atomI >= cSim.atoms) || (atomJ >= cSim.atoms) ) ? 0 : 1;
            index = debugAccumulate( index, debugArray, torqueIPtr, mask, 2.0f );
            index = debugAccumulate( index, debugArray, torqueJPtr, mask, 3.0f );
}
#endif
                }
            }
            else // bExclusion
            {
                unsigned int xi       = x >> GRIDBITS;
                unsigned int cell     = xi + xi*cSim.paddedNumberOfAtoms/GRID-xi*(xi+1)/2;
                int  dScaleMask       = cAmoebaSim.pD_ScaleIndices[cAmoebaSim.pScaleIndicesIndex[cell]+tgx];
                int2 pScaleMask       = cAmoebaSim.pP_ScaleIndices[cAmoebaSim.pScaleIndicesIndex[cell]+tgx];

                for (unsigned int j = 0; j < GRID; j++)
                {

                    unsigned int atomJ      = (y + j);

                    float pScale;
                    float dScale;
                    getMaskedDScaleFactor( j, dScaleMask, &dScale );
                    getMaskedPScaleFactor( j, pScaleMask, &pScale );

                    calculateKirkwoodEDiffPairIxn_kernel( localParticle,                              psA[j],
                                                          pScale,                                     dScale,
                                                          &energy, force, torqueI, torqueJ
#ifdef AMOEBA_DEBUG
                                              , pullBack
#endif
                                            );

                    unsigned int mask       =  ( (atomI == atomJ) || (atomI >= cSim.atoms) || (atomJ >= cSim.atoms) ) ? 0 : 1;

                    // torques include i == j contribution

                    localParticle.torque[0]           += mask ? torqueI[0]  : 0.0f;
                    localParticle.torque[1]           += mask ? torqueI[1]  : 0.0f;
                    localParticle.torque[2]           += mask ? torqueI[2]  : 0.0f;

                    totalEnergy            += mask ? 0.5f*energy : 0.0f;

                    // add to field at atomI the field due atomJ's charge/dipole/quadrupole

                    localParticle.force[0]            += mask ? force[0]    : 0.0f;
                    localParticle.force[1]            += mask ? force[1]    : 0.0f;
                    localParticle.force[2]            += mask ? force[2]    : 0.0f;


#ifdef AMOEBA_DEBUG
if( atomI == targetAtom  || atomJ == targetAtom ){

            unsigned int index                 = (atomI == targetAtom) ? atomJ : atomI;
            float* torqueIPtr                  = (atomI == targetAtom) ? torqueI : torqueJ;
            float* torqueJPtr                  = (atomI == targetAtom) ? torqueJ : torqueI;

            debugArray[index].x                = (float) atomI;
            debugArray[index].y                = (float) atomJ;
            debugArray[index].z                = mask ? tinker_f*energy : 0.0f;

            index = debugAccumulate( index, debugArray, force,          mask, 1.0f );

            //mask  =  ( (atomI >= cSim.atoms) || (atomJ >= cSim.atoms) ) ? 0 : 1;
            index = debugAccumulate( index, debugArray, torqueIPtr, mask, 2.0f );
            index = debugAccumulate( index, debugArray, torqueJPtr, mask, 3.0f );
}
#endif

                } // end of j-loop

            } // end of exclusion

            // scale and write results

            scale3dArray( tinker_f, localParticle.force );
            scale3dArray( tinker_f, localParticle.torque );

#ifdef USE_OUTPUT_BUFFER_PER_WARP
            unsigned int offset                 = x + tgx + warp*cSim.paddedNumberOfAtoms;

            add3dArrayToFloat4(         offset, localParticle.force,  cSim.pForce4 );
            load3dArrayBufferPerWarp( 3*offset, localParticle.torque, outputTorque );

#else
            unsigned int offset                 = x + tgx + (x >> GRIDBITS) * cSim.paddedNumberOfAtoms;
            add3dArrayToFloat4( offset, localParticle.force,  cSim.pForce4 );
            load3dArray(      3*offset, localParticle.torque, outputTorque );

#endif


        }
        else        // 100% utilization
        {
            // Read fixed atom data into registers and GRF

            if (lasty != y)
            {
                // load shared data

                loadKirkwoodEDiffShared( &(sA[threadIdx.x]), (y+tgx),
                                         atomCoord, labFrameDipole, labFrameQuadrupole,
                                         inducedDipole,  inducedDipolePolar,
                                         inducedDipoleS, inducedDipolePolarS );


            }

            // zero j-atom output fields

            zeroKirkwoodEDiffParticleSharedField( &(sA[threadIdx.x]) );

            if (!bExclusionFlag)
            {

                float pScale = 1.0f;
                float dScale = 1.0f;

                for (unsigned int j = 0; j < GRID; j++)
                {

                    unsigned int atomJ    = y + tj;

                    calculateKirkwoodEDiffPairIxn_kernel( localParticle,                              psA[tj],
                                                          pScale,                                     dScale,
                                                          &energy,                                    force, 
                                                          torqueI,                                    torqueJ
#ifdef AMOEBA_DEBUG
                                              ,pullBack
#endif
                                            );


                    unsigned int mask    =  ( (atomI >= cSim.atoms) || ( atomJ >= cSim.atoms) ) ? 0 : 1;

                    // add force and torque to atom I due atom J

                    localParticle.force[0]               += mask ? force[0]    : 0.0f;
                    localParticle.force[1]               += mask ? force[1]    : 0.0f;
                    localParticle.force[2]               += mask ? force[2]    : 0.0f;

                    localParticle.torque[0]              += mask ? torqueI[0]  : 0.0f;
                    localParticle.torque[1]              += mask ? torqueI[1]  : 0.0f;
                    localParticle.torque[2]              += mask ? torqueI[2]  : 0.0f;

                    totalEnergy               += mask ? energy      : 0.0f;

                    // add force and torque to atom J due atom I

                    psA[tj].force[0]          -= mask ? force[0]    : 0.0f;
                    psA[tj].force[1]          -= mask ? force[1]    : 0.0f;
                    psA[tj].force[2]          -= mask ? force[2]    : 0.0f;

                    psA[tj].torque[0]         += mask ? torqueJ[0]  : 0.0f;
                    psA[tj].torque[1]         += mask ? torqueJ[1]  : 0.0f;
                    psA[tj].torque[2]         += mask ? torqueJ[2]  : 0.0f;

#ifdef AMOEBA_DEBUG
if( atomI == targetAtom  || atomJ == targetAtom ){
            unsigned int index                 = (atomI == targetAtom) ? atomJ : atomI;
            float* torqueIPtr                  = (atomI == targetAtom) ? torqueI : torqueJ;
            float* torqueJPtr                  = (atomI == targetAtom) ? torqueJ : torqueI;

            debugArray[index].x                = (float) atomI;
            debugArray[index].y                = (float) atomJ;
            debugArray[index].z                = mask ? tinker_f*energy : 0.0f;

            index = debugAccumulate( index, debugArray, force,      mask, -1.0f );
            index = debugAccumulate( index, debugArray, torqueIPtr, mask, -2.0f );
            index = debugAccumulate( index, debugArray, torqueJPtr, mask, -3.0f );

}
#endif

                    tj                  = (tj + 1) & (GRID - 1);

                } // end of j-loop

            }
            else  // exclusion
            {

                unsigned int xi   = x >> GRIDBITS;
                unsigned int yi   = y >> GRIDBITS;
                unsigned int cell = xi+yi*cSim.paddedNumberOfAtoms/GRID-yi*(yi+1)/2;
                int  dScaleMask   = cAmoebaSim.pD_ScaleIndices[cAmoebaSim.pScaleIndicesIndex[cell]+tgx];
                int2 pScaleMask   = cAmoebaSim.pP_ScaleIndices[cAmoebaSim.pScaleIndicesIndex[cell]+tgx];

                for (unsigned int j = 0; j < GRID; j++)
                {

                    unsigned int atomJ    = y + tj;
                    float dScale;
                    float pScale;
                    getMaskedDScaleFactor( tj, dScaleMask, &dScale );
                    getMaskedPScaleFactor( tj, pScaleMask, &pScale );

                    calculateKirkwoodEDiffPairIxn_kernel( localParticle,                              psA[tj],
                                                          pScale,                                     dScale,
                                                          &energy,                                    force, 
                                                          torqueI,                                    torqueJ
#ifdef AMOEBA_DEBUG
                                              ,pullBack
#endif
                                            );


                    unsigned int mask    =  ( (atomI >= cSim.atoms) || ( atomJ >= cSim.atoms) ) ? 0 : 1;

                    // add force and torque to atom I due atom J

                    localParticle.force[0]               += mask ? force[0]    : 0.0f;
                    localParticle.force[1]               += mask ? force[1]    : 0.0f;
                    localParticle.force[2]               += mask ? force[2]    : 0.0f;

                    localParticle.torque[0]              += mask ? torqueI[0]  : 0.0f;
                    localParticle.torque[1]              += mask ? torqueI[1]  : 0.0f;
                    localParticle.torque[2]              += mask ? torqueI[2]  : 0.0f;

                    totalEnergy               += mask ? energy      : 0.0f;

                    // add force and torque to atom J due atom I

                    psA[tj].force[0]          -= mask ? force[0]    : 0.0f;
                    psA[tj].force[1]          -= mask ? force[1]    : 0.0f;
                    psA[tj].force[2]          -= mask ? force[2]    : 0.0f;

                    psA[tj].torque[0]         += mask ? torqueJ[0]  : 0.0f;
                    psA[tj].torque[1]         += mask ? torqueJ[1]  : 0.0f;
                    psA[tj].torque[2]         += mask ? torqueJ[2]  : 0.0f;

#ifdef AMOEBA_DEBUG
if( atomI == targetAtom  || atomJ == targetAtom ){
            unsigned int index                 = (atomI == targetAtom) ? atomJ : atomI;
            float* torqueIPtr                  = (atomI == targetAtom) ? torqueI : torqueJ;
            float* torqueJPtr                  = (atomI == targetAtom) ? torqueJ : torqueI;

            debugArray[index].x                = (float) atomI;
            debugArray[index].y                = (float) atomJ;
            debugArray[index].z                = mask ? tinker_f*energy : 0.0f;

            index = debugAccumulate( index, debugArray, force,      mask, -1.0f );
            index = debugAccumulate( index, debugArray, torqueIPtr, mask, -2.0f );
            index = debugAccumulate( index, debugArray, torqueJPtr, mask, -3.0f );

}
#endif

                    tj                  = (tj + 1) & (GRID - 1);

                } // end of j-loop

            } // end of exclusion

            // scale and write results

            scale3dArray( tinker_f, localParticle.force );
            scale3dArray( tinker_f, localParticle.torque );

            scale3dArray( tinker_f, sA[threadIdx.x].force );
            scale3dArray( tinker_f, sA[threadIdx.x].torque );

#ifdef USE_OUTPUT_BUFFER_PER_WARP

            unsigned int offset                 = x + tgx + warp*cSim.paddedNumberOfAtoms;

            add3dArrayToFloat4( offset, localParticle.force,  cSim.pForce4 );
            load3dArrayBufferPerWarp( 3*offset, localParticle.torque, outputTorque );

            offset                              = y + tgx + warp*cSim.paddedNumberOfAtoms;

            add3dArrayToFloat4(         offset, sA[threadIdx.x].force,  cSim.pForce4 );
            load3dArrayBufferPerWarp( 3*offset, sA[threadIdx.x].torque, outputTorque );
#else
            unsigned int offset                 = x + tgx + (y >> GRIDBITS) * cSim.paddedNumberOfAtoms;

            add3dArrayToFloat4( offset, localParticle.force,  cSim.pForce4 );
            load3dArray(      3*offset, localParticle.torque, outputTorque );


            offset                              = y + tgx + (x >> GRIDBITS) * cSim.paddedNumberOfAtoms;

            add3dArrayToFloat4( offset, sA[threadIdx.x].force,  cSim.pForce4 );
            load3dArray(      3*offset, sA[threadIdx.x].torque, outputTorque );

#endif
            lasty = y;
        }

        pos++;
    }
    cSim.pEnergy[blockIdx.x*blockDim.x+threadIdx.x] += (tinker_f*totalEnergy);
}
