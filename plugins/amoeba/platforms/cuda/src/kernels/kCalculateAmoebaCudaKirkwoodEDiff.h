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
__launch_bounds__(512, 1)
#elif (__CUDA_ARCH__ >= 120)
__launch_bounds__(96, 1)
#else
__launch_bounds__(32, 1)
#endif
void METHOD_NAME(kCalculateAmoebaCudaKirkwoodEDiff, Forces_kernel)(
                            unsigned int* workUnit, float* outputTorque
){

    extern __shared__ KirkwoodEDiffParticle sA[];

    unsigned int totalWarps      = gridDim.x*blockDim.x/GRID;
    unsigned int warp            = (blockIdx.x*blockDim.x+threadIdx.x)/GRID;
    unsigned int numWorkUnits    = cSim.pInteractionCount[0];
    unsigned int pos             = warp*numWorkUnits/totalWarps;
    unsigned int end             = (warp+1)*numWorkUnits/totalWarps;
    unsigned int lasty           = 0xFFFFFFFF;

    float totalEnergy            = 0.0f;
    float tinker_f               = (cAmoebaSim.electric/cAmoebaSim.dielec);

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

        KirkwoodEDiffParticle* psA    = &sA[tbx];

        unsigned int atomI            = x + tgx;
        KirkwoodEDiffParticle localParticle;
        loadKirkwoodEDiffShared(&localParticle, atomI );
        zeroKirkwoodEDiffParticleSharedField( &localParticle );

        if( x == y ){

            // load shared data

            loadKirkwoodEDiffShared( &(sA[threadIdx.x]), atomI );

            // first force and then torque
            
            unsigned int xi       = x >> GRIDBITS;
            unsigned int cell     = xi + xi*cSim.paddedNumberOfAtoms/GRID-xi*(xi+1)/2;
            int  dScaleMask       = cAmoebaSim.pD_ScaleIndices[cAmoebaSim.pScaleIndicesIndex[cell]+tgx];
            int2 pScaleMask       = cAmoebaSim.pP_ScaleIndices[cAmoebaSim.pScaleIndicesIndex[cell]+tgx];

            for (unsigned int j = 0; j < GRID; j++){

                unsigned int atomJ      = (y + j);

                float pScale;
                float dScale;
                getMaskedDScaleFactor( j, dScaleMask, &dScale );
                getMaskedPScaleFactor( j, pScaleMask, &pScale );

                if( (atomI != atomJ) && (atomI < cSim.atoms) && (atomJ < cSim.atoms) ){
                    float force[3];
                    float energy;
                    calculateKirkwoodEDiffPairIxnF1Scale_kernel( localParticle, psA[j], pScale, dScale, &energy, force);
                    totalEnergy            += 0.25f*energy;
                    localParticle.force[0] += force[0];
                    localParticle.force[1] += force[1];
                    localParticle.force[2] += force[2];
                }

            } // end of j-loop

            // scale and write results

            scale3dArray( tinker_f, localParticle.force );

#ifdef USE_OUTPUT_BUFFER_PER_WARP
            unsigned int offset                 = x + tgx + warp*cSim.paddedNumberOfAtoms;
            add3dArrayToFloat4(         offset, localParticle.force,  cSim.pForce4 );
#else
            unsigned int offset                 = x + tgx + (x >> GRIDBITS) * cSim.paddedNumberOfAtoms;
            add3dArrayToFloat4( offset, localParticle.force,  cSim.pForce4 );
#endif

            zeroKirkwoodEDiffParticleSharedField( &localParticle );
            for (unsigned int j = 0; j < GRID; j++)
            {

                unsigned int atomJ      = (y + j);

                float pScale;
                float dScale;
                getMaskedDScaleFactor( j, dScaleMask, &dScale );
                getMaskedPScaleFactor( j, pScaleMask, &pScale );

                if( (atomI != atomJ) && (atomI < cSim.atoms) && (atomJ < cSim.atoms) ){
                    float force[3];
                    calculateKirkwoodEDiffPairIxnT1Scale_kernel( localParticle, psA[j], pScale, dScale, force);
                    localParticle.force[0] += force[0];
                    localParticle.force[1] += force[1];
                    localParticle.force[2] += force[2];
                }

            } // end of j-loop

            // scale and write results

            scale3dArray( tinker_f, localParticle.force );

#ifdef USE_OUTPUT_BUFFER_PER_WARP
            offset                 = x + tgx + warp*cSim.paddedNumberOfAtoms;
            add3dArray(               3*offset, localParticle.force, outputTorque );
#else
            offset                 = x + tgx + (x >> GRIDBITS) * cSim.paddedNumberOfAtoms;
            add3dArray(       3*offset, localParticle.force, outputTorque );
#endif



        } else {

            // Read fixed atom data into registers and GRF

            if (lasty != y) {
                loadKirkwoodEDiffShared( &(sA[threadIdx.x]), (y+tgx) );
            }

            // zero j-atom output fields

            zeroKirkwoodEDiffParticleSharedField( &(sA[threadIdx.x]) );

            float dScale;
            float pScale;
            int  dScaleMask;
            int2 pScaleMask;
            if( bExclusionFlag ){
                unsigned int xi   = x >> GRIDBITS;
                unsigned int yi   = y >> GRIDBITS;
                unsigned int cell = xi+yi*cSim.paddedNumberOfAtoms/GRID-yi*(yi+1)/2;
                dScaleMask        = cAmoebaSim.pD_ScaleIndices[cAmoebaSim.pScaleIndicesIndex[cell]+tgx];
                pScaleMask        = cAmoebaSim.pP_ScaleIndices[cAmoebaSim.pScaleIndicesIndex[cell]+tgx];
            } else {
                 pScale = dScale = 1.0f;
            }

            for (unsigned int j = 0; j < GRID; j++) {

                unsigned int atomJ    = y + tj;
                if( (atomI < cSim.atoms) && ( atomJ < cSim.atoms) ){

                    float force[3];
                    float energy;
#ifdef APPLY_SCALE
                    if( bExclusionFlag ){
                        getMaskedDScaleFactor( tj, dScaleMask, &dScale );
                        getMaskedPScaleFactor( tj, pScaleMask, &pScale );
                        calculateKirkwoodEDiffPairIxnF1Scale_kernel( localParticle, psA[tj], pScale, dScale, &energy, force);
                    } else {
                        calculateKirkwoodEDiffPairIxnF1_kernel( localParticle, psA[tj], &energy, force);
                    }
#else
                    if( bExclusionFlag ){
                        getMaskedDScaleFactor( tj, dScaleMask, &dScale );
                        getMaskedPScaleFactor( tj, pScaleMask, &pScale );
                    }
                    calculateKirkwoodEDiffPairIxnF1Scale_kernel( localParticle, psA[tj], pScale, dScale, &energy, force);
#endif
           
                    totalEnergy            += 0.5f*energy;
                    localParticle.force[0] += force[0];
                    localParticle.force[1] += force[1];
                    localParticle.force[2] += force[2];
                    psA[tj].force[0]       -= force[0];
                    psA[tj].force[1]       -= force[1];
                    psA[tj].force[2]       -= force[2];
                }
                tj                  = (tj + 1) & (GRID - 1);

            } // end of j-loop

            // scale and write results

            scale3dArray( tinker_f, localParticle.force );
            scale3dArray( tinker_f, sA[threadIdx.x].force );

#ifdef USE_OUTPUT_BUFFER_PER_WARP

            unsigned int offset                 = x + tgx + warp*cSim.paddedNumberOfAtoms;
            add3dArrayToFloat4( offset, localParticle.force,  cSim.pForce4 );

            offset                              = y + tgx + warp*cSim.paddedNumberOfAtoms;
            add3dArrayToFloat4( offset, sA[threadIdx.x].force,  cSim.pForce4 );
#else
            unsigned int offset                 = x + tgx + (y >> GRIDBITS) * cSim.paddedNumberOfAtoms;
            add3dArrayToFloat4( offset, localParticle.force,  cSim.pForce4 );

            offset                              = y + tgx + (x >> GRIDBITS) * cSim.paddedNumberOfAtoms;
            add3dArrayToFloat4( offset, sA[threadIdx.x].force,  cSim.pForce4 );

#endif
            zeroKirkwoodEDiffParticleSharedField( &localParticle );
            zeroKirkwoodEDiffParticleSharedField( &(sA[threadIdx.x]) );
            for (unsigned int j = 0; j < GRID; j++) {

                unsigned int atomJ    = y + tj;
                if( (atomI < cSim.atoms) && ( atomJ < cSim.atoms) ){

                    float force[3];
#ifdef APPLY_SCALE
                    if( bExclusionFlag ){
                        getMaskedDScaleFactor( tj, dScaleMask, &dScale );
                        getMaskedPScaleFactor( tj, pScaleMask, &pScale );
                        calculateKirkwoodEDiffPairIxnT1Scale_kernel( localParticle, psA[tj], pScale, dScale, force);
                    } else {
                        calculateKirkwoodEDiffPairIxnT1_kernel( localParticle, psA[tj], force);
                    }
#else           
                    if( bExclusionFlag ){
                        getMaskedDScaleFactor( tj, dScaleMask, &dScale );
                        getMaskedPScaleFactor( tj, pScaleMask, &pScale );
                    }
                    calculateKirkwoodEDiffPairIxnT1Scale_kernel( localParticle, psA[tj], pScale, dScale, force);
#endif           
                    localParticle.force[0] += force[0];
                    localParticle.force[1] += force[1];
                    localParticle.force[2] += force[2];

#ifdef APPLY_SCALE
                    if( bExclusionFlag ){
                        calculateKirkwoodEDiffPairIxnT3Scale_kernel( localParticle, psA[tj], pScale, dScale, force);
                    } else {
                        calculateKirkwoodEDiffPairIxnT3_kernel( localParticle, psA[tj], force);
                    }
#else
                    calculateKirkwoodEDiffPairIxnT3Scale_kernel( localParticle, psA[tj], pScale, dScale, force);
#endif
                    psA[tj].force[0]       += force[0];
                    psA[tj].force[1]       += force[1];
                    psA[tj].force[2]       += force[2];
                }

                tj = (tj + 1) & (GRID - 1);

            } // end of j-loop

            // scale and write results

            scale3dArray( tinker_f, localParticle.force );
            scale3dArray( tinker_f, sA[threadIdx.x].force );

#ifdef USE_OUTPUT_BUFFER_PER_WARP

            offset                 = x + tgx + warp*cSim.paddedNumberOfAtoms;
            add3dArray(       3*offset, localParticle.force, outputTorque );

            offset                 = y + tgx + warp*cSim.paddedNumberOfAtoms;
            add3dArray(       3*offset, sA[threadIdx.x].force, outputTorque );
#else
            offset                 = x + tgx + (y >> GRIDBITS) * cSim.paddedNumberOfAtoms;
            add3dArray(       3*offset, localParticle.force, outputTorque );

            offset                 = y + tgx + (x >> GRIDBITS) * cSim.paddedNumberOfAtoms;
            add3dArray(       3*offset, sA[threadIdx.x].force, outputTorque );

#endif
            lasty = y;
        }

        pos++;
    }
    cSim.pEnergy[blockIdx.x*blockDim.x+threadIdx.x] += (tinker_f*totalEnergy);
}
