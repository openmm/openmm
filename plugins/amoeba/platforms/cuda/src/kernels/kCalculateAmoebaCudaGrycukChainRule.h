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
#elif (__CUDA_ARCH__ >= 120)
__launch_bounds__(GT2XX_NONBOND_THREADS_PER_BLOCK, 1)
#else
__launch_bounds__(G8X_NONBOND_THREADS_PER_BLOCK, 1)
#endif
void METHOD_NAME(kCalculateAmoebaGrycukChainRule, _kernel)( unsigned int* workUnit 
#ifdef AMOEBA_DEBUG
                           , float4* debugArray, unsigned int targetAtom
#endif
){

    extern __shared__ GrycukChainRuleParticle sAChainRule[];

    unsigned int totalWarps      = gridDim.x*blockDim.x/GRID;
    unsigned int warp            = (blockIdx.x*blockDim.x+threadIdx.x)/GRID;
    unsigned int numWorkUnits    = cSim.pInteractionCount[0];
    unsigned int pos             = warp*numWorkUnits/totalWarps;
    unsigned int end             = (warp+1)*numWorkUnits/totalWarps;
    unsigned int lasty           = 0xFFFFFFFF;

#ifdef AMOEBA_DEBUG
    float4 pullDebug[5];
#endif

    while (pos < end)
    {

        unsigned int x;
        unsigned int y;
        bool bExclusionFlag;

        // Extract cell coordinates

        decodeCell( workUnit[pos], &x, &y, &bExclusionFlag );

        unsigned int tgx                          = threadIdx.x & (GRID - 1);
        unsigned int tbx                          = threadIdx.x - tgx;
        unsigned int tj                           = tgx;

        GrycukChainRuleParticle*  psAChainRule    = &sAChainRule[tbx];
        unsigned int atomI                        = x + tgx;
        GrycukChainRuleParticle localParticle;
        loadGrycukChainRuleParticleShared( &localParticle, atomI );

        zeroGrycukChainRuleParticleSharedField( &localParticle );

        if (x == y){

            // load shared data and zero force

            loadGrycukChainRuleParticleShared( &(sAChainRule[threadIdx.x]), atomI );
            zeroGrycukChainRuleParticleSharedField( &(sAChainRule[threadIdx.x]));

            for (unsigned int j = (tgx+1)&(GRID-1); j != tgx; j = (j+1)&(GRID-1))
            {
                float localForce[3];
                calculateGrycukChainRulePairIxn_kernel( localParticle, psAChainRule[j], localForce
#ifdef AMOEBA_DEBUG
,  pullDebug
#endif
 );
                if( (atomI != (y + j)) && (atomI < cSim.atoms) && ((y+j) < cSim.atoms) ){

                    localParticle.force[0]     -= localForce[0];
                    localParticle.force[1]     -= localForce[1];
                    localParticle.force[2]     -= localForce[2];

                    psAChainRule[j].force[0]   += localForce[0];
                    psAChainRule[j].force[1]   += localForce[1];
                    psAChainRule[j].force[2]   += localForce[2];

#ifdef AMOEBA_DEBUG
if( atomI == targetAtom || (y+j) == targetAtom ){
        unsigned int index                 = (atomI == targetAtom) ? (y + j) : atomI;

        debugArray[index].x                = (float) atomI;
        debugArray[index].y                = (float) (y + j); 
        debugArray[index].z                = -1.0f;
        debugArray[index].w                = -1.0f;

        index                             += cSim.paddedNumberOfAtoms;
        debugArray[index].x                = (float) x;
        debugArray[index].y                = (float) y;
        debugArray[index].z                = (float) tgx;
        debugArray[index].w                = -2.0f;

        index                             += cSim.paddedNumberOfAtoms;
        debugArray[index].x                = pullDebug[0].x;
        debugArray[index].y                = pullDebug[0].y;
        debugArray[index].z                = pullDebug[0].z;
        debugArray[index].w                = pullDebug[0].w;

        index                             += cSim.paddedNumberOfAtoms;
        debugArray[index].x                = pullDebug[1].x;
        debugArray[index].y                = pullDebug[1].y;
        debugArray[index].z                = pullDebug[1].z;
        debugArray[index].w                = pullDebug[1].w;

        index                             += cSim.paddedNumberOfAtoms;
        debugArray[index].x                = localForce[0];
        debugArray[index].y                = localForce[1];
        debugArray[index].z                = localForce[2];
        debugArray[index].w                = -12.0f;

 calculateGrycukChainRulePairIxn_kernel( psAChainRule[j], localParticle, localForce ,  pullDebug );

        index                             += cSim.paddedNumberOfAtoms;
        debugArray[index].x                = pullDebug[0].x;
        debugArray[index].y                = pullDebug[0].y;
        debugArray[index].z                = pullDebug[0].z;
        debugArray[index].w                = -13.0f;

        index                             += cSim.paddedNumberOfAtoms;
        debugArray[index].x                = localForce[0];
        debugArray[index].y                = localForce[1];
        debugArray[index].z                = localForce[2];
        debugArray[index].w                = -14.0f;
}
#endif

                }
            }

            // Write results
            float4 of; 
#ifdef USE_OUTPUT_BUFFER_PER_WARP
            unsigned int offset         = x + tgx + warp*cSim.stride;
#else
            unsigned int offset         = x + tgx + (x >> GRIDBITS) * cSim.stride;
#endif
            of                          = cSim.pForce4[offset];
            of.x                       += localParticle.force[0]  + sAChainRule[threadIdx.x].force[0];
            of.y                       += localParticle.force[1]  + sAChainRule[threadIdx.x].force[1];
            of.z                       += localParticle.force[2]  + sAChainRule[threadIdx.x].force[2];
            cSim.pForce4[offset]       = of; 

        } else {

            if (lasty != y) {
                unsigned int atomJ        = y + tgx;
                loadGrycukChainRuleParticleShared( &(sAChainRule[threadIdx.x]), atomJ );
            }

           // zero shared fields

            zeroGrycukChainRuleParticleSharedField(  &(sAChainRule[threadIdx.x]) );

            for (unsigned int j = 0; j < GRID; j++)
            {

                if( (atomI < cSim.atoms) && ((y+tj) < cSim.atoms) ){
                    float localForce[3];
                    calculateGrycukChainRulePairIxn_kernel( localParticle, psAChainRule[tj], localForce 
#ifdef AMOEBA_DEBUG
,  pullDebug
#endif
);
    
                    localParticle.force[0]     -= localForce[0];
                    localParticle.force[1]     -= localForce[1];
                    localParticle.force[2]     -= localForce[2];
    
                    psAChainRule[tj].force[0]  += localForce[0];
                    psAChainRule[tj].force[1]  += localForce[1];
                    psAChainRule[tj].force[2]  += localForce[2];
    
#ifdef AMOEBA_DEBUG
unsigned int index                 = (atomI == targetAtom) ? (y + tj) : atomI;
if( atomI == targetAtom || (y+tj) == targetAtom ){

        debugArray[index].x                = (float) atomI;
        debugArray[index].y                = (float) (y + tj); 
        debugArray[index].z                = -1.0f;
        debugArray[index].w                = -1.0f;

        index                             += cSim.paddedNumberOfAtoms;
        debugArray[index].x                = (float) x;
        debugArray[index].y                = (float) y;
        debugArray[index].z                = (float) tgx;
        debugArray[index].w                = -2.0f;

        index                             += cSim.paddedNumberOfAtoms;
        debugArray[index].x                = pullDebug[0].x;
        debugArray[index].y                = pullDebug[0].y;
        debugArray[index].z                = pullDebug[0].z;
        debugArray[index].w                = pullDebug[0].w;

        index                             += cSim.paddedNumberOfAtoms;
        debugArray[index].x                = pullDebug[1].x;
        debugArray[index].y                = pullDebug[1].y;
        debugArray[index].z                = pullDebug[1].z;
        debugArray[index].w                = pullDebug[1].w;

        index                             += cSim.paddedNumberOfAtoms;
        debugArray[index].x                = localForce[0];
        debugArray[index].y                = localForce[1];
        debugArray[index].z                = localForce[2];
        debugArray[index].w                = -10.0f;
}
#endif
                    calculateGrycukChainRulePairIxn_kernel( psAChainRule[tj], localParticle, localForce
#ifdef AMOEBA_DEBUG
,  pullDebug
#endif
 );
#ifdef AMOEBA_DEBUG
if( atomI == targetAtom || (y+tj) == targetAtom ){
        index                             += cSim.paddedNumberOfAtoms;
        debugArray[index].x                = pullDebug[0].x;
        debugArray[index].y                = localForce[1];
        debugArray[index].z                = localForce[2];
        debugArray[index].w                = -11.0f;
}
#endif
    
                    localParticle.force[0]     += localForce[0];
                    localParticle.force[1]     += localForce[1];
                    localParticle.force[2]     += localForce[2];
    
                    psAChainRule[tj].force[0]  -= localForce[0];
                    psAChainRule[tj].force[1]  -= localForce[1];
                    psAChainRule[tj].force[2]  -= localForce[2];
                }

                tj  = (tj + 1) & (GRID - 1);

            }

            // Write results

            float4 of;

#ifdef USE_OUTPUT_BUFFER_PER_WARP
            unsigned int offset         = x + tgx + warp*cSim.stride;
#else
            unsigned int offset         = x + tgx + (y >> GRIDBITS) * cSim.stride;
#endif
            of                          = cSim.pForce4[offset];
            of.x                       += localParticle.force[0];
            of.y                       += localParticle.force[1];
            of.z                       += localParticle.force[2];
            cSim.pForce4[offset]       = of;

#ifdef USE_OUTPUT_BUFFER_PER_WARP
            offset                      = y + tgx + warp*cSim.stride;
#else
            offset                      = y + tgx + (x >> GRIDBITS) * cSim.stride;
#endif
            of                          = cSim.pForce4[offset];
            of.x                       += sAChainRule[threadIdx.x].force[0];
            of.y                       += sAChainRule[threadIdx.x].force[1];
            of.z                       += sAChainRule[threadIdx.x].force[2];
            cSim.pForce4[offset]       = of;

            lasty = y;
        }

        pos++;
    }
}
