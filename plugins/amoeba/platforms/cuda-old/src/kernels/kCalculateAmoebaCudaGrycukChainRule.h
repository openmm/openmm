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
void METHOD_NAME(kCalculateAmoebaGrycukChainRule, _kernel)( unsigned int* workUnit ){

    extern __shared__ GrycukChainRuleParticle sAChainRule[];

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
                calculateGrycukChainRulePairIxn_kernel( localParticle, psAChainRule[j], localForce);
                if( (atomI != (y + j)) && (atomI < cSim.atoms) && ((y+j) < cSim.atoms) ){

                    localParticle.force[0]     -= localForce[0];
                    localParticle.force[1]     -= localForce[1];
                    localParticle.force[2]     -= localForce[2];

                    psAChainRule[j].force[0]   += localForce[0];
                    psAChainRule[j].force[1]   += localForce[1];
                    psAChainRule[j].force[2]   += localForce[2];


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
                    calculateGrycukChainRulePairIxn_kernel( localParticle, psAChainRule[tj], localForce );
    
                    localParticle.force[0]     -= localForce[0];
                    localParticle.force[1]     -= localForce[1];
                    localParticle.force[2]     -= localForce[2];
    
                    psAChainRule[tj].force[0]  += localForce[0];
                    psAChainRule[tj].force[1]  += localForce[1];
                    psAChainRule[tj].force[2]  += localForce[2];
    
                    calculateGrycukChainRulePairIxn_kernel( psAChainRule[tj], localParticle, localForce);
    
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
