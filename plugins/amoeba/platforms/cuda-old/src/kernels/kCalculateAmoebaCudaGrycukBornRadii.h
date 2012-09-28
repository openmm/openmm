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
void METHOD_NAME(kCalculateAmoebaGrycukBornRadii, _kernel)( unsigned int* workUnit ){

    extern __shared__ GrycukParticle sA[];

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

        unsigned int tgx                 = threadIdx.x & (GRID - 1);
        unsigned int tbx                 = threadIdx.x - tgx;
        unsigned int tj                  = tgx;

        GrycukParticle*  psA             = &sA[tbx];
        unsigned int atomI               = x + tgx;
        GrycukParticle localParticle;
        loadGrycukShared( &localParticle, atomI );

        float bornSum                    = 0.0f;

        if (x == y) // Handle diagonals uniquely at 50% efficiency
        {

            // load shared data

            loadGrycukShared( &(sA[threadIdx.x]), atomI );

            for (unsigned int j = 0; j < GRID; j++)
            {
                float localBornSum;
                calculateGrycukBornRadiiPairIxn_kernel( localParticle, psA[j], &localBornSum );
                bornSum   +=  ( (atomI == (y + j)) || (atomI >= cSim.atoms) || ((y+j) >= cSim.atoms) ) ? 0.0 : localBornSum;
            }

            // Write results

#ifdef USE_OUTPUT_BUFFER_PER_WARP
            unsigned int offset    = x + tgx + warp*cSim.stride;
            cSim.pBornSum[offset] += bornSum;
#else
            unsigned int offset   = x + tgx + (y >> GRIDBITS) * cSim.stride;
            cSim.pBornSum[offset] = bornSum;
#endif

        } else {

            if (lasty != y) {
                unsigned int atomJ        = y + tgx;
                loadGrycukShared( &(sA[threadIdx.x]), atomJ );
            }

           // zero shared fields

            zeroGrycukParticleSharedField(  &(sA[threadIdx.x]) );

            for (unsigned int j = 0; j < GRID; j++)
            {

                float localBornSum;
                calculateGrycukBornRadiiPairIxn_kernel( localParticle, psA[tj], &localBornSum );
                bornSum           +=  ( (atomI >= cSim.atoms) || ((y+tj) >= cSim.atoms) ) ? 0.0 : localBornSum;

                calculateGrycukBornRadiiPairIxn_kernel( psA[tj], localParticle, &localBornSum );
                psA[tj].bornSum   +=  ( (atomI >= cSim.atoms) || ((y+tj) >= cSim.atoms) ) ? 0.0 : localBornSum;

                tj                 = (tj + 1) & (GRID - 1);

            }

            // Write results

#ifdef USE_OUTPUT_BUFFER_PER_WARP

            unsigned int offset    = x + tgx + warp*cSim.stride;
            cSim.pBornSum[offset] += bornSum;

            offset = y + tgx + warp*cSim.stride;
            cSim.pBornSum[offset] += sA[threadIdx.x].bornSum;
#else

            unsigned int offset   = x + tgx + (y >> GRIDBITS) * cSim.stride;
            cSim.pBornSum[offset] = bornSum;

            offset = y + tgx + (x >> GRIDBITS) * cSim.stride;
            cSim.pBornSum[offset] = sA[threadIdx.x].bornSum;
#endif

            lasty = y;
        }

        pos++;
    }
}
