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
__launch_bounds__(128, 1)
#else
__launch_bounds__(64, 1)
#endif
void METHOD_NAME(kCalculateAmoebaCudaElectrostaticPotential, _kernel)( void ){

    extern __shared__ volatile ElectrostaticPotentialParticle sAPotential[];

    unsigned int* workUnit       = cAmoebaSim.pPotentialWorkUnit;
    unsigned int totalWarps      = gridDim.x*blockDim.x/GRID;
    unsigned int warp            = (blockIdx.x*blockDim.x+threadIdx.x)/GRID;
    unsigned int numWorkUnits    = cAmoebaSim.potentialWorkUnits;
    unsigned int pos             = warp*numWorkUnits/totalWarps;
    unsigned int end             = (warp+1)*numWorkUnits/totalWarps;

    while (pos < end){

        unsigned int x;
        unsigned int y;
        bool bExclusionFlag;

        // Extract cell coordinates

        decodeCell( workUnit[pos], &x, &y, &bExclusionFlag );

        unsigned int tgx              = threadIdx.x & (GRID - 1);
        unsigned int tbx              = threadIdx.x - tgx;
        unsigned int tj               = tgx;

        volatile ElectrostaticPotentialParticle* psA = &sAPotential[tbx];
        unsigned int gridPointIndex   = x + tgx;
        unsigned int particleIndex    = y + tgx;

        // load particle info

        loadElectrostaticPotentialParticle( &(sAPotential[threadIdx.x]), particleIndex );

        float totalPotential  = 0.0f;
        for (unsigned int j = 0; j < GRID; j++){
            unsigned int particleJ = y + tj;
            float potential;
            calculateElectrostaticPotentialForAtomGridPoint_kernel( psA[tj], cAmoebaSim.pPotentialGrid[gridPointIndex], &potential );

            if( particleJ < cSim.atoms && gridPointIndex < cAmoebaSim.potentialGridSize ){
                totalPotential += potential;
            }

            tj = (tj + 1) & (GRID - 1);
        }

        // Write results

#ifdef USE_OUTPUT_BUFFER_PER_WARP
        unsigned int offset            = (x + tgx + warp*cAmoebaSim.paddedPotentialGridSize);
        cAmoebaSim.pPotential[offset] += totalPotential; 
#else
        unsigned int offset            = (x + tgx + (y >> GRIDBITS)*cAmoebaSim.paddedPotentialGridSize);
        cAmoebaSim.pPotential[offset]  = totalPotential; 
#endif
        pos++;
    }
}
