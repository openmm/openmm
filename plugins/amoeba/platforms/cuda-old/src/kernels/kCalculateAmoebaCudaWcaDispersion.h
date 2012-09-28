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

__global__
#if (__CUDA_ARCH__ >= 200)
__launch_bounds__(384, 1)
#elif (__CUDA_ARCH__ >= 120)
__launch_bounds__(192, 1)
#else
__launch_bounds__(64, 1)
#endif

void METHOD_NAME(kCalculateAmoebaWcaDispersion, _kernel)( unsigned int* workUnit ){

    extern __shared__ WcaDispersionParticle sA[];

    unsigned int totalWarps      = gridDim.x*blockDim.x/GRID;
    unsigned int warp            = (blockIdx.x*blockDim.x+threadIdx.x)/GRID;
    unsigned int numWorkUnits    = cSim.pInteractionCount[0];
    unsigned int pos             = warp*numWorkUnits/totalWarps;
    unsigned int end             = (warp+1)*numWorkUnits/totalWarps;
    unsigned int lasty           = 0xFFFFFFFF;

    float4 jCoord;
    float jRadius;
    float jEpsilon;
    float totalEnergy            = 0.0f;

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

        WcaDispersionParticle*  psA      = &sA[tbx];
        unsigned int atomI               = x + tgx;
        float4 iCoord                    = cSim.pPosq[atomI];
        float iRadius                    = cAmoebaSim.pWcaDispersionRadiusEpsilon[atomI].x;
        float iEpsilon                   = cAmoebaSim.pWcaDispersionRadiusEpsilon[atomI].y;

        float forceSum[3];

        float emixo,emixh;
        float rmixo,rmixh;

        float emjxo,emjxh;
        float rmjxo,rmjxh;

        calculateWcaDispersionInit_kernel( iRadius,   iEpsilon,
                                           &rmixo,    &rmixh,
                                           &emixo,    &emixh );

        forceSum[0]                      = 0.0f;
        forceSum[1]                      = 0.0f;
        forceSum[2]                      = 0.0f;

        // load coordinates, charge, ...

        if (lasty != y)
        {
            loadWcaDispersionShared( &(sA[threadIdx.x]), (y+tgx), cSim.pPosq, cAmoebaSim.pWcaDispersionRadiusEpsilon );

        }

       // zero shared fields

        zeroWcaDispersionSharedForce( &(sA[threadIdx.x]) );

        for (unsigned int j = 0; j < GRID; j++)
        {

            float ijForce[3];

            // load coords, charge, ...

            loadWcaDispersionData( &(psA[tj]), &jCoord, &jRadius, &jEpsilon ); 

            // calculate force

            float energy;
            calculateWcaDispersionPairIxn_kernel( iCoord, jCoord,
                                                  iRadius,jRadius,
                                                  rmixo,  rmixh,
                                                  emixo,  emixh,
                                                  ijForce, &energy);

            if( (atomI != (y+tj)) && (atomI < cSim.atoms) && ((y+tj) < cSim.atoms) ){
       
                // add to field at atomI the field due atomJ's dipole

                forceSum[0]              += ijForce[0];
                forceSum[1]              += ijForce[1];
                forceSum[2]              += ijForce[2];
    
                // add to field at atomJ the field due atomI's dipole

                psA[tj].force[0]         -= ijForce[0];
                psA[tj].force[1]         -= ijForce[1];
                psA[tj].force[2]         -= ijForce[2];

                totalEnergy              += (x == y) ? 0.5f*energy : energy;
            }

            calculateWcaDispersionInit_kernel( jRadius,   jEpsilon,
                                               &rmjxo,    &rmjxh,
                                               &emjxo,    &emjxh );

            calculateWcaDispersionPairIxn_kernel( jCoord, iCoord,
                                                  jRadius,iRadius,
                                                  rmjxo,  rmjxh,
                                                  emjxo,  emjxh,
                                                  ijForce, &energy);

            if( (atomI != (y+tj)) && (atomI < cSim.atoms) && ((y+tj) < cSim.atoms) ){
       
                // add to field at atomI the field due atomJ's dipole

                forceSum[0]              -= ijForce[0];
                forceSum[1]              -= ijForce[1];
                forceSum[2]              -= ijForce[2];
    
                // add to field at atomJ the field due atomI's dipole

                psA[tj].force[0]         += ijForce[0];
                psA[tj].force[1]         += ijForce[1];
                psA[tj].force[2]         += ijForce[2];

                totalEnergy              += (x == y) ? 0.5f*energy : energy;
            }

            tj  = (tj + 1) & (GRID - 1);

        }

        // Write results

#ifdef USE_OUTPUT_BUFFER_PER_WARP
        unsigned int offset                 = (x + tgx + warp*cSim.paddedNumberOfAtoms);
        add3dArrayToFloat4( offset, forceSum,  cSim.pForce4);

        // include diagonal only once

        if( x != y ){
            offset                              = (y + tgx + warp*cSim.paddedNumberOfAtoms);
            add3dArrayToFloat4( offset, sA[threadIdx.x].force, cSim.pForce4);
        }

#else
        unsigned int offset                 = (x + tgx + (y >> GRIDBITS) * cSim.paddedNumberOfAtoms);
        add3dArrayToFloat4( offset, forceSum,    cSim.pForce4);

        // include diagonal only once

        if( x != y ){
            offset                              = (y + tgx + (x >> GRIDBITS) * cSim.paddedNumberOfAtoms);
            add3dArrayToFloat4( offset, sA[threadIdx.x].force,    cSim.pForce4 );
        }

#endif
        lasty = y;
        pos++;
    }

    cSim.pEnergy[blockIdx.x * blockDim.x + threadIdx.x] -= cAmoebaSim.awater*totalEnergy;
    if( (blockIdx.x*blockDim.x + threadIdx.x) == 0 ){
        cSim.pEnergy[0] += cAmoebaSim.totalMaxWcaDispersionEnergy;
    }
}
