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
__launch_bounds__(GF1XX_NONBOND_THREADS_PER_BLOCK, 1)
#elif (__CUDA_ARCH__ >= 130)
__launch_bounds__(GT2XX_NONBOND_THREADS_PER_BLOCK, 1)
#else
__launch_bounds__(G8X_NONBOND_THREADS_PER_BLOCK, 1)
#endif
void METHOD_NAME(kCalculateAmoebaWcaDispersion, _kernel)(
                            unsigned int* workUnit,
                            float4* atomCoord,
                            float2*  wcaDispersionParameters,
                            float* outputForce
#ifdef AMOEBA_DEBUG
                           , float4* debugArray, unsigned int targetAtom
#endif
){

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

#ifdef AMOEBA_DEBUG
    float4 pullDebug[3];
#endif
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
        float4 iCoord                    = atomCoord[atomI];
        float iRadius                    = wcaDispersionParameters[atomI].x;
        float iEpsilon                   = wcaDispersionParameters[atomI].y;

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
            loadWcaDispersionShared( &(sA[threadIdx.x]), (y+tgx), atomCoord, wcaDispersionParameters );

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
                                                  ijForce, &energy
#ifdef AMOEBA_DEBUG
,  pullDebug
#endif
   );

#ifdef AMOEBA_DEBUG
if( (atomI == targetAtom) || ( (y+tj) == targetAtom ) ){

unsigned int index                 = (atomI == targetAtom) ? (y + tj) : atomI;

debugArray[index].x                = (float) atomI;
debugArray[index].y                = (float) (y + tj); 
//debugArray[index].z                = (float) cAmoebaSim.numberOfAtoms;
debugArray[index].z                = atomI == (y+tj) ? 0.0f : energy;
energy                             = ( (atomI != (y+tj)) && (atomI < cAmoebaSim.numberOfAtoms) && ((y+tj) < cAmoebaSim.numberOfAtoms) )  ? (energy) : 0.0f;
debugArray[index].w                = energy+totalEnergy;

index                             += cAmoebaSim.paddedNumberOfAtoms;
debugArray[index].x                = iCoord.x;
debugArray[index].y                = iCoord.y;
debugArray[index].z                = iCoord.z;
debugArray[index].w                = (float) (blockIdx.x * blockDim.x + threadIdx.x);

index                             += cAmoebaSim.paddedNumberOfAtoms;
debugArray[index].x                = jCoord.x;
debugArray[index].y                = jCoord.y;
debugArray[index].z                = jCoord.z;
debugArray[index].w                = -4.0f;

index                             += cAmoebaSim.paddedNumberOfAtoms;
debugArray[index].x                = emixo;
debugArray[index].y                = emixh;
debugArray[index].z                = rmixo;
debugArray[index].w                = rmixh;

index                             += cAmoebaSim.paddedNumberOfAtoms;
debugArray[index].x                = pullDebug[0].x;
debugArray[index].y                = pullDebug[0].y;
debugArray[index].z                = pullDebug[0].z;
debugArray[index].w                = pullDebug[0].w;

#if 0
index                             += cAmoebaSim.paddedNumberOfAtoms;
debugArray[index].x                = pullDebug[1].x;
debugArray[index].y                = pullDebug[1].y;
debugArray[index].z                = pullDebug[1].z;
debugArray[index].w                = pullDebug[1].w;

index                             += cAmoebaSim.paddedNumberOfAtoms;
debugArray[index].x                = pullDebug[2].x;
debugArray[index].y                = pullDebug[2].y;
debugArray[index].z                = pullDebug[2].z;
debugArray[index].w                = pullDebug[2].w;
#endif

} else {
//    energy = 0.0f;
}
#endif
            if( (atomI != (y+tj)) && (atomI < cAmoebaSim.numberOfAtoms) && ((y+tj) < cAmoebaSim.numberOfAtoms) ){
       
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
                                                  ijForce, &energy
#ifdef AMOEBA_DEBUG
,  pullDebug
#endif
   );

#ifdef AMOEBA_DEBUG
if( (atomI == targetAtom) || ( (y+tj) == targetAtom ) ){

unsigned int index                 = (atomI == targetAtom) ? (y + tj) : atomI;
index                             += 2*cAmoebaSim.paddedNumberOfAtoms;

debugArray[index].x                = (float) atomI;
debugArray[index].y                = (float) (y + tj); 
debugArray[index].z                = atomI == (y+tj) ? 0.0f : energy;
energy                             = ( (atomI != (y+tj)) && (atomI < cAmoebaSim.numberOfAtoms) && ((y+tj) < cAmoebaSim.numberOfAtoms) )  ? (energy) : 0.0f;
debugArray[index].w                = energy+totalEnergy;
//debugArray[index].w                = -2.0f;

index                             += cAmoebaSim.paddedNumberOfAtoms;
debugArray[index].x                = emjxo;
debugArray[index].y                = emjxh;
debugArray[index].z                = rmjxo;
debugArray[index].w                = rmjxh;

index                             += cAmoebaSim.paddedNumberOfAtoms;
debugArray[index].x                = pullDebug[0].x;
debugArray[index].y                = pullDebug[0].y;
debugArray[index].z                = pullDebug[0].z;
debugArray[index].w                = pullDebug[0].w;
#if 0
index                             += cAmoebaSim.paddedNumberOfAtoms;
debugArray[index].x                = pullDebug[1].x;
debugArray[index].y                = pullDebug[1].y;
debugArray[index].z                = pullDebug[1].z;
debugArray[index].w                = pullDebug[1].w;

index                             += cAmoebaSim.paddedNumberOfAtoms;
debugArray[index].x                = pullDebug[2].x;
debugArray[index].y                = pullDebug[2].y;
debugArray[index].z                = pullDebug[2].z;
debugArray[index].w                = pullDebug[2].w;
#endif

} else {
    //energy = 0.0f;
}
#endif
            if( (atomI != (y+tj)) && (atomI < cAmoebaSim.numberOfAtoms) && ((y+tj) < cAmoebaSim.numberOfAtoms) ){
       
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
        unsigned int offset                 = 3*(x + tgx + warp*cAmoebaSim.paddedNumberOfAtoms);
        load3dArrayBufferPerWarp( offset, forceSum,       outputForce );

        offset                              = 3*(y + tgx + warp*cAmoebaSim.paddedNumberOfAtoms);

        load3dArrayBufferPerWarp( offset, sA[threadIdx.x].force,       outputForce );

#else
        unsigned int offset                 = 3*(x + tgx + (y >> GRIDBITS) * cAmoebaSim.paddedNumberOfAtoms);
        load3dArray( offset, forceSum,       outputForce );

        offset                              = 3*(y + tgx + (x >> GRIDBITS) * cAmoebaSim.paddedNumberOfAtoms);
        load3dArray( offset, sA[threadIdx.x].force,       outputForce );

#endif
        lasty = y;
        pos++;
    }
    cSim.pEnergy[blockIdx.x * blockDim.x + threadIdx.x] -= cAmoebaSim.awater*totalEnergy;
    if( (blockIdx.x*blockDim.x + threadIdx.x) == 0 ){
        cSim.pEnergy[0] += cAmoebaSim.totalMaxWcaDispersionEnergy;
    }
}
