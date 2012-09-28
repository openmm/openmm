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
void METHOD_NAME(kCalculateAmoebaCudaKirkwood, Forces_kernel)(
                            unsigned int* workUnit){

    extern __shared__ KirkwoodParticle sA[];

    unsigned int totalWarps      = gridDim.x*blockDim.x/GRID;
    unsigned int warp            = (blockIdx.x*blockDim.x+threadIdx.x)/GRID;
    unsigned int numWorkUnits    = cSim.pInteractionCount[0];
    unsigned int pos             = warp*numWorkUnits/totalWarps;
    unsigned int end             = (warp+1)*numWorkUnits/totalWarps;
    unsigned int lasty           = 0xFFFFFFFF;

    // pWorkArray_3_1 == torque

    // pWorkArray_1_1 == dBorn
    // pWorkArray_1_2 == dBornPolar

    float energySum = 0.0f;     

    while (pos < end)
    {

        unsigned int x;
        unsigned int y;
        bool bExclusionFlag;

        // extract cell coordinates

        decodeCell( workUnit[pos], &x, &y, &bExclusionFlag );

        unsigned int tgx                         = threadIdx.x & (GRID - 1);
        unsigned int tbx                         = threadIdx.x - tgx;
        unsigned int tj                          = tgx;

        KirkwoodParticle* psA                    = &sA[tbx];
 
        unsigned int atomI                       = x + tgx;
        KirkwoodParticle localParticle;
        loadKirkwoodShared(&localParticle, atomI );

        zeroKirkwoodParticleSharedField( &localParticle );
        if (x == y) 
        {

            loadKirkwoodShared( &(sA[threadIdx.x]), atomI );
            if( atomI < cSim.atoms ){
                for (unsigned int j = 0; j < GRID && (y+j) < cSim.atoms; j++){
                    //calculateKirkwoodPairIxn_kernel( localParticle, psA[j], 0.5f, &energySum);
                    float force[3];
                    float energy;
                    calculateKirkwoodPairIxnF1_kernel( localParticle, psA[j], &energy, force);
                    calculateKirkwoodPairIxnF2_kernel( localParticle, psA[j], &energy, force);

                    localParticle.force[0]  += force[0];
                    localParticle.force[1]  += force[1];
                    localParticle.force[2]  += force[2];
                
                    energySum               += 0.5f*energy;
                
                }
            }

            localParticle.force[0]  *= 0.5f;
            localParticle.force[1]  *= 0.5f;
            localParticle.force[2]  *= 0.5f;

            // write results

#ifdef USE_OUTPUT_BUFFER_PER_WARP
            unsigned int offset                 = x + tgx + warp*cSim.paddedNumberOfAtoms;
            add3dArrayToFloat4(   offset, localParticle.force,  cSim.pForce4 );
#else
            unsigned int offset                 = x + tgx + (x >> GRIDBITS) * cSim.paddedNumberOfAtoms;
            add3dArrayToFloat4(    offset, localParticle.force,  cSim.pForce4);
#endif

            zeroKirkwoodParticleSharedField( &localParticle );
            if( atomI < cSim.atoms ){
                for (unsigned int j = 0; j < GRID && (y+j) < cSim.atoms; j++){
                    //calculateKirkwoodPairIxn_kernel( localParticle, psA[j], 0.5f, &energySum);
#ifdef INCLUDE_TORQUE

                    calculateKirkwoodPairIxnT1_kernel( localParticle, psA[j] );
                    calculateKirkwoodPairIxnT2_kernel( localParticle, psA[j] );

#else

                    float torque[3];
                    calculateKirkwoodPairIxnT1_kernel( localParticle, psA[j], torque );
                    calculateKirkwoodPairIxnT2_kernel( localParticle, psA[j], torque );
                    localParticle.force[0]  += torque[0];
                    localParticle.force[1]  += torque[1];
                    localParticle.force[2]  += torque[2];

#endif
            //        calculateKirkwoodPairIxnB1B2_kernel( localParticle, psA[j] );
                
                }
            }

            // write results

#ifdef USE_OUTPUT_BUFFER_PER_WARP
            float  of;
            offset                              = x + tgx + warp*cSim.paddedNumberOfAtoms;
/*
            of                                  = cAmoebaSim.pWorkArray_1_1[offset];
            of                                 += localParticle.dBornRadius;
            cAmoebaSim.pWorkArray_1_1[offset]   = of;

            of                                  = cAmoebaSim.pWorkArray_1_2[offset];
            of                                 += 0.5f*localParticle.dBornRadiusPolar;
            cAmoebaSim.pWorkArray_1_2[offset]   = of;
*/

#ifdef INCLUDE_TORQUE
            add3dArray( 3*offset, localParticle.torque, cAmoebaSim.pWorkArray_3_1 );
#else
            add3dArray( 3*offset, localParticle.force, cAmoebaSim.pWorkArray_3_1 );
#endif



#else


            offset                              = x + tgx + (x >> GRIDBITS) * cSim.paddedNumberOfAtoms;

/*
            cAmoebaSim.pWorkArray_1_1[offset]   = localParticle.dBornRadius;
            cAmoebaSim.pWorkArray_1_2[offset]   = 0.5f*localParticle.dBornRadiusPolar;
*/
#ifdef INCLUDE_TORQUE
            add3dArray( 3*offset, localParticle.torque, cAmoebaSim.pWorkArray_3_1 );
#else
            add3dArray( 3*offset, localParticle.force, cAmoebaSim.pWorkArray_3_1 );
#endif

#endif
            zeroKirkwoodParticleSharedField( &localParticle );
            if( atomI < cSim.atoms ){
                for (unsigned int j = 0; j < GRID && (y+j) < cSim.atoms; j++){
                    //calculateKirkwoodPairIxn_kernel( localParticle, psA[j], 0.5f, &energySum);
/*
#ifdef INCLUDE_TORQUE

                    calculateKirkwoodPairIxnT1_kernel( localParticle, psA[j] );
                    calculateKirkwoodPairIxnT2_kernel( localParticle, psA[j] );

#else

                    float torque[3];
                    calculateKirkwoodPairIxnT1_kernel( localParticle, psA[j], torque );
                    calculateKirkwoodPairIxnT2_kernel( localParticle, psA[j], torque );
                    localParticle.force[0]  += torque[0];
                    localParticle.force[1]  += torque[1];
                    localParticle.force[2]  += torque[2];

#endif
*/
                    calculateKirkwoodPairIxnB1B2_kernel( localParticle, psA[j] );
                
                }
            }

            // write results

#ifdef USE_OUTPUT_BUFFER_PER_WARP
            //float  of;
            offset                              = x + tgx + warp*cSim.paddedNumberOfAtoms;

            of                                  = cAmoebaSim.pWorkArray_1_1[offset];
            of                                 += localParticle.dBornRadius;
            cAmoebaSim.pWorkArray_1_1[offset]   = of;

            of                                  = cAmoebaSim.pWorkArray_1_2[offset];
            of                                 += 0.5f*localParticle.dBornRadiusPolar;
            cAmoebaSim.pWorkArray_1_2[offset]   = of;
/*
#ifdef INCLUDE_TORQUE
            add3dArray( 3*offset, localParticle.torque, cAmoebaSim.pWorkArray_3_1 );
#else
            add3dArray( 3*offset, localParticle.force, cAmoebaSim.pWorkArray_3_1 );
#endif
*/



#else


            offset                              = x + tgx + (x >> GRIDBITS) * cSim.paddedNumberOfAtoms;

            cAmoebaSim.pWorkArray_1_1[offset]   = localParticle.dBornRadius;
            cAmoebaSim.pWorkArray_1_2[offset]   = 0.5f*localParticle.dBornRadiusPolar;
/*
#ifdef INCLUDE_TORQUE
            add3dArray( 3*offset, localParticle.torque, cAmoebaSim.pWorkArray_3_1 );
#else
            add3dArray( 3*offset, localParticle.force, cAmoebaSim.pWorkArray_3_1 );
#endif
*/

#endif

        } else {

            if (lasty != y){
                loadKirkwoodShared( &(sA[threadIdx.x]), (y+tgx) );

            }
            zeroKirkwoodParticleSharedField( &(sA[threadIdx.x]) );

            for (unsigned int j = 0; j < GRID; j++){
                if( atomI < cSim.atoms && (y+tj) < cSim.atoms ){
                    //calculateKirkwoodPairIxn_kernel( localParticle, psA[tj], 1.0f, &energySum);
                    float force[3];
                    float energy;
                    calculateKirkwoodPairIxnF1_kernel( localParticle, psA[tj], &energy, force);
                    calculateKirkwoodPairIxnF2_kernel( localParticle, psA[tj], &energy, force);

                    localParticle.force[0]  += force[0];
                    localParticle.force[1]  += force[1];
                    localParticle.force[2]  += force[2];
                
                    psA[tj].force[0]        -= force[0];
                    psA[tj].force[1]        -= force[1];
                    psA[tj].force[2]        -= force[2];
                
                    energySum               += energy;
                
                }
                tj  = (tj + 1) & (GRID - 1);
            }

            localParticle.force[0]    *= 0.5f;
            localParticle.force[1]    *= 0.5f;
            localParticle.force[2]    *= 0.5f;

            sA[threadIdx.x].force[0]  *= 0.5f;
            sA[threadIdx.x].force[1]  *= 0.5f;
            sA[threadIdx.x].force[2]  *= 0.5f;

            // Write results

#ifdef USE_OUTPUT_BUFFER_PER_WARP

            unsigned int offset                 = x + tgx + warp*cSim.paddedNumberOfAtoms;
            add3dArrayToFloat4( offset, localParticle.force,  cSim.pForce4 );

            offset                              = y + tgx + warp*cSim.paddedNumberOfAtoms;
            add3dArrayToFloat4(   offset, sA[threadIdx.x].force,  cSim.pForce4 );
#else
            unsigned int offset                 = x + tgx + (y >> GRIDBITS) * cSim.paddedNumberOfAtoms;
            add3dArrayToFloat4(   offset, localParticle.force,  cSim.pForce4 );

            offset                              = y + tgx + (x >> GRIDBITS) * cSim.paddedNumberOfAtoms;
            add3dArrayToFloat4( offset,    sA[threadIdx.x].force,  cSim.pForce4 );

#endif

#ifndef INCLUDE_TORQUE
            zeroKirkwoodParticleSharedField( &localParticle );
            zeroKirkwoodParticleSharedField( &(sA[threadIdx.x]) );
#endif
            for (unsigned int j = 0; j < GRID; j++){
                if( atomI < cSim.atoms && (y+tj) < cSim.atoms ){
                    //calculateKirkwoodPairIxn_kernel( localParticle, psA[tj], 1.0f, &energySum);
#ifdef INCLUDE_TORQUE

                    calculateKirkwoodPairIxnT1_kernel( localParticle, psA[tj] );
                    calculateKirkwoodPairIxnT2_kernel( localParticle, psA[tj] );
                    calculateKirkwoodPairIxnT1_kernel( psA[tj], localParticle );
                    calculateKirkwoodPairIxnT2_kernel( psA[tj], localParticle );
#else

                    float torque[3];
                    calculateKirkwoodPairIxnT1_kernel( localParticle, psA[tj], torque );
                    calculateKirkwoodPairIxnT2_kernel( localParticle, psA[tj], torque );

                    localParticle.force[0]  += torque[0];
                    localParticle.force[1]  += torque[1];
                    localParticle.force[2]  += torque[2];

                    calculateKirkwoodPairIxnT1_kernel( psA[tj], localParticle, torque );
                    calculateKirkwoodPairIxnT2_kernel( psA[tj], localParticle, torque );

                    psA[tj].force[0]        += torque[0];
                    psA[tj].force[1]        += torque[1];
                    psA[tj].force[2]        += torque[2];

#endif

//                    calculateKirkwoodPairIxnB1B2_kernel( localParticle, psA[tj] );
                
                }
                tj  = (tj + 1) & (GRID - 1);
            }

            // Write results

#ifdef USE_OUTPUT_BUFFER_PER_WARP

            //float of;
            offset                              = x + tgx + warp*cSim.paddedNumberOfAtoms;
/*
            of                                  = cAmoebaSim.pWorkArray_1_1[offset];
            of                                 += localParticle.dBornRadius;
            cAmoebaSim.pWorkArray_1_1[offset]   = of;

            of                                  = cAmoebaSim.pWorkArray_1_2[offset];
            of                                 += 0.5f*localParticle.dBornRadiusPolar;
            cAmoebaSim.pWorkArray_1_2[offset]   = of;
*/

            //add3dArrayToFloat4( offset, localParticle.force,  cSim.pForce4 );
#ifdef INCLUDE_TORQUE
            add3dArray(       3*offset, localParticle.torque, cAmoebaSim.pWorkArray_3_1 );
#else
            add3dArray(       3*offset, localParticle.force, cAmoebaSim.pWorkArray_3_1 );
#endif

            offset                              = y + tgx + warp*cSim.paddedNumberOfAtoms;
/*
            of                                  = cAmoebaSim.pWorkArray_1_1[offset];
            of                                 += sA[threadIdx.x].dBornRadius;
            cAmoebaSim.pWorkArray_1_1[offset]   = of;

            of                                  = cAmoebaSim.pWorkArray_1_2[offset];
            of                                 += 0.5f*sA[threadIdx.x].dBornRadiusPolar;
            cAmoebaSim.pWorkArray_1_2[offset]   = of;
*/

            //add3dArrayToFloat4(   offset, sA[threadIdx.x].force,  cSim.pForce4 );
#ifdef INCLUDE_TORQUE
            add3dArray(         3*offset, sA[threadIdx.x].torque, cAmoebaSim.pWorkArray_3_1 );
#else
            add3dArray(         3*offset, sA[threadIdx.x].force, cAmoebaSim.pWorkArray_3_1 );
#endif

#else
            offset                              = x + tgx + (y >> GRIDBITS) * cSim.paddedNumberOfAtoms;
/*
            cAmoebaSim.pWorkArray_1_1[offset]   = localParticle.dBornRadius;
            cAmoebaSim.pWorkArray_1_2[offset]   = 0.5f*localParticle.dBornRadiusPolar;
*/

            //add3dArrayToFloat4(   offset, localParticle.force,  cSim.pForce4 );
#ifdef INCLUDE_TORQUE
            add3dArray(         3*offset, localParticle.torque, cAmoebaSim.pWorkArray_3_1 );
#else
            add3dArray(         3*offset, localParticle.force, cAmoebaSim.pWorkArray_3_1 );
#endif

            offset                              = y + tgx + (x >> GRIDBITS) * cSim.paddedNumberOfAtoms;
/*
            cAmoebaSim.pWorkArray_1_1[offset]   = sA[threadIdx.x].dBornRadius;
            cAmoebaSim.pWorkArray_1_2[offset]   = 0.5f*sA[threadIdx.x].dBornRadiusPolar;
*/

            //add3dArrayToFloat4( offset,    sA[threadIdx.x].force,  cSim.pForce4 );
#ifdef INCLUDE_TORQUE
            add3dArray(       3*offset,    sA[threadIdx.x].torque, cAmoebaSim.pWorkArray_3_1 );
#else
            add3dArray(       3*offset,    sA[threadIdx.x].force, cAmoebaSim.pWorkArray_3_1 );
#endif

#endif

#ifndef INCLUDE_TORQUE
            zeroKirkwoodParticleSharedField( &localParticle );
            zeroKirkwoodParticleSharedField( &(sA[threadIdx.x]) );
#endif
            for (unsigned int j = 0; j < GRID; j++){
                if( atomI < cSim.atoms && (y+tj) < cSim.atoms ){
/*
                    //calculateKirkwoodPairIxn_kernel( localParticle, psA[tj], 1.0f, &energySum);
#ifdef INCLUDE_TORQUE

                    calculateKirkwoodPairIxnT1_kernel( localParticle, psA[tj] );
                    calculateKirkwoodPairIxnT2_kernel( localParticle, psA[tj] );
                    calculateKirkwoodPairIxnT1_kernel( psA[tj], localParticle );
                    calculateKirkwoodPairIxnT2_kernel( psA[tj], localParticle );
#else

                    float torque[3];
                    calculateKirkwoodPairIxnT1_kernel( localParticle, psA[tj], torque );
                    calculateKirkwoodPairIxnT2_kernel( localParticle, psA[tj], torque );

                    localParticle.force[0]  += torque[0];
                    localParticle.force[1]  += torque[1];
                    localParticle.force[2]  += torque[2];

                    calculateKirkwoodPairIxnT1_kernel( psA[tj], localParticle, torque );
                    calculateKirkwoodPairIxnT2_kernel( psA[tj], localParticle, torque );

                    psA[tj].force[0]        += torque[0];
                    psA[tj].force[1]        += torque[1];
                    psA[tj].force[2]        += torque[2];

#endif
*/

                    calculateKirkwoodPairIxnB1B2_kernel( localParticle, psA[tj] );
                
                }
                tj  = (tj + 1) & (GRID - 1);
            }

            // Write results

#ifdef USE_OUTPUT_BUFFER_PER_WARP

            float of;
            offset                              = x + tgx + warp*cSim.paddedNumberOfAtoms;

            of                                  = cAmoebaSim.pWorkArray_1_1[offset];
            of                                 += localParticle.dBornRadius;
            cAmoebaSim.pWorkArray_1_1[offset]   = of;

            of                                  = cAmoebaSim.pWorkArray_1_2[offset];
            of                                 += 0.5f*localParticle.dBornRadiusPolar;
            cAmoebaSim.pWorkArray_1_2[offset]   = of;

            //add3dArrayToFloat4( offset, localParticle.force,  cSim.pForce4 );
/*         
#ifdef INCLUDE_TORQUE
            add3dArray(       3*offset, localParticle.torque, cAmoebaSim.pWorkArray_3_1 );
#else
            add3dArray(       3*offset, localParticle.force, cAmoebaSim.pWorkArray_3_1 );
#endif
*/

            offset                              = y + tgx + warp*cSim.paddedNumberOfAtoms;

            of                                  = cAmoebaSim.pWorkArray_1_1[offset];
            of                                 += sA[threadIdx.x].dBornRadius;
            cAmoebaSim.pWorkArray_1_1[offset]   = of;

            of                                  = cAmoebaSim.pWorkArray_1_2[offset];
            of                                 += 0.5f*sA[threadIdx.x].dBornRadiusPolar;
            cAmoebaSim.pWorkArray_1_2[offset]   = of;

            //add3dArrayToFloat4(   offset, sA[threadIdx.x].force,  cSim.pForce4 );
/*
#ifdef INCLUDE_TORQUE
            add3dArray(         3*offset, sA[threadIdx.x].torque, cAmoebaSim.pWorkArray_3_1 );
#else
            add3dArray(         3*offset, sA[threadIdx.x].force, cAmoebaSim.pWorkArray_3_1 );
#endif
*/

#else

            offset                              = x + tgx + (y >> GRIDBITS) * cSim.paddedNumberOfAtoms;

            cAmoebaSim.pWorkArray_1_1[offset]   = localParticle.dBornRadius;
            cAmoebaSim.pWorkArray_1_2[offset]   = 0.5f*localParticle.dBornRadiusPolar;

/*
            //add3dArrayToFloat4(   offset, localParticle.force,  cSim.pForce4 );
#ifdef INCLUDE_TORQUE
            add3dArray(         3*offset, localParticle.torque, cAmoebaSim.pWorkArray_3_1 );
#else
            add3dArray(         3*offset, localParticle.force, cAmoebaSim.pWorkArray_3_1 );
#endif
*/

            offset                              = y + tgx + (x >> GRIDBITS) * cSim.paddedNumberOfAtoms;

            cAmoebaSim.pWorkArray_1_1[offset]   = sA[threadIdx.x].dBornRadius;
            cAmoebaSim.pWorkArray_1_2[offset]   = 0.5f*sA[threadIdx.x].dBornRadiusPolar;

            //add3dArrayToFloat4( offset,    sA[threadIdx.x].force,  cSim.pForce4 );
/*            
#ifdef INCLUDE_TORQUE
            add3dArray(       3*offset,    sA[threadIdx.x].torque, cAmoebaSim.pWorkArray_3_1 );
#else
            add3dArray(       3*offset,    sA[threadIdx.x].force, cAmoebaSim.pWorkArray_3_1 );
#endif
*/

#endif


            lasty = y;
        }
        pos++;
    }
    cSim.pEnergy[blockIdx.x*blockDim.x+threadIdx.x] += 0.5f*energySum;
}
