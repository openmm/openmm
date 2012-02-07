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
#elif (__CUDA_ARCH__ >= 120)
__launch_bounds__(GT2XX_NONBOND_THREADS_PER_BLOCK, 1)
#else
__launch_bounds__(G8X_NONBOND_THREADS_PER_BLOCK, 1)
#endif
void METHOD_NAME(kCalculateAmoebaVdw14_7, _kernel)(
                            unsigned int* workUnit,
                            float4* atomCoord,
                            float2*  vdwParameters,
                            int sigmaCombiningRule,
                            int epsilonCombiningRule,
                            float* outputForce
){

    extern __shared__ Vdw14_7Particle sA[];

    unsigned int totalWarps      = gridDim.x*blockDim.x/GRID;
    unsigned int warp            = (blockIdx.x*blockDim.x+threadIdx.x)/GRID;
    unsigned int numWorkUnits    = cSim.pInteractionCount[0];
    unsigned int pos             = warp*numWorkUnits/totalWarps;
    unsigned int end             = (warp+1)*numWorkUnits/totalWarps;
    unsigned int lasty           = 0xFFFFFFFF;

    int exclusionIndex;
    int exclusionMask;
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

        Vdw14_7Particle* psA             = &sA[tbx];

        Vdw14_7Particle  localParticle;
        unsigned int atomI               = x + tgx;
        loadVdw14_7Shared( &localParticle, atomI, atomCoord, vdwParameters );
 
        float forceSum[3];

        // forceSum:      field at i due to j
        // fieldPolarSum: field at i due to j polar

        forceSum[0]                      = 0.0f;
        forceSum[1]                      = 0.0f;
        forceSum[2]                      = 0.0f;

        if (x == y) 
        {

            if( bExclusionFlag ){
                unsigned int xi              = x >> GRIDBITS;
                unsigned int cell            = xi + xi*cSim.paddedNumberOfAtoms/GRID-xi*(xi+1)/2;
                exclusionIndex               = cAmoebaSim.pVdwExclusionIndicesIndex[cell]+tgx;
                exclusionMask                = cAmoebaSim.pVdwExclusionIndices[exclusionIndex];
            }

            // load shared data

            loadVdw14_7Shared( &(sA[threadIdx.x]), atomI, atomCoord, vdwParameters );

            for (unsigned int j = 0; j < GRID; j++)
            {

                float ijForce[3];

                // get combined sigma and epsilon

                float combindedSigma;
                float combindedEpsilon;
                getVdw14_7CombindedSigmaEpsilon_kernel( sigmaCombiningRule,   localParticle.sigma,   psA[j].sigma,   &combindedSigma,
                                                        epsilonCombiningRule, localParticle.epsilon, psA[j].epsilon, &combindedEpsilon );
 
                // calculate force

                ijForce[0]    = psA[j].x - localParticle.x;
                ijForce[1]    = psA[j].y - localParticle.y;
                ijForce[2]    = psA[j].z - localParticle.z;
                if( cAmoebaSim.vdwUsePBC )
                {
                    ijForce[0]   -= floorf(ijForce[0]*cSim.invPeriodicBoxSizeX+0.5f)*cSim.periodicBoxSizeX;
                    ijForce[1]   -= floorf(ijForce[1]*cSim.invPeriodicBoxSizeY+0.5f)*cSim.periodicBoxSizeY;
                    ijForce[2]   -= floorf(ijForce[2]*cSim.invPeriodicBoxSizeZ+0.5f)*cSim.periodicBoxSizeZ;
                }

                float energy;
                calculateVdw14_7PairIxn_kernel( combindedSigma, combindedEpsilon, ijForce, &energy);
                // mask out excluded ixns

                unsigned int mask  =  ( (atomI >= cSim.atoms) || ((y+j) >= cSim.atoms) ) ? 0 : 1;
                if( mask && bExclusionFlag ){
                    unsigned int maskIndex  = 1 << j;
                    mask = (exclusionMask & maskIndex) ? 0 : 1;
                }

                // add to field at atomI the field due atomJ's dipole

                forceSum[0]            += mask ? ijForce[0]  : 0.0f;
                forceSum[1]            += mask ? ijForce[1]  : 0.0f;
                forceSum[2]            += mask ? ijForce[2]  : 0.0f;
                totalEnergy            += mask ? 0.5f*energy : 0.0f;
            }

            // Write results

#ifdef USE_OUTPUT_BUFFER_PER_WARP

            unsigned int offset                 = 3*(x + tgx + warp*cSim.paddedNumberOfAtoms);
            load3dArrayBufferPerWarp( offset, forceSum, outputForce );

#else
            unsigned int offset                 = 3*(x + tgx + (x >> GRIDBITS) * cSim.paddedNumberOfAtoms);
            load3dArray( offset, forceSum, outputForce );
#endif

        } else {

            if (lasty != y)
            {
                // load coordinates, charge, ...

                loadVdw14_7Shared( &(sA[threadIdx.x]), (y+tgx), atomCoord, vdwParameters );

            }

#ifdef USE_CUTOFF
            unsigned int flags = cSim.pInteractionFlag[pos];
            if (flags == 0) {
            } else {
#endif
               // zero shared fields
    
                zeroVdw14_7SharedForce( &(sA[threadIdx.x]) );
    
                if( bExclusionFlag )
                {
    
                    unsigned int xi              = x >> GRIDBITS;
                    unsigned int yi              = y >> GRIDBITS;
                    unsigned int cell            = xi+yi*cSim.paddedNumberOfAtoms/GRID-yi*(yi+1)/2;
    
                    exclusionIndex               = cAmoebaSim.pVdwExclusionIndicesIndex[cell]+tgx;
                    exclusionMask                = cAmoebaSim.pVdwExclusionIndices[exclusionIndex];
                }
    
                for (unsigned int j = 0; j < GRID; j++)
                {
        
                    float ijForce[3];
#ifdef USE_CUTOFF
                    if ((flags&(1<<j)) != 0)
                    {
                        unsigned int jIdx  = (flags == 0xFFFFFFFF) ? tj : j;
#else
                        unsigned int jIdx  = tj;
#endif
            
                        // get combined sigma and epsilon

                        float combindedSigma;
                        float combindedEpsilon;
                        getVdw14_7CombindedSigmaEpsilon_kernel( sigmaCombiningRule,   localParticle.sigma,   psA[jIdx].sigma,   &combindedSigma,
                                                                epsilonCombiningRule, localParticle.epsilon, psA[jIdx].epsilon, &combindedEpsilon );
            
                        // calculate force
            
                        float energy;
                        ijForce[0]    = psA[jIdx].x - localParticle.x;
                        ijForce[1]    = psA[jIdx].y - localParticle.y;
                        ijForce[2]    = psA[jIdx].z - localParticle.z;
                        if( cAmoebaSim.vdwUsePBC )
                        {
                            ijForce[0]   -= floorf(ijForce[0]*cSim.invPeriodicBoxSizeX+0.5f)*cSim.periodicBoxSizeX;
                            ijForce[1]   -= floorf(ijForce[1]*cSim.invPeriodicBoxSizeY+0.5f)*cSim.periodicBoxSizeY;
                            ijForce[2]   -= floorf(ijForce[2]*cSim.invPeriodicBoxSizeZ+0.5f)*cSim.periodicBoxSizeZ;
                        }
                        calculateVdw14_7PairIxn_kernel( combindedSigma, combindedEpsilon, ijForce, &energy);
            
                        // mask out excluded ixns
        
                        unsigned int mask       =  ( (atomI >= cSim.atoms) || ((y+jIdx) >= cSim.atoms) ) ? 0 : 1;
                        if( mask && bExclusionFlag ){
                            unsigned int maskIndex  = 1 << jIdx;
                            mask =  (exclusionMask & maskIndex) ? 0 : 1;
                        }
                        if( mask == 0 )
                        {
                            energy = ijForce[0] = ijForce[1] = ijForce[2] = 0.0f;
                        }
                       
                        // accumulate force for atomI
                
                        forceSum[0]        += ijForce[0];
                        forceSum[1]        += ijForce[1];
                        forceSum[2]        += ijForce[2];

                        totalEnergy        += energy;
                    
#ifndef USE_CUTOFF
                        psA[jIdx].force[0]   -= ijForce[0];
                        psA[jIdx].force[1]   -= ijForce[1];
                        psA[jIdx].force[2]   -= ijForce[2];
#else
                        if( flags == 0xFFFFFFFF ){
                
                            psA[jIdx].force[0]   -= ijForce[0];
                            psA[jIdx].force[1]   -= ijForce[1];
                            psA[jIdx].force[2]   -= ijForce[2];

                        } else {

                            sA[threadIdx.x].tempForce[0]     = ijForce[0];
                            sA[threadIdx.x].tempForce[1]     = ijForce[1];
                            sA[threadIdx.x].tempForce[2]     = ijForce[2];

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

                            }

                        }
#endif
        
                
#ifdef USE_CUTOFF
                    }
#endif
                    tj                  = (tj + 1) & (GRID - 1);
        
                } // end of j-loop 
#ifdef USE_CUTOFF
            }
#endif

            // Write results

#ifdef USE_OUTPUT_BUFFER_PER_WARP
            unsigned int offset                 = 3*(x + tgx + warp*cSim.paddedNumberOfAtoms);
            load3dArrayBufferPerWarp( offset, forceSum,       outputForce );

            offset                              = 3*(y + tgx + warp*cSim.paddedNumberOfAtoms);

            load3dArrayBufferPerWarp( offset, sA[threadIdx.x].force,       outputForce );
#else
            unsigned int offset                 = 3*(x + tgx + (y >> GRIDBITS) * cSim.paddedNumberOfAtoms);
            load3dArray( offset, forceSum,       outputForce );

            offset                              = 3*(y + tgx + (x >> GRIDBITS) * cSim.paddedNumberOfAtoms);
            load3dArray( offset, sA[threadIdx.x].force,       outputForce );

#endif
            lasty = y;

        } // x == y block
        pos++;
    }
    cSim.pEnergy[blockIdx.x * blockDim.x + threadIdx.x] += totalEnergy;
}
