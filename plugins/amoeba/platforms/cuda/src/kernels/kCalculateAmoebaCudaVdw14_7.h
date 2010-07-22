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
void METHOD_NAME(kCalculateAmoebaVdw14_7, _kernel)(
                            unsigned int* workUnit,
                            float4* atomCoord,
                            float2*  vdwParameters,
                            int sigmaCombiningRule,
                            int epsilonCombiningRule,
                            float* outputForce
#ifdef AMOEBA_DEBUG
                           , float4* debugArray, unsigned int targetAtom
#endif
){

    extern __shared__ Vdw14_7Particle sA[];

    unsigned int totalWarps      = gridDim.x*blockDim.x/GRID;
    unsigned int warp            = (blockIdx.x*blockDim.x+threadIdx.x)/GRID;
    unsigned int numWorkUnits    = cSim.pInteractionCount[0];
    unsigned int pos             = warp*numWorkUnits/totalWarps;
    unsigned int end             = (warp+1)*numWorkUnits/totalWarps;
    unsigned int lasty           = 0xFFFFFFFF;

    float4 jCoord;
    float jSigma;
    float jEpsilon;
    float totalEnergy            = 0.0f;

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

        unsigned int tgx                 = threadIdx.x & (GRID - 1);
        unsigned int tbx                 = threadIdx.x - tgx;
        unsigned int tj                  = tgx;

        Vdw14_7Particle*  psA            = &sA[tbx];
        unsigned int atomI               = x + tgx;
        float4 iCoord                    = atomCoord[atomI];
        float iSigma                     = vdwParameters[atomI].x;
        float iEpsilon                   = vdwParameters[atomI].y;

        float forceSum[3];

        // forceSum:      field at i due to j
        // fieldPolarSum: field at i due to j polar

        forceSum[0]                      = 0.0f;
        forceSum[1]                      = 0.0f;
        forceSum[2]                      = 0.0f;

        if (x == y) 
        {

            unsigned int xi              = x >> GRIDBITS;
            unsigned int cell            = xi + xi*cAmoebaSim.paddedNumberOfAtoms/GRID-xi*(xi+1)/2;
            int exclusionIndex           = cAmoebaSim.pVdwExclusionIndicesIndex[cell]+tgx;
            int exclusionMask            = cAmoebaSim.pVdwExclusionIndices[exclusionIndex];

            // load shared data

            loadVdw14_7Shared( &(sA[threadIdx.x]), atomI, atomCoord, vdwParameters );

            for (unsigned int j = 0; j < GRID; j++)
            {

                float ijForce[3];

                // load coords, charge, ...

                loadVdw14_7Data( &(psA[j]), &jCoord, &jSigma, &jEpsilon ); 

                // get combined sigma and epsilon

                float combindedSigma;
                float combindedEpsilon;
                getVdw14_7CombindedSigmaEpsilon_kernel( sigmaCombiningRule,   iSigma,   jSigma,   &combindedSigma,
                                                        epsilonCombiningRule, iEpsilon, jEpsilon, &combindedEpsilon );
 
                // calculate force

                float energy;
                calculateVdw14_7PairIxn_kernel( iCoord, jCoord, combindedSigma, combindedEpsilon, ijForce, &energy
#ifdef AMOEBA_DEBUG
,  pullDebug
#endif
);

                // mask out excluded ixns

                unsigned int maskIndex  = 1 << j;
                unsigned int mask       =  ( (exclusionMask & maskIndex) || (atomI >= cAmoebaSim.numberOfAtoms) || ((y+j) >= cAmoebaSim.numberOfAtoms) ) ? 0 : 1;

                // add to field at atomI the field due atomJ's dipole

                forceSum[0]            += mask ? ijForce[0] : 0.0f;
                forceSum[1]            += mask ? ijForce[1] : 0.0f;
                forceSum[2]            += mask ? ijForce[2] : 0.0f;
                totalEnergy            += mask ? 0.5*energy : 0.0f;

#ifdef AMOEBA_DEBUG
if( atomI == targetAtom || (y+j) == targetAtom ){
        unsigned int index                 = (atomI == targetAtom) ? (y + j) : atomI;

        debugArray[index].x                = (float) atomI;
        debugArray[index].y                = (float) (y + j); 
        debugArray[index].z                = -1.0f;
        debugArray[index].w                = (float) (mask + 1); 

        index                             += cAmoebaSim.paddedNumberOfAtoms;
        debugArray[index].x                = (float) x;
        debugArray[index].y                = (float) y;
        debugArray[index].z                = (float) cell+tgx;
        debugArray[index].w                = energy;

        index                             += cAmoebaSim.paddedNumberOfAtoms;
        debugArray[index].x                = pullDebug[0].x;
        debugArray[index].y                = pullDebug[0].y;
        debugArray[index].z                = pullDebug[0].z;
        debugArray[index].w                = pullDebug[0].w;

        index                             += cAmoebaSim.paddedNumberOfAtoms;
        debugArray[index].x                = pullDebug[1].x;
        debugArray[index].y                = pullDebug[1].y;
        debugArray[index].z                = pullDebug[1].z;
        debugArray[index].w                = pullDebug[1].w;

        index                             += cAmoebaSim.paddedNumberOfAtoms;
        debugArray[index].x                = mask ? ijForce[0] : 0.0f;
        debugArray[index].y                = mask ? ijForce[1] : 0.0f;
        debugArray[index].z                = mask ? ijForce[2] : 0.0f;
}
#endif

            }

            // Write results

#ifdef USE_OUTPUT_BUFFER_PER_WARP

            unsigned int offset                 = 3*(x + tgx + warp*cAmoebaSim.paddedNumberOfAtoms);
            load3dArrayBufferPerWarp( offset, forceSum, outputForce );

#else
            unsigned int offset                   = 3*(x + tgx + (x >> GRIDBITS) * cAmoebaSim.paddedNumberOfAtoms);
            load3dArray( offset, forceSum, outputForce );
#endif

        }
        else
        {
            // Read fixed atom data into registers and GRF
            if (lasty != y)
            {
                // load coordinates, charge, ...

                loadVdw14_7Shared( &(sA[threadIdx.x]), (y+tgx), atomCoord, vdwParameters );

            }

           // zero shared fields

            zeroVdw14_7SharedForce(  &(sA[threadIdx.x]) );

            if( !bExclusionFlag )
            {
                for (unsigned int j = 0; j < GRID; j++)
                {
    
                    float ijForce[3];
    
                    // load coords, charge, ...
    
                    loadVdw14_7Data( &(psA[tj]), &jCoord, &jSigma, &jEpsilon ); 
    
                    // get combined sigma and epsilon
    
                    float combindedSigma;
                    float combindedEpsilon;
                    getVdw14_7CombindedSigmaEpsilon_kernel( sigmaCombiningRule,   iSigma,   jSigma,   &combindedSigma,
                                                            epsilonCombiningRule, iEpsilon, jEpsilon, &combindedEpsilon );
    
                    // calculate force
    
                    float energy;
                    calculateVdw14_7PairIxn_kernel( iCoord, jCoord, combindedSigma, combindedEpsilon, ijForce, &energy
#ifdef AMOEBA_DEBUG
    ,  pullDebug
#endif
       );
    
                    if( (atomI < cAmoebaSim.numberOfAtoms) && ((y+tj) < cAmoebaSim.numberOfAtoms) ){
               
                        // add to field at atomI the field due atomJ's dipole
        
                        forceSum[0]              += ijForce[0];
                        forceSum[1]              += ijForce[1];
                        forceSum[2]              += ijForce[2];
            
                        // add to field at atomJ the field due atomI's dipole
        
                        psA[tj].force[0]         -= ijForce[0];
                        psA[tj].force[1]         -= ijForce[1];
                        psA[tj].force[2]         -= ijForce[2];

                        totalEnergy              += energy;
        
                    }
    
#ifdef AMOEBA_DEBUG
if( atomI == targetAtom || (y+tj) == targetAtom ){
        unsigned int index                 = (atomI == targetAtom) ? (y + tj) : atomI;

        debugArray[index].x                = (float) atomI;
        debugArray[index].y                = (float) (y + tj); 
        debugArray[index].z                = -2.0f;
        debugArray[index].w                = -1.0f;

        index                             += cAmoebaSim.paddedNumberOfAtoms;
        debugArray[index].x                = (float) x;
        debugArray[index].y                = (float) y;
        debugArray[index].z                = -1.0f;
        debugArray[index].w                = energy;

        index                             += cAmoebaSim.paddedNumberOfAtoms;
        debugArray[index].x                = pullDebug[0].x;
        debugArray[index].y                = pullDebug[0].y;
        debugArray[index].z                = pullDebug[0].z;
        debugArray[index].w                = pullDebug[0].w;

        index                             += cAmoebaSim.paddedNumberOfAtoms;
        debugArray[index].x                = pullDebug[1].x;
        debugArray[index].y                = pullDebug[1].y;
        debugArray[index].z                = pullDebug[1].z;
        debugArray[index].w                = pullDebug[1].w;

        index                             += cAmoebaSim.paddedNumberOfAtoms;
        debugArray[index].x                = ijForce[0];
        debugArray[index].y                = ijForce[1];
        debugArray[index].z                = ijForce[2];
}
#endif
                    tj                  = (tj + 1) & (GRID - 1);
    
                }
            }
            else 
            {
                unsigned int xi              = x >> GRIDBITS;
                unsigned int yi              = y >> GRIDBITS;
                unsigned int cell            = xi+yi*cSim.paddedNumberOfAtoms/GRID-yi*(yi+1)/2;

                int exclusionIndex           = cAmoebaSim.pVdwExclusionIndicesIndex[cell]+tgx;
                int exclusionMask            = cAmoebaSim.pVdwExclusionIndices[exclusionIndex];
                for (unsigned int j = 0; j < GRID; j++)
                {
    
                    float ijForce[3];
    
                    // load coords, charge, ...
    
                    loadVdw14_7Data( &(psA[tj]), &jCoord, &jSigma, &jEpsilon ); 
    
                    // get combined sigma and epsilon
    
                    float combindedSigma;
                    float combindedEpsilon;
                    getVdw14_7CombindedSigmaEpsilon_kernel( sigmaCombiningRule,   iSigma,   jSigma,   &combindedSigma,
                                                            epsilonCombiningRule, iEpsilon, jEpsilon, &combindedEpsilon );
    
                    // calculate force
    
                    float energy;
                    calculateVdw14_7PairIxn_kernel( iCoord, jCoord, combindedSigma, combindedEpsilon, ijForce, &energy
#ifdef AMOEBA_DEBUG
    ,  pullDebug
#endif
       );
    
                    // mask out excluded ixns

                    unsigned int maskIndex  = 1 << tj;
                    unsigned int mask       =  ( (exclusionMask & maskIndex) || (atomI >= cAmoebaSim.numberOfAtoms) || ((y+tj) >= cAmoebaSim.numberOfAtoms) ) ? 0 : 1;
               
                    // accumulate force for atomI
        
                    forceSum[0]        += mask ? ijForce[0] : 0.0f;
                    forceSum[1]        += mask ? ijForce[1] : 0.0f;
                    forceSum[2]        += mask ? ijForce[2] : 0.0f;
            
                    // accumulate force for atomJ
        
                    psA[tj].force[0]   -= mask ? ijForce[0] : 0.0f;
                    psA[tj].force[1]   -= mask ? ijForce[1] : 0.0f;
                    psA[tj].force[2]   -= mask ? ijForce[2] : 0.0f;

                    totalEnergy        += mask ? energy     : 0.0f;
        
#ifdef AMOEBA_DEBUG
if( atomI == targetAtom || (y+tj) == targetAtom ){
        unsigned int index                 = (atomI == targetAtom) ? (y + tj) : atomI;

        debugArray[index].x                = (float) atomI;
        debugArray[index].y                = (float) (y + tj); 
        debugArray[index].z                = -3.0;
        debugArray[index].w                = (float) (mask + 1); 

        index                             += cAmoebaSim.paddedNumberOfAtoms;
        debugArray[index].x                = (float) x;
        debugArray[index].y                = (float) y;
        debugArray[index].z                = (float) cell+tgx;
        debugArray[index].w                = energy;

        index                             += cAmoebaSim.paddedNumberOfAtoms;
        debugArray[index].x                = pullDebug[0].x;
        debugArray[index].y                = pullDebug[0].y;
        debugArray[index].z                = pullDebug[0].z;
        debugArray[index].w                = pullDebug[0].w;

        index                             += cAmoebaSim.paddedNumberOfAtoms;
        debugArray[index].x                = pullDebug[1].x;
        debugArray[index].y                = pullDebug[1].y;
        debugArray[index].z                = pullDebug[1].z;
        debugArray[index].w                = pullDebug[1].w;

        index                             += cAmoebaSim.paddedNumberOfAtoms;
        debugArray[index].x                = mask ? ijForce[0] : 0.0f;
        debugArray[index].y                = mask ? ijForce[1] : 0.0f;
        debugArray[index].z                = mask ? ijForce[2] : 0.0f;
}
#endif
                    tj                  = (tj + 1) & (GRID - 1);
    
                } 
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
        }
        pos++;
    }
    cSim.pEnergy[blockIdx.x * blockDim.x + threadIdx.x] += totalEnergy;
}
