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
__launch_bounds__(384, 1)
#elif (__CUDA_ARCH__ >= 120)
__launch_bounds__(128, 1)
#else
__launch_bounds__(64, 1)
#endif
void METHOD_NAME(kCalculateAmoebaCudaElectrostatic, Forces_kernel)(
                            unsigned int* workUnit,
                            float4* atomCoord,
                            float* labFrameDipole,
                            float* labFrameQuadrupole,
                            float* inducedDipole,
                            float* inducedDipolePolar,
                            float* outputForce,
                            float* outputTorque

#ifdef AMOEBA_DEBUG
                           , float4* debugArray, unsigned int targetAtom
#endif
){

#ifdef AMOEBA_DEBUG
    float4 pullBack[20];
#endif

    extern __shared__ ElectrostaticParticle sA[];

    unsigned int totalWarps      = gridDim.x*blockDim.x/GRID;
    unsigned int warp            = (blockIdx.x*blockDim.x+threadIdx.x)/GRID;
    unsigned int numWorkUnits    = cSim.pInteractionCount[0];
    unsigned int pos             = warp*numWorkUnits/totalWarps;
    unsigned int end             = (warp+1)*numWorkUnits/totalWarps;
    unsigned int lasty           = 0xFFFFFFFF;
    float totalEnergy            = 0.0f;     
    float conversionFactor       = (cAmoebaSim.electric/cAmoebaSim.dielec);
    float scalingFactors[LastScalingIndex];

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

        ElectrostaticParticle* psA    = &sA[tbx];
        unsigned int atomI            = x + tgx;
        ElectrostaticParticle localParticle;
        loadElectrostaticShared(&localParticle, atomI,
                                atomCoord, labFrameDipole, labFrameQuadrupole,
                                inducedDipole, inducedDipolePolar, cAmoebaSim.pDampingFactorAndThole );

        localParticle.force[0]        = 0.0f;
        localParticle.force[1]        = 0.0f;
        localParticle.force[2]        = 0.0f;

        localParticle.torque[0]       = 0.0f;
        localParticle.torque[1]       = 0.0f;
        localParticle.torque[2]       = 0.0f;

        scalingFactors[PScaleIndex]   = 1.0f;
        scalingFactors[DScaleIndex]   = 1.0f;
        scalingFactors[UScaleIndex]   = 1.0f;
        scalingFactors[MScaleIndex]   = 1.0f;

        if (x == y) // Handle diagonals uniquely at 50% efficiency
        {

            // load shared data

           loadElectrostaticShared( &(sA[threadIdx.x]), atomI,
                                    atomCoord, labFrameDipole, labFrameQuadrupole,
                                    inducedDipole, inducedDipolePolar, cAmoebaSim.pDampingFactorAndThole );

            if (!bExclusionFlag)
            {

                // this branch is never exercised since it includes the
                // interaction between atomI and itself which is always excluded

                for (unsigned int j = 0; j < GRID; j++)
                {

                    float4 force;
                    float4 torque[2];

                    calculateElectrostaticPairIxn_kernel( localParticle, psA[j], scalingFactors, &force, torque
#ifdef AMOEBA_DEBUG
, pullBack
#endif
 );

                    if( (atomI != (y + j)) && (atomI < cAmoebaSim.numberOfAtoms) && ((y+j) < cAmoebaSim.numberOfAtoms) ){

                         localParticle.force[0]            += force.x;
                         localParticle.force[1]            += force.y;
                         localParticle.force[2]            += force.z;

                         totalEnergy                       += 0.5f*force.w;
     
                         localParticle.torque[0]           += torque[0].x;
                         localParticle.torque[1]           += torque[0].y;
                         localParticle.torque[2]           += torque[0].z;
      
                    }
                }

            }
            else  // bExclusion
            {
                unsigned int xi       = x >> GRIDBITS;
                unsigned int cell     = xi + xi*cAmoebaSim.paddedNumberOfAtoms/GRID-xi*(xi+1)/2;
                int  dScaleMask       = cAmoebaSim.pD_ScaleIndices[cAmoebaSim.pScaleIndicesIndex[cell]+tgx];
                int2 pScaleMask       = cAmoebaSim.pP_ScaleIndices[cAmoebaSim.pScaleIndicesIndex[cell]+tgx];
                int2 mScaleMask       = cAmoebaSim.pM_ScaleIndices[cAmoebaSim.pScaleIndicesIndex[cell]+tgx];

                for (unsigned int j = 0; j < GRID; j++)
                {

                    float4 force;
                    float4 torque[2];

                    unsigned int atomJ = y + j;

                    // set scale factors

                    getMaskedDScaleFactor( j, dScaleMask, scalingFactors + DScaleIndex );
                    getMaskedPScaleFactor( j, pScaleMask, scalingFactors + PScaleIndex );
                    getMaskedMScaleFactor( j, mScaleMask, scalingFactors + MScaleIndex );

                    // force

                    calculateElectrostaticPairIxn_kernel( localParticle, psA[j], scalingFactors, &force, torque
#ifdef AMOEBA_DEBUG
, pullBack
#endif
 );

                    // nan*0.0 = nan not 0.0, so explicitly exclude (atomI == atomJ) contribution
                    // by setting match flag

                    if( (atomI != atomJ) && (atomI < cAmoebaSim.numberOfAtoms) && (atomJ < cAmoebaSim.numberOfAtoms) ){

                        localParticle.force[0]            += force.x;
                        localParticle.force[1]            += force.y;
                        localParticle.force[2]            += force.z;

                        totalEnergy                       += 0.5f*force.w;

                        localParticle.torque[0]           += torque[0].x;
                        localParticle.torque[1]           += torque[0].y;
                        localParticle.torque[2]           += torque[0].z;
      

                    }
    
#ifdef AMOEBA_DEBUG
if( atomI == targetAtom ){
            unsigned int index                 = (atomI == targetAtom) ? atomJ : atomI;

            debugArray[index].x                = (float) atomI;
            debugArray[index].y                = (float) atomJ;
            debugArray[index].z                = 1.0f;
            debugArray[index].w                = (float) y;

            unsigned int mask                  = (atomI != atomJ) && (atomI < cAmoebaSim.numberOfAtoms) && (atomJ < cAmoebaSim.numberOfAtoms) ? 1 : 0;
            index                             += cAmoebaSim.paddedNumberOfAtoms;
            debugArray[index].x                = mask ? conversionFactor*force.x : 0.0f;
            debugArray[index].y                = mask ? conversionFactor*force.y : 0.0f;
            debugArray[index].z                = mask ? conversionFactor*force.z : 0.0f;
            debugArray[index].w                = -1.1f;


            index                             += cAmoebaSim.paddedNumberOfAtoms;
            debugArray[index].x                = mask ? conversionFactor*torque[0].x : 0.0f;
            debugArray[index].y                = mask ? conversionFactor*torque[0].y : 0.0f;
            debugArray[index].z                = mask ? conversionFactor*torque[0].z : 0.0f;
            debugArray[index].w                = -2.1f;


            index                             += cAmoebaSim.paddedNumberOfAtoms;
            debugArray[index].x                = mask ? conversionFactor*torque[1].x : 0.0f;
            debugArray[index].y                = mask ? conversionFactor*torque[1].y : 0.0f;
            debugArray[index].z                = mask ? conversionFactor*torque[1].z : 0.0f;
            debugArray[index].w                = -2.2f;


            index                             += cAmoebaSim.paddedNumberOfAtoms;
            int pullIndex                      = 0;
            debugArray[index].x                = pullBack[pullIndex].x;
            debugArray[index].y                = pullBack[pullIndex].y;
            debugArray[index].z                = pullBack[pullIndex].z;
            debugArray[index].w                = pullBack[pullIndex].w;

            index                             += cAmoebaSim.paddedNumberOfAtoms;
            pullIndex++;
            debugArray[index].x                = conversionFactor*pullBack[pullIndex].x;
            debugArray[index].y                = conversionFactor*pullBack[pullIndex].y;
            debugArray[index].z                = conversionFactor*pullBack[pullIndex].z;
            debugArray[index].w                = pullBack[pullIndex].w;

            index                             += cAmoebaSim.paddedNumberOfAtoms;
            pullIndex++;
            debugArray[index].x                = conversionFactor*pullBack[pullIndex].x;
            debugArray[index].y                = conversionFactor*pullBack[pullIndex].y;
            debugArray[index].z                = conversionFactor*pullBack[pullIndex].z;
            debugArray[index].w                = pullBack[pullIndex].w;

            index                             += cAmoebaSim.paddedNumberOfAtoms;
            pullIndex++;
            debugArray[index].x                = conversionFactor*pullBack[pullIndex].x;
            debugArray[index].y                = conversionFactor*pullBack[pullIndex].y;
            debugArray[index].z                = conversionFactor*pullBack[pullIndex].z;
            debugArray[index].w                = pullBack[pullIndex].w;

}
#endif

                }
            }

            // Write results

           localParticle.force[0]  *= conversionFactor;
           localParticle.force[1]  *= conversionFactor;
           localParticle.force[2]  *= conversionFactor;

           localParticle.torque[0] *= conversionFactor;
           localParticle.torque[1] *= conversionFactor;
           localParticle.torque[2] *= conversionFactor;

#ifdef USE_OUTPUT_BUFFER_PER_WARP
            float  of;
            unsigned int offset                 = 3*(x + tgx + warp*cAmoebaSim.paddedNumberOfAtoms);
            of                                  = outputForce[offset];
            of                                 += localParticle.force[0];
            outputForce[offset]                 = of;

            of                                  = outputForce[offset+1];
            of                                 += localParticle.force[1];
            outputForce[offset+1]               = of;

            of                                  = outputForce[offset+2];
            of                                 += localParticle.force[2];
            outputForce[offset+2]               = of;

            of                                  = outputTorque[offset];
            of                                 += localParticle.torque[0];
            outputTorque[offset]                = of;

            of                                  = outputTorque[offset+1];
            of                                 += localParticle.torque[1];
            outputTorque[offset+1]              = of;

            of                                  = outputTorque[offset+2];
            of                                 += localParticle.torque[2];
            outputTorque[offset+2]              = of;

#else
            unsigned int offset                 = 3*(x + tgx + (x >> GRIDBITS) * cAmoebaSim.paddedNumberOfAtoms);
            outputForce[offset]                 = localParticle.force[0];
            outputForce[offset+1]               = localParticle.force[1];
            outputForce[offset+2]               = localParticle.force[2];

            outputTorque[offset]                = localParticle.torque[0];
            outputTorque[offset+1]              = localParticle.torque[1];
            outputTorque[offset+2]              = localParticle.torque[2];
#endif

        }
        else        // 100% utilization
        {
            // Read fixed atom data into registers and GRF

            if (lasty != y)
            {
                // load shared data

               loadElectrostaticShared( &(sA[threadIdx.x]), (y+tgx),
                                        atomCoord, labFrameDipole, labFrameQuadrupole,
                                        inducedDipole, inducedDipolePolar, cAmoebaSim.pDampingFactorAndThole );

            }

            sA[threadIdx.x].force[0]     = 0.0f;
            sA[threadIdx.x].force[1]     = 0.0f;
            sA[threadIdx.x].force[2]     = 0.0f;

            sA[threadIdx.x].torque[0]    = 0.0f;
            sA[threadIdx.x].torque[1]    = 0.0f;
            sA[threadIdx.x].torque[2]    = 0.0f;

            if (!bExclusionFlag)
            {
                for (unsigned int j = 0; j < GRID; j++)
                {

                    float4 force;
                    float4 torque[2];

                    unsigned int atomJ = y + tj;

                    calculateElectrostaticPairIxn_kernel( localParticle, psA[tj], scalingFactors, &force, torque
#ifdef AMOEBA_DEBUG
, pullBack
#endif
 );
       
                    if( (atomI < cAmoebaSim.numberOfAtoms) && (atomJ < cAmoebaSim.numberOfAtoms) ){ 

                        // add force and torque to atom I due atom J
           
                        localParticle.force[0]         += force.x;
                        localParticle.force[1]         += force.y;
                        localParticle.force[2]         += force.z;

                        localParticle.torque[0]        += torque[0].x;
                        localParticle.torque[1]        += torque[0].y;
                        localParticle.torque[2]        += torque[0].z;
    
                        totalEnergy                    += force.w;
               
                        // add force and torque to atom J due atom I
           
                        psA[tj].force[0]               -= force.x;
                        psA[tj].force[1]               -= force.y;
                        psA[tj].force[2]               -= force.z;
           
                        psA[tj].torque[0]              += torque[1].x;
                        psA[tj].torque[1]              += torque[1].y;
                        psA[tj].torque[2]              += torque[1].z;
                    }
       
#ifdef AMOEBA_DEBUG
if( atomI == targetAtom  || atomJ == targetAtom ){
        unsigned int index                 = (atomI == targetAtom) ? atomJ : atomI;
        unsigned int indexI                = (atomI == targetAtom) ? 0 : 1;
        unsigned int indexJ                = (atomI == targetAtom) ? 1 : 0;
        float forceSign                    = (atomI == targetAtom) ? 1.0f : -1.0f;

        debugArray[index].x                = (float) atomI;
        debugArray[index].y                = (float) atomJ;
        debugArray[index].z                = 2.0f;
        debugArray[index].w                = (float) y;

        unsigned int mask                  = (atomI < cAmoebaSim.numberOfAtoms) && (atomJ < cAmoebaSim.numberOfAtoms) ? 1 : 0;
        index                             += cAmoebaSim.paddedNumberOfAtoms;
        debugArray[index].x                = mask ? conversionFactor*forceSign*force.x : 0.0f;
        debugArray[index].y                = mask ? conversionFactor*forceSign*force.y : 0.0f;
        debugArray[index].z                = mask ? conversionFactor*forceSign*force.z : 0.0f;
        debugArray[index].w                = -3.1f;

        index                             += cAmoebaSim.paddedNumberOfAtoms;
        debugArray[index].x                = mask ? conversionFactor*torque[indexI].x  : 0.0f;
        debugArray[index].y                = mask ? conversionFactor*torque[indexI].y  : 0.0f;
        debugArray[index].z                = mask ? conversionFactor*torque[indexI].z  : 0.0f;
        debugArray[index].w                = -4.1f;

        index                             += cAmoebaSim.paddedNumberOfAtoms;
        debugArray[index].x                = mask ? conversionFactor*torque[indexJ].x  : 0.0f;
        debugArray[index].y                = mask ? conversionFactor*torque[indexJ].y  : 0.0f;
        debugArray[index].z                = mask ? conversionFactor*torque[indexJ].z  : 0.0f;
        debugArray[index].w                = -5.1f;

        index                             += cAmoebaSim.paddedNumberOfAtoms;
        int pullIndex                      = 0;
        debugArray[index].x                = pullBack[pullIndex].x;
        debugArray[index].y                = pullBack[pullIndex].y;
        debugArray[index].z                = pullBack[pullIndex].z;
        debugArray[index].w                = pullBack[pullIndex].w;

        index                             += cAmoebaSim.paddedNumberOfAtoms;
        pullIndex++;
        debugArray[index].x                = conversionFactor*pullBack[pullIndex].x;
        debugArray[index].y                = conversionFactor*pullBack[pullIndex].y;
        debugArray[index].z                = conversionFactor*pullBack[pullIndex].z;
        debugArray[index].w                = conversionFactor*pullBack[pullIndex].w;

        index                             += cAmoebaSim.paddedNumberOfAtoms;
        pullIndex++;
        debugArray[index].x                = conversionFactor*pullBack[pullIndex].x;
        debugArray[index].y                = conversionFactor*pullBack[pullIndex].y;
        debugArray[index].z                = conversionFactor*pullBack[pullIndex].z;
        debugArray[index].w                = conversionFactor*pullBack[pullIndex].w;
       
        index                             += cAmoebaSim.paddedNumberOfAtoms;
        pullIndex++;
        debugArray[index].x                = pullBack[pullIndex].x;
        debugArray[index].y                = pullBack[pullIndex].y;
        debugArray[index].z                = pullBack[pullIndex].z;
        debugArray[index].w                = pullBack[pullIndex].w;
       
}
#endif
       
                           tj                  = (tj + 1) & (GRID - 1);
       
                       }
                   }
                   else  // bExclusion
                   {
                       // Read fixed atom data into registers and GRF
       
                       unsigned int xi   = x >> GRIDBITS;
                       unsigned int yi   = y >> GRIDBITS;
                       unsigned int cell = xi+yi*cAmoebaSim.paddedNumberOfAtoms/GRID-yi*(yi+1)/2;
                       int  dScaleMask   = cAmoebaSim.pD_ScaleIndices[cAmoebaSim.pScaleIndicesIndex[cell]+tgx];
                       int2 pScaleMask   = cAmoebaSim.pP_ScaleIndices[cAmoebaSim.pScaleIndicesIndex[cell]+tgx];
                       int2 mScaleMask   = cAmoebaSim.pM_ScaleIndices[cAmoebaSim.pScaleIndicesIndex[cell]+tgx];
       
                       for (unsigned int j = 0; j < GRID; j++)
                       {
       
                           float4 force;
                           float4 torque[2];
       
                           unsigned int atomJ = y + tj;
       
                           // set scale factors
       
                           getMaskedDScaleFactor( tj, dScaleMask, scalingFactors + DScaleIndex );
                           getMaskedPScaleFactor( tj, pScaleMask, scalingFactors + PScaleIndex );
                           getMaskedMScaleFactor( tj, mScaleMask, scalingFactors + MScaleIndex );
       
                           // force
       
                           calculateElectrostaticPairIxn_kernel( localParticle, psA[tj], scalingFactors, &force, torque
#ifdef AMOEBA_DEBUG
, pullBack
#endif
 );

                    // check if atoms out-of-bounds

                    if( (atomI < cAmoebaSim.numberOfAtoms) && (atomJ < cAmoebaSim.numberOfAtoms) ){

                        // add force and torque to atom I due atom J
    
                        localParticle.force[0]         += force.x;
                        localParticle.force[1]         += force.y;
                        localParticle.force[2]         += force.z;
    
                        localParticle.torque[0]        += torque[0].x;
                        localParticle.torque[1]        += torque[0].y;
                        localParticle.torque[2]        += torque[0].z;
        
                        totalEnergy                    += force.w;
    
                        // add force and torque to atom J due atom I
    
                        psA[tj].force[0]               -= force.x;
                        psA[tj].force[1]               -= force.y;
                        psA[tj].force[2]               -= force.z;
    
                        psA[tj].torque[0]              += torque[1].x;
                        psA[tj].torque[1]              += torque[1].y;
                        psA[tj].torque[2]              += torque[1].z;
    
                    }

#ifdef AMOEBA_DEBUG
if( atomI == targetAtom  || atomJ == targetAtom ){
        unsigned int index                 = (atomI == targetAtom) ? atomJ : atomI;
        unsigned int indexI                = (atomI == targetAtom) ? 0 : 1;
        unsigned int indexJ                = (atomI == targetAtom) ? 1 : 0;
        float forceSign                    = (atomI == targetAtom) ? 1.0f : -1.0f;

        debugArray[index].x                   = (float) atomI;
        debugArray[index].y                   = (float) atomJ;
        debugArray[index].z                   = 3.0f;
        debugArray[index].w                   = (float) y;

        unsigned int mask                     = (atomI < cAmoebaSim.numberOfAtoms) && (atomJ < cAmoebaSim.numberOfAtoms) ? 1 : 0;
        index                                += cAmoebaSim.paddedNumberOfAtoms;
        debugArray[index].x                   = mask ? conversionFactor*forceSign*force.x : 0.0f;
        debugArray[index].y                   = mask ? conversionFactor*forceSign*force.y : 0.0f;
        debugArray[index].z                   = mask ? conversionFactor*forceSign*force.z : 0.0f;
        debugArray[index].w                   = -6.1f;

        index                                += cAmoebaSim.paddedNumberOfAtoms;
        debugArray[index].x                   = mask ? conversionFactor*torque[indexI].x  : 0.0f;
        debugArray[index].y                   = mask ? conversionFactor*torque[indexI].y  : 0.0f;
        debugArray[index].z                   = mask ? conversionFactor*torque[indexI].z  : 0.0f;
        debugArray[index].w                   = -7.1f;

        index                                += cAmoebaSim.paddedNumberOfAtoms;
        debugArray[index].x                   = mask ? conversionFactor*torque[indexJ].x  : 0.0f;
        debugArray[index].y                   = mask ? conversionFactor*torque[indexJ].y  : 0.0f;
        debugArray[index].z                   = mask ? conversionFactor*torque[indexJ].z  : 0.0f;
        debugArray[index].w                   = -8.1f;

        index                                += cAmoebaSim.paddedNumberOfAtoms;
        int pullIndex                         = 0;
        debugArray[index].x                   = conversionFactor*pullBack[pullIndex].x;
        debugArray[index].y                   = conversionFactor*pullBack[pullIndex].y;
        debugArray[index].z                   = conversionFactor*pullBack[pullIndex].z;
        debugArray[index].w                   = conversionFactor*pullBack[pullIndex].w;

        index                                += cAmoebaSim.paddedNumberOfAtoms;
        pullIndex++;
        debugArray[index].x                   = conversionFactor*pullBack[pullIndex].x;
        debugArray[index].y                   = conversionFactor*pullBack[pullIndex].y;
        debugArray[index].z                   = conversionFactor*pullBack[pullIndex].z;
        debugArray[index].w                   = conversionFactor*pullBack[pullIndex].w;

        index                                += cAmoebaSim.paddedNumberOfAtoms;
        pullIndex++;
        debugArray[index].x                   = conversionFactor*pullBack[pullIndex].x;
        debugArray[index].y                   = conversionFactor*pullBack[pullIndex].y;
        debugArray[index].z                   = conversionFactor*pullBack[pullIndex].z;
        debugArray[index].w                   = conversionFactor*pullBack[pullIndex].w;
       
        index                                += cAmoebaSim.paddedNumberOfAtoms;
        pullIndex++;
        debugArray[index].x                   = pullBack[pullIndex].x;
        debugArray[index].y                   = pullBack[pullIndex].y;
        debugArray[index].z                   = pullBack[pullIndex].z;
        debugArray[index].w                   = pullBack[pullIndex].w;
       
#if 0
            index                                += cAmoebaSim.paddedNumberOfAtoms;
            debugArray[index].x                   = scalingFactors[DScaleIndex];
            debugArray[index].y                   = dScaleVal;
            debugArray[index].z                   = scalingFactors[PScaleIndex];
            debugArray[index].w                   = pScaleVal;

            index                                += cAmoebaSim.paddedNumberOfAtoms;
            debugArray[index].x                   = scalingFactors[MScaleIndex];
            debugArray[index].y                   = mScaleVal;
            for( int pIndex    = 0; pIndex < 14; pIndex++ ){
                index                                += cAmoebaSim.paddedNumberOfAtoms;
                debugArray[index].x                   = pullBack[pIndex].x;
                debugArray[index].y                   = pullBack[pIndex].y;
                debugArray[index].z                   = pullBack[pIndex].z;
                debugArray[index].w                   = pullBack[pIndex].w;
            }

            index                                += cAmoebaSim.paddedNumberOfAtoms;
            debugArray[index].x                   = labFrameDipole[3*atomI];
            debugArray[index].y                   = labFrameDipole[3*atomI+1];
            debugArray[index].z                   = labFrameDipole[3*atomI+2];
            debugArray[index].w                   = 25.0f;

            index                                += cAmoebaSim.paddedNumberOfAtoms;
            debugArray[index].x                   = labFrameDipole[3*atomJ];
            debugArray[index].y                   = labFrameDipole[3*atomJ+1];
            debugArray[index].z                   = labFrameDipole[3*atomJ+2];
            debugArray[index].w                   = 26.0f;

            index                                += cAmoebaSim.paddedNumberOfAtoms;
            debugArray[index].x                   = jDipole[0];
            debugArray[index].y                   = jDipole[1];
            debugArray[index].z                   = jDipole[2];
            debugArray[index].w                   = 27.0f;
#endif


}
#endif


                    tj                     = (tj + 1) & (GRID - 1);
                }
            }

            // Write results

           localParticle.force[0]     *= conversionFactor;
           localParticle.force[1]     *= conversionFactor;
           localParticle.force[2]     *= conversionFactor;

           localParticle.torque[0]    *= conversionFactor;
           localParticle.torque[1]    *= conversionFactor;
           localParticle.torque[2]    *= conversionFactor;

           sA[threadIdx.x].force[0]   *= conversionFactor;
           sA[threadIdx.x].force[1]   *= conversionFactor;
           sA[threadIdx.x].force[2]   *= conversionFactor;

           sA[threadIdx.x].torque[0]  *= conversionFactor;
           sA[threadIdx.x].torque[1]  *= conversionFactor;
           sA[threadIdx.x].torque[2]  *= conversionFactor;

#ifdef USE_OUTPUT_BUFFER_PER_WARP

            float of;
            unsigned int offset                 = 3*(x + tgx + warp*cAmoebaSim.paddedNumberOfAtoms);
            of                                  = outputForce[offset];
            of                                 += localParticle.force[0];
            outputForce[offset]                 = of;

            of                                  = outputForce[offset+1];
            of                                 += localParticle.force[1];
            outputForce[offset+1]               = of;

            of                                  = outputForce[offset+2];
            of                                 += localParticle.force[2];
            outputForce[offset+2]               = of;

            of                                  = outputTorque[offset];
            of                                 += localParticle.torque[0];
            outputTorque[offset]                 = of;

            of                                  = outputTorque[offset+1];
            of                                 += localParticle.torque[1];
            outputTorque[offset+1]               = of;

            of                                  = outputTorque[offset+2];
            of                                 += localParticle.torque[2];
            outputTorque[offset+2]               = of;

            offset                              = 3*(y + tgx + warp*cAmoebaSim.paddedNumberOfAtoms);

            of                                  = outputForce[offset];
            of                                 += sA[threadIdx.x].force[0];
            outputForce[offset]                 = of;

            of                                  = outputForce[offset+1];
            of                                 += sA[threadIdx.x].force[1];
            outputForce[offset+1]               = of;

            of                                  = outputForce[offset+2];
            of                                 += sA[threadIdx.x].force[2];
            outputForce[offset+2]               = of;

            of                                  = outputTorque[offset];
            of                                 += sA[threadIdx.x].torque[0];
            outputTorque[offset]                = of;

            of                                  = outputTorque[offset+1];
            of                                 += sA[threadIdx.x].torque[1];
            outputTorque[offset+1]              = of;

            of                                  = outputTorque[offset+2];
            of                                 += sA[threadIdx.x].torque[2];
            outputTorque[offset+2]              = of;

#else
            unsigned int offset                 = 3*(x + tgx + (y >> GRIDBITS) * cAmoebaSim.paddedNumberOfAtoms);

            outputForce[offset]                 = localParticle.force[0];
            outputForce[offset+1]               = localParticle.force[1];
            outputForce[offset+2]               = localParticle.force[2];

            outputTorque[offset]                = localParticle.torque[0];
            outputTorque[offset+1]              = localParticle.torque[1];
            outputTorque[offset+2]              = localParticle.torque[2];

            offset                              = 3*(y + tgx + (x >> GRIDBITS) * cAmoebaSim.paddedNumberOfAtoms);

            outputForce[offset]                 = sA[threadIdx.x].force[0];
            outputForce[offset+1]               = sA[threadIdx.x].force[1];
            outputForce[offset+2]               = sA[threadIdx.x].force[2];

            outputTorque[offset]                = sA[threadIdx.x].torque[0];
            outputTorque[offset+1]              = sA[threadIdx.x].torque[1];
            outputTorque[offset+2]              = sA[threadIdx.x].torque[2];


#endif
            lasty = y;
        }

        pos++;
    }
    cSim.pEnergy[blockIdx.x * blockDim.x + threadIdx.x] += (conversionFactor*totalEnergy);
}
