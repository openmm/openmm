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

#include "cudaKernels.h"
#include "amoebaCudaKernels.h"
#include "kCalculateAmoebaCudaUtilities.h"
#include "openmm/OpenMMException.h"

#include <stdio.h>
#include <cuda.h>
#include <cstdlib>
using namespace std; 

#define SQRT sqrtf

static __constant__ cudaGmxSimulation cSim;
static __constant__ cudaAmoebaGmxSimulation cAmoebaSim;
extern __global__ void kFindInteractionsWithinBlocksPeriodic_kernel(unsigned int*);

void SetCalculateAmoebaMultipoleForcesSim(amoebaGpuContext amoebaGpu)
{
    cudaError_t status;
    gpuContext gpu = amoebaGpu->gpuContext;
    status         = cudaMemcpyToSymbol(cSim, &gpu->sim, sizeof(cudaGmxSimulation));     
    RTERROR(status, "SetCalculateAmoebaMultipoleForcesSim: cudaMemcpyToSymbol: SetSim copy to cSim failed");
    status         = cudaMemcpyToSymbol(cAmoebaSim, &amoebaGpu->amoebaSim, sizeof(cudaAmoebaGmxSimulation));     
    RTERROR(status, "SetCalculateAmoebaMultipoleForcesSim: cudaMemcpyToSymbol: SetSim copy to cAmoebaSim failed");
}

void GetCalculateAmoebaMultipoleForcesSim(amoebaGpuContext amoebaGpu)
{
    cudaError_t status;
    gpuContext gpu = amoebaGpu->gpuContext;
    status = cudaMemcpyFromSymbol(&gpu->sim, cSim, sizeof(cudaGmxSimulation));     
    RTERROR(status, "GetCalculateAmoebaMultipoleForcesSim: cudaMemcpyFromSymbol: SetSim copy from cSim failed");
    status = cudaMemcpyFromSymbol(&amoebaGpu->amoebaSim, cAmoebaSim, sizeof(cudaAmoebaGmxSimulation));     
    RTERROR(status, "GetCalculateAmoebaMultipoleForcesSim: cudaMemcpyFromSymbol: SetSim copy from cAmoebaSim failed");
}

__device__ static float normVector3( float* vector )
{

    float norm                    = DOT3( vector, vector );
    float returnNorm              = SQRT( norm );
    norm                          = returnNorm > 0.0f ? 1.0f/returnNorm : 0.0f;

    vector[0]                    *= norm;
    vector[1]                    *= norm;
    vector[2]                    *= norm;

    return returnNorm;
}

// ZThenX     == 0
// Bisector   == 1
// ZBisect    == 2
// ThreeFold  == 3
// ZOnly      == 4
// NoAxisType == 5

__global__
#if (__CUDA_ARCH__ >= 200)
__launch_bounds__(GF1XX_THREADS_PER_BLOCK, 1)
#elif (__CUDA_ARCH__ >= 120)
__launch_bounds__(GT2XX_THREADS_PER_BLOCK, 1)
#else
__launch_bounds__(G8X_THREADS_PER_BLOCK, 1)
#endif
void kCudaComputeCheckChiral_kernel( void )
{

    const int AD          = 0;
    const int BD          = 1;
    const int CD          = 2;
    const int C           = 3;
    float delta[4][3];
 
    float4* particleCoord        = cSim.pPosq;
    int4* multiPoleParticles     = cAmoebaSim.pMultipoleParticlesIdsAndAxisType;
    float* labFrameDipole        = cAmoebaSim.pLabFrameDipole;
    float* labFrameQuadrupole    = cAmoebaSim.pLabFrameQuadrupole;
 
    // ---------------------------------------------------------------------------------------
 
    int particleIndex            = blockIdx.x*blockDim.x + threadIdx.x;
    int numberOfParticles        = cSim.atoms;
    while( particleIndex < numberOfParticles )
    { 
        // skip z-then-x
    
        int axisType             = multiPoleParticles[particleIndex].w; 
        if( axisType != 0 && multiPoleParticles[particleIndex].x >= 0 && multiPoleParticles[particleIndex].y >=0 && multiPoleParticles[particleIndex].z >= 0 )
        {
     
            // ---------------------------------------------------------------------------------------
         
            int particleA                = particleIndex;
            int particleB                = multiPoleParticles[particleIndex].z;
            int particleC                = multiPoleParticles[particleIndex].x;
            int particleD                = multiPoleParticles[particleIndex].y;
        
            delta[AD][0]                 = particleCoord[particleA].x - particleCoord[particleD].x;
            delta[AD][1]                 = particleCoord[particleA].y - particleCoord[particleD].y;
            delta[AD][2]                 = particleCoord[particleA].z - particleCoord[particleD].z;
        
            delta[BD][0]                 = particleCoord[particleB].x - particleCoord[particleD].x;
            delta[BD][1]                 = particleCoord[particleB].y - particleCoord[particleD].y;
            delta[BD][2]                 = particleCoord[particleB].z - particleCoord[particleD].z;
        
            delta[CD][0]                 = particleCoord[particleC].x - particleCoord[particleD].x;
            delta[CD][1]                 = particleCoord[particleC].y - particleCoord[particleD].y;
            delta[CD][2]                 = particleCoord[particleC].z - particleCoord[particleD].z;
        
            delta[C][0]                  = delta[BD][1]*delta[CD][2] - delta[BD][2]*delta[CD][1];
            delta[C][1]                  = delta[CD][1]*delta[AD][2] - delta[CD][2]*delta[AD][1];
            delta[C][2]                  = delta[AD][1]*delta[BD][2] - delta[AD][2]*delta[BD][1];
         
            float volume                 = delta[C][0]*delta[AD][0] + delta[C][1]*delta[BD][0] + delta[C][2]*delta[CD][0];
            if( volume < 0.0 ){
                labFrameDipole[particleIndex*3+1]            *= -1.0f; // pole(3,i)
                labFrameQuadrupole[particleIndex*9+1]        *= -1.0f; // pole(6,i)  && pole(8,i)
                labFrameQuadrupole[particleIndex*9+3]        *= -1.0f; // pole(10,i) && pole(12,i)
                labFrameQuadrupole[particleIndex*9+5]        *= -1.0f; // pole(6,i)  && pole(8,i)
                labFrameQuadrupole[particleIndex*9+7]        *= -1.0f; // pole(10,i) && pole(12,i)
            }
        }
    
        particleIndex                += gridDim.x*blockDim.x;
    }    
}

__global__
#if (__CUDA_ARCH__ >= 200)
__launch_bounds__(GF1XX_THREADS_PER_BLOCK, 1)
#elif (__CUDA_ARCH__ >= 120)
__launch_bounds__(GT2XX_THREADS_PER_BLOCK, 1)
#else
__launch_bounds__(G8X_THREADS_PER_BLOCK, 1)
#endif
void kCudaComputeLabFrameMoments_kernel( void )
{

    float vectorX[3];
    float vectorY[3];
    float vectorZ[3];
 
    int particleIndex            = blockIdx.x*blockDim.x + threadIdx.x;

    float4* particleCoord        = cSim.pPosq;
    int4* multiPoleParticles     = cAmoebaSim.pMultipoleParticlesIdsAndAxisType;
    float* labFrameDipole        = cAmoebaSim.pLabFrameDipole;
    float* labFrameQuadrupole    = cAmoebaSim.pLabFrameQuadrupole;
 
    // get coordinates of this atom and the z & x axis atoms
    // compute the vector between the atoms and 1/sqrt(d2), d2 is distance between
    // this atom and the axis atom
 
    // this atom is referred to as the k-atom in notes below
 
    // code common to ZThenX and Bisector
    
    while( particleIndex < cSim.atoms )
    {

        if( multiPoleParticles[particleIndex].x >= 0 && multiPoleParticles[particleIndex].z >= 0 )
        {
            float4 coordinatesThisParticle   = particleCoord[particleIndex];
         
            int multipoleParticleIndex       = multiPoleParticles[particleIndex].z;
            float4 coordinatesAxisParticle   = particleCoord[multipoleParticleIndex];
         
            vectorZ[0]                       = coordinatesAxisParticle.x - coordinatesThisParticle.x;
            vectorZ[1]                       = coordinatesAxisParticle.y - coordinatesThisParticle.y;
            vectorZ[2]                       = coordinatesAxisParticle.z - coordinatesThisParticle.z;
              
            multipoleParticleIndex           = multiPoleParticles[particleIndex].x; 
            coordinatesAxisParticle          = particleCoord[multipoleParticleIndex];
         
            vectorX[0]                       = coordinatesAxisParticle.x - coordinatesThisParticle.x;
            vectorX[1]                       = coordinatesAxisParticle.y - coordinatesThisParticle.y;
            vectorX[2]                       = coordinatesAxisParticle.z - coordinatesThisParticle.z;
         
            int axisType                     = multiPoleParticles[particleIndex].w; 
    
            /*
                z-only
                   (1) norm z
                   (2) select random x
                   (3) x = x - (x.z)z
                   (4) norm x
        
                z-then-x
                   (1) norm z
                   (2) norm x (not needed)
                   (3) x = x - (x.z)z
                   (4) norm x
        
                bisector
                   (1) norm z
                   (2) norm x 
                   (3) z = x + z
                   (4) norm z
                   (5) x = x - (x.z)z 
                   (6) norm x 
        
                z-bisect
                   (1) norm z
                   (2) norm x 
                   (3) norm y 
                   (3) x = x + y
                   (4) norm x
                   (5) x = x - (x.z)z 
                   (6) norm x 
        
                3-fold
                   (1) norm z
                   (2) norm x 
                   (3) norm y 
                   (4) z = x + y + z
                   (5) norm z
                   (6) x = x - (x.z)z 
                   (7) norm x 
        
            */
        
            // branch based on axis type
             
            float sum                   = normVector3( vectorZ );
        
            if( axisType == 1 ){
        
                // bisector
                
                sum                     = normVector3( vectorX );
                
                vectorZ[0]             += vectorX[0];
                vectorZ[1]             += vectorX[1];
                vectorZ[2]             += vectorX[2];
           
                sum                     = normVector3( vectorZ );
        
            } else if( axisType == 2 || axisType == 3 ){ 
         
                // z-bisect
        
                multipoleParticleIndex  = multiPoleParticles[particleIndex].y; 
                if( multipoleParticleIndex >= 0 && multipoleParticleIndex < cSim.atoms ){
                    coordinatesAxisParticle = particleCoord[multipoleParticleIndex];
                    vectorY[0]              = coordinatesAxisParticle.x - coordinatesThisParticle.x;
                    vectorY[1]              = coordinatesAxisParticle.y - coordinatesThisParticle.y;
                    vectorY[2]              = coordinatesAxisParticle.z - coordinatesThisParticle.z;
            
                    sum                     = normVector3( vectorY );
                    sum                     = normVector3( vectorX );
            
                    if( axisType == 2 ){
            
                        vectorX[0]         += vectorY[0];
                        vectorX[1]         += vectorY[1];
                        vectorX[2]         += vectorY[2];
                        sum                 = normVector3( vectorX );
             
                    } else { 
             
                        // 3-fold
                
                        vectorZ[0]         += vectorX[0] + vectorY[0];
                        vectorZ[1]         += vectorX[1] + vectorY[1];
                        vectorZ[2]         += vectorX[2] + vectorY[2];
                        sum                 = normVector3( vectorZ );
                    }
                }
         
            } else if( axisType >= 4 ){ 
        
                vectorX[0]             = 0.1f;
                vectorX[1]             = 0.1f;
                vectorX[2]             = 0.1f;
            }
            
            // x = x - (x.z)z
        
            float dot         = vectorZ[0]*vectorX[0] + vectorZ[1]*vectorX[1] + vectorZ[2]*vectorX[2];
                
            vectorX[0]       -= dot*vectorZ[0];
            vectorX[1]       -= dot*vectorZ[1];
            vectorX[2]       -= dot*vectorZ[2];
             
            sum               = normVector3( vectorX );
        
            vectorY[0]        = (vectorZ[1]*vectorX[2]) - (vectorZ[2]*vectorX[1]);
            vectorY[1]        = (vectorZ[2]*vectorX[0]) - (vectorZ[0]*vectorX[2]);
            vectorY[2]        = (vectorZ[0]*vectorX[1]) - (vectorZ[1]*vectorX[0]);
         
            // use identity rotation matrix for unrecognized axis types
        
            if( axisType < 0 || axisType > 4 ){
        
                vectorX[0] = 1.0f;
                vectorX[1] = 0.0f;
                vectorX[2] = 0.0f;
        
                vectorY[0] = 0.0f;
                vectorY[1] = 1.0f;
                vectorY[2] = 0.0f;
        
                vectorZ[0] = 0.0f;
                vectorZ[1] = 0.0f;
                vectorZ[2] = 1.0f;
            }
        
            unsigned int offset           = 3*particleIndex;
    
            float molDipole[3];
            molDipole[0]                  = labFrameDipole[offset];
            molDipole[1]                  = labFrameDipole[offset+1];
            molDipole[2]                  = labFrameDipole[offset+2];
            
            // set out-of-range elements to 0.0f
         
            labFrameDipole[offset]        = molDipole[0]*vectorX[0] + molDipole[1]*vectorY[0] + molDipole[2]*vectorZ[0];
            labFrameDipole[offset+1]      = molDipole[0]*vectorX[1] + molDipole[1]*vectorY[1] + molDipole[2]*vectorZ[1];
            labFrameDipole[offset+2]      = molDipole[0]*vectorX[2] + molDipole[1]*vectorY[2] + molDipole[2]*vectorZ[2];
            
            // ---------------------------------------------------------------------------------------
            
            float mPole[3][3];
            offset                        = 9*particleIndex;
            
            mPole[0][0]                   = labFrameQuadrupole[offset];
            mPole[0][1]                   = labFrameQuadrupole[offset+1];
            mPole[0][2]                   = labFrameQuadrupole[offset+2];
        
            mPole[1][0]                   = labFrameQuadrupole[offset+3];
            mPole[1][1]                   = labFrameQuadrupole[offset+4];
            mPole[1][2]                   = labFrameQuadrupole[offset+5];
        
            mPole[2][0]                   = labFrameQuadrupole[offset+6];
            mPole[2][1]                   = labFrameQuadrupole[offset+7];
            mPole[2][2]                   = labFrameQuadrupole[offset+8];
        
            labFrameQuadrupole[offset+8]  = vectorX[2]*(vectorX[2]*mPole[0][0] + vectorY[2]*mPole[0][1] + vectorZ[2]*mPole[0][2]);
            labFrameQuadrupole[offset+8] += vectorY[2]*(vectorX[2]*mPole[1][0] + vectorY[2]*mPole[1][1] + vectorZ[2]*mPole[1][2]);
            labFrameQuadrupole[offset+8] += vectorZ[2]*(vectorX[2]*mPole[2][0] + vectorY[2]*mPole[2][1] + vectorZ[2]*mPole[2][2]);
    
            labFrameQuadrupole[offset+4]  = vectorX[1]*(vectorX[1]*mPole[0][0] + vectorY[1]*mPole[0][1] + vectorZ[1]*mPole[0][2]);
            labFrameQuadrupole[offset+4] += vectorY[1]*(vectorX[1]*mPole[1][0] + vectorY[1]*mPole[1][1] + vectorZ[1]*mPole[1][2]);
            labFrameQuadrupole[offset+4] += vectorZ[1]*(vectorX[1]*mPole[2][0] + vectorY[1]*mPole[2][1] + vectorZ[1]*mPole[2][2]);
    
            labFrameQuadrupole[offset+5]  = vectorX[1]*(vectorX[2]*mPole[0][0] + vectorY[2]*mPole[0][1] + vectorZ[2]*mPole[0][2]);
            labFrameQuadrupole[offset+5] += vectorY[1]*(vectorX[2]*mPole[1][0] + vectorY[2]*mPole[1][1] + vectorZ[2]*mPole[1][2]);
            labFrameQuadrupole[offset+5] += vectorZ[1]*(vectorX[2]*mPole[2][0] + vectorY[2]*mPole[2][1] + vectorZ[2]*mPole[2][2]);
    
            labFrameQuadrupole[offset]    = vectorX[0]*(vectorX[0]*mPole[0][0] + vectorY[0]*mPole[0][1] + vectorZ[0]*mPole[0][2]);
            labFrameQuadrupole[offset]   += vectorY[0]*(vectorX[0]*mPole[1][0] + vectorY[0]*mPole[1][1] + vectorZ[0]*mPole[1][2]);
            labFrameQuadrupole[offset]   += vectorZ[0]*(vectorX[0]*mPole[2][0] + vectorY[0]*mPole[2][1] + vectorZ[0]*mPole[2][2]);
    
            labFrameQuadrupole[offset+1]  = vectorX[0]*(vectorX[1]*mPole[0][0] + vectorY[1]*mPole[0][1] + vectorZ[1]*mPole[0][2]);
            labFrameQuadrupole[offset+1] += vectorY[0]*(vectorX[1]*mPole[1][0] + vectorY[1]*mPole[1][1] + vectorZ[1]*mPole[1][2]);
            labFrameQuadrupole[offset+1] += vectorZ[0]*(vectorX[1]*mPole[2][0] + vectorY[1]*mPole[2][1] + vectorZ[1]*mPole[2][2]);
    
            labFrameQuadrupole[offset+2]  = vectorX[0]*(vectorX[2]*mPole[0][0] + vectorY[2]*mPole[0][1] + vectorZ[2]*mPole[0][2]);
            labFrameQuadrupole[offset+2] += vectorY[0]*(vectorX[2]*mPole[1][0] + vectorY[2]*mPole[1][1] + vectorZ[2]*mPole[1][2]);
            labFrameQuadrupole[offset+2] += vectorZ[0]*(vectorX[2]*mPole[2][0] + vectorY[2]*mPole[2][1] + vectorZ[2]*mPole[2][2]);
     
            labFrameQuadrupole[offset+3]  = labFrameQuadrupole[offset+1];
            labFrameQuadrupole[offset+6]  = labFrameQuadrupole[offset+2];
            labFrameQuadrupole[offset+7]  = labFrameQuadrupole[offset+5];
    
        }

        particleIndex                += gridDim.x*blockDim.x;
    }
     
}

void cudaComputeAmoebaLabFrameMoments( amoebaGpuContext amoebaGpu )
{

   // ---------------------------------------------------------------------------------------

    gpuContext gpu    = amoebaGpu->gpuContext;

    int numBlocks     = gpu->sim.blocks;
    int numThreads    = gpu->sim.threads_per_block;

    // copy molecular moments to lab frame moment arrays
    // check if chiral center requires moments to have sign flipped
    // compute lab frame moments

    cudaMemcpy( amoebaGpu->psLabFrameDipole->_pDevData,      amoebaGpu->psMolecularDipole->_pDevData,    3*gpu->sim.paddedNumberOfAtoms*sizeof( float ), cudaMemcpyDeviceToDevice ); 
    cudaMemcpy( amoebaGpu->psLabFrameQuadrupole->_pDevData, amoebaGpu->psMolecularQuadrupole->_pDevData, 9*gpu->sim.paddedNumberOfAtoms*sizeof( float ), cudaMemcpyDeviceToDevice ); 

    kCudaComputeCheckChiral_kernel<<< numBlocks, numThreads>>> ( );
    LAUNCHERROR("kCudaComputeCheckChiral");

    kCudaComputeLabFrameMoments_kernel<<< numBlocks, numThreads>>> ( );
    LAUNCHERROR("kCudaComputeLabFrameMoments");

}

void kCalculateAmoebaMultipoleForces(amoebaGpuContext amoebaGpu, bool hasAmoebaGeneralizedKirkwood ) 
{
    std::string methodName = "kCalculateAmoebaMultipoleForces";

    // compute lab frame moments

    cudaComputeAmoebaLabFrameMoments( amoebaGpu );

    if( 0 ){
        gpuContext gpu                       = amoebaGpu->gpuContext;
        std::vector<int> fileId;
        //fileId.push_back( 0 );
        VectorOfDoubleVectors outputVector;
        //cudaLoadCudaFloat4Array( gpu->natoms, 3, gpu->psPosq4,              outputVector, gpu->psAtomIndex->_pSysData, 1.0f );
        cudaLoadCudaFloatArray( gpu->natoms,  3, amoebaGpu->psLabFrameDipole,     outputVector, gpu->psAtomIndex->_pSysData, 1.0f );
        cudaLoadCudaFloatArray( gpu->natoms,  9, amoebaGpu->psLabFrameQuadrupole, outputVector, gpu->psAtomIndex->_pSysData, 1.0f );
        cudaWriteVectorOfDoubleVectorsToFile( "CudaLabMoments", fileId, outputVector );
    }   

    // compute fixed E-field and mutual induced field 

    if( hasAmoebaGeneralizedKirkwood ){
        cudaComputeAmoebaFixedEAndGkFields( amoebaGpu );
        cudaComputeAmoebaMutualInducedAndGkField( amoebaGpu );
    } else {

        if( amoebaGpu->multipoleNonbondedMethod == AMOEBA_NO_CUTOFF ){

            cudaComputeAmoebaFixedEField( amoebaGpu );
            cudaComputeAmoebaMutualInducedField( amoebaGpu );

        } else {

            gpuContext gpu = amoebaGpu->gpuContext;
            kFindBlockBoundsPeriodic_kernel<<<(gpu->psGridBoundingBox->_length+63)/64, 64>>>();
            LAUNCHERROR("kFindBlockBoundsPeriodic");
            kFindBlocksWithInteractionsPeriodic_kernel<<<gpu->sim.interaction_blocks, gpu->sim.interaction_threads_per_block>>>();
            LAUNCHERROR("kFindBlocksWithInteractionsPeriodic");

            compactStream(gpu->compactPlan, gpu->sim.pInteractingWorkUnit, gpu->sim.pWorkUnit, gpu->sim.pInteractionFlag, gpu->sim.workUnits, gpu->sim.pInteractionCount);

            //compactStream( gpu->compactPlan, 
            //               gpu->sim.pInteractingWorkUnit, unsigned int* dOut
            //               amoebaGpu->psWorkUnit->_pDevData, const unsigned int* dIn
            //               gpu->sim.pInteractionFlag,        const unsigned int* dValid
            //               gpu->sim.workUnits,               gpu
            //               gpu->sim.pInteractionCount);
            kFindInteractionsWithinBlocksPeriodic_kernel<<<gpu->sim.nonbond_blocks, gpu->sim.nonbond_threads_per_block,
                    sizeof(unsigned int)*gpu->sim.nonbond_threads_per_block>>>(gpu->sim.pInteractingWorkUnit);
            LAUNCHERROR("kFindInteractionsWithinBlocksPeriodic");

            cudaComputeAmoebaPmeFixedEField( amoebaGpu );
            cudaComputeAmoebaPmeMutualInducedField( amoebaGpu );
        }
    }

    // check if induce dipole calculation converged -- abort if it did not

    if( amoebaGpu->mutualInducedDone == 0 ){
       throw OpenMM::OpenMMException("Induced dipole calculation did not converge" );
    }

    // calculate electrostatic forces

    if( amoebaGpu->multipoleNonbondedMethod == AMOEBA_NO_CUTOFF ){
        cudaComputeAmoebaElectrostatic( amoebaGpu, (hasAmoebaGeneralizedKirkwood ? 0 : 1) );
    } else {
        cudaComputeAmoebaPmeElectrostatic( amoebaGpu );
    }

}
