//-----------------------------------------------------------------------------------------

//-----------------------------------------------------------------------------------------

#include "cudaKernels.h"
#include "amoebaCudaKernels.h"
#include "kCalculateAmoebaCudaUtilities.h"

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

#undef AMOEBA_DEBUG

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
 
    float4* atomCoord            = cSim.pPosq;
    int4* multiPoleAtoms         = cAmoebaSim.pMultipoleParticlesIdsAndAxisType;
    float* molecularDipole       = cAmoebaSim.pMolecularDipole;
    float* molecularQuadrupole   = cAmoebaSim.pMolecularQuadrupole;
    float* labFrameDipole        = cAmoebaSim.pLabFrameDipole;
    float* labFrameQuadrupole    = cAmoebaSim.pLabFrameQuadrupole;
 
    // ---------------------------------------------------------------------------------------
 
    int atomIndex                = blockIdx.x;
    if( atomIndex >= cSim.atoms )return;
 
    int axisType                 = multiPoleAtoms[atomIndex].w; 
 
    float* molDipole             = &(molecularDipole[atomIndex*3]);
    float* labDipole             = &(labFrameDipole[atomIndex*3]);
    labDipole[0]                 = molDipole[0];
    labDipole[1]                 = molDipole[1];
    labDipole[2]                 = molDipole[2];
 
    float* molQuadrupole         = &(molecularQuadrupole[atomIndex*9]);
    float* labQuadrupole         = &(labFrameQuadrupole[atomIndex*9]);
    labQuadrupole[0]             = molQuadrupole[0];
    labQuadrupole[1]             = molQuadrupole[1];
    labQuadrupole[2]             = molQuadrupole[2];
    labQuadrupole[3]             = molQuadrupole[3];
    labQuadrupole[4]             = molQuadrupole[4];
    labQuadrupole[5]             = molQuadrupole[5];
    labQuadrupole[6]             = molQuadrupole[6];
    labQuadrupole[7]             = molQuadrupole[7];
    labQuadrupole[8]             = molQuadrupole[8];

    // skip z-then-x

    if( axisType == 0 || multiPoleAtoms[atomIndex].y < 0 )return;
 
    // ---------------------------------------------------------------------------------------
 
    int atomA                    = atomIndex;
    int atomB                    = multiPoleAtoms[atomIndex].z;
    int atomC                    = multiPoleAtoms[atomIndex].x;
    int atomD                    = multiPoleAtoms[atomIndex].y;

    delta[AD][0]                 = atomCoord[atomA].x - atomCoord[atomD].x;
    delta[AD][1]                 = atomCoord[atomA].y - atomCoord[atomD].y;
    delta[AD][2]                 = atomCoord[atomA].z - atomCoord[atomD].z;

    delta[BD][0]                 = atomCoord[atomB].x - atomCoord[atomD].x;
    delta[BD][1]                 = atomCoord[atomB].y - atomCoord[atomD].y;
    delta[BD][2]                 = atomCoord[atomB].z - atomCoord[atomD].z;

    delta[CD][0]                 = atomCoord[atomC].x - atomCoord[atomD].x;
    delta[CD][1]                 = atomCoord[atomC].y - atomCoord[atomD].y;
    delta[CD][2]                 = atomCoord[atomC].z - atomCoord[atomD].z;

    delta[C][0]                  = delta[BD][1]*delta[CD][2] - delta[BD][2]*delta[CD][1];
    delta[C][1]                  = delta[CD][1]*delta[AD][2] - delta[CD][2]*delta[AD][1];
    delta[C][2]                  = delta[AD][1]*delta[BD][2] - delta[AD][2]*delta[BD][1];
 
    float volume                 = delta[C][0]*delta[AD][0] + delta[C][1]*delta[BD][0] + delta[C][2]*delta[CD][0];
    if( volume < 0.0 ){
        labDipole[1]            *= -1.0f; // pole(3,i)
        labQuadrupole[1]        *= -1.0f; // pole(6,i)  && pole(8,i)
        labQuadrupole[3]        *= -1.0f; // pole(10,i) && pole(12,i)
        labQuadrupole[5]        *= -1.0f; // pole(6,i)  && pole(8,i)
        labQuadrupole[7]        *= -1.0f; // pole(10,i) && pole(12,i)
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
 
    int numOfAtoms               = cSim.atoms;
    //float* rotationMatrix        = cAmoebaSim.pRotationMatrix;
    float4* atomCoord            = cSim.pPosq;
    int4* multiPoleAtoms         = cAmoebaSim.pMultipoleParticlesIdsAndAxisType;
    float* labFrameDipole        = cAmoebaSim.pLabFrameDipole;
    float* labFrameQuadrupole    = cAmoebaSim.pLabFrameQuadrupole;
 
    // ---------------------------------------------------------------------------------------
 
    int atomIndex = blockIdx.x;
 
    // ---------------------------------------------------------------------------------------
 
    // get coordinates of this atom and the z & x axis atoms
    // compute the vector between the atoms and 1/sqrt(d2), d2 is distance between
    // this atom and the axis atom
 
    // this atom is referred to as the k-atom in notes below
 
    // code common to ZThenX and Bisector
    
 /*
    vectorX                          = &(rotationMatrix[atomIndex*9]);
    vectorY                          = &(rotationMatrix[atomIndex*9+ 3]);
    vectorZ                          = &(rotationMatrix[atomIndex*9+ 6]);
 */
 
    float4 coordinatesThisAtom       = atomCoord[atomIndex];
 
    int multipoleAtomIndex           = multiPoleAtoms[atomIndex].z;
    float4 coordinatesAxisAtom       = atomCoord[multipoleAtomIndex];
 
    vectorZ[0]                       = coordinatesAxisAtom.x - coordinatesThisAtom.x;
    vectorZ[1]                       = coordinatesAxisAtom.y - coordinatesThisAtom.y;
    vectorZ[2]                       = coordinatesAxisAtom.z - coordinatesThisAtom.z;
      
    multipoleAtomIndex               = multiPoleAtoms[atomIndex].x; 
    coordinatesAxisAtom              = atomCoord[multipoleAtomIndex];
 
    vectorX[0]                       = coordinatesAxisAtom.x - coordinatesThisAtom.x;
    vectorX[1]                       = coordinatesAxisAtom.y - coordinatesThisAtom.y;
    vectorX[2]                       = coordinatesAxisAtom.z - coordinatesThisAtom.z;
 
    int axisType                     = multiPoleAtoms[atomIndex].w; 
      
    
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
     
    float sum                        = normVector3( vectorZ );

    if( axisType == 1 ){

        // bisector
        
        sum                     = normVector3( vectorX );
        
        vectorZ[0]             += vectorX[0];
        vectorZ[1]             += vectorX[1];
        vectorZ[2]             += vectorX[2];
   
        sum                     = normVector3( vectorZ );

    } else if( axisType == 2 || axisType == 3 ){ 
 
        // z-bisect

        multipoleAtomIndex      = multiPoleAtoms[atomIndex].y; 
        coordinatesAxisAtom     = atomCoord[multipoleAtomIndex];
        vectorY[0]              = coordinatesAxisAtom.x - coordinatesThisAtom.x;
        vectorY[1]              = coordinatesAxisAtom.y - coordinatesThisAtom.y;
        vectorY[2]              = coordinatesAxisAtom.z - coordinatesThisAtom.z;

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

    float molDipole[3];
    float* labDipole  = &(labFrameDipole[atomIndex*3]);
    molDipole[0]      = labDipole[0];
    molDipole[1]      = labDipole[1];
    molDipole[2]      = labDipole[2];
    
    // set out-of-range elements to 0.0f
 
    labDipole[0]      = atomIndex >= numOfAtoms ? 0.0f : molDipole[0]*vectorX[0] + molDipole[1]*vectorY[0] + molDipole[2]*vectorZ[0];
    labDipole[1]      = atomIndex >= numOfAtoms ? 0.0f : molDipole[0]*vectorX[1] + molDipole[1]*vectorY[1] + molDipole[2]*vectorZ[1];
    labDipole[2]      = atomIndex >= numOfAtoms ? 0.0f : molDipole[0]*vectorX[2] + molDipole[1]*vectorY[2] + molDipole[2]*vectorZ[2];
    
    // ---------------------------------------------------------------------------------------
    
    float* rPole[3];
    float mPole[3][3];
    float* labQuadrupole       = &(labFrameQuadrupole[atomIndex*9]);
    
    for( int ii = 0; ii < 3; ii++ ){
        mPole[ii][0]   = labQuadrupole[3*ii+0];
        mPole[ii][1]   = labQuadrupole[3*ii+1];
        mPole[ii][2]   = labQuadrupole[3*ii+2];

        rPole[ii]      = labQuadrupole + ii*3;
        rPole[ii][0]   = 0.0f;
        rPole[ii][1]   = 0.0f;
        rPole[ii][2]   = 0.0f;

    }
    
    int ii = threadIdx.x;
    if( ii < 3 ){
        for( int jj = ii; jj < 3; jj++ ){
 
            rPole[ii][jj] += vectorX[ii]*vectorX[jj]*mPole[0][0];
            rPole[ii][jj] += vectorX[ii]*vectorY[jj]*mPole[0][1];
            rPole[ii][jj] += vectorX[ii]*vectorZ[jj]*mPole[0][2];
       	
            rPole[ii][jj] += vectorY[ii]*vectorX[jj]*mPole[1][0];
            rPole[ii][jj] += vectorY[ii]*vectorY[jj]*mPole[1][1];
            rPole[ii][jj] += vectorY[ii]*vectorZ[jj]*mPole[1][2];
       	
            rPole[ii][jj] += vectorZ[ii]*vectorX[jj]*mPole[2][0];
            rPole[ii][jj] += vectorZ[ii]*vectorY[jj]*mPole[2][1];
            rPole[ii][jj] += vectorZ[ii]*vectorZ[jj]*mPole[2][2];
       }
    }
 
    __syncthreads();
 
 
    rPole[1][0] = rPole[0][1];
    rPole[2][0] = rPole[0][2];
    rPole[2][1] = rPole[1][2];
 
    // set out-of-range elements to 0.0f
 
    labQuadrupole[0]   = atomIndex >= numOfAtoms ? 0.0f : labQuadrupole[0];
    labQuadrupole[1]   = atomIndex >= numOfAtoms ? 0.0f : labQuadrupole[1];
    labQuadrupole[2]   = atomIndex >= numOfAtoms ? 0.0f : labQuadrupole[2];
    labQuadrupole[3]   = atomIndex >= numOfAtoms ? 0.0f : labQuadrupole[3];
    labQuadrupole[4]   = atomIndex >= numOfAtoms ? 0.0f : labQuadrupole[4];
    labQuadrupole[5]   = atomIndex >= numOfAtoms ? 0.0f : labQuadrupole[5];
    labQuadrupole[6]   = atomIndex >= numOfAtoms ? 0.0f : labQuadrupole[6];
    labQuadrupole[7]   = atomIndex >= numOfAtoms ? 0.0f : labQuadrupole[7];
    labQuadrupole[8]   = atomIndex >= numOfAtoms ? 0.0f : labQuadrupole[8];
}

void cudaComputeAmoebaLabFrameMoments( amoebaGpuContext amoebaGpu )
{

   // ---------------------------------------------------------------------------------------

   static const char* methodName = "computeCudaAmoebaLabFrameMoments";

   // ---------------------------------------------------------------------------------------

    gpuContext gpu    = amoebaGpu->gpuContext;

    int numBlocks     =  amoebaGpu->paddedNumberOfAtoms;
    int numThreads    =  20;

//#define AMOEBA_DEBUG  
#ifdef AMOEBA_DEBUG
    if( 0 && amoebaGpu->log ){
        (void) fprintf( amoebaGpu->log, "%s: numBlocks/atoms=%d\n", methodName, numBlocks ); (void) fflush( amoebaGpu->log );
        amoebaGpu->psMultipoleParticlesIdsAndAxisType->Download();
        amoebaGpu->psMolecularDipole->Download();
        gpu->psPosq4->Download();
        for( int ii = 0; ii < gpu->natoms; ii++ ){
            int mIndex = 3*ii;
             (void) fprintf( amoebaGpu->log,"%6d [%6d %6d %6d] x[%16.9e %16.9e %16.9e] dpl[%16.9e %16.9e %16.9e]\nRot[%16.9e %16.9e %16.9e] [%16.9e %16.9e %16.9e] [%16.9e %16.9e %16.9e]\n\n", ii,
                             amoebaGpu->psMultipoleParticlesIdsAndAxisType->_pSysData[ii].x,
                             amoebaGpu->psMultipoleParticlesIdsAndAxisType->_pSysData[ii].y,
                             amoebaGpu->psMultipoleParticlesIdsAndAxisType->_pSysData[ii].w,
                             gpu->psPosq4->_pSysData[ii].x,
                             gpu->psPosq4->_pSysData[ii].y,
                             gpu->psPosq4->_pSysData[ii].z,
                             amoebaGpu->psMolecularDipole->_pSysData[mIndex],
                             amoebaGpu->psMolecularDipole->_pSysData[mIndex+1],
                             amoebaGpu->psMolecularDipole->_pSysData[mIndex+2] );
        }
    }
//    int64 kernelTime = AmoebaTiming::getTimeOfDay();
    double kernelTime = 0.0;
#endif

    kCudaComputeCheckChiral_kernel<<< numBlocks, numThreads>>> ( );
    LAUNCHERROR("kCudaComputeCheckChiral");

    kCudaComputeLabFrameMoments_kernel<<< numBlocks, numThreads>>> ( );
    LAUNCHERROR(methodName);

#ifdef AMOEBA_DEBUG
    if( amoebaGpu->log ){
        static int timestep = 0;
        timestep++;
        (void) fprintf( amoebaGpu->log, "Finished rotation kernel execution in %lf us\n", kernelTime ); (void) fflush( amoebaGpu->log );
        (void) fflush( amoebaGpu->log );

        amoebaGpu->psLabFrameDipole->Download();
        (void) fprintf( amoebaGpu->log, "psLabFrameDipole completed\n" );  (void) fflush( amoebaGpu->log );

        amoebaGpu->psLabFrameQuadrupole->Download();
        (void) fprintf( amoebaGpu->log, "psLabFrameQpole completed\n" );  (void) fflush( amoebaGpu->log );

        int maxPrint = 10;
        for( int ii = 0; ii < amoebaGpu->paddedNumberOfAtoms; ii++ ){

             int dipoleOffset     = 3*ii;
             int quadrupoleOffset = 9*ii;

             (void) fprintf( amoebaGpu->log,"\n%6d [%6d %6d %6d] ", ii,
                             amoebaGpu->psMultipoleParticlesIdsAndAxisType->_pSysData[ii].x,
                             amoebaGpu->psMultipoleParticlesIdsAndAxisType->_pSysData[ii].y,
                             amoebaGpu->psMultipoleParticlesIdsAndAxisType->_pSysData[ii].w );
             // coords

             (void) fprintf( amoebaGpu->log,"x[%16.9e %16.9e %16.9e]\n",
                             gpu->psPosq4->_pSysData[ii].x,
                             gpu->psPosq4->_pSysData[ii].y,
                             gpu->psPosq4->_pSysData[ii].z);
/*
             (void) fprintf( amoebaGpu->log,"   R[%16.9e %16.9e %16.9e] [%16.9e %16.9e %16.9e] [%16.9e %16.9e %16.9e]\n",
                             amoebaGpu->psRotationMatrix->_pSysData[quadrupoleOffset],
                             amoebaGpu->psRotationMatrix->_pSysData[quadrupoleOffset+1],
                             amoebaGpu->psRotationMatrix->_pSysData[quadrupoleOffset+2],
                             amoebaGpu->psRotationMatrix->_pSysData[quadrupoleOffset+3],
                             amoebaGpu->psRotationMatrix->_pSysData[quadrupoleOffset+4],
                             amoebaGpu->psRotationMatrix->_pSysData[quadrupoleOffset+5],
                             amoebaGpu->psRotationMatrix->_pSysData[quadrupoleOffset+6],
                             amoebaGpu->psRotationMatrix->_pSysData[quadrupoleOffset+7],
                             amoebaGpu->psRotationMatrix->_pSysData[quadrupoleOffset+8] );
*/
             // dipole

             (void) fprintf( amoebaGpu->log,"   D[%16.9e %16.9e %16.9e]\n",
                             amoebaGpu->psLabFrameDipole->_pSysData[dipoleOffset],
                             amoebaGpu->psLabFrameDipole->_pSysData[dipoleOffset+1],
                             amoebaGpu->psLabFrameDipole->_pSysData[dipoleOffset+2] );
    
             // quadrupole

             (void) fprintf( amoebaGpu->log,"   Q[%16.9e %16.9e %16.9e] [%16.9e %16.9e %16.9e] [%16.9e %16.9e %16.9e]\n",
                             amoebaGpu->psLabFrameQuadrupole->_pSysData[quadrupoleOffset],
                             amoebaGpu->psLabFrameQuadrupole->_pSysData[quadrupoleOffset+1],
                             amoebaGpu->psLabFrameQuadrupole->_pSysData[quadrupoleOffset+2],
                             amoebaGpu->psLabFrameQuadrupole->_pSysData[quadrupoleOffset+3],
                             amoebaGpu->psLabFrameQuadrupole->_pSysData[quadrupoleOffset+4],
                             amoebaGpu->psLabFrameQuadrupole->_pSysData[quadrupoleOffset+5],
                             amoebaGpu->psLabFrameQuadrupole->_pSysData[quadrupoleOffset+6],
                             amoebaGpu->psLabFrameQuadrupole->_pSysData[quadrupoleOffset+7],
                             amoebaGpu->psLabFrameQuadrupole->_pSysData[quadrupoleOffset+8] );

            if( ii == maxPrint && (ii < (gpu->natoms - maxPrint)) ){
                ii = gpu->natoms - maxPrint;
            }
        }
        int nansDetected   = checkForNansAndInfinities( amoebaGpu->paddedNumberOfAtoms*3, amoebaGpu->psLabFrameDipole );
            nansDetected  += checkForNansAndInfinities( amoebaGpu->paddedNumberOfAtoms*9, amoebaGpu->psLabFrameQuadrupole );
        if( nansDetected ){
             (void) fprintf( amoebaGpu->log,"Nans detected in dipole/quadrupoles.\n" );
             exit(0);
        }
        (void) fflush( amoebaGpu->log );
    }
#endif

    if( 0 ){
//        int particles = particles;
        int particles = amoebaGpu->paddedNumberOfAtoms;
        std::vector<int> fileId;
        //fileId.push_back( 0 );
        VectorOfDoubleVectors outputVector;
        cudaLoadCudaFloat4Array( particles, 3, gpu->psPosq4,                     outputVector, gpu->psAtomIndex->_pSysData, 1.0f );
        cudaLoadCudaFloatArray( particles,  9, amoebaGpu->psRotationMatrix,      outputVector, gpu->psAtomIndex->_pSysData, 1.0f );
        cudaWriteVectorOfDoubleVectorsToFile( "CudaRotationMatrices", fileId, outputVector );
    }
    if( 0 ){

        int particles = amoebaGpu->paddedNumberOfAtoms;
        std::vector<int> fileId;
        //fileId.push_back( 0 );

        VectorOfDoubleVectors outputVector;
        cudaLoadCudaFloat4Array( particles, 3, gpu->psPosq4,                     outputVector, gpu->psAtomIndex->_pSysData, 1.0f );
        cudaLoadCudaFloatArray( particles,  3, amoebaGpu->psLabFrameDipole,      outputVector, gpu->psAtomIndex->_pSysData, 1.0f );
        cudaLoadCudaFloatArray( particles,  9, amoebaGpu->psLabFrameQuadrupole,  outputVector, gpu->psAtomIndex->_pSysData, 1.0f );
        cudaWriteVectorOfDoubleVectorsToFile( "CudaRotatedMoments", fileId, outputVector );
    }
  
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
            //compactStream(gpu->compactPlan, gpu->sim.pInteractingWorkUnit, gpu->sim.pWorkUnit, gpu->sim.pInteractionFlag, gpu->sim.workUnits, gpu->sim.pInteractionCount);
            compactStream(gpu->compactPlan, gpu->sim.pInteractingWorkUnit, amoebaGpu->psWorkUnit->_pDevData, gpu->sim.pInteractionFlag, gpu->sim.workUnits, gpu->sim.pInteractionCount);
            kFindInteractionsWithinBlocksPeriodic_kernel<<<gpu->sim.nonbond_blocks, gpu->sim.nonbond_threads_per_block,
                    sizeof(unsigned int)*gpu->sim.nonbond_threads_per_block>>>(gpu->sim.pInteractingWorkUnit);
            LAUNCHERROR("kFindInteractionsWithinBlocksPeriodic");
            cudaComputeAmoebaPmeFixedEField( amoebaGpu );
            cudaComputeAmoebaPmeMutualInducedField( amoebaGpu );
        }
    }

    // check if induce dipole calculation converged -- abort if it did not

    if( amoebaGpu->mutualInducedDone == 0 ){
       (void) fprintf( amoebaGpu->log, "%s induced dipole calculation did not converge -- aborting!\n", methodName.c_str() );
       (void) fflush( amoebaGpu->log );
       exit(-1);
    }

    // calculate electrostatic forces

    if( amoebaGpu->multipoleNonbondedMethod == AMOEBA_NO_CUTOFF ){

        cudaComputeAmoebaElectrostatic( amoebaGpu );

        // map torques to forces

        cudaComputeAmoebaMapTorquesAndAddTotalForce( amoebaGpu, amoebaGpu->psTorque, amoebaGpu->psForce, amoebaGpu->gpuContext->psForce4 );

        if( 0 ){
            gpuContext gpu = amoebaGpu->gpuContext;
            std::vector<int> fileId;
            //fileId.push_back( 0 );
            VectorOfDoubleVectors outputVector;
            //cudaLoadCudaFloat4Array( gpu->natoms, 3, gpu->psPosq4,              outputVector, gpu->psAtomIndex->_pSysData, 1.0f );
            cudaLoadCudaFloat4Array( gpu->natoms,  3, amoebaGpu->gpuContext->psForce4,     outputVector, gpu->psAtomIndex->_pSysData, 1.0f/4.184 );
            cudaWriteVectorOfDoubleVectorsToFile( "CudaMpole", fileId, outputVector );
        }   

    } else {
        cudaComputeAmoebaPmeElectrostatic( amoebaGpu );
    }
}

#undef AMOEBA_DEBUG
