//-----------------------------------------------------------------------------------------

//-----------------------------------------------------------------------------------------

#include "amoebaCudaKernels.h"

#include <stdio.h>
#include <cuda.h>
#include <cstdlib>
using namespace std; 

#define SQRT sqrtf

static __constant__ cudaGmxSimulation cSim;
static __constant__ cudaAmoebaGmxSimulation cAmoebaSim;

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

__device__ float normVector3( float* vector ) 
{

    float norm                    = DOT3( vector, vector );
    float returnNorm              = SQRT( norm );
    norm                          = returnNorm > 0.0f ? 1.0f/returnNorm : 0.0f;

    vector[0]                    *= norm;
    vector[1]                    *= norm;
    vector[2]                    *= norm;

    return returnNorm;
}

__global__
#if (__CUDA_ARCH__ >= 200)
__launch_bounds__(GF1XX_THREADS_PER_BLOCK, 1)
#elif (__CUDA_ARCH__ >= 130)
__launch_bounds__(GT2XX_THREADS_PER_BLOCK, 1)
#else
__launch_bounds__(G8X_THREADS_PER_BLOCK, 1)
#endif
void kCudaComputeLabFrameMoments_kernel(
				   int numOfAtoms,
				   float *rotationMatrix,
				   float4 *atomCoord,
				   int4 *multiPoleAtoms,
				   float *molecularDipole, float *molecularQuadrupole,
				   float *labFrameDipole,  float *labFrameQuadrupole )
{

   float* vectorX;
   float* vectorY;
   float* vectorZ;

   // ---------------------------------------------------------------------------------------

   int atomIndex = blockIdx.x;//__mul24(blockIdx.x,blockDim.x) + threadIdx.x ;

   // ---------------------------------------------------------------------------------------

   // get coordinates of this atom and the z & x axis atoms
   // compute the vector between the atoms and 1/sqrt(d2), d2 is distance between
   // this atom and the axis atom

   // this atom is referred to as the k-atom in notes below

   // code common to ZThenX and Bisector
   
   vectorX                          = &(rotationMatrix[atomIndex*9]);
   vectorY                          = &(rotationMatrix[atomIndex*9+ 3]);
   vectorZ                          = &(rotationMatrix[atomIndex*9+ 6]);

   float4 coordinatesThisAtom       = atomCoord[atomIndex];

   int multipoleAtomIndex           = multiPoleAtoms[atomIndex].x;
   float4 coordinatesAxisAtom       = atomCoord[multipoleAtomIndex];

   vectorZ[0]                       = coordinatesAxisAtom.x - coordinatesThisAtom.x;
   vectorZ[1]                       = coordinatesAxisAtom.y - coordinatesThisAtom.y;
   vectorZ[2]                       = coordinatesAxisAtom.z - coordinatesThisAtom.z;
     
   multipoleAtomIndex               = multiPoleAtoms[atomIndex].y; 
   coordinatesAxisAtom              = atomCoord[multipoleAtomIndex];

   vectorX[0]                       = coordinatesAxisAtom.x - coordinatesThisAtom.x;
   vectorX[1]                       = coordinatesAxisAtom.y - coordinatesThisAtom.y;
   vectorX[2]                       = coordinatesAxisAtom.z - coordinatesThisAtom.z;

   int axisType                     = multiPoleAtoms[atomIndex].w; 
     
   float sum                        = normVector3( vectorZ );
   
   // branch based on axis type
    
   if( axisType == 1 ){

     // bisector
     // dx = dx1 + dx2 (in Tinker code)
     
     sum               = normVector3( vectorX );
     
     vectorZ[0]       += vectorX[0];
     vectorZ[1]       += vectorX[1];
     vectorZ[2]       += vectorX[2];

     sum               = normVector3( vectorZ );

   }
   
   float dot         = vectorZ[0]*vectorX[0] + vectorZ[1]*vectorX[1] + vectorZ[2]*vectorX[2];
   
   vectorX[0]       -= dot*vectorZ[0];
   vectorX[1]       -= dot*vectorZ[1];
   vectorX[2]       -= dot*vectorZ[2];

   sum               = normVector3( vectorX );
   
   vectorY[0]        = (vectorZ[1]*vectorX[2]) - (vectorZ[2]*vectorX[1]);
   vectorY[1]        = (vectorZ[2]*vectorX[0]) - (vectorZ[0]*vectorX[2]);
   vectorY[2]        = (vectorZ[0]*vectorX[1]) - (vectorZ[1]*vectorX[0]);

   float* molDipole  = &(molecularDipole[atomIndex*3]);
   float* labDipole  = &(labFrameDipole[atomIndex*3]);
   
   // set out-of-range elements to 0.0f

   labDipole[0]      = atomIndex >= numOfAtoms ? 0.0f : molDipole[0]*vectorX[0] + molDipole[1]*vectorY[0] + molDipole[2]*vectorZ[0];
   labDipole[1]      = atomIndex >= numOfAtoms ? 0.0f : molDipole[0]*vectorX[1] + molDipole[1]*vectorY[1] + molDipole[2]*vectorZ[1];
   labDipole[2]      = atomIndex >= numOfAtoms ? 0.0f : molDipole[0]*vectorX[2] + molDipole[1]*vectorY[2] + molDipole[2]*vectorZ[2];
   
   // ---------------------------------------------------------------------------------------
   
   const float * mPole[3];
   float* rPole[3];
   
   float* molQuadrupole       = &(molecularQuadrupole[atomIndex*9]);
   float* labQuadrupole       = &(labFrameQuadrupole[atomIndex*9]);
   
   for( int ii = 0; ii < 3; ii++ ){
      mPole[ii]    = molQuadrupole + ii*3;
      rPole[ii]    = labQuadrupole + ii*3;
      rPole[ii][0] = rPole[ii][1] = rPole[ii][2] = 0.0f;
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
                             amoebaGpu->psMultipoleParticlesIdsAndAxisType->_pSysStream[0][ii].x,
                             amoebaGpu->psMultipoleParticlesIdsAndAxisType->_pSysStream[0][ii].y,
                             amoebaGpu->psMultipoleParticlesIdsAndAxisType->_pSysStream[0][ii].w,
                             gpu->psPosq4->_pSysStream[0][ii].x,
                             gpu->psPosq4->_pSysStream[0][ii].y,
                             gpu->psPosq4->_pSysStream[0][ii].z,
                             amoebaGpu->psMolecularDipole->_pSysStream[0][mIndex],
                             amoebaGpu->psMolecularDipole->_pSysStream[0][mIndex+1],
                             amoebaGpu->psMolecularDipole->_pSysStream[0][mIndex+2] );
        }
    }
//    int64 kernelTime = AmoebaTiming::getTimeOfDay();
    double kernelTime = 0.0;
#endif

    kCudaComputeLabFrameMoments_kernel<<< numBlocks, numThreads>>> (
       gpu->natoms,
       amoebaGpu->psRotationMatrix->_pDevStream[0],
       gpu->psPosq4->_pDevStream[0],
       amoebaGpu->psMultipoleParticlesIdsAndAxisType->_pDevStream[0],
       amoebaGpu->psMolecularDipole->_pDevStream[0],
       amoebaGpu->psMolecularQuadrupole->_pDevStream[0],
       amoebaGpu->psLabFrameDipole->_pDevStream[0],
       amoebaGpu->psLabFrameQuadrupole->_pDevStream[0] 
       );
    LAUNCHERROR(methodName);

#ifdef AMOEBA_DEBUG
    if( amoebaGpu->log ){
//        kernelTime          = AmoebaTiming::getTimeOfDay() - kernelTime;
        static int timestep = 0;
        timestep++;
        (void) fprintf( amoebaGpu->log, "Finished rotation kernel execution in %lf us\n", kernelTime ); (void) fflush( amoebaGpu->log );
        (void) fprintf( amoebaGpu->log, "psLabFrameDipole=%p _pSysStream=%p _pSysStream[0]=%p _pDevStream=%p _pDevStream[0]=%p\n",
                        amoebaGpu->psLabFrameDipole,  amoebaGpu->psLabFrameDipole->_pSysStream, 
                        amoebaGpu->psLabFrameDipole->_pSysStream[0], amoebaGpu->psLabFrameDipole->_pDevStream, amoebaGpu->psLabFrameDipole->_pDevStream[0] );
        fflush( amoebaGpu->log );

        amoebaGpu->psRotationMatrix->Download();
        amoebaGpu->psLabFrameDipole->Download();
        (void) fprintf( amoebaGpu->log, "psLabFrameDipole completed\n" );  (void) fflush( amoebaGpu->log );

        amoebaGpu->psLabFrameQuadrupole->Download();
        (void) fprintf( amoebaGpu->log, "psLabFrameQpole completed\n" );  (void) fflush( amoebaGpu->log );

        int maxPrint = 10;
        for( int ii = 0; ii < amoebaGpu->paddedNumberOfAtoms; ii++ ){

             int dipoleOffset     = 3*ii;
             int quadrupoleOffset = 9*ii;

             (void) fprintf( amoebaGpu->log,"\n%6d [%6d %6d %6d] ", ii,
                             amoebaGpu->psMultipoleParticlesIdsAndAxisType->_pSysStream[0][ii].x,
                             amoebaGpu->psMultipoleParticlesIdsAndAxisType->_pSysStream[0][ii].y,
                             amoebaGpu->psMultipoleParticlesIdsAndAxisType->_pSysStream[0][ii].w );
             // coords

             (void) fprintf( amoebaGpu->log,"x[%16.9e %16.9e %16.9e]\n",
                             gpu->psPosq4->_pSysStream[0][ii].x,
                             gpu->psPosq4->_pSysStream[0][ii].y,
                             gpu->psPosq4->_pSysStream[0][ii].z);

             (void) fprintf( amoebaGpu->log,"   R[%16.9e %16.9e %16.9e] [%16.9e %16.9e %16.9e] [%16.9e %16.9e %16.9e]\n",
                             amoebaGpu->psRotationMatrix->_pSysStream[0][quadrupoleOffset],
                             amoebaGpu->psRotationMatrix->_pSysStream[0][quadrupoleOffset+1],
                             amoebaGpu->psRotationMatrix->_pSysStream[0][quadrupoleOffset+2],
                             amoebaGpu->psRotationMatrix->_pSysStream[0][quadrupoleOffset+3],
                             amoebaGpu->psRotationMatrix->_pSysStream[0][quadrupoleOffset+4],
                             amoebaGpu->psRotationMatrix->_pSysStream[0][quadrupoleOffset+5],
                             amoebaGpu->psRotationMatrix->_pSysStream[0][quadrupoleOffset+6],
                             amoebaGpu->psRotationMatrix->_pSysStream[0][quadrupoleOffset+7],
                             amoebaGpu->psRotationMatrix->_pSysStream[0][quadrupoleOffset+8] );

             // dipole

             (void) fprintf( amoebaGpu->log,"   D[%16.9e %16.9e %16.9e]\n",
                             amoebaGpu->psLabFrameDipole->_pSysStream[0][dipoleOffset],
                             amoebaGpu->psLabFrameDipole->_pSysStream[0][dipoleOffset+1],
                             amoebaGpu->psLabFrameDipole->_pSysStream[0][dipoleOffset+2] );
    
             // quadrupole

             (void) fprintf( amoebaGpu->log,"   Q[%16.9e %16.9e %16.9e] [%16.9e %16.9e %16.9e] [%16.9e %16.9e %16.9e]\n",
                             amoebaGpu->psLabFrameQuadrupole->_pSysStream[0][quadrupoleOffset],
                             amoebaGpu->psLabFrameQuadrupole->_pSysStream[0][quadrupoleOffset+1],
                             amoebaGpu->psLabFrameQuadrupole->_pSysStream[0][quadrupoleOffset+2],
                             amoebaGpu->psLabFrameQuadrupole->_pSysStream[0][quadrupoleOffset+3],
                             amoebaGpu->psLabFrameQuadrupole->_pSysStream[0][quadrupoleOffset+4],
                             amoebaGpu->psLabFrameQuadrupole->_pSysStream[0][quadrupoleOffset+5],
                             amoebaGpu->psLabFrameQuadrupole->_pSysStream[0][quadrupoleOffset+6],
                             amoebaGpu->psLabFrameQuadrupole->_pSysStream[0][quadrupoleOffset+7],
                             amoebaGpu->psLabFrameQuadrupole->_pSysStream[0][quadrupoleOffset+8] );

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
        cudaLoadCudaFloat4Array( particles, 3, gpu->psPosq4,                     outputVector );
        cudaLoadCudaFloatArray( particles,  9, amoebaGpu->psRotationMatrix,      outputVector );
        cudaWriteVectorOfDoubleVectorsToFile( "CudaRotationMatrices", fileId, outputVector );
    }
    if( 0 ){

        int particles = amoebaGpu->paddedNumberOfAtoms;
        std::vector<int> fileId;
        //fileId.push_back( 0 );

        VectorOfDoubleVectors outputVector;
        cudaLoadCudaFloat4Array( particles, 3, gpu->psPosq4,                     outputVector );
        cudaLoadCudaFloatArray( particles,  3, amoebaGpu->psLabFrameDipole,      outputVector );
        cudaLoadCudaFloatArray( particles,  9, amoebaGpu->psLabFrameQuadrupole,  outputVector );
        cudaWriteVectorOfDoubleVectorsToFile( "CudaRotatedMoments", fileId, outputVector );
    }
  
}

void kCalculateAmoebaMultipoleForces(amoebaGpuContext amoebaGpu, bool hasAmoebaGeneralizedKirkwood ) 
{
    std::string methodName = "kCalculateAmoebaMultipoleForces";
    //printf("%s \n", methodName.c_str() ); fflush( stdout );

    // compute lab frame moments

    cudaComputeAmoebaLabFrameMoments( amoebaGpu );

    // compute fixed E-field and mutual induced field 

    if( hasAmoebaGeneralizedKirkwood ){
        cudaComputeAmoebaFixedEAndGkFields( amoebaGpu );
        cudaComputeAmoebaMutualInducedAndGkField( amoebaGpu );
    } else {
        cudaComputeAmoebaFixedEField( amoebaGpu );
        cudaComputeAmoebaMutualInducedField( amoebaGpu );
    }

    // check if induce dipole calculation converged -- abort if it did not

    if( amoebaGpu->mutualInducedDone ){
       //cudaComputeAmoebaElectrostatic( amoebaGpuContextGlobal );
    } else {
       (void) fprintf( amoebaGpu->log, "%s induced dipole calculation did not converge -- aborting!\n", methodName.c_str() );
       (void) fflush( amoebaGpu->log );
       exit(-1);
    }

    // calculate electrostatic forces

    cudaComputeAmoebaElectrostatic( amoebaGpu );

    // map torques to forces

    cudaComputeAmoebaMapTorquesAndAddTotalForce( amoebaGpu, amoebaGpu->psTorque, amoebaGpu->psForce, amoebaGpu->gpuContext->psForce4 );
   
    if( 0 && amoebaGpu->log ){
        (void) fprintf( amoebaGpu->log, "Done mapping torques -> forces%s\n", methodName.c_str() ); fflush( NULL );
        (void) fflush( NULL );
    }
}

#undef AMOEBA_DEBUG
