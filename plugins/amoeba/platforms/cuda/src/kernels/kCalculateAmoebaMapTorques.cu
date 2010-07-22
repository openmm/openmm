#include "amoebaCudaKernels.h"

//#define AMOEBA_DEBUG
#define BLOCK_SIZE 128

using namespace std; 

static __constant__ cudaGmxSimulation cSim;
static __constant__ cudaAmoebaGmxSimulation cAmoebaSim;

void SetCalculateAmoebaCudaMapTorquesSim(amoebaGpuContext amoebaGpu)
{
    cudaError_t status;
    gpuContext gpu = amoebaGpu->gpuContext;
    status         = cudaMemcpyToSymbol(cSim, &gpu->sim, sizeof(cudaGmxSimulation));    
    RTERROR(status, "SetCalculateAmoebaCudaMapTorquesSim: cudaMemcpyToSymbol: SetSim copy to cSim failed");
    status         = cudaMemcpyToSymbol(cAmoebaSim, &amoebaGpu->amoebaSim, sizeof(cudaAmoebaGmxSimulation));    
    RTERROR(status, "SetCalculateAmoebaCudaMapTorquesSim: cudaMemcpyToSymbol: SetSim copy to cAmoebaSim failed");
}

void GetCalculateAmoebaCudaMapTorquesSim(amoebaGpuContext amoebaGpu)
{
    cudaError_t status;
    gpuContext gpu = amoebaGpu->gpuContext;
    status = cudaMemcpyFromSymbol(&gpu->sim, cSim, sizeof(cudaGmxSimulation));    
    RTERROR(status, "GetCalculateAmoebaCudaMapTorquesSim: cudaMemcpyFromSymbol: SetSim copy from cSim failed");
    status = cudaMemcpyFromSymbol(&amoebaGpu->amoebaSim, cAmoebaSim, sizeof(cudaAmoebaGmxSimulation));    
    RTERROR(status, "GetCalculateAmoebaCudaMapTorquesSim: cudaMemcpyFromSymbol: SetSim copy from cAmoebaSim failed");
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
void amoebaMapTorqueToForce_kernel( 
				    int numOfAtoms,
				    float4* atomCoord,
				    float* torque,
				    int4* multiPoleAtoms,
				    int maxDiff,
				    float* tempElecForce
				    ){

  // ---------------------------------------------------------------------------------------
  
  int threadId = __mul24(blockIdx.x,blockDim.x) + threadIdx.x ;
  if(threadId >= numOfAtoms)return;
  
  int U     = 0;
  int V     = 1;
  int W     = 2;
  
  int X     = 0;
  int Y     = 1;
  int Z     = 2;
  
  float forces[3][3];
  float norms[3];
  float vector[3][3];
  
  // ---------------------------------------------------------------------------------------
  
    int axisAtom                  = multiPoleAtoms[threadId].x;

    vector[U][0]                  = atomCoord[axisAtom].x - atomCoord[threadId].x;
    vector[U][1]                  = atomCoord[axisAtom].y - atomCoord[threadId].y;
    vector[U][2]                  = atomCoord[axisAtom].z - atomCoord[threadId].z;

    norms[U]                      = normVector3( vector[U] );

    axisAtom                      = multiPoleAtoms[threadId].y;

    vector[V][0]                  = atomCoord[axisAtom].x - atomCoord[threadId].x;
    vector[V][1]                  = atomCoord[axisAtom].y - atomCoord[threadId].y;
    vector[V][2]                  = atomCoord[axisAtom].z - atomCoord[threadId].z;

    norms[V]                      = normVector3( vector[V] );

    // W = UxV

    vector[W][0]                  = vector[U][1]*vector[V][2] - vector[U][2]*vector[V][1];
    vector[W][1]                  = vector[U][2]*vector[V][0] - vector[U][0]*vector[V][2];
    vector[W][2]                  = vector[U][0]*vector[V][1] - vector[U][1]*vector[V][0];

    norms[W]                      = normVector3( vector[W] );

    float diff[3];
    diff[0]                       = vector[V][0] - vector[U][0];
    diff[1]                       = vector[V][1] - vector[U][1];
    diff[2]                       = vector[V][2] - vector[U][2];

    float dotDu                   = DOT3( vector[U], diff );
    float dotDv                   = DOT3( vector[V], diff );

    float up[3], vp[3];

    up[0]                         = diff[0] - dotDu*vector[U][0];
    vp[0]                         = diff[0] - dotDv*vector[V][0];

    up[1]                         = diff[1] - dotDu*vector[U][1];
    vp[1]                         = diff[1] - dotDv*vector[V][1];

    up[2]                         = diff[2] - dotDu*vector[U][2];
    vp[2]                         = diff[2] - dotDv*vector[V][2];

    float norm                    = normVector3( up );
    norm                          = normVector3( vp );

    float dphi[3];
    dphi[0]                       = vector[0][0]*torque[threadId*3] + vector[0][1]*torque[threadId*3+1] + vector[0][2]*torque[threadId*3+2];
    dphi[1]                       = vector[1][0]*torque[threadId*3] + vector[1][1]*torque[threadId*3+1] + vector[1][2]*torque[threadId*3+2];
    dphi[2]                       = vector[2][0]*torque[threadId*3] + vector[2][1]*torque[threadId*3+1] + vector[2][2]*torque[threadId*3+2];

    // clamp c to interval [-1,1]

    float c                       = DOT3( vector[U], vector[V] );
          c                       = c >  1.0f ?  1.0f : c;
          c                       = c < -1.0f ? -1.0f : c;

    float s                       = SQRT( 1.0f - (c*c) );
    float uvdis                   = norms[U]*s;
    float vudis                   = norms[V]*s;

    float factorX;
    float factorZ;
    if( multiPoleAtoms[threadId].w == 1 ){       
        factorX = 0.5f;
        factorZ = 0.5f;
    } else {
        factorX = 1.0f;
        factorZ = 0.0f;
    }

    forces[X][0]                  = -vector[W][0]*dphi[V]/uvdis + factorX*up[0]*dphi[W]/norms[U];
    forces[Z][0]                  =  vector[W][0]*dphi[U]/vudis + factorZ*vp[0]*dphi[W]/norms[V];

    forces[X][1]                  = -vector[W][1]*dphi[V]/uvdis + factorX*up[1]*dphi[W]/norms[U];
    forces[Z][1]                  =  vector[W][1]*dphi[U]/vudis + factorZ*vp[1]*dphi[W]/norms[V];

    forces[X][2]                  = -vector[W][2]*dphi[V]/uvdis + factorX*up[2]*dphi[W]/norms[U];
    forces[Z][2]                  =  vector[W][2]*dphi[U]/vudis + factorZ*vp[2]*dphi[W]/norms[V];

    forces[Y][0]                  = -(forces[X][0] + forces[Z][0]);
    forces[Y][1]                  = -(forces[X][1] + forces[Z][1]);
    forces[Y][2]                  = -(forces[X][2] + forces[Z][2]);

    int temp                      = multiPoleAtoms[threadId].x;
    int min                       = multiPoleAtoms[temp].z;
    int offset                    = 3*(temp*maxDiff + threadId-min);
    tempElecForce[offset  ]       = forces[X][0];
    tempElecForce[offset+1]       = forces[X][1];
    tempElecForce[offset+2]       = forces[X][2];
   
    temp                          = multiPoleAtoms[threadId].y;
    min                           = multiPoleAtoms[temp].z;
    offset                        = 3*(temp*maxDiff + threadId-min);
    tempElecForce[offset  ]       = forces[Z][0];
    tempElecForce[offset+1]       = forces[Z][1];
    tempElecForce[offset+2]       = forces[Z][2];
 
    min                           = multiPoleAtoms[threadId].z;
    offset                        = 3*(threadId*(maxDiff + 1)-min); 
    tempElecForce[offset  ]       = forces[Y][0];
    tempElecForce[offset+1]       = forces[Y][1];
    tempElecForce[offset+2]       = forces[Y][2];
   
}

__global__
#if (__CUDA_ARCH__ >= 200)
__launch_bounds__(GF1XX_THREADS_PER_BLOCK, 1)
#elif (__CUDA_ARCH__ >= 130)
__launch_bounds__(GT2XX_THREADS_PER_BLOCK, 1)
#else
__launch_bounds__(G8X_THREADS_PER_BLOCK, 1)
#endif
void amoebaMapTorqueReduce_kernel( 
					   int numThreads,
					   int numOfAtoms,
					   int maxDiff,
					   float* tempElecForce,
					   float* elecForce
					   ){
  
    unsigned int tid = threadIdx.x;
    
    __shared__ float sfx[BLOCK_SIZE];
    __shared__ float sfy[BLOCK_SIZE];
    __shared__ float sfz[BLOCK_SIZE];
  
    // load values then sum and add results to elecForce

    if( tid < maxDiff ){
  
        sfx[tid] = tempElecForce[blockIdx.x*3*maxDiff + tid*3  ];
        sfy[tid] = tempElecForce[blockIdx.x*3*maxDiff + tid*3+1];
        sfz[tid] = tempElecForce[blockIdx.x*3*maxDiff + tid*3+2];
  
    } else {
        sfx[tid] = sfy[tid] = sfz[tid] = 0.0f;
    }
  
    __syncthreads();
  
    for( unsigned int s = (blockDim.x)/2; s != 0; s >>= 1 ){
        if(tid<s){
            sfx[tid] += sfx[tid+s];
            sfy[tid] += sfy[tid+s];
            sfz[tid] += sfz[tid+s];
        }
        __syncthreads();
    }
  
    if( tid == 0 ){
        elecForce[blockIdx.x*3]    += sfx[0];
        elecForce[blockIdx.x*3+1]  += sfy[0];
        elecForce[blockIdx.x*3+2]  += sfz[0];
    }  

}

__global__
#if (__CUDA_ARCH__ >= 200)
__launch_bounds__(GF1XX_THREADS_PER_BLOCK, 1)
#elif (__CUDA_ARCH__ >= 130)
__launch_bounds__(GT2XX_THREADS_PER_BLOCK, 1)
#else
__launch_bounds__(G8X_THREADS_PER_BLOCK, 1)
#endif
void amoebaMapTorqueReduce_kernel2( 
					   int numThreads,
					   int numOfAtoms,
					   int maxDiff,
					   float* tempElecForce,
					   float* elecForce,
					   float4* outputForce
					   ){
  
    unsigned int tid = threadIdx.x;
    
    __shared__ float sfx[BLOCK_SIZE];
    __shared__ float sfy[BLOCK_SIZE];
    __shared__ float sfz[BLOCK_SIZE];
  
    // load values then sum and add results to elecForce

    if( tid < maxDiff ){
  
        sfx[tid] = tempElecForce[blockIdx.x*3*maxDiff + tid*3  ];
        sfy[tid] = tempElecForce[blockIdx.x*3*maxDiff + tid*3+1];
        sfz[tid] = tempElecForce[blockIdx.x*3*maxDiff + tid*3+2];
  
    } else {
        sfx[tid] = sfy[tid] = sfz[tid] = 0.0f;
    }
  
    __syncthreads();
  
    for( unsigned int s = (blockDim.x)/2; s != 0; s >>= 1 ){
        if(tid<s){
            sfx[tid] += sfx[tid+s];
            sfy[tid] += sfy[tid+s];
            sfz[tid] += sfz[tid+s];
        }
        __syncthreads();
    }
  
    if( tid == 0 ){
        outputForce[blockIdx.x].x  += elecForce[blockIdx.x*3  ] + sfx[0];
        outputForce[blockIdx.x].y  += elecForce[blockIdx.x*3+1] + sfy[0];
        outputForce[blockIdx.x].z  += elecForce[blockIdx.x*3+2] + sfz[0];
    }  

}

void cudaComputeAmoebaMapTorques( amoebaGpuContext amoebaGpu, CUDAStream<float>* psTorque, CUDAStream<float>* psForce )
{
  
   // ---------------------------------------------------------------------------------------

//#define AMOEBA_DEBUG
#ifdef AMOEBA_DEBUG
    static const char* methodName = "cudaMapAmoebaTorqueToForce";
    static int timestep = 0; 
    std::vector<int> fileId;
    timestep++;
    fileId.resize( 2 ); 
    fileId[0] = timestep;
    fileId[1] = 1; 
#endif

    // check that BLOCK_SIZE >= amoebaGpu->maxMapTorqueDifference

    if( amoebaGpu->maxMapTorqueDifference > BLOCK_SIZE ){
        (void) fprintf( amoebaGpu->log, "block size (%d) in amoebaMapTorqueReduce_kernel is too small ( > %d)! -- aborting.\n",
                        BLOCK_SIZE, amoebaGpu->maxMapTorqueDifference );
        exit(-1);  
    }

   // ---------------------------------------------------------------------------------------
    
    gpuContext gpu    = amoebaGpu->gpuContext;
    int numThreads    = min(256, (gpu->natoms));
    int numBlocks     =  1 + (gpu->natoms/numThreads);

//#ifdef AMOEBA_DEBUG
#if 0
    if( amoebaGpu->log ){
        (void) fprintf( amoebaGpu->log, "%s: numBlocks=%d numThreads=%d\n", methodName, numBlocks, numThreads ); (void) fflush( amoebaGpu->log );
    }
    amoebaGpu->psForce->Download();
    amoebaGpu->psTorque->Download();
    int maxPrint        = 20;
    (void) fprintf( amoebaGpu->log,"Pre torqueMap\n" );
    for( int ii = 0; ii < gpu->natoms; ii++ ){
        (void) fprintf( amoebaGpu->log, "%5d ", ii);

        int indexOffset     = ii*3;

        (void) fprintf( amoebaGpu->log,"E[%16.9e %16.9e %16.9e] ",
                        amoebaGpu->psForce->_pSysStream[0][indexOffset],
                        amoebaGpu->psForce->_pSysStream[0][indexOffset+1],
                        amoebaGpu->psForce->_pSysStream[0][indexOffset+2] );
        (void) fprintf( amoebaGpu->log,"T[%16.9e %16.9e %16.9e]\n",
                        amoebaGpu->psTorque->_pSysStream[0][indexOffset],
                        amoebaGpu->psTorque->_pSysStream[0][indexOffset+1],
                        amoebaGpu->psTorque->_pSysStream[0][indexOffset+2] );
        if( ii == maxPrint && (gpu->natoms - maxPrint) > ii )ii = gpu->natoms - maxPrint;
    }
    int nansDetected   = checkForNansAndInfinities( gpu->natoms*3, amoebaGpu->psForce );
    nansDetected      += checkForNansAndInfinities( gpu->natoms*3, amoebaGpu->psTorque );
    if( nansDetected ){
        (void) fprintf( amoebaGpu->log,"WARNING: %d nans/infinities detected force/torques.\n", nansDetected );
    } else {
        (void) fprintf( amoebaGpu->log,"No nans/infinities detected in force/torques.\n" );
    }

// zero forces
#if 0
    for( int ii = 0; ii < 3*gpu->natoms; ii++ ){
        amoebaGpu->psForce->_pSysStream[0][ii] = 0.0f;
    }
    amoebaGpu->psForce->Upload();
#endif

#endif

    // torqueMapForce is zeroed when initialized; should not need to be reinitialized 

/*
    AmoebaTorqueMapZeroKernel<<< numBlocks, numThreads >>>( 
       gpu->natoms, amoebaGpu->torqueMapForce->_pDevStream[0] );
    LAUNCHERROR("AmoebaMapTrqZeroKernel");  
*/


    amoebaMapTorqueToForce_kernel<<< numBlocks, numThreads>>> (
                gpu->natoms,
                gpu->psPosq4->_pDevStream[0],
                psTorque->_pDevStream[0],
                amoebaGpu->psMultipoleParticlesIdsAndAxisType->_pDevStream[0],
                amoebaGpu->maxMapTorqueDifference,
                amoebaGpu->torqueMapForce->_pDevStream[0] );
    LAUNCHERROR("AmoebaMapTrqKernel");
  
//#ifdef AMOEBA_DEBUG
#if 0
    amoebaGpu->torqueMapForce->Download();
    //int maxPrint        = 10;
    (void) fprintf( amoebaGpu->log,"Post AmoebaMapTrqKernel maxMapTorqueDifference=%d\n", amoebaGpu->maxMapTorqueDifference );
    for( int ii = 0; ii < gpu->natoms; ii++ ){
        (void) fprintf( amoebaGpu->log, "\n%5d multi[%d %d %d %d]\n", ii, 
                        amoebaGpu->psMultipoleParticlesIdsAndAxisType->_pSysStream[0][ii].x,
                        amoebaGpu->psMultipoleParticlesIdsAndAxisType->_pSysStream[0][ii].y,
                        amoebaGpu->psMultipoleParticlesIdsAndAxisType->_pSysStream[0][ii].z,
                        amoebaGpu->psMultipoleParticlesIdsAndAxisType->_pSysStream[0][ii].w );

        int indexOffset     = ii*3*amoebaGpu->maxMapTorqueDifference;
        float sum[3]        = { 0.0f,  0.0f,  0.0f };
        for( int jj = 0; jj < amoebaGpu->maxMapTorqueDifference; jj++ ){
            if( amoebaGpu->torqueMapForce->_pSysStream[0][indexOffset] != 0.0f ){ 
                (void) fprintf( amoebaGpu->log,"    %4d %4d Temp[%16.9e %16.9e %16.9e] %d\n",
                                ii, jj + amoebaGpu->psMultipoleParticlesIdsAndAxisType->_pSysStream[0][ii].z,
                                amoebaGpu->torqueMapForce->_pSysStream[0][indexOffset],
                                amoebaGpu->torqueMapForce->_pSysStream[0][indexOffset+1],
                                amoebaGpu->torqueMapForce->_pSysStream[0][indexOffset+2], indexOffset );
                sum[0] += amoebaGpu->torqueMapForce->_pSysStream[0][indexOffset];
                sum[1] += amoebaGpu->torqueMapForce->_pSysStream[0][indexOffset+1];
                sum[2] += amoebaGpu->torqueMapForce->_pSysStream[0][indexOffset+2];
            }
            indexOffset += 3;
        }
        (void) fprintf( amoebaGpu->log,"               Sum[%16.9e %16.9e %16.9e]\n", sum[0], sum[1], sum[2] );
        if( ii == maxPrint && (gpu->natoms - maxPrint) > ii )ii = gpu->natoms - maxPrint;
    }
#endif

    numBlocks  = gpu->natoms;
    numThreads = amoebaGpu->maxMapTorqueDifferencePow2;

    amoebaMapTorqueReduce_kernel<<< numBlocks, numThreads>>>(
            numThreads, gpu->natoms,
            amoebaGpu->maxMapTorqueDifference,
            amoebaGpu->torqueMapForce->_pDevStream[0],
            psForce->_pDevStream[0] );
    LAUNCHERROR("amoebaMapTorqueReduce_kernel");

#ifdef AMOEBA_DEBUG
    if( 0 && amoebaGpu->log ){
        (void) fprintf( amoebaGpu->log, "%s: numBlocks=%d numThreads=%d %d\n", methodName, numBlocks, numThreads, amoebaGpu->maxMapTorqueDifferencePow2); (void) fflush( amoebaGpu->log );
        amoebaGpu->psForce->Download();
        amoebaGpu->psTorque->Download();
        int maxPrint        = 10;
        (void) fprintf( amoebaGpu->log,"Post torqueMap\n" );
        for( int ii = 0; ii < gpu->natoms; ii++ ){
            (void) fprintf( amoebaGpu->log, "%5d ", ii);
    
            int indexOffset     = ii*3;
    
            (void) fprintf( amoebaGpu->log,"E[%16.9e %16.9e %16.9e] ",
                            amoebaGpu->psForce->_pSysStream[0][indexOffset],
                            amoebaGpu->psForce->_pSysStream[0][indexOffset+1],
                            amoebaGpu->psForce->_pSysStream[0][indexOffset+2] );
            (void) fprintf( amoebaGpu->log,"T[%16.9e %16.9e %16.9e]\n",
                            amoebaGpu->psTorque->_pSysStream[0][indexOffset],
                            amoebaGpu->psTorque->_pSysStream[0][indexOffset+1],
                            amoebaGpu->psTorque->_pSysStream[0][indexOffset+2] );
            if( ii == maxPrint && (gpu->natoms - maxPrint) > ii )ii = gpu->natoms - maxPrint;
        }
        (void) fflush( amoebaGpu->log );
    }
    if( 1 ){
        //std::vector<int> fileId;
        //fileId.push_back( 0 );
        VectorOfDoubleVectors outputVector;
        cudaLoadCudaFloat4Array( gpu->natoms, 3, gpu->psPosq4,              outputVector );
        cudaLoadCudaFloatArray( gpu->natoms,  3, amoebaGpu->psForce,      outputVector );
        cudaLoadCudaFloatArray( gpu->natoms,  3, amoebaGpu->psTorque, outputVector);
        cudaWriteVectorOfDoubleVectorsToFile( "CudaVacuumElecForce", fileId, outputVector );
    }
#endif

}

void cudaComputeAmoebaMapTorquesAndAddTotalForce( amoebaGpuContext amoebaGpu,
                                                  CUDAStream<float>* psTorque,
                                                  CUDAStream<float>* psForce,
                                                  CUDAStream<float4>* psCudaForce4 )
{
  
   // ---------------------------------------------------------------------------------------

//#define AMOEBA_DEBUG
#ifdef AMOEBA_DEBUG
    static const char* methodName = "cudaComputeAmoebaMapTorquesAndAddTotalForce";
    static int timestep = 0; 
    std::vector<int> fileId;
    timestep++;
    fileId.resize( 2 ); 
    fileId[0] = timestep;
    fileId[1] = 1; 
#endif

    // check that BLOCK_SIZE >= amoebaGpu->maxMapTorqueDifference

    if( amoebaGpu->maxMapTorqueDifference > BLOCK_SIZE ){
        (void) fprintf( amoebaGpu->log, "block size (%d) in amoebaMapTorqueReduce_kernel is too small ( > %d)! -- aborting.\n",
                        BLOCK_SIZE, amoebaGpu->maxMapTorqueDifference );
        exit(-1);  
    }

   // ---------------------------------------------------------------------------------------
    
    gpuContext gpu    = amoebaGpu->gpuContext;
    int numThreads    = min(256, (gpu->natoms));
    int numBlocks     =  1 + (gpu->natoms/numThreads);

//#ifdef AMOEBA_DEBUG
#if 0
    if( amoebaGpu->log ){
        (void) fprintf( amoebaGpu->log, "%s: numBlocks=%d numThreads=%d\n", methodName, numBlocks, numThreads ); (void) fflush( amoebaGpu->log );
    }
    amoebaGpu->psForce->Download();
    amoebaGpu->psTorque->Download();
    int maxPrint        = 20;
    (void) fprintf( amoebaGpu->log,"Pre torqueMap\n" );
    for( int ii = 0; ii < gpu->natoms; ii++ ){
        (void) fprintf( amoebaGpu->log, "%5d ", ii);

        int indexOffset     = ii*3;

        (void) fprintf( amoebaGpu->log,"E[%16.9e %16.9e %16.9e] ",
                        amoebaGpu->psForce->_pSysStream[0][indexOffset],
                        amoebaGpu->psForce->_pSysStream[0][indexOffset+1],
                        amoebaGpu->psForce->_pSysStream[0][indexOffset+2] );
        (void) fprintf( amoebaGpu->log,"T[%16.9e %16.9e %16.9e]\n",
                        amoebaGpu->psTorque->_pSysStream[0][indexOffset],
                        amoebaGpu->psTorque->_pSysStream[0][indexOffset+1],
                        amoebaGpu->psTorque->_pSysStream[0][indexOffset+2] );
        if( ii == maxPrint && (gpu->natoms - maxPrint) > ii )ii = gpu->natoms - maxPrint;
    }
    int nansDetected   = checkForNansAndInfinities( gpu->natoms*3, amoebaGpu->psForce );
    nansDetected      += checkForNansAndInfinities( gpu->natoms*3, amoebaGpu->psTorque );
    if( nansDetected ){
        (void) fprintf( amoebaGpu->log,"WARNING: %d nans/infinities detected force/torques.\n", nansDetected );
    } else {
        (void) fprintf( amoebaGpu->log,"No nans/infinities detected in force/torques.\n" );
    }

// zero forces
#if 0
    for( int ii = 0; ii < 3*gpu->natoms; ii++ ){
        amoebaGpu->psForce->_pSysStream[0][ii] = 0.0f;
    }
    amoebaGpu->psForce->Upload();
#endif

#endif

    // torqueMapForce is zeroed when initialized; should not need to be reinitialized 

/*
    AmoebaTorqueMapZeroKernel<<< numBlocks, numThreads >>>( 
       gpu->natoms, amoebaGpu->torqueMapForce->_pDevStream[0] );
    LAUNCHERROR("AmoebaMapTrqZeroKernel");  
*/


    amoebaMapTorqueToForce_kernel<<< numBlocks, numThreads>>> (
                gpu->natoms,
                gpu->psPosq4->_pDevStream[0],
                psTorque->_pDevStream[0],
                amoebaGpu->psMultipoleParticlesIdsAndAxisType->_pDevStream[0],
                amoebaGpu->maxMapTorqueDifference,
                amoebaGpu->torqueMapForce->_pDevStream[0] );
    LAUNCHERROR("AmoebaMapTrqKernel");
  
//#ifdef AMOEBA_DEBUG
#if 0
    amoebaGpu->torqueMapForce->Download();
    //int maxPrint        = 10;
    (void) fprintf( amoebaGpu->log,"Post AmoebaMapTrqKernel maxMapTorqueDifference=%d\n", amoebaGpu->maxMapTorqueDifference );
    for( int ii = 0; ii < gpu->natoms; ii++ ){
        (void) fprintf( amoebaGpu->log, "\n%5d multi[%d %d %d %d]\n", ii, 
                        amoebaGpu->psMultipoleParticlesIdsAndAxisType->_pSysStream[0][ii].x,
                        amoebaGpu->psMultipoleParticlesIdsAndAxisType->_pSysStream[0][ii].y,
                        amoebaGpu->psMultipoleParticlesIdsAndAxisType->_pSysStream[0][ii].z,
                        amoebaGpu->psMultipoleParticlesIdsAndAxisType->_pSysStream[0][ii].w );

        int indexOffset     = ii*3*amoebaGpu->maxMapTorqueDifference;
        float sum[3]        = { 0.0f,  0.0f,  0.0f };
        for( int jj = 0; jj < amoebaGpu->maxMapTorqueDifference; jj++ ){
            if( amoebaGpu->torqueMapForce->_pSysStream[0][indexOffset] != 0.0f ){ 
                (void) fprintf( amoebaGpu->log,"    %4d %4d Temp[%16.9e %16.9e %16.9e] %d\n",
                                ii, jj + amoebaGpu->psMultipoleParticlesIdsAndAxisType->_pSysStream[0][ii].z,
                                amoebaGpu->torqueMapForce->_pSysStream[0][indexOffset],
                                amoebaGpu->torqueMapForce->_pSysStream[0][indexOffset+1],
                                amoebaGpu->torqueMapForce->_pSysStream[0][indexOffset+2], indexOffset );
                sum[0] += amoebaGpu->torqueMapForce->_pSysStream[0][indexOffset];
                sum[1] += amoebaGpu->torqueMapForce->_pSysStream[0][indexOffset+1];
                sum[2] += amoebaGpu->torqueMapForce->_pSysStream[0][indexOffset+2];
            }
            indexOffset += 3;
        }
        (void) fprintf( amoebaGpu->log,"               Sum[%16.9e %16.9e %16.9e]\n", sum[0], sum[1], sum[2] );
        if( ii == maxPrint && (gpu->natoms - maxPrint) > ii )ii = gpu->natoms - maxPrint;
    }
#endif

    numBlocks  = gpu->natoms;
    numThreads = amoebaGpu->maxMapTorqueDifferencePow2;

    amoebaMapTorqueReduce_kernel2<<< numBlocks, numThreads>>>(
            numThreads, gpu->natoms,
            amoebaGpu->maxMapTorqueDifference,
            amoebaGpu->torqueMapForce->_pDevStream[0],
            psForce->_pDevStream[0], psCudaForce4->_pDevStream[0] );
    LAUNCHERROR("amoebaMapTorqueReduce_kernel2");

#ifdef AMOEBA_DEBUG
    if( amoebaGpu->log ){
        (void) fprintf( amoebaGpu->log, "%s: numBlocks=%d numThreads=%d %d\n", methodName, numBlocks, numThreads, amoebaGpu->maxMapTorqueDifferencePow2); (void) fflush( amoebaGpu->log );
        amoebaGpu->psForce->Download();
        psCudaForce4->Download();
 
        amoebaGpu->psTorque->Download();
        int maxPrint        = 10;
        (void) fprintf( amoebaGpu->log,"Post torqueMap\n" );
        for( int ii = 0; ii < gpu->natoms; ii++ ){
            (void) fprintf( amoebaGpu->log, "%5d ", ii);
    
            int indexOffset     = ii*3;
    
            (void) fprintf( amoebaGpu->log,"FTtl[%16.9e %16.9e %16.9e] ",
                            psCudaForce4->_pSysStream[0][ii].x,
                            psCudaForce4->_pSysStream[0][ii].y,
                            psCudaForce4->_pSysStream[0][ii].z );
            (void) fprintf( amoebaGpu->log,"F[%16.9e %16.9e %16.9e] ",
                            amoebaGpu->psForce->_pSysStream[0][indexOffset],
                            amoebaGpu->psForce->_pSysStream[0][indexOffset+1],
                            amoebaGpu->psForce->_pSysStream[0][indexOffset+2] );
            (void) fprintf( amoebaGpu->log,"T[%16.9e %16.9e %16.9e]\n",
                            amoebaGpu->psTorque->_pSysStream[0][indexOffset],
                            amoebaGpu->psTorque->_pSysStream[0][indexOffset+1],
                            amoebaGpu->psTorque->_pSysStream[0][indexOffset+2] );
            if( ii == maxPrint && (gpu->natoms - maxPrint) > ii )ii = gpu->natoms - maxPrint;
        }
        (void) fflush( amoebaGpu->log );
    }
    if( 1 ){
        //std::vector<int> fileId;
        //fileId.push_back( 0 );
        VectorOfDoubleVectors outputVector;
        //cudaLoadCudaFloat4Array( gpu->natoms, 3, gpu->psPosq4,              outputVector );
        cudaLoadCudaFloat4Array( gpu->natoms, 4, gpu->psForce4,             outputVector );
        cudaLoadCudaFloatArray( gpu->natoms,  3, amoebaGpu->psForce,        outputVector );
        cudaLoadCudaFloatArray( gpu->natoms,  3, amoebaGpu->psTorque,       outputVector);
        cudaWriteVectorOfDoubleVectorsToFile( "CudaVacuumElecForce", fileId, outputVector );
    }
#endif

}
