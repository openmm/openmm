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

__device__ static void crossVector3( float* vector1, float* vector2, float* vector3 )
{

    vector3[0]                    = vector1[1]*vector2[2] - vector1[2]*vector2[1];
    vector3[1]                    = vector1[2]*vector2[0] - vector1[0]*vector2[2];
    vector3[2]                    = vector1[0]*vector2[1] - vector1[1]*vector2[0];

}

__global__
#if (__CUDA_ARCH__ >= 200)
__launch_bounds__(GF1XX_THREADS_PER_BLOCK, 1)
#elif (__CUDA_ARCH__ >= 120)
__launch_bounds__(GT2XX_THREADS_PER_BLOCK, 1)
#else
__launch_bounds__(G8X_THREADS_PER_BLOCK, 1)
#endif
void amoebaMapTorqueToForce_kernel( float* torque, int maxDiff, float* tempElecForce){

  // ---------------------------------------------------------------------------------------
  
    int ii;
    int threadId                 = __mul24(blockIdx.x,blockDim.x) + threadIdx.x ;
    int numOfAtoms               = cSim.atoms;
    if( threadId >= numOfAtoms )return;
    
    float4* atomCoord            = cSim.pPosq;
    int4* multiPoleAtoms         = cAmoebaSim.pMultipoleParticlesIdsAndAxisType;
    int* multiPoleAtomOffset     = cAmoebaSim.pMultipoleAxisOffset;
  
    const int U                  = 0;
    const int V                  = 1;
    const int W                  = 2;
    const int R                  = 3;
    const int S                  = 4;
    const int UV                 = 5;
    const int UW                 = 6;
    const int VW                 = 7;
    const int UR                 = 8;
    const int US                 = 9;
    const int VS                 = 10;
    const int WS                 = 11;
    const int LastVectorIndex    = 12;
    
    const int X                  = 0;
    const int Y                  = 1;
    const int Z                  = 2;
    const int I                  = 3;
    
    float forces[4][3];
    float norms[LastVectorIndex];
    float vector[LastVectorIndex][3];
    float angles[LastVectorIndex][2];
  
  // ---------------------------------------------------------------------------------------
  
    int axisAtom                  = multiPoleAtoms[threadId].z;
    int axisType                  = multiPoleAtoms[threadId].w;

    vector[U][0]                  = atomCoord[threadId].x - atomCoord[axisAtom].x;
    vector[U][1]                  = atomCoord[threadId].y - atomCoord[axisAtom].y;
    vector[U][2]                  = atomCoord[threadId].z - atomCoord[axisAtom].z;

    norms[U]                      = normVector3( vector[U] );

    if( axisType != 4 ){

        axisAtom                  = multiPoleAtoms[threadId].x;
        vector[V][0]              = atomCoord[threadId].x - atomCoord[axisAtom].x;
        vector[V][1]              = atomCoord[threadId].y - atomCoord[axisAtom].y;
        vector[V][2]              = atomCoord[threadId].z - atomCoord[axisAtom].z;

    } else {
        vector[V][0]              = 0.1f;
        vector[V][1]              = 0.1f;
        vector[V][2]              = 0.1f;
    }

    norms[V]                      = normVector3( vector[V] );

    // W = UxV

    if( axisType < 2 || axisType > 3  ){
        crossVector3( vector[U], vector[V], vector[W] );
    } else { 
        axisAtom                  = multiPoleAtoms[threadId].y;
    
        vector[W][0]              = atomCoord[threadId].x - atomCoord[axisAtom].x;
        vector[W][1]              = atomCoord[threadId].y - atomCoord[axisAtom].y;
        vector[W][2]              = atomCoord[threadId].z - atomCoord[axisAtom].z;
    } 
    norms[W]                      = normVector3( vector[W] );

    crossVector3( vector[V], vector[U], vector[UV] );
    crossVector3( vector[W], vector[U], vector[UW] );
    crossVector3( vector[W], vector[V], vector[VW] );

    norms[UV]                     = normVector3( vector[UV] );
    norms[UW]                     = normVector3( vector[UW] );
    norms[VW]                     = normVector3( vector[VW] );

    angles[UV][0]                 = DOT3( vector[U], vector[V] );
    angles[UV][1]                 = sqrtf( 1.0f - angles[UV][0]*angles[UV][0]);

    angles[UW][0]                 = DOT3( vector[U], vector[W] );
    angles[UW][1]                 = sqrtf( 1.0f - angles[UW][0]*angles[UW][0]);

    angles[VW][0]                 = DOT3( vector[V], vector[W] );
    angles[VW][1]                 = sqrtf( 1.0f - angles[VW][0]*angles[VW][0]);

    float dphi[3];
    dphi[U]                       = DOT3( vector[U], (torque + threadId*3) );
    dphi[V]                       = DOT3( vector[V], (torque + threadId*3) );
    dphi[W]                       = DOT3( vector[W], (torque + threadId*3) );

    dphi[U]                      *= -1.0f;
    dphi[V]                      *= -1.0f;
    dphi[W]                      *= -1.0f;

    // z-then-x and bisector

    if( axisType == 0 || axisType == 1 ){

        float factor1;
        float factor2;
        float factor3;
        float factor4;
    
        factor1                 =  dphi[V]/(norms[U]*angles[UV][1]);
        factor2                 =  dphi[W]/(norms[U]);
        factor3                 = -dphi[U]/(norms[V]*angles[UV][1]);
    
        if( axisType == 1 ){
            factor2    *= 0.5f;
            factor4     = 0.5f*dphi[W]/(norms[V]);
        } else {
            factor4     = 0.0f;
        }
 
        for( ii = 0; ii < 3; ii++ ){
            forces[Z][ii]                  =  vector[UV][ii]*factor1 + factor2*vector[UW][ii];
            forces[X][ii]                  =  vector[UV][ii]*factor3 + factor4*vector[VW][ii];
            forces[I][ii]                  = -(forces[X][ii] + forces[Z][ii]);
            forces[Y][ii]                  = 0.0f;
        }

    } else if( axisType == 2 ){

        // z-bisect

        vector[R][0] = vector[V][0] + vector[W][0]; 
        vector[R][1] = vector[V][1] + vector[W][1]; 
        vector[R][2] = vector[V][2] + vector[W][2]; 

        crossVector3( vector[U], vector[R], vector[S] );

        norms[R]     = normVector3( vector[R] );
        norms[S]     = normVector3( vector[S] );

        crossVector3( vector[R], vector[U], vector[UR] );
        crossVector3( vector[S], vector[U], vector[US] );
        crossVector3( vector[S], vector[V], vector[VS] );
        crossVector3( vector[S], vector[W], vector[WS] );

        norms[UR]     = normVector3( vector[UR] );
        norms[US]     = normVector3( vector[US] );
        norms[VS]     = normVector3( vector[VS] );
        norms[WS]     = normVector3( vector[WS] );

        angles[UR][0]             = DOT3( vector[U], vector[R] );
        angles[UR][1]             = sqrtf( 1.0f - angles[UR][0]*angles[UR][0]);

        angles[US][0]             = DOT3( vector[U], vector[S] );
        angles[US][1]             = sqrtf( 1.0f - angles[US][0]*angles[US][0]);

        angles[VS][0]             = DOT3( vector[V], vector[S] );
        angles[VS][1]             = sqrtf( 1.0f - angles[VS][0]*angles[VS][0]);

        angles[WS][0]             = DOT3( vector[W], vector[S] );
        angles[WS][1]             = sqrtf( 1.0f - angles[WS][0]*angles[WS][0]);
 
        float t1[3];
        float t2[3];
        t1[0]                     = vector[V][0] - vector[S][0]*angles[VS][0];
        t1[1]                     = vector[V][1] - vector[S][1]*angles[VS][0];
        t1[2]                     = vector[V][2] - vector[S][2]*angles[VS][0];

        t2[0]                     = vector[W][0] - vector[S][0]*angles[WS][0];
        t2[1]                     = vector[W][1] - vector[S][1]*angles[WS][0];
        t2[2]                     = vector[W][2] - vector[S][2]*angles[WS][0];
        float notUsed             = normVector3( t1 );
              notUsed             = normVector3( t2 );
        float ut1cos              = DOT3( vector[U], t1 );
        float ut1sin              = sqrtf( 1.0f - ut1cos*ut1cos);
        float ut2cos              = DOT3( vector[U], t2 );
        float ut2sin              = sqrtf( 1.0f - ut2cos*ut2cos);

        float dphiR               = -1.0f*DOT3( vector[R], (torque + threadId*3) );
        float dphiS               = -1.0f*DOT3( vector[S], (torque + threadId*3) );

        float factor1             = dphiR/(norms[U]*angles[UR][1]);
        float factor2             = dphiS/(norms[U]);
        float factor3             = dphi[U]/(norms[V]*(ut1sin+ut2sin));
        float factor4             = dphi[U]/(norms[W]*(ut1sin+ut2sin));
        for( ii = 0; ii < 3; ii++ ){
            forces[Z][ii]         =  vector[UR][ii]*factor1 + factor2*vector[US][ii];
            forces[X][ii]         = (angles[VS][1]*vector[S][ii] - angles[VS][0]*t1[ii])*factor3;
            forces[Y][ii]         = (angles[WS][1]*vector[S][ii] - angles[WS][0]*t2[ii])*factor4;
            forces[I][ii]         = -(forces[X][ii] + forces[Y][ii] + forces[Z][ii]);
        }
 
    } else if( axisType == 3 ){

        // 3-fold

        for( ii = 0; ii < 3; ii++ ){
            float du =  vector[UW][ii]*dphi[W]/(norms[U]*angles[UW][1]) +
                        vector[UV][ii]*dphi[V]/(norms[U]*angles[UV][1]) -
                        vector[UW][ii]*dphi[U]/(norms[U]*angles[UW][1]) -
                        vector[UV][ii]*dphi[U]/(norms[U]*angles[UV][1]);

            float dv =  vector[VW][ii]*dphi[W]/(norms[V]*angles[VW][1]) -
                        vector[UV][ii]*dphi[U]/(norms[V]*angles[UV][1]) -
                        vector[VW][ii]*dphi[V]/(norms[V]*angles[VW][1]) +
                        vector[UV][ii]*dphi[V]/(norms[V]*angles[UV][1]);

            float dw = -vector[UW][ii]*dphi[U]/(norms[W]*angles[UW][1]) -
                        vector[VW][ii]*dphi[V]/(norms[W]*angles[VW][1]) +
                        vector[UW][ii]*dphi[W]/(norms[W]*angles[UW][1]) +
                        vector[VW][ii]*dphi[W]/(norms[W]*angles[VW][1]);

            du      /= 3.0f;
            dv      /= 3.0f;
            dw      /= 3.0f;

            forces[Z][ii]  = du;
            forces[X][ii]  = dv;
            forces[Y][ii]  = dw;
            forces[I][ii]  = -(du + dv + dw);
        } 

    } else if( axisType == 4 ){

        // z-only

        for( ii = 0; ii < 3; ii++ ){
            float du              = vector[UV][ii]*dphi[V]/(norms[U]*angles[UV][1]);
            forces[Z][ii]         = du;
            forces[X][ii]         = 0.0f;
            forces[Y][ii]         = 0.0f;
            forces[I][ii]         = -du;
        }
    } else {

        for( ii = 0; ii < 3; ii++ ){
            forces[Z][ii]         = 0.0f;
            forces[X][ii]         = 0.0f;
            forces[Y][ii]         = 0.0f;
            forces[I][ii]         = 0.0f;
        }
    }

    // load results

    int temp                      = multiPoleAtoms[threadId].z;
    int min                       = multiPoleAtomOffset[temp];
    int offset                    = 3*(temp*maxDiff + threadId-min);
    tempElecForce[offset  ]       = forces[Z][0];
    tempElecForce[offset+1]       = forces[Z][1];
    tempElecForce[offset+2]       = forces[Z][2];
   
   
    if( axisType != 4 ){
        temp                          = multiPoleAtoms[threadId].x;
        min                           = multiPoleAtomOffset[temp];
        offset                        = 3*(temp*maxDiff + threadId-min);
        tempElecForce[offset  ]       = forces[X][0];
        tempElecForce[offset+1]       = forces[X][1];
        tempElecForce[offset+2]       = forces[X][2];
    }
 
    if( axisType == 2 || axisType == 3 ){
        temp                          = multiPoleAtoms[threadId].y;
        if( temp > -1 ){
            min                           = multiPoleAtomOffset[temp];
            offset                        = 3*(temp*maxDiff + threadId-min);
            tempElecForce[offset  ]       = forces[Y][0];
            tempElecForce[offset+1]       = forces[Y][1];
            tempElecForce[offset+2]       = forces[Y][2];
        }
    }
 
    min                           = multiPoleAtomOffset[threadId];
    offset                        = 3*(threadId*(maxDiff + 1)-min); 
    tempElecForce[offset  ]       = forces[I][0];
    tempElecForce[offset+1]       = forces[I][1];
    tempElecForce[offset+2]       = forces[I][2];
   
}

__global__
#if (__CUDA_ARCH__ >= 200)
__launch_bounds__(GF1XX_THREADS_PER_BLOCK, 1)
#elif (__CUDA_ARCH__ >= 120)
__launch_bounds__(GT2XX_THREADS_PER_BLOCK, 1)
#else
__launch_bounds__(G8X_THREADS_PER_BLOCK, 1)
#endif
void amoebaMapTorqueToForceOld_kernel( float* torque, int maxDiff, float* tempElecForce){

  // ---------------------------------------------------------------------------------------
  
    int threadId                 = __mul24(blockIdx.x,blockDim.x) + threadIdx.x ;
    int numOfAtoms               = cSim.atoms;
    if(threadId >= numOfAtoms)return;
    
    float4* atomCoord            = cSim.pPosq;
    int4* multiPoleAtoms         = cAmoebaSim.pMultipoleParticlesIdsAndAxisType;
    int* multiPoleAtomOffset     = cAmoebaSim.pMultipoleAxisOffset;
  
    const int U                  = 0;
    const int V                  = 1;
    const int W                  = 2;
    const int R                  = 3;
    const int S                  = 4;
    
    const int X                  = 0;
    const int Y                  = 1;
    const int Z                  = 2;
    
    float forces[3][3];
    float norms[5];
    float vector[5][3];
  
  // ---------------------------------------------------------------------------------------
  
    int axisAtom                  = multiPoleAtoms[threadId].z;
    int axisType                  = multiPoleAtoms[threadId].w;

    vector[U][0]                  = atomCoord[axisAtom].x - atomCoord[threadId].x;
    vector[U][1]                  = atomCoord[axisAtom].y - atomCoord[threadId].y;
    vector[U][2]                  = atomCoord[axisAtom].z - atomCoord[threadId].z;

    norms[U]                      = normVector3( vector[U] );

    if( axisType != 4 ){
        axisAtom                      = multiPoleAtoms[threadId].x;
    
        vector[V][0]                  = atomCoord[axisAtom].x - atomCoord[threadId].x;
        vector[V][1]                  = atomCoord[axisAtom].y - atomCoord[threadId].y;
        vector[V][2]                  = atomCoord[axisAtom].z - atomCoord[threadId].z;
    } else {
        vector[V][0]                  = 0.1f;
        vector[V][1]                  = 0.1f;
        vector[V][2]                  = 0.1f;
    }

    norms[V]                      = normVector3( vector[V] );

    // W = UxV

    if( axisType < 2 || axisType > 3  ){
        vector[W][0]                  = vector[U][1]*vector[V][2] - vector[U][2]*vector[V][1];
        vector[W][1]                  = vector[U][2]*vector[V][0] - vector[U][0]*vector[V][2];
        vector[W][2]                  = vector[U][0]*vector[V][1] - vector[U][1]*vector[V][0];
   } else { 
        axisAtom                      = multiPoleAtoms[threadId].y;
    
        vector[W][0]                  = atomCoord[axisAtom].x - atomCoord[threadId].x;
        vector[W][1]                  = atomCoord[axisAtom].y - atomCoord[threadId].y;
        vector[W][2]                  = atomCoord[axisAtom].z - atomCoord[threadId].z;
    } 
    norms[W]                      = normVector3( vector[W] );

    if( axisType == 2 ){

        vector[R][0] = vector[V][0] + vector[W][0]; 
        vector[R][1] = vector[V][1] + vector[W][1]; 
        vector[R][2] = vector[V][2] + vector[W][2]; 

        crossVector3( vector[U], vector[R], vector[S] );
        norms[R]     = normVector3( vector[R] );
        norms[S]     = normVector3( vector[S] );
    }


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
    dphi[0]                       = vector[U][0]*torque[threadId*3] + vector[U][1]*torque[threadId*3+1] + vector[U][2]*torque[threadId*3+2];
    dphi[1]                       = vector[V][0]*torque[threadId*3] + vector[V][1]*torque[threadId*3+1] + vector[V][2]*torque[threadId*3+2];
    dphi[2]                       = vector[W][0]*torque[threadId*3] + vector[W][1]*torque[threadId*3+1] + vector[W][2]*torque[threadId*3+2];

    // clamp c to interval [-1,1]

    float c                       = DOT3( vector[U], vector[V] );
          c                       = c >  1.0f ?  1.0f : c;
          c                       = c < -1.0f ? -1.0f : c;

    float s                       = SQRT( 1.0f - (c*c) );
    float uvdis                   = norms[U]*s;
    float vudis                   = norms[V]*s;

    float factorX;
    float factorZ;
    if( axisType == 1 ){       
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

    int temp                      = multiPoleAtoms[threadId].z;
    int min                       = multiPoleAtomOffset[temp];
    int offset                    = 3*(temp*maxDiff + threadId-min);
    tempElecForce[offset  ]       = forces[X][0];
    tempElecForce[offset+1]       = forces[X][1];
    tempElecForce[offset+2]       = forces[X][2];
   
    temp                          = multiPoleAtoms[threadId].x;
    min                           = multiPoleAtomOffset[temp];
    offset                        = 3*(temp*maxDiff + threadId-min);
    tempElecForce[offset  ]       = forces[Z][0];
    tempElecForce[offset+1]       = forces[Z][1];
    tempElecForce[offset+2]       = forces[Z][2];
 
    min                           = multiPoleAtomOffset[threadId];
    offset                        = 3*(threadId*(maxDiff + 1)-min); 
    tempElecForce[offset  ]       = forces[Y][0];
    tempElecForce[offset+1]       = forces[Y][1];
    tempElecForce[offset+2]       = forces[Y][2];
   
}

__global__
#if (__CUDA_ARCH__ >= 200)
__launch_bounds__(GF1XX_THREADS_PER_BLOCK, 1)
#elif (__CUDA_ARCH__ >= 120)
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
#elif (__CUDA_ARCH__ >= 120)
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

__global__
#if (__CUDA_ARCH__ >= 200)
__launch_bounds__(GF1XX_THREADS_PER_BLOCK, 1)
#elif (__CUDA_ARCH__ >= 120)
__launch_bounds__(GT2XX_THREADS_PER_BLOCK, 1)
#else
__launch_bounds__(G8X_THREADS_PER_BLOCK, 1)
#endif
void amoebaMapTorqueReduce_kernel3(
					   int numThreads,
					   int numOfAtoms,
					   int maxDiff,
					   float* tempElecForce,
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
        outputForce[blockIdx.x].x  += sfx[0];
        outputForce[blockIdx.x].y  += sfy[0];
        outputForce[blockIdx.x].z  += sfz[0];
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
                        amoebaGpu->psForce->_pSysData[indexOffset],
                        amoebaGpu->psForce->_pSysData[indexOffset+1],
                        amoebaGpu->psForce->_pSysData[indexOffset+2] );
        (void) fprintf( amoebaGpu->log,"T[%16.9e %16.9e %16.9e]\n",
                        amoebaGpu->psTorque->_pSysData[indexOffset],
                        amoebaGpu->psTorque->_pSysData[indexOffset+1],
                        amoebaGpu->psTorque->_pSysData[indexOffset+2] );
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
        amoebaGpu->psForce->_pSysData[ii] = 0.0f;
    }
    amoebaGpu->psForce->Upload();
#endif

#endif

    // torqueMapForce is zeroed when initialized; should not need to be reinitialized 

/*
    AmoebaTorqueMapZeroKernel<<< numBlocks, numThreads >>>( 
       gpu->natoms, amoebaGpu->torqueMapForce->_pDevData );
    LAUNCHERROR("AmoebaMapTrqZeroKernel");  
*/


    amoebaMapTorqueToForce_kernel<<< numBlocks, numThreads>>> (
                psTorque->_pDevData,
                amoebaGpu->maxMapTorqueDifference,
                amoebaGpu->torqueMapForce->_pDevData );
    LAUNCHERROR("AmoebaMapTrqKernel");
  
//#ifdef AMOEBA_DEBUG
#if 0
    amoebaGpu->torqueMapForce->Download();
    //int maxPrint        = 10;
    (void) fprintf( amoebaGpu->log,"Post AmoebaMapTrqKernel maxMapTorqueDifference=%d\n", amoebaGpu->maxMapTorqueDifference );
    for( int ii = 0; ii < gpu->natoms; ii++ ){
        (void) fprintf( amoebaGpu->log, "\n%5d multi[%d %d %d %d] offset=%d\n", ii, 
                        amoebaGpu->psMultipoleParticlesIdsAndAxisType->_pSysData[ii].x,
                        amoebaGpu->psMultipoleParticlesIdsAndAxisType->_pSysData[ii].y,
                        amoebaGpu->psMultipoleParticlesIdsAndAxisType->_pSysData[ii].z,
                        amoebaGpu->psMultipoleParticlesIdsAndAxisType->_pSysData[ii].w,
                        amoebaGpu->psMultipoleAxisOffset->_pSysData[ii] );

        int indexOffset     = ii*3*amoebaGpu->maxMapTorqueDifference;
        float sum[3]        = { 0.0f,  0.0f,  0.0f };
        for( int jj = 0; jj < amoebaGpu->maxMapTorqueDifference; jj++ ){
            if( amoebaGpu->torqueMapForce->_pSysData[indexOffset] != 0.0f ){ 
                (void) fprintf( amoebaGpu->log,"    %4d %4d Temp[%16.9e %16.9e %16.9e] %d\n",
                                ii, jj + amoebaGpu->psMultipoleAxisOffset->_pSysData[ii],
                                amoebaGpu->torqueMapForce->_pSysData[indexOffset],
                                amoebaGpu->torqueMapForce->_pSysData[indexOffset+1],
                                amoebaGpu->torqueMapForce->_pSysData[indexOffset+2], indexOffset );
                sum[0] += amoebaGpu->torqueMapForce->_pSysData[indexOffset];
                sum[1] += amoebaGpu->torqueMapForce->_pSysData[indexOffset+1];
                sum[2] += amoebaGpu->torqueMapForce->_pSysData[indexOffset+2];
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
            amoebaGpu->torqueMapForce->_pDevData,
            psForce->_pDevData );
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
                            amoebaGpu->psForce->_pSysData[indexOffset],
                            amoebaGpu->psForce->_pSysData[indexOffset+1],
                            amoebaGpu->psForce->_pSysData[indexOffset+2] );
            (void) fprintf( amoebaGpu->log,"T[%16.9e %16.9e %16.9e]\n",
                            amoebaGpu->psTorque->_pSysData[indexOffset],
                            amoebaGpu->psTorque->_pSysData[indexOffset+1],
                            amoebaGpu->psTorque->_pSysData[indexOffset+2] );
            if( ii == maxPrint && (gpu->natoms - maxPrint) > ii )ii = gpu->natoms - maxPrint;
        }
        (void) fflush( amoebaGpu->log );
    }
    if( 1 ){
        //std::vector<int> fileId;
        //fileId.push_back( 0 );
        VectorOfDoubleVectors outputVector;
        cudaLoadCudaFloat4Array( gpu->natoms, 3, gpu->psPosq4,              outputVector );
        cudaLoadCudaFloatArray( gpu->natoms,  3, amoebaGpu->psForce,      outputVector, NULL );
        cudaLoadCudaFloatArray( gpu->natoms,  3, amoebaGpu->psTorque, outputVector, NULL);
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
    psForce->Download();
    psTorque->Download();
    int maxPrint        = 20;
    (void) fprintf( amoebaGpu->log,"Pre torqueMap\n" );
    for( int ii = 0; ii < gpu->natoms; ii++ ){
        (void) fprintf( amoebaGpu->log, "%5d ", ii);

        int indexOffset     = ii*3;

        (void) fprintf( amoebaGpu->log,"E[%16.9e %16.9e %16.9e] ",
                        psForce->_pSysData[indexOffset],
                        psForce->_pSysData[indexOffset+1],
                        psForce->_pSysData[indexOffset+2] );
        (void) fprintf( amoebaGpu->log,"T[%16.9e %16.9e %16.9e]\n",
                        psTorque->_pSysData[indexOffset],
                        psTorque->_pSysData[indexOffset+1],
                        psTorque->_pSysData[indexOffset+2] );
        if( ii == maxPrint && (gpu->natoms - maxPrint) > ii )ii = gpu->natoms - maxPrint;
    }
    int nansDetected   = checkForNansAndInfinities( gpu->natoms*3, psForce );
    nansDetected      += checkForNansAndInfinities( gpu->natoms*3, psTorque );
    if( nansDetected ){
        (void) fprintf( amoebaGpu->log,"WARNING: %d nans/infinities detected force/torques.\n", nansDetected );
    } else {
        (void) fprintf( amoebaGpu->log,"No nans/infinities detected in force/torques.\n" );
    }

// zero forces
#if 0
    for( int ii = 0; ii < 3*gpu->natoms; ii++ ){
        psForce->_pSysData[ii] = 0.0f;
    }
    psForce->Upload();
#endif

#if 0
    (void) fprintf( amoebaGpu->log,"Setting force & torque values.\n" );
    for( int ii = 0; ii < 3*gpu->natoms; ii += 3 ){

        psTorque->_pSysData[ii]   = 1.0f;
        psTorque->_pSysData[ii+1] = 0.0f;
        psTorque->_pSysData[ii+2] = 0.0f;

        psForce->_pSysData[ii]    = 0.0f;
        psForce->_pSysData[ii+1]  = 0.0f;
        psForce->_pSysData[ii+2]  = 0.0f;
    }
    for( int ii = 0; ii < gpu->natoms; ii++ ){
        psCudaForce4->_pSysData[ii].x        = 0.0f;
        psCudaForce4->_pSysData[ii].y        = 0.0f;
        psCudaForce4->_pSysData[ii].z        = 0.0f;
    }
    psForce->Upload();
    psTorque->Upload();
    psCudaForce4->Upload();
#endif

#endif

    // torqueMapForce is zeroed when initialized; should not need to be reinitialized 

/*
    AmoebaTorqueMapZeroKernel<<< numBlocks, numThreads >>>( 
       gpu->natoms, amoebaGpu->torqueMapForce->_pDevData );
    LAUNCHERROR("AmoebaMapTrqZeroKernel");  
*/


    //amoebaMapTorqueToForceOld_kernel<<< numBlocks, numThreads>>> (
    amoebaMapTorqueToForce_kernel<<< numBlocks, numThreads>>> (
                psTorque->_pDevData,
                amoebaGpu->maxMapTorqueDifference,
                amoebaGpu->torqueMapForce->_pDevData );
    LAUNCHERROR("AmoebaMapTrqKernel");
  
//#ifdef AMOEBA_DEBUG
#if 0
    amoebaGpu->torqueMapForce->Download();
    //int maxPrint        = 10;
    (void) fprintf( amoebaGpu->log,"Post AmoebaMapTrqKernel maxMapTorqueDifference=%d\n", amoebaGpu->maxMapTorqueDifference );
    for( int ii = 0; ii < gpu->natoms; ii++ ){
        (void) fprintf( amoebaGpu->log, "\n%5d multi[%d %d %d %d]\n", ii, 
                        amoebaGpu->psMultipoleParticlesIdsAndAxisType->_pSysData[ii].x,
                        amoebaGpu->psMultipoleParticlesIdsAndAxisType->_pSysData[ii].y,
                        amoebaGpu->psMultipoleParticlesIdsAndAxisType->_pSysData[ii].z,
                        amoebaGpu->psMultipoleParticlesIdsAndAxisType->_pSysData[ii].w, amoebaGpu->psMultipoleAxisOffset->_pSysData[ii] );

        int indexOffset     = ii*3*amoebaGpu->maxMapTorqueDifference;
        float sum[3]        = { 0.0f,  0.0f,  0.0f };
        for( int jj = 0; jj < amoebaGpu->maxMapTorqueDifference; jj++ ){
            if( amoebaGpu->torqueMapForce->_pSysData[indexOffset] != 0.0f ){ 
                (void) fprintf( amoebaGpu->log,"    %4d %4d Temp[%16.9e %16.9e %16.9e] %d\n",
                                ii, jj + amoebaGpu->psMultipoleAxisOffset->_pSysData[ii],
                                amoebaGpu->torqueMapForce->_pSysData[indexOffset],
                                amoebaGpu->torqueMapForce->_pSysData[indexOffset+1],
                                amoebaGpu->torqueMapForce->_pSysData[indexOffset+2], indexOffset );
                sum[0] += amoebaGpu->torqueMapForce->_pSysData[indexOffset];
                sum[1] += amoebaGpu->torqueMapForce->_pSysData[indexOffset+1];
                sum[2] += amoebaGpu->torqueMapForce->_pSysData[indexOffset+2];
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
            amoebaGpu->torqueMapForce->_pDevData,
            psForce->_pDevData, psCudaForce4->_pDevData );
    LAUNCHERROR("amoebaMapTorqueReduce_kernel2");

#ifdef AMOEBA_DEBUG
    if( amoebaGpu->log ){
        (void) fprintf( amoebaGpu->log, "%s: numBlocks=%d numThreads=%d %d\n", methodName, numBlocks, numThreads, amoebaGpu->maxMapTorqueDifferencePow2); (void) fflush( amoebaGpu->log );
        psForce->Download();
        psCudaForce4->Download();
        amoebaGpu->torqueMapForce->Download();
        psTorque->Download();
        int maxPrint        = 10;
        (void) fprintf( amoebaGpu->log,"Post torqueMap\n" );
        for( int ii = 0; ii < gpu->natoms; ii++ ){
            (void) fprintf( amoebaGpu->log, "%5d ", ii);
    
            int indexOffset     = ii*3;
    
            (void) fprintf( amoebaGpu->log,"FTtl[%16.9e %16.9e %16.9e] ",
                            psCudaForce4->_pSysData[ii].x,
                            psCudaForce4->_pSysData[ii].y,
                            psCudaForce4->_pSysData[ii].z );
            (void) fprintf( amoebaGpu->log,"F[%16.9e %16.9e %16.9e] ",
                            psForce->_pSysData[indexOffset],
                            psForce->_pSysData[indexOffset+1],
                            psForce->_pSysData[indexOffset+2] );
            (void) fprintf( amoebaGpu->log,"fT[%16.9e %16.9e %16.9e] ",
                            amoebaGpu->torqueMapForce->_pSysData[indexOffset],
                            amoebaGpu->torqueMapForce->_pSysData[indexOffset+1],
                            amoebaGpu->torqueMapForce->_pSysData[indexOffset+2] );
            (void) fprintf( amoebaGpu->log,"T[%16.9e %16.9e %16.9e]\n",
                            psTorque->_pSysData[indexOffset],
                            psTorque->_pSysData[indexOffset+1],
                            psTorque->_pSysData[indexOffset+2] );
            if( ii == maxPrint && (gpu->natoms - maxPrint) > ii )ii = gpu->natoms - maxPrint;
        }
        (void) fflush( amoebaGpu->log );
    }
    if( 1 ){
        //std::vector<int> fileId;
        //fileId.push_back( 0 );
        VectorOfDoubleVectors outputVector;
        //cudaLoadCudaFloat4Array( gpu->natoms, 3, gpu->psPosq4,              outputVector );
        cudaLoadCudaFloat4Array( gpu->natoms, 4, gpu->psForce4,             outputVector, NULL );
        cudaLoadCudaFloatArray( gpu->natoms,  3, psForce,        outputVector, NULL );
        cudaLoadCudaFloatArray( gpu->natoms,  3, psTorque,       outputVector, NULL);
        cudaWriteVectorOfDoubleVectorsToFile( "CudaVacuumElecForce", fileId, outputVector );
    }
#endif

}

void cudaComputeAmoebaMapTorquesAndAddTotalForce2( amoebaGpuContext amoebaGpu,
                                                  CUDAStream<float>* psTorque,
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

    amoebaMapTorqueToForce_kernel<<< numBlocks, numThreads>>> (
                psTorque->_pDevData,
                amoebaGpu->maxMapTorqueDifference,
                amoebaGpu->torqueMapForce->_pDevData );
    LAUNCHERROR("AmoebaMapTrqKernel");

    numBlocks  = gpu->natoms;
    numThreads = amoebaGpu->maxMapTorqueDifferencePow2;

    amoebaMapTorqueReduce_kernel3<<< numBlocks, numThreads>>>(
            numThreads, gpu->natoms,
            amoebaGpu->maxMapTorqueDifference,
            amoebaGpu->torqueMapForce->_pDevData,
            psCudaForce4->_pDevData );
    LAUNCHERROR("amoebaMapTorqueReduce_kernel3");

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
                            psCudaForce4->_pSysData[ii].x,
                            psCudaForce4->_pSysData[ii].y,
                            psCudaForce4->_pSysData[ii].z );
            (void) fprintf( amoebaGpu->log,"F[%16.9e %16.9e %16.9e] ",
                            amoebaGpu->psForce->_pSysData[indexOffset],
                            amoebaGpu->psForce->_pSysData[indexOffset+1],
                            amoebaGpu->psForce->_pSysData[indexOffset+2] );
            (void) fprintf( amoebaGpu->log,"T[%16.9e %16.9e %16.9e]\n",
                            amoebaGpu->psTorque->_pSysData[indexOffset],
                            amoebaGpu->psTorque->_pSysData[indexOffset+1],
                            amoebaGpu->psTorque->_pSysData[indexOffset+2] );
            if( ii == maxPrint && (gpu->natoms - maxPrint) > ii )ii = gpu->natoms - maxPrint;
        }
        (void) fflush( amoebaGpu->log );
    }
    if( 1 ){
        //std::vector<int> fileId;
        //fileId.push_back( 0 );
        VectorOfDoubleVectors outputVector;
        //cudaLoadCudaFloat4Array( gpu->natoms, 3, gpu->psPosq4,              outputVector );
        cudaLoadCudaFloat4Array( gpu->natoms, 4, gpu->psForce4,             outputVector, NULL );
        cudaLoadCudaFloatArray( gpu->natoms,  3, amoebaGpu->psForce,        outputVector, NULL );
        cudaLoadCudaFloatArray( gpu->natoms,  3, amoebaGpu->psTorque,       outputVector, NULL);
        cudaWriteVectorOfDoubleVectorsToFile( "CudaVacuumElecForce", fileId, outputVector );
    }
#endif

}
