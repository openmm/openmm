
#include "cudaKernels.h"
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

__device__ static void loadMappedTorque( int particleId, int bufferIndex, float* forceToAdd )
{

    if( bufferIndex < 0 )return;
    unsigned int offset                     = particleId + bufferIndex*cSim.stride;
    float4 force                            = cAmoebaSim.pTorqueMapForce4[offset];
    force.x                                += forceToAdd[0];
    force.y                                += forceToAdd[1];
    force.z                                += forceToAdd[2];
    cAmoebaSim.pTorqueMapForce4[offset]     = force;
   

}

__global__
#if (__CUDA_ARCH__ >= 200)
__launch_bounds__(GF1XX_THREADS_PER_BLOCK, 1)
#elif (__CUDA_ARCH__ >= 120)
__launch_bounds__(GT2XX_THREADS_PER_BLOCK, 1)
#else
__launch_bounds__(G8X_THREADS_PER_BLOCK, 1)
#endif
void amoebaAddMapTorqueForceToForce_kernel( void )
{

  // ---------------------------------------------------------------------------------------
  
    int pos                      = blockIdx.x*blockDim.x + threadIdx.x;
  
  // ---------------------------------------------------------------------------------------
  
    while (pos < cSim.stride4 )
    {   
        float totalForce = 0.0f;
        float* pFt       = (float*)cAmoebaSim.pTorqueMapForce4 + pos;
        int i            = cAmoebaSim.maxTorqueBufferIndex;
        while (i >= 4)
        {   
            float f1    = *pFt;
            pFt        += cSim.stride4;
            float f2    = *pFt;
            pFt        += cSim.stride4;
            float f3    = *pFt;
            pFt        += cSim.stride4;
            float f4    = *pFt;
            pFt        += cSim.stride4;
            totalForce += f1 + f2 + f3 + f4; 
            i -= 4;
        }   
        if (i >= 2)
        {   
            float f1    = *pFt;
            pFt        += cSim.stride4;
            float f2    = *pFt;
            pFt        += cSim.stride4;
            totalForce += f1 + f2; 
            i -= 2;
        }   
        if (i > 0)
        {   
            totalForce += *pFt;
        }   

        pFt             = (float*)cSim.pForce4 + pos;
        *pFt           += totalForce;
        pos            += gridDim.x * blockDim.x;
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
void amoebaClearMapTorqueForce_kernel( void )
{

    int pos                      = blockIdx.x*blockDim.x + threadIdx.x;
    while (pos < cSim.stride4*cAmoebaSim.maxTorqueBufferIndex )
    {   
        cAmoebaSim.pTorqueMapForce4[pos]  = make_float4(0.0f, 0.0f, 0.0f, 0.0f);
        pos                              += gridDim.x * blockDim.x;
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
void amoebaMapTorqueToForce_kernel( float* torque )
{

  // ---------------------------------------------------------------------------------------
  
    int ii;
    int particleIndex            = blockIdx.x*blockDim.x + threadIdx.x;
    float4* atomCoord            = cSim.pPosq;
    int4* multiPoleAtoms         = cAmoebaSim.pMultipoleParticlesIdsAndAxisType;
  
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
  
    while( particleIndex < cSim.atoms )
    {   

        int axisAtom                  = multiPoleAtoms[particleIndex].z;
        int axisType                  = multiPoleAtoms[particleIndex].w;
    
        // NoAxisType
    
        if( axisType < 5 && multiPoleAtoms[particleIndex].z >= 0 )
        { 
        
            vector[U][0]                  = atomCoord[particleIndex].x - atomCoord[axisAtom].x;
            vector[U][1]                  = atomCoord[particleIndex].y - atomCoord[axisAtom].y;
            vector[U][2]                  = atomCoord[particleIndex].z - atomCoord[axisAtom].z;
        
            norms[U]                      = normVector3( vector[U] );
        
            if( axisType != 4 && multiPoleAtoms[particleIndex].x >= 0 ){
        
                axisAtom                  = multiPoleAtoms[particleIndex].x;
                vector[V][0]              = atomCoord[particleIndex].x - atomCoord[axisAtom].x;
                vector[V][1]              = atomCoord[particleIndex].y - atomCoord[axisAtom].y;
                vector[V][2]              = atomCoord[particleIndex].z - atomCoord[axisAtom].z;
        
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
                axisAtom                  = multiPoleAtoms[particleIndex].y;
            
                vector[W][0]              = atomCoord[particleIndex].x - atomCoord[axisAtom].x;
                vector[W][1]              = atomCoord[particleIndex].y - atomCoord[axisAtom].y;
                vector[W][2]              = atomCoord[particleIndex].z - atomCoord[axisAtom].z;
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
            dphi[U]                       = DOT3( vector[U], (torque + particleIndex*3) );
            dphi[V]                       = DOT3( vector[V], (torque + particleIndex*3) );
            dphi[W]                       = DOT3( vector[W], (torque + particleIndex*3) );
        
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
        
                float dphiR               = -1.0f*DOT3( vector[R], (torque + particleIndex*3) );
                float dphiS               = -1.0f*DOT3( vector[S], (torque + particleIndex*3) );
        
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
        
            // Z
        
            int4  forceBufferIndices      = cAmoebaSim.pMultipoleParticlesTorqueBufferIndices[particleIndex];
            loadMappedTorque( multiPoleAtoms[particleIndex].z, forceBufferIndices.z, forces[Z] );
           
            // X
        
            if( axisType != 4 ){
                loadMappedTorque( multiPoleAtoms[particleIndex].x, forceBufferIndices.x, forces[X] );
            }
         
            // Y
        
            if( axisType == 2 || axisType == 3 ){
                int particleId = multiPoleAtoms[particleIndex].y;
                if( particleId > -1 ){
                    loadMappedTorque( multiPoleAtoms[particleIndex].y, forceBufferIndices.y, forces[Y] );
                }
            }
         
            // put particle force in buffer 0
        
            loadMappedTorque( particleIndex, 0, forces[I] );
        }
        particleIndex   += gridDim.x*blockDim.x;
    }
}

void cudaComputeAmoebaMapTorqueAndAddToForce( amoebaGpuContext amoebaGpu, CUDAStream<float>* psTorque )
{
  
    gpuContext gpu    = amoebaGpu->gpuContext;

    if( amoebaGpu->amoebaSim.maxTorqueBufferIndex > amoebaGpu->gpuContext->sim.outputBuffers && amoebaGpu->psTorqueMapForce4 != amoebaGpu->gpuContext->psForce4 && amoebaGpu->psTorqueMapForce4 ){
        amoebaClearMapTorqueForce_kernel<<< gpu->sim.blocks, gpu->sim.threads_per_block>>> ( );
        LAUNCHERROR("amoebaClearMapTorqueForce");
    }

    if( 0 ){  
        VectorOfDoubleVectors outputVector;

        std::vector<int> fileId;
        static int call = 0;  
        fileId.push_back( call++ );

        int paddedNumberOfAtoms  = amoebaGpu->gpuContext->sim.paddedNumberOfAtoms;
        CUDAStream<float>* temp  = new CUDAStream<float>(3*paddedNumberOfAtoms, 1, "Temp1");

        reduceAndCopyCUDAStreamFloat4( gpu->psForce4, temp, 1.0 );
        cudaLoadCudaFloatArray( gpu->natoms,  3, temp, outputVector, NULL, 1.0f/4.184f );

        reduceAndCopyCUDAStreamFloat( psTorque, temp, 1.0 );
        cudaLoadCudaFloatArray( gpu->natoms,  3, temp, outputVector, NULL, 1.0f/4.184f );

        reduceAndCopyCUDAStreamFloat4( amoebaGpu->psTorqueMapForce4, temp, 1.0 );
        cudaLoadCudaFloatArray( gpu->natoms,  3, temp, outputVector, NULL, 1.0f/4.184f );

        cudaWriteVectorOfDoubleVectorsToFile( "CudaElectrostatiPreTorqueForce", fileId, outputVector );
        delete temp;
    }    

    amoebaMapTorqueToForce_kernel<<< gpu->sim.blocks, gpu->sim.threads_per_block>>> ( psTorque->_pDevData );
    LAUNCHERROR("amoebaMapTorqueToForce");

    if( amoebaGpu->amoebaSim.maxTorqueBufferIndex > amoebaGpu->gpuContext->sim.outputBuffers && amoebaGpu->psTorqueMapForce4 != amoebaGpu->gpuContext->psForce4 && amoebaGpu->psTorqueMapForce4 ){
        amoebaAddMapTorqueForceToForce_kernel<<< gpu->sim.blocks, gpu->sim.threads_per_block>>> ( );
        LAUNCHERROR("amoebaAddMapTorqueForceToForce");
    }

#ifdef AMOEBA_DEBUG
    if( 0 ){  
        VectorOfDoubleVectors outputVector;

        std::vector<int> fileId;
        static int call = 0;  
        fileId.push_back( call++ );

        int paddedNumberOfAtoms  = amoebaGpu->gpuContext->sim.paddedNumberOfAtoms;
        CUDAStream<float>* temp  = new CUDAStream<float>(3*paddedNumberOfAtoms, 1, "Temp1");

        reduceAndCopyCUDAStreamFloat4( gpu->psForce4, temp, 1.0 );
        cudaLoadCudaFloatArray( gpu->natoms,  3, temp, outputVector, NULL, 1.0f/4.184f );

        reduceAndCopyCUDAStreamFloat4( amoebaGpu->psTorqueMapForce4, temp, 1.0 );
        cudaLoadCudaFloatArray( gpu->natoms,  3, temp, outputVector, NULL, 1.0f/4.184f );
        for( int pId = 0; pId < 5; pId++ ){
        float sum[3] = { 0.0f,  0.0f,  0.0f };
        (void) fprintf( stderr, "\n\nTorqueForceToForce for part=%d\n", pId );
        for( int ii = 0; ii < amoebaGpu->amoebaSim.maxTorqueBufferIndex; ii++ ){
            (void) fprintf( stderr, "%4d [%15.7e %15.7e %15.7e]\n", ii,
                            amoebaGpu->psTorqueMapForce4->_pSysStream[ii][pId].x,
                            amoebaGpu->psTorqueMapForce4->_pSysStream[ii][pId].y,
                            amoebaGpu->psTorqueMapForce4->_pSysStream[ii][pId].z );
            sum[0] += amoebaGpu->psTorqueMapForce4->_pSysStream[ii][pId].x;
            sum[1] += amoebaGpu->psTorqueMapForce4->_pSysStream[ii][pId].y;
            sum[2] += amoebaGpu->psTorqueMapForce4->_pSysStream[ii][pId].z;
        }
        (void) fprintf( stderr, "TorqueForceToForce for partcle=%d [%15.7e %15.7e %15.7e]  [%15.7e %15.7e %15.7e]\n", pId, sum[0], sum[1], sum[2], sum[0]/4.184f, sum[1]/4.184f, sum[2]/4.184f );
        }

        cudaWriteVectorOfDoubleVectorsToFile( "CudaElectrostatiPostTorqueForce", fileId, outputVector );
        delete temp;
    }
#endif

}
