
//-----------------------------------------------------------------------------------------

//-----------------------------------------------------------------------------------------

typedef unsigned int uint;

__global__
#if (__CUDA_ARCH__ >= 200)
__launch_bounds__(GF1XX_THREADS_PER_BLOCK, 1)
#elif (__CUDA_ARCH__ >= 130)
__launch_bounds__(GT2XX_THREADS_PER_BLOCK, 1)
#else
__launch_bounds__(G8X_THREADS_PER_BLOCK, 1)
#endif
void METHOD_NAME(kCalculateAmoebaReduce, N2ToNBlockLevel)( float *N2Array, float *partialSum, int num,int numberOfBlocksPerAtom )
{

   uint tid = threadIdx.x;
  
   __shared__ float asx[THREADS_PER_BLOCK];
   asx[tid] = 0.0f;

#ifdef AMOEBA_OFFSET_3
    __shared__ float asy[THREADS_PER_BLOCK];
    __shared__ float asz[THREADS_PER_BLOCK];
    asx[tid] = 0.0f;
    asy[tid] = asz[tid] = 0.0f;
    int offset = 3;
#else
    int offset = 1;
#endif

    int atomI =  blockIdx.x / numberOfBlocksPerAtom;
    int atomJ = (blockIdx.x % numberOfBlocksPerAtom)*blockDim.x+tid;

#ifdef AMOEBA_INCLUDE_DIAGONAL
    if( atomJ < num && atomI < num ){
#else
    if( atomJ < num && atomJ != atomI ){
#endif

      int index = offset*(atomI*num + atomJ);
      asx[tid] = N2Array[index];
#ifdef AMOEBA_OFFSET_3
      asy[tid] = N2Array[index+1];
      asz[tid] = N2Array[index+2];
#endif
    }
    __syncthreads(); //to make sure all the elements are loaded

    for( uint s = (blockDim.x)/2; s != 0; s >>= 1 ){
      if( tid < s ){
        asx[tid] += asx[tid+s];
#ifdef AMOEBA_OFFSET_3
        asy[tid] += asy[tid+s];
        asz[tid] += asz[tid+s];
#endif
      }
      __syncthreads();
    }
  
    if( tid == 0 ){
      partialSum[blockIdx.x*offset] = asx[0];
#ifdef AMOEBA_OFFSET_3
      partialSum[blockIdx.x*3+1]    = asy[0];
      partialSum[blockIdx.x*3+2]    = asz[0];
#endif
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
void METHOD_NAME(kCalculateAmoebaReduce, N2ToNFinal)( float *partialSum, float *final,int num,int numberOfBlocksPerAtom )
{

    uint thread_id = __mul24(blockIdx.x,blockDim.x) + threadIdx.x;
    if( thread_id > num )return;
  
    float3 sum;
  
#ifdef AMOEBA_OFFSET_3
    int offset = 3;
    sum.x = sum.y = sum.z = 0.0f;
#else
    int offset = 1;
    sum.x      = 0.0f;
#endif

    int index = thread_id*offset*numberOfBlocksPerAtom;
    for( int i=0; i < numberOfBlocksPerAtom; i++ ){
      sum.x += partialSum[index + i*offset];
#ifdef AMOEBA_OFFSET_3
      sum.y += partialSum[index + i*offset+1];
      sum.z += partialSum[index + i*offset+2];
#endif
      
    }
    final[thread_id*offset   ] = sum.x;
#ifdef AMOEBA_OFFSET_3
    final[thread_id*3+1 ]      = sum.y;
    final[thread_id*3+2 ]      = sum.z;
#endif
}

