/*
 * Authored by: Chen, Shifu
 * 
 * Email: chen@gmtk.org
 *
 * Website: http://www.gmtk.org/gsort
 *
 * The code is distributed under BSD license, you are allowed to use, modify or sell this code, but a statement is required if you used this code any where.
 * 
 */
#ifndef _BBSORT_KERNEL_H_
#define _BBSORT_KERNEL_H_

#include "bbsort.h"
#include "math_constants.h"

texture<unsigned int, 1, cudaReadModeElementType> tBucketSizes;
texture<unsigned int, 1, cudaReadModeElementType> tBucketOffsets;
texture<unsigned int, 1, cudaReadModeElementType> tBucketOfSlices;
texture<unsigned int, 1, cudaReadModeElementType> tSliceOffsetInBucket;

__device__ int dGetValue(int2 v){
	return v.y;
}

template <typename T>
__device__ T dGetValue(T v){
	return v;
}


__device__ void dPad(int2& v){
	v.x=0x3fffffff;
	v.y=0x4fffffff;
}

template <typename T>
__device__ void dPad(T & v){
	v=0x7fffffff;
}

template <typename T>
__global__ static void reduceMaxD(T * dData,int step,int length)
{
    int index = blockIdx.x * blockDim.x + threadIdx.x;

	if(index + step >=length)
		return ;
	dData[index] = dGetValue(dData[index])>dGetValue(dData[index+step])?dData[index]:dData[index+step];
}

template <typename T>
__global__ static void reduceMinD(T * dData,int step,int length)
{
	
    int index = blockIdx.x * blockDim.x + threadIdx.x;

	if(index + step >=length)
		return ;

	dData[index] = dGetValue(dData[index])<dGetValue(dData[index+step])?dData[index]:dData[index+step];
}

__global__ static void reduceSumD(float * dDiffData,int step,int length)
{

    int index = blockIdx.x * blockDim.x + threadIdx.x;

	if(index + step >=length)
		return ;

	dDiffData[index] += dDiffData[index+step];
}

template <typename T>
__global__ static void calDifferenceD(T * dData,float * dDiffData,int size)
{
    int index = blockIdx.x * blockDim.x + threadIdx.x;

	if(index > size-1)
		return ;

    const unsigned int tid = threadIdx.x;
	
	extern __shared__ T sData[];

	sData[tid]=dData[index];

	__syncthreads();
	
	if(tid < blockDim.x -1)
		dDiffData[index] = abs(dGetValue(sData[tid+1]) - dGetValue(sData[tid]));
	else 
		dDiffData[index] =0;
	
}

template <typename T>
__device__ inline void dSwap(T & a, T & b)
{
    T tmp = a;
    a = b;
    b = tmp;
}


template <typename T>
__global__ static void bitonicSortD(T * datas)
{
    extern __shared__ T shared[];

	const unsigned int bid=blockIdx.x;

    const unsigned int tid = threadIdx.x;

	__shared__ unsigned int count;
	__shared__ unsigned int offset;

	if(tid == 0)
	{
		count=tex1Dfetch(tBucketSizes,bid);
		offset=tex1Dfetch(tBucketOffsets,bid);
	}

	__syncthreads();

    if(tid < count)
		shared[tid] = datas[tid+offset];
	else 
	{
		dPad(shared[tid]);
	}

    __syncthreads();

    for (unsigned int k = 2; k <= BLOCK_SIZE; k *= 2)
    {
        for (unsigned int j = k / 2; j>0; j /= 2)
        {
            unsigned int ixj = tid ^ j;
            

            if (ixj > tid)
            {
                if ((tid & k) == 0)
                {
                    if (dGetValue(shared[tid]) > dGetValue(shared[ixj]))
                    {
                        dSwap(shared[tid], shared[ixj]);
                    }
                }
                else
                {
                    if (dGetValue(shared[tid]) < dGetValue(shared[ixj]))
                    {
                        dSwap(shared[tid], shared[ixj]);
                    }
                }
            }
            
            __syncthreads();
        }
    }
    if(tid < count)
		datas[tid+offset] = shared[tid];
}

template <typename T>

__global__ void assignElementToSlicesD(T* dDatas,int number,unsigned int* dSliceCounts,unsigned int* dOffsetInSlice,float minValue,float step,int sliceSize)
{
	unsigned int index= __mul24(blockIdx.x,blockDim.x) + threadIdx.x;

	if(index > number-1)
		return ;

	unsigned int s=((dGetValue(dDatas[index]) - minValue)/ step);

	unsigned int offset=atomicInc(dSliceCounts + s,0xFFFFFFF);

	dOffsetInSlice[index] = offset;

}

template <typename T>
__global__ void assignElementToSlicesNearlySortedD(T* dDatas,int number,unsigned int* dSliceCounts,unsigned int* dOffsetInSlice,float minValue,float step,int sliceSize,int blockCount)
{
	unsigned int index= blockIdx.x + blockCount * threadIdx.x;

	if(index > number-1)
		return ;

	unsigned int s=((dGetValue(dDatas[index]) - minValue)/ step);

	unsigned int offset=atomicInc(dSliceCounts + s,0xFFFFFFF);

	dOffsetInSlice[index] = offset;

}

template <typename T>
__global__ void assignElementToBucketD(T* dDatas,T*  dNewDatas,int number,unsigned int* dOffsetInSlice,float minValue,float step)
{

	unsigned int index= __mul24(blockIdx.x,blockDim.x) + threadIdx.x;

	if(index > number-1)
		return ;

	unsigned int s=((dGetValue(dDatas[index]) - minValue)/ step);

	unsigned int b=tex1Dfetch(tBucketOfSlices,s);

	unsigned int offset =tex1Dfetch(tBucketOffsets,b) + tex1Dfetch(tSliceOffsetInBucket,s) + dOffsetInSlice[index];

	dNewDatas[offset] =dDatas[index];

}

#endif // _BBSORT_KERNEL_H_
