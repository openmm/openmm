#ifndef CALCULATE_AMOEBA_CUDA_UTILITIES_H
#define CALCULATE_AMOEBA_CUDA_UTILITIES_H

#include "amoebaCudaKernels.h"

__global__ void kReduceFields_kernel( unsigned int fieldComponents, unsigned int outputBuffers, float* fieldIn, float* fieldOut, int addTo );
__global__ void kReduceAndCombineFields_kernel( unsigned int fieldComponents, unsigned int outputBuffers, float* fieldIn1, float* fieldIn2, float* fieldOut );
__global__ void kReduceFieldsToFloat4_kernel( unsigned int fieldComponents, unsigned int outputBuffers, float* fieldIn, float4* fieldOut );

extern __global__ void kFindBlockBoundsPeriodic_kernel();
extern __global__ void kFindBlocksWithInteractionsPeriodic_kernel();
//extern __global__ void kFindInteractionsWithinBlocksPeriodic_kernel(unsigned int*);


extern __global__ void kFindBlocksWithInteractionsVdwPeriodic_kernel();
extern __global__ void kFindInteractionsWithinBlocksVdwPeriodic_kernel(unsigned int*);


#endif
