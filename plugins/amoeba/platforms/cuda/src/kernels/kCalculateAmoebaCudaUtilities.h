#ifndef CALCULATE_AMOEBA_CUDA_UTILITIES_H
#define CALCULATE_AMOEBA_CUDA_UTILITIES_H

#include "amoebaCudaKernels.h"

__global__ void kReduceFields_kernel( unsigned int fieldComponents, unsigned int outputBuffers, float* fieldIn, float* fieldOut );
__global__ void kReduceAndCombineFields_kernel( unsigned int fieldComponents, unsigned int outputBuffers, float* fieldIn1, float* fieldIn2, float* fieldOut );
__global__ void kReduceFieldsToFloat4_kernel( unsigned int fieldComponents, unsigned int outputBuffers, float* fieldIn, float4* fieldOut );

#endif
