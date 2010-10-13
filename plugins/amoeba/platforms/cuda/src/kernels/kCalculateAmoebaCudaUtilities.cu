
//-----------------------------------------------------------------------------------------

//-----------------------------------------------------------------------------------------

#include "amoebaCudaKernels.h"
//#define AMOEBA_DEBUG

static __constant__ cudaGmxSimulation cSim;
static __constant__ cudaAmoebaGmxSimulation cAmoebaSim;

void SetCalculateAmoebaCudaUtilitiesSim(amoebaGpuContext amoebaGpu)
{
    cudaError_t status;
    gpuContext gpu = amoebaGpu->gpuContext;
    status         = cudaMemcpyToSymbol(cSim, &gpu->sim, sizeof(cudaGmxSimulation));    
    RTERROR(status, "SetCalculateAmoebaCudaUtilitiesSim: cudaMemcpyToSymbol: SetSim copy to cSim failed");
    status         = cudaMemcpyToSymbol(cAmoebaSim, &amoebaGpu->amoebaSim, sizeof(cudaAmoebaGmxSimulation));    
    RTERROR(status, "SetCalculateAmoebaCudaUtilitiesSim: cudaMemcpyToSymbol: SetSim copy to cAmoebaSim failed");
}

void GetCalculateAmoebaCudaUtilitiesSim(amoebaGpuContext amoebaGpu)
{
    cudaError_t status;
    gpuContext gpu = amoebaGpu->gpuContext;
    status = cudaMemcpyFromSymbol(&gpu->sim, cSim, sizeof(cudaGmxSimulation));    
    RTERROR(status, "GetCalculateAmoebaCudaUtilitiesSim: cudaMemcpyFromSymbol: SetSim copy from cSim failed");
    status = cudaMemcpyFromSymbol(&amoebaGpu->amoebaSim, cAmoebaSim, sizeof(cudaAmoebaGmxSimulation));    
    RTERROR(status, "GetCalculateAmoebaCudaUtilitiesSim: cudaMemcpyFromSymbol: SetSim copy from cAmoebaSim failed");
}

#undef METHOD_NAME
#define USE_PERIODIC
#define METHOD_NAME(a, b) a##Periodic##b
#include "kFindInteractingBlocks.h"
#undef METHOD_NAME
#undef USE_PERIODIC

__global__ 
#if (__CUDA_ARCH__ >= 200)
__launch_bounds__(GF1XX_THREADS_PER_BLOCK, 1)
#elif (__CUDA_ARCH__ >= 130)
__launch_bounds__(GT2XX_THREADS_PER_BLOCK, 1)
#else
__launch_bounds__(G8X_THREADS_PER_BLOCK, 1)
#endif
void kClearFloat4_kernel( unsigned int bufferLength, float4* fieldToClear )
{
    unsigned int pos = blockIdx.x * blockDim.x + threadIdx.x;
    while (pos < bufferLength )
    {   
        fieldToClear[pos].x       = 0.0f;
        fieldToClear[pos].y       = 0.0f;
        fieldToClear[pos].z       = 0.0f;
        fieldToClear[pos].w       = 0.0f;
        pos                      += gridDim.x * blockDim.x;
    }   
}

void kClearFloat4( amoebaGpuContext amoebaGpu, unsigned int entries, CUDAStream<float4>* fieldToClear )
{
    kClearFloat4_kernel<<<amoebaGpu->gpuContext->blocksPerSM, 384>>>( entries, fieldToClear->_pDevStream[0] );
    LAUNCHERROR("kClearFloat4");
}

__global__ 
#if (__CUDA_ARCH__ >= 200)
__launch_bounds__(GF1XX_THREADS_PER_BLOCK, 1)
#elif (__CUDA_ARCH__ >= 130)
__launch_bounds__(GT2XX_THREADS_PER_BLOCK, 1)
#else
__launch_bounds__(G8X_THREADS_PER_BLOCK, 1)
#endif
void kClearFloat_kernel( unsigned int bufferLength, float* fieldToClear )
{
    unsigned int pos = blockIdx.x * blockDim.x + threadIdx.x;
    while (pos < bufferLength )
    {   
        fieldToClear[pos]       = 0.0f;
        pos                    += gridDim.x * blockDim.x;
    }   
}

void kClearFloat( amoebaGpuContext amoebaGpu, unsigned int entries, CUDAStream<float>* fieldToClear )
{
    kClearFloat_kernel<<<amoebaGpu->gpuContext->blocksPerSM, 384>>>( entries, fieldToClear->_pDevStream[0] );
    LAUNCHERROR("kClearFloat");
}

__global__ 
#if (__CUDA_ARCH__ >= 200)
__launch_bounds__(GF1XX_THREADS_PER_BLOCK, 1)
#elif (__CUDA_ARCH__ >= 130)
__launch_bounds__(GT2XX_THREADS_PER_BLOCK, 1)
#else
__launch_bounds__(G8X_THREADS_PER_BLOCK, 1)
#endif
void kClearFields_kernel( unsigned int bufferLength, float* EField )
{
    unsigned int pos = blockIdx.x * blockDim.x + threadIdx.x;
    while (pos < bufferLength )
    {   
        EField[pos]       = 0.0f;
        pos              += gridDim.x * blockDim.x;
    }   
}

// clear psWorkArray_3_1 & psWorkArray_3_2

void kClearFields_3( amoebaGpuContext amoebaGpu, unsigned int numberToClear )
{

    kClearFields_kernel<<<amoebaGpu->nonbondBlocks, 384>>>( amoebaGpu->paddedNumberOfAtoms*3*amoebaGpu->outputBuffers,
                                                            amoebaGpu->psWorkArray_3_1->_pDevStream[0] );
    LAUNCHERROR("kClearFields_3_1");

    kClearFields_kernel<<<amoebaGpu->nonbondBlocks, 384>>>( amoebaGpu->paddedNumberOfAtoms*3*amoebaGpu->outputBuffers,
                                                            amoebaGpu->psWorkArray_3_2->_pDevStream[0] );
    LAUNCHERROR("kClearFields_3_2");

    if( numberToClear > 2 ){
        kClearFields_kernel<<<amoebaGpu->nonbondBlocks, 384>>>( amoebaGpu->paddedNumberOfAtoms*3*amoebaGpu->outputBuffers,
                                                            amoebaGpu->psWorkArray_3_3->_pDevStream[0] );
        LAUNCHERROR("kClearFields_3_3");
    } else {
        return;
    }
    if( numberToClear > 3 ){
        kClearFields_kernel<<<amoebaGpu->nonbondBlocks, 384>>>( amoebaGpu->paddedNumberOfAtoms*3*amoebaGpu->outputBuffers,
                                                                amoebaGpu->psWorkArray_3_4->_pDevStream[0] );
        LAUNCHERROR("kClearFields_3_4");
    }
}

// clear psWorkArray_1_1 & psWorkArray_1_2

void kClearFields_1( amoebaGpuContext amoebaGpu )
{

    kClearFields_kernel<<<amoebaGpu->nonbondBlocks, 384>>>( amoebaGpu->paddedNumberOfAtoms*amoebaGpu->outputBuffers,
                                                            amoebaGpu->psWorkArray_1_1->_pDevStream[0] );
    LAUNCHERROR("kClearFields_1_1");

    kClearFields_kernel<<<amoebaGpu->nonbondBlocks, 384>>>( amoebaGpu->paddedNumberOfAtoms*amoebaGpu->outputBuffers,
                                                            amoebaGpu->psWorkArray_1_2->_pDevStream[0] );
    LAUNCHERROR("kClearFields_1_2");
}

__global__ 
#if (__CUDA_ARCH__ >= 200)
__launch_bounds__(GF1XX_THREADS_PER_BLOCK, 1)
#elif (__CUDA_ARCH__ >= 130)
__launch_bounds__(GT2XX_THREADS_PER_BLOCK, 1)
#else
__launch_bounds__(G8X_THREADS_PER_BLOCK, 1)
#endif
void kReduceFields_kernel( unsigned int fieldComponents, unsigned int outputBuffers, float* fieldIn, float* fieldOut )
{
    unsigned int pos = blockIdx.x * blockDim.x + threadIdx.x;

    // Reduce field

    while (pos < fieldComponents)
    {   
        float totalField = 0.0f;
        float* pFt       = fieldIn + pos;
        unsigned int i   = outputBuffers;
        while (i >= 4)
        {   
            totalField += pFt[0] + pFt[fieldComponents] + pFt[2*fieldComponents] + pFt[3*fieldComponents];
            pFt        += fieldComponents*4;
            i          -= 4;
        }   

        if (i >= 2)
        {   
            totalField += pFt[0] + pFt[fieldComponents];
            pFt        += fieldComponents*2;
            i          -= 2;
        }   

        if (i > 0)
        {   
            totalField += pFt[0];
        }   

        fieldOut[pos]   = totalField;
        pos            += gridDim.x * blockDim.x;
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
void kReduceAndCombineFields_kernel( unsigned int fieldComponents, unsigned int outputBuffers, float* fieldIn1, float* fieldIn2, float* fieldOut )
{
    unsigned int pos = blockIdx.x * blockDim.x + threadIdx.x;

    // Reduce field

    while (pos < fieldComponents)
    {   
        float totalField = 0.0f;
        float* pFt1      = fieldIn1 + pos;
        float* pFt2      = fieldIn2 + pos;
        unsigned int i   = outputBuffers;
        while (i >= 4)
        {   
            totalField += pFt1[0] + pFt1[fieldComponents] + pFt1[2*fieldComponents] + pFt1[3*fieldComponents];
            totalField += pFt2[0] + pFt2[fieldComponents] + pFt2[2*fieldComponents] + pFt2[3*fieldComponents];
            pFt1       += fieldComponents*4;
            pFt2       += fieldComponents*4;
            i          -= 4;
        }   

        if (i >= 2)
        {   
            totalField += pFt1[0] + pFt1[fieldComponents];
            totalField += pFt2[0] + pFt2[fieldComponents];
            pFt1       += fieldComponents*2;
            pFt2       += fieldComponents*2;
            i          -= 2;
        }   

        if (i > 0)
        {   
            totalField += pFt1[0];
            totalField += pFt2[0];
        }   

        fieldOut[pos]   = totalField;
        pos            += gridDim.x * blockDim.x;
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
void kReduceFieldsToFloat4_kernel( unsigned int fieldComponents, unsigned int outputBuffers, float* fieldIn, float4* field4Out )
{
    unsigned int pos = blockIdx.x * blockDim.x + threadIdx.x;

    // Reduce field

    float* fieldOut = (float*) field4Out;
    while (pos < fieldComponents)
    {   
        float totalField = 0.0f;
        float* pFt       = fieldIn + pos;
        unsigned int i   = outputBuffers;
        while (i >= 4)
        {   
            totalField += pFt[0] + pFt[fieldComponents] + pFt[2*fieldComponents] + pFt[3*fieldComponents];
            pFt        += fieldComponents*4;
            i          -= 4;
        }   

        if (i >= 2)
        {   
            totalField += pFt[0] + pFt[fieldComponents];
            pFt        += fieldComponents*2;
            i          -= 2;
        }   

        if (i > 0)
        {   
            totalField += pFt[0];
        }   

        unsigned int j   = pos/3;
        unsigned int k   = pos - 3*j;
        fieldOut[4*j+k] += totalField;
        pos             += gridDim.x * blockDim.x;
    }   
}

