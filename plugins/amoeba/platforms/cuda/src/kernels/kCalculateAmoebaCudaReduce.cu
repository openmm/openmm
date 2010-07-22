
//-----------------------------------------------------------------------------------------

//-----------------------------------------------------------------------------------------

#include "amoebaGpuTypes.h"
#include "amoebaCudaKernels.h"
#include <stdio.h>

#undef AMOEBA_OFFSET_3
#undef AMOEBA_INCLUDE_DIAGONAL 
#define METHOD_NAME(a, b) a##ExcludeDiagonalOffset1##b
#include "kCalculateAmoebaCudaReduce.h"
#undef METHOD_NAME

#define AMOEBA_OFFSET_3
#define METHOD_NAME(a, b) a##ExcludeDiagonalOffset3##b
#include "kCalculateAmoebaCudaReduce.h"
#undef METHOD_NAME

#undef AMOEBA_OFFSET_3
#define AMOEBA_INCLUDE_DIAGONAL 
#define METHOD_NAME(a, b) a##IncludeDiagonalOffset1##b
#include "kCalculateAmoebaCudaReduce.h"
#undef METHOD_NAME

#define AMOEBA_OFFSET_3
#define METHOD_NAME(a, b) a##IncludeDiagonalOffset3##b
#include "kCalculateAmoebaCudaReduce.h"
#undef METHOD_NAME
#undef AMOEBA_OFFSET_3
#undef AMOEBA_INCLUDE_DIAGONAL 

void cudaReduceN2ToN( float *N2Array, int Nsz, float *NArray, int includeDiagonal, int offset )
{

    int numThreads        = min(THREADS_PER_BLOCK, (Nsz));
    int numBlocksPerAtom  = (Nsz / numThreads);

    if( Nsz % numThreads ){
        numBlocksPerAtom++;
    }

    int numBlocks = numBlocksPerAtom*Nsz;
    float *partialSum1_d;

    // allocate GPU memory 

    cudaMalloc( (void**) &partialSum1_d, numBlocks*offset*sizeof(float) );  

    if( includeDiagonal ){
       if( offset == 3 ){
           kCalculateAmoebaReduceIncludeDiagonalOffset3N2ToNBlockLevel<<< numBlocks, numThreads >>>( N2Array, partialSum1_d, Nsz, numBlocksPerAtom );
           LAUNCHERROR("kCalculateAmoebaReduceN2ToNBlockLevel1");
       } else if( offset == 1 ){
           kCalculateAmoebaReduceIncludeDiagonalOffset1N2ToNBlockLevel<<< numBlocks, numThreads >>>( N2Array, partialSum1_d, Nsz, numBlocksPerAtom );
           LAUNCHERROR("kCalculateAmoebaReduceN2ToNBlockLevel2");
       }
    } else {
       if( offset == 3 ){
           kCalculateAmoebaReduceExcludeDiagonalOffset3N2ToNBlockLevel<<< numBlocks, numThreads >>>( N2Array, partialSum1_d, Nsz, numBlocksPerAtom );
           LAUNCHERROR("kCalculateAmoebaReduceN2ToNBlockLevel3");
       } else if( offset == 1 ){
           kCalculateAmoebaReduceExcludeDiagonalOffset1N2ToNBlockLevel<<< numBlocks, numThreads >>>( N2Array, partialSum1_d, Nsz, numBlocksPerAtom );
           LAUNCHERROR("kCalculateAmoebaReduceN2ToNBlockLevel4");
       }
    }

    int numBlocks2 = numBlocks;
    numBlocks      = numBlocks2*Nsz/numThreads;
    if( (numBlocks2*Nsz) % numThreads ){
        numBlocks++;
    }

    if( offset == 3 ){
        kCalculateAmoebaReduceIncludeDiagonalOffset3N2ToNFinal<<< numBlocks, numThreads >>>(partialSum1_d, NArray, Nsz, numBlocksPerAtom );
        LAUNCHERROR("kCalculateAmoebaReduceN2ToNFinal3");
    } else if( offset == 1 ){
        kCalculateAmoebaReduceIncludeDiagonalOffset1N2ToNFinal<<< numBlocks, numThreads >>>(partialSum1_d, NArray, Nsz, numBlocksPerAtom );
        LAUNCHERROR("kCalculateAmoebaReduceN2ToNFinal1");
    }

    //Free memory

    cudaFree(partialSum1_d);
}
